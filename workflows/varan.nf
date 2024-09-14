// Validate input parameters
//WorkflowRnavar.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.gtf,
    params.gff,
    params.dbsnp,
    params.dbsnp_tbi,
    params.known_indels,
    params.known_indels_tbi,
    params.snpeff_cache,
    params.vep_cache,
    params.star_index]

for (param in checkPathParamList) {if (param) file(param, checkIfExists: true)}
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

include { PREPARE_GENOME                } from '../subworkflows/local/prepare_genome'           // Build the genome index and other reference files
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'              // Validate the input samplesheet.csv and prepare input channels
include { ANNOTATE                      } from '../subworkflows/local/annotate'

include { GATK4_BEDTOINTERVALLIST       } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS       } from '../modules/nf-core/gatk4/intervallisttools/main'
include { MARKDUPLICATES                } from '../subworkflows/nf-core/markduplicates'     // Mark duplicates in the BAM file
include { GATK4_BASERECALIBRATOR        } from '../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_HAPLOTYPECALLER         } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_MERGEVCFS               } from '../modules/nf-core/gatk4/mergevcfs/main'
include { SPLITNCIGAR                   } from '../subworkflows/nf-core/splitncigar'
include { TABIX_TABIX as TABIX          } from '../modules/nf-core/tabix/tabix/main'
include { GATK4_VARIANTFILTRATION       } from '../modules/nf-core/gatk4/variantfiltration/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'

ch_dbsnp                = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()               : Channel.empty()
ch_dbsnp_tbi            = params.dbsnp_tbi         ? Channel.fromPath(params.dbsnp_tbi).collect()           : Channel.empty()
ch_known_indels         = params.known_indels      ? Channel.fromPath(params.known_indels).collect()        : Channel.empty()
ch_known_indels_tbi     = params.known_indels_tbi  ? Channel.fromPath(params.known_indels_tbi).collect()    : Channel.empty()

workflow VARAN {
    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()


	PREPARE_GENOME()
	ch_genome_bed = Channel.from([id:'genome.bed']).combine(PREPARE_GENOME.out.exon_bed)
	ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)


	//
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

	INPUT_CHECK (
    	ch_input,
	)
	.reads
	.set { ch_genome_bam }

	ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

	// Create samplesheet channel (after input check)
	ch_samplesheet = Channel.fromPath(params.input)

    
    // MODULE: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
    //
    ch_interval_list = Channel.empty()
    GATK4_BEDTOINTERVALLIST(
        ch_genome_bed,
        PREPARE_GENOME.out.dict
    )
    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions.first().ifEmpty(null))

    //
    // MODULE: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
	
    ch_interval_list_split = Channel.empty()
    if (!params.skip_intervallisttools) {
        GATK4_INTERVALLISTTOOLS(
            ch_interval_list
        )
        ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ meta, bed -> [bed] }.flatten()
	}
    else ch_interval_list_split = ch_interval_list


		
	MARKDUPLICATES (
		ch_genome_bam
    )
	
	ch_genome_bam             = MARKDUPLICATES.out.bam_bai
	
    //Gather QC reports
    //ch_reports                = ch_reports.mix(MARKDUPLICATES.out.stats.collect{it[1]}.ifEmpty([]))
    ch_reports                = ch_reports.mix(MARKDUPLICATES.out.metrics.collect{it[1]}.ifEmpty([]))
    ch_versions               = ch_versions.mix(MARKDUPLICATES.out.versions.first().ifEmpty(null))

    // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
    // Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
    //
    ch_splitncigar_bam_bai = Channel.empty()
    SPLITNCIGAR (
        ch_genome_bam,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        ch_interval_list_split
    )
    ch_splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
    ch_versions             = ch_versions.mix(SPLITNCIGAR.out.versions.first().ifEmpty(null))


    //
    // MODULE: BaseRecalibrator from GATK4
    // Generates a recalibration table based on various co-variates
    //
    ch_bam_variant_calling = Channel.empty()
    if(!params.skip_baserecalibration) {
        ch_bqsr_table   = Channel.empty()
        // known_sites is made by grouping both the dbsnp and the known indels ressources
        // they can either or both be optional
        ch_known_sites     = ch_dbsnp.concat(ch_known_indels).collect()
        ch_known_sites_tbi = ch_dbsnp_tbi.concat(ch_known_indels_tbi).collect()

        ch_interval_list_recalib = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
        ch_splitncigar_bam_bai.combine(ch_interval_list_recalib)
            .map{ meta, bam, bai, interval -> [ meta, bam, bai, interval]
        }.set{ch_splitncigar_bam_bai_interval}

        GATK4_BASERECALIBRATOR(
            ch_splitncigar_bam_bai_interval,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict,
            ch_known_sites,
            ch_known_sites_tbi
        )
        ch_bqsr_table   = GATK4_BASERECALIBRATOR.out.table
        // Gather QC reports
        ch_reports  = ch_reports.mix(ch_bqsr_table.map{ meta, table -> table})
        ch_versions     = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first().ifEmpty(null))

        ch_bam_applybqsr       = ch_splitncigar_bam_bai.join(ch_bqsr_table, by: [0])
        ch_bam_recalibrated_qc = Channel.empty()

        ch_interval_list_applybqsr = ch_interval_list.map{ meta, bed -> [bed] }.flatten()
        ch_bam_applybqsr.combine(ch_interval_list_applybqsr)
            .map{ meta, bam, bai, table, interval -> [ meta, bam, bai, table, interval]
        }.set{ch_applybqsr_bam_bai_interval}

        //
        // MODULE: ApplyBaseRecalibrator from GATK4
        // Recalibrates the base qualities of the input reads based on the recalibration table produced by the GATK BaseRecalibrator tool.
        //
        RECALIBRATE(
            params.skip_multiqc,
            ch_applybqsr_bam_bai_interval,
            PREPARE_GENOME.out.dict,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.fasta
        )

        ch_bam_variant_calling = RECALIBRATE.out.bam
        ch_bam_recalibrated_qc = RECALIBRATE.out.qc

        // Gather QC reports
        ch_reports  = ch_reports.mix(RECALIBRATE.out.qc.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(RECALIBRATE.out.versions.first().ifEmpty(null))
    } else {
        ch_bam_variant_calling = ch_splitncigar_bam_bai
    }

    interval_flag = params.no_intervals


	// Run haplotyper even in the absence of dbSNP files
    if (!params.dbsnp){
        ch_dbsnp = []
        ch_dbsnp_tbi = []
    }

    ch_haplotypecaller_vcf = Channel.empty()
    ch_haplotypecaller_interval_bam = ch_bam_variant_calling.combine(ch_interval_list_split)

	bam_intervals = ch_haplotypecaller_interval_bam.map { bam_data, bam_file1, bam_file2, intervals_data, intervals_file ->
    	[bam_data, bam_file1, bam_file2, intervals_file]
	}



    //
    // MODULE: HaplotypeCaller from GATK4
    // Calls germline SNPs and indels via local re-assembly of haplotypes.
    //
    GATK4_HAPLOTYPECALLER(
        bam_intervals,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        ch_dbsnp,
        ch_dbsnp_tbi
    )

    ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
        .map{ meta, vcf ->
            meta.id = meta.sample
            [meta, vcf]}
        .groupTuple()

    ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first().ifEmpty(null))


	//
    // MODULE: MergeVCFS from GATK4
    // Merge multiple VCF files into one VCF
    //
    GATK4_MERGEVCFS(
        ch_haplotypecaller_raw,
        PREPARE_GENOME.out.dict
    )
    ch_haplotypecaller_vcf = GATK4_MERGEVCFS.out.vcf
    ch_versions  = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first().ifEmpty(null))

    //
    // MODULE: Index the VCF using TABIX
    //
    TABIX(
        ch_haplotypecaller_vcf
    )

    ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_vcf
        .join(TABIX.out.tbi, by: [0], remainder: true)
        .join(TABIX.out.csi, by: [0], remainder: true)
        .map{meta, vcf, tbi, csi ->
            if (tbi) [meta, vcf, tbi]
            else [meta, vcf, csi]
        }

    ch_versions     = ch_versions.mix(TABIX.out.versions.first().ifEmpty(null))
    ch_final_vcf    = ch_haplotypecaller_vcf

    //
    // MODULE: VariantFiltration from GATK4
    // Filter variant calls based on certain criteria
    //
    if (!params.skip_variantfiltration && !params.bam_csi_index ) {

        GATK4_VARIANTFILTRATION(
            ch_haplotypecaller_vcf_tbi,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.dict
        )

        ch_filtered_vcf = GATK4_VARIANTFILTRATION.out.vcf
        ch_final_vcf    = ch_filtered_vcf
        ch_versions     = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.first().ifEmpty(null))
    }

    //
    // SUBWORKFLOW: Annotate variants using snpEff and Ensembl VEP if enabled.
    //
    if((!params.skip_variantannotation) && (params.annotate_tools) && (params.annotate_tools.contains('merge') || params.annotate_tools.contains('snpeff') || params.annotate_tools.contains('vep'))) {
        ANNOTATE(
            ch_final_vcf,
            params.annotate_tools,
            ch_snpeff_db,
            ch_snpeff_cache,
            ch_vep_genome,
            ch_vep_species,
            ch_vep_cache_version,
            ch_vep_cache)

        // Gather QC reports
        ch_reports  = ch_reports.mix(ANNOTATE.out.reports)
        ch_versions = ch_versions.mix(ANNOTATE.out.versions.first().ifEmpty(null))
    }

    ch_version_yaml = Channel.empty()
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()



}
