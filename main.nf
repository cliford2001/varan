#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARAN } from './workflows/varan.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_varan_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_varan_pipeline'
include { getGenomeAttribute } from './subworkflows/local/utils_nfcore_varan_pipeline/main.nf'

params.fasta                = getGenomeAttribute('fasta')
params.fasta_fai            = getGenomeAttribute('fasta_fai')
params.dict                 = getGenomeAttribute('dict')
params.gtf                  = getGenomeAttribute('gtf')
params.gff                  = getGenomeAttribute('gff')
params.exon_bed             = getGenomeAttribute('exon_bed')
params.star_index           = getGenomeAttribute('star')
params.dbsnp                = getGenomeAttribute('dbsnp')
params.dbsnp_tbi            = getGenomeAttribute('dbsnp_tbi')
params.known_indels         = getGenomeAttribute('known_indels')
params.known_indels_tbi     = getGenomeAttribute('known_indels_tbi')
params.snpeff_db            = getGenomeAttribute('snpeff_db')
params.vep_cache_version    = getGenomeAttribute('vep_cache_version')
params.vep_genome           = getGenomeAttribute('vep_genome')
params.vep_species          = getGenomeAttribute('vep_species')

//
// WORKFLOW: Run main nf-core/rnavar analysis pipeline
//
workflow NFCORE_VARAN {
    VARAN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//

process GATKProcess {
    script:
    """
    gatk --version
    """
	print("hola")
}

workflow {
	main:
    VARAN ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
