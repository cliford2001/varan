g VARAN Pipeline

## Descripción
Este pipeline está diseñado para realizar análisis de variantes y procesos relacionados con el alineamiento de secuencias de RNA-Seq y la anotación de variantes. El flujo de trabajo está implementado usando **Nextflow** y sigue las mejores prácticas para la reproducibilidad y portabilidad. Incluye módulos y sub-workflows de **GATK4**, **snpEff**, **Ensembl VEP**, entre otros.

## Requisitos

- **Nextflow**: Para ejecutar el pipeline.
- **Docker/Singularity**: Para la portabilidad de las dependencias y el entorno de ejecución.
- **Git**: Para clonar el repositorio.
  
Además, este pipeline utiliza contenedores para asegurar la reproducibilidad de todos los pasos del análisis, incluyendo herramientas como:

- GATK4

## Estructura del Pipeline

El flujo de trabajo incluye los siguientes pasos principales:

1. **Preparación del genoma**:
   - Se crea el índice del genoma y otros archivos de referencia necesarios para el análisis.
   
2. **Validación de archivos de entrada**:
   - Se verifica y valida el archivo `samplesheet.csv` que contiene los datos de entrada y se preparan los canales de entrada.

3. **Llamada de variantes**:
   - Procesamiento de las secuencias mediante el algoritmo de GATK HaplotypeCaller para detectar variantes SNPs e indels.

4. **Filtrado de variantes**:
   - Las variantes detectadas se filtran mediante GATK VariantFiltration.

5. **Generación de informes**:
   - Se recogen todos los informes de calidad de control (QC) y versiones de software usadas, que se recopilan para su posterior análisis con **MultiQC**.

## Parámetros de Entrada

El pipeline requiere los siguientes parámetros de entrada:

- `--input`: Archivo de entrada en formato CSV con las muestras a procesar.
- `--fasta`: Archivo FASTA del genoma de referencia.
- `--fasta_fai`: Índice del archivo FASTA.
- `--dict`: Archivo de diccionario del genoma de referencia.
- `--gtf/gff`: Anotaciones de genes en formato GTF o GFF.
- `--dbsnp`: Archivo VCF de variantes conocidas de dbSNP.
- `--known_indels`: Archivo VCF de indels conocidos.
- `--snpeff_cache`: Directorio de caché para snpEff.
- `--star_index`: Índice STAR para alineamiento de RNA-Seq.

## Ejecución del Pipeline

Para ejecutar el pipeline, asegúrate de tener **Nextflow** instalado, y luego usa el siguiente comando:

```bash
nextflow run varan.nf --input samplesheet.csv --fasta genome.fasta --dbsnp dbsnp.vcf --snpeff_cache snpeff_cache/

