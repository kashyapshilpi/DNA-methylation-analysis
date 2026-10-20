#!/usr/bin/env nextflow
/*
 * Bisulfite-seq DNA methylation pipeline
 *
 *   FastQC -> Trimmomatic -> Bismark (Bowtie2) -> deduplicate_bismark
 *          -> bismark_methylation_extractor -> methylKit (DMCs / DMRs) -> MultiQC
 *
 * Usage:
 *   nextflow run kashyapshilpi/DNA-methylation-analysis -profile test,docker
 *   nextflow run kashyapshilpi/DNA-methylation-analysis -profile conda \
 *       --input samplesheet.csv --fasta genome.fa
 */

include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIMMOMATIC              } from './modules/trimmomatic'
include { BISMARK_GENOMEPREP       } from './modules/bismark_genomeprep'
include { BISMARK_ALIGN            } from './modules/bismark_align'
include { BISMARK_DEDUP            } from './modules/bismark_dedup'
include { BISMARK_EXTRACT          } from './modules/bismark_extract'
include { METHYLKIT                } from './modules/methylkit'
include { MULTIQC                  } from './modules/multiqc'

workflow {

    if (!params.input) { error "Missing --input (samplesheet CSV with columns: sample,fastq_1,fastq_2,condition)" }
    if (!params.fasta) { error "Missing --fasta (reference genome FASTA)" }

    log.info """
    D N A   M E T H Y L A T I O N   P I P E L I N E
    ===============================================
    input      : ${params.input}
    fasta      : ${params.fasta}
    contrast   : ${params.contrast_target} vs ${params.contrast_reference}
    outdir     : ${params.outdir}
    """.stripIndent()

    // fastq_2 may be empty for single-end samples
    ch_reads = channel.fromPath(file(params.input, checkIfExists: true))
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.condition) {
                error "Samplesheet rows need sample, fastq_1 and condition (fastq_2 optional): ${row}"
            }
            def single = !row.fastq_2
            def meta   = [id: row.sample, condition: row.condition, single_end: single]
            def reads  = single ? [file(row.fastq_1, checkIfExists: true)]
                                : [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
            [meta, reads]
        }

    fasta       = file(params.fasta, checkIfExists: true)
    adapters_pe = file(params.adapters_pe, checkIfExists: true)
    adapters_se = file(params.adapters_se, checkIfExists: true)
    genes_bed   = params.genes_bed ? file(params.genes_bed, checkIfExists: true) : []
    cpg_bed     = params.cpg_bed   ? file(params.cpg_bed,   checkIfExists: true) : []

    // QC and trimming
    FASTQC_RAW(ch_reads)
    TRIMMOMATIC(ch_reads, adapters_pe, adapters_se)
    FASTQC_TRIMMED(TRIMMOMATIC.out.reads)

    // Bisulfite alignment and methylation calling (Bismark handles the C->T conversion)
    BISMARK_GENOMEPREP(fasta)
    BISMARK_ALIGN(TRIMMOMATIC.out.reads, BISMARK_GENOMEPREP.out.index)
    BISMARK_DEDUP(BISMARK_ALIGN.out.bam)
    BISMARK_EXTRACT(BISMARK_DEDUP.out.bam)

    // Differential methylation: all samples are compared together, so collect them first
    ch_cov = BISMARK_EXTRACT.out.coverage
        .toList()
        .map { items -> [items.collect { it[0] }, items.collect { it[1] }] }
    METHYLKIT(ch_cov, genes_bed, cpg_bed)

    // Combined QC report
    ch_multiqc = FASTQC_RAW.out.reports
        .mix(FASTQC_TRIMMED.out.reports)
        .mix(TRIMMOMATIC.out.log)
        .mix(BISMARK_ALIGN.out.report)
        .mix(BISMARK_DEDUP.out.report)
        .mix(BISMARK_EXTRACT.out.reports)
        .collect()
    MULTIQC(ch_multiqc)
}
