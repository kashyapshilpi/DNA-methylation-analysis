process BISMARK_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda 'bioconda::bismark=0.24'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}*.bam"), emit: bam
    path "${meta.id}_*_report.txt",           emit: report

    script:
    def r        = [reads].flatten()
    def input    = meta.single_end ? "${r[0]}" : "-1 ${r[0]} -2 ${r[1]}"
    def nondir   = params.non_directional ? '--non_directional' : ''
    def parallel = params.bismark_parallel > 1 ? "--parallel ${params.bismark_parallel}" : ''
    """
    bismark --genome ${index} ${input} --basename ${meta.id} ${parallel} ${nondir}
    """

    stub:
    def bam    = meta.single_end ? "${meta.id}.bam" : "${meta.id}_pe.bam"
    def report = meta.single_end ? "${meta.id}_SE_report.txt" : "${meta.id}_PE_report.txt"
    """
    touch ${bam} ${report}
    """
}
