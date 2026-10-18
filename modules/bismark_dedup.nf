process BISMARK_DEDUP {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bismark=0.24'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.deduplicated.bam"), emit: bam
    path "*.deduplication_report.txt",           emit: report

    script:
    def mode = meta.single_end ? '-s' : '-p'
    """
    deduplicate_bismark ${mode} --bam ${bam}
    """

    stub:
    """
    touch ${bam.baseName}.deduplicated.bam ${bam.baseName}.deduplication_report.txt
    """
}
