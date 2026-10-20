process BISMARK_DEDUP {
    // Why: Remove PCR duplicates so one DNA fragment does not count as many methylation calls.
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
