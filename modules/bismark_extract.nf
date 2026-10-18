process BISMARK_EXTRACT {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bismark=0.24'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bismark.cov.gz"),           emit: coverage
    path "*.bedGraph.gz",                                emit: bedgraph
    path "*{_splitting_report.txt,.M-bias.txt}",         emit: reports

    script:
    def mode = meta.single_end ? '-s' : '-p'
    """
    bismark_methylation_extractor ${mode} --gzip --bedGraph ${bam}
    # per-read context files are large and not needed downstream
    rm -f CpG_context_*.txt.gz CHG_context_*.txt.gz CHH_context_*.txt.gz
    """

    stub:
    """
    echo -e "chr\\t1\\t1\\t80\\t8\\t2" | gzip > ${bam.baseName}.bismark.cov.gz
    touch ${bam.baseName}.bedGraph.gz ${bam.baseName}_splitting_report.txt ${bam.baseName}.M-bias.txt
    """
}
