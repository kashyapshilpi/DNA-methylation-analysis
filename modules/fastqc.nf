process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::fastqc=0.12.1'

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    """

    stub:
    """
    touch ${[reads].flatten()[0].simpleName}_fastqc.zip ${[reads].flatten()[0].simpleName}_fastqc.html
    """
}
