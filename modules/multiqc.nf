process MULTIQC {
    label 'process_low'
    conda 'bioconda::multiqc=1.25'

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data',        emit: data

    script:
    """
    multiqc --force --title 'DNA methylation pipeline' .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html
    """
}
