process METHYLKIT {
    label 'process_medium'
    conda 'bioconda::bioconductor-methylkit bioconda::bioconductor-genomation'

    input:
    tuple val(metas), path(covs)
    path genes_bed      // optional: pass [] to skip
    path cpg_bed        // optional: pass [] to skip

    output:
    path 'methylkit/*', emit: results

    script:
    def files = [covs].flatten()
    def rows  = [metas, files].transpose().collect { m, f -> "${m.id},${f.name},${m.condition}" }.join(' ')
    def genes = genes_bed ? "${genes_bed}" : 'NA'
    def cpg   = cpg_bed   ? "${cpg_bed}"   : 'NA'
    """
    printf '%s\\n' sample_id,coverage_file,condition ${rows} > samplesheet.csv

    MIN_COV=${params.min_cov} DIFF=${params.diff_cutoff} QVAL=${params.qvalue} TILE=${params.tile_size} \\
    ASSEMBLY=${params.assembly} CORES=${task.cpus} \\
        methylkit_dmc.R samplesheet.csv ${genes} ${cpg} methylkit \\
        ${params.contrast_reference} ${params.contrast_target}
    """

    stub:
    """
    mkdir methylkit
    touch methylkit/DMC_all.tsv methylkit/DMR_tiles.tsv methylkit/sample_similarity.pdf
    """
}
