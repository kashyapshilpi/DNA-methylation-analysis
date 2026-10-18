process TRIMMOMATIC {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::trimmomatic=0.39'

    input:
    tuple val(meta), path(reads)
    path adapters_pe
    path adapters_se

    output:
    tuple val(meta), path("${meta.id}*.trimmed.fq.gz"), emit: reads
    path "${meta.id}.trimmomatic.log",                   emit: log

    script:
    def r     = [reads].flatten()
    def p     = meta.id
    def steps = "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:${params.min_len}"
    def cmd   = meta.single_end
        ? "trimmomatic SE -threads ${task.cpus} ${r[0]} ${p}.trimmed.fq.gz ILLUMINACLIP:${adapters_se}:2:30:10 ${steps}"
        : "trimmomatic PE -threads ${task.cpus} ${r[0]} ${r[1]} ${p}_R1.trimmed.fq.gz ${p}_R1.unpaired.fq.gz ${p}_R2.trimmed.fq.gz ${p}_R2.unpaired.fq.gz ILLUMINACLIP:${adapters_pe}:2:30:10:2:True ${steps}"
    """
    ${cmd} 2> ${p}.trimmomatic.log
    rm -f ${p}_R1.unpaired.fq.gz ${p}_R2.unpaired.fq.gz
    """

    stub:
    def outs = meta.single_end ? "${meta.id}.trimmed.fq.gz" : "${meta.id}_R1.trimmed.fq.gz ${meta.id}_R2.trimmed.fq.gz"
    """
    touch ${outs} ${meta.id}.trimmomatic.log
    """
}
