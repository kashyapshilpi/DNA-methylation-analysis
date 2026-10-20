process BISMARK_GENOMEPREP {
    // Why: Make C->T and G->A converted copies of the genome, since bisulfite turns unmethylated C into T.
    tag "${fasta.simpleName}"
    label 'process_high'
    conda 'bioconda::bismark=0.24'

    input:
    path fasta

    output:
    path 'BismarkIndex', emit: index

    script:
    """
    mkdir BismarkIndex
    cp -L ${fasta} BismarkIndex/
    bismark_genome_preparation --bowtie2 BismarkIndex
    """

    stub:
    """
    mkdir -p BismarkIndex/Bisulfite_Genome
    """
}
