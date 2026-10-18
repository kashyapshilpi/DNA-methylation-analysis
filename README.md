# DNA Methylation Analysis (Bisulfite-seq)

Paired-end whole-genome bisulfite sequencing workflow, from raw reads to differentially methylated CpGs and regions with gene and CpG-island annotation.

> Developed during my MSc dissertation at CSIR–IHBT. Kept as a reference implementation; not actively maintained.

## Workflow

```
FASTQ ─► FastQC ─► Trimmomatic ─► Bismark (Bowtie2) ─► deduplicate_bismark
      ─► bismark_methylation_extractor ─► methylKit
          ├─ per-sample QC, correlation, clustering, PCA
          ├─ DMCs (overdispersion-corrected χ² test)
          ├─ DMRs (1 kb tiles)
          └─ annotation: gene parts, nearest TSS, CpG islands / shores
```

## Repository contents

| File | Purpose |
|---|---|
| `01_bismark_preprocess.sh` | QC, trimming, alignment, deduplication, methylation extraction, MultiQC report |
| `02_methylkit_dmc.R` | Differential methylation, plots and annotation with methylKit |
| `GeneDMRs.R` | Alternative DMR analysis with GeneDMRs |
| `samplesheet.example.csv` | Example sample sheet for step 2 |

## Requirements

- FastQC, MultiQC
- Trimmomatic 0.39
- Bismark 0.24 with Bowtie2
- R ≥ 4.1 with `methylKit`, `genomation`, `GenomicRanges` (Bioconductor)

Example install with conda:

```bash
conda create -n methylation -c bioconda -c conda-forge \
    fastqc multiqc trimmomatic bismark bowtie2 \
    bioconductor-methylkit bioconductor-genomation
```

## Usage

**Step 1: preprocessing.** List sample IDs in `samples.txt`, one per line, with reads at `data/raw/<ID>/<ID>_1.fastq.gz` and `_2.fastq.gz`:

```bash
GENOME_DIR=ref/genome ADAPTERS=TruSeq3-PE.fa THREADS=16 \
    bash 01_bismark_preprocess.sh samples.txt
```

Outputs land in `results/`. Coverage files for methylKit are written to `results/methylation/*.bismark.cov.gz`, and the combined QC report to `results/multiqc/`.

**Step 2: differential methylation.** Fill in a sample sheet (see `samplesheet.example.csv`), then:

```bash
ASSEMBLY=Sorghum_bicolor_NCBIv3 CORES=4 \
    Rscript 02_methylkit_dmc.R samplesheet.csv genes.bed cpg_islands.bed results/methylkit
```

## Outputs (`results/methylkit/`)

| File | Contents |
|---|---|
| `per_sample_qc.pdf` | Methylation and coverage distributions per sample |
| `sample_similarity.pdf` | Correlation matrix, hierarchical clustering, PCA |
| `DMC_all.tsv`, `DMC_hyper.tsv`, `DMC_hypo.tsv` | Differentially methylated CpGs (≥25% difference, q < 0.01) |
| `DMR_1kb_tiles.tsv` | Differentially methylated 1 kb regions |
| `DMC_nearest_TSS.tsv` | Nearest transcript start site for each DMC |
| `DMC_annotation.pdf`, `DMC_per_chromosome.pdf` | Annotation and per-chromosome summaries |

## Data

Raw data are not included. The analysis was run on public bisulfite-seq data from NCBI BioProject PRJNA289892.

## License

MIT. See [licence.txt](licence.txt).
