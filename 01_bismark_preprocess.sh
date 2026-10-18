#!/usr/bin/env bash
# Bisulfite-seq preprocessing: QC -> trimming -> Bismark alignment -> deduplication
# -> methylation extraction -> MultiQC summary.
#
# Usage:
#   bash 01_bismark_preprocess.sh samples.txt
#
# samples.txt holds one sample ID per line. Paired reads are expected at
#   ${READS_DIR}/<ID>/<ID>_1.fastq.gz and ${READS_DIR}/<ID>/<ID>_2.fastq.gz
#
# Settings can be overridden with environment variables, e.g.
#   THREADS=16 GENOME_DIR=ref/sorghum bash 01_bismark_preprocess.sh samples.txt

set -euo pipefail

READS_DIR=${READS_DIR:-data/raw}
GENOME_DIR=${GENOME_DIR:-ref/genome}     # folder containing the reference FASTA
OUT_DIR=${OUT_DIR:-results}
ADAPTERS=${ADAPTERS:-TruSeq3-PE.fa}
TRIMMOMATIC=${TRIMMOMATIC:-trimmomatic}  # or e.g. "java -jar /path/to/trimmomatic-0.39.jar"
THREADS=${THREADS:-8}
BISMARK_PARALLEL=${BISMARK_PARALLEL:-4}  # each Bismark instance uses several cores

samples_file=${1:?Usage: bash 01_bismark_preprocess.sh samples.txt}
mapfile -t samples < <(grep -v '^[[:space:]]*$' "$samples_file")

mkdir -p "$OUT_DIR"/{qc_raw,trimmed,qc_trimmed,bismark,dedup,methylation}

# 1. QC of raw reads
fastqc -t "$THREADS" -o "$OUT_DIR/qc_raw" "$READS_DIR"/*/*.fastq.gz

# 2. Bisulfite-convert and index the reference (only needed once)
if [[ ! -d "$GENOME_DIR/Bisulfite_Genome" ]]; then
    bismark_genome_preparation --bowtie2 --parallel "$BISMARK_PARALLEL" "$GENOME_DIR"
fi

for s in "${samples[@]}"; do
    echo "=== $s ==="
    r1="$READS_DIR/$s/${s}_1.fastq.gz"
    r2="$READS_DIR/$s/${s}_2.fastq.gz"
    t="$OUT_DIR/trimmed/$s"

    # 3. Adapter and quality trimming
    $TRIMMOMATIC PE -threads "$THREADS" "$r1" "$r2" \
        "${t}_R1.paired.fq.gz" "${t}_R1.unpaired.fq.gz" \
        "${t}_R2.paired.fq.gz" "${t}_R2.unpaired.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    fastqc -t "$THREADS" -o "$OUT_DIR/qc_trimmed" "${t}_R1.paired.fq.gz" "${t}_R2.paired.fq.gz"

    # 4. Alignment (output: <ID>_pe.bam)
    bismark --parallel "$BISMARK_PARALLEL" --genome "$GENOME_DIR" \
        -1 "${t}_R1.paired.fq.gz" -2 "${t}_R2.paired.fq.gz" \
        --basename "$s" -o "$OUT_DIR/bismark"

    # 5. Remove PCR duplicates (output: <ID>_pe.deduplicated.bam)
    deduplicate_bismark --paired --output_dir "$OUT_DIR/dedup" "$OUT_DIR/bismark/${s}_pe.bam"

    # 6. Methylation calls; produces <ID>_pe.deduplicated.bismark.cov.gz for methylKit
    bismark_methylation_extractor --paired-end --gzip --bedGraph \
        --parallel "$BISMARK_PARALLEL" -o "$OUT_DIR/methylation" \
        "$OUT_DIR/dedup/${s}_pe.deduplicated.bam"
done

# 7. One summary report across all samples (MultiQC reads FastQC, Trimmomatic and Bismark reports)
multiqc -f -o "$OUT_DIR/multiqc" "$OUT_DIR"

echo "Done. Coverage files for methylKit: $OUT_DIR/methylation/*.bismark.cov.gz"
