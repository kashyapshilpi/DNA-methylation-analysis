#!/usr/bin/env Rscript
# Differential CpG methylation between two groups with methylKit.
#
# Usage:
#   Rscript 02_methylkit_dmc.R samplesheet.csv genes.bed cpg_islands.bed results/methylkit
#
# samplesheet.csv columns: sample_id,coverage_file,treatment
#   coverage_file = Bismark *.bismark.cov.gz from 01_bismark_preprocess.sh
#   treatment     = 1 for the test group, 0 for the control group
# genes.bed       = 12-column BED of transcripts for the same genome
# cpg_islands.bed = BED of CpG islands (not a whole-genome or gene BED)
#
# Optional environment variables: ASSEMBLY (label only), CORES

suppressPackageStartupMessages({
    library(methylKit)
    library(genomation)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript 02_methylkit_dmc.R samplesheet.csv genes.bed cpg_islands.bed out_dir")
}
ss       <- read.csv(args[1], stringsAsFactors = FALSE)
gene_bed <- args[2]
cpg_bed  <- args[3]
out_dir  <- args[4]

assembly  <- Sys.getenv("ASSEMBLY", "Sorghum_bicolor_NCBIv3")
cores     <- as.integer(Sys.getenv("CORES", "2"))
min_cov   <- 10    # minimum reads per CpG
diff_cut  <- 25    # minimum % methylation difference
q_cut     <- 0.01  # maximum q-value

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(out_dir, x)
write_tsv <- function(x, f) write.table(x, out(f), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Read Bismark coverage files ------------------------------------------
obj <- methRead(as.list(ss$coverage_file),
                sample.id = as.list(ss$sample_id),
                assembly  = assembly,
                treatment = ss$treatment,
                context   = "CpG",
                pipeline  = "bismarkCoverage",
                mincov    = min_cov)

# ---- Per-sample QC ---------------------------------------------------------
pdf(out("per_sample_qc.pdf"), width = 10, height = 8)
for (i in seq_along(obj)) {
    getMethylationStats(obj[[i]], plot = TRUE, both.strands = FALSE)
    getCoverageStats(obj[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

# ---- Filter, normalise, merge ------------------------------------------------
# Drop extreme-coverage CpGs (likely PCR artefacts), then equalise library depth.
filt <- filterByCoverage(obj, lo.count = min_cov, hi.perc = 99.9)
filt <- normalizeCoverage(filt)
meth <- unite(filt, destrand = FALSE)

# ---- Sample similarity -------------------------------------------------------
pdf(out("sample_similarity.pdf"), width = 10, height = 10)
getCorrelation(meth, plot = TRUE)
clusterSamples(meth, dist = "correlation", method = "ward.D", plot = TRUE)
PCASamples(meth, screeplot = TRUE)
PCASamples(meth)
dev.off()

# ---- Differentially methylated CpGs (DMCs) -----------------------------------
# Overdispersion correction is appropriate with biological replicates.
diff <- calculateDiffMeth(meth, overdispersion = "MN", test = "Chisq", mc.cores = cores)

dmc       <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut)
dmc_hyper <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut, type = "hyper")
dmc_hypo  <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut, type = "hypo")

write_tsv(getData(dmc),       "DMC_all.tsv")
write_tsv(getData(dmc_hyper), "DMC_hyper.tsv")
write_tsv(getData(dmc_hypo),  "DMC_hypo.tsv")

pdf(out("DMC_per_chromosome.pdf"), width = 10, height = 8)
diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q_cut, meth.cutoff = diff_cut)
dev.off()

# ---- Differentially methylated regions (1 kb tiles) --------------------------
tiles  <- tileMethylCounts(filt, win.size = 1000, step.size = 1000, cov.bases = 10)
meth_t <- unite(tiles, destrand = FALSE)
diff_t <- calculateDiffMeth(meth_t, overdispersion = "MN", test = "Chisq", mc.cores = cores)
dmr    <- getMethylDiff(diff_t, difference = diff_cut, qvalue = q_cut)
write_tsv(getData(dmr), "DMR_1kb_tiles.tsv")

# ---- Annotation --------------------------------------------------------------
genes  <- readTranscriptFeatures(gene_bed)
cpg    <- readFeatureFlank(cpg_bed, feature.flank.name = c("CpGi", "shores"))
dmc_gr <- as(dmc, "GRanges")

ann_genes <- annotateWithGeneParts(dmc_gr, genes)
ann_cpg   <- annotateWithFeatureFlank(dmc_gr, feature = cpg$CpGi, flank = cpg$shores,
                                      feature.name = "CpGi", flank.name = "shores")

write_tsv(getAssociationWithTSS(ann_genes), "DMC_nearest_TSS.tsv")

pdf(out("DMC_annotation.pdf"), width = 10, height = 10)
plotTargetAnnotation(ann_genes, precedence = TRUE, main = "DMCs by gene part")
plotTargetAnnotation(ann_cpg, col = c("green", "gray", "white"),
                     main = "DMCs by CpG-island context")
dev.off()

message("Done. Results written to ", normalizePath(out_dir))
