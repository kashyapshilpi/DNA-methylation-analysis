#!/usr/bin/env Rscript
# Differential CpG methylation between two conditions with methylKit
# (called by the Nextflow pipeline, also usable on its own).
#
# Usage:
#   methylkit_dmc.R samplesheet.csv genes.bed|NA cpg_islands.bed|NA out_dir [reference] [target]
#
# samplesheet.csv : columns sample_id,coverage_file,condition
#                   coverage_file = Bismark *.bismark.cov.gz
# genes.bed       : optional 12-column transcript BED for gene-part annotation (or NA)
# cpg_islands.bed : optional CpG-island BED for island/shore annotation (or NA)
#
# Environment variables: MIN_COV (10), DIFF (25), QVAL (0.01), TILE (1000), ASSEMBLY, CORES (1)

suppressPackageStartupMessages({
    library(methylKit)
    library(genomation)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: methylkit_dmc.R samplesheet.csv genes.bed|NA cpg_islands.bed|NA out_dir [reference] [target]")
}
ss        <- read.csv(args[1], stringsAsFactors = FALSE, strip.white = TRUE)
gene_bed  <- if (args[2] == "NA") NA else args[2]
cpg_bed   <- if (args[3] == "NA") NA else args[3]
out_dir   <- args[4]
ref_level <- if (length(args) >= 5) args[5] else "control"
tgt_level <- if (length(args) >= 6) args[6] else "treatment"

min_cov  <- as.integer(Sys.getenv("MIN_COV", "10"))
diff_cut <- as.numeric(Sys.getenv("DIFF", "25"))
q_cut    <- as.numeric(Sys.getenv("QVAL", "0.01"))
tile     <- as.integer(Sys.getenv("TILE", "1000"))
assembly <- Sys.getenv("ASSEMBLY", "genome")
cores    <- as.integer(Sys.getenv("CORES", "1"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(out_dir, x)
write_tsv <- function(x, f) write.table(x, out(f), sep = "\t", quote = FALSE, row.names = FALSE)
safely <- function(label, expr) tryCatch(expr, error = function(e) {
    message("Skipped ", label, ": ", conditionMessage(e)); NULL
})

ss <- ss[ss$condition %in% c(ref_level, tgt_level), ]
if (!all(c(ref_level, tgt_level) %in% ss$condition)) {
    stop("Samplesheet must contain both conditions: ", ref_level, " and ", tgt_level)
}
ss$treatment <- ifelse(ss$condition == tgt_level, 1, 0)
group_sizes  <- table(ss$condition)
message("Samples: ", paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "))

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
    safely("methylation stats", getMethylationStats(obj[[i]], plot = TRUE, both.strands = FALSE))
    safely("coverage stats",    getCoverageStats(obj[[i]], plot = TRUE, both.strands = FALSE))
}
invisible(dev.off())

# ---- Filter, normalise, merge ----------------------------------------------
# Very high coverage usually means PCR artefacts or repeats, so drop the top 0.1%.
# Then equalise sequencing depth so samples are comparable.
filt <- filterByCoverage(obj, lo.count = min_cov, hi.perc = 99.9)
filt <- normalizeCoverage(filt)
meth <- unite(filt, destrand = FALSE, mc.cores = cores)
message(nrow(meth), " CpGs covered in all samples")
if (nrow(meth) == 0) stop("No CpG is covered in every sample at MIN_COV=", min_cov, "; lower --min_cov")

# ---- Sample similarity -------------------------------------------------------
if (nrow(ss) >= 3) {
    pdf(out("sample_similarity.pdf"), width = 10, height = 10)
    safely("correlation", getCorrelation(meth, plot = TRUE))
    safely("clustering",  clusterSamples(meth, dist = "correlation", method = "ward.D", plot = TRUE))
    safely("PCA",         PCASamples(meth, screeplot = TRUE))
    safely("PCA",         PCASamples(meth))
    invisible(dev.off())
}

# ---- Differentially methylated CpGs (DMCs) -----------------------------------
# Biological replicates vary more than a simple binomial model expects (overdispersion);
# correcting for it avoids many false positives, but needs replicates in both groups.
od <- if (min(group_sizes) >= 2) "MN" else "none"
message("Overdispersion correction: ", od)
diff <- calculateDiffMeth(meth, overdispersion = od, test = if (od == "MN") "Chisq" else "F",
                          mc.cores = cores)

dmc       <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut)
dmc_hyper <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut, type = "hyper")
dmc_hypo  <- getMethylDiff(diff, difference = diff_cut, qvalue = q_cut, type = "hypo")

write_tsv(getData(diff),      "CpG_all_tested.tsv")
write_tsv(getData(dmc),       "DMC_all.tsv")
write_tsv(getData(dmc_hyper), "DMC_hyper.tsv")
write_tsv(getData(dmc_hypo),  "DMC_hypo.tsv")
message(nrow(dmc), " DMCs (|diff| >= ", diff_cut, "%, q < ", q_cut, ")")

if (nrow(dmc) > 0) {
    pdf(out("DMC_per_chromosome.pdf"), width = 10, height = 8)
    safely("per-chromosome plot", diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q_cut, meth.cutoff = diff_cut))
    invisible(dev.off())
}

# ---- Differentially methylated regions (tiles) -------------------------------
# Single CpGs are noisy; summing them in windows finds regions that change together.
safely("DMR tiling", {
    tiles  <- tileMethylCounts(filt, win.size = tile, step.size = tile, cov.bases = min_cov)
    meth_t <- unite(tiles, destrand = FALSE, mc.cores = cores)
    diff_t <- calculateDiffMeth(meth_t, overdispersion = od, test = if (od == "MN") "Chisq" else "F",
                                mc.cores = cores)
    dmr    <- getMethylDiff(diff_t, difference = diff_cut, qvalue = q_cut)
    write_tsv(getData(dmr), "DMR_tiles.tsv")
    message(nrow(dmr), " DMR tiles of ", tile, " bp")
})

# ---- Annotation (optional) ---------------------------------------------------
if (nrow(dmc) > 0 && (!is.na(gene_bed) || !is.na(cpg_bed))) {
    dmc_gr <- as(dmc, "GRanges")
    pdf(out("DMC_annotation.pdf"), width = 10, height = 10)
    if (!is.na(gene_bed)) safely("gene annotation", {
        genes     <- readTranscriptFeatures(gene_bed)
        ann_genes <- annotateWithGeneParts(dmc_gr, genes)
        write_tsv(getAssociationWithTSS(ann_genes), "DMC_nearest_TSS.tsv")
        plotTargetAnnotation(ann_genes, precedence = TRUE, main = "DMCs by gene part")
    })
    if (!is.na(cpg_bed)) safely("CpG-island annotation", {
        cpg     <- readFeatureFlank(cpg_bed, feature.flank.name = c("CpGi", "shores"))
        ann_cpg <- annotateWithFeatureFlank(dmc_gr, feature = cpg$CpGi, flank = cpg$shores,
                                            feature.name = "CpGi", flank.name = "shores")
        plotTargetAnnotation(ann_cpg, col = c("green", "gray", "white"), main = "DMCs by CpG-island context")
    })
    invisible(dev.off())
}

message("Done. Results written to ", normalizePath(out_dir))
