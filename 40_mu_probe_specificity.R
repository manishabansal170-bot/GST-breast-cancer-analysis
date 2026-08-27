# =============================================================================
#  40_mu_probe_specificity.R
#
#  THE OBJECTION, AND WHY IT IS THE MOST CONSEQUENTIAL ONE LEFT
#    The only clinical association in this manuscript is that GSTM1 and GSTM2
#    expression predicts reduced pathological complete response in
#    receptor-positive disease in GSE25066. Two facts about that result invite
#    a specific doubt.
#
#    First, GSTM1, GSTM2, GSTM4 and GSTM5 are tandem paralogues on chromosome
#    1p13.3 with high sequence identity, and probe cross-hybridisation within
#    this family on Affymetrix arrays is documented.
#
#    Second, the two hits are precisely a paralogue pair, and their effect
#    estimates are nearly identical: odds ratios 0.483 and 0.508.
#
#    If the two probe sets are measuring the same transcripts, the manuscript
#    reports one finding as two. The reviewer is right that this must be
#    excluded rather than assumed away.
#
#  WHAT THIS SCRIPT ESTABLISHES
#    1. Which probe set identifiers are used for each Mu gene, in each cohort
#    2. How strongly those probe values correlate with one another
#    3. Whether the same probe set exists on both platforms, which bears on the
#       separate question of whether the failed replication reflects a platform
#       difference rather than sample size
#
#  HOW TO READ THE RESULT
#    A correlation above roughly 0.8 between the GSTM1 and GSTM2 probe values
#    means the two associations should be reported as one measurement, and the
#    manuscript must say so.
#
#    A correlation comparable to that between GSTM1 and non-adjacent family
#    members, or to the GSTM1-GSTM2 correlation in RNA-seq data where
#    cross-hybridisation cannot occur, means the probes are behaving as
#    intended and the two associations are separable.
#
#    The RNA-seq comparison is the important control: it gives the correlation
#    these genes genuinely have, against which the array correlation can be
#    judged.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleE_drug_response")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- c("GEOquery","hgu133a.db","hgu133plus2.db","AnnotationDbi")
miss <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(GEOquery); library(AnnotationDbi)
  library(hgu133a.db); library(hgu133plus2.db) })

MU <- c("GSTM1","GSTM2","GSTM3","GSTM4","GSTM5")
OTHER <- c("GSTP1","GSTZ1","MGST1")     # non-Mu comparators


# =============================================================================
# 1. WHICH PROBE SETS MAP TO EACH MU GENE
# =============================================================================
probes_for <- function(db, genes) {
  m <- AnnotationDbi::select(db, keys = genes, columns = "PROBEID", keytype = "SYMBOL")
  m[!is.na(m$PROBEID), ]
}

pa <- probes_for(hgu133a.db, c(MU, OTHER))
pp <- probes_for(hgu133plus2.db, c(MU, OTHER))

cat(strrep("=", 76), "\nPROBE SETS PER GENE\n", strrep("=", 76), "\n", sep = "")
cat("\nHG-U133A, used by GSE25066:\n")
print(pa[pa$SYMBOL %in% MU, ], row.names = FALSE)
cat("\nHG-U133 Plus 2, used by GSE32646:\n")
print(pp[pp$SYMBOL %in% MU, ], row.names = FALSE)

shared <- intersect(pa$PROBEID[pa$SYMBOL %in% MU], pp$PROBEID[pp$SYMBOL %in% MU])
cat("\nProbe sets present on BOTH platforms for Mu genes:", length(shared), "\n")
if (length(shared)) print(pa[pa$PROBEID %in% shared, ], row.names = FALSE)
cat("\nThis bears on the replication question: if the same probe set is on both\n")
cat("arrays, the failure to replicate cannot be attributed to platform.\n")


# =============================================================================
# 2. HOW STRONGLY DO THE PROBE VALUES CORRELATE, IN GSE25066
# =============================================================================
f1 <- list.files(CACHE, pattern = "GSE25066.*series_matrix", full.names = TRUE)[1]
if (is.na(f1)) stop("GSE25066 series matrix not in cache.")
e1 <- getGEO(filename = f1, getGPL = FALSE)
ex1 <- Biobase::exprs(e1)
if (max(ex1, na.rm = TRUE) > 100) ex1 <- log2(ex1 + 1)
cat("\nGSE25066:", nrow(ex1), "probe sets x", ncol(ex1), "samples\n")

# the probe set actually used per gene: highest variance, as in script 35
pick <- function(ex, map, genes) {
  out <- character(0)
  for (g in genes) {
    ids <- intersect(map$PROBEID[map$SYMBOL == g], rownames(ex))
    if (!length(ids)) next
    v <- apply(ex[ids, , drop = FALSE], 1, var, na.rm = TRUE)
    out[g] <- ids[which.max(v)]
  }
  out
}
sel <- pick(ex1, pa, c(MU, OTHER))
cat("\nProbe set selected per gene (highest variance):\n")
print(data.frame(gene = names(sel), probe = unname(sel)), row.names = FALSE)

cat("\n", strrep("=", 76), "\nCORRELATION BETWEEN MU PROBE VALUES, GSE25066 ARRAY\n",
    strrep("=", 76), "\n", sep = "")
M <- ex1[sel, , drop = FALSE]; rownames(M) <- names(sel)
cm <- round(cor(t(M), method = "spearman", use = "pairwise.complete.obs"), 3)
print(cm)

key <- cm["GSTM1","GSTM2"]
cat(sprintf("\nGSTM1 against GSTM2 on the array: rho = %.3f\n", key))


# =============================================================================
# 3. THE CONTROL: THE SAME PAIR IN RNA-SEQ, WHERE CROSS-HYBRIDISATION
#    CANNOT OCCUR
# =============================================================================
cat("\n", strrep("=", 76), "\nTHE SAME PAIR IN RNA-SEQ\n", strrep("=", 76), "\n", sep = "")
cat("Short-read RNA-seq assigns reads by sequence, so a high correlation here\n")
cat("reflects genuine co-regulation rather than probe cross-hybridisation.\n\n")

gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^gm(raw) - 0.001; lin[lin < 0] <- 0
E <- log2(lin + 1)
ts <- grep("^TCGA-.*-01$", colnames(E), value = TRUE)
have <- intersect(c(MU, OTHER), rownames(E))
cs <- round(cor(t(E[have, ts, drop = FALSE]), method = "spearman",
                use = "pairwise.complete.obs"), 3)
print(cs)
seq_key <- cs["GSTM1","GSTM2"]
cat(sprintf("\nGSTM1 against GSTM2 in RNA-seq: rho = %.3f\n", seq_key))


# =============================================================================
# 4. VERDICT
# =============================================================================
cat("\n", strrep("=", 76), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 76), "\n", sep = "")
cat(sprintf("Array rho %.3f against RNA-seq rho %.3f, difference %+.3f\n",
            key, seq_key, key - seq_key))

if (key > 0.8 && key - seq_key > 0.2) {
  cat("\nThe two probe sets correlate far more strongly on the array than the\n")
  cat("genes do in sequence-based data. This is the signature of cross-\n")
  cat("hybridisation: the array is not separating the paralogues. The GSTM1 and\n")
  cat("GSTM2 associations should be reported as a single Mu-class measurement,\n")
  cat("not as two findings.\n")
} else if (key > 0.8) {
  cat("\nThe probe values correlate strongly, but so do the genes themselves in\n")
  cat("RNA-seq, so the array correlation is consistent with genuine co-regulation\n")
  cat("rather than cross-hybridisation. Report both values and state this\n")
  cat("comparison, since a reader cannot otherwise distinguish the two causes.\n")
} else {
  cat("\nThe probe values do not correlate strongly enough to suggest they are the\n")
  cat("same measurement. The two associations are separable. Report the\n")
  cat("correlation and the probe set identifiers so the reader need not assume it.\n")
}

write.csv(data.frame(gene = names(sel), probe_set = unname(sel),
                     platform = "HG-U133A"),
          file.path(OUT, "mu_probe_sets_GSE25066.csv"), row.names = FALSE)
write.csv(as.data.frame(cm), file.path(OUT, "mu_probe_correlation_array.csv"))
write.csv(as.data.frame(cs), file.path(OUT, "mu_gene_correlation_rnaseq.csv"))

writeLines(c(paste("Run:", Sys.time()),
             sprintf("GSTM1-GSTM2 array rho: %.3f", key),
             sprintf("GSTM1-GSTM2 RNA-seq rho: %.3f", seq_key),
             sprintf("Mu probe sets shared between U133A and U133 Plus 2: %d", length(shared)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_mu_probes.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
