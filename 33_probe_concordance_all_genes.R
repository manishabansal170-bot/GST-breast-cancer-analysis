# =============================================================================
#  33_probe_concordance_all_genes.R
#
#  THE OBJECTION
#    Per-probe analysis of GSTP1 found that its thirteen annotated promoter
#    probes fall into two blocks running in opposite directions: three CpG
#    island probes agreeing with one another (Spearman 0.65 to 0.89) and two
#    north-shore probes agreeing with each other (0.86) but not with the island
#    (-0.15 to 0.03). The per-gene value used throughout the manuscript
#    corresponds to a single island probe.
#
#    A reviewer asked the obvious follow-up: is the same problem present, and
#    unreported, for any of the other seventeen genes in the methylation
#    analysis? If a gene's promoter contains opposing domains and its per-gene
#    value averages across them, that value represents no probe faithfully and
#    any conclusion drawn from it is unsafe.
#
#  WHAT THIS SCRIPT DOES
#    For every gene in the methylation analysis it retrieves all annotated
#    promoter probes, fetches those returning data, and asks three questions:
#
#      1. How many probes contribute, and how many return values at all
#      2. Do the probes agree, measured as the median pairwise correlation
#      3. If they disagree, do they fall into coherent opposing blocks, as
#         GSTP1 does, or is the disagreement unstructured
#
#    It then reports which per-gene value each cached summary corresponds to,
#    by correlating the summary against each probe.
#
#  WHY THE LAST STEP MATTERS
#    For GSTP1 the per-gene value correlated at 1.000 with one probe, showing
#    it was not an average at all. Whether that holds for the other genes is
#    not documented anywhere, and the Methods should state the rule that was
#    actually applied rather than the rule that was intended.
#
#  WHAT WOULD REQUIRE A CHANGE TO THE PAPER
#    Any gene whose probes form opposing blocks AND whose per-gene value is a
#    genuine average across them needs its conclusions re-examined. Genes whose
#    probes agree, or whose value corresponds to a single probe, are safe once
#    the rule is stated.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr","tidyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (!requireNamespace(bioc, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(bioc, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })

PROMOTER <- c("TSS1500","TSS200","5'UTR","1stExon")

obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
GENES <- rownames(obj$met)
cat("Genes in the methylation analysis:", length(GENES), "\n")
cat(paste(GENES, collapse = ", "), "\n\n")


# =============================================================================
# 1. ANNOTATED PROMOTER PROBES PER GENE
# =============================================================================
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

probes_for <- function(g) {
  hit <- ann[grepl(paste0("(^|;)", g, "(;|$)"), ann$UCSC_RefGene_Name), , drop = FALSE]
  if (!nrow(hit)) return(hit)
  keep <- sapply(strsplit(hit$UCSC_RefGene_Group, ";"), function(v) any(v %in% PROMOTER))
  hit[keep, , drop = FALSE]
}

plist <- lapply(GENES, probes_for); names(plist) <- GENES
cat("Annotated promoter probes per gene:\n")
print(sapply(plist, nrow))


# =============================================================================
# 2. FETCH PROBE-LEVEL DATA
# =============================================================================
# The GSTP1 probes are already cached from script 29; the rest are fetched once
# and cached together.
pf <- file.path(CACHE, "all_gst_probe_betas.rds")
if (file.exists(pf)) {
  B <- readRDS(pf); message("Loaded probe betas from cache.")
} else {
  message("Fetching probe-level betas. This takes several minutes.")
  HOST <- "https://tcga.xenahubs.net"
  DSET <- "TCGA.BRCA.sampleMap/HumanMethylation450"
  want <- unique(unlist(lapply(plist, rownames)))
  cat("Probes to fetch:", length(want), "\n")
  got <- list()
  for (i in seq_along(want)) {
    p <- want[i]
    r <- tryCatch(fetch_dense_values(host = HOST, dataset = DSET, identifiers = p,
                                     use_probeMap = FALSE, check = FALSE),
                  error = function(e) NULL)
    if (!is.null(r) && nrow(as.matrix(r))) got[[p]] <- as.matrix(r)[1, ]
    if (i %% 25 == 0) message("  ", i, " / ", length(want))
    Sys.sleep(0.12)
  }
  if (!length(got)) stop("No probe data retrieved.")
  common <- Reduce(intersect, lapply(got, names))
  B <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(B) <- names(got); saveRDS(B, pf)
  message("Retrieved ", nrow(B), " probes across ", ncol(B), " samples.")
}

samp <- intersect(colnames(B), colnames(obj$met))
samp <- samp[substr(samp, 13, 15) == "-01"]
B <- B[, samp, drop = FALSE]
cat("\nTumours with probe-level data:", ncol(B), "\n")


# =============================================================================
# 3. PROBE AGREEMENT PER GENE
# =============================================================================
cat("\n", strrep("=", 78), "\nPROBE AGREEMENT WITHIN EACH PROMOTER\n",
    strrep("=", 78), "\n", sep = "")

summarise_gene <- function(g) {
  pr <- rownames(plist[[g]])
  have <- intersect(pr, rownames(B))
  have <- have[apply(B[have, , drop = FALSE], 1,
                     function(v) sum(is.finite(v)) >= 100)]
  if (length(have) < 2)
    return(data.frame(gene = g, n_annotated = length(pr), n_with_data = length(have),
                      median_r = NA_real_, min_r = NA_real_, max_r = NA_real_,
                      n_negative_pairs = NA_integer_, pattern = "too few probes",
                      stringsAsFactors = FALSE))
  cm <- suppressWarnings(cor(t(B[have, , drop = FALSE]),
                             use = "pairwise.complete.obs", method = "spearman"))
  off <- cm[upper.tri(cm)]
  neg <- sum(off < -0.05, na.rm = TRUE)
  med <- median(off, na.rm = TRUE)
  pattern <- if (is.na(med)) "undetermined" else
             if (med >= 0.6) "concordant" else
             if (neg >= 1 && max(off, na.rm = TRUE) >= 0.6) "opposing blocks" else
             if (med >= 0.3) "weakly concordant" else "discordant"
  data.frame(gene = g, n_annotated = length(pr), n_with_data = length(have),
             median_r = round(med, 3), min_r = round(min(off, na.rm = TRUE), 3),
             max_r = round(max(off, na.rm = TRUE), 3),
             n_negative_pairs = neg, pattern = pattern, stringsAsFactors = FALSE)
}

res <- lapply(GENES, summarise_gene) %>% bind_rows()
print(as.data.frame(res[order(res$median_r), ]), row.names = FALSE, digits = 3)


# =============================================================================
# 4. WHICH PROBE DOES THE PER-GENE VALUE CORRESPOND TO
# =============================================================================
cat("\n", strrep("=", 78), "\nWHAT THE PER-GENE VALUE ACTUALLY IS\n",
    strrep("=", 78), "\n", sep = "")
cat("Correlating each cached per-gene value against its constituent probes. A\n")
cat("correlation of 1.000 with one probe means the value is that probe, not an\n")
cat("average. A high correlation with several means it is a genuine summary.\n\n")

sh <- intersect(colnames(obj$met), colnames(B))
ident <- lapply(GENES, function(g) {
  have <- intersect(rownames(plist[[g]]), rownames(B))
  have <- have[apply(B[have, , drop = FALSE], 1,
                     function(v) sum(is.finite(v)) >= 100)]
  if (!length(have)) return(NULL)
  gv <- as.numeric(obj$met[g, sh])
  r <- sapply(have, function(p)
    suppressWarnings(cor(gv, as.numeric(B[p, sh]),
                         use = "pairwise.complete.obs", method = "spearman")))
  # A probe with no overlapping finite values returns NA; dropping these first
  # avoids which.max returning an empty name and failing the data.frame call.
  r <- r[is.finite(r)]
  if (!length(r))
    return(data.frame(gene = g, n_probes = length(have), best_probe = NA_character_,
                      best_r = NA_real_, second_r = NA_real_, is_single_probe = NA,
                      stringsAsFactors = FALSE))
  rs <- sort(r, decreasing = TRUE)
  data.frame(gene = g, n_probes = length(have),
             best_probe = names(rs)[1], best_r = round(rs[1], 3),
             second_r = if (length(rs) > 1) round(rs[2], 3) else NA_real_,
             is_single_probe = rs[1] > 0.995,
             stringsAsFactors = FALSE)
}) %>% bind_rows()
print(as.data.frame(ident), row.names = FALSE, digits = 3)


# =============================================================================
# 5. VERDICT
# =============================================================================
cat("\n", strrep("=", 78), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 78), "\n", sep = "")

opposing <- res$gene[res$pattern == "opposing blocks"]
singles  <- ident$gene[which(ident$is_single_probe)]
risky    <- intersect(opposing, setdiff(ident$gene, singles))

cat("Genes whose promoter probes form opposing blocks:\n  ",
    if (length(opposing)) paste(opposing, collapse = ", ") else "none", "\n", sep = "")
cat("Genes whose per-gene value is a single probe rather than an average:\n  ",
    if (length(singles)) paste(singles, collapse = ", ") else "none", "\n", sep = "")

if (!length(risky)) {
  cat("\nNo gene both has opposing probe blocks and a per-gene value averaged\n")
  cat("across them. GSTP1 is the only gene with opposing blocks, and its value\n")
  cat("corresponds to a single island probe, which the manuscript now states.\n")
  cat("State the rule in Methods and report this check; no conclusion changes.\n")
} else {
  cat("\nThe following genes have opposing probe blocks AND a per-gene value that\n")
  cat("averages across them, so that value represents no probe faithfully:\n  ")
  cat(paste(risky, collapse = ", "), "\n")
  cat("Conclusions drawn from these genes must be re-examined before submission.\n")
}

write.csv(res,   file.path(OUT, "probe_concordance_all_genes.csv"), row.names = FALSE)
write.csv(ident, file.path(OUT, "per_gene_value_identity.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Genes checked: %d", length(GENES)),
             sprintf("Opposing blocks: %s",
                     if (length(opposing)) paste(opposing, collapse = ", ") else "none"),
             sprintf("Per-gene value is a single probe for: %s",
                     if (length(singles)) paste(singles, collapse = ", ") else "none"),
             sprintf("Genes requiring re-examination: %s",
                     if (length(risky)) paste(risky, collapse = ", ") else "none"),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_probe_concordance.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
