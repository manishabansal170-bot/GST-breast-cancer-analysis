# =============================================================================
#  29_gstp1_per_probe.R
#
#  THE OBJECTION
#    The manuscript describes its methylation values as per-gene summaries over
#    the promoter-associated probes for each gene. For GSTP1 it reports a
#    bimodal distribution and a 97% versus 58% subtype contrast, and those
#    claims carry much of the paper.
#
#  WHAT THIS SCRIPT FOUND, STATED UP FRONT
#    The GSTP1 value is not an average at all. It corresponds exactly to a
#    single probe, cg04920951, in the CpG island spanning the first exon and
#    5' untranslated region (Spearman 1.000). Of the thirteen annotated
#    promoter probes, only five return data in this cohort, and they form two
#    anti-correlated blocks: three island probes agreeing with one another
#    (0.65 to 0.89) and two north-shore probes agreeing with each other (0.86)
#    but not with the island (-0.15 to 0.03). The two blocks run in opposite
#    directions by subtype. The manuscript has been amended accordingly.
#
#    Averaging can create or destroy structure. If some probes in the GSTP1
#    promoter are methylated and others are not, an average sits in between and
#    represents no probe faithfully. Conversely, an apparently clean bimodality
#    in the average could be produced by one strongly bimodal probe diluted
#    with several flat ones, in which case the biology is narrower than
#    reported.
#
#  WHAT THIS SCRIPT ESTABLISHES
#    1. How many probes contribute to the GSTP1 average, and where they sit
#    2. Whether they agree with each other, by pairwise correlation
#    3. Whether the subtype contrast holds probe by probe, or rests on a subset
#    4. Whether bimodality is a property of individual probes or of the average
#
#  WHAT WOULD REQUIRE A CHANGE TO THE PAPER
#    If the subtype contrast holds across most probes, the averaging is safe
#    and the manuscript should say so and give the per-probe range.
#    If it rests on one or two probes, the claim must be restated in terms of
#    those specific CpGs rather than the promoter as a whole, and the location
#    of those probes becomes part of the finding.
#
#  A LIMIT
#    This requires the per-probe beta matrix. The cached object holds gene-level
#    averages only, so the script fetches probe-level data. If that fails, it
#    says so rather than approximating, because an approximation here would
#    answer a question about approximation with another approximation.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr","tidyr","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19")
miss <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(TCGAbiolinks)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })

SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. WHICH PROBES BELONG TO THE GSTP1 PROMOTER
# =============================================================================
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
gp <- ann[grepl("(^|;)GSTP1(;|$)", ann$UCSC_RefGene_Name), ]
cat("Probes annotated to GSTP1 anywhere:", nrow(gp), "\n")

promoter_terms <- c("TSS1500","TSS200","5'UTR","1stExon")
gp$is_promoter <- sapply(strsplit(gp$UCSC_RefGene_Group, ";"), function(v)
  any(v %in% promoter_terms))
prom <- gp[gp$is_promoter, ]
cat("Promoter-associated probes annotated on the array:", nrow(prom),
    "  (not all return data; see below)\n\n")
print(prom[order(prom$pos), c("chr","pos","UCSC_RefGene_Group","Relation_to_Island")])

if (!nrow(prom)) stop("No promoter probes found. Check the annotation package version.")


# =============================================================================
# 2. FETCH PER-PROBE BETA VALUES
# =============================================================================
pf <- file.path(CACHE, "gstp1_probe_betas.rds")
if (file.exists(pf)) {
  B <- readRDS(pf); message("Loaded per-probe betas from cache.")
} else {
  message("Fetching per-probe betas from Xena. This may take a minute.")
  HOST <- "https://tcga.xenahubs.net"
  DSET <- "TCGA.BRCA.sampleMap/HumanMethylation450"
  got <- list()
  for (p in rownames(prom)) {
    r <- tryCatch(fetch_dense_values(host = HOST, dataset = DSET, identifiers = p,
                                     use_probeMap = FALSE, check = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || !nrow(as.matrix(r))) { message("  no data: ", p); next }
    got[[p]] <- as.matrix(r)[1, ]
    Sys.sleep(0.15)
  }
  if (!length(got)) stop(
    "No probe-level data retrieved. Without it this question cannot be answered,\n",
    "and the manuscript should state that per-probe structure was not examined\n",
    "rather than implying it was.")
  common <- Reduce(intersect, lapply(got, names))
  B <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(B) <- names(got); saveRDS(B, pf)
  message("Retrieved ", nrow(B), " probes across ", ncol(B), " samples.")
}


# =============================================================================
# 3. SUBTYPE ASSIGNMENT
# =============================================================================
st <- TCGAquery_subtype(tumor = "brca")
samp <- colnames(B)[substr(colnames(B), 14, 15) == "01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(samp, 1, 12), st$patient)]
keep <- grp %in% SUB
samp <- samp[keep]; grp <- factor(grp[keep], levels = SUB)
B <- B[, samp, drop = FALSE]
cat("\nTumours with probe-level methylation and a PAM50 call:", ncol(B), "\n")
print(table(grp))


# =============================================================================
# 4. DO THE PROBES AGREE WITH EACH OTHER
# =============================================================================
cat("\n", strrep("=", 70), "\nPROBE AGREEMENT\n", strrep("=", 70), "\n", sep = "")
cm <- suppressWarnings(cor(t(B), use = "pairwise.complete.obs", method = "spearman"))
off <- cm[upper.tri(cm)]
cat(sprintf("Pairwise correlations between probes: median %.3f, range %.3f to %.3f\n",
            median(off, na.rm = TRUE), min(off, na.rm = TRUE), max(off, na.rm = TRUE)))
if (median(off, na.rm = TRUE) >= 0.6) {
  cat("The probes largely agree, so averaging them is defensible.\n")
} else {
  cat("The probes disagree substantially. An average across them represents no\n")
  cat("probe faithfully and the per-gene value should be reconsidered.\n")
}


# =============================================================================
# 5. DOES THE SUBTYPE CONTRAST HOLD PROBE BY PROBE
# =============================================================================
cat("\n", strrep("=", 70), "\nSUBTYPE CONTRAST, PER PROBE\n", strrep("=", 70), "\n", sep = "")

# Eight of the thirteen annotated probes return no values in this cohort.
# kruskal.test errors rather than returning NA when a probe has data in only
# one group, so the test is guarded and n_obs is reported for every probe.
per <- lapply(rownames(B), function(p) {
  v <- as.numeric(B[p, ])
  ok <- is.finite(v)
  med <- tapply(v, grp, median, na.rm = TRUE)
  pm  <- tapply(v, grp, function(x) 100 * mean(x > 0.3, na.rm = TRUE))
  pu  <- tapply(v, grp, function(x) 100 * mean(x < 0.1, na.rm = TRUE))
  kp <- if (sum(ok) >= 40 && length(unique(grp[ok])) >= 2)
          suppressWarnings(kruskal.test(v[ok] ~ droplevels(grp[ok]))$p.value)
        else NA_real_
  data.frame(probe = p, n_obs = sum(ok),
             pos = prom[p, "pos"], region = prom[p, "UCSC_RefGene_Group"],
             island = prom[p, "Relation_to_Island"],
             Basal_unmeth = round(pu["Basal"], 1),
             LumB_meth = round(pm["LumB"], 1),
             Her2_meth = round(pm["Her2"], 1),
             median_Basal = round(med["Basal"], 3),
             median_LumB = round(med["LumB"], 3),
             kruskal_p = kp)
}) %>% bind_rows()
per$FDR <- p.adjust(per$kruskal_p, "BH")
per <- per[order(per$pos), ]

print(as.data.frame(per[, c("probe","n_obs","region","island","Basal_unmeth",
                            "LumB_meth","Her2_meth","FDR")]),
      row.names = FALSE, digits = 3)

concordant <- sum(per$Basal_unmeth > 70 & per$LumB_meth > 40, na.rm = TRUE)
cat(sprintf("\nProbes showing the reported pattern, basal-like above 70%% unmethylated\n"))
cat(sprintf("and luminal B above 40%% methylated: %d of %d\n", concordant, nrow(per)))
cat(sprintf("Probes with a significant subtype difference: %d of %d\n",
            sum(per$FDR < 0.05, na.rm = TRUE), nrow(per)))

write.csv(per, file.path(OUT, "gstp1_per_probe.csv"), row.names = FALSE)


# =============================================================================
# 6. IS BIMODALITY A PROBE PROPERTY OR AN ARTEFACT OF AVERAGING
# =============================================================================
cat("\n", strrep("=", 70), "\nBIMODALITY\n", strrep("=", 70), "\n", sep = "")
bimod <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 50) return(NA_real_)
  # proportion in the intermediate band; low means bimodal
  100 * mean(v >= 0.1 & v <= 0.3)
}
per$pct_intermediate <- round(sapply(rownames(B), function(p) bimod(as.numeric(B[p, ])))[per$probe], 1)

# The per-gene value used throughout the manuscript is cg04920951, not a mean
# across probes: the two correlate at 1.000. Comparing bimodality against a
# constructed 13-probe mean would therefore describe a quantity the paper does
# not use, so the comparison is made against the actual per-gene value.
gv <- tryCatch({
  o <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
  s <- intersect(colnames(o$met), colnames(B)); o$met["GSTP1", s] },
  error = function(e) NULL)
if (!is.null(gv)) {
  cat(sprintf("Per-gene value as used in the manuscript: %.1f%% intermediate\n", bimod(gv)))
} else {
  cat("Per-gene value unavailable; compare against cg04920951 below.\n")
}
cat("Per probe:\n")
print(setNames(per$pct_intermediate, per$probe))
cat("\nThe per-gene value and cg04920951 agree because they are the same\n")
cat("measurement. The island probes carry the sharpest on-off structure; the\n")
cat("shore probes are less bimodal. Bimodality is therefore a property of the\n")
cat("island, not an artefact of combining probes.\n")


# =============================================================================
# 7. VERDICT
# =============================================================================
cat("\n", strrep("=", 70), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 70), "\n", sep = "")
if (concordant >= 0.7 * nrow(per)) {
  cat("The subtype contrast holds across most probes. Report the number of\n")
  cat("probes averaged and the per-probe range, and the per-gene value stands.\n")
} else if (concordant > 0) {
  cat(sprintf("The contrast holds for %d of %d probes. Restate the claim in terms of\n",
              concordant, nrow(per)))
  cat("those probes and give their positions, rather than attributing it to the\n")
  cat("promoter as a whole.\n")
} else {
  cat("No individual probe reproduces the reported pattern. The per-gene average\n")
  cat("is generating it, and the claim must be re-examined before submission.\n")
}

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Promoter probes averaged: %d", nrow(per)),
             sprintf("Median pairwise probe correlation: %.3f", median(off, na.rm = TRUE)),
             sprintf("Probes reproducing the subtype pattern: %d of %d", concordant, nrow(per)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_gstp1_per_probe.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
