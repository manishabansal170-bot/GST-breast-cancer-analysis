# =============================================================================
#  32_variance_partition_within_subtype.R
#
#  THE OBJECTION, AND WHY IT IS THE MOST SERIOUS ONE RAISED
#    The manuscript's central claim is that promoter methylation, rather than
#    copy number, marks the boundary between the two GST programmes. The
#    evidence is a joint regression of rank expression on rank copy number and
#    rank methylation, pooled across 770 tumours of all subtypes.
#
#    Both methylation and expression vary strongly by subtype. A gene that is
#    highly expressed and unmethylated in basal-like disease, and silenced and
#    methylated in luminal B, will produce a large unique methylation R-squared
#    from that contrast alone, even if methylation explains nothing about the
#    variation between tumours WITHIN either subtype.
#
#    This is precisely the confound the manuscript identifies and dismantles
#    for the co-expression analysis, where GSTP1 against EGFR went from +0.475
#    pooled to -0.095 within basal-like, and for the microRNA analysis. Applying
#    that standard to those analyses but not to the variance partition is
#    inconsistent, and the reviewer is right to say so.
#
#  WHAT THIS SCRIPT DOES
#    1. Repeats the variance partition within each PAM50 subtype separately
#    2. Repeats it pooled, with subtype as a covariate
#    3. Repeats it with sample-level mean methylation as a covariate, which
#       addresses the separate objection that luminal B tumours are enriched
#       for a global hypermethylator phenotype
#    4. States plainly whether purity adjustment was applied
#
#  HOW TO READ THE RESULT
#    If unique methylation R-squared remains substantial within basal-like
#    tumours alone, the claim is secure and considerably stronger than at
#    present: methylation explains variation between tumours of the same
#    subtype, which composition cannot produce.
#
#    If it collapses the way the co-expression correlations collapsed, the
#    claim must be restated. Methylation would then be a marker of the
#    programme division rather than a mechanism operating within it. That is a
#    real finding and still publishable, but it is a different sentence and the
#    manuscript must say so.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleBG_alterations_protein")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")

# Genes carried forward for inference: the 18 in the model, less the two
# confounded (GSTM1 germline deletion, GSTT2 paralogue ambiguity) and the two
# below detection (GSTA2, GSTT1).
EXCLUDE <- c("GSTM1","GSTT2","GSTA2","GSTT1")


# =============================================================================
# 1. INPUTS
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; expr <- obj$exp
# Linear copy number, not the discrete GISTIC calls: the published analysis
# ranks its predictors, and discrete calls take too few distinct values to
# rank informatively.
cnv <- NULL
for (nm in c("tcga_cna_linear_gst.rds", "tcga_cna_gst.rds", "tcga_cna_gistic_gst.rds")) {
  fp <- file.path(CACHE, nm)
  if (file.exists(fp)) { cnv <- readRDS(fp); message("Copy number from: ", nm); break }
}
if (is.null(cnv)) {
  cand <- list.files(CACHE, pattern = "cna|cnv|copy", full.names = FALSE)
  cat("No copy number object found. Files matching in cache:\n"); print(cand)
  stop("Set the copy number path manually and rerun.")
}
if (is.list(cnv) && !is.matrix(cnv)) cnv <- cnv[[which(vapply(cnv, is.matrix, logical(1)))[1]]]

st <- TCGAquery_subtype(tumor = "brca")

s <- Reduce(intersect, list(colnames(met), colnames(expr), colnames(cnv)))
s <- s[substr(s, 13, 15) == "-01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(s, 1, 12), st$patient)]
ok <- grp %in% SUB
s <- s[ok]; grp <- factor(grp[ok], levels = SUB)

genes <- Reduce(intersect, list(rownames(met), rownames(expr), rownames(cnv)))
cat("Tumours with all three assays and a PAM50 call:", length(s), "\n")
print(table(grp))
cat("Genes in all three matrices:", length(genes), "\n")

# purity, to state explicitly whether the model is adjusted
pur <- rep(NA_real_, length(s))
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  pur <- num(tp$CPE)[match(substr(s, 1, 12), tp$patient)]
}, silent = TRUE)
cat("Tumours with a purity estimate:", sum(is.finite(pur)), "\n")

# global methylation, for the hypermethylator objection
global_meth <- colMeans(met[, s, drop = FALSE], na.rm = TRUE)
cat(sprintf("Sample-level mean beta: median %.3f, range %.3f to %.3f\n",
            median(global_meth, na.rm = TRUE), min(global_meth, na.rm = TRUE),
            max(global_meth, na.rm = TRUE)))
cat("Mean beta by subtype:\n"); print(round(tapply(global_meth, grp, mean, na.rm = TRUE), 3))


# =============================================================================
# 2. THE PARTITION FUNCTION
# =============================================================================
# Unique contribution of each predictor is the drop in R-squared when it is
# removed from the full model. Ranks are used throughout, matching the
# published analysis.
partition <- function(g, idx, covars = NULL, label = "") {
  e <- rank(as.numeric(expr[g, s[idx]]))
  m <- rank(as.numeric(met [g, s[idx]]))
  c_ <- rank(as.numeric(cnv[g, s[idx]]))
  df <- data.frame(e = e, m = m, c = c_)
  if (!is.null(covars)) for (nm in names(covars)) df[[nm]] <- covars[[nm]][idx]
  df <- df[complete.cases(df), ]
  if (nrow(df) < 40) return(NULL)

  extra <- setdiff(names(df), c("e","m","c"))
  rhs <- function(v) paste(c(v, extra), collapse = " + ")
  r2 <- function(f) summary(lm(as.formula(f), data = df))$r.squared

  full  <- r2(paste("e ~", rhs(c("m","c"))))
  no_m  <- r2(paste("e ~", rhs("c")))
  no_c  <- r2(paste("e ~", rhs("m")))
  base  <- if (length(extra)) r2(paste("e ~", paste(extra, collapse = " + "))) else 0

  fit <- lm(as.formula(paste("e ~", rhs(c("m","c")))), data = df)
  pm <- tryCatch(summary(fit)$coefficients["m", 4], error = function(e) NA_real_)

  data.frame(stratum = label, gene = g, n = nrow(df),
             R2_total  = round(full - base, 3),
             unique_meth = round(full - no_m, 3),
             unique_cn   = round(full - no_c, 3),
             p_meth = pm,
             stringsAsFactors = FALSE)
}


# =============================================================================
# 3. POOLED, AS PUBLISHED
# =============================================================================
cat("\n", strrep("=", 74), "\nA. POOLED, AS PUBLISHED\n", strrep("=", 74), "\n", sep = "")
all_idx <- seq_along(s)
pooled <- lapply(genes, function(g) partition(g, all_idx, NULL, "pooled")) %>% bind_rows()
pooled$FDR <- p.adjust(pooled$p_meth, "BH")
pooled$interpretable <- !(pooled$gene %in% EXCLUDE)
print(as.data.frame(pooled[order(-pooled$unique_meth),
      c("gene","n","R2_total","unique_meth","unique_cn","FDR")]), row.names = FALSE, digits = 3)


# =============================================================================
# 4. WITHIN SUBTYPE
# =============================================================================
cat("\n", strrep("=", 74), "\nB. WITHIN EACH SUBTYPE\n", strrep("=", 74), "\n", sep = "")
within <- lapply(SUB, function(k) {
  idx <- which(grp == k)
  lapply(genes, function(g) partition(g, idx, NULL, k)) %>% bind_rows()
}) %>% bind_rows()
if (nrow(within)) within <- within %>% group_by(stratum) %>%
  mutate(FDR = p.adjust(p_meth, "BH")) %>% ungroup()

focus <- c("GSTP1","GSTM3","GSTM5","GSTM2","GSTM4","GSTO2","GSTA4","GSTA1")
cat("Unique methylation R-squared, by subtype, for the genes carrying the model:\n\n")
w <- within %>% dplyr::filter(gene %in% focus) %>%
  dplyr::select(gene, stratum, unique_meth) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = unique_meth)
p_ <- pooled %>% dplyr::filter(gene %in% focus) %>% dplyr::select(gene, pooled = unique_meth)
print(as.data.frame(dplyr::left_join(p_, w, by = "gene")), row.names = FALSE, digits = 3)


# =============================================================================
# 5. POOLED WITH SUBTYPE AS A COVARIATE
# =============================================================================
cat("\n", strrep("=", 74), "\nC. POOLED, WITH SUBTYPE AS A COVARIATE\n", strrep("=", 74), "\n", sep = "")
adj_sub <- lapply(genes, function(g)
  partition(g, all_idx, list(subtype = grp), "subtype-adjusted")) %>% bind_rows()
adj_sub$FDR <- p.adjust(adj_sub$p_meth, "BH")
cmp <- dplyr::left_join(
  pooled %>% dplyr::select(gene, pooled = unique_meth),
  adj_sub %>% dplyr::select(gene, subtype_adjusted = unique_meth), by = "gene")
cmp$change <- round(cmp$subtype_adjusted - cmp$pooled, 3)
print(as.data.frame(cmp[order(-cmp$pooled), ]), row.names = FALSE, digits = 3)


# =============================================================================
# 6. GLOBAL METHYLATION AS A COVARIATE
# =============================================================================
cat("\n", strrep("=", 74), "\nD. WITH SAMPLE-LEVEL MEAN METHYLATION AS A COVARIATE\n",
    strrep("=", 74), "\n", sep = "")
cat("Addresses the objection that luminal B tumours are enriched for a global\n")
cat("hypermethylator phenotype, so that per-gene methylation may be a proxy for\n")
cat("a genome-wide state rather than locus-specific regulation.\n\n")
adj_glob <- lapply(genes, function(g)
  partition(g, all_idx, list(gmeth = rank(global_meth), subtype = grp),
            "subtype and global methylation adjusted")) %>% bind_rows()
adj_glob$FDR <- p.adjust(adj_glob$p_meth, "BH")
cmp2 <- dplyr::left_join(cmp,
  adj_glob %>% dplyr::select(gene, plus_global = unique_meth), by = "gene")
print(as.data.frame(cmp2[order(-cmp2$pooled), ]), row.names = FALSE, digits = 3)


# =============================================================================
# 7. VERDICT
# =============================================================================
cat("\n", strrep("=", 74), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 74), "\n", sep = "")

bas <- within %>% dplyr::filter(stratum == "Basal", gene %in% focus)
kept <- sum(bas$unique_meth >= 0.10, na.rm = TRUE)
cat(sprintf("Within basal-like tumours alone, %d of %d programme genes retain unique\n",
            kept, nrow(bas)))
cat("methylation R-squared of at least 0.10.\n\n")

gp <- within$unique_meth[within$gene == "GSTP1" & within$stratum == "Basal"]
gpp <- pooled$unique_meth[pooled$gene == "GSTP1"]
if (length(gp) && length(gpp))
  cat(sprintf("GSTP1: pooled %.3f, within basal-like %.3f\n", gpp, gp))

if (kept >= 0.6 * nrow(bas)) {
  cat("\nMethylation retains a substantial unique contribution within subtype.\n")
  cat("The mechanistic claim is secure and is stronger than the pooled analysis\n")
  cat("alone showed, because within-subtype variation cannot be produced by\n")
  cat("composition. Report the within-subtype values alongside the pooled ones.\n")
} else if (kept > 0) {
  cat("\nMethylation retains a contribution for some genes but not most. Restrict\n")
  cat("the mechanistic claim to the genes where it survives, name them, and\n")
  cat("describe methylation as marking the division for the remainder.\n")
} else {
  cat("\nThe unique methylation contribution does not survive within subtype. The\n")
  cat("pooled result reflects the between-subtype contrast. The claim must be\n")
  cat("restated: methylation marks the programme division rather than operating\n")
  cat("as a mechanism within it. This is still a finding, and reporting it\n")
  cat("honestly is consistent with how the co-expression and microRNA analyses\n")
  cat("are already handled in this manuscript.\n")
}

out <- bind_rows(pooled, within, adj_sub, adj_glob)
write.csv(out, file.path(OUT, "variance_partition_within_subtype.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d", length(s)),
             sprintf("Purity estimates available for %d tumours; the published joint model was NOT purity-adjusted and this is now stated in the manuscript", sum(is.finite(pur))),
             sprintf("Programme genes retaining unique methylation R2 >= 0.10 within basal-like: %d of %d", kept, nrow(bas)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_variance_within.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
