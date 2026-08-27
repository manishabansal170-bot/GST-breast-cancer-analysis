# =============================================================================
#  34_matrix_comparison_and_fair_partition.R
#
#  TWO OBJECTIONS ABOUT THE SAME TABLE
#
#  POINT 2. TWO EXPRESSION MATRICES CARRY THE PAPER
#    Table 2, Figures 1 to 3 and the whole expression architecture use the
#    harmonised GRCh38 Toil recompute. The methylation correlations, the joint
#    variance partition and Figure 4 use a differently annotated TCGA matrix
#    that receives none of the annotation care described in Methods. Table 3's
#    own legend warns that its two halves are not paired measurements.
#
#    A reader cannot currently judge how much this matters. This script reports
#    the per-gene correlation between the two matrices across shared samples,
#    and reruns the joint partition on the harmonised matrix so the two can be
#    compared directly.
#
#  POINT 5. THE CONTEST BETWEEN METHYLATION AND COPY NUMBER IS NOT FAIR
#    450k beta is continuous, wide-ranging and strongly bimodal at these loci.
#    Segment-level copy number is coarse and near zero for most tumours at most
#    of these genes. Comparing unique R-squared between predictors of such
#    different resolution favours the finer one regardless of biology.
#
#    This script therefore reports, for every gene: the marginal R-squared for
#    each predictor fitted alone, the shared variance component, both unique
#    terms, and a sensitivity analysis using discrete GISTIC calls in place of
#    linear copy number. If the ranking survives all three, the word "dominant"
#    is defensible; if it does not, "explains substantially more variance in
#    this model" is the honest phrasing.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleBG_alterations_protein")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")
EXCLUDE <- c("GSTM1","GSTT2","GSTA2","GSTT1")


# =============================================================================
# 1. THE TWO EXPRESSION MATRICES
# =============================================================================
obj  <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met  <- obj$met
expr_meth <- obj$exp                    # the matrix used for methylation work

get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^get_matrix(raw) - 0.001; lin[lin < 0] <- 0
expr_toil <- log2(lin + 1)              # the harmonised GRCh38 matrix

cnv_lin <- readRDS(file.path(CACHE, "tcga_cna_linear_gst.rds"))
cnv_dis <- tryCatch(readRDS(file.path(CACHE, "tcga_cna_gistic_gst.rds")),
                    error = function(e) NULL)
if (is.list(cnv_lin) && !is.matrix(cnv_lin))
  cnv_lin <- cnv_lin[[which(vapply(cnv_lin, is.matrix, logical(1)))[1]]]
if (!is.null(cnv_dis) && is.list(cnv_dis) && !is.matrix(cnv_dis))
  cnv_dis <- cnv_dis[[which(vapply(cnv_dis, is.matrix, logical(1)))[1]]]

st <- TCGAquery_subtype(tumor = "brca")

# Toil sample ids lack the trailing letter; harmonise to the 15-character form
norm_id <- function(v) substr(v, 1, 15)
colnames(expr_meth) <- norm_id(colnames(expr_meth))
colnames(met)       <- norm_id(colnames(met))
colnames(cnv_lin)   <- norm_id(colnames(cnv_lin))
if (!is.null(cnv_dis)) colnames(cnv_dis) <- norm_id(colnames(cnv_dis))
colnames(expr_toil) <- norm_id(colnames(expr_toil))

s <- Reduce(intersect, list(colnames(met), colnames(expr_meth),
                            colnames(expr_toil), colnames(cnv_lin)))
s <- s[substr(s, 13, 15) == "-01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(s, 1, 12), st$patient)]
ok <- grp %in% SUB
s <- s[ok]; grp <- factor(grp[ok], levels = SUB)

genes <- Reduce(intersect, list(rownames(met), rownames(expr_meth),
                                rownames(expr_toil), rownames(cnv_lin)))
cat("Tumours present in all four matrices with a PAM50 call:", length(s), "\n")
cat("Genes present in all four:", length(genes), "\n")
if (length(s) < 100) stop("Too few shared samples. Check the identifier formats.")


# =============================================================================
# 2. POINT 2: HOW FAR DO THE TWO EXPRESSION MATRICES AGREE
# =============================================================================
cat("\n", strrep("=", 76), "\nA. AGREEMENT BETWEEN THE TWO EXPRESSION MATRICES\n",
    strrep("=", 76), "\n", sep = "")

agree <- lapply(genes, function(g) {
  a <- as.numeric(expr_meth[g, s]); b <- as.numeric(expr_toil[g, s])
  okk <- is.finite(a) & is.finite(b)
  if (sum(okk) < 50) return(NULL)
  r <- suppressWarnings(cor(a[okk], b[okk], method = "spearman"))
  data.frame(gene = g, n = sum(okk), rho = round(r, 3),
             interpretable = !(g %in% EXCLUDE), stringsAsFactors = FALSE)
}) %>% bind_rows()
agree <- agree[order(agree$rho), ]
print(as.data.frame(agree), row.names = FALSE, digits = 3)

int <- agree[agree$interpretable, ]
cat(sprintf("\nAcross the %d genes used for inference: median rho %.3f, range %.3f to %.3f\n",
            nrow(int), median(int$rho), min(int$rho), max(int$rho)))
cat("A high correlation means the choice of matrix is largely immaterial for\n")
cat("that gene. A low one means the two quantifications disagree and any result\n")
cat("computed on one should not be read across to the other.\n")
write.csv(agree, file.path(OUT, "expression_matrix_agreement.csv"), row.names = FALSE)


# =============================================================================
# 3. THE PARTITION FUNCTION, REPORTING EVERY COMPONENT
# =============================================================================
# Marginal R2 is each predictor fitted alone. Shared variance is what the two
# predictors explain in common: marginal(m) + marginal(c) - full. Reporting it
# alongside the unique terms is what the reviewer asked for.
partition <- function(g, E, C, label) {
  # Index on the samples this matrix actually has. The discrete GISTIC matrix
  # covers far fewer tumours than the linear one, and indexing it with the full
  # sample vector fails rather than subsetting.
  if (!(g %in% rownames(C))) return(NULL)
  ss <- intersect(s, colnames(C))
  if (length(ss) < 100) return(NULL)
  e <- rank(as.numeric(E[g, ss])); m <- rank(as.numeric(met[g, ss]))
  cc <- as.numeric(C[g, ss])
  cc <- if (length(unique(cc[is.finite(cc)])) > 6) rank(cc) else cc
  df <- data.frame(e = e, m = m, c = cc)
  df <- df[complete.cases(df), ]
  if (nrow(df) < 40 || length(unique(df$c)) < 2) return(NULL)
  r2 <- function(f) summary(lm(as.formula(f), data = df))$r.squared
  full <- r2("e ~ m + c"); mm <- r2("e ~ m"); mc <- r2("e ~ c")
  data.frame(matrix = label, gene = g, n = nrow(df),
             R2_full = round(full, 3),
             marginal_meth = round(mm, 3), marginal_cn = round(mc, 3),
             unique_meth = round(full - mc, 3), unique_cn = round(full - mm, 3),
             shared = round(mm + mc - full, 3), stringsAsFactors = FALSE)
}


# =============================================================================
# 4. POINT 2: THE PARTITION ON EACH MATRIX
# =============================================================================
cat("\n", strrep("=", 76), "\nB. JOINT PARTITION ON EACH EXPRESSION MATRIX\n",
    strrep("=", 76), "\n", sep = "")

p_meth <- lapply(genes, function(g) partition(g, expr_meth, cnv_lin, "methylation matrix")) %>% bind_rows()
p_toil <- lapply(genes, function(g) partition(g, expr_toil, cnv_lin, "harmonised GRCh38")) %>% bind_rows()

cmp <- dplyr::left_join(
  p_meth %>% dplyr::select(gene, meth_matrix = unique_meth),
  p_toil %>% dplyr::select(gene, harmonised = unique_meth), by = "gene")
cmp$difference <- round(cmp$harmonised - cmp$meth_matrix, 3)
cmp$interpretable <- !(cmp$gene %in% EXCLUDE)
print(as.data.frame(cmp[order(-cmp$meth_matrix), ]), row.names = FALSE, digits = 3)

ci <- cmp[cmp$interpretable & is.finite(cmp$harmonised), ]
cat(sprintf("\nAcross interpretable genes, unique methylation R-squared differs between\n"))
cat(sprintf("matrices by a median of %.3f (range %.3f to %.3f).\n",
            median(ci$difference), min(ci$difference), max(ci$difference)))


# =============================================================================
# 5. POINT 5: THE FULL VARIANCE DECOMPOSITION, AND DISCRETE COPY NUMBER
# =============================================================================
cat("\n", strrep("=", 76), "\nC. FULL DECOMPOSITION, HARMONISED MATRIX, LINEAR COPY NUMBER\n",
    strrep("=", 76), "\n", sep = "")
show <- p_toil[!(p_toil$gene %in% EXCLUDE), ]
print(as.data.frame(show[order(-show$unique_meth),
      c("gene","n","R2_full","marginal_meth","marginal_cn","unique_meth","unique_cn","shared")]),
      row.names = FALSE, digits = 3)

if (!is.null(cnv_dis)) {
  cat("\n", strrep("=", 76), "\nD. SENSITIVITY: DISCRETE GISTIC CALLS IN PLACE OF LINEAR COPY NUMBER\n",
      strrep("=", 76), "\n", sep = "")
  sd_ <- intersect(s, colnames(cnv_dis))
  cat("Tumours with discrete calls:", length(sd_), "\n")
  cat("Distinct values in the discrete matrix:",
      length(unique(as.numeric(cnv_dis[, sd_]))), "\n\n")
  p_dis <- lapply(intersect(genes, rownames(cnv_dis)), function(g)
    partition(g, expr_toil, cnv_dis, "discrete GISTIC")) %>% bind_rows()
  # An empty result has no columns, so the join must be guarded rather than
  # attempted and allowed to fail.
  if (nrow(p_dis)) {
    cd <- dplyr::left_join(
      p_toil %>% dplyr::select(gene, linear_cn = unique_cn, linear_meth = unique_meth),
      p_dis  %>% dplyr::select(gene, discrete_cn = unique_cn, discrete_meth = unique_meth),
      by = "gene")
    cd <- cd[!(cd$gene %in% EXCLUDE), ]
    print(as.data.frame(cd[order(-cd$linear_meth), ]), row.names = FALSE, digits = 3)
  } else {
    cat("No gene could be modelled on the discrete calls. The discrete matrix\n")
    cat("covers too few tumours at too few distinct levels for a stable fit, so\n")
    cat("this sensitivity analysis cannot be performed with the data available\n")
    cat("and is reported as a limitation rather than omitted.\n")
    p_dis <- NULL
  }
} else {
  cat("\nNo discrete copy number object found; the GISTIC sensitivity analysis\n")
  cat("could not be run and this should be stated as a limitation.\n")
  p_dis <- NULL
}


# =============================================================================
# 6. VERDICT
# =============================================================================
cat("\n", strrep("=", 76), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 76), "\n", sep = "")

pass_marginal <- sum(show$marginal_meth > show$marginal_cn, na.rm = TRUE)
pass_unique   <- sum(show$unique_meth   > show$unique_cn,   na.rm = TRUE)
cat(sprintf("Of %d interpretable genes on the harmonised matrix:\n", nrow(show)))
cat(sprintf("  methylation exceeds copy number on marginal R-squared: %d\n", pass_marginal))
cat(sprintf("  methylation exceeds copy number on unique R-squared:   %d\n", pass_unique))
if (!is.null(p_dis)) {
  cdd <- dplyr::left_join(show %>% dplyr::select(gene, unique_meth),
                          p_dis %>% dplyr::select(gene, unique_cn), by = "gene")
  cat(sprintf("  methylation exceeds discrete copy number on unique R-squared: %d\n",
              sum(cdd$unique_meth > cdd$unique_cn, na.rm = TRUE)))
}
cat("\nIf methylation leads on all three, 'dominant' is defensible and should be\n")
cat("stated as surviving a fairer parameterisation. If it leads on unique but\n")
cat("not marginal R-squared, the resolution objection has force and the phrasing\n")
cat("should become 'explains substantially more variance in this model'.\n")

allout <- bind_rows(p_meth, p_toil, if (!is.null(p_dis)) p_dis)
write.csv(allout, file.path(OUT, "variance_partition_full_decomposition.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Shared tumours: %d", length(s)),
             sprintf("Expression matrix agreement, interpretable genes: median rho %.3f (%.3f to %.3f)",
                     median(int$rho), min(int$rho), max(int$rho)),
             sprintf("Methylation leads on marginal R2 for %d of %d genes", pass_marginal, nrow(show)),
             sprintf("Methylation leads on unique R2 for %d of %d genes", pass_unique, nrow(show)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_matrix_and_partition.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
