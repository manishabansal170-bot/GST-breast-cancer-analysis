# =============================================================================
#  42_methylation_on_harmonised_matrix.R
#
#  THE OBJECTION
#    The manuscript uses two expression matrices. The expression architecture,
#    Table 2 and Figures 1 to 3 use the uniformly requantified GRCh38 Toil
#    recompute. The methylation correlations, the joint variance partition in
#    Table 4 and Figure 4 use the TCGA expression matrix distributed alongside
#    the methylation data, which is differently annotated.
#
#    The manuscript is transparent about this, and a cross-matrix sensitivity
#    analysis showed the classification does not change and unique methylation
#    variance moves by a median of 0.002. But a reader who meets GSTP1 at
#    522.2 TPM in one table, 519.1 in another and an unstated value behind a
#    figure will lose the thread, and the reviewer is right that carrying one
#    matrix is worth the small cost.
#
#  WHAT THIS SCRIPT DOES
#    Recomputes, on the harmonised GRCh38 matrix:
#      1. The methylation-expression correlations, raw and purity-adjusted
#      2. The joint variance partition with p values and driver classification
#      3. The subtype-stratified partition
#    and writes them as the primary output, with the original-matrix values
#    retained alongside so the supplementary sensitivity comparison can be
#    built from one file.
#
#  A CONSEQUENCE THAT MUST BE REPORTED
#    GSTT1 is absent from the GRCh38 primary assembly, so it does not exist in
#    the harmonised matrix. The joint model therefore covers 17 genes rather
#    than 18. This is not a loss: the manuscript already states that GSTT1 is
#    not quantifiable on this assembly and draws no conclusion from it. Moving
#    to one matrix simply makes the gene list consistent with that position
#    instead of carrying a gene the paper elsewhere says cannot be measured.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleBG_alterations_protein")
MOUT    <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB     <- c("Basal","Her2","LumA","LumB")
EXCLUDE <- c("GSTM1","GSTT2","GSTA2","GSTT1")


# =============================================================================
# 1. INPUTS, ON ONE ANNOTATION
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met
expr_old <- obj$exp                       # retained only for the comparison

gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^gm(raw) - 0.001; lin[lin < 0] <- 0
expr <- log2(lin + 1)                     # the harmonised GRCh38 matrix

cnv <- readRDS(file.path(CACHE, "tcga_cna_linear_gst.rds"))
if (is.list(cnv) && !is.matrix(cnv)) cnv <- cnv[[which(vapply(cnv, is.matrix, logical(1)))[1]]]

# Harmonise sample identifiers to the 15-character form. Written out one at a
# time rather than in a loop, because assign() inside a for loop is easy to get
# subtly wrong and a silent failure here produces a subscript error much later.
id <- function(v) substr(v, 1, 15)
colnames(met)      <- id(colnames(met))
colnames(expr_old) <- id(colnames(expr_old))
colnames(expr)     <- id(colnames(expr))
colnames(cnv)      <- id(colnames(cnv))

st <- TCGAquery_subtype(tumor = "brca")
# expr_old is included in the intersection even though it is only used for the
# sensitivity column: one tumour carries methylation, harmonised expression and
# copy number but no value in the original matrix, and omitting it here produces
# a subscript error far downstream rather than at the point of the mismatch.
s <- Reduce(intersect, list(colnames(met), colnames(expr),
                            colnames(expr_old), colnames(cnv)))
s <- s[substr(s, 13, 15) == "-01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(s, 1, 12), st$patient)]
ok <- grp %in% SUB; s <- s[ok]; grp <- factor(grp[ok], levels = SUB)

genes <- Reduce(intersect, list(rownames(met), rownames(expr), rownames(cnv)))
cat("Tumours with methylation, harmonised expression and copy number:", length(s), "\n")
print(table(grp))
cat("Genes present in all three on this annotation:", length(genes), "\n")
absent <- setdiff(rownames(met), genes)
if (length(absent))
  cat("Genes in the methylation set but absent from the harmonised matrix:",
      paste(absent, collapse = ", "), "\n")
cat("\n")

# purity, for the adjusted correlations
pur <- rep(NA_real_, length(s))
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  pur <- num(tp$CPE)[match(substr(s, 1, 12), tp$patient)]
}, silent = TRUE)
cat("Tumours with a purity estimate:", sum(is.finite(pur)), "\n\n")


# =============================================================================
# 2. METHYLATION-EXPRESSION CORRELATIONS
# =============================================================================
cat(strrep("=", 78), "\nA. METHYLATION AGAINST EXPRESSION, HARMONISED MATRIX\n",
    strrep("=", 78), "\n", sep = "")

corr <- lapply(genes, function(g) {
  m <- as.numeric(met[g, s]); e <- as.numeric(expr[g, s])
  eo <- if (g %in% rownames(expr_old)) as.numeric(expr_old[g, s]) else rep(NA_real_, length(s))
  okk <- is.finite(m) & is.finite(e)
  if (sum(okk) < 50) return(NULL)
  raw_r <- suppressWarnings(cor(m[okk], e[okk], method = "spearman"))
  old_r <- if (all(is.na(eo))) NA_real_ else suppressWarnings(cor(m[okk], eo[okk], method = "spearman", use = "pairwise.complete.obs"))
  # purity-adjusted: rank residuals after regressing on purity
  adj_r <- NA_real_
  p2 <- is.finite(pur) & okk
  if (sum(p2) >= 50) {
    pr <- rank(pur[p2])
    mr <- resid(lm(rank(m[p2]) ~ pr)); er <- resid(lm(rank(e[p2]) ~ pr))
    adj_r <- suppressWarnings(cor(mr, er))
  }
  data.frame(gene = g, n = sum(okk),
             rho_harmonised = round(raw_r, 3),
             rho_purity_adjusted = round(adj_r, 3),
             rho_original_matrix = round(old_r, 3),
             interpretable = !(g %in% EXCLUDE), stringsAsFactors = FALSE)
}) %>% bind_rows()
corr <- corr[order(corr$rho_harmonised), ]
print(as.data.frame(corr), row.names = FALSE, digits = 3)

ci <- corr[corr$interpretable, ]
cat(sprintf("\nAcross interpretable genes, the harmonised and original coefficients\n"))
cat(sprintf("differ by a median of %.3f (range %.3f to %.3f).\n",
            median(ci$rho_harmonised - ci$rho_original_matrix),
            min(ci$rho_harmonised - ci$rho_original_matrix),
            max(ci$rho_harmonised - ci$rho_original_matrix)))
write.csv(corr, file.path(MOUT, "methylation_expression_harmonised.csv"), row.names = FALSE)


# =============================================================================
# 3. JOINT VARIANCE PARTITION, WITH P VALUES AND DRIVER CALL
# =============================================================================
partition <- function(g, idx, label) {
  e <- rank(as.numeric(expr[g, s[idx]]))
  m <- rank(as.numeric(met [g, s[idx]]))
  c_ <- rank(as.numeric(cnv [g, s[idx]]))
  df <- data.frame(e = e, m = m, c = c_); df <- df[complete.cases(df), ]
  if (nrow(df) < 40) return(NULL)
  r2 <- function(f) summary(lm(as.formula(f), data = df))$r.squared
  full <- r2("e ~ m + c"); mm <- r2("e ~ m"); mc <- r2("e ~ c")
  fit <- lm(e ~ m + c, data = df); co <- summary(fit)$coefficients
  data.frame(stratum = label, gene = g, n = nrow(df),
             R2_total = round(full, 3),
             unique_cn = round(full - mm, 3),
             unique_meth = round(full - mc, 3),
             shared = round(mm + mc - full, 3),
             p_cn   = tryCatch(co["c", 4], error = function(e) NA_real_),
             p_meth = tryCatch(co["m", 4], error = function(e) NA_real_),
             stringsAsFactors = FALSE)
}

cat("\n", strrep("=", 78), "\nB. JOINT PARTITION, HARMONISED MATRIX\n",
    strrep("=", 78), "\n", sep = "")
pool <- lapply(genes, function(g) partition(g, seq_along(s), "pooled")) %>% bind_rows()

# the published classification rule, restated here so it is visible
pool$driver <- with(pool, ifelse(
  unique_meth > 2 * unique_cn & p_meth < 0.05, "Methylation",
  ifelse(unique_cn > 2 * unique_meth & p_cn < 0.05, "Copy number", "Both")))
pool$interpretable <- !(pool$gene %in% EXCLUDE)

print(as.data.frame(pool[order(-pool$unique_meth),
      c("gene","n","R2_total","unique_cn","unique_meth","shared","p_cn","p_meth","driver")]),
      row.names = FALSE, digits = 3)

int <- pool[pool$interpretable, ]
cat(sprintf("\nInterpretable genes: %d, of which methylation-driven: %d\n",
            nrow(int), sum(int$driver == "Methylation")))
write.csv(pool, file.path(OUT, "variance_partition_harmonised.csv"), row.names = FALSE)


# =============================================================================
# 4. WITHIN SUBTYPE, ON THE SAME MATRIX
# =============================================================================
cat("\n", strrep("=", 78), "\nC. WITHIN SUBTYPE, HARMONISED MATRIX\n",
    strrep("=", 78), "\n", sep = "")
wi <- lapply(SUB, function(k) {
  idx <- which(grp == k)
  lapply(genes, function(g) partition(g, idx, k)) %>% bind_rows()
}) %>% bind_rows()

focus <- c("GSTP1","GSTM3","GSTM5","GSTM2","GSTM4","GSTO2","GSTA4","GSTA1")
w <- wi %>% dplyr::filter(gene %in% focus) %>%
  dplyr::select(gene, stratum, unique_meth) %>%
  tidyr::pivot_wider(names_from = stratum, values_from = unique_meth)
p_ <- pool %>% dplyr::filter(gene %in% focus) %>% dplyr::select(gene, pooled = unique_meth)
print(as.data.frame(dplyr::left_join(p_, w, by = "gene")), row.names = FALSE, digits = 3)

write.csv(bind_rows(pool, wi),
          file.path(OUT, "variance_partition_harmonised_all.csv"), row.names = FALSE)


# =============================================================================
# 5. WHAT TO CHANGE IN THE MANUSCRIPT
# =============================================================================
cat("\n", strrep("=", 78), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("Table 4 should be rebuilt from section B: %d genes, n = %d.\n",
            nrow(pool), pool$n[1]))
if (length(absent))
  cat(sprintf("It will lose %s, which is absent from the GRCh38 assembly and which\n",
              paste(absent, collapse = ", ")),
      "the manuscript already states cannot be quantified there.\n", sep = "")
cat("The methylation-expression correlations in section A replace the current\n")
cat("supplementary table, with the original-matrix column retained as the\n")
cat("sensitivity comparison the reviewer asked to be relegated.\n")
cat("Figure 4's expression row must also be regenerated from this matrix; the\n")
cat("figure script reads obj$exp and should be pointed at the Toil object.\n")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d, genes: %d", length(s), length(genes)),
             sprintf("Absent from harmonised matrix: %s",
                     if (length(absent)) paste(absent, collapse = ", ") else "none"),
             sprintf("Methylation-driven among interpretable: %d of %d",
                     sum(int$driver == "Methylation"), nrow(int)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_harmonised_methylation.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
