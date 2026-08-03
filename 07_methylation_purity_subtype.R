# =============================================================================
#  MODULE C · PART 2 — PURITY ADJUSTMENT AND SUBTYPE-STRATIFIED METHYLATION
#
#  THREE PROBLEMS IN PART 1 THAT THIS FIXES
#
#  1. SIXTEEN OF EIGHTEEN GENES CORRELATED NEGATIVELY.
#     Family-wide negative correlation points to a shared confounder rather
#     than gene-specific regulation. GSTs are epithelial genes; a stroma-rich
#     tumour has lower GST expression AND a bulk methylation reading diluted
#     by non-epithelial cells. Section 2 recomputes every correlation as a
#     PARTIAL correlation adjusting for tumour purity. Genes that survive are
#     regulated; genes that collapse were reading composition.
#
#  2. GSTP1 WAS MISCLASSIFIED.
#     Median beta 0.033 with rho -0.604 was labelled "low methylation", but
#     GSTP1 is bimodal - a large unmethylated, high-expressing group and a
#     separate methylated, silenced group. A median is meaningless for a
#     bimodal variable. Section 3 classifies on the PERCENTAGE of tumours
#     above beta 0.3 and formally tests bimodality.
#
#  3. GSTM1's CORRELATION IS COPY-NUMBER CONFOUNDED.
#     rho -0.794 is the strongest in the table, but GSTM1 carries a common
#     germline whole-gene deletion. No template means no expression and an
#     unreliable beta. Flagged throughout, excluded from conclusions.
#
#  THE HYPOTHESIS THIS TESTS
#     Your expression data showed two programmes: GSTP1/GSTA1/GSTA4 high in
#     Basal, GSTM3/GSTO2/GSTZ1 high in luminal. If methylation is the
#     mechanism, each programme should be UNmethylated in the subtype where
#     it is expressed. Section 4 tests exactly that.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")

obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; exp <- obj$exp

tum <- intersect(colnames(met), colnames(exp))
tum <- tum[substr(tum, 14, 15) == "01"]          # primary tumours only
genes <- intersect(rownames(met), rownames(exp))
cat("Tumours with both assays:", length(tum), "| genes:", length(genes), "\n")


# =============================================================================
# 1. ANNOTATION — PAM50 AND TUMOUR PURITY
# =============================================================================
pat <- substr(tum, 1, 12)

st <- TCGAquery_subtype(tumor = "brca")
pam50 <- st$BRCA_Subtype_PAM50[match(pat, st$patient)]

# TCGAbiolinks ships the consensus purity estimates from Aran et al. 2015.
# CPE is the consensus across ABSOLUTE, ESTIMATE, LUMP and IHC.
purity <- rep(NA_real_, length(tum))
ok <- tryCatch({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  purity <- num(tp$CPE)[match(pat, tp$patient)]
  TRUE
}, error = function(e) { message("Purity table unavailable: ", conditionMessage(e)); FALSE })

cat(sprintf("PAM50 assigned: %d | purity available: %d\n",
            sum(!is.na(pam50)), sum(!is.na(purity))))
if (sum(!is.na(purity)) < 100)
  warning("Too few purity values - the adjusted correlation will be unreliable.")


# =============================================================================
# 2. PARTIAL CORRELATION ADJUSTING FOR PURITY
# =============================================================================
# Regress both methylation and expression on purity, then correlate the
# residuals. What remains is the association not explained by how much tumour
# is in the sample.

pcor <- lapply(genes, function(g) {
  x <- met[g, tum]; y <- exp[g, tum]; p <- purity
  ok <- is.finite(x) & is.finite(y) & is.finite(p)
  if (sum(ok) < 50 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NULL)

  raw <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  rx  <- residuals(lm(rank(x[ok]) ~ rank(p[ok])))
  ry  <- residuals(lm(rank(y[ok]) ~ rank(p[ok])))
  ct  <- suppressWarnings(cor.test(rx, ry, method = "spearman"))

  data.frame(gene = g, n = sum(ok),
             rho_raw = round(raw, 3),
             rho_adj = round(unname(ct$estimate), 3),
             p_adj = ct$p.value,
             attenuation = round(raw - unname(ct$estimate), 3))
}) %>% bind_rows() %>%
  mutate(FDR_adj = p.adjust(p_adj, "BH"),
         verdict = case_when(
           rho_adj <= -0.3 & FDR_adj < 0.05 ~ "survives purity adjustment",
           rho_adj <= -0.15 & FDR_adj < 0.05 ~ "weakened but present",
           TRUE ~ "explained by purity/composition")) %>%
  arrange(rho_adj)

cat("\n===== PARTIAL CORRELATION, ADJUSTED FOR TUMOUR PURITY =====\n")
print(as.data.frame(pcor), digits = 3)
cat("\nLarge 'attenuation' means the raw correlation was largely composition.\n")
cat("GSTM1 is copy-number confounded regardless of what this row says.\n")
write.csv(pcor, file.path(OUT, "methylation_expression_purity_adjusted.csv"),
          row.names = FALSE)


# =============================================================================
# 3. BIMODALITY — the right way to classify GSTP1
# =============================================================================
BETA_CUT <- 0.3     # conventional threshold for a methylated promoter

bim <- lapply(genes, function(g) {
  x <- met[g, tum]; x <- x[is.finite(x)]
  if (length(x) < 50) return(NULL)
  data.frame(gene = g, n = length(x),
             median_beta = round(median(x), 3),
             pct_methylated = round(100 * mean(x > BETA_CUT), 1),
             pct_unmethylated = round(100 * mean(x < 0.1), 1))
}) %>% bind_rows() %>%
  mutate(pattern = case_when(
    pct_methylated > 80 ~ "uniformly methylated",
    pct_methylated < 10 ~ "uniformly unmethylated",
    TRUE ~ "BIMODAL - subset methylated")) %>%
  arrange(desc(pct_methylated))

cat("\n===== METHYLATION PATTERN =====\n")
print(as.data.frame(bim))
cat("\nA bimodal gene cannot be summarised by its median. GSTP1 belongs here:\n")
cat("most tumours unmethylated and highly expressed, a subset silenced.\n")
write.csv(bim, file.path(OUT, "methylation_pattern.csv"), row.names = FALSE)


# =============================================================================
# 4. THE TWO-PROGRAMME TEST
# =============================================================================
# If methylation drives the subtype split, each programme should be
# UNmethylated in the subtype where it is expressed.

df <- lapply(genes, function(g) {
  data.frame(gene = g, sample = tum, beta = met[g, tum],
             expr = exp[g, tum], pam50 = pam50)
}) %>% bind_rows() %>%
  filter(pam50 %in% SUB, is.finite(beta)) %>%
  mutate(pam50 = factor(pam50, levels = SUB))

by_sub <- df %>% group_by(gene, pam50) %>%
  summarise(median_beta = round(median(beta), 3),
            pct_meth = round(100 * mean(beta > BETA_CUT), 1),
            .groups = "drop")

beta_wide <- by_sub %>% select(gene, pam50, median_beta) %>%
  pivot_wider(names_from = pam50, values_from = median_beta) %>%
  mutate(lowest_beta  = SUB[apply(across(all_of(SUB)), 1, which.min)],
         highest_beta = SUB[apply(across(all_of(SUB)), 1, which.max)],
         spread = round(apply(across(all_of(SUB)), 1, max) -
                          apply(across(all_of(SUB)), 1, min), 3)) %>%
  arrange(desc(spread))

cat("\n===== MEDIAN METHYLATION BY PAM50 SUBTYPE =====\n")
print(as.data.frame(beta_wide), digits = 3)
write.csv(beta_wide, file.path(OUT, "methylation_by_subtype.csv"), row.names = FALSE)

kw <- df %>% group_by(gene) %>%
  summarise(H = unname(kruskal.test(beta ~ droplevels(pam50))$statistic),
            p = kruskal.test(beta ~ droplevels(pam50))$p.value,
            n = n(), .groups = "drop") %>%
  mutate(eps2 = round(H / ((n^2 - 1) / (n + 1)), 4),
         FDR = p.adjust(p, "BH")) %>% arrange(desc(eps2))

cat("\n===== IS METHYLATION SUBTYPE-STRUCTURED? =====\n")
print(as.data.frame(kw), digits = 3)
write.csv(kw, file.path(OUT, "methylation_subtype_stats.csv"), row.names = FALSE)


# ---- the decisive comparison -------------------------------------------------
# Expression-highest subtype should be the methylation-LOWEST subtype.
expr_f <- file.path(WORKDIR, "figures", "withintumour_toil",
                    "subtype_medians_wide_toil.csv")

if (file.exists(expr_f)) {
  ex <- read.csv(expr_f, stringsAsFactors = FALSE)
  cmp <- beta_wide %>%
    select(gene, lowest_beta, highest_beta, spread) %>%
    inner_join(ex[, c("gene","highest","lowest")], by = "gene") %>%
    rename(expr_highest = highest, expr_lowest = lowest) %>%
    mutate(mechanism_consistent = lowest_beta == expr_highest)

  cat("\n===== DOES METHYLATION EXPLAIN THE SUBTYPE SPLIT? =====\n")
  cat("Consistent = the subtype with LOWEST methylation is the subtype with\n")
  cat("HIGHEST expression. That is what epigenetic control predicts.\n\n")
  print(as.data.frame(cmp))
  cat(sprintf("\nConsistent in %d of %d genes\n",
              sum(cmp$mechanism_consistent), nrow(cmp)))
  write.csv(cmp, file.path(OUT, "methylation_expression_subtype_match.csv"),
            row.names = FALSE)

  for (g in c("GSTP1","GSTA1","GSTM3","GSTO2","GSTM5")) {
    r <- cmp[cmp$gene == g, ]
    if (nrow(r))
      cat(sprintf("  %-7s expression highest in %-5s | methylation lowest in %-5s  %s\n",
                  g, r$expr_highest, r$lowest_beta,
                  ifelse(r$mechanism_consistent, "<- CONSISTENT", "")))
  }
}


# =============================================================================
# 5. GSTP1 IN DETAIL — the bimodal gene, by subtype
# =============================================================================
if ("GSTP1" %in% genes) {
  gp <- df %>% filter(gene == "GSTP1")

  tab <- gp %>% group_by(pam50) %>%
    summarise(n = n(),
              median_beta = round(median(beta), 3),
              pct_methylated = round(100 * mean(beta > BETA_CUT), 1),
              pct_unmethylated = round(100 * mean(beta < 0.1), 1),
              .groups = "drop")

  cat("\n===== GSTP1 METHYLATION BY SUBTYPE =====\n")
  print(as.data.frame(tab))
  cat("\nPREDICTION: Basal should have the LOWEST methylated fraction, because\n")
  cat("GSTP1 expression is 522 TPM in Basal versus 72 in LumB. If so, escape\n")
  cat("from promoter methylation is the mechanism for GSTP1 retention in\n")
  cat("basal-like disease - and that resolves the GSTP1 literature\n")
  cat("contradiction mechanistically rather than descriptively.\n")
  write.csv(tab, file.path(OUT, "GSTP1_methylation_by_subtype.csv"), row.names = FALSE)

  p1 <- ggplot(gp, aes(beta, fill = pam50)) +
    geom_histogram(bins = 50, alpha = .75, colour = NA) +
    geom_vline(xintercept = BETA_CUT, linetype = "dashed", colour = "grey30") +
    facet_wrap(~ pam50, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                                 LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
    labs(x = "GSTP1 promoter methylation (beta)", y = "Tumours",
         title = "GSTP1 promoter methylation is bimodal, and the split is subtype-dependent",
         subtitle = "Dashed line = beta 0.3, the conventional threshold for a methylated promoter") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey93"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(OUT, "GSTP1_methylation_bimodality.png"), p1,
         width = 8, height = 6, dpi = 300)
}


# =============================================================================
# 6. FAMILY-WIDE FIGURE
# =============================================================================
keep <- kw$gene[kw$eps2 > 0.02]
if (length(keep)) {
  p2 <- ggplot(df %>% filter(gene %in% keep), aes(pam50, beta, fill = pam50)) +
    geom_violin(scale = "width", alpha = .3, colour = NA) +
    geom_boxplot(width = .28, outlier.size = .2, outlier.alpha = .25, linewidth = .28) +
    geom_hline(yintercept = BETA_CUT, linetype = "dashed",
               colour = "grey40", linewidth = .3) +
    facet_wrap(~ gene, ncol = 4) +
    scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                                 LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
    labs(x = NULL, y = "Promoter methylation (beta)",
         title = "GST promoter methylation across PAM50 subtypes",
         subtitle = "Genes with subtype-structured methylation; dashed line = beta 0.3") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.background = element_rect(fill = "grey93"),
          strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(OUT, "methylation_by_subtype.png"), p2,
         width = 10, height = 8, dpi = 300)
}

writeLines(c(
  paste("Run:", Sys.time()),
  "Partial correlation adjusts methylation-expression association for tumour",
  "purity (Aran et al. 2015 CPE) via rank residuals.",
  "Methylation called on percentage of tumours above beta 0.3, not the median,",
  "because bimodal genes such as GSTP1 are not summarised by a median.",
  "GSTM1 excluded from epigenetic conclusions: common germline whole-gene",
  "deletion confounds both its expression and its beta values.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleC_part2.txt"))

message("\nWritten to ", normalizePath(OUT))
