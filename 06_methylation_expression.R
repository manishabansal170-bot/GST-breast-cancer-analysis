# =============================================================================
#  OBJECTIVE 1 · MODULE C — METHYLATION AND EXPRESSION
#
#  THE QUESTION THIS ANSWERS
#    You have shown the Mu class and GSTA1 are lost in tumour. That is a
#    description. This asks WHY: is the loss epigenetic?
#
#    A gene silenced by promoter hypermethylation shows (a) higher methylation
#    in tumour than normal, and (b) a negative correlation between methylation
#    and expression across samples. Both together is a mechanism. Either alone
#    is suggestive.
#
#  WHY THIS MATTERS FOR YOUR THESIS SPECIFICALLY
#    It converts Objective 1 from descriptive to mechanistic, and it is the
#    natural bridge to Objective 2 - you already have a North Indian methylation
#    cohort with GSTP1, p16 and BRCA1. If TCGA shows GSTP1 hypermethylation
#    and your own cohort shows the same, that is public data and wet lab
#    pointing the same way in two populations.
#
#  DATA: cBioPortal API. Small - only the GST genes. No large download.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("cBioPortalData", ask = FALSE, update = FALSE) }
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(cBioPortalData) })

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")


# =============================================================================
# 1. FIND A TCGA-BRCA STUDY WITH BOTH METHYLATION AND EXPRESSION
# =============================================================================
cbio <- cBioPortal()

studies <- getStudies(cbio)
cand <- studies$studyId[grepl("^brca_tcga", studies$studyId)]
cat("\n===== CANDIDATE TCGA-BRCA STUDIES =====\n"); print(cand)

STUDY <- NA; MET <- NA; EXPR <- NA

for (s in cand) {
  mp <- tryCatch(molecularProfiles(cbio, s), error = function(e) NULL)
  if (is.null(mp)) next
  ids <- mp$molecularProfileId
  m <- ids[grepl("methylation", ids, ignore.case = TRUE)]
  # Prefer a raw expression profile over a z-scored one.
  e <- ids[grepl("rna_seq.*mrna|mrna.*rna_seq", ids, ignore.case = TRUE) &
             !grepl("Zscore|zscore", ids)]
  if (!length(e)) e <- ids[grepl("mrna", ids) & !grepl("Zscore|zscore", ids)]
  if (length(m) && length(e)) {
    STUDY <- s
    # hm450 covers more CpGs than hm27; prefer it.
    MET  <- if (any(grepl("hm450", m))) m[grepl("hm450", m)][1] else m[1]
    EXPR <- e[1]
    break
  }
}

if (is.na(STUDY)) stop("No TCGA-BRCA study with both methylation and expression.")
cat(sprintf("\nStudy: %s\nMethylation profile: %s\nExpression profile: %s\n",
            STUDY, MET, EXPR))


# =============================================================================
# 2. FETCH  (cached)
# =============================================================================
cf <- file.path(CACHE, "tcga_methylation_gst.rds")

if (file.exists(cf)) {
  obj <- readRDS(cf); met <- obj$met; exp <- obj$exp
  message("Loaded from cache.")
} else {
  message("Fetching methylation...")
  md <- getDataByGenes(cbio, studyId = STUDY, genes = GST,
                       by = "hugoGeneSymbol", molecularProfileIds = MET)[[1]]
  message("Fetching expression...")
  ed <- getDataByGenes(cbio, studyId = STUDY, genes = GST,
                       by = "hugoGeneSymbol", molecularProfileIds = EXPR)[[1]]

  to_mat <- function(d) {
    w <- d %>% select(hugoGeneSymbol, sampleId, value) %>%
      filter(!is.na(value)) %>%
      group_by(hugoGeneSymbol, sampleId) %>%
      summarise(value = mean(value), .groups = "drop") %>%
      pivot_wider(names_from = sampleId, values_from = value)
    m <- as.matrix(w[, -1]); rownames(m) <- w$hugoGeneSymbol; m
  }
  met <- to_mat(md); exp <- to_mat(ed)
  saveRDS(list(met = met, exp = exp), cf)
  message("Cached.")
}

cat(sprintf("\nMethylation: %d genes x %d samples\n", nrow(met), ncol(met)))
cat(sprintf("Expression : %d genes x %d samples\n", nrow(exp), ncol(exp)))

# A methylation beta value runs 0-1. If the range says otherwise, these are
# M-values and the direction of interpretation still holds but the scale does not.
cat(sprintf("Methylation value range: %.3f to %.3f %s\n",
            min(met, na.rm = TRUE), max(met, na.rm = TRUE),
            if (max(met, na.rm = TRUE) <= 1.05) "(beta values)" else "(NOT beta - check scale)"))


# =============================================================================
# 3. TUMOUR VS NORMAL METHYLATION
# =============================================================================
# TCGA sample barcodes end in a two-digit sample-type code:
#   01 = primary tumour, 11 = solid tissue normal.
stype <- substr(colnames(met), 14, 15)
grp <- ifelse(stype == "11", "Normal", ifelse(stype == "01", "Tumour", NA))

cat("\n===== METHYLATION SAMPLES =====\n"); print(table(grp, useNA = "ifany"))

if (sum(grp == "Normal", na.rm = TRUE) >= 10) {

  tvn <- lapply(rownames(met), function(g) {
    t_ <- met[g, grp == "Tumour" & !is.na(grp)]
    n_ <- met[g, grp == "Normal" & !is.na(grp)]
    t_ <- t_[!is.na(t_)]; n_ <- n_[!is.na(n_)]
    if (length(t_) < 10 || length(n_) < 5) return(NULL)
    data.frame(gene = g,
               beta_normal = round(median(n_), 3),
               beta_tumour = round(median(t_), 3),
               delta_beta  = round(median(t_) - median(n_), 3),
               p = suppressWarnings(wilcox.test(t_, n_, exact = FALSE)$p.value))
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, "BH"),
           # delta-beta above 0.2 is the conventional threshold for a
           # biologically meaningful methylation difference
           call = case_when(delta_beta >=  0.2 & FDR < 0.05 ~ "HYPERmethylated in tumour",
                            delta_beta <= -0.2 & FDR < 0.05 ~ "HYPOmethylated in tumour",
                            TRUE ~ "no meaningful change")) %>%
    arrange(desc(delta_beta))

  cat("\n===== TUMOUR vs NORMAL METHYLATION =====\n")
  print(as.data.frame(tvn), digits = 3)
  cat("\ndelta_beta > 0.2 with FDR < 0.05 is the conventional bar for a real\n")
  cat("methylation difference. Smaller shifts are common and rarely silencing.\n")
  write.csv(tvn, file.path(OUT, "methylation_tumour_vs_normal.csv"), row.names = FALSE)

} else {
  message("Too few normal methylation samples for a tumour-vs-normal test.")
  tvn <- NULL
}


# =============================================================================
# 4. METHYLATION vs EXPRESSION — the correlation that shows silencing
# =============================================================================
shared_s <- intersect(colnames(met), colnames(exp))
shared_g <- intersect(rownames(met), rownames(exp))
cat(sprintf("\nMatched: %d genes across %d samples\n",
            length(shared_g), length(shared_s)))

# Restrict to tumours: including normals can manufacture a correlation purely
# from the tumour/normal difference in both variables.
tum_s <- shared_s[substr(shared_s, 14, 15) == "01"]
cat("Tumour samples used for correlation:", length(tum_s), "\n")

corr <- lapply(shared_g, function(g) {
  x <- met[g, tum_s]; y <- exp[g, tum_s]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 30 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(gene = g, n = sum(ok),
             rho = round(unname(ct$estimate), 3),
             p = ct$p.value,
             median_beta = round(median(x[ok]), 3))
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "BH"),
         # A silenced gene needs BOTH a negative correlation AND appreciable
         # methylation. Strong negative rho at beta 0.02 is not silencing.
         interpretation = case_when(
           rho <= -0.3 & FDR < 0.05 & median_beta > 0.2 ~
             "methylation-associated silencing",
           rho <= -0.3 & FDR < 0.05 ~
             "negative correlation but low methylation",
           rho >= 0.3 & FDR < 0.05 ~ "positive correlation (unusual)",
           TRUE ~ "no clear relationship")) %>%
  arrange(rho)

cat("\n===== METHYLATION-EXPRESSION CORRELATION (tumours only) =====\n")
print(as.data.frame(corr), digits = 3)
write.csv(corr, file.path(OUT, "methylation_expression_correlation.csv"),
          row.names = FALSE)


# =============================================================================
# 5. VERDICT PER GENE — is the loss you found epigenetic?
# =============================================================================
if (!is.null(tvn)) {
  verdict <- full_join(tvn[, c("gene","beta_normal","beta_tumour","delta_beta","call")],
                       corr[, c("gene","rho","median_beta","interpretation")],
                       by = "gene") %>%
    mutate(epigenetic = case_when(
      grepl("HYPER", call) & rho <= -0.3 ~ "YES - hypermethylated AND anticorrelated",
      grepl("HYPER", call)               ~ "partial - hypermethylated, weak correlation",
      rho <= -0.3 & median_beta > 0.2    ~ "partial - anticorrelated, no tumour shift",
      TRUE                               ~ "no evidence")) %>%
    arrange(factor(epigenetic,
                   levels = c("YES - hypermethylated AND anticorrelated",
                              "partial - hypermethylated, weak correlation",
                              "partial - anticorrelated, no tumour shift",
                              "no evidence")))

  cat("\n===== IS THE LOSS EPIGENETIC? =====\n")
  print(as.data.frame(verdict), digits = 3)
  write.csv(verdict, file.path(OUT, "epigenetic_verdict.csv"), row.names = FALSE)

  cat("\n--- The genes your expression work flagged as lost ---\n")
  for (g in c("GSTM5","GSTM2","GSTA1","GSTP1","GSTM3")) {
    r <- verdict[verdict$gene == g, ]
    if (nrow(r)) cat(sprintf("  %-7s delta-beta %+.3f | rho %+.3f  ->  %s\n",
                             g, r$delta_beta, r$rho, r$epigenetic))
  }
  cat("\nGSTP1 is the one to watch. Promoter hypermethylation of GSTP1 in\n")
  cat("breast cancer is well established, so it doubles as a positive control:\n")
  cat("if the method cannot recover GSTP1, treat the other genes with caution.\n")
}


# =============================================================================
# 6. FIGURES
# =============================================================================
if (length(tum_s) > 30) {
  pl <- lapply(intersect(c("GSTP1","GSTM5","GSTM2","GSTA1","GSTM3","GSTM4"), shared_g),
               function(g) data.frame(gene = g, beta = met[g, tum_s],
                                      expr = exp[g, tum_s])) %>% bind_rows() %>%
    filter(is.finite(beta), is.finite(expr))

  if (nrow(pl)) {
    p <- ggplot(pl, aes(beta, expr)) +
      geom_point(alpha = .25, size = .7, colour = "#2E86AB") +
      geom_smooth(method = "lm", se = TRUE, colour = "#C0392B", linewidth = .6) +
      facet_wrap(~ gene, scales = "free_y", ncol = 3) +
      labs(x = "Promoter methylation (beta)", y = "Expression",
           title = "Methylation versus expression in TCGA-BRCA tumours",
           subtitle = "A negative slope with appreciable beta indicates epigenetic silencing") +
      theme_bw(base_size = 10) +
      theme(strip.background = element_rect(fill = "grey93"),
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", size = 12))
    ggsave(file.path(OUT, "methylation_vs_expression.png"), p,
           width = 9, height = 6, dpi = 300)
  }
}

if (!is.null(tvn)) {
  bx <- lapply(rownames(met), function(g) {
    data.frame(gene = g, beta = met[g, ], grp = grp)
  }) %>% bind_rows() %>% filter(!is.na(grp), is.finite(beta))

  q <- ggplot(bx, aes(grp, beta, fill = grp)) +
    geom_violin(scale = "width", alpha = .3, colour = NA) +
    geom_boxplot(width = .3, outlier.size = .3, outlier.alpha = .3, linewidth = .3) +
    facet_wrap(~ gene, ncol = 5) +
    scale_fill_manual(values = c(Normal = "#9AA5AC", Tumour = "#C0392B"), guide = "none") +
    labs(x = NULL, y = "Promoter methylation (beta)",
         title = "GST family promoter methylation, tumour versus adjacent normal",
         subtitle = "TCGA-BRCA HM450") +
    theme_bw(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey93"),
          strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold", size = 11))
  ggsave(file.path(OUT, "methylation_tumour_vs_normal.png"), q,
         width = 10, height = 8, dpi = 300)
}

writeLines(c(
  paste("Run:", Sys.time()),
  paste("Study:", STUDY),
  paste("Methylation profile:", MET),
  paste("Expression profile:", EXPR),
  "Correlation computed on primary tumours only, to avoid a spurious",
  "association driven by the tumour/normal difference in both variables.",
  "Silencing called on delta-beta > 0.2 with FDR < 0.05 AND rho <= -0.3.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleC.txt"))

message("\nWritten to ", normalizePath(OUT))
