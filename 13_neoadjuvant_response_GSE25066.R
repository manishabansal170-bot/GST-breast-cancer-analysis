# =============================================================================
#  MODULE E · PART B — GSE25066 NEOADJUVANT RESPONSE
#  From a locally downloaded series matrix
#
#  Probe-to-gene mapping uses the hgu133a.db annotation package rather than the
#  GPL96 file from GEO. It installs once from Bioconductor, works offline
#  thereafter, and is more reliable than hand-transcribed probe identifiers.
#
#  WHAT THIS TESTS
#    Prognostic = outcome regardless of treatment. Already tested; negative.
#    Predictive = response to a specific treatment. Tested here.
#    Your thesis title claims both. This is the analysis that can support the
#    second one.
#
#  THE SPECIFIC PREDICTION
#    GSTP1 is unmethylated and highly expressed in basal-like/TNBC, and
#    silenced in roughly half of luminal tumours. If it mediates resistance,
#    the effect on pCR should appear in ER-negative disease and be absent in
#    ER-positive. Section B2 tests exactly that.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleE_drug_response")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

# ---- EDIT THIS to wherever the file landed --------------------------------
SERIES <- file.path(CACHE, "GSE25066_series_matrix.txt.gz")

if (!file.exists(SERIES)) {
  cat("Not found at:", SERIES, "\n\nLooking for it...\n")
  hits <- list.files(c(CACHE, WORKDIR, path.expand("~/Downloads"),
                       "C:/Users/Manisha Gupta/Downloads"),
                     pattern = "GSE25066.*series_matrix", full.names = TRUE,
                     recursive = TRUE)
  if (length(hits)) { print(hits); SERIES <- hits[1]
                      message("Using: ", SERIES) }
  else stop("Series matrix not found. Set SERIES to its full path.")
}

pk <- c("dplyr","tidyr","tibble","ggplot2","pROC")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- c("GEOquery","hgu133a.db","AnnotationDbi")
miss <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(GEOquery); library(hgu133a.db); library(AnnotationDbi) })

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTO1","GSTO2",
         "GSTZ1","GSTK1","MGST1","MGST2","MGST3")


# =============================================================================
# 1. LOAD
# =============================================================================
message("Reading series matrix...")
eset <- getGEO(filename = SERIES, getGPL = FALSE)

ex <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)
cat("Expression:", nrow(ex), "probes x", ncol(ex), "samples\n")

# Values should already be log-scale. If the maximum is large, they are not.
cat(sprintf("Value range: %.2f to %.2f\n", min(ex, na.rm = TRUE), max(ex, na.rm = TRUE)))
if (max(ex, na.rm = TRUE) > 100) {
  message("Values look linear rather than log. Applying log2(x + 1).")
  ex <- log2(ex + 1)
}


# =============================================================================
# 2. PROBE TO GENE
# =============================================================================
map <- AnnotationDbi::select(hgu133a.db, keys = rownames(ex),
                             columns = "SYMBOL", keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL), ]

hit <- map[map$SYMBOL %in% GST, ]
cat("\nProbes matching GST genes:", nrow(hit), "\n")
print(table(hit$SYMBOL))

em <- ex[hit$PROBEID, , drop = FALSE]

# One probe per gene: highest variance. A flat probe is background, not a
# second measurement.
v <- apply(em, 1, var, na.rm = TRUE)
idx <- tapply(seq_len(nrow(hit)), hit$SYMBOL, function(i) i[which.max(v[i])])
E25 <- em[unlist(idx), , drop = FALSE]
rownames(E25) <- names(idx)

cat("\nGenes available:", paste(sort(rownames(E25)), collapse = ", "), "\n")
cat("Absent from U133A:", paste(setdiff(GST, rownames(E25)), collapse = ", "), "\n")


# =============================================================================
# 3. PHENOTYPE
# =============================================================================
ch <- pd[, grep("characteristics_ch1", colnames(pd)), drop = FALSE]
cat("\n===== EXAMPLE CHARACTERISTICS (first sample) =====\n")
print(unname(unlist(ch[1, ])))

grab <- function(pat) apply(ch, 1, function(r) {
  h <- grep(pat, r, value = TRUE, ignore.case = TRUE)
  if (!length(h)) return(NA_character_)
  trimws(sub("^[^:]*:\\s*", "", h[1])) })

ph <- data.frame(
  sample   = rownames(pd),
  response = grab("pathologic.?response|pcr.?rd|^response"),
  er       = grab("^er.?status|esr1.?status|er_status"),
  her2     = grab("her2|erbb2"),
  grade    = grab("grade"),
  stringsAsFactors = FALSE)

cat("\n===== RESPONSE LABELS AS SUPPLIED =====\n")
print(table(ph$response, useNA = "ifany"))
cat("\n===== ER LABELS AS SUPPLIED =====\n")
print(table(ph$er, useNA = "ifany"))

ph$pcr <- NA
ph$pcr[grepl("pCR|complete", ph$response, ignore.case = TRUE)] <- 1
ph$pcr[grepl("^RD|residual", ph$response, ignore.case = TRUE)]  <- 0

ph$er_pos <- NA
ph$er_pos[grepl("^P$|^pos|positive", ph$er, ignore.case = TRUE)] <- 1
ph$er_pos[grepl("^N$|^neg|negative", ph$er, ignore.case = TRUE)] <- 0

cat(sprintf("\nCoded: pCR %d | RD %d | ER+ %d | ER- %d\n",
            sum(ph$pcr == 1, na.rm = TRUE), sum(ph$pcr == 0, na.rm = TRUE),
            sum(ph$er_pos == 1, na.rm = TRUE), sum(ph$er_pos == 0, na.rm = TRUE)))

if (sum(!is.na(ph$pcr)) < 100)
  stop("Fewer than 100 responses coded. Look at the label table printed above ",
       "and widen the grepl pattern for ph$pcr to match the actual strings.")


# =============================================================================
# 4. DOES GST EXPRESSION PREDICT pCR?
# =============================================================================
d <- ph[!is.na(ph$pcr) & ph$sample %in% colnames(E25), ]
X <- E25[, d$sample, drop = FALSE]
cat("\nAnalysis set:", nrow(d), "patients |",
    sprintf("pCR rate %.1f%%\n", 100 * mean(d$pcr)))

pred <- lapply(rownames(X), function(g) {
  x <- as.numeric(X[g, ]); y <- d$pcr
  ok <- is.finite(x) & !is.na(y)
  if (sum(ok) < 50 || sd(x[ok]) == 0) return(NULL)
  fit <- glm(y[ok] ~ scale(x[ok]), family = binomial)
  s <- summary(fit)$coefficients
  auc <- suppressMessages(as.numeric(pROC::auc(pROC::roc(y[ok], x[ok], quiet = TRUE))))
  data.frame(gene = g, n = sum(ok),
             OR_perSD = round(exp(s[2, "Estimate"]), 3),
             CI_low  = round(exp(s[2, "Estimate"] - 1.96 * s[2, "Std. Error"]), 3),
             CI_high = round(exp(s[2, "Estimate"] + 1.96 * s[2, "Std. Error"]), 3),
             p = s[2, "Pr(>|z|)"], AUC = round(auc, 3))
}) %>% bind_rows() %>% mutate(FDR = p.adjust(p, "BH")) %>% arrange(p)

cat("\n===== GST EXPRESSION AND pCR, ALL PATIENTS =====\n")
print(as.data.frame(pred), digits = 3)
cat("\nOR below 1 per SD = higher expression predicts LOWER pCR, i.e.\n")
cat("resistance, which is the direction the GST hypothesis predicts.\n")
write.csv(pred, file.path(OUT, "GSE25066_pCR_prediction.csv"), row.names = FALSE)


# =============================================================================
# 5. STRATIFIED BY ER — the decisive test
# =============================================================================
# pCR rates differ substantially between ER strata, so ER confounds any
# unstratified analysis. More importantly, our methylation finding predicts
# the GSTP1 effect should be confined to ER-negative disease.

if (sum(!is.na(d$er_pos)) > 100) {
  strat <- lapply(c(0, 1), function(erv) {
    dd <- d[which(d$er_pos == erv), ]
    if (nrow(dd) < 40) return(NULL)
    XX <- E25[, dd$sample, drop = FALSE]
    lapply(rownames(XX), function(g) {
      x <- as.numeric(XX[g, ]); y <- dd$pcr
      ok <- is.finite(x) & !is.na(y)
      if (sum(ok) < 30 || sd(x[ok]) == 0 || length(unique(y[ok])) < 2) return(NULL)
      fit <- glm(y[ok] ~ scale(x[ok]), family = binomial)
      s <- summary(fit)$coefficients
      data.frame(ER = ifelse(erv == 1, "ER-positive", "ER-negative"),
                 gene = g, n = sum(ok), pCR_pct = round(100 * mean(y[ok]), 1),
                 OR_perSD = round(exp(s[2, "Estimate"]), 3),
                 p = s[2, "Pr(>|z|)"])
    }) %>% bind_rows()
  }) %>% bind_rows()

  if (nrow(strat)) {
    strat <- strat %>% group_by(ER) %>% mutate(FDR = p.adjust(p, "BH")) %>%
      ungroup() %>% arrange(ER, p)
    cat("\n===== STRATIFIED BY ER STATUS =====\n")
    print(as.data.frame(strat), digits = 3)
    write.csv(strat, file.path(OUT, "GSE25066_pCR_by_ER.csv"), row.names = FALSE)

    g1 <- strat[strat$gene == "GSTP1", ]
    if (nrow(g1)) {
      cat("\n--- GSTP1 ---\n"); print(as.data.frame(g1), digits = 3)
      cat("\nPrediction from the methylation data: significant in ER-negative,\n")
      cat("null in ER-positive. If that holds, the chain runs from promoter\n")
      cat("methylation through expression to clinical response.\n")
    }
  }
}


# =============================================================================
# 6. FIGURES
# =============================================================================
pdf_df <- as.data.frame(t(X)) %>% rownames_to_column("sample") %>%
  left_join(d[, c("sample","pcr","er_pos")], by = "sample") %>%
  pivot_longer(-c(sample, pcr, er_pos), names_to = "gene", values_to = "expr") %>%
  dplyr::filter(!is.na(pcr), is.finite(expr)) %>%
  mutate(Response = factor(ifelse(pcr == 1, "pCR", "Residual disease"),
                           levels = c("Residual disease", "pCR")))

p1 <- ggplot(pdf_df, aes(Response, expr, fill = Response)) +
  geom_violin(scale = "width", alpha = .3, colour = NA) +
  geom_boxplot(width = .28, outlier.size = .3, linewidth = .3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Residual disease" = "#C0392B", "pCR" = "#2E86AB"),
                    guide = "none") +
  labs(x = NULL, y = "Expression (log2)",
       title = "GST expression and response to neoadjuvant chemotherapy",
       subtitle = paste0("GSE25066, n = ", nrow(d), ", taxane-anthracycline")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "GSE25066_pCR_boxplots.png"), p1, width = 10, height = 8, dpi = 300)

if ("GSTP1" %in% rownames(X) && sum(!is.na(d$er_pos)) > 100) {
  pg <- pdf_df %>% dplyr::filter(gene == "GSTP1", !is.na(er_pos)) %>%
    mutate(ER = ifelse(er_pos == 1, "ER-positive", "ER-negative"))
  p2 <- ggplot(pg, aes(Response, expr, fill = Response)) +
    geom_violin(scale = "width", alpha = .3, colour = NA) +
    geom_boxplot(width = .3, outlier.size = .4, linewidth = .35) +
    facet_wrap(~ ER) +
    scale_fill_manual(values = c("Residual disease" = "#C0392B", "pCR" = "#2E86AB"),
                      guide = "none") +
    labs(x = NULL, y = "GSTP1 expression (log2)",
         title = "GSTP1 and neoadjuvant response by ER status",
         subtitle = "Prediction: effect in ER-negative disease, absent in ER-positive") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, "GSE25066_GSTP1_by_ER.png"), p2, width = 8, height = 4.5, dpi = 300)
}

saveRDS(list(expr = E25, pheno = ph), file.path(CACHE, "gse25066.rds"))

writeLines(c(
  paste("Run:", Sys.time()),
  "GSE25066 neoadjuvant taxane-anthracycline cohort, pCR vs residual disease.",
  "Probe-to-gene mapping via hgu133a.db; one probe per gene by highest variance.",
  "Logistic regression on expression scaled per SD, stratified by ER status",
  "because pCR rates differ substantially between ER strata.",
  paste("Patients analysed:", nrow(d)),
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleE_partB.txt"))

message("\nWritten to ", normalizePath(OUT))
