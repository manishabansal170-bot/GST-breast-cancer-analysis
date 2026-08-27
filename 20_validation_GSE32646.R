# =============================================================================
#  20_validation_GSE32646.R
#
#  INDEPENDENT VALIDATION OF THE pCR PREDICTION
#
#  WHAT IS BEING VALIDATED
#    Script 13 found, in GSE25066 (n = 488), that GSTM1 and GSTM2 expression
#    predicted reduced pathological complete response in oestrogen-receptor-
#    positive disease (OR 0.483 and 0.508 per SD), with no effect in
#    receptor-negative disease. That rests on a single cohort.
#
#  WHY THIS COHORT AND NOT THE OBVIOUS ALTERNATIVES
#    GSE20194, GSE20271, GSE22093 and GSE23988 all draw substantially on the
#    same MD Anderson neoadjuvant series that contributes to GSE25066. Treating
#    them as independent would inflate apparent replication. GSE32646 is a
#    single-institution Japanese cohort (Osaka) with no overlap, and is
#    therefore a genuine independent test.
#
#  WHAT TO EXPECT, AND HOW TO READ IT
#    GSE32646 has 115 patients against 488. After ER stratification the
#    receptor-positive stratum will be roughly 60 to 70 patients, giving
#    perhaps 50 to 60 percent power for an odds ratio of the size reported.
#    A non-significant result is therefore ambiguous rather than refuting.
#    Read the direction and the confidence interval, not the p value alone:
#    an odds ratio pointing the same way with overlapping intervals is
#    supportive evidence even at p above 0.05.
#
#  A NOTE ON GSTP1
#    Miyake et al. generated this cohort and reported that GSTP1 protein
#    predicted poor response in receptor-NEGATIVE disease. GSTP1 is therefore
#    tested here as a check against that prior finding, reported separately,
#    and is not part of the validation. Testing a hypothesis in the cohort
#    that generated it is not independent evidence.
#
#  INPUT
#    GSE32646 series matrix, downloaded from GEO. Platform is Affymetrix
#    U133 Plus 2.0, so annotation uses hgu133plus2.db rather than hgu133a.db.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleE_drug_response")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

SERIES <- file.path(CACHE, "GSE32646_series_matrix.txt.gz")

if (!file.exists(SERIES)) {
  cat("Not found at:", SERIES, "\n\nLooking for it...\n")
  hits <- list.files(c(CACHE, WORKDIR, path.expand("~/Downloads"),
                       "C:/Users/Manisha Gupta/Downloads"),
                     pattern = "GSE32646.*series_matrix", full.names = TRUE,
                     recursive = TRUE)
  if (length(hits)) { print(hits); SERIES <- hits[1]; message("Using: ", SERIES) }
  else stop("Series matrix not found.\n",
            "Download GSE32646_series_matrix.txt.gz from:\n",
            "  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32646/matrix/\n",
            "and place it in ", CACHE)
}

pk <- c("dplyr","tidyr","tibble","ggplot2","pROC")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- c("GEOquery","hgu133plus2.db","AnnotationDbi")
miss <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(GEOquery); library(hgu133plus2.db); library(AnnotationDbi) })

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTO1","GSTO2",
         "GSTZ1","GSTK1","MGST1","MGST2","MGST3")

# The genes the validation is actually about. Everything else is exploratory
# and is reported with that status made explicit.
PRIMARY <- c("GSTM1","GSTM2")

# Discovery estimates from GSE25066, ER-positive stratum, for comparison.
DISCOVERY <- data.frame(gene = c("GSTM1","GSTM2"), OR_discovery = c(0.483, 0.508))


# =============================================================================
# 1. LOAD
# =============================================================================
message("Reading series matrix...")
eset <- getGEO(filename = SERIES, getGPL = FALSE)
ex <- Biobase::exprs(eset)
pd <- Biobase::pData(eset)
cat("Expression:", nrow(ex), "probes x", ncol(ex), "samples\n")
cat(sprintf("Value range: %.2f to %.2f\n", min(ex, na.rm = TRUE), max(ex, na.rm = TRUE)))
if (max(ex, na.rm = TRUE) > 100) {
  message("Values look linear rather than log. Applying log2(x + 1).")
  ex <- log2(ex + 1)
}


# =============================================================================
# 2. PROBE TO GENE — same rule as script 13, different annotation package
# =============================================================================
map <- AnnotationDbi::select(hgu133plus2.db, keys = rownames(ex),
                             columns = "SYMBOL", keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL), ]
hit <- map[map$SYMBOL %in% GST, ]
cat("\nProbes matching GST genes:", nrow(hit), "\n")
print(table(hit$SYMBOL))

em <- ex[hit$PROBEID, , drop = FALSE]
v <- apply(em, 1, var, na.rm = TRUE)
idx <- tapply(seq_len(nrow(hit)), hit$SYMBOL, function(i) i[which.max(v[i])])
E32 <- em[unlist(idx), , drop = FALSE]
rownames(E32) <- names(idx)

cat("\nGenes available:", paste(sort(rownames(E32)), collapse = ", "), "\n")
cat("Absent:", paste(setdiff(GST, rownames(E32)), collapse = ", "), "\n")

if (!all(PRIMARY %in% rownames(E32)))
  stop("A primary gene is absent from this platform: ",
       paste(setdiff(PRIMARY, rownames(E32)), collapse = ", "),
       ". The validation cannot proceed for that gene.")


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
  response = grab("pathologic.?response|pathological.?response|pcr|response"),
  er       = grab("^er.?status|estrogen|esr1"),
  her2     = grab("her2|erbb2"),
  stringsAsFactors = FALSE)

cat("\n===== RESPONSE LABELS AS SUPPLIED =====\n")
print(table(ph$response, useNA = "ifany"))
cat("\n===== ER LABELS AS SUPPLIED =====\n")
print(table(ph$er, useNA = "ifany"))

ph$pcr <- NA
ph$pcr[grepl("pCR|complete", ph$response, ignore.case = TRUE)] <- 1
ph$pcr[grepl("^RD|residual|non.?pCR|nonpCR", ph$response, ignore.case = TRUE)] <- 0

ph$er_pos <- NA
ph$er_pos[grepl("^P$|^pos|positive", ph$er, ignore.case = TRUE)] <- 1
ph$er_pos[grepl("^N$|^neg|negative", ph$er, ignore.case = TRUE)] <- 0

cat(sprintf("\nCoded: pCR %d | RD %d | ER+ %d | ER- %d\n",
            sum(ph$pcr == 1, na.rm = TRUE), sum(ph$pcr == 0, na.rm = TRUE),
            sum(ph$er_pos == 1, na.rm = TRUE), sum(ph$er_pos == 0, na.rm = TRUE)))

if (sum(!is.na(ph$pcr)) < 60)
  stop("Fewer than 60 responses coded. Inspect the label table printed above ",
       "and widen the grepl pattern to match the actual strings.")


# =============================================================================
# 4. ER-STRATIFIED MODEL — identical specification to script 13
# =============================================================================
d <- ph[!is.na(ph$pcr) & ph$sample %in% colnames(E32), ]
cat("\nAnalysis set:", nrow(d), "patients |",
    sprintf("pCR rate %.1f%%\n", 100 * mean(d$pcr)))

fit_stratum <- function(dd, label) {
  XX <- E32[, dd$sample, drop = FALSE]
  lapply(rownames(XX), function(g) {
    x <- as.numeric(XX[g, ]); y <- dd$pcr
    ok <- is.finite(x) & !is.na(y)
    if (sum(ok) < 25 || sd(x[ok]) == 0 || length(unique(y[ok])) < 2) return(NULL)
    fit <- glm(y[ok] ~ scale(x[ok]), family = binomial)
    s <- summary(fit)$coefficients
    est <- s[2, "Estimate"]; se <- s[2, "Std. Error"]
    data.frame(ER = label, gene = g, n = sum(ok),
               pCR_pct = round(100 * mean(y[ok]), 1),
               OR_perSD = round(exp(est), 3),
               CI_low  = round(exp(est - 1.96 * se), 3),
               CI_high = round(exp(est + 1.96 * se), 3),
               p = s[2, "Pr(>|z|)"])
  }) %>% bind_rows()
}

strat <- bind_rows(
  if (sum(d$er_pos == 1, na.rm = TRUE) >= 25)
    fit_stratum(d[which(d$er_pos == 1), ], "ER-positive"),
  if (sum(d$er_pos == 0, na.rm = TRUE) >= 25)
    fit_stratum(d[which(d$er_pos == 0), ], "ER-negative"))

if (!nrow(strat)) stop("Neither stratum has enough patients to model.")

# FDR within stratum, matching script 13
strat <- strat %>% group_by(ER) %>% mutate(FDR = p.adjust(p, "BH")) %>%
  ungroup() %>% arrange(ER, p)

cat("\n===== FULL ER-STRATIFIED RESULTS =====\n")
print(as.data.frame(strat), digits = 3)
write.csv(strat, file.path(OUT, "GSE32646_pCR_stratified.csv"), row.names = FALSE)


# =============================================================================
# 5. THE VALIDATION ITSELF
# =============================================================================
# Only GSTM1 and GSTM2 in the ER-positive stratum. Everything above is context.

cat("\n", strrep("=", 70), "\nVALIDATION: GSTM1 AND GSTM2, ER-POSITIVE STRATUM\n",
    strrep("=", 70), "\n", sep = "")

val <- strat %>%
  dplyr::filter(ER == "ER-positive", gene %in% PRIMARY) %>%
  dplyr::left_join(DISCOVERY, by = "gene") %>%
  dplyr::mutate(
    same_direction = (OR_perSD < 1) == (OR_discovery < 1),
    discovery_in_CI = OR_discovery >= CI_low & OR_discovery <= CI_high)

if (!nrow(val)) {
  cat("Neither primary gene could be modelled in the ER-positive stratum.\n")
  cat("Check the stratum size reported above.\n")
} else {
  print(as.data.frame(val[, c("gene","n","OR_perSD","CI_low","CI_high","p",
                              "OR_discovery","same_direction","discovery_in_CI")]),
        digits = 3, row.names = FALSE)

  cat("\n----- READING THIS -----\n")
  for (i in seq_len(nrow(val))) {
    r <- val[i, ]
    cat(sprintf("\n%s:\n", r$gene))
    cat(sprintf("  GSE25066 discovery OR %.3f | GSE32646 OR %.3f (%.3f to %.3f), p = %.3g\n",
                r$OR_discovery, r$OR_perSD, r$CI_low, r$CI_high, r$p))
    if (r$same_direction && r$p < 0.05) {
      cat("  REPLICATED. Same direction and significant in an independent cohort.\n")
    } else if (r$same_direction && r$discovery_in_CI) {
      cat("  CONSISTENT. Same direction, and the discovery estimate lies within\n")
      cat("  this cohort's confidence interval. Given the sample size this is\n")
      cat("  supportive, though not independently significant.\n")
    } else if (r$same_direction) {
      cat("  Same direction but the discovery estimate falls outside the interval.\n")
      cat("  Weak support; the effect here is of a different magnitude.\n")
    } else {
      cat("  NOT REPLICATED. The effect points the opposite way in this cohort.\n")
      cat("  This must be reported. The predictive claim should be softened to\n")
      cat("  exploratory, or removed from the title.\n")
    }
  }
  write.csv(val, file.path(OUT, "GSE32646_validation_primary.csv"), row.names = FALSE)
}


# =============================================================================
# 6. GSTP1 — CHECK AGAINST THE PRIOR REPORT, NOT PART OF THE VALIDATION
# =============================================================================
cat("\n", strrep("=", 70), "\nGSTP1, REPORTED SEPARATELY\n", strrep("=", 70), "\n", sep = "")
cat("Miyake et al. generated this cohort and reported that GSTP1 protein\n")
cat("predicted poor response in ER-negative disease. This is a check against\n")
cat("that finding at transcript level, not independent validation, because the\n")
cat("hypothesis and the cohort share an origin.\n\n")

gp <- strat %>% dplyr::filter(gene == "GSTP1")
if (nrow(gp)) {
  print(as.data.frame(gp[, c("ER","n","pCR_pct","OR_perSD","CI_low","CI_high","p")]),
        digits = 3, row.names = FALSE)
  neg <- gp[gp$ER == "ER-negative", ]
  if (nrow(neg)) {
    if (neg$OR_perSD < 1 && neg$p < 0.05) {
      cat("\nConsistent with the prior report: higher GSTP1 predicts lower pCR in\n")
      cat("receptor-negative disease at transcript level.\n")
    } else {
      cat("\nThe prior protein-level finding is not reproduced at transcript level\n")
      cat("here. The two are not directly comparable: transcript-protein\n")
      cat("correlation for GSTP1 is 0.733 in our proteomic analysis, and the\n")
      cat("prior study measured protein by immunohistochemistry.\n")
    }
  }
} else {
  cat("GSTP1 could not be modelled in either stratum.\n")
}


# =============================================================================
# 7. PROVENANCE
# =============================================================================
writeLines(c(
  paste("Run:", Sys.time()),
  paste("Series:", SERIES),
  sprintf("Patients with coded response: %d", nrow(d)),
  sprintf("ER-positive: %d | ER-negative: %d",
          sum(d$er_pos == 1, na.rm = TRUE), sum(d$er_pos == 0, na.rm = TRUE)),
  "Independent of GSE25066: single-institution Japanese cohort, no shared",
  "patients with the MD Anderson series.",
  "Primary validation restricted to GSTM1 and GSTM2 in the ER-positive stratum.",
  "", capture.output(sessionInfo())),
  file.path(OUT, "provenance_GSE32646.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
