# =============================================================================
#  OBJECTIVE 1 · MODULE E — DRUG SENSITIVITY AND TREATMENT RESPONSE
#
#  TWO LEVELS, AND THE SECOND ONE IS THE IMPORTANT ONE
#
#  PART A — CELL LINES (CCLE / GDSC via cBioPortal)
#    Does GST expression correlate with drug IC50 or AUC in breast cancer cell
#    lines? Direct, mechanistic, but in vitro. Cell lines drift from tumours
#    and drug exposure in a dish is not chemotherapy in a patient.
#
#  PART B — NEOADJUVANT PATIENT RESPONSE (GSE25066)
#    Does GST expression predict pathological complete response to
#    taxane-anthracycline chemotherapy in ~500 patients? This is what
#    "predictive" means in the clinical sense, and it is the analysis your
#    thesis title currently promises without evidence.
#
#    PROGNOSTIC = outcome regardless of treatment (tested, negative)
#    PREDICTIVE = response to a specific treatment (tested here)
#    These are not interchangeable and reviewers know the difference.
#
#  THE HYPOTHESIS BEING TESTED
#    GSTs conjugate glutathione to alkylating agents and anthracyclines. If
#    that confers resistance, higher GST expression should predict LOWER pCR.
#    Our subtype finding sharpens this: GSTP1 is retained in basal-like
#    disease and silenced in half of luminal tumours, so any GSTP1 effect
#    should be concentrated in basal-like/TNBC.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleE_drug_response")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2","pROC")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(lapply(pk, library, character.only = TRUE))

GST <- c("GSTA1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTM5","GSTP1",
         "GSTT2B","GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")

# Drugs that are GST substrates or otherwise mechanistically relevant.
GST_SUBSTRATES <- c("Cyclophosphamide","Chlorambucil","Melphalan","Busulfan",
                    "Doxorubicin","Epirubicin","Cisplatin","Carboplatin",
                    "Etoposide","Ifosfamide","Bleomycin","Mitoxantrone",
                    "Thiotepa","Carmustine","BCNU")


# #############################################################################
#  PART A — CELL LINE DRUG SENSITIVITY
# #############################################################################

if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("cBioPortalData", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(cBioPortalData))

cbio <- cBioPortal()

cat("\n===== SEARCHING FOR CELL LINE STUDIES WITH DRUG DATA =====\n")
studies <- getStudies(cbio)
cand <- studies$studyId[grepl("ccle|cellline|gdsc", studies$studyId, ignore.case = TRUE)]
print(cand)

CL_STUDY <- NA; CL_EXPR <- NA; CL_DRUG <- NA
for (s in cand) {
  mp <- tryCatch(molecularProfiles(cbio, s), error = function(e) NULL)
  if (is.null(mp)) next
  ids <- mp$molecularProfileId
  e <- ids[grepl("mrna", ids, ignore.case = TRUE) & !grepl("Zscore|zscore", ids)]
  d <- ids[grepl("treatment|drug|ic50|auc", ids, ignore.case = TRUE)]
  cat(sprintf("  %s: %d expression, %d drug profiles\n", s, length(e), length(d)))
  if (length(e) && length(d) && is.na(CL_STUDY)) {
    CL_STUDY <- s; CL_EXPR <- e[1]; CL_DRUG <- d[1]
  }
}

if (!is.na(CL_STUDY)) {
  cat(sprintf("\nUsing %s\n  expression: %s\n  drug: %s\n",
              CL_STUDY, CL_EXPR, CL_DRUG))

  to_mat <- function(d) {
    vcol <- intersect(c("value","alteration","score"), names(d))[1]
    gcol <- intersect(c("hugoGeneSymbol","gene","treatmentId","entityStableId"),
                      names(d))[1]
    scol <- intersect(c("sampleId","patientId"), names(d))[1]
    w <- d %>% transmute(g = .data[[gcol]], s = .data[[scol]],
                         v = suppressWarnings(as.numeric(.data[[vcol]]))) %>%
      dplyr::filter(!is.na(v)) %>% group_by(g, s) %>%
      summarise(v = mean(v), .groups = "drop") %>%
      pivot_wider(names_from = s, values_from = v)
    m <- as.matrix(w[, -1]); rownames(m) <- w$g; m
  }

  cl_res <- tryCatch({
    ex <- getDataByGenes(cbio, studyId = CL_STUDY, genes = GST,
                         by = "hugoGeneSymbol", molecularProfileIds = CL_EXPR)[[1]]
    E <- to_mat(ex)
    cat("Cell line expression:", nrow(E), "genes x", ncol(E), "lines\n")

    # Restrict to breast lines where the study permits it
    cs <- tryCatch(as.data.frame(clinicalData(cbio, CL_STUDY)), error = function(e) NULL)
    breast <- colnames(E)
    if (!is.null(cs)) {
      tcol <- grep("CANCER_TYPE|TUMOR_TYPE|LINEAGE|PRIMARY_SITE", colnames(cs),
                   value = TRUE)[1]
      if (!is.na(tcol)) {
        bl <- cs$sampleId[grepl("breast", cs[[tcol]], ignore.case = TRUE)]
        if (length(bl) > 10) {
          breast <- intersect(colnames(E), bl)
          cat("Breast cell lines:", length(breast), "\n")
        } else {
          cat("Too few breast lines identified; using all lineages.\n")
        }
      }
    }
    E <- E[, breast, drop = FALSE]

    dr <- getDataByGenes(cbio, studyId = CL_STUDY, genes = NULL,
                         molecularProfileIds = CL_DRUG)[[1]]
    D <- to_mat(dr)
    cat("Drug response:", nrow(D), "compounds x", ncol(D), "lines\n")

    sh <- intersect(colnames(E), colnames(D))
    cat("Lines with both:", length(sh), "\n")
    if (length(sh) < 15) stop("too few overlapping cell lines")

    res <- lapply(rownames(E), function(g) {
      lapply(rownames(D), function(d) {
        x <- E[g, sh]; y <- D[d, sh]
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 15 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NULL)
        ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
        data.frame(gene = g, drug = d, n = sum(ok),
                   rho = round(unname(ct$estimate), 3), p = ct$p.value)
      }) %>% bind_rows()
    }) %>% bind_rows()

    if (nrow(res)) {
      res <- res %>% mutate(FDR = p.adjust(p, "BH")) %>% arrange(p)
      cat("\n===== TOP GST-DRUG ASSOCIATIONS =====\n")
      print(head(as.data.frame(res), 30))
      write.csv(res, file.path(OUT, "celline_gst_drug_correlation.csv"),
                row.names = FALSE)

      sub <- res %>% dplyr::filter(grepl(paste(GST_SUBSTRATES, collapse = "|"),
                                 drug, ignore.case = TRUE))
      if (nrow(sub)) {
        cat("\n===== GST-SUBSTRATE DRUGS SPECIFICALLY =====\n")
        print(head(as.data.frame(sub %>% arrange(p)), 25))
        write.csv(sub, file.path(OUT, "celline_gst_substrate_drugs.csv"),
                  row.names = FALSE)
      }
      cat("\nNOTE on direction: if the drug metric is IC50 or AUC, HIGHER means\n")
      cat("MORE resistant, so a POSITIVE rho supports GST-mediated resistance.\n")
      cat("Confirm which metric this profile reports before interpreting.\n")
    }
    TRUE
  }, error = function(e) { message("Part A failed: ", conditionMessage(e)); FALSE })
} else {
  cat("\nNo cell line study with both expression and drug data found.\n")
  cat("Alternative: download GDSC2 fitted dose response from\n")
  cat("cancerrxgene.org/downloads together with the RMA expression matrix.\n")
}


# #############################################################################
#  PART B — NEOADJUVANT RESPONSE  (the one that matters for your title)
# #############################################################################

cat("\n\n#############################################################\n")
cat("PART B - GSE25066 NEOADJUVANT CHEMOTHERAPY RESPONSE\n")
cat("#############################################################\n")
cat("~500 patients, taxane-anthracycline neoadjuvant therapy, with\n")
cat("pathological complete response (pCR) versus residual disease (RD).\n")
cat("Affymetrix U133A. Download is modest.\n\n")

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(GEOquery))

gf <- file.path(CACHE, "gse25066.rds")
if (file.exists(gf)) {
  obj <- readRDS(gf); E25 <- obj$expr; ph25 <- obj$pheno
  message("Loaded GSE25066 from cache.")
} else {
  message("Fetching GSE25066 (several minutes)...")
  gse <- getGEO("GSE25066", GSEMatrix = TRUE, getGPL = TRUE)
  eset <- gse[[1]]

  ex <- Biobase::exprs(eset)
  fd <- Biobase::fData(eset)
  pd <- Biobase::pData(eset)

  sym_col <- grep("Gene.?symbol|GENE_SYMBOL", colnames(fd),
                  value = TRUE, ignore.case = TRUE)[1]
  cat("Annotation column:", sym_col, "\n")

  sym <- as.character(fd[[sym_col]])
  sym <- sub("\\s*///.*$", "", sym)          # first symbol of multi-mappers

  keep <- sym %in% GST & !is.na(sym)
  cat("Probes matching GST genes:", sum(keep), "\n")

  em <- ex[keep, , drop = FALSE]
  sm <- sym[keep]

  # One probe per gene: highest variance, the standard choice, because a flat
  # probe is usually background rather than a second real measurement.
  v <- apply(em, 1, var, na.rm = TRUE)
  idx <- tapply(seq_along(sm), sm, function(i) i[which.max(v[i])])
  E25 <- em[unlist(idx), , drop = FALSE]
  rownames(E25) <- names(idx)

  ch <- pd[, grep("characteristics_ch1", colnames(pd)), drop = FALSE]
  grab <- function(pat) apply(ch, 1, function(r) {
    h <- grep(pat, r, value = TRUE, ignore.case = TRUE)
    if (!length(h)) return(NA_character_)
    trimws(sub("^[^:]*:\\s*", "", h[1])) })

  ph25 <- data.frame(
    sample   = rownames(pd),
    response = grab("pathologic.?response|pcr"),
    er       = grab("^esr1|er.?status|er_status"),
    her2     = grab("erbb2|her2"),
    subtype  = grab("pam50|subtype"),
    grade    = grab("grade"),
    stringsAsFactors = FALSE)

  saveRDS(list(expr = E25, pheno = ph25), gf)
  message("Cached.")
}

cat("\nGST genes on the array:", paste(sort(rownames(E25)), collapse = ", "), "\n")
cat("Absent (no U133A probe):",
    paste(setdiff(GST, rownames(E25)), collapse = ", "), "\n")

cat("\n===== RESPONSE LABELS AS SUPPLIED =====\n")
print(table(ph25$response, useNA = "ifany"))

ph25$pcr <- NA
ph25$pcr[grepl("^pCR|complete", ph25$response, ignore.case = TRUE)] <- 1
ph25$pcr[grepl("^RD|residual", ph25$response, ignore.case = TRUE)]  <- 0

cat("\nCoded: pCR =", sum(ph25$pcr == 1, na.rm = TRUE),
    "| RD =", sum(ph25$pcr == 0, na.rm = TRUE), "\n")
if (sum(!is.na(ph25$pcr)) < 100)
  warning("Few coded responses - inspect the label strings above and adjust ",
          "the grepl patterns.")

ph25$er_pos <- NA
ph25$er_pos[grepl("^P|pos", ph25$er, ignore.case = TRUE)] <- 1
ph25$er_pos[grepl("^N|neg", ph25$er, ignore.case = TRUE)] <- 0
cat("ER positive:", sum(ph25$er_pos == 1, na.rm = TRUE),
    "| ER negative:", sum(ph25$er_pos == 0, na.rm = TRUE), "\n")


# ---------------------------------------------------------------------------
#  B1. Does GST expression predict pCR?
# ---------------------------------------------------------------------------
d <- ph25[!is.na(ph25$pcr) & ph25$sample %in% colnames(E25), ]
X <- E25[, d$sample, drop = FALSE]
cat("\nAnalysis set:", nrow(d), "patients\n")

pred <- lapply(rownames(X), function(g) {
  x <- as.numeric(X[g, ]); y <- d$pcr
  ok <- is.finite(x) & !is.na(y)
  if (sum(ok) < 50 || sd(x[ok]) == 0) return(NULL)
  fit <- glm(y[ok] ~ scale(x[ok]), family = binomial)
  s <- summary(fit)$coefficients
  auc <- suppressMessages(as.numeric(pROC::auc(pROC::roc(y[ok], x[ok], quiet = TRUE))))
  data.frame(gene = g, n = sum(ok),
             OR_perSD = round(exp(s[2, "Estimate"]), 3),
             CI_low = round(exp(s[2, "Estimate"] - 1.96 * s[2, "Std. Error"]), 3),
             CI_high = round(exp(s[2, "Estimate"] + 1.96 * s[2, "Std. Error"]), 3),
             p = s[2, "Pr(>|z|)"], AUC = round(auc, 3))
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "BH")) %>% arrange(p)

cat("\n===== GST EXPRESSION AND pCR (all patients) =====\n")
print(as.data.frame(pred), digits = 3)
cat("\nOR below 1 per SD means higher expression predicts LOWER pCR, i.e.\n")
cat("resistance - the direction the GST hypothesis predicts.\n")
write.csv(pred, file.path(OUT, "GSE25066_pCR_prediction.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
#  B2. Stratified by ER status — where our subtype finding predicts an effect
# ---------------------------------------------------------------------------
# GSTP1 is retained in basal-like/TNBC and silenced in roughly half of luminal
# tumours. If it mediates resistance, the effect should be concentrated in
# ER-negative disease. This is the specific, falsifiable prediction our
# methylation result generates.

if (sum(!is.na(d$er_pos)) > 100) {
  strat <- lapply(c(0, 1), function(erv) {
    dd <- d[which(d$er_pos == erv), ]
    if (nrow(dd) < 50) return(NULL)
    XX <- E25[, dd$sample, drop = FALSE]
    lapply(rownames(XX), function(g) {
      x <- as.numeric(XX[g, ]); y <- dd$pcr
      ok <- is.finite(x) & !is.na(y)
      if (sum(ok) < 30 || sd(x[ok]) == 0 || length(unique(y[ok])) < 2) return(NULL)
      fit <- glm(y[ok] ~ scale(x[ok]), family = binomial)
      s <- summary(fit)$coefficients
      data.frame(ER = ifelse(erv == 1, "ER-positive", "ER-negative"),
                 gene = g, n = sum(ok),
                 pCR_rate = round(100 * mean(y[ok]), 1),
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
      cat("\n--- GSTP1, the gene carrying the therapeutic argument ---\n")
      print(as.data.frame(g1), digits = 3)
      cat("\nOur prediction: the effect should be present in ER-negative disease,\n")
      cat("where GSTP1 is unmethylated and highly expressed, and absent in\n")
      cat("ER-positive disease, where roughly half of tumours have silenced it.\n")
    }
  }
}


# ---------------------------------------------------------------------------
#  B3. Figures
# ---------------------------------------------------------------------------
plot_df <- as.data.frame(t(X)) %>% rownames_to_column("sample") %>%
  left_join(d[, c("sample","pcr","er_pos")], by = "sample") %>%
  pivot_longer(-c(sample, pcr, er_pos), names_to = "gene", values_to = "expr") %>%
  dplyr::filter(!is.na(pcr), is.finite(expr)) %>%
  mutate(Response = factor(ifelse(pcr == 1, "pCR", "Residual disease"),
                           levels = c("Residual disease", "pCR")))

p1 <- ggplot(plot_df, aes(Response, expr, fill = Response)) +
  geom_violin(scale = "width", alpha = .3, colour = NA) +
  geom_boxplot(width = .28, outlier.size = .3, linewidth = .3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Residual disease" = "#C0392B", "pCR" = "#2E86AB"),
                    guide = "none") +
  labs(x = NULL, y = "Expression (log2)",
       title = "GST expression and response to neoadjuvant chemotherapy",
       subtitle = paste0("GSE25066, n = ", nrow(d),
                         ", taxane-anthracycline; pCR vs residual disease")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "GSE25066_pCR_boxplots.png"), p1,
       width = 10, height = 8, dpi = 300)

if ("GSTP1" %in% rownames(X) && sum(!is.na(d$er_pos)) > 100) {
  pg <- plot_df %>% dplyr::filter(gene == "GSTP1", !is.na(er_pos)) %>%
    mutate(ER = ifelse(er_pos == 1, "ER-positive", "ER-negative"))
  p2 <- ggplot(pg, aes(Response, expr, fill = Response)) +
    geom_violin(scale = "width", alpha = .3, colour = NA) +
    geom_boxplot(width = .3, outlier.size = .4, linewidth = .35) +
    facet_wrap(~ ER) +
    scale_fill_manual(values = c("Residual disease" = "#C0392B", "pCR" = "#2E86AB"),
                      guide = "none") +
    labs(x = NULL, y = "GSTP1 expression (log2)",
         title = "GSTP1 and neoadjuvant response, stratified by ER status",
         subtitle = "Prediction: effect present in ER-negative disease, absent in ER-positive") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, "GSE25066_GSTP1_by_ER.png"), p2,
         width = 8, height = 4.5, dpi = 300)
}


# =============================================================================
#  HOW TO READ THIS
# =============================================================================
cat("\n============================================================\n")
cat("INTERPRETING PART B\n")
cat("============================================================\n")
cat("An odds ratio below 1 per SD means higher expression predicts LOWER pCR,\n")
cat("which is resistance - the direction the GST hypothesis predicts.\n\n")
cat("The result that would matter most: GSTP1 significant in ER-negative and\n")
cat("null in ER-positive. That is the specific prediction generated by the\n")
cat("methylation finding, and confirming it would let you write 'predictive'\n")
cat("in the title with evidence behind it.\n\n")
cat("If nothing reaches significance, that is also reportable. It would mean\n")
cat("GST expression does not predict chemotherapy response at the transcript\n")
cat("level in this cohort - narrowing your claim to mechanism rather than\n")
cat("clinical prediction, which is where the rest of the evidence points.\n\n")
cat("CAVEATS to state either way:\n")
cat("  - U133A is an older array with compressed dynamic range\n")
cat("  - pCR after taxane-anthracycline is far commoner in ER-negative\n")
cat("    disease, so ER status confounds any unstratified analysis\n")
cat("  - one regimen only; results may not generalise to platinum or\n")
cat("    cyclophosphamide-based therapy\n")
cat("============================================================\n")

writeLines(c(
  paste("Run:", Sys.time()),
  "Part A: cell line drug sensitivity via cBioPortal.",
  paste("  Study:", CL_STUDY, "| expression:", CL_EXPR, "| drug:", CL_DRUG),
  "Part B: GSE25066 neoadjuvant taxane-anthracycline cohort, pCR vs residual",
  "disease. One probe per gene selected by highest variance. Logistic",
  "regression on expression scaled per standard deviation, stratified by ER",
  "status because pCR rates differ substantially between ER strata.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleE.txt"))

message("\nWritten to ", normalizePath(OUT))
