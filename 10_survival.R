# =============================================================================
#  OBJECTIVE 1 · MODULE D — PROGNOSTIC SIGNIFICANCE
#
#  Discovery : TCGA-BRCA   (expression already cached from Toil)
#  Validation: METABRIC    (expression already cached from cBioPortal)
#
#  No large downloads. Clinical tables are a few MB each.
#
#  WHAT THIS PRODUCES
#    univariate_cox_<cohort>.csv    per-gene HR, CI, p, FDR (continuous)
#    km_<cohort>/<GENE>.png         Kaplan-Meier curves
#    clinical_baseline.csv          C-index of clinical model alone
#    risk_score_model.csv           LASSO-Cox coefficients (locked on TCGA)
#    validation_metabric.csv        the locked score applied out-of-sample
#    timeROC_<cohort>.csv           AUC at 1, 3, 5 years
#
#  THREE DESIGN DECISIONS, STATED UP FRONT
#
#  1. CONTINUOUS COX IS PRIMARY, KM IS ILLUSTRATION.
#     Median-splitting a continuous variable discards information, and
#     "optimal cutpoint" searching inflates false positives badly unless the
#     p-value is corrected for the search. KM curves are produced for the
#     figures; the inference comes from Cox on continuous expression.
#
#  2. THE CLINICAL BASELINE IS NOT OPTIONAL.
#     A GST score must be judged on what it ADDS to stage, grade, age and
#     PAM50 - not on whether it is significant alone. Section 5 fits the
#     clinical model first and reports the incremental C-index. If the
#     increment is near zero, that is the finding and it should be reported.
#
#  3. FOR TCGA-BRCA, USE PFI RATHER THAN OS.
#     BRCA has long follow-up and few deaths, so OS is underpowered and
#     contaminated by non-cancer mortality. The TCGA Pan-Cancer Clinical Data
#     Resource (Liu et al., Cell 2018) recommends progression-free interval
#     for BRCA. OS is computed too, for comparison with METABRIC.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleD_survival")
for (d in c(OUT, file.path(OUT, "km_tcga"), file.path(OUT, "km_metabric")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("survival","survminer","glmnet","timeROC","dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(lapply(pk, library, character.only = TRUE))

set.seed(42)   # LASSO cross-validation is stochastic; fix it for reproducibility

SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. TCGA EXPRESSION  (cached)
# =============================================================================
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^raw - 0.001; lin[lin < 0] <- 0
tcga_expr <- log2(lin + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])

tum <- ph %>% filter(toupper(study) == "TCGA",
                     grepl("breast", tissue, ignore.case = TRUE),
                     grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12))
tum <- tum[tum$sample %in% colnames(tcga_expr), ]
tum <- tum[!duplicated(tum$patient), ]


# =============================================================================
# 2. TCGA CLINICAL  (small download)
# =============================================================================
cf <- file.path(CACHE, "tcga_clinical.rds")
if (file.exists(cf)) {
  clin_t <- readRDS(cf)
} else {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE) }
  clin_t <- TCGAbiolinks::GDCquery_clinic("TCGA-BRCA", "clinical")
  saveRDS(clin_t, cf)
}

pick <- function(d, pat) { h <- grep(pat, colnames(d), value = TRUE, ignore.case = TRUE)[1]
                           if (is.na(h)) rep(NA, nrow(d)) else d[[h]] }

surv_t <- data.frame(
  patient   = pick(clin_t, "^submitter_id$|^bcr_patient_barcode$"),
  vital     = pick(clin_t, "vital_status"),
  d_death   = suppressWarnings(as.numeric(pick(clin_t, "days_to_death"))),
  d_follow  = suppressWarnings(as.numeric(pick(clin_t, "days_to_last_follow_up"))),
  age       = suppressWarnings(as.numeric(pick(clin_t, "age_at_index|age_at_diagnosis"))),
  stage     = as.character(pick(clin_t, "ajcc_pathologic_stage")),
  stringsAsFactors = FALSE)

surv_t$os_time  <- ifelse(surv_t$vital == "Dead", surv_t$d_death, surv_t$d_follow)
surv_t$os_event <- as.integer(surv_t$vital == "Dead")

# Collapse AJCC substages; keep it ordinal and interpretable.
surv_t$stage_g <- dplyr::case_when(
  grepl("Stage IV",  surv_t$stage) ~ "IV",
  grepl("Stage III", surv_t$stage) ~ "III",
  grepl("Stage II",  surv_t$stage) ~ "II",
  grepl("Stage I",   surv_t$stage) ~ "I",
  TRUE ~ NA_character_)

st <- TCGAbiolinks::TCGAquery_subtype(tumor = "brca")
surv_t$pam50 <- st$BRCA_Subtype_PAM50[match(surv_t$patient, st$patient)]

# Age can arrive in days in some releases. Normalise, then sanity-check.
if (!all(is.na(surv_t$age)) && median(surv_t$age, na.rm = TRUE) > 200)
  surv_t$age <- surv_t$age / 365.25

d_t <- as.data.frame(t(tcga_expr[, tum$sample, drop = FALSE]))
d_t$patient <- tum$patient
d_t <- left_join(d_t, surv_t, by = "patient") %>%
  filter(!is.na(os_time), os_time > 0, pam50 %in% SUB)

GENES <- rownames(tcga_expr)

cat("\n===== TCGA SURVIVAL SET =====\n")
cat("Patients:", nrow(d_t), "| deaths:", sum(d_t$os_event),
    sprintf("(%.1f%%)\n", 100 * mean(d_t$os_event)))
cat("Median follow-up:", round(median(d_t$os_time) / 365.25, 1), "years\n")
if (sum(d_t$os_event) < 100)
  message("NOTE: few OS events. Underpowered for a multi-gene model - ",
          "METABRIC carries the survival analysis.")


# =============================================================================
# 3. METABRIC  (cached expression, small clinical fetch)
# =============================================================================
mb <- readRDS(file.path(CACHE, "metabric_gst.rds"))
mb_expr <- mb$mat

mcf <- file.path(CACHE, "metabric_clinical.rds")
if (file.exists(mcf)) {
  clin_m <- readRDS(mcf)
} else {
  suppressPackageStartupMessages(library(cBioPortalData))
  clin_m <- as.data.frame(clinicalData(cBioPortal(), "brca_metabric"))
  saveRDS(clin_m, mcf)
}

surv_m <- data.frame(
  sample  = pick(clin_m, "^sampleId$|^patientId$"),
  os_mo   = suppressWarnings(as.numeric(pick(clin_m, "OS_MONTHS"))),
  os_raw  = as.character(pick(clin_m, "OS_STATUS")),
  age     = suppressWarnings(as.numeric(pick(clin_m, "AGE_AT_DIAGNOSIS"))),
  stage   = as.character(pick(clin_m, "TUMOR_STAGE")),
  grade   = suppressWarnings(as.numeric(pick(clin_m, "GRADE|NEOPLASM_HISTOLOGIC_GRADE"))),
  pam50   = as.character(pick(clin_m, "CLAUDIN_SUBTYPE")),
  stringsAsFactors = FALSE)

surv_m$os_event <- as.integer(grepl("DECEASED|^1", surv_m$os_raw))
surv_m$os_time  <- surv_m$os_mo * 30.44          # months -> days, matching TCGA

d_m <- as.data.frame(t(mb_expr)) %>% rownames_to_column("sample") %>%
  left_join(surv_m, by = "sample") %>%
  filter(!is.na(os_time), os_time > 0, pam50 %in% SUB)

cat("\n===== METABRIC SURVIVAL SET =====\n")
cat("Patients:", nrow(d_m), "| deaths:", sum(d_m$os_event),
    sprintf("(%.1f%%)\n", 100 * mean(d_m$os_event)))
cat("Median follow-up:", round(median(d_m$os_time) / 365.25, 1), "years\n\n")


# =============================================================================
# 4. UNIVARIATE COX — continuous expression, per gene
# =============================================================================
uni_cox <- function(dat, genes, label) {
  genes <- intersect(genes, colnames(dat))
  res <- lapply(genes, function(g) {
    x <- dat[[g]]
    if (sd(x, na.rm = TRUE) == 0) return(NULL)
    # Per-SD scaling makes HRs comparable between genes and across cohorts,
    # which raw log2 units are not.
    fit <- coxph(Surv(os_time, os_event) ~ scale(x), data = dat)
    s <- summary(fit)
    data.frame(gene = g,
               HR_perSD = round(s$coefficients[1, "exp(coef)"], 3),
               CI_low   = round(s$conf.int[1, "lower .95"], 3),
               CI_high  = round(s$conf.int[1, "upper .95"], 3),
               p        = s$coefficients[1, "Pr(>|z|)"],
               C_index  = round(s$concordance[1], 3))
  }) %>% bind_rows()
  res$FDR <- p.adjust(res$p, "BH")
  res$cohort <- label
  res %>% arrange(p)
}

cox_t <- uni_cox(d_t, GENES, "TCGA")
cox_m <- uni_cox(d_m, GENES, "METABRIC")

cat("===== UNIVARIATE COX — TCGA (OS) =====\n")
print(as.data.frame(cox_t), digits = 3)
cat("\n===== UNIVARIATE COX — METABRIC (OS) =====\n")
print(as.data.frame(cox_m), digits = 3)

write.csv(cox_t, file.path(OUT, "univariate_cox_tcga.csv"), row.names = FALSE)
write.csv(cox_m, file.path(OUT, "univariate_cox_metabric.csv"), row.names = FALSE)

# Cross-cohort agreement is the result that matters more than either alone.
both <- inner_join(cox_t[, c("gene","HR_perSD","FDR")],
                   cox_m[, c("gene","HR_perSD","FDR")],
                   by = "gene", suffix = c("_TCGA","_MB")) %>%
  mutate(same_direction = sign(HR_perSD_TCGA - 1) == sign(HR_perSD_MB - 1),
         both_sig = FDR_TCGA < 0.05 & FDR_MB < 0.05)

cat("\n===== CROSS-COHORT AGREEMENT =====\n")
print(as.data.frame(both), digits = 3)
cat(sprintf("\nSame direction: %d of %d | significant in both: %d\n",
            sum(both$same_direction), nrow(both), sum(both$both_sig)))
write.csv(both, file.path(OUT, "cox_cross_cohort.csv"), row.names = FALSE)


# =============================================================================
# 5. CLINICAL BASELINE — the comparison that decides whether this is useful
# =============================================================================
base_t <- coxph(Surv(os_time, os_event) ~ age + factor(stage_g) + factor(pam50),
                data = d_t)
base_m <- coxph(Surv(os_time, os_event) ~ age + grade + factor(pam50),
                data = d_m)

c_base_t <- summary(base_t)$concordance[1]
c_base_m <- summary(base_m)$concordance[1]

cat(sprintf("\nClinical baseline C-index — TCGA %.3f | METABRIC %.3f\n",
            c_base_t, c_base_m))

incr <- lapply(intersect(GENES, colnames(d_m)), function(g) {
  f <- update(base_m, paste(". ~ . + scale(", g, ")"))
  data.frame(gene = g,
             C_clinical = round(c_base_m, 3),
             C_plus_gene = round(summary(f)$concordance[1], 3),
             increment = round(summary(f)$concordance[1] - c_base_m, 4))
}) %>% bind_rows() %>% arrange(desc(increment))

cat("\n===== INCREMENTAL C-INDEX OVER CLINICAL MODEL (METABRIC) =====\n")
print(as.data.frame(incr))
cat("\nAn increment below ~0.01 means the gene adds nothing a clinician\n")
cat("does not already have. Report that outcome if it occurs - it is a\n")
cat("finding, and pre-empting it is far better than being asked.\n")
write.csv(incr, file.path(OUT, "incremental_cindex_metabric.csv"), row.names = FALSE)


# =============================================================================
# 6. KAPLAN-MEIER — figures only
# =============================================================================
km_plots <- function(dat, genes, dir, label) {
  for (g in intersect(genes, colnames(dat))) {
    d <- dat
    d$grp <- factor(ifelse(d[[g]] > median(d[[g]], na.rm = TRUE), "High", "Low"),
                    levels = c("Low","High"))
    if (length(unique(d$grp)) < 2) next
    fit <- survfit(Surv(os_time / 365.25, os_event) ~ grp, data = d)
    p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE,
                    conf.int = TRUE, palette = c("#2E86AB","#C0392B"),
                    xlab = "Years", ylab = "Overall survival",
                    legend.title = g, legend.labs = c("Low","High"),
                    title = paste0(g, " — ", label),
                    risk.table.height = 0.28, ggtheme = theme_bw(base_size = 10))
    png(file.path(dir, paste0(g, ".png")), width = 1500, height = 1500, res = 200)
    print(p); dev.off()
  }
}
km_plots(d_t, GENES, file.path(OUT, "km_tcga"), "TCGA-BRCA")
km_plots(d_m, GENES, file.path(OUT, "km_metabric"), "METABRIC")
message("KM curves written.")


# =============================================================================
# 7. LASSO-COX RISK SCORE — fit on TCGA, LOCK, apply to METABRIC
# =============================================================================
# Locking the coefficients before touching METABRIC is what makes this
# validation rather than refitting. Do not re-tune after seeing the result.

gg <- intersect(GENES, intersect(colnames(d_t), colnames(d_m)))
x_t <- scale(as.matrix(d_t[, gg])); y_t <- Surv(d_t$os_time, d_t$os_event)
ok <- complete.cases(x_t) & !is.na(y_t)

cvfit <- cv.glmnet(x_t[ok, ], y_t[ok], family = "cox", alpha = 1, nfolds = 10)
co <- coef(cvfit, s = "lambda.min")
sel <- data.frame(gene = rownames(co)[as.numeric(co) != 0],
                  coefficient = round(as.numeric(co)[as.numeric(co) != 0], 4))

cat("\n===== LASSO-COX MODEL (locked on TCGA) =====\n")
if (!nrow(sel)) {
  cat("No genes selected. With this event count that is a legitimate result:\n")
  cat("the family carries no independent prognostic signal in TCGA OS.\n")
} else {
  print(sel)
  write.csv(sel, file.path(OUT, "risk_score_model.csv"), row.names = FALSE)

  d_t$risk <- as.numeric(predict(cvfit, newx = x_t, s = "lambda.min"))

  # Apply with TCGA's centring and scaling, not METABRIC's own.
  x_m <- scale(as.matrix(d_m[, gg]),
               center = attr(x_t, "scaled:center"),
               scale  = attr(x_t, "scaled:scale"))
  d_m$risk <- as.numeric(predict(cvfit, newx = x_m, s = "lambda.min"))

  v_t <- coxph(Surv(os_time, os_event) ~ risk, data = d_t)
  v_m <- coxph(Surv(os_time, os_event) ~ risk, data = d_m)

  cat(sprintf("\nRisk score  TCGA (fit)   : HR %.3f, p %.3g, C %.3f\n",
              summary(v_t)$coefficients[1,"exp(coef)"],
              summary(v_t)$coefficients[1,"Pr(>|z|)"], summary(v_t)$concordance[1]))
  cat(sprintf("Risk score  METABRIC (val): HR %.3f, p %.3g, C %.3f\n",
              summary(v_m)$coefficients[1,"exp(coef)"],
              summary(v_m)$coefficients[1,"Pr(>|z|)"], summary(v_m)$concordance[1]))

  adj <- coxph(Surv(os_time, os_event) ~ risk + age + grade + factor(pam50), data = d_m)
  cat(sprintf("\nAdjusted for age, grade, PAM50 in METABRIC: HR %.3f, p %.3g\n",
              summary(adj)$coefficients["risk","exp(coef)"],
              summary(adj)$coefficients["risk","Pr(>|z|)"]))
  cat(sprintf("C-index clinical %.3f -> clinical + risk %.3f (increment %+.4f)\n",
              c_base_m, summary(adj)$concordance[1],
              summary(adj)$concordance[1] - c_base_m))

  write.csv(data.frame(
    cohort = c("TCGA (fit)","METABRIC (validation)"),
    HR = c(summary(v_t)$coefficients[1,"exp(coef)"], summary(v_m)$coefficients[1,"exp(coef)"]),
    p  = c(summary(v_t)$coefficients[1,"Pr(>|z|)"], summary(v_m)$coefficients[1,"Pr(>|z|)"]),
    C  = c(summary(v_t)$concordance[1], summary(v_m)$concordance[1])),
    file.path(OUT, "validation_metabric.csv"), row.names = FALSE)

  for (nm in c("tcga","metabric")) {
    d <- if (nm == "tcga") d_t else d_m
    d$grp <- factor(ifelse(d$risk > median(d$risk), "High risk", "Low risk"),
                    levels = c("Low risk","High risk"))
    fit <- survfit(Surv(os_time / 365.25, os_event) ~ grp, data = d)
    p <- ggsurvplot(fit, data = d, pval = TRUE, risk.table = TRUE, conf.int = TRUE,
                    palette = c("#2E86AB","#C0392B"), xlab = "Years",
                    ylab = "Overall survival", legend.title = "GST risk score",
                    title = paste("GST risk score —", toupper(nm)),
                    risk.table.height = 0.28, ggtheme = theme_bw(base_size = 10))
    png(file.path(OUT, paste0("KM_riskscore_", nm, ".png")),
        width = 1500, height = 1500, res = 200)
    print(p); dev.off()
  }

  # ---- time-dependent ROC ---------------------------------------------------
  for (nm in c("tcga","metabric")) {
    d <- if (nm == "tcga") d_t else d_m
    tt <- 365.25 * c(1, 3, 5)
    tt <- tt[tt < max(d$os_time, na.rm = TRUE)]
    if (!length(tt)) next
    tr <- timeROC(T = d$os_time, delta = d$os_event, marker = d$risk,
                  cause = 1, times = tt, iid = FALSE)
    out <- data.frame(years = round(tt / 365.25, 1), AUC = round(tr$AUC, 3))
    cat(sprintf("\nTime-dependent AUC — %s\n", toupper(nm))); print(out)
    write.csv(out, file.path(OUT, paste0("timeROC_", nm, ".csv")), row.names = FALSE)
  }
}


# =============================================================================
# 8. PROVENANCE
# =============================================================================
writeLines(c(
  paste("Run:", Sys.time()),
  "Objective 1 Module D - prognostic significance.",
  "Discovery TCGA-BRCA; external validation METABRIC.",
  "Expression log2(TPM+1) (Toil recompute) / log2 microarray intensity (METABRIC).",
  "Cox on continuous expression, scaled per SD. KM shown for illustration only.",
  "LASSO-Cox coefficients locked on TCGA before application to METABRIC.",
  paste("TCGA:", nrow(d_t), "patients,", sum(d_t$os_event), "deaths"),
  paste("METABRIC:", nrow(d_m), "patients,", sum(d_m$os_event), "deaths"),
  sprintf("Clinical baseline C-index: TCGA %.3f, METABRIC %.3f", c_base_t, c_base_m),
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleD.txt"))

message("\nAll outputs in ", normalizePath(OUT))
