# =============================================================================
#  31_verify_reported_numbers.R
#
#  WHAT THIS IS FOR
#    The manuscript claims that every figure and result can be regenerated from
#    the deposited code. A reader cannot check that claim without running the
#    whole pipeline, and a referee asked for evidence rather than assertion.
#
#    This script checks the numbers the manuscript reports against the files
#    the pipeline produced. It does not rerun the analyses; it verifies that
#    what is written in the paper matches what the code wrote to disk.
#
#  WHY THAT IS THE USEFUL CHECK
#    Rerunning everything proves the code executes. It does not prove the
#    manuscript reports what the code produced. Numbers drift during revision:
#    a figure is regenerated, a filter changes, a value in the text is not
#    updated. That gap is what this catches, and it is the failure mode that
#    actually occurs.
#
#  HOW TO READ THE OUTPUT
#    Every check prints PASS, FAIL or ABSENT. A FAIL means the manuscript and
#    the deposit disagree and one of them is wrong. ABSENT means the output
#    file is missing, which is itself worth knowing before submission.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
setwd(WORKDIR)
OUT <- file.path(WORKDIR, "figures")

pass <- 0; fail <- 0; absent <- 0
results <- list()

check <- function(label, claimed, found, tol = 0.02, unit = "") {
  if (is.null(found) || !length(found) || all(is.na(found))) {
    verdict <- "ABSENT"; absent <<- absent + 1
    cat(sprintf("  %-8s %-52s claimed %s, file missing or empty\n",
                verdict, label, paste(claimed, collapse = "/")))
  } else {
    ok <- if (is.numeric(claimed) && is.numeric(found))
            all(abs(found - claimed) <= tol * pmax(abs(claimed), 1))
          else identical(as.character(found), as.character(claimed))
    verdict <- if (ok) "PASS" else "FAIL"
    if (ok) pass <<- pass + 1 else fail <<- fail + 1
    cat(sprintf("  %-8s %-52s claimed %s  found %s %s\n", verdict, label,
                paste(claimed, collapse = "/"), paste(round(found, 4), collapse = "/"), unit))
  }
  results[[length(results) + 1]] <<-
    data.frame(check = label, claimed = paste(claimed, collapse = "/"),
               found = paste(round(as.numeric(found), 4), collapse = "/"),
               verdict = verdict, stringsAsFactors = FALSE)
}

rd <- function(path) {
  p <- file.path(WORKDIR, path)
  if (!file.exists(p)) return(NULL)
  tryCatch(read.csv(p, check.names = FALSE), error = function(e) NULL) }

cat(strrep("=", 78), "\nVERIFYING REPORTED NUMBERS AGAINST DEPOSITED OUTPUT\n",
    strrep("=", 78), "\n\n", sep = "")


# -----------------------------------------------------------------------------
cat("EXPRESSION AND DETECTION\n")
det <- rd("figures/detection_toil.csv")
if (!is.null(det)) {
  check("Genes reaching detection threshold", 16, sum(det$detected))
  check("Genes not detected", 4, sum(!det$detected))
} else { check("Detection table", 16, NULL) }

thr <- rd("figures/sensitivity_detection_threshold.csv")
if (!is.null(thr)) {
  check("Detected at 0.5 TPM", 17, thr$n_detected[thr$threshold_TPM == 0.5])
  check("Detected at 2 TPM",   16, thr$n_detected[thr$threshold_TPM == 2])
} else { check("Detection sensitivity", 17, NULL) }


# -----------------------------------------------------------------------------
cat("\nLOSS AXIS AND EFFECT SIZES\n")
ax <- rd("figures/sensitivity_loss_axis_definition.csv")
if (!is.null(ax)) {
  g <- function(x) ax$pooled_delta[ax$gene == x]
  check("GSTM5 Cliff's delta", -0.905, g("GSTM5"))
  check("GSTM2 Cliff's delta", -0.720, g("GSTM2"))
  check("GSTA1 Cliff's delta", -0.612, g("GSTA1"))
  check("GSTP1 Cliff's delta", -0.434, g("GSTP1"))
  check("GSTO1 rises in tumour", TRUE, g("GSTO1") > 0)
} else { check("Loss axis table", -0.905, NULL) }

pr <- rd("figures/sensitivity_paired_analysis.csv")
if (!is.null(pr)) {
  check("Matched tumour-normal pairs", 109, unique(pr$n_pairs)[1])
  check("Loss genes significant when paired", 5, sum(pr$paired_FDR < 0.05))
} else { check("Paired analysis", 109, NULL) }


# -----------------------------------------------------------------------------
cat("\nAXIS INDEPENDENCE\n")
co <- rd("figures/axis_coordinates.csv")
if (!is.null(co)) {
  ct <- suppressWarnings(cor(abs(co$loss_delta), co$subtype_eps2, method = "spearman"))
  check("Gene-level axis correlation (rho)", 0.226, ct, tol = 0.05)
  check("Genes on both axes", 20, nrow(co))
} else { check("Axis coordinates", 0.226, NULL) }

sl <- rd("figures/axis_independence_sample_level.csv")
if (!is.null(sl)) {
  r <- sl$r[sl$comparison == "All tumours, disjoint gene sets"]
  check("Sample-level axis correlation (r)", -0.066, r, tol = 0.05)
} else { check("Sample-level independence", -0.066, NULL) }

pca <- rd("figures/pca_per_sample_variance.csv")
if (!is.null(pca)) {
  check("PC1 variance (%)", 18.9, pca$variance_pct[1], tol = 0.03)
  check("PC2 variance (%)", 11.3, pca$variance_pct[2], tol = 0.03)
} else { check("Per-sample PCA", 18.9, NULL) }


# -----------------------------------------------------------------------------
cat("\nMETHYLATION\n")
gm <- rd("objective1/moduleC_methylation/GSTP1_methylation_by_subtype.csv")
if (!is.null(gm)) {
  check("Basal-like unmethylated (%)", 97.0, gm$pct_unmethylated[gm$pam50 == "Basal"])
  check("Luminal B methylated (%)",    58.2, gm$pct_methylated[gm$pam50 == "LumB"])
  check("HER2 methylated (%)",         63.0, gm$pct_methylated[gm$pam50 == "Her2"])
  check("Table 3 subtype total",       737,  sum(gm$n))
} else { check("GSTP1 by subtype", 97.0, NULL) }

pp <- rd("objective1/moduleC_methylation/gstp1_per_probe.csv")
if (!is.null(pp)) {
  check("GSTP1 promoter probes with data", 5, sum(pp$n_obs > 0))
} else { check("GSTP1 per-probe", 5, NULL) }

vp <- rd("objective1/moduleBG_alterations_protein/cna_vs_methylation_variance.csv")
if (!is.null(vp)) {
  check("Methylation-dominant genes", 11, sum(grepl("METHYLATION", vp$driver)))
  check("Genes in joint model",       18, nrow(vp))
} else { check("Variance partition", 11, NULL) }


# -----------------------------------------------------------------------------
cat("\nIMMUNE AND PURITY\n")
po <- rd("objective1/moduleG_immune/purity_negative_control_observed.csv")
if (!is.null(po)) {
  check("GST associations surviving (%)", 1.95, po$pct_survive, tol = 0.03)
  check("GST mean attenuation",           0.054, po$mean_atten, tol = 0.05)
} else { check("Purity control", 1.95, NULL) }

pe <- rd("objective1/moduleG_immune/purity_control_epithelial_null.csv")
if (!is.null(pe)) {
  check("Epithelial-matched control median attenuation", -0.011,
        median(pe$mean_atten), tol = 0.5)
} else { check("Epithelial control", -0.011, NULL) }


# -----------------------------------------------------------------------------
cat("\nSURVIVAL\n")
ic <- rd("objective1/moduleD_survival/incremental_cindex_metabric.csv")
if (!is.null(ic)) {
  check("Maximum incremental C-index", 0.0054, max(ic$increment), tol = 0.05)
  check("Clinical model C-index",      0.621,  unique(ic$C_clinical)[1])
} else { check("Incremental C-index", 0.0054, NULL) }

rs <- rd("objective1/moduleD_survival/riskscore_direction_and_scaling.csv")
if (!is.null(rs)) {
  pubc <- rs$C_index[grepl("As published", rs$configuration)]
  revc <- rs$C_index[grepl("METABRIC-trained", rs$configuration)]
  check("Published-direction C-index", 0.560, pubc)
  check("Reverse-direction C-index",   0.587, revc)
} else { check("Risk score directions", 0.560, NULL) }


# -----------------------------------------------------------------------------
cat("\nNEOADJUVANT RESPONSE\n")
vd <- rd("objective1/moduleE_drug_response/GSE32646_validation_primary.csv")
if (!is.null(vd)) {
  check("GSE32646 GSTM2 odds ratio", 1.170, vd$OR_perSD[vd$gene == "GSTM2"], tol = 0.03)
  check("GSE32646 GSTM1 odds ratio", 1.342, vd$OR_perSD[vd$gene == "GSTM1"], tol = 0.03)
} else { check("GSE32646 validation", 1.170, NULL) }


# -----------------------------------------------------------------------------
cat("\nCO-EXPRESSION\n")
cw <- rd("objective1/moduleF_pathway_network/coexpression_within_subtype.csv")
if (!is.null(cw)) {
  check("Pairs tested", 162, nrow(cw))
  check("Pairs surviving within subtype", 3, sum(cw$survives, na.rm = TRUE))
} else { check("Co-expression within subtype", 162, NULL) }


# -----------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\nSUMMARY\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("  PASS   %d\n  FAIL   %d\n  ABSENT %d\n", pass, fail, absent))

res <- do.call(rbind, results)
write.csv(res, file.path(OUT, "verification_report.csv"), row.names = FALSE)

if (fail > 0) {
  cat("\nDisagreements between the manuscript and the deposit:\n")
  print(res[res$verdict == "FAIL", ], row.names = FALSE)
  cat("\nEach must be resolved before submission. Either the text is stale or the\n")
  cat("output was regenerated with different settings; the file is the record of\n")
  cat("what was computed, so it is usually the text that needs correcting.\n")
} else if (absent > 0) {
  cat("\nAll present numbers agree. Missing files should be regenerated so the\n")
  cat("deposit is complete:\n")
  print(res[res$verdict == "ABSENT", c("check")], row.names = FALSE)
} else {
  cat("\nEvery number checked matches the deposited output. This report can be\n")
  cat("cited in the data availability statement as evidence that the manuscript\n")
  cat("reports what the code produced.\n")
}

writeLines(c(paste("Run:", Sys.time()),
             sprintf("PASS %d, FAIL %d, ABSENT %d", pass, fail, absent),
             "", capture.output(print(res, row.names = FALSE)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_verification.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
