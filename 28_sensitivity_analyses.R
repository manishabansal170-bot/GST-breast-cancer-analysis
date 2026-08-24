# =============================================================================
#  28_sensitivity_analyses.R
#
#  Three checks requested at review, each testing whether a choice made during
#  the analysis determines the result reported.
#
#  C1  THE UNIVERSAL LOSS AXIS IS DEFINED POST HOC
#      Membership was defined by direction in every subtype rather than by
#      pooled effect size, and the manuscript states openly that this
#      definition was chosen because GSTP1 has a larger pooled delta than
#      GSTA4 yet does not fit the axis. A definition selected after seeing
#      which genes it admits needs a sensitivity analysis under the
#      alternative. This script reports both memberships side by side.
#
#  C3  THE PAIRED STRUCTURE IS DISCARDED
#      A substantial number of TCGA-BRCA adjacent normals come from patients
#      who also contributed a tumour. Unpaired Wilcoxon throws that away and
#      is less powerful. The paired test is reported here for the loss-axis
#      genes, along with how many pairs actually exist.
#
#  C5  THE 1 TPM DETECTION FLOOR IS UNJUSTIFIED
#      That threshold determines Table 1, the "16 of 22" headline, and which
#      genes enter every downstream analysis. Sensitivity is reported at 0.5,
#      1 and 2 TPM.
#
#  None of these is expected to overturn the paper. The purpose is to show
#  which conclusions depend on a threshold and which do not, so that a reader
#  does not have to take the choices on trust.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")
LOSS_PUBLISHED <- c("GSTM5","GSTM2","GSTA1","MGST1","GSTA4")


# =============================================================================
# SETUP
# =============================================================================
get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
LINEAR <- 2^get_matrix(raw) - 0.001; LINEAR[LINEAR < 0] <- 0   # TPM
E <- log2(LINEAR + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]], study = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

gtex <- ph %>% dplyr::filter(toupper(study) == "GTEX",
                             grepl("breast", tissue, ignore.case = TRUE)) %>%
  mutate(group = "GTEx Normal")
adjn <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                             grepl("breast", tissue, ignore.case = TRUE),
                             grepl("Solid Tissue Normal", stype)) %>%
  mutate(group = "AdjNormal", patient = substr(sample, 1, 12))
tum <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                            grepl("breast", tissue, ignore.case = TRUE),
                            grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

meta <- bind_rows(gtex %>% dplyr::select(sample, group),
                  adjn %>% dplyr::select(sample, group),
                  tum  %>% dplyr::select(sample, group))
meta <- meta[meta$sample %in% colnames(E), ]
genes <- rownames(E)[apply(E[, meta$sample, drop = FALSE], 1,
                           function(v) sd(v, na.rm = TRUE) > 1e-8)]
cat("Samples:", nrow(meta), " Genes with signal:", length(genes), "\n\n")


# =============================================================================
# C5. DETECTION FLOOR SENSITIVITY
# =============================================================================
cat(strrep("=", 70), "\nC5. DETECTION THRESHOLD SENSITIVITY\n", strrep("=", 70), "\n", sep = "")

grp_med_tpm <- sapply(c("GTEx Normal","AdjNormal", SUB), function(g) {
  s <- meta$sample[meta$group == g]
  apply(LINEAR[genes, s, drop = FALSE], 1, median, na.rm = TRUE) })
max_grp <- apply(grp_med_tpm, 1, max, na.rm = TRUE)

thr_tab <- lapply(c(0.5, 1, 2), function(th) {
  det <- names(max_grp)[max_grp >= th]
  data.frame(threshold_TPM = th, n_detected = length(det),
             genes = paste(sort(det), collapse = ", ")) }) %>% bind_rows()
print(thr_tab[, c("threshold_TPM","n_detected")], row.names = FALSE)

d1 <- names(max_grp)[max_grp >= 1]
for (th in c(0.5, 2)) {
  dt <- names(max_grp)[max_grp >= th]
  cat(sprintf("\nAt %.1f TPM against 1 TPM: %+d genes\n", th, length(dt) - length(d1)))
  add <- setdiff(dt, d1); drop <- setdiff(d1, dt)
  if (length(add))  cat("  gained:", paste(sort(add), collapse = ", "), "\n")
  if (length(drop)) cat("  lost  :", paste(sort(drop), collapse = ", "), "\n")
  if (!length(add) && !length(drop)) cat("  membership unchanged\n")
}
cat("\nGenes near the boundary, maximum group median 0.3 to 3 TPM:\n")
near <- sort(max_grp[max_grp >= 0.3 & max_grp <= 3])
print(round(near, 2))
write.csv(thr_tab, file.path(OUT, "sensitivity_detection_threshold.csv"), row.names = FALSE)


# =============================================================================
# C1. LOSS AXIS UNDER AN ALTERNATIVE DEFINITION
# =============================================================================
cat("\n", strrep("=", 70), "\nC1. LOSS AXIS MEMBERSHIP UNDER TWO DEFINITIONS\n",
    strrep("=", 70), "\n", sep = "")

norm_s <- meta$sample[meta$group == "AdjNormal"]
tum_s  <- meta$sample[meta$group %in% SUB]

cliffs <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  r <- rank(c(a, b)); U <- sum(r[seq_along(a)]) - length(a)*(length(a)+1)/2
  2*U/(length(a)*length(b)) - 1 }

info <- lapply(genes, function(g) {
  x <- as.numeric(E[g, ])
  nv <- x[match(norm_s, colnames(E))]
  per <- sapply(SUB, function(k) {
    s <- meta$sample[meta$group == k]
    median(x[match(s, colnames(E))], na.rm = TRUE) })
  data.frame(gene = g,
             pooled_delta = round(cliffs(x[match(tum_s, colnames(E))], nv), 3),
             down_in_all = all(per < median(nv, na.rm = TRUE)),
             n_subtypes_down = sum(per < median(nv, na.rm = TRUE))) }) %>% bind_rows()

info$rank_by_delta <- rank(info$pooled_delta)
by_direction <- info$gene[info$down_in_all]
by_effect    <- info$gene[order(info$pooled_delta)][1:5]

cat("Published definition, direction in every subtype:\n  ",
    paste(sort(by_direction), collapse = ", "), "\n", sep = "")
cat("Alternative definition, five largest negative pooled deltas:\n  ",
    paste(sort(by_effect), collapse = ", "), "\n", sep = "")
cat("\nIn one but not the other:\n")
cat("  direction only:", paste(setdiff(by_direction, by_effect), collapse = ", "), "\n")
cat("  effect only   :", paste(setdiff(by_effect, by_direction), collapse = ", "), "\n")

cat("\nGenes ordered by pooled delta, with how many subtypes fall below normal:\n")
print(as.data.frame(info[order(info$pooled_delta),
                         c("gene","pooled_delta","n_subtypes_down")])[1:10, ],
      row.names = FALSE)
write.csv(info, file.path(OUT, "sensitivity_loss_axis_definition.csv"), row.names = FALSE)

agree <- length(intersect(by_direction, by_effect))
cat(sprintf("\nThe two definitions agree on %d of 5 genes.\n", agree))
if (agree >= 4) {
  cat("Membership is largely robust to the definition. Report the overlap and\n")
  cat("name the gene or genes on which they differ.\n")
} else {
  cat("Membership depends substantially on the definition. This must be stated,\n")
  cat("and the direction-based definition justified on grounds other than fit.\n")
}


# =============================================================================
# C3. PAIRED ANALYSIS
# =============================================================================
cat("\n", strrep("=", 70), "\nC3. PAIRED TUMOUR-NORMAL ANALYSIS\n", strrep("=", 70), "\n", sep = "")

adjn$patient <- substr(adjn$sample, 1, 12)
tum$patient  <- substr(tum$sample, 1, 12)
shared <- intersect(adjn$patient, tum$patient)
cat("Adjacent normals:", nrow(adjn), "| tumours:", nrow(tum),
    "| patients contributing both:", length(shared), "\n")

if (length(shared) >= 20) {
  pn <- adjn$sample[match(shared, adjn$patient)]
  pt <- tum$sample[match(shared, tum$patient)]
  keep <- pn %in% colnames(E) & pt %in% colnames(E)
  pn <- pn[keep]; pt <- pt[keep]
  cat("Usable pairs:", length(pn), "\n\n")

  paired <- lapply(LOSS_PUBLISHED, function(g) {
    if (!g %in% rownames(E)) return(NULL)
    a <- as.numeric(E[g, pt]); b <- as.numeric(E[g, pn])
    ok <- is.finite(a) & is.finite(b)
    w <- suppressWarnings(wilcox.test(a[ok], b[ok], paired = TRUE))
    u <- suppressWarnings(wilcox.test(a[ok], b[ok], paired = FALSE))
    data.frame(gene = g, n_pairs = sum(ok),
               median_diff = round(median(a[ok] - b[ok]), 3),
               paired_p = w$p.value, unpaired_p = u$p.value) }) %>% bind_rows()
  paired$paired_FDR <- p.adjust(paired$paired_p, "BH")
  print(as.data.frame(paired), row.names = FALSE, digits = 3)
  write.csv(paired, file.path(OUT, "sensitivity_paired_analysis.csv"), row.names = FALSE)

  cat("\nAll five remain significant paired:",
      all(paired$paired_FDR < 0.05, na.rm = TRUE), "\n")
  cat("The paired test uses only patients contributing both sample types and is\n")
  cat("the more appropriate comparison where that structure exists. Report it\n")
  cat("alongside the unpaired result rather than in place of it, since the\n")
  cat("unpaired analysis uses the full cohort.\n")
} else {
  cat("Too few matched pairs for a paired analysis. State this in the Methods\n")
  cat("as the reason the unpaired test was used.\n")
}

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Detection: %s",
                     paste(sprintf("%.1f TPM -> %d genes", thr_tab$threshold_TPM,
                                   thr_tab$n_detected), collapse = "; ")),
             sprintf("Loss axis definitions agree on %d of 5", agree),
             sprintf("Matched tumour-normal patients: %d", length(shared)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_sensitivity.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
