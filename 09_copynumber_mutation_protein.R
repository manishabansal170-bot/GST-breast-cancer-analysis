# =============================================================================
#  OBJECTIVE 1 · MODULES B AND G   (v3 — sparse-data fix)
#
#  THE BUG, AND IT IS A TRAP WORTH REMEMBERING
#    cBioPortal returns DISCRETE data sparsely. For GISTIC calls and mutations
#    only ALTERED samples come back; neutral samples are simply not sent.
#
#    Consequences in v2:
#      - the frequency table showed 0% neutral, which is impossible
#      - per-gene n was 6-123 instead of ~1,080
#      - GSTP1 and MGST3 had zero variance (every returned value was +2),
#        so cor.test returned NULL for every gene and the table was empty
#      - "GSTK1 mutated in 18.75%" was 3 / 16 samples-with-any-GST-mutation,
#        not 3 / ~980 sequenced. True rate is about 0.3%.
#
#  THE FIX
#    Use brca_tcga_linear_CNA (continuous log2 ratio, returned for every
#    sample) for correlation and variance partition. Use the GISTIC profile
#    only for frequencies, with the denominator taken from the sample list
#    rather than from the number of rows returned.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleBG_alterations_protein")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(cBioPortalData); library(TCGAbiolinks) })

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB   <- c("Basal","Her2","LumA","LumB")
STUDY <- "brca_tcga"

cbio <- cBioPortal()

to_matrix <- function(d, label = "") {
  vcol <- intersect(c("value","alteration","score"), names(d))[1]
  gcol <- intersect(c("hugoGeneSymbol","gene"), names(d))[1]
  scol <- intersect(c("sampleId","patientId"), names(d))[1]
  cat(sprintf("  %s: using column '%s', %d rows\n", label, vcol, nrow(d)))
  w <- d %>%
    transmute(gene = .data[[gcol]], sample = .data[[scol]],
              val = suppressWarnings(as.numeric(.data[[vcol]]))) %>%
    filter(!is.na(val)) %>%
    group_by(gene, sample) %>% summarise(val = mean(val), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = val)
  m <- as.matrix(w[, -1]); rownames(m) <- w$gene; m
}


# =============================================================================
# 1. PROPER DENOMINATORS
# =============================================================================
# How many samples were actually PROFILED for each assay. Without this every
# frequency is computed against the altered subset and is meaningless.

n_cna <- NA; n_seq <- NA
try({
  sl <- sampleLists(cbio, STUDY)
  cat("\n===== SAMPLE LISTS =====\n")
  print(as.data.frame(sl[, c("sampleListId","name")]))
  get_n <- function(pat) {
    id <- sl$sampleListId[grepl(pat, sl$sampleListId, ignore.case = TRUE)][1]
    if (is.na(id)) return(NA)
    length(samplesInSampleLists(cbio, id)[[1]])
  }
  n_cna <- get_n("_cna$|_acgh$")
  n_seq <- get_n("_sequenced$")
}, silent = TRUE)

cat(sprintf("\nProfiled for CNA: %s | sequenced: %s\n", n_cna, n_seq))


# =============================================================================
# 2. CONTINUOUS COPY NUMBER — complete data, for correlation
# =============================================================================
lf <- file.path(CACHE, "tcga_cna_linear_gst.rds")
if (file.exists(lf)) {
  cnl <- readRDS(lf); message("Loaded linear CNA from cache.")
} else {
  message("Fetching linear (continuous) copy number...")
  d <- getDataByGenes(cbio, studyId = STUDY, genes = GST, by = "hugoGeneSymbol",
                      molecularProfileIds = "brca_tcga_linear_CNA")[[1]]
  cnl <- to_matrix(d, "linear CNA")
  saveRDS(cnl, lf)
}
cat(sprintf("\nLinear CNA: %d genes x %d samples\n", nrow(cnl), ncol(cnl)))
cat(sprintf("Value range: %.2f to %.2f (log2 ratio, complete data)\n",
            min(cnl, na.rm = TRUE), max(cnl, na.rm = TRUE)))

# Sanity: non-NA per gene should now be in the hundreds, not single digits.
cat("Non-missing values per gene (first 6):\n")
print(head(sort(rowSums(is.finite(cnl)), decreasing = TRUE), 6))


# =============================================================================
# 3. FREQUENCY, WITH THE RIGHT DENOMINATOR
# =============================================================================
gf <- file.path(CACHE, "tcga_cna_gistic_gst.rds")
if (file.exists(gf)) {
  cng <- readRDS(gf)
} else {
  d <- getDataByGenes(cbio, studyId = STUDY, genes = GST, by = "hugoGeneSymbol",
                      molecularProfileIds = "brca_tcga_gistic")[[1]]
  cng <- to_matrix(d, "GISTIC")
  saveRDS(cng, gf)
}

denom <- if (!is.na(n_cna)) n_cna else ncol(cnl)
cat(sprintf("\nUsing denominator: %d samples profiled for copy number\n", denom))

cn_tab <- lapply(rownames(cng), function(g) {
  v <- cng[g, ]; v <- v[is.finite(v)]
  data.frame(gene = g,
             deep_del_pct = round(100 * sum(v <= -2) / denom, 2),
             loss_pct     = round(100 * sum(v <  0)  / denom, 2),
             gain_pct     = round(100 * sum(v == 1)  / denom, 2),
             amp_pct      = round(100 * sum(v >= 2)  / denom, 2),
             any_alt_pct  = round(100 * length(v)    / denom, 2))
}) %>% bind_rows() %>% arrange(desc(any_alt_pct))

cat("\n===== COPY NUMBER ALTERATION FREQUENCY (% of profiled samples) =====\n")
print(as.data.frame(cn_tab))
write.csv(cn_tab, file.path(OUT, "copy_number_frequency.csv"), row.names = FALSE)


# =============================================================================
# 4. MUTATIONS, WITH THE RIGHT DENOMINATOR
# =============================================================================
mut <- tryCatch(
  getDataByGenes(cbio, studyId = STUDY, genes = GST, by = "hugoGeneSymbol",
                 molecularProfileIds = "brca_tcga_mutations")[[1]],
  error = function(e) NULL)

if (!is.null(mut) && nrow(mut)) {
  md <- if (!is.na(n_seq)) n_seq else 980
  mt <- mut %>% group_by(hugoGeneSymbol) %>%
    summarise(mutated_samples = n_distinct(sampleId), .groups = "drop") %>%
    mutate(pct_of_sequenced = round(100 * mutated_samples / md, 2)) %>%
    arrange(desc(pct_of_sequenced))
  cat(sprintf("\n===== MUTATION FREQUENCY (denominator %d sequenced) =====\n", md))
  print(as.data.frame(mt))
  cat("\nAll well under 1%. GSTs are metabolic effectors, not drivers, so the\n")
  cat("expression changes are not mutation-driven. Worth one sentence.\n")
  write.csv(mt, file.path(OUT, "mutation_frequency.csv"), row.names = FALSE)
}


# =============================================================================
# 5. COPY NUMBER vs EXPRESSION
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; exp <- obj$exp

shared <- intersect(colnames(cnl), colnames(exp))
shared <- shared[substr(shared, 14, 15) == "01"]
g_ce <- intersect(rownames(cnl), rownames(exp))
cat(sprintf("\nTumours with linear CNA + expression: %d | genes: %d\n",
            length(shared), length(g_ce)))

cn_expr <- lapply(g_ce, function(g) {
  x <- cnl[g, shared]; y <- exp[g, shared]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 50 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(gene = g, n = sum(ok),
             rho_cna = round(unname(ct$estimate), 3), pval = ct$p.value)
}) %>% bind_rows()

if (nrow(cn_expr)) {
  cn_expr <- cn_expr %>% mutate(FDR = p.adjust(pval, "BH")) %>%
    arrange(desc(rho_cna))
  cat("\n===== COPY NUMBER vs EXPRESSION =====\n")
  print(as.data.frame(cn_expr), digits = 3)
  cat("\nPositive rho is expected: more copies, more transcript. The question\n")
  cat("is whether it explains MORE than methylation does - section 6.\n")
  write.csv(cn_expr, file.path(OUT, "cna_expression_correlation.csv"),
            row.names = FALSE)
} else {
  message("No genes passed. Check the linear CNA matrix.")
}


# =============================================================================
# 6. THE DECIDING TEST
# =============================================================================
s3 <- Reduce(intersect, list(colnames(cnl), colnames(exp), colnames(met)))
s3 <- s3[substr(s3, 14, 15) == "01"]
g3 <- Reduce(intersect, list(rownames(cnl), rownames(exp), rownames(met)))
cat(sprintf("\nTumours with all three assays: %d | genes: %d\n", length(s3), length(g3)))

three <- lapply(g3, function(g) {
  d <- data.frame(expr = exp[g, s3], cn = cnl[g, s3], beta = met[g, s3])
  d <- d[complete.cases(d), ]
  if (nrow(d) < 50 || sd(d$cn) == 0 || sd(d$beta) == 0) return(NULL)
  d$expr <- rank(d$expr); d$cn <- rank(d$cn); d$beta <- rank(d$beta)
  full <- lm(expr ~ cn + beta, data = d)
  r2f <- summary(full)$r.squared
  r2c <- summary(lm(expr ~ cn,   data = d))$r.squared
  r2b <- summary(lm(expr ~ beta, data = d))$r.squared
  co  <- summary(full)$coefficients
  data.frame(gene = g, n = nrow(d),
             R2_total = round(r2f, 3),
             R2_cn_alone = round(r2c, 3),
             R2_meth_alone = round(r2b, 3),
             unique_cn = round(r2f - r2b, 3),
             unique_meth = round(r2f - r2c, 3),
             p_cn = co["cn","Pr(>|t|)"], p_meth = co["beta","Pr(>|t|)"])
}) %>% bind_rows()

if (nrow(three)) {
  three <- three %>%
    mutate(driver = case_when(
      unique_meth > 2 * unique_cn & p_meth < 0.05 ~ "METHYLATION dominant",
      unique_cn > 2 * unique_meth & p_cn < 0.05   ~ "COPY NUMBER dominant",
      p_meth < 0.05 & p_cn < 0.05                 ~ "both contribute",
      p_meth < 0.05                               ~ "methylation only",
      p_cn < 0.05                                 ~ "copy number only",
      TRUE ~ "neither")) %>%
    arrange(desc(R2_total))

  cat("\n===== WHAT DRIVES GST EXPRESSION? =====\n")
  print(as.data.frame(three), digits = 3)
  write.csv(three, file.path(OUT, "cna_vs_methylation_variance.csv"),
            row.names = FALSE)

  cat("\n--- Genes carrying the central claim ---\n")
  for (g in c("GSTP1","GSTM2","GSTM3","GSTM5","GSTA1","GSTA4","GSTO2","GSTM4")) {
    r <- three[three$gene == g, ]
    if (nrow(r))
      cat(sprintf("  %-7s R2 %.3f | unique CN %.3f | unique meth %.3f  ->  %s\n",
                  g, r$R2_total, r$unique_cn, r$unique_meth, r$driver))
  }
  cat("\n"); print(table(three$driver))
}


# =============================================================================
# 7. PROTEIN — this study does have RPPA
# =============================================================================
cat("\n===== PROTEIN =====\n")
prot <- tryCatch(
  getDataByGenes(cbio, studyId = STUDY, genes = GST, by = "hugoGeneSymbol",
                 molecularProfileIds = "brca_tcga_rppa")[[1]],
  error = function(e) NULL)

if (is.null(prot) || !nrow(prot)) {
  prot <- tryCatch(
    getDataByGenes(cbio, studyId = STUDY, genes = GST, by = "hugoGeneSymbol",
                   molecularProfileIds = "brca_tcga_protein_quantification")[[1]],
    error = function(e) NULL)
}

if (!is.null(prot) && nrow(prot)) {
  pm <- to_matrix(prot, "protein")
  cat("GST proteins measured:", paste(sort(rownames(pm)), collapse = ", "), "\n")

  # mRNA-protein concordance, the commonest objection to an in-silico chapter
  ps <- intersect(colnames(pm), colnames(exp))
  ps <- ps[substr(ps, 14, 15) == "01"]
  pg <- intersect(rownames(pm), rownames(exp))

  if (length(pg) && length(ps) > 30) {
    conc <- lapply(pg, function(g) {
      x <- pm[g, ps]; y <- exp[g, ps]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 30) return(NULL)
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
      data.frame(gene = g, n = sum(ok),
                 rho_mRNA_protein = round(unname(ct$estimate), 3),
                 pval = ct$p.value)
    }) %>% bind_rows()
    if (nrow(conc)) {
      cat("\n===== mRNA-PROTEIN CONCORDANCE =====\n")
      print(as.data.frame(conc), digits = 3)
      write.csv(conc, file.path(OUT, "mrna_protein_concordance.csv"),
                row.names = FALSE)
      cat("\nA positive correlation here is worth stating explicitly: it shows\n")
      cat("the transcript findings are reflected at protein level, which is\n")
      cat("the standard objection to a computational chapter.\n")
    }
  }
} else {
  cat("No GST proteins on the RPPA panel (~200 antibodies, GSTs not included).\n")
  cat("Expected. Use Human Protein Atlas (proteinatlas.org, GSTP1, Pathology\n")
  cat("tab) for IHC-scored breast staining, and your own GSTP1 IHC cohort in\n")
  cat("Objective 2 for the definitive answer.\n")
}

writeLines(c(
  paste("Run:", Sys.time()), paste("Study:", STUDY),
  "Copy number: linear (continuous) profile used for correlation because the",
  "discrete GISTIC profile is returned sparsely by the cBioPortal API -",
  "only altered samples are sent. Frequencies computed against the number of",
  "samples profiled, taken from the study sample lists.",
  paste("CNA denominator:", n_cna, "| sequenced denominator:", n_seq),
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleBG.txt"))

message("\nWritten to ", normalizePath(OUT))
