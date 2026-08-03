# =============================================================================
#  METABRIC VALIDATION -- via the cBioPortal API
#
#  Replaces the tarball download, which is returning HTTP 403. The API pulls
#  only the 22 GST genes rather than the whole 150 MB study, so it is faster
#  and avoids the broken S3 endpoint entirely.
#
#  SAME CAVEAT AS BEFORE
#  METABRIC is Illumina HT-12 microarray, log2 intensity. There is no TPM.
#  What replicates across platforms is RANK STRUCTURE, not absolute level.
#  Microarrays compress at both ends, so a gene near background in RNA-seq
#  (GSTA1 in luminal tumours) will show a muted effect here. That is a
#  platform ceiling, not a failed replication.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures", "metabric")
for (d in c(CACHE, OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("cBioPortalData", quietly = TRUE))
  BiocManager::install("cBioPortalData", ask = FALSE, update = FALSE)
for (p in c("dplyr","tidyr","tibble","ggplot2"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(cBioPortalData); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2)
})

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB <- c("Basal","Her2","LumA","LumB")


# ---------------------------------------------------------------------------
# 1. CONNECT AND FIND THE RIGHT PROFILE
# ---------------------------------------------------------------------------
cache_f <- file.path(CACHE, "metabric_gst.rds")

if (file.exists(cache_f)) {
  obj <- readRDS(cache_f); mat <- obj$mat; clin <- obj$clin
  message("Loaded METABRIC from cache.")
} else {

  cbio <- cBioPortal()

  mp <- molecularProfiles(cbio, "brca_metabric")
  cat("\n===== MOLECULAR PROFILES AVAILABLE =====\n")
  print(as.data.frame(mp[, c("molecularProfileId","molecularAlterationType")]))

  # Prefer the raw intensity profile over any z-scored one. A z-score profile
  # is already centred within the cohort, which destroys the between-gene
  # hierarchy you want to compare against TCGA.
  ids <- mp$molecularProfileId
  PROF <- ids[grepl("mrna", ids) & !grepl("zscore|Zscore|median_Zscores", ids)][1]
  if (is.na(PROF)) PROF <- ids[grepl("mrna", ids)][1]
  message("\nUsing profile: ", PROF)

  dat <- getDataByGenes(cbio, studyId = "brca_metabric",
                        genes = GST, by = "hugoGeneSymbol",
                        molecularProfileIds = PROF)[[1]]

  message("Returned ", nrow(dat), " measurements")

  mat_df <- dat %>%
    select(hugoGeneSymbol, sampleId, value) %>%
    filter(!is.na(value)) %>%
    group_by(hugoGeneSymbol, sampleId) %>%
    summarise(value = mean(value), .groups = "drop") %>%   # collapse duplicate probes
    pivot_wider(names_from = sampleId, values_from = value)

  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$hugoGeneSymbol

  # -------------------------------------------------------------------------
  clin_raw <- clinicalData(cbio, "brca_metabric")

  sub_col <- grep("CLAUDIN_SUBTYPE|PAM50", colnames(clin_raw),
                  value = TRUE, ignore.case = TRUE)[1]
  message("Subtype column: ", sub_col)

  clin <- data.frame(sample  = clin_raw$sampleId,
                     subtype = as.character(clin_raw[[sub_col]]),
                     stringsAsFactors = FALSE)

  saveRDS(list(mat = mat, clin = clin), cache_f)
  message("Cached to ", cache_f)
}


# ---------------------------------------------------------------------------
# 2. COVERAGE AND SUBTYPES
# ---------------------------------------------------------------------------
cat("\n===== GENE COVERAGE ON THE ARRAY =====\n")
cat("Present:", paste(sort(rownames(mat)), collapse = ", "), "\n")
cat("Absent :", paste(setdiff(GST, rownames(mat)), collapse = ", "), "\n")
cat("\nAbsence means no probe on the HT-12 array. Platform limitation,\n")
cat("not evidence about expression. Report as not assayed.\n")

cat("\n===== ALL SUBTYPE LABELS =====\n")
print(table(clin$subtype, useNA = "ifany"))

# *** TRAP ***
# "Normal" here is PAM50 NORMAL-LIKE TUMOUR, not normal breast tissue.
# Using it as a normal control would be a serious error. Excluded below.
# "claudin-low" is a subtype PAM50 does not assign; excluded from the TCGA
# comparison because there is nothing to compare it against.

clin <- clin[clin$subtype %in% SUB, ]
clin <- clin[clin$sample %in% colnames(mat), ]
clin$subtype <- factor(clin$subtype, levels = SUB)
mat  <- mat[, clin$sample, drop = FALSE]

cat("\n===== ANALYSIS SET =====\n"); print(table(clin$subtype)); cat("\n")

df <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(clin, by = "sample") %>%
  pivot_longer(-c(sample, subtype), names_to = "gene", values_to = "expr") %>%
  filter(!is.na(subtype), is.finite(expr))


# ---------------------------------------------------------------------------
# 3. TABLES
# ---------------------------------------------------------------------------
mb_wide <- df %>% group_by(gene, subtype) %>%
  summarise(m = round(median(expr), 3), .groups = "drop") %>%
  pivot_wider(names_from = subtype, values_from = m) %>%
  mutate(overall = round(apply(across(all_of(SUB)), 1, median), 3),
         highest = SUB[apply(across(all_of(SUB)), 1, which.max)],
         lowest  = SUB[apply(across(all_of(SUB)), 1, which.min)],
         spread  = round(apply(across(all_of(SUB)), 1, max) -
                           apply(across(all_of(SUB)), 1, min), 3)) %>%
  arrange(desc(overall))

cat("===== METABRIC MEDIAN INTENSITY BY SUBTYPE =====\n")
print(as.data.frame(mb_wide), digits = 3)
write.csv(mb_wide, file.path(OUT, "metabric_subtype_medians.csv"), row.names = FALSE)

stats <- df %>% group_by(gene) %>%
  summarise(n = n(),
            H = unname(kruskal.test(expr ~ droplevels(subtype))$statistic),
            p = kruskal.test(expr ~ droplevels(subtype))$p.value, .groups = "drop") %>%
  mutate(eps2 = round(H / ((n^2 - 1) / (n + 1)), 4),
         FDR = p.adjust(p, "BH")) %>%
  arrange(desc(eps2))

cat("\n===== SUBTYPE STRUCTURE IN METABRIC =====\n")
print(as.data.frame(stats), digits = 3)
write.csv(stats, file.path(OUT, "metabric_subtype_stats.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# 4. CONCORDANCE WITH TCGA -- the replication test
# ---------------------------------------------------------------------------
tcga_f <- file.path(WORKDIR, "figures", "withintumour_toil",
                    "subtype_medians_wide_toil.csv")

if (file.exists(tcga_f)) {

  tc <- read.csv(tcga_f, stringsAsFactors = FALSE)
  common <- intersect(tc$gene, mb_wide$gene)
  message("\nMeasurable in both cohorts: ", length(common), " genes")

  a <- tc$max_TPM[match(common, tc$gene)]
  b <- mb_wide$overall[match(common, mb_wide$gene)]
  ct <- suppressWarnings(cor.test(a, b, method = "spearman"))

  cat("\n===== TEST 1: FAMILY EXPRESSION HIERARCHY =====\n")
  cat(sprintf("Spearman rho = %.3f, p = %.3g, n = %d genes\n",
              unname(ct$estimate), ct$p.value, length(common)))
  cat("rho above ~0.7 means the hierarchy replicates across platforms.\n")

  conc <- lapply(common, function(g) {
    x <- unlist(tc[tc$gene == g, SUB]); y <- unlist(mb_wide[mb_wide$gene == g, SUB])
    data.frame(gene = g,
               TCGA_high = SUB[which.max(x)], MB_high = SUB[which.max(y)],
               TCGA_low  = SUB[which.min(x)], MB_low  = SUB[which.min(y)],
               rho_4pt = round(cor(rank(x), rank(y), method = "spearman"), 3))
  }) %>% bind_rows() %>%
    mutate(verdict = case_when(
      TCGA_high == MB_high & TCGA_low == MB_low ~ "replicated",
      TCGA_high == MB_high | TCGA_low == MB_low ~ "partial",
      TRUE ~ "discordant")) %>%
    arrange(desc(rho_4pt))

  cat("\n===== TEST 2: SUBTYPE ORDERING PER GENE =====\n")
  print(as.data.frame(conc))
  write.csv(conc, file.path(OUT, "tcga_metabric_concordance.csv"), row.names = FALSE)

  cat(sprintf("\n%d replicated, %d partial, %d discordant of %d\n",
              sum(conc$verdict == "replicated"), sum(conc$verdict == "partial"),
              sum(conc$verdict == "discordant"), nrow(conc)))

  g <- conc[conc$gene == "GSTA1", ]
  if (nrow(g))
    cat(sprintf("\nGSTA1 -> TCGA high %s / low %s | METABRIC high %s / low %s : %s\n",
                g$TCGA_high, g$TCGA_low, g$MB_high, g$MB_low, g$verdict))

} else {
  message("\nRun GST_within_tumour.R first to enable the concordance test.")
}


# ---------------------------------------------------------------------------
# 5. FIGURE
# ---------------------------------------------------------------------------
df$gene <- factor(df$gene, levels = stats$gene)
n_tab <- df %>% distinct(sample, subtype) %>% count(subtype) %>%
  mutate(lab = paste0(subtype, "\n(n=", n, ")"))
df$glab <- factor(n_tab$lab[match(df$subtype, n_tab$subtype)], levels = n_tab$lab)
PAL <- setNames(c("#C0392B","#8E7CC3","#2E86AB","#E8A33D")[
  match(n_tab$subtype, SUB)], n_tab$lab)

p <- ggplot(df, aes(glab, expr, fill = glab)) +
  geom_violin(scale = "width", width = .85, alpha = .3, colour = NA, linewidth = .2) +
  geom_boxplot(width = .28, outlier.size = .2, outlier.alpha = .25, linewidth = .28) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = "log2 microarray intensity (Illumina HT-12)",
       title = "GST family across PAM50 subtypes in METABRIC",
       subtitle = "Independent cohort; microarray intensity, not comparable in absolute terms to RNA-seq TPM") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.3),
        strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "Panel_METABRIC.png"), p, width = 10.5, height = 13, dpi = 300)
ggsave(file.path(OUT, "Panel_METABRIC.pdf"), p, width = 10.5, height = 13)

write.csv(as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
            left_join(clin, by = "sample"),
          file.path(OUT, "metabric_persample.csv"), row.names = FALSE)

message("\nWritten to ", normalizePath(OUT))


# =============================================================================
#  IF THE API ALSO FAILS
# =============================================================================
#  Manual route: go to https://www.cbioportal.org/datasets, find
#  "Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)", download the
#  study archive, extract it into cache/brca_metabric/, then run the original
#  METABRIC_validation.R which reads from that folder.
#
#  Note the expression file was renamed in a recent datahub release: older
#  copies use data_expression_median.txt, newer ones
#  data_mrna_illumina_microarray.txt. The original script matches both.
# =============================================================================
