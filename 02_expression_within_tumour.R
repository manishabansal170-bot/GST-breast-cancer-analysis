# =============================================================================
#  GST FAMILY -- WITHIN-TUMOUR ANALYSIS ACROSS PAM50 SUBTYPES
#  All genes, log2(TPM + 1), no normal groups
#
#  SOURCE = "gdc"   TCGAbiolinks / GDC STAR-Counts  (needs full download)
#  SOURCE = "toil"  UCSC Toil recompute             (works now)
#
#  Same TCGA tumours either way. GDC quantifies with STAR + GENCODE v36,
#  Toil with RSEM + an earlier GENCODE. Run both and compare - agreement
#  across two pipelines is a result worth reporting.
#
#  OUTPUTS  (all in figures/withintumour_<source>/)
#    subtype_summary_<source>.csv     median / IQR / n per gene per subtype
#    subtype_medians_wide_<source>.csv  compact gene x subtype table
#    persample_matrix_<source>.csv    every sample, for the supplement
#    subtype_stats_<source>.csv       Kruskal-Wallis, spread, epsilon-squared
#    pairwise_<source>.csv            all subtype pairs, BH-adjusted
#    Panel_withintumour_<source>.png  all genes, one grid
# =============================================================================

SOURCE  <- "toil"
WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures", paste0("withintumour_", SOURCE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

for (p in c("dplyr","tidyr","tibble","ggplot2"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")

SUB <- c("Basal","Her2","LumA","LumB")


# ---------------------------------------------------------------------------
#  DATA
# ---------------------------------------------------------------------------
if (SOURCE == "gdc") {

  suppressPackageStartupMessages({
    library(TCGAbiolinks); library(SummarizedExperiment) })

  f <- file.path(CACHE, "brca_se.rds")
  if (!file.exists(f))
    stop("cache/brca_se.rds not found - the GDC download has not finished.\n",
         "  Use SOURCE <- \"toil\" for now; it is the same TCGA tumours.")
  se <- readRDS(f)

  mat <- log2(assay(se, "tpm_unstrand") + 1)
  rownames(mat) <- sub("\\..*$", "", rownames(mat))

  gm <- data.frame(ensembl = sub("\\..*$", "", rowData(se)$gene_id),
                   symbol = rowData(se)$gene_name,
                   type = rowData(se)$gene_type, stringsAsFactors = FALSE)
  pick <- function(s) { c1 <- gm[gm$symbol == s, ]
    if (!nrow(c1)) return(NA_character_)
    if (any(c1$type == "protein_coding")) c1 <- c1[c1$type == "protein_coding", ]
    c1$ensembl[1] }
  hit <- vapply(GST, pick, character(1))
  biotype <- setNames(gm$type[match(hit, gm$ensembl)], names(hit))
  mat <- mat[hit[!is.na(hit)], , drop = FALSE]
  rownames(mat) <- names(hit)[!is.na(hit)]

  clin <- as.data.frame(colData(se))
  clin$sample <- colnames(se); clin$patient <- substr(clin$sample, 1, 12)
  pc <- grep("PAM50", colnames(clin), value = TRUE)[1]
  clin$pam50 <- if (!is.na(pc)) as.character(clin[[pc]]) else {
    st <- TCGAquery_subtype(tumor = "brca")
    st$BRCA_Subtype_PAM50[match(clin$patient, st$patient)] }

  # PRIMARY TUMOURS ONLY. Normals excluded before subtype is applied, so no
  # adjacent normal can inherit a PAM50 label.
  clin <- clin[clin$sample_type == "Primary Tumor", ]
  clin$group <- clin$pam50

} else {

  suppressPackageStartupMessages({
    library(UCSCXenaTools); library(TCGAbiolinks) })

  raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
  lin <- 2^raw - 0.001; lin[lin < 0] <- 0
  mat <- log2(lin + 1)
  biotype <- setNames(rep(NA_character_, nrow(mat)), rownames(mat))

  ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
  cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
  ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                         study = .data[[cl("study")]],
                         tissue = .data[[cl("primary disease or tissue")]],
                         stype = .data[[cl("sample_type")]])
  st <- TCGAquery_subtype(tumor = "brca")

  clin <- ph %>%
    filter(toupper(study) == "TCGA", grepl("breast", tissue, ignore.case = TRUE),
           grepl("Primary Tumor", stype)) %>%
    mutate(patient = substr(sample, 1, 12),
           group   = st$BRCA_Subtype_PAM50[match(patient, st$patient)])
}

clin <- clin[clin$group %in% SUB & !is.na(clin$group), ]
clin <- clin %>% arrange(sample) %>% distinct(patient, .keep_all = TRUE)
clin$group <- factor(clin$group, levels = SUB)
clin <- clin[clin$sample %in% colnames(mat), ]
mat  <- mat[, clin$sample, drop = FALSE]

cat("\n===== TUMOURS PER SUBTYPE =====\n"); print(table(clin$group)); cat("\n")

df <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(clin[, c("sample","group")], by = "sample") %>%
  pivot_longer(-c(sample, group), names_to = "gene", values_to = "expr") %>%
  filter(!is.na(group), is.finite(expr))


# ---------------------------------------------------------------------------
#  TABLE 1 -- long summary, median / IQR / n per gene per subtype
# ---------------------------------------------------------------------------
summ <- df %>%
  group_by(gene, group) %>%
  summarise(n          = n(),
            median_log2 = round(median(expr), 3),
            TPM         = round(2^median(expr) - 1, 3),
            Q1          = round(quantile(expr, .25), 3),
            Q3          = round(quantile(expr, .75), 3),
            pct_zero    = round(100 * mean(expr == 0), 1),
            .groups = "drop")

write.csv(summ, file.path(OUT, paste0("subtype_summary_", SOURCE, ".csv")),
          row.names = FALSE)


# ---------------------------------------------------------------------------
#  TABLE 2 -- compact wide table, the one for your thesis
# ---------------------------------------------------------------------------
wide <- summ %>%
  select(gene, group, TPM) %>%
  pivot_wider(names_from = group, values_from = TPM) %>%
  mutate(max_TPM = apply(across(all_of(SUB)), 1, max, na.rm = TRUE),
         min_TPM = apply(across(all_of(SUB)), 1, min, na.rm = TRUE),
         fold_hi_lo = round(ifelse(min_TPM > 0, max_TPM / min_TPM, NA), 2),
         highest = SUB[apply(across(all_of(SUB)), 1, which.max)],
         lowest  = SUB[apply(across(all_of(SUB)), 1, which.min)],
         # detection judged WITHIN tumours, since a gene can be silenced in
         # tumour yet expressed in normal - pooling would mislabel it
         status = ifelse(max_TPM >= 1, "detected", "not detected")) %>%
  arrange(desc(max_TPM))

cat("===== MEDIAN TPM BY SUBTYPE =====\n")
print(as.data.frame(wide), digits = 3)
write.csv(wide, file.path(OUT, paste0("subtype_medians_wide_", SOURCE, ".csv")),
          row.names = FALSE)


# ---------------------------------------------------------------------------
#  TABLE 3 -- per-sample matrix, for the supplement
# ---------------------------------------------------------------------------
as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(clin[, c("sample","group")], by = "sample") %>%
  relocate(group, .after = sample) %>%
  write.csv(file.path(OUT, paste0("persample_matrix_", SOURCE, ".csv")),
            row.names = FALSE)


# ---------------------------------------------------------------------------
#  TABLE 4 -- across-subtype statistics
# ---------------------------------------------------------------------------
# epsilon-squared is the effect size for Kruskal-Wallis: the proportion of
# rank variance explained by subtype. H / ((n^2 - 1)/(n + 1)).
# Rough reading: <0.01 negligible, <0.08 small, <0.26 medium, else large.
# With n over 1,000 the p-value is nearly always tiny; this is the number
# that tells you whether subtype actually structures the gene.

stats <- df %>%
  group_by(gene) %>%
  summarise(n = n(),
            H = unname(kruskal.test(expr ~ droplevels(group))$statistic),
            p = kruskal.test(expr ~ droplevels(group))$p.value,
            spread_log2 = round(diff(range(tapply(expr, droplevels(group), median))), 3),
            .groups = "drop") %>%
  mutate(eps2 = round(H / ((n^2 - 1) / (n + 1)), 4),
         FDR  = p.adjust(p, "BH"),
         effect = cut(eps2, c(-Inf, .01, .08, .26, Inf),
                      labels = c("negligible","small","medium","large"))) %>%
  left_join(wide[, c("gene","max_TPM","status","highest","lowest")], by = "gene") %>%
  arrange(desc(eps2))

cat("\n===== SUBTYPE STRUCTURE, RANKED BY EFFECT SIZE =====\n")
print(as.data.frame(stats), digits = 3)
write.csv(stats, file.path(OUT, paste0("subtype_stats_", SOURCE, ".csv")),
          row.names = FALSE)


# ---------------------------------------------------------------------------
#  TABLE 5 -- all pairwise comparisons
# ---------------------------------------------------------------------------
pw <- lapply(unique(df$gene), function(g) {
  d <- df[df$gene == g, ]
  m <- pairwise.wilcox.test(d$expr, droplevels(d$group),
                            p.adjust.method = "BH")$p.value
  data.frame(gene = g, as.data.frame(as.table(m)))
}) %>% bind_rows() %>% filter(!is.na(Freq)) %>%
  rename(subtype1 = Var1, subtype2 = Var2, p_adj = Freq) %>%
  mutate(p_adj = signif(p_adj, 3),
         sig = ifelse(p_adj < .0001, "****",
               ifelse(p_adj < .001, "***",
               ifelse(p_adj < .01, "**",
               ifelse(p_adj < .05, "*", "ns")))))

write.csv(pw, file.path(OUT, paste0("pairwise_", SOURCE, ".csv")), row.names = FALSE)


# ---------------------------------------------------------------------------
#  PANEL FIGURE -- genes ordered by how strongly subtype structures them
# ---------------------------------------------------------------------------
df$gene <- factor(df$gene, levels = stats$gene)

n_tab <- df %>% distinct(sample, group) %>% count(group) %>%
  mutate(lab = paste0(group, "\n(n=", n, ")"))
df$glab <- factor(n_tab$lab[match(df$group, n_tab$group)], levels = n_tab$lab)
PAL <- setNames(c("#C0392B","#8E7CC3","#2E86AB","#E8A33D")[
  match(n_tab$group, SUB)], n_tab$lab)

p <- ggplot(df, aes(glab, expr, fill = glab)) +
  geom_violin(scale = "width", width = .85, alpha = .3, colour = NA, linewidth = .2) +
  geom_boxplot(width = .28, outlier.size = .2, outlier.alpha = .25, linewidth = .28) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = expression(log[2](TPM + 1)),
       title = "GST family across TCGA-BRCA PAM50 subtypes (primary tumours only)",
       subtitle = paste0(
         if (SOURCE == "toil") "UCSC Toil recompute" else "GDC STAR-Counts",
         "; panels ordered by epsilon-squared, strongest subtype structure first")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.3),
        strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, paste0("Panel_withintumour_", SOURCE, ".png")), p,
       width = 10.5, height = 13, dpi = 300)
ggsave(file.path(OUT, paste0("Panel_withintumour_", SOURCE, ".pdf")), p,
       width = 10.5, height = 13)

message("\nAll tables and figures written to ", normalizePath(OUT))
message("For individual gene figures with brackets, run ",
        "GST_pergene_publication.R with INCLUDE_NORMAL <- FALSE")
