# =============================================================================
#  GST FAMILY -- MASTER FIGURE SCRIPT
#  Regenerates every figure for all 22 genes from corrected data.
#
#  SOURCE = "gdc"   TCGA only: adjacent normal + 4 PAM50 subtypes
#                   Needs the full GDC download (cache/brca_se.rds)
#  SOURCE = "toil"  GTEx normal + TCGA adjacent normal + 4 PAM50 subtypes
#                   Needs no download. Run this one first.
#
#  Both branches feed the SAME plotting code, so your whole figure suite is
#  visually consistent and one edit changes all of it.
#
#  OUTPUTS
#    Fig_panel_<source>.png/.pdf     all genes, one grid
#    Fig_heatmap_<source>.png        group medians, z-scored
#    pergene_<source>/<GENE>.png     one publication figure per gene
#    stats_<source>.csv              medians, log2FC, Cliff's delta, FDR
#    detection_<source>.csv          which genes are usable and why
#    provenance_<source>.txt         your methods paragraph
#
#  ALL CORRECTIONS FROM THE DIAGNOSTIC ARE BAKED IN:
#    - Ensembl version suffix stripped   (fixes GSTM2, GSTM4)
#    - protein_coding preferred on ties
#    - biotype reported, so GSTT2's pseudogene status is visible
#    - GSTT1 handled explicitly: absent from GRCh38 primary assembly
#    - adjacent normals never inherit a PAM50 label
#    - detection floor applied and shown, not silently applied
# =============================================================================

SOURCE  <- "toil"          # "toil" (no download needed) or "gdc"
WORKDIR <- "~/GST_BRCA"

DETECT_TPM <- 1.0          # TPM at or above this = "detected"
SHOW_VIOLIN <- TRUE        # reveals GSTM1 / GSTT1 bimodality

OUT   <- file.path(WORKDIR, "figures")
CACHE <- file.path(WORKDIR, "cache")
PERG  <- file.path(OUT, paste0("pergene_", SOURCE))
for (d in c(OUT, CACHE, PERG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

options(timeout = 1000000)

for (p in c("dplyr","tidyr","tibble","ggplot2","pheatmap","RColorBrewer"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(pheatmap); library(RColorBrewer)
})

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5",
         "GSTM1","GSTM2","GSTM3","GSTM4","GSTM5",
         "GSTP1",
         "GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1",
         "MGST1","MGST2","MGST3")


# =============================================================================
#  DATA
# =============================================================================
if (SOURCE == "gdc") {

  suppressPackageStartupMessages({
    library(TCGAbiolinks); library(SummarizedExperiment) })

  se_file <- file.path(CACHE, "brca_se.rds")
  if (!file.exists(se_file))
    stop("cache/brca_se.rds not found. The GDC download is not finished.\n",
         "  Run GDC_retry_and_diagnose.R, or set SOURCE <- \"toil\".")
  se <- readRDS(se_file)

  mat <- log2(assay(se, "tpm_unstrand") + 1)
  rownames(mat) <- sub("\\..*$", "", rownames(mat))     # THE GSTM2/GSTM4 FIX

  gm <- data.frame(ensembl = sub("\\..*$", "", rowData(se)$gene_id),
                   symbol  = rowData(se)$gene_name,
                   type    = rowData(se)$gene_type,
                   stringsAsFactors = FALSE)

  pick <- function(s) {
    c1 <- gm[gm$symbol == s, ]
    if (!nrow(c1)) return(NA_character_)
    if (any(c1$type == "protein_coding")) c1 <- c1[c1$type == "protein_coding", ]
    c1$ensembl[1]
  }
  hit <- vapply(GST, pick, character(1))
  biotype <- setNames(gm$type[match(hit, gm$ensembl)], names(hit))

  mat <- mat[hit[!is.na(hit)], , drop = FALSE]
  rownames(mat) <- names(hit)[!is.na(hit)]

  clin <- as.data.frame(colData(se))
  clin$sample  <- colnames(se)
  clin$patient <- substr(clin$sample, 1, 12)

  pcol <- grep("PAM50", colnames(clin), value = TRUE)[1]
  clin$pam50 <- if (!is.na(pcol)) as.character(clin[[pcol]]) else {
    st <- TCGAquery_subtype(tumor = "brca")
    st$BRCA_Subtype_PAM50[match(clin$patient, st$patient)] }

  clin <- clin[clin$sample_type %in% c("Primary Tumor","Solid Tissue Normal"), ]
  # sample type resolved BEFORE subtype, so normals never inherit a PAM50 call
  clin$group <- ifelse(clin$sample_type == "Solid Tissue Normal",
                       "Adjacent Normal", clin$pam50)

  LEV <- c("Adjacent Normal","Basal","Her2","LumA","LumB")
  REF <- "Adjacent Normal"

} else {

  if (!requireNamespace("UCSCXenaTools", quietly = TRUE))
    install.packages("UCSCXenaTools")
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("TCGAbiolinks", ask = FALSE, update = FALSE) }
  suppressPackageStartupMessages({
    library(UCSCXenaTools); library(TCGAbiolinks) })

  HOST  <- "https://toil.xenahubs.net"
  EXPR  <- "TcgaTargetGtex_rsem_gene_tpm"
  PHENO <- "TcgaTargetGTEX_phenotype.txt"

  ph_f <- file.path(CACHE, "toil_pheno.rds")
  ph <- if (file.exists(ph_f)) readRDS(ph_f) else {
    x <- XenaGenerate(subset = XenaDatasets == PHENO) %>% XenaQuery() %>%
      XenaDownload(destdir = CACHE, trans_slash = TRUE) %>% XenaPrepare()
    if (is.list(x) && !is.data.frame(x)) x <- x[[1]]
    x <- as.data.frame(x); saveRDS(x, ph_f); x }

  cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
  ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                         study  = .data[[cl("study")]],
                         tissue = .data[[cl("primary disease or tissue")]],
                         stype  = .data[[cl("sample_type")]])

  ex_f <- file.path(CACHE, "toil_gst.rds")
  raw <- if (file.exists(ex_f)) readRDS(ex_f) else {
    x <- fetch_dense_values(host = HOST, dataset = EXPR, identifiers = GST,
                            use_probeMap = TRUE, check = FALSE)
    saveRDS(x, ex_f); x }

  # Toil ships log2(TPM + 0.001). Back-transform, then re-log on +1.
  lin <- 2^as.matrix(raw) - 0.001; lin[lin < 0] <- 0
  mat <- log2(lin + 1)
  biotype <- setNames(rep(NA_character_, nrow(mat)), rownames(mat))

  br <- ph %>% filter(grepl("breast", tissue, ignore.case = TRUE))
  st <- TCGAquery_subtype(tumor = "brca")

  # The study label in the Toil phenotype file is "GTEX" in upper case. Matching
  # case-sensitively on "GTEx" returns nothing and silently drops all 179 GTEx
  # breast samples, producing a plausible-looking figure with a missing group.
  # Matched case-insensitively here, with an assertion below.
  clin <- bind_rows(
    br %>% filter(toupper(study) == "GTEX") %>% transmute(sample, group = "GTEx Normal"),
    br %>% filter(toupper(study) == "TCGA", grepl("Solid Tissue Normal", stype)) %>%
      transmute(sample, group = "TCGA Adjacent Normal"),
    br %>% filter(toupper(study) == "TCGA", grepl("Primary Tumor", stype)) %>%
      mutate(group = st$BRCA_Subtype_PAM50[match(substr(sample,1,12), st$patient)]) %>%
      filter(group %in% c("Basal","Her2","LumA","LumB")) %>% select(sample, group)
  )
  clin$patient <- ifelse(grepl("^TCGA", clin$sample),
                         substr(clin$sample, 1, 12), clin$sample)

  LEV <- c("GTEx Normal","TCGA Adjacent Normal","Basal","Her2","LumA","LumB")
  REF <- "TCGA Adjacent Normal"
}

clin <- clin[clin$group %in% LEV, ]
clin <- clin %>% arrange(sample) %>% distinct(patient, group, .keep_all = TRUE)
clin$group <- factor(clin$group, levels = LEV)
clin <- clin[clin$sample %in% colnames(mat), ]
mat  <- mat[, clin$sample, drop = FALSE]

cat("\n===== SAMPLES PER GROUP =====\n"); print(table(clin$group)); cat("\n")

# An empty group means a filter has failed to match, not that the data are
# missing. Stop rather than produce a figure with a blank panel column.
if (any(table(clin$group) == 0))
  stop("A sample group is empty. Check the study and sample_type filters ",
       "against the actual values in the phenotype file:\n  study: ",
       paste(unique(ph$study), collapse = ", "))
if (any(table(clin$group) < 3)) warning("A group has n < 3.")


# =============================================================================
#  DETECTION CLASSIFICATION -- applied visibly, never silently
# =============================================================================
missing_genes <- setdiff(GST, rownames(mat))

# Detection is judged PER GROUP, not on a pooled median. A gene silenced in
# tumour but expressed in normal tissue is misclassified as undetected when the
# median is dominated by the tumour fraction - GSTA1 (7.12 TPM in adjacent
# normal, 0.13 TPM in luminal B) is the case in point.
grp_med <- sapply(levels(clin$group), function(g) {
  s <- clin$sample[clin$group == g]
  apply(mat[, s, drop = FALSE], 1, median, na.rm = TRUE)
})
max_tpm <- apply(2^grp_med - 1, 1, max, na.rm = TRUE)

det <- data.frame(
  gene          = rownames(mat),
  biotype       = biotype[rownames(mat)],
  median        = round(apply(mat, 1, median, na.rm = TRUE), 3),
  TPM           = round(2^apply(mat, 1, median, na.rm = TRUE) - 1, 3),
  max_group_TPM = round(max_tpm, 3),
  pct_zero      = round(100 * rowMeans(mat == 0, na.rm = TRUE), 1),
  row.names     = NULL) %>%
  mutate(status = ifelse(max_group_TPM >= DETECT_TPM,
                         "detected", "not detected")) %>%
  arrange(desc(max_group_TPM))

if (length(missing_genes))
  det <- bind_rows(det, data.frame(
    gene = missing_genes, biotype = NA, median = NA, TPM = NA,
    max_group_TPM = NA, pct_zero = NA,
    status = ifelse(missing_genes == "GSTT1",
                    "absent - GRCh38 primary assembly (alt scaffold NT_187633.1 only)",
                    "absent from annotation")))

print(det)
write.csv(det, file.path(OUT, paste0("detection_", SOURCE, ".csv")), row.names = FALSE)


# =============================================================================
#  TIDY
# =============================================================================
df <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(clin[, c("sample","group")], by = "sample") %>%
  pivot_longer(-c(sample, group), names_to = "gene", values_to = "expr")

n_tab <- df %>% distinct(sample, group) %>% count(group) %>%
  mutate(lab = paste0(group, "\n(n=", n, ")"))
df$glab <- factor(n_tab$lab[match(df$group, n_tab$group)], levels = n_tab$lab)

BASE <- c("GTEx Normal"="#4A6572","TCGA Adjacent Normal"="#9AA5AC",
          "Adjacent Normal"="#9AA5AC","Basal"="#C0392B",
          "Her2"="#8E7CC3","LumA"="#2E86AB","LumB"="#E8A33D")
PAL <- setNames(BASE[as.character(n_tab$group)], n_tab$lab)

# Genes ordered by expression, so the panel reads as a hierarchy
ord <- det$gene[!is.na(det$max_group_TPM)][order(-det$max_group_TPM[!is.na(det$max_group_TPM)])]
df$gene <- factor(df$gene, levels = ord)
df <- df[!is.na(df$gene) & is.finite(df$expr), ]


# =============================================================================
#  FIGURE 1 -- master panel, every gene
# =============================================================================
lowlab <- df %>% group_by(gene) %>% summarise(y = max(expr, na.rm = TRUE), .groups="drop") %>%
  left_join(det[, c("gene","status")], by = "gene") %>%
  filter(status != "detected")

p <- ggplot(df, aes(glab, expr, fill = glab))
if (SHOW_VIOLIN)
  p <- p + geom_violin(scale = "width", width = .85, linewidth = .2,
                       alpha = .35, colour = NA)
p <- p +
  geom_boxplot(width = if (SHOW_VIOLIN) .28 else .62, outlier.size = .2,
               outlier.alpha = .25, linewidth = .28) +
  geom_text(data = lowlab, aes(x = 1.5, y = y, label = toupper(status)),
            inherit.aes = FALSE, size = 2.1, colour = "#B03A2E",
            fontface = "bold", hjust = 0, vjust = 1) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = expression(log[2](TPM + 1)),
       title = "Glutathione S-transferase family expression in breast tissue and TCGA-BRCA",
       subtitle = if (SOURCE == "toil")
         "GTEx and TCGA uniformly requantified (UCSC Toil recompute); panels ordered by median expression"
       else "GDC STAR-Counts, TPM; panels ordered by median expression") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.3),
        strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, paste0("Fig_panel_", SOURCE, ".png")), p,
       width = 10.5, height = 13, dpi = 300)
ggsave(file.path(OUT, paste0("Fig_panel_", SOURCE, ".pdf")), p,
       width = 10.5, height = 13)


# =============================================================================
#  FIGURE 2 -- heatmap
# =============================================================================
hm <- df %>% group_by(gene, group) %>%
  summarise(m = median(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = m) %>%
  column_to_rownames("gene") %>% as.matrix()
hm <- hm[rowSums(is.na(hm)) == 0 & apply(hm, 1, sd) > 0, , drop = FALSE]

png(file.path(OUT, paste0("Fig_heatmap_", SOURCE, ".png")),
    width = 1500, height = 2100, res = 300)
pheatmap(hm, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         border_color = "grey85", fontsize = 9,
         main = "GST family\nmedian log2(TPM+1), z-scored within gene")
dev.off()

write.csv(round(hm, 3), file.path(OUT, paste0("group_medians_", SOURCE, ".csv")))


# =============================================================================
#  FIGURE 3 -- one file per gene
# =============================================================================
starz <- function(p) if (is.na(p)) "" else
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else
    if (p < 1e-2) "**" else if (p < .05) "*" else "ns"

for (g in levels(df$gene)) {

  d  <- df[df$gene == g, ]
  st <- det[det$gene == g, ]

  ref <- d$expr[d$group == REF]
  ann <- lapply(setdiff(levels(d$group), REF), function(gr) {
    v <- d$expr[d$group == gr]
    if (length(v) < 3 || length(ref) < 3) return(NULL)
    data.frame(glab = unique(d$glab[d$group == gr]),
               y = max(d$expr, na.rm = TRUE) * 1.04,
               lab = starz(suppressWarnings(
                 wilcox.test(v, ref, exact = FALSE)$p.value)))
  }) %>% bind_rows()

  q <- ggplot(d, aes(glab, expr, fill = glab))
  if (SHOW_VIOLIN)
    q <- q + geom_violin(scale = "width", width = .85, alpha = .35,
                         colour = NA, linewidth = .2)
  q <- q +
    geom_boxplot(width = if (SHOW_VIOLIN) .3 else .6, outlier.size = .5,
                 outlier.alpha = .3, linewidth = .35) +
    { if (nrow(ann)) geom_text(data = ann, aes(glab, y, label = lab),
                               inherit.aes = FALSE, size = 3.2) } +
    scale_fill_manual(values = PAL, guide = "none") +
    labs(x = NULL, y = expression(log[2](TPM + 1)), title = g,
         subtitle = sprintf("median %.2f TPM | zero in %.1f%% of samples | %s%s",
                            st$TPM, st$pct_zero, st$status,
                            ifelse(!is.na(st$biotype) &&
                                     st$biotype != "protein_coding",
                                   paste0(" | ", st$biotype), "")),
         caption = paste0("Wilcoxon vs ", REF,
                          "; * <.05  ** <.01  *** <.001  **** <.0001")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 8.5, colour = "grey35"),
          plot.caption = element_text(size = 7, colour = "grey45"),
          panel.grid.minor = element_blank())

  ggsave(file.path(PERG, paste0(g, ".png")), q, width = 5.2, height = 4.6, dpi = 300)
}
message("Per-gene figures: ", length(list.files(PERG)), " written to ", PERG)


# =============================================================================
#  STATS
# =============================================================================
cliffs <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))$statistic
  as.numeric(2 * w / (length(x) * length(y)) - 1)
}
TUM <- c("Basal","Her2","LumA","LumB")

res <- df %>% group_by(gene) %>%
  summarise(med_ref = median(expr[group == REF], na.rm = TRUE),
            med_tum = median(expr[group %in% TUM], na.rm = TRUE),
            log2FC  = med_tum - med_ref,
            delta   = cliffs(expr[group %in% TUM], expr[group == REF]),
            p_tvn   = suppressWarnings(wilcox.test(expr[group %in% TUM],
                                                   expr[group == REF],
                                                   exact = FALSE)$p.value),
            kruskal_p = kruskal.test(expr[group %in% TUM] ~
                                       droplevels(group[group %in% TUM]))$p.value,
            subtype_spread = diff(range(tapply(expr[group %in% TUM],
                                               droplevels(group[group %in% TUM]),
                                               median))),
            .groups = "drop") %>%
  mutate(FDR_tvn = p.adjust(p_tvn, "BH"),
         FDR_kw  = p.adjust(kruskal_p, "BH"),
         effect  = cut(abs(delta), c(-Inf,.15,.33,.47,Inf),
                       labels = c("negligible","small","medium","large"))) %>%
  left_join(det[, c("gene","status","TPM")], by = "gene") %>%
  arrange(desc(abs(log2FC)))

print(as.data.frame(res), digits = 3)
write.csv(res, file.path(OUT, paste0("stats_", SOURCE, ".csv")), row.names = FALSE)

writeLines(c(
  paste("Run:", Sys.time()), paste("Source:", SOURCE),
  "Expression as log2(TPM + 1).",
  if (SOURCE == "toil")
    "UCSC Toil recompute; back-transformed from log2(TPM+0.001) before re-logging on +1."
  else "GDC STAR-Counts, tpm_unstrand assay.",
  "Ensembl version suffixes stripped before mapping; protein_coding preferred on ties.",
  paste("Genes screened:", length(GST),
        "| detected:", sum(det$status == "detected", na.rm = TRUE),
        "| low:", sum(det$status == "low", na.rm = TRUE),
        "| not detected:", sum(det$status == "not detected", na.rm = TRUE)),
  paste("Absent from annotation:", paste(missing_genes, collapse = ", ")),
  "GSTT1 is present only on GRCh38 alternate scaffold NT_187633.1 and cannot be",
  "quantified against the primary assembly (Genome Biology 2022, 23:265).",
  "", capture.output(print(table(clin$group))),
  "", capture.output(sessionInfo())
), file.path(OUT, paste0("provenance_", SOURCE, ".txt")))

message("\nAll figures written to ", normalizePath(OUT))
