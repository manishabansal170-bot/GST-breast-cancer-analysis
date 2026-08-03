# =============================================================================
#  PER-GENE PUBLICATION FIGURES  --  log2(TPM + 1)
#
#  Reproduces the GSTA1-style figure (significance brackets + global
#  Kruskal-Wallis) for every gene, on a properly normalised axis.
#
#  WHY THE NUMBERS WILL CHANGE
#  ---------------------------
#  Your existing panel is log2(counts + 1). Raw counts scale with library
#  size, so a deeply sequenced sample shows higher counts for every gene
#  regardless of biology. Between-sample comparison on that axis is not
#  valid, which is why this panel disagreed with your TPM panels on subtype
#  ordering. TPM divides out both library size and gene length.
#
#  Expect: smaller y range, smaller apparent differences, and possibly a
#  different subtype ordering. The TPM version is the correct one.
# =============================================================================

SOURCE         <- "toil"    # "toil" or "gdc"
INCLUDE_NORMAL <- TRUE      # add the normal group(s); strongly recommended
ALL_PAIRWISE   <- FALSE     # TRUE = every pair (busy), FALSE = each vs reference
SHOW_VIOLIN    <- TRUE

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures", paste0("pergene_pub_", SOURCE))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

for (p in c("dplyr","tidyr","tibble","ggplot2","ggpubr"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(ggpubr)
})

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")


# ---------------------------------------------------------------------------
#  DATA  (same preparation as the master script)
# ---------------------------------------------------------------------------
if (SOURCE == "gdc") {
  suppressPackageStartupMessages({
    library(TCGAbiolinks); library(SummarizedExperiment) })
  se <- readRDS(file.path(CACHE, "brca_se.rds"))

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
  mat <- mat[hit[!is.na(hit)], , drop = FALSE]
  rownames(mat) <- names(hit)[!is.na(hit)]

  clin <- as.data.frame(colData(se))
  clin$sample <- colnames(se); clin$patient <- substr(clin$sample, 1, 12)
  pc <- grep("PAM50", colnames(clin), value = TRUE)[1]
  clin$pam50 <- if (!is.na(pc)) as.character(clin[[pc]]) else {
    st <- TCGAquery_subtype(tumor = "brca")
    st$BRCA_Subtype_PAM50[match(clin$patient, st$patient)] }
  clin <- clin[clin$sample_type %in% c("Primary Tumor","Solid Tissue Normal"), ]
  clin$group <- ifelse(clin$sample_type == "Solid Tissue Normal",
                       "Adjacent Normal", clin$pam50)
  LEV <- c("Adjacent Normal","Basal","Her2","LumA","LumB")
  REF <- "Adjacent Normal"

} else {
  suppressPackageStartupMessages({
    library(UCSCXenaTools); library(TCGAbiolinks) })
  raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
  lin <- 2^raw - 0.001; lin[lin < 0] <- 0
  mat <- log2(lin + 1)

  ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
  cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
  ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                         study = .data[[cl("study")]],
                         tissue = .data[[cl("primary disease or tissue")]],
                         stype = .data[[cl("sample_type")]])
  br <- ph %>% filter(grepl("breast", tissue, ignore.case = TRUE))
  st <- TCGAquery_subtype(tumor = "brca")
  clin <- bind_rows(
    br %>% filter(toupper(study) == "GTEX") %>% transmute(sample, group = "GTEx Normal"),
    br %>% filter(toupper(study) == "TCGA", grepl("Solid Tissue Normal", stype)) %>%
      transmute(sample, group = "TCGA Adjacent Normal"),
    br %>% filter(toupper(study) == "TCGA", grepl("Primary Tumor", stype)) %>%
      mutate(group = st$BRCA_Subtype_PAM50[match(substr(sample,1,12), st$patient)]) %>%
      filter(group %in% c("Basal","Her2","LumA","LumB")) %>% select(sample, group))
  clin$patient <- ifelse(grepl("^TCGA", clin$sample),
                         substr(clin$sample, 1, 12), clin$sample)
  LEV <- c("GTEx Normal","TCGA Adjacent Normal","Basal","Her2","LumA","LumB")
  REF <- "TCGA Adjacent Normal"
}

if (!INCLUDE_NORMAL) {
  LEV <- setdiff(LEV, c("GTEx Normal","TCGA Adjacent Normal","Adjacent Normal"))
  REF <- "Basal"
}

clin <- clin[clin$group %in% LEV, ] %>% arrange(sample) %>%
  distinct(patient, group, .keep_all = TRUE)
clin$group <- factor(clin$group, levels = LEV)

# A filter that fails to match produces an empty group and a figure with a
# blank column, which looks like real data. Stop instead.
if (any(table(clin$group) == 0))
  stop("A sample group is empty. Check the study and sample_type filters.")
clin <- clin[clin$sample %in% colnames(mat), ]
mat  <- mat[, clin$sample, drop = FALSE]

df <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(clin[, c("sample","group")], by = "sample") %>%
  pivot_longer(-c(sample, group), names_to = "gene", values_to = "expr") %>%
  filter(!is.na(group), is.finite(expr))

n_tab <- df %>% distinct(sample, group) %>% count(group) %>%
  mutate(lab = paste0(group, "\n(n=", n, ")"))
df$glab <- factor(n_tab$lab[match(df$group, n_tab$group)], levels = n_tab$lab)

BASE <- c("GTEx Normal"="#4A6572","TCGA Adjacent Normal"="#9AA5AC",
          "Adjacent Normal"="#9AA5AC","Basal"="#C0392B",
          "Her2"="#8E7CC3","LumA"="#2E86AB","LumB"="#E8A33D")
PAL <- setNames(BASE[as.character(n_tab$group)], n_tab$lab)

cat("\nSamples per group:\n"); print(table(clin$group)); cat("\n")


# ---------------------------------------------------------------------------
#  COMPARISON LIST
# ---------------------------------------------------------------------------
lv     <- levels(df$glab)
ref_lab <- n_tab$lab[n_tab$group == REF]

comps <- if (ALL_PAIRWISE) combn(lv, 2, simplify = FALSE) else
  lapply(setdiff(lv, ref_lab), function(x) c(ref_lab, x))

cliffs <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE))$statistic
  as.numeric(2 * w / (length(x) * length(y)) - 1)
}
TUM <- intersect(c("Basal","Her2","LumA","LumB"), LEV)


# ---------------------------------------------------------------------------
#  PLOT EVERY GENE
# ---------------------------------------------------------------------------
for (g in sort(unique(df$gene))) {

  d <- df[df$gene == g, ]
  if (!nrow(d)) next

  med_tpm <- 2^median(d$expr, na.rm = TRUE) - 1
  pct0    <- 100 * mean(d$expr == 0, na.rm = TRUE)
  dlt     <- cliffs(d$expr[d$group %in% TUM], d$expr[d$group == REF])

  # Bracket ladder above the data. Fixed spacing stops the overlap between
  # brackets and the Kruskal-Wallis label seen in the counts version.
  top  <- max(d$expr, na.rm = TRUE)
  step <- max(top * 0.09, 0.25)
  ys   <- top + step * seq_along(comps)
  kw_y <- top + step * (length(comps) + 1.4)

  p <- ggplot(d, aes(glab, expr, fill = glab))
  if (SHOW_VIOLIN)
    p <- p + geom_violin(scale = "width", width = .85, alpha = .3,
                         colour = NA, linewidth = .2)
  p <- p +
    geom_boxplot(width = if (SHOW_VIOLIN) .3 else .6, outlier.size = .5,
                 outlier.alpha = .3, linewidth = .4) +
    stat_compare_means(comparisons = comps, method = "wilcox.test",
                       label = "p.signif", label.y = ys,
                       tip.length = .012, size = 3.2, bracket.size = .3) +
    stat_compare_means(method = "kruskal.test", label.y = kw_y,
                       size = 3.1, colour = "grey25") +
    scale_fill_manual(values = PAL, guide = "none") +
    coord_cartesian(ylim = c(min(0, min(d$expr, na.rm = TRUE)), kw_y + step)) +
    labs(x = NULL, y = expression(log[2](TPM + 1)), title = g,
         subtitle = sprintf(
           "median %.2f TPM  |  zero in %.1f%% of samples  |  Cliff's d vs %s = %s",
           med_tpm, pct0, REF,
           ifelse(is.na(dlt), "NA", sprintf("%+.2f", dlt))),
         caption = "Wilcoxon rank-sum; ns  * <.05  ** <.01  *** <.001  **** <.0001") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          plot.title = element_text(face = "bold", size = 16, hjust = .5),
          plot.subtitle = element_text(size = 8.5, colour = "grey35", hjust = .5),
          plot.caption = element_text(size = 7, colour = "grey45"),
          axis.line = element_line(linewidth = .4))

  ggsave(file.path(OUT, paste0(g, ".png")), p, width = 5.6, height = 5.4, dpi = 300)
  ggsave(file.path(OUT, paste0(g, ".pdf")), p, width = 5.6, height = 5.4)
}

message(length(list.files(OUT, pattern = "\\.png$")),
        " figures written to ", normalizePath(OUT))


# ---------------------------------------------------------------------------
#  WHY THE STARS ARE NOT THE POINT
# ---------------------------------------------------------------------------
# At n in the hundreds to thousands, **** appears for differences far too
# small to matter. Cliff's delta in each subtitle is the honest summary:
#   |d| < 0.15 negligible | < 0.33 small | < 0.47 medium | else large
# A gene with **** and d = 0.08 is a statistically detectable non-finding.
# Report both, and lead with the effect size.

summ <- df %>% group_by(gene) %>%
  summarise(median_TPM = round(2^median(expr, na.rm = TRUE) - 1, 3),
            pct_zero   = round(100 * mean(expr == 0, na.rm = TRUE), 1),
            delta_vs_ref = round(cliffs(expr[group %in% TUM],
                                        expr[group == REF]), 3),
            kruskal_p  = signif(kruskal.test(expr ~ droplevels(group))$p.value, 3),
            .groups = "drop") %>%
  mutate(kruskal_FDR = signif(p.adjust(kruskal_p, "BH"), 3),
         effect = cut(abs(delta_vs_ref), c(-Inf,.15,.33,.47,Inf),
                      labels = c("negligible","small","medium","large"))) %>%
  arrange(desc(abs(delta_vs_ref)))

print(as.data.frame(summ))
write.csv(summ, file.path(OUT, "pergene_summary.csv"), row.names = FALSE)
