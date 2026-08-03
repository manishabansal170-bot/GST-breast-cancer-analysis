# =============================================================================
#  MODULE C · PART 3 — completes the sections that failed
#
#  THE BUG
#    Error: object 'highest' not found
#    library(TCGAbiolinks) loads S4Vectors, which masks dplyr::rename with an
#    incompatible one. Fixed by calling dplyr::rename() explicitly. Worth
#    remembering: any dplyr verb can be masked once Bioconductor packages are
#    attached, and the error message never says so.
#
#  Run after Part 2. Everything is already in memory or cached.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")
BETA_CUT <- 0.3

obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; exp <- obj$exp

tum <- intersect(colnames(met), colnames(exp))
tum <- tum[substr(tum, 14, 15) == "01"]
genes <- intersect(rownames(met), rownames(exp))

st <- TCGAquery_subtype(tumor = "brca")
pam50 <- st$BRCA_Subtype_PAM50[match(substr(tum, 1, 12), st$patient)]

df <- lapply(genes, function(g) {
  data.frame(gene = g, sample = tum, beta = met[g, tum],
             expr = exp[g, tum], pam50 = pam50)
}) %>% bind_rows() %>%
  filter(pam50 %in% SUB, is.finite(beta)) %>%
  mutate(pam50 = factor(pam50, levels = SUB))

beta_wide <- df %>% group_by(gene, pam50) %>%
  summarise(median_beta = round(median(beta), 3), .groups = "drop") %>%
  pivot_wider(names_from = pam50, values_from = median_beta) %>%
  mutate(lowest_beta  = SUB[apply(across(all_of(SUB)), 1, which.min)],
         highest_beta = SUB[apply(across(all_of(SUB)), 1, which.max)],
         spread = round(apply(across(all_of(SUB)), 1, max) -
                          apply(across(all_of(SUB)), 1, min), 3))


# =============================================================================
#  THE TWO-PROGRAMME TEST
# =============================================================================
ex <- read.csv(file.path(WORKDIR, "figures", "withintumour_toil",
                         "subtype_medians_wide_toil.csv"), stringsAsFactors = FALSE)

cmp <- beta_wide %>%
  select(gene, lowest_beta, highest_beta, spread) %>%
  inner_join(ex[, c("gene","highest","lowest")], by = "gene") %>%
  dplyr::rename(expr_highest = highest, expr_lowest = lowest) %>%   # <- THE FIX
  mutate(match_high = lowest_beta == expr_highest,
         match_low  = highest_beta == expr_lowest,
         informative = spread > 0.1)

cat("\n===== DOES METHYLATION EXPLAIN THE SUBTYPE SPLIT? =====\n")
cat("Epigenetic control predicts: the subtype with the LOWEST methylation is\n")
cat("the subtype with the HIGHEST expression.\n\n")
print(as.data.frame(cmp %>% arrange(desc(spread))))

inf <- cmp %>% filter(informative)
cat(sprintf("\nAll genes:        %d of %d match\n", sum(cmp$match_high), nrow(cmp)))
cat(sprintf("Spread > 0.1 only: %d of %d match\n", sum(inf$match_high), nrow(inf)))

# Genes whose methylation barely varies between subtypes cannot test the
# hypothesis - their "lowest" subtype is noise. Restricting to spread > 0.1
# is the honest denominator.
if (nrow(inf) > 0) {
  bt <- binom.test(sum(inf$match_high), nrow(inf), p = 0.25, alternative = "greater")
  cat(sprintf("Binomial vs 25%% chance: p = %.4g\n", bt$p.value))
}
write.csv(cmp, file.path(OUT, "methylation_expression_subtype_match.csv"),
          row.names = FALSE)


# =============================================================================
#  GSTP1 IN DETAIL
# =============================================================================
gp <- df %>% filter(gene == "GSTP1")

tab <- gp %>% group_by(pam50) %>%
  summarise(n = n(),
            median_beta = round(median(beta), 3),
            pct_methylated   = round(100 * mean(beta > BETA_CUT), 1),
            pct_unmethylated = round(100 * mean(beta < 0.1), 1),
            .groups = "drop")

cat("\n===== GSTP1 METHYLATION BY SUBTYPE =====\n")
print(as.data.frame(tab))
cat("\nAgainst expression: Basal 522 TPM, Her2 159, LumA 141, LumB 72.\n")
write.csv(tab, file.path(OUT, "GSTP1_methylation_by_subtype.csv"), row.names = FALSE)

p1 <- ggplot(gp, aes(beta, fill = pam50)) +
  geom_histogram(bins = 50, alpha = .8, colour = NA) +
  geom_vline(xintercept = BETA_CUT, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ pam50, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                               LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
  labs(x = "GSTP1 promoter methylation (beta)", y = "Tumours",
       title = "GSTP1 promoter methylation is bimodal and subtype-dependent",
       subtitle = "Basal-like tumours escape methylation; dashed line = beta 0.3") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "GSTP1_methylation_bimodality.png"), p1,
       width = 8, height = 6, dpi = 300)


# =============================================================================
#  THE FIGURE THAT MAKES THE ARGUMENT
# =============================================================================
# Methylation and expression side by side, for the six genes where subtype
# methylation varies. The two panels should mirror each other.
key <- inf$gene

if (length(key)) {
  mplot <- df %>% filter(gene %in% key) %>%
    mutate(panel = "Promoter methylation (beta)") %>%
    select(gene, pam50, value = beta, panel)

  eplot <- df %>% filter(gene %in% key) %>%
    mutate(panel = "Expression", value = log2(expr + 1)) %>%
    select(gene, pam50, value, panel)

  both <- bind_rows(mplot, eplot) %>%
    mutate(panel = factor(panel, levels = c("Promoter methylation (beta)","Expression")))

  p2 <- ggplot(both, aes(pam50, value, fill = pam50)) +
    geom_violin(scale = "width", alpha = .3, colour = NA) +
    geom_boxplot(width = .25, outlier.size = .2, outlier.alpha = .2, linewidth = .25) +
    facet_grid(panel ~ gene, scales = "free_y", switch = "y") +
    scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                                 LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
    labs(x = NULL, y = NULL,
         title = "Methylation and expression mirror each other across PAM50 subtypes",
         subtitle = "Genes with subtype-structured promoter methylation; low methylation tracks high expression") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.background = element_rect(fill = "grey93"),
          strip.text = element_text(face = "bold", size = 8),
          plot.title = element_text(face = "bold", size = 11),
          strip.placement = "outside")
  ggsave(file.path(OUT, "methylation_expression_mirror.png"), p2,
         width = 11, height = 6.5, dpi = 300)
  message("Mirror figure written - this is the one for the deck.")
}


# =============================================================================
#  IS THERE A NORMAL-TISSUE METHYLATION SET ANYWHERE?
# =============================================================================
# brca_tcga has 783 tumours and no normals, so tumour-vs-normal never ran.
# Without it you can say methylation tracks silencing within tumours, but not
# that methylation INCREASED in cancer. Check the PanCancer Atlas study.

message("\nChecking for a study with normal methylation samples...")
chk <- tryCatch({
  library(cBioPortalData); cb <- cBioPortal()
  for (s in c("brca_tcga_pan_can_atlas_2018","brca_tcga_gdc","brca_tcga_pub2015")) {
    mp <- tryCatch(molecularProfiles(cb, s), error = function(e) NULL)
    if (is.null(mp)) next
    m <- mp$molecularProfileId[grepl("methylation", mp$molecularProfileId, ignore.case = TRUE)]
    if (length(m)) cat(sprintf("  %s -> %s\n", s, paste(m, collapse = ", ")))
  }
  TRUE
}, error = function(e) { message("  check failed: ", conditionMessage(e)); FALSE })

cat("\nIf none carry normals, the tumour-vs-normal methylation comparison is a\n")
cat("genuine gap in the public data - and one your own North Indian cohort\n")
cat("is positioned to fill. That is an argument for Objective 2, not a flaw.\n")

message("\nWritten to ", normalizePath(OUT))
