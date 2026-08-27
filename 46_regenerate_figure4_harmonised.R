# =============================================================================
#  46_regenerate_figure4_harmonised.R
#
#  THE PROBLEM
#    Table 4's caption states that all expression values in the manuscript now
#    derive from a single quantification. Figure 4's legend says its expression
#    row comes from a differently annotated matrix whose values "should not be
#    read across figures". Both cannot be true, and a reviewer has flagged the
#    contradiction.
#
#    The variance partition was moved to the harmonised GRCh38 matrix in an
#    earlier step; this figure was not, because it is drawn by a different
#    script that reads the expression matrix bundled with the methylation
#    object. Regenerating it on the harmonised matrix removes the contradiction
#    rather than qualifying it away, and makes every expression value in the
#    manuscript comparable to every other.
#
#  WHAT CHANGES AND WHAT DOES NOT
#    The methylation row is untouched: it comes from the same beta values as
#    before. Only the expression row is recomputed, on TPM from the GRCh38
#    requantification, in the same tumours.
#
#    The six genes shown are those with subtype-variable methylation, selected
#    as before by the spread of median beta across subtypes, so the panel
#    membership is not being changed to suit the new values.
#
#  A CHECK WORTH MAKING
#    The legend claims the least-methylated subtype is the most-expressed in
#    five of the six genes. That count is recomputed here on the new values
#    rather than carried over, because the expression row has changed.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr", "tidyr", "ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. THE TWO LAYERS, ON ONE ANNOTATION
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; colnames(met) <- substr(colnames(met), 1, 15)

gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
LIN <- 2^gm(raw) - 0.001; LIN[LIN < 0] <- 0          # TPM, harmonised GRCh38
colnames(LIN) <- substr(colnames(LIN), 1, 15)

st <- TCGAquery_subtype(tumor = "brca")
s  <- intersect(colnames(met), colnames(LIN))
s  <- s[substr(s, 13, 15) == "-01"]
g  <- st$BRCA_Subtype_PAM50[match(substr(s, 1, 12), st$patient)]
ok <- g %in% SUB; s <- s[ok]; grp <- factor(g[ok], levels = SUB)

cat("Tumours with methylation and harmonised expression:", length(s), "\n")
print(table(grp))

genes <- intersect(rownames(met), rownames(LIN))
cat("Genes on both layers:", length(genes), "\n\n")


# =============================================================================
# 2. WHICH SIX GENES, BY THE PUBLISHED RULE
# =============================================================================
# Selection is on the spread of median beta across subtypes, as before, so the
# panel is not being reselected to suit the new expression values.
spread <- sapply(genes, function(gn) {
  m <- tapply(as.numeric(met[gn, s]), grp, median, na.rm = TRUE)
  diff(range(m, na.rm = TRUE))
})
key <- names(sort(spread, decreasing = TRUE))[1:6]
cat("Genes shown, by spread of median beta across subtypes:\n")
print(round(sort(spread, decreasing = TRUE)[1:6], 3))


# =============================================================================
# 3. DOES THE LEAST-METHYLATED SUBTYPE REMAIN THE MOST-EXPRESSED
# =============================================================================
cat("\n", strrep("=", 76), "\nRECOMPUTED ON THE HARMONISED MATRIX\n",
    strrep("=", 76), "\n", sep = "")
chk <- do.call(rbind, lapply(key, function(gn) {
  mb <- tapply(as.numeric(met[gn, s]), grp, median, na.rm = TRUE)
  ex <- tapply(as.numeric(LIN[gn, s]), grp, median, na.rm = TRUE)
  data.frame(gene = gn,
             lowest_beta = names(which.min(mb)),
             highest_expr = names(which.max(ex)),
             match = names(which.min(mb)) == names(which.max(ex)),
             stringsAsFactors = FALSE)
}))
print(chk, row.names = FALSE)
cat(sprintf("\nLeast-methylated is most-expressed in %d of the %d genes shown.\n",
            sum(chk$match), nrow(chk)))
cat("The legend must state whatever this count is, not the previous one.\n")


# =============================================================================
# 4. THE FIGURE
# =============================================================================
mp <- lapply(key, function(gn) data.frame(
  gene = gn, pam50 = grp, value = as.numeric(met[gn, s]),
  panel = "Promoter methylation (beta)")) %>% bind_rows()
ep <- lapply(key, function(gn) data.frame(
  gene = gn, pam50 = grp, value = log2(as.numeric(LIN[gn, s]) + 1),
  panel = "Expression, log2(TPM + 1)")) %>% bind_rows()

both <- bind_rows(mp, ep) %>%
  filter(is.finite(value)) %>%
  mutate(panel = factor(panel, levels = c("Promoter methylation (beta)",
                                          "Expression, log2(TPM + 1)")),
         gene = factor(gene, levels = key))

# Threshold lines belong only on the methylation row; carrying the panel factor
# confines them there, where a bare geom_hline would draw across both.
thr <- data.frame(panel = factor("Promoter methylation (beta)",
                                 levels = levels(both$panel)),
                  y = c(0.1, 0.3))

p <- ggplot(both, aes(pam50, value, fill = pam50)) +
  geom_hline(data = thr, aes(yintercept = y), linetype = "dashed",
             colour = "grey45", linewidth = .35) +
  geom_violin(scale = "width", alpha = .3, colour = NA) +
  geom_boxplot(width = .25, outlier.size = .2, outlier.alpha = .2, linewidth = .25) +
  facet_grid(panel ~ gene, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                               LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Promoter methylation and expression by PAM50 subtype",
       subtitle = sprintf(paste("Expression from the uniformly requantified GRCh38 data, the same",
                                "matrix as Table 2 and Figure 1 (n = %d tumours).",
                                "Dashed lines mark beta 0.1 and 0.3."), length(s))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8),
        strip.placement = "outside")

ggsave(file.path(OUT, "Figure4_harmonised.png"), p, width = 11, height = 6.5,
       dpi = 300, bg = "white")
ggsave(file.path(OUT, "Figure4_harmonised.tiff"), p, width = 11, height = 6.5,
       dpi = 300, bg = "white", compression = "lzw")

cat("\nWritten:\n")
cat("  ", normalizePath(file.path(OUT, "Figure4_harmonised.png")), "\n")
cat("  ", normalizePath(file.path(OUT, "Figure4_harmonised.tiff")), "\n")

write.csv(chk, file.path(OUT, "figure4_direction_check.csv"), row.names = FALSE)
writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d", length(s)),
             sprintf("Genes shown: %s", paste(key, collapse = ", ")),
             sprintf("Least-methylated is most-expressed in %d of %d",
                     sum(chk$match), nrow(chk)),
             "Expression row regenerated on the harmonised GRCh38 matrix.",
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_figure4_harmonised.txt"))

cat("\nReplace Figure 4 in the manuscript with the TIFF above, and rewrite the\n")
cat("legend to state the harmonised matrix and the recomputed count.\n")
