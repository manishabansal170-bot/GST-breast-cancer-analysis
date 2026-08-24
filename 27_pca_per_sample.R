# =============================================================================
#  27_pca_per_sample.R
#
#  THE OBJECTION
#    Hierarchical clustering of a 5 x 20 group-median matrix was used to
#    identify two blocks. A principal component analysis of THE SAME matrix
#    was then described as having "recovered the same structure independently".
#
#    Two unsupervised decompositions of one heavily aggregated matrix are not
#    independent evidence. Aggregating 1,041 tumours into five group medians
#    discards essentially all within-group variation, and any structure the
#    clustering found will reappear in the PCA because it is the same twenty
#    numbers arranged differently.
#
#    The referee also noted that PC3 carried 18.9%, nearly as much as PC2's
#    25.8%, and loaded on a dimension the model does not describe. Dismissing
#    it because HER2-enriched is the smallest group was post hoc.
#
#  WHAT THIS SCRIPT DOES
#    Runs the decomposition on per-sample data: 1,041 tumours by the analysed
#    genes, with no aggregation. This asks whether the two-dimensional
#    structure exists among tumours, which is the claim the manuscript makes,
#    rather than among five averages.
#
#    Component count is assessed against a parallel analysis null rather than
#    by eye, so the number of real dimensions is decided by a criterion
#    stated in advance.
#
#  READING THE RESULT
#    The manuscript's model predicts two interpretable components: one
#    contrasting tumours with adjacent normal, one contrasting basal-like with
#    luminal. If parallel analysis retains more than two, the phenotype has
#    more structure than the model describes, and that must be reported rather
#    than explained away. A third component loading on HER2 would be a real
#    finding, not an inconvenience.
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

set.seed(20260824)
SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. PER-SAMPLE MATRIX, TUMOURS AND ADJACENT NORMAL
# =============================================================================
get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^get_matrix(raw) - 0.001; lin[lin < 0] <- 0
E <- log2(lin + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]], study = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

adjn <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                             grepl("breast", tissue, ignore.case = TRUE),
                             grepl("Solid Tissue Normal", stype)) %>%
  mutate(group = "AdjNormal")
tum <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                            grepl("breast", tissue, ignore.case = TRUE),
                            grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

meta <- bind_rows(adjn %>% dplyr::select(sample, group),
                  tum  %>% dplyr::select(sample, group))
meta <- meta[meta$sample %in% colnames(E), ]
X <- E[, meta$sample, drop = FALSE]
X <- X[apply(X, 1, function(v) sd(v, na.rm = TRUE) > 1e-8), , drop = FALSE]

cat("Samples:", ncol(X), " Genes:", nrow(X), "\n"); print(table(meta$group))

# Samples as rows, genes as columns; scale so abundance does not dominate
M <- t(X)
M <- M[, apply(M, 2, function(v) all(is.finite(v))), drop = FALSE]
M <- scale(M)
cat("Matrix for decomposition:", nrow(M), "samples x", ncol(M), "genes\n")


# =============================================================================
# 2. PRINCIPAL COMPONENTS
# =============================================================================
pc <- prcomp(M, center = FALSE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)

cat("\n", strrep("=", 70), "\nVARIANCE EXPLAINED, PER-SAMPLE DATA\n", strrep("=", 70), "\n", sep = "")
for (i in 1:min(8, length(ve)))
  cat(sprintf("  PC%-2d %5.1f%%   cumulative %5.1f%%\n", i, ve[i], sum(ve[1:i])))


# =============================================================================
# 3. HOW MANY COMPONENTS ARE REAL
# =============================================================================
# Parallel analysis: decompose matrices of the same size built by permuting
# each gene independently, which destroys between-gene structure while keeping
# each gene's distribution. A component is retained if it exceeds the 95th
# centile of the corresponding null component.
cat("\n", strrep("=", 70), "\nPARALLEL ANALYSIS\n", strrep("=", 70), "\n", sep = "")
NPERM <- 200
null <- replicate(NPERM, {
  Mp <- apply(M, 2, sample)
  s <- prcomp(Mp, center = FALSE, scale. = FALSE)$sdev^2
  100 * s / sum(s) })
thr <- apply(null, 1, quantile, 0.95)

keep <- which(ve > thr)
keep <- if (length(keep)) seq_len(max(which(cumprod(ve > thr) == 1))) else integer(0)
cat(sprintf("%-6s %10s %10s %s\n", "", "observed", "null 95%", ""))
for (i in 1:min(6, length(ve)))
  cat(sprintf("  PC%-2d %8.1f%% %9.1f%%   %s\n", i, ve[i], thr[i],
              if (ve[i] > thr[i]) "retained" else "not retained"))
cat(sprintf("\nComponents retained by parallel analysis: %d\n", length(keep)))


# =============================================================================
# 4. WHAT EACH COMPONENT REPRESENTS
# =============================================================================
cat("\n", strrep("=", 70), "\nCOMPONENT MEANING\n", strrep("=", 70), "\n", sep = "")
sc <- as.data.frame(pc$x[, 1:min(4, ncol(pc$x)), drop = FALSE])
sc$group <- factor(meta$group, levels = c("AdjNormal", SUB))

for (i in seq_len(min(4, ncol(pc$x)))) {
  nm <- paste0("PC", i)
  mu <- tapply(sc[[nm]], sc$group, mean, na.rm = TRUE)
  cat(sprintf("\n%s (%.1f%% of variance) group means:\n", nm, ve[i]))
  print(round(mu, 2))
  k <- kruskal.test(sc[[nm]] ~ sc$group)
  cat(sprintf("  Kruskal-Wallis across groups: p = %.3g\n", k$p.value))
  ld <- sort(pc$rotation[, i])
  cat("  most negative loadings:", paste(names(head(ld, 3)), collapse = ", "), "\n")
  cat("  most positive loadings:", paste(names(tail(ld, 3)), collapse = ", "), "\n")
}


# =============================================================================
# 5. VERDICT
# =============================================================================
cat("\n", strrep("=", 70), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 70), "\n", sep = "")
if (length(keep) == 2) {
  cat("Parallel analysis retains exactly two components, matching the model.\n")
  cat("Report these per-sample values in place of the group-median decomposition,\n")
  cat("and state the retention criterion.\n")
} else if (length(keep) > 2) {
  cat(sprintf("Parallel analysis retains %d components, more than the model describes.\n", length(keep)))
  cat("Report this honestly: the two-programme model captures the first two\n")
  cat("dimensions but the phenotype has further structure. Describe what the\n")
  cat("additional components load on rather than dismissing them.\n")
} else {
  cat(sprintf("Parallel analysis retains %d component.\n", length(keep)))
  cat("A single dimension would not support a two-axis model, and the framing\n")
  cat("would need substantial revision.\n")
}
cat("\nIn either case the group-median PCA should be removed from the Results,\n")
cat("since decomposing the same matrix the clustering used is not independent\n")
cat("confirmation and the referee is right to say so.\n")

write.csv(data.frame(component = paste0("PC", seq_along(ve)),
                     variance_pct = round(ve, 2),
                     null_95_pct = round(c(thr, rep(NA, max(0, length(ve) - length(thr)))), 2)),
          file.path(OUT, "pca_per_sample_variance.csv"), row.names = FALSE)
write.csv(round(pc$rotation[, 1:min(4, ncol(pc$rotation))], 4),
          file.path(OUT, "pca_per_sample_loadings.csv"))


# =============================================================================
# 6. FIGURE
# =============================================================================
p <- ggplot(sc, aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.1, alpha = .6) +
  stat_ellipse(level = 0.68, linewidth = .5) +
  scale_colour_manual(values = c(AdjNormal = "#7F8C8D", Basal = "#C0392B",
                                 Her2 = "#8E7CC3", LumA = "#2E86AB",
                                 LumB = "#E8A33D"), name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "GST family expression decomposed at sample level",
       subtitle = sprintf("%d samples, %d genes, no aggregation; %d components retained by parallel analysis",
                          nrow(M), ncol(M), length(keep))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "pca_per_sample.png"), p, width = 7, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT, "pca_per_sample.tiff"), p, width = 7, height = 6, dpi = 300,
       bg = "white", compression = "lzw")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Samples: %d  Genes: %d", nrow(M), ncol(M)),
             sprintf("Variance: %s", paste(sprintf("PC%d %.1f%%", 1:min(5,length(ve)),
                                                   ve[1:min(5,length(ve))]), collapse = "  ")),
             sprintf("Components retained by parallel analysis: %d", length(keep)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_pca_per_sample.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
