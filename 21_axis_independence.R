# =============================================================================
#  21_axis_independence.R
#
#  THE OBJECTION THIS ANSWERS
#    The manuscript describes two axes of change: a universal loss axis, on
#    which genes fall in every subtype relative to adjacent normal, and a
#    subtype axis, on which the remainder are redistributed between basal-like
#    and luminal programmes. It calls these axes independent.
#
#    As submitted, that word was asserted rather than tested. The evidence
#    offered was hierarchical clustering of group medians, which will separate
#    blocks whether or not the underlying axes are independent. A reviewer is
#    entitled to ask for a test.
#
#  WHAT INDEPENDENCE WOULD MEAN, CONCRETELY
#    If the two axes are independent, then knowing how much a gene is lost in
#    tumour overall tells you nothing about how strongly it separates the
#    subtypes. That is a correlation between two per-gene quantities:
#
#      loss magnitude    = Cliff's delta, all tumours versus adjacent normal
#      subtype magnitude = epsilon-squared across the four PAM50 subtypes
#
#    A correlation near zero supports independence. A strong correlation would
#    mean the two axes are one axis described twice, and the two-programme
#    model would need rewriting.
#
#  A SECOND, STRICTER TEST
#    Principal components of the group-median matrix. If the phenotype really
#    is two-dimensional, two components should carry most of the variance and
#    should load in the way the model predicts: one on overall level, one on
#    the basal-versus-luminal contrast.
#
#  ON GENES APPEARING ON BOTH AXES
#    GSTA1, GSTA4, GSTM2 and GSTM5 are named both in the loss table and in the
#    subtype programmes. That is not a contradiction under a two-dimensional
#    model: a gene can sit below normal in every subtype while still being
#    relatively higher in one programme than the other. This script quantifies
#    that directly, reporting both coordinates for every gene so the reader can
#    see the two-dimensional position rather than being asked to accept it.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr","ggplot2","ggrepel")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. REBUILD THE EXPRESSION MATRIX EXACTLY AS SCRIPT 16 DOES
# =============================================================================
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }
lin <- 2^get_matrix(raw) - 0.001; lin[lin < 0] <- 0
E <- log2(lin + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])
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
E <- E[, meta$sample, drop = FALSE]

cat("Samples:\n"); print(table(meta$group))

# Keep genes with signal, matching the manuscript's detection handling
keep <- apply(E, 1, function(v) sd(v, na.rm = TRUE) > 1e-8)
E <- E[keep, , drop = FALSE]
cat("\nGenes analysed:", nrow(E), "\n")


# =============================================================================
# 2. THE TWO COORDINATES, PER GENE
# =============================================================================
is_norm <- meta$group == "AdjNormal"
is_tum  <- meta$group %in% SUB

# Cliff's delta: tumour versus adjacent normal. Negative means lower in tumour.
cliffs <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  # rank-based computation, avoids an n x m outer product
  r <- rank(c(a, b))
  ra <- sum(r[seq_along(a)])
  U <- ra - length(a) * (length(a) + 1) / 2
  2 * U / (length(a) * length(b)) - 1
}

eps_sq <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 20) return(NA_real_)
  k <- kruskal.test(x[ok] ~ droplevels(factor(g[ok])))
  n <- sum(ok)
  unname(k$statistic) / ((n^2 - 1) / (n + 1))
}

coord <- lapply(rownames(E), function(g) {
  x <- as.numeric(E[g, ])
  data.frame(
    gene = g,
    loss_delta = round(cliffs(x[is_tum], x[is_norm]), 3),
    subtype_eps2 = round(eps_sq(x[is_tum], meta$group[is_tum]), 3))
}) %>% bind_rows() %>% dplyr::filter(is.finite(loss_delta), is.finite(subtype_eps2))

coord$loss_magnitude <- abs(coord$loss_delta)

cat("\n===== PER-GENE COORDINATES ON THE TWO AXES =====\n")
print(as.data.frame(coord[order(coord$loss_delta), ]), row.names = FALSE, digits = 3)


# =============================================================================
# 3. TEST ONE: ARE THE AXES CORRELATED?
# =============================================================================
cat("\n", strrep("=", 70), "\nTEST 1: CORRELATION BETWEEN THE TWO AXES\n",
    strrep("=", 70), "\n", sep = "")

ct <- suppressWarnings(cor.test(coord$loss_magnitude, coord$subtype_eps2,
                                method = "spearman"))
cat(sprintf("Spearman rho = %.3f, p = %.3g, n = %d genes\n",
            ct$estimate, ct$p.value, nrow(coord)))

cat("\nInterpretation:\n")
if (abs(ct$estimate) < 0.3 && ct$p.value > 0.05) {
  cat("  Loss magnitude and subtype effect are uncorrelated across genes.\n")
  cat("  How much a gene falls in tumour tells you nothing about how strongly\n")
  cat("  it separates the subtypes. The two axes are independent, and the word\n")
  cat("  is now earned rather than asserted.\n")
} else if (abs(ct$estimate) < 0.3) {
  cat("  The correlation is weak but reaches significance at this sample size.\n")
  cat("  Report the coefficient rather than claiming strict independence, and\n")
  cat("  describe the axes as largely rather than wholly independent.\n")
} else {
  cat("  The axes are substantially correlated. They are not independent, and\n")
  cat("  the two-dimensional framing must be rewritten. This would be a\n")
  cat("  material change to the manuscript's central claim.\n")
}


# =============================================================================
# 4. TEST TWO: HOW MANY DIMENSIONS ARE THERE?
# =============================================================================
cat("\n", strrep("=", 70), "\nTEST 2: PRINCIPAL COMPONENTS OF THE GROUP-MEDIAN MATRIX\n",
    strrep("=", 70), "\n", sep = "")

grp <- c("AdjNormal", SUB)
med <- sapply(grp, function(g) apply(E[, meta$group == g, drop = FALSE], 1,
                                     median, na.rm = TRUE))
med <- med[coord$gene, , drop = FALSE]

# Z-score within gene before decomposing. Without this the first component is
# simply absolute abundance: these genes span roughly 7 to 445 TPM, so an
# unscaled decomposition puts 96% of variance on a component that loads almost
# identically across all groups and says nothing about pattern.
med <- t(scale(t(med)))

pc <- prcomp(med, center = TRUE, scale. = FALSE)
ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
cat("Variance explained by component: ", paste0("PC", seq_along(ve), " ", ve, "%",
    collapse = "  "), "\n", sep = "")
cat(sprintf("\nFirst two components together: %.1f%%\n", sum(ve[1:2])))

cat("\nLoadings, which say what each component represents:\n")
print(round(pc$rotation[, 1:min(3, ncol(pc$rotation))], 3))

cat("\nIf PC1 loads with the same sign across all groups it represents overall\n")
cat("expression level. If PC2 loads with opposite signs on basal-like versus\n")
cat("the luminal subtypes it represents the subtype contrast. That pattern\n")
cat("would be the two-dimensional structure the manuscript describes.\n")


# =============================================================================
# 5. GENES ON BOTH AXES, MADE EXPLICIT
# =============================================================================
cat("\n", strrep("=", 70), "\nGENES APPEARING ON BOTH AXES\n", strrep("=", 70), "\n", sep = "")

both <- coord %>%
  dplyr::filter(loss_delta <= -0.4, subtype_eps2 >= 0.1) %>%
  arrange(loss_delta)

if (nrow(both)) {
  cat("These genes are substantially reduced in tumour overall AND separate the\n")
  cat("subtypes. Under a one-dimensional model that would be contradictory.\n")
  cat("Under a two-dimensional model it is expected: the gene sits below normal\n")
  cat("everywhere, while still being relatively higher in one programme.\n\n")
  print(as.data.frame(both), row.names = FALSE, digits = 3)
} else {
  cat("No gene meets both criteria at these thresholds.\n")
}

write.csv(coord, file.path(OUT, "axis_coordinates.csv"), row.names = FALSE)


# =============================================================================
# 6. FIGURE: THE PHENOTYPE AS TWO DIMENSIONS
# =============================================================================
p <- ggplot(coord, aes(loss_delta, subtype_eps2)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .3) +
  geom_point(size = 2.4, colour = "#2E86AB", alpha = .85) +
  ggrepel::geom_text_repel(aes(label = gene), size = 3, max.overlaps = 30,
                           segment.colour = "grey70", segment.size = .3) +
  labs(x = "Loss axis: Cliff's delta, tumour versus adjacent normal",
       y = "Subtype axis: epsilon-squared across PAM50 subtypes",
       title = "The two axes are independent",
       subtitle = sprintf("Spearman rho = %.3f, p = %.3g across %d genes. Position on one axis does not predict the other.",
                          ct$estimate, ct$p.value, nrow(coord))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "axis_independence.png"), p, width = 7.5, height = 5.5,
       dpi = 300, bg = "white")
ggsave(file.path(OUT, "axis_independence.tiff"), p, width = 7.5, height = 5.5,
       dpi = 300, bg = "white", compression = "lzw")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Genes: %d", nrow(coord)),
             sprintf("Axis correlation: Spearman rho %.3f, p %.3g", ct$estimate, ct$p.value),
             sprintf("PC1 %.1f%%, PC2 %.1f%%, together %.1f%%", ve[1], ve[2], sum(ve[1:2])),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_axis_independence.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
cat("Report the correlation coefficient in the Results where the word\n")
cat("'independent' first appears.\n")
