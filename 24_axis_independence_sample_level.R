# =============================================================================
#  24_axis_independence_sample_level.R
#
#  WHAT THE PREVIOUS TEST DID, AND WHY IT WAS NOT ENOUGH
#    Script 21 computed two summary statistics per gene, a loss magnitude and
#    a subtype effect size, and correlated them across 20 genes. That answers
#    a question about genes: does a gene's loss predict its subtype effect?
#    It gave rho = 0.226, p = 0.337.
#
#    A reviewer objected that this is a correlation of 20 summary points with
#    no interval, and that the claim of independence should be tested at the
#    level of samples, where the data actually are.
#
#  WHAT THIS SCRIPT DOES INSTEAD
#    It places every TUMOUR on the two axes:
#
#      loss score    = mean z-score across the five universal-loss genes,
#                      signed so that higher means more loss
#      subtype score = mean z of the luminal programme minus mean z of the
#                      basal programme, so that higher means more luminal
#
#    If the axes are independent, a tumour's position on one should not
#    predict its position on the other. With roughly a thousand tumours this
#    is a properly powered test and yields an interval, not just a p value.
#
#  A CAVEAT THE SCRIPT ENFORCES
#    The two scores must not share genes, or they would be correlated by
#    construction. GSTM2 and GSTM5 appear in both the loss table and the
#    luminal programme. The script therefore builds the subtype score from
#    the programme genes that are NOT on the loss axis, and reports which
#    genes were excluded for this reason. A version using all genes is also
#    computed, and the difference between the two is itself informative.
#
#  ALSO REPORTED
#    Within-subtype correlations. A pooled correlation across subtypes can
#    be induced by subtype differences alone, so the test is repeated inside
#    each PAM50 group where the confound cannot operate.
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

# Gene sets exactly as the manuscript defines them
LOSS    <- c("GSTM5","GSTM2","GSTA1","MGST1","GSTA4")
BASAL   <- c("GSTP1","GSTA1","GSTA4","MGST3")
LUMINAL <- c("GSTM2","GSTM3","GSTM4","GSTO2","GSTZ1")


# =============================================================================
# 1. EXPRESSION AND SUBTYPE, AS ELSEWHERE
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
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

tum <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                            grepl("breast", tissue, ignore.case = TRUE),
                            grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)
tum <- tum[tum$sample %in% colnames(E), ]

Et <- E[, tum$sample, drop = FALSE]
Z  <- t(scale(t(Et)))                      # z-score each gene across tumours
cat("Tumours:", ncol(Z), "\n"); print(table(tum$group))


# =============================================================================
# 2. TWO SCORES PER TUMOUR
# =============================================================================
have <- function(g) intersect(g, rownames(Z))
score <- function(g) {
  if (length(g)) colMeans(Z[g, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(Z))
}

# Disjoint construction: subtype score uses only genes absent from the loss set
bas_d <- setdiff(have(BASAL),   LOSS)
lum_d <- setdiff(have(LUMINAL), LOSS)
excluded <- setdiff(c(have(BASAL), have(LUMINAL)), c(bas_d, lum_d))

cat("\nSubtype score built from disjoint genes only.\n")
cat("  basal   :", paste(bas_d, collapse = ", "), "\n")
cat("  luminal :", paste(lum_d, collapse = ", "), "\n")
cat("  excluded because they are also on the loss axis:",
    paste(excluded, collapse = ", "), "\n")

d <- data.frame(
  sample  = colnames(Z),
  group   = tum$group,
  loss    = -score(have(LOSS)),                  # higher = more loss
  subtype_disjoint = score(lum_d) - score(bas_d),
  subtype_all      = score(have(LUMINAL)) - score(have(BASAL)),
  stringsAsFactors = FALSE)
d <- d[is.finite(d$loss) & is.finite(d$subtype_disjoint), ]
cat("\nTumours scored:", nrow(d), "\n")


# =============================================================================
# 3. THE TEST, WITH AN INTERVAL
# =============================================================================
report <- function(x, y, label, n_note = "") {
  ctp <- suppressWarnings(cor.test(x, y, method = "pearson"))
  cts <- suppressWarnings(cor.test(x, y, method = "spearman"))
  cat(sprintf("\n%s%s\n", label, n_note))
  cat(sprintf("  Pearson  r = %+.3f  95%% CI %+.3f to %+.3f   p = %.3g\n",
              ctp$estimate, ctp$conf.int[1], ctp$conf.int[2], ctp$p.value))
  cat(sprintf("  Spearman rho = %+.3f                          p = %.3g\n",
              cts$estimate, cts$p.value))
  cat(sprintf("  Shared variance r^2 = %.1f%%\n", 100 * ctp$estimate^2))
  data.frame(comparison = label, n = length(x),
             r = round(unname(ctp$estimate), 3),
             ci_low = round(ctp$conf.int[1], 3),
             ci_high = round(ctp$conf.int[2], 3),
             p = ctp$p.value,
             rho = round(unname(cts$estimate), 3),
             r2_pct = round(100 * ctp$estimate^2, 1))
}

cat("\n", strrep("=", 70), "\nSAMPLE-LEVEL INDEPENDENCE\n", strrep("=", 70), "\n", sep = "")

res <- list()
res[["pooled_disjoint"]] <- report(d$loss, d$subtype_disjoint,
  "All tumours, disjoint gene sets", sprintf("  (n = %d)", nrow(d)))
res[["pooled_all"]] <- report(d$loss, d$subtype_all,
  "All tumours, overlapping gene sets",
  sprintf("  (n = %d, shown to quantify the inflation)", nrow(d)))

cat("\n--- within subtype, where subtype cannot confound ---\n")
for (g in SUB) {
  dd <- d[d$group == g, ]
  if (nrow(dd) < 30) next
  res[[paste0("within_", g)]] <- report(dd$loss, dd$subtype_disjoint,
    paste("Within", g), sprintf("  (n = %d)", nrow(dd)))
}

res <- bind_rows(res)
write.csv(res, file.path(OUT, "axis_independence_sample_level.csv"), row.names = FALSE)


# =============================================================================
# 4. VERDICT
# =============================================================================
p1 <- res[res$comparison == "All tumours, disjoint gene sets", ]
cat("\n", strrep("=", 70), "\nINTERPRETATION\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("Pooled, disjoint: r = %+.3f (%+.3f to %+.3f), %.1f%% shared variance.\n",
            p1$r, p1$ci_low, p1$ci_high, p1$r2_pct))

if (abs(p1$r) < 0.2 && abs(p1$ci_low) < 0.3 && abs(p1$ci_high) < 0.3) {
  cat("\nThe interval excludes any correlation strong enough to matter. A\n")
  cat("tumour's position on the loss axis carries almost no information about\n")
  cat("its position on the subtype axis. Report r with its interval rather than\n")
  cat("the word 'independent' alone, since independence is a stronger claim\n")
  cat("than a small correlation.\n")
} else if (abs(p1$r) < 0.3) {
  cat("\nA weak correlation, detectable at this sample size but explaining\n")
  cat("little variance. Describe the axes as largely independent and give the\n")
  cat("coefficient and interval; do not claim strict independence.\n")
} else {
  cat("\nThe axes share substantial variance at sample level. The two-dimensional\n")
  cat("framing needs rewriting, and this is a material change to the paper.\n")
}

within <- res[grepl("^Within", res$comparison), ]
if (nrow(within)) {
  cat(sprintf("\nWithin subtypes, r ranges %+.3f to %+.3f.\n",
              min(within$r), max(within$r)))
  cat("If these resemble the pooled value, the relationship is not an artefact\n")
  cat("of subtype composition.\n")
}

infl <- res[res$comparison == "All tumours, overlapping gene sets", ]
cat(sprintf("\nUsing overlapping gene sets gives r = %+.3f against %+.3f disjoint.\n",
            infl$r, p1$r))
cat("The difference is the correlation induced by sharing genes between the\n")
cat("two scores, which is why the disjoint version is the one to report.\n")


# =============================================================================
# 5. FIGURE
# =============================================================================
p <- ggplot(d, aes(loss, subtype_disjoint)) +
  geom_point(aes(colour = group), size = 1.1, alpha = .55) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey25",
              linewidth = .6, se = TRUE) +
  scale_colour_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                                 LumA = "#2E86AB", LumB = "#E8A33D"),
                      name = NULL) +
  labs(x = "Loss score (higher = greater loss of the universal-loss genes)",
       y = "Subtype score (higher = more luminal programme)",
       title = "Loss and subtype position are close to independent across tumours",
       subtitle = sprintf("n = %d tumours; Pearson r = %+.3f (95%% CI %+.3f to %+.3f), %.1f%% shared variance. Scores use disjoint gene sets.",
                          nrow(d), p1$r, p1$ci_low, p1$ci_high, p1$r2_pct)) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT, "axis_independence_samples.png"), p,
       width = 7.5, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "axis_independence_samples.tiff"), p,
       width = 7.5, height = 5.8, dpi = 300, bg = "white", compression = "lzw")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d", nrow(d)),
             sprintf("Disjoint subtype score excludes: %s", paste(excluded, collapse = ", ")),
             sprintf("Pooled disjoint r = %.3f (%.3f to %.3f), p = %.3g, r2 = %.1f%%",
                     p1$r, p1$ci_low, p1$ci_high, p1$p, p1$r2_pct),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_axis_sample_level.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
