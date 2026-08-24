# =============================================================================
#  26_coexpression_within_subtype.R
#
#  THE OBJECTION
#    Genes were assigned to a basal or a luminal programme BY their differential
#    expression across PAM50 subtypes. The co-expression analysis then
#    correlated those same genes against ESR1, GATA3, FOXA1, EGFR and the
#    keratins ACROSS ALL SUBTYPES POOLED, and the result was described as
#    reproducing the division "obtained independently".
#
#    It is not independent. Any gene higher in luminal tumours will correlate
#    with ESR1 in a mixed cohort, because ESR1 is itself higher in luminal
#    tumours. The correlation is arithmetically guaranteed by the assignment.
#
#  WHY THIS IS FIXABLE RATHER THAN FATAL
#    The manuscript already applies the right correction elsewhere. In the
#    microRNA analysis, recomputing within subtype collapsed miR-17-5p against
#    GSTM3 from -0.396 pooled to -0.078 within basal-like. The same operation
#    is required here.
#
#  WHAT SURVIVING WITHIN SUBTYPE WOULD MEAN
#    A correlation that persists inside a single subtype cannot be produced by
#    subtype composition, because there is no subtype variation left to produce
#    it. Such a correlation is evidence of a relationship between the two genes.
#    A correlation that vanishes was composition all along.
#
#  THE HONEST EXPECTATION
#    Most of these correlations will shrink substantially. The pooled values
#    around 0.4 to 0.5 are large for tumour co-expression, and the microRNA
#    precedent in this same manuscript suggests attenuation of that magnitude.
#    If they collapse, the section should be rewritten to say that the
#    co-expression pattern reflects subtype membership rather than a
#    transcriptional partnership, or removed. Either is preferable to leaving
#    a self-fulfilling result presented as corroboration.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleF_pathway_network")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")

# The gene-level partners named in the Results, plus the process panel members
# most relied upon. Extended only where the manuscript quotes a coefficient.
PARTNERS <- c("ESR1","GATA3","FOXA1","EGFR","KRT5","KRT17","GREB1","AR","XBP1",
              "MKI67","VIM","CDH1","ERBB2","PGR","TFF1","KRT14","SNAI2","ZEB1")
GST_FOCUS <- c("GSTP1","GSTA1","GSTA4","MGST3",          # basal programme
               "GSTM2","GSTM3","GSTM4","GSTO2","GSTZ1")  # luminal programme


# =============================================================================
# 1. EXPRESSION AND SUBTYPE
# =============================================================================
get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^get_matrix(raw) - 0.001; lin[lin < 0] <- 0
G <- log2(lin + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]], study = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")
clin <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                             grepl("breast", tissue, ignore.case = TRUE),
                             grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)


# =============================================================================
# 2. FETCH THE PARTNER GENES
# =============================================================================
pf <- file.path(CACHE, "coexpression_partners.rds")
if (file.exists(pf)) {
  Praw <- readRDS(pf); message("Loaded partner genes from cache.")
} else {
  HOST <- "https://toil.xenahubs.net"; EXPR <- "TcgaTargetGtex_rsem_gene_tpm"
  got <- list()
  for (g in PARTNERS) {
    r <- tryCatch(fetch_dense_values(host = HOST, dataset = EXPR, identifiers = g,
                                     use_probeMap = TRUE, check = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || !nrow(as.matrix(r))) { message("  failed: ", g); next }
    m <- as.matrix(r)
    got[[g]] <- if (nrow(m) > 1) colMeans(m, na.rm = TRUE) else m[1, ]
    Sys.sleep(0.15)
  }
  common <- Reduce(intersect, lapply(got, names))
  Praw <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(Praw) <- names(got); saveRDS(Praw, pf)
  message("Fetched ", nrow(Praw), " partner genes.")
}
pl <- 2^Praw - 0.001; pl[pl < 0] <- 0
P <- log2(pl + 1)

s <- Reduce(intersect, list(clin$sample, colnames(G), colnames(P)))
clin <- clin[match(s, clin$sample), ]
G <- G[intersect(GST_FOCUS, rownames(G)), s, drop = FALSE]
P <- P[, s, drop = FALSE]
cat("Tumours:", length(s), "\n"); print(table(clin$group))
cat("GST genes:", nrow(G), "| partners:", nrow(P), "\n")


# =============================================================================
# 3. POOLED VERSUS WITHIN-SUBTYPE
# =============================================================================
rho <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 25) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman")) }

res <- lapply(rownames(G), function(g) lapply(rownames(P), function(p) {
  x <- as.numeric(G[g, ]); y <- as.numeric(P[p, ])
  within <- sapply(SUB, function(k) { i <- clin$group == k; rho(x[i], y[i]) })
  data.frame(gst = g, partner = p,
             pooled = round(rho(x, y), 3),
             Basal = round(within["Basal"], 3), Her2 = round(within["Her2"], 3),
             LumA = round(within["LumA"], 3),  LumB = round(within["LumB"], 3),
             within_mean = round(mean(within, na.rm = TRUE), 3),
             max_abs_within = round(max(abs(within), na.rm = TRUE), 3))
}) %>% bind_rows()) %>% bind_rows()

res$attenuation <- round(abs(res$pooled) - abs(res$within_mean), 3)
res$survives <- abs(res$within_mean) >= 0.3

write.csv(res, file.path(OUT, "coexpression_within_subtype.csv"), row.names = FALSE)


# =============================================================================
# 4. THE PAIRS THE MANUSCRIPT QUOTES
# =============================================================================
QUOTED <- list(
  c("GSTP1","GATA3"), c("GSTP1","ESR1"), c("GSTP1","FOXA1"),
  c("GSTP1","EGFR"),  c("GSTP1","KRT17"), c("GSTP1","KRT5"),
  c("GSTM3","GREB1"), c("GSTM3","AR"),   c("GSTM3","XBP1"))

cat("\n", strrep("=", 78), "\nPAIRS QUOTED IN THE RESULTS\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("%-8s %-8s %8s %8s %8s %8s %8s %8s\n",
            "GST","partner","pooled","Basal","Her2","LumA","LumB","mean"))
for (q in QUOTED) {
  r <- res[res$gst == q[1] & res$partner == q[2], ]
  if (!nrow(r)) next
  cat(sprintf("%-8s %-8s %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f\n",
              r$gst, r$partner, r$pooled, r$Basal, r$Her2, r$LumA, r$LumB, r$within_mean))
}

quoted_df <- do.call(rbind, lapply(QUOTED, function(q)
  res[res$gst == q[1] & res$partner == q[2], ]))

cat(sprintf("\nMean absolute correlation: pooled %.3f, within subtype %.3f\n",
            mean(abs(quoted_df$pooled), na.rm = TRUE),
            mean(abs(quoted_df$within_mean), na.rm = TRUE)))
cat(sprintf("Of the %d quoted pairs, %d retain |rho| of at least 0.3 within subtype.\n",
            nrow(quoted_df), sum(quoted_df$survives, na.rm = TRUE)))


# =============================================================================
# 5. ALL PAIRS
# =============================================================================
cat("\n", strrep("=", 78), "\nALL PAIRS\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("Tested: %d\n", nrow(res)))
cat(sprintf("|rho| >= 0.3 pooled        : %d\n", sum(abs(res$pooled) >= 0.3, na.rm = TRUE)))
cat(sprintf("|rho| >= 0.3 within subtype: %d\n", sum(res$survives, na.rm = TRUE)))
cat(sprintf("Median attenuation         : %.3f\n", median(res$attenuation, na.rm = TRUE)))

surv <- res[which(res$survives), ]
if (nrow(surv)) {
  cat("\nPairs surviving within subtype, which cannot be produced by composition:\n")
  print(as.data.frame(surv[order(-abs(surv$within_mean)),
                           c("gst","partner","pooled","within_mean")]),
        row.names = FALSE, digits = 3)
} else {
  cat("\nNo pair retains |rho| of at least 0.3 within subtype.\n")
}


# =============================================================================
# 6. WHAT TO DO WITH THE SECTION
# =============================================================================
cat("\n", strrep("=", 78), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 78), "\n", sep = "")
frac <- mean(abs(quoted_df$within_mean) >= 0.3, na.rm = TRUE)
if (is.finite(frac) && frac >= 0.5) {
  cat("Most quoted correlations persist within subtype. The section can stand,\n")
  cat("but must report the within-subtype values alongside the pooled ones and\n")
  cat("drop the word 'independently', since the gene assignment and the pooled\n")
  cat("correlation share their source of variation.\n")
} else if (is.finite(frac) && frac > 0) {
  cat("A minority persist. Report both sets of values, restrict any claim to the\n")
  cat("pairs that survive, and state plainly that the remainder reflect subtype\n")
  cat("composition rather than a transcriptional relationship.\n")
} else {
  cat("None persists. The pooled correlations are a restatement of the subtype\n")
  cat("assignment. The section should be removed, or rewritten to report the\n")
  cat("attenuation itself as the finding, in the same way the microRNA analysis\n")
  cat("already does in this manuscript.\n")
}

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d", length(s)),
             sprintf("Quoted pairs: mean |rho| pooled %.3f, within subtype %.3f",
                     mean(abs(quoted_df$pooled), na.rm = TRUE),
                     mean(abs(quoted_df$within_mean), na.rm = TRUE)),
             sprintf("Pairs surviving within subtype: %d of %d",
                     sum(res$survives, na.rm = TRUE), nrow(res)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_coexpression_within.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
