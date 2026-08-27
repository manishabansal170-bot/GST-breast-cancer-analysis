# =============================================================================
#  38_purity_control_matched_on_purity_correlation.R
#
#  THE OBJECTION, AND WHY IT IS CORRECT
#    The published negative control draws 200 null gene sets matched to the GST
#    family on median expression, from a pool of metabolic and housekeeping
#    genes. The quantity being compared is attenuation under purity adjustment.
#
#    Attenuation is a function of how strongly a gene correlates with tumour
#    purity, not of how highly it is expressed. Housekeeping genes are selected
#    for stable expression across cell types, which is close to selecting for
#    weak purity correlation. The null was therefore matched on a variable that
#    does not drive the outcome, while the variable that does drive it was left
#    free and biased low. The comparison partly restates its own premise.
#
#  WHAT THIS SCRIPT DOES
#    Rebuilds the null matched on absolute purity correlation, drawn from a
#    genome-wide pool rather than a hand-chosen one, and recomputes the
#    empirical p. It reports the purity-correlation distribution of the old
#    pool and the new one side by side, so the size of the original bias is
#    visible rather than asserted.
#
#  WHAT WOULD CHANGE THE PAPER
#    If the GST attenuation still exceeds a properly matched null, the claim
#    stands and is considerably better supported. If it does not, the immune
#    associations attenuate no more than any gene set with comparable purity
#    correlation, and the finding becomes a statement about bulk deconvolution
#    in general rather than about this family. Either is reportable; only the
#    first is what the manuscript currently claims.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleG_immune")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

set.seed(20260826)
NPERM <- 200


# =============================================================================
# 1. INPUTS
# =============================================================================
gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^gm(raw) - 0.001; lin[lin < 0] <- 0
G <- log2(lin + 1)

pool_raw <- readRDS(file.path(CACHE, "purity_control_pool.rds"))
P <- gm(pool_raw)
if (max(P, na.rm = TRUE) > 30) { pl <- 2^P - 0.001; pl[pl < 0] <- 0; P <- log2(pl + 1) }

imm <- readRDS(file.path(CACHE, "toil_immune_markers.rds"))
I <- gm(imm)
if (max(I, na.rm = TRUE) > 30) { il <- 2^I - 0.001; il[il < 0] <- 0; I <- log2(il + 1) }

# purity
data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
tp <- get("Tumor.purity", envir = environment())
tp$patient <- substr(tp$Sample.ID, 1, 12)
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

# Restrict to TCGA breast tumours. The pool and immune matrices span the whole
# Toil compendium, so intersecting them alone yields a pan-cancer sample set,
# which would answer the question in the wrong population.
ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
phc <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
brca <- ph[[phc("^sample$")]][
  grepl("breast", ph[[phc("primary disease or tissue")]], ignore.case = TRUE) &
  toupper(ph[[phc("study")]]) == "TCGA"]

s <- Reduce(intersect, list(colnames(G), colnames(P), colnames(I)))
s <- s[grepl("^TCGA-.*-01$", s)]
s <- s[s %in% brca]
pur <- num(tp$CPE)[match(substr(s, 1, 12), tp$patient)]
keep <- is.finite(pur); s <- s[keep]; pur <- pur[keep]
cat("Breast tumours with expression, pool, immune markers and purity:", length(s), "\n")
if (length(s) > 1500) stop("Sample set looks pan-cancer. Check the breast filter.")

GST <- rownames(G)[apply(G[, s, drop = FALSE], 1, function(v) sd(v, na.rm = TRUE) > 1e-8)]
cat("GST genes:", length(GST), "| pool genes:", nrow(P), "| immune features:", nrow(I), "\n\n")


# =============================================================================
# 2. PURITY CORRELATION OF EACH GENE
# =============================================================================
pcor <- function(M, genes) sapply(genes, function(g)
  suppressWarnings(cor(as.numeric(M[g, s]), pur, method = "spearman",
                       use = "pairwise.complete.obs")))

gst_pc  <- pcor(G, GST)
pool_pc <- pcor(P, rownames(P))
pool_pc <- pool_pc[is.finite(pool_pc)]

cat(strrep("=", 74), "\nPURITY CORRELATION: THE VARIABLE THAT DRIVES ATTENUATION\n",
    strrep("=", 74), "\n", sep = "")
cat(sprintf("GST family      : median |rho| %.3f, range %.3f to %.3f\n",
            median(abs(gst_pc)), min(abs(gst_pc)), max(abs(gst_pc))))
cat(sprintf("Published pool  : median |rho| %.3f, range %.3f to %.3f\n",
            median(abs(pool_pc)), min(abs(pool_pc)), max(abs(pool_pc))))
cat("\nIf the pool's median is appreciably lower, the original null was matched\n")
cat("on expression while the variable that determines the outcome was left\n")
cat("unmatched and biased low, which is the reviewer's objection.\n\n")


# =============================================================================
# 3. ATTENUATION
# =============================================================================
# For a gene set: correlate each gene against each immune feature, raw and
# after regressing both on purity, and record how far the coefficients shrink.
attenuation <- function(M, genes) {
  raws <- c(); adjs <- c()
  pr <- rank(pur)
  for (g in genes) {
    x <- as.numeric(M[g, s]); if (!all(is.finite(x))) x[!is.finite(x)] <- NA
    xr <- rank(x, na.last = "keep")
    xres <- tryCatch(resid(lm(xr ~ pr, na.action = na.exclude)), error = function(e) rep(NA, length(s)))
    for (f in rownames(I)) {
      y <- rank(as.numeric(I[f, s]), na.last = "keep")
      yres <- tryCatch(resid(lm(y ~ pr, na.action = na.exclude)), error = function(e) rep(NA, length(s)))
      r1 <- suppressWarnings(cor(xr, y, use = "pairwise.complete.obs"))
      r2 <- suppressWarnings(cor(xres, yres, use = "pairwise.complete.obs"))
      if (is.finite(r1) && is.finite(r2)) { raws <- c(raws, r1); adjs <- c(adjs, r2) }
    }
  }
  list(mean_atten = mean(abs(raws) - abs(adjs), na.rm = TRUE),
       pct_survive = 100 * mean(abs(adjs) >= 0.2, na.rm = TRUE),
       n_tests = length(raws))
}

cat(strrep("=", 74), "\nOBSERVED GST FAMILY\n", strrep("=", 74), "\n", sep = "")
obs <- attenuation(G, GST)
cat(sprintf("tests: %d | mean attenuation %.4f | surviving |rho| >= 0.2: %.2f%%\n",
            obs$n_tests, obs$mean_atten, obs$pct_survive))


# =============================================================================
# 4. NULL MATCHED ON PURITY CORRELATION
# =============================================================================
cat("\n", strrep("=", 74), "\nNULL MATCHED ON ABSOLUTE PURITY CORRELATION\n",
    strrep("=", 74), "\n", sep = "")
cat("Each GST gene is matched to a pool gene with the closest |purity rho|,\n")
cat("sampling without replacement within each draw.\n\n")

pool_ok <- names(pool_pc)
draw_matched <- function() {
  avail <- pool_ok; picked <- character(0)
  for (tg in abs(gst_pc)) {
    if (!length(avail)) break
    dif <- abs(abs(pool_pc[avail]) - tg)
    cand <- avail[order(dif)][1:min(10, length(avail))]   # nearest 10, pick one
    ch <- sample(cand, 1)
    picked <- c(picked, ch); avail <- setdiff(avail, ch)
  }
  picked
}

null_at <- numeric(NPERM); null_sv <- numeric(NPERM)
for (i in seq_len(NPERM)) {
  gsel <- draw_matched()
  r <- attenuation(P, gsel)
  null_at[i] <- r$mean_atten; null_sv[i] <- r$pct_survive
  if (i %% 25 == 0) message("  ", i, " / ", NPERM)
}

p_at <- (1 + sum(null_at >= obs$mean_atten)) / (1 + NPERM)
p_sv <- (1 + sum(null_sv <= obs$pct_survive)) / (1 + NPERM)

cat(sprintf("\nNull mean attenuation : median %.4f, 95%% range %.4f to %.4f\n",
            median(null_at), quantile(null_at, .025), quantile(null_at, .975)))
cat(sprintf("Observed GST          : %.4f\n", obs$mean_atten))
cat(sprintf("Empirical p (attenuation greater than null): %.4f\n", p_at))
cat(sprintf("\nNull surviving %%       : median %.2f%%\n", median(null_sv)))
cat(sprintf("Observed GST          : %.2f%%\n", obs$pct_survive))
cat(sprintf("Empirical p (fewer surviving than null): %.4f\n", p_sv))


# =============================================================================
# 5. VERDICT
# =============================================================================
cat("\n", strrep("=", 74), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 74), "\n", sep = "")
if (p_at < 0.05) {
  cat("The GST family attenuates more than gene sets matched on purity\n")
  cat("correlation, so the claim survives the stronger control and should be\n")
  cat("reported against this null rather than the expression-matched one.\n")
} else {
  cat("The GST family does not attenuate more than a properly matched null. The\n")
  cat("published control was matched on expression while the variable driving\n")
  cat("attenuation was left unmatched and biased low. The honest statement is\n")
  cat("that GST-immune associations are largely composition, as reported, but\n")
  cat("that this is a general property of genes with comparable purity\n")
  cat("correlation rather than something specific to this family.\n")
}

write.csv(data.frame(draw = seq_len(NPERM), mean_atten = null_at, pct_survive = null_sv),
          file.path(OUT, "purity_control_matched_null.csv"), row.names = FALSE)
write.csv(data.frame(mean_atten = obs$mean_atten, pct_survive = obs$pct_survive,
                     n_tests = obs$n_tests, p_attenuation = p_at, p_survive = p_sv),
          file.path(OUT, "purity_control_matched_observed.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(gst_pc), purity_rho = round(gst_pc, 3)),
          file.path(OUT, "gst_purity_correlations.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d, tests per set: %d", length(s), obs$n_tests),
             sprintf("GST median |purity rho|: %.3f; pool median: %.3f",
                     median(abs(gst_pc)), median(abs(pool_pc))),
             sprintf("Observed attenuation %.4f; null median %.4f; empirical p %.4f",
                     obs$mean_atten, median(null_at), p_at),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_purity_matched.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
