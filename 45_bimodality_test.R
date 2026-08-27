# =============================================================================
#  45_bimodality_test.R
#
#  The Figure 5 legend states that GSTP1 methylation "is bimodal in every
#  subtype". A reviewer objects that the luminal A panel shows a large mass
#  near zero with a long decaying tail and no visible second mode, and that
#  Table 3 supports that reading: 26.7% methylated, 56.9% unmethylated,
#  median 0.055.
#
#  Bimodality is testable rather than a matter of impression. Hartigan's dip
#  statistic tests the null of unimodality; a small p rejects it. This runs the
#  test per subtype and reports the two-component mixture fit alongside, so
#  that the legend can be written from a result rather than from the shape of
#  a violin plot.
# =============================================================================

WORKDIR <- "~/GST_BRCA"; CACHE <- file.path(WORKDIR,"cache")
OUT <- file.path(WORKDIR,"objective1","moduleC_methylation")
setwd(WORKDIR)
for (p in c("diptest","mclust"))
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
suppressPackageStartupMessages({
  library(diptest); library(mclust); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")
obj <- readRDS(file.path(CACHE,"tcga_methylation_gst.rds"))
met <- obj$met; colnames(met) <- substr(colnames(met),1,15)
st  <- TCGAquery_subtype(tumor="brca")
s   <- colnames(met)[substr(colnames(met),13,15)=="-01"]
g   <- st$BRCA_Subtype_PAM50[match(substr(s,1,12), st$patient)]
ok  <- g %in% SUB; s <- s[ok]; g <- factor(g[ok], levels=SUB)
b   <- as.numeric(met["GSTP1", s])

cat(strrep("=",78), "\nHARTIGAN DIP TEST FOR UNIMODALITY, GSTP1 BETA BY SUBTYPE\n",
    strrep("=",78), "\n", sep="")
cat("The null is that the distribution is unimodal. A small p rejects it and\n")
cat("supports bimodality; a large p means unimodality cannot be rejected.\n\n")

res <- do.call(rbind, lapply(SUB, function(k) {
  v <- b[g==k]; v <- v[is.finite(v)]
  if (length(v) < 20) return(NULL)
  dt <- dip.test(v)
  mc <- tryCatch(Mclust(v, G=1:2, verbose=FALSE), error=function(e) NULL)
  data.frame(subtype=k, n=length(v),
             median=round(median(v),3),
             dip=round(dt$statistic,4), dip_p=dt$p.value,
             best_G = if (is.null(mc)) NA_integer_ else mc$G,
             verdict = if (dt$p.value < 0.05) "bimodal" else "unimodal not rejected",
             stringsAsFactors=FALSE)
}))
print(res, row.names=FALSE, digits=3)

cat("\nMixture fit: best_G is the number of components favoured by BIC. Where the\n")
cat("dip test does not reject unimodality and BIC prefers one component, the\n")
cat("legend should not describe that subtype as bimodal.\n")

cat("\n", strrep("=",78), "\nCONSEQUENCE FOR THE LEGEND\n", strrep("=",78), "\n", sep="")
bi <- res$subtype[res$dip_p < 0.05]
un <- res$subtype[res$dip_p >= 0.05]
cat("Bimodal by the dip test:", if(length(bi)) paste(bi, collapse=", ") else "none", "\n")
cat("Unimodality not rejected:", if(length(un)) paste(un, collapse=", ") else "none", "\n")
if (length(un)) {
  cat("\nThe legend claim of bimodality 'in every subtype' is not supported. Name\n")
  cat("the subtypes in which it holds and describe the others as they are.\n")
} else cat("\nThe legend stands; cite the dip statistics.\n")

write.csv(res, file.path(OUT,"gstp1_bimodality_diptest.csv"), row.names=FALSE)
cat("\nWritten to", normalizePath(OUT), "\n")
