# =============================================================================
#  Table 3, paired: methylation and expression in the same tumours
#
#  As published, Table 3 places a methylation column from 737 array tumours
#  beside an expression column from 1,041 tumours on a different annotation.
#  The caption says plainly that the two should not be read as paired, which is
#  candid but leaves the paper's most-quoted table a juxtaposition rather than
#  a measurement.
#
#  780 tumours carry both assays. Computing the expression column in those
#  same tumours, on the same annotation as the methylation work, makes the
#  table paired. Both annotations are reported here so the difference is
#  visible rather than asserted.
# =============================================================================
WORKDIR <- "~/GST_BRCA"; CACHE <- file.path(WORKDIR,"cache")
OUT <- file.path(WORKDIR,"objective1","moduleC_methylation")
setwd(WORKDIR)
suppressPackageStartupMessages({library(dplyr); library(TCGAbiolinks)})
SUB <- c("Basal","Her2","LumA","LumB")

obj <- readRDS(file.path(CACHE,"tcga_methylation_gst.rds"))
met <- obj$met; exp_m <- obj$exp

gm <- function(x){ if(is.matrix(x)) return(x)
  if(is.list(x)){ if(!is.null(x$mat)) return(x$mat)
    i<-which(vapply(x,is.matrix,logical(1))); if(length(i)) return(x[[i[1]]])}
  stop("no matrix")}
raw <- readRDS(file.path(CACHE,"toil_gst.rds"))
LIN <- 2^gm(raw)-0.001; LIN[LIN<0] <- 0        # TPM on the harmonised annotation

id <- function(v) substr(v,1,15)
colnames(met) <- id(colnames(met)); colnames(exp_m) <- id(colnames(exp_m))
colnames(LIN) <- id(colnames(LIN))

st <- TCGAquery_subtype(tumor="brca")
s <- intersect(colnames(met), colnames(exp_m))
s <- s[substr(s,13,15)=="-01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(s,1,12), st$patient)]
ok <- grp %in% SUB
s <- s[ok]; grp <- factor(grp[ok], levels=SUB)
cat("Tumours with methylation, expression and a PAM50 call:", length(s), "\n")
print(table(grp))

sh <- intersect(s, colnames(LIN))
cat("Of those, also present in the harmonised GRCh38 matrix:", length(sh), "\n\n")

b  <- as.numeric(met["GSTP1", s])
tab <- data.frame(
  subtype = c("Basal-like","HER2-enriched","Luminal A","Luminal B"),
  pam50   = SUB,
  n       = as.integer(table(grp)),
  median_beta      = round(tapply(b, grp, median, na.rm=TRUE), 3),
  pct_methylated   = round(tapply(b, grp, function(x) 100*mean(x>0.3, na.rm=TRUE)), 1),
  pct_unmethylated = round(tapply(b, grp, function(x) 100*mean(x<0.1, na.rm=TRUE)), 1),
  pct_intermediate = round(tapply(b, grp, function(x) 100*mean(x>=0.1 & x<=0.3, na.rm=TRUE)), 1),
  stringsAsFactors = FALSE)

# expression in the SAME tumours, harmonised annotation, TPM
gs <- grp[match(sh, s)]
tpm <- as.numeric(LIN["GSTP1", sh])
tab$expr_TPM_paired <- round(tapply(tpm, gs, median, na.rm=TRUE)[SUB], 1)
tab$n_expr <- as.integer(table(gs)[SUB])

# and on the methylation matrix's own annotation, for comparison
em <- as.numeric(exp_m["GSTP1", s])
tab$expr_methmatrix_log2 <- round(tapply(em, grp, median, na.rm=TRUE)[SUB], 2)

print(tab[, c("subtype","n","median_beta","pct_methylated","pct_unmethylated",
              "pct_intermediate","n_expr","expr_TPM_paired","expr_methmatrix_log2")],
      row.names=FALSE)

cat("\nPreviously published expression column, from 1,041 tumours:\n")
cat("  Basal 522.2, HER2 159.2, LumA 141.3, LumB 72.2 TPM\n")
cat("\nPaired column above is computed in the tumours that also have methylation,\n")
cat("on the same GRCh38 annotation, so the two columns of Table 3 now describe\n")
cat("the same patients.\n")

write.csv(tab, file.path(OUT,"table3_paired.csv"), row.names=FALSE)
cat("\nWritten to", normalizePath(OUT), "\n")
