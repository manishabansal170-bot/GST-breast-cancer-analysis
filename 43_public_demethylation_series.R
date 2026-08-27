# =============================================================================
#  43_public_demethylation_series.R
#
#  WHY THIS IS WORTH DOING
#    The manuscript's regulatory model is associative: promoter methylation and
#    expression covary, alternatives have been excluded, but nothing here shows
#    that removing the methylation restores the transcript. That step needs a
#    demethylation experiment, which this project has not run.
#
#    A reviewer pointed out that the experiment may not need running. Public
#    series exist in which breast cancer cell lines were treated with a DNA
#    methyltransferase inhibitor and profiled for expression. If GSTP1
#    re-expresses in a line where it is silenced, that is a causal step
#    obtained without new data.
#
#  THE SERIES USED
#    GSE74251 profiles MCF7 treated with decitabine (5-aza-2'-deoxycytidine)
#    at 100 nM, by RNA-seq, alongside untreated control. MCF7 is luminal and
#    carries low GSTP1 in the cell line panel (0.94 log2 TPM against 8.54 in
#    MDA-MB-231), so it is the right background: a gene that is already
#    expressed cannot be reactivated.
#
#  WHAT THIS CAN AND CANNOT SHOW
#    These series are small, often one or two samples per condition, and were
#    designed to answer other questions. A clear reactivation is suggestive
#    and worth reporting; an absence of change is uninformative rather than
#    negative, because a single unreplicated comparison has little power.
#
#    This script therefore reports the fold change for the GST family
#    alongside positive controls known to reactivate on demethylation, so that
#    the reader can judge whether the treatment worked at all before judging
#    what it did to the GST genes. Without those controls the result cannot be
#    interpreted in either direction.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleI_demethylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({ lapply(pk, library, character.only = TRUE); library(GEOquery) })

SERIES <- "GSE74251"
GST <- c("GSTP1","GSTM2","GSTM3","GSTM4","GSTM5","GSTA1","GSTA4","GSTO2","GSTZ1","GSTK1",
         "GSTO1","MGST1","MGST2","MGST3")
# Genes repeatedly shown to reactivate on DNMT inhibition; if these do not move,
# the treatment did not work in this series and nothing else can be read from it.
POSCON <- c("MAGEA1","MAGEA3","MAGEA6","GAGE1","XIST","SFRP1","RASSF1","CDKN2A","MLH1","TIMP3")


# =============================================================================
# 1. FETCH
# =============================================================================
cat(strrep("=", 76), "\n", SERIES, "\n", strrep("=", 76), "\n", sep = "")
g <- tryCatch(getGEO(SERIES, GSEMatrix = TRUE, getGPL = FALSE),
              error = function(e) { message("  fetch failed: ", conditionMessage(e)); NULL })
if (is.null(g)) stop("Could not retrieve ", SERIES,
                     ". Check the connection, or download the series matrix by hand.")

# The series matrix for this study carries no expression values; the counts are
# in a supplementary workbook with one sheet per cell line, already holding the
# three control and three treated replicates plus logFC and FDR.
es <- g[[1]]
pd <- Biobase::pData(es)
cat("Sample titles:\n")
print(data.frame(sample = rownames(pd), title = pd$title), row.names = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)
sf <- getGEOSuppFiles(SERIES, baseDir = CACHE, makeDirectory = TRUE)
xl <- rownames(sf)[grepl("\\.xlsx$", rownames(sf))][1]
if (is.na(xl)) stop("No supplementary workbook found for ", SERIES)
cat("\nSupplementary workbook:", basename(xl), "\n")
cat("Sheets:", paste(getSheetNames(xl), collapse = ", "), "\n\n")

# MCF7 is the breast line; the other sheet is a colon line and is not used.
d <- read.xlsx(xl, sheet = "MCF7")
names(d)[1] <- "gene"
cat("Genes in the MCF7 table:", nrow(d), "\n")
ctrl <- grep("Cntrl", names(d), value = TRUE)
dac  <- grep("DAC",   names(d), value = TRUE)
cat("Control replicates:", length(ctrl), "| treated replicates:", length(dac), "\n\n")
cat("Samples:", ncol(ex), " Features:", nrow(ex), "\n\n")

cat("Sample titles:\n")
print(data.frame(sample = rownames(pd), title = pd$title), row.names = FALSE)

if (nrow(ex) == 0) {
  cat("\nThe series matrix carries no expression values. Many RNA-seq series on\n")
  cat("GEO hold counts only in supplementary files rather than in the matrix.\n")
  cat("Retrieve the supplementary counts table with getGEOSuppFiles('", SERIES,
      "') and rerun from section 2.\n", sep = "")
  stop("No expression matrix in this series.")
}


# =============================================================================
# 2. THE COMPARISON
# =============================================================================
report <- function(genes, label) {
  cat("\n", strrep("-", 76), "\n", label, "\n", strrep("-", 76), "\n", sep = "")
  x <- d[d$gene %in% genes, c("gene", ctrl, dac, "logFC", "FDR")]
  if (!nrow(x)) { cat("  none of these genes are in the table\n"); return(NULL) }
  x$ctrl_mean <- round(rowMeans(x[, ctrl], na.rm = TRUE), 1)
  x$dac_mean  <- round(rowMeans(x[, dac],  na.rm = TRUE), 1)
  print(x[order(-x$logFC), c("gene","ctrl_mean","dac_mean","logFC","FDR")],
        row.names = FALSE, digits = 3)
  x$set <- label
  x[, c("gene","ctrl_mean","dac_mean","logFC","FDR","set")]
}

pc <- report(POSCON, "POSITIVE CONTROLS: methylation-silenced genes known to reactivate")
gs <- report(GST,    "GST FAMILY")


# =============================================================================
# 3. VERDICT
# =============================================================================
cat("\n", strrep("=", 76), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 76), "\n", sep = "")

up_pc <- sum(pc$logFC > 1 & pc$FDR < 0.05, na.rm = TRUE)
cat(sprintf("Positive controls rising more than twofold at FDR < 0.05: %d of %d\n",
            up_pc, nrow(pc)))

if (up_pc == 0) {
  cat("\nThe treatment did not reactivate the genes it should have. Nothing can be\n")
  cat("concluded about the GST family from this series in either direction.\n")
} else {
  gp <- gs[gs$gene == "GSTP1", ]
  cat(sprintf("\nGSTP1: control %.1f, treated %.1f, log2 FC %+.2f, FDR %.2g\n",
              gp$ctrl_mean, gp$dac_mean, gp$logFC, gp$FDR))
  risen <- gs$gene[gs$logFC > 1 & gs$FDR < 0.05]
  flat  <- gs$gene[abs(gs$logFC) < 1]
  cat("Rising more than twofold:", paste(risen, collapse = ", "), "\n")
  cat("Essentially unchanged   :", paste(flat, collapse = ", "), "\n")
  cat("\nThe genes that reactivate should be those the manuscript identifies as\n")
  cat("methylation-silenced and that sit near zero untreated. If abundantly\n")
  cat("expressed family members do not rise, the effect is specific rather than a\n")
  cat("general transcriptional response, and that specificity is the argument.\n")
  cat("\nThis is one cell line in one published series. It is a causal step, not a\n")
  cat("demonstration that the mechanism operates in tumours.\n")
}

res <- bind_rows(pc, gs)
write.csv(res, file.path(OUT, "demethylation_GSE74251_MCF7.csv"), row.names = FALSE)
writeLines(c(paste("Run:", Sys.time()), paste("Series:", SERIES),
             sprintf("MCF7 control replicates: %d, treated: %d", length(ctrl), length(dac)),
             sprintf("Positive controls reactivating: %d of %d", up_pc, nrow(pc)),
             sprintf("GSTP1 log2 FC: %+.2f (FDR %.2g)",
                     gs$logFC[gs$gene == "GSTP1"], gs$FDR[gs$gene == "GSTP1"]),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_demethylation.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
