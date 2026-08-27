# =============================================================================
#  39_cell_line_crosscheck.R
#
#  WHAT THIS ANSWERS
#    Every cohort in this manuscript is bulk tumour tissue, so every association
#    is open to the objection that it reflects what surrounds the cancer cells
#    rather than the cells themselves. The purity analysis addresses this
#    statistically. Cell lines address it by construction: they are pure
#    epithelium, with no stroma, no immune infiltrate and no adjacent normal.
#
#    If the GSTP1 subtype contrast holds across a panel of breast cancer cell
#    lines, composition cannot be producing it, because there is no
#    non-epithelial compartment for it to arise from.
#
#  WHAT THIS IS NOT
#    Cell lines are clonal, long-passaged and drift in culture. A pattern that
#    holds in them supports a tumour finding; it does not replace one, and the
#    manuscript reports this as a cross-check rather than as replication.
#
#    Note also that this is expression only. DepMap no longer distributes the
#    CCLE RRBS methylation data, so the methylation-expression relationship
#    cannot be tested in this system with the current release.
#
#  A SECOND USE
#    The same table pre-screens the cell lines named in the demethylation
#    protocol. A line labelled basal-like that does not in fact express GSTP1
#    cannot serve as an unmethylated control, and it is cheaper to discover
#    that here than after three weeks of culture.
#
#  DATA
#    DepMap Public, downloaded by hand from https://depmap.org/portal/download
#      Expression tab  -> OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
#      Model tab       -> Model.csv
#    Both are CC BY 4.0 and must be cited by release. The expression file is
#    around 300 MB, so it is not redeposited here; this script extracts what
#    the manuscript reports.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleH_celllines")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr", "readr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(lapply(pk, library, character.only = TRUE))

FOCUS <- c("GSTP1","GSTM2","GSTM3","GSTM5","GSTA1","GSTO2","MGST1")
PROTOCOL_LINES <- c("MCF7","T47D","MDAMB231","HCC1937")

f_expr <- file.path(CACHE, "ccle_expression.csv")
f_anno <- file.path(CACHE, "ccle_sample_info.csv")
for (f in c(f_expr, f_anno))
  if (!file.exists(f)) stop("Missing ", basename(f),
    ". Download from depmap.org/portal/download and save to the cache with this name.")


# =============================================================================
# 1. READ ONLY THE COLUMNS NEEDED
# =============================================================================
# The expression file is roughly 20,000 columns wide. Reading the header first
# and then selecting turns a twenty-minute read into a one-minute one.
hdr <- names(readr::read_csv(f_expr, n_max = 0, show_col_types = FALSE,
                             progress = FALSE))
sym <- sub(" .*$", "", hdr)                       # "GSTP1 (2950)" -> "GSTP1"
if (!("ModelID" %in% hdr))
  stop("ModelID not in the expression header. It begins: ",
       paste(head(hdr, 6), collapse = ", "))

want <- unique(c("ModelID", hdr[sym %in% FOCUS]))
cat("Reading", length(want), "of", length(hdr), "columns\n")
ex <- readr::read_csv(f_expr, col_select = all_of(want),
                      show_col_types = FALSE, progress = FALSE)
names(ex) <- sub(" .*$", "", names(ex))

# One cell line can contribute several sequencing profiles; keep one row each.
ex <- ex[!duplicated(ex$ModelID), ]

an <- readr::read_csv(f_anno, show_col_types = FALSE, progress = FALSE)
br <- an[grepl("breast", an$OncotreeLineage, ignore.case = TRUE), ]
cat("Breast lines in the annotation:", nrow(br), "\n")

ex <- ex[ex$ModelID %in% br$ModelID, ]
cat("Breast lines with expression:", nrow(ex), "\n\n")

i <- match(ex$ModelID, br$ModelID)
tab <- data.frame(model_id = ex$ModelID,
                  line     = br$StrippedCellLineName[i],
                  subtype_features = br$ModelSubtypeFeatures[i],
                  stringsAsFactors = FALSE)
for (g in intersect(FOCUS, names(ex))) tab[[g]] <- round(ex[[g]], 2)


# =============================================================================
# 2. THE FOUR LINES IN THE DEMETHYLATION PROTOCOL
# =============================================================================
cat(strrep("=", 74), "\nTHE FOUR LINES NAMED IN THE PROTOCOL\n",
    strrep("=", 74), "\n", sep = "")
cat("Values are log2(TPM + 1). The tumour data predict high GSTP1 in basal-like\n")
cat("lines and low GSTP1 in luminal lines.\n\n")

sel <- tab[toupper(gsub("[^A-Za-z0-9]", "", tab$line)) %in% PROTOCOL_LINES, ]
if (nrow(sel)) {
  print(sel[order(-sel$GSTP1), c("line","subtype_features","GSTP1","GSTM2","GSTM3","GSTM5")],
        row.names = FALSE)
} else cat("None of the four matched by name.\n")


# =============================================================================
# 3. DOES THE SUBTYPE CONTRAST HOLD ACROSS THE PANEL
# =============================================================================
cat("\n", strrep("=", 74), "\nGSTP1 ACROSS THE PANEL BY SUBTYPE\n",
    strrep("=", 74), "\n", sep = "")
cat("DepMap annotates breast lines with free-text receptor and subtype features\n")
cat("rather than a PAM50 call, so lines are grouped on that text. Lines carrying\n")
cat("both a basal and a receptor label, such as 'basal_A HER2+', are assigned to\n")
cat("the basal group, since the basal designation is the stronger signal; this\n")
cat("choice is stated because it is a judgement rather than a rule.\n\n")

tab$group <- NA_character_
isb <- grepl("basal|TNBC", tab$subtype_features, ignore.case = TRUE)
isl <- grepl("luminal|ER\\+|HER2\\+", tab$subtype_features, ignore.case = TRUE)
tab$group[isb] <- "Basal-like / TNBC"
tab$group[isl & !isb] <- "Luminal / HER2-positive"

t2 <- tab[!is.na(tab$group) & is.finite(tab$GSTP1), ]
cat("Classifiable lines:", nrow(t2), "of", nrow(tab), "\n\n")

print(t2 %>% group_by(group) %>%
        summarise(n = n(),
                  median_GSTP1 = round(median(GSTP1), 2),
                  Q1 = round(quantile(GSTP1, .25), 2),
                  Q3 = round(quantile(GSTP1, .75), 2), .groups = "drop"))

w <- wilcox.test(GSTP1 ~ group, data = t2)
a <- t2$GSTP1[t2$group == "Basal-like / TNBC"]
b <- t2$GSTP1[t2$group != "Basal-like / TNBC"]
r <- rank(c(a, b)); U <- sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2
delta <- 2 * U / (length(a) * length(b)) - 1
cat(sprintf("\nWilcoxon p = %.4g | Cliff's delta = %+.3f\n", w$p.value, delta))


# =============================================================================
# 4. WHAT THIS SUPPORTS
# =============================================================================
cat("\n", strrep("=", 74), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 74), "\n", sep = "")
if (w$p.value < 0.05 && delta > 0) {
  cat("The subtype contrast holds in pure epithelium. Because cell lines carry no\n")
  cat("stromal or immune compartment, cellular composition cannot generate this\n")
  cat("difference, which the bulk-tumour analyses cannot fully exclude on their\n")
  cat("own. Report as a cross-check in an independent system, not as replication.\n")
} else {
  cat("The contrast does not hold in cell lines. This does not refute the tumour\n")
  cat("finding, since cultured lines drift, but it should be reported rather than\n")
  cat("omitted, and it weakens the case that the difference is cell-intrinsic.\n")
}

write.csv(tab[order(-tab$GSTP1), ], file.path(OUT, "cell_line_gst_expression.csv"),
          row.names = FALSE)
writeLines(c(paste("Run:", Sys.time()),
             "Source: DepMap Public, expression and Model annotation, CC BY 4.0",
             sprintf("Breast lines with expression: %d", nrow(tab)),
             sprintf("Classifiable by subtype text: %d", nrow(t2)),
             sprintf("Median GSTP1 basal %.2f against luminal %.2f",
                     median(a), median(b)),
             sprintf("Wilcoxon p = %.4g, Cliff's delta = %+.3f", w$p.value, delta),
             "Note: DepMap no longer distributes CCLE RRBS methylation, so the",
             "methylation-expression relationship cannot be tested in this system.",
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_cell_lines.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
