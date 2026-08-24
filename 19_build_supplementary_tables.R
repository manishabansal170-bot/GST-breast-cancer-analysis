# =============================================================================
#  19_build_supplementary_tables.R
#
#  Assembles the supplementary tables into a single Excel workbook, one sheet
#  per table, with a contents sheet listing what each contains and which
#  script produced it.
#
#  WHAT THIS DOES NOT DO
#    The six per-sample matrices (S1a to S1f) are left as separate CSV files.
#    They run to thousands of columns and belong in the data deposit, not in
#    a spreadsheet a reviewer opens. The workbook covers S2 to S10, which are
#    the tables a reader actually consults while reading the paper.
#
#  BEFORE COMMITTING TO THE OUTPUT
#    The script reports any source file that is missing, empty, or duplicated
#    between the working folder and zenodo_deposit. Read that report. Several
#    files in the survival module were zero bytes at last check, and a table
#    promised in the manuscript but empty in the deposit is worse than one
#    that was never promised.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
OUT     <- file.path(WORKDIR, "supplementary_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
suppressPackageStartupMessages(library(openxlsx))

# -----------------------------------------------------------------------------
# 1. TABLE DEFINITIONS
#    Each entry: sheet name, caption, and one or more source files. Where two
#    files are listed they are stacked with a 'source' column identifying which
#    is which, because splitting them across sheets would fragment one table.
# -----------------------------------------------------------------------------
TABLES <- list(
  S2 = list(
    caption = "Per-gene summary and detection status by sample group. All 20 genes returning any signal in the GRCh38 requantification; GSTT1 and GSTT4 are absent because they are not quantifiable on that assembly.",
    files = c("figures/pergene_pub_toil/pergene_summary.csv",
              "figures/pergene_pub_toil/detection_bygroup_toil.csv"),
    labels = c("summary", "detection_by_group")),

  S3 = list(
    caption = "Cross-cohort concordance of per-gene subtype ordering. 19 genes shared across all three cohorts.",
    files = c("figures/scanb/three_cohort_concordance.csv"),
    labels = c("three_cohort")),

  S4 = list(
    caption = "Methylation-expression correlations, purity adjusted. 18 genes carrying probes on the methylation array; this is two fewer than the expression panel and two more than the 16 quantifiable by RNA-seq.",
    files = c("objective1/moduleC_methylation/methylation_expression_purity_adjusted.csv"),
    labels = c("purity_adjusted")),

  S5 = list(
    caption = "Promoter methylation by intrinsic subtype. Same 18 array-covered genes as S4.",
    files = c("objective1/moduleC_methylation/methylation_by_subtype.csv",
              "objective1/moduleC_methylation/GSTP1_methylation_by_subtype.csv"),
    labels = c("all_genes", "GSTP1_detail")),

  S6 = list(
    caption = "Copy number and mutation frequency, computed against the number of samples profiled rather than the number returned. Copy number covers 21 loci; mutation covers only the 12 genes with any observed variant.",
    files = c("objective1/moduleBG_alterations_protein/copy_number_frequency.csv",
              "objective1/moduleBG_alterations_protein/mutation_frequency.csv"),
    labels = c("copy_number", "mutation")),

  S7 = list(
    caption = "Transcript-protein concordance in the 74 tumours with matched proteomics. 17 proteins measured.",
    files = c("supplementary_tables/Supplementary_Table_S7_mrna_protein_concordance.csv",
              "objective1/moduleBG_alterations_protein/mrna_protein_concordance.csv"),
    labels = c("concordance", "concordance_alt")),

  S8 = list(
    caption = "Univariate Cox results and incremental C-index against the clinical model. TCGA n = 1,018; METABRIC n = 1,608, both restricted to the four PAM50 subtypes.",
    files = c("objective1/moduleD_survival/univariate_cox_tcga.csv",
              "objective1/moduleD_survival/univariate_cox_metabric.csv",
              "objective1/moduleD_survival/incremental_cindex_metabric.csv"),
    labels = c("cox_TCGA", "cox_METABRIC", "incremental_cindex")),

  S9 = list(
    caption = "Variance partition, methylation versus copy number, in the 770 tumours with all three assays. Same 18 array-covered genes as S4 and S5.",
    files = c("supplementary_tables/Supplementary_Table_S9_variance_partition.csv",
              "objective1/moduleBG_alterations_protein/cna_vs_methylation_variance.csv"),
    labels = c("variance_partition", "variance_partition_alt")),

  S10 = list(
    caption = "Purity negative control. Null distribution across 200 expression-matched random gene sets, with the observed GST family result for comparison.",
    files = c("objective1/moduleG_immune/purity_negative_control_null.csv",
              "objective1/moduleG_immune/purity_negative_control_observed.csv"),
    labels = c("null_distribution", "observed_GST"))
)

# -----------------------------------------------------------------------------
# 2. AUDIT BEFORE BUILDING
#    Report anything missing or empty rather than silently writing a blank
#    sheet. A promised-but-empty table is the kind of thing an editorial
#    office returns a submission over.
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nSOURCE FILE AUDIT\n", strrep("=", 70), "\n", sep = "")

audit <- do.call(rbind, lapply(names(TABLES), function(id) {
  do.call(rbind, lapply(seq_along(TABLES[[id]]$files), function(i) {
    f <- TABLES[[id]]$files[i]
    ex <- file.exists(f)
    sz <- if (ex) file.size(f) else NA_real_
    nr <- NA_integer_
    if (ex && !is.na(sz) && sz > 0) {
      nr <- tryCatch(nrow(read.csv(f, check.names = FALSE)), error = function(e) NA_integer_)
    }
    data.frame(table = id, label = TABLES[[id]]$labels[i], file = f,
               exists = ex, bytes = sz, rows = nr,
               status = if (!ex) "MISSING" else if (is.na(sz) || sz == 0) "EMPTY"
                        else if (is.na(nr) || nr == 0) "NO ROWS" else "ok",
               stringsAsFactors = FALSE)
  }))
}))

print(audit[, c("table","label","status","rows","bytes")], row.names = FALSE)

bad <- audit[audit$status != "ok", ]
if (nrow(bad)) {
  cat("\n--- PROBLEMS ---\n")
  for (i in seq_len(nrow(bad)))
    cat(sprintf("  %-4s %-22s %s\n     %s\n",
                bad$table[i], bad$label[i], bad$status[i], bad$file[i]))
  cat("\nSheets for these will be omitted. If a table is promised in the\n")
  cat("manuscript, either regenerate the source or remove the promise.\n")
}

# -----------------------------------------------------------------------------
# 3. CHECK FOR DUPLICATE COPIES
#    Two folders hold apparently identical supplementary files. If they differ,
#    that matters; if identical, the zenodo_deposit copy is redundant.
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nDUPLICATE CHECK\n", strrep("=", 70), "\n", sep = "")
a_dir <- "supplementary_tables"
b_dir <- "zenodo_deposit/supplementary_tables"
if (dir.exists(a_dir) && dir.exists(b_dir)) {
  common <- intersect(list.files(a_dir, pattern = "\\.csv$"),
                      list.files(b_dir, pattern = "\\.csv$"))
  if (length(common)) {
    for (f in common) {
      m1 <- tools::md5sum(file.path(a_dir, f))
      m2 <- tools::md5sum(file.path(b_dir, f))
      cat(sprintf("  %-58s %s\n", substr(f, 1, 58),
                  if (identical(unname(m1), unname(m2))) "identical" else "DIFFER"))
    }
    cat("\nAny marked DIFFER must be reconciled before deposit; the two copies\n")
    cat("are not interchangeable and only one can be correct.\n")
  }
} else {
  cat("  One or both folders absent; skipped.\n")
}

# -----------------------------------------------------------------------------
# 4. BUILD THE WORKBOOK
# -----------------------------------------------------------------------------
wb <- createWorkbook()

hdr <- createStyle(textDecoration = "bold", fgFill = "#1F3864",
                   fontColour = "#FFFFFF", halign = "left", border = "bottom")
cap <- createStyle(textDecoration = "italic", fontColour = "#555555")

# Contents sheet first
addWorksheet(wb, "Contents")
contents <- data.frame(
  Table = character(0), Contents = character(0),
  Rows = integer(0), Note = character(0), stringsAsFactors = FALSE)

built <- character(0)
for (id in names(TABLES)) {
  spec <- TABLES[[id]]
  parts <- list()
  for (i in seq_along(spec$files)) {
    f <- spec$files[i]
    if (!file.exists(f) || file.size(f) == 0) next
    df <- tryCatch(read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    if (length(spec$files) > 1) df <- cbind(source = spec$labels[i], df)
    parts[[length(parts) + 1]] <- df
  }
  if (!length(parts)) {
    contents <- rbind(contents, data.frame(Table = id, Contents = spec$caption,
      Rows = NA_integer_, Note = "source missing or empty", stringsAsFactors = FALSE))
    next
  }
  # stack only if columns align; otherwise keep the first and note it
  dat <- parts[[1]]
  if (length(parts) > 1) {
    same <- all(vapply(parts, function(x) identical(names(x), names(parts[[1]])), logical(1)))
    if (same) dat <- do.call(rbind, parts)
  }
  addWorksheet(wb, id)
  writeData(wb, id, spec$caption, startRow = 1)
  addStyle(wb, id, cap, rows = 1, cols = 1)
  writeData(wb, id, dat, startRow = 3, headerStyle = hdr)
  freezePane(wb, id, firstActiveRow = 4)
  setColWidths(wb, id, cols = seq_len(ncol(dat)), widths = "auto")
  built <- c(built, id)
  contents <- rbind(contents, data.frame(Table = id, Contents = spec$caption,
    Rows = nrow(dat),
    Note = if (length(parts) > 1) "multiple sources, see 'source' column" else "",
    stringsAsFactors = FALSE))
}

writeData(wb, "Contents",
  "Supplementary tables. Per-sample matrices (S1a to S1f) are provided as separate CSV files.",
  startRow = 1)
addStyle(wb, "Contents", cap, rows = 1, cols = 1)
writeData(wb, "Contents", contents, startRow = 3, headerStyle = hdr)
setColWidths(wb, "Contents", cols = 1:4, widths = c(10, 62, 10, 40))

xl <- file.path(OUT, "Supplementary_Tables_S2_to_S10.xlsx")
saveWorkbook(wb, xl, overwrite = TRUE)

cat("\n", strrep("=", 70), "\nDONE\n", strrep("=", 70), "\n", sep = "")
cat("Sheets written:", paste(built, collapse = ", "), "\n")
cat("Workbook:", normalizePath(xl), "\n")
cat("\nS1a to S1f remain as separate CSVs in", normalizePath(OUT), "\n")
cat("Check the audit above before submitting.\n")
