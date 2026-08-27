# =============================================================================
#  22_assemble_submission_folder.R
#
#  Collects everything a journal needs into one clean folder, renamed so that
#  a stranger can tell what each file is without opening it, and reports
#  anything the manuscript promises but the folder lacks.
#
#  WHY RENAME
#    Your working filenames encode where a file came from in your pipeline,
#    which is useful to you and meaningless to an editorial office. A file
#    called Panel_SCANB.png is Figure 2; a file called axis_independence.png
#    is Supplementary Figure S6. Naming them as the manuscript refers to them
#    removes an entire category of submission error.
#
#  THE ORIGINALS ARE NOT TOUCHED
#    Everything is copied, never moved. Rerun as often as you like.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
SUB     <- file.path(WORKDIR, "submission")
setwd(WORKDIR)

for (dd in c("main_figures", "supplementary_figures", "supplementary_tables"))
  dir.create(file.path(SUB, dd), recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# What the manuscript promises. Source on the left, submission name on the right.
# -----------------------------------------------------------------------------
MAIN <- c(
  "figures/Fig_panel_toil.tiff"                      = "Figure_1.tiff",
  "figures/scanb/Panel_SCANB.tiff"                   = "Figure_2.tiff",
  "figures/Panel_SCANB.tiff"                         = "Figure_2.tiff",
  "figures/metabric/Panel_METABRIC.tiff"             = "Figure_3.tiff",
  "figures/Panel_METABRIC.tiff"                      = "Figure_3.tiff",
  "figures/Panel_METABRIC.png"                       = "Figure_3.png"
)

# Renumbered after the former S1, a group-median heatmap, was withdrawn: no
# script in the deposit produced it, so a figure in the paper had no code
# behind it. METABRIC is now Figure 3 in the main text, not a supplementary
# panel, so it is listed under MAIN.
SUPP_FIG <- c(
  "figures/withintumour_toil/Panel_withintumour_toil.png" = "Supplementary_Figure_S1a.png",
  "figures/Panel_withintumour_toil.png"              = "Supplementary_Figure_S1a.png",
  "figures/metabric/Panel_METABRIC.png"              = "Supplementary_Figure_S1b.png",
  "figures/Panel_METABRIC.png"                       = "Supplementary_Figure_S1b.png",
  "figures/methylation_vs_expression.png"            = "Supplementary_Figure_S2.png",
  "objective1/moduleC_methylation/methylation_vs_expression.png" = "Supplementary_Figure_S2.png",
  "objective1/moduleD_survival/KM_riskscore_tcga.png" = "Supplementary_Figure_S3a.png",
  "objective1/moduleD_survival/KM_riskscore_metabric.png" = "Supplementary_Figure_S3b.png",
  "objective1/moduleE_drug_response/GSE25066_pCR_boxplots.png" = "Supplementary_Figure_S4a.png",
  "objective1/moduleE_drug_response/GSE25066_GSTP1_by_ER.png"  = "Supplementary_Figure_S4b.png",
  "figures/axis_independence.png"                    = "Supplementary_Figure_S5.png"
)

SUPP_TAB <- c(
  "supplementary_tables/Supplementary_Table_S1a_expression_TCGA_GTEx.csv" = "Supplementary_Table_S1a_expression_TCGA_GTEx.csv",
  "supplementary_tables/Supplementary_Table_S1b_expression_METABRIC.csv"  = "Supplementary_Table_S1b_expression_METABRIC.csv",
  "supplementary_tables/Supplementary_Table_S1c_expression_SCANB.csv"     = "Supplementary_Table_S1c_expression_SCANB.csv",
  "supplementary_tables/Supplementary_Table_S1d_methylation_TCGA.csv"     = "Supplementary_Table_S1d_methylation_TCGA.csv",
  "supplementary_tables/Supplementary_Table_S1e_copynumber_TCGA.csv"      = "Supplementary_Table_S1e_copynumber_TCGA.csv",
  "supplementary_tables/Supplementary_Table_S1f_sample_groups.csv"        = "Supplementary_Table_S1f_sample_groups.csv",
  "supplementary_tables/Supplementary_Tables_S2_to_S10.xlsx"              = "Supplementary_Tables_S2_to_S10.xlsx"
)


copy_set <- function(map, dest, label) {
  cat("\n", strrep("-", 70), "\n", label, "\n", strrep("-", 70), "\n", sep = "")
  done <- character(0); missing <- character(0)
  for (i in seq_along(map)) {
    src <- names(map)[i]; dst <- unname(map[i])
    if (dst %in% done) next                     # an earlier alternative worked
    if (file.exists(src)) {
      ok <- file.copy(src, file.path(SUB, dest, dst), overwrite = TRUE)
      if (ok) { done <- c(done, dst)
        cat(sprintf("  ok       %-42s <- %s\n", dst, src)) }
      else cat(sprintf("  FAILED   %-42s <- %s\n", dst, src))
    }
  }
  wanted <- unique(unname(map))
  missing <- setdiff(wanted, done)
  if (length(missing)) {
    cat("\n  NOT FOUND:\n")
    for (m in missing) {
      alts <- names(map)[unname(map) == m]
      cat(sprintf("    %s\n", m))
      for (a in alts) cat(sprintf("       looked in: %s\n", a))
    }
  }
  invisible(list(done = done, missing = missing))
}

r1 <- copy_set(MAIN,     "main_figures",          "MAIN FIGURES")
r2 <- copy_set(SUPP_FIG, "supplementary_figures", "SUPPLEMENTARY FIGURES")
r3 <- copy_set(SUPP_TAB, "supplementary_tables",  "SUPPLEMENTARY TABLES")


# -----------------------------------------------------------------------------
# Search for anything missing, in case it lives somewhere unexpected
# -----------------------------------------------------------------------------
allmiss <- c(r1$missing, r2$missing, r3$missing)
if (length(allmiss)) {
  cat("\n", strrep("=", 70), "\nSEARCHING FOR MISSING FILES\n", strrep("=", 70), "\n", sep = "")
  for (m in allmiss) {
    stem <- sub("^Supplementary_Figure_S[0-9]+[ab]?$", "", tools::file_path_sans_ext(m))
    # search on a distinctive fragment of the expected original name
    frag <- switch(m,
      "Supplementary_Figure_S1a.png"  = "withintumour",
      "Supplementary_Figure_S1b.png"  = "METABRIC",
      "Supplementary_Figure_S2.png"   = "methylation_vs",
      "Supplementary_Figure_S4a.png"  = "pCR_boxplot",
      "Supplementary_Figure_S4b.png"  = "GSTP1_by_ER",
      "Figure_3.tiff"                 = "Panel_METABRIC",
      "Figure_1.tiff"                 = "Fig_panel_toil",
      "Figure_2.tiff"                 = "Panel_SCANB",
      NULL)
    if (is.null(frag)) next
    hits <- list.files(WORKDIR, pattern = frag, recursive = TRUE,
                       full.names = TRUE, ignore.case = TRUE)
    cat(sprintf("\n%s  (looking for '%s')\n", m, frag))
    if (length(hits)) for (h in head(hits, 6)) cat("   found:", sub(".*GST_BRCA/", "", h), "\n")
    else cat("   nothing matching found anywhere under GST_BRCA\n")
  }
}


# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nFOLDER CONTENTS\n", strrep("=", 70), "\n", sep = "")
inv <- list.files(SUB, recursive = TRUE, full.names = TRUE)
if (length(inv)) {
  man <- data.frame(
    file = sub(paste0(normalizePath(SUB, winslash = "/"), "/"), "",
               normalizePath(inv, winslash = "/")),
    kb = round(file.size(inv) / 1024, 1), stringsAsFactors = FALSE)
  print(man, row.names = FALSE)
  write.csv(man, file.path(SUB, "MANIFEST.csv"), row.names = FALSE)
  cat(sprintf("\n%d files, %.1f MB total\n", nrow(man), sum(man$kb) / 1024))
} else cat("Folder is empty.\n")

cat("\nFolder:", normalizePath(SUB), "\n")
cat("\nStill to add by hand, since they are not produced by these scripts:\n")
cat("  the manuscript .docx, the Supplementary Material .docx, and the cover letter\n")
