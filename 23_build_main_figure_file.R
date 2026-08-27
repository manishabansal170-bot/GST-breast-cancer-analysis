# =============================================================================
#  23_build_main_figure_file.R
#
#  Produces one multi-page PDF containing all six main figures, one per page
#  with its caption, plus individual TIFFs for journals that want them
#  separately.
#
#  THE PROBLEM THIS SOLVES
#    Figures 1 and 2 exist as standalone files because scripts 16 and 21 write
#    them. Figures 3, 4, 5A and 5B exist only as images embedded in the
#    manuscript, with no separate copies. A .docx is a zip archive, so those
#    can be extracted rather than regenerated.
#
#  A LIMIT WORTH UNDERSTANDING BEFORE YOU RELY ON THIS
#    Extraction cannot recover resolution that was lost when the image was
#    embedded. If a figure went into Word at 150 dpi, that is what comes out.
#    The script reports the effective resolution of every figure at its printed
#    size and flags any that fall below 300 dpi, which is the usual journal
#    threshold. Anything flagged must be regenerated from the script that made
#    it, not extracted.
#
#  WHERE EACH FIGURE COMES FROM, IF REGENERATION IS NEEDED
#    Figure 1  script 16, Fig_panel_toil
#    Figure 2  script 16, Panel_SCANB
#    Figure 3  methylation module, methylation and expression by subtype
#    Figure 4  methylation module, GSTP1 methylation by subtype
#    Figure 5A pathway module, co-expression with functional processes
#    Figure 5B pathway module, NRF2 activity by subtype
# =============================================================================

WORKDIR <- "~/GST_BRCA"
OUT     <- file.path(WORKDIR, "submission", "main_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

# Point this at whichever manuscript version you are submitting.
DOCX <- file.path(WORKDIR, "GST_manuscript_v28.docx")
if (!file.exists(DOCX)) {
  hits <- list.files(c(WORKDIR, path.expand("~/Downloads"),
                       "C:/Users/Manisha Gupta/Downloads"),
                     pattern = "GST_manuscript.*\\.docx$", full.names = TRUE)
  hits <- hits[!grepl("^~\\$", basename(hits))]     # ignore Word lock files
  if (!length(hits)) stop("Manuscript .docx not found. Set DOCX manually.")
  DOCX <- hits[order(file.mtime(hits), decreasing = TRUE)][1]
  message("Using most recent manuscript found: ", DOCX)
}

pk <- c("png", "magick")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
have_magick <- requireNamespace("magick", quietly = TRUE)
if (!have_magick) stop("The magick package is required. install.packages('magick')")
suppressPackageStartupMessages(library(magick))


# -----------------------------------------------------------------------------
# 1. EXTRACT THE EMBEDDED IMAGES
# -----------------------------------------------------------------------------
tmp <- file.path(tempdir(), "docx_media")
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
utils::unzip(DOCX, exdir = tmp)

media <- list.files(file.path(tmp, "word", "media"), full.names = TRUE)
cat("Images found inside the manuscript:", length(media), "\n")
for (m in sort(media)) {
  inf <- image_info(image_read(m))
  cat(sprintf("  %-14s %5d x %5d\n", basename(m), inf$width, inf$height))
}


# -----------------------------------------------------------------------------
# 2. MAP IMAGES TO FIGURE NUMBERS
#    Order follows the manuscript. image6 is Figure 2 because it was added
#    after the others, so the numbering of the media files does not match the
#    numbering of the figures.
#    Standalone regenerated files are preferred where they exist, because they
#    carry full resolution; the embedded copy is the fallback.
# -----------------------------------------------------------------------------
FIGS <- list(
  list(id = "Figure_1",  media = "image1.png",
       prefer = c("figures/Fig_panel_toil.tiff", "figures/Fig_panel_toil.png"),
       cap = "Figure 1. GST family expression across normal breast and TCGA-BRCA molecular subtypes. GTEx normal breast (n = 179), TCGA adjacent normal (n = 113), and primary tumours by PAM50: basal-like (n = 190), HER2-enriched (n = 81), luminal A (n = 561), luminal B (n = 209). Values are log2(TPM + 1) from the UCSC Toil recompute."),
  list(id = "Figure_2",  media = "image6.png",
       prefer = c("figures/Panel_SCANB.tiff", "figures/scanb/Panel_SCANB.tiff",
                  "figures/Panel_SCANB.png"),
       cap = "Figure 2. GST family expression across PAM50 subtypes in SCAN-B (GSE96058). Independent Swedish cohort, RNA-seq, n = 3,184. Panels ordered by subtype effect size. GSTT1 and GSTT2 are omitted as not quantifiable; GSTA2, GSTA3 and GSTA5 are shown but fall below the detection threshold in every subtype."),
  list(id = "Figure_3",  media = "image2.png", prefer = character(0),
       cap = "Figure 3. Methylation and expression across PAM50 subtypes. Promoter methylation (upper row) and expression (lower row) for the six genes with subtype-variable methylation."),
  list(id = "Figure_4",  media = "image3.png", prefer = character(0),
       cap = "Figure 4. GSTP1 promoter methylation by PAM50 subtype. Dashed line, beta = 0.3, the conventional threshold for a methylated promoter."),
  list(id = "Figure_5A", media = "image4.png", prefer = character(0),
       cap = "Figure 5A. Co-expression of GST genes with ten functional processes in TCGA-BRCA tumours, as mean Spearman correlation with the constituent genes of each process."),
  list(id = "Figure_5B", media = "image5.png", prefer = character(0),
       cap = "Figure 5B. NRF2 pathway activity across PAM50 subtypes, computed as the mean z-score of 20 canonical NRF2 target genes (Kruskal-Wallis p = 1.87e-20).")
)


# -----------------------------------------------------------------------------
# 3. RESOLVE EACH FIGURE TO A FILE, AND ASSESS RESOLUTION
# -----------------------------------------------------------------------------
# Printed width assumed to be a single journal column at 180 mm, roughly 7.1 in.
PRINT_IN <- 7.1

cat("\n", strrep("=", 70), "\nRESOLUTION AUDIT\n", strrep("=", 70), "\n", sep = "")

resolved <- list()
for (f in FIGS) {
  src <- NA_character_; how <- ""
  for (p in f$prefer) if (file.exists(p)) { src <- p; how <- "standalone"; break }
  if (is.na(src)) {
    cand <- file.path(tmp, "word", "media", f$media)
    if (file.exists(cand)) { src <- cand; how <- "extracted" }
  }
  if (is.na(src)) { cat(sprintf("  %-10s NOT FOUND\n", f$id)); next }

  img <- image_read(src)
  inf <- image_info(img)
  dpi <- round(inf$width / PRINT_IN)
  flag <- if (dpi >= 300) "ok" else if (dpi >= 200) "MARGINAL" else "TOO LOW"
  cat(sprintf("  %-10s %-11s %5d x %5d  %4d dpi at %.1f in  %s\n",
              f$id, how, inf$width, inf$height, dpi, PRINT_IN, flag))

  resolved[[f$id]] <- list(id = f$id, img = img, cap = f$cap,
                           dpi = dpi, flag = flag, how = how, src = src)
}

low <- names(Filter(function(x) x$flag != "ok", resolved))
if (length(low)) {
  cat("\nBelow 300 dpi at printed size:", paste(low, collapse = ", "), "\n")
  cat("These cannot be improved by extraction, because the resolution was\n")
  cat("lost when the image was embedded. Regenerate them from the script that\n")
  cat("produced them, saving with dpi = 300, before submitting.\n")
} else {
  cat("\nAll figures meet 300 dpi at the assumed printed width.\n")
}


# -----------------------------------------------------------------------------
# 4. INDIVIDUAL TIFFs
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\nWRITING INDIVIDUAL FILES\n", strrep("=", 70), "\n", sep = "")
for (r in resolved) {
  dst <- file.path(OUT, paste0(r$id, ".tiff"))
  image_write(r$img, dst, format = "tiff", compression = "LZW")
  cat(sprintf("  %s  (%.1f MB)\n", basename(dst), file.size(dst) / 1024^2))
}


# -----------------------------------------------------------------------------
# 5. THE COMBINED FILE
#    One figure per page, caption beneath, in manuscript order.
# -----------------------------------------------------------------------------
PAGE_W <- 2550   # 8.5 in at 300 dpi
PAGE_H <- 3300   # 11 in at 300 dpi
MARGIN <- 225    # 0.75 in

pages <- lapply(resolved, function(r) {
  avail_w <- PAGE_W - 2 * MARGIN
  avail_h <- PAGE_H - 2 * MARGIN - 400        # room for the caption
  im <- image_resize(r$img, paste0(avail_w, "x", avail_h))
  canvas <- image_blank(PAGE_W, PAGE_H, color = "white")
  inf <- image_info(im)
  x <- round((PAGE_W - inf$width) / 2)
  canvas <- image_composite(canvas, im, offset = paste0("+", x, "+", MARGIN))
  wrapped <- paste(strwrap(r$cap, width = 105), collapse = "\n")
  image_annotate(canvas, wrapped, size = 34, color = "black", font = "sans",
                 location = paste0("+", MARGIN, "+", MARGIN + inf$height + 60))
})

combined <- image_join(pages)
pdf_out <- file.path(OUT, "Main_Figures_1_to_5.pdf")
image_write(combined, pdf_out, format = "pdf", density = "300x300")

cat("\n", strrep("=", 70), "\nDONE\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("Combined file: %s  (%.1f MB, %d pages)\n",
            basename(pdf_out), file.size(pdf_out) / 1024^2, length(pages)))
cat("Folder:", normalizePath(OUT), "\n")

writeLines(c(paste("Run:", Sys.time()),
             paste("Source manuscript:", DOCX),
             paste("Assumed printed width (in):", PRINT_IN),
             "",
             "Figure, source, pixels, dpi at printed width, verdict",
             sapply(resolved, function(r) {
               inf <- image_info(r$img)
               sprintf("%s, %s, %dx%d, %d, %s", r$id, r$how, inf$width, inf$height,
                       r$dpi, r$flag) }),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_main_figures.txt"))

cat("\nCheck the resolution audit above before submitting. Journals reject\n")
cat("figure files below their stated threshold, and it is a slow rejection.\n")
