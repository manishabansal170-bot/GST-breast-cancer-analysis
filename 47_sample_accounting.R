# =============================================================================
#  47_sample_accounting.R
#
#  THE OBJECTION
#    The manuscript uses at least eleven distinct denominators: 1,041, 1,083,
#    1,018, 783, 780, 770, 737, 726, 705, 488 and 74. Each is justified
#    somewhere, but no reader can hold them, and one of them (1,083) exceeds
#    the 1,041-tumour expression cohort used everywhere else without
#    explanation.
#
#    A reviewer has asked twice for a single accounting table. This builds it.
#
#  WHY IT IS COMPUTED RATHER THAN TYPED
#    The discrepancy the reviewer found arose because numbers were carried
#    between analyses by hand. Assembling the table the same way would risk
#    the same error. Every figure below is recomputed from the cached objects,
#    so the table states what the data contain rather than what the manuscript
#    remembers.
#
#  WHAT IT SHOULD SHOW
#    For each analysis: the starting set, each filter applied, how many samples
#    each filter removes, and the final denominator. Where two analyses differ
#    by a handful of samples, the table should make the reason visible without
#    the reader having to reconstruct it.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "supplementary_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB <- c("Basal","Her2","LumA","LumB")
id  <- function(v) substr(v, 1, 15)

gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }

rows <- list()
add <- function(analysis, cohort, n, note) {
  rows[[length(rows) + 1]] <<- data.frame(
    analysis = analysis, cohort = cohort, n = n, basis = note,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-42s %6d   %s\n", analysis, n, note))
}


# =============================================================================
# 1. THE TCGA EXPRESSION COHORT AND WHAT DERIVES FROM IT
# =============================================================================
cat(strrep("=", 96), "\nTCGA EXPRESSION\n", strrep("=", 96), "\n", sep = "")

raw  <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin  <- 2^gm(raw) - 0.001; lin[lin < 0] <- 0
E    <- log2(lin + 1); colnames(E) <- id(colnames(E))
ph   <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl   <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
st   <- TCGAquery_subtype(tumor = "brca")

brca_all <- id(ph[[cl("^sample$")]][
  toupper(ph[[cl("study")]]) == "TCGA" &
  grepl("breast", ph[[cl("primary disease or tissue")]], ignore.case = TRUE)])
add("TCGA breast samples in the compendium", "TCGA", length(unique(brca_all)),
    "all sample types")

prim <- id(ph[[cl("^sample$")]][
  toupper(ph[[cl("study")]]) == "TCGA" &
  grepl("breast", ph[[cl("primary disease or tissue")]], ignore.case = TRUE) &
  grepl("Primary Tumor", ph[[cl("sample_type")]])])
prim <- intersect(unique(prim), colnames(E))
add("  primary tumours with expression", "TCGA", length(prim), "excludes normals")

p50 <- st$BRCA_Subtype_PAM50[match(substr(prim, 1, 12), st$patient)]
add("  with any PAM50 call", "TCGA", sum(!is.na(p50)), "")
add("  in the four modelled subtypes", "TCGA", sum(p50 %in% SUB),
    "excludes normal-like and unassigned")

tum <- prim[p50 %in% SUB]
tum <- tum[!duplicated(substr(tum, 1, 12))]
add("EXPRESSION ARCHITECTURE", "TCGA", length(tum),
    "one sample per patient; Tables 1-2, Figures 1 and 4-6")

adjn <- id(ph[[cl("^sample$")]][
  toupper(ph[[cl("study")]]) == "TCGA" &
  grepl("breast", ph[[cl("primary disease or tissue")]], ignore.case = TRUE) &
  grepl("Solid Tissue Normal", ph[[cl("sample_type")]])])
adjn <- intersect(unique(adjn), colnames(E))
add("ADJACENT NORMAL", "TCGA", length(adjn), "the loss-axis comparator")

gtex <- id(ph[[cl("^sample$")]][
  toupper(ph[[cl("study")]]) == "GTEX" &
  grepl("breast", ph[[cl("primary disease or tissue")]], ignore.case = TRUE)])
add("GTEx NORMAL BREAST", "GTEx", length(intersect(unique(gtex), colnames(E))), "Figure 1")


# =============================================================================
# 2. SURVIVAL
# =============================================================================
cat("\n", strrep("=", 96), "\nSURVIVAL\n", strrep("=", 96), "\n", sep = "")
tc <- readRDS(file.path(CACHE, "tcga_clinical.rds"))
i  <- match(substr(tum, 1, 12), tc$submitter_id)
dead <- grepl("Dead|Deceased", tc$vital_status[i], ignore.case = TRUE)
tt <- suppressWarnings(ifelse(dead, as.numeric(tc$days_to_death[i]),
                              as.numeric(tc$days_to_last_follow_up[i]))) / 365.25
add("TCGA survival analysis", "TCGA", sum(is.finite(tt) & tt > 0),
    sprintf("%d excluded for absent or non-positive follow-up",
            length(tum) - sum(is.finite(tt) & tt > 0)))


# =============================================================================
# 3. THE METHYLATION SERIES
# =============================================================================
cat("\n", strrep("=", 96), "\nMETHYLATION\n", strrep("=", 96), "\n", sep = "")
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met; colnames(met) <- id(colnames(met))
exo <- obj$exp; colnames(exo) <- id(colnames(exo))

pb <- tryCatch({ B <- readRDS(file.path(CACHE, "all_gst_probe_betas.rds"))
                 length(intersect(id(colnames(B)), colnames(met))) },
               error = function(e) NA_integer_)
if (!is.na(pb)) add("probe-level beta values retrieved", "TCGA", pb,
                    "script 33, all promoter probes")

m1 <- colnames(met)[substr(colnames(met), 13, 15) == "-01"]
add("primary tumours with array methylation", "TCGA", length(m1), "")

m2 <- intersect(m1, colnames(exo))
add("  and expression on the array's own matrix", "TCGA", length(m2), "")

g2 <- st$BRCA_Subtype_PAM50[match(substr(m2, 1, 12), st$patient)]
add("  in the four modelled subtypes", "TCGA", sum(g2 %in% SUB),
    sprintf("%d normal-like, %d unassigned",
            sum(g2 == "Normal", na.rm = TRUE), sum(is.na(g2))))

cnv <- readRDS(file.path(CACHE, "tcga_cna_linear_gst.rds"))
if (is.list(cnv) && !is.matrix(cnv)) cnv <- cnv[[which(vapply(cnv, is.matrix, logical(1)))[1]]]
colnames(cnv) <- id(colnames(cnv))

m3 <- Reduce(intersect, list(m1, colnames(exo), colnames(cnv)))
g3 <- st$BRCA_Subtype_PAM50[match(substr(m3, 1, 12), st$patient)]
add("VARIANCE PARTITION, original matrix", "TCGA", length(m3),
    "methylation, expression and copy number")

m4 <- Reduce(intersect, list(m1, colnames(E), colnames(cnv)))
g4 <- st$BRCA_Subtype_PAM50[match(substr(m4, 1, 12), st$patient)]
add("VARIANCE PARTITION, harmonised matrix", "TCGA", sum(g4 %in% SUB),
    "as reported in Table 4")

m5 <- m1[st$BRCA_Subtype_PAM50[match(substr(m1, 1, 12), st$patient)] %in% SUB]
add("TABLE 3, methylation by subtype", "TCGA", length(m5), "")


# =============================================================================
# 4. THE ANALYSES WITH THEIR OWN DENOMINATORS
# =============================================================================
cat("\n", strrep("=", 96), "\nOTHER TCGA ANALYSES\n", strrep("=", 96), "\n", sep = "")

# immune: this is the one the reviewer queried
imm <- tryCatch(gm(readRDS(file.path(CACHE, "toil_immune_markers.rds"))), error = function(e) NULL)
pool <- tryCatch(gm(readRDS(file.path(CACHE, "purity_control_pool.rds"))), error = function(e) NULL)
if (!is.null(imm) && !is.null(pool)) {
  colnames(imm) <- id(colnames(imm)); colnames(pool) <- id(colnames(pool))
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  si <- Reduce(intersect, list(colnames(E), colnames(pool), colnames(imm)))
  si <- si[substr(si, 13, 15) == "-01"]
  si <- si[si %in% brca_all]
  pu <- num(tp$CPE)[match(substr(si, 1, 12), tp$patient)]
  si <- si[is.finite(pu)]
  gi <- st$BRCA_Subtype_PAM50[match(substr(si, 1, 12), st$patient)]
  add("IMMUNE AND PURITY", "TCGA", length(si),
      sprintf("not restricted to PAM50 subtypes: %d normal-like, %d unassigned included",
              sum(gi == "Normal", na.rm = TRUE), sum(is.na(gi))))
  add("  restricted to the four subtypes, for comparison", "TCGA", sum(gi %in% SUB),
      "shows why this differs from the expression cohort")
}

mi <- tryCatch({ M <- gm(readRDS(file.path(CACHE, "tcga_brca_mirna.rds")))
                 length(intersect(id(colnames(M)), tum)) },
               error = function(e) NA_integer_)
if (!is.na(mi)) add("MICRORNA", "TCGA", mi, "tumours with matched miRNA")


# =============================================================================
# 5. EXTERNAL COHORTS
# =============================================================================
cat("\n", strrep("=", 96), "\nEXTERNAL COHORTS\n", strrep("=", 96), "\n", sep = "")
# The cached objects hold every sample the source distributes. The manuscript
# reports figures after subtype filtering and, for survival, after requiring
# follow-up. Walking each cohort down its filters here prevents a reader
# meeting a third set of numbers.

# --- METABRIC -----------------------------------------------------------
tryCatch({
  M  <- gm(readRDS(file.path(CACHE, "metabric_gst.rds")))
  mc <- readRDS(file.path(CACHE, "metabric_clinical.rds"))
  sid <- grep("sampleId|SAMPLE_ID", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
  pc  <- grep("CLAUDIN_SUBTYPE", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
  mt  <- grep("OS_MONTHS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
  me  <- grep("OS_STATUS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]

  add("METABRIC, samples in the cached matrix", "METABRIC", ncol(M), "as distributed")
  sh <- intersect(colnames(M), as.character(mc[[sid]]))
  add("  with clinical annotation", "METABRIC", length(sh), "")
  g <- as.character(mc[[pc]][match(sh, mc[[sid]])])
  add("  in the four modelled subtypes", "METABRIC", sum(g %in% SUB),
      sprintf("%d in other classes or unassigned", length(sh) - sum(g %in% SUB)))
  sh2 <- sh[g %in% SUB]
  tm <- suppressWarnings(as.numeric(mc[[mt]][match(sh2, mc[[sid]])])) / 12
  ev <- as.numeric(grepl("^1|DECEASED|Died", as.character(mc[[me]][match(sh2, mc[[sid]])]),
                         ignore.case = TRUE))
  keep <- is.finite(tm) & tm > 0 & is.finite(ev)
  add("METABRIC SURVIVAL", "METABRIC", sum(keep),
      sprintf("%d deaths; excludes %d without positive follow-up",
              sum(ev[keep]), sum(g %in% SUB) - sum(keep)))
}, error = function(e) message("  METABRIC accounting failed: ", conditionMessage(e)))

# --- SCAN-B -------------------------------------------------------------
tryCatch({
  S  <- gm(readRDS(file.path(CACHE, "gse96058_gst.rds")))
  sp <- readRDS(file.path(CACHE, "gse96058_pheno.rds"))
  add("SCAN-B, samples in the cached matrix", "SCAN-B", ncol(S), "GSE96058, as distributed")
  pcol <- grep("pam50|subtype", colnames(sp), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(pcol)) { add("SCAN-B, subtype column not found", "SCAN-B", NA_integer_,
                         paste(head(colnames(sp), 6), collapse = ", ")) } else {
    scol <- grep("^sample|title|geo", colnames(sp), value = TRUE, ignore.case = TRUE)[1]
    g <- as.character(sp[[pcol]][match(colnames(S), as.character(sp[[scol]]))])
    g <- gsub("^\\s+|\\s+$", "", g)
    add("  with a subtype call", "SCAN-B", sum(!is.na(g) & g != ""), "")
    inSUB <- g %in% SUB | tolower(g) %in% tolower(c("Basal","Her2","LumA","LumB"))
    add("SCAN-B VALIDATION", "SCAN-B", sum(inSUB, na.rm = TRUE),
        sprintf("%d in other classes or unassigned", ncol(S) - sum(inSUB, na.rm = TRUE)))
  }
}, error = function(e) message("  SCAN-B accounting failed: ", conditionMessage(e)))

# --- neoadjuvant and proteomics ------------------------------------------
tryCatch({
  f1 <- list.files(CACHE, pattern = "GSE25066.*series_matrix", full.names = TRUE)[1]
  if (!is.na(f1) && requireNamespace("GEOquery", quietly = TRUE)) {
    suppressPackageStartupMessages(library(GEOquery))
    e1 <- getGEO(filename = f1, getGPL = FALSE)
    add("GSE25066, samples in the series", "GSE25066", ncol(Biobase::exprs(e1)), "")
  }
}, error = function(e) message("  GSE25066 accounting skipped: ", conditionMessage(e)))

tryCatch({
  pr <- read.csv(file.path(WORKDIR, "objective1", "moduleBG_alterations_protein",
                           "mrna_protein_concordance.csv"))
  add("PROTEOMICS", "TCGA", max(pr$n, na.rm = TRUE),
      sprintf("range %d to %d across proteins", min(pr$n), max(pr$n)))
}, error = function(e) message("  proteomics accounting skipped: ", conditionMessage(e)))


# =============================================================================
# 6. THE TABLE
# =============================================================================
tab <- bind_rows(rows)
write.csv(tab, file.path(OUT, "Supplementary_Table_S14_sample_accounting.csv"),
          row.names = FALSE)

cat("\n", strrep("=", 96), "\nWHAT TO DO WITH THIS\n", strrep("=", 96), "\n", sep = "")
cat("Deposit as Supplementary Table S14 and cite it once in Methods, at the point\n")
cat("where the first denominator appears. The reader who wants to reconcile two\n")
cat("numbers then has one place to look.\n\n")
cat("Check the immune row against the expression cohort. If it is larger, the\n")
cat("reason is that the immune analysis was not restricted to the four PAM50\n")
cat("subtypes, and the manuscript must say so where that number is reported.\n")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Rows: %d", nrow(tab)), "",
             capture.output(print(as.data.frame(tab), row.names = FALSE)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_sample_accounting.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
