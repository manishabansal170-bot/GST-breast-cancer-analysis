# =============================================================================
#  35_response_multivariable.R
#
#  THE OBJECTION
#    The neoadjuvant response models are stratified by oestrogen receptor
#    status and otherwise unadjusted. Grade, nodal status, HER2 status and age
#    all predict pathological complete response, and all correlate with the
#    subtype structure this paper documents. An association between GSTM2 and
#    response could therefore reflect any of them rather than the gene.
#
#    The reviewer asks for the models to be refitted with these as covariates,
#    in both cohorts.
#
#  WHAT THIS SCRIPT DOES
#    Extracts whatever clinical covariates each GEO series actually carries,
#    reports which were found, and refits the logistic models with them. Where
#    a covariate is absent it says so rather than silently omitting it, because
#    an adjustment that could not be made is a limitation to state, not a gap
#    to leave invisible.
#
#  WHAT WOULD CHANGE THE PAPER
#    If GSTM2 survives adjustment in the discovery cohort, the exploratory
#    claim stands as reported. If it does not, the association was confounded
#    by clinical variables and the honest statement is that it does not survive
#    adjustment, which is a cleaner outcome than the current "tested and not
#    reproduced".
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleE_drug_response")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- c("GEOquery","hgu133a.db","hgu133plus2.db","AnnotationDbi")
miss <- bioc[!sapply(bioc, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(GEOquery); library(AnnotationDbi)
  library(hgu133a.db); library(hgu133plus2.db) })

FOCUS <- c("GSTM1","GSTM2")


# =============================================================================
# HELPERS
# =============================================================================
# Pull a labelled field out of the characteristics columns, whatever its
# position, and return the text after the colon.
grab <- function(pd, pattern) {
  ch <- pd[, grep("characteristics_ch1", colnames(pd)), drop = FALSE]
  apply(ch, 1, function(r) {
    h <- grep(pattern, r, value = TRUE, ignore.case = TRUE)
    if (!length(h)) return(NA_character_)
    trimws(sub("^[^:]*:\\s*", "", h[1]))
  })
}

collapse_probes <- function(ex, db, genes) {
  map <- AnnotationDbi::select(db, keys = rownames(ex), columns = "SYMBOL",
                               keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL %in% genes, ]
  if (!nrow(map)) return(NULL)
  em <- ex[map$PROBEID, , drop = FALSE]
  v <- apply(em, 1, var, na.rm = TRUE)
  idx <- tapply(seq_len(nrow(map)), map$SYMBOL, function(i) i[which.max(v[i])])
  out <- em[unlist(idx), , drop = FALSE]; rownames(out) <- names(idx); out
}

fit_models <- function(E, cl, cohort) {
  cat("\n", strrep("=", 74), "\n", cohort, "\n", strrep("=", 74), "\n", sep = "")

  covars <- c("grade","node","her2","age")
  present <- covars[sapply(covars, function(v)
    sum(!is.na(cl[[v]])) >= 0.7 * nrow(cl) && length(unique(na.omit(cl[[v]]))) > 1)]
  cat("Covariates available:", if (length(present)) paste(present, collapse = ", ") else "none", "\n")
  absent <- setdiff(covars, present)
  if (length(absent)) cat("Not recorded in this series:", paste(absent, collapse = ", "), "\n")

  out <- list()
  for (er in c("ER-positive","ER-negative")) {
    dd <- cl[which(cl$er_pos == (er == "ER-positive")), ]
    if (nrow(dd) < 40) next
    cat(sprintf("\n%s stratum: n = %d, %d complete responses (%.1f%%)\n",
                er, nrow(dd), sum(dd$pcr), 100 * mean(dd$pcr)))
    for (g in intersect(FOCUS, rownames(E))) {
      x <- as.numeric(E[g, dd$sample])
      base <- data.frame(y = dd$pcr, x = scale(x)[, 1])
      for (v in present) base[[v]] <- dd[[v]]
      base <- base[complete.cases(base), ]
      if (nrow(base) < 40 || length(unique(base$y)) < 2) next

      m0 <- glm(y ~ x, data = base, family = binomial)
      f1 <- as.formula(paste("y ~ x", if (length(present))
                             paste("+", paste(present, collapse = " + ")) else ""))
      m1 <- glm(f1, data = base, family = binomial)
      s0 <- summary(m0)$coefficients["x", ]; s1 <- summary(m1)$coefficients["x", ]
      ci <- function(s) c(exp(s[1] - 1.96 * s[2]), exp(s[1] + 1.96 * s[2]))
      c0 <- ci(s0); c1 <- ci(s1)
      out[[length(out) + 1]] <- data.frame(
        cohort = cohort, ER = er, gene = g, n = nrow(base), events = sum(base$y),
        OR_unadjusted = round(exp(s0[1]), 3),
        unadj_low = round(c0[1], 3), unadj_high = round(c0[2], 3), p_unadj = s0[4],
        OR_adjusted = round(exp(s1[1]), 3),
        adj_low = round(c1[1], 3), adj_high = round(c1[2], 3), p_adj = s1[4],
        covariates = if (length(present)) paste(present, collapse = "+") else "none",
        stringsAsFactors = FALSE)
    }
  }
  bind_rows(out)
}


# =============================================================================
# 1. GSE25066, DISCOVERY
# =============================================================================
f1 <- list.files(CACHE, pattern = "GSE25066.*series_matrix", full.names = TRUE)[1]
if (is.na(f1)) stop("GSE25066 series matrix not found in cache.")
e1 <- getGEO(filename = f1, getGPL = FALSE)
ex1 <- Biobase::exprs(e1); pd1 <- Biobase::pData(e1)
if (max(ex1, na.rm = TRUE) > 100) ex1 <- log2(ex1 + 1)

cl1 <- data.frame(sample = rownames(pd1), stringsAsFactors = FALSE)
cl1$resp  <- grab(pd1, "pathologic.?response|pcr|response")
cl1$er    <- grab(pd1, "^er.?status|estrogen|esr1")
cl1$grade <- grab(pd1, "grade")
cl1$node  <- grab(pd1, "node|nodal|^n.?stage")
cl1$her2  <- grab(pd1, "her2|erbb2")
cl1$age   <- suppressWarnings(as.numeric(grab(pd1, "age")))

cat("Raw label examples, GSE25066:\n")
for (v in c("resp","er","grade","node","her2")) {
  tb <- table(cl1[[v]], useNA = "ifany")
  cat(" ", v, ":", paste(head(names(tb), 6), collapse = " | "), "\n")
}

cl1$pcr <- NA
cl1$pcr[grepl("^pCR|complete", cl1$resp, ignore.case = TRUE)] <- 1
cl1$pcr[grepl("^RD|residual|non.?pCR|^nCR$", cl1$resp, ignore.case = TRUE)] <- 0
cl1$er_pos <- NA
cl1$er_pos[grepl("^P$|^pos|positive", cl1$er, ignore.case = TRUE)] <- TRUE
cl1$er_pos[grepl("^N$|^neg|negative", cl1$er, ignore.case = TRUE)] <- FALSE
cl1 <- cl1[!is.na(cl1$pcr) & !is.na(cl1$er_pos) & cl1$sample %in% colnames(ex1), ]
cat("\nGSE25066 analysable:", nrow(cl1), "\n")

E1 <- collapse_probes(ex1, hgu133a.db, FOCUS)
r1 <- if (!is.null(E1)) fit_models(E1, cl1, "GSE25066 (discovery)") else NULL


# =============================================================================
# 2. GSE32646, VALIDATION
# =============================================================================
f2 <- list.files(CACHE, pattern = "GSE32646.*series_matrix", full.names = TRUE)
f2 <- f2[file.size(f2) > 1000][1]
r2 <- NULL
if (!is.na(f2)) {
  e2 <- getGEO(filename = f2, getGPL = FALSE)
  ex2 <- Biobase::exprs(e2); pd2 <- Biobase::pData(e2)
  if (max(ex2, na.rm = TRUE) > 100) ex2 <- log2(ex2 + 1)
  cl2 <- data.frame(sample = rownames(pd2), stringsAsFactors = FALSE)
  cl2$resp  <- grab(pd2, "pathologic.?response|pcr|response")
  cl2$er    <- grab(pd2, "^er.?status|estrogen|esr1")
  cl2$grade <- grab(pd2, "grade")
  cl2$node  <- grab(pd2, "node|nodal")
  cl2$her2  <- grab(pd2, "her2|erbb2")
  cl2$age   <- suppressWarnings(as.numeric(grab(pd2, "age")))
  cl2$pcr <- NA
  cl2$pcr[cl2$resp == "pCR"] <- 1
  cl2$pcr[cl2$resp == "nCR"] <- 0
  cl2$er_pos <- NA
  cl2$er_pos[grepl("^P$|^pos|positive", cl2$er, ignore.case = TRUE)] <- TRUE
  cl2$er_pos[grepl("^N$|^neg|negative", cl2$er, ignore.case = TRUE)] <- FALSE
  cl2 <- cl2[!is.na(cl2$pcr) & !is.na(cl2$er_pos) & cl2$sample %in% colnames(ex2), ]
  cat("\nGSE32646 analysable:", nrow(cl2), "\n")
  E2 <- collapse_probes(ex2, hgu133plus2.db, FOCUS)
  r2 <- if (!is.null(E2)) fit_models(E2, cl2, "GSE32646 (validation)") else NULL
} else {
  cat("\nGSE32646 series matrix not available at usable size; validation skipped.\n")
}


# =============================================================================
# 3. RESULT
# =============================================================================
res <- bind_rows(r1, r2)
if (!nrow(res)) stop("No models fitted.")

cat("\n", strrep("=", 74), "\nUNADJUSTED AGAINST ADJUSTED\n", strrep("=", 74), "\n", sep = "")
print(as.data.frame(res[, c("cohort","ER","gene","n","events",
                            "OR_unadjusted","OR_adjusted","adj_low","adj_high","p_adj")]),
      row.names = FALSE, digits = 3)

key <- res[res$ER == "ER-positive" & res$gene == "GSTM2" &
             grepl("discovery", res$cohort), ]
if (nrow(key)) {
  cat("\n", strrep("=", 74), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 74), "\n", sep = "")
  cat(sprintf("GSTM2, receptor-positive, discovery cohort:\n"))
  cat(sprintf("  unadjusted OR %.3f (%.3f to %.3f)\n",
              key$OR_unadjusted, key$unadj_low, key$unadj_high))
  cat(sprintf("  adjusted   OR %.3f (%.3f to %.3f), p = %.4g\n",
              key$OR_adjusted, key$adj_low, key$adj_high, key$p_adj))
  cat(sprintf("  covariates: %s\n\n", key$covariates))
  if (key$p_adj < 0.05 && key$OR_adjusted < 1) {
    cat("The association survives adjustment. Report both estimates and name the\n")
    cat("covariates; the exploratory claim stands as written.\n")
  } else if (key$OR_adjusted < 1) {
    cat("The direction holds but significance does not survive adjustment. Report\n")
    cat("both estimates and state that the association is not independent of the\n")
    cat("clinical covariates available in this cohort.\n")
  } else {
    cat("The association does not survive adjustment. This is a cleaner negative\n")
    cat("than the failed replication and should replace it as the primary reason\n")
    cat("the claim is not carried forward.\n")
  }
}

write.csv(res, file.path(OUT, "response_multivariable.csv"), row.names = FALSE)
writeLines(c(paste("Run:", Sys.time()),
             sprintf("Models fitted: %d", nrow(res)),
             sprintf("Covariates used: %s", paste(unique(res$covariates), collapse = "; ")),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_response_multivariable.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
