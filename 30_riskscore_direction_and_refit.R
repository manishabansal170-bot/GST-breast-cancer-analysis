# =============================================================================
#  30_riskscore_direction_and_refit.R
#
#  TWO OBJECTIONS ABOUT THE SAME RISK SCORE
#
#  C4  THE COEFFICIENT IS LOCKED BUT THE SCALING IS NOT
#      The manuscript describes applying TCGA-derived coefficients unchanged to
#      METABRIC, which is the right thing to do. But the score standardises
#      GSTM4 using the TCGA mean and standard deviation, and METABRIC is a
#      microarray cohort on a different scale entirely. Transferring a
#      standardisation between platforms with different dynamic ranges is not
#      the same as transferring a coefficient, and the referee is right that it
#      needs justifying or testing.
#
#  C6  THE DIRECTION IS ARBITRARY
#      TCGA has 141 deaths and 2.3 years median follow-up. METABRIC has 965
#      deaths and 9.7 years. Training on the smaller, shorter-followed cohort
#      and testing on the larger is the harder direction, and defensible, but
#      the paper does not say why that direction was chosen. Reversing it is a
#      one-line change and either supports the result or does not.
#
#  WHAT IS TESTED HERE
#      1. The published direction, as a baseline
#      2. Within-METABRIC standardisation instead of transferred TCGA scaling
#      3. Rank-based scoring, which is scale-free and sidesteps the issue
#      4. The reverse direction: LASSO on METABRIC, applied to TCGA
#
#  WHAT WOULD MATTER
#      The manuscript's claim is that the family carries no clinically useful
#      prognostic signal. If the score is significant in every configuration
#      but with a small C-index throughout, that claim strengthens: the result
#      is not an artefact of one arbitrary choice. If significance appears and
#      disappears depending on the scaling, the analysis is fragile and the
#      negative claim needs to be stated more carefully, not less.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleD_survival")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","survival","glmnet")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

set.seed(20260825)
SUB <- c("Basal","Her2","LumA","LumB")
PUBLISHED_GENE <- "GSTM4"
PUBLISHED_COEF <- -0.0437


# =============================================================================
# 1. BUILD BOTH COHORTS, MATCHING SCRIPT 16
# =============================================================================
get_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("No matrix in cached object") }

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^get_matrix(raw) - 0.001; lin[lin < 0] <- 0
E <- log2(lin + 1)

st <- TCGAquery_subtype(tumor = "brca")
tc <- readRDS(file.path(CACHE, "tcga_clinical.rds"))
samp <- grep("^TCGA-.*-01$", colnames(E), value = TRUE)
d_t <- data.frame(sample = samp, patient = substr(samp, 1, 12),
                  stringsAsFactors = FALSE)
d_t <- d_t[!duplicated(d_t$patient), ]
i <- match(d_t$patient, tc$submitter_id)
d_t <- d_t[!is.na(i), ]; i <- i[!is.na(i)]
dead <- grepl("Dead|Deceased", tc$vital_status[i], ignore.case = TRUE)
d_t$time <- suppressWarnings(ifelse(dead, as.numeric(tc$days_to_death[i]),
                                    as.numeric(tc$days_to_last_follow_up[i]))) / 365.25
d_t$event <- as.numeric(dead)
d_t$pam50 <- st$BRCA_Subtype_PAM50[match(d_t$patient, st$patient)]
d_t <- d_t[is.finite(d_t$time) & d_t$time > 0 & d_t$pam50 %in% SUB, ]

M  <- get_matrix(readRDS(file.path(CACHE, "metabric_gst.rds")))
mc <- readRDS(file.path(CACHE, "metabric_clinical.rds"))
sid <- grep("sampleId|SAMPLE_ID", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
mt  <- grep("OS_MONTHS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
me  <- grep("OS_STATUS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
pc  <- grep("CLAUDIN_SUBTYPE", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
sh <- intersect(colnames(M), as.character(mc[[sid]]))
d_m <- data.frame(sample = sh, stringsAsFactors = FALSE)
j <- match(d_m$sample, mc[[sid]])
d_m$time  <- suppressWarnings(as.numeric(mc[[mt]][j])) / 12
d_m$event <- as.numeric(grepl("^1|DECEASED|Died", as.character(mc[[me]][j]), ignore.case = TRUE))
d_m$pam50 <- as.character(mc[[pc]][j])
d_m <- d_m[is.finite(d_m$time) & d_m$time > 0 & d_m$pam50 %in% SUB, ]

genes <- intersect(rownames(E), rownames(M))
cat("TCGA:", nrow(d_t), "patients,", sum(d_t$event), "deaths\n")
cat("METABRIC:", nrow(d_m), "patients,", sum(d_m$event), "deaths\n")
cat("Genes in both cohorts:", length(genes), "\n")

Xt <- t(E[genes, d_t$sample, drop = FALSE])
Xm <- t(M[genes, d_m$sample, drop = FALSE])
cat(sprintf("\nScale check, %s: TCGA mean %.2f sd %.2f | METABRIC mean %.2f sd %.2f\n",
            PUBLISHED_GENE, mean(Xt[, PUBLISHED_GENE]), sd(Xt[, PUBLISHED_GENE]),
            mean(Xm[, PUBLISHED_GENE]), sd(Xm[, PUBLISHED_GENE])))
cat("These are different measurement scales. Transferring a standardisation\n")
cat("between them is the choice under test.\n")


# =============================================================================
# 2. EVALUATION HELPER
# =============================================================================
evaluate <- function(score, dd, label) {
  ok <- is.finite(score) & is.finite(dd$time) & is.finite(dd$event)
  cx <- coxph(Surv(dd$time[ok], dd$event[ok]) ~ scale(score[ok]))
  s <- summary(cx)
  data.frame(configuration = label, n = sum(ok), events = sum(dd$event[ok]),
             HR_perSD = round(exp(coef(cx)), 3),
             CI_low  = round(exp(confint(cx))[1], 3),
             CI_high = round(exp(confint(cx))[2], 3),
             p = s$coefficients[1, 5],
             C_index = round(s$concordance[1], 3))
}


# =============================================================================
# 3. C4: THREE SCALINGS OF THE SAME LOCKED COEFFICIENT
# =============================================================================
cat("\n", strrep("=", 70), "\nC4. SCALING OF THE LOCKED COEFFICIENT\n", strrep("=", 70), "\n", sep = "")

g <- PUBLISHED_GENE
res <- list()

# as published: TCGA mean and sd applied to METABRIC
s_pub <- PUBLISHED_COEF * ((Xm[, g] - mean(Xt[, g])) / sd(Xt[, g]))
res[["published"]] <- evaluate(s_pub, d_m, "As published, TCGA scaling transferred")

# within-cohort standardisation
s_own <- PUBLISHED_COEF * scale(Xm[, g])[, 1]
res[["within"]] <- evaluate(s_own, d_m, "Within-METABRIC standardisation")

# rank-based, scale free
s_rank <- PUBLISHED_COEF * scale(rank(Xm[, g]))[, 1]
res[["rank"]] <- evaluate(s_rank, d_m, "Rank-based, scale free")

out <- bind_rows(res)
print(as.data.frame(out), row.names = FALSE, digits = 3)

cat("\nThe hazard ratios and C-indices should be near-identical across the three,\n")
cat("because a linear rescaling of a single-gene score cannot change the ordering\n")
cat("of patients, and Cox regression on a standardised score depends only on that\n")
cat("ordering. If they differ materially, something else is wrong.\n")


# =============================================================================
# 4. C6: THE REVERSE DIRECTION
# =============================================================================
cat("\n", strrep("=", 70), "\nC6. TRAINING ON METABRIC, TESTING ON TCGA\n", strrep("=", 70), "\n", sep = "")

fit_lasso <- function(X, dd, label) {
  y <- Surv(dd$time, dd$event)
  Z <- scale(X)
  cvf <- cv.glmnet(Z, y, family = "cox", alpha = 1, nfolds = 10)
  co <- as.matrix(coef(cvf, s = "lambda.min"))
  sel <- co[co[, 1] != 0, , drop = FALSE]
  cat(sprintf("\nLASSO on %s selected %d gene(s) at lambda.min:\n", label, nrow(sel)))
  if (nrow(sel)) print(round(sel, 4)) else cat("  none\n")
  list(coef = setNames(co[, 1], rownames(co)),
       centre = attr(Z, "scaled:center"), scale = attr(Z, "scaled:scale"))
}

# reproduce the published direction first, to show the pipeline agrees
fit_t <- fit_lasso(Xt, d_t, "TCGA")
sel_t <- names(fit_t$coef)[fit_t$coef != 0]

# reverse
fit_m <- fit_lasso(Xm, d_m, "METABRIC")
sel_m <- names(fit_m$coef)[fit_m$coef != 0]

cat("\nGenes selected in TCGA:   ", if (length(sel_t)) paste(sel_t, collapse = ", ") else "none", "\n")
cat("Genes selected in METABRIC:", if (length(sel_m)) paste(sel_m, collapse = ", ") else "none", "\n")
cat("In both:", paste(intersect(sel_t, sel_m), collapse = ", "), "\n")

rev_out <- NULL
if (length(sel_m)) {
  # apply METABRIC-derived model to TCGA, coefficients and scaling locked
  sc <- as.numeric(scale(Xt[, names(fit_m$coef), drop = FALSE],
                         center = fit_m$centre, scale = fit_m$scale) %*% fit_m$coef)
  rev_out <- evaluate(sc, d_t, "METABRIC-trained, tested in TCGA")
  print(as.data.frame(rev_out), row.names = FALSE, digits = 3)
} else {
  cat("\nLASSO selected no gene in METABRIC, so there is no model to transfer.\n")
  cat("That is itself an answer: in the larger, longer-followed cohort the\n")
  cat("penalised model retains nothing, which supports the negative claim more\n")
  cat("directly than the published direction does.\n")
}


# =============================================================================
# 5. WHAT THIS MEANS
# =============================================================================
cat("\n", strrep("=", 70), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 70), "\n", sep = "")
spread <- diff(range(out$C_index, na.rm = TRUE))
cat(sprintf("C-index across the three scalings: %.3f to %.3f, spread %.3f\n",
            min(out$C_index), max(out$C_index), spread))
if (spread < 0.01) {
  cat("The result does not depend on the scaling choice. State in the Methods\n")
  cat("that the score is a single standardised gene, so any linear rescaling\n")
  cat("leaves the patient ordering and therefore the Cox model unchanged.\n")
} else {
  cat("The result varies with the scaling. This must be reported and the\n")
  cat("published configuration justified.\n")
}

all_out <- bind_rows(out, rev_out)
write.csv(all_out, file.path(OUT, "riskscore_direction_and_scaling.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("TCGA %d patients %d deaths; METABRIC %d patients %d deaths",
                     nrow(d_t), sum(d_t$event), nrow(d_m), sum(d_m$event)),
             sprintf("TCGA LASSO selected: %s", paste(sel_t, collapse = ", ")),
             sprintf("METABRIC LASSO selected: %s",
                     if (length(sel_m)) paste(sel_m, collapse = ", ") else "none"),
             sprintf("C-index spread across scalings: %.3f", spread),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_riskscore_direction.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
