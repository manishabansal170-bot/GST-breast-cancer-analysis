# =============================================================================
#  44_methylation_purity_confounding.R
#
#  THE OBJECTION, AND WHY IT IS THE SHARPEST ONE LEFT
#    This manuscript adjusts expression for tumour purity throughout and argues
#    the case at length. The methylation side receives no equivalent treatment,
#    and a reviewer points out that the same argument applies with more force.
#
#    Basal-like tumours are systematically less pure and more heavily
#    infiltrated. Normal breast epithelium and leukocytes are unmethylated at
#    the GSTP1 CpG island. A bulk beta value from a basal-like tumour is
#    therefore diluted toward zero by its own non-tumour content, and the
#    97.0% unmethylated figure could in principle be a composition artefact
#    rather than a property of the cancer cells.
#
#    This is exactly the dilution argument the manuscript rebuts for
#    expression. Applying it there and not here is inconsistent, and the
#    reviewer is right to say the linear model on rank expression does not
#    address it, because there the confounded variable is the outcome rather
#    than the predictor.
#
#  WHAT THIS SCRIPT DOES
#    1. Reports purity by subtype, so the size of the imbalance is visible
#    2. Correlates GSTP1 beta against purity directly
#    3. Repeats the Table 3 subtype comparison restricted to high-purity
#       tumours, where dilution cannot be doing the work
#    4. Repeats it with beta adjusted for purity by rank residualisation
#    5. Extends the same treatment to the five universal-loss genes, comparing
#       adjacent normal against tumour before and after adjustment
#
#  HOW TO READ IT
#    If the basal-versus-luminal contrast holds among high-purity tumours
#    alone, dilution is not generating it and the claim is secure. If it
#    collapses, the manuscript must say that the subtype difference in bulk
#    beta partly reflects differing non-tumour content, which would be a
#    genuine and reportable limitation of every 450k tumour study of this kind.
#
#    For the loss axis the reasoning is the same but the comparison is
#    adjacent normal against tumour, and the manuscript already makes the
#    stromal-origin argument for GSTM5. The question is whether the other four
#    behave as GSTM5 does.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB  <- c("Basal","Her2","LumA","LumB")
LOSS <- c("GSTM5","GSTM2","GSTA1","MGST1","GSTA4")


# =============================================================================
# 1. INPUTS
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
met <- obj$met
colnames(met) <- substr(colnames(met), 1, 15)

gm <- function(x) { if (is.matrix(x)) return(x)
  if (is.list(x)) { if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1))); if (length(i)) return(x[[i[1]]]) }
  stop("no matrix") }
raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^gm(raw) - 0.001; lin[lin < 0] <- 0
expr <- log2(lin + 1); colnames(expr) <- substr(colnames(expr), 1, 15)

st <- TCGAquery_subtype(tumor = "brca")
data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
tp <- get("Tumor.purity", envir = environment())
tp$patient <- substr(tp$Sample.ID, 1, 12)
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

s <- colnames(met)[substr(colnames(met), 13, 15) == "-01"]
grp <- st$BRCA_Subtype_PAM50[match(substr(s, 1, 12), st$patient)]
ok <- grp %in% SUB
s <- s[ok]; grp <- factor(grp[ok], levels = SUB)
pur <- num(tp$CPE)[match(substr(s, 1, 12), tp$patient)]

cat("Tumours with methylation and a PAM50 call:", length(s), "\n")
cat("Of those, with a purity estimate:", sum(is.finite(pur)), "\n\n")


# =============================================================================
# 2. IS PURITY IMBALANCED BY SUBTYPE
# =============================================================================
cat(strrep("=", 78), "\nA. TUMOUR PURITY BY SUBTYPE\n", strrep("=", 78), "\n", sep = "")
pt <- data.frame(subtype = grp, purity = pur)
pt <- pt[is.finite(pt$purity), ]
print(pt %>% group_by(subtype) %>%
        summarise(n = n(), median_purity = round(median(purity), 3),
                  Q1 = round(quantile(purity, .25), 3),
                  Q3 = round(quantile(purity, .75), 3), .groups = "drop"))
kw <- kruskal.test(purity ~ subtype, data = pt)
cat(sprintf("\nKruskal-Wallis across subtypes: p = %.3g\n", kw$p.value))
cat("If basal-like tumours are systematically less pure, the dilution concern is\n")
cat("real and must be addressed rather than assumed away.\n")


# =============================================================================
# 3. DOES GSTP1 BETA TRACK PURITY
# =============================================================================
cat("\n", strrep("=", 78), "\nB. GSTP1 BETA AGAINST PURITY\n", strrep("=", 78), "\n", sep = "")
b <- as.numeric(met["GSTP1", s])
okp <- is.finite(b) & is.finite(pur)
cat(sprintf("Overall: Spearman rho = %.3f (n = %d)\n",
            cor(b[okp], pur[okp], method = "spearman"), sum(okp)))
cat("\nWithin each subtype, where the subtype contrast cannot contribute:\n")
for (k in SUB) {
  i <- okp & grp == k
  if (sum(i) < 20) next
  cat(sprintf("  %-6s rho = %+.3f  (n = %d)\n", k,
              cor(b[i], pur[i], method = "spearman"), sum(i)))
}
cat("\nA positive correlation means purer tumours are more methylated, which is\n")
cat("what dilution by unmethylated normal tissue would produce.\n")


# =============================================================================
# 4. TABLE 3 RESTRICTED TO HIGH-PURITY TUMOURS
# =============================================================================
cat("\n", strrep("=", 78), "\nC. THE SUBTYPE CONTRAST IN HIGH-PURITY TUMOURS ONLY\n",
    strrep("=", 78), "\n", sep = "")

summarise_beta <- function(idx, label) {
  data.frame(
    analysis = label,
    subtype  = SUB,
    n = as.integer(table(droplevels(grp[idx]))[SUB]),
    median_beta = round(tapply(b[idx], droplevels(grp[idx]), median, na.rm = TRUE)[SUB], 3),
    pct_methylated = round(tapply(b[idx], droplevels(grp[idx]),
                          function(x) 100 * mean(x > 0.3, na.rm = TRUE))[SUB], 1),
    pct_unmethylated = round(tapply(b[idx], droplevels(grp[idx]),
                          function(x) 100 * mean(x < 0.1, na.rm = TRUE))[SUB], 1),
    stringsAsFactors = FALSE)
}

all_t <- summarise_beta(okp, "all tumours")
print(all_t, row.names = FALSE)

for (thr in c(0.6, 0.7)) {
  hi <- okp & pur >= thr
  cat(sprintf("\nRestricted to purity >= %.1f (%d tumours):\n", thr, sum(hi)))
  if (sum(hi) < 80) { cat("  too few tumours for a stable comparison\n"); next }
  print(summarise_beta(hi, sprintf("purity >= %.1f", thr)), row.names = FALSE)
}


# =============================================================================
# 5. BETA ADJUSTED FOR PURITY
# =============================================================================
cat("\n", strrep("=", 78), "\nD. SUBTYPE CONTRAST WITH BETA ADJUSTED FOR PURITY\n",
    strrep("=", 78), "\n", sep = "")
cat("Beta is residualised on purity by rank regression, then compared across\n")
cat("subtypes. This asks whether the contrast survives once the part of beta\n")
cat("explained by non-tumour content is removed.\n\n")

br <- resid(lm(rank(b[okp]) ~ rank(pur[okp])))
g2 <- droplevels(grp[okp])
print(data.frame(subtype = levels(g2),
                 n = as.integer(table(g2)),
                 mean_residual_rank = round(tapply(br, g2, mean), 1)),
      row.names = FALSE)
k2 <- kruskal.test(br ~ g2)
cat(sprintf("\nKruskal-Wallis on purity-adjusted beta: p = %.3g\n", k2$p.value))
k1 <- kruskal.test(b[okp] ~ g2)
cat(sprintf("Unadjusted, for comparison            : p = %.3g\n", k1$p.value))


# =============================================================================
# 6. THE LOSS AXIS, ADJACENT NORMAL AGAINST TUMOUR
# =============================================================================
cat("\n", strrep("=", 78), "\nE. THE LOSS AXIS, ADJUSTED FOR PURITY\n",
    strrep("=", 78), "\n", sep = "")
cat("Adjacent normal is epithelium-rich and tumour is roughly 70% pure, so the\n")
cat("same dilution argument applies to every gene on the loss axis, not only to\n")
cat("GSTM5 where the manuscript already makes it.\n\n")

# expression-side check: does the tumour-normal difference survive adjustment
ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
adjn <- ph[[cl("^sample$")]][
  toupper(ph[[cl("study")]]) == "TCGA" &
  grepl("breast", ph[[cl("primary disease or tissue")]], ignore.case = TRUE) &
  grepl("Solid Tissue Normal", ph[[cl("sample_type")]])]
adjn <- substr(adjn, 1, 15)
adjn <- intersect(adjn, colnames(expr))
tum  <- intersect(s, colnames(expr))
cat("Adjacent normal:", length(adjn), "| tumours:", length(tum), "\n\n")

cliffs <- function(a, b) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  r <- rank(c(a, b)); U <- sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2
  2 * U / (length(a) * length(b)) - 1
}

out <- lapply(LOSS, function(g) {
  if (!(g %in% rownames(expr))) return(NULL)
  tv <- as.numeric(expr[g, tum]); nv <- as.numeric(expr[g, adjn])
  d_all <- cliffs(tv, nv)
  # high-purity tumours only, against the same normals
  tp2 <- pur[match(tum, s)]
  hi <- is.finite(tp2) & tp2 >= 0.7
  d_hi <- if (sum(hi) >= 50) cliffs(tv[hi], nv) else NA_real_
  data.frame(gene = g,
             delta_all = round(d_all, 3),
             delta_high_purity = round(d_hi, 3),
             n_high_purity = sum(hi),
             change = round(d_hi - d_all, 3), stringsAsFactors = FALSE)
}) %>% bind_rows()
print(as.data.frame(out), row.names = FALSE, digits = 3)

cat("\nIf the effect sizes are similar in high-purity tumours, dilution is not\n")
cat("producing the loss and Table 2 stands as written. If they shrink toward\n")
cat("zero, the loss is partly a composition effect and the manuscript should\n")
cat("say so for those genes, as it already does for GSTM5.\n")


# =============================================================================
# 7. WHAT TO WRITE
# =============================================================================
cat("\n", strrep("=", 78), "\nCONSEQUENCE FOR THE MANUSCRIPT\n", strrep("=", 78), "\n", sep = "")
cat("Report section A regardless: if purity differs by subtype, the reader needs\n")
cat("to know before reading Table 3.\n\n")
cat("If the contrast in section C holds among high-purity tumours, state that\n")
cat("and cite the restricted numbers. That is the cleanest possible answer to\n")
cat("this objection and costs nothing.\n\n")
cat("If it weakens, report the restricted figures alongside the full ones and\n")
cat("qualify the 97.0% accordingly. A bulk beta value cannot distinguish an\n")
cat("unmethylated cancer cell from an unmethylated lymphocyte, and no 450k\n")
cat("tumour study can fully escape that.\n")

write.csv(out, file.path(OUT, "loss_axis_purity_sensitivity.csv"), row.names = FALSE)
writeLines(c(paste("Run:", Sys.time()),
             sprintf("Tumours: %d, with purity: %d", length(s), sum(is.finite(pur))),
             sprintf("GSTP1 beta against purity, overall rho: %.3f",
                     cor(b[okp], pur[okp], method = "spearman")),
             sprintf("Subtype contrast p, unadjusted %.3g, purity-adjusted %.3g",
                     k1$p.value, k2$p.value),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_methylation_purity.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
