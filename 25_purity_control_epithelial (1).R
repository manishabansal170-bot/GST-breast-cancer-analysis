# =============================================================================
#  25_purity_control_epithelial.R
#
#  THE WEAKNESS IN THE FIRST NEGATIVE CONTROL
#    Script 18 compared the GST family against 200 random gene sets matched
#    on expression level, and found that purity adjustment attenuated GST
#    immune correlations (mean |rho| 0.128 to 0.074) while leaving the random
#    sets untouched (0.126 to 0.127).
#
#    That comparison is softer than it looks. The control pool was drawn from
#    metabolic and housekeeping genes, which are expressed across all cell
#    types including immune cells. Such genes are not expected to be purity
#    sensitive, so the contrast partly restates the design rather than testing
#    it. A reviewer is right to ask for a control that is epithelial-restricted
#    like the GSTs, where purity adjustment SHOULD attenuate correlations.
#
#  THE QUESTION THIS SCRIPT ASKS
#    Among genes that are similarly epithelial-restricted, is the GST family
#    unusual? If GSTs attenuate like other epithelial genes, the finding is
#    about epithelial genes in bulk tissue generally, which is a real but
#    different claim. If GSTs attenuate more, the family is distinctive.
#
#    Either answer is publishable. The first requires rewording the manuscript
#    sentence; the second strengthens it. What is not defensible is leaving
#    the question unasked.
#
#  HOW EPITHELIAL RESTRICTION IS DEFINED HERE
#    Rather than importing an external list, restriction is measured in the
#    same data: a gene's correlation with tumour purity. Genes whose expression
#    rises with purity are, by construction, concentrated in the tumour
#    epithelial compartment. Control genes are matched to the GST family on
#    BOTH expression level and purity correlation, so the comparison holds
#    constant the property that drives the attenuation.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleG_immune")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

set.seed(20260824)
N_SETS <- 200
SUB <- c("Basal","Her2","LumA","LumB")

GST <- c("GSTA1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTM5","GSTP1",
         "GSTT2B","GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")


# =============================================================================
# 1. INPUTS, IDENTICAL TO SCRIPTS 14 AND 18
# =============================================================================
IMMUNE <- list(
  `T cells` = c("CD3D","CD3E","CD3G","CD2","CD6","SH2D1A","TRAT1","CD28","LCK"),
  `CD8 T cells` = c("CD8A","CD8B"),
  `Cytotoxic cells` = c("CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","NKG7","PRF1"),
  `Exhausted CD8` = c("CD244","EOMES","LAG3","PTGER4"),
  `Tregs` = c("FOXP3","IL2RA","IKZF2"),
  `Th1 cells` = c("TBX21","IFNG"),
  `B cells` = c("BLK","CD19","MS4A1","TNFRSF17","FCRL2","PNOC","SPIB","TCL1A"),
  `NK cells` = c("NCR1","XCL1","XCL2","KIR2DL3","KIR3DL1"),
  `Macrophages` = c("CD163","CD68","CD84","MS4A4A","MRC1"),
  `Dendritic cells` = c("CCL13","CD209","HSD11B1"),
  `Mast cells` = c("CPA3","HDC","MS4A2","TPSAB1","TPSB2"),
  `Neutrophils` = c("CSF3R","S100A12","CEACAM3","FCAR","FCGR3B","FPR1"),
  `Total immune` = c("PTPRC"))
CHECKPOINT <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
                "IDO1","BTLA","VSIR","TNFRSF9","ICOS")
STROMA <- c("FAP","PDGFRB","PDGFRA","ACTA2","COL1A1","COL1A2","COL3A1","THY1","POSTN","VIM")
CYTOKINE <- c("IFNG","TNF","IL1B","IL6","IL10","TGFB1","CXCL9","CXCL10","CXCL13","CCL5","GZMK")

mk <- readRDS(file.path(CACHE, "toil_immune_markers.rds"))
lin <- 2^mk - 0.001; lin[lin < 0] <- 0
M <- log2(lin + 1); M <- M[apply(M, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]

gr <- readRDS(file.path(CACHE, "toil_gst.rds"))
gl <- 2^gr - 0.001; gl[gl < 0] <- 0
G_all <- log2(gl + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]], study = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")
clin <- ph %>% dplyr::filter(toupper(study) == "TCGA",
                             grepl("breast", tissue, ignore.case = TRUE),
                             grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

pool_raw <- readRDS(file.path(CACHE, "purity_control_pool.rds"))   # from script 18
pl <- 2^pool_raw - 0.001; pl[pl < 0] <- 0
P_all <- log2(pl + 1)

s <- Reduce(intersect, list(clin$sample, colnames(M), colnames(G_all), colnames(P_all)))
clin <- clin[match(s, clin$sample), ]
M <- M[, s, drop = FALSE]
G <- G_all[intersect(GST, rownames(G_all)), s, drop = FALSE]
P <- P_all[, s, drop = FALSE]
P <- P[apply(P, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]

clin$purity <- NA_real_
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment()); tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  clin$purity <- num(tp$CPE)[match(clin$patient, tp$patient)] }, silent = TRUE)

cat("Tumours:", length(s), "| with purity:", sum(!is.na(clin$purity)),
    "| control pool:", nrow(P), "\n")

score_of <- function(g, mat) {
  gg <- intersect(g, rownames(mat)); if (!length(gg)) return(rep(NA_real_, ncol(mat)))
  colMeans(t(scale(t(mat[gg, , drop = FALSE]))), na.rm = TRUE) }
S <- sapply(names(IMMUNE), function(k) score_of(IMMUNE[[k]], M))
S <- cbind(S, Checkpoint = score_of(CHECKPOINT, M),
           `Stromal / CAF` = score_of(STROMA, M),
           `Inflammatory cytokines` = score_of(CYTOKINE, M))
rownames(S) <- s; S <- S[, colSums(is.na(S)) < nrow(S), drop = FALSE]


# =============================================================================
# 2. MEASURE EPITHELIAL RESTRICTION AS CORRELATION WITH PURITY
# =============================================================================
pur_cor <- function(v) {
  ok <- is.finite(v) & is.finite(clin$purity)
  if (sum(ok) < 50) return(NA_real_)
  suppressWarnings(cor(v[ok], clin$purity[ok], method = "spearman")) }

gst_pc  <- apply(G, 1, pur_cor)
pool_pc <- apply(P, 1, pur_cor)
gst_med <- apply(G, 1, median, na.rm = TRUE)
pool_med <- apply(P, 1, median, na.rm = TRUE)

cat("\n===== EPITHELIAL RESTRICTION, AS CORRELATION WITH TUMOUR PURITY =====\n")
cat(sprintf("GST family      : median rho %+.3f  (range %+.3f to %+.3f)\n",
            median(gst_pc, na.rm = TRUE), min(gst_pc, na.rm = TRUE), max(gst_pc, na.rm = TRUE)))
cat(sprintf("Control pool    : median rho %+.3f  (range %+.3f to %+.3f)\n",
            median(pool_pc, na.rm = TRUE), min(pool_pc, na.rm = TRUE), max(pool_pc, na.rm = TRUE)))

ok_pool <- names(pool_pc)[is.finite(pool_pc) & is.finite(pool_med)]
cat(sprintf("\nPool genes with purity correlation at least as high as the GST median: %d\n",
            sum(pool_pc[ok_pool] >= median(gst_pc, na.rm = TRUE))))


# =============================================================================
# 3. THE PIPELINE, UNCHANGED FROM SCRIPT 14
# =============================================================================
partial_rho <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 50) return(c(NA, NA, NA))
  raw <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  rx <- residuals(lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- residuals(lm(rank(y[ok]) ~ rank(z[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman"))
  c(raw, unname(ct$estimate), ct$p.value) }

run_set <- function(mat) {
  res <- lapply(rownames(mat), function(g) lapply(colnames(S), function(p) {
    v <- partial_rho(mat[g, ], S[, p], clin$purity)
    if (all(is.na(v))) return(NULL)
    data.frame(rho_raw = v[1], rho_adj = v[2], p_adj = v[3]) }) %>% bind_rows()) %>% bind_rows()
  if (!nrow(res)) return(NULL)
  res$FDR <- p.adjust(res$p_adj, "BH")
  data.frame(n_tested = nrow(res),
             pct_survive = 100 * mean(abs(res$rho_adj) >= 0.3 & res$FDR < 0.05, na.rm = TRUE),
             mean_raw = mean(abs(res$rho_raw), na.rm = TRUE),
             mean_adj = mean(abs(res$rho_adj), na.rm = TRUE),
             mean_atten = mean(abs(res$rho_raw) - abs(res$rho_adj), na.rm = TRUE)) }

obs <- run_set(G)
cat("\n===== OBSERVED, GST FAMILY =====\n"); print(obs, digits = 3)


# =============================================================================
# 4. MATCH ON BOTH EXPRESSION AND PURITY CORRELATION
# =============================================================================
# Each GST is matched to a pool gene close in BOTH properties, using a
# standardised distance so neither dominates.
z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
all_med <- c(gst_med, pool_med[ok_pool]); all_pc <- c(gst_pc, pool_pc[ok_pool])
zmed <- z(all_med); zpc <- z(all_pc)

draw_matched <- function() {
  avail <- ok_pool; chosen <- character(0)
  for (g in names(gst_med)) {
    if (!length(avail)) break
    dist <- sqrt((zmed[avail] - zmed[g])^2 + (zpc[avail] - zpc[g])^2)
    cand <- avail[order(dist)][seq_len(min(10, length(avail)))]
    pick <- sample(cand, 1); chosen <- c(chosen, pick); avail <- setdiff(avail, pick) }
  P[chosen, , drop = FALSE] }

chk <- draw_matched()
cat(sprintf("\nExample matched set: median purity rho %+.3f against GST %+.3f\n",
            median(pur_cor2 <- apply(chk, 1, pur_cor), na.rm = TRUE),
            median(gst_pc, na.rm = TRUE)))

cat("\nRunning", N_SETS, "purity-matched sets. A few minutes.\n")
null <- vector("list", N_SETS)
for (i in seq_len(N_SETS)) {
  null[[i]] <- run_set(draw_matched())
  if (i %% 20 == 0) message("  ", i, " / ", N_SETS) }
null <- bind_rows(null)


# =============================================================================
# 5. RESULT
# =============================================================================
q  <- quantile(null$pct_survive, c(0.025, 0.5, 0.975), na.rm = TRUE)
qa <- quantile(null$mean_atten,  c(0.025, 0.5, 0.975), na.rm = TRUE)

cat("\n", strrep("=", 70), "\nGST FAMILY AGAINST PURITY-MATCHED CONTROLS\n",
    strrep("=", 70), "\n", sep = "")
cat(sprintf("Surviving associations   GST %.2f%%   controls %.2f%% (95%% %.2f to %.2f)\n",
            obs$pct_survive, q[2], q[1], q[3]))
cat(sprintf("Mean attenuation         GST %.3f    controls %.3f (95%% %.3f to %.3f)\n",
            obs$mean_atten, qa[2], qa[1], qa[3]))
cat(sprintf("Empirical p, attenuation at least as large: %.3f\n",
            mean(null$mean_atten >= obs$mean_atten, na.rm = TRUE)))

cat("\n----- INTERPRETATION -----\n")
if (obs$mean_atten >= qa[1] && obs$mean_atten <= qa[3]) {
  cat("GST attenuation lies within the range for equally epithelial-restricted\n")
  cat("genes. The loss of immune correlations under purity adjustment is a\n")
  cat("property of epithelial genes in bulk tissue, not something particular to\n")
  cat("this family. The manuscript should say so: the correlations reflect\n")
  cat("cellular composition, as would be expected for any epithelial gene set.\n")
} else if (obs$mean_atten > qa[3]) {
  cat("GST attenuation exceeds that of equally epithelial-restricted controls.\n")
  cat("The family is unusually composition-driven even by that standard, and\n")
  cat("the claim of specificity survives a stricter test than before.\n")
} else {
  cat("GST correlations attenuate LESS than epithelial-matched controls, so\n")
  cat("more of their immune association survives adjustment than expected.\n")
  cat("That inverts the current framing and must be reported.\n")
}

write.csv(null, file.path(OUT, "purity_control_epithelial_null.csv"), row.names = FALSE)
write.csv(cbind(obs, gst_median_purity_rho = median(gst_pc, na.rm = TRUE)),
          file.path(OUT, "purity_control_epithelial_observed.csv"), row.names = FALSE)

p1 <- ggplot(null, aes(mean_atten)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40", linewidth = .2) +
  geom_vline(xintercept = obs$mean_atten, colour = "#C0392B", linewidth = 1) +
  annotate("text", x = obs$mean_atten, y = Inf, label = "  GST family",
           hjust = 0, vjust = 2, colour = "#C0392B", size = 3.4, fontface = "bold") +
  labs(x = "Mean attenuation of immune correlation after purity adjustment",
       y = "Matched control gene sets",
       title = "GST family against equally epithelial-restricted controls",
       subtitle = sprintf("%d sets matched on expression level and on correlation with tumour purity",
                          N_SETS)) +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "purity_control_epithelial.png"), p1,
       width = 7.5, height = 4.5, dpi = 300, bg = "white")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("GST median purity rho: %.3f", median(gst_pc, na.rm = TRUE)),
             sprintf("GST attenuation %.3f; control median %.3f (95%% %.3f to %.3f)",
                     obs$mean_atten, qa[2], qa[1], qa[3]),
             sprintf("Empirical p: %.3f", mean(null$mean_atten >= obs$mean_atten, na.rm = TRUE)),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_purity_epithelial.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
