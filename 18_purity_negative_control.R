# =============================================================================
#  18_purity_negative_control.R
#
#  THE QUESTION THIS ANSWERS
#    Script 14 reports that 242 of 256 GST-immune associations (94.5%) do not
#    survive adjustment for tumour purity. That number is only interpretable
#    against a baseline. Three readings are possible and the manuscript cannot
#    currently distinguish them:
#
#      (a) GSTs specifically are not immune-associated
#      (b) bulk immune deconvolution is broadly confounded by composition, so
#          almost any gene set would behave this way
#      (c) the consensus purity estimate over-corrects
#
#    This script tests (a) against (b) by running the identical pipeline on
#    random gene sets matched to the GST family for expression level, and
#    asking what fraction of THEIR associations survive.
#
#  WHY NOT A CROSS-COHORT REPLICATION INSTEAD
#    METABRIC and SCAN-B have no consensus purity estimate and no methylation
#    data. The only available substitute is ESTIMATE, which derives purity
#    from immune and stromal signatures - so adjusting immune correlations for
#    it removes the signal by construction. That analysis would be circular
#    and is deliberately not attempted here.
#
#  INTERPRETING THE OUTPUT
#    If random sets survive at a similar rate to GSTs, the 94.5% figure is a
#    property of bulk deconvolution, not of GST biology, and the manuscript
#    sentence must be reworded.
#    If random sets survive at a substantially higher rate, the GST result is
#    specific and the claim stands as written.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleG_immune")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(TCGAbiolinks) })

set.seed(20260823)          # reproducible sampling
N_SETS   <- 200             # random gene sets to draw
POOL_N   <- 400             # candidate genes to fetch for the pool

GST <- c("GSTA1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTM5","GSTP1",
         "GSTT2B","GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. REBUILD THE EXACT INPUTS USED BY SCRIPT 14
# =============================================================================
IMMUNE <- list(
  `T cells`          = c("CD3D","CD3E","CD3G","CD2","CD6","SH2D1A","TRAT1","CD28","LCK"),
  `CD8 T cells`      = c("CD8A","CD8B"),
  `Cytotoxic cells`  = c("CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","NKG7","PRF1"),
  `Exhausted CD8`    = c("CD244","EOMES","LAG3","PTGER4"),
  `Tregs`            = c("FOXP3","IL2RA","IKZF2"),
  `Th1 cells`        = c("TBX21","IFNG"),
  `B cells`          = c("BLK","CD19","MS4A1","TNFRSF17","FCRL2","PNOC","SPIB","TCL1A"),
  `NK cells`         = c("NCR1","XCL1","XCL2","KIR2DL3","KIR3DL1"),
  `Macrophages`      = c("CD163","CD68","CD84","MS4A4A","MRC1"),
  `Dendritic cells`  = c("CCL13","CD209","HSD11B1"),
  `Mast cells`       = c("CPA3","HDC","MS4A2","TPSAB1","TPSB2"),
  `Neutrophils`      = c("CSF3R","S100A12","CEACAM3","FCAR","FCGR3B","FPR1"),
  `Total immune`     = c("PTPRC")
)
CHECKPOINT <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
                "IDO1","BTLA","VSIR","TNFRSF9","ICOS")
STROMA <- c("FAP","PDGFRB","PDGFRA","ACTA2","COL1A1","COL1A2","COL3A1","THY1",
            "POSTN","VIM")
CYTOKINE <- c("IFNG","TNF","IL1B","IL6","IL10","TGFB1","CXCL9","CXCL10",
              "CXCL13","CCL5","GZMK")

marker_raw <- readRDS(file.path(CACHE, "toil_immune_markers.rds"))
lin <- 2^marker_raw - 0.001; lin[lin < 0] <- 0
M <- log2(lin + 1)
M <- M[apply(M, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]

gst_raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
gl <- 2^gst_raw - 0.001; gl[gl < 0] <- 0
G_all <- log2(gl + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

clin <- ph %>%
  dplyr::filter(toupper(study) == "TCGA",
                grepl("breast", tissue, ignore.case = TRUE),
                grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

s <- Reduce(intersect, list(clin$sample, colnames(M), colnames(G_all)))
clin <- clin[match(s, clin$sample), ]
M <- M[, s, drop = FALSE]
G <- G_all[intersect(GST, rownames(G_all)), s, drop = FALSE]

clin$purity <- NA_real_
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  clin$purity <- num(tp$CPE)[match(clin$patient, tp$patient)]
}, silent = TRUE)

cat("Tumours:", length(s), "| with purity:", sum(!is.na(clin$purity)), "\n")

score_of <- function(genes, mat) {
  gg <- intersect(genes, rownames(mat))
  if (!length(gg)) return(rep(NA_real_, ncol(mat)))
  colMeans(t(scale(t(mat[gg, , drop = FALSE]))), na.rm = TRUE)
}
S <- sapply(names(IMMUNE), function(k) score_of(IMMUNE[[k]], M))
S <- cbind(S,
           `Checkpoint` = score_of(CHECKPOINT, M),
           `Stromal / CAF` = score_of(STROMA, M),
           `Inflammatory cytokines` = score_of(CYTOKINE, M))
rownames(S) <- s
S <- S[, colSums(is.na(S)) < nrow(S), drop = FALSE]
cat("Populations scored:", ncol(S), "\n")


# =============================================================================
# 2. THE SAME PARTIAL CORRELATION FUNCTION AS SCRIPT 14
# =============================================================================
partial_rho <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 50) return(c(NA, NA, NA))
  raw <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  rx <- residuals(lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- residuals(lm(rank(y[ok]) ~ rank(z[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman"))
  c(raw, unname(ct$estimate), ct$p.value)
}

# Run one gene set through the full pipeline and return summary counts.
run_set <- function(mat) {
  res <- lapply(rownames(mat), function(g) {
    lapply(colnames(S), function(p) {
      v <- partial_rho(mat[g, ], S[, p], clin$purity)
      if (all(is.na(v))) return(NULL)
      data.frame(rho_raw = v[1], rho_adj = v[2], p_adj = v[3])
    }) %>% bind_rows()
  }) %>% bind_rows()
  if (!nrow(res)) return(NULL)
  res$FDR <- p.adjust(res$p_adj, "BH")
  data.frame(
    n_tested    = nrow(res),
    n_survive   = sum(abs(res$rho_adj) >= 0.3 & res$FDR < 0.05, na.rm = TRUE),
    pct_survive = 100 * mean(abs(res$rho_adj) >= 0.3 & res$FDR < 0.05, na.rm = TRUE),
    mean_atten  = mean(abs(res$rho_raw) - abs(res$rho_adj), na.rm = TRUE),
    mean_raw    = mean(abs(res$rho_raw), na.rm = TRUE),
    mean_adj    = mean(abs(res$rho_adj), na.rm = TRUE))
}

obs <- run_set(G)
cat("\n===== OBSERVED, GST FAMILY =====\n")
print(obs, digits = 3)
cat(sprintf("Not surviving: %.1f%%\n", 100 - obs$pct_survive))


# =============================================================================
# 3. BUILD A MATCHED RANDOM GENE POOL
# =============================================================================
# Random genes must be matched on expression level. A set of highly expressed
# genes and a set of near-silent genes have very different correlation
# behaviour, so an unmatched comparison would be uninformative.

poolf <- file.path(CACHE, "purity_control_pool.rds")
if (file.exists(poolf)) {
  pool_raw <- readRDS(poolf); message("Loaded control gene pool from cache.")
} else {
  # A broad, deliberately arbitrary set of protein-coding genes spanning the
  # expression range. Not curated for any biological property.
  CAND <- c(
    "AARS1","ABCB6","ABCC5","ACAA2","ACADM","ACAT1","ACLY","ACO2","ACTR2","ADSL",
    "AK2","ALDH3A2","ALG3","ANAPC5","AP2M1","APEX1","ARF4","ARPC2","ASNS","ATG3",
    "ATIC","ATP5F1B","ATP6V1A","AURKB","BAX","BCAP31","BLVRA","BRIX1","BUB3","CAD",
    "CANX","CAPZA1","CBX3","CCT2","CCT5","CDC20","CDK4","CENPA","CHCHD2","CHMP4B",
    "CKAP4","CLTB","CNOT1","COPB1","COX5A","CPSF1","CS","CSNK1A1","CTPS1","CUL1",
    "DAD1","DARS1","DCTN2","DDOST","DDX21","DHX15","DLAT","DNAJA1","DPM1","DUT",
    "EEF1B2","EIF2S1","EIF3E","EIF4A1","ENO1","EPRS1","ERH","ETF1","EXOSC7","FARSA",
    "FDPS","FEN1","FH","FKBP4","G6PD","GAPDH","GART","GLA","GLRX3","GMPS",
    "GNB1","GOT2","GPI","GSR","GTF2B","HADHA","HAT1","HDAC1","HMGCS1","HNRNPA1",
    "HPRT1","HSD17B10","HSPA4","HSPD1","IDH1","IDH3A","ILF2","IMPDH2","IPO4","ITGB1",
    "KARS1","KIF11","KPNB1","LAMP1","LDHA","LSM4","MAD2L1","MAPK1","MCM2","MDH1",
    "ME1","MRPL13","MRPS12","MTHFD1","MYC","NAMPT","NCAPD2","NDUFA4","NME1","NOP56",
    "NPM1","NUDT1","OAT","P4HB","PAICS","PCNA","PDHA1","PFKP","PGAM1","PGD",
    "PGK1","PHB1","PLK1","PMM2","POLD1","POLR2A","PPA1","PPIA","PRDX1","PRDX4",
    "PSMA1","PSMB5","PSMD1","PTGES3","PYCR1","RAD51","RAN","RANBP1","RBM8A","RFC4",
    "RPA1","RPL13A","RPN1","RRM1","RRM2","RUVBL1","SDHA","SEC61A1","SERPINH1","SF3B1",
    "SHMT2","SLC1A5","SLC2A1","SMC2","SNRPB","SOD1","SRM","SRSF1","SSRP1","STIP1",
    "SUPT16H","TALDO1","TARS1","TCP1","TFRC","THBS1","TIMM50","TK1","TKT","TMED2",
    "TOMM40","TPI1","TRAP1","TUBB","TXN","TXNRD1","TYMS","UBE2C","UBE2L3","UMPS",
    "UNG","UQCRC1","USP1","VARS1","VCP","VDAC1","XPO1","XRCC5","YARS1","YWHAE",
    "AASS","ABAT","ACADS","ACSS2","ADH5","AGPS","AKR1A1","ALAS1","ALDH2","AOX1",
    "ARG2","ASS1","BCKDHA","BDH1","CAT","CBR1","CES1","CHDH","CPT1A","CYB5A",
    "DAO","DBT","DECR1","DHRS4","ECHS1","EPHX1","ESD","ETFA","FMO2","GAMT",
    "GATM","GCAT","GCDH","GLUD1","GLUL","GPX3","GSTCD","HAAO","HGD","HIBCH",
    "HMGCL","HPD","IVD","KMO","MAOA","MCEE","MGLL","MTHFD2","NNMT","NQO1"
  )
  CAND <- setdiff(unique(CAND), c(GST, rownames(M)))
  CAND <- head(CAND, POOL_N)
  HOST <- "https://toil.xenahubs.net"
  EXPR <- "TcgaTargetGtex_rsem_gene_tpm"
  got <- list(); failed <- character(0)
  for (g in CAND) {
    r <- tryCatch(
      fetch_dense_values(host = HOST, dataset = EXPR, identifiers = g,
                         use_probeMap = TRUE, check = FALSE),
      error = function(e) NULL)
    if (is.null(r) || !nrow(as.matrix(r))) { failed <- c(failed, g); next }
    m <- as.matrix(r)
    got[[g]] <- if (nrow(m) > 1) colMeans(m, na.rm = TRUE) else m[1, ]
    Sys.sleep(0.15)
    if (length(got) %% 25 == 0) message("  fetched ", length(got), " / ", length(CAND))
  }
  common <- Reduce(intersect, lapply(got, names))
  pool_raw <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(pool_raw) <- names(got)
  saveRDS(pool_raw, poolf)
  message("Pool fetched: ", nrow(pool_raw), ". Failed: ", length(failed))
}

pl <- 2^pool_raw - 0.001; pl[pl < 0] <- 0
P <- log2(pl + 1)
P <- P[, intersect(colnames(P), s), drop = FALSE]
P <- P[, match(s, colnames(P)), drop = FALSE]
P <- P[apply(P, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]
cat("\nControl pool genes available:", nrow(P), "\n")


# =============================================================================
# 4. EXPRESSION-MATCHED SAMPLING
# =============================================================================
# Match each GST to a pool gene of similar median expression, sampling without
# replacement within the set. Repeat N_SETS times.

gst_med  <- apply(G, 1, median, na.rm = TRUE)
pool_med <- apply(P, 1, median, na.rm = TRUE)

draw_matched <- function() {
  avail <- names(pool_med); chosen <- character(0)
  for (target in gst_med) {
    if (!length(avail)) break
    d <- abs(pool_med[avail] - target)
    # sample from the 10 nearest to avoid always picking the same gene
    cand <- avail[order(d)][seq_len(min(10, length(avail)))]
    pick <- sample(cand, 1)
    chosen <- c(chosen, pick); avail <- setdiff(avail, pick)
  }
  P[chosen, , drop = FALSE]
}

cat("\nRunning", N_SETS, "matched random sets. This takes a few minutes.\n")
null <- vector("list", N_SETS)
for (i in seq_len(N_SETS)) {
  null[[i]] <- run_set(draw_matched())
  if (i %% 20 == 0) message("  ", i, " / ", N_SETS)
}
null <- bind_rows(null)


# =============================================================================
# 5. RESULT
# =============================================================================
cat("\n===== NULL DISTRIBUTION, EXPRESSION-MATCHED RANDOM SETS =====\n")
print(summary(null$pct_survive))
q <- quantile(null$pct_survive, c(0.025, 0.5, 0.975), na.rm = TRUE)
cat(sprintf("\nRandom sets, %% of associations surviving purity adjustment:\n"))
cat(sprintf("  median %.2f%%   95%% range %.2f%% to %.2f%%\n", q[2], q[1], q[3]))
cat(sprintf("Observed GST family: %.2f%%\n", obs$pct_survive))

emp_p <- mean(null$pct_survive <= obs$pct_survive, na.rm = TRUE)
cat(sprintf("\nProportion of random sets with survival at or below GST: %.3f\n", emp_p))

cat("\nMean attenuation (raw minus adjusted, absolute rho):\n")
cat(sprintf("  random sets median %.3f   GST family %.3f\n",
            median(null$mean_atten, na.rm = TRUE), obs$mean_atten))

cat("\n----- INTERPRETATION -----\n")
if (obs$pct_survive >= q[1] && obs$pct_survive <= q[3]) {
  cat("The GST result falls INSIDE the range expected for arbitrary\n")
  cat("expression-matched gene sets. The 94.5% figure therefore describes\n")
  cat("bulk deconvolution in general, not GST biology specifically. The\n")
  cat("manuscript sentence should be reworded to say so.\n")
} else if (obs$pct_survive < q[1]) {
  cat("The GST family survives purity adjustment LESS often than arbitrary\n")
  cat("gene sets. The association with immune content is weaker than\n")
  cat("baseline, which strengthens the negative claim as written.\n")
} else {
  cat("The GST family survives purity adjustment MORE often than arbitrary\n")
  cat("gene sets, so the surviving associations are not merely residual\n")
  cat("composition. Report the comparison alongside the 94.5% figure.\n")
}

write.csv(null, file.path(OUT, "purity_negative_control_null.csv"), row.names = FALSE)
write.csv(obs,  file.path(OUT, "purity_negative_control_observed.csv"), row.names = FALSE)

p1 <- ggplot(null, aes(pct_survive)) +
  geom_histogram(bins = 30, fill = "grey80", colour = "grey40", linewidth = .2) +
  geom_vline(xintercept = obs$pct_survive, colour = "#C0392B", linewidth = 1) +
  annotate("text", x = obs$pct_survive, y = Inf, label = "  GST family",
           hjust = 0, vjust = 2, colour = "#C0392B", size = 3.4, fontface = "bold") +
  labs(x = "% of gene-population associations surviving purity adjustment",
       y = "Random matched gene sets",
       title = "Is the GST purity attenuation specific to the family?",
       subtitle = sprintf("%d expression-matched random sets of %d genes, identical pipeline",
                          N_SETS, nrow(G))) +
  theme_bw(base_size = 10)
ggsave(file.path(OUT, "purity_negative_control.png"), p1,
       width = 7.5, height = 4.5, dpi = 300, bg = "white")

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Observed GST survival: %.2f%%", obs$pct_survive),
             sprintf("Random median: %.2f%% (95%% %.2f to %.2f)", q[2], q[1], q[3]),
             sprintf("Empirical p: %.3f", emp_p),
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_purity_control.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
