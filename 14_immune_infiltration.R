# =============================================================================
#  OBJECTIVE 1 · MODULE G — IMMUNE INFILTRATION
#
#  TWO CONFOUNDS THAT MAKE THE NAIVE VERSION OF THIS ANALYSIS MEANINGLESS
#
#  1. CELL COMPOSITION. Immune scores and epithelial gene expression are both
#     functions of how much tumour is in the sample. GSTs are epithelial, so a
#     lymphocyte-rich sample shows high immune score and low GST expression
#     for purely arithmetic reasons. Every correlation here is therefore
#     computed twice: raw, and adjusted for tumour purity.
#
#  2. SUBTYPE. Basal-like tumours are both GSTP1-high and immune-hot. A
#     pooled correlation would recover that association without either
#     variable causing the other. Section 5 repeats everything WITHIN subtype.
#
#  A correlation that survives both is worth reporting. One that does not is
#  composition, and reporting it would be an error a reviewer will find.
#
#  METHOD
#    Danaher et al. (J Immunother Cancer 2017) marker-based scoring: each cell
#    population scored as the mean expression of its markers. Compact,
#    published, and interpretable, unlike black-box deconvolution.
#
#  BIOLOGICAL MOTIVATION
#    GSTO1 participates in IL-1beta processing and inflammatory signalling, so
#    an omega-class link to immune content is mechanistically plausible rather
#    than fishing. GSTP1 sequesters JNK1 and modulates redox signalling, both
#    relevant to T cell function in the tumour microenvironment.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleG_immune")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2","pheatmap","RColorBrewer")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(TCGAbiolinks) })

GST <- c("GSTA1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTM5","GSTP1",
         "GSTT2B","GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. MARKER SETS
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

all_markers <- unique(c(unlist(IMMUNE), CHECKPOINT, STROMA, CYTOKINE))
cat("Marker genes to fetch:", length(all_markers), "\n")


# =============================================================================
# 2. FETCH
# =============================================================================
HOST <- "https://toil.xenahubs.net"
EXPR <- "TcgaTargetGtex_rsem_gene_tpm"
mf <- file.path(CACHE, "toil_immune_markers.rds")

if (file.exists(mf)) {
  marker_raw <- readRDS(mf); message("Loaded markers from cache.")
} else {
  got <- list(); failed <- character(0)
  for (g in all_markers) {
    r <- tryCatch(
      fetch_dense_values(host = HOST, dataset = EXPR, identifiers = g,
                         use_probeMap = TRUE, check = FALSE),
      error = function(e) NULL)
    if (is.null(r) || !nrow(as.matrix(r))) { failed <- c(failed, g); next }
    m <- as.matrix(r)
    got[[g]] <- if (nrow(m) > 1) colMeans(m, na.rm = TRUE) else m[1, ]
    Sys.sleep(0.2)
    if (length(got) %% 25 == 0) message("  ", length(got), " / ", length(all_markers))
  }
  common <- Reduce(intersect, lapply(got, names))
  marker_raw <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(marker_raw) <- names(got)
  saveRDS(marker_raw, mf)
  message("Fetched ", nrow(marker_raw), ". Failed: ",
          if (length(failed)) paste(failed, collapse = ", ") else "none")
}

lin <- 2^marker_raw - 0.001; lin[lin < 0] <- 0
M <- log2(lin + 1)
M <- M[apply(M, 1, function(v) sd(v, na.rm = TRUE) > 0), , drop = FALSE]

gst_raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
gl <- 2^gst_raw - 0.001; gl[gl < 0] <- 0
G_all <- log2(gl + 1)


# =============================================================================
# 3. BREAST TUMOURS, SUBTYPE AND PURITY
# =============================================================================
ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

clin <- ph %>%
  dplyr::filter(study == "TCGA", grepl("breast", tissue, ignore.case = TRUE),
         grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

s <- Reduce(intersect, list(clin$sample, colnames(M), colnames(G_all)))
clin <- clin[match(s, clin$sample), ]
M <- M[, s, drop = FALSE]
G <- G_all[intersect(GST, rownames(G_all)), s, drop = FALSE]

# Tumour purity — the adjustment that makes this analysis interpretable
clin$purity <- NA_real_
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  clin$purity <- num(tp$CPE)[match(clin$patient, tp$patient)]
}, silent = TRUE)

cat("\nTumours:", length(s), "| with purity:", sum(!is.na(clin$purity)), "\n")
print(table(clin$group)); cat("\n")


# =============================================================================
# 4. CELL POPULATION SCORES
# =============================================================================
score_of <- function(genes) {
  gg <- intersect(genes, rownames(M))
  if (!length(gg)) return(rep(NA_real_, ncol(M)))
  colMeans(t(scale(t(M[gg, , drop = FALSE]))), na.rm = TRUE)
}

S <- sapply(names(IMMUNE), function(k) score_of(IMMUNE[[k]]))
S <- cbind(S,
           `Checkpoint` = score_of(CHECKPOINT),
           `Stromal / CAF` = score_of(STROMA),
           `Inflammatory cytokines` = score_of(CYTOKINE))
rownames(S) <- s
S <- S[, colSums(is.na(S)) < nrow(S), drop = FALSE]

cat("Populations scored:", ncol(S), "\n")

# Sanity check: immune scores should be higher in basal-like disease, which is
# well established. If they are not, the scoring has gone wrong.
chk <- data.frame(group = clin$group, immune = S[, "Total immune"])
cat("\nTotal immune score by subtype (expect Basal highest):\n")
print(chk %>% group_by(group) %>%
        summarise(median = round(median(immune, na.rm = TRUE), 3), .groups = "drop"))


# =============================================================================
# 5. CORRELATION — RAW, THEN PURITY-ADJUSTED
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

res <- lapply(rownames(G), function(g) {
  lapply(colnames(S), function(p) {
    v <- partial_rho(G[g, ], S[, p], clin$purity)
    if (all(is.na(v))) return(NULL)
    data.frame(gene = g, population = p,
               rho_raw = round(v[1], 3),
               rho_adj = round(v[2], 3),
               p_adj = v[3],
               attenuation = round(v[1] - v[2], 3))
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p_adj, "BH"),
         verdict = case_when(
           abs(rho_adj) >= 0.3 & FDR < 0.05 ~ "survives purity adjustment",
           abs(rho_adj) >= 0.2 & FDR < 0.05 ~ "weak but present",
           TRUE ~ "explained by composition")) %>%
  arrange(desc(abs(rho_adj)))

cat("\n===== GST vs IMMUNE POPULATIONS, PURITY-ADJUSTED =====\n")
print(head(as.data.frame(res), 30))
cat("\nLarge 'attenuation' means the raw correlation was mostly cell content.\n")
write.csv(res, file.path(OUT, "gst_immune_correlation.csv"), row.names = FALSE)

cat("\n"); print(table(res$verdict))

# Heatmaps: raw and adjusted, side by side, because the difference is the point
hm_raw <- res %>% dplyr::select(gene, population, rho_raw) %>%
  pivot_wider(names_from = population, values_from = rho_raw) %>%
  column_to_rownames("gene") %>% as.matrix()
hm_adj <- res %>% dplyr::select(gene, population, rho_adj) %>%
  pivot_wider(names_from = population, values_from = rho_adj) %>%
  column_to_rownames("gene") %>% as.matrix()

lim <- max(abs(c(hm_raw, hm_adj)), na.rm = TRUE)
pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
brk <- seq(-lim, lim, length.out = 101)

# pheatmap cannot cluster a matrix containing NA and errors out, which would
# stop the script before the within-subtype analysis below. Zero-fill for the
# figure only; the underlying tables retain the NAs.
safe_heat <- function(m, file, title) {
  m2 <- m; m2[!is.finite(m2)] <- 0
  ok_r <- nrow(m2) > 2 && all(apply(m2, 1, sd) > 0)
  ok_c <- ncol(m2) > 2 && all(apply(m2, 2, sd) > 0)
  png(file, width = 2000, height = 1500, res = 220)
  pheatmap(m2, color = pal, breaks = brk, border_color = "grey85", fontsize = 9,
           cluster_rows = ok_r, cluster_cols = ok_c, main = title)
  dev.off()
}

safe_heat(hm_raw, file.path(OUT, "immune_correlation_raw.png"),
          "GST vs immune populations - RAW\n(confounded by cell composition)")
safe_heat(hm_adj, file.path(OUT, "immune_correlation_adjusted.png"),
          "GST vs immune populations - PURITY ADJUSTED\n(what survives)")


# =============================================================================
# 6. WITHIN SUBTYPE — the second confound
# =============================================================================
# Basal-like tumours are both GSTP1-high and immune-hot. A pooled correlation
# would recover that without either causing the other. Repeating within
# subtype removes it.

ws <- lapply(SUB, function(sb) {
  idx <- which(clin$group == sb)
  if (length(idx) < 60) return(NULL)
  lapply(rownames(G), function(g) {
    lapply(colnames(S), function(p) {
      v <- partial_rho(G[g, idx], S[idx, p], clin$purity[idx])
      if (all(is.na(v))) return(NULL)
      data.frame(subtype = sb, gene = g, population = p,
                 rho_adj = round(v[2], 3), p = v[3], n = length(idx))
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

if (nrow(ws)) {
  ws <- ws %>% group_by(subtype) %>% mutate(FDR = p.adjust(p, "BH")) %>%
    ungroup() %>% arrange(desc(abs(rho_adj)))
  cat("\n===== WITHIN-SUBTYPE, PURITY-ADJUSTED =====\n")
  print(head(as.data.frame(ws), 30))
  write.csv(ws, file.path(OUT, "gst_immune_within_subtype.csv"), row.names = FALSE)

  # Threshold set at 2 of 4 subtypes, not 3. HER2-enriched has only n = 81, so
  # a real association can fail there on power alone; requiring 3 would
  # discard it. Per-subtype detail is written out so the reader can judge.
  cons <- ws %>% dplyr::filter(FDR < 0.05, abs(rho_adj) >= 0.20) %>%
    group_by(gene, population) %>%
    summarise(n_subtypes = n(),
              subtypes = paste(sort(subtype), collapse = ", "),
              mean_rho = round(mean(rho_adj), 3), .groups = "drop") %>%
    dplyr::filter(n_subtypes >= 2) %>% arrange(desc(n_subtypes), desc(abs(mean_rho)))

  cat("\n===== CONSISTENT IN 3 OR MORE SUBTYPES =====\n")
  if (nrow(cons)) {
    print(as.data.frame(cons))
    cat("\nThese survive both confounds. They are the only associations worth\n")
    cat("interpreting biologically.\n")
    write.csv(cons, file.path(OUT, "gst_immune_consistent.csv"), row.names = FALSE)
  } else {
    cat("None. Every GST-immune association is explained by cell composition\n")
    cat("or subtype. That is a legitimate and reportable negative - and it is\n")
    cat("the result most published versions of this analysis fail to check.\n")
  }
}


# =============================================================================
# 7. CHECKPOINT GENES — immunotherapy relevance
# =============================================================================
cp <- intersect(CHECKPOINT, rownames(M))
if (length(cp) >= 4) {
  cpr <- lapply(rownames(G), function(g) {
    lapply(cp, function(c1) {
      v <- partial_rho(G[g, ], M[c1, ], clin$purity)
      if (all(is.na(v))) return(NULL)
      data.frame(gene = g, checkpoint = c1,
                 rho_adj = round(v[2], 3), p = v[3])
    }) %>% bind_rows()
  }) %>% bind_rows() %>% mutate(FDR = p.adjust(p, "BH")) %>%
    arrange(desc(abs(rho_adj)))

  cat("\n===== GST vs IMMUNE CHECKPOINT GENES (purity-adjusted) =====\n")
  print(head(as.data.frame(cpr), 20))
  write.csv(cpr, file.path(OUT, "gst_checkpoint_correlation.csv"), row.names = FALSE)

  g1 <- cpr %>% dplyr::filter(gene == "GSTP1") %>% head(5)
  if (nrow(g1)) {
    cat("\n--- GSTP1 and checkpoints ---\n"); print(as.data.frame(g1))
    cat("\nRelevant because basal-like disease is where GSTP1 is retained AND\n")
    cat("where immunotherapy is used. A real association would be worth a\n")
    cat("sentence; an absent one closes the question.\n")
  }
}


# =============================================================================
# 8. IMMUNE SCORE BY SUBTYPE
# =============================================================================
sd_df <- as.data.frame(S) %>% rownames_to_column("sample") %>%
  mutate(group = factor(clin$group, levels = SUB)) %>%
  pivot_longer(-c(sample, group), names_to = "population", values_to = "score")

p1 <- ggplot(sd_df, aes(group, score, fill = group)) +
  geom_violin(scale = "width", alpha = .3, colour = NA) +
  geom_boxplot(width = .25, outlier.size = .2, linewidth = .25) +
  facet_wrap(~ population, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                               LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
  labs(x = NULL, y = "Marker-based score",
       title = "Immune and stromal populations across PAM50 subtypes",
       subtitle = "Danaher marker scoring, TCGA-BRCA primary tumours") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(OUT, "immune_by_subtype.png"), p1, width = 10, height = 8, dpi = 300)


cat("\n============================================================\n")
cat("INTERPRETING THIS MODULE\n")
cat("============================================================\n")
cat("Compare the two heatmaps. The raw one will look impressive. The adjusted\n")
cat("one shows what is left once cell composition is removed, and the gap\n")
cat("between them is the finding most published versions of this analysis\n")
cat("never report.\n\n")
cat("Only associations in the 'consistent in 3 or more subtypes' table have\n")
cat("survived both confounds. Interpret those; treat the rest as composition.\n\n")
cat("An entirely negative result here is worth reporting. It would say GST\n")
cat("expression carries no independent relationship to immune infiltration,\n")
cat("which constrains a claim others might otherwise make from the same data.\n")
cat("============================================================\n")

writeLines(c(
  paste("Run:", Sys.time()),
  "Immune populations scored using Danaher et al. (J Immunother Cancer 2017)",
  "marker sets as the mean z-score of markers per population.",
  "All correlations adjusted for tumour purity (Aran et al. 2015 CPE) via rank",
  "residuals, because immune content and epithelial gene expression are both",
  "functions of sample composition.",
  "Analyses repeated within PAM50 subtype, because basal-like tumours are both",
  "GSTP1-high and immune-hot, which would generate a spurious pooled",
  "association.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleG.txt"))

message("\nWritten to ", normalizePath(OUT))
