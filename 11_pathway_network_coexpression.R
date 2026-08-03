# =============================================================================
#  OBJECTIVE 1 · MODULE F — PATHWAY, NETWORK AND CO-EXPRESSION
#
#  WHY THE OBVIOUS ANALYSIS IS WORTHLESS
#    Running GO/KEGG enrichment on 22 GST genes returns "glutathione metabolism"
#    and "glutathione transferase activity". That is the definition of the
#    input, not a finding. Any reviewer will say so.
#
#  WHAT THIS DOES INSTEAD — three questions the data can actually answer
#
#    1. CO-EXPRESSION. What functional processes do GSTs track with across
#       tumours? Tested against a curated panel of ~90 genes spanning
#       oxidative stress, drug metabolism and transport, ER signalling,
#       proliferation and the basal/EMT programme. Non-circular because the
#       panel is chosen independently of the GST family.
#
#    2. THE NRF2 HYPOTHESIS. GST transcription is canonically driven by
#       NFE2L2 (NRF2) through antioxidant response elements. This scores each
#       tumour for NRF2 target activity and asks whether GST expression tracks
#       it - a directional test of a specific mechanism, not a list overlap.
#
#    3. PROTEIN NETWORK. STRING interaction network for the family plus its
#       first-shell partners, with enrichment computed on the EXPANDED set.
#       Including partners is what makes the enrichment informative rather
#       than tautological.
#
#  DATA: Toil hub (already working) plus STRING API. No large downloads.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleF_pathway_network")
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
# 1. THE COMPARISON PANEL
# =============================================================================
# Chosen a priori from the literature, not from anything in our data. This is
# what keeps the co-expression analysis non-circular.

PANEL <- list(
  `NRF2 targets` = c("NQO1","HMOX1","GCLC","GCLM","TXNRD1","SLC7A11","G6PD",
                     "PGD","ME1","SRXN1","PRDX1","GSR","FTL","FTH1","AKR1C1",
                     "AKR1B10","CBR1","EPHX1","ABCC1","ABCC2"),
  `NRF2 axis`    = c("NFE2L2","KEAP1","CUL3","MAF","MAFF","MAFG","MAFK"),
  `Drug transport` = c("ABCB1","ABCC3","ABCC4","ABCC5","ABCG2","SLC22A1",
                       "SLC47A1"),
  `Phase I metabolism` = c("CYP1A1","CYP1B1","CYP2E1","CYP3A4","CYP2D6",
                           "ALDH1A1","ALDH3A1","NAT1","NAT2"),
  `Phase II other` = c("UGT1A1","UGT1A6","SULT1A1","SULT1E1","COMT","TPMT",
                       "NQO2"),
  `Oxidative stress` = c("SOD1","SOD2","SOD3","CAT","GPX1","GPX2","GPX3",
                         "GPX4","PRDX2","PRDX3","PRDX6","TXN","TXN2"),
  `ER signalling` = c("ESR1","PGR","GATA3","FOXA1","XBP1","TFF1","GREB1",
                      "AR","SPDEF"),
  `Proliferation` = c("MKI67","CCNB1","AURKA","TOP2A","BIRC5","PLK1","BUB1",
                      "CDK1","MYBL2"),
  `Basal / EMT`  = c("KRT5","KRT14","KRT17","EGFR","VIM","CDH2","SNAI2",
                     "TWIST1","ZEB1","FN1"),
  `Apoptosis`    = c("BAX","BCL2","BCL2L1","CASP3","CASP8","CASP9","TP53",
                     "MAPK8","MAPK9")
)

panel_genes <- unique(unlist(PANEL))
cat("Comparison panel:", length(panel_genes), "genes across",
    length(PANEL), "processes\n")


# =============================================================================
# 2. FETCH  (batched, cached, resumable)
# =============================================================================
HOST <- "https://toil.xenahubs.net"
EXPR <- "TcgaTargetGtex_rsem_gene_tpm"
pf <- file.path(CACHE, "toil_panel.rds")

if (file.exists(pf)) {
  panel_raw <- readRDS(pf); message("Loaded panel from cache.")
} else {
  got <- list(); failed <- character(0)
  for (g in panel_genes) {
    r <- tryCatch(
      fetch_dense_values(host = HOST, dataset = EXPR, identifiers = g,
                         use_probeMap = TRUE, check = FALSE),
      error = function(e) NULL)
    if (is.null(r)) { failed <- c(failed, g); next }
    m <- as.matrix(r)
    if (!nrow(m)) { failed <- c(failed, g); next }
    got[[g]] <- if (nrow(m) > 1) colMeans(m, na.rm = TRUE) else m[1, ]
    Sys.sleep(0.2)
    if (length(got) %% 20 == 0) message("  ", length(got), " / ", length(panel_genes))
  }
  common <- Reduce(intersect, lapply(got, names))
  panel_raw <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(panel_raw) <- names(got)
  saveRDS(panel_raw, pf)
  message("Fetched ", nrow(panel_raw), " genes. Failed: ",
          if (length(failed)) paste(failed, collapse = ", ") else "none")
}

lin <- 2^panel_raw - 0.001; lin[lin < 0] <- 0
panel_mat <- log2(lin + 1)
panel_mat <- panel_mat[apply(panel_mat, 1, function(v) sd(v, na.rm = TRUE) > 0), ,
                       drop = FALSE]

gst_raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
gl <- 2^gst_raw - 0.001; gl[gl < 0] <- 0
gst_mat <- log2(gl + 1)


# =============================================================================
# 3. BREAST TUMOURS ONLY
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

s <- Reduce(intersect, list(clin$sample, colnames(panel_mat), colnames(gst_mat)))
clin <- clin[match(s, clin$sample), ]
P <- panel_mat[, s, drop = FALSE]
G <- gst_mat[intersect(GST, rownames(gst_mat)), s, drop = FALSE]

cat("\nTumours analysed:", length(s), "\n")
print(table(clin$group)); cat("\n")


# =============================================================================
# 4. CO-EXPRESSION — what do GSTs track with?
# =============================================================================
cormat <- matrix(NA_real_, nrow(G), nrow(P),
                 dimnames = list(rownames(G), rownames(P)))
pmat <- cormat
for (i in rownames(G)) for (j in rownames(P)) {
  x <- G[i, ]; y <- P[j, ]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 50 || sd(x[ok]) == 0 || sd(y[ok]) == 0) next
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  cormat[i, j] <- unname(ct$estimate); pmat[i, j] <- ct$p.value
}
fdr <- matrix(p.adjust(pmat, "BH"), nrow(pmat), dimnames = dimnames(pmat))

# Mean correlation of each GST with each process
proc <- sapply(names(PANEL), function(k) {
  gg <- intersect(PANEL[[k]], colnames(cormat))
  if (!length(gg)) return(rep(NA, nrow(cormat)))
  rowMeans(cormat[, gg, drop = FALSE], na.rm = TRUE)
})
proc <- proc[, colSums(is.na(proc)) < nrow(proc), drop = FALSE]

cat("===== MEAN CORRELATION OF EACH GST WITH EACH PROCESS =====\n")
print(round(proc, 3))
write.csv(round(proc, 3), file.path(OUT, "gst_process_correlation.csv"))

png(file.path(OUT, "coexpression_heatmap.png"), width = 2000, height = 1600, res = 220)
pheatmap(proc, cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-max(abs(proc), na.rm = TRUE), max(abs(proc), na.rm = TRUE),
                      length.out = 101),
         border_color = "grey85", fontsize = 10,
         main = "GST co-expression with functional processes\nSpearman rho, TCGA-BRCA tumours")
dev.off()

# Strongest individual partners
top <- which(!is.na(cormat) & abs(cormat) > 0.4 & fdr < 0.05, arr.ind = TRUE)
if (nrow(top)) {
  tp <- data.frame(GST = rownames(cormat)[top[, 1]],
                   partner = colnames(cormat)[top[, 2]],
                   rho = round(cormat[top], 3),
                   FDR = signif(fdr[top], 3)) %>% arrange(desc(abs(rho)))
  tp$process <- sapply(tp$partner, function(g) {
    h <- names(PANEL)[sapply(PANEL, function(v) g %in% v)]
    if (length(h)) h[1] else NA })
  cat("\n===== STRONGEST CO-EXPRESSION PARTNERS (|rho| > 0.4, FDR < 0.05) =====\n")
  print(head(as.data.frame(tp), 40))
  write.csv(tp, file.path(OUT, "top_coexpression_partners.csv"), row.names = FALSE)
}


# =============================================================================
# 5. THE NRF2 HYPOTHESIS
# =============================================================================
# GST transcription is canonically NRF2-driven via antioxidant response
# elements. If that holds here, GST expression should track a per-tumour NRF2
# target score. Scoring on TARGETS rather than on NFE2L2 mRNA matters: NRF2 is
# regulated post-translationally by KEAP1, so its own transcript level is a
# poor proxy for pathway activity.

nrf2_t <- intersect(PANEL$`NRF2 targets`, rownames(P))
cat("\nNRF2 target genes available:", length(nrf2_t), "\n")

z <- t(scale(t(P[nrf2_t, , drop = FALSE])))
nrf2_score <- colMeans(z, na.rm = TRUE)

nrf2_res <- lapply(rownames(G), function(g) {
  x <- G[g, ]; ok <- is.finite(x) & is.finite(nrf2_score)
  if (sum(ok) < 50 || sd(x[ok]) == 0) return(NULL)
  ct <- suppressWarnings(cor.test(x[ok], nrf2_score[ok], method = "spearman"))
  data.frame(gene = g, rho = round(unname(ct$estimate), 3), p = ct$p.value)
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "BH")) %>% arrange(desc(rho))

cat("\n===== GST EXPRESSION vs NRF2 TARGET SCORE =====\n")
print(as.data.frame(nrf2_res), digits = 3)
cat("\nPositive rho supports NRF2-driven transcription. Genes near zero are\n")
cat("regulated by something else - and the methylation result gives a\n")
cat("candidate for what.\n")
write.csv(nrf2_res, file.path(OUT, "nrf2_correlation.csv"), row.names = FALSE)

# Does the NRF2 score itself differ by subtype?
nd <- data.frame(sample = s, score = nrf2_score,
                 group = factor(clin$group, levels = SUB))
kw <- kruskal.test(score ~ group, data = nd)
cat(sprintf("\nNRF2 score by subtype: Kruskal-Wallis p = %.3g\n", kw$p.value))
print(nd %>% group_by(group) %>%
        summarise(median_score = round(median(score), 3), n = n(), .groups = "drop"))

p1 <- ggplot(nd, aes(group, score, fill = group)) +
  geom_violin(scale = "width", alpha = .3, colour = NA) +
  geom_boxplot(width = .28, outlier.size = .3, linewidth = .3) +
  scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                               LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
  labs(x = NULL, y = "NRF2 target activity score",
       title = "NRF2 pathway activity across PAM50 subtypes",
       subtitle = "Mean z-score of 20 canonical NRF2 target genes") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT, "nrf2_score_by_subtype.png"), p1, width = 6, height = 4.5, dpi = 300)


# =============================================================================
# 6. STRING PROTEIN NETWORK
# =============================================================================
# Enrichment is computed on the family PLUS its first-shell interaction
# partners. Enriching on the family alone would only recover "glutathione
# metabolism"; the partners are what make the result informative.

if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("STRINGdb", ask = FALSE, update = FALSE)
}

ok <- tryCatch({
  suppressPackageStartupMessages(library(STRINGdb))
  sdb <- STRINGdb$new(version = "12.0", species = 9606,
                      score_threshold = 700,   # high confidence only
                      input_directory = CACHE)

  mapped <- sdb$map(data.frame(gene = rownames(G), stringsAsFactors = FALSE),
                    "gene", removeUnmappedRows = TRUE)
  cat("\nMapped to STRING:", nrow(mapped), "of", nrow(G), "\n")

  png(file.path(OUT, "STRING_network_family.png"), width = 2000, height = 2000, res = 200)
  sdb$plot_network(mapped$STRING_id)
  dev.off()

  # First-shell partners
  nb <- sdb$get_neighbors(mapped$STRING_id)
  cat("First-shell interaction partners:", length(nb), "\n")

  expanded <- unique(c(mapped$STRING_id, nb))
  enr <- sdb$get_enrichment(expanded)

  if (!is.null(enr) && nrow(enr)) {
    enr <- enr %>% arrange(fdr)
    cat("\n===== STRING ENRICHMENT, FAMILY + PARTNERS =====\n")
    print(head(as.data.frame(enr[, intersect(c("category","term","description",
                                               "number_of_genes","fdr"),
                                             colnames(enr))]), 30))
    write.csv(enr, file.path(OUT, "STRING_enrichment_expanded.csv"), row.names = FALSE)

    top20 <- head(enr, 20)
    dcol <- intersect(c("description","term"), colnames(top20))[1]
    p2 <- ggplot(top20, aes(x = reorder(.data[[dcol]], -log10(fdr)),
                            y = -log10(fdr))) +
      geom_col(fill = "#2E86AB", width = .7) +
      coord_flip() +
      labs(x = NULL, y = expression(-log[10](FDR)),
           title = "Functional enrichment: GST family and first-shell partners",
           subtitle = "STRING v12, high-confidence interactions (score > 700)") +
      theme_bw(base_size = 10)
    ggsave(file.path(OUT, "STRING_enrichment.png"), p2, width = 9, height = 6, dpi = 300)
  }
  TRUE
}, error = function(e) { message("STRING step failed: ", conditionMessage(e)); FALSE })

if (!ok) {
  cat("\nSTRING via R failed. Manual alternative, two minutes:\n")
  cat("  string-db.org -> Multiple proteins -> paste the 16 gene symbols\n")
  cat("  -> Homo sapiens -> Search. Then Analysis tab for enrichment, and\n")
  cat("  '+ More' to add first-shell partners before re-running enrichment.\n")
}


# =============================================================================
# 7. GO AND KEGG ON THE EXPANDED SET
# =============================================================================
need <- c("clusterProfiler","org.Hs.eg.db")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}

ok2 <- tryCatch({
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db) })

  # Universe = family + strongest co-expression partners. Enriching the family
  # alone is circular; this asks what the GST-associated module does.
  partners <- if (exists("tp") && nrow(tp)) unique(tp$partner) else character(0)
  gene_set <- unique(c(rownames(G), partners))
  cat("\nEnrichment input:", length(gene_set), "genes (family +",
      length(partners), "co-expression partners)\n")

  eg <- bitr(gene_set, fromType = "SYMBOL", toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)

  go <- enrichGO(eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                 pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    cat("\n===== GO BIOLOGICAL PROCESS =====\n")
    print(head(as.data.frame(go)[, c("Description","GeneRatio","p.adjust")], 20))
    write.csv(as.data.frame(go), file.path(OUT, "GO_BP_enrichment.csv"), row.names = FALSE)
    ggsave(file.path(OUT, "GO_BP_dotplot.png"),
           dotplot(go, showCategory = 20) + theme_bw(base_size = 9),
           width = 9, height = 7, dpi = 300)
  }

  kk <- enrichKEGG(eg$ENTREZID, organism = "hsa", pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)
  if (!is.null(kk) && nrow(as.data.frame(kk))) {
    cat("\n===== KEGG PATHWAYS =====\n")
    print(head(as.data.frame(kk)[, c("Description","GeneRatio","p.adjust")], 20))
    write.csv(as.data.frame(kk), file.path(OUT, "KEGG_enrichment.csv"), row.names = FALSE)
    ggsave(file.path(OUT, "KEGG_dotplot.png"),
           dotplot(kk, showCategory = 20) + theme_bw(base_size = 9),
           width = 9, height = 7, dpi = 300)
  }
  TRUE
}, error = function(e) { message("Enrichment failed: ", conditionMessage(e)); FALSE })


# =============================================================================
# 8. HOW TO READ THIS
# =============================================================================
cat("\n============================================================\n")
cat("INTERPRETING THE OUTPUT\n")
cat("============================================================\n")
cat("1. Co-expression heatmap. Do the basal-programme genes (GSTP1, GSTA1,\n")
cat("   GSTA4) and the luminal-programme genes (GSTM2, GSTM3, GSTM4, GSTO2)\n")
cat("   separate into different clusters? If so, the two programmes are\n")
cat("   embedded in different biology, not merely differently methylated.\n\n")
cat("2. NRF2 correlation. A positive rho supports canonical NRF2-driven\n")
cat("   transcription. Genes near zero are regulated by something else -\n")
cat("   and your methylation data supplies a candidate mechanism.\n\n")
cat("3. If the NRF2 score is highest in Basal, that is a second, independent\n")
cat("   explanation for GSTP1 retention there, complementary to the escape\n")
cat("   from promoter methylation.\n\n")
cat("4. STRING and GO/KEGG on the expanded set. Terms beyond glutathione\n")
cat("   metabolism are the informative ones - they show what the GST module\n")
cat("   is connected to. Glutathione terms recovering is a positive control.\n")
cat("============================================================\n")

writeLines(c(
  paste("Run:", Sys.time()),
  "Co-expression tested against an a priori curated panel of 10 processes,",
  "chosen independently of the GST family to avoid circularity.",
  "NRF2 activity scored on 20 canonical target genes rather than NFE2L2 mRNA,",
  "because NRF2 is regulated post-translationally by KEAP1.",
  "STRING and GO/KEGG enrichment computed on the family PLUS first-shell",
  "interaction partners; enrichment on the family alone would be tautological.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleF.txt"))

message("\nWritten to ", normalizePath(OUT))
