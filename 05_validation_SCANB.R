# =============================================================================
#  GSE96058 (SCAN-B) — PROCESS + ANALYSE + THREE-COHORT CONCORDANCE
#
#  Self-contained. Needs only:
#    cache/GSE96058/...csv.gz          (complete, 591.7 MB — you have this)
#    cache/gse96058_pheno.rds          (already cached)
#    cache/toil_gst.rds                (TCGA, already cached)
#    cache/metabric_gst.rds            (METABRIC, already cached)
#
#  PRODUCES
#    scanb_subtype_medians.csv      median TPM per gene per subtype
#    scanb_subtype_stats.csv        Kruskal-Wallis + epsilon-squared
#    three_cohort_concordance.csv   TCGA vs METABRIC vs SCAN-B, per gene
#    tcga_vs_scanb_absolute.csv     absolute TPM agreement, RNA-seq vs RNA-seq
#    Panel_SCANB.png/.pdf           all genes across subtypes
#    scanb_survival_input.csv       expression + OS, ready for Cox
#
#  DISK: decompresses to a few GB, deleted at the end.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "figures", "scanb")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("data.table","R.utils","dplyr","tidyr","tibble","ggplot2")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages(lapply(pk, library, character.only = TRUE))

GST <- c("GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTM1","GSTM2","GSTM3",
         "GSTM4","GSTM5","GSTP1","GSTT1","GSTT2","GSTT2B","GSTT4",
         "GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. PROCESS THE COMPLETED DOWNLOAD
# =============================================================================
ex_f <- file.path(CACHE, "gse96058_gst.rds")

if (file.exists(ex_f)) {
  obj <- readRDS(ex_f); mat <- obj$mat
  message("Loaded SCAN-B from cache.")
} else {

  supp  <- file.path(CACHE, "GSE96058")
  gz    <- list.files(supp, pattern = "gene_expression.*csv\\.gz$", full.names = TRUE)[1]
  if (is.na(gz)) stop("Expression .gz not found in ", supp)
  plain <- sub("\\.gz$", "", gz)

  cat(sprintf("Compressed file: %.1f MB\n", file.size(gz) / 1e6))
  if (file.size(gz) < 5.8e8)
    stop("The .gz is only ", round(file.size(gz)/1e6, 1),
         " MB. Expected ~592 MB. Re-run the resume loop.")

  if (file.exists(plain)) unlink(plain)
  message("Decompressing...")
  gunzip(gz, destname = plain, remove = FALSE, overwrite = TRUE)
  cat(sprintf("Decompressed: %.2f GB\n", file.size(plain) / 1e9))

  # ---- integrity: the last line must have as many fields as the header ----
  hdr <- readLines(plain, n = 1)
  n_hdr <- length(strsplit(hdr, ",")[[1]])
  lastl <- tail(readLines(plain, warn = FALSE), 1)
  n_last <- length(strsplit(lastl, ",")[[1]])
  cat(sprintf("Header fields: %d | last line fields: %d\n", n_hdr, n_last))
  if (abs(n_last - n_hdr) > 1)
    stop("Last line still truncated. The file is incomplete.")
  message("Integrity check passed.")

  message("Reading (several minutes)...")
  dt <- fread(plain, showProgress = TRUE)
  genes <- dt[[1]]
  message("Rows read: ", length(genes))
  if (length(genes) < 15000)
    stop("Only ", length(genes), " rows. Still incomplete.")

  found <- intersect(GST, genes)
  cat("\nGST genes present:", length(found), "of", length(GST), "\n")
  cat(paste(sort(found), collapse = ", "), "\n")
  cat("Absent:", paste(setdiff(GST, found), collapse = ", "), "\n\n")
  if (!length(found)) stop("No GST genes found. Check ID format: ",
                           paste(head(genes, 5), collapse = ", "))

  m <- as.matrix(dt[, -1, with = FALSE]); rownames(m) <- genes
  rm(dt); gc()

  # ---- confirm the transform before inverting it --------------------------
  # log2(FPKM + 0.1) has a floor at log2(0.1) = -3.3219 when FPKM = 0.
  obs_min <- min(m, na.rm = TRUE)
  cat(sprintf("Observed minimum: %.4f  (log2(0.1) = %.4f)\n", obs_min, log2(0.1)))
  OFFSET <- 0.1
  if (abs(obs_min - log2(OFFSET)) > 0.05)
    warning("Minimum does not match log2(0.1). Verify the GEO transform ",
            "before trusting absolute TPM. Rank results remain valid.")

  # ---- per-sample FPKM totals, chunked to keep memory flat ----------------
  message("Computing per-sample FPKM totals...")
  tot <- numeric(ncol(m)); step <- 2000
  for (i in seq(1, nrow(m), by = step)) {
    j <- min(i + step - 1, nrow(m))
    blk <- 2^m[i:j, , drop = FALSE] - OFFSET
    blk[blk < 0] <- 0
    tot <- tot + colSums(blk, na.rm = TRUE)
    rm(blk)
  }
  gc()
  cat(sprintf("Total FPKM per sample: median %.0f (range %.0f-%.0f)\n",
              median(tot), min(tot), max(tot)))

  sub_m <- m[found, , drop = FALSE]; rm(m); gc()
  fpkm <- 2^sub_m - OFFSET; fpkm[fpkm < 0] <- 0
  tpm  <- t(t(fpkm) / tot) * 1e6
  mat  <- log2(tpm + 1)

  saveRDS(list(mat = mat, unit = "log2(TPM+1)"), ex_f)
  message("Cached to ", ex_f)
  unlink(plain)
  message("Removed the decompressed csv.")
}


# =============================================================================
# 2. MATCH TO PHENOTYPE
# =============================================================================
ph <- readRDS(file.path(CACHE, "gse96058_pheno.rds"))
ph$subtype <- c(basal = "Basal", her2 = "Her2", luma = "LumA",
                lumb = "LumB", normal = NA)[tolower(ph$pam50)]

idx <- match(colnames(mat), ph$title)
if (sum(!is.na(idx)) < 0.5 * ncol(mat)) idx <- match(colnames(mat), ph$gsm)
if (sum(!is.na(idx)) < 0.5 * ncol(mat))
  stop("Could not match expression columns to phenotype rows.")

ann <- ph[idx, ]
ann$sample <- colnames(mat)
ann <- ann[!duplicated(ann$title), ]      # drop the 136 technical replicates
ann <- ann[!is.na(ann$subtype), ]
mat <- mat[, ann$sample, drop = FALSE]
ann$subtype <- factor(ann$subtype, levels = SUB)

cat("\n===== SCAN-B ANALYSIS SET =====\n"); print(table(ann$subtype)); cat("\n")

# Drop genes with no measurable variance. GSTT1 and GSTT2 are not quantified in
# this dataset and would otherwise render as flat zeros on a 1e-16 axis.
keep_var <- apply(mat, 1, function(v) {
  v <- v[is.finite(v)]
  length(v) > 0 && sd(v) > 1e-10
})
if (any(!keep_var)) {
  message("Dropping genes with no variance in SCAN-B: ",
          paste(rownames(mat)[!keep_var], collapse = ", "))
  mat <- mat[keep_var, , drop = FALSE]
}

df <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(ann[, c("sample","subtype")], by = "sample") %>%
  pivot_longer(-c(sample, subtype), names_to = "gene", values_to = "expr") %>%
  filter(!is.na(subtype), is.finite(expr))


# =============================================================================
# 3. TABLES
# =============================================================================
sb_wide <- df %>% group_by(gene, subtype) %>%
  summarise(TPM = round(2^median(expr) - 1, 3), .groups = "drop") %>%
  pivot_wider(names_from = subtype, values_from = TPM) %>%
  mutate(max_TPM = apply(across(all_of(SUB)), 1, max, na.rm = TRUE),
         highest = SUB[apply(across(all_of(SUB)), 1, which.max)],
         lowest  = SUB[apply(across(all_of(SUB)), 1, which.min)],
         status  = ifelse(max_TPM >= 1, "detected", "not detected")) %>%
  arrange(desc(max_TPM))

cat("===== SCAN-B MEDIAN TPM BY SUBTYPE =====\n")
print(as.data.frame(sb_wide), digits = 3)
write.csv(sb_wide, file.path(OUT, "scanb_subtype_medians.csv"), row.names = FALSE)

sb_stats <- df %>% group_by(gene) %>%
  summarise(n = n(),
            H = unname(kruskal.test(expr ~ droplevels(subtype))$statistic),
            p = kruskal.test(expr ~ droplevels(subtype))$p.value, .groups = "drop") %>%
  mutate(eps2 = round(H / ((n^2 - 1) / (n + 1)), 4),
         FDR  = p.adjust(p, "BH"),
         effect = cut(eps2, c(-Inf,.01,.08,.26,Inf),
                      labels = c("negligible","small","medium","large"))) %>%
  arrange(desc(eps2))

cat("\n===== SUBTYPE STRUCTURE IN SCAN-B =====\n")
print(as.data.frame(sb_stats), digits = 3)
write.csv(sb_stats, file.path(OUT, "scanb_subtype_stats.csv"), row.names = FALSE)


# =============================================================================
# 4. THREE-COHORT CONCORDANCE — the payoff
# =============================================================================
tc_f <- file.path(WORKDIR, "figures", "withintumour_toil",
                  "subtype_medians_wide_toil.csv")
mb_f <- file.path(WORKDIR, "figures", "metabric", "metabric_subtype_medians.csv")

if (file.exists(tc_f) && file.exists(mb_f)) {

  tc <- read.csv(tc_f, stringsAsFactors = FALSE)
  mb <- read.csv(mb_f, stringsAsFactors = FALSE)
  common <- Reduce(intersect, list(tc$gene, mb$gene, sb_wide$gene))
  message("\nMeasurable in all three cohorts: ", length(common), " genes")

  # ---- Test 1: does the expression hierarchy replicate? -------------------
  a <- tc$max_TPM[match(common, tc$gene)]
  b <- mb$overall[match(common, mb$gene)]
  s <- sb_wide$max_TPM[match(common, sb_wide$gene)]

  r_ts <- suppressWarnings(cor.test(a, s, method = "spearman"))
  r_tm <- suppressWarnings(cor.test(a, b, method = "spearman"))
  r_ms <- suppressWarnings(cor.test(b, s, method = "spearman"))

  cat("\n===== TEST 1: FAMILY EXPRESSION HIERARCHY =====\n")
  cat(sprintf("TCGA vs SCAN-B     rho = %.3f  (both RNA-seq)\n", unname(r_ts$estimate)))
  cat(sprintf("TCGA vs METABRIC   rho = %.3f\n", unname(r_tm$estimate)))
  cat(sprintf("METABRIC vs SCAN-B rho = %.3f\n", unname(r_ms$estimate)))

  # ---- Test 1b: ABSOLUTE agreement, only possible RNA-seq vs RNA-seq ------
  cat("\n===== TEST 1b: ABSOLUTE TPM AGREEMENT (TCGA vs SCAN-B) =====\n")
  cmp <- data.frame(gene = common,
                    TCGA_TPM  = round(a, 2),
                    SCANB_TPM = round(s, 2),
                    ratio     = round(s / pmax(a, 1e-6), 2)) %>%
    arrange(desc(TCGA_TPM))
  print(cmp, row.names = FALSE)
  cat(sprintf("\nMedian ratio %.2f | %d of %d genes within 2-fold\n",
              median(cmp$ratio, na.rm = TRUE),
              sum(cmp$ratio > 0.5 & cmp$ratio < 2, na.rm = TRUE), nrow(cmp)))
  cat("METABRIC cannot support this test - microarray has no TPM. Two RNA-seq\n")
  cat("cohorts agreeing on absolute values is the strongest public-data claim.\n")
  write.csv(cmp, file.path(OUT, "tcga_vs_scanb_absolute.csv"), row.names = FALSE)

  # ---- Test 2: per-gene subtype ordering across all three -----------------
  conc <- lapply(common, function(g) {
    x <- unlist(tc[tc$gene == g, SUB])
    y <- unlist(mb[mb$gene == g, SUB])
    z <- unlist(sb_wide[sb_wide$gene == g, SUB])
    data.frame(gene = g,
               TCGA_high  = SUB[which.max(x)], TCGA_low  = SUB[which.min(x)],
               MB_high    = SUB[which.max(y)], MB_low    = SUB[which.min(y)],
               SCANB_high = SUB[which.max(z)], SCANB_low = SUB[which.min(z)],
               rho_TCGA_SCANB = round(cor(rank(x), rank(z), method = "spearman"), 3))
  }) %>% bind_rows() %>%
    mutate(all_three_high = TCGA_high == MB_high & TCGA_high == SCANB_high,
           all_three_low  = TCGA_low  == MB_low  & TCGA_low  == SCANB_low,
           verdict = case_when(
             all_three_high & all_three_low ~ "replicated in all 3",
             all_three_high | all_three_low ~ "partial",
             TCGA_high == SCANB_high        ~ "RNA-seq cohorts agree",
             TRUE ~ "discordant")) %>%
    arrange(desc(rho_TCGA_SCANB))

  cat("\n===== TEST 2: SUBTYPE ORDERING ACROSS THREE COHORTS =====\n")
  print(as.data.frame(conc[, c("gene","TCGA_high","MB_high","SCANB_high",
                               "TCGA_low","MB_low","SCANB_low",
                               "rho_TCGA_SCANB","verdict")]))
  write.csv(conc, file.path(OUT, "three_cohort_concordance.csv"), row.names = FALSE)

  cat("\n"); print(table(conc$verdict))

  for (g in c("GSTA1","GSTM5","GSTP1","GSTM2")) {
    r <- conc[conc$gene == g, ]
    if (nrow(r))
      cat(sprintf("\n%-7s highest: TCGA %s | METABRIC %s | SCAN-B %s  ->  %s",
                  g, r$TCGA_high, r$MB_high, r$SCANB_high, r$verdict))
  }
  cat("\n")

} else {
  message("\nTCGA or METABRIC tables not found - skipping concordance.")
}


# =============================================================================
# 5. FIGURE
# =============================================================================
df$gene <- factor(df$gene, levels = sb_stats$gene)
n_tab <- df %>% distinct(sample, subtype) %>% count(subtype) %>%
  mutate(lab = paste0(subtype, "\n(n=", n, ")"))
df$glab <- factor(n_tab$lab[match(df$subtype, n_tab$subtype)], levels = n_tab$lab)
PAL <- setNames(c("#C0392B","#8E7CC3","#2E86AB","#E8A33D")[
  match(n_tab$subtype, SUB)], n_tab$lab)

p <- ggplot(df, aes(glab, expr, fill = glab)) +
  geom_violin(scale = "width", width = .85, alpha = .3, colour = NA, linewidth = .2) +
  geom_boxplot(width = .28, outlier.size = .2, outlier.alpha = .25, linewidth = .28) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = expression(log[2](TPM + 1)),
       title = "GST family across PAM50 subtypes in SCAN-B (GSE96058)",
       subtitle = paste0("Independent Swedish cohort, RNA-seq, n = ",
                         nrow(ann),
                         "; FPKM converted to TPM across the full transcriptome")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.3),
        strip.background = element_rect(fill = "grey93", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "Panel_SCANB.png"), p, width = 10.5, height = 13, dpi = 300)
ggsave(file.path(OUT, "Panel_SCANB.pdf"), p, width = 10.5, height = 13)


# =============================================================================
# 6. SURVIVAL TABLE  (prepared, not analysed here)
# =============================================================================
surv <- as.data.frame(t(mat)) %>% rownames_to_column("sample") %>%
  left_join(ann[, c("sample","subtype","os_days","os_event")], by = "sample") %>%
  filter(!is.na(os_days), !is.na(os_event))

write.csv(surv, file.path(OUT, "scanb_survival_input.csv"), row.names = FALSE)
cat(sprintf("\nSurvival table: %d patients, %d events (%.1f%%)\n",
            nrow(surv), sum(surv$os_event), 100 * mean(surv$os_event)))

writeLines(c(
  paste("Run:", Sys.time()),
  "SCAN-B GSE96058. log2(FPKM+0.1) as published, back-transformed and",
  "converted to TPM across the full transcriptome, then log2(TPM+1).",
  paste("Samples analysed:", nrow(ann), "(technical replicates removed)"),
  "", capture.output(print(table(ann$subtype))),
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_scanb.txt"))

message("\nWritten to ", normalizePath(OUT))
