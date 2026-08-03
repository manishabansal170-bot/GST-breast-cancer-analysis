# =============================================================================
#  OBJECTIVE 1 · MODULE H — miRNA REGULATION
#
#  WHY THE OBVIOUS VERSION FAILS
#    Correlating all ~1,900 miRNAs against 16 GST genes tests 30,000 pairs.
#    After correction almost nothing survives, and whatever does is as likely
#    to be chance as biology. Worse, a negative correlation between any two
#    genes in bulk tumour tissue is weak evidence of anything.
#
#  THE DESIGN THAT WORKS
#    1. Retrieve miRNAs PREDICTED or VALIDATED to target GST genes
#       (multiMiR: miRTarBase, TarBase, TargetScan, miRDB, DIANA).
#    2. Test only those pairs. This cuts the multiple-testing burden by two
#       orders of magnitude and every test has a prior reason to exist.
#    3. Require NEGATIVE correlation, since miRNAs repress their targets.
#    4. Adjust for tumour purity - miRNA and mRNA abundance are both
#       composition-sensitive, as Module G demonstrated emphatically.
#    5. Check within subtype.
#
#  WHAT THIS CAN AND CANNOT SHOW
#    A predicted target that is anticorrelated in 1,000 tumours is a credible
#    candidate. It is not proof of targeting - that needs a luciferase assay
#    or miRNA overexpression. State it as candidate regulation, nothing more.
#
#  RELEVANCE TO YOUR MAIN FINDING
#    You have shown promoter methylation explains 28-43% of expression
#    variance for the programme genes. Most of the variance is unexplained.
#    miRNA regulation is a plausible second layer, and if a miRNA tracks the
#    subtype split it becomes part of the mechanism rather than a footnote.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleH_mirna")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)
options(timeout = 1000000, download.file.method = "libcurl")

pk <- c("dplyr","tidyr","tibble","ggplot2","pheatmap","RColorBrewer")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("multiMiR", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(UCSCXenaTools); library(TCGAbiolinks); library(multiMiR) })

GST <- c("GSTA1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTM5","GSTP1",
         "GSTT2B","GSTO1","GSTO2","GSTZ1","GSTK1","MGST1","MGST2","MGST3")
SUB <- c("Basal","Her2","LumA","LumB")


# =============================================================================
# 1. WHICH miRNAs TARGET GST GENES
# =============================================================================
tf <- file.path(CACHE, "mirna_targets_gst.rds")

if (file.exists(tf)) {
  targets <- readRDS(tf); message("Loaded target predictions from cache.")
} else {
  message("Querying multiMiR - this takes a few minutes...")

  val <- tryCatch({
    r <- get_multimir(org = "hsa", target = GST, table = "validated",
                      summary = FALSE)
    as.data.frame(r@data)
  }, error = function(e) { message("  validated query failed: ",
                                   conditionMessage(e)); NULL })

  # Top 10% of predictions only. The full set runs to tens of thousands and is
  # dominated by weak seed matches.
  pre <- tryCatch({
    r <- get_multimir(org = "hsa", target = GST, table = "predicted",
                      summary = FALSE, predicted.cutoff.type = "p",
                      predicted.cutoff = 10)
    as.data.frame(r@data)
  }, error = function(e) { message("  predicted query failed: ",
                                   conditionMessage(e)); NULL })

  targets <- bind_rows(
    if (!is.null(val) && nrow(val)) val %>% mutate(evidence = "validated") else NULL,
    if (!is.null(pre) && nrow(pre)) pre %>% mutate(evidence = "predicted") else NULL)

  if (!nrow(targets))
    stop("multiMiR returned nothing. The service may be down - retry later, ",
         "or use ENCORI/starBase manually at rnasysu.com/encori")

  saveRDS(targets, tf)
}

mcol <- grep("mature_mirna_id|mature_mirna", colnames(targets), value = TRUE)[1]
gcol <- grep("target_symbol", colnames(targets), value = TRUE)[1]

pairs <- targets %>%
  transmute(mirna = .data[[mcol]], gene = .data[[gcol]], evidence, database) %>%
  dplyr::filter(!is.na(mirna), gene %in% GST) %>%
  group_by(mirna, gene) %>%
  summarise(evidence = ifelse(any(evidence == "validated"), "validated", "predicted"),
            n_db = n_distinct(database), .groups = "drop")

cat("\n===== TARGET PREDICTIONS =====\n")
cat("Unique miRNA-GST pairs:", nrow(pairs), "\n")
cat("Distinct miRNAs:", n_distinct(pairs$mirna), "\n")
print(table(pairs$evidence))
cat("\nPairs per gene:\n"); print(sort(table(pairs$gene), decreasing = TRUE))

# Require support from at least two databases for predicted pairs. Single-
# database predictions are mostly seed-match noise.
pairs <- pairs %>% dplyr::filter(evidence == "validated" | n_db >= 2)
cat("\nAfter requiring validated OR >=2 databases:", nrow(pairs), "pairs\n")
write.csv(pairs, file.path(OUT, "mirna_gst_target_pairs.csv"), row.names = FALSE)


# =============================================================================
# 2. FETCH miRNA EXPRESSION
# =============================================================================
# Xena TCGA hub. Names there are of the form "hsa-mir-21"; multiMiR returns
# "hsa-miR-21-5p". Normalise both to a common form before matching.

norm_mir <- function(x) {
  x <- tolower(x)
  x <- gsub("-[35]p$", "", x)
  x <- gsub("^hsa-", "", x)
  x
}

HOST_T <- "https://tcga.xenahubs.net"

av <- XenaData %>%
  dplyr::filter(XenaHostNames == "tcgaHub", grepl("BRCA", XenaCohorts),
         grepl("miRNA|mirna", XenaDatasets))
cat("\n===== miRNA DATASETS AVAILABLE =====\n")
print(as.data.frame(av[, c("XenaDatasets","Unit")]))

MIR_DS <- av$XenaDatasets[grepl("miRNA_HiSeq_gene|mirna", av$XenaDatasets)][1]
if (is.na(MIR_DS)) stop("No BRCA miRNA dataset found on the TCGA hub.")
cat("Using:", MIR_DS, "\n")

mf <- file.path(CACHE, "tcga_brca_mirna.rds")
if (file.exists(mf)) {
  MIR <- readRDS(mf)
} else {
  ids <- fetch_dataset_identifiers(HOST_T, MIR_DS)
  cat("miRNAs in dataset:", length(ids), "\n")

  want <- unique(pairs$mirna)
  lookup <- setNames(ids, norm_mir(ids))
  matched <- unique(na.omit(lookup[norm_mir(want)]))
  cat("Matched to dataset:", length(matched), "of", length(want), "\n")
  if (!length(matched)) stop("No miRNA names matched. Inspect head(ids).")

  got <- list()
  for (m in matched) {
    r <- tryCatch(fetch_dense_values(HOST_T, MIR_DS, identifiers = m,
                                     check = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || !nrow(as.matrix(r))) next
    got[[m]] <- as.matrix(r)[1, ]
    if (length(got) %% 25 == 0) message("  ", length(got), " / ", length(matched))
    Sys.sleep(0.15)
  }
  common <- Reduce(intersect, lapply(got, names))
  MIR <- do.call(rbind, lapply(got, function(v) v[common]))
  rownames(MIR) <- names(got)
  saveRDS(MIR, mf)
}
cat("miRNA matrix:", nrow(MIR), "x", ncol(MIR), "\n")


# =============================================================================
# 3. MATCH SAMPLES
# =============================================================================
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
  dplyr::filter(study == "TCGA", grepl("breast", tissue, ignore.case = TRUE),
         grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

# Xena TCGA hub uses 15-character barcodes; align on the first 15 characters.
colnames(MIR) <- substr(colnames(MIR), 1, 15)
clin$short <- substr(clin$sample, 1, 15)

s <- intersect(clin$short, colnames(MIR))
clin <- clin[match(s, clin$short), ]
M <- MIR[, s, drop = FALSE]
G <- G_all[intersect(GST, rownames(G_all)), clin$sample, drop = FALSE]
colnames(G) <- s

cat("\nTumours with both miRNA and GST expression:", length(s), "\n")
print(table(clin$group))

clin$purity <- NA_real_
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- get("Tumor.purity", envir = environment())
  tp$patient <- substr(tp$Sample.ID, 1, 12)
  num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
  clin$purity <- num(tp$CPE)[match(clin$patient, tp$patient)]
}, silent = TRUE)
cat("With purity estimate:", sum(!is.na(clin$purity)), "\n")


# =============================================================================
# 4. TEST THE PREDICTED PAIRS
# =============================================================================
lookup <- setNames(rownames(M), norm_mir(rownames(M)))
pairs$mir_id <- lookup[norm_mir(pairs$mirna)]
test <- pairs %>% dplyr::filter(!is.na(mir_id), gene %in% rownames(G))
cat("\nTestable pairs:", nrow(test), "\n")

partial_rho <- function(x, y, z) {
  ok <- is.finite(x) & is.finite(y) & is.finite(z)
  if (sum(ok) < 50) return(c(NA, NA, NA))
  raw <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
  rx <- residuals(lm(rank(x[ok]) ~ rank(z[ok])))
  ry <- residuals(lm(rank(y[ok]) ~ rank(z[ok])))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman"))
  c(raw, unname(ct$estimate), ct$p.value)
}

res <- lapply(seq_len(nrow(test)), function(i) {
  x <- M[test$mir_id[i], ]; y <- G[test$gene[i], ]
  v <- partial_rho(x, y, clin$purity)
  if (all(is.na(v))) return(NULL)
  data.frame(mirna = test$mirna[i], gene = test$gene[i],
             evidence = test$evidence[i], n_db = test$n_db[i],
             rho_raw = round(v[1], 3), rho_adj = round(v[2], 3),
             p = v[3])
}) %>% bind_rows()

if (!nrow(res)) stop("No pairs testable. Check sample matching above.")

res <- res %>% mutate(FDR = p.adjust(p, "BH"),
                      call = case_when(
                        rho_adj <= -0.3 & FDR < 0.05 ~ "candidate repressor",
                        rho_adj <= -0.15 & FDR < 0.05 ~ "weak negative",
                        rho_adj >= 0.3 & FDR < 0.05 ~ "positive (not repression)",
                        TRUE ~ "no association")) %>%
  arrange(rho_adj)

cat("\n===== miRNA-GST CORRELATION, PURITY-ADJUSTED =====\n")
print(head(as.data.frame(res), 30))
cat("\n"); print(table(res$call))
write.csv(res, file.path(OUT, "mirna_gst_correlation.csv"), row.names = FALSE)

cand <- res %>% dplyr::filter(call == "candidate repressor")

# ---------------------------------------------------------------------------
#  ACCESSION -> NAME
# ---------------------------------------------------------------------------
# multiMiR returns MIMAT accessions, which are unreadable. Converting them is
# not cosmetic: it is what reveals whether the hits are scattered individual
# miRNAs or members of one polycistronic cluster. In this dataset the GSTM3
# repressors resolve largely to miR-17~92 and its paralogues, which is the
# finding.
if (nrow(cand)) {
  if (!requireNamespace("miRBaseConverter", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("miRBaseConverter", ask = FALSE, update = FALSE)
  }
  nm <- tryCatch({
    suppressPackageStartupMessages(library(miRBaseConverter))
    miRNA_AccessionToName(cand$mirna, targetVersion = "v22")$TargetName
  }, error = function(e) { message("Name conversion failed: ",
                                   conditionMessage(e)); NA_character_ })
  cand$mirna_name <- nm
  cat("\n===== CANDIDATE REPRESSORS, NAMED =====\n")
  print(cand[, c("mirna_name","gene","rho_adj","FDR")], row.names = FALSE)
  write.csv(cand, file.path(OUT, "candidate_repressors_named.csv"),
            row.names = FALSE)
}
if (nrow(cand)) {
  cat("\n===== CANDIDATE REPRESSORS =====\n")
  cat("Predicted or validated targets, anticorrelated after purity adjustment.\n\n")
  print(as.data.frame(cand))
  write.csv(cand, file.path(OUT, "candidate_repressors.csv"), row.names = FALSE)
} else {
  cat("\nNo pair reached the candidate threshold. Reportable: predicted miRNA\n")
  cat("targeting of GST genes is not reflected in bulk tumour expression.\n")
}


# =============================================================================
# 5. WITHIN SUBTYPE — does any miRNA track the programme split?
# =============================================================================
# If a miRNA drives the basal/luminal GST split, its correlation should hold
# within subtype and its own expression should differ between subtypes.

if (nrow(cand)) {
  ws <- lapply(SUB, function(sb) {
    idx <- which(clin$group == sb)
    if (length(idx) < 60) return(NULL)
    lapply(seq_len(nrow(cand)), function(i) {
      mid <- lookup[norm_mir(cand$mirna[i])]
      if (is.na(mid)) return(NULL)
      v <- partial_rho(M[mid, idx], G[cand$gene[i], idx], clin$purity[idx])
      if (all(is.na(v))) return(NULL)
      data.frame(subtype = sb, mirna = cand$mirna[i], gene = cand$gene[i],
                 n = length(idx), rho_adj = round(v[2], 3), p = v[3])
    }) %>% bind_rows()
  }) %>% bind_rows()

  if (nrow(ws)) {
    ws <- ws %>% group_by(subtype) %>% mutate(FDR = p.adjust(p, "BH")) %>%
      ungroup()
    cat("\n===== CANDIDATES WITHIN SUBTYPE =====\n")
    print(as.data.frame(ws %>% arrange(mirna, gene, subtype)), digits = 3)
    write.csv(ws, file.path(OUT, "candidates_within_subtype.csv"), row.names = FALSE)

    robust <- ws %>% dplyr::filter(FDR < 0.05, rho_adj <= -0.2) %>%
      group_by(mirna, gene) %>%
      summarise(n_subtypes = n(), mean_rho = round(mean(rho_adj), 3),
                .groups = "drop") %>% arrange(desc(n_subtypes))
    cat("\n===== HOLDING IN MULTIPLE SUBTYPES =====\n")
    if (nrow(robust)) print(as.data.frame(robust)) else
      cat("None - the pooled associations are subtype-driven.\n")
  }

  # Do the candidate miRNAs themselves differ by subtype?
  cm <- unique(na.omit(lookup[norm_mir(cand$mirna)]))
  if (length(cm)) {
    md <- lapply(cm, function(m) data.frame(
      mirna = m, expr = M[m, ], group = factor(clin$group, levels = SUB))
    ) %>% bind_rows() %>% dplyr::filter(is.finite(expr))

    ks <- md %>% group_by(mirna) %>%
      summarise(kw_p = kruskal.test(expr ~ droplevels(group))$p.value,
                Basal = round(median(expr[group == "Basal"]), 2),
                Her2  = round(median(expr[group == "Her2"]), 2),
                LumA  = round(median(expr[group == "LumA"]), 2),
                LumB  = round(median(expr[group == "LumB"]), 2),
                .groups = "drop") %>%
      mutate(FDR = p.adjust(kw_p, "BH")) %>% arrange(kw_p)

    cat("\n===== CANDIDATE miRNA EXPRESSION BY SUBTYPE =====\n")
    print(as.data.frame(ks), digits = 3)
    cat("\nA miRNA that both represses a GST and differs by subtype is a\n")
    cat("candidate driver of the programme split, alongside methylation.\n")
    write.csv(ks, file.path(OUT, "candidate_mirna_by_subtype.csv"), row.names = FALSE)

    p1 <- ggplot(md, aes(group, expr, fill = group)) +
      geom_violin(scale = "width", alpha = .3, colour = NA) +
      geom_boxplot(width = .25, outlier.size = .2, linewidth = .25) +
      facet_wrap(~ mirna, scales = "free_y") +
      scale_fill_manual(values = c(Basal = "#C0392B", Her2 = "#8E7CC3",
                                   LumA = "#2E86AB", LumB = "#E8A33D"), guide = "none") +
      labs(x = NULL, y = "miRNA expression (log2)",
           title = "Candidate GST-repressing miRNAs across PAM50 subtypes") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
    ggsave(file.path(OUT, "candidate_mirna_by_subtype.png"), p1,
           width = 9, height = 6, dpi = 300)
  }
}


# =============================================================================
# 6. FIGURE — the full tested set
# =============================================================================
if (nrow(res) > 5) {
  # Base R throughout. dplyr verbs are masked once multiMiR loads AnnotationDbi,
  # and a figure is not worth a namespace fight.
  hm <- res[, c("mirna", "gene", "rho_adj")]
  hm <- aggregate(rho_adj ~ mirna + gene, data = hm, FUN = mean)
  hm <- reshape(hm, idvar = "mirna", timevar = "gene", direction = "wide")
  rownames(hm) <- hm$mirna; hm$mirna <- NULL
  colnames(hm) <- sub("^rho_adj\\.", "", colnames(hm))
  hm <- as.matrix(hm)

  keep <- apply(hm, 1, function(v) any(abs(v) > 0.15, na.rm = TRUE))
  hm <- hm[keep, , drop = FALSE]

  if (nrow(hm) > 2) {
    hm[!is.finite(hm)] <- 0
    lim <- max(abs(hm))
    png(file.path(OUT, "mirna_gst_heatmap.png"), width = 1800,
        height = max(1200, 40 * nrow(hm)), res = 200)
    pheatmap::pheatmap(hm,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      breaks = seq(-lim, lim, length.out = 101), border_color = "grey85",
      fontsize = 8,
      main = "Predicted miRNA-GST pairs\nSpearman rho, purity-adjusted\n(blue = anticorrelated, consistent with repression)")
    dev.off()
  }
}


cat("\n============================================================\n")
cat("INTERPRETING THIS MODULE\n")
cat("============================================================\n")
cat("Only pairs with a prior reason to exist were tested - predicted or\n")
cat("validated targets - which is what makes the multiple testing tractable.\n\n")
cat("A candidate repressor is a predicted target that is anticorrelated in\n")
cat("about a thousand tumours after purity adjustment. That is a credible\n")
cat("candidate, not proof: demonstrating targeting needs a luciferase assay\n")
cat("or miRNA overexpression, and the text should say so.\n\n")
cat("The most valuable outcome would be a miRNA that both represses a\n")
cat("programme gene AND differs by subtype - that would place miRNA\n")
cat("regulation alongside methylation in the mechanism rather than beside it.\n\n")
cat("A null result is reportable: it would mean predicted miRNA targeting of\n")
cat("GST genes is not reflected in bulk tumour expression, leaving promoter\n")
cat("methylation as the dominant regulatory layer you have identified.\n")
cat("============================================================\n")

writeLines(c(
  paste("Run:", Sys.time()),
  "miRNA-target pairs from multiMiR (validated: miRTarBase, TarBase;",
  "predicted: TargetScan, miRDB, DIANA - top 10% only).",
  "Predicted pairs required support from at least two databases.",
  "miRNA expression from UCSC Xena TCGA hub; GST expression from the Toil",
  "recompute; matched on 15-character TCGA barcodes.",
  "Correlations adjusted for tumour purity (Aran et al. 2015 CPE) via rank",
  "residuals, and repeated within PAM50 subtype.",
  "Correlation indicates candidate regulation only; targeting requires",
  "experimental validation.",
  "", capture.output(sessionInfo())
), file.path(OUT, "provenance_moduleH.txt"))

message("\nWritten to ", normalizePath(OUT))
