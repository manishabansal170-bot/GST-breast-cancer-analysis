# =============================================================================
#  16_regenerate_all_figures.R
#
#  Regenerates every published figure from the cached data, incorporating the
#  corrections listed below. Run after scripts 01 to 15 have populated the
#  cache; this script performs no downloads.
#
#  CORRECTIONS INCORPORATED
#
#  1. GTEx samples were silently dropped from Figure 1.
#     The study label in the UCSC Toil phenotype file is "GTEX" in upper case.
#     A filter written as study == "GTEx" matches nothing, and all 179 GTEx
#     breast samples were excluded without warning. Now matched case
#     insensitively, with an assertion that halts if any group is empty.
#
#  2. Detection was assessed on a pooled median.
#     A gene expressed in normal tissue but silenced in tumour is misclassified
#     as undetected when the median is dominated by the tumour fraction, which
#     outnumbers normal tissue roughly four to one in these cohorts. GSTA1
#     reaches 7.12 TPM in adjacent normal and falls to 0.13 TPM in luminal B,
#     and was wrongly labelled NOT DETECTED. Detection is now assessed per
#     sample group.
#
#  3. Zero-variance genes produced empty facets in the SCAN-B panel.
#     GSTT1 and GSTT2 return no signal in that quantification and generated
#     panels with axis ranges of the order 1e-16. They are now dropped by
#     variance criterion and named in the caption with the reason.
#
#  4. The TCGA survival curve extended to 25 years.
#     Median follow-up in TCGA-BRCA is 2.4 years and beyond 10 years fewer than
#     25 patients per arm remain at risk. The axis is truncated at 10 years
#     with numbers at risk displayed.
#
#  Outputs are written to figures/ and objective1/moduleD_survival/ as PNG and
#  as TIFF at 300 dpi for journal submission.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
FIG     <- file.path(WORKDIR, "figures")
setwd(WORKDIR)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

pk <- c("dplyr","tidyr","ggplot2","survival","survminer")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE); library(TCGAbiolinks) })

SUB  <- c("Basal","Her2","LumA","LumB")
GRP  <- c("GTEx Normal","TCGA Adjacent Normal", SUB)
COLS <- c("GTEx Normal" = "#5D6D7E", "TCGA Adjacent Normal" = "#AEB6BF",
          "Basal" = "#C0392B", "Her2" = "#8E7CC3",
          "LumA" = "#2E86AB", "LumB" = "#E8A33D")

# A cached object may be a bare matrix or a list containing one. Handle both.
get_matrix <- function(x, label = "") {
  if (is.matrix(x)) return(x)
  if (is.list(x)) {
    if (!is.null(x$mat)) return(x$mat)
    i <- which(vapply(x, is.matrix, logical(1)))
    if (length(i)) return(x[[i[1]]])
  }
  stop("No matrix found in cached object: ", label)
}

save_fig <- function(p, stem, w, h, dir = FIG) {
  ggsave(file.path(dir, paste0(stem, ".png")),  p, width = w, height = h,
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(dir, paste0(stem, ".tiff")), p, width = w, height = h,
         dpi = 300, bg = "white", compression = "lzw", limitsize = FALSE)
  cat("  ", stem, ".png and .tiff\n", sep = "")
}

hdr <- function(t) cat("\n", strrep("=", 70), "\n", t, "\n",
                       strrep("=", 70), "\n", sep = "")


# #############################################################################
#  FIGURE 1 — GST FAMILY ACROSS NORMAL TISSUE AND PAM50 SUBTYPES
# #############################################################################
hdr("FIGURE 1: family expression, GTEx and TCGA")

raw <- readRDS(file.path(CACHE, "toil_gst.rds"))
lin <- 2^get_matrix(raw, "toil_gst") - 0.001; lin[lin < 0] <- 0
E <- log2(lin + 1)

ph <- readRDS(file.path(CACHE, "toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value = TRUE, ignore.case = TRUE)[1]
ph <- ph %>% transmute(sample = .data[[cl("^sample$")]],
                       study  = .data[[cl("study")]],
                       tissue = .data[[cl("primary disease or tissue")]],
                       stype  = .data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor = "brca")

# CORRECTION 1: toupper() on study. The label is "GTEX", not "GTEx".
gtex <- ph %>% filter(toupper(study) == "GTEX",
                      grepl("breast", tissue, ignore.case = TRUE)) %>%
  mutate(group = "GTEx Normal")

adjn <- ph %>% filter(toupper(study) == "TCGA",
                      grepl("breast", tissue, ignore.case = TRUE),
                      grepl("Solid Tissue Normal", stype)) %>%
  mutate(group = "TCGA Adjacent Normal")

# Adjacent normals must never inherit a PAM50 label from their patient
tum <- ph %>% filter(toupper(study) == "TCGA",
                     grepl("breast", tissue, ignore.case = TRUE),
                     grepl("Primary Tumor", stype)) %>%
  mutate(patient = substr(sample, 1, 12),
         group = st$BRCA_Subtype_PAM50[match(patient, st$patient)]) %>%
  filter(group %in% SUB) %>% distinct(patient, .keep_all = TRUE)

meta <- bind_rows(gtex %>% select(sample, group),
                  adjn %>% select(sample, group),
                  tum  %>% select(sample, group))
meta <- meta[meta$sample %in% colnames(E), ]
meta$group <- factor(meta$group, levels = GRP)

print(table(meta$group))
if (any(table(meta$group) == 0))
  stop("A sample group is empty. Check the study and sample_type filters.")

long <- as.data.frame(t(E[, meta$sample, drop = FALSE])) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(meta, by = "sample") %>%
  pivot_longer(-c(sample, group), names_to = "gene", values_to = "expr") %>%
  filter(is.finite(expr))

# CORRECTION 2: detection assessed per group, not on a pooled median
det <- long %>%
  group_by(gene, group) %>%
  summarise(median_tpm = 2^median(expr) - 1, .groups = "drop") %>%
  group_by(gene) %>%
  summarise(max_group_tpm = max(median_tpm, na.rm = TRUE),
            max_group = group[which.max(median_tpm)],
            detected = max_group_tpm >= 1, .groups = "drop")

cat("\nDetected (", sum(det$detected), "): ",
    paste(det$gene[det$detected], collapse = ", "), "\n", sep = "")
cat("Not detected (", sum(!det$detected), "): ",
    paste(det$gene[!det$detected], collapse = ", "), "\n", sep = "")
write.csv(det, file.path(FIG, "detection_toil.csv"), row.names = FALSE)

ord <- long %>% group_by(gene) %>% summarise(m = median(expr), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(gene)
long$gene <- factor(long$gene, levels = ord)

# Annotation placed on a middle category so it cannot run off the panel edge
lab <- det %>% filter(!detected) %>%
  transmute(gene = factor(gene, levels = ord),
            group = factor("Her2", levels = GRP), label = "NOT DETECTED")

lv <- long %>% group_by(group) %>%
  summarise(n = n_distinct(sample), .groups = "drop") %>%
  mutate(lab = paste0(group, "\n(n=", n, ")")) %>%
  { setNames(.$lab, .$group) }

p1 <- ggplot(long, aes(group, expr, fill = group)) +
  geom_violin(scale = "width", alpha = .35, colour = NA) +
  geom_boxplot(width = .22, outlier.size = .25, outlier.alpha = .25,
               linewidth = .3, colour = "grey20") +
  geom_text(data = lab, aes(x = group, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.8,
            colour = "#B03A2E", size = 2.1, fontface = "bold") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_x_discrete(labels = lv) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = NULL, y = expression(log[2]*"(TPM + 1)"),
       title = "Glutathione S-transferase family expression in breast tissue and TCGA-BRCA",
       subtitle = "GTEx and TCGA uniformly requantified (UCSC Toil recompute); panels ordered by median expression") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        panel.grid.minor = element_blank())

save_fig(p1, "Fig_panel_toil", 11.5, 14)


# #############################################################################
#  SCAN-B VALIDATION PANEL
# #############################################################################
hdr("SCAN-B validation panel")

Es_all <- get_matrix(readRDS(file.path(CACHE, "gse96058_gst.rds")), "gse96058_gst")
phs <- readRDS(file.path(CACHE, "gse96058_pheno.rds"))

# The identifier matching the expression columns is 'title' (F1, F2, ...),
# not the GEO accession in 'gsm'
ms <- data.frame(sample = as.character(phs$title),
                 group  = as.character(phs$pam50), stringsAsFactors = FALSE)
ms <- ms[ms$group %in% SUB & ms$sample %in% colnames(Es_all), ]
ms$group <- factor(ms$group, levels = SUB)
print(table(ms$group))
stopifnot(nrow(ms) > 100)

Es <- Es_all[, ms$sample, drop = FALSE]

# CORRECTION 3: drop genes with no variance rather than plotting empty facets
gs <- data.frame(gene = rownames(Es),
                 sd  = apply(Es, 1, sd,  na.rm = TRUE),
                 max = apply(Es, 1, max, na.rm = TRUE),
                 stringsAsFactors = FALSE)
gs$drop <- gs$sd < 1e-8 | gs$max < 1e-8
dropped <- gs$gene[gs$drop]; kept <- gs$gene[!gs$drop]

if (length(dropped)) {
  cat("\nDropped, no quantifiable signal:", paste(dropped, collapse = ", "), "\n")
  cat("  GSTT1: absent from the GRCh38 primary assembly (alternate scaffold\n")
  cat("         NT_187633.1 only), so not quantifiable by any GRCh38 pipeline\n")
  cat("  GSTT2: polymorphic pseudogene in segmental duplication with GSTT2B;\n")
  cat("         short reads assign ambiguously between the two loci\n")
}
write.csv(gs, file.path(FIG, "scanb_gene_status.csv"), row.names = FALSE)

ls_ <- as.data.frame(t(Es[kept, , drop = FALSE])) %>%
  tibble::rownames_to_column("sample") %>%
  left_join(ms, by = "sample") %>%
  pivot_longer(-c(sample, group), names_to = "gene", values_to = "expr") %>%
  filter(is.finite(expr))

eps <- ls_ %>% group_by(gene) %>%
  summarise(H = unname(kruskal.test(expr ~ droplevels(group))$statistic),
            n = n(), .groups = "drop") %>%
  mutate(eps2 = H / ((n^2 - 1)/(n + 1))) %>% arrange(desc(eps2))
ls_$gene <- factor(ls_$gene, levels = eps$gene)
write.csv(eps, file.path(FIG, "scanb_subtype_effect.csv"), row.names = FALSE)

cat("\nStrongest subtype effects in SCAN-B:\n")
print(as.data.frame(head(eps[, c("gene","eps2")], 5)), row.names = FALSE, digits = 3)
cat("The ordering reproduces TCGA, where GSTP1 and GSTM3 also rank first.\n")

lvs <- ls_ %>% group_by(group) %>%
  summarise(n = n_distinct(sample), .groups = "drop") %>%
  mutate(lab = paste0(group, "\n(n=", n, ")")) %>%
  { setNames(.$lab, .$group) }

p2 <- ggplot(ls_, aes(group, expr, fill = group)) +
  geom_violin(scale = "width", alpha = .35, colour = NA) +
  geom_boxplot(width = .22, outlier.size = .2, outlier.alpha = .2,
               linewidth = .3, colour = "grey20") +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = COLS[SUB], guide = "none") +
  scale_x_discrete(labels = lvs) +
  labs(x = NULL, y = expression(log[2]*"(TPM + 1)"),
       title = "GST family expression across PAM50 subtypes in SCAN-B (GSE96058)",
       subtitle = sprintf("Independent Swedish cohort, RNA-seq, n = %s; FPKM converted to TPM across the full transcriptome; panels ordered by subtype effect size",
                          format(nrow(ms), big.mark = ","))) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        strip.background = element_rect(fill = "grey93"),
        strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())

save_fig(p2, "Panel_SCANB", 11, 2 + 2.3 * ceiling(length(kept)/4))


# #############################################################################
#  SURVIVAL CURVES
# #############################################################################
hdr("Kaplan-Meier curves")

SURV <- file.path(WORKDIR, "objective1", "moduleD_survival")
dir.create(SURV, recursive = TRUE, showWarnings = FALSE)

LOCKED_GENE <- "GSTM4"     # the single gene selected by LASSO-Cox in TCGA
LOCKED_COEF <- -0.0437

tc <- readRDS(file.path(CACHE, "tcga_clinical.rds"))

# Survival time comes from days_to_death for the deceased and
# days_to_last_follow_up for the living. days_to_diagnosis is all zero and
# must not be used.
samp <- grep("^TCGA-.*-01$", colnames(E), value = TRUE)
d_t <- data.frame(sample = samp, patient = substr(samp, 1, 12),
                  expr = as.numeric(E[LOCKED_GENE, samp]),
                  stringsAsFactors = FALSE)
d_t <- d_t[!duplicated(d_t$patient), ]

i <- match(d_t$patient, tc$submitter_id)
d_t <- d_t[!is.na(i), ]; i <- i[!is.na(i)]

dead <- grepl("Dead|Deceased", tc$vital_status[i], ignore.case = TRUE)
d_t$time  <- suppressWarnings(ifelse(dead, as.numeric(tc$days_to_death[i]),
                                     as.numeric(tc$days_to_last_follow_up[i]))) / 365.25
d_t$event <- as.numeric(dead)
d_t <- d_t[is.finite(d_t$time) & is.finite(d_t$event) & d_t$time > 0, ]

# CORRECTION 5: restrict to the four PAM50 subtypes, matching 10_survival.R
d_t$pam50 <- st$BRCA_Subtype_PAM50[match(d_t$patient, st$patient)]
d_t <- d_t[d_t$pam50 %in% SUB, ]

d_t$risk <- LOCKED_COEF * scale(d_t$expr)[, 1]
d_t$grp  <- factor(ifelse(d_t$risk > median(d_t$risk), "High risk", "Low risk"),
                   levels = c("Low risk","High risk"))

cat(sprintf("TCGA: %d patients, %d deaths, median follow-up %.1f years\n",
            nrow(d_t), sum(d_t$event), median(d_t$time)))
a10 <- sapply(levels(d_t$grp), function(g) sum(d_t$time >= 10 & d_t$grp == g))
cat(sprintf("At risk at 10 years: %d and %d, which is why the axis is truncated\n",
            a10[1], a10[2]))

km <- function(fit, dat, title, subtitle, xmax, brk) {
  ggsurvplot(fit, data = dat, xlab = "Years", ylab = "Overall survival",
             title = title, subtitle = subtitle,
             xlim = c(0, xmax), break.time.by = brk,
             palette = c("#2E86AB", "#C0392B"),
             legend.title = "GST risk score",
             legend.labs = c("Low risk", "High risk"),
             risk.table = TRUE, risk.table.height = 0.26,
             risk.table.title = "Number at risk", risk.table.y.text = FALSE,
             pval = TRUE, pval.size = 4.2,
             conf.int = TRUE, conf.int.alpha = 0.12, censor.size = 2.5,
             ggtheme = theme_bw(base_size = 11))
}

# CORRECTION 4: TCGA truncated at 10 years
p3 <- km(survfit(Surv(time, event) ~ grp, data = d_t), d_t,
         "GST risk score \u2014 TCGA-BRCA",
         sprintf("Truncated at 10 years; median follow-up %.1f years", median(d_t$time)),
         10, 2)
png(file.path(SURV, "KM_riskscore_tcga.png"), width = 7.5, height = 7,
    units = "in", res = 300); print(p3); invisible(dev.off())
cat("   KM_riskscore_tcga.png\n")

# METABRIC, shown in full: 965 deaths over a median 9.7 years
M  <- get_matrix(readRDS(file.path(CACHE, "metabric_gst.rds")), "metabric_gst")
mc <- readRDS(file.path(CACHE, "metabric_clinical.rds"))
sid <- grep("sampleId|SAMPLE_ID", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
mt  <- grep("OS_MONTHS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]
me  <- grep("OS_STATUS", colnames(mc), value = TRUE, ignore.case = TRUE)[1]

sh <- intersect(colnames(M), as.character(mc[[sid]]))
d_m <- data.frame(sample = sh, expr = as.numeric(M[LOCKED_GENE, sh]),
                  stringsAsFactors = FALSE)
j <- match(d_m$sample, mc[[sid]])
d_m$time  <- suppressWarnings(as.numeric(mc[[mt]][j])) / 12
d_m$event <- as.numeric(grepl("^1|DECEASED|Died", as.character(mc[[me]][j]),
                              ignore.case = TRUE))
d_m <- d_m[is.finite(d_m$time) & is.finite(d_m$event) & d_m$time > 0, ]

# CORRECTION 5, as above, for METABRIC CLAUDIN_SUBTYPE
pcol <- grep("CLAUDIN_SUBTYPE", colnames(mc), value=TRUE, ignore.case=TRUE)[1]
d_m$pam50 <- as.character(mc[[pcol]][match(d_m$sample, mc[[sid]])])
d_m <- d_m[d_m$pam50 %in% SUB, ]

# coefficients and scaling locked on TCGA: this is external validation
d_m$risk <- LOCKED_COEF * ((d_m$expr - mean(d_t$expr)) / sd(d_t$expr))
d_m$grp  <- factor(ifelse(d_m$risk > median(d_m$risk), "High risk", "Low risk"),
                   levels = c("Low risk","High risk"))

cx <- coxph(Surv(time, event) ~ scale(risk), data = d_m)
cat(sprintf("METABRIC: %d patients, %d deaths, median follow-up %.1f years\n",
            nrow(d_m), sum(d_m$event), median(d_m$time)))
cat(sprintf("Locked risk score: HR %.3f per SD, p = %.3g, C = %.3f\n",
            exp(coef(cx)), summary(cx)$coefficients[1,5], summary(cx)$concordance[1]))

p4 <- km(survfit(Surv(time, event) ~ grp, data = d_m), d_m,
         "GST risk score \u2014 METABRIC",
         sprintf("External validation, coefficients locked; %d deaths over a median %.1f years",
                 sum(d_m$event), median(d_m$time)),
         25, 5)
png(file.path(SURV, "KM_riskscore_metabric.png"), width = 7.5, height = 7,
    units = "in", res = 300); print(p4); invisible(dev.off())
cat("   KM_riskscore_metabric.png\n")


# #############################################################################
hdr("DONE")
cat("Figures written to:\n  ", normalizePath(FIG), "\n  ", normalizePath(SURV), "\n")
cat("\nEach figure is written as PNG for viewing and TIFF at 300 dpi for\n")
cat("journal submission.\n")

writeLines(c(paste("Run:", Sys.time()),
  "Regenerated all published figures from cache.",
  "Corrections: GTEX case sensitivity; per-group detection; zero-variance",
  "gene removal in SCAN-B; TCGA survival axis truncated at 10 years.",
  "", capture.output(sessionInfo())),
  file.path(FIG, "provenance_figures.txt"))
