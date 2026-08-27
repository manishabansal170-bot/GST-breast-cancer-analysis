# =============================================================================
#  41_probe_quality_control.R
#
#  THE OBJECTION
#    The methylation analysis is the mechanistic core of this manuscript, and
#    the whole GSTP1 island result rests on two probes: cg04920951 and
#    cg02659086. The Methods describe how probes were selected by genomic
#    annotation but say nothing about probe quality.
#
#    Two failure modes matter on the 450k array:
#
#    CROSS-REACTIVE PROBES map to more than one genomic location, so their
#    beta values are a blend of several loci. Chen and colleagues (2013)
#    identified roughly 30,000 such probes.
#
#    SNP-OVERLAPPING PROBES sit on a common polymorphism, so beta reflects
#    genotype rather than methylation. This matters particularly for a family
#    with common germline deletions.
#
#    Neither is checked anywhere in the current analysis.
#
#  WHAT THIS SCRIPT DOES
#    Uses the SNP annotation distributed with the Illumina annotation package,
#    which gives the distance from each probe and its extension base to the
#    nearest common SNP, and reports every probe used in this manuscript
#    against it. It also flags probes on the X and Y chromosomes, and reports
#    the probe design type, since Type I and Type II probes have different
#    beta distributions and the 0.1 and 0.3 thresholds are sensitive to that.
#
#  A LIMIT WORTH STATING
#    The definitive cross-reactive list is a published supplementary file, not
#    part of any package. This script checks what the annotation provides and
#    names the two probes that carry the manuscript, so that a reader or
#    reviewer can look them up against Chen et al. directly. Where the check
#    cannot be made from the data to hand, the script says so rather than
#    implying a clean result.
# =============================================================================

WORKDIR <- "~/GST_BRCA"
CACHE   <- file.path(WORKDIR, "cache")
OUT     <- file.path(WORKDIR, "objective1", "moduleC_methylation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(WORKDIR)

pk <- c("dplyr")
for (p in pk) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (!requireNamespace(bioc, quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(bioc, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  lapply(pk, library, character.only = TRUE)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) })

CARRY <- c("cg04920951", "cg02659086")   # the two the GSTP1 result rests on


# =============================================================================
# 1. THE PROBES ACTUALLY USED
# =============================================================================
obj <- readRDS(file.path(CACHE, "tcga_methylation_gst.rds"))
GENES <- rownames(obj$met)

ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
PROMOTER <- c("TSS1500","TSS200","5'UTR","1stExon")

probes_for <- function(g) {
  hit <- ann[grepl(paste0("(^|;)", g, "(;|$)"), ann$UCSC_RefGene_Name), , drop = FALSE]
  if (!nrow(hit)) return(character(0))
  keep <- sapply(strsplit(hit$UCSC_RefGene_Group, ";"), function(v) any(v %in% PROMOTER))
  rownames(hit)[keep]
}

used <- unique(unlist(lapply(GENES, probes_for)))
cat("Promoter probes annotated across the", length(GENES), "genes:", length(used), "\n")

# restrict to those that actually returned data, since the rest are irrelevant
pf <- file.path(CACHE, "all_gst_probe_betas.rds")
if (file.exists(pf)) {
  B <- readRDS(pf)
  with_data <- intersect(used, rownames(B))
  with_data <- with_data[apply(B[with_data, , drop = FALSE], 1,
                               function(v) sum(is.finite(v)) >= 100)]
  cat("Of those, returning data in this cohort:", length(with_data), "\n\n")
} else {
  with_data <- used
  cat("Probe beta cache not found; checking all annotated probes.\n\n")
}


# =============================================================================
# 2. SNP OVERLAP
# =============================================================================
# Probe_rs and CpG_rs give the nearest common SNP to the probe body and to the
# interrogated CpG; CpG_maf and Probe_maf give its minor allele frequency. A
# SNP at the CpG itself is the serious case, because beta then reports
# genotype.
cat(strrep("=", 76), "\nSNP OVERLAP\n", strrep("=", 76), "\n", sep = "")

sn <- ann[with_data, grep("_rs$|_maf$|SBE", colnames(ann), value = TRUE), drop = FALSE]
sn$probe <- rownames(sn)

has_cpg_snp <- !is.na(sn$CpG_rs)
has_sbe_snp <- if ("SBE_rs" %in% names(sn)) !is.na(sn$SBE_rs) else rep(FALSE, nrow(sn))
has_probe_snp <- !is.na(sn$Probe_rs)

cat(sprintf("Probes with a SNP at the interrogated CpG      : %d\n", sum(has_cpg_snp)))
cat(sprintf("Probes with a SNP at the single-base extension : %d\n", sum(has_sbe_snp)))
cat(sprintf("Probes with a SNP anywhere in the probe body   : %d\n", sum(has_probe_snp)))

flagged <- sn$probe[has_cpg_snp | has_sbe_snp]
if (length(flagged)) {
  cat("\nProbes with a SNP at the CpG or extension base, which are the positions\n")
  cat("that make beta unreliable:\n")
  fl <- ann[flagged, c("chr","pos","UCSC_RefGene_Name","Relation_to_Island"), drop = FALSE]
  fl$CpG_rs <- ann[flagged, "CpG_rs"]; fl$CpG_maf <- round(ann[flagged, "CpG_maf"], 3)
  print(fl)
} else cat("\nNo probe used here carries a SNP at the CpG or extension base.\n")


# =============================================================================
# 3. THE TWO PROBES THE GSTP1 RESULT RESTS ON
# =============================================================================
cat("\n", strrep("=", 76), "\nTHE TWO PROBES CARRYING THE GSTP1 RESULT\n",
    strrep("=", 76), "\n", sep = "")

for (p in CARRY) {
  if (!(p %in% rownames(ann))) { cat(p, ": not in the annotation\n"); next }
  a <- ann[p, ]
  cat(sprintf("\n%s\n", p))
  cat(sprintf("  position        : %s:%s\n", a$chr, format(a$pos, big.mark = ",")))
  cat(sprintf("  region          : %s, %s\n", a$UCSC_RefGene_Group, a$Relation_to_Island))
  cat(sprintf("  design type     : %s\n", a$Type))
  cat(sprintf("  SNP at CpG      : %s\n", ifelse(is.na(a$CpG_rs), "none", a$CpG_rs)))
  cat(sprintf("  SNP at extension: %s\n",
              if ("SBE_rs" %in% names(a)) ifelse(is.na(a$SBE_rs), "none", a$SBE_rs) else "not annotated"))
  cat(sprintf("  SNP in probe    : %s\n", ifelse(is.na(a$Probe_rs), "none", a$Probe_rs)))
}


# =============================================================================
# 4. PROBE DESIGN TYPE, WHICH AFFECTS THE BETA THRESHOLDS
# =============================================================================
cat("\n", strrep("=", 76), "\nPROBE DESIGN TYPE\n", strrep("=", 76), "\n", sep = "")
cat("Type I and Type II probes have different beta distributions. Thresholds of\n")
cat("0.1 and 0.3 are therefore not exactly equivalent between them, and a set\n")
cat("mixing the two should say so.\n\n")
print(table(ann[with_data, "Type"]))

cat("\nSex chromosome probes among those used:",
    sum(ann[with_data, "chr"] %in% c("chrX","chrY")), "\n")


# =============================================================================
# 5. CROSS-REACTIVE PROBES
# =============================================================================
# The Chen et al. (2013) list of 29,233 co-hybridising probes is published as a
# supplementary file rather than in a package. Place it in the cache as
# chen_crossreactive.csv, with a TargetID column, and this section runs.
cat("\n", strrep("=", 76), "\nCROSS-REACTIVE PROBES\n", strrep("=", 76), "\n", sep = "")
cf <- file.path(CACHE, "chen_crossreactive.csv")
if (file.exists(cf)) {
  xr <- read.csv(cf, stringsAsFactors = FALSE)
  idc <- grep("TargetID|probe|IlmnID", names(xr), value = TRUE, ignore.case = TRUE)[1]
  xrp <- unique(xr[[idc]])
  cat("Cross-reactive probes in the list:", length(xrp), "\n")
  hit <- intersect(with_data, xrp)
  cat("Of the", length(with_data), "probes used here,", length(hit), "are cross-reactive.\n")
  if (length(hit)) {
    h <- ann[hit, c("chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island")]
    print(h)
  }
  cat("\nThe two probes carrying the GSTP1 result:\n")
  for (p in CARRY)
    cat(sprintf("  %s: %s\n", p, ifelse(p %in% xrp, "CROSS-REACTIVE", "not listed")))
} else {
  cat("chen_crossreactive.csv not in the cache; this check was not run.\n")
  cat("Obtain it from the Weksberg lab supplementary files and rerun.\n")
}


# =============================================================================
# 6. WHAT CANNOT BE CHECKED HERE
# =============================================================================
cat("\n", strrep("=", 76), "\nWHAT THIS CHECK DOES NOT COVER\n", strrep("=", 76), "\n", sep = "")
cat("Whether the cBioPortal betas are normalised for probe design is a property\n")
cat("of the distributed file rather than something derivable from it, and should\n")
cat("be stated from the data portal's own documentation.\n")

out <- data.frame(probe = with_data,
                  chr = ann[with_data, "chr"],
                  pos = ann[with_data, "pos"],
                  gene = ann[with_data, "UCSC_RefGene_Name"],
                  region = ann[with_data, "UCSC_RefGene_Group"],
                  island = ann[with_data, "Relation_to_Island"],
                  type = ann[with_data, "Type"],
                  CpG_rs = ann[with_data, "CpG_rs"],
                  CpG_maf = ann[with_data, "CpG_maf"],
                  Probe_rs = ann[with_data, "Probe_rs"],
                  stringsAsFactors = FALSE)
write.csv(out, file.path(OUT, "probe_quality_control.csv"), row.names = FALSE)

writeLines(c(paste("Run:", Sys.time()),
             sprintf("Probes used: %d", length(with_data)),
             sprintf("With SNP at CpG: %d", sum(has_cpg_snp)),
             sprintf("With SNP at extension base: %d", sum(has_sbe_snp)),
             sprintf("Type I: %d, Type II: %d",
                     sum(ann[with_data,"Type"] == "I"), sum(ann[with_data,"Type"] == "II")),
             "Cross-reactive list of Chen et al. 2013 not applied; see script notes.",
             "", capture.output(sessionInfo())),
           file.path(OUT, "provenance_probe_qc.txt"))

cat("\nWritten to", normalizePath(OUT), "\n")
