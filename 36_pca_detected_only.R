# =============================================================================
#  Rerun the per-sample PCA on detected genes only
#
#  The manuscript states that PC1's positive loadings are dominated by GSTA2
#  and GSTA5, and that both fall below the detection threshold so their
#  contribution is not distinguishable from noise. A reviewer points out the
#  obvious consequence: if that is true they should not be in the decomposition
#  at all. Scaling a gene whose signal is noise gives it unit variance and lets
#  it load freely on any component.
#
#  Every figure the abstract quotes from this analysis is conditional on that
#  choice. This rerun restricts the matrix to the 16 genes reaching detection.
# =============================================================================
WORKDIR <- "~/GST_BRCA"; CACHE <- file.path(WORKDIR,"cache"); OUT <- file.path(WORKDIR,"figures")
setwd(WORKDIR)
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(TCGAbiolinks)})
set.seed(20260824)
SUB <- c("Basal","Her2","LumA","LumB")

gm <- function(x){ if(is.matrix(x)) return(x)
  if(is.list(x)){ if(!is.null(x$mat)) return(x$mat)
    i<-which(vapply(x,is.matrix,logical(1))); if(length(i)) return(x[[i[1]]])}
  stop("no matrix")}
raw <- readRDS(file.path(CACHE,"toil_gst.rds"))
LIN <- 2^gm(raw)-0.001; LIN[LIN<0] <- 0
E <- log2(LIN+1)

ph <- readRDS(file.path(CACHE,"toil_pheno.rds"))
cl <- function(p) grep(p, colnames(ph), value=TRUE, ignore.case=TRUE)[1]
ph <- ph %>% transmute(sample=.data[[cl("^sample$")]], study=.data[[cl("study")]],
                       tissue=.data[[cl("primary disease or tissue")]], stype=.data[[cl("sample_type")]])
st <- TCGAquery_subtype(tumor="brca")
adjn <- ph %>% dplyr::filter(toupper(study)=="TCGA", grepl("breast",tissue,ignore.case=TRUE),
                             grepl("Solid Tissue Normal",stype)) %>% mutate(group="AdjNormal")
tum <- ph %>% dplyr::filter(toupper(study)=="TCGA", grepl("breast",tissue,ignore.case=TRUE),
                            grepl("Primary Tumor",stype)) %>%
  mutate(patient=substr(sample,1,12), group=st$BRCA_Subtype_PAM50[match(patient,st$patient)]) %>%
  dplyr::filter(group %in% SUB) %>% distinct(patient,.keep_all=TRUE)
meta <- bind_rows(adjn%>%dplyr::select(sample,group), tum%>%dplyr::select(sample,group))
meta <- meta[meta$sample %in% colnames(E),]

# detection: median >= 1 TPM in at least one group, matching Table 1
grp_med <- sapply(unique(meta$group), function(g)
  apply(LIN[, meta$sample[meta$group==g], drop=FALSE], 1, median, na.rm=TRUE))
detected <- rownames(LIN)[apply(grp_med,1,max,na.rm=TRUE) >= 1]
cat("Detected genes:", length(detected), "\n"); cat(paste(sort(detected), collapse=", "), "\n\n")

M <- t(E[detected, meta$sample, drop=FALSE])
M <- M[, apply(M,2,function(v) all(is.finite(v)) && sd(v)>1e-8), drop=FALSE]
M <- scale(M)
cat("Matrix:", nrow(M), "samples x", ncol(M), "genes\n\n")

pc <- prcomp(M, center=FALSE, scale.=FALSE)
ve <- 100*pc$sdev^2/sum(pc$sdev^2)
for (i in 1:6) cat(sprintf("  PC%-2d %5.1f%%   cumulative %5.1f%%\n", i, ve[i], sum(ve[1:i])))

null <- replicate(200, { Mp <- apply(M,2,sample)
  s <- prcomp(Mp,center=FALSE,scale.=FALSE)$sdev^2; 100*s/sum(s) })
thr <- apply(null,1,quantile,0.95)
keep <- which(cumprod(ve>thr)==1)
cat(sprintf("\nComponents retained by parallel analysis: %d\n", length(keep)))
for (i in 1:6) cat(sprintf("  PC%-2d observed %5.1f%%  null 95%% %5.1f%%  %s\n",
                           i, ve[i], thr[i], if(ve[i]>thr[i]) "retained" else "not retained"))

sc <- as.data.frame(pc$x[,1:3]); sc$group <- factor(meta$group, levels=c("AdjNormal",SUB))
for (i in 1:3) {
  cat(sprintf("\nPC%d (%.1f%%) group means:\n", i, ve[i]))
  print(round(tapply(sc[[paste0("PC",i)]], sc$group, mean, na.rm=TRUE),2))
  ld <- sort(pc$rotation[,i])
  cat("  most negative:", paste(names(head(ld,3)),collapse=", "),
      "| most positive:", paste(names(tail(ld,3)),collapse=", "), "\n")
}

write.csv(data.frame(component=paste0("PC",seq_along(ve)), variance_pct=round(ve,2),
                     null_95_pct=round(c(thr,rep(NA,max(0,length(ve)-length(thr)))),2)),
          file.path(OUT,"pca_detected_only_variance.csv"), row.names=FALSE)
write.csv(round(pc$rotation[,1:3],4), file.path(OUT,"pca_detected_only_loadings.csv"))

cat(sprintf("\nFOR THE MANUSCRIPT: PC1 %.1f%%, PC2 %.1f%%, together %.1f%%, %d components retained\n",
            ve[1], ve[2], ve[1]+ve[2], length(keep)))
