# Table 3 is 738 by the accounting script but 737 in the manuscript.
# The difference is whether expression is required. Check which.
CACHE <- path.expand("~/GST_BRCA/cache")
suppressPackageStartupMessages(library(TCGAbiolinks))
SUB <- c("Basal","Her2","LumA","LumB"); id <- function(v) substr(v,1,15)
obj <- readRDS(file.path(CACHE,"tcga_methylation_gst.rds"))
met <- obj$met; colnames(met) <- id(colnames(met))
exo <- obj$exp; colnames(exo) <- id(colnames(exo))
st <- TCGAquery_subtype(tumor="brca")
m1 <- colnames(met)[substr(colnames(met),13,15)=="-01"]
g  <- st$BRCA_Subtype_PAM50[match(substr(m1,1,12), st$patient)]
cat("methylation only, four subtypes      :", sum(g %in% SUB), "\n")
m2 <- intersect(m1, colnames(exo))
g2 <- st$BRCA_Subtype_PAM50[match(substr(m2,1,12), st$patient)]
cat("methylation AND expression, four subs:", sum(g2 %in% SUB), "\n")
cat("\nTable 3 reports methylation percentages, which need only methylation,\n")
cat("but its expression column needs both. State which n applies to which column.\n")
