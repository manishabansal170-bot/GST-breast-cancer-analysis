d <- read.csv("~/GST_BRCA/objective1/moduleBG_alterations_protein/mrna_protein_concordance.csv")
z  <- atanh(d$rho_mRNA_protein)
se <- 1/sqrt(d$n - 3)
d$CI_low  <- round(tanh(z - 1.96*se), 3)
d$CI_high <- round(tanh(z + 1.96*se), 3)
d$FDR     <- p.adjust(d$pval, "BH")
d$sig     <- d$FDR < 0.05
d <- d[order(-d$rho_mRNA_protein), ]
print(d[, c("gene","n","rho_mRNA_protein","CI_low","CI_high","FDR","sig")],
      row.names = FALSE, digits = 3)
cat("\nProteins with FDR < 0.05:", sum(d$sig), "of", nrow(d), "\n")
cat("n range:", min(d$n), "to", max(d$n), "\n")
cat("Non-significant:", paste(d$gene[!d$sig], collapse = ", "), "\n")
write.csv(d, "~/GST_BRCA/objective1/moduleBG_alterations_protein/mrna_protein_concordance_CI.csv",
          row.names = FALSE)
