GST breast cancer analysis — scripts 32 to 47
=============================================

These are the scripts written during the R3 to R6 revisions, in the state that
actually ran. Where a script failed on first execution and was corrected at the
console, the correction is in the file: the sample-set mismatch in 42, the
kruskal.test guard in 29 and 33, the empty-result guard in 34, and the
supplementary-workbook reader in 43.

RUN ORDER
  Each script is self-contained and reads from ~/GST_BRCA/cache. They do not
  depend on one another except where noted below, so any can be run alone.

  33 must run before 41  — 33 populates cache/all_gst_probe_betas.rds, which 41
                            uses to restrict the QC to probes returning data.
  42 before 46           — 42 establishes the harmonised-matrix result that 46's
                            figure must agree with.

WHAT EACH ONE ANSWERS
  32  Does the methylation contribution survive within subtype, or is it the
      between-subtype contrast? (It survives: 7 of 8 programme genes.)
  33  Do other genes carry the opposing-probe-block structure found in GSTP1?
      (GSTA4, GSTM4 and GSTZ1 do; two of them average across the blocks.)
  34  Do the two expression matrices agree, and is the methylation-versus-copy-
      number comparison fair? (Median rho 0.917; methylation leads on both
      marginal and unique variance.)
  35  Does the GSTM2 response association survive clinical covariates?
      (Yes: adjusted OR 0.569, p = 0.026.)
  36  PCA restricted to detected genes only. (PC1 21.7%, PC2 12.4%, 4 retained.)
  37  Table 3 with expression computed in the same tumours as methylation.
  38  Purity null matched on purity correlation rather than expression.
      NOTE: must be run on breast tumours only. The first run was pan-cancer
      and gave the opposite conclusion.
  39  DepMap breast cell lines: basal versus luminal GSTP1.
  40  Are the GSTM1 and GSTM2 array probes measuring the same thing?
      (Yes: array rho 0.900 against RNA-seq 0.206. Cross-hybridisation.)
  41  450k probe quality: SNPs, design type, sex chromosomes, cross-reactivity.
      Place chen_crossreactive.csv (included) in cache/ for the full check.
  42  All methylation analyses recomputed on the harmonised GRCh38 matrix.
  43  Public decitabine series GSE74251: does GSTP1 re-express?
      (Yes: log2FC 4.08, FDR 6.8e-68, with positive controls confirming the
      treatment worked. Values in the workbook are normalised, not counts.)
  44  Is the methylation subtype contrast a purity artefact? (No: unchanged at
      purity >= 0.7, and beta does not track purity within basal-like.)
  45  Hartigan dip test for bimodality by subtype. (HER2 and LumB only.)
  46  Figure 4 regenerated on the harmonised matrix.
  47  Sample accounting for every analysis. Produces Supplementary Table S14.

  check_table3_n.R        resolves the 737 versus 738 denominator
  check_GSE74251_scale.R  establishes that the decitabine values are normalised
  compute_proteomic_CIs.R adds confidence intervals to the protein concordance

A NOTE ON 39
  Two versions exist. 39_cell_line_check.R attempts an automatic DepMap
  download, which failed against the current release. 39_cell_line_crosscheck.R
  is the one that ran, reading files downloaded by hand. Only the latter is
  included here.
