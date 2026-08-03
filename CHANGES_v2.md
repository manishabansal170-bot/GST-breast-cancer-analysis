# Changes in v2.0

## Scripts modified

### 01_expression_TCGA_GTEx.R
- Study filters made case-insensitive with `toupper(study)`. The label in the
  Toil phenotype file is "GTEX" in upper case; a case-sensitive match on "GTEx"
  returns nothing and silently drops all 179 GTEx breast samples, producing a
  plausible figure with a missing group.
- Empty sample groups now stop execution rather than warn. An empty group means
  a filter has failed, not that data are missing, and the resulting figure looks
  like real data.

### 02_expression_within_tumour.R
- Study filter made case-insensitive.

### 03_pergene_figures.R
- Study filters made case-insensitive.
- Empty-group guard added.

### 10_survival.R
- Study filter made case-insensitive.

## Scripts added

### 16_regenerate_all_figures.R
Regenerates every published figure from the cache with no downloads. Contains
four corrections:

1. GTEX case sensitivity, as above
2. Detection assessed per sample group rather than on a pooled median, which had
   misclassified GSTA1 as not detected. GSTA1 reaches 7.12 TPM in adjacent
   normal tissue and falls to 0.13 TPM in luminal B, so a pooled median is
   dominated by the tumour fraction
3. Zero-variance genes removed from the SCAN-B panel rather than plotted as
   empty facets with axis ranges of the order 1e-16
4. TCGA survival axis truncated at 10 years; median follow-up is 2.4 years and
   fewer than 25 patients per arm remain at risk beyond that point

### 17_objective2_indian_cohort.R
Objective 2 validation in 104 Indian invasive ductal carcinomas. Halts with an
explanatory message if the patient data file is absent, since that cohort
comprises identifiable clinical material and is not deposited.

## Not changed

Scripts 04 to 09 and 11 to 15 required no correction. Scripts 04, 05, 06, 07,
08, 09 do not filter on the study column; 11 to 15 were verified during the
v1.1 release.

## Effect on published results

None. The figures in the manuscript were generated before the ad-hoc scripts
that introduced the case-sensitivity error, and the deposited scripts 01 and 03
already used the correct upper-case label. These changes harden the code against
a failure mode that would otherwise be silent, and add the two new scripts.
