# GST-breast-cancer-analysis

Analysis code for subtype-specific epigenetic regulation of the glutathione
S-transferase family in breast cancer.

**Zenodo (concept DOI, resolves to the latest version):**
[10.5281/zenodo.21676592](https://doi.org/10.5281/zenodo.21676592)

---

## What this repository contains

Code for a two-part study.

**Objective 1** profiles 22 GST and GST-like genes across 5,833 breast tumours
in three independent cohorts, integrating expression with promoter methylation,
somatic copy number, mutation, proteomics, microRNA and immune deconvolution,
and assessing prognostic and predictive clinical significance.

**Objective 2** validates the principal finding in an independent cohort of 104
Indian patients assayed by methylation-specific PCR.

All computational analyses run from public data. Scripts download and cache what
they need on first run.

---

## Principal findings

Breast tumours alter GST expression along two independent axes: a universal loss
of mu- and theta-class capacity, and a subtype-specific redistribution of what
remains.

Subtype-specific promoter methylation determines which programme survives. GSTP1
is unmethylated in 97% of basal-like tumours and methylated in 58–63% of luminal
B and HER2-enriched tumours, with corresponding differences in expression. Copy
number, mutation, NRF2 activity and microRNA targeting were each tested and
excluded as explanations.

The family carries no clinically meaningful prognostic information, with a
maximum incremental concordance index of 0.0054 over a clinical baseline. Mu
class expression predicts reduced pathological complete response to neoadjuvant
chemotherapy in oestrogen-receptor-positive disease.

In an independent Indian cohort, GSTP1 promoter methylation was absent from all
16 triple-negative tumours and present in 20.5% of the remainder. HIC1
methylation was, by contrast, more frequent in triple-negative disease,
establishing that the deficit is gene-specific.

---

## Running order

Scripts are numbered in dependency order. Scripts 01 and 02 populate the cache
that later scripts rely on.

| Script | Produces |
|---|---|
| `01_expression_TCGA_GTEx.R` | Family expression across normal tissue and PAM50 subtypes; per-group detection |
| `02_expression_within_tumour.R` | Subtype comparison, effect sizes |
| `03_pergene_figures.R` | Individual gene figures with significance brackets |
| `04_validation_METABRIC.R` | METABRIC replication |
| `05_validation_SCANB.R` | SCAN-B replication, FPKM to TPM conversion |
| `06_methylation_expression.R` | Methylation–expression correlation |
| `07_methylation_purity_subtype.R` | Purity-adjusted partial correlation, bimodality |
| `08_methylation_two_programme_test.R` | The subtype-matching test |
| `09_copynumber_mutation_protein.R` | Variance partition, mutation frequency, proteomics |
| `10_survival.R` | Cox models, incremental C-index, LASSO with external validation |
| `11_pathway_network_coexpression.R` | Functional co-expression, NRF2 target scoring |
| `12_drug_sensitivity.R` | Cell line drug response |
| `13_neoadjuvant_response_GSE25066.R` | Pathological complete response, stratified by receptor status |
| `14_immune_infiltration.R` | Immune and stromal scoring, purity-adjusted |
| `15_mirna_regulation.R` | Predicted and validated miRNA targeting |
| `16_regenerate_all_figures.R` | **Regenerates all published figures from cache** |
| `17_objective2_indian_cohort.R` | Objective 2 validation cohort |

Script 16 is the one to run if you want the published figures. It performs no
downloads.

---

## Notes on the analysis

Several decisions in these scripts are not obvious and were made in response to
specific failure modes. They are recorded here because the reasons are not
apparent from the code alone.

**Expression is taken from the UCSC Toil recompute.** Comparing tumour with
healthy donor tissue requires that both be quantified through an identical
pipeline; otherwise differences in aligner, annotation and library preparation
dominate the comparison. Values are supplied as log2(TPM + 0.001) and are
back-transformed before re-expression as log2(TPM + 1).

**The study label in the Toil phenotype file is "GTEX" in upper case.** A filter
written as `study == "GTEx"` matches nothing and silently drops all 179 GTEx
breast samples. Script 16 matches case-insensitively and asserts that no sample
group is empty.

**Detection is assessed per sample group, not on a pooled median.** A gene
expressed in normal tissue but silenced in tumour is misclassified as undetected
when the median is dominated by the tumour fraction. GSTA1, at 7.12 TPM in
adjacent normal and 0.13 TPM in luminal B, is the case in point.

**Ensembl version suffixes are stripped before identifier mapping.** Omitting
this causes affected genes to return zero across all samples without any
warning.

**Tumour-adjacent normal samples never inherit a PAM50 label.** Subtype is
assigned per patient in TCGA, so a naive join labels the adjacent normal with
its patient's tumour subtype.

**GSTT1 is absent from the GRCh38 primary assembly.** It is annotated only on
the alternate locus scaffold NT_187633.1, the primary assembly having been built
from a haplotype carrying the common germline deletion. No GRCh38-based pipeline
can quantify it. This is a limitation of the reference and not evidence of
absent expression, and the distinction is preserved throughout.

**GSTT2 is a polymorphic pseudogene** forming a segmental duplication with
GSTT2B. Short reads assign ambiguously between the two loci and they are not
interpreted independently.

**cBioPortal returns discrete data sparsely.** For GISTIC copy number calls and
mutations, only altered samples are transmitted. Alteration frequencies are
computed against the number of samples profiled, taken from the study sample
lists, rather than against the number of records returned. Correlation analyses
use the continuous linear copy number profile.

**All composition-sensitive correlations are purity-adjusted.** GSTs are
epithelial, so a stroma-rich specimen shows both lower expression and diluted
bulk methylation. In this dataset 242 of 256 raw GST–immune associations (95%)
did not survive adjustment.

**Enrichment is computed on the family plus its interaction partners.** GO or
KEGG enrichment of 16 GST genes returns "glutathione metabolism", which restates
the input. Co-expression is tested against a panel of processes selected a
priori and independently of the family.

**MicroRNA testing is restricted to predicted or validated target pairs.**
Correlating all microRNAs against all family members entails roughly 30,000
tests and returns noise. Only pairs with prior database support are tested, and
predicted pairs require support from at least two databases.

**Survival inference is from continuous models, not from Kaplan-Meier curves.**
Dichotomisation discards information and optimal cut-point selection inflates
type I error. Curves are presented for illustration, and specifically to show
that a visually convincing separation is compatible with a hazard ratio of 1.16
and a concordance index of 0.560.

**The TCGA survival axis is truncated at 10 years.** Median follow-up is 2.4
years and fewer than 25 patients per arm remain at risk beyond that point.

**Predictive analysis is stratified by receptor status.** Complete response
rates differ approximately threefold between strata, so any marker associated
with receptor status appears associated with response in an unstratified
analysis irrespective of causation.

---

## Data availability

**Public data used by scripts 01 to 16:**

| Source | Access |
|---|---|
| TCGA-BRCA and GTEx expression | UCSC Xena Toil hub |
| METABRIC, TCGA methylation, copy number, mutation, proteomics | cBioPortal API |
| SCAN-B | GEO GSE96058 |
| Neoadjuvant cohort | GEO GSE25066 |
| TCGA miRNA | UCSC Xena TCGA hub |

Scripts download and cache these automatically.

**Objective 2 patient data are not deposited.** The Indian validation cohort
comprises identifiable clinical material collected under institutional ethics
approval, and patient-level data cannot be made public. Script 17 will halt with
an explanatory message if the data file is absent. Aggregate results are in the
manuscript, and enquiries may be directed to the corresponding author.

---

## Requirements

R 4.4 or later. Packages install automatically on first run:

```
dplyr tidyr ggplot2 UCSCXenaTools cBioPortalData TCGAbiolinks GEOquery
survival survminer glmnet timeROC multiMiR miRBaseConverter pROC
hgu133a.db pheatmap RColorBrewer readxl meta metafor
```

On Windows, before installing Bioconductor packages:

```r
Sys.setenv(TAR = "internal")
options(timeout = 1000000)
```

---

## Reproducing the figures

Place the scripts from this repository in a `scripts/` subfolder of your
working directory, so the layout is `~/GST_BRCA/scripts/`, then:

```r
setwd("~/GST_BRCA")
source("scripts/16_regenerate_all_figures.R")
```

Requires the cache populated by scripts 01 to 15. Each figure is written as PNG
for viewing and as TIFF at 300 dpi for submission.

---

## Citation

If you use this code, please cite the Zenodo record:

> Gupta M. GST-breast-cancer-analysis. Zenodo.
> https://doi.org/10.5281/zenodo.21676592

---

## Licence

MIT
