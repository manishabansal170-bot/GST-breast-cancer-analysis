Files missing from the repository as of this check
==================================================

29_gstp1_per_probe.R
    Produces the per-probe GSTP1 result the manuscript relies on: the CpG
    island and north shore form two anti-correlated blocks, and the per-gene
    value corresponds to a single island probe (cg04920951, rho = 1.000).
    Cited in the Results, the Discussion and Supplementary Table S11.
    Contains the kruskal.test guard: eight of the thirteen annotated probes
    return no data, and the unguarded test errors rather than returning NA.

20_validation_GSE32646.R
    The neoadjuvant validation cohort.

25_purity_control_epithelial.R
    The epithelial-matched purity control, reported in the immune section.

22_assemble_submission_folder.R
23_build_main_figure_file.R
    Submission-assembly utilities. Not analytical, but they document how the
    figure files were produced and at what resolution.

check_table3_n.R
    Establishes that Table 3's methylation columns rest on 738 tumours and its
    expression column on 737.

check_GSE74251_scale.R
    Establishes that the decitabine values are normalised, not raw counts.
    The manuscript states this; this is the check behind it.

compute_proteomic_CIs.R
    Adds confidence intervals and FDR to the transcript-protein concordance.
    Produces the 11-of-14 interpretable figure quoted in the Abstract.

chen_crossreactive.csv
    The Chen et al. (2013) list of 29,233 cross-reactive 450k probes, from the
    Weksberg lab supplementary files. Place in cache/ so that script 41 runs
    its cross-reactivity check rather than reporting that it cannot.

ALSO WORTH CHECKING
    51_scale.R is in the repository but is not part of this script series.
    Confirm what it is before the Zenodo deposit freezes it.
