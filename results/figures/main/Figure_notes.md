# Figure panel guide

Panel notes for the Nature-style data figures, written by scripts/39_main_figures.py.

## Fig1

- **a**: Datasets used in the study. Counts come from the existing study description.
- **b**: WNN coordinates from the reanalysis. Colours and labelled counts show eight broad cell types. Black points mark M-6.
- **c**: Cell-type make-up of the nuclei in each donor, with the yield printed above each bar. Squares below mark RA (vermilion) and OA (grey).
- **d**: RNA weights by subtype. Points are subtype summaries, black ticks are cell-type medians and lines run from the lowest to the highest subtype value.
- **e**: Raw benchmark metrics on a 0–1 scale. All columns use the same colour scale (white = 0, saturated blue = 1). The overall score is taken as given by the analysis.
- **f**: Percentage of variance explained by each of the first eight MOFA+ factors.
- **g**: Median genes versus median peaks for each donor. Point area is proportional to the number of retained nuclei. RNA and ATAC QC counts may differ slightly.

## Fig2

- **a**: Colour shows S-LDSC tau z and point area shows nominal −log10 P, capped at 6. A black ring marks FDR < 0.05.
- **b**: Top 10 subtypes by mean RA S-LDSC statistic in the two European RA cohorts.
- **c**: Fine-mapped PIP mass relative to a permutation null. Shape and colour show the GWAS, point area shows nominal P and a black outline marks FDR < 0.05.
- **d**: Mean RA S-LDSC versus mean RA scDRS for states that have both estimates, with the Spearman correlation across states. M-6 has no S-LDSC estimate and is left out.
- **e**: Top variants, chosen by the highest supplied PIP with one variant per nearest protein-coding gene. Labels include the peak state. Grey dashes mark missing study records, not zero PIP. Exact variant IDs are in the source table saved with this figure.
- **f**: Paired B-cell neighbourhood means within each section, with the positive/negative fold ratio. Lines join values from the same section.
- **g**: Xenium 5K section CT_21045-A-15, read from a processed plotting table made by step 20. The outlined square is the ROI shown in h. A cell counts as CXCL13-positive when its raw count is >0.
- **h**: 1.2-mm ROI centred on the densest 500-µm B-cell bin. It is chosen by a reproducible rule to show the structure of the aggregate, not to measure enrichment. Same cells and coordinates as g.
- **i**: Empirical permutation P for each section, capped at 1/1001, the smallest value possible with 1,000 permutations. Global label shuffling does not keep the spatial clustering of labels or the tissue compartments. These associations do not show causal organisation.

## Fig3

- **a**: DESeq2 active/remission effects and Wald statistics from the original analysis. P is recomputed from the stored Wald statistic so that no P value is rounded to zero. FDR status uses the supplied adjusted P.
- **b**: S1 is defined by 15 discovery genes. The 13 measurable genes give the RNA and ATAC gene-activity z scores shown for all 14 myeloid states, on one shared centred scale.
- **c**: S1 discovery-signature scores for each donor from RNA and ATAC gene activity (13 of 15 genes measurable), with the Spearman correlation and the supplied random-gene-set test. The line is a descriptive least-squares fit.
- **d**: The eight Hallmark pathways with the highest disease NES, with treatment NES for the same pathways. * marks supplied GSEA FDR < 0.05. Disease = active/remission macrophages; treatment = post/pre bulk in good responders.
- **e**: Paired regression coefficients for each section, for IFNG-positive and IFNG-negative T/NK neighbours. P values are one-sided spatial permutation tests from the existing table.
- **f**: S2 is the predefined 10-gene IFN programme. Each dot is a section and black bars are medians. The comparison shown is the supplied donor-level RA/healthy test, so section counts are not independent donor counts.
- **g**: Hallmark IFN-gamma protein score. Sublining joins immune_sublining, fibrous_sublining and sublining before the average within each donor, as in the statistical test (RA n=6, NH n=3 in each compartment). Black bars are medians. Two-sided Mann-Whitney P values are recomputed and checked against the source table. Protein-set specificity is supported in lining but not in sublining (random-set P=0.339).
- **h**: S2 is the predefined 10-gene IFN programme and is different from S1. Points are single biopsy scores. Boxes show the IQR and median, and whiskers 1.5×IQR. Counts are biopsies. The supplied random-gene-set comparison is shown.
- **i**: S1 is the 15-gene discovery signature (13 measured in this cohort). Points are the paired change for each patient after technical residualisation. Boxes show the IQR and median, and whiskers 1.5×IQR. The P values shown are supplied paired Wilcoxon tests, not comparisons between response groups.

## Fig4

- **a**: Differences in transcription-factor RNA regulon activity and chromVAR motif activity. Candidates pass BH q <0.05 across all tested RNA gene-label-permutation P values and rank in the top motif decile. Data come from two donors, and gene-label permutations are not donor replication. These are candidate regulators.
- **b**: Top 10 factors by summed RNA and motif rank. Point area = 15 + 2 sqrt(regulon target count). Each column has its own labelled colour scale in its own units (RNA 0–7, ATAC 0–6). A black outline marks RNA gene-label-permutation BH q <0.05. Printed numbers are unscaled effects.
- **c**: Highest IRF5 H4 for each project/biosample/condition and RA lead signal. The top 10 conditions by highest H4 are shown. Grey dashes mark missing records.
- **d**: Existing MR odds ratios and supplied 95% CI for six genes. All supplied estimates are shown, including more than one lead signal for some gene/GWAS pairs; numbers in brackets tell these apart, in source-table order. The right track gives same-direction/total datasets. Colour shows direction and grey marks agreement <75%. These are exploratory estimates, not measures of drug efficacy.
- **e**: Fraction of immune QTL datasets in which the risk allele raises expression. Point area and the number on the right give the number of datasets. A black ring marks the supplied consistency flag.
- **f**: Top 24 genes across the 11 original evidence layers, on a shared 0–1 scale. The original evidence score is printed above. White = 0, mid blue = 0.5, saturated blue = 1. A green marker below a gene shows low prior RA evidence, based on the supplied novelty label. The evidence layers are not assumed to be statistically independent.
- **g**: All prioritised genes: original integrated score versus Open Targets association score. Point area = 5 + 2 × number of evidence layers. Green marks genes with the supplied low-prior-evidence label.

## ExtData_Fig1

- **a**: Mean AUC ± fold SD from the revised repeated nested CV, with all fitting done on training data only. Error bars are SD, not CI. For random controls, the error is the SD across the means of the draws. Observed and label-permuted F1 use the same 5 × 5 outer CV, 3-fold tuning and preprocessing. These are internal estimates from one cohort.
- **b**: AUC for good response versus nonresponse at baseline, with the supplied bootstrap 95% CI. The predictive signal seen at baseline is confounded by technical factors.
- **c**: Fraction of perturbations that meet each original rule. The right column gives the number of perturbations. Metric definitions and ranges are in robustness_summary.tsv.
- **d**: Observed share of anti-inflammatory compound classes versus the random median. P values are supplied empirical null comparisons.
- **e**: Existing RNA and ATAC PAGA connectivity to M-6. Closeness in the graph does not show developmental or temporal direction.
- **f**: Boundary findings. The blood P and the revised classifier random-set P are read from their result tables. The other entries come from the named source analyses.

## ExtData_Fig2

- **a**: M-6 programme score (30 genes defined in RA synovial tissue) in synovial-fluid macrophage pseudobulk for each patient, treatment-naive and on treatment. Lines join the same patient, colour = drug and black bars are medians. The paired Wilcoxon test and the comparison with 200 expression-matched random gene sets come from ra_fluid_tests.tsv.
- **b**: Share of macrophages above the 90th percentile of treatment-naive M-6 cell scores, for each patient and time point.
- **c**: Tph programme score in T-cell pseudobulk. Eight patients have enough T cells at both time points.
- **d**: Top markers of the Leiden cluster with the highest M-6 score among 39,260 primary-tumour macrophages (LuCA core atlas, Harmony by study). Red marks genes in the RA M-6 programme. P compares this cluster with the best cluster of 200 random gene sets.
- **e**: Macrophage pseudobulk M-6 score for each donor in adjacent normal lung and primary tumour, for the 29 donors with >= 50 macrophages in both.
- **f**: Macrophage M-6 versus CD4 T-cell Tph programme score for each tumour, each centred within study. Spearman correlation on the study-adjusted residuals.
- **g**: Cox proportional hazards per SD of the bulk programme score in TCGA primary tumours, adjusted for age and stage, with and without adjustment for leukocyte abundance (PTPRC). Bars are 95% CI.
- **h**: Hypotheses as set in scripts/40_ra_fluid_treatment.py and 41_nsclc_extension.py before scoring. Verdicts use the stated two-sided tests.
