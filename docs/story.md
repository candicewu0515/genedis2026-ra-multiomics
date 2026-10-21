# RA genetic and macrophage-state evidence across public multi-omics resources

**GeNeDis 2026 · Xia Wu, PhD student, University of Iowa · xwu76@uiowa.edu**

## Summary

This workflow links complementary RA data layers and tests where candidate signatures fail. We reanalyse a published synovial RNA/ATAC atlas, map a sorted-macrophage discovery signature, compare bulk cohorts and assess spatial, protein and genetic evidence. CD4 genetic enrichment and IFN-associated macrophage features in tissue occur together. These data do not establish a directional causal chain from lymphocytes to macrophages.

Our contribution is a reproducible evidence map across resources, with integration and model baselines, clear score definitions and negative controls. The atlas states and their original genetic enrichment come from [Gupta et al., Nature Genetics 2023](https://doi.org/10.1038/s41588-023-01577-7). We keep the authors' fine-state annotations, including M-6 STAT1+ CXCL10+ macrophages. We do not present these states as new discoveries.

## Results

### Atlas and integration (Fig. 1)

The joint RNA/ATAC analysis contains 31,520 nuclei from 12 donors (11 RA, 1 OA). The author-defined M-6 state contains 65 nuclei. The benchmark compares RNA, ATAC, concatenation, WNN and one MOFA+ configuration. It uses shared cells and labels derived from both RNA and ATAC (Fig. 1e). Concatenation and WNN score best overall, but no method isolates M-6 as a cluster. Rare-state neighbourhood purity is highest with RNA alone. This evaluates only the tested configuration. It is not a general ranking of integration methods.

### Genetic association and tissue niches (Fig. 2)

Conditional S-LDSC links broad T-cell chromatin to RA heritability in Ishigaki EUR and FinnGen (tau z about 4.7). CD4 Tph/Treg/OX40-associated states are among the leading subtypes. T2D is a trait control. scDRS uses expression and the same GWAS traits. Its agreement with S-LDSC across states is modest (rho 0.26, P 0.047; Fig. 2d), and it is complementary methodological evidence. It has not been established that all GWAS participants are independent across resources.

Fine-mapping supports an ANKRD55-region variant across ancestries (Fig. 2e). This does not extend to all T-cell enrichment, because the East Asian analysis shows no broad T-cell enrichment. A non-positive myeloid conditional coefficient cannot establish that macrophages are downstream. M-6 is too small for its own S-LDSC annotation.

CXCL13-positive T cells have more B-cell neighbours in all four RA 5K sections (Fig. 2f–h). Each section reaches the empirical permutation resolution of about 0.001. Fig. 2i reports empirical P and does not extrapolate large z scores into extreme significance. Global label shuffling keeps cell positions but not label clustering or tissue compartments. This spatial association does not show that CXCL13 created the aggregates. CXCL13/CXCR5 communication ranks point to a compatible interaction for follow-up.

### Distinct macrophage measurements (Fig. 3)

**S1: sorted-cell discovery signature.** Fifteen genes are higher in active than in remission sorted macrophages at FDR <0.05. They include metallothioneins, glycolytic genes, APOL6 and GBP5 (Fig. 3a). Thirteen of these genes are measurable in the multiome comparison. S1 scores highest in M-6, and donor RNA and ATAC gene-activity scores agree (rho 0.84; random-set P 0.008; Fig. 3b–c). Both modalities measure the same nuclei. No S1 gene has an FDR-significant peak-gene link in this analysis. The independent multiome association with donor inflammation is not significant (10 RA donors, rho 0.43, P 0.21). In bulk RA/OA tissue, S1 fails its association and matched-random specificity tests (P 0.327; random-set P 0.525). When discovery donor SA139 is dropped, only 1/15 genes stays significant, although the positive directions remain. S1 remains exploratory.

**S2: predefined IFN programme.** CXCL9, CXCL10, CXCL11, STAT1, GBP1, GBP5, IDO1, IRF1, CD274 and APOL6 form a separate ten-gene score. S1 and S2 share only APOL6 and GBP5. S2 is associated with IFNG-positive lymphocyte neighbourhoods in 3/4 RA 5K sections (Fig. 3e). It is also associated with more IFN-high Visium spots in RA (Fig. 3f). In bulk tissue, its RA/OA difference goes beyond matched random sets (random-set P 0.002; P 0.001 after adjustment for macrophage abundance; Fig. 3h). These are S2 associations, not a successful S1 bulk validation. The limits of global spatial shuffling also apply here.

**Hallmark pathway/protein evidence.** Disease and treatment pathway NES values are inversely related (Fig. 3d). The protein panel uses Hallmark IFN-gamma response genes, not the exact S1/S2 set. RA and healthy donor scores differ in both lining and sublining (6 versus 3 donors per compartment, Mann-Whitney P 0.024 each; Fig. 3g). Fibrous and immune sublining regions are averaged within each donor before testing and plotting. Matched random protein sets support specificity in lining (0/1,000 as extreme; resolution about 0.001) but not in sublining (P 0.339). Xenium, Visium and proteomics share blocks and cohort, so they are not three independent patient replications.

**Treatment.** In the anti-TNF cohort, 13/15 S1 genes are measured. The technically residualised S1 score decreases within good responders (median change -0.70, paired P 0.014; Fig. 3i). However, the change does not differ significantly between good responders and nonresponders (P 0.131). This does not establish response specificity or prediction. Apparent baseline prediction is confounded by sequencing depth.

### Candidate regulators and targets (Fig. 4)

RNA regulon and motif comparisons point to NF-kB, STAT1 and IRF-related candidates (Fig. 4a–b). The M-6-bearing units come from two donors. The candidate flag combines BH-adjusted gene-label-permutation q <0.05 with a top-decile motif ranking. Gene permutations do not give donor replication or causal evidence.

Two IRF5 GWAS signals colocalise with myeloid expression QTLs (Fig. 4c). Under the model assumptions, this supports shared association signals. Single-instrument Wald ratios link higher CTLA4 expression to lower risk and higher IL6ST expression to higher risk (Fig. 4d). Pleiotropy cannot be tested. These estimates do not recover drug efficacy. [Abatacept binds CD80/CD86](https://www.ema.europa.eu/en/medicines/human/EPAR/orencia) and [tocilizumab binds IL6R](https://www.ema.europa.eu/en/medicines/human/EPAR/roactemra), while IL6ST encodes gp130. IL6ST and ANKRD55 share an instrument and cannot be separated here. IRF5 directions conflict across datasets and signals. Drug annotations give biological context only.

The transparent 11-layer gene ranking (Fig. 4f–g) sets priorities among follow-up hypotheses. The layers overlap, and the scores are not causal probabilities. A low Open Targets association describes that database release only. It does not prove a lack of literature or guarantee novelty.

### Extension and boundaries (Extended Data Figs. 1–2)

The extensions use separate 30-gene M-6/Tph programmes derived from the author-defined atlas states. In synovial fluid (9 paired patients), the M-6-like fraction decreases in 7/9 (P 0.027). The paired change in M-6 pseudobulk score is not significant (P 0.129), even though the matched-random comparison gives P 0.005. Cluster evidence is borderline (P 0.055), response labels are not available, and two treatments are pooled.

In LuCA, macrophage M-6 scores are higher in tumour than in adjacent lung in 20/29 paired donors (P 0.0038). The macrophage/CD4 Tph correlation is modest and not conventionally significant (rho 0.24, P 0.092), although the random-set P is 0.030. The TCGA primary survival tests adjusted for age and stage are negative. The leukocyte-adjusted LUAD result is a sensitivity analysis. It can be affected by tumour-cell expression and by multiple testing. These are cross-disease associations. They do not show a shared treatment mechanism.

Classifier evaluation uses training-only repeated nested CV (25 outer folds) and permutations that repeat the same procedure. In RA versus OA (174 biopsies), the eight-programme model has a mean AUC of 0.900, against 0.979 for the genome-wide baseline. The approximate corrected resampled comparison favours the genome-wide model (P=0.00596). In ACPA-positive versus ACPA-negative RA (139 biopsies), the matching means are 0.797 and 0.726. The approximate corrected comparison finds no significant difference (P=0.149). The prioritised-gene model does not do significantly better than 100 random gene sets matched for feature count, in either task (empirical P=0.2475 and 0.5149). Both programme-model label tests reach P=0.005, the smallest value possible with 199 permutations. This shows internal association under this procedure, not external prediction or biological specificity. Full results are in QC report section 33 and classifier_baselines.tsv. These CV and permutation estimates are conditional on the given feature lists. They do not correct for earlier exploratory choices of tasks or programmes. These are exploratory tests for each task. The SD across correlated folds is not a confidence interval, and no independent clinical model validation is available. Blood transfer, deconvolution sanity, L1000 specificity and chromatin-before-RNA timing fail or remain uninterpretable. They are kept in Extended Data Fig. 1.

## Figure legends

Panel letters refer to the final figures. Gene sets and measurement units are kept separate.

### Fig1 | Study design, published atlas states and integration baselines

- **a**, Dataset inventory, with counts from the existing study description.
- **b**, WNN coordinates from the reanalysis. Colour and labelled counts show eight broad cell types. Black points mark M-6.
- **c**, Nuclear composition by donor, with yield printed above each bar. Squares below mark RA (vermilion) and OA (grey).
- **d**, Subtype RNA weights. Points are subtype summaries, black ticks are cell-type medians and ranges show the subtype minimum to maximum.
- **e**, Raw 0–1 benchmark metrics, with the same colour scale in every column (white = 0, saturated blue = 1). The overall score is as supplied by the analysis.
- **f**, Percent variance explained by each of the first eight MOFA+ factors.
- **g**, Donor median genes versus peaks, with point area proportional to retained nuclei. QC RNA and ATAC counts may differ slightly.

### Fig2 | Genetic association and spatial neighbourhoods

- **a**, Colour shows S-LDSC tau z, and area shows nominal −log10 P capped at 6. A black ring marks FDR < 0.05.
- **b**, Top 10 mean-RA S-LDSC subtype statistics in the two European RA cohorts.
- **c**, Fine-mapped PIP mass relative to the permutation null. Shape and colour identify the GWAS, area shows nominal P and a black outline marks FDR < 0.05.
- **d**, Mean RA S-LDSC versus mean RA scDRS for states with both estimates, with Spearman correlation across states. M-6 has no S-LDSC estimate and is excluded.
- **e**, Top variants, selected by maximum supplied PIP and unique nearest protein-coding gene. Labels include the peak state. Missing study records are grey dashes, not zero PIP. Exact variant IDs are in the source table that goes with the figure.
- **f**, Paired B-cell neighbourhood means within each section, with positive/negative fold ratios. Lines connect the same section.
- **g**, Xenium 5K section CT_21045-A-15, supplied as a processed plotting table generated by step 20. The outlined square is the ROI in h. CXCL13 positivity means a raw count >0.
- **h**, A 1.2-mm ROI centred on the densest 500-µm B-cell bin. It is selected reproducibly to show aggregate architecture, not to quantify enrichment. Same cells and coordinates as g.
- **i**, Section-level empirical permutation P, capped at the 1/1001 resolution of the 1,000 permutations. Global label shuffling does not keep spatial label clustering or tissue compartments. These associations do not establish causal organisation.

### Fig3 | Distinct discovery, IFN-programme and protein measurements

- **a**, Original DESeq2 active/remission effects and Wald statistics. P is recomputed from the stored Wald statistic to avoid P values rounded to zero. FDR status uses the supplied adjusted P.
- **b**, S1: 15-gene discovery definition. The 13 measurable genes are used in the RNA and ATAC gene-activity z scores shown for all 14 myeloid states, on a common centred scale.
- **c**, S1: donor-level RNA/ATAC gene-activity scores for the discovery signature (13 of 15 genes measurable), with Spearman correlation and the supplied random-gene-set test. The line is a descriptive least-squares fit.
- **d**, Top eight disease Hallmark NES values, with treatment NES for the same pathways; * supplied GSEA FDR < 0.05. Disease = active/remission macrophages; treatment = post/pre bulk in good responders.
- **e**, Paired section-level regression coefficients for IFNG-positive and IFNG-negative T/NK neighbours. P values are one-sided spatial permutation tests from the existing table.
- **f**, S2 is the predefined 10-gene IFN programme. Each dot is a section, and black bars are medians. The comparison shown is the supplied donor-level RA/healthy test, so section counts are not independent donor counts.
- **g**, Hallmark IFN-gamma protein score. Sublining combines immune_sublining, fibrous_sublining and sublining before averaging within each donor, as in the statistical test (RA n=6, NH n=3 in each compartment). Black bars are medians. Two-sided Mann-Whitney P values are recomputed and checked against the source table. Protein-set specificity is supported in lining, not in sublining (random-set P=0.339).
- **h**, S2: predefined 10-gene IFN programme, separate from S1. Individual biopsy scores, with box IQR, median and 1.5×IQR whiskers. Counts are biopsies. The supplied random-gene-set comparison is shown.
- **i**, S1: 15-gene discovery signature (13 measured in this cohort). Individual paired patient changes after technical residualisation, with box IQR, median and 1.5×IQR whiskers. The P values shown are supplied paired Wilcoxon tests, not comparisons between response groups.

### Fig4 | Candidate regulatory and expression-associated genetic evidence

- **a**, Transcription-factor RNA regulon and chromVAR motif differences. Candidates need BH q <0.05 across all tested RNA gene-label-permutation P values plus a top-decile motif ranking. Data come from two donors, and gene-label permutations are not donor replication. These are candidate regulators.
- **b**, Top 10 summed RNA/motif ranks. Point area = 15 + 2 sqrt(regulon target count). Each colour column uses its labelled scale in its own units (RNA 0–7, ATAC 0–6). A black outline marks RNA gene-label-permutation BH q <0.05. Numbers are unscaled effects.
- **c**, Maximum IRF5 H4 within each project/biosample/condition and RA lead signal. Top 10 conditions by maximum H4. Grey dashes mark missing records.
- **d**, Existing MR odds ratios and supplied 95% CI for six genes. All supplied estimates are shown, including several lead signals per gene/GWAS combination. Bracketed indices tell these apart in source-table order. The right track gives same-direction/total datasets. Colour shows direction, and grey marks agreement <75%. The estimates are exploratory, not drug efficacy.
- **e**, Fraction of immune QTL datasets in which the risk allele raises expression. Area and the right annotation give the number of datasets. A black ring marks the supplied consistency flag.
- **f**, Top 24 genes and 11 original evidence layers on a common 0–1 scale. The original evidence score is shown above. White = 0, intermediate blue = 0.5, saturated blue = 1. A green marker below flags low prior RA evidence, based on the supplied novelty label. Evidence layers are not assumed to be statistically independent.
- **g**, All prioritised genes: original integrated score versus Open Targets association score. Point area = 5 + 2 × evidence-layer count. Green uses the supplied low-prior-evidence label.

### ExtData_Fig1 | Baselines, robustness and boundary findings

- **a**, Revised training-only repeated nested-CV AUC, mean ± fold SD (not a CI). Random-control error is the SD across draw means. Observed and label-permuted F1 use the same 5 × 5 outer CV, 3-fold tuning and preprocessing. These are internal estimates from one cohort.
- **b**, Baseline good-versus-nonresponse AUC with the supplied bootstrap 95% CI. The observed baseline predictive signal is technically confounded.
- **c**, Fraction meeting each original perturbation rule. The right column gives the number of perturbations. The underlying metric definitions and ranges are kept in robustness_summary.tsv.
- **d**, Observed versus random-median share of anti-inflammatory compound classes. P values are supplied empirical null comparisons.
- **e**, Existing RNA and ATAC PAGA connectivity to M-6. Graph proximity does not establish developmental or temporal direction.
- **f**, Boundary findings. Blood P and the revised classifier random-set P are read from the matching result tables. Other boundaries refer to the named source analyses.

### ExtData_Fig2 | Exploratory synovial-fluid and NSCLC extensions

- **a**, Patient-level synovial-fluid macrophage pseudobulk M-6 programme score (30 genes defined in RA synovial tissue), treatment-naive and on treatment. Lines join patients, colour = drug and black bars are medians. Paired Wilcoxon test and comparison with 200 expression-matched random gene sets, from ra_fluid_tests.tsv.
- **b**, Share of macrophages above the 90th percentile of treatment-naive M-6 cell scores, per patient and time point.
- **c**, T-cell pseudobulk Tph programme score. Eight patients have enough T cells at both time points.
- **d**, Top markers of the Leiden cluster with the highest M-6 score among 39,260 primary-tumour macrophages (LuCA core atlas, Harmony by study). Red marks genes in the RA M-6 programme. P compares the cluster with the best cluster of 200 random gene sets.
- **e**, Donor-level macrophage pseudobulk M-6 score in adjacent normal lung and primary tumour for the 29 donors with >= 50 macrophages in both.
- **f**, Tumour-level macrophage M-6 versus CD4 T-cell Tph programme scores, each centred within study. Spearman correlation on study-adjusted residuals.
- **g**, Cox proportional hazards per SD of the bulk programme score in TCGA primary tumours, adjusted for age and stage, with and without leukocyte abundance (PTPRC). Bars are 95% CI.
- **h**, Hypotheses as fixed in scripts/40_ra_fluid_treatment.py and 41_nsclc_extension.py before scoring. Verdicts use the stated two-sided tests.

## Reproducibility and interpretation

README.md gives the commands. docs/REVISION_VALIDATION.md records the final clean run and package verification. docs/evidence_map.tsv lists the exact score definitions. Analysis script headers state the questions, and this work was not preregistered. The figures are data plots with supplied source hashes and panel notes.
