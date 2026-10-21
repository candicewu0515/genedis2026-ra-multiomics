# QC and model report: RA synovium single-cell module

**Reading guide:** This is the retained detailed analysis history. Read `docs/methods_and_QC.md` or `presentation/GeNeDis2026_methods_QC.pdf` for the current consolidated methods/QC appendix, and `presentation/GeNeDis2026_submission_report.pdf` for the competition technical report.
Generated 2026-10-18 and extended 2026-10-19/20 (sections 19–35) from the pipeline in `scripts/`; logs are in `logs/`.
Section 36 item 16 describes an earlier run. The current revision and a fresh validation are recorded in docs/REVISION_VALIDATION.md. Where that record differs from the older runtime or classifier statements below, the record takes precedence. The current figure set is the Nature-inspired redesign; older diagnostic images have been removed. The integrated narrative and main figures are in `docs/story.md` and `results/figures/main/`.
Each result has one of three status tags: **ready** (analysis available; not a claim of clinical validity), **exploratory** (direction-consistent but
not significant, or confounded), **unavailable** (not done, or not possible with these data).

## 1. GSE243917: whole synovial tissue, snRNA (Multiome RNA)

### QC (`results/tables/gse243917_qc_per_sample.tsv`, `results/figures/main/` (current figures; panel map in `Figure_notes.md`))

- 31,547 nuclei x 36,601 genes as deposited (the authors had already removed low-quality nuclei).
  Floor QC (>= 200 genes, < 20 % mito) removed 0 nuclei. Scrublet flagged 27 (0.09 %), and these were removed.
- Two inferred nuclei-prep batches: median mito 0.5–1.1 % (BRI-1586, -1588, -1590, -1592, -1594 = "lowMT")
  vs 6.7–11.3 % (the other seven). Three lowMT samples (1590, 1592, 1594) captured very few fibroblasts.
  Batch is kept as a covariate, and every feature is tested against it.
- BRI-1588 gave only 137 nuclei (4 macrophages). It is excluded wherever a per-donor pseudobulk is needed.

### Annotation (**ready**)

- Harmony (donor) + Leiden. Broad types were assigned from marker scores, without using author labels.
- Agreement with author cell types: **ARI 0.952, 97.3 % of nuclei identical**
  (`results/tables/gse243917_annotation_vs_author.tsv`). Main disagreements: ~300 author-NK nuclei called T
  (NK/T boundary) and ~300 author-mural nuclei called fibroblast. The result did not change after the marker sets were revised (section 4).

### Composition vs histology (**exploratory**; `results/tables/gse243917_feature_vs_histology.tsv`)

- n = 11 RA donors. 8 of 11 share lymphocytic-inflammation grade 2, so power is low.
- No feature passes FDR < 0.05 (all FDR >= 0.86). The trends point the expected way. The myeloid share of immune cells
  falls with inflammation (Spearman rho = -0.56, p = 0.075), while plasma cells (rho = 0.51, p = 0.11) and T cells
  (rho = 0.49, p = 0.13) rise. This is the expected lymphoid-rich pattern.
- Macrophage-state features (inflammatory share, homeostasis index, homeostatic score) differ by prep batch
  (Mann-Whitney p = 0.004–0.009; MERTK+ subset share p = 0.052). They are therefore **not interpretable** in this
  cohort without a batch-balanced design.

## 2. E-MTAB-8322: sorted synovial macrophages, scRNA

### QC (`results/tables/emtab8322_qc_per_sample.tsv`)

- 99,512 barcodes. Keeping cells with >= 300 genes and < 15 % mito left 81,601 cells after Scrublet doublet removal.
- SA149P was merged with SA149 (same patient) -> 26 donors: 4 healthy, 3 UPA, 6 naive-active RA,
  6 treatment-resistant active RA, 7 remission RA.

### Non-macrophage cells in a "macrophage" dataset

- The Atlas labels every barcode "macrophage" at sample level. But stage-1 clustering finds fibroblasts,
  T/NK, B/plasma, endothelial, mural and mast cells: **only 42,963 / 81,601 cells (52.7 %) are myeloid**.
  The contamination is mostly in extended-sort-panel samples (myeloid share per sample: median 97 %, min 6 %;
  `results/tables/emtab8322_broad_composition.tsv`). All downstream analysis uses myeloid cells only, and the
  pseudobulk model adjusts for sort panel.

### Macrophage states (`results/tables/emtab8322_clusters.tsv`)

12 myeloid subclusters. 3 DC-like clusters (CLEC10A/FCER1A/HLA-class II) are left out of the macrophage
features. Homeostatic clusters: C0 SELENOP+ RNASE1+ LGMN+ (LYVE1-high), C1 NUPR1+ MARCO+, C2 CCL3+ CCL4+,
C3 SELENOP+ C1QC+ FOLR2+. Inflammatory clusters: C4 SPP1+ FN1+ FABP5+, C5 FCN1+ LYZ+ S100A12+ (monocyte-like),
C9 GBP1+ CXCL10+ (IFN), C11 APOC1+.

### Donor-level group tests (`results/tables/emtab8322_group_tests.tsv`)

| Feature | Healthy | RA naive-active | RA resistant-active | RA remission | Active vs healthy (FDR) | Active vs remission (FDR) |
|---|---|---|---|---|---|---|
| Inflammatory share of macrophages | 0.09 | 0.42 | 0.51 | 0.32 | 0.028 | 0.64 |
| Homeostatic program score | 0.43 | 0.03 | -0.08 | 0.19 | 0.028 | 0.64 |
| SPP1+ FN1+ cluster (C4) share | 0.03 | 0.10 | 0.12 | 0.04 | 0.028 | 0.64 (raw p 0.056) |
| Homeostasis index, log2(homeostatic/inflammatory) | 3.40 | 0.47 | -0.05 | 1.08 | 0.028 | 0.64 |

- Active RA vs healthy: a **directionally consistent reanalysis** of Alivernini et al. 2020 (loss of resident homeostatic
  macrophages, gain of inflammatory/monocyte-derived states), **but confounded**. The healthy donors are younger
  (22–48 y) and 3 of 4 are male, while the RA donors are 42–77 y and all female.
- Active vs remission (**exploratory**): remission sits between healthy and active on most features. But no
  difference survives FDR with 12 vs 7 donors.

### Pseudobulk DE, active vs remission, macrophages only (**exploratory discovery**)

`results/tables/emtab8322_DE_active_vs_remission.tsv`; PyDESeq2, `~ sort_panel + condition`, 12 vs 7 donors.
18 genes at FDR < 0.05 (17 up in active): metallothioneins (MT1H, MT1G, MT2A, MT1F, MT1X, MT1E, MT1M, MT1L),
glycolysis (PKM, PFKP) and the IFN-gamma/STAT1 axis (GBP5, APOL6, SLAMF7, IL4I1; STAT1 and ANKRD22 at FDR 0.088,
CD40 at 0.10). This fits an IFN-gamma-activated, glycolytic inflammatory macrophage state in active disease.
Caveat: metallothioneins also respond to dissociation stress. The available analysis cannot separate out the stress-related part.

## 3. Models

### Signature transfer to an independent cohort (**exploratory**)

`results/tables/cross_dataset_signature_stats.tsv`. The 15 FDR genes up in active RA (E-MTAB-8322, scRNA)
were scored in GSE243917 macrophage pseudobulk (snRNA, different donors, different lab). Only donors with >= 50
macrophages were used (n = 10 RA).

| Signature | Genes found | rho vs inflammation grade | p | Batch p |
|---|---|---|---|---|
| Primary: FDR genes up in active (15) | 13/15 | 0.43 | 0.21 | 0.26 |
| Sensitivity: nominal top-50 each direction | 48/50 + 47/50 | 0.32 | 0.38 | 0.76 |

The direction is consistent (more inflamed tissue -> higher active-RA macrophage score), but the result is not significant.

### Supervised baseline: active vs remission (**exploratory (negative)**)

`results/tables/classifier_loocv.tsv`. Donor-level, leave-one-donor-out, standardisation inside each fold,
L2 logistic regression; label-permutation null (500 permutations).

| Model | AUC | Balanced accuracy | Permutation p |
|---|---|---|---|
| Majority class baseline | 0.50 | 0.50 | n/a |
| Single feature: homeostatic program score | 0.60 | 0.60 | 0.12 |
| Three a priori features | 0.50 | 0.56 | 0.25 |

Macrophage-composition features alone do not separate active from remission in 19 donors. They do not beat
the baseline.

**Leakage controls.** Labels: active = naive-active + treatment-resistant active; remission = sustained
remission (as annotated by the depositors). Features are gene-program scores fixed from the literature and
label-free cluster assignments. No feature was selected using the labels it predicts. One unit is one donor, so
no donor's cells appear in both training and test folds. Caveat: the feature list was fixed after the group
tests had been inspected, so the classifier is exploratory.

## 4. Chromatin layer and RNA+ATAC integration (GSE243917, same nuclei)

### ATAC QC (`results/tables/atac_qc_per_sample.tsv`, `results/figures/main/` (current figures; panel map in `Figure_notes.md`)) (**ready**)

- 132,520 fixed-width (200 bp) peaks. All 31,520 RNA-retained nuclei have ATAC. Median accessible peaks per
  nucleus: 5,278–11,474 by sample. Every nucleus passes the >= 1,000-peak floor. Peaks open in < 0.5 % of nuclei
  were dropped (132,479 kept).
- TF-IDF + LSI on the 50,000 most frequently open peaks. LSI component 1 tracks depth (|r| = 0.99), so it is
  dropped. Harmony by donor.

### What chromatin alone recovers (`results/tables/atac_label_transfer.tsv`) (**ready**)

- kNN label transfer in ATAC space recovers author cell types for **97.8 %** of nuclei (mural 59 %, NK 80 %,
  all others >= 98 %), but author subtypes for only **60.3 %**. Best-recovered subtypes: mural (96 %), IgG plasma
  (95 %), PRG4+ CLIC5+ lining fibroblasts (85 %), CXCL12+ sublining fibroblasts (83 %). Small or
  activation-defined states are poorly recovered (proliferating T 1 %, STAT1+ CXCL10+ macrophages 3 %). kNN
  majority voting also penalises small groups.

### Joint embedding (WNN) (`results/tables/wnn_recovery.tsv`, `wnn_modality_weights.tsv`) (**ready**)

| Graph | Neighbour agreement, cell type | Neighbour agreement, subtype | Clustering ARI vs subtypes |
|---|---|---|---|
| RNA | 0.964 | 0.515 | 0.382 (Leiden, step 02) |
| ATAC | 0.965 | 0.450 | n/a |
| WNN (RNA+ATAC) | 0.968 | 0.503 | 0.414 |

Per-nucleus modality weights fall in a narrow range (median RNA weight 0.46–0.52 by subtype). Fibroblast subtypes rely most
on chromatin (0.46–0.47). Proliferating and dendritic-cell states rely more on RNA (0.50–0.52).

### RNA vs gene activity and peak-gene links (**exploratory**)

- Gene activity = ATAC counts in gene body + 2 kb upstream (Ensembl 98). Across 238 pseudobulk units
  (subtype x donor, >= 30 nuclei), the median per-gene RNA-activity correlation is r = 0.21 (8,979 expressed
  genes; 22 % with r > 0.5). This is in the usual range for multiome data.
- Peak-gene links (+/-250 kb, 421,253 cis pairs, 9,104 genes): the trans-chromosome null is wide (99th
  percentile r = 0.60), because units differ mainly by cell type. So only **277 links (151 genes)** pass
  r > 0.45 and FDR < 0.05 (median distance to TSS 15 kb). No signature gene has a link at this threshold.

### The active-RA macrophage signature on both layers (**ready (descriptive) / exploratory (clinical)**)

`results/tables/signature_chromatin.tsv`, `signature_by_myeloid_subset.tsv`, `multiome_signature_stats.tsv`,
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). Signature = 15 FDR genes up in active RA macrophages from
E-MTAB-8322 (13 present in the Ensembl gene set used here).

- **Which cell state:** in GSE243917, on RNA, the signature is highest in the **M-6 STAT1+ CXCL10+ (IFN-activated)**
  macrophage subset (mean z 1.49; next highest 0.34). It is lowest in resident M-0 MERTK+ SELENOP+ LYVE1+
  (-0.36). On chromatin the ranking is less specific (top: M-14 LAMP3+ DC 0.88, M-3 MERTK+ HBEGF+ 0.42,
  M-4 SPP1+ 0.41, M-6 0.33). Across myeloid subsets, RNA and gene activity agree significantly for 7 of 13 genes (Spearman
  rho 0.55–0.90, p < 0.05: MT1H, MT1G, MT2A, PKM, IL4I1, SLAMF7, GBP5); APOL6 rho 0.50 (p = 0.07).
- **Donor level (RA, >= 50 macrophages, n = 10):** ATAC gene-activity signature score vs RNA signature score,
  **rho = 0.84 (p = 0.002)**. Control: 1,000 random gene sets matched on expression decile give median rho 0.19
  (95th percentile 0.69); **empirical p = 0.008**. So the agreement is specific to the signature and is not a global
  nucleus-quality effect. Neither layer is associated with batch (p = 0.26 RNA, 0.61 ATAC).
- vs inflammation grade: RNA rho 0.43 (p = 0.21), ATAC rho 0.27 (p = 0.45). The direction is consistent, but neither is significant.

Interpretation: the IFN-gamma / metallothionein macrophage programme was found in active RA in the scRNA cohort. Its
donor-to-donor variation is also present in an independent snMultiome cohort. There it is concentrated in STAT1+ CXCL10+ macrophages and
mirrored in chromatin accessibility of the same genes. A larger, clinically labelled cohort is needed to test whether it
tracks clinical activity.

## 5. Genetics layer: RA risk loci near cell-state-specific chromatin

`scripts/09_gwas_download.py`, `scripts/10_gwas_peak_enrichment.py`; tables `results/tables/gwas_*.tsv`;
figure `results/figures/main/` (current figures; panel map in `Figure_notes.md`).

### Inputs (**ready**)

- GWAS Catalog REST API v2 (the v1 API is no longer available), queried 2026-10-18, GRCh38. RA: 3,221 associations from
  95 studies. Associations from cross-disease scans (e.g. hypothyroidism + RA, 1,035 rows; RA + T1D, IBD, MS)
  were removed by a rule fixed before any enrichment was run. The rule keeps rows whose mapped traits are all RA terms
  (RA, ACPA+ RA, RF/ACPA seropositivity). After p <= 5e-8 and MHC exclusion: **494 lead SNPs from 40 studies ->
  240 loci at +/-25 kb** (top loci PTPN22, CCR6, ANKRD55, TNFAIP3, PADI4). Negative control, LDL cholesterol:
  4,302 lead SNPs -> 1,210 loci.
- Peak sets: top 4,000 specific peaks per state (8 cell types vs all others; 57 subtypes with >= 100 nuclei vs
  their own lineage). Universe: 131,644 non-MHC peaks. Null: 1,000 random sets matched on accessibility decile x
  promoter/distal. The statistic counts loci, not peaks, so clustered peaks at one locus count once.

### Results (`results/tables/gwas_enrichment.tsv`)

| State (cell type vs all) | RA loci hit / null mean | Enrichment | FDR, 10 kb | FDR, 25 kb | FDR, 100 kb |
|---|---|---|---|---|---|
| T cell | 76 / 47.9 | 1.59 | 0.004 | **0.008** | 1.00 |
| B cell | 59 / 47.3 | 1.25 | 0.004 | 0.084 | 1.00 |
| NK | 58 / 47.8 | 1.21 | 0.012 | 0.115 | 1.00 |
| Plasma | 51 / 46.1 | 1.11 | 0.005 | 0.45 | 1.00 |
| Myeloid | 39 / 47.7 | 0.82 | n/a | 1.00 | n/a |
| Fibroblast | 20 / 47.9 | 0.42 | n/a | 1.00 | n/a |
| Mural | 18 / 48.5 | 0.37 | n/a | 1.00 | n/a |

- **Lymphocyte enrichment (ready).** RA risk loci sit near T-, B- and NK-specific chromatin and away from
  stromal chromatin. This reproduces the established T-cell-centred architecture of RA heritability. At 100 kb
  almost every locus contains some peak (210/217), so the test saturates.
- **Control (ready).** LDL loci show no enrichment in any cell type or subtype (all FDR = 1; enrichment 0.58–1.11).
- **Subtypes within lineage (exploratory).** Top: T-10 CD4+ OX40+ NR3C1+ (1.43x), NK-1 CD56dim IFNG+ CD160+
  (1.42x), T-9 CD25-low Treg (1.31x), with FDR 0.057–0.095 at 25 kb and 0.028–0.10 at 10 kb. T-8 CD25-high Treg has FDR 0.038
  at 10 kb. Tph subsets (T-3, T-7) are not enriched at locus level (1.04–1.07x). Resident MERTK+ macrophage states
  are depleted (M-0 0.74x, M-3 0.67x, M-2 0.52x).
- **The STAT1+ CXCL10+ macrophage state that carries the active-RA signature (section 4) has 65 nuclei. This is
  below the 100-nucleus minimum, so its genetic enrichment cannot be tested in this dataset.**

### Candidate regulatory variants (**exploratory**; `results/tables/gwas_candidate_variants.tsv`)

37 RA lead SNPs lie within 250 bp of a state-specific peak (108 SNP-state pairs). Examples with a coherent
biological reading: rs3087243 (CTLA4) in Treg and Tph peaks; rs117605225 (CXCL13, the hallmark Tph chemokine)
and rs7731626 (ANKRD55) in Tph/Tfh peaks; rs9927316 (IRF8 region) in IFNG+ NK and myeloid peaks; rs1883832 (CD40)
in naive CD4 peaks; rs3807306 (IRF5) in B/plasma peaks. None of these peaks has a significant peak-gene link at
the strict threshold of section 4, so target genes are nearest-gene assignments. Lead SNPs are not necessarily
causal. Fine-mapped credible sets would be needed to prioritise variants.

## 6. Clinical-response test: anti-TNF cohort (GSE198520)

`scripts/11_gse198520_response.py`; tables `results/tables/gse198520_tests.tsv`, `gse198520_auc.tsv`;
figure `results/figures/main/` (current figures; panel map in `Figure_notes.md`).

**Cohort choice.** R4RA (E-MTAB-11611; rituximab vs tocilizumab) was the first choice, but it cannot be used. Only raw
FASTQ files are public, the authors' scripts load internal processed files, and the public sample table has CDAI50 response
but no randomised arm. GSE198520 (Wang et al. 2022, *Arthritis Rheumatol*) provides raw gene counts for 46 RA patients.
They were biopsied before and 12 weeks after anti-TNF (certolizumab 27, etanercept 19), and EULAR response is given (good 19, moderate
13, non 14). Sequencing depth is low (median 0.84 M reads/sample), and some genes are not covered (e.g. CD68).
13 of the 15 signature genes are measured (the same 13 as in GSE243917).

**Declared hypotheses** (in the script header, written before outcome data were analysed): H1, baseline score vs
disease activity; H2, baseline score vs response (good vs non; AUC vs majority class, baseline DAS28 and random gene
sets); H3, on-treatment change by response group.

### Technical confounding (**ready (data-quality finding)**)

- Baseline PC1 explains 57 % of expression variance. It tracks log library size (rho -0.84), mitochondrial-read
  fraction (-0.81), RIN (-0.64), EULAR response (-0.64) and the signature score (-0.71).
- **Library size alone separates good from non-responders with AUC 0.84** (mito fraction 0.80). Random
  expression-matched gene sets give AUC ~0.15, because a global shift lowers most genes in good responders.
- The signature genes depend unusually strongly on library size (mean per-gene rho 0.42 vs -0.24 for matched random
  genes, p = 0.0005). Library size + mito fraction explain 43 % of signature-score variance.
- Library size also tracks the histological inflammatory score (rho 0.46, p = 0.002). So part of this axis is
  probably biological (more cellular, inflamed tissue yields more RNA). The technical adjustment below is therefore
  conservative and may remove real signal. The two cannot be separated in this cohort.

### Results

| Test | Raw score | After technical adjustment |
|---|---|---|
| H2 good vs non, AUC (n = 33) | 0.73 (95 % CI 0.53–0.89; permutation p 0.026) | 0.52 (p 0.90) |
| H2 logistic, OR per SD (adjusted for biologic / library size + mito) | 4.8 (p 0.019) | 1.14 (p 0.87) |
| H2 ordinal Spearman (non < moderate < good, n = 46) | rho 0.35 (p 0.018) | rho 0.02 (p 0.88) |
| H1 vs histological inflammatory score (n = 42) | rho 0.40 (p 0.009) | rho 0.24 (p 0.12) |
| H1 vs baseline DAS28 / CRP | rho -0.19 / -0.14 (ns) | n/a |
| H3 change pre -> week 12, good responders | median -0.60 (p 0.002) | median -0.70 (p 0.014) |
| H3 change, non-responders | median -0.19 (p 0.09) | median -0.27 (p 0.10) |
| H3 change, good vs non | p 0.025 | p 0.13 |

Comparators for H2: baseline DAS28 AUC 0.29 (good responders start lower, as EULAR "good" requires reaching DAS28
<= 3.2); majority class 0.50. Library size did not change significantly with treatment in either group. So the on-treatment
fall in the signature is not a library-size artefact.

**Interpretation: negative for prediction, exploratory for monitoring.** In this cohort the active-RA macrophage
signature does not predict anti-TNF response beyond a sample-quality / cellularity axis, and that axis itself
predicts response. The signature does fall after anti-TNF, most clearly in good responders. This fits with it marking
treatment-sensitive inflammation, but the between-group difference does not survive adjustment. It also means that
expression-based response predictors built on GSE198520 should be checked against library size.

## 7. Pathways (GSEA) (**ready**)

`scripts/12_pathway_gsea.py`; `results/tables/gsea_results.tsv`, `gse198520_DE_post_vs_pre_good.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). Pre-ranked GSEA (gseapy, 1,000 permutations, set size 15–500), MSigDB v2024.1.Hs.

- A: active RA vs remission macrophages (E-MTAB-8322, Wald statistic). B: anti-TNF good responders, week 12 vs
  baseline (GSE198520, paired PyDESeq2 `~ patient + timepoint`, 19 pairs, 85 genes FDR < 0.05; library size did not
  change in this group).
- Hallmark: A has 20 pathways up (FDR < 0.05) and none down. B has 31 down and 2 up. 17 are up in A and down in B.
  NES across all 49 pathways correlate negatively (Spearman rho = -0.61, p = 3e-6).
- Top shared: IFN-gamma response (A NES +2.67, B -3.21), IFN-alpha response, IL6-JAK-STAT3, hypoxia, mTORC1, unfolded
  protein response, PI3K-AKT-mTOR, IL2-STAT5, complement, ROS, inflammatory response, allograft rejection. Glycolysis
  is up in A (+2.17) but not significantly down in B.
- Reactome: 133 up in A, 393 down in B, 100 shared (rho -0.37 across 872). Proteasome / antigen-processing terms lead the list,
  in line with an IFN-gamma-induced immunoproteasome.
- Caveat: pathways share genes, so the p-value of the NES correlation is optimistic.

## 8. Cell-cell communication (LIANA) (**exploratory**)

`scripts/13_cell_communication.py`; `results/tables/liana_results.tsv.gz`, `liana_to_M6_specific.tsv`,
`ifng_senders.tsv`, `ifng_donor_link.tsv`, `cell_communication_tests.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
LIANA rank_aggregate (consensus resource, 1,000 permutations) over 68 author subtypes with >= 30 nuclei (31,082
nuclei; 968,510 scored interactions).

- By magnitude alone, the signals received by M-6 are dominated by abundant fibroblast ECM ligands to CD44 (FN1, PRG4,
  collagens). These are uninformative. M-6-specific received signals include autocrine CXCL10 -> TLR4 and CXCL9 -> FCGR2A,
  venular CCL14 -> CCR1/CCR5 (monocyte recruitment), VWF -> SIRPA and fibroblast BGN -> TLR2 (damage signal).
- IFNG -> IFNGR1/IFNGR2 on M-6: the strongest and most specific sender is NK-1 CD56dim CD16+ IFNG+ CD160+
  (lr_mean 0.85, specificity rank 0.004). This is the NK state that is also enriched for RA risk loci (section 5, 1.42x).
  IFNG is low in nuclear RNA, so this axis ranks mid-table on magnitude (0.53–0.62).
- M-6 chemokines: CXCL9/10/11 rank highest towards ACKR1 on venular endothelium (magnitude rank 0.001). This means chemokine
  presentation for leukocyte recruitment. Towards CXCR3 the top targets are CD4+ Tph, Vdelta1 and CD25-high Treg
  (magnitude ~0.38–0.40).
- Negative tests: IFNG expression does not track RA GWAS enrichment across 23 T/NK states (rho 0.06, p 0.78). Across
  RA donors, lymphocyte IFNG does not track the macrophage IFN-gamma programme (rho 0.10, p 0.78, n = 10) or the M-6
  share (rho 0.00).
- Ligand-receptor co-expression points to possible communication. It does not prove it.

## 9. Transcription-factor motifs (chromVAR) (**sanity ready; biology exploratory**)

`scripts/14_tf_motif_chromvar.py`; `results/tables/chromvar_*.tsv*`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
pychromvar on 514 pseudobulk units (subtype x donor with >= 10 nuclei, plus 11 donor macrophage pools), 132,479 peaks,
UCSC hg38 sequence, 1,019 JASPAR 2026 CORE vertebrate motifs and 50 GC/accessibility-matched background iterations.
Single-nucleus chromVAR was not feasible (dense 31k x 132k matrix > 16 GB).

- Sanity: lineage motifs are recovered. PU.1 family (SPI1/SPIB/SPIC) in myeloid and B cells, OCT (POU2F) in B and
  plasma cells, BCL11B / RUNX / T-box / ZEB in T and NK, NFI and TEAD in fibroblasts, NFAT / STAT in endothelium.
- M-6 STAT1+ CXCL10+ vs the other macrophage subtypes of the same donor: only 2 donors have an M-6 unit of >= 10 nuclei.
  Both show higher IRF (IRF7 +5.5 z, IRF4 +3.8, IRF5 +3.7, IRF8, IRF2, IRF9 +3.1–3.3) and NF-kB (REL +4.2, RELA +3.9,
  RELB +3.3, NFKB2 +2.0) motif accessibility. This is descriptive only.
- Donor level (RA, n = 10): the single motifs most correlated with the macrophage RNA signature are noise-prone
  (GLI/SOX, 1,019 motifs). But the declared STAT/IRF set correlates more than other motifs (median rho 0.37 vs
  -0.04; one-sided Mann-Whitney p = 5e-5; motif redundancy makes this p optimistic).

## 10. Knowledge mining: targets and drugs (Open Targets) (**controls ready; leads exploratory**)

`scripts/15_target_drug_mining.py`; `results/tables/target_drug_mining.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`);
cached API responses in `data/raw/opentargets/` (data release 26.09; RA = MONDO_0008383, 8,420 associated targets).

- 115 genes were taken programmatically from steps 05, 10, 12 and 13 (PALM2AKAP2 has no Ensembl 98 ID). The positive controls
  IL6R, JAK1, TNF and MS4A1 are all returned as "approved in RA" (tocilizumab/sarilumab, JAK inhibitors, anti-TNFs,
  rituximab).
- RA stage is taken from each drug's indication list, because a drug's overall maximum stage can come from another disease
  (e.g. maraviroc is phase 4 for HIV but phase 2 in RA). Drug action on the target (inhibitory / activating) is taken
  from mechanism-of-action records.
- Novelty: the macrophage signature genes are largely absent from RA knowledge in Open Targets (median association
  score 0.015, none >= 0.3). The GWAS-derived genes, by contrast, are well known (median 0.42; 21/32 >= 0.3).
- Already tested in RA, not approved: IFNG (fontolizumab, phase 2), CXCL10 (eldelumab, phase 2), CD40 (iscalimab,
  phase 2), CSF2 (otilimab, phase 3), CCR1 and CCR5 (phase 2), ITGB2 (efalizumab, phase 2). So the IFN-gamma / CXCL10 axis
  that the data highlight has already been probed clinically. Emapalumab (anti-IFN-gamma) is approved in
  another indication.
- Repurposing leads (late-stage inhibitor of a gene that is up in disease): SLAMF7 (elotuzumab; in the signature and
  the IFN-gamma leading edge), CD38 (daratumumab, isatuximab), complement C1S (sutimlimab), TLR4 (eritoran), VWF
  (caplacizumab), VEGFC (conbercept).
- GWAS genes with late-stage inhibitors need direction checks. CTLA4 blockers (ipilimumab) act in the wrong direction
  (RA benefits from CTLA4-Ig, abatacept). OX40L (TNFSF4) blockade with amlitelimab is direction-plausible. It fits the
  enrichment of RA loci in OX40+ activated CD4 T cells (section 5).
- These are hypotheses for follow-up, not recommendations. Efficacy and safety in RA are not assessed here.

## 11. Variant-to-motif effects of RA lead SNPs (**exploratory (sequence-based predictions)**)

`scripts/16_variant_motif_effects.py`; `results/tables/variant_alleles.tsv`, `variant_motif_effects.tsv`.
37 RA lead SNPs lie inside state-specific peaks (section 5). Of these, 36 are mappable SNVs whose Ensembl reference allele matches
hg38. The GWAS Catalog risk allele is resolved for 32 (26 alternative, 6 reference). All 1,019 JASPAR motifs are scored
on both alleles (+/-30 bp). An effect is called when only one allele passes the p < 1e-4 match threshold and the
log-odds change by >= 2. The TF must also be detected in >= 10 % of nuclei of the peak's cell state.

- 275 motif-changing allele/TF pairs across 33 SNPs; 93 across 22 SNPs with the TF expressed in the state; 60 across
  18 SNPs after setting aside long zinc-finger motifs (ZNF/ZSCAN/ZBTB), which have low information content.
- Examples with a coherent reading (risk-allele effect, state of the peak):
  rs6859219 (ANKRD55): the risk allele gains STAT5B, BCL6B and RBPJ motifs in CD4+ Tph/Tfh peaks.
  rs2192008 (STAT1/STAT4 region): gains STAT2, VDR and TCF7L2 motifs (proliferating T, plasma, MERTK+ HBEGF+ macrophages).
  rs3087243 (CTLA4): loses an ARNT::HIF1A motif in Treg, OX40+ CD4 and Tph peaks.
  rs9927316 (IRF8 region): loses KLF6/KLF12/SP4 motifs in IFNG+ CD56dim NK and myeloid peaks.
  rs7117261 (CXCR5): loses ELF1/ETV6/IKZF3 motifs in plasma-cell and SPP1+ macrophage peaks.
  rs10872479 (TNFAIP3 region): loses E-box (TCF3/TCF4/TCF12, ZEB1/2) motifs in plasmablast peaks.
- These are predictions from sequence. Allelic imbalance (e.g. allele-specific accessibility in heterozygous donors) or
  reporter assays would be needed to confirm them. Lead SNPs may be in LD with the causal variant.


## 12. Robustness and ablation (**ready**)

`scripts/17_robustness.py`; `results/tables/robustness_summary.tsv`, `robustness_R1_loo_genes.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`).

| Result | Full data | Range under perturbation | Holds? |
|---|---|---|---|
| R1 signature genes still up in active RA (19 leave-one-donor-out DESeq2 refits) | 15/15 | 15/15 in every refit | yes |
| R1 genome-wide Wald statistics vs full data (Spearman) | 1 | 0.93–0.99 | yes |
| R1 share of the 15 genes still FDR < 0.05 | 1.00 | 0.07–1.00 (>= 0.5 in 17/19) | partly |
| R2 active vs healthy, inflammatory macrophage share (Mann-Whitney p) | 0.013 | 0.006–0.048 | yes |
| R2 active vs healthy, homeostatic score (p) | 0.013 | 0.001–0.031 | yes |
| R3 donor-level RNA vs ATAC signature rho (GSE243917) | 0.84 | 0.78–0.91 (all above random-set 95th percentile 0.69) | yes |
| R4 signature transfer vs inflammation grade rho | 0.43 | 0.31–0.59 (always positive, never significant) | direction only |
| R5 RA GWAS enrichment in T-cell peaks, 2,000 / 4,000 / 8,000 peaks | 1.59 (p 0.001) | 1.60 (p 0.002) / 1.24 (p 0.004) | yes |
| R6 GSEA cross-contrast rho, rank by log2FC instead of Wald | -0.61 | -0.55 | yes |

The signature's direction is stable, but its FDR significance rests largely on one donor. Dropping SA139 (naive active
RA, largest sample, 4,199 macrophages, extended sort panel) leaves 1 of 15 genes at FDR < 0.05.

## 13. Blood validation (CELLxGENE PBMC, RA vs matched controls) (**ready (negative)**)

`scripts/18_blood_validation.py`; `results/tables/blood_tests.tsv`; `data/processed/blood_donor_scores.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). 30 donors with >= 100 monocytes (15 RA, 13 complete RA-control pairs);
donor monocyte pseudobulks (classical, non-classical, IL1B+, IFITM3+).

- Only 3 of the 15 signature genes are measurable in monocytes (MT2A, PKM, APOL6). The Hallmark IFN-gamma response
  (105 genes) is the more informative readout.
- Neither score differs between RA and matched controls (paired Wilcoxon p 0.46 and 0.41). Neither tracks DAS28-CRP / CDAI
  within RA (|rho| <= 0.26, p >= 0.39). The IFITM3+ (IFN-like) monocyte share is, if anything, lower in RA (p 0.057).
  Sequencing batch is borderline (Kruskal p 0.07–0.08). Methotrexate use is not associated (p 0.23).
- Interpretation: the IFN-gamma macrophage programme is a tissue state that is not mirrored in circulating monocytes.
  Monitoring it would therefore need synovial tissue rather than blood.

## 14. Cross-indication: the CXCL9 vs SPP1 macrophage axis in NSCLC (TCGA) (**ready**)

`scripts/19_nsclc_macrophage_polarity.py`; `results/tables/nsclc_ra_axis_check.tsv`, `nsclc_cox.tsv`,
`nsclc_tumour_vs_normal.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).

- The axis is first checked in RA. In GSE243917 myeloid cells, CXCL9 is highest in M-6 STAT1+ CXCL10+ (mean 1.29;
  CXCL10 2.30) and SPP1 in M-4 SPP1+ (1.80). These are the IFN-activated and SPP1+ states of the RA analysis.
- Declared H1 (polarity z(CXCL9) - z(SPP1) in primary tumours; Cox OS adjusted for age and stage):

| Cancer | Patients (events) | HR per SD of polarity (95 % CI) | p | Median-split log-rank p |
|---|---|---|---|---|
| LUAD | 491 (179) | 0.85 (0.73–0.995) | 0.043 | 0.12 |
| LUSC | 489 (212) | 0.86 (0.75–0.99) | 0.033 | 0.030 |

- H2: the full RA active-macrophage signature does not predict survival (HR 1.05 LUAD, 1.09 LUSC; ns). What
  transfers is the macrophage polarity axis, not the RA-specific gene set.
- H3: tumours are shifted towards SPP1 relative to adjacent normal lung (median paired difference -1.66 in 58 LUAD
  pairs, p = 5e-10; -2.10 in 51 LUSC pairs, p = 1e-9).
- Interpretation: the same macrophage axis has opposite clinical meaning in the two diseases. In RA the CXCL9/IFN pole marks active joint
  inflammation. In NSCLC it marks anti-tumour immunity and longer survival. Caveats: these are bulk tumour data, and SPP1 is
  also expressed by tumour cells (especially in adenocarcinoma), so the score is not purely macrophage-derived. The
  prognostic association replicates published pan-cancer work; it is not new.

## 15. Spatial transcriptomics (Xenium, GSE347958) (**exploratory spatial association**)

`scripts/20_spatial_xenium.py`; `results/tables/xenium_sections.tsv`, `xenium_ifng_neighbourhood.tsv`,
`xenium_cxcl9_macrophages.tsv`, `xenium_tests.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
17 FFPE sections (Sisó et al. 2026): 11 on a 480-gene panel (no IFNG/STAT1/GBP) and 6 on the 5K panel (4 RA, 1 healthy,
1 OA). Cell types come from per-section Leiden clusters named by marker scores. Macrophage IFN score = mean z of CXCL9,
CXCL10, CXCL11, STAT1, GBP1, GBP5, IDO1, IRF1, CD274, APOL6.

- H1 (declared; IFN score ~ IFNG+ lymphocytes within 30 um + IFNG- lymphocytes within 30 um; IFNG labels
  permuted 1,000x):

| Section | Condition | Macrophages | IFNG+ T/NK | Coefficient per IFNG+ neighbour | per IFNG- neighbour | Permutation p |
|---|---|---|---|---|---|---|
| CT_21045-A-15 | RA | 19,376 | 0.8 % | 0.125 | 0.009 | 0.001 |
| CT_19060-A-4 | RA | 12,955 | 1.5 % | 0.092 | 0.012 | 0.005 |
| CT_513-A-1 | RA | 54,184 | 0.5 % | 0.131 | 0.008 | 0.001 |
| CT_514-A-1 | RA | 3,498 | 1.5 % | -0.015 | 0.006 | 0.61 |
| CT_23031-A-8 | OA | 15,934 | 1.0 % | 0.128 | 0.015 | 0.030 |
| CT_19060-D-8 | healthy | 1,159 | 0 % | n/a | 0.021 | n/a |

  3 of 4 RA sections show a positive coefficient with empirical permutation P <=0.005. Global label shuffling does not keep the spatial clustering of labels or the tissue compartments. The combined z statistic is not used to claim extreme significance.
  Macrophages with an IFNG+ neighbour have a mean IFN score of ~0.16–0.21, vs ~-0.01 for those without one.
- The same neighbourhood relation is present in the OA section, and healthy synovium contains no IFNG+ lymphocytes.
  So the neighbourhood association is not shown to be RA-specific, and it does not identify a causal interaction. For the healthy section, the positive-neighbour coefficient and permutation P cannot be estimated (stored as missing), because there are no IFNG-positive neighbours. A structural zero is not a negative statistical result.
- Secondary: CXCL9/10-high macrophage share, RA median 10 % vs healthy/OA 2.7 % (480 panel; Mann-Whitney p = 0.66).
  One OA block reaches 34 %, and the OA blocks come from one donor.
- This is a complementary spatial association for S2. It does not rescue the non-significant S1 donor-inflammation test in section 8, because the scores differ.

## 16. Spatial proteomics (LCM-MS, PRIDE PXD062696) (**exploratory**)

`scripts/21_proteomics.py`; `results/tables/proteomics_tests.tsv`, `proteomics_detected_genes.tsv`;
`data/processed/proteomics_region_scores.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
22 NH/RA regions (3 NH donors: lining, sublining; 6 RA donors: lining, fibrous sublining, and immune sublining in 4);
2,767 proteins quantified in >= 70 % of regions; per-region median normalisation.

- Coverage: 65 Hallmark IFN-gamma response proteins, but only 3 signature proteins (PKM, PFKP, IL4I1). Metallothioneins
  are small and rarely detected by MS.
- IFN-gamma response score, RA vs NH (donor means, 6 vs 3): lining +1.11 (p = 0.024), sublining +1.31 (p = 0.024, the
  smallest p possible). Abundance-matched random protein sets: 0 % reach the lining difference (specific), but 34 %
  reach the sublining difference (part of a broad RA-vs-NH proteome shift).
- Signature proteins: sublining p = 0.048, not specific (73 % of random sets as large). Inflammatory macrophage proteins
  (S100A8/9, ICAM1) are higher in RA lining and sublining (6 vs 2, p = 0.071) and specific in lining.
- Immune vs fibrous sublining within RA donors: no difference (2 of 4 donors higher).

## 17. Integrated gene prioritisation (**exploratory (ranking)**)

`scripts/22_gene_prioritisation.py`; `results/tables/gene_prioritisation.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
Eleven layers: scRNA disease, M-6 state, treatment response, chromatin, genetics, variant motif, pathway, communication,
protein, eQTL colocalisation and fine-mapping (the last two were added after sections 20 and 24). Each layer scores 1 / 0.5 / 0 by fixed
rules, and the scores are summed (max 11). Druggability and novelty are shown alongside but not summed.

- Top: CD40 (6.0; genetics, colocalisation, fine-mapping, M-6, pathway; tested in RA), STAT4 (5; genetics,
  T/NK-cell colocalisation, M-6, pathway, chromatin) and SOD2 (5; M-6, treatment, pathway, protein; new to RA in Open
  Targets). Next come TNFAIP3, IRF8 and CD2 (4.5; genetics, variant motif, fine-mapping plus functional layers). Then at 4 come SLAMF7
  (new to RA, approved inhibitor elotuzumab), ANKRD55 and FCGR2A (colocalisation), and CXCL9/CXCL10, STAT1, APOL6, LAP3 (the
  IFN-gamma axis). L10 also counts immune-cell colocalisation from section 27.
- SLAMF7 is the clearest repurposing candidate (support from several layers, little prior RA evidence, approved antibody).
- Caveats: equal layer weights are a choice. Nearest-gene assignment can promote bystanders (e.g. RPLP1). The score
  ranks hypotheses; it is not a probability of causality.

## 18. Knowledge graph (**ready (deliverable)**)

`scripts/23_knowledge_graph.py`; `results/knowledge_graph/nodes.tsv`, `edges.tsv`, `ra_multiomics_kg.graphml`,
`kg_viewer.html`. 460 nodes (175 genes, 68 TFs, 58 drugs, 67 cell states, 53 variants, 36 pathways, 3 disease /
treatment nodes) and 731 typed edges. These include colocalisation and risk-allele direction, Visium, fine-mapping, S-LDSC,
disease-stage, trajectory, TF-activity, Tph spatial, scDRS and MR edges. Every edge records the step that produced it, its statistic and its value. The HTML viewer
(vis-network) supports search and highlights a node's neighbours. Isolated nodes are hidden in the viewer only.

## 19. Visium (GSE347797) (**ready (donor level); exploratory (within tissue)**)

`scripts/24_visium.py`; `results/tables/visium_sections.tsv`, `visium_within_section.tsv`, `visium_tests.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). 12 FFPE sections (6 RA, 3 healthy, 3 OA blocks from one donor); spots with >= 200 UMIs.

- H1 (donor level): share of spots above the healthy 95th percentile of the IFN programme score. RA median 32 %
  (range 10.5–75.6 %) vs healthy 5.5 % (4.4–7.3 %). Every RA section exceeds every healthy section (Mann-Whitney
  p = 0.024, the minimum for 6 vs 3). The single OA donor's three blocks range 5–56 %.
- H2 (within RA sections; partial Spearman adjusted for spot depth; vs 200 expression-matched random gene sets): the
  IFN score co-localises with the lymphoid-aggregate score beyond random in 2/6 sections (CT_510, CT_19060-A-4). It co-localises
  with IFNG in the spot and its six neighbours in 2/6 (CT_21045-A-15, CT_19060-A-4).
- Read together with section 15: Visium confirms that RA tissue carries far more IFN-programme signal than healthy tissue.
  The 55 um multicellular spots and sparse IFNG detection make within-tissue coupling weaker than in single-cell Xenium.

## 20. eQTL colocalisation (RA GWAS x myeloid eQTLs) (**ready**)

`scripts/25_eqtl_coloc.py`; `results/tables/eqtl_coloc.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
RA GWAS: Ishigaki 2022 European ancestry (GCST90132223). eQTL Catalogue gene-expression datasets: Alasoo 2018 macrophage
naive and IFN-gamma 18 h (n = 84), Nedelec 2016 macrophage naive (n = 163), Quach 2016 monocyte naive (n = 200). Both were read
by remote tabix (the eQTL Catalogue REST API is no longer available). coloc.abf with default priors; TSS +/- 500 kb; MHC excluded.
37 genes with genetic evidence; 136 gene x dataset tests with >= 50 shared variants.

| Gene | Dataset | Shared variants | GWAS min p | eQTL min p | PP.H4 |
|---|---|---|---|---|---|
| IRF5 | Nedelec 2016 macrophage naive | 3,265 | 4e-17 | 3e-32 | **0.96** |
| IRF5 | Quach 2016 monocyte naive | 3,268 | 4e-17 | 2e-14 | **0.94** |
| CD40 | Nedelec 2016 macrophage naive | 3,166 | 6e-13 | 5e-5 | **0.87** |
| TDRD10 (IL6R locus) | Alasoo 2018 macrophage IFN-gamma | 1,844 | 5e-7 | 1e-5 | 0.70 (suggestive) |
| TSPAN32 | Nedelec 2016 macrophage naive | 3,530 | 1e-3 | 1e-11 | 0.57 (suggestive) |

- Distinct causal variants (PP.H3 > 0.8): STAT4 in all four myeloid datasets. These tests do not support a shared expression/RA signal in the tested myeloid datasets. Absence of support does not exclude myeloid mediation. CD40 is H3 in monocytes but H4 in macrophages, which suggests cell-state
  specificity. Also ETV7, PFKFB3, PLGRKT.
- Cross-layer link: IRF5 is a master regulator of inflammatory macrophage polarisation. Its motif is among the most
  accessible in M-6 STAT1+ CXCL10+ macrophages (section 9, +3.7 z). The shared IRF5 GWAS/expression-QTL signals nominate myeloid expression as a candidate mediator. They do not establish that mediation, or any regulation of M-6. (The prioritisation score
  does not count TF-motif evidence, so IRF5 ranks low there. This is stated here rather than fixed by re-weighting.)
- Caveats: single-causal-variant assumption; European GWAS with mixed-sample eQTLs. The IFN-gamma macrophage dataset is
  small (84 donors), so absence of colocalisation there is weak evidence.


## 21. Signature-reversal compound search (L1000CDS2) (**exploratory**)

`scripts/26_l1000_reversal.py`; `results/tables/l1000_reversal_signatures.tsv`, `l1000_reversal_compounds.tsv`,
`l1000_reversal_class_enrichment.tsv`. Top-50 reversing L1000 signatures for Q1 (active-RA macrophage up 50 / remission
up 50) and Q2 (IFN-gamma leading edge up / homeostatic programme down). Classes come from ChEMBL mechanisms plus a documented
compound-name dictionary (many tool compounds have no ChEMBL mechanism record).

- Q1: HSP90 inhibitors (geldanamycin #1, radicicol), mTOR inhibitors (temsirolimus #2, torin-2), BET inhibitors
  (PFI-1 #5, I-BET151 #6, JQ1) and JAK inhibitors (ruxolitinib #15, TG-101348). Together these are 9 of 35 named compounds.
- Q2: NF-kB pathway inhibitors (withaferin A #1, BAY 11-7821 #5), PI3K/mTOR and BET inhibitors.
- Baseline: 20 random gene-list queries give a median of 14 % (3–31 %) named hits in these classes. Q1 gives 26 %
  (empirical p = 0.095) and Q2 14 % (p = 0.52). The classes are mechanistically plausible (they suppress JAK-STAT / NF-kB
  inflammatory programmes). But they are common in the L1000 library, so the enrichment is not significant. L1000 profiles
  come from cancer and other cell lines, not synovial macrophages.

## 22. Bulk deconvolution of the anti-TNF cohort (**negative (method), not interpreted**)

`scripts/27_bulk_deconvolution.py`; `data/processed/gse198520_deconvolution.tsv`, `results/tables/deconvolution_tests.tsv`.
NNLS with a GSE243917 single-nucleus reference (13 states including M-6, SPP1+, MERTK+ resident, monocyte-derived).
- Sanity checks failed in both versions (with and without per-gene scaling). The estimated lymphoid share did not differ
  between lymphoid and fibroid pathotypes (p 0.10 / 0.26) or track the histological inflammatory score (rho 0.20 / 0.12).
  The estimated compositions were implausible (e.g. fibroblasts 3–7 %, one state absorbing up to ~50 %).
- Apparent M-6 associations with response track library size (rho 0.63), and the adjusted model separates. They are
  not interpreted. Likely causes: a nuclear reference vs whole-cell bulk RNA, a low-depth targeted panel, and an M-6
  profile built from 65 nuclei. No further tuning was tried, to avoid searching for a desired result.

## 23. Integration baselines: WNN vs single-modality, concatenation and MOFA+ (**ready (benchmark)**)

`scripts/28_integration_benchmark.py`; `results/tables/integration_benchmark.tsv`, `mofa_factors_top_features.tsv`,
`mofa_variance_explained.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`). Every method uses the same 31,520 nuclei and input features.
All are Harmony-corrected by sample except the uncorrected RNA reference. MOFA+ used 20 factors on 3,000 HVG
and 20,000 peaks (settings fixed beforehand; ~30 min on CPU). Metrics are on 20-NN graphs; Leiden over 8 resolutions.

| Method | Cell types ARI | Subtypes ARI | ATAC clusters ARI | M-6 neighbour purity | Sample mixing | Overall (0.6 bio + 0.4 batch) |
|---|---|---|---|---|---|---|
| Concatenation + Harmony | **0.95** | **0.42** | 0.75 | 0.19 | 0.81 | **0.84** |
| WNN (step 07) | 0.83 | 0.41 | 0.73 | 0.09 | 0.82 | 0.71 |
| ATAC LSI + Harmony | 0.88 | 0.37 | 0.83* | 0.10 | 0.87 | 0.66 |
| RNA PCA + Harmony | 0.81 | 0.41 | 0.65 | **0.24** | 0.83 | 0.54 |
| RNA PCA, uncorrected | 0.79 | 0.29 | 0.60 | 0.19 | 0.42 | 0.29 |
| MOFA+ + Harmony | 0.83 | 0.24 | 0.61 | 0.07 | 0.86 | 0.26 |

\* ATAC clusters were defined on this graph (step 06), so ATAC alone is favoured there, as RNA is for author labels.

- **A simple baseline is competitive.** Concatenating the two Harmony-corrected embeddings scores highest overall. It is
  best or tied on every biology metric. WNN comes second overall. Both joint methods beat either modality alone. Harmony
  matters most (uncorrected RNA mixes samples half as well).
- **Run-to-run variation.** The ARI values depend on the ATAC LSI start and on the Leiden clusterings. An earlier, unseeded
  LSI run put WNN slightly ahead of concatenation on subtypes (0.43 vs 0.41) and ATAC clusters (0.78 vs 0.72), and
  ATAC alone at 0.74 on cell types. Differences of this size (≈ 0.05 ARI) between the top methods are therefore not
  meaningful. The ranking concatenation ≥ WNN > single modalities > MOFA+ held in both runs. The LSI is now seeded and
  reproduces exactly (section 36, item 16).
- **Modality balance (label-free):** WNN neighbourhoods overlap the ATAC graph more than the RNA graph (0.27 vs 0.16).
  Concatenation is even (0.20 vs 0.21).
- **Rare state:** no method isolates M-6 STAT1+ CXCL10+ as a cluster (best F1 0.04). Its neighbourhoods are purest on
  RNA alone (0.24) and are diluted by joint methods (WNN 0.09). This fits M-6 being defined by an IFN RNA programme
  that chromatin separates only weakly (section 4). For rare activated states, RNA-only analyses are still needed.
- **MOFA+** ranks last in this single configuration. Its factors follow cell identity (Factor 1: 15 % of RNA, 1.6 % of
  ATAC variance). The factor that best tracks the IFN programme in myeloid nuclei (Factor 20, r 0.27) is highest in M-6
  (median 0.73 vs 0.01 in other myeloid nuclei). But it explains 0.25 % of RNA variance, and its top weights are T-cell
  genes. MOFA+ is designed for interpretable sample-level factors rather than single-nucleus neighbourhoods. We did not
  tune it, to keep the comparison fair.
- Caveats: one dataset; min-max scaling across six methods exaggerates differences; the overall score depends on the
  scIB weights; MultiVI was not run (CPU only).

## 24. Fine-mapped RA variants in cell-state chromatin (**T-cell signal replicates in FinnGen; carried by few loci (S-LDSC, section 28, shows it is genome-wide)**)

`scripts/29_finemap_enrichment.py`; `results/tables/finemap_variants.tsv`, `finemap_enrichment.tsv`, `finemap_tests.tsv`,
`finemap_candidates.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
95 % credible sets from Open Targets 26.09, MHC excluded: Ishigaki 2022 European ancestry (SuSiE-inf; 64 sets, 1,048
variants, 20 with PIP >= 0.5), FinnGen R12 RA (SuSiE; 41 sets, 1,057 variants, 15 with PIP >= 0.5), and LDL cholesterol as control (Sakaue 2021; 165 sets, 4,728 variants). It has not been established here that FinnGen participants are independent of the Ishigaki meta-analysis. The state-specific peaks and
matched null are the same as in section 5. Statistic = PIP mass of credible-set variants within 250 bp of a specific-peak midpoint.

| State | Ishigaki 2022 | FinnGen R12 |
|---|---|---|
| T cells | 6.2x (p 0.020, FDR 0.08) | **11.8x (p 0.001, FDR 0.008)** |
| Plasma cells | 6.4x (p 0.012, FDR 0.08) | 6.8x (p 0.016, FDR 0.06) |
| T-10 CD4+ OX40+ NR3C1+ | 12.9x (p 0.001, FDR 0.06) | 8.8x (p 0.006, FDR 0.17) |
| T-7 CD4+ Tph | 6.7x (p 0.019) | 8.5x (p 0.005) |
| T-3 CD4+ Tfh/Tph | 6.4x (p 0.020) | 6.8x (p 0.015) |
| Myeloid | 1.9x (p 0.14) | 0.7x (p 0.22) |

- F1 (T cells, Ishigaki): nominal (p 0.02, FDR 0.08). F2 (replication): FinnGen T cells FDR 0.008. The enrichments of all
  65 states agree between the two studies (Spearman rho 0.43, p 3e-4). F3 (control): LDL has no state at FDR < 0.05.
- **Carried by few loci.** 82 % of the Ishigaki T-cell PIP mass comes from one credible set (CD2, rs798000). The
  Tph/Tfh/OX40+ T signal comes mostly from ANKRD55 rs7731626 (PIP 1.00 in both studies). The LDL control reaches a nominal 6x
  in endothelium from a single variant, which shows what one high-PIP variant can do. The replication is therefore of
  the same loci across resources, not a broad polygenic pattern. Participant independence is not established here. The result agrees with the window-based result of
  section 5 (T cells 1.59x) and with published T-cell heritability enrichment.
- **Candidate causal variants** (PIP >= 0.1 inside a state-specific peak): 19 variants, 6 of them step-10 lead SNPs.
  ANKRD55 rs7731626 (PIP 1.00, both studies; Tph/Tfh/OX40+ T, plasma), CD2 rs798000 (1.00, both; T cells, lining
  fibroblasts), ARID5B rs77509998 (1.00; CD11c+ ABC B cells, OX40+ T), CLIC6/RUNX1 rs143365361 (0.84), IRF8 rs9927316
  (0.27; IFNG+ NK and myeloid), TNFAIP3, RUNX1, JAK1, CTLA4 rs3087243 (0.13) and CD40 rs1883832 (0.10).
- Consequences for earlier sections: the ANKRD55 motif gain of section 11 (rs6859219) sits on a variant outside both
  credible sets, so rs7731626 is the better candidate there. The IRF5 locus fine-maps to two independent single-variant
  sets (rs4731532 and rs3778754, PIP 1.00 each), and neither lies inside a state-specific peak. Two causal signals break coloc's
  single-variant assumption. The multi-signal colocalisation of section 27 confirms that both colocalise with myeloid
  IRF5 expression.
- Myeloid states are not enriched in either study, as in section 5.
- **East Asian ancestry (F4, added in a later round; Ishigaki 2022 EAS, 51 sets).** Broad T cells are not enriched (1.26x,
  p 0.22; one contributing set). But CD4+ Tfh/Tph peaks are (T-3 14.0x, p 0.001, FDR 0.06; T-7 Tph 8.4x, p 0.011).
  This is carried by the same ANKRD55 variant rs7731626 (PIP ~1.0 in East Asians as well) and by an ARID5B signal with a different
  lead variant. State-level enrichments correlate with the European analysis (rho 0.42, p 6e-4). The CD2 signal that
  dominates the European T-cell result is absent. So the cross-ancestry evidence is for the ANKRD55 causal variant in
  Tph/Tfh chromatin, not for a general T-cell enrichment.

## 25. Disease-stage validation in bulk synovium (GSE89408) (**IFN programme ready; primary signature negative**)

`scripts/30_gse89408_disease_stage.py`; `data/processed/gse89408_sample_scores.tsv`, `results/tables/gse89408_tests.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). 218 biopsies (28 healthy, 10 arthralgia, 6 undifferentiated arthritis,
57 early RA, 95 established RA, 22 OA). OLS with age, sex and log library size (HC3 errors). Every test was repeated on
1,000 random expression-matched gene sets of the same size.

| Score | RA vs healthy | RA vs OA | Trend healthy < arthralgia < early RA | RA vs OA, + macrophage abundance |
|---|---|---|---|---|
| S1 active-RA signature (primary, 13/15 genes) | p 0.12; vs random p 0.52 | p 0.33; vs random 0.52 | p 0.26; vs random 0.57 | p 0.22; vs random 0.48 |
| S2 IFN macrophage programme | **t 8.7, p 3e-18; vs random 0.011** | **t 7.7, p 1e-14; vs random 0.002** | t 6.0, p 2e-9; vs random 0.055 | **t 9.1; vs random 0.001** |
| S3 homeostatic MERTK+ | p 0.03; vs random 0.64 | p 1.0 | p 0.01; vs random 0.57 | p 0.95 |

- **Random gene sets also differ between groups** (median t 1.1–3.5), even after median-of-ratios normalisation. Group
  differences in RNA quality or composition shift all genes. Only differences beyond that background are interpreted.
- The declared primary signature (15 genes from sorted macrophages, active vs remission RA) does not separate stages
  in bulk tissue beyond random genes. The a priori IFN macrophage programme does. It is higher in RA than in healthy tissue and than in
  OA (similar ages, 53–55 years), independent of macrophage abundance. It is also higher against healthy donors aged >= 40 only (vs
  random p 0.05). Its rise from healthy through arthralgia to early RA is borderline beyond random (p 0.055; 0.13 with
  macrophage adjustment). Early vs established RA and ACPA+ vs ACPA- do not differ.
- The homeostatic programme's apparent rise in RA vs healthy disappears once macrophage abundance is adjusted for (healthy
  synovium holds far fewer macrophages). This adjustment was added after seeing that result and is reported as such.
- Interpretation: in bulk tissue, the IFN-activated macrophage axis (CXCL9/10/11, STAT1, GBP, IDO1) is the RA-specific signal that holds up.
  This agrees with Visium (section 19) and proteomics (section 16). The metallothionein/glycolysis
  part of the sorted-cell signature does not carry over to bulk.

## 26. Myeloid trajectory and chromatin priming (**negative / not interpretable**)

`scripts/31_myeloid_trajectory.py`; `results/tables/myeloid_paga.tsv`, `myeloid_trajectory_tests.tsv`,
`myeloid_priming_by_state.tsv`, `myeloid_lag_sensitivity.tsv`, `myeloid_lag_controls.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). 5,849 myeloid nuclei without DC-like states; RNA graph (2,000 HVG, Harmony)
and chromatin graph (ATAC LSI, Harmony); M-6 programme = 50 RNA markers (CXCL10, IDO1, IFIT1, GBP1, CD274, ...).

- **T1 (where M-6 sits): not supported on RNA.** On the RNA graph, M-6's strongest PAGA neighbour is the resident
  M-1 MERTK+ SELENOP+ LYVE1- state (0.83 vs 0.73 for monocyte-derived M-7), in all 12 leave-one-sample-out runs. On the
  chromatin graph, M-7 and M-8 PLCG2+ tie at the top (1.00), with M-1 close behind (0.98). The data fit an IFN-activated state that can be reached from both
  resident and monocyte-derived macrophages. They do not show monocyte origin.
- **T2 (chromatin priming in M-7): negative.** M-7 nuclei carry more M-6 programme chromatin relative to RNA than
  resident nuclei in 11/11 samples (Wilcoxon p 0.001). But random matched gene sets show a larger gap (median 0.72 vs
  0.39 observed; p vs random 0.84). This is a general difference between the states' RNA and chromatin layers, not an M-6 effect.
- **T3 (does chromatin rise before RNA towards M-6?): not interpretable.** Pseudotime from the M-7 root fails its sanity
  check. M-7 nuclei come out latest, not first (median 0.21 vs 0.10–0.19 for other states). This holds whether the root is taken
  as the most monocyte-like nucleus or from smoothed scores. So the nominal lag (chromatin half-max 0.16 quantiles earlier;
  p vs random 0.03; 7/9 bin x threshold settings) has no M-7 -> M-6 direction. The a priori IFN gene list
  shows no lead beyond random (p 0.17). Other end states with the same marker selection show no lag. So the selection
  artefact is not what drives it; the ordering is. Not interpreted.
- From the homeostatic M-0 root (which passes its sanity check), M-6 is the second-latest state after M-7. This fits
  an activated end state rather than a precursor.

## 27. Multi-signal colocalisation and direction of effect (**ready (IRF5 resolved; direction table exploratory)**)

`scripts/32_multisignal_coloc.py`; `results/tables/coloc_multisignal_all.tsv`, `coloc_multisignal_genes.tsv`,
`coloc_direction.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`). Each fine-mapped RA signal (105 credible sets: Ishigaki
EUR 64, FinnGen 41) was matched to fine-mapped QTL credible sets by Open Targets 26.09 (coloc on SuSiE sets plus eCAVIAR;
eQTL Catalogue incl. single-cell eQTL, GTEx v10, UKB-PPP pQTL). 8,016 records; 5,666 with H4 >= 0.8, of which 4,773 are cis.
893 trans or unplaced records were dropped, because UKB-PPP trans-pQTLs had attached hundreds of distant proteins to immune loci.

- **C1 IRF5: both signals colocalise with myeloid IRF5.** Ishigaki rs3778754: 11 myeloid IRF5 QTLs with H4 0.93–1.00
  (macrophages naive and IFN-gamma 18 h, Alasoo; classical monocytes naive / LPS / influenza, Quach; macrophages after
  Listeria / Salmonella, Nedelec). rs4731532: CD14+ monocytes (BLUEPRINT, H4 0.99). FinnGen rs13246321: 9 myeloid QTLs.
  This resolves the single-variant caveat of section 20. The direction in myeloid cells is not consistent (the risk
  allele raises IRF5 in 62 % of myeloid datasets), while in T/NK cells it raises IRF5 in all 11. Whether to inhibit IRF5
  cannot be read from eQTL direction.
- **C2:** 70/105 RA signals colocalise with at least one cis-QTL, and 58/105 with an immune-cell QTL (66 genes for
  Ishigaki, 50 for FinnGen). 20 genes colocalise in both studies (ALDH2, ANKRD55, CCR6, COG6, CTLA4, FDX2, ICAM4/5,
  IL6ST, IRF5, KRI1, PDE4A, POU3F1, PVT1, RAVER1, RBPJ, TNPO3, TYK2, ZGLP1).
- **C3 direction** (one call per QTL dataset; expression and cis-protein QTLs only): of 91 genes with immune-cell
  colocalisation, the risk allele raises the level of 26, lowers 30, is mixed for 6 and rests on one dataset for 29.
  - **CTLA4: risk lowers expression in T/NK cells (21/21 datasets).** This is biologically compatible with CTLA4-mediated inhibitory signalling, but it is not pharmacological validation. Abatacept binds CD80/CD86, and it cannot be equated with changing CTLA4 expression.
  - **ANKRD55 locus (rs7731626): risk raises ANKRD55 (62 datasets) and IL6ST (gp130, 18 T/NK and B datasets).** IL6ST is the
    signal-transducing chain of the IL-6 receptor family. The shared IL6ST/ANKRD55 instrument cannot identify the mediator. Tocilizumab/sarilumab target IL6R, not IL6ST, and no drug efficacy follows from this association.
  - CCR6 and NAB1: risk raises expression (inhibition fits). FCGR2A, IFNGR2, PHF19: risk lowers expression. CD40 is mixed
    (4 B/LCL datasets). STAT4 rests on one T/NK dataset (risk lowers). TRAF1 is cell-type specific (lowers in myeloid,
    raises in T/NK).
- Caveats: QTL tissues are not synovium. Direction from expression level does not fix drug direction for receptors,
  ligands or loss-of-function coding variants (e.g. TYK2 has no immune level-QTL colocalisation here). Results depend on
  Open Targets release 26.09.

## 28. S-LDSC heritability enrichment (**ready (T-cell signal genome-wide, replicated, control clean)**)

`scripts/33_sldsc.py`; `results/tables/sldsc_results.tsv`, `sldsc_allpeaks_summary.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
Annotations on GRCh38 (peaks +/- 500 bp): all 131,644 synovial peaks (6.5 % of reference SNPs), top 10 % specific peaks
per broad type (~0.45 %), top 4,000 per subtype (~0.15 %). Model per state: baselineLD v2.2 + all peaks + state. 1000G EUR
LD (1 cM), HapMap3 regression SNPs, MHC excluded; polyfun's Python 3 LDSC. Total observed-scale h2: 0.20 (Ishigaki,
intercept 1.05), 0.14 (FinnGen RA, 1.09), 0.17 (T2D, 1.39).

| Annotation | Ishigaki 2022 EUR: tau z (FDR), enrichment | FinnGen R12 RA | T2D control |
|---|---|---|---|
| All synovial peaks | 3.1, 9.3x (p 5e-8) | 4.6, 13.5x (p 3e-10) | -0.3, 2.3x (p 0.14) |
| **T cells** | **4.7 (1e-5), 46x** | **4.7 (1e-5), 53x** | -1.3 |
| NK | 3.0 (0.005), 28x | 2.5 (0.018), 32x | -1.8 |
| Plasma | 2.8 (0.008), 29x | 3.6 (6e-4), 41x | 0.8 |
| B cells | 2.6 (0.009), 29x | 2.3 (0.022), 36x | 0.5 |
| Myeloid | -2.2 (0.99) | -1.8 (1.0) | 0.9 |
| Fibroblast / endothelial / mural | all negative | all negative | n.s. |

- **S1 supported and replicated:** T-cell-specific chromatin carries RA heritability beyond the baseline model and beyond
  synovial accessible chromatin in general, in both GWAS (tau z 4.7). This answers the caveat of section 24: the T-cell
  signal is genome-wide, not only the CD2 and ANKRD55 loci. The lymphocyte lineages (NK, plasma, B) follow.
- **S2 subtypes:** CD4+ Tph (T-7), Treg (T-8, T-9) and OX40+ NR3C1+ (T-10) reach FDR 0.02 in FinnGen. Tfh/Tph (T-3), Tph,
  Tregs and HLA-DR+ IgG+ plasmablasts reach FDR 0.055 in Ishigaki. These are the same CD4 states as in sections 5 and 24. Resident and
  inflammatory macrophage subtypes are not enriched (M-2, M-3, M-4 tau z -1.5 to -2.4). M-6 cannot be tested (65 nuclei,
  below the 100-nucleus minimum for a peak set).
- **S3 control:** T2D shows no enrichment in any of 65 annotations (FDR >= 0.80) and none for synovial peaks overall.
- Interpretation: conditional RA heritability is enriched in lymphocyte-associated regulatory DNA. Non-positive myeloid coefficients cannot show that macrophages are downstream or unimportant, and M-6 itself was not tested. IRF5/CD40 provide locus-level candidate myeloid associations (sections 20, 27).
- Caveats: European LD panel; 200-bp peaks extended by 500 bp; subtype annotations are small (low power); the
  synovial snATAC comes from 12 samples.

## 29. TF activity from RNA vs motif accessibility (**sanity available; M-6 candidate regulators exploratory (2 donors)**)

`scripts/34_tf_activity_crossmodal.py`; `results/tables/tf_activity_units.tsv.gz`, `tf_lineage_sanity.tsv`,
`tf_crossmodal_concordance.tsv`, `tf_M6_drivers.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`). The same 514 pseudobulk
units as chromVAR (section 9). Univariate linear model (as in decoupler) on CollecTRI regulons (464 TFs with >= 10 targets).
336 TFs also have a single-TF JASPAR motif.

- **T1 sanity (5 of 6 lineages):** SPI1 #1 in myeloid, EBF1 #4 (PAX5 #11) in B cells, XBP1 #1 in plasma cells, TWIST1 #6
  in fibroblasts, ERG #3 in endothelium. T cells miss (TCF7 #419, LEF1 #446). Synovial T cells are memory/effector
  cells whose top TFs are STAT4, BATF and MAF.
- **T2 cross-layer agreement is weak overall:** within 70 macrophage units, a TF's RNA-inferred activity tracks its own
  motif slightly more than other families' motifs (median rho 0.022 vs -0.002, p 0.024). Agreement is concentrated in AP-1
  (JUN, JUNB, JUND, FOS; rho 0.52–0.61) and NF-kB (RELA, NFKB1, NFKB2, REL; 0.54–0.56). These factors have strong,
  well-defined motifs and large regulons.
- **T3 M-6 candidate regulators:** the submission revision applies Benjamini-Hochberg adjustment to all tested RNA gene-label-permutation P values. It then combines q <0.05 with a positive top-decile motif ranking. The nominal screen is kept separately in both_layers_nominal. IRF5's nominal P~0.036 does not survive this adjustment (q~0.139). Its colocalisation and motif evidence remain separate exploratory observations. Final corrected candidates are in tf_M6_drivers.tsv. This adjustment does not provide donor-level replication. It also does not fix the correlation assumptions of the gene-label null.
- Caveats: M-6 units exist in 2 donors. Gene-label permutation ignores gene-gene correlation, so permutation p-values are
  optimistic. NF-kB dimer regulons overlap. RNA and ATAC come from the same nuclei.

## 30. The CD4 Tph / Treg axis across layers (**spatial and communication ready; bulk and composition negative**)

`scripts/35_tph_treg_axis.py`; `results/tables/tph_treg_tests.tsv`, `tph_treg_composition.tsv`, `tph_spatial.tsv`;
`results/figures/main/` (current figures; panel map in `Figure_notes.md`). This follows section 28: the states that carry RA heritability (CD4 Tph T-7, Tfh/Tph
T-3; Treg T-8, T-9).

- **P1 composition (10 RA donors with >= 50 T cells):** Tph share of T cells vs lymphocytic inflammation grade, rho 0.24
  (p 0.50); Treg, -0.12 (p 0.74). Tph share vs M-6 share of macrophages, rho 0.45 (p 0.19). No prep-batch difference.
  This is underpowered rather than informative.
- **P2 bulk (GSE89408):** adjusted for age, sex, depth and T-cell abundance, neither the Tph programme (t -1.6; vs random
  p 0.81) nor the Treg programme (t 1.0; p 0.38) is higher in RA than OA. Within RA, the IFN macrophage programme does
  not track the Tph programme beyond random sets (p 0.17). In bulk tissue these programmes follow T-cell infiltration.
- **P3 communication:** CXCL13 -> CXCR5 from both Tph states to B cells (unswitched memory B-1) is in the top 0.05 % of
  each Tph state's interactions by specificity. IL21 is too sparse in snRNA to be scored. Tph ranks #10 among IFNG senders
  to M-6 (top: NK-1 IFNG+ CD56dim NK).
- **P4 spatial (4 RA sections, Xenium 5K):** B cells within 30 um are far more frequent around CXCL13+ T cells (5–9 % of
  T cells) than around CXCL13- T cells, in 4/4 sections. The section-level empirical permutation p is about 0.001 each; the former combined normal-tail P is removed from the final summary. Higher macrophage IFN scores next to CXCL13+ T cells are not supported (0/4 sections with empirical P<0.05).
- **P5 anti-TNF (exploratory):** the Tph programme falls after anti-TNF in good responders (median change -0.26, Wilcoxon
  p 0.045). Good minus non-responders: -0.33 (Mann-Whitney p 0.15; beyond random gene sets p 0.035).
- Interpretation: CXCL13-positive T cells are associated with B-cell-rich neighbourhoods, consistent with a candidate interaction. This does not establish aggregate organisation or a causal Tph-to-macrophage chain. Bulk Tph scores cannot be separated from T-cell abundance.
- Caveats: CXCL13 alone marks Tph in Xenium (Tfh-like cells share it); 4 sections; donor-level tests have 10 donors.

## 31. scDRS single-cell polygenic risk (**complementary expression-based genetic association; M-6 modest**)

`scripts/36_scdrs.py`; `results/tables/magma_genes_<trait>.tsv`, `scdrs_groups.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`). MAGMA
(SNP-wise mean, 10-kb windows, 1000G EUR GRCh38, MHC excluded) was run on the full summary statistics. scDRS was run on all 31,520 nuclei,
with the top 1,000 MAGMA genes outside the extended MHC and 1,000 matched control gene sets.

| Group | Ishigaki 2022 EUR (MC z) | FinnGen R12 RA | T2D control |
|---|---|---|---|
| T cells | **10.5** (928 nuclei FDR < 0.1) | **8.5** (430) | -1.9 |
| NK | 5.6 | 2.5 | -0.8 |
| B cells | 4.1 | 1.8 | -0.7 |
| Myeloid | 2.9 | 2.7 | 1.2 |
| Plasma | 0.5 | 0.7 | 0.2 |
| Fibroblast / endothelial / mural | -3.6 / -2.7 / -4.0 | -3.6 / -3.4 / -2.8 | 1.6 / **4.4 / 6.7** |
| T-8 CD25-high Treg | **11.5** | **9.6** | -1.9 |
| T-7 CD4+ Tph | **11.1** | **8.5** | -0.4 |
| T-10 CD4+ OX40+ NR3C1+ | 9.4 | 7.4 | -1.2 |
| T-3 CD4+ Tfh/Tph | 8.2 | 6.8 | -0.8 |
| M-6 STAT1+ CXCL10+ | 3.0 (MC p 0.006) | 1.8 (p 0.041) | -1.2 (p 0.90) |

- **D1 supported:** T cells, and above all CD25-high Treg, Tph, OX40+ and Tfh/Tph states, express RA risk genes far
  more than matched controls in both GWAS. These are the same states as in S-LDSC (section 28). Here they appear at single-nucleus resolution and
  with a method based on gene expression rather than chromatin. T2D maps to vascular cells instead (mural z 6.7,
  endothelial 4.4) and not to T cells, so the result is trait-specific.
- **D2 M-6:** nominally associated in both RA GWAS (p 0.006, 0.041) and not with T2D, but not distinctive within the
  myeloid compartment. Dendritic cells (DC1 z 4.7, LAMP3+ 4.2) and several macrophage states (MERTK+ HBEGF+ 3.8, SPP1+ 3.3)
  score as high or higher. Myeloid cells overall show a weak association (z 2.7–2.9) that S-LDSC did not detect.
- **D3:** RA scores are heterogeneous within T cells (heterogeneity MC p 0.001), consistent with subtype-level differences.
- Caveats: snRNA detection is sparse. scDRS tests expression of risk genes, not that risk variants act in these cells. The
  MAGMA gene windows (10 kb) assign some variants to neighbouring genes.

## 32. Drug-target Mendelian randomisation (**exploratory expression-associated estimates; no drug-efficacy validation**)

`scripts/37_drug_target_mr.py`; `results/tables/mr_wald_all.tsv`, `mr_primary.tsv`; `results/figures/main/` (current figures; panel map in `Figure_notes.md`).
906 colocalising level-QTL pairs (95 genes) from section 27 with a QTL lead variant. GWAS effects were found for 366/453
(Ishigaki) and 332/453 (FinnGen) instruments. Wald ratio per signal x gene; primary dataset = highest CLPP in immune cells.

| Gene | GWAS | Primary dataset | OR per SD (95 % CI) | Datasets with same sign |
|---|---|---|---|---|
| **CTLA4** | Ishigaki | OneK1K CD4+ T (sc-eQTL) | **0.71 (0.66–0.77)** | 6/6 |
| **CTLA4** | FinnGen | CEDAR CD4+ T | **0.86 (0.83–0.89)** | 14/14 |
| **IL6ST** | Ishigaki | IBDverse regulatory T (sc-eQTL) | **1.48 (1.39–1.58)** | 9/9 |
| **IL6ST** | FinnGen | IBDverse regulatory T | **1.27 (1.20–1.34)** | 9/9 |
| ANKRD55 | Ishigaki / FinnGen | IBDverse central memory CD4+ T | 1.30 / 1.17 | 35/35 |
| CD40 | Ishigaki | IBDverse plasma cell | 1.26 (1.18–1.34) | 3/4 |
| FCGR2A | Ishigaki | Alasoo macrophage IFN-gamma | 0.93 (0.91–0.95) | 7/7 |
| CCR6 | Ishigaki / FinnGen | Fairfax B cell | 1.35 / 1.25 | 1/2 |
| IRF5 | Ishigaki rs3778754 | Alasoo macrophage IFN-gamma | 0.88 (0.85–0.90) | 4/18 |
| IRF5 | Ishigaki rs4731532 | GTEx LCL | 1.32 (1.21–1.42) | 3/5 |

- **M1 expression association:** higher CTLA4 expression is associated with lower RA risk in both GWAS. Abatacept binds CD80/CD86. This is biological context, not evidence that genetic expression reproduces drug action.
- **M2 expression association:** higher IL6ST/gp130 expression is associated with higher risk. IL6ST is not IL6R, the tocilizumab/sarilumab target. The shared rs7731626 cannot separate IL6ST from ANKRD55.
- Exploratory: higher CD40 in plasma/B cells raises risk (this fits CD40/CD40L blockade, which is in RA trials). Higher FCGR2A in
  IFN-gamma macrophages protects. IRF5 estimates conflict between datasets and signals, so no drug direction is supported.
- Caveats: one instrument per signal (pleiotropy not testable). QTL effects were measured in blood and cultured cells, per SD
  of normalised expression. About 20 % of instruments are missing or allele-mismatched in the GWAS files.

## 33. Revised classifiers against baselines (internal validation only)

Previous classifier values and permutation P values are superseded. GEO specifies one biopsy per subject without replicates in GSE89408, and all 218 GEO sample IDs are unique. So for this deposited design, biopsy-level splitting is subject-level (logs/gse89408_design_verification.json). The revised script takes counts into every fold. Filtering, normalisation reference, variable-gene selection, means/SDs, imputation, scaling and tuning use only the relevant training fold. Observed F1 and 199 permuted labels use the same 5-fold outer CV repeated 5 times and 3-fold inner tuning. Each of 100 random controls uses a fixed gene ranking. It uses the same number of training-expressed genes as were retained from the fixed F2 list. Fold SD is not a confidence interval. Final values are saved in classifier_baselines.tsv, classifier_permutations.tsv and classifier_protocol.json.

These CV and permutation estimates are conditional on the specified feature lists. They do not correct for earlier exploratory choices of tasks or programmes.

The separate classifier_F1_coefficients.tsv table is a descriptive full-sample F1 fit at fixed C=1, after the transformations were trained on the complete task cohort. It is not a held-out importance estimate, an average of outer-fold coefficients or the source of the reported CV AUCs.

| Task | Model | Mean AUC | SD |
|---|---|---:|---:|
| RA vs OA | F0 covariates | 0.574 | 0.143 |
| RA vs OA | F1 eight programme scores | 0.900 | 0.048 |
| RA vs OA | F2 top-30 prioritised genes | 0.933 | 0.042 |
| RA vs OA | F3 genome-wide elastic net | 0.979 | 0.029 |
| RA vs OA | F4 random 30-gene sets (median of 100) | 0.909 | 0.046 |
| ACPA+ vs ACPA- (RA) | F0 covariates | 0.332 | 0.079 |
| ACPA+ vs ACPA- (RA) | F1 eight programme scores | 0.797 | 0.062 |
| ACPA+ vs ACPA- (RA) | F2 top-30 prioritised genes | 0.647 | 0.083 |
| ACPA+ vs ACPA- (RA) | F3 genome-wide elastic net | 0.726 | 0.089 |
| ACPA+ vs ACPA- (RA) | F4 random 30-gene sets (median of 100) | 0.649 | 0.054 |

AUC SD is across correlated outer folds; random-set SD is across draw means. Neither is a confidence interval.

- **RA vs OA**: F1 permutation p 0.005; same repeated nested-CV procedure. F2 vs random, empirical p 0.2475; SD across random-set mean AUCs. F1 − genome-wide: -0.079; t -3.017, p 0.00596; approximate corrected resampled test.
- **ACPA+ vs ACPA- (RA)**: F1 permutation p 0.005; same repeated nested-CV procedure. F2 vs random, empirical p 0.5149; SD across random-set mean AUCs. F1 − genome-wide: 0.072; t 1.49, p 0.1492; approximate corrected resampled test.

All estimates are internal to GSE89408. A low label-permutation P shows an association under this modelling procedure. It does not show external predictive performance, or that the candidate biology is specific.


## 34. Second RA single-cell cohort: synovial fluid before and on treatment (GSE296117) (**exploratory**)

`scripts/40_ra_fluid_treatment.py`; programmes from `scripts/state_markers.py` (`results/tables/state_markers_gse243917.tsv`);
`results/tables/ra_fluid_tests.tsv`, `ra_fluid_m6_cluster_markers.tsv`; `data/processed/ra_fluid_unit_scores.tsv`;
Extended Data Fig. 2a–c.

- **Programmes, fixed before any new cohort was scored.**
  - M-6 programme: M-6 STAT1+ CXCL10+ vs the other macrophage and monocyte subtypes of GSE243917 (65 vs 5,784 nuclei).
    Its 30 genes include GBP1, CXCL10, STAT1, PARP14, IRF1, GBP5, APOL6, SOD2, IL4I1 and CXCL9.
  - Tph programme: T-7 Tph vs the other CD4 T-cell subtypes (T-3 Tfh/Tph excluded; 285 vs 4,367 nuclei). Its 30 genes
    include TOX, MAF, CTLA4, IL6ST, CXCL13 and TOX2.
  - Rule: log2 fold change > 1, detected in >= 10 % of the state's nuclei, no mitochondrial or ribosomal genes, top 30
    by Wilcoxon score.
- **Cohort.** 100,387 synovial-fluid cells from 9 RA patients, each sampled treatment-naive and on treatment (5 adalimumab,
  4 tofacitinib); author cell types (36,820 macrophages). The samples are fluid, not tissue, and there are no clinical response labels.

| Test (fixed in the script header) | Result | Verdict |
|---|---|---|
| V1 M-6-like cluster among naive macrophages | best cluster 4.3 % of 22,091 cells, present in 7/9 patients, 9/30 top markers are M-6 genes; vs 200 random sets P = 0.055 | borderline |
| V2 macrophage M-6 score on treatment (primary) | fell in 6/9, median change −0.51, paired P = 0.13; larger drop than 200 matched random gene sets (P = 0.005) | specific but not consistent across patients |
| V2 active-RA signature / Hallmark IFN-γ | paired P = 0.82 / 0.20 | not supported |
| V3 share of M-6-like macrophages | median 5.6 % → 2.9 %, fell in 7/9, P = 0.027 | supported |
| V4 Tph score of T cells | fell in 7/8, P = 0.055 | borderline |

- **Reading.** The M-6 state recurs in RA synovial fluid and becomes rarer on treatment. This is in line with the anti-TNF
  reversal in bulk tissue (sections 6–7). By drug, the median M-6 change was −0.83 with adalimumab (n = 5) and −0.41 with
  tofacitinib (n = 4). Four patients are too few to test the stronger JAK-inhibitor prediction.
- **Caveats.** Nine patients; fluid cells are not tissue cells; time on treatment and clinical response are not
  deposited; author annotations are used as deposited.

## 35. NSCLC extension of the RA programmes (LuCA atlas, TCGA) (**ready for the state (N1–N2); coupling exploratory; bulk survival negative**)

`scripts/41_nsclc_extension.py`; `results/tables/nsclc_extension_tests.tsv`, `nsclc_extension_cox.tsv`,
`luca_m6_cluster_markers.tsv`; `data/processed/luca_unit_scores.tsv`; Extended Data Fig. 2d–h.

- **Data.** LuCA core atlas (Salcher et al. 2022; 19 studies, raw counts), lung tissue only: 153,743 non-alveolar
  macrophages and CD4 T cells from primary tumour, adjacent normal and normal lung. TCGA LUAD and LUSC bulk RNA-seq and
  curated survival, as in section 14. The programmes are the same as in section 34 and were not re-fitted.

| Test (fixed in the script header) | Result | Verdict |
|---|---|---|
| N1 M-6-like state among primary-tumour macrophages | best cluster 7.2 % of 39,260 cells, >= 1 % of macrophages in 74/131 donors; top markers CXCL9, GBP1, STAT1, CALHM6, CXCL10 (14/30 are M-6 genes); vs 200 random sets P = 0.005 | supported |
| N2 macrophage M-6 score, tumour vs adjacent normal (29 paired donors) | higher in tumour in 20/29, paired P = 0.004; vs random sets P = 0.005 | supported |
| N3 macrophage M-6 vs CD4 T-cell Tph across tumours (study-adjusted; 52 tumours, 10 studies) | ρ = 0.24, P = 0.09; higher than 97 % of random macrophage gene sets (P = 0.03) | suggestive |
| N4 TCGA survival (Cox per SD; age, stage) | M-6: LUAD HR 0.99 (P = 0.84), LUSC 1.04 (P = 0.57); Tph: 0.92 (P = 0.24), 1.07 (P = 0.32). With leukocyte adjustment (PTPRC): LUAD M-6 HR 1.59 (1.25–2.03, P = 2e-4), others P >= 0.07 | no association; adjusted LUAD signal opposite to expectation |
| N4 TCGA tumour vs adjacent normal (bulk) | M-6 lower in LUSC tumours (P < 1e-4; LUAD P = 0.10); Tph lower in both (P <= 1e-4) | bulk dilution; not comparable with N2 |

- **Reading.** The IFN-γ-activated macrophage state defined in RA synovium is a recognisable state in NSCLC tumours. It is
  enriched over adjacent lung. Its conventional correlation with Tph-like CD4 T-cell scores is only modest and not significant (P=0.092). Its prognostic meaning in bulk is unclear. Alone it is not associated with survival. After leukocyte adjustment,
  the LUAD association points to worse survival, opposite to the CXCL9:SPP1 polarity of section 14. Several M-6 genes
  (STAT1, B2M, HLA class I) are also expressed by tumour cells, so the bulk score is not macrophage-specific.
- **Caveats.** The atlas integrates 19 studies and platforms (handled by within-donor pairing for N2 and within-study
  centring for N3). Adjacent normal lung may carry field effects. TCGA is bulk.

## 36. Changes made during the analysis (for transparency)

1. The HVG flavour was switched from `seurat_v3` to `seurat`, because the LOESS fit fails on the 137-nucleus sample.
2. Marker sets were revised after the first E-MTAB-8322 pass mis-called FCN1+ LYZ+ monocytes as endothelial
   (PECAM1 is expressed by monocytes). Myeloid markers now cover monocytes (CD14, LYZ, TYROBP, FCER1G, AIF1), and
   PECAM1 was removed from the endothelial set. The GSE243917 annotation was re-run and did not change (ARI 0.952).
3. The E-MTAB-8322 analysis was split into two stages after ~47 % of cells were found to be non-myeloid.
4. "MERTK+" clusters were first defined by a MERTK-detection threshold. This mislabelled the LYVE1-high
   SELENOP+ resident cluster as MERTK-negative, because MERTK is detected in <= 21 % of cells. It was replaced by the
   dominant a priori program (homeostatic vs inflammatory). The change was made after the first group tests had
   been seen, but before the classifier was run. The active-vs-healthy result holds under both definitions.
5. Signature transfer: the declared fallback (nominal top-50 per direction) was kept as a sensitivity
   analysis. The primary signature uses only the FDR-significant direction. Donors with < 50 macrophages were excluded.
6. ATAC step: muon's `tl.leiden` clusters each modality separately, so the joint WNN graph is clustered with scanpy directly.
   The random-gene-set control for the donor-level RNA-ATAC agreement was added after the first run showed
   rho = 0.84. Its purpose is to rule out a shared nucleus-quality effect.
7. GWAS step: the old GWAS Catalog API is no longer available, so v2 is used. The LDL download timed out once, and the
   script was made resumable (per-page cache). The RA-only trait filter and all enrichment parameters (N = 4,000
   peaks, >= 100 nuclei, windows, null matching) were fixed before the enrichment was run.
8. Clinical-response step: R4RA was replaced by GSE198520 (reasons in section 6). The first run showed a dominant
   library-size axis (random gene sets AUC ~0.15). After that, CPM was switched to nuclear genes only. Technical predictors,
   a technical-residual score and covariate-adjusted models were also added. The first-run numbers (raw score AUC 0.70) are in
   the run log history and lead to the same conclusion.
9. Steps 12–15: single-nucleus chromVAR was replaced by pseudobulk units (memory). LIANA's magnitude ranking was
   supplemented with a specificity view after the first output was dominated by ECM -> CD44 pairs. The knowledge-mining
   step was corrected twice before interpretation. RA stage is now indication-specific (maraviroc was first mislabelled
   phase 4 in RA), and drug action direction was added (ipilimumab and interferon gamma-1b were first listed as leads).
10. Steps 16–19: step 09 was extended to keep GWAS risk alleles (RA was re-queried with 50-row pages because the API timed
   out on larger pages). In step 18 only 3 signature genes passed the monocyte expression filter, so the Hallmark
   IFN-gamma response was used alongside it. Step 19 p-values are written in full precision. The RA re-query on
   2026-10-19 returned the same 3,221 associations as on 2026-10-18 (identical association IDs), so step 10 is unaffected.
11. Steps 20–23: only per-sample matrices and cell tables were downloaded from GSE347958 (the 57 GB archive was skipped).
   Proteomics donor means drop missing values (one NH donor lacks S100A8/9/ICAM1). The knowledge-graph viewer was fixed
   after a first test (search did not focus while the layout was still moving).
12. Steps 24–27: the eQTL Catalogue REST API is no longer available, so remote tabix reads of its FTP files are used (one dataset
   lacks an rsid column; this is handled). L1000 hits needed a compound-name class dictionary, because ChEMBL has no mechanism
   for many tool compounds. A random-gene-list baseline was added before the classes were interpreted. Deconvolution was
   tried once more with per-gene scaling after the first version failed sanity checks, and then stopped.
13. Steps 28–31: step 25 now caches remote tabix indexes in `data/raw/eqtl/` (they had been written to the project root).
   Integration: methods, inputs, MOFA+ settings and metrics were fixed before the run. The heatmap text colour was changed
   for legibility and the step was re-run (MOFA+ cached; numbers unchanged).
   Fine-mapping: Open Targets rejects 25 credible sets per request as too expensive, so 10 are fetched per page. The
   per-state count of contributing credible sets and the share of the largest were added after the first run showed one
   variant dominating. Disease stage: CPM was replaced by median-of-ratios size factors after random gene sets came out
   higher in RA than OA. The random-set comparison was then extended to every test. Macrophage-abundance adjustment was extended to
   every comparison after the homeostatic programme rose in RA vs healthy. Trajectory: the chromatin score uses gene
   activity, because only 277 peak-gene links pass FDR and none fall on IFN genes (decided before any test). The
   M-7 -> M-6 path first used shortest paths on -log(connectivity), which failed on connectivities of 1. It is now taken
   from the PAGA spanning tree (same path). The half-max crossing was first defined as the first bin >= 0.5. That picks
   the starting bin when a curve begins high. It was redefined as the crossing of the final rise to the maximum after
   the first result had been seen (lag -0.32 -> +0.26), and a bins x threshold grid plus selection controls were added.
   A check then showed that M-7 nuclei come out latest in M-7-rooted pseudotime. Roots from smoothed scores did not fix
   it, so T3 was declared not interpretable and no further variants were tried.
14. Steps 32–34 and the East Asian replication: superseded backups and exploration probes from earlier rounds were
   deleted (no stale project outputs existed). Colocalisation: the first run counted trans-pQTL colocalisations
   (UKB-PPP) and mixed splicing / transcript-usage directions. The cis restriction and one direction call per dataset
   from expression and protein QTLs were added before any result was interpreted. TF activity: OmniPath's UniProt
   columns collided with the renamed gene-symbol columns (fixed). NF-kB dimer regulons are hidden in the figure only.
   S-LDSC: polyfun's LDSC has no --h2-cts, so one --h2 model is fitted per state. It needed a second environment
   (Python 3.10, numpy 1.26, pandas 1.5). The first LD-score run was stopped because parallel jobs oversubscribed CPU
   threads. Regressions first failed because polyfun requires CHR and BP columns in the summary statistics. The all-peaks
   summary was added after the first regression output was seen. The East Asian study was added to step 29 as F4.
15. Steps 35–38: MAGMA first failed on chromosome 6 (an MHC p-value of 2.6e-319). The MHC was excluded and p-values were
   floored at 1e-300. Because 261 non-MHC SNPs were affected by the floor, all chromosomes were re-run. The extended
   MHC (chr6:25–34 Mb) was removed from the scDRS gene set after flanking genes (ITPR3, ZSCAN31) topped the MAGMA list.
   scDRS 1.0.2 needed a NumPy 2 alias (np.float_). MR instruments were matched to the step-32 records through the
   rounded H4/CLPP values written by step 32. The classifier's coefficient table and the MR plot's consistency labels
   were added after the first run; no model or test was changed.
16. Full clean rerun for submission (2026-10-19/20). All intermediate files and results were deleted. Every step was then
   re-run from the downloaded raw data with `scripts/run_all.sh`. Each table was compared with the previous run using
   `scripts/compare_results.py` (`results/reproducibility_check.tsv`). The rerun found and fixed:
   - **Three sources of run-to-run randomness.**
     - muon's ATAC LSI calls an unseeded truncated SVD. It is replaced by a seeded SVD with a fixed sign convention.
       Steps 06–07 then gave identical outputs in two consecutive runs.
     - MAGMA's gene-level permutations now use `--seed 1`. Two runs of chromosome 16 gave byte-identical output.
     - The saga solver in step 38 is now seeded.
   - **Step 25.** It failed three times on interrupted remote reads of eQTL Catalogue files. It now:
     - reads the GWAS from the local copy;
     - retries remote reads with back-off and caches every region it reads (`data/raw/eqtl/regions/`);
     - stops instead of skipping a dataset after repeated errors.
     A merge also lost the shared-variant label. It now reports rsIDs where the dataset has them (previously
     `chr_pos_ref_alt`).
   - **Run order.** Step 25 needs the first pass of step 22, and step 38 needs the final one. Steps 22 and 32 run twice
     (22 reads 32's colocalisations, and 32 is annotated with 22's ranking).
   - **Environment paths and caches.** The LDSC interpreter is found next to the active environment or via `LDSC_PYTHON`.
     The cached MOFA+ model and LIANA results were deleted and recomputed. Both reproduced exactly.

   **Outcome: 93 of 104 tables are identical; 11 differ.**

   | Tables that differ | What changed | Why |
   |---|---|---|
   | integration_benchmark | WNN no longer leads on subtypes and ATAC clusters; overall ranking held (section 23) | LSI was unseeded |
   | atac_label_transfer, wnn_recovery, wnn_modality_weights | ≤ 0.03 per row (cell types: ≤ 0.002); WNN ARI 0.412 → 0.414 | LSI was unseeded |
   | magma_genes (Ishigaki, T2D), scdrs_groups | a few gene p-values; top-1,000 gene sets unchanged; scDRS group z ≤ 0.013, one quoted value −1.3 → −1.2 | MAGMA was unseeded |
   | classifier_baselines, classifier_folds | ACPA genome-wide model AUC 0.7156 → 0.7161; no quoted value changed | saga was unseeded |
   | eqtl_coloc | variant label format only; all statistics and lead positions identical | label fix above |
   | coloc_direction | evidence_score only | 32 now runs after the final 22 |

   No conclusion in this report changed.

17. Additions on 2026-10-20 after the full rerun:
   - **Reproducibility.** Exact package locks of both environments (`env/`). `scripts/verify_quick.sh` checks the manifest
     checksums and input checksums, and re-runs steps 11, 12, 23 and 30 to compare them with the packaged tables. It passes on a fresh copy of the archive. Step 17
     was in the first version, but it reads an intermediate ATAC file that the archive does not contain. It was replaced by
     step 12 after a test on a fresh copy. A
     `GENEDIS_FIG_DIR` option was added to `config.py` so that the check does not write into the figure folder. `macOS grep` does
     not support `\|`, which stopped the first test of the check. It now uses awk.
   - **Figures.** The main figures were redesigned (Nature-style layout; `39_main_figures.py`
     rewritten; the old per-step images were dropped and left out of the package; this version was kept).
     Extended Data Fig. 2 was then added for steps 40–41 in the same style. Three text overlaps were fixed after the
     first render. Earlier the same day, a panel letter overlapping a title in Fig. 2 had been moved.
   - **New cohorts (sections 34–35).** Hypotheses were written in each script header before scoring, and the
     programmes were derived once in GSE243917. A second RA tissue cohort was searched for first. CELLxGENE has no RA
     synovium dataset. The only open RA synovium scRNA-seq dataset on GEO with macrophage clusters (GSE246416)
     deposits only its healthy donors. So the synovial-fluid cohort GSE296117 was used. Step 40 failed once on an integer
     index name (code fix). It was re-run after the per-patient share was added to its output table (test results
     identical). `data/raw/SHA256SUMS.txt` was recomputed: no earlier checksum changed, and the new entries are the two datasets
     and the cached eQTL regions of step 25.

## 37. Limitations

- Small donor numbers (11 RA with histology; 19 RA with treatment status). Most results are exploratory.
- The healthy-vs-RA comparison is confounded by age and sex.
- GSE243917 has an inferred nuclei-prep batch that aligns with macrophage-state features.
- snRNA (GSE243917) and scRNA (E-MTAB-8322) differ in transcript capture (nuclear vs whole-cell).
- Treatment-status labels are used as deposited. No DAS28 or drug-level detail is available in the public metadata.
- ATAC: only a peak x cell matrix is deposited (no fragments), so TSS enrichment, FRiP and footprinting are not
  possible. RNA and ATAC come from the same nuclei, so cross-layer agreement is not independent replication.
- GWAS enrichment (section 5) uses lead SNPs with fixed windows. The fine-mapped version (section 24) rests on a few
  high-PIP loci (CD2, ANKRD55, ARID5B). S-LDSC (section 28) uses a 1000 Genomes European LD panel only.
- Direction of effect (section 27) comes from QTLs measured in blood and cultured immune cells, not synovium.
- Clinical-response evidence comes from one anti-TNF cohort (n = 33 for good vs non), and its response labels are
  entangled with library size. Rituximab/tocilizumab response (R4RA) could not be tested because the data are not public.
- Cell-cell communication and TF-motif results for the STAT1+ CXCL10+ state rest on 65 nuclei (2 donors with enough
  nuclei for paired motif analysis). Open Targets content changes between releases (26.09 recorded).
- Blood validation uses one PBMC cohort (13 complete RA-control pairs). Only 3 signature genes are measurable in monocytes.
- NSCLC uses bulk TCGA data. SPP1 is also expressed by tumour cells, and the survival association replicates published
  pan-cancer work.
- Variant-to-motif effects are sequence predictions on lead SNPs without allelic validation. The ANKRD55 SNP of section 11
  (rs6859219) lies outside the fine-mapped credible sets, whose single causal candidate is rs7731626.
- The active-RA signature's FDR significance depends largely on one donor (SA139). Its direction is stable.
- Xenium: 4 RA sections on the 5K panel (one does not support H1); the 480-gene panel lacks IFNG; OA blocks are from one
  donor. Proteomics: 6 vs 3 donors, and few signature proteins are detectable.
- Gene prioritisation uses equal weights and nearest-gene assignment. The knowledge graph only contains results that
  passed each step's thresholds.
- Visium spots are multicellular; L1000 uses cell lines; deconvolution failed; coloc assumes one causal variant per
  locus; the IFN-gamma macrophage eQTL dataset has 84 donors.
- **Unavailable here:** a technically clean responder cohort with public processed data (R4RA would qualify if the
  authors share processed counts and arms);
  genetic testing of the STAT1+ CXCL10+ macrophage state (too few nuclei); MultiVI (no GPU).
- Bulk synovium (GSE89408): groups differ on random genes too, so only differences beyond random gene sets are interpreted.
  Healthy biopsies come from younger donors (median 31 vs 54 years).
- Tph/Treg: Xenium marks Tph by CXCL13 alone (Tfh-like cells share it). IL21 is too sparse in snRNA to test.
- scDRS tests expression of risk genes, not where risk variants act. MR uses one instrument per signal (no pleiotropy
  test) and QTL effects from blood and cultured cells. IL6ST and ANKRD55 share their instrument.
- Classifier: random genes also separate RA from OA (group-wide shift), so only comparisons with baselines are meaningful.
- Trajectory: PAGA and pseudotime order states but do not prove lineage. M-6 has 65 nuclei. The deposited count matrix
  has no spliced / unspliced split, so RNA velocity is not possible.
- Second RA cohort (section 34): synovial fluid from 9 patients without response labels. NSCLC (section 35) uses an
  integrated multi-study atlas and bulk TCGA. GSE296117 needs R with Seurat for a one-off conversion.
