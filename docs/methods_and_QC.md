# Methods and quality control

**Supplement to the GeNeDis 2026 RA multi-omics technical report**

Xia Wu | PhD student, University of Iowa | xwu76@uiowa.edu

## A1. Cohort structure and preprocessing

### Analysis units and sample handling

| Component | Processing / analysis unit | QC finding or restriction |
|---|---|---|
| GSE243917 RNA/ATAC | Nucleus for embedding; donor for pseudobulk and clinical association | 31,547 deposited nuclei; 27 doublet flags removed, leaving 31,520. Twelve donors include 11 RA and one OA. |
| E-MTAB-8322 scRNA | Donor macrophage pseudobulk for active/remission DE | 81,601 cells pass cell QC; 42,963 are broadly myeloid. Non-myeloid cells are excluded; DC-like clusters are excluded from macrophage features. |
| GSE89408 bulk | One biopsy per subject; no repeated samples in model tasks | 218 biopsies across healthy, arthralgia, UA, early/established RA and OA. Labels determine task-specific inclusion. |
| GSE198520 bulk | Paired patient before/after treatment | 46 patients, 92 biopsies. Low depth and variable mitochondrial fractions affect score coverage and baseline prediction. |
| Xenium / Visium / protein | Section-specific spatial tests; donor summaries for disease comparisons | Serial/shared blocks and cohort overlap are tracked. Regions are pooled within donor before protein tests. |

The scRNA discovery cohort contains 27 samples from 26 donors after merging SA149P with SA149. The active/remission contrast has 12 active and seven remission donors. PyDESeq2 models donor pseudobulk using `~ sort_panel + condition`. The primary signature contains 15 symbol-mapped, active-up FDR genes; the full DE table retains additional Ensembl identifiers. Discovery results are sensitive to removing donor SA139: only one of the 15 signature genes remains significant.

Broad cell annotation is checked against author annotations (ARI 0.952; 97.3% agreement). Inferred nuclei-preparation batches have different mitochondrial fractions and cell composition. They are carried into relevant comparisons. BRI-1588 has only 137 nuclei and four macrophages and is excluded where a donor pseudobulk requires sufficient cells. Healthy/active comparisons in the sorted-cell cohort also have age/sex imbalance.

Gene identifiers and measured coverage are recorded by cohort. Unmeasured genes are not treated as biological zeros or imputed as a full signature. Raw references, source accessions, selection criteria, labels and access terms are listed in `docs/sources.tsv`; QC tables remain available at sample level.

<!-- PAGEBREAK -->

## A2. Score definitions and assay coverage

| Score | Definition | Uses | Interpretation restriction |
|---|---|---|---|
| S1 | 15 named active-up discovery genes | Macrophage transfer, bulk disease tests and treatment-context analysis | Only 13/15 measured in atlas and anti-TNF RNA; independent inflammation/bulk tests are negative. |
| S2 | Ten-gene IFN programme | Xenium, Visium and bulk-stage analyses | Shares only APOL6 and GBP5 with S1; cannot validate S1 by substitution. |
| Hallmark IFNG | Measurable genes from HALLMARK_INTERFERON_GAMMA_RESPONSE | Pathway and spatial-protein analyses | Assay-specific coverage; distinct from exact S1/S2 gene sets. |
| M-6 programme | Thirty top markers of the author-defined atlas state | Synovial-fluid and NSCLC transfer | Exploratory state score; not a new cell type or the S1 discovery set. |
| Tph programme | Thirty atlas-state markers | Synovial-fluid / lung coupling tests | Distinct from the short Tph feature in the bulk classifier. |

S1: MT1H, MT1G, MT2A, MT1F, MT1L, MT1X, PKM, APOL6, MT1E, PFKP, MT1M, PALM2AKAP2, IL4I1, SLAMF7 and GBP5.

S2: CXCL9, CXCL10, CXCL11, STAT1, GBP1, GBP5, IDO1, IRF1, CD274 and APOL6. All programme lists, their roles and limitations are supplied in `docs/evidence_map.tsv` and `results/tables/state_markers_gse243917.tsv`.

For the independent atlas transfer, macrophage pseudobulk requires at least 50 macrophages per donor, leaving ten RA donors. RNA/ATAC score concordance describes two measurements of the same nuclei. No S1 gene has an FDR-significant peak-gene link under this analysis's cis-link criteria. In bulk, S1's RA/OA association and matched-random specificity tests are negative. S2's bulk comparison remains supported against matched random genes and after macrophage-abundance adjustment.

In the anti-TNF cohort, RNA scores use log2(CPM+1) and baseline-derived gene standardisation. A descriptive technical-residual score is fitted label-free against library size and mitochondrial fraction across all 92 samples. This is not a prospective held-out predictor. The good-responder within-group change and good/nonresponder difference are reported separately. Apparent baseline prediction is technically confounded and is not presented as a validated model.

For proteomics, lining and all immune/fibrous/other sublining regions are first averaged within donor. Each compartment comparison has six RA and three healthy donor means. The Hallmark score differs at P=0.02381 in each compartment. Random protein-set comparisons support lining specificity (0/1,000 as extreme; resolution about 0.001) but not sublining specificity (P=0.339).

<!-- PAGEBREAK -->

## A3. Supervised model protocol

### Tasks and baselines

The GSE89408 metadata specify one biopsy per subject. Task 1 compares RA with OA (n=174; 152 RA). Task 2 compares ACPA-positive with ACPA-negative RA (n=139; 96 positive). Samples with unavailable task labels are omitted. No external test cohort is claimed.

F0 uses age, sex and log library size as covariates. F1 uses eight specified programme scores: S1, IFN, homeostatic, macrophage, short Tph, Treg, T-cell and B/plasma features. F2 uses a fixed 30-gene prioritised list, retaining the genes that pass training-fold filtering. F3 is elastic-net logistic regression using up to 2,000 training-selected variable genes. All models use class balancing.

### Nested training and evaluation

- Outer evaluation: stratified five-fold CV, repeated five times with a recorded seed. Each model has 25 held-out AUCs.
- Inner selection: three stratified folds. Inverse regularisation C is chosen from 0.01, 0.1 and 1.0. The elastic-net comparator also tunes its L1 ratio over 0.1, 0.5 and 0.9.
- Count filtering: at least 1 CPM in at least 20% of training observations. Median-of-ratios reference genes and geometric means are fitted using training counts.
- Transformation: held-out counts use the fitted training reference. Gene scaling, programme aggregation, age imputation and feature standardisation are fitted within the relevant training split, including inner splits.
- Variable genes: ranked by variance in the training data after normalisation. No held-out observation selects variable genes.
- Convergence: logistic-regression warnings trigger explicit higher-iteration retries; unresolved convergence stops the analysis.

### Null comparisons and uncertainty

For F1, each of 199 label permutations repeats the same preprocessing, outer repeats, inner tuning and summary statistic as the observed model. With a plus-one empirical correction, the minimum P is 0.005. Tests are exploratory and per task.

F2 is compared with 100 fixed random gene rankings. Within every inner and outer training split, random controls match the number of prioritised genes retained after filtering. The random-set comparison uses mean AUC across the same outer folds. The F2 random-set P values are 0.2475 and 0.5149; neither supports prioritised-gene specificity.

Fold SD describes correlated folds and is not a confidence interval. F1/F3 comparisons use an approximate corrected resampled t-test and do not substitute for independent validation. Feature lists and tasks reflect earlier exploratory choices; nested CV does not correct those earlier choices. Full-sample F1 coefficients at fixed C=1 are descriptive and separate from held-out performance.

Primary records: `classifier_protocol.json`, `classifier_baselines.tsv`, `classifier_folds.tsv`, `classifier_permutations.tsv` and `classifier_random_sets.tsv` in `results/tables/`.

<!-- PAGEBREAK -->

## A4. Integration and association methods

### Paired RNA/ATAC benchmark

Shared cells are compared using RNA, ATAC, concatenation plus Harmony, WNN and one MOFA+ configuration. Metrics cover broad/fine author labels, ATAC-derived clustering, neighbour-label agreement, rare M-6 purity, batch entropy/connectivity and agreement with each unimodal graph. The supplied weighted summary is a benchmark of these configurations, not a universal ranking of methods. Original author-state labels remain the reference.

### Genetic and regulatory evidence

Conditional S-LDSC uses baselineLD annotations and the frozen full summary statistics. The final family contains 195 models: 65 annotations for each of Ishigaki European RA, FinnGen RA and the T2D control. Any failed, missing, duplicated or non-finite regression stops the step before BH adjustment. Adjustment is checked within GWAS-by-annotation-level families. Native Parquet caches preserve exact baseline annotation values/dtypes; original compressed inputs are retained.

scDRS and S-LDSC use complementary expression/chromatin measurements but share GWAS traits. GWAS participant independence across resources has not been established. Credible-set overlap, colocalisation and gene ranking are separate evidence types. Overlapping ranking layers are not independent probabilities.

RNA TF screening applies BH to all 464 gene-label-permutation tests, followed by the positive/top-decile motif criterion. M-6 regulatory units come from two donors. Gene-label permutations do not create donor replication. IRF5 has genetic/colocalisation evidence but is not a corrected RNA candidate in this screen (RNA q about 0.139).

Single-instrument Wald ratios cannot assess horizontal pleiotropy. Shared IL6ST/ANKRD55 instruments do not distinguish their effects. IRF5 directions differ across signals/datasets. CTLA4 expression is not equivalent to CTLA4-Ig pharmacology, and IL6ST/gp130 is not IL6R. Drug annotations give biological context only.

### Spatial empirical tests

Within each 5K Xenium section, the macrophage IFN score is regressed on the numbers of IFNG-positive and IFNG-negative T/NK neighbours within 30 micrometres. IFNG labels are shuffled among lymphocytes 1,000 times at fixed positions. Effects and empirical P values are reported by section; there is no normal-tail extrapolation or combined extreme Stouffer P. A healthy section with no IFNG-positive exposure is marked non-estimable.

CXCL13/B-cell and IFNG/macrophage spatial results are neighbourhood associations. Global label shuffles do not preserve expression-label clustering or tissue compartments. Shared blocks, targeted panels and segmentation/annotation uncertainty limit both independence and generalisation.

<!-- PAGEBREAK -->

## A5. Reproduction and result availability

### Re-run and integrity record

The tested full-analysis platform is macOS arm64 on an Apple M4 Max with 48 GB RAM. Python environments and package versions are pinned; the R/Seurat conversion session and external-tool hashes are recorded. MAGMA uses the tested macOS x86_64 build through Rosetta 2. The supplied downloader/tool-hash check is macOS-specific; Linux requires an appropriate MAGMA binary/download adaptation. Cross-platform numerical identity has not been certified.

The scientific validation began with empty analysis intermediates/results, byte-verified public primary inputs and frozen reference inputs. All 43 ordered invocations completed. Independent branches were scheduled within the same clean project; any reuse was checked against exact code/input/output hashes. The recorded 2.78-hour elapsed run is platform-specific. 

The frozen reference archive preserves API/GMT/eQTL inputs, including Open Targets 26.09, with an archive hash and 1,904 per-file checks. Primary raw inputs have a 391-file hash record. The complete submission ZIP includes the frozen references but omits large primary raw files and intermediate models. A Git checkout restores the same pinned reference asset; large public inputs are fetched by the documented downloader.

Install the two analysis environments, then use `bash scripts/verify_quick.sh` for a non-mutating extracted-package check. Full reconstruction uses `bash scripts/00_download.sh` followed by `bash scripts/run_all.sh`. Quick verification previously reproduced 112 result and 15 processed tables and regenerated all six figures in an isolated directory. The manifest records the size and checksum of each archived file.

### Available results

| Output | Status | Main qualification |
|---|---|---|
| Processed cohorts, integration benchmarks and source-linked figures | Available | Source labels, QC exclusions and assay coverage retained |
| Cross-cohort S1 validation | Exploratory / negative | Inflammation and bulk specificity tests do not pass |
| Spatial IFN/CXCL13 associations | Exploratory | Section-level empirical null; compartment/cluster structure omitted |
| GSE89408 classifier comparison | Available internal benchmark | No external clinical validation; conditional on explored features/tasks |
| Anti-TNF prediction, blood transfer and deconvolution | Unsuccessful or uninterpretable | Depth/technical confounding, failed transfer or sanity checks |
| Chromatin-before-RNA timing and drug efficacy | Unavailable | Cross-sectional associations and single instruments cannot establish them |
| Synovial-fluid and NSCLC transfer | Exploratory | Pooled therapies, unavailable response labels and negative primary survival tests |

### Data and code locations

`docs/sources.tsv` holds provenance/access terms. `docs/evidence_map.tsv` defines scores. `data/processed/` and `results/tables/` contain features/results. `results/figures/main/` contains six original figure sets and source hashes; `presentation/GeNeDis2026_figure_atlas.pdf` is the separate figure atlas. Environments, raw/reference manifests and `logs/` support reconstruction. `docs/REVISION_VALIDATION.md` and the retained historical QC record document analysis history; they record the original analyses and subsequent corrections.
