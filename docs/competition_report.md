# A reproducible multi-omics workflow for rheumatoid arthritis

**Technical report | GeNeDis 2026 Hackathon**

Xia Wu | PhD student, University of Iowa | xwu76@uiowa.edu

## 1. Project overview

### Summary

Rheumatoid arthritis (RA) datasets describe genetic susceptibility, synovial cell states and tissue inflammation at different biological scales. This project builds a reproducible workflow that connects those measurements through documented cohort features, explicit gene-set definitions and baseline comparisons. It combines a published RNA/ATAC atlas with macrophage single-cell RNA sequencing, bulk synovial cohorts, spatial transcriptomics, tissue proteomics and public genetic resources. The clearest convergent pattern is CD4-associated genetic enrichment alongside an interferon-associated tissue signal. A separate macrophage discovery signature transfers only partly: its independent inflammation and bulk disease tests are negative. Internally evaluated classifiers separate disease groups, but prioritised genes do not outperform matched random gene sets. The deliverable is a reusable analysis and evidence resource for follow-up studies, with the limits of each association retained.

### Question and scope

Which findings remain supported when RA genetic risk, macrophage state, tissue neighbourhood and clinical-context measurements are examined together? The primary scope is RA synovial tissue. Synovial-fluid treatment and non-small cell lung cancer (NSCLC) analyses test whether specified atlas-state programmes can be transferred to additional datasets.

### What this submission contributes

- **A documented data structure:** source/accession and access-term records, cohort labels, processed feature tables and reference checksums.
- **Comparable analyses:** RNA/ATAC integration benchmarks, donor-level signature tests and nested cross-validation against simple and genome-wide models.
- **Traceable evidence:** separate score definitions, assay coverage, negative controls and links from results to figures and a knowledge graph.

Published cell-state annotations follow Gupta et al. [1], and the macrophage discovery cohort follows Alivernini et al. [2]. This work's contribution is the implemented integration, comparison and validation workflow built on those resources.

{{WORKFLOW}}

### Material to read with this report

The separate **Methods and QC appendix** provides preprocessing, sample-unit and model details. The **Figure atlas** contains four main and two Extended Data figures, with panel captions and source data. Code, processed features and full results are in the reproducibility archive. The slide deck is optional supporting material.

<!-- PAGEBREAK -->

## 2. Data and analytical design

### Core resources

| Resource | Measurements | Analysis units | Main role |
|---|---|---|---|
| GSE243917 [1] | Paired nuclear RNA and ATAC | 31,520 retained nuclei; 12 donors, including 11 RA and 1 OA | Cell-state reference and integration benchmark |
| E-MTAB-8322 [2] | Sorted-tissue scRNA-seq | 26 donors; 12 active and 7 remission RA in discovery contrast | Macrophage features and S1 signature |
| GSE89408 [4] | Bulk synovial RNA-seq | 218 biopsies; 174 in RA/OA and 139 in ACPA tasks | Disease-stage tests and model comparison |
| GSE198520 [3] | Bulk synovial RNA-seq | 46 patients sampled before and after anti-TNF | Paired treatment-context analysis |
| GSE347958 / GSE347797 [5] | Xenium / Visium | Four RA 5K sections for neighbourhood tests; donor-level tissue comparisons | Spatial distribution of cells and IFN scores |
| PXD062696 [6] | Region-resolved proteomics | Six RA and three healthy donors in compartment comparisons | Protein-level pathway evidence |
| Public GWAS and QTL resources | Association, credible sets and expression QTLs | Resource-specific variants, cell annotations and biosamples | Genetic enrichment and candidate-gene evidence |

RNA and ATAC in the reference atlas measure the same nuclei. Spatial RNA and protein measurements share tissue blocks/cohort. They provide complementary measurements, not independent patient replications. Possible participant overlap between GWAS resources is unresolved. The complete source and access-term inventory is supplied in `docs/sources.tsv`.

### Processing and integration

Cell QC, doublet screening, broad-type annotation and donor-aware analysis precede signature testing. Low-yield samples and non-myeloid contaminants are handled explicitly. Macrophage differential expression uses donor pseudobulk and adjusts for sorting panel. Technical preparation effects are recorded rather than interpreted as disease effects.

RNA, ATAC, concatenation, weighted nearest neighbours (WNN) and one MOFA+ configuration are compared on shared cells using biological-label and batch metrics. Cross-cohort analysis uses named gene sets and records which genes each assay measures. No all-patient, all-modality matrix is assumed where samples are unmatched.

Genetic evidence combines locus/credible-set summaries, conditional stratified LD-score regression and expression-based scDRS. Spatial tests use within-section empirical label permutations. Protein comparisons aggregate regions within donor before testing. Supervised modelling is confined to the bulk cohort and uses training-only preprocessing inside nested cross-validation.

<!-- PAGEBREAK -->

## 3. Integrated biological findings

### Cell states and genetic association

The RNA/ATAC benchmark favours concatenation and WNN overall in the tested configuration. RNA alone better preserves rare M-6 neighbourhoods, and none of the tested integrations isolates M-6 as a cluster. The author-defined M-6 state has only 65 nuclei. These observations make rare-state resolution a practical constraint on interpretation (Figure atlas, Fig. 1).

Conditional S-LDSC associates T-cell chromatin with RA heritability in the Ishigaki European and FinnGen analyses, with Tph/Treg-related states among the leading subtypes. scDRS provides complementary expression-based evidence using the same genetic traits; cross-state agreement is modest (rho 0.26, P=0.047). This reanalysis supports the relevance of published immune states without establishing a directional T-cell-to-macrophage mechanism (Fig. 2).

### Distinct scores give different answers

| Measurement | Supported result | Boundary |
|---|---|---|
| S1: 15-gene macrophage discovery signature | Highest scores in M-6; donor RNA/ATAC concordance rho 0.84 | Independent inflammation association P=0.210; bulk RA/OA P=0.327; random-set P=0.525 |
| S2: separate 10-gene IFN programme | IFNG-positive neighbourhood association in 3/4 RA sections; bulk specificity P=0.002 | Spatial label shuffles do not preserve compartments or expression-label clustering |
| Hallmark IFN-gamma protein score | RA/healthy donor comparison P=0.024 in both lining and sublining | Random-set specificity supported in lining; sublining P=0.339 |
| S1 during anti-TNF treatment | Decrease within good responders, paired P=0.014 | Good/nonresponder changes differ at P=0.131; response specificity is unproven |

S1 and S2 overlap at only APOL6 and GBP5. Thirteen of the 15 S1 genes are measurable in the atlas and anti-TNF cohort. The protein pathway measurement uses Hallmark genes available in the protein assay. The S2 and protein results therefore cannot substitute for successful validation of S1 (Fig. 3).

CXCL13-positive T cells have more B-cell neighbours in all four RA 5K sections, each at the empirical permutation resolution of about P=0.001. These section-level associations identify tissue patterns worth investigating. The global-shuffle null does not establish that a ligand caused the observed aggregates.

<!-- PAGEBREAK -->

## 4. Model benchmarks and follow-up hypotheses

### Internal model evaluation

The RA/OA task uses 174 biopsies (152 RA, 22 OA). The ACPA-positive/negative RA task uses 139 biopsies (96 positive, 43 negative). Both use five outer folds repeated five times and three-fold inner tuning. Filtering, normalisation reference, scaling, variable-gene selection and imputation are fitted within the relevant training fold.

{{CLASSIFIER_TABLE}}

Values are mean AUC ± SD across 25 correlated outer folds; SD is not a confidence interval. The genome-wide comparator is an elastic-net model using up to 2,000 training-selected variable genes.

In RA/OA, the genome-wide model performs better than the programme model (approximate corrected resampled comparison P=0.00596). The programme model's higher ACPA point estimate is not significant in that comparison (P=0.149). Prioritised genes do not outperform 100 feature-count-matched random sets in either task (empirical P=0.2475 and 0.5149). Programme-model label tests reach the 199-permutation resolution of P=0.005 using the same training and tuning procedure. These tests establish internal association conditional on the specified feature lists; they do not adjust for earlier exploratory task or programme choices. External clinical validation is unavailable (Extended Data Fig. 1).

### Candidate regulatory and genetic evidence

RNA regulon and chromatin-motif comparisons nominate candidate regulators, with RNA BH adjustment across 464 TF tests. M-6-bearing units come from only two donors, so this is not replicated donor-level regulator discovery. IRF5 colocalisation supports shared genetic association signals. Single-instrument CTLA4/IL6ST Wald ratios provide expression-associated hypotheses, with unresolved pleiotropy and shared-instrument attribution. Abatacept binds CD80/CD86; IL6ST/gp130 is distinct from the IL6R target of tocilizumab. These analyses do not estimate drug efficacy (Fig. 4).

### Transfer beyond synovial tissue

In nine paired synovial-fluid patients, the M-6-like fraction decreases in 7/9, while the paired M-6 pseudobulk score test is negative (P=0.129). In LuCA, macrophage M-6 scores are higher in tumour in 20/29 paired donors (P=0.0038), but macrophage/Tph coupling is not conventionally significant (P=0.092), and primary TCGA survival tests are negative. These are exploratory transfers of separate 30-gene programmes (Extended Data Fig. 2).

<!-- PAGEBREAK -->

## 5. Reproducibility, limitations and conclusion

### Reproducibility and delivery

The ordered pipeline was rerun in a clean analysis project with pinned environments, tool versions and frozen reference/API inputs. It completed all 43 scheduled invocations. An extracted-package check independently reproduced 112 result and 15 processed tables and regenerated all six figures without large single-cell intermediates. Raw-input and reference integrity, model-procedure invariants and figure/source consistency were checked. Detailed commands, platform requirements and evidence are in the Methods and QC appendix and archive, rather than this report's narrative.

### Limits on interpretation

The macrophage discovery signature is sensitive to a single donor, has negative independent inflammation/bulk tests and has no FDR-significant peak-gene link in this analysis. Rare-state and regulatory comparisons have limited donor support. Spatial shuffles omit compartment and label-clustering structure. Shared nuclei or tissue blocks cannot be counted as independent validation cohorts. Anti-TNF changes are descriptive and do not demonstrate responder-specific prediction. All classifiers are internally evaluated and conditional on exploratory feature choices. Single-instrument genetic estimates cannot resolve pleiotropy, drug effects or all genes sharing an instrument.

### Conclusion

The submission supplies a repeatable route from open multi-omics resources to comparable features, baseline models and a documented evidence map. It identifies IFN-associated tissue patterns and candidate follow-up genes while retaining unsuccessful signature, transfer and prediction checks. Its immediate use is reproducible research and study design; mechanistic and clinical claims require further data.

### Selected sources

[1] Gupta A, Weinand K, et al. Dynamic regulatory elements in single-cell multimodal data implicate key immune cell states enriched for autoimmune disease heritability. Nat Genet. 2023;55:2200–2210. [doi:10.1038/s41588-023-01577-7](https://doi.org/10.1038/s41588-023-01577-7).

[2] Alivernini S, et al. Distinct synovial tissue macrophage subsets regulate inflammation and remission in rheumatoid arthritis. Nat Med. 2020;26:1295–1306. [doi:10.1038/s41591-020-0939-8](https://doi.org/10.1038/s41591-020-0939-8).

[3] Wang J, et al. Synovial inflammatory pathways characterize anti-TNF-responsive rheumatoid arthritis patients. Arthritis Rheumatol. 2022;74:1916–1927. [doi:10.1002/art.42295](https://doi.org/10.1002/art.42295).

[4] Guo Y, Walsh AM, et al. J Immunol. 2017;198:4490. [doi:10.4049/jimmunol.1601988](https://doi.org/10.4049/jimmunol.1601988). Synovial bulk cohort: GSE89408.

[5] Sisó S, et al. Spatial multiomics identifies rheumatoid arthritis synovial niches and de-novo laminar histological patterns of fibroblasts subtypes. Front Immunol. 2026. [doi:10.3389/fimmu.2026.1883263](https://doi.org/10.3389/fimmu.2026.1883263). GSE347958/GSE347797.

[6] Wang et al. Proteomes. 2025;13:17. [doi:10.3390/proteomes13020017](https://doi.org/10.3390/proteomes13020017). Spatial proteomics: PXD062696. Full citations and accessions for all resources are in `docs/sources.tsv`.
