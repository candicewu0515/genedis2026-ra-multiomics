# A reproducible multi-omics workflow for rheumatoid arthritis

**Technical report | GeNeDis 2026 Hackathon**

Xia Wu | PhD student, University of Iowa | xwu76@uiowa.edu

## 1. Study overview

### Abstract

Public rheumatoid arthritis (RA) datasets differ in tissue sampling, assay coverage and clinical labels. I analysed a paired RNA/ATAC atlas, sorted synovial macrophages, bulk RNA cohorts, spatial transcriptomics, region-resolved proteomics and genetic association data to assess how cell-state findings transfer between these settings. The atlas contained 31,520 nuclei from 12 donors. T-cell chromatin was enriched for RA heritability, while spatial and protein analyses supported an IFN-associated tissue signal. A 15-gene macrophage signature, derived from an active/remission contrast, scored highest in the published M-6 state but failed independent inflammation and bulk RA/OA tests. In nested cross-validation, RA/OA AUC was 0.900 for programme scores and 0.979 for a genome-wide comparator; prioritised genes did not outperform matched random sets. Shared nuclei and tissue blocks limit the independence of some comparisons. The accompanying code and processed tables allow the analyses, including unsuccessful transfer tests, to be reproduced.

### Aim

The main question was whether genetic association, macrophage expression and tissue neighbourhood measurements supported the same RA cell-state patterns. Synovial tissue formed the primary analysis. Synovial-fluid and non-small cell lung cancer (NSCLC) cohorts were used to examine transfer of specified atlas-state programmes.

### Analysis

The RNA/ATAC comparison used published cell-state annotations from Gupta et al. [1]. The macrophage discovery contrast used the active/remission cohort of Alivernini et al. [2]. Donor pseudobulk, measured gene coverage and matched random sets were used to assess transfer between assays. Classification was restricted to the bulk cohort, with simple covariate and genome-wide baselines. Source accessions, processing steps, score definitions and figure data are supplied with the code.

{{WORKFLOW}}

### Supporting files

The **Methods and QC appendix** describes preprocessing, sample units and model fitting. The **Figure atlas** contains four main and two Extended Data figures with panel captions. The archive contains code, processed features, result tables and source records. A slide deck is also available.

<!-- PAGEBREAK -->

## 2. Data and methods

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

Cell filtering, doublet screening and broad-type annotation preceded signature testing. Low-yield samples and non-myeloid cells were excluded according to the criteria in the Methods appendix. Macrophage differential expression used donor pseudobulk with sorting panel as a covariate. Preparation batch was included in relevant comparisons.

RNA, ATAC, concatenation, weighted nearest neighbours (WNN) and one MOFA+ configuration were compared on shared cells using cell-label and batch metrics. Cross-cohort scores used predefined gene lists with assay-specific coverage. Unmatched cohorts were analysed separately.

Genetic analyses included locus and credible-set summaries, conditional stratified LD-score regression (S-LDSC) and expression-based scDRS. Spatial tests used label permutations within each section. Protein regions were averaged within donor before testing. Bulk classification used nested cross-validation, with preprocessing fitted on training data.

<!-- PAGEBREAK -->

## 3. Cell-state and tissue associations

### Cell states and genetic association

Concatenation and WNN ranked highest in the RNA/ATAC benchmark. RNA alone better preserved M-6 neighbourhoods. This state comprised 65 nuclei, and none of the tested integrations recovered it as a separate cluster (Figure atlas, Fig. 1).

Conditional S-LDSC associated T-cell chromatin with RA heritability in the Ishigaki European and FinnGen analyses; Tph/Treg-related states were among the leading subtypes. Agreement with expression-based scDRS across states was modest (rho 0.26, P=0.047). Both analyses used the same GWAS traits. They do not establish a directional interaction between T cells and macrophages (Fig. 2).

### Signature transfer and tissue associations

| Measurement | Result | Limitation |
|---|---|---|
| S1: 15-gene macrophage discovery signature | Highest scores in M-6; donor RNA/ATAC concordance rho 0.84 | Independent inflammation association P=0.210; bulk RA/OA P=0.327; random-set P=0.525 |
| S2: separate 10-gene IFN programme | IFNG-positive neighbourhood association in 3/4 RA sections; bulk specificity P=0.002 | Spatial label shuffles do not preserve compartments or expression-label clustering |
| Hallmark IFN-gamma protein score | RA/healthy donor comparison P=0.024 in both lining and sublining | Random-set specificity supported in lining; sublining P=0.339 |
| S1 during anti-TNF treatment | Decrease within good responders, paired P=0.014 | Between-group change was not significant (P=0.131) |

S1 and S2 share APOL6 and GBP5. Thirteen of the 15 S1 genes were measured in the atlas and anti-TNF cohort. The protein score used measured Hallmark IFN-gamma genes. These are three different scores: the S2 and protein associations do not validate the S1 discovery signature (Fig. 3).

CXCL13-positive T cells had more B-cell neighbours in all four RA 5K sections (each empirical P about 0.001). Global label shuffles did not preserve tissue compartments or expression-label clustering, so these tests cannot determine whether CXCL13 caused the aggregates.

<!-- PAGEBREAK -->

## 4. Classification and follow-up analyses

### Internal model evaluation

The RA/OA task uses 174 biopsies (152 RA, 22 OA). The ACPA-positive/negative RA task uses 139 biopsies (96 positive, 43 negative). Both use five outer folds repeated five times and three-fold inner tuning. Filtering, normalisation reference, scaling, variable-gene selection and imputation are fitted within the relevant training fold.

{{CLASSIFIER_TABLE}}

Values are mean AUC ± SD across 25 correlated outer folds; SD is not a confidence interval. The genome-wide comparator is an elastic-net model using up to 2,000 training-selected variable genes.

In RA/OA, the genome-wide model performs better than the programme model (approximate corrected resampled comparison P=0.00596). The programme model's higher ACPA point estimate is not significant in that comparison (P=0.149). Prioritised genes do not outperform 100 feature-count-matched random sets in either task (empirical P=0.2475 and 0.5149). Programme-model label tests reach the 199-permutation resolution of P=0.005 using the same training and tuning procedure. These tests establish internal association conditional on the specified feature lists; they do not adjust for earlier exploratory task or programme choices. External clinical validation is unavailable (Extended Data Fig. 1).

### Regulatory and genetic candidates

RNA regulon tests used BH adjustment across 464 transcription factors, followed by a chromatin-motif criterion. M-6-bearing units came from two donors, limiting donor-level replication. IRF5 colocalisation supported shared genetic association signals. CTLA4/IL6ST Wald ratios used single instruments; pleiotropy and shared-instrument attribution remain unresolved. Abatacept binds CD80/CD86, while IL6ST/gp130 differs from the IL6R target of tocilizumab. The genetic estimates therefore cannot be interpreted as drug effects (Fig. 4).

### Transfer beyond synovial tissue

In nine paired synovial-fluid patients, the M-6-like fraction decreases in 7/9, while the paired M-6 pseudobulk score test is negative (P=0.129). In LuCA, macrophage M-6 scores are higher in tumour in 20/29 paired donors (P=0.0038), but macrophage/Tph coupling is not conventionally significant (P=0.092), and primary TCGA survival tests are negative. These are exploratory transfers of separate 30-gene programmes (Extended Data Fig. 2).

<!-- PAGEBREAK -->

## 5. Reproducibility and limitations

### Reproduction

The pipeline completed 43 ordered invocations in a clean analysis directory using pinned environments, tool versions and reference inputs. A separate check of the extracted archive reproduced 112 result tables and 15 processed tables and regenerated all six figures without large single-cell intermediates. Input checksums, model-preprocessing checks and figure/source comparisons are recorded in the archive. Commands and platform requirements are given in the Methods and QC appendix.

### Limits on interpretation

The macrophage discovery signature is sensitive to a single donor, has negative independent inflammation/bulk tests and has no FDR-significant peak-gene link in this analysis. Rare-state and regulatory comparisons have limited donor support. Spatial shuffles omit compartment and label-clustering structure. Shared nuclei or tissue blocks cannot be counted as independent validation cohorts. Anti-TNF changes are descriptive and do not demonstrate responder-specific prediction. All classifiers are internally evaluated and conditional on exploratory feature choices. Single-instrument genetic estimates cannot resolve pleiotropy, drug effects or all genes sharing an instrument.

### Conclusion

T-cell genetic enrichment and IFN-associated tissue patterns were observed, but support for the macrophage discovery signature was limited outside its discovery setting. Prioritised genes also failed the matched-random classifier comparisons. These results favour follow-up of specific tissue associations over claims of a general RA biomarker. The supplied tables and code provide a basis for testing those associations in independent donors.

### Selected sources

[1] Gupta A, Weinand K, et al. Dynamic regulatory elements in single-cell multimodal data implicate key immune cell states enriched for autoimmune disease heritability. Nat Genet. 2023;55:2200–2210. [doi:10.1038/s41588-023-01577-7](https://doi.org/10.1038/s41588-023-01577-7).

[2] Alivernini S, et al. Distinct synovial tissue macrophage subsets regulate inflammation and remission in rheumatoid arthritis. Nat Med. 2020;26:1295–1306. [doi:10.1038/s41591-020-0939-8](https://doi.org/10.1038/s41591-020-0939-8).

[3] Wang J, et al. Synovial inflammatory pathways characterize anti-TNF-responsive rheumatoid arthritis patients. Arthritis Rheumatol. 2022;74:1916–1927. [doi:10.1002/art.42295](https://doi.org/10.1002/art.42295).

[4] Guo Y, Walsh AM, et al. J Immunol. 2017;198:4490. [doi:10.4049/jimmunol.1601988](https://doi.org/10.4049/jimmunol.1601988). Synovial bulk cohort: GSE89408.

[5] Sisó S, et al. Spatial multiomics identifies rheumatoid arthritis synovial niches and de-novo laminar histological patterns of fibroblasts subtypes. Front Immunol. 2026. [doi:10.3389/fimmu.2026.1883263](https://doi.org/10.3389/fimmu.2026.1883263). GSE347958/GSE347797.

[6] Wang et al. Proteomes. 2025;13:17. [doi:10.3390/proteomes13020017](https://doi.org/10.3390/proteomes13020017). Spatial proteomics: PXD062696. Full citations and accessions for all resources are in `docs/sources.tsv`.
