# RA multi-omics evidence map (GeNeDis 2026)

**Individual participant:** Xia Wu · PhD student, University of Iowa · xwu76@uiowa.edu  
**Scope:** rheumatoid arthritis, with an exploratory NSCLC extension.

This is a reproducible workflow that combines open genetic, single-cell, chromatin, bulk, spatial and protein data. It produces an evidence map to guide follow-up experiments. Associations, internal classifiers and single-variant genetic estimates do not establish a causal pathway or clinical utility.

## Start here

- [Study narrative and figure legends](docs/story.md)
- [Detailed QC and model report](docs/QC_model_report.md)
- [Score definitions and the limits of each score](docs/evidence_map.tsv)
- [Source/accession list](docs/sources.tsv)
- [Revision and validation record](docs/REVISION_VALIDATION.md)
- Four main figures and two Extended Data figures: [results/figures/main](results/figures/main)
- [Presentation](presentation/GeNeDis2026_RA_results.pptx)

The CD4 genetic enrichment comes from a reanalysis of a published atlas. M-6 is a macrophage state defined by the atlas authors. These are different measurements: the discovery signature (S1, 15 genes), the predefined IFN programme (S2, 10 genes), the Hallmark IFN-gamma genes and the 30-gene M-6/Tph programmes. RNA and ATAC come from the same nuclei, and the spatial and protein data come from shared tissue blocks. These modalities complement each other. They are not independent patient replications.

## Findings and limits

1. Regulatory DNA linked to CD4 Tph/Treg states shows RA genetic enrichment in two GWAS. scDRS gives a complementary, expression-based analysis of the same genetic traits. We reuse the fine-state annotations from Gupta et al. We do not claim to have discovered these cell states or their original heritability enrichment.
2. S1 scores highest in M-6, but its association with donor inflammation in the multiome cohort is not significant. S1 fails its bulk RA/OA specificity test. S2 shows bulk and spatial associations, but these cannot replace S1 validation.
3. Hallmark IFN-gamma protein scores differ between RA and healthy donors in each compartment. Random-set tests support specificity in lining, not in sublining. The corrected plot averages all sublining regions within each donor (6 RA, 3 healthy).
4. S1 decreases after anti-TNF within good responders. This does not show response specificity, because the change does not differ significantly between good responders and nonresponders. Candidate TF evidence comes from only two M-6-bearing donors and uses gene-label permutations with BH adjustment.
5. The CTLA4/IL6ST Wald ratios are exploratory, expression-associated evidence. Abatacept binds CD80/CD86, and IL6ST is not the IL6R target of tocilizumab. Genetic effects on expression do not estimate drug efficacy.
6. The revised eight-programme classifier has an RA/OA AUC of 0.900, against 0.979 for the genome-wide baseline. The prioritised genes do not do significantly better than matched random gene sets in either task (empirical P=0.2475 and 0.5149).
7. The synovial-fluid and NSCLC analyses are exploratory extensions that use separate 30-gene programmes. All classifier performance is internal to one bulk cohort, and no clinical classifier is validated.

## Data (open, no controlled access)

| Accession | Content | Donors | Labels used |
|---|---|---|---|
| GSE243917 (GEO) | Synovial tissue, 10x single-nucleus Multiome (RNA + ATAC, same nuclei) | 11 RA + 1 OA | disease, lymphocytic-inflammation grade, lining-hyperplasia grade |
| E-MTAB-8322 (EBI SCEA) | Sorted synovial tissue macrophages, 10x 3' scRNA-seq | 26 (4 healthy, 3 UPA, 6 naive-active RA, 6 treatment-resistant RA, 7 remission RA) | disease, treatment status, sort panel |
| GSE198520 (GEO) | Synovial bulk RNA-seq before and 12 weeks after anti-TNF | 46 RA (92 samples) | EULAR response, biologic, DAS28, pathotype, histology, RIN |
| GSE89408 (GEO) | Synovial bulk RNA-seq across disease stages | 28 healthy, 10 arthralgia, 6 UA, 57 early RA, 95 established RA, 22 OA | stage, age, sex, ACPA |
| GWAS Catalog | RA and LDL (control) genome-wide significant associations, GRCh38 | n/a | trait |
| GSE347958 (GEO) | Synovium Xenium in situ (480-gene and 5K panels), FFPE | 5 RA, 3 healthy, 1 OA donor (3 blocks); 17 sections | condition |
| PXD062696 (PRIDE) | LCM spatial proteomics of the same blocks (DIA-NN matrices) | 6 RA, 3 healthy, 1 OA (24 regions) | disease, region |
| GSE347797 (GEO) | Synovium Visium (FFPE), same cohort | 6 RA, 3 healthy, 1 OA donor (3 blocks) | condition |
| GCST90132223 + eQTL Catalogue | RA GWAS summary statistics (EUR) and myeloid eQTLs (tabix) | 22,350 cases; 84–200 eQTL donors | n/a |
| Open Targets credible sets and colocalisation | Fine-mapped RA GWAS (Ishigaki 2022 EUR and EAS, FinnGen R12), LDL control; QTL colocalisation | 64 + 51 + 41 + 165 credible sets | n/a |
| S-LDSC reference + full summary statistics | 1000G EUR GRCh38 baselineLD v2.2 (Zenodo 10515792); Ishigaki 2022 EUR RA, FinnGen R12 RA and T2D (control) | 489 reference genomes | n/a |
| CollecTRI (OmniPath) | Signed TF-target regulons | 464 TFs used | n/a |
| L1000CDS2 / ChEMBL | Perturbation signatures and compound mechanisms (APIs) | n/a | n/a |
| CELLxGENE PBMC (JCI Insight 2024) | RA disease-activity blood scRNA-seq | 18 RA + 18 matched controls | disease, DAS28, CDAI |
| TCGA LUAD / LUSC (UCSC Xena) | Bulk RNA-seq + curated survival | 491 + 489 patients | OS, stage, age |
| GSE296117 (GEO) | RA synovial fluid scRNA-seq, treatment-naive and on adalimumab or tofacitinib | 9 RA (18 samples) | treatment, drug, author cell type |
| LuCA core atlas (CELLxGENE) | NSCLC and normal lung scRNA-seq atlas (19 studies), lung tissue only | 151 tumour donors; 29 with paired adjacent normal macrophages | origin, disease, study |

Full source details (repository, licence, publication, files, role) are in [`docs/sources.tsv`](docs/sources.tsv).
Raw files are not stored in the repository. `scripts/00_download.sh` downloads them and writes
`data/raw/SHA256SUMS.txt`.

## Reproduce

The tested platform is macOS arm64 (Apple M4 Max, 48 GB RAM). The exact Python environments are in env/. Other platforms need compatible environments and a suitable MAGMA binary. We do not claim identical numerical results across platforms. The tested MAGMA build is x86_64 and needs Rosetta 2 on Apple Silicon. The full downloader and tool-hash check currently works only on macOS. On Linux, the MAGMA download and platform hash must be adapted. Windows needs a POSIX environment such as WSL2. A full cross-platform run has not been certified.

~~~bash
mamba env create -f environment.yml
mamba env create -f environment_ldsc.yml
mamba activate genedis_ra
bash scripts/verify_quick.sh
bash scripts/00_download.sh
bash scripts/run_all.sh
~~~

The complete ZIP manifest covers every packaged file. The Git-checkout manifest covers the versioned files. The large reference archive is left out of Git and pinned separately by env/reference_snapshot.sha256. Quick verification needs only the versioned files, including the small plotting inputs and offline viewer assets.

Quick verification checks package hashes and reruns the small bulk analyses in a temporary directory. It then compares tables and renders all six figures from the packaged plotting inputs. It does not change the submission tables. The full analysis starts from public raw inputs (about 24 GB plus extracted references) and rebuilds the intermediate objects. Final runtime, completion status and table comparisons are in docs/REVISION_VALIDATION.md and logs/.

GSE296117 is deposited as a Seurat object. The first conversion needs R 4.5.2 with Seurat 5.5.0 and Matrix. The script scripts/install_r_dependencies.R installs these packages. The tested R session is recorded in env/R_sessionInfo.txt. The conversion exports sparse raw counts and metadata, not analysis results.

The archive includes references/reference_snapshot.tar.gz with a per-file manifest. It holds the exact API responses (Open Targets 26.09, GWAS Catalog, Ensembl, CollecTRI, L1000/ChEMBL), GMT gene sets and eQTL regions used. The download script restores and verifies it before analysis. If the snapshot archive is missing, a Git checkout downloads the unchanged reference asset from v1.4-submission. Its expected SHA-256 is pinned in env/reference_snapshot.sha256. The saved snapshot is the reference version for the submission. Live APIs may change.

Rendering the figures alone needs no raw or interim data:

~~~bash
python scripts/39_main_figures.py
~~~

The small plotting tables are outputs of steps 07/20 and are supplied under data/processed/. Every figure comes with source hashes and panel descriptions.

## Analysis structure

- 01–08: QC, annotations, macrophage discovery, RNA/ATAC integration and signature mapping.
- 09–17: genetic overlap, response/pathway tests, interaction hypotheses, motifs, target annotations and robustness.
- 18–27: blood/lung context, Xenium, proteomics, gene ranking, Visium, colocalisation and negative controls.
- 28–38: integration baselines, fine-mapping, bulk stages, graph proximity, QTL directions, S-LDSC, TF candidates, Tph neighbourhoods, scDRS, exploratory Wald ratios and revised nested-CV classifiers.
- 40–41: synovial-fluid treatment and NSCLC extensions. Steps 23/39 build the graph and figures.

Step 33 needs all 195 unique finite regression results from the fixed atlas version before it computes the BH adjustment. A failed or incomplete family stops the step. Step 33 also builds native Parquet caches of the 22 baseline annotation files with the pinned LDSC parser. It checks an exact value/dtype round trip and keeps the original gzip inputs. LDSC natively prefers this storage format. Model inputs and parameters are unchanged. Cache checksums and conversion checks are in data/interim/ldsc/annotation_parquet_cache.json (a copy for the final saved run is in logs).

Steps 22/32 run twice to resolve dependencies. Script headers state the analysis questions. This public-data reanalysis was not preregistered. The QC report separates scores fixed before transfer from post-hoc changes. Random generators use recorded seeds.

The revised classifier fits these steps only within each training fold, including inner folds: filtering, the count-normalisation reference, gene scaling, age imputation, variable-gene selection and hyperparameters. Observed and permuted F1 use the same 5-fold outer CV repeated 5 times, 3-fold tuning and the same grid. There are 199 label permutations and 100 random controls. With the plus-one correction, the smallest possible empirical permutation P is 0.005. These are exploratory P values for each task. Random controls use fixed random rankings. In each training fold, they match the retained feature count of the 30-gene prioritised list. These CV and permutation estimates are conditional on the given feature lists. They do not correct for earlier exploratory choices of tasks or programmes. The SD across correlated folds is not a confidence interval.

The package contains code, exact environments, documentation, processed tables, six composite figures (PDF/SVG/600-dpi PNG), figure source data, the knowledge graph, the presentation, run logs and a SHA-256 manifest. It leaves out primary raw inputs, interim models and external tools, but it includes the saved reference/API inputs. The revised submission is release [v1.4-submission](https://github.com/candicewu0515/genedis2026-ra-multiomics/releases/tag/v1.4-submission). The older v1.0 to v1.3 releases are kept as a record and have been superseded.

## Rebuilding the presentation and readable report

The final PDFs/PPTX are included, so you can view them without a build environment. To rebuild them after the analysis and validation record are complete, use a separate Python environment with `pip install -r env/presentation_requirements.txt`. Then run `python scripts/build_presentation.py` and `python scripts/build_report.py`. The report PDF is built directly with ReportLab. Slide PDF conversion needs LibreOffice (`soffice` or `libreoffice` on PATH). These document-build dependencies are separate from the pinned analysis environments. Arial is used when available, and the report falls back to DejaVu Sans. LibreOffice converts the slides to PDF, and the editable PPTX remains the source presentation.

After rebuilding the PPTX, convert it with `soffice --headless --convert-to pdf --outdir presentation presentation/GeNeDis2026_RA_results.pptx`. Then run `python scripts/check_artifacts.py`. It checks the 12-slide bounds, embedded figure hashes, PDF page counts and contact text, and the 600-dpi figure exports.
