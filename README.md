# RA multi-omics workflow | GeNeDis 2026

**Xia Wu**, PhD student, University of Iowa · **xwu76@uiowa.edu** · Individual entry

A reproducible workflow linking open genetic, RNA/chromatin, bulk, spatial and protein resources for rheumatoid arthritis (RA). It produces documented cohort features, integration/model benchmarks and a source-linked evidence map. RA synovial tissue is the primary scope; synovial-fluid and NSCLC analyses are exploratory extensions.

## Read first

1. [Technical report](presentation/GeNeDis2026_submission_report.pdf) — five pages covering the question, contribution, data/methods, main findings, model comparisons and limitations. [Editable source](docs/competition_report.md).
2. [Methods and QC appendix](presentation/GeNeDis2026_methods_QC.pdf) — five pages covering preprocessing, sample units, score definitions, nested CV and reproducibility. [Editable source](docs/methods_and_QC.md).
3. [Figure atlas](presentation/GeNeDis2026_figure_atlas.pdf) — four main and two Extended Data figures with panel captions. Individual vector/600-dpi figures and source hashes: [results/figures/main](results/figures/main).

The [PPTX](presentation/GeNeDis2026_RA_results.pptx) and [slide PDF](presentation/GeNeDis2026_RA_results.pdf) are optional supporting material. The official public deliverable list does not require an initial PPT submission; shortlisted entrants present their work.

## What the workflow shows

Published atlas states provide the reference for paired RNA/ATAC integration and genetic enrichment. S1, the 15-gene macrophage discovery signature, maps to M-6 but has negative independent inflammation/bulk tests. The separate S2 IFN programme and Hallmark protein measurements support tissue associations within their assay-specific limits. Prioritised-gene classifiers do not significantly outperform matched random gene sets; all model evaluation is internal.

RNA/ATAC from the same nuclei and spatial/protein measurements from shared blocks are complementary measurements, not independent patient replications. Regulatory comparisons use two M-6-bearing donors. Spatial associations and single-instrument Wald ratios do not establish a causal mechanism or drug efficacy. These limitations and failed transfer checks remain in the report.

## Data, code and outputs

| Deliverable | Location |
|---|---|
| Source/accession, access terms and provenance | [docs/sources.tsv](docs/sources.tsv) |
| Score definitions and measured coverage | [docs/evidence_map.tsv](docs/evidence_map.tsv) |
| Pipeline and environments | scripts/; environment.yml; environment_ldsc.yml; env/ |
| Processed cohort/feature tables | data/processed/ |
| Result tables and statistical protocols | results/tables/ |
| Figures and knowledge graph | results/figures/main/; results/knowledge_graph/ |
| Current QC/model report and limits | Technical report and Methods/QC appendix linked above |
| Detailed validation and analysis history | docs/REVISION_VALIDATION.md; docs/QC_model_report.md; logs/ |

The detailed historical records remain available for audit. They are separate from the competition reading report. The workflow uses secondary open data; large primary raw inputs are not committed to Git.

## Reproduce

The full scientific workflow was validated on macOS arm64 (Apple M4 Max, 48 GB RAM). Python environments, the R/Seurat session and tool hashes are recorded in env/. MAGMA uses the tested x86_64 macOS build through Rosetta 2. Linux requires adapting the MAGMA download/binary; Windows requires a POSIX environment such as WSL2. Cross-platform numerical identity has not been certified.

~~~bash
mamba env create -f environment.yml
mamba env create -f environment_ldsc.yml
mamba activate genedis_ra
bash scripts/verify_quick.sh
~~~

Quick verification checks the frozen manifest, recreates small analyses in an isolated directory, compares tables and renders all six figures. It leaves the submission unchanged. Full reconstruction uses:

~~~bash
bash scripts/00_download.sh
bash scripts/run_all.sh
~~~

The complete ZIP includes the saved reference archive (Open Targets 26.09 and exact API/GMT/eQTL inputs), with archive/per-file checksums. Git checkouts restore the identical pinned asset from retained **v1.4-submission**. Primary raw files and intermediate models are omitted from the archive; the downloader and documented conversion steps reconstruct them. Full primary inputs require about 24 GB plus reference/extraction space. GSE296117 conversion requires the recorded R/Seurat environment.

The scientific validation completed 43 ordered invocations. Its extracted-package check reproduced 112 result and 15 processed tables and regenerated all six figures. Report restructuring does not rerun or alter those calculations. Version-specific integrity and public-download checks accompany each release.

## Documents and submission

Current competition package: [v1.5-submission](https://github.com/candicewu0515/genedis2026-ra-multiomics/releases/tag/v1.5-submission). The earlier v1.4 package is retained as history and as the source of the unchanged reference asset.

The [official requirements](https://2026.genedis.org/hackathon/) request a README, source list, code/pipeline, processed features, key outputs, QC/model report and limitations. This package follows that public list; no additional organiser template was provided by the entrant. Registration is separate and is not asserted completed. Submission email has not been sent.

To rebuild the report/appendix after editing their Markdown sources, install `env/presentation_requirements.txt` in a document-build environment and run `python scripts/build_report.py`. The separate figure atlas is rebuilt from the accepted figure sources if absent. `python scripts/check_artifacts.py` checks the PDFs, slide bounds and embedded figure hashes. The tested presentation builder and LibreOffice PDF conversion remain available in the source scripts.
