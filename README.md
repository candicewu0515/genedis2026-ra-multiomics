# RA multi-omics analysis | GeNeDis 2026

Xia Wu · PhD student, University of Iowa · xwu76@uiowa.edu

This repository contains analyses of RA genetic association, synovial cell states and tissue inflammation. It combines a published paired RNA/ATAC atlas with sorted macrophage RNA, bulk RNA, spatial transcriptomics and region-resolved proteomics. Synovial-fluid and NSCLC cohorts are included as exploratory comparisons.

## Reports

- [Technical report](presentation/GeNeDis2026_submission_report.pdf) — study question, datasets, methods and results (5 pages).
- [Methods and QC](presentation/GeNeDis2026_methods_QC.pdf) — sample handling, score definitions, nested cross-validation and reproduction (5 pages).
- [Figure atlas](presentation/GeNeDis2026_figure_atlas.pdf) — four main and two Extended Data figures with captions (6 pages).
- [Slides](presentation/GeNeDis2026_RA_results.pptx) / [slide PDF](presentation/GeNeDis2026_RA_results.pdf) — optional presentation.

The report sources are in `docs/competition_report.md` and `docs/methods_and_QC.md`.

## Results

T-cell chromatin was enriched for RA heritability, and spatial/protein measurements supported IFN-associated tissue patterns. S1, the 15-gene macrophage discovery signature, scored highest in the published M-6 state but failed independent inflammation and bulk RA/OA tests. The separate S2 IFN programme and Hallmark protein score use different gene sets. Prioritised-gene classifiers did not outperform matched random sets; model evaluation was internal.

Paired RNA/ATAC measures the same nuclei, and spatial/protein measurements share tissue blocks. These comparisons are not independent patient replications. Regulatory analyses include only two M-6-bearing donors. Spatial associations and single-instrument genetic estimates do not establish mechanisms or drug effects.

## Repository contents

| Files | Location |
|---|---|
| Accessions, access terms and provenance | docs/sources.tsv |
| Gene sets and measured coverage | docs/evidence_map.tsv |
| Analysis scripts and environments | scripts/; environment.yml; environment_ldsc.yml; env/ |
| Processed features | data/processed/ |
| Statistical results and model protocol | results/tables/ |
| Figures, panel data and hashes | results/figures/main/ |
| Knowledge graph | results/knowledge_graph/ |
| Analysis history and run records | docs/QC_model_report.md; docs/REVISION_VALIDATION.md; logs/ |

## Installation and reproduction

The full workflow was tested on macOS arm64 (Apple M4 Max, 48 GB RAM). Python versions, the R/Seurat session and tool hashes are recorded in `env/`. MAGMA uses the x86_64 macOS build through Rosetta 2. Linux requires an appropriate MAGMA binary and download adaptation; Windows requires a POSIX environment such as WSL2. Cross-platform numerical identity has not been tested.

~~~bash
mamba env create -f environment.yml
mamba env create -f environment_ldsc.yml
mamba activate genedis_ra
bash scripts/verify_quick.sh
~~~

The quick check verifies the manifest, downloads four small raw inputs, reruns three analyses in an isolated directory, compares tables and renders six figures. It does not modify the submission. Full reconstruction uses:

~~~bash
bash scripts/00_download.sh
bash scripts/run_all.sh
~~~

The full run completed 43 ordered invocations. Extracted-package verification reproduced 112 result tables and 15 processed tables and rendered all six figures. Primary inputs require about 24 GB, plus reference/extraction space. The GSE296117 conversion requires the recorded R/Seurat environment.

The complete ZIP includes frozen API/GMT/eQTL references, including Open Targets 26.09. Git checkouts restore the same checksum-pinned archive from retained v1.4-submission. Large primary raw files and intermediate models are omitted; the downloader and conversion steps reconstruct them. Archive and per-file checksums are supplied.

## Submission

Download: [v1.6-submission](https://github.com/candicewu0515/genedis2026-ra-multiomics/releases/tag/v1.6-submission).

The [official deliverables](https://2026.genedis.org/hackathon/) are a README, source list, pipeline, processed features, key outputs, QC/model report and limitations. These files are included. PPT is optional for the initial submission. Registration is separate; this repository does not record registration or email delivery.

To rebuild the report, install `env/presentation_requirements.txt` and run `python scripts/build_report.py`. Run `python scripts/build_presentation.py` to rebuild the deck, then convert it to PDF with LibreOffice. `python scripts/check_artifacts.py` checks page counts, slide bounds and embedded figure hashes.
