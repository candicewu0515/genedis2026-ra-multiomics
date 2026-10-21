# Revision and validation record (GeNeDis 2026)

Participant: Xia Wu, PhD student, University of Iowa · xwu76@uiowa.edu.
Scope: an exploratory multi-omics integration workflow that follows the official public deliverable list. No separate organiser template was supplied.

## Revision traceability

| Issue | Correction | Evidence |
|---|---|---|
| Protein plot and test used different sublining groups | All immune/fibrous/sublining regions are pooled within each donor before plotting and testing | Fig3g source values; sample/P assertions; submission checks |
| Hidden figure cache | Steps 07/20 export small plotting inputs, and the six figures render from packaged tables | Figures rebuilt in isolation |
| Classifier preprocessing and null mismatch | Every fitted transform stays inside the inner/outer training folds. The 199 null labels repeat the same tuning/CV. The 100 random controls match the retained F2 gene count | Protocol JSON, fold/null tables, invariance tests |
| Nominal TF selection | BH adjustment over all 464 RNA TF tests; the nominal flag is kept for transparency | Corrected q values and candidate flags; IRF5 is not a corrected RNA candidate |
| Zero positive-neighbour exposure | The healthy Xenium section has no IFNG-positive exposure, so its coefficient and P are marked non-estimable, not zero/null | Extra step-20 rerun; test_status and submission assertion |
| Extreme spatial combined P | Normal-tail extrapolation from effect z removed; section-level effects and empirical P kept | Xenium/Tph tables; figure panel 2i |
| Mixed signature claims | Clear S1/S2/Hallmark/M-6/Tph definitions and a score registry | Narrative, figure titles/legends and evidence_map.tsv |
| Causal/drug overclaims | Association/candidate wording; molecular distinctions for CTLA4-Ig and IL6ST/IL6R; two-donor and pleiotropy limits stated | README, narrative, QC report, presentation and email |
| Independence/novelty | Published states and original genetic findings credited; same-nucleus and shared-block measurements marked as such; GWAS independence not claimed | Sources and narrative |
| Knowledge-graph semantics | S2 has its own score node; spatial nodes reflect measured cells; edges describe associations or predicted interactions | Graph source, nodes/edges and semantic checks |
| Reference APIs could change | Saved copies of Open Targets 26.09 and an API/GMT/eQTL snapshot, with per-file checksums and a verified restore | Reference manifest and restore log |
| Platform/conversion dependencies | Python locks, a documented R inventory and conversion, Rosetta/MAGMA constraints, pinned tool hashes | Environments, R conversion comparison and tool checks |
| Partial S-LDSC output could silently shrink the FDR family | The step stops on any failed/missing model and needs all 195 unique finite results before BH | Completeness failure-mode tests; final step-33 entry-point check and submission checks |
| Old captions/slides/contact fields | Panel notes, the six-figure set and the slide deck rebuilt; individual details filled in | Final artifact checks |
| Submission package could change itself during checking | Quick verification runs in a temporary project and rechecks the original hashes | Quick verification of the extracted package |
| Report layout and wording | The report and slides use black and grey headings, text, links and table lines. The README, report, figure legends and slides were reworded in plain English. Numbers, tables and figures did not change | Rebuilt report and slides; number and name check of every edited text |

## Validation status

**Ordered analysis: PASS.** All 43 runs completed successfully (41 numbered steps, with steps 22/32 repeated to resolve dependencies). The analysis began 2026-10-20 13:30:36, and the final figure export ended 2026-10-20 16:17:35 (America/Chicago; 2.78 h elapsed).

The project started with empty analysis intermediates and results. Primary raw files were byte-identical local copies of the public input files, checked by hashes. The API/GMT/eQTL reference snapshots were kept fixed on purpose. Independent branches were scheduled after their data dependencies. The ordered driver reused only derived files generated inside this same clean validation project. The early step-36 branch was stopped after MAGMA to control memory. A resumed step-36 branch then completed scDRS after MOFA released memory. The ordered step 36 checked and reused these results by their exact input and output hashes. These scheduling actions are recorded in the logs, and no old analysis model was substituted. During validation, the 22 baseline annotation files were converted with the pinned LDSC parser to Parquet, a format LDSC supports natively. Exact value/dtype round-trip checks were run, and the original gzip inputs were kept. Later regressions used this storage cache. Model definitions and numeric input values did not change. The conversion evidence is in logs/annotation_parquet_cache.json. An extra step-20 rerun marks the healthy section with no IFNG-positive neighbours as non-estimable. It keeps the RA estimator and the permutation RNG sequence unchanged, and its completion is recorded separately. A further step-33 run checked the final cache-preparation entry point. It reused only the regressions freshly computed in this project and is logged separately.

Table comparison against the pre-revision baseline: {"identical": 103, "expected revision": 7, "new output": 2}. Differences in the classifier, corrected TF selection and spatial summaries are expected and documented below. The full per-table comparison is in results/reproducibility_check.tsv. New null-distribution tables are marked as new outputs.

Checks completed before the final package was saved:
- Four Seurat-converted input files were reproduced with identical decompressed contents.
- Invariant checks passed for training reference, query independence, nested holdout, matched permutation and matched random-feature count.
- The first isolated quick check reproduced 110 result tables, 15 processed tables and all six figures, without the large single-cell raw/interim objects.
- Final submission consistency checks are in results/submission_checks.json.
- Verification of the extracted archive and the final internal re-review are recorded in the release verification note, after the archive is finalised.

## Interpretation boundaries

M-6 has 65 discovery nuclei, and regulatory comparisons use two donors. Global spatial label shuffles do not keep tissue compartments or spatial label clustering. The independent inflammation and bulk tests of signature S1 are negative. The anti-TNF result is a within-group change, not a significant difference between response groups. Wald ratios are single-instrument, expression-associated estimates, not drug effects. All classifiers are internally validated, and no clinical utility is established. Modalities from a shared cohort give complementary measurements, not independent human replication.

## Revised classifier values

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

AUC SD is taken across correlated outer folds, and random-set SD across draw means. Neither is a confidence interval.

- **RA vs OA**: F1 permutation p 0.005, with the same repeated nested-CV procedure. F2 vs random empirical p 0.2475; SD across random-set mean AUCs. F1 − genome-wide: -0.079; t -3.017, p 0.00596; approximate corrected resampled test.
- **ACPA+ vs ACPA- (RA)**: F1 permutation p 0.005, with the same repeated nested-CV procedure. F2 vs random empirical p 0.5149; SD across random-set mean AUCs. F1 − genome-wide: 0.072; t 1.49, p 0.1492; approximate corrected resampled test.

All estimates are internal to GSE89408. A low label-permutation P shows association under this modelling procedure. It does not show external predictive performance or specificity of the candidate biology.

