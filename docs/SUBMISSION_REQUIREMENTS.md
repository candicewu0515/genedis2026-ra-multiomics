# Official submission requirements and package map

Basis: the [official GeNeDis 2026 hackathon page](https://2026.genedis.org/hackathon/), checked 20 October 2026. No additional organiser template was supplied. This package follows the public requirements. It does not claim to match a template that was not available.

Participant: **Xia Wu**, PhD student, University of Iowa. Contact: **xwu76@uiowa.edu**. Individual entry. Scope: **RA**, with an exploratory NSCLC extension.

| Public deliverable | Included material |
|---|---|
| Short README | README.md |
| Source/accession list and documentation | docs/sources.tsv; source-selection and dataset sections in README.md and docs/QC_model_report.md; raw-input and reference manifests |
| Code or pipeline | scripts/; environment.yml; environment_ldsc.yml; exact environment inventories in env/ |
| Processed feature tables | data/processed/; cohort labels and source identifiers kept |
| Key outputs or representative sample | results/tables/; results/figures/main/; results/knowledge_graph/; presentation/ |
| QC/model report | docs/QC_model_report.md; docs/REVISION_VALIDATION.md; results/submission_checks.json; logs/ |
| Clearly stated limitations | README.md; docs/QC_model_report.md; docs/story.md; presentation |

The release also contains a readable report, a presentation, six composite figures with captions and an integrity manifest. The figures use the final regenerated results. The archive leaves out large primary raw inputs and intermediate models. The download and conversion steps rebuild these from documented public sources. Saved API/GMT/eQTL reference inputs are included.

## Delivery

The public page asks for submission by email to **aris.vrahatis@ionio.gr** by **31 October 2026**, as a compressed archive or a link to a repository or storage location. Registration is a separate form linked from the official page.

The official page gives inconsistent hackathon-day dates in different sections. The submission files do not state an event-day date. The submission deadline is the same throughout the page.
