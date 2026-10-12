# PIC sodium and age association study: analysis code

This repository contains the analysis and plotting implementations used for a retrospective study of early sodium measurements and subsequent in-hospital mortality in the Pediatric Intensive Care (PIC) database. The accepted analysis cohort contains 3,459 children and 306 subsequent in-hospital deaths at a 72-hour landmark.

The repository shares code and fixed method parameters. It contains no patient-level records, database extracts, imputed datasets, fitted-model objects, individual prediction scores, manuscript files or figure/data outputs. Users must obtain PIC access independently and comply with the applicable data-use agreement.

## Statistical implementation

* `data_completion_all_variables_20261009_corrected/code/run_clean_MI.R`: the actual accepted clinical imputation implementation (50 imputations, 30 iterations, five PMM donors; QR representation of the fixed auxiliary design). `preserve_export_literals.py` preserves observed CSV literals without changing originally missing cells.
* `analysis_daily_mean_20261009_v1/code/outcome_engine.R`: primary two-segment sodium-change model, continuous-age interactions, daily-level association models and CBC/measurement-frequency supporting models; Rubin pooling and fixed joint tests.
* `prediction_day1_day3_20261009_v1/code/prediction_engine.R`: exploratory day-1/day-3 clinical models, repeated five-fold validation and train-only clinical imputation. Prediction uses 10 imputations within each training fold; it does not use the 50 full-cohort association imputations as cross-validation input.
* `km_full_followup/draw_km.R`: full hospital follow-up KM estimation and log-rank calculations. `sodium_day3_CIF_preview_20261009_v1/code/render.R` supplies the supporting Aalen–Johansen implementation with live discharge as a competing event.

## Current manuscript figure mapping

Historical source filenames are retained to preserve provenance. Their numeric suffixes do not define the current manuscript numbering.

| Current display | Code source |
|---|---|
| Figure 1: cohort flow | `figure_annotation_revision_20261011_v17/code/figure_worker/Figure1_revised.R` |
| Figure 2: observed daily changes and adjusted dynamic ORs | `side_figure3_dynamic_composite_20261010_v1/code/Figure3_dynamic_composite.R`; bounded title removal in `evidence_reporting_revision_20261010_v9/figure_title_cleanup/cleanup_vector_text.py` |
| Figure 3: KM | `figure_annotation_revision_20261011_v17/code/figure_worker/Figure4_revised.R` |
| eFigure 1: competing-risk estimates | same directory, `Supplementary_Figure1_revised.R` |
| eFigure 2: daily-level RCS associations | same directory, `Figure3_revised.R` |
| eFigure 3: day-3 status associations | `evidence_reporting_revision_20261010_v9/figure_title_cleanup/source_code/Supplementary_Figure2_original.R` |
| eFigure 4: low/high sodium histories | same `source_code` directory, `Supplementary_Figure3_original.R` (history panel g1) |
| eFigure 5: day-1/day-3 ROC | `figure_annotation_revision_20261011_v17/code/figure_worker/Figure5_revised.R` |
| eFigure 6: calibration | same `source_code` directory, `Supplementary_Figure3_original.R` (calibration panel g3) |

## Reproduction boundary

The original scripts intentionally require authorized source bindings, independent-release metadata and input hashes. Those internal records and the study data are not distributed. Therefore this code archive is not a data-free, one-command replication package. The release and no-overwrite gates remain in the code. Sanitized path changes require a new local configuration/review binding before scientific execution; original acceptance hashes must not be fabricated or reused for changed code.

See `docs/RUNNING.md`, `docs/INPUT_SCHEMA.md`, `docs/DEPENDENCIES.md` and `docs/FIGURE_LINEAGE.md`. Original study data, accepted results and submitted materials were not altered when this repository was prepared. Syntax parsing and privacy checks do not establish fresh numerical or visual reproduction.
