# Revised analysis (November 2026)

Corrected clone-censor-weighting (CCW) analysis for the revised manuscript. Scripts were run by the authors on a local MIMIC-IV v3.1 PostgreSQL installation with R 4.5 (packages: dplyr, tidyr, lubridate, survival, splines, parallel, DBI, RPostgres, zoo, readr, ggplot2, ggsurvfit, patchwork, forestploter, ragg, bit64) and Python 3 (standard library; reportlab for Figure 1). They are provided for transparency and reproduction with credentialed data access; they are not a one-command pipeline.

## Specification (primary analysis)

* Eligibility at the first physiologic KDIGO stage 3 (RRT-based staging excluded); t0 = max(ICU admission, KDIGO stage 3 time); SOFA at t0 from the latest 24-hour window ending at or before t0; all baseline values measured at or before t0.
* Strategies: CRRT start in [0, 24) hours after t0 versus no CRRT within 24 hours; each patient cloned into both strategies.
* Deferred strategy: hourly pooled logistic initiation model among clones alive, in the ICU and untreated at the start of each hour, conditional on baseline covariates (L0), covariates from the preceding hour (Lk) and a natural spline in hour; stabilising numerator depends on hour only; cumulative weight up to the start of each interval applied in counting-process format.
* Early strategy: weight 1 during the grace period; after 24 hours, initiators weighted by P(initiation by 24 h) / P(initiation by 24 h | L0) among patients alive at 24 hours.
* Cumulative stabilised weights truncated at the strategy-specific 99th percentile, as prespecified in the registered protocol (OSF 10.17605/OSF.IO/ZHP7C); untruncated weights are a sensitivity analysis.
* IPCW-weighted Cox model (counting process, strategy only), weighted Kaplan-Meier risks, risk difference and 28-day RMST; bias-corrected and accelerated (BCa) 95% CIs from 1,000 patient-level bootstrap resamples refitting both weight models, acceleration from a 200-group jackknife.

## Run order

| Step | Script | Produces |
|---|---|---|
| 1 | `analysis/step1_sofa_boundary.R` | strict SOFA at t0, early definition [0, 24) -> `private_data/cohort_v3_step1.rds` |
| 2 | `analysis/00_add_death_dates.R` | dates of death, removal of 2 inconsistent records -> `private_data/cohort_v3_step3.rds` |
| 3 | `analysis/round1_info_timing.R` | baseline values at or before t0, strictly lagged hourly SOFA/KDIGO -> `cohort_round1.rds`, `tv_hourly_round1.rds` |
| 4 | `analysis/round3_final.R` | untruncated-weight inference (BCa), used for the "no weight truncation" sensitivity row |
| 5 | `analysis/round7_truncated_primary.R` | **primary analysis** (99th-percentile truncation): runs `round3_final.R`, `round4_sensitivity.R`, `fig3_rows_2026-11.R` and `round6_results_extras.R` with `options(cp_trunc = TRUE)` |
| 6 | `analysis/round5_extra.R` | flow-diagram counts, mortality-ascertainment counts |
| 7 | `figures_tables/run_table1_adapter.R`, `run_table1_backfill.R`, `run_S14_rich.R`, `run_S5_fill.R` | Table 1, Supplemental Tables S4, S5, S8, S10, S14 (aggregate CSVs) |
| 8 | `figures_tables/bayes_evalue_qba_2026-11.py` | Bayesian reanalysis (7 priors), E-values, quantitative bias analysis (Supplemental Tables S12-S13) |
| 9 | `figures_tables/figure1_flow_2026-11_reportlab.py`, `figure2_2026-11.R`, `figure3_forest_2026-11.R`, `figureS2S3S4_2026-11.R`, `figureS6_bayes_2026-11.R` | Figures 1-3 and Supplemental Figures S2-S4, S6 |

`analysis/cp_core.R` holds the shared functions (`prep_pop`, `build`, `est`). `figures_tables/fig3_forest_rows.csv` contains the aggregate estimates plotted in Figure 3. `figures_tables/dependencies/` contains the unchanged layout scripts reused for the supplemental figures and the read-only Table 1 query adapter.

## Inputs (private, not distributed)

Built by the scripts in `../upstream/` from MIMIC-IV v3.1: `ch24_q4_corrected_inputs/{baseline,tv_hourly,tv_hour23,ccw_final}.rds`; the alternative-eligibility and urgent-indication branches under `ch24_q4_corrected_other_weights/<branch>/NEW/data/`; `ch24_q4_corrected_tv48/tv_hourly_0_47.rds` (48-hour grace period); `ch24_q4_corrected_secondary/NEW/run13/ANALYSIS_DF_AUDIT.rds` (renal recovery and RRT-free days); the cohort funnel CSV.

## Environment variables

| Variable | Meaning |
|---|---|
| `CRRT_REPO_DIR` | local clone of this repository |
| `CRRT_SCRIPT_DIR` | directory containing `cp_core.R` and the round scripts (`revision_2026_11/analysis`) |
| `CRRT_REVISION_DIR` | working directory with `private_data/` (patient-level, local only) and `out/` (aggregate reports) |
| `CRRT_INPUT_DIR` | directory holding the private input datasets listed above |
| `CRRT_OUTPUT_DIR` | output directory with `tables/` and `figures/` |
| `CRRT_FUNNEL_CSV` | cohort funnel written by `upstream/analysis/strictK3_primary_extract_v2.R` |
| `CRRT_FIG2_TEMPLATE_DIR`, `CRRT_FIG2_RLIB` | editable ggsurvfit objects and pinned library used for the submitted Figure 2 layout (not distributed) |
| `CRRT_TABLE1_REFERENCE_CSV` | earlier Table 1 output used only by the adapter's optional old-replay comparison |
| `PGHOST`, `PGPORT`, `PGDATABASE`, `PGUSER`, `PGPASSWORD` | local MIMIC-IV PostgreSQL connection |

## Changes from the October 2026 submission

See Supplemental Table S15. In brief: time-varying instead of time-fixed weights with at-risk rows ending at death or ICU discharge; early-strategy weights no longer condition on post-initiation covariates; stabilising numerators depend on time only; protocol-specified strategy-specific truncation; baseline windows ending at t0 and strictly lagged hourly covariates; APS III removed from L0 (retained as sensitivity); two records with inconsistent timestamps removed (N = 2,461; primary population 1,469); BCa bootstrap CIs refitting the weight models; renal recovery as weighted day-28 proportions. Analyses that depended on the superseded implementation (eICU feasibility, covariate-stratified analyses, SOFA cut-point scan, SOFA-2 re-ascertainment, pre-exposure falsification endpoints, multiple-imputation diagnostics, 12:00 rule, tipping-point contour) were not repeated.
