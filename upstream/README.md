# Upstream input preparation

Scripts used to extract the cohort from MIMIC-IV v3.1 and build the private input datasets read by `../revision_2026_11/` (see that README for the list). They are the October 2026 versions with the vasopressin unit correction (vasopressin units/hour divided by 60 before the 2.5 norepinephrine-equivalent multiplier). Comments not in English were translated or removed, and two output folder names were renamed (`table1_output/results`); code is otherwise unchanged.

| Folder | Role |
|---|---|
| `analysis/` | cohort extraction (`strictK3_primary_extract_v2.R`), baseline covariates (`strictK3_baseline_covariates_v2.R`), hourly covariates (`03a_timevarying_covariates.R`). `03_ccw_ipcw.R` and `strictK3_run_timevarying_ccw*.R` contain the **superseded** weighting implementation; they are kept only because their saved clone file is read by the input checks in `step1_sofa_boundary.R`. |
| `anchor/`, `rescue/` | alternative eligibility anchors and urgent-indication exclusion variants (branch cohorts). Their estimation steps are superseded by `revision_2026_11/analysis/round4_sensitivity.R`. |
| `manuscript/` | hourly covariates for hours 0-47 (48-hour grace period). |
| `secondary/` | day-28 renal recovery and RRT-free days dataset (run13). |
| `table1/` | Table 1 / Supplemental Table S14 queries and generator. |
