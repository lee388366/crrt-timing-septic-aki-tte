# CRRT timing TTE code dependency and revision candidate

LOCAL REVIEW CANDIDATE — NOT A RELEASE. This candidate completes source dependencies and records the NEE revision plan. It has not been executed as an end-to-end portable package and has not been uploaded or pushed. The existing submitted artifacts and public repository are not changed by this build.

`legacy_snapshot/` contains the previously prepared, hash-tracked code snapshot, including its historical estimates and Q5 implementation. It is retained for provenance, not designated as the current NEE-corrected production workflow. `current_sources/` adds the original missing core dependencies, validated current branch adapters and current summary/figure inputs. `PATCH_LOG.json` records source IDs, source/payload hashes and changes. `PATCHES/` contains credential-neutral, path-neutral minimal formula/input diffs. `MANIFEST.tsv` includes original source and packaged hashes; source IDs are logical environment-root paths rather than personal filesystem paths.

## Data access and confidentiality

No patient-level data, DOD lookup, clinical IDs, hourly caches, model RDS, database dumps or private logs are distributed. MIMIC-IV v3.1 and eICU-CRD require PhysioNet credentialing and the applicable data use agreement. Generated private inputs remain outside this package. The package contains code and editable diagram sources; no data access or database extraction was performed to build it.

Use `CONFIG.env.example` as a neutral local template. Do not add passwords, tokens or private data to this directory or a release. R scripts read database passwords from the existing `PGPASSWORD` adaptation. Other database connection fields retain their original local defaults; exporting PGHOST/PGPORT/PGUSER alone does not override explicit original fields. A fully portable database-configuration adapter remains a stated task, not a tested feature.

## What has changed scientifically

The verified correction is the vasopressin term divided by 60 before applying the original 2.5 coefficient, with all other terms, ordinary segment averaging, windows, LOCF/coalescing, cohorts and estimation choices retained. The supplied hourly core and 48h builder have this one dose-term patch. The v2 wrapper already contained its unit patch; the patched core no longer matches its old search string, so the term is not divided twice. Keep branch-specific cohorts and weights: the primary, RRT-inclusive, anchor and rescue branches cannot share weights by joining primary patients.

Original run15 era-stratified analysis uses CH01 `corrected_follow_days/corrected_event`; run14 PH QA uses its original `follow_days/event` implementation. They are separate recorded branches. Original full-precision L0 and era interface bridge calls must be used for Figure3 rather than rounded display CSVs. Figure3 retains 16 original analyses, order, labels, sample counts and plot parameters; CH26 adds omitted leading zeros without changing precision.

Current Bayesian HR functions retain the seven original priors, the original 2*1.959964 SE denominator and fixed clinical thresholds. `current_sources/bias/hr_catalog_current_ast.py` and the current Bayesian figure load HR/CI through `CH01_PRIMARY_CSV`, rather than stale numeric constants. This entry does not execute a new RMST prior catalog. Current QBA uses the already-approved common-outcome conversion, applied to observed and target HRs separately, with the true inverse for approximate HR equivalents. Set `QBA_MODE=common_outcome` explicitly; retain legacy mode only for old regression. For current weighted falsification results, do not reuse historical weighted rows through `QBA_ONLY=1`; the summary-only entry must point to the current versioned negative-control source after that branch is accepted.

## Current local numerical anchors for review

These are candidate summary facts, not fixed constants in the new summary/figure computation. They must be re-read from the same accepted result version when the package becomes a release.

| Item | Current NEE-corrected local candidate | Logical source |
|---|---|---|
| Primary population | N1502; HR1.05144237278742, CI0.92985279471207–1.18893126909961; RMST difference−1.17247324148331days | `CRRT_WORKSPACE_ROOT/09_工作稿/out/ch24_q4_corrected_primary/ch01_primary_dod_CORRECTED_CANDIDATE.csv` |
| Full cohort | N2463; HR1.11186346301827, CI1.00869995484273–1.22557788811216; RMST difference−2.10161682166659days | same two-population summary |
| Primary maximum post-weighting SMD | 0.023836968666115364, all nine variables below0.10 | `CRRT_WORKSPACE_ROOT/09_工作稿/out/ch24_q4_corrected_diagnostics/NEW/strictK3_v2_ipcw_balance_smd.csv` |
| Primary/full common-outcome E-values | 1.2267772147179932/1.362741730850597; CI-bound1/1.0838592916445193 | `CRRT_WORKSPACE_ROOT/09_工作稿/out/ch24_q4_corrected_qba/NEW/ch24_q5_evalues_summary.csv` |
| Seven HR priors | original three primary, four sensitivity; probabilities depend on authoritative HR/CI | `CRRT_WORKSPACE_ROOT/09_工作稿/out/ch24_q4_corrected_bayes/NEW/current_hr_prior_catalog_results.csv` |

Old values in the legacy code are deliberately historical. They must not be used to populate current figures or a manuscript. Patient/event counts and Table1 NEE are not changed by a blanket numeric replacement: the Table1 display NEE row independently reproduces the corrected baseline field and already agrees with its historical rich/FULL producer.

## Stage B completed offline interfaces

Primary CH01/Cox/RMST, original B1000 bootstrap/noon sensitivity, run13 FineGray/RFD, run18 extended RFD, both acid-base implementations, independent anchor outcomes and five-rescue outcomes now have actual portable inputs/outputs. See PORTABLE_IO_CONTRACT.md and STAGE_B.env.example for execution order and required external fields/layouts. All source functions and statistical expressions in the seven primary entry files pass the static AST gate. The rescue driver evaluates original saved AST rather than reconstructing a model. These exact portable entries remain unexecuted; static equivalence does not establish numerical replication.

## Inputs and outputs

Configure private directory roots, then inspect each source's input/output assignments before any execution. The current adapters retain the original branch directory layout beneath `CRRT_WORKSPACE_ROOT/09_工作稿/out/`; they are not a one-command fresh-install runner. The wrappers resolve included03a and03 via CRRT_CODE_ROOT; their path-token redirect is adapted to neutral env paths, with no double application of the /60 term. Generated clinical data must remain in private data directories, never this code directory.

Figure2 requires the current corrected CH01 clone RDS and matching summary. Figure3 requires the three current export interfaces specified by FIG3_FOREST_CSV, FIG3_QAPH_CSV and FIG3_PATH_RESCUE_CSV; its defaults are not substitutes for those interfaces. The actual original assembly and era/L0 bridge entry are now in `current_sources/figure3_interface/`; configure STAGE_C.env.example. All21 main/supplement display objects and exceptions are mapped in OBJECT_SOURCE_MAP.json. environment/ contains actual sessionInfo and installed versions, explicitly not a historical renv lockfile. Bayesian/QBA entries require an authoritative two-population HR/CI summary with `population`, `hr`, `ci_lower`, `ci_upper`; the primary key is `Equipoise SOFA 8-14`. Outputs use isolated branch directories. A fresh package does not contain those private inputs or OLD-regression artifacts.

## Table1 and S1 boundary

The basic `strictK3_v2_table1_mmr_spec.R` is not the complete rich/FULL Table1 generator. The supplied rich source generates the original rich table with nine missing rows; the original nine-row backfill SQL is also supplied. Patient-level augmentation/lab snapshots needed for exact full-table offline replay have not been located. The NEE display row alone passed its original SQL/function audit. CH28 approved S1's current model six-hour baseline NEE mean display (quantiles unchanged; missing0 with zero-filling limitation disclosed). Its exact original overall statistical producer remains unlocated. This package does not claim those full tables are end-to-end reproducible.

## Environment and verification boundary

Original computations used R4.5-family and the recorded R package families (survival, glmnet, dplyr/readr/tidyr, DBI/RPostgres, ggplot2/forestploter/ragg, sandwich/lmtest, zoo/lubridate); Python HR summaries use standard-library math/statistics/csv. Exact original versions/session information remain in private validation records. Static parse/hash/confidentiality checks on this portable candidate are distinct from the previously completed local branch regressions. No dependency lockfile, clean-room database run, container or end-to-end execution is claimed.

Before release: review the completed source/object map and explicit generator exceptions, then verify the remaining runtime/configuration boundaries, synchronize final accepted manuscript estimates and graph versions, rerun the package checks, bind an immutable code version to the actual public repository and verify public accessibility. Publication requires separate authorization.

CH29 adds Supplemental Table S9 (six existing pre-exposure diagnostic rows) to the object map, making21 objects with9 supplement tables. Its original validated prefix and source/result hashes are supplied, with the data CSV excluded. Figure2 current source includes only the approved leading-zero P-label format patch; historical sources remain unchanged.
