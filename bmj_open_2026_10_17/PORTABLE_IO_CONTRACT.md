# Stage B portable offline execution contract (candidate only)

No command below was executed during packaging. These entries preserve original statistical expressions and use externally supplied, credentialed, version-matched inputs. Copy/export the two neutral env templates; check every required variable is nonempty and every input exists before execution. All outputs must be outside the code directory, isolated from authoritative artifacts. The ZIP contains no data or example patient records.

| Entry under current_sources | Required external inputs | Output / execution order |
|---|---|---|
| primary/ch01_dod_primary.R | CRRT_V2_CCW; CRRT_DOD_RDS | CRRT_PRIMARY_OUT/ch01_primary_dod.csv, ccw_ch01_dod.rds and mortality ascertainment CSV. Original CH01 date/censoring assembly and Cox/RMST point estimates. |
| primary/ch01_rmst_sensitivity.R | CRRT_CH01_CCW produced by accepted primary version | Original stdout RMST B1000 full then equipoise, seed20260602; also original noon-date sensitivity. Capture stdout externally. No changes to bootstrap merge/resampling/CI. |
| secondary/run13/ORIGINAL_OFFLINE_PATH_ADAPTER.R | CRRT_V2_CCW; CRRT_SECONDARY_CSV; CRRT_DOD_RDS | CRRT_SECONDARY_OUT/run13; original FineGray and RFD outputs. |
| secondary/run18/ORIGINAL_OFFLINE_PATH_ADAPTER.R | same three inputs, plus same-version run13 RFD result | CRRT_SECONDARY_OUT/run18; retains original run13 numerical cross-check. Run13 first. |
| acidbase/median/NEW_adapter.R | CRRT_V2_DATA_DIR/{baseline.rds,tv_hourly.rds}; CRRT_ACIDBASE_RDS; CRRT_DOD_RDS | CRRT_ACIDBASE_OUT/median; pass original --grace=24 to reproduce the reported branch. Original default 12/24 and --tv48 argument retained. This entry refits original weights and model only if explicitly run. |
| acidbase/missing_indicator/NEW_adapter.R | same as median | CRRT_ACIDBASE_OUT/missing_indicator; --grace=24; original zero+indicator implementation retained. |
| anchor/NEW_selected_original_AST.R | independent CRRT_ANCHOR_MINUS6_CCW, CRRT_ANCHOR_STAGE2_CCW; CRRT_DOD_RDS | CRRT_ANCHOR_OUT/anchor_original_eight_rows.csv. Eight original old-outcome/current-CH01 calls remain labelled explicitly; select CH01 rows for current display. |
| rescue/run_rescue_existing_weights.R | CRRT_CODE_ROOT; CRRT_RESCUE_WEIGHT_ROOT; CRRT_DOD_RDS | CRRT_RESCUE_OUT/strictK3_v2_rescue_sensitivity_model_results.csv; only original Cox functions/final bind_rows/rounding evaluated. Does not run make_variant or regenerate weights. |

Rescue layout is `${CRRT_RESCUE_WEIGHT_ROOT}/rescue_${name}/NEW/ccw_final.rds`, for the five original names: no_rescue_exclusion; sensitivity_R1_R2_R2b_R3_R4; sensitivity_R1_confirm_R2_R2b_R3; sensitivity_liberal_thresholds_R1_R2_R2b_R3; sensitivity_conservative_thresholds_R1_R2_R2b_R3. Preserve each branch's own cohort, protocol censoring and weights. The main branch is not substituted into this loop.

## External schema and provenance

These are required column contracts, not fabricated example records. Additional original columns must be retained. Private RDS must preserve original R classes, factor levels, row order and numeric precision. `PRIVATE_INPUT_SCHEMA.json` (code-only schema, if present) records complete names/classes from the existing accepted local snapshots without IDs or values.

- DOD RDS: unique `subject_id` and `dod` (Date). Derive from credentialed MIMIC patients, no manual filling. Every branch subject must be covered by the applicable lookup. Do not ship this lookup.
- Main ccw: original `clone_id`, `stay_id`, `subject_id`, `strategy`, `t0`, `deathtime`, `protocol_censor_h`, `ipcw_w`, `ipcw_w_trim`, `equipoise`, `early_crrt`, `died_28d`, `follow_time_28d`, plus original covariates and `days_to_death` required by secondary. Current NEE correction must already be applied in this supplied ccw. Primary original assembly derives CH01 corrected fields itself.
- CH01 ccw: prior output plus `strategy_f`, `corrected_follow_days`, `corrected_event`, `death_days_inhospital`, `dod`; preserve all original date and censoring columns for noon sensitivity. Primary full/equipoise outputs share this same input version.
- Secondary CSV: original unique `stay_id`; `early_crrt`, `died_28d`, `renal_recovery_28d`, `crrt_total_days`, `rrt_free_days_28d`. Original boolean coercion is retained; original CH01 secondary_NEW.sql is included for credentialed extraction rather than manufactured replacement endpoints.
- Baseline/hourly: original stage-v2 baseline and hourly schemas consumed by included 03/03a and acidbase formulas; preserve all baseline eligibility/date columns. Acidbase RDS unique `stay_id`, `ph_t0`, `k_t0`, `hco3_t0` from the original closest-to-t0 extraction. This entry does not extract new labs.
- Anchor ccw: respective original anchor cohort with `follow_days`, `event`, `protocol_censor_h`, `equipoise`, date/strategy/weight identifiers; main-cohort ccw is not interchangeable.
- Rescue ccw: original five independent branch ccw with `died_28d`, `follow_time_28d`, `follow_days`, `event`, original covariates and counts. The portable driver replays original CH01 date recoding as an exact assertion before Cox; it stops on mismatch.

## Verification and unresolved work

`STATIC_AST_GATE.R` parses only. It compares full original function ASTs and the ordered non-I/O top-level expression sequence after explicit path normalization/configuration exclusions. Rescue definitions are exact original source copies; the driver evaluates only the original final bind_rows expression and original rounding. This demonstrates source equivalence, not successful runtime or numerical reproducibility of the portable package. Original local branch OLD→NEW evidence remains separate from this packaging check. A clean configured run still needs authorized credentialed data, the original R environment/dependencies, compatible original snapshots, and branch-specific original weight builders. No database, weight fit, Cox/FineGray, RMST bootstrap or clinical analysis was run while building this candidate.

Remaining release gaps: portable Figure3 full-precision assembly/era-L0 bridge; complete branch input-production/DB configuration wiring; rich/FULL Table1 augmentation and unlocated S1 producer; runtime regression of this exact portable tree; dependency lock; approved immutable public version. Public repository availability of a previous snapshot does not establish publication of this candidate. No upload/push/repository creation is authorized here.
