# =============================================================================
# Strict clean first-KDIGO3 completion audit and manuscript-facing summary
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
fig_dir <- file.path(strict_dir, "figures")

primary <- read_csv(file.path(results_dir, "strictK3_primary_v2_mortality_results.csv"), show_col_types = FALSE)
finegray <- read_csv(file.path(results_dir, "strictK3_v2_weighted_finegray_results.csv"), show_col_types = FALSE)
rfd <- read_csv(file.path(results_dir, "strictK3_v2_rfd_results.csv"), show_col_types = FALSE)
rfd_qa <- read_csv(file.path(results_dir, "strictK3_v2_QA_rfd_model_comparison.csv"), show_col_types = FALSE)
ph_qa <- read_csv(file.path(results_dir, "strictK3_v2_QA_ph_treatment_results.csv"), show_col_types = FALSE)
rescue <- read_csv(file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.csv"), show_col_types = FALSE)
anchor <- read_csv(file.path(results_dir, "strictK3_v2_anchor_sensitivity_comparison_results.csv"), show_col_types = FALSE)
era_adj <- read_csv(file.path(results_dir, "strictK3_v2_calendar_era_adjusted_results.csv"), show_col_types = FALSE)
funnel <- read_csv(file.path(data_dir, "strictK3_primary_v2_funnel.csv"), show_col_types = FALSE)
weights <- read_csv(file.path(results_dir, "strictK3_v2_ipcw_weight_summary.csv"), show_col_types = FALSE)
balance <- read_csv(file.path(results_dir, "strictK3_v2_ipcw_balance_smd.csv"), show_col_types = FALSE)
path_sens <- read_csv(file.path(results_dir, "strictK3_v2_path_adjusted_ipcw_results.csv"), show_col_types = FALSE)
v1v2 <- read_csv(file.path(results_dir, "strictK3_v2_vs_v1_primary_sofa8_14.csv"), show_col_types = FALSE)

fmt_hr <- function(hr, lo, hi, p) {
  sprintf("%.3f (%.3f to %.3f), p=%s", hr, lo, hi, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

primary_eq <- primary %>% filter(analysis == "Primary SOFA 8-14 equipoise")
primary_full <- primary %>% filter(analysis == "Full strictK3_v2 primary cohort")
fg_eq <- finegray %>% filter(Analysis == "Equipoise SOFA 8-14")
rfd_eq <- rfd %>% filter(Analysis == "Equipoise SOFA 8-14")
rfd_corrected <- rfd_qa %>% filter(grepl("^Corrected", analysis))
ph_eq_strat <- ph_qa %>% filter(model == "SOFA 8-14 era-stratified")
anchor_sofa <- anchor %>% filter(grepl("SOFA 8-14$", analysis))
rescue_sofa <- rescue %>% filter(grepl("SOFA 8-14$", analysis))
era_adj_eq <- era_adj %>% filter(grepl("SOFA 8-14", analysis))
balance_eq <- balance %>% filter(population == "SOFA 8-14")
weights_eq <- weights %>% filter(equipoise %in% TRUE | equipoise == "TRUE")

audit <- tibble::tribble(
  ~requirement, ~status, ~evidence,
  "Strict clean first-KDIGO3 cohort extraction with R1+R2+R2b+R3 rescue exclusion", "complete", paste0("strictK3_primary_extract_v2.R; strictK3_primary_v2_extraction_report.md; cohort_base n=", primary_full$n_patients, "; SOFA 8-14 n=", primary_eq$n_patients, "."),
  "Funnel counts with separate rescue exclusion step", "complete", "strictK3_primary_v2_funnel.csv and extraction report.",
  "Table 1 early vs deferred with SMD flags", "complete", "strictK3_v2_table1_sofa8_14.csv/md; SMD>=0.10 flags included.",
  "IPCW weight diagnostics, ESS, balance and positivity", "complete", "strictK3_v2_ipcw_weight_summary.csv; strictK3_v2_ipcw_balance_smd.csv; strictK3_v2_ipcw_diagnostics_report.md; Love plot.",
  "Primary weighted Cox and RMST", "complete", "strictK3_primary_v2_mortality_results.csv/md.",
  "E-values", "complete", "strictK3_primary_v2_mortality_results.csv includes point and CI-bound E-values.",
  "Weighted Fine-Gray and RRT-free days", "complete", "strictK3_v2_weighted_finegray_results.csv; strictK3_v2_rfd_results.csv; strictK3_v2_QA_rfd_report.md documents corrected strategy-compatible RFD.",
  "SOFA strata and continuous/band scan", "complete", "strictK3_v2_sofa_strata_results.csv; strictK3_v2_sofa_band_scan_results.csv; heatmap figure.",
  "Calendar-era strata, era-adjusted Cox, and PH QA", "complete", "strictK3_v2_calendar_era_results.csv; strictK3_v2_calendar_era_adjusted_results.csv; strictK3_v2_QA_ph_report.md.",
  "Anchor sensitivity: KDIGO3-6h and KDIGO2", "complete", "strictK3_v2_anchor_sensitivity_models.R; strictK3_v2_anchor_sensitivity_report.md; strictK3_v2_anchor_sensitivity_diagnostics_report.md; independent CCW/IPCW directories.",
  "Rescue sensitivity: no rescue, R1-R4, R1-confirm, threshold variants", "complete", "strictK3_v2_rescue_sensitivity_models.R; strictK3_v2_rescue_sensitivity_model_results.csv.",
  "L0-only denominator sensitivity", "complete", "strictK3_v2_l0_only_sensitivity_results.csv; strictK3_v2_l0_only_sensitivity_report.md.",
  "AKI staging-path sensitivity and path-adjusted IPCW", "complete", "strictK3_path_sensitivity_v2.R; strictK3_v2_path_sensitivity_report.md; stage3_path included in alternative IPCW numerator/denominator.",
  "v1 contaminated vs v2 clean-t0 comparison", "complete", "strictK3_v1_v2_comparison.R; strictK3_v2_vs_v1_comparison_report.md.",
  "Figures: forest, cumulative mortality, Love plot, SOFA scan, RFD", "complete", "strictK3_primary_v2/figures/*.png/*.pdf.",
  "Assumptions/ambiguities documented", "complete", "This report; anchor sensitivity report notes no-rescue guardrail for KDIGO2/-6h."
)

write_csv(audit, file.path(results_dir, "strictK3_v2_completion_audit_checklist.csv"))

report <- c(
  "# Strict clean first-KDIGO3 TTE/CCW/IPCW completion audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This audit covers the requested MIMIC-IV v3.1 strict clean first-KDIGO3 target-trial emulation. Existing legacy scripts were not edited; strict outputs are under `analysis/strictK3_primary_v2/` and companion scripts are named `strictK3_v2_*`.",
  "",
  "## Assumptions / ambiguities",
  "",
  "- Clean t0 is defined as first physiologic KDIGO stage 3 using `aki_stage_creat = 3 OR aki_stage_uo = 3`; `aki_stage_crrt` is explicitly not used to define t0.",
  "- Prior CRRT exclusion is strict: `crrt_first_time < clean_t0` is excluded; CRRT at or after clean t0 remains eligible as the exposure process.",
  "- Rescue dictionary v1.0 is applied to the strict clean first-KDIGO3 primary analysis: R1 OR (R2 OR R2b) OR R3. R4/R5/R6 are not primary exclusions.",
  "- Anchor sensitivity models for clean KDIGO3-minus-6h and clean KDIGO2 are run as no-rescue anchor comparisons. This is intentional because `rescue_dictionary_v1.0_frozen.md` is defined for strict clean first-KDIGO3 t0 and should not be silently imposed on KDIGO2 or the -6h anchor.",
  "- RRT-free days use distinct CRRT calendar days from `mimiciv_derived.rrt` within 28 days, preserving the prior project definition.",
  "- Fixed-28d RFD is reported only with the QA-corrected strategy-compatible clone rows; the naive all-clone fixed-time RFD estimate is documented as non-causal/artifactual.",
  "- Table 1 reports observed early vs deferred baseline imbalance; IPCW-corrected balance is shown separately in the Love plot and SMD diagnostics.",
  "",
  "## Database/schema evidence",
  "",
  "- Source database: local PostgreSQL `mimiciv31`.",
  "- Main source schemas/tables: `mimiciv_icu`, `mimiciv_hosp`, `mimiciv_derived`.",
  "- Strict derived schema: `data_extract_crrt_strictk3_v2` (`cohort_pre_rescue_all`, `cohort_base`, `baseline`, `secondary_outcomes`).",
  "",
  "## Key results",
  "",
  paste0("- Primary SOFA 8-14 cohort: n=", primary_eq$n_patients, ", early n=", primary_eq$early_n, "."),
  paste0("- Primary weighted Cox HR: ", fmt_hr(primary_eq$weighted_hr, primary_eq$ci_lower, primary_eq$ci_upper, primary_eq$p_value), "."),
  paste0("- Primary RMST difference (Early minus Deferred): ", sprintf("%.2f days (%.2f to %.2f)", primary_eq$rmst_diff, primary_eq$rmst_diff_ci_lower, primary_eq$rmst_diff_ci_upper), "."),
  paste0("- E-value point / CI-bound: ", sprintf("%.2f / %.2f", primary_eq$e_value_point, primary_eq$e_value_ci_bound), "."),
  paste0("- Renal recovery Fine-Gray sHR in SOFA 8-14: ", fmt_hr(fg_eq$HR[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$HR_lo[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$HR_hi[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$p_value[fg_eq$Endpoint == "Renal recovery (death competing)"]), "."),
  paste0("- Corrected strategy-compatible RFD weighted mean difference in SOFA 8-14: ", sprintf("%.3f days (%.3f to %.3f), p=%s", rfd_corrected$mean_diff, rfd_corrected$ci_lower, rfd_corrected$ci_upper, rfd_corrected$p_display), "."),
  paste0("- Era-stratified Cox PH-QA SOFA 8-14 HR: ", sprintf("%.3f (%.3f to %.3f), p=%s", ph_eq_strat$HR, ph_eq_strat$CI_lo, ph_eq_strat$CI_hi, ph_eq_strat$p_display), "."),
  paste0("- Max post-IPCW SMD in SOFA 8-14: ", sprintf("%.3f", max(balance_eq$post_weighted_smd, na.rm = TRUE)), "."),
  paste0("- Path-adjusted IPCW SOFA 8-14 HR: ", {
    x <- path_sens %>% filter(grepl("SOFA 8-14", analysis))
    fmt_hr(x$HR, x$CI_lo, x$CI_hi, x$p)
  }, "."),
  "",
  "## Difference versus clean KDIGO3-minus-6h and KDIGO2 anchors",
  "",
  paste0("- Strict primary SOFA 8-14: ", fmt_hr(primary_eq$weighted_hr, primary_eq$ci_lower, primary_eq$ci_upper, primary_eq$p_value), "."),
  paste0("- Clean KDIGO3-minus-6h SOFA 8-14 anchor sensitivity: ", {
    x <- anchor_sofa %>% filter(grepl("KDIGO3_minus6h", analysis))
    fmt_hr(x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
  }, "."),
  paste0("- KDIGO2 exploratory SOFA 8-14 anchor sensitivity: ", {
    x <- anchor_sofa %>% filter(grepl("^KDIGO2", analysis))
    fmt_hr(x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
  }, "."),
  "",
  "Interpretation: the strict clean first-KDIGO3 primary analysis remains null in the prespecified SOFA 8-14 equipoise band. The clean -6h anchor is also null in SOFA 8-14, while KDIGO2 becomes a different, exploratory question and shows harm-direction association rather than support for early CRRT.",
  "",
  "## v1 contaminated vs v2 clean-t0 primary shift",
  "",
  paste0("- v1 contaminated SOFA 8-14: n=", v1v2$n_patients[v1v2$version == "v1 contaminated t0"], ", early n=", v1v2$early_n[v1v2$version == "v1 contaminated t0"], ", HR=", sprintf("%.3f", v1v2$HR[v1v2$version == "v1 contaminated t0"]), "."),
  paste0("- v2 clean SOFA 8-14: n=", v1v2$n_patients[v1v2$version == "v2 clean Cr/UO-only t0"], ", early n=", v1v2$early_n[v1v2$version == "v2 clean Cr/UO-only t0"], ", HR=", sprintf("%.3f", v1v2$HR[v1v2$version == "v2 clean Cr/UO-only t0"]), "."),
  "- Full comparison report: `strictK3_v2_vs_v1_comparison_report.md`.",
  "",
  "## Requirement checklist",
  "",
  "| requirement | status | evidence |",
  "|---|---|---|",
  apply(audit, 1, function(x) sprintf("| %s | %s | %s |", x[["requirement"]], x[["status"]], x[["evidence"]])),
  "",
  "## Core output index",
  "",
  paste0("- Cohort extraction: `", file.path(strict_dir, "strictK3_primary_v2_extraction_report.md"), "`"),
  paste0("- Primary mortality: `", file.path(results_dir, "strictK3_primary_v2_mortality_report.md"), "`"),
  paste0("- IPCW diagnostics: `", file.path(results_dir, "strictK3_v2_ipcw_diagnostics_report.md"), "`"),
  paste0("- Secondary outcomes: `", file.path(results_dir, "strictK3_v2_weighted_finegray_rfd_report.md"), "`"),
  paste0("- Corrected RFD QA: `", file.path(results_dir, "strictK3_v2_QA_rfd_report.md"), "`"),
  paste0("- Cox PH QA: `", file.path(results_dir, "strictK3_v2_QA_ph_report.md"), "`"),
  paste0("- SOFA/era: `", file.path(results_dir, "strictK3_v2_sofa_era_sensitivity_report.md"), "`"),
  paste0("- Anchor sensitivity: `", file.path(results_dir, "strictK3_v2_anchor_sensitivity_report.md"), "`"),
  paste0("- Anchor overlap/balance diagnostics: `", file.path(results_dir, "strictK3_v2_anchor_sensitivity_diagnostics_report.md"), "`"),
  paste0("- Rescue sensitivity: `", file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_report.md"), "`"),
  paste0("- L0-only sensitivity: `", file.path(results_dir, "strictK3_v2_l0_only_sensitivity_report.md"), "`"),
  paste0("- Path sensitivity: `", file.path(results_dir, "strictK3_v2_path_sensitivity_report.md"), "`"),
  paste0("- v1-v2 comparison: `", file.path(results_dir, "strictK3_v2_vs_v1_comparison_report.md"), "`"),
  paste0("- Tables/figures: `", file.path(results_dir, "strictK3_v2_tables_figures_summary.md"), "`")
)

writeLines(report, file.path(results_dir, "strictK3_v2_completion_audit_report.md"))
cat("Wrote completion audit:", file.path(results_dir, "strictK3_v2_completion_audit_report.md"), "\n")
