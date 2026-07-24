# =============================================================================
# Strict first-KDIGO3 completion audit and manuscript-facing summary
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
fig_dir <- file.path(strict_dir, "figures")

primary <- read_csv(file.path(results_dir, "strictK3_primary_mortality_results.csv"), show_col_types = FALSE)
finegray <- read_csv(file.path(results_dir, "strictK3_weighted_finegray_results.csv"), show_col_types = FALSE)
rfd <- read_csv(file.path(results_dir, "strictK3_rfd_results.csv"), show_col_types = FALSE)
rescue <- read_csv(file.path(results_dir, "strictK3_rescue_sensitivity_model_results.csv"), show_col_types = FALSE)
anchor <- read_csv(file.path(results_dir, "strictK3_anchor_sensitivity_comparison_results.csv"), show_col_types = FALSE)
era_adj <- read_csv(file.path(results_dir, "strictK3_calendar_era_adjusted_results.csv"), show_col_types = FALSE)
funnel <- read_csv(file.path(data_dir, "strictK3_primary_funnel.csv"), show_col_types = FALSE)
weights <- read_csv(file.path(results_dir, "strictK3_ipcw_weight_summary.csv"), show_col_types = FALSE)
balance <- read_csv(file.path(results_dir, "strictK3_ipcw_balance_smd.csv"), show_col_types = FALSE)

fmt_hr <- function(hr, lo, hi, p) {
  sprintf("%.3f (%.3f to %.3f), p=%s", hr, lo, hi, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

primary_eq <- primary %>% filter(analysis == "Primary SOFA 8-14 equipoise")
primary_full <- primary %>% filter(analysis == "Full strictK3 primary cohort")
fg_eq <- finegray %>% filter(Analysis == "Equipoise SOFA 8-14")
rfd_eq <- rfd %>% filter(Analysis == "Equipoise SOFA 8-14")
anchor_sofa <- anchor %>% filter(grepl("SOFA 8-14$", analysis))
rescue_sofa <- rescue %>% filter(grepl("SOFA 8-14$", analysis))
era_adj_eq <- era_adj %>% filter(grepl("SOFA 8-14", analysis))
balance_eq <- balance %>% filter(population == "SOFA 8-14")
weights_eq <- weights %>% filter(equipoise %in% TRUE | equipoise == "TRUE")

audit <- tibble::tribble(
  ~requirement, ~status, ~evidence,
  "Strict first-KDIGO3 cohort extraction with R1+R2+R2b+R3 rescue exclusion", "complete", "strictK3_primary_extract.R; strictK3_primary_extraction_report.md; cohort_base n=2566; SOFA 8-14 n=1568.",
  "Funnel counts with separate rescue exclusion step", "complete", "strictK3_primary_funnel.csv and extraction report.",
  "Table 1 early vs deferred with SMD flags", "complete", "strictK3_table1_sofa8_14.csv/md; SMD>=0.10 flags included.",
  "IPCW weight diagnostics, ESS, balance and positivity", "complete", "strictK3_ipcw_weight_summary.csv; strictK3_ipcw_balance_smd.csv; strictK3_ipcw_diagnostics_report.md; Love plot.",
  "Primary weighted Cox and RMST", "complete", "strictK3_primary_mortality_results.csv/md.",
  "E-values", "complete", "strictK3_primary_mortality_results.csv includes point and CI-bound E-values.",
  "Weighted Fine-Gray and RRT-free days", "complete", "strictK3_weighted_finegray_results.csv; strictK3_rfd_results.csv; strictK3_weighted_finegray_rfd_report.md.",
  "SOFA strata and continuous/band scan", "complete", "strictK3_sofa_strata_results.csv; strictK3_sofa_band_scan_results.csv; heatmap figure.",
  "Calendar-era strata and era-adjusted Cox", "complete", "strictK3_calendar_era_results.csv; strictK3_calendar_era_adjusted_results.csv.",
  "Anchor sensitivity: KDIGO3-6h and KDIGO2", "complete", "strictK3_anchor_sensitivity_models.R; strictK3_anchor_sensitivity_report.md; strictK3_anchor_sensitivity_diagnostics_report.md; independent CCW/IPCW directories.",
  "Rescue sensitivity: no rescue, R1-R4, R1-confirm, threshold variants", "complete", "strictK3_rescue_sensitivity_models.R; strictK3_rescue_sensitivity_model_results.csv.",
  "L0-only denominator sensitivity", "complete", "strictK3_l0_only_sensitivity_results.csv; strictK3_l0_only_sensitivity_report.md.",
  "Figures: forest, cumulative mortality, Love plot, SOFA scan, RFD", "complete", "strictK3_primary/figures/*.png/*.pdf.",
  "Assumptions/ambiguities documented", "complete", "This report; anchor sensitivity report notes no-rescue guardrail for KDIGO2/-6h."
)

write_csv(audit, file.path(results_dir, "strictK3_completion_audit_checklist.csv"))

report <- c(
  "# Strict first-KDIGO3 TTE/CCW/IPCW completion audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This audit covers the requested MIMIC-IV v3.1 strict first-KDIGO3 target-trial emulation. Existing legacy scripts were not edited; strict outputs are under `analysis/strictK3_primary/` and companion scripts are named `strictK3_*`.",
  "",
  "## Assumptions / ambiguities",
  "",
  "- Rescue dictionary v1.0 is applied to the strict first-KDIGO3 primary analysis: R1 OR (R2 OR R2b) OR R3. R4/R5/R6 are not primary exclusions.",
  "- Anchor sensitivity models for KDIGO3-minus-6h and KDIGO2 are run as no-rescue anchor comparisons. This is intentional because `rescue_dictionary_v1.0_frozen.md` is defined for strict first-KDIGO3 t0 and should not be silently imposed on KDIGO2 or the legacy -6h anchor.",
  "- RRT-free days use distinct CRRT calendar days from `mimiciv_derived.rrt` within 28 days, preserving the prior project definition.",
  "- Table 1 reports observed early vs deferred baseline imbalance; IPCW-corrected balance is shown separately in the Love plot and SMD diagnostics.",
  "",
  "## Database/schema evidence",
  "",
  "- Source database: local PostgreSQL `mimiciv31`.",
  "- Main source schemas/tables: `mimiciv_icu`, `mimiciv_hosp`, `mimiciv_derived`.",
  "- Strict derived schema: `data_extract_crrt_strictk3` (`cohort_pre_rescue_all`, `cohort_base`, `baseline`, `secondary_outcomes`).",
  "",
  "## Key results",
  "",
  paste0("- Primary SOFA 8-14 cohort: n=", primary_eq$n_patients, ", early n=", primary_eq$early_n, "."),
  paste0("- Primary weighted Cox HR: ", fmt_hr(primary_eq$weighted_hr, primary_eq$ci_lower, primary_eq$ci_upper, primary_eq$p_value), "."),
  paste0("- Primary RMST difference (Early minus Deferred): ", sprintf("%.2f days (%.2f to %.2f)", primary_eq$rmst_diff, primary_eq$rmst_diff_ci_lower, primary_eq$rmst_diff_ci_upper), "."),
  paste0("- E-value point / CI-bound: ", sprintf("%.2f / %.2f", primary_eq$e_value_point, primary_eq$e_value_ci_bound), "."),
  paste0("- Renal recovery Fine-Gray sHR in SOFA 8-14: ", fmt_hr(fg_eq$HR[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$HR_lo[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$HR_hi[fg_eq$Endpoint == "Renal recovery (death competing)"], fg_eq$p_value[fg_eq$Endpoint == "Renal recovery (death competing)"]), "."),
  paste0("- RFD weighted mean difference in SOFA 8-14: ", sprintf("%.3f days (%.3f to %.3f), p=%s", rfd_eq$Weighted_mean_diff, rfd_eq$Diff_lo, rfd_eq$Diff_hi, rfd_eq$Weighted_p_display), "."),
  paste0("- Max post-IPCW SMD in SOFA 8-14: ", sprintf("%.3f", max(balance_eq$post_weighted_smd, na.rm = TRUE)), "."),
  "",
  "## Difference versus legacy KDIGO3-minus-6h",
  "",
  paste0("- Strict primary SOFA 8-14: ", fmt_hr(primary_eq$weighted_hr, primary_eq$ci_lower, primary_eq$ci_upper, primary_eq$p_value), "."),
  paste0("- Legacy KDIGO3-minus-6h SOFA 8-14 anchor sensitivity: ", {
    x <- anchor_sofa %>% filter(grepl("KDIGO3_minus6h", analysis))
    fmt_hr(x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
  }, "."),
  paste0("- KDIGO2 exploratory SOFA 8-14 anchor sensitivity: ", {
    x <- anchor_sofa %>% filter(grepl("^KDIGO2", analysis))
    fmt_hr(x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
  }, "."),
  "",
  "Interpretation: the strict first-KDIGO3 primary analysis remains null in the prespecified SOFA 8-14 equipoise band. The legacy -6h anchor is also null in SOFA 8-14, while KDIGO2 becomes a different, exploratory question and shows harm-direction association rather than support for early CRRT.",
  "",
  "## Requirement checklist",
  "",
  "| requirement | status | evidence |",
  "|---|---|---|",
  apply(audit, 1, function(x) sprintf("| %s | %s | %s |", x[["requirement"]], x[["status"]], x[["evidence"]])),
  "",
  "## Core output index",
  "",
  paste0("- Cohort extraction: `", file.path(strict_dir, "strictK3_primary_extraction_report.md"), "`"),
  paste0("- Primary mortality: `", file.path(results_dir, "strictK3_primary_mortality_report.md"), "`"),
  paste0("- IPCW diagnostics: `", file.path(results_dir, "strictK3_ipcw_diagnostics_report.md"), "`"),
  paste0("- Secondary outcomes: `", file.path(results_dir, "strictK3_weighted_finegray_rfd_report.md"), "`"),
  paste0("- SOFA/era: `", file.path(results_dir, "strictK3_sofa_era_sensitivity_report.md"), "`"),
  paste0("- Anchor sensitivity: `", file.path(results_dir, "strictK3_anchor_sensitivity_report.md"), "`"),
  paste0("- Anchor overlap/balance diagnostics: `", file.path(results_dir, "strictK3_anchor_sensitivity_diagnostics_report.md"), "`"),
  paste0("- Rescue sensitivity: `", file.path(results_dir, "strictK3_rescue_sensitivity_model_report.md"), "`"),
  paste0("- L0-only sensitivity: `", file.path(results_dir, "strictK3_l0_only_sensitivity_report.md"), "`"),
  paste0("- Tables/figures: `", file.path(results_dir, "strictK3_tables_figures_summary.md"), "`")
)

writeLines(report, file.path(results_dir, "strictK3_completion_audit_report.md"))
cat("Wrote completion audit:", file.path(results_dir, "strictK3_completion_audit_report.md"), "\n")
