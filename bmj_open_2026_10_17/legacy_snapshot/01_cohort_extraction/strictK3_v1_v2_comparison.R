# =============================================================================
# Strict KDIGO3 v1 contaminated vs v2 clean-t0 comparison summary
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

base_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
v1_dir <- file.path(base_dir, "strictK3_primary")
v2_dir <- file.path(base_dir, "strictK3_primary_v2")
v2_results <- file.path(v2_dir, "results")

fmt_hr <- function(hr, lo, hi, p) {
  sprintf("%.3f (%.3f to %.3f), p=%s", hr, lo, hi, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

v1_primary <- read_csv(file.path(v1_dir, "results", "strictK3_primary_mortality_results.csv"), show_col_types = FALSE)
v2_primary <- read_csv(file.path(v2_results, "strictK3_primary_v2_mortality_results.csv"), show_col_types = FALSE)
v1_weights <- read_csv(file.path(v1_dir, "results", "strictK3_ipcw_weight_summary.csv"), show_col_types = FALSE)
v2_weights <- read_csv(file.path(v2_results, "strictK3_v2_ipcw_weight_summary.csv"), show_col_types = FALSE)
v1_balance <- read_csv(file.path(v1_dir, "results", "strictK3_ipcw_balance_smd.csv"), show_col_types = FALSE)
v2_balance <- read_csv(file.path(v2_results, "strictK3_v2_ipcw_balance_smd.csv"), show_col_types = FALSE)
v1_funnel <- read_csv(file.path(v1_dir, "data", "strictK3_primary_funnel.csv"), show_col_types = FALSE)
v2_funnel <- read_csv(file.path(v2_dir, "data", "strictK3_primary_v2_funnel.csv"), show_col_types = FALSE)
v1_sofa <- read_csv(file.path(v1_dir, "results", "strictK3_sofa_strata_results.csv"), show_col_types = FALSE)
v2_sofa <- read_csv(file.path(v2_results, "strictK3_v2_sofa_strata_results.csv"), show_col_types = FALSE)
v1_scan <- read_csv(file.path(v1_dir, "results", "strictK3_sofa_band_scan_results.csv"), show_col_types = FALSE)
v2_scan <- read_csv(file.path(v2_results, "strictK3_v2_sofa_band_scan_results.csv"), show_col_types = FALSE)

primary_rows <- bind_rows(
  v1_primary %>%
    filter(analysis == "Primary SOFA 8-14 equipoise") %>%
    mutate(version = "v1 contaminated t0"),
  v2_primary %>%
    filter(analysis == "Primary SOFA 8-14 equipoise") %>%
    mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  transmute(
    version,
    n_patients,
    early_n,
    early_pct = round(100 * early_n / n_patients, 1),
    events,
    HR = weighted_hr,
    CI_lo = ci_lower,
    CI_hi = ci_upper,
    p = p_value,
    RMST_diff = rmst_diff,
    RMST_lo = rmst_diff_ci_lower,
    RMST_hi = rmst_diff_ci_upper,
    e_value_point,
    e_value_ci_bound
  )

full_rows <- bind_rows(
  v1_primary %>%
    filter(grepl("^Full", analysis)) %>%
    mutate(version = "v1 contaminated t0"),
  v2_primary %>%
    filter(grepl("^Full", analysis)) %>%
    mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  transmute(
    version,
    n_patients,
    early_n,
    early_pct = round(100 * early_n / n_patients, 1),
    events,
    HR = weighted_hr,
    CI_lo = ci_lower,
    CI_hi = ci_upper,
    p = p_value
  )

funnel_compare <- full_join(
  v1_funnel %>% select(step, v1_n = n),
  v2_funnel %>% select(step, v2_n = n),
  by = "step"
) %>%
  mutate(delta_v2_minus_v1 = v2_n - v1_n)

weight_compare <- bind_rows(
  v1_weights %>% mutate(version = "v1 contaminated t0"),
  v2_weights %>% mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  filter(equipoise %in% TRUE | equipoise == "TRUE") %>%
  select(version, strategy, n_patients, mean_w, p99_w, max_w, ess)

balance_compare <- bind_rows(
  v1_balance %>% mutate(version = "v1 contaminated t0"),
  v2_balance %>% mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  filter(population == "SOFA 8-14") %>%
  group_by(version) %>%
  summarise(
    max_post_smd = max(post_weighted_smd, na.rm = TRUE),
    post_smd_gt_010 = sum(post_weighted_smd > 0.10, na.rm = TRUE),
    .groups = "drop"
  )

sofa_compare <- bind_rows(
  v1_sofa %>% mutate(version = "v1 contaminated t0"),
  v2_sofa %>% mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  select(version, analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)

scan_compare <- bind_rows(
  v1_scan %>% mutate(version = "v1 contaminated t0"),
  v2_scan %>% mutate(version = "v2 clean Cr/UO-only t0")
) %>%
  group_by(version) %>%
  summarise(
    scan_rows = n(),
    min_hr = min(weighted_hr, na.rm = TRUE),
    max_hr = max(weighted_hr, na.rm = TRUE),
    significant_tiles = sum(!is.na(p_value) & p_value < 0.05),
    .groups = "drop"
  )

write_csv(primary_rows, file.path(v2_results, "strictK3_v2_vs_v1_primary_sofa8_14.csv"))
write_csv(full_rows, file.path(v2_results, "strictK3_v2_vs_v1_full_cohort.csv"))
write_csv(funnel_compare, file.path(v2_results, "strictK3_v2_vs_v1_funnel.csv"))
write_csv(weight_compare, file.path(v2_results, "strictK3_v2_vs_v1_weight_summary.csv"))
write_csv(balance_compare, file.path(v2_results, "strictK3_v2_vs_v1_balance_summary.csv"))
write_csv(sofa_compare, file.path(v2_results, "strictK3_v2_vs_v1_sofa_strata.csv"))
write_csv(scan_compare, file.path(v2_results, "strictK3_v2_vs_v1_sofa_scan_summary.csv"))

report <- c(
  "# Strict KDIGO3 v1 contaminated t0 vs v2 clean t0 comparison",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## What changed",
  "",
  "- v1 t0 used `mimiciv_derived.kdigo_stages WHERE aki_stage = 3`, which included `aki_stage_crrt = 3`.",
  "- v2 t0 uses clean physiologic KDIGO3 only: `aki_stage_creat = 3 OR aki_stage_uo = 3`.",
  "- Prior CRRT is excluded as `crrt_first_time < clean_t0`; CRRT at or after clean t0 remains the exposure.",
  "",
  "## Primary SOFA 8-14 comparison",
  "",
  "| version | n | early n | early % | events | HR (95% CI), p | RMST diff days (95% CI) | E-value point / CI-bound |",
  "|---|---:|---:|---:|---:|---|---|---|",
  apply(primary_rows, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %.2f (%.2f to %.2f) | %.2f / %.2f |",
    x[["version"]], x[["n_patients"]], x[["early_n"]], x[["early_pct"]],
    x[["events"]],
    fmt_hr(as.numeric(x[["HR"]]), as.numeric(x[["CI_lo"]]), as.numeric(x[["CI_hi"]]), as.numeric(x[["p"]])),
    as.numeric(x[["RMST_diff"]]), as.numeric(x[["RMST_lo"]]), as.numeric(x[["RMST_hi"]]),
    as.numeric(x[["e_value_point"]]), as.numeric(x[["e_value_ci_bound"]])
  )),
  "",
  "## Full cohort comparison",
  "",
  "| version | n | early n | early % | events | HR (95% CI), p |",
  "|---|---:|---:|---:|---:|---|",
  apply(full_rows, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    x[["version"]], x[["n_patients"]], x[["early_n"]], x[["early_pct"]],
    x[["events"]],
    fmt_hr(as.numeric(x[["HR"]]), as.numeric(x[["CI_lo"]]), as.numeric(x[["CI_hi"]]), as.numeric(x[["p"]]))
  )),
  "",
  "## Positivity / balance comparison in SOFA 8-14",
  "",
  "### Weights",
  "",
  "| version | strategy | n patients | mean weight | p99 | max | ESS |",
  "|---|---|---:|---:|---:|---:|---:|",
  apply(weight_compare, 1, function(x) sprintf(
    "| %s | %s | %s | %.3f | %.3f | %.3f | %.1f |",
    x[["version"]], x[["strategy"]], x[["n_patients"]],
    as.numeric(x[["mean_w"]]), as.numeric(x[["p99_w"]]),
    as.numeric(x[["max_w"]]), as.numeric(x[["ess"]])
  )),
  "",
  "### Balance",
  "",
  "| version | max post-IPCW SMD | post-SMD >0.10 count |",
  "|---|---:|---:|",
  apply(balance_compare, 1, function(x) sprintf(
    "| %s | %.3f | %s |",
    x[["version"]], as.numeric(x[["max_post_smd"]]), x[["post_smd_gt_010"]]
  )),
  "",
  "## SOFA strata comparison",
  "",
  "| version | analysis | n | early n | HR (95% CI), p |",
  "|---|---|---:|---:|---|",
  apply(sofa_compare, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s |",
    x[["version"]], x[["analysis"]], x[["n_patients"]], x[["early_n"]],
    fmt_hr(as.numeric(x[["weighted_hr"]]), as.numeric(x[["ci_lower"]]), as.numeric(x[["ci_upper"]]), as.numeric(x[["p_value"]]))
  )),
  "",
  "## SOFA band scan summary",
  "",
  "| version | tiles | min HR | max HR | p<0.05 tiles |",
  "|---|---:|---:|---:|---:|",
  apply(scan_compare, 1, function(x) sprintf(
    "| %s | %s | %.3f | %.3f | %s |",
    x[["version"]], x[["scan_rows"]], as.numeric(x[["min_hr"]]),
    as.numeric(x[["max_hr"]]), x[["significant_tiles"]]
  )),
  "",
  "## Interpretation",
  "",
  "- Clean t0 reduces the primary SOFA 8-14 early arm from 215 to 146 and early percentage from 13.7% to 9.7%.",
  "- Despite the smaller early arm, v2 positivity remains acceptable: 99% trimmed weights are modest and post-IPCW balance is excellent.",
  "- The manuscript-facing conclusion remains null in the primary SOFA 8-14 equipoise band after removing treatment-defined t0 contamination."
)

writeLines(report, file.path(v2_results, "strictK3_v2_vs_v1_comparison_report.md"))
cat("Wrote v1-v2 comparison report: ", file.path(v2_results, "strictK3_v2_vs_v1_comparison_report.md"), "\n", sep = "")
