# Strict clean first-KDIGO3 v2: negative-control falsification and QBA.
# Scope: SOFA 8-14 primary equipoise cohort only. Reuses existing ccw_final
# clone-censor-weighted data and IPCW weights. Does not modify earlier outputs.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sandwich)
})

set.seed(20260716)

base_dir <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
data_dir <- file.path(base_dir, "data")
# CH-01 改动（仅路径）：输入改读本目录下的 CH-01 版 primary results；产物落本目录，不覆盖原 results/
source_dir <- Sys.getenv("QBA_SOURCE_DIR", file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "10_图件/out/fig07/qba_ch01"))
qba_only <- identical(Sys.getenv("QBA_ONLY", "0"), "1")
qba_mode <- Sys.getenv("QBA_MODE", "legacy")
stopifnot(qba_mode %in% c("legacy", "common_outcome"))
results_dir <- Sys.getenv("QBA_RESULTS_DIR", source_dir)
if (qba_only && normalizePath(results_dir, mustWork = FALSE) == normalizePath(source_dir)) {
  stop("QBA_ONLY requires an isolated results directory; production source is read-only")
}
# Extract ONLY the three summary functions, never source/evaluate the model script.
function_source <- Sys.getenv("QBA_FUNCTION_SOURCE", "11_docx回填/GITHUB_CODE/10_Analysis_Code/04_all_cause_28d_reanalysis/run20_NEW_v1_outcomes.R")
summary_function_names <- c("ch24_hr_to_rr", "ch24_rr_to_hr", "evalue_from_hr")
summary_expressions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  is.symbol(x[[2]]) && as.character(x[[2]]) %in% summary_function_names, parse(function_source))
stopifnot(length(summary_expressions) == 3L)
for (expr in summary_expressions) eval(expr)
figures_dir <- results_dir
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

safe_fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}


# QBA_ONLY never loads patient-level data or calls a fitting function.
if (!qba_only) {
ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    strategy_name = if_else(trt == 1L, "A_early", "B_deferred")
  )

eq_ccw <- ccw %>%
  filter(equipoise %in% TRUE) %>%
  mutate(
    nc_high_sofa12 = as.integer(sofa_t0 >= 12),
    nc_high_creat3 = as.integer(creatinine_t0 >= 3)
  )

eq_patient <- eq_ccw %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(early_crrt),
    strategy_name = if_else(trt == 1L, "observed_early", "observed_deferred")
  )

candidate_definitions <- tibble(
  candidate = c("nc_high_sofa12", "nc_high_creat3"),
  label = c("Pre-exposure high organ-failure burden", "Pre-exposure creatinine-dominant AKI severity"),
  operational_definition = c(
    "SOFA score at t0 at least 12, within the prespecified SOFA 8-14 equipoise cohort.",
    "Serum creatinine at t0 at least 3 mg/dL, before assignment to the early or deferred CRRT strategy."
  ),
  causal_path_argument = c(
    "CRRT timing after t0 cannot cause the SOFA value already observed at t0. The variable can still share confounding-by-indication because higher baseline organ failure may influence clinicians toward earlier CRRT and is prognostic for death.",
    "CRRT timing after t0 cannot cause the creatinine value already observed at t0. The variable can still share the renal-severity indication structure that drives early CRRT and mortality risk."
  ),
  limitation = c(
    "This is a pre-exposure falsification endpoint, not a post-treatment clinical negative-control outcome.",
    "This is a pre-exposure falsification endpoint, not a post-treatment clinical negative-control outcome."
  )
)


robust_poisson_rr <- function(dat, outcome, weight_col = NULL, analysis_label) {
  dat <- dat %>%
    filter(!is.na(.data[[outcome]]), !is.na(trt)) %>%
    mutate(
      y = .data[[outcome]],
      w = if (is.null(weight_col)) 1 else .data[[weight_col]]
    ) %>%
    filter(!is.na(w), is.finite(w), w > 0)

  if (nrow(dat) == 0 || length(unique(dat$trt)) < 2 || length(unique(dat$y)) < 2) {
    return(tibble(
      analysis = analysis_label,
      candidate = outcome,
      n_rows = nrow(dat),
      n_patients = n_distinct(dat$stay_id),
      events_early = NA_integer_,
      events_deferred = NA_integer_,
      risk_early = NA_real_,
      risk_deferred = NA_real_,
      risk_difference = NA_real_,
      rr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_
    ))
  }

  fit <- glm(y ~ trt, data = dat, family = poisson(link = "log"), weights = w)
  vc <- vcovCL(fit, cluster = dat$stay_id, type = "HC0")
  beta <- unname(coef(fit)["trt"])
  se <- sqrt(unname(vc["trt", "trt"]))
  if (!is.finite(se) || se < 1e-12) {
    se <- 0
  }
  z <- if (se == 0) NA_real_ else beta / se

  risk_summary <- dat %>%
    group_by(trt) %>%
    summarise(
      events = sum(y == 1, na.rm = TRUE),
      risk = weighted.mean(y, w, na.rm = TRUE),
      .groups = "drop"
    )
  risk_early <- risk_summary$risk[risk_summary$trt == 1]
  risk_deferred <- risk_summary$risk[risk_summary$trt == 0]
  events_early <- risk_summary$events[risk_summary$trt == 1]
  events_deferred <- risk_summary$events[risk_summary$trt == 0]

  tibble(
    analysis = analysis_label,
    candidate = outcome,
    n_rows = nrow(dat),
    n_patients = n_distinct(dat$stay_id),
    events_early = events_early,
    events_deferred = events_deferred,
    risk_early = risk_early,
    risk_deferred = risk_deferred,
    risk_difference = risk_early - risk_deferred,
    rr = exp(beta),
    ci_lower = exp(beta - 1.96 * se),
    ci_upper = exp(beta + 1.96 * se),
    p_value = if (is.na(z)) NA_real_ else 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

nc_results <- bind_rows(
  robust_poisson_rr(eq_patient, "nc_high_sofa12", NULL, "Observed actual early versus deferred, unweighted"),
  robust_poisson_rr(eq_patient, "nc_high_creat3", NULL, "Observed actual early versus deferred, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_sofa12", NULL, "Clone assigned strategies, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_creat3", NULL, "Clone assigned strategies, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_sofa12", "ipcw_w_trim", "Clone assigned strategies, IPCW weighted"),
  robust_poisson_rr(eq_ccw, "nc_high_creat3", "ipcw_w_trim", "Clone assigned strategies, IPCW weighted")
) %>%
  left_join(candidate_definitions %>% select(candidate, label), by = "candidate") %>%
  relocate(label, .after = candidate)

write_csv(candidate_definitions, file.path(results_dir, "strictK3_v2_negctrl_candidate_definitions.csv"))
write_csv(nc_results, file.path(results_dir, "strictK3_v2_negctrl_results.csv"))
} else {
  candidate_definitions <- read_csv(file.path(source_dir, "strictK3_v2_negctrl_candidate_definitions.csv"), show_col_types = FALSE)
  nc_results <- read_csv(file.path(source_dir, "strictK3_v2_negctrl_results.csv"), show_col_types = FALSE)
}


primary_results <- read_csv(
  file.path(source_dir, "strictK3_primary_v2_mortality_results.csv"),
  show_col_types = FALSE
)

primary <- primary_results %>%
  filter(analysis == "Primary SOFA 8-14 equipoise") %>%
  slice(1)

obs_hr <- primary$weighted_hr
obs_ci_lo <- primary$ci_lower
obs_ci_hi <- primary$ci_upper
e_value_point <- primary$e_value_point
e_value_ci_bound <- primary$e_value_ci_bound
if (qba_mode == "common_outcome") {
  # Authoritative CH01 HR/CI summaries; do not rerun Cox, weights, or outcomes.
  summary_source <- Sys.getenv("QBA_PRIMARY_SOURCE", "08_CH01/out/ch01_primary_dod.csv")
  all_cause_summary <- read_csv(summary_source, show_col_types = FALSE)
  stopifnot(setequal(all_cause_summary$population, c("Full cohort", "Equipoise SOFA 8-14")))
  summary_evalues <- bind_rows(lapply(seq_len(nrow(all_cause_summary)), function(i) {
    row <- all_cause_summary[i, ]
    ev <- evalue_from_hr(row$hr, row$ci_lower, row$ci_upper, mode = qba_mode)
    tibble(population = row$population, hr = row$hr, ci_lower = row$ci_lower,
           ci_upper = row$ci_upper, approximate_rr = ch24_hr_to_rr(row$hr),
           e_value_point = unname(ev[1]), e_value_ci_bound = unname(ev[2]),
           mode = qba_mode)
  }))
  write_csv(summary_evalues, file.path(results_dir, "ch24_q5_evalues_summary.csv"))
  primary <- summary_evalues %>% filter(population == "Equipoise SOFA 8-14")
  stopifnot(nrow(primary) == 1L)
  obs_hr <- primary$hr
  obs_ci_lo <- primary$ci_lower
  obs_ci_hi <- primary$ci_upper
  e_value_point <- primary$e_value_point
  e_value_ci_bound <- primary$e_value_ci_bound
}
obs_qba <- ch24_hr_to_rr(obs_hr, qba_mode)
target_qba <- ch24_hr_to_rr(c(1, 0.90, 1.10), qba_mode)


evalue_one <- function(rr) {
  ifelse(rr <= 1, 1, rr + sqrt(rr * (rr - 1)))
}

bias_factor <- function(rr_eu, rr_ud) {
  (rr_eu * rr_ud) / (rr_eu + rr_ud - 1)
}

equal_strength <- function(required_bf) {
  evalue_one(required_bf)
}

target_table <- tibble(
  target = c("Move point estimate to null HR 1.00", "Move point estimate to clinical benefit HR 0.90", "Move point estimate to clinical harm HR 1.10"),
  direction = c("upward confounding for observed harm", "upward confounding masking benefit", "downward confounding masking harm"),
  required_bias_factor = c(obs_qba / target_qba[1], obs_qba / target_qba[2], target_qba[3] / obs_qba),
  equal_rr_strength_for_exposure_and_outcome = equal_strength(required_bias_factor)
)

rr_grid <- seq(1, 3, by = 0.01)
qba_grid <- expand.grid(rr_exposure_confounder = rr_grid, rr_confounder_outcome = rr_grid) %>%
  as_tibble() %>%
  mutate(
    bias_factor = bias_factor(rr_exposure_confounder, rr_confounder_outcome),
    corrected_hr_if_upward_confounding = ch24_rr_to_hr(obs_qba / bias_factor, qba_mode),
    corrected_hr_if_downward_confounding = ch24_rr_to_hr(obs_qba * bias_factor, qba_mode),
    reaches_clinical_benefit_0_90 = corrected_hr_if_upward_confounding <= 0.90,
    reaches_clinical_harm_1_10 = corrected_hr_if_downward_confounding >= 1.10
  )

selected_rr <- c(1.10, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00)
solve_rr_ud <- function(required_bf, rr_eu) {
  ifelse(required_bf <= 1, 1, ifelse(rr_eu <= required_bf, Inf, required_bf * (rr_eu - 1) / (rr_eu - required_bf)))
}

tipping_array <- expand.grid(
  target = target_table$target,
  rr_exposure_confounder = selected_rr
) %>%
  as_tibble() %>%
  left_join(target_table %>% select(target, required_bias_factor), by = "target") %>%
  mutate(
    minimum_rr_confounder_outcome = solve_rr_ud(required_bias_factor, rr_exposure_confounder),
    minimum_rr_confounder_outcome = ifelse(
      is.infinite(minimum_rr_confounder_outcome),
      NA_real_,
      minimum_rr_confounder_outcome
    )
  )

scenario_table <- tibble(
  scenario = c("mild equal associations", "moderate equal associations", "stronger equal associations"),
  rr_exposure_confounder = c(1.25, 1.50, 1.75),
  rr_confounder_outcome = c(1.25, 1.50, 1.75)
) %>%
  mutate(
    bias_factor = bias_factor(rr_exposure_confounder, rr_confounder_outcome),
    corrected_hr_if_upward_confounding = ch24_rr_to_hr(obs_qba / bias_factor, qba_mode),
    corrected_hr_if_downward_confounding = ch24_rr_to_hr(obs_qba * bias_factor, qba_mode)
  )

write_csv(target_table, file.path(results_dir, "strictK3_v2_qba_target_thresholds.csv"))
write_csv(qba_grid, file.path(results_dir, "strictK3_v2_qba_grid.csv"))
write_csv(tipping_array, file.path(results_dir, "strictK3_v2_qba_tipping_array.csv"))
write_csv(scenario_table, file.path(results_dir, "strictK3_v2_qba_scenarios.csv"))

contour_breaks <- sort(unique(c(1, target_table$required_bias_factor, 2, 3)))

contour_plot <- ggplot(qba_grid, aes(rr_exposure_confounder, rr_confounder_outcome, z = bias_factor)) +
  geom_contour_filled(breaks = contour_breaks) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.35) +
  labs(
    title = "Tipping-point QBA for the SOFA 8-14 primary estimate",
    subtitle = if (qba_mode == "common_outcome") sprintf("Observed HR %.3f converted to approximate RR; contours are required bias factors", obs_hr) else sprintf("LEGACY: observed HR %.3f; contours are required bias factors", obs_hr),
    x = "RR: unmeasured confounder with early strategy",
    y = "RR: unmeasured confounder with 28-day mortality",
    fill = "Bias factor"
  ) +
  theme_minimal(base_size = 11)

ggsave(
  file.path(figures_dir, "strictK3_v2_qba_tipping_point_contour.png"),
  contour_plot,
  width = 7,
  height = 5,
  dpi = 300
)

fmt_rr <- function(x) {
  sprintf("%.3f (%.3f to %.3f)", x[["rr"]], x[["ci_lower"]], x[["ci_upper"]])
}

nc_table_lines <- nc_results %>%
  mutate(
    risk_early_text = safe_fmt(risk_early),
    risk_deferred_text = safe_fmt(risk_deferred),
    risk_difference_text = safe_fmt(risk_difference),
    rr_ci = sprintf("%.3f (%.3f to %.3f)", rr, ci_lower, ci_upper),
    p_text = if_else(is.na(p_value), "NA", if_else(p_value < 0.001, "less than 0.001", sprintf("%.3f", p_value)))
  ) %>%
  select(label, analysis, n_patients, risk_early_text, risk_deferred_text, risk_difference_text, rr_ci, p_text) %>%
  apply(1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    x[["label"]], x[["analysis"]], x[["n_patients"]],
    x[["risk_early_text"]], x[["risk_deferred_text"]], x[["risk_difference_text"]],
    x[["rr_ci"]], x[["p_text"]]
  ))

target_lines <- target_table %>%
  mutate(
    required_bias_factor = safe_fmt(required_bias_factor),
    equal_rr_strength_for_exposure_and_outcome = safe_fmt(equal_rr_strength_for_exposure_and_outcome)
  ) %>%
  apply(1, function(x) sprintf(
    "| %s | %s | %s | %s |",
    x[["target"]], x[["direction"]], x[["required_bias_factor"]],
    x[["equal_rr_strength_for_exposure_and_outcome"]]
  ))

scenario_lines <- scenario_table %>%
  mutate(
    bias_factor = safe_fmt(bias_factor),
    corrected_hr_if_upward_confounding = safe_fmt(corrected_hr_if_upward_confounding),
    corrected_hr_if_downward_confounding = safe_fmt(corrected_hr_if_downward_confounding)
  ) %>%
  apply(1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    x[["scenario"]], x[["rr_exposure_confounder"]], x[["rr_confounder_outcome"]],
    x[["bias_factor"]], x[["corrected_hr_if_upward_confounding"]],
    x[["corrected_hr_if_downward_confounding"]]
  ))

report <- c(
  "# strictK3 v2 negative-control falsification and QBA report",
  "",
  "## ASSUMPTIONS",
  "",
  "- Analysis is restricted to the prespecified SOFA 8-14 primary equipoise cohort.",
  if (qba_only) "- QBA_ONLY: existing negative-control CSVs were loaded unchanged; no patient-level data or model fitting was used." else "- The script reuses the existing clone-censor-weighted dataset and trimmed IPCW weights.",
  "- No clean post-treatment clinical negative-control outcome was available in the existing derived analysis files. Therefore, Task 1 is implemented as a pre-exposure falsification analysis rather than as a conventional post-treatment negative-control outcome analysis.",
  "- The two falsification endpoints are measured at t0, before assignment to the early or deferred CRRT strategy in the emulated trial. CRRT timing after t0 cannot cause these t0 measurements, but both variables plausibly share the confounding-by-indication structure.",
  "- The observed actual early-versus-deferred row is included to show that each candidate tracks clinical indication. The clone assigned-strategy rows are the actual CCW/IPCW falsification checks.",
  "- For clone assigned-strategy falsification rows, p values are secondary because each pre-exposure endpoint is duplicated across strategy clones; absolute risk differences and RR magnitude are more informative than null-hypothesis testing.",
  if (qba_mode == "common_outcome") "- QBA applies RR bias factors after the common-outcome HR-to-approximate-RR conversion. Corrected HR fields are numerical inverse-mapped approximate HR equivalents, not refitted or identified causal HRs. Targets are separately converted. These sensitivity parameters do not establish that such a confounder exists." else "- LEGACY REGRESSION ONLY: direct HR-as-RR calculations are reproduced without claiming common-outcome validity.",
  "",
  "## Variables and source tables",
  "",
  "- `ccw_final.rds`: `stay_id`, `strategy`, `early_crrt`, `equipoise`, `follow_days`, `event`, `ipcw_w_trim`, `sofa_t0`, `creatinine_t0`.",
  if (qba_mode == "common_outcome") "- Primary HR/CI source: 08_CH01/out/ch01_primary_dod.csv; no model rerun." else "- LEGACY primary estimate source: production mortality-summary CSV, unchanged.",
  "",
  "## Task 1. Candidate definitions and rationale",
  "",
  "| candidate | operational definition | no-causal-path argument | limitation |",
  "|---|---|---|---|",
  apply(candidate_definitions, 1, function(x) sprintf(
    "| %s | %s | %s | %s |",
    x[["label"]], x[["operational_definition"]],
    x[["causal_path_argument"]], x[["limitation"]]
  )),
  "",
  "## Task 1. Results",
  "",
  "| candidate | analysis | n patients | early risk | deferred risk | risk difference | RR (95 percent CI) | p |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  nc_table_lines,
  "",
  "Interpretation for manuscript:",
  "",
  "The observed actual treatment groups showed clear indication-related imbalance in pre-exposure severity markers. After CCW/IPCW, assigned-strategy absolute differences were small, although statistically detectable for the creatinine marker because of the cloned design and high precision. This supports, but cannot prove, that the primary IPCW analysis materially reduced measured severity-driven confounding; residual unmeasured clinician-judgement confounding remains possible.",
  "",
  "## Task 2. Tipping-point QBA",
  "",
  sprintf(
    "Primary SOFA 8-14 estimate: HR %.3f with 95 percent CI %.3f to %.3f; E-value point %.3f and CI-bound %.3f.",
    obs_hr, obs_ci_lo, obs_ci_hi, e_value_point, e_value_ci_bound
  ),
  "",
  "| target | confounding direction | required bias factor | equal RR strength for exposure and outcome |",
  "|---|---|---:|---:|",
  target_lines,
  "",
  "Scenario bounds:",
  "",
  "| scenario | RR exposure-confounder | RR confounder-outcome | bias factor | approximate HR equivalent if upward confounding | approximate HR equivalent if downward confounding |",
  "|---|---:|---:|---:|---:|---:|",
  scenario_lines,
  "",
  "Interpretation for manuscript:",
  "",
  "The point-estimate E-value was small and the confidence-interval-bound E-value was 1 because the primary CI crossed the null. A single unmeasured confounder with modest-to-moderate joint associations with early CRRT selection and 28-day mortality could shift the point estimate across clinically important benefit or harm thresholds. This should be reported as a limitation rather than downplayed.",
  "",
  "## Output files",
  "",
  if (qba_only) paste0("- Existing read-only negative-control definitions: ", file.path(source_dir, "strictK3_v2_negctrl_candidate_definitions.csv")) else "- `results/strictK3_v2_negctrl_candidate_definitions.csv`",
  if (qba_only) paste0("- Existing read-only negative-control results: ", file.path(source_dir, "strictK3_v2_negctrl_results.csv")) else "- `results/strictK3_v2_negctrl_results.csv`",
  "- `results/strictK3_v2_qba_target_thresholds.csv`",
  "- `results/strictK3_v2_qba_grid.csv`",
  "- `results/strictK3_v2_qba_tipping_array.csv`",
  "- `results/strictK3_v2_qba_scenarios.csv`",
  "- `figures/strictK3_v2_qba_tipping_point_contour.png`"
)

writeLines(report, file.path(results_dir, "strictK3_v2_negctrl_qba_report.md"))

message("Wrote strictK3_v2_negctrl_qba_report.md and supporting CSV/figure outputs.")

if (qba_mode == "common_outcome") {
  write_csv(tibble(mode = qba_mode, observed_hr = obs_hr, observed_approximate_rr = obs_qba,
                   target_hr = c(1, 0.9, 1.1), target_approximate_rr = target_qba,
                   corrected_hr_fields = "inverse-mapped approximate HR equivalents"),
            file.path(results_dir, "ch24_q5_scale_provenance.csv"))
}
