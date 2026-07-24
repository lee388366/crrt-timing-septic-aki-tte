# =============================================================================
# Strict clean first-KDIGO3 v2: physiologic KDIGO3 path sensitivity
#
# Adds:
#   1) SOFA 8-14 stratified HR by clean KDIGO3 path:
#      creatinine_only / urine_output_only / both_creatinine_and_uo
#   2) IPCW sensitivity with stage3_path included in numerator and denominator
#      weight models. Does not overwrite ccw_final.rds.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmnet)
  library(survival)
  library(tidyr)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260716)

baseline <- readRDS(file.path(data_dir, "baseline.rds")) %>%
  mutate(
    stage3_path = factor(stage3_path, levels = c(
      "urine_output_only", "creatinine_only", "both_creatinine_and_uo"
    ))
  )

ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    stage3_path = factor(stage3_path, levels = c(
      "urine_output_only", "creatinine_only", "both_creatinine_and_uo"
    ))
  )

tv_hourly <- readRDS(file.path(data_dir, "tv_hourly.rds"))
tv_hour23 <- readRDS(file.path(data_dir, "tv_hour23.rds"))

L0_vars <- c("age", "sofa_t0", "apsiii", "creatinine_baseline",
             "nee_at_t0", "lactate_t0", "uo_rate_t0", "map_t0", "cci")
Lk_vars <- c("map_h", "nee_h", "lactate_h", "sofa_h",
             "kdigo_h", "creatinine_h", "creatinine_delta", "fluid_balance_h")

impute_median <- function(df, vars) {
  for (v in vars) {
    med <- median(df[[v]], na.rm = TRUE)
    df[[v]][is.na(df[[v]])] <- med
  }
  df
}

run_cox <- function(dat, label) {
  dat <- dat %>%
    filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim), !is.na(trt))
  if (n_distinct(dat$stay_id) < 50 || length(unique(dat$trt)) < 2) {
    return(tibble(
      analysis = label,
      n_patients = n_distinct(dat$stay_id),
      clone_rows = nrow(dat),
      early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
      events = sum(dat$event, na.rm = TRUE),
      HR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_, p = NA_real_
    ))
  }
  fit <- coxph(
    Surv(follow_days, event) ~ trt,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  s <- summary(fit)
  tibble(
    analysis = label,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat),
    early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
    events = sum(dat$event, na.rm = TRUE),
    HR = unname(s$coefficients["trt", "exp(coef)"]),
    CI_lo = unname(s$conf.int["trt", "lower .95"]),
    CI_hi = unname(s$conf.int["trt", "upper .95"]),
    p = unname(s$coefficients["trt", "Pr(>|z|)"])
  )
}

path_counts <- baseline %>%
  filter(sofa_equipoise %in% TRUE) %>%
  count(stage3_path, early_crrt, name = "n") %>%
  group_by(stage3_path) %>%
  mutate(path_n = sum(n), pct_within_path = round(100 * n / path_n, 1)) %>%
  ungroup()

path_summary <- baseline %>%
  filter(sofa_equipoise %in% TRUE) %>%
  group_by(stage3_path) %>%
  summarise(
    n = n(),
    early_n = sum(early_crrt, na.rm = TRUE),
    early_pct = round(100 * mean(early_crrt, na.rm = TRUE), 1),
    death28_pct = round(100 * mean(died_28d, na.rm = TRUE), 1),
    .groups = "drop"
  )

path_hr <- bind_rows(lapply(levels(baseline$stage3_path), function(pth) {
  run_cox(
    ccw %>% filter(equipoise %in% TRUE, stage3_path == pth),
    paste0("SOFA 8-14 path: ", pth)
  )
}))

path_weight_summary <- ccw %>%
  filter(equipoise %in% TRUE) %>%
  group_by(stage3_path, strategy) %>%
  summarise(
    n_patients = n_distinct(stay_id),
    mean_w = mean(ipcw_w_trim, na.rm = TRUE),
    p99_w = unname(quantile(ipcw_w_trim, 0.99, na.rm = TRUE)),
    max_w = max(ipcw_w_trim, na.rm = TRUE),
    ess = sum(ipcw_w_trim, na.rm = TRUE)^2 / sum(ipcw_w_trim^2, na.rm = TRUE),
    .groups = "drop"
  )

# Refit IPCW with stage3_path included ----------------------------------------

cohort <- baseline %>%
  group_by(subject_id) %>%
  slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(sofa_t0)) %>%
  impute_median(L0_vars)

GRACE <- 24

ccw_base <- cohort %>%
  crossing(strategy = c("A_early", "B_deferred")) %>%
  mutate(
    clone_id = paste(stay_id, strategy, sep = "_"),
    h_crrt = crrt_hours_from_t0,
    protocol_censor_h = case_when(
      strategy == "A_early" & (is.na(h_crrt) | h_crrt > GRACE) ~ as.double(GRACE),
      strategy == "B_deferred" & !is.na(h_crrt) & h_crrt <= GRACE ~ h_crrt,
      TRUE ~ NA_real_
    ),
    protocol_censored = !is.na(protocol_censor_h),
    follow_h = pmin(follow_time_28d * 24, coalesce(protocol_censor_h, follow_time_28d * 24)),
    follow_days = follow_h / 24,
    event = died_28d & (is.na(protocol_censor_h) | follow_time_28d * 24 < protocol_censor_h)
  )

tv_hour23_uniq <- tv_hour23 %>% distinct(stay_id, .keep_all = TRUE)

A_data <- ccw_base %>%
  filter(strategy == "A_early") %>%
  left_join(tv_hour23_uniq, by = "stay_id") %>%
  mutate(y_A = as.integer(early_crrt)) %>%
  impute_median(Lk_vars)

X_A_denom <- model.matrix(
  ~ age + sofa_t0 + apsiii + creatinine_baseline + nee_at_t0 +
    lactate_t0 + uo_rate_t0 + map_t0 + cci + stage3_path +
    map_h + nee_h + lactate_h + sofa_h + kdigo_h +
    creatinine_h + creatinine_delta + fluid_balance_h,
  data = A_data
)[, -1]

cv_A <- cv.glmnet(X_A_denom, A_data$y_A, family = "binomial", alpha = 0,
                  nfolds = 10, type.measure = "deviance")
p_denom_A <- as.vector(predict(cv_A, X_A_denom, s = "lambda.1se", type = "response"))

lm_num_A <- glm(
  y_A ~ age + sofa_t0 + apsiii + creatinine_baseline +
    nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci + stage3_path,
  data = A_data,
  family = binomial()
)
p_num_A <- predict(lm_num_A, type = "response")

A_data <- A_data %>%
  mutate(
    p_denom = p_denom_A,
    p_num = p_num_A,
    ipcw_w = if_else(
      early_crrt,
      p_num / pmax(p_denom, 0.01),
      (1 - p_num) / pmax(1 - p_denom, 0.01)
    )
  )

B_base <- ccw_base %>% filter(strategy == "B_deferred")

B_long <- B_base %>%
  select(stay_id, early_crrt, h_crrt, stage3_path, all_of(L0_vars)) %>%
  left_join(tv_hourly, by = "stay_id", relationship = "many-to-many") %>%
  mutate(
    crrt_this_hour = !is.na(h_crrt) & h_crrt >= hour_index &
      h_crrt < (hour_index + 1) & h_crrt <= GRACE,
    at_risk = is.na(h_crrt) | h_crrt > hour_index
  ) %>%
  filter(at_risk) %>%
  impute_median(Lk_vars)

X_B_denom <- model.matrix(
  ~ age + sofa_t0 + apsiii + creatinine_baseline + nee_at_t0 +
    lactate_t0 + uo_rate_t0 + map_t0 + cci + stage3_path + hour_index +
    map_h + nee_h + lactate_h + sofa_h + kdigo_h +
    creatinine_h + creatinine_delta + fluid_balance_h,
  data = B_long
)[, -1]
y_B <- as.integer(B_long$crrt_this_hour)

cv_B <- cv.glmnet(X_B_denom, y_B, family = "binomial", alpha = 0,
                  nfolds = 10, type.measure = "deviance")
p_denom_B <- as.vector(predict(cv_B, X_B_denom, s = "lambda.1se", type = "response"))

lm_num_B <- glm(
  crrt_this_hour ~ age + sofa_t0 + apsiii + creatinine_baseline +
    nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci + stage3_path + hour_index,
  data = B_long,
  family = binomial()
)
p_num_B <- predict(lm_num_B, type = "response")

B_long <- B_long %>%
  mutate(
    p_denom_h = p_denom_B,
    p_num_h = p_num_B,
    w_contrib = (1 - p_num_h) / pmax(1 - p_denom_h, 0.01)
  )

B_weights <- B_long %>%
  group_by(stay_id) %>%
  summarise(ipcw_w = prod(w_contrib, na.rm = TRUE), .groups = "drop")

w99_A <- quantile(A_data$ipcw_w, 0.99, na.rm = TRUE)
w99_B <- quantile(B_weights$ipcw_w, 0.99, na.rm = TRUE)

A_weights <- A_data %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  transmute(stay_id, ipcw_w, ipcw_w_trim = pmin(ipcw_w, w99_A))
B_weights <- B_weights %>%
  mutate(ipcw_w_trim = pmin(ipcw_w, w99_B))

ccw_A <- ccw_base %>%
  filter(strategy == "A_early") %>%
  left_join(A_weights, by = "stay_id", relationship = "many-to-one")
ccw_B <- ccw_base %>%
  filter(strategy == "B_deferred") %>%
  left_join(B_weights, by = "stay_id", relationship = "many-to-one") %>%
  mutate(ipcw_w = coalesce(ipcw_w, 1), ipcw_w_trim = coalesce(ipcw_w_trim, 1))

ccw_path <- bind_rows(ccw_A, ccw_B) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    equipoise = !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14
  )

saveRDS(ccw_path, file.path(data_dir, "ccw_final_path_weighted.rds"))

path_weight_model_results <- bind_rows(
  run_cox(ccw_path, "Full cohort: IPCW includes stage3_path"),
  run_cox(ccw_path %>% filter(equipoise %in% TRUE), "SOFA 8-14: IPCW includes stage3_path")
)

path_weight_model_summary <- ccw_path %>%
  group_by(strategy, equipoise) %>%
  summarise(
    n_patients = n_distinct(stay_id),
    mean_w = mean(ipcw_w_trim, na.rm = TRUE),
    p99_w = unname(quantile(ipcw_w_trim, 0.99, na.rm = TRUE)),
    max_w = max(ipcw_w_trim, na.rm = TRUE),
    ess = sum(ipcw_w_trim, na.rm = TRUE)^2 / sum(ipcw_w_trim^2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(path_summary, file.path(results_dir, "strictK3_v2_path_summary_sofa8_14.csv"))
write_csv(path_counts, file.path(results_dir, "strictK3_v2_path_counts_sofa8_14.csv"))
write_csv(path_hr, file.path(results_dir, "strictK3_v2_path_stratified_hr.csv"))
write_csv(path_weight_summary, file.path(results_dir, "strictK3_v2_path_weight_summary_current_ipcw.csv"))
write_csv(path_weight_model_results, file.path(results_dir, "strictK3_v2_path_adjusted_ipcw_results.csv"))
write_csv(path_weight_model_summary, file.path(results_dir, "strictK3_v2_path_adjusted_ipcw_weight_summary.csv"))

fmt_hr <- function(hr, lo, hi, p) {
  ifelse(
    is.na(hr),
    "NA",
    sprintf("%.3f (%.3f to %.3f), p=%s", hr, lo, hi, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  )
}

report <- c(
  "# Strict clean KDIGO3 v2: AKI staging-path sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Assumptions",
  "",
  "- `stage3_path` is defined at clean physiologic KDIGO3 t0 only: `creatinine_only`, `urine_output_only`, or `both_creatinine_and_uo`.",
  "- Main stratified HRs use the primary v2 IPCW weights.",
  "- The path-adjusted IPCW sensitivity adds `stage3_path` to both numerator and denominator models, preserving the stabilized-weight structure.",
  "",
  "## SOFA 8-14 path distribution",
  "",
  "| stage3 path | n | early n | early % | 28d death % |",
  "|---|---:|---:|---:|---:|",
  apply(path_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s |",
    x[["stage3_path"]], x[["n"]], x[["early_n"]], x[["early_pct"]], x[["death28_pct"]]
  )),
  "",
  "## SOFA 8-14 stratified weighted Cox by path",
  "",
  "| analysis | n | early n | events | HR (95% CI), p |",
  "|---|---:|---:|---:|---:|",
  apply(path_hr, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s |",
    x[["analysis"]], x[["n_patients"]], x[["early_n"]], x[["events"]],
    fmt_hr(as.numeric(x[["HR"]]), as.numeric(x[["CI_lo"]]), as.numeric(x[["CI_hi"]]), as.numeric(x[["p"]]))
  )),
  "",
  "## IPCW sensitivity including stage3_path",
  "",
  sprintf("- Weight trim caps with stage3_path: A %.3f; B %.3f.", w99_A, w99_B),
  "",
  "| analysis | n | early n | events | HR (95% CI), p |",
  "|---|---:|---:|---:|---:|",
  apply(path_weight_model_results, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s |",
    x[["analysis"]], x[["n_patients"]], x[["early_n"]], x[["events"]],
    fmt_hr(as.numeric(x[["HR"]]), as.numeric(x[["CI_lo"]]), as.numeric(x[["CI_hi"]]), as.numeric(x[["p"]]))
  )),
  "",
  "## Path-adjusted IPCW weight summary",
  "",
  "| strategy | equipoise | n patients | mean weight | p99 | max | ESS |",
  "|---|---:|---:|---:|---:|---:|---:|",
  apply(path_weight_model_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %.3f | %.3f | %.3f | %.1f |",
    x[["strategy"]], x[["equipoise"]], x[["n_patients"]],
    as.numeric(x[["mean_w"]]), as.numeric(x[["p99_w"]]),
    as.numeric(x[["max_w"]]), as.numeric(x[["ess"]])
  ))
)

writeLines(report, file.path(results_dir, "strictK3_v2_path_sensitivity_report.md"))
print(path_summary)
print(path_hr)
print(path_weight_model_results)
