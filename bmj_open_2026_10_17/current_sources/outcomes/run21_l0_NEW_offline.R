# =============================================================================
# Strict clean first-KDIGO3: IPCW denominator L0-only sensitivity
# Mirrors the existing project sensitivity in 22_fix_c1_c3_sensitivity.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(DBI)
  library(RPostgres)
})

base_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_inputs")
results_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_l0/NEW/results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

get_cox_res <- function(fit) {
  s <- summary(fit)
  coef_name <- rownames(s$coefficients)[1]
  tibble(
    HR = unname(s$coefficients[coef_name, "exp(coef)"]),
    CI_lo = unname(s$conf.int[coef_name, "lower .95"]),
    CI_hi = unname(s$conf.int[coef_name, "upper .95"]),
    p = unname(s$coefficients[coef_name, "Pr(>|z|)"])
  )
}

run_cox <- function(dat, label, weight_col = "ipcw_w_trim") {
  fit <- coxph(
    as.formula(paste0("Surv(follow_days, event) ~ strategy_f")),
    data = dat,
    weights = dat[[weight_col]],
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  get_cox_res(fit) %>%
    mutate(
      analysis = label,
      n_patients = n_distinct(dat$stay_id),
      clone_rows = nrow(dat),
      early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
      events = sum(dat$event, na.rm = TRUE),
      .before = HR
    )
}

impute_median <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df)) {
      med <- median(df[[v]], na.rm = TRUE)
      df[[v]][is.na(df[[v]])] <- med
    }
  }
  df
}

baseline <- readRDS(file.path(data_dir, "baseline.rds")) %>%
  group_by(subject_id) %>%
  slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(sofa_t0))

# ===== RUN-21 · CH-01 最小改动：只重定义结局两列（died_28d / follow_time_28d）=====
# 本脚本模型消费 follow_days/event（L83-85 由两列导出）⇒ 必须同式重导（run14 教训）。
.dod <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_l0/dod_current_CH01_unique.rds"))
.dod$dod <- as.Date(.dod$dod)
ch01_outcome <- function(df) {
  df %>% left_join(.dod, by = "subject_id") %>%
    mutate(.t0 = as.POSIXct(t0, tz = "UTC"), .dt = as.POSIXct(deathtime, tz = "UTC"),
      .dc = coalesce(as.numeric(difftime(.dt, .t0, units = "days")), as.numeric(dod - as.Date(.t0))),
      died_28d = !is.na(.dc) & .dc >= 0 & .dc <= 28,
      follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))) %>%
    select(-.dc, -.t0, -.dt, -dod)
}
ch01_derive <- function(df) {
  df %>% mutate(
    follow_h = pmin(follow_time_28d * 24, coalesce(protocol_censor_h, follow_time_28d * 24)),
    follow_days = follow_h / 24,
    event = died_28d & (is.na(protocol_censor_h) | follow_time_28d * 24 < protocol_censor_h))
}

L0_vars <- c(
  "age", "sofa_t0", "apsiii", "creatinine_baseline",
  "nee_at_t0", "lactate_t0", "uo_rate_t0", "map_t0", "cci"
)
baseline <- impute_median(baseline, L0_vars)
baseline <- ch01_outcome(baseline)   # RUN-21: 只重定义两列（下游 L83-85 自动重导）
cat(sprintf("[RUN-21 check] baseline died_28d = %d / %d\n", sum(baseline$died_28d), nrow(baseline)))

GRACE <- 24
ccw_sens <- baseline %>%
  tidyr::crossing(strategy = c("A_early", "B_deferred")) %>%
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

A_sens <- ccw_sens %>%
  filter(strategy == "A_early") %>%
  mutate(y_A = as.integer(early_crrt))

denom_A <- glm(
  y_A ~ age + sofa_t0 + apsiii + creatinine_baseline +
    nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci,
  data = A_sens,
  family = binomial()
)
p_denom_A <- predict(denom_A, type = "response")
num_A <- glm(y_A ~ 1, data = A_sens, family = binomial())
p_num_A <- predict(num_A, type = "response")

A_sens <- A_sens %>%
  mutate(
    p_denom = p_denom_A,
    p_num = p_num_A,
    ipcw_w = if_else(
      early_crrt,
      p_num / pmax(p_denom, 0.01),
      (1 - p_num) / pmax(1 - p_denom, 0.01)
    )
  )

w99_A <- quantile(A_sens$ipcw_w, 0.99, na.rm = TRUE)
A_sens <- A_sens %>% mutate(ipcw_w_trim = pmin(ipcw_w, w99_A))

ccw_orig <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)
ccw_orig <- ch01_derive(ch01_outcome(ccw_orig))   # RUN-21: 两列 + 同式重导 follow_days/event
cat(sprintf("[RUN-21 check] ccw_orig died_28d (unique stays) = %d\n", length(unique(ccw_orig$stay_id[ccw_orig$died_28d %in% TRUE]))))

B_orig_weights <- ccw_orig %>%
  filter(strategy == "B_deferred") %>%
  select(stay_id, ipcw_w_orig = ipcw_w, ipcw_w_trim_orig = ipcw_w_trim)

ccw_B_sens <- ccw_sens %>%
  filter(strategy == "B_deferred") %>%
  left_join(B_orig_weights, by = "stay_id") %>%
  mutate(
    ipcw_w = coalesce(ipcw_w_orig, 1),
    ipcw_w_trim = coalesce(ipcw_w_trim_orig, 1)
  )

ccw_l0 <- bind_rows(
  A_sens %>% select(any_of(names(ccw_B_sens))),
  ccw_B_sens %>% select(any_of(names(A_sens)))
) %>%
  mutate(
    strategy_f = factor(strategy, levels = c("B_deferred", "A_early")),
    equipoise = !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14
  )

ccw_orig2 <- ccw_orig %>%
  mutate(
    strategy_f = factor(strategy, levels = c("B_deferred", "A_early")),
    equipoise = !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14
  )

results <- bind_rows(
  run_cox(ccw_orig2, "Original denominator: full cohort"),
  run_cox(filter(ccw_orig2, equipoise), "Original denominator: SOFA 8-14"),
  run_cox(ccw_l0, "L0-only denominator: full cohort"),
  run_cox(filter(ccw_l0, equipoise), "L0-only denominator: SOFA 8-14")
) %>%
  mutate(
    HR = round(HR, 3),
    CI_lo = round(CI_lo, 3),
    CI_hi = round(CI_hi, 3),
    p = round(p, 4)
  )

weight_summary <- bind_rows(
  ccw_orig2 %>%
    group_by(strategy) %>%
    summarise(version = "original", mean_w = mean(ipcw_w_trim), median_w = median(ipcw_w_trim), max_w = max(ipcw_w_trim), .groups = "drop"),
  ccw_l0 %>%
    group_by(strategy) %>%
    summarise(version = "L0_only_A", mean_w = mean(ipcw_w_trim), median_w = median(ipcw_w_trim), max_w = max(ipcw_w_trim), .groups = "drop")
)

write_csv(results, file.path(results_dir, "strictK3_v2_l0_only_sensitivity_results.csv"))
write_csv(weight_summary, file.path(results_dir, "strictK3_v2_l0_only_weight_summary.csv"))
saveRDS(list(results = results, weights = weight_summary, ccw_l0 = ccw_l0), file.path(results_dir, "strictK3_v2_l0_only_sensitivity.rds"))

report <- c(
  "# Strict clean first-KDIGO3 L0-only IPCW denominator sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Scope: mirrors the existing project sensitivity in `22_fix_c1_c3_sensitivity.R`: Strategy A denominator is re-estimated using L0 variables only; Strategy B weights are retained from the primary pooled hourly model.",
  "",
  "## Cox results",
  "",
  paste(capture.output(print(results, n = nrow(results))), collapse = "\n"),
  "",
  "## Weight summary",
  "",
  paste(capture.output(print(weight_summary, n = nrow(weight_summary))), collapse = "\n"),
  "",
  "Outputs:",
  "",
  paste0("- `", file.path(results_dir, "strictK3_v2_l0_only_sensitivity_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_l0_only_weight_summary.csv"), "`")
)
writeLines(report, file.path(results_dir, "strictK3_v2_l0_only_sensitivity_report.md"))

cat("\n=== strictK3_v2 L0-only sensitivity ===\n")
print(results, n = nrow(results))
cat("\nSaved report:", file.path(results_dir, "strictK3_v2_l0_only_sensitivity_report.md"), "\n")