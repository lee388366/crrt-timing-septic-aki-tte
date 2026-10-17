# =============================================================================
# Strict clean first-KDIGO3: rescue-exclusion sensitivity models
#
# Variants:
#   1. no_rescue_exclusion
#   2. primary_R1_R2_R2b_R3 (reference, already run)
#   3. sensitivity_R1_R2_R2b_R3_R4
#   4. sensitivity_R1_confirm_R2_R2b_R3
#
# Each non-reference variant gets its own data directory. The legacy 03a/03 CCW
# scripts are evaluated with path redirection and tie fixes, so primary outputs
# are not overwritten.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(fs)
  library(survival)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
strict_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run17/NEW/rescue/strictK3_primary_v2")
primary_data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimiciv31",
  host = "localhost",
  port = 5432,
  user = "postgres",
  password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

# ===== RUN-17 · CH-01 最小改动：只重定义结局两列（died_28d / follow_time_28d）=====
# 03_ccw_ipcw.R L61-66 由这两列导出 follow_days/event ⇒ 本补丁落在 build_baseline 内即可贯通口径。
ch01_dod <- dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients") %>% mutate(dod = as.Date(dod))
ch01_recode <- function(df) {
  df %>%
    left_join(ch01_dod, by = "subject_id") %>%
    mutate(
      .t0 = as.POSIXct(t0, tz = "UTC"),
      .dt = as.POSIXct(deathtime, tz = "UTC"),
      .dc = coalesce(as.numeric(difftime(.dt, .t0, units = "days")), as.numeric(dod - as.Date(.t0))),
      died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
      follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
    ) %>%
    select(-.dc, -.t0, -.dt, -dod)
}

pre_rescue <- readRDS(file.path(primary_data_dir, "cohort_pre_rescue_all.rds"))
primary_cohort <- readRDS(file.path(primary_data_dir, "cohort_base.rds"))
cat("[RUN-17 check] primary cohort CH-01 died_28d:", sum(ch01_recode(primary_cohort)$died_28d, na.rm = TRUE), " (expect 1248)\n")

threshold_flags <- dbGetQuery(con, "
  SELECT
    pr.stay_id,
    coalesce(k.k_ge_6, false) AS r1_k_ge_6,
    coalesce(k.k_ge_7, false) AS r1_k_ge_7,
    coalesce(ph.ph_le_720, false) AS r2_ph_le_720,
    coalesce(ph.ph_le_710, false) AS r2_ph_le_710,
    coalesce(hco3.hco3_le_15, false) AS r2b_hco3_le_15,
    coalesce(hco3.hco3_le_10, false) AS r2b_hco3_le_10,
    coalesce(bun.bun_ge_90, false) AS r3_bun_ge_90,
    coalesce(bun.bun_ge_112, false) AS r3_bun_ge_112
  FROM data_extract_crrt_strictk3_v2.cohort_pre_rescue_all pr
  LEFT JOIN LATERAL (
    SELECT
      bool_or(c.potassium >= 6.0) AS k_ge_6,
      bool_or(c.potassium >= 7.0) AS k_ge_7
    FROM mimiciv_derived.chemistry c
    WHERE c.hadm_id = pr.hadm_id
      AND c.charttime > pr.t0 - interval '6 hours'
      AND c.charttime <= pr.t0
      AND c.potassium IS NOT NULL
  ) k ON true
  LEFT JOIN LATERAL (
    SELECT
      bool_or(bg.ph <= 7.20) AS ph_le_720,
      bool_or(bg.ph <= 7.10) AS ph_le_710
    FROM mimiciv_derived.bg bg
    WHERE bg.hadm_id = pr.hadm_id
      AND bg.charttime > pr.t0 - interval '6 hours'
      AND bg.charttime <= pr.t0
      AND bg.ph IS NOT NULL
      AND bg.specimen = 'ART.'
  ) ph ON true
  LEFT JOIN LATERAL (
    SELECT
      bool_or(c.bicarbonate <= 15) AS hco3_le_15,
      bool_or(c.bicarbonate <= 10) AS hco3_le_10
    FROM mimiciv_derived.chemistry c
    WHERE c.hadm_id = pr.hadm_id
      AND c.charttime > pr.t0 - interval '6 hours'
      AND c.charttime <= pr.t0
      AND c.bicarbonate IS NOT NULL
  ) hco3 ON true
  LEFT JOIN LATERAL (
    SELECT
      bool_or(c.bun >= 90) AS bun_ge_90,
      bool_or(c.bun >= 112) AS bun_ge_112
    FROM mimiciv_derived.chemistry c
    WHERE c.hadm_id = pr.hadm_id
      AND c.charttime > pr.t0 - interval '24 hours'
      AND c.charttime <= pr.t0
      AND c.bun IS NOT NULL
  ) bun ON true
")

pre_rescue_ext <- pre_rescue %>%
  left_join(threshold_flags, by = "stay_id") %>%
  mutate(
    rescue_threshold_liberal = r1_k_ge_6 | r2_ph_le_720 | r2b_hco3_le_15 | r3_bun_ge_90,
    rescue_threshold_conservative = r1_k_ge_7 | r2_ph_le_710 | r2b_hco3_le_10 | r3_bun_ge_112
  )

impute_median <- function(df, vars) {
  for (v in vars) {
    med <- median(df[[v]], na.rm = TRUE)
    df[[v]][is.na(df[[v]])] <- med
  }
  df
}

build_baseline <- function(cohort_raw, data_dir) {
  cohort <- cohort_raw %>%
    group_by(subject_id) %>%
    slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
    ungroup()

  # RUN-17：只重定义结局两列（下游 03_ccw_ipcw.R 重导 follow_days/event，口径自动贯通）
  cohort <- ch01_recode(cohort)

  sofa_hourly <- dbGetQuery(con, "
    SELECT stay_id, starttime, sofa_24hours AS sofa_total
    FROM mimiciv_derived.sofa
  ")
  sofa_at_t0 <- sofa_hourly %>%
    inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
    filter(starttime <= t0 + hours(1)) %>%
    group_by(stay_id) %>%
    slice_max(starttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(stay_id, sofa_t0 = sofa_total)

  apsiii <- dbGetQuery(con, "SELECT stay_id, apsiii FROM mimiciv_derived.apsiii")
  creat_ts <- dbGetQuery(con, "
    SELECT hadm_id, charttime, creatinine
    FROM mimiciv_derived.chemistry
    WHERE creatinine IS NOT NULL
  ")
  creat_at_t0 <- creat_ts %>%
    inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
    filter(charttime <= t0 + hours(2)) %>%
    group_by(stay_id) %>%
    slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(stay_id, creatinine_t0 = creatinine)

  lactate_ts <- dbGetQuery(con, "
    SELECT hadm_id, charttime, lactate
    FROM mimiciv_derived.bg
    WHERE lactate IS NOT NULL
  ")
  lactate_at_t0 <- lactate_ts %>%
    inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
    filter(charttime <= t0 + hours(4), charttime >= t0 - hours(12)) %>%
    group_by(stay_id) %>%
    slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(stay_id, lactate_t0 = lactate)

  uo_ts <- dbGetQuery(con, "
    SELECT stay_id, charttime, uo_mlkghr_6hr AS urineoutput_rate
    FROM mimiciv_derived.urine_output_rate
    WHERE uo_mlkghr_6hr IS NOT NULL
  ")
  uo_at_t0 <- uo_ts %>%
    inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
    filter(charttime <= t0 + hours(2)) %>%
    group_by(stay_id) %>%
    slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(stay_id, uo_rate_t0 = urineoutput_rate)

  charlson <- dbGetQuery(con, "
    SELECT hadm_id, charlson_comorbidity_index AS cci
    FROM mimiciv_derived.charlson
  ")
  vitals_ts <- dbGetQuery(con, "
    SELECT stay_id, charttime, mbp AS map
    FROM mimiciv_derived.vitalsign
    WHERE mbp IS NOT NULL
  ")
  map_at_t0 <- vitals_ts %>%
    inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
    filter(charttime <= t0 + hours(1), charttime >= t0 - hours(2)) %>%
    group_by(stay_id) %>%
    slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(stay_id, map_t0 = map)

  baseline <- cohort %>%
    left_join(sofa_at_t0, by = "stay_id") %>%
    left_join(apsiii, by = "stay_id") %>%
    left_join(creat_at_t0, by = "stay_id") %>%
    left_join(lactate_at_t0, by = "stay_id") %>%
    left_join(uo_at_t0, by = "stay_id") %>%
    left_join(charlson, by = "hadm_id") %>%
    left_join(map_at_t0, by = "stay_id") %>%
    mutate(
      sofa_equipoise = sofa_t0 >= 8 & sofa_t0 <= 14,
      sofa_low = sofa_t0 <= 7,
      sofa_high = sofa_t0 >= 15,
      sofa_band_prespecified = case_when(
        sofa_t0 < 9 ~ "SOFA_lt9",
        sofa_t0 >= 9 & sofa_t0 <= 13 ~ "SOFA_9_13",
        sofa_t0 > 13 ~ "SOFA_gt13",
        TRUE ~ NA_character_
      ),
      creat_3x_base = !is.na(creatinine_t0) & !is.na(creatinine_baseline) &
        creatinine_t0 >= 3 * creatinine_baseline,
      creat_4mg_rise = !is.na(creatinine_t0) &
        creatinine_t0 >= 4.0 &
        (creatinine_t0 - coalesce(creatinine_baseline, 0)) >= 0.5,
      crrt_free_days = if_else(died_28d, 0, pmax(0, 28 - coalesce(crrt_hours_from_t0 / 24, 0)))
    )

  saveRDS(baseline, file.path(data_dir, "baseline.rds"))
  write_csv(baseline, file.path(data_dir, "baseline.csv"))
  invisible(baseline)
}

run_legacy_step_with_data_dir <- function(script_path, data_dir) {
  legacy_data_dir <- file.path(root_dir, "data")
  code <- readLines(script_path, warn = FALSE)
  code <- gsub(legacy_data_dir, data_dir, code, fixed = TRUE)
  # RUN-16 沙箱重定向：legacy 脚本内另有工作区外绝对路径（广西医科大学/analysis/data=软链，文件存在）
  # ⚠️ 2026-10-16 更正：该路径非失效；仅因文件沙箱禁止写工作区外 ⇒ 重定向到本 variant 的 data_dir
  code <- gsub(file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "data"), data_dir, code, fixed = TRUE)
  code <- gsub("slice_max(charttime, n=1)", "slice_max(charttime, n=1, with_ties=FALSE)", code, fixed = TRUE)
  code <- gsub("slice_max(starttime, n=1)", "slice_max(starttime, n=1, with_ties=FALSE)", code, fixed = TRUE)
  step_env <- new.env(parent = globalenv())
  eval(parse(text = code), envir = step_env)
}

run_cox <- function(ccw, label) {
  dat <- ccw %>%
    distinct(clone_id, .keep_all = TRUE) %>%
    mutate(trt = as.integer(strategy == "A_early")) %>%
    filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim))

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
    sofa8_14_n = n_distinct(dat$stay_id[dat$equipoise %in% TRUE]),
    sofa8_14_early_n = n_distinct(dat$stay_id[dat$equipoise %in% TRUE & dat$early_crrt %in% TRUE]),
    events = sum(dat$event, na.rm = TRUE),
    weighted_hr = unname(s$coefficients["trt", "exp(coef)"]),
    ci_lower = unname(s$conf.int["trt", "lower .95"]),
    ci_upper = unname(s$conf.int["trt", "upper .95"]),
    p_value = unname(s$coefficients["trt", "Pr(>|z|)"])
  )
}

make_variant <- function(name, cohort) {
  variant_dir <- file.path(strict_dir, paste0("rescue_sensitivity_", name))
  data_dir <- file.path(variant_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(cohort, file.path(data_dir, "cohort_base.rds"))
  write_csv(cohort, file.path(data_dir, "cohort_base.csv"))
  build_baseline(cohort, data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03a_timevarying_covariates.R"), data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03_ccw_ipcw.R"), data_dir)

  ccw <- readRDS(file.path(data_dir, "ccw_final.rds"))
  bind_rows(
    run_cox(ccw, paste0(name, ": full cohort")),
    run_cox(filter(ccw, equipoise %in% TRUE), paste0(name, ": SOFA 8-14"))
  )
}

variant_cohorts <- list(
  no_rescue_exclusion = pre_rescue_ext,
  primary_R1_R2_R2b_R3 = primary_cohort,
  sensitivity_R1_R2_R2b_R3_R4 = pre_rescue_ext %>% filter(!rescue_r1_r4_sensitivity),
  sensitivity_R1_confirm_R2_R2b_R3 = pre_rescue_ext %>%
    mutate(r1_confirm_trigger = r1_trigger & !r1_single_unrechecked_k_high) %>%
    filter(!(r1_confirm_trigger | r2_trigger | r2b_trigger | r3_trigger)),
  sensitivity_liberal_thresholds_R1_R2_R2b_R3 = pre_rescue_ext %>%
    filter(!rescue_threshold_liberal),
  sensitivity_conservative_thresholds_R1_R2_R2b_R3 = pre_rescue_ext %>%
    filter(!rescue_threshold_conservative)
)

primary_old <- readRDS(file.path(results_dir, "strictK3_primary_v2_mortality_results.rds"))
primary_ch01 <- read_csv(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ch01_primary_dod.csv"), show_col_types = FALSE)
primary_results <- primary_old %>%
  filter(analysis %in% c("Full strictK3_v2 primary cohort", "Primary SOFA 8-14 equipoise")) %>%
  mutate(pop = recode(
    analysis,
    "Full strictK3_v2 primary cohort" = "Full cohort",
    "Primary SOFA 8-14 equipoise" = "Equipoise SOFA 8-14"
  )) %>%
  left_join(
    primary_ch01 %>% select(pop = population, events_ch01 = events, hr_ch01 = hr, lo_ch01 = ci_lower, hi_ch01 = ci_upper, p_ch01 = p_value),
    by = "pop"
  ) %>%
  transmute(
    analysis = recode(
      pop,
      "Full cohort" = "primary_R1_R2_R2b_R3: full cohort",
      "Equipoise SOFA 8-14" = "primary_R1_R2_R2b_R3: SOFA 8-14"
    ),
    n_patients, clone_rows, early_n,
    sofa8_14_n = NA_integer_,
    sofa8_14_early_n = NA_integer_,
    events = events_ch01,
    weighted_hr = hr_ch01, ci_lower = lo_ch01, ci_upper = hi_ch01, p_value = p_ch01
  )

results <- list()
for (nm in names(variant_cohorts)) {
  if (nm == "primary_R1_R2_R2b_R3") {
    results[[nm]] <- primary_results
  } else {
    message("\n=== Running rescue sensitivity variant: ", nm, " ===")
    results[[nm]] <- make_variant(nm, variant_cohorts[[nm]])
  }
}

out <- bind_rows(results) %>%
  mutate(
    weighted_hr = round(weighted_hr, 3),
    ci_lower = round(ci_lower, 3),
    ci_upper = round(ci_upper, 3),
    p_value = round(p_value, 4)
  )

write_csv(out, file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.csv"))
saveRDS(out, file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.rds"))

report <- c(
  "# Strict clean first-KDIGO3 rescue sensitivity models",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Each non-reference rescue variant was re-run through baseline extraction, time-varying covariates, CCW/IPCW, and weighted Cox in a separate data directory.",
  "",
  paste(capture.output(print(out, n = nrow(out))), collapse = "\n"),
  "",
  "Outputs:",
  "",
  paste0("- `", file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.csv"), "`"),
  paste0("- Variant data directories under `", strict_dir, "/rescue_sensitivity_*`.")
)
writeLines(report, file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_report.md"))

cat("\n=== Rescue sensitivity model results ===\n")
print(out, n = nrow(out))
cat("\nSaved report:", file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_report.md"), "\n")