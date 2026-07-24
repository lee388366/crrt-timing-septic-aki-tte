# =============================================================================
# Strict clean first-KDIGO3 project: anchor sensitivity models
#
# Anchors:
#   - KDIGO3_minus6h: legacy sensitivity, t0 = max(ICU intime, first KDIGO3 - 6h)
#   - KDIGO2: exploratory secondary, t0 = first KDIGO2; simultaneous stage3 excluded
#
# These anchor comparisons are intentionally no-rescue sensitivity analyses:
# rescue_dictionary_v1.0_frozen.md is defined for strict clean first-KDIGO3 t0 and is
# not silently imposed on KDIGO2 or the legacy -6h anchor.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(survival)
})

root_dir <- "analysis"
strict_dir <- file.path(root_dir, "strictK3_primary_v2")
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

run <- function(sql) invisible(dbExecute(con, sql))

extract_anchor_cohorts <- function() {
  run("DROP TABLE IF EXISTS anchor_sens_a")
  run("DROP TABLE IF EXISTS anchor_sens_source")
  run("DROP TABLE IF EXISTS anchor_sens_eligible")

  run("
    CREATE TEMP TABLE anchor_sens_a AS
    WITH icu AS (
      SELECT stay_id, intime
      FROM mimiciv_icu.icustays
    ),
    k3 AS (
      SELECT DISTINCT ON (stay_id)
        stay_id,
        charttime AS kdigo3_time,
        aki_stage_creat,
        aki_stage_uo,
        aki_stage_crrt,
        uo_rt_12hr,
        uo_rt_24hr
      FROM mimiciv_derived.kdigo_stages
      WHERE aki_stage_creat = 3 OR aki_stage_uo = 3
      ORDER BY stay_id, charttime
    ),
    k2 AS (
      SELECT DISTINCT ON (stay_id)
        stay_id,
        charttime AS kdigo2_time,
        aki_stage_creat,
        aki_stage_uo,
        aki_stage_crrt,
        uo_rt_12hr,
        uo_rt_24hr
      FROM mimiciv_derived.kdigo_stages
      WHERE aki_stage_creat = 2 OR aki_stage_uo = 2
      ORDER BY stay_id, charttime
    )
    SELECT
      'KDIGO3_minus6h'::text AS anchor,
      k3.stay_id,
      k3.kdigo3_time,
      NULL::timestamp AS kdigo2_time,
      greatest(icu.intime, k3.kdigo3_time - interval '6 hours') AS t0,
      k3.aki_stage_creat,
      k3.aki_stage_uo,
      k3.aki_stage_crrt,
      k3.uo_rt_12hr,
      k3.uo_rt_24hr
    FROM k3
    JOIN icu USING (stay_id)
    UNION ALL
    SELECT
      'KDIGO2'::text AS anchor,
      k2.stay_id,
      NULL::timestamp AS kdigo3_time,
      k2.kdigo2_time,
      k2.kdigo2_time AS t0,
      k2.aki_stage_creat,
      k2.aki_stage_uo,
      k2.aki_stage_crrt,
      k2.uo_rt_12hr,
      k2.uo_rt_24hr
    FROM k2
  ")
  run("CREATE INDEX ON anchor_sens_a(anchor, stay_id)")
  run("CREATE INDEX ON anchor_sens_a(stay_id)")

  run("
    CREATE TEMP TABLE anchor_sens_source AS
    WITH sepsis_at_t0 AS (
      SELECT
        a.anchor,
        a.stay_id,
        min(s.sofa_time) AS sofa_time,
        (array_agg(s.sofa_score ORDER BY s.sofa_time))[1] AS sofa_score
      FROM anchor_sens_a a
      JOIN mimiciv_derived.sepsis3 s
        ON s.stay_id = a.stay_id
       AND s.sepsis3 = true
       AND s.sofa_time <= a.t0
      GROUP BY a.anchor, a.stay_id
    ),
    vaso_window AS (
      SELECT
        a.anchor,
        a.stay_id,
        true AS vaso_in_window,
        avg(
          coalesce(v.norepinephrine, 0)
          + coalesce(v.epinephrine, 0)
          + coalesce(v.dopamine, 0) / 100.0
          + coalesce(v.phenylephrine, 0) * 0.1
          + (coalesce(v.vasopressin, 0) / 60.0) * 2.5
        ) AS nee_at_t0
      FROM anchor_sens_a a
      JOIN mimiciv_derived.vasoactive_agent v
        ON v.stay_id = a.stay_id
       AND v.starttime <= a.t0
       AND v.endtime >= a.t0 - interval '6 hours'
      GROUP BY a.anchor, a.stay_id
    ),
    crrt AS (
      SELECT stay_id, min(charttime) AS crrt_first_time
      FROM mimiciv_derived.rrt
      WHERE dialysis_active = 1
        AND dialysis_type = 'CRRT'
      GROUP BY stay_id
    )
    SELECT
      i.subject_id,
      i.hadm_id,
      i.stay_id,
      i.intime AS icu_intime,
      i.outtime AS icu_outtime,
      extract(epoch from (i.outtime - i.intime)) / 3600 AS icu_los_hours,
      p.anchor_age AS age,
      p.gender,
      p.anchor_year,
      p.anchor_year_group,
      ad.hospital_expire_flag,
      ad.admittime,
      ad.dischtime,
      ad.deathtime,
      a.anchor,
      a.kdigo3_time,
      a.kdigo2_time,
      a.t0,
      (a.aki_stage_creat = 3) AS stage3_by_creat,
      (a.aki_stage_uo = 3) AS stage3_by_uo,
      (a.aki_stage_crrt = 3) AS stage3_by_crrt,
      a.uo_rt_12hr,
      a.uo_rt_24hr,
      cb.scr_baseline AS creatinine_baseline,
      c.crrt_first_time,
      s.sofa_time,
      s.sofa_score,
      v.vaso_in_window,
      v.nee_at_t0
    FROM anchor_sens_a a
    JOIN mimiciv_icu.icustays i USING (stay_id)
    JOIN mimiciv_hosp.patients p USING (subject_id)
    JOIN mimiciv_hosp.admissions ad USING (hadm_id)
    JOIN sepsis_at_t0 s USING (anchor, stay_id)
    JOIN vaso_window v USING (anchor, stay_id)
    LEFT JOIN mimiciv_derived.creatinine_baseline cb USING (hadm_id)
    LEFT JOIN crrt c USING (stay_id)
    WHERE p.anchor_age >= 18
      AND NOT (
        a.anchor = 'KDIGO2'
        AND EXISTS (
          SELECT 1
          FROM mimiciv_derived.kdigo_stages k
          WHERE k.stay_id = a.stay_id
            AND k.charttime = a.t0
            AND (k.aki_stage_creat = 3 OR k.aki_stage_uo = 3)
        )
      )
  ")
  run("CREATE INDEX ON anchor_sens_source(anchor, stay_id)")

  run("
    CREATE TEMP TABLE anchor_sens_eligible AS
    SELECT src.*
    FROM anchor_sens_source src
    WHERE (crrt_first_time IS NULL OR crrt_first_time >= t0)
      AND NOT EXISTS (
        SELECT 1
        FROM mimiciv_hosp.diagnoses_icd d
        WHERE d.hadm_id = src.hadm_id
          AND (
            (d.icd_version = 10 AND d.icd_code IN ('N186', 'N185'))
            OR (d.icd_version = 9 AND d.icd_code IN ('5856', 'V4511', 'V4512'))
          )
      )
      AND icu_los_hours >= 24
      AND creatinine_baseline IS NOT NULL
      AND nee_at_t0 IS NOT NULL
      AND t0 IS NOT NULL
      AND greatest(
        0,
        least(
          28,
          coalesce(
            extract(epoch from (deathtime - t0)) / 86400,
            extract(epoch from (dischtime - t0)) / 86400,
            28
          )
        )
      ) > 0
  ")

  dbGetQuery(con, "
    SELECT DISTINCT ON (anchor, subject_id)
      *,
      extract(epoch from (crrt_first_time - t0)) / 3600 AS crrt_hours_from_t0,
      coalesce(crrt_first_time >= t0 AND crrt_first_time <= t0 + interval '24 hours', false) AS early_crrt,
      coalesce(crrt_first_time < t0, false) AS crrt_before_t0,
      extract(epoch from (deathtime - t0)) / 86400 AS days_to_death,
      coalesce(deathtime <= t0 + interval '28 days', false) AS died_28d,
      greatest(
        0,
        least(
          28,
          coalesce(
            extract(epoch from (deathtime - t0)) / 86400,
            extract(epoch from (dischtime - t0)) / 86400,
            28
          )
        )
      ) AS follow_time_28d,
      extract(epoch from (icu_outtime - t0)) / 86400 AS icu_los_from_t0,
      anchor AS t0_anchor,
      'no_rescue_anchor_sensitivity'::text AS rescue_definition
    FROM anchor_sens_eligible
    ORDER BY anchor, subject_id, icu_intime
  ")
}

build_baseline <- function(cohort_raw, data_dir) {
  cohort <- cohort_raw %>%
    group_by(subject_id) %>%
    slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
    ungroup()

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
  code <- gsub("slice_min(icu_intime, n=1)", "slice_min(icu_intime, n=1, with_ties=FALSE)", code, fixed = TRUE)
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
    events = sum(dat$event, na.rm = TRUE),
    weighted_hr = unname(s$coefficients["trt", "exp(coef)"]),
    ci_lower = unname(s$conf.int["trt", "lower .95"]),
    ci_upper = unname(s$conf.int["trt", "upper .95"]),
    p_value = unname(s$coefficients["trt", "Pr(>|z|)"])
  )
}

weight_diag <- function(ccw, label) {
  ccw %>%
    distinct(clone_id, .keep_all = TRUE) %>%
    group_by(strategy, equipoise) %>%
    summarise(
      anchor = label,
      n_patients = n_distinct(stay_id),
      mean_w = mean(ipcw_w_trim, na.rm = TRUE),
      p99_w = unname(quantile(ipcw_w_trim, 0.99, na.rm = TRUE)),
      max_w = max(ipcw_w_trim, na.rm = TRUE),
      ess = sum(ipcw_w_trim, na.rm = TRUE)^2 / sum(ipcw_w_trim^2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(anchor, everything())
}

cohorts <- extract_anchor_cohorts()
write_csv(cohorts, file.path(results_dir, "strictK3_v2_anchor_sensitivity_cohort_rows.csv"))

anchor_names <- c("KDIGO3_minus6h", "KDIGO2")
results <- list()
weights <- list()
summaries <- list()

for (anchor_name in anchor_names) {
  message("\n=== Running anchor sensitivity: ", anchor_name, " ===")
  cohort <- cohorts %>% filter(anchor == anchor_name)
  variant_dir <- file.path(strict_dir, paste0("anchor_sensitivity_", anchor_name))
  data_dir <- file.path(variant_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(cohort, file.path(data_dir, "cohort_base.rds"))
  write_csv(cohort, file.path(data_dir, "cohort_base.csv"))
  build_baseline(cohort, data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03a_timevarying_covariates.R"), data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03_ccw_ipcw.R"), data_dir)

  ccw <- readRDS(file.path(data_dir, "ccw_final.rds"))
  results[[anchor_name]] <- bind_rows(
    run_cox(ccw, paste0(anchor_name, ": full cohort")),
    run_cox(filter(ccw, equipoise %in% TRUE), paste0(anchor_name, ": SOFA 8-14"))
  )
  weights[[anchor_name]] <- weight_diag(ccw, anchor_name)
  summaries[[anchor_name]] <- ccw %>%
    distinct(stay_id, .keep_all = TRUE) %>%
    summarise(
      anchor = anchor_name,
      n = n(),
      early_n = sum(early_crrt, na.rm = TRUE),
      early_pct = round(100 * mean(early_crrt, na.rm = TRUE), 1),
      sofa8_14_n = sum(equipoise, na.rm = TRUE),
      sofa8_14_early_n = sum(equipoise & early_crrt, na.rm = TRUE),
      sofa8_14_early_pct = round(100 * mean(early_crrt[equipoise], na.rm = TRUE), 1)
    )
}

strict_primary <- readRDS(file.path(results_dir, "strictK3_primary_v2_mortality_results.rds")) %>%
  filter(analysis %in% c("Full strictK3_v2 primary cohort", "Primary SOFA 8-14 equipoise")) %>%
  transmute(
    analysis = recode(
      analysis,
      "Full strictK3_v2 primary cohort" = "strict_KDIGO3_primary_R1_R2_R2b_R3: full cohort",
      "Primary SOFA 8-14 equipoise" = "strict_KDIGO3_primary_R1_R2_R2b_R3: SOFA 8-14"
    ),
    n_patients, clone_rows, early_n, events, weighted_hr, ci_lower, ci_upper, p_value
  )

anchor_results <- bind_rows(results) %>%
  mutate(
    weighted_hr = round(weighted_hr, 3),
    ci_lower = round(ci_lower, 3),
    ci_upper = round(ci_upper, 3),
    p_value = round(p_value, 4)
  )

comparison_results <- bind_rows(
  strict_primary %>%
    mutate(
      weighted_hr = round(weighted_hr, 3),
      ci_lower = round(ci_lower, 3),
      ci_upper = round(ci_upper, 3),
      p_value = round(p_value, 4)
    ),
  anchor_results
)

anchor_summary <- bind_rows(summaries)
weight_summary <- bind_rows(weights)

write_csv(anchor_results, file.path(results_dir, "strictK3_v2_anchor_sensitivity_model_results.csv"))
write_csv(comparison_results, file.path(results_dir, "strictK3_v2_anchor_sensitivity_comparison_results.csv"))
write_csv(anchor_summary, file.path(results_dir, "strictK3_v2_anchor_sensitivity_cohort_summary.csv"))
write_csv(weight_summary, file.path(results_dir, "strictK3_v2_anchor_sensitivity_weight_summary.csv"))
saveRDS(
  list(results = anchor_results, comparison = comparison_results, cohort_summary = anchor_summary, weights = weight_summary),
  file.path(results_dir, "strictK3_v2_anchor_sensitivity_results.rds")
)

report <- c(
  "# Strict clean first-KDIGO3 anchor sensitivity models",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Interpretation guardrail",
  "",
  "- `KDIGO3_minus6h` is a legacy sensitivity anchor.",
  "- `KDIGO2` is exploratory secondary.",
  "- Both anchor sensitivity models are no-rescue comparisons because `rescue_dictionary_v1.0_frozen.md` is defined for strict clean first-KDIGO3 t0 and was not silently imposed on other anchors.",
  "- Baseline/time-varying covariates, CCW/IPCW, ridge models, stabilization, and trimming reuse the existing 03a/03 pipeline via path redirection.",
  "",
  "## Cohort summary",
  "",
  paste(capture.output(print(anchor_summary, n = nrow(anchor_summary))), collapse = "\n"),
  "",
  "## Cox comparison",
  "",
  paste(capture.output(print(comparison_results, n = nrow(comparison_results))), collapse = "\n"),
  "",
  "## Weight summary",
  "",
  paste(capture.output(print(weight_summary, n = nrow(weight_summary))), collapse = "\n"),
  "",
  "Outputs:",
  "",
  paste0("- `", file.path(results_dir, "strictK3_v2_anchor_sensitivity_model_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_anchor_sensitivity_comparison_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_anchor_sensitivity_cohort_summary.csv"), "`"),
  paste0("- Anchor data directories under `", strict_dir, "/anchor_sensitivity_*`.")
)
writeLines(report, file.path(results_dir, "strictK3_v2_anchor_sensitivity_report.md"))

cat("\n=== Anchor sensitivity comparison ===\n")
print(comparison_results, n = nrow(comparison_results))
cat("\nSaved report:", file.path(results_dir, "strictK3_v2_anchor_sensitivity_report.md"), "\n")
