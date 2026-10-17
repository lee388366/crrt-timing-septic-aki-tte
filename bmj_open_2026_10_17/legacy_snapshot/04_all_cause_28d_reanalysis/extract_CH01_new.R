# Strict clean first-KDIGO3 cohort extraction for the primary TTE analysis.
# Scope: cohort, rescue flags, and funnel counts only. Existing analysis scripts
# are not modified and legacy analysis/data/*.rds files are not overwritten.
suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(tidyr)
  library(readr)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
out_dir <- file.path(root_dir, "strictK3_primary_v2")
data_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/extract_ch01")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

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

run("
CREATE TEMP TABLE k3 AS
SELECT DISTINCT ON (stay_id)
  stay_id,
  charttime AS kdigo3_time,
  charttime AS t0,
  aki_stage_creat,
  aki_stage_uo,
  aki_stage_crrt,
  uo_rt_12hr,
  uo_rt_24hr
FROM mimiciv_derived.kdigo_stages
WHERE aki_stage_creat = 3 OR aki_stage_uo = 3
ORDER BY stay_id, charttime
")
run("CREATE INDEX ON k3(stay_id)")

run("
CREATE TEMP TABLE adult_k3 AS
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
  k3.kdigo3_time,
  k3.t0,
  bool_or(k3.aki_stage_creat = 3) AS stage3_by_creat,
  bool_or(k3.aki_stage_uo = 3) AS stage3_by_uo,
  bool_or(k3.aki_stage_crrt = 3) AS stage3_by_crrt,
  max(k3.uo_rt_12hr) AS uo_rt_12hr,
  max(k3.uo_rt_24hr) AS uo_rt_24hr
FROM k3
JOIN mimiciv_icu.icustays i USING (stay_id)
JOIN mimiciv_hosp.patients p USING (subject_id)
JOIN mimiciv_hosp.admissions ad USING (hadm_id)
WHERE p.anchor_age >= 18
GROUP BY
  i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime,
  p.anchor_age, p.gender, p.anchor_year, p.anchor_year_group,
  ad.hospital_expire_flag, ad.admittime, ad.dischtime, ad.deathtime,
  k3.kdigo3_time, k3.t0
")
run("CREATE INDEX ON adult_k3(stay_id)")
run("CREATE INDEX ON adult_k3(hadm_id, t0)")

run("
CREATE TEMP TABLE sepsis_at_t0 AS
SELECT
  a.stay_id,
  min(s.sofa_time) AS sofa_time,
  (array_agg(s.sofa_score ORDER BY s.sofa_time))[1] AS sofa_score
FROM adult_k3 a
JOIN mimiciv_derived.sepsis3 s
  ON s.stay_id = a.stay_id
 AND s.sepsis3 = true
 AND s.sofa_time <= a.t0
GROUP BY a.stay_id
")
run("CREATE INDEX ON sepsis_at_t0(stay_id)")

run("
CREATE TEMP TABLE vaso_window AS
SELECT
  a.stay_id,
  true AS vaso_in_window,
  avg(
    coalesce(v.norepinephrine, 0)
    + coalesce(v.epinephrine, 0)
    + coalesce(v.dopamine, 0) / 100.0
    + coalesce(v.phenylephrine, 0) * 0.1
    + (coalesce(v.vasopressin, 0) / 60.0) * 2.5
  ) AS nee_at_t0
FROM adult_k3 a
JOIN mimiciv_derived.vasoactive_agent v
  ON v.stay_id = a.stay_id
 AND v.starttime <= a.t0
 AND v.endtime >= a.t0 - interval '6 hours'
GROUP BY a.stay_id
")
run("CREATE INDEX ON vaso_window(stay_id)")

run("
CREATE TEMP TABLE crrt AS
SELECT stay_id, min(charttime) AS crrt_first_time
FROM mimiciv_derived.rrt
WHERE dialysis_active = 1
  AND dialysis_type = 'CRRT'
GROUP BY stay_id
")
run("CREATE INDEX ON crrt(stay_id)")

run("
CREATE TEMP TABLE source AS
SELECT
  a.*,
  cb.scr_baseline AS creatinine_baseline,
  c.crrt_first_time,
  s.sofa_time,
  s.sofa_score,
  v.vaso_in_window,
  v.nee_at_t0
FROM adult_k3 a
JOIN sepsis_at_t0 s USING (stay_id)
JOIN vaso_window v USING (stay_id)
LEFT JOIN mimiciv_derived.creatinine_baseline cb USING (hadm_id)
LEFT JOIN crrt c USING (stay_id)
")
run("CREATE INDEX ON source(stay_id)")
run("CREATE INDEX ON source(hadm_id, t0)")

source_counts <- dbGetQuery(con, "
WITH
  n_k3 AS (SELECT count(*)::integer AS n FROM k3),
  n_adult AS (SELECT count(*)::integer AS n FROM adult_k3),
  n_sepsis AS (
    SELECT count(*)::integer AS n
    FROM adult_k3 a
    JOIN sepsis_at_t0 s USING (stay_id)
  ),
  n_vaso AS (SELECT count(*)::integer AS n FROM source),
  n_no_prior_rrt AS (
    SELECT count(*)::integer AS n
    FROM source
    WHERE crrt_first_time IS NULL OR crrt_first_time >= t0
  ),
  n_no_esrd AS (
    SELECT count(*)::integer AS n
    FROM source src
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
  ),
  n_los AS (
    SELECT count(*)::integer AS n
    FROM source src
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
  ),
  n_baseline AS (
    SELECT count(*)::integer AS n
    FROM source src
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
  ),
  n_followup AS (
    SELECT count(*)::integer AS n
    FROM source src
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
          coalesce(extract(epoch from (deathtime - t0)) / 86400,
                   extract(epoch from (dischtime - t0)) / 86400,
                   28)
        )
      ) > 0
  )
SELECT 'first_KDIGO3_stays' AS step, n FROM n_k3
UNION ALL SELECT 'adult_KDIGO3_stays', n FROM n_adult
UNION ALL SELECT 'sepsis3_at_or_before_t0', n FROM n_sepsis
UNION ALL SELECT 'vasopressor_window_existing_definition', n FROM n_vaso
UNION ALL SELECT 'exclude_RRT_before_t0', n FROM n_no_prior_rrt
UNION ALL SELECT 'exclude_ESRD_chronic_dialysis', n FROM n_no_esrd
UNION ALL SELECT 'exclude_ICU_LOS_lt_24h', n FROM n_los
UNION ALL SELECT 'exclude_missing_baseline_essentials', n FROM n_baseline
UNION ALL SELECT 'exclude_zero_28d_followup', n FROM n_followup
")

run("
CREATE TEMP TABLE pre_rescue_all_stays AS
SELECT src.*
FROM source src
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
      coalesce(extract(epoch from (deathtime - t0)) / 86400,
               extract(epoch from (dischtime - t0)) / 86400,
               28)
    )
  ) > 0
")
run("CREATE INDEX ON pre_rescue_all_stays(subject_id, icu_intime)")
run("CREATE INDEX ON pre_rescue_all_stays(hadm_id, t0)")

run("
CREATE TEMP TABLE pre_rescue AS
SELECT DISTINCT ON (subject_id) *
FROM pre_rescue_all_stays
ORDER BY subject_id, icu_intime
")
run("CREATE INDEX ON pre_rescue(stay_id)")
run("CREATE INDEX ON pre_rescue(hadm_id, t0)")

run("
CREATE TEMP TABLE rescue_flags AS
SELECT
  pr.*,
  coalesce(kstats.k_n, 0) AS r1_k_n,
  coalesce(kstats.k_high_n, 0) AS r1_k_high_n,
  coalesce(kstats.r1_trigger, false) AS r1_trigger,
  coalesce(abg.abg_ph_n, 0) AS r2_abg_ph_n,
  coalesce(abg.r2_trigger, false) AS r2_trigger,
  coalesce(nonart.nonarterial_ph_n, 0) AS r2_nonarterial_ph_n,
  coalesce(hco3.hco3_n, 0) AS r2b_hco3_n,
  coalesce(hco3.r2b_trigger, false) AS r2b_trigger,
  coalesce(bun.bun_n, 0) AS r3_bun_n,
  coalesce(bun.r3_trigger, false) AS r3_trigger,
  (pr.uo_rt_12hr IS NOT NULL OR pr.uo_rt_24hr IS NOT NULL) AS r4_uo_measured,
  (coalesce(pr.uo_rt_24hr < 0.3, false) OR coalesce(pr.uo_rt_12hr = 0, false)) AS r4_trigger,
  coalesce(pr.uo_rt_12hr = 0, false) AS r4_anuria12_proxy,
  coalesce(pr.uo_rt_24hr < 0.3, false) AS r4_severe_oliguria24
FROM pre_rescue pr
LEFT JOIN LATERAL (
  SELECT
    count(*) FILTER (WHERE c.potassium IS NOT NULL) AS k_n,
    count(*) FILTER (WHERE c.potassium >= 6.5) AS k_high_n,
    bool_or(c.potassium >= 6.5) AS r1_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '6 hours'
    AND c.charttime <= pr.t0
    AND c.potassium IS NOT NULL
) kstats ON true
LEFT JOIN LATERAL (
  SELECT
    count(*) FILTER (WHERE bg.ph IS NOT NULL) AS abg_ph_n,
    bool_or(bg.ph <= 7.15) AS r2_trigger
  FROM mimiciv_derived.bg bg
  WHERE bg.hadm_id = pr.hadm_id
    AND bg.charttime > pr.t0 - interval '6 hours'
    AND bg.charttime <= pr.t0
    AND bg.specimen = 'ART.'
    AND bg.ph IS NOT NULL
) abg ON true
LEFT JOIN LATERAL (
  SELECT count(*) FILTER (WHERE bg.ph IS NOT NULL) AS nonarterial_ph_n
  FROM mimiciv_derived.bg bg
  WHERE bg.hadm_id = pr.hadm_id
    AND bg.charttime > pr.t0 - interval '6 hours'
    AND bg.charttime <= pr.t0
    AND bg.specimen <> 'ART.'
    AND bg.ph IS NOT NULL
) nonart ON true
LEFT JOIN LATERAL (
  SELECT
    count(*) FILTER (WHERE c.bicarbonate IS NOT NULL) AS hco3_n,
    bool_or(c.bicarbonate <= 12) AS r2b_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '6 hours'
    AND c.charttime <= pr.t0
    AND c.bicarbonate IS NOT NULL
) hco3 ON true
LEFT JOIN LATERAL (
  SELECT
    count(*) FILTER (WHERE c.bun IS NOT NULL) AS bun_n,
    bool_or(c.bun >= 100) AS r3_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '24 hours'
    AND c.charttime <= pr.t0
    AND c.bun IS NOT NULL
) bun ON true
")

run("
CREATE TEMP TABLE cohort_flags AS
SELECT
  rf.*,
  (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger) AS rescue_primary_r1_r2_r2b_r3,
  (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger OR r4_trigger) AS rescue_r1_r4_sensitivity,
  (r1_trigger AND r1_k_n = 1) AS r1_single_unrechecked_k_high,
  coalesce(crrt_first_time BETWEEN t0 AND t0 + interval '24 hours', false) AS early_crrt,
  (crrt_first_time IS NOT NULL AND crrt_first_time < t0) AS crrt_before_t0,
  extract(epoch from (crrt_first_time - t0)) / 3600 AS crrt_hours_from_t0,
  extract(epoch from (deathtime - t0)) / 86400 AS days_to_death,
  (deathtime IS NOT NULL AND extract(epoch from (deathtime - t0)) / 86400 <= 28) AS died_28d,
  greatest(
    0,
    least(
      28,
      coalesce(extract(epoch from (deathtime - t0)) / 86400,
               extract(epoch from (dischtime - t0)) / 86400,
               28)
    )
  ) AS follow_time_28d,
  extract(epoch from (icu_outtime - t0)) / 86400 AS icu_los_from_t0
FROM rescue_flags rf
")

cohort_all <- dbGetQuery(con, "SELECT * FROM cohort_flags ORDER BY subject_id")

# ===== CH-01 最小改动：仅重定义 died_28d / follow_time_28d 两列 =====
.dod <- dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients")
cohort_all <- cohort_all %>%
  left_join(.dod, by = "subject_id") %>%
  mutate(
    .dc = coalesce(days_to_death, as.numeric(dod - as.Date(t0))),
    died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
    follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
  ) %>%
  select(-.dc, -dod)
cat(sprintf("[CH-01 patch] died_28d = %d / %d\n", sum(cohort_all$died_28d, na.rm=TRUE), nrow(cohort_all)))
# ===== 补丁结束 =====
cohort_primary <- cohort_all %>%
  filter(!rescue_primary_r1_r2_r2b_r3) %>%
  mutate(
    stage3_path = case_when(
      stage3_by_creat %in% TRUE & stage3_by_uo %in% TRUE ~ "both_creatinine_and_uo",
      stage3_by_creat %in% TRUE ~ "creatinine_only",
      stage3_by_uo %in% TRUE ~ "urine_output_only",
      TRUE ~ "no_creat_or_uo_stage3_ERROR"
    ),
    t0_anchor = "strict_first_KDIGO3_clean_CrUO",
    rescue_definition = "primary_R1_R2_R2b_R3"
  ) %>%
  select(
    subject_id, hadm_id, stay_id, icu_intime, icu_outtime, icu_los_hours,
    age, gender, anchor_year, anchor_year_group, hospital_expire_flag,
    admittime, dischtime, deathtime, kdigo3_time, stage3_by_creat,
    stage3_by_uo, stage3_by_crrt, stage3_path, uo_rt_12hr, uo_rt_24hr,
    creatinine_baseline, crrt_first_time, t0, sofa_time, sofa_score,
    vaso_in_window, nee_at_t0, crrt_hours_from_t0, early_crrt,
    crrt_before_t0, days_to_death, died_28d, follow_time_28d,
    icu_los_from_t0, starts_with("r1_"), starts_with("r2_"),
    starts_with("r2b_"), starts_with("r3_"), starts_with("r4_"),
    rescue_primary_r1_r2_r2b_r3, rescue_r1_r4_sensitivity,
    t0_anchor, rescue_definition
  )

saveRDS(cohort_all, file.path(data_dir, "cohort_pre_rescue_all.rds"))
saveRDS(cohort_primary, file.path(data_dir, "cohort_base.rds"))
write_csv(cohort_all, file.path(data_dir, "cohort_pre_rescue_all.csv"))
write_csv(cohort_primary, file.path(data_dir, "cohort_base.csv"))

run("CREATE SCHEMA IF NOT EXISTS data_extract_crrt_strictk3_v2")
dbWriteTable(con, Id(schema = "data_extract_crrt_strictk3_v2", table = "cohort_pre_rescue_all_ch01"),
             cohort_all, overwrite = TRUE)
dbWriteTable(con, Id(schema = "data_extract_crrt_strictk3_v2", table = "cohort_base_ch01"),
             cohort_primary, overwrite = TRUE)

source_counts <- source_counts %>%
  mutate(
    previous_n = lag(n),
    excluded_from_previous = if_else(is.na(previous_n), NA_integer_, previous_n - n)
  )

rescue_funnel <- tibble(
  step = c("pre_rescue_first_eligible_ICU", "exclude_primary_rescue_R1_R2_R2b_R3"),
  n = c(nrow(cohort_all), nrow(cohort_primary))
) %>%
  mutate(
    previous_n = c(tail(source_counts$n, 1), lag(n)[-1]),
    excluded_from_previous = if_else(is.na(previous_n), NA_integer_, previous_n - n)
  )
funnel <- bind_rows(source_counts, rescue_funnel)
write_csv(funnel, file.path(data_dir, "strictK3_primary_v2_funnel.csv"))

indicator_summary <- tibble(
  indicator = c("R1_hyperkalemia", "R2_arterial_acidemia", "R2b_bicarbonate", "R3_azotemia", "R4_UO_sensitivity"),
  main_or_sensitivity = c("main", "main", "main", "main", "sensitivity_only"),
  lookback_window = c("6h", "6h", "6h", "24h", "12h/24h derived UO rates at t0"),
  denominator_n = nrow(cohort_all),
  measured_n = c(
    sum(cohort_all$r1_k_n > 0, na.rm = TRUE),
    sum(cohort_all$r2_abg_ph_n > 0, na.rm = TRUE),
    sum(cohort_all$r2b_hco3_n > 0, na.rm = TRUE),
    sum(cohort_all$r3_bun_n > 0, na.rm = TRUE),
    sum(cohort_all$r4_uo_measured, na.rm = TRUE)
  ),
  trigger_n = c(
    sum(cohort_all$r1_trigger, na.rm = TRUE),
    sum(cohort_all$r2_trigger, na.rm = TRUE),
    sum(cohort_all$r2b_trigger, na.rm = TRUE),
    sum(cohort_all$r3_trigger, na.rm = TRUE),
    sum(cohort_all$r4_trigger, na.rm = TRUE)
  )
) %>%
  mutate(
    missing_n = denominator_n - measured_n,
    missing_pct = round(100 * missing_n / denominator_n, 1),
    trigger_pct_all = round(100 * trigger_n / denominator_n, 1),
    trigger_pct_measured = round(100 * trigger_n / measured_n, 1)
  )
write_csv(indicator_summary, file.path(data_dir, "strictK3_v2_rescue_indicator_summary.csv"))

version_summary <- bind_rows(
  cohort_all %>%
    summarise(version = "no_rescue_exclusion", n = n(), early_n = sum(early_crrt), early_pct = round(100 * mean(early_crrt), 1)),
  cohort_all %>%
    filter(!rescue_primary_r1_r2_r2b_r3) %>%
    summarise(version = "primary_exclude_R1_R2_R2b_R3", n = n(), early_n = sum(early_crrt), early_pct = round(100 * mean(early_crrt), 1)),
  cohort_all %>%
    filter(!rescue_r1_r4_sensitivity) %>%
    summarise(version = "sensitivity_exclude_R1_R4", n = n(), early_n = sum(early_crrt), early_pct = round(100 * mean(early_crrt), 1))
)
write_csv(version_summary, file.path(data_dir, "strictK3_v2_rescue_version_summary.csv"))

report <- c(
  "# Strict First-KDIGO3 Primary Cohort Extraction",
  "",
  "Scope: cohort extraction and rescue diagnostics only. CCW/IPCW, Cox, RMST, and competing-risk analyses are not run by this script.",
  "",
  "## Assumptions / Ambiguities",
  "",
  "- The existing pipeline's clinical definitions are reused: Sepsis-3 at or before t0 from `mimiciv_derived.sepsis3`, CRRT from `mimiciv_derived.rrt`, NEE from `mimiciv_derived.vasoactive_agent`. Vasopressin in `mimiciv_derived.vasoactive_agent` is stored as U/hour and is converted to U/min before applying the 2.5 norepinephrine-equivalent multiplier.",
  "- Clean t0 is changed to the first physiologic KDIGO3 charttime defined by `aki_stage_creat = 3 OR aki_stage_uo = 3`; `aki_stage_crrt` is explicitly not allowed to define t0.",
  "- This clean t0 replaces the legacy `max(icu_intime, KDIGO3_time - 6h)` anchor and removes CRRT-defined t0 circularity.",
  "- First qualifying ICU stay per patient is fixed after source eligibility/non-rescue exclusions and before rescue exclusion; rescue-positive first stays are not replaced by later stays.",
  "- Primary rescue exclusion follows `rescue_dictionary_v1.0_frozen.md`: R1 OR (R2 OR R2b) OR R3.",
  "- R4/UO, R5, and R6 are not used in the primary exclusion. R4 flags are retained only for sensitivity diagnostics.",
  "- R2 uses arterial pH only (`specimen = 'ART.'`). Nonarterial pH is counted but cannot trigger R2.",
  "- Lab triggers use any qualifying measurement in the lookback window, matching the frozen dictionary's '达标即触发' rule; missing measurements are treated as not triggered and reported.",
  "",
  "## Funnel",
  "",
  "| step | n | excluded from previous |",
  "|---|---:|---:|",
  apply(funnel, 1, function(x) sprintf("| %s | %s | %s |", x[["step"]], x[["n"]], ifelse(is.na(x[["excluded_from_previous"]]), "", x[["excluded_from_previous"]]))),
  "",
  "## Rescue Versions",
  "",
  "| version | n | early n | early % |",
  "|---|---:|---:|---:|",
  apply(version_summary, 1, function(x) sprintf("| %s | %s | %s | %s |", x[["version"]], x[["n"]], x[["early_n"]], x[["early_pct"]])),
  "",
  "## Per-Indicator Rescue Diagnostics",
  "",
  "| indicator | role | window | measured n | missing n | missing % | trigger n | trigger % all | trigger % measured |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|",
  apply(indicator_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    x[["indicator"]], x[["main_or_sensitivity"]], x[["lookback_window"]],
    x[["measured_n"]], x[["missing_n"]], x[["missing_pct"]],
    x[["trigger_n"]], x[["trigger_pct_all"]], x[["trigger_pct_measured"]]
  )),
  "",
  "## Outputs",
  "",
  "- `data/cohort_base.rds`: strict-K3-v2 primary cohort after R1/R2/R2b/R3 rescue exclusion.",
  "- `data/cohort_pre_rescue_all.rds`: first eligible strict-K3-v2 cohort before rescue exclusion.",
  "- `data/strictK3_primary_v2_funnel.csv`",
  "- `data/strictK3_v2_rescue_indicator_summary.csv`",
  "- `data/strictK3_v2_rescue_version_summary.csv`",
  "- PostgreSQL schema: `data_extract_crrt_strictk3_v2`."
)
writeLines(report, file.path(data_dir, "strictK3_primary_v2_extraction_report.md"))

print(funnel)
print(version_summary)
print(indicator_summary)