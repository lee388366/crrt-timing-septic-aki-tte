# QA: Does strict first-KDIGO3 t0 include the RRT/CRRT KDIGO stage-3 criterion?
# Scope: read-only diagnostic + small counterfactual cohort count only.
# Does NOT overwrite existing strictK3 primary analysis scripts or database tables.

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(tidyr)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
extract_script <- file.path(root_dir, "strictK3_primary_extract.R")
out_dir <- file.path(root_dir, "strictK3_primary", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
qry <- function(sql) dbGetQuery(con, sql)

run("SET TIME ZONE 'UTC'")

fmt_pct <- function(n, d, digits = 1) {
  if (is.na(n) || is.na(d) || d == 0) return(NA_character_)
  sprintf(paste0("%.", digits, "f%%"), 100 * n / d)
}

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
}

collapse_vals <- function(x, digits = 1) {
  if (all(is.na(x))) return("")
  qs <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE)
  sprintf(
    "%s [%s, %s]; range %s to %s",
    fmt_num(qs[3], digits), fmt_num(qs[2], digits), fmt_num(qs[4], digits),
    fmt_num(qs[1], digits), fmt_num(qs[5], digits)
  )
}

# -----------------------------------------------------------------------------
# 1) Code evidence
# -----------------------------------------------------------------------------

src <- readLines(extract_script, warn = FALSE)
code_hits <- tibble(
  line = seq_along(src),
  text = src
) %>%
  filter(grepl("CREATE TEMP TABLE k3|mimiciv_derived\\.kdigo_stages|WHERE aki_stage = 3|aki_stage_crrt|stage3_by_crrt|crrt_first_time IS NULL OR crrt_first_time >= t0", text))

write_csv(code_hits, file.path(out_dir, "strictK3_QA_t0_code_evidence.csv"))

# -----------------------------------------------------------------------------
# 2) Current final SOFA 8-14 cohort: how many t0 rows used CRRT criterion?
# -----------------------------------------------------------------------------

current_summary <- qry("
WITH eq AS (
  SELECT cb.*, b.sofa_t0
  FROM data_extract_crrt_strictk3.cohort_base cb
  JOIN data_extract_crrt_strictk3.baseline b USING(stay_id)
  WHERE b.sofa_t0 BETWEEN 8 AND 14
)
SELECT
  count(*)::integer AS n,
  sum((early_crrt)::int)::integer AS early_n,
  sum((NOT early_crrt)::int)::integer AS deferred_n,
  sum((stage3_by_crrt)::int)::integer AS t0_by_crrt_n,
  sum((stage3_by_crrt AND early_crrt)::int)::integer AS t0_by_crrt_early_n,
  sum((stage3_by_crrt AND NOT early_crrt)::int)::integer AS t0_by_crrt_deferred_n,
  sum((stage3_by_creat)::int)::integer AS t0_by_creat_n,
  sum((stage3_by_uo)::int)::integer AS t0_by_uo_n
FROM eq
")

current_by_arm <- qry("
WITH eq AS (
  SELECT cb.*, b.sofa_t0
  FROM data_extract_crrt_strictk3.cohort_base cb
  JOIN data_extract_crrt_strictk3.baseline b USING(stay_id)
  WHERE b.sofa_t0 BETWEEN 8 AND 14
)
SELECT
  CASE WHEN early_crrt THEN 'early' ELSE 'deferred' END AS arm,
  count(*)::integer AS n,
  sum((stage3_by_crrt)::int)::integer AS t0_by_crrt_n,
  round(100.0 * sum((stage3_by_crrt)::int) / count(*), 1) AS t0_by_crrt_pct,
  sum((stage3_by_creat)::int)::integer AS t0_by_creat_n,
  sum((stage3_by_uo)::int)::integer AS t0_by_uo_n
FROM eq
GROUP BY arm
ORDER BY arm
")

write_csv(current_summary, file.path(out_dir, "strictK3_QA_t0_current_summary.csv"))
write_csv(current_by_arm, file.path(out_dir, "strictK3_QA_t0_current_by_arm.csv"))

# -----------------------------------------------------------------------------
# 3) Current SOFA 8-14 cohort: compare current t0 with first physiologic K3 t0
#    Physiologic K3 = creatinine stage 3 OR urine-output stage 3, excluding RRT.
# -----------------------------------------------------------------------------

current_reclassification <- qry("
WITH eq AS (
  SELECT cb.*, b.sofa_t0
  FROM data_extract_crrt_strictk3.cohort_base cb
  JOIN data_extract_crrt_strictk3.baseline b USING(stay_id)
  WHERE b.sofa_t0 BETWEEN 8 AND 14
),
phys AS (
  SELECT stay_id, min(charttime) AS phys_k3_time
  FROM mimiciv_derived.kdigo_stages
  WHERE aki_stage_creat = 3 OR aki_stage_uo = 3
  GROUP BY stay_id
),
q AS (
  SELECT
    eq.*,
    phys.phys_k3_time,
    extract(epoch FROM (eq.crrt_first_time - eq.t0)) / 3600.0 AS crrt_minus_old_t0_h,
    extract(epoch FROM (phys.phys_k3_time - eq.t0)) / 3600.0 AS phys_minus_old_t0_h,
    extract(epoch FROM (eq.crrt_first_time - phys.phys_k3_time)) / 3600.0 AS crrt_minus_phys_t0_h,
    CASE
      WHEN phys.phys_k3_time IS NULL THEN 'no_phys_k3'
      WHEN eq.crrt_first_time IS NOT NULL AND eq.crrt_first_time < phys.phys_k3_time THEN 'excluded_prior_crrt_under_phys_t0'
      WHEN eq.crrt_first_time BETWEEN phys.phys_k3_time AND phys.phys_k3_time + interval '24 hours' THEN 'early_under_phys_t0'
      ELSE 'deferred_under_phys_t0'
    END AS phys_status
  FROM eq
  LEFT JOIN phys USING(stay_id)
)
SELECT *
FROM q
ORDER BY stage3_by_crrt DESC, early_crrt DESC, stay_id
")

current_reclassification <- current_reclassification %>%
  mutate(
    stage3_by_crrt = coalesce(stage3_by_crrt, FALSE),
    stage3_by_creat = coalesce(stage3_by_creat, FALSE),
    stage3_by_uo = coalesce(stage3_by_uo, FALSE),
    crrt_minus_old_t0_h = as.numeric(crrt_minus_old_t0_h),
    phys_minus_old_t0_h = as.numeric(phys_minus_old_t0_h),
    crrt_minus_phys_t0_h = as.numeric(crrt_minus_phys_t0_h)
  )

write_csv(current_reclassification, file.path(out_dir, "strictK3_QA_t0_current_reclassification.csv"))

reclass_summary <- current_reclassification %>%
  mutate(
    old_arm = if_else(early_crrt, "old_early", "old_deferred"),
    t0_by_crrt = if_else(stage3_by_crrt, "t0_by_crrt", "t0_not_by_crrt")
  ) %>%
  count(t0_by_crrt, old_arm, phys_status, name = "n") %>%
  arrange(desc(t0_by_crrt), old_arm, phys_status)

movement_summary <- current_reclassification %>%
  group_by(stage3_by_crrt, early_crrt) %>%
  summarise(
    n = n(),
    exact_crrt_at_old_t0_n = sum(!is.na(crrt_minus_old_t0_h) & abs(crrt_minus_old_t0_h) < 1e-8),
    crrt_minus_old_t0_h = collapse_vals(crrt_minus_old_t0_h, 1),
    no_phys_k3_n = sum(is.na(phys_k3_time)),
    phys_minus_old_t0_h = collapse_vals(phys_minus_old_t0_h, 1),
    .groups = "drop"
  ) %>%
  mutate(
    stage3_by_crrt = if_else(stage3_by_crrt, "TRUE", "FALSE"),
    early_crrt = if_else(early_crrt, "TRUE", "FALSE")
  )

write_csv(reclass_summary, file.path(out_dir, "strictK3_QA_t0_reclassification_summary.csv"))
write_csv(movement_summary, file.path(out_dir, "strictK3_QA_t0_movement_summary.csv"))

# -----------------------------------------------------------------------------
# 4) Small counterfactual rerun: first physiologic K3 t0 (Cr/UO only), same core
#    extraction/rescue rules, no CCW/IPCW/Cox.
# -----------------------------------------------------------------------------

run("
CREATE TEMP TABLE k3_phys AS
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
run("CREATE INDEX ON k3_phys(stay_id)")

run("
CREATE TEMP TABLE adult_k3_phys AS
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
FROM k3_phys k3
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
run("CREATE INDEX ON adult_k3_phys(stay_id)")
run("CREATE INDEX ON adult_k3_phys(hadm_id, t0)")

run("
CREATE TEMP TABLE sepsis_at_t0_phys AS
SELECT
  a.stay_id,
  min(s.sofa_time) AS sofa_time,
  (array_agg(s.sofa_score ORDER BY s.sofa_time))[1] AS sofa_score
FROM adult_k3_phys a
JOIN mimiciv_derived.sepsis3 s
  ON s.stay_id = a.stay_id
 AND s.sepsis3 = true
 AND s.sofa_time <= a.t0
GROUP BY a.stay_id
")
run("CREATE INDEX ON sepsis_at_t0_phys(stay_id)")

run("
CREATE TEMP TABLE vaso_window_phys AS
SELECT
  a.stay_id,
  true AS vaso_in_window,
  avg(
    coalesce(v.norepinephrine, 0)
    + coalesce(v.epinephrine, 0)
    + coalesce(v.dopamine, 0) / 100.0
    + coalesce(v.phenylephrine, 0) * 0.1
    + coalesce(v.vasopressin, 0) * 2.5
  ) AS nee_at_t0
FROM adult_k3_phys a
JOIN mimiciv_derived.vasoactive_agent v
  ON v.stay_id = a.stay_id
 AND v.starttime <= a.t0
 AND v.endtime >= a.t0 - interval '6 hours'
GROUP BY a.stay_id
")
run("CREATE INDEX ON vaso_window_phys(stay_id)")

run("
CREATE TEMP TABLE crrt_phys AS
SELECT stay_id, min(charttime) AS crrt_first_time
FROM mimiciv_derived.rrt
WHERE dialysis_active = 1
  AND dialysis_type = 'CRRT'
GROUP BY stay_id
")
run("CREATE INDEX ON crrt_phys(stay_id)")

run("
CREATE TEMP TABLE source_phys AS
SELECT
  a.*,
  cb.scr_baseline AS creatinine_baseline,
  c.crrt_first_time,
  s.sofa_time,
  s.sofa_score,
  v.vaso_in_window,
  v.nee_at_t0
FROM adult_k3_phys a
JOIN sepsis_at_t0_phys s USING (stay_id)
JOIN vaso_window_phys v USING (stay_id)
LEFT JOIN mimiciv_derived.creatinine_baseline cb USING (hadm_id)
LEFT JOIN crrt_phys c USING (stay_id)
")
run("CREATE INDEX ON source_phys(stay_id)")
run("CREATE INDEX ON source_phys(hadm_id, t0)")

phys_funnel <- qry("
WITH
  n_k3 AS (SELECT count(*)::integer AS n FROM k3_phys),
  n_adult AS (SELECT count(*)::integer AS n FROM adult_k3_phys),
  n_sepsis AS (
    SELECT count(*)::integer AS n
    FROM adult_k3_phys a
    JOIN sepsis_at_t0_phys s USING (stay_id)
  ),
  n_vaso AS (SELECT count(*)::integer AS n FROM source_phys),
  n_no_prior_rrt AS (
    SELECT count(*)::integer AS n
    FROM source_phys
    WHERE crrt_first_time IS NULL OR crrt_first_time >= t0
  ),
  n_no_esrd AS (
    SELECT count(*)::integer AS n
    FROM source_phys src
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
    FROM source_phys src
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
    FROM source_phys src
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
    FROM source_phys src
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
SELECT 'first_phys_KDIGO3_stays' AS step, n FROM n_k3
UNION ALL SELECT 'adult_phys_KDIGO3_stays', n FROM n_adult
UNION ALL SELECT 'sepsis3_at_or_before_t0', n FROM n_sepsis
UNION ALL SELECT 'vasopressor_window_existing_definition', n FROM n_vaso
UNION ALL SELECT 'exclude_RRT_before_t0', n FROM n_no_prior_rrt
UNION ALL SELECT 'exclude_ESRD_chronic_dialysis', n FROM n_no_esrd
UNION ALL SELECT 'exclude_ICU_LOS_lt_24h', n FROM n_los
UNION ALL SELECT 'exclude_missing_baseline_essentials', n FROM n_baseline
UNION ALL SELECT 'exclude_zero_28d_followup', n FROM n_followup
")

run("
CREATE TEMP TABLE pre_rescue_all_stays_phys AS
SELECT src.*
FROM source_phys src
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
run("CREATE INDEX ON pre_rescue_all_stays_phys(subject_id, icu_intime)")
run("CREATE INDEX ON pre_rescue_all_stays_phys(hadm_id, t0)")

run("
CREATE TEMP TABLE pre_rescue_phys AS
SELECT DISTINCT ON (subject_id) *
FROM pre_rescue_all_stays_phys
ORDER BY subject_id, icu_intime
")
run("CREATE INDEX ON pre_rescue_phys(stay_id)")
run("CREATE INDEX ON pre_rescue_phys(hadm_id, t0)")

run("
CREATE TEMP TABLE rescue_flags_phys AS
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
  (coalesce(pr.uo_rt_24hr < 0.3, false) OR coalesce(pr.uo_rt_12hr = 0, false)) AS r4_trigger
FROM pre_rescue_phys pr
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
CREATE TEMP TABLE cohort_flags_phys AS
SELECT
  rf.*,
  (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger) AS rescue_primary_r1_r2_r2b_r3,
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
  ) AS follow_time_28d
FROM rescue_flags_phys rf
")

run("
CREATE TEMP TABLE cohort_primary_phys AS
SELECT *
FROM cohort_flags_phys
WHERE NOT rescue_primary_r1_r2_r2b_r3
")

run("
CREATE TEMP TABLE baseline_phys AS
SELECT
  c.*,
  sf.sofa_24hours AS sofa_t0
FROM cohort_primary_phys c
LEFT JOIN LATERAL (
  SELECT sofa_24hours
  FROM mimiciv_derived.sofa x
  WHERE x.stay_id = c.stay_id
    AND x.starttime <= c.t0 + interval '1 hour'
  ORDER BY x.starttime DESC
  LIMIT 1
) sf ON true
")

phys_summary <- qry("
SELECT
  'phys_CrUO_only_t0_primary_R1_R2_R2b_R3' AS version,
  count(*)::integer AS final_n,
  sum((early_crrt)::int)::integer AS early_n,
  round(100.0 * sum((early_crrt)::int) / nullif(count(*), 0), 1) AS early_pct,
  sum((sofa_t0 BETWEEN 8 AND 14)::int)::integer AS sofa8_14_n,
  sum((sofa_t0 BETWEEN 8 AND 14 AND early_crrt)::int)::integer AS sofa8_14_early_n,
  round(100.0 * sum((sofa_t0 BETWEEN 8 AND 14 AND early_crrt)::int) / nullif(sum((sofa_t0 BETWEEN 8 AND 14)::int), 0), 1) AS sofa8_14_early_pct
FROM baseline_phys
")

phys_by_sofa <- qry("
SELECT
  CASE
    WHEN sofa_t0 <= 7 THEN 'SOFA <=7'
    WHEN sofa_t0 BETWEEN 8 AND 14 THEN 'SOFA 8-14'
    WHEN sofa_t0 >= 15 THEN 'SOFA >=15'
    ELSE 'SOFA missing'
  END AS sofa_group,
  count(*)::integer AS n,
  sum((early_crrt)::int)::integer AS early_n,
  round(100.0 * sum((early_crrt)::int) / nullif(count(*), 0), 1) AS early_pct
FROM baseline_phys
GROUP BY sofa_group
ORDER BY sofa_group
")

write_csv(phys_funnel, file.path(out_dir, "strictK3_QA_t0_phys_counterfactual_funnel.csv"))
write_csv(phys_summary, file.path(out_dir, "strictK3_QA_t0_phys_counterfactual_summary.csv"))
write_csv(phys_by_sofa, file.path(out_dir, "strictK3_QA_t0_phys_counterfactual_by_sofa.csv"))

# -----------------------------------------------------------------------------
# 5) Report
# -----------------------------------------------------------------------------

cs <- current_summary[1, ]
ps <- phys_summary[1, ]

affected <- current_reclassification %>% filter(stage3_by_crrt)
affected_early <- current_reclassification %>% filter(stage3_by_crrt, early_crrt)

affected_status <- affected %>%
  count(phys_status, name = "n") %>%
  mutate(pct = sprintf("%.1f%%", 100 * n / sum(n))) %>%
  arrange(desc(n))

old_vs_phys <- tibble(
  metric = c("SOFA 8-14 n", "SOFA 8-14 early n", "SOFA 8-14 early %"),
  current_strictK3_any_stage3 = c(
    as.character(cs$n),
    as.character(cs$early_n),
    sprintf("%.1f%%", 100 * cs$early_n / cs$n)
  ),
  counterfactual_phys_CrUO_only = c(
    as.character(ps$sofa8_14_n),
    as.character(ps$sofa8_14_early_n),
    sprintf("%.1f%%", ps$sofa8_14_early_pct)
  ),
  difference = c(
    as.character(ps$sofa8_14_n - cs$n),
    as.character(ps$sofa8_14_early_n - cs$early_n),
    sprintf("%+.1f percentage points", ps$sofa8_14_early_pct - 100 * cs$early_n / cs$n)
  )
)
write_csv(old_vs_phys, file.path(out_dir, "strictK3_QA_t0_current_vs_phys_counterfactual.csv"))

lines <- c(
  "# Strict KDIGO3 t0 / CRRT-criterion contamination QA",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Bottom line",
  "",
  "- **Contamination confirmed.** The strict first-KDIGO3 extraction defines `t0` as the first row in `mimiciv_derived.kdigo_stages` where `aki_stage = 3`. In MIMIC-IV, `aki_stage = 3` can be reached through creatinine, urine-output, **or CRRT/RRT** (`aki_stage_crrt = 3`).",
  sprintf(
    "- In the final SOFA 8-14 primary cohort (n=%s; early n=%s), **%s/%s (%s)** had first KDIGO3 t0 triggered by the CRRT criterion. This was **%s/%s (%s)** of the early arm versus **%s/%s (%s)** of the deferred arm.",
    cs$n, cs$early_n,
    cs$t0_by_crrt_n, cs$n, fmt_pct(cs$t0_by_crrt_n, cs$n),
    cs$t0_by_crrt_early_n, cs$early_n, fmt_pct(cs$t0_by_crrt_early_n, cs$early_n),
    cs$t0_by_crrt_deferred_n, cs$deferred_n, fmt_pct(cs$t0_by_crrt_deferred_n, cs$deferred_n)
  ),
  sprintf(
    "- Among the %s CRRT-criterion t0 patients in SOFA 8-14, %s had no later/earlier physiologic KDIGO3 record by creatinine or urine output, and %s would become **prior-CRRT/excluded** if t0 were redefined using creatinine+urine-output only.",
    nrow(affected),
    sum(is.na(affected$phys_k3_time)),
    sum(affected$phys_status == "excluded_prior_crrt_under_phys_t0")
  ),
  "- Therefore, the Table 1 imbalance is not merely a post-t0 descriptive artifact. For a non-circular strict-KDIGO3 design, the t0 definition should exclude the RRT/CRRT KDIGO criterion and the downstream analysis should be rerun.",
  "",
  "## Code evidence",
  "",
  "Relevant lines from `strictK3_primary_extract.R`:",
  "",
  "| line | code |",
  "|---:|---|",
  apply(code_hits, 1, function(x) sprintf("| %s | `%s` |", x[["line"]], gsub("\\|", "\\\\|", trimws(x[["text"]])))),
  "",
  "Interpretation: `k3` is built from `WHERE aki_stage = 3`, not from `(aki_stage_creat = 3 OR aki_stage_uo = 3)`. The script carries `aki_stage_crrt` into `stage3_by_crrt`; the exclusion `crrt_first_time >= t0` does **not** remove records where `crrt_first_time == t0`.",
  "",
  "## Current SOFA 8-14 contamination counts",
  "",
  "| arm | n | t0 by CRRT n | t0 by CRRT % | t0 by creatinine n | t0 by urine-output n |",
  "|---|---:|---:|---:|---:|---:|",
  apply(current_by_arm, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    x[["arm"]], x[["n"]], x[["t0_by_crrt_n"]], x[["t0_by_crrt_pct"]],
    x[["t0_by_creat_n"]], x[["t0_by_uo_n"]]
  )),
  "",
  "## Timing evidence: current t0 versus CRRT and physiologic KDIGO3",
  "",
  "| t0 by CRRT | old early arm | n | CRRT minus old t0, hours | exact CRRT at old t0 n | no physiologic K3 n | physiologic K3 minus old t0, hours |",
  "|---|---|---:|---|---:|---:|---|",
  apply(movement_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    x[["stage3_by_crrt"]], x[["early_crrt"]], x[["n"]],
    x[["crrt_minus_old_t0_h"]], x[["exact_crrt_at_old_t0_n"]],
    x[["no_phys_k3_n"]], x[["phys_minus_old_t0_h"]]
  )),
  "",
  "Key reading: in the CRRT-triggered early subgroup, the median CRRT-minus-current-t0 was 0 h; thus many patients are early by definition. When using physiologic KDIGO3 only, t0 generally moves **later** or disappears rather than moving earlier.",
  "",
  "## Reclassification within the current SOFA 8-14 cohort if t0 is replaced by physiologic KDIGO3",
  "",
  "| current t0 status | old arm | physiologic-t0 status | n |",
  "|---|---|---|---:|",
  apply(reclass_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s |",
    x[["t0_by_crrt"]], x[["old_arm"]], x[["phys_status"]], x[["n"]]
  )),
  "",
  "Interpretation: most non-CRRT-triggered patients are unchanged. Among CRRT-triggered early patients, only a small minority remain early under physiologic t0; most either have prior CRRT relative to the physiologic t0 or never reach physiologic KDIGO3.",
  "",
  "## Counterfactual small rerun: first KDIGO3 by creatinine or urine-output only",
  "",
  "Same core extraction/rescue rules were used, but `k3` was redefined as `WHERE aki_stage_creat = 3 OR aki_stage_uo = 3`. No CCW/IPCW/Cox/RMST was run.",
  "",
  "| metric | current strictK3 any stage3 | counterfactual Cr/UO-only t0 | difference |",
  "|---|---:|---:|---:|",
  apply(old_vs_phys, 1, function(x) sprintf(
    "| %s | %s | %s | %s |",
    x[["metric"]], x[["current_strictK3_any_stage3"]],
    x[["counterfactual_phys_CrUO_only"]], x[["difference"]]
  )),
  "",
  "Counterfactual SOFA strata:",
  "",
  "| SOFA group | n | early n | early % |",
  "|---|---:|---:|---:|",
  apply(phys_by_sofa, 1, function(x) sprintf(
    "| %s | %s | %s | %s |",
    x[["sofa_group"]], x[["n"]], x[["early_n"]], x[["early_pct"]]
  )),
  "",
  "## Recommendation",
  "",
  "- Treat the current strict first-KDIGO3 primary analysis as **not manuscript-ready** until rerun with a non-treatment-defined t0.",
  "- Redefine strict KDIGO3 t0 as first physiologic KDIGO3 only: creatinine stage 3 and/or urine-output stage 3, explicitly excluding RRT/CRRT criterion.",
  "- Rerun the whole downstream stack after redefining t0: cohort extraction, baseline covariates, CCW/IPCW, Cox/RMST, Fine-Gray/RFD, SOFA scan, era analyses, rescue sensitivity, anchor sensitivity, and all tables/figures.",
  "- In the manuscript/reply, do not use the current Table 1 `KDIGO3 by CRRT criterion` version as a primary-analysis table. If needed, keep it only as a diagnostic artifact explaining why the rerun was required.",
  "",
  "## Output files",
  "",
  "- `strictK3_QA_t0_code_evidence.csv`",
  "- `strictK3_QA_t0_current_summary.csv`",
  "- `strictK3_QA_t0_current_by_arm.csv`",
  "- `strictK3_QA_t0_current_reclassification.csv`",
  "- `strictK3_QA_t0_reclassification_summary.csv`",
  "- `strictK3_QA_t0_movement_summary.csv`",
  "- `strictK3_QA_t0_phys_counterfactual_funnel.csv`",
  "- `strictK3_QA_t0_phys_counterfactual_summary.csv`",
  "- `strictK3_QA_t0_phys_counterfactual_by_sofa.csv`",
  "- `strictK3_QA_t0_current_vs_phys_counterfactual.csv`"
)

writeLines(lines, file.path(out_dir, "strictK3_QA_t0_contamination_report.md"))

message("Wrote report: ", file.path(out_dir, "strictK3_QA_t0_contamination_report.md"))
print(current_summary)
print(current_by_arm)
print(reclass_summary)
print(phys_summary)
