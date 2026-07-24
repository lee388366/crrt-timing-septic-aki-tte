-- =============================================================================
-- strictK3 v2 — Figure 1 funnel verification (clean Cr/UO-only t0)
-- Source of truth (logic): analysis/strictK3_primary_extract_v2.R
-- DB: mimiciv31 @ localhost
-- Purpose: reproduce funnel step counts + primary rescue cohort + SOFA 8–14
--          equipoise (sofa_t0 from mimiciv_derived.sofa, same rule as
--          strictK3_baseline_covariates_v2.R: last sofa_24hours with
--          starttime <= t0 + 1 hour).
-- Does NOT write patient-level tables; TEMP only.
-- =============================================================================

\timing on
\pset footer off

-- -----------------------------------------------------------------------------
-- 1) Clean physiologic KDIGO3 t0 (Cr or UO stage 3; RRT criterion NOT used)
-- -----------------------------------------------------------------------------
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
ORDER BY stay_id, charttime;
CREATE INDEX ON k3(stay_id);

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
  k3.kdigo3_time, k3.t0;
CREATE INDEX ON adult_k3(stay_id);
CREATE INDEX ON adult_k3(hadm_id, t0);

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
GROUP BY a.stay_id;
CREATE INDEX ON sepsis_at_t0(stay_id);

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
GROUP BY a.stay_id;
CREATE INDEX ON vaso_window(stay_id);

CREATE TEMP TABLE crrt AS
SELECT stay_id, min(charttime) AS crrt_first_time
FROM mimiciv_derived.rrt
WHERE dialysis_active = 1
  AND dialysis_type = 'CRRT'
GROUP BY stay_id;
CREATE INDEX ON crrt(stay_id);

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
LEFT JOIN crrt c USING (stay_id);
CREATE INDEX ON source(stay_id);
CREATE INDEX ON source(hadm_id, t0);

-- -----------------------------------------------------------------------------
-- 2) Funnel steps through follow-up (matches extract_v2 source_counts CTE)
-- -----------------------------------------------------------------------------
CREATE TEMP TABLE funnel_core AS
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
SELECT * FROM (
  SELECT 1 AS ord, 'first_KDIGO3_stays' AS step, n FROM n_k3
  UNION ALL SELECT 2, 'adult_KDIGO3_stays', n FROM n_adult
  UNION ALL SELECT 3, 'sepsis3_at_or_before_t0', n FROM n_sepsis
  UNION ALL SELECT 4, 'vasopressor_window_existing_definition', n FROM n_vaso
  UNION ALL SELECT 5, 'exclude_RRT_before_t0', n FROM n_no_prior_rrt
  UNION ALL SELECT 6, 'exclude_ESRD_chronic_dialysis', n FROM n_no_esrd
  UNION ALL SELECT 7, 'exclude_ICU_LOS_lt_24h', n FROM n_los
  UNION ALL SELECT 8, 'exclude_missing_baseline_essentials', n FROM n_baseline
  UNION ALL SELECT 9, 'exclude_zero_28d_followup', n FROM n_followup
) s
ORDER BY ord;

-- -----------------------------------------------------------------------------
-- 3) Pre-rescue (first eligible ICU per subject) + primary rescue exclusion
-- -----------------------------------------------------------------------------
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
  ) > 0;
CREATE INDEX ON pre_rescue_all_stays(subject_id, icu_intime);
CREATE INDEX ON pre_rescue_all_stays(hadm_id, t0);

CREATE TEMP TABLE pre_rescue AS
SELECT DISTINCT ON (subject_id) *
FROM pre_rescue_all_stays
ORDER BY subject_id, icu_intime;
CREATE INDEX ON pre_rescue(stay_id);
CREATE INDEX ON pre_rescue(hadm_id, t0);

CREATE TEMP TABLE rescue_flags AS
SELECT
  pr.*,
  coalesce(kstats.r1_trigger, false) AS r1_trigger,
  coalesce(abg.r2_trigger, false) AS r2_trigger,
  coalesce(hco3.r2b_trigger, false) AS r2b_trigger,
  coalesce(bun.r3_trigger, false) AS r3_trigger,
  (coalesce(pr.uo_rt_24hr < 0.3, false) OR coalesce(pr.uo_rt_12hr = 0, false)) AS r4_trigger
FROM pre_rescue pr
LEFT JOIN LATERAL (
  SELECT bool_or(c.potassium >= 6.5) AS r1_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '6 hours'
    AND c.charttime <= pr.t0
    AND c.potassium IS NOT NULL
) kstats ON true
LEFT JOIN LATERAL (
  SELECT bool_or(bg.ph <= 7.15) AS r2_trigger
  FROM mimiciv_derived.bg bg
  WHERE bg.hadm_id = pr.hadm_id
    AND bg.charttime > pr.t0 - interval '6 hours'
    AND bg.charttime <= pr.t0
    AND bg.specimen = 'ART.'
    AND bg.ph IS NOT NULL
) abg ON true
LEFT JOIN LATERAL (
  SELECT bool_or(c.bicarbonate <= 12) AS r2b_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '6 hours'
    AND c.charttime <= pr.t0
    AND c.bicarbonate IS NOT NULL
) hco3 ON true
LEFT JOIN LATERAL (
  SELECT bool_or(c.bun >= 100) AS r3_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '24 hours'
    AND c.charttime <= pr.t0
    AND c.bun IS NOT NULL
) bun ON true;

CREATE TEMP TABLE cohort_primary AS
SELECT
  rf.*,
  (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger) AS rescue_primary,
  coalesce(crrt_first_time BETWEEN t0 AND t0 + interval '24 hours', false) AS early_crrt
FROM rescue_flags rf
WHERE NOT (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger);
CREATE INDEX ON cohort_primary(stay_id);

-- -----------------------------------------------------------------------------
-- 4) SOFA at t0 (equipoise) — same rule as baseline_covariates_v2.R
-- -----------------------------------------------------------------------------
CREATE TEMP TABLE sofa_t0 AS
SELECT DISTINCT ON (c.stay_id)
  c.stay_id,
  s.sofa_24hours AS sofa_t0
FROM cohort_primary c
JOIN mimiciv_derived.sofa s
  ON s.stay_id = c.stay_id
 AND s.starttime <= c.t0 + interval '1 hour'
ORDER BY c.stay_id, s.starttime DESC;
CREATE INDEX ON sofa_t0(stay_id);

-- -----------------------------------------------------------------------------
-- OUTPUT A: funnel core
-- -----------------------------------------------------------------------------
\echo '=== FUNNEL_CORE ==='
SELECT step, n FROM funnel_core ORDER BY ord;

-- -----------------------------------------------------------------------------
-- OUTPUT B: rescue funnel + early
-- -----------------------------------------------------------------------------
\echo '=== RESCUE_AND_PRIMARY ==='
SELECT 'pre_rescue_first_eligible_ICU' AS step, count(*)::int AS n FROM pre_rescue
UNION ALL
SELECT 'exclude_primary_rescue_R1_R2_R2b_R3', count(*)::int FROM cohort_primary
UNION ALL
SELECT 'primary_early_crrt', count(*)::int FROM cohort_primary WHERE early_crrt
UNION ALL
SELECT 'primary_deferred', count(*)::int FROM cohort_primary WHERE NOT early_crrt;

-- -----------------------------------------------------------------------------
-- OUTPUT C: SOFA 8-14 equipoise (Fig1 bottom box)
-- -----------------------------------------------------------------------------
\echo '=== SOFA_EQUIPOISE ==='
SELECT
  count(*)::int AS n_sofa_8_14,
  count(*) FILTER (WHERE early_crrt)::int AS early_n,
  count(*) FILTER (WHERE NOT early_crrt)::int AS deferred_n
FROM cohort_primary c
JOIN sofa_t0 s USING (stay_id)
WHERE s.sofa_t0 >= 8 AND s.sofa_t0 <= 14;
