-- =============================================================================
-- Strict first-KDIGO3 secondary outcomes extraction
-- Study: Early vs Deferred CRRT in septic shock + strict first KDIGO stage 3 AKI
-- Data:  MIMIC-IV v3.1
-- Source cohort: data_extract_crrt_strictk3.cohort_base
--
-- IMPORTANT:
--   This script intentionally does NOT reuse analysis/data/analysis_master.csv,
--   because that file belongs to the legacy KDIGO3-minus-6h anchor. All windows
--   below are recalculated from the strict first-KDIGO3 t0.
--
-- Secondary Outcome 1: 28-day renal recovery
--   Survival to day 28 AND last creatinine in [t0+21d, t0+28d]
--   <= 1.5 x pre-AKI baseline AND no ongoing CRRT in [t0+27d, t0+28d].
--   If the patient survived but no day-21-to-day-28 creatinine is available,
--   renal_recovery_28d is NULL.
--
-- Secondary Outcome 2: RRT-free days alive at day 28
--   If died within 28 days: 0.
--   If survived: 28 - number of distinct calendar dates with active CRRT
--   records in [t0, t0+28d].
-- =============================================================================

CREATE SCHEMA IF NOT EXISTS data_extract_crrt_strictk3;

DROP TABLE IF EXISTS data_extract_crrt_strictk3.secondary_outcomes;

CREATE TABLE data_extract_crrt_strictk3.secondary_outcomes AS
WITH cohort AS (
    SELECT
        stay_id,
        hadm_id,
        subject_id,
        t0,
        creatinine_baseline,
        died_28d,
        follow_time_28d,
        early_crrt
    FROM data_extract_crrt_strictk3.cohort_base
),

scr_day28 AS (
    SELECT DISTINCT ON (c.stay_id)
        c.stay_id,
        ch.creatinine AS scr_day28,
        ch.charttime  AS scr_day28_time
    FROM cohort c
    JOIN mimiciv_derived.chemistry ch
      ON ch.hadm_id = c.hadm_id
     AND ch.creatinine IS NOT NULL
     AND ch.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '21 days'
     AND ch.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
    ORDER BY c.stay_id, ch.charttime DESC
),

rrt_at_day28 AS (
    SELECT DISTINCT
        c.stay_id,
        TRUE AS rrt_ongoing_day28
    FROM cohort c
    JOIN mimiciv_derived.rrt r USING (stay_id)
    WHERE r.dialysis_active = 1
      AND r.dialysis_type = 'CRRT'
      AND r.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '27 days'
      AND r.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
),

crrt_day_count AS (
    SELECT
        c.stay_id,
        COUNT(DISTINCT DATE(r.charttime)) AS crrt_total_days
    FROM cohort c
    JOIN mimiciv_derived.rrt r USING (stay_id)
    WHERE r.dialysis_active = 1
      AND r.dialysis_type = 'CRRT'
      AND r.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp
      AND r.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
    GROUP BY c.stay_id
)

SELECT
    c.stay_id,
    c.hadm_id,
    c.subject_id,
    c.early_crrt,
    c.died_28d,
    c.follow_time_28d,
    s.scr_day28,
    s.scr_day28_time,
    c.creatinine_baseline,
    ROUND((1.5 * c.creatinine_baseline)::numeric, 3) AS scr_recovery_threshold,
    COALESCE(r.rrt_ongoing_day28, FALSE) AS rrt_ongoing_day28,
    CASE
        WHEN c.died_28d = TRUE THEN FALSE
        WHEN s.scr_day28 IS NULL THEN NULL
        WHEN s.scr_day28 > 1.5 * c.creatinine_baseline THEN FALSE
        WHEN COALESCE(r.rrt_ongoing_day28, FALSE) = TRUE THEN FALSE
        ELSE TRUE
    END AS renal_recovery_28d,
    COALESCE(d.crrt_total_days, 0) AS crrt_total_days,
    CASE
        WHEN c.died_28d = TRUE THEN 0.0
        ELSE GREATEST(0.0, 28.0 - COALESCE(d.crrt_total_days, 0))
    END AS rrt_free_days_28d
FROM cohort c
LEFT JOIN scr_day28      s USING (stay_id)
LEFT JOIN rrt_at_day28   r USING (stay_id)
LEFT JOIN crrt_day_count d USING (stay_id);

CREATE INDEX IF NOT EXISTS idx_strictk3_secondary_outcomes_stay
    ON data_extract_crrt_strictk3.secondary_outcomes (stay_id);

ANALYZE data_extract_crrt_strictk3.secondary_outcomes;

-- Quick sanity check
SELECT
    COUNT(*) AS n_total,
    SUM(CASE WHEN died_28d THEN 1 ELSE 0 END) AS n_died_28d,
    SUM(CASE WHEN renal_recovery_28d = TRUE THEN 1 ELSE 0 END) AS n_renal_recovery,
    SUM(CASE WHEN renal_recovery_28d IS NULL THEN 1 ELSE 0 END) AS n_missing_scr_day28,
    ROUND(AVG(rrt_free_days_28d)::numeric, 1) AS mean_rrt_free_days,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY rrt_free_days_28d) AS median_rrt_free_days,
    ROUND(
        100.0 * SUM(CASE WHEN renal_recovery_28d = TRUE THEN 1 ELSE 0 END)
        / COUNT(*)::numeric,
        1
    ) AS pct_renal_recovery
FROM data_extract_crrt_strictk3.secondary_outcomes;

-- Navicat-friendly result view
SELECT *
FROM data_extract_crrt_strictk3.secondary_outcomes
LIMIT 100;
