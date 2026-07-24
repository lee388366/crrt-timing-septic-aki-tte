-- =============================================================================
-- Step 6: Secondary Outcomes Extraction
-- Study: Early vs Deferred CRRT in Septic Shock + KDIGO Stage 3 AKI
-- Data:  MIMIC-IV v3.1 | Source: data_extract_crrt.cohort_base (N=3,112)
--
-- Secondary Outcome 1: 28-day Renal Recovery
--   Definition (per methods_section.docx, Line 30):
--     Survival to day 28 WITH (serum creatinine ≤ 1.5× pre-AKI baseline
--     value) AND (no ongoing RRT requirement at day 28)
--   Operationalization:
--     a. Survived to day 28: died_28d = FALSE
--     b. SCr at day 28: last available creatinine in [t₀+21d, t₀+28d]
--        from mimiciv_derived.chemistry, keyed on hadm_id
--     c. No ongoing RRT: no CRRT record (dialysis_active=1, dialysis_type='CRRT')
--        in [t₀+27d, t₀+28d] from mimiciv_derived.rrt, keyed on stay_id
--     Renal recovery = TRUE iff (a AND b_value ≤ 1.5×baseline AND NOT c)
--     → NULL when died_28d=FALSE but scr_day28 is missing (informatively censored)
--
-- Secondary Outcome 2: RRT-free days alive at day 28
--   Definition (per methods_section.docx, Line 30 & 70):
--     For patients who died within 28 days: = 0
--     For survivors: = 28 − total number of days on CRRT within [t₀, t₀+28d]
--   Operationalization:
--     Total CRRT days = count of distinct calendar dates in [t₀, t₀+28d]
--     where at least one charttime record with dialysis_active=1, dialysis_type='CRRT'
--     exists in mimiciv_derived.rrt (point-in-time charting, not continuous intervals)
-- =============================================================================

-- Drop and recreate secondary outcomes table
DROP TABLE IF EXISTS data_extract_crrt.secondary_outcomes;

CREATE TABLE data_extract_crrt.secondary_outcomes AS
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
    FROM data_extract_crrt.cohort_base
),

-- ------------------------------------------------------------------
-- CTE A: Last available SCr in [t₀+21d, t₀+28d]
--   → Proxy for "SCr at day 28" (within 7-day look-back window)
--   Keyed on hadm_id (chemistry table has no stay_id)
-- ------------------------------------------------------------------
scr_day28 AS (
    SELECT DISTINCT ON (c.stay_id)
        c.stay_id,
        ch.creatinine                          AS scr_day28,
        ch.charttime                           AS scr_day28_time
    FROM cohort c
    JOIN mimiciv_derived.chemistry ch
        ON ch.hadm_id = c.hadm_id
       AND ch.creatinine IS NOT NULL
       AND ch.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '21 days'
       AND ch.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
    ORDER BY c.stay_id, ch.charttime DESC   -- take the LAST value closest to day 28
),

-- ------------------------------------------------------------------
-- CTE B: Any CRRT charted on day 28 window [t₀+27d, t₀+28d]
--   → "Ongoing RRT requirement at day 28"
-- ------------------------------------------------------------------
rrt_at_day28 AS (
    SELECT DISTINCT c.stay_id,
        TRUE AS rrt_ongoing_day28
    FROM cohort c
    JOIN mimiciv_derived.rrt r USING (stay_id)
    WHERE r.dialysis_active = 1
      AND r.dialysis_type   = 'CRRT'
      AND r.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '27 days'
      AND r.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
),

-- ------------------------------------------------------------------
-- CTE C: Count distinct CRRT days in [t₀, t₀+28d]
--   → Used for RRT-free days alive calculation
-- ------------------------------------------------------------------
crrt_day_count AS (
    SELECT
        c.stay_id,
        COUNT(DISTINCT DATE(r.charttime)) AS crrt_total_days
    FROM cohort c
    JOIN mimiciv_derived.rrt r USING (stay_id)
    WHERE r.dialysis_active = 1
      AND r.dialysis_type   = 'CRRT'
      AND r.charttime >= (c.t0 AT TIME ZONE 'UTC')::timestamp
      AND r.charttime <= (c.t0 AT TIME ZONE 'UTC')::timestamp + INTERVAL '28 days'
    GROUP BY c.stay_id
)

-- ------------------------------------------------------------------
-- Final: Combine all components
-- ------------------------------------------------------------------
SELECT
    c.stay_id,
    c.hadm_id,
    c.subject_id,
    c.early_crrt,
    c.died_28d,
    c.follow_time_28d,

    -- ── Secondary Outcome 1: 28-day Renal Recovery ──────────────────
    s.scr_day28,
    s.scr_day28_time,
    c.creatinine_baseline,
    -- threshold: 1.5 × baseline
    ROUND((1.5 * c.creatinine_baseline)::numeric, 3)  AS scr_recovery_threshold,
    -- sub-components for transparency
    COALESCE(r.rrt_ongoing_day28, FALSE)               AS rrt_ongoing_day28,
    -- final binary flag
    CASE
        WHEN c.died_28d = TRUE   THEN FALSE        -- must survive to day 28
        WHEN s.scr_day28 IS NULL THEN NULL          -- missing SCr → informative missingness, flag as NULL
        WHEN s.scr_day28 > 1.5 * c.creatinine_baseline THEN FALSE  -- creatinine criterion not met
        WHEN COALESCE(r.rrt_ongoing_day28, FALSE) = TRUE THEN FALSE -- still on RRT
        ELSE TRUE                                   -- all criteria met → recovered
    END                                             AS renal_recovery_28d,

    -- ── Secondary Outcome 2: RRT-free days alive at day 28 ──────────
    COALESCE(d.crrt_total_days, 0)                  AS crrt_total_days,
    CASE
        WHEN c.died_28d = TRUE
            THEN 0.0                                -- died → 0 (per protocol)
        ELSE GREATEST(0.0,
             28.0 - COALESCE(d.crrt_total_days, 0))-- survived → 28 minus CRRT days
    END                                             AS rrt_free_days_28d

FROM cohort c
LEFT JOIN scr_day28       s USING (stay_id)
LEFT JOIN rrt_at_day28    r USING (stay_id)
LEFT JOIN crrt_day_count  d USING (stay_id);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_sec_outcomes_stay
    ON data_extract_crrt.secondary_outcomes (stay_id);

ANALYZE data_extract_crrt.secondary_outcomes;

-- ------------------------------------------------------------------
-- Quick sanity check
-- ------------------------------------------------------------------
SELECT
    COUNT(*)                                                      AS n_total,
    SUM(CASE WHEN died_28d THEN 1 ELSE 0 END)                   AS n_died_28d,
    SUM(CASE WHEN renal_recovery_28d = TRUE  THEN 1 ELSE 0 END) AS n_renal_recovery,
    SUM(CASE WHEN renal_recovery_28d IS NULL THEN 1 ELSE 0 END) AS n_missing_scr_day28,
    ROUND(AVG(rrt_free_days_28d)::numeric, 1)                   AS mean_rrt_free_days,
    PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY rrt_free_days_28d)
                                                                  AS median_rrt_free_days,
    ROUND(100.0 * SUM(CASE WHEN renal_recovery_28d = TRUE THEN 1 ELSE 0 END)
          / COUNT(*)::numeric, 1)                                 AS pct_renal_recovery
FROM data_extract_crrt.secondary_outcomes;
