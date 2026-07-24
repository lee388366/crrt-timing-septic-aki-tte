-- =============================================================================
-- strictK3 v2 — Fig2 / CCW clone-skeleton verification (SQL layer)
-- Design source: analysis/03_ccw_ipcw.R (clone + protocol censor + follow/event)
--               via strictK3_run_timevarying_ccw_v2.R → data/ccw_final.rds
--
-- What this SQL DOES:
--   1) Recompute sofa_t0 from mimiciv_derived.sofa (same rule as baseline R)
--      and compare to data_extract_crrt_strictk3_v2.baseline.sofa_t0
--   2) Build A_early / B_deferred clones with grace=24h protocol censoring
--   3) Emit summary counts for alignment vs ccw_final.rds
--
-- What this SQL does NOT do (by design — stays in R/glmnet):
--   Ridge IPCW weights (ipcw_w / ipcw_w_trim). Those require tv_hourly + glmnet.
--   Fig2 weighted n.risk uses R survfit on ccw_final.rds (see fig2_step2_*).
--
-- DB: mimiciv31 @ localhost
-- Requires: data_extract_crrt_strictk3_v2.baseline (from extract + baseline R)
-- =============================================================================

\timing on
\pset footer off

-- -----------------------------------------------------------------------------
-- 1) sofa_t0 from MIMIC derived (last sofa_24hours with starttime <= t0 + 1h)
-- -----------------------------------------------------------------------------
CREATE TEMP TABLE sofa_t0_recompute AS
SELECT DISTINCT ON (b.stay_id)
  b.stay_id,
  s.sofa_24hours AS sofa_t0_sql
FROM data_extract_crrt_strictk3_v2.baseline b
JOIN mimiciv_derived.sofa s
  ON s.stay_id = b.stay_id
 AND s.starttime <= b.t0 + interval '1 hour'
ORDER BY b.stay_id, s.starttime DESC;

CREATE TEMP TABLE sofa_t0_check AS
SELECT
  b.stay_id,
  b.sofa_t0 AS sofa_t0_baseline,
  r.sofa_t0_sql,
  (b.sofa_t0 IS NOT DISTINCT FROM r.sofa_t0_sql) AS match_sofa
FROM data_extract_crrt_strictk3_v2.baseline b
LEFT JOIN sofa_t0_recompute r USING (stay_id);

\echo '=== SOFA_T0_RECOMPUTE_VS_BASELINE ==='
SELECT
  count(*)::int AS n,
  count(*) FILTER (WHERE match_sofa)::int AS n_match,
  count(*) FILTER (WHERE NOT match_sofa OR sofa_t0_sql IS NULL)::int AS n_mismatch_or_null
FROM sofa_t0_check;

-- -----------------------------------------------------------------------------
-- 2) Clone-censor skeleton (03_ccw_ipcw.R; GRACE = 24 h)
-- -----------------------------------------------------------------------------
CREATE TEMP TABLE ccw_clones_sql AS
WITH cohort AS (
  SELECT
    b.*,
    b.crrt_hours_from_t0 AS h_crrt,
    (b.sofa_t0 >= 8 AND b.sofa_t0 <= 14) AS equipoise
  FROM data_extract_crrt_strictk3_v2.baseline b
  WHERE b.sofa_t0 IS NOT NULL
),
strategies AS (
  SELECT 'A_early'::text AS strategy
  UNION ALL
  SELECT 'B_deferred'::text
),
clones AS (
  SELECT
    c.*,
    s.strategy,
    c.stay_id::text || '_' || s.strategy AS clone_id,
    CASE
      WHEN s.strategy = 'A_early'
           AND (c.h_crrt IS NULL OR c.h_crrt > 24) THEN 24.0
      WHEN s.strategy = 'B_deferred'
           AND c.h_crrt IS NOT NULL AND c.h_crrt <= 24 THEN c.h_crrt
      ELSE NULL
    END AS protocol_censor_h
  FROM cohort c
  CROSS JOIN strategies s
)
SELECT
  stay_id,
  subject_id,
  strategy,
  clone_id,
  early_crrt,
  h_crrt,
  follow_time_28d,
  died_28d,
  sofa_t0,
  equipoise,
  protocol_censor_h,
  (protocol_censor_h IS NOT NULL) AS protocol_censored,
  LEAST(follow_time_28d * 24.0, COALESCE(protocol_censor_h, follow_time_28d * 24.0)) AS follow_h,
  LEAST(follow_time_28d * 24.0, COALESCE(protocol_censor_h, follow_time_28d * 24.0)) / 24.0 AS follow_days,
  (
    died_28d
    AND (
      protocol_censor_h IS NULL
      OR follow_time_28d * 24.0 < protocol_censor_h
    )
  ) AS event
FROM clones;

CREATE INDEX ON ccw_clones_sql(stay_id, strategy);

\echo '=== CCW_CLONE_SKELETON_SUMMARY ==='
SELECT
  count(*)::int AS n_clone_rows,
  count(DISTINCT stay_id)::int AS n_stays,
  count(*) FILTER (WHERE strategy = 'A_early')::int AS n_A,
  count(*) FILTER (WHERE strategy = 'B_deferred')::int AS n_B,
  count(*) FILTER (WHERE equipoise)::int AS n_equipoise_clones,
  count(DISTINCT stay_id) FILTER (WHERE equipoise)::int AS n_equipoise_stays,
  count(*) FILTER (WHERE strategy = 'A_early' AND protocol_censored)::int AS n_censored_A,
  count(*) FILTER (WHERE strategy = 'B_deferred' AND protocol_censored)::int AS n_censored_B,
  count(*) FILTER (WHERE equipoise AND strategy = 'A_early' AND early_crrt)::int AS n_eq_early_true_on_A
FROM ccw_clones_sql;

\echo '=== EQUIPOSE_AT_RISK_UNWEIGHTED_NOTE ==='
-- Unweighted clone counts at t0 (NOT Fig2 weighted n.risk; weights are R-only)
SELECT
  strategy,
  count(*)::int AS n_clones_equipoise,
  count(*) FILTER (WHERE NOT protocol_censored OR protocol_censor_h > 0)::int AS n_enter_followup
FROM ccw_clones_sql
WHERE equipoise
GROUP BY strategy
ORDER BY strategy;

-- -----------------------------------------------------------------------------
-- 3) Export clone keys for RDS join (run via psql \copy from shell wrapper)
-- -----------------------------------------------------------------------------
-- Columns chosen for exact compare vs ccw_final.rds (excluding ipcw weights)
\echo '=== EXPORT_HINT ==='
\echo 'Use companion shell/R to \\copy ccw_clones_sql key columns to CSV and diff vs RDS.'
