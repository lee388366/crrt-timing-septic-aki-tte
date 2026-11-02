-- Baseline characteristics query used for Table 1 and Supplemental Table S14 (comments in the original language removed;
-- the SELECT body is extracted by revision_2026_11/figures_tables/dependencies/TABLE1_READONLY_ORIGINAL_ADAPTER.R).
-- ============================================================
-- 05) cohort_baseline_v2 — Baseline characteristics (Table 1)
-- Project  : MIMIC-IV CCW — CRRT timing in Septic Shock ∩ AKI (KDIGO3)
-- Output   : data_extract_crrt.cohort_baseline_v2
-- Changes vs v1:
--
-- ─────────────────────────────────────────────────────────────
--
-- ─────────────────────────────────────────────────────────────
--
-- ─────────────────────────────────────────────────────────────
-- mimiciv_icu.icustays              → first_careunit
-- mimiciv_derived.weight_durations  → weight_kg
-- mimiciv_derived.vitalsign         → hr/sbp/dbp/mbp/rr/spo2/temperature
-- mimiciv_derived.chemistry         → na/cl/k/bicarbonate/creatinine/bun
-- mimiciv_derived.bg                → ph/lactate/po2/fio2
-- mimiciv_derived.coagulation       → inr
-- mimiciv_derived.complete_blood_count → wbc/hemoglobin/platelet
-- mimiciv_derived.gcs               → gcs_motor/gcs_verbal/gcs_eyes
-- mimiciv_derived.ventilation       → ventilation_status = 'InvasiveVent'
-- mimiciv_derived.urine_output      → urineoutput（mL/event）
-- mimiciv_derived.kdigo_stages      → aki_stage_creat / aki_stage_uo
--
-- ─────────────────────────────────────────────────────────────
-- ============================================================

DROP TABLE IF EXISTS data_extract_crrt.cohort_baseline_v2;

CREATE TABLE data_extract_crrt.cohort_baseline_v2 AS

WITH cohort AS (
  SELECT
    cb.stay_id,
    cb.subject_id,
    cb.hadm_id,
    cb.icu_intime                          AS intime,
    cb.icu_outtime                         AS outtime,
    cb.age,
    cb.gender,
    adm.race,
    adm.insurance,
    cb.icu_los_hours / 24.0                AS los_icu,
    cb.kdigo3_time,
    cb.t0                                  AS t0_time,
    cb.t0 - INTERVAL '6 hour'             AS baseline_start,
    adm.admittime
  FROM data_extract_crrt.cohort_base cb
  JOIN mimiciv_hosp.admissions adm ON adm.hadm_id = cb.hadm_id
  WHERE cb.stay_id IN (
    SELECT DISTINCT ON (subject_id) stay_id
    FROM data_extract_crrt.cohort_base
    ORDER BY subject_id, icu_intime
  )
)

SELECT
  -- ── identifiers ──────────────────────────────────────────────────────────────────
  c.subject_id,
  c.hadm_id,
  c.stay_id,
  c.intime,
  c.outtime,
  c.t0_time,

  icu_det.first_careunit,
  EXTRACT(YEAR FROM c.intime)::INT AS admission_year,

  c.age,
  c.gender,
  c.race,
  c.insurance,
  w.weight_kg,

  cc.charlson_comorbidity_index AS cci,
  sf.sofa_baseline,
  gcs.gcs_baseline,
  ap.apsiii_score,

  v.hr,
  v.sbp,
  v.dbp,
  v.map,
  v.rr,
  v.spo2,
  v.temp_c,

  cb.wbc,
  cb.hb,
  cb.plt,

  ch.na,
  ch.cl,
  ch.k,
  ch.hco3,
  ch.scr,
  ch.bun,

  enz.ast,
  enz.alt,

  b.ph,
  b.lactate,
  b.pao2,
  b.fio2,
  CASE
    WHEN b.pao2 IS NOT NULL AND b.fio2 IS NOT NULL AND b.fio2 > 0
    THEN b.pao2 / NULLIF(b.fio2 / 100.0, 0)
    ELSE NULL
  END AS pf_ratio,

  cg.inr,

  GREATEST(COALESCE(cc.diabetes_without_cc, 0),
           COALESCE(cc.diabetes_with_cc,    0)) AS dm,
  COALESCE(htn.htn, 0)                          AS htn,
  COALESCE(cc.congestive_heart_failure,     0)  AS chf,
  COALESCE(cc.renal_disease,               0)   AS ckd,
  GREATEST(COALESCE(cc.mild_liver_disease,  0),
           COALESCE(cc.severe_liver_disease, 0)) AS liver_disease,
  GREATEST(COALESCE(cc.malignant_cancer,    0),
           COALESCE(cc.metastatic_solid_tumor, 0)) AS malignancy,
  COALESCE(cc.chronic_pulmonary_disease,   0)   AS copd,

  sp.mv_baseline,
  sp.vaso_baseline,
  sp.ne_eq_baseline,

  CASE
    WHEN w.weight_kg > 0 THEN cr.crystalloid_ml / w.weight_kg
    ELSE NULL
  END AS crystalloid_ml_kg,
  co.colloid_flag,
  nt.nephrotoxic_flag,

  uo_t0.uo_mlkgph_t0_last6h,
  akip.aki_pathway,
  creat_bl.scr_baseline

FROM cohort c

LEFT JOIN mimiciv_icu.icustays icu_det
  ON icu_det.stay_id = c.stay_id

LEFT JOIN LATERAL (
  SELECT wd.weight AS weight_kg
  FROM mimiciv_derived.weight_durations wd
  WHERE wd.stay_id = c.stay_id
    AND wd.starttime <= c.t0_time
    AND wd.weight IS NOT NULL
  ORDER BY wd.starttime DESC
  LIMIT 1
) w ON true

LEFT JOIN mimiciv_derived.charlson cc
  ON cc.hadm_id = c.hadm_id

LEFT JOIN LATERAL (
  SELECT
    v.heart_rate AS hr, v.sbp, v.dbp, v.mbp AS map,
    v.resp_rate AS rr, v.spo2, v.temperature AS temp_c
  FROM mimiciv_derived.vitalsign v
  WHERE v.stay_id = c.stay_id
    AND v.charttime >  c.baseline_start
    AND v.charttime <= c.t0_time
  ORDER BY v.charttime DESC
  LIMIT 1
) v ON true

LEFT JOIN LATERAL (
  SELECT
    ch.sodium AS na, ch.chloride AS cl, ch.potassium AS k,
    ch.bicarbonate AS hco3, ch.creatinine AS scr, ch.bun
  FROM mimiciv_derived.chemistry ch
  WHERE ch.hadm_id = c.hadm_id
    AND ch.charttime >  c.baseline_start
    AND ch.charttime <= c.t0_time
  ORDER BY ch.charttime DESC
  LIMIT 1
) ch ON true

LEFT JOIN LATERAL (
  SELECT b.ph, b.lactate, b.po2 AS pao2, b.fio2
  FROM mimiciv_derived.bg b
  WHERE b.hadm_id = c.hadm_id
    AND b.charttime >  c.baseline_start
    AND b.charttime <= c.t0_time
  ORDER BY b.charttime DESC
  LIMIT 1
) b ON true

LEFT JOIN LATERAL (
  SELECT e.ast, e.alt
  FROM mimiciv_derived.enzyme e
  WHERE e.hadm_id = c.hadm_id
    AND e.charttime >  c.baseline_start
    AND e.charttime <= c.t0_time
  ORDER BY e.charttime DESC
  LIMIT 1
) enz ON true

LEFT JOIN LATERAL (
  SELECT cg.inr
  FROM mimiciv_derived.coagulation cg
  WHERE cg.hadm_id = c.hadm_id
    AND cg.charttime >  c.baseline_start
    AND cg.charttime <= c.t0_time
  ORDER BY cg.charttime DESC
  LIMIT 1
) cg ON true

LEFT JOIN LATERAL (
  SELECT cb.wbc, cb.hemoglobin AS hb, cb.platelet AS plt
  FROM mimiciv_derived.complete_blood_count cb
  WHERE cb.hadm_id = c.hadm_id
    AND cb.charttime >  c.baseline_start
    AND cb.charttime <= c.t0_time
  ORDER BY cb.charttime DESC
  LIMIT 1
) cb ON true

-- ── JOIN 11: GCS ───────────────────────────────────────────────────────────────
LEFT JOIN LATERAL (
  SELECT (g.gcs_motor + g.gcs_verbal + g.gcs_eyes) AS gcs_baseline
  FROM mimiciv_derived.gcs g
  WHERE g.stay_id = c.stay_id
    AND g.charttime >  c.baseline_start
    AND g.charttime <= c.t0_time
    AND g.gcs_motor IS NOT NULL AND g.gcs_verbal IS NOT NULL AND g.gcs_eyes IS NOT NULL
  ORDER BY g.charttime DESC
  LIMIT 1
) gcs ON true

-- ── JOIN 12: APSIII ────────────────────────────────────────────────────────────
LEFT JOIN LATERAL (
  SELECT a.apsiii AS apsiii_score
  FROM mimiciv_derived.apsiii a
  WHERE a.stay_id = c.stay_id
  LIMIT 1
) ap ON true

LEFT JOIN LATERAL (
  SELECT
    EXISTS (
      SELECT 1 FROM mimiciv_derived.ventilation vv
      WHERE vv.stay_id = c.stay_id
        AND vv.starttime < c.t0_time
        AND COALESCE(vv.endtime, c.t0_time) > c.baseline_start
        AND vv.ventilation_status = 'InvasiveVent'
    )::int AS mv_baseline,
    EXISTS (
      SELECT 1 FROM mimiciv_derived.vasoactive_agent va
      WHERE va.stay_id = c.stay_id
        AND va.starttime < c.t0_time
        AND COALESCE(va.endtime, c.t0_time) > c.baseline_start
        AND (  COALESCE(va.norepinephrine, 0) > 0
            OR COALESCE(va.epinephrine,    0) > 0
            OR COALESCE(va.phenylephrine,  0) > 0
            OR COALESCE(va.vasopressin,    0) > 0
            OR COALESCE(va.dopamine,       0) > 0)
    )::int AS vaso_baseline,
    (SELECT ne.norepinephrine_equivalent_dose
     FROM mimiciv_derived.norepinephrine_equivalent_dose ne
     WHERE ne.stay_id = c.stay_id
       AND ne.starttime < c.t0_time
       AND COALESCE(ne.endtime, c.t0_time) > c.baseline_start
       AND ne.norepinephrine_equivalent_dose IS NOT NULL
     ORDER BY ne.endtime DESC
     LIMIT 1) AS ne_eq_baseline
  FROM (SELECT 1) _
) sp ON true

LEFT JOIN LATERAL (
  SELECT s.sofa_24hours AS sofa_baseline
  FROM mimiciv_derived.sofa s
  WHERE s.stay_id = c.stay_id
    AND s.endtime >  c.baseline_start
    AND s.endtime <= c.t0_time
    AND s.sofa_24hours IS NOT NULL
  ORDER BY s.endtime DESC
  LIMIT 1
) sf ON true

LEFT JOIN LATERAL (
  SELECT EXISTS (
    SELECT 1 FROM mimiciv_hosp.diagnoses_icd d
    WHERE d.hadm_id = c.hadm_id
      AND (
        (d.icd_version = 9  AND d.icd_code ~ '^40[1-5]')
        OR
        (d.icd_version = 10 AND (d.icd_code ~ '^I1[0-3]' OR d.icd_code ~ '^I15'))
      )
  )::int AS htn
  FROM (SELECT 1) _
) htn ON true

--   225827 D5+0.5%NaCl     225941 D5+0.25%NaCl
LEFT JOIN LATERAL (
  SELECT COALESCE(SUM(ie.amount), 0) AS crystalloid_ml
  FROM mimiciv_icu.inputevents ie
  WHERE ie.stay_id = c.stay_id
    AND ie.starttime >= c.intime
    AND ie.starttime <  c.t0_time
    AND ie.itemid IN (
        225158, 225159, 225161, 225828, 225944,
        226364, 226375, 227533, 220949, 220950,
        225797, 225825, 225827, 225941
    )
    AND ie.amountuom = 'mL'
) cr ON true

LEFT JOIN LATERAL (
  SELECT EXISTS (
    SELECT 1 FROM mimiciv_icu.inputevents ie
    WHERE ie.stay_id = c.stay_id
      AND ie.starttime >= c.intime
      AND ie.starttime <  c.t0_time
      AND ie.itemid IN (225174, 225795, 225796)
  )::int AS colloid_flag
  FROM (SELECT 1) _
) co ON true

LEFT JOIN LATERAL (
  SELECT EXISTS (
    SELECT 1 FROM mimiciv_hosp.prescriptions p
    WHERE p.hadm_id = c.hadm_id
      AND p.starttime >= c.intime
      AND p.starttime <  c.t0_time
      AND p.drug ~* 'gentamicin|tobramycin|amikacin|streptomycin|vancomycin|ketorolac|ibuprofen|indomethacin|naproxen|amphotericin|colistin|polymyxin'
  )::int AS nephrotoxic_flag
  FROM (SELECT 1) _
) nt ON true

LEFT JOIN LATERAL (
  SELECT
    COALESCE(SUM(uo.urineoutput), 0)
      / NULLIF((
          SELECT wd.weight
          FROM mimiciv_derived.weight_durations wd
          WHERE wd.stay_id = c.stay_id
            AND wd.starttime <= c.t0_time
            AND COALESCE(wd.endtime, c.t0_time) >= c.t0_time
            AND wd.weight IS NOT NULL
          ORDER BY wd.starttime DESC
          LIMIT 1
        ), 0)
      / 6.0 AS uo_mlkgph_t0_last6h
  FROM mimiciv_derived.urine_output uo
  WHERE uo.stay_id = c.stay_id
    AND uo.charttime >  c.baseline_start
    AND uo.charttime <= c.t0_time
) uo_t0 ON true

LEFT JOIN LATERAL (
  SELECT
    CASE
      WHEN MAX(ks.aki_stage_creat) >= 3
       AND MAX(ks.aki_stage_uo)    >= 3 THEN 'both'
      WHEN MAX(ks.aki_stage_uo)    >= 3 THEN 'oliguria'
      WHEN MAX(ks.aki_stage_creat) >= 3 THEN 'creatinine'
      ELSE                                   'unknown'
    END AS aki_pathway
  FROM mimiciv_derived.kdigo_stages ks
  WHERE ks.stay_id  = c.stay_id
    AND ks.charttime >  c.kdigo3_time - INTERVAL '12 hour'
    AND ks.charttime <= c.kdigo3_time + INTERVAL '6 hour'
) akip ON true

LEFT JOIN mimiciv_derived.creatinine_baseline creat_bl
  ON creat_bl.hadm_id = c.hadm_id
;

CREATE INDEX IF NOT EXISTS idx_baseline_v2_stay
  ON data_extract_crrt.cohort_baseline_v2 (stay_id);
CREATE INDEX IF NOT EXISTS idx_baseline_v2_hadm
  ON data_extract_crrt.cohort_baseline_v2 (hadm_id);
CREATE INDEX IF NOT EXISTS idx_baseline_v2_subject
  ON data_extract_crrt.cohort_baseline_v2 (subject_id);

ANALYZE data_extract_crrt.cohort_baseline_v2;

SELECT * FROM data_extract_crrt.cohort_baseline_v2;
