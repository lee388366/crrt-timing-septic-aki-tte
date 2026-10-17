# Diagnose and rescue the eICU external-validation cohort.
# This script writes only to diagnose_rescue/ and a new PostgreSQL schema.

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(tibble)
  library(survival)
})

set.seed(20260716)

root_dir <- "CONFIGURE_OTHER_SOURCE_ROOT/eICU数据提取"
out_dir <- file.path(root_dir, "diagnose_rescue")
data_dir <- file.path(out_dir, "data")
report_dir <- file.path(out_dir, "reports")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

pg_password <- Sys.getenv("PGPASSWORD")
if (!nzchar(pg_password)) {
  stop("PGPASSWORD is not set. Export the local PostgreSQL password before running.")
}

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = Sys.getenv("EICU_DBNAME", "eicu"),
  host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")),
  user = Sys.getenv("PGUSER", "postgres"),
  password = pg_password
)
on.exit(dbDisconnect(con), add = TRUE)

exec <- function(sql) invisible(dbExecute(con, sql))
qry <- function(sql) dbGetQuery(con, sql)

schema_name <- "data_extract_crrt_eicu_diagnose"
message("Creating diagnosis schema ", schema_name, " ...")
exec(sprintf("DROP SCHEMA IF EXISTS %s CASCADE", schema_name))
exec(sprintf("CREATE SCHEMA %s", schema_name))

message("Building clean physiologic KDIGO3 feature base ...")

exec("
CREATE TEMP TABLE patient_base AS
SELECT
  patientunitstayid,
  patienthealthsystemstayid,
  uniquepid,
  hospitalid,
  wardid,
  unitvisitnumber,
  unittype,
  unitstaytype,
  gender,
  CASE
    WHEN age = '> 89' THEN 90
    WHEN age ~ '^[0-9]+$' THEN age::integer
    ELSE NULL
  END AS age,
  ethnicity,
  admissionweight,
  unitdischargeoffset,
  hospitaldischargeoffset,
  hospitaldischargestatus,
  unitdischargestatus,
  apacheadmissiondx
FROM eicu_crd.patient
WHERE unitdischargeoffset IS NOT NULL
")
exec("CREATE INDEX ON patient_base(patientunitstayid)")
exec("CREATE INDEX ON patient_base(hospitalid)")

exec("
CREATE TEMP TABLE adult_patients AS
SELECT *
FROM patient_base
WHERE age >= 18
")
exec("CREATE INDEX ON adult_patients(patientunitstayid)")

exec("
CREATE TEMP TABLE sepsis_candidates AS
SELECT DISTINCT patientunitstayid
FROM (
  SELECT d.patientunitstayid
  FROM eicu_crd.diagnosis d
  JOIN adult_patients p USING (patientunitstayid)
  WHERE d.diagnosisstring ILIKE '%sepsis%'
     OR d.diagnosisstring ILIKE '%septic shock%'
  UNION
  SELECT a.patientunitstayid
  FROM eicu_crd.apachepredvar a
  JOIN adult_patients p USING (patientunitstayid)
  WHERE a.admitdiagnosis ILIKE '%sepsis%'
     OR a.admitdiagnosis ILIKE '%septic%'
) x
")
exec("CREATE INDEX ON sepsis_candidates(patientunitstayid)")

exec("
CREATE TEMP TABLE lab_clean AS
SELECT
  l.patientunitstayid,
  l.labresultoffset AS offset_min,
  lower(l.labname) AS labname,
  CASE
    WHEN lower(l.labname) = 'ph' AND l.labresult > 100 THEN l.labresult / 1000.0
    ELSE l.labresult
  END AS value
FROM eicu_crd.lab l
JOIN sepsis_candidates sc USING (patientunitstayid)
WHERE l.labresult IS NOT NULL
  AND l.labresultoffset IS NOT NULL
  AND lower(l.labname) IN ('creatinine','bun','potassium','bicarbonate','hco3','ph','lactate')
")
exec("CREATE INDEX ON lab_clean(patientunitstayid, offset_min)")
exec("CREATE INDEX ON lab_clean(patientunitstayid, labname, offset_min)")

exec("
CREATE TEMP TABLE creat_labs AS
SELECT patientunitstayid, offset_min, value AS creatinine
FROM lab_clean
WHERE labname = 'creatinine'
  AND value BETWEEN 0.1 AND 30
")
exec("CREATE INDEX ON creat_labs(patientunitstayid, offset_min)")

exec("
CREATE TEMP TABLE creat_baseline AS
SELECT
  patientunitstayid,
  COALESCE(
    MIN(creatinine) FILTER (WHERE offset_min BETWEEN -1440 AND 1440),
    MIN(creatinine)
  ) AS creatinine_baseline,
  CASE
    WHEN COUNT(*) FILTER (WHERE offset_min BETWEEN -1440 AND 1440) > 0
      THEN 'min_-24h_to_+24h_from_ICU_admit'
    ELSE 'fallback_min_observed_ICU'
  END AS creatinine_baseline_source
FROM creat_labs
GROUP BY patientunitstayid
HAVING COALESCE(
    MIN(creatinine) FILTER (WHERE offset_min BETWEEN -1440 AND 1440),
    MIN(creatinine)
  ) BETWEEN 0.3 AND 20
")
exec("CREATE INDEX ON creat_baseline(patientunitstayid)")

exec("
CREATE TEMP TABLE k3_creat AS
SELECT DISTINCT ON (c.patientunitstayid)
  c.patientunitstayid,
  c.offset_min AS kdigo3_creat_offset,
  c.creatinine AS creatinine_at_kdigo3,
  b.creatinine_baseline,
  b.creatinine_baseline_source
FROM creat_labs c
JOIN creat_baseline b USING (patientunitstayid)
WHERE c.offset_min >= 0
  AND (
    c.creatinine >= 3.0 * b.creatinine_baseline
    OR (c.creatinine >= 4.0 AND c.creatinine - b.creatinine_baseline >= 0.5)
  )
ORDER BY c.patientunitstayid, c.offset_min
")
exec("CREATE INDEX ON k3_creat(patientunitstayid)")

exec("
CREATE TEMP TABLE io_weight AS
SELECT
  io.patientunitstayid,
  MIN(CASE
    WHEN lower(io.celllabel) = 'bodyweight (kg)'
      AND io.cellvaluenumeric BETWEEN 30 AND 300 THEN io.cellvaluenumeric
    WHEN lower(io.celllabel) = 'bodyweight (lb)'
      AND io.cellvaluenumeric BETWEEN 66 AND 660 THEN io.cellvaluenumeric * 0.45359237
    ELSE NULL
  END) AS io_weight_kg
FROM eicu_crd.intakeoutput io
JOIN sepsis_candidates sc USING (patientunitstayid)
WHERE io.cellvaluenumeric IS NOT NULL
GROUP BY io.patientunitstayid
")
exec("CREATE INDEX ON io_weight(patientunitstayid)")

exec("
CREATE TEMP TABLE weight_base AS
SELECT
  p.patientunitstayid,
  COALESCE(
    CASE WHEN p.admissionweight BETWEEN 30 AND 300 THEN p.admissionweight END,
    iw.io_weight_kg
  ) AS weight_kg,
  CASE
    WHEN p.admissionweight BETWEEN 30 AND 300 THEN 'patient.admissionweight'
    WHEN iw.io_weight_kg IS NOT NULL THEN 'intakeoutput.bodyweight'
    ELSE 'missing'
  END AS weight_source
FROM adult_patients p
LEFT JOIN io_weight iw USING (patientunitstayid)
")
exec("CREATE INDEX ON weight_base(patientunitstayid)")

exec("
CREATE TEMP TABLE urine_events AS
SELECT
  io.patientunitstayid,
  io.intakeoutputoffset AS offset_min,
  SUM(io.outputtotal) AS urine_ml
FROM eicu_crd.intakeoutput io
JOIN sepsis_candidates sc USING (patientunitstayid)
WHERE io.intakeoutputoffset >= 0
  AND io.outputtotal IS NOT NULL
  AND io.outputtotal BETWEEN 0 AND 2500
  AND lower(io.cellpath) LIKE '%output (ml)%'
  AND (
    lower(io.cellpath) LIKE '%urine%'
    OR lower(io.celllabel) LIKE '%urine%'
    OR lower(io.celllabel) LIKE '%foley%'
    OR lower(io.celllabel) LIKE '%void%'
  )
  AND lower(io.celllabel) NOT LIKE '%count%'
  AND lower(io.celllabel) NOT LIKE '%occurrence%'
GROUP BY io.patientunitstayid, io.intakeoutputoffset
")
exec("CREATE INDEX ON urine_events(patientunitstayid, offset_min)")

exec("
CREATE TEMP TABLE urine_6h AS
SELECT
  patientunitstayid,
  floor(offset_min / 360.0)::integer AS bin6,
  SUM(urine_ml) AS urine_ml_6h,
  COUNT(*) AS urine_entries_6h
FROM urine_events
GROUP BY patientunitstayid, floor(offset_min / 360.0)::integer
")
exec("CREATE INDEX ON urine_6h(patientunitstayid, bin6)")

exec("
CREATE TEMP TABLE uo_roll AS
SELECT
  u.patientunitstayid,
  (u.bin6 + 1) * 360 AS offset_min,
  SUM(u.urine_ml_6h) OVER (
    PARTITION BY u.patientunitstayid
    ORDER BY u.bin6
    ROWS BETWEEN 3 PRECEDING AND CURRENT ROW
  ) AS urine_24h_ml,
  COUNT(*) OVER (
    PARTITION BY u.patientunitstayid
    ORDER BY u.bin6
    ROWS BETWEEN 3 PRECEDING AND CURRENT ROW
  ) AS urine_bins_24h,
  SUM(u.urine_ml_6h) OVER (
    PARTITION BY u.patientunitstayid
    ORDER BY u.bin6
    ROWS BETWEEN 1 PRECEDING AND CURRENT ROW
  ) AS urine_12h_ml,
  COUNT(*) OVER (
    PARTITION BY u.patientunitstayid
    ORDER BY u.bin6
    ROWS BETWEEN 1 PRECEDING AND CURRENT ROW
  ) AS urine_bins_12h
FROM urine_6h u
")
exec("CREATE INDEX ON uo_roll(patientunitstayid, offset_min)")

exec("
CREATE TEMP TABLE k3_uo AS
SELECT DISTINCT ON (u.patientunitstayid)
  u.patientunitstayid,
  u.offset_min AS kdigo3_uo_offset,
  u.urine_24h_ml,
  u.urine_bins_24h,
  u.urine_12h_ml,
  u.urine_bins_12h,
  w.weight_kg
FROM uo_roll u
JOIN weight_base w USING (patientunitstayid)
WHERE w.weight_kg BETWEEN 30 AND 300
  AND u.offset_min >= 0
  AND (
    (u.urine_bins_24h >= 4 AND u.urine_24h_ml / w.weight_kg / 24.0 < 0.3)
    OR (u.urine_bins_12h >= 2 AND u.urine_12h_ml = 0)
  )
ORDER BY u.patientunitstayid, u.offset_min
")
exec("CREATE INDEX ON k3_uo(patientunitstayid)")

exec("
CREATE TEMP TABLE clean_kdigo3 AS
SELECT
  COALESCE(c.patientunitstayid, u.patientunitstayid) AS patientunitstayid,
  CASE
    WHEN c.kdigo3_creat_offset IS NOT NULL
      AND (u.kdigo3_uo_offset IS NULL OR c.kdigo3_creat_offset <= u.kdigo3_uo_offset)
      THEN c.kdigo3_creat_offset
    ELSE u.kdigo3_uo_offset
  END AS t0_offset,
  c.kdigo3_creat_offset,
  u.kdigo3_uo_offset,
  (c.kdigo3_creat_offset IS NOT NULL) AS stage3_by_creat,
  (u.kdigo3_uo_offset IS NOT NULL) AS stage3_by_uo,
  c.creatinine_at_kdigo3,
  c.creatinine_baseline,
  c.creatinine_baseline_source,
  u.urine_24h_ml,
  u.urine_bins_24h,
  u.urine_12h_ml,
  u.urine_bins_12h,
  u.weight_kg AS uo_weight_kg
FROM k3_creat c
FULL OUTER JOIN k3_uo u USING (patientunitstayid)
")
exec("CREATE INDEX ON clean_kdigo3(patientunitstayid)")
exec("CREATE INDEX ON clean_kdigo3(t0_offset)")

exec("
CREATE TEMP TABLE crrt_first AS
SELECT
  patientunitstayid,
  MIN(treatmentoffset) AS crrt_first_offset
FROM eicu_crd.treatment
WHERE treatmentoffset IS NOT NULL
  AND (
    treatmentstring ILIKE '%C V V H%'
    OR treatmentstring ILIKE '%CVVH%'
    OR treatmentstring ILIKE '%CVVHD%'
    OR treatmentstring ILIKE '%CRRT%'
    OR treatmentstring ILIKE '%continuous%renal%'
  )
GROUP BY patientunitstayid
")
exec("CREATE INDEX ON crrt_first(patientunitstayid)")

exec("
CREATE TEMP TABLE chronic_dialysis AS
SELECT DISTINCT patientunitstayid
FROM (
  SELECT patientunitstayid
  FROM eicu_crd.diagnosis
  WHERE diagnosisstring ILIKE '%end stage renal%'
     OR diagnosisstring ILIKE '%ESRD%'
     OR diagnosisstring ILIKE '%chronic renal failure%'
  UNION
  SELECT patientunitstayid
  FROM eicu_crd.treatment
  WHERE treatmentstring ILIKE '%for chronic renal failure%'
     OR treatmentstring ILIKE '%peritoneal dialysis%for chronic%'
) x
")
exec("CREATE INDEX ON chronic_dialysis(patientunitstayid)")

exec("
CREATE TEMP TABLE sepsis_diag_features AS
SELECT
  k.patientunitstayid,
  bool_or(d.diagnosisstring ILIKE '%sepsis%' OR d.diagnosisstring ILIKE '%septic shock%')
    FILTER (WHERE d.diagnosisoffset IS NULL OR d.diagnosisoffset <= k.t0_offset + 360) AS diagnosis_sepsis_before_t0,
  bool_or(d.diagnosisstring ILIKE '%septic shock%')
    FILTER (WHERE d.diagnosisoffset IS NULL OR d.diagnosisoffset <= k.t0_offset + 360) AS diagnosis_septic_shock_before_t0
FROM clean_kdigo3 k
LEFT JOIN eicu_crd.diagnosis d USING (patientunitstayid)
GROUP BY k.patientunitstayid
")
exec("CREATE INDEX ON sepsis_diag_features(patientunitstayid)")

exec("
CREATE TEMP TABLE apache_features AS
SELECT
  k.patientunitstayid,
  bool_or(a.admitdiagnosis ILIKE '%sepsis%' OR a.admitdiagnosis ILIKE '%septic%') AS apache_sepsis,
  bool_or(a.admitdiagnosis ILIKE '%septic shock%') AS apache_septic_shock
FROM clean_kdigo3 k
LEFT JOIN eicu_crd.apachepredvar a USING (patientunitstayid)
GROUP BY k.patientunitstayid
")
exec("CREATE INDEX ON apache_features(patientunitstayid)")

exec("
CREATE TEMP TABLE vaso_events AS
SELECT
  i.patientunitstayid,
  i.infusionoffset AS offset_min,
  lower(i.drugname) AS drugname,
  NULLIF(substring(coalesce(i.drugrate, i.infusionrate, '') from '[-+]?[0-9]*\\.?[0-9]+'), '')::numeric AS rate_numeric
FROM eicu_crd.infusiondrug i
JOIN clean_kdigo3 k USING (patientunitstayid)
WHERE i.infusionoffset IS NOT NULL
  AND (
    i.drugname ILIKE '%norepi%'
    OR i.drugname ILIKE '%levophed%'
    OR i.drugname ILIKE '%epi%'
    OR i.drugname ILIKE '%dopamine%'
    OR i.drugname ILIKE '%phenyl%'
    OR i.drugname ILIKE '%neo%'
    OR i.drugname ILIKE '%vasopressin%'
  )
  AND i.drugname NOT ILIKE '%epidural%'
")
exec("CREATE INDEX ON vaso_events(patientunitstayid, offset_min)")

exec("
CREATE TEMP TABLE shock_features AS
SELECT
  k.patientunitstayid,
  COUNT(v.*) FILTER (WHERE v.offset_min BETWEEN k.t0_offset - 360 AND k.t0_offset + 60)::integer AS vaso_rows_near_t0,
  COUNT(v.*) FILTER (
    WHERE v.offset_min BETWEEN k.t0_offset - 360 AND k.t0_offset + 60
      AND v.rate_numeric IS NOT NULL
  )::integer AS vaso_numeric_rows_near_t0,
  COUNT(l.*) FILTER (
    WHERE l.labname = 'lactate'
      AND l.value BETWEEN 0.1 AND 30
      AND l.offset_min BETWEEN k.t0_offset - 720 AND k.t0_offset + 360
  )::integer AS lactate_n_near_t0,
  MAX(l.value) FILTER (
    WHERE l.labname = 'lactate'
      AND l.value BETWEEN 0.1 AND 30
      AND l.offset_min BETWEEN k.t0_offset - 720 AND k.t0_offset + 360
  ) AS lactate_near_t0
FROM clean_kdigo3 k
LEFT JOIN vaso_events v USING (patientunitstayid)
LEFT JOIN lab_clean l USING (patientunitstayid)
GROUP BY k.patientunitstayid
")
exec("CREATE INDEX ON shock_features(patientunitstayid)")

exec("
CREATE TEMP TABLE hospital_vaso_quality AS
SELECT
  p.hospitalid,
  COUNT(v.*) AS vaso_rows,
  COUNT(v.*) FILTER (WHERE v.rate_numeric IS NOT NULL) AS vaso_rate_numeric_rows,
  COALESCE(COUNT(v.*) FILTER (WHERE v.rate_numeric IS NOT NULL)::numeric / NULLIF(COUNT(v.*), 0), 0) AS vaso_rate_numeric_prop
FROM adult_patients p
LEFT JOIN vaso_events v USING (patientunitstayid)
GROUP BY p.hospitalid
")

exec("
CREATE TEMP TABLE feature_pre_rescue AS
SELECT
  p.patientunitstayid,
  p.patienthealthsystemstayid,
  p.uniquepid,
  p.hospitalid,
  p.wardid,
  p.unitvisitnumber,
  p.unittype,
  p.unitstaytype,
  p.gender,
  p.age,
  p.ethnicity,
  p.admissionweight,
  wb.weight_kg,
  wb.weight_source,
  p.unitdischargeoffset,
  p.hospitaldischargeoffset,
  p.hospitaldischargestatus,
  p.unitdischargestatus,
  p.apacheadmissiondx,
  k.t0_offset,
  k.kdigo3_creat_offset,
  k.kdigo3_uo_offset,
  k.stage3_by_creat,
  k.stage3_by_uo,
  k.creatinine_at_kdigo3,
  k.creatinine_baseline,
  k.creatinine_baseline_source,
  k.urine_24h_ml,
  k.urine_bins_24h,
  k.urine_12h_ml,
  k.urine_bins_12h,
  c.crrt_first_offset,
  (cd.patientunitstayid IS NOT NULL) AS chronic_dialysis_or_esrd,
  COALESCE(sd.diagnosis_sepsis_before_t0, false) AS diagnosis_sepsis_before_t0,
  COALESCE(af.apache_sepsis, false) AS apache_sepsis,
  COALESCE(sd.diagnosis_septic_shock_before_t0, false) AS diagnosis_septic_shock_before_t0,
  COALESCE(af.apache_septic_shock, false) AS apache_septic_shock,
  COALESCE(sf.vaso_rows_near_t0, 0) AS vaso_rows_near_t0,
  COALESCE(sf.vaso_numeric_rows_near_t0, 0) AS vaso_numeric_rows_near_t0,
  COALESCE(sf.lactate_n_near_t0, 0) AS lactate_n_near_t0,
  sf.lactate_near_t0,
  apr.acutephysiologyscore,
  apr.apachescore,
  apv.admitdiagnosis AS apache_admitdiagnosis,
  aps.meanbp AS apache_meanbp,
  aps.ph AS apache_ph,
  aps.bun AS apache_bun,
  aps.creatinine AS apache_creatinine,
  aps.urine AS apache_urine,
  aps.eyes AS apache_eyes,
  aps.motor AS apache_motor,
  aps.verbal AS apache_verbal
FROM adult_patients p
JOIN clean_kdigo3 k USING (patientunitstayid)
LEFT JOIN weight_base wb USING (patientunitstayid)
LEFT JOIN crrt_first c USING (patientunitstayid)
LEFT JOIN chronic_dialysis cd USING (patientunitstayid)
LEFT JOIN sepsis_diag_features sd USING (patientunitstayid)
LEFT JOIN apache_features af USING (patientunitstayid)
LEFT JOIN shock_features sf USING (patientunitstayid)
LEFT JOIN eicu_crd.apachepatientresult apr USING (patientunitstayid)
LEFT JOIN eicu_crd.apachepredvar apv USING (patientunitstayid)
LEFT JOIN eicu_crd.apacheapsvar aps USING (patientunitstayid)
")
exec("CREATE INDEX ON feature_pre_rescue(patientunitstayid)")
exec("CREATE INDEX ON feature_pre_rescue(t0_offset)")

exec("
CREATE TEMP TABLE feature_pre_rescue_dedup AS
SELECT DISTINCT ON (patientunitstayid) *
FROM feature_pre_rescue
ORDER BY patientunitstayid
")
exec("CREATE INDEX ON feature_pre_rescue_dedup(patientunitstayid)")
exec("CREATE INDEX ON feature_pre_rescue_dedup(t0_offset)")

exec("
CREATE TEMP TABLE feature_rescue AS
SELECT
  f.*,
  COALESCE(kstats.k_n, 0) AS r1_k_n,
  COALESCE(kstats.k_high_n, 0) AS r1_k_high_n,
  COALESCE(kstats.r1_trigger, false) AS r1_trigger,
  COALESCE(phstats.ph_n, 0) AS r2_ph_n,
  COALESCE(phstats.r2_trigger, false) AS r2_trigger,
  COALESCE(hco3.hco3_n, 0) AS r2b_hco3_n,
  COALESCE(hco3.r2b_trigger, false) AS r2b_trigger,
  COALESCE(bun.bun_n, 0) AS r3_bun_n,
  COALESCE(bun.r3_trigger, false) AS r3_trigger
FROM feature_pre_rescue_dedup f
LEFT JOIN LATERAL (
  SELECT
    COUNT(*) FILTER (WHERE l.value IS NOT NULL) AS k_n,
    COUNT(*) FILTER (WHERE l.value >= 6.5) AS k_high_n,
    bool_or(l.value >= 6.5) AS r1_trigger
  FROM lab_clean l
  WHERE l.patientunitstayid = f.patientunitstayid
    AND l.labname = 'potassium'
    AND l.value BETWEEN 1 AND 12
    AND l.offset_min > f.t0_offset - 360
    AND l.offset_min <= f.t0_offset
) kstats ON true
LEFT JOIN LATERAL (
  SELECT
    COUNT(*) FILTER (WHERE l.value IS NOT NULL) AS ph_n,
    bool_or(l.value <= 7.15) AS r2_trigger
  FROM lab_clean l
  WHERE l.patientunitstayid = f.patientunitstayid
    AND l.labname = 'ph'
    AND l.value BETWEEN 6.5 AND 7.8
    AND l.offset_min > f.t0_offset - 360
    AND l.offset_min <= f.t0_offset
) phstats ON true
LEFT JOIN LATERAL (
  SELECT
    COUNT(*) FILTER (WHERE l.value IS NOT NULL) AS hco3_n,
    bool_or(l.value <= 12) AS r2b_trigger
  FROM lab_clean l
  WHERE l.patientunitstayid = f.patientunitstayid
    AND l.labname IN ('bicarbonate','hco3')
    AND l.value BETWEEN 1 AND 80
    AND l.offset_min > f.t0_offset - 360
    AND l.offset_min <= f.t0_offset
) hco3 ON true
LEFT JOIN LATERAL (
  SELECT
    COUNT(*) FILTER (WHERE l.value IS NOT NULL) AS bun_n,
    bool_or(l.value >= 100) AS r3_trigger
  FROM lab_clean l
  WHERE l.patientunitstayid = f.patientunitstayid
    AND l.labname = 'bun'
    AND l.value BETWEEN 1 AND 300
    AND l.offset_min > f.t0_offset - 1440
    AND l.offset_min <= f.t0_offset
) bun ON true
")

exec("
CREATE TEMP TABLE uo_doc_24h AS
SELECT
  f.patientunitstayid,
  COUNT(ue.*) AS urine_entries_0_24h
FROM feature_rescue f
LEFT JOIN urine_events ue
  ON ue.patientunitstayid = f.patientunitstayid
 AND ue.offset_min >= f.t0_offset
 AND ue.offset_min < f.t0_offset + 1440
GROUP BY f.patientunitstayid
")
exec("CREATE INDEX ON uo_doc_24h(patientunitstayid)")

exec(sprintf("
CREATE TABLE %s.kdigo3_feature_base AS
SELECT
  f.*,
  COALESCE(u.urine_entries_0_24h, 0) AS urine_entries_0_24h,
  COALESCE(hv.vaso_rate_numeric_prop, 0) AS hospital_vaso_rate_numeric_prop,
  (NOT f.chronic_dialysis_or_esrd
    AND (f.crrt_first_offset IS NULL OR f.crrt_first_offset >= f.t0_offset)
    AND f.unitdischargeoffset > f.t0_offset
    AND COALESCE(f.hospitaldischargeoffset, f.unitdischargeoffset) > f.t0_offset
  ) AS core_eligible,
  (COALESCE(f.diagnosis_sepsis_before_t0, false) OR COALESCE(f.apache_sepsis, false)) AS sepsis_evidence_at_t0,
  (COALESCE(f.diagnosis_septic_shock_before_t0, false) OR COALESCE(f.apache_septic_shock, false)) AS septic_shock_diagnosis_proxy,
  (COALESCE(f.vaso_rows_near_t0, 0) > 0) AS vaso_any_near_t0,
  (COALESCE(f.vaso_numeric_rows_near_t0, 0) > 0) AS vaso_numeric_near_t0,
  (COALESCE(f.vaso_rows_near_t0, 0) > 0 AND COALESCE(f.vaso_numeric_rows_near_t0, 0) = 0) AS vaso_record_no_numeric_rate,
  (COALESCE(f.vaso_rows_near_t0, 0) = 0) AS no_vaso_record_near_t0,
  (COALESCE(f.lactate_n_near_t0, 0) > 0) AS lactate_measured_near_t0,
  (COALESCE(f.lactate_near_t0, -999) > 2) AS lactate_gt2_near_t0,
  (COALESCE(f.lactate_n_near_t0, 0) > 0 AND COALESCE(f.lactate_near_t0, -999) <= 2) AS lactate_le2_near_t0,
  (COALESCE(f.lactate_n_near_t0, 0) = 0) AS lactate_missing_near_t0,
  (r1_trigger OR r2_trigger OR r2b_trigger OR r3_trigger) AS rescue_primary_r1_r2_r2b_r3,
  (r1_trigger OR r2b_trigger OR r3_trigger) AS rescue_no_r2,
  COALESCE(f.crrt_first_offset BETWEEN f.t0_offset AND f.t0_offset + 1440, false) AS early_crrt,
  CASE
    WHEN f.stage3_by_creat AND f.stage3_by_uo THEN 'both_creatinine_and_uo'
    WHEN f.stage3_by_creat THEN 'creatinine_only'
    WHEN f.stage3_by_uo THEN 'urine_output_only'
    ELSE 'no_creat_or_uo_stage3_ERROR'
  END AS stage3_path,
  (f.hospitaldischargestatus ILIKE '%%Expired%%'
    AND f.hospitaldischargeoffset IS NOT NULL
    AND f.hospitaldischargeoffset > f.t0_offset
    AND f.hospitaldischargeoffset <= f.t0_offset + 40320) AS died_28d,
  CASE
    WHEN f.hospitaldischargestatus ILIKE '%%Expired%%'
      AND f.hospitaldischargeoffset IS NOT NULL
      THEN (f.hospitaldischargeoffset - f.t0_offset) / 1440.0
    WHEN f.hospitaldischargeoffset IS NOT NULL
      THEN LEAST(28.0, GREATEST(0.0, (f.hospitaldischargeoffset - f.t0_offset) / 1440.0))
    ELSE LEAST(28.0, GREATEST(0.0, (f.unitdischargeoffset - f.t0_offset) / 1440.0))
  END AS follow_time_28d,
  CASE
    WHEN f.hospitaldischargestatus ILIKE '%%Expired%%'
      AND f.hospitaldischargeoffset IS NOT NULL
      THEN (f.hospitaldischargeoffset - f.t0_offset) / 1440.0
    ELSE NULL
  END AS days_to_death,
  (f.hospitaldischargestatus NOT ILIKE '%%Expired%%'
    AND f.hospitaldischargeoffset IS NOT NULL
    AND f.hospitaldischargeoffset > f.t0_offset
    AND f.hospitaldischargeoffset <= f.t0_offset + 40320) AS discharged_alive_before_28d,
  (f.crrt_first_offset - f.t0_offset) / 60.0 AS crrt_hours_from_t0
FROM feature_rescue f
LEFT JOIN uo_doc_24h u USING (patientunitstayid)
LEFT JOIN hospital_vaso_quality hv USING (hospitalid)
", schema_name))
exec(sprintf("CREATE INDEX ON %s.kdigo3_feature_base(patientunitstayid)", schema_name))
exec(sprintf("CREATE INDEX ON %s.kdigo3_feature_base(hospitalid)", schema_name))

feature <- qry(sprintf("SELECT * FROM %s.kdigo3_feature_base ORDER BY patientunitstayid", schema_name))
write_csv(feature, file.path(data_dir, "task0_kdigo3_feature_base.csv"))

message("Feature base n: ", nrow(feature))

truthy <- function(x) x %in% TRUE

coalesce_unitvisit <- function(x) ifelse(is.na(x), 999999L, x)

build_cohort <- function(df, def_flag) {
  src <- df %>%
    filter(.data[[def_flag]] %in% TRUE) %>%
    arrange(uniquepid, coalesce_unitvisit(unitvisitnumber), patientunitstayid)
  pre <- src %>%
    group_by(uniquepid) %>%
    slice(1) %>%
    ungroup()
  base <- pre %>% filter(!rescue_primary_r1_r2_r2b_r3)
  no_r2 <- pre %>% filter(!rescue_no_r2)
  list(source_all_stays = src, pre_rescue = pre, cohort_base = base, cohort_no_r2 = no_r2)
}

add_definition_flags <- function(df) {
  df %>%
    mutate(
      def_old_any_vaso_lactate = core_eligible & sepsis_evidence_at_t0 & vaso_any_near_t0 & lactate_gt2_near_t0,
      def_a_strict_numeric_lactate = core_eligible & sepsis_evidence_at_t0 & vaso_numeric_near_t0 & lactate_gt2_near_t0,
      def_b_vaso_or_lactate = core_eligible & sepsis_evidence_at_t0 & (vaso_any_near_t0 | lactate_gt2_near_t0),
      def_c_any_vaso_no_lactate_requirement = core_eligible & sepsis_evidence_at_t0 & vaso_any_near_t0,
      def_d_septic_shock_diagnosis_proxy = core_eligible & septic_shock_diagnosis_proxy
    )
}

feature <- add_definition_flags(feature)

definition_specs <- tibble(
  definition = c(
    "Old current script: sepsis + any vasopressor + lactate above 2",
    "Def A: sepsis + numeric-rate vasopressor + lactate above 2",
    "Def B: sepsis + any vasopressor or lactate above 2",
    "Def C: sepsis + any vasopressor; lactate not required",
    "Def D: diagnosis or APACHE septic shock proxy"
  ),
  flag = c(
    "def_old_any_vaso_lactate",
    "def_a_strict_numeric_lactate",
    "def_b_vaso_or_lactate",
    "def_c_any_vaso_no_lactate_requirement",
    "def_d_septic_shock_diagnosis_proxy"
  ),
  defensibility = c(
    "Reproduces the previous executable pipeline, but lactate documentation can make it too specific in eICU.",
    "Most documentation-strict shock definition; useful as sensitivity, but likely excludes true vasopressor-treated shock when rates are missing.",
    "Maximizes sensitivity but may include hyperlactatemia without vasopressor-supported shock; better as broad sensitivity.",
    "Clinically defensible for septic shock treatment need in eICU because vasopressor record is less dependent on numeric-rate completeness and lactate is not forced.",
    "Diagnosis/APACHE label is transparent but vulnerable to coding heterogeneity; useful as a label-based sensitivity."
  )
)

definitions <- setNames(lapply(definition_specs$flag, function(flag) build_cohort(feature, flag)), definition_specs$definition)

task1_crosstab <- feature %>%
  mutate(
    vaso_category = case_when(
      vaso_numeric_near_t0 ~ "vasopressor record with numeric rate",
      vaso_any_near_t0 ~ "vasopressor record without numeric rate",
      TRUE ~ "no vasopressor record"
    ),
    lactate_category = case_when(
      lactate_missing_near_t0 ~ "no lactate measurement",
      lactate_gt2_near_t0 ~ "lactate measured above 2",
      TRUE ~ "lactate measured 2 or below"
    )
  ) %>%
  count(vaso_category, lactate_category, name = "n") %>%
  mutate(pct_of_5932 = round(100 * n / sum(n), 1))

task1_loss_priority <- feature %>%
  mutate(old_current_proxy = def_old_any_vaso_lactate) %>%
  filter(!old_current_proxy) %>%
  mutate(
    loss_reason_priority = case_when(
      !core_eligible ~ "not core eligible: chronic dialysis, prior CRRT, or no follow-up after t0",
      !sepsis_evidence_at_t0 ~ "no sepsis evidence at or before t0 plus 6h",
      no_vaso_record_near_t0 & lactate_measured_near_t0 & lactate_le2_near_t0 ~ "measured non-shock: no vasopressor and lactate 2 or below",
      vaso_record_no_numeric_rate ~ "vasopressor record present but no numeric rate",
      no_vaso_record_near_t0 ~ "no vasopressor record near t0",
      lactate_missing_near_t0 ~ "lactate missing near t0",
      TRUE ~ "other"
    )
  ) %>%
  count(loss_reason_priority, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(pct_of_excluded = round(100 * n / sum(n), 1))

task1_overlap <- feature %>%
  filter(!def_old_any_vaso_lactate) %>%
  summarise(
    excluded_from_old_current_n = n(),
    not_core_eligible_n = sum(!core_eligible, na.rm = TRUE),
    no_sepsis_evidence_n = sum(core_eligible & !sepsis_evidence_at_t0, na.rm = TRUE),
    no_vaso_record_n = sum(core_eligible & sepsis_evidence_at_t0 & no_vaso_record_near_t0, na.rm = TRUE),
    vaso_record_no_numeric_rate_n = sum(core_eligible & sepsis_evidence_at_t0 & vaso_record_no_numeric_rate, na.rm = TRUE),
    lactate_missing_n = sum(core_eligible & sepsis_evidence_at_t0 & lactate_missing_near_t0, na.rm = TRUE),
    measured_nonshock_no_vaso_lactate_le2_n = sum(core_eligible & sepsis_evidence_at_t0 & no_vaso_record_near_t0 & lactate_measured_near_t0 & lactate_le2_near_t0, na.rm = TRUE)
  )

sequential_counts <- c(
  "first physiologic KDIGO3 among adult sepsis candidates" = nrow(feature),
  "after chronic dialysis or ESRD exclusion" = sum(!feature$chronic_dialysis_or_esrd, na.rm = TRUE),
  "after prior CRRT exclusion" = sum(!feature$chronic_dialysis_or_esrd &
                                      (is.na(feature$crrt_first_offset) | feature$crrt_first_offset >= feature$t0_offset), na.rm = TRUE),
  "after requiring positive follow-up after t0" = sum(feature$core_eligible, na.rm = TRUE),
  "after requiring sepsis evidence at or before t0 plus 6h" = sum(feature$core_eligible & feature$sepsis_evidence_at_t0, na.rm = TRUE),
  "after old shock ascertainment: any vasopressor plus lactate above 2" = sum(feature$def_old_any_vaso_lactate, na.rm = TRUE),
  "after first stay per patient under old definition" = nrow(definitions[["Old current script: sepsis + any vasopressor + lactate above 2"]]$pre_rescue),
  "after primary rescue exclusion under old definition" = nrow(definitions[["Old current script: sepsis + any vasopressor + lactate above 2"]]$cohort_base)
)

task1_sequential_funnel <- tibble(
  step = names(sequential_counts),
  n = as.integer(sequential_counts),
  previous_n = dplyr::lag(n),
  excluded_from_previous = previous_n - n,
  retained_pct_from_previous = round(100 * n / previous_n, 1)
)

old_funnel_path <- file.path(root_dir, "data", "eicu_extraction_funnel.csv")
old_reported_shock_rows <- NA_integer_
if (file.exists(old_funnel_path)) {
  old_funnel <- read_csv(old_funnel_path, show_col_types = FALSE)
  old_reported_shock_rows <- old_funnel %>%
    filter(step == "septic_shock_proxy_at_t0") %>%
    pull(n) %>%
    as.integer()
}
task0_old_row_reconciliation <- tibble(
  item = c(
    "old reported septic_shock_proxy_at_t0 rows",
    "deduplicated unique stays meeting old executable shock definition",
    "first stay per patient under old executable shock definition"
  ),
  n = c(
    old_reported_shock_rows,
    sum(feature$def_old_any_vaso_lactate, na.rm = TRUE),
    nrow(definitions[["Old current script: sepsis + any vasopressor + lactate above 2"]]$pre_rescue)
  ),
  interpretation = c(
    "Old scaffold row count before correcting join duplication.",
    "Correct unique patientunitstayid count after deduplication.",
    "Patient-level first-stay count used for downstream cohort construction."
  )
)

compute_hospital_quality <- function(pre_df, min_n, uo_prop, vaso_prop) {
  if (nrow(pre_df) == 0) {
    return(tibble(
      hospitalid = integer(),
      pre_rescue_n = integer(),
      urine_doc_24h_prop = numeric(),
      vaso_rate_numeric_prop = numeric(),
      high_fidelity_hospital = logical()
    ))
  }
  pre_df %>%
    group_by(hospitalid) %>%
    summarise(
      pre_rescue_n = n(),
      urine_doc_24h_prop = mean(urine_entries_0_24h >= 4, na.rm = TRUE),
      vaso_rate_numeric_prop = max(hospital_vaso_rate_numeric_prop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      high_fidelity_hospital = pre_rescue_n >= min_n &
        urine_doc_24h_prop >= uo_prop &
        vaso_rate_numeric_prop >= vaso_prop
    )
}

summarise_definition <- function(def_name, obj) {
  hq <- compute_hospital_quality(obj$pre_rescue, 20, 0.70, 0.70)
  base <- obj$cohort_base %>% left_join(select(hq, hospitalid, high_fidelity_hospital), by = "hospitalid")
  stage <- base %>%
    count(stage3_path, early_crrt, name = "n") %>%
    mutate(label = paste0(stage3_path, " / ", ifelse(early_crrt, "early", "deferred"), ": ", n)) %>%
    pull(label) %>%
    paste(collapse = "; ")
  tibble(
    definition = def_name,
    source_all_stays_n = nrow(obj$source_all_stays),
    pre_rescue_first_stay_n = nrow(obj$pre_rescue),
    cohort_base_n = nrow(base),
    early_n = sum(base$early_crrt, na.rm = TRUE),
    early_pct = round(100 * mean(base$early_crrt, na.rm = TRUE), 1),
    death_28d_n = sum(base$died_28d, na.rm = TRUE),
    discharged_alive_before_28d_n = sum(base$discharged_alive_before_28d, na.rm = TRUE),
    landmark_n = sum(base$follow_time_28d > 1, na.rm = TRUE),
    landmark_early_n = sum(base$follow_time_28d > 1 & base$early_crrt, na.rm = TRUE),
    high_fidelity_hospital_n_default = sum(hq$high_fidelity_hospital, na.rm = TRUE),
    high_fidelity_cohort_n_default = sum(base$high_fidelity_hospital, na.rm = TRUE),
    high_fidelity_early_n_default = sum(base$high_fidelity_hospital & base$early_crrt, na.rm = TRUE),
    stage3_path_composition = stage
  )
}

task2_definition_comparison <- bind_rows(Map(summarise_definition, names(definitions), definitions)) %>%
  left_join(definition_specs, by = "definition") %>%
  select(definition, defensibility, source_all_stays_n, pre_rescue_first_stay_n, cohort_base_n, early_n, early_pct,
         death_28d_n, discharged_alive_before_28d_n, landmark_n, landmark_early_n,
         high_fidelity_hospital_n_default, high_fidelity_cohort_n_default, high_fidelity_early_n_default,
         stage3_path_composition)

threshold_grid <- expand.grid(
  min_hospital_n = c(10L, 20L),
  min_uo_doc_prop = c(0.50, 0.60, 0.70),
  min_vaso_rate_prop = c(0.50, 0.60, 0.70),
  KEEP.OUT.ATTRS = FALSE
) %>% as_tibble()

summarise_grid_one <- function(def_name, obj, min_n, uo_prop, vaso_prop) {
  hq <- compute_hospital_quality(obj$pre_rescue, min_n, uo_prop, vaso_prop)
  base <- obj$cohort_base %>% left_join(select(hq, hospitalid, high_fidelity_hospital), by = "hospitalid")
  tibble(
    definition = def_name,
    min_hospital_n = min_n,
    min_uo_doc_prop = uo_prop,
    min_vaso_rate_prop = vaso_prop,
    high_fidelity_hospital_n = sum(hq$high_fidelity_hospital, na.rm = TRUE),
    high_fidelity_cohort_n = sum(base$high_fidelity_hospital, na.rm = TRUE),
    high_fidelity_early_n = sum(base$high_fidelity_hospital & base$early_crrt, na.rm = TRUE),
    high_fidelity_event_n = sum(base$high_fidelity_hospital & base$died_28d, na.rm = TRUE),
    high_fidelity_landmark_n = sum(base$high_fidelity_hospital & base$follow_time_28d > 1, na.rm = TRUE),
    high_fidelity_landmark_early_n = sum(base$high_fidelity_hospital & base$follow_time_28d > 1 & base$early_crrt, na.rm = TRUE),
    all_hospital_cohort_n = nrow(base),
    all_hospital_early_n = sum(base$early_crrt, na.rm = TRUE),
    all_hospital_event_n = sum(base$died_28d, na.rm = TRUE)
  )
}

task3_high_fidelity_grid <- bind_rows(lapply(names(definitions), function(def_name) {
  obj <- definitions[[def_name]]
  bind_rows(lapply(seq_len(nrow(threshold_grid)), function(i) {
    summarise_grid_one(
      def_name, obj,
      threshold_grid$min_hospital_n[i],
      threshold_grid$min_uo_doc_prop[i],
      threshold_grid$min_vaso_rate_prop[i]
    )
  }))
}))

recommended_definition <- "Def C: sepsis + any vasopressor; lactate not required"
def_c_grid <- task3_high_fidelity_grid %>%
  filter(definition == recommended_definition)

eligible_grid <- def_c_grid %>%
  filter(high_fidelity_early_n >= 30) %>%
  arrange(desc(min_uo_doc_prop + min_vaso_rate_prop), desc(min_hospital_n), desc(high_fidelity_early_n))

if (nrow(eligible_grid) > 0) {
  recommended_grid <- eligible_grid %>% slice(1)
} else {
  recommended_grid <- def_c_grid %>%
    arrange(desc(high_fidelity_early_n), desc(high_fidelity_cohort_n), desc(min_uo_doc_prop + min_vaso_rate_prop), desc(min_hospital_n)) %>%
    slice(1)
}

recommended_obj <- definitions[[recommended_definition]]
recommended_hq <- compute_hospital_quality(
  recommended_obj$pre_rescue,
  recommended_grid$min_hospital_n,
  recommended_grid$min_uo_doc_prop,
  recommended_grid$min_vaso_rate_prop
)
recommended_cohort <- recommended_obj$cohort_base %>%
  left_join(select(recommended_hq, hospitalid, high_fidelity_hospital), by = "hospitalid")

write_csv(recommended_cohort, file.path(data_dir, "task4_recommended_def_c_cohort_base.csv"))

numeric_covariates <- c(
  "age",
  "apachescore",
  "acutephysiologyscore",
  "creatinine_baseline",
  "creatinine_at_kdigo3",
  "lactate_near_t0",
  "weight_kg",
  "apache_meanbp",
  "apache_ph",
  "apache_bun",
  "apache_creatinine",
  "apache_urine"
)
available_numeric_covariates <- intersect(numeric_covariates, names(recommended_cohort))
MIMIC_HR_REFERENCE <- as.numeric(Sys.getenv("MIMIC_STRICTK3_V2_HR", "1.044"))
TRIM_Q <- as.numeric(Sys.getenv("EICU_WEIGHT_TRIM_Q", "0.99"))

impute_median <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    df[[v]][df[[v]] %in% c(-1, -99, -999)] <- NA
    med <- median(df[[v]], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    df[[v]][is.na(df[[v]])] <- med
  }
  df
}

standardize_numeric <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    s <- sd(df[[v]], na.rm = TRUE)
    m <- mean(df[[v]], na.rm = TRUE)
    if (is.finite(s) && s > 0) {
      df[[v]] <- (df[[v]] - m) / s
    }
  }
  df
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

smd_numeric <- function(df, var, weight_var = NULL, treat_var = "early_crrt") {
  a <- df[df[[treat_var]] %in% TRUE & is.finite(df[[var]]), ]
  b <- df[df[[treat_var]] %in% FALSE & is.finite(df[[var]]), ]
  if (nrow(a) < 2 || nrow(b) < 2) return(NA_real_)
  if (is.null(weight_var)) {
    ma <- mean(a[[var]], na.rm = TRUE)
    mb <- mean(b[[var]], na.rm = TRUE)
    va <- stats::var(a[[var]], na.rm = TRUE)
    vb <- stats::var(b[[var]], na.rm = TRUE)
  } else {
    ma <- weighted.mean(a[[var]], a[[weight_var]], na.rm = TRUE)
    mb <- weighted.mean(b[[var]], b[[weight_var]], na.rm = TRUE)
    va <- weighted_var(a[[var]], a[[weight_var]])
    vb <- weighted_var(b[[var]], b[[weight_var]])
  }
  denom <- sqrt((va + vb) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  abs(ma - mb) / denom
}

make_balance_frame <- function(df) {
  out <- df
  if ("stage3_path" %in% names(out)) {
    out$bal_stage3_creatinine_only <- as.numeric(out$stage3_path == "creatinine_only")
    out$bal_stage3_urine_output_only <- as.numeric(out$stage3_path == "urine_output_only")
    out$bal_stage3_both <- as.numeric(out$stage3_path == "both_creatinine_and_uo")
  }
  if ("gender" %in% names(out)) {
    out$bal_gender_male <- as.numeric(tolower(as.character(out$gender)) == "male")
  }
  out
}

balance_covariates <- function(df) {
  c(
    intersect(available_numeric_covariates, names(df)),
    grep("^bal_", names(df), value = TRUE)
  )
}

calc_rmst_one <- function(time, event, weight, tau) {
  ok <- is.finite(time) & !is.na(event) & is.finite(weight) & weight > 0
  if (sum(ok) < 2 || length(unique(event[ok])) == 0) return(NA_real_)
  fit <- survfit(Surv(time[ok], event[ok]) ~ 1, weights = weight[ok])
  tt <- fit$time[fit$time < tau]
  ss <- fit$surv[fit$time < tau]
  times <- c(0, tt, tau)
  surv <- c(1, ss)
  sum(diff(times) * surv)
}

run_landmark_analysis <- function(df, label) {
  landmark <- df %>%
    filter(follow_time_28d > 1) %>%
    mutate(
      early_crrt = as.logical(early_crrt),
      event_after_landmark = died_28d %in% TRUE & !is.na(days_to_death) & days_to_death > 1,
      follow_after_landmark = pmax(0.001, pmin(27, follow_time_28d - 1)),
      stage3_path = factor(stage3_path),
      gender = factor(gender)
    ) %>%
    impute_median(available_numeric_covariates)

  if (nrow(landmark) < 30 || length(unique(landmark$early_crrt)) < 2) {
    return(list(
      label = label,
      data = landmark,
      diagnostics = tibble(analysis = label, status = "insufficient_landmark_contrast", n = nrow(landmark),
                           early_n = sum(landmark$early_crrt, na.rm = TRUE)),
      result = tibble(analysis = label, status = "insufficient_landmark_contrast"),
      smd = tibble(),
      weight_distribution = tibble()
    ))
  }

  early_n <- sum(landmark$early_crrt, na.rm = TRUE)
  if (early_n < 10) {
    denom_terms <- intersect(c("age", "acutephysiologyscore", "lactate_near_t0", "creatinine_baseline"), names(landmark))
    num_terms <- "1"
  } else if (early_n < 25) {
    denom_terms <- intersect(c("age", "apachescore", "acutephysiologyscore", "creatinine_baseline", "lactate_near_t0", "stage3_path"), names(landmark))
    num_terms <- intersect(c("age", "acutephysiologyscore"), names(landmark))
  } else {
    denom_terms <- intersect(
      c("age", "apachescore", "acutephysiologyscore", "creatinine_baseline",
        "creatinine_at_kdigo3", "lactate_near_t0", "weight_kg",
        "apache_meanbp", "apache_ph", "apache_bun", "apache_creatinine",
        "apache_urine", "stage3_path"),
      names(landmark)
    )
    num_terms <- intersect(c("age", "apachescore", "acutephysiologyscore"), names(landmark))
  }

  landmark_model <- standardize_numeric(landmark, intersect(available_numeric_covariates, denom_terms))
  denom_formula <- as.formula(paste("early_crrt ~", paste(denom_terms, collapse = " + ")))
  num_formula <- if (identical(num_terms, "1") || length(num_terms) == 0) {
    early_crrt ~ 1
  } else {
    as.formula(paste("early_crrt ~", paste(num_terms, collapse = " + ")))
  }

  fit_status <- "ok"
  denom_fit <- tryCatch(
    suppressWarnings(glm(denom_formula, data = landmark_model, family = binomial())),
    error = function(e) e
  )
  num_fit <- tryCatch(
    suppressWarnings(glm(num_formula, data = landmark_model, family = binomial())),
    error = function(e) e
  )
  if (inherits(denom_fit, "error") || inherits(num_fit, "error")) {
    fit_status <- "weight_model_fallback_intercept_only"
    denom_fit <- glm(early_crrt ~ 1, data = landmark_model, family = binomial())
    num_fit <- glm(early_crrt ~ 1, data = landmark_model, family = binomial())
    denom_formula <- early_crrt ~ 1
    num_formula <- early_crrt ~ 1
  }

  p_denom <- pmin(pmax(as.numeric(predict(denom_fit, newdata = landmark_model, type = "response")), 0.01), 0.99)
  p_num <- pmin(pmax(as.numeric(predict(num_fit, newdata = landmark_model, type = "response")), 0.01), 0.99)

  analysis_df <- landmark %>%
    mutate(
      p_denom = p_denom,
      p_num = p_num,
      sw = if_else(early_crrt, p_num / p_denom, (1 - p_num) / (1 - p_denom)),
      sw_trim = pmin(sw, as.numeric(quantile(sw, TRIM_Q, na.rm = TRUE)))
    )

  weight_distribution <- bind_rows(lapply(c("early", "deferred", "overall"), function(arm) {
    d <- analysis_df
    if (arm == "early") d <- d %>% filter(early_crrt)
    if (arm == "deferred") d <- d %>% filter(!early_crrt)
    tibble(
      analysis = label,
      arm = arm,
      n = nrow(d),
      p_denom_min = min(d$p_denom, na.rm = TRUE),
      p_denom_p01 = as.numeric(quantile(d$p_denom, 0.01, na.rm = TRUE)),
      p_denom_p05 = as.numeric(quantile(d$p_denom, 0.05, na.rm = TRUE)),
      p_denom_median = median(d$p_denom, na.rm = TRUE),
      p_denom_p95 = as.numeric(quantile(d$p_denom, 0.95, na.rm = TRUE)),
      p_denom_p99 = as.numeric(quantile(d$p_denom, 0.99, na.rm = TRUE)),
      p_denom_max = max(d$p_denom, na.rm = TRUE),
      sw_p95 = as.numeric(quantile(d$sw, 0.95, na.rm = TRUE)),
      sw_p99 = as.numeric(quantile(d$sw, 0.99, na.rm = TRUE)),
      sw_max = max(d$sw, na.rm = TRUE),
      sw_trim_max = max(d$sw_trim, na.rm = TRUE)
    )
  }))

  cox_fit <- tryCatch(
    coxph(Surv(follow_after_landmark, event_after_landmark) ~ early_crrt,
          data = analysis_df, weights = sw_trim, robust = TRUE),
    error = function(e) e
  )
  if (inherits(cox_fit, "error")) {
    result <- tibble(analysis = label, status = "cox_model_failed")
    diagnostics <- tibble(
      analysis = label,
      status = "cox_model_failed",
      n = nrow(analysis_df),
      early_n = sum(analysis_df$early_crrt, na.rm = TRUE)
    )
    return(list(label = label, data = analysis_df, diagnostics = diagnostics, result = result, smd = tibble(), weight_distribution = weight_distribution))
  }

  beta <- unname(coef(cox_fit)[1])
  se <- sqrt(vcov(cox_fit)[1, 1])
  hr <- exp(beta)
  ci_low <- exp(beta - 1.96 * se)
  ci_high <- exp(beta + 1.96 * se)
  p_value <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)

  rmst_early <- calc_rmst_one(
    analysis_df$follow_after_landmark[analysis_df$early_crrt],
    analysis_df$event_after_landmark[analysis_df$early_crrt],
    analysis_df$sw_trim[analysis_df$early_crrt],
    tau = 27
  )
  rmst_deferred <- calc_rmst_one(
    analysis_df$follow_after_landmark[!analysis_df$early_crrt],
    analysis_df$event_after_landmark[!analysis_df$early_crrt],
    analysis_df$sw_trim[!analysis_df$early_crrt],
    tau = 27
  )

  balance_df <- make_balance_frame(analysis_df)
  bal_vars <- balance_covariates(balance_df)
  smd_tbl <- tibble(covariate = bal_vars) %>%
    mutate(
      smd_unweighted = vapply(covariate, function(v) smd_numeric(balance_df, v), numeric(1)),
      smd_weighted = vapply(covariate, function(v) smd_numeric(balance_df, v, "sw_trim"), numeric(1))
    )

  diagnostics <- tibble(
    analysis = label,
    status = fit_status,
    n = nrow(analysis_df),
    early_n = sum(analysis_df$early_crrt, na.rm = TRUE),
    deferred_n = sum(!analysis_df$early_crrt, na.rm = TRUE),
    event_n = sum(analysis_df$event_after_landmark, na.rm = TRUE),
    untrimmed_weight_p99 = as.numeric(quantile(analysis_df$sw, 0.99, na.rm = TRUE)),
    untrimmed_weight_max = max(analysis_df$sw, na.rm = TRUE),
    trim_q = TRIM_Q,
    trim_cap = max(analysis_df$sw_trim, na.rm = TRUE),
    trimmed_fraction = mean(analysis_df$sw > analysis_df$sw_trim, na.rm = TRUE),
    ess = sum(analysis_df$sw_trim, na.rm = TRUE)^2 / sum(analysis_df$sw_trim^2, na.rm = TRUE),
    ess_early = with(analysis_df[analysis_df$early_crrt %in% TRUE, ],
                     sum(sw_trim, na.rm = TRUE)^2 / sum(sw_trim^2, na.rm = TRUE)),
    ess_deferred = with(analysis_df[analysis_df$early_crrt %in% FALSE, ],
                        sum(sw_trim, na.rm = TRUE)^2 / sum(sw_trim^2, na.rm = TRUE)),
    max_weighted_smd = max(smd_tbl$smd_weighted, na.rm = TRUE),
    denominator_formula = paste(deparse(denom_formula), collapse = " "),
    numerator_formula = paste(deparse(num_formula), collapse = " ")
  )

  result <- tibble(
    analysis = label,
    status = "ok",
    hr = hr,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_ratio_high_to_low = ci_high / ci_low,
    p_value = p_value,
    mimic_reference_hr = MIMIC_HR_REFERENCE,
    hr_ratio_vs_mimic = hr / MIMIC_HR_REFERENCE,
    hr_ratio_inside_0_80_to_1_25 = (hr / MIMIC_HR_REFERENCE) >= 0.80 & (hr / MIMIC_HR_REFERENCE) <= 1.25,
    rmst_early_27d_after_landmark = rmst_early,
    rmst_deferred_27d_after_landmark = rmst_deferred,
    rmst_diff_early_minus_deferred = rmst_early - rmst_deferred
  )

  list(label = label, data = analysis_df, diagnostics = diagnostics, result = result, smd = smd_tbl, weight_distribution = weight_distribution)
}

analysis_high_fidelity <- recommended_cohort %>% filter(high_fidelity_hospital %in% TRUE)
analysis_all_hospitals <- recommended_cohort

landmark_high <- run_landmark_analysis(analysis_high_fidelity, "recommended_def_c_high_fidelity")
landmark_all <- run_landmark_analysis(analysis_all_hospitals, "recommended_def_c_all_hospitals")

task4_results <- bind_rows(landmark_high$result, landmark_all$result)
task4_diagnostics <- bind_rows(landmark_high$diagnostics, landmark_all$diagnostics)
task4_smd <- bind_rows(
  if (nrow(landmark_high$smd) > 0) mutate(landmark_high$smd, analysis = landmark_high$label),
  if (nrow(landmark_all$smd) > 0) mutate(landmark_all$smd, analysis = landmark_all$label)
)
task4_weight_distribution <- bind_rows(landmark_high$weight_distribution, landmark_all$weight_distribution)

if (!is.null(landmark_high$data)) {
  write_csv(landmark_high$data, file.path(data_dir, "task4_landmark_high_fidelity.csv"))
  saveRDS(landmark_high$data, file.path(data_dir, "task4_landmark_high_fidelity.rds"))
}
if (!is.null(landmark_all$data)) {
  write_csv(landmark_all$data, file.path(data_dir, "task4_landmark_all_hospitals.csv"))
  saveRDS(landmark_all$data, file.path(data_dir, "task4_landmark_all_hospitals.rds"))
}

task5_precision_gate <- task4_results %>%
  left_join(select(task4_diagnostics, analysis, n, early_n, ess, ess_early, ess_deferred, max_weighted_smd), by = "analysis") %>%
  mutate(
    early_n_gate = early_n >= 30,
    early_ess_gate = is.finite(ess_early) & ess_early >= 30 & ess_early >= 0.70 * early_n,
    balance_gate = is.finite(max_weighted_smd) & max_weighted_smd <= 0.10,
    ci_width_gate = is.finite(ci_ratio_high_to_low) & ci_ratio_high_to_low <= 4,
    point_similarity_gate = hr_ratio_inside_0_80_to_1_25 %in% TRUE,
    usable_external_validation = early_n_gate & early_ess_gate & balance_gate & ci_width_gate & point_similarity_gate,
    interpretation = case_when(
      usable_external_validation ~ "usable external validation",
      point_similarity_gate & !usable_external_validation ~ "directionally consistent but insufficient precision or balance",
      TRUE ~ "not replicated"
    )
  )

recommended_no_r2 <- recommended_obj$cohort_no_r2
secondary_outcome_censoring <- bind_rows(
  recommended_cohort %>%
    summarise(
      analysis_population = "Def C primary rescue exclusion",
      n = n(),
      early_n = sum(early_crrt, na.rm = TRUE),
      death_28d_n = sum(died_28d, na.rm = TRUE),
      discharged_alive_before_28d_n = sum(discharged_alive_before_28d, na.rm = TRUE),
      discharged_alive_before_28d_pct = round(100 * mean(discharged_alive_before_28d, na.rm = TRUE), 1)
    ),
  recommended_no_r2 %>%
    summarise(
      analysis_population = "Def C rescue exclusion without pH R2",
      n = n(),
      early_n = sum(early_crrt, na.rm = TRUE),
      death_28d_n = sum(died_28d, na.rm = TRUE),
      discharged_alive_before_28d_n = sum(discharged_alive_before_28d, na.rm = TRUE),
      discharged_alive_before_28d_pct = round(100 * mean(discharged_alive_before_28d, na.rm = TRUE), 1)
    )
)

write_csv(task1_crosstab, file.path(data_dir, "task1_vaso_lactate_crosstab.csv"))
write_csv(task0_old_row_reconciliation, file.path(data_dir, "task0_old_row_reconciliation.csv"))
write_csv(task1_sequential_funnel, file.path(data_dir, "task1_sequential_funnel.csv"))
write_csv(task1_loss_priority, file.path(data_dir, "task1_flow_loss_priority.csv"))
write_csv(task1_overlap, file.path(data_dir, "task1_excluded_overlap_indicators.csv"))
write_csv(task2_definition_comparison, file.path(data_dir, "task2_shock_definition_comparison.csv"))
write_csv(task3_high_fidelity_grid, file.path(data_dir, "task3_high_fidelity_grid.csv"))
write_csv(recommended_grid, file.path(data_dir, "task3_recommended_grid_choice.csv"))
write_csv(task4_results, file.path(data_dir, "task4_recommended_landmark_results.csv"))
write_csv(task4_diagnostics, file.path(data_dir, "task4_recommended_diagnostics.csv"))
write_csv(task4_smd, file.path(data_dir, "task4_recommended_smd.csv"))
write_csv(task4_weight_distribution, file.path(data_dir, "task4_propensity_weight_distribution.csv"))
write_csv(task5_precision_gate, file.path(data_dir, "task5_precision_gate.csv"))
write_csv(secondary_outcome_censoring, file.path(data_dir, "secondary_outcome_censoring_and_no_r2.csv"))

fmt_table <- function(df, digits = 3) {
  if (is.null(df) || nrow(df) == 0) return(character())
  x <- df
  for (nm in names(x)) {
    if (is.numeric(x[[nm]])) {
      x[[nm]] <- format(round(x[[nm]], digits), trim = TRUE, scientific = FALSE)
    } else {
      x[[nm]] <- as.character(x[[nm]])
    }
    x[[nm]][is.na(x[[nm]])] <- ""
  }
  header <- paste(names(x), collapse = " | ")
  sep <- paste(rep("---", ncol(x)), collapse = " | ")
  rows <- apply(x, 1, function(row) paste(row, collapse = " | "))
  c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |"))
}

top_gate <- task5_precision_gate %>%
  filter(analysis == "recommended_def_c_high_fidelity")
if (nrow(top_gate) == 0) top_gate <- task5_precision_gate %>% slice(1)

top_conclusion <- if (nrow(top_gate) == 0 || !isTRUE(top_gate$usable_external_validation[1])) {
  if (nrow(top_gate) > 0 && !isTRUE(top_gate$point_similarity_gate[1])) {
    paste0(
      "eICU 不能作为确认性外部验证，也不能作为方向一致的复现证据；推荐 Def C 加高保真阈值 ",
      "每院至少 ", recommended_grid$min_hospital_n,
      " 例、UO 文档比例至少 ", recommended_grid$min_uo_doc_prop,
      "、升压药速率文档比例至少 ", recommended_grid$min_vaso_rate_prop,
      " 后，landmark early n 为 ", top_gate$early_n[1],
      "，加权后最大 SMD 为 ", round(top_gate$max_weighted_smd[1], 3),
      "，HR ratio 为 ", round(top_gate$hr_ratio_vs_mimic[1], 3),
      "，只能作为外部数据可行性/迁移性探索。"
    )
  } else {
  paste0(
    "eICU 仍不能作为确认性外部验证；推荐 Def C 加高保真阈值 ",
    "每院至少 ", recommended_grid$min_hospital_n,
    " 例、UO 文档比例至少 ", recommended_grid$min_uo_doc_prop,
    "、升压药速率文档比例至少 ", recommended_grid$min_vaso_rate_prop,
    "，但 landmark early n 为 ", top_gate$early_n[1],
    "，加权后最大 SMD 为 ", round(top_gate$max_weighted_smd[1], 3),
    "，只能写作样本和精度不足的支持性分析。"
  )
  }
} else {
  paste0(
    "eICU 可作为有限外部验证；推荐 Def C 加高保真阈值 ",
    "每院至少 ", recommended_grid$min_hospital_n,
    " 例、UO 文档比例至少 ", recommended_grid$min_uo_doc_prop,
    "、升压药速率文档比例至少 ", recommended_grid$min_vaso_rate_prop,
    "，landmark early n 为 ", top_gate$early_n[1],
    "，加权后最大 SMD 为 ", round(top_gate$max_weighted_smd[1], 3),
    "。"
  )
}

report <- c(
  "# eICU Diagnose and Rescue Report",
  "",
  "## Top-Line Conclusion",
  "",
  top_conclusion,
  "",
  "Important correction: the previous scaffold's `septic_shock_proxy_at_t0 = 819` was a row count inflated by one-to-many joins. After deduplication by `patientunitstayid`, the old executable shock definition contains 425 unique ICU stays and 422 first stays per patient.",
  "",
  "## What This Script Did",
  "",
  "- Rebuilt the eICU clean physiologic KDIGO3 base without using CRRT as a t0 criterion.",
  "- Decomposed the first physiologic KDIGO3 to previous septic-shock-proxy loss.",
  "- Tested shock definitions A to D and the previous executable definition.",
  "- Tested high-fidelity hospital thresholds across UO documentation and vasopressor numeric-rate documentation.",
  "- Re-ran a 24-hour landmark stabilized-weight analysis for the recommended Def C option.",
  "- Added a precision gate before calling any eICU result external validation.",
  "",
  "## Task 1. Loss From First Physiologic KDIGO3 to Old Shock Proxy",
  "",
  "Old row-count reconciliation:",
  "",
  fmt_table(task0_old_row_reconciliation),
  "",
  "Vasopressor and lactate cross-tabulation in the 5932 first physiologic KDIGO3 patients:",
  "",
  fmt_table(task1_crosstab),
  "",
  "Sequential eligibility and old shock-proxy funnel:",
  "",
  fmt_table(task1_sequential_funnel),
  "",
  "Priority loss decomposition among patients excluded by the previous executable shock proxy:",
  "",
  fmt_table(task1_loss_priority),
  "",
  "Non-mutually-exclusive indicators among excluded patients:",
  "",
  fmt_table(task1_overlap),
  "",
  "Interpretation: the priority table is mutually exclusive; the overlap table shows how much missingness coexists with other exclusion mechanisms.",
  "",
  "## Task 2. Shock Definition Comparison",
  "",
  fmt_table(task2_definition_comparison %>% select(-stage3_path_composition), digits = 2),
  "",
  "Stage-3 path composition is in `task2_shock_definition_comparison.csv` because it is too wide for the main report.",
  "",
  "Recommended main option: Def C. It keeps the vasopressor-supported septic-shock phenotype but does not force lactate measurement or numeric-rate documentation, both of which are known weak points in eICU.",
  "",
  "## Task 3. High-Fidelity Threshold Sensitivity",
  "",
  "Recommended threshold row:",
  "",
  fmt_table(recommended_grid, digits = 2),
  "",
  "Full threshold grid is saved as `task3_high_fidelity_grid.csv`.",
  "",
  "## Task 4. Landmark/IPCW Re-Run Under Recommended Def C",
  "",
  "Results:",
  "",
  fmt_table(task4_results, digits = 3),
  "",
  "Weight and positivity diagnostics:",
  "",
  fmt_table(task4_diagnostics, digits = 3),
  "",
  "Propensity and stabilized-weight distribution by arm:",
  "",
  fmt_table(task4_weight_distribution, digits = 3),
  "",
  "Maximum SMD after weighting must be read before interpreting the HR. If balance fails, the HR is descriptive only.",
  "",
  "## Task 5. Precision Gate",
  "",
  fmt_table(task5_precision_gate %>% select(analysis, early_n, ess_early, max_weighted_smd, ci_ratio_high_to_low, hr_ratio_vs_mimic, early_n_gate, early_ess_gate, balance_gate, ci_width_gate, point_similarity_gate, usable_external_validation, interpretation), digits = 3),
  "",
  "Manuscript sentence if the similarity gate and precision/balance gates fail: The eICU analysis did not meet the prespecified information, balance, precision, or point-estimate similarity requirements for confirmatory external validation; therefore, it was interpreted as an exploratory feasibility and transportability analysis rather than successful replication.",
  "",
  "Manuscript sentence if point estimates are compatible but precision/balance gates fail: Although the eICU point estimate was directionally compatible with the MIMIC estimate, the validation analysis failed prespecified precision and balance criteria; therefore, the eICU analysis was considered insufficiently powered rather than confirmatory replication.",
  "",
  "## Secondary Fixes",
  "",
  fmt_table(secondary_outcome_censoring, digits = 2),
  "",
  "The eICU outcome is hospital death within 28 days from t0 with discharge censoring, whereas MIMIC uses true 28-day all-cause mortality from death records. This difference must be stated in the manuscript.",
  "",
  "## Assumptions and Exact Logic",
  "",
  "- t0: first physiologic KDIGO3 by creatinine or urine-output reconstruction only; CRRT/RRT does not define t0.",
  "- Vasopressor near t0 window: t0 minus 360 minutes to t0 plus 60 minutes.",
  "- Lactate near t0 window: t0 minus 720 minutes to t0 plus 360 minutes.",
  "- Early CRRT: first continuous renal replacement therapy string from t0 through 24 hours after t0.",
  "- Primary rescue exclusion: R1 hyperkalemia, R2 any-specimen pH, R2b bicarbonate/HCO3, or R3 BUN trigger at or before t0.",
  "- No-R2 sensitivity: drops the pH rescue component because eICU pH lacks reliable specimen type.",
  "- High-fidelity hospitals: hospital-level minimum cohort size, urine documentation in the first 24 hours, and vasopressor numeric-rate documentation.",
  "- External validation gate: point-estimate similarity is not sufficient. The estimate must also have enough early-arm information, acceptable balance, and acceptable CI width.",
  "",
  "## Outputs",
  "",
  paste0("- PostgreSQL schema: `", schema_name, "`"),
  "- `data/task1_vaso_lactate_crosstab.csv`",
  "- `data/task0_old_row_reconciliation.csv`",
  "- `data/task1_sequential_funnel.csv`",
  "- `data/task1_flow_loss_priority.csv`",
  "- `data/task2_shock_definition_comparison.csv`",
  "- `data/task3_high_fidelity_grid.csv`",
  "- `data/task4_recommended_landmark_results.csv`",
  "- `data/task4_recommended_diagnostics.csv`",
  "- `data/task4_recommended_smd.csv`",
  "- `data/task4_propensity_weight_distribution.csv`",
  "- `data/task5_precision_gate.csv`",
  "- `data/secondary_outcome_censoring_and_no_r2.csv`"
)

writeLines(report, file.path(report_dir, "eicu_diagnose_rescue_report.md"))

message("\n=== Top-line conclusion ===")
message(top_conclusion)
message("\nSaved diagnosis/rescue outputs under: ", out_dir)
