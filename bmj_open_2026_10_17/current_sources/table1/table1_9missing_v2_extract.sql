-- 回补母板 Table1 的 9 个缺失变量（HR/RR/SpO2/温度/pH/MV/HTN/晶体/肾毒素）到 v2 SOFA8-14
-- 忠实照搬母板 05_baseline_table1_v2.sql 的 LATERAL 逻辑；唯一改动=队列源换成 v2 clean-t0 + 按 early_crrt 分臂。
-- 只读聚合（percentile_cont + 计数/均值/标准差），不建表、不改科学逻辑。
WITH cohort AS (
  SELECT b.stay_id, b.subject_id, b.hadm_id,
         b.icu_intime AS intime, b.kdigo3_time,
         b.t0 AS t0_time, b.t0 - INTERVAL '6 hour' AS baseline_start,
         b.early_crrt
  FROM data_extract_crrt_strictk3_v2.baseline b
  WHERE b.sofa_equipoise IS TRUE          -- SOFA 8-14，n=1502
),
pat AS (
  SELECT c.stay_id, c.early_crrt,
    v.hr, v.rr, v.spo2, v.temp_c, bg.ph, sp.mv_baseline, htn.htn,
    CASE WHEN w.weight_kg > 0 THEN cr.crystalloid_ml / w.weight_kg ELSE NULL END AS crystalloid_ml_kg,
    nt.nephrotoxic_flag
  FROM cohort c
  -- JOIN3 weight（母板逐字）
  LEFT JOIN LATERAL (SELECT wd.weight AS weight_kg FROM mimiciv_derived.weight_durations wd
    WHERE wd.stay_id=c.stay_id AND wd.starttime<=c.t0_time AND wd.weight IS NOT NULL
    ORDER BY wd.starttime DESC LIMIT 1) w ON true
  -- JOIN5 vitalsign（母板逐字）
  LEFT JOIN LATERAL (SELECT vs.heart_rate AS hr, vs.resp_rate AS rr, vs.spo2, vs.temperature AS temp_c
    FROM mimiciv_derived.vitalsign vs
    WHERE vs.stay_id=c.stay_id AND vs.charttime>c.baseline_start AND vs.charttime<=c.t0_time
    ORDER BY vs.charttime DESC LIMIT 1) v ON true
  -- JOIN7 bg（母板逐字，仅取 ph）
  LEFT JOIN LATERAL (SELECT bb.ph FROM mimiciv_derived.bg bb
    WHERE bb.hadm_id=c.hadm_id AND bb.charttime>c.baseline_start AND bb.charttime<=c.t0_time
    ORDER BY bb.charttime DESC LIMIT 1) bg ON true
  -- JOIN13 ventilation（母板逐字）
  LEFT JOIN LATERAL (SELECT EXISTS(SELECT 1 FROM mimiciv_derived.ventilation vv
    WHERE vv.stay_id=c.stay_id AND vv.starttime<c.t0_time AND COALESCE(vv.endtime,c.t0_time)>c.baseline_start
      AND vv.ventilation_status='InvasiveVent')::int AS mv_baseline FROM (SELECT 1) _) sp ON true
  -- JOIN15 htn ICD（母板逐字）
  LEFT JOIN LATERAL (SELECT EXISTS(SELECT 1 FROM mimiciv_hosp.diagnoses_icd d
    WHERE d.hadm_id=c.hadm_id AND ((d.icd_version=9 AND d.icd_code ~ '^40[1-5]')
      OR (d.icd_version=10 AND (d.icd_code ~ '^I1[0-3]' OR d.icd_code ~ '^I15'))))::int AS htn
    FROM (SELECT 1) _) htn ON true
  -- JOIN16 crystalloid（母板逐字 itemid 列表）
  LEFT JOIN LATERAL (SELECT COALESCE(SUM(ie.amount),0) AS crystalloid_ml FROM mimiciv_icu.inputevents ie
    WHERE ie.stay_id=c.stay_id AND ie.starttime>=c.intime AND ie.starttime<c.t0_time
      AND ie.itemid IN (225158,225159,225161,225828,225944,226364,226375,227533,220949,220950,225797,225825,225827,225941)
      AND ie.amountuom='mL') cr ON true
  -- JOIN18 nephrotoxic（母板逐字药名正则）
  LEFT JOIN LATERAL (SELECT EXISTS(SELECT 1 FROM mimiciv_hosp.prescriptions p
    WHERE p.hadm_id=c.hadm_id AND p.starttime>=c.intime AND p.starttime<c.t0_time
      AND p.drug ~* 'gentamicin|tobramycin|amikacin|streptomycin|vancomycin|ketorolac|ibuprofen|indomethacin|naproxen|amphotericin|colistin|polymyxin')::int AS nephrotoxic_flag
    FROM (SELECT 1) _) nt ON true
)
SELECT (CASE WHEN early_crrt THEN 'early' WHEN NOT early_crrt THEN 'nocrrt' ELSE 'all' END) AS arm,
  count(*) n, count(hr) hr_n, percentile_cont(0.5) within group(order by hr) hr_med, /* ...q1/q3/mean/sd for hr,rr,spo2,temp,ph,crystalloid... */
  sum(mv_baseline) mv_n, sum(htn) htn_n, sum(nephrotoxic_flag) neph_n
FROM pat GROUP BY GROUPING SETS ((early_crrt),());
-- 说明：实际执行版含各连续变量 q1/q3/mean/sd（见运行记录）。母板对应块：05_baseline_table1_v2.sql JOIN 3/5/7/13/15/16/18。
