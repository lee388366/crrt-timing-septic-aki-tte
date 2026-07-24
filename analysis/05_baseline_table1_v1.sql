-- ============================================================
-- 05) cohort_baseline_v1 — Baseline characteristics (Table 1)
-- Project  : MIMIC-IV CCW — CRRT timing in Septic Shock ∩ AKI (KDIGO3)
-- Output   : data_extract_crrt.cohort_baseline_v1
-- Source   : data_extract_crrt."301_ss_aki3_t0_no_rrt_cohort"  (n = 789)
--
-- 设计说明
-- ─────────────────────────────────────────────────────────────
-- t0        : max(ICU intime, 首次达到KDIGO 3期 − 6h)
-- 基线窗口   : (t0 − 6h, t0]，用于生命体征、实验室、评分
-- 治疗窗口   : [ICU intime, t0)，用于累积液体/药物暴露
-- 观测策略   : LATERAL + ORDER BY charttime DESC LIMIT 1
--             取窗口内最后一次观测值（非均值），与原始管线一致
--
-- 变量分组
-- ─────────────────────────────────────────────────────────────
-- 标识符     : subject_id, hadm_id, stay_id, intime, outtime, t0_time
-- ICU背景    : first_careunit (ICU类型), admission_year
-- 人口学     : age, gender, race, insurance, weight_kg
-- 严重程度   : cci, sofa_baseline, gcs_baseline, apsiii_score
-- 生命体征   : hr, sbp, dbp, map, rr, spo2, temp_c
-- 血常规     : wbc, hb, plt
-- 生化       : na, cl, k, hco3, scr, bun
-- 肝功能     : ast, alt
-- 血气/乳酸  : ph, lactate, pao2, fio2, pf_ratio
-- 凝血       : inr
-- 合并症     : dm, htn, chf, ckd, liver_disease, malignancy, copd
-- 基线治疗   : mv_baseline, vaso_baseline, ne_eq_baseline
-- 液体/药物  : crystalloid_ml_kg, colloid_flag, nephrotoxic_flag
-- AKI描述    : uo_mlkgph_t0_last6h, aki_pathway, scr_baseline
--
-- 数据来源（MIMIC-IV v3.1）
-- ─────────────────────────────────────────────────────────────
-- mimiciv_icu.icustays              → first_careunit
-- mimiciv_hosp.admissions           → admittime（用于scr_baseline窗口）
-- mimiciv_derived.weight_durations  → weight_kg
-- mimiciv_derived.charlson          → cci及各合并症组件
-- mimiciv_derived.vitalsign         → hr/sbp/dbp/mbp/rr/spo2/temperature
-- mimiciv_derived.chemistry         → na/cl/k/bicarbonate/creatinine/bun
-- mimiciv_derived.bg                → ph/lactate/po2/fio2
-- mimiciv_derived.enzyme            → ast/alt（独立于chemistry表）
-- mimiciv_derived.coagulation       → inr
-- mimiciv_derived.complete_blood_count → wbc/hemoglobin/platelet
-- mimiciv_derived.gcs               → gcs_motor/gcs_verbal/gcs_eyes
-- mimiciv_derived.apsiii            → apsiii（ICU入室第一天评分）
-- mimiciv_derived.sofa              → sofa_24hours（滚动24h窗口）
-- mimiciv_derived.ventilation       → ventilation_status = 'InvasiveVent'
-- mimiciv_derived.vasoactive_agent  → 各升压药剂量字段
-- mimiciv_derived.norepinephrine_equivalent_dose → ne当量（无charttime）
-- mimiciv_derived.urine_output      → urineoutput（mL/event）
-- mimiciv_derived.kdigo_stages      → aki_stage_creat / aki_stage_uo
-- mimiciv_icu.inputevents           → 液体输入（itemid列表见下）
-- mimiciv_hosp.prescriptions        → 药物处方（drug字段文本匹配）
-- mimiciv_hosp.diagnoses_icd        → ICD诊断码（用于高血压）
--
-- Navicat运行方式
-- ─────────────────────────────────────────────────────────────
-- 1. 在已连接的数据库下新建查询
-- 2. Ctrl+A 全选 → Ctrl+R 运行整个文件
-- 3. 跑完后单独选中最后一行 SELECT 再运行查看结果
-- 4. 索引见 05_baseline_indexes.sql（可选，提升后续查询速度）
-- ============================================================

DROP TABLE IF EXISTS data_extract_crrt.cohort_baseline_v1;

CREATE TABLE data_extract_crrt.cohort_baseline_v1 AS

-- ── CTE: 基础队列 ──────────────────────────────────────────────────────────────
-- 来源：步骤03最终队列（脓毒性休克 ∩ KDIGO 3期，排除t0前已行RRT，n=789）
-- baseline_start = t0 - 6h，作为基线窗口下界，避免每个LATERAL重复计算
WITH cohort AS (
  SELECT
    stay_id, subject_id, hadm_id,
    intime, outtime,
    age, gender, race, insurance,
    los_icu, kdigo3_time, t0_time,
    t0_time - INTERVAL '6 hour' AS baseline_start
  FROM data_extract_crrt."301_ss_aki3_t0_no_rrt_cohort"
)

SELECT
  -- ── 标识符 ──────────────────────────────────────────────────────────────────
  c.subject_id,    -- 患者ID（跨住院唯一）
  c.hadm_id,       -- 住院ID
  c.stay_id,       -- ICU住院ID
  c.intime,        -- ICU入室时间
  c.outtime,       -- ICU出室时间
  c.t0_time,       -- 研究时间零点 = max(ICU intime, KDIGO3时间 − 6h)

  -- ── ICU背景 ─────────────────────────────────────────────────────────────────
  icu_det.first_careunit,                        -- ICU类型（MICU/SICU/CCU等）
  EXTRACT(YEAR FROM c.intime)::INT AS admission_year,  -- 入ICU年份

  -- ── 人口学 ──────────────────────────────────────────────────────────────────
  c.age,           -- 入院年龄（岁）
  c.gender,        -- 性别（M/F）
  c.race,          -- 种族
  c.insurance,     -- 保险类型（Medicare/Medicaid/Private等）
  w.weight_kg,     -- 体重（kg），取t0前最近一次记录

  -- ── 严重程度评分 ─────────────────────────────────────────────────────────────
  -- CCI：Charlson合并症指数（范围0-37，越高合并症越多）
  cc.charlson_comorbidity_index AS cci,
  -- SOFA：取含t0的24h滚动窗口内最后一个有效评分（范围0-24）
  sf.sofa_baseline,
  -- GCS：取基线窗口内三分项（运动+言语+睁眼）均非空时的最后一次总分（范围3-15）
  gcs.gcs_baseline,
  -- APSIII：ICU入室第一天计算，作为入ICU初始严重程度的辅助描述
  -- 注意：大多数患者的APSIII时间早于t0（t0为KDIGO3发生后），仅供参考
  ap.apsiii_score,

  -- ── 生命体征（基线窗口内最后一次观测值） ───────────────────────────────────
  v.hr,      -- 心率（次/分）
  v.sbp,     -- 收缩压（mmHg）
  v.dbp,     -- 舒张压（mmHg）
  v.map,     -- 平均动脉压（mmHg）；源字段为vitalsign.mbp
  v.rr,      -- 呼吸频率（次/分）
  v.spo2,    -- 脉氧饱和度（%）
  v.temp_c,  -- 体温（℃）；mimiciv_derived.vitalsign已转换为摄氏度

  -- ── 血常规（基线窗口内最后一次） ────────────────────────────────────────────
  cb.wbc,    -- 白细胞计数（×10⁹/L）
  cb.hb,     -- 血红蛋白（g/dL）；源字段为hemoglobin
  cb.plt,    -- 血小板计数（×10⁹/L）；源字段为platelet

  -- ── 生化（基线窗口内最后一次） ──────────────────────────────────────────────
  ch.na,     -- 血钠（mEq/L）
  ch.cl,     -- 血氯（mEq/L）
  ch.k,      -- 血钾（mEq/L）
  ch.hco3,   -- 碳酸氢根（mEq/L）；源字段为bicarbonate
  ch.scr,    -- 血肌酐（mg/dL）；注意这是t0时值，非基线值
  ch.bun,    -- 血尿素氮（mg/dL）

  -- ── 肝功能（基线窗口内最后一次，来自enzyme表，独立于chemistry） ────────────
  enz.ast,   -- 天冬氨酸转氨酶（U/L）
  enz.alt,   -- 丙氨酸转氨酶（U/L）

  -- ── 血气 & 乳酸（基线窗口内最后一次动脉血气） ──────────────────────────────
  b.ph,      -- 动脉血pH
  b.lactate, -- 血乳酸（mmol/L）
  b.pao2,    -- 动脉氧分压（mmHg）；源字段为bg.po2
  b.fio2,    -- 吸入氧浓度（百分比格式，25–100；mimiciv_derived.bg以百分比存储）
  -- P/F比值 = PaO2(mmHg) / (FiO2%/100)；先将百分比转小数再计算
  -- Berlin ARDS分层：<200重度，200-300中度，>300轻度/正常
  -- 注：bg.fio2在MIMIC-IV中以百分比存储（25-100），必须除以100.0转为小数
  CASE
    WHEN b.pao2 IS NOT NULL AND b.fio2 IS NOT NULL AND b.fio2 > 0
    THEN b.pao2 / NULLIF(b.fio2 / 100.0, 0)
    ELSE NULL
  END AS pf_ratio,

  -- ── 凝血（基线窗口内最后一次） ──────────────────────────────────────────────
  cg.inr,    -- 国际标准化比值

  -- ── 合并症 ──────────────────────────────────────────────────────────────────
  -- 来源1：Charlson衍生表的各组件（二值变量，1=有，0=无）
  -- GREATEST处理有/无并发症的糖尿病分类，取更严重的编码
  GREATEST(COALESCE(cc.diabetes_without_cc, 0),
           COALESCE(cc.diabetes_with_cc,    0)) AS dm,           -- 糖尿病
  -- 来源2：ICD诊断码（高血压不在Charlson指数中，需单独提取）
  COALESCE(htn.htn, 0)                          AS htn,          -- 高血压
  COALESCE(cc.congestive_heart_failure,     0)  AS chf,          -- 充血性心衰
  -- 注：步骤01已排除ESRD，故此处renal_disease=1代表非终末期CKD
  COALESCE(cc.renal_disease,               0)   AS ckd,          -- 慢性肾脏病（非ESRD）
  -- 取轻/重度肝病中更严重的编码
  GREATEST(COALESCE(cc.mild_liver_disease,  0),
           COALESCE(cc.severe_liver_disease, 0)) AS liver_disease, -- 慢性肝病
  -- 取局限性/转移性肿瘤中更严重的编码
  GREATEST(COALESCE(cc.malignant_cancer,    0),
           COALESCE(cc.metastatic_solid_tumor, 0)) AS malignancy, -- 恶性肿瘤
  COALESCE(cc.chronic_pulmonary_disease,   0)   AS copd,          -- 慢阻肺

  -- ── 基线治疗状态（t0时刻） ───────────────────────────────────────────────────
  -- mv_baseline    : 1=t0时正在有创机械通气，0=否
  -- vaso_baseline  : 1=t0时正在使用任意升压药，0=否（设计上所有患者均为1）
  -- ne_eq_baseline : t0时去甲肾上腺素当量（μg/kg/min），取基线窗口内最后有效值
  sp.mv_baseline,
  sp.vaso_baseline,
  sp.ne_eq_baseline,

  -- ── 液体与药物暴露（ICU入室到t0，即[intime, t0)） ───────────────────────────
  -- crystalloid_ml_kg : 累积等/低张晶体液总量（mL）÷ 体重（kg）
  --   itemid覆盖：NaCl 0.9%/0.45%/3%、乳酸林格液、5%葡萄糖、碳酸氢钠等
  --   仅统计amountuom='mL'的记录，排除速率记录
  CASE
    WHEN w.weight_kg > 0 THEN cr.crystalloid_ml / w.weight_kg
    ELSE NULL
  END AS crystalloid_ml_kg,
  -- colloid_flag : 1=曾使用人工胶体（羟乙基淀粉/右旋糖酐），0=否
  --   注：白蛋白为天然胶体，已排除（itemid 220862/220864）
  co.colloid_flag,
  -- nephrotoxic_flag : 1=曾使用肾毒性药物，0=否
  --   覆盖：氨基糖苷类、万古霉素、NSAIDs（酮咯酸/布洛芬等）、两性霉素B、多黏菌素
  nt.nephrotoxic_flag,

  -- ── AKI描述性变量 ────────────────────────────────────────────────────────────
  -- uo_mlkgph_t0_last6h : t0前6h尿量（mL/kg/h），KDIGO少尿判断关键指标
  --   KDIGO 3期少尿定义：<0.3 mL/kg/h持续≥24h 或 <0.5 mL/kg/h持续≥12h
  uo_t0.uo_mlkgph_t0_last6h,
  -- aki_pathway : AKI达到KDIGO 3期的驱动路径
  --   'oliguria'   = 仅少尿标准达到3期
  --   'creatinine' = 仅肌酐标准达到3期
  --   'both'       = 两个标准均达到3期
  --   'unknown'    = kdigo_stages在kdigo3_time附近无有效记录（极少见）
  --   注：此处取kdigo3_time前后12h窗口的MAX分期（非t0），修复了原始版本
  --       因t0比kdigo3_time早6h导致全部显示'unknown'的bug
  akip.aki_pathway,
  -- scr_baseline : [入院时间, t0]期间最低肌酐值（mg/dL）
  --   作为AKI前基础肾功能的代理指标，同时是28天肾功能恢复终点的分母
  --   公式：renal_recovery = alive_28d AND scr_28d ≤ 1.5 × scr_baseline
  scr_bl.scr_baseline

FROM cohort c

-- ── JOIN 1: ICU类型 ────────────────────────────────────────────────────────────
-- mimiciv_icu.icustays是first_careunit的正确来源（非icustay_detail）
LEFT JOIN mimiciv_icu.icustays icu_det
  ON icu_det.stay_id = c.stay_id

-- ── JOIN 2: 住院信息（提供admittime，用于scr_baseline的时间窗口下界） ──────────
LEFT JOIN mimiciv_hosp.admissions adm
  ON adm.hadm_id = c.hadm_id

-- ── JOIN 3: 体重 ───────────────────────────────────────────────────────────────
-- 取t0时刻之前最近一次有效体重记录
-- weight_durations中每行代表一个体重有效时间段（starttime到endtime）
LEFT JOIN LATERAL (
  SELECT wd.weight AS weight_kg
  FROM mimiciv_derived.weight_durations wd
  WHERE wd.stay_id = c.stay_id
    AND wd.starttime <= c.t0_time  -- 体重记录时间在t0之前
    AND wd.weight IS NOT NULL
  ORDER BY wd.starttime DESC       -- 取最接近t0的记录
  LIMIT 1
) w ON true

-- ── JOIN 4: Charlson合并症指数及各组件 ────────────────────────────────────────
-- charlson表按hadm_id存储，一次住院对应一行
-- 同时提供总指数和各个合并症的二值编码（用于SELECT中展开个体合并症）
LEFT JOIN mimiciv_derived.charlson cc
  ON cc.hadm_id = c.hadm_id

-- ── JOIN 5: 生命体征 ───────────────────────────────────────────────────────────
-- 基线窗口(t0-6h, t0]内最后一次记录
-- vitalsign.mbp = 平均动脉压，输出时重命名为map
-- temperature已由MIMIC衍生表转换为摄氏度
LEFT JOIN LATERAL (
  SELECT
    v.heart_rate AS hr, v.sbp, v.dbp, v.mbp AS map,
    v.resp_rate AS rr, v.spo2, v.temperature AS temp_c
  FROM mimiciv_derived.vitalsign v
  WHERE v.stay_id = c.stay_id
    AND v.charttime >  c.baseline_start  -- 严格大于（开区间下界）
    AND v.charttime <= c.t0_time         -- 小于等于（闭区间上界）
  ORDER BY v.charttime DESC
  LIMIT 1
) v ON true

-- ── JOIN 6: 生化 ───────────────────────────────────────────────────────────────
-- 注：chemistry表用hadm_id关联（非stay_id）
-- ast/alt不在此表，来自独立的enzyme表（见JOIN 8）
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

-- ── JOIN 7: 血气 & 乳酸 ────────────────────────────────────────────────────────
-- bg表用hadm_id关联
-- PaO2在bg表中字段名为po2（非pao2），此处别名为pao2便于下游使用
LEFT JOIN LATERAL (
  SELECT b.ph, b.lactate, b.po2 AS pao2, b.fio2
  FROM mimiciv_derived.bg b
  WHERE b.hadm_id = c.hadm_id
    AND b.charttime >  c.baseline_start
    AND b.charttime <= c.t0_time
  ORDER BY b.charttime DESC
  LIMIT 1
) b ON true

-- ── JOIN 8: 肝功能酶学 ─────────────────────────────────────────────────────────
-- enzyme是独立于chemistry的衍生表，hadm_id关联
LEFT JOIN LATERAL (
  SELECT e.ast, e.alt
  FROM mimiciv_derived.enzyme e
  WHERE e.hadm_id = c.hadm_id
    AND e.charttime >  c.baseline_start
    AND e.charttime <= c.t0_time
  ORDER BY e.charttime DESC
  LIMIT 1
) enz ON true

-- ── JOIN 9: 凝血 ───────────────────────────────────────────────────────────────
LEFT JOIN LATERAL (
  SELECT cg.inr
  FROM mimiciv_derived.coagulation cg
  WHERE cg.hadm_id = c.hadm_id
    AND cg.charttime >  c.baseline_start
    AND cg.charttime <= c.t0_time
  ORDER BY cg.charttime DESC
  LIMIT 1
) cg ON true

-- ── JOIN 10: 血常规 ────────────────────────────────────────────────────────────
-- 字段名：hemoglobin（非hgb）、platelet（非plt），此处别名对齐输出
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
-- 仅取三分项（运动/言语/睁眼）均有记录的行，避免部分缺失导致计算偏低
-- GCS总分 = gcs_motor（1-6）+ gcs_verbal（1-5）+ gcs_eyes（1-4），范围3-15
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
-- apsiii表每个stay_id通常只有一行（ICU第一天），LATERAL+LIMIT 1防止极少数多行情况
-- 注意：APSIII基于ICU入室第一天，对于t0较晚的患者（KDIGO3出现较迟），
--       APSIII反映的是入ICU时的严重程度，不完全代表t0时状态
LEFT JOIN LATERAL (
  SELECT a.apsiii AS apsiii_score
  FROM mimiciv_derived.apsiii a
  WHERE a.stay_id = c.stay_id
  LIMIT 1
) ap ON true

-- ── JOIN 13: 机械通气 / 升压药 / NE当量（合并为一个LATERAL减少扫描次数） ────────
-- mv_baseline    : EXISTS查询检测t0时是否有活跃的有创通气（InvasiveVent）
-- vaso_baseline  : EXISTS查询检测t0时是否有任意升压药活跃（设计上全部为1）
-- ne_eq_baseline : 子查询取基线窗口内最后一个有效NE当量值
--   重要：norepinephrine_equivalent_dose表只有starttime/endtime，无charttime
--   时间重叠判断：starttime < t0 AND endtime > baseline_start（区间相交）
--   COALESCE(endtime, t0_time) 处理endtime为NULL（即仍在持续输注）的情况
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

-- ── JOIN 14: SOFA评分 ──────────────────────────────────────────────────────────
-- sofa表每行为一个24h滚动窗口（starttime到endtime）
-- 使用endtime进行窗口匹配：取endtime落在基线窗口内的最后一个有效评分
-- 注：sofa_24hours = 该24h窗口内的SOFA总分（0-24）
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

-- ── JOIN 15: 高血压（ICD诊断码） ───────────────────────────────────────────────
-- 高血压不在Charlson指数中，需从入院诊断ICD码单独提取
-- ICD-9: 401（原发性高血压）、402-405（继发/合并症高血压）
-- ICD-10: I10（原发性高血压）、I11-I13（高血压心脏/肾脏病）、I15（继发性高血压）
-- 使用正则表达式替代多条LIKE语句，性能更优
-- EXISTS返回布尔值，::int 转换为 1/0
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

-- ── JOIN 16: 晶体液累积量 [intime, t0) ────────────────────────────────────────
-- 统计ICU入室到t0之间所有晶体液输入总量（mL），SELECT中除以体重得mL/kg
-- 仅统计 amountuom = 'mL' 的记录（排除速率/浓度等单位的行）
-- ItemID覆盖范围（MIMIC-IV MetaVision）：
--   225158 生理盐水(0.9%)    225828 乳酸林格液
--   225159 0.45%NaCl        226364 手术室晶体液（汇总项）
--   225161 3%NaCl           226375 PACU晶体液（汇总项）
--   225944 无菌水            227533 8.4%碳酸氢钠IV
--   220949 5%葡萄糖(D5W)    220950 D5+0.45%NaCl
--   225797 自由水            225825 D5+乳酸林格液
--   225827 D5+0.5%NaCl     225941 D5+0.25%NaCl
LEFT JOIN LATERAL (
  SELECT COALESCE(SUM(ie.amount), 0) AS crystalloid_ml
  FROM mimiciv_icu.inputevents ie
  WHERE ie.stay_id = c.stay_id
    AND ie.starttime >= c.intime   -- 闭区间，从ICU入室开始
    AND ie.starttime <  c.t0_time  -- 开区间，不含t0时刻
    AND ie.itemid IN (
        225158, 225159, 225161, 225828, 225944,
        226364, 226375, 227533, 220949, 220950,
        225797, 225825, 225827, 225941
    )
    AND ie.amountuom = 'mL'
) cr ON true

-- ── JOIN 17: 人工胶体暴露标志 [intime, t0) ────────────────────────────────────
-- 人工胶体：羟乙基淀粉（HES）和右旋糖酐，不含白蛋白（天然胶体）
--   225174 Hetastarch/Hespan/Hextend 6%（羟乙基淀粉）
--   225795 Dextran-40（右旋糖酐40）
--   225796 Dextran-70（右旋糖酐70）
-- 注：现代ICU中HES使用极少（2013年FDA警告后），该指标在MIMIC中流行率可能很低
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

-- ── JOIN 18: 肾毒性药物暴露标志 [intime, t0) ──────────────────────────────────
-- 在住院处方表（prescriptions，hadm_id级别）中匹配drug字段文本
-- ~* = 不区分大小写的正则匹配，等价于原先12条 ILIKE 语句，性能更优
-- 覆盖药物类别：
--   氨基糖苷类：gentamicin（庆大霉素）、tobramycin（妥布霉素）、
--              amikacin（阿米卡星）、streptomycin（链霉素）
--   糖肽类：   vancomycin（万古霉素）
--   NSAIDs：   ketorolac（酮咯酸）、ibuprofen（布洛芬）、
--              indomethacin（吲哚美辛）、naproxen（萘普生）
--   抗真菌：   amphotericin（两性霉素B）
--   多肽类：   colistin（多黏菌素E）、polymyxin（多黏菌素B）
-- 未纳入：IV造影剂（需链接procedureevents/radiology，可按需补充）
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

-- ── JOIN 19: t0前6h尿量 ────────────────────────────────────────────────────────
-- 统计(t0-6h, t0]内所有尿量事件之和（mL），除以体重（kg）和6小时，得mL/kg/h
-- 内嵌子查询获取体重（避免LATERAL之间的依赖问题）
-- NULLIF(weight, 0) 防止体重为0时除零错误
-- 若6h内无尿量记录，SUM为NULL，COALESCE转为0（代表无尿）
-- 临床意义：KDIGO少尿标准 <0.5 mL/kg/h；严重少尿 <0.3 mL/kg/h
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

-- ── JOIN 20: AKI达到KDIGO 3期的驱动路径 ───────────────────────────────────────
-- 取kdigo_stages在kdigo3_time前后各12h/6h窗口内的最大分期
-- 之所以用kdigo3_time（而非t0）：t0 = kdigo3_time - 6h，在t0时患者尚未达到KDIGO 3
-- 所以若用t0查分期，全部为<3（即原始版本的'unknown' bug）
-- ±12h窗口的理由：基于尿量的KDIGO 3需要12-24h累积少尿时间，其在kdigo_stages
--   表中的记录时间可能比肌酐标准滞后数小时，宽窗口确保同时捕捉到两种标准
LEFT JOIN LATERAL (
  SELECT
    CASE
      WHEN MAX(ks.aki_stage_creat) >= 3
       AND MAX(ks.aki_stage_uo)    >= 3 THEN 'both'        -- 肌酐+尿量双标准
      WHEN MAX(ks.aki_stage_uo)    >= 3 THEN 'oliguria'    -- 仅尿量标准
      WHEN MAX(ks.aki_stage_creat) >= 3 THEN 'creatinine'  -- 仅肌酐标准
      ELSE                                   'unknown'     -- 极少见，需排查
    END AS aki_pathway
  FROM mimiciv_derived.kdigo_stages ks
  WHERE ks.stay_id  = c.stay_id
    AND ks.charttime >  c.kdigo3_time - INTERVAL '12 hour'
    AND ks.charttime <= c.kdigo3_time + INTERVAL '6 hour'
) akip ON true

-- ── JOIN 21: 基线肌酐（入院至t0的最低值） ─────────────────────────────────────
-- 定义：[admittime, t0]期间所有有效肌酐记录的最小值
-- 作为AKI前基础肾功能的代理指标（true baseline SCr通常无直接记录）
-- 临床意义：28天肾功能恢复的判断标准 = 存活 AND 28天SCr ≤ 1.5 × scr_baseline
-- 使用admittime（住院入院时间）而非intime（ICU入室时间），
--   因为患者可能在普通病房期间已有肌酐记录，可提供更早的基线值
LEFT JOIN LATERAL (
  SELECT MIN(ch2.creatinine) AS scr_baseline
  FROM mimiciv_derived.chemistry ch2
  WHERE ch2.hadm_id   = c.hadm_id
    AND ch2.charttime >= adm.admittime  -- 从住院入院时间开始
    AND ch2.charttime <= c.t0_time      -- 到t0为止
    AND ch2.creatinine IS NOT NULL
) scr_bl ON true
;

-- ── 结果表索引（提升后续R分析中的JOIN效率） ───────────────────────────────────
CREATE INDEX IF NOT EXISTS idx_baseline_v1_stay
  ON data_extract_crrt.cohort_baseline_v1 (stay_id);
CREATE INDEX IF NOT EXISTS idx_baseline_v1_hadm
  ON data_extract_crrt.cohort_baseline_v1 (hadm_id);
CREATE INDEX IF NOT EXISTS idx_baseline_v1_subject
  ON data_extract_crrt.cohort_baseline_v1 (subject_id);

ANALYZE data_extract_crrt.cohort_baseline_v1;

-- ====== 查看结果：若结果栏无数据，单独选中下面一行再运行 ======
SELECT * FROM data_extract_crrt.cohort_baseline_v1;
