# Rescue-indication diagnostics for the strict first-KDIGO3 anchor.
# This script is intentionally limited to cohort counts and rescue diagnostics:
# no CCW, IPCW, Cox, RMST, or causal modeling is run here.
suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(tidyr)
  library(readr)
})

out_dir <- "analysis/strictK3_rescue_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
  charttime AS t0,
  aki_stage_creat AS aki_stage_creat_t0,
  aki_stage_uo AS aki_stage_uo_t0,
  aki_stage_crrt AS aki_stage_crrt_t0,
  uo_rt_12hr,
  uo_rt_24hr
FROM mimiciv_derived.kdigo_stages
WHERE aki_stage = 3
ORDER BY stay_id, charttime
")
run("CREATE INDEX ON k3(stay_id)")

run("
CREATE TEMP TABLE sepsis AS
SELECT
  k3.stay_id,
  bool_or(s.sepsis3 AND s.sofa_time <= k3.t0) AS sepsis_ok
FROM k3
LEFT JOIN mimiciv_derived.sepsis3 s USING (stay_id)
GROUP BY k3.stay_id
")
run("CREATE INDEX ON sepsis(stay_id)")

run("
CREATE TEMP TABLE vaso AS
SELECT
  k3.stay_id,
  count(v.*) > 0 AS vaso_ok,
  avg(
    coalesce(v.norepinephrine, 0)
    + coalesce(v.epinephrine, 0)
    + coalesce(v.dopamine, 0) / 100.0
    + coalesce(v.phenylephrine, 0) * 0.1
    + coalesce(v.vasopressin, 0) * 2.5
  ) AS nee_at_t0
FROM k3
LEFT JOIN mimiciv_derived.vasoactive_agent v
  ON v.stay_id = k3.stay_id
 AND v.starttime <= k3.t0
 AND v.endtime >= k3.t0 - interval '6 hours'
GROUP BY k3.stay_id
")
run("CREATE INDEX ON vaso(stay_id)")

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
  k3.t0,
  i.subject_id,
  i.hadm_id,
  i.stay_id,
  i.intime,
  i.outtime,
  extract(epoch from (i.outtime - i.intime)) / 3600 AS icu_los_hours,
  p.anchor_age AS age,
  ad.dischtime,
  ad.deathtime,
  cb.scr_baseline AS creatinine_baseline,
  cr.crrt_first_time,
  v.nee_at_t0,
  k3.aki_stage_creat_t0,
  k3.aki_stage_uo_t0,
  k3.aki_stage_crrt_t0,
  k3.uo_rt_12hr,
  k3.uo_rt_24hr
FROM k3
JOIN mimiciv_icu.icustays i USING (stay_id)
JOIN mimiciv_hosp.patients p USING (subject_id)
JOIN mimiciv_hosp.admissions ad USING (hadm_id)
LEFT JOIN mimiciv_derived.creatinine_baseline cb USING (hadm_id)
LEFT JOIN crrt cr USING (stay_id)
JOIN sepsis s USING (stay_id)
JOIN vaso v USING (stay_id)
WHERE p.anchor_age >= 18
  AND s.sepsis_ok
  AND v.vaso_ok
")
run("CREATE INDEX ON source(stay_id)")
run("CREATE INDEX ON source(hadm_id, t0)")

source_counts <- dbGetQuery(con, "
WITH
  n_k3 AS (SELECT count(*) AS n FROM k3),
  n_adult AS (
    SELECT count(*) AS n
    FROM k3
    JOIN mimiciv_icu.icustays i USING (stay_id)
    JOIN mimiciv_hosp.patients p USING (subject_id)
    WHERE p.anchor_age >= 18
  ),
  n_sepsis_vaso AS (SELECT count(*) AS n FROM source),
  n_no_prior_rrt AS (
    SELECT count(*) AS n FROM source
    WHERE crrt_first_time IS NULL OR crrt_first_time >= t0
  ),
  n_no_esrd AS (
    SELECT count(*) AS n FROM source src
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
    SELECT count(*) AS n FROM source src
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
    SELECT count(*) AS n FROM source src
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
    SELECT count(*) AS n FROM source src
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
SELECT 'K3 stays' AS step, n FROM n_k3
UNION ALL SELECT 'Adult K3 stays', n FROM n_adult
UNION ALL SELECT 'Sepsis-3 by t0 + vasopressor window', n FROM n_sepsis_vaso
UNION ALL SELECT 'No RRT before t0', n FROM n_no_prior_rrt
UNION ALL SELECT 'No ESRD/chronic dialysis', n FROM n_no_esrd
UNION ALL SELECT 'ICU LOS >=24h', n FROM n_los
UNION ALL SELECT 'Nonmissing baseline essentials', n FROM n_baseline
UNION ALL SELECT 'Positive 28d follow-up', n FROM n_followup
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
run("CREATE INDEX ON pre_rescue_all_stays(stay_id)")
run("CREATE INDEX ON pre_rescue_all_stays(hadm_id, t0)")

run("
CREATE TEMP TABLE pre_rescue AS
SELECT DISTINCT ON (subject_id) *
FROM pre_rescue_all_stays
ORDER BY subject_id, intime
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
  coalesce(bun.bun_n, 0) AS r3_bun_n,
  coalesce(bun.r3_trigger, false) AS r3_trigger,
  (pr.uo_rt_12hr IS NOT NULL OR pr.uo_rt_24hr IS NOT NULL) AS r4_uo_measured,
  (coalesce(pr.uo_rt_24hr < 0.3, false) OR coalesce(pr.uo_rt_12hr = 0, false)) AS r4_trigger,
  coalesce(pr.uo_rt_12hr = 0, false) AS r4_anuria12_proxy,
  coalesce(pr.uo_rt_24hr < 0.3, false) AS r4_severe_oliguria24,
  coalesce(pr.crrt_first_time BETWEEN pr.t0 AND pr.t0 + interval '24 hours', false) AS early_crrt
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
    count(*) FILTER (WHERE c.bun IS NOT NULL) AS bun_n,
    bool_or(c.bun >= 100) AS r3_trigger
  FROM mimiciv_derived.chemistry c
  WHERE c.hadm_id = pr.hadm_id
    AND c.charttime > pr.t0 - interval '24 hours'
    AND c.charttime <= pr.t0
    AND c.bun IS NOT NULL
) bun ON true
")
run("CREATE INDEX ON rescue_flags(stay_id)")

run("
CREATE TEMP TABLE rescue_final AS
SELECT
  rf.*,
  sf.sofa_24hours AS sofa_t0,
  (sf.sofa_24hours BETWEEN 8 AND 14) AS sofa8_14,
  (r1_trigger OR r2_trigger OR r3_trigger) AS rescue_r1_r3,
  (r1_trigger OR r2_trigger OR r3_trigger OR r4_trigger) AS rescue_r1_r4,
  (r1_trigger AND r1_k_n = 1) AS r1_single_unrechecked_k_high
FROM rescue_flags rf
LEFT JOIN LATERAL (
  SELECT sofa_24hours
  FROM mimiciv_derived.sofa s
  WHERE s.stay_id = rf.stay_id
    AND s.starttime <= rf.t0 + interval '1 hour'
  ORDER BY s.starttime DESC
  LIMIT 1
) sf ON true
")

pt <- dbGetQuery(con, "SELECT * FROM rescue_final")
write_csv(pt, file.path(out_dir, "strictK3_rescue_patient_level.csv"))

source_counts <- source_counts %>%
  mutate(n = as.numeric(n)) %>%
  mutate(
    previous_n = lag(n),
    excluded_from_previous = if_else(is.na(previous_n), NA_real_, previous_n - n)
  )
write_csv(source_counts, file.path(out_dir, "strictK3_funnel_pre_rescue.csv"))

summarise_version <- function(dat, version, exclude_var) {
  keep <- if (is.null(exclude_var)) rep(TRUE, nrow(dat)) else !dat[[exclude_var]]
  x <- dat[keep, , drop = FALSE]
  eq <- x$sofa8_14 %in% TRUE
  tibble(
    version = version,
    final_n = nrow(x),
    early_n = sum(x$early_crrt, na.rm = TRUE),
    early_pct = round(100 * mean(x$early_crrt, na.rm = TRUE), 1),
    sofa8_14_n = sum(eq, na.rm = TRUE),
    sofa8_14_early_n = sum(eq & x$early_crrt, na.rm = TRUE),
    sofa8_14_early_pct = round(100 * mean(x$early_crrt[eq], na.rm = TRUE), 1),
    rescue_excluded_total = nrow(dat) - nrow(x),
    rescue_excluded_sofa8_14 = sum((dat$sofa8_14 %in% TRUE) & !keep, na.rm = TRUE)
  )
}

version_summary <- bind_rows(
  summarise_version(pt, "no_rescue_exclusion", NULL),
  summarise_version(pt, "exclude_R1_R3", "rescue_r1_r3"),
  summarise_version(pt, "exclude_R1_R4_primary", "rescue_r1_r4")
) %>%
  mutate(
    incremental_vs_R1_R3_total = if_else(
      version == "exclude_R1_R4_primary",
      rescue_excluded_total - rescue_excluded_total[version == "exclude_R1_R3"],
      NA_integer_
    ),
    incremental_vs_R1_R3_sofa8_14 = if_else(
      version == "exclude_R1_R4_primary",
      rescue_excluded_sofa8_14 - rescue_excluded_sofa8_14[version == "exclude_R1_R3"],
      NA_integer_
    )
  )
write_csv(version_summary, file.path(out_dir, "strictK3_rescue_version_summary.csv"))

indicator_summary <- tibble(
  indicator = c("R1_hyperkalemia", "R2_arterial_acidemia", "R3_azotemia", "R4_anuria_or_severe_oliguria"),
  lookback_window = c("6h", "6h", "24h", "12h/24h derived UO rates at t0"),
  denominator_n = nrow(pt),
  measured_n = c(
    sum(pt$r1_k_n > 0, na.rm = TRUE),
    sum(pt$r2_abg_ph_n > 0, na.rm = TRUE),
    sum(pt$r3_bun_n > 0, na.rm = TRUE),
    sum(pt$r4_uo_measured, na.rm = TRUE)
  ),
  trigger_n = c(
    sum(pt$r1_trigger, na.rm = TRUE),
    sum(pt$r2_trigger, na.rm = TRUE),
    sum(pt$r3_trigger, na.rm = TRUE),
    sum(pt$r4_trigger, na.rm = TRUE)
  )
) %>%
  mutate(
    missing_n = denominator_n - measured_n,
    missing_pct = round(100 * missing_n / denominator_n, 1),
    trigger_pct_all = round(100 * trigger_n / denominator_n, 1),
    trigger_pct_measured = round(100 * trigger_n / measured_n, 1)
  )
write_csv(indicator_summary, file.path(out_dir, "strictK3_rescue_indicator_summary.csv"))

uo_loop_summary <- bind_rows(
  pt %>%
    summarise(
      population = "pre_rescue_total",
      n = n(),
      kdigo3_by_uo_n = sum(aki_stage_uo_t0 == 3, na.rm = TRUE),
      kdigo3_by_uo_and_r4_excluded_n = sum(aki_stage_uo_t0 == 3 & r4_trigger, na.rm = TRUE),
      r4_trigger_n = sum(r4_trigger, na.rm = TRUE),
      r4_anuria12_proxy_n = sum(r4_anuria12_proxy, na.rm = TRUE),
      r4_severe_oliguria24_n = sum(r4_severe_oliguria24, na.rm = TRUE)
    ),
  pt %>%
    filter(sofa8_14 %in% TRUE) %>%
    summarise(
      population = "SOFA_8_14_pre_rescue",
      n = n(),
      kdigo3_by_uo_n = sum(aki_stage_uo_t0 == 3, na.rm = TRUE),
      kdigo3_by_uo_and_r4_excluded_n = sum(aki_stage_uo_t0 == 3 & r4_trigger, na.rm = TRUE),
      r4_trigger_n = sum(r4_trigger, na.rm = TRUE),
      r4_anuria12_proxy_n = sum(r4_anuria12_proxy, na.rm = TRUE),
      r4_severe_oliguria24_n = sum(r4_severe_oliguria24, na.rm = TRUE)
    )
)
write_csv(uo_loop_summary, file.path(out_dir, "strictK3_R4_uo_loop_summary.csv"))

k_single <- tibble(
  metric = c(
    "R1 triggers",
    "R1 triggers based on a single K measurement in the 6h window",
    "R1 triggers with >=2 K measurements in the 6h window",
    "Patients with nonarterial pH records ignored for R2"
  ),
  n = c(
    sum(pt$r1_trigger, na.rm = TRUE),
    sum(pt$r1_single_unrechecked_k_high, na.rm = TRUE),
    sum(pt$r1_trigger & pt$r1_k_n >= 2, na.rm = TRUE),
    sum(pt$r2_nonarterial_ph_n > 0, na.rm = TRUE)
  )
)
write_csv(k_single, file.path(out_dir, "strictK3_special_rescue_checks.csv"))

trigger_overlap <- pt %>%
  transmute(
    R1 = r1_trigger,
    R2 = r2_trigger,
    R3 = r3_trigger,
    R4 = r4_trigger
  ) %>%
  mutate(pattern = paste0(if_else(R1, "R1", ""),
                          if_else(R2, "+R2", ""),
                          if_else(R3, "+R3", ""),
                          if_else(R4, "+R4", "")),
         pattern = if_else(pattern == "", "none", sub("^\\+", "", pattern))) %>%
  count(pattern, name = "n") %>%
  arrange(desc(n))
write_csv(trigger_overlap, file.path(out_dir, "strictK3_rescue_trigger_overlap.csv"))

report <- c(
  "# Strict KDIGO3 Rescue Diagnostics",
  "",
  "Scope: counts and rescue/positivity diagnostics only. No CCW, IPCW, Cox, RMST, Fine-Gray, or causal modeling was run.",
  "",
  "## Implementation notes",
  "",
  "- Anchor: first observed KDIGO stage 3 charttime (`mimiciv_derived.kdigo_stages.aki_stage = 3`).",
  "- R1: potassium >= 6.5 mmol/L in `(t0 - 6h, t0]` from `mimiciv_derived.chemistry`.",
  "- R2: arterial blood gas pH <= 7.15 in `(t0 - 6h, t0]` from `mimiciv_derived.bg` with `specimen = 'ART.'`; venous/mixed/central venous pH is not used for triggering.",
  "- R3: BUN >= 100 mg/dL in `(t0 - 24h, t0]` from `mimiciv_derived.chemistry`.",
  "- R4: derived urine-output proxy at t0: `uo_rt_12hr == 0` for anuria >=12h proxy, or `uo_rt_24hr < 0.3` mL/kg/h for severe oliguria >=24h.",
  "- Missing measurements in the lookback window are treated as not meeting that rescue indication, and missingness is reported per indication.",
  "- First qualifying ICU stay is fixed by `subject_id` after source eligibility and before rescue exclusion.",
  "",
  "## Pre-rescue funnel",
  "",
  "| step | n | excluded from previous |",
  "|---|---:|---:|",
  apply(source_counts, 1, function(x) sprintf("| %s | %s | %s |", x[["step"]], x[["n"]], ifelse(is.na(x[["excluded_from_previous"]]), "", x[["excluded_from_previous"]]))),
  "",
  "## Rescue versions",
  "",
  "| version | final n | early n | early % | SOFA 8-14 n | SOFA 8-14 early n | SOFA 8-14 early % | rescue excluded total | rescue excluded SOFA 8-14 | incremental R4 total | incremental R4 SOFA 8-14 |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
  apply(version_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    x[["version"]], x[["final_n"]], x[["early_n"]], x[["early_pct"]],
    x[["sofa8_14_n"]], x[["sofa8_14_early_n"]], x[["sofa8_14_early_pct"]],
    x[["rescue_excluded_total"]], x[["rescue_excluded_sofa8_14"]],
    x[["incremental_vs_R1_R3_total"]], x[["incremental_vs_R1_R3_sofa8_14"]]
  )),
  "",
  "## Per-indication triggers and missingness",
  "",
  "| indicator | window | measured n | missing n | missing % | trigger n | trigger % all | trigger % measured |",
  "|---|---|---:|---:|---:|---:|---:|---:|",
  apply(indicator_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    x[["indicator"]], x[["lookback_window"]], x[["measured_n"]],
    x[["missing_n"]], x[["missing_pct"]], x[["trigger_n"]],
    x[["trigger_pct_all"]], x[["trigger_pct_measured"]]
  )),
  "",
  "## R4 / urine-output loop check",
  "",
  "| population | n | KDIGO3 by UO n | KDIGO3 by UO and R4 excluded n | R4 trigger n | R4 anuria12 proxy n | R4 severe oliguria24 n |",
  "|---|---:|---:|---:|---:|---:|---:|",
  apply(uo_loop_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s |",
    x[["population"]], x[["n"]], x[["kdigo3_by_uo_n"]],
    x[["kdigo3_by_uo_and_r4_excluded_n"]], x[["r4_trigger_n"]],
    x[["r4_anuria12_proxy_n"]], x[["r4_severe_oliguria24_n"]]
  )),
  "",
  "## Special checks",
  "",
  "| metric | n |",
  "|---|---:|",
  apply(k_single, 1, function(x) sprintf("| %s | %s |", x[["metric"]], x[["n"]])),
  "",
  "## Files",
  "",
  "- `strictK3_rescue_patient_level.csv`",
  "- `strictK3_funnel_pre_rescue.csv`",
  "- `strictK3_rescue_version_summary.csv`",
  "- `strictK3_rescue_indicator_summary.csv`",
  "- `strictK3_R4_uo_loop_summary.csv`",
  "- `strictK3_special_rescue_checks.csv`",
  "- `strictK3_rescue_trigger_overlap.csv`"
)
writeLines(report, file.path(out_dir, "strictK3_rescue_diagnostics_report.md"))

print(version_summary)
print(indicator_summary)
print(uo_loop_summary)
print(k_single)
