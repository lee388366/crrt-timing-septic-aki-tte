# =============================================================================
# Rich baseline Table 1 for the MMR master template layout, rebuilt on strictK3 v2.
# Population: primary SOFA 8-14 equipoise cohort (N = 1,502).
# Arms:  Early CRRT within 24 h (n = 146)  vs  No CRRT within 24 h (n = 1,356).
# Columns: Characteristic | Overall (N=1502) | Early CRRT (n=146) | No CRRT (n=1356) | SMD | p-value
#
# METHOD POLICY (per task): reuse EXISTING v2 data + old table1 code; never invent
# a new variable-extraction method; variables that cannot be obtained from ready v2
# data / documented existing code are reported honestly in a "missing" manifest and
# are NOT fabricated in the main table.
#
# Reused code:
#   - strictK3_v2_table1_mmr_spec.R : NEE vasopressin U/h -> U/min correction,
#     first_day_weight join, admissions.race grouping, fmt/SMD/p helper functions.
#   - 10_table1_results.R / 11_write_result_table1.R : master Table-1 variable set,
#     Cohen's d (continuous) + pooled-proportion (binary) SMD, Wilcoxon / chi2 / Fisher.
#
# Anchor cross-checks (must align):
#   - Overall column vs results/strictK3_v2_table1_overall_sofa8_14.csv
#   - shared by-arm rows vs results/strictK3_v2_table1_sofa8_14.csv
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
})

base_dir    <- Sys.getenv("CRRT_ANALYSIS_ROOT")
strict_dir  <- file.path(base_dir, "strictK3_primary_v2")
data_dir    <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")

out_dir <- file.path(Sys.getenv("CRRT_MANUSCRIPT_ROOT"), "table1_output")
csv_out     <- file.path(out_dir, "results",  "table1_rich_sofa8_14_v2.csv")
audit_out   <- file.path(out_dir, "results",  "table1_rich_sofa8_14_v2_bestbeffort_audit.csv")
missing_out <- file.path(out_dir, "results",  "table1_rich_sofa8_14_v2_MISSING.md")

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimiciv31", host = "localhost", port = 5432,
  user = "postgres", password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

# ---------------------------------------------------------------------------
# 1. Cohort + augmentation (reuses mmr_spec join logic verbatim in spirit)
# ---------------------------------------------------------------------------
baseline <- read_csv(file.path(data_dir, "baseline.csv"), show_col_types = FALSE)

augment <- dbGetQuery(con, "
WITH base AS (
  SELECT stay_id, hadm_id, t0
  FROM data_extract_crrt_strictk3_v2.baseline
  WHERE sofa_t0 BETWEEN 8 AND 14
),
nee AS (
  SELECT
    b.stay_id,
    avg(
      coalesce(v.norepinephrine, 0)
      + coalesce(v.epinephrine, 0)
      + coalesce(v.dopamine, 0) / 100.0
      + coalesce(v.phenylephrine, 0) * 0.1
      + (coalesce(v.vasopressin, 0) / 60.0) * 2.5   -- vasopressin U/h -> U/min before 2.5x
    ) AS nee_at_t0_corrected
  FROM base b
  JOIN mimiciv_derived.vasoactive_agent v
    ON v.stay_id = b.stay_id
   AND v.starttime <= b.t0
   AND v.endtime   >= b.t0 - interval '6 hours'
  GROUP BY b.stay_id
)
SELECT
  b.stay_id,
  ad.race,
  fdw.weight   AS weight_kg,
  fds.sofa     AS sofa_admission,
  fdg.gcs_min  AS gcs_min,
  ch.congestive_heart_failure,
  ch.renal_disease,
  ch.chronic_pulmonary_disease,
  (CASE WHEN ch.diabetes_without_cc = 1 OR ch.diabetes_with_cc = 1 THEN 1 ELSE 0 END) AS dm_flag,
  (CASE WHEN ch.mild_liver_disease  = 1 OR ch.severe_liver_disease = 1 THEN 1 ELSE 0 END) AS liver_flag,
  (CASE WHEN ch.malignant_cancer    = 1 OR ch.metastatic_solid_tumor = 1 THEN 1 ELSE 0 END) AS malig_flag,
  n.nee_at_t0_corrected
FROM base b
LEFT JOIN mimiciv_hosp.admissions        ad  ON ad.hadm_id  = b.hadm_id
LEFT JOIN mimiciv_derived.first_day_weight fdw ON fdw.stay_id = b.stay_id
LEFT JOIN mimiciv_derived.first_day_sofa   fds ON fds.stay_id = b.stay_id
LEFT JOIN mimiciv_derived.first_day_gcs    fdg ON fdg.stay_id = b.stay_id
LEFT JOIN mimiciv_derived.charlson         ch  ON ch.hadm_id  = b.hadm_id
LEFT JOIN nee n ON n.stay_id = b.stay_id
")

# Labs at t0: closest (last) non-null value within (t0 - 6 h, t0], per analyte.
# Standard MIMIC derived tables; verified to reproduce the rich overall anchor exactly
# (WBC/Hb/Plt/INR nonmissing counts identical; BUN/HCO3/K within 1).
labs <- dbGetQuery(con, "
WITH c AS (SELECT stay_id, hadm_id, t0 FROM data_extract_crrt_strictk3_v2.baseline WHERE sofa_t0 BETWEEN 8 AND 14)
SELECT c.stay_id,
  (SELECT x.bun         FROM mimiciv_derived.chemistry x           WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.bun         IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) bun,
  (SELECT x.bicarbonate FROM mimiciv_derived.chemistry x           WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.bicarbonate IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) bicarbonate,
  (SELECT x.potassium   FROM mimiciv_derived.chemistry x           WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.potassium   IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) potassium,
  (SELECT x.wbc         FROM mimiciv_derived.complete_blood_count x WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.wbc         IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) wbc,
  (SELECT x.hemoglobin  FROM mimiciv_derived.complete_blood_count x WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.hemoglobin  IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) hemoglobin,
  (SELECT x.platelet    FROM mimiciv_derived.complete_blood_count x WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.platelet    IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) platelet,
  (SELECT x.inr         FROM mimiciv_derived.coagulation x          WHERE x.hadm_id=c.hadm_id AND x.charttime> c.t0-interval '6 hours' AND x.charttime<=c.t0 AND x.inr         IS NOT NULL ORDER BY x.charttime DESC LIMIT 1) inr
FROM c")

race_group <- function(x) {
  x <- toupper(ifelse(is.na(x), "", x))
  dplyr::case_when(
    grepl("^WHITE", x) ~ "White",
    grepl("^BLACK", x) ~ "Black/African American",
    grepl("^ASIAN", x) ~ "Asian",
    grepl("HISPANIC|LATINO", x) ~ "Hispanic/Latino",
    TRUE ~ "Other/Unknown"
  )
}

df <- baseline %>%
  filter(sofa_t0 >= 8, sofa_t0 <= 14) %>%
  left_join(augment, by = "stay_id") %>%
  left_join(labs, by = "stay_id") %>%
  mutate(
    early_crrt = early_crrt %in% TRUE,
    died_28d   = died_28d %in% TRUE,
    female     = gender == "F",
    race_grp   = factor(race_group(race),
                        levels = c("White","Black/African American","Asian",
                                   "Hispanic/Latino","Other/Unknown")),
    dm         = dm_flag == 1,
    htn        = NA,                 # not a Charlson component; ICD def unrecoverable -> missing
    chf        = congestive_heart_failure == 1,
    ckd        = renal_disease == 1,
    liver      = liver_flag == 1,
    malig      = malig_flag == 1,
    copd       = chronic_pulmonary_disease == 1
  )

stopifnot(nrow(df) == 1502, sum(df$early_crrt) == 146, sum(!df$early_crrt) == 1356)

# ---------------------------------------------------------------------------
# 2. Formatting / SMD / p helpers  (verbatim from strictK3_v2_table1_mmr_spec.R)
# ---------------------------------------------------------------------------
fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

fmt_p <- function(p, suffix) {
  if (is.na(p)) return("")
  label <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  paste0(label, suffix)
}

fmt_cont <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), q[2], q[1], q[3])
}

fmt_cat <- function(x) {
  denom <- sum(!is.na(x)); if (denom == 0) return("")
  n <- sum(x %in% TRUE, na.rm = TRUE)
  sprintf("%s (%.1f)", fmt_n(n), 100 * n / denom)
}

smd_cont <- function(x, g) {
  a <- x[g %in% TRUE]; b <- x[g %in% FALSE]
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  pooled <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b)) / (length(a)+length(b)-2))
  if (!is.finite(pooled) || pooled == 0) return(0)
  abs(mean(a) - mean(b)) / pooled
}

smd_bin <- function(x, g) {
  x <- x %in% TRUE; a <- x[g %in% TRUE]; b <- x[g %in% FALSE]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  p1 <- mean(a, na.rm = TRUE); p0 <- mean(b, na.rm = TRUE)
  p <- (p1 + p0) / 2; denom <- sqrt(p * (1 - p))
  if (!is.finite(denom) || denom == 0) return(0)
  abs(p1 - p0) / denom
}

p_cont <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (length(unique(g[ok])) < 2) return("")
  suppressWarnings(fmt_p(wilcox.test(x[ok] ~ g[ok])$p.value, "\u1d43"))  # a = Wilcoxon
}

p_cat <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  tab <- table(x[ok], g[ok])
  if (nrow(tab) < 2 || ncol(tab) < 2) return("")
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    p <- if (all(dim(tab) == c(2, 2))) fisher.test(tab)$p.value
         else { set.seed(20260717); fisher.test(tab, simulate.p.value = TRUE, B = 20000)$p.value }
    fmt_p(p, "\u1d9c")   # c = Fisher exact
  } else {
    suppressWarnings(fmt_p(chisq.test(tab, correct = FALSE)$p.value, "\u1d47"))  # b = chi-square
  }
}

overall <- df; early <- df %>% filter(early_crrt); deferred <- df %>% filter(!early_crrt)
NA_CELL  <- "NA (\u672a\u56de\u8865\uff0c\u89c1\u7f3a\u5931\u6e05\u5355)"

crow <- function(label, var, digits) {
  tibble(Characteristic = label,
         Overall  = fmt_cont(overall[[var]], digits),
         Early    = fmt_cont(early[[var]], digits),
         `No CRRT`= fmt_cont(deferred[[var]], digits),
         SMD      = sprintf("%.3f", smd_cont(overall[[var]], overall$early_crrt)),
         `p-value`= p_cont(overall[[var]], overall$early_crrt))
}
brow <- function(label, vall, vearly, vdef, pvec) {
  tibble(Characteristic = label,
         Overall  = fmt_cat(vall),
         Early    = fmt_cat(vearly),
         `No CRRT`= fmt_cat(vdef),
         SMD      = sprintf("%.3f", smd_bin(vall, overall$early_crrt)),
         `p-value`= p_cat(pvec, overall$early_crrt))
}
srow <- function(label) tibble(Characteristic = label, Overall="", Early="", `No CRRT`="", SMD="", `p-value`="")
miss_row <- function(label) tibble(Characteristic = label, Overall=NA_CELL, Early=NA_CELL,
                                   `No CRRT`=NA_CELL, SMD="\u2014", `p-value`="\u2014")

race_header <- tibble(Characteristic = "Race/ethnicity, n (%)", Overall="", Early="",
                      `No CRRT`="", SMD="", `p-value`= p_cat(overall$race_grp, overall$early_crrt))
race_rows <- bind_rows(lapply(levels(overall$race_grp), function(lv)
  brow(paste0("  ", lv), overall$race_grp==lv, early$race_grp==lv, deferred$race_grp==lv, rep(lv, nrow(overall))) %>%
    mutate(`p-value` = "")))

# ---------------------------------------------------------------------------
# 3. Assemble master-template row set (Tier-1 verified + Tier-2 missing)
# ---------------------------------------------------------------------------
table1 <- bind_rows(
  srow("Demographics"),
  crow("Age, years", "age", 1),
  brow("Female sex, n (%)", overall$female, early$female, deferred$female, overall$female),
  race_header, race_rows,
  crow("Body weight, kg", "weight_kg", 1),

  srow("Illness severity"),
  crow("SOFA at ICU admission", "sofa_admission", 0),
  crow("SOFA at t0", "sofa_t0", 0),
  crow("APACHE III (APS III)", "apsiii", 0),
  crow("GCS (minimum)", "gcs_min", 0),
  crow("Charlson comorbidity index", "cci", 0),

  srow("Vital signs at t0"),
  crow("MAP, mmHg", "map_t0", 1),           # baseline map_t0 (matches by-arm slim; see note on overall)
  miss_row("Heart rate, beats/min"),
  miss_row("Respiratory rate, breaths/min"),
  miss_row("Temperature, \u00b0C"),
  miss_row("SpO2, %"),

  srow("Laboratory values at t0"),
  crow("Serum creatinine (t0), mg/dL", "creatinine_t0", 2),
  crow("Baseline SCr, mg/dL", "creatinine_baseline", 2),
  crow("BUN, mg/dL", "bun", 0),
  crow("Lactate, mmol/L", "lactate_t0", 2),
  miss_row("Arterial pH"),
  crow("Bicarbonate, mmol/L", "bicarbonate", 0),
  crow("Potassium, mmol/L", "potassium", 2),
  crow("White blood cell count, 10^9/L", "wbc", 1),
  crow("Hemoglobin, g/dL", "hemoglobin", 1),
  crow("Platelets, 10^9/L", "platelet", 0),
  crow("INR", "inr", 2),

  srow("Organ support at t0"),
  miss_row("Mechanical ventilation, n (%)"),
  crow("Norepinephrine-equivalent dose, \u00b5g/kg/min", "nee_at_t0_corrected", 2),
  crow("Urine output (last 6 h), mL/kg/h", "uo_rate_t0", 2),
  miss_row("Cumulative crystalloid, mL/kg"),
  miss_row("Nephrotoxin exposure, n (%)"),

  srow("Comorbidities"),
  brow("Diabetes mellitus, n (%)", overall$dm, early$dm, deferred$dm, overall$dm),
  miss_row("Hypertension, n (%)"),
  brow("Congestive heart failure, n (%)", overall$chf, early$chf, deferred$chf, overall$chf),
  brow("Chronic kidney disease, n (%)", overall$ckd, early$ckd, deferred$ckd, overall$ckd),
  brow("Liver disease, n (%)", overall$liver, early$liver, deferred$liver, overall$liver),
  brow("Malignancy, n (%)", overall$malig, early$malig, deferred$malig, overall$malig),
  brow("COPD, n (%)", overall$copd, early$copd, deferred$copd, overall$copd),

  srow("Outcome (unadjusted)"),
  brow("Hospital-recorded 28-day mortality, n (%)",
       overall$died_28d, early$died_28d, deferred$died_28d, overall$died_28d)
)

table1_export <- table1 %>%
  rename(`Overall (N=1502)`        = Overall,
         `Early CRRT (n=146)`      = Early,
         `No CRRT (n=1356)`        = `No CRRT`)
write_csv(table1_export, csv_out)

# ---------------------------------------------------------------------------
# 4. Best-effort audit for Tier-2 vars (standard MIMIC derived tables, t0-6h window).
#    NOT anchor-matched -> reference only, kept OUT of the main table on purpose.
# ---------------------------------------------------------------------------
be <- dbGetQuery(con, "
WITH c AS (SELECT stay_id, hadm_id, t0 FROM data_extract_crrt_strictk3_v2.baseline WHERE sofa_t0 BETWEEN 8 AND 14),
vs AS (SELECT c.stay_id, avg(v.heart_rate) hr, avg(v.resp_rate) rr, avg(v.spo2) spo2, avg(v.temperature) tc
       FROM c JOIN mimiciv_derived.vitalsign v ON v.stay_id=c.stay_id
         AND v.charttime> c.t0 - interval '6 hours' AND v.charttime<=c.t0 GROUP BY c.stay_id),
chem AS (SELECT c.stay_id, ch.bun, ch.potassium, ch.bicarbonate,
                row_number() over (partition by c.stay_id order by ch.charttime desc) rn
         FROM c JOIN mimiciv_derived.chemistry ch ON ch.hadm_id=c.hadm_id
           AND ch.charttime> c.t0 - interval '6 hours' AND ch.charttime<=c.t0),
cbc AS (SELECT c.stay_id, b.wbc, b.hemoglobin, b.platelet,
               row_number() over (partition by c.stay_id order by b.charttime desc) rn
        FROM c JOIN mimiciv_derived.complete_blood_count b ON b.hadm_id=c.hadm_id
          AND b.charttime> c.t0 - interval '6 hours' AND b.charttime<=c.t0),
cg AS (SELECT c.stay_id, co.inr, row_number() over (partition by c.stay_id order by co.charttime desc) rn
       FROM c JOIN mimiciv_derived.coagulation co ON co.hadm_id=c.hadm_id
         AND co.charttime> c.t0 - interval '6 hours' AND co.charttime<=c.t0),
gas AS (SELECT c.stay_id, g.ph, row_number() over (partition by c.stay_id order by g.charttime desc) rn
        FROM c JOIN mimiciv_derived.bg g ON g.hadm_id=c.hadm_id
          AND g.charttime> c.t0 - interval '6 hours' AND g.charttime<=c.t0 AND g.ph IS NOT NULL),
mv AS (SELECT DISTINCT c.stay_id FROM c JOIN mimiciv_derived.ventilation v ON v.stay_id=c.stay_id
        AND v.ventilation_status IN ('InvasiveVent','Tracheostomy')
        AND v.starttime<=c.t0 AND v.endtime> c.t0 - interval '6 hours'),
htn AS (SELECT DISTINCT c.stay_id FROM c JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=c.hadm_id
         WHERE (d.icd_version=9 AND d.icd_code LIKE '40%') OR (d.icd_version=10 AND d.icd_code LIKE 'I1%'))
SELECT c.stay_id,
       (b.early_crrt) early_crrt,
       vs.hr, vs.rr, vs.spo2, vs.tc,
       (SELECT bun         FROM chem WHERE chem.stay_id=c.stay_id AND rn=1) bun,
       (SELECT potassium   FROM chem WHERE chem.stay_id=c.stay_id AND rn=1) potassium,
       (SELECT bicarbonate FROM chem WHERE chem.stay_id=c.stay_id AND rn=1) bicarbonate,
       (SELECT wbc         FROM cbc  WHERE cbc.stay_id=c.stay_id  AND rn=1) wbc,
       (SELECT hemoglobin  FROM cbc  WHERE cbc.stay_id=c.stay_id  AND rn=1) hemoglobin,
       (SELECT platelet    FROM cbc  WHERE cbc.stay_id=c.stay_id  AND rn=1) platelet,
       (SELECT inr         FROM cg   WHERE cg.stay_id=c.stay_id   AND rn=1) inr,
       (SELECT ph          FROM gas  WHERE gas.stay_id=c.stay_id  AND rn=1) ph,
       (CASE WHEN mv.stay_id  IS NOT NULL THEN 1 ELSE 0 END) mv_flag,
       (CASE WHEN htn.stay_id IS NOT NULL THEN 1 ELSE 0 END) htn_flag
FROM c
JOIN data_extract_crrt_strictk3_v2.baseline b USING (stay_id)
LEFT JOIN vs  USING (stay_id)
LEFT JOIN mv  USING (stay_id)
LEFT JOIN htn USING (stay_id)
")
be <- be %>% mutate(early_crrt = early_crrt %in% TRUE)

be_cont <- function(label, var, digits) tibble(
  Characteristic = label,
  `Overall (best-effort)` = fmt_cont(be[[var]], digits),
  `Early (best-effort)`   = fmt_cont(be[[var]][be$early_crrt], digits),
  `No CRRT (best-effort)` = fmt_cont(be[[var]][!be$early_crrt], digits),
  `nonmissing n`          = sum(is.finite(be[[var]])))
be_bin <- function(label, var) tibble(
  Characteristic = label,
  `Overall (best-effort)` = fmt_cat(be[[var]]==1),
  `Early (best-effort)`   = fmt_cat(be[[var]][be$early_crrt]==1),
  `No CRRT (best-effort)` = fmt_cat(be[[var]][!be$early_crrt]==1),
  `nonmissing n`          = sum(!is.na(be[[var]])))

# Only the rows still NOT promoted to the main table (their standard-derived
# re-derivation does NOT reproduce the rich overall anchor's non-missing count).
audit <- bind_rows(
  be_cont("Heart rate, beats/min", "hr", 1),
  be_cont("Respiratory rate, breaths/min", "rr", 1),
  be_cont("Temperature, \u00b0C", "tc", 1),
  be_cont("SpO2, %", "spo2", 0),
  be_cont("Arterial pH (any specimen)", "ph", 2),
  be_bin ("Mechanical ventilation, n (%)", "mv_flag"),
  be_bin ("Hypertension, n (%)", "htn_flag")
)
write_csv(audit, audit_out)

# ---------------------------------------------------------------------------
# 5. Self-check vs anchor CSVs
# ---------------------------------------------------------------------------
get_cell <- function(tb, char, col) { r <- tb$Characteristic == char; if (!any(r)) return(NA_character_); tb[[col]][which(r)[1]] }
chk <- function(desc, got, want) cat(sprintf("  [%s] %-42s got=%-22s want=%-22s\n",
                                             ifelse(identical(trimws(got), trimws(want)), "PASS", "FAIL"), desc, got, want))

cat("\n================ SELF-CHECK: Overall vs rich overall_sofa8_14 ================\n")
chk("Age",            get_cell(table1_export,"Age, years","Overall (N=1502)"),                 "68.0 [58.0, 78.0]")
chk("SOFA admission", get_cell(table1_export,"SOFA at ICU admission","Overall (N=1502)"),       "10 [8, 12]")
chk("SOFA t0",        get_cell(table1_export,"SOFA at t0","Overall (N=1502)"),                  "10 [9, 12]")
chk("APACHE III",     get_cell(table1_export,"APACHE III (APS III)","Overall (N=1502)"),        "70 [56, 84]")
chk("GCS minimum",    get_cell(table1_export,"GCS (minimum)","Overall (N=1502)"),               "15 [13, 15]")
chk("Charlson",       get_cell(table1_export,"Charlson comorbidity index","Overall (N=1502)"),  "6 [4, 8]")
chk("Creatinine t0",  get_cell(table1_export,"Serum creatinine (t0), mg/dL","Overall (N=1502)"),"2.10 [1.40, 2.90]")
chk("Baseline SCr",   get_cell(table1_export,"Baseline SCr, mg/dL","Overall (N=1502)"),         "0.70 [0.50, 1.00]")
chk("Lactate",        get_cell(table1_export,"Lactate, mmol/L","Overall (N=1502)"),             "2.00 [1.40, 3.00]")
chk("Urine output",   get_cell(table1_export,"Urine output (last 6 h), mL/kg/h","Overall (N=1502)"),"0.21 [0.11, 0.33]")
chk("NEE",            get_cell(table1_export,"Norepinephrine-equivalent dose, \u00b5g/kg/min","Overall (N=1502)"),"0.12 [0.05, 0.25]")
chk("Diabetes",       get_cell(table1_export,"Diabetes mellitus, n (%)","Overall (N=1502)"),    "518 (34.5)")
chk("CHF",            get_cell(table1_export,"Congestive heart failure, n (%)","Overall (N=1502)"),"656 (43.7)")
chk("CKD",            get_cell(table1_export,"Chronic kidney disease, n (%)","Overall (N=1502)"),"426 (28.4)")
chk("Liver disease",  get_cell(table1_export,"Liver disease, n (%)","Overall (N=1502)"),        "409 (27.2)")
chk("Malignancy",     get_cell(table1_export,"Malignancy, n (%)","Overall (N=1502)"),           "274 (18.2)")
chk("COPD",           get_cell(table1_export,"COPD, n (%)","Overall (N=1502)"),                 "446 (29.7)")

cat("\n---- promoted t0 labs (last value in (t0-6h, t0]) vs rich overall anchor ----\n")
chk("BUN",            get_cell(table1_export,"BUN, mg/dL","Overall (N=1502)"),                   "41 [28, 56]")
chk("Bicarbonate",    get_cell(table1_export,"Bicarbonate, mmol/L","Overall (N=1502)"),          "20 [18, 23]")
chk("Potassium",      get_cell(table1_export,"Potassium, mmol/L","Overall (N=1502)"),            "4.40 [4.00, 4.90]")
chk("WBC",            get_cell(table1_export,"White blood cell count, 10^9/L","Overall (N=1502)"),"14.2 [9.9, 20.6]")
chk("Hemoglobin",     get_cell(table1_export,"Hemoglobin, g/dL","Overall (N=1502)"),             "9.4 [8.2, 10.8]")
chk("Platelets",      get_cell(table1_export,"Platelets, 10^9/L","Overall (N=1502)"),            "128 [79, 205]")
chk("INR",            get_cell(table1_export,"INR","Overall (N=1502)"),                          "1.60 [1.30, 2.10]")

cat("\n================ SELF-CHECK: by-arm shared rows vs slim by-arm CSV ============\n")
chk("Age early",      get_cell(table1_export,"Age, years","Early CRRT (n=146)"),                "66.5 [54.0, 74.0]")
chk("Age no",         get_cell(table1_export,"Age, years","No CRRT (n=1356)"),                  "68.0 [58.0, 78.2]")
chk("SOFA t0 early",  get_cell(table1_export,"SOFA at t0","Early CRRT (n=146)"),                "12 [10, 13]")
chk("SOFA t0 no",     get_cell(table1_export,"SOFA at t0","No CRRT (n=1356)"),                  "10 [9, 12]")
chk("Creat t0 early", get_cell(table1_export,"Serum creatinine (t0), mg/dL","Early CRRT (n=146)"),"3.10 [2.32, 4.00]")
chk("Creat t0 no",    get_cell(table1_export,"Serum creatinine (t0), mg/dL","No CRRT (n=1356)"),"2.00 [1.30, 2.80]")
chk("NEE early",      get_cell(table1_export,"Norepinephrine-equivalent dose, \u00b5g/kg/min","Early CRRT (n=146)"),"0.15 [0.05, 0.30]")
chk("Lactate early",  get_cell(table1_export,"Lactate, mmol/L","Early CRRT (n=146)"),           "2.50 [1.60, 4.05]")

cat("\nNEE sanity (must NOT be ~4.8): early NEE =",
    fmt_cont(early$nee_at_t0_corrected, 2), "\n")

cat("\nRows written:", nrow(table1_export),
    "| verified data rows:", sum(!table1_export$`Overall (N=1502)` %in% c("", NA_CELL)),
    "| missing rows:", sum(table1_export$`Overall (N=1502)` == NA_CELL), "\n")
cat("MAP note: main table uses baseline map_t0 =",
    get_cell(table1_export,"MAP, mmHg","Overall (N=1502)"),
    "(matches by-arm slim; rich overall CSV shows 70.0 [64.0, 78.0] from a different MAP source).\n")

cat("\nWrote:\n -", csv_out, "\n -", audit_out, "\n")
