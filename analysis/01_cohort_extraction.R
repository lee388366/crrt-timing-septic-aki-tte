# =============================================================================
# Step 1: Cohort Extraction
# Study: Early vs Deferred CRRT in Septic Shock + KDIGO Stage 3 AKI
# Data:  MIMIC-IV v3.1 (local PostgreSQL: mimiciv31)
# =============================================================================

library(DBI)
library(RPostgres)
library(tidyverse)
library(lubridate)

con <- dbConnect(RPostgres::Postgres(),
                 dbname = "mimiciv31", host = "localhost",
                 port = 5432, user = "postgres", password = Sys.getenv("PGPASSWORD"))

# -----------------------------------------------------------------------------
# 1. ICU stay details (age, admission/discharge times, mortality)
# -----------------------------------------------------------------------------
icu_detail <- dbGetQuery(con, "
  SELECT
    icu.subject_id,
    icu.hadm_id,
    icu.stay_id,
    icu.intime                              AS icu_intime,
    icu.outtime                             AS icu_outtime,
    EXTRACT(EPOCH FROM (icu.outtime - icu.intime))/3600 AS icu_los_hours,
    pat.anchor_age                          AS age,
    pat.gender,
    adm.hospital_expire_flag,
    adm.admittime,
    adm.dischtime,
    adm.deathtime
  FROM mimiciv_icu.icustays icu
  JOIN mimiciv_hosp.patients pat USING (subject_id)
  JOIN mimiciv_hosp.admissions adm USING (hadm_id)
  WHERE pat.anchor_age >= 18
")

cat("ICU stays (adults):", nrow(icu_detail), "\n")

# -----------------------------------------------------------------------------
# 2. KDIGO Stage 3 AKI — first occurrence per stay
# -----------------------------------------------------------------------------
kdigo_stage3 <- dbGetQuery(con, "
  SELECT
    stay_id,
    MIN(charttime) AS kdigo3_time,
    -- which criterion triggered stage 3
    BOOL_OR(aki_stage_creat = 3) AS stage3_by_creat,
    BOOL_OR(aki_stage_uo    = 3) AS stage3_by_uo
  FROM mimiciv_derived.kdigo_stages
  WHERE aki_stage = 3
  GROUP BY stay_id
")

cat("Stays with KDIGO Stage 3:", nrow(kdigo_stage3), "\n")

# -----------------------------------------------------------------------------
# 3a. ESRD / chronic dialysis exclusion list (by hadm_id)
#     ICD-10: N186 (ESRD), N185 (CKD stage 5)
#     ICD-9:  5856 (ESRD), V4511/V4512 (dialysis dependence)
# -----------------------------------------------------------------------------
esrd_hadm <- dbGetQuery(con, "
  SELECT DISTINCT hadm_id
  FROM mimiciv_hosp.diagnoses_icd
  WHERE (icd_version = 10 AND icd_code IN ('N186','N185'))
     OR (icd_version = 9  AND icd_code IN ('5856','V4511','V4512'))
")
cat("ESRD / chronic dialysis admissions:", nrow(esrd_hadm), "\n")

# -----------------------------------------------------------------------------
# 3b. Sepsis-3 — must be confirmed AT OR BEFORE t₀
#     (t₀ = kdigo3_time − 6h, so sofa_time must be ≤ t₀)
#     We join after computing t₀ in Step 6
# -----------------------------------------------------------------------------
sepsis3_stays <- dbGetQuery(con, "
  SELECT stay_id, sofa_time, sofa_score
  FROM mimiciv_derived.sepsis3
  WHERE sepsis3 = TRUE
")

# Full vasopressor time-series — compute NEE from individual agents
# NEE formula (weights per MIMIC convention): NE*1 + Epi*1 + Dopa/100 + Phenyl*0.1 + Vaso*2.5
vasoactive <- dbGetQuery(con, "
  SELECT stay_id, starttime, endtime,
         COALESCE(norepinephrine, 0)
           + COALESCE(epinephrine, 0)
           + COALESCE(dopamine, 0) / 100.0
           + COALESCE(phenylephrine, 0) * 0.1
           + COALESCE(vasopressin, 0) * 2.5   AS nee_dose
  FROM mimiciv_derived.vasoactive_agent
")

cat("Sepsis-3 stays:", length(unique(sepsis3_stays$stay_id)), "\n")
cat("Vasopressor records:", nrow(vasoactive), "\n")

# -----------------------------------------------------------------------------
# 4. CRRT timing per stay
#    Use rrt table: dialysis_active=1 AND dialysis_type='CRRT'
# -----------------------------------------------------------------------------
crrt_times <- dbGetQuery(con, "
  SELECT
    stay_id,
    MIN(charttime) AS crrt_first_time
  FROM mimiciv_derived.rrt
  WHERE dialysis_active = 1
    AND dialysis_type = 'CRRT'
  GROUP BY stay_id
")

cat("Stays with CRRT:", nrow(crrt_times), "\n")

# -----------------------------------------------------------------------------
# 5. Baseline creatinine (pre-AKI) — keyed on hadm_id
# -----------------------------------------------------------------------------
creat_baseline <- dbGetQuery(con, "
  SELECT hadm_id, scr_baseline AS creatinine_baseline
  FROM mimiciv_derived.creatinine_baseline
")

# -----------------------------------------------------------------------------
# 6. MERGE and compute t₀, then apply precise time-window checks in R
# -----------------------------------------------------------------------------

# Step A: base merge (KDIGO3 + ICU details + baseline creatinine + CRRT)
base <- icu_detail %>%
  inner_join(kdigo_stage3,   by = "stay_id") %>%
  left_join(creat_baseline,  by = "hadm_id") %>%
  left_join(crrt_times,      by = "stay_id") %>%
  mutate(
    t0 = pmax(icu_intime, kdigo3_time - hours(6))
  )

cat("\nKDIGO3 + ICU detail:", nrow(base), "\n")

# Step B: Sepsis-3 must be confirmed AT OR BEFORE t₀ (methods: "at a prespecified index time point t₀")
#         t₀ is already computed in base; sofa_time must be ≤ t₀
sepsis3_valid <- sepsis3_stays %>%
  inner_join(base %>% select(stay_id, t0), by = "stay_id") %>%
  filter(sofa_time <= t0) %>%        # ← fixed: use t0, not kdigo3_time
  group_by(stay_id) %>%
  slice_min(sofa_time, n = 1) %>%
  ungroup() %>%
  select(stay_id, sofa_time, sofa_score)

cat("Stays with Sepsis-3 at or before t0:", nrow(sepsis3_valid), "\n")

base <- base %>% inner_join(sepsis3_valid, by = "stay_id")
cat("After sepsis3 filter:", nrow(base), "\n")

# Step C: Vasopressor active during [t0-6h, t0] = [kdigo3-12h, kdigo3-6h]
#         i.e., vaso starttime <= t0 AND endtime >= t0 - 6h
vaso_window <- vasoactive %>%
  inner_join(base %>% select(stay_id, t0), by = "stay_id") %>%
  filter(starttime <= t0,
         endtime   >= t0 - hours(6)) %>%
  group_by(stay_id) %>%
  summarise(
    vaso_in_window = TRUE,
    nee_at_t0 = mean(nee_dose, na.rm = TRUE),  # mean NEE in window
    .groups = "drop"
  )

cat("Stays with vasopressor active at t0:", nrow(vaso_window), "\n")

cohort <- base %>%
  inner_join(vaso_window, by = "stay_id") %>%
  mutate(
    # Hours from t₀ to CRRT initiation (NA if no CRRT)
    crrt_hours_from_t0 = as.numeric(
      difftime(crrt_first_time, t0, units = "hours")),

    # Treatment assignment: CRRT within 24h of t₀
    early_crrt = !is.na(crrt_hours_from_t0) & crrt_hours_from_t0 <= 24 &
                  crrt_hours_from_t0 >= 0,

    # CRRT before t₀ (exclusion criterion)
    crrt_before_t0 = !is.na(crrt_first_time) & crrt_first_time < t0,

    # 28-day mortality from t₀
    days_to_death = as.numeric(difftime(deathtime, t0, units = "days")),
    died_28d      = !is.na(days_to_death) & days_to_death <= 28,

    # Follow-up time (capped at 28 days)
    follow_time_28d = pmin(
      28,
      coalesce(
        as.numeric(difftime(deathtime, t0, units = "days")),
        as.numeric(difftime(dischtime, t0, units = "days")),
        28
      )
    ),
    follow_time_28d = pmax(follow_time_28d, 0),

    # ICU LOS from t₀ (days)
    icu_los_from_t0 = as.numeric(difftime(icu_outtime, t0, units = "days"))
  )

cat("\nAfter merging (septic shock + KDIGO3 + vasopressor at t0):", nrow(cohort), "\n")

# -----------------------------------------------------------------------------
# 7. Apply exclusion criteria
# -----------------------------------------------------------------------------
n0 <- nrow(cohort)
excl <- list()

# E1: CRRT before t₀
cohort <- cohort %>% filter(!crrt_before_t0)
excl[["E1_crrt_before_t0"]] <- n0 - nrow(cohort); n0 <- nrow(cohort)

# E2: ESRD or chronic dialysis dependence (new — per methods exclusion criterion #2)
cohort <- cohort %>% filter(!hadm_id %in% esrd_hadm$hadm_id)
excl[["E2_esrd_chronic_dialysis"]] <- n0 - nrow(cohort); n0 <- nrow(cohort)

# E3: ICU stay < 24h from admission
cohort <- cohort %>% filter(icu_los_hours >= 24)
excl[["E3_icu_los_lt_24h"]] <- n0 - nrow(cohort); n0 <- nrow(cohort)

# E4: Missing key baseline covariates: creatinine_baseline, nee_at_t0, t0
#     (methods: "serum creatinine at t₀, vasopressor dose at t₀, or ICU admission time")
cohort <- cohort %>% filter(!is.na(creatinine_baseline), !is.na(nee_at_t0), !is.na(t0))
excl[["E4_missing_baseline"]] <- n0 - nrow(cohort); n0 <- nrow(cohort)

# E5: follow-up time = 0 (t₀ after discharge/death — data anomaly)
cohort <- cohort %>% filter(follow_time_28d > 0)
excl[["E5_zero_followup"]] <- n0 - nrow(cohort); n0 <- nrow(cohort)

# Keep one row per stay_id (in case of duplicates from multiple SOFA records)
cohort <- cohort %>%
  arrange(stay_id, sofa_time) %>%
  distinct(stay_id, .keep_all = TRUE)
excl[["E6_dedup"]] <- n0 - nrow(cohort)

cat("\n=== EXCLUSION FLOW ===\n")
for(nm in names(excl)) cat(sprintf("  %-30s excluded: %d\n", nm, excl[[nm]]))
cat(sprintf("\nFinal cohort: %d stays\n", nrow(cohort)))

# -----------------------------------------------------------------------------
# 8. Descriptive summary
# -----------------------------------------------------------------------------
cat("\n=== COHORT SUMMARY ===\n")
cat(sprintf("  Age (median [IQR]):        %.0f [%.0f–%.0f]\n",
    median(cohort$age), quantile(cohort$age,.25), quantile(cohort$age,.75)))
cat(sprintf("  Female (%%):               %.1f%%\n",
    100*mean(cohort$gender=="F")))
cat(sprintf("  28-day mortality:         %.1f%%\n",
    100*mean(cohort$died_28d)))
cat(sprintf("  Received CRRT (any):      %.1f%%\n",
    100*mean(!is.na(cohort$crrt_first_time))))
cat(sprintf("  Early CRRT (≤24h of t₀): %.1f%%\n",
    100*mean(cohort$early_crrt)))
cat(sprintf("  KDIGO3 by creatinine:     %.1f%%\n",
    100*mean(cohort$stage3_by_creat)))
cat(sprintf("  KDIGO3 by UO:             %.1f%%\n",
    100*mean(cohort$stage3_by_uo)))

# -----------------------------------------------------------------------------
# 9. Save to PostgreSQL and RDS
# -----------------------------------------------------------------------------
dir.create("analysis/data", 
           showWarnings = FALSE, recursive = TRUE)

saveRDS(cohort, "analysis/data/cohort_base.rds")

# Write to DB schema for subsequent steps
dbExecute(con, "CREATE SCHEMA IF NOT EXISTS data_extract_crrt")
dbWriteTable(con, Id(schema="data_extract_crrt", table="cohort_base"),
             cohort, overwrite=TRUE)

cat("\n✅ Saved: analysis/data/cohort_base.rds\n")
cat("✅ Saved: data_extract_crrt.cohort_base (PostgreSQL)\n")

dbDisconnect(con)
