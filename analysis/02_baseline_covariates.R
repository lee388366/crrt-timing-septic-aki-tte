# =============================================================================
# Step 2: Baseline Covariates at t₀
# =============================================================================

library(DBI); library(RPostgres); library(tidyverse); library(lubridate)

con <- dbConnect(RPostgres::Postgres(),
                 dbname="mimiciv31", host="localhost",
                 port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))

cohort_raw <- readRDS("analysis/data/cohort_base.rds")

# Keep first ICU stay per patient
cohort <- cohort_raw %>%
  group_by(subject_id) %>%
  slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("First ICU stays:", nrow(cohort), "\n")

# -----------------------------------------------------------------------------
# 2a. Absolute SOFA score at/near t₀ (closest hourly record before t₀ + 1h)
# -----------------------------------------------------------------------------
sofa_hourly <- dbGetQuery(con, "
  SELECT stay_id, starttime, sofa_24hours AS sofa_total
  FROM mimiciv_derived.sofa
")

sofa_at_t0 <- sofa_hourly %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(starttime <= t0 + hours(1)) %>%          # allow 1h slack
  group_by(stay_id) %>%
  slice_max(starttime, n = 1) %>%                 # closest record before t₀
  ungroup() %>%
  select(stay_id, sofa_t0 = sofa_total)

cat("SOFA available at t0:", nrow(sofa_at_t0), "\n")
cat("SOFA summary:\n"); print(summary(sofa_at_t0$sofa_t0))

# -----------------------------------------------------------------------------
# 2b. APACHE III score (first day)
# -----------------------------------------------------------------------------
apsiii <- dbGetQuery(con, "SELECT stay_id, apsiii FROM mimiciv_derived.apsiii")

# -----------------------------------------------------------------------------
# 2c. Serum creatinine nearest to t₀
# -----------------------------------------------------------------------------
creat_ts <- dbGetQuery(con, "
  SELECT c.hadm_id, c.charttime, c.creatinine
  FROM mimiciv_derived.chemistry c
  WHERE c.creatinine IS NOT NULL
")

# Join via hadm_id (chemistry keyed on hadm_id, not stay_id)
creat_at_t0 <- creat_ts %>%
  inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
  filter(charttime <= t0 + hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1) %>%
  ungroup() %>%
  select(stay_id, creatinine_t0 = creatinine)

cat("Creatinine at t0:", nrow(creat_at_t0), "\n")

# -----------------------------------------------------------------------------
# 2d. Vasopressor NEE dose at t₀ (already have from cohort build)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2e. Lactate nearest to t₀
# -----------------------------------------------------------------------------
lactate_ts <- dbGetQuery(con, "
  SELECT hadm_id, charttime, lactate
  FROM mimiciv_derived.bg
  WHERE lactate IS NOT NULL
")

lactate_at_t0 <- lactate_ts %>%
  inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
  filter(charttime <= t0 + hours(4), charttime >= t0 - hours(12)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1) %>%
  ungroup() %>%
  select(stay_id, lactate_t0 = lactate)

cat("Lactate at t0:", nrow(lactate_at_t0), "\n")

# -----------------------------------------------------------------------------
# 2f. Urine output rate nearest t₀ (ml/kg/h)
# -----------------------------------------------------------------------------
uo_ts <- dbGetQuery(con, "
  SELECT stay_id, charttime, uo_mlkghr_6hr AS urineoutput_rate
  FROM mimiciv_derived.urine_output_rate
  WHERE uo_mlkghr_6hr IS NOT NULL
")

uo_at_t0 <- uo_ts %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(charttime <= t0 + hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1) %>%
  ungroup() %>%
  select(stay_id, uo_rate_t0 = urineoutput_rate)

# -----------------------------------------------------------------------------
# 2g. Charlson Comorbidity Index (baseline severity)
# -----------------------------------------------------------------------------
charlson <- dbGetQuery(con, "SELECT hadm_id, charlson_comorbidity_index AS cci FROM mimiciv_derived.charlson")

# -----------------------------------------------------------------------------
# 2h. MAP at t₀ from vitalsign
# -----------------------------------------------------------------------------
vitals_ts <- dbGetQuery(con, "
  SELECT stay_id, charttime, mbp AS map
  FROM mimiciv_derived.vitalsign
  WHERE mbp IS NOT NULL
")

map_at_t0 <- vitals_ts %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(charttime <= t0 + hours(1), charttime >= t0 - hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1) %>%
  ungroup() %>%
  select(stay_id, map_t0 = map)

# -----------------------------------------------------------------------------
# 3. ASSEMBLE baseline dataset
# -----------------------------------------------------------------------------
baseline <- cohort %>%
  left_join(sofa_at_t0,    by = "stay_id") %>%
  left_join(apsiii,        by = "stay_id") %>%
  left_join(creat_at_t0,   by = "stay_id") %>%
  left_join(lactate_at_t0, by = "stay_id") %>%
  left_join(uo_at_t0,      by = "stay_id") %>%
  left_join(charlson,      by = "hadm_id") %>%
  left_join(map_at_t0,     by = "stay_id") %>%
  mutate(
    # Equipoise zone: SOFA 8-14
    sofa_equipoise  = sofa_t0 >= 8 & sofa_t0 <= 14,
    sofa_low        = sofa_t0 <= 7,
    sofa_high       = sofa_t0 >= 15,

    # Creatinine ≥ 3× baseline OR ≥ 4.0 with ≥0.5 rise (KDIGO3 creat check)
    creat_3x_base   = !is.na(creatinine_t0) & !is.na(creatinine_baseline) &
                       creatinine_t0 >= 3 * creatinine_baseline,
    creat_4mg_rise  = !is.na(creatinine_t0) &
                       creatinine_t0 >= 4.0 &
                       (creatinine_t0 - coalesce(creatinine_baseline, 0)) >= 0.5,

    # Outcome: CRRT-free days alive at day 28
    crrt_free_days  = if_else(died_28d, 0,
                        pmax(0, 28 - coalesce(crrt_hours_from_t0/24, 0)))
  )

# -----------------------------------------------------------------------------
# 4. Summary table: full cohort + equipoise subgroup
# -----------------------------------------------------------------------------
summarise_grp <- function(df, label) {
  tibble(
    Group          = label,
    N              = nrow(df),
    Age_med        = median(df$age, na.rm=T),
    Female_pct     = round(100*mean(df$gender=="F", na.rm=T),1),
    SOFA_med       = median(df$sofa_t0, na.rm=T),
    APSIII_med     = median(df$apsiii, na.rm=T),
    Creat_t0_med   = median(df$creatinine_t0, na.rm=T),
    NEE_t0_med     = median(df$nee_at_t0, na.rm=T),
    Lactate_t0_med = median(df$lactate_t0, na.rm=T),
    Mortality_28d  = round(100*mean(df$died_28d, na.rm=T),1),
    EarlyCRRT_pct  = round(100*mean(df$early_crrt, na.rm=T),1),
    AnyCRRT_pct    = round(100*mean(!is.na(df$crrt_first_time), na.rm=T),1)
  )
}

tbl <- bind_rows(
  summarise_grp(baseline,                            "Full cohort"),
  summarise_grp(filter(baseline, sofa_equipoise),    "SOFA 8-14 (equipoise)"),
  summarise_grp(filter(baseline, sofa_low),          "SOFA ≤7"),
  summarise_grp(filter(baseline, sofa_high),         "SOFA ≥15")
)

cat("\n=== BASELINE SUMMARY TABLE ===\n")
print(t(tbl), quote=FALSE)

# -----------------------------------------------------------------------------
# 5. CRRT initiation rate by SOFA group (to verify equipoise rationale)
# -----------------------------------------------------------------------------
cat("\n=== CRRT INITIATION RATE BY SOFA ===\n")
baseline %>%
  mutate(sofa_grp = case_when(
    sofa_t0 <= 7  ~ "SOFA ≤7",
    sofa_t0 <= 14 ~ "SOFA 8-14",
    TRUE          ~ "SOFA ≥15"
  )) %>%
  group_by(sofa_grp) %>%
  summarise(
    n           = n(),
    crrt_any    = round(100*mean(!is.na(crrt_first_time)),1),
    crrt_early  = round(100*mean(early_crrt),1),
    mort_28d    = round(100*mean(died_28d),1),
    .groups = "drop"
  ) %>%
  print()

# -----------------------------------------------------------------------------
# 6. Save
# -----------------------------------------------------------------------------
saveRDS(baseline, "analysis/data/baseline.rds")
dbWriteTable(con, Id(schema="data_extract_crrt", table="baseline"),
             baseline, overwrite=TRUE)

cat("\n✅ Saved: analysis/data/baseline.rds\n")
dbDisconnect(con)
