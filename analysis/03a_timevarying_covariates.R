# =============================================================================
# Step 3a: Extract Time-Varying Covariates for IPCW (hours 0–23 from t₀)
#
# Per methods: denominator model requires L_k at each decision hour k:
#   - Last 1h MAP (vitalsign.mbp)
#   - NEE dose (µg/kg/min) at hour k (vasoactive_agent)
#   - Last 6h lactate (bg.lactate)
#   - SOFA score trajectory (sofa.sofa_24hours)
#   - KDIGO stage trajectory (kdigo_stages.aki_stage)
#   - Cumulative fluid balance from t₀ (inputevents - outputevents)
#   - Creatinine change from t₀ (chemistry.creatinine delta)
#
# Strategy A denominator: covariates at hour 23 (cross-sectional)
# Strategy B denominator: covariates at each hour 0–23 (pooled longitudinal)
# =============================================================================

library(DBI); library(RPostgres)
library(tidyverse); library(lubridate); library(zoo)

con <- dbConnect(RPostgres::Postgres(),
                 dbname="mimiciv31", host="localhost",
                 port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))

cohort <- readRDS("analysis/data/cohort_base.rds") %>%
  group_by(subject_id) %>% slice_min(icu_intime, n=1) %>% ungroup()

cat("Cohort size:", nrow(cohort), "\n")

# Grace period: 24 hours
GRACE <- 24L

# Build hourly time grid for each patient: t₀ + 0h, t₀ + 1h, ..., t₀ + 23h
time_grid <- cohort %>%
  select(stay_id, hadm_id, subject_id, t0) %>%
  crossing(hour_index = 0:(GRACE - 1)) %>%
  mutate(
    window_end   = t0 + hours(hour_index + 1),   # end of this 1h window
    window_start = t0 + hours(hour_index)         # start of this 1h window
  )

cat("Time grid rows:", nrow(time_grid), "\n")

# =============================================================================
# 1. MAP — last value in [window_start, window_end]
# =============================================================================
cat("Extracting MAP...\n")
vitals_raw <- dbGetQuery(con, "
  SELECT stay_id, charttime, mbp AS map
  FROM mimiciv_derived.vitalsign
  WHERE mbp IS NOT NULL
    AND mbp BETWEEN 20 AND 200
")

map_hourly <- time_grid %>%
  select(stay_id, hour_index, window_start, window_end) %>%
  left_join(vitals_raw, by="stay_id", relationship="many-to-many") %>%
  filter(charttime >= window_start, charttime < window_end) %>%
  group_by(stay_id, hour_index) %>%
  summarise(map_h = mean(map, na.rm=TRUE), .groups="drop")

cat("  MAP hourly records:", nrow(map_hourly), "\n")

# =============================================================================
# 2. NEE dose — mean dose in [window_start, window_end]
#    NEE = NE + Epi + Dopa/100 + Phenyl*0.1 + Vaso*2.5
# =============================================================================
cat("Extracting NEE...\n")
vaso_raw <- dbGetQuery(con, "
  SELECT stay_id, starttime, endtime,
         COALESCE(norepinephrine, 0)
           + COALESCE(epinephrine, 0)
           + COALESCE(dopamine, 0) / 100.0
           + COALESCE(phenylephrine, 0) * 0.1
           + COALESCE(vasopressin, 0) * 2.5  AS nee_dose
  FROM mimiciv_derived.vasoactive_agent
")

nee_hourly <- time_grid %>%
  select(stay_id, hour_index, window_start, window_end) %>%
  left_join(vaso_raw, by="stay_id", relationship="many-to-many") %>%
  filter(starttime < window_end, endtime > window_start) %>%
  group_by(stay_id, hour_index) %>%
  summarise(nee_h = mean(nee_dose, na.rm=TRUE), .groups="drop")

cat("  NEE hourly records:", nrow(nee_hourly), "\n")

# =============================================================================
# 3. Lactate — last value in preceding 6h window [window_end-6h, window_end]
# =============================================================================
cat("Extracting lactate...\n")
lactate_raw <- dbGetQuery(con, "
  SELECT hadm_id, charttime, lactate
  FROM mimiciv_derived.bg
  WHERE lactate IS NOT NULL AND lactate < 30
")

lactate_hourly <- time_grid %>%
  select(stay_id, hadm_id, hour_index, window_end) %>%
  left_join(lactate_raw, by="hadm_id", relationship="many-to-many") %>%
  filter(charttime >= window_end - hours(6), charttime < window_end) %>%
  group_by(stay_id, hour_index) %>%
  summarise(lactate_h = mean(lactate, na.rm=TRUE), .groups="drop")

cat("  Lactate hourly records:", nrow(lactate_hourly), "\n")

# =============================================================================
# 4. SOFA score — last value at or before window_end
# =============================================================================
cat("Extracting SOFA...\n")
sofa_raw <- dbGetQuery(con, "
  SELECT stay_id, starttime, sofa_24hours AS sofa_total
  FROM mimiciv_derived.sofa
  WHERE sofa_24hours IS NOT NULL
")

sofa_hourly <- time_grid %>%
  select(stay_id, hour_index, window_end) %>%
  left_join(sofa_raw, by="stay_id", relationship="many-to-many") %>%
  filter(starttime <= window_end) %>%
  group_by(stay_id, hour_index) %>%
  slice_max(starttime, n=1) %>%
  summarise(sofa_h = first(sofa_total), .groups="drop")

cat("  SOFA hourly records:", nrow(sofa_hourly), "\n")

# =============================================================================
# 5. KDIGO stage — at window_end
# =============================================================================
cat("Extracting KDIGO stage...\n")
kdigo_raw <- dbGetQuery(con, "
  SELECT stay_id, charttime, aki_stage
  FROM mimiciv_derived.kdigo_stages
  WHERE aki_stage IS NOT NULL
")

kdigo_hourly <- time_grid %>%
  select(stay_id, hour_index, window_end) %>%
  left_join(kdigo_raw, by="stay_id", relationship="many-to-many") %>%
  filter(charttime <= window_end) %>%
  group_by(stay_id, hour_index) %>%
  slice_max(charttime, n=1) %>%
  summarise(kdigo_h = first(aki_stage), .groups="drop")

cat("  KDIGO hourly records:", nrow(kdigo_hourly), "\n")

# =============================================================================
# 6. Creatinine — last value in [window_end-6h, window_end]
#    + delta from t₀ creatinine
# =============================================================================
cat("Extracting creatinine...\n")
creat_raw <- dbGetQuery(con, "
  SELECT hadm_id, charttime, creatinine
  FROM mimiciv_derived.chemistry
  WHERE creatinine IS NOT NULL AND creatinine < 25
")

# Creatinine at t₀ (for delta calculation)
creat_t0 <- creat_raw %>%
  inner_join(cohort %>% select(stay_id, hadm_id, t0), by="hadm_id") %>%
  filter(charttime <= t0 + hours(2), charttime >= t0 - hours(6)) %>%
  group_by(stay_id) %>% slice_max(charttime, n=1) %>% ungroup() %>%
  select(stay_id, creatinine_t0 = creatinine)

creat_hourly <- time_grid %>%
  select(stay_id, hadm_id, hour_index, window_end) %>%
  left_join(creat_raw, by="hadm_id", relationship="many-to-many") %>%
  filter(charttime >= window_end - hours(6), charttime < window_end) %>%
  group_by(stay_id, hour_index) %>%
  summarise(creatinine_h = mean(creatinine, na.rm=TRUE), .groups="drop") %>%
  left_join(creat_t0, by="stay_id") %>%
  mutate(creatinine_delta = creatinine_h - coalesce(creatinine_t0, creatinine_h))

cat("  Creatinine hourly records:", nrow(creat_hourly), "\n")

# =============================================================================
# 7. Fluid balance — cumulative net from t₀ to window_end
#    net = total_input - total_output in [t₀, window_end]
# =============================================================================
cat("Extracting fluid balance...\n")
input_raw <- dbGetQuery(con, "
  SELECT stay_id, starttime, amount
  FROM mimiciv_icu.inputevents
  WHERE amount > 0 AND amount < 50000
    AND amountuom IN ('mL','ml')
")

output_raw <- dbGetQuery(con, "
  SELECT stay_id, charttime, value AS amount
  FROM mimiciv_icu.outputevents
  WHERE value > 0 AND value < 20000
")

# Cumulative input up to window_end from t₀
input_cumul <- time_grid %>%
  select(stay_id, hour_index, t0, window_end) %>%
  left_join(input_raw, by="stay_id", relationship="many-to-many") %>%
  filter(starttime >= t0, starttime < window_end) %>%
  group_by(stay_id, hour_index) %>%
  summarise(total_input = sum(amount, na.rm=TRUE), .groups="drop")

output_cumul <- time_grid %>%
  select(stay_id, hour_index, t0, window_end) %>%
  left_join(output_raw, by="stay_id", relationship="many-to-many") %>%
  filter(charttime >= t0, charttime < window_end) %>%
  group_by(stay_id, hour_index) %>%
  summarise(total_output = sum(amount, na.rm=TRUE), .groups="drop")

fluid_hourly <- full_join(input_cumul, output_cumul, by=c("stay_id","hour_index")) %>%
  mutate(fluid_balance_h = coalesce(total_input, 0) - coalesce(total_output, 0))

cat("  Fluid balance hourly records:", nrow(fluid_hourly), "\n")

# =============================================================================
# 8. ASSEMBLE hourly time-varying dataset with LOCF imputation
# =============================================================================
cat("Assembling hourly dataset...\n")
tv_hourly <- time_grid %>%
  select(stay_id, hour_index) %>%
  left_join(map_hourly,      by=c("stay_id","hour_index")) %>%
  left_join(nee_hourly,      by=c("stay_id","hour_index")) %>%
  left_join(lactate_hourly,  by=c("stay_id","hour_index")) %>%
  left_join(sofa_hourly,     by=c("stay_id","hour_index")) %>%
  left_join(kdigo_hourly,    by=c("stay_id","hour_index")) %>%
  left_join(creat_hourly %>% select(stay_id, hour_index, creatinine_h, creatinine_delta),
            by=c("stay_id","hour_index")) %>%
  left_join(fluid_hourly %>% select(stay_id, hour_index, fluid_balance_h),
            by=c("stay_id","hour_index")) %>%
  # LOCF within each patient (per methods: "carried forward within each clone")
  group_by(stay_id) %>%
  mutate(across(c(map_h, nee_h, lactate_h, sofa_h, kdigo_h,
                  creatinine_h, creatinine_delta, fluid_balance_h),
                ~zoo::na.locf(., na.rm=FALSE))) %>%
  # Fill remaining NAs (first hours with no prior obs) with 0 for doses, median for others
  mutate(
    nee_h           = coalesce(nee_h, 0),
    fluid_balance_h = coalesce(fluid_balance_h, 0)
  ) %>%
  ungroup()

cat("Total hourly rows:", nrow(tv_hourly), "\n")
cat("Missing rates:\n")
cat(sprintf("  MAP:             %.1f%%\n", 100*mean(is.na(tv_hourly$map_h))))
cat(sprintf("  NEE:             %.1f%%\n", 100*mean(is.na(tv_hourly$nee_h))))
cat(sprintf("  Lactate:         %.1f%%\n", 100*mean(is.na(tv_hourly$lactate_h))))
cat(sprintf("  SOFA:            %.1f%%\n", 100*mean(is.na(tv_hourly$sofa_h))))
cat(sprintf("  Creatinine:      %.1f%%\n", 100*mean(is.na(tv_hourly$creatinine_h))))
cat(sprintf("  Fluid balance:   %.1f%%\n", 100*mean(is.na(tv_hourly$fluid_balance_h))))

# =============================================================================
# 9. Covariates at hour 23 (for Strategy A denominator)
# =============================================================================
tv_hour23 <- tv_hourly %>% filter(hour_index == 23)
cat("\nHour-23 cross-section: n =", nrow(tv_hour23), "\n")

# =============================================================================
# 10. Save
# =============================================================================
saveRDS(tv_hourly,  "analysis/data/tv_hourly.rds")
saveRDS(tv_hour23,  "analysis/data/tv_hour23.rds")
saveRDS(creat_t0,   "analysis/data/creatinine_t0.rds")

cat("\n✅ Saved tv_hourly, tv_hour23, creatinine_t0\n")
dbDisconnect(con)
