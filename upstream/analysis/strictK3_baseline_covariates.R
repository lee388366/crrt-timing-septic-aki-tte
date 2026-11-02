# Strict first-KDIGO3 baseline covariates at t0.
# Reuses the legacy 02_baseline_covariates.R definitions but reads/writes the
# strictK3_primary data directory so legacy analysis/data outputs stay intact.
suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(tidyverse)
  library(lubridate)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
data_dir <- file.path(root_dir, "strictK3_primary", "data")

con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimiciv31",
  host = "localhost",
  port = 5432,
  user = "postgres",
  password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

cohort_raw <- readRDS(file.path(data_dir, "cohort_base.rds"))

cohort <- cohort_raw %>%
  group_by(subject_id) %>%
  slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("Strict K3 first ICU stays:", nrow(cohort), "\n")

sofa_hourly <- dbGetQuery(con, "
  SELECT stay_id, starttime, sofa_24hours AS sofa_total
  FROM mimiciv_derived.sofa
")

sofa_at_t0 <- sofa_hourly %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(starttime <= t0 + hours(1)) %>%
  group_by(stay_id) %>%
  slice_max(starttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stay_id, sofa_t0 = sofa_total)

cat("SOFA available at t0:", nrow(sofa_at_t0), "\n")

apsiii <- dbGetQuery(con, "SELECT stay_id, apsiii FROM mimiciv_derived.apsiii")

creat_ts <- dbGetQuery(con, "
  SELECT hadm_id, charttime, creatinine
  FROM mimiciv_derived.chemistry
  WHERE creatinine IS NOT NULL
")

creat_at_t0 <- creat_ts %>%
  inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
  filter(charttime <= t0 + hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stay_id, creatinine_t0 = creatinine)

cat("Creatinine at t0:", nrow(creat_at_t0), "\n")

lactate_ts <- dbGetQuery(con, "
  SELECT hadm_id, charttime, lactate
  FROM mimiciv_derived.bg
  WHERE lactate IS NOT NULL
")

lactate_at_t0 <- lactate_ts %>%
  inner_join(cohort %>% select(stay_id, hadm_id, t0), by = "hadm_id") %>%
  filter(charttime <= t0 + hours(4), charttime >= t0 - hours(12)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stay_id, lactate_t0 = lactate)

cat("Lactate at t0:", nrow(lactate_at_t0), "\n")

uo_ts <- dbGetQuery(con, "
  SELECT stay_id, charttime, uo_mlkghr_6hr AS urineoutput_rate
  FROM mimiciv_derived.urine_output_rate
  WHERE uo_mlkghr_6hr IS NOT NULL
")

uo_at_t0 <- uo_ts %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(charttime <= t0 + hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stay_id, uo_rate_t0 = urineoutput_rate)

charlson <- dbGetQuery(con, "
  SELECT hadm_id, charlson_comorbidity_index AS cci
  FROM mimiciv_derived.charlson
")

vitals_ts <- dbGetQuery(con, "
  SELECT stay_id, charttime, mbp AS map
  FROM mimiciv_derived.vitalsign
  WHERE mbp IS NOT NULL
")

map_at_t0 <- vitals_ts %>%
  inner_join(cohort %>% select(stay_id, t0), by = "stay_id") %>%
  filter(charttime <= t0 + hours(1), charttime >= t0 - hours(2)) %>%
  group_by(stay_id) %>%
  slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(stay_id, map_t0 = map)

baseline <- cohort %>%
  left_join(sofa_at_t0, by = "stay_id") %>%
  left_join(apsiii, by = "stay_id") %>%
  left_join(creat_at_t0, by = "stay_id") %>%
  left_join(lactate_at_t0, by = "stay_id") %>%
  left_join(uo_at_t0, by = "stay_id") %>%
  left_join(charlson, by = "hadm_id") %>%
  left_join(map_at_t0, by = "stay_id") %>%
  mutate(
    sofa_equipoise = sofa_t0 >= 8 & sofa_t0 <= 14,
    sofa_low = sofa_t0 <= 7,
    sofa_high = sofa_t0 >= 15,
    sofa_band_prespecified = case_when(
      sofa_t0 < 9 ~ "SOFA_lt9",
      sofa_t0 >= 9 & sofa_t0 <= 13 ~ "SOFA_9_13",
      sofa_t0 > 13 ~ "SOFA_gt13",
      TRUE ~ NA_character_
    ),
    creat_3x_base = !is.na(creatinine_t0) & !is.na(creatinine_baseline) &
      creatinine_t0 >= 3 * creatinine_baseline,
    creat_4mg_rise = !is.na(creatinine_t0) &
      creatinine_t0 >= 4.0 &
      (creatinine_t0 - coalesce(creatinine_baseline, 0)) >= 0.5,
    crrt_free_days = if_else(
      died_28d,
      0,
      pmax(0, 28 - coalesce(crrt_hours_from_t0 / 24, 0))
    )
  )

summarise_grp <- function(df, label) {
  tibble(
    Group = label,
    N = nrow(df),
    Age_med = median(df$age, na.rm = TRUE),
    Female_pct = round(100 * mean(df$gender == "F", na.rm = TRUE), 1),
    SOFA_med = median(df$sofa_t0, na.rm = TRUE),
    APSIII_med = median(df$apsiii, na.rm = TRUE),
    Creat_t0_med = median(df$creatinine_t0, na.rm = TRUE),
    NEE_t0_med = median(df$nee_at_t0, na.rm = TRUE),
    Lactate_t0_med = median(df$lactate_t0, na.rm = TRUE),
    Mortality_28d = round(100 * mean(df$died_28d, na.rm = TRUE), 1),
    EarlyCRRT_pct = round(100 * mean(df$early_crrt, na.rm = TRUE), 1),
    AnyCRRT_pct = round(100 * mean(!is.na(df$crrt_first_time), na.rm = TRUE), 1)
  )
}

baseline_summary <- bind_rows(
  summarise_grp(baseline, "Full strictK3 primary cohort"),
  summarise_grp(filter(baseline, sofa_equipoise), "SOFA 8-14 equipoise"),
  summarise_grp(filter(baseline, sofa_low), "SOFA <=7"),
  summarise_grp(filter(baseline, sofa_high), "SOFA >=15")
)

sofa_group_summary <- baseline %>%
  mutate(sofa_grp = case_when(
    sofa_t0 <= 7 ~ "SOFA <=7",
    sofa_t0 <= 14 ~ "SOFA 8-14",
    TRUE ~ "SOFA >=15"
  )) %>%
  group_by(sofa_grp) %>%
  summarise(
    n = n(),
    crrt_any = round(100 * mean(!is.na(crrt_first_time)), 1),
    crrt_early = round(100 * mean(early_crrt), 1),
    mort_28d = round(100 * mean(died_28d), 1),
    .groups = "drop"
  )

saveRDS(baseline, file.path(data_dir, "baseline.rds"))
write_csv(baseline, file.path(data_dir, "baseline.csv"))
write_csv(baseline_summary, file.path(data_dir, "baseline_summary.csv"))
write_csv(sofa_group_summary, file.path(data_dir, "sofa_group_summary.csv"))

dbWriteTable(
  con,
  Id(schema = "data_extract_crrt_strictk3", table = "baseline"),
  baseline,
  overwrite = TRUE
)

cat("\n=== STRICT K3 BASELINE SUMMARY ===\n")
print(t(baseline_summary), quote = FALSE)
cat("\n=== STRICT K3 SOFA GROUP SUMMARY ===\n")
print(sofa_group_summary)
cat("\nSaved strictK3 baseline outputs under:", data_dir, "\n")
