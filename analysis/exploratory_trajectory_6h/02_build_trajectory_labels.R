# Build composite "clinical deterioration" vs "relatively stable" (exploratory).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

SCRIPT_DIR <- "analysis/exploratory_trajectory_6h"
source(file.path(SCRIPT_DIR, "00_constants.R"))

baseline <- readRDS(PATH_BASELINE) %>%
  distinct(stay_id, .keep_all = TRUE)
tv <- readRDS(PATH_TV)

sql_path <- file.path(DIR_DATA, "sql_uo_lac_6h.rds")
if (!file.exists(sql_path)) {
  stop("Run 01_extract_uo_lac_sql.R first (creates sql_uo_lac_6h.rds).")
}
sql_pack <- readRDS(sql_path)
uo_6h <- sql_pack$uo_6h
lac_sql <- sql_pack$lac_first_6h
db_ok <- isTRUE(sql_pack$db_ok)

# NEE: max in hours 0–5 vs nee_at_t0
nee_agg <- tv %>%
  filter(hour_index %in% HOURS_POST_T0) %>%
  group_by(stay_id) %>%
  summarise(
    nee_max_0_6h = {
      x <- nee_h[is.finite(nee_h)]
      if (!length(x)) NA_real_ else max(x)
    },
    .groups = "drop"
  )

# Lactate fallback: max rolling lactate over hours 1–5 (proxy if no first BG)
lac_tv <- tv %>%
  filter(hour_index %in% 1:5) %>%
  group_by(stay_id) %>%
  summarise(
    lac_max_h1_5 = {
      x <- lactate_h[is.finite(lactate_h)]
      if (!length(x)) NA_real_ else max(x)
    },
    .groups = "drop"
  )

feat <- baseline %>%
  dplyr::select(stay_id, hadm_id, t0, nee_at_t0, sofa_equipoise) %>%
  left_join(nee_agg, by = "stay_id") %>%
  left_join(lac_tv, by = "stay_id") %>%
  left_join(
    uo_6h %>% dplyr::select(stay_id, uo_min_6h, uo_low_6h, uo_n_meas),
    by = "stay_id"
  ) %>%
  left_join(
    lac_sql %>% dplyr::select(stay_id, lac_first_6h, lac_high_first),
    by = "stay_id"
  ) %>%
  mutate(
    nee_at_t0 = replace_na(nee_at_t0, 0),
    nee_rise = !is.na(nee_max_0_6h) & (nee_max_0_6h > nee_at_t0 + NEE_INCREASE_EPS),
    # Prefer first ABG in window; else tv proxy
    lac_high = if_else(
      !is.na(lac_first_6h),
      lac_first_6h > LACTATE_THRESH,
      if_else(!is.na(lac_max_h1_5), lac_max_h1_5 > LACTATE_THRESH, NA)
    ),
    uo_low = if_else(!is.na(uo_low_6h), uo_low_6h, NA)
  )

feat <- feat %>%
  mutate(
    deterior = case_when(
      (uo_low %in% TRUE) | nee_rise | (lac_high %in% TRUE) ~ TRUE,
      (uo_low %in% FALSE | is.na(uo_low)) & !nee_rise & (lac_high %in% FALSE) ~ FALSE,
      TRUE ~ NA
    ),
    # Tight "stable": measured UO ≥ threshold, no NEE rise, lactate clearly not high (not NA).
    tight_stable = is.finite(uo_min_6h) &
      (uo_min_6h >= UO_THRESH_MLKGH) &
      !nee_rise &
      (lac_high %in% FALSE),
    clearly_deterior_tight = nee_rise |
      (lac_high %in% TRUE) |
      (is.finite(uo_min_6h) & (uo_min_6h < UO_THRESH_MLKGH)),
    deterior_tight = case_when(
      tight_stable ~ FALSE,
      clearly_deterior_tight ~ TRUE,
      TRUE ~ NA
    )
  )

saveRDS(feat, file.path(DIR_DATA, "trajectory_labels.rds"))
write_csv(feat, file.path(DIR_DATA, "trajectory_labels.csv"))

sink(file.path(DIR_LOG, "02_trajectory_summary.txt"))
cat("DB UO/lactate used:", db_ok, "\n")
cat("N patients:", nrow(feat), "\n\n")
cat("deterior (TRUE=bad trajectory, FALSE=stable, NA=indeterminate):\n")
print(table(feat$deterior, useNA = "ifany"))
cat("\nAmong non-NA deterior:\n")
print(prop.table(table(feat$deterior)))
cat("\nEquipoise x deterior:\n")
print(table(equipoise = feat$sofa_equipoise, deterior = feat$deterior, useNA = "ifany"))
cat("\ndeterior_tight (tight stable rule; NA = indeterminate):\n")
print(table(feat$deterior_tight, useNA = "ifany"))
cat("\nEquipoise x deterior_tight:\n")
print(table(equipoise = feat$sofa_equipoise, deterior_tight = feat$deterior_tight, useNA = "ifany"))
sink()

message("Saved trajectory_labels.rds")
