# Round 1: information timing (independent audit items M4, m1, m2)
#  (a) Weight convention: interval [k, k+1) uses the cumulative weight up to the start of hour k
#  (b) Hourly SOFA: endtime <= window end (previously starttime <= window end, which overlapped the decision hour); KDIGO: charttime < window end
#      creatinine_delta referenced to creatinine at or before t0; within-patient LOCF as in 03a
#  (c) Baseline laboratory values, MAP and urine output use values at or before t0 only (previous windows extended 1-4 h after t0)
#  (d) Variant: APS III removed from L0 (the MIMIC first-day APS III can include information after t0)
# Reports HR/RMST after each cumulative change and saves the corrected inputs for later rounds
source("cp_core.R")
suppressPackageStartupMessages({library(DBI); library(RPostgres); library(zoo)})
PKG <- Sys.getenv("CRRT_REVISION_DIR")
WS  <- Sys.getenv("CRRT_INPUT_DIR")
OUT <- file.path(PKG, "out"); PRIV <- file.path(PKG, "private_data")
sink(file.path(OUT, "round1_report.txt"), split = TRUE); cat("Round1 run at", format(Sys.time()), "\n\n")

coh <- readRDS(file.path(PRIV, "cohort_v3_step3.rds")) %>% mutate(era = factor(anchor_year_group))
tv0 <- readRDS(file.path(WS, "ch24_q4_corrected_inputs/tv_hourly.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE)
ids <- paste(unique(coh$stay_id), collapse = ","); hids <- paste(unique(coh$hadm_id), collapse = ",")

con <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("PGDATABASE", "mimiciv31"), host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")), user = Sys.getenv("PGUSER", "postgres"), password = Sys.getenv("PGPASSWORD"))
sofa_raw  <- dbGetQuery(con, paste0("SELECT stay_id, starttime, endtime, sofa_24hours AS s FROM mimiciv_derived.sofa WHERE sofa_24hours IS NOT NULL AND stay_id IN (", ids, ")"))
kdigo_raw <- dbGetQuery(con, paste0("SELECT stay_id, charttime, aki_stage FROM mimiciv_derived.kdigo_stages WHERE aki_stage IS NOT NULL AND stay_id IN (", ids, ")"))
creat_raw <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, creatinine FROM mimiciv_derived.chemistry WHERE creatinine IS NOT NULL AND creatinine < 25 AND hadm_id IN (", hids, ")"))
lact_raw  <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, lactate FROM mimiciv_derived.bg WHERE lactate IS NOT NULL AND hadm_id IN (", hids, ")"))
uo_raw    <- dbGetQuery(con, paste0("SELECT stay_id, charttime, uo_mlkghr_6hr AS uo FROM mimiciv_derived.urine_output_rate WHERE uo_mlkghr_6hr IS NOT NULL AND stay_id IN (", ids, ")"))
map_raw   <- dbGetQuery(con, paste0("SELECT stay_id, charttime, mbp AS map FROM mimiciv_derived.vitalsign WHERE mbp IS NOT NULL AND stay_id IN (", ids, ")"))
dbDisconnect(con)

# ---------- (c) Baseline values at or before t0 ----------
last_before <- function(raw, key, val, lower_h) {
  raw %>% inner_join(coh %>% select(stay_id, hadm_id, t0), by = key, relationship = "many-to-many") %>%
    filter(charttime <= t0, charttime >= t0 - hours(lower_h)) %>% group_by(stay_id) %>%
    slice_max(charttime, n = 1, with_ties = FALSE) %>% ungroup() %>% select(stay_id, all_of(val))
}
b_cr  <- last_before(creat_raw, "hadm_id", "creatinine", 1e5) %>% rename(creatinine_t0_v2 = creatinine)
b_la  <- last_before(lact_raw,  "hadm_id", "lactate", 12)     %>% rename(lactate_t0_v2 = lactate)
b_uo  <- last_before(uo_raw %>% left_join(coh %>% select(stay_id, hadm_id), by = "stay_id"), "stay_id", "uo", 1e5) %>% rename(uo_rate_t0_v2 = uo)
b_map <- last_before(map_raw %>% left_join(coh %>% select(stay_id, hadm_id), by = "stay_id"), "stay_id", "map", 2) %>% rename(map_t0_v2 = map)
cohv2 <- coh %>% left_join(b_cr, by = "stay_id") %>% left_join(b_la, by = "stay_id") %>% left_join(b_uo, by = "stay_id") %>% left_join(b_map, by = "stay_id")
cat("Missing baseline values (previous -> at or before t0):\n")
for (v in c("creatinine_t0", "lactate_t0", "uo_rate_t0", "map_t0"))
  cat(sprintf("  %-14s %4d -> %4d | proportion unchanged %.3f\n", v, sum(is.na(cohv2[[v]])), sum(is.na(cohv2[[paste0(v, "_v2")]])),
      mean(cohv2[[v]] == cohv2[[paste0(v, "_v2")]], na.rm = TRUE)))

# ---------- (b) Hourly SOFA / KDIGO / creatinine change ----------
grid <- coh %>% select(stay_id, hadm_id, t0) %>% crossing(hour_index = 0:23) %>% mutate(window_end = t0 + hours(hour_index + 1))
sofa_v2 <- grid %>% select(stay_id, hour_index, window_end) %>% inner_join(sofa_raw, by = "stay_id", relationship = "many-to-many") %>%
  filter(endtime <= window_end) %>% group_by(stay_id, hour_index) %>% slice_max(starttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(stay_id, hour_index, sofa_h_v2 = s)
kdigo_v2 <- grid %>% select(stay_id, hour_index, window_end) %>% inner_join(kdigo_raw, by = "stay_id", relationship = "many-to-many") %>%
  filter(charttime < window_end) %>% group_by(stay_id, hour_index) %>% slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(stay_id, hour_index, kdigo_h_v2 = aki_stage)
cr_h <- grid %>% select(stay_id, hadm_id, hour_index, window_end) %>% inner_join(creat_raw, by = "hadm_id", relationship = "many-to-many") %>%
  filter(charttime >= window_end - hours(6), charttime < window_end) %>% group_by(stay_id, hour_index) %>%
  summarise(cr = mean(creatinine), .groups = "drop")
tv1 <- tv0 %>% left_join(sofa_v2, by = c("stay_id", "hour_index")) %>% left_join(kdigo_v2, by = c("stay_id", "hour_index")) %>%
  left_join(cr_h, by = c("stay_id", "hour_index")) %>% left_join(b_cr, by = "stay_id") %>%
  mutate(cd_v2 = cr - coalesce(creatinine_t0_v2, cr)) %>%
  group_by(stay_id) %>% arrange(hour_index, .by_group = TRUE) %>%
  mutate(across(c(sofa_h_v2, kdigo_h_v2, cd_v2), ~ zoo::na.locf(., na.rm = FALSE))) %>% ungroup()
cat("\nHourly covariates (previous -> corrected): proportion unchanged SOFA %.3f | KDIGO %.3f | creatinine change %.3f\n" |> sprintf(
  mean(tv1$sofa_h == tv1$sofa_h_v2, na.rm = TRUE), mean(tv1$kdigo_h == tv1$kdigo_h_v2, na.rm = TRUE),
  mean(abs(tv1$creatinine_delta - tv1$cd_v2) < 1e-9, na.rm = TRUE)))
cat(sprintf("Missing (previous -> corrected): SOFA %d -> %d | KDIGO %d -> %d | creatinine change %d -> %d\n",
  sum(is.na(tv1$sofa_h)), sum(is.na(tv1$sofa_h_v2)), sum(is.na(tv1$kdigo_h)), sum(is.na(tv1$kdigo_h_v2)),
  sum(is.na(tv1$creatinine_delta)), sum(is.na(tv1$cd_v2))))
tv_v2 <- tv1 %>% transmute(stay_id, hour_index, map_h, nee_h, lactate_h, sofa_h = sofa_h_v2, kdigo_h = kdigo_h_v2,
                           creatinine_h, creatinine_delta = cd_v2, fluid_balance_h)
cohv2b <- cohv2 %>% mutate(creatinine_t0 = creatinine_t0_v2, lactate_t0 = lactate_t0_v2, uo_rate_t0 = uo_rate_t0_v2, map_t0 = map_t0_v2)
saveRDS(tv_v2, file.path(PRIV, "tv_hourly_round1.rds")); saveRDS(cohv2b, file.path(PRIV, "cohort_round1.rds"))

# ---------- Cumulative estimates ----------
L0_full <- L0
run <- function(lab, pop, tvl, wlag) {
  pop$id <- seq_len(nrow(pop)); cp <- build(pop, tvl, wlag = wlag); e <- est(cp); d <- attr(cp, "diag")
  cat(sprintf("  %-40s HR %.3f (%.3f-%.3f, robust) | RMST diff %5.2f | ESS_A %.1f/%d | max|SMD| %.3f | max w A %.2f B %.2f\n",
      lab, e["HR"], e["lo"], e["hi"], e["rmst_diff"], d["ESS_A"], d["n_A_post"], d["max_absSMD_A"], d["max_wA"], d["max_wB"]))
  data.frame(step = lab, t(e), t(d))
}
res <- list()
for (sub in c("primary", "full")) {
  P0 <- if (sub == "primary") filter(coh, eq_strict) else coh
  P2 <- if (sub == "primary") filter(cohv2b, eq_strict) else cohv2b
  cat("\n== ", sub, " (n=", nrow(P0), ")\n", sep = "")
  L0 <- L0_full
  res[[length(res) + 1]] <- cbind(pop = sub, run("0 previous method (reproduces step 2)", P0, tv0, FALSE))
  res[[length(res) + 1]] <- cbind(pop = sub, run("(a) + weight at interval start", P0, tv0, TRUE))
  res[[length(res) + 1]] <- cbind(pop = sub, run("(b) + strictly lagged hourly SOFA/KDIGO", P0, tv_v2, TRUE))
  res[[length(res) + 1]] <- cbind(pop = sub, run("(c) + baseline values at or before t0", P2, tv_v2, TRUE))
  L0 <- setdiff(L0_full, "apsiii")
  res[[length(res) + 1]] <- cbind(pop = sub, run("(d) variant: (c) without APS III in L0", P2, tv_v2, TRUE))
  L0 <- L0_full
}
write.csv(bind_rows(res), file.path(OUT, "round1_estimates.csv"), row.names = FALSE)
cat("\nCorrected inputs saved: private_data/tv_hourly_round1.rds, cohort_round1.rds\nDONE\n"); sink()
