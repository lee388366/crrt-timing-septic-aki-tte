# Step 2 of input preparation: adds patient-level dates of death (mimiciv_hosp.patients.dod), removes the two records with
# inconsistent timestamps (CRRT start after death; t0 before ICU admission) and defines follow-up, early initiation [0, 24 h)
# and the SOFA 8-14 primary population. Code is verbatim from the cohort section of the original step3_bootstrap.R
# (its bootstrap section belonged to an intermediate specification and is not used). Output: private_data/cohort_v3_step3.rds
# Usage: Rscript 00_add_death_dates.R   (after step1_sofa_boundary.R)
suppressPackageStartupMessages({library(dplyr); library(lubridate); library(DBI); library(RPostgres)})
PKG  <- Sys.getenv("CRRT_REVISION_DIR"); PRIV <- file.path(PKG, "private_data")
coh_f <- file.path(PRIV, "cohort_v3_step3.rds")
if (!file.exists(coh_f)) {
  coh <- readRDS(file.path(PRIV, "cohort_v3_step1.rds"))
  con <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("PGDATABASE", "mimiciv31"),
    host = Sys.getenv("PGHOST", "localhost"), port = as.integer(Sys.getenv("PGPORT", "5432")),
    user = Sys.getenv("PGUSER", "postgres"), password = Sys.getenv("PGPASSWORD"))
  dod <- dbGetQuery(con, paste0("SELECT subject_id, dod FROM mimiciv_hosp.patients WHERE subject_id IN (",
                                paste(unique(coh$subject_id), collapse=","), ")")) %>% mutate(dod = as.Date(dod))
  dbDisconnect(con)
  bad <- (!is.na(coh$crrt_first_time) & !is.na(coh$deathtime) & coh$crrt_first_time > coh$deathtime) | coh$t0 < coh$icu_intime
  coh <- coh[!bad, ] %>% left_join(dod, by = "subject_id") %>%
    mutate(d_in = as.numeric(difftime(deathtime, t0, units = "hours")), d_dod = as.numeric(dod - as.Date(t0)) * 24,
      death_h = case_when(!is.na(d_in) ~ d_in, !is.na(d_dod) ~ d_dod, TRUE ~ NA_real_),
      died28 = !is.na(death_h) & death_h >= 0 & death_h <= 672,
      death_h = if_else(died28, pmax(death_h, 0.5), NA_real_), end_h = if_else(died28, death_h, 672),
      icu_out_h = as.numeric(difftime(icu_outtime, t0, units = "hours")),
      h_crrt = if_else(!is.na(crrt_hours_from_t0) & crrt_hours_from_t0 >= 0 & crrt_hours_from_t0 < end_h, crrt_hours_from_t0, NA_real_),
      early = !is.na(h_crrt) & h_crrt < 24,
      eq_strict = sofa_t0_strict >= 8 & sofa_t0_strict <= 14)
  saveRDS(coh, coh_f)
}
cat("cohort_v3_step3.rds:", nrow(readRDS(coh_f)), "patients\n")
