# Step 1: input checks, corrected SOFA-at-t0 window, 24-hour boundary and inconsistent records
# Usage: Rscript step1_sofa_boundary.R
# Writes only aggregate statistics to out/; patient-level intermediate data are written to private_data/ only (never shared)
suppressPackageStartupMessages({library(dplyr); library(lubridate); library(DBI); library(RPostgres)})

PKG  <- Sys.getenv("CRRT_REVISION_DIR")
IN   <- file.path(Sys.getenv("CRRT_INPUT_DIR"), "ch24_q4_corrected_inputs")
OUT  <- file.path(PKG, "out"); PRIV <- file.path(PKG, "private_data")
log_f <- file.path(OUT, "step1_report.txt"); sink(log_f, split = TRUE)
cat("Step1 run at", format(Sys.time()), "\n\n")

# ---------- 0. Input structure (column names / row counts only) ----------
baseline  <- readRDS(file.path(IN, "baseline.rds"))
tv_hourly <- readRDS(file.path(IN, "tv_hourly.rds"))
tv_hour23 <- readRDS(file.path(IN, "tv_hour23.rds"))
ccw_old   <- readRDS(file.path(IN, "ccw_final.rds"))
for (nm in c("baseline","tv_hourly","tv_hour23","ccw_old")) {
  x <- get(nm); cat("==", nm, ":", nrow(x), "rows x", ncol(x), "cols\n"); cat(paste(names(x), collapse=", "), "\n\n")
}

# Cohort of the original analysis (as in 03_ccw_ipcw.R)
cohort <- baseline %>% group_by(subject_id) %>% slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup() %>% filter(!is.na(sofa_t0))
cat("Original analysis cohort n =", nrow(cohort), "| early =", sum(cohort$early_crrt, na.rm=TRUE),
    "| SOFA 8-14 =", sum(cohort$sofa_t0 >= 8 & cohort$sofa_t0 <= 14), "\n\n")

# ---------- 1. SOFA window: previously starttime <= t0 + 1 h (window end t0 + 2 h) -> corrected to endtime <= t0 ----------
con <- dbConnect(RPostgres::Postgres(),
  dbname = Sys.getenv("PGDATABASE", "mimiciv31"), host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")), user = Sys.getenv("PGUSER", "postgres"),
  password = Sys.getenv("PGPASSWORD"))
on.exit(dbDisconnect(con), add = TRUE)
ids <- paste(unique(cohort$stay_id), collapse = ",")
sofa_h <- dbGetQuery(con, paste0("SELECT stay_id, starttime, endtime, sofa_24hours FROM mimiciv_derived.sofa WHERE stay_id IN (", ids, ")"))

pick <- function(rule) {
  s <- sofa_h %>% inner_join(cohort %>% select(stay_id, t0), by = "stay_id")
  s <- if (rule == "orig") filter(s, starttime <= t0 + hours(1)) else filter(s, endtime <= t0)
  s %>% group_by(stay_id) %>% slice_max(starttime, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(stay_id, v = sofa_24hours)
}
s_orig <- pick("orig") %>% rename(sofa_orig_requery = v)
s_strict <- pick("strict") %>% rename(sofa_t0_strict = v)
cmp <- cohort %>% select(stay_id, sofa_t0, early_crrt) %>% left_join(s_orig, by="stay_id") %>% left_join(s_strict, by="stay_id")
cat("Proportion of re-queried original-rule values equal to stored sofa_t0:", round(mean(cmp$sofa_orig_requery == cmp$sofa_t0, na.rm=TRUE), 4),
    "| re-query missing:", sum(is.na(cmp$sofa_orig_requery)), "\n")
cat("Missing under strict rule (endtime <= t0):", sum(is.na(cmp$sofa_t0_strict)), "\n")
cat("SOFA change (strict - original):\n"); print(table(cmp$sofa_t0_strict - cmp$sofa_t0, useNA = "ifany"))
cmp <- cmp %>% mutate(eq_old = sofa_t0 >= 8 & sofa_t0 <= 14,
                      eq_new = !is.na(sofa_t0_strict) & sofa_t0_strict >= 8 & sofa_t0_strict <= 14)
cat("\nSOFA 8-14 membership (rows = original, columns = strict):\n"); print(table(old = cmp$eq_old, strict = cmp$eq_new))
cat("SOFA 8-14 under strict rule: n =", sum(cmp$eq_new), "| early =", sum(cmp$eq_new & cmp$early_crrt), "\n\n")

# ---------- 2. 24-hour boundary and inconsistent records ----------
h <- cohort$crrt_hours_from_t0
cat("CRRT within [23,24] h:", sum(h >= 23 & h <= 24, na.rm=TRUE), "| exactly 24 h:", sum(h == 24, na.rm=TRUE),
    "| early group with h >= 24:", sum(cohort$early_crrt & h >= 24, na.rm=TRUE), "\n")
if ("deathtime" %in% names(cohort) && "crrt_first_time" %in% names(cohort)) {
  cat("CRRT start after deathtime:", sum(cohort$crrt_first_time > cohort$deathtime, na.rm=TRUE), "\n")
} else cat("baseline lacks deathtime/crrt_first_time; CRRT-death order check skipped\n")
cat("t0 before ICU admission:", sum(cohort$t0 < cohort$icu_intime, na.rm=TRUE), "\n")
cat("tv_hourly hours per patient:\n"); print(summary(as.numeric(table(tv_hourly$stay_id))))

# ---------- Save corrected cohort (local only) ----------
cohort_v3 <- cohort %>% left_join(s_strict, by = "stay_id") %>%
  mutate(early_crrt_strict = !is.na(crrt_hours_from_t0) & crrt_hours_from_t0 >= 0 & crrt_hours_from_t0 < 24)
saveRDS(cohort_v3, file.path(PRIV, "cohort_v3_step1.rds"))
cat("\nEarly initiators under [0, 24):", sum(cohort_v3$early_crrt_strict), "(original", sum(cohort$early_crrt), ")\n")
cat("\nDONE\n"); sink()
