# Optional: pull urine-output + first ABG lactate in (t0, t0+6h] from local MIMIC-IV Postgres.
# If DB unavailable, writes placeholder NAs — 02_build_trajectory_labels.R uses tv_hourly fallbacks.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

SCRIPT_DIR <- "analysis/exploratory_trajectory_6h"
source(file.path(SCRIPT_DIR, "00_constants.R"))

logf <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(..., collapse = ""), "\n")
  cat(msg)
  write(msg, file.path(DIR_LOG, "01_sql_extract.log"), append = TRUE)
}

baseline <- readRDS(PATH_BASELINE)
cohort_key <- baseline %>%
  dplyr::select(stay_id, hadm_id, t0) %>%
  distinct()

uo_6h <- NULL
lac_first_6h <- NULL
db_ok <- FALSE

tryCatch(
  {
    if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RPostgres", quietly = TRUE)) {
      stop("DBI/RPostgres not installed")
    }
    con <- DBI::dbConnect(
      RPostgres::Postgres(),
      dbname = "mimiciv31",
      host = "localhost",
      port = 5432,
      user = "postgres",
      password = Sys.getenv("PGPASSWORD")
    )
    on.exit(DBI::dbDisconnect(con), add = TRUE)

    # --- Urine: any measurement in (t0, t0+6h]; flag if minimum 6h rate < threshold
    uo_raw <- DBI::dbGetQuery(con, "
      SELECT stay_id, charttime, uo_mlkghr_6hr AS uo_rate
      FROM mimiciv_derived.urine_output_rate
      WHERE uo_mlkghr_6hr IS NOT NULL
    ")
    uo_6h <- uo_raw %>%
      inner_join(cohort_key, by = "stay_id", relationship = "many-to-many") %>%
      filter(charttime > t0, charttime <= t0 + as.difftime(6, units = "hours")) %>%
      group_by(stay_id) %>%
      summarise(
        uo_min_6h = min(uo_rate, na.rm = TRUE),
        uo_n_meas = dplyr::n(),
        .groups = "drop"
      ) %>%
      mutate(uo_low_6h = uo_min_6h < !!UO_THRESH_MLKGH)

    # --- Lactate: first BG strictly after t0 within 6h
    lac_raw <- DBI::dbGetQuery(con, "
      SELECT hadm_id, charttime, lactate
      FROM mimiciv_derived.bg
      WHERE lactate IS NOT NULL AND lactate < 30
    ")
    lac_first_6h <- lac_raw %>%
      inner_join(cohort_key, by = "hadm_id", relationship = "many-to-many") %>%
      filter(charttime > t0, charttime <= t0 + as.difftime(6, units = "hours")) %>%
      group_by(stay_id) %>%
      arrange(charttime, .by_group = TRUE) %>%
      slice(1L) %>%
      ungroup() %>%
      dplyr::select(stay_id, lac_first_time = charttime, lac_first_6h = lactate) %>%
      mutate(lac_high_first = lac_first_6h > !!LACTATE_THRESH)

    db_ok <- TRUE
    logf("Postgres OK: uo rows ", nrow(uo_6h), ", lac rows ", nrow(lac_first_6h))
  },
  error = function(e) {
    logf("DB extract skipped: ", conditionMessage(e))
    uo_6h <<- tibble(stay_id = integer(), uo_min_6h = double(), uo_n_meas = integer(), uo_low_6h = logical())
    lac_first_6h <<- tibble(stay_id = integer(), lac_first_time = as.POSIXct(character()), lac_first_6h = double(), lac_high_first = logical())
  }
)

saveRDS(
  list(
    uo_6h = uo_6h,
    lac_first_6h = lac_first_6h,
    db_ok = db_ok,
    uo_thresh = UO_THRESH_MLKGH,
    lac_thresh = LACTATE_THRESH
  ),
  file.path(DIR_DATA, "sql_uo_lac_6h.rds")
)

write_csv(
  if (nrow(uo_6h) > 0) uo_6h else tibble(stay_id = NA_integer_, note = "empty"),
  file.path(DIR_DATA, "uo_6h_summary.csv")
)
write_csv(
  if (nrow(lac_first_6h) > 0) lac_first_6h else tibble(stay_id = NA_integer_, note = "empty"),
  file.path(DIR_DATA, "lac_first_6h_summary.csv")
)

logf("Saved ", file.path(DIR_DATA, "sql_uo_lac_6h.rds"))
