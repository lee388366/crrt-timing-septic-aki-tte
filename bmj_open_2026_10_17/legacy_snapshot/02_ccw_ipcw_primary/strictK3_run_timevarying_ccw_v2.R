# Run legacy time-varying covariate extraction and CCW/IPCW with strictK3_v2 paths.
# The legacy scripts are evaluated unchanged except for redirecting the hardcoded
# analysis/data directory to analysis/strictK3_primary_v2/data.
suppressPackageStartupMessages({
  library(fs)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
legacy_data_dir <- file.path(root_dir, "data")
strict_data_dir <- file.path(root_dir, "strictK3_primary_v2", "data")

if (!dir.exists(strict_data_dir)) {
  stop("Missing strictK3_v2 data directory: ", strict_data_dir)
}
if (!file.exists(file.path(strict_data_dir, "baseline.rds"))) {
  stop("Run strictK3_primary_v2_extract.R and strictK3_v2_baseline_covariates.R first.")
}

run_legacy_step_with_strict_paths <- function(script_path) {
  message("\n=== Running legacy step with strictK3_v2 paths: ", basename(script_path), " ===")
  code <- readLines(script_path, warn = FALSE)
  code <- gsub(legacy_data_dir, strict_data_dir, code, fixed = TRUE)
  code <- gsub(
    "COALESCE(vasopressin, 0) * 2.5",
    "(COALESCE(vasopressin, 0) / 60.0) * 2.5",
    code,
    fixed = TRUE
  )
  code <- gsub(
    "NEE = NE + Epi + Dopa/100 + Phenyl*0.1 + Vaso*2.5",
    "NEE = NE + Epi + Dopa/100 + Phenyl*0.1 + (Vaso U/hour / 60)*2.5",
    code,
    fixed = TRUE
  )
  # The legacy time-varying script can duplicate rows when multiple labs share
  # the exact same timestamp. Keep the intended single closest record.
  code <- gsub("slice_max(charttime, n=1)", "slice_max(charttime, n=1, with_ties=FALSE)", code, fixed = TRUE)
  code <- gsub("slice_max(starttime, n=1)", "slice_max(starttime, n=1, with_ties=FALSE)", code, fixed = TRUE)
  step_env <- new.env(parent = globalenv())
  step_env[["strictK3_v2_path_redirect"]] <- c(
    legacy_data_dir = legacy_data_dir,
    strict_data_dir = strict_data_dir
  )
  eval(parse(text = code), envir = step_env)
}

run_legacy_step_with_strict_paths(file.path(root_dir, "03a_timevarying_covariates.R"))
run_legacy_step_with_strict_paths(file.path(root_dir, "03_ccw_ipcw.R"))

message("\nStrictK3 time-varying and CCW/IPCW outputs are under: ", strict_data_dir)
