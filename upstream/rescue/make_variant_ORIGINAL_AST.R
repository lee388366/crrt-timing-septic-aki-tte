make_variant <- function(name, cohort) {
  variant_dir <- file.path(strict_dir, paste0("rescue_sensitivity_", name))
  data_dir <- file.path(variant_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(cohort, file.path(data_dir, "cohort_base.rds"))
  write_csv(cohort, file.path(data_dir, "cohort_base.csv"))
  build_baseline(cohort, data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03a_timevarying_covariates.R"), data_dir)
  run_legacy_step_with_data_dir(file.path(root_dir, "03_ccw_ipcw.R"), data_dir)

  ccw <- readRDS(file.path(data_dir, "ccw_final.rds"))
  bind_rows(
    run_cox(ccw, paste0(name, ": full cohort")),
    run_cox(filter(ccw, equipoise %in% TRUE), paste0(name, ": SOFA 8-14"))
  )
}
