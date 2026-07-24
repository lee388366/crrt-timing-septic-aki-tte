# IPCW-weighted Cox: pooled interaction + stratified HRs; exploratory figures.

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

SCRIPT_DIR <- "analysis/exploratory_trajectory_6h"
source(file.path(SCRIPT_DIR, "00_constants.R"))
source(file.path(SCRIPT_DIR, "analysis_utils.R"))

feat <- readRDS(file.path(DIR_DATA, "trajectory_labels.rds")) %>%
  distinct(stay_id, .keep_all = TRUE)
baseline <- readRDS(PATH_BASELINE) %>% distinct(stay_id, .keep_all = TRUE)
ccw <- readRDS(PATH_CCW) %>% distinct(clone_id, .keep_all = TRUE)

res <- run_trajectory_analysis(
  feat,
  baseline,
  ccw,
  traj_col = "deterior",
  suffix = "main",
  equipoise_only = FALSE,
  legacy_main_03 = TRUE
)

# Equipoise-only interaction on same trajectory labels (legacy auxiliary log)
ccw_eq <- res$ccw_tr %>% filter(equipoise %in% TRUE)
if (nrow(ccw_eq) > 100 && n_distinct(ccw_eq$stay_id) > 50) {
  fit_eq <- coxph(
    Surv(follow_days, event) ~ trt * det,
    data = ccw_eq,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE
  )
  capture.output(
    print(summary(fit_eq)),
    file = file.path(DIR_LOG, "03_cox_interaction_equipoise.txt")
  )
}
