# Variants: tight trajectory definition + equipoise-only IPCW Cox (exploratory).

suppressPackageStartupMessages({
  library(dplyr)
})

SCRIPT_DIR <- "analysis/exploratory_trajectory_6h"
source(file.path(SCRIPT_DIR, "00_constants.R"))
source(file.path(SCRIPT_DIR, "analysis_utils.R"))

feat <- readRDS(file.path(DIR_DATA, "trajectory_labels.rds")) %>%
  distinct(stay_id, .keep_all = TRUE)
baseline <- readRDS(PATH_BASELINE) %>% distinct(stay_id, .keep_all = TRUE)
ccw <- readRDS(PATH_CCW) %>% distinct(clone_id, .keep_all = TRUE)

if (!"deterior_tight" %in% names(feat)) {
  stop("Re-run 02_build_trajectory_labels.R: column deterior_tight missing.")
}

run_trajectory_analysis(
  feat, baseline, ccw,
  traj_col = "deterior_tight",
  suffix = "tight",
  equipoise_only = FALSE,
  title_short = "| tight stable rule"
)

run_trajectory_analysis(
  feat, baseline, ccw,
  traj_col = "deterior",
  suffix = "equipoise",
  equipoise_only = TRUE,
  title_short = "| equipoise only"
)

run_trajectory_analysis(
  feat, baseline, ccw,
  traj_col = "deterior_tight",
  suffix = "equipoise_tight",
  equipoise_only = TRUE,
  title_short = "| equipoise + tight rule"
)

message("Variant analyses written under ", DIR_FIG, " and ", DIR_DATA)
