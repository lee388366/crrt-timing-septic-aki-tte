# Exploratory trajectory pipeline (run from repo root or this folder).
# source("exploratory_trajectory_6h/run_all.R")

SCRIPT_DIR <- "analysis/exploratory_trajectory_6h"

message("=== 01 SQL extract (optional) ===")
source(file.path(SCRIPT_DIR, "01_extract_uo_lac_sql.R"))

message("=== 02 Build trajectory labels ===")
source(file.path(SCRIPT_DIR, "02_build_trajectory_labels.R"))

message("=== 03 IPCW Cox + figures (main exploratory) ===")
source(file.path(SCRIPT_DIR, "03_ipcw_subgroup_cox_and_plots.R"))

message("=== 04 Variants: tight stable + equipoise-only ===")
source(file.path(SCRIPT_DIR, "04_variants_tight_equipoise.R"))

message("Done. Outputs under: ", file.path(SCRIPT_DIR, "output"))
