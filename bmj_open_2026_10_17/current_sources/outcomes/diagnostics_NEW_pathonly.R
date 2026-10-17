# Strict clean first-KDIGO3 IPCW weight and balance diagnostics.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

root_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
data_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_inputs")
out_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_diagnostics/NEW")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

baseline <- readRDS(file.path(data_dir, "baseline.rds"))
ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

vars <- c(
  "age", "sofa_t0", "apsiii", "creatinine_baseline", "nee_at_t0",
  "lactate_t0", "uo_rate_t0", "map_t0", "cci"
)

smd_unweighted <- function(x, g) {
  a <- x[g %in% TRUE]
  b <- x[g %in% FALSE]
  if (sum(is.finite(a)) < 2 || sum(is.finite(b)) < 2) return(NA_real_)
  abs(mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) /
    sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
}

smd_weighted <- function(df, var) {
  a <- df %>% filter(strategy == "A_early", !is.na(.data[[var]]), !is.na(ipcw_w_trim))
  b <- df %>% filter(strategy == "B_deferred", !is.na(.data[[var]]), !is.na(ipcw_w_trim))
  if (nrow(a) < 2 || nrow(b) < 2) return(NA_real_)
  ma <- weighted.mean(a[[var]], a$ipcw_w_trim, na.rm = TRUE)
  mb <- weighted.mean(b[[var]], b$ipcw_w_trim, na.rm = TRUE)
  abs(ma - mb) / sqrt((var(a[[var]], na.rm = TRUE) + var(b[[var]], na.rm = TRUE)) / 2)
}

balance <- bind_rows(lapply(c("Full cohort", "SOFA 8-14"), function(pop) {
  base_pop <- if (pop == "SOFA 8-14") baseline %>% filter(sofa_equipoise %in% TRUE) else baseline
  ccw_pop <- if (pop == "SOFA 8-14") ccw %>% filter(equipoise %in% TRUE) else ccw
  tibble(
    population = pop,
    variable = vars,
    pre_smd = sapply(vars, function(v) smd_unweighted(base_pop[[v]], base_pop$early_crrt)),
    post_weighted_smd = sapply(vars, function(v) smd_weighted(ccw_pop, v))
  )
})) %>%
  mutate(
    pre_flag_gt_0_10 = pre_smd > 0.10,
    post_flag_gt_0_10 = post_weighted_smd > 0.10
  )

weight_summary <- ccw %>%
  group_by(strategy, equipoise) %>%
  summarise(
    n_rows = n(),
    n_patients = n_distinct(stay_id),
    mean_w = mean(ipcw_w_trim, na.rm = TRUE),
    sd_w = sd(ipcw_w_trim, na.rm = TRUE),
    min_w = min(ipcw_w_trim, na.rm = TRUE),
    p01_w = unname(quantile(ipcw_w_trim, 0.01, na.rm = TRUE)),
    p25_w = unname(quantile(ipcw_w_trim, 0.25, na.rm = TRUE)),
    median_w = median(ipcw_w_trim, na.rm = TRUE),
    p75_w = unname(quantile(ipcw_w_trim, 0.75, na.rm = TRUE)),
    p99_w = unname(quantile(ipcw_w_trim, 0.99, na.rm = TRUE)),
    max_w = max(ipcw_w_trim, na.rm = TRUE),
    ess = sum(ipcw_w_trim, na.rm = TRUE)^2 / sum(ipcw_w_trim^2, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(balance, file.path(out_dir, "strictK3_v2_ipcw_balance_smd.csv"))
write_csv(weight_summary, file.path(out_dir, "strictK3_v2_ipcw_weight_summary.csv"))

report <- c(
  "# Strict KDIGO3 IPCW Diagnostics",
  "",
  "## Weight Summary",
  "",
  "| strategy | equipoise | n patients | mean | p99 | max | ESS |",
  "|---|---:|---:|---:|---:|---:|---:|",
  apply(weight_summary, 1, function(x) sprintf(
    "| %s | %s | %s | %.3f | %.3f | %.3f | %.1f |",
    x[["strategy"]], x[["equipoise"]], x[["n_patients"]],
    as.numeric(x[["mean_w"]]), as.numeric(x[["p99_w"]]),
    as.numeric(x[["max_w"]]), as.numeric(x[["ess"]])
  )),
  "",
  "## Balance Summary",
  "",
  "| population | variables with pre SMD >0.10 | variables with post SMD >0.10 | max post SMD |",
  "|---|---:|---:|---:|",
  apply(balance %>% group_by(population) %>% summarise(
    pre_gt_010 = sum(pre_flag_gt_0_10, na.rm = TRUE),
    post_gt_010 = sum(post_flag_gt_0_10, na.rm = TRUE),
    max_post_smd = max(post_weighted_smd, na.rm = TRUE),
    .groups = "drop"
  ), 1, function(x) sprintf(
    "| %s | %s | %s | %.3f |",
    x[["population"]], x[["pre_gt_010"]], x[["post_gt_010"]],
    as.numeric(x[["max_post_smd"]])
  )),
  "",
  "Detailed balance table: `strictK3_v2_ipcw_balance_smd.csv`."
)
writeLines(report, file.path(out_dir, "strictK3_v2_ipcw_diagnostics_report.md"))

print(weight_summary)
print(balance)
