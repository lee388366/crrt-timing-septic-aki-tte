# =============================================================================
# Anchor sensitivity overlap/balance diagnostics
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary")
results_dir <- file.path(strict_dir, "results")

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

one_anchor <- function(anchor) {
  data_dir <- file.path(strict_dir, paste0("anchor_sensitivity_", anchor), "data")
  baseline <- readRDS(file.path(data_dir, "baseline.rds"))
  ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
    distinct(clone_id, .keep_all = TRUE)

  bind_rows(lapply(c("Full cohort", "SOFA 8-14"), function(pop) {
    base_pop <- if (pop == "SOFA 8-14") baseline %>% filter(sofa_equipoise %in% TRUE) else baseline
    ccw_pop <- if (pop == "SOFA 8-14") ccw %>% filter(equipoise %in% TRUE) else ccw
    tibble(
      anchor = anchor,
      population = pop,
      variable = vars,
      pre_smd = sapply(vars, function(v) smd_unweighted(base_pop[[v]], base_pop$early_crrt)),
      post_weighted_smd = sapply(vars, function(v) smd_weighted(ccw_pop, v))
    )
  }))
}

balance <- bind_rows(lapply(c("KDIGO3_minus6h", "KDIGO2"), one_anchor)) %>%
  mutate(
    pre_flag_gt_0_10 = pre_smd > 0.10,
    post_flag_gt_0_10 = post_weighted_smd > 0.10
  )

summary <- balance %>%
  group_by(anchor, population) %>%
  summarise(
    variables_pre_smd_gt_0_10 = sum(pre_flag_gt_0_10, na.rm = TRUE),
    variables_post_smd_gt_0_10 = sum(post_flag_gt_0_10, na.rm = TRUE),
    max_post_smd = max(post_weighted_smd, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(balance, file.path(results_dir, "strictK3_anchor_sensitivity_balance_smd.csv"))
write_csv(summary, file.path(results_dir, "strictK3_anchor_sensitivity_balance_summary.csv"))

report <- c(
  "# Anchor sensitivity overlap and balance diagnostics",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "These diagnostics use each anchor-specific CCW/IPCW output directory and report pre/post-weighting SMDs for the same core baseline covariates as the strictK3 primary diagnostics.",
  "",
  paste(capture.output(print(summary, n = nrow(summary))), collapse = "\n"),
  "",
  "Outputs:",
  "",
  paste0("- `", file.path(results_dir, "strictK3_anchor_sensitivity_balance_smd.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_anchor_sensitivity_balance_summary.csv"), "`")
)

writeLines(report, file.path(results_dir, "strictK3_anchor_sensitivity_diagnostics_report.md"))
cat("Saved anchor diagnostics:", file.path(results_dir, "strictK3_anchor_sensitivity_diagnostics_report.md"), "\n")
