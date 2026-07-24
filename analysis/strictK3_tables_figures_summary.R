# =============================================================================
# Strict first-KDIGO3: manuscript-facing tables and figures
#   - Table 1 with SMD flags
#   - Love plot for pre/post IPCW balance
#   - Weighted cumulative mortality curve for primary SOFA 8-14
#   - Forest plot across primary/sensitivity/anchor analyses
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(tibble)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
fig_dir <- file.path(strict_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(file.path(data_dir, "baseline.rds"))
ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

smd_cont <- function(x, g) {
  a <- x[g %in% TRUE]
  b <- x[g %in% FALSE]
  if (sum(is.finite(a)) < 2 || sum(is.finite(b)) < 2) return(NA_real_)
  abs(mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) /
    sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
}

smd_bin <- function(x, g) {
  x <- as.numeric(x %in% TRUE)
  a <- x[g %in% TRUE]
  b <- x[g %in% FALSE]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  p1 <- mean(a, na.rm = TRUE)
  p0 <- mean(b, na.rm = TRUE)
  p <- (p1 + p0) / 2
  abs(p1 - p0) / sqrt(p * (1 - p))
}

fmt_cont <- function(x) {
  sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE))
}

fmt_bin <- function(x) {
  sprintf("%d (%.1f%%)", sum(x %in% TRUE, na.rm = TRUE), 100 * mean(x %in% TRUE, na.rm = TRUE))
}

table1_vars <- tribble(
  ~variable, ~label, ~type,
  "age", "Age, years", "continuous",
  "gender_female", "Female sex", "binary",
  "sofa_t0", "SOFA at t0", "continuous",
  "apsiii", "APS III", "continuous",
  "creatinine_baseline", "Baseline creatinine, mg/dL", "continuous",
  "creatinine_t0", "Creatinine at t0, mg/dL", "continuous",
  "lactate_t0", "Lactate at t0, mmol/L", "continuous",
  "nee_at_t0", "Norepinephrine-equivalent dose at t0", "continuous",
  "uo_rate_t0", "Urine-output rate at t0, mL/kg/h", "continuous",
  "map_t0", "MAP at t0, mmHg", "continuous",
  "cci", "Charlson comorbidity index", "continuous",
  "stage3_by_creat", "KDIGO3 by creatinine", "binary",
  "stage3_by_uo", "KDIGO3 by urine output", "binary",
  "stage3_by_crrt", "KDIGO3 by CRRT criterion", "binary"
)

table1_df <- baseline %>%
  filter(sofa_equipoise %in% TRUE) %>%
  mutate(gender_female = gender == "F")

table1 <- bind_rows(lapply(seq_len(nrow(table1_vars)), function(i) {
  v <- table1_vars$variable[i]
  label <- table1_vars$label[i]
  type <- table1_vars$type[i]
  early <- table1_df %>% filter(early_crrt %in% TRUE)
  deferred <- table1_df %>% filter(!(early_crrt %in% TRUE))

  if (type == "continuous") {
    smd <- smd_cont(table1_df[[v]], table1_df$early_crrt)
    early_text <- fmt_cont(early[[v]])
    deferred_text <- fmt_cont(deferred[[v]])
  } else {
    smd <- smd_bin(table1_df[[v]], table1_df$early_crrt)
    early_text <- fmt_bin(early[[v]])
    deferred_text <- fmt_bin(deferred[[v]])
  }

  tibble(
    variable = label,
    early = early_text,
    deferred = deferred_text,
    smd = round(smd, 3),
    smd_flag_ge_0_10 = smd >= 0.10
  )
}))

write_csv(table1, file.path(results_dir, "strictK3_table1_sofa8_14.csv"))

table1_md <- c(
  "# Strict first-KDIGO3 Table 1: SOFA 8-14 primary equipoise",
  "",
  paste0("Population: n=", nrow(table1_df), "; early n=", sum(table1_df$early_crrt), "; deferred n=", sum(!table1_df$early_crrt), "."),
  "",
  "| variable | Early CRRT | Deferred CRRT | SMD | SMD >=0.10 |",
  "|---|---:|---:|---:|---:|",
  apply(table1, 1, function(x) sprintf(
    "| %s | %s | %s | %.3f | %s |",
    x[["variable"]], x[["early"]], x[["deferred"]], as.numeric(x[["smd"]]),
    ifelse(x[["smd_flag_ge_0_10"]] %in% c(TRUE, "TRUE"), "YES", "")
  ))
)
writeLines(table1_md, file.path(results_dir, "strictK3_table1_sofa8_14.md"))

# Love plot --------------------------------------------------------------------
balance <- read_csv(file.path(results_dir, "strictK3_ipcw_balance_smd.csv"), show_col_types = FALSE) %>%
  filter(population == "SOFA 8-14") %>%
  mutate(variable = factor(variable, levels = rev(variable)))

love_df <- balance %>%
  select(variable, pre_smd, post_weighted_smd) %>%
  tidyr::pivot_longer(c(pre_smd, post_weighted_smd), names_to = "stage", values_to = "smd") %>%
  mutate(stage = recode(stage, pre_smd = "Before weighting", post_weighted_smd = "After IPCW"))

p_love <- ggplot(love_df, aes(x = smd, y = variable, color = stage, shape = stage)) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.7) +
  scale_color_manual(values = c("Before weighting" = "#C0392B", "After IPCW" = "#2E86C1")) +
  labs(
    title = "Strict first-KDIGO3 primary SOFA 8-14: covariate balance",
    subtitle = "Dashed line marks SMD = 0.10",
    x = "Absolute standardized mean difference",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

png(file.path(fig_dir, "strictK3_love_plot_sofa8_14.png"), width = 2100, height = 1500, res = 300)
print(p_love)
dev.off()
pdf(file.path(fig_dir, "strictK3_love_plot_sofa8_14.pdf"), width = 7, height = 5, onefile = FALSE)
print(p_love)
dev.off()

# Weighted cumulative mortality ------------------------------------------------
km_df <- ccw %>%
  filter(equipoise %in% TRUE) %>%
  mutate(strategy_label = factor(ifelse(strategy == "A_early", "Early CRRT strategy", "Deferred CRRT strategy")))

sf <- survfit(Surv(follow_days, event) ~ strategy_label, data = km_df, weights = ipcw_w_trim)
sf_sum <- summary(sf, times = seq(0, 28, by = 1))
curve_df <- tibble(
  time = sf_sum$time,
  survival = sf_sum$surv,
  lower = sf_sum$lower,
  upper = sf_sum$upper,
  strata = sub("^strategy_label=", "", sf_sum$strata)
) %>%
  mutate(
    cum_mortality = 1 - survival,
    lower_cm = 1 - upper,
    upper_cm = 1 - lower
  )

p_curve <- ggplot(curve_df, aes(x = time, y = cum_mortality, color = strata, fill = strata)) +
  geom_ribbon(aes(ymin = lower_cm, ymax = upper_cm), alpha = 0.14, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(0, 28, by = 7)) +
  scale_color_manual(values = c("Early CRRT strategy" = "#C0392B", "Deferred CRRT strategy" = "#2E86C1")) +
  scale_fill_manual(values = c("Early CRRT strategy" = "#C0392B", "Deferred CRRT strategy" = "#2E86C1")) +
  labs(
    title = "Strict first-KDIGO3 primary SOFA 8-14: IPCW-weighted cumulative mortality",
    subtitle = "Clone-censor-weighted target-trial strategies over 28 days",
    x = "Days since t0",
    y = "Cumulative mortality",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

png(file.path(fig_dir, "strictK3_weighted_cumulative_mortality_sofa8_14.png"), width = 2100, height = 1500, res = 300)
print(p_curve)
dev.off()
pdf(file.path(fig_dir, "strictK3_weighted_cumulative_mortality_sofa8_14.pdf"), width = 7, height = 5, onefile = FALSE)
print(p_curve)
dev.off()

# Forest plot ------------------------------------------------------------------
primary <- read_csv(file.path(results_dir, "strictK3_primary_mortality_results.csv"), show_col_types = FALSE) %>%
  transmute(group = "Primary / SOFA", analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)

sofa <- read_csv(file.path(results_dir, "strictK3_sofa_strata_results.csv"), show_col_types = FALSE) %>%
  transmute(group = "SOFA strata", analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)

rescue <- read_csv(file.path(results_dir, "strictK3_rescue_sensitivity_model_results.csv"), show_col_types = FALSE) %>%
  filter(grepl("SOFA 8-14$", analysis)) %>%
  transmute(group = "Rescue sensitivity", analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)

anchor <- read_csv(file.path(results_dir, "strictK3_anchor_sensitivity_comparison_results.csv"), show_col_types = FALSE) %>%
  filter(grepl("SOFA 8-14$", analysis)) %>%
  transmute(group = "Anchor sensitivity", analysis, n_patients, early_n, weighted_hr, ci_lower, ci_upper, p_value)

forest <- bind_rows(primary, sofa, rescue, anchor) %>%
  mutate(
    label = paste0(analysis, "  (n=", n_patients, ", early=", early_n, ")"),
    label = factor(label, levels = rev(label)),
    p_text = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  )

write_csv(forest, file.path(results_dir, "strictK3_forest_plot_source.csv"))

p_forest <- ggplot(forest, aes(x = weighted_hr, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.18, color = "#2C3E50") +
  geom_point(shape = 18, size = 3, color = "#2C3E50") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(0.75, 1, 1.25, 1.5, 2), limits = c(0.7, 2.05)) +
  labs(
    title = "Strict first-KDIGO3: weighted Cox sensitivity summary",
    subtitle = "Hazard ratio for 28-day mortality, Early vs Deferred CRRT strategy",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(strip.text.y = element_text(angle = 0, face = "bold"), axis.text.y = element_text(size = 7.5))

png(file.path(fig_dir, "strictK3_forest_plot_sensitivity_summary.png"), width = 3000, height = 2600, res = 300)
print(p_forest)
dev.off()
pdf(file.path(fig_dir, "strictK3_forest_plot_sensitivity_summary.pdf"), width = 10, height = 8.7, onefile = FALSE)
print(p_forest)
dev.off()

summary_report <- c(
  "# Strict first-KDIGO3 tables and figures",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Table 1",
  "",
  paste0("- `", file.path(results_dir, "strictK3_table1_sofa8_14.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_table1_sofa8_14.md"), "`"),
  paste0("- Variables with SMD >=0.10: ", sum(table1$smd_flag_ge_0_10, na.rm = TRUE), "."),
  "",
  "## Figures",
  "",
  paste0("- Love plot: `", file.path(fig_dir, "strictK3_love_plot_sofa8_14.png"), "`"),
  paste0("- Weighted cumulative mortality: `", file.path(fig_dir, "strictK3_weighted_cumulative_mortality_sofa8_14.png"), "`"),
  paste0("- Sensitivity forest plot: `", file.path(fig_dir, "strictK3_forest_plot_sensitivity_summary.png"), "`"),
  paste0("- SOFA band scan heatmap: `", file.path(fig_dir, "strictK3_sofa_band_scan_heatmap.png"), "`"),
  paste0("- RRT-free days figure: `", file.path(fig_dir, "strictK3_Figure_RRTfree_days.png"), "`")
)
writeLines(summary_report, file.path(results_dir, "strictK3_tables_figures_summary.md"))

cat("\nSaved Table 1 and figures under:\n")
cat("  ", results_dir, "\n", sep = "")
cat("  ", fig_dir, "\n", sep = "")
