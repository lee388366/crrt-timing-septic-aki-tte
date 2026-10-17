# =============================================================================
# Strict clean first-KDIGO3: IPCW-weighted Fine-Gray + RRT-free days
#
# This script is the strict-KDIGO3 counterpart of 18_weighted_finegray_rfd.R.
# It deliberately reads secondary outcomes from
# data_extract_crrt_strictk3_v2.secondary_outcomes, not legacy analysis_master.csv.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(readr)
  library(tibble)
})

set.seed(20260715)

base_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(Sys.getenv("CRRT_SECONDARY_OUT"),"run13")
fig_dir <- file.path(results_dir,"figs")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ccw_path <- Sys.getenv("CRRT_V2_CCW")
secondary_csv_path <- file.path(results_dir,"strictK3_v2_secondary_outcomes.csv")
fg_csv_path <- file.path(results_dir, "strictK3_v2_weighted_finegray_results.csv")
fg_rds_path <- file.path(results_dir, "strictK3_v2_weighted_finegray_results.rds")
rfd_csv_path <- file.path(results_dir, "strictK3_v2_rfd_results.csv")
rfd_rds_path <- file.path(results_dir, "strictK3_v2_rfd_results.rds")
all_rds_path <- file.path(results_dir, "strictK3_v2_weighted_finegray_rfd_results.rds")
report_path <- file.path(results_dir, "strictK3_v2_weighted_finegray_rfd_report.md")
figure_png_path <- file.path(fig_dir, "strictK3_v2_Figure_RRTfree_days.png")

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(
    x %in% c("t", "TRUE", "true", "1", 1),
    TRUE,
    ifelse(x %in% c("f", "FALSE", "false", "0", 0), FALSE, NA)
  )
}

fmt_p_chr <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

weighted_quantile <- function(x, w, probs) {
  keep <- !(is.na(x) | is.na(w))
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0 || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

fmt_miq <- function(x) {
  sprintf(
    "%.1f [%.1f, %.1f]",
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

extract_fg_res <- function(fit, label, endpoint, n_patients, n_event, n_competing) {
  s <- summary(fit)
  coef_name <- rownames(s$coefficients)[1]
  beta <- s$coefficients[coef_name, "coef"]
  se_col <- if ("robust se" %in% colnames(s$coefficients)) "robust se" else "se(coef)"
  se <- s$coefficients[coef_name, se_col]
  pval <- s$coefficients[coef_name, "Pr(>|z|)"]

  tibble(
    Endpoint = endpoint,
    Analysis = label,
    N_patients = n_patients,
    N_event = n_event,
    N_competing = n_competing,
    HR = round(exp(beta), 3),
    HR_lo = round(exp(beta - 1.96 * se), 3),
    HR_hi = round(exp(beta + 1.96 * se), 3),
    p_value = round(pval, 4),
    p_value_display = fmt_p_chr(pval)
  )
}

run_weighted_finegray <- function(dat, label, endpoint = c("renal_recovery", "mortality")) {
  endpoint <- match.arg(endpoint)

  if (endpoint == "renal_recovery") {
    status_factor <- factor(
      dat$fg_status_recovery,
      levels = c("censor", "renal_recovery", "death"),
      ordered = FALSE
    )
    etype_label <- "renal_recovery"
    endpoint_label <- "Renal recovery (death competing)"
    n_event <- length(unique(dat$stay_id[dat$fg_status_recovery == "renal_recovery"]))
    n_comp <- length(unique(dat$stay_id[dat$fg_status_recovery == "death"]))
  } else {
    status_factor <- factor(
      dat$fg_status_mortality,
      levels = c("censor", "death", "renal_recovery"),
      ordered = FALSE
    )
    etype_label <- "death"
    endpoint_label <- "Mortality (renal recovery competing)"
    n_event <- length(unique(dat$stay_id[dat$fg_status_mortality == "death"]))
    n_comp <- length(unique(dat$stay_id[dat$fg_status_mortality == "renal_recovery"]))
  }

  fg_in <- dat %>%
    transmute(
      stay_id = stay_id,
      trt = as.integer(strategy == "A_early"),
      time = pmin(follow_time_28d, 28),
      status = status_factor,
      ipcw_w_trim = ipcw_w_trim
    ) %>%
    filter(!is.na(time), !is.na(status), !is.na(ipcw_w_trim))

  fg_dat <- finegray(
    Surv(time, status) ~ trt,
    data = fg_in,
    weights = ipcw_w_trim,
    etype = etype_label
  )

  weight_col <- grep("weight", names(fg_dat), ignore.case = TRUE, value = TRUE)
  fg_dat <- fg_dat %>%
    mutate(fg_case_w = if (length(weight_col) > 0) fgwt * .data[[weight_col[1]]] else fgwt)

  fg_fit <- coxph(
    Surv(fgstart, fgstop, fgstatus) ~ trt,
    data = fg_dat,
    weights = fg_case_w,
    ties = "breslow"
  )

  extract_fg_res(
    fit = fg_fit,
    label = label,
    endpoint = endpoint_label,
    n_patients = length(unique(fg_in$stay_id)),
    n_event = n_event,
    n_competing = n_comp
  )
}

run_rfd_analysis <- function(dat, label) {
  dat <- dat %>% filter(!is.na(rfd_28d), !is.na(ipcw_w_trim))

  fit <- lm(rfd_28d ~ trt, data = dat, weights = ipcw_w_trim)
  vc <- vcovCL(fit, cluster = dat$stay_id, type = "HC3")
  ct <- coeftest(fit, vcov. = vc)

  beta <- unname(ct["trt", "Estimate"])
  se <- unname(ct["trt", "Std. Error"])
  pval <- unname(ct["trt", "Pr(>|t|)"])

  patient_level <- dat %>% distinct(stay_id, .keep_all = TRUE)
  early_obs <- patient_level %>% filter(early_crrt == TRUE)
  def_obs <- patient_level %>% filter(early_crrt == FALSE)

  wilcox_p <- tryCatch(
    wilcox.test(rfd_28d ~ early_crrt, data = patient_level)$p.value,
    error = function(e) NA_real_
  )

  early_wq <- weighted_quantile(dat$rfd_28d[dat$trt == 1], dat$ipcw_w_trim[dat$trt == 1], c(0.25, 0.5, 0.75))
  def_wq <- weighted_quantile(dat$rfd_28d[dat$trt == 0], dat$ipcw_w_trim[dat$trt == 0], c(0.25, 0.5, 0.75))

  tibble(
    Analysis = label,
    N_patients = length(unique(dat$stay_id)),
    Early_observed_n = nrow(early_obs),
    Deferred_observed_n = nrow(def_obs),
    Early_observed_median_IQR = fmt_miq(early_obs$rfd_28d),
    Deferred_observed_median_IQR = fmt_miq(def_obs$rfd_28d),
    Early_weighted_median_IQR = sprintf("%.1f [%.1f, %.1f]", early_wq[2], early_wq[1], early_wq[3]),
    Deferred_weighted_median_IQR = sprintf("%.1f [%.1f, %.1f]", def_wq[2], def_wq[1], def_wq[3]),
    Weighted_mean_diff = round(beta, 3),
    Diff_lo = round(beta - 1.96 * se, 3),
    Diff_hi = round(beta + 1.96 * se, 3),
    Weighted_p_value = round(pval, 4),
    Weighted_p_display = fmt_p_chr(pval),
    Wilcoxon_p_observed_groups = round(wilcox_p, 4),
    Wilcoxon_p_observed_display = fmt_p_chr(wilcox_p)
  )
}

# Offline adapter: connection block not evaluated; original local snapshots only.

secondary <- readr::read_csv(Sys.getenv("CRRT_SECONDARY_CSV"), show_col_types=FALSE) %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    early_crrt = fix_bool(early_crrt),
    died_28d = fix_bool(died_28d),
    renal_recovery_28d = fix_bool(renal_recovery_28d),
    crrt_total_days = as.numeric(crrt_total_days),
    rrt_free_days_28d = as.numeric(rrt_free_days_28d)
  )

write_csv(secondary, secondary_csv_path)

ccw <- readRDS(ccw_path) %>%
  distinct(stay_id, strategy, .keep_all = TRUE)

analysis_df <- ccw %>%
  left_join(
    secondary %>%
      select(stay_id, renal_recovery_28d, rrt_free_days_28d, crrt_total_days),
    by = "stay_id"
  ) %>%
  # ===== CH-01 最小改动：仅重定义 died_28d / follow_time_28d（deathtime 或 patients.dod，午夜日界）=====
  left_join(
    readRDS(Sys.getenv("CRRT_DOD_RDS")),
    by = "subject_id"
  ) %>%
  mutate(
    .dc = coalesce(days_to_death, as.numeric(dod - as.Date(t0))),
    died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
    follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
  ) %>%
  select(-.dc, -dod) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    rrt_hours_28d = case_when(
      !is.na(crrt_total_days) ~ crrt_total_days * 24,
      !is.na(rrt_free_days_28d) ~ pmax(0, 28 - rrt_free_days_28d) * 24,
      TRUE ~ NA_real_
    ),
    rfd_28d = case_when(
      !is.na(died_28d) & died_28d ~ 0,
      !is.na(rrt_hours_28d) ~ pmax(0, 28 - rrt_hours_28d / 24),
      TRUE ~ NA_real_
    ),
    fg_status_recovery = case_when(
      renal_recovery_28d == TRUE ~ "renal_recovery",
      died_28d == TRUE ~ "death",
      TRUE ~ "censor"
    ),
    fg_status_mortality = case_when(
      died_28d == TRUE ~ "death",
      renal_recovery_28d == TRUE ~ "renal_recovery",
      TRUE ~ "censor"
    )
  )

cat("Merged strictK3_v2 Fine-Gray/RFD dataset\n")
cat("  Clone rows:", nrow(analysis_df), "\n")
cat("  Unique stays:", length(unique(analysis_df$stay_id)), "\n")
cat("  Equipoise clone rows:", sum(analysis_df$equipoise, na.rm = TRUE), "\n")
cat("  Missing renal recovery:", sum(is.na(analysis_df$renal_recovery_28d)), "\n")
cat("  Missing RFD:", sum(is.na(analysis_df$rfd_28d)), "\n")
cat("  [CH-01 check] unique stays died_28d:", length(unique(analysis_df$stay_id[analysis_df$died_28d %in% TRUE])), "\n")

full_df <- analysis_df
eq_df <- analysis_df %>% filter(equipoise)

fg_results <- bind_rows(
  run_weighted_finegray(full_df, "Full cohort", endpoint = "renal_recovery"),
  run_weighted_finegray(eq_df, "Equipoise SOFA 8-14", endpoint = "renal_recovery"),
  run_weighted_finegray(full_df, "Full cohort", endpoint = "mortality"),
  run_weighted_finegray(eq_df, "Equipoise SOFA 8-14", endpoint = "mortality")
)

rfd_results <- bind_rows(
  run_rfd_analysis(full_df, "Full cohort"),
  run_rfd_analysis(eq_df, "Equipoise SOFA 8-14")
)

write_csv(fg_results, fg_csv_path)
saveRDS(fg_results, fg_rds_path)
write_csv(rfd_results, rfd_csv_path)
saveRDS(rfd_results, rfd_rds_path)
saveRDS(list(finegray_results = fg_results, rfd_results = rfd_results), all_rds_path)

plot_df <- analysis_df %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    group = factor(
      ifelse(early_crrt, "Early CRRT", "Deferred CRRT"),
      levels = c("Early CRRT", "Deferred CRRT")
    ),
    rfd_cat = factor(
      case_when(
        rfd_28d == 0 ~ "0 days",
        rfd_28d == 28 ~ "28 days",
        TRUE ~ "1-27 days"
      ),
      levels = c("0 days", "1-27 days", "28 days")
    )
  )

bar_df <- plot_df %>%
  count(group, rfd_cat, name = "n") %>%
  group_by(group) %>%
  mutate(
    total = sum(n),
    pct = 100 * n / total,
    label = ifelse(pct >= 4, sprintf("%.1f%%", pct), "")
  ) %>%
  ungroup()

full_rfd_model <- rfd_results %>% filter(Analysis == "Full cohort")
subtitle_text <- sprintf(
  "IPCW-weighted linear model: mean diff %.2f days (95%% CI %.2f to %.2f), p %s; descriptive Wilcoxon p %s",
  full_rfd_model$Weighted_mean_diff,
  full_rfd_model$Diff_lo,
  full_rfd_model$Diff_hi,
  full_rfd_model$Weighted_p_display,
  full_rfd_model$Wilcoxon_p_observed_display
)

p_rfd <- ggplot(bar_df, aes(x = group, y = pct, fill = rfd_cat)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.4,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("0 days" = "#C0392B", "1-27 days" = "#F0B27A", "28 days" = "#2E86C1")
  ) +
  scale_x_discrete(
    labels = c(
      "Early CRRT" = sprintf(
        "Early CRRT\n(n = %d, median [IQR] %s)",
        sum(plot_df$early_crrt, na.rm = TRUE),
        fmt_miq(plot_df$rfd_28d[plot_df$early_crrt])
      ),
      "Deferred CRRT" = sprintf(
        "Deferred CRRT\n(n = %d, median [IQR] %s)",
        sum(!plot_df$early_crrt, na.rm = TRUE),
        fmt_miq(plot_df$rfd_28d[!plot_df$early_crrt])
      )
    )
  ) +
  scale_y_continuous(
    limits = c(0, 112),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    title = "Strict clean first-KDIGO3: RRT-Free Days Alive at Day 28",
    subtitle = subtitle_text,
    x = NULL,
    y = "Proportion of patients (%)",
    fill = "RFD category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.2),
    axis.text.x = element_text(size = 9.5, color = "black", lineheight = 1.2),
    axis.title.y = element_text(size = 10.5, face = "bold", margin = margin(r = 8)),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35)
  )

png(figure_png_path, width = 2250, height = 1650, res = 300)
print(p_rfd)
dev.off()
figure_pdf_path <- sub("\\.png$", ".pdf", figure_png_path)
pdf(figure_pdf_path, width = 7.5, height = 5.5, onefile = FALSE)
print(p_rfd)
dev.off()

report <- c(
  "# Strict clean first-KDIGO3 weighted Fine-Gray and RRT-free-days results",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Data source",
  "",
  "- Cohort/weights: `analysis/strictK3_primary_v2/data/ccw_final.rds`.",
  "- Secondary outcomes: `data_extract_crrt_strictk3_v2.secondary_outcomes_ch01` (CH-01口径), generated by `analysis/strictK3_v2_secondary_outcomes.sql`.",
  "- Legacy `analysis_master.csv` was not used.",
  "",
  "## Dataset checks",
  "",
  paste0("- Clone rows: ", nrow(analysis_df), "."),
  paste0("- Unique stays: ", length(unique(analysis_df$stay_id)), "."),
  paste0("- Equipoise clone rows: ", sum(analysis_df$equipoise, na.rm = TRUE), "."),
  paste0("- Missing renal recovery clone rows: ", sum(is.na(analysis_df$renal_recovery_28d)), "."),
  paste0("- Missing RFD clone rows: ", sum(is.na(analysis_df$rfd_28d)), "."),
  "",
  "## IPCW-weighted Fine-Gray results",
  "",
  paste(capture.output(print(fg_results, n = nrow(fg_results))), collapse = "\n"),
  "",
  "## RRT-free days alive at day 28",
  "",
  paste(capture.output(print(rfd_results, n = nrow(rfd_results))), collapse = "\n"),
  "",
  "## Outputs",
  "",
  paste0("- `", fg_csv_path, "`"),
  paste0("- `", rfd_csv_path, "`"),
  paste0("- `", figure_png_path, "`"),
  paste0("- `", figure_pdf_path, "`")
)

writeLines(report, report_path)

cat("\n=== Fine-Gray ===\n")
print(fg_results, n = nrow(fg_results))
cat("\n=== RFD ===\n")
print(rfd_results, n = nrow(rfd_results))
cat("\nSaved report:", report_path, "\n")