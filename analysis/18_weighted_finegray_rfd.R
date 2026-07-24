# =============================================================================
# Step 18: IPCW-Weighted Fine-Gray Models + RRT-Free Days Alive at Day 28
#
# Analyses:
#   1. Weighted Fine-Gray subdistribution hazard models
#      - Event of interest: 28-day renal recovery
#      - Competing event:   death before day 28
#      - Reciprocal analysis: 28-day mortality with renal recovery as competing event
#      - Weighting: stabilized IPCW weights from the CCW pipeline
#      - Populations: full cohort + equipoise subgroup (SOFA 8-14 at t0)
#
#   2. RRT-free days alive (RFD) as a continuous outcome
#      - Death before day 28 => RFD = 0
#      - Survived to day 28 => RFD = 28 - rrt_hours/24
#      - Primary inference: weighted linear model with IPCW weights
#      - Wilcoxon rank-sum test retained as descriptive support only
#      - Median [IQR] per observed treatment group
#      - Figure saved to analysis/figures/Figure3_RRTfree_days.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(tibble)
})

set.seed(20260403)

# ── Paths ────────────────────────────────────────────────────────────────────
base_dir   <- "analysis"
data_dir   <- file.path(base_dir, "data")
fig_dir    <- file.path(base_dir, "figures")

ccw_path        <- file.path(data_dir, "ccw_final.rds")
cohort_path     <- file.path(data_dir, "cohort_base.rds")
tv_hourly_path  <- file.path(data_dir, "tv_hourly.rds")
tv_hour23_path  <- file.path(data_dir, "tv_hour23.rds")
master_csv_path <- file.path(data_dir, "analysis_master.csv")

fg_rds_path     <- file.path(data_dir, "weighted_finegray_results.rds")
fg_csv_path     <- file.path(data_dir, "weighted_finegray_results.csv")
rfd_rds_path    <- file.path(data_dir, "rfd_results.rds")
rfd_csv_path    <- file.path(data_dir, "rfd_results.csv")
all_rds_path    <- file.path(data_dir, "weighted_finegray_rfd_results.rds")
figure_png_path <- file.path(fig_dir, "Figure3_RRTfree_days.png")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper functions ─────────────────────────────────────────────────────────
fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(
    x %in% c("t", "TRUE", "true", "1"),
    TRUE,
    ifelse(x %in% c("f", "FALSE", "false", "0"), FALSE, NA)
  )
}

fmt_p_chr <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

weighted_quantile <- function(x, w, probs) {
  keep <- !(is.na(x) | is.na(w))
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0 || sum(w) <= 0) {
    return(rep(NA_real_, length(probs)))
  }
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
  # use "se(coef)" if "robust se" not present (non-cluster fit)
  se_col <- if ("robust se" %in% colnames(s$coefficients)) "robust se" else "se(coef)"
  se   <- s$coefficients[coef_name, se_col]
  pval <- s$coefficients[coef_name, "Pr(>|z|)"]

  tibble(
    Endpoint         = endpoint,
    Analysis         = label,
    N_patients       = n_patients,
    N_event          = n_event,
    N_competing      = n_competing,
    HR               = round(exp(beta), 3),
    HR_lo            = round(exp(beta - 1.96 * se), 3),
    HR_hi            = round(exp(beta + 1.96 * se), 3),
    p_value          = round(pval, 4),
    p_value_display  = fmt_p_chr(pval)
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
    n_event <- sum(dat$fg_status_recovery == "renal_recovery", na.rm = TRUE)
    n_comp  <- sum(dat$fg_status_recovery == "death", na.rm = TRUE)
  } else {
    status_factor <- factor(
      dat$fg_status_mortality,
      levels = c("censor", "death", "renal_recovery"),
      ordered = FALSE
    )
    etype_label <- "death"
    endpoint_label <- "Mortality (renal recovery competing)"
    n_event <- sum(dat$fg_status_mortality == "death", na.rm = TRUE)
    n_comp  <- sum(dat$fg_status_mortality == "renal_recovery", na.rm = TRUE)
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

  # finegray only preserves formula variables; extra columns (incl. stay_id)
  # are dropped, so we use robust=TRUE without cluster for sandwich SEs.
  weight_col <- grep("weight", names(fg_dat), ignore.case = TRUE, value = TRUE)
  if (length(weight_col) == 0) weight_col <- character(0)

  fg_dat <- fg_dat %>%
    mutate(fg_case_w = if (length(weight_col) > 0)
                         fgwt * .data[[weight_col[1]]]
                       else fgwt)

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
  dat <- dat %>% filter(!is.na(rfd_28d))

  # Primary inference for RFD uses the CCW clone-level dataset with stabilized
  # IPCW weights from ccw_final.rds.
  fit <- lm(rfd_28d ~ trt, data = dat, weights = ipcw_w_trim)
  vc  <- vcovCL(fit, cluster = dat$stay_id, type = "HC3")
  ct  <- coeftest(fit, vcov. = vc)

  beta <- unname(ct["trt", "Estimate"])
  se   <- unname(ct["trt", "Std. Error"])
  pval <- unname(ct["trt", "Pr(>|t|)"])

  early_obs <- dat %>% filter(early_crrt == TRUE)
  def_obs   <- dat %>% filter(early_crrt == FALSE)

  # Wilcoxon test uses observed treatment group (early_crrt), not clone-level
  # strategy — descriptive only, interpreted as observed-group comparison
  wilcox_p <- tryCatch(
    wilcox.test(rfd_28d ~ early_crrt, data = dat)$p.value,
    error = function(e) NA_real_
  )

  early_wq <- weighted_quantile(dat$rfd_28d[dat$trt == 1], dat$ipcw_w_trim[dat$trt == 1], c(0.25, 0.5, 0.75))
  def_wq   <- weighted_quantile(dat$rfd_28d[dat$trt == 0], dat$ipcw_w_trim[dat$trt == 0], c(0.25, 0.5, 0.75))

  tibble(
    Analysis                 = label,
    N_patients               = nrow(dat),
    Early_observed_median_IQR = fmt_miq(early_obs$rfd_28d),
    Deferred_observed_median_IQR = fmt_miq(def_obs$rfd_28d),
    Early_weighted_median_IQR = sprintf("%.1f [%.1f, %.1f]", early_wq[2], early_wq[1], early_wq[3]),
    Deferred_weighted_median_IQR = sprintf("%.1f [%.1f, %.1f]", def_wq[2], def_wq[1], def_wq[3]),
    Weighted_mean_diff       = round(beta, 3),
    Diff_lo                  = round(beta - 1.96 * se, 3),
    Diff_hi                  = round(beta + 1.96 * se, 3),
    Weighted_p_value         = round(pval, 4),
    Weighted_p_display       = fmt_p_chr(pval),
    Wilcoxon_p_observed_groups = round(wilcox_p, 4),
    Wilcoxon_p_observed_display = fmt_p_chr(wilcox_p)
  )
}

# ── Load data ────────────────────────────────────────────────────────────────
ccw <- readRDS(ccw_path) %>%
  distinct(stay_id, strategy, .keep_all = TRUE)

cohort_base <- readRDS(cohort_path) %>%
  distinct(stay_id, .keep_all = TRUE)

tv_hourly <- readRDS(tv_hourly_path)
tv_hour23 <- readRDS(tv_hour23_path)

cat("CCW rows after de-duplication:", nrow(ccw), "\n")
cat("Unique stays:", length(unique(ccw$stay_id)), "\n")
cat("tv_hourly hours:", paste(range(tv_hourly$hour_index, na.rm = TRUE), collapse = " to "), "\n")
cat("tv_hour23 hour index:", paste(unique(tv_hour23$hour_index), collapse = ", "), "\n")

# Secondary outcomes were extracted upstream in analysis/06_secondary_outcomes.sql.
# The hourly RDS files cover the 24-hour grace period for IPCW modeling only and
# do not contain 28-day CRRT exposure, so the authoritative 28-day renal recovery
# and RFD fields are taken from analysis_master.csv when available.
secondary <- read.csv(master_csv_path, stringsAsFactors = FALSE) %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    early_crrt         = fix_bool(early_crrt),
    died_28d           = fix_bool(died_28d),
    renal_recovery_28d = fix_bool(renal_recovery_28d)
  )

# ── Build patient-level analysis dataset ─────────────────────────────────────
# ccw already contains: early_crrt, died_28d, follow_time_28d, days_to_death,
# sofa_t0, equipoise. Only pull the secondary-outcome columns that are absent.
analysis_df <- ccw %>%
  left_join(
    secondary %>%
      select(stay_id, renal_recovery_28d, rrt_free_days_28d, crrt_total_days),
    by = "stay_id"
  ) %>%
  mutate(
    trt = as.integer(strategy == "A_early")
  )

# Derive RFD on an hours scale: RFD = max(0, 28 - rrt_hours/24).
# The upstream secondary-outcomes table stores `crrt_total_days` rather than a
# direct cumulative-hour field, so when available we convert it to an hour-scale
# equivalent (`rrt_hours = crrt_total_days * 24`) before applying the formula.
# This preserves the requested continuous form while keeping consistency with the
# extracted 28-day CRRT exposure table used elsewhere in the project.
analysis_df <- analysis_df %>%
  mutate(
    rrt_hours_28d = case_when(
      !is.na(crrt_total_days)   ~ crrt_total_days * 24,
      !is.na(rrt_free_days_28d) ~ pmax(0, 28 - rrt_free_days_28d) * 24,
      TRUE                      ~ NA_real_
    ),
    rfd_28d = case_when(
      !is.na(died_28d) & died_28d ~ 0,
      !is.na(rrt_hours_28d)       ~ pmax(0, 28 - rrt_hours_28d / 24),
      TRUE                        ~ NA_real_
    ),
    fg_status_recovery = case_when(
      renal_recovery_28d == TRUE ~ "renal_recovery",
      died_28d == TRUE           ~ "death",
      TRUE                       ~ "censor"
    ),
    fg_status_mortality = case_when(
      died_28d == TRUE           ~ "death",
      renal_recovery_28d == TRUE ~ "renal_recovery",
      TRUE                       ~ "censor"
    )
  )

cat("\nMerged analysis dataset:\n")
cat("  Rows:", nrow(analysis_df), "\n")
cat("  Unique stays:", length(unique(analysis_df$stay_id)), "\n")
cat("  Equipoise rows:", sum(analysis_df$equipoise, na.rm = TRUE), "\n")
cat("  Missing renal recovery:", sum(is.na(analysis_df$renal_recovery_28d)), "\n")
cat("  Missing RFD:", sum(is.na(analysis_df$rfd_28d)), "\n")

# ── Fine-Gray analyses ───────────────────────────────────────────────────────
full_df <- analysis_df
eq_df   <- analysis_df %>% filter(equipoise)

fg_results <- bind_rows(
  run_weighted_finegray(full_df, "Full cohort", endpoint = "renal_recovery"),
  run_weighted_finegray(eq_df,   "Equipoise SOFA 8-14", endpoint = "renal_recovery"),
  run_weighted_finegray(full_df, "Full cohort", endpoint = "mortality"),
  run_weighted_finegray(eq_df,   "Equipoise SOFA 8-14", endpoint = "mortality")
)

cat("\n=== IPCW-Weighted Fine-Gray Results ===\n")
print(fg_results, n = nrow(fg_results))

# ── RFD analyses ─────────────────────────────────────────────────────────────
rfd_desc_full <- analysis_df %>%
  distinct(stay_id, .keep_all = TRUE)
rfd_desc_eq <- rfd_desc_full %>%
  filter(equipoise)

rfd_weight_df <- analysis_df %>%
  filter(!is.na(rfd_28d), !is.na(ipcw_w_trim))
rfd_weight_eq <- rfd_weight_df %>%
  filter(equipoise)

rfd_results <- bind_rows(
  run_rfd_analysis(rfd_weight_df, "Full cohort"),
  run_rfd_analysis(rfd_weight_eq, "Equipoise SOFA 8-14")
)

cat("\n=== RFD Results ===\n")
print(rfd_results, n = nrow(rfd_results))

# ── Figure: stacked RFD distribution in the full cohort ─────────────────────
plot_df <- rfd_desc_full %>%
  mutate(
    group = factor(
      ifelse(early_crrt, "Early CRRT", "Deferred CRRT"),
      levels = c("Early CRRT", "Deferred CRRT")
    ),
    rfd_cat = factor(
      case_when(
        rfd_28d == 0  ~ "0 days",
        rfd_28d == 28 ~ "28 days",
        TRUE          ~ "1-27 days"
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
  "Primary IPCW-weighted linear model: mean diff %.2f days (95%% CI %.2f to %.2f), p %s; descriptive Wilcoxon p %s",
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
    values = c(
      "0 days" = "#C0392B",
      "1-27 days" = "#F0B27A",
      "28 days" = "#2E86C1"
    )
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
  annotate("segment", x = 1, xend = 2, y = 104, yend = 104, linewidth = 0.7) +
  annotate("segment", x = 1, xend = 1, y = 104, yend = 101, linewidth = 0.7) +
  annotate("segment", x = 2, xend = 2, y = 104, yend = 101, linewidth = 0.7) +
  annotate(
    "text",
    x = 1.5, y = 108,
    label = paste0("Wilcoxon p ", full_rfd_model$Wilcoxon_p_observed_display),
    size = 3.3
  ) +
  labs(
    title = "Figure 2. RRT-Free Days Alive at Day 28",
    subtitle = subtitle_text,
    x = NULL,
    y = "Proportion of patients (%)",
    fill = "RFD category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 8.5, color = "grey35", hjust = 0, lineheight = 1.2),
    axis.text.x = element_text(size = 9.5, color = "black", lineheight = 1.2),
    axis.text.y = element_text(size = 9.5),
    axis.title.y = element_text(size = 10.5, face = "bold", margin = margin(r = 8)),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 10)
  )

png(figure_png_path, width = 2250, height = 1650, res = 300)
print(p_rfd)
dev.off()

figure_pdf_path <- sub("\\.png$", ".pdf", figure_png_path)
pdf(figure_pdf_path, width = 7.5, height = 5.5, onefile = FALSE)
print(p_rfd)
dev.off()
cat("\nSaved figure:", figure_png_path, "\n")
cat("Saved figure:", figure_pdf_path, "\n")

# ── Save outputs ─────────────────────────────────────────────────────────────
saveRDS(fg_results, fg_rds_path)
write.csv(fg_results, fg_csv_path, row.names = FALSE)

saveRDS(rfd_results, rfd_rds_path)
write.csv(rfd_results, rfd_csv_path, row.names = FALSE)

saveRDS(
  list(
    finegray_results = fg_results,
    rfd_results = rfd_results
  ),
  all_rds_path
)

cat("\nSaved outputs:\n")
cat("  ", fg_rds_path, "\n", sep = "")
cat("  ", fg_csv_path, "\n", sep = "")
cat("  ", rfd_rds_path, "\n", sep = "")
cat("  ", rfd_csv_path, "\n", sep = "")
cat("  ", all_rds_path, "\n", sep = "")
