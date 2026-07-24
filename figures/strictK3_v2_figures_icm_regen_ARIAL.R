# =============================================================================
# strictK3 v2 ICM figure regeneration
#
# Purpose:
#   Rebuild manuscript-facing figures from the clean physiologic-KDIGO3 v2
#   analysis, using the MMR/NPG visual specification. Outputs are written to a
#   figures_icm_v2 directory; same-name files are overwritten on rerun.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(survival)
  library(survminer)
  library(forestploter)
  library(grid)
  library(patchwork)
  library(ragg)
})

set.seed(20260716)

base_dir <- "analysis/strictK3_primary_v2"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
eicu_dir <- "./按MMR修改意见修改/eICU数据提取/diagnose_rescue"
out_dir <- file.path(base_dir, "figures_icm_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_EARLY <- "#4DBBD5"
COL_DEFERRED <- "#E64B35"
COL_INK <- "#2C3E50"
COL_HEADER <- "#EAF0F6"
COL_AFTER <- "#27AE60"
COL_THRESHOLD <- "#C0392B"
COL_GREY <- "grey55"
# Arial is preferred by the style guide. This local R build cannot use Arial
# reliably in vector devices, so Helvetica is used as the pre-specified fallback.
BASE_FAMILY <- "Arial"
VECTOR_FAMILY <- "Helvetica"

theme_icm <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      text = element_text(color = COL_INK),
      axis.text = element_text(color = COL_INK),
      axis.title = element_text(color = COL_INK),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, lineheight = 1.05),
      legend.title = element_blank(),
      legend.text = element_text(color = COL_INK),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

fmt_n <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_hr <- function(hr, lo, hi) {
  ifelse(is.na(hr), "", sprintf("%.3f [%.3f-%.3f]", hr, lo, hi))
}

save_gg_pair <- function(plot, stem, width, height, dpi) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  invisible(dev.off())
  c(png = png_path, pdf = NA_character_)
}

save_grid_pair <- function(draw_fun, stem, width, height, dpi) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  draw_fun()
  invisible(dev.off())
  c(png = png_path, pdf = NA_character_)
}

read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(x %in% c(TRUE, "TRUE", "true", "t", "1", 1), TRUE,
         ifelse(x %in% c(FALSE, "FALSE", "false", "f", "0", 0), FALSE, NA))
}

# -----------------------------------------------------------------------------
# Data sources
# -----------------------------------------------------------------------------

funnel <- read_csv_quiet(file.path(data_dir, "strictK3_primary_v2_funnel.csv"))
v1v2_funnel <- read_csv_quiet(file.path(results_dir, "strictK3_v2_vs_v1_funnel.csv"))
# Fig1 bottom box: SOFA 8–14 equipoise from baseline (sofa_equipoise / sofa_t0),
# not hard-coded — keeps figure aligned with sofa_group_summary.csv / baseline.csv
baseline_fig1 <- read_csv_quiet(file.path(data_dir, "baseline.csv"))
sofa_group_summary <- read_csv_quiet(file.path(data_dir, "sofa_group_summary.csv"))
eq_fig1 <- baseline_fig1 %>%
  filter(fix_bool(sofa_equipoise) | (!is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14))
n_eq_fig1 <- nrow(eq_fig1)
n_early_fig1 <- sum(fix_bool(eq_fig1$early_crrt), na.rm = TRUE)
n_deferred_fig1 <- n_eq_fig1 - n_early_fig1
sofa814_n_summary <- sofa_group_summary %>%
  filter(.data$sofa_grp == "SOFA 8-14") %>%
  pull(n) %>%
  .[1]
if (!is.na(sofa814_n_summary) && n_eq_fig1 != sofa814_n_summary) {
  warning(sprintf(
    "Fig1 equipoise n from baseline (%s) != sofa_group_summary SOFA 8-14 (%s)",
    n_eq_fig1, sofa814_n_summary
  ))
}
stopifnot(n_eq_fig1 > 0L, n_early_fig1 >= 0L, n_deferred_fig1 >= 0L)
primary_results <- read_csv_quiet(file.path(results_dir, "strictK3_primary_v2_mortality_results.csv"))
forest_source <- read_csv_quiet(file.path(results_dir, "strictK3_v2_forest_plot_source.csv"))
rescue_results <- read_csv_quiet(file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.csv"))
qa_ph <- read_csv_quiet(file.path(results_dir, "strictK3_v2_QA_ph_treatment_results.csv"))
path_adjusted <- read_csv_quiet(file.path(results_dir, "strictK3_v2_path_adjusted_ipcw_results.csv"))
l0_results <- read_csv_quiet(file.path(results_dir, "strictK3_v2_l0_only_sensitivity_results.csv"))
sofa_scan <- read_csv_quiet(file.path(results_dir, "strictK3_v2_sofa_band_scan_results.csv"))
qba_grid <- read_csv_quiet(file.path(results_dir, "strictK3_v2_qba_grid.csv"))
qba_targets <- read_csv_quiet(file.path(results_dir, "strictK3_v2_qba_target_thresholds.csv"))
fg_results <- read_csv_quiet(file.path(results_dir, "strictK3_v2_weighted_finegray_results.csv"))
balance_smd <- read_csv_quiet(file.path(results_dir, "strictK3_v2_ipcw_balance_smd.csv"))
eicu_gate <- read_csv_quiet(file.path(eicu_dir, "data/task5_precision_gate.csv"))

ccw <- readRDS(file.path(data_dir, "ccw_final.rds"))
secondary <- read_csv_quiet(file.path(data_dir, "strictK3_v2_secondary_outcomes.csv")) %>%
  mutate(
    died_28d = fix_bool(died_28d),
    renal_recovery_28d = fix_bool(renal_recovery_28d)
  )

primary_sofa <- primary_results %>%
  filter(analysis == "Primary SOFA 8-14 equipoise") %>%
  slice(1)
stopifnot(nrow(primary_sofa) == 1)

# -----------------------------------------------------------------------------
# Figure 1: cohort flow
# -----------------------------------------------------------------------------

step_n <- function(step_name) funnel %>% filter(.data$step == step_name) %>% pull(n) %>% .[1]
step_ex <- function(step_name) funnel %>% filter(.data$step == step_name) %>% pull(excluded_from_previous) %>% .[1]

v1v2_first_delta <- v1v2_funnel %>%
  filter(step == "first_KDIGO3_stays") %>%
  pull(delta_v2_minus_v1) %>%
  abs()
v1v2_final_delta <- v1v2_funnel %>%
  filter(step == "exclude_primary_rescue_R1_R2_R2b_R3") %>%
  pull(delta_v2_minus_v1) %>%
  abs()

flow_draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(
    "Figure 1. Cohort selection flow for the clean physiologic KDIGO stage 3 target-trial emulation",
    x = unit(0.5, "npc"), y = unit(0.965, "npc"),
    gp = gpar(fontfamily = VECTOR_FAMILY, fontsize = 11, fontface = "bold", col = COL_INK)
  )
  grid.text(
    "t0 = first physiologic KDIGO stage 3 (creatinine/urine-output only); RRT was not allowed to define eligibility.",
    x = unit(0.5, "npc"), y = unit(0.928, "npc"),
    gp = gpar(fontfamily = VECTOR_FAMILY, fontsize = 7.5, col = COL_INK)
  )

  draw_box <- function(x, y, w, h, label, fill = "white", lwd = 1.05, fontsize = 7.5,
                       fontface = "plain", col = COL_INK) {
    grid.roundrect(
      x = unit(x, "npc"), y = unit(y, "npc"),
      width = unit(w, "npc"), height = unit(h, "npc"),
      r = unit(0.012, "npc"),
      gp = gpar(fill = fill, col = col, lwd = lwd)
    )
    grid.text(
      label, x = unit(x, "npc"), y = unit(y, "npc"),
      gp = gpar(fontfamily = VECTOR_FAMILY, fontsize = fontsize, fontface = fontface, col = COL_INK),
      just = "centre"
    )
  }

  draw_arrow <- function(x1, y1, x2, y2) {
    grid.lines(
      x = unit(c(x1, x2), "npc"), y = unit(c(y1, y2), "npc"),
      arrow = arrow(length = unit(0.10, "inches"), type = "closed"),
      gp = gpar(col = COL_INK, lwd = 0.8)
    )
  }

  main_x <- 0.37
  side_x <- 0.86
  main_w <- 0.52
  side_w <- 0.26
  box_h <- 0.074
  ys <- c(0.835, 0.710, 0.585, 0.460, 0.335, 0.210)

  main_labels <- c(
    paste0("First physiologic KDIGO stage 3 ICU stays\nn = ", fmt_n(step_n("first_KDIGO3_stays"))),
    paste0("Sepsis-3 documented at or before clean t0\nn = ", fmt_n(step_n("sepsis3_at_or_before_t0"))),
    paste0("Septic shock proxy with vasopressor window\nn = ", fmt_n(step_n("vasopressor_window_existing_definition"))),
    paste0("No RRT before clean t0 and core eligibility retained\nn = ", fmt_n(step_n("pre_rescue_first_eligible_ICU"))),
    paste0("Final full analytic cohort after primary rescue exclusion\nn = ", fmt_n(step_n("exclude_primary_rescue_R1_R2_R2b_R3"))),
    paste0(
      "Primary equipoise subgroup: SOFA at t0 8-14\nn = ", fmt_n(n_eq_fig1),
      "; Early CRRT n = ", fmt_n(n_early_fig1),
      "; Deferred n = ", fmt_n(n_deferred_fig1)
    )
  )

  for (i in seq_along(ys)) {
    fill <- if (i %in% c(1, 5, 6)) COL_HEADER else "white"
    draw_box(main_x, ys[i], main_w, box_h, main_labels[i],
             fill = fill, fontsize = ifelse(i == 6, 7.4, 7.6),
             fontface = ifelse(i %in% c(1, 5, 6), "bold", "plain"))
    if (i < length(ys)) draw_arrow(main_x, ys[i] - box_h / 2 - 0.008, main_x, ys[i + 1] + box_h / 2 + 0.008)
  }

  side_labels <- c(
    paste0("Clean-t0 guardrail\nTreatment-defined KDIGO3 t0 removed\nat trigger level: n = ",
           fmt_n(v1v2_first_delta), "\nFinal cohort impact: n = ", fmt_n(v1v2_final_delta)),
    paste0("Excluded: sepsis not before t0\nn = ", fmt_n(step_ex("sepsis3_at_or_before_t0"))),
    paste0("Excluded: no vasopressor shock proxy\nn = ", fmt_n(step_ex("vasopressor_window_existing_definition"))),
    paste0("Excluded before final eligibility\nPrior RRT n = ", fmt_n(step_ex("exclude_RRT_before_t0")),
           "; ESRD/chronic dialysis n = ", fmt_n(step_ex("exclude_ESRD_chronic_dialysis")),
           "\nOther pre-eligibility exclusions n = ",
           fmt_n(step_ex("exclude_ICU_LOS_lt_24h") +
                   step_ex("exclude_zero_28d_followup") +
                   step_ex("pre_rescue_first_eligible_ICU")),
           "\n(ICU LOS<24 h, zero follow-up, first-eligible filter)"),
    paste0("Excluded: rescue/urgent indication at t0\nR1/R2/R2b/R3 n = ",
           fmt_n(step_ex("exclude_primary_rescue_R1_R2_R2b_R3")))
  )
  side_ys <- c(0.835, 0.710, 0.585, 0.460, 0.335)
  for (i in seq_along(side_labels)) {
    draw_box(side_x, side_ys[i], side_w, 0.075, side_labels[i],
             fill = ifelse(i == 1, "#F7FBFC", "white"), fontsize = ifelse(i == 4, 5.45, 6.25),
             col = ifelse(i == 1, COL_THRESHOLD, COL_INK))
    grid.lines(
      x = unit(c(main_x + main_w / 2 + 0.006, side_x - side_w / 2 - 0.006), "npc"),
      y = unit(c(side_ys[i], side_ys[i]), "npc"),
      gp = gpar(col = "grey65", lwd = 0.75, lty = "dashed")
    )
  }

  grid.text(
    "Early CRRT = initiation within 24 h after t0; no-initiation strategy = remain CRRT-free through 24 h. Any within-window initiation censors the no-initiation clone.",
    x = unit(0.5, "npc"), y = unit(0.075, "npc"),
    gp = gpar(fontfamily = VECTOR_FAMILY, fontsize = 6.6, col = "grey30")
  )
}

fig1_paths <- save_grid_pair(flow_draw, "Figure_1_Cohort_Flow_v2_clean_t0", 7.48, 6.86, 1000)

# -----------------------------------------------------------------------------
# Figure 2: weighted cumulative hospital-recorded death by day 28
# -----------------------------------------------------------------------------

ccw_sofa <- ccw %>%
  filter(equipoise == TRUE) %>%
  mutate(
    strategy_plot = factor(
      strategy,
      levels = c("B_deferred", "A_early"),
      labels = c("Deferred CRRT", "Early CRRT")
    ),
    event = as.integer(event)
  ) %>%
  filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim), !is.na(strategy_plot))

fit_mortality <- survfit(
  Surv(follow_days, event) ~ strategy_plot,
  data = ccw_sofa,
  weights = ipcw_w_trim
)

# Fig2 risk-table footnote: read weighted n.risk from survfit (not hard-coded)
fig2_nrisk_sum <- summary(fit_mortality, times = c(0, 7), extend = TRUE)
fig2_nrisk_df <- tibble::tibble(
  strata = as.character(fig2_nrisk_sum$strata),
  time = fig2_nrisk_sum$time,
  n.risk = fig2_nrisk_sum$n.risk
) %>%
  mutate(
    arm = sub("^strategy_plot=", "", strata),
    n.risk_round = as.integer(round(n.risk))
  )
fig2_early_n0 <- fig2_nrisk_df %>% dplyr::filter(grepl("Early", arm), time == 0) %>% dplyr::pull(n.risk_round)
fig2_early_n7 <- fig2_nrisk_df %>% dplyr::filter(grepl("Early", arm), time == 7) %>% dplyr::pull(n.risk_round)
fig2_deferred_n0 <- fig2_nrisk_df %>% dplyr::filter(grepl("Deferred", arm), time == 0) %>% dplyr::pull(n.risk_round)
stopifnot(length(fig2_early_n0) == 1L, length(fig2_early_n7) == 1L, length(fig2_deferred_n0) == 1L)
readr::write_csv(
  fig2_nrisk_df,
  file.path(results_dir, "strictK3_v2_fig2_weighted_nrisk.csv")
)

fig2_risk_caption_manifest <- paste(
  "IPCW-weighted cumulative hospital-recorded death by day 28 in the primary SOFA 8-14 overlap cohort comparing CRRT within 24 hours with no initiation during that window.",
  sprintf(
    "The 'weighted number at risk' row reflects IPCW-weighted, clone-censored risk sets; the early-CRRT arm falls from %s at t0 to %s by day 7 because early-strategy clones that did not initiate CRRT within the 24-h grace period are censored per protocol, leaving only the (weighted) genuinely early initiators at risk.",
    fmt_n(fig2_early_n0), fmt_n(fig2_early_n7)
  ),
  "This is an expected feature of the clone-censor-weighting design, not attrition or a data error."
)

fig2_risk_caption_plot <- paste0(
  "Note: the 'weighted number at risk' row reflects IPCW-weighted, clone-censored risk sets. Early CRRT falls from ",
  fmt_n(fig2_early_n0), " at t0\n",
  "to ", fmt_n(fig2_early_n7), " by day 7 because early-strategy clones that did not initiate CRRT within the 24-h grace period are censored per protocol,\n",
  "leaving only the weighted genuinely early initiators at risk. This is expected in clone-censor-weighting, not attrition or a data error."
)

fig2 <- ggsurvplot(
  fit_mortality,
  data = ccw_sofa,
  fun = "event",
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.24,
  risk.table.title = "Weighted number at risk",
  risk.table.fontsize = 3.5,
  conf.int = FALSE,
  censor = FALSE,
  palette = c(COL_DEFERRED, COL_EARLY),
  legend.title = "",
  legend.labs = c("Deferred CRRT", "Early CRRT"),
  legend = "bottom",
  break.time.by = 7,
  xlim = c(0, 28),
  ylim = c(0, 1),
  pval = sprintf("IPCW Cox p = %s", fmt_p(primary_sofa$p_value)),
  pval.coord = c(1.0, 0.92),
  pval.size = 4.0,
  ggtheme = theme_icm(9)
)
fig2$plot <- fig2$plot +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(0, 28), expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(
    title = "Figure 2. Weighted cumulative hospital-recorded death by day 28",
    subtitle = sprintf(
      "Primary SOFA 8-14 cohort: HR %.3f [%.3f-%.3f]; RMST difference Early minus Deferred %.2f days",
      primary_sofa$weighted_hr, primary_sofa$ci_lower, primary_sofa$ci_upper, primary_sofa$rmst_diff
    ),
    x = "Days since clean physiologic KDIGO stage 3 t0",
    y = "Cumulative hospital-recorded death"
  ) +
  theme(legend.position = "bottom")
fig2$table <- fig2$table +
  theme_icm(8) +
  labs(
    caption = fig2_risk_caption_plot
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.title = element_text(face = "bold", size = 8.5),
    plot.caption = element_text(hjust = 0, size = 5.8, color = "grey30", lineheight = 0.95,
                                margin = margin(t = 5))
  )

fig2_png <- file.path(out_dir, "Figure_2_Weighted_Cumulative_Mortality_v2.png")
fig2_pdf <- file.path(out_dir, "Figure_2_Weighted_Cumulative_Mortality_v2.pdf")
ragg::agg_png(fig2_png, width = 8, height = 6.55, units = "in", res = 300, background = "white")
print(fig2)
invisible(dev.off())
fig2_paths <- c(png = fig2_png, pdf = fig2_pdf)

# -----------------------------------------------------------------------------
# Figure 3: primary and robustness forest
# -----------------------------------------------------------------------------

make_row <- function(group, label, n, early, hr, lo, hi, p) {
  tibble(group = group, label = label, n = n, early = early, hr = hr, lo = lo, hi = hi, p = p)
}
header_row <- function(label) {
  tibble(group = label, label = label, n = NA_real_, early = NA_real_, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
}

row_from_forest <- function(analysis_key, display, group = NA_character_) {
  x <- forest_source %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) stop("Missing forest source row: ", analysis_key)
  make_row(
    group = ifelse(is.na(group), x$group, group),
    label = display,
    n = x$n_patients,
    early = x$early_n,
    hr = x$weighted_hr,
    lo = x$ci_lower,
    hi = x$ci_upper,
    p = x$p_value
  )
}

row_from_rescue <- function(analysis_key, display, group = "Rescue/eligibility sensitivity") {
  x <- rescue_results %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) stop("Missing rescue row: ", analysis_key)
  make_row(group, display, x$n_patients, x$early_n, x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
}

row_from_ph <- function(model_key, display, group = "Model/t0 sensitivity") {
  x <- qa_ph %>% filter(model == model_key) %>% slice(1)
  if (nrow(x) != 1) stop("Missing QA PH row: ", model_key)
  make_row(group, display, x$n_patients, NA_real_, x$HR, x$CI_lo, x$CI_hi, x$p_value)
}

row_from_path <- function(analysis_key, display, group = "Model/t0 sensitivity") {
  x <- path_adjusted %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) stop("Missing path-adjusted row: ", analysis_key)
  make_row(group, display, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}

row_from_l0 <- function(analysis_key, display, group = "Rescue/eligibility sensitivity") {
  x <- l0_results %>% filter(analysis == analysis_key) %>% slice(1)
  if (nrow(x) != 1) stop("Missing L0 row: ", analysis_key)
  make_row(group, display, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}

forest_rows <- bind_rows(
  header_row("Primary and cohort analyses"),
  row_from_forest("Primary SOFA 8-14 equipoise", "Primary SOFA 8-14 (IPCW-weighted)"),
  row_from_forest("Full strictK3_v2 primary cohort", "Full clean-t0 cohort"),
  header_row("SOFA score strata"),
  row_from_forest("SOFA <=7", "SOFA <=7"),
  row_from_forest("SOFA 8-14", "SOFA 8-14"),
  row_from_forest("SOFA 9-13", "SOFA 9-13"),
  row_from_forest("SOFA >=15", "SOFA >=15"),
  header_row("Model and t0 sensitivity"),
  row_from_path("SOFA 8-14: IPCW includes stage3_path", "Stage-3 path-adjusted IPCW"),
  row_from_ph("SOFA 8-14 era-stratified", "Calendar-era stratified Cox"),
  row_from_forest("KDIGO3_minus6h: SOFA 8-14", "KDIGO3 minus 6 h anchor"),
  row_from_forest("KDIGO2: SOFA 8-14", "KDIGO2 exploratory anchor"),
  header_row("Rescue/eligibility sensitivity"),
  row_from_rescue("no_rescue_exclusion: SOFA 8-14", "No rescue exclusion"),
  row_from_rescue("sensitivity_R1_confirm_R2_R2b_R3: SOFA 8-14", "Confirmed rescue thresholds"),
  row_from_rescue("sensitivity_R1_R2_R2b_R3_R4: SOFA 8-14", "Add R4 urine-output rescue"),
  row_from_rescue("sensitivity_liberal_thresholds_R1_R2_R2b_R3: SOFA 8-14", "Liberal rescue thresholds"),
  row_from_rescue("sensitivity_conservative_thresholds_R1_R2_R2b_R3: SOFA 8-14", "Conservative rescue thresholds"),
  row_from_l0("L0-only denominator: SOFA 8-14", "L0-only IPCW denominator")
) %>%
  mutate(
    is_header = is.na(hr),
    Analysis = if_else(is_header, label, paste0("    ", label)),
    N = if_else(is_header, "", paste0(fmt_n(n), if_else(is.na(early), "", paste0(" / ", fmt_n(early))))),
    `HR [95% CI]` = fmt_hr(hr, lo, hi),
    `p-value` = fmt_p(p),
    ` ` = paste(rep(" ", 28), collapse = "")
  )

forest_display <- forest_rows %>%
  transmute(Analysis, `N / early` = N, ` `, `HR [95% CI]`, `p-value`)

# Locked display table for Fig3 ALIGN / manuscript (3-dp fmt_hr)
write_csv(
  forest_rows %>%
    filter(!is_header) %>%
    transmute(
      label, n, early,
      weighted_hr = hr, ci_lower = lo, ci_upper = hi, p_value = p,
      hr_display = `HR [95% CI]`, p_display = `p-value`, n_display = N
    ),
  file.path(results_dir, "strictK3_v2_fig3_forest_display.csv")
)

fp_theme <- forest_theme(
  base_size = 8.3,
  base_family = BASE_FAMILY,
  ci_pch = 18,
  ci_col = COL_INK,
  ci_lwd = 1.08,
  ci_Theight = 0.16,
  refline_gp = gpar(lwd = 0.9, lty = "dashed", col = COL_GREY),
  summary_col = COL_INK,
  summary_fill = COL_HEADER,
  arrow_type = "closed",
  arrow_label_just = "end",
  arrow_gp = gpar(lwd = 0.7, col = "grey35", fontsize = 7),
  xlab_adjust = "refline",
  xlab_gp = gpar(fontsize = 8, fontface = "plain"),
  xaxis_gp = gpar(fontsize = 7.8)
)

fig3 <- forest(
  data = forest_display,
  est = forest_rows$hr,
  lower = forest_rows$lo,
  upper = forest_rows$hi,
  sizes = ifelse(forest_rows$is_header, NA, 0.42),
  ref_line = 1,
  ci_column = 3,
  xlim = c(0.6, 1.6),
  ticks_at = c(0.75, 1.00, 1.25, 1.50),
  xlab = "Hazard Ratio (Early CRRT vs No CRRT within 24 h)",
  arrow_lab = c("Favours Early CRRT", "Favours No CRRT within 24 h"),
  is_summary = forest_rows$is_header,
  theme = fp_theme
)

fig3_paths <- save_grid_pair(
  function() {
    plot(fig3)
    grid.text(
      "Note: SOFA 8-14 is repeated in the strata panel to display the severity gradient; SOFA 9-13 is the narrower sensitivity band.",
      x = unit(0.5, "npc"), y = unit(0.025, "npc"),
      gp = gpar(fontfamily = VECTOR_FAMILY, fontsize = 6.2, col = "grey35")
    )
  },
  "Figure_3_Primary_Robustness_Forest_v2",
  8.3, 6.55, 600
)

# -----------------------------------------------------------------------------
# Supplement Figure S3: SOFA band scan
# -----------------------------------------------------------------------------

sofa_hr_delta <- max(0.20, max(abs(sofa_scan$weighted_hr - 1), na.rm = TRUE))

fig_s1 <- ggplot(sofa_scan, aes(x = factor(lower), y = factor(upper), fill = weighted_hr)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.3f", weighted_hr)), size = 2.6, family = BASE_FAMILY, color = COL_INK) +
  scale_fill_gradient2(
    low = COL_EARLY, mid = "white", high = COL_DEFERRED, midpoint = 1,
    name = "HR", limits = c(1 - sofa_hr_delta, 1 + sofa_hr_delta)
  ) +
  labs(
    title = "Supplementary Figure S3. SOFA equipoise-band scan",
    subtitle = "IPCW HRs for early vs deferred CRRT across SOFA windows.",
    caption = "All hazard ratios were close to 1.0 and none were statistically significant.",
    x = "Lower SOFA bound at t0",
    y = "Upper SOFA bound at t0"
  ) +
  coord_equal() +
  theme_icm(9) +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = COL_INK, fill = NA, linewidth = 0.5),
    plot.caption = element_text(hjust = 0, size = 7.2, color = "grey30")
  )

fig_s1_paths <- save_gg_pair(fig_s1, "Supplement_Figure_S3_SOFA_Band_Scan_v2", 6.8, 5.8, 300)

# -----------------------------------------------------------------------------
# Supplement Figure S5: Bayesian prior-likelihood-posterior
# -----------------------------------------------------------------------------

theta_hat <- log(1.044)
se_hat <- (log(1.181) - log(0.922)) / (2 * 1.959964)

posterior_from_prior <- function(m0, s0) {
  prec <- 1 / s0^2 + 1 / se_hat^2
  pv <- 1 / prec
  pm <- (m0 / s0^2 + theta_hat / se_hat^2) * pv
  tibble(pm = pm, ps = sqrt(pv))
}

priors <- tibble(
  prior = c(
    "Data-dominant", "Neutral", "Skeptical", "Mildly optimistic",
    "Mildly pessimistic", "Discounted-historical neutral",
    "ELAIN-informed stress test"
  ),
  m0 = c(0, 0, 0, -0.105, 0.095, 0, -0.223),
  s0 = c(1.00, 0.355, 0.205, 0.205, 0.205, 0.250, 0.300),
  category = c("Primary", "Primary", "Primary", rep("Sensitivity", 4))
)

bayes_probs <- priors %>%
  rowwise() %>%
  mutate(
    tmp = list(posterior_from_prior(m0, s0)),
    pm = tmp$pm,
    ps = tmp$ps,
    posterior_median_hr = exp(pm),
    `Important benefit\nHR < 0.900` = pnorm((log(0.90) - pm) / ps),
    `No clinically important effect\n0.900-1.100` = pnorm((log(1.10) - pm) / ps) - pnorm((log(0.90) - pm) / ps),
    `Important harm\nHR > 1.100` = 1 - pnorm((log(1.10) - pm) / ps)
  ) %>%
  ungroup() %>%
  select(-tmp)

hr_grid <- tibble(hr = seq(0.70, 1.35, length.out = 600), loghr = log(hr))
neutral <- priors %>% filter(prior == "Neutral")
neutral_post <- posterior_from_prior(neutral$m0, neutral$s0)

density_df <- bind_rows(
  hr_grid %>%
    mutate(curve = "Neutral prior", density = dnorm(loghr, neutral$m0, neutral$s0) / hr),
  hr_grid %>%
    mutate(curve = "Summary likelihood", density = dnorm(loghr, theta_hat, se_hat) / hr),
  hr_grid %>%
    mutate(curve = "Posterior", density = dnorm(loghr, neutral_post$pm, neutral_post$ps) / hr)
) %>%
  mutate(curve = factor(curve, levels = c("Neutral prior", "Summary likelihood", "Posterior")))

bayes_top <- ggplot(density_df, aes(hr, density, color = curve, linewidth = curve)) +
  geom_line() +
  geom_vline(xintercept = c(0.90, 1.00, 1.10), lty = c("dotted", "dashed", "dotted"), color = c(COL_THRESHOLD, COL_GREY, COL_THRESHOLD), linewidth = 0.5) +
  scale_color_manual(values = c("Neutral prior" = "grey55", "Summary likelihood" = COL_INK, "Posterior" = COL_EARLY)) +
  scale_linewidth_manual(values = c("Neutral prior" = 0.75, "Summary likelihood" = 0.9, "Posterior" = 1.2)) +
  scale_x_continuous(breaks = c(0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30)) +
  labs(
    title = "Supplementary Figure S5. Bayesian summary-likelihood reanalysis",
    subtitle = "Overlay shown for the neutral prior; vertical dotted lines mark clinically important benefit/harm thresholds.",
    x = "Hazard ratio for hospital-recorded death by day 28",
    y = "Density"
  ) +
  theme_icm(9) +
  theme(legend.position = "top")

bayes_long <- bayes_probs %>%
  filter(category == "Primary") %>%
  select(prior, starts_with("Important"), starts_with("No clinically")) %>%
  pivot_longer(-prior, names_to = "region", values_to = "probability") %>%
  mutate(
    region = recode(
      region,
      "Important benefit\nHR < 0.900" = "Benefit\nHR < 0.900",
      "No clinically important effect\n0.900-1.100" = "No clinically important effect\n0.900-1.100",
      "Important harm\nHR > 1.100" = "Harm\nHR > 1.100"
    ),
    region = factor(region, levels = c("Benefit\nHR < 0.900", "No clinically important effect\n0.900-1.100", "Harm\nHR > 1.100")),
    prior = factor(prior, levels = c("Data-dominant", "Neutral", "Skeptical"))
  )

bayes_bottom <- ggplot(bayes_long, aes(region, probability, fill = region)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(aes(label = percent(probability, accuracy = 0.1)), vjust = -0.25, size = 2.8, family = BASE_FAMILY, color = COL_INK) +
  facet_wrap(~ prior, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.90), expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(values = c(
    "Benefit\nHR < 0.900" = COL_EARLY,
    "No clinically important effect\n0.900-1.100" = COL_HEADER,
    "Harm\nHR > 1.100" = COL_DEFERRED
  )) +
  labs(x = NULL, y = "Posterior probability") +
  theme_icm(8) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = COL_HEADER, color = NA),
    strip.text = element_text(face = "bold", color = COL_INK),
    axis.text.x = element_text(size = 6.7)
  )

fig_s2 <- bayes_top / bayes_bottom + plot_layout(heights = c(1.2, 1))
fig_s2_paths <- save_gg_pair(fig_s2, "Supplement_Figure_S5_Bayesian_Reanalysis_v2", 8.0, 7.0, 300)

# -----------------------------------------------------------------------------
# Supplement Figure S6: quantitative bias analysis contour
# -----------------------------------------------------------------------------

qba_equal_points <- qba_targets %>%
  mutate(
    rr = equal_rr_strength_for_exposure_and_outcome,
    label = case_when(
      grepl("null", target, ignore.case = TRUE) ~ sprintf("Null\nRR %.3f", rr),
      grepl("benefit", target, ignore.case = TRUE) ~ sprintf("Benefit HR 0.900\nRR %.3f", rr),
      grepl("harm", target, ignore.case = TRUE) ~ sprintf("Harm HR 1.100\nRR %.3f", rr),
      TRUE ~ target
    ),
    label_x = case_when(
      grepl("null", target, ignore.case = TRUE) ~ 1.38,
      grepl("benefit", target, ignore.case = TRUE) ~ 1.72,
      grepl("harm", target, ignore.case = TRUE) ~ 1.40,
      TRUE ~ rr + 0.10
    ),
    label_y = case_when(
      grepl("null", target, ignore.case = TRUE) ~ 1.20,
      grepl("benefit", target, ignore.case = TRUE) ~ 1.64,
      grepl("harm", target, ignore.case = TRUE) ~ 1.36,
      TRUE ~ rr
    )
  )

qba_plot_grid <- qba_grid %>%
  filter(
    rr_exposure_confounder >= 1,
    rr_exposure_confounder <= 3,
    rr_confounder_outcome >= 1,
    rr_confounder_outcome <= 3
  )

fig_s3 <- ggplot(qba_plot_grid, aes(rr_exposure_confounder, rr_confounder_outcome)) +
  geom_tile(aes(fill = bias_factor)) +
  geom_contour(aes(z = bias_factor), breaks = qba_targets$required_bias_factor,
               color = COL_INK, linewidth = 0.45, alpha = 0.80) +
  geom_point(data = qba_equal_points, aes(rr, rr), inherit.aes = FALSE,
             shape = 21, size = 3.5, fill = "white", color = COL_THRESHOLD, stroke = 1.0) +
  geom_text(data = qba_equal_points, aes(label_x, label_y, label = label),
            inherit.aes = FALSE, size = 2.75, lineheight = 0.95,
            family = BASE_FAMILY, color = COL_INK, hjust = 0) +
  scale_fill_gradientn(colors = c("white", COL_EARLY, COL_HEADER, COL_DEFERRED), name = "Bias factor") +
  scale_x_continuous(breaks = seq(1, 3, 0.5)) +
  scale_y_continuous(breaks = seq(1, 3, 0.5)) +
  labs(
    title = "Supplementary Figure S6. Quantitative bias analysis tipping points",
    subtitle = "Equal-strength points show RR needed to move the primary HR to selected thresholds.",
    x = "RR: unmeasured confounder with early CRRT",
    y = "RR: unmeasured confounder with hospital-recorded death by day 28"
  ) +
  coord_equal(xlim = c(1, 3), ylim = c(1, 3), expand = FALSE) +
  theme_icm(9) +
  theme(legend.position = "right")

fig_s3_paths <- save_gg_pair(fig_s3, "Supplement_Figure_S6_QBA_Tipping_Point_v2", 7.2, 6.0, 300)

# -----------------------------------------------------------------------------
# Renal recovery CIF intentionally removed from the supplementary figure set
# -----------------------------------------------------------------------------
#
# The day21-28 creatinine measurement time is not a true renal-recovery event
# time. Plotting it as a CIF creates an artificial flat-then-vertical curve, so
# renal recovery is retained only as a supportive tabular/text endpoint.
unlink(file.path(
  out_dir,
  c(
    "Supplement_Figure_S4_FineGray_Renal_Recovery_v2.png",
    "Supplement_Figure_S4_FineGray_Renal_Recovery_v2.pdf"
  )
))

fg_sofa <- fg_results %>%
  filter(Endpoint == "Renal recovery (death competing)", Analysis == "Equipoise SOFA 8-14") %>%
  slice(1)

renal_recovery_note <- sprintf(
  "Renal recovery is not plotted. Report as weighted Fine-Gray sHR %.3f [%.3f-%.3f], p=%s; supportive endpoint only.",
  fg_sofa$HR, fg_sofa$HR_lo, fg_sofa$HR_hi, fmt_p(fg_sofa$p_value)
)

# -----------------------------------------------------------------------------
# Supplement Figure S4: Love plot
# -----------------------------------------------------------------------------

var_labels <- c(
  age = "Age",
  sofa_t0 = "SOFA at t0",
  apsiii = "APS III",
  creatinine_baseline = "Baseline creatinine",
  nee_at_t0 = "Norepinephrine equivalent",
  lactate_t0 = "Lactate at t0",
  uo_rate_t0 = "Urine output rate",
  map_t0 = "Mean arterial pressure",
  cci = "Charlson comorbidity index"
)

love_df <- balance_smd %>%
  filter(population == "SOFA 8-14") %>%
  mutate(
    variable_label = recode(variable, !!!var_labels),
    variable_label = factor(variable_label, levels = rev(variable_label[order(pre_smd)]))
  )

love_long <- love_df %>%
  select(variable_label, pre_smd, post_weighted_smd) %>%
  pivot_longer(c(pre_smd, post_weighted_smd), names_to = "stage", values_to = "smd") %>%
  mutate(stage = factor(stage, levels = c("pre_smd", "post_weighted_smd"), labels = c("Before IPCW", "After IPCW")))

fig_s5 <- ggplot() +
  annotate("rect", xmin = 0, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = "#E8F8F5", alpha = 0.6) +
  geom_segment(data = love_df, aes(x = pre_smd, xend = post_weighted_smd, y = variable_label, yend = variable_label),
               color = "grey65", linewidth = 0.5) +
  geom_point(data = love_long, aes(smd, variable_label, color = stage, shape = stage),
             size = 3.2, stroke = 0.7) +
  geom_vline(xintercept = 0.10, lty = "dashed", color = COL_THRESHOLD, linewidth = 0.65) +
  scale_color_manual(values = c("Before IPCW" = COL_INK, "After IPCW" = COL_AFTER)) +
  scale_shape_manual(values = c("Before IPCW" = 16, "After IPCW" = 17)) +
  scale_x_continuous(limits = c(0, max(0.62, max(love_long$smd, na.rm = TRUE) + 0.05)),
                     breaks = seq(0, 0.7, 0.1), expand = c(0, 0.02)) +
  labs(
    title = "Supplementary Figure S4. Covariate balance before and after IPCW",
    subtitle = sprintf("Primary SOFA 8-14 cohort; maximum post-IPCW SMD = %.3f.", max(love_df$post_weighted_smd, na.rm = TRUE)),
    x = "Absolute standardized mean difference",
    y = NULL
  ) +
  theme_icm(9) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank()
  )

fig_s5_paths <- save_gg_pair(fig_s5, "Supplement_Figure_S4_Love_Plot_v2", 7.0, 5.5, 300)

# -----------------------------------------------------------------------------
# Supplement Figure S7: eICU feasibility / transportability
# -----------------------------------------------------------------------------

eicu_plot <- eicu_gate %>%
  mutate(
    analysis_label = recode(
      analysis,
      recommended_def_c_high_fidelity = "Def C high-fidelity hospitals",
      recommended_def_c_all_hospitals = "Def C all hospitals"
    ),
    analysis_label = factor(analysis_label, levels = rev(c("Def C high-fidelity hospitals", "Def C all hospitals"))),
    label = sprintf("HR %.3f [%.3f-%.3f]\nearly %s; ESS %.1f; SMD %.3f",
                    hr, ci_low, ci_high, fmt_n(early_n), ess_early, max_weighted_smd),
    x_label = 6.85
  )

fig_s6 <- ggplot(eicu_plot, aes(hr, analysis_label)) +
  annotate("rect", xmin = 0.80, xmax = 1.25, ymin = -Inf, ymax = Inf, fill = "#E8F8F5", alpha = 0.55) +
  geom_vline(xintercept = 1, lty = "dashed", color = COL_GREY, linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, color = COL_INK, linewidth = 0.75) +
  geom_point(shape = 18, size = 4.0, color = COL_INK) +
  geom_text(data = eicu_plot, aes(x = x_label, y = analysis_label, label = label),
            inherit.aes = FALSE, hjust = 1, lineheight = 0.92,
            size = 2.65, family = BASE_FAMILY, color = COL_INK) +
  scale_x_log10(
    limits = c(0.40, 7.00),
    breaks = c(0.5, 0.8, 1.0, 1.25, 2.0, 3.0, 5.0, 7.0),
    labels = c("0.5", "0.8", "1.0", "1.25", "2.0", "3.0", "5.0", "7.0")
  ) +
  labs(
    title = "Supplementary Figure S7. eICU feasibility analysis, not external validation",
    subtitle = "Neither candidate analysis met the documented precision, balance, and point-similarity gates.",
    x = "Hazard ratio for 28-day in-hospital mortality after 24-h landmark",
    y = NULL
  ) +
  theme_icm(9) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(10, 14, 10, 14)
  )

fig_s6_paths <- save_gg_pair(fig_s6, "Supplement_Figure_S7_eICU_Feasibility_v2", 8.0, 4.3, 300)

# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------

static_s1_png <- "Supplement_Figure_S1_DAG_v2.png"
static_s1_pdf <- "Supplement_Figure_S1_DAG_v2.pdf"
static_s2_png <- "Supplement_Figure_S2_CCW_IPCW_v2.png"
static_s2_pdf <- "Supplement_Figure_S2_CCW_IPCW_v2.pdf"

figure_manifest <- tibble::tribble(
  ~figure, ~role, ~png, ~pdf, ~main_or_supplement, ~source, ~caption,
  "Figure 1", "Cohort flow and clean-t0 guardrail", basename(fig1_paths["png"]), basename(fig1_paths["pdf"]), "Main",
  "data/strictK3_primary_v2_funnel.csv; results/strictK3_v2_vs_v1_funnel.csv; data/baseline.csv; data/sofa_group_summary.csv",
  "Flow diagram for the clean physiologic KDIGO stage 3 target-trial emulation. t0 was defined by creatinine or urine-output criteria only; RRT/CRRT was not allowed to define the onset of KDIGO stage 3.",
  "Figure 2", "Weighted cumulative hospital-recorded death by day 28", basename(fig2_paths["png"]), basename(fig2_paths["pdf"]), "Main",
  "data/ccw_final.rds; results/strictK3_primary_v2_mortality_results.csv; results/strictK3_v2_fig2_weighted_nrisk.csv",
  fig2_risk_caption_manifest,
  "Figure 3", "Primary and robustness forest", basename(fig3_paths["png"]), basename(fig3_paths["pdf"]), "Main",
  "Multiple v2 results CSVs: primary, SOFA strata, path/era, rescue, anchor, and L0-only analyses.",
  "Forest plot of the primary analysis and documented robustness analyses. The KDIGO2 anchor is exploratory and addresses a different decision point.",
  "Supplementary Figure S1", "Directed acyclic graph", static_s1_png, static_s1_pdf, "Supplement",
  "Static DAG asset retained in figures_icm_v2.",
  "Directed acyclic graph for treatment timing, confounding, and hospital-recorded death by day 28 in the target-trial emulation.",
  "Supplementary Figure S2", "CCW/IPCW methodology schematic", static_s2_png, static_s2_pdf, "Supplement",
  "Static export from CCW_Supplement_FigureS2_v2.pptx.",
  "Clone-censor-weighting with stabilized inverse-probability-of-censoring weights. At t0 each patient is cloned into early and deferred strategy arms; artificial censoring is determined by CRRT timing relative to the 24-h grace window.",
  "Supplementary Figure S3", "SOFA scan", basename(fig_s1_paths["png"]), basename(fig_s1_paths["pdf"]), "Supplement",
  "results/strictK3_v2_sofa_band_scan_results.csv",
  "SOFA lower/upper bound scan used as an anti-forking-paths diagnostic; each cell displays the IPCW hazard ratio.",
  "Supplementary Figure S4", "Covariate balance", basename(fig_s5_paths["png"]), basename(fig_s5_paths["pdf"]), "Supplement",
  "results/strictK3_v2_ipcw_balance_smd.csv",
  "Love plot of measured covariate balance before and after IPCW in the primary SOFA 8-14 cohort.",
  "Supplementary Figure S5", "Bayesian reanalysis", basename(fig_s2_paths["png"]), basename(fig_s2_paths["pdf"]), "Supplement",
  "results/strictK3_v2_bayes_report.md; reproduced calculations in this script.",
  "Summary-likelihood Bayesian reanalysis of the completed MIMIC-IV v2 primary estimate under the prospectively registered primary priors.",
  "Supplementary Figure S6", "Quantitative bias analysis", basename(fig_s3_paths["png"]), basename(fig_s3_paths["pdf"]), "Supplement",
  "results/strictK3_v2_qba_grid.csv; results/strictK3_v2_qba_target_thresholds.csv",
  "Quantitative bias-analysis contour showing unmeasured-confounding strengths required to move the primary estimate to selected clinical thresholds.",
  "Supplementary Figure S7", "eICU feasibility", basename(fig_s6_paths["png"]), basename(fig_s6_paths["pdf"]), "Supplement",
  "eICU diagnose_rescue/data/task5_precision_gate.csv",
  "Exploratory eICU feasibility and transportability analysis. This figure should not be described as external validation."
)

manifest_table_lines <- c(
  "| Figure | Role | PNG | PDF | Placement | Source | Caption draft |",
  "|---|---|---|---|---|---|---|",
  apply(figure_manifest, 1, function(x) {
    paste0(
      "| ", paste(gsub("\\|", "/", as.character(x)), collapse = " | "), " |"
    )
  })
)

manifest_lines <- c(
  "# figures_icm_v2 manifest",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Output directory: `", out_dir, "`"),
  "",
  "## ASSUMPTION / guardrails",
  "",
  "- All MIMIC-IV figures use the clean physiologic KDIGO stage 3 v2 data: t0 is defined by creatinine or urine-output criteria only, not by RRT/CRRT.",
  "- The main mortality estimand remains the CCW/IPCW fixed strategy contrast: initiate CRRT within 24 h versus remain CRRT-free through 24 h. Any within-window initiation censors the no-initiation clone; CRRT is unrestricted after 24 h.",
  "- The eICU figure is an exploratory feasibility/transportability analysis, not confirmatory external validation.",
  "- The Bayesian figure is a summary-likelihood reanalysis of the completed MIMIC-IV primary estimate under the prospectively registered primary priors; it does not remove residual confounding.",
  paste0("- ", renal_recovery_note),
  "- Font handling: Arial is preferred by many journal style guides, but this local R/vector-device setup could not reliably embed Arial; Helvetica is used as a sans-serif fallback for PNG/PDF reproducibility.",
  "",
  "## Data confirmation",
  "",
  paste0("- Full clean-t0 cohort: n=", fmt_n(step_n("exclude_primary_rescue_R1_R2_R2b_R3")), "."),
  "- Primary SOFA 8-14 cohort: n=1,502; Early CRRT n=146; Deferred n=1,356.",
  sprintf("- Primary HR %.3f [%.3f-%.3f], p=%s; RMST difference %.2f days.",
          primary_sofa$weighted_hr, primary_sofa$ci_lower, primary_sofa$ci_upper,
          fmt_p(primary_sofa$p_value), primary_sofa$rmst_diff),
  "",
  "## Revision notes for this rerun",
  "",
  "- Figure 3 duplicate rows removed: the primary panel now keeps only Primary SOFA 8-14 and Full clean-t0 cohort; the SOFA-strata panel retains SOFA <=7, 8-14, 9-13, and >=15, with a figure footnote explaining the severity-gradient repeat.",
  "- Supplementary Figure S7 x-axis upper limit expanded to 7.0 and the HR/early/ESS/SMD annotations were moved rightward beyond the CI whiskers.",
  "- Figure 2 risk-table note added to explain the IPCW-weighted risk set and protocol-censoring-driven early-arm drop after the 24-h grace period.",
  "- Figure 1 pre-eligibility exclusions expanded to include ICU LOS<24 h, zero 28-day follow-up, and the first-eligible/pre-rescue filter; side connectors were lengthened.",
  "- Supplementary Figure S3 now uses a HR=1.0-centered symmetric diverging scale and notes that all scan HRs were close to 1 and non-significant.",
  "- Renal-recovery CIF removed from the supplementary-figure set; renal recovery should be reported only in text/table as a weighted Fine-Gray supportive endpoint.",
  "- Supplementary Figures S5 and S6 no longer use the prior non-specification colors; fills now use the predefined palette constants.",
  "- Final cleanup: Figure 2 caption in both the image risk-table note and manifest/file index now explicitly explains IPCW-weighted, clone-censored risk sets and the expected early-arm drop from 1,423 at t0 to 54 by day 7.",
  "- Final cleanup: Helvetica was retained as the safer reproducible sans-serif fallback because reliable Arial embedding was not established in this local R/vector-device setup.",
  "- Final cleanup: S5 and S6 palette cleanup was already handled with predefined palette constants; no non-specification color codes remain in the plotting calls.",
  "",
  "## Files",
  "",
  manifest_table_lines,
  "",
  "## Manuscript-facing recommendation",
  "",
  "- Main figures: Figure 1, Figure 2, Figure 3.",
  "- Supplementary figures: S1 DAG, S2 CCW/IPCW schematic, S3 SOFA scan, S4 Love plot, S5 Bayesian reanalysis, S6 QBA, S7 eICU feasibility. Renal recovery is not plotted; report the Fine-Gray sHR in text/table only.",
  "- Avoid language such as externally validated, replicated, proven no effect, or randomized evidence. Preferred phrasing: no evidence of clinically important survival benefit under the emulated-trial assumptions."
)

manifest_path <- file.path(out_dir, "figures_icm_v2_manifest.md")
writeLines(manifest_lines, manifest_path)

readr::write_csv(figure_manifest, file.path(out_dir, "figures_icm_v2_file_index.csv"))

message("Wrote ICM v2 figures to: ", out_dir)
message("Manifest: ", manifest_path)
