# =============================================================================
# Figure 3 (v2) forest — mother-template forestploter style + v2 data ("改自源件")
#
# The v2 forest in strictK3_v2_figures_icm_regen.R already inherits the mother
# restyle theme (MMR/scripts/mmr_restyle_figure3_s7_style.R): same forestploter
# engine, diamond CI markers, #2C3E50 ink, #EAF0F6 header fill, dashed grey
# refline, closed arrows. This standalone rebuild keeps that identical style and
# the v2 robustness row-set/data, but harmonises two things:
#   (1) terminology -> "deferred CRRT" (matches Fig1/Fig2/body/Tables 1-4),
#       replacing the legacy "No CRRT within 24 h" x-axis/arrow labels;
#   (2) font -> Arial for the raster PNG (matches mother + regenerated Fig1);
#       Helvetica retained for the vector PDF (Arial unreliable in vector device).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(forestploter); library(grid); library(ragg)
})

results_dir <- "analysis/strictK3_primary_v2/results"
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_INK <- "#2C3E50"; COL_HEADER <- "#EAF0F6"; COL_GREY <- "grey55"

rq <- function(p) read_csv(p, show_col_types = FALSE, progress = FALSE)
fmt_n  <- function(x) format(round(as.numeric(x)), big.mark = ",", scientific = FALSE)
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_hr <- function(hr, lo, hi) ifelse(is.na(hr), "", sprintf("%.3f [%.3f-%.3f]", hr, lo, hi))

forest_source <- rq(file.path(results_dir, "strictK3_v2_forest_plot_source.csv"))
rescue_results <- rq(file.path(results_dir, "strictK3_v2_rescue_sensitivity_model_results.csv"))
qa_ph <- rq(file.path(results_dir, "strictK3_v2_QA_ph_treatment_results.csv"))
path_adjusted <- rq(file.path(results_dir, "strictK3_v2_path_adjusted_ipcw_results.csv"))
l0_results <- rq(file.path(results_dir, "strictK3_v2_l0_only_sensitivity_results.csv"))
grace_results <- rq("analysis/_provenance/数据溯源_grace_period/grace_period_sensitivity_results.csv")

make_row <- function(group, label, n, early, hr, lo, hi, p)
  tibble(group = group, label = label, n = n, early = early, hr = hr, lo = lo, hi = hi, p = p)
header_row <- function(label)
  tibble(group = label, label = label, n = NA_real_, early = NA_real_, hr = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_)
row_from_forest <- function(k, disp, group = NA_character_) {
  x <- forest_source %>% filter(analysis == k) %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(ifelse(is.na(group), x$group, group), disp, x$n_patients, x$early_n, x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
}
row_from_rescue <- function(k, disp, group = "Rescue/eligibility sensitivity") {
  x <- rescue_results %>% filter(analysis == k) %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(group, disp, x$n_patients, x$early_n, x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
}
row_from_ph <- function(k, disp, group = "Model/t0 sensitivity") {
  x <- qa_ph %>% filter(model == k) %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(group, disp, x$n_patients, NA_real_, x$HR, x$CI_lo, x$CI_hi, x$p_value)
}
row_from_path <- function(k, disp, group = "Model/t0 sensitivity") {
  x <- path_adjusted %>% filter(analysis == k) %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(group, disp, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}
row_from_l0 <- function(k, disp, group = "Rescue/eligibility sensitivity") {
  x <- l0_results %>% filter(analysis == k) %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(group, disp, x$n_patients, x$early_n, x$HR, x$CI_lo, x$CI_hi, x$p)
}
row_from_grace <- function(gh, disp, group = "Grace-period window sensitivity") {
  x <- grace_results %>% filter(grace_h == gh, population == "Primary SOFA 8-14") %>% slice(1); stopifnot(nrow(x) == 1)
  make_row(group, disp, x$n_patients, x$early_n, x$weighted_hr, x$ci_lower, x$ci_upper, x$p_value)
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
  row_from_l0("L0-only denominator: SOFA 8-14", "L0-only IPCW denominator"),
  header_row("Grace-period window sensitivity"),
  row_from_grace(12, "Grace window 12 h"),
  row_from_grace(48, "Grace window 48 h")
) %>%
  mutate(
    is_header = is.na(hr),
    Analysis = if_else(is_header, label, paste0("    ", label)),
    N = if_else(is_header, "", paste0(fmt_n(n), if_else(is.na(early), "", paste0(" / ", fmt_n(early))))),
    `HR [95% CI]` = fmt_hr(hr, lo, hi),
    `p-value` = fmt_p(p),
    ` ` = paste(rep(" ", 28), collapse = "")
  )

forest_display <- forest_rows %>% transmute(Analysis, `N / early` = N, ` `, `HR [95% CI]`, `p-value`)

build_theme <- function(fam) forest_theme(
  base_size = 8.3, base_family = fam, ci_pch = 18, ci_col = COL_INK, ci_lwd = 1.08, ci_Theight = 0.16,
  refline_gp = gpar(lwd = 0.9, lty = "dashed", col = COL_GREY),
  summary_col = COL_INK, summary_fill = COL_HEADER,
  arrow_type = "closed", arrow_label_just = "end", arrow_gp = gpar(lwd = 0.7, col = "black", fontsize = 8.3, fontfamily = fam),
  xlab_adjust = "refline", xlab_gp = gpar(fontsize = 8, fontface = "plain"), xaxis_gp = gpar(fontsize = 7.8)
)

build_fig <- function(fam) forest(
  data = forest_display, est = forest_rows$hr, lower = forest_rows$lo, upper = forest_rows$hi,
  sizes = ifelse(forest_rows$is_header, NA, 0.42), ref_line = 1, ci_column = 3,
  xlim = c(0.6, 1.6), ticks_at = c(0.75, 1.00, 1.25, 1.50),
  xlab = "Hazard Ratio (Early CRRT vs deferred CRRT)",
  arrow_lab = c("Favours early CRRT", "Favours deferred CRRT"),
  is_summary = forest_rows$is_header, theme = build_theme(fam)
)

png_path <- file.path(out_dir, "Figure_3_Primary_Robustness_Forest_v2_MMRstyle.png")
pdf_path <- file.path(out_dir, "Figure_3_Primary_Robustness_Forest_v2_MMRstyle.pdf")
agg_png(png_path, width = 8.3, height = 7.45, units = "in", res = 600, background = "white")
plot(build_fig("Arial")); invisible(dev.off())
grDevices::pdf(pdf_path, width = 8.3, height = 7.45, family = "Helvetica", bg = "white", useDingbats = FALSE)
plot(build_fig("Helvetica")); invisible(dev.off())
message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
