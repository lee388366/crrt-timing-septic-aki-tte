# Fig. S2/S3-type: mean covariate values over follow-up among uncensored, alive clones,
# unweighted vs IPCW-weighted, by clone arm. S2 = primary SOFA 8-14 cohort; S3 = full cohort.
# Original version (October 2026 submission). Layout functions are reused by figureS2S3S4_2026-11.R;
# the data loading and fixed-weight calculation below belong to the superseded analysis and are not used.
#
#
source("/workspace/tte_bmj/scripts/00_common.R")
ccw <- readRDS(file.path(ROOT, "data/corrected_primary/ccw_ch01_dod.rds")) %>%
  mutate(arm = factor(ifelse(strategy == "A_early", "Early CRRT", "Deferred CRRT"), c("Deferred CRRT", "Early CRRT")),
         fu_h = corrected_follow_h, w = unname(ipcw_w_trim))
stopifnot(nrow(ccw) == 4926, !anyNA(ccw$w), !anyNA(ccw$fu_h))
tv <- readRDS(file.path(ROOT, "data/corrected_inputs/tv_hourly.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE)
L0 <- c(age = "Age (years)", sofa_t0 = "SOFA at t\u2080 (score)", apsiii = "APS III (score)", creatinine_baseline = "Baseline Cr (mg/dL)",
        nee_at_t0 = "NEE at t\u2080 (\u00b5g/kg/min)", lactate_t0 = "Lactate at t\u2080 (mmol/L)",
        uo_rate_t0 = "UO at t\u2080 (mL/kg/h)", map_t0 = "MAP at t\u2080 (mmHg)", cci = "CCI (score)")
LK <- c(sofa_h = "SOFA (score)", map_h = "MAP (mmHg)", nee_h = "NEE (\u00b5g/kg/min)", lactate_h = "Lactate (mmol/L)",
        creatinine_h = "Cr (mg/dL)", kdigo_h = "KDIGO stage", fluid_balance_h = "Fluid balance (mL)")
HOURS_L0 <- 0:167   # 7 days, as in the reference layout
at_risk <- function(d, h) d[d$fu_h > h, ]
wmean <- function(x, w) { ok <- !is.na(x); sum(x[ok] * w[ok]) / sum(w[ok]) }
calc <- function(d) {
  l0 <- lapply(HOURS_L0, function(h) {
    r <- at_risk(d, h)
    bind_rows(lapply(names(L0), function(v) r %>% group_by(arm) %>%
      summarise(Unweighted = mean(.data[[v]], na.rm = TRUE), Weighted = wmean(.data[[v]], w), n = n(), .groups = "drop") %>%
      mutate(var = v, hour = h)))
  }) %>% bind_rows()
  lk <- lapply(0:23, function(h) {
    r <- at_risk(d, h) %>% left_join(filter(tv, hour_index == h), by = "stay_id")
    bind_rows(lapply(names(LK), function(v) r %>% group_by(arm) %>%
      summarise(Unweighted = mean(.data[[v]], na.rm = TRUE), Weighted = wmean(.data[[v]], w), n = sum(!is.na(.data[[v]])), .groups = "drop") %>%
      mutate(var = v, hour = h)))
  }) %>% bind_rows()
  bind_rows(l0, lk) %>% pivot_longer(c(Unweighted, Weighted), names_to = "weighting", values_to = "mean")
}
mini <- function(df, v, lab, xmax) {
  ggplot(filter(df, var == v), aes(hour, mean, colour = arm)) +
    geom_vline(xintercept = 24, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_step(linewidth = 0.6) +
    facet_wrap(~weighting, ncol = 1) +
    scale_colour_manual(values = ARM_COLS) +
    scale_x_continuous(breaks = if (xmax > 30) seq(0, 168, 48) else seq(0, 24, 6), limits = c(0, xmax), expand = c(0.01, 0)) +
    labs(x = "Time (h)", y = lab) +
    theme_bmj(7) + theme(strip.background = element_rect(fill = "grey92", colour = NA), strip.text = element_text(size = 6.5, margin = margin(1, 0, 1, 0)),
      axis.title.y = element_text(size = 6.5), legend.position = "none", panel.spacing = unit(1.5, "mm"))
}
make_fig <- function(d, pop_label, n_pat, stem, figlab) {
  df <- calc(d)
  write.csv(df, file.path(OUT, paste0(stem, "_plot_data.csv")), row.names = FALSE)
  pl <- c(Map(function(v, l) mini(df, v, l, 168), names(L0), L0), Map(function(v, l) mini(df, v, l, 23), names(LK), LK))
  leg <- ggplot(data.frame(arm = factor(names(ARM_COLS), names(ARM_COLS)), x = 1, y = 1), aes(x, y, colour = arm)) +
    geom_line() + scale_colour_manual(values = ARM_COLS, name = "Clone arm") + theme_void(8) +
    theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5))
  p <- wrap_plots(pl, ncol = 4) + plot_layout(guides = "collect") &
    theme(legend.position = "right")
  p <- p & theme(legend.position = "right", legend.title = element_text(size = 8)) & labs(colour = "Clone")
  save_fig(p, stem, 10.5, 11)
}
make_fig(filter(ccw, equipoise), "primary SOFA 8\u201314 cohort", n_distinct(ccw$stay_id[ccw$equipoise]), "FigS2_covariate_means_primary_SOFA8_14", "Fig. S2")
make_fig(ccw, "full cohort", n_distinct(ccw$stay_id), "FigS3_covariate_means_full_cohort", "Fig. S3")
