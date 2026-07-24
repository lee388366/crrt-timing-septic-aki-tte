suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
})

set.seed(20260329)

base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")

d <- readRDS(file.path(base_dir, "data/ccw_final.rds"))
d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))

d1 <- d[!duplicated(d$stay_id), ]
med_sofa  <- median(d1$sofa_t0,      na.rm = TRUE)
med_aps   <- median(d1$apsiii,        na.rm = TRUE)
med_creat <- median(d1$creatinine_t0, na.rm = TRUE)
cat(sprintf("Medians — SOFA: %.0f | APACHE III: %.0f | Creatinine: %.2f\n",
            med_sofa, med_aps, med_creat))

run_cox <- function(dat, label) {
  n_stays <- length(unique(dat$stay_id))
  ev      <- sum(dat$event[!duplicated(dat$stay_id)])
  if (n_stays < 30 || ev < 10) {
    cat(sprintf("  [SKIP] %s: n=%d events=%d\n", label, n_stays, ev))
    return(NULL)
  }
  fit <- tryCatch(
    coxph(Surv(follow_days, event) ~ strategy,
          data = dat, weights = ipcw_w_trim,
          cluster = stay_id, ties = "efron"),
    error = function(e) { cat("Cox error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(fit)) return(NULL)
  s  <- summary(fit)
  hr <- s$conf.int[1, 1]
  lo <- s$conf.int[1, 3]
  hi <- s$conf.int[1, 4]
  pv <- s$coefficients[1, "Pr(>|z|)"]
  cat(sprintf("  %-50s HR=%.3f [%.3f, %.3f] p=%.3f  n=%d ev=%d\n",
              label, hr, lo, hi, pv, n_stays, ev))
  data.frame(Label=label, N_stays=n_stays, N_events=ev,
             HR=hr, CI_lo=lo, CI_hi=hi, p_value=pv)
}

# ── Equipoise (SOFA 8-14)
d_eq <- d[!is.na(d$equipoise) & d$equipoise == TRUE, ]
res_main <- run_cox(d_eq, "Primary: Equipoise SOFA 8-14 (IPCW-weighted)")

# ── Full cohort
res_full <- run_cox(d, "Full cohort (IPCW-weighted)")

# ── SA-1: Conservative equipoise SOFA 9-13
d_sa1 <- d[!is.na(d$sofa_t0) & d$sofa_t0 >= 9 & d$sofa_t0 <= 13, ]
res_sa1 <- run_cox(d_sa1, "SA-1: Conservative equipoise SOFA 9-13")

# ── SA-3: SOFA stratification
d_sofa_lo <- d[!is.na(d$sofa_t0) & d$sofa_t0 < med_sofa, ]
d_sofa_hi <- d[!is.na(d$sofa_t0) & d$sofa_t0 >= med_sofa, ]
res_sofa_lo <- run_cox(d_sofa_lo, sprintf("SA-3: SOFA < %.0f (lower severity)", med_sofa))
res_sofa_hi <- run_cox(d_sofa_hi, sprintf("SA-3: SOFA >= %.0f (higher severity)", med_sofa))

# ── SA-3: Creatinine
d_cr_lo <- d[!is.na(d$creatinine_t0) & d$creatinine_t0 < med_creat, ]
d_cr_hi <- d[!is.na(d$creatinine_t0) & d$creatinine_t0 >= med_creat, ]
res_creat_lo <- run_cox(d_cr_lo, sprintf("SA-3: Creatinine < %.2f mg/dL", med_creat))
res_creat_hi <- run_cox(d_cr_hi, sprintf("SA-3: Creatinine >= %.2f mg/dL", med_creat))

# ── SA-3: APACHE III
d_aps_lo <- d[!is.na(d$apsiii) & d$apsiii < med_aps, ]
d_aps_hi <- d[!is.na(d$apsiii) & d$apsiii >= med_aps, ]
res_aps_lo <- run_cox(d_aps_lo, sprintf("SA-3: APACHE III < %.0f", med_aps))
res_aps_hi <- run_cox(d_aps_hi, sprintf("SA-3: APACHE III >= %.0f", med_aps))

# ── Assemble
all_res <- bind_rows(
  res_main, res_full, res_sa1,
  data.frame(Label = "-- SA-3: SOFA stratification --",
             N_stays=NA, N_events=NA, HR=NA, CI_lo=NA, CI_hi=NA, p_value=NA),
  res_sofa_lo, res_sofa_hi,
  data.frame(Label = "-- SA-3: Creatinine stratification --",
             N_stays=NA, N_events=NA, HR=NA, CI_lo=NA, CI_hi=NA, p_value=NA),
  res_creat_lo, res_creat_hi,
  data.frame(Label = "-- SA-3: APACHE III stratification --",
             N_stays=NA, N_events=NA, HR=NA, CI_lo=NA, CI_hi=NA, p_value=NA),
  res_aps_lo, res_aps_hi
)

all_res$display_label <- all_res$Label
all_res$is_header     <- grepl("^--", all_res$Label)
all_res$p_fmt <- ifelse(is.na(all_res$p_value), "",
                  ifelse(all_res$p_value < 0.001, "< 0.001",
                  sprintf("%.3f", all_res$p_value)))
n_rows        <- nrow(all_res)
all_res$y_pos <- n_rows:1

# ── X-axis zones (linear scale):
#   [-5, -0.5]  row labels
#   [0.4, 2.5]  forest CI
#   [3.0, 8.5]  stats table
X_LABEL_ANCHOR <- -4.2
X_FOREST_LEFT  <-  0.40
X_FOREST_RIGHT <-  2.50
X_NULL         <-  1.00
X_N            <-  3.2
X_HR           <-  5.4
X_P            <-  7.5
X_RIGHT_EDGE   <-  8.4

fp <- ggplot() +
  # ── Forest zone shading
  annotate("rect",
           xmin = X_FOREST_LEFT - 0.12, xmax = X_FOREST_RIGHT + 0.12,
           ymin = 0.25, ymax = n_rows + 0.75,
           fill = "grey97", color = NA) +
  # ── Null line
  geom_vline(xintercept = X_NULL, linetype = "dashed",
             color = "grey55", linewidth = 0.45) +
  # ── CI whiskers (clip to forest zone for visual neatness)
  geom_errorbarh(
    data = all_res %>% filter(!is_header & !is.na(CI_lo)),
    aes(xmin = pmax(CI_lo, X_FOREST_LEFT - 0.08),
        xmax = pmin(CI_hi, X_FOREST_RIGHT + 0.08),
        y = y_pos),
    height = 0.22, linewidth = 0.78, color = "#2C3E50"
  ) +
  # ── Point estimates
  geom_point(
    data = all_res %>% filter(!is_header & !is.na(HR)),
    aes(x = HR, y = y_pos, size = 1 / log(N_stays + 2)),
    shape = 18, color = "#2C3E50", show.legend = FALSE
  ) +
  scale_size_continuous(range = c(2.5, 5.5)) +
  # ── Row labels
  geom_text(
    data = all_res,
    aes(x = X_LABEL_ANCHOR, y = y_pos, label = display_label,
        fontface = if_else(is_header, "bold.italic", "plain"),
        color    = if_else(is_header, "grey45", "black")),
    hjust = 0, size = 2.8
  ) +
  scale_color_identity() +
  # ── Column headers (above separator line)
  annotate("text", x = X_N,  y = n_rows + 1.05, label = "N",
           hjust = 0.5, size = 2.9, fontface = "bold", color = "grey20") +
  annotate("text", x = X_HR, y = n_rows + 1.05, label = "HR [95% CI]",
           hjust = 0.5, size = 2.9, fontface = "bold", color = "grey20") +
  annotate("text", x = X_P,  y = n_rows + 1.05, label = "p-value",
           hjust = 0.5, size = 2.9, fontface = "bold", color = "grey20") +
  # ── Stats data
  geom_text(
    data = all_res %>% filter(!is_header),
    aes(x = X_N, y = y_pos,
        label = ifelse(!is.na(N_stays), as.character(N_stays), "")),
    hjust = 0.5, size = 2.65, color = "grey30"
  ) +
  geom_text(
    data = all_res %>% filter(!is_header),
    aes(x = X_HR, y = y_pos,
        label = ifelse(!is.na(HR),
                       sprintf("%.3f [%.3f, %.3f]", HR, CI_lo, CI_hi), "")),
    hjust = 0.5, size = 2.65, color = "#2C3E50"
  ) +
  geom_text(
    data = all_res %>% filter(!is_header),
    aes(x = X_P, y = y_pos, label = p_fmt),
    hjust = 0.5, size = 2.65, color = "#2C3E50"
  ) +
  # ── Horizontal rules
  annotate("segment",
           x = X_LABEL_ANCHOR - 0.1, xend = X_RIGHT_EDGE,
           y = n_rows + 0.6, yend = n_rows + 0.6,
           color = "grey55", linewidth = 0.4) +
  annotate("segment",
           x = X_LABEL_ANCHOR - 0.1, xend = X_RIGHT_EDGE,
           y = 0.35, yend = 0.35,
           color = "grey55", linewidth = 0.4) +
  # ── Axis ticks only over the forest zone
  scale_x_continuous(
    limits = c(X_LABEL_ANCHOR - 0.3, X_RIGHT_EDGE),
    breaks = c(0.5, 1.0, 1.5, 2.0),
    labels = c("0.50", "1.00", "1.50", "2.00"),
    name   = paste0(strrep(" ", 60),
                    "Hazard Ratio (Early vs Deferred CRRT)")
  ) +
  scale_y_continuous(limits = c(0.15, n_rows + 1.6), breaks = NULL) +
  labs(
    title    = "Figure 6. Sensitivity Analyses: IPCW-Weighted HR for 28-Day All-Cause Mortality",
    subtitle = "Early vs Deferred CRRT  |  Diamond = HR point estimate; whiskers = 95% CI; dashed line = null (HR = 1.0)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 8.5, color = "grey40", hjust = 0),
    axis.line.y        = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.title.x       = element_text(size = 9),
    axis.text.x        = element_text(size = 8.5),
    axis.line.x        = element_line(color = "grey65"),
    axis.ticks.x       = element_line(color = "grey65"),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.28),
    plot.margin        = margin(10, 14, 10, 8)
  )

# ── Save
fig4_pdf <- file.path(fig_dir, "Figure6_ForestPlot_Sensitivity.pdf")
fig4_png <- file.path(fig_dir, "Figure6_ForestPlot_Sensitivity.png")

tryCatch({
  pdf(fig4_pdf, width = 12, height = 5.5, onefile = FALSE)
  print(fp); dev.off()
  cat("\nFigure 6 PDF saved.\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  cat("PDF error:", conditionMessage(e), "\n")
})

tryCatch({
  png(fig4_png, width = 3600, height = 1650, res = 300, type = "quartz")
  print(fp); dev.off()
  cat("Figure 6 PNG saved.\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  png(fig4_png, width = 3600, height = 1650, res = 300)
  print(fp); dev.off()
  cat("Figure 6 PNG saved (fallback).\n")
})
