# =============================================================================
# Figure 2B: Cumulative Incidence of 28-Day All-Cause Death (1 - KM)
#
# Rationale: Complements the KM survival curve (Figure 2) by showing the
# cumulative probability of death accumulating upward over time — an
# intuitive format increasingly adopted by high-impact critical care journals.
#
# Colour palette (unified across all study figures):
#   Early CRRT    → #4DBBD5  (teal/steel blue)
#   Deferred CRRT → #E64B35  (vermillion / brick red)
#
# Output:
#   Figure2B_CumInc_Death.pdf   — 2-panel (Equipoise | Full cohort)
#   Figure2B_CumInc_Equipoise.png
#   FigS1b_CumInc_Full.png
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(survminer)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir  <- "analysis"
data_path <- file.path(base_dir, "data/ccw_final.rds")
fig_dir   <- file.path(base_dir, "figures")

# ── Load data ─────────────────────────────────────────────────────────────────
d <- readRDS(data_path)
d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))

d_eq <- d[!is.na(d$equipoise) & d$equipoise == TRUE, ]

cat(sprintf("Full cohort: %d rows | Equipoise: %d rows\n",
            nrow(d), nrow(d_eq)))

# ── Unified visual theme ──────────────────────────────────────────────────────
COL_EARLY    <- "#4DBBD5"   # teal blue  – Early CRRT
COL_DEFERRED <- "#E64B35"   # vermillion – Deferred CRRT
PALETTE      <- c(COL_DEFERRED, COL_EARLY)   # factor order: Deferred, Early

ci_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(size = 11, face = "bold", hjust = 0),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 9),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.title     = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

# ── Helper: build one cumulative-incidence panel ──────────────────────────────
make_ci_panel <- function(dat, title_str, show_pval = FALSE, pval_str = "") {
  fit <- survfit(Surv(follow_days, event) ~ strategy,
                 data = dat, weights = ipcw_w_trim)

  p <- ggsurvplot(
    fit,
    data              = dat,
    fun               = "event",          # 1 - KM = cumulative incidence
    conf.int          = TRUE,
    conf.int.alpha    = 0.14,
    risk.table        = TRUE,
    risk.table.col    = "strata",
    risk.table.y.text = FALSE,
    palette           = PALETTE,
    legend.labs       = c("Deferred CRRT", "Early CRRT"),
    legend.title      = "Strategy",
    xlab              = "Days from t0",
    ylab              = "Cumulative Incidence of Death (%)",
    ylim              = c(0, 1),
    title             = title_str,
    pval              = if (show_pval) pval_str else FALSE,
    pval.coord        = c(2, 0.92),
    pval.size         = 3.5,
    xlim              = c(0, 28),
    break.time.by     = 7,
    surv.scale        = "percent",
    ggtheme           = ci_theme,
    fontsize          = 3.5,
    tables.theme      = theme_cleantable()
  )

  # Axis formatting: percent labels on y
  p$plot <- p$plot +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1)
    )
  p
}

# ── Cox adjusted p-values (re-run quickly) ────────────────────────────────────
fit_cox_eq <- coxph(Surv(follow_days, event) ~ strategy,
                    data = d_eq, weights = ipcw_w_trim,
                    cluster = stay_id, robust = TRUE)
pv_eq_raw <- 2 * pnorm(-abs(coef(fit_cox_eq) / sqrt(diag(vcov(fit_cox_eq)))))
pv_eq_str <- if (pv_eq_raw < 0.001) "p\u2009<\u20090.001" else sprintf("p\u2009=\u2009%.3f", pv_eq_raw)

fit_cox_full <- coxph(Surv(follow_days, event) ~ strategy,
                      data = d, weights = ipcw_w_trim,
                      cluster = stay_id, robust = TRUE)
pv_full_raw <- 2 * pnorm(-abs(coef(fit_cox_full) / sqrt(diag(vcov(fit_cox_full)))))
pv_full_str <- if (pv_full_raw < 0.001) "p < 0.001" else sprintf("p = %.3f", pv_full_raw)

# ── Build panels ──────────────────────────────────────────────────────────────
p_eq   <- make_ci_panel(d_eq, "Figure 2B. Cumulative Incidence of Death (Equipoise: SOFA 8-14)",
                         show_pval = TRUE, pval_str = pv_eq_str)
p_full <- make_ci_panel(d,    "Supplementary. Cumulative Incidence of Death (Full Cohort)",
                         show_pval = TRUE, pval_str = pv_full_str)

# ── Save Figure 2B (primary: equipoise panel only) ────────────────────────────
f1b_eq_pdf <- file.path(fig_dir, "Figure2B_CumInc_Equipoise.pdf")
f1b_eq_png <- file.path(fig_dir, "Figure2B_CumInc_Equipoise.png")

tryCatch({
  cairo_pdf(f1b_eq_pdf, width = 8, height = 6, onefile = FALSE)
  print(p_eq); dev.off()
  cat("Figure 2B PDF saved:", f1b_eq_pdf, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  cat("Figure 2B PDF error:", conditionMessage(e), "\n")
})

tryCatch({
  png(f1b_eq_png, width = 2400, height = 1800, res = 300, type = "quartz")
  print(p_eq); dev.off()
  cat("Figure 2B PNG saved:", f1b_eq_png, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  png(f1b_eq_png, width = 2400, height = 1800, res = 300)
  print(p_eq); dev.off()
  cat("Figure 2B PNG saved (fallback):", f1b_eq_png, "\n")
})

# ── Save Full cohort cumulative incidence (supplementary) ─────────────────────
f1b_full_pdf <- file.path(fig_dir, "FigS1b_CumInc_Full.pdf")
f1b_full_png <- file.path(fig_dir, "FigS1b_CumInc_Full.png")

tryCatch({
  cairo_pdf(f1b_full_pdf, width = 8, height = 6, onefile = FALSE)
  print(p_full); dev.off()
  cat("Figure S1b PDF saved:", f1b_full_pdf, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  cat("Figure S1b PDF error:", conditionMessage(e), "\n")
})

tryCatch({
  png(f1b_full_png, width = 2400, height = 1800, res = 300, type = "quartz")
  print(p_full); dev.off()
  cat("Figure S1b PNG saved:", f1b_full_png, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  png(f1b_full_png, width = 2400, height = 1800, res = 300)
  print(p_full); dev.off()
  cat("Figure S1b PNG saved (fallback):", f1b_full_png, "\n")
})

cat("\nFigure 2B complete.\n")
