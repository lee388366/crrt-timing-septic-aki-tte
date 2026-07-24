# =============================================================================
# Rebuild KM / Cumulative-Incidence figures with reference-code risk table style
#
# Reference style:
#   - risk.table.y.text = TRUE   → show strategy labels on y-axis
#   - tables.theme keeps axis lines, ticks, and "Days from t0" x-label
#   - colour palette unified: Deferred=#E64B35, Early=#4DBBD5
#
# Outputs (PNG only — cairo_pdf has macOS issues):
#   Figure2_KM_Equipoise.png
#   FigS1_KM_FullCohort.png
#   Figure2B_CumInc_Equipoise.png
#   FigS1b_CumInc_Full.png
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")

d <- readRDS(file.path(base_dir, "data/ccw_final.rds"))
d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))
d_eq <- d[!is.na(d$equipoise) & d$equipoise == TRUE, ]

COL_DEFERRED <- "#E64B35"
COL_EARLY    <- "#4DBBD5"
PALETTE      <- c(COL_DEFERRED, COL_EARLY)

# ── Shared themes ─────────────────────────────────────────────────────────────
main_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(size = 11, face = "bold", hjust = 0),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9, face = "bold"),
    panel.grid.major   = element_blank(),   # remove all grid lines
    panel.grid.minor   = element_blank()
  )

# Risk table theme: matches reference style — y-axis line visible, clean white background
risk_theme <- theme(
  plot.title       = element_text(size = 9, face = "bold", hjust = 0),
  axis.line.x      = element_line(color = "black", linewidth = 0.5),
  axis.line.y      = element_line(color = "black", linewidth = 0.5),   # ← y-axis line
  axis.ticks.x     = element_line(color = "black", linewidth = 0.4),
  axis.ticks.y     = element_blank(),
  axis.text.x      = element_text(size = 8.5),
  axis.text.y      = element_text(size = 8.5, hjust = 1),
  axis.title.x     = element_text(size = 9),
  axis.title.y     = element_text(size = 9, face = "bold", angle = 90),
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.grid       = element_blank()
)

save_png <- function(p, path, w = 2400, h = 1900) {
  tryCatch({
    png(path, width = w, height = h, res = 300, type = "quartz")
    print(p); dev.off()
    cat("Saved:", basename(path), "\n")
  }, error = function(e) {
    tryCatch(dev.off(), error = function(e2) {})
    png(path, width = w, height = h, res = 300)
    print(p); dev.off()
    cat("Saved (fallback):", basename(path), "\n")
  })
}

# ═════════════════════════════════════════════════════════════════════════════
# SECTION 1: Weighted KM curves
# ═════════════════════════════════════════════════════════════════════════════

# ── Quick Cox p-values ────────────────────────────────────────────────────────
get_pval <- function(dat) {
  fit <- tryCatch(
    coxph(Surv(follow_days, event) ~ strategy,
          data = dat, weights = ipcw_w_trim, cluster = stay_id, ties = "efron"),
    error = function(e) NULL)
  if (is.null(fit)) return("p = NA")
  pv <- summary(fit)$coefficients[1, "Pr(>|z|)"]
  if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
}

pv_eq   <- get_pval(d_eq)
pv_full <- get_pval(d)
cat(sprintf("Equipoise p = %s  |  Full cohort p = %s\n", pv_eq, pv_full))

# 1a. Figure 2 — KM, Equipoise subgroup
fit_eq <- survfit(Surv(follow_days, event) ~ strategy,
                  data = d_eq, weights = ipcw_w_trim)

km_eq <- ggsurvplot(
  fit_eq, data = d_eq,
  conf.int          = TRUE,
  conf.int.alpha    = 0.14,
  risk.table        = TRUE,
  risk.table.col    = "strata",
  risk.table.y.text = TRUE,        # ← show group labels on y-axis
  risk.table.fontsize = 3.2,
  palette           = PALETTE,
  legend.labs       = c("Deferred CRRT", "Early CRRT"),
  legend.title      = "Strategy",
  xlab              = "Days from t0",
  ylab              = "28-Day Survival Probability",
  title             = "Figure 2.  Weighted KM Survival Curves (Equipoise: SOFA 8\u201314)",
  pval              = pv_eq,
  pval.coord        = c(1, 0.12),
  pval.size         = 3.8,
  xlim              = c(0, 28),
  break.time.by     = 7,
  ggtheme           = main_theme,
  tables.theme      = risk_theme   # ← reference-code style
)

save_png(km_eq, file.path(fig_dir, "Figure2_KM_Equipoise.png"))

# 1b. Suppl S1 — KM, Full cohort
fit_full <- survfit(Surv(follow_days, event) ~ strategy,
                    data = d, weights = ipcw_w_trim)

km_full <- ggsurvplot(
  fit_full, data = d,
  conf.int          = TRUE,
  conf.int.alpha    = 0.14,
  risk.table        = TRUE,
  risk.table.col    = "strata",
  risk.table.y.text = TRUE,
  risk.table.fontsize = 3.2,
  palette           = PALETTE,
  legend.labs       = c("Deferred CRRT", "Early CRRT"),
  legend.title      = "Strategy",
  xlab              = "Days from t0",
  ylab              = "28-Day Survival Probability",
  title             = "Supplementary Figure S1.  Weighted KM Curves (Full Cohort; N\u2009=\u20093,112)",
  pval              = pv_full,
  pval.coord        = c(1, 0.12),
  pval.size         = 3.8,
  xlim              = c(0, 28),
  break.time.by     = 7,
  ggtheme           = main_theme,
  tables.theme      = risk_theme
)

save_png(km_full, file.path(fig_dir, "FigS1_KM_FullCohort.png"))

# ═════════════════════════════════════════════════════════════════════════════
# SECTION 2: Cumulative Incidence of Death (1 – KM)
# ═════════════════════════════════════════════════════════════════════════════

make_ci <- function(dat, title_str, pval_str, show_pval = TRUE) {
  fit <- survfit(Surv(follow_days, event) ~ strategy,
                 data = dat, weights = ipcw_w_trim)

  p <- ggsurvplot(
    fit, data = dat,
    fun               = "event",
    conf.int          = TRUE,
    conf.int.alpha    = 0.14,
    risk.table        = TRUE,
    risk.table.col    = "strata",
    risk.table.y.text = TRUE,      # ← show group labels on y-axis
    risk.table.fontsize = 3.2,
    palette           = PALETTE,
    legend.labs       = c("Deferred CRRT", "Early CRRT"),
    legend.title      = "Strategy",
    xlab              = "Days from t0",
    ylab              = "Cumulative Incidence of Death (%)",
    ylim              = c(0, 1),
    title             = title_str,
    pval              = if (show_pval) pval_str else FALSE,
    pval.coord        = c(1, 0.92),
    pval.size         = 3.8,
    xlim              = c(0, 28),
    break.time.by     = 7,
    surv.scale        = "percent",
    ggtheme           = main_theme,
    tables.theme      = risk_theme  # ← reference-code style
  )

  # Keep percent formatting on y-axis
  p$plot <- p$plot +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1)
    )
  p
}

pv_ci_eq   <- get_pval(d_eq)
pv_ci_full <- get_pval(d)

# 2a. Figure 2B — CumInc, Equipoise
ci_eq <- make_ci(
  d_eq,
  "Figure 2B.  Cumulative Incidence of 28-Day All-Cause Death (Equipoise: SOFA 8\u201314)",
  pv_ci_eq
)
save_png(ci_eq, file.path(fig_dir, "Figure2B_CumInc_Equipoise.png"))

# 2b. Suppl S1b — CumInc, Full cohort
ci_full <- make_ci(
  d,
  "Supplementary Figure S1b.  Cumulative Incidence of Death (Full Cohort; N\u2009=\u20093,112)",
  pv_ci_full
)
save_png(ci_full, file.path(fig_dir, "FigS1b_CumInc_Full.png"))

cat("\nAll four figures rebuilt with reference-code risk table style.\n")
