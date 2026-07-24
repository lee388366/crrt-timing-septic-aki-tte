# =============================================================================
# Figure 6 — Sensitivity Analysis Forest Plot (output: Figure6_ForestPlot_Sensitivity)
# Built with forestploter (journal-standard layout)
# Colour palette: unified with other figures
#   Early CRRT   = #4DBBD5 (teal)
#   Deferred     = #E64B35 (vermillion)
#   Point/CI     = #2C3E50 (dark navy)
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(forestploter)
  library(dplyr)
  library(grid)
})

set.seed(20260329)

base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")

d <- readRDS(file.path(base_dir, "data/ccw_final.rds"))
d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))

d1        <- d[!duplicated(d$stay_id), ]
med_sofa  <- median(d1$sofa_t0,      na.rm = TRUE)   # 9
med_aps   <- median(d1$apsiii,        na.rm = TRUE)   # 69
med_creat <- median(d1$creatinine_t0, na.rm = TRUE)   # 1.90

# ── Cox helper ───────────────────────────────────────────────────────────────
run_cox <- function(dat, label) {
  n_stays <- length(unique(dat$stay_id))
  ev      <- sum(dat$event[!duplicated(dat$stay_id)])
  if (n_stays < 30 || ev < 10) return(NULL)
  fit <- tryCatch(
    coxph(Surv(follow_days, event) ~ strategy,
          data = dat, weights = ipcw_w_trim,
          cluster = stay_id, ties = "efron"),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s  <- summary(fit)
  hr <- s$conf.int[1, 1]
  lo <- s$conf.int[1, 3]
  hi <- s$conf.int[1, 4]
  pv <- s$coefficients[1, "Pr(>|z|)"]
  p_fmt <- if (pv < 0.001) "< 0.001" else sprintf("%.3f", pv)
  data.frame(label    = label,
             n        = n_stays,
             events   = ev,
             est      = hr,
             low      = lo,
             high     = hi,
             p_fmt    = p_fmt,
             stringsAsFactors = FALSE)
}

# ── Run analyses ─────────────────────────────────────────────────────────────
d_eq  <- d[!is.na(d$equipoise) & d$equipoise == TRUE, ]
d_sa1 <- d[!is.na(d$sofa_t0) & d$sofa_t0 >= 9 & d$sofa_t0 <= 13, ]

res <- list(
  run_cox(d_eq,  "Primary: Equipoise SOFA 8\u201314 (IPCW-weighted)"),
  run_cox(d,     "Full cohort (IPCW-weighted)"),
  run_cox(d_sa1, "SA-1: Conservative equipoise SOFA 9\u201313"),
  NULL,  # subgroup header
  run_cox(d[!is.na(d$sofa_t0) & d$sofa_t0 < med_sofa,  ],
          sprintf("SOFA < %g (lower severity)", med_sofa)),
  run_cox(d[!is.na(d$sofa_t0) & d$sofa_t0 >= med_sofa, ],
          sprintf("SOFA \u2265 %g (higher severity)", med_sofa)),
  NULL,
  run_cox(d[!is.na(d$creatinine_t0) & d$creatinine_t0 < med_creat,  ],
          sprintf("Creatinine < %.2f mg/dL", med_creat)),
  run_cox(d[!is.na(d$creatinine_t0) & d$creatinine_t0 >= med_creat, ],
          sprintf("Creatinine \u2265 %.2f mg/dL", med_creat)),
  NULL,
  run_cox(d[!is.na(d$apsiii) & d$apsiii < med_aps,  ],
          sprintf("APACHE III < %g", med_aps)),
  run_cox(d[!is.na(d$apsiii) & d$apsiii >= med_aps, ],
          sprintf("APACHE III \u2265 %g", med_aps))
)

# ── Assemble display table ───────────────────────────────────────────────────
header_labels <- c(
  "SA-3: SOFA stratification",
  "SA-3: Creatinine stratification",
  "SA-3: APACHE III stratification"
)
header_positions <- c(4, 7, 10)  # positions in 'res' list (1-indexed)

rows <- list()
h_idx <- 1
for (i in seq_along(res)) {
  if (i %in% header_positions) {
    rows[[length(rows) + 1]] <- data.frame(
      label   = header_labels[h_idx],
      n       = NA_integer_,
      events  = NA_integer_,
      est     = NA_real_,
      low     = NA_real_,
      high    = NA_real_,
      p_fmt   = NA_character_,
      is_header = TRUE,
      stringsAsFactors = FALSE
    )
    h_idx <- h_idx + 1
  }
  if (!is.null(res[[i]])) {
    rows[[length(rows) + 1]] <- cbind(res[[i]], is_header = FALSE)
  }
}

dt <- do.call(rbind, rows)
dt$is_header <- as.logical(dt$is_header)

# ── Build forestploter display columns ───────────────────────────────────────
# Indent subgroup rows
dt$Analysis <- ifelse(
  dt$is_header,
  dt$label,
  ifelse(grepl("^(SOFA|Creatinine|APACHE)", dt$label),
         paste0("    ", dt$label),
         dt$label)
)

dt$N <- ifelse(is.na(dt$n),      "", as.character(dt$n))
dt$E <- ifelse(is.na(dt$events), "", as.character(dt$events))

dt$HR_CI <- ifelse(
  is.na(dt$est), "",
  sprintf("%.3f [%.3f\u2013%.3f]", dt$est, dt$low, dt$high)
)

dt$P <- ifelse(is.na(dt$p_fmt), "", dt$p_fmt)

# Blank column — wider space to hold the CI plot
dt$` ` <- paste(rep(" ", 26), collapse = "")

# Rename for display headers
names(dt)[names(dt) == "E"]      <- "Events"
names(dt)[names(dt) == "HR_CI"]  <- "HR [95% CI]"
names(dt)[names(dt) == "P"]      <- "p-value"

# Column order for display
display_cols <- c("Analysis", "N", "Events", " ", "HR [95% CI]", "p-value")

# ── Theme / aesthetics ────────────────────────────────────────────────────────
ci_col   <- "#2C3E50"
ref_col  <- "grey55"

fp_theme <- forest_theme(
  base_size    = 9,               # slightly smaller base to prevent crowding
  base_family  = "",
  # CI point and line
  ci_pch       = 18,              # diamond shape
  ci_col       = ci_col,
  ci_lwd       = 1.3,
  ci_Theight   = 0.16,
  # Reference line
  refline_gp   = gpar(lwd = 1.1, lty = "dashed", col = ref_col),
  # Summary rows (section headers)
  summary_col  = "#2C3E50",
  summary_fill = "#EAF0F6",
  # Footnote — smaller and contained
  footnote_gp  = gpar(fontsize = 7.5, fontface = "italic", col = "grey40"),
  # Arrow below x-axis
  arrow_type         = "closed",
  arrow_label_just   = "end",
  arrow_gp           = gpar(lwd = 0.8, col = "grey30", fontsize = 7.5),
  # x-axis
  xlab_adjust  = "refline",
  xlab_gp      = gpar(fontsize = 8.5, fontface = "plain"),
  xaxis_gp     = gpar(fontsize = 8)
)

# ── Draw forest plot ──────────────────────────────────────────────────────────
fp <- forest(
  data       = dt[, display_cols],
  est        = dt$est,
  lower      = dt$low,
  upper      = dt$high,
  sizes      = ifelse(is.na(dt$est), NA,
                      pmax(0.35, pmin(0.65, 300 / dt$n))),
  ref_line   = 1,
  ci_column  = 4,
  xlim       = c(0.4, 2.2),
  ticks_at   = c(0.5, 1.0, 1.5, 2.0),
  xlab       = "Hazard Ratio (Early vs Deferred CRRT)",
  arrow_lab  = c("Favours Early CRRT", "Favours Deferred CRRT"),
  is_summary = dt$is_header,
  footnote   = NULL,              # drawn outside the plot via grid.text
  title      = paste0(
    "Figure 6.  Sensitivity Analyses: IPCW-Weighted HR for 28-Day All-Cause Mortality\n",
    "Early vs Deferred CRRT  |  Diamond = point estimate; whiskers = 95% CI; dashed = null (HR = 1.0)"
  ),
  theme      = fp_theme
)

# ── Save ─────────────────────────────────────────────────────────────────────
fig4_pdf <- file.path(fig_dir, "Figure6_ForestPlot_Sensitivity.pdf")
fig4_png <- file.path(fig_dir, "Figure6_ForestPlot_Sensitivity.png")

footnote_txt <- paste0(
  "Abbreviations: IPCW, inverse probability of censoring weighting; HR, hazard ratio; CI, confidence interval.  ",
  sprintf("Cohort medians: SOFA = %g, creatinine = %.2f mg/dL, APACHE III = %g.  ",
          med_sofa, med_creat, med_aps),
  "All Cox models: cluster-robust sandwich SE, clustered on ICU stay."
)

draw_footnote <- function() {
  grid.text(
    footnote_txt,
    x    = unit(0.01, "npc"),
    y    = unit(0.012, "npc"),
    just = "left",
    gp   = gpar(fontsize = 7, fontface = "italic", col = "grey40")
  )
}

tryCatch({
  pdf(fig4_pdf, width = 14, height = 7.2, onefile = FALSE)
  plot(fp)
  draw_footnote()
  dev.off()
  cat("\nFigure 6 (forestploter) PDF saved.\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  cat("PDF error:", conditionMessage(e), "\n")
})

tryCatch({
  png(fig4_png, width = 4200, height = 2160, res = 300, type = "quartz")
  plot(fp)
  draw_footnote()
  dev.off()
  cat("Figure 6 (forestploter) PNG saved.\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  png(fig4_png, width = 4200, height = 2160, res = 300)
  plot(fp)
  draw_footnote()
  dev.off()
  cat("Figure 6 (forestploter) PNG saved (fallback).\n")
})
