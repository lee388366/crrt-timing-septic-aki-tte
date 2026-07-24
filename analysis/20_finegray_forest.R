# =============================================================================
# Figure 5 — IPCW-Weighted Fine-Gray sHR Forest Plot
# Built with forestploter — style identical to Figure 6 (sensitivity forest)
#
# Colour palette:
#   Renal recovery endpoint  = #2E86C1 (blue)
#   Mortality endpoint       = #C0392B (red)
#   Section headers          = #EAF0F6 fill
# =============================================================================

suppressPackageStartupMessages({
  library(forestploter)
  library(dplyr)
  library(grid)
})

base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")
csv_path <- file.path(base_dir, "data/weighted_finegray_results.csv")

if (!file.exists(csv_path)) {
  stop("Input CSV not found: ", csv_path,
       "\nRun analysis/18_weighted_finegray_rfd.R first.")
}

res <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("Loaded Fine-Gray results:\n")
print(res[, c("Endpoint","Analysis","HR","HR_lo","HR_hi","p_value_display")])

# ── Build display table ───────────────────────────────────────────────────────
recovery_full <- res[res$Endpoint == "Renal recovery (death competing)" &
                       res$Analysis == "Full cohort", ]
recovery_eq   <- res[res$Endpoint == "Renal recovery (death competing)" &
                       grepl("Equipoise", res$Analysis), ]
mortality_full <- res[res$Endpoint == "Mortality (renal recovery competing)" &
                        res$Analysis == "Full cohort", ]
mortality_eq   <- res[res$Endpoint == "Mortality (renal recovery competing)" &
                        grepl("Equipoise", res$Analysis), ]

make_row <- function(label, n, events, est, lo, hi, p_fmt, is_header, endpoint_col) {
  data.frame(
    label        = label,
    n            = n,
    events       = events,
    est          = est,
    low          = lo,
    high         = hi,
    p_fmt        = p_fmt,
    is_header    = is_header,
    endpoint_col = endpoint_col,
    stringsAsFactors = FALSE
  )
}

p_chr <- function(p) ifelse(is.na(p), NA_character_,
                             ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p)))

dt <- rbind(
  make_row("Renal Recovery  (death as competing event)",
           NA, NA, NA, NA, NA, NA, TRUE, "#2E86C1"),
  make_row(paste0("    Full cohort (n = ", recovery_full$N_patients, ")"),
           recovery_full$N_patients, recovery_full$N_event,
           recovery_full$HR, recovery_full$HR_lo, recovery_full$HR_hi,
           p_chr(recovery_full$p_value), FALSE, "#2E86C1"),
  make_row(paste0("    Equipoise subgroup, SOFA 8\u201314 (n = ", recovery_eq$N_patients, ")"),
           recovery_eq$N_patients, recovery_eq$N_event,
           recovery_eq$HR, recovery_eq$HR_lo, recovery_eq$HR_hi,
           p_chr(recovery_eq$p_value), FALSE, "#2E86C1"),
  make_row("Mortality  (renal recovery as competing event)",
           NA, NA, NA, NA, NA, NA, TRUE, "#C0392B"),
  make_row(paste0("    Full cohort (n = ", mortality_full$N_patients, ")"),
           mortality_full$N_patients, mortality_full$N_event,
           mortality_full$HR, mortality_full$HR_lo, mortality_full$HR_hi,
           p_chr(mortality_full$p_value), FALSE, "#C0392B"),
  make_row(paste0("    Equipoise subgroup, SOFA 8\u201314 (n = ", mortality_eq$N_patients, ")"),
           mortality_eq$N_patients, mortality_eq$N_event,
           mortality_eq$HR, mortality_eq$HR_lo, mortality_eq$HR_hi,
           p_chr(mortality_eq$p_value), FALSE, "#C0392B")
)

dt$N      <- ifelse(is.na(dt$n),      "", as.character(dt$n))
dt$Events <- ifelse(is.na(dt$events), "", as.character(dt$events))
dt$`sHR [95% CI]` <- ifelse(
  is.na(dt$est), "",
  sprintf("%.3f [%.3f\u2013%.3f]", dt$est, dt$low, dt$high)
)
dt$`p-value` <- ifelse(is.na(dt$p_fmt), "", dt$p_fmt)
dt$` ` <- paste(rep(" ", 26), collapse = "")

names(dt)[names(dt) == "label"] <- "Analysis"
display_cols <- c("Analysis", "N", "Events", " ", "sHR [95% CI]", "p-value")

# ── Theme (identical to Figure 6 forest) ─────────────────────────────────────
fp_theme <- forest_theme(
  base_size      = 9,
  base_family    = "",
  ci_pch         = 18,
  ci_col         = "#2C3E50",
  ci_lwd         = 1.3,
  ci_Theight     = 0.16,
  refline_gp     = gpar(lwd = 1.1, lty = "dashed", col = "grey55"),
  summary_col    = "#2C3E50",
  summary_fill   = "#EAF0F6",
  footnote_gp    = gpar(fontsize = 7.5, fontface = "italic", col = "grey40"),
  arrow_type         = "closed",
  arrow_label_just   = "end",
  arrow_gp           = gpar(lwd = 0.8, col = "grey30", fontsize = 7.5),
  xlab_adjust    = "refline",
  xlab_gp        = gpar(fontsize = 8.5, fontface = "plain"),
  xaxis_gp       = gpar(fontsize = 8)
)

# ── Draw forest plot ──────────────────────────────────────────────────────────
fp <- forest(
  data       = dt[, display_cols],
  est        = dt$est,
  lower      = dt$low,
  upper      = dt$high,
  sizes      = ifelse(is.na(dt$est), NA, 0.5),
  ref_line   = 1,
  ci_column  = 4,
  xlim       = c(0.4, 2.2),
  ticks_at   = c(0.5, 1.0, 1.5, 2.0),
  xlab       = "Subdistribution Hazard Ratio (Early vs Deferred CRRT)",
  arrow_lab  = c("Favours Early CRRT", "Favours Deferred CRRT"),
  is_summary = dt$is_header,
  footnote   = NULL,
  title      = paste0(
    "Figure 5.  Fine-Gray Competing Risks: IPCW-Weighted sHR\n",
    "Early vs Deferred CRRT  |  Diamond = point estimate; whiskers = 95% CI; dashed = null (sHR = 1.0)"
  ),
  theme      = fp_theme
)

# ── Footnote ─────────────────────────────────────────────────────────────────
footnote_txt <- paste0(
  "Abbreviations: sHR, subdistribution hazard ratio; IPCW, inverse probability of censoring weighting; CI, confidence interval.  ",
  "Models: survival::finegray() expanded dataset + weighted coxph (weights: ipcw_w_trim from CCW pipeline); ",
  "cluster-robust sandwich SE, clustered on ICU stay.  Blue = renal recovery endpoint; Red = mortality endpoint."
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

# ── Save ─────────────────────────────────────────────────────────────────────
fig5_pdf <- file.path(fig_dir, "Figure5_FineGray_Forest.pdf")
fig5_png <- file.path(fig_dir, "Figure5_FineGray_Forest.png")

tryCatch({
  pdf(fig5_pdf, width = 14, height = 5, onefile = FALSE)
  plot(fp); draw_footnote()
  dev.off()
  cat("\nFigure 5 PDF saved:", fig5_pdf, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  cat("PDF error:", conditionMessage(e), "\n")
})

tryCatch({
  png(fig5_png, width = 4200, height = 1500, res = 300, type = "quartz")
  plot(fp); draw_footnote()
  dev.off()
  cat("Figure 5 PNG saved:", fig5_png, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  png(fig5_png, width = 4200, height = 1500, res = 300)
  plot(fp); draw_footnote()
  dev.off()
  cat("Figure 5 PNG saved (fallback):", fig5_png, "\n")
})
