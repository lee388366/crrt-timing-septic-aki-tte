# =============================================================================
# Figure 4 (v2): CIF Plots — 28-Day Renal Recovery (Competing Risks)
#
# Colour palette aligned with Figure 2 / Figure 2B:
#   Early CRRT    → #4DBBD5  (teal blue)
#   Deferred CRRT → #E64B35  (vermillion red)
#
# Improvements over v1:
#   • Unified colour palette
#   • 95% CI bands from cuminc() variance
#   • theme_classic + grey gridlines consistent with other figures
#   • Side-by-side PDF (Panel A + B) + individual PNGs
# =============================================================================

suppressPackageStartupMessages({
  library(cmprsk)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir <- "analysis"
fig_dir  <- file.path(base_dir, "figures")

# ── Load data ─────────────────────────────────────────────────────────────────
master <- read.csv(file.path(base_dir, "data/analysis_master.csv"),
                   stringsAsFactors = FALSE)
fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"),  TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))
master$early_crrt         <- fix_bool(master$early_crrt)
master$died_28d           <- fix_bool(master$died_28d)
master$renal_recovery_28d <- fix_bool(master$renal_recovery_28d)

master <- master %>%
  mutate(
    cr_status = case_when(
      renal_recovery_28d == TRUE ~ 1L,
      died_28d == TRUE           ~ 2L,
      TRUE                       ~ 0L
    ),
    cr_time   = follow_time_28d,
    trt       = as.integer(early_crrt),
    equipoise = !is.na(sofa_baseline) & sofa_baseline >= 8 & sofa_baseline <= 14
  )

# ── Unified colour palette ────────────────────────────────────────────────────
COL_EARLY    <- "#4DBBD5"   # teal blue  – Early CRRT
COL_DEFERRED <- "#E64B35"   # vermillion – Deferred CRRT

# ── Unified visual theme ──────────────────────────────────────────────────────
cif_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title         = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle      = element_text(size = 9, hjust = 0, color = "grey40"),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    legend.title       = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin        = margin(8, 16, 8, 8)
  )

# ── Helper: extract CIF + 95% CI from cuminc() ───────────────────────────────
extract_cif_ci <- function(ci_obj, grp_event_str, group_label) {
  g <- ci_obj[[grp_event_str]]
  se <- sqrt(g$var)
  data.frame(
    time  = g$time,
    cif   = g$est,
    ci_lo = pmax(0, g$est - 1.96 * se),
    ci_hi = pmin(1, g$est + 1.96 * se),
    group = group_label
  )
}

# ── Helper: build one CIF ggplot panel ───────────────────────────────────────
make_cif_panel <- function(dat, title_str, subtitle_str = "") {
  cif_obj <- cuminc(
    ftime   = dat$cr_time,
    fstatus = dat$cr_status,
    group   = dat$trt
  )

  # Gray test p-value
  pval <- cif_obj$Tests["1", "pv"]
  pval_label <- if (pval < 0.001) "Gray test p < 0.001" else sprintf("Gray test p = %.3f", pval)

  d1 <- extract_cif_ci(cif_obj, "1 1", "Early CRRT")
  d0 <- extract_cif_ci(cif_obj, "0 1", "Deferred CRRT")

  cif_df <- bind_rows(d1, d0) %>%
    mutate(group = factor(group, levels = c("Deferred CRRT", "Early CRRT")))

  color_map   <- c("Early CRRT" = COL_EARLY, "Deferred CRRT" = COL_DEFERRED)
  ltype_map   <- c("Early CRRT" = "solid",   "Deferred CRRT" = "dashed")
  fill_map    <- c("Early CRRT" = COL_EARLY, "Deferred CRRT" = COL_DEFERRED)

  # N at risk annotation
  n_early    <- sum(dat$trt == 1)
  n_deferred <- sum(dat$trt == 0)
  n_label    <- sprintf("Early CRRT n=%d  |  Deferred CRRT n=%d", n_early, n_deferred)

  # Max y-axis (ceiling to nearest 5%)
  max_cif <- max(cif_df$ci_hi, na.rm = TRUE)
  y_max   <- min(ceiling(max_cif * 100 / 5) * 5 / 100 + 0.03, 0.35)

  ggplot(cif_df, aes(x = time, color = group, fill = group, linetype = group)) +
    geom_ribbon(aes(ymin = ci_lo * 100, ymax = ci_hi * 100),
                alpha = 0.13, color = NA) +
    geom_step(aes(y = cif * 100), linewidth = 1.0) +
    scale_color_manual(values = color_map) +
    scale_fill_manual(values  = fill_map) +
    scale_linetype_manual(values = ltype_map) +
    scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(0, 28.5)) +
    scale_y_continuous(
      limits = c(0, y_max * 100),
      breaks = seq(0, floor(y_max * 100 / 5) * 5, 5),
      labels = function(x) paste0(x, "%")
    ) +
    annotate("text", x = 28.2, y = y_max * 100 * 0.95,
             label = pval_label, hjust = 1, size = 3.0,
             color = "grey30", fontface = "italic") +
    labs(
      title    = title_str,
      subtitle = n_label,
      x        = "Days from t0",
      y        = "Cumulative Incidence of Renal Recovery (%)"
    ) +
    cif_theme
}

# ── Build panels ──────────────────────────────────────────────────────────────
p_full <- make_cif_panel(master,
  "Figure 4A. CIF: 28-Day Renal Recovery (Full Cohort, n = 3,112)")
p_eq   <- make_cif_panel(master %>% filter(equipoise),
  "Figure 4B. CIF: 28-Day Renal Recovery (Equipoise: SOFA 8-14)")

# ── Save combined PDF ─────────────────────────────────────────────────────────
pdf_path <- file.path(fig_dir, "Figure4_CIF_RenalRecovery.pdf")
tryCatch({
  cairo_pdf(pdf_path, width = 10, height = 4.8, onefile = FALSE)
  gridExtra::grid.arrange(p_full, p_eq, ncol = 2)
  dev.off()
  cat("Figure 4 combined PDF saved:", pdf_path, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(e2) {})
  pdf(pdf_path, width = 10, height = 4.8, onefile = FALSE)
  gridExtra::grid.arrange(p_full, p_eq, ncol = 2)
  dev.off()
  cat("Figure 4 combined PDF (fallback):", pdf_path, "\n")
})

# ── Save individual PNGs ──────────────────────────────────────────────────────
save_png <- function(p, path) {
  tryCatch({
    png(path, width = 1600, height = 1400, res = 300, type = "quartz")
    print(p); dev.off()
    cat("PNG saved:", path, "\n")
  }, error = function(e) {
    tryCatch(dev.off(), error = function(e2) {})
    png(path, width = 1600, height = 1400, res = 300)
    print(p); dev.off()
    cat("PNG saved (fallback):", path, "\n")
  })
}

save_png(p_full, file.path(fig_dir, "Figure4A_CIF_Full.png"))
save_png(p_eq,   file.path(fig_dir, "Figure4B_CIF_Equipoise.png"))

cat("\nFigure 4 (v2) complete — colours and CI bands unified.\n")
