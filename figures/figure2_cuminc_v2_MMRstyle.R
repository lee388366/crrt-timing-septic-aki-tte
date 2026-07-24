# =============================================================================
# Figure 2 (v2) restyled to the mother-template Figure 2 layout ("改自源件")
#
# Same survminer/ggsurvplot engine + same v2 data object (ccw_final.rds) used by
# strictK3_v2_figures_icm_regen.R, but the LAYOUT is harmonised to the mother
# template's Figure_2_Cumulative_Mortality.png conventions:
#   - centered bold in-figure title (mother hjust=0.5), no subtitle
#   - legend at TOP, titled "Strategy" (mother)
#   - 95% CI shaded bands ON (mother had bands; also makes docx caption correct)
#   - plain "p = ..." annotation (mother style)
#   - concise x-axis "Days from t0"; no in-figure footnote (moved to caption)
# v2 SCIENCE preserved: hospital-recorded death, Deferred/Early strategy names,
#   weighted (clone-censored) risk table, IPCW Cox p = 0.497, red/teal palette.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
  library(ggplot2); library(scales); library(survival); library(survminer); library(ragg)
})

base_dir <- "analysis/strictK3_primary_v2"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
out_dir <- "figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_EARLY <- "#4DBBD5"; COL_DEFERRED <- "#E64B35"; COL_INK <- "#2C3E50"
# PNG uses Arial (matches mother + regenerated Fig1/Fig3); PDF keeps Helvetica
# (Arial is unreliable in the vector device on this R build).
RASTER_FAMILY <- "Arial"; VECTOR_FAMILY <- "Helvetica"

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# theme mirrors strictK3_v2 theme_icm but with CENTERED title (mother convention)
theme_icm_centered <- function(base_size = 9, fam = "Helvetica") {
  theme_classic(base_size = base_size, base_family = fam) +
    theme(
      text = element_text(color = COL_INK),
      axis.text = element_text(color = COL_INK),
      axis.title = element_text(color = COL_INK),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      legend.title = element_text(color = COL_INK),
      legend.text = element_text(color = COL_INK),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

ccw <- readRDS(file.path(data_dir, "ccw_final.rds"))
primary_results <- read_csv(file.path(results_dir, "strictK3_primary_v2_mortality_results.csv"),
                            show_col_types = FALSE, progress = FALSE)
primary_sofa <- primary_results %>% filter(analysis == "Primary SOFA 8-14 equipoise") %>% slice(1)
stopifnot(nrow(primary_sofa) == 1)

ccw_sofa <- ccw %>%
  filter(equipoise == TRUE) %>%
  mutate(
    strategy_plot = factor(strategy, levels = c("B_deferred", "A_early"),
                           labels = c("Deferred CRRT", "Early CRRT")),
    event = as.integer(event)
  ) %>%
  filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim), !is.na(strategy_plot))

fit_mortality <- survfit(Surv(follow_days, event) ~ strategy_plot,
                         data = ccw_sofa, weights = ipcw_w_trim)

build_fig2 <- function(fam) {
  f <- ggsurvplot(
    fit_mortality, data = ccw_sofa, fun = "event",
    risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.22,
    risk.table.title = "Weighted number at risk", risk.table.fontsize = 3.5,
    conf.int = TRUE, conf.int.alpha = 0.15, censor = FALSE,
    palette = c(COL_DEFERRED, COL_EARLY),
    legend.title = "Strategy",
    legend.labs = c("Deferred CRRT", "Early CRRT"),
    legend = "top",
    break.time.by = 7, xlim = c(0, 28), ylim = c(0, 1),
    pval = sprintf("p = %s", fmt_p(primary_sofa$p_value)),
    pval.coord = c(1.0, 0.92), pval.size = 4.2,
    ggtheme = theme_icm_centered(9, fam)
  )
  f$plot <- f$plot +
    scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(0, 28), expand = c(0.01, 0.01)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(
      title = "Figure 2. Weighted cumulative hospital-recorded death by day 28",
      x = "Days from t0",
      y = "Cumulative incidence of hospital-recorded death (%)"
    ) +
    theme(legend.position = "top", text = element_text(family = fam))
  f$table <- f$table +
    theme_icm_centered(8, fam) +
    labs(y = "Strategy") +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
      axis.title.x = element_blank()
    )
  f
}

png_path <- file.path(out_dir, "Figure_2_Cumulative_Mortality_v2_MMRstyle.png")
pdf_path <- file.path(out_dir, "Figure_2_Cumulative_Mortality_v2_MMRstyle.pdf")
agg_png(png_path, width = 8, height = 6.55, units = "in", res = 600, background = "white")
print(build_fig2(RASTER_FAMILY)); invisible(dev.off())
grDevices::pdf(pdf_path, width = 8, height = 6.55, family = VECTOR_FAMILY, bg = "white", useDingbats = FALSE)
print(build_fig2(VECTOR_FAMILY)); invisible(dev.off())
message("Wrote: ", png_path)
message("Wrote: ", pdf_path)
