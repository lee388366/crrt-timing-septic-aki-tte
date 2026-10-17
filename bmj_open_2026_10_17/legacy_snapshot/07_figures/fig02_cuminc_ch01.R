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


# ── CH-01 改动 1/4：路径（原 data/results 目录只读；新增 CH-01 产物目录与本次输出目录）──
base_dir    <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
ch01_dir    <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out")
out_dir     <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "10_图件/out/fig02")
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
BASE_FAMILY <- "Helvetica"
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

# P values follow the manuscript's AMA convention: capital P, no leading zero, 2 decimals.
# Three decimals are kept where 2 would give .00 or .05, and for the three values the
# manuscript itself quotes to 3 decimals (P = .046, .364, .604), so that the figure and the
# text read identically digit for digit.
LEGACY_P3 <- c(0.046, 0.364, 0.604)
fmt_p1 <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("< .001")
  for (v in LEGACY_P3) if (round(p, 3) == v) return(sprintf(".%03d", round(p * 1000)))
  r <- round(p, 2)
  if (abs(r - 0.05) < 1e-9 || r == 0 || r >= 0.995) return(sprintf(".%03d", round(p * 1000)))
  sprintf(".%02d", round(p * 100))
}
fmt_p <- function(p) vapply(p, fmt_p1, character(1))

fmt_hr <- function(hr, lo, hi) {
  ifelse(is.na(hr), "", sprintf("%.3f [%.3f-%.3f]", hr, lo, hi))
}

save_gg_pair <- function(plot, stem, width, height, dpi) {
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  invisible(dev.off())
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = VECTOR_FAMILY,
                 bg = "white")
  print(plot)
  invisible(dev.off())
  c(png = png_path, pdf = pdf_path)
}
read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(x %in% c(TRUE, "TRUE", "true", "t", "1", 1), TRUE,
         ifelse(x %in% c(FALSE, "FALSE", "false", "f", "0", 0), FALSE, NA))
}

# ── CH-01 改动 2/4：主分析结果改读 CH-01 产物（列名映射到原变量名，下游零改动）──────
primary_results <- read_csv_quiet(file.path(ch01_dir, "ch01_primary_dod.csv")) %>%
  transmute(
    analysis    = if_else(population == "Equipoise SOFA 8-14",
                          "Primary SOFA 8-14 equipoise", "Full strictK3_v2 primary cohort"),
    weighted_hr = hr,
    ci_lower    = ci_lower,
    ci_upper    = ci_upper,
    p_value     = p_value,
    rmst_diff   = rmst_difference
  )

# ── CH-01 改动 3/4：结局口径 = 28 天全因死亡（deathtime 或 dod，不做出院删失）────────
#    只换随访/事件两列；权重、克隆、协变量、模型、种子一律不动
ccw <- readRDS(file.path(ch01_dir, "ccw_ch01_dod.rds")) %>%
  mutate(
    follow_days = corrected_follow_days,
    event       = as.integer(corrected_event)
  )

# ── CH-01 改动：DPI 可调（默认 300，与原 pipeline 一致；600 时输出新文件、不覆盖 300dpi 版）──
FIG2_DPI <- as.integer(Sys.getenv("FIG2_DPI", "300"))
FIG2_SUFFIX <- if (FIG2_DPI == 300L) "" else paste0("_", FIG2_DPI, "dpi")

primary_sofa <- primary_results %>%
  filter(analysis == "Primary SOFA 8-14 equipoise") %>%
  slice(1)
stopifnot(nrow(primary_sofa) == 1)
# -----------------------------------------------------------------------------
# Figure 2: weighted cumulative 28-day all-cause death
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
if (FIG2_DPI == 300L) {
  readr::write_csv(
    fig2_nrisk_df,
    file.path(out_dir, "strictK3_v2_fig2_weighted_nrisk.csv")
  )
}

fig2_risk_caption_manifest <- paste(
  "IPCW-weighted cumulative 28-day all-cause death in the primary SOFA 8-14 overlap cohort comparing CRRT within 24 hours with no initiation during that window.",
  sprintf(
    "The 'weighted number at risk' row reflects IPCW-weighted, clone-censored risk sets; the early-CRRT arm falls from %s at t₀ to %s by day 7 because early-strategy clones that did not initiate CRRT within the 24-h grace period are censored per protocol, leaving only the (weighted) genuinely early initiators at risk.",
    fmt_n(fig2_early_n0), fmt_n(fig2_early_n7)
  ),
  "This is an expected feature of the clone-censor-weighting design, not attrition or a data error."
)

fig2_risk_caption_plot <- paste0(
  "Note: the 'weighted number at risk' row reflects IPCW-weighted, clone-censored risk sets. Early CRRT falls from ",
  fmt_n(fig2_early_n0), " at t₀\n",
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
  pval = sprintf("IPCW Cox P = %s", fmt_p(primary_sofa$p_value)),
  pval.coord = c(1.0, 0.92),
  pval.size = 4.0,
  ggtheme = theme_icm(9)
)
fig2$plot <- fig2$plot +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(0, 28), expand = c(0.01, 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(
    title = "Figure 2. Weighted cumulative 28-day all-cause death",
    subtitle = sprintf(
      "Primary SOFA 8-14 cohort: HR %.3f [%.3f-%.3f]; RMST difference Early minus Deferred %.2f days",
      primary_sofa$weighted_hr, primary_sofa$ci_lower, primary_sofa$ci_upper, primary_sofa$rmst_diff
    ),
    x = "Days since first physiologic KDIGO stage 3 AKI (t₀)",
    y = "Cumulative all-cause death"
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

fig2_png <- file.path(out_dir, paste0("Figure_2_Weighted_Cumulative_Mortality_v2", FIG2_SUFFIX, ".png"))
fig2_pdf <- file.path(out_dir, paste0("Figure_2_Weighted_Cumulative_Mortality_v2", FIG2_SUFFIX, ".pdf"))
ragg::agg_png(fig2_png, width = 8, height = 6.55, units = "in", res = FIG2_DPI, background = "white")
print(fig2)
invisible(dev.off())
grDevices::cairo_pdf(fig2_pdf, width = 8, height = 6.55, family = VECTOR_FAMILY,
               bg = "white")
print(fig2)
invisible(dev.off())
fig2_paths <- c(png = fig2_png, pdf = fig2_pdf)

# ── CH-01 改动 4/4：校验打印（不写回原 results/）──────────────────────────────────
cat("\n[FIG-02] CH-01 口径 Figure 2 写出：", out_dir, "\n")
cat(sprintf("[FIG-02] 主分析 HR %.4f [%.4f-%.4f] p=%.4f RMST diff=%.4f\n",
            primary_sofa$weighted_hr, primary_sofa$ci_lower, primary_sofa$ci_upper,
            primary_sofa$p_value, primary_sofa$rmst_diff))
