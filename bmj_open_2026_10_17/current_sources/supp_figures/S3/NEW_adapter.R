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


# ── CH-01 改动 1/3：路径（原 data/results 只读）───────────────────────────────────────
base_dir    <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
results_dir <- file.path(base_dir, "results")
work_dir    <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out")
bandscan_csv <- Sys.getenv("FIG5_BANDSCAN_CSV", file.path(work_dir, "run15", "行值_bandscan_接口.csv"))
out_dir <- Sys.getenv("CRRT_FIG_S3_OUT")
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

fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

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

# ── CH-01 改动 2/3：数据源 = run15 全量重跑（25 行；旧口径 25 行来自 results/）──────────
# ── CH-01 改动：DPI 可调（默认 300；600 时输出新文件、不覆盖 300dpi 版）─────────────
FIG5_DPI <- as.integer(Sys.getenv("FIG5_DPI", "300"))
FIG5_SUFFIX <- if (FIG5_DPI == 300L) "" else paste0("_", FIG5_DPI, "dpi")

sofa_scan <- readr::read_csv(bandscan_csv, show_col_types = FALSE, progress = FALSE)
cat("[FIG-05] band scan 行数:", nrow(sofa_scan), "| HR 范围", round(min(sofa_scan$weighted_hr),4),
    "-", round(max(sofa_scan$weighted_hr),4), "| 最小 p", signif(min(sofa_scan$p_value),3), "\n")
# -----------------------------------------------------------------------------
# Supplement Figure S3: SOFA band scan (CH-01 = 28-day all-cause mortality)
# -----------------------------------------------------------------------------

sofa_hr_delta <- max(0.20, max(abs(sofa_scan$weighted_hr - 1), na.rm = TRUE))

fig_s1 <- ggplot(sofa_scan, aes(x = factor(lower), y = factor(upper), fill = weighted_hr)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.3f", weighted_hr)), size = 2.6, family = BASE_FAMILY, color = COL_INK) +
  scale_fill_gradient2(
    low = COL_EARLY, mid = "white", high = COL_DEFERRED, midpoint = 1,
    name = "HR", limits = c(1 - sofa_hr_delta, 1 + sofa_hr_delta)
  ) +
  labs(
    title = "Supplementary Figure S3. SOFA equipoise-band scan",
    subtitle = "IPCW HRs for early vs deferred CRRT across SOFA windows.",
    caption = "All hazard ratios were close to 1.0 and none were statistically significant.",
    x = "Lower SOFA bound at t₀",
    y = "Upper SOFA bound at t₀"
  ) +
  coord_equal() +
  theme_icm(9) +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = COL_INK, fill = NA, linewidth = 0.5),
    plot.caption = element_text(hjust = 0, size = 7.2, color = "grey30")
  )

fig_s1_paths <- save_gg_pair(fig_s1, paste0("Supplement_Figure_S3_SOFA_Band_Scan_v3", FIG5_SUFFIX), 6.8, 5.8, FIG5_DPI)

# ── CH-01 改动 3/3：落盘源数据副本，便于核对 ─────────────────────────────────────────
if (FIG5_DPI == 300L) readr::write_csv(sofa_scan, file.path(out_dir, "Supplement_Figure_S3_SOFA_Band_Scan_v3_source.csv"))
cat("[FIG-05] 写出:", out_dir, "\n")
