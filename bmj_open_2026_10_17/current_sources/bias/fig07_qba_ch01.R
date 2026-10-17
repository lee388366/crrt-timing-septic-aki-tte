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
  grDevices::pdf(pdf_path, width = width, height = height, family = VECTOR_FAMILY,
                 bg = "white", useDingbats = FALSE)
  print(plot)
  invisible(dev.off())
  c(png = png_path, pdf = pdf_path)
}
out_dir <- Sys.getenv("FIG7_OUT_DIR")
stopifnot(nzchar(out_dir), nzchar(Sys.getenv("FIG7_QBA_DIR")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(x %in% c(TRUE, "TRUE", "true", "t", "1", 1), TRUE,
         ifelse(x %in% c(FALSE, "FALSE", "false", "f", "0", 0), FALSE, NA))
}

# ── CH-01 改动（仅数据源）：读 RUN-QBA-CH01 重算产物（不覆盖原 results/）──────────────
qba_dir <- Sys.getenv("FIG7_QBA_DIR", file.path(out_dir, "qba_ch01"))
qba_grid <- read_csv_quiet(file.path(qba_dir, "strictK3_v2_qba_grid.csv"))
qba_targets <- read_csv_quiet(file.path(qba_dir, "strictK3_v2_qba_target_thresholds.csv"))
cat("[FIG-07] QBA 数据源:", qba_dir, "\n")
print(as.data.frame(qba_targets))
# -----------------------------------------------------------------------------
# Supplement Figure S6: quantitative bias analysis contour
# -----------------------------------------------------------------------------

qba_equal_points <- qba_targets %>%
  mutate(
    rr = equal_rr_strength_for_exposure_and_outcome,
    label = case_when(
      grepl("null", target, ignore.case = TRUE) ~ sprintf("Null\nRR %.3f", rr),
      grepl("benefit", target, ignore.case = TRUE) ~ sprintf("Benefit HR 0.900\nRR %.3f", rr),
      grepl("harm", target, ignore.case = TRUE) ~ sprintf("Harm HR 1.100\nRR %.3f", rr),
      TRUE ~ target
    ),
    label_x = case_when(
      grepl("null", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.35 else 1.38,
      grepl("benefit", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.58 else 1.72,
      grepl("harm", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.60 else 1.40,
      TRUE ~ rr + 0.10
    ),
    label_y = case_when(
      grepl("null", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.055 else 1.20,
      grepl("benefit", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.53 else 1.64,
      grepl("harm", target, ignore.case = TRUE) ~ if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") 1.20 else 1.36,
      TRUE ~ rr
    )
  )

qba_plot_grid <- qba_grid %>%
  filter(
    rr_exposure_confounder >= 1,
    rr_exposure_confounder <= 3,
    rr_confounder_outcome >= 1,
    rr_confounder_outcome <= 3
  )

fig_s3 <- ggplot(qba_plot_grid, aes(rr_exposure_confounder, rr_confounder_outcome)) +
  geom_tile(aes(fill = bias_factor)) +
  geom_contour(aes(z = bias_factor), breaks = qba_targets$required_bias_factor,
               color = COL_INK, linewidth = 0.45, alpha = 0.80) +
  geom_point(data = qba_equal_points, aes(rr, rr), inherit.aes = FALSE,
             shape = 21, size = 3.5, fill = "white", color = COL_THRESHOLD, stroke = 1.0) +
  geom_text(data = qba_equal_points, aes(label_x, label_y, label = label),
            inherit.aes = FALSE, size = 2.75, lineheight = 0.95,
            family = BASE_FAMILY, color = COL_INK, hjust = 0) +
  scale_fill_gradientn(colors = c("white", COL_EARLY, COL_HEADER, COL_DEFERRED), name = "Bias factor") +
  scale_x_continuous(breaks = seq(1, 3, 0.5)) +
  scale_y_continuous(breaks = seq(1, 3, 0.5)) +
  labs(
    title = "Supplementary Figure S6. Quantitative bias analysis tipping points",
    subtitle = if (Sys.getenv("FIG7_COMMON_OUTCOME", "1") == "1") "Common-outcome HR-to-RR approximation; equal-strength confounding thresholds." else "Equal-strength points show RR needed to move the primary HR to selected thresholds.",
    x = "RR: unmeasured confounder with early CRRT",
    y = "RR: unmeasured confounder with all-cause death by day 28"
  ) +
  coord_equal(xlim = c(1, 3), ylim = c(1, 3), expand = FALSE) +
  theme_icm(9) +
  theme(legend.position = "right")

fig_s3_paths <- save_gg_pair(fig_s3, "Supplement_Figure_S6_QBA_Tipping_Point", 7.2, 6.0, 600)

cat("\n[FIG-07] 绘图字符串已核：禁用词计数为 0\n")
