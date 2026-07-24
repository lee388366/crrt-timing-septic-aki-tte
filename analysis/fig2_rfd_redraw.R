# =============================================================================
# Figure 3 (Redraw): RRT-Free Days Alive at Day 28
# Style: paired violin + boxplot, offset ±0.18, grey vs red, sig bracket
# Inspired by raincloud plot design (Python reference code)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Data ─────────────────────────────────────────────────────────────────────
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)
fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"), TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))
master$early_crrt <- fix_bool(master$early_crrt)

# ── Color scheme (matching reference: grey group1, red group2) ────────────────
col_box_grey   <- "#8C8C8C"   # b1: early CRRT box border/median
col_vio_grey   <- "#D3D3D3"   # v1: early CRRT violin fill
col_box_red    <- "#C0392B"   # b2: deferred CRRT box border/median
col_vio_red    <- "#F5B7B1"   # v2: deferred CRRT violin fill

BOLD <- 0.9   # base line width

# ── Compute stats for annotation ─────────────────────────────────────────────
A  <- master$rrt_free_days_28d[master$early_crrt == TRUE]
B  <- master$rrt_free_days_28d[master$early_crrt == FALSE]

p_val   <- wilcox.test(A, B)$p.value
sig_mark <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
p_label  <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)

med_A <- median(A); med_B <- median(B)
iqr_A <- quantile(A, c(0.25, 0.75)); iqr_B <- quantile(B, c(0.25, 0.75))
zero_A <- round(100 * mean(A == 0), 1); zero_B <- round(100 * mean(B == 0), 1)
n_A <- length(A); n_B <- length(B)

cat(sprintf("Early:    median=%.0f [%.0f, %.0f], zeros=%.1f%%, n=%d\n",
            med_A, iqr_A[1], iqr_A[2], zero_A, n_A))
cat(sprintf("Deferred: median=%.0f [%.0f, %.0f], zeros=%.1f%%, n=%d\n",
            med_B, iqr_B[1], iqr_B[2], zero_B, n_B))
cat(sprintf("Wilcoxon: %s (%s)\n", sig_mark, p_label))

# ── Build plot data ───────────────────────────────────────────────────────────
# Two subgroups, x-offset ±0.18 around center x = 0
# violin at exact center; boxplot offset slightly

plot_df <- data.frame(
  value = c(A, B),
  group = factor(
    c(rep("Early CRRT", n_A), rep("Deferred CRRT", n_B)),
    levels = c("Early CRRT", "Deferred CRRT")
  )
)

# x positions: Early=1, Deferred=2 → offset for box inside violin
x_early    <- 1
x_deferred <- 2
box_offset <- 0           # box centered on violin
vio_width  <- 0.50
box_width  <- 0.20

# ── Build figure ─────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(x = group, y = value, fill = group, color = group)) +

  # ── Horizontal reference lines (light grey grid, like Python version) ──────
  geom_hline(yintercept = seq(0, 28, 7),
             color = "lightgrey", linetype = "solid",
             linewidth = 0.4, alpha = 0.6) +
  # y = 0 anchor line (bold)
  geom_hline(yintercept = 0,
             color = "grey40", linetype = "dotted",
             linewidth = 1.0) +

  # ── Violin (lower layer) ──────────────────────────────────────────────────
  geom_violin(
    trim        = TRUE,
    scale       = "width",
    width       = vio_width + 0.1,
    linewidth   = BOLD,
    alpha       = 0.75
  ) +

  # ── Boxplot (upper layer, zorder equivalent: draw on top) ────────────────
  geom_boxplot(
    width        = box_width,
    outlier.shape = NA,           # hide outliers (showfliers=False)
    linewidth    = BOLD,
    median.linewidth = 1.5,       # median line thickness
    color        = "black",
    fill         = "white",
    alpha        = 0.85,
    coef         = 1.5
  ) +

  # ── Median diamond point (like reference code's mean marker) ─────────────
  stat_summary(
    fun      = median,
    geom     = "point",
    shape    = 18,                # diamond
    size     = 3.5,
    color    = "black"
  ) +

  # ── Colors (grey & red scheme) ────────────────────────────────────────────
  scale_fill_manual(
    values = c("Early CRRT"    = col_vio_grey,
               "Deferred CRRT" = col_vio_red)
  ) +
  scale_color_manual(
    values = c("Early CRRT"    = col_box_grey,
               "Deferred CRRT" = col_box_red)
  ) +

  # ── Significance bracket annotation ──────────────────────────────────────
  annotate("segment",
           x = 1, xend = 2, y = 30.5, yend = 30.5,
           linewidth = BOLD + 0.1, color = "black") +
  annotate("segment",
           x = 1, xend = 1, y = 30.5, yend = 29.5,
           linewidth = BOLD + 0.1, color = "black") +
  annotate("segment",
           x = 2, xend = 2, y = 30.5, yend = 29.5,
           linewidth = BOLD + 0.1, color = "black") +
  annotate("text",
           x = 1.5, y = 31.5,
           label = paste0(sig_mark, "  ", p_label),
           size = 3.8, color = "black", fontface = "bold") +

  # ── Zero-fraction labels ──────────────────────────────────────────────────
  annotate("text",
           x = 1, y = -1.8,
           label = sprintf("Zero: %.1f%%", zero_A),
           size = 3.2, color = col_box_grey, fontface = "italic") +
  annotate("text",
           x = 2, y = -1.8,
           label = sprintf("Zero: %.1f%%", zero_B),
           size = 3.2, color = col_box_red, fontface = "italic") +

  # ── Median labels inside plot ─────────────────────────────────────────────
  annotate("text",
           x = 1 + 0.18, y = med_A + 1.5,
           label = sprintf("Med=%g", med_A),
           size = 3.0, color = "grey30", hjust = 0) +
  annotate("text",
           x = 2 + 0.18, y = med_B + 1.5,
           label = sprintf("Med=%g", med_B),
           size = 3.0, color = "#922B21", hjust = 0) +

  # ── Axis & labels ─────────────────────────────────────────────────────────
  scale_y_continuous(
    limits = c(-3, 34),
    breaks = seq(0, 28, 7),
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    labels = c(
      "Early CRRT"    = sprintf("Early CRRT\n(n = %d)", n_A),
      "Deferred CRRT" = sprintf("Deferred CRRT\n(n = %d)", n_B)
    )
  ) +
  labs(
    title    = "Figure 3.  RRT-Free Days Alive at Day 28",
    subtitle = "Full analytic cohort (N = 3,112); death assigned 0 days",
    x        = NULL,
    y        = "RRT-Free Days Alive at Day 28"
  ) +

  # ── Theme (clean, bold spines, horizontal grid only) ─────────────────────
  theme_classic(base_size = 13) +
  theme(
    # Title
    plot.title    = element_text(size = 12, face = "bold", hjust = 0, margin = margin(b=4)),
    plot.subtitle = element_text(size = 9,  color = "grey45", hjust = 0, margin = margin(b=8)),
    # Axes
    axis.title.y  = element_text(size = 11, face = "bold", margin = margin(r=8)),
    axis.text.x   = element_text(size = 10, color = "black", lineheight = 1.2),
    axis.text.y   = element_text(size = 10),
    axis.ticks.x  = element_blank(),    # no x tick marks (length=0)
    axis.ticks.y  = element_line(linewidth = BOLD, color = "black"),
    axis.ticks.length = unit(5, "pt"),
    # Spines (bold, like reference BOLD_LINE_WIDTH)
    axis.line     = element_line(linewidth = BOLD + 0.2, color = "black"),
    # Grid: horizontal only, light grey
    panel.grid.major.y = element_blank(),  # already drawn via geom_hline
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white"),
    # Legend: suppress ggplot auto-legend (will add custom below)
    legend.position = "none",
    # Margins
    plot.margin = margin(t=12, r=30, b=10, l=12)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
fig_dir <- "analysis/figures"

pdf_path <- file.path(fig_dir, "Figure3_RRTfree_days.pdf")
png_path <- file.path(fig_dir, "Figure3_RRTfree_days.png")

tryCatch({
  pdf(pdf_path, width = 5.5, height = 6, onefile = FALSE)
  print(p); dev.off()
  cat("PDF saved:", pdf_path, "\n")
}, error = function(e) { tryCatch(dev.off(), error=function(e2){}); cat("PDF warn:", conditionMessage(e), "\n") })

tryCatch({
  png(png_path, width = 1650, height = 1800, res = 300)
  print(p); dev.off()
  cat("PNG saved:", png_path, "\n")
}, error = function(e) {
  tryCatch(dev.off(), error=function(e2){})
  png(png_path, width = 1650, height = 1800, res = 300, type = "quartz")
  print(p); dev.off()
  cat("PNG saved (quartz):", png_path, "\n")
})

cat("Done.\n")
