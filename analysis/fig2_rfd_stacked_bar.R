# =============================================================================
# Figure 3 (Final): RRT-Free Days at Day 28 — Stacked Bar Chart
# Three categories: Died/ongoing RRT (=0) | Partial (1-27) | Full survival (=28)
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

# ── Categorise ────────────────────────────────────────────────────────────────
master <- master %>%
  mutate(
    group = factor(
      ifelse(early_crrt, "Early CRRT", "Deferred CRRT"),
      levels = c("Early CRRT", "Deferred CRRT")
    ),
    rfd_cat = factor(
      case_when(
        rrt_free_days_28d == 0  ~ "Died or ongoing RRT (= 0 days)",
        rrt_free_days_28d == 28 ~ "Full 28-day RRT-free survival (= 28 days)",
        TRUE                    ~ "Partial RRT-free survival (1\u201327 days)"
      ),
      levels = c(
        "Died or ongoing RRT (= 0 days)",
        "Partial RRT-free survival (1\u201327 days)",
        "Full 28-day RRT-free survival (= 28 days)"
      )
    )
  )

# ── Compute proportions ───────────────────────────────────────────────────────
bar_df <- master %>%
  group_by(group, rfd_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(
    total = sum(n),
    pct   = 100 * n / total,
    label = sprintf("%.1f%%", pct)
  ) %>%
  ungroup()

cat("=== Stacked bar proportions ===\n")
print(bar_df %>% select(group, rfd_cat, n, total, pct))

# Wilcoxon p-value for annotation
A  <- master$rrt_free_days_28d[master$early_crrt == TRUE]
B  <- master$rrt_free_days_28d[master$early_crrt == FALSE]
p_val    <- wilcox.test(A, B)$p.value
p_label  <- if (p_val < 0.001) "Wilcoxon p < 0.001" else sprintf("Wilcoxon p = %.3f", p_val)

# ── Colour palette (matches study colour scheme) ──────────────────────────────
cat_colors <- c(
  "Died or ongoing RRT (= 0 days)"           = "#C0392B",  # red
  "Partial RRT-free survival (1\u201327 days)"     = "#F0B27A",  # orange
  "Full 28-day RRT-free survival (= 28 days)" = "#2E86C1"   # blue
)

# ── Plot ─────────────────────────────────────────────────────────────────────
p_bar <- ggplot(bar_df, aes(x = group, y = pct, fill = rfd_cat)) +

  geom_col(
    width    = 0.55,
    color    = "white",
    linewidth = 0.5
  ) +

  # Percentage labels inside bars (only if segment is wide enough)
  geom_text(
    aes(label = label),
    position  = position_stack(vjust = 0.5),
    size      = 3.6,
    color     = "white",
    fontface  = "bold",
    data      = bar_df %>% filter(pct > 4)
  ) +

  # n= labels below x axis (via subtitle in axis tick)
  scale_x_discrete(
    labels = c(
      "Early CRRT"    = sprintf("Early CRRT\n(n = %d)", sum(master$early_crrt == TRUE, na.rm=TRUE)),
      "Deferred CRRT" = sprintf("Deferred CRRT\n(n = %d)", sum(master$early_crrt == FALSE, na.rm=TRUE))
    )
  ) +

  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(
    limits = c(0, 114),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +

  # Significance bracket — sits clearly above 100% bar top
  # horizontal bar
  annotate("segment", x = 1, xend = 2, y = 105, yend = 105,
           linewidth = 0.7, color = "black") +
  # left tick down
  annotate("segment", x = 1, xend = 1, y = 105, yend = 102,
           linewidth = 0.7, color = "black") +
  # right tick down
  annotate("segment", x = 2, xend = 2, y = 105, yend = 102,
           linewidth = 0.7, color = "black") +
  # p-value text above bracket
  annotate("text", x = 1.5, y = 109,
           label = p_label,
           size = 3.3, color = "black", fontface = "plain") +

  labs(
    title    = "Figure 3.  RRT-Free Days Alive at Day 28 by Treatment Group",
    subtitle = "Patients categorised by outcome: died or ongoing RRT (0 days), partial (1-27 days), full RRT-free survival (28 days)",
    x        = NULL,
    y        = "Proportion of Patients (%)",
    fill     = "Patient category"
  ) +

  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle   = element_text(size = 8.5, color = "grey40", hjust = 0,
                                   margin = margin(b = 8), lineheight = 1.3),
    axis.text.x     = element_text(size = 10.5, color = "black", lineheight = 1.3),
    axis.text.y     = element_text(size = 10),
    axis.title.y    = element_text(size = 10.5, face = "bold", margin = margin(r = 8)),
    axis.ticks.x    = element_blank(),
    axis.line       = element_line(linewidth = 0.7, color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8.5),
    legend.key.size = unit(11, "pt"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    plot.margin     = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ── Save ─────────────────────────────────────────────────────────────────────
fig_dir <- "analysis/figures"

pdf(file.path(fig_dir, "Figure3_RRTfree_days.pdf"), width = 7.5, height = 5.5, onefile = FALSE)
print(p_bar); dev.off()
png(file.path(fig_dir, "Figure3_RRTfree_days.png"), width = 2250, height = 1650, res = 300)
print(p_bar); dev.off()
cat("Figure 3 (stacked bar) saved.\n")
