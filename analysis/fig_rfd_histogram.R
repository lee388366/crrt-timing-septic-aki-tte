# =============================================================================
# Diagnostic: histogram of rrt_free_days_28d to inspect raw distribution
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# Use hour-precision rfd_28d derived in script 18, sourced from analysis_master.csv
# rfd_28d = 0 if died_28d; = pmax(0, 28 - crrt_total_days) otherwise (hour-precision)
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)
fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"), TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))
master$early_crrt <- fix_bool(master$early_crrt)
master$died_28d   <- fix_bool(master$died_28d)

# Hour-precision RFD (consistent with 18_weighted_finegray_rfd.R)
master$rrt_hours_28d <- ifelse(
  !is.na(master$crrt_total_days),
  master$crrt_total_days * 24,
  ifelse(!is.na(master$rrt_free_days_28d),
         pmax(0, 28 - master$rrt_free_days_28d) * 24,
         NA_real_)
)
master$rfd_28d <- ifelse(
  !is.na(master$died_28d) & master$died_28d, 0,
  ifelse(!is.na(master$rrt_hours_28d),
         pmax(0, 28 - master$rrt_hours_28d / 24),
         NA_real_)
)

plot_df <- data.frame(
  value = master$rfd_28d,
  group = factor(
    ifelse(master$early_crrt, "Early CRRT (n=500)", "Deferred CRRT (n=2,612)"),
    levels = c("Early CRRT (n=500)", "Deferred CRRT (n=2,612)")
  )
)

# Compute zero% and 28-day% per group
ann <- plot_df %>%
  group_by(group) %>%
  summarise(
    n      = n(),
    pct0   = round(100 * mean(value == 0), 1),
    pct28  = round(100 * mean(value == 28), 1),
    .groups = "drop"
  )
cat("Group summary:\n"); print(ann)

p_hist <- ggplot(plot_df, aes(x = value, fill = group)) +

  geom_histogram(
    binwidth  = 1,
    color     = "white",
    linewidth = 0.25,
    alpha     = 0.85
  ) +

  # vertical reference lines at the two spikes
  geom_vline(xintercept = 0,  linetype = "dashed", color = "#922B21", linewidth = 0.8) +
  geom_vline(xintercept = 28, linetype = "dashed", color = "#1A5276", linewidth = 0.8) +

  annotate("text", x = 0.5,  y = Inf, label = "Spike at 0\n(death / ongoing RRT)",
           color = "#922B21", size = 3, hjust = 0, vjust = 1.4) +
  annotate("text", x = 27.5, y = Inf, label = "Spike at 28\n(survived RRT-free)",
           color = "#1A5276", size = 3, hjust = 1, vjust = 1.4) +

  scale_fill_manual(values = c("Early CRRT (n=500)"     = "#8C8C8C",
                                "Deferred CRRT (n=2,612)" = "#E06666")) +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28), limits = c(-0.5, 28.5)) +

  facet_wrap(~ group, ncol = 1, scales = "free_y") +

  labs(
    title    = "Supplementary Figure S3. Distribution of RRT-Free Days Alive at Day 28",
    subtitle = "Bimodal structure: large spikes at 0 (death/RRT) and 28 (full survival)",
    x        = "RRT-Free Days Alive at Day 28",
    y        = "Count",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, color = "grey45"),
    strip.text      = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none",
    plot.margin     = margin(10, 20, 10, 10)
  )

fig_dir <- "analysis/figures"
pdf(file.path(fig_dir, "FigS3_RFD_Distribution_Diagnostic.pdf"), width = 6.5, height = 6, onefile = FALSE)
print(p_hist); dev.off()
png(file.path(fig_dir, "FigS3_RFD_Distribution_Diagnostic.png"), width = 1950, height = 1800, res = 300)
print(p_hist); dev.off()
cat("FigS3 saved.\n")
