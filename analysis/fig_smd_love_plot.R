# =============================================================================
# IPCW Weight Diagnostics: Love Plot (Pre vs Post Weighting SMD)
# Pre:  analysis_master.csv, by early_crrt (unweighted)
# Post: ccw_final.rds, by strategy, weighted by ipcw_w_trim
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Load data ─────────────────────────────────────────────────────────────────
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)
fix_bool <- function(x)
  ifelse(x %in% c("t","TRUE","true","1"), TRUE,
  ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))
master$early_crrt <- fix_bool(master$early_crrt)

ccw <- readRDS("analysis/data/ccw_final.rds")

# ── SMD helpers ──────────────────────────────────────────────────────────────
smd_cont_w <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  keep <- !is.na(x) & !is.na(g) & !is.na(w)
  x <- x[keep]; g <- g[keep]; w <- w[keep]
  m1 <- sum(w[g==1] * x[g==1]) / sum(w[g==1])
  m0 <- sum(w[g==0] * x[g==0]) / sum(w[g==0])
  s  <- sqrt((var(x[g==1]) + var(x[g==0])) / 2)
  if (is.na(s) || s == 0) return(NA)
  abs((m1 - m0) / s)
}

smd_bin_w <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  x <- as.integer(as.logical(x))
  keep <- !is.na(x) & !is.na(g) & !is.na(w)
  x <- x[keep]; g <- g[keep]; w <- w[keep]
  p1 <- sum(w[g==1] * x[g==1]) / sum(w[g==1])
  p0 <- sum(w[g==0] * x[g==0]) / sum(w[g==0])
  denom <- sqrt((p1*(1-p1) + p0*(1-p0)) / 2)
  if (is.na(denom) || denom == 0) return(NA)
  abs((p1 - p0) / denom)
}

# ── Pre-weighting: prepare master ─────────────────────────────────────────────
master$gender_female <- as.integer(master$gender == "F")
master$vasopressor   <- as.integer(!is.na(master$vaso_baseline) & master$vaso_baseline > 0)
g_pre <- as.integer(master$early_crrt)

# Variable mapping: name in master → display label → type
pre_vars <- list(
  list(col="age",               label="Age (years)",              type="cont"),
  list(col="sofa_baseline",     label="SOFA score at t0",         type="cont"),
  list(col="apsiii_score",      label="APACHE III score",         type="cont"),
  list(col="scr",               label="Creatinine at t0 (mg/dL)", type="cont"),
  list(col="lactate",           label="Lactate at t0 (mmol/L)",   type="cont"),
  list(col="map",               label="MAP at t0 (mmHg)",         type="cont"),
  list(col="ne_eq_baseline",    label="NEE at t0 (ug/kg/min)",    type="cont"),
  list(col="uo_mlkgph_t0_last6h",label="Urine output (mL/kg/h)", type="cont"),
  list(col="cci",               label="Charlson Comorbidity Index",type="cont"),
  list(col="gender_female",     label="Female sex",               type="bin"),
  list(col="vasopressor",       label="Vasopressor use",          type="bin")
)

pre_res <- lapply(pre_vars, function(v) {
  if (!v$col %in% names(master)) return(list(label=v$label, smd=NA))
  x <- master[[v$col]]
  s <- if (v$type == "cont") smd_cont_w(x, g_pre) else smd_bin_w(x, g_pre)
  list(label=v$label, smd=s)
})

# ── Post-weighting: prepare CCW baseline snapshot ────────────────────────────
# One row per (stay_id, strategy) = baseline characteristics
ccw_base <- ccw %>%
  group_by(stay_id, strategy) %>%
  slice(1) %>%
  ungroup()

ccw_base$gender_female <- as.integer(ccw_base$gender == "F")
ccw_base$vasopressor   <- as.integer(ccw_base$vaso_in_window == TRUE)

g_post <- as.integer(ccw_base$strategy == "A_early")
w_post <- ccw_base$ipcw_w_trim

# Variable mapping: name in ccw → label → type  (same order as pre_vars)
post_vars <- list(
  list(col="age",            label="Age (years)",              type="cont"),
  list(col="sofa_t0",        label="SOFA score at t0",         type="cont"),
  list(col="apsiii",         label="APACHE III score",         type="cont"),
  list(col="creatinine_t0",  label="Creatinine at t0 (mg/dL)", type="cont"),
  list(col="lactate_t0",     label="Lactate at t0 (mmol/L)",   type="cont"),
  list(col="map_t0",         label="MAP at t0 (mmHg)",         type="cont"),
  list(col="nee_at_t0",      label="NEE at t0 (ug/kg/min)",    type="cont"),
  list(col="uo_rate_t0",     label="Urine output (mL/kg/h)",   type="cont"),
  list(col="cci",            label="Charlson Comorbidity Index",type="cont"),
  list(col="gender_female",  label="Female sex",               type="bin"),
  list(col="vasopressor",    label="Vasopressor use",          type="bin")
)

post_res <- lapply(post_vars, function(v) {
  if (!v$col %in% names(ccw_base)) return(list(label=v$label, smd=NA))
  x <- ccw_base[[v$col]]
  s <- if (v$type == "cont") smd_cont_w(x, g_post, w_post) else smd_bin_w(x, g_post, w_post)
  list(label=v$label, smd=s)
})

# ── Assemble data frame ───────────────────────────────────────────────────────
labels_pre  <- sapply(pre_res,  `[[`, "label")
labels_post <- sapply(post_res, `[[`, "label")
smd_pre     <- as.numeric(sapply(pre_res,  `[[`, "smd"))
smd_post    <- as.numeric(sapply(post_res, `[[`, "smd"))

love_df <- data.frame(
  label = labels_pre,
  pre   = smd_pre,
  post  = smd_post,
  stringsAsFactors = FALSE
)
love_df <- love_df[!is.na(love_df$pre) & !is.na(love_df$post), ]
love_df <- love_df[order(love_df$pre, decreasing = FALSE), ]
love_df$label <- factor(love_df$label, levels = love_df$label)

cat("\nBalance summary:\n")
cat(sprintf("  Pre-weighting  SMD < 0.1: %d/%d (%.0f%%)\n",
    sum(love_df$pre < 0.1), nrow(love_df), 100*mean(love_df$pre < 0.1)))
cat(sprintf("  Post-weighting SMD < 0.1: %d/%d (%.0f%%)\n",
    sum(love_df$post < 0.1), nrow(love_df), 100*mean(love_df$post < 0.1)))
cat(sprintf("  Max post-weighting SMD: %.3f\n", max(love_df$post, na.rm=TRUE)))
cat("\nDetailed SMDs:\n")
print(data.frame(
  Variable = as.character(love_df$label),
  Pre      = round(love_df$pre, 3),
  Post     = round(love_df$post, 3)
), row.names = FALSE)

# ── Love Plot ─────────────────────────────────────────────────────────────────
love_long <- rbind(
  data.frame(label=love_df$label, smd=love_df$pre,  timing="Before IPCW weighting",
             stringsAsFactors=FALSE),
  data.frame(label=love_df$label, smd=love_df$post, timing="After IPCW weighting",
             stringsAsFactors=FALSE)
)
love_long$timing <- factor(love_long$timing,
                            levels=c("Before IPCW weighting","After IPCW weighting"))
love_long$label  <- factor(love_long$label, levels=levels(love_df$label))

x_max <- max(love_long$smd, na.rm=TRUE) * 1.12 + 0.02

p_love <- ggplot(love_long, aes(x = smd, y = label, color = timing, shape = timing)) +

  # background shading for "well balanced" zone
  annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf,
           fill = "#E8F8F5", alpha = 0.6) +

  # Reference lines
  geom_vline(xintercept = 0,   color = "black",   linewidth = 0.5) +
  geom_vline(xintercept = 0.1, color = "#C0392B", linetype = "dashed",
             linewidth = 0.9, alpha = 0.9) +

  # Connector segments (pre → post per variable)
  geom_segment(
    data    = love_df,
    inherit.aes = FALSE,
    mapping = aes(x = pre, xend = post, y = label, yend = label),
    color = "grey65", linewidth = 0.5, linetype = "solid"
  ) +

  geom_point(size = 3.2, stroke = 0.7) +

  # SMD = 0.10 label anchored to top variable row (avoids bottom clipping)
  annotate("text",
           x     = 0.102,
           y     = nrow(love_df),          # top variable row
           label = "SMD = 0.10\n(balance threshold)",
           color = "#C0392B", size = 2.7,
           hjust = 0, vjust = 1.3, lineheight = 1.15) +

  scale_color_manual(
    values = c("Before IPCW weighting" = "#2C3E50",
               "After IPCW weighting"  = "#27AE60"),
    guide  = guide_legend(reverse = FALSE)
  ) +
  scale_shape_manual(
    values = c("Before IPCW weighting" = 16,
               "After IPCW weighting"  = 17),
    guide  = guide_legend(reverse = FALSE)
  ) +

  scale_x_continuous(
    limits = c(0, x_max),
    breaks = seq(0, ceiling(x_max*10)/10, 0.1),
    expand = c(0, 0.005)
  ) +
  labs(
    title    = "Supplementary Figure S2. Covariate Balance Before and After IPCW Weighting",
    subtitle = "Love plot: absolute standardized mean differences (SMD) by treatment strategy",
    x        = "Absolute Standardized Mean Difference",
    y        = NULL,
    color    = NULL,
    shape    = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle   = element_text(size = 9, color = "grey45", hjust = 0, margin = margin(b=8)),
    axis.text.y     = element_text(size = 9.5, color = "black"),
    axis.text.x     = element_text(size = 9.5),
    axis.title.x    = element_text(size = 10, face = "bold", margin = margin(t=6)),
    legend.position = "top",
    legend.text     = element_text(size = 9),
    legend.key.size = unit(12, "pt"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    plot.margin     = margin(t=8, r=25, b=10, l=10)
  )

# ── Save ──────────────────────────────────────────────────────────────────────
fig_dir <- "analysis/figures"
pdf(file.path(fig_dir, "FigS2_LovePlot_IPCW_Balance.pdf"), width=7, height=5.5, onefile=FALSE)
print(p_love); dev.off()
png(file.path(fig_dir, "FigS2_LovePlot_IPCW_Balance.png"), width=2100, height=1650, res=300)
print(p_love); dev.off()
cat("FigS2 Love Plot saved.\n")
