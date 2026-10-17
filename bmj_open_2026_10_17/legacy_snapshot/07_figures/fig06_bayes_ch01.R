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
out_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "10_图件/out/fig06")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_csv_quiet <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(x %in% c(TRUE, "TRUE", "true", "t", "1", 1), TRUE,
         ifelse(x %in% c(FALSE, "FALSE", "false", "f", "0", 0), FALSE, NA))
}
# -----------------------------------------------------------------------------
# Supplement Figure S5: Bayesian prior-likelihood-posterior
# -----------------------------------------------------------------------------

theta_hat <- log(1.04438401026141)   # CH-01 主人群点估计
se_hat <- (log(1.18184543930523) - log(0.922910834711954)) / (2 * 1.959964)   # CH-01 95% CI

posterior_from_prior <- function(m0, s0) {
  prec <- 1 / s0^2 + 1 / se_hat^2
  pv <- 1 / prec
  pm <- (m0 / s0^2 + theta_hat / se_hat^2) * pv
  tibble(pm = pm, ps = sqrt(pv))
}

priors <- tibble(
  prior = c(
    "Data-dominant", "Neutral", "Skeptical", "Mildly optimistic",
    "Mildly pessimistic", "Discounted-historical neutral",
    "ELAIN-informed stress test"
  ),
  m0 = c(0, 0, 0, -0.105, 0.095, 0, -0.223),
  s0 = c(1.00, 0.355, 0.205, 0.205, 0.205, 0.250, 0.300),
  category = c("Primary", "Primary", "Primary", rep("Sensitivity", 4))
)

bayes_probs <- priors %>%
  rowwise() %>%
  mutate(
    tmp = list(posterior_from_prior(m0, s0)),
    pm = tmp$pm,
    ps = tmp$ps,
    posterior_median_hr = exp(pm),
    `Important benefit\nHR < 0.900` = pnorm((log(0.90) - pm) / ps),
    `No clinically important effect\n0.900-1.100` = pnorm((log(1.10) - pm) / ps) - pnorm((log(0.90) - pm) / ps),
    `Important harm\nHR > 1.100` = 1 - pnorm((log(1.10) - pm) / ps)
  ) %>%
  ungroup() %>%
  select(-tmp)

hr_grid <- tibble(hr = seq(0.70, 1.35, length.out = 600), loghr = log(hr))
neutral <- priors %>% filter(prior == "Neutral")
neutral_post <- posterior_from_prior(neutral$m0, neutral$s0)

density_df <- bind_rows(
  hr_grid %>%
    mutate(curve = "Neutral prior", density = dnorm(loghr, neutral$m0, neutral$s0) / hr),
  hr_grid %>%
    mutate(curve = "Summary likelihood", density = dnorm(loghr, theta_hat, se_hat) / hr),
  hr_grid %>%
    mutate(curve = "Posterior", density = dnorm(loghr, neutral_post$pm, neutral_post$ps) / hr)
) %>%
  mutate(curve = factor(curve, levels = c("Neutral prior", "Summary likelihood", "Posterior")))

bayes_top <- ggplot(density_df, aes(hr, density, color = curve, linewidth = curve)) +
  geom_line() +
  geom_vline(xintercept = c(0.90, 1.00, 1.10), lty = c("dotted", "dashed", "dotted"), color = c(COL_THRESHOLD, COL_GREY, COL_THRESHOLD), linewidth = 0.5) +
  scale_color_manual(values = c("Neutral prior" = "grey55", "Summary likelihood" = COL_INK, "Posterior" = COL_EARLY)) +
  scale_linewidth_manual(values = c("Neutral prior" = 0.75, "Summary likelihood" = 0.9, "Posterior" = 1.2)) +
  scale_x_continuous(breaks = c(0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30)) +
  labs(
    title = "Supplementary Figure S5. Bayesian summary-likelihood reanalysis",
    subtitle = "Overlay shown for the neutral prior; vertical dotted lines mark clinically important benefit/harm thresholds.",
    x = "Hazard ratio for all-cause death by day 28",
    y = "Density"
  ) +
  theme_icm(9) +
  theme(legend.position = "top")

bayes_long <- bayes_probs %>%
  filter(category == "Primary") %>%
  select(prior, starts_with("Important"), starts_with("No clinically")) %>%
  pivot_longer(-prior, names_to = "region", values_to = "probability") %>%
  mutate(
    region = recode(
      region,
      "Important benefit\nHR < 0.900" = "Benefit\nHR < 0.900",
      "No clinically important effect\n0.900-1.100" = "No clinically important effect\n0.900-1.100",
      "Important harm\nHR > 1.100" = "Harm\nHR > 1.100"
    ),
    region = factor(region, levels = c("Benefit\nHR < 0.900", "No clinically important effect\n0.900-1.100", "Harm\nHR > 1.100")),
    prior = factor(prior, levels = c("Data-dominant", "Neutral", "Skeptical"))
  )

bayes_bottom <- ggplot(bayes_long, aes(region, probability, fill = region)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  geom_text(aes(label = percent(probability, accuracy = 0.1)), vjust = -0.25, size = 2.8, family = BASE_FAMILY, color = COL_INK) +
  facet_wrap(~ prior, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.90), expand = expansion(mult = c(0, 0.06))) +
  scale_fill_manual(values = c(
    "Benefit\nHR < 0.900" = COL_EARLY,
    "No clinically important effect\n0.900-1.100" = COL_HEADER,
    "Harm\nHR > 1.100" = COL_DEFERRED
  )) +
  labs(x = NULL, y = "Posterior probability") +
  theme_icm(8) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = COL_HEADER, color = NA),
    strip.text = element_text(face = "bold", color = COL_INK),
    axis.text.x = element_text(size = 6.7)
  )

fig_s2 <- bayes_top / bayes_bottom + plot_layout(heights = c(1.2, 1))
fig_s2_paths <- save_gg_pair(fig_s2, "Supplement_Figure_S5_Bayesian_Reanalysis", 8.0, 7.0, 600)

cat("\n[FIG-06] 贝叶斯先验/后验概率表（CH-01 输入）:\n")
print(as.data.frame(bayes_probs %>% select(prior, posterior_median_hr, starts_with("Important"), starts_with("No clinically"))))
cat("[FIG-06] theta_hat =", theta_hat, " se_hat =", se_hat, "\n")
