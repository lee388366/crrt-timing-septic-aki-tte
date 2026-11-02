# Figure S12 (old S5): Bayesian reanalysis under the three primary priors
# Minimal modification of the earlier Bayesian figure adapter
# Posterior locations/SDs and probabilities are READ from current_hr_prior_catalog_results.csv
# (7-row catalog; only the 3 primary priors are shown). The summary likelihood is taken from
# ch01_primary_dod_CORRECTED_CANDIDATE.csv (Equipoise SOFA 8-14), SE = (log CI_hi - log CI_lo)/(2*1.959964)
# exactly as in the adapter. Density transform (dnorm(logHR)/HR) and thresholds 0.90/1.10 unchanged.
# No priors redefined; no model fitted. Probabilities x100 only for display.
SCRIPT_DIR <- local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()})
source(file.path(SCRIPT_DIR, "00_style.R"))

catalog <- read_csv(inp("current_hr_prior_catalog_results.csv"), show_col_types = FALSE)
primary <- read_csv(inp("ch01_primary_dod_CORRECTED_CANDIDATE.csv"), show_col_types = FALSE)
stopifnot(nrow(catalog) == 7, sum(catalog$role == "primary") == 3)

main <- filter(primary, population == "Equipoise SOFA 8-14")
stopifnot(nrow(main) == 1)
theta_hat <- log(main$hr)
se_hat    <- (log(main$ci_upper) - log(main$ci_lower)) / (2 * 1.959964)

# NOTE: in the saved catalog, posterior_location is on the HR scale (= exp(posterior mean logHR));
# posterior_sd is on the logHR scale. Verified below against the probabilities stored in the same table.
cat3 <- catalog %>% filter(role == "primary") %>%
  mutate(pm_log = log(posterior_location), ps_log = posterior_sd,
         benefit_chk = pnorm((log(0.90) - pm_log) / ps_log),
         harm_chk    = 1 - pnorm((log(1.10) - pm_log) / ps_log))
stopifnot(max(abs(cat3$benefit_chk - cat3$benefit)) < 1e-6,
          max(abs(cat3$harm_chk - cat3$harm)) < 1e-6)

neutral  <- filter(cat3, prior == "Neutral")
hr_grid  <- tibble(hr = seq(0.70, 1.35, length.out = 600), loghr = log(hr))
density_df <- bind_rows(
  mutate(hr_grid, curve = "Neutral prior",      density = dnorm(loghr, neutral$prior_mean, neutral$prior_sd) / hr),
  mutate(hr_grid, curve = "Summary likelihood", density = dnorm(loghr, theta_hat, se_hat) / hr),
  mutate(hr_grid, curve = "Posterior",          density = dnorm(loghr, neutral$pm_log, neutral$ps_log) / hr)
) %>% mutate(curve = factor(curve, levels = c("Neutral prior", "Summary likelihood", "Posterior")))

pa <- ggplot(density_df, aes(hr, density, colour = curve, linewidth = curve)) +
  geom_vline(xintercept = c(0.90, 1.10), linetype = "dotted", colour = COL_MID, linewidth = 0.45) +
  geom_vline(xintercept = 1.00, linetype = "dashed", colour = COL_MID, linewidth = 0.45) +
  geom_line() +
  scale_colour_manual(values = c("Neutral prior" = "grey65", "Summary likelihood" = COL_INK,
                                 "Posterior" = COL_BLUE), name = NULL) +
  scale_linewidth_manual(values = c("Neutral prior" = 0.6, "Summary likelihood" = 0.6, "Posterior" = 1.0),
                         name = NULL) +
  scale_x_continuous(breaks = seq(0.70, 1.30, 0.10), labels = sprintf("%.2f", seq(0.70, 1.30, 0.10)),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 6, 2), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Hazard ratio for all-cause death by day 28", y = "Density") +
  theme_ch55() +
  theme(legend.position = c(0.01, 0.99), legend.justification = c(0, 1),
        legend.key.width = unit(14, "pt"))

reg_lv <- c("Benefit (HR < 0.900)", "No clinically important effect (0.900–1.100)", "Harm (HR > 1.100)")
bayes_long <- cat3 %>%
  select(prior, benefit, equivalence, harm) %>%
  pivot_longer(-prior, names_to = "region_key", values_to = "probability") %>%
  mutate(region = factor(recode(region_key, benefit = reg_lv[1], equivalence = reg_lv[2], harm = reg_lv[3]),
                         levels = reg_lv),
         # display label only (UK spelling); the catalog key "Skeptical" is kept for CSV/numeric checks
         prior_label = factor(recode(prior, Skeptical = "Sceptical"), levels = c("Data-dominant", "Neutral", "Sceptical")),
         prior = factor(prior, levels = c("Data-dominant", "Neutral", "Skeptical")),
         label = percent(probability, accuracy = 0.1))

pb <- ggplot(bayes_long, aes(region, probability, fill = region)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = label), vjust = -0.35, size = 2.5, family = FONT, colour = COL_INK) +
  facet_wrap(~ prior_label, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.90),
                     breaks = seq(0, 0.75, 0.25), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = setNames(c(COL_BLUE, COL_LIGHT, COL_VERM), reg_lv)) +
  labs(x = NULL, y = "Posterior probability", fill = NULL) +
  theme_ch55() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.spacing.x = unit(12, "pt"))

p <- (pa / pb) + plot_layout(heights = c(1, 0.95)) +
  plot_annotation(title = "Figure S12", tag_levels = "a",
                  theme = theme_ch55() + theme(plot.margin = margin(4, 6, 2, 4))) &
  theme(plot.tag = element_text(face = "bold", size = BASE_SIZE + 2, family = FONT))

write_csv(bind_rows(
  density_df %>% transmute(panel = "a", series = as.character(curve), x_hr = hr, y = density),
  bayes_long %>% transmute(panel = "b", series = paste(prior, region_key, sep = " | "),
                           x_hr = NA_real_, y = probability)
), file.path(OUT_DIR, "FigS12_bayesian_three_priors_plotdata.csv"))
write_csv(cat3 %>% transmute(prior, prior_mean, prior_sd, posterior_location_hr = posterior_location,
                             posterior_sd_log = posterior_sd, benefit, equivalence, harm,
                             likelihood_theta_log = theta_hat, likelihood_se_log = se_hat),
          file.path(OUT_DIR, "FigS12_bayesian_three_priors_parameters.csv"))
save_fig(p, "FigS12_bayesian_three_priors", 6.0, 5.0)
