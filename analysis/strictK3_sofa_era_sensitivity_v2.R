# =============================================================================
# Strict clean first-KDIGO3: SOFA strata, SOFA band scan, and calendar-era sensitivity
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(tibble)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
fig_dir <- file.path(strict_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    sofa_stratum = case_when(
      sofa_t0 <= 7 ~ "SOFA <=7",
      sofa_t0 >= 8 & sofa_t0 <= 14 ~ "SOFA 8-14",
      sofa_t0 >= 15 ~ "SOFA >=15",
      TRUE ~ NA_character_
    )
  )

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

run_weighted_cox <- function(dat, label) {
  dat <- dat %>%
    filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim), !is.na(trt))

  if (n_distinct(dat$stay_id) < 50 || length(unique(dat$trt)) < 2) {
    return(tibble(
      analysis = label,
      n_patients = n_distinct(dat$stay_id),
      clone_rows = nrow(dat),
      early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
      events = sum(dat$event, na.rm = TRUE),
      weighted_hr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      ph_global_p = NA_real_
    ))
  }

  fit <- coxph(
    Surv(follow_days, event) ~ trt,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  s <- summary(fit)
  ph_p <- tryCatch(cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_)

  tibble(
    analysis = label,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat),
    early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
    events = sum(dat$event, na.rm = TRUE),
    weighted_hr = unname(s$coefficients["trt", "exp(coef)"]),
    ci_lower = unname(s$conf.int["trt", "lower .95"]),
    ci_upper = unname(s$conf.int["trt", "upper .95"]),
    p_value = unname(s$coefficients["trt", "Pr(>|z|)"]),
    ph_global_p = ph_p
  )
}

interaction_p <- function(dat, group_var) {
  dat <- dat %>%
    filter(!is.na(.data[[group_var]]), !is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim)) %>%
    mutate(group_tmp = factor(.data[[group_var]]))
  if (n_distinct(dat$group_tmp) < 2) return(NA_real_)
  fit1 <- coxph(
    Surv(follow_days, event) ~ trt * group_tmp,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  beta <- coef(fit1)
  interaction_terms <- grep("^trt:group_tmp", names(beta), value = TRUE)
  if (length(interaction_terms) == 0) return(NA_real_)
  V <- vcov(fit1)[interaction_terms, interaction_terms, drop = FALSE]
  b <- beta[interaction_terms]
  tryCatch({
    chisq <- as.numeric(t(b) %*% solve(V) %*% b)
    pchisq(chisq, df = length(b), lower.tail = FALSE)
  }, error = function(e) NA_real_)
}

run_era_adjusted_cox <- function(dat, label) {
  dat <- dat %>%
    filter(
      !is.na(follow_days),
      !is.na(event),
      !is.na(ipcw_w_trim),
      !is.na(trt),
      !is.na(anchor_year_group)
    ) %>%
    mutate(anchor_year_group = factor(anchor_year_group))

  fit <- coxph(
    Surv(follow_days, event) ~ trt + anchor_year_group,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  s <- summary(fit)
  ph_p <- tryCatch(cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_)

  tibble(
    analysis = label,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat),
    early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
    events = sum(dat$event, na.rm = TRUE),
    weighted_hr = unname(s$coefficients["trt", "exp(coef)"]),
    ci_lower = unname(s$conf.int["trt", "lower .95"]),
    ci_upper = unname(s$conf.int["trt", "upper .95"]),
    p_value = unname(s$coefficients["trt", "Pr(>|z|)"]),
    ph_global_p = ph_p
  )
}

strata_results <- bind_rows(
  run_weighted_cox(filter(ccw, sofa_t0 <= 7), "SOFA <=7"),
  run_weighted_cox(filter(ccw, sofa_t0 >= 8 & sofa_t0 <= 14), "SOFA 8-14"),
  run_weighted_cox(filter(ccw, sofa_t0 >= 9 & sofa_t0 <= 13), "SOFA 9-13"),
  run_weighted_cox(filter(ccw, sofa_t0 >= 15), "SOFA >=15")
)

sofa_interaction_p <- interaction_p(
  ccw %>% filter(!is.na(sofa_stratum)),
  "sofa_stratum"
)

scan_grid <- expand.grid(
  lower = 6:10,
  upper = 12:16
) %>%
  as_tibble() %>%
  filter(lower < upper) %>%
  mutate(label = paste0("SOFA ", lower, "-", upper))

scan_results <- bind_rows(lapply(seq_len(nrow(scan_grid)), function(i) {
  lo <- scan_grid$lower[i]
  hi <- scan_grid$upper[i]
  run_weighted_cox(
    filter(ccw, sofa_t0 >= lo & sofa_t0 <= hi),
    scan_grid$label[i]
  ) %>%
    mutate(lower = lo, upper = hi, .before = analysis)
}))

era_results <- ccw %>%
  filter(equipoise %in% TRUE, !is.na(anchor_year_group)) %>%
  group_split(anchor_year_group) %>%
  lapply(function(x) {
    run_weighted_cox(x, paste0("SOFA 8-14 era ", unique(x$anchor_year_group)))
  }) %>%
  bind_rows() %>%
  mutate(anchor_year_group = sub("^SOFA 8-14 era ", "", analysis), .before = analysis)

era_interaction_p <- interaction_p(
  ccw %>% filter(equipoise %in% TRUE, !is.na(anchor_year_group)),
  "anchor_year_group"
)

era_adjusted_results <- bind_rows(
  run_era_adjusted_cox(ccw, "Full cohort adjusted for anchor_year_group"),
  run_era_adjusted_cox(filter(ccw, equipoise %in% TRUE), "SOFA 8-14 adjusted for anchor_year_group")
)

write_csv(strata_results, file.path(results_dir, "strictK3_v2_sofa_strata_results.csv"))
write_csv(scan_results, file.path(results_dir, "strictK3_v2_sofa_band_scan_results.csv"))
write_csv(era_results, file.path(results_dir, "strictK3_v2_calendar_era_results.csv"))
write_csv(era_adjusted_results, file.path(results_dir, "strictK3_v2_calendar_era_adjusted_results.csv"))
saveRDS(
  list(
    sofa_strata = strata_results,
    sofa_scan = scan_results,
    calendar_era = era_results,
    calendar_era_adjusted = era_adjusted_results,
    sofa_interaction_p = sofa_interaction_p,
    era_interaction_p = era_interaction_p
  ),
  file.path(results_dir, "strictK3_v2_sofa_era_sensitivity_results.rds")
)

scan_plot <- scan_results %>%
  mutate(
    lower = factor(lower),
    upper = factor(upper),
    sig = !is.na(p_value) & p_value < 0.05
  )

p <- ggplot(scan_plot, aes(x = upper, y = lower, fill = weighted_hr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(weighted_hr), "NA", sprintf("%.2f\nn=%d", weighted_hr, n_patients))), size = 3) +
  scale_fill_gradient2(
    low = "#2E86C1",
    mid = "white",
    high = "#C0392B",
    midpoint = 1,
    na.value = "grey85"
  ) +
  labs(
    title = "Strict clean first-KDIGO3 SOFA band scan",
    subtitle = "IPCW-weighted Cox HR for Early vs Deferred CRRT; each tile is SOFA lower-upper band",
    x = "Upper SOFA cutpoint",
    y = "Lower SOFA cutpoint",
    fill = "HR"
  ) +
  theme_minimal(base_size = 12)

png(file.path(fig_dir, "strictK3_v2_sofa_band_scan_heatmap.png"), width = 2100, height = 1650, res = 300)
print(p)
dev.off()

pdf(file.path(fig_dir, "strictK3_v2_sofa_band_scan_heatmap.pdf"), width = 7, height = 5.5, onefile = FALSE)
print(p)
dev.off()

fmt_table <- function(df, include_prefix_cols = character()) {
  df %>%
    mutate(
      hr_ci = ifelse(is.na(weighted_hr), "NA", sprintf("%.3f (%.3f to %.3f)", weighted_hr, ci_lower, ci_upper)),
      p_text = fmt_p(p_value)
    ) %>%
    select(any_of(include_prefix_cols), analysis, n_patients, early_n, events, hr_ci, p_text, ph_global_p)
}

report <- c(
  "# Strict clean first-KDIGO3 SOFA and calendar-era sensitivity",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## SOFA strata",
  "",
  paste(capture.output(print(fmt_table(strata_results), n = nrow(strata_results))), collapse = "\n"),
  "",
  paste0("SOFA stratum interaction p: ", fmt_p(sofa_interaction_p)),
  "",
  "## Calendar-era strata within primary SOFA 8-14",
  "",
  paste(capture.output(print(fmt_table(era_results, "anchor_year_group"), n = nrow(era_results))), collapse = "\n"),
  "",
  paste0("Calendar-era interaction p: ", fmt_p(era_interaction_p)),
  "",
  "## Calendar-era adjusted Cox models",
  "",
  paste(capture.output(print(fmt_table(era_adjusted_results), n = nrow(era_adjusted_results))), collapse = "\n"),
  "",
  "## SOFA band scan",
  "",
  paste(capture.output(print(fmt_table(scan_results, c("lower", "upper")), n = nrow(scan_results))), collapse = "\n"),
  "",
  "Outputs:",
  "",
  paste0("- `", file.path(results_dir, "strictK3_v2_sofa_strata_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_sofa_band_scan_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_calendar_era_results.csv"), "`"),
  paste0("- `", file.path(results_dir, "strictK3_v2_calendar_era_adjusted_results.csv"), "`"),
  paste0("- `", file.path(fig_dir, "strictK3_v2_sofa_band_scan_heatmap.png"), "`")
)

writeLines(report, file.path(results_dir, "strictK3_v2_sofa_era_sensitivity_report.md"))

cat("\n=== SOFA strata ===\n")
print(strata_results, n = nrow(strata_results))
cat("\nSOFA interaction p:", sofa_interaction_p, "\n")
cat("\n=== Calendar era ===\n")
print(era_results, n = nrow(era_results))
cat("\nEra interaction p:", era_interaction_p, "\n")
cat("\n=== Calendar-era adjusted ===\n")
print(era_adjusted_results, n = nrow(era_adjusted_results))
cat("\nSaved report:", file.path(results_dir, "strictK3_v2_sofa_era_sensitivity_report.md"), "\n")
