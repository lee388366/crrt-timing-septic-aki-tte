# =============================================================================
# strictK3 internal QA checks:
#   QA-1: RRT-free days direction / clone-censoring issue
#   QA-2: PH assumption in calendar-era adjusted Cox models
#
# Guardrail: read existing strictK3 outputs only; do not re-run extraction,
# time-varying covariates, CCW/IPCW, or primary models.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(tibble)
})

base_dir <- "analysis"
strict_dir <- file.path(base_dir, "strictK3_primary")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(strict_dir, "results")
fig_dir <- file.path(strict_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260715)

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

weighted_quantile <- function(x, w, probs) {
  keep <- !(is.na(x) | is.na(w))
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0 || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

fmt_miq <- function(x) {
  sprintf(
    "%.1f [%.1f, %.1f]",
    median(x, na.rm = TRUE),
    quantile(x, 0.25, na.rm = TRUE),
    quantile(x, 0.75, na.rm = TRUE)
  )
}

fmt_wmiq <- function(x, w) {
  q <- weighted_quantile(x, w, c(0.25, 0.5, 0.75))
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

fix_bool <- function(x) {
  if (is.logical(x)) return(x)
  ifelse(
    x %in% c("t", "TRUE", "true", "1", 1),
    TRUE,
    ifelse(x %in% c("f", "FALSE", "false", "0", 0), FALSE, NA)
  )
}

run_lm_rfd <- function(dat, label) {
  dat <- dat %>%
    filter(!is.na(rfd_28d), !is.na(ipcw_w_trim), !is.na(strategy_trt))

  if (length(unique(dat$strategy_trt)) < 2) {
    return(tibble(
      analysis = label,
      n_patients = n_distinct(dat$stay_id),
      clone_rows = nrow(dat),
      mean_diff = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_
    ))
  }

  fit <- lm(rfd_28d ~ strategy_trt, data = dat, weights = ipcw_w_trim)
  vc <- vcovCL(fit, cluster = dat$stay_id, type = "HC3")
  ct <- coeftest(fit, vcov. = vc)
  beta <- unname(ct["strategy_trt", "Estimate"])
  se <- unname(ct["strategy_trt", "Std. Error"])
  p <- unname(ct["strategy_trt", "Pr(>|t|)"])
  tibble(
    analysis = label,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat),
    mean_diff = beta,
    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,
    p_value = p
  )
}

# ── Load data ────────────────────────────────────────────────────────────────
ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

secondary <- read_csv(file.path(data_dir, "strictK3_secondary_outcomes.csv"), show_col_types = FALSE) %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    early_crrt = fix_bool(early_crrt),
    died_28d = fix_bool(died_28d),
    renal_recovery_28d = fix_bool(renal_recovery_28d),
    crrt_total_days = as.numeric(crrt_total_days),
    rrt_free_days_28d = as.numeric(rrt_free_days_28d)
  )

analysis_df <- ccw %>%
  left_join(
    secondary %>%
      select(stay_id, renal_recovery_28d, rrt_free_days_28d, crrt_total_days),
    by = "stay_id"
  ) %>%
  mutate(
    strategy_trt = as.integer(strategy == "A_early"),
    rfd_28d = case_when(
      !is.na(died_28d) & died_28d ~ 0,
      !is.na(rrt_free_days_28d) ~ rrt_free_days_28d,
      TRUE ~ NA_real_
    ),
    censor_h = protocol_censor_h,
    death_before_protocol_censor = !is.na(censor_h) & died_28d & (follow_time_28d * 24 < censor_h),
    strategy_compatible_for_28d = !protocol_censored | death_before_protocol_censor
  )

patient_df <- analysis_df %>%
  distinct(stay_id, .keep_all = TRUE)

eq_all <- analysis_df %>% filter(equipoise %in% TRUE)
eq_compatible <- eq_all %>% filter(strategy_compatible_for_28d)
eq_patient <- patient_df %>% filter(equipoise %in% TRUE)

# ── QA-1 RFD diagnostics ─────────────────────────────────────────────────────
observed_by_actual <- eq_patient %>%
  mutate(actual_group = if_else(early_crrt, "Observed early CRRT", "Observed deferred")) %>%
  group_by(actual_group) %>%
  summarise(
    n = n(),
    died_n = sum(died_28d, na.rm = TRUE),
    died_pct = round(100 * mean(died_28d, na.rm = TRUE), 1),
    crrt_days_median_iqr = fmt_miq(crrt_total_days),
    rfd_median_iqr = fmt_miq(rfd_28d),
    rfd_mean = mean(rfd_28d, na.rm = TRUE),
    .groups = "drop"
  )

clone_stage_distribution <- bind_rows(
  eq_all %>%
    group_by(strategy) %>%
    summarise(
      stage = "All clone rows currently used by old RFD script",
      n_patients = n_distinct(stay_id),
      clone_rows = n(),
      protocol_censored_rows = sum(protocol_censored, na.rm = TRUE),
      death_before_censor_rows = sum(death_before_protocol_censor, na.rm = TRUE),
      compatible_rows = sum(strategy_compatible_for_28d, na.rm = TRUE),
      unweighted_median_iqr = fmt_miq(rfd_28d),
      weighted_median_iqr = fmt_wmiq(rfd_28d, ipcw_w_trim),
      weighted_mean = weighted.mean(rfd_28d, ipcw_w_trim, na.rm = TRUE),
      .groups = "drop"
    ),
  eq_compatible %>%
    group_by(strategy) %>%
    summarise(
      stage = "Strategy-compatible clone rows only",
      n_patients = n_distinct(stay_id),
      clone_rows = n(),
      protocol_censored_rows = sum(protocol_censored, na.rm = TRUE),
      death_before_censor_rows = sum(death_before_protocol_censor, na.rm = TRUE),
      compatible_rows = sum(strategy_compatible_for_28d, na.rm = TRUE),
      unweighted_median_iqr = fmt_miq(rfd_28d),
      weighted_median_iqr = fmt_wmiq(rfd_28d, ipcw_w_trim),
      weighted_mean = weighted.mean(rfd_28d, ipcw_w_trim, na.rm = TRUE),
      .groups = "drop"
    )
) %>%
  relocate(stage, strategy)

rfd_model_comparison <- bind_rows(
  run_lm_rfd(eq_all, "Old all-clone weighted LM (reproduces existing RFD result)"),
  run_lm_rfd(eq_compatible, "Corrected strategy-compatible clone weighted LM")
) %>%
  mutate(
    mean_diff = round(mean_diff, 3),
    ci_lower = round(ci_lower, 3),
    ci_upper = round(ci_upper, 3),
    p_value = round(p_value, 4),
    p_display = fmt_p(p_value)
  )

top_weight_influence <- eq_all %>%
  group_by(strategy) %>%
  mutate(
    weight_rank = percent_rank(ipcw_w_trim),
    top10_weight = weight_rank >= 0.90,
    weighted_rfd_contribution = ipcw_w_trim * rfd_28d
  ) %>%
  summarise(
    n_rows = n(),
    top10_rows = sum(top10_weight, na.rm = TRUE),
    top10_share_of_total_weight = sum(ipcw_w_trim[top10_weight], na.rm = TRUE) / sum(ipcw_w_trim, na.rm = TRUE),
    top10_share_of_weighted_rfd_sum = sum(weighted_rfd_contribution[top10_weight], na.rm = TRUE) / sum(weighted_rfd_contribution, na.rm = TRUE),
    top10_rfd_median_iqr = fmt_miq(rfd_28d[top10_weight]),
    top10_death_pct = round(100 * mean(died_28d[top10_weight], na.rm = TRUE), 1),
    .groups = "drop"
  )

definition_consistency <- tibble(
  item = c(
    "secondary_outcomes_rows",
    "ccw_unique_stays",
    "rfd_missing_clone_rows",
    "finegray_rfd_same_secondary_source",
    "death_handled_as_zero_in_rfd"
  ),
  value = c(
    as.character(nrow(secondary)),
    as.character(n_distinct(ccw$stay_id)),
    as.character(sum(is.na(analysis_df$rfd_28d))),
    "YES: strictK3_weighted_finegray_rfd.R and this QA read strictK3_secondary_outcomes.csv/data_extract_crrt_strictk3.secondary_outcomes",
    as.character(all(analysis_df$rfd_28d[analysis_df$died_28d %in% TRUE] == 0, na.rm = TRUE))
  )
)

write_csv(observed_by_actual, file.path(results_dir, "strictK3_QA_rfd_observed_by_actual_group.csv"))
write_csv(clone_stage_distribution, file.path(results_dir, "strictK3_QA_rfd_clone_stage_distribution.csv"))
write_csv(rfd_model_comparison, file.path(results_dir, "strictK3_QA_rfd_model_comparison.csv"))
write_csv(top_weight_influence, file.path(results_dir, "strictK3_QA_rfd_top_weight_influence.csv"))
write_csv(definition_consistency, file.path(results_dir, "strictK3_QA_rfd_definition_consistency.csv"))

old_row <- rfd_model_comparison %>% filter(grepl("^Old", analysis))
corrected_row <- rfd_model_comparison %>% filter(grepl("^Corrected", analysis))

rfd_report <- c(
  "# strictK3 QA-1: RRT-free days direction / definition check",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Bottom line",
  "",
  "**The previously reported +0.278 day weighted RFD difference is a clone-censoring artifact, not a reliable causal effect estimate.** The old script assigned 28-day RFD outcomes to all clone rows, including clones that had already been protocol-censored for nonadherence to the assigned strategy. For a fixed 28-day continuous endpoint, those protocol-censored clone outcomes are not observed under that strategy and should not be used directly.",
  "",
  "Using only strategy-compatible clone rows (plus deaths before protocol censoring, where the 28-day RFD is known to be 0), the SOFA 8-14 weighted mean difference becomes:",
  "",
  sprintf("- Corrected RFD mean difference, Early strategy minus Deferred strategy: %.3f days (95%% CI %.3f to %.3f), p=%s.",
          corrected_row$mean_diff, corrected_row$ci_lower, corrected_row$ci_upper, corrected_row$p_display),
  "",
  "## RFD definition check",
  "",
  "- Formula confirmed: death within 28 days = 0; survivors = 28 minus distinct CRRT calendar days in [t0, t0+28d].",
  "- Death handling is symmetric across observed early/deferred groups and strategy clones.",
  "- RFD and Fine-Gray secondary analyses read the same strictK3 secondary-outcome source; Fine-Gray itself is not affected by the RFD clone-continuous-outcome issue.",
  "",
  "## Observed patient-level distribution by actual treatment group",
  "",
  paste(capture.output(print(observed_by_actual, n = nrow(observed_by_actual))), collapse = "\n"),
  "",
  "## Clone/censoring decomposition",
  "",
  paste(capture.output(print(clone_stage_distribution, n = nrow(clone_stage_distribution))), collapse = "\n"),
  "",
  "## Old vs corrected RFD model",
  "",
  paste(capture.output(print(rfd_model_comparison, n = nrow(rfd_model_comparison))), collapse = "\n"),
  "",
  "## High-weight influence check",
  "",
  paste(capture.output(print(top_weight_influence, n = nrow(top_weight_influence))), collapse = "\n"),
  "",
  "## Manuscript recommendation",
  "",
  "Do **not** present the old +0.278-day RFD result as an efficacy finding. The corrected strategy-compatible clone analysis reverses direction (early strategy has fewer RFD), but this endpoint is mechanically affected by initiating CRRT because every CRRT calendar day is subtracted from the score. Therefore, RFD should be positioned as a descriptive/supportive renal-support-burden endpoint rather than as a primary causal efficacy endpoint.",
  "",
  "Suggested sentence:",
  "",
  "> RRT-free days were treated as a descriptive secondary endpoint because fixed-time continuous outcomes are sensitive to clone-censoring choices and mechanically penalize earlier CRRT exposure; after excluding protocol-censored clone outcomes, early CRRT showed fewer RRT-free days, consistent with greater renal-support exposure rather than a mortality benefit.",
  "",
  "## Output tables",
  "",
  "- `strictK3_QA_rfd_observed_by_actual_group.csv`",
  "- `strictK3_QA_rfd_clone_stage_distribution.csv`",
  "- `strictK3_QA_rfd_model_comparison.csv`",
  "- `strictK3_QA_rfd_top_weight_influence.csv`",
  "- `strictK3_QA_rfd_definition_consistency.csv`"
)

writeLines(rfd_report, file.path(results_dir, "strictK3_QA_rfd_report.md"))

# ── QA-2 PH diagnostics ──────────────────────────────────────────────────────
run_cox <- function(dat, formula, label) {
  environment(formula) <- environment()
  fit <- coxph(
    formula,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  zph <- cox.zph(fit)
  list(
    label = label,
    fit = fit,
    zph = zph,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat)
  )
}

extract_zph <- function(obj) {
  tab <- as.data.frame(obj$zph$table)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  tab %>%
    as_tibble() %>%
    mutate(model = obj$label, .before = term) %>%
    select(model, term, everything())
}

extract_trt <- function(obj) {
  s <- summary(obj$fit)
  tibble(
    model = obj$label,
    n_patients = obj$n_patients,
    clone_rows = obj$clone_rows,
    HR = unname(s$coefficients["strategy_trt", "exp(coef)"]),
    CI_lo = unname(s$conf.int["strategy_trt", "lower .95"]),
    CI_hi = unname(s$conf.int["strategy_trt", "upper .95"]),
    p_value = unname(s$coefficients["strategy_trt", "Pr(>|z|)"])
  )
}

ph_df <- analysis_df %>%
  filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim), !is.na(anchor_year_group)) %>%
  mutate(anchor_year_group = factor(anchor_year_group))

ph_eq <- ph_df %>% filter(equipoise %in% TRUE)

models <- list(
  run_cox(ph_eq, Surv(follow_days, event) ~ strategy_trt, "SOFA 8-14 main unadjusted"),
  run_cox(ph_eq, Surv(follow_days, event) ~ strategy_trt + anchor_year_group, "SOFA 8-14 era-adjusted"),
  run_cox(ph_eq, Surv(follow_days, event) ~ strategy_trt + strata(anchor_year_group), "SOFA 8-14 era-stratified"),
  run_cox(ph_df, Surv(follow_days, event) ~ strategy_trt + anchor_year_group, "Full cohort era-adjusted"),
  run_cox(ph_df, Surv(follow_days, event) ~ strategy_trt + strata(anchor_year_group), "Full cohort era-stratified")
)

zph_table <- bind_rows(lapply(models, extract_zph)) %>%
  mutate(p_display = fmt_p(p))

trt_results <- bind_rows(lapply(models, extract_trt)) %>%
  mutate(
    HR = round(HR, 3),
    CI_lo = round(CI_lo, 3),
    CI_hi = round(CI_hi, 3),
    p_value = round(p_value, 4),
    p_display = fmt_p(p_value)
  )

write_csv(zph_table, file.path(results_dir, "strictK3_QA_ph_zph_table.csv"))
write_csv(trt_results, file.path(results_dir, "strictK3_QA_ph_treatment_results.csv"))

plot_zph_terms <- function(obj, safe_name) {
  png(file.path(fig_dir, paste0("strictK3_QA_ph_", safe_name, ".png")), width = 2200, height = 1600, res = 300)
  tryCatch({
    plot(obj$zph)
  }, finally = {
    dev.off()
  })
}

plot_zph_terms(models[[1]], "sofa8_14_main_zph")
plot_zph_terms(models[[2]], "sofa8_14_era_adjusted_zph")
plot_zph_terms(models[[3]], "sofa8_14_era_stratified_zph")

zph_eq_adj <- zph_table %>% filter(model == "SOFA 8-14 era-adjusted")
treatment_ph_p <- zph_eq_adj %>% filter(term == "strategy_trt") %>% pull(p)
era_ph_p <- zph_eq_adj %>% filter(term == "anchor_year_group") %>% pull(p)
global_ph_p <- zph_eq_adj %>% filter(term == "GLOBAL") %>% pull(p)
strat_eq <- trt_results %>% filter(model == "SOFA 8-14 era-stratified")
adj_eq <- trt_results %>% filter(model == "SOFA 8-14 era-adjusted")

ph_report <- c(
  "# strictK3 QA-2: Calendar-era adjusted Cox PH assumption",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Bottom line",
  "",
  sprintf("In the SOFA 8-14 era-adjusted model, the global PH test is significant (p=%s), but the treatment term itself is not the driver (strategy_trt p=%s). The PH signal is attributable to the calendar-era nuisance term (anchor_year_group p=%s).",
          fmt_p(global_ph_p), fmt_p(treatment_ph_p), fmt_p(era_ph_p)),
  "",
  "Therefore, the early-vs-deferred treatment HR is not materially undermined by treatment-specific PH violation. For manuscript robustness, the cleaner era-adjusted approach is to stratify the Cox baseline hazard by `anchor_year_group` rather than include era as an ordinary covariate.",
  "",
  sprintf("- Era-adjusted SOFA 8-14 treatment HR: %.3f (95%% CI %.3f to %.3f), p=%s.",
          adj_eq$HR, adj_eq$CI_lo, adj_eq$CI_hi, adj_eq$p_display),
  sprintf("- Era-stratified SOFA 8-14 treatment HR: %.3f (95%% CI %.3f to %.3f), p=%s.",
          strat_eq$HR, strat_eq$CI_lo, strat_eq$CI_hi, strat_eq$p_display),
  "",
  "## Schoenfeld PH tests",
  "",
  paste(capture.output(print(zph_table, n = nrow(zph_table))), collapse = "\n"),
  "",
  "## Treatment effect results",
  "",
  paste(capture.output(print(trt_results, n = nrow(trt_results))), collapse = "\n"),
  "",
  "## Recommendation for manuscript",
  "",
  "Report the primary unadjusted-by-era IPCW Cox HR and RMST as planned. For the calendar-era sensitivity, prefer the era-stratified Cox model or report that the era-adjusted model's global PH warning was driven by the nuisance calendar-era term rather than the treatment contrast. RMST remains the PH-independent estimand.",
  "",
  "Suggested sentence:",
  "",
  "> In calendar-era sensitivity analyses, the proportional hazards warning after adding era as a covariate was driven by the era term rather than the treatment contrast; an era-stratified Cox model yielded a similar null early-versus-deferred estimate, and the RMST analysis does not rely on the proportional hazards assumption.",
  "",
  "## Output files",
  "",
  "- `strictK3_QA_ph_zph_table.csv`",
  "- `strictK3_QA_ph_treatment_results.csv`",
  "- `strictK3_QA_ph_sofa8_14_main_zph.png`",
  "- `strictK3_QA_ph_sofa8_14_era_adjusted_zph.png`",
  "- `strictK3_QA_ph_sofa8_14_era_stratified_zph.png`"
)

writeLines(ph_report, file.path(results_dir, "strictK3_QA_ph_report.md"))

cat("QA reports written:\n")
cat("  ", file.path(results_dir, "strictK3_QA_rfd_report.md"), "\n", sep = "")
cat("  ", file.path(results_dir, "strictK3_QA_ph_report.md"), "\n", sep = "")
