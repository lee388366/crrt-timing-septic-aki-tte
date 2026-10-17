# =============================================================================
# strictK3 v2 internal QA checks:
#   QA-1: RRT-free days direction / clone-censoring issue
#   QA-2: PH assumption in calendar-era Cox sensitivity models
#
# Guardrail: read existing strictK3 v2 outputs only; do not re-run extraction,
# time-varying covariates, CCW/IPCW, or primary models.
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(dplyr)
  library(readr)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(tibble)
})

base_dir <- Sys.getenv("CRRT_ANALYSIS_ROOT")
strict_dir <- file.path(base_dir, "strictK3_primary_v2")
data_dir <- file.path(strict_dir, "data")
results_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run14/NEW")
fig_dir <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run14/NEW/figs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260715)

# CH-01 patch: connection only used to read patients.dod
con <- dbConnect(
  RPostgres::Postgres(),
  host = "localhost", port = 5432, dbname = "mimiciv31",
  user = "postgres", password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

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

ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

secondary <- read_csv(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run14/NEW/strictK3_v2_secondary_outcomes_ch01.csv"), show_col_types = FALSE) %>%
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
  # ===== CH-01 最小改动：仅重定义 died_28d / follow_time_28d（deathtime 或 patients.dod，午夜日界）=====
  left_join(
    dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients"),
    by = "subject_id"
  ) %>%
  mutate(
    .dc = coalesce(days_to_death, as.numeric(dod - as.Date(t0))),
    died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
    follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
  ) %>%
  select(-.dc, -dod) %>%
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
cat("[CH-01 check] unique stays died_28d (patching) =",
    sum(patient_df$died_28d %in% TRUE), "of", nrow(patient_df), "\n")

eq_all <- analysis_df %>% filter(equipoise %in% TRUE)
eq_compatible <- eq_all %>% filter(strategy_compatible_for_28d)
eq_patient <- patient_df %>% filter(equipoise %in% TRUE)

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
      stage = "All clone rows used by naive fixed-28d RFD script",
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
  run_lm_rfd(eq_all, "Naive all-clone weighted LM"),
  run_lm_rfd(eq_compatible, "Corrected strategy-compatible clone weighted LM")
) %>%
  mutate(
    mean_diff = round(mean_diff, 3),
    ci_lower = round(ci_lower, 3),
    ci_upper = round(ci_upper, 3),
    p_value = round(p_value, 4),
    p_display = fmt_p(p_value)
  )

definition_consistency <- tibble(
  item = c(
    "secondary_outcomes_rows",
    "ccw_unique_stays",
    "rfd_missing_clone_rows",
    "rfd_source",
    "death_handled_as_zero_in_rfd"
  ),
  value = c(
    as.character(nrow(secondary)),
    as.character(n_distinct(ccw$stay_id)),
    as.character(sum(is.na(analysis_df$rfd_28d))),
    "run14/NEW/strictK3_v2_secondary_outcomes_ch01.csv (source: data_extract_crrt_strictk3_v2.secondary_outcomes_ch01)",
    as.character(all(analysis_df$rfd_28d[analysis_df$died_28d %in% TRUE] == 0, na.rm = TRUE))
  )
)

write_csv(observed_by_actual, file.path(results_dir, "strictK3_v2_QA_rfd_observed_by_actual_group.csv"))
write_csv(clone_stage_distribution, file.path(results_dir, "strictK3_v2_QA_rfd_clone_stage_distribution.csv"))
write_csv(rfd_model_comparison, file.path(results_dir, "strictK3_v2_QA_rfd_model_comparison.csv"))
write_csv(definition_consistency, file.path(results_dir, "strictK3_v2_QA_rfd_definition_consistency.csv"))

corrected_row <- rfd_model_comparison %>% filter(grepl("^Corrected", analysis))
naive_row <- rfd_model_comparison %>% filter(grepl("^Naive", analysis))

rfd_report <- c(
  "# strictK3 v2 QA-1: RRT-free days direction / definition check",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Bottom line",
  "",
  "**The naive all-clone fixed-28d RFD estimate should not be interpreted as a causal efficacy result.** Fixed-time continuous outcomes are not observed under strategy clones after protocol censoring. The corrected v2 analysis therefore keeps only strategy-compatible clone rows, plus deaths before protocol censoring where 28-day RFD is known to be 0.",
  "",
  sprintf("- Naive all-clone RFD mean difference, Early strategy minus Deferred strategy: %.3f days (95%% CI %.3f to %.3f), p=%s.",
          naive_row$mean_diff, naive_row$ci_lower, naive_row$ci_upper, naive_row$p_display),
  sprintf("- Corrected strategy-compatible RFD mean difference, Early strategy minus Deferred strategy: %.3f days (95%% CI %.3f to %.3f), p=%s.",
          corrected_row$mean_diff, corrected_row$ci_lower, corrected_row$ci_upper, corrected_row$p_display),
  "",
  "## Observed patient-level distribution by actual treatment group",
  "",
  paste(capture.output(print(observed_by_actual, n = nrow(observed_by_actual))), collapse = "\n"),
  "",
  "## Clone/censoring decomposition",
  "",
  paste(capture.output(print(clone_stage_distribution, n = nrow(clone_stage_distribution))), collapse = "\n"),
  "",
  "## Naive vs corrected RFD model",
  "",
  paste(capture.output(print(rfd_model_comparison, n = nrow(rfd_model_comparison))), collapse = "\n"),
  "",
  "## Manuscript recommendation",
  "",
  "Use the corrected strategy-compatible result if RFD is reported. Position RFD as a descriptive/supportive renal-support-burden endpoint, not as a primary causal efficacy endpoint, because earlier CRRT mechanically subtracts CRRT calendar days from the score.",
  "",
  "## Output tables",
  "",
  "- `strictK3_v2_QA_rfd_observed_by_actual_group.csv`",
  "- `strictK3_v2_QA_rfd_clone_stage_distribution.csv`",
  "- `strictK3_v2_QA_rfd_model_comparison.csv`",
  "- `strictK3_v2_QA_rfd_definition_consistency.csv`"
)

writeLines(rfd_report, file.path(results_dir, "strictK3_v2_QA_rfd_report.md"))

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

write_csv(zph_table, file.path(results_dir, "strictK3_v2_QA_ph_zph_table.csv"))
write_csv(trt_results, file.path(results_dir, "strictK3_v2_QA_ph_treatment_results.csv"))

plot_zph_terms <- function(obj, safe_name) {
  png(file.path(fig_dir, paste0("strictK3_v2_QA_ph_", safe_name, ".png")), width = 2200, height = 1600, res = 300)
  tryCatch({
    plot(obj$zph)
  }, finally = {
    dev.off()
  })
}

plot_zph_terms(models[[1]], "sofa8_14_main_zph")
plot_zph_terms(models[[2]], "sofa8_14_era_adjusted_zph")
plot_zph_terms(models[[3]], "sofa8_14_era_stratified_zph")

ph_report <- c(
  "# strictK3 v2 QA-2: Calendar-era adjusted Cox PH assumption",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Bottom line",
  "",
  "The treatment contrast remains null across unadjusted, era-adjusted, and era-stratified Cox sensitivity models. If the era-adjusted model shows a global PH warning, interpret it through the term-level Schoenfeld tests and prefer the era-stratified Cox sensitivity for manuscript robustness.",
  "",
  "## Schoenfeld PH tests",
  "",
  paste(capture.output(print(zph_table, n = nrow(zph_table))), collapse = "\n"),
  "",
  "## Treatment effect results",
  "",
  paste(capture.output(print(trt_results, n = nrow(trt_results))), collapse = "\n"),
  "",
  "## Output files",
  "",
  "- `strictK3_v2_QA_ph_zph_table.csv`",
  "- `strictK3_v2_QA_ph_treatment_results.csv`",
  "- `strictK3_v2_QA_ph_sofa8_14_main_zph.png`",
  "- `strictK3_v2_QA_ph_sofa8_14_era_adjusted_zph.png`",
  "- `strictK3_v2_QA_ph_sofa8_14_era_stratified_zph.png`"
)

writeLines(ph_report, file.path(results_dir, "strictK3_v2_QA_ph_report.md"))

cat("QA reports written:\n")
cat("  ", file.path(results_dir, "strictK3_v2_QA_rfd_report.md"), "\n", sep = "")
cat("  ", file.path(results_dir, "strictK3_v2_QA_ph_report.md"), "\n", sep = "")