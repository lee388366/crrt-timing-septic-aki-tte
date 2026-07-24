# =============================================================================
# Script 22: Fix C1 (Bootstrap) + C3 (LRT→Wald) + C2 Sensitivity + C4 Check
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(DBI)
  library(RPostgres)
})

set.seed(20260412)

base_dir  <- "analysis"
data_dir  <- file.path(base_dir, "data")

# ── Load CCW data ─────────────────────────────────────────────────────────────
d <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))

cat(sprintf("CCW dataset: %d rows, %d unique stays\n",
            nrow(d), length(unique(d$stay_id))))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 1: FIX C1 — Bootstrap RMST with proper resampling
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("PART 1: FIX C1 — Corrected Bootstrap RMST\n")
cat(strrep("=", 60), "\n")

rmst_wkm <- function(times, events, weights, tau = 28) {
  fit <- survfit(Surv(times, events) ~ 1, weights = weights)
  t   <- c(0, fit$time)
  s   <- c(1, fit$surv)
  if (max(t) < tau) {
    t <- c(t, tau); s <- c(s, s[length(s)])
  } else {
    idx <- max(which(t <= tau))
    t <- c(t[seq_len(idx)], tau)
    s <- c(s[seq_len(idx)], s[idx])
  }
  dup <- duplicated(t); t <- t[!dup]; s <- s[!dup]
  sum(diff(t) * s[-length(s)])
}

# FIXED bootstrap: use merge to preserve duplicate draws
boot_rmst_diff_fixed <- function(data, B = 1000, tau = 28) {
  ids   <- unique(data$stay_id)
  N     <- length(ids)
  diffs <- numeric(B)
  for (i in seq_len(B)) {
    samp_df  <- data.frame(stay_id = sample(ids, N, replace = TRUE))
    boot_df  <- merge(samp_df, data, by = "stay_id", all.x = TRUE)
    dA <- boot_df[boot_df$strategy == "Early CRRT",    ]
    dB <- boot_df[boot_df$strategy == "Deferred CRRT", ]
    r_A <- tryCatch(rmst_wkm(dA$follow_days, dA$event, dA$ipcw_w_trim, tau),
                    error = function(e) NA)
    r_B <- tryCatch(rmst_wkm(dB$follow_days, dB$event, dB$ipcw_w_trim, tau),
                    error = function(e) NA)
    diffs[i] <- r_A - r_B
  }
  diffs <- diffs[!is.na(diffs)]
  c(diff = mean(diffs),
    lo   = unname(quantile(diffs, 0.025)),
    hi   = unname(quantile(diffs, 0.975)))
}

get_cox_res <- function(fit) {
  co  <- coef(fit)
  se  <- sqrt(diag(vcov(fit)))
  idx <- which(names(co) == "strategyEarly CRRT")
  hr  <- exp(co[idx])
  lo  <- exp(co[idx] - 1.96 * se[idx])
  hi  <- exp(co[idx] + 1.96 * se[idx])
  pv  <- 2 * pnorm(-abs(co[idx] / se[idx]))
  list(HR = hr, CI_lo = lo, CI_hi = hi, p = pv)
}

# ── Full cohort ─────────────────────────────────────────
cat("\n── Full cohort ──\n")
fit_full <- coxph(Surv(follow_days, event) ~ strategy,
                  data = d, weights = ipcw_w_trim,
                  cluster = stay_id, robust = TRUE)
cox_full <- get_cox_res(fit_full)
cat(sprintf("IPCW Cox HR: %.3f [%.3f, %.3f], p=%.4f\n",
            cox_full$HR, cox_full$CI_lo, cox_full$CI_hi, cox_full$p))

dA_full <- d[d$strategy == "Early CRRT",    ]
dB_full <- d[d$strategy == "Deferred CRRT", ]
rmst_A_full <- rmst_wkm(dA_full$follow_days, dA_full$event, dA_full$ipcw_w_trim)
rmst_B_full <- rmst_wkm(dB_full$follow_days, dB_full$event, dB_full$ipcw_w_trim)
cat(sprintf("RMST: Early=%.2f d, Deferred=%.2f d\n", rmst_A_full, rmst_B_full))

cat("Running FIXED Bootstrap RMST (full, B=1000)...\n")
boot_full <- boot_rmst_diff_fixed(d, B = 1000)
cat(sprintf("RMST diff = %.2f d [%.2f, %.2f]\n",
            boot_full["diff"], boot_full["lo"], boot_full["hi"]))

# ── Equipoise subgroup ──────────────────────────────────
cat("\n── Equipoise SOFA 8-14 ──\n")
d_eq <- d[d$equipoise == TRUE, ]
fit_eq <- coxph(Surv(follow_days, event) ~ strategy,
                data = d_eq, weights = ipcw_w_trim,
                cluster = stay_id, robust = TRUE)
cox_eq <- get_cox_res(fit_eq)
cat(sprintf("IPCW Cox HR: %.3f [%.3f, %.3f], p=%.4f\n",
            cox_eq$HR, cox_eq$CI_lo, cox_eq$CI_hi, cox_eq$p))

dA_eq <- d_eq[d_eq$strategy == "Early CRRT",    ]
dB_eq <- d_eq[d_eq$strategy == "Deferred CRRT", ]
rmst_A_eq <- rmst_wkm(dA_eq$follow_days, dA_eq$event, dA_eq$ipcw_w_trim)
rmst_B_eq <- rmst_wkm(dB_eq$follow_days, dB_eq$event, dB_eq$ipcw_w_trim)
cat(sprintf("RMST: Early=%.2f d, Deferred=%.2f d\n", rmst_A_eq, rmst_B_eq))

cat("Running FIXED Bootstrap RMST (equipoise, B=1000)...\n")
boot_eq <- boot_rmst_diff_fixed(d_eq, B = 1000)
cat(sprintf("RMST diff = %.2f d [%.2f, %.2f]\n",
            boot_eq["diff"], boot_eq["lo"], boot_eq["hi"]))

# Save C1 results
c1_results <- list(
  full_cohort = list(
    cox = cox_full,
    rmst_early = rmst_A_full, rmst_deferred = rmst_B_full,
    rmst_diff = boot_full["diff"],
    rmst_ci_lo = boot_full["lo"], rmst_ci_hi = boot_full["hi"]
  ),
  equipoise = list(
    cox = cox_eq,
    rmst_early = rmst_A_eq, rmst_deferred = rmst_B_eq,
    rmst_diff = boot_eq["diff"],
    rmst_ci_lo = boot_eq["lo"], rmst_ci_hi = boot_eq["hi"]
  )
)

cat("\n── C1 FIX: Comparison of old vs new RMST CIs ──\n")
old_results <- readRDS(file.path(data_dir, "main_results.rds"))
cat("OLD (buggy bootstrap) RMST CI:\n")
print(old_results[, c("Analysis", "RMST_diff", "CI_lo", "CI_hi")])
cat("\nNEW (corrected bootstrap) RMST CI:\n")
cat(sprintf("  Full cohort:   diff=%.3f [%.3f, %.3f]\n",
            boot_full["diff"], boot_full["lo"], boot_full["hi"]))
cat(sprintf("  Equipoise:     diff=%.3f [%.3f, %.3f]\n",
            boot_eq["diff"], boot_eq["lo"], boot_eq["hi"]))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 2: FIX C3 — Interaction test: use Wald instead of LRT
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("PART 2: FIX C3 — Interaction Test (Wald, not LRT)\n")
cat(strrep("=", 60), "\n")

d_int <- d %>%
  mutate(
    trt    = as.integer(strategy == "Early CRRT"),
    sofa_c = sofa_t0 - median(sofa_t0, na.rm = TRUE)
  ) %>%
  filter(!is.na(sofa_t0), !is.na(ipcw_w_trim))

d_int_eq <- d_int %>% filter(equipoise == TRUE)

run_interaction_fixed <- function(dat, label) {
  cat(sprintf("\n── %s ──\n", label))

  fit1 <- coxph(Surv(follow_days, event) ~ trt * sofa_c,
                data = dat, weights = ipcw_w_trim,
                cluster = stay_id, ties = "efron")

  s1 <- summary(fit1)
  cf <- s1$coefficients

  int_row <- grep("trt:sofa_c|sofa_c:trt", rownames(cf), value = TRUE)
  if (length(int_row) == 0) {
    cat("WARNING: interaction term not found\n")
    return(list(p_wald = NA))
  }

  hr_int <- exp(cf[int_row, "coef"])
  se_int <- cf[int_row, "robust se"]
  p_wald <- cf[int_row, "Pr(>|z|)"]

  ci_lo <- exp(cf[int_row, "coef"] - 1.96 * se_int)
  ci_hi <- exp(cf[int_row, "coef"] + 1.96 * se_int)

  cat(sprintf("  Interaction HR = %.4f [%.4f, %.4f]\n", hr_int, ci_lo, ci_hi))
  cat(sprintf("  Wald p = %.4f  (CORRECT — robust SE)\n", p_wald))

  list(hr_int = hr_int, ci_lo = ci_lo, ci_hi = ci_hi, p_wald = p_wald)
}

# Compare with old LRT results
old_int <- readRDS(file.path(data_dir, "interaction_sofa_results.rds"))
cat("\nOLD (invalid LRT):\n")
cat(sprintf("  Full cohort    LRT p = %.4f\n", old_int$full_cohort$p_int_lrt))
cat(sprintf("  Equipoise      LRT p = %.4f\n", old_int$equipoise_subgroup$p_int_lrt))
cat(sprintf("  Full cohort    Wald p = %.4f\n", old_int$full_cohort$p_int_wald))
cat(sprintf("  Equipoise      Wald p = %.4f\n", old_int$equipoise_subgroup$p_int_wald))

cat("\nNEW (corrected — Wald test with robust SE):\n")
int_full <- run_interaction_fixed(d_int,    "Full cohort")
int_eq   <- run_interaction_fixed(d_int_eq, "Equipoise SOFA 8-14")

c3_results <- list(
  full_cohort = int_full,
  equipoise   = int_eq,
  conclusion  = if (any(c(int_full$p_wald, int_eq$p_wald) < 0.10, na.rm = TRUE))
    "Evidence of effect modification" else "No evidence of effect modification"
)
cat(sprintf("\nConclusion: %s\n", c3_results$conclusion))

# ═══════════════════════════════════════════════════════════════════════════════
# PART 3: C2 SENSITIVITY — Strategy A with L₀-only IPCW denominator
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("PART 3: C2 SENSITIVITY — Strategy A L₀-only IPCW\n")
cat(strrep("=", 60), "\n")

baseline <- readRDS(file.path(data_dir, "baseline.rds"))

cohort_sens <- baseline %>%
  group_by(subject_id) %>%
  slice_min(icu_intime, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(sofa_t0))

L0_vars <- c("age", "sofa_t0", "apsiii", "creatinine_baseline",
             "nee_at_t0", "lactate_t0", "uo_rate_t0", "map_t0", "cci")

impute_median <- function(df, vars) {
  for (v in vars) {
    if (v %in% names(df)) {
      med <- median(df[[v]], na.rm = TRUE)
      df[[v]][is.na(df[[v]])] <- med
    }
  }
  df
}
cohort_sens <- impute_median(cohort_sens, L0_vars)

GRACE <- 24

ccw_sens <- cohort_sens %>%
  tidyr::crossing(strategy = c("A_early", "B_deferred")) %>%
  mutate(
    clone_id = paste(stay_id, strategy, sep = "_"),
    h_crrt   = crrt_hours_from_t0,
    protocol_censor_h = case_when(
      strategy == "A_early"    & (is.na(h_crrt) | h_crrt > GRACE) ~ as.double(GRACE),
      strategy == "B_deferred" & !is.na(h_crrt) & h_crrt <= GRACE ~ h_crrt,
      TRUE ~ NA_real_
    ),
    protocol_censored = !is.na(protocol_censor_h),
    follow_h    = pmin(follow_time_28d * 24,
                       coalesce(protocol_censor_h, follow_time_28d * 24)),
    follow_days = follow_h / 24,
    event = died_28d & (is.na(protocol_censor_h) |
                        follow_time_28d * 24 < protocol_censor_h)
  )

# Strategy A: L₀-only denominator (sensitivity)
A_sens <- ccw_sens %>%
  filter(strategy == "A_early") %>%
  mutate(y_A = as.integer(early_crrt))

# Denominator: L₀ only (NO time-varying L₂₃)
lm_denom_A_sens <- glm(
  y_A ~ age + sofa_t0 + apsiii + creatinine_baseline +
    nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci,
  data = A_sens, family = binomial()
)
p_denom_A_sens <- predict(lm_denom_A_sens, type = "response")

# Numerator: also L₀ only (same as original — stabilized weight)
p_num_A_sens <- p_denom_A_sens  # same model → weight ≈ 1

# With identical numerator and denominator, stabilized weights = 1
# So we use a simpler approach: numerator = marginal probability only
lm_num_A_simple <- glm(y_A ~ 1, data = A_sens, family = binomial())
p_num_A_simple  <- predict(lm_num_A_simple, type = "response")

A_sens <- A_sens %>%
  mutate(
    p_denom = p_denom_A_sens,
    p_num   = p_num_A_simple,
    ipcw_w  = if_else(early_crrt,
                      p_num   / pmax(p_denom,   0.01),
                      (1 - p_num) / pmax(1 - p_denom, 0.01))
  )

w99_A_sens <- quantile(A_sens$ipcw_w, 0.99, na.rm = TRUE)
A_sens <- A_sens %>%
  mutate(ipcw_w_trim = pmin(ipcw_w, w99_A_sens))

cat(sprintf("Strategy A L0-only weight summary:\n"))
cat(sprintf("  Mean: %.3f  Median: %.3f  Max: %.3f  99th: %.3f\n",
            mean(A_sens$ipcw_w_trim), median(A_sens$ipcw_w_trim),
            max(A_sens$ipcw_w_trim), w99_A_sens))

# Strategy B: keep original weights from ccw_final
ccw_orig <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)

B_orig_weights <- ccw_orig %>%
  filter(strategy == "B_deferred") %>%
  select(stay_id, ipcw_w_orig = ipcw_w, ipcw_w_trim_orig = ipcw_w_trim)

ccw_B_sens <- ccw_sens %>%
  filter(strategy == "B_deferred") %>%
  left_join(B_orig_weights, by = "stay_id") %>%
  mutate(
    ipcw_w      = coalesce(ipcw_w_orig, 1),
    ipcw_w_trim = coalesce(ipcw_w_trim_orig, 1)
  )

# Combine
ccw_c2 <- bind_rows(
  A_sens %>% select(any_of(names(ccw_B_sens))),
  ccw_B_sens %>% select(any_of(names(A_sens)))
)

ccw_c2 <- ccw_c2 %>%
  mutate(
    strategy_f = factor(strategy,
                        levels = c("B_deferred", "A_early"),
                        labels = c("Deferred CRRT", "Early CRRT")),
    equipoise = !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14
  )

# Cox model — full cohort
fit_c2_full <- coxph(Surv(follow_days, event) ~ strategy_f,
                     data = ccw_c2, weights = ipcw_w_trim,
                     cluster = stay_id, robust = TRUE)
cox_c2_full <- get_cox_res(fit_c2_full)

# Cox model — equipoise
ccw_c2_eq <- ccw_c2 %>% filter(equipoise)
fit_c2_eq <- coxph(Surv(follow_days, event) ~ strategy_f,
                   data = ccw_c2_eq, weights = ipcw_w_trim,
                   cluster = stay_id, robust = TRUE)
cox_c2_eq <- get_cox_res(fit_c2_eq)

cat("\n── C2 Sensitivity: Strategy A with L₀-only denominator ──\n")
cat(sprintf("  Full cohort:  HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_c2_full$HR, cox_c2_full$CI_lo, cox_c2_full$CI_hi, cox_c2_full$p))
cat(sprintf("  Equipoise:    HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_c2_eq$HR, cox_c2_eq$CI_lo, cox_c2_eq$CI_hi, cox_c2_eq$p))

cat("\n── Comparison: Original vs L₀-only ──\n")
cat(sprintf("  Original full:   HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_full$HR, cox_full$CI_lo, cox_full$CI_hi, cox_full$p))
cat(sprintf("  L0-only full:    HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_c2_full$HR, cox_c2_full$CI_lo, cox_c2_full$CI_hi, cox_c2_full$p))
cat(sprintf("  Original equip:  HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_eq$HR, cox_eq$CI_lo, cox_eq$CI_hi, cox_eq$p))
cat(sprintf("  L0-only equip:   HR=%.3f [%.3f, %.3f]  p=%.4f\n",
            cox_c2_eq$HR, cox_c2_eq$CI_lo, cox_c2_eq$CI_hi, cox_c2_eq$p))

hr_diff_full <- abs(cox_full$HR - cox_c2_full$HR)
hr_diff_eq   <- abs(cox_eq$HR - cox_c2_eq$HR)
cat(sprintf("\n  HR absolute difference (full): %.4f\n", hr_diff_full))
cat(sprintf("  HR absolute difference (equip): %.4f\n", hr_diff_eq))

if (hr_diff_full < 0.05 && hr_diff_eq < 0.05) {
  cat("\n  → CONCLUSION: L₀-only sensitivity yields very similar results.\n")
  cat("    The post-treatment covariate concern (C2) does NOT materially\n")
  cat("    affect the primary conclusions. No full re-run needed for C2.\n")
} else {
  cat("\n  → CONCLUSION: Meaningful difference detected. Consider full re-run.\n")
}

c2_results <- list(
  original_full = cox_full, original_eq = cox_eq,
  l0only_full = cox_c2_full, l0only_eq = cox_c2_eq,
  hr_diff_full = hr_diff_full, hr_diff_eq = hr_diff_eq
)

# ═══════════════════════════════════════════════════════════════════════════════
# PART 4: C4 CHECK — patients.dod coverage
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("PART 4: C4 CHECK — patients.dod coverage\n")
cat(strrep("=", 60), "\n")

con <- tryCatch(
  dbConnect(RPostgres::Postgres(),
            dbname = "mimiciv31", host = "localhost",
            port = 5432, user = "postgres", password = Sys.getenv("PGPASSWORD")),
  error = function(e) { cat("DB connection failed:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(con)) {
  c4_check <- dbGetQuery(con, "
    WITH cohort AS (
      SELECT subject_id, hadm_id, stay_id,
             t0, deathtime, follow_time_28d, died_28d
      FROM data_extract_crrt.cohort_base
    )
    SELECT
      COUNT(*)                                             AS total_patients,
      SUM(CASE WHEN p.dod IS NOT NULL THEN 1 ELSE 0 END)  AS has_dod,
      SUM(CASE WHEN c.deathtime IS NOT NULL THEN 1 ELSE 0 END) AS has_deathtime,
      -- Key: patients with dod but no deathtime (died out-of-hospital)
      SUM(CASE WHEN p.dod IS NOT NULL AND c.deathtime IS NULL THEN 1 ELSE 0 END)
        AS dod_only_no_deathtime,
      -- Among those: died within 28 days of t0
      SUM(CASE WHEN p.dod IS NOT NULL
                AND c.deathtime IS NULL
                AND p.dod <= (c.t0 + INTERVAL '28 day')::date
           THEN 1 ELSE 0 END)
        AS dod_within_28d_no_deathtime,
      -- Discharged alive but died before day 28
      SUM(CASE WHEN p.dod IS NOT NULL
                AND c.deathtime IS NULL
                AND p.dod <= (c.t0 + INTERVAL '28 day')::date
                AND c.died_28d = FALSE
           THEN 1 ELSE 0 END)
        AS missed_deaths_28d
    FROM cohort c
    JOIN mimiciv_hosp.patients p USING (subject_id)
  ")

  cat("\n── patients.dod Coverage Analysis ──\n")
  cat(sprintf("  Total patients in cohort:           %d\n", c4_check$total_patients))
  cat(sprintf("  Has dod (any death record):         %d (%.1f%%)\n",
              c4_check$has_dod, 100 * c4_check$has_dod / c4_check$total_patients))
  cat(sprintf("  Has deathtime (in-hospital):        %d (%.1f%%)\n",
              c4_check$has_deathtime, 100 * c4_check$has_deathtime / c4_check$total_patients))
  cat(sprintf("  dod only (no deathtime):            %d\n", c4_check$dod_only_no_deathtime))
  cat(sprintf("  dod within 28d (no deathtime):      %d\n", c4_check$dod_within_28d_no_deathtime))
  cat(sprintf("  MISSED deaths (alive per current):  %d\n", c4_check$missed_deaths_28d))
  cat(sprintf("  Impact: %.1f%% of cohort affected\n",
              100 * c4_check$missed_deaths_28d / c4_check$total_patients))

  if (c4_check$missed_deaths_28d > 0) {
    # Check by treatment group
    c4_by_group <- dbGetQuery(con, "
      WITH cohort AS (
        SELECT subject_id, hadm_id, stay_id, t0, deathtime,
               died_28d, early_crrt
        FROM data_extract_crrt.cohort_base
      )
      SELECT
        c.early_crrt,
        COUNT(*) AS n,
        SUM(CASE WHEN p.dod IS NOT NULL
                  AND c.deathtime IS NULL
                  AND p.dod <= (c.t0 + INTERVAL '28 day')::date
                  AND c.died_28d = FALSE
             THEN 1 ELSE 0 END) AS missed_deaths
      FROM cohort c
      JOIN mimiciv_hosp.patients p USING (subject_id)
      GROUP BY c.early_crrt
    ")
    cat("\n── Missed deaths by treatment group ──\n")
    print(c4_by_group)
  }

  dbDisconnect(con)
} else {
  cat("Skipping C4 check — database not available\n")
  c4_check <- NULL
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART 5: Save all results
# ═══════════════════════════════════════════════════════════════════════════════
all_fix_results <- list(
  c1_bootstrap = c1_results,
  c3_interaction = c3_results,
  c2_sensitivity = c2_results,
  c4_dod_check = c4_check,
  timestamp = Sys.time()
)

saveRDS(all_fix_results, file.path(data_dir, "fix_c1_c3_c2_c4_results.rds"))
cat("\n✅ All results saved to: analysis/data/fix_c1_c3_c2_c4_results.rds\n")

# ═══════════════════════════════════════════════════════════════════════════════
# GRAND SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("GRAND SUMMARY\n")
cat(strrep("=", 60), "\n")

cat("\nC1 (Bootstrap): FIXED — new RMST CIs computed with correct resampling\n")
cat(sprintf("  Equipoise RMST diff: %.2f [%.2f, %.2f]\n",
            boot_eq["diff"], boot_eq["lo"], boot_eq["hi"]))

cat(sprintf("\nC3 (Interaction): FIXED — Wald p used instead of invalid LRT\n"))
cat(sprintf("  Full cohort Wald p = %.4f\n", int_full$p_wald))
cat(sprintf("  Equipoise   Wald p = %.4f\n", int_eq$p_wald))

cat(sprintf("\nC2 (Strategy A IPCW): Sensitivity analysis completed\n"))
cat(sprintf("  HR difference (full): %.4f\n", hr_diff_full))
cat(sprintf("  HR difference (equip): %.4f\n", hr_diff_eq))

if (!is.null(c4_check)) {
  cat(sprintf("\nC4 (patients.dod): %d missed deaths (%.1f%% of cohort)\n",
              c4_check$missed_deaths_28d,
              100 * c4_check$missed_deaths_28d / c4_check$total_patients))
}
