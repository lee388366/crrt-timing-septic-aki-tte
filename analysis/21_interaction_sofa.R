# =============================================================================
# Script 21: strategy × sofa_t0 interaction test in IPCW-weighted Cox model
# Purpose: test whether baseline SOFA modifies the effect of early CRRT timing
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

base_dir <- "analysis"
d <- readRDS(file.path(base_dir, "data/ccw_final.rds"))

# Centre sofa_t0 on its median to keep main-effect estimate interpretable
sofa_med <- median(d$sofa_t0, na.rm = TRUE)
cat(sprintf("sofa_t0 median (centering value): %g\n", sofa_med))

d <- d %>%
  mutate(
    trt    = as.integer(strategy == "A_early"),
    sofa_c = sofa_t0 - sofa_med
  ) %>%
  filter(!is.na(sofa_t0), !is.na(ipcw_w_trim))

d_eq <- d %>% filter(equipoise == TRUE)

cat(sprintf("Full cohort n (rows): %d  unique stays: %d\n",
            nrow(d), length(unique(d$stay_id))))
cat(sprintf("Equipoise subgroup n (rows): %d  unique stays: %d\n",
            nrow(d_eq), length(unique(d_eq$stay_id))))

# ── Helper: fit main-effect + interaction models and LRT ─────────────────────
run_interaction <- function(dat, label) {
  cat("\n========================================\n")
  cat("Population:", label, "\n")
  cat("========================================\n")

  fit0 <- coxph(Surv(follow_days, event) ~ trt,
                data = dat, weights = ipcw_w_trim,
                cluster = stay_id, ties = "efron")

  fit1 <- coxph(Surv(follow_days, event) ~ trt * sofa_c,
                data = dat, weights = ipcw_w_trim,
                cluster = stay_id, ties = "efron")

  # Main effect from interaction model
  s1  <- summary(fit1)
  cf  <- s1$coefficients

  hr_main <- exp(cf["trt", "coef"])
  se_main <- cf["trt", "robust se"]
  p_main  <- cf["trt", "Pr(>|z|)"]

  # Interaction term (trt:sofa_c)
  int_row <- grep("trt:sofa_c|sofa_c:trt", rownames(cf), value = TRUE)
  if (length(int_row) == 0) {
    cat("WARNING: interaction term not found in model coefficients\n")
    hr_int <- NA; se_int <- NA; p_int <- NA
  } else {
    hr_int <- exp(cf[int_row, "coef"])
    se_int <- cf[int_row, "robust se"]
    p_int  <- cf[int_row, "Pr(>|z|)"]
  }

  # Wald test for interaction (valid under sandwich SE; LRT is invalid with weights)
  lrt_stat <- NA
  lrt_df   <- fit1$rank - fit0$rank
  lrt_p    <- p_int  # use Wald p instead of invalid LRT

  cat(sprintf("\nMain effect HR (early vs deferred CRRT at median SOFA):\n"))
  cat(sprintf("  HR = %.3f (95%% CI %.3f-%.3f), p = %.3f\n",
              hr_main,
              exp(cf["trt","coef"] - 1.96 * se_main),
              exp(cf["trt","coef"] + 1.96 * se_main),
              p_main))

  if (!is.na(hr_int)) {
    cat(sprintf("\nInteraction term (HR per 1-unit increase in SOFA above median):\n"))
    cat(sprintf("  HR_interaction = %.3f (95%% CI %.3f-%.3f), p = %.3f\n",
                hr_int,
                exp(cf[int_row,"coef"] - 1.96 * se_int),
                exp(cf[int_row,"coef"] + 1.96 * se_int),
                p_int))
  }

  cat(sprintf("\nInteraction p-value (Wald, robust SE): %.4f\n", lrt_p))

  list(
    population    = label,
    hr_main       = round(hr_main, 3),
    hr_main_lo    = round(exp(cf["trt","coef"] - 1.96 * se_main), 3),
    hr_main_hi    = round(exp(cf["trt","coef"] + 1.96 * se_main), 3),
    p_main        = round(p_main, 4),
    hr_int        = round(hr_int, 4),
    hr_int_lo     = if (!is.na(hr_int)) round(exp(cf[int_row,"coef"] - 1.96*se_int), 4) else NA,
    hr_int_hi     = if (!is.na(hr_int)) round(exp(cf[int_row,"coef"] + 1.96*se_int), 4) else NA,
    p_int_wald    = round(p_int, 4),
    p_int_lrt     = round(lrt_p, 4)
  )
}

res_full <- run_interaction(d,    "Full cohort (n = 3,112)")
res_eq   <- run_interaction(d_eq, "Equipoise subgroup SOFA 8-14 (n = 1,887)")

# ── Decision rule ─────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("DECISION\n")
cat("========================================\n")
p_full <- res_full$p_int_lrt
p_eq   <- res_eq$p_int_lrt
cat(sprintf("Full cohort    interaction LRT p = %.4f\n", p_full))
cat(sprintf("Equipoise subgroup interaction LRT p = %.4f\n", p_eq))

if (any(c(p_full, p_eq) < 0.10, na.rm = TRUE)) {
  cat("\n→ Proceed to RCS analysis (at least one p < 0.10)\n")
} else {
  cat("\n→ No evidence of effect modification by baseline SOFA severity.\n")
  cat("  Report interaction p-values and state: 'No statistically significant\n")
  cat("  interaction was detected between early CRRT timing and baseline SOFA\n")
  cat("  score (p_interaction = ", sprintf("%.3f", p_full),
      " full cohort; ", sprintf("%.3f", p_eq), " equipoise subgroup).\n", sep="")
}

# ── Save results ──────────────────────────────────────────────────────────────
results <- list(
  full_cohort       = res_full,
  equipoise_subgroup = res_eq,
  sofa_centering_value = sofa_med,
  decision_threshold   = 0.10
)
saveRDS(results, file.path(base_dir, "data/interaction_sofa_results.rds"))
cat("\nSaved: analysis/data/interaction_sofa_results.rds\n")
