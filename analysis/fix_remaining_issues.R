# =============================================================================
# fix_remaining_issues.R
# Three fixes:
#   G1 — Schoenfeld residuals test (PH assumption) → report results
#   G2 — IPCW-weighted logistic regression for renal recovery (pre-specified)
#   G3 — ICU LOS: update methods document (handled in Python script)
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(sandwich)
  library(lmtest)
  library(dplyr)
})

base_dir <- "analysis"
ccw_path <- file.path(base_dir, "data/ccw_final.rds")
master_path <- file.path(base_dir, "data/analysis_master.csv")

d      <- readRDS(ccw_path)
master <- read.csv(master_path)

d$strategy <- factor(d$strategy,
                     levels = c("B_deferred", "A_early"),
                     labels = c("Deferred CRRT", "Early CRRT"))

# ── G1: Schoenfeld residuals for PH assumption ──────────────────────────────
cat("\n====== G1: Schoenfeld Residuals Test ======\n")

# Primary equipoise population
d_eq <- d[d$equipoise == TRUE, ]
cat(sprintf("Equipoise (SOFA 8-14): %d rows, %d unique stays\n",
            nrow(d_eq), length(unique(d_eq$stay_id))))

fit_adj_eq <- coxph(
  Surv(follow_days, event) ~ strategy,
  data    = d_eq,
  weights = ipcw_w_trim,
  cluster = stay_id,
  ties    = "efron"
)

schoen_eq <- cox.zph(fit_adj_eq)
cat("\nEquipoise — Schoenfeld residuals (scaled):\n")
print(schoen_eq)

# Full cohort
fit_adj_full <- coxph(
  Surv(follow_days, event) ~ strategy,
  data    = d,
  weights = ipcw_w_trim,
  cluster = stay_id,
  ties    = "efron"
)

schoen_full <- cox.zph(fit_adj_full)
cat("\nFull cohort — Schoenfeld residuals:\n")
print(schoen_full)

# Store p-values for docx update
ph_p_eq   <- schoen_eq$table["GLOBAL", "p"]
ph_p_full <- schoen_full$table["GLOBAL", "p"]
cat(sprintf("\nPH p-value: Equipoise=%.3f, Full=%.3f\n", ph_p_eq, ph_p_full))

# ── G2: IPCW-weighted logistic regression for renal recovery ────────────────
cat("\n====== G2: IPCW-Weighted Logistic (Renal Recovery) ======\n")

# Merge renal_recovery_28d from analysis_master into CCW dataset
# Use one row per stay (patient-level, since recovery is a single binary outcome)
master_sub <- master[, c("stay_id", "renal_recovery_28d", "rrt_free_days_28d",
                          "icu_los_from_t0", "sofa_baseline")]

# Join to CCW (keep one row per clone, but recovery is at patient level)
# For IPCW logistic, we need patient-level data with IPCW weights from CCW
# Strategy: take one clone per patient (e.g., Strategy A clone) and attach recovery
# Actually for logistic, we want to use the full CCW framework analogously to Cox:
# Each clone contributes to the outcome with IPCW weights.
# The outcome for each clone is the same (recovery_bin), independent of strategy.
# Logistic: P(recovery | strategy, IPCW weight)

d_merge <- d %>%
  left_join(master_sub, by = "stay_id") %>%
  mutate(
    recovery_bin = as.integer(renal_recovery_28d == 1),
    strategy_bin = as.integer(strategy == "Early CRRT")
  )

# Full cohort IPCW logistic
fit_logit_full <- glm(
  recovery_bin ~ strategy_bin,
  data   = d_merge,
  family = binomial(),
  weights = ipcw_w_trim
)
# HC3 sandwich SEs clustered on stay_id
vcov_full  <- vcovCL(fit_logit_full, cluster = ~stay_id, type = "HC3")
coef_full  <- coeftest(fit_logit_full, vcov = vcov_full)
or_full    <- exp(coef_full["strategy_bin", "Estimate"])
ci_lo_full <- exp(coef_full["strategy_bin", "Estimate"] - 1.96 * coef_full["strategy_bin", "Std. Error"])
ci_hi_full <- exp(coef_full["strategy_bin", "Estimate"] + 1.96 * coef_full["strategy_bin", "Std. Error"])
pv_full    <- coef_full["strategy_bin", "Pr(>|z|)"]
cat(sprintf("\nIPCW logistic — Full cohort:\n  OR = %.3f [%.3f, %.3f]; p = %.3f\n",
            or_full, ci_lo_full, ci_hi_full, pv_full))

# Equipoise IPCW logistic
d_eq_merge <- d_merge[d_merge$equipoise == TRUE, ]
fit_logit_eq <- glm(
  recovery_bin ~ strategy_bin,
  data   = d_eq_merge,
  family = binomial(),
  weights = ipcw_w_trim
)
vcov_eq  <- vcovCL(fit_logit_eq, cluster = ~stay_id, type = "HC3")
coef_eq  <- coeftest(fit_logit_eq, vcov = vcov_eq)
or_eq    <- exp(coef_eq["strategy_bin", "Estimate"])
ci_lo_eq <- exp(coef_eq["strategy_bin", "Estimate"] - 1.96 * coef_eq["strategy_bin", "Std. Error"])
ci_hi_eq <- exp(coef_eq["strategy_bin", "Estimate"] + 1.96 * coef_eq["strategy_bin", "Std. Error"])
pv_eq    <- coef_eq["strategy_bin", "Pr(>|z|)"]
cat(sprintf("IPCW logistic — Equipoise (SOFA 8-14, n=1887):\n  OR = %.3f [%.3f, %.3f]; p = %.3f\n",
            or_eq, ci_lo_eq, ci_hi_eq, pv_eq))

# Save results to a simple text file for Python to read
results_txt <- sprintf(
  "SCHOENFELD_EQ_P=%.4f\nSCHOENFELD_FULL_P=%.4f\nIPCW_OR_FULL=%.3f\nIPCW_CI_LO_FULL=%.3f\nIPCW_CI_HI_FULL=%.3f\nIPCW_PV_FULL=%.4f\nIPCW_OR_EQ=%.3f\nIPCW_CI_LO_EQ=%.3f\nIPCW_CI_HI_EQ=%.3f\nIPCW_PV_EQ=%.4f\n",
  ph_p_eq, ph_p_full,
  or_full, ci_lo_full, ci_hi_full, pv_full,
  or_eq, ci_lo_eq, ci_hi_eq, pv_eq
)
writeLines(results_txt, file.path(base_dir, "data/fix_remaining_results.txt"))
cat("\nResults saved to data/fix_remaining_results.txt\n")
