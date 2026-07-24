# =============================================================================
# Step 5: IPCW-Weighted Cox Model + RMST + Kaplan-Meier Curves
# =============================================================================

library(tidyverse); library(survival); library(sandwich); library(lmtest)
library(boot); library(ggplot2)

ccw <- readRDS("analysis/data/ccw_final.rds")

# Remove any duplicated clone_ids from many-to-many join issue
ccw <- ccw %>% distinct(clone_id, .keep_all = TRUE)

cat("CCW dataset:", nrow(ccw), "rows,", nrow(ccw)/2, "patients\n")

run_analysis <- function(data, label) {
  cat("\n", strrep("=", 60), "\n")
  cat("ANALYSIS:", label, "(n =", nrow(data)/2, "patients)\n")
  cat(strrep("=", 60), "\n")

  # -------------------------------------------------------------------
  # 1. IPCW-Weighted Cox proportional hazards model
  #    Cluster on stay_id to account for paired clone structure
  # -------------------------------------------------------------------
  data <- data %>% mutate(
    trt = as.integer(strategy == "A_early"),  # 1=early, 0=deferred
    w   = ipcw_w_trim
  )

  cox_fit <- coxph(
    Surv(follow_days, event) ~ trt,
    data    = data,
    weights = w,
    cluster = stay_id,   # cluster-robust SE for paired clones
    robust  = TRUE
  )

  cox_s <- summary(cox_fit)
  hr    <- cox_s$coefficients["trt", "exp(coef)"]
  hr_lo <- cox_s$conf.int["trt", "lower .95"]
  hr_hi <- cox_s$conf.int["trt", "upper .95"]
  pval  <- cox_s$coefficients["trt", "Pr(>|z|)"]

  cat(sprintf("Cox HR (Early vs Deferred): %.3f (95%% CI: %.3f–%.3f), p=%.4f\n",
              hr, hr_lo, hr_hi, pval))

  # Proportional hazards test
  ph_test <- cox.zph(cox_fit)
  cat(sprintf("PH test (global p): %.4f\n", ph_test$table["GLOBAL", "p"]))

  # -------------------------------------------------------------------
  # 2. Weighted Kaplan-Meier curves
  # -------------------------------------------------------------------
  km_A <- survfit(Surv(follow_days, event) ~ 1,
                  data    = data %>% filter(strategy == "A_early"),
                  weights = w)
  km_B <- survfit(Surv(follow_days, event) ~ 1,
                  data    = data %>% filter(strategy == "B_deferred"),
                  weights = w)

  # -------------------------------------------------------------------
  # 3. RMST at 28 days (trapezoidal rule from KM)
  # -------------------------------------------------------------------
  rmst_from_km <- function(km_fit, tau = 28) {
    t <- c(0, km_fit$time)
    s <- c(1, km_fit$surv)
    t_trunc <- pmin(t, tau)
    idx <- t <= tau
    t_use <- c(t[idx], tau)
    s_use <- c(s[idx], s[max(which(t <= tau))])
    sum(diff(t_use) * s_use[-length(s_use)])
  }

  rmst_A <- rmst_from_km(km_A)
  rmst_B <- rmst_from_km(km_B)
  rmst_diff <- rmst_A - rmst_B
  cat(sprintf("RMST (28d): Early=%.2f days, Deferred=%.2f days, Diff=%.2f days\n",
              rmst_A, rmst_B, rmst_diff))

  # Bootstrap 95% CI for RMST difference (nonparametric, resample on stay_id)
  set.seed(42)
  unique_stays <- unique(data$stay_id)
  n_stays      <- length(unique_stays)
  B            <- 200   # 200 replicates, percentile method

  # Pre-split data by strategy for speed
  dA <- data[data$strategy == "A_early", c("stay_id","follow_days","event","w")]
  dB <- data[data$strategy == "B_deferred", c("stay_id","follow_days","event","w")]

  boot_diffs <- replicate(B, {
    samp_df <- data.frame(stay_id = sample(unique_stays, n_stays, replace = TRUE))
    sA  <- merge(samp_df, dA, by = "stay_id", all.x = TRUE)
    sB  <- merge(samp_df, dB, by = "stay_id", all.x = TRUE)
    tryCatch({
      rmst_from_km(survfit(Surv(follow_days, event) ~ 1, data=sA, weights=sA$w)) -
      rmst_from_km(survfit(Surv(follow_days, event) ~ 1, data=sB, weights=sB$w))
    }, error = function(e) NA_real_)
  })

  ci_lo <- quantile(boot_diffs, 0.025, na.rm = TRUE)
  ci_hi <- quantile(boot_diffs, 0.975, na.rm = TRUE)
  cat(sprintf("RMST diff 95%% CI (bootstrap): %.2f to %.2f days\n", ci_lo, ci_hi))

  # -------------------------------------------------------------------
  # 4. KM plot
  # -------------------------------------------------------------------
  km_df <- bind_rows(
    tibble(time=c(0, km_A$time), surv=c(1, km_A$surv), strategy="Early CRRT"),
    tibble(time=c(0, km_B$time), surv=c(1, km_B$surv), strategy="Deferred")
  )

  p <- ggplot(km_df, aes(x=time, y=1-surv, color=strategy)) +
    geom_step(linewidth=1) +
    scale_color_manual(values=c("Early CRRT"="#E41A1C","Deferred"="#377EB8")) +
    scale_x_continuous(limits=c(0,28), breaks=seq(0,28,7)) +
    scale_y_continuous(limits=c(0,1), labels=scales::percent) +
    labs(
      title    = paste("28-day Cumulative Mortality:", label),
      subtitle = sprintf("HR=%.2f (95%% CI %.2f–%.2f), p=%.3f | RMST diff=%.2f days",
                         hr, hr_lo, hr_hi, pval, rmst_diff),
      x="Days from t₀", y="Cumulative Mortality",
      color="Strategy"
    ) +
    theme_classic(base_size=13) +
    theme(legend.position="bottom")

  fname <- paste0("analysis/figures/km_",
                  gsub("[^a-zA-Z0-9]","_", label), ".pdf")
  dir.create("analysis/figures", showWarnings=FALSE, recursive=TRUE)
  ggsave(fname, p, width=7, height=5)
  cat(sprintf("KM plot saved: %s\n", fname))

  list(hr=hr, hr_lo=hr_lo, hr_hi=hr_hi, pval=pval,
       rmst_A=rmst_A, rmst_B=rmst_B, rmst_diff=rmst_diff,
       ci_lo=ci_lo, ci_hi=ci_hi)
}

# =============================================================================
# Main analyses
# =============================================================================

# 1. Full cohort
res_full <- run_analysis(ccw, "Full cohort (n=3523)")

# 2. Primary: equipoise zone SOFA 8-14
res_eq <- run_analysis(
  filter(ccw, !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14),
  "Equipoise SOFA 8-14"
)

# 3. Sensitivity: conservative SOFA 9-13
res_cons <- run_analysis(
  filter(ccw, !is.na(sofa_t0) & sofa_t0 >= 9 & sofa_t0 <= 13),
  "Conservative SOFA 9-13"
)

# =============================================================================
# Summary table
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("RESULTS SUMMARY\n")
cat(strrep("=", 60), "\n")
results <- tibble(
  Analysis   = c("Full cohort", "Equipoise SOFA 8-14", "Conservative SOFA 9-13"),
  N_patients = c(nrow(filter(ccw, !duplicated(stay_id))),
                 nrow(filter(ccw, !duplicated(paste(stay_id, sofa_t0 >= 8 & sofa_t0 <= 14)),
                              sofa_t0 >= 8 & sofa_t0 <= 14)),
                 nrow(filter(ccw, !duplicated(paste(stay_id, sofa_t0 >= 9 & sofa_t0 <= 13)),
                              sofa_t0 >= 9 & sofa_t0 <= 13))),
  HR         = c(res_full$hr, res_eq$hr, res_cons$hr),
  HR_lo      = c(res_full$hr_lo, res_eq$hr_lo, res_cons$hr_lo),
  HR_hi      = c(res_full$hr_hi, res_eq$hr_hi, res_cons$hr_hi),
  p_value    = c(res_full$pval, res_eq$pval, res_cons$pval),
  RMST_diff  = c(res_full$rmst_diff, res_eq$rmst_diff, res_cons$rmst_diff),
  CI_lo      = c(res_full$ci_lo, res_eq$ci_lo, res_cons$ci_lo),
  CI_hi      = c(res_full$ci_hi, res_eq$ci_hi, res_cons$ci_hi)
) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

print(results, width=Inf)

saveRDS(results, "analysis/data/main_results.rds")
write_csv(results, "analysis/data/main_results.csv")
cat("\n✅ Results saved to analysis/data/main_results.csv\n")
