# =============================================================================
# Step 3+4: Clone-Censor-Weighting (CCW) + IPCW with Time-Varying Covariates
#
# Per methods:
#   Strategy A denominator: cross-sectional at hour 23, covariates L₀ + L₂₃
#   Strategy B denominator: pooled person-hours 0–23, covariates L₀ + L_k + hour_index
#   Numerator (both):       L₀ only (+ hour_index for B)
#   Regularization:         Ridge (alpha=0), 10-fold CV, lambda.1se
#   Missing time-varying:   LOCF (already applied in 03a)
# =============================================================================

library(tidyverse); library(lubridate); library(glmnet)

baseline  <- readRDS("analysis/data/baseline.rds")
tv_hourly <- readRDS("analysis/data/tv_hourly.rds")
tv_hour23 <- readRDS("analysis/data/tv_hour23.rds")

# First ICU stay per patient, SOFA available
cohort <- baseline %>%
  group_by(subject_id) %>% slice_min(icu_intime, n=1) %>% ungroup() %>%
  filter(!is.na(sofa_t0))

cat("Analytic cohort:", nrow(cohort), "\n")

# Baseline covariate names (L₀)
L0_vars <- c("age", "sofa_t0", "apsiii", "creatinine_baseline",
             "nee_at_t0", "lactate_t0", "uo_rate_t0", "map_t0", "cci")

# Time-varying covariate names (L_k)
Lk_vars <- c("map_h", "nee_h", "lactate_h", "sofa_h",
             "kdigo_h", "creatinine_h", "creatinine_delta", "fluid_balance_h")

# Impute missing L₀ with cohort median
impute_median <- function(df, vars) {
  for (v in vars) {
    med <- median(df[[v]], na.rm=TRUE)
    df[[v]][is.na(df[[v]])] <- med
  }
  df
}
cohort <- impute_median(cohort, L0_vars)

# =============================================================================
# STEP 3: Clone-Censor-Weighting
# =============================================================================
GRACE <- 24

ccw <- cohort %>%
  crossing(strategy = c("A_early", "B_deferred")) %>%
  mutate(
    clone_id = paste(stay_id, strategy, sep="_"),
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

cat("\nCCW clones:", nrow(ccw),
    "| Censored A:", sum(ccw$protocol_censored & ccw$strategy=="A_early"),
    "| Censored B:", sum(ccw$protocol_censored & ccw$strategy=="B_deferred"), "\n")

# =============================================================================
# STEP 4A: IPCW — Strategy A (cross-sectional at hour 23)
#
# Per methods: "denominator model fitted cross-sectionally on at-risk observations
#               at hour 23—the time point immediately preceding the single
#               censoring event at hour 24"
# Outcome:     CRRT initiated within 24h of t₀ (early_crrt)
# Denominator: Ridge logistic ~ L₀ + L₂₃
# Numerator:   Plain logistic ~ L₀ only
# =============================================================================

tv_hour23_uniq <- tv_hour23 %>% distinct(stay_id, .keep_all=TRUE)

A_data <- ccw %>%
  filter(strategy == "A_early") %>%
  left_join(tv_hour23_uniq, by="stay_id") %>%
  mutate(y_A = as.integer(early_crrt))

# Impute time-varying covariates at hour 23 with cohort median
A_data <- impute_median(A_data, Lk_vars)

# --- Denominator: L₀ + L₂₃ ---
X_A_denom <- model.matrix(
  ~ age + sofa_t0 + apsiii + creatinine_baseline + nee_at_t0 +
    lactate_t0 + uo_rate_t0 + map_t0 + cci +
    map_h + nee_h + lactate_h + sofa_h + kdigo_h +
    creatinine_h + creatinine_delta + fluid_balance_h,
  data = A_data)[, -1]

set.seed(42)
cv_A <- cv.glmnet(X_A_denom, A_data$y_A, family="binomial", alpha=0,
                  nfolds=10, type.measure="deviance")
p_denom_A <- as.vector(predict(cv_A, X_A_denom, s="lambda.1se", type="response"))
cat(sprintf("Strategy A: Ridge lambda=%.4f | mean P(CRRT|L0,L23)=%.3f\n",
            cv_A$lambda.1se, mean(p_denom_A)))

# --- Numerator: L₀ only ---
lm_num_A <- glm(y_A ~ age + sofa_t0 + apsiii + creatinine_baseline +
                  nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci,
                data=A_data, family=binomial())
p_num_A <- predict(lm_num_A, type="response")

# Stabilized weights
A_data <- A_data %>%
  mutate(
    p_denom = p_denom_A,
    p_num   = p_num_A,
    ipcw_w  = if_else(early_crrt,
                      p_num   / pmax(p_denom,   0.01),
                      (1-p_num) / pmax(1-p_denom, 0.01))
  )

cat("Strategy A weight summary:\n"); print(summary(A_data$ipcw_w))

# =============================================================================
# STEP 4B: IPCW — Strategy B (pooled longitudinal, hours 0–23)
#
# Per methods: "denominator model pooled all at-risk person-hours (hours 0–23)
#               with hour_index included as an additional predictor"
# Outcome:     CRRT initiated THIS hour (within-hour event)
# Denominator: Ridge logistic ~ L₀ + L_k + hour_index
# Numerator:   Plain logistic ~ L₀ + hour_index
# =============================================================================

B_base <- ccw %>% filter(strategy == "B_deferred")

# Build person-hour rows for the grace period
B_long <- B_base %>%
  select(stay_id, early_crrt, h_crrt, all_of(L0_vars)) %>%
  left_join(tv_hourly, by="stay_id", relationship="many-to-many") %>%
  mutate(
    crrt_this_hour = !is.na(h_crrt) & h_crrt >= hour_index &
                     h_crrt <  (hour_index + 1) & h_crrt <= GRACE,
    at_risk = is.na(h_crrt) | h_crrt > hour_index
  ) %>%
  filter(at_risk)

# Impute missing Lk with cohort median within B_long
B_long <- impute_median(B_long, Lk_vars)

cat("Strategy B long format:", nrow(B_long), "person-hours\n")

# --- Denominator: L₀ + L_k + hour_index ---
X_B_denom <- model.matrix(
  ~ age + sofa_t0 + apsiii + creatinine_baseline + nee_at_t0 +
    lactate_t0 + uo_rate_t0 + map_t0 + cci + hour_index +
    map_h + nee_h + lactate_h + sofa_h + kdigo_h +
    creatinine_h + creatinine_delta + fluid_balance_h,
  data = B_long)[, -1]

y_B <- as.integer(B_long$crrt_this_hour)

set.seed(42)
cv_B <- cv.glmnet(X_B_denom, y_B, family="binomial", alpha=0,
                  nfolds=10, type.measure="deviance")
p_denom_B <- as.vector(predict(cv_B, X_B_denom, s="lambda.1se", type="response"))
cat(sprintf("Strategy B: Ridge lambda=%.4f | mean P(CRRT@hour|L0,Lk)=%.4f\n",
            cv_B$lambda.1se, mean(p_denom_B)))

# --- Numerator: L₀ + hour_index ---
lm_num_B <- glm(crrt_this_hour ~ age + sofa_t0 + apsiii + creatinine_baseline +
                  nee_at_t0 + lactate_t0 + uo_rate_t0 + map_t0 + cci + hour_index,
                data=B_long, family=binomial())
p_num_B <- predict(lm_num_B, type="response")

B_long <- B_long %>%
  mutate(
    p_denom_h = p_denom_B,
    p_num_h   = p_num_B,
    w_contrib = (1 - p_num_h) / pmax(1 - p_denom_h, 0.01)
  )

# Cumulative product per patient
B_weights <- B_long %>%
  group_by(stay_id) %>%
  summarise(ipcw_w = prod(w_contrib, na.rm=TRUE), .groups="drop")

cat("Strategy B weight summary:\n"); print(summary(B_weights$ipcw_w))

# =============================================================================
# Trim extreme weights at 99th percentile (per methods)
# =============================================================================
w99_A <- quantile(A_data$ipcw_w, 0.99, na.rm=TRUE)
w99_B <- quantile(B_weights$ipcw_w, 0.99, na.rm=TRUE)

A_data    <- A_data    %>% mutate(ipcw_w_trim = pmin(ipcw_w, w99_A))
B_weights <- B_weights %>% mutate(ipcw_w_trim = pmin(ipcw_w, w99_B))
cat(sprintf("\nWeight trim caps: A=%.3f, B=%.3f\n", w99_A, w99_B))

# =============================================================================
# Assemble final CCW dataset with weights
# =============================================================================
A_weights <- A_data %>%
  distinct(stay_id, .keep_all=TRUE) %>%
  select(stay_id, ipcw_w, ipcw_w_trim)

ccw_A <- ccw %>%
  filter(strategy == "A_early") %>%
  left_join(A_weights, by="stay_id", relationship="many-to-one")

ccw_B <- ccw %>%
  filter(strategy == "B_deferred") %>%
  left_join(B_weights, by="stay_id", relationship="many-to-one") %>%
  mutate(ipcw_w      = coalesce(ipcw_w, 1),
         ipcw_w_trim = coalesce(ipcw_w_trim, 1))

ccw_final <- bind_rows(ccw_A, ccw_B) %>%
  mutate(
    strategy_label = if_else(strategy=="A_early", "Early CRRT", "Deferred"),
    equipoise      = !is.na(sofa_t0) & sofa_t0 >= 8 & sofa_t0 <= 14
  )

cat("\nFinal CCW dataset:", nrow(ccw_final), "rows |",
    nrow(ccw_final)/2, "patients\n")
cat("Mean weight A (trimmed):", round(mean(ccw_final$ipcw_w_trim[ccw_final$strategy=="A_early"], na.rm=T), 3), "\n")
cat("Mean weight B (trimmed):", round(mean(ccw_final$ipcw_w_trim[ccw_final$strategy=="B_deferred"], na.rm=T), 3), "\n")

# =============================================================================
# Weight diagnostics: weighted SMD for key covariates (per methods)
# =============================================================================
smd_simple <- function(ccw_df, var) {
  grpA <- ccw_df %>% filter(strategy=="A_early",    !is.na(.data[[var]]))
  grpB <- ccw_df %>% filter(strategy=="B_deferred", !is.na(.data[[var]]))
  mA <- weighted.mean(grpA[[var]], grpA$ipcw_w_trim)
  mB <- weighted.mean(grpB[[var]], grpB$ipcw_w_trim)
  sdA <- sd(grpA[[var]])
  sdB <- sd(grpB[[var]])
  abs(mA - mB) / sqrt((sdA^2 + sdB^2) / 2)
}

cat("\n=== Weighted SMD (target <0.10) ===\n")
for (v in c("age","sofa_t0","apsiii","creatinine_baseline","nee_at_t0","lactate_t0")) {
  if (v %in% names(ccw_final)) {
    s <- tryCatch({
      grpA <- ccw_final %>% filter(strategy=="A_early",    !is.na(.data[[v]]))
      grpB <- ccw_final %>% filter(strategy=="B_deferred", !is.na(.data[[v]]))
      mA <- weighted.mean(grpA[[v]], grpA$ipcw_w_trim)
      mB <- weighted.mean(grpB[[v]], grpB$ipcw_w_trim)
      sdA <- sd(grpA[[v]]); sdB <- sd(grpB[[v]])
      abs(mA - mB) / sqrt((sdA^2 + sdB^2) / 2)
    }, error=function(e) NA_real_)
    cat(sprintf("  %-25s SMD = %.3f %s\n", v, s,
                if(!is.na(s) && s > 0.10) "⚠" else ""))
  }
}

# =============================================================================
# Save
# =============================================================================
saveRDS(ccw_final, "analysis/data/ccw_final.rds")
cat("\n✅ Saved: analysis/data/ccw_final.rds\n")
