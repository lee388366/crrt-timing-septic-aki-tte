# Weight-model sensitivity for the recommended Def C eICU rescue cohort.
# Reads outputs from 04_eicu_diagnose_rescue.R; writes diagnosis-only sensitivity tables.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(survival)
})

out_dir <- "CONFIGURE_OTHER_SOURCE_ROOT/diagnose_rescue"
data_dir <- file.path(out_dir, "data")
report_dir <- file.path(out_dir, "reports")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

cohort_path <- file.path(data_dir, "task4_recommended_def_c_cohort_base.csv")
if (!file.exists(cohort_path)) {
  stop("Missing recommended cohort. Run 04_eicu_diagnose_rescue.R first.")
}

cohort <- read_csv(cohort_path, show_col_types = FALSE)
TRIM_Q <- as.numeric(Sys.getenv("EICU_WEIGHT_TRIM_Q", "0.99"))
MIMIC_HR_REFERENCE <- as.numeric(Sys.getenv("MIMIC_STRICTK3_V2_HR", "1.044"))

numeric_covariates <- c(
  "age", "apachescore", "acutephysiologyscore", "creatinine_baseline",
  "creatinine_at_kdigo3", "lactate_near_t0", "weight_kg",
  "apache_meanbp", "apache_ph", "apache_bun", "apache_creatinine",
  "apache_urine"
)
available_numeric_covariates <- intersect(numeric_covariates, names(cohort))

impute_median <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    df[[v]][df[[v]] %in% c(-1, -99, -999)] <- NA
    med <- median(df[[v]], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    df[[v]][is.na(df[[v]])] <- med
  }
  df
}

standardize_numeric <- function(df, vars) {
  for (v in vars) {
    if (!v %in% names(df)) next
    s <- sd(df[[v]], na.rm = TRUE)
    m <- mean(df[[v]], na.rm = TRUE)
    if (is.finite(s) && s > 0) df[[v]] <- (df[[v]] - m) / s
  }
  df
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- weighted.mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

smd_numeric <- function(df, var, weight_var = NULL) {
  a <- df[df$early_crrt %in% TRUE & is.finite(df[[var]]), ]
  b <- df[df$early_crrt %in% FALSE & is.finite(df[[var]]), ]
  if (nrow(a) < 2 || nrow(b) < 2) return(NA_real_)
  if (is.null(weight_var)) {
    ma <- mean(a[[var]], na.rm = TRUE)
    mb <- mean(b[[var]], na.rm = TRUE)
    va <- stats::var(a[[var]], na.rm = TRUE)
    vb <- stats::var(b[[var]], na.rm = TRUE)
  } else {
    ma <- weighted.mean(a[[var]], a[[weight_var]], na.rm = TRUE)
    mb <- weighted.mean(b[[var]], b[[weight_var]], na.rm = TRUE)
    va <- weighted_var(a[[var]], a[[weight_var]])
    vb <- weighted_var(b[[var]], b[[weight_var]])
  }
  denom <- sqrt((va + vb) / 2)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  abs(ma - mb) / denom
}

make_balance_frame <- function(df) {
  out <- df
  out$bal_stage3_creatinine_only <- as.numeric(out$stage3_path == "creatinine_only")
  out$bal_stage3_urine_output_only <- as.numeric(out$stage3_path == "urine_output_only")
  out$bal_stage3_both <- as.numeric(out$stage3_path == "both_creatinine_and_uo")
  out$bal_gender_male <- as.numeric(tolower(as.character(out$gender)) == "male")
  out
}

balance_covariates <- function(df) {
  c(intersect(available_numeric_covariates, names(df)), grep("^bal_", names(df), value = TRUE))
}

safe_formula <- function(terms) {
  terms <- terms[terms %in% names(cohort) | terms %in% c("stage3_path")]
  if (length(terms) == 0) return(early_crrt ~ 1)
  as.formula(paste("early_crrt ~", paste(terms, collapse = " + ")))
}

model_specs <- tibble(
  model_name = c("intercept_only", "minimal_core", "severity_stage3", "full_current"),
  denom_terms = list(
    character(),
    c("age", "acutephysiologyscore", "creatinine_baseline", "lactate_near_t0", "stage3_path"),
    c("age", "apachescore", "acutephysiologyscore", "stage3_path"),
    c("age", "apachescore", "acutephysiologyscore", "creatinine_baseline",
      "creatinine_at_kdigo3", "lactate_near_t0", "weight_kg", "apache_meanbp",
      "apache_ph", "apache_bun", "apache_creatinine", "apache_urine", "stage3_path")
  ),
  num_terms = list(
    character(),
    character(),
    c("age"),
    c("age", "apachescore", "acutephysiologyscore")
  )
)

run_model <- function(df, subset_label, model_name, denom_terms, num_terms) {
  landmark <- df %>%
    filter(follow_time_28d > 1) %>%
    mutate(
      early_crrt = as.logical(early_crrt),
      event_after_landmark = died_28d %in% TRUE & !is.na(days_to_death) & days_to_death > 1,
      follow_after_landmark = pmax(0.001, pmin(27, follow_time_28d - 1)),
      stage3_path = factor(stage3_path),
      gender = factor(gender)
    ) %>%
    impute_median(available_numeric_covariates)

  if (nrow(landmark) < 30 || length(unique(landmark$early_crrt)) < 2) {
    return(tibble(
      subset = subset_label,
      model_name = model_name,
      status = "insufficient_contrast",
      n = nrow(landmark),
      early_n = sum(landmark$early_crrt, na.rm = TRUE)
    ))
  }

  denom_formula <- safe_formula(denom_terms)
  num_formula <- safe_formula(num_terms)
  model_df <- standardize_numeric(landmark, intersect(available_numeric_covariates, unique(c(denom_terms, num_terms))))
  denom_fit <- tryCatch(suppressWarnings(glm(denom_formula, data = model_df, family = binomial())), error = function(e) e)
  num_fit <- tryCatch(suppressWarnings(glm(num_formula, data = model_df, family = binomial())), error = function(e) e)
  status <- "ok"
  if (inherits(denom_fit, "error") || inherits(num_fit, "error")) {
    status <- "fallback_intercept_only"
    denom_fit <- glm(early_crrt ~ 1, data = model_df, family = binomial())
    num_fit <- glm(early_crrt ~ 1, data = model_df, family = binomial())
    denom_formula <- early_crrt ~ 1
    num_formula <- early_crrt ~ 1
  }

  p_denom <- pmin(pmax(as.numeric(predict(denom_fit, newdata = model_df, type = "response")), 0.01), 0.99)
  p_num <- pmin(pmax(as.numeric(predict(num_fit, newdata = model_df, type = "response")), 0.01), 0.99)
  analysis_df <- landmark %>%
    mutate(
      p_denom = p_denom,
      p_num = p_num,
      sw = if_else(early_crrt, p_num / p_denom, (1 - p_num) / (1 - p_denom)),
      sw_trim = pmin(sw, as.numeric(quantile(sw, TRIM_Q, na.rm = TRUE)))
    )
  cox_fit <- tryCatch(
    coxph(Surv(follow_after_landmark, event_after_landmark) ~ early_crrt,
          data = analysis_df, weights = sw_trim, robust = TRUE),
    error = function(e) e
  )
  if (inherits(cox_fit, "error")) {
    hr <- ci_low <- ci_high <- p_value <- NA_real_
  } else {
    beta <- unname(coef(cox_fit)[1])
    se <- sqrt(vcov(cox_fit)[1, 1])
    hr <- exp(beta)
    ci_low <- exp(beta - 1.96 * se)
    ci_high <- exp(beta + 1.96 * se)
    p_value <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  }
  balance_df <- make_balance_frame(analysis_df)
  smd_tbl <- tibble(covariate = balance_covariates(balance_df)) %>%
    mutate(
      smd_unweighted = vapply(covariate, function(v) smd_numeric(balance_df, v), numeric(1)),
      smd_weighted = vapply(covariate, function(v) smd_numeric(balance_df, v, "sw_trim"), numeric(1))
    )

  tibble(
    subset = subset_label,
    model_name = model_name,
    status = status,
    n = nrow(analysis_df),
    early_n = sum(analysis_df$early_crrt, na.rm = TRUE),
    ess = sum(analysis_df$sw_trim, na.rm = TRUE)^2 / sum(analysis_df$sw_trim^2, na.rm = TRUE),
    ess_early = with(analysis_df[analysis_df$early_crrt %in% TRUE, ],
                     sum(sw_trim, na.rm = TRUE)^2 / sum(sw_trim^2, na.rm = TRUE)),
    max_weighted_smd = max(smd_tbl$smd_weighted, na.rm = TRUE),
    worst_covariate = smd_tbl$covariate[which.max(smd_tbl$smd_weighted)],
    hr = hr,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_ratio_high_to_low = ci_high / ci_low,
    hr_ratio_vs_mimic = hr / MIMIC_HR_REFERENCE,
    usable_gate = early_n >= 30 &
      (sum(analysis_df$sw_trim[analysis_df$early_crrt], na.rm = TRUE)^2 / sum(analysis_df$sw_trim[analysis_df$early_crrt]^2, na.rm = TRUE)) >= 30 &
      max(smd_tbl$smd_weighted, na.rm = TRUE) <= 0.10 &
      is.finite(ci_high / ci_low) & (ci_high / ci_low) <= 4 &
      (hr / MIMIC_HR_REFERENCE) >= 0.80 & (hr / MIMIC_HR_REFERENCE) <= 1.25,
    denominator_formula = paste(deparse(denom_formula), collapse = " "),
    numerator_formula = paste(deparse(num_formula), collapse = " ")
  )
}

subsets <- list(
  high_fidelity = cohort %>% filter(high_fidelity_hospital %in% TRUE),
  all_hospitals = cohort
)

results <- bind_rows(lapply(names(subsets), function(snm) {
  bind_rows(lapply(seq_len(nrow(model_specs)), function(i) {
    run_model(subsets[[snm]], snm, model_specs$model_name[i], model_specs$denom_terms[[i]], model_specs$num_terms[[i]])
  }))
}))

write_csv(results, file.path(data_dir, "task4_weight_model_sensitivity.csv"))

fmt_table <- function(df, digits = 3) {
  x <- df
  for (nm in names(x)) {
    if (is.numeric(x[[nm]])) x[[nm]] <- format(round(x[[nm]], digits), trim = TRUE, scientific = FALSE)
    else x[[nm]] <- as.character(x[[nm]])
    x[[nm]][is.na(x[[nm]])] <- ""
  }
  c(
    paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  )
}

report <- c(
  "# eICU Weight-Model Sensitivity",
  "",
  "This sensitivity asks whether the Def C eICU result can be rescued by changing the stabilized-weight model.",
  "",
  fmt_table(results %>% select(subset, model_name, n, early_n, ess_early, max_weighted_smd, worst_covariate, hr, ci_low, ci_high, ci_ratio_high_to_low, hr_ratio_vs_mimic, usable_gate), digits = 3),
  "",
  "Conclusion: if all usable gates remain false, the failure is not a single-formula artifact; it reflects sparse early treatment support and persistent imbalance."
)
writeLines(report, file.path(report_dir, "eicu_weight_model_sensitivity_report.md"))

message("Saved weight-model sensitivity outputs under: ", out_dir)
