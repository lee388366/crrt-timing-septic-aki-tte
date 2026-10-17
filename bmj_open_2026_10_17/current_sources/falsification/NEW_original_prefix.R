# Strict clean first-KDIGO3 v2: negative-control falsification and QBA.
# Scope: SOFA 8-14 primary equipoise cohort only. Reuses existing ccw_final
# clone-censor-weighted data and IPCW weights. Does not modify earlier outputs.
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sandwich)
})

set.seed(20260716)

base_dir <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
data_dir <- Sys.getenv("CRRT_V2_DATA_DIR")
results_dir <- Sys.getenv("CRRT_FALSIFICATION_OUT")
figures_dir <- Sys.getenv("CRRT_FALSIFICATION_OUT")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

ccw <- readRDS(file.path(data_dir, "ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(strategy == "A_early"),
    strategy_name = if_else(trt == 1L, "A_early", "B_deferred")
  )

eq_ccw <- ccw %>%
  filter(equipoise %in% TRUE) %>%
  mutate(
    nc_high_sofa12 = as.integer(sofa_t0 >= 12),
    nc_high_creat3 = as.integer(creatinine_t0 >= 3)
  )

eq_patient <- eq_ccw %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  mutate(
    trt = as.integer(early_crrt),
    strategy_name = if_else(trt == 1L, "observed_early", "observed_deferred")
  )

candidate_definitions <- tibble(
  candidate = c("nc_high_sofa12", "nc_high_creat3"),
  label = c("Pre-exposure high organ-failure burden", "Pre-exposure creatinine-dominant AKI severity"),
  operational_definition = c(
    "SOFA score at t0 at least 12, within the prespecified SOFA 8-14 equipoise cohort.",
    "Serum creatinine at t0 at least 3 mg/dL, before assignment to the early or deferred CRRT strategy."
  ),
  causal_path_argument = c(
    "CRRT timing after t0 cannot cause the SOFA value already observed at t0. The variable can still share confounding-by-indication because higher baseline organ failure may influence clinicians toward earlier CRRT and is prognostic for death.",
    "CRRT timing after t0 cannot cause the creatinine value already observed at t0. The variable can still share the renal-severity indication structure that drives early CRRT and mortality risk."
  ),
  limitation = c(
    "This is a pre-exposure falsification endpoint, not a post-treatment clinical negative-control outcome.",
    "This is a pre-exposure falsification endpoint, not a post-treatment clinical negative-control outcome."
  )
)

safe_fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

robust_poisson_rr <- function(dat, outcome, weight_col = NULL, analysis_label) {
  dat <- dat %>%
    filter(!is.na(.data[[outcome]]), !is.na(trt)) %>%
    mutate(
      y = .data[[outcome]],
      w = if (is.null(weight_col)) 1 else .data[[weight_col]]
    ) %>%
    filter(!is.na(w), is.finite(w), w > 0)

  if (nrow(dat) == 0 || length(unique(dat$trt)) < 2 || length(unique(dat$y)) < 2) {
    return(tibble(
      analysis = analysis_label,
      candidate = outcome,
      n_rows = nrow(dat),
      n_patients = n_distinct(dat$stay_id),
      events_early = NA_integer_,
      events_deferred = NA_integer_,
      risk_early = NA_real_,
      risk_deferred = NA_real_,
      risk_difference = NA_real_,
      rr = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_
    ))
  }

  fit <- glm(y ~ trt, data = dat, family = poisson(link = "log"), weights = w)
  vc <- vcovCL(fit, cluster = dat$stay_id, type = "HC0")
  beta <- unname(coef(fit)["trt"])
  se <- sqrt(unname(vc["trt", "trt"]))
  if (!is.finite(se) || se < 1e-12) {
    se <- 0
  }
  z <- if (se == 0) NA_real_ else beta / se

  risk_summary <- dat %>%
    group_by(trt) %>%
    summarise(
      events = sum(y == 1, na.rm = TRUE),
      risk = weighted.mean(y, w, na.rm = TRUE),
      .groups = "drop"
    )
  risk_early <- risk_summary$risk[risk_summary$trt == 1]
  risk_deferred <- risk_summary$risk[risk_summary$trt == 0]
  events_early <- risk_summary$events[risk_summary$trt == 1]
  events_deferred <- risk_summary$events[risk_summary$trt == 0]

  tibble(
    analysis = analysis_label,
    candidate = outcome,
    n_rows = nrow(dat),
    n_patients = n_distinct(dat$stay_id),
    events_early = events_early,
    events_deferred = events_deferred,
    risk_early = risk_early,
    risk_deferred = risk_deferred,
    risk_difference = risk_early - risk_deferred,
    rr = exp(beta),
    ci_lower = exp(beta - 1.96 * se),
    ci_upper = exp(beta + 1.96 * se),
    p_value = if (is.na(z)) NA_real_ else 2 * pnorm(abs(z), lower.tail = FALSE)
  )
}

nc_results <- bind_rows(
  robust_poisson_rr(eq_patient, "nc_high_sofa12", NULL, "Observed actual early versus deferred, unweighted"),
  robust_poisson_rr(eq_patient, "nc_high_creat3", NULL, "Observed actual early versus deferred, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_sofa12", NULL, "Clone assigned strategies, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_creat3", NULL, "Clone assigned strategies, unweighted"),
  robust_poisson_rr(eq_ccw, "nc_high_sofa12", "ipcw_w_trim", "Clone assigned strategies, IPCW weighted"),
  robust_poisson_rr(eq_ccw, "nc_high_creat3", "ipcw_w_trim", "Clone assigned strategies, IPCW weighted")
) %>%
  left_join(candidate_definitions %>% select(candidate, label), by = "candidate") %>%
  relocate(label, .after = candidate)

write_csv(candidate_definitions, file.path(results_dir, "strictK3_v2_negctrl_candidate_definitions.csv"))
write_csv(nc_results, file.path(results_dir, "strictK3_v2_negctrl_results.csv"))

