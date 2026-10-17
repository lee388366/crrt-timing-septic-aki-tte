run_cox <- function(ccw, label) {
  dat <- ccw %>%
    distinct(clone_id, .keep_all = TRUE) %>%
    mutate(trt = as.integer(strategy == "A_early")) %>%
    filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim))

  fit <- coxph(
    Surv(follow_days, event) ~ trt,
    data = dat,
    weights = ipcw_w_trim,
    cluster = stay_id,
    robust = TRUE,
    ties = "efron"
  )
  s <- summary(fit)
  tibble(
    analysis = label,
    n_patients = n_distinct(dat$stay_id),
    clone_rows = nrow(dat),
    early_n = n_distinct(dat$stay_id[dat$early_crrt %in% TRUE]),
    sofa8_14_n = n_distinct(dat$stay_id[dat$equipoise %in% TRUE]),
    sofa8_14_early_n = n_distinct(dat$stay_id[dat$equipoise %in% TRUE & dat$early_crrt %in% TRUE]),
    events = sum(dat$event, na.rm = TRUE),
    weighted_hr = unname(s$coefficients["trt", "exp(coef)"]),
    ci_lower = unname(s$conf.int["trt", "lower .95"]),
    ci_upper = unname(s$conf.int["trt", "upper .95"]),
    p_value = unname(s$coefficients["trt", "Pr(>|z|)"])
  )
}
