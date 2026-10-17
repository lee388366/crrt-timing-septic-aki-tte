mk <- function(x) {
    x %>% distinct(clone_id, .keep_all = TRUE) %>% left_join(dod, 
        by = "subject_id") %>% mutate(t0 = as.POSIXct(t0, tz = "UTC"), 
        deathtime = as.POSIXct(deathtime, tz = "UTC"), di = as.numeric(difftime(deathtime, 
            t0, units = "days")), dd = as.numeric(dod - as.Date(t0)), 
        dc = coalesce(di, dd), ev = !is.na(dc) & dc >= 0 & dc <= 
            28, dh = dc * 24, fo = pmin(28 * 24, coalesce(protocol_censor_h, 
            28 * 24), if_else(ev, dh, 28 * 24), na.rm = TRUE), 
        ev = ev & (is.na(protocol_censor_h) | dh < protocol_censor_h), 
        fo = fo/24, sf = factor(strategy, levels = c("B_deferred", 
            "A_early"), labels = c("Deferred CRRT", "Early CRRT"))) %>% 
        as.data.frame()
}
cox1 <- function(dat, lab, pop, w = "ipcw_w_trim", strata_var = NULL, 
    tcol = "fo", ecol = "ev") {
    f <- if (is.null(strata_var)) {
        coxph(as.formula(sprintf("Surv(%s,%s) ~ sf", tcol, ecol)), 
            data = dat, weights = dat[[w]], cluster = stay_id, 
            robust = TRUE, ties = "efron")
    }
    else {
        coxph(as.formula(sprintf("Surv(%s,%s) ~ sf + strata(%s)", 
            tcol, ecol, strata_var)), data = dat, weights = dat[[w]], 
            cluster = stay_id, robust = TRUE, ties = "efron")
    }
    co <- coef(f)[["sfEarly CRRT"]]
    se <- sqrt(diag(vcov(f)))[["sfEarly CRRT"]]
    php <- tryCatch(cox.zph(f)$table["GLOBAL", "p"], error = function(e) NA_real_)
    tibble(analysis = lab, population = pop, n_patients = length(unique(dat$stay_id)), 
        HR = exp(co), ci_lower = exp(co - 1.96 * se), ci_upper = exp(co + 
            1.96 * se), p_value = 2 * pnorm(-abs(co/se)), ph_global_p = php)
}
rows[["l0_new_eq"]] <- cox1(l0[l0$equipoise %in% TRUE, ], "L0-only denominator (NEW CH-01)", 
    "SOFA 8-14")
rows[["l0_new_all"]] <- cox1(l0, "L0-only denominator (NEW CH-01)", 
    "Full cohort")
rows[["era_strat_eq"]] <- cox1(as.data.frame(dc)[dc$equipoise %in% 
    TRUE, ], "Calendar-era stratified Cox (NEW CH-01)", "SOFA 8-14", 
    strata_var = "anchor_year_group", tcol = "corrected_follow_days", 
    ecol = "corrected_event")
rows[["era_strat_all"]] <- cox1(as.data.frame(dc), "Calendar-era stratified Cox (NEW CH-01)", 
    "Full cohort", strata_var = "anchor_year_group", tcol = "corrected_follow_days", 
    ecol = "corrected_event")
