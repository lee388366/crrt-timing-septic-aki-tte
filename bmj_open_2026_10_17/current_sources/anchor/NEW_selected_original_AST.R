suppressPackageStartupMessages({library(dplyr);library(survival);library(readr)})
B <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
dod <- readRDS(Sys.getenv("CRRT_DOD_RDS")) %>% mutate(dod = as.Date(dod))
mk <- function(x){
  x %>% distinct(clone_id, .keep_all = TRUE) %>% left_join(dod, by = "subject_id") %>%
  mutate(t0 = as.POSIXct(t0, tz = "UTC"), deathtime = as.POSIXct(deathtime, tz = "UTC"),
    di = as.numeric(difftime(deathtime, t0, units = "days")), dd = as.numeric(dod - as.Date(t0)),
    dc = coalesce(di, dd),
    ev = !is.na(dc) & dc >= 0 & dc <= 28, dh = dc * 24,
    fo = pmin(28 * 24, coalesce(protocol_censor_h, 28 * 24), if_else(ev, dh, 28 * 24), na.rm = TRUE),
    ev = ev & (is.na(protocol_censor_h) | dh < protocol_censor_h), fo = fo / 24,
    sf = factor(strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))) %>% as.data.frame()
}

cox1 <- function(dat, lab, pop, w = "ipcw_w_trim", strata_var = NULL, tcol = "fo", ecol = "ev") {
  f <- if (is.null(strata_var)) {
    coxph(as.formula(sprintf("Surv(%s,%s) ~ sf", tcol, ecol)), data = dat, weights = dat[[w]],
          cluster = stay_id, robust = TRUE, ties = "efron")
  } else {
    coxph(as.formula(sprintf("Surv(%s,%s) ~ sf + strata(%s)", tcol, ecol, strata_var)), data = dat,
          weights = dat[[w]], cluster = stay_id, robust = TRUE, ties = "efron")
  }
  co <- coef(f)[["sfEarly CRRT"]]; se <- sqrt(diag(vcov(f)))[["sfEarly CRRT"]]
  php <- tryCatch(cox.zph(f)$table["GLOBAL", "p"], error = function(e) NA_real_)
  tibble(analysis = lab, population = pop, n_patients = length(unique(dat$stay_id)),
         HR = exp(co), ci_lower = exp(co - 1.96 * se), ci_upper = exp(co + 1.96 * se),
         p_value = 2 * pnorm(-abs(co / se)), ph_global_p = php)
}


rows <- list()
a_raw <- readRDS(Sys.getenv("CRRT_ANCHOR_MINUS6_CCW"))
a_raw$sf <- factor(a_raw$strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))
rows[["6h_old_eq"]]  <- cox1(as.data.frame(a_raw)[a_raw$equipoise %in% TRUE, ], "legacy K3-6h anchor (OLD dao)", "SOFA 8-14", tcol = "follow_days", ecol = "event")
rows[["6h_old_all"]] <- cox1(as.data.frame(a_raw), "legacy K3-6h anchor (OLD dao)", "Full cohort", tcol = "follow_days", ecol = "event")
a <- mk(a_raw)
rows[["6h_new_eq"]]  <- cox1(a[a$equipoise %in% TRUE, ], "legacy K3-6h anchor (NEW CH-01)", "SOFA 8-14")
rows[["6h_new_all"]] <- cox1(a, "legacy K3-6h anchor (NEW CH-01)", "Full cohort")

## ---- 2) KDIGO2 anchor ----
k_raw <- readRDS(Sys.getenv("CRRT_ANCHOR_STAGE2_CCW"))
k_raw$sf <- factor(k_raw$strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))
rows[["k2_old_eq"]]  <- cox1(as.data.frame(k_raw)[k_raw$equipoise %in% TRUE, ], "KDIGO2 anchor (OLD dao)", "SOFA 8-14", tcol = "follow_days", ecol = "event")
rows[["k2_old_all"]] <- cox1(as.data.frame(k_raw), "KDIGO2 anchor (OLD dao)", "Full cohort", tcol = "follow_days", ecol = "event")
k <- mk(k_raw)
rows[["k2_new_eq"]]  <- cox1(k[k$equipoise %in% TRUE, ], "KDIGO2 anchor (NEW CH-01)", "SOFA 8-14")
rows[["k2_new_all"]] <- cox1(k, "KDIGO2 anchor (NEW CH-01)", "Full cohort")


out <- bind_rows(rows)

dir.create(Sys.getenv("CRRT_ANCHOR_OUT"),recursive=TRUE,showWarnings=FALSE)
readr::write_csv(out,file.path(Sys.getenv("CRRT_ANCHOR_OUT"),"anchor_original_eight_rows.csv"))
