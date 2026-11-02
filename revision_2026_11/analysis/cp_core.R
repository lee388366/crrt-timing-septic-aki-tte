# Core functions for the corrected analysis: counting-process clone-censor-weighting with time-varying IPCW
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lubridate); library(survival); library(splines)})
L0 <- c("age", "sofa_t0_strict", "apsiii", "creatinine_baseline", "nee_at_t0", "lactate_t0", "uo_rate_t0", "map_t0", "cci")
Lk <- c("map_h", "nee_h", "lactate_h", "sofa_h", "kdigo_h", "creatinine_h", "creatinine_delta", "fluid_balance_h")
med_imp <- function(df, vars) { for (v in vars) df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE); df }
ess <- function(w) sum(w)^2 / sum(w^2)

# baseline + strict SOFA + date of death -> analysis cohort; G = grace period (hours)
prep_pop <- function(baseline, sofa_strict, dod, G = 24) {
  b <- baseline %>% group_by(subject_id) %>% slice_min(icu_intime, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(-any_of(c("sofa_t0_strict", "dod"))) %>% left_join(sofa_strict, by = "stay_id") %>% left_join(dod, by = "subject_id")
  bad <- (!is.na(b$crrt_first_time) & !is.na(b$deathtime) & b$crrt_first_time > b$deathtime) | b$t0 < b$icu_intime
  b <- b[!bad & !is.na(b$sofa_t0_strict), ] %>%
    mutate(d_in = as.numeric(difftime(deathtime, t0, units = "hours")), d_dod = as.numeric(dod - as.Date(t0)) * 24,
      death_h = case_when(!is.na(d_in) ~ d_in, !is.na(d_dod) ~ d_dod, TRUE ~ NA_real_),
      died28 = !is.na(death_h) & death_h >= 0 & death_h <= 672,
      death_h = if_else(died28, pmax(death_h, 0.5), NA_real_), end_h = if_else(died28, death_h, 672),
      icu_out_h = as.numeric(difftime(icu_outtime, t0, units = "hours")),
      h_crrt = if_else(!is.na(crrt_hours_from_t0) & crrt_hours_from_t0 >= 0 & crrt_hours_from_t0 < end_h, crrt_hours_from_t0, NA_real_),
      early = !is.na(h_crrt) & h_crrt < G,
      eq_strict = sofa_t0_strict >= 8 & sofa_t0_strict <= 14,
      era = factor(anchor_year_group))
  attr(b, "n_bad") <- sum(bad); b$id <- seq_len(nrow(b)); b
}

# tvl: stay_id, hour_index, Lk (the 1-hour lag, k = hour_index + 1, is applied inside build())
build <- function(pop, tvl, G = 24, kdigo0 = 3, era = FALSE, use_tv = TRUE, wlag = FALSE, amode = "L0") {
  # amode = "L0": early-arm weight = P(initiation) / P(initiation | L0) among patients alive at 24 h (landmark);
  # amode = "haz": uses the deferred-arm hourly initiation model with pre-initiation history: P(initiation by 24 h | history) = 1 - prod(1 - h_k)
  # wlag = TRUE: interval [k, k+1) uses the cumulative weight up to the start of hour k (prod_{j<k}); after 24 h the full product is used
  if (!"era" %in% names(pop)) pop$era <- factor(pop$anchor_year_group)
  pop <- med_imp(pop, L0)
  tvlag <- tvl %>% distinct(stay_id, hour_index, .keep_all = TRUE) %>% mutate(k = as.integer(hour_index) + 1L) %>% select(-hour_index)
  ph <- expand_grid(pop %>% select(id, stay_id, all_of(L0), creatinine_t0, end_h, icu_out_h, h_crrt, era), k = 0:(G - 1)) %>%
    filter(k < end_h, k < icu_out_h, is.na(h_crrt) | h_crrt >= k) %>%
    left_join(tvlag, by = c("stay_id", "k")) %>%
    mutate(map_h = if_else(k == 0, map_t0, map_h), nee_h = if_else(k == 0, nee_at_t0, nee_h),
           lactate_h = if_else(k == 0, lactate_t0, lactate_h), sofa_h = if_else(k == 0, sofa_t0_strict, sofa_h),
           kdigo_h = if_else(k == 0, as.numeric(kdigo0), as.numeric(kdigo_h)), creatinine_h = if_else(k == 0, creatinine_t0, creatinine_h),
           creatinine_delta = if_else(k == 0, 0, creatinine_delta), fluid_balance_h = if_else(k == 0, 0, fluid_balance_h),
           y = as.integer(!is.na(h_crrt) & h_crrt < k + 1))
  ph <- med_imp(ph, Lk)
  ex <- if (era) "era" else NULL
  mB <- suppressWarnings(glm(reformulate(c(L0, if (use_tv) Lk, "ns(k, 3)", ex), "y"), data = ph, family = binomial()))
  pd <- predict(mB, type = "response")
  pn <- predict(glm(reformulate(c("ns(k, 3)", ex), "y"), data = ph, family = binomial()), type = "response")
  ph <- ph %>% mutate(r = (1 - pn) / (1 - pd)) %>% group_by(id) %>% arrange(k, .by_group = TRUE) %>% mutate(Wfull = cumprod(r), W = if (wlag) lag(Wfull, default = 1) else Wfull) %>% ungroup()
  aG <- pop %>% filter(end_h > G)
  pA <- suppressWarnings(predict(glm(reformulate(c(L0, ex), "early"), data = aG, family = binomial()), type = "response"))
  nA <- if (era) predict(glm(early ~ era, data = aG, family = binomial()), type = "response") else mean(aG$early)
  aG$wA <- nA / pA
  if (amode == "haz") {
    ini <- aG$id[aG$early]
    hist <- ph %>% filter(id %in% ini)
    lastr <- hist %>% group_by(id) %>% slice_max(k, n = 1, with_ties = FALSE) %>% ungroup()
    ext <- lastr %>% rename(kmax = k) %>% crossing(k = 0:(G - 1)) %>% filter(k > kmax) %>% select(-kmax)
    allr <- bind_rows(hist, ext)
    allr$h <- predict(mB, newdata = allr, type = "response")
    PA <- allr %>% group_by(id) %>% summarise(PA = 1 - prod(1 - h), .groups = "drop")
    aG <- aG %>% left_join(PA, by = "id") %>% mutate(wA = if_else(early, mean(early) / PA, wA))
    attr(aG, "n_last_hour") <- sum(!is.na(aG$h_crrt) & aG$h_crrt >= G - 1 & aG$h_crrt < G)
  }
  key <- pop %>% select(id, end_h, died28, h_crrt, era)
  Bh <- ph %>% transmute(id, tstart = k, tstop = pmin(k + 1, h_crrt, end_h, na.rm = TRUE), w = W) %>% left_join(key, by = "id")
  lastW <- ph %>% group_by(id) %>% summarise(Wlast = last(Wfull), klast = max(k), .groups = "drop")
  Bp <- key %>% filter(is.na(h_crrt) | h_crrt >= G) %>% left_join(lastW, by = "id") %>%
    mutate(tstart = pmin(coalesce(klast, -1L) + 1, end_h), tstop = end_h, w = coalesce(Wlast, 1))
  Bcp <- bind_rows(Bh, Bp %>% select(id, tstart, tstop, w, end_h, died28, h_crrt, era)) %>% filter(tstop > tstart) %>%
    mutate(event = died28 & abs(tstop - end_h) < 1e-9 & (is.na(h_crrt) | h_crrt >= G), strategy = "Deferred")
  Acp <- bind_rows(key %>% transmute(id, tstart = 0, tstop = pmin(G, end_h), w = 1, end_h, died28, h_crrt, era),
                   aG %>% filter(early) %>% transmute(id, tstart = G, tstop = end_h, w = wA, end_h, died28, h_crrt, era)) %>%
    filter(tstop > tstart) %>% mutate(event = died28 & abs(tstop - end_h) < 1e-9, strategy = "Early")
  cp <- bind_rows(Acp, Bcp) %>% mutate(strategy = factor(strategy, c("Deferred", "Early")), tstart = tstart / 24, tstop = tstop / 24)
  # Protocol-specified truncation (primary analysis): cumulative stabilised weights truncated at the strategy-specific 99th percentile
  # (early strategy: weights after 24 h; deferred strategy: all intervals). Active when options(cp_trunc = TRUE); default is no truncation.
  cA <- NA_real_; cB <- NA_real_
  if (isTRUE(getOption("cp_trunc", FALSE))) {
    cA <- quantile(cp$w[cp$strategy == "Early" & cp$tstart > 0], .99, names = FALSE); cB <- quantile(cp$w[cp$strategy == "Deferred"], .99, names = FALSE)
    cp <- cp %>% mutate(w = if_else(strategy == "Early", if_else(tstart > 0, pmin(w, cA), w), pmin(w, cB)))
    aG$wA <- pmin(aG$wA, cA); ph$Wfull <- pmin(ph$Wfull, cB)
  }
  # Diagnostics: early-strategy clones still under follow-up at the end of the grace period
  pa <- aG %>% filter(early)
  smd <- sapply(L0, function(v) (weighted.mean(pa[[v]], pa$wA) - mean(aG[[v]])) / sd(aG[[v]]))
  attr(cp, "diag") <- c(n_A_post = nrow(pa), ESS_A = ess(pa$wA), max_absSMD_A = max(abs(smd)), max_wA = max(pa$wA), max_wB = max(ph$Wfull), cap_A = cA, cap_B = cB, n_init_last_hour = sum(!is.na(aG$h_crrt) & aG$h_crrt >= G - 1 & aG$h_crrt < G))
  cp
}

est <- function(cp, era = FALSE, adjL0 = FALSE, pop = NULL) {
  rhs <- c("strategy", if (adjL0) L0, if (era) "strata(era)")
  if (adjL0) cp <- cp %>% left_join(med_imp(pop, L0) %>% select(id, all_of(L0)), by = "id")
  m <- coxph(reformulate(rhs, "Surv(tstart, tstop, event)"), data = cp, weights = w, cluster = id, ties = "efron")
  b <- coef(m)[["strategyEarly"]]; se <- sqrt(vcov(m)["strategyEarly", "strategyEarly"])
  s <- summary(survfit(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w), rmean = 28)$table[, "rmean"]
  c(HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), logHR = b, SE = se,
    rmst_diff = unname(s[2] - s[1]))
}
