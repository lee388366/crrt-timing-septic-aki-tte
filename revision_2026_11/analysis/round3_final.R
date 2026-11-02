# Round 3: final inference (independent audit items M2, m5); with options(cp_trunc = TRUE) via round7 this gives the primary analysis
# Specification: round-1 inputs (strictly lagged hourly SOFA/KDIGO, baseline values at or before t0) + weight at interval start + L0 without APS III + early-arm L0 landmark weights
# Output: HR, weighted 28-day risks and risk difference, RMST, secondary outcomes; 1,000 bootstrap resamples (both weight models refitted) + 200-group jackknife for BCa
# Usage: Rscript round3_final.R    (resumable: completed batches are skipped)
source("cp_core.R")
suppressPackageStartupMessages({library(parallel)})
PKG <- Sys.getenv("CRRT_REVISION_DIR")
WS  <- Sys.getenv("CRRT_INPUT_DIR")
OUT <- file.path(PKG, "out"); PRIV <- file.path(PKG, "private_data"); BOOT <- file.path(PRIV, "boot"); dir.create(BOOT, showWarnings = FALSE)
sink(file.path(OUT, "round3_report.txt"), split = TRUE); cat("Round3 run at", format(Sys.time()), "\n\n")
NREP <- 1000; NJK <- 200; CORES <- 4; RNGkind("L'Ecuyer-CMRG")

L0  <- setdiff(L0, "apsiii")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds")); tvl <- readRDS(file.path(PRIV, "tv_hourly_round1.rds"))
aud <- readRDS(file.path(WS, "ch24_q4_corrected_secondary/NEW/run13/ANALYSIS_DF_AUDIT.rds")) %>%
  distinct(stay_id, .keep_all = TRUE) %>% select(stay_id, renal_recovery_28d, rfd_28d)

stat <- function(pop, keep = FALSE) {
  cp <- build(pop, tvl, wlag = TRUE, amode = "L0")
  m  <- coxph(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w, ties = "efron")
  sf <- survfit(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w)
  rm <- summary(sf, rmean = 28)$table[, "rmean"]; s28 <- summary(sf, times = 28, extend = TRUE)$surv
  fin <- cp %>% group_by(strategy, id) %>% slice_max(tstop, n = 1, with_ties = FALSE) %>% ungroup() %>%
    filter(abs(tstop - end_h / 24) < 1e-9) %>% left_join(pop %>% select(id, stay_id), by = "id") %>% left_join(aud, by = "stay_id")
  sm <- fin %>% group_by(strategy) %>% summarise(rec = weighted.mean(renal_recovery_28d %in% TRUE, w), rfd = weighted.mean(rfd_28d, w, na.rm = TRUE), .groups = "drop")
  r <- c(logHR = unname(coef(m)), rmst_def = unname(rm[1]), rmst_early = unname(rm[2]), rmst_diff = unname(rm[2] - rm[1]),
    risk_def = 1 - s28[1], risk_early = 1 - s28[2], rd = (1 - s28[2]) - (1 - s28[1]),
    rec_def = sm$rec[1], rec_early = sm$rec[2], rec_diff = sm$rec[2] - sm$rec[1],
    rfd_def = sm$rfd[1], rfd_early = sm$rfd[2], rfd_diff = sm$rfd[2] - sm$rfd[1])
  if (keep) attr(r, "cp") <- cp
  r
}
reid <- function(p) { p$id <- seq_len(nrow(p)); p }
bca <- function(th, bt, jk, alpha = c(.025, .975)) {
  bt <- bt[is.finite(bt)]; jk <- jk[is.finite(jk)]
  z0 <- qnorm(mean(bt < th) + 0.5 * mean(bt == th)); jm <- mean(jk)
  a <- sum((jm - jk)^3) / (6 * sum((jm - jk)^2)^1.5); z <- qnorm(alpha)
  quantile(bt, pnorm(z0 + (z0 + z) / (1 - a * (z0 + z))), names = FALSE)
}
run_chunks <- function(f, idx, fun) {
  done <- if (file.exists(f)) readRDS(f) else NULL
  todo <- setdiff(idx, if (is.null(done)) integer(0) else done[, "rep"])
  for (ch in split(todo, ceiling(seq_along(todo) / 100))) {
    r <- do.call(rbind, mclapply(ch, function(b) c(rep = b, tryCatch(fun(b), error = function(e) rep(NA_real_, 13))), mc.cores = CORES))
    done <- rbind(done, r); saveRDS(done, f); cat(format(Sys.time(), "%H:%M:%S"), basename(f), nrow(done), "/", length(idx), "\n")
  }
  done
}
priors <- list(c("Data-dominant", 0, 1, "primary"), c("Neutral", 0, .355, "primary"), c("Skeptical", 0, .205, "primary"),
  c("Mildly optimistic", -.105, .205, "sensitivity"), c("Mildly pessimistic", .095, .205, "sensitivity"),
  c("Discounted-historical neutral", 0, .25, "sensitivity"), c("ELAIN-informed stress test", -.223, .3, "sensitivity"))
evalue <- function(hr) { if (hr < 1) hr <- 1 / hr; rr <- (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr)); rr + sqrt(rr * (rr - 1)) }

allrows <- list()
for (sub in c("primary", "full")) {
  pop <- reid(if (sub == "primary") filter(coh, eq_strict) else coh)
  cat("==================== ", sub, " (n=", nrow(pop), ", early=", sum(pop$early), ", 28-day deaths=", sum(pop$died28), ") ====================\n", sep = "")
  pe <- stat(pop, keep = TRUE); cp <- attr(pe, "cp"); d <- attr(cp, "diag")
  rob <- coxph(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w, cluster = id, ties = "efron")
  seR <- sqrt(vcov(rob)[1, 1])
  postB <- cp %>% filter(strategy == "Deferred", tstart <= 1, tstop > 1) %>% distinct(id, .keep_all = TRUE)
  cat(sprintf("Diagnostics: early arm after 24 h n=%d, ESS %.1f, max|SMD| %.3f, max w %.2f | deferred arm at 24 h n=%d, ESS %.1f, max w %.2f\n",
      d["n_A_post"], d["ESS_A"], d["max_absSMD_A"], d["max_wA"], nrow(postB), ess(postB$w), d["max_wB"]))
  sf <- survfit(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w)
  ss <- summary(sf, times = seq(0, 28, by = 0.25), extend = TRUE)
  write.csv(data.frame(strategy = sub("strategy=", "", ss$strata), time = ss$time, surv = ss$surv, n_risk_weighted = ss$n.risk),
            file.path(OUT, paste0("round3_km_curve_", sub, ".csv")), row.names = FALSE)

  bt <- run_chunks(file.path(BOOT, paste0("round3_boot_", sub, ".rds")), 1:NREP, function(b) { set.seed(20261102 + b)
    stat(reid(pop[sample.int(nrow(pop), replace = TRUE), ])) })
  grp <- { set.seed(7); sample(rep(seq_len(NJK), length.out = nrow(pop))) }
  jk <- run_chunks(file.path(BOOT, paste0("round3_jk_", sub, ".rds")), 1:NJK, function(g) stat(reid(pop[grp != g, ])))

  rep_line <- function(lab, key, tf = identity, fmt = "%.3f") {
    th <- pe[[key]]; x <- bt[, key]; j <- jk[, key]; ok <- is.finite(x)
    ci_b <- bca(th, x, j); ci_p <- quantile(x[ok], c(.025, .975), names = FALSE); se <- sd(x[ok])
    ci_n <- c(th - 1.96 * se, th + 1.96 * se)
    f <- function(v) sprintf(fmt, tf(v))
    cat(sprintf("  %-28s %s | BCa %s to %s | percentile %s to %s | normal %s to %s | valid %d\n", lab, f(th), f(ci_b[1]), f(ci_b[2]),
        f(ci_p[1]), f(ci_p[2]), f(ci_n[1]), f(ci_n[2]), sum(ok)))
    allrows[[length(allrows) + 1]] <<- data.frame(population = sub, measure = lab, estimate = tf(th), bca_lo = tf(ci_b[1]), bca_hi = tf(ci_b[2]),
      pct_lo = tf(ci_p[1]), pct_hi = tf(ci_p[2]), norm_lo = tf(ci_n[1]), norm_hi = tf(ci_n[2]), n_boot = sum(ok))
    ci_b
  }
  cat("Primary outcome (28-day all-cause mortality)\n")
  ciHR <- rep_line("HR (early vs deferred)", "logHR", exp)
  cat(sprintf("  (robust-variance CI: %.3f-%.3f, for comparison only)\n", exp(pe[["logHR"]] - 1.96 * seR), exp(pe[["logHR"]] + 1.96 * seR)))
  rep_line("28-day risk, deferred", "risk_def", function(v) 100 * v, "%.1f%%"); rep_line("28-day risk, early", "risk_early", function(v) 100 * v, "%.1f%%")
  rep_line("28-day risk difference (percentage points)", "rd", function(v) 100 * v, "%.1f")
  rep_line("RMST deferred (days)", "rmst_def", identity, "%.2f"); rep_line("RMST early (days)", "rmst_early", identity, "%.2f")
  rep_line("RMST difference (days)", "rmst_diff", identity, "%.2f")
  cat("Secondary outcomes (fixed at day 28)\n")
  rep_line("Renal recovery, deferred", "rec_def", function(v) 100 * v, "%.1f%%"); rep_line("Renal recovery, early", "rec_early", function(v) 100 * v, "%.1f%%")
  rep_line("Renal recovery difference (percentage points)", "rec_diff", function(v) 100 * v, "%.1f")
  rep_line("RRT-free days, deferred", "rfd_def", identity, "%.2f"); rep_line("RRT-free days, early", "rfd_early", identity, "%.2f")
  rep_line("RRT-free days difference", "rfd_diff", identity, "%.2f")
  x <- bt[is.finite(bt[, "logHR"]), "logHR"]
  cat(sprintf("bootstrap logHR: median HR %.3f, skewness %.2f, SD %.4f (robust SE %.4f)\n", exp(median(x)), mean((x - mean(x))^3) / sd(x)^3, sd(x), seR))

  # Bayesian (normal likelihood; SE derived from the BCa interval width) and E-values
  se_b <- (ciHR[2] - ciHR[1]) / (2 * 1.959964); est_l <- pe[["logHR"]]   # ciHR is the BCa interval on the log-HR scale
  cat(sprintf("Bayesian (logHR %.4f, SE from BCa %.4f):\n", est_l, se_b))
  for (p in priors) { m0 <- as.numeric(p[2]); s0 <- as.numeric(p[3]); v <- 1 / (1 / se_b^2 + 1 / s0^2); mu <- v * (est_l / se_b^2 + m0 / s0^2)
    b <- pnorm(log(.9), mu, sqrt(v)); h <- 1 - pnorm(log(1.1), mu, sqrt(v))
    cat(sprintf("  %-30s %-11s posterior median HR %.3f | benefit %4.1f%% | equivalence %4.1f%% | harm %4.1f%%\n", p[1], p[4], exp(mu), 100 * b, 100 * (1 - b - h), 100 * h))
    allrows[[length(allrows) + 1]] <- data.frame(population = sub, measure = paste("Bayes", p[1]), estimate = exp(mu), bca_lo = b, bca_hi = h,
      pct_lo = NA, pct_hi = NA, norm_lo = NA, norm_hi = NA, n_boot = NA) }
  cat(sprintf("E-value (common-outcome approximation): point %.2f; CI limit %s\n\n", evalue(exp(est_l)),
      if (ciHR[1] <= 0 && ciHR[2] >= 0) "1.00 (CI includes 1)" else sprintf("%.2f", evalue(if (est_l < 0) exp(ciHR[2]) else exp(ciHR[1])))))
}
write.csv(bind_rows(allrows), file.path(OUT, "round3_final_estimates.csv"), row.names = FALSE)
cat("In Bayes rows, bca_lo / bca_hi hold the posterior probabilities of benefit / harm\nDONE\n"); sink()
