suppressPackageStartupMessages({library(dplyr); library(survival)})
set.seed(20260602)
OUT <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out")
d <- readRDS(file.path(OUT,"ccw_ch01_dod.rds"))
rmst_wkm <- function(times, events, weights, tau=28){
  fit <- survfit(Surv(times, events) ~ 1, weights=weights)
  tt <- c(0, fit$time); ss <- c(1, fit$surv)
  if (max(tt) < tau) { tt <- c(tt, tau); ss <- c(ss, ss[length(ss)]) }
  else { idx <- max(which(tt <= tau)); tt <- c(tt[seq_len(idx)], tau); ss <- c(ss[seq_len(idx)], ss[idx]) }
  keep <- !duplicated(tt); tt <- tt[keep]; ss <- ss[keep]
  sum(diff(tt) * ss[-length(ss)])
}
boot_rmst <- function(dat, b=1000, tau=28){
  ids <- unique(dat$stay_id); n <- length(ids); diffs <- numeric(b)
  for (i in seq_len(b)) {
    s <- data.frame(stay_id=sample(ids,n,replace=TRUE))
    bo <- merge(s, dat, by="stay_id", all.x=TRUE)
    e <- bo[bo$strategy_f=="Early CRRT",]; df <- bo[bo$strategy_f=="Deferred CRRT",]
    diffs[i] <- rmst_wkm(e$corrected_follow_days,e$corrected_event,e$ipcw_w_trim,tau) -
                rmst_wkm(df$corrected_follow_days,df$corrected_event,df$ipcw_w_trim,tau)
  }
  c(mean=mean(diffs,na.rm=TRUE), lower=unname(quantile(diffs,0.025,na.rm=TRUE)),
    upper=unname(quantile(diffs,0.975,na.rm=TRUE)))
}
cat("=== CH-01 RMST 差（bootstrap B=1000）===\n")
for (pop in c("full","equipoise")) {
  dt <- if (pop=="full") d else d[d$equipoise,]
  b <- boot_rmst(dt)
  cat(sprintf("%-20s RMST差 = %.2f 天 (95%% CI %.2f to %.2f)\n", pop, b[["mean"]], b[["lower"]], b[["upper"]]))
}

# ---- 敏感性：dod 仅有时点的规则（午夜日界 vs 12:00 规则）----
cat("\n=== 敏感性：dod 时点规则 ===\n")
d2 <- d
d2$death_days_dod_noon <- as.numeric(as.Date(d2$dod) - as.Date(d2$t0)) + 0.5
d2$dd_noon <- ifelse(!is.na(d2$death_days_inhospital), d2$death_days_inhospital,
              ifelse(!is.na(d2$death_days_dod_noon), d2$death_days_dod_noon, NA))
d2$ev_noon <- !is.na(d2$dd_noon) & d2$dd_noon >= 0 & d2$dd_noon <= 28
d2$death_h_noon <- d2$dd_noon * 24
d2$fo_noon <- pmin(28*24, ifelse(is.na(d2$protocol_censor_h), 28*24, d2$protocol_censor_h),
                   ifelse(d2$ev_noon, d2$death_h_noon, 28*24), na.rm=TRUE)
d2$ev_noon <- d2$ev_noon & (is.na(d2$protocol_censor_h) | d2$death_h_noon < d2$protocol_censor_h)
d2$fo_noon_days <- d2$fo_noon/24
for (pop in c("full","equipoise")) {
  dt <- if (pop=="full") d2 else d2[d2$equipoise,]
  fit <- coxph(Surv(fo_noon_days, ev_noon) ~ strategy_f, data=dt, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE, ties="efron")
  co <- coef(fit)[["strategy_fEarly CRRT"]]; se <- sqrt(diag(vcov(fit)))[["strategy_fEarly CRRT"]]
  cat(sprintf("%-20s 12:00 规则 HR=%.4f (%.4f-%.4f) p=%.4f | 事件=%d\n", pop, exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se)), sum(dt$ev_noon)))
}
