suppressPackageStartupMessages({library(dplyr); library(survival)})
d <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ccw_ch01_dod.rds"))
rmst_wkm <- function(times, events, weights, tau=28){
  fit <- survfit(Surv(times, events) ~ 1, weights=weights)
  tt <- c(0, fit$time); ss <- c(1, fit$surv)
  if (max(tt) < tau) { tt <- c(tt, tau); ss <- c(ss, ss[length(ss)]) }
  else { idx <- max(which(tt <= tau)); tt <- c(tt[seq_len(idx)], tau); ss <- c(ss[seq_len(idx)], ss[idx]) }
  keep <- !duplicated(tt); tt <- tt[keep]; ss <- ss[keep]
  sum(diff(tt) * ss[-length(ss)])
}
for (pop in c("full","equipoise")) {
  dt <- if (pop=="full") d else d[d$equipoise,]
  e <- dt[dt$strategy=="A_early",]; df <- dt[dt$strategy=="B_deferred",]
  re <- rmst_wkm(e$corrected_follow_days, e$corrected_event, e$ipcw_w_trim)
  rd <- rmst_wkm(df$corrected_follow_days, df$corrected_event, df$ipcw_w_trim)
  cat(sprintf("%-10s Early %.2f  Deferred %.2f  Diff %.2f\n", pop, re, rd, re-rd))
}
cat("\n旧值（稿件）：全队列 15.61 / 17.24 / -1.63；主人群 16.07 / 16.74 / -0.66\n")
