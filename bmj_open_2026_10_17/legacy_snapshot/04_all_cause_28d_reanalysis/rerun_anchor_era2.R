suppressPackageStartupMessages({library(dplyr); library(survival)})
B <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
cox2 <- function(dat, lab, tcol, ecol, wcol="ipcw_w_trim", strata_var=NULL){
  f <- if (is.null(strata_var)) coxph(Surv(dat[[tcol]], dat[[ecol]]) ~ strategy_f, data=dat, weights=dat[[wcol]], cluster=stay_id, robust=TRUE, ties="efron")
       else coxph(as.formula(sprintf("Surv(%s,%s) ~ strategy_f + strata(%s)", tcol, ecol, strata_var)), data=dat, weights=dat[[wcol]], cluster=stay_id, robust=TRUE, ties="efron")
  co<-coef(f)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f)))[["strategy_fEarly CRRT"]]
  sprintf("%-30s N=%4d  HR=%.3f [%.3f, %.3f] p=%.3f", lab, length(unique(dat$stay_id)), exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se)))
}
cat("=== A) legacy 6h anchor：先复现旧口径 1.053（主人群）===\n")
a0 <- readRDS(file.path(B,"anchor_sensitivity_KDIGO3_minus6h/data/ccw_final.rds")) %>% as.data.frame()
a0$strategy_f <- factor(a0$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
cat(" ", cox2(a0[a0$equipoise,], "legacy6h 主人群（旧口径）","follow_days","event"), "\n")
cat(" ", cox2(a0, "legacy6h 全队列（旧口径）","follow_days","event"), "\n")
cat("  → 稿中锁定：主人群 1.053 (0.929-1.193)\n")

cat("\n=== B) calendar era（新口径）===\n")
d <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ccw_ch01_dod.rds"))
d$strategy_f <- factor(d$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
cat(" ", cox2(d, "era-stratified 全队列", "corrected_follow_days","corrected_event", strata_var="anchor_year_group"), "\n")
cat(" ", cox2(d, "era-adjusted 全队列 (covariate)", "corrected_follow_days","corrected_event"), "\n")
f2 <- coxph(Surv(corrected_follow_days,corrected_event) ~ strategy_f + anchor_year_group, data=d, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE)
co<-coef(f2)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f2)))[["strategy_fEarly CRRT"]]
cat(sprintf("  %-30s HR=%.3f [%.3f, %.3f] p=%.3f\n","era-as-covariate（字符串）", exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))))
d$era_f <- factor(d$anchor_year_group)
f3 <- coxph(Surv(corrected_follow_days,corrected_event) ~ strategy_f + era_f, data=d, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE)
co<-coef(f3)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f3)))[["strategy_fEarly CRRT"]]
cat(sprintf("  %-30s HR=%.3f [%.3f, %.3f] p=%.3f\n","era-adjusted（因子）", exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))))
cat("  旧口径参照：era-stratified 1.041 / era-adjusted 1.036，交互 p=0.85\n")
