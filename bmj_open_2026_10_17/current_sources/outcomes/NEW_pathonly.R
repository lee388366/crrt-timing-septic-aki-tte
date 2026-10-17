suppressPackageStartupMessages({library(dplyr); library(survival)})
d <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_primary/ccw_ch01_dod.rds"))
d$strategy_f <- factor(d$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
run <- function(dat, lab){
  f <- coxph(Surv(corrected_follow_days, corrected_event) ~ strategy_f, data=dat,
             weights=ipcw_w_trim, cluster=stay_id, robust=TRUE, ties="efron")
  co <- coef(f)[["strategy_fEarly CRRT"]]; se <- sqrt(diag(vcov(f)))[["strategy_fEarly CRRT"]]
  data.frame(analysis=lab, N=length(unique(dat$stay_id)),
             early_events=sum(dat$corrected_event[dat$strategy_f=="Early CRRT"]),
             HR=exp(co), lo=exp(co-1.96*se), hi=exp(co+1.96*se), p=2*pnorm(-abs(co/se)), row.names=NULL)
}
mc <- median(d$creatinine_t0, na.rm=TRUE); ma <- median(d$apsiii, na.rm=TRUE)
cat("中位数: creatinine =", round(mc,3), " APS III =", ma, "\n\n")
res <- bind_rows(
 run(filter(d, equipoise), "Primary: SOFA at t0 8-14"),
 run(d, "Full cohort"),
 run(filter(d, sofa_band_prespecified=="SOFA_9_13"), "Conservative SOFA 9-13"),
 run(filter(d, sofa_band_prespecified=="SOFA_lt9"), "SOFA <9"),
 run(filter(d, sofa_band_prespecified %in% c("SOFA_9_13","SOFA_gt13")), "SOFA >=9"),
 run(filter(d, creatinine_t0 <  mc), "Creatinine < median"),
 run(filter(d, creatinine_t0 >= mc), "Creatinine >= median"),
 run(filter(d, apsiii <  ma), "APS III < median"),
 run(filter(d, apsiii >= ma), "APS III >= median")
)
res$HR <- round(res$HR,3); res$lo <- round(res$lo,3); res$hi <- round(res$hi,3); res$p <- round(res$p,3)
print(as.data.frame(res), row.names=FALSE)
write.csv(res, file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/ch24_q4_corrected_tableS4/NEW/tableS4_ch01.csv"), row.names=FALSE)
cat("\n（旧口径对照：主人群 1.044[0.922,1.181] p=.497 / 全队列 1.095[0.993,1.207] p=.069 / 9-13 1.044[0.900,1.212] p=.568 / <9 1.069 / >=9 1.043 / Cr< 1.121 / Cr>= 1.063 / APS< 1.034 / APS>= 1.065）\n")
