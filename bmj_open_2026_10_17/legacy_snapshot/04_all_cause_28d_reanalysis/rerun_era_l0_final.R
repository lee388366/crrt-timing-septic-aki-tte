suppressPackageStartupMessages({library(dplyr); library(survival)})
B <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
con <- DBI::dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost", port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
dod <- DBI::dbGetQuery(con,"SELECT subject_id, dod FROM mimiciv_hosp.patients") %>% mutate(dod=as.Date(dod)); DBI::dbDisconnect(con)
mk <- function(x) x %>% distinct(clone_id,.keep_all=TRUE) %>% left_join(dod,by="subject_id") %>%
  mutate(t0=as.POSIXct(t0,tz="UTC"), deathtime=as.POSIXct(deathtime,tz="UTC"),
    di=as.numeric(difftime(deathtime,t0,units="days")), dd=as.numeric(dod-as.Date(t0)), dc=coalesce(di,dd),
    ev=!is.na(dc)&dc>=0&dc<=28, dh=dc*24,
    fo=pmin(28*24, coalesce(protocol_censor_h,28*24), if_else(ev,dh,28*24), na.rm=TRUE),
    ev=ev&(is.na(protocol_censor_h)|dh<protocol_censor_h), fo=fo/24,
    strategy_f=factor(strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))) %>% as.data.frame()

cat("=== 1) era（主人群，新口径）===\n")
d <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ccw_ch01_dod.rds"))
d$strategy_f <- factor(d$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
eq <- d[d$equipoise,]
f1 <- coxph(Surv(corrected_follow_days,corrected_event)~strategy_f+strata(anchor_year_group), data=eq, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE)
co<-coef(f1)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f1)))[["strategy_fEarly CRRT"]]
cat(sprintf("  %-26s HR=%.3f [%.3f, %.3f] p=%.3f\n","era-stratified 主人群", exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))))
eq$era_f<-factor(eq$anchor_year_group)
f2 <- coxph(Surv(corrected_follow_days,corrected_event)~strategy_f+era_f, data=eq, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE)
co<-coef(f2)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f2)))[["strategy_fEarly CRRT"]]
cat(sprintf("  %-26s HR=%.3f [%.3f, %.3f] p=%.3f\n","era-adjusted 主人群", exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))))
cat("  旧口径参照 1.041 / 1.036\n")

cat("\n=== 2) L0-only denominator（新口径）===\n")
L <- readRDS(file.path(B,"results/strictK3_v2_l0_only_sensitivity.rds"))
l0 <- mk(L$ccw_l0)
cat("  ccw_l0 列含 ipcw_w_trim:", "ipcw_w_trim" %in% colnames(l0), "\n")
run <- function(dat,lab){ f<-coxph(Surv(fo,ev)~strategy_f,data=dat,weights=ipcw_w_trim,cluster=stay_id,robust=TRUE,ties="efron")
  co<-coef(f)[["strategy_fEarly CRRT"]]; se<-sqrt(diag(vcov(f)))[["strategy_fEarly CRRT"]]
  sprintf("  %-30s N=%4d HR=%.3f [%.3f, %.3f] p=%.3f", lab, length(unique(dat$stay_id)), exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))) }
cat(run(l0,"L0-only 全队列"),"\n"); cat(run(l0[l0$equipoise,],"L0-only 主人群"),"\n")
cat("  旧口径参照：全队列 1.117 [1.002,1.245] p=.045 / 主人群 1.073 [0.936,1.230] p=.310\n")
