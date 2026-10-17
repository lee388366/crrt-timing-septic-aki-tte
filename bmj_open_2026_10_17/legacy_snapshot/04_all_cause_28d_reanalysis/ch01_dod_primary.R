#!/usr/bin/env Rscript
# CH-01：主结局改为 28 天全因死亡（deathtime 或 patients.dod），不做出院删失
# 输入（只读）：strictK3_primary_v2/data/ccw_final.rds（干净队列 2,463，4,926 克隆行）
# 输出：本目录 out/
# 方法学来源：04_数据与真值/mmr_reanalysis/01_rerun_primary_with_dod.R（2026-06-02）
#   ⛔ 其数字跑在污染队列（3,112）上，属禁用清单；本脚本只复用其口径。
suppressPackageStartupMessages({ library(DBI); library(RPostgres); library(dplyr); library(survival) })
set.seed(20260602)

OUT <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ccw <- readRDS(file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data/ccw_final.rds")) %>%
  distinct(clone_id, .keep_all = TRUE)
cat("输入克隆行（clone_id 去重后）:", nrow(ccw), " 唯一 ICU stay:", length(unique(ccw$stay_id)), "\n")

con <- dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost", port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
on.exit(dbDisconnect(con), add = TRUE)
dod <- dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients") %>% mutate(dod = as.Date(dod))

ccw2 <- ccw %>% left_join(dod, by="subject_id") %>%
  mutate(
    deathtime = as.POSIXct(deathtime, tz="UTC"),
    t0 = as.POSIXct(t0, tz="UTC"),
    death_days_inhospital = as.numeric(difftime(deathtime, t0, units="days")),
    death_days_dod        = as.numeric(dod - as.Date(t0)),
    death_days_corrected  = case_when(!is.na(death_days_inhospital) ~ death_days_inhospital,
                                      !is.na(death_days_dod) ~ death_days_dod, TRUE ~ NA_real_),
    died_28d_corrected = !is.na(death_days_corrected) & death_days_corrected >= 0 & death_days_corrected <= 28,
    death_h_corrected  = death_days_corrected * 24,
    admin_h = 28 * 24,
    corrected_follow_h = pmin(admin_h, coalesce(protocol_censor_h, admin_h),
                              if_else(died_28d_corrected, death_h_corrected, admin_h), na.rm = TRUE),
    corrected_event = died_28d_corrected & (is.na(protocol_censor_h) | death_h_corrected < protocol_censor_h),
    corrected_follow_days = corrected_follow_h / 24,
    strategy_f = factor(strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
  ) %>% as.data.frame()

# ---- K3：dod 覆盖与补获 ----
cat("\ndod 覆盖（患者级，n=", nrow(ccw2)/2, "）：\n")
pt <- ccw2 %>% distinct(stay_id, .keep_all=TRUE)
cat("  有 deathtime:", sum(!is.na(pt$death_days_inhospital)), "\n")
cat("  有 dod      :", sum(!is.na(pt$death_days_dod)), "\n")
cat("  两者都无    :", sum(is.na(pt$death_days_inhospital) & is.na(pt$death_days_dod)), "\n")
rep_tab <- pt %>% mutate(old=as.logical(died_28d), new=died_28d_corrected, newly = new & !old) %>%
  group_by(early_crrt) %>% summarise(n=n(), old_deaths=sum(old,na.rm=TRUE),
    new_deaths=sum(new,na.rm=TRUE), newly_captured=sum(newly,na.rm=TRUE), .groups="drop")
cat("\n按实际治疗组的死亡补获：\n"); print(as.data.frame(rep_tab))
write.csv(rep_tab, file.path(OUT,"mortality_ascertainment_by_group.csv"), row.names=FALSE)

# ---- 主分析：Cox（IPCW 加权，按 stay_id 聚类）----
rmst_wkm <- function(times, events, weights, tau=28){
  fit <- survfit(Surv(times, events) ~ 1, weights=weights)
  tt <- c(0, fit$time); ss <- c(1, fit$surv)
  if (max(tt) < tau) { tt <- c(tt, tau); ss <- c(ss, ss[length(ss)]) }
  else { idx <- max(which(tt <= tau)); tt <- c(tt[seq_len(idx)], tau); ss <- c(ss[seq_len(idx)], ss[idx]) }
  keep <- !duplicated(tt); tt <- tt[keep]; ss <- ss[keep]
  sum(diff(tt) * ss[-length(ss)])
}
analyse <- function(dat, label){
  fit <- coxph(Surv(corrected_follow_days, corrected_event) ~ strategy_f, data=dat,
               weights=ipcw_w_trim, cluster=stay_id, robust=TRUE, ties="efron")
  co <- coef(fit)[["strategy_fEarly CRRT"]]; se <- sqrt(diag(vcov(fit)))[["strategy_fEarly CRRT"]]
  e <- dat[dat$strategy_f=="Early CRRT",]; d <- dat[dat$strategy_f=="Deferred CRRT",]
  data.frame(population=label, n_patients=length(unique(dat$stay_id)), clone_rows=nrow(dat),
    events=sum(dat$corrected_event), hr=exp(co), ci_lower=exp(co-1.96*se), ci_upper=exp(co+1.96*se),
    p_value=2*pnorm(-abs(co/se)),
    rmst_early=rmst_wkm(e$corrected_follow_days,e$corrected_event,e$ipcw_w_trim),
    rmst_deferred=rmst_wkm(d$corrected_follow_days,d$corrected_event,d$ipcw_w_trim),
    schoenfeld_global_p=cox.zph(fit)$table["GLOBAL","p"], row.names=NULL)
}
res <- bind_rows(analyse(ccw2, "Full cohort"), analyse(filter(ccw2, equipoise), "Equipoise SOFA 8-14"))
res$rmst_difference <- res$rmst_early - res$rmst_deferred
write.csv(res, file.path(OUT,"ch01_primary_dod.csv"), row.names=FALSE)
saveRDS(ccw2, file.path(OUT,"ccw_ch01_dod.rds"))
cat("\n=== CH-01 主结果（28 天全因死亡，不做出院删失）===\n"); print(res)
