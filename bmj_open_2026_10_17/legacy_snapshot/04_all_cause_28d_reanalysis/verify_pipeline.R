suppressPackageStartupMessages({library(dplyr); library(survival)})
ccw <- readRDS(file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data/ccw_final.rds")) %>% distinct(clone_id, .keep_all=TRUE) %>% as.data.frame()
ccw$strategy_f <- factor(ccw$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
cat("=== 管线校验：用 ccw_final.rds 里已有的 follow_days/event 复现锁定值 ===\n")
cat("（预期：全队列 1.095 (0.993-1.207) p=0.069；主人群 1.044 (0.922-1.181) p=0.497）\n\n")
f <- function(d, lab){
  fit <- coxph(Surv(follow_days, event) ~ strategy_f, data=d, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE, ties="efron")
  co <- coef(fit)[["strategy_fEarly CRRT"]]; se <- sqrt(diag(vcov(fit)))[["strategy_fEarly CRRT"]]
  cat(sprintf("%-24s n=%4d  events=%4d  HR=%.4f (%.4f-%.4f)  p=%.4f\n",
      lab, length(unique(d$stay_id)), sum(d$event), exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se))))
}
f(ccw, "全队列 follow_days")
f(ccw[ccw$equipoise,], "主人群 follow_days")

cat("\n=== 逐患者：主人群，旧口径 vs dod 补正 ===\n")
con <- DBI::dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost", port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
dod <- DBI::dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients") %>% mutate(dod=as.Date(dod))
DBI::dbDisconnect(con)
d <- ccw %>% left_join(dod, by="subject_id") %>%
  mutate(t0=as.POSIXct(t0,tz="UTC"), deathtime=as.POSIXct(deathtime,tz="UTC"),
         dd_inh=as.numeric(difftime(deathtime,t0,units="days")), dd_dod=as.numeric(dod-as.Date(t0)),
         dd=coalesce(dd_inh, dd_dod),
         new28 = !is.na(dd) & dd>=0 & dd<=28) %>%
  distinct(stay_id, .keep_all=TRUE)
for (pop in c("full","equipoise")) {
  x <- if (pop=="full") d else d[d$sofa_equipoise,]
  tab <- x %>% mutate(old=as.logical(died_28d), new=new28) %>%
    group_by(early_crrt) %>% summarise(n=n(), old=sum(old), new=sum(new), added=sum(new & !old), .groups="drop")
  cat("\n--- ", pop, " (n=", nrow(x), ") ---\n", sep=""); print(as.data.frame(tab))
  cat("  合计 old=", sum(tab$old), " new=", sum(tab$new), " 新增=", sum(tab$added),
      " 旧死亡率=", sprintf("%.1f%%",100*sum(tab$old)/nrow(x)),
      " 新死亡率=", sprintf("%.1f%%",100*sum(tab$new)/nrow(x)), "\n", sep="")
}
