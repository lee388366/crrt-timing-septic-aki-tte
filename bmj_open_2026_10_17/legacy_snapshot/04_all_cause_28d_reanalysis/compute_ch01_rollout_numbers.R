suppressPackageStartupMessages({library(dplyr); library(survival)})
d <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ccw_ch01_dod.rds"))
d$strategy_f <- factor(d$strategy, levels=c("B_deferred","A_early"), labels=c("Deferred CRRT","Early CRRT"))
cat("=== 1) 粗（未加权）Cox —— 新主结局 ===\n")
for (pop in c("full","equipoise")) {
  dt <- if (pop=="full") d else d[d$equipoise,]
  fit <- coxph(Surv(corrected_follow_days, corrected_event) ~ strategy_f, data=dt, cluster=stay_id, robust=TRUE, ties="efron")
  co <- coef(fit)[["strategy_fEarly CRRT"]]; se <- sqrt(diag(vcov(fit)))[["strategy_fEarly CRRT"]]
  cat(sprintf("  %-12s crude HR=%.3f (%.3f, %.3f)  p=%.3f  events=%d  n=%d\n",
      pop, exp(co), exp(co-1.96*se), exp(co+1.96*se), 2*pnorm(-abs(co/se)), sum(dt$corrected_event), length(unique(dt$stay_id))))
}
cat("\n=== 2) Table 1 行：患者级，按实际治疗组（主人群）===\n")
pt <- d %>% distinct(stay_id, .keep_all=TRUE) %>% filter(sofa_equipoise)
tab <- pt %>% group_by(early_crrt) %>%
  summarise(n=n(), deaths=sum(died_28d_corrected), pct=100*sum(died_28d_corrected)/n(), .groups="drop")
print(as.data.frame(tab))
p_early <- tab$pct[tab$early_crrt]; p_def <- tab$pct[!tab$early_crrt]
smd <- abs(p_early-p_def)/sqrt((p_early*(100-p_early)+p_def*(100-p_def))/2)
cat(sprintf("\n  合计 %d (%.1f%%) | SMD=%.3f\n", sum(tab$deaths), 100*sum(tab$deaths)/nrow(pt), smd))
m <- matrix(c(sum(tab$deaths[tab$early_crrt]), tab$n[tab$early_crrt]-sum(tab$deaths[tab$early_crrt]),
              sum(tab$deaths[!tab$early_crrt]), tab$n[!tab$early_crrt]-sum(tab$deaths[!tab$early_crrt])), nrow=2, byrow=TRUE)
cat("  chi-square p =", sprintf("%.3f", chisq.test(m)$p.value), "\n")
cat("\n  旧口径同法（对照，应得 720 / 47.9% / SMD 0.106 / p .22）：\n")
tab2 <- pt %>% group_by(early_crrt) %>% summarise(n=n(), d=sum(as.logical(died_28d)), pct=100*sum(as.logical(died_28d))/n(), .groups="drop")
print(as.data.frame(tab2))
p1<-tab2$pct[tab2$early_crrt]; p0<-tab2$pct[!tab2$early_crrt]
cat(sprintf("  合计 %d (%.1f%%) SMD=%.3f\n", sum(tab2$d), 100*sum(tab2$d)/nrow(pt),
    abs(p1-p0)/sqrt((p1*(100-p1)+p0*(100-p0))/2)))
cat("\n=== 3) 全队列患者级（补充 Table S6 用）===\n")
ptf <- d %>% distinct(stay_id, .keep_all=TRUE)
cat(sprintf("  旧 %d (%.1f%%) → 新 %d (%.1f%%)\n", sum(as.logical(ptf$died_28d)), 100*mean(as.logical(ptf$died_28d)),
    sum(ptf$died_28d_corrected), 100*mean(ptf$died_28d_corrected)))
