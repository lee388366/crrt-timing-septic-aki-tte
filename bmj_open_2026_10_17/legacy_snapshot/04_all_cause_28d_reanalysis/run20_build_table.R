# RUN-20 组装：v1/v2 并排表（旧口径 vs 新口径）+ 污染数 + 权重诊断
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble); library(DBI); library(RPostgres)})
SA <- Sys.getenv("CRRT_ANALYSIS_ROOT"); R20 <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run20")
old_v1 <- read_csv(file.path(R20, "OLD/v1/results/strictK3_primary_mortality_results.csv"), show_col_types = FALSE)
new_v1 <- read_csv(file.path(R20, "NEW/v1/results/strictK3_primary_mortality_results.csv"), show_col_types = FALSE)
old_v2 <- read_csv(file.path(SA, "strictK3_primary_v2/results/strictK3_primary_v2_mortality_results.csv"), show_col_types = FALSE)
ch01   <- read_csv(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ch01_primary_dod.csv"), show_col_types = FALSE)

t0sum <- read_csv(file.path(SA, "strictK3_primary/results/strictK3_QA_t0_current_summary.csv"), show_col_types = FALSE)
t0arm <- read_csv(file.path(SA, "strictK3_primary/results/strictK3_QA_t0_current_by_arm.csv"), show_col_types = FALSE)
contam <- tibble(
  metric = c("v1 污染版主人群 early 数", "v1 污染版主人群 N", "v1 early 占比 %",
             "v1 early 中 treatment-defined t0 数", "v1 early 中 treatment-defined t0 占比 %",
             "v1 主人群 treatment-defined t0 合计（early+deferred）",
             "v2 干净版主人群 early 数", "v2 干净版主人群 N", "v2 early 占比 %"),
  value = c(t0sum$early_n, t0sum$n, round(100*t0sum$early_n/t0sum$n, 2),
            t0arm$t0_by_crrt_n[t0arm$arm=="early"],
            round(100*t0arm$t0_by_crrt_n[t0arm$arm=="early"]/t0arm$n[t0arm$arm=="early"], 2),
            t0sum$t0_by_crrt_n, 146, 1502, round(100*146/1502, 2))
)
write_csv(contam, file.path(R20, "污染数_v1v2.csv"))

wdiag <- function(ccw, version) {
  caps <- ccw %>% group_by(strategy) %>% summarise(cap = unname(quantile(ipcw_w, 0.99, na.rm = TRUE)), .groups="drop")
  d <- ccw %>% left_join(caps, by = "strategy")
  d %>% group_by(strategy) %>%
    summarise(version = version, cap_p99_full_arm = first(cap),
      trimmed_pct_full_arm = 100*mean(ipcw_w > cap),
      n_equipoise = sum(equipoise %in% TRUE),
      trimmed_pct_equipoise = 100*mean(ipcw_w[equipoise %in% TRUE] > cap[equipoise %in% TRUE]),
      ess_equipoise = { w <- ipcw_w_trim[equipoise %in% TRUE]; sum(w)^2/sum(w^2) },
      mean_w_equipoise = mean(ipcw_w_trim[equipoise %in% TRUE]),
      max_w_equipoise = max(ipcw_w_trim[equipoise %in% TRUE]),
      .groups = "drop")
}
w1 <- readRDS(file.path(SA, "strictK3_primary/data/ccw_final.rds")) %>% distinct(clone_id, .keep_all = TRUE)
w2 <- readRDS(file.path(SA, "strictK3_primary_v2/data/ccw_final.rds")) %>% distinct(clone_id, .keep_all = TRUE)
wtab <- bind_rows(wdiag(w1, "v1 contaminated t0"), wdiag(w2, "v2 clean t0"))
write_csv(wtab, file.path(R20, "权重诊断_ESS_trimming.csv"))

pick <- function(df, a, idcol = "analysis") df %>% filter(.data[[idcol]] == a) %>% slice(1)
mk <- function(pop, oldrow, newrow, n, early_n, ev_old, ev_new) {
  tibble(population = pop, N = n, early_n = early_n, events_old = ev_old, events_new = ev_new,
    HR_old = oldrow$weighted_hr, HR_old_lo = oldrow$ci_lower, HR_old_hi = oldrow$ci_upper, HR_old_p = oldrow$p_value,
    HR_new = newrow$weighted_hr, HR_new_lo = newrow$ci_lower, HR_new_hi = newrow$ci_upper, HR_new_p = newrow$p_value,
    RMST_old = oldrow$rmst_diff, RMST_old_lo = oldrow$rmst_diff_ci_lower, RMST_old_hi = oldrow$rmst_diff_ci_upper,
    RMST_new = newrow$rmst_diff, RMST_new_lo = newrow$rmst_diff_ci_lower, RMST_new_hi = newrow$rmst_diff_ci_upper,
    dao_old = "旧口径（院内死亡 + 出院删失）", dao_new = "新口径 CH-01（28 天全因死亡，不做出院删失）")
}
tab <- bind_rows(
  mk("污染版主人群 (v1, SOFA 8-14)", pick(old_v1,"Primary SOFA 8-14 equipoise"), pick(new_v1,"Primary SOFA 8-14 equipoise"), 1568, 215, 873, 944),
  mk("污染版全队列 (v1)", pick(old_v1,"Full strictK3 primary cohort"), pick(new_v1,"Full strictK3 primary cohort"), 2566, 299, 1365, 1486)
)
write_csv(tab, file.path(R20, "行值_v1v2并排表.csv"))

# 干净版参考行：HR/CI/p 与 RMST 差取自 CH-01 canonical 产物；新口径 RMST CI 来自 08_CH01/out/rmst_sens.txt（bootstrap B=1000）
ce <- pick(ch01, "Equipoise SOFA 8-14", "population"); cf <- pick(ch01, "Full cohort", "population")
oe <- pick(old_v2, "Primary SOFA 8-14 equipoise"); of <- old_v2 %>% filter(grepl("^Full", analysis)) %>% slice(1)
clean <- bind_rows(
  tibble(population="干净版主人群 (v2, SOFA 8-14)", N=1502, early_n=146, events_old=836, events_new=ce$events,
    HR_old=oe$weighted_hr, HR_old_lo=oe$ci_lower, HR_old_hi=oe$ci_upper, HR_old_p=oe$p_value,
    HR_new=ce$hr, HR_new_lo=ce$ci_lower, HR_new_hi=ce$ci_upper, HR_new_p=ce$p_value,
    RMST_old=oe$rmst_diff, RMST_old_lo=oe$rmst_diff_ci_lower, RMST_old_hi=oe$rmst_diff_ci_upper,
    RMST_new=-1.05, RMST_new_lo=-3.23, RMST_new_hi=1.06,
    dao_old="旧口径", dao_new="新口径 CH-01（RMST CI 出处 08_CH01/out/rmst_sens.txt）"),
  tibble(population="干净版全队列 (v2)", N=2463, early_n=195, events_old=1309, events_new=cf$events,
    HR_old=of$weighted_hr, HR_old_lo=of$ci_lower, HR_old_hi=of$ci_upper, HR_old_p=of$p_value,
    HR_new=cf$hr, HR_new_lo=cf$ci_lower, HR_new_hi=cf$ci_upper, HR_new_p=cf$p_value,
    RMST_old=of$rmst_diff, RMST_old_lo=of$rmst_diff_ci_lower, RMST_old_hi=of$rmst_diff_ci_upper,
    RMST_new=-1.90, RMST_new_lo=-3.74, RMST_new_hi=-0.13,
    dao_old="旧口径", dao_new="新口径 CH-01（RMST CI 出处 08_CH01/out/rmst_sens.txt）")
)
write_csv(clean, file.path(R20, "行值_v1v2并排表_干净版参考.csv"))

cat("=== 污染数 ===\n"); print(as.data.frame(contam), row.names = FALSE)
cat("\n=== 权重诊断 ===\n"); print(as.data.frame(wtab %>% mutate(across(where(is.numeric), ~round(.x,4)))), row.names = FALSE)
cat("\n=== 并排表（污染版）===\n"); print(as.data.frame(tab %>% transmute(population, N, early_n, events_old, events_new,
  HR_old = sprintf("%.3f (%.3f-%.3f) p=%.3f", HR_old, HR_old_lo, HR_old_hi, HR_old_p),
  HR_new = sprintf("%.3f (%.3f-%.3f) p=%.3f", HR_new, HR_new_lo, HR_new_hi, HR_new_p),
  RMST_old = sprintf("%.2f (%.2f to %.2f)", RMST_old, RMST_old_lo, RMST_old_hi),
  RMST_new = sprintf("%.2f (%.2f to %.2f)", RMST_new, RMST_new_lo, RMST_new_hi))), row.names = FALSE)
cat("\n=== 并排表（干净版参考）===\n"); print(as.data.frame(clean %>% transmute(population, N, early_n, events_old, events_new,
  HR_old = sprintf("%.3f (%.3f-%.3f) p=%.3f", HR_old, HR_old_lo, HR_old_hi, HR_old_p),
  HR_new = sprintf("%.3f (%.3f-%.3f) p=%.3f", HR_new, HR_new_lo, HR_new_hi, HR_new_p),
  RMST_old = sprintf("%.2f (%.2f to %.2f)", RMST_old, RMST_old_lo, RMST_old_hi),
  RMST_new = sprintf("%.2f (%.2f to %.2f)", RMST_new, RMST_new_lo, RMST_new_hi))), row.names = FALSE)