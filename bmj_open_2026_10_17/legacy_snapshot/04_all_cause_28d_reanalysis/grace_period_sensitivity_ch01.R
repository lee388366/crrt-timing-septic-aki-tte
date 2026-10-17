# =============================================================================
# Q6 (reviewer) — Grace-period length sensitivity analysis
#
# Reproduces the locked v2 CCW/IPCW pipeline (03a + 03_ccw_ipcw.R, GRACE=24)
# and re-runs it under alternative grace-period windows (12 h, 24 h, 48 h).
#
# Authoritative v2 inputs (clean first-KDIGO3 t0 cohort, N=2,463):
#   strictK3_primary_v2/data/{baseline.rds, tv_hourly.rds, cohort_base.rds}
#
# GRACE<=24 : uses the existing tv_hourly (hours 0-23) subset -> no DB.
# GRACE>24  : requires an hours 0-(GRACE-1) tv_hourly rebuilt from the DB
#             (built by grace_period_build_tv48.R and passed via --tv48=...).
#
# The clone-censoring, IPCW denominators/numerators, ridge specification,
# 99th-pct weight trimming, and weighted Cox mirror 03_ccw_ipcw.R exactly;
# only the grace window (and the hour-(GRACE-1) cross-section for Strategy A,
# and the "early = CRRT within GRACE h" definition) are parameterized.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lubridate); library(glmnet); library(survival); library(readr); library(DBI); library(RPostgres)
})

V2_DIR   <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data")
OUT_DIR  <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out")

L0_vars <- c("age","sofa_t0","apsiii","creatinine_baseline",
             "nee_at_t0","lactate_t0","uo_rate_t0","map_t0","cci")
Lk_vars <- c("map_h","nee_h","lactate_h","sofa_h",
             "kdigo_h","creatinine_h","creatinine_delta","fluid_balance_h")

impute_median <- function(df, vars) {
  for (v in vars) { med <- median(df[[v]], na.rm=TRUE); df[[v]][is.na(df[[v]])] <- med }
  df
}

# ------------------------------------------------------------------ core CCW/IPCW
build_ccw <- function(GRACE, baseline, tv_hourly) {
  stopifnot(max(tv_hourly$hour_index) >= (GRACE - 1))
  cohort <- baseline %>%
    group_by(subject_id) %>% slice_min(icu_intime, n=1, with_ties=FALSE) %>% ungroup() %>%
    filter(!is.na(sofa_t0))
  cohort <- impute_median(cohort, L0_vars)

  # Clone-censor (identical structure to 03_ccw_ipcw.R; early redefined as <=GRACE)
  ccw <- cohort %>%
    crossing(strategy = c("A_early","B_deferred")) %>%
    mutate(
      clone_id = paste(stay_id, strategy, sep="_"),
      h_crrt   = crrt_hours_from_t0,
      early_g  = !is.na(h_crrt) & h_crrt <= GRACE,
      protocol_censor_h = case_when(
        strategy=="A_early"    & (is.na(h_crrt) | h_crrt > GRACE) ~ as.double(GRACE),
        strategy=="B_deferred" & !is.na(h_crrt) & h_crrt <= GRACE ~ h_crrt,
        TRUE ~ NA_real_),
      protocol_censored = !is.na(protocol_censor_h),
      follow_h    = pmin(follow_time_28d*24, coalesce(protocol_censor_h, follow_time_28d*24)),
      follow_days = follow_h/24,
      event = died_28d & (is.na(protocol_censor_h) | follow_time_28d*24 < protocol_censor_h)
    )

  # ---- Strategy A: cross-sectional at hour (GRACE-1) ----
  tv_hourA <- tv_hourly %>% filter(hour_index == (GRACE-1)) %>% distinct(stay_id, .keep_all=TRUE)
  A_data <- ccw %>% filter(strategy=="A_early") %>%
    left_join(tv_hourA, by="stay_id") %>%
    mutate(y_A = as.integer(early_g))
  A_data <- impute_median(A_data, Lk_vars)
  X_A <- model.matrix(~ age+sofa_t0+apsiii+creatinine_baseline+nee_at_t0+lactate_t0+uo_rate_t0+map_t0+cci+
                        map_h+nee_h+lactate_h+sofa_h+kdigo_h+creatinine_h+creatinine_delta+fluid_balance_h,
                      data=A_data)[,-1]
  set.seed(42)
  cv_A <- cv.glmnet(X_A, A_data$y_A, family="binomial", alpha=0, nfolds=10, type.measure="deviance")
  p_denom_A <- as.vector(predict(cv_A, X_A, s="lambda.1se", type="response"))
  lm_num_A <- glm(y_A ~ age+sofa_t0+apsiii+creatinine_baseline+nee_at_t0+lactate_t0+uo_rate_t0+map_t0+cci,
                  data=A_data, family=binomial())
  p_num_A <- predict(lm_num_A, type="response")
  A_data <- A_data %>% mutate(
    ipcw_w = if_else(early_g, p_num_A/pmax(p_denom_A,0.01), (1-p_num_A)/pmax(1-p_denom_A,0.01)))

  # ---- Strategy B: pooled person-hours 0..(GRACE-1) ----
  B_base <- ccw %>% filter(strategy=="B_deferred")
  B_long <- B_base %>% select(stay_id, early_g, h_crrt, all_of(L0_vars)) %>%
    left_join(tv_hourly %>% filter(hour_index <= (GRACE-1)), by="stay_id", relationship="many-to-many") %>%
    mutate(crrt_this_hour = !is.na(h_crrt) & h_crrt>=hour_index & h_crrt<(hour_index+1) & h_crrt<=GRACE,
           at_risk = is.na(h_crrt) | h_crrt > hour_index) %>%
    filter(at_risk)
  B_long <- impute_median(B_long, Lk_vars)
  X_B <- model.matrix(~ age+sofa_t0+apsiii+creatinine_baseline+nee_at_t0+lactate_t0+uo_rate_t0+map_t0+cci+hour_index+
                        map_h+nee_h+lactate_h+sofa_h+kdigo_h+creatinine_h+creatinine_delta+fluid_balance_h,
                      data=B_long)[,-1]
  y_B <- as.integer(B_long$crrt_this_hour)
  set.seed(42)
  cv_B <- cv.glmnet(X_B, y_B, family="binomial", alpha=0, nfolds=10, type.measure="deviance")
  p_denom_B <- as.vector(predict(cv_B, X_B, s="lambda.1se", type="response"))
  lm_num_B <- glm(crrt_this_hour ~ age+sofa_t0+apsiii+creatinine_baseline+nee_at_t0+lactate_t0+uo_rate_t0+map_t0+cci+hour_index,
                  data=B_long, family=binomial())
  p_num_B <- predict(lm_num_B, type="response")
  B_long <- B_long %>% mutate(w_contrib=(1-p_num_B)/pmax(1-p_denom_B,0.01))
  B_weights <- B_long %>% group_by(stay_id) %>% summarise(ipcw_w=prod(w_contrib,na.rm=TRUE), .groups="drop")

  # ---- trim @ 99th pct ----
  w99_A <- quantile(A_data$ipcw_w,0.99,na.rm=TRUE); w99_B <- quantile(B_weights$ipcw_w,0.99,na.rm=TRUE)
  A_data <- A_data %>% mutate(ipcw_w_trim=pmin(ipcw_w,w99_A))
  B_weights <- B_weights %>% mutate(ipcw_w_trim=pmin(ipcw_w,w99_B))

  A_w <- A_data %>% distinct(stay_id,.keep_all=TRUE) %>% select(stay_id,ipcw_w,ipcw_w_trim)
  ccw_A <- ccw %>% filter(strategy=="A_early") %>% left_join(A_w, by="stay_id", relationship="many-to-one")
  ccw_B <- ccw %>% filter(strategy=="B_deferred") %>% left_join(B_weights, by="stay_id", relationship="many-to-one") %>%
    mutate(ipcw_w=coalesce(ipcw_w,1), ipcw_w_trim=coalesce(ipcw_w_trim,1))
  bind_rows(ccw_A, ccw_B) %>%
    mutate(equipoise = !is.na(sofa_t0) & sofa_t0>=8 & sofa_t0<=14) %>%
    distinct(clone_id, .keep_all=TRUE) %>% mutate(trt=as.integer(strategy=="A_early"))
}

run_cox <- function(dat, GRACE, pop) {
  dat <- dat %>% filter(!is.na(follow_days), !is.na(event), !is.na(ipcw_w_trim))
  f <- coxph(Surv(follow_days,event)~trt, data=dat, weights=ipcw_w_trim, cluster=stay_id, robust=TRUE, ties="efron")
  s <- summary(f)
  tibble(grace_h=GRACE, population=pop,
         n_patients=n_distinct(dat$stay_id),
         early_n=n_distinct(dat$stay_id[dat$early_g %in% TRUE]),
         events=sum(dat$event,na.rm=TRUE),
         weighted_hr=unname(s$coefficients["trt","exp(coef)"]),
         ci_lower=unname(s$conf.int["trt","lower .95"]),
         ci_upper=unname(s$conf.int["trt","upper .95"]),
         p_value=unname(s$coefficients["trt","Pr(>|z|)"]),
         ph_global_p=tryCatch(cox.zph(f)$table["GLOBAL","p"],error=function(e)NA_real_))
}

run_grace <- function(GRACE, baseline, tv_hourly) {
  ccw <- build_ccw(GRACE, baseline, tv_hourly)
  bind_rows(run_cox(ccw, GRACE, "Full cohort"),
            run_cox(filter(ccw, equipoise %in% TRUE), GRACE, "Primary SOFA 8-14"))
}

# ------------------------------------------------------------------ driver
args <- commandArgs(trailingOnly=TRUE)
tv48_path <- sub("^--tv48=","", grep("^--tv48=", args, value=TRUE))
graces <- as.integer(sub("^--grace=","", grep("^--grace=", args, value=TRUE)))
if (length(graces)==0) graces <- c(12L, 24L)

baseline <- readRDS(file.path(V2_DIR,"baseline.rds"))

# ===== CH-01 补丁（2026-10-15）：主结局改为 28 天全因死亡 =====
# 原口径：died_28d = 院内死亡（出院即删失）；原 follow_time_28d = min(至死亡, 至出院, 28)
# 新口径：28 天内全因死亡 = deathtime 或 patients.dod（不做出院删失）；随访至死亡或 28 天
# dod 仅到日粒度 ⇒ 采用午夜日界规则（与 CH-01 主分析一致）
.con <- DBI::dbConnect(RPostgres::Postgres(), dbname="mimiciv31", host="localhost",
                       port=5432, user="postgres", password = Sys.getenv("PGPASSWORD"))
.dod <- DBI::dbGetQuery(.con, "SELECT subject_id, dod FROM mimiciv_hosp.patients")
DBI::dbDisconnect(.con)
.dod$dod <- as.Date(.dod$dod)
baseline <- baseline %>%
  left_join(.dod, by="subject_id") %>%
  mutate(
    t0        = as.POSIXct(t0, tz="UTC"),
    deathtime = as.POSIXct(deathtime, tz="UTC"),
    .di = as.numeric(difftime(deathtime, t0, units="days")),
    .dd = as.numeric(dod - as.Date(t0)),
    .dc = coalesce(.di, .dd),
    died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
    follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
  ) %>%
  select(-.di, -.dd, -.dc)
cat(sprintf("[CH-01 patch] 新口径死亡数 = %d / %d\n", sum(baseline$died_28d), nrow(baseline)))
# ===== 补丁结束 =====
tv_hourly <- readRDS(file.path(V2_DIR,"tv_hourly.rds"))

res <- list()
for (g in graces) {
  if (g <= 24) {
    res[[as.character(g)]] <- run_grace(g, baseline, tv_hourly)
  } else {
    if (length(tv48_path)==0 || !file.exists(tv48_path))
      stop(sprintf("GRACE=%d needs --tv48=<tv_hourly_0_%d.rds>", g, g-1))
    tv_big <- readRDS(tv48_path)
    res[[as.character(g)]] <- run_grace(g, baseline, tv_big)
  }
  cat("done grace =", g, "\n"); flush.console()
}
out <- bind_rows(res) %>% arrange(population, grace_h)
# RUN-22b（lead 2026-10-16 批准）：唯一输出名，防止被 ab_ipcw_*/ab2_*/ab3_* 系列（共用旧名路径）覆盖
write_csv(out, file.path(OUT_DIR, "grace_period_sensitivity_ch01_results.csv"))
print(as.data.frame(out), digits=4)