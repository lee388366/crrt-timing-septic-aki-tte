# =============================================================================
# RUN-15B: anchors (legacy K3-6h / KDIGO2) + L0-only + era-stratified
# 口径：CH-01 = 28 天全因死亡（deathtime 或 patients.dod，午夜日界）
# §1D：mk() 逐字取自 CH-01 rollout 的 rerun_anchor_era.R（同一 code path），
#      先用旧口径复现已报数值，再输出新口径；模型公式/权重/ties 一律不动。
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(survival); library(readr)})
B <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2")
OUT <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run15/NEW")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
canon <- readRDS(file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "08_CH01/out/ccw_ch01_dod.rds"))

con <- DBI::dbConnect(RPostgres::Postgres(), dbname = "mimiciv31", host = "localhost",
                      port = 5432, user = "postgres", password = Sys.getenv("PGPASSWORD"))
dod <- DBI::dbGetQuery(con, "SELECT subject_id, dod FROM mimiciv_hosp.patients") %>% mutate(dod = as.Date(dod))
DBI::dbDisconnect(con)

## ---- canonical CH-01 patch (verbatim from rerun_anchor_era.R) ----
mk <- function(x){
  x %>% distinct(clone_id, .keep_all = TRUE) %>% left_join(dod, by = "subject_id") %>%
  mutate(t0 = as.POSIXct(t0, tz = "UTC"), deathtime = as.POSIXct(deathtime, tz = "UTC"),
    di = as.numeric(difftime(deathtime, t0, units = "days")), dd = as.numeric(dod - as.Date(t0)),
    dc = coalesce(di, dd),
    ev = !is.na(dc) & dc >= 0 & dc <= 28, dh = dc * 24,
    fo = pmin(28 * 24, coalesce(protocol_censor_h, 28 * 24), if_else(ev, dh, 28 * 24), na.rm = TRUE),
    ev = ev & (is.na(protocol_censor_h) | dh < protocol_censor_h), fo = fo / 24,
    sf = factor(strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))) %>% as.data.frame()
}

cox1 <- function(dat, lab, pop, w = "ipcw_w_trim", strata_var = NULL, tcol = "fo", ecol = "ev") {
  f <- if (is.null(strata_var)) {
    coxph(as.formula(sprintf("Surv(%s,%s) ~ sf", tcol, ecol)), data = dat, weights = dat[[w]],
          cluster = stay_id, robust = TRUE, ties = "efron")
  } else {
    coxph(as.formula(sprintf("Surv(%s,%s) ~ sf + strata(%s)", tcol, ecol, strata_var)), data = dat,
          weights = dat[[w]], cluster = stay_id, robust = TRUE, ties = "efron")
  }
  co <- coef(f)[["sfEarly CRRT"]]; se <- sqrt(diag(vcov(f)))[["sfEarly CRRT"]]
  php <- tryCatch(cox.zph(f)$table["GLOBAL", "p"], error = function(e) NA_real_)
  tibble(analysis = lab, population = pop, n_patients = length(unique(dat$stay_id)),
         HR = exp(co), ci_lower = exp(co - 1.96 * se), ci_upper = exp(co + 1.96 * se),
         p_value = 2 * pnorm(-abs(co / se)), ph_global_p = php)
}

## ---- 0) 校验：mk() 是否等于 CH-01 canonical 数据集 ----
chk <- mk(readRDS(file.path(B, "data/ccw_final.rds")))
c2 <- canon %>% distinct(clone_id, .keep_all = TRUE) %>% arrange(clone_id)
chk <- chk %>% arrange(clone_id)
cat(sprintf("[check] mk() ev=%d (expect 1248) | max|fo-corrected_follow_days|=%.3g | ev==corrected_event: %s\n",
  sum(chk$ev), max(abs(chk$fo - c2$corrected_follow_days)),
  all(chk$ev == as.logical(c2$corrected_event))))

rows <- list()

## ---- 1) legacy K3-6h anchor ----
a_raw <- readRDS(file.path(B, "anchor_sensitivity_KDIGO3_minus6h/data/ccw_final.rds"))
a_raw$sf <- factor(a_raw$strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))
rows[["6h_old_eq"]]  <- cox1(as.data.frame(a_raw)[a_raw$equipoise %in% TRUE, ], "legacy K3-6h anchor (OLD dao)", "SOFA 8-14", tcol = "follow_days", ecol = "event")
rows[["6h_old_all"]] <- cox1(as.data.frame(a_raw), "legacy K3-6h anchor (OLD dao)", "Full cohort", tcol = "follow_days", ecol = "event")
a <- mk(a_raw)
rows[["6h_new_eq"]]  <- cox1(a[a$equipoise %in% TRUE, ], "legacy K3-6h anchor (NEW CH-01)", "SOFA 8-14")
rows[["6h_new_all"]] <- cox1(a, "legacy K3-6h anchor (NEW CH-01)", "Full cohort")

## ---- 2) KDIGO2 anchor ----
k_raw <- readRDS(file.path(B, "anchor_sensitivity_KDIGO2/data/ccw_final.rds"))
k_raw$sf <- factor(k_raw$strategy, levels = c("B_deferred", "A_early"), labels = c("Deferred CRRT", "Early CRRT"))
rows[["k2_old_eq"]]  <- cox1(as.data.frame(k_raw)[k_raw$equipoise %in% TRUE, ], "KDIGO2 anchor (OLD dao)", "SOFA 8-14", tcol = "follow_days", ecol = "event")
rows[["k2_old_all"]] <- cox1(as.data.frame(k_raw), "KDIGO2 anchor (OLD dao)", "Full cohort", tcol = "follow_days", ecol = "event")
k <- mk(k_raw)
rows[["k2_new_eq"]]  <- cox1(k[k$equipoise %in% TRUE, ], "KDIGO2 anchor (NEW CH-01)", "SOFA 8-14")
rows[["k2_new_all"]] <- cox1(k, "KDIGO2 anchor (NEW CH-01)", "Full cohort")

## ---- 3) L0-only denominator ----
L <- readRDS(file.path(B, "results/strictK3_v2_l0_only_sensitivity.rds"))
l0_old <- L$ccw_l0; l0_old$sf <- factor(l0_old$strategy, levels = c("B_deferred","A_early"), labels = c("Deferred CRRT","Early CRRT"))
rows[["l0_old_eq"]]  <- cox1(as.data.frame(l0_old)[l0_old$equipoise %in% TRUE, ], "L0-only denominator (OLD dao)", "SOFA 8-14", tcol = "follow_days", ecol = "event")
rows[["l0_old_all"]] <- cox1(as.data.frame(l0_old), "L0-only denominator (OLD dao)", "Full cohort", tcol = "follow_days", ecol = "event")
l0 <- mk(L$ccw_l0)
rows[["l0_new_eq"]]  <- cox1(l0[l0$equipoise %in% TRUE, ], "L0-only denominator (NEW CH-01)", "SOFA 8-14")
rows[["l0_new_all"]] <- cox1(l0, "L0-only denominator (NEW CH-01)", "Full cohort")

## ---- 4) calendar-era stratified (CH-01 canonical dataset) ----
dc <- canon %>% distinct(clone_id, .keep_all = TRUE) %>% as.data.frame()
dc$sf <- factor(dc$strategy, levels = c("B_deferred","A_early"), labels = c("Deferred CRRT","Early CRRT"))
rows[["era_strat_eq"]]  <- cox1(as.data.frame(dc)[dc$equipoise %in% TRUE, ], "Calendar-era stratified Cox (NEW CH-01)", "SOFA 8-14", strata_var = "anchor_year_group", tcol = "corrected_follow_days", ecol = "corrected_event")
rows[["era_strat_all"]] <- cox1(as.data.frame(dc), "Calendar-era stratified Cox (NEW CH-01)", "Full cohort", strata_var = "anchor_year_group", tcol = "corrected_follow_days", ecol = "corrected_event")

out <- bind_rows(rows)
write_csv(out, file.path(OUT, "anchors_L0_era_rows.csv"))
cat("\n=== rows ===\n")
print(as.data.frame(out %>% mutate(across(c(HR, ci_lower, ci_upper, p_value, ph_global_p), ~round(.x, 4)))), row.names = FALSE)
cat("\nSaved:", file.path(OUT, "anchors_L0_era_rows.csv"), "\n")
