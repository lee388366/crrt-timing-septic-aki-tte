# =============================================================================
# crude_cox_v2.R
#   把旧母板 (12_primary_outcome.R) 的 crude / IPCW Cox 写法，
#   原样适配到 strictK3 v2 clone 数据 (ccw_final.rds) 上重跑。
#   不重新原创方法，仅把旧代码指向 v2 数据。
#
#   4 个人群 × 2 个模型 (crude 未加权 / IPCW 加权)：
#     a. 全队列 (整份数据)
#     b. 主分析 SOFA 8-14 (equipoise %in% TRUE)
#     c. SOFA<9 (sofa_t0 < 9)
#     d. SOFA>=9 (sofa_t0 >= 9)
#   全部 cluster = stay_id, robust = TRUE
#   HR 取 strategy 的 A_early (参照 B_deferred)
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
})

# ── Paths ────────────────────────────────────────────────────────────────────
v2_data <- file.path(Sys.getenv("CRRT_ANALYSIS_ROOT"), "strictK3_primary_v2/data/ccw_final.rds")
out_dir <- file.path(Sys.getenv("CRRT_MANUSCRIPT_ROOT"), "11_母板重建_数据回补/结果")
out_csv <- file.path(out_dir, "crude_cox_v2_results.csv")

# ── Load v2 CCW data ──────────────────────────────────────────────────────────
d <- readRDS(v2_data)

# strategy 以 B_deferred 为参照 (relevel)；保留 v2 原始水平名 A_early / B_deferred
d$strategy <- factor(d$strategy, levels = c("B_deferred", "A_early"))

cat(sprintf("v2 CCW dataset: %d rows, %d unique stay_ids\n",
            nrow(d), length(unique(d$stay_id))))
cat(sprintf("Events: %d / %d (%.1f%%)\n",
            sum(d$event), nrow(d), 100 * mean(d$event)))
cat(sprintf("equipoise TRUE: %d rows\n\n", sum(d$equipoise %in% TRUE)))

# ── Helper: extract Cox results (get_cox_res 风格) ────────────────────────────
# HR = exp(coef); 95%CI = exp(coef ± 1.96*se); p = 2*pnorm(-|z|)
get_cox_res <- function(fit) {
  co  <- coef(fit)
  se  <- sqrt(diag(vcov(fit)))
  idx <- which(names(co) == "strategyA_early")
  hr  <- exp(co[idx])
  lo  <- exp(co[idx] - 1.96 * se[idx])
  hi  <- exp(co[idx] + 1.96 * se[idx])
  pv  <- 2 * pnorm(-abs(co[idx] / se[idx]))
  c(HR = unname(hr), CI_lower = unname(lo), CI_upper = unname(hi), p_value = unname(pv))
}

# ── Helper: 对某一人群跑 crude + IPCW 两版，返回 2 行 data.frame ──────────────
run_population <- function(dat, pop_label) {
  n_stays <- length(unique(dat$stay_id))
  events  <- sum(dat$event)

  # crude 未加权
  fit_crude <- coxph(Surv(follow_days, event) ~ strategy,
                     data = dat, cluster = stay_id, robust = TRUE)
  r_crude <- get_cox_res(fit_crude)

  # IPCW 加权
  fit_ipcw <- coxph(Surv(follow_days, event) ~ strategy,
                    data = dat, weights = ipcw_w_trim,
                    cluster = stay_id, robust = TRUE)
  r_ipcw <- get_cox_res(fit_ipcw)

  cat(sprintf("── %s (n_stays=%d, events=%d) ──\n", pop_label, n_stays, events))
  cat(sprintf("   crude HR = %.3f [%.3f, %.3f]  p=%.4f\n",
              r_crude["HR"], r_crude["CI_lower"], r_crude["CI_upper"], r_crude["p_value"]))
  cat(sprintf("   IPCW  HR = %.3f [%.3f, %.3f]  p=%.4f\n\n",
              r_ipcw["HR"], r_ipcw["CI_lower"], r_ipcw["CI_upper"], r_ipcw["p_value"]))

  data.frame(
    population = pop_label,
    model      = c("crude", "IPCW"),
    n_stays    = n_stays,
    events     = events,
    HR         = c(r_crude["HR"], r_ipcw["HR"]),
    CI_lower   = c(r_crude["CI_lower"], r_ipcw["CI_lower"]),
    CI_upper   = c(r_crude["CI_upper"], r_ipcw["CI_upper"]),
    p_value    = c(r_crude["p_value"], r_ipcw["p_value"]),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# ── 4 个人群 ──────────────────────────────────────────────────────────────────
d_full  <- d
d_eq    <- d[d$equipoise %in% TRUE, ]
d_lt9   <- d[d$sofa_t0 <  9, ]
d_ge9   <- d[d$sofa_t0 >= 9, ]

res <- rbind(
  run_population(d_full, "Full cohort"),
  run_population(d_eq,   "Primary SOFA 8-14 (equipoise)"),
  run_population(d_lt9,  "SOFA<9"),
  run_population(d_ge9,  "SOFA>=9")
)

# ── 输出 CSV ──────────────────────────────────────────────────────────────────
res_out <- res
res_out$HR       <- round(res_out$HR, 3)
res_out$CI_lower <- round(res_out$CI_lower, 3)
res_out$CI_upper <- round(res_out$CI_upper, 3)
res_out$p_value  <- round(res_out$p_value, 4)
write.csv(res_out, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("✓ Results CSV saved: %s\n\n", out_csv))

# ── 自检: IPCW 锚点复现 ───────────────────────────────────────────────────────
anchor_check <- function(pop, expect_hr, expect_lo, expect_hi, tol = 0.01) {
  row <- res[res$population == pop & res$model == "IPCW", ]
  ok  <- abs(row$HR - expect_hr) <= tol &&
         abs(row$CI_lower - expect_lo) <= tol &&
         abs(row$CI_upper - expect_hi) <= tol
  cat(sprintf("[anchor] %-32s IPCW HR=%.3f [%.3f, %.3f]  expect %.3f [%.3f, %.3f]  -> %s\n",
              pop, row$HR, row$CI_lower, row$CI_upper,
              expect_hr, expect_lo, expect_hi, ifelse(ok, "PASS", "FAIL")))
  ok
}

cat(strrep("=", 60), "\n")
cat("IPCW anchor self-check\n")
cat(strrep("=", 60), "\n")
ok1 <- anchor_check("Primary SOFA 8-14 (equipoise)", 1.044, 0.922, 1.181)
ok2 <- anchor_check("Full cohort",                   1.095, 0.993, 1.207)
if (ok1 && ok2) {
  cat("\n✓ ALL ANCHORS REPRODUCED\n")
} else {
  cat("\n✗ ANCHOR MISMATCH — check script/data\n")
}
