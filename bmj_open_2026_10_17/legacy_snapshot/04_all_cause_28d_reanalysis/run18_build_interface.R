# RUN-18 接口：把加强指标整理成可直接入表的形式
suppressPackageStartupMessages({library(dplyr); library(readr)})
R18 <- file.path(Sys.getenv("CRRT_WORKSPACE_ROOT"), "09_工作稿/out/run18")
x <- read_csv(file.path(R18, "NEW/strictK3_v2_rfd_extended_metrics.csv"), show_col_types = FALSE)
d <- read_csv(file.path(R18, "NEW/strictK3_v2_rfd_restricted_mean_diff.csv"), show_col_types = FALSE)

tbl <- x %>%
  transmute(
    population, framing, arm, n = n_patients,
    zero_n, zero_pct = round(zero_pct_unweighted, 1),
    zero_pct_weighted = round(zero_pct_weighted, 1),
    ess = round(ess, 1),
    restricted_mean = round(restricted_mean, 2),
    rmean_lo = round(rmean_lo, 2), rmean_hi = round(rmean_hi, 2),
    rmean_display = sprintf("%.2f (%.2f–%.2f)", restricted_mean, rmean_lo, rmean_hi),
    p25, p50, p75, p90
  )
diff_tbl <- d %>% transmute(population, n = n_patients, weighted_mean_diff, diff_lo, diff_hi, p_value,
                            diff_display = sprintf("%.3f (%.3f–%.3f), p=%s", weighted_mean_diff, diff_lo, diff_hi,
                              ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))))
write_csv(tbl, file.path(R18, "行值_RFD加强指标_入表.csv"))
write_csv(diff_tbl, file.path(R18, "行值_RFD受限均值差_入表.csv"))
cat("=== 入表 CSV 1（RFD 加强指标）===\n"); print(as.data.frame(tbl), row.names = FALSE)
cat("\n=== 入表 CSV 2（受限均值差 = 加权均值差）===\n"); print(as.data.frame(diff_tbl), row.names = FALSE)
