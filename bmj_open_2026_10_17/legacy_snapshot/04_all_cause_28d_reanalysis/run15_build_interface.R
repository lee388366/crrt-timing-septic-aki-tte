# RUN-15D: figure-runner 接口 CSV（列名/key 逐字沿用原图件脚本）
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble)})
W <- Sys.getenv("CRRT_WORKSPACE_ROOT"); R15 <- file.path(W, "09_工作稿/out/run15"); NE <- file.path(R15, "NEW")
strata <- read_csv(file.path(NE, "sofa_era/strictK3_v2_sofa_strata_results.csv"), show_col_types = FALSE)
scan   <- read_csv(file.path(NE, "sofa_era/strictK3_v2_sofa_band_scan_results.csv"), show_col_types = FALSE)
anch   <- read_csv(file.path(NE, "anchors_L0_era_rows.csv"), show_col_types = FALSE)
prim   <- read_csv(file.path(W, "08_CH01/out/ch01_primary_dod.csv"), show_col_types = FALSE)
# CH-10 改动（作者 2026-10-16 裁决）：era-stratified 行需要 early_n，改由此件取
#   口径来源: 09_工作稿/out/run15/NEW/sofa_era/strictK3_v2_calendar_era_results.csv
#             （CH-01 口径；逐年代含 early_n；5 层合计 early=146、n=1,502）
#   为什么按年代层求和: 合并(pooled)的 stratified 估计取自 anchors_L0_era_rows.csv，该件无 early_n 列
era    <- read_csv(file.path(NE, "sofa_era/strictK3_v2_calendar_era_results.csv"), show_col_types = FALSE)
era814 <- era %>% filter(grepl("^SOFA 8-14 era ", analysis))
stopifnot(nrow(era814) == 5, sum(era814$n_patients) == 1502)   # 与主行同人群，不成立即停
g <- function(x,a) x[[a]]
S <- function(l) strata %>% filter(analysis == l)
A <- function(l, pop = "SOFA 8-14") anch %>% filter(analysis == l, population == pop)

pe <- prim %>% filter(population == "Equipoise SOFA 8-14")
pa <- prim %>% filter(population == "Full cohort")

forest <- bind_rows(
  tibble(group="Primary / SOFA", analysis="Primary SOFA 8-14 equipoise", n_patients=1502, early_n=146,
         weighted_hr=g(pe,"hr"), ci_lower=g(pe,"ci_lower"), ci_upper=g(pe,"ci_upper"), p_value=g(pe,"p_value")),
  tibble(group="Primary / SOFA", analysis="Full strictK3_v2 primary cohort", n_patients=2463, early_n=195,
         weighted_hr=g(pa,"hr"), ci_lower=g(pa,"ci_lower"), ci_upper=g(pa,"ci_upper"), p_value=g(pa,"p_value")),
  tibble(group="SOFA strata", analysis="SOFA <=7", n_patients=g(S("SOFA <=7"),"n_patients"), early_n=g(S("SOFA <=7"),"early_n"),
         weighted_hr=g(S("SOFA <=7"),"weighted_hr"), ci_lower=g(S("SOFA <=7"),"ci_lower"), ci_upper=g(S("SOFA <=7"),"ci_upper"), p_value=g(S("SOFA <=7"),"p_value")),
  tibble(group="SOFA strata", analysis="SOFA 8-14", n_patients=g(S("SOFA 8-14"),"n_patients"), early_n=g(S("SOFA 8-14"),"early_n"),
         weighted_hr=g(S("SOFA 8-14"),"weighted_hr"), ci_lower=g(S("SOFA 8-14"),"ci_lower"), ci_upper=g(S("SOFA 8-14"),"ci_upper"), p_value=g(S("SOFA 8-14"),"p_value")),
  tibble(group="SOFA strata", analysis="SOFA 9-13", n_patients=g(S("SOFA 9-13"),"n_patients"), early_n=g(S("SOFA 9-13"),"early_n"),
         weighted_hr=g(S("SOFA 9-13"),"weighted_hr"), ci_lower=g(S("SOFA 9-13"),"ci_lower"), ci_upper=g(S("SOFA 9-13"),"ci_upper"), p_value=g(S("SOFA 9-13"),"p_value")),
  tibble(group="SOFA strata", analysis="SOFA >=15", n_patients=g(S("SOFA >=15"),"n_patients"), early_n=g(S("SOFA >=15"),"early_n"),
         weighted_hr=g(S("SOFA >=15"),"weighted_hr"), ci_lower=g(S("SOFA >=15"),"ci_lower"), ci_upper=g(S("SOFA >=15"),"ci_upper"), p_value=g(S("SOFA >=15"),"p_value")),
  tibble(group="Anchor sensitivity", analysis="KDIGO3_minus6h: SOFA 8-14", n_patients=1701, early_n=192,
         weighted_hr=g(A("legacy K3-6h anchor (NEW CH-01)"),"HR"), ci_lower=g(A("legacy K3-6h anchor (NEW CH-01)"),"ci_lower"),
         ci_upper=g(A("legacy K3-6h anchor (NEW CH-01)"),"ci_upper"), p_value=g(A("legacy K3-6h anchor (NEW CH-01)"),"p_value")),
  tibble(group="Anchor sensitivity", analysis="KDIGO2: SOFA 8-14", n_patients=2874, early_n=187,
         weighted_hr=g(A("KDIGO2 anchor (NEW CH-01)"),"HR"), ci_lower=g(A("KDIGO2 anchor (NEW CH-01)"),"ci_lower"),
         ci_upper=g(A("KDIGO2 anchor (NEW CH-01)"),"ci_upper"), p_value=g(A("KDIGO2 anchor (NEW CH-01)"),"p_value")),
  tibble(group="L0-only denominator", analysis="L0-only denominator: SOFA 8-14", n_patients=1502, early_n=146,
         weighted_hr=g(A("L0-only denominator (NEW CH-01)"),"HR"), ci_lower=g(A("L0-only denominator (NEW CH-01)"),"ci_lower"),
         ci_upper=g(A("L0-only denominator (NEW CH-01)"),"ci_upper"), p_value=g(A("L0-only denominator (NEW CH-01)"),"p_value"))
)
write_csv(forest, file.path(R15, "行值_forest_接口.csv"))

band <- scan %>% transmute(lower = lower, upper = upper, weighted_hr = weighted_hr,
                           n_patients = n_patients, early_n = early_n, p_value = p_value, ph_global_p = ph_global_p)
write_csv(band, file.path(R15, "行值_bandscan_接口.csv"))

era_strat <- tibble(
  model = "SOFA 8-14 era-stratified", n_patients = 1502,
  early_n = sum(era814$early_n),   # CH-10：由逐年代层求和（= 146），见文件头注释
  HR = g(A("Calendar-era stratified Cox (NEW CH-01)"),"HR"),
  CI_lo = g(A("Calendar-era stratified Cox (NEW CH-01)"),"ci_lower"),
  CI_hi = g(A("Calendar-era stratified Cox (NEW CH-01)"),"ci_upper"),
  p_value = g(A("Calendar-era stratified Cox (NEW CH-01)"),"p_value"),
  ph_global_p = g(A("Calendar-era stratified Cox (NEW CH-01)"),"ph_global_p"))
write_csv(era_strat, file.path(R15, "行值_era_stratified_接口.csv"))

cat("forest rows:", nrow(forest), "\n"); print(as.data.frame(forest %>% mutate(across(where(is.numeric), ~round(.x,4)))))
cat("band rows:", nrow(band), "\n"); cat("era-strat rows:", nrow(era_strat), "\n")
