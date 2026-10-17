# RUN-15C: 组装 Fig3 / Fig S3 行值 CSV（每行含 HR/95%CI/p + 来源脚本 + 口径）
suppressPackageStartupMessages({library(dplyr); library(readr); library(tibble)})
W  <- Sys.getenv("CRRT_WORKSPACE_ROOT")
R15 <- file.path(W, "09_工作稿/out/run15")
QA  <- file.path(W, "09_工作稿/out/qa_validation")
NE  <- file.path(R15, "NEW")

old3   <- read_csv(file.path(QA, "strictK3_v2_fig3_forest_display.csv"), show_col_types = FALSE)
strata <- read_csv(file.path(NE, "sofa_era/strictK3_v2_sofa_strata_results.csv"), show_col_types = FALSE)
scan   <- read_csv(file.path(NE, "sofa_era/strictK3_v2_sofa_band_scan_results.csv"), show_col_types = FALSE)
scan_o <- read_csv(file.path(QA, "strictK3_v2_sofa_band_scan_results.csv"), show_col_types = FALSE)
era    <- read_csv(file.path(NE, "sofa_era/strictK3_v2_calendar_era_results.csv"), show_col_types = FALSE)
erao   <- read_csv(file.path(QA, "strictK3_v2_calendar_era_results.csv"), show_col_types = FALSE)
eraadj <- read_csv(file.path(NE, "sofa_era/strictK3_v2_calendar_era_adjusted_results.csv"), show_col_types = FALSE)
eraadjo<- read_csv(file.path(QA, "strictK3_v2_calendar_era_adjusted_results.csv"), show_col_types = FALSE)
anch   <- read_csv(file.path(NE, "anchors_L0_era_rows.csv"), show_col_types = FALSE)
prim   <- read_csv(file.path(W, "08_CH01/out/ch01_primary_dod.csv"), show_col_types = FALSE)

g <- function(x, a) x[[a]]
pick_strata <- function(lab) strata %>% filter(analysis == lab)
pick_anch <- function(lab, pop) anch %>% filter(analysis == lab, population == pop)

mkrow <- function(fig_row, population, n, early_n, hr, lo, hi, p, status, src_script, src_file, note = "") {
  tibble(fig_row = fig_row, population = population, n = n, early_n = early_n,
         hr = hr, ci_lower = lo, ci_upper = hi, p_value = p,
         status = status, source_script = src_script, source_file = src_file, note = note)
}
oldrow <- function(lab) old3 %>% filter(label == lab)

# --- Fig 3: 16 行 ---
p_eq  <- prim %>% filter(population == "Equipoise SOFA 8-14")
p_all <- prim %>% filter(population == "Full cohort")
s7 <- pick_strata("SOFA <=7"); s814 <- pick_strata("SOFA 8-14"); s913 <- pick_strata("SOFA 9-13"); s15 <- pick_strata("SOFA >=15")

fig3 <- bind_rows(
  mkrow("Primary SOFA 8-14 (IPCW-weighted)", "SOFA 8-14", 1502, 146, g(p_eq,"hr"), g(p_eq,"ci_lower"), g(p_eq,"ci_upper"), g(p_eq,"p_value"),
        "新口径 CH-01", "08_CH01/ch01_dod_primary.R", "08_CH01/out/ch01_primary_dod.csv", "与 RUN-15 NEW SOFA 8-14 逐位一致"),
  mkrow("Full clean-t0 cohort", "Full cohort", 2463, 195, g(p_all,"hr"), g(p_all,"ci_lower"), g(p_all,"ci_upper"), g(p_all,"p_value"),
        "新口径 CH-01", "08_CH01/ch01_dod_primary.R", "08_CH01/out/ch01_primary_dod.csv", "Schoenfeld global p=0.3394"),
  mkrow("SOFA <=7", "SOFA <=7", g(s7,"n_patients"), g(s7,"early_n"), g(s7,"weighted_hr"), g(s7,"ci_lower"), g(s7,"ci_upper"), g(s7,"p_value"),
        "新口径 CH-01", "run15_NEW_sofa_era.R", "run15/NEW/sofa_era/strictK3_v2_sofa_strata_results.csv", "本次新算"),
  mkrow("SOFA 8-14", "SOFA 8-14", g(s814,"n_patients"), g(s814,"early_n"), g(s814,"weighted_hr"), g(s814,"ci_lower"), g(s814,"ci_upper"), g(s814,"p_value"),
        "新口径 CH-01", "run15_NEW_sofa_era.R", "run15/NEW/sofa_era/strictK3_v2_sofa_strata_results.csv", "= CH-01 主分析"),
  mkrow("SOFA 9-13", "SOFA 9-13", g(s913,"n_patients"), g(s913,"early_n"), g(s913,"weighted_hr"), g(s913,"ci_lower"), g(s913,"ci_upper"), g(s913,"p_value"),
        "新口径 CH-01", "run15_NEW_sofa_era.R", "run15/NEW/sofa_era/strictK3_v2_sofa_strata_results.csv", "= CH-01 记录 1.040"),
  mkrow("SOFA >=15", "SOFA >=15", g(s15,"n_patients"), g(s15,"early_n"), g(s15,"weighted_hr"), g(s15,"ci_lower"), g(s15,"ci_upper"), g(s15,"p_value"),
        "新口径 CH-01", "run15_NEW_sofa_era.R", "run15/NEW/sofa_era/strictK3_v2_sofa_strata_results.csv", "本次新算"),
  mkrow("Stage-3 path-adjusted IPCW", "SOFA 8-14", 1502, 146, g(oldrow("Stage-3 path-adjusted IPCW"),"weighted_hr"),
        g(oldrow("Stage-3 path-adjusted IPCW"),"ci_lower"), g(oldrow("Stage-3 path-adjusted IPCW"),"ci_upper"),
        g(oldrow("Stage-3 path-adjusted IPCW"),"p_value"), "⛔ 无源（CH-05 已删）", "(none)", "out/qa_validation/strictK3_v2_fig3_forest_display.csv",
        "旧值，不得当新口径；见 无源行清单.csv"),
  mkrow("Calendar-era stratified Cox", "SOFA 8-14", 1502, 146, g(pick_anch("Calendar-era stratified Cox (NEW CH-01)","SOFA 8-14"),"HR"),
        g(pick_anch("Calendar-era stratified Cox (NEW CH-01)","SOFA 8-14"),"ci_lower"), g(pick_anch("Calendar-era stratified Cox (NEW CH-01)","SOFA 8-14"),"ci_upper"),
        g(pick_anch("Calendar-era stratified Cox (NEW CH-01)","SOFA 8-14"),"p_value"), "新口径 CH-01", "run15_NEW_anchors_L0_era.R",
        "run15/NEW/anchors_L0_era_rows.csv", "PH 全局 p=0.7174"),
  mkrow("KDIGO3 minus 6 h anchor", "SOFA 8-14", 1701, 192, g(pick_anch("legacy K3-6h anchor (NEW CH-01)","SOFA 8-14"),"HR"),
        g(pick_anch("legacy K3-6h anchor (NEW CH-01)","SOFA 8-14"),"ci_lower"), g(pick_anch("legacy K3-6h anchor (NEW CH-01)","SOFA 8-14"),"ci_upper"),
        g(pick_anch("legacy K3-6h anchor (NEW CH-01)","SOFA 8-14"),"p_value"), "新口径 CH-01", "run15_NEW_anchors_L0_era.R",
        "run15/NEW/anchors_L0_era_rows.csv", "= CH-01 记录 1.065；旧口径 1.053 已复现"),
  mkrow("KDIGO2 exploratory anchor", "SOFA 8-14", 2874, 187, g(pick_anch("KDIGO2 anchor (NEW CH-01)","SOFA 8-14"),"HR"),
        g(pick_anch("KDIGO2 anchor (NEW CH-01)","SOFA 8-14"),"ci_lower"), g(pick_anch("KDIGO2 anchor (NEW CH-01)","SOFA 8-14"),"ci_upper"),
        g(pick_anch("KDIGO2 anchor (NEW CH-01)","SOFA 8-14"),"p_value"), "新口径 CH-01", "run15_NEW_anchors_L0_era.R",
        "run15/NEW/anchors_L0_era_rows.csv", "旧口径 1.305 (1.157-1.472) 已逐位复现；本研究行由本次新算"),
  mkrow("No rescue exclusion", "SOFA 8-14", 1766, 208, g(oldrow("No rescue exclusion"),"weighted_hr"), g(oldrow("No rescue exclusion"),"ci_lower"),
        g(oldrow("No rescue exclusion"),"ci_upper"), g(oldrow("No rescue exclusion"),"p_value"), "⛔ 无源（CH-06 已删）", "(none)",
        "out/qa_validation/strictK3_v2_fig3_forest_display.csv", "旧值"),
  mkrow("Confirmed rescue thresholds", "SOFA 8-14", 1508, 148, g(oldrow("Confirmed rescue thresholds"),"weighted_hr"), g(oldrow("Confirmed rescue thresholds"),"ci_lower"),
        g(oldrow("Confirmed rescue thresholds"),"ci_upper"), g(oldrow("Confirmed rescue thresholds"),"p_value"), "⛔ 无源（CH-06 已删）", "(none)",
        "out/qa_validation/strictK3_v2_fig3_forest_display.csv", "旧值"),
  mkrow("Add R4 urine-output rescue", "SOFA 8-14", 343, 64, g(oldrow("Add R4 urine-output rescue"),"weighted_hr"), g(oldrow("Add R4 urine-output rescue"),"ci_lower"),
        g(oldrow("Add R4 urine-output rescue"),"ci_upper"), g(oldrow("Add R4 urine-output rescue"),"p_value"), "⛔ 无源（CH-06 已删）", "(none)",
        "out/qa_validation/strictK3_v2_fig3_forest_display.csv", "旧值"),
  mkrow("Liberal rescue thresholds", "SOFA 8-14", 1316, 100, g(oldrow("Liberal rescue thresholds"),"weighted_hr"), g(oldrow("Liberal rescue thresholds"),"ci_lower"),
        g(oldrow("Liberal rescue thresholds"),"ci_upper"), g(oldrow("Liberal rescue thresholds"),"p_value"), "⛔ 无源（CH-06 已删）", "(none)",
        "out/qa_validation/strictK3_v2_fig3_forest_display.csv", "旧值"),
  mkrow("Conservative rescue thresholds", "SOFA 8-14", 1631, 169, g(oldrow("Conservative rescue thresholds"),"weighted_hr"), g(oldrow("Conservative rescue thresholds"),"ci_lower"),
        g(oldrow("Conservative rescue thresholds"),"ci_upper"), g(oldrow("Conservative rescue thresholds"),"p_value"), "⛔ 无源（CH-06 已删）", "(none)",
        "out/qa_validation/strictK3_v2_fig3_forest_display.csv", "旧值"),
  mkrow("L0-only IPCW denominator", "SOFA 8-14", 1502, 146, g(pick_anch("L0-only denominator (NEW CH-01)","SOFA 8-14"),"HR"),
        g(pick_anch("L0-only denominator (NEW CH-01)","SOFA 8-14"),"ci_lower"), g(pick_anch("L0-only denominator (NEW CH-01)","SOFA 8-14"),"ci_upper"),
        g(pick_anch("L0-only denominator (NEW CH-01)","SOFA 8-14"),"p_value"), "新口径 CH-01", "run15_NEW_anchors_L0_era.R",
        "run15/NEW/anchors_L0_era_rows.csv", "= CH-01 记录 1.066；全队列 1.111 (0.996-1.240) p=.059")
)
write_csv(fig3, file.path(R15, "行值_Fig3.csv"))

# --- Fig S3: band scan 25 行（新旧并列）---
scan_cmp <- scan %>%
  select(lower, upper, analysis, n_patients, early_n, hr_new = weighted_hr, lo_new = ci_lower, hi_new = ci_upper,
         p_new = p_value, ph_new = ph_global_p) %>%
  left_join(scan_o %>% select(lower, upper, hr_old = weighted_hr, lo_old = ci_lower, hi_old = ci_upper, p_old = p_value, ph_old = ph_global_p),
            by = c("lower","upper")) %>%
  mutate(status = "新口径 CH-01", source_script = "run15_NEW_sofa_era.R",
         source_file = "run15/NEW/sofa_era/strictK3_v2_sofa_band_scan_results.csv")
write_csv(scan_cmp, file.path(R15, "行值_FigS3_bandscan.csv"))

# --- calendar-era（per-era + adjusted，含 PH）---
era_cmp <- era %>%
  select(anchor_year_group, analysis, n_patients, early_n, hr_new = weighted_hr, lo_new = ci_lower, hi_new = ci_upper, p_new = p_value, ph_new = ph_global_p) %>%
  left_join(erao %>% select(anchor_year_group, hr_old = weighted_hr, lo_old = ci_lower, hi_old = ci_upper, p_old = p_value, ph_old = ph_global_p), by = "anchor_year_group") %>%
  mutate(model = "per-era (SOFA 8-14)", status = "新口径 CH-01")
eraadj_cmp <- eraadj %>%
  select(analysis, n_patients, early_n, hr_new = weighted_hr, lo_new = ci_lower, hi_new = ci_upper, p_new = p_value, ph_new = ph_global_p) %>%
  left_join(eraadjo %>% select(analysis, hr_old = weighted_hr, lo_old = ci_lower, hi_old = ci_upper, p_old = p_value, ph_old = ph_global_p), by = "analysis") %>%
  mutate(model = "era-adjusted", status = "新口径 CH-01")
write_csv(bind_rows(era_cmp, eraadj_cmp), file.path(R15, "行值_calendar_era.csv"))

# --- SOFA strata + anchors + L0（全列 + 口径）---
strat_cmp <- strata %>%
  select(analysis, n_patients, early_n, hr_new = weighted_hr, lo_new = ci_lower, hi_new = ci_upper, p_new = p_value, ph_new = ph_global_p) %>%
  mutate(model_family = "SOFA strata", status = "新口径 CH-01", source_script = "run15_NEW_sofa_era.R")
write_csv(strat_cmp, file.path(R15, "行值_SOFA_strata.csv"))

# --- 无源行清单 ---
nosrc <- old3 %>% filter(label %in% c("Stage-3 path-adjusted IPCW","No rescue exclusion","Confirmed rescue thresholds",
                                      "Add R4 urine-output rescue","Liberal rescue thresholds","Conservative rescue thresholds")) %>%
  transmute(fig_row = label, n = n, early_n = early, hr_old = weighted_hr, ci_lower_old = ci_lower, ci_upper_old = ci_upper, p_old = p_value,
            deleted_as = ifelse(grepl("path-adjusted", label), "CH-05", "CH-06"),
            basis = "07_计划/马斯克5步法审阅_20261015.md §②：不挂 OSF 后业主消失 ⇒ 删",
            action = "不恢复跑、不删图行；等作者裁决")
write_csv(nosrc, file.path(R15, "无源行清单.csv"))

cat("Fig3 rows:", nrow(fig3), "| 新口径:", sum(fig3$status != "⛔ 无源（CH-05 已删）" & fig3$status != "⛔ 无源（CH-06 已删）"), "| 无源:", sum(grepl("无源", fig3$status)), "\n")
cat("band scan rows:", nrow(scan_cmp), "\n")
cat("era rows:", nrow(era_cmp) + nrow(eraadj_cmp), "\n")
cat("no-source list rows:", nrow(nosrc), "\n")
