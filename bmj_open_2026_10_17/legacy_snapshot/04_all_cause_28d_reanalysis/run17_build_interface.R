# RUN-17 · 合并 figure-runner 接口（9 行原有 + 6 行恢复 = 15 行）
suppressPackageStartupMessages({library(dplyr); library(readr)})
W <- Sys.getenv("CRRT_WORKSPACE_ROOT"); R15 <- file.path(W, "09_工作稿/out/run15"); R17 <- file.path(W, "09_工作稿/out/run17")
old9 <- read_csv(file.path(R15, "行值_forest_接口.csv"), show_col_types = FALSE)
pathn <- read_csv(file.path(R17, "NEW/path/results/strictK3_v2_path_adjusted_ipcw_results.csv"), show_col_types = FALSE)
rescn <- read_csv(file.path(R17, "NEW/rescue/strictK3_primary_v2/results/strictK3_v2_rescue_sensitivity_model_results.csv"), show_col_types = FALSE)
oldp  <- read_csv(file.path(W, "09_工作稿/out/qa_validation/strictK3_v2_path_adjusted_ipcw_results.csv"), show_col_types = FALSE)
oldr  <- read_csv(file.path(W, "09_工作稿/out/qa_validation/strictK3_v2_rescue_sensitivity_model_results.csv"), show_col_types = FALSE)

pr <- pathn %>% filter(analysis == "SOFA 8-14: IPCW includes stage3_path")
new6 <- bind_rows(
  tibble(group="Robustness", analysis="SOFA 8-14: IPCW includes stage3_path",
         n_patients=pr$n_patients, early_n=pr$early_n, weighted_hr=pr$HR, ci_lower=pr$CI_lo, ci_upper=pr$CI_hi, p_value=pr$p),
  rescn %>% filter(grepl(": SOFA 8-14$", analysis), !grepl("^primary_R1_R2_R2b_R3", analysis)) %>%
    transmute(group="Rescue sensitivity", analysis=analysis, n_patients=n_patients, early_n=early_n,
              weighted_hr=weighted_hr, ci_lower=ci_lower, ci_upper=ci_upper, p_value=p_value)
)
# 顺序：按 figure-runner 给的 key 顺序（path 行在前，然后 5 条 rescue 按原图序）
ord <- c("SOFA 8-14: IPCW includes stage3_path","no_rescue_exclusion: SOFA 8-14",
         "sensitivity_R1_confirm_R2_R2b_R3: SOFA 8-14","sensitivity_R1_R2_R2b_R3_R4: SOFA 8-14",
         "sensitivity_liberal_thresholds_R1_R2_R2b_R3: SOFA 8-14","sensitivity_conservative_thresholds_R1_R2_R2b_R3: SOFA 8-14")
new6 <- new6 %>% mutate(.o = match(analysis, ord)) %>% arrange(.o) %>% select(-.o)
all15 <- bind_rows(old9, new6)
write_csv(all15, file.path(R17, "行值_forest_接口_15行.csv"))

detail <- bind_rows(
  tibble(fig_row="Stage-3 path-adjusted IPCW", analysis_key=pr$analysis, population="SOFA 8-14",
         script="strictK3_path_sensitivity_v2.R", run17_script="run17_NEW_path.R",
         hr=pr$HR, ci_lo=pr$CI_lo, ci_hi=pr$CI_hi, p=pr$p, n=pr$n_patients, events=pr$events,
         old_hr=oldp$HR[oldp$analysis==pr$analysis], old_ci_lo=oldp$CI_lo[oldp$analysis==pr$analysis],
         old_ci_hi=oldp$CI_hi[oldp$analysis==pr$analysis], old_p=oldp$p[oldp$analysis==pr$analysis], dao="新口径 CH-01"),
  rescn %>% filter(grepl(": SOFA 8-14$", analysis), !grepl("^primary_R1_R2_R2b_R3", analysis)) %>%
    transmute(fig_row=recode(analysis, "no_rescue_exclusion: SOFA 8-14"="No rescue exclusion",
      "sensitivity_R1_confirm_R2_R2b_R3: SOFA 8-14"="Confirmed rescue thresholds",
      "sensitivity_R1_R2_R2b_R3_R4: SOFA 8-14"="Add R4 urine-output rescue",
      "sensitivity_liberal_thresholds_R1_R2_R2b_R3: SOFA 8-14"="Liberal rescue thresholds",
      "sensitivity_conservative_thresholds_R1_R2_R2b_R3: SOFA 8-14"="Conservative rescue thresholds"),
      analysis_key=analysis, population="SOFA 8-14", script="strictK3_rescue_sensitivity_models_v2.R",
      run17_script="run17_NEW_rescue.R", hr=weighted_hr, ci_lo=ci_lower, ci_hi=ci_upper, p=p_value,
      n=n_patients, events=events,
      old_hr=oldr$weighted_hr[match(analysis, oldr$analysis)], old_ci_lo=oldr$ci_lower[match(analysis, oldr$analysis)],
      old_ci_hi=oldr$ci_upper[match(analysis, oldr$analysis)], old_p=oldr$p_value[match(analysis, oldr$analysis)],
      dao="新口径 CH-01")
)
write_csv(detail, file.path(R17, "行值_6行_新口径.csv"))
cat("all15 rows:", nrow(all15), " | new6 rows:", nrow(new6), " | detail rows:", nrow(detail), "\n")
print(as.data.frame(all15 %>% mutate(across(where(is.numeric), ~round(.x, 4)))), row.names = FALSE)
