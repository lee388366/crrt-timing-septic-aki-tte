variant_cohorts <- list(
  no_rescue_exclusion = pre_rescue_ext,
  primary_R1_R2_R2b_R3 = primary_cohort,
  sensitivity_R1_R2_R2b_R3_R4 = pre_rescue_ext %>% filter(!rescue_r1_r4_sensitivity),
  sensitivity_R1_confirm_R2_R2b_R3 = pre_rescue_ext %>%
    mutate(r1_confirm_trigger = r1_trigger & !r1_single_unrechecked_k_high) %>%
    filter(!(r1_confirm_trigger | r2_trigger | r2b_trigger | r3_trigger)),
  sensitivity_liberal_thresholds_R1_R2_R2b_R3 = pre_rescue_ext %>%
    filter(!rescue_threshold_liberal),
  sensitivity_conservative_thresholds_R1_R2_R2b_R3 = pre_rescue_ext %>%
    filter(!rescue_threshold_conservative)
)
