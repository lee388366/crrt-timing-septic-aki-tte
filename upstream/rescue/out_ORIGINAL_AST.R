out <- bind_rows(results) %>%
  mutate(
    weighted_hr = round(weighted_hr, 3),
    ci_lower = round(ci_lower, 3),
    ci_upper = round(ci_upper, 3),
    p_value = round(p_value, 4)
  )
