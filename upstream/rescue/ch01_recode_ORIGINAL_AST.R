ch01_recode <- function(df) {
  df %>%
    left_join(ch01_dod, by = "subject_id") %>%
    mutate(
      .t0 = as.POSIXct(t0, tz = "UTC"),
      .dt = as.POSIXct(deathtime, tz = "UTC"),
      .dc = coalesce(as.numeric(difftime(.dt, .t0, units = "days")), as.numeric(dod - as.Date(.t0))),
      died_28d        = !is.na(.dc) & .dc >= 0 & .dc <= 28,
      follow_time_28d = pmin(28, ifelse(is.na(.dc), 28, .dc))
    ) %>%
    select(-.dc, -.t0, -.dt, -dod)
}
