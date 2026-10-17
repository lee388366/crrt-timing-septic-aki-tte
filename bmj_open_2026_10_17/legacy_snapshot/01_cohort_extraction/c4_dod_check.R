library(DBI)
library(RPostgres)

con <- dbConnect(RPostgres::Postgres(), dbname = "mimiciv31", host = "localhost",
                 port = 5432, user = "postgres", password = Sys.getenv("PGPASSWORD"))

c4 <- dbGetQuery(con, "
WITH cohort AS (
  SELECT subject_id, hadm_id, stay_id, t0, deathtime,
         follow_time_28d, died_28d, early_crrt
  FROM data_extract_crrt.cohort_base
)
SELECT
  COUNT(*)::int                                        AS total,
  SUM((p.dod IS NOT NULL)::int)::int                   AS has_dod,
  SUM((c.deathtime IS NOT NULL)::int)::int              AS has_deathtime,
  SUM((p.dod IS NOT NULL AND c.deathtime IS NULL)::int)::int AS dod_only,
  SUM(CASE WHEN p.dod IS NOT NULL AND c.deathtime IS NULL
            AND p.dod <= (c.t0 + INTERVAL '28 day')::date
       THEN 1 ELSE 0 END)::int                         AS dod_28d_no_dt,
  SUM(CASE WHEN p.dod IS NOT NULL AND c.deathtime IS NULL
            AND p.dod <= (c.t0 + INTERVAL '28 day')::date
            AND c.died_28d = FALSE
       THEN 1 ELSE 0 END)::int                         AS missed_deaths
FROM cohort c
JOIN mimiciv_hosp.patients p USING (subject_id)
")

cat("=== C4: patients.dod Coverage ===\n")
cat(sprintf("Total patients:          %d\n", c4$total))
cat(sprintf("Has dod record:          %d (%.1f%%)\n", c4$has_dod, 100 * c4$has_dod / c4$total))
cat(sprintf("Has in-hospital death:   %d (%.1f%%)\n", c4$has_deathtime, 100 * c4$has_deathtime / c4$total))
cat(sprintf("dod only (no in-hosp):   %d\n", c4$dod_only))
cat(sprintf("dod within 28d:          %d\n", c4$dod_28d_no_dt))
cat(sprintf("MISSED deaths (current): %d (%.2f%%)\n", c4$missed_deaths, 100 * c4$missed_deaths / c4$total))

by_grp <- dbGetQuery(con, "
WITH cohort AS (
  SELECT subject_id, stay_id, t0, deathtime, died_28d, early_crrt
  FROM data_extract_crrt.cohort_base
)
SELECT
  c.early_crrt,
  COUNT(*)::int AS n,
  SUM(CASE WHEN p.dod IS NOT NULL AND c.deathtime IS NULL
            AND p.dod <= (c.t0 + INTERVAL '28 day')::date
            AND c.died_28d = FALSE
       THEN 1 ELSE 0 END)::int AS missed
FROM cohort c
JOIN mimiciv_hosp.patients p USING (subject_id)
GROUP BY c.early_crrt ORDER BY c.early_crrt
")

cat("\nBy treatment group:\n")
print(by_grp)

dbDisconnect(con)
