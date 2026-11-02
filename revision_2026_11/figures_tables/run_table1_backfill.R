# Table 1 variables not produced by the generator (weight, arterial pH, invasive ventilation); JOIN logic verbatim from upstream/table1/table1_9missing_v2_extract.sql
suppressPackageStartupMessages({library(dplyr); library(DBI); library(RPostgres)})
ROOT <- Sys.getenv("CRRT_REVISION_DIR")
coh <- readRDS(file.path(Sys.getenv("CRRT_REVISION_DIR"), "private_data/cohort_round1.rds"))
con <- dbConnect(RPostgres::Postgres(), dbname = "mimiciv31", host = "localhost", port = 5432, user = "postgres", password = Sys.getenv("PGPASSWORD"))
out <- list()
for (pop in c("full", "primary")) {
  p <- if (pop == "primary") filter(coh, eq_strict) else coh
  vals <- paste(sprintf("(%d,%d,'%s'::timestamp,%s)", p$stay_id, p$hadm_id, format(p$t0, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                        ifelse(p$early, "TRUE", "FALSE")), collapse = ",")
  sql <- paste0("
WITH cohort AS (SELECT stay_id, hadm_id, t0 AS t0_time, t0 - INTERVAL '6 hour' AS baseline_start, early_crrt
                FROM (VALUES ", vals, ") AS vv(stay_id,hadm_id,t0,early_crrt)),
pat AS (
  SELECT c.stay_id, c.early_crrt, w.weight_kg, bg.ph, sp.mv_baseline
  FROM cohort c
  LEFT JOIN LATERAL (SELECT wd.weight AS weight_kg FROM mimiciv_derived.weight_durations wd
    WHERE wd.stay_id=c.stay_id AND wd.starttime<=c.t0_time AND wd.weight IS NOT NULL
    ORDER BY wd.starttime DESC LIMIT 1) w ON true
  LEFT JOIN LATERAL (SELECT bb.ph FROM mimiciv_derived.bg bb
    WHERE bb.hadm_id=c.hadm_id AND bb.charttime>c.baseline_start AND bb.charttime<=c.t0_time
    ORDER BY bb.charttime DESC LIMIT 1) bg ON true
  LEFT JOIN LATERAL (SELECT EXISTS(SELECT 1 FROM mimiciv_derived.ventilation vv
    WHERE vv.stay_id=c.stay_id AND vv.starttime<c.t0_time AND COALESCE(vv.endtime,c.t0_time)>c.baseline_start
      AND vv.ventilation_status='InvasiveVent')::int AS mv_baseline FROM (SELECT 1) _) sp ON true)
SELECT (CASE WHEN early_crrt THEN 'early' WHEN NOT early_crrt THEN 'deferred' ELSE 'overall' END) AS arm, count(*) n,
  count(weight_kg) wt_n, percentile_cont(ARRAY[0.25,0.5,0.75]) within group(order by weight_kg) wt_q,
  count(ph) ph_n, percentile_cont(ARRAY[0.25,0.5,0.75]) within group(order by ph) ph_q,
  sum(mv_baseline) mv_n
FROM pat GROUP BY GROUPING SETS ((early_crrt),())")
  r <- dbGetQuery(con, sql); r$pop <- pop; out[[pop]] <- r
}
dbDisconnect(con)
res <- bind_rows(out); res$wt_q <- sapply(res$wt_q, function(x) paste(round(as.numeric(gsub("[{}]", "", strsplit(x, ",")[[1]])), 2), collapse = " / "))
res$ph_q <- sapply(res$ph_q, function(x) paste(round(as.numeric(gsub("[{}]", "", strsplit(x, ",")[[1]])), 3), collapse = " / "))
res$mv_pct <- round(100 * as.numeric(res$mv_n) / as.numeric(res$n), 1)
print(res, row.names = FALSE)
write.csv(res, file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables/table1_backfill_weight_ph_mv.csv"), row.names = FALSE)
# For Table 2: ICU length of stay from t0 (days), by observed initiation within 24 h (descriptive)
q <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
for (pop in c("primary", "full")) {
  p <- if (pop == "primary") filter(coh, eq_strict) else coh; los <- p$icu_los_from_t0
  if (max(los, na.rm = TRUE) > 200) los <- los / 24      # convert hours to days if needed
  cat(sprintf("ICU LOS from t0 (days) [%s]: early %s | deferred %s | Wilcoxon P = %.3f\n", pop, q(los[p$early]), q(los[!p$early]),
      wilcox.test(los[p$early], los[!p$early])$p.value))
}
