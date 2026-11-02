# Round 5: supplementary analyses and counts reported in the manuscript
#  1) flow-diagram and text counts  2) mortality ascertainment (in-hospital death, censoring at discharge)  3) acid-base/electrolyte adjustment
source("cp_core.R")
suppressPackageStartupMessages({library(DBI); library(RPostgres)})
PKG <- Sys.getenv("CRRT_REVISION_DIR")
OUT <- file.path(PKG, "out"); PRIV <- file.path(PKG, "private_data")
sink(file.path(OUT, "round5_report.txt"), split = TRUE); cat("Round5 run at", format(Sys.time()), "\n\n")
L0 <- setdiff(L0, "apsiii")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds")); tvl <- readRDS(file.path(PRIV, "tv_hourly_round1.rds"))
raw <- readRDS(file.path(PRIV, "cohort_v3_step1.rds"))          # before removal of inconsistent records (2,463)
fun <- read.csv(Sys.getenv("CRRT_FUNNEL_CSV"))   # strictK3_primary_v2_funnel.csv written by upstream/analysis/strictK3_primary_extract_v2.R

# ---------- 1. Counts ----------
cat("== 1. Flow diagram\n"); print(fun[, c("step", "n", "excluded_from_previous")], row.names = FALSE)
cat(sprintf("After urgent-indication exclusion %d -> inconsistent timestamps removed %d (CRRT after death / t0 before ICU admission) -> analysis cohort %d (early %d, deferred %d)\n",
  nrow(raw), nrow(raw) - nrow(coh), nrow(coh), sum(coh$early), sum(!coh$early)))
sb <- cut(coh$sofa_t0_strict, c(-Inf, 7, 14, Inf), labels = c("SOFA ≤7", "SOFA 8–14", "SOFA ≥15"))
print(table(SOFA = sb, early = coh$early))
for (sub in c("primary", "full")) {
  p <- if (sub == "primary") filter(coh, eq_strict) else coh
  dfr <- p[!p$early, ]; never <- is.na(dfr$crrt_first_time) | is.na(dfr$h_crrt)
  inh <- !is.na(p$d_in) & p$d_in >= 0 & p$d_in <= 672; dodonly <- p$died28 & !inh
  cat(sprintf("\n[%s] n=%d | early %d, observed deaths %d/%d (%.1f%%) | deferred %d, observed deaths %d/%d (%.1f%%)\n", sub, nrow(p), sum(p$early),
    sum(p$died28 & p$early), sum(p$early), 100 * mean(p$died28[p$early]), sum(!p$early), sum(p$died28 & !p$early), sum(!p$early), 100 * mean(p$died28[!p$early])))
  cat(sprintf("  Deferred group never receiving CRRT within 28 days: %d/%d (%.1f%%); initiated after 24 h: %d\n", sum(never), nrow(dfr), 100 * mean(never), sum(!never)))
  cat(sprintf("  28-day deaths identified only by date of death: %d (early %d, deferred %d); date of death without in-hospital death time: %d\n", sum(dodonly),
    sum(dodonly & p$early), sum(dodonly & !p$early), sum(!is.na(p$dod) & is.na(p$deathtime))))
}

# ---------- 2. Mortality ascertainment ----------
cat("\n== 2. Mortality ascertainment: in-hospital death only, censored at discharge\n")
for (sub in c("primary", "full")) {
  p <- if (sub == "primary") filter(coh, eq_strict) else coh
  p <- p %>% mutate(disch_h = as.numeric(difftime(dischtime, t0, units = "hours")),
    inh = !is.na(d_in) & d_in >= 0 & d_in <= 672,
    died28 = inh, end_h = case_when(inh ~ pmax(d_in, 0.5), !is.na(disch_h) & disch_h < 672 ~ pmax(disch_h, 0.5), TRUE ~ 672),
    h_crrt = if_else(!is.na(h_crrt) & h_crrt < end_h, h_crrt, NA_real_), early = !is.na(h_crrt) & h_crrt < 24)
  p$id <- seq_len(nrow(p)); e <- est(build(p, tvl, wlag = TRUE))
  cat(sprintf("  %-8s n=%d early=%d events=%d | HR %.3f (%.3f-%.3f) | RMST diff %.2f\n", sub, nrow(p), sum(p$early), sum(p$died28), e["HR"], e["lo"], e["hi"], e["rmst_diff"]))
}

# ---------- 3. Acid-base and electrolyte adjustment ----------
cat("\n== 3. Acid-base/electrolytes: latest pH, potassium and bicarbonate within 24 h before t0 added to L0\n")
con <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("PGDATABASE", "mimiciv31"), host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")), user = Sys.getenv("PGUSER", "postgres"), password = Sys.getenv("PGPASSWORD"))
H <- paste(unique(coh$hadm_id), collapse = ",")
ph <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, ph AS v FROM mimiciv_derived.bg WHERE ph IS NOT NULL AND hadm_id IN (", H, ")"))
ch <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, potassium, bicarbonate FROM mimiciv_derived.chemistry WHERE hadm_id IN (", H, ")"))
dbDisconnect(con)
lastv <- function(r, nm) r %>% filter(!is.na(v)) %>% inner_join(coh %>% select(stay_id, hadm_id, t0), by = "hadm_id", relationship = "many-to-many") %>%
  filter(charttime <= t0, charttime >= t0 - hours(24)) %>% group_by(stay_id) %>% slice_max(charttime, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(stay_id, v) %>% rename(!!nm := v)
ab <- coh %>% left_join(lastv(ph, "ph_t0"), by = "stay_id") %>% left_join(lastv(ch %>% transmute(hadm_id, charttime, v = potassium), "k_t0"), by = "stay_id") %>%
  left_join(lastv(ch %>% transmute(hadm_id, charttime, v = bicarbonate), "hco3_t0"), by = "stay_id")
cat(sprintf("  Missing: pH %.1f%% | potassium %.1f%% | bicarbonate %.1f%%\n", 100 * mean(is.na(ab$ph_t0)), 100 * mean(is.na(ab$k_t0)), 100 * mean(is.na(ab$hco3_t0))))
L0b <- L0; L0 <- c(L0b, "ph_t0", "k_t0", "hco3_t0")
for (sub in c("primary", "full")) {
  p <- if (sub == "primary") filter(ab, eq_strict) else ab; p$id <- seq_len(nrow(p))
  cp <- build(p, tvl, wlag = TRUE); e <- est(cp); d <- attr(cp, "diag")
  cat(sprintf("  %-8s HR %.3f (%.3f-%.3f) | RMST diff %.2f | ESS_A %.1f/%d | max|SMD| %.3f\n", sub, e["HR"], e["lo"], e["hi"], e["rmst_diff"], d["ESS_A"], d["n_A_post"], d["max_absSMD_A"]))
}
L0 <- L0b
cat("\nDONE\n"); sink()
