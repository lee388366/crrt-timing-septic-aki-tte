# Composition of the deferred group (no CRRT within 24 h of t0), by observed treatment before cloning.
# Produces the counts in Supplemental Table S16 and the corresponding Results/Discussion sentences.
# Aggregate output only; no patient-level data are written.
# Definitions as in cp_core.R::prep_pop(): died28 = death within 28 days of t0 (in-hospital death time or
# date of death); h_crrt = hours from t0 to CRRT start (NA if no CRRT before death or day 28).
# Usage: CRRT_REVISION_DIR=<folder containing private_data/> Rscript counts_deferred_2026-11-03.R
suppressPackageStartupMessages(library(dplyr))
PKG <- Sys.getenv("CRRT_REVISION_DIR", ".."); PRIV <- file.path(PKG, "private_data")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds"))
for (pop in c("primary", "full")) {
  p <- if (pop == "primary") filter(coh, eq_strict) else coh       # primary: SOFA at t0 8-14
  d <- filter(p, !early)
  never <- is.na(d$h_crrt); late <- !never                          # late = CRRT started at or after 24 h
  inh <- !is.na(d$d_in) & d$d_in >= 0 & d$d_in <= 672               # death time from hospital record (hourly precision)
  dn <- d$died28 & never
  cat(sprintf("\n[%s] deferred n = %d | never CRRT within 28 d = %d | CRRT started >= 24 h = %d\n", pop, nrow(d), sum(never), sum(late)))
  cat(sprintf("  survived 28 d and never CRRT:            %d\n", sum(!d$died28 & never)))
  cat(sprintf("  died within 28 d without CRRT:           %d\n", sum(dn)))
  cat(sprintf("    of whom died < 24 h after t0:          %d  (hospital death time %d; date-only deaths %d)\n",
      sum(dn & d$death_h < 24), sum(dn & inh & d$death_h < 24), sum(dn & !inh & d$death_h < 24)))
  cat(sprintf("    died 24 h - day 28:                    %d\n", sum(dn & d$death_h >= 24)))
  cat(sprintf("  CRRT started >= 24 h: died within 28 d:  %d / %d (%.1f%%)\n", sum(d$died28 & late), sum(late), 100 * mean(d$died28[late])))
  cat(sprintf("    median start of late CRRT, h [IQR]:    %.1f [%.1f, %.1f]\n", median(d$h_crrt[late]), quantile(d$h_crrt[late], .25), quantile(d$h_crrt[late], .75)))
}
