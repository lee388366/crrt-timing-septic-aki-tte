# Round 6: crude HR (unweighted clones) and global Schoenfeld test of the weighted model
setwd(Sys.getenv("CRRT_SCRIPT_DIR", getwd()))
source("cp_core.R"); L0 <- setdiff(L0, "apsiii")
coh <- readRDS("../private_data/cohort_round1.rds"); tvl <- readRDS("../private_data/tv_hourly_round1.rds")
for (sub in c("primary", "full")) {
  p <- if (sub == "primary") filter(coh, eq_strict) else coh; p$id <- seq_len(nrow(p))
  cp <- build(p, tvl, wlag = TRUE)
  u <- est(cp %>% mutate(w = 1)); w <- est(cp)
  m <- coxph(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w, ties = "efron")
  cat(sprintf("[%s] crude HR (unweighted) %.3f (%.3f-%.3f), P = %.3f | weighted HR %.3f | global Schoenfeld P = %.2f\n",
      sub, u[["HR"]], u[["lo"]], u[["hi"]], u[["p"]], w[["HR"]], cox.zph(m)$table["GLOBAL", "p"]))
}
