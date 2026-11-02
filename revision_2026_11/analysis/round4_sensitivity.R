# Round 4: all sensitivity analyses under the final specification (run with options(cp_trunc = TRUE) via round7 for the primary specification)
# Final specification = strict SOFA(t0) + baseline values at or before t0 + strictly lagged hourly SOFA/KDIGO + weight at interval start + L0 without APS III + early-arm L0 landmark
# Key contrasts (RRT-inclusive staging; creatinine-only vs urine-output-only) use 1,000 bootstrap resamples + grouped jackknife BCa; others use robust CIs
# Usage: Rscript round4_sensitivity.R   (bootstrap is resumable)
source("cp_core.R")
suppressPackageStartupMessages({library(parallel); library(DBI); library(RPostgres); library(zoo)})
PKG <- Sys.getenv("CRRT_REVISION_DIR")
WS  <- Sys.getenv("CRRT_INPUT_DIR")
OUT <- file.path(PKG, "out"); PRIV <- file.path(PKG, "private_data"); BOOT <- file.path(PRIV, "boot")
sink(file.path(OUT, "round4_report.txt"), split = TRUE); cat("Round4 run at", format(Sys.time()), "\n\n")
NREP <- 1000; NJK <- 200; CORES <- 4; RNGkind("L'Ecuyer-CMRG")
L0_full <- L0; L0 <- setdiff(L0_full, "apsiii")

ow <- function(b) { for (d in file.path(WS, "ch24_q4_corrected_other_weights", b, c("NEW/data", "NEW"))) if (file.exists(file.path(d, "baseline.rds"))) return(d); NA }
branches <- list(
  list(key = "primary", name = "Primary (physiologic KDIGO 3)", dir = file.path(WS, "ch24_q4_corrected_inputs"), kdigo0 = 3),
  list(key = "rrt", name = "RRT-inclusive staging", dir = ow("RRT_inclusive_v1"), kdigo0 = 3),
  list(key = "m6h", name = "Anchor: 6 h before KDIGO 3", dir = ow("anchor_minus6h"), kdigo0 = NA),
  list(key = "st2", name = "Anchor: KDIGO stage 2", dir = ow("anchor_stage2"), kdigo0 = 2),
  list(key = "nores", name = "Urgent exclusion: none", dir = ow("rescue_no_rescue_exclusion"), kdigo0 = 3),
  list(key = "r1c", name = "Urgent exclusion: hyperkalaemia reconfirmed", dir = ow("rescue_sensitivity_R1_confirm_R2_R2b_R3"), kdigo0 = 3),
  list(key = "cons", name = "Urgent exclusion: cut-offs raised", dir = ow("rescue_sensitivity_conservative_thresholds_R1_R2_R2b_R3"), kdigo0 = 3),
  list(key = "lib", name = "Urgent exclusion: cut-offs lowered", dir = ow("rescue_sensitivity_liberal_thresholds_R1_R2_R2b_R3"), kdigo0 = 3))
branches <- branches[!is.na(sapply(branches, `[[`, "dir"))]
bl  <- lapply(branches, function(b) readRDS(file.path(b$dir, "baseline.rds")))
tvs <- lapply(branches, function(b) readRDS(file.path(b$dir, "tv_hourly.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE))
tv48 <- readRDS(file.path(WS, "ch24_q4_corrected_tv48", "tv_hourly_0_47.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE)

# ---------- One-off queries (all branches) ----------
allb <- bind_rows(lapply(bl, function(x) x %>% select(stay_id, hadm_id, subject_id, t0))) %>% distinct()
S <- paste(unique(allb$stay_id), collapse = ","); H <- paste(unique(allb$hadm_id), collapse = ","); P <- paste(unique(allb$subject_id), collapse = ",")
con <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("PGDATABASE", "mimiciv31"), host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")), user = Sys.getenv("PGUSER", "postgres"), password = Sys.getenv("PGPASSWORD"))
q <- function(sql) dbGetQuery(con, sql)
sofa_raw  <- q(paste0("SELECT stay_id, starttime, endtime, sofa_24hours AS s FROM mimiciv_derived.sofa WHERE sofa_24hours IS NOT NULL AND stay_id IN (", S, ")"))
kdigo_raw <- q(paste0("SELECT stay_id, charttime, aki_stage FROM mimiciv_derived.kdigo_stages WHERE aki_stage IS NOT NULL AND stay_id IN (", S, ")"))
creat_raw <- q(paste0("SELECT hadm_id, charttime, creatinine FROM mimiciv_derived.chemistry WHERE creatinine IS NOT NULL AND creatinine < 25 AND hadm_id IN (", H, ")"))
lact_raw  <- q(paste0("SELECT hadm_id, charttime, lactate FROM mimiciv_derived.bg WHERE lactate IS NOT NULL AND hadm_id IN (", H, ")"))
uo_raw    <- q(paste0("SELECT stay_id, charttime, uo_mlkghr_6hr AS uo FROM mimiciv_derived.urine_output_rate WHERE uo_mlkghr_6hr IS NOT NULL AND stay_id IN (", S, ")"))
map_raw   <- q(paste0("SELECT stay_id, charttime, mbp AS map FROM mimiciv_derived.vitalsign WHERE mbp IS NOT NULL AND stay_id IN (", S, ")"))
dod <- q(paste0("SELECT subject_id, dod FROM mimiciv_hosp.patients WHERE subject_id IN (", P, ")")) %>% mutate(dod = as.Date(dod))
dbDisconnect(con)

# ---------- Common timing corrections ----------
fix_pop <- function(b) {
  ss <- sofa_raw %>% inner_join(b %>% select(stay_id, t0) %>% distinct(), by = "stay_id", relationship = "many-to-many") %>%
    filter(endtime <= t0, endtime >= t0 - hours(48)) %>% group_by(stay_id) %>% slice_max(starttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(stay_id, sofa_t0_strict = s)
  p <- prep_pop(b, ss, dod); nb <- attr(p, "n_bad")
  lb <- function(raw, key, val, lower_h) raw %>% inner_join(p %>% select(stay_id, hadm_id, t0), by = key, relationship = "many-to-many") %>%
    filter(charttime <= t0, charttime >= t0 - hours(lower_h)) %>% group_by(stay_id) %>% slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(stay_id, v = all_of(val))
  p <- p %>% select(-any_of(c("creatinine_t0", "lactate_t0", "uo_rate_t0", "map_t0"))) %>%
    left_join(lb(creat_raw, "hadm_id", "creatinine", 1e5) %>% rename(creatinine_t0 = v), by = "stay_id") %>%
    left_join(lb(lact_raw, "hadm_id", "lactate", 12) %>% rename(lactate_t0 = v), by = "stay_id") %>%
    left_join(lb(uo_raw, "stay_id", "uo", 1e5) %>% rename(uo_rate_t0 = v), by = "stay_id") %>%
    left_join(lb(map_raw, "stay_id", "map", 2) %>% rename(map_t0 = v), by = "stay_id")
  p$id <- seq_len(nrow(p)); attr(p, "n_bad") <- nb; p
}
fix_tv <- function(p, tv, G = 24) {
  grid <- p %>% select(stay_id, hadm_id, t0, creatinine_t0) %>% crossing(hour_index = 0:(G - 1)) %>% mutate(we = t0 + hours(hour_index + 1))
  sv <- grid %>% select(stay_id, hour_index, t0, we) %>% inner_join(sofa_raw, by = "stay_id", relationship = "many-to-many") %>%
    filter(endtime <= we, endtime >= t0 - hours(48)) %>% group_by(stay_id, hour_index) %>% slice_max(starttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(stay_id, hour_index, sofa_h = s)
  kv <- grid %>% select(stay_id, hour_index, t0, we) %>% inner_join(kdigo_raw, by = "stay_id", relationship = "many-to-many") %>%
    filter(charttime < we, charttime >= t0 - hours(72)) %>% group_by(stay_id, hour_index) %>% slice_max(charttime, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(stay_id, hour_index, kdigo_h = aki_stage)
  cv <- grid %>% select(stay_id, hadm_id, hour_index, we, creatinine_t0) %>% inner_join(creat_raw, by = "hadm_id", relationship = "many-to-many") %>%
    filter(charttime >= we - hours(6), charttime < we) %>% group_by(stay_id, hour_index, creatinine_t0) %>%
    summarise(cr = mean(creatinine), .groups = "drop") %>% mutate(creatinine_delta = cr - coalesce(creatinine_t0, cr)) %>%
    select(stay_id, hour_index, creatinine_delta)
  tv %>% filter(stay_id %in% p$stay_id, hour_index < G) %>% select(-any_of(c("sofa_h", "kdigo_h", "creatinine_delta"))) %>%
    left_join(sv, by = c("stay_id", "hour_index")) %>% left_join(kv, by = c("stay_id", "hour_index")) %>% left_join(cv, by = c("stay_id", "hour_index")) %>%
    group_by(stay_id) %>% arrange(hour_index, .by_group = TRUE) %>%
    mutate(across(c(sofa_h, kdigo_h, creatinine_delta), ~ zoo::na.locf(., na.rm = FALSE))) %>% ungroup()
}

# ---------- Estimation and bootstrap ----------
rows <- list()
run <- function(lab, pop, tvl, kdigo0 = 3, G = 24, era = FALSE, use_tv = TRUE, adjL0 = FALSE, amode = "L0", trunc = FALSE) {
  r <- tryCatch({
    pop$id <- seq_len(nrow(pop)); cp <- build(pop, tvl, G = G, kdigo0 = kdigo0, era = era, use_tv = use_tv, wlag = TRUE, amode = amode)
    if (trunc) { cA <- quantile(cp$w[cp$strategy == "Early" & cp$tstart > 0], .99); cB <- quantile(cp$w[cp$strategy == "Deferred"], .99)
      cp <- cp %>% mutate(w = if_else(strategy == "Early", if_else(tstart > 0, pmin(w, cA), w), pmin(w, cB))) }
    e <- est(cp, era = era, adjL0 = adjL0, pop = pop); d <- attr(cp, "diag")
    data.frame(analysis = lab, n = nrow(pop), early = sum(pop$early), deaths = sum(pop$died28), t(e), t(d))
  }, error = function(e) { cat("  !! failed:", lab, conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) { cat(sprintf("%-44s n=%4d early=%3d | HR %.3f (%.3f-%.3f) | RMST diff %5.2f | ESS_A %.0f/%d | max|SMD| %.2f | maxw A %.1f B %.1f\n",
    lab, r$n, r$early, r$HR, r$lo, r$hi, r$rmst_diff, r$ESS_A, r$n_A_post, r$max_absSMD_A, r$max_wA, r$max_wB)); rows[[length(rows) + 1]] <<- r }
}
statH <- function(pop, tvl) {
  cp <- build(pop, tvl, wlag = TRUE); m <- coxph(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w, ties = "efron")
  sf <- survfit(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w)
  rm <- summary(sf, rmean = 28)$table[, "rmean"]; s28 <- summary(sf, times = 28, extend = TRUE)$surv
  c(logHR = unname(coef(m)), rmst_diff = unname(rm[2] - rm[1]), rd = s28[1] - s28[2])
}
reid <- function(p) { p$id <- seq_len(nrow(p)); p }
bca <- function(th, bt, jk, alpha = c(.025, .975)) { bt <- bt[is.finite(bt)]; jk <- jk[is.finite(jk)]
  z0 <- qnorm(mean(bt < th) + 0.5 * mean(bt == th)); jm <- mean(jk); a <- sum((jm - jk)^3) / (6 * sum((jm - jk)^2)^1.5)
  z <- qnorm(alpha); quantile(bt, pnorm(z0 + (z0 + z) / (1 - a * (z0 + z))), names = FALSE) }
run_chunks <- function(f, idx, fun) { done <- if (file.exists(f)) readRDS(f) else NULL
  todo <- setdiff(idx, if (is.null(done)) integer(0) else done[, "rep"])
  for (ch in split(todo, ceiling(seq_along(todo) / 100))) {
    r <- do.call(rbind, mclapply(ch, function(b) c(rep = b, tryCatch(fun(b), error = function(e) c(logHR = NA, rmst_diff = NA, rd = NA))), mc.cores = CORES))
    done <- rbind(done, r); saveRDS(done, f); cat("   ", format(Sys.time(), "%H:%M:%S"), basename(f), nrow(done), "/", length(idx), "\n") }
  done }
bootBCa <- function(lab, pop, tvl, key) {
  pop <- reid(pop); pe <- statH(pop, tvl)
  bt <- run_chunks(file.path(BOOT, paste0("round4_boot_", key, ".rds")), 1:NREP, function(b) { set.seed(20261102 + b); statH(reid(pop[sample.int(nrow(pop), replace = TRUE), ]), tvl) })
  grp <- { set.seed(7); sample(rep(seq_len(NJK), length.out = nrow(pop))) }
  jk <- run_chunks(file.path(BOOT, paste0("round4_jk_", key, ".rds")), 1:NJK, function(g) statH(reid(pop[grp != g, ]), tvl))
  ci <- bca(pe[["logHR"]], bt[, "logHR"], jk[, "logHR"]); cr <- bca(pe[["rmst_diff"]], bt[, "rmst_diff"], jk[, "rmst_diff"]); cd <- bca(pe[["rd"]], bt[, "rd"], jk[, "rd"])
  cat(sprintf("  ★ %-36s HR %.3f (BCa %.3f–%.3f) | risk difference %.1f (%.1f to %.1f) points | RMST diff %.2f (%.2f to %.2f)\n", lab, exp(pe[["logHR"]]), exp(ci[1]), exp(ci[2]),
      100 * pe[["rd"]], 100 * cd[1], 100 * cd[2], pe[["rmst_diff"]], cr[1], cr[2]))
  rows[[length(rows) + 1]] <<- data.frame(analysis = paste0(lab, " [BCa]"), HR = exp(pe[["logHR"]]), lo = exp(ci[1]), hi = exp(ci[2]), rmst_diff = pe[["rmst_diff"]], rd = pe[["rd"]])
  c(logHR = pe[["logHR"]], lo = ci[1], hi = ci[2])
}

# ---------- Main loop ----------
for (i in seq_along(branches)) {
  br <- branches[[i]]; pop <- fix_pop(bl[[i]]); tvf <- fix_tv(pop, tvs[[i]])
  cat("\n== ", br$name, " (inconsistent records removed ", attr(pop, "n_bad"), "; missing hourly SOFA/KDIGO ", sum(is.na(tvf$sofa_h)), "/", sum(is.na(tvf$kdigo_h)), ")\n", sep = "")
  run(paste0(br$name, " | SOFA 8–14"), filter(pop, eq_strict), tvf, br$kdigo0)
  run(paste0(br$name, " | full cohort"), pop, tvf, br$kdigo0)
  if (br$key %in% c("m6h", "st2", "nores")) {
    run(paste0(br$name, " | SOFA 8-14 | 99% truncation"), filter(pop, eq_strict), tvf, br$kdigo0, trunc = TRUE)
    run(paste0(br$name, " | full cohort | 99% truncation"), pop, tvf, br$kdigo0, trunc = TRUE)
  }
  if (br$key == "rrt") {
    e <- pop[pop$early & pop$eq_strict, ]
    cat(sprintf("  RRT-inclusive SOFA 8-14 early initiators %d: stage3_by_crrt=TRUE %d (%.1f%%); not in physiologic cohort %d (%.1f%%)\n", nrow(e),
      sum(e$stage3_by_crrt %in% TRUE), 100 * mean(e$stage3_by_crrt %in% TRUE), sum(!(e$stay_id %in% P0$stay_id)), 100 * mean(!(e$stay_id %in% P0$stay_id))))
    bootBCa("RRT-inclusive | SOFA 8–14", filter(pop, eq_strict), tvf, "rrt_sofa")
  }
  if (br$key == "primary") {
    P0 <- filter(pop, eq_strict); TV0 <- tvf
    cat("\n== Other sensitivity analyses in the primary population\n")
    run("Deferred-arm weights, baseline only", P0, TV0, use_tv = FALSE)
    run("Calendar era", P0, TV0, era = TRUE)
    run("Outcome model baseline-adjusted", P0, TV0, adjL0 = TRUE)
    run("Early-arm weights, hourly initiation model", P0, TV0, amode = "haz")
    run("Weights truncated at 99th percentile", P0, TV0, trunc = TRUE)
    L0 <- L0_full; run("APS III retained", P0, TV0); L0 <- setdiff(L0_full, "apsiii")
    run("grace 12h", P0, TV0, G = 12)
    p48 <- fix_pop(bl[[i]]) %>% mutate(early = !is.na(h_crrt) & h_crrt < 48) %>% filter(eq_strict)
    run("grace 48h", p48, fix_tv(p48, tv48, G = 48), G = 48)
    cat("\n== KDIGO stage 3 criterion (post hoc, full cohort)\n")
    pc <- filter(pop, stage3_path == "creatinine_only"); pu <- filter(pop, stage3_path == "urine_output_only")
    run("Creatinine criterion only", pc, tvf); run("Urine-output criterion only", pu, tvf)
    bc <- bootBCa("Creatinine criterion only", pc, tvf, "creat_only"); bu <- bootBCa("Urine-output criterion only", pu, tvf, "uo_only")
    xc <- readRDS(file.path(BOOT, "round4_boot_creat_only.rds"))[, "logHR"]; xu <- readRDS(file.path(BOOT, "round4_boot_uo_only.rds"))[, "logHR"]
    dlt <- xc - xu; dlt <- dlt[is.finite(dlt)]; d0 <- bc[["logHR"]] - bu[["logHR"]]
    cat(sprintf("  Interaction (percentiles of independent bootstrap differences): ratio of HRs %.3f (%.3f-%.3f), two-sided P approx %.3f\n", exp(d0), exp(quantile(dlt, .025)), exp(quantile(dlt, .975)),
        2 * min(mean(dlt <= 0), mean(dlt >= 0))))
  }
}
L0 <- L0_full
write.csv(bind_rows(rows), file.path(OUT, "round4_estimates.csv"), row.names = FALSE)
cat("\nDONE\n"); sink()
