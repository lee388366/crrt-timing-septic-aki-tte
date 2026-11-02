# Figure 3 row data (final specification): recomputes each row and writes tables/fig3_forest_rows.csv
# Primary, full cohort, RRT-inclusive and untruncated rows: BCa CIs from saved bootstrap/jackknife replicates, bootstrap P values; other rows: robust-variance CIs and Wald P
# Timing corrections and branch handling are taken verbatim from round4_sensitivity.R (extracted between section markers)
setwd(Sys.getenv("CRRT_SCRIPT_DIR", getwd()))
source("cp_core.R"); suppressPackageStartupMessages({library(DBI); library(RPostgres); library(zoo)})
PKG <- Sys.getenv("CRRT_REVISION_DIR"); WS <- Sys.getenv("CRRT_INPUT_DIR")
PRIV <- file.path(PKG, "private_data"); BOOT <- file.path(PRIV, "boot")
OUTCSV <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables/fig3_forest_rows.csv")
L0_full <- L0; L0 <- setdiff(L0_full, "apsiii")
ow <- function(b) { for (d in file.path(WS, "ch24_q4_corrected_other_weights", b, c("NEW/data", "NEW"))) if (file.exists(file.path(d, "baseline.rds"))) return(d); NA }
branches <- list(
  list(key = "primary", dir = file.path(WS, "ch24_q4_corrected_inputs")), list(key = "rrt", dir = ow("RRT_inclusive_v1")),
  list(key = "r1c", dir = ow("rescue_sensitivity_R1_confirm_R2_R2b_R3")), list(key = "r14", dir = ow("rescue_sensitivity_R1_R2_R2b_R3_R4")),
  list(key = "lib", dir = ow("rescue_sensitivity_liberal_thresholds_R1_R2_R2b_R3")), list(key = "cons", dir = ow("rescue_sensitivity_conservative_thresholds_R1_R2_R2b_R3")))
stopifnot(!any(is.na(sapply(branches, `[[`, "dir"))))
bl  <- lapply(branches, function(b) readRDS(file.path(b$dir, "baseline.rds"))); names(bl) <- sapply(branches, `[[`, "key")
tvs <- lapply(branches, function(b) readRDS(file.path(b$dir, "tv_hourly.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE)); names(tvs) <- names(bl)
tv48 <- readRDS(file.path(WS, "ch24_q4_corrected_tv48", "tv_hourly_0_47.rds")) %>% distinct(stay_id, hour_index, .keep_all = TRUE)
r4 <- readLines("round4_sensitivity.R", warn = FALSE)
seg <- function(a, b) r4[(grep(a, r4, fixed = TRUE)[1]):(grep(b, r4, fixed = TRUE)[1] - 1)]
eval(parse(text = seg("# ---------- One-off queries", "# ---------- Estimation and bootstrap")))   # queries + fix_pop + fix_tv
con <- dbConnect(RPostgres::Postgres(), dbname = Sys.getenv("PGDATABASE", "mimiciv31"), host = Sys.getenv("PGHOST", "localhost"),
  port = as.integer(Sys.getenv("PGPORT", "5432")), user = Sys.getenv("PGUSER", "postgres"), password = Sys.getenv("PGPASSWORD"))
P <- lapply(bl, fix_pop); TV <- mapply(function(p, t) fix_tv(p, t), P, tvs, SIMPLIFY = FALSE)
P0 <- filter(P$primary, eq_strict); TV0 <- TV$primary
H <- paste(unique(P0$hadm_id), collapse = ",")
ph <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, ph AS v FROM mimiciv_derived.bg WHERE ph IS NOT NULL AND hadm_id IN (", H, ")"))
ch <- dbGetQuery(con, paste0("SELECT hadm_id, charttime, potassium, bicarbonate FROM mimiciv_derived.chemistry WHERE hadm_id IN (", H, ")"))
dbDisconnect(con)

rows <- list()
add <- function(group, label, pop, ci, hr, lo, hi, p) rows[[length(rows) + 1]] <<- data.frame(group, label, n = nrow(pop), early = sum(pop$early), hr, lo, hi, p, ci)
rob <- function(group, label, pop, tvl, G = 24, era = FALSE, use_tv = TRUE, adjL0 = FALSE, amode = "L0", trunc = FALSE) {
  pop$id <- seq_len(nrow(pop)); cp <- build(pop, tvl, G = G, era = era, use_tv = use_tv, wlag = TRUE, amode = amode)
  if (trunc) { cA <- quantile(cp$w[cp$strategy == "Early" & cp$tstart > 0], .99); cB <- quantile(cp$w[cp$strategy == "Deferred"], .99)
    cp <- cp %>% mutate(w = if_else(strategy == "Early", if_else(tstart > 0, pmin(w, cA), w), pmin(w, cB))) }
  e <- est(cp, era = era, adjL0 = adjL0, pop = pop); add(group, label, pop, "robust", e[["HR"]], e[["lo"]], e[["hi"]], e[["p"]])
  cat(sprintf("%-46s %.3f (%.3f–%.3f)\n", label, e[["HR"]], e[["lo"]], e[["hi"]]))
}
bca <- function(th, bt, jk) { bt <- bt[is.finite(bt)]; jk <- jk[is.finite(jk)]; z0 <- qnorm(mean(bt < th) + 0.5 * mean(bt == th))
  jm <- mean(jk); a <- sum((jm - jk)^3) / (6 * sum((jm - jk)^2)^1.5); z <- qnorm(c(.025, .975))
  quantile(bt, pnorm(z0 + (z0 + z) / (1 - a * (z0 + z))), names = FALSE) }
bootrow <- function(group, label, pop, tvl, bf, jf) {
  pop$id <- seq_len(nrow(pop)); e <- est(build(pop, tvl, wlag = TRUE)); bt <- readRDS(file.path(BOOT, bf))[, "logHR"]; jk <- readRDS(file.path(BOOT, jf))[, "logHR"]
  ci <- exp(bca(log(e[["HR"]]), bt, jk)); bt <- bt[is.finite(bt)]; p <- min(1, 2 * min(mean(bt <= 0), mean(bt >= 0)))
  add(group, label, pop, "BCa", e[["HR"]], ci[1], ci[2], p); cat(sprintf("%-46s %.3f (%.3f–%.3f) BCa\n", label, e[["HR"]], ci[1], ci[2]))
}
G1 <- "Primary and cohort analyses"; G2 <- "Eligibility definition"; G3 <- "Model specification"; G4 <- "Grace period and outcome"; G5 <- "Urgent-indication criterion"
bootrow(G1, "Primary: SOFA 8–14", P0, TV0, "round3_boot_primary.rds", "round3_jk_primary.rds")
bootrow(G1, "Full cohort", P$primary, TV0, "round3_boot_full.rds", "round3_jk_full.rds")
bootrow(G2, "RRT-inclusive KDIGO staging", filter(P$rrt, eq_strict), TV$rrt, "round4_boot_rrt_sofa.rds", "round4_jk_rrt_sofa.rds")
rob(G3, "Deferred-arm weights, baseline only", P0, TV0, use_tv = FALSE)
rob(G3, "Early-arm weights, hourly initiation model", P0, TV0, amode = "haz")
rob(G3, "Weights truncated at 99th percentile", P0, TV0, trunc = TRUE)
rob(G3, "Calendar era", P0, TV0, era = TRUE)
rob(G3, "Outcome model baseline-adjusted", P0, TV0, adjL0 = TRUE)
lastv <- function(r, nm) r %>% filter(!is.na(v)) %>% inner_join(P0 %>% select(stay_id, hadm_id, t0), by = "hadm_id", relationship = "many-to-many") %>%
  filter(charttime <= t0, charttime >= t0 - hours(24)) %>% group_by(stay_id) %>% slice_max(charttime, n = 1, with_ties = FALSE) %>% ungroup() %>% select(stay_id, v) %>% rename(!!nm := v)
PAB <- P0 %>% left_join(lastv(ph, "ph_t0"), by = "stay_id") %>% left_join(lastv(ch %>% transmute(hadm_id, charttime, v = potassium), "k_t0"), by = "stay_id") %>%
  left_join(lastv(ch %>% transmute(hadm_id, charttime, v = bicarbonate), "hco3_t0"), by = "stay_id")
L0 <- c(setdiff(L0_full, "apsiii"), "ph_t0", "k_t0", "hco3_t0"); rob(G3, "Acid–base and electrolytes added", PAB, TV0)
L0 <- L0_full; rob(G3, "APS III retained", P0, TV0); L0 <- setdiff(L0_full, "apsiii")
rob(G4, "Grace period 12 h", P0, TV0, G = 12)
p48 <- P$primary %>% mutate(early = !is.na(h_crrt) & h_crrt < 48) %>% filter(eq_strict); rob(G4, "Grace period 48 h", p48, fix_tv(p48, tv48, G = 48), G = 48)
pin <- P0 %>% mutate(disch_h = as.numeric(difftime(dischtime, t0, units = "hours")), inh = !is.na(d_in) & d_in >= 0 & d_in <= 672, died28 = inh,
  end_h = case_when(inh ~ pmax(d_in, 0.5), !is.na(disch_h) & disch_h < 672 ~ pmax(disch_h, 0.5), TRUE ~ 672),
  h_crrt = if_else(!is.na(h_crrt) & h_crrt < end_h, h_crrt, NA_real_), early = !is.na(h_crrt) & h_crrt < 24)
rob(G4, "In-hospital death, censored at discharge", pin, TV0)
rob(G5, "Hyperkalaemia reconfirmed", filter(P$r1c, eq_strict), TV$r1c)
rob(G5, "Urine-output criterion added", filter(P$r14, eq_strict), TV$r14)
rob(G5, "Cut-offs lowered", filter(P$lib, eq_strict), TV$lib)
rob(G5, "Cut-offs raised", filter(P$cons, eq_strict), TV$cons)
out <- bind_rows(rows); write.csv(out, OUTCSV, row.names = FALSE, fileEncoding = "UTF-8"); cat("\nwrote", OUTCSV, "rows:", nrow(out), "\n")
