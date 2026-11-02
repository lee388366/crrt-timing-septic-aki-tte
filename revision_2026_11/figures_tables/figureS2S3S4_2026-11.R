# Supplemental Figures S2/S3 (covariate means over time, unweighted and weighted) and S4 (missingness), using the layout functions in dependencies/
# S2/S3 changes: (1) data -> corrected cohort and final weights (cp_core::build, wlag = TRUE, amode = "L0", L0 without APS III);
#   (2) "Weighted" uses the time-varying weight in effect for each clone at that hour; (3) APS III panel removed.
# S4 changes: inputs -> corrected cohort (N = 2,461) and round-1 hourly data; APS III removed. Plotting functions unchanged.
GS    <- file.path(Sys.getenv("CRRT_REPO_DIR"), "revision_2026_11/figures_tables/dependencies")
source(file.path(Sys.getenv("CRRT_REPO_DIR"), "revision_2026_11/analysis/cp_core.R"))
options(cp_trunc = TRUE)   # primary analysis: weights truncated at the 99th percentile
L0_MODEL <- setdiff(L0, "apsiii")
cm <- readLines(file.path(GS, "00_common.R"), warn = FALSE); eval(parse(text = cm))
OUT <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "figures")                         # output directory
PRIV <- file.path(Sys.getenv("CRRT_REVISION_DIR"), "private_data")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds")); tvl <- readRDS(file.path(PRIV, "tv_hourly_round1.rds"))
stopifnot(nrow(coh) == 2461)

# ---------------- S2 / S3 ----------------
src <- readLines(file.path(GS, "figS2S3_covariate_means_over_time.R"), warn = FALSE)
keep <- function(name) { i <- grep(paste0("^", name, " <- function"), src); j <- i
  while (!grepl("^}", src[j]) && !(j == i && grepl("}\\s*$", src[j]))) j <- j + 1; src[i:j] }
for (f in c("wmean", "mini", "make_fig")) eval(parse(text = keep(f)))
L0 <- c(age = "Age (years)", sofa_t0_strict = "SOFA at t₀ (score)", creatinine_baseline = "Baseline Cr (mg/dL)",
        nee_at_t0 = "NEE at t₀ (µg/kg/min)", lactate_t0 = "Lactate at t₀ (mmol/L)",
        uo_rate_t0 = "UO at t₀ (mL/kg/h)", map_t0 = "MAP at t₀ (mmHg)", cci = "CCI (score)")
LK <- c(sofa_h = "SOFA (score)", map_h = "MAP (mmHg)", nee_h = "NEE (µg/kg/min)", lactate_h = "Lactate (mmol/L)",
        creatinine_h = "Cr (mg/dL)", kdigo_h = "KDIGO stage", fluid_balance_h = "Fluid balance (mL)")
HOURS_L0 <- 0:167
tv <- tvl %>% distinct(stay_id, hour_index, .keep_all = TRUE)
mk <- function(pop) {
  Lsave <- L0; L0 <<- L0_MODEL; cp <- build(pop, tvl, wlag = TRUE, amode = "L0"); L0 <<- Lsave
  cp %>% mutate(arm = factor(ifelse(strategy == "Early", "Early CRRT", "Deferred CRRT"), c("Deferred CRRT", "Early CRRT"))) %>%
    left_join(pop %>% select(id, stay_id, all_of(names(L0))), by = "id")
}
calc <- function(d) {     # time-varying weights: at hour h use the counting-process row covering h and its weight
  at_h <- function(h) d[d$tstart * 24 <= h & d$tstop * 24 > h, ]
  l0 <- lapply(HOURS_L0, function(h) { r <- at_h(h)
    bind_rows(lapply(names(L0), function(v) r %>% group_by(arm) %>%
      summarise(Unweighted = mean(.data[[v]], na.rm = TRUE), Weighted = wmean(.data[[v]], w), n = n(), .groups = "drop") %>% mutate(var = v, hour = h))) }) %>% bind_rows()
  lk <- lapply(0:23, function(h) { r <- at_h(h) %>% left_join(filter(tv, hour_index == h), by = "stay_id")
    bind_rows(lapply(names(LK), function(v) r %>% group_by(arm) %>%
      summarise(Unweighted = mean(.data[[v]], na.rm = TRUE), Weighted = wmean(.data[[v]], w), n = sum(!is.na(.data[[v]])), .groups = "drop") %>% mutate(var = v, hour = h))) }) %>% bind_rows()
  bind_rows(l0, lk) %>% pivot_longer(c(Unweighted, Weighted), names_to = "weighting", values_to = "mean")
}
pP <- filter(coh, eq_strict); pP$id <- seq_len(nrow(pP)); pF <- coh; pF$id <- seq_len(nrow(pF))
make_fig(mk(pP), "primary SOFA 8–14 cohort", nrow(pP), "FigS2_covariate_means_primary_SOFA8_14", "Fig. S2")
make_fig(mk(pF), "full cohort", nrow(pF), "FigS3_covariate_means_full_cohort", "Fig. S3")

# ---------------- S4 ----------------
s4 <- readLines(file.path(GS, "figS4_missingness.R"), warn = FALSE)
s4 <- s4[!grepl("^source\\(", s4)]
s4 <- sub('^b  <- readRDS\\(.*$', 'b <- coh %>% mutate(sofa_t0 = sofa_t0_strict)', s4)
s4 <- sub('^tv <- readRDS\\(.*$', 'tv <- tvl %>% filter(stay_id %in% coh$stay_id, hour_index %in% 0:23) %>% distinct(stay_id, hour_index, .keep_all = TRUE)', s4)
s4 <- sub('^stopifnot\\(nrow\\(b\\) == 2463, nrow\\(tv\\) == 59112\\)', 'stopifnot(nrow(b) == 2461); cat("S4 person-hours:", nrow(tv), "\\n")', s4)
s4 <- sub('apsiii = "APS III", ', '', s4, fixed = TRUE)
stopifnot(!any(grepl("apsiii", s4)), any(grepl("^b <- coh", s4)), any(grepl("^tv <- tvl", s4)))
eval(parse(text = s4))
tab <- read.csv(file.path(OUT, "FigS4_missingness_table.csv")); print(tab[, c("panel", "label", "n_miss", "n", "pct")])
