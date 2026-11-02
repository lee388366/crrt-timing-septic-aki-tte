# Supplemental Table S5 (primary population by observed treatment): the nine rows not produced by the Table 1 script (heart rate, respiratory rate, temperature, SpO2, pH, mechanical ventilation, crystalloid, nephrotoxins, hypertension)
# Data: rich SQL output from run_S14_rich.R (windows at or before t0), summarised by observed group; aggregate output only. Also produces counts for Tables S8, S10 and S14.
suppressPackageStartupMessages({library(dplyr); library(bit64)})   # stay_id from the SQL output is integer64
PRIV <- file.path(Sys.getenv("CRRT_REVISION_DIR"), "private_data"); OUTD <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "tables")
d   <- readRDS(file.path(PRIV, "S14_rich", "PRIVATE_original_select.rds"))
coh <- readRDS(file.path(PRIV, "cohort_round1.rds")) %>% select(stay_id, eq_strict, early)
d <- inner_join(mutate(d, stay_id = as.character(stay_id)), mutate(coh, stay_id = as.character(stay_id)), by = "stay_id")
stopifnot(nrow(d) == 2461)
f1 <- function(x, k) formatC(x, format = "f", digits = k, big.mark = ",")
cont <- function(x, k) sprintf("%s [%s, %s]", f1(median(x, na.rm = TRUE), k), f1(quantile(x, .25, na.rm = TRUE, type = 7), k), f1(quantile(x, .75, na.rm = TRUE, type = 7), k))
binr <- function(x) { x <- x %in% c(1, TRUE); sprintf("%s (%.1f)", formatC(sum(x), big.mark = ","), 100 * mean(x)) }
pfmt <- function(p) ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, "<0.01", sprintf("%.2f", p)))
smd_c <- function(a, b) abs(mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
smd_b <- function(a, b) { p1 <- mean(a %in% c(1, TRUE)); p2 <- mean(b %in% c(1, TRUE)); abs(p1 - p2) / sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2) }
spec <- list(c("Heart rate, bpm", "hr", "c", 1), c("Respiratory rate, /min", "rr", "c", 1), c("Temperature, °C", "temp_c", "c", 1),
             c("SpO₂, %", "spo2", "c", 1), c("pH", "ph", "c", 2), c("Mechanical ventilation, n (%)", "mv_baseline", "b", 0),
             c("Cumulative crystalloid, mL/kg", "crystalloid_ml_kg", "c", 1), c("Nephrotoxin exposure, n (%)", "nephrotoxic_flag", "b", 0),
             c("Hypertension, n (%)", "htn", "b", 0))
out <- list()
for (pop in c("primary", "full")) {
  p <- if (pop == "primary") filter(d, eq_strict) else d
  E <- filter(p, early); D <- filter(p, !early)
  for (s in spec) {
    x <- p[[s[2]]]; a <- E[[s[2]]]; b <- D[[s[2]]]
    if (s[3] == "c") { k <- as.integer(s[4]); row <- c(cont(x, k), cont(a, k), cont(b, k), sprintf("%.3f", smd_c(a, b)), pfmt(wilcox.test(a, b)$p.value)) }
    else row <- c(binr(x), binr(a), binr(b), sprintf("%.3f", smd_b(a, b)), pfmt(suppressWarnings(chisq.test(table(c(a, b) %in% c(1, TRUE), rep(1:2, c(length(a), length(b)))))$p.value)))
    out[[length(out) + 1]] <- data.frame(pop = pop, label = s[1], overall = row[1], early = row[2], deferred = row[3], smd = row[4], p = row[5],
                                         missing = sprintf("%d (%.1f%%)", sum(is.na(x)), 100 * mean(is.na(x))))
  }
}
r <- do.call(rbind, out); write.csv(r, file.path(OUTD, "S5_fill_9rows.csv"), row.names = FALSE); print(r, right = FALSE)

# ---- Supplemental Table S10 (mortality ascertainment) counts ----
c2 <- readRDS(file.path(PRIV, "cohort_round1.rds"))
stopifnot(all(c("dod", "deathtime", "died28") %in% names(c2)))
inh <- !is.na(c2$deathtime) & as.numeric(difftime(c2$deathtime, c2$t0, units = "hours")) <= 672
s10 <- data.frame(metric = c("total", "any_dod", "inhosp_deathtime", "dod_no_inhosp", "dod_only_death_28d", "died28"),
  value = c(nrow(c2), sum(!is.na(c2$dod)), sum(!is.na(c2$deathtime)), sum(!is.na(c2$dod) & is.na(c2$deathtime)),
            sum(c2$died28 & !inh), sum(c2$died28)))
write.csv(s10, file.path(OUTD, "S10_ascertainment_counts.csv"), row.names = FALSE); print(s10)

# ---- Supplemental Table S8 (observed-group ICU length of stay from t0, days) ----
stopifnot("icu_los_from_t0" %in% names(c2))
los <- do.call(rbind, lapply(c("primary", "full"), function(pop) {
  p <- if (pop == "primary") filter(c2, eq_strict) else c2
  a <- p$icu_los_from_t0[p$early]; b <- p$icu_los_from_t0[!p$early]
  q <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
  data.frame(pop = pop, early = q(a), deferred = q(b), p = sprintf("%.3f", wilcox.test(a, b)$p.value), units_check_max = round(max(p$icu_los_from_t0, na.rm = TRUE), 1))
}))
write.csv(los, file.path(OUTD, "S8_icu_los.csv"), row.names = FALSE); print(los)

# ---- Supplemental Table S14: missingness of model-input variables (as in Table 1) and KDIGO stage 3 pathway ----
bl <- readRDS(file.path(Sys.getenv("CRRT_INPUT_DIR"), "ch24_q4_corrected_inputs/baseline.rds"))
c3 <- left_join(c2 %>% select(-any_of("stage3_path")), bl %>% select(stay_id, stage3_path), by = "stay_id")
mv <- function(v) sprintf("%s (%.1f%%)", formatC(sum(is.na(c3[[v]])), big.mark = ","), 100 * mean(is.na(c3[[v]])))
s14 <- data.frame(item = c("miss_creatinine_t0", "miss_lactate_t0", "miss_uo_rate_t0", "miss_map_t0", paste0("path_", names(table(c3$stage3_path, useNA = "ifany")))),
                  value = c(mv("creatinine_t0"), mv("lactate_t0"), mv("uo_rate_t0"), mv("map_t0"), as.character(as.vector(table(c3$stage3_path, useNA = "ifany")))))
write.csv(s14, file.path(OUTD, "S14_analysis_vars.csv"), row.names = FALSE); print(s14)
