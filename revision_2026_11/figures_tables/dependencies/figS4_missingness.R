# Fig. S4-type: missing-data profile (DataExplorer::plot_missing style, re-themed).
# Original version (October 2026 submission); figureS2S3S4_2026-11.R replaces its inputs with the corrected cohort (n = 2,461).
# Hourly panel: person-hours 0-23 after LOCF, before median imputation.
# Only covariates entering the IPCW models; bands follow DataExplorer defaults: Good <5%, OK 5-40%, Bad 40-80%, Remove >80%.
source("/workspace/tte_bmj/scripts/00_common.R")
b  <- readRDS(file.path(ROOT, "data/corrected_inputs/baseline.rds"))
tv <- readRDS(file.path(ROOT, "data/corrected_inputs/tv_hourly.rds"))
stopifnot(nrow(b) == 2463, nrow(tv) == 59112)
L0 <- c(age = "Age", sofa_t0 = "SOFA at t\u2080", apsiii = "APS III", creatinine_baseline = "Baseline creatinine",
        nee_at_t0 = "NEE at t\u2080", lactate_t0 = "Lactate at t\u2080", uo_rate_t0 = "Urine output rate at t\u2080",
        map_t0 = "MAP at t\u2080", cci = "Charlson index")
LK <- c(map_h = "MAP", nee_h = "NEE", lactate_h = "Lactate", sofa_h = "SOFA", kdigo_h = "KDIGO stage",
        creatinine_h = "Creatinine", creatinine_delta = "Creatinine change from baseline", fluid_balance_h = "Cumulative fluid balance")
BAND_COLS <- c("Good" = "#009E73", "OK" = "#E69F00", "Bad" = "#D55E00", "Remove" = "#000000")
prof <- function(d, vars) {
  tibble(var = names(vars), label = unname(vars), n_miss = sapply(names(vars), function(v) sum(is.na(d[[v]]))), n = nrow(d)) %>%
    mutate(pct = 100 * n_miss / n,
           band = factor(cut(pct, c(-Inf, 5, 40, 80, Inf), right = FALSE, labels = names(BAND_COLS)), levels = names(BAND_COLS)),
           label = reorder(label, pct))
}
pA <- prof(b, L0); pB <- prof(tv, LK)
write.csv(bind_rows(mutate(pA, panel = "A_baseline_patients"), mutate(pB, panel = "B_person_hours_0_23")),
          file.path(OUT, "FigS4_missingness_table.csv"), row.names = FALSE)
bar <- function(p, ttl) {
  ggplot(p, aes(pct, label, fill = band)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.2f%%", pct)), hjust = -0.1, size = 2.7, colour = COL_INK) +
    scale_fill_manual(values = BAND_COLS, drop = TRUE, name = "Band") +
    scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 10), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
    labs(title = ttl, x = "Missing", y = NULL) +
    theme_bmj(9) + theme(plot.title = element_text(size = 9.5), plot.title.position = "plot", legend.title = element_text(size = 8.5))
}
p <- (bar(pA, "a  Baseline covariates") / bar(pB, "b  Time-varying covariates (0\u201323 h)")) +
  plot_layout(guides = "collect", heights = c(9, 8)) & theme(legend.position = "right")
save_fig(p, "FigS4_missingness", 6.5, 5.8)
