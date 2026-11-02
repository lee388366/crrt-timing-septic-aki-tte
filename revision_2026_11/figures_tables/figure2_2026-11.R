# Figure 2: loads the editable ggsurvfit objects used for the submitted figure (layout template, not distributed),
# replaces their data with the corrected IPCW-weighted Kaplan-Meier estimates (primary population, n = 1,469),
# and runs the original layout code verbatim. Colours, bands, risk table, theme and size are unchanged.
project_dir <- Sys.getenv("CRRT_PROJECT_ROOT")
if (nzchar(Sys.getenv("CRRT_FIG2_RLIB"))) .libPaths(c(Sys.getenv("CRRT_FIG2_RLIB"), .libPaths()))   # pinned ggsurvfit version used for the submitted figure
SCR <- Sys.getenv("CRRT_SCRIPT_DIR", getwd())
owd <- setwd(SCR); source("cp_core.R"); setwd(owd)
options(cp_trunc = TRUE)   # primary analysis: weights truncated at the 99th percentile (registered protocol)
library(ggsurvfit); library(patchwork)
PRIV <- file.path(dirname(SCR), "private_data")
out_dir <- file.path(Sys.getenv("CRRT_OUTPUT_DIR"), "figures")

# ---- Corrected primary analysis (same specification as round3/round7) ----
L0 <- setdiff(L0, "apsiii")
coh <- readRDS(file.path(PRIV, "cohort_round1.rds")); tvl <- readRDS(file.path(PRIV, "tv_hourly_round1.rds"))
pop <- filter(coh, eq_strict); pop$id <- seq_len(nrow(pop))
cp <- build(pop, tvl, wlag = TRUE)
cp$strategy <- factor(ifelse(cp$strategy == "Early", "Early CRRT", "Deferred CRRT"), levels = c("Deferred CRRT", "Early CRRT"))
fit <- survfit(Surv(tstart, tstop, event) ~ strategy, data = cp, weights = w)
cat("Check: weighted 28-day mortality risk", round(100 * (1 - summary(fit, times = 28, extend = TRUE)$surv), 1), "(see out/round7_report.txt)\n")

# ---- Load the original plot objects and replace their data ----
input_dir <- Sys.getenv("CRRT_FIG2_TEMPLATE_DIR")   # SA_AKI_KM_survival_editable.rds / SA_AKI_cumulative_risk_editable.rds
s <- readRDS(file.path(input_dir, "SA_AKI_KM_survival_editable.rds"))
r <- readRDS(file.path(input_dir, "SA_AKI_cumulative_risk_editable.rds"))
ds <- ggsurvfit(fit, type = "survival")$data; dr <- ggsurvfit(fit, type = "risk")$data   # same structure as the original objects (incl. survfit column used by the risk table)
cat("Original columns:", paste(names(s$data), collapse = ", "), "\nNew columns:", paste(names(ds), collapse = ", "), "\n")
cat("Original survfit column class:", class(s$data$survfit), "\n")
if (!"survfit" %in% names(ds)) {          # newer ggsurvfit versions omit this column: add it as in the original object
  if (is.list(s$data$survfit)) { ds$survfit <- rep(list(fit), nrow(ds)); dr$survfit <- rep(list(fit), nrow(dr)) }
  else { ds$survfit <- s$data$survfit[1]; dr$survfit <- r$data$survfit[1] }
}
cat("Original strata:", paste0("[", unique(as.character(s$data$strata)), "]"), "| new strata:", paste0("[", unique(as.character(ds$strata)), "]"), "\n")
relab <- function(d, ref) {               # map Early/Deferred to the original labels and factor levels
  old <- unique(as.character(ref$strata)); new <- as.character(d$strata)
  key <- function(x) ifelse(grepl("Early", x, ignore.case = TRUE), "E", "D")
  m <- setNames(old, key(old)); lab <- unname(m[key(new)])
  d$strata <- if (is.factor(ref$strata)) factor(lab, levels = levels(ref$strata)) else lab
  if ("strata_label" %in% names(ref)) d$strata_label <- ref$strata_label[1]
  d
}
ds <- relab(ds, s$data); dr <- relab(dr, r$data)
stopifnot(setequal(names(s$data), names(ds)), setequal(unique(as.character(s$data$strata)), unique(as.character(ds$strata))))
s$data <- ds[, names(s$data)]; r$data <- dr[, names(r$data)]

# ---- Layout code below is verbatim from the submitted-figure script (only output path/file name differ) ----
stopifnot(identical(s$data$time, r$data$time),
          identical(s$data$strata, r$data$strata),
          identical(s$data$n.risk, r$data$n.risk),
          max(abs(s$data$estimate + r$data$estimate - 1)) < 1e-14)
s <- s + labs(title = "A  Kaplan-Meier survival", x = "Days since t0",
              y = "Survival probability") +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28),
                     expand = expansion(mult = c(0.08, 0.05))) +
  theme(plot.title = element_text(size = 12, face = "bold"))
r <- r + labs(title = "B  Cumulative mortality risk", x = "Days since t0",
              y = "Cumulative all-cause mortality risk") +
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28),
                     expand = expansion(mult = c(0.08, 0.05))) +
  theme(plot.title = element_text(size = 12, face = "bold"))
panel_a <- ggsurvfit_build(s, combine_plots = FALSE)
panel_b <- ggsurvfit_build(r, combine_plots = FALSE)
panel_a[[2]]$layout$clip <- "off"
panel_b[[2]]$layout$clip <- "off"
combined <- wrap_plots(
  wrap_elements(full = panel_a[[1]], clip = FALSE),
  wrap_elements(full = panel_b[[1]], clip = FALSE),
  wrap_elements(full = panel_a[[2]], clip = FALSE),
  wrap_elements(full = panel_b[[2]], clip = FALSE),
  ncol = 2, heights = c(1, 0.16)
)
ggsave(file.path(out_dir, "03_Figure_2_Cumulative_Mortality_2026-11.pdf"), combined,
       device = grDevices::cairo_pdf, width = 11.5, height = 5.8,
       units = "in", bg = "white")
ggsave(file.path(out_dir, "03_Figure_2_Cumulative_Mortality_2026-11.png"), combined,
       device = ragg::agg_png, width = 11.5, height = 5.8,
       units = "in", dpi = 600, bg = "white")
cat("Saved to:", out_dir, "\n")
