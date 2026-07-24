# =============================================================================
# Step 10: Table 1 — Baseline Characteristics by Treatment Group
# Study: Early vs Deferred CRRT in Septic Shock + KDIGO Stage 3 AKI
# Output: analysis/result.docx (Section 3.1)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tableone)
  library(flextable)
  library(officer)
})

# ─── 1. Load and fix data types ───────────────────────────────────────────────
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)

fix_bool <- function(x) ifelse(x %in% c("t","TRUE","true","1"), TRUE,
                         ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))

master <- master %>%
  mutate(
    early_crrt         = fix_bool(early_crrt),
    died_28d           = fix_bool(died_28d),
    mv_baseline        = fix_bool(mv_baseline),
    vaso_baseline      = fix_bool(vaso_baseline),
    colloid_flag       = fix_bool(colloid_flag),
    nephrotoxic_flag   = fix_bool(nephrotoxic_flag),
    rrt_ongoing_day28  = fix_bool(rrt_ongoing_day28),
    renal_recovery_28d = fix_bool(renal_recovery_28d),
    dm                 = fix_bool(dm),
    htn                = fix_bool(htn),
    chf                = fix_bool(chf),
    ckd                = fix_bool(ckd),
    liver_disease      = fix_bool(liver_disease),
    malignancy         = fix_bool(malignancy),
    copd               = fix_bool(copd),
    # AKI pathway as factor
    aki_pathway        = factor(aki_pathway,
                                levels = c("creatinine_only","uo_only","both"),
                                labels = c("Creatinine only","Urine output only","Both")),
    # treatment label
    trt_label = factor(ifelse(early_crrt, "Early CRRT (≤24h)", "Deferred CRRT (>24h)"),
                       levels = c("Early CRRT (≤24h)", "Deferred CRRT (>24h)"))
  )

cat(sprintf("N = %d  |  Early CRRT: %d (%.1f%%)  |  Deferred: %d (%.1f%%)\n",
    nrow(master),
    sum(master$early_crrt, na.rm=TRUE),
    100*mean(master$early_crrt, na.rm=TRUE),
    sum(!master$early_crrt, na.rm=TRUE),
    100*mean(!master$early_crrt, na.rm=TRUE)))

cat(sprintf("28-day mortality: %.1f%%\n", 100*mean(master$died_28d, na.rm=TRUE)))

# ─── 2. Define Table 1 variables ─────────────────────────────────────────────
cont_vars <- c(
  "age", "weight_kg",
  "sofa_baseline", "gcs_baseline", "apsiii_score",
  "map", "hr", "rr", "temp_c", "spo2",
  "scr_baseline", "scr",
  "lactate", "ph", "hco3", "k", "bun", "wbc", "hb", "plt", "inr",
  "uo_mlkgph_t0_last6h", "crystalloid_ml_kg", "ne_eq_baseline",
  "cci"
)

cat_vars <- c(
  "gender", "race",
  "dm", "htn", "chf", "ckd", "liver_disease", "malignancy", "copd",
  "mv_baseline", "colloid_flag", "nephrotoxic_flag",
  "aki_pathway"
)

all_vars <- c(cont_vars, cat_vars)

# ─── 3. Normality check (Shapiro-Wilk, n≤5000) ───────────────────────────────
nonnormal <- c()
for (v in cont_vars) {
  x <- master[[v]][!is.na(master[[v]])]
  if (length(x) >= 3 && length(x) <= 5000) {
    p <- shapiro.test(x[1:min(5000, length(x))])$p.value
    if (p < 0.05) nonnormal <- c(nonnormal, v)
  }
}
cat("\nNon-normal continuous vars (", length(nonnormal), "):", paste(nonnormal, collapse=", "), "\n")

# ─── 4. Sparsity check for chi-square vs Fisher ───────────────────────────────
exact_vars <- c()
for (v in cat_vars) {
  tbl <- table(master[[v]], master$early_crrt, useNA="no")
  if (nrow(tbl) >= 2 && ncol(tbl) >= 2) {
    suppressWarnings(exp <- chisq.test(tbl)$expected)
    if (any(exp < 5)) exact_vars <- c(exact_vars, v)
  }
}
cat("Fisher's exact vars:", paste(exact_vars, collapse=", "), "\n\n")

# ─── 5. Create TableOne ───────────────────────────────────────────────────────
tab1 <- CreateTableOne(
  vars      = all_vars,
  strata    = "trt_label",
  data      = master,
  factorVars = cat_vars,
  test      = TRUE,
  addOverall = TRUE
)

# Print to console
print(tab1, nonnormal = nonnormal, exact = exact_vars, quote = FALSE, noSpaces = TRUE)

# Extract as data frame for export
tab1_df <- print(tab1,
                 nonnormal   = nonnormal,
                 exact       = exact_vars,
                 quote       = FALSE,
                 noSpaces    = TRUE,
                 printToggle = FALSE,
                 showAllLevels = FALSE)

# Convert row names to column
tab1_df <- as.data.frame(tab1_df, stringsAsFactors = FALSE)
tab1_df <- cbind(Variable = rownames(tab1_df), tab1_df)
rownames(tab1_df) <- NULL

# Rename columns
names(tab1_df)[1] <- "Variable"
# Columns: Variable, Overall (n=3112), Early CRRT (<=24h) (n=500), Deferred (>24h) (n=2612), p, test
cat("Table columns:", names(tab1_df), "\n")

# ─── 6. SMD calculation ───────────────────────────────────────────────────────
smd_vals <- ExtractSmd(tab1)
smd_df   <- data.frame(
  Variable = rownames(smd_vals),
  SMD      = round(as.numeric(smd_vals[,1]), 3),
  stringsAsFactors = FALSE
)

# ─── 7. Save table objects for use in result.docx generation ─────────────────
saveRDS(list(
  tab1_df  = tab1_df,
  smd_df   = smd_df,
  nonnormal = nonnormal,
  exact_vars = exact_vars,
  n_total   = nrow(master),
  n_early   = sum(master$early_crrt, na.rm=TRUE),
  n_deferred = sum(!master$early_crrt, na.rm=TRUE),
  pct_early = round(100*mean(master$early_crrt, na.rm=TRUE),1),
  mort_overall = round(100*mean(master$died_28d, na.rm=TRUE),1),
  mort_early   = round(100*mean(master$died_28d[master$early_crrt], na.rm=TRUE),1),
  mort_deferred= round(100*mean(master$died_28d[!master$early_crrt], na.rm=TRUE),1)
), "analysis/data/table1_results.rds")

cat("\n✓ Table 1 results saved to data/table1_results.rds\n")
