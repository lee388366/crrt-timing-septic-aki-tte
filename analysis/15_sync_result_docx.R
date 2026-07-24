# =============================================================================
# Script 15: Sync result.docx — Final content reconciliation
# Actions:
#   1. Replace Figure 3 placeholder with proper stacked-bar figure caption
#   2. Add Figure 4 reference (CIF) after Table 4 / Fine-Gray section
#   3. Add Section 3.4 IPCW Weight Diagnostics (weight stats + Love Plot)
#   4. Append Supplementary Materials inventory at document end
# =============================================================================

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

doc_path <- "analysis/result.docx"
doc      <- read_docx(doc_path)

TNR <- "Times New Roman"
fp_h2  <- fp_text(font.family = TNR, font.size = 12, bold = TRUE)
fp_h3  <- fp_text(font.family = TNR, font.size = 11, bold = TRUE)
fp_bd  <- fp_text(font.family = TNR, font.size = 11, bold = TRUE)
fp_nm  <- fp_text(font.family = TNR, font.size = 11, bold = FALSE)
fp_it  <- fp_text(font.family = TNR, font.size = 11, bold = FALSE, italic = TRUE)
fp_sm  <- fp_text(font.family = TNR, font.size = 9,  bold = FALSE, color = "grey30")
fp_smb <- fp_text(font.family = TNR, font.size = 9,  bold = TRUE)

# Helper ──────────────────────────────────────────────────────────────────────
add_nm <- function(doc, txt, space_before = 0, space_after = 4) {
  body_add_fpar(doc,
    fpar(ftext(txt, fp_nm),
         fp_p = fp_par(padding.top = space_before, padding.bottom = space_after)),
    style = "Normal")
}
add_it <- function(doc, txt) {
  body_add_fpar(doc,
    fpar(ftext(txt, fp_it),
         fp_p = fp_par(padding.top = 0, padding.bottom = 3)),
    style = "Normal")
}
add_sm <- function(doc, txt) {
  body_add_fpar(doc,
    fpar(ftext(txt, fp_sm),
         fp_p = fp_par(padding.top = 0, padding.bottom = 3)),
    style = "Normal")
}

# ── 1. Replace Figure 3 placeholder with proper caption ──────────────────────
# Current text: "[Figure 3 — see analysis/figures/Figure3_RRTfree_days.pdf]"
# Replace with a proper italic figure caption + figure note

doc <- cursor_reach(doc, keyword = "Figure 3 — see analysis/figures/Figure3_RRTfree_days")
doc <- body_remove(doc)   # remove the old placeholder

doc <- body_add_fpar(doc,
  fpar(
    ftext("Figure 3.  ", fp_bd),
    ftext(paste0(
      "Distribution of RRT-Free Days Alive at Day 28 by Treatment Group. ",
      "Stacked bar chart showing the proportion of patients in each of three mutually ",
      "exclusive outcome categories: died or ongoing RRT (= 0 days, red), partial ",
      "RRT-free survival (1-27 days, orange), and full 28-day RRT-free survival ",
      "(= 28 days, blue). Early CRRT group: n = 500; Deferred CRRT group: n = 2,612. ",
      "Wilcoxon p < 0.001. Note that no patient in the early CRRT group achieved ",
      "full 28-day RRT-free survival, as all had initiated CRRT before day 24."
    ), fp_it),
    fp_p = fp_par(padding.top = 4, padding.bottom = 6)
  ),
  style = "Normal"
)
cat("Figure 3 caption updated.\n")

# ── 2. Add Figure 4 reference after Fine-Gray footnote ───────────────────────
doc <- cursor_reach(doc,
  keyword = "Event 1 = 28-day renal recovery; Event 2 = death within 28 days")
doc <- cursor_forward(doc)   # move past the footnote

doc <- body_add_fpar(doc,
  fpar(
    ftext("Figure 4.  ", fp_bd),
    ftext(paste0(
      "Cumulative Incidence Functions (CIF) for 28-Day Renal Recovery. ",
      "Panels A and B show the cumulative incidence of renal recovery (Event 1) ",
      "in the presence of the competing risk of death (Event 2), estimated by the ",
      "nonparametric estimator of Aalen et al. (cmprsk package). ",
      "Panel A: full analytic cohort (n = 3,112). ",
      "Panel B: equipoise subgroup (SOFA 8-14, n = 1,952). ",
      "Grey solid line = Early CRRT; Red dashed line = Deferred CRRT."
    ), fp_it),
    fp_p = fp_par(padding.top = 4, padding.bottom = 8)
  ),
  style = "Normal"
)
cat("Figure 4 reference added.\n")

# ── 3. Add Section 3.4 IPCW Weight Diagnostics ───────────────────────────────
doc <- body_add_fpar(doc,
  fpar(ftext("3.4  IPCW Weight Diagnostics and Covariate Balance", fp_h2),
       fp_p = fp_par(padding.top = 10, padding.bottom = 4)),
  style = "heading 2"
)

doc <- add_nm(doc, paste0(
  "To assess the adequacy of the clone-censor-weighting procedure, we examined the ",
  "distribution of the cumulative stabilized IPCW weights and evaluated covariate ",
  "balance before and after weighting."
), space_before = 2)

doc <- add_nm(doc, paste0(
  "3.4.1  Weight distribution. The cumulative stabilized IPCW weights (ipcw_w_trim, ",
  "after 99th-percentile truncation) were well-behaved: median 0.981 (mean 0.959; ",
  "range 0.036-1.976). The narrow distribution and mean close to 1.0 indicate that ",
  "the weighting procedure introduced minimal variance inflation and that no extreme ",
  "weights required further stabilization beyond the applied trimming."
))

doc <- add_nm(doc, paste0(
  "3.4.2  Covariate balance (Love plot). Prior to IPCW weighting, substantial ",
  "imbalance was observed between the early and deferred CRRT groups across all ",
  "10 assessed baseline covariates (only 2/10 had SMD < 0.10), reflecting strong ",
  "confounding by indication. The largest pre-weighting imbalances were seen for ",
  "SOFA score at t0 (SMD = 0.817), lactate at t0 (SMD = 0.692), APACHE III score ",
  "(SMD = 0.687), and norepinephrine-equivalent dose (SMD = 0.644). After IPCW ",
  "weighting, all 10 covariates achieved SMD < 0.10 (maximum post-weighting ",
  "SMD = 0.076 for creatinine at t0), confirming that the CCW procedure effectively ",
  "balanced baseline characteristics between treatment arms. The Love plot displaying ",
  "pre- and post-weighting SMDs for all assessed covariates is presented in ",
  "Supplementary Figure S2."
))
cat("Section 3.4 added.\n")

# ── 4. Append Supplementary Materials inventory ───────────────────────────────
doc <- body_add_fpar(doc,
  fpar(ftext("Supplementary Materials", fp_h2),
       fp_p = fp_par(padding.top = 14, padding.bottom = 6)),
  style = "heading 2"
)

doc <- add_nm(doc, paste0(
  "The following supplementary materials are available online:"
), space_before = 2)

# Supplementary table
doc <- body_add_fpar(doc,
  fpar(
    ftext("Supplementary Figure S1.  ", fp_smb),
    ftext(paste0(
      "Weighted Kaplan-Meier survival curves for 28-day all-cause mortality, ",
      "full analytic cohort (n = 3,112). IPCW-weighted survival estimates for ",
      "early CRRT (Strategy A) versus deferred CRRT (Strategy B), with 95% ",
      "confidence bands. [See: analysis/figures/FigS1_KM_FullCohort.pdf]"
    ), fp_sm),
    fp_p = fp_par(padding.top = 3, padding.bottom = 3)
  ),
  style = "Normal"
)

doc <- body_add_fpar(doc,
  fpar(
    ftext("Supplementary Figure S2.  ", fp_smb),
    ftext(paste0(
      "Love plot: absolute standardized mean differences (SMD) for 10 baseline ",
      "covariates before and after IPCW weighting. The vertical dashed red line ",
      "marks the SMD = 0.10 balance threshold. Green triangles represent ",
      "post-weighting SMDs; dark circles represent pre-weighting SMDs. ",
      "All post-weighting SMDs fall below 0.10, confirming adequate balance. ",
      "[See: analysis/figures/FigS2_LovePlot_IPCW_Balance.pdf]"
    ), fp_sm),
    fp_p = fp_par(padding.top = 3, padding.bottom = 3)
  ),
  style = "Normal"
)

doc <- body_add_fpar(doc,
  fpar(
    ftext("Supplementary Figure S3.  ", fp_smb),
    ftext(paste0(
      "Diagnostic histogram of RRT-free days alive at day 28 (full cohort, n = 3,112), ",
      "illustrating the bimodal distribution with spikes at 0 (death or ongoing RRT) ",
      "and 28 (full RRT-free survival). ",
      "[See: analysis/figures/FigS3_RFD_Distribution_Diagnostic.pdf]"
    ), fp_sm),
    fp_p = fp_par(padding.top = 3, padding.bottom = 6)
  ),
  style = "Normal"
)

cat("Supplementary Materials section added.\n")

# ── Save ─────────────────────────────────────────────────────────────────────
print(doc, target = doc_path)
cat("result.docx saved.\n")
