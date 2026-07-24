# =============================================================================
# Step 11: Write result.docx — Section 3.1 Baseline Characteristics (Table 1)
# Font: Times New Roman (matches methods_section.docx)
# Table: 6 columns — Characteristic | Overall | Early | Deferred | SMD | p-value
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(flextable)
  library(officer)
})

# ── Global flextable defaults (Times New Roman) ───────────────────────────────
set_flextable_defaults(
  font.family = "Times New Roman",
  font.size   = 10,
  padding     = 3
)

# ─── Load & clean data ────────────────────────────────────────────────────────
master <- read.csv(
  "analysis/data/analysis_master.csv",
  stringsAsFactors = FALSE
)

fix_bool <- function(x) ifelse(x %in% c("t","TRUE","true","1"), TRUE,
                         ifelse(x %in% c("f","FALSE","false","0"), FALSE, NA))

master <- master %>%
  mutate(
    early_crrt       = fix_bool(early_crrt),
    died_28d         = fix_bool(died_28d),
    mv_baseline      = fix_bool(mv_baseline),
    dm               = fix_bool(dm),
    htn              = fix_bool(htn),
    chf              = fix_bool(chf),
    ckd              = fix_bool(ckd),
    liver_disease    = fix_bool(liver_disease),
    malignancy       = fix_bool(malignancy),
    copd             = fix_bool(copd),
    nephrotoxic_flag = fix_bool(nephrotoxic_flag),
    race_grp = case_when(
      grepl("^WHITE",    race) ~ "White",
      grepl("^BLACK",    race) ~ "Black/African American",
      grepl("^ASIAN",    race) ~ "Asian",
      grepl("^HISPANIC", race) ~ "Hispanic/Latino",
      TRUE                     ~ "Other/Unknown"
    ),
    race_grp = factor(race_grp,
      levels = c("White","Black/African American","Asian","Hispanic/Latino","Other/Unknown")),
    aki_path = case_when(
      aki_pathway == "creatinine_only" ~ "Creatinine criterion only",
      aki_pathway == "uo_only"         ~ "Urine output criterion only",
      aki_pathway == "both"            ~ "Both criteria",
      TRUE                             ~ NA_character_
    ),
    aki_path = factor(aki_path,
      levels = c("Creatinine criterion only","Urine output criterion only","Both criteria"))
  )

A  <- master %>% filter(early_crrt == TRUE)
B  <- master %>% filter(early_crrt == FALSE)
AL <- master

# ─── SMD helpers ─────────────────────────────────────────────────────────────
# Continuous: mean-based Cohen's d (standard in clinical epidemiology)
smd_cont <- function(xA, xB) {
  xA <- xA[!is.na(xA)]; xB <- xB[!is.na(xB)]
  if (length(xA) < 2 || length(xB) < 2) return("")
  pooled_sd <- sqrt(((length(xA)-1)*var(xA) + (length(xB)-1)*var(xB)) /
                    (length(xA)+length(xB)-2))
  if (pooled_sd == 0) return("0.000")
  sprintf("%.3f", abs(mean(xA) - mean(xB)) / pooled_sd)
}

# Binary (proportion): pooled-proportion formula
smd_bin <- function(xA, xB) {
  xA <- xA[!is.na(xA)]; xB <- xB[!is.na(xB)]
  pA <- mean(xA); pB <- mean(xB)
  denom <- sqrt((pA*(1-pA) + pB*(1-pB)) / 2)
  if (is.nan(denom) || denom == 0) return("0.000")
  sprintf("%.3f", abs(pA - pB) / denom)
}

# Section header rows: no SMD
smd_na <- ""

# ─── Cell value helpers ────────────────────────────────────────────────────────
miq <- function(x, d=1) {
  x <- x[!is.na(x)]
  sprintf("%s [%s, %s]",
          formatC(median(x), format="f", digits=d),
          formatC(quantile(x, 0.25), format="f", digits=d),
          formatC(quantile(x, 0.75), format="f", digits=d))
}

npct <- function(x, d=1) {
  n <- sum(x, na.rm=TRUE); N <- sum(!is.na(x))
  sprintf("%d (%.1f)", n, 100*n/N)
}

# p-values — superscript 'a' = Wilcoxon, 'b' = Fisher's exact, 'c' = chi-square
wp <- function(x, g) {
  p <- wilcox.test(x ~ g)$p.value
  lab <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  paste0(lab, "\u1d43")   # superscript 'a'
}

cp <- function(v, g, exact=FALSE) {
  tbl <- table(v, g, useNA="no")
  if (nrow(tbl) < 2) return("\u2014")
  if (exact) {
    p <- fisher.test(tbl, simulate.p.value=TRUE, B=2000)$p.value
    lab <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
    paste0(lab, "\u1d47")  # superscript 'b'
  } else {
    suppressWarnings(p <- chisq.test(tbl)$p.value)
    lab <- if (p < 0.001) "<0.001" else sprintf("%.3f", p)
    paste0(lab, "\u1d9c")  # superscript 'c'
  }
}

# ─── Build curated Table 1 — 6 columns ───────────────────────────────────────
# col order: Characteristic | Overall | Early | Deferred | SMD | p-value
rows <- list(
  # ── Demographics ────────────────────────────────────────────────────────────
  c("**Demographics**",              "",                "",               "",               smd_na, ""),
  c("Age, years",
    miq(AL$age),  miq(A$age),  miq(B$age),
    smd_cont(A$age, B$age),   wp(master$age, master$early_crrt)),

  c("Female sex, n (%)",
    npct(AL$gender=="F"), npct(A$gender=="F"), npct(B$gender=="F"),
    smd_bin(A$gender=="F", B$gender=="F"),  cp(master$gender, master$early_crrt)),

  c("Body weight, kg",
    miq(AL$weight_kg,1), miq(A$weight_kg,1), miq(B$weight_kg,1),
    smd_cont(A$weight_kg, B$weight_kg),   wp(master$weight_kg, master$early_crrt)),

  c("Race/ethnicity, n (%)", "", "", "", smd_na,
    cp(master$race_grp, master$early_crrt, exact=TRUE)),
  c("  White",
    npct(AL$race_grp=="White"), npct(A$race_grp=="White"), npct(B$race_grp=="White"),
    smd_bin(A$race_grp=="White", B$race_grp=="White"), ""),
  c("  Black/African American",
    npct(AL$race_grp=="Black/African American"),
    npct(A$race_grp=="Black/African American"),
    npct(B$race_grp=="Black/African American"),
    smd_bin(A$race_grp=="Black/African American", B$race_grp=="Black/African American"), ""),
  c("  Asian",
    npct(AL$race_grp=="Asian"), npct(A$race_grp=="Asian"), npct(B$race_grp=="Asian"),
    smd_bin(A$race_grp=="Asian", B$race_grp=="Asian"), ""),
  c("  Hispanic/Latino",
    npct(AL$race_grp=="Hispanic/Latino"), npct(A$race_grp=="Hispanic/Latino"), npct(B$race_grp=="Hispanic/Latino"),
    smd_bin(A$race_grp=="Hispanic/Latino", B$race_grp=="Hispanic/Latino"), ""),
  c("  Other/Unknown",
    npct(AL$race_grp=="Other/Unknown"), npct(A$race_grp=="Other/Unknown"), npct(B$race_grp=="Other/Unknown"),
    smd_bin(A$race_grp=="Other/Unknown", B$race_grp=="Other/Unknown"), ""),

  # ── Illness Severity ────────────────────────────────────────────────────────
  c("**Illness Severity at t\u2080**",  "",  "",  "",  smd_na, ""),
  c("SOFA score",
    miq(AL$sofa_baseline,0), miq(A$sofa_baseline,0), miq(B$sofa_baseline,0),
    smd_cont(A$sofa_baseline, B$sofa_baseline),  wp(master$sofa_baseline, master$early_crrt)),
  c("APACHE III score",
    miq(AL$apsiii_score,0), miq(A$apsiii_score,0), miq(B$apsiii_score,0),
    smd_cont(A$apsiii_score, B$apsiii_score),    wp(master$apsiii_score, master$early_crrt)),
  c("GCS (minimum)",
    miq(AL$gcs_baseline,0), miq(A$gcs_baseline,0), miq(B$gcs_baseline,0),
    smd_cont(A$gcs_baseline, B$gcs_baseline),    wp(master$gcs_baseline, master$early_crrt)),
  c("Charlson comorbidity index",
    miq(AL$cci,0), miq(A$cci,0), miq(B$cci,0),
    smd_cont(A$cci, B$cci),                      wp(master$cci, master$early_crrt)),

  # ── Vital Signs ─────────────────────────────────────────────────────────────
  c("**Vital Signs at t\u2080**",  "",  "",  "",  smd_na, ""),
  c("MAP, mmHg",
    miq(AL$map,1), miq(A$map,1), miq(B$map,1),
    smd_cont(A$map, B$map),   wp(master$map, master$early_crrt)),
  c("Heart rate, bpm",
    miq(AL$hr,0), miq(A$hr,0), miq(B$hr,0),
    smd_cont(A$hr, B$hr),     wp(master$hr, master$early_crrt)),
  c("Respiratory rate, /min",
    miq(AL$rr,0), miq(A$rr,0), miq(B$rr,0),
    smd_cont(A$rr, B$rr),     wp(master$rr, master$early_crrt)),
  c("Temperature, \u00b0C",
    miq(AL$temp_c,1), miq(A$temp_c,1), miq(B$temp_c,1),
    smd_cont(A$temp_c, B$temp_c),  wp(master$temp_c, master$early_crrt)),
  c("SpO\u2082, %",
    miq(AL$spo2,0), miq(A$spo2,0), miq(B$spo2,0),
    smd_cont(A$spo2, B$spo2),  wp(master$spo2, master$early_crrt)),

  # ── Laboratory Values ───────────────────────────────────────────────────────
  c("**Laboratory Values at t\u2080**",  "",  "",  "",  smd_na, ""),
  c("Serum creatinine (t\u2080), mg/dL",
    miq(AL$scr,2), miq(A$scr,2), miq(B$scr,2),
    smd_cont(A$scr, B$scr),  wp(master$scr, master$early_crrt)),
  c("Baseline SCr (pre-AKI), mg/dL",
    miq(AL$scr_baseline,2), miq(A$scr_baseline,2), miq(B$scr_baseline,2),
    smd_cont(A$scr_baseline, B$scr_baseline),  wp(master$scr_baseline, master$early_crrt)),
  c("BUN, mg/dL",
    miq(AL$bun,0), miq(A$bun,0), miq(B$bun,0),
    smd_cont(A$bun, B$bun),  wp(master$bun, master$early_crrt)),
  c("Lactate, mmol/L",
    miq(AL$lactate,2), miq(A$lactate,2), miq(B$lactate,2),
    smd_cont(A$lactate, B$lactate),  wp(master$lactate, master$early_crrt)),
  c("pH",
    miq(AL$ph,2), miq(A$ph,2), miq(B$ph,2),
    smd_cont(A$ph, B$ph),    wp(master$ph, master$early_crrt)),
  c("Bicarbonate, mEq/L",
    miq(AL$hco3,0), miq(A$hco3,0), miq(B$hco3,0),
    smd_cont(A$hco3, B$hco3),  wp(master$hco3, master$early_crrt)),
  c("Potassium, mEq/L",
    miq(AL$k,1), miq(A$k,1), miq(B$k,1),
    smd_cont(A$k, B$k),      wp(master$k, master$early_crrt)),
  c("WBC, \u00d710\u2079/L",
    miq(AL$wbc,1), miq(A$wbc,1), miq(B$wbc,1),
    smd_cont(A$wbc, B$wbc),  wp(master$wbc, master$early_crrt)),
  c("Hemoglobin, g/dL",
    miq(AL$hb,1), miq(A$hb,1), miq(B$hb,1),
    smd_cont(A$hb, B$hb),    wp(master$hb, master$early_crrt)),
  c("Platelets, \u00d710\u2079/L",
    miq(AL$plt,0), miq(A$plt,0), miq(B$plt,0),
    smd_cont(A$plt, B$plt),  wp(master$plt, master$early_crrt)),
  c("INR",
    miq(AL$inr,2), miq(A$inr,2), miq(B$inr,2),
    smd_cont(A$inr, B$inr),  wp(master$inr, master$early_crrt)),

  # ── Organ Support ───────────────────────────────────────────────────────────
  c("**Organ Support at t\u2080**",  "",  "",  "",  smd_na, ""),
  c("Mechanical ventilation, n (%)",
    npct(AL$mv_baseline), npct(A$mv_baseline), npct(B$mv_baseline),
    smd_bin(A$mv_baseline, B$mv_baseline),  cp(master$mv_baseline, master$early_crrt)),
  c("NEE dose, \u00b5g/kg/min",
    miq(AL$ne_eq_baseline,2), miq(A$ne_eq_baseline,2), miq(B$ne_eq_baseline,2),
    smd_cont(A$ne_eq_baseline, B$ne_eq_baseline),  wp(master$ne_eq_baseline, master$early_crrt)),
  c("Urine output (last 6h), mL/kg/h",
    miq(AL$uo_mlkgph_t0_last6h,2), miq(A$uo_mlkgph_t0_last6h,2), miq(B$uo_mlkgph_t0_last6h,2),
    smd_cont(A$uo_mlkgph_t0_last6h, B$uo_mlkgph_t0_last6h),  wp(master$uo_mlkgph_t0_last6h, master$early_crrt)),
  c("Cumulative crystalloid, mL/kg",
    miq(AL$crystalloid_ml_kg,1), miq(A$crystalloid_ml_kg,1), miq(B$crystalloid_ml_kg,1),
    smd_cont(A$crystalloid_ml_kg, B$crystalloid_ml_kg),  wp(master$crystalloid_ml_kg, master$early_crrt)),
  c("Nephrotoxin exposure, n (%)",
    npct(AL$nephrotoxic_flag), npct(A$nephrotoxic_flag), npct(B$nephrotoxic_flag),
    smd_bin(A$nephrotoxic_flag, B$nephrotoxic_flag),  cp(master$nephrotoxic_flag, master$early_crrt)),

  # ── Comorbidities ───────────────────────────────────────────────────────────
  c("**Comorbidities**",  "",  "",  "",  smd_na, ""),
  c("Diabetes mellitus, n (%)",
    npct(AL$dm),  npct(A$dm),  npct(B$dm),
    smd_bin(A$dm, B$dm),  cp(master$dm, master$early_crrt)),
  c("Hypertension, n (%)",
    npct(AL$htn), npct(A$htn), npct(B$htn),
    smd_bin(A$htn, B$htn),  cp(master$htn, master$early_crrt)),
  c("Congestive heart failure, n (%)",
    npct(AL$chf), npct(A$chf), npct(B$chf),
    smd_bin(A$chf, B$chf),  cp(master$chf, master$early_crrt)),
  c("Chronic kidney disease, n (%)",
    npct(AL$ckd), npct(A$ckd), npct(B$ckd),
    smd_bin(A$ckd, B$ckd),  cp(master$ckd, master$early_crrt)),
  c("Liver disease, n (%)",
    npct(AL$liver_disease), npct(A$liver_disease), npct(B$liver_disease),
    smd_bin(A$liver_disease, B$liver_disease),  cp(master$liver_disease, master$early_crrt)),
  c("Malignancy, n (%)",
    npct(AL$malignancy), npct(A$malignancy), npct(B$malignancy),
    smd_bin(A$malignancy, B$malignancy),  cp(master$malignancy, master$early_crrt)),
  c("COPD, n (%)",
    npct(AL$copd), npct(A$copd), npct(B$copd),
    smd_bin(A$copd, B$copd),  cp(master$copd, master$early_crrt)),

  # ── Primary Outcome (unadjusted) ─────────────────────────────────────────────
  c("**Primary Outcome (unadjusted)**",  "",  "",  "",  smd_na, ""),
  c("28-day all-cause mortality, n (%)",
    npct(AL$died_28d), npct(A$died_28d), npct(B$died_28d),
    smd_bin(A$died_28d, B$died_28d),  cp(master$died_28d, master$early_crrt))
)

tbl <- as.data.frame(do.call(rbind, rows), stringsAsFactors=FALSE)
names(tbl) <- c("Characteristic",
                sprintf("Overall\n(N = %d)", nrow(master)),
                "Early CRRT\n(\u224824h, n = 500)",
                "Deferred CRRT\n(>24h, n = 2,612)",
                "SMD",
                "p-value")

# ─── Build flextable ─────────────────────────────────────────────────────────
header_rows   <- which(grepl("^\\*\\*", tbl$Characteristic))
indent_rows   <- which(grepl("^  ", tbl$Characteristic))  # race sub-rows

# Clean display text
tbl_display <- tbl
tbl_display$Characteristic <- gsub("^\\*\\*|\\*\\*$", "", tbl_display$Characteristic)
tbl_display$Characteristic <- trimws(tbl_display$Characteristic)

# Restore race indentation for display (we'll handle via padding)
race_label_rows <- which(tbl_display$Characteristic %in%
  c("White","Black/African American","Asian","Hispanic/Latino","Other/Unknown"))

ft <- flextable(tbl_display) %>%
  # ── Font: Times New Roman throughout ────────────────────────────────────────
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  # ── Header formatting ────────────────────────────────────────────────────────
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  # ── Section header rows: bold, left-aligned, light grey background ───────────
  bold(i = header_rows, j = 1) %>%
  bg(i = header_rows, bg = "#EEEEEE", part = "body") %>%
  # ── Indent race sub-rows ─────────────────────────────────────────────────────
  padding(i = race_label_rows, j = 1, padding.left = 20) %>%
  # ── Column alignment ─────────────────────────────────────────────────────────
  align(j = 1, align = "left",   part = "body") %>%
  align(j = 2:6, align = "center", part = "body") %>%
  # ── SMD column: bold if ≥0.1 (imbalanced) ────────────────────────────────────
  # (applied after table build — highlight notable imbalance)
  color(i = which(suppressWarnings(as.numeric(tbl_display$SMD)) >= 0.1),
        j = "SMD", color = "#B22222", part = "body") %>%
  # ── Column widths (total ≈ 7.4 inches for A4 with 1-inch margins) ────────────
  width(j = 1, width = 2.30) %>%   # Characteristic
  width(j = 2, width = 1.15) %>%   # Overall
  width(j = 3, width = 1.25) %>%   # Early
  width(j = 4, width = 1.30) %>%   # Deferred
  width(j = 5, width = 0.60) %>%   # SMD
  width(j = 6, width = 0.75) %>%   # p-value
  # ── Borders ──────────────────────────────────────────────────────────────────
  border_outer(part = "all",    border = fp_border(width = 1.2)) %>%
  hline(part = "header",         border = fp_border(width = 1.2)) %>%
  hline(i = header_rows,         border = fp_border(width = 0.5, color = "#AAAAAA")) %>%
  # ── Caption ──────────────────────────────────────────────────────────────────
  set_caption(
    caption = as_paragraph(
      as_chunk("Table 1. ",
               props = fp_text(font.family="Times New Roman", font.size=10, bold=TRUE)),
      as_chunk("Baseline Characteristics of the Study Cohort, by Treatment Group",
               props = fp_text(font.family="Times New Roman", font.size=10, bold=FALSE))
    )
  )

# ─── Text properties for officer paragraphs ──────────────────────────────────
tp_h1   <- fp_text(font.family="Times New Roman", font.size=14, bold=TRUE)
tp_h2   <- fp_text(font.family="Times New Roman", font.size=12, bold=TRUE)
tp_body <- fp_text(font.family="Times New Roman", font.size=11)
tp_fn   <- fp_text(font.family="Times New Roman", font.size=9)

# fp_par only supports: text.align, padding.*, line_spacing, border.*, word_style
pp_body <- fp_par(line_spacing=1.15, text.align="left")
pp_head <- fp_par(line_spacing=1.15, text.align="left")
pp_fn   <- fp_par(line_spacing=1.0,  text.align="left")

# ─── Build narrative ─────────────────────────────────────────────────────────
mort_pct <- round(100*mean(master$died_28d, na.rm=TRUE), 1)
mort_n   <- sum(master$died_28d, na.rm=TRUE)

narrative <- paste0(
  "Between 2008 and 2022, a total of 3,112 ICU admissions at Beth Israel Deaconess ",
  "Medical Center met the eligibility criteria for inclusion in the analytic cohort. ",
  "Of these, 500 patients (16.1%) initiated CRRT within 24 hours of t\u2080 ",
  "(early CRRT group), and 2,612 patients (83.9%) did not receive CRRT within the ",
  "grace period (deferred CRRT group). Overall 28-day all-cause mortality was ",
  sprintf("%.1f%% (%d/3,112).", mort_pct, mort_n),

  "\n\nAt the index time t\u2080, patients who subsequently received early CRRT exhibited ",
  "markedly higher acute illness severity compared with those managed with watchful ",
  "waiting. The early CRRT group had higher SOFA scores (median 12 [IQR 10\u201314] vs ",
  "9 [7\u201312]; p\u2009<\u20090.001), higher APACHE III scores (84 [IQR 69\u201398] vs ",
  "66 [IQR 52\u201383]; p\u2009<\u20090.001), greater lactate concentrations ",
  "(4.45 [2.40\u20139.83] vs 2.30 [1.50\u20133.90] mmol/L; p\u2009<\u20090.001), ",
  "more severe metabolic acidosis (pH 7.26 [7.20\u20137.35] vs 7.34 [7.28\u20137.39]; ",
  "p\u2009<\u20090.001), lower urine output (0.09 [0.02\u20130.25] vs 0.22 [0.12\u20130.35] ",
  "mL/kg/h; p\u2009<\u20090.001), and higher norepinephrine-equivalent vasopressor doses ",
  "(0.32 [0.15\u20130.53] vs 0.12 [0.05\u20130.26] \u00b5g/kg/min; p\u2009<\u20090.001) at t\u2080. ",
  "The early CRRT group also had a higher prevalence of mechanical ventilation ",
  "(88.4% vs 70.7%; p\u2009<\u20090.001). Pre-AKI baseline serum creatinine ",
  "(median 0.70 mg/dL in both groups; p\u2009=\u20090.305) and the prevalence of chronic ",
  "kidney disease (30.0% vs 27.8%; p\u2009=\u20090.351) were similar between groups, ",
  "confirming that the acute severity differential was not attributable to pre-existing ",
  "renal impairment. Large standardised mean differences (SMD\u2009\u22650.10) were observed ",
  "for SOFA score, APACHE III score, lactate, vasopressor dose, pH, and urine output, ",
  "reflecting substantial baseline imbalance that motivates the application of ",
  "Clone-Censor-Weighting to achieve causal comparability. ",
  "Detailed baseline characteristics are presented in Table\u00a01."
)

footnote_txt <- paste0(
  "\u1d43Wilcoxon rank-sum test (all continuous variables failed the Shapiro\u2013Wilk ",
  "normality test, \u03b1\u2009=\u20090.05). ",
  "\u1d47Fisher\u2019s exact test (simulated p-value, B\u2009=\u20092,000). ",
  "\u1d9cChi-square test. ",
  "SMD: standardised mean difference (Cohen\u2019s d for continuous, pooled-proportion ",
  "formula for binary variables); SMD\u2009\u22650.10 (shown in red) indicates meaningful ",
  "imbalance. ",
  "Abbreviations: APACHE\u00a0III, Acute Physiology Score\u00a0III; CKD, chronic kidney disease; ",
  "COPD, chronic obstructive pulmonary disease; GCS, Glasgow Coma Scale; INR, ",
  "international normalised ratio; IQR, interquartile range; MAP, mean arterial pressure; ",
  "NEE, norepinephrine-equivalent dose; SCr, serum creatinine; SMD, standardised mean ",
  "difference; SOFA, Sequential Organ Failure Assessment; t\u2080, index time; WBC, white blood ",
  "cell count."
)

# ─── Assemble document using template (Times New Roman styles) ────────────────
template_path <- "analysis/result_template.docx"
out_path      <- "analysis/result.docx"

doc <- read_docx(template_path)

doc <- doc %>%
  # Heading 1 — explicit Times New Roman bold 14pt
  body_add_fpar(fpar(ftext("Results", tp_h1), fp_p = pp_head)) %>%
  # Heading 2 — explicit Times New Roman bold 12pt
  body_add_fpar(fpar(ftext("3.1  Study Population and Baseline Characteristics", tp_h2),
                     fp_p = pp_head)) %>%
  # Narrative — Normal style from template (inherits TNR 11pt)
  body_add_par(narrative, style = "Normal") %>%
  # Spacer
  body_add_par("", style = "Normal") %>%
  # Table
  body_add_flextable(ft) %>%
  # Footnote — explicit 9pt TNR
  body_add_fpar(fpar(ftext(footnote_txt, tp_fn), fp_p = pp_fn))

print(doc, target = out_path)

cat("\n\u2713 result.docx saved:", out_path, "\n")
cat("  Font:    Times New Roman 11pt body / 10pt table\n")
cat("  Table:   ", nrow(tbl_display), "rows \u00d7 6 columns (+ SMD column)\n")
cat("  Notable imbalance (SMD\u22650.10) in",
    sum(suppressWarnings(as.numeric(tbl_display$SMD)) >= 0.1, na.rm=TRUE),
    "variables\n")
