#!/usr/bin/env python3
"""
Comprehensive final review fix for result.docx.
Fixes: code/variable names, figure numbering, remaining Discussion content,
spelling consistency, supplementary duplicates, verbose annotations.
"""
import re
from docx import Document
from docx.shared import Pt

RESULT_PATH = "analysis/result.docx"
THIN = '\u2009'

doc = Document(RESULT_PATH)
changes = []

def norm(text):
    return text.replace(THIN, ' ')

def rewrite(para, new_text):
    if not para.runs:
        para.text = new_text
        return
    fn = para.runs[0].font.name
    fs = para.runs[0].font.size
    fb = para.runs[0].font.bold
    fi = para.runs[0].font.italic
    for r in para.runs:
        r.text = ""
    para.runs[0].text = new_text
    for r in para.runs:
        r.font.name = fn or "Times New Roman"
        r.font.size = fs or Pt(11)
        if fb is not None:
            r.font.bold = fb
        if fi is not None:
            r.font.italic = fi


for i, para in enumerate(doc.paragraphs):
    t = para.text.strip()
    tn = norm(t)
    if not t:
        continue
    modified = False
    new_t = para.text

    # ===================================================================
    # FIX GROUP 1: Remove code/variable/database references
    # ===================================================================

    # P006: Table 1 footnote — remove database-specific references
    if "mimiciv_derived.creatinine_baseline" in tn:
        new_t = (
            "\u00b9Continuous variables: median (IQR); categorical variables: n (%). "
            "\u00b2Vasopressor use is 100% by design (inclusion criterion). "
            "\u00b3Baseline serum creatinine defined as the lowest value within the 7 days preceding ICU admission. "
            "ESRD patients are excluded; CKD here represents non-ESRD chronic kidney disease. "
            "Data source: MIMIC-IV v3.1 (Beth Israel Deaconess Medical Center, 2008\u20132022). "
            "\u2074Norepinephrine-equivalent dose (NEE) reflects the last recorded value within the baseline window. "
            "Although vasopressor use was an inclusion requirement (hence 100% vasopressor use), "
            "quantitative dose data were unavailable for 17.5% of patients. "
            "Patients with missing NEE dose remained in the analytic cohort, "
            "as the exclusion criterion applied to the presence of vasopressor activity at t\u2080 "
            "(not to the quantitative dose value)."
        )
        modified = True
        changes.append(f"P{i:03d}: Removed database table names and ICD codes from Table 1 footnote")

    # P029: Table 4 footnote — remove code variable names
    if "sofa_baseline" in tn and "sofa_t" in tn and "Table 4" in tn or \
       ("CCW IPCW weights are derived from a time-varying censoring model" in tn):
        new_t = (
            "\u1d43Multivariable logistic regression with covariate adjustment "
            "(age, baseline SOFA score, APACHE III score, mechanical ventilation at t\u2080, "
            "Charlson comorbidity index); HC3 sandwich standard errors. "
            "Covariate-adjusted logistic regression was used for this binary secondary outcome "
            "because the IPCW weights, derived from a time-varying censoring model for the primary "
            "time-to-event analysis, are not directly applicable to a single-point binary outcome. "
            "Variables with substantial missingness (lactate, serum creatinine at t\u2080, vasopressor dose) "
            "were excluded from the adjustment model to preserve sample size. "
            "Renal recovery denominator: all patients (N = 3,112); death and ongoing RRT-dependence "
            "at day 28 classified as no recovery. "
            "RRT-free days reported as median [IQR]; a supplementary Wilcoxon rank-sum test was performed "
            "as an exploratory descriptive comparison (not used for causal inference). "
            "ICU LOS reported as median [IQR] and compared using Wilcoxon rank-sum test. "
            "Abbreviations: CI, confidence interval; ICU, intensive care unit; IQR, interquartile range; "
            "OR, odds ratio; RRT, renal replacement therapy; SCr, serum creatinine; "
            "SOFA, Sequential Organ Failure Assessment; t\u2080, index time. "
            "Note: the equipoise subgroup for secondary outcomes (renal recovery, RRT-free days) "
            "is defined using the SOFA score at ICU admission (n \u2248 1,952), "
            "which differs slightly from the primary analysis equipoise defined using SOFA at t\u2080 (n = 1,887)."
        )
        modified = True
        changes.append(f"P{i:03d}: Removed code variable names from Table 4 footnote")

    # P036: Table 5 footnote — remove R function names
    if "survival::finegray" in tn:
        new_t = (
            "sHR = subdistribution hazard ratio from IPCW-weighted Fine-Gray model "
            "(cluster-robust sandwich standard errors clustered on ICU stay). "
            "No covariate adjustment; confounding addressed through IPCW weighting. "
            "CIF = cumulative incidence function. "
            "Event 1 = 28-day renal recovery; Event 2 = death within 28 days (competing event). "
            "The equipoise subgroup (SOFA 8\u201314) is defined using SOFA at index time t\u2080, "
            "consistent with the pre-specified methods definition (n = 1,887). "
            "Event counts in this table are reported at the clone level "
            "(each patient is represented by two clones in the CCW framework), "
            "and are therefore approximately twice the patient-level event counts reported in the text "
            "(e.g., 434 clone-level recovery events \u2248 217 patient-level events \u00d7 2)."
        )
        modified = True
        changes.append(f"P{i:03d}: Removed R function names from Table 5 footnote")

    # P038: Figure 5 caption — remove code references + renumber to Figure 4
    if "Figure 5." in tn and "subdistribution hazard ratios" in tn:
        new_t = (
            "Figure 4.  Forest plot of IPCW-weighted subdistribution hazard ratios (sHR) "
            "from Fine-Gray competing risks models. Four estimates are shown: "
            "28-day renal recovery (full cohort and equipoise subgroup, SOFA 8\u201314) "
            "and 28-day mortality (full cohort and equipoise subgroup), "
            "each with the opposing event as the competing risk. "
            "Point estimates and 95% confidence intervals are derived from the IPCW-weighted "
            "Fine-Gray model. The \u2018Events\u2019 column reports events in the expanded clone dataset "
            "(each patient duplicated into two clones), which is approximately twice the "
            "patient-level event count reported in the text. "
            "The vertical dashed line marks sHR = 1 (no effect). "
            "Blue = renal recovery endpoint; Red = mortality endpoint."
        )
        modified = True
        changes.append(f"P{i:03d}: Cleaned Figure 5 caption (removed code refs, renumbered to Figure 4)")

    # P041: Weight diagnostics — remove variable name
    if "ipcw_w_trim" in tn and "well-behaved" in tn:
        new_t = new_t.replace("ipcw_w_trim, after", "after")
        new_t = re.sub(r'\(ipcw_w_trim,?\s*', '(', new_t)
        modified = True
        changes.append(f"P{i:03d}: Removed variable name from weight diagnostics")

    # P052: Table 6 footnote — remove code variable names
    if "sofa_t0 from CCW dataset" in tn or "A_early" in tn and "Table 6" in tn:
        if "IPCW, inverse probability" in tn and "Medians computed" in tn:
            new_t = (
                "IPCW, inverse probability of censoring weighting. HR, hazard ratio. "
                "All analyses use IPCW-weighted Cox proportional hazards models with cluster-robust "
                "sandwich standard errors (clustered on ICU stay). "
                "Medians computed in the full analytic cohort (n = 3,112): "
                "SOFA at t\u2080 = 9, creatinine = 1.90 mg/dL, APACHE III = 69. "
                "Note: the SOFA median used for stratification (9) reflects the SOFA score at t\u2080, "
                "which differs from the SOFA at ICU admission reported in Table 2 (median = 10). "
                "\u2018Events\u2019 refers to the number of outcome events observed in the early CRRT "
                "strategy arm within the CCW framework; total model events (both strategy arms combined) "
                "are larger."
            )
            modified = True
            changes.append(f"P{i:03d}: Removed code variable names from Table 6 footnote")

    # P054: Remove database table name
    if "MIMIC-IV chartevents table" in tn:
        new_t = new_t.replace("MIMIC-IV chartevents table", "MIMIC-IV database")
        new_t = new_t.replace(THIN + "MIMIC-IV chartevents table", THIN + "MIMIC-IV database")
        modified = True
        changes.append(f"P{i:03d}: Replaced database table name with generic reference")

    # ===================================================================
    # FIX GROUP 2: Figure renumbering (swap 4↔5)
    # ===================================================================

    # P034: Cross-reference "Figure 5" → "Figure 4"
    if "summarised visually in Figure 5" in tn:
        new_t = new_t.replace("Figure 5", "Figure 4")
        modified = True
        changes.append(f"P{i:03d}: Updated cross-reference Figure 5 → Figure 4")

    # P053: "Figure 4." caption → "Figure 5."
    if tn.startswith("Figure 4.") and "sensitivity" in tn.lower():
        new_t = new_t.replace("Figure 4.", "Figure 5.", 1)
        # Also remove the explanation about omitting crude HR
        omit_phrase = "The unadjusted (crude) Cox analysis is reported separately in Table 3"
        omit_idx = norm(new_t).find(omit_phrase)
        if omit_idx >= 0:
            new_t = new_t[:omit_idx].rstrip()
            if not new_t.endswith('.'):
                new_t += '.'
        modified = True
        changes.append(f"P{i:03d}: Renumbered Figure 4 → Figure 5, removed crude HR omission explanation")

    # P063: Supplementary S4 is duplicate of main text Figure 5 → remove
    if "Supplementary Figure S4" in tn and "Forest plot" in tn:
        new_t = ""
        modified = True
        changes.append(f"P{i:03d}: Removed duplicate Supplementary Figure S4 (= main text Figure 5)")

    # P062: Supplementary S3 is duplicate of main text Figure 3 → remove
    if "Supplementary Figure S3" in tn and "renal recovery" in tn.lower():
        new_t = ""
        modified = True
        changes.append(f"P{i:03d}: Removed duplicate Supplementary Figure S3 (= main text Figure 3)")

    # P064: Supplementary S5 → renumber to S3
    if "Supplementary Figure S5" in tn:
        new_t = new_t.replace("Supplementary Figure S5", "Supplementary Figure S3")
        # Also remove "This distribution motivated..." sentence
        moti_phrase = "This distribution motivated the adoption of the stacked bar chart format (Figure 2)."
        new_t_clean = norm(new_t)
        moti_idx = new_t_clean.find(moti_phrase)
        if moti_idx >= 0:
            new_t = new_t[:moti_idx].rstrip()
            if not new_t.endswith('.'):
                new_t += '.'
        modified = True
        changes.append(f"P{i:03d}: Renumbered S5 → S3, removed methodological justification")

    # ===================================================================
    # FIX GROUP 3: Remaining Discussion content
    # ===================================================================

    # P020: Remove "mathematically equivalent" explanation
    if "The two representations are mathematically equivalent" in tn:
        math_phrase = "The two representations are mathematically equivalent; the cumulative incidence format is included to facilitate direct visual comparison with the CIF curves for the secondary outcome (Figure 3)."
        new_t_n = norm(new_t)
        idx = new_t_n.find("The two representations are mathematically equivalent")
        if idx >= 0:
            end_marker = "(Figure 3)."
            end_idx = new_t_n.find(end_marker, idx)
            if end_idx >= 0:
                actual_end = end_idx + len(end_marker)
                new_t = (new_t[:idx].rstrip() + " " + new_t[actual_end:].lstrip()).strip()
                modified = True
                changes.append(f"P{i:03d}: Removed 'mathematically equivalent' explanation from Figure 1B caption")

    # P034: Remove "These estimates reflect the net effect..."
    if "These estimates reflect the net effect" in tn:
        idx = norm(new_t).find("These estimates reflect the net effect")
        if idx >= 0:
            end_marker = "competing risk of death."
            end_idx = norm(new_t).find(end_marker, idx)
            if end_idx >= 0:
                actual_end = end_idx + len(end_marker)
                new_t = (new_t[:idx].rstrip() + " " + new_t[actual_end:].lstrip()).strip()
                modified = True
                changes.append(f"P{i:03d}: Removed 'net effect' interpretation from Fine-Gray results")

    # ===================================================================
    # FIX GROUP 4: Spelling consistency (standardized → standardised)
    # ===================================================================
    if "standardized" in new_t.lower() and "standardised" not in new_t.lower():
        new_t = new_t.replace("standardized", "standardised")
        new_t = new_t.replace("Standardized", "Standardised")
        modified = True
        changes.append(f"P{i:03d}: Standardized → Standardised (spelling consistency)")

    # ===================================================================
    # FIX GROUP 5: Simplify verbose P046
    # ===================================================================
    if "sofa_t0 in the CCW dataset" in tn or \
       ("reflecting sofa_t0" in tn and "CCW time-varying covariate pipeline" in tn):
        new_t = (
            "3.5.3  Baseline severity-stratified analyses. "
            "To assess directional consistency of the treatment-outcome association across "
            "clinically distinct subgroups, we stratified the full analytic cohort at the cohort "
            "median for three baseline severity indicators: SOFA score at t\u2080 (median = 9), "
            "serum creatinine at t\u2080 (median = 1.90 mg/dL; <1% missing), "
            "and APACHE III score (median = 69)."
        )
        modified = True
        changes.append(f"P{i:03d}: Simplified severity-stratified intro (removed verbose variable explanations)")

    # Apply changes
    if modified:
        new_t = re.sub(r'\s{2,}', ' ', new_t).strip()
        rewrite(para, new_t)

doc.save(RESULT_PATH)

print("=" * 70)
print("ALL CHANGES:")
for c in changes:
    print(f"  \u2713 {c}")
print(f"\nTotal: {len(changes)} fixes applied")
