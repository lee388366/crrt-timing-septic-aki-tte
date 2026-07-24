"""
Fix all 13 identified inconsistencies between methods_section.docx and result.docx.

CRITICAL:
  #1  M-76: LRT → Wald test for interaction
  #2  R-13: RMST values update to corrected bootstrap CIs
  #3  R-45: E-value recalculation (CI crosses null → E=1.0; add full cohort)
  #4  R-81/82, M-89: C4 numbers to analytic cohort (N=3112)

MAJOR:
  #6  M-78: SMD formula clarification
  #8  M-88: "pre-specified" → "post-hoc"
  #9  R-02: SOFA median clarification (at ICU admission)
  #10 Merge errata into main text, remove duplication

MINOR:
  #13 R-50: interaction p-value → Wald, consistent precision
"""

import shutil
from docx import Document
from docx.oxml.ns import qn


def delete_paragraph(paragraph):
    p = paragraph._element
    p.getparent().remove(p)


def replace_in_para(para, old, new):
    full = para.text
    if old not in full:
        return False
    new_text = full.replace(old, new, 1)
    if not para.runs:
        return False
    para.runs[0].text = new_text
    for r in para.runs[1:]:
        r.text = ""
    return True


# ── Backup originals ──────────────────────────────────────────────────────────
shutil.copy2("methods_section.docx", "methods_section.docx.prefinal")
shutil.copy2("analysis/result.docx", "analysis/result.docx.prefinal")


# ═══════════════════════════════════════════════════════════════════════════════
#  PART 1: methods_section.docx
# ═══════════════════════════════════════════════════════════════════════════════
mdoc = Document("methods_section.docx")
mp = mdoc.paragraphs
changes_m = 0

# ── #1  M-76: LRT → Wald ─────────────────────────────────────────────────────
for p in mp:
    if "A likelihood ratio test comparing models with and without the interaction term" in p.text:
        ok = replace_in_para(
            p,
            "A likelihood ratio test comparing models with and without the interaction term "
            "was used to evaluate whether the treatment effect varied across the severity "
            "spectrum; a threshold of p < 0.10 was applied given the well-documented low "
            "power of interaction tests in observational studies of this size.",

            "A Wald test based on the cluster-robust (sandwich) standard error of the "
            "interaction coefficient was used to evaluate whether the treatment effect "
            "varied across the severity spectrum; the Wald test was preferred over the "
            "likelihood ratio test because the chi-squared approximation for the weighted "
            "partial log-likelihood does not hold under the pseudo-population framework "
            "created by IPCW weights. A threshold of p < 0.10 was applied given the "
            "well-documented low power of interaction tests in observational studies of "
            "this size."
        )
        if ok:
            changes_m += 1
            print("[M] #1 Fixed: LRT → Wald in M-76")
        break

# ── #6  M-78: SMD formula clarification ──────────────────────────────────────
for p in mp:
    if "SMDs of <0.10 were considered indicative of adequate balance" in p.text:
        ok = replace_in_para(
            p,
            "SMDs of <0.10 were considered indicative of adequate balance.",
            "Weighted SMDs were computed as the difference in IPCW-weighted means "
            "divided by the pooled unweighted standard deviation, following the convention "
            "for balance diagnostics in inverse probability weighting analyses; "
            "SMDs of <0.10 were considered indicative of adequate balance."
        )
        if ok:
            changes_m += 1
            print("[M] #6 Fixed: SMD formula clarification in M-78")
        break

# ── #3 partial: M-72 E-value CI-crosses-null convention ──────────────────────
for p in mp:
    if "This approximation assumes rare outcomes" in p.text:
        ok = replace_in_para(
            p,
            "This approximation assumes rare outcomes",
            "When the 95% confidence interval includes the null (HR = 1), the E-value "
            "for the confidence interval is defined as 1.0, because no unmeasured "
            "confounding is needed to shift the interval to include the null. "
            "This approximation assumes rare outcomes"
        )
        if ok:
            changes_m += 1
            print("[M] #3 partial: E-value CI convention added in M-72")
        break

# ── M-68: Add bootstrap merge detail ─────────────────────────────────────────
for p in mp:
    if "B = 1,000 replicates; percentile method)." in p.text:
        ok = replace_in_para(
            p,
            "B = 1,000 replicates; percentile method).",
            "B = 1,000 replicates; percentile method). Each bootstrap iteration used a "
            "merge-based resampling approach to ensure that patients drawn multiple "
            "times were represented with the appropriate frequency, preserving the "
            "with-replacement property essential for valid variance estimation."
        )
        if ok:
            changes_m += 1
            print("[M] Bootstrap merge detail added to M-68")
        break

# ── #8  M-88: pre-specified → post-hoc ───────────────────────────────────────
for p in mp:
    if "As a pre-specified sensitivity analysis, the Strategy A" in p.text:
        ok = replace_in_para(p, "As a pre-specified sensitivity analysis", "As a post-hoc sensitivity analysis")
        if ok:
            changes_m += 1
            print("[M] #8 Fixed: pre-specified → post-hoc in M-88")
        break

# ── #4 partial: M-89 C4 numbers ──────────────────────────────────────────────
for p in mp:
    if "identified 144 additional patients (4.5%)" in p.text:
        ok = replace_in_para(p, "identified 144 additional patients (4.5%)", "identified 138 additional patients (4.4%)")
        if ok:
            replace_in_para(p, "1.7% Early vs 5.0% Deferred", "1.8% Early vs 4.9% Deferred")
            changes_m += 1
            print("[M] #4 Fixed: C4 numbers in M-89")

# ── #10 partial: Remove addendum header and merge markers ─────────────────────
for p in mp:
    if p.text.strip() == "Methodological Addendum (2026-04-12)":
        delete_paragraph(p)
        changes_m += 1
        print("[M] #10 Deleted addendum header M-85")
        break

# Remove M-86 (bootstrap addendum) - now redundant with M-68 update
for p in mp:
    if p.text.startswith("Bootstrap resampling for RMST confidence intervals."):
        delete_paragraph(p)
        changes_m += 1
        print("[M] #10 Deleted bootstrap addendum M-86 (merged into M-68)")
        break

# Remove M-87 (interaction testing addendum) - now redundant with M-76 update
for p in mp:
    if p.text.startswith("Interaction testing. Effect modification by baseline SOFA"):
        delete_paragraph(p)
        changes_m += 1
        print("[M] #10 Deleted interaction addendum M-87 (merged into M-76)")
        break

# Remove addendum label from M-88 (IPCW sensitivity)
for p in mp:
    if p.text.startswith("IPCW denominator sensitivity analysis. As a post-hoc"):
        ok = replace_in_para(
            p,
            "IPCW denominator sensitivity analysis. As a post-hoc sensitivity analysis, the Strategy A",
            "Sensitivity analysis — IPCW denominator specification. As a post-hoc sensitivity analysis, the Strategy A"
        )
        if ok:
            changes_m += 1
            print("[M] #10 Reformatted M-88 as regular sensitivity section")
        break

# Reformat M-89 (mortality ascertainment) into proper limitation paragraph
for p in mp:
    if p.text.startswith("Mortality ascertainment. The primary analysis used"):
        ok = replace_in_para(
            p,
            "Mortality ascertainment. The primary analysis used",
            "Sensitivity analysis — Mortality ascertainment. The primary analysis used"
        )
        if ok:
            changes_m += 1
            print("[M] #10 Reformatted M-89 as sensitivity paragraph")
        break

mdoc.save("methods_section.docx")
print(f"\n✓ methods_section.docx: {changes_m} changes saved.\n")


# ═══════════════════════════════════════════════════════════════════════════════
#  PART 2: analysis/result.docx
# ═══════════════════════════════════════════════════════════════════════════════
rdoc = Document("analysis/result.docx")
rp = rdoc.paragraphs
changes_r = 0

# ── #9  R-02: SOFA at ICU admission clarification ────────────────────────────
for p in rp:
    if "the median SOFA score was 10 (IQR 7" in p.text:
        ok = replace_in_para(p, "the median SOFA score was 10", "the median SOFA score at ICU admission was 10")
        if ok:
            changes_r += 1
            print("[R] #9 Fixed: SOFA 'at ICU admission' clarification in R-02")
        break

# ── #2  R-13: RMST corrected values ──────────────────────────────────────────
for p in rp:
    if "-0.62 days (Bootstrap 95% CI -1.80 to 0.50 days)" in p.text:
        ok1 = replace_in_para(p,
            "-0.62 days (Bootstrap 95% CI -1.80 to 0.50 days)",
            "-0.63 days (Bootstrap 95% CI -2.07 to 0.84 days)")
        ok2 = replace_in_para(p,
            "-2.13 days (Bootstrap 95% CI -3.16 to -1.17)",
            "-2.12 days (Bootstrap 95% CI -3.30 to -0.79)")
        if ok1 or ok2:
            changes_r += 1
            print("[R] #2 Fixed: RMST values updated in R-13")
        break

# ── #3  R-45: E-value fix ────────────────────────────────────────────────────
for p in rp:
    if "The corresponding E-value for the point estimate was 1.236" in p.text:
        ok = replace_in_para(
            p,
            "The corresponding E-value for the point estimate was 1.236, and the E-value "
            "for the confidence interval bound closest to the null was 1.413.",

            "The corresponding E-value for the point estimate was 1.237. Because the "
            "95% confidence interval crossed the null (lower bound 0.915 < 1.0), the "
            "E-value for the confidence interval is 1.0, indicating that no unmeasured "
            "confounding is needed to shift the confidence interval to include the null. "
            "In the full analytic cohort (HR = 1.164, 95% CI 1.047 to 1.293), the "
            "E-value for the point estimate was 1.601 and the E-value for the "
            "confidence interval lower bound was 1.269, indicating that an unmeasured "
            "confounder would need to be associated with both early CRRT and 28-day "
            "mortality by a risk ratio of at least 1.27 to explain away the statistical "
            "significance of the full-cohort result."
        )
        if ok:
            changes_r += 1
            print("[R] #3 Fixed: E-value corrected in R-45")
        break

# ── #13 R-50: Interaction p-value precision + Wald label ─────────────────────
for p in rp:
    if "interaction HR per 1-unit SOFA above median = 1.018" in p.text:
        ok = replace_in_para(p, "p = 0.215)", "Wald p = 0.214)")
        if ok:
            changes_r += 1
            print("[R] #13 Fixed: interaction p → 'Wald p = 0.214'")
        break

# ── #4  R-81: C4 numbers ─────────────────────────────────────────────────────
for p in rp:
    if "9/517 = 1.7%" in p.text:
        ok = replace_in_para(p, "9 missed deaths in the Early CRRT group (9/517 = 1.7%) and "
            "135 missed deaths in the Deferred CRRT group (135/2,716 = 5.0%)",
            "9 missed deaths in the Early CRRT group (9/500 = 1.8%) and "
            "129 missed deaths in the Deferred CRRT group (129/2,612 = 4.9%)")
        if ok:
            changes_r += 1
            print("[R] #4 Fixed: C4 numbers in R-81")
        break

# ── #4  R-82: C4 total ───────────────────────────────────────────────────────
for p in rp:
    if "A total of 144 patients (4.5% of the cohort)" in p.text:
        ok = replace_in_para(p, "A total of 144 patients (4.5% of the cohort)",
                            "A total of 138 patients (4.4% of the cohort)")
        if ok:
            replace_in_para(p, "(1.7% vs 5.0%)", "(1.8% vs 4.9%)")
            changes_r += 1
            print("[R] #4 Fixed: C4 total in R-82")

# ── #10 Merge errata into main text ──────────────────────────────────────────
# E1 and E2 corrections are now reflected in main text.
# Convert errata section: keep S1/S2 as proper subsections, remove E1/E2.

# Delete errata title
for p in rp:
    if p.text.strip() == "ERRATA & UPDATED RESULTS (2026-04-12)":
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted errata title R-66")
        break

# Delete errata intro paragraph
for p in rp:
    if "The following corrections and supplementary analyses were performed on 2026-04-12" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted errata intro R-67")
        break

# Delete E1 header
for p in rp:
    if p.text.strip() == "E1. Corrected Bootstrap RMST Confidence Intervals":
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E1 header R-68")
        break

# Delete E1 issue
for p in rp:
    if "The original bootstrap resampling for RMST confidence intervals used R" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E1 issue R-69")
        break

# Delete E1 fix
for p in rp:
    if "The bootstrap was re-implemented using merge()" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E1 fix R-70")
        break

# Delete E1 impact
for p in rp:
    if "The corrected confidence intervals are slightly wider" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E1 impact R-71")
        break

# Delete E2 header
for p in rp:
    if p.text.strip() == "E2. Corrected Interaction Test (SOFA x Treatment)":
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E2 header R-72")
        break

# Delete E2 issue
for p in rp:
    if "The original interaction test between baseline SOFA score and treatment strategy used a likelihood ratio test" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E2 issue R-73")
        break

# Delete E2 fix
for p in rp:
    if "The Wald test based on the cluster-robust (sandwich) standard error was used instead" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E2 fix R-74")
        break

# Delete E2 impact
for p in rp:
    if "The corrected Wald p-values are virtually identical" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted E2 impact R-75")
        break

# Relabel S1 as proper subsection
for p in rp:
    if p.text.startswith("S1. Sensitivity Analysis: Strategy A IPCW"):
        ok = replace_in_para(p,
            "S1. Sensitivity Analysis: Strategy A IPCW with Baseline-Only Denominator",
            "3.5.5 IPCW denominator sensitivity analysis")
        if ok:
            changes_r += 1
            print("[R] #10 Relabeled S1 → 3.5.5")
        break

# Fix "pre-specified" in S1 rationale
for p in rp:
    if "Rationale: The primary IPCW model for Strategy A" in p.text:
        replace_in_para(p, "Rationale: ", "")
        changes_r += 1
        print("[R] #8 Removed 'Rationale:' prefix from S1")
        break

# Relabel S2 as proper subsection
for p in rp:
    if p.text.startswith("S2. Assessment of Out-of-Hospital Deaths"):
        ok = replace_in_para(p,
            "S2. Assessment of Out-of-Hospital Deaths (patients.dod)",
            "3.5.6 Mortality ascertainment assessment")
        if ok:
            changes_r += 1
            print("[R] #10 Relabeled S2 → 3.5.6")
        break

# Clean up S2 rationale prefix
for p in rp:
    if "Rationale: The primary analysis defined 28-day mortality using admissions.deathtime" in p.text:
        replace_in_para(p, "Rationale: ", "")
        changes_r += 1
        print("[R] Removed 'Rationale:' prefix from S2")
        break

# Relabel S1 impact prefix
for p in rp:
    if p.text.startswith("Impact: The L0-only sensitivity analysis"):
        replace_in_para(p, "Impact: ", "")
        changes_r += 1
        print("[R] Removed 'Impact:' prefix from S1")
        break

# Relabel S2 impact prefix
for p in rp:
    if p.text.startswith("Impact assessment: A total of 138"):
        replace_in_para(p, "Impact assessment: ", "")
        changes_r += 1
        print("[R] Removed 'Impact assessment:' prefix from S2")
        break

# Delete summary header
for p in rp:
    if p.text.strip() == "Summary of Corrections and Sensitivity Analyses":
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted summary header R-83")
        break

# Delete summary text
for p in rp:
    if "Two code-level corrections (E1: bootstrap resampling; E2: interaction test statistic)" in p.text:
        delete_paragraph(p)
        changes_r += 1
        print("[R] #10 Deleted summary text R-84")
        break

rdoc.save("analysis/result.docx")
print(f"\n✓ analysis/result.docx: {changes_r} changes saved.")
print(f"\nTotal changes: methods={changes_m}, results={changes_r}")
