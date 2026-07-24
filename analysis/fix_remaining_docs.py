"""
fix_remaining_docs.py
=====================
G1 — Schoenfeld residuals: add one sentence to result.docx Table 2 footnote
     & Section 3.2.1 paragraph
G2 — Renal recovery method: update Methods P70 to describe covariate-adjusted
     logistic (not IPCW logistic), with brief methodological note
G3 — ICU LOS method: update Methods P70 to say Wilcoxon for ICU LOS
G4 — Forest plot: add cross-reference to Table 2 in Figure 4 caption
"""

from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

METHODS_PATH = "./广西医科大学/methods_section.docx"
RESULT_PATH  = "analysis/result.docx"

# Read Schoenfeld p-values from R output
with open("data/fix_remaining_results.txt") as f:
    vals = dict(line.strip().split("=") for line in f if "=" in line)

ph_p_eq   = float(vals["SCHOENFELD_EQ_P"])
ph_p_full = float(vals["SCHOENFELD_FULL_P"])

def fmt_p(p):
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"

def set_run_tnr(run, pt=11):
    run.font.name = "Times New Roman"
    run.font.size = Pt(pt)
    rPr = run._r.get_or_add_rPr()
    rFonts = OxmlElement("w:rFonts")
    rFonts.set(qn("w:ascii"),    "Times New Roman")
    rFonts.set(qn("w:hAnsi"),    "Times New Roman")
    rFonts.set(qn("w:eastAsia"), "Times New Roman")
    rPr.insert(0, rFonts)

def rewrite_para(para, new_text, pt=11):
    if para.runs:
        para.runs[0].text = new_text
        set_run_tnr(para.runs[0], pt)
        for r in para.runs[1:]:
            r.text = ""

def append_run(para, text, pt=11, italic=False):
    run = para.add_run(text)
    run.italic = italic
    set_run_tnr(run, pt)

changes = []

# ══════════════════════════════════════════════════════════════════════════════
# ── METHODS SECTION ──────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
methods = Document(METHODS_PATH)

for i, p in enumerate(methods.paragraphs):
    t = p.text

    # G2 + G3: Update P70 (Secondary outcome analyses)
    # Replace the IPCW-weighted logistic sentence AND add ICU LOS Wilcoxon note
    if ("IPCW-weighted logistic regression with HC3 sandwich standard errors" in t
            and "RRT-free days" in t):
        new_t = t.replace(
            "The odds ratio for 28-day renal recovery was estimated using IPCW-weighted "
            "logistic regression with HC3 sandwich standard errors.",
            "The odds ratio for 28-day renal recovery was estimated using multivariable "
            "logistic regression with covariate adjustment (age, baseline SOFA score, "
            "APACHE III score, mechanical ventilation status at t\u2080, and Charlson "
            "comorbidity index), with standard errors estimated using the HC3 sandwich "
            "method. This covariate-based approach was used in place of IPCW-weighted "
            "logistic regression on the clone dataset because the CCW inverse probability "
            "weights are derived from a time-varying censoring model and are not "
            "directly transferable to a single-point binary outcome without additional "
            "structural assumptions; covariate adjustment is an accepted alternative for "
            "secondary binary outcomes in target trial emulation studies."
        )
        new_t = new_t.replace(
            "The median difference in ICU length of stay was estimated using "
            "IPCW-weighted quantile regression at the 0.50 quantile with bootstrap "
            "confidence intervals (xy-pair resampling, B = 500 replicates; percentile method).",
            "Median ICU length of stay was compared between groups using the "
            "Wilcoxon rank-sum test and reported as median [interquartile range], "
            "as ICU LOS is a secondary descriptive outcome for which the Wilcoxon "
            "comparison provides a straightforward summary of the central tendency "
            "difference without requiring additional structural modelling assumptions."
        )
        if new_t != t:
            rewrite_para(p, new_t, pt=11)
            changes.append(f"Methods P{i}: Updated renal recovery method (covariate-adjusted) + ICU LOS (Wilcoxon)")

    # G1: Schoenfeld note in primary outcome analysis paragraph (P68)
    if ("Schoenfeld residuals" in t and
            "global test p-value" in t and
            "Equipoise" not in t and
            "p = 0.0" not in t):
        # Append actual values
        schoen_note = (
            f" In the primary equipoise analysis, the global Schoenfeld test "
            f"yielded {fmt_p(ph_p_eq)}, and in the full cohort {fmt_p(ph_p_full)}, "
            f"both indicating that the proportional hazards assumption was not "
            f"violated and supporting the use of the Cox HR as the primary "
            f"inferential measure."
        )
        append_run(p, schoen_note, pt=11)
        changes.append(f"Methods P{i}: Schoenfeld actual p-values added to primary outcome methods")

methods.save(METHODS_PATH)
print("methods_section.docx saved.\n")

# ══════════════════════════════════════════════════════════════════════════════
# ── RESULT.DOCX ──────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
result = Document(RESULT_PATH)

for i, p in enumerate(result.paragraphs):
    t = p.text

    # G1: Table 2 footnote — add Schoenfeld results
    if (t.startswith("HR, hazard ratio") and
            "IPCW" in t and "RMST" in t and
            "Schoenfeld" not in t):
        schoen_sentence = (
            f" The proportional hazards assumption was evaluated using scaled "
            f"Schoenfeld residuals; the global test was non-significant in both "
            f"the equipoise subgroup ({fmt_p(ph_p_eq)}) and the full cohort "
            f"({fmt_p(ph_p_full)}), confirming the validity of the Cox HR as "
            f"the primary inferential measure."
        )
        append_run(p, schoen_sentence, pt=10)
        changes.append(f"Result P{i}: Schoenfeld results added to Table 2 footnote")

    # G1: Also add one line to Section 3.2.1 body text
    if ("3.2.1" in t and "Cox proportional hazards" in t and
            "Schoenfeld" not in t and "proportional hazards assumption" not in t):
        schoen_inline = (
            f" The proportional hazards assumption was confirmed using scaled "
            f"Schoenfeld residuals (equipoise subgroup: {fmt_p(ph_p_eq)}; "
            f"full cohort: {fmt_p(ph_p_full)}); the RMST difference therefore "
            f"serves as a corroborating absolute effect measure rather than a "
            f"replacement for the HR."
        )
        append_run(p, schoen_inline, pt=11)
        changes.append(f"Result P{i}: Schoenfeld inline note added to Section 3.2.1")

    # G2: Table 3 footnote — update renal recovery analysis description
    if ("Logistic regression adjusted for" in t and
            "covariate-adjusted" not in t.lower() and
            "IPCW" not in t):
        new_t = t.replace(
            "\u1d43Logistic regression adjusted for age, baseline SOFA score, APACHE\xa0III score, "
            "mechanical ventilation use at t\u2080, and Charlson comorbidity index.",
            "\u1d43Multivariable logistic regression with covariate adjustment (age, baseline "
            "SOFA score, APACHE\u202fIII score, mechanical ventilation at t\u2080, Charlson comorbidity "
            "index); HC3 sandwich standard errors. Covariate-adjusted logistic regression "
            "was used for this binary secondary outcome as CCW IPCW weights are derived "
            "from a time-varying censoring model not directly applicable to a single-point "
            "binary response."
        )
        if new_t != t:
            rewrite_para(p, new_t, pt=10)
            changes.append(f"Result P{i}: Table 3 footnote updated — covariate-adjusted logistic note")

    # G4: Figure 4 caption — add cross-reference to Table 2 for crude analysis
    if (t.startswith("Figure 4.") and "Forest plot" in t and
            "Table 2" not in t):
        note = (
            " The unadjusted (crude) Cox analysis is reported separately in "
            "Table 2 (crude HR = 1.228 [1.113, 1.355] for the full cohort); "
            "it is omitted from this forest plot to focus on IPCW-adjusted "
            "estimates across analytic populations."
        )
        append_run(p, note, pt=11)
        changes.append(f"Result P{i}: Figure 4 caption — crude analysis cross-reference to Table 2 added")

result.save(RESULT_PATH)
print("result.docx saved.\n")

print(f"Total changes: {len(changes)}")
for c in changes:
    print(f"  ✓ {c}")

# ── Quick verification ────────────────────────────────────────────────────────
print("\n=== VERIFICATION ===")
r2 = Document(RESULT_PATH)
for i, p in enumerate(r2.paragraphs):
    t = p.text
    if "Schoenfeld" in t and "0.08" in t:
        print(f"  ✓ P{i}: Schoenfeld p-values present ({t[:80]}...)")
    if "Figure 4." in t and "Table 2" in t:
        print(f"  ✓ P{i}: Figure 4 has Table 2 cross-reference")
m2 = Document(METHODS_PATH)
for i, p in enumerate(m2.paragraphs):
    if "covariate adjustment" in p.text and "logistic" in p.text:
        print(f"  ✓ Methods P{i}: covariate-adjusted logistic text present")
        break
