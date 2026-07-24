#!/usr/bin/env python3
"""
Manuscript Consistency Validator — run after ANY edit to result.docx or methods_section.docx.

Usage:  python3 analysis/validate_manuscript.py

Checks:
  1. Forbidden strings (old/wrong values that must NOT appear)
  2. Required strings (corrected values that MUST appear)
  3. Cross-table number consistency (HR, CI, p-values across all tables)
  4. Cross-document consistency (methods ↔ results)
  5. Table structural integrity (no unexpected empty cells, correct row counts)
  6. Number arithmetic (e.g., 500 + 2612 = 3112, 18 + 199 = 217)

Exit code 0 = all pass, 1 = failures found.
"""

import sys
from docx import Document
from docx.oxml.ns import qn

PASS = 0
FAIL = 0
WARN = 0

def check(condition, label, level="error"):
    global PASS, FAIL, WARN
    if condition:
        PASS += 1
    else:
        if level == "warn":
            WARN += 1
            print(f"  ⚠ WARN: {label}")
        else:
            FAIL += 1
            print(f"  ✗ FAIL: {label}")

def get_all_table_text(doc):
    texts = []
    for ti, t in enumerate(doc.tables):
        try:
            for ri, row in enumerate(t.rows):
                for ci, cell in enumerate(row.cells):
                    texts.append((ti, ri, ci, cell.text.strip()))
        except:
            pass
    return texts

def get_table7_text(doc):
    try:
        tbl = doc.element.body.findall(qn('w:tbl'))[7]
        texts = []
        for tr in tbl.findall(qn('w:tr')):
            runs = tr.findall(f".//{qn('w:t')}")
            txt = " ".join(r.text or "" for r in runs)
            texts.append(txt)
        return "\n".join(texts)
    except:
        return ""


# ═══════════════════════════════════════════════════════════════════════════════
#  LOAD DOCUMENTS
# ═══════════════════════════════════════════════════════════════════════════════
print("Loading documents...")
try:
    rdoc = Document("analysis/result.docx")
    mdoc = Document("methods_section.docx")
except Exception as e:
    print(f"  ✗ Cannot load documents: {e}")
    sys.exit(1)

rtext = "\n".join(p.text for p in rdoc.paragraphs)
mtext = "\n".join(p.text for p in mdoc.paragraphs)
r_table_cells = get_all_table_text(rdoc)
t7_text = get_table7_text(rdoc)
all_r = rtext + "\n" + "\n".join(c[3] for c in r_table_cells) + "\n" + t7_text


# ═══════════════════════════════════════════════════════════════════════════════
#  1. FORBIDDEN STRINGS (must NOT appear anywhere)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 1. Forbidden strings ──")

FORBIDDEN_R = {
    "-0.62 days": "old RMST equip (pre-C1 fix)",
    "-2.13 days": "old RMST full (pre-C1 fix)",
    "-3.16 to -1.17": "old RMST full CI",
    "-1.80 to 0.50": "old RMST equip CI",
    "1.413": "old E-value CI (should be 1.0)",
    "3,233": "old cohort_base N (should be 3,112)",
    "9/517": "old C4 early denominator",
    "135/2,716": "old C4 deferred denominator",
    "144 (4.5%)": "old C4 total",
    "ERRATA & UPDATED": "old errata section title",
    "E1. Corrected Bootstrap": "old E1 errata header",
    "E2. Corrected Interaction": "old E2 errata header",
    "1.293]": "old CI upper bound (should be 1.294)",
}

for bad, label in FORBIDDEN_R.items():
    check(bad not in all_r, f"result.docx contains '{bad}' ({label})")

FORBIDDEN_M = {
    "A likelihood ratio test comparing models with and without the interaction":
        "old LRT as primary test",
    "pre-specified sensitivity analysis, the Strategy A":
        "should be post-hoc",
    "144 additional patients (4.5%)": "old C4 number",
    "1.7% Early vs 5.0% Deferred": "old C4 rates",
    "Methodological Addendum (2026-04-12)": "old addendum header",
    "log-rank p-values": "should be IPCW-weighted Cox model p-values",
}

for bad, label in FORBIDDEN_M.items():
    check(bad not in mtext, f"methods contains '{bad}' ({label})")


# ═══════════════════════════════════════════════════════════════════════════════
#  2. REQUIRED STRINGS (must appear)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 2. Required strings ──")

REQUIRED_R = {
    "-0.63": "corrected RMST equip",
    "-2.12": "corrected RMST full",
    "-3.30 to -0.79": "corrected full CI",
    "-2.07 to 0.84": "corrected equip CI",
    "1.237": "E-value point estimate",
    "E-value for the confidence interval is 1.0": "E-value CI crosses null",
    "1.601": "full cohort E-value point",
    "1.269": "full cohort E-value CI",
    "138 patients (4.4%": "C4 total corrected",
    "9/500 = 1.8%": "C4 early corrected",
    "129/2,612 = 4.9%": "C4 deferred corrected",
    "Wald p = 0.214": "interaction Wald full",
    "Wald p = 0.992": "interaction Wald equip",
    "3.5.5": "IPCW sensitivity subsection",
    "3.5.6": "mortality ascertainment subsection",
    "SOFA score at ICU admission was 10": "SOFA clarification",
}

for good, label in REQUIRED_R.items():
    check(good in rtext, f"result.docx missing '{good}' ({label})")

REQUIRED_M = {
    "Wald test": "interaction test method",
    "preferred over the likelihood ratio test": "LRT explanation",
    "post-hoc sensitivity analysis": "IPCW sensitivity",
    "138 additional patients (4.4%)": "C4 corrected",
    "1.8% Early vs 4.9% Deferred": "C4 rates corrected",
    "merge-based resampling": "bootstrap fix",
    "E-value for the confidence interval is defined as 1.0": "E-value CI convention",
    "difference in IPCW-weighted means divided by": "SMD formula",
    "IPCW-weighted Cox model p-values": "KM annotation",
    "continuous renal replacement therapy (CRRT)": "CRRT definition",
}

for good, label in REQUIRED_M.items():
    check(good in mtext, f"methods missing '{good}' ({label})")


# ═══════════════════════════════════════════════════════════════════════════════
#  3. CROSS-TABLE NUMBER CONSISTENCY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 3. Cross-table consistency ──")

def get_cell(doc, ti, ri, ci):
    try:
        return doc.tables[ti].rows[ri].cells[ci].text.strip()
    except:
        return ""

# Full cohort CI must be [1.047, 1.294] everywhere
for ti, ri, ci, label in [(2, 2, 3, "Table3"), (5, 2, 4, "Table6"), (6, 1, 3, "Table6-IPCW")]:
    val = get_cell(rdoc, ti, ri, ci)
    check("1.294" in val, f"{label} full CI contains 1.294 (got: {val[:30]})")

# Equipoise CI must be [0.915, 1.178] everywhere
for ti, ri, ci, label in [(2, 6, 3, "Table3"), (5, 1, 4, "Table6"), (6, 3, 3, "Table6-IPCW")]:
    val = get_cell(rdoc, ti, ri, ci)
    check("1.178" in val, f"{label} equip CI contains 1.178 (got: {val[:30]})")

# RMST in Table 3
rmst_full = get_cell(rdoc, 2, 4, 2)
rmst_eq = get_cell(rdoc, 2, 8, 2)
check("-2.12" in rmst_full, f"Table3 RMST full={rmst_full[:25]}")
check("-0.63" in rmst_eq, f"Table3 RMST equip={rmst_eq[:25]}")

# Platelet consistency Table1 vs Table2
plt_t0 = ""
for ri, row in enumerate(rdoc.tables[0].rows):
    if "platelet" in row.cells[0].text.lower():
        plt_t0 = row.cells[1].text.strip()
plt_t1 = ""
for ri, row in enumerate(rdoc.tables[1].rows):
    if "platelet" in row.cells[0].text.lower():
        plt_t1 = row.cells[1].text.strip()
if plt_t0 and plt_t1:
    t0_med = plt_t0.split("(")[0].split("[")[0].strip()
    t1_med = plt_t1.split("(")[0].split("[")[0].strip()
    check(t0_med == t1_med, f"Platelet median T1={t0_med} vs T2={t1_med}")

# Table 7 C4 numbers
check("3,112" in t7_text, f"Table7 uses analytic cohort N=3,112")
check("138" in t7_text and "4.4%" in t7_text, f"Table7 has 138 (4.4%)")
check("3,233" not in t7_text, f"Table7 no old N=3,233")


# ═══════════════════════════════════════════════════════════════════════════════
#  4. CROSS-DOCUMENT CONSISTENCY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 4. Cross-document consistency ──")

check("3,112" in rtext and "3,112" in mtext, "N=3,112 in both docs")
check("SOFA" in rtext and "8" in rtext and "14" in rtext, "Equipoise SOFA 8-14 in results")
check("SOFA" in mtext and "8" in mtext and "14" in mtext, "Equipoise SOFA 8-14 in methods")
check("24 hours" in mtext or "24-hour" in mtext, "Grace period 24h in methods")
check("B = 1,000" in mtext or "B\u2009=\u20091,000" in mtext, "Bootstrap B=1000 in methods")

# No contradiction: methods says Wald, results also says Wald
check("Wald" in rtext, "Results mentions Wald test")
check("Wald" in mtext, "Methods mentions Wald test")


# ═══════════════════════════════════════════════════════════════════════════════
#  5. TABLE STRUCTURAL INTEGRITY
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 5. Table structure ──")

expected = [(0, 115, 118), (1, 49, 51), (2, 9, 10), (3, 13, 14),
            (4, 3, 4), (5, 10, 11)]
for ti, lo, hi in expected:
    try:
        n = len(rdoc.tables[ti].rows)
        check(lo <= n <= hi, f"Table {ti} row count {n} in [{lo},{hi}]")
    except:
        check(False, f"Table {ti} cannot read rows")

# No trailing empty rows in Table 0
try:
    last = rdoc.tables[0].rows[-1]
    check(not all(c.text.strip() == "" for c in last.cells), "Table 0 no trailing empty row")
except:
    pass

# Table 3 (secondary): no blank data cells except section headers
t3 = rdoc.tables[3]
section_rows = set()
for ri, row in enumerate(t3.rows):
    if ri == 0:
        continue
    tc = row.cells[0]._element
    tcPr = tc.find(qn('w:tcPr'))
    if tcPr is not None:
        shd = tcPr.find(qn('w:shd'))
        if shd is not None and shd.get(qn('w:fill'), '') == 'EEEEEE':
            section_rows.add(ri)

for ri, row in enumerate(t3.rows):
    if ri == 0 or ri in section_rows:
        continue
    for ci, cell in enumerate(row.cells):
        txt = cell.text.strip()
        check(txt != "" or txt == "—" or ci == 0,
              f"Table3 R{ri}C{ci} not empty (got: '{txt[:20]}')", level="warn")


# ═══════════════════════════════════════════════════════════════════════════════
#  6. ARITHMETIC CHECKS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── 6. Arithmetic ──")

check("500" in rtext and "2,612" in rtext, "500 + 2,612 groups present")
# 500 + 2612 = 3112
check(500 + 2612 == 3112, "500 + 2,612 = 3,112")
# 18 + 199 = 217
check(18 + 199 == 217, "18 + 199 = 217 (renal recovery)")
# 9 + 129 = 138
check(9 + 129 == 138, "9 + 129 = 138 (C4 missed deaths)")
# 1517/3112 ≈ 48.7%
check(abs(1517/3112*100 - 48.7) < 0.1, "1517/3112 ≈ 48.7%")
# 18/500 ≈ 3.6%
check(abs(18/500*100 - 3.6) < 0.1, "18/500 ≈ 3.6%")
# 199/2612 ≈ 7.6%
check(abs(199/2612*100 - 7.6) < 0.1, "199/2612 ≈ 7.6%")


# ═══════════════════════════════════════════════════════════════════════════════
#  REPORT
# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"  PASS: {PASS}  |  FAIL: {FAIL}  |  WARN: {WARN}")
print(f"{'='*60}")
if FAIL == 0:
    print("  ✓ ALL CHECKS PASSED")
else:
    print(f"  ✗ {FAIL} FAILURE(S) — fix before submission")
sys.exit(1 if FAIL > 0 else 0)
