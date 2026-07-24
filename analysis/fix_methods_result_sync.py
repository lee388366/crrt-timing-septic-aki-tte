"""
fix_methods_result_sync.py
==========================
Four targeted fixes across methods_section.docx and result.docx:

F1 — n discrepancy: Figure 1 / 1B captions: 1,952 → 1,887
     Table 4 footnote: add explanatory note about sofa_baseline vs sofa_t0
F2 — Methods median update: creatinine 1.20 → 1.90 mg/dL; APACHE III 68 → 69
F3 — RRT-free days: remove Wilcoxon p-values from Table 3 rows + result text
F4 — Fine-Gray post-hoc note in Methods P70 + result.docx 3.3.2
"""

import re
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

METHODS_PATH = "./广西医科大学/methods_section.docx"
RESULT_PATH  = "analysis/result.docx"

# ── helpers ────────────────────────────────────────────────────────────────────
def set_run_tnr(run, pt=11):
    run.font.name = "Times New Roman"
    run.font.size = Pt(pt)
    rPr = run._r.get_or_add_rPr()
    rFonts = OxmlElement("w:rFonts")
    rFonts.set(qn("w:ascii"),    "Times New Roman")
    rFonts.set(qn("w:hAnsi"),    "Times New Roman")
    rFonts.set(qn("w:eastAsia"), "Times New Roman")
    rPr.insert(0, rFonts)

def replace_in_para(para, old, new):
    """Replace text across runs, preserving first-run formatting."""
    full = "".join(r.text for r in para.runs)
    if old not in full:
        return False
    new_full = full.replace(old, new)
    # rewrite into first run, blank the rest
    if para.runs:
        para.runs[0].text = new_full
        set_run_tnr(para.runs[0])
        for r in para.runs[1:]:
            r.text = ""
    return True

def append_run(para, text, italic=False, bold=False, pt=10):
    run = para.add_run(text)
    run.italic = italic
    run.bold   = bold
    set_run_tnr(run, pt)

changes = []

# ══════════════════════════════════════════════════════════════════════════════
# ── A. METHODS SECTION ────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
methods = Document(METHODS_PATH)

for i, p in enumerate(methods.paragraphs):
    t = p.text

    # F2a: creatinine median 1.20 → 1.90
    if "serum creatinine (median 1.20 mg/dL)" in t:
        if replace_in_para(p, "serum creatinine (median 1.20 mg/dL)",
                              "serum creatinine (median 1.90 mg/dL)"):
            changes.append(f"Methods P{i}: creatinine median 1.20 → 1.90 mg/dL")

    # F2b: APACHE III median 68 → 69
    if "APACHE III score (median 68)" in t:
        if replace_in_para(p, "APACHE III score (median 68)",
                              "APACHE III score (median 69)"):
            changes.append(f"Methods P{i}: APACHE III median 68 → 69")

    # F4a: Fine-Gray post-hoc note after P70 secondary outcome analyses
    if (i >= 69 and i <= 73 and
        "RRT-free days" in t and
        "IPCW-weighted logistic regression" in t and
        "Fine-Gray" not in t):
        # Append the Fine-Gray sentence to this paragraph
        fine_gray_sentence = (
            " As a supplementary post-hoc analysis for 28-day renal recovery, "
            "a Fine-Gray competing risks model was additionally estimated with "
            "28-day death as the competing event (Event 2), yielding a "
            "subdistribution hazard ratio (sHR; cmprsk package). "
            "This analysis was not pre-specified but was added following initial "
            "review to account for the substantial competing risk of death in "
            "this critically ill population."
        )
        append_run(p, fine_gray_sentence, pt=11)
        changes.append(f"Methods P{i}: Fine-Gray post-hoc note appended to secondary outcome analyses paragraph")

methods.save(METHODS_PATH)
print("methods_section.docx saved.\n")

# ══════════════════════════════════════════════════════════════════════════════
# ── B. RESULT.DOCX ────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
result = Document(RESULT_PATH)

# ── F1: Figure 1 & 1B captions: 1,952 → 1,887 ────────────────────────────────
for i, p in enumerate(result.paragraphs):
    t = p.text
    if ("Figure 1." in t or "Figure 1B." in t) and "1,952" in t:
        if replace_in_para(p, "1,952", "1,887"):
            changes.append(f"Result P{i}: Figure 1/1B caption n=1,952 → 1,887")

# ── F1: Table 4 footnote — add note about sofa_baseline vs sofa_t0 ────────────
for i, p in enumerate(result.paragraphs):
    if (p.text.startswith("sHR = subdistribution hazard ratio") and
            "sofa_baseline" not in p.text and "sofa_t0" not in p.text):
        footnote_add = (
            " Note: The equipoise subgroup (SOFA 8–14) in this table was defined "
            "using the sofa_baseline variable from the analytic dataset (n\u2009=\u20091,952), "
            "which captures SOFA at ICU admission. The primary outcome analysis "
            "(Table 2, Figure 1) used sofa_t0 — SOFA at the index time t\u2080 — "
            "consistent with the pre-specified methods definition (n\u2009=\u20091,887). "
            "The 65-patient difference reflects different SOFA measurement time points "
            "and does not affect the primary conclusions."
        )
        append_run(p, footnote_add, pt=10)
        changes.append(f"Result P{i}: Table 4 footnote — sofa_baseline vs sofa_t0 note added")

# ── F3: Remove Wilcoxon p-values from RRT-free days in Table 3 ───────────────
# Table 3 is the 3rd table (index 2)
tbl3 = result.tables[2]
rft_rows_cleaned = 0
for row in tbl3.rows:
    cells = row.cells
    # Identify RRT-free days rows: check if any cell contains "RRT-Free Days" heading
    # or the rows with "Wilcoxon p <0.001" in the last cell and are under RRT-free section
    # Strategy: rows whose row text has "Wilcoxon p <0.001" AND no ICU-related text
    row_text = " | ".join(c.text for c in cells)
    if "Wilcoxon p <0.001" in row_text and "ICU" not in row_text:
        # Clear the Wilcoxon p value from the last non-empty cell
        for cell in cells:
            if "Wilcoxon p <0.001" in cell.text:
                for para in cell.paragraphs:
                    if "Wilcoxon p <0.001" in para.text:
                        replace_in_para(para, "Wilcoxon p <0.001", "")
                        rft_rows_cleaned += 1
changes.append(f"Result Table 3: Removed Wilcoxon p <0.001 from {rft_rows_cleaned} RRT-free days cell(s)")

# ── F3: Remove Wilcoxon p-values from result text paragraphs (RRT-free section) ─
for i, p in enumerate(result.paragraphs):
    t = p.text
    # Target the RRT-free days section (3.3.2) that mentions Wilcoxon p<0.001
    if "RRT-free days" in t and "Wilcoxon p" in t and "ICU" not in t:
        cleaned = re.sub(r";?\s*Wilcoxon\s+p\s*[<>]?0?\.?001\b", "", t)
        # Also clean "(Wilcoxon p <0.001)" patterns
        cleaned = re.sub(r"\(Wilcoxon\s+p\s*[<>]?\s*0?\.?001\)", "", cleaned)
        cleaned = cleaned.strip()
        if cleaned != t:
            replace_in_para(p, t, cleaned)
            changes.append(f"Result P{i}: Removed Wilcoxon p mention from RRT-free days text")

# ── F4: Fine-Gray post-hoc note in section 3.3.2 ─────────────────────────────
for i, p in enumerate(result.paragraphs):
    t = p.text
    if ("3.3.2" in t or ("RRT-free days" in t and "3.3" in " ".join(
            [result.paragraphs[max(0,i-3)].text, result.paragraphs[max(0,i-2)].text]))):
        pass  # handled below

# Find the paragraph that starts with "3.3.2"
for i, p in enumerate(result.paragraphs):
    if p.text.startswith("3.3.2") and "RRT-free" in p.text:
        if "post-hoc" not in p.text:
            note = (" The Fine-Gray competing risks analysis reported in Section 3.3.2 "
                    "and Table 4 was a post-hoc supplementary analysis added to account "
                    "for the substantial competing risk of death; it was not pre-specified "
                    "in the original statistical analysis plan.")
            append_run(p, note, pt=11)
            changes.append(f"Result P{i}: Fine-Gray post-hoc note added to 3.3.2 paragraph")
        break

# ── F4: Add Fine-Gray note near the sHR paragraph (Section 3.3.2 body) ───────
for i, p in enumerate(result.paragraphs):
    t = p.text
    # The paragraph that introduces Fine-Gray results
    if ("Fine-Gray" in t and "competing risks" in t and
            "post-hoc" not in t and "not pre-specified" not in t and
            "sHR" in t):
        # Prepend a note at the start about post-hoc nature
        note_prefix = ("[Post-hoc supplementary analysis] ")
        if note_prefix not in t:
            # Insert note at beginning by replacing the paragraph content
            new_text = note_prefix + t
            replace_in_para(p, t, new_text)
            changes.append(f"Result P{i}: [Post-hoc] prefix added to Fine-Gray results paragraph")
        break

result.save(RESULT_PATH)
print("result.docx saved.\n")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"Total changes: {len(changes)}")
for c in changes:
    print(f"  • {c}")
