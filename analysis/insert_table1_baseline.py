"""
insert_table1_baseline.py
─────────────────────────
Fixes and inserts table1_baseline.docx into result.docx as the new Table 1.
- Fixes: GCS label, Colloid display, footnote merge, adds table caption
- Inserts as Table 1 before the existing Table 1
- Renames existing Tables 1→2, 2→3, 3→4, 4→5, 5→6
- Updates all Table references in body text
"""

import copy
import re
from docx import Document
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

BASE   = "./广西医科大学"
SRC    = f"{BASE}/table1_baseline.docx"
RESULT = f"{BASE}/analysis/result.docx"

# ─────────────────────────────────────────────────────────────────────────────
# Helper: copy run formatting from a template run
# ─────────────────────────────────────────────────────────────────────────────
def copy_run_fmt(run, bold=False, font_name="Times New Roman", font_size=9):
    run.bold = bold
    run.font.name = font_name
    run.font.size = Pt(font_size)


def set_cell_font(cell, font_name="Times New Roman", font_size=9, bold=False, italic=False):
    for para in cell.paragraphs:
        for run in para.runs:
            run.font.name = font_name
            run.font.size = Pt(font_size)
            run.bold = bold
            run.italic = italic


# ─────────────────────────────────────────────────────────────────────────────
# Step 1: Load & fix table1_baseline.docx content
# ─────────────────────────────────────────────────────────────────────────────
src = Document(SRC)
src_table = src.tables[0]

# Fix GCS label: "GCS total score" → "GCS score (minimum)"
for row in src_table.rows:
    cell0 = row.cells[0]
    for para in cell0.paragraphs:
        for run in para.runs:
            if "GCS total score" in run.text:
                run.text = run.text.replace("GCS total score", "GCS score (minimum)")
        if "GCS total score" in para.text:
            full_text = para.text
            new_text = full_text.replace("GCS total score", "GCS score (minimum)")
            for run in para.runs:
                run.text = ""
            if para.runs:
                para.runs[0].text = new_text

# Fix Colloid: "1 (0.0%)" → "1 (<0.1%)"
for row in src_table.rows:
    cell0_text = row.cells[0].text.strip()
    cell1 = row.cells[1]
    if "Artificial colloid" in cell0_text:
        for para in cell1.paragraphs:
            for run in para.runs:
                if "0.0%" in run.text:
                    run.text = run.text.replace("0.0%", "<0.1%")
            if "0.0%" in para.text:
                for run in para.runs:
                    run.text = run.text.replace("0.0%", "<0.1%")

# Fix footnote: clear cell 1 of the footnote row (cell 0 has the full footnote already)
footnote_row = src_table.rows[116]
# Clear cell 1 to avoid duplication
for para in footnote_row.cells[1].paragraphs:
    for run in para.runs:
        run.text = ""
# Also merge cells in the footnote row using gridSpan
tc0 = footnote_row.cells[0]._tc
# Set gridSpan = 2 on cell 0 so it spans both columns
tcPr0 = tc0.find(qn('w:tcPr'))
if tcPr0 is None:
    tcPr0 = OxmlElement('w:tcPr')
    tc0.insert(0, tcPr0)
gridSpan = tcPr0.find(qn('w:gridSpan'))
if gridSpan is None:
    gridSpan = OxmlElement('w:gridSpan')
    tcPr0.append(gridSpan)
gridSpan.set(qn('w:val'), '2')

# Mark cell 1 as continuation of merge (vMerge or just make it empty with no gridSpan)
tc1 = footnote_row.cells[1]._tc
# Remove cell 1 from the row entirely if possible, or mark as merged continuation
# Simplest: set its content to empty and shrink it
tcPr1 = tc1.find(qn('w:tcPr'))
if tcPr1 is None:
    tcPr1 = OxmlElement('w:tcPr')
    tc1.insert(0, tcPr1)

print("[OK] Content fixes applied to src table")

# ─────────────────────────────────────────────────────────────────────────────
# Step 2: Load result.docx and find insertion point
# ─────────────────────────────────────────────────────────────────────────────
result = Document(RESULT)

# Find paragraph index of "Table 1. Baseline Characteristics..."
target_para_idx = None
for i, p in enumerate(result.paragraphs):
    if p.text.strip().startswith("Table 1."):
        target_para_idx = i
        break

print(f"[OK] Target insertion point: before paragraph index {target_para_idx}")
print(f"     Text: {result.paragraphs[target_para_idx].text[:60]}")

# ─────────────────────────────────────────────────────────────────────────────
# Step 3: Rename existing Tables 1→2, 2→3, 3→4, 4→5, 5→6 in paragraph titles
# ─────────────────────────────────────────────────────────────────────────────
def renumber_table_refs_in_para(para, direction="bump"):
    """Replace 'Table N' with 'Table N+1' in paragraph text (preserving runs)."""
    full_text = para.text
    # Check if this paragraph contains a Table reference to renumber
    if not re.search(r'\bTable\s+[1-5]\b', full_text):
        return False
    # Build new text by replacing Table N → Table N+1 (in reverse order to avoid double-replacement)
    new_text = re.sub(r'\bTable\s+5\b', 'TABLE_6_TEMP', full_text)
    new_text = re.sub(r'\bTable\s+4\b', 'TABLE_5_TEMP', new_text)
    new_text = re.sub(r'\bTable\s+3\b', 'TABLE_4_TEMP', new_text)
    new_text = re.sub(r'\bTable\s+2\b', 'TABLE_3_TEMP', new_text)
    new_text = re.sub(r'\bTable\s+1\b', 'TABLE_2_TEMP', new_text)
    # Restore TEMP markers
    new_text = new_text.replace('TABLE_6_TEMP', 'Table 6')
    new_text = new_text.replace('TABLE_5_TEMP', 'Table 5')
    new_text = new_text.replace('TABLE_4_TEMP', 'Table 4')
    new_text = new_text.replace('TABLE_3_TEMP', 'Table 3')
    new_text = new_text.replace('TABLE_2_TEMP', 'Table 2')

    if new_text == full_text:
        return False

    # Rewrite using existing runs (preserve formatting)
    if para.runs:
        # Clear all runs then put into first run
        first_run = para.runs[0]
        first_run_fmt = (first_run.bold, first_run.font.name, first_run.font.size)
        for run in para.runs:
            run.text = ""
        first_run.text = new_text
        first_run.bold = first_run_fmt[0]
    return True

changes = []

# Update all paragraph text references
for i, p in enumerate(result.paragraphs):
    t = p.text
    if re.search(r'\bTable\s+[1-5]\b', t):
        old_text = t
        if renumber_table_refs_in_para(p):
            changes.append(f"P{i:03d}: '{old_text[:80]}' → '{p.text[:80]}'")

print(f"\n[OK] Renumbered {len(changes)} paragraph(s):")
for c in changes:
    print(f"  {c}")

# ─────────────────────────────────────────────────────────────────────────────
# Step 4: Insert new Table 1 caption paragraph + table XML before target para
# ─────────────────────────────────────────────────────────────────────────────
target_para_elem = result.paragraphs[target_para_idx]._element

# 4a: Insert blank paragraph after the new table (before old Table 2 title)
blank_para_elem = OxmlElement('w:p')
target_para_elem.addprevious(blank_para_elem)

# 4b: Insert the footnote text as a plain paragraph
footnote_text = (
    "\u00b9Continuous variables: median (IQR); categorical variables: n (%). "
    "\u00b2Vasopressor use is 100% by design (inclusion criterion). "
    "\u00b3Baseline SCr derived from mimiciv_derived.creatinine_baseline "
    "(lowest SCr within the 7 days preceding ICU admission). "
    "ESRD patients are excluded (ICD-10: N185/N186; ICD-9: 5856/V4511/V4512); "
    "CKD here represents non-ESRD chronic kidney disease. "
    "Data source: MIMIC-IV v3.1 (Beth Israel Deaconess Medical Center, 2008\u20132019; dates shifted to 2110\u20132210)."
)

fn_para_elem = OxmlElement('w:p')
fn_run_elem  = OxmlElement('w:r')
fn_rpr_elem  = OxmlElement('w:rPr')
fn_rFonts   = OxmlElement('w:rFonts')
fn_rFonts.set(qn('w:ascii'), 'Times New Roman')
fn_rFonts.set(qn('w:hAnsi'), 'Times New Roman')
fn_sz_elem   = OxmlElement('w:sz')
fn_sz_elem.set(qn('w:val'), '18')  # 9pt = 18 half-points
fn_szCs_elem = OxmlElement('w:szCs')
fn_szCs_elem.set(qn('w:val'), '18')
fn_rpr_elem.append(fn_rFonts)
fn_rpr_elem.append(fn_sz_elem)
fn_rpr_elem.append(fn_szCs_elem)
fn_run_elem.append(fn_rpr_elem)
fn_t_elem = OxmlElement('w:t')
fn_t_elem.text = footnote_text
fn_run_elem.append(fn_t_elem)
fn_para_elem.append(fn_run_elem)
blank_para_elem.addprevious(fn_para_elem)

# 4c: Deep copy the table XML from src and insert before the footnote paragraph
src_tbl_elem = copy.deepcopy(src_table._tbl)
fn_para_elem.addprevious(src_tbl_elem)

# 4d: Insert the new Table 1 caption paragraph before the table
caption_text = "Table 1.  Baseline Characteristics of the Study Cohort (Overall, N\u2009=\u20093,112)"
cap_para_elem = OxmlElement('w:p')
cap_pPr_elem  = OxmlElement('w:pPr')
cap_jc        = OxmlElement('w:jc')
cap_jc.set(qn('w:val'), 'left')
cap_pPr_elem.append(cap_jc)
cap_para_elem.append(cap_pPr_elem)

cap_run_elem  = OxmlElement('w:r')
cap_rpr_elem  = OxmlElement('w:rPr')
cap_rFonts    = OxmlElement('w:rFonts')
cap_rFonts.set(qn('w:ascii'), 'Times New Roman')
cap_rFonts.set(qn('w:hAnsi'), 'Times New Roman')
cap_rpr_elem.append(cap_rFonts)
cap_bold      = OxmlElement('w:b')
cap_rpr_elem.append(cap_bold)
cap_sz        = OxmlElement('w:sz')
cap_sz.set(qn('w:val'), '22')  # 11pt
cap_szCs      = OxmlElement('w:szCs')
cap_szCs.set(qn('w:val'), '22')
cap_rpr_elem.append(cap_sz)
cap_rpr_elem.append(cap_szCs)
cap_run_elem.append(cap_rpr_elem)
cap_t_elem    = OxmlElement('w:t')
cap_t_elem.text = caption_text
cap_run_elem.append(cap_t_elem)
cap_para_elem.append(cap_run_elem)
src_tbl_elem.addprevious(cap_para_elem)

print("\n[OK] New Table 1 inserted into result.docx")

# ─────────────────────────────────────────────────────────────────────────────
# Step 5: Update P002 text (cohort description) to reference both Table 1 and 2
# ─────────────────────────────────────────────────────────────────────────────
# P002 likely ends with "... are shown in Table 1" or similar
for i, p in enumerate(result.paragraphs):
    t = p.text
    if "3,112 ICU admissions" in t and "Table 2" in t:
        # Already updated
        break
    if "3,112 ICU admissions" in t and ("Table 1" not in t and "Table 2" not in t):
        # Add reference
        if p.runs:
            last_run = p.runs[-1]
            if not last_run.text.endswith('.'):
                last_run.text += " Baseline characteristics are summarized in Table 1 (overall cohort) and Table 2 (by treatment group)."
        break

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
result.save(RESULT)
print("\n[DONE] result.docx saved successfully.")
print("       New table structure: Table 1 (overall), Table 2 (stratified), Table 3-6 (outcomes/sensitivity)")
