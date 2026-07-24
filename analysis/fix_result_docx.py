"""
fix_result_docx.py
==================
One-pass fix for result.docx:
  Fix 1 — P13: replace [Figure 1 bracket placeholder] with formal caption
  Fix 2 — P14: update [Figure S1 ...] text to proper Suppl reference (.pdf→.png)
  Fix 3 — P16: update [Supplementary Figure S1b ...] .pdf → .png
  Fix 4 — P18: move Suppl. Figure S1b description to Supplementary Materials block
  Fix 5 — P26: move Figure 2 caption from inside 3.3.2 to after Table 3 (3.3.1)
  Fix 6 — P36-P40: move entire Supplementary Materials section to end of document
"""

from docx import Document
from docx.shared import Pt, RGBColor
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from lxml import etree
import copy

DOC_PATH = "analysis/result.docx"

doc = Document(DOC_PATH)
body = doc.element.body

# ── helpers ───────────────────────────────────────────────────────────────────
def para_text(p):
    return "".join(r.text or "" for r in p.iter(qn("w:r")))

def set_para_text(para, text, bold=False, font_size_pt=11):
    """Replace all runs in a paragraph with a single formatted run."""
    for r in list(para.runs):
        r._r.getparent().remove(r._r)
    run = para.add_run(text)
    run.bold = bold
    run.font.name = "Times New Roman"
    run.font.size = Pt(font_size_pt)
    # Ensure font is applied via rPr
    rPr = run._r.get_or_add_rPr()
    rFonts = OxmlElement('w:rFonts')
    rFonts.set(qn('w:ascii'),    'Times New Roman')
    rFonts.set(qn('w:hAnsi'),   'Times New Roman')
    rFonts.set(qn('w:eastAsia'),'Times New Roman')
    rPr.insert(0, rFonts)

# Collect all <w:p> elements (paragraphs) from the body
# (There may also be <w:tbl> elements for tables — we must keep those in place)
all_children = list(body)

# Build a lookup: paragraph object → index in body children
para_to_body_idx = {}
para_idx = 0
for ci, child in enumerate(all_children):
    if child.tag == qn("w:p"):
        if para_idx < len(doc.paragraphs):
            para_to_body_idx[para_idx] = ci
        para_idx += 1

print("Total paragraphs:", len(doc.paragraphs))
print("Total body children:", len(all_children))

# Print current paragraph map
for pi, para in enumerate(doc.paragraphs):
    print(f"  P{pi} (body[{para_to_body_idx.get(pi,'?')}]): {repr(para.text[:80])}")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 1 — P13: replace placeholder with formal Figure 1 caption
# ─────────────────────────────────────────────────────────────────────────────
p13 = doc.paragraphs[13]
assert "[Figure 1" in p13.text, f"Unexpected P13: {p13.text}"

fig1_text = (
    "Figure 2.  Weighted Kaplan\u2013Meier survival curves for 28-day all-cause mortality "
    "by treatment group, equipoise subgroup (SOFA 8\u201314; n\u2009=\u20091,952). "
    "IPCW-weighted curves are displayed with 95% confidence bands. "
    "Number at risk is shown below the x-axis at 0, 7, 14, 21, and 28 days. "
    "The p-value shown is derived from the IPCW-weighted Cox proportional hazards model "
    "(p\u2009=\u20090.565). "
    "See: analysis/figures/Figure2_KM_Equipoise.png."
)
set_para_text(p13, fig1_text)
print("Fix 1 done: P13 Figure 1 caption updated.")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 2 — P14: update Suppl S1 reference (.pdf → .png, proper format)
# ─────────────────────────────────────────────────────────────────────────────
p14 = doc.paragraphs[14]
assert "Figure S1" in p14.text or "FigS1" in p14.text, f"Unexpected P14: {p14.text}"

s1_text = (
    "Supplementary Figure S1.  Weighted Kaplan\u2013Meier survival curves for 28-day "
    "all-cause mortality, full analytic cohort (n\u2009=\u20093,112). "
    "See: analysis/figures/FigS1_KM_FullCohort.png."
)
set_para_text(p14, s1_text)
print("Fix 2 done: P14 Suppl S1 reference updated.")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 3 — P16: update S1b reference .pdf → .png
# ─────────────────────────────────────────────────────────────────────────────
p16 = doc.paragraphs[16]
assert "S1b" in p16.text, f"Unexpected P16: {p16.text}"

s1b_ref_text = (
    "Supplementary Figure S1b.  Cumulative incidence of 28-day all-cause death, "
    "full analytic cohort (n\u2009=\u20093,112). "
    "See: analysis/figures/FigS1b_CumInc_Full.png."
)
set_para_text(p16, s1b_ref_text)
print("Fix 3 done: P16 Suppl S1b reference updated.")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 4 — P18: this is a duplicate S1b description sitting in main results;
#          replace with blank (it will be covered by updated P16 above)
# ─────────────────────────────────────────────────────────────────────────────
p18 = doc.paragraphs[18]
if "Supplementary Figure S1b" in p18.text:
    set_para_text(p18, "")
    print("Fix 4 done: P18 duplicate S1b description cleared.")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 5 — Move P26 (Figure 2 caption) to BEFORE P25 (3.3.2 heading)
#          i.e., it goes after Table 3 footnote / blank → before 3.3.2
# ─────────────────────────────────────────────────────────────────────────────
# P25 = "3.3.2  Competing Risks Analysis..." (Heading 2)
# P26 = "Figure 2.  Distribution of RRT-Free Days..." — should be BEFORE P25
p25_elem = doc.paragraphs[25]._p
p26_elem = doc.paragraphs[26]._p

assert "3.3.2" in para_text(p25_elem), f"P25 mismatch: {para_text(p25_elem)[:60]}"
assert "Figure 2" in para_text(p26_elem), f"P26 mismatch: {para_text(p26_elem)[:60]}"

# Move P26 to just before P25 in the body XML
fig2_xml = copy.deepcopy(p26_elem)    # deep copy of Figure 2 paragraph
p26_elem.getparent().remove(p26_elem) # remove from current position
p25_elem.addprevious(fig2_xml)        # insert before 3.3.2 heading
print("Fix 5 done: Figure 2 caption moved before 3.3.2 heading.")

# ─────────────────────────────────────────────────────────────────────────────
# FIX 6 — Move Supplementary Materials block (currently between 3.4 and 3.5)
#          to AFTER the last paragraph of the document (after 3.5.4 + SA-5)
# ─────────────────────────────────────────────────────────────────────────────
# Reload paragraph references after DOM changes
def get_paras():
    return doc.paragraphs

paras = get_paras()

# Find Supplementary Materials heading
suppl_start_idx = None
for i, p in enumerate(paras):
    if p.text.strip() == "Supplementary Materials" and p.style.name.startswith("Heading"):
        suppl_start_idx = i
        break

assert suppl_start_idx is not None, "Could not find Supplementary Materials heading"

# Find 3.5 heading (it comes after Supplementary Materials in current layout)
sa_idx = None
for i, p in enumerate(paras):
    if "3.5" in p.text and "Sensitivity" in p.text:
        sa_idx = i
        break

assert sa_idx is not None, "Could not find 3.5 Sensitivity Analyses heading"

print(f"Supplementary Materials block starts at P{suppl_start_idx}")
print(f"3.5 Sensitivity Analyses at P{sa_idx}")

# Collect the block: from suppl_start_idx up to (but not including) sa_idx
suppl_elems = [copy.deepcopy(p._p) for p in paras[suppl_start_idx:sa_idx]]
print(f"Collecting {len(suppl_elems)} paragraphs for Supplementary Materials block")

# Remove original paragraphs from body
for p in paras[suppl_start_idx:sa_idx]:
    p._p.getparent().remove(p._p)

# Find the new last paragraph of the document (after removal)
paras = get_paras()
last_para = paras[-1]._p

# Also add a blank separator before Supplementary Materials
blank = OxmlElement("w:p")
last_para.addnext(blank)

# Append the supplementary block after the last paragraph
insert_after = blank
for elem in suppl_elems:
    insert_after.addnext(elem)
    insert_after = elem

# ── Also update references inside the moved Supplementary block to .png ──────
paras = get_paras()
for p in paras:
    t = p.text
    # Update any remaining .pdf references to .png
    if ".pdf" in t and ("Figure" in t or "figure" in t):
        new_t = t.replace("FigS1_KM_FullCohort.pdf", "FigS1_KM_FullCohort.png")
        new_t = new_t.replace("Figure2_KM_Equipoise.pdf", "Figure2_KM_Equipoise.png")
        new_t = new_t.replace("FigS1b_CumInc_Full.pdf", "FigS1b_CumInc_Full.png")
        new_t = new_t.replace("Figure2B_CumInc_Equipoise.pdf", "Figure2B_CumInc_Equipoise.png")
        new_t = new_t.replace("FigS2_LovePlot_IPCW_Balance.pdf", "FigS2_LovePlot_IPCW_Balance.png")
        new_t = new_t.replace("Figure4_CIF_RenalRecovery.pdf", "Figure4_CIF_RenalRecovery.pdf")  # keep pdf for Figure 4
        if new_t != t:
            set_para_text(p, new_t)

print("Fix 6 done: Supplementary Materials block moved to end of document.")

# ─────────────────────────────────────────────────────────────────────────────
# SAVE
# ─────────────────────────────────────────────────────────────────────────────
doc.save(DOC_PATH)
print("\nAll fixes applied. result.docx saved.")

# ── Quick verification ────────────────────────────────────────────────────────
doc2 = Document(DOC_PATH)
print("\n=== FINAL PARAGRAPH STRUCTURE ===")
for i, p in enumerate(doc2.paragraphs):
    if p.text.strip():
        print(f"  P{i}: [{p.style.name[:12]}] {repr(p.text[:90])}")
