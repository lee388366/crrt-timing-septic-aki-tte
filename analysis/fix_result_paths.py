"""
fix_result_paths.py
===================
Append file-path notes to main-text / supplementary figure captions (paths match disk names).
"""

from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

DOC_PATH = "analysis/result.docx"

doc = Document(DOC_PATH)


def append_path_note(para, path_note):
    """Append a new run to existing paragraph with the file path note in italic TNR 10pt."""
    run = para.add_run("  " + path_note)
    run.italic = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(10)
    rPr = run._r.get_or_add_rPr()
    rFonts = OxmlElement("w:rFonts")
    rFonts.set(qn("w:ascii"), "Times New Roman")
    rFonts.set(qn("w:hAnsi"), "Times New Roman")
    rFonts.set(qn("w:eastAsia"), "Times New Roman")
    rPr.insert(0, rFonts)

# ── Track what we changed ─────────────────────────────────────────────────────
changes = []

for i, p in enumerate(doc.paragraphs):
    t = p.text

    # Figure 2B — cumulative incidence (equipoise)
    if t.startswith("Figure 2B.") and "Cumulative Incidence" in t:
        if "Figure2B_CumInc_Equipoise.png" not in t:
            append_path_note(p, "[File: analysis/figures/Figure2B_CumInc_Equipoise.png]")
            changes.append(f"P{i}: added path to Figure 2B caption")

    # Figure 3 — RRT-free days stacked bar
    elif t.startswith("Figure 3.") and "RRT-Free Days" in t:
        if "Figure3_RRTfree_days" not in t:
            append_path_note(
                p,
                "[Files: analysis/figures/Figure3_RRTfree_days.png  |  "
                "Figure3_RRTfree_days.pdf]",
            )
            changes.append(f"P{i}: added path to Figure 3 caption")

    # Figure 4 — CIF renal recovery
    elif t.startswith("Figure 4.") and "CIF" in t:
        if "Figure4_CIF_RenalRecovery" not in t and "Figure4A_CIF" not in t:
            append_path_note(
                p,
                "[Files: analysis/figures/Figure4_CIF_RenalRecovery.pdf (combined)  |  "
                "Figure4A_CIF_Full.png (Panel A)  |  "
                "Figure4B_CIF_Equipoise.png (Panel B)]",
            )
            changes.append(f"P{i}: added path to Figure 4 caption")

    # Supplementary Figure S3 — RFD diagnostic histogram
    elif t.startswith("Supplementary Figure S3.") and "histogram" in t.lower():
        if "FigS3_RFD" not in t:
            append_path_note(
                p,
                "[Files: analysis/figures/FigS3_RFD_Distribution_Diagnostic.png  |  "
                "FigS3_RFD_Distribution_Diagnostic.pdf]",
            )
            changes.append(f"P{i}: added paths to Supplementary Figure S3 (histogram)")


doc.save(DOC_PATH)
print("Saved result.docx\n")
print("Changes made:")
for c in changes:
    print(" •", c)

# ── Quick verification ────────────────────────────────────────────────────────
doc2 = Document(DOC_PATH)
print("\n=== FINAL FIGURE/SUPPL PARAGRAPHS (prefix) ===")
for p in doc2.paragraphs:
    if any(
        kw in p.text
        for kw in ["Figure 2B", "Figure 3.", "Figure 4.", "Supplementary Figure S3"]
    ):
        print(f"  {repr(p.text[:120])}")
