from pathlib import Path

from docx import Document


DOC_PATH = Path("analysis/result.docx")
APPEND_TEXT = " These results are summarized visually in Figure 5."
INSERT_TEXT = (
    "Figure 5. Forest plot of IPCW-weighted subdistribution hazard ratios "
    "(sHR) from Fine-Gray competing risks models. Four estimates are shown: "
    "28-day renal recovery (full cohort and equipoise subgroup, SOFA 8–14) "
    "and 28-day mortality (full cohort and equipoise subgroup), each with "
    "the opposing event as the competing risk. Point estimates and 95% "
    "confidence intervals are derived from IPCW-weighted coxph on the "
    "finegray-expanded dataset (weights: ipcw_w_trim). The vertical dashed "
    "line marks sHR = 1 (no effect). Blue = renal recovery endpoint; Red = "
    "mortality endpoint."
)


def normalize(text: str) -> str:
    return " ".join(text.split())


def insert_paragraph_after(paragraph, text: str):
    new_p = paragraph._element.addnext(paragraph._element.__class__("w:p"))
    inserted = paragraph._parent.paragraphs
    # Refresh document paragraph list by reloading after save.
    return new_p


def main() -> None:
    doc = Document(DOC_PATH)
    numeric_idx = None
    figure3_idx = None

    for idx, para in enumerate(doc.paragraphs):
        text = normalize(para.text)
        text_l = text.lower()
        if numeric_idx is None and all(
            key in text_l for key in ["subdistribution hazard ratio", "full cohort", "equipoise subgroup"]
        ):
            numeric_idx = idx
        if figure3_idx is None and text.startswith("Figure 3."):
            figure3_idx = idx

    if numeric_idx is None:
        for idx, para in enumerate(doc.paragraphs):
            text_l = normalize(para.text).lower()
            if "shr" in text_l and "equipoise subgroup" in text_l:
                numeric_idx = idx
                break

    if figure3_idx is None:
        for idx, para in enumerate(doc.paragraphs):
            text_l = normalize(para.text).lower()
            if "figure 3" in text_l and "cumulative incidence" in text_l:
                figure3_idx = idx
                break

    if numeric_idx is None:
        raise RuntimeError("Could not locate the sHR results paragraph.")
    if figure3_idx is None:
        raise RuntimeError("Could not locate the Figure 3 paragraph.")

    numeric_para = doc.paragraphs[numeric_idx]
    if APPEND_TEXT.strip() not in normalize(numeric_para.text):
        numeric_para.add_run(APPEND_TEXT)

    next_text = ""
    if figure3_idx + 1 < len(doc.paragraphs):
        next_text = normalize(doc.paragraphs[figure3_idx + 1].text)

    if normalize(INSERT_TEXT) != next_text:
        from docx.oxml import OxmlElement
        from docx.text.paragraph import Paragraph

        new_p = OxmlElement("w:p")
        doc.paragraphs[figure3_idx]._element.addnext(new_p)
        new_para = Paragraph(new_p, doc.paragraphs[figure3_idx]._parent)
        new_para.add_run(INSERT_TEXT)

    doc.save(DOC_PATH)

    verify_doc = Document(DOC_PATH)
    verify_numeric = any(APPEND_TEXT.strip() in normalize(p.text) for p in verify_doc.paragraphs)
    verify_insert = any(normalize(p.text) == normalize(INSERT_TEXT) for p in verify_doc.paragraphs)
    if not verify_numeric or not verify_insert:
      raise RuntimeError("Verification failed for result.docx edits.")

    print(f"Updated {DOC_PATH}")
    print(APPEND_TEXT.strip())
    print(INSERT_TEXT)


if __name__ == "__main__":
    main()
