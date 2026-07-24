from pathlib import Path

from docx import Document


DOC_PATH = Path("./广西医科大学/methods_section.docx")
APPEND_TEXT = (
    " Results are visualized as a forest plot (Figure 5) displaying the sHR "
    "with 95% CI for each analytic population and endpoint."
)


def normalize(text: str) -> str:
    return " ".join(text.split())


def main() -> None:
    doc = Document(DOC_PATH)
    target_idx = None

    for idx, para in enumerate(doc.paragraphs):
        text = normalize(para.text)
        if (
            "secondary outcome analyses" in text.lower()
            and "fine-gray" in text.lower()
            and "subdistribution hazard ratio" in text.lower()
        ):
            target_idx = idx
            break

    if target_idx is None:
        for idx, para in enumerate(doc.paragraphs):
            text = normalize(para.text)
            if "fine-gray" in text.lower() and "secondary outcome analyses" in text.lower():
                target_idx = idx
                break

    if target_idx is None:
        raise RuntimeError("Could not locate the Fine-Gray methods paragraph.")

    para = doc.paragraphs[target_idx]
    current_text = para.text
    if APPEND_TEXT.strip() not in normalize(current_text):
        para.add_run(APPEND_TEXT)
        doc.save(DOC_PATH)
    else:
        doc.save(DOC_PATH)

    verify_doc = Document(DOC_PATH)
    verify_text = normalize(verify_doc.paragraphs[target_idx].text)
    if APPEND_TEXT.strip() not in verify_text:
        raise RuntimeError("Verification failed: appended sentence not found.")

    print(f"Updated paragraph {target_idx + 1} in {DOC_PATH}")
    print(APPEND_TEXT.strip())


if __name__ == "__main__":
    main()
