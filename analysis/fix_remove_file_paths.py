#!/usr/bin/env python3
"""Remove all internal file path references from result.docx figure captions."""
import re
from docx import Document
from docx.shared import Pt

RESULT_PATH = "analysis/result.docx"
THIN = '\u2009'

doc = Document(RESULT_PATH)
changes = []

for i, para in enumerate(doc.paragraphs):
    original = para.text
    if not original.strip():
        continue

    new_text = original

    # Pattern 1: "See: analysis/figures/xxx.png." or "See: analysis/figures/xxx.png"
    new_text = re.sub(r'\s*See:\s*analysis/figures/[^\.\s]+\.(png|pdf)\.?', '', new_text)

    # Pattern 2: "[See: analysis/figures/xxx.png]"
    new_text = re.sub(r'\s*\[See:\s*analysis/figures/[^\]]+\]', '', new_text)

    # Pattern 3: "[Files: analysis/figures/xxx.png | xxx.pdf]" (various formats)
    new_text = re.sub(r'\s*\[Files?:\s*[^\]]*\]', '', new_text)

    # Clean up trailing whitespace and double spaces
    new_text = re.sub(r'\s{2,}', ' ', new_text).strip()

    if new_text != original.strip():
        # Preserve formatting from first run
        if para.runs:
            first_run = para.runs[0]
            font_name = first_run.font.name
            font_size = first_run.font.size
            font_bold = first_run.font.bold
            font_italic = first_run.font.italic

            for run in para.runs:
                run.text = ""
            para.runs[0].text = new_text

            for run in para.runs:
                run.font.name = font_name or "Times New Roman"
                run.font.size = font_size or Pt(11)
                if font_bold is not None:
                    run.font.bold = font_bold
                if font_italic is not None:
                    run.font.italic = font_italic
        else:
            para.text = new_text

        changes.append(f"P{i:03d}: Removed file path reference")

doc.save(RESULT_PATH)

print("=" * 60)
print("CHANGES MADE:")
for c in changes:
    print(f"  ✓ {c}")
print(f"\nTotal: {len(changes)} paragraphs cleaned")
