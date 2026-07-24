from __future__ import annotations

import csv
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from docx import Document


BASE_DIR = Path("./广西医科大学")
CSV_PATH = BASE_DIR / "analysis/data/weighted_finegray_results.csv"
DOCX_PATH = BASE_DIR / "analysis/result.docx"
BACKUP_PATH = BASE_DIR / "analysis/result.docx.bak2"

RENAL_ENDPOINT = "Renal recovery (death competing)"
MORTALITY_ENDPOINT = "Mortality (renal recovery competing)"
TBD_PATTERN = r"\[TBD[^\]]*\]"


@dataclass(frozen=True)
class ResultRow:
    endpoint: str
    analysis: str
    n_patients: str
    n_event: str
    n_competing: str
    hr: float
    hr_lo: float
    hr_hi: float
    p_value_display: str

    @property
    def shr_full(self) -> str:
        return (
            f"sHR = {self.hr:.3f} "
            f"(95% CI [{self.hr_lo:.3f}, {self.hr_hi:.3f}]; p {self.p_value_display})"
        )

    @property
    def shr_value_only(self) -> str:
        return (
            f"{self.hr:.3f} "
            f"(95% CI [{self.hr_lo:.3f}, {self.hr_hi:.3f}]; p {self.p_value_display})"
        )

    @property
    def ci_text(self) -> str:
        return f"[{self.hr_lo:.3f}, {self.hr_hi:.3f}]"

    @property
    def hr_text(self) -> str:
        return f"{self.hr:.3f}"


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"Input CSV not found: {path}")

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    required = {
        "Endpoint",
        "Analysis",
        "N_patients",
        "N_event",
        "N_competing",
        "HR",
        "HR_lo",
        "HR_hi",
        "p_value_display",
    }
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ValueError(f"CSV missing required columns: {missing}")
    if not rows:
        raise ValueError(f"CSV has no data rows: {path}")
    return rows


def select_row(
    rows: Iterable[dict[str, str]],
    *,
    endpoint: str,
    analysis_contains: str | None = None,
    analysis_exact: str | None = None,
) -> ResultRow:
    for row in rows:
        if row["Endpoint"] != endpoint:
            continue
        analysis = row["Analysis"]
        if analysis_exact is not None and analysis != analysis_exact:
            continue
        if analysis_contains is not None and analysis_contains not in analysis:
            continue
        return ResultRow(
            endpoint=row["Endpoint"],
            analysis=analysis,
            n_patients=row["N_patients"],
            n_event=row["N_event"],
            n_competing=row["N_competing"],
            hr=float(row["HR"]),
            hr_lo=float(row["HR_lo"]),
            hr_hi=float(row["HR_hi"]),
            p_value_display=row["p_value_display"],
        )
    selector = analysis_exact if analysis_exact is not None else f"contains '{analysis_contains}'"
    raise LookupError(f"No CSV row matched endpoint={endpoint!r}, analysis={selector}")


def clear_paragraph(paragraph) -> None:
    element = paragraph._p
    for child in list(element):
        if child.tag.endswith("}pPr"):
            continue
        element.remove(child)


def set_paragraph_text(paragraph, new_text: str) -> None:
    style = paragraph.style
    alignment = paragraph.alignment
    clear_paragraph(paragraph)
    paragraph.style = style
    paragraph.alignment = alignment
    paragraph.add_run(new_text)


def set_cell_text(cell, new_text: str) -> None:
    if not cell.paragraphs:
        cell.text = new_text
        return
    first = cell.paragraphs[0]
    style = first.style
    alignment = first.alignment
    cell.text = ""
    first = cell.paragraphs[0]
    first.style = style
    first.alignment = alignment
    first.add_run(new_text)


def replace_nth_tbd(text: str, replacement: str, occurrence_index: int) -> tuple[str, str]:
    matches = list(re.finditer(TBD_PATTERN, text))
    if occurrence_index >= len(matches):
        raise ValueError(f"Paragraph has only {len(matches)} TBD placeholders; need index {occurrence_index}")
    match = matches[occurrence_index]
    before = match.group(0)
    after = replacement
    new_text = text[: match.start()] + replacement + text[match.end() :]
    return new_text, f"{before} -> {after}"


def backup_docx(src: Path, dest: Path) -> None:
    shutil.copy2(src, dest)


def update_tbd_paragraphs(doc: Document, renal_full: ResultRow, renal_eq: ResultRow, mortality: list[ResultRow]) -> list[str]:
    logs: list[str] = []
    mortality_iter = iter(mortality)

    for idx, paragraph in enumerate(doc.paragraphs, start=1):
        text = paragraph.text
        if "[TBD" not in text:
            continue

        original = text
        new_text = text
        replacements: list[str] = []

        if "renal recovery compared with deferred CRRT" in new_text:
            new_text, log_item = replace_nth_tbd(new_text, renal_full.shr_value_only, 0)
            replacements.append(f"full renal recovery: {log_item}")
        if "the corresponding sHR was" in new_text and "[TBD" in new_text:
            new_text, log_item = replace_nth_tbd(new_text, renal_eq.shr_value_only, 0)
            replacements.append(f"equipoise renal recovery: {log_item}")

        while "[TBD" in new_text:
            try:
                next_result = next(mortality_iter)
            except StopIteration as exc:
                raise RuntimeError(
                    f"Paragraph {idx} still has unresolved TBD placeholder(s), but no mortality results remain to map."
                ) from exc
            new_text, log_item = replace_nth_tbd(new_text, next_result.shr_full, 0)
            replacements.append(f"{next_result.analysis} / {next_result.endpoint}: {log_item}")

        if new_text != original:
            set_paragraph_text(paragraph, new_text)
            logs.append(
                f"Paragraph {idx}\nBEFORE: {original}\nAFTER:  {new_text}\nDETAILS: " + "; ".join(replacements)
            )

    return logs


def update_table5(doc: Document, renal_full: ResultRow, renal_eq: ResultRow) -> list[str]:
    logs: list[str] = []

    for table_idx, table in enumerate(doc.tables, start=1):
        header = [cell.text.strip() for cell in table.rows[0].cells] if table.rows else []
        if "sHR" not in " ".join(header):
            continue

        for row_idx, row in enumerate(table.rows[1:], start=2):
            cells = row.cells
            if not cells:
                continue
            label = cells[0].text.strip()
            if label == "Full cohort":
                result = renal_full
            elif "Equipoise" in label:
                result = renal_eq
            else:
                continue

            old_values = [cell.text for cell in cells[:7]]
            new_values = [
                label,
                format_int_like(result.n_patients),
                format_int_like(result.n_event),
                format_int_like(result.n_competing),
                result.hr_text,
                result.ci_text,
                result.p_value_display,
            ]

            for cell, new_value in zip(cells[:7], new_values):
                set_cell_text(cell, new_value)

            logs.append(
                f"Table {table_idx} row {row_idx}\nBEFORE: {old_values}\nAFTER:  {new_values}"
            )

    return logs


def format_int_like(value: str) -> str:
    stripped = value.replace(",", "").strip()
    if stripped.isdigit():
        return f"{int(stripped):,}"
    return value


def main() -> int:
    rows = load_csv_rows(CSV_PATH)

    renal_full = select_row(rows, endpoint=RENAL_ENDPOINT, analysis_exact="Full cohort")
    renal_eq = select_row(rows, endpoint=RENAL_ENDPOINT, analysis_contains="Equipoise")
    mortality_full = select_row(rows, endpoint=MORTALITY_ENDPOINT, analysis_exact="Full cohort")
    mortality_eq = select_row(rows, endpoint=MORTALITY_ENDPOINT, analysis_contains="Equipoise")

    if not DOCX_PATH.exists():
        raise FileNotFoundError(f"Word document not found: {DOCX_PATH}")

    backup_docx(DOCX_PATH, BACKUP_PATH)
    print(f"Backup created: {BACKUP_PATH}")

    doc = Document(DOCX_PATH)

    paragraph_logs = update_tbd_paragraphs(
        doc,
        renal_full=renal_full,
        renal_eq=renal_eq,
        mortality=[mortality_full, mortality_eq],
    )
    table_logs = update_table5(doc, renal_full=renal_full, renal_eq=renal_eq)

    doc.save(DOCX_PATH)
    print(f"Saved updated document: {DOCX_PATH}")

    if paragraph_logs:
        print("\nParagraph replacements:")
        for log in paragraph_logs:
            print(log)
            print("-" * 80)
    else:
        print("\nParagraph replacements: none")

    if table_logs:
        print("\nTable replacements:")
        for log in table_logs:
            print(log)
            print("-" * 80)
    else:
        print("\nTable replacements: none")

    print("\nExtracted results:")
    for result in (renal_full, renal_eq, mortality_full, mortality_eq):
        print(
            f"{result.analysis} | {result.endpoint} | "
            f"N={result.n_patients}, event={result.n_event}, competing={result.n_competing}, {result.shr_full}"
        )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
