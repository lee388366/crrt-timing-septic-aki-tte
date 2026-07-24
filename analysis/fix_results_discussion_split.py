#!/usr/bin/env python3
"""
Strip Discussion content from Results section of result.docx.
Save extracted Discussion material to discussion_notes.md.
Principle: Results = WHAT was found. Discussion = WHY.
"""
import copy
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_LINE_SPACING

RESULT_PATH = "analysis/result.docx"
THIN = '\u2009'

doc = Document(RESULT_PATH)

discussion_notes = []
changes_made = []

def get_text(para):
    return para.text.strip()

def replace_para_text(para, new_text):
    """Replace paragraph text while preserving first run's formatting."""
    if not para.runs:
        para.text = new_text
        return
    first_run = para.runs[0]
    font_name = first_run.font.name
    font_size = first_run.font.size
    font_bold = first_run.font.bold
    font_italic = first_run.font.italic
    font_color = first_run.font.color.rgb if first_run.font.color and first_run.font.color.rgb else None

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
        if font_color:
            run.font.color.rgb = font_color

def safe_replace(text, old, new):
    """Replace handling thin spaces."""
    result = text.replace(old, new)
    if result == text:
        old_thin = old.replace(' ', THIN)
        result = text.replace(old_thin, new)
    return result

for i, para in enumerate(doc.paragraphs):
    t = get_text(para)
    if not t:
        continue
    original = t

    # === P003: Remove CCW methodological justification ===
    if "reflecting substantial baseline imbalance that underscores the need for Clone-Censor-Weighting" in t or \
       "reflecting substantial baseline imbalance that underscores the need for Clone" in t.replace(THIN, ' '):
        old_clause = "reflecting substantial baseline imbalance that underscores the need for Clone-Censor-Weighting to achieve causal comparability"
        new_clause = "reflecting substantial baseline imbalance between treatment groups"

        disc_text = "The large baseline imbalance underscores the need for Clone-Censor-Weighting (CCW) to achieve causal comparability between treatment groups."
        discussion_notes.append(("P003 — Baseline imbalance → CCW justification", disc_text))

        new_t = safe_replace(para.text, old_clause, new_clause)
        if new_t != para.text:
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed CCW methodological justification")

    # === P012: Major rewrite — the most problematic paragraph ===
    if "3.2.1" in t and "Cox proportional hazards" in t:
        full_t = para.text

        # Extract all Discussion content before rewriting
        disc_pieces = []
        if "This pattern reflects the expected confounding by indication" in full_t.replace(THIN, ' '):
            disc_pieces.append("The nominally elevated crude HR reflects the expected confounding by indication: patients selected for early CRRT were substantially sicker at t₀ (Table 2), generating a spurious apparent harm signal in the unadjusted analysis.")
        if "The magnitude of attenuation from the crude to the weighted estimate" in full_t.replace(THIN, ' '):
            disc_pieces.append("The magnitude of attenuation from the crude to the weighted estimate (HR 1.098 → 1.038) is consistent with partial correction for the strong severity-driven confounding documented in Table 2.")
        if "residual structural confounding at the severity extremes" in full_t.replace(THIN, ' '):
            disc_pieces.append("The residual statistical significance in the full-cohort analysis reflects contributions from patients at SOFA extremes (SOFA ≤7 and ≥15), where CRRT decisions are less clinically ambiguous and residual structural confounding after IPCW may be more difficult to eliminate — suggesting that residual confounding, rather than a genuine treatment effect, may drive the full-cohort estimate.")
        if "This underscores the importance of the equipoise-restricted" in full_t.replace(THIN, ' '):
            disc_pieces.append("This underscores the importance of the equipoise-restricted primary analysis.")

        if disc_pieces:
            discussion_notes.append(("P012 — Confounding interpretation & equipoise justification", "\n".join(disc_pieces)))

        # Rewrite: data only
        new_p012 = (
            "3.2.1  Cox proportional hazards analysis. "
            "In the unadjusted (crude) analysis of the primary analytic population\u2014"
            "patients in the clinical equipoise zone (SOFA 8\u201314 at t\u2080; n = 1,887)\u2014"
            "early CRRT was associated with a nominally higher hazard of 28-day all-cause "
            "mortality (crude HR = 1.098, 95% CI 0.977 to 1.233; p = 0.117; Table 3). "
            "After application of IPCW weights, the hazard ratio attenuated toward the null "
            "(IPCW-weighted HR = 1.038, 95% CI 0.915 to 1.178; p = 0.565; Table 3).\n"
            "In the full analytic cohort (n = 3,112), the unadjusted Cox model yielded "
            "HR = 1.228 (95% CI 1.113 to 1.355; p < 0.001); after IPCW adjustment, "
            "HR = 1.164 (95% CI 1.047 to 1.294; p = 0.005). "
            "Results for all analytic strata are presented in Table 3."
        )
        replace_para_text(para, new_p012)
        changes_made.append(f"P{i:03d}: Major rewrite — removed ~6 sentences of Discussion content from 3.2.1")

    # === P014: Remove RMST methodological commentary ===
    if "the RMST difference therefore serves as a corroborating absolute effect measure" in t.replace(THIN, ' '):
        old_rmst = "; the RMST difference therefore serves as a corroborating absolute effect measure rather than a replacement for the HR."
        disc_text = "Because the proportional hazards assumption held, the RMST difference serves as a corroborating absolute effect measure rather than a replacement for the HR."
        discussion_notes.append(("P014 — RMST as corroborating measure", disc_text))

        new_t = safe_replace(para.text, old_rmst, ".")
        if new_t != para.text:
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed RMST methodological commentary")

    # === P024-area: Renal recovery interpretation ===
    if "warrants careful interpretation" in t.replace(THIN, ' ') and "renal recovery" in t.lower():
        old_interp = "The lower renal recovery rate in the early CRRT group warrants careful interpretation: patients selected for early CRRT had substantially higher baseline severity (Table 2), and the composite nature of this outcome (incorporating 28-day mortality) means that group differences in recovery largely reflect the higher mortality burden rather than a differential effect on renal repair per se."

        disc_text = "The lower renal recovery rate in the early CRRT group warrants careful interpretation: patients selected for early CRRT had substantially higher baseline severity (Table 2), and the composite nature of this outcome (incorporating 28-day mortality) means that group differences largely reflect the higher mortality burden rather than a differential effect on renal repair per se."
        discussion_notes.append(("Secondary — Renal recovery interpretation", disc_text))

        # Try with and without thin spaces
        new_t = para.text
        for variant in [old_interp, old_interp.replace(' ', THIN)]:
            if variant in new_t:
                new_t = new_t.replace(variant, "")
                break
        # Also try partial match
        if new_t == para.text:
            # Fallback: search for the key phrase
            idx = new_t.find("The lower renal recovery rate")
            if idx == -1:
                idx = new_t.replace(THIN, ' ').find("The lower renal recovery rate")
            if idx >= 0:
                end_idx = new_t.find("per se.", idx)
                if end_idx == -1:
                    end_idx_nts = new_t.replace(THIN, ' ').find("per se.", idx)
                    if end_idx_nts >= 0:
                        end_idx = end_idx_nts
                if end_idx >= 0:
                    removal = new_t[idx:end_idx + len("per se.")]
                    new_t = new_t[:idx].rstrip() + new_t[end_idx + len("per se."):].lstrip()

        if new_t != para.text:
            new_t = new_t.strip()
            if new_t:
                replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed renal recovery interpretation (→ Discussion)")

    # === P025-area: RFD causal interpretation ===
    if "primarily reflects the longer duration of CRRT exposure" in t.replace(THIN, ' '):
        old_rfd = "The difference in this outcome primarily reflects the longer duration of CRRT exposure in the early group and the higher mortality burden documented in the unadjusted analysis, rather than a differential effect on CRRT weaning."

        disc_text = old_rfd
        discussion_notes.append(("Secondary — RFD causal interpretation", disc_text))

        new_t = para.text
        for variant in [old_rfd, old_rfd.replace(' ', THIN)]:
            if variant in new_t:
                new_t = new_t.replace(variant, "").strip()
                break
        if new_t == para.text.strip():
            idx = new_t.replace(THIN, ' ').find("The difference in this outcome primarily reflects")
            if idx >= 0:
                end_marker = "CRRT weaning."
                end_idx = new_t.replace(THIN, ' ').find(end_marker, idx)
                if end_idx >= 0:
                    actual_end = end_idx + len(end_marker)
                    search_text = new_t.replace(THIN, ' ')
                    actual_segment = new_t[idx:actual_end]
                    new_t = (new_t[:idx].rstrip() + " " + new_t[actual_end:].lstrip()).strip()

        if new_t != para.text.strip():
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed RFD causal interpretation (→ Discussion)")

    # === P026-area: ICU LOS clinical significance judgment ===
    if "unlikely to represent a clinically meaningful treatment effect" in t.replace(THIN, ' '):
        old_los = "The modest difference in ICU LOS between groups, while statistically significant in the full cohort, is unlikely to represent a clinically meaningful treatment effect given the large baseline severity imbalance."

        disc_text = old_los
        discussion_notes.append(("Secondary — ICU LOS clinical significance", disc_text))

        new_t = para.text
        for variant in [old_los, old_los.replace(' ', THIN)]:
            if variant in new_t:
                new_t = new_t.replace(variant, "").strip()
                break
        if new_t == para.text.strip():
            idx = new_t.replace(THIN, ' ').find("The modest difference in ICU LOS")
            if idx >= 0:
                end_marker = "severity imbalance."
                end_idx = new_t.replace(THIN, ' ').find(end_marker, idx)
                if end_idx >= 0:
                    actual_end = end_idx + len(end_marker)
                    new_t = (new_t[:idx].rstrip() + " " + new_t[actual_end:].lstrip()).strip()

        if new_t != para.text.strip():
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed ICU LOS clinical significance judgment (→ Discussion)")

    # === P033-area: Fine-Gray methodological justification ===
    if "Because death is a competing event for renal recovery" in t.replace(THIN, ' ') and "standard logistic regression may conflate" in t.replace(THIN, ' '):
        disc_text = "Standard logistic regression may conflate distinct failure modes — patients who die and those whose kidneys fail to recover — because it codes all deaths as 'no recovery'. The Fine-Gray subdistribution hazard model separates these processes, providing a more nuanced estimate of the treatment effect on renal recovery in the presence of competing mortality."
        discussion_notes.append(("Secondary — Fine-Gray methodological justification", disc_text))

        # Rewrite to be brief
        new_p033 = para.text
        # Find the full sentence and replace with concise version
        old_long = "Because death is a competing event for renal recovery"
        idx = new_p033.replace(THIN, ' ').find(old_long)
        if idx >= 0:
            # Find end of the justification paragraph (up to "We therefore applied" or the post-hoc note)
            we_idx = new_p033.replace(THIN, ' ').find("We therefore applied")
            if we_idx < 0:
                we_idx = new_p033.replace(THIN, ' ').find("we therefore applied")

            if we_idx >= 0:
                # Replace the methodological justification with a brief version
                brief = "To account for the competing risk of death, we applied "
                # Find what comes after "We therefore applied"
                after_we = new_p033[we_idx + len("We therefore applied"):]
                after_we_clean = after_we.replace(THIN, ' ')
                new_p033 = new_p033[:idx] + brief + after_we.lstrip()
                replace_para_text(para, new_p033)
                changes_made.append(f"P{i:03d}: Condensed Fine-Gray intro to one sentence (methodology → Discussion)")

    # === P034-area: OR vs sHR discordance interpretation ===
    if "discordance between the adjusted OR" in t.replace(THIN, ' ') or \
       "Standard logistic regression codes all deaths" in t.replace(THIN, ' '):

        disc_pieces_34 = []
        plain = t.replace(THIN, ' ')
        if "The discordance between the adjusted OR" in plain:
            disc_pieces_34.append("The discordance between the adjusted OR (0.514 in logistic regression, indicating lower odds of recovery with early CRRT) and the sHR (1.096, indicating a comparable subdistribution hazard of recovery) underscores the importance of the competing risks framework.")
        if "Standard logistic regression codes all deaths" in plain:
            disc_pieces_34.append("Standard logistic regression codes all deaths as 'no recovery', conflating true impairment of renal repair with the higher competing mortality in the early CRRT group. The Fine-Gray model properly separates these distinct processes; after accounting for the competing risk of death, there was no significant impairment of renal recovery associated with early CRRT.")

        if disc_pieces_34:
            discussion_notes.append(("Secondary — OR vs sHR discordance interpretation", "\n".join(disc_pieces_34)))

        # Remove the interpretation sentences, keep only data
        new_t = para.text
        # Remove "The discordance between..." sentence
        markers = [
            ("The discordance between the adjusted OR", "competing risks framework."),
            ("Standard logistic regression codes all deaths", "early CRRT."),
        ]
        for start_phrase, end_phrase in markers:
            plain_t = new_t.replace(THIN, ' ')
            s_idx = plain_t.find(start_phrase)
            if s_idx >= 0:
                e_idx = plain_t.find(end_phrase, s_idx)
                if e_idx >= 0:
                    actual_end = e_idx + len(end_phrase)
                    new_t = (new_t[:s_idx].rstrip() + " " + new_t[actual_end:].lstrip()).strip()

        if new_t != para.text.strip():
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed OR vs sHR interpretation (→ Discussion)")

    # === P045: E-value interpretation ===
    if "E-value" in t and "consistent with a true null" in t.replace(THIN, ' '):
        disc_text = "The small E-value (close to 1.0) is consistent with a true null or near-null treatment effect: even weak unmeasured confounders would suffice to shift the estimate to null. This supports the interpretation that, after rigorous IPCW-based causal adjustment, early CRRT does not confer a detectable mortality benefit in this high-risk subgroup, rather than a protective effect masked by residual confounding."
        discussion_notes.append(("Sensitivity — E-value interpretation", disc_text))

        new_t = para.text
        # Remove from "The small E-value" to end of that interpretation
        plain_t = new_t.replace(THIN, ' ')
        s_idx = plain_t.find("The small E-value")
        if s_idx >= 0:
            # Find end: "residual confounding."
            e_idx = plain_t.find("residual confounding.", s_idx)
            if e_idx >= 0:
                actual_end = e_idx + len("residual confounding.")
                new_t = (new_t[:s_idx].rstrip() + new_t[actual_end:].lstrip()).strip()

        if new_t != para.text.strip():
            replace_para_text(para, new_t)
            changes_made.append(f"P{i:03d}: Removed E-value interpretation (→ Discussion)")

doc.save(RESULT_PATH)

# Save Discussion notes
disc_path = "analysis/discussion_notes_extracted.md"
with open(disc_path, "w", encoding="utf-8") as f:
    f.write("# Discussion Content Extracted from Results Section\n\n")
    f.write("These sentences/paragraphs were removed from `result.docx` because they contain\n")
    f.write("interpretive, causal, or methodological content that belongs in the Discussion section,\n")
    f.write("not in Results.\n\n")
    f.write("**Principle**: Results = WHAT was found (data, numbers, facts).\n")
    f.write("Discussion = WHY it was found (causal interpretation, clinical significance, methodological justification).\n\n")
    f.write("---\n\n")
    for title, content in discussion_notes:
        f.write(f"## {title}\n\n")
        f.write(f"{content}\n\n")
    f.write("---\n\n")
    f.write("*Extracted on 2026-03-28 during Results section audit.*\n")

print("=" * 60)
print("CHANGES MADE:")
for c in changes_made:
    print(f"  ✓ {c}")
print(f"\nTotal changes: {len(changes_made)}")
print(f"\nDiscussion notes saved to: {disc_path}")
print(f"  ({len(discussion_notes)} items extracted)")
