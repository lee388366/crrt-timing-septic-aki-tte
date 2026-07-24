"""
fix_f1_to_f11.py  —  Apply F1–F11 fixes to result.docx
"""
from docx import Document
from docx.shared import Pt
import re, copy

RESULT = "analysis/result.docx"

doc = Document(RESULT)
paras = doc.paragraphs
changes = []

def rewrite_para(para, new_text, pt=11, font="Times New Roman"):
    """Replace all run text with new_text, preserving first-run formatting."""
    if not para.runs:
        return
    saved = {
        "bold": para.runs[0].bold,
        "italic": para.runs[0].italic,
        "size": para.runs[0].font.size,
        "name": para.runs[0].font.name,
    }
    for r in para.runs:
        r.text = ""
    para.runs[0].text = new_text
    para.runs[0].bold = saved["bold"]
    para.runs[0].italic = saved["italic"]
    para.runs[0].font.name = saved["name"] or font
    para.runs[0].font.size = saved["size"] or Pt(pt)

# ─── F1: [17] Table 3 footnote ───────────────────────────────────────────────
# "confirming the validity" → "indicating no statistically significant violation"
p17 = paras[17]
old = p17.text
if "confirming the validity of the Cox HR as the primary inferential measure" in old:
    new = old.replace(
        "confirming the validity of the Cox HR as the primary inferential measure",
        "indicating no statistically significant violation of the proportional hazards assumption, "
        "and supporting the use of the Cox HR as the primary inferential measure"
    )
    rewrite_para(p17, new)
    changes.append("F1 [17]: Schoenfeld footnote wording updated")

# ─── F2: [31] Figure 2 caption ───────────────────────────────────────────────
p31 = paras[31]
old = p31.text
if "as all had initiated CRRT before day 24" in old:
    new = old.replace(
        "as all had initiated CRRT before day 24",
        "as all patients who received early CRRT were on RRT for at least one day, "
        "precluding full 28-day RRT-free survival"
    )
    rewrite_para(p31, new)
    changes.append("F2 [31]: Figure 2 caption — corrected explanation for 0 full-RFD in early group")

# ─── F3: [42] Love plot note ─────────────────────────────────────────────────
p42 = paras[42]
old = p42.text
note = (" Note: SMDs in the Love plot reflect the IPCW-weighted CCW population "
        "(using the clone-level dataset with stabilised weights), whereas SMDs in Table 2 "
        "were computed in the original unweighted cohort; the two sets of SMDs are therefore "
        "not directly comparable.")
if "Love plot" in old and note not in old:
    new = old.rstrip() + note
    rewrite_para(p42, new)
    changes.append("F3 [42]: Love plot dataset-difference note added")

# ─── F4: [46] creatinine variable label ──────────────────────────────────────
p46 = paras[46]
old = p46.text
if "serum creatinine at t0 (median = 1.90 mg/dL)" in old:
    new = old.replace(
        "serum creatinine at t0 (median = 1.90 mg/dL)",
        "serum creatinine at t₀ as derived from the CCW time-varying covariate pipeline "
        "(median = 1.90 mg/dL; <1% missing; this variable differs from the raw "
        "laboratory-recorded creatinine at t₀ in Table 1 [median = 2.30 mg/dL, 36.3% missing], "
        "which reflects point-in-time laboratory values without interpolation)"
    )
    rewrite_para(p46, new)
    changes.append("F4 [46]: Creatinine variable label clarified (CCW pipeline vs raw lab)")

# ─── F5: [46] SOFA median clarification ──────────────────────────────────────
p46 = paras[46]  # re-read after potential F4 modification
old = p46.text
if "SOFA score (median = 9)" in old and "sofa_t0" not in old:
    new = old.replace(
        "SOFA score (median = 9)",
        "SOFA score at t₀ (median = 9; reflecting sofa_t0 in the CCW dataset, "
        "which may differ from the sofa_baseline reported in Table 2 [median = 10], "
        "as the two variables capture SOFA at slightly different time windows)"
    )
    rewrite_para(p46, new)
    changes.append("F5 [46]: SOFA median = 9 vs 10 discrepancy explained (sofa_t0 vs sofa_baseline)")

# Also update Table 6 footnote [52]
p52 = paras[52]
old52 = p52.text
if "SOFA = 9" in old52 and "sofa_t0" not in old52:
    new52 = old52.replace(
        "SOFA = 9",
        "SOFA at t₀ = 9 (sofa_t0 from CCW dataset; sofa_baseline in Table 2 = 10)"
    )
    rewrite_para(p52, new52)
    changes.append("F5 [52]: Table 6 footnote — SOFA median clarified")

# ─── F6: [3] Add age difference, fix SMD list ────────────────────────────────
p3 = paras[3]
old = p3.text
if "Large standardised mean differences (SMD ≥0.10) were observed for SOFA score" in old \
        and "age" not in old.split("Large")[1][:200]:
    new = old.replace(
        "Large standardised mean differences (SMD ≥0.10) were observed for SOFA score, "
        "APACHE III score, lactate, vasopressor dose, pH, and urine output,",
        "Large standardised mean differences (SMD ≥0.10) were observed for multiple "
        "covariates, including age (median 62 [IQR 52–72] years in the early CRRT group "
        "vs 68 [58–79] years in the deferred group; SMD = 0.440), SOFA score, APACHE III "
        "score, lactate, vasopressor dose, pH, and urine output,"
    )
    rewrite_para(p3, new)
    changes.append("F6 [3]: Age difference (62 vs 68 yrs, SMD=0.440) added to SMD list; 'including' qualifier added")

# ─── F7: [3] Remove duplicate Table reference at end ─────────────────────────
p3 = paras[3]
old = p3.text
dup = ("Detailed baseline characteristics of the overall cohort are presented in Table 1; "
       "a comparison of baseline characteristics by treatment group, including standardised mean differences, "
       "is provided in Table 2.")
if dup in old and old.count("Table 1") >= 2:
    new = old.replace(" " + dup, "").replace(dup, "")
    rewrite_para(p3, new)
    changes.append("F7 [3]: Duplicate Table 1/2 reference at end of para removed")

# ─── F8: [13] Add Bootstrap to full-cohort RMST ──────────────────────────────
p13 = paras[13]
old = p13.text
if "full cohort, the RMST difference was -2.13 days (95% CI" in old:
    new = old.replace(
        "full cohort, the RMST difference was -2.13 days (95% CI",
        "full cohort, the RMST difference was -2.13 days (Bootstrap 95% CI"
    )
    rewrite_para(p13, new)
    changes.append("F8 [13]: 'Bootstrap' added before full-cohort RMST CI")

# ─── F9: [25] "pattern was similar" → more precise; fix "p <0.001" spacing ──
p25 = paras[25]
old = p25.text
new = old
if "the pattern was similar" in new:
    new = new.replace("the pattern was similar", "a similar distributional pattern was observed")
# Fix space before 0.001
new = re.sub(r'p <0\.001', 'p < 0.001', new)
if new != old:
    rewrite_para(p25, new)
    changes.append("F9 [25]: 'pattern was similar' → 'similar distributional pattern'; 'p <0.001' → 'p < 0.001'")

# ─── F10: [34] "summarized" → "summarised" ───────────────────────────────────
p34 = paras[34]
old = p34.text
if "summarized visually" in old:
    new = old.replace("summarized visually", "summarised visually")
    rewrite_para(p34, new)
    changes.append("F10 [34]: 'summarized' → 'summarised'")

# ─── F11: [12] "N = 3,112" → "n = 3,112" ────────────────────────────────────
p12 = paras[12]
old = p12.text
if "(N = 3,112)" in old:
    new = old.replace("(N = 3,112)", "(n = 3,112)")
    rewrite_para(p12, new)
    changes.append("F11 [12]: 'N = 3,112' → 'n = 3,112'")

# ─── Save ─────────────────────────────────────────────────────────────────────
doc.save(RESULT)
print(f"[DONE] result.docx saved — {len(changes)} changes applied:\n")
for c in changes:
    print(f"  ✓ {c}")
