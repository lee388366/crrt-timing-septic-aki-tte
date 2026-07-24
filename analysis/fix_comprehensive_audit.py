"""
fix_comprehensive_audit.py — Opus 4.6 逐句审计修复
"""
from docx import Document
from docx.shared import Pt
import re

RESULT = "analysis/result.docx"
doc = Document(RESULT)
paras = doc.paragraphs
changes = []
THIN = '\u2009'  # thin space used throughout doc

def safe_replace(para, old, new, label):
    """Replace text in paragraph preserving first-run formatting."""
    t = para.text
    if old not in t:
        changes.append(f"SKIP {label}: pattern not found")
        return False
    new_t = t.replace(old, new)
    if not para.runs:
        return False
    saved = {k: getattr(para.runs[0].font, k) for k in ['name', 'size', 'bold', 'italic']}
    for r in para.runs:
        r.text = ""
    para.runs[0].text = new_t
    for k, v in saved.items():
        setattr(para.runs[0].font, k, v)
    changes.append(f"OK {label}")
    return True


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 1 (HIGH): P012 — "baseline SOFA 8–14" → "SOFA 8–14 at t₀"
#   And add Table 3 cross-reference for the HR results
# ════════════════════════════════════════════════════════════════════════════
p12 = paras[12]
t12 = p12.text

# 1a: Fix "baseline SOFA" → "SOFA at t₀"
safe_replace(p12,
    "clinical equipoise zone (baseline SOFA 8\u201314",
    "clinical equipoise zone (SOFA 8\u201314 at t\u2080",
    "I1a [12]: baseline SOFA → SOFA at t₀")

# 1b: Add Table 3 reference after equipoise IPCW result
p12 = paras[12]  # re-read after modification
safe_replace(p12,
    "there was no statistically detectable difference in 28-day mortality between early and deferred CRRT strategies in this equipoise population.",
    "there was no statistically detectable difference in 28-day mortality between early and deferred CRRT strategies in this equipoise population (Table 3).",
    "I1b [12]: Add Table 3 reference")

# 1c: Add Table 3 reference for full-cohort result
p12 = paras[12]
safe_replace(p12,
    "This underscores the importance of the equipoise-restricted primary analysis.",
    "These results are summarised in Table 3. This underscores the importance of the equipoise-restricted primary analysis.",
    "I1c [12]: Add Table 3 ref for full cohort")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 2 (HIGH): P025 — Remove Wilcoxon p-values from RFD text
#   Methodology: RFD is "descriptive only, no formal inference"
#   Table 4 footnote already labels Wilcoxon as "exploratory"
# ════════════════════════════════════════════════════════════════════════════
p25 = paras[25]
t25 = p25.text

# Remove "Wilcoxon p <0.001" (with thin space) from full cohort
# Pattern: "; Wilcoxon p <0.001" or "; Wilcoxon p<0.001"
new25 = t25
new25 = re.sub(r';\s*Wilcoxon\s+p\s*[<>≤≥]\s*0\.001', '', new25)
# Remove "; p <0.001" or "; p<0.001" from equipoise (without "Wilcoxon" prefix)
# But be careful not to remove legitimate p-values — only after RFD medians
# The equipoise sentence: "...vs 1.5 [0.0, 28.0] days; p <0.001)"
new25 = re.sub(r'days;\s*p\s*[<>]\s*0\.001\)', 'days).', new25)
# Fix double period if any
new25 = new25.replace(')..', ').').replace(').).', ').')
if new25 != t25:
    if not p25.runs:
        changes.append("SKIP I2: no runs")
    else:
        saved = {k: getattr(p25.runs[0].font, k) for k in ['name','size','bold','italic']}
        for r in p25.runs:
            r.text = ""
        p25.runs[0].text = new25
        for k, v in saved.items():
            setattr(p25.runs[0].font, k, v)
        changes.append("OK I2 [25]: Removed Wilcoxon p-values from RFD text (descriptive only)")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 3 (HIGH): P050 — "p_LRT" → "p" (Wald test, not LRT)
#   Can't compute LRT with robust sandwich SEs
# ════════════════════════════════════════════════════════════════════════════
p50 = paras[50]
safe_replace(p50,
    "p_LRT = 0.215",
    "p = 0.215",
    "I3a [50]: p_LRT → p (full cohort interaction)")
p50 = paras[50]
safe_replace(p50,
    "p_LRT = 0.992",
    "p = 0.992",
    "I3b [50]: p_LRT → p (equipoise interaction)")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 4 (HIGH): P029 — Add equipoise definition note for secondary outcomes
#   Renal recovery uses sofa_baseline (n≈1,952), not sofa_t0 (n=1,887)
# ════════════════════════════════════════════════════════════════════════════
p29 = paras[29]
t29 = p29.text
note_equip = (
    " Note: the equipoise subgroup for secondary outcomes "
    "(renal recovery, RRT-free days) is defined using sofa_baseline "
    "(SOFA at ICU admission; n \u2248 1,952), which differs slightly from "
    "the primary analysis equipoise defined using sofa_t\u2080 "
    "(SOFA at index time t\u2080; n = 1,887)."
)
if "sofa_baseline" not in t29 and "1,952" not in t29:
    new29 = t29.rstrip() + note_equip
    if p29.runs:
        saved = {k: getattr(p29.runs[0].font, k) for k in ['name','size','bold','italic']}
        for r in p29.runs:
            r.text = ""
        p29.runs[0].text = new29
        for k, v in saved.items():
            setattr(p29.runs[0].font, k, v)
        changes.append("OK I4 [29]: Added equipoise definition note (sofa_baseline vs sofa_t0)")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 5 (MEDIUM): P036 — Table 5 (Fine-Gray) footnote: clarify clone-level events
# ════════════════════════════════════════════════════════════════════════════
p36 = paras[36]
t36 = p36.text
clone_note = (
    " Event counts in this table are reported at the clone level "
    "(each patient is represented by two clones in the CCW framework), "
    "and are therefore approximately twice the patient-level event "
    "counts reported in the text (e.g., 434 clone-level recovery events "
    "\u2248 217 patient-level events \u00d7 2)."
)
if "clone level" not in t36 and "clone-level" not in t36:
    new36 = t36.rstrip() + clone_note
    if p36.runs:
        saved = {k: getattr(p36.runs[0].font, k) for k in ['name','size','bold','italic']}
        for r in p36.runs:
            r.text = ""
        p36.runs[0].text = new36
        for k, v in saved.items():
            setattr(p36.runs[0].font, k, v)
        changes.append("OK I5 [36]: Fine-Gray Table footnote — clone-level events clarified")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 6 (MEDIUM): P052 — Table 6 footnote: clarify Events column
#   Events = A_early (treatment strategy) clone events, verified against data:
#   Full: 502, Equip: 349, SA-1: 265, etc.
# ════════════════════════════════════════════════════════════════════════════
p52 = paras[52]
t52 = p52.text
events_note = (
    " 'Events' refers to the number of outcome events observed in the "
    "early CRRT strategy (A_early) clones within the CCW framework; "
    "total model events (both strategy arms combined) are larger."
)
if "A_early" not in t52:
    new52 = t52.rstrip() + events_note
    if p52.runs:
        saved = {k: getattr(p52.runs[0].font, k) for k in ['name','size','bold','italic']}
        for r in p52.runs:
            r.text = ""
        p52.runs[0].text = new52
        for k, v in saved.items():
            setattr(p52.runs[0].font, k, v)
        changes.append("OK I6 [52]: Table 6 footnote — Events column clarified (A_early clone events)")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 7 (MEDIUM): P044 — "death events = 265" clarify source
# ════════════════════════════════════════════════════════════════════════════
p44 = paras[44]
safe_replace(p44,
    "death events = 265",
    "death events = 265 in the early CRRT strategy arm",
    "I7 [44]: SA-1 events clarified as A_early clone events")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 8 (MEDIUM): Consistent SOFA terminology
#   P014 "baseline SOFA 8–14" already fixed in P012.
#   Check P019 (Figure 1 caption) and P044 (SA-1)
# ════════════════════════════════════════════════════════════════════════════
# P019 already says "SOFA 8-14" without "baseline" → OK
# P044 says "baseline SOFA scores of 9-13" → should clarify
p44 = paras[44]
safe_replace(p44,
    "baseline SOFA scores of 9-13",
    "SOFA scores at t\u2080 of 9\u201313",
    "I8 [44]: baseline SOFA → SOFA at t₀ for SA-1")


# ════════════════════════════════════════════════════════════════════════════
# ISSUE 9 (LOW): P012 — mention full-cohort crude HR for completeness
#   Currently only equipoise crude and full IPCW are in text;
#   full crude (HR=1.228) is only in Table 3 and Figure 4 caption
# ════════════════════════════════════════════════════════════════════════════
p12 = paras[12]
t12_cur = p12.text
target = f"In the full analytic cohort (n{THIN}={THIN}3,112), the IPCW-weighted Cox model yielded HR{THIN}={THIN}1.164"
if target in t12_cur:
    replacement = (
        f"In the full analytic cohort (n{THIN}={THIN}3,112), "
        f"the unadjusted Cox model showed a statistically significant association "
        f"(crude HR{THIN}={THIN}1.228, 95%{THIN}CI 1.113 to 1.355; p{THIN}<{THIN}0.001). "
        f"After IPCW adjustment, the IPCW-weighted Cox model yielded HR{THIN}={THIN}1.164"
    )
    new_t12 = t12_cur.replace(target, replacement)
    if p12.runs:
        saved = {k: getattr(p12.runs[0].font, k) for k in ['name','size','bold','italic']}
        for r in p12.runs:
            r.text = ""
        p12.runs[0].text = new_t12
        for k, v in saved.items():
            setattr(p12.runs[0].font, k, v)
        changes.append("OK I9 [12]: Added full-cohort crude HR for narrative completeness")


# ════════════════════════════════════════════════════════════════════════════
# Save and report
# ════════════════════════════════════════════════════════════════════════════
doc.save(RESULT)
ok_count = sum(1 for c in changes if c.startswith("OK"))
skip_count = sum(1 for c in changes if c.startswith("SKIP"))
print(f"\n{'='*60}")
print(f"AUDIT FIX COMPLETE: {ok_count} applied, {skip_count} skipped")
print(f"{'='*60}\n")
for c in changes:
    print(f"  {'✓' if c.startswith('OK') else '⚠'} {c}")
