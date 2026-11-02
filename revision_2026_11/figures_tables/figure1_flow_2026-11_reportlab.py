"""Figure 1 cohort flow (2026-11 revision). Reproduces the layout of the submitted
02_Figure_1_Cohort_Flow.pdf (ReportLab; 479.2 pt wide; main column left, exclusion boxes right).
Counts updated; wording shortened."""
import os, sys
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT

FONTS = [("/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
         ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"),
         ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf")]
for reg, bold in FONTS:
    if os.path.exists(reg) and os.path.exists(bold):
        pdfmetrics.registerFont(TTFont("Sans", reg)); pdfmetrics.registerFont(TTFont("Sans-Bold", bold)); break
from reportlab.lib.fonts import addMapping
addMapping("Sans", 0, 0, "Sans"); addMapping("Sans", 1, 0, "Sans-Bold")

FS, LEAD, PAD = 8.6, 10.6, 5
MAIN = ParagraphStyle("m", fontName="Sans", fontSize=FS, leading=LEAD, alignment=TA_CENTER)
SIDE = ParagraphStyle("s", fontName="Sans", fontSize=FS, leading=LEAD, alignment=TA_LEFT, leftIndent=11, firstLineIndent=-11)
SIDEH = ParagraphStyle("h", parent=SIDE, leftIndent=0, firstLineIndent=0)
SIDEC = ParagraphStyle("c", parent=SIDE, leftIndent=11, firstLineIndent=0)
T0 = "t<sub>0</sub>"

main = [
    ["ICU stays with first physiologic", "KDIGO stage 3 AKI", "MIMIC-IV v3.1 (2008–2022)", "<b>n = 14,625</b>"],
    ["Sepsis-3 with vasopressor support", "<b>n = 3,728</b>"],
    ["First eligible ICU stay per patient", "<b>n = 2,839</b>"],
    ["Analytic cohort", "<b>n = 2,461</b>"],
    [f"Primary population: SOFA 8–14 at {T0}", "<b>n = 1,469</b>", "Early CRRT 139; no CRRT within 24 h 1,330"],
]
side = [
    ["<b>Excluded (n = 10,897)</b>", f"(1) No Sepsis-3 at or before {T0} (n = 5,872)", "(2) No vasopressor in the preceding 6 h (n = 5,025)"],
    ["<b>Excluded (n = 889)</b>", f"(1) RRT before {T0} (n = 419)", "(2) End-stage kidney disease or chronic dialysis (n = 350)",
     "(3) ICU stay &lt;24 h (n = 14)", "(4) Missing baseline data (n = 0)", "(5) No follow-up time (n = 8)", "(6) Repeat ICU stays (n = 98)"],
    ["<b>Excluded (n = 378)</b>", f"(1) Urgent dialysis indication at {T0} (n = 376)",
     ">Potassium ≥6.5 mmol/L, pH ≤7.15, bicarbonate ≤12 mmol/L or urea nitrogen ≥100 mg/dL", "(2) Inconsistent timestamps (n = 2)"],
    ["<b>Observed CRRT initiation</b>", f"Within 24 h of {T0}: 192 (7.8%)", "None within 24 h: 2,269 (92.2%)"],
]

W = 479.2; MX, MW = 12, 196; SX, SW = 236, 232; GAP_MIN = 22
def para(lines, style, w, header_style=None):
    ps = [Paragraph(t[1:], SIDEC) if t.startswith(">") else Paragraph(t, header_style if (header_style and i == 0) else style) for i, t in enumerate(lines)]
    hs = [p.wrap(w - 2 * PAD, 1000)[1] for p in ps]
    return ps, hs, sum(hs) + 2 * PAD
M = [para(l, MAIN, MW) for l in main]; S = [para(l, SIDE, SW, SIDEH) for l in side]
gaps = [max(GAP_MIN, s[2] + 14) for s in S]
H = 14 + sum(m[2] for m in M) + sum(gaps) + 14

def box(c, x, ytop, w, h, ps, hs):
    c.setLineWidth(0.8); c.rect(x, ytop - h, w, h)
    y = ytop - PAD
    for p, hh in zip(ps, hs):
        p.drawOn(c, x + PAD, y - hh); y -= hh
def arrow(c, x1, y1, x2, y2):
    c.setLineWidth(0.8); c.line(x1, y1, x2, y2)
    import math
    a = math.atan2(y2 - y1, x2 - x1); L = 5.5
    p = c.beginPath(); p.moveTo(x2, y2)
    p.lineTo(x2 - L * math.cos(a - 0.38), y2 - L * math.sin(a - 0.38)); p.lineTo(x2 - L * math.cos(a + 0.38), y2 - L * math.sin(a + 0.38)); p.close()
    c.drawPath(p, fill=1, stroke=0)

def draw(path):
    c = canvas.Canvas(path, pagesize=(W, H)); c.setTitle("Figure 1. Cohort flow"); c.setAuthor("")
    y = H - 14; cx = MX + MW / 2
    for i, (ps, hs, h) in enumerate(M):
        box(c, MX, y, MW, h, ps, hs); y -= h
        if i < len(S):
            g = gaps[i]; sps, shs, sh = S[i]
            arrow(c, cx, y, cx, y - g)
            mid = y - g / 2
            box(c, SX, mid + sh / 2, SW, sh, sps, shs)
            c.setLineWidth(0.8); arrow(c, cx, mid, SX, mid)
            y -= g
    c.showPage(); c.save()

out = sys.argv[1] if len(sys.argv) > 1 else "figures/02_Figure_1_Cohort_Flow_2026-11.pdf"
draw(out); print("wrote", out, "page", W, "x", round(H, 1))
