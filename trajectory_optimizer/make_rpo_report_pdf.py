"""
Render RPO_report.md content into a polished PDF (RPO_report.pdf).

Self-contained: builds the report with ReportLab Platypus and embeds the four
trajectory figures.  Uses the DejaVuSans font (shipped with matplotlib) so that
maths glyphs (Delta, x, micro, sqrt, super/subscripts) render correctly.

    cd c:\\Projects\\DSE\\REAVER
    python trajectory_optimizer/make_rpo_report_pdf.py
"""
import os
import matplotlib
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, ListFlowable, ListItem)
from reportlab.lib.utils import ImageReader

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "RPO_report.pdf")

# ── Fonts (DejaVu from matplotlib data dir → full glyph coverage) ────────────
_ttf = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
pdfmetrics.registerFont(TTFont("DejaVu",      os.path.join(_ttf, "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-Bold", os.path.join(_ttf, "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-Obl",  os.path.join(_ttf, "DejaVuSans-Oblique.ttf")))
pdfmetrics.registerFontFamily("DejaVu", normal="DejaVu", bold="DejaVu-Bold",
                              italic="DejaVu-Obl", boldItalic="DejaVu-Bold")

ACCENT = colors.HexColor("#1f4e79")
HEADER = colors.HexColor("#dce6f1")
GRID   = colors.HexColor("#9bb8d3")

styles = getSampleStyleSheet()
def _st(name, **kw):
    base = dict(fontName="DejaVu", fontSize=9.5, leading=13, textColor=colors.black)
    base.update(kw)
    return ParagraphStyle(name, **base)

S = {
    "title":   _st("title",   fontName="DejaVu-Bold", fontSize=18, leading=22, textColor=ACCENT, spaceAfter=2),
    "subtitle":_st("subtitle",fontName="DejaVu-Obl",  fontSize=10, leading=13, textColor=colors.HexColor("#555555"), spaceAfter=10),
    "h1":      _st("h1",       fontName="DejaVu-Bold", fontSize=13, leading=16, textColor=ACCENT, spaceBefore=12, spaceAfter=5),
    "h2":      _st("h2",       fontName="DejaVu-Bold", fontSize=11, leading=14, textColor=colors.HexColor("#244b6b"), spaceBefore=9, spaceAfter=3),
    "body":    _st("body",     alignment=TA_LEFT, spaceAfter=6),
    "bullet":  _st("bullet",   leftIndent=10, spaceAfter=2),
    "cap":     _st("cap",      fontName="DejaVu-Obl", fontSize=8.5, leading=11, textColor=colors.HexColor("#444444"), spaceBefore=2, spaceAfter=10),
    "cell":    _st("cell",     fontSize=8.8, leading=11),
    "cellb":   _st("cellb",    fontName="DejaVu-Bold", fontSize=8.8, leading=11),
    "cellr":   _st("cellr",    fontSize=8.8, leading=11, alignment=2),
    "cellbr":  _st("cellbr",   fontName="DejaVu-Bold", fontSize=8.8, leading=11, alignment=2),
    "foot":    _st("foot",     fontName="DejaVu-Obl", fontSize=8, leading=11, textColor=colors.HexColor("#555555"), spaceBefore=6),
}

PAGE_W, PAGE_H = A4
MARGIN = 18 * mm
USABLE = PAGE_W - 2 * MARGIN

story = []

def P(txt, s="body"):     story.append(Paragraph(txt, S[s]))
def H1(txt):              story.append(Paragraph(txt, S["h1"]))
def H2(txt):              story.append(Paragraph(txt, S["h2"]))
def GAP(h=4):             story.append(Spacer(1, h))

def bullets(items):
    story.append(ListFlowable(
        [ListItem(Paragraph(t, S["bullet"]), leftIndent=12, value="•") for t in items],
        bulletType="bullet", start="•", leftIndent=10, spaceAfter=6))

def figure(fname, caption, width=120 * mm):
    path = os.path.join(HERE, fname)
    iw, ih = ImageReader(path).getSize()
    h = width * ih / iw
    img = Image(path, width=width, height=h)
    img.hAlign = "CENTER"
    story.append(img)
    story.append(Paragraph(caption, S["cap"]))

def result_table(rows, widths):
    """rows: list of (cells, is_total). cells already Paragraphs."""
    data = [[Paragraph("Phase", S["cellb"]),
             Paragraph("ΔV [m/s]", S["cellbr"]),
             Paragraph("Prop [kg]", S["cellbr"]),
             Paragraph("Time", S["cellbr"])]]
    style = [
        ("BACKGROUND", (0, 0), (-1, 0), HEADER),
        ("GRID", (0, 0), (-1, -1), 0.5, GRID),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
    ]
    for i, (cells, total) in enumerate(rows, start=1):
        data.append(cells)
        if total:
            style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#eef3f9")))
            style.append(("LINEABOVE", (0, i), (-1, i), 0.8, ACCENT))
    t = Table(data, colWidths=widths)
    t.setStyle(TableStyle(style))
    story.append(t)
    GAP(10)

def rrow(phase, dv, prop, time, total=False):
    ps, ns = ("cellb", "cellbr") if total else ("cell", "cellr")
    return ([Paragraph(phase, S[ps]), Paragraph(dv, S[ns]),
             Paragraph(prop, S[ns]), Paragraph(time, S[ns])], total)

# =============================================================================
# TITLE
# =============================================================================
P("REAVER — Rendezvous &amp; Proximity Operations (RPO) Analysis", "title")
P("Close-range ΔV budget and AOCS sizing — debris capture and Recycling-Hub handover", "subtitle")

# =============================================================================
# 1. INTRODUCTION
# =============================================================================
H1("1.  Introduction")
P("The REAVER mothership (MS) performs close-range RPO twice per debris: once to "
  "<b>capture</b> a tumbling debris object, and once per cycle to <b>hand it over</b> "
  "at the cooperative Recycling Hub (RH). These maneuvers set the close-range ΔV "
  "budget — the long-range orbit transfers are sized separately by the mission "
  "optimizer — and they drive the AOCS sizing.")
P("This report walks through the two campaigns — <b>Part 1: Debris RPO</b> and "
  "<b>Part 2: RH docking cycle</b> — using four trajectory figures. Each maneuver is a "
  "sequence of short translations between hold points, plus the attitude work "
  "(spin-up, momentum dump, re-point). The figures are the primary explanation aid; "
  "the tables give the per-phase numbers and the equations close the loop. All "
  "trajectories use the target <b>Hill (LVLH) frame</b>: R-bar (radial), V-bar "
  "(along-track), H-bar (cross-track).")

# =============================================================================
# 2. ASSUMPTIONS
# =============================================================================
H1("2.  Assumptions")
_assume = [
    "Relative motion is the linearised <b>Clohessy–Wiltshire</b> model about a circular GEO orbit, mean motion n = √(µ/a³) ≈ 7.11×10<super>-5</super> rad/s.",
    "Each translation is a <b>rest-to-rest two-impulse CW transfer</b>; transfer time keeps the mean speed ≈ ½·V<sub>MAX</sub>, with <b>V<sub>MAX</sub> = 0.05 m/s</b> (low-impact docking limit).",
    "Burns are <b>impulsive</b>; propellant follows <b>Tsiolkovsky</b> with MS monopropellant Isp = 253 s → v<sub>e</sub> = 2481 m/s.",
    "Attitude uses <b>24 RCS thrusters</b> (4 N each); the slew/detumble couple is τ = 2·(3·4 N)·1.05 m = 25.2 N·m. Momentum-dump propellant = H/(arm·v<sub>e</sub>) is thruster-count independent.",
    "Reference debris is <b>COMSATBW-1</b> (2440 kg, 17.2 m span) tumbling at <b>6 °/s (1 rpm)</b> — the REQ-MIS-06 ceiling and AOCS-dimensioning case.",
    "Inertia tensors are <b>box models</b>: MS 3.0×3.0×3.5 m, tug 0.75×0.75×1.2 m, debris bus 2.5×2.5×3.0 m + panels.",
    "<b>Keep-out spheres (KOS)</b> are centred on the target: KOS1 = (panel span/2)·1.20; KOS2 from hardware geometry (see §5).",
    "<b>Slew rates</b> are bounded by V<sub>MAX</sub> at the payload tip (ω ≈ V<sub>MAX</sub>/r): 0.5 °/s for the dock re-point, 0.2 °/s for the clamp settle. Costed as a single momentum build.",
    "The RH is <b>cooperative</b> (known state, fixed docking axis): no spin-synchronisation, only a ~30 min go/no-go dwell.",
    "A <b>+50 % abort margin</b> is applied to the nominal debris-RPO ΔV (not to the detumble line).",
    "<b>TBD (flagged):</b> robotic-arm length (5.0 m placeholder), RH structural size, and debris panel orientation during the RH dock.",
]
adata = [[Paragraph("#", S["cellb"]), Paragraph("Assumption", S["cellb"])]]
for i, a in enumerate(_assume, 1):
    adata.append([Paragraph(str(i), S["cell"]), Paragraph(a, S["cell"])])
at = Table(adata, colWidths=[10 * mm, USABLE - 10 * mm])
at.setStyle(TableStyle([
    ("BACKGROUND", (0, 0), (-1, 0), HEADER),
    ("GRID", (0, 0), (-1, -1), 0.5, GRID),
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("TOPPADDING", (0, 0), (-1, -1), 2.5),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
]))
story.append(at)
GAP(6)

# =============================================================================
# 3. RPO PHASES
# =============================================================================
H1("3.  RPO Phases")

H2("Part 1 — Debris RPO  (Fig 1)")
P("The debris is <b>uncooperative</b> — silent, tumbling, no fixed docking axis — so the "
  "MS cannot just pick a corridor. It first <b>observes from a stable hold</b>, then "
  "approaches along the tumble axis it has measured. The hold is on <b>V-bar</b>, the only "
  "naturally fixed CW hold point (a body at rest at an along-track offset stays put; an "
  "H-bar offset would oscillate back through the plane). The capture corridor is then flown "
  "along the <b>identified tumble axis</b> (modelled as H-bar).")
figure("rpo_fig1_debris_traj.png",
       "<b>Fig 1.</b> Two-stage debris RPO: (blue) approach along −V-bar to the V-bar inspection "
       "hold (20 m); (green) transition around the keep-out sphere onto the tumble axis; (orange) "
       "capture corridor down the tumble axis through KOS1 (10.3 m) → KOS2 (5.6 m, MEV lock), inside the 5° cone.")
bullets([
    "<b>P1 Approach:</b> cut-off (150 m) → 20 m <b>V-bar</b> inspection hold.",
    "<b>P2 Inspection:</b> 8 h dwell on the stable V-bar hold — measure spin rate + tumble axis, ground go/no-go (you must know the axis before approaching it).",
    "<b>P3 Transition:</b> reposition around the keep-out sphere onto the identified tumble axis (H-bar), then close to KOS1.",
    "<b>P4 Arm extension:</b> deploy the arm + tug at KOS1 (before spin-up).",
    "<b>P5 Spin-sync:</b> spin the deployed assembly up to match the tumble (spin-up + inertia growth sum to the same total either way).",
    "<b>P6 MEV lock + detumble:</b> final KOS1→KOS2 approach, then the combined-body momentum dump — the largest ΔV line and the AOCS sizing case.",
    "<b>P7 Retreat:</b> MS backs off to a safe standoff after tug release.",
])

H2("Part 2, Phase A — Tug Meet  (Fig 2a)")
P("At the RH the MS first <b>meets the cooperative tug+debris</b> and clamps it with the "
  "bare robotic arm. Being cooperative, the approach is a direct <b>V-bar</b> corridor — "
  "no H-bar detour, no spin-sync.")
figure("rpo_fig2a_rh_meet.png",
       "<b>Fig 2a.</b> Tug-meet close approach: 30 m standoff → KOS1 (10.3 m) → KOS2 "
       "(5.0 m, <b>bare arm length</b> — no tug on the MS arm yet, so tug/2 is not added).")

H2("Part 2, Phase B — RH Dock  (Fig 2c then Fig 2b)")
P("The MS now carries the combined body (MS + <b>extended arm + tug + debris</b>) back to the "
  "RH and hard-docks via a dedicated port on <b>top</b> of the bus — the arm stays extended, "
  "holding the payload to the side. Split into a far-range carry and a close approach.")
figure("rpo_fig2c_rh_carry.png",
       "<b>Fig 2c.</b> Far-range carry (B1): combined body from the 500 m capture hold to the "
       "30 m standoff. The R-bar excursion is the natural CW coasting arc; green = handoff.")
figure("rpo_fig2b_rh_dock.png",
       "<b>Fig 2b.</b> Close approach + hard-dock: 30 m → KOS1_DOCK (11.0 m, clears the extended "
       "payload) → KOS2_DOCK (2.15 m, top-port contact). The off-axis payload forces the B2 re-point.")

# =============================================================================
# 4. RESULTS
# =============================================================================
H1("4.  Results")
P("ΔV is geometry-driven (mass-independent for translations); propellant is at each part's "
  "reference chaser mass. <b>Mission constants (per event):</b> DV_RPO_DEBRIS = 0.65, "
  "DV_RPO_DETUMBLE = 1.33, DV_RPO_TUG_MEET = 0.14, DV_RPO_RH_DOCK = 0.27, DV_RPO_RH = 0.41 m/s.")

W = [USABLE - 3 * 26 * mm, 26 * mm, 26 * mm, 26 * mm]

H2("Table 1 — Part 1: Debris RPO  (chaser 3228 kg, 6 °/s)")
result_table([
    rrow("P1 Approach (cut-off → V-bar hold)", "0.051", "0.066", "86.7 min"),
    rrow("P2 Inspection (V-bar hold: spin + axis)", "0.005", "0.006", "8.0 h"),
    rrow("P3 Transition to tumble axis (→ KOS1)", "0.100", "0.130", "25.3 min"),
    rrow("P4 Arm extension (deploy + lock onto tug)", "0.011", "0.015", "1.5 min"),
    rrow("P5 Spin-sync (match tumble)", "0.168", "0.219", "0.4 min"),
    rrow("P6a Final approach (KOS1 → KOS2)", "0.050", "0.065", "3.2 min"),
    rrow("P6b Momentum dump (detumble)", "1.334", "2.017", "3.5 min"),
    rrow("P7 Retreat (KOS2 → safe)", "0.050", "0.065", "9.6 min"),
    rrow("DV_RPO_DEBRIS (nominal ×1.50 abort)", "0.653", "0.850", "—", total=True),
    rrow("DV_RPO_DETUMBLE (no abort margin)", "1.334", "2.017", "—", total=True),
], W)

H2("Table 2 — Part 2 Phase A: Tug Meet  (MS 1884 kg)")
result_table([
    rrow("A1 Close approach (30 m → KOS1)", "0.050", "0.038", "13.1 min"),
    rrow("A2 KOS1 hold (synch + comm)", "0.000", "0.000", "30.0 min"),
    rrow("A3 Capture run (KOS1 → KOS2 = arm)", "0.050", "0.038", "3.6 min"),
    rrow("A4 Arm extend + dock (CoM shift)", "0.044", "0.067", "2.0 min"),
    rrow("DV_RPO_TUG_MEET", "0.145", "0.143", "—", total=True),
], W)

H2("Table 3 — Part 2 Phase B: RH Dock  (MS 1884 kg + payload)")
result_table([
    rrow("B1 Carry (500 m hold → 30 m) — Fig 2c", "0.054", "0.099", "313.3 min"),
    rrow("B2 Re-point (combined inertia slew)", "0.111", "0.168", "2.0 min"),
    rrow("B3 Close approach (30 m → KOS1) — Fig 2b", "0.050", "0.091", "12.6 min"),
    rrow("B4 KOS1 hold (align + go/no-go)", "0.000", "0.001", "30.0 min"),
    rrow("B5 Hard-dock (KOS1 → KOS2, top port)", "0.050", "0.091", "5.9 min"),
    rrow("DV_RPO_RH_DOCK", "0.266", "0.449", "—", total=True),
], W)

# =============================================================================
# 5. EQUATIONS
# =============================================================================
H1("5.  Relevant Equations")
bullets([
    "<b>Mean motion (CW):</b>  n = √(µ/a³)  — at GEO ≈ 7.11×10<super>-5</super> rad/s.",
    "<b>Translation (2-impulse CW):</b>  solve Φ<sub>rv</sub>·v<sub>0</sub><super>+</super> = r<sub>f</sub> − Φ<sub>rr</sub>·r<sub>0</sub>, then ΔV = |Δv<sub>1</sub>| + |Δv<sub>2</sub>|, with t = chord/(0.5·V<sub>MAX</sub>).",
    "<b>Stationkeeping:</b>  ΔV = a<sub>dist</sub>·t, with a<sub>dist</sub> = C<sub>R</sub>·(A/m)·P<sub>sol</sub> + a<sub>triax</sub>.",
    "<b>Attitude (detumble / spin-sync / re-point):</b>  H = |J·ω|; m<sub>p</sub> = H/(arm·v<sub>e</sub>); ΔV = v<sub>e</sub>·ln(m/(m−m<sub>p</sub>)); slew time = H/(2·F<sub>edge</sub>·arm); ω ≈ V<sub>MAX</sub>/r.",
    "<b>Propellant (Tsiolkovsky):</b>  m<sub>p</sub> = m<sub>0</sub>·(1 − e<super>−ΔV/v_e</super>),  v<sub>e</sub> = Isp·g<sub>0</sub>.",
    "<b>Keep-out spheres:</b>  KOS1 = (span/2)·1.20;  KOS2(debris) = arm + tug/2 = 5.60 m;  KOS2(tug-meet) = arm = 5.00 m;  KOS2(dock) = MS_h/2 + port = 2.15 m;  KOS1(dock) = (arm + tug + bus)·1.20 ≈ 11.0 m.",
])

P("Figures included: <b>Fig 1</b> (§3 Part 1, Table 1), <b>Fig 2a</b> (§3 Phase A, Table 2), "
  "<b>Fig 2c</b> and <b>Fig 2b</b> (§3 Phase B, Table 3). Figures generated by rpo_plots.py; "
  "numbers by rpo_run.py.", "foot")

doc = SimpleDocTemplate(OUT, pagesize=A4,
                        leftMargin=MARGIN, rightMargin=MARGIN,
                        topMargin=16 * mm, bottomMargin=16 * mm,
                        title="REAVER RPO Analysis")
doc.build(story)
print(f"Saved -> {OUT}")
