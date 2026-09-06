"""Revision figures via frontend construction: self-contained HTML/SVG -> Chrome headless.

Visual language: Nature-style NPG palette on the original submission's clean
grid language (hairline grids, dashed reference lines, right-aligned value
columns, no text-on-marker overlap). Compact panels: Fig1 A/B same width with
the timeline truncated at +5y (post-diagnosis follow-up is not the study
focus); Fig2 and Fig3 panels share identical plot heights and plot widths.
Grey footnote strips removed; only the minimal marker legend is kept in-figure
(full explanations live in the manuscript figure legends).

Content: the domain-correct v4.7.0 figure set (CSV estimates asserted
before drawing).

Outputs: figN_rev1.html (tmp/fig_html), figN_rev1.png (scale 3) and
figN_rev1.pdf (vector) into revision/submission_ready_v4.7.0_20260902/.
"""

from __future__ import annotations

import hashlib
import json
from rev5_domain_results import load_results
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MASTER_PATH = ROOT / "outputs" / "revision_round1_v4_r461_sens2007" / "master_results.csv"
COHORT_PATH = ROOT / "outputs" / "cohort" / "analytic_cohort_v4inc2007.parquet"
OUT_DIR = ROOT / "revision" / "submission_ready_v4.7.0_20260902"
HTML_DIR = ROOT / "tmp" / "fig_html"

EXPECTED_MASTER_SHA256 = "a77d58f3a7221693107527026a729c2d35989604e33179dea527e7c16c568ecc"
EXPECTED_COHORT_SHA256 = "c6fa6d882b00e8af60d12ba5e98afd75578927a76f490637064f046e00a8bf03"

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"

INK = "#1a1a1a"
MUTED = "#555555"
GRID = "#e3e3e3"
BORDER = "#2b2b2b"
# NPG (Nature) palette
BLUE = "#3C5488"        # established pre-cancer DM (principal)
SALMON = "#F39B7F"      # peri-diagnostic
RED = "#E64B35"         # time-zero emphasis
TEAL = "#00A087"        # DM >10y before
LAV = "#8491B4"         # post-cancer
GREY = "#999999"        # order unknown
SLATE = "#4D4D4D"       # no diabetes (reference)

FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def esc(text) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def load_data() -> dict:
    assert sha256(MASTER_PATH) == EXPECTED_MASTER_SHA256, "master hash mismatch"
    assert sha256(COHORT_PATH) == EXPECTED_COHORT_SHA256, "cohort hash mismatch"
    master = load_results()
    meta = master["meta"]
    assert meta.get("status") == "VALIDATED_V4_R461_SUBMISSION"
    assert int(meta["cohort_n"]) == 5123 and int(meta["n_analysis"]) == 4910 and int(meta["events_10y"]) == 2001

    principal = {row["term"]: row for row in master["A_principal_5cat"]}
    expected = {'established_pre_cancer_dm': (1.6868, 1.4355, 1.9822), 'peri_diagnostic': (1.6819, 1.2763, 2.2164), 'post_cancer_dm': (0.9409, 0.7595, 1.1656), 'dm_order_unknown': (1.052, 0.6078, 1.8208)}
    for term, (hr, lo, hi) in expected.items():
        row = principal[term]
        assert abs(row["HR"] - hr) < 0.01 and abs(row["ci_lo"] - lo) < 0.01 and abs(row["ci_hi"] - hi) < 0.01, term

    import pandas as pd
    risk = pd.DataFrame(master["G_absolute_risk_bootstrap"])
    est = risk[(risk.trajectory == "established_pre_cancer_dm") & (risk.horizon_months == 120)].iloc[0]
    nodm = risk[(risk.trajectory == "no_diabetes") & (risk.horizon_months == 120)].iloc[0]
    assert abs(est.mortality - 0.598) < 0.005 and abs(nodm.mortality - 0.449) < 0.005
    assert abs(est.abs_risk_diff_per_100 - 14.84) < 0.05
    return {
        "meta": meta,
        "principal": principal,
        "risk": {"est": est, "nodm": nodm},
        "timing": master["C_timing_distributions"],
        "burden": {r["variable"]: r for r in master["F_burden_comparison"]},
    }


# ----------------------------------------------------------------- svg atoms

def text(x, y, content, size=12.5, fill=INK, anchor="start", weight="400", style="", spacing=""):
    extra = f' font-style="{style}"' if style else ""
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    return (f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" '
            f'text-anchor="{anchor}" font-weight="{weight}"{extra}{ls}>{esc(content)}</text>')


def line(x1, y1, x2, y2, stroke=INK, w=1, dash="", cap="butt"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{w}" stroke-linecap="{cap}"{d}/>'


def rect(x, y, w, h, fill, stroke=BORDER, sw=1.2, rx=0):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" rx="{rx}"/>'


def circle(cx, cy, r, fill, stroke="none", sw=0):
    return f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'


def square(cx, cy, s, fill):
    return f'<rect x="{cx - s / 2}" y="{cy - s / 2}" width="{s}" height="{s}" fill="{fill}"/>'


def ci(x_lo, x_hi, x, y, color, lw=2.6, r=6.5, cap=4.5):
    parts = [
        line(x_lo, y, x_hi, y, color, lw),
        line(x_lo, y - cap, x_lo, y + cap, color, lw),
        line(x_hi, y - cap, x_hi, y + cap, color, lw),
        circle(x, y, r, color),
    ]
    return "".join(parts)


def panel_header(x, y, letter, title):
    return (text(x, y, letter, size=17, weight="700") +
            text(x + 22, y, title, size=15, weight="700"))


def svg_page(width, height, body) -> str:
    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><style>
  html, body {{ margin: 0; padding: 0; background: #ffffff; }}
  @page {{ size: {width}px {height}px; margin: 0; }}
  svg {{ display: block; font-family: {FONT}; }}
</style></head>
<body>
<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg"
     font-family="{FONT}">
{body}
</svg>
</body></html>"""


# ----------------------------------------------------------------- figure 1

def fig1(data) -> tuple[str, int, int]:
    Wp, Hp = 1120, 520
    meta = data["meta"]
    b = []

    b.append(panel_header(40, 44, "A", "Analytic cohort"))
    b.append(panel_header(540, 44, "B", "Time origin and survivor selection"))

    # ---- panel a: cohort flow (white boxes, exclusion notes right of arrows) ----
    bx, bw, bh = 60, 360, 80
    boxes = [
        (80, "Digestive cancer survivors", "NHIS 1997\u20132018, public-use LMF", f"n = {int(meta['cohort_n']):,}"),
        (232, "Positive follow-up", "public-use quarter timing", "n = 5,082"),
        (384, "Principal Cox sample", "10-year deaths = 2,001", f"n = {int(meta['n_analysis']):,}"),
    ]
    arrow_x = bx + bw / 2
    note_lines = [
        (160, 232, "41 records with 0-month", "public-use follow-up excluded"),
        (312, 384, "172 positive-follow-up records", "with incomplete covariates"),
    ]
    for y, title, sub, n in boxes:
        b.append(rect(bx, y, bw, bh, "#ffffff", BORDER, 1.4))
        b.append(text(bx + 16, y + 31, title, size=14.5, weight="700"))
        b.append(text(bx + 16, y + 53, sub, size=12, fill=MUTED))
        b.append(text(bx + bw - 16, y + 53, n, size=13, fill=INK, anchor="end", weight="700"))
    for y1, y2, l1, l2 in note_lines:
        mid = (y1 + y2) / 2
        b.append(line(arrow_x, y1 + 5, arrow_x, y2 - 13, SLATE, 1.6))
        b.append(f'<path d="M {arrow_x - 5} {y2 - 13} L {arrow_x + 5} {y2 - 13} L {arrow_x} {y2 - 3} Z" fill="{SLATE}"/>')
        b.append(text(arrow_x + 16, mid - 4, l1, size=11.5, fill=MUTED))
        b.append(text(arrow_x + 16, mid + 13, l2, size=11.5, fill=MUTED))

    # ---- panel b: trajectory timeline, years -12 .. +5 ----
    px0, px1 = 540, 1070          # plot x range (same visual weight as panel a)
    py0, py1 = 120, 440           # plot y range
    span = 17.0                   # -12 .. +5

    def X(year):
        return px0 + (year + 12) / span * (px1 - px0)

    # minimal marker legend (replaces the former footnote strip)
    ly, lx = 84, 542
    b.append(circle(lx, ly - 4, 5, BLUE))
    b.append(text(lx + 12, ly, "Diabetes diagnosis", size=12, fill=MUTED))
    b.append(square(lx + 140, ly - 4, 9, "#1F2937"))
    b.append(text(lx + 152, ly, "Cancer diagnosis", size=12, fill=MUTED))
    b.append(line(lx + 262, ly - 4, lx + 296, ly - 4, GREY, 2))
    b.append(f'<path d="M {lx + 296} {ly - 8.5} L {lx + 296} {ly + 0.5} L {lx + 308} {ly - 4} Z" fill="{GREY}"/>')
    b.append(text(lx + 314, ly, "Follow-up", size=12, fill=MUTED))

    # grid + axes
    for yr in (-10, -5, 0, 5):
        b.append(line(X(yr), py0, X(yr), py1, GRID, 1))
    b.append(line(px0, py1, px1, py1, INK, 1.4))
    for yr, lab in [(-10, "\u221210"), (-5, "\u22125"), (0, "Interview"), (5, "+5")]:
        b.append(line(X(yr), py1, X(yr), py1 + 6, INK, 1.2))
        b.append(text(X(yr), py1 + 24, lab, size=12, fill=MUTED, anchor="middle"))
    b.append(text((px0 + px1) / 2, py1 + 50, "Years relative to NHIS interview", size=12.5, fill=INK, anchor="middle"))
    # time-zero emphasis
    b.append(line(X(0), py0 - 4, X(0), py1, RED, 2))
    b.append(text(X(0) + 8, py0 - 12, "time zero", size=12, fill=RED, weight="700"))

    # (label, color, dm_year, cancer_year, dm_label_dy)
    rows = [
        ("Established pre-cancer DM", BLUE, -8.0, -3.0, 24),
        ("Peri-diagnostic DM", SALMON, -3.2, -3.0, -24),
        ("No diabetes", SLATE, None, -3.0, 0),
        ("Post-cancer DM", LAV, -4.0, -8.0, -24),
    ]
    rh = (py1 - py0) / 4
    for i, (label, color, dm_x, ca_x, dm_dy) in enumerate(rows):
        y = py0 + rh * (i + 0.55)
        b.append(line(px0, y, px1, y, "#dcdcdc", 1))
        b.append(text(px0 + 2, y - 20, label, size=12.5, fill=color, weight="700"))
        # follow-up arrow from time zero (truncated at +5y)
        b.append(line(X(0), y, px1 - 10, y, "#c9c9c9", 1.6))
        b.append(f'<path d="M {px1 - 10} {y - 4.5} L {px1 - 10} {y + 4.5} L {px1} {y} Z" fill="#c9c9c9"/>')
        # cancer marker (label always below)
        b.append(square(X(ca_x), y, 10, "#1F2937"))
        b.append(text(X(ca_x), y + 24, "Cancer", size=10.5, fill=MUTED, anchor="middle"))
        if dm_x is not None:
            b.append(circle(X(dm_x), y, 6, color))
            b.append(text(X(dm_x), y + dm_dy, "DM", size=10.5, fill=color, anchor="middle"))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- figure 2

def fig2(data) -> tuple[str, int, int]:
    Wp, Hp = 1180, 420
    b = []
    b.append(panel_header(40, 44, "A", "Adjusted all-cause mortality"))
    b.append(panel_header(740, 44, "B", "Standardized 10-year mortality"))

    import math
    ay0, ay1 = 100, 350       # shared plot height

    # ---- panel a: forest (log scale), plot 280 wide ----
    ax0, ax1 = 280, 560
    vmin, vmax = 0.55, 2.55

    def X(v):
        return ax0 + (math.log(v) - math.log(vmin)) / (math.log(vmax) - math.log(vmin)) * (ax1 - ax0)

    for tv in (0.6, 1, 1.5, 2, 2.5):
        b.append(line(X(tv), ay0, X(tv), ay1, GRID, 1))
    b.append(line(X(1), ay0 - 6, X(1), ay1, "#9a9a9a", 1.4, dash="5,4"))
    b.append(line(ax0, ay1, ax1, ay1, INK, 1.4))
    for tv in (0.6, 1, 1.5, 2, 2.5):
        b.append(line(X(tv), ay1, X(tv), ay1 + 6, INK, 1.2))
        b.append(text(X(tv), ay1 + 24, f"{tv:g}", size=12, fill=MUTED, anchor="middle"))
    b.append(text((ax0 + ax1) / 2, ay1 + 50, "Hazard ratio (95% CI), reference = no diabetes",
                  size=12.5, fill=INK, anchor="middle"))

    p = data["principal"]
    rows = [
        ("Established pre-cancer DM", p["established_pre_cancer_dm"], BLUE, 122),
        ("Peri-diagnostic DM", p["peri_diagnostic"], SALMON, 186),
        None,
        ("Post-cancer DM", p["post_cancer_dm"], LAV, 276),
        ("Order unknown", p["dm_order_unknown"], GREY, 328),
    ]
    for row in rows:
        if row is None:
            b.append(line(ax0, 232, ax1, 232, "#c9c9c9", 1))
            b.append(text(ax1, 252, "conditional descriptive", size=11.5, fill=MUTED,
                          anchor="end", style="italic"))
            continue
        label, r, color, y = row
        b.append(text(ax0 - 16, y + 4, label, size=12.5, fill=INK, anchor="end"))
        b.append(ci(X(r["ci_lo"]), X(r["ci_hi"]), X(r["HR"]), y, color))
        b.append(text(700, y + 4, f"{r['HR']:.2f} ({r['ci_lo']:.2f}\u2013{r['ci_hi']:.2f})",
                      size=12.5, fill=MUTED, anchor="end"))

    # ---- panel b: standardized mortality, same plot height and width ----
    bx0, bx1 = 872, 1152
    lo, hi = 38, 68

    def BX(v):
        return bx0 + (v - lo) / (hi - lo) * (bx1 - bx0)

    for tv in (40, 45, 50, 55, 60, 65):
        b.append(line(BX(tv), ay0, BX(tv), ay1, GRID, 1))
    b.append(line(bx0, ay1, bx1, ay1, INK, 1.4))
    for tv in (40, 45, 50, 55, 60, 65):
        b.append(line(BX(tv), ay1, BX(tv), ay1 + 6, INK, 1.2))
        b.append(text(BX(tv), ay1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    b.append(text((bx0 + bx1) / 2, ay1 + 50, "Mortality probability (%)", size=12.5, fill=INK, anchor="middle"))

    risk = data["risk"]
    b.append(text(bx0 - 16, 154, "No diabetes", size=12.5, fill=INK, anchor="end"))
    b.append(text(bx0 - 16, 302, "Established", size=12.5, fill=INK, anchor="end"))
    b.append(text(bx0 - 16, 318, "pre-cancer DM", size=12.5, fill=INK, anchor="end"))
    rows_b = [(risk["nodm"], SLATE, 150), (risk["est"], BLUE, 310)]
    for r, color, y in rows_b:
        m = float(r.mortality) * 100
        lo_v, hi_v = float(r.mortality_ci_lo) * 100, float(r.mortality_ci_hi) * 100
        b.append(ci(BX(lo_v), BX(hi_v), BX(m), y, color))
        b.append(text(BX(hi_v) + 10, y + 4, f"{m:.1f}%", size=13, fill=color, weight="700"))

    est = risk["est"]
    b.append(text(bx1, 330, f"Risk difference +{float(est.abs_risk_diff_per_100):.2f} per 100",
                  size=12, fill=BLUE, weight="700", anchor="end"))
    b.append(text(bx1, 346, f"95% CI {float(est.abs_risk_diff_per_100_ci_lo):.2f} to {float(est.abs_risk_diff_per_100_ci_hi):.2f}",
                  size=11.5, fill=MUTED, anchor="end"))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- figure 3

def fig3(data) -> tuple[str, int, int]:
    Wp, Hp = 1260, 560
    b = []
    b.append(panel_header(40, 44, "A", "Reported diagnostic intervals"))
    b.append(panel_header(700, 44, "B", "Standardized burden differences"))

    timing = {(r["trajectory"], r["variable"]): r for r in data["timing"]}
    ay0, ay1 = 100, 480       # shared plot height

    # ---- panel a: interval plot, 300 wide ----
    ax0, ax1 = 280, 580

    def X(v):
        return ax0 + v / 50 * (ax1 - ax0)

    for tv in (0, 10, 20, 30, 40, 50):
        b.append(line(X(tv), ay0, X(tv), ay1, GRID, 1))
    b.append(line(ax0, ay1, ax1, ay1, INK, 1.4))
    for tv in (0, 10, 20, 30, 40, 50):
        b.append(line(X(tv), ay1, X(tv), ay1 + 6, INK, 1.2))
        b.append(text(X(tv), ay1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    b.append(text((ax0 + ax1) / 2, ay1 + 48, "Years", size=12.5, fill=INK, anchor="middle"))

    groups = [
        ("Cancer-to-interview", [
            ("No diabetes", "gi_only", SLATE),
            ("DM 2\u201310 years before", "dm_to_gi_2_10y", BLUE),
            ("DM >10 years before", "dm_to_gi_gt10y", TEAL),
            ("Peri-diagnostic", "peri_diagnostic", SALMON),
            ("Post-cancer", "gi_to_dm", LAV),
        ], 134),
        ("Diabetes-to-interview", [
            ("DM 2\u201310 years before", "dm_to_gi_2_10y", BLUE),
            ("DM >10 years before", "dm_to_gi_gt10y", TEAL),
            ("Peri-diagnostic", "peri_diagnostic", SALMON),
            ("Post-cancer", "gi_to_dm", LAV),
        ], 346),
    ]
    for gtitle, rows, y_start in groups:
        b.append(text(ax0 + 2, y_start - 26, gtitle, size=11.5, fill=MUTED, weight="700"))
        for i, (label, traj, color) in enumerate(rows):
            y = y_start + i * 34
            r = timing[(traj, "cancer_to_interview" if gtitle.startswith("Cancer") else "dm_to_interview")]
            b.append(text(ax0 - 16, y + 4, label, size=12.5, fill=INK, anchor="end"))
            b.append(line(X(r["p10"]), y, X(r["p90"]), y, color, 2))
            b.append(line(X(r["p25"]), y, X(r["p75"]), y, color, 6.5, cap="round"))
            b.append(circle(X(r["median"]), y, 5.5, "#ffffff", color, 2))
    b.append(line(ax0, 296, ax1, 296, "#c9c9c9", 1))

    # mini encoding legend (replaces the long axis note)
    gy = ay1 + 74
    b.append(line(292, gy - 4, 318, gy - 4, MUTED, 2))
    b.append(text(324, gy, "P10\u2013P90", size=11.5, fill=MUTED))
    b.append(line(404, gy - 4, 430, gy - 4, MUTED, 6.5, cap="round"))
    b.append(text(436, gy, "IQR", size=11.5, fill=MUTED))
    b.append(circle(488, gy - 4, 5, "#ffffff", MUTED, 2))
    b.append(text(498, gy, "median", size=11.5, fill=MUTED))

    # ---- panel b: burden differences, same plot height and width ----
    burden = data["burden"]
    bx0, bx1 = 830, 1130
    vmax = 30

    def BX2(v):
        return bx0 + v / vmax * (bx1 - bx0)

    for tv in (0, 10, 20, 30):
        b.append(line(BX2(tv), ay0, BX2(tv), ay1, GRID, 1))
    b.append(line(BX2(0), ay0 - 6, BX2(0), ay1, "#9a9a9a", 1.4, dash="5,4"))
    b.append(line(bx0, ay1, bx1, ay1, INK, 1.4))
    for tv in (0, 10, 20, 30):
        b.append(line(BX2(tv), ay1, BX2(tv), ay1 + 6, INK, 1.2))
        b.append(text(BX2(tv), ay1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    b.append(text((bx0 + bx1) / 2, ay1 + 48, "Established minus no diabetes (percentage points)",
                  size=12, fill=INK, anchor="middle"))

    rows_b = [
        ("Hypertension", "hypertension_b"),
        ("Obesity", "obesity"),
        ("Coronary heart disease", "chd_b"),
        ("Stroke", "stroke_b"),
        ("Fair/poor health", "srh_fairpoor"),
    ]
    LIGHTBLUE = "#2F6F9F"  # reference-image steel blue (per owner 2026-09-03)
    for i, (label, var) in enumerate(rows_b):
        y = 134 + i * 72
        r = burden[var]
        d, dlo, dhi = r["diff"] * 100, r["diff_ci_lo"] * 100, r["diff_ci_hi"] * 100
        b.append(text(bx0 - 16, y + 4, label, size=12.5, fill=INK, anchor="end"))
        b.append(ci(BX2(dlo), BX2(dhi), BX2(d), y, LIGHTBLUE))
        b.append(text(1245, y + 4, f"+{d:.1f} ({dlo:.1f}\u2013{dhi:.1f})", size=12, fill=MUTED, anchor="end"))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- render

def render(html_path: Path, stem: str, width: int, height: int) -> None:
    png = OUT_DIR / f"{stem}.png"
    pdf = OUT_DIR / f"{stem}.pdf"
    url = html_path.resolve().as_uri()
    subprocess.run([
        CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars",
        f"--window-size={width},{height}", "--force-device-scale-factor=3",
        f"--screenshot={png}", url,
    ], check=True, capture_output=True, timeout=120)
    subprocess.run([
        CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
        f"--print-to-pdf={pdf}", url,
    ], check=True, capture_output=True, timeout=120)
    print(f"{stem}: png={png.stat().st_size} pdf={pdf.stat().st_size}")


def main() -> None:
    render_only = "--render-only" in sys.argv
    data = None if render_only else load_data()
    HTML_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    specs = []
    if render_only:
        specs = [("fig1_rev1", 1240, 520), ("fig2_rev1", 1180, 420), ("fig3_rev1", 1260, 560)]
    else:
        for stem, maker in [("fig1_rev1", fig1), ("fig2_rev1", fig2), ("fig3_rev1", fig3)]:
            html, w, h = maker(data)
            path = HTML_DIR / f"{stem}.html"
            path.write_text(html, encoding="utf-8")
            specs.append((stem, w, h))

    for stem, w, h in specs:
        render(HTML_DIR / f"{stem}.html", stem, w, h)
    print("FRONTEND FIGURES OK")


if __name__ == "__main__":
    main()
