"""Supplementary figures sfig1/sfig2 via frontend construction (HTML/SVG -> Chrome headless).

Same visual language as rev2_figures_frontend.py (NPG palette, hairline grids,
dashed reference lines, right-aligned value columns, compact panels).

  sfig1  Two-panel forest on a unified log-HR scale:
         A = established pre-cancer DM split by diagnostic interval (2-10y vs >10y)
             with the Wald heterogeneity P value;
         B = sensitivity analyses for the established pre-cancer DM contrast.
  sfig2  Six-state trajectory distribution of the full cohort (unweighted n and
         survey-weighted percentage), horizontal bars.

Outputs: sfigN_rev1.html (tmp/fig_html), sfigN_rev1.png (scale 3) and
sfigN_rev1.pdf (vector) into revision/submission_ready_v4.7.0_20260902/.
"""

from __future__ import annotations

import hashlib
import json
from rev5_domain_results import load_results
import math
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
    assert int(meta["cohort_n"]) == 5123

    lag = {row["term"]: row for row in master["B_lag_6cat"]}
    assert abs(lag["dm_to_gi_2_10y"]["HR"] - 1.6395) < 0.001
    assert abs(lag["dm_to_gi_gt10y"]["HR"] - 1.7421) < 0.001
    wald = master["B_wald"][0]
    assert abs(wald["p_heterogeneity"] - 0.6708) < 0.001

    principal = {row["term"]: row for row in master["A_principal_5cat"]}["established_pre_cancer_dm"]
    assert abs(principal["HR"] - 1.687) < 0.01

    def est_row(key):
        for row in master[key]:
            if row.get("term") == "established_pre_cancer_dm":
                return row
        raise KeyError(key)

    sens = [
        ("Principal model", principal, BLUE),
        ("+ Cancer-to-interview interval", est_row("C_plus_cancer_to_interview"), SLATE),
        ("\u2212 Type-1-proxy diagnoses", est_row("D_t1d_proxy"), SLATE),
        ("\u2212 Income adjustment", est_row("A2_no_income_sensitivity"), SLATE),
        # Domain-correct zero-month sensitivity results
        ("+ Zero-month records (0.5/1-month)", est_row("G_zero_time_time_z05"), SLATE),
        ("+ Zero-month records (3-month)", est_row("G_zero_time_time_z3"), SLATE),
    ]

    # Six-state distribution from the full 22-year cohort.
    import pandas as pd
    cohort = pd.read_parquet(COHORT_PATH, columns=["trajectory_6cat", "sa_wgt_pool"])
    summary = cohort.groupby("trajectory_6cat", observed=True).sa_wgt_pool.agg(["size", "sum"])
    df = [(key, row["size"], row["sum"]) for key, row in summary.iterrows()]
    total_n = sum(r[1] for r in df)
    total_w = sum(r[2] for r in df)
    assert total_n == 5123, total_n
    dist = {r[0]: (int(r[1]), 100.0 * r[2] / total_w) for r in df}
    expected = {'gi_only': (4034, 79.164), 'peri_diagnostic': (151, 2.881), 'dm_to_gi_2_10y': (305, 5.952), 'dm_to_gi_gt10y': (275, 5.018), 'gi_to_dm': (315, 5.976), 'dm_order_unknown': (43, 1.009)}
    for k, (n, pct) in expected.items():
        assert dist[k][0] == n, (k, dist[k])
        assert abs(dist[k][1] - pct) < 0.05, (k, dist[k])

    return {"lag": lag, "wald": wald, "sens": sens, "dist": dist}


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


# ----------------------------------------------------------------- sfig 1

def sfig1(data) -> tuple[str, int, int]:
    Wp, Hp = 1180, 420
    b = []
    b.append(panel_header(40, 44, "A", "Pre-cancer DM by diagnostic interval"))
    b.append(panel_header(630, 44, "B", "Sensitivity analyses, established pre-cancer DM"))

    ay0, ay1 = 96, 340       # shared plot height
    vmin, vmax = 0.9, 2.3    # unified log-HR scale
    ticks = (1.0, 1.25, 1.5, 1.75, 2.0)

    def draw_forest(px0, px1, val_x):
        def X(v):
            return px0 + (math.log(v) - math.log(vmin)) / (math.log(vmax) - math.log(vmin)) * (px1 - px0)
        for tv in ticks:
            b.append(line(X(tv), ay0, X(tv), ay1, GRID, 1))
        b.append(line(X(1), ay0 - 6, X(1), ay1, "#9a9a9a", 1.4, dash="5,4"))
        b.append(line(px0, ay1, px1, ay1, INK, 1.4))
        for tv in ticks:
            b.append(line(X(tv), ay1, X(tv), ay1 + 6, INK, 1.2))
            b.append(text(X(tv), ay1 + 24, f"{tv:g}", size=12, fill=MUTED, anchor="middle"))
        b.append(text((px0 + px1) / 2, ay1 + 50, "Hazard ratio (95% CI), reference = no diabetes",
                      size=12.5, fill=INK, anchor="middle"))
        return X

    # ---- panel a (plot width identical to panel b; bottom row kept near axis) ----
    XA = draw_forest(270, 500, 596)
    lag, wald = data["lag"], data["wald"]
    rows_a = [
        ("DM 2\u201310 years before", lag["dm_to_gi_2_10y"], BLUE, 190),
        ("DM >10 years before", lag["dm_to_gi_gt10y"], TEAL, 310),
    ]
    for label, r, color, y in rows_a:
        b.append(text(254, y + 4, label, size=12.5, fill=INK, anchor="end"))
        b.append(ci(XA(r["ci_lo"]), XA(r["ci_hi"]), XA(r["HR"]), y, color))
        b.append(text(596, y + 4, f"{r['HR']:.2f} ({r['ci_lo']:.2f}\u2013{r['ci_hi']:.2f})",
                      size=12.5, fill=MUTED, anchor="end"))
    b.append(text(500, 120, f"P for heterogeneity = {wald['p_heterogeneity']:.3f}",
                  size=12, fill=MUTED, anchor="end", style="italic"))

    # ---- panel b (same plot width; principal flush-left, delta rows flush-right) ----
    XB = draw_forest(856, 1086, 1170)
    pitch = 40
    y0 = 112
    for i, (label, r, color) in enumerate(data["sens"]):
        y = y0 + i * pitch
        if i == 0:
            b.append(text(630, y + 4, label, size=12, fill=INK, anchor="start", weight="700"))
        else:
            b.append(text(844, y + 4, label, size=12, fill=INK, anchor="end"))
        b.append(ci(XB(r["ci_lo"]), XB(r["ci_hi"]), XB(r["HR"]), y, color, r=5.5, lw=2.2, cap=4))
        b.append(text(1170, y + 4, f"{r['HR']:.2f} ({r['ci_lo']:.2f}\u2013{r['ci_hi']:.2f})",
                      size=11.5, fill=MUTED, anchor="end"))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- sfig 2

def sfig2(data) -> tuple[str, int, int]:
    Wp, Hp = 1080, 400
    b = []
    dist = data["dist"]

    px0, px1 = 320, 900
    py0, py1 = 76, 320
    vmax = 85.0

    def X(v):
        return px0 + v / vmax * (px1 - px0)

    for tv in (0, 20, 40, 60, 80):
        b.append(line(X(tv), py0, X(tv), py1, GRID, 1))
    b.append(line(px0, py1, px1, py1, INK, 1.4))
    for tv in (0, 20, 40, 60, 80):
        b.append(line(X(tv), py1, X(tv), py1 + 6, INK, 1.2))
        b.append(text(X(tv), py1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    b.append(text((px0 + px1) / 2, py1 + 50,
                  "Survey-weighted share of cohort (%)", size=12.5, fill=INK, anchor="middle"))

    rows = [
        ("No diabetes", "gi_only", SLATE),
        ("DM 2\u201310 years before", "dm_to_gi_2_10y", BLUE),
        ("DM >10 years before", "dm_to_gi_gt10y", TEAL),
        ("Peri-diagnostic DM", "peri_diagnostic", SALMON),
        ("Post-cancer DM", "gi_to_dm", LAV),
        ("Order unknown", "dm_order_unknown", GREY),
    ]
    pitch = (py1 - py0) / len(rows)
    bh = 26
    for i, (label, key, color) in enumerate(rows):
        y = py0 + pitch * (i + 0.5)
        n, pct = dist[key]
        b.append(text(px0 - 16, y + 4, label, size=12.5, fill=INK, anchor="end"))
        b.append(rect(px0, y - bh / 2, max(X(pct) - px0, 1.5), bh, color, stroke="none", sw=0))
        b.append(text(X(pct) + 10, y + 4, f"n = {n:,} \u00b7 {pct:.1f}%", size=12, fill=MUTED))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- sfig 3

def path_d(points, close_y=None):
    d = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    d += [f"L {x:.1f} {y:.1f}" for x, y in points[1:]]
    return " ".join(d)


def step_d(points):
    """Step-post path (KM convention): horizontal hold, vertical jump at each x."""
    d = [f"M {points[0][0]:.1f} {points[0][1]:.1f}"]
    for x, y in points[1:]:
        d.append(f"H {x:.1f}")
        d.append(f"V {y:.1f}")
    return " ".join(d)


def sfig3(data) -> tuple[str, int, int]:
    """Model-standardized cumulative mortality, no diabetes vs established pre-cancer DM.

    Reads revision/sfig3_stdcurve/sfig3_curves.json produced by rev2_stdcurve.py
    (500-rep stratified PSU bootstrap; point QC replicated to locked values).
    """
    curves = json.loads((ROOT / "revision" / "sfig3_stdcurve" / "sfig3_curves.json")
                        .read_text(encoding="utf-8"))
    assert curves["n_boot"] == 500 and curves["boot_qc_pass"]
    assert curves["bootstrap_method"] == "rao_wu_nminus1_period_psu_replicate_weights_v1"
    grid = curves["grid"]
    assert len(grid) == 120 and grid[-1] == 120
    # locked point values
    assert abs(curves["no_diabetes"]["point"][119] - 0.449094060) < 1e-6
    assert abs(curves["established"]["point"][119] - 0.597508438) < 1e-6

    Wp, Hp = 980, 480
    b = []
    px0, px1 = 110, 900
    py0, py1 = 60, 400
    ymax = 70.0

    def X(month):
        return px0 + month / 120.0 * (px1 - px0)

    def Y(pct):
        return py1 - pct / ymax * (py1 - py0)

    for tv in range(0, 11, 2):
        b.append(line(X(tv * 12), py0, X(tv * 12), py1, GRID, 1))
    for tv in range(0, 71, 10):
        b.append(line(px0, Y(tv), px1, Y(tv), GRID, 1))
    b.append(line(px0, py1, px1, py1, INK, 1.4))
    b.append(line(px0, py0, px0, py1, INK, 1.4))
    for tv in range(0, 11, 2):
        b.append(line(X(tv * 12), py1, X(tv * 12), py1 + 6, INK, 1.2))
        b.append(text(X(tv * 12), py1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    for tv in range(0, 71, 10):
        b.append(line(px0 - 6, Y(tv), px0, Y(tv), INK, 1.2))
        b.append(text(px0 - 12, Y(tv) + 4, str(tv), size=12, fill=MUTED, anchor="end"))
    b.append(text((px0 + px1) / 2, py1 + 50, "Years since interview", size=12.5, fill=INK, anchor="middle"))
    b.append(f'<text x="{px0 - 66}" y="{(py0 + py1) / 2}" font-size="12.5" fill="{INK}" '
             f'text-anchor="middle" transform="rotate(-90 {px0 - 66} {(py0 + py1) / 2})">'
             f'Standardized cumulative mortality (%)</text>')

    def band(key, color):
        hi = [100 * v for v in curves[key]["hi"]]
        lo = [100 * v for v in curves[key]["lo"]]
        hi_pts = [(X(m), Y(v)) for m, v in zip(grid, hi)]
        lo_pts = [(X(m), Y(v)) for m, v in zip(grid, lo)]
        # upper edge: step-post forward; lower edge: step-consistent backward (V then H)
        d = step_d(hi_pts)
        d += f" L {lo_pts[-1][0]:.1f} {lo_pts[-1][1]:.1f}"
        for x, y in reversed(lo_pts[:-1]):
            d += f" V {y:.1f} H {x:.1f}"
        d += " Z"
        return f'<path d="{d}" fill="{color}" fill-opacity="0.18" stroke="none"/>'

    def curve(key, color):
        pt = [100 * v for v in curves[key]["point"]]
        pts = [(X(0), Y(0.0))] + [(X(m), Y(v)) for m, v in zip(grid, pt)]
        return f'<path d="{step_d(pts)}" fill="none" stroke="{color}" stroke-width="2.6"/>'

    b.append(band("no_diabetes", SLATE))
    b.append(band("established", BLUE))
    b.append(curve("no_diabetes", SLATE))
    b.append(curve("established", BLUE))

    # end labels
    m0 = 100 * curves["no_diabetes"]["point"][119]
    m1 = 100 * curves["established"]["point"][119]
    b.append(text(px1 + 12, Y(m0) + 4, f"{m0:.1f}%", size=12.5, fill=SLATE, weight="700"))
    b.append(text(px1 + 12, Y(m1) + 4, f"{m1:.1f}%", size=12.5, fill=BLUE, weight="700"))

    # in-plot legend (upper left, plot is empty there)
    lx, ly = px0 + 26, py0 + 34
    b.append(line(lx, ly - 4, lx + 34, ly - 4, SLATE, 2.6))
    b.append(text(lx + 42, ly, "No diabetes", size=12.5, fill=INK))
    b.append(line(lx, ly + 22, lx + 34, ly + 22, BLUE, 2.6))
    b.append(text(lx + 42, ly + 26, "Established pre-cancer DM", size=12.5, fill=INK))

    # RD annotation — report the LOCKED R bootstrap CI (manuscript text: 10.17-19.14);
    # the Python replication band is visually identical (QC-passed, diff ~1e-9 on points).
    master = load_results()
    locked = [r for r in master["G_absolute_risk_bootstrap"]
              if r.get("trajectory") == "established_pre_cancer_dm"
              and r.get("horizon_months") == 120]
    assert len(locked) == 1
    rd = locked[0]["abs_risk_diff_per_100"]
    rdc = [locked[0]["abs_risk_diff_per_100_ci_lo"], locked[0]["abs_risk_diff_per_100_ci_hi"]]
    assert abs(rd - 14.8414) < 1e-3 and abs(rdc[0] - 10.1701) < 1e-3 and abs(rdc[1] - 19.1380) < 1e-3
    b.append(text(px1 - 10, py1 - 64, f"Risk difference at 10 years: +{rd:.1f} per 100",
                  size=12, fill=BLUE, weight="700", anchor="end"))
    b.append(text(px1 - 10, py1 - 46, f"95% CI {rdc[0]:.1f} to {rdc[1]:.1f} (Rao–Wu bootstrap)",
                  size=11.5, fill=MUTED, anchor="end"))

    return svg_page(Wp, Hp, "".join(b)), Wp, Hp


# ----------------------------------------------------------------- sfig 4

def sfig4(data) -> tuple[str, int, int]:
    """Exploratory spline lag–mortality association (survey-weighted Cox, 3-knot RCS).

    Reads spline_lag_curve.csv / spline_lag_peri.csv / spline_lag_rug.csv /
    spline_lag_diagnostics.csv produced by scripts/rev1_p3_spline_lag.R
    (500-rep stratified PSU bootstrap pointwise CI).
    """
    import csv

    out = ROOT / "outputs" / "revision_round1_v4_r461_sens2007"

    def read_csv(name):
        with open(out / name, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))

    diag = {r["item"]: float(r["value"]) for r in read_csv("spline_lag_diagnostics.csv")}
    assert diag["n_analysis"] == 4910 and diag["events"] == 2001
    assert diag["boot_converged"] == 500
    assert abs(diag["p_nonlinearity"] - 0.6868) < 0.01
    assert [diag["knot3_1"], diag["knot3_2"], diag["knot3_3"]] == [3.0, 10.0, 29.0]

    curve = [(float(r["lag_years"]), float(r["hr"]), float(r["ci_lo"]), float(r["ci_hi"]))
             for r in read_csv("spline_lag_curve.csv")]
    peri = read_csv("spline_lag_peri.csv")[0]
    peri_hr, peri_lo, peri_hi = (float(peri["HR"]), float(peri["ci_lo_boot"]),
                                 float(peri["ci_hi_boot"]))
    rug = [float(r["lag_years"]) for r in read_csv("spline_lag_rug.csv")]
    lmax = curve[-1][0]
    assert lmax == 46.0

    Wp, Hp = 980, 470
    b = []
    px0, px1 = 110, 920              # continuous lag axis, starts at 2 y
    py0, py1 = 50, 386
    vmin = 0.8
    vmax = math.ceil(max(hi for _, _, _, hi in curve) * 1.2 * 2) / 2
    y_ticks = tuple(v for v in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5) if v < vmax)
    assert all(vmin <= lo <= hr <= hi <= vmax for _, hr, lo, hi in curve)

    def X(lag):
        return px0 + (lag - 2.0) / (lmax - 2.0) * (px1 - px0)

    def Y(hr):
        return py1 - (math.log(hr) - math.log(vmin)) / (math.log(vmax) - math.log(vmin)) * (py1 - py0)

    # light background band for the 2-10 y classification window
    b.append(f'<rect x="{X(2):.1f}" y="{py0}" width="{X(10) - X(2):.1f}" height="{py1 - py0}" '
             f'fill="#eef2f7" stroke="none"/>')

    # grids + axes
    for tv in y_ticks:
        b.append(line(px0, Y(tv), px1, Y(tv), GRID, 1))
    for tv in (2, 10, 20, 30, 40):
        b.append(line(X(tv), py0, X(tv), py1, GRID, 1))
    b.append(line(px0, py1, px1, py1, INK, 1.4))
    b.append(line(px0, py0, px0, py1, INK, 1.4))
    b.append(line(px0, Y(1), px1, Y(1), "#9a9a9a", 1.4, dash="5,4"))   # HR = 1
    b.append(line(X(10), py0, X(10), py1, SLATE, 1.2, dash="2,3"))     # 10-y boundary
    b.append(text(X(10) + 7, py0 + 12, "10 y", size=11, fill=MUTED))
    for tv in (2, 10, 20, 30, 40):
        b.append(line(X(tv), py1, X(tv), py1 + 6, INK, 1.2))
        b.append(text(X(tv), py1 + 24, str(tv), size=12, fill=MUTED, anchor="middle"))
    for tv in y_ticks:
        b.append(line(px0 - 6, Y(tv), px0, Y(tv), INK, 1.2))
        b.append(text(px0 - 12, Y(tv) + 4, f"{tv:g}", size=12, fill=MUTED, anchor="end"))
    b.append(text((px0 + px1) / 2, py1 + 50, "Reported diabetes-to-cancer lag (years)",
                  size=12.5, fill=INK, anchor="middle"))
    b.append(f'<text x="{px0 - 62}" y="{(py0 + py1) / 2}" font-size="12.5" fill="{INK}" '
             f'text-anchor="middle" transform="rotate(-90 {px0 - 62} {(py0 + py1) / 2})">'
             f'Adjusted hazard ratio vs no diabetes</text>')

    # peri-diagnostic estimate as a numeric annotation (not joined to the curve)
    PERI_BLUE = "#7C8EB1"            # lighter tint of the curve blue (#3C5488)
    b.append(text(X(2) + 10, py0 + 20, "Peri-diagnostic DM:", size=12, fill=PERI_BLUE, weight="700"))
    b.append(text(X(2) + 10, py0 + 37, f"HR {peri_hr:.2f} ({peri_lo:.2f}\u2013{peri_hi:.2f})",
                  size=12, fill=PERI_BLUE))

    # bootstrap CI band + spline curve
    hi_pts = [(X(l), Y(hi)) for l, _, _, hi in curve]
    lo_pts = [(X(l), Y(lo)) for l, _, lo, _ in curve]
    d = path_d(hi_pts) + f" L {lo_pts[-1][0]:.1f} {lo_pts[-1][1]:.1f} " + \
        " ".join(f"L {x:.1f} {y:.1f}" for x, y in reversed(lo_pts[:-1])) + " Z"
    b.append(f'<path d="{d}" fill="{BLUE}" fill-opacity="0.15" stroke="none"/>')
    b.append(f'<path d="{path_d([(X(l), Y(h)) for l, h, _, _ in curve])}" fill="none" '
             f'stroke="{BLUE}" stroke-width="2.6"/>')

    # unweighted rug of observed lag values (within display range)
    rug_in = [v for v in rug if v <= lmax]
    b.append(f'<g stroke="{BLUE}" stroke-width="1" opacity="0.28">')
    b.append("".join(line(X(v), py1, X(v), py1 - 7, BLUE, 1) for v in rug_in))
    b.append("</g>")

    b.append(text(px1, py0 + 14, f"P for nonlinearity = {diag['p_nonlinearity']:.3f}",
                  size=12, fill=MUTED, anchor="end", style="italic"))

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
    data = load_data()
    HTML_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    specs = []
    for stem, maker in [("sfig1_rev1", sfig1), ("sfig2_rev1", sfig2), ("sfig3_rev1", sfig3),
                        ("sfig4_rev1", sfig4)]:
        html, w, h = maker(data)
        path = HTML_DIR / f"{stem}.html"
        path.write_text(html, encoding="utf-8")
        specs.append((stem, w, h))

    for stem, w, h in specs:
        render(HTML_DIR / f"{stem}.html", stem, w, h)
    print("SFIGS 1-4 OK")


if __name__ == "__main__":
    main()
