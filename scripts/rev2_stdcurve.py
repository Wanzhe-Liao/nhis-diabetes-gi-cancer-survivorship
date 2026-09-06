"""sfig3 computation: model-standardized cumulative mortality curves.

Replicates the locked R point-model standardization in pure Python
(scipy/numpy), with a full-population domain Rao--Wu PSU bootstrap.  Gates:

  * weighted Cox with Efron ties, case weights = sa_wgt_pool
    QC: |HR(established) - 1.68683| < 0.005 (and the other three trajectory HRs)
  * Breslow baseline cumulative hazard (survival::basehaz centered=FALSE)
  * standardized mortality M_g(t) = weighted mean of 1 - exp(-H0(t) exp(lp_i^g))
    over the principal-model analysis sample (target population), g in
    {no_diabetes, established_pre_cancer_dm}
    QC: M_no_dm(120) ~ 0.449, M_est(120) ~ 0.598 (+-0.005), RD ~ 14.84/100
  * Rao--Wu n_h-1 bootstrap replicates: draw PSUs with replacement within
    pooled design-period x original-stratum groups.  Selected PSU rows are
    retained once and receive their draw multiplicity times n_h/(n_h-1) in
    the replicate weights.  PSUs are drawn from the full Sample Adult frame
    before subsetting to GI rows.

Usage:
  python scripts/rev2_stdcurve.py --qc-only            # point fit + gates + curves
  python scripts/rev2_stdcurve.py --boot --start 1 --end 100
  python scripts/rev2_stdcurve.py --finalize           # merge chunks -> bands + QC
"""

from __future__ import annotations

import argparse
import json
from rev5_domain_results import load_results, COHORT_PATH, FULLSAMPLE_PATH, POVERTY_PATH, SOURCE
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "revision" / "sfig3_stdcurve"
# Method-specific location prevents resuming the former year-stratified,
# unscaled PSU chunks after changing the bootstrap design.
CHUNK_DIR = ROOT / "tmp" / "sfig3_rao_wu_period_chunks"

BOOT_SEED = 20260829
BOOT_REPS = 500
GRID = np.arange(1, 121)  # months 1..120

# Reference estimates come from the current R primary model and bootstrap export.
RISK_REFERENCE = pd.read_csv(SOURCE / "rev1_absolute_risk_with_rd.csv")
LOCKED_MORT = {g: {int(r.horizon_months): float(r.mortality) for r in RISK_REFERENCE.loc[RISK_REFERENCE.trajectory.eq(g)].itertuples()}
               for g in ["no_diabetes", "established_pre_cancer_dm"]}
_rd = RISK_REFERENCE.loc[RISK_REFERENCE.trajectory.eq("established_pre_cancer_dm") & RISK_REFERENCE.horizon_months.eq(120)].iloc[0]
LOCKED_RD_120 = float(_rd["abs_risk_diff_per_100"])

DESIGN_PERIODS = ((1997, 2005, "1997-2005"),
                  (2006, 2015, "2006-2015"),
                  (2016, 2018, "2016-2018"))
BOOT_METHOD_ID = "rao_wu_nminus1_period_psu_replicate_weights_v1"

SITE_FLAGS = ["colon_flag", "esoph_flag", "gallbladder_flag", "liver_flag",
              "pancreas_flag", "rectum_flag", "stomach_flag"]


# ------------------------------------------------------------------ data prep

def design_period_id(year: pd.Series | np.ndarray) -> pd.Series:
    """Return the pooled NHIS design-period ID for each survey year.

    The original design strata and PSU identifiers are reused across survey
    years.  Pooling them within the three design periods preserves that
    identifier continuity while keeping the 2006 and 2016 redesign boundaries.
    """
    years = pd.Series(year)
    out = pd.Series(pd.NA, index=years.index, dtype="string")
    years_num = pd.to_numeric(years, errors="raise")
    for lo, hi, label in DESIGN_PERIODS:
        out.loc[years_num.between(lo, hi, inclusive="both")] = label
    if out.isna().any():
        bad = sorted(pd.unique(years_num[out.isna()]).tolist())
        raise ValueError(f"Survey years outside the pooled design periods: {bad}")
    return out


def add_period_design_ids(frame: pd.DataFrame) -> pd.DataFrame:
    """Add period-level strata/PSU IDs using the original design fields."""
    frame = frame.copy()
    frame["design_period"] = design_period_id(frame["year"])
    strata = frame["design_strata"].astype(int).astype(str)
    psu = frame["design_psu"].astype(int).astype(str)
    frame["design_strata_prefixed"] = frame["design_period"] + "." + strata
    frame["design_psu_prefixed"] = (frame["design_strata_prefixed"] + "." + psu)
    return frame

def load_analysis_frame() -> pd.DataFrame:
    cols = ["trajectory_6cat", "followup_years", "mortstat", "design_psu",
            "design_strata", "sa_wgt_pool", "year", "hhx", "fmx", "age", "sex",
            "race", "region", "bmi", "smoking_3cat", "survey_year",
            "education_4cat"] + SITE_FLAGS
    df = pd.read_parquet(COHORT_PATH, columns=cols)

    pov = pd.read_csv(POVERTY_PATH, dtype={"hhx": str, "fmx": str})
    df["hhx_k"] = df["hhx"].astype(float).astype(int).astype(str)
    df["fmx_k"] = df["fmx"].astype(float).astype(int).astype(str)
    pov["hhx_k"] = pov["hhx"].astype(float).astype(int).astype(str)
    pov["fmx_k"] = pov["fmx"].astype(float).astype(int).astype(str)
    pov["year"] = pov["year"].astype(int)
    df["year"] = df["year"].astype(int)
    n0 = len(df)
    df = df.merge(pov[["year", "hhx_k", "fmx_k", "poverty_3cat_rep"]],
                  on=["year", "hhx_k", "fmx_k"], how="left", validate="m:1")
    assert len(df) == n0, "poverty join changed row count"

    # survival variables (round4 rule, identical to rev1_p1_main.R)
    df["time_months"] = np.minimum(df["followup_years"] * 12.0, 120.0)
    df["event"] = ((df["mortstat"] == 1) & (df["followup_years"] <= 10)).astype(float)
    df = df[df["time_months"].notna() & (df["time_months"] > 0)]
    df = df[df["design_psu"].notna() & df["design_strata"].notna() & df["sa_wgt_pool"].notna()]

    # five-state exposure
    traj5 = df["trajectory_6cat"].astype(str).replace({
        "dm_to_gi_2_10y": "established_pre_cancer_dm",
        "dm_to_gi_gt10y": "established_pre_cancer_dm",
        "gi_to_dm": "post_cancer_dm",
        "gi_only": "no_diabetes",
    })
    df["trajectory_5cat_rev"] = traj5

    # R make_factor: "" -> NA; poverty NA -> explicit "missing" level
    df["poverty_3cat_rep"] = df["poverty_3cat_rep"].fillna("missing")

    # model-frame complete cases (poverty excluded: missing is a level)
    mask = (
        df["bmi"].notna()
        & df["smoking_3cat"].notna() & (df["smoking_3cat"] != "")
        & df["education_4cat"].notna() & (df["education_4cat"] != "")
        & df["age"].notna() & df["sex"].notna() & df["race"].notna()
        & df["region"].notna() & df["survey_year"].notna()
        & df["trajectory_5cat_rev"].notna()
    )
    for v in SITE_FLAGS:
        mask &= df[v].notna()
    df = df[mask].reset_index(drop=True)

    n, ev = len(df), int(df["event"].sum())
    assert n == 4910 and ev == 2001, (n, ev)

    return add_period_design_ids(df)


# ------------------------------------------------------------- design matrix

TRAJ_LEVELS = ["no_diabetes", "established_pre_cancer_dm", "peri_diagnostic",
               "post_cancer_dm", "dm_order_unknown"]  # ref = no_diabetes
CAT_SPEC = {
    "sex": ("1", ["2"]),
    "race": ("1", ["2", "97"]),
    "region": ("1", ["2", "3", "4"]),
    "smoking_3cat": ("current", ["former", "never"]),
    "education_4cat": ("college_grad", ["high_school", "lt_high_school", "some_college"]),
    "poverty_3cat_rep": ("2_0_to_3_99", ["ge_4_0", "lt_2_0", "missing"]),
}
SURVEY_YEARS = list(range(1997, 2019))


def build_design(df: pd.DataFrame) -> tuple[np.ndarray, list[str], dict[str, slice]]:
    """Fixed-column design matrix; returns X, names, trajectory column slice."""
    n = len(df)
    blocks, names = [], []

    def add(mat, cols):
        blocks.append(np.asarray(mat, dtype=float).reshape(n, -1))
        names.extend(cols)

    add((df["age"].to_numpy() - 70.0) / 12.0, ["age"])
    for var, (ref, dums) in CAT_SPEC.items():
        s = df[var].astype("Int64").astype(str) if var in ("sex", "race", "region") else df[var].astype(str)
        for lev in dums:
            add((s == lev).to_numpy().astype(float), [f"{var}{lev}"])
    add((df["bmi"].to_numpy() - 27.0) / 5.0, ["bmi"])
    for yr in SURVEY_YEARS[1:]:
        add((df["survey_year"] == yr).to_numpy().astype(float), [f"survey_year{yr}"])
    for v in SITE_FLAGS:
        add(df[v].astype(float).to_numpy(), [v])
    # trajectory dummies last, keep slice for counterfactual replacement
    t0 = len(names)
    for lev in TRAJ_LEVELS[1:]:
        add((df["trajectory_5cat_rev"] == lev).to_numpy().astype(float),
            [f"trajectory_5cat_rev{lev}"])
    X = np.hstack(blocks)
    return X, names, slice(t0, len(names))


# ------------------------------------------------- weighted Cox (Efron ties)

def make_groups(time: np.ndarray, event: np.ndarray):
    """Sort descending by time; per unique event time -> (risk_end, death_idx)."""
    order = np.argsort(-time, kind="stable")
    t_s, e_s = time[order], event[order]
    groups = []
    i, n = 0, len(t_s)
    while i < n:
        j = i
        while j + 1 < n and t_s[j + 1] == t_s[i]:
            j += 1
        d = np.where(e_s[i:j + 1] > 0)[0]
        if len(d):
            groups.append((j, i + d))  # risk set = rows 0..j (time >= t)
        i = j + 1
    return order, groups


def cox_negloglik_grad(beta, X, w, groups, p):
    eta = np.exp(np.clip(X @ beta, -30, 30))
    weta = w * eta
    cA0 = np.cumsum(weta)
    cA1 = np.cumsum(weta[:, None] * X, axis=0)
    if not (np.all(np.isfinite(eta)) and np.all(np.isfinite(weta))
            and np.all(np.isfinite(cA0)) and np.all(np.isfinite(cA1))):
        return np.inf, np.full(p, np.nan)
    ll = 0.0
    g = np.zeros(p)
    for r_end, d_idx in groups:
        Xd = X[d_idx]
        wd = w[d_idx]
        nd = len(d_idx)
        A0 = cA0[r_end]
        A1 = cA1[r_end]
        ll += (wd * (Xd @ beta)).sum()
        g += (wd[:, None] * Xd).sum(axis=0)
        B0 = weta[d_idx].sum()
        B1 = (weta[d_idx, None] * Xd).sum(axis=0)
        fr = np.arange(nd) / nd
        den = A0 - fr * B0
        if not np.all(np.isfinite(den)) or np.any(den <= 0):
            return np.inf, np.full(p, np.nan)
        wsum = wd.sum()
        ll -= (wsum / nd) * np.log(den).sum()
        g -= (wsum / nd) * ((A1[None, :] - fr[:, None] * B1[None, :])
                            / den[:, None]).sum(axis=0)
    value, gradient = -ll, -g
    if not (np.isfinite(value) and np.all(np.isfinite(gradient))):
        return np.inf, np.full(p, np.nan)
    return value, gradient


def fit_cox(X, w, time, event, beta0=None):
    from scipy.optimize import minimize
    X = np.asarray(X, dtype=float)
    w = np.asarray(w, dtype=float)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    if X.ndim != 2 or len(X) != len(w) or len(time) != len(w) or len(event) != len(w):
        raise ValueError("Cox inputs have inconsistent lengths")
    if (len(w) == 0 or not np.all(np.isfinite(X)) or not np.all(np.isfinite(w))
            or not np.all(np.isfinite(time)) or not np.all(np.isfinite(event))
            or np.any(w <= 0)):
        raise ValueError("Cox inputs contain non-finite values or non-positive weights")
    w = w / w.mean()  # scale-invariant; keeps loglik/gradient magnitudes sane
    order, groups = make_groups(time, event)
    Xs, ws = X[order], w[order]
    p = X.shape[1]
    beta0 = np.zeros(p) if beta0 is None else beta0.copy()

    def fg(beta):
        return cox_negloglik_grad(beta, Xs, ws, groups, p)

    res = minimize(fg, beta0, method="BFGS", jac=True,
                   options={"gtol": 1e-8, "maxiter": 500})
    if (not np.all(np.isfinite(res.x)) or not np.isfinite(res.fun)
            or not np.all(np.isfinite(res.jac))):
        raise RuntimeError("non-finite Cox fit result")
    return res.x, res


def breslow_h0(beta, X, w, time, event):
    """Baseline cumulative hazard at unique event times.

    Mirrors survival::basehaz(fit, centered=FALSE) for an Efron-ties coxph fit:
    survfit.coxph defaults ctype=2 for method='efron', so the increment at a
    tied event time is (D0w/nd) * sum_j 1/(S0 - (j/nd)*B0r)  [agsurv5], where
    S0 = sum(w*eta) at risk, B0r = sum(w*eta) of deaths, D0w = sum(w) of deaths.
    Weight-scale invariant.  Returns ascending (times, cumulative H0).
    """
    w = w / w.mean()
    order, groups = make_groups(time, event)
    Xs, ws = X[order], w[order]
    eta = np.exp(np.clip(Xs @ beta, -30, 30))
    weta = ws * eta
    cA0 = np.cumsum(weta)
    t_sorted = time[order]
    pairs = []
    for r_end, d_idx in groups:
        S0 = cA0[r_end]
        wd = ws[d_idx]                    # sorted-order weights (d_idx are sorted positions)
        nd = len(d_idx)
        B0r = weta[d_idx].sum()
        fr = np.arange(nd) / nd
        den = S0 - fr * B0r
        if not np.all(np.isfinite(den)) or np.any(den <= 0):
            raise RuntimeError("non-positive or non-finite Efron risk denominator")
        inc = (wd.sum() / nd) * (1.0 / den).sum()
        if not np.isfinite(inc):
            raise RuntimeError("non-finite Breslow hazard increment")
        pairs.append((t_sorted[r_end], inc))
    pairs.sort(key=lambda p: p[0])
    times = np.array([p[0] for p in pairs])
    h0 = np.cumsum([p[1] for p in pairs])
    if not (np.all(np.isfinite(times)) and np.all(np.isfinite(h0))):
        raise RuntimeError("non-finite Breslow baseline hazard")
    return times, h0


def h0_at(times, h0, grid):
    idx = np.searchsorted(times, grid, side="right") - 1
    out = np.zeros(len(grid))
    ok = idx >= 0
    out[ok] = h0[idx[ok]]
    return out


def standardized_mortality(beta, X, w, traj_slice, h0_times, h0, grid):
    """M_g(t) for g in {no_diabetes, established} over the given target rows."""
    lp_base = X @ beta
    out = {}
    wsum = w.sum()
    if not (np.isfinite(wsum) and wsum > 0):
        raise RuntimeError("non-finite or non-positive standardization weight sum")
    H = h0_at(h0_times, h0, grid)
    if not np.all(np.isfinite(H)):
        raise RuntimeError("non-finite baseline hazard on standardization grid")
    for g, cols in [("no_diabetes", None), ("established_pre_cancer_dm", 0)]:
        lp = lp_base.copy()
        if cols is None:
            lp = lp - X[:, traj_slice] @ beta[traj_slice]  # all dummies to 0
        else:
            lp = lp - X[:, traj_slice] @ beta[traj_slice] + beta[traj_slice][0]
        eta = np.exp(np.clip(lp, -30, 30))
        surv = np.exp(-np.clip(H[:, None] * eta[None, :], 0, 50))
        out[g] = ((1 - surv) * w[None, :]).sum(axis=1) / wsum
        if not np.all(np.isfinite(out[g])):
            raise RuntimeError(f"non-finite standardized mortality curve: {g}")
    return out


# ------------------------------------------------------------------ bootstrap

def prepare_domain_sampler(df):
    """Build the full-frame pooled-period PSU sampler.

    The frame contains every eligible Sample Adult PSU, including PSUs with no
    analytic GI rows.  Strata and PSU IDs use the original NHIS identifiers
    within each design period, so a PSU's GI rows from every year in that
    period are kept together.
    """
    full = pd.read_parquet(FULLSAMPLE_PATH,
                           columns=["year", "design_strata", "design_psu"])
    assert len(full) == 646201
    full = add_period_design_ids(full)
    full["stratum"] = full["design_strata_prefixed"]
    full["psu"] = full["design_psu_prefixed"]
    frame = full[["stratum", "psu"]].drop_duplicates()
    rows = df.groupby("design_psu_prefixed", sort=False).indices
    assert set(rows).issubset(set(frame.psu))
    sampler = []
    for stratum, group in frame.groupby("stratum", sort=True):
        psus = [rows.get(p, np.array([], dtype=int)) for p in group.psu]
        if len(psus) < 2:
            raise ValueError(
                f"Unexpected singleton full-frame stratum {stratum!r}; "
                "Rao-Wu n_h-1 bootstrap requires n_h >= 2"
            )
        sampler.append(psus)
    print(f"Domain bootstrap frame: {len(full):,} adults, {len(sampler):,} strata, {len(frame):,} PSUs", flush=True)
    return sampler


def stratified_psu_sample(sampler, rng):
    """Draw Rao--Wu n_h-1 PSUs per stratum with replacement.

    Returns analysis-row indices and the corresponding replicate-weight
    multipliers.  A selected PSU's rows are represented once with its draw
    multiplicity; this is the weighted-row representation used by the shared
    R bootstrap and keeps Efron tie handling identical to that fit.  Empty
    domain PSUs remain valid frame members and simply contribute no rows.
    """
    idx_chunks = []
    multiplier_chunks = []
    for h, psus in enumerate(sampler):
        n_h = len(psus)
        if n_h < 2:
            raise ValueError(
                f"Unexpected singleton bootstrap stratum at position {h}; "
                "Rao-Wu n_h-1 bootstrap requires n_h >= 2"
            )
        draws = rng.integers(0, n_h, size=n_h - 1)
        counts = np.bincount(draws, minlength=n_h)
        scale = n_h / (n_h - 1)
        for j, multiplicity in enumerate(counts):
            if multiplicity == 0:
                continue
            rows = np.asarray(psus[j], dtype=int)
            if rows.size == 0:
                continue
            idx_chunks.append(rows)
            multiplier_chunks.append(
                np.full(rows.size, multiplicity * scale, dtype=float)
            )
    if not idx_chunks:
        return np.empty(0, dtype=int), np.empty(0, dtype=float)
    return np.concatenate(idx_chunks), np.concatenate(multiplier_chunks)


def run_bootstrap(df, X, traj_slice, beta_point, start, end):
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    w = df["sa_wgt_pool"].to_numpy()
    time_v = df["time_months"].to_numpy()
    event_v = df["event"].to_numpy()
    path = CHUNK_DIR / f"boot_{start:03d}_{end:03d}.npz"
    results = {}
    if path.exists():
        z = np.load(path, allow_pickle=True)
        results = {int(k): z[k] for k in z.files}
    sampler = prepare_domain_sampler(df)
    t0 = time.time()
    for b in range(start, end + 1):
        if b in results:
            continue
        rng = np.random.default_rng([BOOT_SEED, b])
        try:
            ridx, wb_multiplier = stratified_psu_sample(sampler, rng)
            if ridx.size == 0:
                raise RuntimeError("bootstrap draw contains no analytic GI rows")
            Xb, wb = X[ridx], w[ridx] * wb_multiplier
            tb, eb = time_v[ridx], event_v[ridx]
            beta_b, fit = fit_cox(Xb, wb, tb, eb, beta0=beta_point)
            if (not np.all(np.isfinite(beta_b)) or not np.all(np.isfinite(fit.jac))
                    or np.max(np.abs(fit.jac)) > 1e-3):
                raise RuntimeError("Cox gradient did not converge")
            ht, h0 = breslow_h0(beta_b, Xb, wb, tb, eb)
            m = standardized_mortality(beta_b, Xb, wb, traj_slice, ht, h0, GRID)
            curve = np.vstack([m["no_diabetes"], m["established_pre_cancer_dm"]])
            if not np.all(np.isfinite(curve)):
                raise RuntimeError("non-finite standardized mortality curve")
            results[b] = curve
        except Exception as exc:  # noqa: BLE001 - failed reps are skipped like R tryCatch
            print(f"rep {b}: FAILED ({exc})", flush=True)
        if b % 10 == 0:
            np.savez_compressed(path, **{str(k): v for k, v in results.items()})
            rate = (time.time() - t0) / max(1, b - start + 1)
            print(f"rep {b}/{end}  {rate:.1f}s/rep  converged={len(results)}", flush=True)
    np.savez_compressed(path, **{str(k): v for k, v in results.items()})
    print(f"chunk {start}-{end} done: {len(results)} reps -> {path}")


# ------------------------------------------------------------------ qc + main

def qc_point(beta, names, curves):
    master = load_results()
    principal = {r["term"]: r for r in master["A_principal_5cat"]}
    checks = []
    ok = True
    for term in ["established_pre_cancer_dm", "peri_diagnostic",
                 "post_cancer_dm", "dm_order_unknown"]:
        col = f"trajectory_5cat_rev{term}"
        b = beta[names.index(col)]
        hr = float(np.exp(b))
        ref = principal[term]
        good = abs(hr - ref["HR"]) < 0.005
        ok &= good
        checks.append({"term": term, "hr_python": hr, "hr_locked": ref["HR"],
                       "diff": hr - ref["HR"], "pass": good})
    for g, per_h in LOCKED_MORT.items():
        for horizon, ref in per_h.items():
            m_h = float(curves[g][GRID == horizon][0])
            good = abs(m_h - ref) < 0.002
            ok &= good
            checks.append({"term": f"mortality_{horizon}_{g}", "python": m_h,
                           "locked": ref, "diff": m_h - ref, "pass": good})
    rd = float((curves["established_pre_cancer_dm"] - curves["no_diabetes"])[GRID == 120][0]) * 100
    good = abs(rd - LOCKED_RD_120) < 0.1
    ok &= good
    checks.append({"term": "rd_120_per100", "python": rd,
                   "locked": LOCKED_RD_120, "diff": rd - LOCKED_RD_120,
                   "pass": good})
    return ok, checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-only", action="store_true")
    ap.add_argument("--boot", action="store_true")
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=BOOT_REPS)
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)

    df = load_analysis_frame()
    X, names, traj_slice = build_design(df)
    w = df["sa_wgt_pool"].to_numpy()
    time_v = df["time_months"].to_numpy()
    event_v = df["event"].to_numpy()

    t0 = time.time()
    beta, res = fit_cox(X, w, time_v, event_v)
    print(f"point fit: {time.time() - t0:.1f}s, nit={res.nit}, success={res.success}")

    if args.qc_only or args.boot or args.finalize:
        ht, h0 = breslow_h0(beta, X, w, time_v, event_v)
        curves = standardized_mortality(beta, X, w, traj_slice, ht, h0, GRID)
        ok, checks = qc_point(beta, names, curves)
        report = {"point_qc_pass": ok, "checks": checks,
                  "fit_success": bool(res.success), "nit": int(res.nit)}
        for c in checks:
            print(("PASS" if c["pass"] else "FAIL"), c)
        if not ok:
            print("POINT QC FAILED - stop before drawing", file=sys.stderr)
            sys.exit(2)
        point_path = CHUNK_DIR / "point_curves.npz"
        if not point_path.exists():
            np.savez_compressed(point_path,
                                grid=GRID, no_diabetes=curves["no_diabetes"],
                                established=curves["established_pre_cancer_dm"],
                                beta=beta)
            print("point curves saved")
        else:
            print("point curves already present; reusing point fit")

    if args.boot:
        run_bootstrap(df, X, traj_slice, beta, args.start, args.end)

    if args.finalize:
        reps = []
        for f in sorted(CHUNK_DIR.glob("boot_*.npz")):
            z = np.load(f, allow_pickle=True)
            for k in z.files:
                if not k.isdigit():
                    raise RuntimeError(f"Unexpected non-replicate entry {k!r} in {f}")
                replicate_id = int(k)
                if not 1 <= replicate_id <= BOOT_REPS:
                    raise RuntimeError(f"Bootstrap replicate {replicate_id} is outside 1..{BOOT_REPS}")
                value = np.asarray(z[k], dtype=float)
                if value.shape != (2, len(GRID)) or not np.all(np.isfinite(value)):
                    raise RuntimeError(f"Invalid bootstrap curves for replicate {k} in {f}")
                reps.append((replicate_id, value))
        replicate_ids = [replicate_id for replicate_id, _ in reps]
        expected_ids = set(range(1, BOOT_REPS + 1))
        observed_ids = set(replicate_ids)
        duplicate_ids = sorted({replicate_id for replicate_id in replicate_ids
                                if replicate_ids.count(replicate_id) > 1})
        missing_ids = sorted(expected_ids - observed_ids)
        if duplicate_ids or missing_ids or observed_ids - expected_ids:
            raise RuntimeError(
                "Bootstrap replicate set must contain each ID 1..500 exactly once; "
                f"duplicates={duplicate_ids}, missing={missing_ids[:10]}"
            )
        reps = dict(reps)
        print(f"loaded {len(reps)} bootstrap reps")
        if len(reps) != BOOT_REPS:
            raise RuntimeError(f"Expected {BOOT_REPS} bootstrap reps, found {len(reps)}")
        arr0 = np.array([v[0] for v in reps.values()])   # no_diabetes
        arr1 = np.array([v[1] for v in reps.values()])   # established
        pt = np.load(CHUNK_DIR / "point_curves.npz")
        bands = {
            "grid": GRID.tolist(),
            "no_diabetes": {
                "point": pt["no_diabetes"].tolist(),
                "lo": np.percentile(arr0, 2.5, axis=0).tolist(),
                "hi": np.percentile(arr0, 97.5, axis=0).tolist(),
            },
            "established": {
                "point": pt["established"].tolist(),
                "lo": np.percentile(arr1, 2.5, axis=0).tolist(),
                "hi": np.percentile(arr1, 97.5, axis=0).tolist(),
            },
            "n_boot": len(reps),
            "bootstrap_method": BOOT_METHOD_ID,
            "bootstrap_frame": "full pooled ELIGSTAT=1 Sample Adults; subset to GI after PSU resampling",
            "bootstrap_periods": [f"{lo}-{hi}" for lo, hi, _ in DESIGN_PERIODS],
            "bootstrap_weighting": "Rao-Wu n_h-1; replicate weights = original analysis weights * n_h/(n_h-1) * PSU draw multiplicity",
            "full_population_n": 646201,
            "seed": BOOT_SEED,
        }
        rd = (arr1 - arr0) * 100
        rd120 = rd[:, GRID.tolist().index(120)]
        ci = np.percentile(rd120, [2.5, 97.5])
        bands["rd_120_per100_ci"] = ci.tolist()
        bands["boot_qc_pass"] = bool(np.all(np.isfinite(ci)))
        (OUT_DIR / "sfig3_curves.json").write_text(json.dumps(bands), encoding="utf-8")
        print(f"RD120 CI python=({ci[0]:.2f},{ci[1]:.2f})"
              f"  boot_qc={'PASS' if bands['boot_qc_pass'] else 'FAIL'}  n_boot={len(reps)}")


if __name__ == "__main__":
    main()
