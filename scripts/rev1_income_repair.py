# rev1_income_repair.py
# 修复 poverty_3cat/poverty_4cat：
#  - 1997-2006: RAT_CAT/FRAT_CAT 类别码（01-14）曾被误当数值比率 → 正确类别映射
#  - 2009,2011-2014: 外层 csv 丢失键列 → 从官方 .dat 定宽文件按码本位置解析 RAT_CAT3/RAT_CAT5
#  - 2015: zip 内完整 csv；2008/2010/2016-2018: 外层 csv（含键）
# Output: outputs/revision_round1_v4_r461_sens2007/income_repair/poverty_repaired.csv plus validation tables.
# The cohort is read-only; the repaired fields are linked by (year, hhx, fmx).

import csv
import argparse
import io
import json
import os
import zipfile
from collections import Counter
from pathlib import Path

import duckdb
import pandas as pd

PROJECT = Path(__file__).resolve().parent.parent


def configured_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    path = Path(value).expanduser() if value else default
    if not path.is_absolute():
        path = PROJECT / path
    return path.resolve()


RAW = configured_path("NHIS_RAW_ROOT", PROJECT / "data" / "NHIS_raw")
DEFAULT_OUT = configured_path(
    "REV1_INCOME_REPAIR_DIR",
    PROJECT / "outputs" / "revision_round1_v4_r461_sens2007" / "income_repair",
)
COHORT = configured_path(
    "REV1_INCOME_COHORT_PATH",
    PROJECT / "outputs" / "cohort" / "analytic_cohort_v4inc2007.parquet",
)

# year -> (source, variable)
# source: 'csv' = 外层 csv; 'zipcsv' = zip 内 csv; 'dat' = zip 内定宽 dat（位置 1-based 闭区间）
YEAR_SPEC = {
    1997: ("csv", "FRAT_CAT"),
    1998: ("csv", "FRAT_CAT"), 1999: ("csv", "FRAT_CAT"), 2000: ("csv", "FRAT_CAT"),
    2001: ("csv", "FRAT_CAT"), 2002: ("csv", "FRAT_CAT"), 2003: ("csv", "FRAT_CAT"),
    2004: ("csv", "RAT_CAT"), 2005: ("csv", "RAT_CAT"), 2006: ("csv", "RAT_CAT"),
    2007: ("csv", "RAT_CAT3"),
    2008: ("csv", "RAT_CAT3"),
    2009: ("dat", "RAT_CAT3", (158, 159)),
    2010: ("csv", "RAT_CAT3"),
    2011: ("dat", "RAT_CAT3", (204, 205)),
    2012: ("dat", "RAT_CAT3", (204, 205)),
    2013: ("dat", "RAT_CAT3", (206, 207)),
    2014: ("dat", "RAT_CAT5", (196, 197)),
    2015: ("zipcsv", "RAT_CAT5"),
    2016: ("csv", "RAT_CAT5"), 2017: ("csv", "RAT_CAT5"), 2018: ("csv", "RAT_CAT5"),
}
DAT_KEY_POS = {"srvy_yr": (3, 6), "hhx": (7, 12), "fmx": (13, 14)}  # 已用原始字节验证 (2011)

# 类别映射（1997-2006 仅 01-14；2007+ RAT_CAT3/5 含 15-18 追问粗码）
def map_band(code):
    if code is None:
        return None
    c = code.strip()
    if not c or not c.isdigit():
        return None
    v = int(c)
    if v in (96, 97, 98, 99) or v == 0:
        return None
    if 1 <= v <= 7 or v in (15, 16):
        return "lt_2_0"
    if 8 <= v <= 11 or v == 17:
        return "2_0_to_3_99"
    if 12 <= v <= 14 or v == 18:
        return "ge_4_0"
    return None

def map_band4(code):
    if code is None:
        return None
    c = code.strip()
    if not c or not c.isdigit():
        return None
    v = int(c)
    if v in (96, 97, 98, 99) or v == 0:
        return None
    if 1 <= v <= 3 or v == 15:
        return "lt_1_0"
    if 4 <= v <= 7 or v == 16:
        return "1_0_to_1_99"
    if 8 <= v <= 11 or v == 17:
        return "2_0_to_3_99"
    if 12 <= v <= 14 or v == 18:
        return "ge_4_0"
    return None


def norm_key(s):
    if s is None:
        return None
    t = str(s).strip()
    if t == "" or t.upper() == "NA":
        return None
    try:
        return str(int(float(t)))
    except ValueError:
        return t


def read_csv_year(year, var, fh=None):
    """外层 csv 或 zip 内 csv → {(hhx,fmx): code}"""
    if fh is None:
        fh = open(RAW / str(year) / ("ratcat.csv" if year == 2007 else "familyxx.csv"), encoding="utf-8-sig", newline="")
        close = True
    else:
        close = False
    out = {}
    with fh if close else fh:
        r = csv.DictReader(fh)
        fieldnames = {f.strip().lstrip("﻿") for f in r.fieldnames}
        assert var in fieldnames, f"{year}: {var} 不在列中"
        hhx_col = "HHX" if "HHX" in fieldnames else ("HH" if "HH" in fieldnames else None)
        assert hhx_col and "FMX" in fieldnames, f"{year}: 缺键列"
        for row in r:
            h, f = norm_key(row.get(hhx_col)), norm_key(row.get("FMX"))
            if h is None:
                continue
            out[(h, f)] = (row.get(var) or "").strip()
    return out


def read_dat_year(year, var, pos):
    """zip 内定宽 dat → {(hhx,fmx): code}；附带健全性检查"""
    zp = zipfile.ZipFile(RAW / str(year) / "familyxx.zip")
    name = [n for n in zp.namelist() if n.lower().endswith(".dat")][0]
    out = {}
    rectype_bad = year_bad = 0
    with zp.open(name) as fh:
        for raw in io.TextIOWrapper(fh, encoding="latin-1"):
            line = raw.rstrip("\r\n")
            if len(line) < pos[1]:
                continue
            if line[0:2].strip() != "60":
                rectype_bad += 1
                continue
            yr = line[DAT_KEY_POS["srvy_yr"][0]-1:DAT_KEY_POS["srvy_yr"][1]].strip()
            if yr != str(year):
                year_bad += 1
            h = norm_key(line[DAT_KEY_POS["hhx"][0]-1:DAT_KEY_POS["hhx"][1]])
            f = norm_key(line[DAT_KEY_POS["fmx"][0]-1:DAT_KEY_POS["fmx"][1]])
            out[(h, f)] = line[pos[0]-1:pos[1]].strip()
    if rectype_bad or year_bad:
        print(f"  !! {year}: rectype_bad={rectype_bad} year_bad={year_bad}")
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        help="Explicit output directory. Required unless REV1_INCOME_REPAIR_DIR is set.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow writing into an existing non-protected output directory.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_value = args.output_dir or os.environ.get("REV1_INCOME_REPAIR_DIR")
    if not output_value:
        raise SystemExit("Refusing implicit output: pass --output-dir or set REV1_INCOME_REPAIR_DIR.")
    out = Path(output_value).expanduser()
    if not out.is_absolute():
        out = PROJECT / out
    out = out.resolve()
    protected = [
        (PROJECT / "revision" / "submission_ready_v4.6.1_20260830").resolve(),
        (PROJECT / "revision" / "submission_ready_v4.6.2_20260902").resolve(),
        (PROJECT / "revision" / "submission_ready_v4.6.3_20260902").resolve(),
    ]
    if any(out == path or path in out.parents for path in protected):
        raise SystemExit(f"Refusing to write inside protected formal directory: {out}")
    if out.exists() and any(out.iterdir()) and not args.force:
        raise SystemExit(f"Refusing to overwrite non-empty output directory without --force: {out}")
    out.mkdir(parents=True, exist_ok=True)

    frames = []
    dist_rows = []
    for year, spec in YEAR_SPEC.items():
        kind = spec[0]
        var = spec[1]
        if kind == "csv":
            d = read_csv_year(year, var)
        elif kind == "zipcsv":
            zp = zipfile.ZipFile(RAW / str(year) / "familyxx.zip")
            name = [n for n in zp.namelist() if n.lower().endswith(".csv")][0]
            d = read_csv_year(year, var, fh=io.TextIOWrapper(zp.open(name), encoding="utf-8-sig", newline=""))
        else:
            d = read_dat_year(year, var, spec[2])
        dist = Counter(d.values())
        dist_rows.append({"year": year, "var": var, "families": len(d),
                          "dist": dict(sorted(dist.items()))})
        df = pd.DataFrame({"hhx": [k[0] for k in d], "fmx": [k[1] for k in d],
                           "rat_code": list(d.values())})
        df["year"] = year
        frames.append(df)
        print(f"{year}: {var} families={len(d)} codes={len(dist)}")

    rat = pd.concat(frames, ignore_index=True)
    rat["poverty_3cat_rep"] = rat["rat_code"].map(map_band)
    rat["poverty_4cat_rep"] = rat["rat_code"].map(map_band4)

    # 冻结队列键
    co = duckdb.sql(f"select year, hhx, fmx, person_key, poverty_3cat as old from '{COHORT}'").fetchdf()
    co["hhx"] = co["hhx"].map(norm_key)
    co["fmx"] = co["fmx"].map(norm_key)
    co["year"] = co["year"].astype(int)

    m = co.merge(rat, on=["year", "hhx", "fmx"], how="left", validate="many_to_one")
    merge_ok = m["rat_code"].notna().mean()
    print(f"\n队列行数={len(m)}  合并命中率={merge_ok:.4f}")

    # 逐年命中率与新旧分布
    by_year = (m.assign(hit=m["rat_code"].notna())
                .groupby("year")
                .agg(n=("hit", "size"), merge_rate=("hit", "mean"),
                     old_missing=("old", lambda s: (s.isna() | (s == "missing")).mean()),
                     new_missing=("poverty_3cat_rep", lambda s: s.isna().mean()),
                     new_lt2=("poverty_3cat_rep", lambda s: (s == "lt_2_0").mean()))
                .round(3).reset_index())
    print(by_year.to_string(index=False))

    m[["person_key", "year", "hhx", "fmx", "rat_code", "poverty_3cat_rep", "poverty_4cat_rep"]].to_csv(
        out / "poverty_repaired.csv", index=False)
    by_year.to_csv(out / "validation_by_year.csv", index=False)
    with open(out / "raw_code_distributions.json", "w") as fh:
        json.dump(dist_rows, fh, indent=1)

    # 全体新旧对比
    overall = {
        "old_lt2_share": float(((m["old"] == "lt_2_0")).mean()),
        "old_missing_share": float((m["old"].isna() | (m["old"] == "missing")).mean()),
        "new_lt2_share": float((m["poverty_3cat_rep"] == "lt_2_0").mean()),
        "new_missing_share": float(m["poverty_3cat_rep"].isna().mean()),
        "merge_rate": float(merge_ok),
    }
    with open(out / "validation_overall.json", "w") as fh:
        json.dump(overall, fh, indent=1)
    print(json.dumps(overall, indent=1))


if __name__ == "__main__":
    main()
