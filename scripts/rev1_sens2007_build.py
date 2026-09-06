"""Build the complete 22-year primary inputs from existing 21-year inputs and 2007 public data.
Existing v4 files remain unchanged; outputs retain their v4inc2007 filenames.
"""
import os
import sys
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
ARCHIVE = Path(r"E:/Liao Lab/_cleanup_archive_20260829/data/NHIS 数据库")

# 必须在导入 project_paths 之前设置（模块导入时读取环境变量）
os.environ.setdefault("NHIS_RAW_ROOT", str(ARCHIVE / "Raw_data"))
os.environ.setdefault("LMF_RAW_ROOT", str(ARCHIVE / "Linked_Mortality_Files"))
os.environ["NHIS_GI_PIPELINE_VERSION"] = "v4inc2007"

sys.path.insert(0, str(PROJECT / "scripts"))

import pandas as pd

import project_paths

YEARS22 = list(range(1997, 2019))


def step1_harmonize() -> None:
    import harmonize_nhis_gi_diabetes as H

    frame07 = H.harmonize_year(2007)
    print(f"2007 harmonized rows: {len(frame07)}")
    v4_path = H.HARMONIZED_DIR / "nhis_gi_diabetes_harmonized_v4.parquet"
    frame_v4 = pd.read_parquet(v4_path)
    combined = pd.concat([frame_v4, frame07], ignore_index=True)
    H.write_parquet_atomic(combined, H.HARMONIZED_PATH)
    print(f"Wrote {H.HARMONIZED_PATH} rows={len(combined)}")
    metadata = H.build_metadata(combined)
    import json

    with H.METADATA_PATH.open("w", encoding="utf-8") as fh:
        json.dump(metadata, fh, indent=2, ensure_ascii=False)
    gi07 = int(frame07["gi_any"].fillna(False).sum())
    print(f"2007 GI any: {gi07}")


def step2_lmf() -> None:
    import parse_nhis_lmf as L

    L.ANALYSIS_YEARS = YEARS22
    L.main()


def step3_cohort() -> None:
    import build_analytic_cohort as B

    B.ANALYSIS_YEARS = YEARS22
    B.main()


def step4_domain_fullsample() -> None:
    import build_domain_fullsample as D

    D.ANALYSIS_YEARS = YEARS22
    D.OUT_PATH = D.OUTPUTS_DIR / "cohort" / "analytic_cohort_v4inc2007_domain_fullsample.parquet"
    D.main()


def _map_band(code):
    if code is None:
        return None
    c = str(code).strip()
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


def _map_band4(code):
    if code is None:
        return None
    c = str(code).strip()
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


def step5_income() -> None:
    """既有 21 年 poverty_repaired.csv + 2007（ratcat.csv 内 RAT_CAT3）。"""
    base = pd.read_csv(
        PROJECT / "outputs" / "revision_round1" / "income_repair" / "poverty_repaired.csv",
        dtype="string",
    )
    rat = pd.read_csv(project_paths.NHIS_RAW_ROOT / "2007" / "ratcat.csv", dtype="string")
    rat.columns = [c.strip().lstrip("﻿") for c in rat.columns]
    rows = pd.DataFrame(
        {
            "person_key": pd.NA,
            "year": "2007",
            "hhx": rat["HHX"].str.strip().str.lstrip("0").replace("", "0"),
            "fmx": rat["FMX"].str.strip(),
            "rat_code": rat["RAT_CAT3"].str.strip(),
        }
    )
    rows["poverty_3cat_rep"] = rows["rat_code"].map(_map_band)
    rows["poverty_4cat_rep"] = rows["rat_code"].map(_map_band4)
    # hhx 与既有文件保持同一规范（整数字符串）
    base_hhx = base["hhx"].str.strip()
    base = base.assign(hhx=base_hhx)
    combined = pd.concat([base, rows[base.columns]], ignore_index=True)
    out_dir = PROJECT / "outputs" / "revision_round1_v4_r461_sens2007" / "income_repair"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "poverty_repaired.csv"
    combined.to_csv(out_path, index=False)
    print(f"Wrote {out_path} rows={len(combined)} (2007 rows={len(rows)}, "
          f"2007 mapped={rows['poverty_3cat_rep'].notna().sum()})")


STEPS = {
    "harmonize": step1_harmonize,
    "lmf": step2_lmf,
    "cohort": step3_cohort,
    "domain": step4_domain_fullsample,
    "income": step5_income,
}

if __name__ == "__main__":
    selected = sys.argv[1:] or list(STEPS)
    for name in selected:
        print(f"\n===== step: {name} =====")
        STEPS[name]()
