import json
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import *

LMF_OUTPUT_DIR = OUTPUTS_DIR / "lmf"
PIPELINE_VERSION = os.environ.get("NHIS_GI_PIPELINE_VERSION", "v4inc2007")
LMF_OUTPUT_PATH = LMF_OUTPUT_DIR / f"nhis_lmf_parsed_{PIPELINE_VERSION}.parquet"
LMF_AUDIT_PATH = LMF_OUTPUT_DIR / f"lmf_parse_audit_{PIPELINE_VERSION}.json"

INTEGER_COLUMNS = ["eligstat", "mortstat", "ucod_leading", "diabetes", "hyperten", "dodqtr", "dodyear"]
FLOAT_COLUMNS = ["wgt_new", "sa_wgt_new"]
VALID_DODQTR_VALUES = {1, 2, 3, 4}


def normalize_publicid(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip()
    return normalized.replace({"": pd.NA})


def normalize_text_values(frame: pd.DataFrame) -> pd.DataFrame:
    for column in frame.columns:
        if column == "publicid":
            continue
        if frame[column].dtype == "object" or pd.api.types.is_string_dtype(frame[column]):
            normalized = frame[column].astype("string").str.strip()
            frame[column] = normalized.replace({"": pd.NA, ".": pd.NA})
    return frame


def sanitize_dodqtr_series(series: pd.Series) -> pd.Series:
    dodqtr_values = pd.to_numeric(series, errors="coerce")
    return dodqtr_values.where(dodqtr_values.isin(VALID_DODQTR_VALUES))


def parse_year_file(year: int) -> pd.DataFrame:
    file_path = LMF_RAW_ROOT / f"NHIS_{year}_MORT_2019_PUBLIC.dat"
    if not file_path.exists():
        raise FileNotFoundError(f"Missing LMF file: {file_path}")
    frame = pd.read_fwf(file_path, colspecs=LMF_COLSPECS, names=LMF_COLNAMES, dtype={"publicid": str})
    frame["publicid"] = normalize_publicid(frame["publicid"]).astype("string")
    frame = normalize_text_values(frame)
    for column in INTEGER_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("float64")
    frame["dodqtr_raw"] = frame["dodqtr"]
    frame["dodqtr"] = sanitize_dodqtr_series(frame["dodqtr"])
    for column in FLOAT_COLUMNS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce").astype("float64")
    frame["survey_year"] = pd.Series(year, index=frame.index, dtype="int64")
    return frame[LMF_COLNAMES + ["survey_year", "dodqtr_raw"]]


def build_year_audit(frame: pd.DataFrame) -> dict[str, int]:
    eligible_frame = frame[frame["eligstat"].eq(1)]
    dodqtr_raw = pd.to_numeric(frame["dodqtr_raw"], errors="coerce") if "dodqtr_raw" in frame.columns else pd.Series(dtype="float64")
    invalid_dodqtr_count = int((dodqtr_raw.notna() & ~dodqtr_raw.isin(VALID_DODQTR_VALUES)).sum())
    return {
        "raw_rows": int(len(frame)),
        "eligible_rows": int(len(eligible_frame)),
        "deaths": int(eligible_frame["mortstat"].eq(1).sum()),
        "cancer_deaths": int(eligible_frame["ucod_leading"].eq(2).sum()),
        "diabetes_flag": int(eligible_frame["diabetes"].eq(1).sum()),
        "invalid_or_reserved_dodqtr_masked": invalid_dodqtr_count,
    }


def build_audit_payload(year_frames: dict[int, pd.DataFrame], eligible_frame: pd.DataFrame) -> dict[str, object]:
    raw_rows_per_year = {}
    eligible_rows_per_year = {}
    deaths_per_year = {}
    cancer_deaths_per_year = {}
    diabetes_flag_per_year = {}
    invalid_dodqtr_per_year = {}
    publicid_unique_check = True
    for year, frame in year_frames.items():
        year_audit = build_year_audit(frame)
        raw_rows_per_year[str(year)] = year_audit["raw_rows"]
        eligible_rows_per_year[str(year)] = year_audit["eligible_rows"]
        deaths_per_year[str(year)] = year_audit["deaths"]
        cancer_deaths_per_year[str(year)] = year_audit["cancer_deaths"]
        diabetes_flag_per_year[str(year)] = year_audit["diabetes_flag"]
        invalid_dodqtr_per_year[str(year)] = year_audit["invalid_or_reserved_dodqtr_masked"]
        publicid_unique_check = publicid_unique_check and (not frame["publicid"].duplicated().any())
    return {
        "pipeline_version": PIPELINE_VERSION,
        "analysis_years": ANALYSIS_YEARS,
        "raw_rows_per_year": raw_rows_per_year,
        "eligible_rows_per_year": eligible_rows_per_year,
        "deaths_per_year": deaths_per_year,
        "cancer_deaths_per_year": cancer_deaths_per_year,
        "diabetes_flag_per_year": diabetes_flag_per_year,
        "total_eligible_rows": int(len(eligible_frame)),
        "total_deaths": int(eligible_frame["mortstat"].eq(1).sum()),
        "total_cancer_deaths": int(eligible_frame["ucod_leading"].eq(2).sum()),
        "invalid_or_reserved_dodqtr_masked_per_year": invalid_dodqtr_per_year,
        "total_invalid_or_reserved_dodqtr_masked": int(sum(invalid_dodqtr_per_year.values())),
        "publicid_unique_check": publicid_unique_check,
    }


def write_parquet_atomic(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".parquet", delete=False) as handle:
        temporary_path = Path(handle.name)
    try:
        frame.to_parquet(temporary_path, index=False, engine="pyarrow")
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def main() -> None:
    LMF_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    year_frames = {year: parse_year_file(year) for year in ANALYSIS_YEARS}
    combined_frame = pd.concat([year_frames[year] for year in ANALYSIS_YEARS], ignore_index=True)
    eligible_frame = combined_frame[combined_frame["eligstat"].eq(1)].reset_index(drop=True)
    audit_payload = build_audit_payload(year_frames, eligible_frame)
    audit_payload["total_lmf_rows"] = int(len(combined_frame))
    write_parquet_atomic(combined_frame, LMF_OUTPUT_PATH)
    with LMF_AUDIT_PATH.open("w", encoding="utf-8") as handle:
        json.dump(audit_payload, handle, indent=2, ensure_ascii=False)
    print(f"Wrote {LMF_OUTPUT_PATH}")
    print(f"Wrote {LMF_AUDIT_PATH}")
    print(f"Total eligible rows: {audit_payload['total_eligible_rows']}")
    print(f"Total deaths: {audit_payload['total_deaths']}")
    print(f"Total cancer deaths: {audit_payload['total_cancer_deaths']}")


if __name__ == "__main__":
    main()
