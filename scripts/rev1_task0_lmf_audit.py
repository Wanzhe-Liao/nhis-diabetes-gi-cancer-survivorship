"""Task 0 audit for the three NHIS-LMF flow denominators.

This script reads the raw Sample Adult and public-use 2019 LMF files only.  It
does not rebuild the analytic cohort and does not modify any locked result.
The three reported layers are kept distinct:

1. Sample Adult records whose 14-character publicid is present in the LMF;
2. those matched records with ELIGSTAT == 1;
3. those eligible records with MORTSTAT == 1.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]


def configured_path(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    path = Path(value).expanduser() if value else default
    if not path.is_absolute():
        path = PROJECT / path
    return path.resolve()


RAW_ROOT = configured_path("NHIS_RAW_ROOT", PROJECT / "data" / "NHIS_raw")
LMF_ROOT = configured_path("LMF_RAW_ROOT", PROJECT / "data" / "LMF_raw")
LOCKED_OUT = (PROJECT / "outputs" / "revision_round1").resolve()
OUT = configured_path(
    "REV1_OUTPUT_DIR",
    PROJECT / "outputs" / "revision_round1_v4_r461_sens2007",
)
if OUT == LOCKED_OUT:
    raise RuntimeError(
        "The historical outputs/revision_round1 directory is read-only; "
        "choose the versioned R 4.6.1 output directory."
    )

YEARS = list(range(1997, 2019))

# The harmonization code uses these year-specific identifiers to construct the
# 14-character publicid = survey year + HHX + FMX + person key.
ID_COLUMNS = {
    1997: ("SRVY_YR", "HHX", "FMX", "PX"),
    1998: ("SRVY_YR", "HHX", "FMX", "PX"),
    1999: ("SRVY_YR", "HHX", "FMX", "PX"),
    2000: ("SRVY_YR", "HHX", "FMX", "PX"),
    2001: ("SRVY_YR", "HHX", "FMX", "PX"),
    2002: ("SRVY_YR", "HHX", "FMX", "PX"),
    2003: ("SRVY_YR", "HHX", "FMX", "PX"),
    2004: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2005: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2006: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2007: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2008: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2009: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2010: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2011: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2012: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2013: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2014: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2015: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2016: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2017: ("SRVY_YR", "HHX", "FMX", "FPX"),
    2018: ("SRVY_YR", "HHX", "FMX", "FPX"),
}


def fmt(value: object, width: int) -> str:
    """Format an identifier component as the harmonization code does."""
    text = str(value).strip()
    if text in {"", "nan", "None", "<NA>"}:
        return ""
    try:
        text = str(int(float(text)))
    except ValueError:
        pass
    return text.zfill(width)


def sample_adult_ids(year: int) -> list[str]:
    year_col, hhx_col, fmx_col, person_col = ID_COLUMNS[year]
    path = RAW_ROOT / str(year) / "samadult.csv"
    frame = pd.read_csv(
        path,
        usecols=[year_col, hhx_col, fmx_col, person_col],
        dtype=str,
        low_memory=False,
    )
    ids = []
    for _, row in frame.iterrows():
        survey_year = fmt(row[year_col], 4)
        # The 2005 public-use file stores SRVY_YR as the short code "5";
        # harmonization replaces such short codes with the directory year.
        try:
            if int(float(str(row[year_col]).strip())) < 1000:
                survey_year = str(year)
        except ValueError:
            pass
        ids.append(
            survey_year
            + fmt(row[hhx_col], 6)
            + fmt(row[fmx_col], 2)
            + fmt(row[person_col], 2)
        )
    if any(not value for value in ids):
        raise ValueError(f"{year}: missing publicid component")
    if len(ids) != len(set(ids)):
        raise ValueError(f"{year}: duplicate Sample Adult publicid")
    return ids


def lmf_status_by_id(year: int, selected_ids: set[str]) -> dict[str, tuple[str, str]]:
    """Return status fields for selected publicids from the fixed-width LMF."""
    path = LMF_ROOT / f"NHIS_{year}_MORT_2019_PUBLIC.dat"
    found: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="latin-1", newline="") as handle:
        for line in handle:
            if len(line) < 16:
                continue
            publicid = line[:14]
            if publicid in selected_ids:
                found[publicid] = (line[14:15], line[15:16])
    return found


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    total_sample = 0
    total_matched = 0
    total_eligible = 0
    total_deaths = 0

    for year in YEARS:
        ids = sample_adult_ids(year)
        status = lmf_status_by_id(year, set(ids))
        matched = sum(publicid in status for publicid in ids)
        eligible = sum(status.get(publicid, ("", ""))[0] == "1" for publicid in ids)
        deaths = sum(
            status.get(publicid, ("", ""))[0] == "1"
            and status.get(publicid, ("", ""))[1] == "1"
            for publicid in ids
        )
        rows.append(
            {
                "survey_year": year,
                "sample_adult_n": len(ids),
                "matched_to_lmf_n": matched,
                "eligstat_1_n": eligible,
                "mortstat_1_among_eligible_n": deaths,
                "matched_rate": matched / len(ids),
                "eligibility_rate_among_matched": eligible / matched if matched else None,
                "eligibility_rate_of_sample_adult": eligible / len(ids),
                "death_rate_among_eligible": deaths / eligible if eligible else None,
            }
        )
        total_sample += len(ids)
        total_matched += matched
        total_eligible += eligible
        total_deaths += deaths

    summary = {
        "survey_years": YEARS,
        "sample_adult_n": total_sample,
        "matched_to_lmf_n": total_matched,
        "eligstat_1_n": total_eligible,
        "mortstat_1_among_eligible_n": total_deaths,
        "matched_rate": total_matched / total_sample,
        "eligibility_rate_among_matched": total_eligible / total_matched,
        "eligibility_rate_of_sample_adult": total_eligible / total_sample,
        "death_rate_among_eligible": total_deaths / total_eligible,
        "interpretation": {
            "merge_layer": "Sample Adult publicid present in the 2019 LMF",
            "eligibility_layer": "Matched record with ELIGSTAT == 1",
            "death_layer": "ELIGSTAT == 1 record with MORTSTAT == 1",
        },
    }
    pd.DataFrame(rows).to_csv(OUT / "lmf_linkage_audit_rev1.csv", index=False)
    with (OUT / "lmf_linkage_audit_rev1.json").open("w", encoding="utf-8") as handle:
        json.dump({"summary": summary, "by_year": rows}, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
