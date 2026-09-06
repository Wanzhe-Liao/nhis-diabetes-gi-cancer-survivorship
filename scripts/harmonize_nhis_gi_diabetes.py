import io
import hashlib
import json
import os
import sys
import tempfile
import zipfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import *

HARMONIZED_DIR = OUTPUTS_DIR / "harmonized"
HARMONIZATION_VERSION = os.environ.get("NHIS_GI_PIPELINE_VERSION", "v4inc2007")
HARMONIZED_PATH = HARMONIZED_DIR / f"nhis_gi_diabetes_harmonized_{HARMONIZATION_VERSION}.parquet"
METADATA_PATH = HARMONIZED_DIR / f"harmonization_metadata_{HARMONIZATION_VERSION}.json"

SITE_FLAG_COLUMNS = {
    "colon": "colon_flag",
    "esophagus": "esoph_flag",
    "gallbladder": "gallbladder_flag",
    "liver": "liver_flag",
    "pancreas": "pancreas_flag",
    "rectum": "rectum_flag",
    "stomach": "stomach_flag",
}

RENAMED_CANONICAL_COLUMNS = [
    "person_key",
    "hhx",
    "fmx",
    "srvy_yr",
    "intv_qrt",
    "age",
    "sex",
    "race",
    "hispanic",
    "education",
    "region",
    "poverty",
    "marital",
    "canev",
    "cnkind7", "canage7",
    "cnkind8", "canage8",
    "cnkind9", "canage9",
    "cnkind13", "canage13",
    "cnkind19", "canage19",
    "cnkind21", "canage21",
    "cnkind25", "canage25",
    "dm_ever",
    "dm_age",
    "dm_duration",
    "dm_med_pills",
    "dm_med_insulin",
    "smoking_ever",
    "smoking_status",
    "bmi",
    "afford_rx",
    "wtfa_sa",
    "design_strata",
    "design_psu",
    "hypertension_ever",
    "cholesterol_high_ever",
    "chd_ever",
    "stroke_ever",
    "alcohol_status",
    "alcohol_freq_past_year",
    "phys_mod_freq",
    "phys_vig_freq",
    "insurance_notcov",
    "insurance_cover",
    "private_flag",
    "other_public_flag",
    "military_flag",
    "insurance_type",
    "srh",
    "insurance_cover65",
    "medicare_flag",
    "medicaid_flag",
]

DERIVED_COLUMNS = [
    "year",
    "publicid",
    "colon_flag",
    "esoph_flag",
    "gallbladder_flag",
    "liver_flag",
    "pancreas_flag",
    "rectum_flag",
    "stomach_flag",
    "gi_any",
    "smoking_3cat",
    "dm_med_intensity",
    "source_dm_var",
    "source_dm_age_var",
    "source_design_strata_var",
    "source_design_psu_var",
    "source_person_key_var",
    "race_structural_missing",
    "race_ethnicity_4cat",
    "phys_active_any",
    "cost_barrier_rx",
    "alcohol_status_raw",
    "education_5cat",
    "education_4cat",
    "poverty_4cat",
    "poverty_3cat",
    "source_poverty_var",
]

FINAL_COLUMNS = ["year", "publicid", "srvy_yr", "hhx", "fmx", "person_key"] + [column for column in RENAMED_CANONICAL_COLUMNS if column not in {"person_key", "hhx", "fmx", "srvy_yr"}] + [column for column in DERIVED_COLUMNS if column != "year" and column != "publicid"]

# Round 3: personsx variables available for all years
PERSONSX_SOURCE_COLUMNS = {
    "srh": "PHSTAT",
    "education": "EDUC",
    "insurance_notcov": "NOTCOV",
    "age": "AGE_P",
    "sex": "SEX",
    "design_strata": "STRAT_P",
    "design_psu": "PSU",
}

# Some NHIS years keep demographics/design variables outside samadult. Use
# source-file fallbacks rather than treating whole years as design missing.
PERSONSX_FALLBACK_CANDIDATES = {
    "age": ["AGE_P"],
    "sex": ["SEX"],
    "race": ["MRACRPI2", "RACERPI2", "RACEIMP2", "RACRECI2"],
    "hispanic": ["HISPAN_I"],
    "education": ["EDUC1", "EDUC"],
    "region": ["REGION"],
    "marital": ["R_MARITL"],
    "design_strata": ["STRAT_P", "STRATUM"],
    "design_psu": ["PSU_P", "PSU"],
}

HOUSEHLD_FALLBACK_CANDIDATES = {
    "intv_qrt": ["INTV_QRT"],
    "region": ["REGION"],
    "design_strata": ["STRAT_P", "STRATUM"],
    "design_psu": ["PSU_P", "PSU"],
}

INSURANCE_RECODE_CANDIDATES = {
    "medicare_flag": ["MEDICARE"],
    "medicaid_flag": ["MEDICAID"],
    "private_flag": ["PRIVATE"],
    "other_public_flag": ["OTHPUB", "OTHERPUB"],
    "military_flag": ["MILITARN", "MILCARE", "MILITARY"],
}

# Official NHIS Sample Adult files for 2011 and 2015 include race/Hispanic
# fields. Earlier local extracts omitted them; the repaired pipeline no longer
# treats any analysis year as structurally missing for race.
STRUCTURAL_UNKNOWN_RACE_YEARS: set[int] = set()


def normalize_key_series(series: pd.Series) -> pd.Series:
    """Canonicalize NHIS merge keys without losing nonnumeric identifiers."""
    normalized = series.astype("string").str.strip().replace({"": pd.NA})
    numeric = pd.to_numeric(normalized, errors="coerce")
    integer_like = numeric.notna() & numeric.eq(numeric.round())
    normalized.loc[integer_like] = numeric.loc[integer_like].astype("Int64").astype("string")
    return normalized


def format_publicid_component(series: pd.Series, width: int) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    formatted = normalize_key_series(series)
    numeric_mask = numeric.notna()
    formatted.loc[numeric_mask] = numeric.loc[numeric_mask].astype("Int64").astype("string").str.zfill(width)
    return formatted.replace({"": pd.NA})


def build_rename_map(year_variable_map: dict[str, str | None]) -> dict[str, str]:
    rename_map = {}
    for canonical_name, source_name in year_variable_map.items():
        if source_name is None:
            continue
        if canonical_name == "survey_year_col":
            rename_map[source_name] = "srvy_yr"
        elif canonical_name == "weight_sa":
            rename_map[source_name] = "wtfa_sa"
        elif canonical_name == "dm_pills":
            rename_map[source_name] = "dm_med_pills"
        elif canonical_name == "dm_insulin":
            rename_map[source_name] = "dm_med_insulin"
        else:
            rename_map[source_name] = canonical_name
    return rename_map


def binary_with_missing(series: pd.Series, true_values: set[int], false_values: set[int]) -> pd.Series:
    numeric_series = pd.to_numeric(series, errors="coerce")
    output = pd.Series(pd.NA, index=series.index, dtype="boolean")
    output.loc[numeric_series.isin(true_values)] = True
    output.loc[numeric_series.isin(false_values)] = False
    return output


def numeric_with_special_missing(series: pd.Series) -> pd.Series:
    numeric_series = pd.to_numeric(series, errors="coerce")
    return numeric_series.mask(numeric_series.isin([96, 97, 98, 99]))


def clean_activity_frequency(series: pd.Series) -> pd.Series:
    """Return interpretable weekly frequency for MODFREQW/VIGFREQW.

    NHIS codes 00-28 as weekly frequency, 95 as never, 96 as unable, and
    97-99 as nonresponse. Never/unable contribute zero sessions; nonresponse
    and out-of-range values remain missing.
    """
    raw = pd.to_numeric(series, errors="coerce")
    cleaned = raw.where(raw.between(0, 28))
    cleaned.loc[raw.isin([95, 96])] = 0
    return cleaned


def insurance_plan_recode_boolean(series: pd.Series) -> pd.Series:
    """Decode MEDICARE/MEDICAID/PRIVATE recodes (1/2=yes, 3=no)."""
    return binary_with_missing(series, {1, 2}, {3})


def other_public_recode_boolean(series: pd.Series) -> pd.Series:
    """Decode OTHPUB/OTHERPUB recodes (1=yes, 2=no; other codes unknown)."""
    return binary_with_missing(series, {1}, {2})


def military_recode_boolean(series: pd.Series, year: int) -> pd.Series:
    """Decode the year-specific NHIS military coverage recode."""
    if year <= 2003:
        return binary_with_missing(series, {1, 2, 3, 4}, {5})
    return binary_with_missing(series, {1, 2}, {3})


def derive_alcohol_status(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    """Harmonize the two materially different ALCSTAT coding schemes."""
    raw = clean_numeric(frame["alcohol_status"])
    frame["alcohol_status_raw"] = raw
    frame["alcohol_status"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    frame.loc[raw.eq(1), "alcohol_status"] = "never"
    if year <= 2003:
        # ALCSTAT1: 1 lifetime abstainer, 2 former, 3 current, 9 unknown.
        frame.loc[raw.eq(2), "alcohol_status"] = "former"
        frame.loc[raw.eq(3), "alcohol_status"] = "current"
    else:
        # ALCSTAT: 2-4 former; 5-9 current (frequency/intensity subtypes).
        frame.loc[raw.isin([2, 3, 4]), "alcohol_status"] = "former"
        frame.loc[raw.isin([5, 6, 7, 8, 9]), "alcohol_status"] = "current"
    return frame


def clean_bmi(series: pd.Series) -> pd.Series:
    numeric_series = pd.to_numeric(series, errors="coerce")
    numeric_series = numeric_series.mask(numeric_series >= 9997)
    scaled_series = numeric_series.where(numeric_series.isna() | (numeric_series < 100), numeric_series / 100)
    return scaled_series.mask((scaled_series >= 99) & (scaled_series < 100))


def clean_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_publicid(frame: pd.DataFrame) -> pd.Series:
    publicid_parts = pd.DataFrame(index=frame.index)
    publicid_parts["srvy_yr"] = format_publicid_component(frame["srvy_yr"], 4)
    publicid_parts["hhx"] = format_publicid_component(frame["hhx"], 6)
    publicid_parts["fmx"] = format_publicid_component(frame["fmx"], 2)
    publicid_parts["person_key"] = format_publicid_component(frame["person_key"], 2)
    srvy_numeric = pd.to_numeric(publicid_parts["srvy_yr"], errors="coerce")
    year_numeric = pd.to_numeric(frame.get("year"), errors="coerce")
    short_year_mask = srvy_numeric.notna() & srvy_numeric.lt(1000) & year_numeric.notna()
    publicid_parts.loc[short_year_mask, "srvy_yr"] = (
        year_numeric.loc[short_year_mask].astype("Int64").astype("string")
    )
    complete_mask = publicid_parts.notna().all(axis=1)
    publicid = publicid_parts["srvy_yr"] + publicid_parts["hhx"] + publicid_parts["fmx"] + publicid_parts["person_key"]
    return publicid.where(complete_mask, pd.NA)


def add_missing_columns(frame: pd.DataFrame) -> pd.DataFrame:
    for column in RENAMED_CANONICAL_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA
    return frame


def derive_education_5cat(frame: pd.DataFrame) -> pd.DataFrame:
    """Harmonize EDUC/EDUC1 credentials using the official 00-21 codes."""
    frame["education_5cat"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    educ = pd.to_numeric(frame["education"], errors="coerce")
    frame.loc[educ.between(0, 12), "education_5cat"] = "lt_high_school"
    frame.loc[educ.between(13, 14), "education_5cat"] = "high_school"
    frame.loc[educ.between(15, 17), "education_5cat"] = "some_college"
    frame.loc[educ.between(18, 21), "education_5cat"] = "college_grad"
    return frame


def derive_education_4cat(frame: pd.DataFrame) -> pd.DataFrame:
    """Expose the four credential groups used by the analysis."""
    frame["education_4cat"] = frame["education_5cat"].astype("string")
    return frame


def derive_poverty_4cat(frame: pd.DataFrame) -> pd.DataFrame:
    frame["poverty_4cat"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    # FRAT_CAT/RAT_CAT/RAT_CAT3/RAT_CAT5 are category codes, not raw ratios.
    pov = pd.to_numeric(frame["poverty"], errors="coerce")
    frame.loc[pov.isin([1, 2, 3, 15]), "poverty_4cat"] = "lt_1_0"
    frame.loc[pov.isin([4, 5, 6, 7, 16]), "poverty_4cat"] = "1_0_to_1_99"
    frame.loc[pov.isin([8, 9, 10, 11, 17]), "poverty_4cat"] = "2_0_to_3_99"
    frame.loc[pov.isin([12, 13, 14, 18]), "poverty_4cat"] = "ge_4_0"
    return frame


def derive_poverty_3cat(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse the official four poverty bands at a ratio of 2.0."""
    frame["poverty_3cat"] = frame["poverty_4cat"].astype("string")
    frame.loc[frame["poverty_3cat"].isin(["lt_1_0", "1_0_to_1_99"]), "poverty_3cat"] = "lt_2_0"
    return frame


def derive_race_ethnicity_4cat(frame: pd.DataFrame) -> pd.DataFrame:
    """Create mutually exclusive Hispanic/race groups across survey years."""
    race = pd.to_numeric(frame["race"], errors="coerce")
    hispanic = pd.to_numeric(frame["hispanic"], errors="coerce")
    is_hispanic = hispanic.between(0, 11)
    is_non_hispanic = hispanic.eq(12)
    out = pd.Series(pd.NA, index=frame.index, dtype="string")
    out.loc[is_hispanic] = "hispanic"
    out.loc[is_non_hispanic & race.eq(1)] = "non_hispanic_white"
    out.loc[is_non_hispanic & race.eq(2)] = "non_hispanic_black"
    out.loc[is_non_hispanic & race.notna() & ~race.isin([1, 2, 97, 98, 99])] = "non_hispanic_other"
    frame["race_ethnicity_4cat"] = out
    return frame


def read_year_frame(year: int) -> pd.DataFrame:
    year_variable_map = YEAR_VARIABLE_MAP[year]
    samadult_ext = YEAR_SAMADULT_EXTENSION_MAP.get(year, {})
    samadult_path = NHIS_RAW_ROOT / str(year) / "samadult.csv"
    # Collect all samadult columns
    base_cols = {c for c in year_variable_map.values() if c is not None}
    ext_cols = {c for c in samadult_ext.values() if c is not None}
    requested_columns = sorted(base_cols | ext_cols)
    samadult_header = set(pd.read_csv(samadult_path, nrows=0).columns)
    fallback_columns = {
        source
        for candidates in list(PERSONSX_FALLBACK_CANDIDATES.values()) + list(HOUSEHLD_FALLBACK_CANDIDATES.values())
        for source in candidates
    }
    cancerxx_path = NHIS_RAW_ROOT / str(year) / "cancerxx.csv"
    if cancerxx_path.exists():
        cancerxx_header = set(pd.read_csv(cancerxx_path, nrows=0).columns)
        if "WTFA_SA" in cancerxx_header:
            fallback_columns.add("WTFA_SA")
    fwf_columns = set(YEAR_SAMADULT_FWF_SUPPLEMENT_MAP.get(year, {}).get("fields", {}))
    zip_csv_columns = set(YEAR_SAMADULT_ZIP_CSV_SUPPLEMENT_MAP.get(year, {}).get("fields", []))
    missing_requested = set(requested_columns) - samadult_header
    unexpected_missing = missing_requested - fallback_columns - fwf_columns - zip_csv_columns
    if unexpected_missing:
        missing = ", ".join(sorted(unexpected_missing))
        raise ValueError(f"Required samadult columns missing for {year}: {missing}")
    use_columns = [column for column in requested_columns if column in samadult_header]
    dtype_map = {
        year_variable_map["survey_year_col"]: "string",
        year_variable_map["hhx"]: "string",
        year_variable_map["fmx"]: "string",
        year_variable_map["person_key"]: "string",
    }
    frame = pd.read_csv(samadult_path, usecols=use_columns, dtype=dtype_map, low_memory=False)
    # Rename using combined map
    combined_map = {**year_variable_map, **samadult_ext}
    frame = frame.rename(columns=build_rename_map(combined_map))
    frame = add_missing_columns(frame)
    frame["hhx"] = normalize_key_series(frame["hhx"])
    frame["fmx"] = normalize_key_series(frame["fmx"])
    frame["person_key"] = normalize_key_series(frame["person_key"])
    frame["srvy_yr"] = normalize_key_series(frame["srvy_yr"])
    frame["year"] = year
    frame["source_dm_var"] = year_variable_map["dm_ever"]
    frame["source_dm_age_var"] = year_variable_map["dm_age"]
    frame["source_design_strata_var"] = year_variable_map["design_strata"]
    frame["source_design_psu_var"] = year_variable_map["design_psu"]
    frame["source_person_key_var"] = year_variable_map["person_key"]
    return frame


def read_samadult_fwf_supplement(year: int) -> pd.DataFrame | None:
    """Read fields omitted by a local CSV extract from the official archive."""
    csv_spec = YEAR_SAMADULT_ZIP_CSV_SUPPLEMENT_MAP.get(year)
    spec = YEAR_SAMADULT_FWF_SUPPLEMENT_MAP.get(year)
    if spec is None and csv_spec is None:
        return None
    zip_path = NHIS_RAW_ROOT / str(year) / "samadult.zip"
    if not zip_path.exists():
        raise FileNotFoundError(f"Required official Sample Adult archive missing for {year}: {zip_path}")
    with zipfile.ZipFile(zip_path) as archive:
        if csv_spec is not None:
            members = {name.lower(): name for name in archive.namelist()}
            member = members.get(str(csv_spec["member"]).lower())
            if member is None:
                raise ValueError(f"Official Sample Adult CSV member missing for {year}")
            usecols = list(csv_spec["keys"]) + list(csv_spec["fields"])
            with archive.open(member) as raw_handle:
                supplement = pd.read_csv(raw_handle, usecols=usecols, dtype="string", low_memory=False)
            source_to_canonical = build_rename_map(YEAR_VARIABLE_MAP[year])
            supplement = supplement.rename(columns=source_to_canonical)
            for key in ["srvy_yr", "hhx", "fmx", "person_key"]:
                supplement[key] = normalize_key_series(supplement[key])
            keys = ["hhx", "fmx", "person_key"]
            if supplement.duplicated(subset=keys).any():
                raise ValueError(f"Duplicate Sample Adult archive merge keys detected for {year}")
            return supplement

        members = {name.lower(): name for name in archive.namelist()}
        requested_member = str(spec["member"]).lower()
        if requested_member not in members:
            dat_members = [name for name in archive.namelist() if name.lower().endswith(".dat")]
            if len(dat_members) != 1:
                raise ValueError(f"Could not identify Sample Adult DAT member for {year}")
            member = dat_members[0]
        else:
            member = members[requested_member]
        names = list(spec["keys"]) + list(spec["fields"])
        positions = list(spec["keys"].values()) + list(spec["fields"].values())
        colspecs = [(start - 1, end) for start, end in positions]
        with archive.open(member) as raw_handle:
            text_handle = io.TextIOWrapper(raw_handle, encoding="latin-1")
            supplement = pd.read_fwf(text_handle, colspecs=colspecs, names=names, dtype="string")

    source_to_canonical = build_rename_map(YEAR_VARIABLE_MAP[year])
    supplement = supplement.rename(columns=source_to_canonical)
    for key in ["srvy_yr", "hhx", "fmx", "person_key"]:
        supplement[key] = normalize_key_series(supplement[key])
    keys = ["hhx", "fmx", "person_key"]
    if supplement.duplicated(subset=keys).any():
        raise ValueError(f"Duplicate Sample Adult DAT merge keys detected for {year}")
    return supplement


def merge_samadult_fwf_supplement(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    supplement = read_samadult_fwf_supplement(year)
    if supplement is None:
        return frame
    original_rows = len(frame)
    keys = ["hhx", "fmx", "person_key"]
    value_columns = [column for column in supplement.columns if column not in keys + ["srvy_yr"]]
    supplement = supplement.rename(columns={name: f"{name}_official" for name in value_columns})
    frame = frame.merge(supplement[keys + [f"{name}_official" for name in value_columns]],
                        on=keys, how="left", validate="one_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after Sample Adult DAT merge for {year}")
    for column in value_columns:
        official = f"{column}_official"
        frame[column] = frame[column].where(frame[column].notna(), frame[official])
        frame = frame.drop(columns=official)
    missing_all = [column for column in value_columns if frame[column].isna().all()]
    if missing_all:
        raise ValueError(f"Sample Adult DAT fields failed to merge for {year}: {missing_all}")
    return frame


def read_personsx_frame(year: int) -> tuple[pd.DataFrame | None, list[str]]:
    personsx_path = NHIS_RAW_ROOT / str(year) / "personsx.csv"
    if not personsx_path.exists():
        return None, []
    year_variable_map = YEAR_VARIABLE_MAP[year]
    personsx_map = YEAR_PERSONSX_MAP.get(year, {})
    header = pd.read_csv(personsx_path, nrows=0).columns.tolist()
    key_cols = [year_variable_map["hhx"], year_variable_map["fmx"], year_variable_map["person_key"]]
    available = []
    rename_map = {
        year_variable_map["hhx"]: "hhx",
        year_variable_map["fmx"]: "fmx",
        year_variable_map["person_key"]: "person_key",
    }
    for canonical_name, source_name in personsx_map.items():
        if source_name is not None and source_name in header:
            available.append(canonical_name)
            rename_map[source_name] = canonical_name
    # Also look for basic demographic/design vars in personsx as fallback
    # (for years where samadult omits these fields).
    for canonical_name, source_names in PERSONSX_FALLBACK_CANDIDATES.items():
        for source_name in source_names:
            if source_name in header and canonical_name not in rename_map.values():
                rename_map[source_name] = f"{canonical_name}_personsx"
                break
    for canonical_name, source_name in PERSONSX_SOURCE_COLUMNS.items():
        if source_name in header and canonical_name not in rename_map.values() and source_name not in rename_map:
            rename_map[source_name] = f"{canonical_name}_personsx"
    # Prefer official coverage recodes over hand-built HIKIND combinations.
    for canonical_name, candidates in INSURANCE_RECODE_CANDIDATES.items():
        for source_name in candidates:
            if source_name in header and source_name not in rename_map:
                rename_map[source_name] = canonical_name
                available.append(canonical_name)
                break
    if len(available) == 0 and len(rename_map) == 3:
        return None, []
    usecols = list(rename_map.keys())
    dtype_map = {
        year_variable_map["hhx"]: "string",
        year_variable_map["fmx"]: "string",
        year_variable_map["person_key"]: "string",
    }
    pframe = pd.read_csv(personsx_path, usecols=usecols, dtype=dtype_map, low_memory=False)
    pframe = pframe.rename(columns=rename_map)
    pframe["hhx"] = normalize_key_series(pframe["hhx"])
    pframe["fmx"] = normalize_key_series(pframe["fmx"])
    pframe["person_key"] = normalize_key_series(pframe["person_key"])
    if pframe.duplicated(subset=["hhx", "fmx", "person_key"]).any():
        raise ValueError(f"Duplicate personsx merge keys detected for {year}")
    return pframe, available


# Round 3c B2a: supplemental file readers for year-specific cholesterol variables
SUPPLEMENTAL_CHOLESTEROL_MAP: dict[int, tuple[str, str]] = {
    1998: ("prevadlt.csv", "CHLHIGH"),
}


def read_supplemental_cholesterol(year: int) -> pd.DataFrame | None:
    """Read cholesterol variable from supplemental files (prevadlt, cancerxx) when not in samadult."""
    if year not in SUPPLEMENTAL_CHOLESTEROL_MAP:
        return None
    filename, source_col = SUPPLEMENTAL_CHOLESTEROL_MAP[year]
    filepath = NHIS_RAW_ROOT / str(year) / filename
    if not filepath.exists():
        return None
    header = pd.read_csv(filepath, nrows=0).columns.tolist()
    if source_col not in header:
        return None
    year_variable_map = YEAR_VARIABLE_MAP[year]
    hhx_col = year_variable_map["hhx"]
    fmx_col = year_variable_map["fmx"]
    px_col = year_variable_map["person_key"]
    usecols = [hhx_col, fmx_col, px_col, source_col]
    dtype_map = {hhx_col: "string", fmx_col: "string", px_col: "string"}
    df = pd.read_csv(filepath, usecols=usecols, dtype=dtype_map, low_memory=False)
    df = df.rename(columns={
        hhx_col: "hhx",
        fmx_col: "fmx",
        px_col: "person_key",
        source_col: "cholesterol_high_ever",
    })
    df["hhx"] = normalize_key_series(df["hhx"])
    df["fmx"] = normalize_key_series(df["fmx"])
    df["person_key"] = normalize_key_series(df["person_key"])
    if df.duplicated(subset=["hhx", "fmx", "person_key"]).any():
        df = df.drop_duplicates(subset=["hhx", "fmx", "person_key"])
    return df


def merge_supplemental_cholesterol(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    """Merge supplemental cholesterol data for years where samadult lacks the variable."""
    supp = read_supplemental_cholesterol(year)
    if supp is None:
        return frame
    original_rows = len(frame)
    # Only merge if frame doesn't already have cholesterol data
    if "cholesterol_high_ever" in frame.columns and frame["cholesterol_high_ever"].notna().any():
        return frame
    # Drop placeholder column to avoid merge suffix conflicts
    if "cholesterol_high_ever" in frame.columns:
        frame = frame.drop(columns=["cholesterol_high_ever"])
    frame = frame.merge(supp, on=["hhx", "fmx", "person_key"], how="left", validate="one_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after supplemental cholesterol merge for {year}")
    return frame


def _read_family_csv_candidate(path: Path, source_poverty: str) -> pd.DataFrame | None:
    if not path.exists():
        return None
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            csv_members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
            if len(csv_members) != 1:
                return None
            with archive.open(csv_members[0]) as raw_handle:
                payload = raw_handle.read()
        source = io.BytesIO(payload)
    else:
        source = path
    header = pd.read_csv(source, nrows=0).columns.tolist()
    if hasattr(source, "seek"):
        source.seek(0)
    if source_poverty not in header:
        return None
    hhx_col = "HHX" if "HHX" in header else ("HH" if "HH" in header else None)
    fmx_col = "FMX" if "FMX" in header else None
    if hhx_col is None or fmx_col is None:
        return None
    return pd.read_csv(
        source,
        usecols=[hhx_col, fmx_col, source_poverty],
        dtype={hhx_col: "string", fmx_col: "string", source_poverty: "string"},
        low_memory=False,
    ).rename(columns={hhx_col: "hhx", fmx_col: "fmx", source_poverty: "poverty"})


def _read_family_fwf_poverty(year: int, source_poverty: str) -> pd.DataFrame | None:
    fwf_spec = YEAR_FAMILYXX_FWF_POVERTY_MAP.get(year)
    if fwf_spec is None or fwf_spec[0] != source_poverty:
        return None
    zip_path = NHIS_RAW_ROOT / str(year) / "familyxx.zip"
    if not zip_path.exists():
        raise FileNotFoundError(f"Required official Family archive missing for {year}: {zip_path}")
    with zipfile.ZipFile(zip_path) as archive:
        dat_members = [name for name in archive.namelist() if name.lower().endswith(".dat")]
        if len(dat_members) != 1:
            raise ValueError(f"Could not identify Family DAT member for {year}")
        positions = [(7, 12), (13, 14), fwf_spec[1]]
        colspecs = [(start - 1, end) for start, end in positions]
        with archive.open(dat_members[0]) as raw_handle:
            text_handle = io.TextIOWrapper(raw_handle, encoding="latin-1")
            return pd.read_fwf(
                text_handle,
                colspecs=colspecs,
                names=["hhx", "fmx", "poverty"],
                dtype="string",
            )


def read_familyxx_frame(year: int) -> pd.DataFrame:
    source_poverty = YEAR_FAMILYXX_MAP.get(year, {}).get("poverty")
    if source_poverty is None:
        raise ValueError(f"No poverty source mapping configured for {year}")
    year_dir = NHIS_RAW_ROOT / str(year)
    candidates = [year_dir / "familyxx.csv", year_dir / "ratcat.csv", year_dir / "familyxx.zip"]
    fframe = None
    for candidate in candidates:
        fframe = _read_family_csv_candidate(candidate, source_poverty)
        if fframe is not None:
            break
    if fframe is None:
        fframe = _read_family_fwf_poverty(year, source_poverty)
    if fframe is None:
        raise ValueError(f"Could not reconstruct keyed {source_poverty} poverty data for {year}")
    fframe["hhx"] = normalize_key_series(fframe["hhx"])
    fframe["fmx"] = normalize_key_series(fframe["fmx"])
    duplicate_keys = fframe.duplicated(subset=["hhx", "fmx"], keep=False)
    if duplicate_keys.any():
        conflicts = fframe.loc[duplicate_keys].groupby(["hhx", "fmx"])["poverty"].nunique(dropna=False)
        if conflicts.gt(1).any():
            raise ValueError(f"Conflicting duplicate family poverty keys detected for {year}")
        fframe = fframe.drop_duplicates(subset=["hhx", "fmx"])
    fframe["source_poverty_var"] = source_poverty
    return fframe


def read_househld_frame(year: int) -> tuple[pd.DataFrame | None, dict[str, str]]:
    househld_path = NHIS_RAW_ROOT / str(year) / "househld.csv"
    if not househld_path.exists():
        return None, {}
    header = pd.read_csv(househld_path, nrows=0).columns.tolist()
    year_variable_map = YEAR_VARIABLE_MAP[year]
    hhx_col = year_variable_map.get("hhx")
    if hhx_col not in header:
        hhx_col = "HHX" if "HHX" in header else None
    if hhx_col is None:
        return None, {}

    rename_map = {hhx_col: "hhx"}
    source_map: dict[str, str] = {}
    for canonical_name, source_names in HOUSEHLD_FALLBACK_CANDIDATES.items():
        for source_name in source_names:
            if source_name in header:
                rename_map[source_name] = f"{canonical_name}_househld"
                source_map[canonical_name] = source_name
                break
    if len(rename_map) == 1:
        return None, {}

    dtype_map = {hhx_col: "string"}
    hframe = pd.read_csv(househld_path, usecols=list(rename_map.keys()), dtype=dtype_map, low_memory=False)
    hframe = hframe.rename(columns=rename_map)
    hframe["hhx"] = normalize_key_series(hframe["hhx"])
    if hframe.duplicated(subset=["hhx"]).any():
        hframe = hframe.drop_duplicates(subset=["hhx"])
    return hframe, source_map


def merge_personsx(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    pframe, available = read_personsx_frame(year)
    if pframe is None:
        return frame
    original_rows = len(frame)
    # Determine overlap columns to avoid suffixes
    overlap = [
        c
        for c in [
            "education",
            "srh",
            "insurance_notcov",
            "insurance_cover",
            "insurance_cover65",
            "medicare_flag",
            "medicaid_flag",
            "private_flag",
            "other_public_flag",
            "military_flag",
        ]
        if c in frame.columns and c in pframe.columns
    ]
    frame = frame.drop(columns=overlap, errors="ignore")
    frame = frame.merge(pframe, on=["hhx", "fmx", "person_key"], how="left", validate="one_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after personsx merge for {year}")
    # Fallback demographic variables from personsx for supplement years
    for demo in ["age", "sex", "race", "hispanic", "education", "region", "marital", "design_strata", "design_psu"]:
        psx_col = f"{demo}_personsx"
        if psx_col in frame.columns:
            frame[demo] = frame[demo].where(frame[demo].notna(), frame[psx_col])
            frame = frame.drop(columns=[psx_col])
    return frame


def merge_househld(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    hframe, source_map = read_househld_frame(year)
    if hframe is None:
        return frame
    original_rows = len(frame)
    frame = frame.merge(hframe, on="hhx", how="left", validate="many_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after househld merge for {year}")
    for demo in ["intv_qrt", "region", "design_strata", "design_psu"]:
        hh_col = f"{demo}_househld"
        if hh_col in frame.columns:
            frame[demo] = frame[demo].where(frame[demo].notna(), frame[hh_col])
            frame = frame.drop(columns=[hh_col])
    if "design_strata" in source_map and "source_design_strata_var" in frame.columns:
        frame["source_design_strata_var"] = frame["source_design_strata_var"].where(
            frame["source_design_strata_var"].notna(), source_map["design_strata"]
        )
    if "design_psu" in source_map and "source_design_psu_var" in frame.columns:
        frame["source_design_psu_var"] = frame["source_design_psu_var"].where(
            frame["source_design_psu_var"].notna(), source_map["design_psu"]
        )
    return frame


def merge_cancerxx_weight(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    """Recover sample-adult final weight when it is stored in cancerxx.

    Some local NHIS CSV extracts for later years keep the true sample-adult
    final annual weight outside ``samadult.csv`` while ``samadult.csv`` has
    ``AWEIGHTP`` (self-reported body weight in pounds). Never use AWEIGHTP as
    a survey weight; fill ``wtfa_sa`` only from a real ``WTFA_SA`` column.
    """
    cancerxx_path = NHIS_RAW_ROOT / str(year) / "cancerxx.csv"
    if not cancerxx_path.exists():
        return frame
    header = pd.read_csv(cancerxx_path, nrows=0).columns.tolist()
    if "WTFA_SA" not in header:
        return frame
    year_variable_map = YEAR_VARIABLE_MAP[year]
    key_cols = [year_variable_map["hhx"], year_variable_map["fmx"], year_variable_map["person_key"]]
    if any(col not in header for col in key_cols):
        return frame
    dtype_map = {
        year_variable_map["hhx"]: "string",
        year_variable_map["fmx"]: "string",
        year_variable_map["person_key"]: "string",
    }
    wframe = pd.read_csv(cancerxx_path, usecols=key_cols + ["WTFA_SA"], dtype=dtype_map, low_memory=False)
    wframe = wframe.rename(
        columns={
            year_variable_map["hhx"]: "hhx",
            year_variable_map["fmx"]: "fmx",
            year_variable_map["person_key"]: "person_key",
            "WTFA_SA": "wtfa_sa_cancerxx",
        }
    )
    wframe["hhx"] = normalize_key_series(wframe["hhx"])
    wframe["fmx"] = normalize_key_series(wframe["fmx"])
    wframe["person_key"] = normalize_key_series(wframe["person_key"])
    if wframe.duplicated(subset=["hhx", "fmx", "person_key"]).any():
        raise ValueError(f"Duplicate cancerxx weight merge keys detected for {year}")
    original_rows = len(frame)
    frame = frame.merge(wframe, on=["hhx", "fmx", "person_key"], how="left", validate="one_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after cancerxx weight merge for {year}")
    frame["wtfa_sa"] = frame["wtfa_sa"].where(frame["wtfa_sa"].notna(), frame["wtfa_sa_cancerxx"])
    frame = frame.drop(columns=["wtfa_sa_cancerxx"])
    return frame


def merge_familyxx(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    fframe = read_familyxx_frame(year)
    original_rows = len(frame)
    frame = frame.drop(columns=["poverty", "source_poverty_var"], errors="ignore")
    frame = frame.merge(fframe, on=["hhx", "fmx"], how="left", validate="many_to_one")
    if len(frame) != original_rows:
        raise ValueError(f"Row count changed after familyxx merge for {year}")
    if frame["source_poverty_var"].isna().any():
        missing = int(frame["source_poverty_var"].isna().sum())
        raise ValueError(f"Family poverty key failed to merge for {missing} Sample Adult rows in {year}")
    return frame


def derive_insurance_type(frame: pd.DataFrame) -> pd.DataFrame:
    """Derive coverage hierarchy without imputing unobserved plan type."""
    frame["insurance_type"] = pd.Series(pd.NA, index=frame.index, dtype="string")

    age = pd.to_numeric(frame["age"], errors="coerce")
    year = pd.to_numeric(frame["year"], errors="coerce")
    notcov = pd.to_numeric(frame.get("insurance_notcov", pd.Series(pd.NA, index=frame.index)), errors="coerce")
    medicare = frame.get("medicare_flag", pd.Series(pd.NA, index=frame.index, dtype="boolean"))
    medicaid = frame.get("medicaid_flag", pd.Series(pd.NA, index=frame.index, dtype="boolean"))
    private = frame.get("private_flag", pd.Series(pd.NA, index=frame.index, dtype="boolean"))
    other_public = frame.get("other_public_flag", pd.Series(pd.NA, index=frame.index, dtype="boolean"))
    military = frame.get("military_flag", pd.Series(pd.NA, index=frame.index, dtype="boolean"))

    # The official hierarchy is available from 2015 onward.
    cover65 = pd.to_numeric(frame.get("insurance_cover65", pd.Series(pd.NA, index=frame.index)), errors="coerce")
    cover = pd.to_numeric(frame.get("insurance_cover", pd.Series(pd.NA, index=frame.index)), errors="coerce")
    m65 = year.ge(2015) & age.ge(65)
    mlt65 = year.ge(2015) & age.lt(65)

    frame.loc[m65 & cover65.eq(1), "insurance_type"] = "private_only"
    frame.loc[m65 & cover65.eq(2), "insurance_type"] = "medicare_dual"
    frame.loc[m65 & cover65.eq(3), "insurance_type"] = "medicare_private"
    frame.loc[m65 & cover65.eq(4), "insurance_type"] = "medicare_only"
    frame.loc[m65 & cover65.eq(5), "insurance_type"] = "other"
    frame.loc[m65 & cover65.eq(6), "insurance_type"] = "uninsured"

    frame.loc[mlt65 & cover.eq(1), "insurance_type"] = "private_only"
    frame.loc[mlt65 & cover.eq(2), "insurance_type"] = "medicaid_only"
    frame.loc[mlt65 & cover.eq(3), "insurance_type"] = "other"
    frame.loc[mlt65 & cover.eq(4), "insurance_type"] = "uninsured"

    # Explicit official unknown codes must remain missing, not be overwritten
    # by demographic assumptions or a modal-category backstop.
    official_unknown = ((m65 & cover65.eq(7)) | (mlt65 & cover.eq(5))).fillna(False)
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & notcov.eq(1), "insurance_type"] = "uninsured"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & medicare.eq(True) & medicaid.eq(True), "insurance_type"] = "medicare_dual"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & medicare.eq(True) & private.eq(True), "insurance_type"] = "medicare_private"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & medicare.eq(True), "insurance_type"] = "medicare_only"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & private.eq(True), "insurance_type"] = "private_only"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & medicaid.eq(True), "insurance_type"] = "medicaid_only"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & (other_public.eq(True) | military.eq(True)), "insurance_type"] = "other"
    manual = frame["insurance_type"].isna() & ~official_unknown
    frame.loc[manual & notcov.eq(2), "insurance_type"] = "other"

    return frame


def clean_year_frame(frame: pd.DataFrame, year: int) -> pd.DataFrame:
    frame["canev"] = binary_with_missing(frame["canev"], {1}, {2})
    frame["dm_ever"] = binary_with_missing(frame["dm_ever"], {1}, {2, 3})
    frame["dm_med_pills"] = binary_with_missing(frame["dm_med_pills"], {1}, {2})
    frame["dm_med_insulin"] = binary_with_missing(frame["dm_med_insulin"], {1}, {2})
    for age_column in ["canage7", "canage8", "canage9", "canage13", "canage19", "canage21", "canage25", "dm_age"]:
        frame[age_column] = numeric_with_special_missing(frame[age_column])
    frame["dm_duration"] = clean_numeric(frame["dm_duration"])
    frame["bmi"] = clean_bmi(frame["bmi"])
    frame["smoking_status"] = clean_numeric(frame["smoking_status"])
    frame["wtfa_sa"] = clean_numeric(frame["wtfa_sa"])
    frame["design_strata"] = clean_numeric(frame["design_strata"])
    frame["design_psu"] = clean_numeric(frame["design_psu"])
    for column in ["age", "sex", "race", "hispanic", "education", "region", "poverty", "marital", "smoking_ever", "afford_rx", "intv_qrt"]:
        frame[column] = clean_numeric(frame[column])
    frame["race_structural_missing"] = False
    if year in STRUCTURAL_UNKNOWN_RACE_YEARS:
        frame["race_structural_missing"] = frame["race"].isna() | frame["hispanic"].isna()
        frame["race"] = frame["race"].where(frame["race"].notna(), 99)
        frame["hispanic"] = frame["hispanic"].where(frame["hispanic"].notna(), 99)
    frame = derive_race_ethnicity_4cat(frame)
    # Stabilize the legacy model adjustment: rare public-use race codes share
    # an Other/unknown category. The mutually exclusive race/ethnicity field is
    # retained separately for sensitivity analysis.
    # structurally unavailable race years share an Other/unknown category.
    frame["race"] = frame["race"].where(frame["race"].isin([1, 2]), 97)
    # New comorbidity variables
    frame["hypertension_ever"] = binary_with_missing(frame["hypertension_ever"], {1}, {2})
    frame["chd_ever"] = binary_with_missing(frame["chd_ever"], {1}, {2})
    frame["stroke_ever"] = binary_with_missing(frame["stroke_ever"], {1}, {2})
    frame["cholesterol_high_ever"] = binary_with_missing(frame["cholesterol_high_ever"], {1}, {2})
    # Physical activity — 95/96 mean never/unable (zero sessions), whereas
    # 97-99 and out-of-range values are nonresponse and stay missing.
    frame["phys_mod_freq"] = clean_activity_frequency(frame["phys_mod_freq"])
    frame["phys_vig_freq"] = clean_activity_frequency(frame["phys_vig_freq"])
    mod_ok = frame["phys_mod_freq"].notna()
    vig_ok = frame["phys_vig_freq"].notna()
    phys_any = frame["phys_mod_freq"].ge(3) | frame["phys_vig_freq"].ge(3)
    frame["phys_active_any"] = phys_any.where(mod_ok | vig_ok)
    frame = derive_alcohol_status(frame, year)
    frame["cost_barrier_rx"] = binary_with_missing(frame["afford_rx"], {1}, {2})
    for column in ["medicare_flag", "medicaid_flag", "private_flag"]:
        frame[column] = insurance_plan_recode_boolean(frame[column])
    frame["other_public_flag"] = other_public_recode_boolean(frame["other_public_flag"])
    frame["military_flag"] = military_recode_boolean(frame["military_flag"], year)
    # Self-rated health
    frame["srh"] = clean_numeric(frame["srh"])
    frame = derive_insurance_type(frame)
    # Site flags
    for site_name, source_column in GI_SITE_FLAGS.items():
        flag_column = SITE_FLAG_COLUMNS[site_name]
        frame[flag_column] = pd.to_numeric(frame[source_column.lower()], errors="coerce").eq(1)
        frame[source_column.lower()] = frame[flag_column]
    frame["gi_any"] = frame[list(SITE_FLAG_COLUMNS.values())].any(axis=1)
    frame["smoking_3cat"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    frame.loc[frame["smoking_status"].eq(4), "smoking_3cat"] = "never"
    frame.loc[frame["smoking_status"].isin([2, 3]), "smoking_3cat"] = "former"
    frame.loc[frame["smoking_status"].eq(1), "smoking_3cat"] = "current"
    frame["dm_med_intensity"] = pd.Series("unknown", index=frame.index, dtype="string")
    frame.loc[frame["dm_med_insulin"].eq(True), "dm_med_intensity"] = "insulin_containing"
    frame.loc[frame["dm_med_pills"].eq(True) & frame["dm_med_insulin"].eq(False), "dm_med_intensity"] = "pills_only"
    frame.loc[frame["dm_med_pills"].eq(False) & frame["dm_med_insulin"].eq(False), "dm_med_intensity"] = "none"
    frame["publicid"] = build_publicid(frame)
    # Derived categorizations
    frame = derive_education_5cat(frame)
    frame = derive_education_4cat(frame)
    frame = derive_poverty_4cat(frame)
    frame = derive_poverty_3cat(frame)
    return frame



def harmonize_year(year: int) -> pd.DataFrame:
    frame = read_year_frame(year)
    frame = merge_samadult_fwf_supplement(frame, year)
    frame = merge_personsx(frame, year)
    frame = merge_househld(frame, year)
    frame = merge_cancerxx_weight(frame, year)
    frame = merge_supplemental_cholesterol(frame, year)
    frame = merge_familyxx(frame, year)
    frame = clean_year_frame(frame, year)
    frame = add_missing_columns(frame)
    for column in DERIVED_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA
    return frame[FINAL_COLUMNS]


def build_metadata(frame: pd.DataFrame) -> dict[str, object]:
    row_count_by_year = {str(year): int(count) for year, count in frame.groupby("year").size().sort_index().items()}
    non_missing_rate_by_column = {
        column: round(float(rate), 6)
        for column, rate in frame.notna().mean().items()
    }
    gi_positive_count_by_site = {
        "gi_any": int(frame["gi_any"].fillna(False).sum()),
        "colon_flag": int(frame["colon_flag"].fillna(False).sum()),
        "esoph_flag": int(frame["esoph_flag"].fillna(False).sum()),
        "gallbladder_flag": int(frame["gallbladder_flag"].fillna(False).sum()),
        "liver_flag": int(frame["liver_flag"].fillna(False).sum()),
        "pancreas_flag": int(frame["pancreas_flag"].fillna(False).sum()),
        "rectum_flag": int(frame["rectum_flag"].fillna(False).sum()),
        "stomach_flag": int(frame["stomach_flag"].fillna(False).sum()),
    }
    new_variable_coverage = {}
    for var in ["hypertension_ever", "cholesterol_high_ever", "chd_ever", "stroke_ever",
                "alcohol_status", "phys_active_any", "insurance_type", "srh",
                "education_5cat", "poverty_4cat"]:
        if var in frame.columns:
            by_year = frame.groupby("year")[var].apply(lambda s: s.notna().mean())
            new_variable_coverage[var] = {str(y): round(float(v), 4) for y, v in by_year.items()}
    return {
        "pipeline_version": HARMONIZATION_VERSION,
        "analysis_years": ANALYSIS_YEARS,
        "explicitly_excluded_primary_years": OPTIONAL_SENSITIVITY_YEARS,
        "row_count_by_year": row_count_by_year,
        "total_rows": int(len(frame)),
        "non_missing_rate_by_column": non_missing_rate_by_column,
        "gi_positive_count_by_site": gi_positive_count_by_site,
        "new_variable_coverage": new_variable_coverage,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_parquet_atomic(frame: pd.DataFrame, path: Path) -> None:
    """Write a complete parquet file before replacing any prior v4 output."""
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
    HARMONIZED_DIR.mkdir(parents=True, exist_ok=True)
    yearly_frames = [harmonize_year(year) for year in ANALYSIS_YEARS]
    harmonized_frame = pd.concat(yearly_frames, ignore_index=True)
    metadata = build_metadata(harmonized_frame)
    write_parquet_atomic(harmonized_frame, HARMONIZED_PATH)
    metadata["harmonized_parquet_sha256"] = sha256_file(HARMONIZED_PATH)
    with METADATA_PATH.open("w", encoding="utf-8") as metadata_file:
        json.dump(metadata, metadata_file, indent=2, ensure_ascii=False)
    print(f"Wrote {HARMONIZED_PATH}")
    print(f"Wrote {METADATA_PATH}")
    print(f"Total rows: {metadata['total_rows']}")
    print(f"GI cancer count: {metadata['gi_positive_count_by_site']['gi_any']}")
    print(f"Stomach count: {metadata['gi_positive_count_by_site']['stomach_flag']}")
    print("Row counts by year:")
    for year, count in metadata["row_count_by_year"].items():
        print(f"  {year}: {count}")
    print("New variable coverage (fraction non-missing) by year:")
    for var, cov in metadata.get("new_variable_coverage", {}).items():
        covered_years = sum(1 for v in cov.values() if v >= 0.8)
        print(f"  {var}: {covered_years} years with >=80% coverage")


if __name__ == "__main__":
    main()
