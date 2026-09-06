import json
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import *
from parse_nhis_lmf import sanitize_dodqtr_series

COHORT_OUTPUT_DIR = OUTPUTS_DIR / "cohort"
PIPELINE_VERSION = os.environ.get("NHIS_GI_PIPELINE_VERSION", "v4inc2007")
COHORT_OUTPUT_PATH = COHORT_OUTPUT_DIR / f"analytic_cohort_{PIPELINE_VERSION}.parquet"
VARIABLE_COVERAGE_PATH = COHORT_OUTPUT_DIR / f"variable_coverage_{PIPELINE_VERSION}.md"
COHORT_FLOW_PATH = COHORT_OUTPUT_DIR / f"cohort_flow_{PIPELINE_VERSION}.json"
TRAJECTORY_DISTRIBUTION_PATH = COHORT_OUTPUT_DIR / f"trajectory_distribution_{PIPELINE_VERSION}.csv"
TRAJECTORY_DISTRIBUTION_6CAT_PATH = COHORT_OUTPUT_DIR / f"trajectory_distribution_6cat_{PIPELINE_VERSION}.csv"
NHIS_HARMONIZED_PATH = OUTPUTS_DIR / "harmonized" / f"nhis_gi_diabetes_harmonized_{PIPELINE_VERSION}.parquet"
LMF_PARSED_PATH = OUTPUTS_DIR / "lmf" / f"nhis_lmf_parsed_{PIPELINE_VERSION}.parquet"

SITE_FLAG_COLUMNS = [
    "colon_flag",
    "esoph_flag",
    "gallbladder_flag",
    "liver_flag",
    "pancreas_flag",
    "rectum_flag",
    "stomach_flag",
]

GI_AGE_COLUMNS = ["canage7", "canage8", "canage9", "canage13", "canage19", "canage21", "canage25"]

TRAJECTORY_ORDER = [
    "gi_only",
    "dm_to_gi_lt2y",
    "dm_to_gi_2_10y",
    "dm_to_gi_gt10y",
    "gi_to_dm",
    "same_year",
    "dm_order_unknown",
]

TRAJECTORY_6CAT_ORDER = [
    "gi_only",
    "peri_diagnostic",
    "dm_to_gi_2_10y",
    "dm_to_gi_gt10y",
    "gi_to_dm",
    "dm_order_unknown",
]

OUTCOME_COLUMNS = [
    "death_allcause_5y",
    "death_allcause_10y",
    "death_cancer_10y",
    "death_dm_contributing_10y",
]


def normalize_publicid(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip()
    return normalized.replace({"": pd.NA})


def normalize_boolean_series(series: pd.Series) -> pd.Series:
    numeric_values = pd.to_numeric(series, errors="coerce")
    string_values = series.astype("string").str.strip().str.lower()
    normalized = pd.Series(pd.NA, index=series.index, dtype="boolean")
    normalized.loc[numeric_values.eq(1)] = True
    normalized.loc[numeric_values.eq(0)] = False
    normalized.loc[string_values.eq("true")] = True
    normalized.loc[string_values.eq("false")] = False
    return normalized


def _format_identifier_value(value):
    if pd.isna(value):
        return pd.NA
    numeric_value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.notna(numeric_value):
        if float(numeric_value).is_integer():
            return str(int(numeric_value))
        return str(float(numeric_value))
    text_value = str(value).strip()
    return text_value if text_value else pd.NA


def derive_design_identifiers(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    years = pd.to_numeric(cohort_frame["year"], errors="raise")
    if not years.between(1997, 2018).all():
        raise ValueError("NHIS design-period identifiers support 1997–2018 only")
    year_values = pd.Series("2016-2018", index=cohort_frame.index, dtype="string")
    year_values.loc[years <= 2015] = "2006-2015"
    year_values.loc[years <= 2005] = "1997-2005"
    strata_values = cohort_frame["design_strata"].map(_format_identifier_value)
    psu_values = cohort_frame["design_psu"].map(_format_identifier_value)
    cohort_frame["design_strata_prefixed"] = pd.Series(pd.NA, index=cohort_frame.index, dtype="string")
    cohort_frame["design_psu_prefixed"] = pd.Series(pd.NA, index=cohort_frame.index, dtype="string")
    strata_mask = year_values.notna() & strata_values.notna()
    psu_mask = year_values.notna() & strata_values.notna() & psu_values.notna()
    cohort_frame.loc[strata_mask, "design_strata_prefixed"] = year_values[strata_mask] + "." + strata_values[strata_mask]
    cohort_frame.loc[psu_mask, "design_psu_prefixed"] = (
        year_values[psu_mask] + "." + strata_values[psu_mask] + "." + psu_values[psu_mask]
    )
    return cohort_frame


def load_input_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    nhis_frame = pd.read_parquet(NHIS_HARMONIZED_PATH)
    lmf_frame = pd.read_parquet(LMF_PARSED_PATH)
    nhis_frame["publicid"] = normalize_publicid(nhis_frame["publicid"])
    lmf_frame["publicid"] = normalize_publicid(lmf_frame["publicid"])
    if nhis_frame["publicid"].isna().any():
        raise ValueError("Missing publicid values detected in harmonized NHIS data")
    if nhis_frame["publicid"].str.len().ne(14).any():
        raise ValueError("Harmonized NHIS publicid values are not all 14 characters")
    if nhis_frame["publicid"].duplicated().any():
        raise ValueError("Duplicate publicid values detected in harmonized NHIS data")
    if lmf_frame["publicid"].duplicated().any():
        raise ValueError("Duplicate publicid values detected in parsed LMF data")
    return nhis_frame, lmf_frame


def coerce_input_types(nhis_frame: pd.DataFrame, lmf_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    nhis_numeric_columns = ["year", "intv_qrt", "age", "dm_age", "dm_duration", "wtfa_sa"] + GI_AGE_COLUMNS
    lmf_numeric_columns = ["eligstat", "mortstat", "ucod_leading", "diabetes", "dodqtr", "dodyear", "sa_wgt_new", "survey_year"]
    for column in nhis_numeric_columns:
        nhis_frame[column] = pd.to_numeric(nhis_frame[column], errors="coerce")
    for column in lmf_numeric_columns:
        lmf_frame[column] = pd.to_numeric(lmf_frame[column], errors="coerce")
    lmf_frame["dodqtr"] = sanitize_dodqtr_series(lmf_frame["dodqtr"])
    nhis_boolean_columns = ["gi_any", "dm_ever"] + SITE_FLAG_COLUMNS
    for column in nhis_boolean_columns:
        nhis_frame[column] = normalize_boolean_series(nhis_frame[column])
    # Round 3: coerce new numeric columns
    for column in ["phys_mod_freq", "phys_vig_freq", "srh"]:
        if column in nhis_frame.columns:
            nhis_frame[column] = pd.to_numeric(nhis_frame[column], errors="coerce")
    return nhis_frame, lmf_frame


def merge_nhis_lmf(nhis_frame: pd.DataFrame, lmf_frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int | float]]:
    merged_frame = nhis_frame.merge(lmf_frame, on="publicid", how="left", indicator=True, validate="one_to_one")
    id_present_mask = merged_frame["_merge"].eq("both")
    eligible_mask = id_present_mask & merged_frame["eligstat"].eq(1)
    death_mask = eligible_mask & merged_frame["mortstat"].eq(1)
    layers: dict[str, int | float] = {
        "sample_adult_n": int(len(nhis_frame)),
        "lmf_publicid_present_n": int(id_present_mask.sum()),
        "eligstat_1_n": int(eligible_mask.sum()),
        "mortstat_1_among_eligible_n": int(death_mask.sum()),
        "lmf_publicid_present_rate": float(id_present_mask.mean()),
        "eligstat_1_rate_of_sample_adult": float(eligible_mask.mean()),
        "eligstat_1_rate_among_id_present": float(eligible_mask.sum() / id_present_mask.sum()),
        "death_rate_among_eligible": float(death_mask.sum() / eligible_mask.sum()),
    }
    merged_frame = merged_frame.loc[eligible_mask].drop(columns=["_merge"]).reset_index(drop=True)
    if not merged_frame["eligstat"].eq(1).all():
        raise ValueError("Merged analytic frame contains non-eligible LMF rows")
    year_mismatch_count = int(merged_frame["year"].ne(merged_frame["survey_year"]).sum())
    if year_mismatch_count != 0:
        raise ValueError(f"Year mismatch detected between NHIS and LMF for {year_mismatch_count} rows")
    return merged_frame, layers


def derive_trajectory_fields(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    interview_age = pd.to_numeric(cohort_frame["age"], errors="coerce")
    dm_dx_age_raw = pd.to_numeric(cohort_frame["dm_age"], errors="coerce")
    impossible_dm_age_mask = dm_dx_age_raw.lt(0) | (interview_age.notna() & dm_dx_age_raw.gt(interview_age))
    gi_dx_age_frame = cohort_frame[GI_AGE_COLUMNS].apply(lambda column: pd.to_numeric(column, errors="coerce"))
    impossible_gi_age_frame = gi_dx_age_frame.lt(0) | gi_dx_age_frame.gt(interview_age, axis=0)
    cohort_frame["flag_impossible_dm_age"] = normalize_boolean_series(impossible_dm_age_mask)
    cohort_frame["flag_impossible_gi_age"] = normalize_boolean_series(impossible_gi_age_frame.any(axis=1))
    cohort_frame["gi_first_dx_age"] = gi_dx_age_frame.mask(impossible_gi_age_frame).min(axis=1, skipna=True)
    cohort_frame["dm_dx_age"] = dm_dx_age_raw.mask(impossible_dm_age_mask)
    cohort_frame["gi_dm_lag_years"] = cohort_frame["gi_first_dx_age"] - cohort_frame["dm_dx_age"]
    cohort_frame["trajectory_7cat"] = pd.Series("dm_order_unknown", index=cohort_frame.index, dtype="string")
    gi_only_mask = cohort_frame["dm_ever"].eq(False)
    orderable_mask = cohort_frame["dm_ever"].eq(True) & cohort_frame["dm_dx_age"].notna() & cohort_frame["gi_first_dx_age"].notna()
    lag_values = cohort_frame["gi_dm_lag_years"]
    cohort_frame.loc[gi_only_mask, "trajectory_7cat"] = "gi_only"
    cohort_frame.loc[orderable_mask & lag_values.gt(0) & lag_values.lt(2), "trajectory_7cat"] = "dm_to_gi_lt2y"
    cohort_frame.loc[orderable_mask & lag_values.ge(2) & lag_values.le(10), "trajectory_7cat"] = "dm_to_gi_2_10y"
    cohort_frame.loc[orderable_mask & lag_values.gt(10), "trajectory_7cat"] = "dm_to_gi_gt10y"
    cohort_frame.loc[orderable_mask & lag_values.lt(0), "trajectory_7cat"] = "gi_to_dm"
    cohort_frame.loc[orderable_mask & lag_values.eq(0), "trajectory_7cat"] = "same_year"
    cohort_frame["trajectory_7cat"] = pd.Categorical(cohort_frame["trajectory_7cat"], categories=TRAJECTORY_ORDER, ordered=True)
    # Round 3: derive 6-category trajectory (merge <2y and same_year into peri_diagnostic)
    cohort_frame["trajectory_6cat"] = cohort_frame["trajectory_7cat"].astype(str)
    cohort_frame.loc[cohort_frame["trajectory_6cat"].isin(["dm_to_gi_lt2y", "same_year"]), "trajectory_6cat"] = "peri_diagnostic"
    cohort_frame["trajectory_6cat"] = pd.Categorical(cohort_frame["trajectory_6cat"], categories=TRAJECTORY_6CAT_ORDER, ordered=True)
    return cohort_frame


def derive_followup_and_outcomes(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    cohort_frame["interview_time"] = cohort_frame["year"] + ((cohort_frame["intv_qrt"] - 1) * 0.25)
    cohort_frame.loc[cohort_frame["intv_qrt"].isna(), "interview_time"] = cohort_frame.loc[cohort_frame["intv_qrt"].isna(), "year"] + 0.5
    cohort_frame["dodqtr"] = sanitize_dodqtr_series(cohort_frame["dodqtr"])
    cohort_frame["death_time"] = cohort_frame["dodyear"] + ((cohort_frame["dodqtr"] - 1) * 0.25)
    death_mask = cohort_frame["mortstat"].eq(1)
    cohort_frame.loc[~death_mask, "death_time"] = pd.NA
    cohort_frame["censor_time"] = 2019.99
    cohort_frame["followup_years"] = cohort_frame["censor_time"] - cohort_frame["interview_time"]
    cohort_frame.loc[death_mask, "followup_years"] = cohort_frame.loc[death_mask, "death_time"] - cohort_frame.loc[death_mask, "interview_time"]
    invalid_followup_mask = cohort_frame["followup_years"].isna() | cohort_frame["followup_years"].lt(0)
    censor_before_5y_mask = ~death_mask & cohort_frame["followup_years"].lt(5)
    censor_before_10y_mask = ~death_mask & cohort_frame["followup_years"].lt(10)

    cohort_frame["death_allcause_5y"] = pd.Series(0.0, index=cohort_frame.index)
    cohort_frame.loc[death_mask & cohort_frame["followup_years"].le(5), "death_allcause_5y"] = 1.0
    cohort_frame.loc[censor_before_5y_mask | invalid_followup_mask, "death_allcause_5y"] = pd.NA

    cohort_frame["death_allcause_10y"] = pd.Series(0.0, index=cohort_frame.index)
    cohort_frame.loc[death_mask & cohort_frame["followup_years"].le(10), "death_allcause_10y"] = 1.0
    cohort_frame.loc[censor_before_10y_mask | invalid_followup_mask, "death_allcause_10y"] = pd.NA

    cohort_frame["death_cancer_10y"] = pd.Series(0.0, index=cohort_frame.index)
    cohort_frame.loc[
        death_mask & cohort_frame["ucod_leading"].eq(2) & cohort_frame["followup_years"].le(10),
        "death_cancer_10y",
    ] = 1.0
    cohort_frame.loc[censor_before_10y_mask | invalid_followup_mask, "death_cancer_10y"] = pd.NA

    cohort_frame["death_dm_contributing_10y"] = pd.Series(0.0, index=cohort_frame.index)
    cohort_frame.loc[
        death_mask & cohort_frame["diabetes"].eq(1) & cohort_frame["followup_years"].le(10),
        "death_dm_contributing_10y",
    ] = 1.0
    cohort_frame.loc[censor_before_10y_mask | invalid_followup_mask, "death_dm_contributing_10y"] = pd.NA

    return cohort_frame


def derive_weights_and_flags(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    pooled_years = len(ANALYSIS_YEARS)
    cohort_frame["sa_wgt_pool"] = cohort_frame["sa_wgt_new"] / pooled_years
    cohort_frame["wtfa_sa_pool"] = cohort_frame["wtfa_sa"] / pooled_years
    if "flag_impossible_dm_age" not in cohort_frame.columns:
        cohort_frame["flag_impossible_dm_age"] = normalize_boolean_series(cohort_frame["dm_dx_age"].gt(cohort_frame["age"]))
    if "flag_impossible_gi_age" not in cohort_frame.columns:
        cohort_frame["flag_impossible_gi_age"] = normalize_boolean_series(cohort_frame["gi_first_dx_age"].gt(cohort_frame["age"]))
    cohort_frame["flag_negative_followup"] = normalize_boolean_series(cohort_frame["followup_years"].lt(0))
    return cohort_frame


def build_trajectory_distribution(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    distribution_series = cohort_frame["trajectory_7cat"].value_counts(dropna=False).reindex(TRAJECTORY_ORDER, fill_value=0)
    distribution_frame = distribution_series.rename_axis("trajectory_7cat").reset_index(name="n")
    distribution_frame["n"] = distribution_frame["n"].astype(int)
    return distribution_frame


def build_trajectory_distribution_6cat(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    distribution_series = cohort_frame["trajectory_6cat"].value_counts(dropna=False).reindex(TRAJECTORY_6CAT_ORDER, fill_value=0)
    distribution_frame = distribution_series.rename_axis("trajectory_6cat").reset_index(name="n")
    distribution_frame["n"] = distribution_frame["n"].astype(int)
    return distribution_frame


def derive_missing_indicators(cohort_frame: pd.DataFrame) -> pd.DataFrame:
    # Round 3: encode high-missingness vars with explicit "missing" category to avoid listwise deletion collapse
    new_vars = ["hypertension_ever", "cholesterol_high_ever", "chd_ever", "stroke_ever",
                "alcohol_status", "phys_active_any", "insurance_type", "srh",
                "education_5cat", "education_4cat", "poverty_4cat", "poverty_3cat"]
    for var in new_vars:
        if var not in cohort_frame.columns:
            continue
        miss_rate = cohort_frame[var].isna().mean()
        if miss_rate > 0.15:
            cohort_frame[var] = cohort_frame[var].astype("string").fillna("missing")
    return cohort_frame


def build_variable_coverage_report(cohort_frame: pd.DataFrame) -> str:
    lines = ["# Variable Coverage Report (Round 3b)", ""]
    lines.append("This report shows **raw missing rate** (before missing-indicator encoding) and **effective category rate** (non-missing, excluding explicit 'missing' category).")
    lines.append("")
    new_vars = ["hypertension_ever", "cholesterol_high_ever", "chd_ever", "stroke_ever",
                "alcohol_status", "phys_active_any", "insurance_type", "srh",
                "education_5cat", "education_4cat", "poverty_4cat", "poverty_3cat"]
    for var in new_vars:
        if var not in cohort_frame.columns:
            continue
        lines.append(f"## {var}")
        lines.append("")
        # Compute raw missing rate (before missing-indicator encoding)
        raw_series = cohort_frame[var]
        if raw_series.dtype == "object":
            raw_missing = ((raw_series.isna()) | (raw_series == "missing")).mean()
            raw_valid = 1 - raw_missing
        else:
            raw_valid = raw_series.notna().mean()
        lines.append(f"**Overall raw valid (non-missing) rate: {raw_valid:.2%}**")
        lines.append("")
        by_year = cohort_frame.groupby("year")[var].apply(lambda s: s.notna().mean() if s.dtype != "object" else ((s.notna()) & (s != "missing")).mean())
        lines.append("| Year | Raw Valid Rate | Missing-Indicator Effective |")
        lines.append("|------|----------------|----------------------------|")
        for year, cov in by_year.items():
            eff = cohort_frame.loc[cohort_frame["year"] == year, var].notna().mean()
            lines.append(f"| {year} | {cov:.2%} | {eff:.2%} |")
        # Stratified by trajectory
        lines.append("")
        lines.append(f"### Stratified by trajectory_6cat")
        lines.append("")
        strat_raw = cohort_frame.groupby("trajectory_6cat", observed=True).apply(
            lambda df: ((df[var].notna()) & (df[var] != "missing")).mean() if df[var].dtype == "object" else df[var].notna().mean(),
            include_groups=False
        ).reset_index()
        strat_raw.columns = ["trajectory", f"{var}_raw_valid"]
        lines.append(strat_raw.to_markdown(index=False))
        lines.append("")
    return "\n".join(lines)


def build_outcome_event_counts(cohort_frame: pd.DataFrame) -> dict[str, dict[str, int]]:
    outcome_event_counts = {}
    for column in OUTCOME_COLUMNS:
        outcome_event_counts[column] = {
            "count_1": int(cohort_frame[column].eq(1).sum()),
            "count_0": int(cohort_frame[column].eq(0).sum()),
            "count_nan": int(cohort_frame[column].isna().sum()),
        }
    return outcome_event_counts


def build_flow_payload(
    linkage_layers: dict[str, int | float],
    cohort_frame: pd.DataFrame,
    trajectory_distribution_frame: pd.DataFrame,
    trajectory_distribution_6cat_frame: pd.DataFrame,
) -> dict[str, object]:
    site_counts = {column: int(cohort_frame[column].eq(True).sum()) for column in SITE_FLAG_COLUMNS}
    impossible_pattern_counts = {
        "flag_impossible_dm_age": int(cohort_frame["flag_impossible_dm_age"].eq(True).sum()),
        "flag_impossible_gi_age": int(cohort_frame["flag_impossible_gi_age"].eq(True).sum()),
        "flag_negative_followup": int(cohort_frame["flag_negative_followup"].eq(True).sum()),
    }
    trajectory_distribution = {row["trajectory_7cat"]: int(row["n"]) for _, row in trajectory_distribution_frame.iterrows()}
    trajectory_distribution_6cat = {row["trajectory_6cat"]: int(row["n"]) for _, row in trajectory_distribution_6cat_frame.iterrows()}
    return {
        "pipeline_version": PIPELINE_VERSION,
        "analysis_years": ANALYSIS_YEARS,
        "linkage_layers": linkage_layers,
        "gi_cohort_n": int(len(cohort_frame)),
        "stomach_subgroup_n": int(cohort_frame["stomach_flag"].eq(True).sum()),
        "site_counts": site_counts,
        "trajectory_distribution": trajectory_distribution,
        "trajectory_distribution_6cat": trajectory_distribution_6cat,
        "outcome_event_counts": build_outcome_event_counts(cohort_frame),
        "impossible_pattern_counts": impossible_pattern_counts,
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
    COHORT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    nhis_frame, lmf_frame = load_input_frames()
    nhis_frame, lmf_frame = coerce_input_types(nhis_frame, lmf_frame)
    merged_frame, linkage_layers = merge_nhis_lmf(nhis_frame, lmf_frame)
    cohort_frame = merged_frame.loc[merged_frame["gi_any"].eq(True)].copy().reset_index(drop=True)
    cohort_frame = derive_trajectory_fields(cohort_frame)
    cohort_frame = derive_followup_and_outcomes(cohort_frame)
    cohort_frame = derive_weights_and_flags(cohort_frame)
    cohort_frame = derive_design_identifiers(cohort_frame)
    cohort_frame = derive_missing_indicators(cohort_frame)
    trajectory_distribution_frame = build_trajectory_distribution(cohort_frame)
    trajectory_distribution_6cat_frame = build_trajectory_distribution_6cat(cohort_frame)
    flow_payload = build_flow_payload(
        linkage_layers,
        cohort_frame,
        trajectory_distribution_frame,
        trajectory_distribution_6cat_frame,
    )
    write_parquet_atomic(cohort_frame, COHORT_OUTPUT_PATH)
    trajectory_distribution_frame.to_csv(TRAJECTORY_DISTRIBUTION_PATH, index=False)
    trajectory_distribution_6cat_frame.to_csv(TRAJECTORY_DISTRIBUTION_6CAT_PATH, index=False)
    with COHORT_FLOW_PATH.open("w", encoding="utf-8") as handle:
        json.dump(flow_payload, handle, indent=2, ensure_ascii=False)
    coverage_md = build_variable_coverage_report(cohort_frame)
    with VARIABLE_COVERAGE_PATH.open("w", encoding="utf-8") as handle:
        handle.write(coverage_md)
    print(f"Wrote {COHORT_OUTPUT_PATH}")
    print(f"Wrote {COHORT_FLOW_PATH}")
    print(f"Wrote {VARIABLE_COVERAGE_PATH}")
    print(f"Wrote {TRAJECTORY_DISTRIBUTION_PATH}")
    print(f"Wrote {TRAJECTORY_DISTRIBUTION_6CAT_PATH}")
    print(f"Cohort n: {flow_payload['gi_cohort_n']}")
    print(f"Stomach n: {flow_payload['stomach_subgroup_n']}")
    print("Trajectory distribution:")
    for trajectory_name, count in flow_payload["trajectory_distribution"].items():
        print(f"  {trajectory_name}: {count}")
    print("Trajectory distribution (6-category):")
    for trajectory_name, count in flow_payload["trajectory_distribution_6cat"].items():
        print(f"  {trajectory_name}: {count}")
    print("Outcome counts:")
    for outcome_name, outcome_counts in flow_payload["outcome_event_counts"].items():
        print(f"  {outcome_name}: 1={outcome_counts['count_1']} 0={outcome_counts['count_0']} NaN={outcome_counts['count_nan']}")
    print("Impossible pattern counts:")
    for flag_name, count in flow_payload["impossible_pattern_counts"].items():
        print(f"  {flag_name}: {count}")


if __name__ == "__main__":
    main()
