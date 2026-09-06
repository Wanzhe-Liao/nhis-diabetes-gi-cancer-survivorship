import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from harmonize_nhis_gi_diabetes import (
    clean_activity_frequency,
    derive_alcohol_status,
    derive_poverty_3cat,
    derive_poverty_4cat,
    derive_race_ethnicity_4cat,
    insurance_plan_recode_boolean,
    military_recode_boolean,
    normalize_key_series,
    other_public_recode_boolean,
)
from project_paths import ANALYSIS_YEARS, OPTIONAL_SENSITIVITY_YEARS


def test_primary_year_set_includes_all_22_years():
    assert ANALYSIS_YEARS == list(range(1997, 2019))
    assert 2007 in ANALYSIS_YEARS
    assert OPTIONAL_SENSITIVITY_YEARS == []


def test_numeric_merge_keys_are_canonicalized_without_destroying_text():
    result = normalize_key_series(pd.Series(["000123", "01", "2.0", " A1 ", ""]))
    assert result.iloc[:4].tolist() == ["123", "1", "2", "A1"]
    assert pd.isna(result.iloc[4])


def test_activity_frequency_distinguishes_zero_from_nonresponse():
    result = clean_activity_frequency(pd.Series([0, 3, 28, 95, 96, 97, 98, 99, 29, -1]))
    assert result.iloc[:5].tolist() == [0, 3, 28, 0, 0]
    assert result.iloc[5:].isna().all()


def test_alcohol_status_uses_year_specific_codebooks():
    early = derive_alcohol_status(pd.DataFrame({"alcohol_status": [1, 2, 3, 9]}), 2003)
    late = derive_alcohol_status(pd.DataFrame({"alcohol_status": list(range(1, 11))}), 2004)
    assert early["alcohol_status"].tolist()[:3] == ["never", "former", "current"]
    assert pd.isna(early.loc[3, "alcohol_status"])
    assert late["alcohol_status"].tolist()[:9] == [
        "never", "former", "former", "former",
        "current", "current", "current", "current", "current",
    ]
    assert pd.isna(late.loc[9, "alcohol_status"])


def test_poverty_category_codes_are_not_treated_as_continuous_ratios():
    frame = pd.DataFrame({"poverty": [1, 3, 4, 7, 8, 11, 12, 14, 15, 16, 17, 18, 99]})
    result = derive_poverty_3cat(derive_poverty_4cat(frame))
    assert result["poverty_3cat"].tolist()[:12] == [
        "lt_2_0", "lt_2_0", "lt_2_0", "lt_2_0",
        "2_0_to_3_99", "2_0_to_3_99", "ge_4_0", "ge_4_0",
        "lt_2_0", "lt_2_0", "2_0_to_3_99", "ge_4_0",
    ]
    assert pd.isna(result.loc[12, "poverty_3cat"])


def test_race_ethnicity_groups_are_mutually_exclusive():
    frame = pd.DataFrame({
        "race": [1, 1, 2, 3, 97, 2],
        "hispanic": [1, 12, 12, 12, 12, 99],
    })
    result = derive_race_ethnicity_4cat(frame)
    assert result["race_ethnicity_4cat"].tolist()[:4] == [
        "hispanic", "non_hispanic_white", "non_hispanic_black",
        "non_hispanic_other",
    ]
    assert pd.isna(result.loc[4, "race_ethnicity_4cat"])
    assert pd.isna(result.loc[5, "race_ethnicity_4cat"])


def test_insurance_recode_families_keep_unknown_codes_missing():
    plan = insurance_plan_recode_boolean(pd.Series([1, 2, 3, 7, 8, 9]))
    other = other_public_recode_boolean(pd.Series([1, 2, 3, 7, 8, 9]))
    military_old = military_recode_boolean(pd.Series([1, 2, 3, 4, 5, 7]), 2003)
    military_new = military_recode_boolean(pd.Series([1, 2, 3, 7]), 2004)
    assert plan.iloc[:3].tolist() == [True, True, False]
    assert plan.iloc[3:].isna().all()
    assert other.iloc[:2].tolist() == [True, False]
    assert other.iloc[2:].isna().all()
    assert military_old.iloc[:5].tolist() == [True, True, True, True, False]
    assert pd.isna(military_old.iloc[5])
    assert military_new.iloc[:3].tolist() == [True, True, False]
    assert pd.isna(military_new.iloc[3])
