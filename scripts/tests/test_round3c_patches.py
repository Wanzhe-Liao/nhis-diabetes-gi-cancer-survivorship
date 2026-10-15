"""Unit tests for Round 3c patches (B1-B4)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from harmonize_nhis_gi_diabetes import (
    derive_education_5cat,
    derive_education_4cat,
    derive_insurance_type,
    derive_poverty_3cat,
    merge_supplemental_cholesterol,
)


class TestCollapsedCategories:
    def test_derive_education_4cat_preserves_official_four_groups(self):
        frame = pd.DataFrame({
            "education_5cat": ["lt_high_school", "high_school", "some_college", "college_grad", pd.NA],
        })
        frame = derive_education_4cat(frame)
        assert frame["education_4cat"].iloc[:4].tolist() == [
            "lt_high_school", "high_school", "some_college", "college_grad"
        ]
        assert pd.isna(frame.loc[4, "education_4cat"])

    def test_education_credentials_follow_official_00_to_21_recode(self):
        frame = pd.DataFrame({"education": [0, 12, 13, 14, 15, 17, 18, 21, 97]})
        result = derive_education_5cat(frame)
        assert result["education_5cat"].tolist()[:8] == [
            "lt_high_school", "lt_high_school", "high_school", "high_school",
            "some_college", "some_college", "college_grad", "college_grad",
        ]
        assert pd.isna(result.loc[8, "education_5cat"])

    def test_derive_poverty_3cat_collapses_lt_2_0(self):
        frame = pd.DataFrame({
            "poverty_4cat": ["lt_1_0", "1_0_to_1_99", "2_0_to_3_99", "ge_4_0"],
        })
        frame = derive_poverty_3cat(frame)
        assert frame["poverty_3cat"].tolist() == ["lt_2_0", "lt_2_0", "2_0_to_3_99", "ge_4_0"]


class TestInsuranceDerivation:
    def test_insurance_cover65_medicare_private(self):
        frame = pd.DataFrame({
            "insurance_notcov": pd.Series([2, 2, 2, 1], dtype="Int64"),
            "insurance_cover": pd.Series([pd.NA, pd.NA, pd.NA, pd.NA], dtype="Int64"),
            "insurance_cover65": pd.Series([1, 2, 5, pd.NA], dtype="Int64"),
            "age": pd.Series([70, 70, 70, 70], dtype="float"),
            "year": [2015, 2015, 2015, 2015],
        })
        frame = derive_insurance_type(frame)
        assert frame["insurance_type"].tolist() == ["private_only", "medicare_dual", "other", "uninsured"]

    def test_insurance_cover_private_only_under65(self):
        frame = pd.DataFrame({
            "insurance_notcov": pd.Series([2, 2, 2, 1], dtype="Int64"),
            "insurance_cover": pd.Series([1, 2, 3, pd.NA], dtype="Int64"),
            "insurance_cover65": pd.Series([pd.NA, pd.NA, pd.NA, pd.NA], dtype="Int64"),
            "age": pd.Series([50, 50, 50, 50], dtype="float"),
            "year": [2015, 2015, 2015, 2015],
        })
        frame = derive_insurance_type(frame)
        assert frame["insurance_type"].tolist() == ["private_only", "medicaid_only", "other", "uninsured"]

    def test_insurance_medicare_medicaid_priority(self):
        frame = pd.DataFrame({
            "insurance_notcov": pd.Series([2, 2, 2, 2], dtype="Int64"),
            "medicare_flag": pd.Series([True, True, pd.NA, pd.NA], dtype="boolean"),
            "medicaid_flag": pd.Series([pd.NA, True, True, pd.NA], dtype="boolean"),
            "private_flag": pd.Series([pd.NA, pd.NA, pd.NA, True], dtype="boolean"),
            "age": pd.Series([70, 70, 70, 70], dtype="float"),
            "year": [2014, 2014, 2014, 2014],
        })
        frame = derive_insurance_type(frame)
        assert frame["insurance_type"].tolist() == [
            "medicare_only", "medicare_dual", "medicaid_only", "private_only"
        ]


class TestSupplementalCholesterolMerge:
    def test_merge_supplemental_drops_placeholder_and_merges(self):
        frame = pd.DataFrame({
            "hhx": ["1", "2", "3"],
            "fmx": ["01", "01", "01"],
            "person_key": ["01", "02", "03"],
            "cholesterol_high_ever": pd.Series([pd.NA, pd.NA, pd.NA], dtype="Int64"),
        })
        supp = pd.DataFrame({
            "hhx": ["1", "2"],
            "fmx": ["01", "01"],
            "person_key": ["01", "02"],
            "cholesterol_high_ever": pd.Series([1, 2], dtype="Int64"),
        })
        # Manually inject supplemental data by temporarily overriding the reader
        import harmonize_nhis_gi_diabetes as hmod
        original_reader = hmod.read_supplemental_cholesterol
        hmod.read_supplemental_cholesterol = lambda year: supp
        try:
            result = merge_supplemental_cholesterol(frame, 1998)
            assert result.loc[0, "cholesterol_high_ever"] == 1
            assert result.loc[1, "cholesterol_high_ever"] == 2
            assert pd.isna(result.loc[2, "cholesterol_high_ever"])
        finally:
            hmod.read_supplemental_cholesterol = original_reader


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
