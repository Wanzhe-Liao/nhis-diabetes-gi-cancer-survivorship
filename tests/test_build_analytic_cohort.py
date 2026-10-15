import sys
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT_OR_CODE = Path(__file__).resolve().parent.parent
for _candidate in (_ROOT_OR_CODE / "scripts", _ROOT_OR_CODE):
    if (_candidate / "build_analytic_cohort.py").exists():
        sys.path.insert(0, str(_candidate))
        break
from build_analytic_cohort import (
    GI_AGE_COLUMNS,
    derive_design_identifiers,
    derive_followup_and_outcomes,
    derive_trajectory_fields,
)


def test_cause_specific_non_target_deaths_are_zero_and_invalid_followup_is_na():
    frame = pd.DataFrame(
        {
            "year": [2010, 2010, 2017, 2010, 2010],
            "intv_qrt": [1, 1, 4, 1, 1],
            "mortstat": [1, 1, 0, 1, 1],
            "dodyear": [2012, 2012, np.nan, 2009, 2012],
            "dodqtr": [1, 1, np.nan, 4, 9],
            "ucod_leading": [2, 1, np.nan, 2, 2],
            "diabetes": [0, 0, np.nan, 1, 1],
        }
    )

    result = derive_followup_and_outcomes(frame)

    assert result.loc[0, "death_cancer_10y"] == 1.0
    assert result.loc[0, "death_dm_contributing_10y"] == 0.0
    assert result.loc[1, "death_cancer_10y"] == 0.0
    assert result.loc[1, "death_dm_contributing_10y"] == 0.0
    assert pd.isna(result.loc[2, "death_cancer_10y"])
    assert pd.isna(result.loc[3, "death_allcause_10y"])
    assert pd.isna(result.loc[4, "dodqtr"])
    assert pd.isna(result.loc[4, "death_time"])
    assert pd.isna(result.loc[4, "death_cancer_10y"])


def test_design_identifiers_use_period_and_original_fields():
    frame = pd.DataFrame(
        {
            "year": [2010, 2011, 2016],
            "design_strata": [1.0, "02", 3],
            "design_psu": [10.0, "003", 30],
        }
    )

    result = derive_design_identifiers(frame.copy())
    assert result["design_strata_prefixed"].tolist() == [
        "2006-2015.1", "2006-2015.2", "2016-2018.3"
    ]
    assert result["design_psu_prefixed"].tolist() == [
        "2006-2015.1.10", "2006-2015.2.3", "2016-2018.3.30"
    ]


def test_impossible_diagnosis_ages_are_removed_before_trajectory_derivation():
    frame = pd.DataFrame(
        {
            "age": [60, 60, 60],
            "dm_age": [70, 40, 40],
            "dm_ever": [True, True, True],
            **{column: [np.nan, np.nan, np.nan] for column in GI_AGE_COLUMNS},
        }
    )
    frame.loc[0, "canage7"] = 50
    frame.loc[1, "canage7"] = 70
    frame.loc[2, "canage7"] = 45
    frame.loc[2, "canage8"] = 70

    result = derive_trajectory_fields(frame)

    assert pd.isna(result.loc[0, "dm_dx_age"])
    assert bool(result.loc[0, "flag_impossible_dm_age"])
    assert result.loc[0, "trajectory_7cat"] == "dm_order_unknown"
    assert pd.isna(result.loc[1, "gi_first_dx_age"])
    assert bool(result.loc[1, "flag_impossible_gi_age"])
    assert result.loc[1, "trajectory_7cat"] == "dm_order_unknown"
    assert result.loc[2, "gi_first_dx_age"] == 45
    assert bool(result.loc[2, "flag_impossible_gi_age"])
    assert result.loc[2, "trajectory_7cat"] == "dm_to_gi_2_10y"
