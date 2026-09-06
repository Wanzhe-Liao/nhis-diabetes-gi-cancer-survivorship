import sys
from pathlib import Path

import pandas as pd

_ROOT_OR_CODE = Path(__file__).resolve().parent.parent
for _candidate in (_ROOT_OR_CODE / "scripts", _ROOT_OR_CODE):
    if (_candidate / "parse_nhis_lmf.py").exists():
        sys.path.insert(0, str(_candidate))
        break
from parse_nhis_lmf import build_year_audit, sanitize_dodqtr_series


def test_sanitize_dodqtr_series_masks_reserved_and_invalid_quarters():
    result = sanitize_dodqtr_series(pd.Series([1, 2, 3, 4, 0, 5, 9, "", ".", "2"]))

    assert result.iloc[:4].tolist() == [1, 2, 3, 4]
    assert result.iloc[9] == 2
    assert result.iloc[4:9].isna().all()


def test_year_audit_counts_masked_invalid_dodqtr_values():
    frame = pd.DataFrame(
        {
            "eligstat": [1, 1, 2, 1],
            "mortstat": [1, 0, 1, 1],
            "ucod_leading": [2, 1, 2, 3],
            "diabetes": [0, 1, 1, 0],
            "dodqtr_raw": [1, 9, 0, 5],
        }
    )

    audit = build_year_audit(frame)

    assert audit["eligible_rows"] == 3
    assert audit["deaths"] == 2
    assert audit["cancer_deaths"] == 1
    assert audit["invalid_or_reserved_dodqtr_masked"] == 3
