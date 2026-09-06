# build_domain_fullsample.py —— 为 domain/subpopulation 正确方差估计构建全样本骨架。
# 输出: outputs/cohort/analytic_cohort_v4inc2007_domain_fullsample.parquet
# 内容: 所有 ELIGSTAT=1 的 pooled Sample Adults (646,201 行)，含设计变量与 gi_any 标记。
# 非 GI 行不携带分析协变量——它们只向 survey design 提供 PSU/strata 设计信息。
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import ANALYSIS_YEARS, OUTPUTS_DIR
from build_analytic_cohort import (
    coerce_input_types,
    derive_design_identifiers,
    load_input_frames,
    merge_nhis_lmf,
    write_parquet_atomic,
)

OUT_PATH = OUTPUTS_DIR / "cohort" / "analytic_cohort_v4inc2007_domain_fullsample.parquet"

KEEP_COLUMNS = [
    "publicid",
    "year",
    "design_strata",
    "design_psu",
    "design_strata_prefixed",
    "design_psu_prefixed",
    "sa_wgt_new",
    "sa_wgt_pool",
    "gi_any",
    "mortstat",
]


def main() -> None:
    nhis_frame, lmf_frame = load_input_frames()
    nhis_frame, lmf_frame = coerce_input_types(nhis_frame, lmf_frame)
    merged_frame, layers = merge_nhis_lmf(nhis_frame, lmf_frame)
    pooled_years = len(ANALYSIS_YEARS)
    merged_frame["sa_wgt_pool"] = merged_frame["sa_wgt_new"] / pooled_years
    merged_frame = derive_design_identifiers(merged_frame)
    skeleton = merged_frame[KEEP_COLUMNS].copy()
    write_parquet_atomic(skeleton, OUT_PATH)
    print(f"Wrote {OUT_PATH}")
    print(f"rows: {len(skeleton)} (expected eligstat_1_n = {layers['eligstat_1_n']})")
    print(f"gi_any TRUE: {int(skeleton['gi_any'].eq(True).sum())}")
    print(f"design strata: {skeleton['design_strata_prefixed'].nunique()}, "
          f"psu: {skeleton['design_psu_prefixed'].nunique()}")
    missing_design = skeleton["design_psu_prefixed"].isna().sum()
    print(f"rows missing design ids: {missing_design}")


if __name__ == "__main__":
    main()
