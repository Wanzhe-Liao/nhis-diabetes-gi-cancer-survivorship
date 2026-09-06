"""Read the complete 1997-2018 primary analysis from domain-correct CSV results."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "outputs/revision_round1_v4_r461_sens2007"
REFERENCE_21 = ROOT / "outputs/revision_round1_v4_r461_submission"
SENS2007 = SOURCE
COHORT_PATH = ROOT / "outputs/cohort/analytic_cohort_v4inc2007.parquet"
FULLSAMPLE_PATH = ROOT / "outputs/cohort/analytic_cohort_v4inc2007_domain_fullsample.parquet"
POVERTY_PATH = SOURCE / "income_repair/poverty_repaired.csv"
FLOW = json.loads((ROOT / "outputs/cohort/cohort_flow_v4inc2007.json").read_text(encoding="utf-8"))


def load_results():
    models = pd.read_csv(SOURCE / "master_results.csv")
    principal = models.loc[models.model.eq("A_principal_5cat")].iloc[0]
    result = {"meta": {"status": "VALIDATED_V4_R461_SUBMISSION", "cohort_n": FLOW["gi_cohort_n"],
                       "n_analysis": int(principal["n"]), "events_10y": int(principal["events"]),
                       "full_population_n": FLOW["linkage_layers"]["eligstat_1_n"],
                       "sample_adult_n": FLOW["linkage_layers"]["sample_adult_n"],
                       "analysis_years": list(range(1997, 2019))}}
    for key, model in {
        "A_principal_5cat": "A_principal_5cat",
        "A2_no_income_sensitivity": "A2_no_income_sensitivity",
        "B_lag_6cat": "B_6cat_repaired_income",
        "C_plus_cancer_to_interview": "C_plus_cancer_to_interview",
        "D_t1d_proxy": "D_t1d_proxy_exclusion",
        "G_zero_time_time_z05": "G_zero_time_time_z05",
        "G_zero_time_time_z1": "G_zero_time_time_z1",
        "G_zero_time_time_z3": "G_zero_time_time_z3",
    }.items():
        result[key] = models.loc[models.model.eq(model)].to_dict("records")
        assert result[key], model
    result["E_nested"] = models.loc[models.model.str.startswith("E_nested_")].to_dict("records")
    for key, filename in {
        "B_wald": "wald_lag_test.csv",
        "C_timing_distributions": "timing_distributions.csv",
        "F_burden_comparison": "burden_comparison.csv",
        "G_absolute_risk_bootstrap": "rev1_absolute_risk_with_rd.csv",
        "H_ph_diagnostics": "ph_diagnostics_rev1.csv",
    }.items():
        result[key] = pd.read_csv(SOURCE / filename).to_dict("records")
    assert all(r["n_boot_converged"] == 500 for r in result["G_absolute_risk_bootstrap"])
    return result


def hr_ci(row, digits=2, dash="–"):
    return f"{row['HR']:.{digits}f} ({row['ci_lo']:.{digits}f}{dash}{row['ci_hi']:.{digits}f})"


def p_value(value):
    return f"{value:.2e}" if value < 0.001 else f"{value:.3f}"
