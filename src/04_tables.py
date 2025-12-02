#!/usr/bin/env python3
"""
04_tables.py -- csv tables for the manuscript.
"""

import pandas as pd
from .utils import build_paths

DOSE_ORDER = {
    "StachysDose_High": 0,
    "StachysDose_Mid": 1,
    "StachysDose_Low": 2,
}


def main():
    paths = build_paths()
    df = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # sort
    df["_s"] = df["condition"].map(DOSE_ORDER)
    df = df.sort_values("_s").drop(columns="_s").reset_index(drop=True)

    # table 1: SCFA inputs
    t1 = df[["condition", "acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr",
             "butyrate_mmol_gDW_hr"]].copy()
    t1.columns = ["Condition", "Acetate (mmol/gDW/hr)",
                   "Propionate (mmol/gDW/hr)", "Butyrate (mmol/gDW/hr)"]
    t1.to_csv(paths.tables_dir / "table_scfa_inputs.csv", index=False)
    print("  table_scfa_inputs.csv")

    # table 2: host fluxes
    want_cols = [
        "condition", "objective_id", "objective_value", "baseline_objective",
        "objective_delta", "objective_pct_change",
        "glucose_flux", "oxygen_flux", "co2_flux",
        "acetate_flux", "propionate_flux", "butyrate_flux",
    ]

    want_cols = [c for c in want_cols if c in df.columns]
    t2 = df[want_cols]
    t2.to_csv(paths.tables_dir / "table_host_fluxes.csv", index=False)
    print("  table_host_fluxes.csv")

    # table 3: summary
    summary_cols = ["condition", "objective_value"]
    for scfa in ["acetate", "propionate", "butyrate"]:
        for suf in ["_mmol_gDW_hr", "_flux"]:
            c = scfa + suf
            if c in df.columns:
                summary_cols.append(c)

    t3 = df[[c for c in summary_cols if c in df.columns]]
    t3.to_csv(paths.tables_dir / "table_summary.csv", index=False)
    print("  table_summary.csv")

    print("done")


if __name__ == "__main__":
    main()
