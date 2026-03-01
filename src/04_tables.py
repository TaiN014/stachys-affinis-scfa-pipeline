"""
Step 04 - make tables for the paper.
"""

import pandas as pd
from rich.console import Console
from .utils import build_paths

console = Console()

def main():
    paths = build_paths()
    merged = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # Table 1: SCFA inputs
    t1 = merged[["condition", "acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr",
                  "butyrate_mmol_gDW_hr"]].copy()
    t1.columns = ["Condition", "Acetate (mmol/gDW/hr)", "Propionate (mmol/gDW/hr)",
                   "Butyrate (mmol/gDW/hr)"]
    t1.to_csv(paths.tables_dir / "table_scfa_inputs.csv", index=False)
    console.log("table_scfa_inputs.csv")

    # Table 2: host fluxes
    cols = ["condition", "objective_id", "objective_value", "baseline_objective",
            "objective_delta", "objective_pct_change",
            "glucose_flux", "oxygen_flux", "co2_flux",
            "acetate_flux", "propionate_flux", "butyrate_flux"]
    cols = [c for c in cols if c in merged.columns]
    t2 = merged[cols]
    t2.to_csv(paths.tables_dir / "table_host_fluxes.csv", index=False)
    console.log("table_host_fluxes.csv")

    # Table 3: simple summary
    scols = ["condition", "objective_value"]
    for s in ["acetate", "propionate", "butyrate"]:
        for suffix in ["_mmol_gDW_hr", "_flux"]:
            c = s + suffix
            if c in merged.columns:
                scols.append(c)
    t3 = merged[[c for c in scols if c in merged.columns]]
    t3.to_csv(paths.tables_dir / "table_summary.csv", index=False)
    console.log("table_summary.csv")

    console.log("Done.")

if __name__ == "__main__":
    main()
