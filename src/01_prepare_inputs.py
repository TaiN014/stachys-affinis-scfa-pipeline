"""
Step 01 - validate SCFA inputs and write a canonical copy.
"""

import pandas as pd
from rich.console import Console
from .utils import build_paths, load_config, read_scfa_inputs

console = Console()

def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]

    console.log("Step 01: validating SCFA inputs")
    df = read_scfa_inputs(paths.scfa_csv, conditions)

    out = paths.results / "scfa_inputs_canonical.csv"
    df.to_csv(out, index=False)
    console.log(f"Wrote {out}")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
