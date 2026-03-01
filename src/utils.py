"""
Utility functions -- paths, config, input validation, etc.
"""

import gzip
import os
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml
from rich.console import Console

console = Console()

@dataclass
class Paths:
    root: Path
    config_path: Path
    scfa_csv: Path
    sbml_path: Path
    results: Path
    figs_dir: Path
    tables_dir: Path

def get_root():
    return Path(__file__).resolve().parents[1]

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def build_paths():
    """Build all the paths we need. Creates output dirs if missing."""
    root = get_root()
    cfg_path = root / "data" / "inputs" / "project_config.yml"
    cfg = load_config(cfg_path)

    scfa_csv = root / "data" / "inputs" / "scfa_inputs.csv"
    sbml = root / cfg["human_model"]["sbml_path"]
    results = root / "results"
    figs = root / "outputs" / "figs"
    tables = root / "outputs" / "tables"

    for d in [results, figs, tables]:
        d.mkdir(parents=True, exist_ok=True)

    return Paths(
        root=root, config_path=cfg_path, scfa_csv=scfa_csv,
        sbml_path=sbml, results=results, figs_dir=figs, tables_dir=tables
    )


def read_scfa_inputs(path, expected_conditions):
    """Read and validate the SCFA input CSV."""
    df = pd.read_csv(path)

    # check columns
    needed = {"condition", "acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr", "butyrate_mmol_gDW_hr"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # check conditions match
    df["condition"] = df["condition"].astype(str)
    if set(df["condition"]) != set(expected_conditions):
        raise ValueError(
            f"Conditions don't match.\nExpected: {sorted(expected_conditions)}\n"
            f"Got: {sorted(df['condition'].unique())}"
        )

    # no negatives
    for c in ["acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr", "butyrate_mmol_gDW_hr"]:
        if (df[c] < 0).any():
            raise ValueError(f"Negative values in {c}")

    return df.sort_values("condition").reset_index(drop=True)


def decompress_gz(sbml_path):
    """Decompress .xml.gz to cache/ if needed. Returns path to .xml file."""
    if sbml_path.suffix != ".gz":
        return sbml_path

    cache = sbml_path.parent / "cache"
    cache.mkdir(exist_ok=True)
    out = cache / sbml_path.name.replace(".gz", "")

    if out.exists() and out.stat().st_size > 0:
        return out

    console.log(f"Decompressing {sbml_path.name}...")
    with gzip.open(sbml_path, "rb") as fin:
        with open(out, "wb") as fout:
            fout.write(fin.read())
    return out
