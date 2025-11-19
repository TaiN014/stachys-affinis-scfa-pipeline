# utils.py -- shared helpers

import gzip
import os
from pathlib import Path

import pandas as pd
import yaml

class Paths:
    """Holds all the directory paths the pipeline needs."""
    def __init__(self, root, config_path, scfa_csv, sbml_path, results, figs_dir, tables_dir):
        self.root = root
        self.config_path = config_path
        self.scfa_csv = scfa_csv
        self.sbml_path = sbml_path
        self.results = results
        self.figs_dir = figs_dir
        self.tables_dir = tables_dir


def get_root():
    return Path(__file__).resolve().parents[1]


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def build_paths():
    root = get_root()
    cfg_path = root / "data" / "inputs" / "project_config.yml"
    cfg = load_config(cfg_path)

    scfa_csv = root / "data" / "inputs" / "scfa_inputs.csv"
    sbml = root / cfg["human_model"]["sbml_path"]
    results_dir = root / "results"
    figs = root / "outputs" / "figs"
    tables = root / "outputs" / "tables"

    for d in [results_dir, figs, tables]:
        d.mkdir(parents=True, exist_ok=True)

    return Paths(root, cfg_path, scfa_csv, sbml, results_dir, figs, tables)


def read_scfa_inputs(path, expected_conditions):
    """Load and validate SCFA input csv."""
    df = pd.read_csv(path)

    required_cols = {"condition", "acetate_mmol_gDW_hr",
                     "propionate_mmol_gDW_hr", "butyrate_mmol_gDW_hr"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"SCFA input is missing columns: {missing}")

    # conditions should match config
    df["condition"] = df["condition"].astype(str)
    got = set(df["condition"])
    want = set(expected_conditions)
    if got != want:
        raise ValueError(
            f"Condition mismatch!\n  config says: {sorted(want)}\n  csv has: {sorted(got)}"
        )

    # no negatives
    for col in ["acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr", "butyrate_mmol_gDW_hr"]:
        if (df[col] < 0).any():
            raise ValueError(f"Negative values in {col}")

    # sort high->low
    order = {"StachysDose_High": 0, "StachysDose_Mid": 1, "StachysDose_Low": 2}
    df["_sort"] = df["condition"].map(order)
    df = df.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    return df


def decompress_gz(sbml_path):
    """Decompress .xml.gz to cache/ if needed."""
    if sbml_path.suffix != ".gz":
        return sbml_path

    cache_dir = sbml_path.parent / "cache"
    cache_dir.mkdir(exist_ok=True)
    out_path = cache_dir / sbml_path.name.replace(".gz", "")

    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    print(f"  decompressing {sbml_path.name}...")
    with gzip.open(sbml_path, "rb") as fin:
        with open(out_path, "wb") as fout:
            fout.write(fin.read())
    return out_path
