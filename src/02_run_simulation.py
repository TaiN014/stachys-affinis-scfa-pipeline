#!/usr/bin/env python3
"""
02_run_simulation.py

Load Recon3D, set up hepatocyte-like medium, run ATPM maximization
for each SCFA dose condition.
"""

import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model

from .utils import build_paths, load_config, decompress_gz, read_scfa_inputs


# Recon3D uses inconsistent naming for exchange rxns
SCFA_EXCHANGE_IDS = {
    "acetate":    ["EX_ac_e", "EX_ac[u]", "EX_ac(e)"],
    "propionate": ["EX_ppa_e", "EX_propn_e", "EX_ppn_e"],
    "butyrate":   ["EX_but_e", "EX_btn_e"],
}
GLUCOSE_IDS = ["EX_glc__D_e", "EX_glc_D_e", "EX_glc[e]"]
O2_IDS      = ["EX_o2_e", "EX_o2[e]"]
ATPM_IDS    = ["ATPM", "DM_atp_c_"]


# medium components
FREE_EXCHANGE = ["EX_h2o_e", "EX_h_e"]

# inorganic ions
IONS = {
    "EX_pi_e": 10, "EX_so4_e": 10, "EX_k_e": 10, "EX_na1_e": 10,
    "EX_ca2_e": 10, "EX_cl_e": 10, "EX_mg2_e": 10, "EX_fe2_e": 10,
}

# secretion
SECRETION_RXNS = {"EX_co2_e": 1000, "EX_nh4_e": 100}

# essential AAs
ESSENTIAL_AA = [
    "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_thr__L_e", "EX_trp__L_e",
    "EX_val__L_e",
]

VITAMINS = [
    "EX_thm_e", "EX_ribflv_e", "EX_ncam_e", "EX_pnto__R_e",
    "EX_pydxn_e", "EX_fol_e", "EX_cbl1_e", "EX_chol_e", "EX_inost_e",
]

# pathway reactions to read out
PATHWAY_RXNS = [
    ("PYK",         "pyruvate kinase"),
    ("PDHm",        "pyruvate dehydrogenase"),
    ("CSm",         "citrate synthase"),
    ("ACONTm",      "aconitase"),
    ("ICDHxm",      "isocitrate DH"),
    ("AKGDm",       "alpha-KG DH"),
    ("SUCDi",       "succinate DH"),
    ("FUMm",        "fumarase"),
    ("MDHm",        "malate DH"),
    ("PCm",         "pyruvate carboxylase"),
    ("PEPCK",       "PEPCK"),
    ("G6PDH2r",     "G6PD"),
    ("FBA",         "aldolase"),
    ("PFK",         "PFK"),
    ("ATPS4mi",     "ATP synthase"),
    ("NADH2_u10mi", "complex I"),
    ("CYOOm3i",     "complex IV"),
]


def _find_rxn(model, id_list):
    """Return first matching reaction or None."""
    for rid in id_list:
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    return None


def setup_medium(model, cfg):
    """
    Hepatocyte-like medium: cap internals, close all boundary rxns,
    reopen a curated set. Prevents thermodynamic loops.
    """
    sim_cfg = cfg["host_simulation"]

    # cap internal reactions
    n_capped = 0
    for rxn in model.reactions:
        if rxn.id.startswith(("EX_", "DM_", "sink_", "SK_")):
            continue
        if rxn.lower_bound < -500:
            rxn.lower_bound = -500
            n_capped += 1
        if rxn.upper_bound > 500:
            rxn.upper_bound = 500
            n_capped += 1

    # close all boundary rxns
    n_closed = 0
    for rxn in model.reactions:
        if rxn.id.startswith(("EX_", "DM_", "sink_", "SK_")):
            rxn.bounds = (0, 0)
            n_closed += 1
    print(f"  capped {n_capped} internal bounds, closed {n_closed} boundary rxns")

    # reopen essentials
    for rid in FREE_EXCHANGE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-1000, 1000)

    # ions
    for rid, ub in IONS.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-ub, 100)

    # secretion
    for rid, ub in SECRETION_RXNS.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (0, ub)
    # small NH4 uptake
    if "EX_nh4_e" in model.reactions:
        model.reactions.get_by_id("EX_nh4_e").lower_bound = -0.5

    # O2
    o2_uptake = float(sim_cfg["oxygen_uptake"])
    o2_rxn = _find_rxn(model, O2_IDS)
    if o2_rxn:
        o2_rxn.bounds = (-o2_uptake, 0)

    # glucose (scarce on purpose)
    glc_uptake = float(sim_cfg["glucose_uptake"])
    glc_rxn = _find_rxn(model, GLUCOSE_IDS)
    if glc_rxn:
        glc_rxn.bounds = (-glc_uptake, 0)

    # AAs
    aa_rate = float(sim_cfg.get("amino_acid_uptake", 0.01))
    for rid in ESSENTIAL_AA:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-aa_rate, 0)

    # vitamins
    vit_rate = float(sim_cfg.get("vitamin_uptake", 0.01))
    for rid in VITAMINS:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-vit_rate, 0)


def collect_pathway_fluxes(model, solution):
    result = {}
    for rid, _label in PATHWAY_RXNS:
        if rid in model.reactions:
            val = solution.fluxes.get(rid, float("nan"))
            result[f"pathway_{rid}"] = float(val)
        else:
            result[f"pathway_{rid}"] = float("nan")
    return result


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]


    scfa_df = read_scfa_inputs(
        paths.results / "scfa_inputs_canonical.csv", conditions
    )

    # load model
    sbml_file = decompress_gz(paths.sbml_path)
    print(f"loading model from {sbml_file}...")
    model = read_sbml_model(str(sbml_file))
    print(f"  loaded: {len(model.reactions)} rxns, {len(model.metabolites)} mets")

    # set up medium
    setup_medium(model, cfg)

    # objective = ATPM
    atpm_rxn = _find_rxn(model, ATPM_IDS)
    if atpm_rxn is None:
        raise RuntimeError("can't find ATPM reaction")
    atpm_rxn.bounds = (0, 500)
    model.objective = atpm_rxn
    print(f"objective: {atpm_rxn.id}")

    # SCFA exchange rxns
    scfa_rxns = {}
    for name, candidates in SCFA_EXCHANGE_IDS.items():
        scfa_rxns[name] = _find_rxn(model, candidates)
        found_id = scfa_rxns[name].id if scfa_rxns[name] else "NOT FOUND"
        print(f"  {name}: {found_id}")

    glc_rxn = _find_rxn(model, GLUCOSE_IDS)
    o2_rxn = _find_rxn(model, O2_IDS)

    # baseline
    print("\nrunning baseline (no SCFAs)...")
    baseline = model.optimize()
    baseline_atpm = float(baseline.objective_value or 0)
    print(f"  baseline ATPM = {baseline_atpm:.4f}")

    # dose conditions
    rows = []
    for _, row in scfa_df.iterrows():
        cond = row["condition"]
        ac_dose  = float(row["acetate_mmol_gDW_hr"])
        ppa_dose = float(row["propionate_mmol_gDW_hr"])
        but_dose = float(row["butyrate_mmol_gDW_hr"])

        print(f"\nrunning {cond}...")

        # context manager resets bounds after
        with model:
            for scfa_name, dose in [("acetate", ac_dose),
                                     ("propionate", ppa_dose),
                                     ("butyrate", but_dose)]:
                rxn = scfa_rxns.get(scfa_name)
                if rxn and dose > 0:
                    rxn.bounds = (-dose, 0)

            sol = model.optimize()
            atpm_val = float(sol.objective_value or 0)
            delta = atpm_val - baseline_atpm
            pct_change = 100.0 * delta / baseline_atpm if baseline_atpm else float("nan")

            def get_flux(r):
                return float(sol.fluxes.get(r.id, 0)) if r else None

            result = {
                "condition": cond,
                "objective_id": atpm_rxn.id,
                "objective_value": atpm_val,
                "baseline_objective": baseline_atpm,
                "objective_delta": delta,
                "objective_pct_change": pct_change,
                "glucose_flux": get_flux(glc_rxn),
                "oxygen_flux": get_flux(o2_rxn),
                "co2_flux": float(sol.fluxes.get("EX_co2_e", 0)),
                "acetate_flux": get_flux(scfa_rxns.get("acetate")),
                "propionate_flux": get_flux(scfa_rxns.get("propionate")),
                "butyrate_flux": get_flux(scfa_rxns.get("butyrate")),
            }
            result.update(collect_pathway_fluxes(model, sol))
            rows.append(result)

            print(f"  ATPM = {atpm_val:.2f}  ({pct_change:+.1f}% vs baseline)")

    # save
    host_df = pd.DataFrame(rows)
    host_df.to_csv(paths.results / "host_fluxes_by_condition.csv", index=False)

    merged = scfa_df.merge(host_df, on="condition", how="left")
    merged.to_csv(paths.results / "merged_dose_scfa_host.csv", index=False)

    print("\n--- done ---")
    for _, r in host_df.iterrows():
        print(f"  {r['condition']:25s} ATPM={r['objective_value']:.2f}  "
              f"delta={r['objective_pct_change']:+.1f}%")


if __name__ == "__main__":
    main()
