"""
Step 02 - Run FBA with COBRApy on Recon3D.

Loads model, sets up a strict hepatocyte medium (close all boundary rxns,
reopen curated set), then runs ATPM maximization for each SCFA dose condition.
"""

import pandas as pd
import numpy as np
from rich.console import Console
import cobra
from cobra.io import read_sbml_model

from .utils import build_paths, load_config, decompress_gz, read_scfa_inputs

console = Console()

# --- reaction ID candidates (Recon3D naming is inconsistent) ---
SCFA_RXN_IDS = {
    "acetate":    ["EX_ac_e", "EX_ac[u]", "EX_ac(e)"],
    "propionate": ["EX_ppa_e", "EX_propn_e", "EX_ppn_e"],
    "butyrate":   ["EX_but_e", "EX_btn_e"],
}
GLUCOSE_IDS = ["EX_glc__D_e", "EX_glc_D_e", "EX_glc[e]"]
O2_IDS = ["EX_o2_e", "EX_o2[e]"]
ATPM_IDS = ["ATPM", "DM_atp_c_"]

# --- medium components ---
# water and protons can flow freely
FREE = ["EX_h2o_e", "EX_h_e"]

# inorganic ions - moderate uptake allowed
IONS = {
    "EX_pi_e": 10, "EX_so4_e": 10, "EX_k_e": 10, "EX_na1_e": 10,
    "EX_ca2_e": 10, "EX_cl_e": 10, "EX_mg2_e": 10, "EX_fe2_e": 10,
}

# secretion products
SECRETION = {"EX_co2_e": 1000, "EX_nh4_e": 100}

ESSENTIAL_AA = [
    "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_thr__L_e", "EX_trp__L_e", "EX_val__L_e",
]

VITAMINS = [
    "EX_thm_e", "EX_ribflv_e", "EX_ncam_e", "EX_pnto__R_e",
    "EX_pydxn_e", "EX_fol_e", "EX_cbl1_e", "EX_chol_e", "EX_inost_e",
]

# pathway reactions we want to read out
PATHWAYS = [
    ("PYK", "pyruvate kinase"),
    ("PDHm", "pyruvate dehydrogenase"),
    ("CSm", "citrate synthase"),
    ("ACONTm", "aconitase"),
    ("ICDHxm", "isocitrate DH"),
    ("AKGDm", "alpha-KG DH"),
    ("SUCDi", "succinate DH"),
    ("FUMm", "fumarase"),
    ("MDHm", "malate DH"),
    ("PCm", "pyruvate carboxylase"),
    ("PEPCK", "PEPCK"),
    ("G6PDH2r", "G6PD"),
    ("FBA", "aldolase"),
    ("PFK", "PFK"),
    ("ATPS4mi", "ATP synthase"),
    ("NADH2_u10mi", "complex I"),
    ("CYOOm3i", "complex IV"),
]


def find_rxn(model, candidates):
    """Try a list of possible IDs, return the first one that exists."""
    for rid in candidates:
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    return None


def setup_medium(model, cfg):
    """
    Set up hepatocyte-like medium. The key idea is to close EVERYTHING
    first, then selectively reopen what we want. This prevents the
    thermodynamic loops that plague Recon3D.
    """
    sim = cfg["host_simulation"]

    # cap internal rxns
    capped = 0
    for r in model.reactions:
        if not r.id.startswith(("EX_", "DM_", "sink_", "SK_")):
            if r.lower_bound < -500:
                r.lower_bound = -500
                capped += 1
            if r.upper_bound > 500:
                r.upper_bound = 500
                capped += 1

    # close ALL boundary rxns
    closed = 0
    for r in model.reactions:
        if r.id.startswith(("EX_", "DM_", "sink_", "SK_")):
            r.bounds = (0, 0)
            closed += 1

    console.log(f"  capped {capped} internal rxns, closed {closed} boundary rxns")

    # reopen water/protons
    for rid in FREE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-1000, 1000)

    # ions
    for rid, ub in IONS.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-ub, 100)

    # secretion
    for rid, ub in SECRETION.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (0, ub)
    # also allow small NH4 uptake
    if "EX_nh4_e" in model.reactions:
        model.reactions.get_by_id("EX_nh4_e").lower_bound = -0.5

    # O2
    o2_bound = float(sim["oxygen_uptake"])
    rxn = find_rxn(model, O2_IDS)
    if rxn:
        rxn.bounds = (-o2_bound, 0)

    # glucose (scarce on purpose)
    glc = float(sim["glucose_uptake"])
    rxn = find_rxn(model, GLUCOSE_IDS)
    if rxn:
        rxn.bounds = (-glc, 0)

    # amino acids - trace
    aa = float(sim.get("amino_acid_uptake", 0.01))
    for rid in ESSENTIAL_AA:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-aa, 0)

    # vitamins - trace
    vit = float(sim.get("vitamin_uptake", 0.01))
    for rid in VITAMINS:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds = (-vit, 0)


def get_pathway_fluxes(model, solution):
    """Pull out fluxes for the pathway reactions we care about."""
    out = {}
    for rid, label in PATHWAYS:
        if rid in model.reactions:
            out[f"pathway_{rid}"] = float(solution.fluxes.get(rid, float("nan")))
        else:
            out[f"pathway_{rid}"] = float("nan")
    return out


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]

    scfa_df = read_scfa_inputs(paths.results / "scfa_inputs_canonical.csv", conditions)

    # load model
    sbml = decompress_gz(paths.sbml_path)
    console.log(f"Loading {sbml} ...")
    model = read_sbml_model(str(sbml))
    console.log(f"Loaded: {len(model.reactions)} rxns, {len(model.metabolites)} mets")

    # set up medium
    setup_medium(model, cfg)

    # set ATPM as objective
    atpm = find_rxn(model, ATPM_IDS)
    if not atpm:
        raise RuntimeError("Can't find ATPM reaction")
    atpm.bounds = (0, 500)
    model.objective = atpm
    console.log(f"Objective: {atpm.id}")

    # find SCFA exchange rxns
    scfa_rxns = {}
    for name, cands in SCFA_RXN_IDS.items():
        scfa_rxns[name] = find_rxn(model, cands)
        found = scfa_rxns[name].id if scfa_rxns[name] else "NOT FOUND"
        console.log(f"  {name}: {found}")

    rxn_glc = find_rxn(model, GLUCOSE_IDS)
    rxn_o2 = find_rxn(model, O2_IDS)

    # baseline (no SCFAs)
    console.log("Running baseline...")
    baseline_sol = model.optimize()
    baseline_atpm = float(baseline_sol.objective_value or 0)
    console.log(f"  Baseline ATPM = {baseline_atpm:.4f}")

    # run each condition
    rows = []
    for _, row in scfa_df.iterrows():
        cond = row["condition"]
        ac = float(row["acetate_mmol_gDW_hr"])
        ppa = float(row["propionate_mmol_gDW_hr"])
        but = float(row["butyrate_mmol_gDW_hr"])

        console.log(f"Running {cond}...")

        with model:  # context manager resets bounds after
            for scfa_name, dose in [("acetate", ac), ("propionate", ppa), ("butyrate", but)]:
                rxn = scfa_rxns.get(scfa_name)
                if rxn and dose > 0:
                    rxn.bounds = (-dose, 0)

            sol = model.optimize()
            obj = float(sol.objective_value or 0)
            delta = obj - baseline_atpm
            pct = 100 * delta / baseline_atpm if baseline_atpm else float("nan")

            def flux(r):
                if r is None: return None
                return float(sol.fluxes.get(r.id, 0))

            result = {
                "condition": cond,
                "objective_id": atpm.id,
                "objective_value": obj,
                "baseline_objective": baseline_atpm,
                "objective_delta": delta,
                "objective_pct_change": pct,
                "glucose_flux": flux(rxn_glc),
                "oxygen_flux": flux(rxn_o2),
                "co2_flux": float(sol.fluxes.get("EX_co2_e", 0)),
                "acetate_flux": flux(scfa_rxns.get("acetate")),
                "propionate_flux": flux(scfa_rxns.get("propionate")),
                "butyrate_flux": flux(scfa_rxns.get("butyrate")),
            }
            result.update(get_pathway_fluxes(model, sol))
            rows.append(result)

            console.log(f"  ATPM={obj:.2f} ({pct:+.1f}%)")

    # save
    host_df = pd.DataFrame(rows)
    host_df.to_csv(paths.results / "host_fluxes_by_condition.csv", index=False)
    console.log(f"Saved host_fluxes_by_condition.csv")

    merged = scfa_df.merge(host_df, on="condition", how="left")
    merged.to_csv(paths.results / "merged_dose_scfa_host.csv", index=False)
    console.log(f"Saved merged_dose_scfa_host.csv")

    # print summary
    print("\n--- Summary ---")
    for _, r in host_df.iterrows():
        print(f"  {r['condition']:25s} ATPM={r['objective_value']:.2f}  delta={r['objective_pct_change']:+.1f}%")


if __name__ == "__main__":
    main()
