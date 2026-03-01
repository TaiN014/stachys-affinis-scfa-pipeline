"""
Step 03 - generate figures from simulation results.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rich.console import Console

from .utils import build_paths

console = Console()

# colors (roughly colorbrewer dark2)
COLORS = {
    "acetate": "#1b9e77",
    "propionate": "#d95f02",
    "butyrate": "#7570b3",
    "glucose": "#e7298a",
    "oxygen": "#66a61e",
    "co2": "#a65628",
    "objective": "#e6ab02",
    "baseline": "#999999",
    "delta": "#d62728",
}

COND_ORDER = ["StachysDose_Low", "StachysDose_Mid", "StachysDose_High"]
LABELS = ["Low", "Mid", "High"]
FIGSIZE = (7, 4.5)
DPI = 300

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
})


def main():
    paths = build_paths()
    merged = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # sort by dose level
    merged["_ord"] = merged["condition"].map({c: i for i, c in enumerate(COND_ORDER)})
    merged = merged.sort_values("_ord").reset_index(drop=True)
    x = np.arange(len(merged))

    baseline = merged["baseline_objective"].iloc[0]

    # ---- Fig 1: SCFA inputs ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for scfa, col in [("acetate", "acetate_mmol_gDW_hr"),
                       ("propionate", "propionate_mmol_gDW_hr"),
                       ("butyrate", "butyrate_mmol_gDW_hr")]:
        ax.plot(x, merged[col], "o-", color=COLORS[scfa], label=scfa.capitalize(),
                lw=2, markersize=7)
    ax.set_xticks(x); ax.set_xticklabels(LABELS)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("SCFA Availability (mmol/gDW/hr)")
    ax.set_title("SCFA Availability by Dose Condition")
    ax.legend()
    fig.savefig(paths.figs_dir / "fig_scfa_inputs.png")
    plt.close(fig)
    console.log("fig_scfa_inputs.png")

    # ---- Fig 2: SCFA ratios (stacked bar) ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    totals = merged["acetate_mmol_gDW_hr"] + merged["propionate_mmol_gDW_hr"] + merged["butyrate_mmol_gDW_hr"]
    fa = merged["acetate_mmol_gDW_hr"] / totals
    fp = merged["propionate_mmol_gDW_hr"] / totals
    fb = merged["butyrate_mmol_gDW_hr"] / totals
    ax.bar(x, fa, 0.55, color=COLORS["acetate"], label="Acetate")
    ax.bar(x, fp, 0.55, bottom=fa, color=COLORS["propionate"], label="Propionate")
    ax.bar(x, fb, 0.55, bottom=fa+fp, color=COLORS["butyrate"], label="Butyrate")
    ax.set_xticks(x); ax.set_xticklabels(LABELS)
    ax.set_ylabel("Molar Fraction")
    ax.set_title("SCFA Molar Ratio by Condition")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper right")
    fig.savefig(paths.figs_dir / "fig_scfa_ratios.png")
    plt.close(fig)
    console.log("fig_scfa_ratios.png")

    # ---- Fig 3: ATPM bar chart ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.bar(x, merged["objective_value"], 0.55, color=COLORS["objective"])
    ax.axhline(baseline, color=COLORS["baseline"], ls="--", lw=1.5,
               label=f"Baseline ({baseline:.1f})")
    for i, v in enumerate(merged["objective_value"]):
        ax.text(i, v + 1, f"{v:.1f}", ha="center", fontweight="bold", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(LABELS)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("ATPM (mmol/gDW/hr)")
    ax.set_title("ATP Maintenance by Dose Condition")
    ax.legend()
    fig.savefig(paths.figs_dir / "fig_host_objective.png")
    plt.close(fig)
    console.log("fig_host_objective.png")

    # ---- Fig 4: %change ----
    if "objective_pct_change" in merged.columns:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        ax.bar(x, merged["objective_pct_change"], 0.55, color=COLORS["delta"])
        for i, v in enumerate(merged["objective_pct_change"]):
            ax.text(i, v + 2, f"+{v:.1f}%", ha="center", fontweight="bold", fontsize=10)
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(LABELS)
        ax.set_ylabel("ATPM Change (%)")
        ax.set_title("% Change in ATPM vs Baseline")
        fig.savefig(paths.figs_dir / "fig_objective_delta_pct.png")
        plt.close(fig)
        console.log("fig_objective_delta_pct.png")

    # ---- Fig 5: exchange fluxes ----
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for label, col, color in [
        ("Glucose", "glucose_flux", COLORS["glucose"]),
        ("Acetate", "acetate_flux", COLORS["acetate"]),
        ("Propionate", "propionate_flux", COLORS["propionate"]),
        ("Butyrate", "butyrate_flux", COLORS["butyrate"]),
        ("O₂", "oxygen_flux", COLORS["oxygen"]),
        ("CO₂", "co2_flux", COLORS["co2"]),
    ]:
        if col in merged.columns:
            ax.plot(x, merged[col], "o-", color=color, label=label, lw=2, markersize=6)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(LABELS)
    ax.set_ylabel("Flux (mmol/gDW/hr)")
    ax.set_title("Exchange Fluxes by Condition")
    ax.legend(ncol=2, loc="lower left", fontsize=9)
    fig.savefig(paths.figs_dir / "fig_host_exchange_fluxes.png")
    plt.close(fig)
    console.log("fig_host_exchange_fluxes.png")

    # ---- Fig 6: pathway heatmap ----
    pw_cols = [c for c in merged.columns if c.startswith("pathway_")]
    if pw_cols:
        pw = merged[pw_cols].copy()
        pw.index = LABELS
        pw.columns = [c.replace("pathway_", "") for c in pw_cols]
        # drop all-zero or all-NaN columns
        pw = pw.loc[:, (pw != 0).any() & pw.notna().any()]

        if not pw.empty:
            fig, ax = plt.subplots(figsize=(max(7, pw.shape[1] * 0.9), 4.5))
            im = ax.imshow(pw.values.astype(float), aspect="auto", cmap="RdYlBu_r")
            ax.set_xticks(np.arange(pw.shape[1]))
            ax.set_xticklabels(pw.columns, rotation=45, ha="right")
            ax.set_yticks(np.arange(len(pw)))
            ax.set_yticklabels(pw.index)
            ax.set_title("Pathway Fluxes by Condition")

            # annotate cells
            for i in range(len(pw)):
                for j in range(pw.shape[1]):
                    v = pw.iloc[i, j]
                    if pd.notna(v) and v != 0:
                        color = "white" if v > pw.values.max() * 0.6 else "black"
                        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                                fontsize=9, color=color)

            plt.colorbar(im, ax=ax, shrink=0.8)
            fig.tight_layout()
            fig.savefig(paths.figs_dir / "fig_pathway_heatmap.png")
            plt.close(fig)
            console.log("fig_pathway_heatmap.png")

    console.log("Done with figures.")


if __name__ == "__main__":
    main()
