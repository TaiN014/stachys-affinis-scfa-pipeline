#!/usr/bin/env python3
"""
03_figures.py -- generate figures from simulation results.
"""

import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt

from .utils import build_paths

# colors
COLORS = {
    "acetate":    "#1b9e77",
    "propionate": "#d95f02",
    "butyrate":   "#7570b3",
    "glucose":    "#e7298a",
    "oxygen":     "#66a61e",
    "co2":        "#a65628",
    "objective":  "#e6ab02",
    "baseline":   "#999999",
    "delta":      "#d62728",
}


COND_ORDER = ["StachysDose_High", "StachysDose_Mid", "StachysDose_Low"]
XLABELS = ["High", "Mid", "Low"]

FIG_W, FIG_H = 7, 4.5
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
    df = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # sort
    order_map = {c: i for i, c in enumerate(COND_ORDER)}
    df["_ord"] = df["condition"].map(order_map)
    df = df.sort_values("_ord").reset_index(drop=True)
    x = np.arange(len(df))

    baseline_val = df["baseline_objective"].iloc[0]

    # fig 1: SCFA inputs
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    for scfa, col_name in [("acetate", "acetate_mmol_gDW_hr"),
                            ("propionate", "propionate_mmol_gDW_hr"),
                            ("butyrate", "butyrate_mmol_gDW_hr")]:
        ax.plot(x, df[col_name], "o-", color=COLORS[scfa],
                label=scfa.capitalize(), lw=2, markersize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("SCFA Availability (mmol/gDW/hr)")
    ax.set_title("SCFA Availability by Dose Condition")
    ax.legend()
    fig.savefig(paths.figs_dir / "fig_scfa_inputs.png")
    plt.close(fig)
    print("  fig_scfa_inputs.png")

    # fig 2: SCFA ratios
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    total = (df["acetate_mmol_gDW_hr"] + df["propionate_mmol_gDW_hr"]
             + df["butyrate_mmol_gDW_hr"])
    frac_ac  = df["acetate_mmol_gDW_hr"] / total
    frac_ppa = df["propionate_mmol_gDW_hr"] / total
    frac_but = df["butyrate_mmol_gDW_hr"] / total

    bar_w = 0.55
    ax.bar(x, frac_ac, bar_w, color=COLORS["acetate"], label="Acetate")
    ax.bar(x, frac_ppa, bar_w, bottom=frac_ac,
           color=COLORS["propionate"], label="Propionate")
    ax.bar(x, frac_but, bar_w, bottom=frac_ac + frac_ppa,
           color=COLORS["butyrate"], label="Butyrate")
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_ylabel("Molar Fraction")
    ax.set_title("SCFA Molar Ratio by Condition")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper right")
    fig.savefig(paths.figs_dir / "fig_scfa_ratios.png")
    plt.close(fig)
    print("  fig_scfa_ratios.png")

    # fig 3: ATPM
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.bar(x, df["objective_value"], 0.55, color=COLORS["objective"])
    ax.axhline(baseline_val, color=COLORS["baseline"], ls="--", lw=1.5,
               label=f"Baseline ({baseline_val:.1f})")

    for i, val in enumerate(df["objective_value"]):
        ax.text(i, val + 1, f"{val:.1f}", ha="center", fontweight="bold",
                fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_xlabel("Stachyose Dose")
    ax.set_ylabel("ATPM (mmol/gDW/hr)")
    ax.set_title("ATP Maintenance by Dose Condition")
    ax.legend()
    fig.savefig(paths.figs_dir / "fig_host_objective.png")
    plt.close(fig)
    print("  fig_host_objective.png")

    # fig 4: % change
    if "objective_pct_change" in df.columns:
        fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
        ax.bar(x, df["objective_pct_change"], 0.55, color=COLORS["delta"])
        for i, val in enumerate(df["objective_pct_change"]):
            ax.text(i, val + 2, f"+{val:.1f}%", ha="center",
                    fontweight="bold", fontsize=10)
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(XLABELS)
        ax.set_ylabel("ATPM Change (%)")
        ax.set_title("% Change in ATPM vs Baseline")
        fig.savefig(paths.figs_dir / "fig_objective_delta_pct.png")
        plt.close(fig)
        print("  fig_objective_delta_pct.png")

    # fig 5: exchange fluxes
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    flux_series = [
        ("Glucose",    "glucose_flux",    COLORS["glucose"]),
        ("Acetate",    "acetate_flux",    COLORS["acetate"]),
        ("Propionate", "propionate_flux", COLORS["propionate"]),
        ("Butyrate",   "butyrate_flux",   COLORS["butyrate"]),
        ("O\u2082",    "oxygen_flux",     COLORS["oxygen"]),
        ("CO\u2082",   "co2_flux",        COLORS["co2"]),
    ]
    for label, col, color in flux_series:
        if col in df.columns:
            ax.plot(x, df[col], "o-", color=color, label=label,
                    lw=2, markersize=6)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(XLABELS)
    ax.set_ylabel("Flux (mmol/gDW/hr)")
    ax.set_title("Exchange Fluxes by Condition")
    ax.legend(ncol=2, loc="lower left", fontsize=9)
    fig.savefig(paths.figs_dir / "fig_host_exchange_fluxes.png")
    plt.close(fig)
    print("  fig_host_exchange_fluxes.png")

    # fig 6: pathway heatmap
    pw_cols = [c for c in df.columns if c.startswith("pathway_")]
    if pw_cols:
        pw = df[pw_cols].copy()
        pw.index = XLABELS
        pw.columns = [c.replace("pathway_", "") for c in pw_cols]

        # drop all-zero or all-NaN
        pw = pw.loc[:, (pw != 0).any() & pw.notna().any()]

        if not pw.empty:
            ncols = pw.shape[1]
            fig_width = max(7, ncols * 0.9)
            fig, ax = plt.subplots(figsize=(fig_width, 4.5))

            im = ax.imshow(pw.values.astype(float), aspect="auto",
                           cmap="RdYlBu_r")
            ax.set_xticks(np.arange(ncols))
            ax.set_xticklabels(pw.columns, rotation=45, ha="right")
            ax.set_yticks(np.arange(len(pw)))
            ax.set_yticklabels(pw.index)
            ax.set_title("Pathway Fluxes by Condition")

            # annotate
            vmax = pw.values.max()
            for i in range(len(pw)):
                for j in range(ncols):
                    v = pw.iloc[i, j]
                    if pd.notna(v) and v != 0:
                        txt_color = "white" if v > vmax * 0.6 else "black"
                        ax.text(j, i, f"{v:.1f}", ha="center", va="center",
                                fontsize=9, color=txt_color)

            plt.colorbar(im, ax=ax, shrink=0.8)
            fig.tight_layout()
            fig.savefig(paths.figs_dir / "fig_pathway_heatmap.png")
            plt.close(fig)
            print("  fig_pathway_heatmap.png")

    print("done with figures")


if __name__ == "__main__":
    main()
