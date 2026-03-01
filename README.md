# Stachys affinis SCFA modeling - reproducibility

Code and data to reproduce the constraint-based modeling results from the manuscript.

We used COBRApy + Recon3D to simulate how short-chain fatty acids derived from
stachyose fermentation affect hepatic energy metabolism (ATP maintenance flux).

## Setup

You need conda (or mamba) and Python 3.11.

```
conda env create -f environment.yml
conda activate stachys-scfa
```

You also need the Recon3D SBML file. Download from BiGG:
https://bigg.ucsd.edu/models/Recon3D

Put `Recon3D.xml.gz` in `data/models/`.

## How to run

```
make all
```

Or run each step separately:
```
python -m src.01_prepare_inputs
python -m src.02_run_simulation
python -m src.03_figures
python -m src.04_tables
```

Step 01 validates the SCFA input file. Step 02 does the actual FBA (takes a couple
minutes to load the model). Steps 03 and 04 generate the figures and tables from the
simulation output.

## Outputs

- `results/` -- simulation CSVs (fluxes, merged data)
- `outputs/figs/` -- the 6 figures as PNGs
- `outputs/tables/` -- 3 tables as CSVs

## SCFA doses

We defined three dose conditions based on estimated S. affinis tuber intake:

| Condition | Acetate | Propionate | Butyrate |
|-----------|---------|------------|----------|
| Low (~25g tuber)  | 2.0 | 0.7 | 0.4 |
| Mid (~50g tuber)  | 4.0 | 1.4 | 0.8 |
| High (~100g tuber) | 8.0 | 2.8 | 1.6 |

(all in mmol/gDW/hr)

## Notes

- The medium constraints are pretty strict -- we close all 1806 boundary rxns and
  only reopen a small curated set. This is necessary to avoid thermodynamic loops
  in Recon3D that give unrealistic ATP yields.
- Propionate shows zero flux under all conditions. This is a known Recon3D issue
  (the methylmalonyl-CoA mutase pathway doesn't carry flux with these constraints).
- ATPM is the objective function, not biomass, since hepatocytes don't divide.
