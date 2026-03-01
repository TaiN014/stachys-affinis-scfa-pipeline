SHELL := /bin/bash

.PHONY: all clean

all:
	python -m src.01_prepare_inputs
	python -m src.02_run_simulation
	python -m src.03_figures
	python -m src.04_tables

clean:
	rm -rf results/*.csv outputs/figs/*.png outputs/tables/*.csv data/models/cache/
