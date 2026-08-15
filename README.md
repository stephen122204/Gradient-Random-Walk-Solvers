# Gradient Random Walk Methods for the Heat, FitzHugh–Nagumo, and Viscous Burgers Equations

Python software and reproducible numerical examples for the paper
*Controlled Error Attribution for Gradient Random Walk Methods Applied to the
Heat, FitzHugh–Nagumo, and Viscous Burgers Equations* by Stephen Abkin and
Prabir Daripa.

The repository supports two uses:

1. reproduce the reported tables and figures; and
2. modify the supplied configurations or study scripts to run new cases.

## Install

```bash
git clone https://github.com/stephen122204/heat_burgers_fhn.git
cd heat_burgers_fhn
git checkout grw-solvers-v1
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The pinned environment uses Python 3.11.4. Generated files are written under
`output/` or `outputs/`; both directories are ignored by Git.

## Reproduce the Paper

Generate the ten figures used in the paper directly from the committed data:

```bash
python reproduce.py paper
```

Re-run the representative single-seed studies or all ensemble and controlled
studies:

```bash
python reproduce.py studies
python reproduce.py ensembles
```

The individual ensemble targets are `t4` (heat), `t7` (paired heat-grid
control), `t5` (FitzHugh–Nagumo), `t3` (original Cole–Hopf diagnostics), and
`t8` (controlled Burgers attribution). For example:

```bash
python reproduce.py t8
```

The committed reference values can be checked with:

```bash
python reproduce.py verify
python reproduce.py verify-ensembles
```

Use `python reproduce.py verify --deep` to re-run the archived representative
simulations as well as the tabulated studies. Run `python reproduce.py` with no
target to display every available command.

## Run a Modified Case

Copy a JSON file from `configs/`, change its parameters, and pass it to the
solver:

```bash
python main.py configs/heat_step_dirichlet.json
python main.py configs/fhn_grw_steady.json
python main.py configs/burgers_stationary_shock.json
```

`config_template.jsonc` documents the available fields. Custom comparison
figures are saved below `outputs/`. When an exact solution is available, a
modified case can also be checked with:

```bash
python verify_solver.py --equation heat --config configs/heat_step_dirichlet.json
```

The files in `studies/` are complete examples of parameter sweeps, multi-seed
experiments, error decompositions, bootstrap intervals, and controlled
comparisons. They can be copied and edited for new studies.

## Repository Layout

- `simulation.py`, `config.py`: solver and configuration handling.
- `configs/`, `config_template.jsonc`: editable example inputs.
- `studies/`, `study_paper_refinement.py`: paper experiments and reusable study
  examples.
- `reproduce.py`: paper reproduction and verification entry point.
- `figure_data/`, `pinned_ensembles/`, `expected_values.json`: committed data
  behind the reported values and figures.
- `figure_scripts/`: figure generation.
- `tests/`: small smoke tests for core formulas and boundary operations.

## Optional Smoke Test

```bash
python -m unittest discover -s tests
```

## Citation

See `CITATION.cff`.
