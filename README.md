# Gradient Random Walk Solvers for the Heat, FitzHugh–Nagumo, and Burgers' Equations

Python software and reproducible numerical examples for the paper
*On the Accuracy of Gradient Random Walk Methods for the Heat,
FitzHugh–Nagumo, and Burgers' Equations* by Stephen Abkin and
Prabir Daripa ([arXiv:2608.22592](https://doi.org/10.48550/arXiv.2608.22592)).

The repository supports two uses:

1. reproduce the reported tables and figures, and
2. modify the supplied configurations or study scripts to run new cases.

## Install

For the archived version 1.0.0 release, download and extract the ZIP from
[Zenodo](https://doi.org/10.5281/zenodo.22050659), then open a terminal in the
extracted directory:

```bash
cd Gradient-Random-Walk-Solvers-1.0.0
```

Alternatively, clone the development repository:

```bash
git clone https://github.com/stephen122204/Gradient-Random-Walk-Solvers.git
cd Gradient-Random-Walk-Solvers
git checkout grw-solvers-v3
```

Create a Python 3.11 environment:

```bash
python -m venv .venv
```

Activate it on macOS or Linux:

```bash
source .venv/bin/activate
```

On Windows PowerShell, use:

```powershell
.\.venv\Scripts\Activate.ps1
```

On Windows Command Prompt, use:

```bat
.venv\Scripts\activate.bat
```

Then install the pinned dependencies on any platform:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The pinned environment uses Python 3.11.4. Generated files are written under
`output/` or `outputs/`. Both directories are ignored by Git.

## Check the Installation

Confirm the environment works before the longer runs (takes a few seconds):

```bash
python -m unittest discover -s tests
```

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
control), `t5` (FitzHugh–Nagumo), `t3` (Cole–Hopf plateau controls), and
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
simulations as well as the tabulated studies. One command covers everything,
the deep representative checks followed by a full re-run and comparison of the
five ensemble studies:

```bash
python reproduce.py verify-all
```

Measured wall-clock times (Apple-Silicon laptop, pinned environment):

| Command | What it covers | Wall clock |
|---|---|---|
| `python reproduce.py paper` | the paper's ten figures, from committed data | ~6 s |
| `python reproduce.py verify --deep` | single-seed studies plus archived representative arrays (179 checks) | ~17 s |
| `python reproduce.py t3` | Cole–Hopf plateau controls | ~2 s |
| `python reproduce.py t8` | controlled Burgers attribution | ~4 s |
| `python reproduce.py t5` | FitzHugh–Nagumo thirty-seed ensemble | ~15 s |
| `python reproduce.py t4` | heat thirty-seed ensemble | ~1.5 min |
| `python reproduce.py t7` | paired heat output-grid study | ~1.5 min |
| `python reproduce.py verify-all` | release gate: deep checks plus all five ensemble studies, re-run and compared | ~4 min |

Allow longer on older hardware.

Run `python reproduce.py` with no target to display every available command.

## Additional Checks

Two standalone scripts in `checks/` verify specific parts of the code and the
archived data. They are not `reproduce.py` targets; run them from the
repository root.

| Command | What it covers | Wall clock |
|---|---|---|
| `python checks/heat_cdf_identity_check.py` | evaluates the deterministic bin-and-sum control, the empirical-CDF identity for the expected squared errors, the Gaussian/delta approximations to their sampling spread, and the paired identity, and compares them with `pinned_ensembles/heat_grid_paired/` (no solver runs) | <1 s |
| `python checks/heat_cdf_identity_check.py --fresh` | additionally re-runs the heat walk for four independent seed blocks to assess seed-block variability | ~15 s |
| `python checks/interface_regression_checks.py` | confirms the FitzHugh–Nagumo default path is unchanged by the optional `reaction_derivative` callback and that `utils.reconstruct_cumulative` reproduces the pinned paired-heat statistics | ~1 min |

`interface_regression_checks.py --pristine DIR` also compares the default
path against an earlier checkout at `DIR`.

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
- `tests/`: quick installation checks of core formulas and boundary operations.
- `checks/`: standalone identity and regression checks (see Additional Checks).


## Citation

Please cite the paper:

> S. Abkin and P. Daripa, *On the Accuracy of Gradient Random Walk Methods
> for the Heat, FitzHugh–Nagumo, and Burgers' Equations*, arXiv preprint
> [arXiv:2608.22592](https://doi.org/10.48550/arXiv.2608.22592), 2026.

Version 1.0.0 of this software is archived on Zenodo:
[https://doi.org/10.5281/zenodo.22050659](https://doi.org/10.5281/zenodo.22050659).
See `CITATION.cff` for complete citation metadata.

## Acknowledgments

The authors thank Oliver Stalker for providing an early version of the
Python code.

**Principal Investigator:** [Professor Prabir Daripa](https://artsci.tamu.edu/mathematics/contact/profiles/prabir-daripa.html) — Texas A&M University, Department of Mathematics

Other projects from the Daripa Research Group are available on the
[group's GitHub page](https://github.com/Daripa-Research-Group).
