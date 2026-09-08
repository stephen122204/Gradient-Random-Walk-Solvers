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

Two scripts in `checks/` evaluate the paper's analytical identities against the
archived data and confirm that two solver entry points reproduce it. Run them
from the repository root.

| Command | What it covers | Wall clock |
|---|---|---|
| `python checks/heat_cdf_identity_check.py` | the deterministic bin-and-sum control, the empirical-CDF identity for the expected squared errors, their Gaussian/delta sampling approximations, and the paired identity, compared with `pinned_ensembles/heat_grid_paired/` (no solver runs) | <1 s |
| `python checks/heat_cdf_identity_check.py --fresh` | also runs the heat walk for four independent seed blocks to show seed-block variability of the spread | ~15 s |
| `python checks/interface_regression_checks.py` | the `reaction_derivative` callback reproduces the default reaction–diffusion path exactly, and `utils.reconstruct_cumulative` reproduces the pinned paired-heat statistics | ~1 min |

`interface_regression_checks.py --reference DIR` additionally compares the
default path with another checkout of the code at `DIR`.

## Add Your Own Case

Three levels of customization are available.

**A new parameter set.** Copy a JSON file from `configs/`, edit it, and pass it
to the solver:

```bash
python main.py configs/heat_step_dirichlet.json
python main.py configs/fhn_grw_steady.json
python main.py configs/burgers_stationary_shock.json
```

`config_template.jsonc` documents every field. Figures are written below
`outputs/`. When an exact solution is available, the run can be checked with:

```bash
python verify_solver.py --equation heat --config configs/heat_step_dirichlet.json
```

**A new scalar reaction law.** The reaction–diffusion solver evolves the
gradient of `u_t = D u_xx + f(u)` and needs only `f'(u)`. Set
`reaction_derivative` on the configuration to a function of a NumPy array;
leaving it `None` selects the built-in Nagumo-type polynomial used in the paper.

```python
from config import SimulationConfig
cfg = SimulationConfig(equation_type='fitzhugh-nagumo', ...)
cfg.reaction_derivative = lambda u: 1.0 - 2.0 * u      # Fisher–KPP: f(u) = u(1 - u)
```

**A new study.** The files in `studies/` are complete examples of parameter
sweeps, multi-seed ensembles, bias–spread decompositions, bootstrap intervals,
and controlled comparisons. Copy one and edit its parameter block. To compare
a binned cumulative reconstruction with a reference, use
`utils.reconstruct_cumulative`, which returns the coordinates the
reconstruction represents along with the values.

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
