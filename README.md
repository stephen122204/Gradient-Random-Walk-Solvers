# Gradient Random Walk Solvers

Python software and reproducible numerical examples for the revised paper
*Error attribution in gradient random walk methods for parabolic equations*
by Stephen Abkin and Prabir Daripa. The earlier preprint is
[arXiv:2608.22592](https://doi.org/10.48550/arXiv.2608.22592), under its original
heat, FitzHugh-Nagumo, and Burgers title.

The reaction example is a scalar reaction-diffusion equation with a
FitzHugh-Nagumo-type cubic term. The existing `fitzhugh-nagumo` configuration
key is retained for compatibility.

The repository supports two uses:

1. reproduce the reported tables and figures, and
2. modify the supplied configurations or study scripts to run new cases.

## Install

**For the revised paper, use this updated checkout or the updated submission
supplement.** The original Zenodo 1.0.0 archive predates the two-step predictive
extension. A new archive version has not yet been deposited.


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

Generate all eleven figures used in the paper directly from the committed data:

```bash
python reproduce.py paper
```

Re-run the representative single-seed studies or the five original ensemble
and controlled studies:

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

The new heat prediction study is included in `verify-all` and can also be
reproduced on its own:

```bash
python reproduce.py verify-two-step
```

This generates predictions and the validation specification before running the
production and validation ensembles, then compares the outputs with the
committed data. `python reproduce.py t9` is an alias. Fresh files go to
`output/heat_two_step_reproduction/`, preserving the original records. The
particle count is recomputed from the rule and rounded to an even allocation.

For a fast check of archived predictions, historical hashes, the count rule,
and error-moment consistency without solver runs:

```bash
python reproduce.py verify-two-step --pinned-only
```

To regenerate the extension's figure and comparison summary from the fresh run:

```bash
python studies/analyze_t9_heat_two_step.py --data-dir output/heat_two_step_reproduction --output-dir output/heat_two_step_reproduction
```

See [Methods and reproduction](docs/METHODS_AND_REPRODUCIBILITY.md) for the
mathematical sources, assumptions, and complete mapping from paper results to
commands. The predictive study uses classical statistics and reflected heat
solutions to diagnose comparison errors and select particle effort. The
reaction multiplier follows Ghoniem and Sherman (1985), equation (41).

Use `python reproduce.py verify --deep` to re-run the archived representative
simulations as well as the tabulated studies. The combined numerical check
runs the deep representative checks, the five original ensemble studies, and
the new two-step heat extension, comparing the results with committed data:

```bash
python reproduce.py verify-all
```

Measured wall-clock times (Apple-Silicon laptop, pinned environment):

| Command | What it covers | Wall clock |
|---|---|---|
| `python reproduce.py paper` | the paper's eleven figures, from committed data | ~6 s |
| `python reproduce.py verify --deep` | single-seed studies plus archived representative arrays (179 checks) | ~17 s |
| `python reproduce.py t3` | Cole–Hopf plateau controls | ~2 s |
| `python reproduce.py t8` | controlled Burgers attribution | ~4 s |
| `python reproduce.py t5` | FitzHugh–Nagumo thirty-seed ensemble | ~15 s |
| `python reproduce.py t4` | heat thirty-seed ensemble | ~1.5 min |
| `python reproduce.py t7` | paired heat output-grid study | ~1.5 min |
| `python reproduce.py verify-all` | deep checks, five original ensemble studies, and two-step heat, re-run and compared | several minutes; hardware dependent |

Allow longer on older hardware.

Run `python reproduce.py` with no target to display every available command.

## Additional Checks

The following scripts reproduce the analytical, uncertainty, and interface
checks. Run them from the repository root. These are separate from the
numerical rerun gate `verify-all`.

| Command | What it covers | Wall clock |
|---|---|---|
| `python checks/heat_cdf_identity_check.py` | the deterministic bin-and-sum control, the empirical-CDF identity for the expected squared errors, their Gaussian/delta sampling approximations, and the paired identity, compared with `pinned_ensembles/heat_grid_paired/` (no solver runs) | <1 s |
| `python checks/heat_cdf_identity_check.py --fresh` | also runs the heat walk for four independent seed blocks to show seed-block variability of the spread | ~15 s |
| `python checks/interface_regression_checks.py` | the `reaction_derivative` callback reproduces the default reaction–diffusion path exactly, and `utils.reconstruct_cumulative` reproduces the pinned paired-heat statistics | ~1 min |

The revision also uses these checks:

```bash
python checks/two_step_reference_and_seed_blocks.py
python checks/bootstrap_seed_grouping_check.py --heat-profiles --json output/bootstrap-seed-grouping.json
```

The first reproduces the direct-distribution controls and six independent seed
blocks. Its default reference control uses 20 million samples and can use
substantial memory. The second reproduces the bootstrap comparison table,
including heat profile statistics. The manuscript uses joint-by-seed intervals
as primary. Both scripts print their results, and the second also writes JSON.

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
gradient of `u_t = D u_xx + f(u)` and evaluates `f'(u)` in the weight update.
A new case also requires compatible initial and boundary data and a suitable
time step. Set
`reaction_derivative` on the configuration to a function of a NumPy array;
leaving it `None` selects the built-in Nagumo-type polynomial used in the paper.

```python
from config import load_config_from_json
cfg = load_config_from_json('configs/fhn_grw_steady.json')
# Replace the initial/boundary data and reference to match the new problem.
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
- `checks/`: identity, bootstrap, interface, and two-step reproduction checks.
- `docs/METHODS_AND_REPRODUCIBILITY.md`: mathematical sources and paper-to-command map.
- `provenance/heat_two_step/`: preserved original design hashes and production log.


## Build a Reproduction Package

```bash
python tools/build_code_supplement.py --output output/code-supplement.zip
```

The ZIP contains a complete `source/` tree, the committed data and design
records, and a SHA-256 manifest. Optional `--check-outputs DIR` includes
verification logs saved in that directory. This creates a local package only.

## Citation

The original preprint citation is:

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
