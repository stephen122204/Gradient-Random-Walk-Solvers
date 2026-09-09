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

Use the current branch for the revised paper:

```bash
git clone --branch grw-solvers-v3 https://github.com/stephen122204/Gradient-Random-Walk-Solvers.git
cd Gradient-Random-Walk-Solvers
```

If you downloaded the updated submission supplement instead, extract it and
open a terminal inside its `source/` directory. The commands below work from
that directory as well.

The original [Zenodo 1.0.0 archive](https://doi.org/10.5281/zenodo.22050659)
supports the original paper. The two-step predictive extension requires this
updated checkout or supplement. A new archive version has not yet been deposited.

Check that `python --version` reports Python 3.11, then create an environment:

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
control), `t5` (scalar reaction–diffusion), `t3` (Cole–Hopf plateau controls), and
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
| `python reproduce.py t5` | scalar reaction–diffusion thirty-seed ensemble | ~15 s |
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

The first reproduces the direct-distribution and finite-ensemble moment controls
and six independent seed blocks. Its default reference control uses 20 million samples and can use
substantial memory. The second reproduces the bootstrap comparison table,
including heat profile statistics, and creates the requested JSON parent directory. The manuscript uses joint-by-seed intervals
as primary. Both scripts print their results, and the second also writes JSON.

`interface_regression_checks.py --reference DIR` additionally compares the
default path with another checkout of the code at `DIR`.

## Outputs, Interrupted Runs, and Random Seeds

Paper studies use fixed output names. Rerunning a target replaces its generated
files. The `paper` target recreates its figure directory from committed inputs.
The pinned inputs in `figure_data/`, `pinned_ensembles/`, and `provenance/` remain
unchanged during the usual reproduction commands.

| Run | Destination and overwrite behavior |
|---|---|
| `studies` or `verify` | `output/paper_refinement_original_grw/`; fixed study files |
| `t3`, `t4`, `t5`, `t7`, `t8` or `verify-ensembles` | One directory per study under `output/final_prepublication_tests/` |
| `paper` | Recreates `output/final_prepublication_tests/paper_figures/` |
| `verify-two-step` or `t9` | `output/heat_two_step_reproduction/`; fixed prediction and ensemble files |
| Direct two-step script | `output/heat_two_step/` by default, or `--output-dir PATH` |
| `main.py CONFIG` | Plots in `outputs/YYYY-MM-DD_HH-MM-SS/`; the folder identifies the run time, not its parameters |

There is no checkpoint resume within an ensemble. After an interruption, rerun
the affected study. Completed original studies can be checked together without
rerunning them using `python reproduce.py verify-ensembles --no-rerun`.
For the two-step study, if `predict` finished, rerun the interrupted `run` or
`validate` stage with the same output directory. Once all stages finish,
`python reproduce.py verify-two-step --no-rerun` checks the standard reproduction
directory. A failed or interrupted command does not establish a completed check.

Each study fixes its seed list and reuses it across refinement levels. Changing
the particle count can change how draws are assigned to particles, so seed reuse
is not the same as identical trajectories. Paired heat reconstructions reuse the
same final particle set across comparison points and references. The two-step
production uses seeds 5000–5029, validation uses 7000–7029, and the six additional
blocks use separate seed ranges. The manuscript's primary bootstrap resamples
jointly by seed across refinement levels. Historical study outputs retain their
earlier intervals; `bootstrap_seed_grouping_check.py` produces the current
`joint` intervals and the comparison table.

The exploratory `main.py` command does not set a seed. To make a custom run
repeatable, set NumPy's seed before calling a solver and save the configuration
with the results. Matching a seed alone does not pair runs of different equations
or guarantee matching particle trajectories after changing the design.

## Add Your Own Case

Three levels of customization are available. Pinned verification checks apply
to the published configuration. Changing physical parameters, particle counts,
seeds, or reconstruction conventions changes the expected results and requires
new figures and summaries for that case. Keep custom outputs in a separate
directory. Study filenames do not automatically encode modified parameters.
Larger particle counts, more seeds, or more time steps increase the work.

**A new parameter set.** Copy a JSON file from `configs/`, edit it, and pass it
to the solver:

```bash
python main.py configs/heat_step_dirichlet.json
python main.py configs/fhn_grw_steady.json
python main.py configs/burgers_stationary_shock.json
```

`config_template.jsonc` documents every field. Figures are written below
`outputs/`. Save your edited configuration with the output because the timestamp
folder does not record the parameter values. When an exact solution is
available, the run can be checked with:

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
and controlled comparisons. Copy one and edit its parameter block and output
destination. The two-step study requires an even particle count of at least
two for its equal group allocation. To compare
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
See `CITATION.cff` for the original release metadata. Update its version, date,
and DOI when the revised software archive is deposited.

## Acknowledgments

The authors thank Oliver Stalker for providing an early version of the
Python code.

**Principal Investigator:** [Professor Prabir Daripa](https://artsci.tamu.edu/mathematics/contact/profiles/prabir-daripa.html) — Texas A&M University, Department of Mathematics

Other projects from the Daripa Research Group are available on the
[group's GitHub page](https://github.com/Daripa-Research-Group).
