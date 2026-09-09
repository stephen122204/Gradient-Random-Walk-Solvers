# Gradient Random Walk Solvers for the Heat, FitzHugh–Nagumo, and Burgers' Equations

Python software and reproducible numerical examples for the paper
*Error attribution in gradient random walk methods for parabolic equations*
by Stephen Abkin and Prabir Daripa.

The repository supports two uses:

1. reproduce the reported tables and figures, and
2. modify the supplied configurations or study scripts to run new cases.

## Install

Clone the revised-paper branch:

```bash
git clone --branch grw-solvers-v3 https://github.com/stephen122204/Gradient-Random-Walk-Solvers.git
cd Gradient-Random-Walk-Solvers
```

Alternatively, extract the code supplement ZIP and open a terminal in its
`source/` folder.

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

Then install the pinned dependencies:

```bash
python -m pip install -r requirements.txt
```

The pinned environment uses Python 3.11.4. Generated files are written under
`output/` or `outputs/`. Both directories are ignored by Git.

## Reproduce the Paper

Generate all eleven figures directly from the committed data:

```bash
python reproduce.py paper
```

Rerun the representative simulations, refinement studies, ensembles, and
two-step heat study, and compare their outputs with the committed results:

```bash
python reproduce.py verify-all
```

Reproduce the additional heat identity controls, independent seed-block
controls, and fitted-rate bootstrap comparison:

```bash
python checks/heat_cdf_identity_check.py --fresh
python checks/two_step_reference_and_seed_blocks.py
python checks/bootstrap_seed_grouping_check.py --heat-profiles --json output/bootstrap.json
```

The bootstrap output's `joint` intervals are those reported in the manuscript.
Figures are saved under `output/final_prepublication_tests/paper_figures/`.
Rerunning replaces generated outputs. Allow several minutes for the full
studies and several GB of available memory for the direct-distribution control.
Runtime depends on hardware.

To run one study, use `t4` for heat, `t7` for paired heat reconstructions,
`t5` for the scalar reaction–diffusion front, `t3` for Cole–Hopf plateau
controls, `t8` for Burgers controls, or `t9` for two-step heat predictions
and validation. For example:

```bash
python reproduce.py t9
```

Run `python reproduce.py` with no target to display the available commands.

## Run a Modified Case

Copy a JSON file from `configs/`, change its parameters, and pass it to the
solver:

```bash
cp configs/heat_step_dirichlet.json configs/my_heat.json
python main.py configs/my_heat.json
```

Reaction–diffusion and Burgers examples are `fhn_grw_steady.json` and
`burgers_stationary_shock.json` in `configs/`. `config_template.jsonc`
documents the available fields. Plots are saved below `outputs/`.
Save your edited input alongside them.

For a case matching one of the supplied exact references, also compute and
print error metrics with:

```bash
python verify_solver.py --equation heat --config configs/my_heat.json
```

The files in `studies/` are complete examples of parameter sweeps, multi-seed
experiments, error decompositions, and controlled comparisons. Adapt their
parameters, seed lists, and output directories for new studies.

## Repository Layout

- `simulation.py`, `config.py`, `utils.py`: solvers, configuration, and shared helpers.
- `main.py`, `verify_solver.py`: run a case and compare with an exact reference.
- `configs/`, `config_template.jsonc`: editable example inputs.
- `studies/`, `study_paper_refinement.py`: paper experiments and reusable study examples.
- `reproduce.py`, `verify_ensembles.py`, `checks/`: reproduce and check reported results.
- `figure_data/`, `pinned_ensembles/`, `expected_values.json`: reference data for the reported values and figures.
- `figure_scripts/`: figure generation.

## Citation

The original version 1.0.0 is archived on
[Zenodo](https://doi.org/10.5281/zenodo.22050659).
This branch also includes the added two-step heat study and updated checks.
See `CITATION.cff` for citation metadata.

## Acknowledgments

The authors thank Oliver Stalker for providing an early version of the
Python code.

**Principal Investigator:** [Professor Prabir Daripa](https://artsci.tamu.edu/mathematics/contact/profiles/prabir-daripa.html) — Texas A&M University, Department of Mathematics

Other projects from the Daripa Research Group are available on the
[group's GitHub page](https://github.com/Daripa-Research-Group).
