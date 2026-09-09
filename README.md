# Gradient Random Walk Solvers

Code and data for *Error attribution in gradient random walk methods for
parabolic equations* by Stephen Abkin and Prabir Daripa. See the paper for
methods and analysis, and `CITATION.cff` for citation metadata.

## Setup

Use Python 3.11. Clone the revised-paper branch, or extract the supplement
and work inside its `source/` folder:

```bash
git clone --branch grw-solvers-v3 https://github.com/stephen122204/Gradient-Random-Walk-Solvers.git
cd Gradient-Random-Walk-Solvers
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows PowerShell, activate with `.\.venv\Scripts\Activate.ps1` instead.

## Reproduce the paper

Generate all eleven figures from the committed data, then rerun and compare
the numerical studies with the reference results:

```bash
python reproduce.py paper
python reproduce.py verify-all
```

Reproduce the paper’s heat identity controls, seed-block controls, and
bootstrap intervals:

```bash
python checks/heat_cdf_identity_check.py --fresh
python checks/two_step_reference_and_seed_blocks.py
python checks/bootstrap_seed_grouping_check.py --heat-profiles --json output/bootstrap.json
```

The bootstrap output’s `joint` intervals are those used in the manuscript.
Results go under `output/`. The eleven figures are in
`output/final_prepublication_tests/paper_figures/`. Rerunning replaces generated
outputs. Allow several minutes for the studies and several GB of available
memory for the direct-distribution control. Runtime depends on hardware.

## Use your own inputs

Copy a JSON file from `configs/`, edit its physical and numerical parameters,
and run the edited file:

```bash
cp configs/heat_step_dirichlet.json configs/my_heat.json
python main.py configs/my_heat.json
```

Reaction–diffusion and Burgers examples are `fhn_grw_steady.json` and
`burgers_stationary_shock.json` in `configs/`. Fields are described in
`config_template.jsonc`. Plots go under `outputs/`. Save your edited input
alongside them. For a sweep or ensemble, adapt the relevant script in
`studies/`, including its parameters, seed list, and output directory.
