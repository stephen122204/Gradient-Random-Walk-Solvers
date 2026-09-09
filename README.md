# Gradient Random Walk Solvers

Code and data for *Error attribution in gradient random walk methods for
parabolic equations* by Stephen Abkin and Prabir Daripa. The examples cover
heat, scalar reaction–diffusion with a FitzHugh–Nagumo-type cubic term, and
Burgers’ equation.

## Setup

Clone the revised-paper branch:

```bash
git clone --branch grw-solvers-v3 https://github.com/stephen122204/Gradient-Random-Walk-Solvers.git
cd Gradient-Random-Walk-Solvers
```

For the ZIP supplement, extract it and work inside `source/` instead.
Use Python 3.11 and install the pinned dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, replace the activation command with
`.\.venv\Scripts\Activate.ps1` in PowerShell or
`.venv\Scripts\activate.bat` in Command Prompt.
The reference environment is Python 3.11.4.

## Reproduce the paper

Generate all eleven paper figures from the committed data:

```bash
python reproduce.py paper
```

Rerun the representative calculations, five original ensemble/control studies,
and two-step heat prediction study, then compare their results with the
committed reference values:

```bash
python reproduce.py verify-all
```

The paper’s analytical controls and bootstrap comparison table use three
additional commands:

```bash
python checks/heat_cdf_identity_check.py --fresh
python checks/two_step_reference_and_seed_blocks.py
python checks/bootstrap_seed_grouping_check.py --heat-profiles --json output/bootstrap.json
```

These reproduce the heat error identities and seed-block checks, the two-step
direct-distribution and finite-ensemble controls, and the fitted-rate intervals,
respectively. The bootstrap output’s `joint` field gives the manuscript’s
primary intervals. Its `reported` field retains the earlier independent
intervals for comparison. The direct-distribution control uses 20 million
samples and needs several GB of available memory.

To rerun individual parts:

| Paper calculation | Command |
|---|---|
| Representative results and saved particle arrays | `python reproduce.py verify --deep` |
| Heat ensemble | `python reproduce.py t4` |
| Paired heat reconstructions | `python reproduce.py t7` |
| Scalar reaction–diffusion ensemble | `python reproduce.py t5` |
| Cole–Hopf plateau controls | `python reproduce.py t3` |
| Burgers refinement, boundary, and perturbation controls | `python reproduce.py t8` |
| Two-step heat predictions, production, and validation | `python reproduce.py verify-two-step` |

The two-step command calculates predictions and selects the validation particle
count before running the ensembles. For a quick check of its committed
predictions and error statistics, use
`python reproduce.py verify-two-step --pinned-only`.
To plot and summarize a fresh two-step run:

```bash
python studies/analyze_t9_heat_two_step.py --data-dir output/heat_two_step_reproduction --output-dir output/heat_two_step_reproduction
```

## Outputs and running time

| Output | Directory |
|---|---|
| Eleven paper figures | `output/final_prepublication_tests/paper_figures/` |
| Representative refinement tables | `output/paper_refinement_original_grw/` |
| Original ensemble/control studies | Separate study folders in `output/final_prepublication_tests/` |
| Two-step predictions and ensembles | `output/heat_two_step_reproduction/` |
| Custom-case plots | `outputs/YYYY-MM-DD_HH-MM-SS/` |

On an Apple-Silicon laptop, representative verification takes about 20 seconds,
heat ensemble studies take about 1–2 minutes each, and `verify-all` takes several
minutes. Figure generation takes under a minute. Times depend on hardware.

Rerunning a study replaces its generated files. Committed reference data in
`figure_data/`, `pinned_ensembles/`, and `expected_values.json` are retained.
After an interruption, rerun the affected study. To compare completed ensemble
outputs without rerunning, use `python reproduce.py verify-ensembles --no-rerun`
or `python reproduce.py verify-two-step --no-rerun`.

Paper studies specify their parameters and seed lists in the study scripts.
Paired reconstructions share particle trajectories, and ensemble refinement
levels reuse seed lists. The original two-step design timestamp and run log
are retained in `provenance/heat_two_step/`.

## Run your own case

Copy an example from `configs/`, edit its parameters, and run it:

```bash
cp configs/heat_step_dirichlet.json configs/my_heat.json
python main.py configs/my_heat.json
```

`config_template.jsonc` describes the fields. Reaction–diffusion and Burgers
examples are `configs/fhn_grw_steady.json` and
`configs/burgers_stationary_shock.json`. Save the edited configuration with
your output. Custom runs use new random draws; set NumPy’s seed before calling
a solver in a script if you want repeatable results.

For a parameter sweep or ensemble, copy the relevant script from `studies/`
and change its parameters and output directory. The two-step study requires
an even particle count of at least two. For another scalar reaction law,
`config.reaction_derivative` accepts a function acting on a NumPy array;
adapt the initial conditions, boundary data, and reference to your problem.

## Citation and packaging

Citation metadata is in `CITATION.cff`. The original
[Zenodo 1.0.0 release](https://doi.org/10.5281/zenodo.22050659) predates the
revised paper’s two-step extension. Use this branch or the updated supplement
for the revision until its archive version is deposited.

To make a ZIP of the source and committed data:

```bash
python tools/build_code_supplement.py --output output/code-supplement.zip
```

The authors thank Oliver Stalker for providing an early version of the code.
