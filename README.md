# Gradient Random Walk Methods: Heat, FitzHugh–Nagumo, and Burgers

This repository contains computer code for reproducing the numerical results
described in the manuscript *A Numerical Study of Gradient Random Walk
Methods for Heat, FitzHugh–Nagumo, and Burgers' Equations* by Stephen Abkin
and Prabir Daripa.

**Paper:** link to be added (arXiv preprint forthcoming).

## Getting Started

```bash
git clone https://github.com/stephen122204/heat_burgers_fhn.git
cd heat_burgers_fhn
git checkout grw-solvers-v1
```

## Reproducing Numerical Results

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Regenerate every study and figure in the paper (fixed seed 42, ~16 s), and
check the results against the paper's reported values:

```bash
python reproduce.py all       # all studies + all 8 figures
python reproduce.py verify    # PASS/FAIL against the paper's values
```

Figures 1, 3, 4, 6, 7 and Table 2 come from the archived representative
arrays in `figure_data/`; Figures 2, 5, 8 and Table 1 come from the
refinement and domain studies (`reproduce.py studies`). The paper's pinned
configurations live in `study_paper_refinement.py` and
`figure_scripts/regenerate_manuscript_figures.py`, which `reproduce.py`
invokes.

## Running Your Own Experiments

The section above is for convenience in reproducing the paper. To run the
solvers with your own inputs, copy a JSON config from `configs/`, edit it,
and run:

```bash
python main.py configs/heat_step_dirichlet.json       # heat GRW
python main.py configs/fhn_grw_steady.json            # scalar FHN GRW
python main.py configs/burgers_stationary_shock.json  # Cole--Hopf Burgers
```

The main fields are `diff_constant` (α, D, or ν), `time_step`, `total_time`,
`num_points` (particle count), `domain_size`, and `boundary_conditions`
(Dirichlet reflects particles and keeps their weight; Neumann reflects and
negates it). Each equation adds an initial-condition block: heat takes a
step, uniform-gradient, or Gaussian-cloud profile; FHN takes the logistic
front (`steady_solution`), a linear ramp, or a Heaviside step; Burgers takes
a stationary shock, traveling wave, or step, solved through the Cole--Hopf
transformation. `config_template.jsonc` documents every field with comments;
comparison figures are saved under `outputs/<timestamp>/`. To check a custom
run against an exact solution where one exists, use
`python verify_solver.py --equation heat --config <your_config>.json`.

## Citation

See `CITATION.cff`.
