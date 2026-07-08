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

Every config field is documented with comments in `config_template.jsonc`.
Comparison figures are saved under `outputs/<timestamp>/`, and a custom run
can be checked against an exact solution where one exists via
`python verify_solver.py --equation heat --config <your_config>.json`.

## Citation

See `CITATION.cff`.
