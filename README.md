# GRW Solvers — Heat, FitzHugh–Nagumo, Burgers

Companion code for *A Numerical Study of Gradient Random Walk Methods for
Heat, FitzHugh–Nagumo, and Burgers' Equations* (Abkin & Daripa).
Paper link: to be added.

Gradient random walk (GRW) particle solvers for three benchmark PDEs: the
heat equation (direct GRW), a scalar FitzHugh–Nagumo traveling front (GRW
with reaction-driven weight updates), and viscous Burgers' equation through
the Cole–Hopf transformation. Every number, table, and figure in the paper
regenerates from this branch with one command (fixed seed 42).

## How to run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt   # numpy, matplotlib

python reproduce.py all           # every study + every figure   (~16 s)
python reproduce.py verify        # PASS/FAIL vs the paper's values (~12 s)
python reproduce.py verify --deep # + bit-identical array regeneration

# Single exploratory simulations from JSON configs:
python main.py configs/heat_step_dirichlet.json
```

Note: `main.py` and the JSON files under `configs/` are exploration tools;
the paper's pinned configurations live in `study_paper_refinement.py` and
`figure_scripts/regenerate_manuscript_figures.py`, which `reproduce.py`
invokes. All paper computations use NumPy seed 42, re-seeded at the start
of each configuration.

## What was used in the paper

| Paper artifact | Command | Data |
|---|---|---|
| Figures 1, 3, 4, 6, 7 + Table 2 | `python reproduce.py figures` | archived arrays in `figure_data/representative_figure_arrays.npz` (seed 42) |
| Figure 2 (heat fixed-grid diagnostic) | `python reproduce.py studies` | `figure_data/heat_refinement_summary.csv` |
| Figure 5 (FHN refinement) | `python reproduce.py studies` | `figure_data/fhn_refinement_summary.csv` |
| Figure 8 + Table 1 (Burgers domain study) | `python reproduce.py studies` | `figure_data/burgers_domain_sensitivity_summary.csv` |

Two grid details, both stated in the paper: the representative heat run
(Figure 1, Table 2) reconstructs on a 400-bin output grid while the heat
refinement study (Figure 2) uses 300 bins, and the FHN profile/front runs
use Neumann boundaries while the FHN refinement study uses Dirichlet.
`expected_values.json` pins the values `reproduce.py verify` checks.
