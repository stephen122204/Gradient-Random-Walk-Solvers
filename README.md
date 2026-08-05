# Gradient Random Walk Methods: Heat, FitzHugh–Nagumo, and Burgers

This repository contains the code and data accompanying the paper *A Numerical
Study of Gradient Random Walk Methods for Heat, FitzHugh–Nagumo, and Burgers'
Equations* by Stephen Abkin and Prabir Daripa. A clean checkout of this branch
independently verifies every number reported in the paper and regenerates all
eight paper figures, using only files in the repository.

**Paper:** link to be added (arXiv preprint forthcoming).

## Getting Started

```bash
git clone https://github.com/stephen122204/heat_burgers_fhn.git
cd heat_burgers_fhn
git checkout grw-solvers-v1
python -m venv .venv && source .venv/bin/activate   # tested with Python 3.11
pip install -r requirements.txt
```

## Verifying the Paper's Numbers

Every computation in the paper uses a single fixed seed (42). The paper's
design is deliberately single-seed: it reports no fitted convergence
exponents, no confidence intervals, and no convergence theorems, and the
verification below checks exactly the values the paper reports, nothing more.

```bash
python reproduce.py verify        # re-runs the studies, checks every reported value
python reproduce.py verify --deep # also re-runs the three representative simulations
```

`verify` compares against `expected_values.json` using documented tolerances
(relative/absolute `1e-12`): far above cross-platform floating-point noise,
far below the precision of any value reported in the paper, so a genuine
change in any paper value cannot pass. Integers, row counts, and seeds are
compared exactly. Under the pinned environment in `requirements.txt` the
regeneration is bit-for-bit identical; under other NumPy builds, last-digit
noise up to about `2.4e-14` relative has been observed, which the tolerance
absorbs. The summary line prints the actual check count (138 checks, or 175
with `--deep`).

## Regenerating the Paper's Figures

```bash
python reproduce.py figures       # all 8 figures (PDF + PNG) -> outputs/<timestamp>/
python reproduce.py studies       # re-run the refinement/domain studies from scratch
python reproduce.py all           # studies + figures
```

Two kinds of figure sources, matching the paper:

- **Archived representative arrays** (`figure_data/representative_figure_arrays.npz`,
  seed 42) drive Figures 1, 3, 4, 6, 7 and Table 2. `reproduce.py figures`
  re-plots this archived data without rerunning the simulations; use
  `python reproduce.py figures --rerun-arrays` to re-simulate them fresh.
- **Checked-in study CSVs** (`figure_data/*.csv`, the paper's data of record)
  drive Figures 2, 5, 8 and Table 1. `reproduce.py studies` re-runs those
  studies and `verify` confirms the results match the data of record.

The eight figures:

1. **Heat comparison** — GRW reconstruction versus the exact error-function
   profile, with an inset showing the residual `u_N - u_ex` concentrated at
   the transition region.
2. **Heat fixed-grid diagnostic** — `L2_h` error versus particle count `N` on
   a fixed 300-bin grid (the grid-controlled error floor; not a convergence
   rate).
3. **FHN comparison** — traveling front versus the exact wave at
   `t = 0, 3, 6, 9`.
4. **FHN diagnostics** — front location, front-location error, signed glob
   weights (with a zero line making the two Neumann-negated weights near
   `x ~ 26` visible), and the reconstructed fields.
5. **FHN refinement** — profile and front-location error versus `N`, with a
   black dashed `O(N^{-1/2})` reference guide (a visual guide, not a fitted
   exponent).
6. **Burgers comparison** — recovered field versus the exact stationary
   shock.
7. **Burgers diagnostics** — panels (a)–(d); panel (d) shows the pointwise
   components `u_FD - u_ex` and `u_GRW - u_FD`, whose RMSEs are `E_det` and
   `E_GRW`.
8. **Burgers domain sensitivity** — `E_det`, `E_GRW`, and `E_total` versus
   domain size `L`. The study uses `N = 100L` (fixed particle density and
   grid spacing), so it does not isolate the particle count. The decline of
   `E_det` with `L` follows the whole-domain RMSE normalization
   (`E_det * sqrt(L)` is nearly constant — RMSE dilution), as discussed in
   the paper.

Terminology follows the paper: `E_det` is the deterministic component
(deterministic-pipeline mismatch), `E_GRW` the particle-reconstruction
component, `E_total` the total; the components are not orthogonal and need
not sum to the total.

Each figure run writes `figure_regeneration_metadata.json` alongside the
outputs, recording the seed and the source of every figure with
repository-relative paths. The pinned configurations live in
`study_paper_refinement.py` and `figure_scripts/regenerate_paper_figures.py`,
which `reproduce.py` invokes.

## Running Your Own Experiments

The sections above reproduce the paper. To run the solvers with your own
inputs, copy a JSON config from `configs/`, edit it, and run:

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
