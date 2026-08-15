# Gradient Random Walk Methods for the Heat, FitzHugh-Nagumo, and Viscous Burgers Equations

Code and data for the paper *Controlled Error Attribution for Gradient Random
Walk Methods Applied to the Heat, FitzHugh-Nagumo, and Viscous Burgers'
Equations* by Stephen Abkin and Prabir Daripa. A clean
checkout of this branch reproduces and verifies every number and figure in
the paper using only files in the repository.

## Reproduce and Verify

```bash
git clone https://github.com/stephen122204/heat_burgers_fhn.git
cd heat_burgers_fhn
git checkout grw-solvers-v1
python -m venv .venv && source .venv/bin/activate   # tested with Python 3.11.4
pip install -r requirements.txt

python reproduce.py verify             # single-seed pipeline, 138 checks
python reproduce.py verify --deep      # adds fresh seed-42 array comparisons, 175 checks
python reproduce.py verify-ensembles   # re-runs the five ensemble studies, 18 pinned file comparisons
python -m unittest discover -s tests   # fast solver/control invariant checks
```

`verify` re-runs the single-seed refinement and domain studies and compares
every reported value against `expected_values.json` with relative and
absolute tolerance `1e-12`. That tolerance sits far above cross-platform
floating-point noise and far below the precision of any reported value, so a
genuine change in any paper value cannot pass. Under the pinned environment
in `requirements.txt` the regeneration is bit-identical.

`verify-ensembles` re-runs the multi-seed studies and compares every pinned
numeric field against `pinned_ensembles/`, which holds committed copies of
the study outputs down to individual realizations. Wall-clock fields and
order `1e-16` identity residuals are excluded. Tolerances are relative
`1e-9` and absolute `1e-12`.

## The Studies

The single-seed studies (seed 42) live in `study_paper_refinement.py`. The
ensemble studies live under `studies/` and use fixed seed lists documented in
each script. Every study writes its tables, summaries, and figures under
`output/final_prepublication_tests/`.

```bash
python reproduce.py ensembles   # all five ensemble studies (t4, t7, t5, t3, t8)
python reproduce.py t4          # heat thirty-seed study
python reproduce.py t7          # paired heat output-grid study with the aligned correction
python reproduce.py t5          # scalar FHN thirty-seed study
python reproduce.py t3          # Cole-Hopf plateau evidence
python reproduce.py t8          # Burgers controlled attribution (decoupled, boundary, response, domain)
```

* **t4:** heat bias-spread-total decomposition and its ensemble figure
  (spread exponent near the Monte Carlo reference).
* **t7:** the paired output-grid study. The same thirty realizations per
  particle count are reconstructed on the coupled grid and on fixed 300 and
  400 bin grids. The coupled rows reproduce t4 digit for digit, the built-in
  cross-check that the two studies share one solver and one seed list.
* **t5:** scalar FHN profile, center, speed, and aligned-profile convergence
  with realization-level bootstrap intervals, plus the time-step quartet and
  the deterministic Neumann boundary diagnostic.
* **t3:** the four original Cole-Hopf plateau experiments (domain sensitivity,
  deterministic-transform control, perturbed transformed field, coupled
  particle and output-grid refinement), kept as raw evidence.
* **t8:** the controlled attribution behind the paper's Burgers conclusions.
  A validated parameterized pipeline (bit-identical to the packaged solver at
  the paper configuration) decouples glob count, output grid, and smoothing
  bandwidth over twenty-seed ensembles, adds the deterministic
  boundary-consistent control (pinned versus exact transformed endpoint
  data), the white and kernel-smoothed perturbation response curves, and the
  thirty-seed domain decomposition. Its pinned data include every
  realization-level scalar used in the reported intervals and trends.

## Figures

```bash
python reproduce.py figures     # the eight representative figures (PDF and PNG)
python reproduce.py studies     # re-run the single-seed refinement and domain studies
python reproduce.py all         # studies then figures
python reproduce.py paper1-figures  # the ten figures used by the combined paper
```

Two figure sources, matching the paper:

* **Archived representative arrays** (`figure_data/representative_figure_arrays.npz`,
  seed 42) drive the comparison and diagnostic figures. `reproduce.py figures`
  re-plots the archived data without re-running any simulation. Use
  `--rerun-arrays` to re-simulate the arrays fresh.
* **Checked-in study CSVs** (`figure_data/*.csv`, the data of record) drive
  the refinement and domain-sensitivity figures. `reproduce.py studies`
  re-runs those studies and `verify` confirms the results match the data of
  record.

`paper1-figures` is the canonical combined-paper figure target. It draws four
representative figures from the archived seed-42 arrays and six
ensemble/control figures from committed data in `pinned_ensembles/`, then
writes a SHA-256 source manifest beside the ten figures under
`output/final_prepublication_tests/paper_figures/`. It does not rely on an
earlier study run or select a timestamped output directory.

## Run Your Own Cases

The sections above reproduce the paper. The solvers also run on
user-supplied problems. Copy a JSON config from `configs/`, edit it, and run
one of

```bash
python main.py configs/heat_step_dirichlet.json       # heat GRW
python main.py configs/fhn_grw_steady.json            # scalar FHN GRW
python main.py configs/burgers_stationary_shock.json  # Cole-Hopf Burgers
python run_config.py                                  # same, with a file picker
```

Every config field is documented with comments in `config_template.jsonc`.
Comparison figures are saved under `outputs/<timestamp>/`. Where an exact
solution exists, a custom run can be checked against it with
`python verify_solver.py --equation heat --config <your_config>.json`. These
tools are for exploration. They are not the source of any number in the
paper.

## Repository Layout

* **Solvers and studies:** `simulation.py`, `config.py`, `verify_solver.py`
  (exact solutions and error metrics), `study_paper_refinement.py`,
  `studies/`.
* **Entry points:** `reproduce.py` (all targets above) and
  `verify_ensembles.py`.
* **Interactive tools:** `main.py`, `run_config.py`, `utils.py`, `configs/`,
  and `config_template.jsonc` (own-case exploration, not paper inputs).
* **Data of record:** `figure_data/`, `expected_values.json`, and
  `pinned_ensembles/`.
* **Release checks:** `tests/` and `PROVENANCE.md`.

## Citation

See `CITATION.cff`.
