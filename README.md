# GRW Solvers — Heat, FitzHugh–Nagumo, Burgers (companion code)

Companion code and data for:

> Stephen Abkin and Prabir Daripa, *A Numerical Study of Gradient Random Walk
> Methods for Heat, FitzHugh–Nagumo, and Burgers' Equations*, 2026.
> arXiv id to be assigned.

Every number, table, and figure in the paper regenerates from this branch
with one command. All paper computations use NumPy seed 42, re-seeded at the
start of each configuration, so runs are identically reproducible.

The method follows the Gradient Random Walk approach.

<!-- TODO(license): choose and add a LICENSE file (MIT or BSD-3 recommended
     for research code) before making this public. See RELEASE_NOTES.md. -->

---

## Quickstart

```bash
python -m venv .venv
source .venv/bin/activate        # macOS / Linux
pip install -r requirements.txt

python reproduce.py all          # every study + every figure   (~16 s)
python reproduce.py verify       # PASS/FAIL vs the paper's values (~12 s)
```

`MANIFEST.md` maps each table and figure of the paper to its command, output
file, and measured runtime. `expected_values.json` pins every reported value;
`python reproduce.py verify --deep` additionally re-runs the representative
simulations and requires bit-identical arrays against the archived
`figure_data/representative_figure_arrays.npz`.

---

## Important: paper configurations vs. exploration tools

The paper's pinned configurations live in `study_paper_refinement.py` and
`figure_scripts/regenerate_manuscript_figures.py`, which `reproduce.py`
invokes. **`main.py` and the JSON files under `configs/` are interactive
exploration tools; their values are starting points, not the source of the
paper's numbers.** Two details a reader should know (both stated in the
paper, Sec. 6.1–6.2):

- the representative heat run (Figure 1, Table 2) reconstructs on a 400-bin
  output grid, while the heat refinement study (Figure 2) uses 300 bins;
- the FHN profile/front-tracking runs use homogeneous Neumann boundaries,
  while the FHN refinement study uses Dirichlet boundaries.

---

## Directory layout

| Path | Contents |
|---|---|
| `simulation.py` | GRW solvers: heat, scalar FHN (reaction-weight update), Cole–Hopf Burgers |
| `config.py`, `configs/`, `main.py`, `run_config.py` | JSON-config exploration interface |
| `study_paper_refinement.py` | The three paper studies: heat N-sweep, FHN N-sweep, Burgers domain sensitivity (seed 42) |
| `figure_scripts/` | Manuscript figure pipeline (`regenerate_manuscript_figures.py`, `plot_style.py`) |
| `figure_data/` | Data of record: study CSVs + archived representative arrays (seed 42) |
| `reproduce.py` | One command per paper artifact; `verify` target |
| `expected_values.json` | Every reported value, pinned for `verify` |
| `MANIFEST.md` | Table/figure → command → output → runtime map |
| `verify_solver.py`, `verify_grw.py` | Standalone accuracy checks vs exact solutions |
| `make_final_plots_original_grw.py`, `make_publication_plots_original_grw.py` | Supplementary plotting used during drafting |

Generated outputs go to `output/` and `outputs/<timestamp>/` (gitignored).

---

## Run a single simulation (exploration interface)

```bash
python main.py configs/heat_step_dirichlet.json
```

The figure is saved to `outputs/YYYY-MM-DD_HH-MM-SS/`. Pass any config file
from `configs/` or one you wrote yourself.

All eight ready-to-run configs:

```bash
python main.py configs/heat_step_dirichlet.json
python main.py configs/heat_step_neumann.json
python main.py configs/burgers_stationary_shock.json
python main.py configs/burgers_shock.json
python main.py configs/burgers_traveling_wave.json
python main.py configs/fhn_grw_steady.json
python main.py configs/fhn_grw_nonsmooth.json
python main.py configs/fhn_grw_discontinuous.json
```

---

## Config file guide

Each JSON config maps directly to the parameters used in the paper. Copy one
of the files in `configs/` and edit it.

### Fields common to all equations

| JSON field | Paper symbol | Description |
|---|---|---|
| `equation_type` | — | `"heat"`, `"burgers"`, or `"fitzhugh-nagumo"` |
| `domain_type` | — | `"Finite"` (all paper benchmarks use this) |
| `domain_size` | L | Right endpoint; domain is [0, L] |
| `diff_constant` | α / ν / D | Diffusivity for heat, viscosity for Burgers, D for FHN |
| `time_step` | Δt | Time step size |
| `total_time` | T | Final simulation time |
| `num_points` | N | Number of GRW globs (particles) |
| `boundary_conditions` | — | Left and right BC types (see below) |
| `reaction_term` | — | Set `false` (reserved field, not used) |

**Boundary conditions:**

```json
"boundary_conditions": {
  "LEFT":  { "type": "Dirichlet", "value": 0.0 },
  "RIGHT": { "type": "Dirichlet", "value": 1.0 }
}
```

- `"Dirichlet"` — fixes the solution value at the wall; globs that cross the
  boundary are reflected symmetrically (weight preserved).
- `"Neumann"` — zero-flux wall (`"value"` must be `0.0`); globs are reflected
  anti-symmetrically (weight negated).

### Heat equation — u_t = α u_xx

Paper benchmark: α = 0.1, L = 10, x₀ = 5, T = 0.5, Δt = 0.001, N = 50 000.
Exact solution is the error-function profile.

With Dirichlet–Dirichlet BCs and a step IC, the output overlays the exact
error-function solution (Eq. (4.3) in the paper). Increasing N reduces
stochastic noise; the paper uses N = 50 000.

### Burgers equation — u_t + u u_x = ν u_xx (Cole–Hopf GRW)

Paper benchmark: ν = 0.5, A = 1, x_c = 2, L = 4, T = 0.5, Δt = 0.005,
N = 400. Exact solution is the stationary shock u = −A tanh(A(x−x_c)/(2ν)).

The solver applies the Cole–Hopf transformation u = −2ν φ_x/φ to convert
Burgers into a heat equation for φ, runs GRW on φ, then reconstructs u. This
avoids the ill-conditioned direct GRW formulation discussed in the paper
(Sec. 4.3). The paper's diagnostic separates the total error into E_det
(finite-domain/deterministic pipeline) and E_GRW (particle reconstruction);
at L = 4 with N = 400, E_det = 0.172 > E_GRW = 0.081, and the two move in
opposite directions as L changes — see Table 1 of the paper.

`nu` inside `burgers_initial_condition` must match `diff_constant`.

### FitzHugh–Nagumo equation — u_t = D u_xx + f(u)

Paper benchmark: D = 0.5, a = 0.25, L = 30, x_c = 15, T = 9, Δt = 0.01,
N = 500. Wave speed θ = √2(0.5 − a) ≈ 0.354. Exact solution is the traveling
wave u(x,t) = 1/(1 + exp(−(x + θt − x_c)/2)).

The solver updates glob weights at each step using R(u) = f′(u) derived from
the exact wave (paper Secs. 3.5 and 4.2). This is the scalar traveling-wave
reduction used in the paper — not the full two-variable FHN system. Domain
L = 30 is chosen so the front never reaches the walls over the full run.

**`fhn_initial_condition` options:** `"steady_solution"` (logistic quantile
positions — the paper IC), `"nonsmooth"` (linear ramp), `"discontinuous"`
(Heaviside).

---

## Standalone verification against exact solutions

```bash
python verify_solver.py --equation all
python verify_solver.py --equation heat
python verify_solver.py --equation burgers --config configs/burgers_stationary_shock.json
python verify_solver.py --equation fhn
python verify_grw.py                     # heat-only glob statistics checks
```

Outputs (in `output/verify/<equation>/`): `comparison_plot.png`,
`metrics.json`, and for Burgers `cole_hopf_diagnostics.png`. These runs are
stochastic (no fixed seed) — they reproduce the experiments, not the exact
paper numbers; use `reproduce.py` for the paper's pinned values.

---

## Measured runtimes (Apple Silicon, Python 3.11.4)

| Command | Time |
|---|---|
| `python reproduce.py studies` | ~11 s |
| `python reproduce.py figures` | ~5 s |
| `python reproduce.py verify` | ~12 s |
| `python reproduce.py verify --deep` | ~13 s |
