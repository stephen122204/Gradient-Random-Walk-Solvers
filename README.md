# [GRW Solvers] Heat, Burgers, FitzHugh-Nagumo

Gradient Random Walk (GRW) solver for three benchmark PDEs. Each run reads a JSON config file, runs the solver, and saves comparison figures to `outputs/<timestamp>/`.

---

## Setup

```bash
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
pip install numpy matplotlib
```

---

## Run a single simulation

```bash
python main.py configs/heat_step_dirichlet.json
```

The figure is saved to `outputs/YYYY-MM-DD_HH-MM-SS/`. Pass any config file from `configs/` or one you wrote yourself.

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

## Recreate all paper figures

```bash
python figure_scripts/regenerate_manuscript_figures.py
```

This regenerates all eight figures from the paper using a fixed random seed (SEED = 42) and saves them to `outputs/<timestamp>/`. It uses cached simulation arrays from `figure_data/` by default. To force a fresh simulation run:

```bash
python figure_scripts/regenerate_manuscript_figures.py --rerun
```

Terminal output shows the figure filenames and the verification error table (L2, Linf, rel L2) matching the paper's Table 1.

---

## Config file guide

Each JSON config maps directly to the parameters used in the paper. Copy one of the files in `configs/` and edit it.

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

- `"Dirichlet"` — fixes the solution value at the wall; globs that cross the boundary are reflected symmetrically (weight preserved).
- `"Neumann"` — zero-flux wall (`"value"` must be `0.0`); globs are reflected anti-symmetrically (weight negated).

---

### Heat equation — u_t = α u_xx

Paper benchmark: α = 0.1, L = 10, x₀ = 5, T = 0.5, Δt = 0.001, N = 50 000. Exact solution is the error-function profile.

```json
{
  "equation_type": "heat",
  "domain_type": "Finite",
  "domain_size": 10.0,
  "boundary_conditions": {
    "LEFT":  { "type": "Dirichlet", "value": 0.0 },
    "RIGHT": { "type": "Dirichlet", "value": 1.0 }
  },
  "diff_constant": 0.1,
  "time_step": 0.001,
  "total_time": 0.5,
  "num_points": 50000,
  "reaction_term": false,
  "heat_initial_condition": {
    "type": "step",
    "jump_position": 5.0,
    "jump_height": 1.0
  }
}
```

**`heat_initial_condition` options:**

| `type` | Description | Extra fields |
|---|---|---|
| `"step"` | All globs start at `jump_position` | `jump_position`, `jump_height` |
| `"uniform_gradient"` | Globs spread uniformly across [0, L] | `gradient_value` |
| `"gaussian_cloud"` | Globs drawn from a Gaussian | `center`, `sigma`, `jump_height` |

With Dirichlet–Dirichlet BCs and a step IC, the output overlays the exact error-function solution (equation 7 in the paper). Increasing N reduces stochastic noise; the paper uses N = 50 000.

---

### Burgers equation — u_t + u u_x = ν u_xx (Cole-Hopf GRW)

Paper benchmark: ν = 0.5, A = 1, x_c = 2, L = 4, T = 0.5, Δt = 0.005, N = 400. Exact solution is the stationary shock u = −A tanh(A(x−x_c)/(2ν)).

The solver applies the Cole-Hopf transformation u = −2ν φ_x/φ to convert Burgers into a heat equation for φ, runs GRW on φ, then reconstructs u. This avoids the numerically unstable direct GRW formulation discussed in the paper (Section 3.3).

```json
{
  "equation_type": "burgers",
  "domain_type": "Finite",
  "domain_size": 4.0,
  "boundary_conditions": {
    "LEFT":  { "type": "Dirichlet", "value": 0.0 },
    "RIGHT": { "type": "Dirichlet", "value": 0.0 }
  },
  "diff_constant": 0.5,
  "time_step": 0.005,
  "total_time": 0.5,
  "num_points": 400,
  "reaction_term": false,
  "burgers_mode": "cole_hopf_grw",
  "burgers_initial_condition": {
    "type": "stationary_shock",
    "nu": 0.5,
    "x_center": 2.0,
    "amplitude": 1.0
  }
}
```

**`burgers_initial_condition` options:**

| `type` | Description | Extra fields |
|---|---|---|
| `"stationary_shock"` | Exact stationary shock −A tanh(A(x−x_c)/(2ν)) | `nu`, `x_center`, `amplitude` |
| `"traveling_wave"` | Exact traveling wave at unit speed | `nu`, `x_center` |
| `"shock"` | Step profile | `shock_position`, `left_value`, `right_value` |

`nu` inside `burgers_initial_condition` must match `diff_constant`. The paper diagnostic separates the total error into E_det (finite-domain boundary mismatch) and E_GRW (GRW particle error). At L = 4 with N = 400, the paper finds E_det = 0.172 > E_GRW = 0.081. To reduce E_det, increase `domain_size`; to reduce E_GRW, increase `num_points` (these two errors move in opposite directions as L changes — see Table 2 of the paper).

---

### FitzHugh-Nagumo equation — u_t = D u_xx + f(u)

Paper benchmark: D = 0.5, a = 0.25, L = 30, x_c = 15, T = 9, Δt = 0.01, N = 500. Wave speed θ = √2(0.5 − a) ≈ 0.354. Exact solution is the traveling wave u(x,t) = 1/(1 + exp(−(x + θt − x_c)/2)).

The solver updates glob weights at each step using the reaction statistic R(u) = f′(u) derived from the exact wave. This is the scalar traveling-wave reduction used in the paper — not the full two-variable FHN system.

```json
{
  "equation_type": "fitzhugh-nagumo",
  "domain_type": "Finite",
  "domain_size": 30.0,
  "boundary_conditions": {
    "LEFT":  { "type": "Neumann", "value": 0.0 },
    "RIGHT": { "type": "Neumann", "value": 0.0 }
  },
  "diff_constant": 0.5,
  "time_step": 0.01,
  "total_time": 9.0,
  "num_points": 500,
  "reaction_term": false,
  "a": 0.25,
  "b": 1.0,
  "tau": 1.0,
  "fhn_initial_condition": {
    "type": "steady_solution",
    "a": 0.25,
    "x_center": 15.0
  }
}
```

**FHN-specific fields:**

| JSON field | Paper symbol | Description |
|---|---|---|
| `a` | a | Threshold; controls wave speed θ = √2(0.5 − a) |
| `b` | b | Recovery coupling (set 1.0 for scalar benchmark; not used) |
| `tau` | τ | Temporal scaling (set 1.0 for scalar benchmark; not used) |

The `a` value inside `fhn_initial_condition` must match the outer `a` field. `x_center` is x_c, the initial location of the u = 0.5 level.

**`fhn_initial_condition` options:**

| `type` | Description | Extra fields |
|---|---|---|
| `"steady_solution"` | Globs at exact logistic quantile positions (paper IC) | `a`, `x_center` |
| `"nonsmooth"` | Linear-ramp IC — C0 but not C1 | `a`, `x_center`, `half_width` |
| `"discontinuous"` | Heaviside IC — single Dirac-delta glob at x_c | `a`, `x_center` |

The output is a 2×2 panel comparing GRW to the exact traveling wave at t = 0, T/3, 2T/3, T. All three IC types converge to the same traveling wave by T = 9 for the paper parameters. Domain L = 30 is chosen so the front never reaches the walls over the full run.

---

## Verify against exact solutions

```bash
python verify_solver.py --equation heat
python verify_solver.py --equation burgers --config configs/burgers_stationary_shock.json
python verify_solver.py --equation fhn
```

Saves `comparison_plot.png`, `metrics.json`, and (for Burgers) `cole_hopf_diagnostics.png` to `output/verify/<equation>/`.
