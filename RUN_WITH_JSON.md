# Running With JSON Configs

This repo is config-first: edit a JSON file in `configs/`, then run one command.

## 1. Edit Parameters

Start with one of these files:

```text
configs/heat_step_dirichlet.json
configs/heat_step_neumann.json
configs/fhn_grw_steady.json
configs/fhn_grw_nonsmooth.json
configs/fhn_grw_discontinuous.json
configs/burgers_stationary_shock.json
configs/burgers_shock.json
configs/burgers_traveling_wave.json
```

Common parameters to change:

```json
"domain_size": 10.0,
"diff_constant": 0.1,
"time_step": 0.001,
"total_time": 0.5,
"num_points": 50000
```

Equation-specific blocks are:

```json
"heat_initial_condition": { ... }
"fhn_initial_condition": { ... }
"burgers_initial_condition": { ... }
```

## 2. Easiest Launcher

Run:

```bash
python run_config.py
```

If your system supports it, a file picker opens. Choose a JSON file from
`configs/`. If a file picker is not available, the terminal asks for the path.

The launcher then asks for a mode:

```text
verify
simulate
```

Use `verify` when you want comparison plots, exact/reference overlays, and
metrics. This is usually the best choice for paper-style evidence.

Use `simulate` when you only want the raw solver output plot.

## 3. Direct Commands

Verification/comparison plots:

```bash
python run_config.py configs/heat_step_dirichlet.json --mode verify
python run_config.py configs/fhn_grw_steady.json --mode verify
python run_config.py configs/burgers_stationary_shock.json --mode verify
```

Raw solver plots:

```bash
python run_config.py configs/heat_step_dirichlet.json --mode simulate
python run_config.py configs/fhn_grw_steady.json --mode simulate
python run_config.py configs/burgers_stationary_shock.json --mode simulate
```

## 4. Output Locations

Verification outputs:

```text
output/verify/heat/comparison_plot.png
output/verify/heat/metrics.json

output/verify/fhn/comparison_plot.png
output/verify/fhn/fhn_grw_diagnostics.png
output/verify/fhn/metrics.json

output/verify/burgers/comparison_plot.png
output/verify/burgers/cole_hopf_diagnostics.png
output/verify/burgers/metrics.json
```

Raw simulation outputs:

```text
output/heat_density.png
output/heat_field_dirichlet.png
output/fhn_uv.png
output/burgers_u.png
```

## 5. Paper-Style Figures

To generate the paper-style figures from the manuscript plotting pipeline:

```bash
MPLCONFIGDIR=/tmp/grw-figure-mpl \
XDG_CACHE_HOME=/tmp/grw-figure-cache \
python3 figure_scripts/regenerate_manuscript_figures.py
```

Outputs:

```text
_DRAFT__GRW_Methods (1)/figures/
```

That paper-style generator uses fixed seed `42`, archived representative arrays
in `figure_data/representative_figure_arrays.npz`, and CSV summaries in
`output/paper_refinement_original_grw/`.
