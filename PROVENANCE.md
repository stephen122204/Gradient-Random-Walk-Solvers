# Paper 1 reproduction and provenance

This branch is the code and data release for *Controlled Error Attribution for
Gradient Random Walk Methods Applied to the Heat, FitzHugh–Nagumo, and Viscous
Burgers' Equations*. The paper uses representative calculations to show the
spatial error, multi-seed ensembles to measure sampling behavior, and paired
or deterministic controls to identify the responsible numerical operation.

## Release gates

Run from the repository root in the environment pinned by `requirements.txt`:

```bash
python reproduce.py verify --deep
python reproduce.py verify-ensembles
python -m unittest discover -s tests -v
python reproduce.py paper1-figures
```

The first command re-runs the representative studies and compares 175 scalar
and array quantities. The second re-runs the five ensemble/control studies
and compares 18 committed CSV/JSON artifacts. The unit tests check nine
solver and experimental invariants. The final command redraws the ten paper
figures from committed data and writes `paper1_figure_manifest.json` with a
source and SHA-256 digest for each PDF.

## Study map

| Target | Design and paper role | Committed data |
|---|---|---|
| `verify` / `studies` | Seed-42 heat, FHN, and Cole–Hopf representative/refinement calculations | `expected_values.json`, `figure_data/` |
| `t4` | Heat, 30 realizations at each of seven particle counts; bias–spread–total decomposition | `pinned_ensembles/heat_extended/` |
| `t7` | Heat, the identical trajectory for a given `(N, seed)` reconstructed under five grid/evaluation treatments; deterministic reflected-law operator control | `pinned_ensembles/heat_grid_paired/` |
| `t5` | FHN, 30 realizations at each of six particle counts; profile, center, speed, and aligned-profile errors; deterministic boundary diagnostic | `pinned_ensembles/fhn_extended/` |
| `t3` | Original Cole–Hopf plateau evidence retained without automatic causal classification | `pinned_ensembles/cole_hopf_plateau/` |
| `t8` | Burgers attribution: independent variation of initialization density, output grid, and physical bandwidth; exact-boundary and inversion controls; perturbation response; 30-realization domain study | `pinned_ensembles/burgers_controls/` |

The repeated seed identifiers at different particle counts are for
reproducibility; they are not treated as a strict common-random-number
coupling because the random-vector dimensions change. The heat `t7` arms are
genuinely paired because all five treatments reuse the same final particle
set for each `(N, seed)`.

For the Burgers controls, `P` is the number of initialization-grid points and
the run evolves `P-1` transformed-gradient globs; `M` is the number of output
points. `decoupled_per_run.csv` and `domain_per_run.csv` retain the scalar
result for every realization used in the reported slope intervals and domain
trend. The latter trend is a realization-level bootstrap interval for the
linear slope of ensemble-mean total RMSE versus domain size; its interval
contains zero.

## Paper figures

`python reproduce.py paper1-figures` produces:

- representative arrays: `heat_comparison`, `fhn_comparison`,
  `fhn_diagnostics`, `burgers_diagnostics`;
- ensemble/control data: `heat_bias_spread_total_vs_N`, `heat_grid_paired`,
  `fhn_convergence`, `burgers_decoupled`, `burgers_boundary_domain`, and
  `burgers_perturbation_response`.

The representative archive is regenerated with
`python reproduce.py figures --rerun-arrays`. The canonical paper-figure
target deliberately plots the ensemble figures from the committed data of
record, so figure reproduction does not require rerunning the simulations.

## Numerical comparison policy

The representative verifier uses relative and absolute tolerances of
`1e-12`. The ensemble verifier uses relative `1e-9` and absolute `1e-12`.
Machine-dependent timing fields and floating-point identity residuals of
order `1e-16` are excluded. Integer quantities, seed identifiers, row counts,
and strings are exact. These tolerances are far below the printed precision
of the paper and above the observed cross-platform last-digit variation.

The checked-in output is evidence for the tested problems, parameter ranges,
and implementations. It does not substitute for a convergence theorem, and
the scripts do not encode causal conclusions by threshold rules: the paper's
attribution follows from the explicitly matched controls listed above.
