# Methods and reproduction for the revised GRW paper

The revised manuscript is *Error attribution in gradient random walk methods
for parabolic equations*. The original preprint title remains in the arXiv
citation. The repository preserves its numerical results and existing equation
identifiers so that old configurations continue to run.

## Mathematical sources and implementation

| Ingredient | Source and manuscript role | Implementation |
|---|---|---|
| Weighted gradient particles and cumulative field recovery | Ghoniem and Sherman (1985), equations (15)-(18), Sections II.2-II.3 | `simulation.py`, `utils.reconstruct_cumulative` |
| Finite-interval reflection and sign rules | Ghoniem and Sherman, Sections III.1-III.2; classical reflected heat kernel | `simulation.py`; `studies/study_t9_heat_two_step.py` evaluates the image-sum reference |
| Deterministic reaction multiplier `1 + dt*f'(u)` | Ghoniem and Sherman, equation (41); our implementation evaluates the reconstructed field after diffusion | `simulation.simulate_fitzhugh_nagumo_grw`, reaction step |
| Exact pointwise empirical-CDF sampling variance | Classical Bernoulli variance; Bertaglia, Pareschi, and Caflisch (2024), Appendix A.2 | `checks/heat_cdf_identity_check.py`, `studies/study_t9_heat_two_step.py` |
| Error at actual comparison points, floor, and particle-count prediction | Application of those established ingredients in the revised paper's heat study | `predictions`, `covariance`, `validation_spec`, and `ensemble` in the two-step study |
| Reaction callback and direct cumulative interface | Reproduction checks against the existing implementation and pinned heat statistics | `checks/interface_regression_checks.py` |

The deterministic reaction multiplier also has the conditional-mean
interpretation discussed in the manuscript. That is a one-step interpretation
given the current field. The finite-time nonlinear particle calculation is
verified by its own experiments.

The `fitzhugh-nagumo` configuration key names the existing scalar solver. The
paper tests a scalar reaction-diffusion equation with a FitzHugh-Nagumo-type
cubic term, with no recovery variable. The exact front, cubic coefficients,
and tested parameter regime are given in the manuscript and configurations.

References:

- A. F. Ghoniem and F. S. Sherman, *Grid-free Simulation of Diffusion Using
  Random Walk Methods*, Journal of Computational Physics 61 (1985), 1-37.
  [DOI](https://doi.org/10.1016/0021-9991(85)90058-0).
- G. Bertaglia, L. Pareschi, and R. E. Caflisch, *Gradient-based Monte Carlo
  methods for relaxation approximations of hyperbolic conservation laws*,
  Journal of Scientific Computing 100, 60 (2024).
  [DOI](https://doi.org/10.1007/s10915-024-02614-1).

## What the predictive calculation assumes

The two-step heat experiment uses independent particles, fixed deterministic
weights, constant diffusivity, and known reflected position distributions.
With fixed half/half allocation, group weights are `A/(N/2)` and `B/(N/2)`.
The variance contributions add with squared weights. The predictions are

```text
E[E_total^2]  = B_h^2 + V_h/N
E[E_bias^2]   = B_h^2 + V_h/(N*S)
E[E_spread^2] = (S-1)/S * V_h/N
N_star       = V_h/B_h^2                         when B_h > 0
N_required  >= V_h/(target^2 - B_h^2)            when target > B_h
```

Here `B_h` includes the comparison-location and reference discrepancies in
the reported discrete norm, and `V_h` is the sampling coefficient. The selected
count is rounded upward to an even number to retain the prescribed allocation.
For zero floor the count bound becomes `V_h/target^2`. At or below a nonzero
floor, no finite count attains the requested expected RMS accuracy.

The rule concerns expected squared error. The held-out ensemble measured
0.008070141 against a prediction of 0.007999948 and a target of 0.008. This is
agreement within 1% of the prediction. A finite ensemble can fall on either
side of a target chosen at the threshold.

The Gaussian and delta-method uncertainty scales used in the analysis are
approximations. They are distinguished from the exact expected squared
moments. The 160 production configurations share trajectories, and the
independent seed-block control assesses correlated sampling fluctuations.

Reaction weights depend on the evolving reconstructed field, and Cole-Hopf
recovery is nonlinear. The heat prediction's assumptions therefore determine
where it applies. The Burgers result is a separate deterministic boundary
control. The particle Burgers runs retain the original fixed endpoints.

## Paper-to-command map

Run all commands from the repository root with the pinned dependencies.

| Material | Command | Inputs and generated outputs |
|---|---|---|
| All eleven paper figures | `python reproduce.py paper` | Reads `figure_data/` and `pinned_ensembles/`; writes `output/final_prepublication_tests/paper_figures/` and its SHA-256 manifest |
| Representative results and five original ensemble/control studies | `python reproduce.py verify --deep` and `python reproduce.py verify-ensembles` | Re-runs solvers and compares the archived arrays/tables |
| Two-step prediction table, design, production, validation | `python reproduce.py verify-two-step` | Recomputes the design, runs seeds 5000-5029 and 7000-7029, and compares the four archived JSON files; writes `output/heat_two_step_reproduction/` |
| Archived two-step design and error-moment consistency | `python reproduce.py verify-two-step --pinned-only` | Reads committed data, checks historical hashes and recomputed predictions; no solver runs |
| Analysis/figure from the fresh two-step rerun | `python studies/analyze_t9_heat_two_step.py --data-dir output/heat_two_step_reproduction --output-dir output/heat_two_step_reproduction` | Prints comparison ratios and the held-out result and writes `Figure_11.pdf` |
| Heat CDF, paired identity, and uncertainty calculation | `python checks/heat_cdf_identity_check.py --fresh` | Uses pinned paired data and additional seed blocks; prints its checks |
| Direct-distribution reference, finite-ensemble moments, six seed blocks | `python checks/two_step_reference_and_seed_blocks.py` | Fixed defaults reproduce the controls, including 4000 ensembles of 12 realizations; prints the results |
| Joint/independent bootstrap table, including heat profiles | `python checks/bootstrap_seed_grouping_check.py --heat-profiles --json output/bootstrap-seed-grouping.json` | Reruns heat profiles and verifies against pinned summaries before fitting; writes intervals for both resampling schemes |
| Reaction callback and cumulative interface | `python checks/interface_regression_checks.py` | Confirms callback/default agreement and pinned paired-heat statistics |

`python reproduce.py verify-all` combines the representative, original ensemble,
and two-step rerun comparisons. The four standalone analytical/bootstrap/interface
checks in the table are separate commands. Figure generation reads data and
does not rerun the solvers.

The bootstrap script's `reported` field refers to the earlier independent
intervals used as regression references. The current manuscript uses its
`joint` field as primary. Small differences in the independently resampled
intervals are identified as such in the manuscript comparison table.

## Original archive and present reproduction

Version 1.0.0 on Zenodo is the original software archive. The revision's
two-step extension and later checks require this updated checkout or the
updated submission supplement. They have not yet been deposited as a new
Zenodo version. Retain the original published DOI and add a new version under
the same concept when depositing the revision.

Historical design records are in `provenance/heat_two_step/`. A rerun writes a
new output directory and preserves the pinned data and original records. The
new outputs reproduce a previously selected design; their timestamps do not
constitute new prospective validation.

## Editing and reproduction conventions

Keep solver formulas in `simulation.py`, initial data and configuration handling
in `config.py`, and reconstruction/plot helpers in `utils.py`. A study script
should state its mathematical reference, parameters, seeds, output destination,
and comparison statistic before its runner. Use ordinary functions, descriptive
comments and docstrings, and one statement per line. Preserve numerical
operation order and random-number consumption during editorial refactors.

Command-line parsing and file/figure generation belong under a `main` entry
point so importing an analysis module does not run it. Use explicit input/output
paths, close files after reading or writing, and keep the committed reference
data separate from fresh runs. The two-step study rejects odd or undersized
particle counts because its fixed allocation requires two equal-sized groups.

The four heat-extension and bootstrap scripts received this formatting pass.
The earlier solver code retains its established layout. `README.md` documents
output replacement, interrupted-run recovery, seed pairing, and custom-case
limits. The usual `verify-all` command covers numerical reruns; the analytical,
bootstrap and interface commands remain separately listed.
