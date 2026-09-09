# Reproduction verification, 2026-09-08

The Ghoniem-Sherman citation/reproduction pass was checked with Python 3.11.4,
NumPy 2.4.2, SciPy 1.17.1, and Matplotlib 3.10.8.

| Check | Result |
|---|---|
| Unit tests | 14 passed |
| Representative numerical checks | 179 passed, including 40 bit-identical arrays |
| Five original ensemble/control studies | Fresh runs matched all 18 pinned file comparisons |
| Two-step heat study | Fresh predictions, specification, production, and validation matched all four archived JSON files within the documented tolerances |
| Two-step configuration coverage | 160 production and 16 validation configurations, 30 realizations each |
| Historical design | Both recorded prediction/specification hashes matched the pinned files |
| Analytical and supplementary checks | Heat identity/fresh blocks, two-step reference/moment/six-block controls, full joint/independent bootstrap including heat profiles, and interface checks completed successfully |
| Figure generation | 11 of 11 figures regenerated |
| Clean supplement extraction | Unit tests, archived two-step checks, and all eleven figures passed without pre-existing generated output |

The held-out count is unchanged at N = 29682. The rerun measured 0.008070141
against the prediction 0.007999948 and target 0.008. The total-error ratio
across the 160 production configurations ranges from 0.817916 to 1.098724,
with median 0.976829. These reproduce the manuscript's reported results.

All pre-existing pinned data, figure data, expected_values.json, and the solver
files simulation.py, config.py, and utils.py remain byte-identical to the
baseline. Changes affect attribution documentation, reproduction entry points,
the deterministic validation-specification generator, verification, and packaging.
The covariance comment was corrected without changing its calculation.

Full logs are included under check-outputs/ in the updated submission
code-supplement.zip. Run the commands in METHODS_AND_REPRODUCIBILITY.md to
repeat these checks. This is a reproduction of the existing study and does
not constitute a new prospective experiment or a new Zenodo deposition.
