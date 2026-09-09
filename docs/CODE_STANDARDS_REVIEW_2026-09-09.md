# Code and README review — 2026-09-09

## Scope

Reviewed the current `grw-solvers-v3` checkout against the established separation
of solver, configuration, study, figure, and verification code. The main work
was in the recent two-step heat extension and its associated checks. Existing
solver interfaces and numerical operation order were preserved.

## Changes

- Expanded compressed statements and grouped imports in the four extension and
  bootstrap scripts. Added short docstrings describing mathematical quantities,
  their references, and returned values.
- Put the two-step analysis behind a main entry point. Importing it now leaves
  command-line arguments and files untouched. Figure handles are closed after use.
- Replaced unclosed JSON reads/writes in the extension with scoped file operations.
  Bootstrap CSV readers also close their files explicitly.
- Added normal argument parsing to the bootstrap check. Its JSON output creates
  missing parent directories, so the documented command works in a fresh checkout.
- Reject odd, undersized, and noninteger counts in the two-step study. Previously
  `N // 2` in both groups silently used fewer particles for an odd requested count.
  The published counts are all valid and unchanged.
- Corrected the direct-reference control's printed sampling scale to identify its
  pointwise upper bound, and labeled Gaussian/delta fluctuation scales approximate.
- Replaced two console convergence assertions with fitted-exponent descriptions.
  Clarified the separate tolerances used by representative and ensemble checks.
- Reordered the README installation instructions around the current branch and
  updated supplement. Added output/overwrite behavior, interrupted-run recovery,
  seed pairing, custom-case limits, and the status of the original release citation.

No new solver framework, dependencies, random-number generator, physical model,
or study design was introduced. Core solver files, committed numerical data,
requirements, and historical prediction/provenance records retain their hashes.
The archive citation still identifies the original 1.0.0 release; the revised
Zenodo version remains to be deposited.

## Validation

The review uses Python 3.11.4, NumPy 2.4.2, SciPy 1.17.1, and Matplotlib 3.10.8.

- All 16 unit tests pass, including invalid-allocation rejection and analysis
  import safety.
- All 179 representative checks pass, including 40 bit-identical archived arrays.
- Archived design hashes, predictions, count rules, configuration coverage, and
  finite-ensemble identities pass.
- Full bootstrap output reproduces the previous 16-trend JSON exactly. The newly
  parsed JSON output also works with a previously nonexistent parent directory.
- Direct-reference, moment and six-block reports reproduce the previous outputs
  exactly after accounting for the corrected sampling-scale label.
- The analysis report is unchanged. Figure 11's PDF drawing instructions match
  the pre-refactor figure exactly.
- A clean supplement extraction passes the unit tests and pinned-data checks,
  and regenerates all eleven figures with no pre-existing output directory.
- Compilation of Python modules, targeted formatter checks, local documentation
  links, and Git whitespace checks pass.

The full `verify-all` rerun passes. All 18 original ensemble file comparisons
and all four two-step JSON comparisons match the committed values within their
documented tolerances. The two-step study covers 160 production and 16 validation
configurations with 30 realizations each. At N = 29682, the held-out result remains
0.008070141 against prediction 0.007999948.

Logs are stored in `output/code-standards-review/`; selected final records
accompany the rebuilt submission supplement.
