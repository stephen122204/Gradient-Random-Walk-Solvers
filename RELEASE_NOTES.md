# Release-prep notes — Paper 1 branch (`paper1-release-prep`, from `grw-solvers-v1`)

Packaging pass only. **No solver or study code was modified** — the changes
are wrapper/documentation additions (`reproduce.py`, `expected_values.json`,
`MANIFEST.md`, `requirements.txt`, `CITATION.cff`, this file) and rewrites of
the two READMEs. This file collects the evidence and the decisions that
remain with the authors. The repo has NOT been made public and nothing has
been pushed.

## Regression evidence (cardinal rule: bit-identity)

- Baseline established before any change; re-checked after all changes.
- `python reproduce.py verify --deep`: **163 checks, all PASS** —
  - all rows/fields of the three study CSVs identical to the checked-in
    `figure_data/*.csv` (machine-timing column excluded);
  - all 14 representative metrics (Table 2 norms, Table 1 decomposition at
    L = 4, Figure 4 front locations) exact;
  - all 37 archived arrays regenerate **bit-identically** from fresh seed-42
    simulations.
- Verified again inside a fresh venv built from `requirements.txt`
  (see README runtimes; Python 3.11.4, numpy 2.4.2, matplotlib 3.10.8).

## Sensitive-content scan results

- No absolute paths (`/Users/...`), no emails, no tokens/credentials in the
  branch tree. (One stale absolute-path reference in
  `figure_scripts/README.md` was corrected.)
- No `__pycache__`/`.pyc` tracked; `.gitignore` already covers Python
  caches, venvs, outputs, and editor folders.
- `grw-solvers-v1` history is clean: no committed PDFs, zips, or draft
  folders.

## Decisions that remain with the authors (flagged, not decided)

1. **LICENSE** — none added (institution/coauthor decision). Recommend MIT
   or BSD-3 for research code. `TODO(license)` markers are in `README.md`
   and `CITATION.cff`.
2. **How to go public.** This repository also hosts the *unpublished* Paper 2
   (`rb-gbmc-paper2` branch), whose history contains committed output zips,
   draft-manuscript folders, and internal study outputs. Flipping this repo
   to public exposes all branches. Two options:
   - **(Recommended)** create a fresh public repository containing only this
     branch's release tree (history here is clean, but a fresh repo also
     keeps Paper 2 fully private), and point the paper's Code and Data
     Availability statement at it; or
   - keep the single-repo plan and make it public only when Paper 2 is also
     ready for exposure.
   The paper currently cites the `grw-solvers-v1` branch of this repo —
   update the availability statement if option 1 is chosen.
3. **Legacy scripts kept by default** — `make_final_plots_original_grw.py`,
   `make_publication_plots_original_grw.py`, `RUN_WITH_JSON.md` are drafting
   aids, harmless but not needed for reproduction. Confirm keep-or-drop.
4. **Tag** — after review, tag this state (suggest `paper1-code-v1`). Not
   created here so the tag lands on the commit you approve.
5. **CITATION.cff placeholders** — fill in the arXiv id when assigned.

## Known-trap notes carried into the docs (referee-facing)

- Paper configurations are pinned in the entry points; `main.py`/`configs/`
  are exploration tools (stated prominently in README).
- 400-bin representative heat grid vs 300-bin refinement grid: documented in
  README + MANIFEST, matching the paper's Sec. 6.1 disclosure.
- `verify_solver.py`/`verify_grw.py` runs are unseeded (stochastic) — the
  docs now say to use `reproduce.py` for the paper's exact numbers.
