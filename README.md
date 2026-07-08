# Gradient Random Walk Methods: Heat, FitzHugh–Nagumo, and Burgers

This repository contains computer code for reproducing the numerical results
described in the manuscript *A Numerical Study of Gradient Random Walk
Methods for Heat, FitzHugh–Nagumo, and Burgers' Equations* by Stephen Abkin
and Prabir Daripa.

**Paper:** link to be added (arXiv preprint forthcoming).

## Getting Started

```bash
git clone https://github.com/stephen122204/heat_burgers_fhn.git
cd heat_burgers_fhn
git checkout grw-solvers-v1
```

## Reproducing Numerical Results

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Regenerate every study and figure in the paper (fixed seed 42, ~16 s), and
check the results against the paper's reported values:

```bash
python reproduce.py all       # all studies + all 8 figures
python reproduce.py verify    # PASS/FAIL against the paper's values
```

Figures 1, 3, 4, 6, 7 and Table 2 come from the archived representative
arrays in `figure_data/`; Figures 2, 5, 8 and Table 1 come from the
refinement and domain studies (`reproduce.py studies`). The paper's pinned
configurations live in `study_paper_refinement.py` and
`figure_scripts/regenerate_manuscript_figures.py`, which `reproduce.py`
invokes; `main.py` with the JSON files under `configs/` is for exploratory
runs only.

## Citation

See `CITATION.cff`.
