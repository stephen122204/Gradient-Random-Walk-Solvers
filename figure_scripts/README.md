# Manuscript Figure Generator

This folder contains the paper-style plotting pipeline copied from the local
Desktop draft paper and adapted to run from this repository.

Run from the repository root:

```bash
MPLCONFIGDIR=/tmp/grw-figure-mpl \
XDG_CACHE_HOME=/tmp/grw-figure-cache \
python3 figure_scripts/regenerate_manuscript_figures.py
```

Default mode loads archived representative arrays from:

```text
figure_data/representative_figure_arrays.npz
```

and writes PDF/PNG figures to:

```text
_DRAFT__GRW_Methods (1)/figures/
```

To rerun the representative simulations with fixed seed `42` and update the
archive:

```bash
MPLCONFIGDIR=/tmp/grw-figure-mpl \
XDG_CACHE_HOME=/tmp/grw-figure-cache \
python3 figure_scripts/regenerate_manuscript_figures.py --rerun
```

The refinement/domain figures read CSV summaries from:

```text
output/paper_refinement_original_grw/
```
