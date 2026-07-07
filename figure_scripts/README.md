# Manuscript Figure Generator

Paper-style plotting pipeline for all eight manuscript figures.

Run from the repository root (or use `python reproduce.py figures`):

```bash
python figure_scripts/regenerate_manuscript_figures.py
```

Default mode loads the archived representative arrays from:

```text
figure_data/representative_figure_arrays.npz
```

and writes PDF/PNG figures to a timestamped folder:

```text
outputs/<YYYY-MM-DD_HH-MM-SS>/
```

To rerun the representative simulations with fixed seed `42` and replace the
archive (verified bit-identical on the tested platform):

```bash
python figure_scripts/regenerate_manuscript_figures.py --rerun
```

The CSV-backed figures (heat fixed-grid diagnostic, FHN refinement, Burgers
domain sensitivity) read the checked-in summaries from:

```text
figure_data/*.csv
```

Re-running `study_paper_refinement.py` writes fresh copies of the same
summaries to `output/paper_refinement_original_grw/`; `reproduce.py verify`
confirms the fresh copies are identical to the checked-in ones.
