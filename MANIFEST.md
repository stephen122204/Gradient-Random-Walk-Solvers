# Reproduction Manifest

Every table and figure in the paper ("A Numerical Study of Gradient Random
Walk Methods for Heat, FitzHugh--Nagumo, and Burgers' Equations") maps to one
command below. All commands are deterministic: every paper computation uses
NumPy seed 42, re-seeded at the start of each configuration.

Runtimes were measured on the tested machine (Apple Silicon, Python 3.11,
NumPy 2.4.2); everything is seconds-scale.

| Paper artifact | Command | Output | Runtime |
|---|---|---|---|
| Figure 1 (heat profile) | `python reproduce.py figures` | `outputs/<stamp>/heat_comparison.{pdf,png}` | ~5 s (all 8 figures) |
| Figure 2 (heat fixed-grid diagnostic) | `python reproduce.py figures` | `outputs/<stamp>/heat_error_fixed_grid_diagnostic_final.{pdf,png}` | included above |
| Figure 3 (FHN profiles) | `python reproduce.py figures` | `outputs/<stamp>/fhn_comparison.{pdf,png}` | included above |
| Figure 4 (FHN diagnostics) | `python reproduce.py figures` | `outputs/<stamp>/fhn_diagnostics.{pdf,png}` | included above |
| Figure 5 (FHN refinement) | `python reproduce.py figures` | `outputs/<stamp>/fhn_error_vs_N_final.{pdf,png}` | included above |
| Figure 6 (Burgers shock) | `python reproduce.py figures` | `outputs/<stamp>/burgers_comparison.{pdf,png}` | included above |
| Figure 7 (Burgers pipeline diagnostics) | `python reproduce.py figures` | `outputs/<stamp>/burgers_diagnostics.{pdf,png}` | included above |
| Figure 8 (Burgers domain sensitivity) | `python reproduce.py figures` | `outputs/<stamp>/burgers_domain_sensitivity_final.{pdf,png}` | included above |
| Table 1 (Burgers error decomposition, L = 4,6,8,10) | `python reproduce.py studies` | `output/paper_refinement_original_grw/burgers_domain_sensitivity_summary.csv` | < 1 s |
| Table 2 (representative error norms) | `python reproduce.py verify` (printed), or `python figure_scripts/regenerate_manuscript_figures.py` (printed at end) | terminal table | ~5 s |
| Fig. 2 underlying data (heat N-sweep) | `python reproduce.py studies` | `output/paper_refinement_original_grw/heat_refinement_summary.csv` | ~10 s |
| Fig. 5 underlying data (FHN N-sweep + front errors) | `python reproduce.py studies` | `output/paper_refinement_original_grw/fhn_refinement_summary.csv` | < 1 s |
| Front-location numbers quoted in the Fig. 4 discussion (12.0135 / 11.8100 / 0.203) | `python reproduce.py verify` | printed as `fhn_front.*` checks | included in verify |
| Everything | `python reproduce.py all` | all of the above | ~16 s |

## Verification

```bash
python reproduce.py verify          # ~12 s
python reproduce.py verify --deep   # ~13 s, additionally re-runs the three
                                    # representative simulations and requires
                                    # bit-identical arrays vs figure_data/*.npz
```

`verify` re-runs the three studies and compares **every reported field**
(all rows of all three CSVs, machine-specific `runtime_s` excluded) plus the
Table 2 norms, the Table 1 decomposition at L = 4, and the Figure 4 front
locations against `expected_values.json`, printing PASS/FAIL per item.
Current status on the tested platform: **163 checks, all PASS, including
bit-identical regeneration of all 37 archived arrays.**

## Data provenance notes

- The figure pipeline plots Figures 2, 5, 8 from the **checked-in**
  `figure_data/*.csv` (the paper's data of record). `reproduce.py studies`
  writes fresh copies to `output/paper_refinement_original_grw/`; `verify`
  confirms the two are identical.
- The representative arrays behind Figures 1, 3, 4, 6, 7 and Table 2 are
  archived in `figure_data/representative_figure_arrays.npz` (seed 42).
  `python reproduce.py figures --rerun-arrays` regenerates them from fresh
  simulations; `verify --deep` proves the regeneration is bit-identical.
- **Output-grid note (matches paper Sec. 6.1):** the representative heat run
  (Figure 1, Table 2 heat row) uses a **400-bin** output grid; the heat
  refinement study (Figure 2) uses a **300-bin** grid. Both values sit on
  their grid-controlled floors (0.0119 vs 0.0158 at N = 50,000), which is the
  h-scaling discussed in the paper.
- Artifacts of the second paper (relaxation-based GBMC) are **not** in this
  branch; they live on the `rb-gbmc-paper2` branch and are not regenerable
  from this tree.
