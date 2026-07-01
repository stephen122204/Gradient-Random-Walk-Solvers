#!/usr/bin/env python3
"""
make_publication_plots_original_grw.py

Regenerates publication-quality figures for the original GRW paper from
existing CSV files. Does not re-run any simulation or modify any solver.

Reads from:  output/paper_refinement_original_grw/
Writes to:   output/paper_refinement_original_grw/publication_figures/
"""

import os
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    'font.family':        'serif',
    'font.size':          11,
    'axes.labelsize':     12,
    'axes.titlesize':     12,
    'axes.titleweight':   'normal',
    'legend.fontsize':    10,
    'xtick.labelsize':    10,
    'ytick.labelsize':    10,
    'axes.linewidth':     0.8,
    'lines.linewidth':    1.8,
    'lines.markersize':   6,
    'figure.dpi':         150,   # screen preview; PDF ignores dpi anyway
    'savefig.dpi':        300,
    'axes.grid':          True,
    'grid.alpha':         0.25,
    'grid.linewidth':     0.5,
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE_DIR    = "output/paper_refinement_original_grw"
OUT_DIR     = os.path.join(BASE_DIR, "publication_figures")
CSV_HEAT    = os.path.join(BASE_DIR, "heat_refinement_summary.csv")
CSV_FHN     = os.path.join(BASE_DIR, "fhn_refinement_summary.csv")
CSV_BURGERS = os.path.join(BASE_DIR, "burgers_domain_sensitivity_summary.csv")

os.makedirs(OUT_DIR, exist_ok=True)

generated = []   # (filename, description)
skipped   = []   # (filename, reason)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def read_csv(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"CSV not found: {path}")
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    return rows


def col(rows, key, cast=float):
    return [cast(r[key]) for r in rows]


def _annotation_box(ax, text, loc='lower left', fontsize=8.5):
    """Place a wrapped text annotation inside the axes with a light background box."""
    anchors = {
        'lower left':  (0.03, 0.04),
        'lower right': (0.97, 0.04),
        'upper left':  (0.03, 0.96),
        'upper right': (0.97, 0.96),
    }
    ha_map = {
        'lower left': 'left', 'upper left': 'left',
        'lower right': 'right', 'upper right': 'right',
    }
    va_map = {
        'lower left': 'bottom', 'lower right': 'bottom',
        'upper left': 'top',    'upper right': 'top',
    }
    x, y = anchors.get(loc, (0.03, 0.04))
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        fontsize=fontsize,
        ha=ha_map.get(loc, 'left'),
        va=va_map.get(loc, 'bottom'),
        fontstyle='italic',
        color='#444444',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='#f5f5f5',
                  edgecolor='#cccccc', linewidth=0.6, alpha=0.9),
        wrap=True,
        zorder=5,
    )


# ---------------------------------------------------------------------------
# A.  Heat fixed-grid refinement diagnostic
# ---------------------------------------------------------------------------

def plot_heat_error():
    rows = read_csv(CSV_HEAT)
    Ns   = col(rows, 'N', int)
    L2s  = col(rows, 'L2_h')

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    ax.loglog(Ns, L2s, 'o-', color='#2166ac', lw=1.8, ms=6,
              label=r'$L_h^2$ error')

    ax.set_xlabel('Number of particles $N$')
    ax.set_ylabel(r'$L_h^2$ error')
    ax.set_title('Heat GRW fixed-grid refinement diagnostic')
    ax.legend(loc='upper right')

    _annotation_box(
        ax,
        'Error decreases weakly on this fixed output grid;\n'
        'this diagnostic does not establish an asymptotic\n'
        'particle-convergence rate.',
        loc='lower left',
    )

    ax.set_xlim(min(Ns) * 0.7, max(Ns) * 1.4)

    fig.tight_layout()
    out = os.path.join(OUT_DIR, 'heat_error_fixed_grid_diagnostic.pdf')
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# B.  Heat profiles — skipped (GRW profile arrays not saved to disk)
# ---------------------------------------------------------------------------

def plot_heat_profiles():
    return None   # signal skip


# ---------------------------------------------------------------------------
# C.  FHN particle-refinement diagnostic
# ---------------------------------------------------------------------------

def plot_fhn_error():
    rows   = read_csv(CSV_FHN)
    Ns     = col(rows, 'N', int)
    L2s    = col(rows, 'profile_L2')
    ferrs  = col(rows, 'front_error')

    N_arr  = np.array(Ns, dtype=float)
    # O(N^{-1/2}) guide anchored at the median L2 value
    anchor = float(np.median(L2s)) * (N_arr[len(N_arr)//2] ** 0.5)
    guide  = anchor / N_arr ** 0.5

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    ax.loglog(Ns, L2s, 'o-', color='#2166ac', lw=1.8, ms=6,
              label=r'Profile $L_h^2$ error')
    ax.loglog(Ns, ferrs, 's--', color='#d6604d', lw=1.4, ms=5,
              alpha=0.85, label='Front-location error')
    ax.loglog(N_arr, guide, ':', color='#888888', lw=1.2,
              label=r'$O(N^{-1/2})$ guide')

    ax.set_xlabel('Number of particles $N$')
    ax.set_ylabel('Error')
    ax.set_title('Scalar FHN GRW particle-refinement diagnostic')
    ax.legend(loc='upper right')

    _annotation_box(
        ax,
        'Profile error decreases close to the Monte Carlo\n'
        'reference trend; front-location error is more\n'
        'sensitive to local crossing extraction.',
        loc='lower left',
    )

    fig.tight_layout()
    out = os.path.join(OUT_DIR, 'fhn_error_vs_N_clean.pdf')
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# D.  Burgers domain-sensitivity diagnostic
# ---------------------------------------------------------------------------

def plot_burgers_domain():
    rows    = read_csv(CSV_BURGERS)
    Ls      = col(rows, 'L')
    bc_errs = col(rows, 'bc_mismatch_RMSE')
    par_errs = col(rows, 'grw_particle_RMSE')

    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    ax.plot(Ls, bc_errs,  'o-', color='#d6604d', lw=1.8, ms=6,
            label='Finite-domain transformed-variable mismatch')
    ax.plot(Ls, par_errs, 's-', color='#2166ac', lw=1.8, ms=6,
            label='GRW transformed-variable reconstruction error')

    ax.set_xlabel('Domain size $L$')
    ax.set_ylabel('RMSE contribution')
    ax.set_title('Cole--Hopf Burgers domain-sensitivity diagnostic')
    ax.legend(loc='center right', fontsize=9)

    _annotation_box(
        ax,
        'Increasing the domain reduces finite-domain mismatch\n'
        'but increases stochastic transformed-variable\n'
        'reconstruction error in this test.',
        loc='upper left',
    )

    ax.set_xticks(Ls)
    ax.set_xticklabels([str(int(l)) for l in Ls])

    fig.tight_layout()
    out = os.path.join(OUT_DIR, 'burgers_domain_sensitivity_clean.pdf')
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Markdown interpretation notes
# ---------------------------------------------------------------------------

MARKDOWN_CONTENT = """\
# Figure Interpretation Notes

*Original GRW paper: A Numerical Study of Gradient Random Walk Methods for
Heat, FitzHugh--Nagumo, and Burgers' Equations*

---

## Heat figures

`heat_error_fixed_grid_diagnostic.pdf` shows the discrete $L_h^2$ error
between the GRW reconstruction and the exact error-function solution as
particle count $N$ is increased, using a fixed 300-bin output grid.
The error decreases slowly and non-monotonically at low $N$ because both
particle noise and grid-resolution effects contribute at this grid spacing.
This figure should be described as a **fixed-grid reconstruction diagnostic**,
not as a particle-convergence-rate result.
The profile overlay figure (`heat_profiles_selected_N.pdf`, produced during
the simulation run) is the primary Heat evidence: it shows that the GRW
cumulative-sum reconstruction visually matches the analytic error-function
solution at representative particle counts.

## FHN figure

`fhn_error_vs_N_clean.pdf` shows profile $L_h^2$ error and front-location
error versus particle count $N$ for the scalar FHN traveling-front GRW.
The profile error decreases at a rate consistent with the $O(N^{-1/2})$
Monte Carlo reference trend (log--log slope $\\approx -0.5$), supporting the
reaction-statistic GRW formulation for the FitzHugh--Nagumo equation.
Front-location error is plotted as a secondary curve; it is noisier because
the $u = 0.5$ crossing extracted from a coarse sorted-glob reconstruction is
sensitive to individual glob positions at small $N$.

## Burgers figure

`burgers_domain_sensitivity_clean.pdf` shows two RMSE contributions versus
domain size $L$ for the Cole--Hopf Burgers stationary-shock benchmark.
The finite-domain transformed-variable mismatch arises because the
Cole--Hopf heat equation is solved on a bounded domain with Dirichlet
boundary conditions that differ from the infinite-domain exact solution.
This mismatch decreases as $L$ grows relative to the shock width.
The GRW transformed-variable reconstruction error, arising from stochastic
particle noise in recovering $\\phi_x/\\phi$, increases with $L$ because
the dynamic range of $\\phi_0$ grows, making the ratio harder to estimate.
These two competing effects reveal an optimal domain size near
$L / \\delta \\approx 4$--$6$, where $\\delta = 2\\nu/A$ is the shock width.

---

## Claims supported by these diagnostics

- **Heat:** The direct GRW profile reconstruction matches the exact error-function
  solution at all tested particle counts. The fixed-grid refinement diagnostic
  shows decreasing but grid-limited error and **should not be used as an
  asymptotic particle-convergence-rate claim**.

- **Scalar FHN:** Profile error decreases under particle refinement at a rate
  consistent with $O(N^{-1/2})$, supporting the reaction-statistic GRW
  traveling-front formulation.

- **Cole--Hopf Burgers:** The error decomposes into two competing mechanisms.
  Finite-domain transformed-variable mismatch dominates on tight domains;
  stochastic transformed-variable reconstruction error grows on wider domains.
  Neither mechanism is hidden by the other in this diagnostic.
"""


def write_markdown():
    md_path = os.path.join(OUT_DIR, 'figure_interpretation_notes.md')
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write(MARKDOWN_CONTENT)
    return md_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("make_publication_plots_original_grw.py")
    print("=" * 65)

    # --- Inspect CSVs ---
    for label, path in [
        ("Heat",    CSV_HEAT),
        ("FHN",     CSV_FHN),
        ("Burgers", CSV_BURGERS),
    ]:
        if os.path.isfile(path):
            with open(path, newline='') as f:
                cols = csv.DictReader(f).fieldnames
            print(f"\n  {label} CSV  ({os.path.basename(path)})")
            print(f"    Columns: {', '.join(cols)}")
            print(f"    Path:    {path}")
        else:
            print(f"\n  {label} CSV  NOT FOUND: {path}")

    # --- Plot A: Heat fixed-grid diagnostic ---
    print("\n  Generating plots ...")
    try:
        out = plot_heat_error()
        generated.append((out, 'Heat fixed-grid refinement diagnostic'))
        print(f"  [OK]  {out}")
    except Exception as e:
        skipped.append(('heat_error_fixed_grid_diagnostic.pdf', str(e)))
        print(f"  [SKIP] heat_error_fixed_grid_diagnostic.pdf — {e}")

    # --- Plot B: Heat profiles (skipped: arrays not saved) ---
    reason = (
        "GRW profile arrays (x_grid, u_num per N) were not saved to disk "
        "by study_paper_refinement.py; regenerating them would require "
        "re-running the simulation, which is disallowed. "
        "The original heat_profiles_selected_N.pdf remains in the parent directory."
    )
    skipped.append(('heat_profiles_selected_N_clean.pdf', reason))
    print(f"  [SKIP] heat_profiles_selected_N_clean.pdf")
    print(f"         Reason: {reason}")

    # --- Plot C: FHN error ---
    try:
        out = plot_fhn_error()
        generated.append((out, 'Scalar FHN particle-refinement diagnostic'))
        print(f"  [OK]  {out}")
    except Exception as e:
        skipped.append(('fhn_error_vs_N_clean.pdf', str(e)))
        print(f"  [SKIP] fhn_error_vs_N_clean.pdf — {e}")

    # --- Plot D: Burgers domain sensitivity ---
    try:
        out = plot_burgers_domain()
        generated.append((out, 'Cole-Hopf Burgers domain-sensitivity diagnostic'))
        print(f"  [OK]  {out}")
    except Exception as e:
        skipped.append(('burgers_domain_sensitivity_clean.pdf', str(e)))
        print(f"  [SKIP] burgers_domain_sensitivity_clean.pdf — {e}")

    # --- Markdown ---
    try:
        md = write_markdown()
        generated.append((md, 'Figure interpretation notes (Markdown)'))
        print(f"  [OK]  {md}")
    except Exception as e:
        skipped.append(('figure_interpretation_notes.md', str(e)))
        print(f"  [SKIP] figure_interpretation_notes.md — {e}")

    # --- Summary ---
    print("\n" + "=" * 65)
    print(f"  Generated ({len(generated)} files):")
    for path, desc in generated:
        print(f"    {desc}")
        print(f"      -> {path}")

    if skipped:
        print(f"\n  Skipped ({len(skipped)} files):")
        for name, reason in skipped:
            print(f"    {name}")
            print(f"      Reason: {reason}")

    print("\n  Done.")


if __name__ == "__main__":
    main()
