#!/usr/bin/env python3
"""
make_final_plots_original_grw.py

Final publication-quality figures for the original GRW paper.
Reads from existing CSV files only; does not re-run any simulation.

Input:   output/paper_refinement_original_grw/
Output:  output/paper_refinement_original_grw/publication_figures/
"""

import os
import csv
from datetime import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE = "figure_data"
OUT  = os.path.join("outputs", datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
os.makedirs(OUT, exist_ok=True)

CSV_HEAT    = os.path.join(BASE, "heat_refinement_summary.csv")
CSV_FHN     = os.path.join(BASE, "fhn_refinement_summary.csv")
CSV_BURGERS = os.path.join(BASE, "burgers_domain_sensitivity_summary.csv")

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family":       "serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    12,
    "axes.titlepad":     7,
    "legend.fontsize":   10,
    "legend.framealpha": 0.9,
    "legend.edgecolor":  "#cccccc",
    "xtick.labelsize":   10,
    "ytick.labelsize":   10,
    "axes.linewidth":    0.8,
    "lines.linewidth":   1.8,
    "lines.markersize":  6,
    "axes.grid":         True,
    "grid.alpha":        0.22,
    "grid.linewidth":    0.5,
    "figure.dpi":        150,
    "savefig.dpi":       300,
})

C_BLUE   = "#2166ac"
C_RED    = "#d6604d"
C_GREY   = "#888888"

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    cols = list(rows[0].keys())
    return rows, cols


def col(rows, key, cast=float):
    return [cast(r[key]) for r in rows]


def log_slope(xs, ys):
    """Ordinary least-squares slope of log y vs log x."""
    lx = np.log(np.array(xs, dtype=float))
    ly = np.log(np.array(ys, dtype=float))
    return float(np.polyfit(lx, ly, 1)[0])


# ---------------------------------------------------------------------------
# Figure 1 — Heat fixed-grid diagnostic
# ---------------------------------------------------------------------------

def fig_heat(rows):
    Ns  = col(rows, "N", int)
    L2s = col(rows, "L2_h")

    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    ax.loglog(Ns, L2s, "o-", color=C_BLUE, label=r"$L_h^2$ error")

    ax.set_title("Heat GRW fixed-grid diagnostic")
    ax.set_xlabel("Particle count $N$")
    ax.set_ylabel(r"$L_h^2$ error")
    ax.legend(loc="upper right")
    ax.set_xlim(min(Ns) * 0.6, max(Ns) * 1.6)

    fig.tight_layout()
    out = os.path.join(OUT, "heat_error_fixed_grid_diagnostic_final.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out, log_slope(Ns, L2s)


# ---------------------------------------------------------------------------
# Figure 2 — Scalar FHN refinement
# ---------------------------------------------------------------------------

def fig_fhn(rows):
    Ns    = col(rows, "N", int)
    L2s   = col(rows, "profile_L2")
    ferrs = col(rows, "front_error")

    # O(N^{-1/2}) guide anchored at first profile-L2 point
    N_arr = np.array(Ns, dtype=float)
    guide = L2s[0] * (N_arr[0] / N_arr) ** 0.5

    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    ax.loglog(Ns, L2s,  "o-",  color=C_BLUE, lw=1.8, label="Profile error")
    ax.loglog(Ns, ferrs, "s--", color=C_RED,  lw=1.4, alpha=0.85,
              label="Front error")
    ax.loglog(N_arr, guide, ":", color=C_GREY, lw=1.2,
              label=r"$N^{-1/2}$ guide")

    ax.set_title("Scalar FHN GRW refinement")
    ax.set_xlabel("Particle count $N$")
    ax.set_ylabel("Error")
    ax.legend(loc="lower left")

    fig.tight_layout()
    out = os.path.join(OUT, "fhn_error_vs_N_final.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out, log_slope(Ns, L2s), log_slope(Ns, ferrs)


# ---------------------------------------------------------------------------
# Figure 3 — Burgers domain sensitivity
# ---------------------------------------------------------------------------

def fig_burgers(rows):
    Ls       = col(rows, "L")
    bc_errs  = col(rows, "bc_mismatch_RMSE")
    par_errs = col(rows, "grw_particle_RMSE")

    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    ax.plot(Ls, bc_errs,  "o-", color=C_RED,  label="Finite-domain mismatch")
    ax.plot(Ls, par_errs, "s-", color=C_BLUE, label="GRW reconstruction error")

    ax.set_title("Cole--Hopf Burgers domain sensitivity")
    ax.set_xlabel("Domain size $L$")
    ax.set_ylabel("RMSE")
    ax.set_xticks(Ls)
    ax.set_xticklabels([str(int(v)) for v in Ls])
    ax.legend(loc="upper center")

    fig.tight_layout()
    out = os.path.join(OUT, "burgers_domain_sensitivity_final.pdf")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    heat_rows,    _ = load(CSV_HEAT)
    fhn_rows,     _ = load(CSV_FHN)
    burgers_rows, _ = load(CSV_BURGERS)

    path_heat, slope_heat                     = fig_heat(heat_rows)
    path_fhn,  slope_fhn_L2, slope_fhn_front = fig_fhn(fhn_rows)
    path_burg                                  = fig_burgers(burgers_rows)

    print(f"\nFigures → {OUT}/")
    for p in [path_heat, path_fhn, path_burg]:
        print(f"  {os.path.basename(p)}")

    print("\nSlopes (log-log fit)")
    print(f"  Heat  L2_h vs N:        {slope_heat:+.3f}")
    print(f"  FHN   profile_L2 vs N:  {slope_fhn_L2:+.3f}")
    print(f"  FHN   front_error vs N: {slope_fhn_front:+.3f}")


if __name__ == "__main__":
    main()
