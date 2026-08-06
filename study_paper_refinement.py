#!/usr/bin/env python3
"""
study_paper_refinement.py

Publication-oriented numerical evidence for the original GRW paper.
Three sequential studies run to completion, writing all outputs to
output/paper_refinement_original_grw/.

Study 1 – Heat particle refinement (N convergence)
Study 2 – Scalar FHN traveling-wave refinement (N convergence)
Study 3 – Burgers Cole-Hopf domain sensitivity (L sweep)
"""

import sys
import os
import time
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from simulation import (
    simulate_heat_equation,
    simulate_fitzhugh_nagumo_grw,
    simulate_burgers_cole_hopf_grw,
    _reference_phi_heat_fd,
)
from verify_solver import (
    exact_heat_step,
    exact_fhn_traveling_wave,
)
from config import SimulationConfig, generate_fhn_steady_ic

OUT_DIR = "output/paper_refinement_original_grw"
os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# STUDY 1: Heat particle refinement
# ---------------------------------------------------------------------------

def study1_heat_refinement():
    print("\n" + "=" * 65)
    print("STUDY 1: Heat particle refinement")
    print("=" * 65)

    alpha = 0.1
    T = 0.5
    dt = 0.001
    L = 10.0
    x0 = 5.0
    uL = 0.0
    uR = 1.0
    N_values = [1000, 2000, 5000, 10000, 20000, 50000]
    nbins = 300

    bc = {
        'LEFT':  {'type': 'dirichlet', 'value': uL},
        'RIGHT': {'type': 'dirichlet', 'value': uR},
    }
    jump_height = uR - uL

    rows = []
    selected_profiles = {}  # N -> (x_grid, u_num, u_exact)

    for N in N_values:
        np.random.seed(42)

        globs = [{'position': x0, 'value': jump_height / N} for _ in range(N)]

        cfg = SimulationConfig(
            equation_type='heat',
            domain_type='finite',
            domain_size=L,
            boundary_conditions=bc,
            diff_constant=alpha,
            time_step=dt,
            total_time=T,
            num_points=N,
            initial_conditions=[(x0, jump_height / N)] * N,
            reaction_term=False,
        )

        t0 = time.perf_counter()
        final_globs = simulate_heat_equation(globs, cfg)
        runtime = time.perf_counter() - t0

        positions = np.array([g['position'] for g in final_globs])
        values   = np.array([g['value']    for g in final_globs])

        total_weight = float(np.sum(values))
        glob_std     = float(np.std(positions))

        # Reconstruct u on 300-bin grid via histogram + cumsum
        edges  = np.linspace(0.0, L, nbins + 1)
        bin_w, _ = np.histogram(positions, bins=edges, weights=values)
        u_num  = uL + np.cumsum(bin_w)
        x_grid = 0.5 * (edges[:-1] + edges[1:])
        dx     = float(edges[1] - edges[0])

        u_exact = exact_heat_step(x_grid, T, x0, uL, uR, alpha)

        diff    = u_num - u_exact
        L2_h    = float(np.sqrt(np.sum(diff ** 2) * dx))
        Linf    = float(np.max(np.abs(diff)))
        ref_norm = float(np.sqrt(np.sum(u_exact ** 2) * dx))
        rel_L2  = L2_h / ref_norm if ref_norm > 1e-12 else float('nan')

        rows.append({
            'N': N,
            'L2_h': L2_h,
            'Linf': Linf,
            'rel_L2': rel_L2,
            'total_weight': total_weight,
            'glob_std': glob_std,
            'runtime_s': runtime,
        })

        if N in (1000, 5000, 50000):
            selected_profiles[N] = (x_grid.copy(), u_num.copy(), u_exact.copy())

        print(f"  N={N:6d}: L2={L2_h:.4e}  Linf={Linf:.4e}  "
              f"rel_L2={rel_L2:.4e}  rt={runtime:.2f}s")

    # --- CSV ---
    csv_path = os.path.join(OUT_DIR, 'heat_refinement_summary.csv')
    fields = ['N', 'L2_h', 'Linf', 'rel_L2', 'total_weight', 'glob_std', 'runtime_s']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\n  Saved: {csv_path}")

    # --- Plot 1: log-log L2 vs N ---
    Ns    = [r['N']    for r in rows]
    L2s   = [r['L2_h'] for r in rows]
    N_arr = np.array(Ns, dtype=float)
    ref   = L2s[0] * (N_arr[0] / N_arr) ** 0.5

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(Ns, L2s, 'o-', color='steelblue', lw=2, ms=7, label='GRW $L^2$ error')
    ax.loglog(N_arr, ref, '--', color='gray', lw=1.5, label=r'$O(N^{-1/2})$')
    ax.set_xlabel('Number of particles N')
    ax.set_ylabel('$L^2$ error')
    ax.set_title('Heat GRW: error vs particle count')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    pdf1 = os.path.join(OUT_DIR, 'heat_error_vs_N.pdf')
    fig.savefig(pdf1, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {pdf1}")

    # --- Plot 2: profiles at N=1000, 5000, 50000 vs exact ---
    fig, ax = plt.subplots(figsize=(7, 4))
    palette = {1000: '#e41a1c', 5000: '#4daf4a', 50000: '#377eb8'}
    for N_sel in (1000, 5000, 50000):
        xg, un, _ = selected_profiles[N_sel]
        ax.plot(xg, un, '-', color=palette[N_sel], lw=1.5, label=f'GRW N={N_sel}')
    xg_ref, _, ue_ref = selected_profiles[50000]
    ax.plot(xg_ref, ue_ref, 'k--', lw=2.0, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x, T)')
    ax.set_title(f'Heat GRW profiles: N = 1000, 5000, 50000  (T={T})')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    pdf2 = os.path.join(OUT_DIR, 'heat_profiles_selected_N.pdf')
    fig.savefig(pdf2, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {pdf2}")

    # --- Stdout summary table ---
    print("\n  Summary table:")
    hdr = f"  {'N':>8}  {'L2_h':>12}  {'Linf':>12}  {'rel_L2':>12}  {'runtime_s':>10}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['N']:>8d}  {r['L2_h']:>12.4e}  {r['Linf']:>12.4e}  "
              f"{r['rel_L2']:>12.4e}  {r['runtime_s']:>10.2f}")

    return rows


# ---------------------------------------------------------------------------
# STUDY 2: Scalar FHN traveling-wave refinement
# ---------------------------------------------------------------------------

def study2_fhn_refinement():
    print("\n" + "=" * 65)
    print("STUDY 2: Scalar FHN refinement")
    print("=" * 65)

    a        = 0.25
    D        = 0.5
    T        = 9.0
    dt       = 0.01
    L        = 30.0
    x_center = 15.0
    N_values = [100, 200, 500, 1000, 2000]
    n_grid   = 500

    theta         = np.sqrt(2.0) * (0.5 - a)
    x_front_exact = x_center - theta * T

    bc = {
        'LEFT':  {'type': 'dirichlet', 'value': 0},
        'RIGHT': {'type': 'dirichlet', 'value': 0},
    }

    x_grid = np.linspace(0.0, L, n_grid)
    dx     = float(x_grid[1] - x_grid[0])
    u_exact = exact_fhn_traveling_wave(x_grid, T, a, x_center=x_center)

    rows = []

    for N in N_values:
        np.random.seed(42)

        ic    = generate_fhn_steady_ic(N, a, x_center=x_center)
        globs = [{'position': float(pos), 'value': float(w)} for pos, w in ic]

        cfg = SimulationConfig(
            equation_type='fitzhugh-nagumo',
            domain_type='finite',
            domain_size=L,
            boundary_conditions=bc,
            diff_constant=D,
            time_step=dt,
            total_time=T,
            num_points=N,
            initial_conditions=ic,
            reaction_term=True,
            a=a, b=0.5, tau=10.0,
            fhn_ic_type='steady_solution',
        )

        t0 = time.perf_counter()
        final_globs = simulate_fitzhugh_nagumo_grw(globs, cfg)
        runtime = time.perf_counter() - t0

        # Reconstruct u via sorted cumsum on 500-point grid
        x_pos = np.array([g['position'] for g in final_globs])
        w_val = np.array([float(g['value']) for g in final_globs])
        order = np.argsort(x_pos)
        x_s   = x_pos[order]
        w_s   = w_val[order]

        total_weight = float(np.sum(w_s))

        cumw      = np.cumsum(w_s)
        full_cumw = np.concatenate([[0.0], cumw])
        idx_grid  = np.searchsorted(x_s, x_grid, side='right')
        u_num     = full_cumw[idx_grid]

        # Profile error metrics
        diff          = u_num - u_exact
        profile_L2    = float(np.sqrt(np.sum(diff ** 2) * dx))
        ref_norm      = float(np.sqrt(np.sum(u_exact ** 2) * dx))
        profile_rel_L2 = profile_L2 / ref_norm if ref_norm > 1e-12 else float('nan')
        profile_Linf  = float(np.max(np.abs(diff)))

        # Front location: u=0.5 crossing via linear interpolation
        diff05    = u_num - 0.5
        crossings = np.where(np.diff(np.sign(diff05)))[0]
        if len(crossings) > 0:
            i = crossings[0]
            denom = u_num[i + 1] - u_num[i]
            if abs(denom) > 1e-15:
                x_front_num = float(x_grid[i] + (0.5 - u_num[i]) / denom
                                    * (x_grid[i + 1] - x_grid[i]))
            else:
                x_front_num = float(x_grid[i])
        else:
            x_front_num = float(x_grid[np.argmin(np.abs(diff05))])

        front_error = abs(x_front_num - x_front_exact)

        rows.append({
            'N': N,
            'profile_L2': profile_L2,
            'profile_rel_L2': profile_rel_L2,
            'profile_Linf': profile_Linf,
            'front_location_num': x_front_num,
            'front_location_exact': x_front_exact,
            'front_error': front_error,
            'total_weight': total_weight,
            'runtime_s': runtime,
        })

        print(f"  N={N:5d}: L2={profile_L2:.4e}  Linf={profile_Linf:.4e}  "
              f"front_err={front_error:.4e}  rt={runtime:.2f}s")

    # --- CSV ---
    csv_path = os.path.join(OUT_DIR, 'fhn_refinement_summary.csv')
    fields = [
        'N', 'profile_L2', 'profile_rel_L2', 'profile_Linf',
        'front_location_num', 'front_location_exact', 'front_error',
        'total_weight', 'runtime_s',
    ]
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\n  Saved: {csv_path}")

    # --- Plot: two curves on log-log ---
    Ns        = [r['N']           for r in rows]
    L2s       = [r['profile_L2']  for r in rows]
    ferrs     = [r['front_error'] for r in rows]
    N_arr     = np.array(Ns, dtype=float)
    ref_line  = L2s[0] * (N_arr[0] / N_arr) ** 0.5

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.loglog(Ns, L2s,  'o-', color='steelblue',  lw=2, ms=7, label='Profile $L^2$ error')
    ax.loglog(Ns, ferrs, 's-', color='darkorange', lw=2, ms=7, label='Front-location error')
    ax.loglog(N_arr, ref_line, '--', color='gray', lw=1.5, label=r'$O(N^{-1/2})$')
    ax.set_xlabel('Number of particles N')
    ax.set_ylabel('Error')
    ax.set_title('Scalar FHN GRW: error vs particle count')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    pdf = os.path.join(OUT_DIR, 'fhn_error_vs_N.pdf')
    fig.savefig(pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {pdf}")

    # --- Stdout summary table ---
    print(f"\n  Summary table  (theta={theta:.6f}, exact front at T={T}: x={x_front_exact:.4f})")
    hdr = (f"  {'N':>6}  {'profile_L2':>12}  {'profile_Linf':>13}  "
           f"{'front_err':>12}  {'runtime_s':>10}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['N']:>6d}  {r['profile_L2']:>12.4e}  {r['profile_Linf']:>13.4e}  "
              f"{r['front_error']:>12.4e}  {r['runtime_s']:>10.2f}")

    return rows


# ---------------------------------------------------------------------------
# STUDY 3: Burgers Cole-Hopf domain sensitivity
# ---------------------------------------------------------------------------

def _compute_phi0_stationary_shock(L, N, nu, A):
    """
    Compute phi0 on a uniform N-point grid for the stationary-shock IC,
    replicating the normalization used inside simulate_burgers_cole_hopf_grw.
    Returns (x_grid, phi0, phi0_0, phi0_L).
    """
    xc  = L / 2.0
    x   = np.linspace(0.0, L, N)
    u0  = -A * np.tanh(A * (x - xc) / (2.0 * nu))

    # Psi0(x) = integral_0^x u0(s) ds  (trapezoidal rule)
    Psi0 = np.zeros(N)
    for i in range(1, N):
        Psi0[i] = Psi0[i - 1] + 0.5 * (u0[i - 1] + u0[i]) * (x[i] - x[i - 1])

    log_phi0 = -Psi0 / (2.0 * nu)
    log_phi0 -= log_phi0.max()          # normalize: max(log_phi0) = 0  =>  phi0_max = 1
    phi0 = np.exp(np.clip(log_phi0, -700.0, 0.0))

    return x, phi0, float(phi0[0]), float(phi0[-1])


def study3_burgers_domain_sensitivity():
    print("\n" + "=" * 65)
    print("STUDY 3: Burgers Cole-Hopf domain sensitivity")
    print("=" * 65)

    A        = 1.0
    nu       = 0.5
    dt       = 0.005
    T        = 0.5
    L_values = [4.0, 6.0, 8.0, 10.0]

    bc = {
        'LEFT':  {'type': 'dirichlet', 'value': 0},
        'RIGHT': {'type': 'dirichlet', 'value': 0},
    }

    rows = []

    for L in L_values:
        N           = int(100 * L)
        xc          = L / 2.0
        shock_width = 2.0 * nu / A

        print(f"\n  L={L:.1f}  N={N}  shock_width={shock_width:.3f}  "
              f"L/shock_width={L / shock_width:.1f}")

        x_out  = np.linspace(0.0, L, N)
        dx_out = float(x_out[1] - x_out[0])

        # Exact stationary shock (time-independent)
        u_exact = -A * np.tanh(A * (x_out - xc) / (2.0 * nu))

        # phi0 for FD reference (deterministic, computed before RNG seed)
        _, phi0, phi0_0, phi0_L = _compute_phi0_stationary_shock(L, N, nu, A)
        phi_ref_fd = _reference_phi_heat_fd(phi0.copy(), x_out, nu, T, phi0_0, phi0_L)
        phi_x_fd   = np.gradient(phi_ref_fd, dx_out)
        phi_safe_fd = np.where(np.abs(phi_ref_fd) < 1e-12, 1e-12, phi_ref_fd)
        u_fd = -2.0 * nu * phi_x_fd / phi_safe_fd

        # GRW run
        u0_ic = -A * np.tanh(A * (x_out - xc) / (2.0 * nu))
        globs = [{'position': float(x_out[i]), 'value': [float(u0_ic[i])]}
                 for i in range(N)]

        cfg = SimulationConfig(
            equation_type='burgers',
            domain_type='finite',
            domain_size=L,
            boundary_conditions=bc,
            diff_constant=nu,
            time_step=dt,
            total_time=T,
            num_points=N,
            initial_conditions=list(zip(x_out.tolist(), u0_ic.tolist())),
            reaction_term=False,
            burgers_mode='cole_hopf_grw',
            burgers_ic_type='stationary_shock',
            burgers_ic_amplitude=A,
        )

        np.random.seed(42)
        t0 = time.perf_counter()
        result_globs = simulate_burgers_cole_hopf_grw(globs, cfg, _diag_dir=None)
        runtime = time.perf_counter() - t0

        # GRW output is already on x_out = np.linspace(0, L, N)
        u_grw = np.array([float(g['value'][0]) for g in result_globs])

        # Error decomposition on the N-point output grid
        diff_total    = u_grw - u_exact
        diff_bc       = u_fd  - u_exact
        diff_particle = u_grw - u_fd

        total_RMSE       = float(np.sqrt(np.mean(diff_total    ** 2)))
        bc_mismatch_RMSE = float(np.sqrt(np.mean(diff_bc       ** 2)))
        grw_particle_RMSE = float(np.sqrt(np.mean(diff_particle ** 2)))

        L2_h     = float(np.sqrt(np.sum(diff_total ** 2) * dx_out))
        ref_norm = float(np.sqrt(np.sum(u_exact   ** 2) * dx_out))
        rel_L2   = L2_h / ref_norm if ref_norm > 1e-12 else float('nan')

        rows.append({
            'L': L,
            'N': N,
            'shock_width': shock_width,
            'L_over_shock_width': L / shock_width,
            'total_RMSE': total_RMSE,
            'bc_mismatch_RMSE': bc_mismatch_RMSE,
            'grw_particle_RMSE': grw_particle_RMSE,
            'L2_h': L2_h,
            'rel_L2': rel_L2,
            'runtime_s': runtime,
        })

        print(f"  -> total_RMSE={total_RMSE:.4e}  bc_mismatch={bc_mismatch_RMSE:.4e}  "
              f"particle={grw_particle_RMSE:.4e}  rt={runtime:.2f}s")

    # --- CSV ---
    csv_path = os.path.join(OUT_DIR, 'burgers_domain_sensitivity_summary.csv')
    fields = [
        'L', 'N', 'shock_width', 'L_over_shock_width',
        'total_RMSE', 'bc_mismatch_RMSE', 'grw_particle_RMSE',
        'L2_h', 'rel_L2', 'runtime_s',
    ]
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\n  Saved: {csv_path}")

    # --- Plot: bc_mismatch and grw_particle vs L ---
    L_vals       = [r['L']                for r in rows]
    bc_errs      = [r['bc_mismatch_RMSE'] for r in rows]
    particle_errs = [r['grw_particle_RMSE'] for r in rows]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(L_vals, bc_errs,       'o-', color='darkorange', lw=2, ms=7,
            label='deterministic component RMSE')
    ax.plot(L_vals, particle_errs, 's-', color='steelblue',  lw=2, ms=7,
            label='GRW reconstruction component RMSE')
    ax.set_xlabel('Domain size L')
    ax.set_ylabel('RMSE contribution')
    ax.set_title('Cole-Hopf Burgers: error decomposition vs domain size')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    pdf = os.path.join(OUT_DIR, 'burgers_domain_sensitivity.pdf')
    fig.savefig(pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {pdf}")

    # --- Stdout summary table ---
    print("\n  Summary table:")
    hdr = (f"  {'L':>5}  {'N':>5}  {'total_RMSE':>12}  "
           f"{'bc_mismatch':>12}  {'particle':>12}  {'runtime_s':>10}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in rows:
        print(f"  {r['L']:>5.1f}  {r['N']:>5d}  {r['total_RMSE']:>12.4e}  "
              f"{r['bc_mismatch_RMSE']:>12.4e}  {r['grw_particle_RMSE']:>12.4e}  "
              f"{r['runtime_s']:>10.2f}")

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    t_start = time.perf_counter()

    rows1 = study1_heat_refinement()
    rows2 = study2_fhn_refinement()
    rows3 = study3_burgers_domain_sensitivity()

    wall = time.perf_counter() - t_start

    print("\n" + "=" * 65)
    print("ALL STUDIES COMPLETE")
    print("=" * 65)

    files = [
        'heat_refinement_summary.csv',
        'heat_error_vs_N.pdf',
        'heat_profiles_selected_N.pdf',
        'fhn_refinement_summary.csv',
        'fhn_error_vs_N.pdf',
        'burgers_domain_sensitivity_summary.csv',
        'burgers_domain_sensitivity.pdf',
    ]
    print("\n  Output files:")
    for f in files:
        p = os.path.join(OUT_DIR, f)
        exists = os.path.isfile(p)
        print(f"    {'OK' if exists else 'MISSING':6s}  {p}")

    # Study 1 convergence check: slope of log-log fit
    Ns1  = np.log([r['N']    for r in rows1])
    L2s1 = np.log([r['L2_h'] for r in rows1])
    slope1 = float(np.polyfit(Ns1, L2s1, 1)[0])

    # Study 2 profile slope
    Ns2  = np.log([r['N']          for r in rows2])
    L2s2 = np.log([r['profile_L2'] for r in rows2])
    slope2 = float(np.polyfit(Ns2, L2s2, 1)[0])

    print(f"\n  Convergence summary:")
    print(f"    Study 1 (heat)  – log-log slope: {slope1:.3f}  "
          f"(expected ≈ -0.5 for O(N^{{-1/2}}))")
    print(f"    Study 2 (FHN)   – log-log slope: {slope2:.3f}  "
          f"(expected ≈ -0.5 for O(N^{{-1/2}}))")

    print(f"\n  Study 3 error decomposition:")
    for r in rows3:
        frac_bc  = r['bc_mismatch_RMSE']  / max(r['total_RMSE'], 1e-30) * 100
        frac_par = r['grw_particle_RMSE'] / max(r['total_RMSE'], 1e-30) * 100
        print(f"    L={r['L']:.1f}: bc_mismatch={frac_bc:.0f}%  particle={frac_par:.0f}%  "
              f"of total RMSE={r['total_RMSE']:.3e}")

    print(f"\n  Total wall time: {wall:.1f}s")

    # Paper-claim assessment
    print("\n  Paper-claim assessment:")
    if -0.7 <= slope1 <= -0.3:
        print("    [PASS] Heat GRW converges at O(N^{-1/2}) rate.")
    else:
        print(f"   [NOTE] Heat slope {slope1:.3f} deviates from -0.5; "
              f"check particle count range or dt.")
    if -0.7 <= slope2 <= -0.3:
        print("    [PASS] Scalar FHN GRW converges at O(N^{-1/2}) rate.")
    else:
        print(f"   [NOTE] FHN slope {slope2:.3f} deviates from -0.5.")
    bc_dom = all(r['bc_mismatch_RMSE'] > r['grw_particle_RMSE'] for r in rows3)
    if bc_dom:
        print("    [PASS] deterministic component dominates the GRW reconstruction component across all L values.")
    else:
        print("    [NOTE] GRW reconstruction component rivals the deterministic component for some L values.")


if __name__ == "__main__":
    main()
