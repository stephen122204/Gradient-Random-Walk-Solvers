"""Task 7: Paired heat output-grid study.

Same physical and numerical parameters, particle counts, and 30 seeds as the
heat ensemble study (t4): alpha=0.5, domain [0,4], x0=2, increasing unit step,
T=0.5, dt=0.005, N in {500,...,50000}, seeds 42..71.

For each (N, seed) ONE simulation is run; the identical final particle set is
then reconstructed under three output-grid treatments:

  coupled  : sorted cumulative sum interpolated to a uniform N-point grid
             (exactly the t4 ensemble treatment; the coupled rows and fitted
             slopes/CIs reproduce t4's published table digit-for-digit)
  fixed300 : glob weights binned on a fixed 300-bin grid, cumulative sum,
             compared with the exact profile at bin centers
             (exactly the Paper-1 fixed-grid diagnostic treatment)
  fixed400 : the same with 400 bins

Because all three treatments share the particle trajectories seed-for-seed,
differences between the curves isolate the output-grid treatment alone.
The bias--spread--total decomposition is computed per treatment: the fixed-grid
floor is expected to appear in E_bias (deterministic binning offset, common to
all seeds) while E_spread continues to decrease near N^{-1/2}.

Output: output/final_prepublication_tests/heat_grid_paired/
"""
import csv
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import SimulationConfig
from verify_solver import exact_heat_step

OUT_BASE = 'output/final_prepublication_tests/heat_grid_paired'

FIXED_BINS = (300, 400)


def _mk(*parts):
    p = os.path.join(*parts)
    os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    return p


def _savefig(fig, path_noext):
    for ext in ('png', 'pdf'):
        p = path_noext + '.' + ext
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
        fig.savefig(p, dpi=150, bbox_inches='tight')
    plt.close(fig)


def _run_heat_one(N, alpha, T, dt, L, x0, uL, uR, seed):
    """Single Heat GRW run; identical to t4's runner (same RNG consumption).

    Returns the raw particle data (positions, weights) in addition to the
    sorted staircase, so every reconstruction arm works from one simulation.
    """
    from simulation import simulate_heat_equation
    np.random.seed(seed)

    weight = (uL - uR) / N
    ic = [(x0, weight)] * N

    cfg = SimulationConfig(
        equation_type='heat',
        domain_type='Finite',
        domain_size=L,
        boundary_conditions={'LEFT': {'type': 'Dirichlet', 'value': float(uR)},
                             'RIGHT': {'type': 'Dirichlet', 'value': float(uL)}},
        diff_constant=alpha,
        time_step=dt,
        total_time=T,
        num_points=N,
        initial_conditions=ic,
        reaction_term=False,
    )
    globs = [{'position': float(p), 'value': float(w)} for p, w in ic]
    t0 = time.perf_counter()
    result = simulate_heat_equation(globs, cfg)
    elapsed = time.perf_counter() - t0

    x_pos = np.array([g['position'] for g in result])
    w_arr = np.array([g['value']    for g in result])
    order    = np.argsort(x_pos)
    x_sorted = x_pos[order]
    u_sorted = float(uR) + np.cumsum(w_arr[order])
    return x_pos, w_arr, x_sorted, u_sorted, elapsed


def _bootstrap_slope_ci(x_arr, y_arr, n_boot=2000, ci=0.95, rng=None):
    if rng is None:
        rng = np.random.default_rng(77)
    valid = np.isfinite(x_arr) & np.isfinite(y_arr) & (x_arr > 0) & (y_arr > 0)
    n = int(valid.sum())
    if n < 3:
        return float('nan'), float('nan')
    lx = np.log10(x_arr[valid])
    ly = np.log10(y_arr[valid])
    slopes = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            c = np.polyfit(lx[idx], ly[idx], 1)
            slopes.append(c[0])
        except Exception:
            pass
    if not slopes:
        return float('nan'), float('nan')
    lo = float(np.percentile(slopes, 100*(1-ci)/2))
    hi = float(np.percentile(slopes, 100*(1+ci)/2))
    return lo, hi


def _decompose(u_arr, u_exact, dx):
    """Bias--spread--total decomposition over stacked realization profiles."""
    u_mean = u_arr.mean(axis=0)
    E_bias   = float(np.sqrt(np.sum((u_mean - u_exact)**2 * dx)))
    E_spread = float(np.sqrt(np.mean(np.sum((u_arr - u_mean[None, :])**2 * dx, axis=1))))
    E_total  = float(np.sqrt(np.mean(np.sum((u_arr - u_exact[None, :])**2 * dx, axis=1))))
    E_total_check = float(np.sqrt(E_bias**2 + E_spread**2))
    identity_err  = abs(E_total - E_total_check) / max(E_total, 1e-15)
    return E_bias, E_spread, E_total, identity_err


def run_task7(N_seq=None, S=30, base_seed=42,
              alpha=0.5, T=0.5, L=4.0, x0=2.0, uL=1.0, uR=0.0,
              dt=0.005):
    os.makedirs(OUT_BASE, exist_ok=True)
    if N_seq is None:
        N_seq = [500, 1000, 2000, 5000, 10000, 20000, 50000]

    arms = ['coupled'] + [f'fixed{M}' for M in FIXED_BINS]

    print(f"{'='*60}\n  Task 7: Paired heat output-grid study\n"
          f"  N_seq={N_seq}  S={S}  arms={arms}\n{'='*60}")

    # Fixed-grid machinery (Paper-1 diagnostic treatment): bin edges, centers,
    # exact profile at centers. u rises from uR(=0) at the left boundary.
    fixed = {}
    for M in FIXED_BINS:
        edges   = np.linspace(0.0, L, M + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        dxM     = float(edges[1] - edges[0])
        u_ex_M  = exact_heat_step(centers, T, x0, uR, uL, alpha)
        fixed[M] = {'edges': edges, 'centers': centers, 'dx': dxM, 'u_exact': u_ex_M}

    results = {arm: [] for arm in arms}

    for N in N_seq:
        x_grid = np.linspace(0.0, L, N)
        u_exact = exact_heat_step(x_grid, T, x0, uR, uL, alpha)
        dx = float(x_grid[1] - x_grid[0])

        print(f"\n  N={N}")
        seeds = [base_seed + i for i in range(S)]
        runs = {arm: [] for arm in arms}
        runtimes = []

        for s_idx, seed in enumerate(seeds):
            x_pos, w_arr, x_out, u_out, elapsed = _run_heat_one(
                N, alpha, T, dt, L, x0, uL, uR, seed)
            runtimes.append(elapsed)

            # coupled treatment: identical operations to t4
            if len(x_out) != N or not np.allclose(x_out, x_grid, atol=1e-10):
                u_coupled = np.interp(x_grid, x_out, u_out,
                                      left=float(uR), right=float(uL))
            else:
                u_coupled = u_out
            runs['coupled'].append(u_coupled)

            # fixed-grid treatments: identical operations to the Paper-1
            # refinement diagnostic (histogram of weights, cumulative sum,
            # offset by the left boundary value uR=0)
            for M in FIXED_BINS:
                bin_w, _ = np.histogram(x_pos, bins=fixed[M]['edges'], weights=w_arr)
                u_fixed = float(uR) + np.cumsum(bin_w)
                runs[f'fixed{M}'].append(u_fixed)

            if (s_idx + 1) % 10 == 0:
                print(f"    {s_idx+1}/{S} seeds done")

        for arm in arms:
            if arm == 'coupled':
                dxa, uexa = dx, u_exact
            else:
                M = int(arm.replace('fixed', ''))
                dxa, uexa = fixed[M]['dx'], fixed[M]['u_exact']
            E_bias, E_spread, E_total, ident = _decompose(
                np.array(runs[arm]), uexa, dxa)
            results[arm].append({
                'N': N, 'S': len(runs[arm]),
                'E_bias': E_bias, 'E_spread': E_spread, 'E_total': E_total,
                'identity_residual': ident,
            })
            print(f"    {arm:9s} E_bias={E_bias:.5f}  E_spread={E_spread:.5f}  "
                  f"E_total={E_total:.5f}")

    # ---- Fit slopes per arm ----
    fits = {}
    for arm in arms:
        N_arr  = np.array([r['N']        for r in results[arm]], dtype=float)
        bias_a = np.array([r['E_bias']   for r in results[arm]])
        spr_a  = np.array([r['E_spread'] for r in results[arm]])
        tot_a  = np.array([r['E_total']  for r in results[arm]])

        # The coupled arm reuses t4's exact rng sequence (303) and fit order so
        # its slopes AND bootstrap CIs reproduce the published ensemble values.
        rng_ci = np.random.default_rng(303 if arm == 'coupled' else
                                       310 + FIXED_BINS.index(int(arm[5:])))
        tot_lo, tot_hi = _bootstrap_slope_ci(N_arr, tot_a, rng=rng_ci)
        tot_fit = np.polyfit(np.log10(N_arr[tot_a > 0]), np.log10(tot_a[tot_a > 0]), 1)
        tot_slope = float(tot_fit[0])

        valid_b = bias_a > 1e-8
        if valid_b.sum() >= 2:
            bias_fit = np.polyfit(np.log10(N_arr[valid_b]), np.log10(bias_a[valid_b]), 1)
            bias_slope = float(bias_fit[0])
            bias_lo, bias_hi = _bootstrap_slope_ci(N_arr[valid_b], bias_a[valid_b], rng=rng_ci)
        else:
            bias_slope = float('nan'); bias_lo = float('nan'); bias_hi = float('nan')

        spr_fit = np.polyfit(np.log10(N_arr), np.log10(np.maximum(spr_a, 1e-12)), 1)
        spr_slope = float(spr_fit[0])
        spr_lo, spr_hi = _bootstrap_slope_ci(N_arr, spr_a, rng=rng_ci)

        fits[arm] = {
            'total_slope': tot_slope, 'total_ci': [tot_lo, tot_hi],
            'bias_slope': bias_slope, 'bias_ci': [bias_lo, bias_hi],
            'spread_slope': spr_slope, 'spread_ci': [spr_lo, spr_hi],
        }
        print(f"\n  [{arm}] slopes: total={tot_slope:.3f}  bias={bias_slope:.3f}  "
              f"spread={spr_slope:.3f}")

    # ---- Floor diagnostics for the fixed arms ----
    N_max = max(N_seq)
    floor = {}
    for M in FIXED_BINS:
        arm = f'fixed{M}'
        row_max = [r for r in results[arm] if r['N'] == N_max][0]
        # crossover: first N at which the deterministic component exceeds spread
        crossover = None
        for r in results[arm]:
            if r['E_bias'] > r['E_spread']:
                crossover = r['N']
                break
        floor[arm] = {
            'E_total_at_Nmax': row_max['E_total'],
            'E_bias_at_Nmax': row_max['E_bias'],
            'E_spread_at_Nmax': row_max['E_spread'],
            'bias_crossover_N': crossover,
            'bin_width': L / M,
        }
    r300 = floor['fixed300']['E_bias_at_Nmax']
    r400 = floor['fixed400']['E_bias_at_Nmax']
    floor['bias_floor_ratio_300_over_400'] = r300 / r400 if r400 > 0 else float('nan')
    floor['bin_width_ratio_400_over_300'] = 400.0 / 300.0
    print(f"\n  Fixed-grid E_bias floors at N={N_max}: "
          f"300-bin={r300:.5f}  400-bin={r400:.5f}  "
          f"ratio={floor['bias_floor_ratio_300_over_400']:.3f} "
          f"(bin-width ratio {400/300:.3f})")

    # ---- Save outputs ----
    csv_path = _mk(OUT_BASE, 'summary_by_N_paired.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['arm', 'N', 'S', 'E_bias', 'E_spread', 'E_total', 'identity_residual'])
        for arm in arms:
            for r in results[arm]:
                w.writerow([arm, r['N'], r['S'],
                            f"{r['E_bias']:.8f}", f"{r['E_spread']:.8f}",
                            f"{r['E_total']:.8f}", f"{r['identity_residual']:.2e}"])

    full_results = {
        'params': {'alpha': alpha, 'T': T, 'L': L, 'x0': x0, 'uL': uL, 'uR': uR,
                   'dt': dt, 'S': S, 'N_seq': N_seq, 'fixed_bins': list(FIXED_BINS),
                   'seeds': f'{base_seed}..{base_seed+S-1}'},
        'design_note': (
            'One simulation per (N, seed); three reconstructions of the same '
            'particles. Coupled arm reproduces the heat ensemble study (t4); '
            'fixed arms apply the Paper-1 300/400-bin diagnostic treatment.'),
        'fits': fits,
        'floor': floor,
        'per_N': results,
    }
    with open(_mk(OUT_BASE, 'summary.json'), 'w') as f:
        json.dump(full_results, f, indent=2)

    # ---- Paper figure: two panels, no in-figure titles ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))

    N_arr = np.array(N_seq, dtype=float)
    styles = {'coupled':  dict(color='tab:green', marker='^', ls='-'),
              'fixed300': dict(color='tab:blue',  marker='s', ls='-'),
              'fixed400': dict(color='tab:purple', marker='d', ls='--')}
    labels = {'coupled': 'coupled grid ($M=N$)',
              'fixed300': 'fixed grid ($M=300$)',
              'fixed400': 'fixed grid ($M=400$)'}
    for arm in arms:
        tot_a = np.array([r['E_total'] for r in results[arm]])
        axL.loglog(N_arr, tot_a, lw=1.6, ms=5, label=labels[arm], **styles[arm])
    guide_ref = np.array([N_arr.min(), N_arr.max()])
    c0 = np.array([r['E_total'] for r in results['coupled']])[0] * N_arr[0]**0.5
    axL.loglog(guide_ref, c0 * guide_ref**(-0.5), 'k:', lw=1.2, label=r'$N^{-1/2}$ guide')
    axL.set_xlabel(r'$N$'); axL.set_ylabel(r'$E_{\mathrm{total}}$')
    axL.legend(fontsize=8); axL.grid(True, which='both', alpha=0.3)
    axL.text(0.02, 0.02, '(a)', transform=axL.transAxes, fontsize=11)

    arm = 'fixed300'
    bias_a = np.array([r['E_bias']   for r in results[arm]])
    spr_a  = np.array([r['E_spread'] for r in results[arm]])
    tot_a  = np.array([r['E_total']  for r in results[arm]])
    axR.loglog(N_arr, bias_a, color='tab:blue', marker='s', ls='-', lw=1.6, ms=5,
               label=r'$E_{\mathrm{bias}}$')
    axR.loglog(N_arr, spr_a, color='tab:red', marker='o', ls='-', lw=1.6, ms=5,
               label=r'$E_{\mathrm{spread}}$')
    axR.loglog(N_arr, tot_a, color='k', marker='^', ls='--', lw=1.2, ms=4,
               label=r'$E_{\mathrm{total}}$')
    axR.set_xlabel(r'$N$'); axR.set_ylabel('error (fixed 300-bin grid)')
    axR.legend(fontsize=8); axR.grid(True, which='both', alpha=0.3)
    axR.text(0.02, 0.02, '(b)', transform=axR.transAxes, fontsize=11)

    fig.tight_layout()
    _savefig(fig, _mk(OUT_BASE, 'heat_grid_paired'))

    print(f"\n  [Task 7] Done. Outputs in {OUT_BASE}/")
    return full_results


if __name__ == '__main__':
    run_task7()
