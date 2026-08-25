"""Paired heat output-grid study with aligned-evaluation correction and a
deterministic operator control (manuscript `sec:heat-grid-paired`,
`tab:heat-grid-paired`, `fig:heat-grid-paired`, arXiv:2608.22592; reproduce
target: t7).

Same physical and numerical parameters, particle counts, and 30 seeds as the
heat ensemble study (t4): alpha=0.5, domain [0,4], x0=2, increasing unit step,
T=0.5, dt=0.005, N in {500,...,50000}, seeds 42..71.

For each (N, seed) ONE simulation is run; the identical final particle set is
then reconstructed under five output-grid treatments:

  coupled   : sorted cumulative sum interpolated to a uniform N-point grid
              (exactly the t4 ensemble treatment; the coupled rows and the
              rng(303) fitted slopes/CIs reproduce t4's published values
              digit-for-digit)
  fixed300  : weights binned on a fixed 300-bin grid, cumulative sum, compared
              with the exact profile at bin CENTERS (the convention of the
              original Paper-1 diagnostic, retained under that label as the
              diagnosed convention)
  fixed400  : the same with 400 bins
  fixed300e : the same 300-bin reconstruction compared with the exact profile
              at bin RIGHT EDGES (the aligned-evaluation correction: the
              cumulative sum over bins 0..k is the reconstruction at the right
              edge of bin k, so that is where the reference is evaluated)
  fixed400e : the same with 400 bins

The center-compare and edge-compare treatments share the identical reconstruction
vector, so their stochastic spreads agree; only the deterministic comparison
convention differs.

Deterministic operator control (no particles): the exact reflected position
law on [0,4] (method of images) is passed through the same bin-and-sum
operator and compared under both conventions. This predicts the fixed-grid
bias floor of the diagnosed convention (0.00435 / 0.00334 at 300/400 bins)
and the grid-independent finite-domain residual of the aligned convention
(0.00109), with no stochastic input.

Uncertainty: realization-level bootstrap confidence intervals for the fitted
total/spread/bias slopes of every treatment, obtained by resampling the 30
realizations within each particle count (Gram-matrix formulation, n_boot
5000, rng seed 12345, matching the FHN study's procedure).

Output: output/final_prepublication_tests/heat_grid_paired/
"""
import csv
import json
import os
import sys
import time
from math import erf, sqrt

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import SimulationConfig
from verify_solver import exact_heat_step

OUT_BASE = 'output/final_prepublication_tests/heat_grid_paired'

FIXED_BINS = (300, 400)
N_BOOT = 5000
BOOT_SEED = 12345


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
    """Single Heat GRW run; identical to t4's runner (same RNG consumption)."""
    from simulation import simulate_heat_equation
    np.random.seed(seed)

    weight = (uR - uL) / N
    ic = [(x0, weight)] * N

    cfg = SimulationConfig(
        equation_type='heat',
        domain_type='Finite',
        domain_size=L,
        boundary_conditions={'LEFT': {'type': 'Dirichlet', 'value': float(uL)},
                             'RIGHT': {'type': 'Dirichlet', 'value': float(uR)}},
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
    u_sorted = float(uL) + np.cumsum(w_arr[order])
    return x_pos, w_arr, x_sorted, u_sorted, elapsed


def _bootstrap_slope_ci(x_arr, y_arr, n_boot=2000, ci=0.95, rng=None):
    """Design-point (level-resampling) CI; retained so the coupled treatment's
    rng(303) fit sequence reproduces t4's published values exactly."""
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
    u_mean = u_arr.mean(axis=0)
    E_bias   = float(np.sqrt(np.sum((u_mean - u_exact)**2 * dx)))
    E_spread = float(np.sqrt(np.mean(np.sum((u_arr - u_mean[None, :])**2 * dx, axis=1))))
    E_total  = float(np.sqrt(np.mean(np.sum((u_arr - u_exact[None, :])**2 * dx, axis=1))))
    E_total_check = float(np.sqrt(E_bias**2 + E_spread**2))
    identity_err  = abs(E_total - E_total_check) / max(E_total, 1e-15)
    return E_bias, E_spread, E_total, identity_err


def _gram(u_arr, u_exact, dx):
    """Gram matrix of realization deviations d_s = u_s - u_exact under the
    L_h^2 inner product; every bootstrap statistic follows from it."""
    d = u_arr - u_exact[None, :]
    return (d @ d.T) * dx


def _realization_boot_slopes(grams, N_seq, n_boot=N_BOOT, seed=BOOT_SEED):
    """Realization-level bootstrap CIs for the total/spread/bias log-log
    slopes: resample the S realizations within each N (multinomial counts),
    recompute each statistic from the Gram matrix, refit the slope."""
    rng = np.random.default_rng(seed)
    S = grams[0].shape[0]
    lN = np.log10(np.array(N_seq, dtype=float))
    out = {k: [] for k in ('total', 'spread', 'bias')}
    for _ in range(n_boot):
        et, es, eb = [], [], []
        for G in grams:
            idx = rng.integers(0, S, size=S)
            c = np.bincount(idx, minlength=S).astype(float)
            Et2 = float(c @ np.diag(G)) / S
            Eb2 = float(c @ G @ c) / S**2
            Es2 = max(Et2 - Eb2, 0.0)
            et.append(Et2); es.append(Es2); eb.append(Eb2)
        for key, vals in (('total', et), ('spread', es), ('bias', eb)):
            v = np.sqrt(np.maximum(np.array(vals), 1e-30))
            out[key].append(float(np.polyfit(lN, np.log10(v), 1)[0]))
    ci = {}
    for key, slopes in out.items():
        ci[key] = [float(np.percentile(slopes, 2.5)), float(np.percentile(slopes, 97.5))]
    return ci


def _operator_control(L, x0, alpha, T, fixed):
    """Exact reflected position law through the bin-and-sum operator.

    The reflecting walls give the position density by the method of images.
    Binning its exact per-bin masses and summing reproduces the reconstruction
    the particle treatments approximate, with no stochastic input.
    """
    s = sqrt(2 * alpha * T)
    erfv = np.vectorize(erf)

    def img_cdf(x):
        total = np.zeros_like(np.asarray(x, dtype=float))
        for n in range(-4, 5):
            for mu in (2*n*L + x0, 2*n*L - x0):
                total += 0.5*(1 + erfv((np.asarray(x, dtype=float) - mu)/(s*sqrt(2))))
        return total

    control = {}
    for M in FIXED_BINS:
        edges = fixed[M]['edges']
        mass = img_cdf(edges[1:]) - img_cdf(edges[:-1])
        u_det = np.cumsum(mass)
        dxM = fixed[M]['dx']
        err_center = float(np.sqrt(np.sum((u_det - fixed[M]['u_exact_center'])**2 * dxM)))
        err_edge   = float(np.sqrt(np.sum((u_det - fixed[M]['u_exact_edge'])**2 * dxM)))
        control[f'M{M}'] = {'center_compare': err_center, 'edge_compare': err_edge}
    # grid-independent finite-domain reference gap on a fine grid
    Mf = 20000
    edges_f = np.linspace(0.0, L, Mf + 1)
    mass_f = img_cdf(edges_f[1:]) - img_cdf(edges_f[:-1])
    u_f = np.cumsum(mass_f)
    dxf = float(edges_f[1] - edges_f[0])
    uex_f = exact_heat_step(edges_f[1:], T, x0, 0.0, 1.0, alpha)
    control['finite_domain_reference_gap'] = float(np.sqrt(np.sum((u_f - uex_f)**2 * dxf)))
    return control


def run_task7(N_seq=None, S=30, base_seed=42,
              alpha=0.5, T=0.5, L=4.0, x0=2.0, uL=0.0, uR=1.0,
              dt=0.005):
    os.makedirs(OUT_BASE, exist_ok=True)
    if N_seq is None:
        N_seq = [500, 1000, 2000, 5000, 10000, 20000, 50000]

    treatments = ['coupled'] + [f'fixed{M}' for M in FIXED_BINS] \
                       + [f'fixed{M}e' for M in FIXED_BINS]

    print(f"{'='*60}\n  Task 7 (v2): Paired heat output-grid study\n"
          f"  N_seq={N_seq}  S={S}  treatments={treatments}\n{'='*60}")

    fixed = {}
    for M in FIXED_BINS:
        edges   = np.linspace(0.0, L, M + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        dxM     = float(edges[1] - edges[0])
        fixed[M] = {
            'edges': edges, 'centers': centers, 'dx': dxM,
            'u_exact_center': exact_heat_step(centers, T, x0, uL, uR, alpha),
            'u_exact_edge':   exact_heat_step(edges[1:], T, x0, uL, uR, alpha),
        }

    control = _operator_control(L, x0, alpha, T, fixed)
    print("  Operator control (deterministic, no particles):")
    for M in FIXED_BINS:
        c = control[f'M{M}']
        print(f"    M={M}: center-compare={c['center_compare']:.5f}  "
              f"edge-compare={c['edge_compare']:.5f}")
    print(f"    finite-domain reference gap (fine grid): "
          f"{control['finite_domain_reference_gap']:.5f}")

    results = {treatment: [] for treatment in treatments}
    grams   = {treatment: [] for treatment in treatments}

    for N in N_seq:
        x_grid = np.linspace(0.0, L, N)
        u_exact = exact_heat_step(x_grid, T, x0, uL, uR, alpha)
        dx = float(x_grid[1] - x_grid[0])

        print(f"\n  N={N}")
        seeds = [base_seed + i for i in range(S)]
        runs = {treatment: [] for treatment in treatments}

        for s_idx, seed in enumerate(seeds):
            x_pos, w_arr, x_out, u_out, elapsed = _run_heat_one(
                N, alpha, T, dt, L, x0, uL, uR, seed)

            if len(x_out) != N or not np.allclose(x_out, x_grid, atol=1e-10):
                u_coupled = np.interp(x_grid, x_out, u_out,
                                      left=float(uL), right=float(uR))
            else:
                u_coupled = u_out
            runs['coupled'].append(u_coupled)

            for M in FIXED_BINS:
                bin_w, _ = np.histogram(x_pos, bins=fixed[M]['edges'], weights=w_arr)
                u_fixed = float(uL) + np.cumsum(bin_w)
                runs[f'fixed{M}'].append(u_fixed)
                runs[f'fixed{M}e'].append(u_fixed)  # same vector, different reference

            if (s_idx + 1) % 10 == 0:
                print(f"    {s_idx+1}/{S} seeds done")

        for treatment in treatments:
            if treatment == 'coupled':
                dxa, uexa = dx, u_exact
            else:
                M = int(treatment.replace('fixed', '').replace('e', ''))
                dxa = fixed[M]['dx']
                uexa = fixed[M]['u_exact_edge'] if treatment.endswith('e') \
                    else fixed[M]['u_exact_center']
            u_arr = np.array(runs[treatment])
            E_bias, E_spread, E_total, ident = _decompose(u_arr, uexa, dxa)
            grams[treatment].append(_gram(u_arr, uexa, dxa))
            results[treatment].append({
                'N': N, 'S': len(runs[treatment]),
                'E_bias': E_bias, 'E_spread': E_spread, 'E_total': E_total,
                'identity_residual': ident,
            })
            print(f"    {treatment:10s} E_bias={E_bias:.5f}  E_spread={E_spread:.5f}  "
                  f"E_total={E_total:.5f}")

    # ---- Fits ----
    fits = {}
    for treatment in treatments:
        N_arr  = np.array([r['N']        for r in results[treatment]], dtype=float)
        bias_a = np.array([r['E_bias']   for r in results[treatment]])
        spr_a  = np.array([r['E_spread'] for r in results[treatment]])
        tot_a  = np.array([r['E_total']  for r in results[treatment]])

        # coupled treatment reuses t4's exact rng sequence so its slopes AND
        # design-point CIs reproduce the published ensemble values
        rng_ci = np.random.default_rng(303 if treatment == 'coupled' else 400 + treatments.index(treatment))
        tot_lo, tot_hi = _bootstrap_slope_ci(N_arr, tot_a, rng=rng_ci)
        tot_slope = float(np.polyfit(np.log10(N_arr[tot_a > 0]),
                                     np.log10(tot_a[tot_a > 0]), 1)[0])
        valid_b = bias_a > 1e-8
        if valid_b.sum() >= 2:
            bias_slope = float(np.polyfit(np.log10(N_arr[valid_b]),
                                          np.log10(bias_a[valid_b]), 1)[0])
            bias_lo, bias_hi = _bootstrap_slope_ci(N_arr[valid_b], bias_a[valid_b], rng=rng_ci)
        else:
            bias_slope = float('nan'); bias_lo = float('nan'); bias_hi = float('nan')
        spr_slope = float(np.polyfit(np.log10(N_arr),
                                     np.log10(np.maximum(spr_a, 1e-12)), 1)[0])
        spr_lo, spr_hi = _bootstrap_slope_ci(N_arr, spr_a, rng=rng_ci)

        real_ci = _realization_boot_slopes(grams[treatment], N_seq)

        fits[treatment] = {
            'total_slope': tot_slope, 'total_ci_design': [tot_lo, tot_hi],
            'bias_slope': bias_slope, 'bias_ci_design': [bias_lo, bias_hi],
            'spread_slope': spr_slope, 'spread_ci_design': [spr_lo, spr_hi],
            'total_ci_realization': real_ci['total'],
            'spread_ci_realization': real_ci['spread'],
            'bias_ci_realization': real_ci['bias'],
        }
        print(f"\n  [{treatment}] slopes: total={tot_slope:.3f} "
              f"CI_real={real_ci['total']}  spread={spr_slope:.3f} "
              f"CI_real={real_ci['spread']}  bias={bias_slope:.3f}")

    # ---- Floor diagnostics ----
    N_max = max(N_seq)
    floor = {}
    for treatment in treatments:
        if treatment == 'coupled':
            continue
        row_max = [r for r in results[treatment] if r['N'] == N_max][0]
        crossover = None
        for r in results[treatment]:
            if r['E_bias'] > r['E_spread']:
                crossover = r['N']
                break
        floor[treatment] = {
            'E_total_at_Nmax': row_max['E_total'],
            'E_bias_at_Nmax': row_max['E_bias'],
            'E_spread_at_Nmax': row_max['E_spread'],
            'bias_crossover_N': crossover,
        }
    print("\n  Bias at N=%d by treatment:" % N_max)
    for treatment, f in floor.items():
        print(f"    {treatment:10s} E_bias={f['E_bias_at_Nmax']:.5f}  "
              f"crossover_N={f['bias_crossover_N']}")

    # ---- Save ----
    with open(_mk(OUT_BASE, 'summary_by_N_paired.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['treatment', 'N', 'S', 'E_bias', 'E_spread', 'E_total', 'identity_residual'])
        for treatment in treatments:
            for r in results[treatment]:
                w.writerow([treatment, r['N'], r['S'],
                            f"{r['E_bias']:.8f}", f"{r['E_spread']:.8f}",
                            f"{r['E_total']:.8f}", f"{r['identity_residual']:.2e}"])

    full_results = {
        'params': {'alpha': alpha, 'T': T, 'L': L, 'x0': x0, 'uL': uL, 'uR': uR,
                   'dt': dt, 'S': S, 'N_seq': N_seq, 'fixed_bins': list(FIXED_BINS),
                   'seeds': f'{base_seed}..{base_seed+S-1}',
                   'n_boot_realization': N_BOOT, 'boot_seed': BOOT_SEED},
        'design_note': (
            'One simulation per (N, seed); five reconstructions of the same '
            'particles. coupled reproduces t4. fixedM = bin-and-sum compared '
            'at bin centers (diagnosed convention). fixedMe = identical '
            'reconstruction compared at bin right edges (aligned-evaluation '
            'correction). The operator control passes the exact reflected '
            'position law through the same operator.'),
        'operator_control': control,
        'fits': fits,
        'floor': floor,
        'per_N': results,
    }
    with open(_mk(OUT_BASE, 'summary.json'), 'w') as f:
        json.dump(full_results, f, indent=2)

    # ---- Figure: (a) totals, (b) bias vs operator-control predictions ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    N_arr = np.array(N_seq, dtype=float)
    styles = {'coupled':   dict(color='tab:green',  marker='^', ls='-'),
              'fixed300':  dict(color='tab:blue',   marker='s', ls='-'),
              'fixed300e': dict(color='tab:cyan',   marker='o', ls='--'),
              'fixed400':  dict(color='tab:purple', marker='d', ls='-'),
              'fixed400e': dict(color='tab:pink',   marker='v', ls='--')}
    labels = {'coupled':   'coupled grid ($M=N$)',
              'fixed300':  'fixed $M=300$, bin center',
              'fixed300e': 'fixed $M=300$, right edge',
              'fixed400':  'fixed $M=400$, bin center',
              'fixed400e': 'fixed $M=400$, right edge'}
    for treatment in ('coupled', 'fixed300', 'fixed300e'):
        tot_a = np.array([r['E_total'] for r in results[treatment]])
        axL.loglog(N_arr, tot_a, lw=1.6, ms=5, label=labels[treatment], **styles[treatment])
    guide_ref = np.array([N_arr.min(), N_arr.max()])
    c0 = np.array([r['E_total'] for r in results['coupled']])[0] * N_arr[0]**0.5
    axL.loglog(guide_ref, c0 * guide_ref**(-0.5), 'k:', lw=1.2, label=r'$N^{-1/2}$ guide')
    axL.set_xlabel(r'$N$'); axL.set_ylabel(r'$E_{\mathrm{total}}$')
    axL.legend(fontsize=8); axL.grid(True, which='both', alpha=0.3)
    axL.text(0.02, 0.02, '(a)', transform=axL.transAxes, fontsize=11)

    for treatment in ('fixed300', 'fixed400', 'fixed300e', 'fixed400e', 'coupled'):
        bias_a = np.array([r['E_bias'] for r in results[treatment]])
        axR.loglog(N_arr, bias_a, lw=1.4, ms=4, label=labels[treatment], **styles[treatment])
    axR.axhline(control['M300']['center_compare'], color='tab:blue', ls=':', lw=1.1,
                label='deterministic, bin center')
    axR.axhline(control['M400']['center_compare'], color='tab:purple', ls=':', lw=1.1)
    axR.axhline(control['M300']['edge_compare'], color='k', ls=':', lw=1.1,
                label='deterministic, right edge')
    axR.set_xlabel(r'$N$'); axR.set_ylabel(r'$E_{\mathrm{bias}}$')
    axR.legend(fontsize=7); axR.grid(True, which='both', alpha=0.3)
    axR.text(0.02, 0.02, '(b)', transform=axR.transAxes, fontsize=11)

    fig.tight_layout()
    _savefig(fig, _mk(OUT_BASE, 'heat_grid_paired'))

    print(f"\n  [Task 7 v2] Done. Outputs in {OUT_BASE}/")
    return full_results


if __name__ == '__main__':
    run_task7()
