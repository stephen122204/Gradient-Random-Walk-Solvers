"""Controlled attribution for the Cole-Hopf Burgers pipeline (manuscript
`sec:colehopf-diagnosis`, `tab:burgers-decoupled-design`,
`fig:burgers-decoupled`, `fig:burgers-boundary-domain`,
`fig:burgers-perturbation-response`, arXiv:2608.22592; reproduce target: t8).

The paper's coupled refinement varies three numerical choices together: the number of
phi_x globs, the output-grid size (M = N), and the physical smoothing
bandwidth (the kernel standard deviation is fixed at 12 bins, so its physical
width shrinks with the grid). This study decouples them, adds the missing
deterministic boundary control, and turns the one-point perturbation
experiment into a response curve.

Parts
  V  validation: the parameterized pipeline reproduces the packaged solver
     bit-for-bit at the paper configuration (P=M=400, sigma_bins=12, seed 42)
  A  particle refinement at fixed output grid and fixed physical bandwidth
  B  output-grid refinement at fixed glob count and fixed physical bandwidth
  C  bandwidth sweep at fixed glob count and output grid
  D  deterministic boundary controls (no particles), for L in {4,6,8,10}:
     pinned-endpoint FD solve (the pipeline's boundary model) versus an FD
     solve with the exact time-scaled Dirichlet data phi0 * exp(nu k^2 t),
     k = A/(2 nu), which the exact transformed solution satisfies; plus the
     exact-phi-through-recovery control
  E  perturbation response curve: additive noise of several amplitudes on the
     exact normalized transformed field, many realizations each, with the
     measured particle reconstruction error of phi marked for comparison
  F  multi-seed domain study: the paper's N = 100 L design at L in {4,6,8,10}
     with 30 seeds per domain (E_det is deterministic; E_GRW and E_total get
     ensemble means and spreads)

Error conventions: RMSE for the u-space decomposition (the paper's Burgers
convention) with the L_h^2 value also recorded (norms differ by sqrt(h M)).
phi-space errors are RMSE against the pinned-endpoint FD reference, which is
the boundary model the pipeline itself enforces.

Output: output/final_prepublication_tests/burgers_controls/
"""
import contextlib
import csv
import io
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
from simulation import (
    simulate_burgers_cole_hopf_grw,
    _reference_phi_heat_fd,
    _reflect_arrays,
    _validated_step_count,
)

OUT_BASE = 'output/final_prepublication_tests/burgers_controls'

A_AMP = 1.0
NU = 0.5
DT = 0.005
T_END = 0.5
S_DECOUPLED = 20     # seeds for the three decoupled one-factor sweeps
S_DOMAIN = 30        # seeds for the multi-seed domain study
BASE_SEED = 42


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


def _rmse(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


# ---------------------------------------------------------------- pipeline --
def _phi0_from_grid(x, u0, nu):
    """Trapezoidal Psi0 and max-normalized phi0, exactly as the solver builds
    them."""
    inc = 0.5 * (u0[:-1] + u0[1:]) * np.diff(x)
    Psi0 = np.concatenate([[0.0], np.cumsum(inc)])
    log_phi0 = -Psi0 / (2.0 * nu)
    log_phi0 -= log_phi0.max()
    phi0 = np.exp(np.clip(log_phi0, -700.0, 0.0))
    return phi0


def pipeline(P, M, sigma_bins, L, seed, nu=NU, A=A_AMP, dt=DT, T=T_END):
    """Parameterized Cole-Hopf particle pipeline.

    P: number of initialization grid points (P-1 phi_x globs at midpoints)
    M: number of output-grid points (binning, smoothing, recovery)
    sigma_bins: Gaussian kernel standard deviation in units of OUTPUT bins
                (may be fractional); the physical bandwidth is
                sigma_x = sigma_bins * L / (M - 1)

    At P == M and sigma_bins == 12 this reproduces the packaged solver
    operation-for-operation, including the RNG stream.
    """
    xc = L / 2.0
    x_init = np.linspace(0.0, L, P)
    u0 = -A * np.tanh(A * (x_init - xc) / (2.0 * nu))
    phi0 = _phi0_from_grid(x_init, u0, nu)
    phi0_0, phi0_L = float(phi0[0]), float(phi0[-1])
    exact_integral = phi0_L - phi0_0

    w_ph = np.diff(phi0)
    x_ph = 0.5 * (x_init[:-1] + x_init[1:]).copy()

    np.random.seed(seed)
    sigma_step = np.sqrt(2.0 * nu * dt)
    n_steps = _validated_step_count(T, dt)
    for _ in range(n_steps):
        x_ph = x_ph + np.random.normal(0.0, sigma_step, size=x_ph.shape)
        x_ph, w_ph = _reflect_arrays(
            x_ph, w_ph, L, "dirichlet", "dirichlet"
        )

    x_out = np.linspace(0.0, L, M)
    dx_out = float(x_out[1] - x_out[0])
    bin_sums = np.zeros(M)
    idx = np.clip(np.floor(x_ph / dx_out).astype(int), 0, M - 1)
    np.add.at(bin_sums, idx, w_ph)

    kw = int(4 * sigma_bins) + 1
    kernel_x = np.arange(-kw, kw + 1, dtype=float)
    kernel = np.exp(-0.5 * (kernel_x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    raw_conv = np.convolve(bin_sums, kernel, mode='same')
    kernel_norm = np.convolve(np.ones(M, dtype=float), kernel, mode='same')
    bin_sums_s = raw_conv / np.maximum(kernel_norm, 1e-12)

    near_zero = 1e-6 * max(float(np.abs(bin_sums).max()), 1e-30)
    if abs(exact_integral) < near_zero:
        bin_sums_s -= bin_sums_s.mean()
    elif abs(bin_sums_s.sum()) > 1e-30:
        bin_sums_s *= exact_integral / bin_sums_s.sum()

    phi_out = phi0_0 + np.cumsum(bin_sums_s)
    phi_x_out = np.gradient(phi_out, dx_out)

    phi0_min = float(phi0.min())
    phi_floor = max(phi0_min / 2.0, 1e-10)
    clipped = phi_out < phi_floor
    phi_safe = np.where(clipped, phi_floor, phi_out)
    u_out = -2.0 * nu * phi_x_out / phi_safe
    u_out = np.where(clipped, 0.0, u_out)

    return {
        'x_out': x_out, 'u': u_out, 'phi': phi_out, 'phi_x': phi_x_out,
        'clip_count': int(clipped.sum()), 'min_phi': float(phi_out.min()),
        'endpoint_residual': float(abs(phi_out[-1] - phi0_L)),
        'sigma_x': sigma_bins * dx_out,
    }


def _fd_reference(M, L, nu=NU, A=A_AMP, T=T_END, boundary='pinned'):
    """Deterministic FD solve of phi_t = nu phi_xx on the M-point grid.

    boundary='pinned': endpoints held at their initial values (the pipeline's
    boundary model, identical to _reference_phi_heat_fd).
    boundary='exact' : endpoints follow phi0 * exp(nu k^2 t), k = A/(2 nu),
    the values the exact transformed solution takes on the domain edge.
    """
    xc = L / 2.0
    x_out = np.linspace(0.0, L, M)
    u0 = -A * np.tanh(A * (x_out - xc) / (2.0 * nu))
    phi0 = _phi0_from_grid(x_out, u0, nu)
    phi0_0, phi0_L = float(phi0[0]), float(phi0[-1])
    if boundary == 'pinned':
        phi = _reference_phi_heat_fd(phi0.copy(), x_out, nu, T, phi0_0, phi0_L)
    else:
        k = A / (2.0 * nu)
        dx = float(x_out[1] - x_out[0])
        dt_fd = 0.4 * dx**2 / nu
        n_steps = int(np.ceil(T / dt_fd))
        dt_fd = T / n_steps
        r = nu * dt_fd / dx**2
        phi = phi0.copy()
        for n in range(n_steps):
            t_new = (n + 1) * dt_fd
            growth = float(np.exp(nu * k**2 * t_new))
            phi_new = phi.copy()
            phi_new[1:-1] = phi[1:-1] + r * (phi[2:] - 2.0 * phi[1:-1] + phi[:-2])
            phi_new[0] = phi0_0 * growth
            phi_new[-1] = phi0_L * growth
            phi = phi_new
    dx = float(x_out[1] - x_out[0])
    phi_x = np.gradient(phi, dx)
    phi_safe = np.where(np.abs(phi) < 1e-12, 1e-12, phi)
    u = -2.0 * nu * phi_x / phi_safe
    u_exact = -A * np.tanh(A * (x_out - xc) / (2.0 * nu))
    return {'x': x_out, 'phi': phi, 'u': u, 'u_exact': u_exact,
            'rmse_u': _rmse(u, u_exact),
            'l2_u': float(np.sqrt(np.sum((u - u_exact)**2 * dx)))}


# -------------------------------------------------------------- validation --
def validate_pipeline():
    L, M = 4.0, 400
    xc = L / 2.0
    x_out = np.linspace(0.0, L, M)
    u0 = -A_AMP * np.tanh(A_AMP * (x_out - xc) / (2.0 * NU))
    globs = [{'position': float(x_out[i]), 'value': [float(u0[i])]}
             for i in range(M)]
    cfg = SimulationConfig(
        equation_type='burgers', domain_type='finite', domain_size=L,
        boundary_conditions={'LEFT': {'type': 'dirichlet', 'value': 0},
                             'RIGHT': {'type': 'dirichlet', 'value': 0}},
        diff_constant=NU, time_step=DT, total_time=T_END, num_points=M,
        initial_conditions=list(zip(x_out.tolist(), u0.tolist())),
        reaction_term=False, burgers_mode='cole_hopf_grw',
        burgers_ic_type='stationary_shock', burgers_ic_amplitude=A_AMP,
    )
    np.random.seed(42)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        packaged = simulate_burgers_cole_hopf_grw(globs, cfg)
    u_packaged = np.array([float(g['value'][0]) for g in packaged])
    mine = pipeline(P=M, M=M, sigma_bins=12, L=L, seed=42)
    max_dev = float(np.max(np.abs(mine['u'] - u_packaged)))
    print(f"  [V] pipeline vs packaged solver: max|du| = {max_dev:.3e}")
    assert max_dev < 1e-12, 'parameterized pipeline does not match the solver'
    return max_dev


# ----------------------------------------------- decoupled one-factor sweeps --
def _ensemble(configs, part_name, fd_refs, per_seed_store=None,
              per_run_rows=None):
    rows = []
    for tag, (P, M, sigma_bins, L) in configs.items():
        fd = fd_refs[(M, L)]
        u_fd, u_exact = fd['u'], fd['u_exact']
        phi_fd = fd['phi']
        dx = float(fd['x'][1] - fd['x'][0])
        per = {k: [] for k in ('rmse_total', 'rmse_particle', 'l2_total',
                               'rmse_phi', 'rmse_phi_x', 'rmse_ratio',
                               'clip', 'min_phi', 'endpoint_residual')}
        u_runs = []
        for s in range(S_DECOUPLED):
            seed = BASE_SEED + s
            r = pipeline(P=P, M=M, sigma_bins=sigma_bins, L=L,
                         seed=seed)
            u_runs.append(r['u'])
            per['rmse_total'].append(_rmse(r['u'], u_exact))
            per['rmse_particle'].append(_rmse(r['u'], u_fd))
            per['l2_total'].append(
                float(np.sqrt(np.sum((r['u'] - u_exact)**2 * dx))))
            per['rmse_phi'].append(_rmse(r['phi'], phi_fd))
            per['rmse_phi_x'].append(
                _rmse(r['phi_x'], np.gradient(phi_fd, dx)))
            ratio_n = r['phi_x'] / np.where(np.abs(r['phi']) < 1e-12, 1e-12, r['phi'])
            ratio_f = np.gradient(phi_fd, dx) / np.where(
                np.abs(phi_fd) < 1e-12, 1e-12, phi_fd)
            per['rmse_ratio'].append(_rmse(ratio_n, ratio_f))
            per['clip'].append(r['clip_count'])
            per['min_phi'].append(r['min_phi'])
            per['endpoint_residual'].append(r['endpoint_residual'])
            if per_run_rows is not None:
                per_run_rows.append({
                    'part': part_name, 'tag': tag, 'seed': seed,
                    'P': P, 'M': M, 'sigma_bins': sigma_bins,
                    'sigma_x': sigma_bins * L / (M - 1), 'L': L,
                    **{key: values[-1] for key, values in per.items()},
                })
        # u-space decomposition against the deterministic reference: the mean
        # deviation profile is the smoothing/reconstruction bias, and the
        # spread about it is the sampling part (RMSE convention).
        u_arr = np.array(u_runs)
        dev = u_arr - u_fd[None, :]
        dev_mean = dev.mean(axis=0)
        bias_u = float(np.sqrt(np.mean(dev_mean**2)))
        spread_u = float(np.sqrt(np.mean(np.mean((dev - dev_mean[None, :])**2,
                                                 axis=1))))
        row = {'part': part_name, 'tag': tag, 'P': P, 'M': M,
               'sigma_bins': sigma_bins,
               'sigma_x': sigma_bins * L / (M - 1), 'L': L,
               'S': S_DECOUPLED,
               'particle_bias_rmse': bias_u, 'particle_spread_rmse': spread_u}
        for k, v in per.items():
            row[k + '_mean'] = float(np.mean(v))
            row[k + '_std'] = float(np.std(v))
        rows.append(row)
        if per_seed_store is not None:
            per_seed_store[(part_name, tag)] = list(per['rmse_particle'])
        print(f"  [{part_name}] {tag}: u_total={row['rmse_total_mean']:.4f}"
              f"±{row['rmse_total_std']:.4f}  u_particle="
              f"{row['rmse_particle_mean']:.4f} (bias {bias_u:.4f}/spread "
              f"{spread_u:.4f})  phi={row['rmse_phi_mean']:.2e}")
    return rows


def run_decoupled_controls():
    L = 4.0
    sigma_x_ref = 12 * L / 399          # the paper configuration's bandwidth
    fd_refs = {}
    for M in (100, 200, 400, 800, 1600):
        fd_refs[(M, L)] = _fd_reference(M, L, boundary='pinned')

    per_seed = {}
    per_run_rows = []
    partA = {f'P{P}': (P, 400, sigma_x_ref * 399 / L, L)
             for P in (100, 200, 400, 800, 1600, 3200)}
    print('\n  Part A: initialization-grid refinement at fixed M=400, fixed sigma_x')
    rows = _ensemble(partA, 'A', fd_refs, per_seed_store=per_seed,
                     per_run_rows=per_run_rows)

    print('\n  Part B: output-grid refinement at fixed P=400, fixed sigma_x')
    partB = {f'M{M}': (400, M, sigma_x_ref * (M - 1) / L, L)
             for M in (100, 200, 400, 800, 1600)}
    rows += _ensemble(partB, 'B', fd_refs, per_run_rows=per_run_rows)

    print('\n  Part C: bandwidth sweep at fixed P=M=400')
    partC = {f's{sx}': (400, 400, sx * 399 / L, L)
             for sx in (0.03, 0.06, 0.12, 0.24, 0.48)}
    rows += _ensemble(partC, 'C', fd_refs, per_run_rows=per_run_rows)

    # realization-level bootstrap CI for the part-A particle-error slope:
    # resample the 20 per-seed particle errors within each P, recompute the
    # mean at each level, refit the log-log slope
    rng = np.random.default_rng(12345)
    P_vals = sorted({r['P'] for r in rows if r['part'] == 'A'})
    seed_lists = [np.array(per_seed[('A', f'P{P}')]) for P in P_vals]
    lP = np.log10(np.array(P_vals, dtype=float))
    slopes = []
    for _ in range(5000):
        means = [float(np.mean(v[rng.integers(0, len(v), size=len(v))]))
                 for v in seed_lists]
        slopes.append(float(np.polyfit(lP, np.log10(means), 1)[0]))
    partA_ci = [float(np.percentile(slopes, 2.5)),
                float(np.percentile(slopes, 97.5))]
    print(f"  [A] particle-error slope realization-level 95% CI: "
          f"[{partA_ci[0]:.3f}, {partA_ci[1]:.3f}]")
    return rows, partA_ci, per_run_rows


# ---------------------------------------------------- boundary controls (D) --
def run_boundary_controls():
    out = []
    for L in (4.0, 6.0, 8.0, 10.0):
        M = int(100 * L)
        pinned = _fd_reference(M, L, boundary='pinned')
        exact_bc = _fd_reference(M, L, boundary='exact')
        # exact-phi-through-recovery control (no FD solve, no particles)
        xc = L / 2.0
        x = pinned['x']; dx = float(x[1] - x[0])
        k = A_AMP / (2.0 * NU)
        phi_ex = np.cosh(k * (x - xc)) / np.cosh(k * xc)
        u_from_exact_phi = -2.0 * NU * np.gradient(phi_ex, dx) / phi_ex
        u_exact = pinned['u_exact']
        row = {
            'L': L, 'M': M,
            'E_det_pinned_rmse': pinned['rmse_u'],
            'E_det_pinned_l2': pinned['l2_u'],
            'E_det_exactBC_rmse': exact_bc['rmse_u'],
            'E_det_exactBC_l2': exact_bc['l2_u'],
            'exact_phi_recovery_rmse': _rmse(u_from_exact_phi, u_exact),
            'exact_phi_recovery_l2': float(
                np.sqrt(np.sum((u_from_exact_phi - u_exact)**2 * dx))),
        }
        out.append(row)
        print(f"  [D] L={L}: E_det pinned={row['E_det_pinned_rmse']:.4f} RMSE, "
              f"exact-BC={row['E_det_exactBC_rmse']:.2e} RMSE, "
              f"exact-phi recovery={row['exact_phi_recovery_rmse']:.2e} RMSE")
    return out


# ------------------------------------------------- perturbation response (E) --
def run_perturbation_response(measured_phi_rmse, measured_u_l2):
    """Response of the recovery map to perturbations of the transformed field.

    Two perturbation families with the same root-mean-square amplitude:
    white (independent per grid point) and kernel-smoothed (white noise passed
    through the pipeline's own Gaussian kernel, then rescaled). The exact
    transformed-field shape is scaled to the pipeline normalization (unit
    endpoint maximum), so additive amplitudes are comparable to the measured
    reconstruction error. The particle
    reconstruction error is kernel-smooth, so the smoothed family is the one
    the measured pipeline point should land on; the white family shows how
    much stronger the amplification is for uncorrelated errors.
    """
    L, M = 4.0, 400
    xc = L / 2.0
    x = np.linspace(0.0, L, M)
    dx = float(x[1] - x[0])
    k = A_AMP / (2.0 * NU)
    phi_ex = np.cosh(k * (x - xc)) / np.cosh(k * xc)
    u_exact = -A_AMP * np.tanh(A_AMP * (x - xc) / (2.0 * NU))
    sigma_bins = 12
    kw = int(4 * sigma_bins) + 1
    kernel_x = np.arange(-kw, kw + 1, dtype=float)
    kernel = np.exp(-0.5 * (kernel_x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    kernel_norm = np.convolve(np.ones(M), kernel, mode='same')

    amps = (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 2e-2)
    R = 20
    rows = []
    for amp in amps:
        errs = {'white': [], 'smooth': []}
        rng = np.random.default_rng(9000 + int(round(-np.log10(amp) * 10)))
        for _ in range(R):
            w = rng.standard_normal(M)
            fams = {'white': w / max(float(np.sqrt(np.mean(w**2))), 1e-30)}
            ws = np.convolve(w, kernel, mode='same') / np.maximum(kernel_norm, 1e-12)
            fams['smooth'] = ws / max(float(np.sqrt(np.mean(ws**2))), 1e-30)
            for fam, noise in fams.items():
                phi_n = np.clip(phi_ex + amp * noise, 1e-10, None)
                u_n = -2.0 * NU * np.gradient(phi_n, dx) / phi_n
                errs[fam].append(
                    float(np.sqrt(np.sum((u_n - u_exact)**2 * dx))))
        rows.append({'amplitude': amp, 'R': R,
                     'white_l2_mean': float(np.mean(errs['white'])),
                     'white_l2_std': float(np.std(errs['white'])),
                     'smooth_l2_mean': float(np.mean(errs['smooth'])),
                     'smooth_l2_std': float(np.std(errs['smooth']))})
        print(f"  [E] amp={amp:.0e}: white L2 = {rows[-1]['white_l2_mean']:.4f}"
              f" ± {rows[-1]['white_l2_std']:.4f}   smooth L2 = "
              f"{rows[-1]['smooth_l2_mean']:.4f} ± {rows[-1]['smooth_l2_std']:.4f}")
    return {'rows': rows, 'measured_phi_rmse': measured_phi_rmse,
            'measured_u_l2': measured_u_l2}


# --------------------------------------------- multi-seed domain study (F) --
def run_domain_multiseed():
    out = []
    per_run_rows = []
    total_by_L = {}
    for L in (4.0, 6.0, 8.0, 10.0):
        N = int(100 * L)
        xc = L / 2.0
        x_out = np.linspace(0.0, L, N)
        dx = float(x_out[1] - x_out[0])
        u_exact = -A_AMP * np.tanh(A_AMP * (x_out - xc) / (2.0 * NU))
        pinned = _fd_reference(N, L, boundary='pinned')
        u_fd = pinned['u']
        E_det = _rmse(u_fd, u_exact)
        totals, particles = [], []
        for s in range(S_DOMAIN):
            seed = BASE_SEED + s
            r = pipeline(P=N, M=N, sigma_bins=12, L=L, seed=seed)
            total = _rmse(r['u'], u_exact)
            particle = _rmse(r['u'], u_fd)
            totals.append(total)
            particles.append(particle)
            per_run_rows.append({
                'L': L, 'N': N, 'seed': seed,
                'E_det_rmse': E_det,
                'E_grw_rmse': particle,
                'E_total_rmse': total,
            })
        row = {'L': L, 'N': N, 'S': S_DOMAIN,
               'E_det_rmse': E_det,
               'E_grw_mean': float(np.mean(particles)),
               'E_grw_std': float(np.std(particles)),
               'E_total_mean': float(np.mean(totals)),
               'E_total_std': float(np.std(totals))}
        out.append(row)
        total_by_L[L] = np.asarray(totals, dtype=float)
        print(f"  [F] L={L}: E_det={E_det:.4f}  "
              f"E_GRW={row['E_grw_mean']:.4f}±{row['E_grw_std']:.4f}  "
              f"E_total={row['E_total_mean']:.4f}±{row['E_total_std']:.4f}")
    # The paper describes the total as statistically compatible with a flat
    # trend across the four tested domains. Quantify that statement directly:
    # independently resample realizations within each L, refit the linear
    # slope dE_total/dL, and form a percentile interval. The repeated seed
    # identifiers are not treated as common random numbers because N changes.
    L_vals = np.asarray(sorted(total_by_L), dtype=float)
    means = np.asarray([total_by_L[L].mean() for L in L_vals])
    slope = float(np.polyfit(L_vals, means, 1)[0])
    rng = np.random.default_rng(24680)
    slopes = []
    for _ in range(5000):
        boot_means = []
        for L in L_vals:
            values = total_by_L[float(L)]
            draw = values[rng.integers(0, len(values), size=len(values))]
            boot_means.append(float(draw.mean()))
        slopes.append(float(np.polyfit(L_vals, boot_means, 1)[0]))
    trend = {
        'metric': 'linear slope of ensemble-mean E_total RMSE versus L',
        'slope_per_unit_L': slope,
        'slope_ci_realization': [float(np.percentile(slopes, 2.5)),
                                 float(np.percentile(slopes, 97.5))],
        'n_boot': 5000,
        'boot_seed': 24680,
        'ci_contains_zero': bool(np.percentile(slopes, 2.5) <= 0.0
                                 <= np.percentile(slopes, 97.5)),
    }
    print(f"  [F] total-error trend dE/dL={slope:+.4e}, realization-level "
          f"95% CI=[{trend['slope_ci_realization'][0]:+.4e}, "
          f"{trend['slope_ci_realization'][1]:+.4e}]")
    return out, trend, per_run_rows


# ----------------------------------------------------------------- figures --
def make_figures(decoupled_rows, boundary_rows, pert, domain_rows):
    # decoupled one-factor sweeps: three panels
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    A_rows = [r for r in decoupled_rows if r['part'] == 'A']
    B_rows = [r for r in decoupled_rows if r['part'] == 'B']
    C_rows = [r for r in decoupled_rows if r['part'] == 'C']
    ax = axes[0]
    P = [r['P'] for r in A_rows]
    ax.errorbar(P, [r['rmse_total_mean'] for r in A_rows],
                yerr=[r['rmse_total_std'] for r in A_rows],
                fmt='o-', color='tab:blue', ms=4, lw=1.4, capsize=2,
                label='vs exact shock')
    ax.errorbar(P, [r['rmse_particle_mean'] for r in A_rows],
                yerr=[r['rmse_particle_std'] for r in A_rows],
                fmt='s--', color='tab:red', ms=4, lw=1.4, capsize=2,
                label='vs deterministic reference')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('initialization-grid points $P$'); ax.set_ylabel('RMSE')
    ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
    ax.text(0.03, 0.03, '(a)', transform=ax.transAxes, fontsize=11)
    ax = axes[1]
    Mv = [r['M'] for r in B_rows]
    ax.errorbar(Mv, [r['rmse_total_mean'] for r in B_rows],
                yerr=[r['rmse_total_std'] for r in B_rows],
                fmt='o-', color='tab:blue', ms=4, lw=1.4, capsize=2,
                label='vs exact shock')
    ax.errorbar(Mv, [r['rmse_particle_mean'] for r in B_rows],
                yerr=[r['rmse_particle_std'] for r in B_rows],
                fmt='s--', color='tab:red', ms=4, lw=1.4, capsize=2,
                label='vs deterministic reference')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('output-grid points $M$')
    ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
    ax.text(0.03, 0.03, '(b)', transform=ax.transAxes, fontsize=11)
    ax = axes[2]
    sx = [r['sigma_x'] for r in C_rows]
    ax.errorbar(sx, [r['rmse_total_mean'] for r in C_rows],
                yerr=[r['rmse_total_std'] for r in C_rows],
                fmt='o-', color='tab:blue', ms=4, lw=1.4, capsize=2,
                label='vs exact shock')
    ax.errorbar(sx, [r['rmse_particle_mean'] for r in C_rows],
                yerr=[r['rmse_particle_std'] for r in C_rows],
                fmt='s--', color='tab:red', ms=4, lw=1.4, capsize=2,
                label='vs deterministic reference')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'physical bandwidth $\sigma_x$')
    ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
    ax.text(0.03, 0.03, '(c)', transform=ax.transAxes, fontsize=11)
    fig.tight_layout()
    _savefig(fig, _mk(OUT_BASE, 'burgers_decoupled'))

    # boundary controls + multi-seed domain study
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.0))
    Lv = [r['L'] for r in boundary_rows]
    axL.semilogy(Lv, [r['E_det_pinned_rmse'] for r in boundary_rows],
                 's-', color='tab:red', ms=5, lw=1.6, label='pinned endpoints')
    axL.semilogy(Lv, [r['E_det_exactBC_rmse'] for r in boundary_rows],
                 'o-', color='tab:green', ms=5, lw=1.6,
                 label='exact transformed boundary data')
    axL.semilogy(Lv, [r['exact_phi_recovery_rmse'] for r in boundary_rows],
                 'd--', color='tab:gray', ms=5, lw=1.2,
                 label='exact transformed field')
    axL.set_xlabel('domain size $L$'); axL.set_ylabel('deterministic RMSE')
    axL.legend(fontsize=8); axL.grid(True, which='both', alpha=0.3)
    axL.text(0.03, 0.03, '(a)', transform=axL.transAxes, fontsize=11)
    Lv = [r['L'] for r in domain_rows]
    axR.errorbar(Lv, [r['E_total_mean'] for r in domain_rows],
                 yerr=[r['E_total_std'] for r in domain_rows],
                 fmt='^-', color='k', ms=5, lw=1.6, capsize=3,
                 label=r'$E_{\mathrm{total}}$ (mean $\pm$ std)')
    axR.errorbar(Lv, [r['E_grw_mean'] for r in domain_rows],
                 yerr=[r['E_grw_std'] for r in domain_rows],
                 fmt='o-', color='tab:blue', ms=5, lw=1.6, capsize=3,
                 label=r'$E_{\mathrm{GRW}}$ (mean $\pm$ std)')
    axR.plot(Lv, [r['E_det_rmse'] for r in domain_rows],
             's--', color='tab:red', ms=5, lw=1.4,
             label=r'$E_{\mathrm{det}}$ (deterministic)')
    axR.set_xlabel('domain size $L$'); axR.set_ylabel('RMSE')
    axR.legend(fontsize=8); axR.grid(True, alpha=0.3)
    axR.text(0.03, 0.03, '(b)', transform=axR.transAxes, fontsize=11)
    fig.tight_layout()
    _savefig(fig, _mk(OUT_BASE, 'burgers_boundary_domain'))

    # perturbation response curves (white vs kernel-smoothed)
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    amps = [r['amplitude'] for r in pert['rows']]
    ax.errorbar(amps, [r['white_l2_mean'] for r in pert['rows']],
                yerr=[r['white_l2_std'] for r in pert['rows']],
                fmt='o-', color='tab:blue', ms=5, lw=1.6, capsize=3,
                label='white perturbation')
    ax.errorbar(amps, [r['smooth_l2_mean'] for r in pert['rows']],
                yerr=[r['smooth_l2_std'] for r in pert['rows']],
                fmt='s-', color='tab:green', ms=5, lw=1.6, capsize=3,
                label='kernel-smoothed perturbation')
    ax.plot([pert['measured_phi_rmse']], [pert['measured_u_l2']],
            marker='*', color='tab:red', ms=14, ls='none',
            label='GRW computation')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'perturbation amplitude on $\phi$ (RMS)')
    ax.set_ylabel(r'$L^2$ error of recovered $u$')
    ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    _savefig(fig, _mk(OUT_BASE, 'burgers_perturbation_response'))


# -------------------------------------------------------------------- main --
def run_task8():
    os.makedirs(OUT_BASE, exist_ok=True)
    t0 = time.perf_counter()
    print(f"{'='*60}\n  Task 8: Burgers controlled attribution\n{'='*60}")

    print('\n  Part V: pipeline validation against the packaged solver')
    max_dev = validate_pipeline()

    decoupled_rows, partA_ci, decoupled_per_run = run_decoupled_controls()

    print('\n  Part D: deterministic boundary controls')
    boundary_rows = run_boundary_controls()

    # measured particle pipeline point at the paper configuration (part A,
    # P=400): x = phi reconstruction RMSE, y = particle-induced u error in
    # the L_h^2 norm (RMSE times sqrt(h M))
    ref_row = [r for r in decoupled_rows
               if r['part'] == 'A' and r['P'] == 400][0]
    measured_phi = ref_row['rmse_phi_mean']
    conv = float(np.sqrt((4.0 / 399) * 400))
    measured_u_l2 = ref_row['rmse_particle_mean'] * conv
    print(f"\n  Part E: perturbation response (measured phi RMSE = "
          f"{measured_phi:.2e}, particle u error = {measured_u_l2:.3f} L2)")
    pert = run_perturbation_response(measured_phi, measured_u_l2)

    # verification: exact-BC FD solve against the analytic transformed
    # solution, and the decoupled particle-refinement slope
    kk = A_AMP / (2.0 * NU)
    exact_fd = _fd_reference(400, 4.0, boundary='exact')
    x_v = exact_fd['x']
    phi_analytic = (np.exp(NU * kk**2 * T_END)
                    * np.cosh(kk * (x_v - 2.0)) / np.cosh(kk * 2.0))
    phi_bc_dev = float(np.max(np.abs(exact_fd['phi'] - phi_analytic)))
    P_v = np.array([r['P'] for r in decoupled_rows
                    if r['part'] == 'A'], dtype=float)
    e_v = np.array([r['rmse_particle_mean'] for r in decoupled_rows
                    if r['part'] == 'A'])
    slopeA = float(np.polyfit(np.log10(P_v), np.log10(e_v), 1)[0])
    print(f"  [check] exact-BC FD vs analytic phi: max dev = {phi_bc_dev:.2e}")
    print(f"  [check] part-A particle-error slope vs P: {slopeA:.3f}")

    print('\n  Part F: multi-seed domain study')
    domain_rows, domain_trend, domain_per_run = run_domain_multiseed()

    with open(_mk(OUT_BASE, 'decoupled.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(decoupled_rows[0].keys()))
        w.writeheader()
        for r in decoupled_rows:
            w.writerow(r)
    with open(_mk(OUT_BASE, 'domain_multiseed.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(domain_rows[0].keys()))
        w.writeheader()
        for r in domain_rows:
            w.writerow(r)
    with open(_mk(OUT_BASE, 'decoupled_per_run.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(decoupled_per_run[0].keys()))
        w.writeheader()
        for r in decoupled_per_run:
            w.writerow(r)
    with open(_mk(OUT_BASE, 'domain_per_run.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(domain_per_run[0].keys()))
        w.writeheader()
        for r in domain_per_run:
            w.writerow(r)

    summary = {
        'params': {'A': A_AMP, 'nu': NU, 'dt': DT, 'T': T_END,
                   'S_decoupled': S_DECOUPLED, 'S_domain': S_DOMAIN,
                   'seeds': f'{BASE_SEED}..'},
        'validation_max_dev': max_dev,
        'verification': {'exact_bc_fd_vs_analytic_max_dev': phi_bc_dev,
                         'partA_particle_error_slope_vs_P': slopeA,
                         'partA_slope_ci_realization': partA_ci},
        'decoupled': decoupled_rows,
        'boundary_controls': boundary_rows,
        'perturbation_response': pert,
        'domain_multiseed': domain_rows,
        'domain_total_trend': domain_trend,
    }
    with open(_mk(OUT_BASE, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    make_figures(decoupled_rows, boundary_rows, pert, domain_rows)
    print(f"\n  [Task 8] Done in {time.perf_counter()-t0:.1f}s. "
          f"Outputs in {OUT_BASE}/")
    return summary


if __name__ == '__main__':
    run_task8()
