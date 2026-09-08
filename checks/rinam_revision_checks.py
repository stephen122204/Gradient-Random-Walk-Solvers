"""Focused checks for the RINAM-revision code changes.

(a) FHN default path is bit-identical before/after the reaction-derivative
    interface, and an explicit callback equal to the built-in polynomial gives
    identical output.
(b) The coordinate-aware reconstruct_cumulative path reproduces the pinned
    paired-heat statistics (fixed300 / fixed300e) for the requested N.
Run from the repository root:  python checks/rinam_revision_checks.py [--pristine DIR]
"""
import csv, os, subprocess, sys, tempfile, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'studies'))

FHN_SNIPPET = r'''
import sys, numpy as np
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[1] + "/studies")
from study_t5_fhn_extended import _run_fhn_one
import study_t5_fhn_extended as t5
from simulation import simulate_fitzhugh_nagumo_grw
from config import SimulationConfig, generate_fhn_steady_ic
N, a, nu, T, dt, L, xc, seed = 500, 0.25, 0.5, 5.0, 0.01, 30.0, 15.0, 42
np.random.seed(seed)
ic = generate_fhn_steady_ic(N, a, x_center=xc)
cfg = SimulationConfig(equation_type='fitzhugh-nagumo', domain_type='Finite', domain_size=L,
    boundary_conditions={'LEFT': {'type': 'Neumann', 'value': 0.0}, 'RIGHT': {'type': 'Neumann', 'value': 0.0}},
    diff_constant=nu, time_step=dt, total_time=T, num_points=N, initial_conditions=ic, reaction_term=True, a=a)
if len(sys.argv) > 3 and sys.argv[3] == 'callback':
    theta = np.sqrt(2.0) * (0.5 - a)
    c2, c1, c0 = -1.5 * nu, 1.5 * nu - theta, 0.5 * theta - 0.25 * nu
    cfg.reaction_derivative = lambda u: c2 * u**2 + c1 * u + c0
globs = [{'position': float(p), 'value': float(w)} for p, w in ic]
res = simulate_fitzhugh_nagumo_grw(globs, cfg)
np.savez(sys.argv[2], x=np.array([g['position'] for g in res]), w=np.array([g['value'] for g in res]))
'''

def run_fhn(repo, out, mode='default'):
    with tempfile.NamedTemporaryFile('w', suffix='.py', delete=False) as f:
        f.write(FHN_SNIPPET); path = f.name
    subprocess.run([sys.executable, path, repo, out, mode], check=True)
    return np.load(out)

def check_fhn(pristine):
    t0 = time.perf_counter()
    d = tempfile.mkdtemp()
    new = run_fhn(ROOT, os.path.join(d, 'new.npz'))
    cb = run_fhn(ROOT, os.path.join(d, 'cb.npz'), 'callback')
    ident_cb = np.array_equal(new['x'], cb['x']) and np.array_equal(new['w'], cb['w'])
    print(f"(a) callback vs default identical: {ident_cb}  "
          f"max|dx|={np.max(np.abs(new['x']-cb['x'])):.3e} max|dw|={np.max(np.abs(new['w']-cb['w'])):.3e}")
    if pristine:
        old = run_fhn(pristine, os.path.join(d, 'old.npz'))
        ident = np.array_equal(new['x'], old['x']) and np.array_equal(new['w'], old['w'])
        print(f"(a) pristine vs new identical: {ident}  "
              f"max|dx|={np.max(np.abs(new['x']-old['x'])):.3e} max|dw|={np.max(np.abs(new['w']-old['w'])):.3e}")
    print(f"(a) elapsed {time.perf_counter()-t0:.1f}s")

def check_paired_heat(N_list=(5000, 50000)):
    from study_t7_heat_grid_paired import _run_heat_one, _decompose
    from verify_solver import exact_heat_step
    from utils import reconstruct_cumulative
    alpha, T, dt, L, x0, uL, uR, S, base = 0.5, 0.5, 0.005, 4.0, 2.0, 0.0, 1.0, 30, 42
    pinned = {}
    with open(os.path.join(ROOT, 'pinned_ensembles/heat_grid_paired/summary_by_N_paired.csv')) as f:
        for r in csv.DictReader(f):
            pinned[(r['treatment'], int(r['N']))] = r
    M = 300
    edges = np.linspace(0.0, L, M + 1); centers = 0.5 * (edges[:-1] + edges[1:]); dxM = float(edges[1] - edges[0])
    ref_c = exact_heat_step(centers, T, x0, uL, uR, alpha)
    ref_e = exact_heat_step(edges[1:], T, x0, uL, uR, alpha)
    worst = 0.0
    for N in N_list:
        t0 = time.perf_counter(); runs = []
        for s in range(S):
            x_pos, w_arr, _, _, _ = _run_heat_one(N, alpha, T, dt, L, x0, uL, uR, base + s)
            x_eval, u = reconstruct_cumulative(x_pos, w_arr, edges, uL)
            assert np.array_equal(x_eval, edges[1:])
            runs.append(u)
        u_arr = np.array(runs)
        for tr, ref in (('fixed300', ref_c), ('fixed300e', ref_e)):
            got = _decompose(u_arr, ref, dxM)[:3]
            row = pinned[(tr, N)]
            exp = [float(row['E_bias']), float(row['E_spread']), float(row['E_total'])]
            rel = max(abs(g - e) / e for g, e in zip(got, exp))
            worst = max(worst, rel)
            print(f"(b) N={N} {tr:10s} got={['%.8f'%g for g in got]} pinned={['%.8f'%e for e in exp]} max_rel={rel:.2e}")
        print(f"(b) N={N} elapsed {time.perf_counter()-t0:.1f}s")
    print(f"(b) worst relative deviation from pinned CSV: {worst:.2e} (CSV has 8 decimals)")

if __name__ == '__main__':
    pristine = None
    if '--pristine' in sys.argv:
        pristine = sys.argv[sys.argv.index('--pristine') + 1]
    check_fhn(pristine)
    check_paired_heat()
