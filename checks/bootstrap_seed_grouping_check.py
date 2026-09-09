"""Sensitivity of the fitted-trend confidence intervals to how seeds are resampled.

Every ensemble study reuses the same seed list at each parameter value. The
reported realization-level bootstrap resamples the S realizations independently
within each parameter value. This check recomputes those intervals and compares
them with intervals obtained by resampling the seed labels jointly across all
parameter values, which preserves any dependence created by the shared seeds.

Coverage
  from pinned per-realization norms (no solver runs):
    heat E_total(N)                          root mean square of per-realization L_h^2 errors
    reaction profile, location, speed, aligned   means of per-realization errors
    Burgers initialization-point refinement  mean particle RMSE versus P (log-log slope)
    Burgers domain study                     mean total RMSE versus L (linear slope)
  with --heat-profiles (reruns the paired heat ensemble, about two minutes):
    heat total, spread, and bias slopes for the coupled, fixed-300 bin-center,
    and fixed-300 right-edge treatments, from realization profiles. The rerun
    is first verified against pinned_ensembles/heat_grid_paired/summary_by_N_paired.csv.

Run from the repository root:
    python checks/bootstrap_seed_grouping_check.py [--heat-profiles] [--json PATH]
"""
import csv, json, os, sys
from collections import defaultdict
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'studies'))
B, SEED = 5000, 303
REPORTED = {
    'heat E_total': (-0.493, -0.423), 'reaction profile': (-0.496, -0.424),
    'reaction front location': (-0.607, -0.441), 'reaction speed': (-0.607, -0.441),
    'reaction aligned profile': (-0.475, -0.414),
    'Burgers particle RMSE vs P': (-0.524, -0.466),
    'Burgers total RMSE vs L (linear)': (-2.88e-3, 5.44e-4),
    'heat coupled total': (-0.493, -0.423), 'heat coupled spread': (-0.507, -0.430),
}
RESULTS = []

def load(path, key, cols, seedcol='seed'):
    d = {c: defaultdict(dict) for c in cols}
    for r in csv.DictReader(open(path)):
        for c in cols:
            d[c][float(r[key])][int(r[seedcol])] = float(r[c])
    return d

def matrix(byK):
    Ks = sorted(byK); seeds = sorted(byK[Ks[0]])
    return np.array(Ks, float), np.array([[byK[K][s] for s in seeds] for K in Ks])

def fit(Ks, stat, loglog=True):
    return float(np.polyfit(np.log(Ks), np.log(stat), 1)[0]) if loglog else float(np.polyfit(Ks, stat, 1)[0])

def interval(Ks, A, agg, joint, loglog=True):
    rng = np.random.default_rng(SEED); S = A.shape[1]; out = []
    for _ in range(B):
        if joint:
            idx = rng.integers(0, S, S); stat = np.array([agg(A[k, idx]) for k in range(len(Ks))])
        else:
            stat = np.array([agg(A[k, rng.integers(0, S, S)]) for k in range(len(Ks))])
        out.append(fit(Ks, stat, loglog))
    lo, hi = np.percentile(out, [2.5, 97.5]); return float(lo), float(hi)

def report(name, Ks, A, agg, loglog=True, tol=0.005):
    est = fit(Ks, np.array([agg(A[k]) for k in range(len(Ks))]), loglog)
    ind = interval(Ks, A, agg, False, loglog); jnt = interval(Ks, A, agg, True, loglog)
    rep = REPORTED.get(name)
    ok = rep is None or (abs(ind[0]-rep[0]) <= tol*max(1, abs(rep[0])) and abs(ind[1]-rep[1]) <= tol*max(1, abs(rep[1])))
    RESULTS.append(dict(statistic=name, fit=est, independent=ind, joint=jnt, reported=rep, reproduced=ok))
    print(f"{name:34s} fit {est:+.4g}  independent [{ind[0]:+.4g},{ind[1]:+.4g}]  joint [{jnt[0]:+.4g},{jnt[1]:+.4g}]"
          + (f"  reported [{rep[0]:+.4g},{rep[1]:+.4g}] {'reproduced' if ok else 'DIFFERS'}" if rep else ""))
    assert ok, f"{name}: independent interval does not reproduce the reported one"

def heat_profiles():
    """Rerun the paired heat ensemble with study t7's runner and compare total/spread/bias slopes."""
    from study_t7_heat_grid_paired import _run_heat_one, _decompose
    from verify_solver import exact_heat_step
    from utils import reconstruct_cumulative
    alpha, T, dt, L, x0, uL, uR, S, base = 0.5, 0.5, 0.005, 4.0, 2.0, 0.0, 1.0, 30, 42
    Ns = [500, 1000, 2000, 5000, 10000, 20000, 50000]
    pinned = {}
    for r in csv.DictReader(open(os.path.join(ROOT, 'pinned_ensembles/heat_grid_paired/summary_by_N_paired.csv'))):
        pinned[(r['treatment'], int(r['N']))] = r
    M = 300; edges = np.linspace(0, L, M+1); centers = 0.5*(edges[:-1]+edges[1:]); h = float(edges[1]-edges[0])
    ref_c = exact_heat_step(centers, T, x0, uL, uR, alpha); ref_e = exact_heat_step(edges[1:], T, x0, uL, uR, alpha)
    grams = {'coupled': [], 'fixed300': [], 'fixed300e': []}
    for N in Ns:
        fixed_runs, coupled_runs = [], []
        xg = np.linspace(0, L, N); ref_cpl = exact_heat_step(xg, T, x0, uL, uR, alpha); h_cpl = float(xg[1]-xg[0])
        for s in range(S):
            x_pos, w_arr, _, _, _ = _run_heat_one(N, alpha, T, dt, L, x0, uL, uR, base+s)
            _, u = reconstruct_cumulative(x_pos, w_arr, edges, uL); fixed_runs.append(u)
            order = np.argsort(x_pos); xs = x_pos[order]; us = uL + np.cumsum(w_arr[order])
            coupled_runs.append(np.interp(xg, xs, us, left=uL, right=uR))
        u_f = np.array(fixed_runs); u_c = np.array(coupled_runs)
        for tr, ref, arr, hh in (('fixed300', ref_c, u_f, h), ('fixed300e', ref_e, u_f, h), ('coupled', ref_cpl, u_c, h_cpl)):
            got = _decompose(arr, ref, hh)[:3]; row = pinned[(tr, N)]
            exp = [float(row['E_bias']), float(row['E_spread']), float(row['E_total'])]
            assert all(abs(g-e) <= 5.01e-9 for g, e in zip(got, exp)), f'{tr} N={N} differs from pinned'
            d = arr - ref[None, :]; grams[tr].append((d @ d.T) * hh)
        print(f"  heat profiles N={N} reproduce pinned statistics")
    lN = np.log(np.array(Ns, float))
    def stats_from(G, c):
        Sn = len(c); Et2 = float(c @ np.diag(G))/Sn; Eb2 = float(c @ G @ c)/Sn**2; return Et2, max(Et2-Eb2, 0.0), Eb2
    for tr in ('coupled', 'fixed300', 'fixed300e'):
        for key, k in (('total', 0), ('spread', 1), ('bias', 2)):
            def run(joint):
                rng = np.random.default_rng(SEED); out = []
                for _ in range(B):
                    vals = []
                    idx_joint = rng.integers(0, S, S)
                    for G in grams[tr]:
                        idx = idx_joint if joint else rng.integers(0, S, S)
                        c = np.bincount(idx, minlength=S).astype(float); vals.append(stats_from(G, c)[k])
                    v = np.sqrt(np.maximum(np.array(vals), 1e-30)); out.append(float(np.polyfit(lN, np.log(v), 1)[0]))
                lo, hi = np.percentile(out, [2.5, 97.5]); return float(lo), float(hi)
            est = float(np.polyfit(lN, np.log(np.sqrt([stats_from(G, np.ones(S))[k] for G in grams[tr]])), 1)[0])
            ind, jnt = run(False), run(True); name = f'heat {tr} {key}'; rep = REPORTED.get(name)
            ok = rep is None or (abs(ind[0]-rep[0]) <= 0.005 and abs(ind[1]-rep[1]) <= 0.005)
            RESULTS.append(dict(statistic=name, fit=est, independent=ind, joint=jnt, reported=rep, reproduced=ok))
            print(f"{name:34s} fit {est:+.4g}  independent [{ind[0]:+.4g},{ind[1]:+.4g}]  joint [{jnt[0]:+.4g},{jnt[1]:+.4g}]"
                  + (f"  reported [{rep[0]:+.4g},{rep[1]:+.4g}] {'reproduced' if ok else 'DIFFERS'}" if rep else ""))
            assert ok, f'{name}: independent interval does not reproduce the reported one'

if __name__ == '__main__':
    rms = lambda v: float(np.sqrt(np.mean(v**2))); mean = lambda v: float(np.mean(v))
    h = load(os.path.join(ROOT, 'pinned_ensembles/heat_extended/per_run.csv'), 'N', ['l2'])
    Ks, A = matrix(h['l2']); report('heat E_total', Ks, A, rms)
    f = load(os.path.join(ROOT, 'pinned_ensembles/fhn_extended/per_run.csv'), 'N', ['l2', 'center_err', 'speed_err', 'aligned_err'])
    for col, name in [('l2', 'reaction profile'), ('center_err', 'reaction front location'),
                      ('speed_err', 'reaction speed'), ('aligned_err', 'reaction aligned profile')]:
        Ks, A = matrix(f[col]); report(name, Ks, A, mean)
    bp = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(ROOT, 'pinned_ensembles/burgers_controls/decoupled_per_run.csv'))):
        if r['part'] == 'A': bp[float(r['P'])][int(r['seed'])] = float(r['rmse_particle'])
    Ks, A = matrix(bp); report('Burgers particle RMSE vs P', Ks, A, mean)
    bd = load(os.path.join(ROOT, 'pinned_ensembles/burgers_controls/domain_per_run.csv'), 'L', ['E_total_rmse'])
    Ks, A = matrix(bd['E_total_rmse']); report('Burgers total RMSE vs L (linear)', Ks, A, mean, loglog=False, tol=0.02)
    if '--heat-profiles' in sys.argv: heat_profiles()
    if '--json' in sys.argv:
        json.dump(RESULTS, open(sys.argv[sys.argv.index('--json')+1], 'w'), indent=1)
