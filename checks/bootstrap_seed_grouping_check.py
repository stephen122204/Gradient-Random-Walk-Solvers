"""Sensitivity of the fitted-exponent confidence intervals to how seeds are resampled.

The studies reuse the same seed list at every particle count N. The reported
realization-level bootstrap resamples the S realizations independently within
each N. This check recomputes those intervals from the pinned per-realization
norms and compares them with intervals obtained by resampling the seed labels
jointly across all N, which preserves any dependence between counts.

Statistics covered (those formed from per-realization norms):
  heat E_total(N)  = root mean square of the per-realization L_h^2 errors
  reaction profile, front-location, speed, and aligned-profile errors = means

Run from the repository root:
    python checks/bootstrap_seed_grouping_check.py
"""
import csv, os
from collections import defaultdict
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B, SEED = 5000, 303
REPORTED = {  # independent-within-count intervals stated in the manuscript
    'heat E_total': (-0.493, -0.423),
    'reaction profile': (-0.496, -0.424),
    'reaction front location': (-0.607, -0.441),
    'reaction speed': (-0.607, -0.441),
    'reaction aligned profile': (-0.475, -0.414),
}

def load(path, cols):
    d = {c: defaultdict(dict) for c in cols}
    for r in csv.DictReader(open(path)):
        for c in cols:
            d[c][int(r['N'])][int(r['seed'])] = float(r[c])
    return d

def matrix(byN):
    Ns = sorted(byN); seeds = sorted(byN[Ns[0]])
    return np.array(Ns, float), np.array([[byN[N][s] for s in seeds] for N in Ns])

def slope(Ns, stat):
    return float(np.polyfit(np.log(Ns), np.log(stat), 1)[0])

def interval(Ns, A, agg, joint):
    rng = np.random.default_rng(SEED); S = A.shape[1]; out = []
    for _ in range(B):
        if joint:
            idx = rng.integers(0, S, S)
            stat = np.array([agg(A[k, idx]) for k in range(len(Ns))])
        else:
            stat = np.array([agg(A[k, rng.integers(0, S, S)]) for k in range(len(Ns))])
        out.append(slope(Ns, stat))
    lo, hi = np.percentile(out, [2.5, 97.5])
    return float(lo), float(hi)

def report(name, Ns, A, agg):
    fit = slope(Ns, np.array([agg(A[k]) for k in range(len(Ns))]))
    ind = interval(Ns, A, agg, joint=False)
    jnt = interval(Ns, A, agg, joint=True)
    rep = REPORTED[name]
    ok = abs(ind[0] - rep[0]) <= 0.005 and abs(ind[1] - rep[1]) <= 0.005
    contains = lambda iv: iv[0] <= -0.5 <= iv[1]
    print(f"{name:26s} fit {fit:.3f}  independent [{ind[0]:.3f},{ind[1]:.3f}] (reported [{rep[0]:.3f},{rep[1]:.3f}], "
          f"{'reproduced' if ok else 'DIFFERS'})  joint-by-seed [{jnt[0]:.3f},{jnt[1]:.3f}]  "
          f"contains -1/2: independent {contains(ind)}, joint {contains(jnt)}")
    assert ok, f"{name}: independent-scheme interval does not reproduce the reported one"

if __name__ == '__main__':
    rms = lambda v: float(np.sqrt(np.mean(v ** 2)))
    mean = lambda v: float(np.mean(v))
    h = load(os.path.join(ROOT, 'pinned_ensembles/heat_extended/per_run.csv'), ['l2'])
    Ns, A = matrix(h['l2']); report('heat E_total', Ns, A, rms)
    f = load(os.path.join(ROOT, 'pinned_ensembles/fhn_extended/per_run.csv'),
             ['l2', 'center_err', 'speed_err', 'aligned_err'])
    for col, name in [('l2', 'reaction profile'), ('center_err', 'reaction front location'),
                      ('speed_err', 'reaction speed'), ('aligned_err', 'reaction aligned profile')]:
        Ns, A = matrix(f[col]); report(name, Ns, A, mean)
