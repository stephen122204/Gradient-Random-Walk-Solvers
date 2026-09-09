"""Two controls for the two-step heat prediction.

(a) Reference control. Samples the exact final-time reflected distribution directly,
    with no time stepping, and compares the resulting reconstruction mean with the
    closed-form finite-interval solution m. This isolates the reference formula from
    the solver.
(b) Seed-block control. Repeats the solver ensemble on independent blocks of thirty
    seeds at one particle count and records the deviation of the measured statistics
    from the prediction, together with the alignment term whose sign is at issue.
    The production ensembles reuse one seed list at every particle count, so their
    deviations are correlated across counts and cannot show whether a systematic
    effect is present.

Run from the repository root:
    python checks/two_step_reference_and_seed_blocks.py [--blocks K] [--N N]
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'studies'))
from study_t9_heat_two_step import (ALPHA, L, T, UL, A, B, XA, XB, grid, F_R,   # noqa: E402
                                    m_finite, covariance, run_one)
from utils import reconstruct_cumulative                                         # noqa: E402

def fold(y):
    """Reflect a real position into [0, L] repeatedly."""
    y = np.abs(y) % (2 * L)
    return np.where(y > L, 2 * L - y, y)

def reference_control(M, n=20_000_000, seed=11):
    """Direct sampling of the folded Gaussian at time T, no time stepping."""
    edges, be, ce, h = grid(M)
    rng = np.random.default_rng(seed)
    sig = np.sqrt(2 * ALPHA * T)
    Na = Nb = n // 2
    xa = fold(XA + sig * rng.standard_normal(Na))
    xb = fold(XB + sig * rng.standard_normal(Nb))
    x = np.concatenate([xa, xb]); w = np.concatenate([np.full(Na, A / Na), np.full(Nb, B / Nb)])
    _, u = reconstruct_cumulative(x, w, edges, UL)
    m = m_finite(be)
    err = np.abs(u - m)
    sd_point = np.sqrt(A**2 / Na * 0.25 + B**2 / Nb * 0.25)
    print(f"(a) M={M} direct sampling of the final-time law, n={n:,}")
    print(f"    max |mean reconstruction - m| = {err.max():.3e}, pointwise sampling sd = {sd_point:.3e}, "
          f"ratio {err.max()/sd_point:.2f}")
    return err.max() / sd_point

def seed_blocks(N, M_list, blocks, S=30, start=90000):
    print(f"\n(b) independent seed blocks at N={N}, S={S}, {blocks} blocks")
    print(f"    {'block':>7} " + " ".join(f"M={M}: E_tot/pred(center)  align" for M in M_list))
    out = {M: [] for M in M_list}
    for k in range(blocks):
        seeds = range(start + 1000 * k, start + 1000 * k + S)
        recon = {M: [] for M in M_list}
        for s in seeds:
            x, w, _ = run_one(N, s)
            for M in M_list:
                _, u = reconstruct_cumulative(x, w, grid(M)[0], UL)
                recon[M].append(u)
        line = f"    {k:>7} "
        for M in M_list:
            edges, be, ce, h = grid(M)
            u_arr = np.array(recon[M]); m_e = m_finite(be); g_c = m_finite(ce)
            Et2 = h * float(np.mean(np.sum((u_arr - g_c[None, :])**2, axis=1)))
            V = float(np.trace(covariance(M, N)))
            B2 = h * float(np.sum((m_e - g_c)**2))
            ratio = np.sqrt(Et2) / np.sqrt(B2 + V)
            align = h * float(np.sum((u_arr.mean(axis=0) - m_e) * (m_e - g_c)))
            out[M].append((ratio, align))
            line += f"  {ratio:.4f}  {align:+.3e}"
        print(line)
    print()
    for M in M_list:
        r = np.array([v[0] for v in out[M]]); a = np.array([v[1] for v in out[M]])
        print(f"    M={M}: ratio mean {r.mean():.4f} sd {r.std(ddof=1):.4f} | "
              f"alignment term {a.mean():+.3e} sd {a.std(ddof=1):.3e}, negative in {int((a<0).sum())}/{len(a)} blocks")
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--blocks', type=int, default=6)
    ap.add_argument('--N', type=int, default=16000)
    ap.add_argument('--n-direct', type=int, default=20_000_000)
    a = ap.parse_args()
    for M in (50, 200):
        reference_control(M, n=a.n_direct)
    seed_blocks(a.N, [50, 200], a.blocks)
