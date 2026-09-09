"""Check the paired heat error identities and their finite-ensemble uncertainty scales."""
import argparse, csv, sys
from math import erf, sqrt
from pathlib import Path
import numpy as np

ALPHA, T, X0, L, S, DT = 0.5, 0.5, 2.0, 4.0, 30, 0.005
SIG = sqrt(4 * ALPHA * T)
Phi = np.vectorize(lambda z: 0.5 * (1 + erf(z)))

def F_reflected(x):
    tot = np.zeros_like(x, dtype=float)
    for n in range(-4, 5):
        for mu in (2 * n * L + X0, 2 * n * L - X0):
            tot += Phi((x - mu) / SIG) - Phi((0 - mu) / SIG)
    return tot

def ref_inf(x):
    return Phi((x - X0) / SIG)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pinned", default=str(Path(__file__).resolve().parents[1] / "pinned_ensembles/heat_grid_paired/summary_by_N_paired.csv"))
    ap.add_argument("--fresh", action="store_true", help="also run the heat walk for four independent seed blocks")
    a = ap.parse_args()
    print(f"NumPy {np.__version__}; exact second moments; Gaussian/delta SD approximations")
    print("RMS = sqrt(E[statistic**2]), not E[statistic]; approx_sd is not exact")
    rows = list(csv.DictReader(open(a.pinned)))
    for M in (300, 400):
        h = L / M; b = (np.arange(M) + 1) * h; c = (np.arange(M) + 0.5) * h
        Fb = F_reflected(b)
        Fi, Fj = np.meshgrid(Fb, Fb, indexing="ij"); Sig1 = np.minimum(Fi, Fj) - Fi * Fj
        d = ref_inf(b) - ref_inf(c); dn = sqrt(h * np.sum(d ** 2))
        print(f"M={M}: control center={sqrt(h*np.sum((Fb-ref_inf(c))**2)):.5f} edge={sqrt(h*np.sum((Fb-ref_inf(b))**2)):.5f}  ||d||={dn:.5f}")
        for conv, ref, tag in (("center", ref_inf(c), f"fixed{M}"), ("edge", ref_inf(b), f"fixed{M}e")):
            D = Fb - ref; pred2 = h * np.sum(D ** 2)
            for r in rows:
                if r["treatment"] != tag: continue
                N = int(r["N"]); A = h * Sig1 / N
                Etot = sqrt(pred2 + np.trace(A))
                Esp2 = (S - 1) / S * np.trace(A); Vsp2 = (S - 1) / S ** 2 * 2 * np.sum(A * A)
                C = A / S; Eb2 = pred2 + np.trace(C); Vb2 = 4 * h * (D @ C @ D) + 2 * np.sum(C * C)
                print(f"  {conv:6s} N={N:6d} total meas={float(r['E_total']):.5f} RMS={Etot:.5f} | spread meas={float(r['E_spread']):.5f} RMS={sqrt(Esp2):.5f} approx_relsd={0.5*sqrt(Vsp2)/Esp2:.3f} | bias meas={float(r['E_bias']):.5f} RMS={sqrt(Eb2):.5f} approx_sd={sqrt(Vb2)/(2*sqrt(Eb2)):.5f}")
        bc = [float(x["E_bias"]) for x in rows if x["treatment"] == f"fixed{M}" and x["N"] == "50000"][0]
        be = [float(x["E_bias"]) for x in rows if x["treatment"] == f"fixed{M}e" and x["N"] == "50000"][0]
        print(f"  paired identity at N=50000: sqrt(bc^2-be^2)={sqrt(bc**2-be**2):.5f}  ||d||={dn:.5f}  cross term={bc**2-be**2-dn**2:.3e}")
    if a.fresh:
        nsteps = int(round(T / DT)); sig = sqrt(2 * ALPHA * DT); M = 300; h = L / M
        def walk(N, seed):
            rng = np.random.RandomState(seed); x = np.full(N, X0)
            for _ in range(nsteps):
                x = x + rng.normal(0.0, sig, size=N)
                while np.any((x < 0) | (x > L)):
                    x = np.where(x < 0, -x, x); x = np.where(x > L, 2 * L - x, x)
            return x
        for N in (5000, 50000):
            for blk in (42, 1000, 2000, 3000):
                U = np.array([np.cumsum(np.histogram(walk(N, s), bins=np.linspace(0, L, M + 1))[0]) / N for s in range(blk, blk + S)])
                ub = U.mean(0); spread = sqrt(np.mean([h * np.sum((u - ub) ** 2) for u in U]))
                print(f"fresh N={N} seeds {blk}..{blk+S-1}: spread={spread:.5f}")

if __name__ == "__main__":
    sys.exit(main())
