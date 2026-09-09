"""Two-step heat profile: predicted error floors, transition particle counts,
and a particle-count rule, tested against the solver (manuscript Section 5.4).

Initial data u(x,0) = u_L + A H(x-a) + B H(x-b) on [0,L] with 0<a<b<L and
A != B. N_a = N_b = N/2 globs are placed at a and b with weights A/N_a and
B/N_b, evolved by the heat solver with weight-preserving reflection, and
reconstructed by bin-and-sum on M bins.

References
  finite-interval  m(x,T) = u_L + A F_R(x,T;a) + B F_R(x,T;b), F_R the reflected CDF
  infinite-line    u_inf(x,T) = u_L + A Phi(x;a) + B Phi(x;b), Phi the error-function CDF
Comparison conventions: reference at bin right edges (aligned) or at bin centers.

Predictions (independent globs, deterministic weights, fixed allocation fractions)
  E[E_total^2]  = B_h^2 + V_h/N
  E[E_bias^2]   = B_h^2 + V_h/(N S)
  E[E_spread^2] = (S-1)/S * V_h/N
  B_h^2 = h sum_j (m(b_j) - g_j)^2,   V_h = h sum_j [A^2/f_a F_a(1-F_a) + B^2/f_b F_b(1-F_b)](b_j)
  N_* = V_h / B_h^2 (B_h > 0),   N >= V_h / (eps^2 - B_h^2) for a target RMS error eps > B_h.

Usage (from the repository root):
  python studies/study_t9_heat_two_step.py predict           # deterministic predictions only
  python studies/study_t9_heat_two_step.py pilot             # time one N=50000 run
  python studies/study_t9_heat_two_step.py run               # production ensembles (seeds 5000-5029)
  python studies/study_t9_heat_two_step.py validate          # held-out target demonstration (seeds 7000-7029)
Outputs go to output/heat_two_step/.
Use --output-dir PATH to write a separate reproduction. The predict command
also generates validation_spec.json before either ensemble is run.

The cumulative construction and finite-interval reflection are classical
(Ghoniem and Sherman, 1985, Sections II-III, DOI 10.1016/0021-9991(85)90058-0).
The weighted Bernoulli variance is established mathematics; see also BPC
(2024), Appendix A.2, DOI 10.1007/s10915-024-02614-1. This study applies it
to the reported comparison points and tests the resulting predictions.
"""
import argparse, json, os, sys, time
from math import ceil, erf, sqrt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from config import SimulationConfig            # noqa: E402
from simulation import simulate_heat_equation  # noqa: E402
from utils import reconstruct_cumulative       # noqa: E402

OUT = os.path.join(ROOT, 'output', 'heat_two_step')

# ---- problem and design -----------------------------------------------------
ALPHA, L, T, DT = 0.5, 4.0, 0.5, 0.005
UL, A, B, XA, XB = 0.0, 0.3, 0.7, 1.4, 2.6
FA = FB = 0.5                       # allocation fractions N_a = N_b = N/2
M_LIST = [50, 100, 200, 400]
N_LIST = [250, 500, 1000, 2000, 4000, 8000, 16000, 32000, 64000, 128000]
S = 30
SEEDS_RUN = list(range(5000, 5030))
SEEDS_VAL = list(range(7000, 7030))
N_IMAGES = 4

Phi = np.vectorize(lambda z: 0.5 * (1.0 + erf(z)))
SIG = sqrt(4.0 * ALPHA * T)         # erf argument scale: (x-mu)/sqrt(4 alpha T)

def F_inf(x, x0):
    return Phi((np.asarray(x, float) - x0) / SIG)

def G_images(x, x0, n_images=N_IMAGES):
    x = np.asarray(x, float); tot = np.zeros_like(x)
    for n in range(-n_images, n_images + 1):
        for mu in (2 * n * L + x0, 2 * n * L - x0):
            tot += Phi((x - mu) / SIG)
    return tot

def F_R(x, x0, n_images=N_IMAGES):
    return G_images(x, x0, n_images) - G_images(0.0, x0, n_images)

def m_finite(x):  return UL + A * F_R(x, XA) + B * F_R(x, XB)
def u_inf(x):     return UL + A * F_inf(x, XA) + B * F_inf(x, XB)

def grid(M):
    edges = np.linspace(0.0, L, M + 1)
    return edges, edges[1:], 0.5 * (edges[:-1] + edges[1:]), float(edges[1] - edges[0])

def predictions(M):
    edges, be, ce, h = grid(M)
    Fa, Fb = F_R(be, XA), F_R(be, XB)
    Vh = h * float(np.sum(A**2 / FA * Fa * (1 - Fa) + B**2 / FB * Fb * (1 - Fb)))
    m_e = m_finite(be)
    out = {'M': M, 'h': h, 'V_h': Vh}
    for ref_name, ref in (('finite', m_finite), ('infinite', u_inf)):
        for conv, pts in (('edge', be), ('center', ce)):
            B2 = h * float(np.sum((m_e - ref(pts))**2))
            key = f'{ref_name}_{conv}'
            out[key] = {'B2': B2, 'B': sqrt(B2), 'N_star': (Vh / B2 if B2 > 0 else None)}
    return out

def covariance(M, N):
    """Exact covariance of the binned reconstruction at the bin right edges.

    Independent groups add, and within a group the indicators satisfy
    I(x)I(y) = I(min(x,y)), so for group weight w and count n the covariance is
    w^2 n [F(min) - F(x)F(y)]. The factor h converts to the L_h^2 inner product.
    """
    edges, be, ce, h = grid(M)
    Na = Nb = N // 2
    out = np.zeros((len(be), len(be)))
    for x0, wt, n in ((XA, A / Na, Na), (XB, B / Nb, Nb)):
        F = F_R(be, x0)
        Fmin = F_R(np.minimum.outer(be, be).ravel(), x0).reshape(len(be), len(be))
        out += wt**2 * n * (Fmin - np.outer(F, F))
    return h * out


def image_truncation_check():
    x = np.linspace(0, L, 4001)
    d = max(np.max(np.abs(F_R(x, XA, 4) - F_R(x, XA, 8))), np.max(np.abs(F_R(x, XB, 4) - F_R(x, XB, 8))))
    return {'max_abs_diff_4_vs_8_images': float(d),
            'F_R(L;a)-1': float(F_R(L, XA) - 1), 'F_R(L;b)-1': float(F_R(L, XB) - 1)}


def validation_spec():
    """Recompute the paper's fixed target design without running particles.

    Round the bound upward to an even count to preserve the half/half allocation.
    This reproduces the historical design; it is not a new prospective experiment.
    """
    M, target = 200, 0.008
    p = predictions(M)
    edge, center = p['infinite_edge'], p['infinite_center']
    if target <= edge['B']:
        raise ValueError('The edge floor must be below the fixed validation target')
    N = 2 * ceil(p['V_h'] / (target**2 - edge['B2']) / 2)
    return {
        'purpose': 'held-out accuracy-target demonstration, specified before any validation run',
        'reference': 'infinite-line (the reference a user would normally have)',
        'M': M, 'target_rms_error': target,
        'center_convention': {'B': center['B'], 'attainable': target > center['B']},
        'edge_convention': {'B': edge['B'], 'attainable': True, 'N_required': N,
                            'predicted_E_total_at_N': sqrt(edge['B2'] + p['V_h']/N)},
        'predicted_E_total_center_at_same_N': sqrt(center['B2'] + p['V_h']/N),
        'N': N, 'seeds': f'{SEEDS_VAL[0]}-{SEEDS_VAL[-1]}',
        'criterion': 'measured RMS realization error (E_total) against the infinite-line reference at the bin edges compared with eps and with the predicted value; center convention reported at the same N',
    }

# ---- solver path --------------------------------------------------------------
def run_one(N, seed):
    np.random.seed(seed)
    Na = Nb = N // 2
    ic = [(XA, A / Na)] * Na + [(XB, B / Nb)] * Nb
    cfg = SimulationConfig(equation_type='heat', domain_type='Finite', domain_size=L,
                           boundary_conditions={'LEFT': {'type': 'Dirichlet', 'value': UL},
                                                'RIGHT': {'type': 'Dirichlet', 'value': UL + A + B}},
                           diff_constant=ALPHA, time_step=DT, total_time=T, num_points=N,
                           initial_conditions=ic, reaction_term=False)
    globs = [{'position': float(p), 'value': float(w)} for p, w in ic]
    t0 = time.perf_counter()
    res = simulate_heat_equation(globs, cfg)
    el = time.perf_counter() - t0
    return (np.array([g['position'] for g in res]), np.array([g['value'] for g in res]), el)

def decompose(u_arr, ref, h):
    ubar = u_arr.mean(axis=0)
    Eb2 = h * float(np.sum((ubar - ref)**2))
    Et2 = h * float(np.mean(np.sum((u_arr - ref)**2, axis=1)))
    Es2 = h * float(np.mean(np.sum((u_arr - ubar)**2, axis=1)))
    return Eb2, Es2, Et2

def ensemble(N, seeds):
    recon = {M: [] for M in M_LIST}; runtimes = []
    for s in seeds:
        x, w, el = run_one(N, s); runtimes.append(el)
        for M in M_LIST:
            edges = grid(M)[0]
            _, u = reconstruct_cumulative(x, w, edges, UL)
            recon[M].append(u)
    rows = []
    for M in M_LIST:
        edges, be, ce, h = grid(M)
        u_arr = np.array(recon[M])
        for ref_name, ref in (('finite', m_finite), ('infinite', u_inf)):
            for conv, pts in (('edge', be), ('center', ce)):
                g = ref(pts)
                Eb2, Es2, Et2 = decompose(u_arr, g, h)
                per_real = (h * np.sum((u_arr - g[None, :])**2, axis=1)).tolist()
                rows.append({'N': N, 'M': M, 'reference': ref_name, 'convention': conv,
                             'E_bias': sqrt(Eb2), 'E_spread': sqrt(Es2), 'E_total': sqrt(Et2),
                             'S': len(seeds), 'per_realization_sq': per_real})
    return rows, runtimes

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', nargs='?', default='predict', choices=['predict', 'pilot', 'run', 'validate'])
    parser.add_argument('--output-dir', default=OUT)
    args = parser.parse_args()
    mode, OUT = args.mode, os.path.abspath(args.output_dir)
    os.makedirs(OUT, exist_ok=True)
    if mode == 'predict':
        pred = {'design': dict(alpha=ALPHA, L=L, T=T, dt=DT, u_L=UL, A=A, B=B, a=XA, b=XB, f_a=FA, f_b=FB,
                               M_list=M_LIST, N_list=N_LIST, S=S, seeds_run=[SEEDS_RUN[0], SEEDS_RUN[-1]],
                               seeds_validation=[SEEDS_VAL[0], SEEDS_VAL[-1]]),
                'image_truncation': image_truncation_check(),
                'per_M': [predictions(M) for M in M_LIST]}
        json.dump(pred, open(os.path.join(OUT, 'predictions.json'), 'w'), indent=1)
        json.dump(validation_spec(), open(os.path.join(OUT, 'validation_spec.json'), 'w'), indent=1)
        for p in pred['per_M']:
            print(f"M={p['M']:4d} h={p['h']:.4f} V_h={p['V_h']:.4e} | "
                  f"B(fin,center)={p['finite_center']['B']:.3e} N*={p['finite_center']['N_star']:.0f} | "
                  f"B(inf,edge)={p['infinite_edge']['B']:.3e} N*={p['infinite_edge']['N_star']:.0f} | "
                  f"B(inf,center)={p['infinite_center']['B']:.3e} N*={p['infinite_center']['N_star']:.0f} | "
                  f"B(fin,edge)={p['finite_edge']['B']:.1e}")
        print("image truncation:", pred['image_truncation'])
    elif mode == 'pilot':
        x, w, el = run_one(50000, 5000)
        print(f"pilot N=50000: {el:.2f} s, sum w = {w.sum():.6f}, positions in [0,L]: {x.min():.3f}..{x.max():.3f}")
    elif mode == 'run':
        all_rows = []; rt = {}
        for N in N_LIST:
            t0 = time.perf_counter(); rows, runtimes = ensemble(N, SEEDS_RUN); all_rows += rows
            rt[N] = {'mean_run_s': float(np.mean(runtimes)), 'wall_s': time.perf_counter() - t0}
            print(f"N={N:6d} done in {rt[N]['wall_s']:.1f}s")
        json.dump({'rows': all_rows, 'runtime': rt}, open(os.path.join(OUT, 'ensembles.json'), 'w'), indent=1)
    elif mode == 'validate':
        spec = json.load(open(os.path.join(OUT, 'validation_spec.json')))
        rows, runtimes = ensemble(spec['N'], SEEDS_VAL)
        json.dump({'spec': spec, 'rows': rows, 'mean_run_s': float(np.mean(runtimes))},
                  open(os.path.join(OUT, 'validation.json'), 'w'), indent=1)
        for r in rows:
            if r['M'] == spec['M']: print(r)
