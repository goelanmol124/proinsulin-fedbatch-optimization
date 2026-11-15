import time
import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C, WhiteKernel
from wetlab_sim import _rhs, P as TRUE_P, V0, Xv0, Xd0, S0_init, A0, P_sol0, P_ib0, DO0, BOUNDS, run_batch, simulate
warnings.filterwarnings('ignore')


_oracle_cache = {}
def oracle_score(x):
    key = tuple(np.round(np.asarray(x, dtype=float), 6))
    v = _oracle_cache.get(key)
    if v is None:
        tr = simulate(x, seed=0, scatter=False)
        v = float(tr['P_sol_final']) if tr['ok'] else 0.0
        _oracle_cache[key] = v
    return v

def _normalize(x):
    return (np.asarray(x) - BOUNDS[:, 0]) / (BOUNDS[:, 1] - BOUNDS[:, 0])

def _denormalize(z):
    return BOUNDS[:, 0] + np.asarray(z) * (BOUNDS[:, 1] - BOUNDS[:, 0])

def _biased_params():
    p = dict(TRUE_P)
    p['T_opt_p'] = 32.0
    p['alpha'] = 0.06
    p['k_tox_IPTG'] = 0.05
    p['kLa_c'] = 0.027 * 1.25
    p['Y_AS'] = 0.2
    return p
_BIASED = _biased_params()

def low_fidelity(x, t_max=None):
    T_prod, IPTG, X_ind, F_s, t_f = x
    u = (T_prod, IPTG, X_ind, F_s, 500.0, 1000.0, 2.0)
    y0 = [V0, Xv0, Xd0, S0_init, A0, P_sol0, P_ib0, DO0, _BIASED['T_opt_g']]
    if t_max is None:
        t_max = float(t_f)
    sol = solve_ivp(_rhs, (0.0, t_max), y0, args=(u, _BIASED), method='LSODA', t_eval=[t_max], rtol=0.0001, atol=1e-06, max_step=1.0)
    if not sol.success or np.any(np.isnan(sol.y)):
        return 0.0
    return float(max(0.0, sol.y[5, -1]))

class MFSurrogate:

    def __init__(self):
        kernel = C(0.3, (0.001, 5.0)) * Matern(length_scale=[0.3] * 5, length_scale_bounds=(0.05, 5.0), nu=2.5) + WhiteKernel(noise_level=0.05 ** 2, noise_level_bounds=(1e-05, 0.5))
        self.gp = GaussianProcessRegressor(kernel=kernel, normalize_y=False, n_restarts_optimizer=2, alpha=1e-06)
        self.X = np.zeros((0, 5))
        self.y = np.zeros(0)
        self.lf_cache = {}

    def _lf(self, x):
        key = tuple(np.round(x, 6))
        v = self.lf_cache.get(key)
        if v is None:
            v = low_fidelity(x)
            self.lf_cache[key] = v
        return v

    def lf_only(self, X):
        X = np.atleast_2d(X)
        return np.array([self._lf(x) for x in X])

    def fit(self, X, y):
        self.X = np.asarray(X, dtype=float)
        self.y = np.asarray(y, dtype=float)
        lf = self.lf_only(self.X)
        residual = self.y - lf
        Z = _normalize(self.X)
        if residual.std() < 1e-09:
            residual = residual + 1e-06 * np.arange(len(residual))
        self.gp.fit(Z, residual)

    def predict(self, X):
        X = np.atleast_2d(X)
        lf = self.lf_only(X)
        Z = _normalize(X)
        if len(self.X) == 0:
            return (lf, np.full(len(X), self.gp.kernel.k1.k1.constant_value ** 0.5))
        mu_d, std_d = self.gp.predict(Z, return_std=True)
        return (np.maximum(0.0, lf + mu_d), std_d)

def _lhs(n, seed):
    sampler = qmc.LatinHypercube(d=BOUNDS.shape[0], seed=int(seed))
    Z = sampler.random(n)
    return (_denormalize(Z), Z)

def method_trial_error(seed, budget=100, n_scout=25, scatter=True, sigma_rel=0.08):
    rng = np.random.default_rng(int(seed))
    X_scout, _ = _lhs(n_scout, seed)
    history = []
    best = (-np.inf, None, -np.inf)
    for i, x in enumerate(X_scout):
        Pt, Po, _ = run_batch(x, seed=int(seed) * 1000 + i, scatter=scatter)
        history.append((x, Po, Pt))
        if Po > best[0]:
            best = (Po, x, Pt)
    span = BOUNDS[:, 1] - BOUNDS[:, 0]
    n_refine = budget - n_scout
    for k in range(n_refine):
        x_next = np.clip(best[1] + rng.normal(0, sigma_rel * span), BOUNDS[:, 0], BOUNDS[:, 1])
        Pt, Po, _ = run_batch(x_next, seed=int(seed) * 1000 + n_scout + k, scatter=scatter)
        history.append((x_next, Po, Pt))
        if Po > best[0]:
            best = (Po, x_next, Pt)
    return dict(method='trial_error', best_x=best[1], best_P_obs=best[0],
                best_P_true=oracle_score(best[1]), history=history, n_wetlab=len(history))

def _box_behnken_5():
    pairs = [(i, j) for i in range(5) for j in range(i + 1, 5)]
    rows = []
    for i, j in pairs:
        for si in (-1, 1):
            for sj in (-1, 1):
                row = [0] * 5
                row[i] = si
                row[j] = sj
                rows.append(row)
    rows += [[0] * 5] * 6
    return np.array(rows, dtype=float)

def _quad_features(z):
    z = np.atleast_2d(z)
    n, d = z.shape
    feats = [np.ones((n, 1)), z, z * z]
    for i in range(d):
        for j in range(i + 1, d):
            feats.append((z[:, i] * z[:, j]).reshape(-1, 1))
    return np.hstack(feats)

def _bb_Z01():
    return (_box_behnken_5() + 1.0) / 2.0

def _rsm_shrink_box(x_center, shrink):
    lo_g, hi_g = (BOUNDS[:, 0], BOUNDS[:, 1])
    span_g = hi_g - lo_g
    c = np.asarray(x_center, dtype=float).ravel()
    L = shrink * span_g
    lo = np.maximum(lo_g, c - 0.5 * L)
    hi = np.minimum(hi_g, c + 0.5 * L)
    min_span = 0.15 * span_g
    deficit = np.maximum(0.0, min_span - (hi - lo))
    lo = lo - 0.5 * deficit
    hi = hi + 0.5 * deficit
    lo = np.maximum(lo, lo_g)
    hi = np.minimum(hi, hi_g)
    lo = np.minimum(lo, hi - 1e-09 * span_g)
    return np.column_stack([lo, hi])

def _box_behnken_5_compact(n_target):
    base = _box_behnken_5()
    base = base[:40]
    if n_target >= 40:
        return base
    return base[:max(21, n_target)]


def method_doe_rsm(seed, budget=100, scatter=True, shrink=0.45):
    history = []
    sid = int(seed)

    def one_stage(box, Z01_design, run_base, de_add):
        lo, hi = box[:, 0], box[:, 1]
        X = lo + Z01_design * (hi - lo)
        Y = np.zeros(len(X))
        for i, x in enumerate(X):
            Pt, Po, _ = run_batch(x, seed=sid * 1000 + run_base + i, scatter=scatter)
            history.append((x, Po, Pt))
            Y[i] = Po
        F = _quad_features(Z01_design)
        coefs, *_ = np.linalg.lstsq(F, Y, rcond=None)

        def neg_surf(zv):
            return -float(_quad_features(np.asarray(zv, dtype=float).reshape(1, -1)) @ coefs)
        res = differential_evolution(neg_surf, [(0.0, 1.0)] * 5, seed=sid + de_add,
                                     popsize=10, maxiter=30, tol=1e-4)
        x_de = lo + res.x * (hi - lo)
        return x_de, X, Y

    if budget >= 92:
        Z01 = _bb_Z01()
        _, X1, Y1 = one_stage(BOUNDS, Z01, 0, 0)
        x0 = X1[int(np.argmax(Y1))]
        box_loc = _rsm_shrink_box(x0, shrink)
        x_de, X2, Y2 = one_stage(box_loc, Z01, 6000, 31)
        x_obs = X2[int(np.argmax(Y2))] if Y2.max() >= Y1.max() else X1[int(np.argmax(Y1))]
        Pt_de, Po_de, _ = run_batch(x_de, seed=sid * 1000 + 9800, scatter=scatter)
        history.append((x_de, Po_de, Pt_de))
        Pt_o, Po_o, _ = run_batch(x_obs, seed=sid * 1000 + 9810, scatter=scatter)
        history.append((x_obs, Po_o, Pt_o))
        x_rec = x_de if Po_de >= Po_o else x_obs
    else:
        n_design = min(40, budget - 1)
        Z01_c = (_box_behnken_5_compact(n_design) + 1.0) / 2.0
        x_de, X1, Y1 = one_stage(BOUNDS, Z01_c, 0, 0)
        x_obs = X1[int(np.argmax(Y1))]
        if len(history) < budget - 1:
            Pt_de, Po_de, _ = run_batch(x_de, seed=sid * 1000 + 9800, scatter=scatter)
            history.append((x_de, Po_de, Pt_de))
            x_rec = x_de if Po_de >= np.max(Y1) else x_obs
        else:
            x_rec = x_obs

    return dict(method='doe_rsm', best_x=x_rec, best_P_obs=float(np.max([h[1] for h in history])),
                best_P_true=oracle_score(x_rec), history=history, n_wetlab=len(history))

def method_de_sim(seed, scatter=True):

    def f(x):
        Pt, _, _ = run_batch(x, seed=0, scatter=False, noisy=False)
        return -Pt
    bounds = [tuple(b) for b in BOUNDS]
    res = differential_evolution(f, bounds, seed=int(seed), popsize=10, maxiter=12, tol=0.001, polish=True)
    x_star = res.x
    Pt, Po, _ = run_batch(x_star, seed=int(seed) * 1000, scatter=scatter)
    history = [(x_star, Po, Pt)]
    return dict(method='de_sim', best_x=x_star, best_P_obs=Po, best_P_true=Pt, history=history, n_wetlab=1, n_surrogate=res.nfev)

def _fit_gp(Z, y):
    kernel = C(1.0, (0.01, 100.0)) * Matern(length_scale=[0.3] * Z.shape[1], length_scale_bounds=(0.01, 10.0), nu=2.5) + WhiteKernel(noise_level=0.05, noise_level_bounds=(0.0001, 1.0))
    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=2, alpha=1e-06)
    gp.fit(Z, y)
    return gp

def _ei(gp, z_cands, y_best, xi=0.01):
    mu, std = gp.predict(z_cands, return_std=True)
    std = np.maximum(std, 1e-09)
    from scipy.stats import norm
    imp = mu - y_best - xi
    Z = imp / std
    return imp * norm.cdf(Z) + std * norm.pdf(Z)

def method_gp_bo(seed, budget=100, n_init=5, scatter=True):
    n_iter = budget - n_init - 1
    rng = np.random.default_rng(int(seed))
    X_init, Z_init = _lhs(n_init, seed)
    history = []
    Zs, ys, y_trues = (list(Z_init), [], [])
    for i, x in enumerate(X_init):
        Pt, Po, _ = run_batch(x, seed=int(seed) * 1000 + i, scatter=scatter)
        ys.append(Po)
        y_trues.append(Pt)
        history.append((x, Po, Pt))
    Zs = list(np.asarray(Zs))
    for k in range(n_iter):
        Z_arr = np.asarray(Zs)
        y_arr = np.asarray(ys)
        gp = _fit_gp(Z_arr, y_arr)
        z_best = Z_arr[int(np.argmax(y_arr))]
        cand_u = rng.random((1500, BOUNDS.shape[0]))
        cand_l = np.clip(z_best + rng.normal(0, 0.08, size=(800, BOUNDS.shape[0])), 0, 1)
        cand = np.vstack([cand_u, cand_l])
        ei = _ei(gp, cand, y_best=float(np.max(y_arr)))
        z0 = cand[np.argmax(ei)]
        res = minimize(lambda z: -float(_ei(gp, z.reshape(1, -1), y_best=float(np.max(y_arr)))[0]), z0, bounds=[(0.0, 1.0)] * BOUNDS.shape[0], method='L-BFGS-B')
        z_next = np.clip(res.x, 0.0, 1.0)
        x_next = _denormalize(z_next)
        Pt, Po, _ = run_batch(x_next, seed=int(seed) * 1000 + n_init + k, scatter=scatter)
        Zs.append(z_next)
        ys.append(Po)
        y_trues.append(Pt)
        history.append((x_next, Po, Pt))
    Z_arr = np.asarray(Zs)
    y_arr = np.asarray(ys)
    gp = _fit_gp(Z_arr, y_arr)
    cand = rng.random((5000, BOUNDS.shape[0]))
    mu = gp.predict(cand)
    z_rec = cand[np.argmax(mu)]
    x_rec = _denormalize(z_rec)
    Pt_rec, Po_rec, _ = run_batch(x_rec, seed=int(seed) * 1000 + 9999, scatter=scatter)
    history.append((x_rec, Po_rec, Pt_rec))
    return dict(method='gp_bo', best_x=x_rec, best_P_obs=Po_rec,
                best_P_true=oracle_score(x_rec), history=history, n_wetlab=len(history))

def method_mfbo_ude(seed, budget=100, n_init=5, scatter=True):
    n_iter = budget - (n_init + 1) - 1
    rng = np.random.default_rng(int(seed))
    history = []

    def lf_neg(x):
        return -low_fidelity(x)
    res_warm = differential_evolution(lf_neg, [tuple(b) for b in BOUNDS], seed=int(seed), popsize=8, maxiter=10, polish=False)
    x_warm = res_warm.x
    X_init, _ = _lhs(n_init, seed)
    X_train = np.vstack([X_init, x_warm[None, :]])
    y_train = []
    for i, x in enumerate(X_train):
        Pt, Po, _ = run_batch(x, seed=int(seed) * 1000 + i, scatter=scatter)
        history.append((x, Po, Pt))
        y_train.append(Po)
    y_train = np.array(y_train)
    surr = MFSurrogate()
    surr.fit(X_train, y_train)
    cand_z = rng.random((600, 5))
    cand_x = _denormalize(cand_z)
    cand_lf = surr.lf_only(cand_x)
    for k in range(n_iter):
        beta = max(0.05, 1.0 - 0.02 * k)
        y_best = float(np.max(y_train))
        if k % 5 == 0 and len(y_train) > 0:
            x_best_now = X_train[int(np.argmax(y_train))]
            jitter = rng.normal(0, 0.05, size=(80, 5)) * (BOUNDS[:, 1] - BOUNDS[:, 0])
            x_local = np.clip(x_best_now + jitter, BOUNDS[:, 0], BOUNDS[:, 1])
            cand_z = np.vstack([cand_z, _normalize(x_local)])
            cand_lf = np.concatenate([cand_lf, surr.lf_only(x_local)])
        if len(surr.X) == 0:
            mu = cand_lf
            std = np.full(len(cand_lf), np.sqrt(surr.gp.kernel.k1.k1.constant_value))
        else:
            mu_d, std_d = surr.gp.predict(cand_z, return_std=True)
            mu = np.maximum(0.0, cand_lf + mu_d)
            std = std_d
        from scipy.stats import norm
        std_safe = np.maximum(std, 1e-09)
        imp = mu - y_best - 0.01
        Zsc = imp / std_safe
        ei = imp * norm.cdf(Zsc) + std_safe * norm.pdf(Zsc)
        score = ei + beta * std
        z_next = cand_z[np.argmax(score)]
        x_next = _denormalize(z_next.reshape(1, -1))[0]
        Pt, Po, _ = run_batch(x_next, seed=int(seed) * 1000 + n_init + 1 + k, scatter=scatter)
        history.append((x_next, Po, Pt))
        cand_lf_new = surr.lf_only(np.atleast_2d(x_next))
        cand_z = np.vstack([cand_z, _normalize(x_next)])
        cand_lf = np.concatenate([cand_lf, cand_lf_new])
        X_train = np.vstack([X_train, x_next[None, :]])
        y_train = np.append(y_train, Po)
        surr.fit(X_train, y_train)
    if len(surr.X) > 0:
        mu_d, _ = surr.gp.predict(cand_z, return_std=True)
        mu_final = np.maximum(0.0, cand_lf + mu_d)
    else:
        mu_final = cand_lf
    x_rec = _denormalize(cand_z[np.argmax(mu_final)].reshape(1, -1))[0]
    Pt_rec, Po_rec, _ = run_batch(x_rec, seed=int(seed) * 1000 + 9999, scatter=scatter)
    history.append((x_rec, Po_rec, Pt_rec))
    return dict(method='mfbo_ude', best_x=x_rec, best_P_obs=Po_rec,
                best_P_true=oracle_score(x_rec), history=history, n_wetlab=len(history))
METHODS = {'trial_error': method_trial_error, 'doe_rsm': method_doe_rsm, 'gp_bo': method_gp_bo, 'mfbo_ude': method_mfbo_ude}

def run_all(seeds=range(5), budget=40):
    rows = []
    histories = {m: [] for m in METHODS}
    for s in seeds:
        for name, fn in METHODS.items():
            t0 = time.time()
            r = fn(seed=int(s), budget=budget)
            wall = time.time() - t0
            rows.append((name, int(s), r['best_P_true'], r['best_P_obs'], r['n_wetlab'], wall, *r['best_x']))
            histories[name].append(r['history'])
            print(f"seed={s:2d} {name:12s} P_true={r['best_P_true']:.3f}  P_obs={r['best_P_obs']:.3f}  n_wetlab={r['n_wetlab']:3d}  wall={wall:5.1f}s  x={np.round(r['best_x'], 3)}")
    return (rows, histories)
if __name__ == '__main__':
    import sys
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    budget = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    t0 = time.time()
    rows, histories = run_all(range(n_seeds), budget=budget)
    print(f'\ntotal wall-clock: {time.time() - t0:.1f}s')
    dtype = [('method', 'U16'), ('seed', 'i4'), ('P_true', 'f8'), ('P_obs', 'f8'), ('n_wetlab', 'i4'), ('wall_s', 'f8'), ('T_prod', 'f8'), ('IPTG', 'f8'), ('X_ind', 'f8'), ('F_s', 'f8'), ('t_f', 'f8')]
    arr = np.array(rows, dtype=dtype)
    np.savez('results.npz', table=arr, **{f'hist_{m}': np.array(h, dtype=object) for m, h in histories.items()})
    print('wrote results.npz')
    print('\n=== Summary (median across seeds) ===')
    print(f"{'method':15s}{'med P_true':>12s}{'IQR':>10s}{'med n_wetlab':>14s}{'med wall(s)':>14s}")
    for m in METHODS:
        sub = arr[arr['method'] == m]
        q1, med, q3 = np.percentile(sub['P_true'], [25, 50, 75])
        print(f"{m:15s}{med:12.3f}{q3 - q1:10.3f}{int(np.median(sub['n_wetlab'])):14d}{np.median(sub['wall_s']):14.1f}")
