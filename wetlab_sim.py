import numpy as np
from scipy.integrate import solve_ivp
P = dict(mu_max=0.6, Ks=0.05, Ki_S=80.0, Kdo=0.2, Ki_A=8.0, kd_0=0.02, k_tox_A=0.005, k_tox_IPTG=0.1, Yxs=0.45, ms=0.03, Yxo=1.1, mo=0.4, Ypo=0.5, alpha=0.12, beta=0.0, K_IPTG=0.05, n_Hill=2, T_opt_g=37.0, sigma_g=5.0, T_opt_p=28.0, sigma_p=4.0, qs_crit=0.4, Y_AS=0.4, k_au=0.1, K_au=0.5, Yps=0.15, T_IB_50=32.0, sigma_IB=2.0, mu_crit_IB=0.3, phi_min=0.4, kLa_c=0.027, kLa_aN=1.5, kLa_bF=0.5, DO_star=7.5, rho=1.0, Cp=4.18, K_jacket=200.0, Q_ag_ref=10.0, H_O2=0.46, T_feed=25.0)
V0, Xv0, S0_init, A0, P_sol0, P_ib0, DO0, Xd0 = (2.5, 0.1, 10.0, 0.0, 0.0, 0.0, 7.5, 0.0)

def _sig(x, c=0.0, w=1.0):
    return 1.0 / (1.0 + np.exp(-(x - c) / w))

def _rhs(t, y, u, p):
    V, Xv, Xd, S, A, P_sol, P_ib, DO, T = y
    T_set_prod, IPTG_dose, X_ind, F_s, S_f, N_rpm, F_air_vvm = u
    ind = _sig(Xv - X_ind, 0.0, 0.5)
    IPTG = IPTG_dose * ind
    hill = IPTG ** p['n_Hill'] / (p['K_IPTG'] ** p['n_Hill'] + IPTG ** p['n_Hill'] + 1e-12)
    F_active = F_s * _sig(1.0 - S, 0.0, 0.3)
    D = F_active / V
    T_set = p['T_opt_g'] + (T_set_prod - p['T_opt_g']) * ind
    kLa = p['kLa_c'] * N_rpm ** p['kLa_aN'] * F_air_vvm ** p['kLa_bF']
    fT_g = np.exp(-0.5 * ((T - p['T_opt_g']) / p['sigma_g']) ** 2)
    mu = p['mu_max'] * fT_g * S / (p['Ks'] + S + S * S / p['Ki_S']) * DO / (p['Kdo'] + DO) * p['Ki_A'] / (p['Ki_A'] + A)
    phi = 1.0 - (1.0 - p['phi_min']) * ind * hill
    fT_p = np.exp(-0.5 * ((T - p['T_opt_p']) / p['sigma_p']) ** 2)
    qp = (p['alpha'] * mu + p['beta']) * hill * (p['Ki_A'] / (p['Ki_A'] + A)) * fT_p
    qs = mu / p['Yxs'] + p['ms'] + qp / max(p['Yps'], 1e-06)
    q_ac_over = p['Y_AS'] * max(0.0, qs - p['qs_crit'])
    q_ac_up = p['k_au'] * A / (p['K_au'] + A) * (1.0 - _sig(S, 1.0, 0.3))
    f_sol = 1.0 / (1.0 + np.exp((T - p['T_IB_50']) / p['sigma_IB'])) * 1.0 / (1.0 + (mu / p['mu_crit_IB']) ** 2)
    kd = p['kd_0'] + p['k_tox_A'] * A + p['k_tox_IPTG'] * IPTG
    q_O2 = mu / p['Yxo'] + p['mo'] + qp / max(p['Ypo'], 1e-06)
    OUR_vol = q_O2 * Xv * 1000.0
    Q_met = p['H_O2'] * (OUR_vol / 32.0) * V
    Q_cool = max(0.0, p['K_jacket'] * (T - T_set))
    Q_ag = p['Q_ag_ref'] * (N_rpm / 600.0) ** 3
    H_feed = F_active * p['rho'] * p['Cp'] * (p['T_feed'] - T)
    dV = F_active
    dXv = (mu * phi - kd) * Xv - D * Xv
    dXd = kd * Xv - D * Xd
    dS = D * (S_f - S) - Xv * qs
    dA = (q_ac_over - q_ac_up) * Xv - D * A
    dP_sol = qp * f_sol * Xv - D * P_sol
    dP_ib = qp * (1.0 - f_sol) * Xv - D * P_ib
    dDO = kLa * (p['DO_star'] - DO) - OUR_vol - D * DO
    dT = (Q_met + Q_ag - Q_cool + H_feed) / (V * p['rho'] * p['Cp'])
    return [dV, dXv, dXd, dS, dA, dP_sol, dP_ib, dDO, dT]

def _draw_scatter(rng):
    return dict(mu_max=float(rng.lognormal(0.0, 0.05)), Yxs=float(rng.lognormal(0.0, 0.04)), alpha=float(rng.lognormal(0.0, 0.08)), kLa_c=float(rng.lognormal(0.0, 0.1)), Yx_OD=float(rng.lognormal(0.0, 0.08)))

def simulate(x, seed=0, scatter=True, t_max=None):
    T_prod, IPTG, X_ind, F_s, t_f = x
    S_f, N_rpm, F_air_vvm = (500.0, 1000.0, 2.0)
    rng = np.random.default_rng(int(seed))
    p = dict(P)
    if scatter:
        s = _draw_scatter(rng)
        p['mu_max'] *= s['mu_max']
        p['Yxs'] *= s['Yxs']
        p['alpha'] *= s['alpha']
        p['kLa_c'] *= s['kLa_c']
        Yx_OD = s['Yx_OD']
    else:
        Yx_OD = 1.0
    u = (T_prod, IPTG, X_ind, F_s, S_f, N_rpm, F_air_vvm)
    y0 = [V0, Xv0, Xd0, S0_init, A0, P_sol0, P_ib0, DO0, p['T_opt_g']]
    if t_max is None:
        t_max = float(t_f)
    t_eval = np.linspace(0.0, t_max, 200)
    sol = solve_ivp(_rhs, (0.0, t_max), y0, args=(u, p), method='BDF', t_eval=t_eval, rtol=1e-06, atol=1e-08, max_step=0.5)
    if not sol.success or np.any(np.isnan(sol.y)):
        return dict(t=t_eval, y=np.zeros((9, len(t_eval))), P_sol_final=0.0, P_ib_final=0.0, ok=False, Yx_OD=Yx_OD, scatter=scatter)
    return dict(t=sol.t, y=sol.y, P_sol_final=float(sol.y[5, -1]), P_ib_final=float(sol.y[6, -1]), Xv_final=float(sol.y[1, -1]), A_final=float(sol.y[4, -1]), DO_min=float(np.min(sol.y[7])), T_max=float(np.max(sol.y[8])), ok=True, Yx_OD=Yx_OD, scatter=scatter)

def measure(traj, rng):
    if not traj['ok']:
        return 0.0
    P_total_true = traj['P_sol_final'] + traj['P_ib_final']
    f_sol_true = traj['P_sol_final'] / max(P_total_true, 1e-09)
    P_total_obs = P_total_true * max(0.0, 1.0 + rng.normal(0.0, 0.07))
    f_sol_obs = float(np.clip(f_sol_true * (1.0 + rng.normal(0.0, 0.2)), 0.0, 1.0))
    return float(P_total_obs * f_sol_obs)
BOUNDS = np.array([[25.0, 34.0], [0.01, 0.5], [5.0, 30.0], [0.02, 0.2], [16.0, 36.0]])

def run_batch(x, seed=0, scatter=True, noisy=True, csv_path=None, run_id=None, batch_id=1, instant_feedback=False):
    traj = simulate(x, seed=seed, scatter=scatter)
    rng = np.random.default_rng(int(seed) + 10000)
    P_sol_reported = measure(traj, rng) if noisy else traj['P_sol_final']
    if csv_path is not None:
        _write_batch_csv(csv_path, x, traj, rng, run_id or 'run', batch_id, instant_feedback)
    return (traj['P_sol_final'], P_sol_reported, traj)
ONLINE_CHANNELS = ['T_broth', 'pH', 'DO', 'RPM', 'F_air', 'weight', 'OUR', 'CER', 'RQ']
OFFLINE_ASSAYS = [('OD600', 0.5, 0.1, 2 / 60, 10 / 60, 'offline_OD'), ('DCW', 2.0, 0.04, 5 / 60, 26.0, 'offline_DCW'), ('glucose', 1.0, 0.03, 5 / 60, 50 / 60, 'offline_HPLC'), ('acetate', 1.0, 0.06, 5 / 60, 50 / 60, 'offline_HPLC'), ('proinsulin_total', 2.0, 0.07, 5 / 60, 75 / 60, 'offline_HPLC'), ('soluble_fraction', 4.0, 0.2, 10 / 60, 5.0, 'offline_PAGE')]

def _write_batch_csv(path, x, traj, rng, run_id, batch_id, instant_feedback):
    import csv
    if not traj['ok']:
        return
    t = traj['t']
    y = traj['y']
    T_prod, IPTG_dose, X_ind, F_s, t_f = x
    controls = f'{{"T_set_prod":{T_prod:.2f},"IPTG":{IPTG_dose:.3f},"X_ind":{X_ind:.2f},"F_s":{F_s:.3f},"t_f":{t_f:.1f}}}'
    rows = []
    t_online = np.arange(0.0, float(t_f), 10 / 60)
    V = np.interp(t_online, t, y[0])
    Xv = np.interp(t_online, t, y[1])
    S = np.interp(t_online, t, y[3])
    A = np.interp(t_online, t, y[4])
    DO = np.interp(t_online, t, y[7])
    T_broth = np.interp(t_online, t, y[8])
    for i, tt in enumerate(t_online):
        mu_est = max(0.0, 0.6 * S[i] / (0.05 + S[i]) * DO[i] / (0.2 + DO[i]))
        OUR_est = (mu_est / 1.1 + 0.4) * Xv[i] * 1000.0
        CER_est = OUR_est * (1.0 + 0.1 * (A[i] > 1.0))
        RQ_est = CER_est / max(OUR_est, 1.0)
        vals = dict(T_broth=T_broth[i] + rng.normal(0, 0.1), pH=7.0 + rng.normal(0, 0.02), DO=DO[i] / 7.5 * 100 + rng.normal(0, 1.5), RPM=1000.0, F_air=2.0 + rng.normal(0, 0.02), weight=V[i] + rng.normal(0, 0.005), OUR=OUR_est * (1 + rng.normal(0, 0.05)), CER=CER_est * (1 + rng.normal(0, 0.05)), RQ=RQ_est)
        units = dict(T_broth='C', pH='-', DO='%sat', RPM='rpm', F_air='vvm', weight='kg', OUR='mgO2/L/h', CER='mgCO2/L/h', RQ='-')
        for ch in ONLINE_CHANNELS:
            rows.append((run_id, batch_id, tt, tt, ch, vals[ch], units[ch], 'online', False, 0.0, controls))
    for name, dt_assay, sigma_y, sigma_t, turnaround, source in OFFLINE_ASSAYS:
        t_grid = np.arange(dt_assay, float(t_f) + 1e-06, dt_assay)
        for tt in t_grid:
            t_reported = tt + rng.normal(0, sigma_t)
            if name == 'OD600':
                Xv_t = np.interp(tt, t, y[1])
                v_true = 0.5 * Xv_t / max(traj.get('Yx_OD', 1.0), 0.5)
                unit = 'AU'
            elif name == 'DCW':
                v_true = np.interp(tt, t, y[1])
                unit = 'g/L'
            elif name == 'glucose':
                v_true = max(0.0, np.interp(tt, t, y[3]))
                unit = 'g/L'
            elif name == 'acetate':
                v_true = max(0.0, np.interp(tt, t, y[4]))
                unit = 'g/L'
            elif name == 'proinsulin_total':
                v_true = max(0.0, np.interp(tt, t, y[5]) + np.interp(tt, t, y[6]))
                unit = 'g/L'
            elif name == 'soluble_fraction':
                p_sol_t = np.interp(tt, t, y[5])
                p_ib_t = np.interp(tt, t, y[6])
                v_true = float(p_sol_t / max(p_sol_t + p_ib_t, 1e-06))
                unit = '-'
            else:
                continue
            v_obs = v_true * max(0.0, 1.0 + rng.normal(0, sigma_y))
            if name == 'soluble_fraction':
                v_obs = float(np.clip(v_obs, 0.0, 1.0))
            delay_applied = not instant_feedback and turnaround > 0
            rows.append((run_id, batch_id, tt, t_reported, name, v_obs, unit, source, delay_applied, float(turnaround), controls))
    import os
    new = not os.path.exists(path)
    with open(path, 'a', newline='') as f:
        w = csv.writer(f)
        if new:
            w.writerow(['run_id', 'batch_id', 't_true_h', 't_reported_h', 'channel', 'value', 'unit', 'source', 'delay_applied', 'turnaround_h', 'controls_at_sample_json'])
        for r in rows:
            w.writerow(r)
if __name__ == '__main__':
    x0 = (29.0, 0.1, 10.0, 0.1, 30.0)
    P_true, P_obs, traj = run_batch(x0, seed=0, scatter=False)
    print(f'baseline naive  P_sol_true={P_true:.3f}  P_obs={P_obs:.3f}')
    print(f"  Xv_final={traj['Xv_final']:.2f}  A_final={traj['A_final']:.2f}  DO_min={traj['DO_min']:.2f}  T_max={traj['T_max']:.2f}")
    x_de = (25.0, 0.39, 14.1, 0.034, 30.0)
    P_true, P_obs, traj = run_batch(x_de, seed=0, scatter=False)
    print(f'prior DE-like   P_sol_true={P_true:.3f}  P_obs={P_obs:.3f}')
