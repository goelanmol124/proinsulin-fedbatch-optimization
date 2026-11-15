# Optimisation strategy

Goal: maximise final **soluble** proinsulin titre `P_sol(t_end)` over a 7-D operating envelope of a 5 L lab bioreactor, subject to the model in `MODEL.md`. The objective is `P_sol`, not total titre — inclusion-body recovery costs yield and capital, so the optimiser is run against the recoverable fraction.

## 1. Decision variables (7)

`T_set_prod` (°C), `IPTG_dose` (mM), `X_ind` (g/L), `F_s` (L/h), `S_f` (g/L), `N_rpm`, `F_air_vvm`. Bounds and biochemical justification in `MODEL.md` §5.

**Why exactly these seven (was 5):**

- *(a) Biochemical leverage.* The first five are the same biological levers as before: `T_set_prod` decouples growth-vs-folding (only way to exploit the µ/qp temperature-optimum gap and IB-partition); `IPTG_dose` is the Hill knob, now with a real interior optimum once IPTG toxicity entered `kd`; `X_ind` sets pre-vs-post-induction biomass split; `F_s` and `S_f` together control glucose availability and the overflow-acetate threshold.
- *(b) New physical levers — N and F_air.* These were absent in the 5-D problem because DO was a collapsed constant. Now DO is a state, and `kLa = c·N^a·F_air^b` is the operator's real lever against oxygen limitation at high biomass. **N is doubly coupled**: it raises kLa (good, more O2 supply) AND raises mechanical heat ∝ N³ (bad, fights the jacket). The optimiser navigates this trade-off, which is the headline reason we promoted them to decision variables instead of pinning at 600 rpm / 1 vvm.
- *(c) Methodology.* 7-D is still tractable with DE — adds ~2× the simulator calls vs 5-D (popsize=15, maxiter=40 ≈ 600 evals × ~0.2 s each ≈ 2 min wall-clock). Heatmap slices remain interpretable in 2-D pairs once the other five are pinned at the DE optimum.

## 2. Strategy: SciPy `differential_evolution` (DE)

Global, gradient-free, population-based. Default for this problem because:

- The titre surface is **non-convex**: induction-timing creates a ridge; acetate inhibition creates a cliff; oxygen limitation creates a wall in the (Xv, kLa) plane; IPTG toxicity creates an interior optimum that local methods step past.
- **Cheap evaluations** (~150–250 ms per fed-batch sim with `solve_ivp`/BDF for the now-stiffer 9-state system): a population method that needs ~10²–10³ evals is fine.
- **No analytic gradient.** ODE-derived sensitivities (CVODES) exist but adding that machinery violates the dense-code rule and buys little for 7-D box-constrained.
- **Box constraints only.** DE handles them natively.

Settings: `popsize=15`, `maxiter=40`, `tol=1e-3`, `seed=0`, `polish=True` (the final L-BFGS polish recovers the last few % once DE has localised the basin).

## 3. Alternatives considered

| Method | Why not the default | When to switch |
|---|---|---|
| Grid search | Curse of dimensionality: 10⁷ evals at 10 levels per dim. Used only for 2-D visualisation. | Never as primary. |
| Bayesian optimisation (GP) | Sample-efficient but our simulator is cheap; GP overhead competitive with the sim itself. | When wet-lab evaluations replace simulator calls. |
| Nelder–Mead / L-BFGS multistart | Local; misses ridges where DE finds them. Now also misses the IPTG-toxicity interior optimum from naive starts. | Use only as DE's polish phase (already wired). |
| CMA-ES | Comparable to DE on rotated ill-conditioned landscapes. Marginal gain not worth the dependency at 7-D. | If DE plateaus before reaching a known-better titre. |
| RL / dynamic optimisation (open-loop trajectories) | Necessary if `F_s(t)`, `T(t)`, or `N(t)` become time-varying functions instead of scalars. | Phase 2: piecewise-constant feed/temperature profiles. The temperature-state extension makes a real-time T trajectory now meaningful. |

## 4. Visualisation

- **2-D heatmaps** of `P_sol_final` over `(T_prod, IPTG)`, `(X_ind, F_s)`, and the new `(N_rpm, F_air)` slice. The third pair reveals the oxygen-limited wall — runs at low N, low F_air collapse to near-zero titre because DO can't keep up with OUR.
- **Convergence trace** of DE's best `P_sol` vs generation. Monotone curve confirms convergence.

## 5. Current best (auto-updated by notebook)

Last DE run (model as of 2026-05-11, 7-D, popsize=15, maxiter=40):

| variable | optimum | comment |
|---|---|---|
| T_set_prod (°C) | **25.0** | hits lower bound — at low T the IB→soluble shift dominates, and the Gaussian production envelope (sigma_p=4) doesn't drop fast enough to compensate. *Worth a sensitivity check*: tightening sigma_p to 3 or shifting T_opt_p down to 28 would push the optimum interior |
| IPTG_dose (mM) | **0.39** | **interior optimum** — IPTG toxicity (`k_tox_IPTG=0.10`) now balances the Hill saturation. Was at the upper bound 0.5 in the previous 5-D model. This is the single biggest qualitative improvement from the model upgrade |
| X_ind (g/L) | 14.1 | interior — DO-limited growth at 5 L scale flattens the late-phase Xv curve, so "induce later" no longer pays off the way it did at 1 L scale (where it hit 30) |
| F_s (L/h) | 0.034 | interior — high enough to sustain Xv post-induction, low enough to stay below overflow at high biomass |
| S_f (g/L) | 572 | near upper bound (600); concentrated feed minimises dilution |
| N_rpm | 994 | near upper bound (1000); kLa needs to be high to keep DO above ~1 mg/L at peak Xv ≈ 15 g/L. Heat penalty (Q_ag ∝ N³) is real but cooling capacity holds T to ~38.5 °C |
| F_air_vvm | 2.0 | upper bound; foaming-limited in practice — would re-tune downward and add antifoam if model included foaming |

**Final P_sol = 1.41 g/L** (with P_ib = 0.72 g/L → ~66% soluble fraction). Naive setpoint baseline (T=29, IPTG=0.1, X_ind=10, F_s=0.10, S_f=500, N=800, F_air=1.5) gives P_sol = 0.38 g/L → DE finds a **~3.7× improvement** over the naive setpoint.

**Three of seven optima are at bounds** (S_f, N, F_air, all upper). Same pattern as the 5-D version: optimiser rides the physical limits — solubility for S_f, mechanical/foaming limits for N and F_air. T at the lower bound is informative — it says the IB-partition gain at low T is steeper than the µ-and-qp loss within the modelled window. **Two are now genuinely interior** (IPTG, X_ind), which is the qualitative payoff of the upgrade: the previous model's "push IPTG and X_ind to the wall" answers were artefacts of missing IPTG toxicity and unmodelled DO limitation.

Re-run after any model change. If the basin shifts, update this table and record the run.
