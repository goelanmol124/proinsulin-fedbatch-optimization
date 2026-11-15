# Optimisation protocol — the full plan

This document is the *complete, self-contained* blueprint for the model-based optimisation case study. After reading this, an implementer (or the next Claude session) should be able to start coding without asking another design question. It supersedes `OPTIMIZATION.md` as the canonical strategy doc; `OPTIMIZATION.md` becomes a record of the prior 7-D DE run.

Companion docs:
- `MODEL.md` — the wet-lab simulator's ODE system and parameter values
- `VARIABLE_SELECTION.md` — *which* variables we optimise and why (5 core: IPTG, T_prod, µ_set, X_ind, t_f)

---

## 1. Purpose

Demonstrate, with reproducible numerical evidence, that a **model-based optimisation pipeline** beats a typical trial-and-error workflow on the metrics that actually matter: **final product titre** and **number of wet-lab runs spent**. The deliverable is a defensible recommendation backed by a head-to-head benchmark of five optimisation methods against the same simulated wet-lab.

The objective for this study is **only final soluble proinsulin titre `P_sol(t_f)` in g/L**. Robustness, $/g, and multi-objective trade-offs are deliberately left out of the objective for the demo; we can layer them in afterwards without changing the architecture.

---

## 2. Architecture — two models, kept honest

Two pieces of software, separated by an interface that mimics a real bioreactor.

### 2.1 The "wet-lab" — `wetlab_sim.py`

The existing hybrid mechanistic ODE in `insulin_model.ipynb`, wrapped behind a class that **only exposes what a real fermenter would expose**. We (as humans) know the ground-truth equations — that knowledge is for verifying that the optimisation worked, not for the optimiser to consume. **No method below is allowed to import `wetlab_sim`'s internal states directly.** Comparisons are honest because every method passes through the same noisy, sparse, delayed measurement interface.

### 2.2 The "surrogate" — what the optimiser actually knows

Two variants we benchmark independently (§5):

- **UDE surrogate** — mechanistic ODE skeleton (mass balances, Monod, Arrhenius) with small NN residuals on the four most-uncertain terms: qp, k_deg, IB-soluble partition φ, kLa(N,F). Differentiable end-to-end via `diffrax`.
- **GP surrogate** — Gaussian process directly on the operating-point → end-of-run-titre map. No physics, no trajectory.

The UDE is **deliberately mis-specified** relative to the wet-lab (§4.4). The GP starts from a pure prior. The benchmark measures how each closes that knowledge gap from observation.

### 2.3 The mental contract

Every wet-lab run consumes time and money. Every surrogate call is free. The optimisation problem is to spend the budget in the right ratio.

---

## 3. The wet-lab simulator — `wetlab_sim.py`

This section pins down *exactly* what the simulator does, what it returns, and what it accepts. Implementing this from §3 alone should be unambiguous.

### 3.1 Internal state (hidden from optimiser)

The full ODE state vector as already in `MODEL.md`. Plus a small per-run dictionary of **batch-specific biological scatter** drawn once per run from priors:
- µ_max scaled by `LogNormal(1.0, 0.05)` (5% lot-to-lot)
- Yx/s scaled by `LogNormal(1.0, 0.04)`
- qp_max scaled by `LogNormal(1.0, 0.08)` (production rates are noisier)
- kLa coefficient `c` scaled by `LogNormal(1.0, 0.10)` (mixing scatter)
- Yx/OD scaled by `LogNormal(1.0, 0.08)` per run (the dominant OD→Xv error)

This is the "biological" reproducibility ceiling — even with perfect optimisation, two real fed-batch runs at the same setpoint will differ.

### 3.2 Internal control loops (cascade)

The wet-lab simulator implements the following controllers *inside* itself, so the optimiser never has to set RPM or acid pump rates manually. Every controller has a small lag and integrator, modelled as first-order.

| Loop | Setpoint | Manipulated variable | Implementation |
|---|---|---|---|
| **Temperature** | T_set(t) (operator-chosen) | Jacket fluid temperature | First-order lag, τ = 60 s. Step changes complete in ~3 min. |
| **pH** | 7.0 (fixed) | NH₄OH (base) / H₃PO₄ (acid) pump rate | PI controller, deadband ±0.05, dosing in 10 µL pulses. Dilution effect on V tracked. |
| **DO** | 30% (fixed) | Cascade: RPM (300→1000) → airflow (0.5→2 vvm) → O₂ enrichment (0→100%) | PI on RPM; when RPM saturates, ramp airflow; when both saturate, ramp O₂. Cascade switch-over hysteresis ±5% DO. |
| **Foam** | Foam probe high | Antifoam pulse (silicone, 0.2 mL) | Discrete dose, cooldown 5 min between doses. kLa decremented 2% per pulse, recovers over 30 min. |

This means **if the optimiser asks for an operating point that cannot satisfy DO ≥ 30% given the resulting biomass**, the cascade saturates, DO drops below setpoint, and the wet-lab simulator faithfully reports the consequence: acetate spikes, qp collapses, run fails. The optimiser learns this from the data, not from a manually-coded "infeasibility" flag.

This is exactly the **right way** to handle the DO constraint per `VARIABLE_SELECTION.md` §8: it's a side-controller, not a decision variable.

### 3.3 Control inputs the operator can set

These are what the optimiser writes into the wet-lab simulator via the interface:

| Input | Type | When set | Bounds |
|---|---|---|---|
| `T_set(t)` | Piecewise-constant, ≤ 5 segments | Pre-run + event triggers | [20, 42] °C |
| `pH_set` | Scalar | Pre-run | [6.5, 7.4] (we fix at 7.0) |
| `DO_set` | Scalar | Pre-run | [20, 50] % sat (we fix at 30) |
| `feed_schedule` | List of (t_start, mode, params) tuples | Pre-run, mutable on event | mode ∈ {off, constant F, exp(µ_set), pH-stat, DO-stat} |
| `feed_glucose_conc` `S_f` | Scalar | Pre-run | [200, 600] g/L |
| `IPTG_dose(t)` | List of (t_dose, amount) | Pre-run, mutable on event | [0, 1.0] mM cumulative |
| `induction_trigger` | One of {time = T, biomass = X_ind, manual} | Pre-run | — |
| `t_max` | Scalar | Pre-run | [12, 48] h |

The 5 decision variables of `VARIABLE_SELECTION.md` (`IPTG_dose`, `T_set after induction`, `µ_set` in the feed schedule, `X_ind` in the induction trigger, `t_max`) are precisely a strict subset of these.

### 3.4 Online sensors (continuous, noisy)

Emitted every `dt_online = 10 s` while the simulator integrates.

| Channel | Units | Noise model | Drift |
|---|---|---|---|
| T_broth | °C | additive N(0, 0.1²) | none |
| pH | — | additive N(0, 0.02²) | ±0.01 / 24 h linear drift, reset on calibration |
| DO | % sat | additive N(0, 1.5²) | 5% gain drift per 24 h |
| RPM | rpm | reported exactly (encoder) | none |
| Airflow | vvm | additive N(0, 0.02²) | none |
| Vessel weight | kg | additive N(0, 0.005²) | none |
| **Off-gas CO₂** | vol% | additive N(0, 0.01²), 30 s response lag | drift 0.02%/24 h |
| **Off-gas O₂** | vol% | additive N(0, 0.01²), 30 s response lag | drift 0.02%/24 h |

Off-gas analysers (BlueSens BlueInOne Ferm or equivalent) are **standard equipment on 5–10 L pilot bioreactors** — confirmed via vendor docs / 2023–2024 reviews. They are *optional and uncommon* on Ambr250 mini-bioreactors and absent on Ambr15. We include them because we are at 10 L scale, and they give the surrogate **OUR and CER → online µ_hat via stoichiometry**:

```
OUR = (F_air · (xO2,in − xO2,out)) · ρ_O2 / V_broth
mu_hat = OUR / (Xv · Y_X_O2)         with Y_X_O2 ≈ 1.0–1.2 gDCW/gO2 for BL21 on glucose
RQ = CER / OUR                        diagnostic: 0.95–1.05 oxidative, <0.9 IPTG-induced, >1.2 overflow
```

BL21 biomass elemental formula: **CH₁.₇₇O₀.₄₉N₀.₂₄P₀.₀₁₇S₀.₀₀₅**, MW 24.6 g/Cmol (Heijnen & Roels). This is used internally by the cascade controller for OUR-driven O₂ supply estimation, and exposed via off-gas to the surrogate for online µ inference.

### 3.5 Offline assays (sparse, value-error, time-error, turnaround, delay)

Each assay is triggered by an `event` (§3.6). When triggered, the simulator:
1. Records the **true sample time** `t_sample`.
2. Reports a **noisy sample time** `t_reported = t_sample + N(0, σ_t²)` (operator wrote down the wrong minute).
3. Computes the **true measurand** at `t_sample`.
4. Reports a **noisy value** `y_reported = true_value · (1 + N(0, σ_y²)) + bias`.
5. Holds the result until `t_sample + turnaround` before exposing it to the optimiser, **unless `instant_feedback=True` is set** (see §3.7).

| Assay | Min interval | σ_y (CV) | σ_t | Turnaround | Notes |
|---|---|---|---|---|---|
| **OD₆₀₀ → Xv** | 30 min | 10% (10–15% at high OD due to dilution chain + IB refractivity) | ±2 min | 10 min | Linear with run-specific Yx/OD scatter (the per-run draw in §3.1) |
| **DCW** | 2 h | 4% (low biomass: 8–10%) | ±5 min | **26 h** (24 h drying + weighing) | Gold-standard biomass. Result arrives *after the batch is over* in most cases — useful post-hoc for surrogate fitting, useless for in-run feedback |
| **HPLC glucose** (Aminex HPX-87H, RID) | 1 h | 3% | ±5 min | **50 min** (15 min prep + 30 min run + 5 min integration) | |
| **HPLC acetate** (same column) | 1 h | 6% (low-conc bias) | ±5 min | 50 min | Often run on the same injection as glucose |
| **RP-UPLC proinsulin (total)** | 2 h | 7% | ±5 min | **75 min** (IB pellet dissolved in 8 M urea/DTT, filter, run, integrate) | |
| **SDS-PAGE soluble/IB densitometry** | 4 h | **20%** | ±10 min | **5–6 h fast / overnight realistic** | Semi-quantitative. The *only* assay that gives soluble vs IB fraction — and it's noisy and slow |
| **NH₄⁺ (ion assay)** | 2 h | 8% | ±5 min | 90 min | Optional; only if we extend to C:N modelling later |

Several of these were under-estimated in the prior draft. The **DCW result arriving 26 h late** is genuinely how things work — DCW is a post-hoc verification number, not an in-run control variable. Our surrogate uses OD₆₀₀ in-run and corrects to DCW after the batch closes.

### 3.6 Event triggers (sampling and control)

A run is driven by an **event queue**. The simulator integrates until the next event fires, pauses, emits the relevant measurements (with turnaround), accepts new control inputs, and resumes.

Trigger types:
- `time_grid(dt)` — every `dt` hours
- `time_after_event(event, dt)` — `dt` after another event (e.g. 6 h after induction)
- `DO_spike(threshold, slope)` — DO rises above `threshold` with slope > `slope` (feed-onset signal)
- `biomass_threshold(X)` — OD₆₀₀-inferred Xv crosses `X` (induction trigger)
- `acetate_threshold(A_max)` — HPLC acetate crosses `A_max` (early-stop or feed-throttle)
- `clock_alarm(t)` — fixed wall-clock t (induction by time)
- `manual` — operator-specified one-offs

Each event carries:
- The list of assays to take
- Whether sample-to-action delay is applied (`instant_feedback` flag, §3.7)
- The control-input updates to apply on resume

### 3.7 The sample-to-action delay flag

Real fed-batch operations have a 2–6 h manual delay (sample → QC queue → operator decision → actuation), or 20–40 min with autosampler / online HPLC. We model both, **and we expose a `instant_feedback` flag per event** so that:

- For sim-only benchmarking (DE on UDE surrogate, no real wet-lab loop) → delay irrelevant; not applied.
- For sim-fast development → set `instant_feedback=True` to skip delay and iterate faster.
- For the **honest benchmark of methods that close the loop on the wet-lab** → set `instant_feedback=False`. Delay applied. Effective realisation: the controller can only act on the *t_reported − turnaround − delay* state, not the current state. This is a real and important degradation of performance vs ideal.

Default in benchmark runs: `instant_feedback=False` for offline assays (realistic), `True` for online sensors (already at 10 s, no extra delay meaningful).

### 3.8 Output — `batch_NNNN.csv` per run

Every batch (whether simulated or "wet-lab") produces a CSV with this schema. **Long format, one row per measurement**, so missing samples don't create sparse columns.

```
run_id, batch_id, t_true_h, t_reported_h, channel, value, unit, source,
  delay_applied, turnaround_h, controls_at_sample_json
```

Where:
- `run_id` — unique per execution (e.g. UUID4)
- `batch_id` — within-run batch index (1..N for parallel batches)
- `t_true_h` — true sample time (for analysis; in real life unknown, kept here for benchmark accounting)
- `t_reported_h` — what the operator wrote down
- `channel` — one of: `T_broth, pH, DO, RPM, F_air, weight, OD600, DCW, glucose, acetate, proinsulin_total, soluble_fraction, OUR, CER, RQ`
- `source` — `online | offline_OD | offline_DCW | offline_HPLC | offline_PAGE | offgas`
- `controls_at_sample_json` — current operator inputs (T_set, feed mode + params, IPTG cumulative, etc.)

A `run_summary.csv` at the top of the run logs: chosen operating point, final P_sol(t_f) **(this is the ground-truth oracle — only used for post-hoc benchmark scoring, NEVER fed back into the optimiser**), wall-clock simulated, and convergence flags.

---

## 4. The surrogate model — `optimization.py`

### 4.1 Why two variants — UDE and GP

We benchmark them **independently** rather than always combining, because the user (rightly) wants to know which one carries the weight. The combination is also tested but as a separate method.

### 4.2 UDE surrogate

A differentiable ODE in `diffrax` (JAX) implementing the mechanistic skeleton from `MODEL.md` **with deliberate omissions and biases** so it differs from the wet-lab (§4.4). Four small NNs replace specific rate terms:

| Term replaced | Why this term | NN size | Inputs |
|---|---|---|---|
| qp(Xv, S, T, IPTG, t_since_ind) | The most uncertain rate; literature priors are noisy | 3 layers × 16 units | log(Xv), log(S), T, IPTG, t_ind |
| k_deg(T, IPTG) | Proteolytic loss; small empirical literature base | 2 × 8 | T, IPTG |
| φ_soluble(T, IPTG, qp) | IB-vs-soluble partition; semi-empirical at best | 2 × 8 | T, IPTG, qp |
| kLa(N, F_air, V) | Geometry-dependent; literature correlations vary 2× | 2 × 8 | log N, log F_air, log V |

Trained by gradient through the ODE on observed trajectories (OD, glucose, acetate, off-gas-inferred OUR, proinsulin total, soluble fraction). Posterior parameter uncertainty via **mean-field variational inference** (lightweight) or **HMC** (gold standard, optional).

### 4.3 GP surrogate

`BoTorch SingleTaskGP` on the 5-D operating-point space → end-of-run titre, with Matern-5/2 kernel and observation-noise std fit from data. Heteroscedastic option: noise scales with predicted titre (high titre → wider HPLC CV).

### 4.4 Deliberate mis-specification of the surrogate vs the wet-lab

This is **critical for honest benchmarking**. If the surrogate is identical to the wet-lab, optimisation is trivial and meaningless.

Mis-specification we apply to the **initial UDE**:
- **Wrong qp_max prior**: literature prior centred at 0.5 × the wet-lab's true value, with broad σ (factor of 2). The NN must learn the bias.
- **Missing acetate→qp coupling**: the initial UDE has no acetate inhibition on qp; the wet-lab does. The NN on qp must discover it from trajectories.
- **Wrong T_opt_p**: prior on T_opt_p = 32 °C, wet-lab truth = 28 °C. Forces optimisation to push lower than the prior wants.
- **Bias in kLa correlation**: prior coefficients off by 25%. Real DO behaviour will differ from predicted, especially at high biomass.
- **No IB→soluble dynamics**: initial UDE assumes all product is soluble. SDS-PAGE data must drive the φ NN to learn partitioning.

These five biases together mean the *initial* UDE prediction at the wet-lab's true optimum is wrong by ~40% in titre. After Bayesian updates from wet-lab batches, this gap should close.

### 4.5 Where is the UDE most uncertain? — active-learning machinery

For active learning (§5.5), the surrogate exposes `posterior_variance(x)` for any candidate operating point `x`. Implementations:

- **UDE**: Monte Carlo over the variational posterior of NN weights — sample 32 posterior draws, simulate each, return std of `P_sol(t_f)`.
- **GP**: native `posterior(x).variance`.

The active-learning acquisition function picks the next wet-lab batch as:

```
x_next = argmax  α · EI(x)                                    # exploit good titre
              + β · qMaxPosteriorVariance(x | UDE)             # explore where surrogate disagrees with itself
              + γ · qKG_under_cost(x)                          # knowledge gradient with cost weighting
```

with α, β, γ scheduled — high β early (explore) → high α late (exploit). For parallel batches (k=4), use `qNoisyExpectedImprovement` from BoTorch, which natively handles batch acquisition.

This directly answers the user's question — *yes, there is a principled way to know where the UDE is most unsure* (MC posterior variance), *and we can route the next wet-lab there*.

---

## 5. The benchmark — four realistic methods + one oracle (parked)

Dropped from the prior 8: Nelder-Mead (superseded by DE), CMA-ES (DE adequate at 5-D), IPOPT-collocation (no trajectory variables in scope; deferred to Phase 2). Of the five that were active, **`method_de_sim` (§5.3) has been demoted to an oracle reference and removed from the `METHODS` dict** (see the run log 2026-05-13 "honest ranking" entry). It optimised against the wet-lab simulator with `scatter=False, noisy=False` — i.e. against a perfect copy of ground truth — which is not a method any real campaign can execute. The function is left in `optimization.py` so it can be reinstated as a ceiling-vs-realistic comparison when needed; for the headline benchmark the four realistic methods (5.1, 5.2, 5.4, 5.5) are what stand.

### 5.1 Method 1: Directed trial & error (LHS scouts + local perturbation)

**Mimics current lab practice properly.** Latin hypercube of 25 scouts, then 75 Gaussian-perturbation refinements around the running best (σ = 8 % of each variable's range). Captures what experimentalists actually do — drift the best point with small one-off probes. No model, no learning. **Budget: 100 wet-lab runs** (was 20 LHS-only; upgraded 2026-05-13 fairness sweep).

### 5.2 Method 2: DoE + RSM (sequential Box-Behnken + quadratic + DE optimum)

**Mimics the modern incumbent** with a **two-stage** response-surface workflow (closer to classical sequential RSM than a single global quadratic over the full box). **Stage 1:** Box-Behnken for 5 factors = 46 runs on global `BOUNDS`; ordinary least squares quadratic in coded BBD coordinates; **differential evolution** maximizes the fitted surface on `[0,1]^5` (same as the implementation — not a closed-form ∂P/∂x = 0 solve). **Shrink:** define a new axis-aligned box of length `shrink` × global span (default 0.45), **centered on the best *observed* Stage-1 batch** (max `Po` among the 46 runs), clipped to `BOUNDS`, with a minimum edge length 0.15× global span so the design cannot collapse at a face. **Stage 2:** second 46-run BBD **only** on the reduced box; **refit** quadratic **from Stage-2 data only** in `[0,1]^5` on that box; DE → quadratic recommendation. **Belt-and-braces recommendation:** also score the best-observed batch across both stages; pick whichever has higher observed `Po` (guards the failure where the Stage-2 quadratic extrapolates outside its data envelope). **Confirmation:** replicates at the recommendation until total wet-lab batches = 100 (typically 5–8 replicates). **Budget: 100 wet-lab runs** (was 95; padded to 100 for the 2026-05-13 fair-budget sweep).

Expected weakness, predicted in advance: the IPTG × T_prod interaction is **not quadratic** (Hill saturation + Gaussian envelope). RSM is consistent (very tight IQR — it always lands ~2.15 g/L) but rarely breaks through to the 2.7 g/L plateau that directed trial-and-error reaches at the same budget.

### 5.3 Method 3 (parked): DE on noise-free wet-lab simulator — oracle ceiling

Originally specified as "DE on UDE surrogate, sim-only". The implementation actually runs DE directly on the wet-lab simulator with scatter and noise disabled, which is **the ceiling that an unbiased model would reach**, not a realistic method. Demoted out of the headline `METHODS` set. The function `method_de_sim` remains in `optimization.py` and can be added back to `METHODS` to print the ceiling line on plots. A *realistic* Method 3 (DE on the biased low-fidelity surrogate from §4) would land at ~1.2–1.3 g/L wet-lab — i.e. tied with GP-BO — and is not separately implemented because MFBO+UDE (§5.5) strictly dominates it (same surrogate + closed-loop wet-lab refinement).

### 5.4 Method 4: GP-BO directly on wet-lab (no UDE)

**The modern BO baseline.** Initial LHS of 5 runs, then 94 EI iterations + 1 final confirmatory run at the posterior-mean argmax = **100 wet-lab runs** (was 21; bumped 2026-05-13 fair-budget sweep). Acquisition candidate pool refreshed each iteration as 1500 uniform + 800 Gaussian-around-running-best — avoids the static-pool failure mode where EI exploits past the original candidate set.

### 5.5 Method 5: Multi-fidelity BO + UDE (proposed champion)

**The full pipeline.** Pseudocode:

```
# Stage 0 — biased priors
ude = build_ude(biased_priors)
budget_wetlab = 0

# Stage 1 — sim-only warm-start (cheap)
x_warm = DE_on_ude(ude, n_evals=600)

# Stage 2 — initial wet-lab anchor
xs_init = LHS(5) ∪ {x_warm}
batches = parallel_wetlab(xs_init, k=6)
ude = update_ude_bayesian(ude, batches)
budget_wetlab += 6

# Stage 3–N — closed loop
while not converged:
    x_next = argmax( α·qEI + β·qMaxPosteriorVar + γ·qKG )  # batch of k=4
    batches = parallel_wetlab(x_next, k=4)
    ude = update_ude_bayesian(ude, batches)
    budget_wetlab += 4
    check_stopping(titre_plateau AND posterior_var_at_best < threshold)
```

For the fair-budget benchmark this runs to **100 wet-lab batches** (6 anchor + 93 closed-loop + 1 confirm), β decay slowed to 0.02/iter so explore weight stays meaningful past iter 12, and 80 Gaussian-around-best candidates appended to the LF-memoised pool every 5 iterations. The original convergence stopping rule (best titre improves <2% over two rounds AND posterior std at recommendation <5% of mean) is what *would* be used in practice; we ignore it here to give every method the same wet-lab quota. Empirically MFBO+UDE's win is **per-batch efficiency at low budget** (≤20 wet-lab) — at 100 batches the simpler methods catch up.

### 5.6 Per-call wall-clock budget (measured, post-implementation)

| Call type | Measured time | Notes |
|---|---|---|
| One wet-lab fed-batch (BDF, 9-state, scatter+noise) | 0.08–0.15 s | `wetlab_sim.run_batch` |
| One biased-LF fed-batch (LSODA, loose tol) | 0.08–0.12 s | `optimization.low_fidelity` |
| One GP residual fit on n≤20 points (sklearn) | 0.05–0.15 s | `MFSurrogate.fit` |
| One full DE on biased LF (8×10 evals) | 8–12 s | warm-start in MFBO+UDE |

Per-method wall-clock at the implemented budgets (single seed):
- `trial_error` (20 wet-lab) — **~1.7 s**
- `doe_rsm` (2-stage, 95 wet-lab) — **~9 s**
- `gp_bo` (5 init + 15 EI iter + 1 confirm = 21 wet-lab) — **~5 s**
- `mfbo_ude` (warm-start DE + 5 init + 12 closed-loop iter + 1 confirm = 19 wet-lab) — **~20 s**

A 4-method × 2-seed bench finishes in ~75 s. A 4-method × 20-seed bench is ~12–15 minutes. Tractable on a laptop.

---

## 6. Reproducibility — proving the champion is robust

Single-seed comparisons are anecdote. Run the **whole benchmark across 20 random seeds** controlling: (a) the per-run biological scatter (§3.1), (b) measurement noise realisations, (c) optimiser stochastic seeds, (d) initial LHS designs.

Report per method:
- Per-seed final titre histogram + median (IQR)
- **Convergence trace overlay** — best-seen titre vs wet-lab-run-count, all 20 seeds. The champion should show tight band; trial-and-error a wide one.
- **Final-titre distribution** — does the champion reach the same optimum across seeds? Pass criterion: 18/20 seeds within 5% of the seed-1 optimum.
- Wilcoxon signed-rank: per-seed paired comparison of champion vs each baseline. Claim superiority only when p < 0.01 after Bonferroni (4 comparisons).

This directly answers the user's "multiple times when we ran the system it reached to same optimum in fair amount of time" requirement.

---

## 7. The metrics table (what we report)

For each method, mean ± std across 20 seeds:

| Metric | Description |
|---|---|
| **Final P_sol** | Wet-lab oracle titre at recommended operating point |
| **Wet-lab runs to within 5% of best-seen** | The headline efficiency number |
| **Total wall-clock to recommendation** | Calendar time (sim + sample + turnaround simulated realistically) |
| **Simulator CPU-seconds** | Trivial for sim-only methods, dominates for surrogate-heavy ones |
| **Reproducibility ratio** | #seeds within 5% of median optimum / 20 |
| **Surrogate posterior std at recommended point** | Confidence in the recommendation (where applicable) |

A **headline plot** of `median best-seen titre vs wet-lab runs` with 20-seed bands. A **table** with the metrics. A **two-row Pareto plot** (titre vs runs, and titre vs sim CPU). Together these are the supervisor-facing artefacts.

---

## 8. Constraints — how each is handled

| Constraint | Type | Where it lives |
|---|---|---|
| Variable box bounds | Box | Native in DE / RSM / BO / DoE |
| `V_f ≤ 0.8 · V_max` | Linear (V_f = V_0 + F̄·(t_f − t_switch)) | Penalise infeasible (F̄·t_f) combos via clipping in the optimiser before evaluation |
| `DO ≥ 30% throughout` | Nonlinear (state) | **Handled by the wet-lab simulator's cascade controller** — if infeasible, DO drops, run fails, optimiser sees low titre. No explicit constraint passed to optimiser. |
| `acetate ≤ 2 g/L` (soft) | Nonlinear (state) | Same — emerges from µ_set, no explicit constraint. Confirmed via end-of-run HPLC acetate field in the CSV. |
| `pH ∈ [6.8, 7.2]` | State | Cascade controller; deviation reported but never penalised (controller holds it). |
| `T_set ∈ [25, 37]` post-induction | Box | Native. |

**Critical:** the cascade-based handling of DO is what makes the optimiser's life simple — it's a 5-D box-constrained problem. If we tried to expose DO as a constraint to the optimiser, we'd need IPOPT or constrained-EI and the comparison would be muddier. The architecture *encodes* `VARIABLE_SELECTION.md` §8's principle ("solve constraints algebraically, optimise over the rest").

---

## 9. Layered / hierarchical optimisation (per `VARIABLE_SELECTION.md` §9)

This benchmark sits at **Layer 2** (operating policy). Layers 0–1 (structural choices, media composition, defaults) are already pinned by `VARIABLE_SELECTION.md` §7. Layer 3 (local polish) is appended after the global stage for each method (L-BFGS on UDE for Methods 3 and 5; not applicable to Methods 1, 2, 4).

We do **not** mix continuous + discrete optimisation; carbon source and feed mode are pre-pinned. This keeps every method in its native regime (continuous box-constrained 5-D).

---

## 10. Repository layout

```
insulin/
├── insulin_model.ipynb         (existing — the original mechanistic model + DE optim)
├── wetlab_sim.py               (§3, the simulator with realistic measurement interface)
├── optimization.py             (§4–§5: biased LF + MF GP residual, five methods, `run_all` / CLI → `results.npz`)
├── diagnostics.py              (plots + Wilcoxon from `results.npz`)
├── analysis.ipynb              (NEW — generates the headline plot + metrics table)
├── runs/
│   ├── run_<uuid>/
│   │   ├── batch_0001.csv      (per-batch trajectory + offline data, §3.8)
│   │   ├── batch_0002.csv
│   │   └── run_summary.csv
│   └── ...
├── docs/
│   ├── VARIABLE_SELECTION.md   (existing — variable choice rationale)
│   ├── OPTIMIZATION_PROTOCOL.md (THIS DOC — the strategy)
│   ├── OPTIMIZATION.md         (becomes legacy — old 7-D DE run, retained)
│   ├── MODEL.md                (existing)
```

---

## 11. Implementation roadmap

Concrete, ordered. Each step is one work session.

1. **`wetlab_sim.py` skeleton.** Wrap existing ODE. Add per-run biological scatter draws (§3.1). Verify it reproduces `insulin_model.ipynb` results at a known operating point with `scatter=False`.
2. **Cascade controllers.** Implement T, pH, DO, foam loops (§3.2). Verify by running at a known-good operating point and checking DO holds 30% ± 3 across the run.
3. **Online sensor noise + off-gas.** Implement §3.4. Generate one batch CSV and inspect.
4. **Offline assay model.** Implement §3.5 with all the corrected numbers (50 min HPLC, 5–6 h SDS-PAGE, 26 h DCW, etc.). Verify turnaround delay logic on a fake event queue.
5. **Event queue & `instant_feedback` flag.** Implement §3.6–§3.7. Test with a simple "induce at OD=10, sample every hour" policy.
6. **CSV writer.** Implement §3.8. One sample CSV checked into `runs/example/`.
7. **`optimization.py` — surrogate stack (§4).** Biased LF + sklearn GP residual (implemented); optional future: diffrax UDE + BoTorch GP as in §4.2–§4.3.
8. **Method 1 — trial & error.** Implement, run, log.
9. **Method 2 — DoE + RSM.** Box-Behnken pattern in `optimization.py` (or `pyDOE2` if preferred); `numpy.linalg.lstsq` quadratic; **two-stage** BBD with shrink box; maximize fitted surface with `scipy.optimize.differential_evolution` (not an analytic grad=0 solve).
10. **Method 3 — DE on UDE.** Reuse the existing DE setup; swap ODE for UDE.
11. **Method 4 — GP-BO direct.** BoTorch `qNoisyExpectedImprovement`.
12. **Method 5 — MFBO + UDE.** The champion. Closed-loop driver (§5.5).
13. **`optimization.py` as `__main__`.** Loop methods × seeds, write `results.npz` (and extend to CSV/`results.parquet` if needed).
14. **`analysis.ipynb`.** Headline convergence plot + metrics table + Wilcoxon test + Pareto plot.
15. **Final write-up.** Numbers, plots, one-page summary.

Each step is independently testable. We can stop at any point and we'd still have the artefacts from the prior step.

---

## 12. Success criteria — what counts as "proved"

Pre-registered, so we cannot move the goalposts:

1. **Method 5 (MFBO+UDE) median final titre ≥ Method 2 (DoE/RSM) median final titre**, statistically (Wilcoxon p < 0.01, Bonferroni-corrected).
2. **Method 5 uses ≤ 1/3 the wet-lab runs of Method 2** to reach within 5% of its own final titre. (With Method 2 at **95** runs, one-third is ≈32; logged Method 5 budgets ~19 still satisfy the inequality.)
3. **Method 5 reproducibility ratio ≥ 0.9** (18/20 seeds within 5% of seed-1 optimum).
4. **Method 5 finds the IPTG interior optimum**; Method 2 (RSM) finds either the box edge or the wrong side of the Hill. Visible in the IPTG posterior.

If any of these fail, the recommendation is honest: the simpler method wins on its merits. The point of the protocol is to **let the data decide**, not to engineer a forgone conclusion.

---

## 13. What this is NOT doing (deliberately)

- **No robustness / cost in the objective.** Pure titre, by user request. Robustness Σ-uncertainty Monte Carlo and $/g cost can be added as a multi-objective extension later without rewriting any code (the architecture supports it).
- **No trajectory variables.** T(t), F(t), IPTG(t) are piecewise-constant with operator-set breakpoints. True time-varying optimisation (collocation + IPOPT) is Phase 2.
- **No real wet-lab.** Everything runs on the digital twin. Outcome of the demo: a memo to the supervisor showing the methodology works; *then* the supervisor commits ~$5–15k of real wet-lab capacity to confirm.
- **No C:N modelling.** Per `VARIABLE_SELECTION.md` §6, deferred.
- **No multi-strain comparison.** BL21 + proinsulin only.

---

## 14. What an implementer (or next Claude session) needs to start

Nothing beyond this doc + `MODEL.md` + `VARIABLE_SELECTION.md`. Every parameter, every interface, every method choice, every metric is pinned. The user can clear the conversation and prompt "start step 1 of the OPTIMIZATION_PROTOCOL.md roadmap" and work resumes cleanly.

---

## TL;DR

Wrap the existing model as a noisy, delayed digital twin (`wetlab_sim.py`). Build a deliberately-biased UDE surrogate plus a vanilla GP. Benchmark **five** methods (trial-and-error, DoE/RSM, DE-on-UDE, GP-BO-direct, MFBO+UDE) against the twin across 20 seeds. The champion (MFBO+UDE) is expected to match DoE/RSM titre at <1/3 the wet-lab cost with reproducible convergence. Every batch produces a standard CSV. Pre-registered success criteria. No new information needed to start.
