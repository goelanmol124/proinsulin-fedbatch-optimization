# Proinsulin fed-batch model

Mechanistic ODE model of recombinant *E. coli* (BL21-class) producing proinsulin under IPTG induction in a **5 L lab-scale stirred-tank fed-batch bioreactor**. Predicts **soluble product titre P_sol(t_end) (g/L)** as a function of operating inputs, and is wrapped by an optimiser (`OPTIMIZATION.md`).

## 1. Why mechanistic, not data-driven

We have no in-house dataset. Literature gives well-constrained kinetic and stoichiometric constants (see §4), so a Monod–Luedeking–Piret framework with overflow-metabolism, dissolved-oxygen mass transfer, energy balance, and induction kinetics extrapolates beyond a fit window — important because the optimiser will probe the corners of the input space, where a black-box ML model trained on a narrow dataset would lie. A hybrid layer (e.g. a Gaussian-process residual on `qp`) is a natural extension once experimental data exist; the current `qp` is left isolated so this is a one-line swap.

## 2. State variables (9)

The previous 6-state version assumed perfect DO and T control. That assumption breaks in two regimes the optimiser actively explores: high biomass (where kLa cannot meet OUR and DO crashes regardless of setpoint) and high N (where mechanical heat exceeds jacket capacity and T overshoots). Both states are now real.

| Symbol | Meaning | Units |
|---|---|---|
| V     | Reactor working volume                       | L |
| Xv    | Viable biomass                                | g/L |
| Xd    | Dead biomass                                  | g/L |
| S     | Glucose                                       | g/L |
| A     | Acetate                                       | g/L |
| P_sol | Soluble (correctly folded) proinsulin         | g/L |
| P_ib  | Inclusion-body (misfolded) proinsulin         | g/L |
| DO    | Dissolved oxygen                              | mg/L |
| T     | Broth temperature                             | °C |

Soluble vs IB is split because a 24 h run that produces 4 g/L total at 37 °C with 90 % IB is operationally worse than 2 g/L at 28 °C with 80 % soluble — IB recovery requires solubilisation/refolding, which costs yield and capital. The optimiser maximises **P_sol**, not total titre.

IPTG is *not* a state; the bolus instantly mixes and dilution by feed (<3 %/h) is small over the post-induction window. It is computed algebraically as `IPTG = IPTG_dose · σ(Xv − X_ind)`. Re-introduce as a state if you model intracellular vs extracellular inducer.

pH is held at 7.0 by base addition (controller folded out of the model).

## 3. Equations

Dilution `D = F_active/V` (only feed adds volume; evaporation and base flow are <1 %/h, folded into the controller).

**Feed activation** is glucose-on-demand: `F_active = F_s · σ(1−S, 0, 0.3)` — feed turns on when S drops below ~1 g/L. This is the industry-standard *E. coli* protocol; constant feed from t=0 lets glucose accumulate past the Haldane peak before consumption catches up, stalling growth indefinitely.

**Discrete switches → smoothed sigmoids.** Induction trigger and feed switch are realised as `σ(x, c, w) = 1/(1+exp(−(x−c)/w))` with small `w` rather than hard `if`-conditions. Hard switches make the ODE non-Lipschitz; BDF/LSODA bog down or fail. Mathematically equivalent in the limit `w → 0`; numerically essential.

```
dV/dt     = F_active
dXv/dt    = (mu·phi − kd)·Xv − D·Xv
dXd/dt    =  kd·Xv − D·Xd
dS/dt     =  D·(S_f − S) − Xv·qs
dA/dt     =  (q_ac_over − q_ac_up)·Xv − D·A
dP_sol/dt =  qp·f_sol·Xv      − D·P_sol
dP_ib/dt  =  qp·(1 − f_sol)·Xv − D·P_ib
dDO/dt    =  kLa·(DO* − DO) − OUR_vol − D·DO
dT/dt     =  (Q_met + Q_ag − Q_cool + H_feed) / (V·rho·Cp)
```

### Specific growth rate `mu`

```
mu = mu_max · fT_g(T) · S/(Ks + S + S²/Ki_S) · DO/(Kdo + DO) · Ki_A/(Ki_A + A)
fT_g(T) = exp(-0.5·((T − T_opt_g)/sigma_g)²)
```

- **Haldane** in S, not pure Monod: substrate inhibition matters at the start of the feed when residual glucose can spike.
- **Monod in DO** is now real, not collapsed — DO is a state, and `Kdo = 0.2 mg/L` means growth crashes once kLa fails to keep DO above ~1 mg/L.
- **Non-competitive acetate**: same form reused in `qp`.
- **Gaussian temperature envelope** evaluated on the *actual* T (not setpoint), so jacket-tracking error feeds back into growth.

### Specific production rate `qp`

```
qp = (alpha·mu + beta) · IPTG^n/(K_IPTG^n + IPTG^n) · Ki_A/(Ki_A + A) · fT_p(T)
fT_p(T) = exp(-0.5·((T − T_opt_p)/sigma_p)²)
```

- `beta ≈ 0` for plasmid-based BL21; kept as a knob for stronger constitutive promoters.
- **Hill** with n=2 captures cooperative LacI–IPTG response; n=1 over-predicts low-IPTG titre.
- **Crucially `mu` and `qp` have different temperature optima** — that's the entire reason a two-phase temperature profile beats isothermal operation.

### Metabolic burden `phi`

```
phi = 1 − (1 − phi_min) · σ(Xv − X_ind) · IPTG^n/(K_IPTG^n + IPTG^n)
```

Pre-induction `phi = 1` (no recombinant load). Post-induction it drops toward `phi_min = 0.4` as the Hill term saturates with IPTG dose. Was a constant 0.5 before; now dynamic so weak induction yields weak burden — necessary to capture the IPTG-dose response correctly. Ties three pieces together: induction biomass, inducer dose, and growth penalty.

### Overflow metabolism / acetate

```
qs       = mu/Yxs + ms + qp/Yps                      # total specific glucose uptake
q_ac_over = Y_AS · max(0, qs − qs_crit)              # only the overflow above oxidative capacity
q_ac_up   = k_au · A/(K_au + A) · (1 − σ(S, 1, 0.3)) # re-uptake when glucose is depleted
```

`qs_crit = 0.4 g/g/h` is the oxidative capacity of BL21 — above it, glucose is shunted to acetate even with full DO. **The single biggest operational lever** between a 1 g/L and a 5 g/L proinsulin campaign. Acetate re-uptake is load-bearing: without it the model predicts >30 g/L acetate, which never occurs experimentally. Set `k_au = 0` to model an acetate-reuptake-deficient mutant.

### Inclusion-body partitioning

```
f_sol = 1/(1 + exp((T − T_IB_50)/sigma_IB)) · 1/(1 + (mu/mu_crit_IB)²)
```

Two independent bottlenecks for correct folding: temperature (chaperone activity drops above ~32 °C, with a steep sigmoid) and growth rate (high `mu` overwhelms folding capacity even at low T). They are multiplicative — relieving one doesn't help if the other is saturated. This is what makes the temperature downshift at induction worth doing operationally: it doesn't increase total `qp` (the production term still drops with T below 29 °C), it shifts the *partition* toward soluble.

### Death rate

```
kd = kd_0 + k_tox_A·A + k_tox_IPTG·IPTG
```

Three terms: basal turnover, acetate-driven (lowers intracellular pH, dissipates PMF), and **IPTG-driven**. The IPTG term is new and closes the model gap that previously sent the optimiser to IPTG = 0.5 mM (hits Hill saturation, no penalty). With `k_tox_IPTG = 0.10 1/h/mM`, 0.5 mM doubles `kd`, restoring the literature-typical 0.05–0.1 mM optimum. Arrhenius denaturation dropped because we operate within ±10 °C of optimum; re-add if scale-up exposes hot spots.

### Oxygen uptake & mass transfer

```
q_O2     = mu/Yxo + mo + qp/Ypo                  # specific O2 demand, g_O2/g_X/h
OUR_vol  = q_O2 · Xv · 1000                       # mg_O2/L/h
kLa      = kLa_c · N_rpm^a · F_air_vvm^b          # 1/h, lab-scale Rushton correlation
```

`a = 1.5, b = 0.5, kLa_c = 0.027` is tuned so `N=600 rpm, F_air=1 vvm` gives `kLa ≈ 400 1/h` (typical 5 L Rushton). The exponents are mid-range of published correlations; different impellers (axial-flow, marine) shift `a` toward 2.0. DO* = 7.5 mg/L (water/broth saturation at 37 °C, 1 atm).

This is what makes N and F_air decision variables: at high biomass, OUR ≈ 5 g_O2/L/h ≈ 5000 mg/L/h, and `kLa·DO* ≈ 3000 mg/L/h` at 600 rpm — so the run becomes oxygen-limited unless N or F_air is increased. The optimiser navigates this boundary.

### Energy balance

```
Q_met   = H_O2 · (OUR_vol/32) · V                       # kJ/h, OUR mg/L/h → mmol/L/h via /32
Q_cool  = max(0, K_jacket · (T − T_set))                # proportional jacket, cooling-only
Q_ag    = Q_ag_ref · (N_rpm/600)^3                       # mechanical heat ∝ N³ (impeller power law)
H_feed  = F_active · rho · Cp · (T_feed − T)             # sensible heat of cooler feed
dT/dt   = (Q_met + Q_ag − Q_cool + H_feed) / (V·rho·Cp)
```

`K_jacket = 200 kJ/h/K` is a finite proportional cooling capacity — large enough to track setpoint within ~0.3 °C under typical load, but finite, so a high-OUR run sees real overshoot. This makes `N` a coupled lever: high `N` raises `kLa` (good for DO) AND raises `Q_ag` ∝ N³ (bad for T tracking), and high biomass raises `Q_met` even more — both feed back into `mu`, `qp`, and `f_sol` via the Gaussian envelopes. The two-phase setpoint (37 °C → T_set_prod) is encoded as `T_set = T_opt_g + (T_set_prod − T_opt_g)·σ(Xv − X_ind)`.

## 4. Parameter table

| Param | Value | Units | Source / why |
|---|---|---|---|
| mu_max | 0.6 | 1/h | mid of 0.4–0.7 reliable BL21 range; 0.9 only under non-burdened conditions |
| Ks | 0.05 | g/L | standard BL21 value; do not raise above 0.1 |
| Ki_S | 80 | g/L | substrate inhibition, only active during feed spikes |
| Kdo | 0.2 | mg/L | DO half-sat; now load-bearing because DO is a real state |
| Ki_A | 8 | g/L | within the 5–10 g/L consensus; never lower below 5 |
| kd_0 | 0.02 | 1/h | basal death |
| k_tox_A | 0.005 | L/g/h | tuned so kd doubles at A=4 g/L |
| k_tox_IPTG | 0.10 | 1/h/mM | linear IPTG toxicity; 0.5 mM doubles kd. Calibrated to recover the 0.05–0.1 mM literature optimum |
| Yxs | 0.45 | g/g | glucose-grown E. coli on minimal+complex media |
| ms | 0.03 | g/g/h | maintenance |
| Yxo | 1.1 | g_X/g_O2 | mid 1.0–1.2 |
| mo | 0.4 | g_O2/g_X/h | ≈12 mmol/g/h × 32/1000; baseline maintenance O2 |
| Ypo | 0.5 | g_P/g_O2 | order-of-magnitude; small term in OUR |
| alpha | 0.12 | g/g | mid of optimised 0.02–0.15 g/g range |
| beta | 0 | g/g/h | growth-coupled production |
| K_IPTG | 0.05 | mM | within 0.01–0.1 toxicity-balanced range |
| n_Hill | 2 | – | LacI cooperativity |
| T_opt_g | 37 | °C | growth optimum |
| sigma_g | 5 | °C | width of Gaussian envelope |
| T_opt_p | 29 | °C | folding-favoured production optimum |
| sigma_p | 4 | °C | tighter envelope: production drops faster off-target than growth |
| qs_crit | 0.4 | g/g/h | oxidative capacity threshold |
| Y_AS | 0.4 | g/g | acetate yield on overflowed glucose |
| k_au | 0.10 | 1/h | acetate re-uptake when glucose is depleted; load-bearing |
| K_au | 0.5 | g/L | acetate uptake half-sat |
| Yps | 0.15 | g/g | proinsulin yield on glucose |
| T_IB_50 | 32 | °C | T at which 50% of nascent product misfolds (chaperone capacity) |
| sigma_IB | 2 | °C | steepness of T-driven IB onset |
| mu_crit_IB | 0.30 | 1/h | growth rate at which folding capacity is half-saturated |
| phi_min | 0.4 | – | post-induction floor of the burden multiplier |
| kLa_c | 0.027 | – | tuned so N=600 rpm, F_air=1 vvm → kLa≈400 1/h |
| kLa_aN | 1.5 | – | agitation exponent; mid of 1.5–2.0 lab-scale range |
| kLa_bF | 0.5 | – | aeration exponent; mid of 0.5–1.0 |
| DO_star | 7.5 | mg/L | O2 saturation in broth |
| rho | 1.0 | kg/L | broth density (water-like) |
| Cp | 4.18 | kJ/kg/K | broth heat capacity |
| K_jacket | 200 | kJ/h/K | finite proportional cooling — high-OUR runs see real overshoot |
| Q_ag_ref | 10 | kJ/h | gassed Rushton at 5 L scale, 600 rpm reference; scales as (N/600)³ |
| H_O2 | 0.46 | kJ/mmol | universal metabolic heat constant |
| T_feed | 25 | °C | feed reservoir temperature |
| V0 | 2.5 | L | 50% fill of 5 L vessel; leaves room for ~2.5 L of feed |
| Xv0 | 0.1 | g/L | typical 1–2% v/v inoculum |
| S0 | 10 | g/L | gives ~3 h batch phase before feed-on-demand activates |

## 5. Inputs (operator-controlled)

Seven decision variables — up from five — because DO and T are now real states with their own actuators (N, F_air for kLa; N also drives Q_ag).

| Input | Range | Decision var? | Reason |
|---|---|---|---|
| T_set_prod (°C) | 25–34 | ✓ | trades growth vs folding; biggest single-variable lever |
| IPTG (mM) | 0.01–0.5 | ✓ | now has a real interior optimum (Hill saturation vs k_tox_IPTG) |
| X_ind (g/L) | 5–30 | ✓ | sets pre-vs-post-induction biomass split |
| F_s (L/h) | 0.02–0.20 | ✓ | controls overflow vs starvation; range raised vs the 1 L baseline |
| S_f (g/L feed) | 200–600 | ✓ | bounded by glucose solubility (~700) |
| N_rpm | 200–1000 | ✓ | drives kLa AND Q_ag — coupled DO/T trade-off |
| F_air_vvm | 0.5–2.0 | ✓ | drives kLa; upper bound foaming-limited |
| pH | 7.0 (controlled) | ✗ | universal optimum; no operator value re-tuning |
| T_growth (pre-induction) | 37 °C | ✗ | universal pre-induction optimum |
| V0, S0 | 2.5 L, 10 g/L | ✗ | initial-charge convention. Was 1 L / 5 g/L — bumped for 5 L vessel |

## 6. What was removed / parked

- **`Kdo` collapsed to a constant `0.97`** — gone. DO is a state now and the Monod term is real.
- **Constant `phi = 0.5` post-induction** — gone. φ now Hill-saturated by IPTG dose so it tracks induction strength.
- **Single product `P`** — replaced by `(P_sol, P_ib)` so the optimiser sees the folding trade-off explicitly.
- **OUR / kLa / energy balance equations** — were parked in §6 of the previous version "for future scale-up extension"; now live in the model.
- **Arrhenius denaturation in `kd`** — still parked. Negligible within the operating envelope; re-add if exploring T > 40 °C scenarios.
- **Non-producer subpopulation (segregated model)** — not added. Relevant for long fed-batch (>48 h) where mutant escape dominates; over a 24 h horizon the homogeneous model is adequate.
- **Trace-metal limitation, formate accumulation** — not added. Strain-specific and depend on media recipe; not justified without a specific medium specification.
