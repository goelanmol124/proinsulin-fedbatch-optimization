# Choosing optimisation variables — rationale

This document explains **which inputs of the fed-batch process we optimise, which we set from biochemistry, and why**, without running any wet-lab experiment. It is the upstream justification for the 5–7 decision variables that appear in `OPTIMIZATION.md`. The supervisor's question — *"why this variable and not that one?"* — should be answerable by pointing at one row of one of the tables below.

The framing is: a bioprocess has ~30+ tweakable knobs. We cannot optimise all of them — both because the optimisation cost explodes, and because most of those knobs are biochemically *flat* near their conventional values. The job is to identify the subset where the model has real curvature, the controller has real authority, and the wet-lab cost of being wrong is real.

---

## 1. A decision rubric we can apply pairwise

For every candidate variable we score five things. All five are derivable from the model + textbook biochemistry; none requires an experiment.

| Symbol | Criterion | How we estimate it |
|---|---|---|
| **S** | Local sensitivity ∂(titre)/∂(var) at nominal point | Numerical derivative on the ODE; literature dose-response when the model lumps it away |
| **R** | Width of the *biologically feasible* window | Strain physiology — e.g. BL21 viable 20–42 °C, useful 25–37 °C |
| **C** | Hardware controllability (precise, cheap, low lag) | Engineering — RPM yes, intracellular ATP no |
| **∂$** | Cost gradient (how fast does "more" turn into "more money"?) | Reagent/utility prices; linear vs threshold |
| **K** | Coupling with other decision variables | Shared terms in the µ, qp, qs equations |

**Inclusion rule.** Include a variable if `S × R` is large *and* either `C` is high or `K` with another included variable is low. Exclude if `∂$` dominates `S × R`, or if the model literally cannot resolve the variable (sensitivity through the model is zero by construction).

**Pairwise comparison.** Given two candidates A and B competing for one slot, prefer the one with higher `S × R`; break ties with `C`, then by lower `K`. The whole document below is essentially this rule applied 25 times.

A separate consideration on top of the rubric — **shape of the titre vs variable curve.** This dictates the *optimiser* (gradient vs global vs grid), not whether to *include* the variable, but it informs the inclusion cost.

| Shape | Examples | Optimiser implication |
|---|---|---|
| Unimodal smooth, broad | µ_set, T_prod | Cheap. Any global method finds it in <50 evals. |
| Unimodal sharp (Hill / interior max) | IPTG, X_ind | Needs reasonable resolution, but no multi-start. |
| Plateau with cliffs | DO, pH | Don't optimise — set above the cliff with margin. Constraint, not variable. |
| Monotone with diminishing returns | trace metals | Don't optimise — pick the saturating value once. |
| Step / discrete | carbon source, feed mode | Don't optimise as continuous — choose by physics, then optimise inside the chosen branch. |
| Multi-modal | rare here | Would need DE/CMA-ES; not encountered in this model. |

---

## 2. Exhaustive variable list — what could we touch?

Grouped by phase. The verdict column is one of: **OPT** (decision variable), **SET** (fix from biochemistry / literature, justified below), **CONSTRAINT** (handled by a controller, not the optimiser), or **STRUCTURAL** (architectural choice, not a knob).

### Media design (pre-sterilisation)

| Variable | Verdict | One-line reason |
|---|---|---|
| Carbon source identity (glucose / glycerol / mixed) | STRUCTURAL | Discrete; optimisation can't compare apples to oranges in one continuous space. Glucose chosen for BL21 + T7 baseline. |
| Initial glucose **S₀** | SET (defensible default 5–10 g/L) | Real biochemical effect but with a wide flat window; see §4 below. |
| Complex carbon (yeast extract, tryptone) | SET | Not represented in model (Yx/s is lumped). Optimisation over invisible parameters = optimising noise. |
| Nitrogen source identity (NH₄⁺ / urea) | SET | Discrete; NH₄⁺ standard. |
| **C:N ratio** | SET | See §6 — could be modelled but not worth the parameter-identifiability cost on this model. |
| Phosphate / sulfate / Mg | SET | Korz/Riesenberg recipe; saturating dose, monotone-with-diminishing-returns shape. |
| Trace metals & vitamins | SET | Standard cocktail; model has no separate state for any of them. |
| Antifoam prophylactic dose | CONSTRAINT | Reactive dose-on-demand. Excess kills kLa, so it's a side-constraint, not a knob. |
| Buffer capacity / initial pH | SET (7.0) | One flat optimum, set once. |

### Inoculum

| Variable | Verdict | One-line reason |
|---|---|---|
| Inoculum OD at transfer | SET (mid-exponential, OD ~3) | Lag duration ≪ batch length; effect on integrated titre <5%. |
| Inoculum volume fraction | SET (5–10% v/v) | Same. |
| Seed media match to production | SET (matched) | Mismatch causes lag — handled by SOP, not optimised. |
| Number of seed stages | STRUCTURAL | Discrete; chosen by scale ratio. |

### Initial reactor state

| Variable | Verdict | One-line reason |
|---|---|---|
| Initial fill volume **V₀** | SET (~60% of working vol) | Hardware constraint (must leave headroom for feed); no metabolic effect. |
| Initial T / pH at inoculation | SET (37 °C / 7.0) | Same as control setpoints. |

### Growth-phase control

| Variable | Verdict | One-line reason |
|---|---|---|
| Growth-phase temperature **T_growth** | SET (37 °C) | Flat near peak of µ_max(T); ~0.5 °C off costs <2% biomass. |
| **Induction temperature T_prod** | **OPT** | Decouples µ and qp through a two-optimum system — the canonical reason to optimise T at all. |
| Temperature ramp rate at induction | SET (fast — minutes) | Step is fine; cells don't care about the slope, only the destination. |
| pH setpoint | SET (7.0) | One flat optimum. |
| DO setpoint | CONSTRAINT (≥30%) | Plateau-with-cliff. Set above the cliff with margin, don't optimise. |
| RPM, airflow, O₂ enrichment | CONSTRAINT (cascade controller) | These are how the controller *achieves* the DO setpoint; not independent knobs. |
| Back-pressure | SET / hardware-limited | Marginal kLa gain; vessel-rated. |

### Feed strategy

| Variable | Verdict | One-line reason |
|---|---|---|
| Feed-onset trigger | SET (DO spike) | Physics chooses it: the spike *is* the signal that batch glucose ran out. |
| Feed mode (constant / linear / exp / DO-stat) | STRUCTURAL (exponential) | Exponential is the only mode that holds qs below qs_crit through the entire feed phase. Chosen on first principles, then we optimise *inside* this branch. |
| **Target growth rate µ_set** | **OPT** | The continuous knob inside the exponential branch. Controls acetate overflow over the longest phase of the run. |
| Feed glucose concentration `S_f` | OPT (secondary) | Couples to feed rate to set total carbon flux; can be jointly optimised with µ_set. |
| Supplemental N/Mg feeds | SET | Tied to the main feed at a fixed stoichiometric ratio. |

### Induction & production

| Variable | Verdict | One-line reason |
|---|---|---|
| **IPTG concentration** | **OPT** | Hill curve with interior maximum (toxicity vs Hill saturation). Strongest single qp lever for proinsulin. |
| Inducer identity (IPTG vs lactose) | STRUCTURAL | Discrete. IPTG default for tight T7 control; lactose only if cost-driven. |
| **Biomass at induction X_ind** | **OPT** | Sets the pre/post phase split; integral ∫qp·Xv dt has interior optimum in X_ind. |
| Induction profile (single step vs interspersed) | STRUCTURAL (single step) | IPTG isn't metabolised → cells can't clear it → pulse-induction not operationally feasible. See §5. |
| Post-induction µ (independent of growth µ?) | OPT (optional 6th variable) | The cell's metabolic state changes at induction; a separate µ_set_prod ≠ µ_set_growth can lift titre 10–20%. Worth including if compute budget allows. |
| **Batch time t_f** | OPT (if objective is productivity/profit) | Productivity = titre / time; clear interior optimum once cost of time is in the objective. |

### Hardware / scale (10 L vessel)

| Variable | Verdict | One-line reason |
|---|---|---|
| Vessel working volume | STRUCTURAL | Fixed for the case study. |
| Impeller geometry / sparger | STRUCTURAL | Fixed by the vessel. Not a process knob at this scale. |

That's the catalogue. The rest of the document defends each non-trivial decision.

---

## 3. Why so few variables survive — the orthogonality story

If the rubric in §1 looks aggressive, it's because **most bioprocess "knobs" fall into one of three categories that don't deserve to be in the decision space**:

1. **Plateau-with-cliff variables** (DO, pH, trace metals). One value works for a wide range; another value catastrophically fails. There's no continuous optimum to find — there's a *minimum acceptable value plus margin*. These belong in the constraint set, not the decision vector. Putting them in the optimiser wastes evaluations exploring an essentially flat surface and lets the optimiser accidentally drift onto the cliff.

2. **Model-invisible variables** (C:N ratio, trace metals, vitamins, antifoam). Our model has no state for nitrogen, no state for sulfur, no state for individual metals. The model's sensitivity to these is identically zero. Optimising them through this model is meaningless. They are real biological effects that should be set once via separate DoE on the wet-lab strain — not coupled into the operating-policy optimisation.

3. **Variables whose effect is upstream of a stronger variable.** Inoculum OD affects lag, which affects time-to-feed, which affects total biomass at induction — but we *directly* optimise X_ind. Adding inoculum OD on top is redundant; X_ind already captures the integrated effect. This is the same reason we drop S₀ in favour of feed-trigger biomass.

The five surviving variables are mutually low-coupling on the model — they each move a *different* physical quantity (folding bias, overflow margin, phase split, induction strength, productivity-vs-time trade-off) that no other variable can capture.

---

## 4. The S₀ and batch-phase question — answered properly

The user raised this directly. It deserves its own section.

### Does S₀ actually control acetate?

**Yes** — but it's threshold-like, not continuous. The mechanism:

For S ≫ Ks (Ks ≈ 0.05 g/L for glucose), the cell takes glucose at qs_max independent of S. qs_max in BL21 on glucose is ~1.0–1.3 g/g/h, which is *above* the overflow threshold qs_crit ≈ 0.3 g/g/h. So during the batch phase, **every cell is in overflow regardless of S₀ value**, as long as S > ~1 g/L. The cell doesn't see "more substrate", it sees "above the threshold" or "below".

Where S₀ does matter: **total time spent in overflow** = (S₀ − S_residual) / (X·qs_max). Higher S₀ → longer batch → more biomass accumulated at overflow rate → more cumulative acetate.

Approximate scaling: cumulative acetate from batch ≈ X_batch_end · (qs_max − qs_crit) · Y_ac · t_batch ≈ const · ln(X_batch / X₀) · Y_ac. The variation across S₀ ∈ [3, 20] g/L is roughly 0.5 to 2.5 g/L acetate at end-of-batch, **straddling the BL21 inhibition threshold of ~2 g/L**.

So S₀ matters — but it has a wide flat zone (5–10 g/L) where the answer is "fine" and sharp cliffs outside it. **This is exactly a SET variable, not an OPT variable**: pick the centre of the flat zone, gain reliability from being far from both cliffs.

### Why have a batch phase at all? Can't we go fed-batch from inoculation?

In principle, yes. In practice, three reasons not to:

1. **Feeding lag-phase cells is open-loop control with no observability.** During lag, you don't know µ — the population is acclimating, and qs is rising from near-zero to qs_max. Any exponential feed you commit to will either undershoot (substrate accumulates → overflow once cells wake up) or overshoot (carbon-limit cells before they've gained mass). The DO-spike heuristic that triggers feed onset *requires* the cell to have already reached qs_max in batch.

2. **A batch phase is a free identifier of µ_max-in-this-medium.** OD vs t during batch fits a clean exponential whose slope IS µ_max for this strain + media + lot. That slope tunes the exponential feed. Without a batch phase, you'd have to either trust literature µ_max (lot-dependent, unreliable) or run a calibration shake-flask alongside.

3. **Operational robustness.** Batch → DO-spike → feed is one of the most robust transitions in bioprocessing because the trigger is unambiguous. Substrate-sensor-driven feed from inoculation is fragile (sensors drift, foul, lag).

The cost we pay for these three benefits is the modest acetate accumulation in batch, which we manage by keeping S₀ in the 5–10 g/L window. **The architecture is "batch then fed-batch" not because S₀ has high titre leverage, but because the architecture is more reliable and more observable than the alternative.**

### What does acetate actually block?

This matters because the "keep acetate low" rule is otherwise just hand-wave. Mechanism:

- Acetate exists in equilibrium between **acetate⁻** (charged, can't cross membrane) and **acetic acid** (uncharged, lipophilic, crosses freely).
- Outside the cell at pH 7, the equilibrium sits far towards acetate⁻ (pKa = 4.76), so the membrane-permeable fraction is small but non-zero.
- Acetic acid that crosses the membrane finds the cytoplasm near pH 7.6, dissociates, dumps a **proton inside the cell**.
- The cell has to pump that proton out against the gradient. ATP cost per proton ≈ 1/3 ATP. At 2 g/L extracellular acetate, the steady-state intracellular proton load consumes ~10–15% of cellular ATP budget — pure overhead.
- Beyond ATP drain, intracellular acetate-CoA pools elevate, which: (a) inhibits methionine biosynthesis (Roe et al. — Met is rate-limiting for many recombinant proteins), (b) triggers the stringent response above ~3 g/L, down-regulating ribosomes, **directly killing qp**.

So acetate inhibition is **stronger on qp than on µ**. For a producer strain, the acetate ceiling is set by *product* loss, not biomass loss. This is why holding qs below qs_crit in the feed phase via µ_set is the highest-leverage continuous knob in the whole process.

---

## 5. Growth-then-induction vs interspersed phase architecture

The user asked whether we should interleave growth and induction. Short answer: **no, single-step induction is correct for this strain + inducer**. Long answer:

For interspersing to pay off, two conditions must hold:
- qp must **decline meaningfully** during sustained induction (otherwise why pulse?).
- The cell must be able to **clear the inducer** between pulses (otherwise pulsing is just delayed continuous).

For BL21 + IPTG + proinsulin:
- qp does decline post-induction (proteolysis + ribosome diversion + IB overload), but on a 10–20 h timescale. For a 30 h batch with ~8–12 h production phase, the decline within one phase is tolerable — pulse cycles can't fit.
- **IPTG isn't metabolised by E. coli.** It's a gratuitous (non-metabolisable) lactose analogue — that's the point of IPTG. Cells cannot clear it. The only way to "pulse off" is to dilute, which means feeding fresh medium fast enough to wash out IPTG, which dilutes biomass and product too — net negative.

If the inducer were lactose (metabolised, cleared by the cell), pulse-induction would be on the table. With IPTG, it isn't. This is a **structural** choice driven by the inducer's mechanism, not an optimisation outcome.

---

## 6. C:N ratio — model it, or skip it?

The user is right that we *could* extend the model to represent N as a separate state, with its own yield Yx/n and uptake kinetics. The question is whether it's worth the parameter-identifiability cost.

**Arguments for adding N to the model:**
- C:N ratio has a real, well-documented ~10–20% titre effect in literature (Korz et al. 1995, Riesenberg & Guthke 1999).
- N-limitation has a distinctive metabolic signature (carbon channelled to TCA / acetate rather than biomass) that the current lumped Yx/s can't predict.
- It opens the door to optimising N feed separately, which can lift titre on cells with high amino-acid-rich products.

**Arguments against:**
- Adds at least 3 new parameters (Yx/n, Kn, qn_max) that are weakly identifiable from typical online measurements (OD + HPLC glucose + HPLC product). Without offline NH₄⁺ assays, the fit is under-determined.
- The 10–20% titre effect is **smaller** than the leverage of IPTG (3–5×), T_prod (~2×), µ_set (~2×), or X_ind (~2×). Adding it to the decision vector means trading optimiser resolution on the big knobs for resolution on a smaller one.
- The mechanistically-correct way to set C:N for a recombinant strain is **one-time DoE on the wet-lab strain at a single operating point**, not coupled into operating-policy optimisation that has to be re-run every time T_prod or IPTG moves.

**Verdict: skip in this model.** If the supervisor asks "how do we know N isn't limiting?" the defensible answer is: (a) the Korz medium puts N in stoichiometric excess for the expected biomass, (b) if we see Yx/s drop systematically across runs in residuals against the lumped model, that's the signal to extend the model, not before.

This is a textbook case of *not letting model complexity outrun data*.

---

## 7. How to set the variables we are NOT optimising — defensible defaults

A common failure mode is to defend an optimisation by saying "we focused on these 5 variables" without specifying what we did with the other 25. Here's the full set of SET variables with how we'd justify each value to the supervisor, in the order an operator would set up a run:

| Variable | Default | Source / reasoning |
|---|---|---|
| Carbon source | Glucose | BL21 + T7 standard; glycerol only if acetate problem persists after feed-phase optimisation |
| **S₀** | 8 g/L | Centre of the 5–10 g/L flat zone; ~4 generations of batch growth, ~1.5 g/L cumulative acetate (below 2 g/L threshold with margin) |
| Complex N (yeast extract) | 5 g/L | Korz recipe; covers vitamin / co-factor gaps not in defined media |
| Inorganic N (NH₄Cl) | 2.5 g/L | Stoichiometric excess for expected 40 g/L DCW final biomass at Yx/n ≈ 8 |
| Phosphate (KH₂PO₄ / K₂HPO₄) | Korz | Buffering + macronutrient; saturating |
| Mg, Fe, trace metals | Riesenberg cocktail | Standard, saturating |
| Vitamins (thiamine) | 4.5 mg/L | BL21 thiamine auxotrophy; standard |
| Antifoam | Silicone, reactive dose | Silicone over PPG (lower kLa penalty per dose); reactive over scheduled |
| pH setpoint | 7.0 | BL21 cytoplasmic optimum; flat ±0.2 |
| DO setpoint | 30% | Above the 20% overflow cliff with 10% margin for cascade lag |
| Growth T | 37 °C | µ_max peak; flat ±2 °C |
| RPM range (cascade) | 300–1000 | Vessel limit; cascade controller chooses |
| Airflow range (cascade) | 0.5–2 vvm | Vessel limit; cascade controller chooses |
| O₂ enrichment | Off by default, enable above 1.5 vvm | Cost; only if mechanical+air can't hold DO |
| V₀ | 6 L (60% of 10 L) | Leaves headroom for ~3 L feed |
| Inoculum OD | OD₆₀₀ ≈ 3 (mid-exp) | Lag minimisation without acetate carry-over |
| Inoculum % | 5% v/v | Standard; below 1% extends lag, above 10% adds seed-train cost |
| Feed-onset trigger | DO spike (>20% rise in <2 min) | Physics |
| Feed mode | Exponential | Only mode that holds qs<qs_crit through the feed phase |
| Inducer | IPTG | Tight T7 control; lactose only on cost grounds |
| Temperature ramp at induction | Step (cascade-limited) | Cells don't care about slope |
| Induction profile | Single step | IPTG isn't clearable |

The supervisor's "why this value?" question for any of these is answered by the source/reasoning column.

---

## 8. Hard constraints that come from solving equations, not optimising

Several variables get fixed not by literature but by **directly solving the model for the value that keeps a constraint satisfied with a safety margin.** This is much more reliable than including them in the optimisation:

**DO setpoint.** The OUR (oxygen uptake rate) is OUR = qO₂·Xv ≈ (qs · Y_O2/S) · Xv. The OTR (oxygen transfer rate) is kLa·(C* − C_DO). Setting OUR = OTR and solving for C_DO given the peak expected Xv·qs gives the minimum DO that the cascade can actually hold. We pick the setpoint **above** that, with margin. At T = 37 °C, C* ≈ 6.7 mg/L; the cascade can hold ≥30% (~2 mg/L) for Xv up to ~40 g/L at full RPM + airflow.

**Stirrer / airflow cascade.** Given a DO setpoint, we *solve* kLa = OUR/(C*−C_DO) and back out the required RPM from the kLa correlation kLa = c·N^a·(F_air)^b. This is algebra, not optimisation.

**Feed flow at switch.** At the moment of feed onset, F_0 = (µ_set / Yx/s) · X(t_switch) · V(t_switch) / S_f. Closed-form, not optimised.

**Maximum acetate.** µ_set is bounded above by the µ at which qs hits qs_crit: µ_set < µ_max · qs_crit / qs_max. With qs_max ≈ 1.2 g/g/h, qs_crit ≈ 0.3 g/g/h, the cap is ~0.25·µ_max ≈ 0.18 h⁻¹. The optimiser searches *inside* the (0.08, 0.18) box; the box itself is set algebraically with margin.

For reliability we apply margins **on both sides** of each constraint:
- DO ≥ 30%, not the cliff at 20% — 10% buffer.
- µ_set capped at 0.85·µ_overflow — 15% buffer for parameter uncertainty.
- T_prod ≥ 25 °C — below this folding becomes kinetically limited and growth stalls.
- IPTG ≤ 1 mM — well above the Hill saturation; further raises only toxicity.

The optimiser sees a box-constrained problem; the constraints' values come from explicit equation-solving, not guesswork.

---

## 9. Layered optimisation — coarse-to-fine in practice

The user is right that optimisation can and should be layered. The structure:

**Layer 0 (offline, once per strain):** structural choices. Carbon source, inducer, feed mode, induction architecture. Discrete; chosen by physics + strain knowledge. Done in this document.

**Layer 1 (offline, once per medium):** SET variables — media composition, inoculum SOP, V₀, pH, DO, T_growth. One-time DoE in shake flasks if data exists, literature defaults otherwise. Done in §7.

**Layer 2 (this study — operating policy):** the 5–7 OPT variables. Run DE globally on the model. ~600 evaluations. Identifies a basin.

**Layer 3 (refinement):** local polish around the DE basin (L-BFGS with finite-difference gradients). Tightens the optimum by the last few %.

**Layer 4 (wet-lab — out of scope here):** confirmatory runs at the model optimum and at ±1 standard-deviation perturbations to test reliability. Update the model if residuals are systematic. Re-run layer 2.

The reason for this layering is that **layers 0–1 change the model structure or pin its inputs**, so they invalidate any layer-2 optimisation that's already run. Doing them first means layer 2 only runs once. Doing them in the optimiser (mixed continuous + discrete, full media composition, all knobs) would mean ~30 dimensions and combinatorial discrete branches — intractable, even if cheap simulator calls suggest otherwise.

Within layer 2, we could in principle layer further: optimise the most-leverage variables coarsely first (IPTG + T_prod on a grid), then fix those and optimise µ_set + X_ind + t_f. Empirically this is rarely needed for 5 dims with DE — the global method handles it. We'd only resort to it if DE convergence stalled.

---

## 10. Comparing the final shortlist — pairwise, again

This is the deliverable: a table the supervisor can read top-to-bottom and ask "why this and not the alternative right next to it?"

| Rank | Variable | Why included | Closest excluded alternative | Why we picked the included one |
|---|---|---|---|---|
| 1 | **IPTG concentration** | Hill curve with interior max in qp; ~3–5× titre swing; model resolves it explicitly | Trace metals | Model can't resolve trace metals (no state). Optimising what the model doesn't see is optimising noise. |
| 2 | **Induction T_prod** | Two-optimum system (µ_max peak ≠ qp peak); 2× titre swing | Growth-phase T_growth | T_growth is on a flat plateau near 37 °C; T_prod sits on a real gradient. |
| 3 | **µ_set during feed** | Single biggest acetate-overflow knob; acts over the longest phase (~20 h) | S₀ | µ_set acts over 20 h, S₀ over 6 h; same biochemical mechanism (overflow), 3× longer leverage arm. |
| 4 | **X_ind (biomass at induction)** | Sets pre/post phase split; integral ∫qp·Xv has interior max | Inoculum OD | Inoculum OD shifts the time-axis by ~hours; X_ind shifts the phase boundary itself. ~10× the leverage. |
| 5 | **Batch time t_f** | Productivity / cost lever; broad optimum; cheap dim | Hard-coded stopping rule at fixed t | Explicit optimisation lets the objective (titre vs productivity vs profit) flow through to the stopping decision. |
| 6 (opt.) | **µ_set_prod** (post-induction µ, separate from growth) | Cell physiology changes at induction; one µ for both phases is a compromise | Single µ_set | Adds 10–20% titre at the cost of 1 dim. Include if compute budget ≥ 1 min. |
| 7 (opt.) | **Feed glucose S_f** | Couples to feed rate; jointly with µ_set sets total carbon flux | Fix S_f = 500 g/L | Including it lets the optimiser balance dilution vs osmotic shock at high biomass; usually finds 300–500 g/L. |

The 5-variable core is the defensible minimum. The two optional adds (µ_set_prod and S_f) are what we promote to a 7-dim problem if the compute budget allows and the supervisor wants tighter productivity numbers. Anything past 7 dims is diminishing returns: the next-best candidate (e.g. post-induction T ramp rate, or feed-step at IPTG addition) brings <5% titre at the cost of an extra optimiser dimension and worse identifiability.

---

## 11. Three-axis justification of the final shortlist

The shortlist is defended along three axes: **(a) biochemical implications**, **(b) methodological limitations**, **(c) resource constraints (time, cost, compute)**. The rubric in §1 collapses these into one score; here we expand them back out so the reader sees each axis named.

### IPTG concentration

- **(a) Biochemical implications.** IPTG is a gratuitous inducer of the lac operator; binding lacI relieves T7 repression. The dose-response of *expression* is Hill-saturating, but the dose-response of *soluble titre* is non-monotone because (i) high IPTG → high mRNA → ribosome overload → misfolding into inclusion bodies, (ii) IPTG above ~1 mM is metabolically toxic (induces stress response, slows µ), and (iii) the qp peak in this model sits at 0.3–0.6 mM. This shape is precisely the case where optimisation pays off — neither "max" nor "min" is correct.
- **(b) Methodological limitations.** Strong, model-resolvable derivative (IPTG enters qp via an explicit Hill term and µ via a toxicity term in `MODEL.md`). One-shot dose at induction — a scalar, not a profile — fits cleanly into a box-constrained optimiser. No state-estimation problem: IPTG is set by the operator, not inferred.
- **(c) Resource constraints.** Reagent cost trivial (~$1/g, mg-scale dose, ~$0.01 per run). Stable to autoclave. No additional sensor or hardware. Compute cost: one dimension in the DE search, ~50 extra evals to localise.

### Induction temperature T_prod

- **(a) Biochemical implications.** µ_max(T) and qp_soluble(T) have **different optima**: cardinal-temperature µ peak near 37 °C (BL21), but folding kinetics + inclusion-body partitioning push the soluble-titre peak down to 25–30 °C. This decoupling is the unique biological justification for *changing* T during a run, and the only reason T appears in the decision space at all. Lower T also slows proteolysis (k_deg is Arrhenius), which compounds the soluble-fraction gain.
- **(b) Methodological limitations.** Cleanly modelled in `MODEL.md` as a Gaussian production envelope around T_opt_p with σ_p = 4 °C. Identifiable from titre-vs-T scan data, no hidden state. The step at induction is fast (minutes via jacket); we don't have to optimise a ramp profile, just a destination scalar.
- **(c) Resource constraints.** Free to vary — the jacket and chiller already exist. No extra reagents. Compute: one dim. The only real cost is run-time: lower T_prod slows growth and product accumulation, which lengthens batch time and ties up the reactor. This trade-off is captured by including t_f in the same optimisation.

### Target growth rate µ_set during fed-batch

- **(a) Biochemical implications.** The single biggest **acetate-overflow** lever. µ_set determines qs (via Yx/s); when qs > qs_crit the cell short-circuits TCA and excretes acetate, which then poisons qp (see §4). µ_set acts over the longest phase of the run (~20 h), so its integrated influence on titre dominates anything that acts on the batch phase. There is a real interior optimum: too low wastes time (under-fed cells, low biomass), too high hits overflow.
- **(b) Methodological limitations.** Implemented as the exponent in an exponential feed F(t) = F₀·exp(µ_set·(t−t_switch)) — a clean scalar parameter, not a profile. Model resolves the qs → acetate → qp pathway explicitly, so sensitivity is well-posed. Open-loop µ-stat (set-and-forget exponential) is simpler than feedback µ-control, doesn't need an online biomass sensor.
- **(c) Resource constraints.** No extra hardware (one peristaltic pump). No reagent change. Compute: one dim. Cheap.

### Biomass at induction X_ind

- **(a) Biochemical implications.** Defines the pre/post boundary of the run. The objective ∫qp·Xv dt over the production phase has an interior optimum because: induce-too-early → long production window but low Xv; induce-too-late → high Xv but production window truncated by acetate accumulation, oxygen limit, proteolysis, and run-time cap. Couples weakly to IPTG (high IPTG → growth arrests sooner → optimal X_ind shifts up).
- **(b) Methodological limitations.** Implemented as a *trigger* (induce when Xv reaches X_ind), not a clock time — so it's robust to run-to-run variation in lag and µ_max. Requires an online biomass measurement (OD or capacitance probe); both are standard at 10 L.
- **(c) Resource constraints.** No extra reagent. Standard OD probe (~$2k capex, already on vessel). Compute: one dim. The main real-world cost is the operator overhead of switching to induction phase at a non-fixed clock time — manageable in supervised runs, trivially automated in PLC.

### Batch time t_f

- **(a) Biochemical implications.** After qp drops below the rate of proteolytic loss + maintenance carbon consumption, the cell is *net-degrading* product. There is a clear interior optimum once any cost-of-time enters the objective (utilities, depreciation, labour, downstream window). Without t_f in the objective, the optimiser will run the reactor forever for the last 1% of titre.
- **(b) Methodological limitations.** Trivial — it's the simulation endpoint. The only methodological subtlety is whether to optimise titre or productivity: optimising titre alone implicitly picks t_f at the asymptote, which is unrealistic. The objective `g(t_f) = P(t_f)/t_f` (productivity) or `P(t_f) − c·t_f` (profit) handles this directly.
- **(c) Resource constraints.** No reagent or hardware impact. Compute: one dim. **This variable is unique in that including it *reduces* total resource use** — the optimiser will pick a shorter run if the marginal titre doesn't justify the marginal time.

### Joint justification — why exactly these 5 (not 4, not 7) for the demo

- **(a) Biochemical.** They span the four biologically distinct levers: folding bias (T_prod), translation rate (IPTG), carbon flux (µ_set), phase split (X_ind), and stopping rule (t_f). No two of them are redundant; each addresses a mechanism no other can substitute for. A reasonable 4-variable cut (drop t_f, fix at a reasonable wall-clock) is acceptable if the objective is pure titre, but loses productivity / cost handling. A 7-variable extension (add µ_set_prod, S_f) is biochemically justified but methodologically marginal at this model fidelity.
- **(b) Methodological.** All 5 are scalars with clean model derivatives, box-constrained, no profiles or feedback policies needed. This keeps the optimiser inside the regime where DE + L-BFGS polish is provably fine. The next candidate variables (post-induction T ramp, IPTG re-dose, time-varying feed profile) would push us into trajectory optimisation (collocation / RL), which is methodological overkill for the demonstration's claim ("model-based beats trial-and-error").
- **(c) Resource constraints.** DE on 5 dims with popsize=15, maxiter=40 ≈ 600 simulator calls; at ~0.2 s per fed-batch ODE solve ≈ 2 minutes wall-clock on a laptop. Each added dimension scales linearly in DE evals at fixed convergence quality, so 7 dims ≈ 3 minutes — still trivial. *Beyond ~10 dims*, DE convergence quality degrades faster than linearly and we'd want sequential / surrogate methods. We're well below that. Wet-lab cost: zero (model-based). Hardware cost: zero (all 5 are knobs on standard equipment). The only meaningful resource cost is the **wet-lab confirmation** of the model optimum afterwards, which is one fed-batch run (~$2–5k, ~3 days) — orders of magnitude less than the ~30 trial-and-error runs the old SOP would have used to find the same point.

The framing the supervisor should walk away with: *we picked the smallest variable set that (a) spans every distinct biological mechanism in our model, (b) keeps the optimisation inside a regime where the chosen algorithm has theoretical guarantees, and (c) finishes in minutes on commodity compute while telling us where to commit ~$5k of wet-lab capacity. Adding more variables doesn't improve any of these axes; dropping any of these five loses a biological mechanism the model can resolve.*

---

## 12. What would change the conclusion?

For completeness — the cases that would justify adding variables we currently exclude:

- **DO becomes a real variable**, not a constraint, if we move to a vessel where kLa is hard-limited (large scale, viscous broth). At 10 L this doesn't happen; at 1000 L it does. Then DO joins the decision set and the cascade controller can't paper over it.
- **C:N enters the decision set** if we move to a high-density (>80 g DCW/L) process where N stoichiometric excess is no longer cheap to maintain.
- **pH enters** if we switch to a strain or product where the production-phase pH optimum diverges from the growth-phase optimum (some yeast / mammalian systems). Not relevant to BL21.
- **Antifoam** enters if foam-out probability becomes a binding constraint at the optimum. Solved as a side-controller for now.
- **Inoculum OD / phase** enter for short batches (<10 h) where lag is a meaningful fraction. Not relevant here.

Each of these is a clear, scope-bounded trigger. The supervisor's follow-up question — "what would make you change the variable list?" — has an answer at this level of specificity for every excluded variable.

---

## TL;DR

We optimise **5 variables** (IPTG, T_prod, µ_set, X_ind, t_f), and 7 if we add µ_set_prod and S_f.

- **Biochemically (§11a):** each of the 5 addresses a mechanism no other can substitute for — folding bias, translation rate, carbon flux, phase split, stopping rule.
- **Methodologically (§11b):** all 5 are scalars with clean model derivatives and box constraints — DE + L-BFGS polish is sound. The next-tier candidates would push into trajectory optimisation, which exceeds the demo's claim.
- **Resource-wise (§11c):** 5 dims ≈ 2 min wall-clock; 7 ≈ 3 min. All knobs are on standard equipment. Wet-lab cost displaced: ~30 trial-and-error runs (~$60–150k) collapses to ~1–3 confirmatory runs (~$5–15k).

Everything else is either: (i) on a plateau-with-cliff where a constraint + margin is the correct answer, (ii) invisible to this model, (iii) upstream of an already-included variable, or (iv) a structural / discrete choice fixed by first principles. Each inclusion has a model-resolvable derivative; each exclusion has a written reason. The 5-variable surface is unimodal-or-broadly-multimodal with interior optima, which is why DE with L-BFGS polish is the chosen optimiser (see `OPTIMIZATION.md`).
