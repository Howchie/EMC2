# Plan: cross-subrace correlated-target capacity for LogicalRules

Date: 2026-07-17  
Related architecture: `bawl_corr_exact_kernel_plan.md` (implemented and
verified on `bawl-correlated-race`; the concrete reuse points below name the
code that actually exists, not the plan's original sketches)
Primary targets: `LogicalRulesLBA` with `OR_DETECTION_ANALYTIC`,
`OR_DETECTION_GNG`, and the OR/AND/XOR/ID choice rules

## Outcome

Add a redundant-target capacity model in which the latent drift draws for
target A and target B are correlated when both target stimuli are present,
while the two targets remain in separate logical-rule subraces.

The model supports:

- a mean capacity effect (`kappa`),
- a between-trial capacity-variability effect (`tau`),
- either effect alone or both together,
- detection rules with only A/B target detectors (`OR_DETECTION_ANALYTIC`),
- go/nogo detection (`OR_DETECTION_GNG`): the A/B capacity math is identical;
  the nogo accumulator is a zero-loading independent racer,
- choice rules with A/n_A and B/n_B subraces,
- independent n_A, n_B, and nogo racers,
- unchanged single-target trials.

The central likelihood operation is:

```text
draw/integrate one shared capacity factor z
    -> evaluate A subrace conditional on z
    -> evaluate B subrace conditional on z
    -> assemble the logical rule conditional on z
    -> integrate the complete trial likelihood over z
```

This is not the same as treating A and B as competitors in one race. The
shared-factor integration machinery from the completed BAwLCorr architecture
is reusable, but the single-race pair-density evaluator is not the logical-rule
integrand for choice models.

## Scope and non-goals

In scope:

1. Correlation between the target drifts A and B only.
2. A/B correlation only in the dual-target (`AB`) condition.
3. Independent residual drift variation for all racers.
4. Independent n_A and n_B racers in choice models, and the independent nogo
   racer in `OR_DETECTION_GNG`.
5. Correct positive-drift normalization for the active set, with the A/B
   joint term computed closed-form (bivariate normal orthant), not by
   quadrature.
6. Finite RTs, censoring, truncation, omissions, contaminants, and expansion
   through the normal logical-rule likelihood contracts.
7. A reference R simulator and a production C++ likelihood.
8. Single-target trials remain on the existing batched fast path untouched;
   only AB trials pay any factor-integration cost.

Out of scope:

- correlation between A and n_A, B and n_B, n_A and n_B, or any other pair;
- a general multi-factor covariance model;
- correlation on single-target trials;
- treating the cross-subrace pair as an ordinary within-race winner/loser
  pair;
- changing the public BAwLCorr model or its correct/error semantics.

## Mathematical specification

### 1. Trial activation

Let

```text
pair_active_j = 1(S_j == "AB")
```

For detection models (`OR_DETECTION_ANALYTIC`), the active target set is:

```text
S == "A"   -> A
S == "B"   -> B
S == "AB"  -> A and B
```

For go/nogo detection (`OR_DETECTION_GNG`), the nogo accumulator is active on
every trial in addition to the set above. It never carries a capacity
loading; conditional on `z` it enters the rule assembly exactly as in the
ordinary GNG likelihood.

For choice models, the active subraces are:

```text
S == "A"   -> A and n_A
S == "B"   -> B and n_B
S == "AB"  -> A/n_A and B/n_B
```

The fixed accumulator roles are retained on every expanded trial. Activity is
controlled by the condition/RACE mask; an inactive fixed role must not enter
the race likelihood, positivity normalizer, or logical-rule state.

### 2. Shared capacity factor

For an AB trial, define one trial-level factor:

\[
G_j = \kappa_j + \tau_j Z_j,
\qquad Z_j \sim N(0,1).
\]

For target base means \(\mu_A\) and \(\mu_B\):

\[
V_A = G_j\mu_A + \epsilon_A,
\qquad
V_B = G_j\mu_B + \epsilon_B,
\]

where

\[
\epsilon_A \sim N(0,sv_A^2),
\qquad
\epsilon_B \sim N(0,sv_B^2),
\]

and the residuals are independent conditional on \(Z_j\).

Thus the target drift covariance is

\[
\operatorname{Cov}(V_A,V_B)=\mu_A\mu_B\tau_j^2,
\]

with marginal variances

\[
\operatorname{Var}(V_i)=sv_i^2+\mu_i^2\tau_j^2.
\]

The nontarget racers remain independent:

\[
V_{n_A}=\mu_{n_A}+\epsilon_{n_A},
\qquad
V_{n_B}=\mu_{n_B}+\epsilon_{n_B}.
\]

Their shared-factor loadings are exactly zero.

For single-target trials, the capacity loading is zero. Equivalently, use
`kappa = 1`, `tau = 0` for those trials. This ensures that the single-target
RT distributions are unchanged by the redundant-target capacity model.

### 3. Parameterization choices

The implementation must distinguish two related but different uses of a mean
capacity effect:

1. A fixed design effect on the target means, for example
   `v ~ 0 + lR + lR:DT`.
2. The multiplicative capacity factor above, where `kappa` multiplies the
   target means in the AB cell.

They should not both be added accidentally. The supported model variants are:

```text
baseline:       kappa = 1, tau = 0
kappa only:     kappa free, tau = 0
tau only:       kappa = 1, tau free
both:           kappa free, tau free
```

If an additive pair effect is desired instead, it should be represented in the
ordinary `v` design and the correlated capacity factor should retain
`kappa = 1` unless a second multiplicative effect is explicitly intended.

Use a positive natural-scale parameterization for `kappa`, with default 1,
and a nonnegative parameterization for `tau`, with default 0. A zero-valued
`tau` must be a supported fixed/exception value so the ordinary logical-rule
path is recovered exactly.

### 4. Conditional logical likelihood

The fundamental identity is:

\[
L_j = \int L_j(\text{data}\mid z)\phi(z)\,dz.
\]

The logical-rule assembly must occur inside the integral. It is incorrect to
integrate the A and B subraces separately and then assemble their marginal
probabilities, because that removes their shared trial-level dependence.

#### Detection-only model

For an AB detection trial, conditional on `z`:

\[
L_{AB}(t\mid z)
 = f_A(t\mid z)S_B(t\mid z)
 + f_B(t\mid z)S_A(t\mid z).
\]

Therefore:

\[
L_{AB}(t)
 = \int\left[f_A(t\mid z)S_B(t\mid z)
              +f_B(t\mid z)S_A(t\mid z)\right]\phi(z)\,dz.
\]

For A-only and B-only trials, use the ordinary single-detector likelihood.

#### Choice/logical-rule model

For each quadrature node, evaluate the A and B subraces conditionally on the
same `z`. A subrace evaluator returns the state quantities required by the
logical rule, such as:

```text
G_yes(t | z)
G_no(t | z)
S_undecided(t | z)
```

for each active subrace. The existing logical-rule assembler then returns the
conditional event density or censoring mass:

\[
L_{\text{rule}}(t,R\mid z).
\]

The production likelihood integrates that complete conditional result over
the shared factor.

Conditional on `z`, A/n_A and B/n_B are independent subraces. Marginally,
their rule outcomes are dependent because A and B share `z`.

### 5. Positive-drift semantics

Positive-drift models must use joint active-set conditioning, not an
independent marginal `sv` approximation and not a `G > 0` shortcut.

For an AB trial, the joint target positivity mass is closed-form. Marginally,
\((V_A,V_B)\) is jointly bivariate normal with

\[
\mathrm{E}[V_i] = \kappa\mu_i,
\qquad
\mathrm{sd}_i = \sqrt{sv_i^2+\tau^2\mu_i^2},
\qquad
\rho_{AB} = \frac{\tau^2\mu_A\mu_B}{\mathrm{sd}_A\,\mathrm{sd}_B},
\]

so

\[
q_{AB} = P(V_A>0,\,V_B>0)
\]

is one bivariate-normal orthant probability. Compute it with the existing
`norm_cdf_2d()` / `bawl_corr_pair_positive_normalizer()` machinery from
`src/bawl_corr_exact.h` (including its nearly-singular \(\rho\) limit); do
not reimplement it and do not integrate \(\int q_A(z)q_B(z)\phi(z)\,dz\) by
quadrature in production. The quadrature identity is retained only as a test
cross-check of the closed form.

For an AB choice trial, independent n_A/n_B (and nogo) positivity terms are
univariate normal probabilities that multiply the conditional node
contribution; only A/B share the joint term above.

For A-only trials, only the active A subrace participates in the positivity
normalizer. B must not be included merely because its fixed role exists.

The final per-trial contract is:

```text
log numerator = log ∫ conditional_event_mass(z) phi(z) dz
log denominator = log(observation normalizer)
trial log likelihood = log numerator - log denominator
```

The observation normalizer factors as the closed-form `q_AB` (above) times
the independent univariate positivity terms; truncation windows additionally
require the z-integrated total survivor evaluated at both endpoints,
following the BAwLCorr `log S_all(t)` contract (every survivor factor inside
the endpoint subtraction; only time-independent `q` constants cancel).

Each normalizer is applied exactly once. Follow the convention the BAwLCorr
implementation actually adopted: component evaluators receive their
normalizer explicitly as an argument and divide exactly once at assembly.
(The `LogComponentValue` struct sketched in the BAwLCorr plan was not
implemented; do not introduce it here.)

## Architecture

### 1. Reuse from the completed BAwLCorr architecture

The BAwLCorr implementation exists and has been verified; reuse these
concrete pieces rather than the plan-era abstractions:

- `src/gh_quad.h`: `GHRule`/`gh_rule()` cached rules,
  `gh_standard_normal_weight()`, `AGHCenter`, `agh_center_from_scan()`,
  `agh_log_weight()`, and the environment node-count override. The override
  helper is currently named `bawl_corr_quad_nodes()`; rename/move it to a
  shared name rather than duplicating it.
- `src/bawl_geometry.h`: `BAwLTimeGeometry` and `BAwLPreparedRow`.
  `BAwLPreparedRow` already has exactly the conditional form this model
  needs — `v_q = v0 + slope*z` with residual `sv` — and its prepared
  endpoint modes (`f0`, `F0`, `q`, `q-F0`) are precisely the terms the
  detection integrand and its survivors require. `LogicalRulesLBA` is plain
  LBA (`k = 0`), which the geometry covers.
- `src/bawl_corr_exact.h`: the guarded univariate incomplete-normal moments,
  `norm_cdf_2d()` usage patterns, and
  `bawl_corr_pair_positive_normalizer()` for the closed-form `q_AB`.
- Shared-state lifecycle: build once outside the particle loop, direct
  `ParamTable` column pointers, per-particle layout classification, stable
  log-sum-exp node reduction, per-trial finalisation
  (contaminant/truncation/floor/expansion applied once), and env-gated route
  counters — all following the `BAwLCorrSharedState` /
  `bawl_corr_classify_particle()` pattern in `src/particle_ll.cpp`.

Do not reuse the BAwLCorr assumption that `lM` identifies the correlated pair.
For this model, pair membership is explicit: the correlated pair is A/B.

Do not reuse BAwLcorr's `rho -> loading` Ttransform or its loading
convention. BAwLcorr loads the factor on the sv scale
(`mu + sign(rho)*sv*loading*z`); the capacity factor loads on the drift-mean
scale (`slope = tau*v_base`). The two simulators will look similar — do not
copy the transform.

The final BAwLCorr exact pair kernel is for a pair of racers within a common
race and should remain unchanged. It is not the conditional logical-rule
integrand for the choice model. The reusable numerical abstraction is the
shared-factor integrator, the prepared-row conditional geometry, and the
explicit component/normalizer convention.

### 2. New logical-rule correlated state

Add a model-specific shared state, following the final BAwLCorr shared-state
lifecycle:

```text
LogicalRulesCapacitySharedState
```

It should contain data-fixed information:

- number of rows, roles, and unique trials;
- condition/rule/response codes;
- fixed role indices for A, B, n_A, and n_B;
- active-role masks for A-only, B-only, and AB trials;
- pair-active mask;
- censoring/truncation bounds;
- winner/response and expansion mappings;
- parameter-column indices for `v`, `sv`, `kappa`, and `tau`;
- optional parameter-column indices for deterministic capacity effects;
- direct parameter-column pointers;
- reusable per-trial quadrature and likelihood buffers.

The state must not infer the pair from `lM`. `lM` may remain available for
ordinary models, but it is not required for capacity correlation.

### 3. New trial layout

Use a logical-rule-specific layout, following the final
`BAwLCorrTrialLayout` design but distinguishing subraces:

```text
LogicalRulesCapacityTrialLayout
    rule_code
    condition_code
    pair_active
    row_A, row_B
    row_nA, row_nB, row_nogo
    active_A, active_B, active_nA, active_nB, active_nogo
    route
    positive_drift_route
    valid
```

Routes should at least distinguish:

```text
ordinary_single_or_tau_zero
capacity_detection_pair      (OR_DETECTION_ANALYTIC and OR_DETECTION_GNG;
                              GNG adds the independent nogo racer to the
                              same conditional detection integrand)
capacity_logical_pair
invalid
```

The layout must classify activity before positivity dimension or quadrature
work is assigned. In particular, a fixed inactive role contributes neither a
correlated loading nor a positivity dimension.

Validate that `kappa` and `tau` are shared across A and B within a trial when
they define one common capacity factor. If different target-specific values
are eventually needed, that should be a separate multi-loading model rather
than an accidental row-wise mismatch.

### 4. Conditional parameter view

At factor node `z`, construct a conditional view without changing the base
particle parameters:

```text
A:   v_q = kappa * v_base + tau * v_base * z,  sv_q = sv_base
B:   v_q = kappa * v_base + tau * v_base * z,  sv_q = sv_base
n_A: v_q = v_base,                              sv_q = sv_base
n_B: v_q = v_base,                              sv_q = sv_base
```

This maps directly onto the existing `BAwLPreparedRow`:

```text
v0    = kappa * v_base
slope = tau * v_base
sv    = sv_base          (residual SD, unchanged)
```

Build one `BAwLTimeGeometry`/`BAwLPreparedRow` per loaded row and evaluation
time, and update only the conditional mean per node — the same
node-invariant preparation discipline the BAwLCorr fused fallback uses. Do
not invent a second conditional-view representation.

Do not replace the conditional representation with only:

```text
sv_eff = sqrt(sv^2 + v^2 * tau^2)
```

That expression is a marginal variance identity. It loses the shared A/B
dependence if used in an ordinary independent logical-rule likelihood.

### 5. Existing fast-path structure and where the real cost is

The current LogicalRules fast path in `c_log_likelihood_logicalrules()` is
aggressively batched and is nowhere near per-trial factor-neutral. It
performs:

1. one raw `dfun`/`pfun` batch sweep over all rows at the observed RTs
   (`model_dfun_raw`/`model_pfun_raw` with direct column pointers), and
2. a 31-point Gauss-Legendre batch pre-pass (`use_gl_pass`) computing the
   channel quantities `GA_no = P(n_A wins channel A before RT)` and `GB_no`
   across all unique trials in vectorised sweeps.

Making this factor-neutral is the dominant work item of the whole plan, and
its cost lands only on AB trials:

- **Detection (`OR_DETECTION_ANALYTIC`/`OR_DETECTION_GNG`)**: no time
  integral is needed for finite RT. The conditional integrand is assembled
  from prepared-row endpoint evaluations at the observed RT, so the per-node
  cost is a handful of endpoint calls — comparable to one BAwLCorr fused-GH
  node.
- **Choice rules**: conditional on `z`, the channel quantity becomes

  ```text
  G_no(t | z) = integral_lo^t f_nA(s) * S_A(s | z) ds,
  ```

  which must be re-evaluated per GH node because `S_A(s|z)` depends on `z`.
  That is `n_nodes x n_GL` endpoint evaluations per channel per AB trial.
  Two facts keep this tractable and must be built in from the start, not as
  a later optimization:

  - `f_nA(s)` at the fixed GL abscissae is node-invariant — compute it once
    per trial and reuse it across all nodes. Only `S_A(s|z)` varies.
  - `S_A(s|z)` at each GL abscissa should come from a `BAwLPreparedRow`
    prepared once per abscissa (time geometry fixed), with only the
    conditional mean updated per node.

  With both in place the per-node incremental cost is `n_GL` prepared
  survivor evaluations per channel, which is the same order as the BAwLCorr
  generic-GH fallback and acceptable for AB trials only.

Single-target trials (and all trials when `tau = 0`) stay on the existing
batched sweeps unchanged; the factor-neutral evaluator is entered per trial
by route, exactly as BAwLCorr routes `exact_pair` trials off the GH path.

## Simulator implementation

### 1. Reference simulator

Implement a reference R simulator that samples one factor per unique trial,
not one factor per accumulator row:

```text
if AB:
    z <- rnorm(1)
    target drift means <- (kappa + tau*z) * base target means
else:
    target drift means <- base target means
```

Draw residual drifts independently for A, B, n_A, n_B, and (for GNG) nogo.
Then pass the actual drift draws to the appropriate subrace/rule simulator.

`rBAwL_corr()` in `R/model_LBA.R` is the template for one-factor-per-trial
sampling with joint positive-set reweighting — but its loading convention is
different and must not be copied: BAwLcorr loads on the sv scale
(`mu + sign(rho)*sv*loading*z`), while capacity loads on the drift-mean
scale (`(kappa + tau*z) * mu`). Reuse the structure (per-trial factor,
active-set-only conditioning), not the loading arithmetic or the
`rho -> loading` Ttransform.

For positive-drift models, condition only on active racers. A reference
implementation may use joint rejection or a multivariate truncated-normal
draw, but it must not condition on inactive fixed roles.

### 2. Choice-task and GNG activity

On an A-only choice trial, draw A and n_A only. On a B-only trial, draw B and
n_B only. On an AB choice trial, draw all four active racers, with the shared
factor loading only on A and B. For `OR_DETECTION_GNG`, additionally draw the
nogo racer on every trial; it never receives the factor loading.

### 3. Simulation/likelihood parity

The simulator and likelihood must agree on all of the following:

- whether capacity is active only in AB;
- whether `kappa` multiplies target base means or is represented as a design
  effect;
- whether drift positivity is jointly conditioned;
- which rows are active under each condition;
- whether `tau` is a residual/common-factor effect rather than a marginal
  independent `sv` effect.

## Likelihood implementation phases

### Phase 1: parameter and design contract

1. Add a model constructor/variant for LogicalRules capacity correlation.
2. Define `kappa` and `tau` transforms, defaults, bounds, and exceptions.
   Use the same exception mechanism the correlated BAwL constructor uses for
   `rho = 0` (`exception = c(rho = 0)` in `R/model_LBA.R`) so `tau = 0` is an
   exactly representable value and the ordinary route is recovered
   bit-for-bit; `kappa` defaults to exactly 1 the same way.
3. Require stable fixed accumulator roles for the relevant model variants.
4. Validate condition coding and A/B role availability (and nogo for GNG).
5. Validate that capacity parameters are shared within each trial.
6. Define whether an additive `DT` drift effect is used instead of, or in
   addition to, multiplicative `kappa`.

### Phase 2: factor-neutral logical-rule evaluator

Refactor the logical-rule likelihood so that the core per-trial evaluator can
consume conditional parameter views and return an unnormalized event mass.
This is the dominant work item — see architecture section 5 for what the
current batched fast path does and what must change. The refactor must leave
single-target and `tau = 0` trials on the existing batched sweeps; the
factor-neutral evaluator is entered per trial by route.

It must support:

- detection analytical finite-RT density (with and without nogo);
- choice/logical-rule finite-RT density;
- lower and upper censoring;
- truncation-normalizer endpoint probabilities;
- positive-drift node normalizers (closed-form `q_AB` for the pair,
  univariate terms for independents);
- response identities and missing responses.

The evaluator must not perform the outer capacity integration itself.

### Phase 3: detection implementation

Implement the AB detection route first because it has the simplest integrand:

```text
conditional node:
    fA, FA <- target A at z      (prepared-row endpoint modes f0 / F0)
    fB, FB <- target B at z
    mass   <- fA * (1 - FB) + fB * (1 - FA)
integrate mass over z
```

Assemble the node terms from `BAwLPreparedRow` endpoint evaluators; do not
write new scalar LBA endpoint code. Add corresponding conditional survivor
and censoring functions (the conditional total survivor is
`S_A(t|z) * S_B(t|z)` times independent survivors). `OR_DETECTION_GNG` uses
the same A/B math with the nogo racer's density/survivor entering the rule
assembly at each node exactly as in the ordinary GNG likelihood. Validate
against direct numerical integration and Monte Carlo simulation.

### Phase 4: choice/logical-rule implementation

For each factor node:

1. Evaluate the A/n_A subrace conditionally.
2. Evaluate the B/n_B subrace conditionally.
3. Construct the logical-rule channel states.
4. Evaluate the logical-rule event mass.
5. Accumulate the node contribution.

The per-node channel quantities require the conditional Gauss-Legendre
integral from architecture section 5. Build the two mandatory caches in this
phase, not later: node-invariant `f_nA`/`f_nB` values at the fixed GL
abscissae computed once per trial, and per-abscissa `BAwLPreparedRow`
geometry with only the conditional mean updated per node.

Do not route the A/B target rows through the ordinary single-race correlated
pair-cause formula. They are correlated across the two subrace evaluators.

### Phase 5: shared-factor integration

Use the BAwLCorr shared-factor machinery from `src/gh_quad.h` and the
patterns in `src/particle_ll.cpp`:

- standard-normal factor nodes and weights (`gh_rule()`,
  `gh_standard_normal_weight()`), with the env node-count override moved to
  a shared name;
- scan/recentre (`AGHCenter`, `agh_center_from_scan()`) where needed; the
  analytic loaded-winner centre does not transfer — logical-rule integrands
  are not single-race winner integrands — so use scan-based centring
  initially;
- stable max-then-sum log-sum-exp accumulation;
- the closed-form `q_AB` positivity denominator (mathematical specification
  section 5) — no denominator quadrature for the pair;
- truncation through the z-integrated total survivor at both endpoints;
- env-gated route counters and node diagnostics, following
  `src/bawl_corr_counters.h` (and actually increment every counter that is
  exported — BAwLCorr's `bvn_corner_evaluations` is currently declared but
  never incremented; do not repeat that).

For `tau == 0`, route exactly to the ordinary logical-rule evaluator so the
baseline is bit-for-bit equivalent wherever the existing contract permits.

### Phase 6: optimization

Only after correctness is established (the `f_nA` GL-abscissa cache and
per-abscissa prepared geometry are Phase 4 requirements, not optional
optimizations):

- batch A and B conditional scalar-kernel evaluations;
- use specialized detection formulas for `OR_DETECTION_ANALYTIC`;
- reuse equal-parameter subrace shortcuts;
- add a fast all-finite route analogous to the completed BAwLCorr route;
- preserve the generic censoring/truncation fallback.

## Validation plan

### Mathematical and simulator checks

1. Empirical A/B drift covariance matches `mu_A * mu_B * tau^2`.
2. Empirical target marginal variances match
   `sv^2 + mu^2 * tau^2`.
3. n_A and n_B remain independent of the shared factor and of each other.
4. A-only and B-only distributions do not change when AB capacity parameters
   are varied.
5. The AB drift covariance disappears when `tau = 0`.
6. `kappa = 1, tau = 0` reproduces the ordinary model.

### Likelihood checks

1. Detection AB likelihood agrees with direct one-dimensional quadrature.
2. Detection likelihood agrees with Monte Carlo simulation over a grid of
   `kappa`, `tau`, drift means, and thresholds.
3. Choice likelihood agrees with simulation for OR, AND, XOR, and ID rules.
4. `OR_DETECTION_GNG` with capacity agrees with simulation, and its A/B
   conditional terms match `OR_DETECTION_ANALYTIC` values when the nogo
   racer is disabled.
5. The closed-form `q_AB` matches the quadrature identity
   `integral qA(z) qB(z) phi(z) dz` across the `kappa`/`tau`/drift grid,
   including near-singular `rho_AB`.
6. A/B target correlation changes AB trials but not single-target trials.
7. Changing n_A/n_B (or nogo) parameters does not alter the A/B covariance
   structure.
8. Role-order permutations produce identical likelihoods.
9. Fixed inactive roles are absent from race and positivity calculations.
10. Positive-drift likelihood normalization agrees with the simulator’s
    jointly truncated active-set semantics.
11. Censoring and truncation agree with finite-RT likelihoods in limiting
    cases.
12. Contaminant and expansion corrections are applied once after full trial
    assembly.
13. Single-target and `tau = 0` trials produce bit-for-bit the existing
    batched fast-path values, and route counters confirm they never enter
    the factor-neutral evaluator.

### Regression and routing checks

Track route counts for:

```text
ordinary single/tau-zero
capacity detection pair
capacity logical pair
invalid
```

Verify that only A and B carry capacity loadings, while n_A, n_B, and the
GNG nogo racer remain independent singleton racers throughout the new model.

## Final design principle

The BAwLCorr plan supplies the mature shared-factor integration, numerical
stability, normalization, layout, and diagnostics architecture. This plan
supplies a different conditional event evaluator:

```text
BAwLCorr:       correlated pair -> one race likelihood -> integrate
LogicalRules:   correlated A/B subraces -> logical rule -> integrate
```

The latent-factor mathematics is shared; the race/logical-rule assembly is
not.
