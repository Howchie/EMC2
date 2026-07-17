# Plan: cross-subrace correlated-target capacity for LogicalRules

Date: 2026-07-17  
Related architecture: `bawl_corr_exact_kernel_plan.md`  
Primary targets: `LogicalRulesLBA` and `OR_DETECTION_ANALYTIC`

## Outcome

Add a redundant-target capacity model in which the latent drift draws for
target A and target B are correlated when both target stimuli are present,
while the two targets remain in separate logical-rule subraces.

The model supports:

- a mean capacity effect (`kappa`),
- a between-trial capacity-variability effect (`tau`),
- either effect alone or both together,
- detection rules with only A/B target detectors,
- choice rules with A/n_A and B/n_B subraces,
- independent n_A and n_B racers,
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
4. Independent n_A and n_B racers in choice models.
5. Correct positive-drift normalization for the active set.
6. Finite RTs, censoring, truncation, omissions, contaminants, and expansion
   through the normal logical-rule likelihood contracts.
7. A reference R simulator and a production C++ likelihood.

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

For detection models, the active target set is:

```text
S == "A"   -> A
S == "B"   -> B
S == "AB"  -> A and B
```

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

For an AB detection trial, the factor-level positivity mass is:

\[
q_{AB} = \int q_A(z)q_B(z)\phi(z)\,dz,
\]

where \(q_i(z)\) is the conditional probability that the relevant target
drift is positive.

For an AB choice trial, independent n_A/n_B positivity terms may multiply the
conditional node contribution, but A/B must retain their joint factor
integration.

For A-only trials, only the active A subrace participates in the positivity
normalizer. B must not be included merely because its fixed role exists.

The final per-trial contract is:

```text
log numerator = log ∫ conditional_event_mass(z) phi(z) dz
log denominator = log ∫ conditional_observation_mass(z) phi(z) dz
trial log likelihood = log numerator - log denominator
```

Each normalizer is applied exactly once.

## Architecture

### 1. Reuse from the completed BAwLCorr architecture

Assume the final implementation of `bawl_corr_exact_kernel_plan.md` exists.
Reuse its final architecture for:

- shared state constructed once outside the particle loop;
- direct parameter-column access through the final parameter table;
- canonical per-trial layout and route classification;
- shared-factor GH rules, weights, centering, refinement, and log-sum-exp
  accumulation;
- positive-drift node normalizers and truncation-normalizer contracts;
- probability/tail primitives and numerical-status handling;
- component-level likelihood values that preserve their own normalizer;
- per-trial output and expansion handling;
- route counters and diagnostic instrumentation.

Do not reuse the BAwLCorr assumption that `lM` identifies the correlated pair.
For this model, pair membership is explicit: the correlated pair is A/B.

The final BAwLCorr exact pair kernel is for a pair of racers within a common
race and should remain unchanged. It is not the conditional logical-rule
integrand for the choice model. The reusable numerical abstraction is the
shared-factor integrator and its component/normalizer interfaces.

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
    row_nA, row_nB
    active_A, active_B, active_nA, active_nB
    route
    positive_drift_route
    valid
```

Routes should at least distinguish:

```text
ordinary_single_or_tau_zero
capacity_detection_pair
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

The exact implementation may use direct column views, prepared row objects,
or a compact per-node buffer, depending on the completed shared-factor API.

Do not replace the conditional representation with only:

```text
sv_eff = sqrt(sv^2 + v^2 * tau^2)
```

That expression is a marginal variance identity. It loses the shared A/B
dependence if used in an ordinary independent logical-rule likelihood.

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

Draw residual drifts independently for A, B, n_A, and n_B. Then pass the
actual drift draws to the appropriate subrace/rule simulator.

For positive-drift models, condition only on active racers. A reference
implementation may use joint rejection or a multivariate truncated-normal
draw, but it must not condition on inactive fixed roles.

### 2. Choice-task activity

On an A-only choice trial, draw A and n_A only. On a B-only trial, draw B and
n_B only. On an AB choice trial, draw all four active racers, with the shared
factor loading only on A and B.

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
3. Require stable fixed accumulator roles for the relevant model variants.
4. Validate condition coding and A/B role availability.
5. Validate that capacity parameters are shared within each trial.
6. Define whether an additive `DT` drift effect is used instead of, or in
   addition to, multiplicative `kappa`.

### Phase 2: factor-neutral logical-rule evaluator

Refactor the logical-rule likelihood so that the core per-trial evaluator can
consume conditional parameter views and return an unnormalized event mass.

It must support:

- detection analytical finite-RT density;
- choice/logical-rule finite-RT density;
- lower and upper censoring;
- truncation-normalizer endpoint probabilities;
- positive-drift node normalizers;
- response identities and missing responses.

The evaluator must not perform the outer capacity integration itself.

### Phase 3: detection implementation

Implement the AB detection route first because it has the simplest integrand:

```text
conditional node:
    fA, FA <- target A at z
    fB, FB <- target B at z
    mass   <- fA * (1 - FB) + fB * (1 - FA)
integrate mass over z
```

Add corresponding conditional survivor and censoring functions. Validate it
against direct numerical integration and Monte Carlo simulation.

### Phase 4: choice/logical-rule implementation

For each factor node:

1. Evaluate the A/n_A subrace conditionally.
2. Evaluate the B/n_B subrace conditionally.
3. Construct the logical-rule channel states.
4. Evaluate the logical-rule event mass.
5. Accumulate the node contribution.

Do not route the A/B target rows through the ordinary single-race correlated
pair-cause formula. They are correlated across the two subrace evaluators.

### Phase 5: shared-factor integration

Use the completed BAwLCorr shared-factor integration service:

- standard-normal factor nodes and weights;
- adaptive scan/refinement where needed;
- stable log-sum-exp accumulation;
- separate numerator and positivity/truncation-denominator integrations;
- route counters and node diagnostics.

For `tau == 0`, route exactly to the ordinary logical-rule evaluator so the
baseline is bit-for-bit equivalent wherever the existing contract permits.

### Phase 6: optimization

Only after correctness is established:

- cache n_A/n_B singleton quantities that are invariant across factor nodes;
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
4. A/B target correlation changes AB trials but not single-target trials.
5. Changing n_A/n_B parameters does not alter the A/B covariance structure.
6. Role-order permutations produce identical likelihoods.
7. Fixed inactive roles are absent from race and positivity calculations.
8. Positive-drift likelihood normalization agrees with the simulator’s jointly
   truncated active-set semantics.
9. Censoring and truncation agree with finite-RT likelihoods in limiting cases.
10. Contaminant and expansion corrections are applied once after full trial
    assembly.

### Regression and routing checks

Track route counts for:

```text
ordinary single/tau-zero
capacity detection pair
capacity logical pair
invalid
```

Verify that only A and B carry capacity loadings, while n_A and n_B remain
independent singleton racers throughout the new model.

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
