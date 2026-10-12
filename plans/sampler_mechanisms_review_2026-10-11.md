# Sampler mechanisms review, 2026-10-11

Review of HEAD `55b29f4d`, concentrating on the October 8–11 optional-move,
allocation, transport, and scheduling commits. The uncommitted worktree changes
are the DDM variability/bin likelihood rewrite and its tests. No sampler code
was changed by this review.

The main finding is that several allocation and feature-gating errors can
explain unexpected activation and cost before invoking difficult posterior
geometry. Beyond those errors, lag-one jump statistics are insufficient to
predict long-run mixing in a hierarchy with several timescales. Fix the
accounting and gating first, then assess that prediction model separately.

Evidence is source inspection, existing experiment reports/CSV files, and the
companion arithmetic-only script. No fits, continuations, likelihood probes,
or sampler transitions were launched. An isolated package build begun before
the request to avoid competing with in-flight benchmarks was interrupted;
package installation and test suites were not completed. The existing graph
does not contain the new allocator/shear/HMC implementation, so findings were
checked directly against source. Arithmetic reproductions do not establish
the size of an effect in a real fit.

## Findings requiring fixes

### 1. P1: reachability ignores collapsed-only movement

`R/move_allocation.R:438` and `:828`, introduced by `55b29f4d`.

`.emc_alloc_range()` considers a fixed group unit unreachable. It does not
consider the joint collapsed move or the group credit from an HMC regime.
`.emc_alloc_solve_units()` freezes this incomplete mask before solving any of
the alternatives. A group coordinate with no effective shear is therefore
excluded even when the collapsed move can improve it. Within a direct solve,
lowering a fixed unit's rho with a collapsed candidate still does not make that
unit movable: its best rho remains equal to its candidate level.

The companion probe uses exactly the fixed-group setup of the existing test
at `tests/testthat/test-collapsed-move.R:442`: rho `(0.95, 0.9)`, collapsed
factors `(0.3, 0.4)`, and round cost `0.05`. Current code buys zero collapsed
rounds and reports `tau0 = 1`, although `tau0_all = 39`. Including the two
reachable units buys four rounds. The existing test's positive-round and
reached-rho assertions are contradicted by this arithmetic reproduction;
the complete test file was not run.

Compute reachability from the union of available kernels, including the
joint move and regime effects, and use that same coordinate set for every
candidate. Retain the distinction between a genuinely unreachable coordinate
and one excluded merely because a particular kernel cannot move it.

### 2. P1: zero measured group movement is turned into positive movement

`R/move_allocation.R:602` and `:638`; interacts with `55b29f4d`.

Shear and collapsed factors are capped at `1 - 1e-6`, whereas particle and
subject-HMC factors can equal one. Thus a measured zero jump becomes
`g = 0.999999`, passes the new reachability threshold, and is treated as a
kernel capable of improving the coordinate. This is consequential when the
rest of the iteration also has zero movement against positive pooled variance.

With baseline rho one, zero optional-move ESJD, and shear cost `0.05`, the
current solver buys the maximum four rounds and predicts a gain of
`99.99997%`: its baseline tau is approximately `2e12`, while the fabricated
movement produces tau approximately `5e5`. This is numerical regularization
being interpreted as statistical evidence, rather than a small harmless clamp.

Preserve an exact factor of one for measured zero movement. Keep numerical
guards inside tau/inverse calculations and make exploration an explicit policy
for unmeasured kernels. Add coverage through the measurement path, not only
solver tests that supply `g = 1` directly.

### 3. P1: HMC regimes predict a different schedule from execution

`R/move_allocation.R:797`, `:841`, and `:934`, restored in `416191b2`.

The regime first multiplies subject rho by `1 - r(1-h)` and charges `r*c_h`.
The solver then retains the subject's HMC factor and can buy another rate `x`,
multiplying rho again by `1 - x(1-h)` and charging another `x*c_h`. Execution
instead assigns a single rate `min(1, r+x)` and runs at most one transition
in the final HMC round.

A deterministic probe selects regime one plus supplemental rate about
`0.430`. It predicts HMC factor `0.1145`, but execution supplies `0.1772`.
It charges approximately `0.1787` HMC work units for work actually costed at
`0.1250` in the same model. Particle substitution can consequently be chosen
against movement that execution never supplies. This affects the default
regime-enabled allocator; extra interleaving is not required.

Use one total rate with a regime lower bound and price its incremental
movement consistently, or explicitly schedule the two modeled transitions.
The choice requires a corresponding adjustment to the group-credit model.
Assert agreement between solver movement/cost and the actual schedule.

Two related defects need separate coverage:

* At `:884–890`, calibration reads only `sol$x`, losing `sol$regime`. The
  optimistic probe selects regime one with supplemental rate zero, so its
  calibration gate is false despite an intended positive HMC rate. A newly
  eligible subject outside burn can therefore miss bracketing/pilots and be
  disabled for lack of a usable epsilon. Select pilots from total scheduled
  rates and the actual subjects covered by the regime.
* At `:811`, the regime numerator averages only subjects with usable
  deviation factors, whereas `r_bar` at `:679` averages all finite subject
  deviation coordinates. Mixed eligibility changes the population. In the
  probe, even calling the projection with rate zero lowers group rho from
  `0.9` to `0.6333`. Keep unsupported subjects' unchanged deviations in the
  numerator and assert that a zero-rate regime is identity.

### 4. P1: the curve switch does not disable a checkpoint's existing atlas

`R/shear_move.R:488`, `:598`, `:667`, and `:952`; switch added in `416191b2`.

`emc2.shear_curve = FALSE` prevents the initial pilot from requesting an
atlas. But `.emc_shear_atlas_add()` immediately includes any chain that already
has one, `.emc_shear_build_curves()` still builds the regression, and the
transport context still accepts the resulting curve. The six-node arithmetic
probe builds an active curve with the option false.

This can invalidate a resumed "curves off" arm if its checkpoint contains an
atlas. It can also keep paying for atlas top-ups and curve validation. Gate
atlas fitting, curve construction, and curve consumption consistently while
preserving cached atlas data if desired. Before interpreting in-flight arms,
inspect the restored checkpoint and resulting active frames. No running arm
was changed by this review.

### 5. P2: interleaved HMC is charged to shear and omitted from HMC evidence

`R/sampling.R:1082–1096`, compared with `:1150–1201`, added in `416191b2`.

The intermediate `.emc_hmc_round()` results update proposals, streams, and
the worker plan but discard `hr$tasks`. Only the final round calls
`.emc_hmc_record_tasks()` and passes tasks into allocation cost recording.
Intermediate HMC time sits inside `shear_el`, inflating the inferred shear
cost while leaving HMC attempts, acceptance, jump evidence, and reported HMC
work incomplete. Warmup epsilon changes made in intermediate rounds are
also not recorded through the normal evidence path.

With interleave K and fixed subject rate r, actual expected attempts are
`K*r`, but the final-round counter sees only `r`. Splitting rounds also
changes the per-call move-count histories used for shear saturation; they
cease to describe the original one-call conditional kernel. Interleaved HMC
even runs between parts when the allocated shear rates are all zero.

Record every HMC round and separate its cost before deriving shear cost.
Either make the allocator explicitly model interleaving or label such arms
as controlled schedules outside its current predictive assumptions. This
finding applies when `emc2.shear_interleave > 1`; the default is one.

### 6. P2: shared particle rates undo locally priced HMC substitution

`R/move_allocation.R:984–992`, originating in `509cfac2`.

Each chain optimizes its own particle/HMC pair, then every chain's particle
rate is replaced by the cross-chain maximum. HMC rates and recorded objective
and cost are not reoptimized. If one chain selects particle rate zero with
HMC replacing it, while another selects particle rate one, the first chain
actually pays for both kernels. The cost that justified substitution is no
longer the cost of its executed plan.

The recently removed collapsed-rate synchronization had the same type of
problem. If common particle rates are required operationally, optimize them
jointly and recompute each chain's predictions. Otherwise preserve local
rates. Existing fields distinguish `particle_rate` from `particle_shared`,
which permits checking saved allocation histories without new sampling.

## Trigger and prediction assumptions to assess

### Curve triggering uses a different variance reference from allocation

`R/shear_move.R:473–477` uses the current chain's proposal-window variance;
`R/move_allocation.R:590–604` uses pooled-chain variance and, in production,
can use all retained production draws. A chain confined to a narrow part of
the posterior can therefore look productive to the curve trigger but
unproductive to the allocator.

The same per-round ESJD `0.01` gives factor `0.625` in a short own-chain
window of variance `0.01333`, but factor `0.995` against pooled variance one.
The first fails the curve trigger's `g >= 0.9`; the second clearly describes
ineffective movement. The companion script verifies the failed trigger.
The trigger also requires five measured rounds, `lam >= 0.5`, a cadence
check, and a subsampled pilot gain at least two. A useful diagnostic should
record the result of each gate rather than report only "curve present".

Align the intended variance reference, or document why these two detectors
deliberately answer different questions. A pooled denominator can be inflated
by nonconvergence, so alignment alone is not a guarantee of useful curvature.

### Lag-one ESJD is not integrated autocorrelation time

`R/move_allocation.R:245–253` converts a jump-based lag-one estimate to
`tau = (1+rho)/(1-rho)`. This is the AR(1) formula. General ESS depends on
the sum over all autocorrelation lags; see the primary-source
[Stan reference manual](https://mc-stan.org/docs/2_40/reference-manual/analysis.html).

For a stationary coordinate containing 20% variance from an AR mode with
correlation `0.99` and 80% from a mode with correlation `0.1`, lag one is
`0.278`. The allocator's tau is `1.77`; the true integrated tau is `40.78`.
Removing only the fast mode's autocorrelation and increasing cost 5% gives a
predicted tau-cost improvement exceeding the allocator's 10% gain threshold,
while the actual tau-cost worsens from `40.78` to `42.63`.

This is an exact arithmetic counterexample, not a diagnosis of a particular
saved fit. It demonstrates why better jumps can coexist with worse ESS per
second. Multiplying separately measured rho factors for different kernels
requires further assumptions about shared slow modes; kernel invariance
does not supply those assumptions. The current dense-grid tests validate
optimization of that model, not its predictive accuracy on composed kernels.

Group allocation also monitors mean directions and marginal log SDs rather
than every free covariance entry or combined mean/scale/subject slow mode.
The continuation driver already scores fixed combined group directions;
comparing predictions with those diagnostics is preferable to adding another
local jump heuristic without validation.

### The optional refresh floor guarantees an expected count, not a refresh

`R/move_allocation.R:493–501`, particle selection at `R/sampling.R:930`,
and HMC selection at `R/subject_hmc.R:302`.

Under independent Bernoulli selection, `q+r >= 1` does not guarantee a
transition is attempted every iteration. At `q=r=0.5`, neither runs with
probability `0.25`. The current factor product matches independent selection;
coupling selection to guarantee coverage changes that formula too. Clarify
whether the switch intends an expected-work floor or actual coverage, then
make scheduling, prediction, and documentation agree. This option is not
enabled by default.

## What the existing experiments establish

The chronological report in
`plans/collapsed_transport_adaptation_2026-10-09/report.md` explicitly says
its exact-pilot scale searches and several proposed policies were reverted.
Those sections must not be interpreted as the current implementation.

Its later original-scale, support-conditioned LBA screen still found expected
exact acceptance approximately `4.9e-34`. Its subject-block prototype
recovered acceptance around `0.468` but had worse ESS per second in the short
continuation and was removed. These observations support accumulated
surrogate error and inadequate long-range movement as unresolved problems;
they do not justify simply reinstating global scale calibration or the block
prototype. The large earlier support failures are documented separately in
`plans/lba_sampling_audit_2026-10-09/report.md`.

The historical T30/corr15/Forstmann comparisons support useful movement in
some targets and substantial per-parameter tradeoffs. They predate the latest
reachability change, and the benchmark README mixes several older kernel
variants in its initial tables. No general current-HEAD efficiency ranking
follows from those tables. In-flight results need source/configuration
provenance and effective feature state before being combined with them.

## Worktree likelihood changes

The DDM rewrite removes an adaptive cubature and introduces fixed rules for
start-point/onset variability and CDFs. Smooth likelihood values matter to
finite-difference HMC forces, so this can address a different source of poor
HMC efficiency. It does not repair the allocator/schedule issues above.

The added tests check smoothness and bin/CDF consistency against integration
of the new unbinned density. They do not independently validate that density
over the full parameter domain. Document the need for a high-precision
independent density/CDF comparison, limit checks as variability approaches
zero, and checks across the small/large-time series and bin-rule switches.
HMC gradient-scale checks should include `h/2`, `h`, and `2h` on fixed stored
states; likelihood error differences are divided by `h`. No independent
numerical validation or DDM efficiency benchmark was run in this review.

## Required follow-up after the in-flight benchmarks

First add arithmetic/regression coverage for collapsed-only reachability,
zero-ESJD measurement, total-rate regime accounting, regime-only pilot
selection, mixed eligibility, and curves-off restored checkpoints. Add a
small mocked scheduler check that all interleaved tasks and costs are
recorded, and a plan/execution check for shared particle rates. Run the
existing collapsed/allocation suites after fixing the latest regression.

From existing saved histories, compare local and executed rates, active
curve/atlas state, total attempt counters, and allocation predictions with
observed long-lag and fixed-direction ESS. Record excluded coordinates and
every curve gate. Resumed curves-off arms with an active curve, or interleaved
arms relying on current HMC counters, need qualification before interpretation.

Only after these checks, run the existing fixed-length matched-checkpoint
protocol, one job at a time, with several streams and explicit source hashes,
options, inherited checkpoint configuration, and active feature counts.
Compare CPU/wall cost, bulk/tail ESS, and R-hat in slow fixed directions as
well as per-parameter families. Keep exploration/calibration work in the
cost. Include a frozen-plan control to separate adaptive-plan errors from
kernel performance. These are documented requirements, not jobs launched by
this review.
