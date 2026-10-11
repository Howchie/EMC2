# Occasional subject HMC with adaptive allocation

Implementation plan, 2026-10-07. Baseline: branch `playground`, commit
`49bfe8d0`. This document consolidates the original proposal and its reviews;
all implementation decisions and acceptance criteria are stated here.

## 1. Objective and evidence

Add occasional subject-conditional HMC to rescue chains trapped along curved
subject posterior ridges, including traps that develop during production.
Allocate its cost using a deterministic mixing-efficiency surrogate, and judge
success by convergence and minimum ESS per total CPU on a fixed test panel.
HMC rates and tuning remain adaptive during production, on the same doubling
schedule as the other proposals.

Distinguish two failure mechanisms when diagnosing results:

- **Coupling between hierarchy levels.** Group means and weakly informed
  subjects move slowly together. The reported T30 case has identified subjects
  and known correlations: default fits without a group move and with ensemble
  updates both exhausted `max_tries`, with mu R-hat about 1.32. Group moves
  address this mechanism. Subject HMC may influence subsequent group updates,
  but a subject-only transition is not a direct group transport.
- **Poor exploration of a subject conditional.** In the reported lba71sim
  checkpoint, subject 11, chain 1 stayed in a tail region through production.
  Multi-start conditional optimisation found one mode. With the group fixed,
  the particle step produced sd(A) about 0.05 over 1,500 iterations against
  conditional sd about 0.167; finite-difference HMC reached the bulk in about
  50 iterations. The A/B/t0 tradeoff motivates following a curved trajectory
  within the subject conditional.

These are existing observations to reproduce, rather than new benchmark
results. Preserve the separation throughout the experiments. Weak
identification, a second mode, and conditional trapping require different
interpretations; low R-hat alone does not distinguish them.

## 2. Scope and sampler contract

The first implementation supports `standard` and `diagonal-gamma` hierarchies
with a common Gaussian subject prior, including covariance blocks represented
by `standard`. Initially exclude group designs, latent-factor/SEM hierarchies,
nuisance samplers, marginalized parameters, pooled coefficients, and joint
likelihood bookkeeping. Use an explicit support predicate covering these
features and parameter ordering. Extending support later requires deriving
and testing each configuration's exact conditional.

For eligible subject s, hold the group state (mu, V) fixed throughout the
particle update and following HMC transition. The target in the sampler's
existing alpha coordinates is

    log pi_s(a) = ll_s(a) - 0.5 * (a - mu)' V^-1 (a - mu) + constant.

Use the same likelihood, transforms, support rules and subject prior as the
particle kernel. Audit any transform Jacobians already included in this target;
add no duplicated Jacobian. A constant affine whitening Jacobian cancels in the
HMC acceptance ratio.

Each subject takes its existing particle update, then takes HMC before the
ensemble group update and any separately enabled group move. The first regular
subject task after plan initialization or a rebuild makes one scheduled HMC
measurement attempt for every eligible subject. On the remaining iterations,
HMC is selected with probability r_s. The adaptive rate controls additional
attempts; a zero rate does not remove the scheduled measurement.
Return the final alpha and its exact likelihood together so the ensemble pool
includes the post-HMC current state. Measure the particle and HMC jumps
separately, before ensemble reselection can change the subject again.

A fixed, state-independent rate mixes identity with a valid HMC kernel, and
composition with a valid particle kernel preserves the conditional target.
A scheduled measurement applies that same valid kernel at a predetermined
iteration; it also preserves the conditional target for a fixed plan.
The group-dependent metric below is also valid because the group is held fixed
for the conditional update. History-based tuning is adaptation; its separate
production contract is in section 5.

## 3. The conditional HMC kernel

### 3.1 Metric and coordinates

Use the finite canonical likelihood precision H_s stored in `lik_laplace`.
`R/mode_start.R` already projects this precision to positive semidefinite in
`.emc_lik_canonical()`. Do not require a positive-definite likelihood precision
for routine HMC eligibility: a valid prior precision makes

    Q_s = H_s + V^-1
    Sigma_s = Q_s^-1

positive definite even when H_s is singular. Validate dimensions, parameter
names, finiteness, symmetry and factorization. A failed fit or factorization
makes that subject ineligible for the current plan version, with a recorded
reason. Reconsider eligibility at subsequent proposal rebuilds. Do not
silently substitute a prior-only metric for an informative failed likelihood
fit; a fallback needs its own performance evidence.

Freeze H_s between proposal rebuilds. For each transition compute a factor S
such that S S' = Sigma_s under the current group state, then keep it fixed
throughout the trajectory. Implement

    a = mu + S z
    momentum v ~ N(0, I)
    K(v) = 0.5 * v' v.

Cholesky solves avoid forming inverses in the loop. This is equivalent to
M^-1 = Sigma_s in alpha coordinates. Evolve z with unit mass; no matrix
factorization or metric adaptation occurs inside leapfrog.

### 3.2 Deterministic finite-difference force

For each coordinate i, evaluate the likelihood at a +/- h_s S e_i:

    g_ll,z[i] = {ll_s(a + h_s S e_i) - ll_s(a - h_s S e_i)} / (2 h_s)
    g_prior,z = -S' V^-1 (a - mu)
    g_z = g_ll,z + g_prior,z.

Batch the 2p probes in one likelihood call per gradient. Include the centre
likelihood at each new position for support checks and reuse it for the final
acceptance calculation. Reuse the particle step's starting likelihood, but
compute the starting HMC gradient explicitly.

h_s is dimensionless in whitened coordinates and fixed for a plan version.
Use 1e-5 and 1e-4 as initial sensitivity candidates, not as universal accuracy
guarantees. Before accepting a default, compare force stability over adjacent
scales on the Gaussian control, curved toy, representative likelihood sizes,
and states near support boundaries. Check for cancellation and model-level
numerical discontinuities. Choose the smallest stable error region supported
by those checks, and record the choice. There is no position-dependent probe
shrinkage, random likelihood approximation, or history-dependent force within
a trajectory.

Leapfrog momentum kicks are volume-preserving shears even with a deterministic
approximate force; symmetric kicks make the map reversible. Therefore an
approximate force can be corrected with the exact Hamiltonian acceptance
ratio. This construction follows the reversible, volume-preserving proposal
principle in [Neal, section 5.5.3](https://arxiv.org/abs/1206.1901); finite
differences affect efficiency rather than substituting an approximate target.

### 3.3 Trajectory, acceptance and failure handling

Draw eps = eps_s * U, U uniform on [0.9, 1.1], independently of the current
subject position. Set L = max(1, ceil((pi/2)/eps)). Hold eps and L fixed for the
entire trajectory. Integration time is approximately a quarter period of the
whitened Gaussian surrogate, a starting design choice to validate on curved
posteriors. It is not a derived optimum for every target.

Apply the usual half momentum kick, L position steps with intervening kicks,
and final half kick. Include momentum reversal in the proposal definition.
Accept with

    log acceptance = min(0, log pi_s(a_new) - log pi_s(a_old)
                            + K(v_old) - K(v_new)).

Use exact endpoint likelihoods and the analytic Gaussian prior. Record energy
error and acceptance probability as diagnostics, even though acceptance rate
is not the tuning objective.

Abort and reject the entire trajectory when a checked centre or probe lies
outside support or has a non-finite likelihood. Apply the same position-domain
check at every gradient location, including the start and endpoint: a completed
forward path and its reverse visit the same checked positions. Verify this
symmetry with boundary tests. Never accept a partial path, project onto a
boundary, or retry a smaller eps inside the same transition. Numerical failures
retain the post-particle state and its likelihood; programming errors follow
the existing failure policy rather than being hidden as rejections.

### 3.4 Execution, RNG and cost

Run the whole particle-plus-HMC subject task in its assigned worker. Keep all
leapfrog batches local to that task. There is one subject task round trip,
rather than one synchronous master/worker round trip per leapfrog step.

Use the subject's `rng$subjects[[s]]` stream for selection, jitter, momentum
and acceptance, in a documented order shared by serial, static and queued
paths. Transport the final RNG state with alpha, likelihood and diagnostics.
Fallback recomputation starts from the original task input and RNG state; it
must not duplicate an accepted transition or advance the stream twice. With
HMC disabled or the hierarchy unsupported, consume no additional RNG draws.

Count actual evaluated parameter rows, likelihood batches, completed position
steps, early aborts and calibration work. With the implementation above, a
complete uncached trajectory costs

    c_H = 2p * (L + 1) + L

likelihood rows: L + 1 gradient batches and L new centre evaluations. The
starting centre likelihood is cached. For p = 71 and L = 30 this is 4,432 rows,
not 4,261. This example is a cost illustration, not a predicted trajectory
length; L is determined by the calibrated eps. Count shortened paths and any
extra evaluations directly rather than charging this formula unconditionally.

## 4. Measurement and allocation

### 4.1 What to record

All allocation state is chain-local. At each rebuild construct a reference
Sigma_ref,s = (H_s + V_bar^-1)^-1 from the window-mean group covariance and use
its eigenbasis to whiten jumps. Keep this diagnostic frame distinct from the
current group's integration frame. A Laplace width is a scale surrogate, not
proof of the true conditional variance.

Record raw alpha displacement vectors so a window's particle and HMC jumps can
be expressed in the same reference frame. For each kernel k:

    ESJD_k,s = min_i mean_t { (S_ref,s^-1 delta_a_k,s,t)[i]^2 }
    rho_k,s  = clip(1 - ESJD_k,s / 2, 0, 1).

Use an eigenfactor for S_ref so coordinates correspond to its eigenvectors;
set eigenvector conventions deterministically. The particle mean includes
zero jumps. For HMC tuning and allocation use acceptance probability times
the proposed squared displacement, including zero for aborted trajectories;
the mean includes scheduled and Bernoulli-selected HMC attempts, rather than
counting skipped iterations as HMC failures. Also
retain realised accepted jumps to check the acceptance-weighted estimate.

Whitening with the Laplace frame avoids scaling a trapped chain by its own
narrow empirical variance. It does not turn ESJD into an exact autocorrelation
estimate. The formula rho = 1 - ESJD/2 assumes unit marginal variance and a
stationary lag-one interpretation; clipping deliberately forgoes claimed
benefits from negative correlation. Keep per-coordinate measurements as well
as their minimum, and compare predictions against measured multi-lag behaviour
and `posterior::ess_bulk`.

Record the final group state of every complete sweep independently of stored
sample thinning and trimming. Include mu, log standard deviations, and every
free correlation; include fixed window-based linear directions of mu where
coordinate correlations hide slow movement. Exclude structurally fixed entries.
For the initial allocator use

    rho_G = maximum nonnegative lag-one correlation over these group coordinates
    tau_G = (1 + rho_G) / (1 - rho_G).

This is an empirical AR(1) surrogate for group mixing, not a measurement of
integrated autocorrelation time. Use final-to-final sweep pairs, including all
enabled group updates, rather than reusing the shear's pre-shear pairs without
changing their meaning. Subject HMC may change the following ensemble update,
so recompute this estimate at every rebuild and do not treat it as an immutable
physical floor.

C_base counts all non-HMC likelihood work, including particle retries,
ensemble pools and enabled group moves. Keep raw rows as the initial
reproducible cost unit; also record each subject's data size and batch sizes.
Raw row cost is a proxy for CPU, especially with heterogeneous subjects or
models. Introduce fixed data-based row weights only if profiling shows that
the raw proxy systematically chooses poor allocations, and evaluate that
change separately. Timings never drive a scheduling decision.

### 4.2 The allocation model

For ordinary iterations with the composed particle-then-optional-HMC update,
approximate

    rho_s(r) = rho_P,s * (1 - r + r * rho_H,s)
    tau_s(r) = (1 + rho_s(r)) / (1 - rho_s(r))
    C(r) = C_fixed + sum_s r_s * c_optional,s
    efficiency(r) = 1 / { max(tau_G, max_s tau_s(r_s)) * C(r) }.

If the next rebuild interval has B iterations, charge the mandatory measurement
work as sum_s c_probe,s / B, and let c_optional,s = (B - 1) * c_H,s / B.
C_fixed is C_base plus that measurement cost and any additional calibration
cost amortized over B. Thus the expected number of regular HMC attempts for
subject s is 1 + (B - 1) * r_s, or average rate

    r_effective,s = {1 + (B - 1) * r_s} / B.

An abort counts as an attempt, with its actual cost. For B = 1 there are only
scheduled attempts, so additional rates are zero. For an explicitly frozen
rebuild schedule, count the initialization measurement as a one-time cost;
there are no recurring rebuild measurements.

The ordinary-iteration rho model deliberately credits no mixing gain to the
scheduled measurement; its cost is still charged. It is a conservative proxy
for that part of the design, not an exact autocorrelation formula for the
periodic schedule. Do not substitute r_effective into the Bernoulli formula
and claim the schedules are equivalent. Validate actual block mixing and ESS
in the panel, including measurement-only allocations.

Costs and rho values are window estimates. This model treats slow directions
as shared scalar AR(1) modes and costs as additive; actual kernels need not
share modes or commute. It proposes an allocation, while measured minimum
ESS per CPU decides whether the method succeeds.

For each valid, calibrated subject with rho_P > 0 and rho_H < 1, and feasible
target tau (rho_star = (tau - 1)/(tau + 1)), the cheapest rate is

    r_s(tau) = max(0, (1 - rho_star/rho_P,s) / (1 - rho_H,s)).

Handle these cases explicitly:

- rho_P = 0: r = 0 under the surrogate; no division by zero.
- rho_H = 1: HMC adds no measured contraction; r = 0 until recalibration.
- 0 <= rho_H < 1: HMC can help even if rho_H >= rho_P. For example,
  rho_P = 0.2 and rho_H = 0.5 compose to 0.1 at r = 1. Comparing the two
  standalone rho values is not a valid exclusion rule.
- Unsupported or unusable subjects: r = 0 and their unassisted tau remains
  in the objective. Report when an unaddressable subject limits the fit.
- Missing or sparse HMC measurements: use the calibration policy in section 5,
  rather than interpreting missing evidence as proven failure.
- A measured rho of 1 gives infinite tau. Do not silently cap it to make
  progress look finite. Use the optimistic calibration test to investigate a
  supported trapped subject; return a documented unresolved bottleneck if
  the calibrated kernel still offers no finite allocation.

Define tau_s(0) and tau_s(1) for useful eligible subjects. For subjects that
cannot benefit, their attainable minimum is tau_s(0). Then

    tau_min = max(tau_G, max_s attainable_minimum_tau_s)
    tau_max = max(tau_G, max_s tau_s(0)).

Search feasible targets in [tau_min, tau_max] at the endpoints and the
subject breakpoints tau_s(0); redundant tau_s(1) candidates may be included.
Compute the actual rates, max tau and cost at each candidate. Reject a target
requiring r > 1 rather than clipping an infeasible target and pretending it
was attained. Handle infinite endpoints through their limits or explicit
infeasibility. Include the all-zero additional-rate allocation, which still
contains scheduled measurements, and compare its total cost with the original
sampler without HMC. Break ties toward lower HMC cost, then a fixed subject
ordering.

For finite inputs, on any interval with a fixed active subject set, write
C(tau) = A - b * (tau - 1)/(tau + 1), b >= 0. Then
(tau * C(tau))'' = -4b/(tau + 1)^3 <= 0, so its minimum is at an interval
endpoint. This justifies the breakpoint search for this surrogate. It does
not prove an optimal allocation for actual ESS, CPU, a stochastic trajectory
cost, or an unknown group response to HMC.

## 5. Calibration and adaptation through production

### 5.1 A shared rebuild schedule

Use the existing proposal rebuild opportunities in burn and adapt. In
production update the HMC plan at the same transition counts used by
`add_proposals()` and `.emc_next_sample_rebuild()`: by default 100, 200, 400,
800, ... . Honour explicitly configured rebuild schedules as well. Counts
are sampler transitions, independent of thinning and stored-window trimming.
Do not introduce a separate production freeze or a second adaptation clock.
Plan initialization and every rebuild arm one measurement attempt per eligible
subject in the next regular subject task. This is an actual Metropolis-corrected
transition, so its accepted state can immediately contribute to rescue. Retain
pending-measurement flags in checkpoints and consume each flag once, including
on rejection, so a resume or worker retry does not repeat it.

At each rebuild, version the likelihood metric, reference frame,
finite-difference scale, eps, cost estimates, calibration status and rates.
Compute the new snapshot from the preceding window and calibration results,
then activate it for subsequent transitions. Within the production block the
learned tuning and rates stay fixed; the scheduled measurement, group-conditioned
integration frame, independent jitter and Bernoulli selection follow the rules
above.
Statistics collected inside the block do not immediately change its kernel.

Every eligible subject is reconsidered, including r = 0 subjects. One fresh
attempt per interval is the measurement minimum, not enough by itself for a
reliable HMC performance estimate. Retain compatible earlier HMC measurements
with their counts and plan versions; invalidate estimates when changed tuning
makes them incomparable, and use additional calibration when needed. Production
must be able to activate HMC after a subject that was mixing well becomes
trapped, and to retune a previously ineffective trajectory at a later rebuild.
The resulting detection delay can be as long as the next doubling interval;
report it in the late-trap test rather than claiming immediate detection.

### 5.2 Calibration without starvation

Start eps near p^(-1/4), then bracket a numerically usable scale with bounded
doubling/halving trials. The one-step search in
[Hoffman and Gelman, algorithm 4](https://www.jmlr.org/papers/v15/hoffman14a.html)
is a useful initialization procedure. Its acceptance comparison establishes
a starting bracket; it is not the final efficiency objective or a requirement
to use dual averaging.

Use acceptance-weighted ESJD per actual likelihood cost to choose eps. Compare
the two independent jitter halves, retain their sample counts and uncertainty,
and make bounded log-eps updates at rebuilds using a gain decreasing with
informative tuning windows. Validate this stochastic optimization against a
small explicit eps grid before trusting its optimum. If the window cannot
resolve a slope, keep the previous eps and collect more calibration evidence.
Changes in L and early-abort cost must be included in the score.

An unmeasured HMC kernel has optimistic rho_H = 0 for assessing whether
learning its performance could improve allocation. Never replace this by
rho_H = 1 after a handful of failed trajectories. When the optimistic
allocation assigns positive rate to an uncalibrated subject, allocate a
bounded extra calibration batch independent of that subject's incidental
Bernoulli HMC count. Scheduled measurements supply ongoing evidence even for
zero-rate subjects; they do not replace initial step-size bracketing or adequate
calibration. Start with 20 pilot transitions, check sensitivity at 10 and 40,
and retain this provisional batch size only if it reliably avoids starvation
without dominating cost. This is a calibration parameter to validate, not a
claimed universal optimum or a subject-specific rescue threshold.

At production rebuilds run pilots on private copies of the subject state with
the group fixed; discard their states, retain tuning evidence, and charge all
likelihood work. Calibration does not silently replace stored posterior draws.
Only subsequent scheduled or Bernoulli-selected HMC transitions move the
production chain. Preserve and return the RNG stream advanced by those pilots
for reproducible resume.

Reassess dormant subjects under the optimistic allocation at later rebuilds,
especially when their particle movement has deteriorated or their plan has
changed. If learning HMC performance could improve the objective, recalibrate
instead of treating an old rho_H = 1 as permanent exclusion. Sparse or
obsolete evidence must not lock the kernel off. Include calibration costs in
the total efficiency report and amortize planned pilot work over the next
known proposal block when comparing candidate allocations. A group-limited
fit may need no additional subject HMC. Measure this allocation outcome and
check the resulting group behaviour.

Use ordinary Bernoulli selection for additional attempts, plus the explicitly
costed measurement at each rebuild interval. A schedule taking a step every
floor(1/r) iterations has the wrong frequency for general r and does not have
the same mixture/autocorrelation model. The mandatory measurement is a declared
experimental-design cost, not a replacement frequency rule for r_s.

### 5.3 What validity is claimed

Each fixed snapshot defines an HMC kernel targeting the exact conditional;
finite differences do not change the target. The whole production run is
adaptive because its snapshots depend on chain history. A doubling schedule
makes adaptation increasingly rare; it does not by itself prove diminishing
kernel changes, containment, or finite-sample stationarity.

Use an explicit convergence argument for the adaptive sampler. The
[AIR MCMC framework](https://arxiv.org/abs/1801.09309) provides a route through
increasingly rare adaptation under simultaneous or local drift assumptions;
[Roberts and Rosenthal](https://www2.stat.duke.edu/~scs/Courses/Stat376/Papers/AdaptiveMC/RobRosCouplingErgod2007.pdf)
provide another through diminishing adaptation and containment. Establish
which result and assumptions apply to the supported full sampler, including
its existing proposal refreshes. Bound admissible tuning and pilot work,
retain the particle kernel, and state what those restrictions establish.
Neither bounds nor toy invariance checks alone prove containment.

Produce this argument with the prototype. If assumptions remain unproved,
state the limitation explicitly and keep the adaptive method experimental;
do not resolve it by silently freezing production rescue. A fixed-snapshot
mode is useful as a validation control, while the intended implementation
adapts throughout production on the shared schedule.

## 6. Evaluation panel and decision rules

| Fixture | Purpose | Required comparison |
|---|---|---|
| Identified free-A LBA or RDM, to build | Curved subject-conditional trapping | Baseline must exhibit the diagnosed failure; HMC must improve escape and fit efficiency |
| lba71sim, saved subject 11 / chain 1 trap | Known tail-state replay | Fixed-group escape, then full-hierarchy production rescue |
| Production late-trap replay | Reactivation after initially good mixing | A zero-rate subject is reconsidered at a later doubling rebuild and rescued |
| T30 correlated Wald RDM, 40 subjects, 30 trials/cell | Group coupling and correlation diagnostics | No convergence or minimum ESS/CPU regression; record actual HMC allocation |
| corr15 | Neutral simulated cost control | No material convergence or minimum ESS/CPU regression; low sustained HMC cost |
| forstmann LNR under `fit()` | Neutral real-data cost control | Same cost and convergence assessment on real data |

`gtest/fx_T30_s1.0.rds` and the old `gtest/screen.R` are absent in this checkout.
Before benchmarking, recover their generator and specification or build a
versioned equivalent; do not claim to have run the original T30 fixture from
that path. Existing starting points include
`benchmarks/mode-start/simulate_fixtures.R`,
`benchmarks/mode-start/plausible_fixtures.R`, and
`benchmarks/group-move/fixture_fits.R`; audit their dependencies too. Identify
the exact lba71sim checkpoint and confirm its alpha, likelihood, group state,
RNG streams and tuning state. A newly generated lba71sim dataset is not a
replacement for the saved stuck-state replay.

Before implementation, verify free-A identification at each subject's own
multi-start likelihood MLE: the negative Hessian of the log likelihood
(observed information) must be positive definite and numerically stable.
Checking the Hessian at simulation truth is insufficient. Check simulation
plausibility and correct subject/parameter alignment. Establish the baseline
failure on 2-3 preregistered fit seeds using dispersed starts, exact likelihood
checks, per-subject chain diagnostics and conditional-width comparisons.
Revise an unsuitable fixture transparently, then lock its generator and seeds
before evaluating HMC. Do not search fixtures merely to favour the method.

Run these arms with identical data, starts, group-move configuration, library
build and core budget:

1. Existing default, without subject HMC.
2. A calibrated fixed-rate/fixed-tuning HMC control on the conditional escape
   and curved-toy tests, isolating the kernel from allocation failure.
3. Adaptive HMC using the shared rebuild schedule, for the whole panel.

Use default `fit()` settings for the full-fit panel. Fixed-group, imposed-tail
and fixed-snapshot checks are separately labelled diagnostic experiments. For
the late-trap replay first establish good production movement and r = 0, then
resume a labelled test checkpoint with one subject in the documented tail;
verify exact likelihood bookkeeping. Do not treat this imposed perturbation
as a production posterior sample. Also seek naturally occurring late traps
across held-out seeds.

Two seeds per arm are an initial screen. Confirm apparent gains and ambiguous
regressions with held-out seeds and an independently generated identified
fixture before choosing a default. Record seed selection and stopping rules
before confirmation runs; never discard an inconvenient chain or seed.

Measure parent and worker CPU together using cgroup `cpu.stat` around each
run. Use an isolated `systemd-run --user --scope` or equivalent cgroup, verify
its process coverage and permissions, and wait for actual completion/status
files. Parent `proc.time()` and an immediate scope-launch return are not
completed-fit CPU measurements. Report wall time and pool imbalance as well.
Run comparisons only on an otherwise idle machine; avoid oversubscription if
running two arms concurrently, and preserve immutable scripts during runs.

Report:

- `fit()`'s own stopping outcome, including its alpha/mu R-hat < 1.1 rule or
  exhaustion of tries. Keep this distinct from external diagnostic criteria.
- `posterior::ess_bulk`, tail ESS and rank-normalized R-hat for alpha, mu,
  log group SDs and free correlations. For an inference-quality conclusion,
  seek R-hat <= 1.01 and adequate ESS/Monte Carlo error on every reported
  quantity; a default-fit success alone is insufficient.
- Minimum ESS per total CPU, including calibration, rebuild and warm-up
  costs, plus production efficiency separately. Do not interpret an ESS/CPU
  gain from an unconverged arm as a successful sampler.
- Truth recovery with uncertainty on simulated fixtures, and predictive
  agreement with a trusted reference. Separate Monte Carlo failure from
  weak identification or model misspecification.
- Escape-time distributions, not just the fastest replay; r_s and eps_s
  histories, selected/attempted HMC counts, predicted versus measured mixing,
  group behaviour, failed-fit coverage and likelihood-row totals.

Proceed only if identified subject traps are rescued consistently, including
after production begins, with improved minimum ESS/CPU; neutral fixtures have
no material cost or diagnostic regression; and T30 does not regress. Define
cost uncertainty from paired repeats before judging an increase to be within
noise. If fixed HMC works but adaptive HMC does not, fix calibration or
allocation rather than changing the target or tuning a fixture-specific rate.

## 7. Implementation sequence and deliverables

1. **Lock the baseline and fixtures.** Record commit, dirty patch, installed
   library, options, seeds and data/checkpoint hashes. The existing uncommitted
   shear edits are a separate experiment; compare HMC arms against the same
   chosen baseline. Recover T30 and the stuck checkpoint, build and screen the
   identified free-A fixture, and reproduce the conditional failure. Deliver
   runnable generators/replays and a fixture manifest.
2. **Build the scratch kernel.** Implement the support predicate, whitening,
   deterministic force, reversible rejection rules, exact acceptance and row
   accounting outside the package. Test Gaussian force/sign/scale and O(eps^2)
   energy error at fixed integration time, reversibility, curved-target
   invariance from exact draws, and support-boundary rejection. Include a
   deliberately broken control to establish test power. Deliver the kernel
   and a correctness report before adding a rate controller.
3. **Validate escape and calibration.** Replay the subject 11 conditional
   from dispersed positions, compare against the particle step and fixed HMC,
   verify h sensitivity, and compare ESJD tuning against an eps grid. Test
   calibration sizes and early failure recovery. Deliver calibrated costs,
   escape distributions and supported numerical bounds.
4. **Add measurement and allocation.** Implement shared final-sweep group
   recording including correlations, per-kernel jump records, pilot accounting
   and the breakpoint solver. Test rho_P = 0/1, rho_H = 0/1, rho_H >= rho_P,
   infinite/unaddressable bottlenecks, group-limited cases, ties, missing fits,
   sparse evidence, measurement-only allocations, B = 1, and heterogeneous
   costs. Compare the solver with a dense rate-grid oracle on small synthetic
   problems. Deliver a validated
   allocator and the adaptive-convergence argument.
5. **Exercise production adaptation.** Reuse the existing doubling rebuild
   clock for all HMC snapshots. Test the scheduled-attempt minimum, r = 0
   reactivation, formerly ineffective HMC recalibration, late-trap rescue,
   changing group state, proposal trimming,
   and checkpoint continuation across a rebuild. Run the initial whole panel
   and then held-out confirmation. Deliver an auditable decision report using
   section 6's criteria.
6. **Port only after the gates pass.** Add `R/subject_hmc.R` for plan state,
   calibration and kernel helpers. Integrate rebuilds in `R/fitting.R`, final
   group recording and serial subject tasks in `R/sampling.R`, and the same
   particle-plus-HMC task in `.emc_wpool_compute_particle()` /
   `.emc_wpool_iter()` in `R/chain_pool.R`. Extend `R/profile_schema.R` and
   checkpoint metadata with counts, plan versions and calibration status.
   Keep transport and result assembly consistent across all routes.
7. **Complete package validation and handoff.** Add meaningful kernel,
   allocation, supported/unsupported-path, serial/static/queued identity,
   worker fallback, trimming and exact-resume checks; run the relevant tests
   and full suite in an isolated install. An old checkpoint without HMC
   metadata resumes with HMC disabled until an explicitly recorded rebuild
   initializes it; distinguish that from exact continuation of a new-format
   checkpoint. Document the selected behaviour and experimental/default
   status. Prepare changes for review; commit or install into the user's
   active library only when requested.

## 8. Remaining research decisions

The conditional HMC construction, one measurement attempt per rebuild interval,
Bernoulli selection for additional attempts, worker-local execution and shared
production rebuild schedule are decided. The following need prototype evidence
before a package default is selected:

- Whether finite-difference forces remain efficient across the supported
  model sizes and boundary geometries, and which h/eps bounds the evidence
  supports.
- Whether the scalar ESJD/AR(1) allocation predicts useful whole-fit gains.
  If it fails, revise its mixing or cost model and validate the revision;
  do not add arbitrary subject thresholds to mask the discrepancy.
- Whether a simultaneous-drift or other explicit adaptive-MCMC argument
  covers the full supported sampler and its growing-block schedule.
- Whether residual free-A slowness is group coupling requiring a separate
  group-level change. Evaluate that after the subject mechanism is isolated;
  do not bundle a joint HMC or shear redesign into this implementation.

No acceptance-rate target governs the final eps or rate allocation. The
quarter-period trajectory, jitter width, finite-difference scale and pilot
batch size are explicit choices requiring sensitivity evidence. Minimum ESS
per likelihood cost is the allocation surrogate; convergence and minimum ESS
per measured total CPU determine whether this plan earns a package default.
