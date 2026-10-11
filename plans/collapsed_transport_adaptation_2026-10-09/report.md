# Exact-correction calibration for collapsed transport

> **Current status (2026-10-10):** This file is a chronological experiment
> record, not a description of the current sampler. The scalar and
> coordinate-specific exact-pilot calibration explored below was reverted.
> Production keeps the original covariance-metric scale and has no exact-pilot
> scale search. The retained changes are the flat-curvature surrogate
> safeguard, cached sanitized fits, chain-local collapsed-round allocation,
> and conditional transport for subject coefficients with hard bounds and
> `t0`. A zero-round fractional allocation skips collapsed setup. The LBA ESS
> bottleneck remains open. See “Package decision” and the follow-up sections
> for the current evidence.

## Initial proposal: exact-correction calibration (reverted)

The initial experimental version tuned the collapsed integrator against the
Gaussian surrogate. High inner acceptance could coexist with essentially zero
exact acceptance when transported subjects crossed likelihood support.
Allocation could then suppress the move without learning a usable distance.

That version used `scale * C`, where `C` is the existing learned covariance
square root and `scale` is dimensionless. At proposal rebuilds, discarded
pilots from the current state evaluated the real likelihood and the existing
exact Metropolis correction. The objective was

`E[exact_acceptance * ||solve(C, proposed_phi - current_phi)||^2] / cost`.

The intended objective balanced distance, acceptance and reserved core time.
It did not target a universal acceptance percentage or prescribe a distance
in parameter units. The exact correction included the original prior and
transport terms; it added no likelihood clipping, support relaxation or
proposal repair.

The search bracketed the objective in log distance, expanding or shrinking by two,
then refining adjacent intervals. Eight pilots per candidate and at most ten
candidates bound calibration effort. An approximate one-sided 95% Student-t
lower estimate of expected squared jump discounts noisy candidates; this is
a selection heuristic, not a finite-sample guarantee. Sixteen independent
validation pilots estimate the selected kernel's movement and cost for
allocation, avoiding direct reuse of the winning selection score.

Warmup applied the selected scale directly. Production interpolated log scale
with gain `k^(-2/3)` at calibration number `k`; the counter persists through
proposal rebuilds and checkpoints. The first update is unrestricted, allowing
old checkpoints to recover. Scale is fixed between rebuilds. Pilots never
replace retained states, and consume the chain's Gibbs RNG stream. The caller's
RNG is preserved by the allocator. Worker results return metadata and seed,
without replacing process-local sampler pointers.

Rebuilding a surrogate reset obsolete collapsed rejection records and requested
fresh calibration even if the preceding allocation was zero. Allocation can
still choose zero when the calibrated move offers insufficient benefit for its
cost. Subsequent rebuilds reconsider it. Calibration diagnostics and scale are
saved in allocation history.

## Validation

- Collapsed-move regression and fixed-kernel invariance tests passed.
- New adaptation tests passed, including hard support, diagonal-gamma and
  diagonal-HW hierarchies, unit changes from 0.001 to 1000, zero-allocation
  recovery, underflow, uncertainty scoring, unchanged retained draws and
  checkpointed diminishing updates.
- Allocation and ordinary sample proposal schedule tests passed. The optional
  full model validation in the schedule suite was skipped; see the saved log.
- A worker-pool/allocator integration probe used two cores, preserved retained
  draws, advanced the chain stream and recovered 65% mean exact acceptance.
  Allocation nevertheless chose zero rounds on that probe. A separate forced
  single transition had stored likelihoods matching fresh evaluations exactly.

Independent checks on the frozen LBA checkpoint used 20 new proposals per
chain at the scale learned by calibration:

| Chain | Learned scale | Mean exact acceptance probability |
|---|---:|---:|
| 1 | 0.044194 | 0.328 |
| 2 | 0.062500 | 0.343 |
| 3 | 0.088388 | 0.123 |

Previously the unscaled transport had effectively zero exact acceptance.
Calibration took about 8–10 seconds per chain with one core on this snapshot.
These are local proposal checks, not a full-chain ESS or convergence benchmark.
The final diminishing-update addition leaves these first calibration updates
unchanged. Scripts and console results are alongside this report; reproducing
the checkpoint probes requires the frozen RDS files under
`/tmp/lba_sampling_audit/`, which are not copied into the repository.

## Limits

The search uses finite local evidence and is not guaranteed to find a global
efficiency optimum. Its pilot budgets, confidence level and search resolution
control numerical effort and uncertainty, rather than encode LBA parameter
scales. A scalar multiplier respects the learned covariance geometry but does
not fix every directional defect in a surrogate. Actual sampling efficiency
still needs assessment over full runs.

Diminishing this scale update alone is not a proof of ergodicity for the
sampler's entire existing adaptive proposal and allocation system. For any
fixed plan the original exact acceptance correction remains in use. Neither
this implementation nor its checks establish that the LBA ridge is resolved.

Changes are in the source checkout only. The running fit, checkpoint and
installed package have not been modified.

## ESS and core-utilization check

I compared the previous sampler at `HEAD` with the adapted working tree using
the same frozen three-chain LBA checkpoint. Each arm used the same data,
proposal settings, 200 retained draws per chain and 12 total cores (three chains
with four likelihood cores each); subject HMC was disabled in both arms to
isolate the collapsed move. The baseline's first 200 draws reproduced exactly
in a fresh control run. Calibration pilots consume the chain RNG, so subsequent
adapted transitions are not common-random-number paired. The adapted arm was
then continued to 400 draws per chain and compared with the saved 400-draw
baseline.

| Bulk ESS summary | Baseline, 200 | Adapted, 200 | Change | Baseline, 400 | Adapted, 400 | Change |
|---|---:|---:|---:|---:|---:|---:|
| Median across group means | 44.8 | 39.4 | -12% | 124.1 | 65.1 | -48% |
| Median across participant 52091 parameters | 75.2 | 109.8 | +46% | 95.3 | 53.8 | -44% |

For the matched 200-draw runs, elapsed time including rebuild and calibration
was 323.6 seconds for baseline and 376.9 seconds for adapted. This makes the
adapted arm's median ESS per second about 24% lower for group means and 25%
higher for participant 52091. The baseline 400-draw elapsed time was not
recorded, so the 400-draw ESS comparison is not cost-adjusted. At 400 draws,
individual parameters remain mixed: for example, participant 52091's `A`,
`B`, and `t0_LoadHigh` ESS improved, while the group means and other participant
parameters generally declined. Several adapted R-hats remain above 1.1. This
short experiment does not establish a general ESS improvement or adequate
convergence.

The adapted kernel selected scale `0.044194` and allocated four collapsed
rounds per transition. Production exact acceptance was 0.431, 0.457 and 0.393
across the three chains, confirming the bound/support failure was recoverable.
The observed ESS cost shows that recovering acceptance alone does not ensure
better end-to-end mixing; allocation and the LBA ridge need further study.

The earlier 8–10 second-per-chain timing above was a deliberately single-core
diagnostic. In this 12-core run, allocation gave each of three concurrent
chain calibrations four cores. Their recorded pilot times were 5.0, 7.3 and
5.1 seconds; rebuild wall time was 16.5 seconds. The calibration starts a
likelihood worker pool using its assigned core budget. During the ESS run, the
12 likelihood workers were active on distinct CPUs in the process affinity
set (0–31). The separate subject-HMC pilot path parallelizes subject jobs up
to the supplied core count and uses one likelihood core per job, avoiding
nested oversubscription. Thus the single-core note described that diagnostic
configuration, not an intentional limit in the production allocator.

Compact ESS summaries, run logs and scripts are stored beside this report.
The 400-draw baseline elapsed time and a longer replicated ESS comparison
remain unavailable, so these results should guide the next benchmark rather
than serve as a performance guarantee.

## Follow-up: chain synchronization interaction

The earlier ESS comparison exposed a concrete interaction that was not
accounted for in its interpretation. The allocator solves each chain's move
budget from that chain's state, pilot movement and measured cost. At the end
of rebuilding, however, the sampler replaced every active chain's collapsed
rate with the largest rate selected by any chain. On the frozen checkpoint,
local choices included `1.125`, `4`, and `0`; synchronizing to the maximum made
all chains execute four collapsed rounds per transition. A chain whose own
allocation rejected the move could therefore be made to spend substantial
time on it. Chains are independent Markov chains, so equal collapsed rates
are not required for target invariance.

I removed that cross-chain maximum and added a regression test asserting that
the allocator's distinct chain rates survive proposal rebuilding. The
`move-allocation` test file passed (including the new regression test).

To isolate the effect on ESS, I calibrated once and cloned the same
post-calibration state, scale and per-chain RNG state into two 200-draw arms.
The only policy difference was whether chains retained their local rates or
all used the largest local rate. The local choices in this calibration were
`1.125`, `0`, and `0`, with scale `0.044194` for each chain; the shared arm
used `1.125` in all three chains.

| Summary across parameters | Chain-local rates | Shared maximum | Difference |
|---|---:|---:|---:|
| Median bulk ESS, group means | 79.5 | 59.1 | +34% |
| Median bulk ESS, participant 52091 | 140.2 | 111.0 | +26% |
| Sampling elapsed time | 381.5 s | 537.9 s | -29% |

The local arm had exact collapsed acceptance `0.393` in chain 1 and no
collapsed attempts in chains 2–3. The shared arm's exact acceptance was
`0.393`, `0.464`, and `0.445`. Thus acceptance remained healthy in chains
where the move was already selected, while forcing it on the other chains
increased runtime and coincided with lower pooled ESS in this paired run.
The elapsed times include sequential arms under a busy host, so the runtime
difference is suggestive rather than a clean throughput estimate. Although
both arms started from identical states and RNG seeds, their trajectories
diverged as soon as the move policies differed.

This supports cross-chain rate synchronization as a real contributor to the
observed slowdown and ESS regression. It does not prove that extra collapsed
moves lower expected ESS: this is one 200-draw paired comparison, and
finite-run ESS is variable. The earlier 400-draw chain-local run also had
higher median ESS than the four-round shared run (77.1 vs 65.1 for group
means; 62.7 vs 53.8 for participant 52091), but remained below the 400-draw
baseline (124.1 and 95.3). These mixed, parameter-specific outcomes still
point to the difficult LBA ridge and need longer replicated comparisons.

Calibration/allocation choices also varied between rebuilds on this
checkpoint: observed local rates included `(1.125, 4, 0)`, `(1.75, 4, 0)`,
and `(1.125, 0, 0)`. The learned scale stayed near `0.044`, while measured
pilot costs and resulting move budgets changed. That variability is a separate
reason to avoid treating a single ESS run as decisive; it merits checking the
stability of the cost estimate under host load in the next benchmark.

## Follow-up, 2026-10-10: saturation, setup cost and incremental allocation

Reading the code and profiling the frozen checkpoint exposed three additional
problems beyond the already removed cross-chain maximum. None requires
changing the likelihood, target distribution or exact acceptance correction.

First, every surrogate rebuild resets the collapsed move records. Exact
calibration then populated the allocation record using sixteen independent
single-round proposals from the same retained state, assigning `N = E` and
putting all counts in column one of `K`. Consequently the conditional
saturation estimator always returned `a = g`, irrespective of the preceding
production block's evidence. Repeating this at every rebuild permanently
restored the optimistic assumption that repeated collapsed rounds compound
geometrically. Keeping old records would mix evidence from different kernels,
so the fix measures consecutive rounds under the newly calibrated kernel.
Validation now runs sixteen independent private batches of two rounds. Each
batch starts at the retained state, carries accepted states into its second
round, records the actual net jump in `N`, and contributes one count in column
two of `K`. Selection remains eight pilots per candidate; validation therefore
costs thirty-two rounds rather than sixteen. No pilot state is retained.

Second, calibration priced only its round loop. Production first reconstructs
the collapsed context and the current subjects' conditional frames on every
iteration using the move. A twenty-call, single-process profile of the frozen
LBA state measured approximately 3.1 ms for context construction, 20.0 ms for
the current frame state, 1.6 ms for the inner integrator, and 137.9 ms for an
entire one-round move including serial likelihood evaluation. These are local
measurements, not portable throughput estimates, but show that the integrator
is a small part of this workload. Validation now separately measures context
and frame setup. Allocation charges this once on each iteration running at
least one round: `min(1, expected_rounds) * fixed_cost`, in addition to the
per-round cost. Fractional allocations also return before constructing the
context when their Bernoulli draw selects zero rounds; counters and the chain
RNG advance as before. The fixed cost and the objective without collapsed
moves are saved in allocation history.

Third, the existing ten-percent minimum gain applied only to the combined
optional-move plan against an iteration with none of those moves. Effective
shears could pass this check and thereby let collapsed rounds be purchased
for a negligible extra predicted gain. On the fresh paired-round calibration,
chain one's objective was 185.737 with collapsed rounds versus 188.572 with
the other moves optimized alone: only a 1.50% improvement. Both comfortably
beat the no-optional-move objective of 350.287. The allocator selected 1.25
collapsed rounds despite their small incremental benefit. It now applies the
same existing ten-percent threshold to the best collapsed allocation against
the optimized allocation without collapsed moves. This does not prescribe a
new acceptance target or transport distance.

The same calibration found single-round factors of approximately 0.993–0.998
in the slow group-mean coordinates. Thus restored acceptance translates to
only small movement relative to the posterior spread. Two-round validation
identified saturation in 8, 30 and 35 of the 36 group coordinates across the
three chains. Setup and saturation corrections alone left the selected rates
at `(1.25, 0, 0)` in this probe, with small changes in shear rates; their
presence is not evidence of an ESS improvement. The incremental check is
needed to distinguish a useful combined plan from a useful collapsed move.

The lag-one correlation allocation model is still an approximation to ESS.
Local jump statistics cannot establish movement along the LBA ridge, and the
conditional saturation model does not establish independence between shear
and collapsed movement. The paired-round saturation estimate is also noisy
with sixteen starts. Short continuation benchmarks must therefore report
ESS, elapsed time and convergence diagnostics without interpreting healthy
acceptance as a performance guarantee.

A cost sensitivity probe multiplied both calibrated marginal and setup costs
by `0.25`, `0.5`, `1`, `2`, and `4`. The incremental-gain policy selected
`(0, 0, 0)` collapsed rates in every case. On this state the rejection of
extra collapsed work is therefore stable over substantial timing variation;
this probe does not vary movement estimates or establish stability on other
states. Results are in `cost_sensitivity.txt`.

Regression checks passed for collapsed adaptation, move allocation, the
collapsed kernel's fixed-plan invariance, and the production proposal
schedule. The schedule suite skipped its opt-in full model validation
(`EMC2_TEST_LEVEL=full`). Added regressions check net cancellation of accepted
private moves, fixed setup charges above and below one expected round,
skipping setup on zero-round fractional draws, rejecting a negligible joint
improvement over an effective shear plan, and retaining a substantial joint
improvement. The worker-pool integration check passed, preserved retained
draws, advanced the chain stream, and found zero discrepancy between stored
and freshly evaluated likelihoods after a forced retained move. Logs are
`setup_saturation_tests.txt` and `setup_saturation_integration.txt`.

A separate zero-round profile compared the old and new function directly on
the frozen chain, forcing the fractional allocation to draw zero. The outputs,
including RNG state and counters, were identical. Over 100 calls it measured
4.75 ms per call before versus 0.69 ms after. This is evidence for the removed
setup work, not an estimate of complete-chain speedup. The script and output
are `profile_zero_rounds.R` and `zero_round_profile.txt`; its previous-source
backup is in `/tmp/collapsed_followup_before/`.

When HMC is available, the incremental comparison also considers the best
allocation with and without its fixed round cost. A collapsed plan must beat
the best policy without collapsed moves, even when that alternative uses a
different HMC policy. A regression covers a collapsed plan that wins against
its own HMC policy but gives less than ten-percent improvement against the
best alternative. The final allocation test file passed 706 assertions;
see `incremental_budget_tests.txt`.

### Matched continuation assessment

The continuation comparison isolates the incremental-gain check. Both arms
already use the setup and saturation corrections described above. They start
at the same calibrated retained states, scales and chain RNG states; the
reference uses the former combined-plan threshold and rates `(1.25, 0, 0)`,
while the incremental check selects `(0, 0, 0)` and reoptimizes the shear
allocation. Scales are `(0.044194, 0.044194, 0.088388)`. HMC is disabled to
isolate this decision. Each arm has three chains, four likelihood cores per
chain and 200 new retained draws per chain. The second repetition advances
all chain streams by three full L'Ecuyer streams from the same calibrated
starting state, and reverses arm order. Proposals and allocations stay fixed
through these continuations. This is a comparison of allocation policies at
one checkpoint, rather than a new comparison with the pre-adaptation `HEAD`
sampler.

The initial proposal rebuild and calibration took 15.794 seconds, shared by
both policies. The table gives sampling wall time separately; adding that
common cost to each arm gives elapsed time including the initial rebuild.
Timing excludes writing the benchmark artifacts.

| Repeat | Allocation | Group mean ESS | Group log sd ESS | All participant ESS | Participant 52091 ESS | Sample seconds |
|---|---|---:|---:|---:|---:|---:|
| 1 | Former threshold | 28.0 | 207.0 | 294.9 | 87.4 | 351.7 |
| 1 | Incremental check | 64.0 | 185.4 | 320.1 | 156.4 | 358.0 |
| 2 | Former threshold | 39.4 | 212.5 | 264.1 | 126.1 | 514.1 |
| 2 | Incremental check | 54.8 | 176.3 | 267.2 | 125.9 | 395.1 |

ESS values are medians of rank-normalized bulk ESS across 18 group means,
18 group log standard deviations, all 1,674 participant parameters, and the
18 parameters of participant 52091, respectively. Removing the marginal
collapsed allocation increased median group mean ESS by 128% and 39% in the
two repetitions. Across all participants the increases were only 9% and 1%;
the much larger improvement for participant 52091 in the first repetition
did not recur in the second. Median group log sd ESS declined by 10% and 17%.
This is consistent with an allocation focused on the slower group means
sacrificing some movement in the already faster variance coordinates. It is
not evidence that the scaled transport cannot benefit those coordinates.

Sampling time increased by 2% in the first pair and declined by 23% in the
second. Group mean ESS per second including the common rebuild cost increased
by about 124% and 79%. The sequential runs used a busy host; reversing order
helps expose timing variation but does not control it. These two short
repetitions support rejecting the marginal collapsed work on this checkpoint,
with mixed consequences for individual parameter families. They do not
establish a general speedup or an ESS guarantee. Exact collapsed acceptance
in the reference's active chain was 0.400 and 0.438, again demonstrating that
healthy acceptance does not by itself justify that allocation.

Convergence remains incomplete: maximum R-hat across the 36 group coordinates
was 1.160 and 1.144 under the reference, and 1.116 and 1.147 under the new
check. Of 1,674 participant coordinates, 920 and 839 had R-hat above 1.01 in
the reference versus 706 and 753 with the new check. These continuations are
sampler diagnostics; they do not support substantive posterior conclusions
or establish that the LBA ridge has been resolved.

The scripts are `compare_setup_saturation.R`, `previous_budget_solver.R`,
`compare_incremental_budget.R` and `summarise_incremental_budget.R`.
The latter two produce `incremental_budget_ess.rds`,
`incremental_budget_summary.csv`, and four CSV files with bulk ESS, tail ESS
and R-hat for every sampled group and participant coordinate. Console logs
are beside them. The initial calibration checkpoint and full continuation
fits are saved in `/tmp/collapsed_transport_followup/`; the original frozen
checkpoint remains under `/tmp/lba_sampling_audit/`. The initial diagnostic
script was stopped after its first complete arm when the incremental-gain
issue was identified; that saved arm is reused as the first reference in the
final paired comparison. To reproduce the former allocation with current
source, `previous_budget_solver.R` removes the incremental check from local
copies of the solvers without modifying the package namespace.

Changes remain in the source checkout. The original frozen checkpoint and
installed package have not been modified.

## Follow-up, 2026-10-10: shared scale, t0 support and ESS

The shared multiplier is a real source of cross-parameter coupling. The
proposal first builds `C` from the covariance of the unconstrained group chart
across the proposal window, then the C++ Hamiltonian trajectory uses
`phi = C x`. `C` carries the chart coordinates' different units, marginal
scales and correlations. The additional scalar `scale` multiplies every row
of `C`, so a support limit in one direction can force every chart direction
to use a smaller fraction of its learned scale. Raw `t0` units do not directly
set the other parameters' scales; the common multiplier does.

In the current proposal the effective factor matrix is `D C`: each non-t0
chart row in `C` is multiplied by the global scale `s`; t0 mean rows are
multiplied by `s * r_mean`; and t0 variance/covariance rows by
`s * r_hyper`. Thus at `s = 0.03125`, `r_mean = 0.5`, and
`r_hyper = 0.25`, non-t0 chart directions use 3.125% of their learned row
scale, t0 means use 1.5625%, and t0 hierarchy directions use 0.78125%. The
learned rows of `C` still encode each chart coordinate's own scale and
correlation. These are proposal displacement factors, not the leapfrog step
size `eps`: the C++ position update is `x += eps * C * v`, so shrinking `C`
shrinks movement in every non-t0 direction as well as t0. `eps` controls
integration resolution and trajectory length controls total integration time
(`ceil(length / eps)` leapfrog steps). In these runs `eps` was at its cap of
1.5, so step-size tuning could not undo a small `s`; a longer path changes
integration time and is scored separately.

For this LBA checkpoint, the subject-level model bound is `t0 >= 0.05`, and a
finite observed response also requires `t0 < RT`. The model maps each `t0_*`
coefficient through `t0 = exp(alpha)`, so these subject-level constraints are
`alpha >= log(0.05)` and `alpha < log(RT)`. The group chart's `t0_*` means are
on this log-t0 scale and are not themselves clipped at those limits. A
proposed group change is transported to all subjects, and the exact likelihood
correction evaluates the resulting subject parameters together. One subject
crossing support can therefore reject the whole group move. The bounds probe
found, averaged over the saved proposals and chains:

| Global scale | Subjects with proposed t0 below 0.05 | Subjects with t0 at/above an observed RT |
|---:|---:|---:|
| 1.00 | 10.75 | 8.25 |
| 0.30 | 3.33 | 0.00 |
| 0.10 | 0.42 | 0.00 |
| 0.03 | 0.00 | 0.00 |

The parameter-by-parameter audit corrected the earlier conclusion that other
bounded parameters also contributed: every explicit model-bound crossing in
these probes was `t0 < 0.05`; no proposal crossed the bounds for `A`, `B`,
`sv`, `k`, `mG`, `mK`, `eta`, `pContaminant`, `pGuess`, or `v`. The second
column is a separate likelihood-support failure, `t0 >= RT` for at least one
observed response. Counts are averages per proposal over the three chains.
Both limits apply to subject values, while the transported chart coordinates
are group means and variance/covariance terms. One subject outside support
rejects the entire joint group move.

The prior implementation used one global multiplier on every row of `C`,
then a second factor on t0 means only. On the earlier chain-1 plan, the three
t0 mean row norms of `C` were 0.070–0.081 on the `log(t0 in seconds)` scale.
At global scale 0.0625 and mean factor 0.5, their initial velocity scales
were about 0.0022–0.0025 log-scale units, while every other chart row received
the same 0.0625 factor. These are chart velocities, not final proposal jumps.

The source keeps the learned covariance factor and exact Metropolis
correction, and applies separate effective scales to t0 means and
variance/covariance chart directions. In diagonal-gamma hierarchies this
includes the log-variance coordinates; in inverse-Wishart blocks it includes
the Cholesky entries involving t0 as a correlated block. The ratios are
stored relative to the non-t0 scale, but calibration tunes the t0 effective
scales first and then holds them fixed while searching the non-t0 scale. The
first joint search is only a bootstrap; later rebuilds use the independent
parameterization. This removes t0 support as a reason to shrink the
non-t0 multiplier. A single scale still ties all non-t0 chart rows together,
and the exact correction can penalize larger joint moves for other reasons.

The calibration also previously compared different metric scales using the
same Hamiltonian step size, even though changing `scale` changes the
surrogate Hamiltonian's curvature. It now dual-averages a step size on the
surrogate for each candidate metric and path length, retains the selected
step size, and includes a half-scale candidate in the t0 joint search. It
also searches trajectory lengths using exact accepted squared movement per
measured cost. The surrogate step-size tuning itself does not use likelihood
calls. On fresh calibration of the frozen checkpoint, all three chains
selected the maximum allowed step size (1.5), so the cap is active for these
small metrics.

Before separating t0 means from their variance/covariance coordinates, fresh
full-grid calibration with candidate-specific step size and path-length search
gave:

| Chain | Global scale | t0 factor | Effective t0 scale | Path length | Step size | Exact validation acceptance | Calibration seconds | Exact pilot proposals |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.03125 | 0.25 | 0.00781 | 3.00 | 1.50 | 0.443 | 12.63 | 240 |
| 2 | 0.03125 | 0.25 | 0.00781 | 3.00 | 1.50 | 0.579 | 12.03 | 248 |
| 3 | 0.04419 | 1.00 | 0.04419 | 0.75 | 1.50 | 0.582 | 13.73 | 256 |

A 24-proposal length probe at the previously calibrated scales compared
lengths 1.5, 3 and 6. Accepted squared chart movement per second was 3.01,
8.22 and 6.21 on chain 1; 3.44, 8.23 and 9.48 on chain 2; and 6.38, 7.36
and 1.30 on chain 3. Exact acceptance at length 6 fell to 0.144, 0.176 and
0.017, respectively. This supported trying longer paths than 1.5, while also
showing that long paths can be wasteful for some chains. The final calibration
selected lengths 3, 3 and 0.75, respectively. These are short local
diagnostics, not ESS measurements.

The global multiplier remains small even after t0 directions are separated.
The search is sequential over global/t0 scale and path length, and its local
score does not optimize long-run ESS or guarantee movement along the slow LBA
ridge. The earlier means-only full t0 ratio search used roughly
70–77% more exact pilot proposals per chain than the scalar-only search (about
176–192 versus 104–112). Its three-chain proposal rebuild took 24.2 seconds;
an earlier scalar-only rebuild took 15.8 seconds. Rebuild timing is
host-sensitive, but the extra pilot count is real and must be repaid by better
sampling efficiency.

To reduce repeated cost when that full search was followed by zero allocated
rounds, later same-stage rebuilds now screen a smaller neighborhood around the
last calibrated scale and t0 factor. At the time of the experiments below,
stage changes still used the full grid. A follow-up now also screens the first
sample rebuild when a burn/adapt stage already produced a calibrated plan;
un-calibrated starts retain the full grid. The earlier
means-only, fixed-step screen on chain 1 evaluated 112 exact pilots, found
validation acceptance 0.556, and took 8.1 seconds. The next screen added t0
variance directions and candidate-specific step tuning but predated path
length search; it also used 112 exact pilots, selected scale 0.0404 and t0
factor 0.323, reached validation acceptance 0.527, and took 8.45 seconds on
four cores. These screens are local neighborhoods with independent
validation, not full searches for a new optimum.

The path-length-aware screen on chain 1 used 136 exact pilot proposals,
selected scale 0.03125, t0 factor 0.646, path length 0.969 and step size
1.493, with validation acceptance 0.668. It took 8.96 seconds on four cores.
That is 24 more pilots and about 6% more elapsed time than the previous
112-pilot screen on this host; the small single-chain timing difference is
not a stable runtime estimate. Its CSV and reproduction script are
`length_adaptive_calibration_screen.csv` and
`benchmark_length_adaptive_calibration_screen.R`.

To test the allocated move itself, I compared a no-collapsed arm with an arm
that ran the t0-adapted move at the allocator's no-gate candidate rate, 0.875
rounds on chain 1 and zero on chains 2–3. At that earlier means-only,
fixed-step checkpoint, the allocator chose zero on all three chains: chain
1's predicted objective improved by only
3.3%, below its existing 10% marginal-gain threshold; chains 2–3 have
essentially no predicted gain. Therefore this continuation is a counterfactual
test of that marginal move, not the automatic policy's current behavior.
Both arms share the same calibration rebuild cost (24.179 seconds). Each arm
ran three chains for 200 retained draws per chain, with four likelihood cores
per chain and subject HMC disabled. The starting states and all other move
rates matched within each pair.

| Repeat | Arm | Total seconds including rebuild | Median group-mean ESS/sec | Median t0 group-mean ESS/sec | Median participant ESS/sec | Participant 52091 ESS/sec | Sum ESS/sec across 1,710 variables |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | No collapsed move | 376.3 | 0.0671 | 0.0801 | 0.6931 | 0.1237 | reference |
| 1 | t0-adapted move | 417.0 | 0.0859 (+28%) | 0.2142 (+168%) | 0.6873 (-1%) | 0.3211 (+160%) | 0.928x |
| 2 | No collapsed move | 414.4 | 0.0353 | 0.0297 | 0.3000 | 0.0666 | reference |
| 2 | t0-adapted move | 378.7 | 0.0414 (+17%) | 0.0323 (+9%) | 0.3819 (+27%) | 0.0901 (+35%) | 1.114x |

The second pair advanced the chain RNG streams and reversed arm order. These
ESS comparisons used the earlier means-only, fixed-candidate-step
implementation; they do not test the revised t0 variance scaling or
candidate-specific step tuning below. The observed time difference reversed
too: the t0 arm was 11% slower in repeat 1
and 9% faster in repeat 2, so timing is not decisive at this run length. The
sum of per-variable ESS per second also flipped from 7% worse to 11% better.
Across individual variables, 598 of 1,710 improved in repeat 1 and 1,218 did
in repeat 2. Median group-mean ESS/sec improved in both pairs, including the
non-t0 group means, but other family and individual results remain mixed.
The active chain's exact acceptance estimates were 0.294 and 0.311. Maximum
R-hat was 1.23–1.30 in these short continuations, so their ESS differences
are noisy and convergence remains incomplete.

I am keeping the 10% marginal allocation gate. These two short comparisons do
not establish that the t0-adjusted move improves overall wall-time per ESS,
and the current pilot search adds measurable rebuild work. The zero-acceptance
proposals in the original diagnostic were not useful for sampling. If
production repeatedly invoked that kernel, the likelihood calls spent on
rejected proposals were wasteful; preserving the target distribution does not
make those calls efficient. In the saved frozen baseline continuation the
allocator assigned zero collapsed rounds on all chains, so it avoided
repeating that production cost there. That distinction is why “as good” was
the wrong interpretation: the unscaled kernel remained valid as a Markov
transition, but was a poor movement proposal. Raising its acceptance repaired
that failure mode; it did not by itself establish a better overall sampler.
The row-specific calibration remains a candidate fix, with broader throughput
benefit still unproven.

The paired fits and per-variable summaries are under
`/tmp/collapsed_transport_followup/` and in
`t0_policy_allocation_family_summary.csv`,
`t0_policy_allocation_per_variable.csv`, and their `_repeat2` versions beside
this report. `benchmark_t0_policy_allocation.R`,
`benchmark_t0_policy_allocation_repeat2.R`, and
`summarise_t0_policy_allocation.R` reproduce the comparisons.

The parameter-by-parameter bound audit, revised calibration, and focused
adaptation, collapsed-move, and move-allocation tests pass. These results
confirm that t0 support drove the earlier hard-bound failures and that the
global scale ties the other chart directions to it. They do not yet show that
the current scale objective maximizes long-run ESS per second. The
path-length paired continuations below add short-run evidence, but do not
establish a general ESS-per-second improvement.
The attribution is reproducible with `bounds_by_parameter.R`; fresh and
screen calibration runs are in `benchmark_fresh_retuned_calibration.R` and
`benchmark_zero_rate_calibration_screen.R`, with their CSV/RDS outputs beside
this report.

### Path-length policy paired continuation

The current allocator's prepared plan selected four collapsed rounds only on
chain 2; chains 1 and 3 received zero. Its chain-local settings were
`scale/t0 factor/path length/step size` = `0.03125/1/1.5/1.5`,
`0.03125/0.25/3/1.496`, and `0.04419/0.5/3/1.5`. I compared this allocated
policy with the same saved starts and RNG streams but no collapsed rounds.
Both arms included the same 24.15-second calibration/rebuild cost and ran 200
retained draws per chain, with three chain workers and four likelihood cores
per chain. The second repeat advanced the RNG streams and reversed arm order.

| Repeat | Arm | Total seconds incl. rebuild | Median group-mean ESS/sec | Change | Median participant ESS/sec | Change | Max key-parameter R-hat |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | No collapsed rounds | 369.3 | 0.1340 | — | 0.2408 | — | 1.152 |
| 1 | Allocated policy | 429.4 | 0.0997 | -25.6% | 0.2233 | -7.3% | 1.128 |
| 2 | No collapsed rounds | 394.3 | 0.1530 | — | 0.3326 | — | 1.096 |
| 2 | Allocated policy | 387.4 | 0.1607 | +5.0% | 0.4078 | +22.6% | 1.086 |

The allocated chain's exact collapsed acceptance was 0.453 in repeat 1 and
0.463 in repeat 2. Its runtime was 16% slower in the first pair and 2% faster
in the second; group-mean ESS/sec declined in one pair and improved slightly
in the other. Participant ESS/sec moved in the same direction as group means
in each repeat, with larger swings. This is mixed evidence from short,
200-draw continuations. It does not justify claiming the calibrated move
improves wall-time per ESS. The current marginal allocation gate remains
important: it chose only four rounds on one of three chains, and neither
zero-acceptance proposals nor a healthy acceptance estimate alone decide the
efficiency question.

The comparison is reproduced by `benchmark_length_policy.R` and
`benchmark_length_policy_repeat2.R`; compact results are in
`length_policy_ess_comparison.rds` and
`length_policy_repeat2_ess_comparison.rds`.

### Split t0 mean/hierarchy scale counterfactual

The combined t0 factor reduced ESS/sec for t0 group means in both of its
200-draw paired repeats. The source now gives t0 group means and t0
variance/covariance chart rows distinct relative factors. On a frozen chain-2
state, a 100-proposal grid found a useful tradeoff: changing the mean factor
from 0.25 to 0.5 while holding the global scale at 0.03125 and the hierarchy
factor at 0.25 raised t0-mean accepted movement per second from 0.0170 to
0.0895 (5.3x), while total accepted movement per second changed from 3.592 to
3.506 (-2.4%). This is a local pilot score, not ESS.

I then compared a forced four-round split-factor move on chain 2 with no
collapsed rounds on any chain, using the same saved starts and RNG streams.
Both arms included the same 25.511-second rebuild; each sampled 200 draws per
chain with subject HMC disabled. The forced chain-2 plan used global scale
0.03125, t0 mean factor 0.5, t0 hierarchy factor 0.25, path length 3 and
`eps = 1.5`. Its exact acceptance was 0.426. This was a counterfactual: the
actual allocator assigned zero collapsed rounds to all three chains.

| Summary | No collapsed rounds | Forced split t0 | Change |
|---|---:|---:|---:|
| Total seconds, including rebuild | 399.7 | 393.2 | -1.6% |
| Median group-mean ESS/sec | 0.1111 | 0.1323 | +19% |
| Median t0 group-mean ESS/sec | 0.0340 | 0.0477 | +40% |
| Median non-t0 group-mean ESS/sec | 0.1220 | 0.1422 | +17% |
| Median group log-SD ESS/sec | 0.4936 | 0.4778 | -3% |
| Median participant ESS/sec | 0.6284 | 0.7049 | +12% |
| Sum of variable ESS/sec | 1,143.3 | 1,213.7 | +6% |
| Maximum R-hat | 1.339 | 1.311 | — |
| Variables with R-hat above 1.05 | 197 | 206 | — |

The allocator's zero decision is explainable from its current bottleneck
metric. In the prepared state, the slowest raw group coordinate was
`mu:v_LoadHigh` in all three chains (rho about 0.959–0.961). Its measured
one-round collapsed factors were 0.995, 0.995 and 0.999, compared with
0.956, 0.951 and 0.955 for shear. The minimax allocation model therefore
predicts that collapsed transport barely changes the current bottleneck,
while shear changes it more. The short forced run nevertheless showed
`mu:v_LoadHigh` ESS/sec up about 21% and median non-t0 group-mean ESS/sec up
about 17%. This is a concrete model-versus-run discrepancy to investigate;
with 200 draws and R-hat above 1.3, it does not yet prove the allocation model
wrong.

Participant 52091 moved the other way: its median ESS/sec was 0.1209 with no
collapsed rounds and 0.1069 with the forced move (-12%). The small elapsed-time
difference is not evidence of a speedup; the forced arm ran first, and the
same rebuild cost was charged to both arms. Several R-hats are far above 1.01
and hundreds of variables exceed 1.05, so these 200-draw ESS estimates are
noisy and neither arm establishes convergence. Taken together, the counterfactual
is mildly positive in aggregate ESS/sec but mixed by parameter family and
participant. It does not establish a material end-to-end improvement, and it
does not show that the actual allocator benefits because that allocator chose
zero rounds.

The source-level conclusion is more specific: the earlier shared scale did
couple t0 support to every other chart direction. Calibration now removes
that coupling when it refines the non-t0 scale, while one multiplier still
ties all non-t0 rows together. A short forced run with the earlier
parameterization showed mixed ESS evidence; the local probes do not establish
that the revised scale search improves long-run ESS per second. The
calibration result, controlled run, all-variable ESS table and reproduction
script are `fresh_retuned_calibration.csv`,
`split_t0_counterfactual_family_summary.csv`,
`split_t0_counterfactual_per_variable.csv`, and
`summarise_split_t0_counterfactual.R` beside this report. Full fits remain in
`/tmp/collapsed_transport_followup/`.

### Follow-up, 2026-10-10: independent evaluation of the surrogate fit

The review's diagnosis of the fitted quadratic is supported by the saved
chain-1 fits. Across the 93 subjects, the fitted precision ranks ranged from
13 to 17 in 18 dimensions, and all 93 had a nonzero fitted slope in at least
one flat direction. A linear term in a flat likelihood direction makes the
Gaussian surrogate unbounded; in the collapsed integral this can create an
exponential incentive to expand group covariance. This distorts the proposal
Hamiltonian, although the unchanged exact Metropolis correction still targets
the original posterior and prevents that surrogate error from changing the
target distribution.

The fit now projects the refitted linear term into the positive-curvature
subspace after the PSD projection. Saved and Laplace fallback fits receive
the same sanitization when needed. A generic Laplace fallback was not used:
the saved Laplace fits can also have weak or flat directions, so their label
alone does not guarantee a bounded quadratic. The sampler caches the
sanitized precision, linear term and pooled cross block, then reconditions
those base fits when pooled coefficients change. In a 20-context profile on
the same checkpoint, context construction took 0.359 seconds before caching
and 0.034 seconds after caching. A collapsed-move invariance test caught and
prevented caching fallback Laplace fits when a caller recomputes them from
changing data. The focused collapsed-adaptation, collapsed-move and
move-allocation suites passed 127, 64 and 706 assertions, respectively.

I recalibrated one frozen chain with the sanitized fits and independent
effective t0 scales. Its stored calibration started from scale 0.0625, t0
mean scale 0.125 and t0 hierarchy scale 0.015625. The new search selected a
non-t0 scale of 0.02210, t0 mean scale 0.005524 and t0 hierarchy scale
0.002762; the resulting exact validation acceptance was 0.546. It used 296
pilot proposals and took 16.94 seconds. On this run the t0 effective scales
were held fixed during the non-t0 refinement. A synthetic regression with
support restricted only on t0 likewise selected a non-t0 scale more than
twice its effective t0 mean scale.

This does not show that proposals became more efficient. The new calibration
had higher acceptance than the saved calibration's 0.130, but its validation
record's unscaled group-coordinate squared movement per proposal was lower
(0.00441 versus 0.00738). Those runs did not use common random numbers, and
this raw displacement is not a long-run ESS measure. The pilot metric
preferred scale 0.0221 over 0.0442; at 0.0442 the pilot acceptance was about
0.29, but the one-sided lower-confidence score returned negative infinity.
That is a conservative, noisy selection outcome, and it is another reason the
selected scale should not be treated as an established efficiency optimum.
The independent-scale change fixes the t0-to-non-t0 coupling, but the real
checkpoint still selects a small common scale for the non-t0 coordinates.

The review's support explanation is directionally right, with a bound detail
corrected: the subject support is `0.05 <= t0 < min(RT_s)` for finite
responses. The review's quoted inequality reverses both limits. Existing
support probes show that a single subject crossing can zero the joint exact
acceptance. They do not establish that excluding t0 is better than tuning its
effective scale, so t0 remains in the move. The earlier short split-t0 ESS
counterfactual was mixed and preceded the null-slope fix; no full ESS/sec
comparison has yet been run with the sanitized surrogate. The review also
overgeneralizes the ESS loss: the matched 200-draw comparison had lower group
mean ESS/sec but higher participant-52091 ESS/sec, and the 400-draw baseline
had no recorded runtime. The corrected report above keeps those limits
explicit.

The compact current-chain calibration and rank audit are
`safe_fit_calibration_summary.csv` and `safe_fit_rank_audit.csv`; the script
`benchmark_safe_fit_independent_scale.R` reproduces them from
`/tmp/collapsed_transport_followup/fresh_retuned_calibrated.rds`.

### Matched continuation with the sanitized surrogate

To check sampling efficiency, I recalibrated all three frozen chains with the
new fit safeguard and independent t0 effective scales. I cloned the resulting
states and RNG streams into two 200-draw arms, with subject HMC disabled. One
arm forced four collapsed rounds on chain 2; the other used no collapsed
rounds. Both include the same 45.6-second recalibration cost. The selected
chain-2 settings were non-t0 scale 0.04419, t0 mean scale 0.04419, t0
hierarchy scale 0.01105, path length 3, and step size 1.5. Thus the allocator
was allowed to decouple the t0 mean but did not select a smaller mean scale
for this chain.

| Arm | Sample seconds | Total seconds | Median group-mean ESS/sec | t0 group-mean ESS/sec | Non-t0 group-mean ESS/sec | Participant ESS/sec | Key max R-hat |
|---|---:|---:|---:|---:|---:|---:|---:|
| No collapsed rounds | 415.8 | 461.4 | 0.0645 | 0.0461 | 0.0660 | 0.0711 | 1.216 |
| Safe fit, four rounds on chain 2 | 399.1 | 444.7 | 0.0664 | 0.0481 | 0.0687 | 0.0734 | 1.195 |

This is not a material ESS improvement. The per-second gains of 3–4% came
with nearly unchanged raw ESS: the median group-mean ESS was 29.8 without the
move and 29.5 with it; the median participant ESS was 32.8 and 32.6. Ten of
the ten selected key variables still had R-hat above 1.05 in both arms. The
forced arm's production acceptance was 0.171, below its 0.184 calibration
validation estimate. The runs were sequential on a shared host, so the 3.6%
wall-time difference can reflect timing noise. The comparison is one short
pair, with poor convergence, and chain 2 did not choose a smaller t0 mean
scale. It shows that sanitizing the surrogate and isolating t0 scales did not
produce a clear mixing gain in this continuation; it does not establish that
either change is ineffective in longer runs or other states.

The continuation and compact summary are reproducible with
`benchmark_safe_fit_ess.R` and `summarise_safe_fit_ess.R`; their saved
summaries are `safe_fit_ess_comparison.rds` and
`safe_fit_ess_summary.csv`. Full fits are under
`/tmp/collapsed_transport_followup/`.

### Original global scale, t0, and the earlier fixtures

I compared three proposal configurations at the same frozen chain-2 state,
with the sanitized surrogate fits, path length 3, and candidate-specific
inner step-size tuning. Each arm used 64 exact proposals. The t0 mean and
hierarchy rows were assigned effective scales 0.04419 and 0.01105 in every
arm. Results are local expected acceptance and C-whitened expected squared
jump, not ESS:

| Proposal | Common scale | v.g row scale | Exact acceptance | Total C-metric ESJD |
|---|---:|---:|---:|---:|
| Unit common scale, t0 row scales held | 1.000 | 1.000 | 1.2e-11 | 1.4e-10 |
| Prior common scale, t0 row scales held | 0.04419 | 0.04419 | 0.168 | 0.379 |
| Targeted v.g scale | 0.01845 | 0.35355 | 0.436 | 0.720 |

Holding the t0 row multipliers fixed was not enough to keep subject t0 values
above 0.05: 61/64 unit-scale proposals put at least one trial's t0 below that
floor. I then froze all t0 group-chart rows while leaving the remaining rows
at unit scale. This still put at least one trial's t0 below 0.05 in 60/64
proposals. None of the 256 proposals across these four arms had a trial with
`t0 >= RT`. Among the four proposals in the frozen-t0 arm with no floor
crossing, mean exact acceptance was only `6.8e-13`. So the physical lower
bound contributes to rejection, but it does not explain why the unit common
scale fails. This is a small local probe, and it is not the same as
marginalising t0 from the model; that LBA comparison remains open.

Correction to the first version of this later safe-fit audit:
`get_pars_oo()` expands condition-specific t0 values across trial rows. Its
first counter checked only the first trial's t0 against every RT, which
misclassified the RT ceiling and could miss floor crossings in other
conditions. The corrected audit compares each trial's t0 with that same
trial's RT and checks the 0.05 floor across all rows. This changes the
64-proposal counts for this safe-fit checkpoint. The earlier 12-proposal
audit in `plans/lba_sampling_audit_2026-10-09/bounds_probe.R` already used
rowwise comparisons and remains evidence of RT-ceiling crossings in its
different saved LBA state. The revised script and CSV are
`audit_t0_bounds_original_scale.R` and `t0_bounds_original_scale_audit.csv`.

The targeted arm's total local ESJD was about 1.9 times the 0.04419 arm's;
its v.g movement divided by recent posterior variance was about 7.3 times
larger. The full directional calibration that selected it used 328 exact
pilot proposals and took 118 seconds serially on this checkpoint. A four-core
rerun failed in the worker likelihood transport (`empty (zero-byte) input
file`), so the production-pool calibration cost is unknown. The saved
four-core calibration results are not evidence of a speedup. A fixed 0.04419
scale avoids search cost and is far better than scale 1 here. The matched
ESS/sec check below finds that the measured serial directional search does
not repay its cost in this run.

The published corr15 and Forstmann table above still does not validate this
collapsed kernel: its `full` arm sets `emc2.collapsed_move = FALSE`, so those
acceptance figures belong to the ensemble, shear, and subject-HMC moves.
Separate `fullc` runs do exist in the scratchpad archive and are summarized
below. The collapsed unit tests use small synthetic Gaussian-likelihood
targets (2–5 subjects); they check invariance and algebra, not acceptance on
the 93-subject nonlinear LBA snapshot.

I also updated the fixture harness so its summary script recognizes `fullc`
outputs and reports attempted and accepted exact corrections and expected
acceptance. The historical one-line logs reported accepted iterations, which
is not the per-proposal rate. The retained EMC checkpoints did preserve the
last exact-correction block, allowing the historical rates below to be
reconstructed. A synthetic `fullc` summary smoke check verified the updated
columns.

The proposal and bound-audit scripts and outputs are
`benchmark_original_global_scale_t0_fixed.R`,
`original_global_scale_t0_fixed.csv`, `audit_t0_bounds_original_scale.R`,
and `t0_bounds_original_scale_audit.csv`.

### Matched fixed-scale versus directional ESS check

I cloned the same frozen checkpoint and RNG streams into two 200-draw arms,
with subject HMC disabled and four collapsed rounds on chain 2. The fixed arm
used common scale 0.04419. The directional arm used common scale 0.01845 and
raised the `v.g` row to effective scale 0.35355; the two t0 row scales stayed
fixed. The observed chain-2 production acceptance rose from 0.171 to 0.500,
but this did not produce a comparable gain in the target row's ESS.

| Arm | Sample seconds | Rebuild seconds charged | Total seconds | Median group ESS | Median group ESS/sec | Median participant ESS | Median participant ESS/sec | `v.g` ESS | `v.g` ESS/sec | Max key R-hat |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Fixed scale | 399.1 | 45.6 | 444.7 | 29.53 | 0.0664 | 32.64 | 0.0734 | 15.46 | 0.0348 | 1.195 |
| Directional | 410.6 | 150.1 | 560.6 | 31.95 | 0.0570 | 35.55 | 0.0634 | 15.51 | 0.0277 | 1.195 |

For the fixed arm, the 45.6 seconds are the measured three-chain rebuild
cost. To charge the directional arm fairly, I replaced chain 2's 13.7-second
fixed calibration with the 118.2-second serial directional calibration; the
other two chain costs are shared. Before setup cost, median group ESS/sec was
0.0778 for the directional arm versus 0.0740 for fixed scale (+5.2%). Including
the observed calibration cost reverses that result: median group ESS/sec is
14.2% lower, median participant ESS/sec is 13.6% lower, and `v.g` ESS/sec is
20.4% lower. The higher production acceptance therefore did not pay for its
pilot search in this run. With the measured sampling times and ESS, the
directional chain-2 calibration would need to cost about 39 seconds or less
to break even on median group ESS/sec; that is a calculated threshold, not a
measured four-core runtime.

The result is still weak evidence about long-run mixing: there are only 200
new transitions per chain, and key R-hats remain around 1.195. A planned
continuation to 400 draws was stopped before producing an extended fit. The
host was under severe memory pressure (about 61 of 63 GiB of swap in use),
and the sampling workers were blocked on page reads. No conclusion about the
400-draw ESS follows from that aborted run. The 200-draw outputs are in
`directional_ess_summary.csv` and `directional_ess_comparison.rds`; the
reproduction and continuation scripts are
`benchmark_directional_ess.R` and `extend_directional_ess_400.R`.

### Lower-cost calibration on repeated rebuilds

I changed proposal rebuilds so uncalibrated starts use the full exact search
and later rebuilds use a smaller exact screen near the last calibrated
settings, whether or not the previous allocation assigned rounds. The first
sample rebuild also warm-starts from a calibrated burn/adapt plan. This avoids
repeating the full pilot grid at a stage boundary where the target is
unchanged, while still screening the refreshed proposal chart with exact
likelihoods.
The initial comparison on a saved checkpoint was not paired and suggested the
screen retained the same scale; the paired runs below supersede that reading.

I compared the full search with the screen three times from the same saved
single-chain checkpoint. Each arm ran 200 new transitions with four collapsed
rounds per transition and the same RNG stream within its pair; the second pair
advanced the stream and reversed arm order. HMC was disabled. Each calibration
used 368 pilots for the full search and 156 for the screen. Calibration took
about 36 seconds versus 15 seconds. The table reports total elapsed time and
ESS/sec for group means, group log-SDs, all 1,674 subject-parameter series,
and the previously tracked subject 52091 summary.

| Repeat | Arm | Total seconds | Acceptance | Median group-mean ESS/sec | Median group log-SD ESS/sec | Median all-subject ESS/sec | Subject 52091 ESS/sec | `v.g` ESS/sec |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | Full | 430.8 | 0.520 | 0.0376 | 0.1741 | 0.2128 | 0.1140 | 0.0328 |
| 1 | Screen | 426.7 | 0.492 | 0.0784 | 0.1917 | 0.2827 | 0.2287 | 0.0348 |
| 2 | Full | 438.7 | 0.547 | 0.0710 | 0.1789 | 0.2848 | 0.1306 | 0.0580 |
| 2 | Screen | 405.9 | 0.500 | 0.1211 | 0.1664 | 0.3292 | 0.0134 | 0.0418 |
| 3 | Full | 430.5 | 0.570 | 0.0113 | 0.1323 | 0.1298 | 0.0168 | 0.0045 |
| 3 | Screen | 417.8 | 0.530 | 0.0525 | 0.1757 | 0.2219 | 0.1589 | 0.0516 |

The screen used nearly half the calibration time and had similar or lower
total wall time in both pairs. It improved median group-mean ESS/sec by 109%
and 70% in the first two pairs; repeat 3 added a 363% gain. Median
all-subject ESS/sec improved by 33%, 16%, and 71%. Repeat 3 also reversed the
earlier `v.g` ambiguity: its ESS/sec was 11.5 times higher under the screen.
Group log-SD ESS/sec rose 10%, fell 7%, and rose 33% across the three pairs.
Subject 52091 remained highly variable: the screen helped in repeats 1 and 3
but reduced its ESS/sec by 90% in repeat 2. The aggregate result favors the
screen in these three starts, while individual coordinates and subjects still
show substantial regressions.

The screen changed the selected proposal as well as its cost. Its common
scales were 0.01845, 0.02575, and 0.01322; the full search selected 0.01322,
0.01322, and 0.00487. Effective t0 scales also varied. Trajectory length and
the extra `v.g` row scale were the same within each pair. These runs cannot
isolate the screen policy from the selected scales. They are single-chain
200-draw experiments, so they provide no R-hat and do not establish
convergence or a general ESS/sec improvement. The third pair was completed
before the per-rebuild scale cap below. Per-variable summaries are in
`rebuild_calibration_ess_per_variable.csv`,
`rebuild_calibration_ess_repeat2_per_variable.csv`, and
`rebuild_calibration_ess_repeat3_uncapped_per_variable.csv`; paired outputs
are
`rebuild_calibration_ess_summary.csv`,
`rebuild_calibration_ess_repeat2_summary.csv`, and
`rebuild_calibration_ess_repeat3_uncapped_summary.csv`. Scripts:
`benchmark_paired_rebuild_calibration.R`,
`benchmark_rebuild_calibration_ess.R`, and
`summarise_rebuild_calibration_ess.R`.

### Cap on later production scale changes

The third uncapped pair showed that one full-search rebuild could move the
common multiplier from 0.01845 to 0.00487, while the same stream's screen
rebuild selected 0.01322. The full arm's higher acceptance (0.570 versus
0.530) came with much lower median group-mean and all-subject ESS/sec. I
replayed that full/screen pair from the same checkpoint and RNG streams with
a factor-of-two cap on later production scale updates. The first calibration
in a stage remains unrestricted so it can recover from an unusable starting
scale.

| Full-search arm | Total seconds | Acceptance | Median group-mean ESS/sec | Median group log-SD ESS/sec | Median subject 52091 ESS/sec | Median all-subject ESS/sec | `v.g` ESS/sec | Common scale |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Uncapped | 430.5 | 0.570 | 0.0113 | 0.1323 | 0.0168 | 0.1298 | 0.0045 | 0.00487 |
| Capped | 422.8 | 0.585 | 0.0810 | 0.1925 | 0.1207 | 0.2719 | 0.0542 | 0.00923 |

The capped scale was exactly half the starting 0.01845, the largest allowed
change in one rebuild. The held t0 effective scales, targeted mean scale,
trajectory length, and production RNG start were unchanged. Median
all-subject ESS/sec doubled; group-mean ESS/sec was 7.2 times higher. The
screen arm stayed near 0.01322 and had median all-subject ESS/sec 0.2239 in
this replay, so the capped full arm was also higher on that aggregate. These
are strong short-run differences from one chain and one 200-draw stream, with
no R-hat; they validate the cap against this observed overshrink but do not
establish its effect across chains or longer runs. The uncapped outputs are
saved with `_uncapped` suffixes; current capped results use the usual
`repeat3` filenames. The adaptation regression tests cover the cap and retain
unrestricted first-stage recovery.

### What the archived fullc fixture runs establish

I previously said no `fullc` fixture fits had been run. That was wrong: I
treated the published `base`/`full` table and the obvious benchmark output
directories as the complete run history. The table's `full` arm does disable
collapsed transport, but there are separate `fullc` fits in the R scratchpad
run archive. Their saved EMC checkpoints show that collapsed proposals were
attempted and accepted on several fixture targets. This corrects the earlier
claim; it does not make those runs equivalent to the current 93-subject LBA
checkpoint.

The archived final suite has three seeds for each case. I reconstructed the
exact per-proposal acceptance from each checkpoint's last measured correction
block (`plan$last$acc / plan$last$n` and `moved / n`). The harness's printed
“accepted iterations” count means iterations with at least one accepted
collapsed round; it is not the per-round acceptance rate.

| Case | Subjects / group coefficients / collapsed chart dimension | Total exact rounds across 3 seeds × 3 chains | Last measured exact acceptance, expected / observed | Notes |
|---|---:|---:|---:|---|
| T30 (RDM) | 40 / 13 / 104 | 74,477 | 0.230 / 0.230 across all 9 chains | `t0 ~ 1`; all chains had collapsed proposals |
| corr15 (RDM) | 30 / 13 / 104 | 1,711 | 0.604 / 0.616 across 6 nonempty measured blocks | The third seed had no proposals in its last measured block |
| Forstmann (LNR) | 19 / 4 / 14 | 693 | 0.364 / 0.356 across 6 chains | The first seed allocated no collapsed proposals |

These are real full-sampler runs, not only unit tests. The T30 and corr15
artifacts therefore disprove the broad idea that collapsed exact acceptance
was absent from all fixture tests. They also show why the published table
looked healthy for another reason: its `full` arm's acceptance is still HMC
acceptance, while the separate fullc artifacts were not included in that
table. The fullc targets had 19–40 subjects and at most 13 group coefficients;
their collapsed chart dimensions were 14–104. The current LBA checkpoint has
93 subjects, 18 group coefficients, three population `t0` coefficients, and a
36-dimensional diagonal-gamma chart. T30's shared `t0` had a minimum observed
subject slack (`minRT_s - t0_s`) of about 0.034 s. Thus `t0` support alone does
not force every collapsed proposal to fail; subject count, support margins,
model, hierarchy, and proposal state all matter. The archived summary is
`fullc_fixture_historical_summary.csv`; the retained raw fits were under the
temporary `.Rtmp/.../scratchpad/runs/final/` directory.

The checkpoint plans answer the scale question directly: none of the 27
chain plans has a scalar `scale` or `t0_ratio` field. These are pre-calibration
plans using the full saved metric, with fixture-specific tuned step sizes
(T30 0.046–0.115, corr15 0.212–0.284, Forstmann 0.346–0.492). That is evidence
that the successful fullc runs did not rely on the later 0.044 global shrink.
It does not mean they had the same displacement as a unit-scale LBA proposal:
their metrics, step sizes, chart dimensions, and targets differ.

I replayed 32 proposals from T30 seed 1, chain 1 at its saved final state
using the current source, the full saved metric, and its archived step size.
Mean expected exact acceptance was 0.280, and none of the proposals crossed
either measured subject-level `t0` bound. This independently reproduces
healthy collapsed acceptance on a fixture under the current fit safeguards;
it is not just the old log's HMC acceptance. The replay script and output are
`replay_fullc_unit_scale.R` and `fullc_current_source_replay.csv`.

To separate support rejection from surrogate error, I also decomposed the
exact correction by subject for 24 proposals with no trial-level t0 crossing
on each target. T30 used its full archived metric. LBA used unit scale on the
non-t0 chart rows, with t0 group-chart rows frozen and the inner step size
retuned; proposals were screened before exact evaluation only when an actual
trial crossed the `t0` floor or RT ceiling.

| Target and proposal | Subjects | Chart dimension | Attempts / t0-valid proposals | Median total log correction | Median signed correction per subject | Median absolute subject correction | Median expected acceptance |
|---|---:|---:|---:|---:|---:|---:|---:|
| T30, full metric | 40 | 104 | 25 / 24 | -1.20 | -0.030 | 0.146 | 0.301 |
| LBA, unit non-t0 rows | 93 | 36 | 394 / 24 | -26.86 | -0.289 | 0.406 | 2.8e-12 |

The typical LBA subject contributed a larger exact-minus-surrogate error, and
those errors were summed over more than twice as many subjects. This gives a
concrete reason that full-metric T30 proposals can accept while unit-scale LBA
proposals fail even after removing t0-crossing proposals. It does not isolate
the subject-count effect from the different likelihood, hierarchy, and
posterior state. The proposal-level decomposition and setup are saved in
`fullc_lba_correction_decomposition.csv` and
`fullc_lba_correction_decomposition_setup.csv`; the reproduction script is
`benchmark_fullc_lba_correction_decomposition.R`.

The saved subject fits also differ materially. In corr15 and Forstmann, none
of the 441 non-null subject fits was rank deficient. In T30, 63 of 360 were
rank deficient, but the largest fraction of fitted linear slope lying in a
flat precision direction was about 1.0e-6. By comparison, the stored LBA fits
are rank 13–17 in an 18-parameter subject space for all 93 subjects, with
unsupported-slope norms from 0.28 to 7.8. The rank audit is
`fullc_fixture_implementation_audit.csv`; the LBA audit is
`safe_fit_rank_audit.csv`. This points to a real difference in surrogate
geometry, alongside the larger subject count and different `t0` structure.
It is not a controlled attribution: the fixture and LBA targets differ in
several ways at once.

The fixture evidence does not explain away the current failure. On the 93-
subject LBA state, 400 forced unit-scale exact corrections had zero
acceptance. In the corrected local probe that froze the t0 group-chart rows,
60/64 proposals put at least one trial's t0 below 0.05, and none put a t0 at
or above that trial's RT. The four proposals without a floor crossing still
had mean exact acceptance only `6.8e-13`. Thus lower-bound crossings and poor
exact correction among the remaining proposals both contribute in this
state. The fixture suite did not reproduce the same subject count, model,
three-coordinate t0 structure, or proposal state.

The direct collapsed unit tests remain small correctness checks: their main
hierarchy fixture has four parameters and five subjects, and stationary-law
checks use two to three parameters and four to six subjects with mocked
Gaussian likelihoods. They verify invariance for those targets, not the
93-subject LBA acceptance rate. The published benchmark table and these unit
tests should not be described as the whole `fullc` history.

### Direct cost of zero-acceptance rounds

The forced full-scale sampler arm made 400 collapsed attempts on the 93-
subject checkpoint and accepted none. Each attempted round sends one proposed
row for every subject through the exact likelihood path before it can form
the joint correction. Thus this arm issued 37,200 subject-likelihood rows.
The collapsed round currently waits for the complete worker result; it does
not stop after the first subject returns `-Inf`.

I timed that path directly from the same checkpoint, using the serial fallback
and repeating the same rejected unit-scale candidate 12 times. All 12 exact
acceptance probabilities were zero. The 1,116 subject-likelihood rows took
3.70 CPU seconds in total, or 0.308 seconds per round; including proposal
integration, the full round took 0.316 seconds on average. At that serial rate,
400 rounds represent about 123 CPU seconds of exact-evaluation work. This is
a measured kernel cost, not a four-worker wall-time estimate. The script and
CSV are `benchmark_zero_acceptance_round_cost.R` and
`zero_acceptance_round_cost.csv`.

The paired 100-draw sampler comparison does not show a reliable end-to-end
slowdown: the forced unit-scale arm took 180.3 seconds wall and 1,887 cgroup
CPU seconds; the no-collapsed arm took 192.8 wall and 2,119 cgroup CPU
seconds. This one sequential pair goes opposite to the direct kernel-cost
measurement, so it cannot support a percent runtime penalty or speedup. It
does establish that the forced move did futile exact work. The benchmark
disabled subject HMC and forced the move, so it does not measure allocator
crowd-out. After a rebuild has enough rejected-round evidence, the allocator
prices the move's measured cost against its zero acceptance-weighted
displacement; the current tests and allocation probes send such a move to
zero rounds. A stale schedule can still spend rounds until that rebuild.

### Allocator direction versus calibrated direction

The capped allocator probe still selected zero collapsed rounds in all three
chains. Its worst group coordinate was shear coordinate 1 (`rho` 0.959–0.961)
in every chain, while calibration selected raw `v.g` as its slow mean. These
are different directions: allocator coordinates are `D^{-1} mu`, but the
current directional scale multiplies one raw parameter row. On chain 2, the
worst shear direction had only 0.7% of its squared norm on the t0 mean rows,
so it can be scaled separately while leaving t0 means fixed.

A matched 24-proposal local probe compared that current raw-row setting with
rank-one scaling along the support-safe projection of the worst shear
direction. Exact acceptance and accepted squared displacement in that
direction were:

| Candidate | Exact acceptance | Accepted squared jump / direction variance |
|---|---:|---:|
| Current raw-row scale | 0.633 | 0.01465 |
| No extra mean scale | 0.631 | 0.00254 |
| Shear direction x2 | 0.568 | 0.00821 |
| Shear direction x4 | 0.366 | 0.01644 |
| Shear direction x8 | 0.087 | 0.01570 |
| Shear direction x16 | 0.0002 | 0.00015 |

The factor-four probe increased movement in the allocator's exact shear
coordinate by 12% relative to the current raw-row scale. Its estimated
one-round autocorrelation factor is still about 0.992, versus about 0.95 for
the separately allocated shear round. The shared allocator therefore has a
reason to prefer shear on this checkpoint. The pilot shows a modest coordinate
alignment gain, but it does not justify the extra calibration search or show
that directional collapsed moves improve ESS per second. Results and the
reproduction script are `allocator_direction_scale_probe.csv` and
`probe_allocator_direction_scale.R`.

I then ran a matched 200-transition continuation from the same saved
three-chain checkpoint to see whether the x2 shear-direction candidate
improved ESS. Subject HMC was disabled in both arms. The candidate used four
collapsed rounds on chain 2 and none on chains 1 or 3; its 22.9-second
single-core pilot cost is charged below.

| Measure | Shear direction x2 | No collapsed rounds |
|---|---:|---:|
| Sampling seconds | 358.2 | 382.5 |
| Total seconds charged | 381.1 | 382.5 |
| Exact collapsed acceptance | 0.590 | — |
| Slow shear-coordinate ESS/sec | 0.124 | 0.137 |
| Median group-mean ESS/sec | 0.148 | 0.138 |
| Median group log-SD ESS/sec | 0.553 | 0.420 |
| Median subject ESS/sec | 0.759 | 0.757 |

The short run is mixed: the allocated shear coordinate did not gain ESS/sec,
median group-mean and group log-SD ESS/sec increased, and subject ESS/sec was
unchanged. The two elapsed times were nearly equal after charging the pilot,
but one pair cannot distinguish a runtime effect from noise. The summary
script initially failed to extract R-hat because `posterior::summarise_draws`
names that column `posterior::rhat`; these results therefore have no valid
R-hat check. I fixed the extraction for future runs. This probe does not
justify adding a directional scale or claiming a sampler-wide efficiency
gain. The CSV and reproduction script are `allocator_direction_ess.csv` and
`benchmark_allocator_direction_ess.R`.

### T30 funnel continuations and the default

The LBA result is not representative of every target: the shared allocator
still assigns it zero collapsed rounds, but a T30 full-covariance funnel uses
the move and shows better group and subject ESS per second. I compared the
calibrated collapsed allocation with a no-collapsed control, starting both
arms from the same saved three-chain checkpoint and charging the collapsed
arm for its proposal rebuild and calibration. The first two rows disable
subject HMC; the last row retains it. The two 200-draw checkpoints are
different saved configurations from the same simulated T30 seed family.

| Saved checkpoint and continuation | Setup charged to collapsed | Total seconds, collapsed / control | Median group log-SD ESS/sec, collapsed / control | Exact collapsed acceptance |
|---|---:|---:|---:|---:|
| `T30_seed1_emc`, 200 draws, HMC off | 15.4 | 60.0 / 42.2 | 1.155 / 0.667 | 0.224–0.426 |
| `T30_seed1_panel5_emc`, 200 draws, HMC off | 17.6 | 59.8 / 39.3 | 0.807 / 0.455 | 0.224–0.332 |
| `T30_seed1_emc`, 500 draws, HMC on | 15.3 | 129.2 / 103.1 | 0.862 / 0.337 | 0.163–0.324 |

Across these continuations, median group log-SD ESS per second was 73–156%
higher after charging calibration, and the 500-draw HMC-enabled run also had
2.1 times the median subject ESS per second. These are short runs from saved
checkpoints, with only one simulated target family; they do not establish the
same gain for LBA or other models. The LBA calibration cost and zero-round
allocation remain a real downside for fits that do not have this funnel.

### Charging calibration setup to the allocator

The allocator previously counted a collapsed move's per-iteration context and
frame setup, but omitted the one-time exact-pilot calibration that can run at
the same proposal rebuild. The rebuild now measures that work in reserved
core-seconds and adds `calibration cost / B` to `C_rest`, where `B` is the
number of transitions before the next rebuild. This lets the round allocation
compare ESS against the actual configured block cost. The one-time charge is
logged separately from context/frame setup in allocation history.

I reran the matched T30 stage-transition comparison with this accounting. The
full search used 1,224 pilots and took 17.0 seconds for sample setup; the
screen used 508 pilots and took 11.1 seconds. Both allocated zero collapsed
rounds and produced identical group and subject ESS because the production
streams and non-collapsed state were shared. The setup measurement therefore
shows a 5.9-second saving from screening, with no sampling difference in this
checkpoint. Total elapsed time was 119.9 seconds for full search and 89.5 for
screen, but the measured sampling phase itself differed by 24.5 seconds
despite identical draws. Treat that single-run end-to-end difference as timing
noise, not a demonstrated production speedup. The row-level results are in
`stage_transition_t30_calibration_charged.csv`; allocation cost accounting is
covered by `test-move-allocation.R`.

Given the substantial end-to-end gain on the target the move was designed to
help, it is enabled by default, and users can disable it with
`options(emc2.collapsed_move = FALSE)`. The decision does not claim universal
gains. Checkpoints continue to carry the recorded production option so
resumed runs retain their kernel.

I also probed splitting one LBA subject's HMC into separate `t0` and non-`t0`
blocks. It reduced likelihood rows per transition but lowered ESS for the slow
`t0_LoadLow` coordinate in both tested chains, so that simple block split was
not adopted. The frozen-checkpoint probe is recorded in
`/tmp/lba_sampling_audit/hmc_split_conditional_ess.rds`.

### Coupling `t0` to its correlated subject parameters

Because the raw `t0` split may sever posterior correlations, I tested a
two-block conditional HMC kernel on subject 52091. The first block contains
the three `t0` coordinates and the parameters most correlated with them in
this checkpoint (`B`, `B_lR_target`, `B_LoadLow`, `B.w_prac`, `sv_lMTRUE`,
and `A`); the second block contains the other nine parameters. Each update
uses the exact conditional Gaussian group prior and the exact subject
likelihood, so this is an invariant kernel for the same subject conditional.
The benchmark fixes the group state and all other subjects, and is not a
whole-sampler comparison.

I ran 400 retained transitions from each of two saved chains, tuning each
block with either the production 20-move warmup or 100 warmup moves. Step
sizes were then held fixed during retained transitions. The table reports
bulk ESS summed over the 18 coordinates per second and the median ESS per
second for the three `t0` coordinates.

| Chain | Warmup moves | Blocks | `eps` | Acceptance | Aborts | Seconds | Sum ESS/sec | Median `t0` ESS/sec |
|---|---:|---|---|---|---|---:|---:|---:|
| 2 | 20 | full | 0.4685 | 0.125 | 178 | 9.64 | 122.2 | 2.31 |
| 2 | 20 | coupled 9+9 | 0.3605+0.8721 | 0.775+0.770 | 3+0 | 10.20 | 252.3 | 3.91 |
| 3 | 20 | full | 0.5170 | 0.630 | 7 | 8.93 | 256.1 | 6.63 |
| 3 | 20 | coupled 9+9 | 0.5279+1.0640 | 0.738+0.700 | 0+0 | 8.69 | 404.2 | 12.85 |
| 2 | 100 | full | 0.3332 | 0.762 | 2 | 14.29 | 194.2 | 1.00 |
| 2 | 100 | coupled 9+9 | 0.4142+1.0880 | 0.700+0.672 | 7+0 | 11.27 | 261.2 | 5.90 |
| 3 | 100 | full | 0.5039 | 0.652 | 4 | 11.08 | 278.8 | 11.05 |
| 3 | 100 | coupled 9+9 | 0.4833+0.8027 | 0.797+0.830 | 0+0 | 11.18 | 238.4 | 6.87 |

The 20-move full-block result on chain 2 is poorly tuned: extending the same
warmup rule raised its acceptance from 0.125 to 0.762 and reduced support
aborts from 178 to 2. Production does not freeze `eps` forever after this
warmup: at proposal rebuilds it uses accumulated transition evidence and
halves `eps` when mean acceptance is below 0.2. Thus the 20-move full-block
result exaggerates the persistent failure of the production kernel.

The coupled block improves median `t0` ESS/sec in chain 2 under both warmup
lengths. In chain 3, however, it improves `t0` ESS/sec only with the short
warmup, while the 100-move run is worse than full HMC on both median `t0` and
summed ESS/sec. It also reduces ESS/sec for some drift coordinates, which
remain correlated with `B` and `sv`. This chain- and tuning-sensitive result
does not support adding a fixed block partition to the production sampler.
The short single-subject runs also have no R-hat or whole-fit ESS/sec check.
The benchmark script and compact outputs are
`benchmark_hmc_t0_coupled_block.R`, `hmc_t0_coupled_block_summary.csv`, and
`hmc_t0_coupled_block_summary_warmup100.csv`.

### Calibration validation budget probe

I tested whether a local calibration screen could halve its independent
two-round validation batches from 16 to 8. This saves 16 exact pilot proposals
per chain. On the saved LBA allocation state, both budgets assigned zero
collapsed rounds to all chains. The 8-batch run used 444 total pilots across
the chains, compared with 492 for 16 batches; it therefore showed no sampling
gain to weigh against the smaller validation sample.

I also tested the two settings on the saved T30 checkpoint. Both left chains 1
and 3 at zero rounds, while chain 2 received 0.50 rounds with 8 batches and
0.75 with 16. I continued each resulting plan for 200 draws with subject HMC
off. The following ESS/sec values include the measured local-screen time
(chain calibrations were run sequentially); the shared initial rebuild is
excluded from both totals.

| Validation batches | Pilot proposals | Screen seconds | Sample seconds | Group-mean ESS/sec | Group log-SD ESS/sec | t0 group-mean ESS/sec | All-subject ESS/sec | Max R-hat, group / subject |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 436 | 14.52 | 38.65 | 0.449 | 0.443 | 1.304 | 0.883 | 1.50 / 1.63 |
| 16 | 484 | 14.18 | 39.91 | 0.543 | 0.231 | 1.374 | 0.684 | 1.70 / 1.60 |

The 8-batch arm gained group log-SD and all-subject ESS/sec; the 16-batch arm
gained group-mean ESS/sec and was slightly higher for t0 group means. Screen
time did not fall in the T30 comparison, and the plans also selected different
scales and path lengths, so this does not isolate the batch count. Both runs
are only 200 draws and have high R-hat. I therefore kept 16 validation batches
as the production default. The private `validation_batches` argument exists
for this diagnostic only. Scripts and outputs are
`benchmark_validation_batches.R`, `benchmark_validation_batch_ess.R`,
`lba_validation_batch_comparison.csv`, `t30_validation_batch_comparison.csv`,
and `validation_batch_ess_comparison.csv`.

### Subject-HMC trajectory duration probe

I exposed a fixed trajectory-duration option, `emc2.subject_hmc_length`, so a
longer path could be compared without changing the default `pi/2` kernel. The
duration is stored with production options and carried to HMC workers; old
checkpoints inherit `pi/2`. A focused subject-HMC test verifies the leapfrog
count, worker routes, and legacy-checkpoint default.

The whole-sampler check continued the same frozen 93-subject LBA checkpoint
for 200 sample transitions per chain at `pi/2` and `3*pi/4`. Both arms forced
the same sample proposal rebuild, disabled collapsed transport, enabled HMC,
and used four likelihood cores per chain. ESS below uses only the 200 new
draws per chain. Rebuild/calibration and sampling are included in total CPU and
wall time.

| Path | Wall seconds | CPU seconds | Group-mean min / median ESS | Group log-SD min / median ESS | Subject 52091 min / median ESS | Key min / median ESS | Max R-hat, key | HMC attempts / acceptance / aborts |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `pi/2` | 315.0 | 3444.1 | 19.7 / 48.1 | 86.2 / 204.5 | 47.1 / 172.4 | 19.7 / 133.4 | 1.138 | 200 / 0.725 / 0.020 |
| `3*pi/4` | 320.0 | 3476.8 | 13.1 / 74.0 | 73.6 / 202.9 | 31.6 / 284.9 | 39.2 / 130.2 | 1.075 | 200 / 0.680 / 0.065 |

The longer path raised the median ESS for group means and participant 52091,
but lowered their minimum ESS; group log-SD ESS was nearly unchanged. For the
five selected group means and participant parameters, changes also went in
both directions. It used about 1% more CPU, with HMC's own CPU rising from
8.1 to 10.9 seconds and HMC calibration CPU rising from 27.2 to 34.2 seconds;
only chain 3 received HMC attempts under the shared allocator in these
continuations. This single short comparison does not support changing the
production default. It does show why overall HMC
acceptance alone is insufficient: the longer path had slightly lower
acceptance and more support aborts, while its ESS effects depended on the
coordinate.

An earlier paired 100-draw screen favored the longer path on its worst group
mean and selected subject parameters, but that advantage did not consistently
hold at 200 draws. Isolated subject-conditional runs also favored `3*pi/4` on
the worst coordinate across four tested subject-chain states, while a 7- or
20-transition pilot score did not reliably select that duration. These
conditional and short whole-fit results are exploratory; the default remains
unchanged. Scripts and compact summaries are `benchmark_hmc_path_fit.R`,
`summarise_hmc_path_fit.R`, `hmc_path_fit_comparison.csv`, and
`hmc_path_fit_key_parameters.csv`.

### Three-chain full-search versus warm-start calibration

I continued the same frozen three-chain LBA checkpoint for 200 sample draws per
chain, comparing a full calibration search with the local neighborhood screen
used when a calibrated plan is rebuilt. Both arms forced a sample-stage
proposal rebuild, disabled subject HMC, used four likelihood cores per chain,
and reset each retained chain to its pre-calibration Gibbs RNG stream. The
calibration work is included in total wall and cgroup CPU; pilot draws are not
retained.

| Calibration | Pilot proposals, chains 1 / 2 / 3 | Rebuild wall | Sample wall | Total wall | Cgroup CPU | Collapsed rounds per iteration, chains 1 / 2 / 3 |
|---|---:|---:|---:|---:|---:|---:|
| Full search | 376 / 384 / 392 | 30.1 s | 298.0 s | 328.2 s | 3528.3 s | 0 / 0 / 0 |
| Local screen | 164 / 164 / 164 | 18.2 s | 298.1 s | 316.3 s | 3436.1 s | 0 / 0 / 0 |

The local screen cut pilot count by 58% and rebuild wall time by 11.9 seconds
(about 40%). End-to-end wall time fell by 11.9 seconds, or 3.6%; cgroup CPU
fell by 2.6%. Both plans assigned zero collapsed rounds, so neither arm gained
sampling moves. Their short-run ESS summaries varied in both directions: the
full search had higher median ESS for group means (60.8 vs 55.1), group
log-SDs (217.8 vs 201.9), and all subject coefficients (327.4 vs 307.3),
while the screen had higher median ESS for subject 52091 (85.5 vs 76.4).
Minimum ESS and R-hat also differed modestly. Since both arms assigned no
collapsed transitions and this was only 200 draws per chain, these differences
do not establish a sampling-quality effect.

This comparison confirms that warm-start screening reduces calibration
overhead, but calibration is not the dominant runtime cost in this continuation:
sampling accounts for 298 seconds, and the screen saves only 12 seconds of
end-to-end wall time. The bigger inefficiency in this state is that the
allocator spends pilot work on a collapsed proposal and then assigns it zero
rounds. The current screen already limits this cost at later rebuilds; I did
not reduce its pilot budget further because the earlier validation-budget
probe changed ESS results in mixed directions and did not isolate the effect.
The benchmark and output tables are `benchmark_multichain_calibration_screen.R`,
`multichain_calibration_screen_comparison.csv`, and
`multichain_calibration_screen_chains.csv`.

### Factor-four scaling on the allocator's slow shear direction

The earlier rank-one pilot found its largest exact accepted squared movement
on the current slow shear direction at a factor of four, but the whole-sampler
continuation had tested only factor two. I therefore continued the same frozen
three-chain LBA state for 200 draws per chain with four collapsed rounds on
chain 2 only, scaling the t0-safe projection of that shear direction by four.
The other chains received no collapsed rounds. A 64-proposal private exact
pilot cost was measured and charged to the candidate arm; subject HMC was off
in both arms.

| Measure | Slow-direction x4 | No collapsed rounds |
|---|---:|---:|
| Sampling seconds | 337.8 | 358.5 |
| Exact pilot seconds charged | 22.1 | 0 |
| Total seconds charged | 359.9 | 358.5 |
| Exact collapsed acceptance | 0.523 | — |
| Slow shear-coordinate ESS/sec | 0.1454 | 0.1463 |
| Median group-mean ESS/sec | 0.1571 | 0.1468 |
| Median group log-SD ESS/sec | 0.4778 | 0.4485 |
| Median subject 52091 ESS/sec | 0.433 | 0.372 |
| Median all-subject ESS/sec | 0.852 | 0.808 |

The scaled proposal was essentially tied with the control in charged wall time
and did not improve ESS/sec for its intended slow coordinate. Group-level and
subject summaries moved in both directions, and the minimum group log-SD
ESS/sec was lower with the move. This single 200-draw continuation is noisy;
it does not support promoting factor-four directional scaling into the
sampler. The pilot's estimated target movement also varied with its seed: the
64-proposal batch estimated 0.0228 of the target coordinate variance, versus
0.0164 in the earlier 24-proposal batch. The continuation's target-coordinate
ESS/sec remained unchanged despite realized acceptance above one half. This
reinforces that exact acceptance and pilot ESJD are not sufficient evidence
of improved slow-direction mixing. The reproducible benchmark and outputs are
`benchmark_allocator_direction_x4.R`, `allocator_direction_x4_summary.csv`,
`allocator_direction_x4_timing.csv`, and `allocator_direction_x4.rds`.

### Whole-sampler subject-HMC block probe and warmup-state persistence

The single-subject block experiments above motivated a whole-sampler check.
The LBA continuation used 93 subjects and 200 new draws per chain, with the
same saved starting checkpoint, subject HMC enabled, collapsed transport
disabled, four likelihood cores per chain, a forced sample proposal rebuild,
and the default `pi/2` duration. The experimental arm updated a first block
containing the three `t0` coordinates and six correlated parameters, then
updated the remaining nine parameters as a second block. Each block has its
own step size and conditional Gaussian prior, so narrow `t0` support no longer
forces one step size on every coordinate. These were opt-in research arms;
the block-partition feature is not being retained in the package because the
three matched pairs did not show a consistent package-wide gain.

| Same-seed pair | Arm | Wall / CPU seconds | Group mean min / median ESS/sec | Group log-SD min / median ESS/sec | Subject 52091 min / median ESS/sec | Key min / median ESS/sec | Max R-hat: group mean / subject / key |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | Full | 313.3 / 3428.4 | 0.020 / 0.071 | 0.145 / 0.572 | 0.038 / 0.664 | 0.032 / 0.191 | 1.447 / 1.189 / 1.231 |
| 1 | Coupled blocks | 316.6 / 3473.4 | 0.155 / 0.257 | 0.172 / 0.609 | 0.160 / 0.441 | 0.160 / 0.318 | 1.042 / 1.059 / 1.059 |
| 2 | Full | 309.9 / 3402.4 | 0.079 / 0.178 | 0.215 / 0.722 | 0.137 / 0.428 | 0.102 / 0.325 | 1.108 / 1.053 / 1.103 |
| 2 | Coupled blocks | 313.2 / 3434.3 | 0.137 / 0.194 | 0.402 / 0.823 | 0.089 / 0.414 | 0.089 / 0.366 | 1.069 / 1.098 / 1.098 |
| 3 | Full | 308.8 / 3386.0 | 0.047 / 0.121 | 0.158 / 0.732 | 0.130 / 0.353 | 0.082 / 0.269 | 1.206 / 1.109 / 1.113 |
| 3 | Coupled blocks | 322.1 / 3509.8 | 0.030 / 0.061 | 0.103 / 0.502 | 0.081 / 0.405 | 0.033 / 0.263 | 1.287 / 1.092 / 1.241 |

The block arm cost about 1% more wall time and CPU in both pairs. It improved
the weakest group-mean and group log-SD ESS/sec in both pairs, and improved
the group-mean R-hat. Participant-level results were mixed: the first pair's
median participant ESS/sec fell by 34%, while the fresh-seed pair's minimum
participant ESS/sec fell by 35%. The key-coordinate median improved in both
pairs, but its minimum fell in the repeat. Thus the apparent large first-pair
gain shrank substantially on a fresh seed. Pair 3 was worse on the group
summaries and key minimum, while participant median ESS/sec was higher. The
mixed results and small wall-time penalty do not justify shipping a block
partition option.

Acceptance also needs to be read per block. In the first pair, the selected
subject's first-block mean acceptance probability was 0.441, compared with
0.895 for the remaining coordinates; combining them into the reported
whole-transition acceptance hid this difference. In the fresh-seed pair the
block arm had 200 HMC attempts, 185 accepted transitions and 20 aborts; the
full arm had 386 attempts, 302 accepted and 12 aborts. These aggregate counts
are affected by the shared move allocator and are not a controlled measure of
block quality. HMC was allocated only on chain 3 in these short runs. Neither
run included collapsed transport, so they do not resolve the LBA collapsed
move's acceptance or cost.

While wiring the experiment, I also tried persisting unfinished HMC warmup
state across calls. The candidate change did not pass the focused split-run
check, so I removed it; no HMC warmup behavior is retained from that probe.
The original subject-HMC suite passes with its original test setup. The
benchmark outputs are
`hmc_full_200.rds`, `hmc_t0_coupled_200.rds`,
`hmc_full_200_repeat2.rds`, `hmc_t0_coupled_200_repeat2.rds`,
`hmc_full_200_repeat3.rds`, and `hmc_t0_coupled_200_repeat3.rds`.

### Fullc fixtures versus the LBA fit: surrogate comparison

The fullc fixture results remain important evidence: the saved T30 chain had
40 subjects and a 104-coordinate group chart, and its unshrunk full metric had
no `t0` support crossings in the current-source replay. On eight proposals,
the empirical regression surrogate gave mean exact acceptance 0.280 and
median log correction -1.72. Replacing only its likelihood fits with the
saved mode-centered Laplace fits gave mean exact acceptance 0.324 and median
log correction -2.05. Both variants were usable at the original scale. This
small probe agrees with the historical fixture acceptance range; it does not
replace the recorded multi-seed fullc runs.

I ran the same fit comparison on the 93-subject LBA checkpoint, holding the
learned non-t0 metric at scale one and freezing the t0 chart rows. The
candidate still transports each subject through the exact conditional frame,
so this isolates support-safe proposals rather than changing the target.

| Fit used to define the transport | Attempts | Proposals within t0 floor | Median log correction on safe proposals | Mean exact acceptance on safe proposals |
|---|---:|---:|---:|---:|
| Empirical regression | 244 | 8 | -28.3 | 7.2e-7 |
| Mode-centered Laplace | 500 | 2 | -92.2 | 5.4e-36 |

Neither arm crossed an observed-RT ceiling; every screened proposal failed at
the hard 0.05-second floor. The safe-proposal counts are small, especially
for Laplace, so the numerical acceptance estimates are local diagnostics.
The difference is nevertheless large enough to reject “replace the regression
fits with Laplace fits” as a fix for this checkpoint. T30 and LBA are not a
controlled model comparison: LBA has 93 subjects, 18 group parameters, and
three t0 coefficients, while the T30 fixture has 40 subjects and 13 group
parameters. The fullc fixture history does not imply that every target has
the same surrogate error.

This also clarifies the role of the scalar. The proposal integrator advances
`phi` by `eps * C * momentum`; multiplying `C` by the selected scalar directly
shrinks every chart direction relative to the learned covariance geometry,
while dual averaging retunes `eps`. It does not tie raw parameters such as
`t0` and `B` to a common physical unit. The LBA safe-proposal correction
remains poor at scale one, so removing only t0-floor collisions is not enough
to justify restoring the full scale. The scripts and outputs are
`benchmark_lba_regression_vs_laplace_surrogate.R`,
`lba_regression_vs_laplace_summary.csv`,
`t30_regression_vs_laplace_summary.csv`, and the corresponding proposal and
subject-correction CSVs.

### t0-floor preflight experiment (discarded)

The historical microbenchmarks measured a large cost reduction by skipping
all exact subject-likelihood calls whenever any transported t0 fell below the
model's lower bound. That shortcut was not target-preserving. The sampler's
`is_ok` verdict is per data row, and the compiled LBA kernel separately checks
`rt - t0 <= 0` before evaluating that row's density. The subject likelihood
still sums its other trials; a support crossing on one trial does not make the
whole subject likelihood zero. The implementation also applies a finite floor,
and a nonzero guess component can assign likelihood to responses before t0.
So the measured bypass could not safely return a zero joint acceptance.

I removed the whole-subject shortcut from production calibration and
sampling, together with its state flags and tests. Exact likelihood dispatch
continues for every transported subject, while the existing kernels skip
unsupported trial densities. The old benchmark scripts and tables remain in
this directory only as an audit trail for a discarded optimization; their
speed measurements do not describe the current package.

The data is available during `design(data = ...)`, but the model bounds it
builds are global and apply to mapped parameter values. A per-subject upper
cap could be expressed for a simple intercept-only t0 model by carrying
subject-specific data bounds into each likelihood call. With condition or
covariate designs, support becomes a collection of row-wise constraints on
linear predictors, not one scalar parameter bound. Guess mixtures and other
race models also change whether `t0 < min(RT)` is a hard support condition.
Since the compiled LBA path already avoids the density call for each
`rt <= t0` row, a general dynamic cap would add a second support mechanism
without fixing the collapsed transport. I did not add that package-wide
mechanism.

### Package decision, 2026-10-10

I reverted the exact-pilot global scale, separate t0 and slow-mean multipliers,
trajectory-length search, diminishing scale updates, and calibration-cost
allocation changes. The retained LBA continuation showed no material ESS gain
after calibration cost, and the measured directional ESS/sec was effectively
tied with its no-collapsed control. The unit-scale LBA probe also remained
poor even when t0 was frozen, so neither “just bound t0” nor a smaller shared
scale is established as the general fix. In contrast, the archived `fullc`
fixtures used the original unshrunk transport successfully. These findings do
not support changing the package-wide scale based on the current LBA fit.

The package retains the collapsed-surrogate safeguard, which removes linear
slopes in flat fitted-curvature directions from empirical, fallback, and
saved fits. The HMC warmup-continuation attempt, manual HMC block and
trajectory-length options, and unsafe whole-subject t0-floor preflight were
removed after the focused checks or whole-sampler comparisons did not support
shipping them. No dynamic per-subject upper-bound mechanism was added: the
current compiled LBA path already avoids density work for each row with
`rt <= t0`, while enforcing a new hard cap would need to account for the full
t0 design and guess-mixture semantics.

Concretely, `LBA()` declares the global natural-scale bound `t0` in
`(0.05, Inf)` (with the existing zero exception). The `is_ok` mask comes from
those supplied parameter bounds; it does not inspect reaction times. The
compiled race path separately checks each row's `rt - t0` and skips that
row's race density when the difference is nonpositive. It then evaluates the
other rows and applies any enabled guess or contaminant mixture. So this path
does not impose a subject-level `t0 < min(rt)` cap, and a blanket cap would
change the target for supported mixture models. A correct data-derived bound
would need to be attached to mapped subject/trial values only in models where
pre-`t0` responses have exactly zero likelihood.

The package changes therefore clean up the unsupported scale and support
shortcuts without claiming that the current 93-subject LBA collapsed move has
been fixed. Its scale-one exact correction remains a separate open problem;
the current evidence says that replacing scale calibration alone cannot solve
it. The archived calibration, preflight, and block experiments remain here
for reproducibility, but they are not production features.

### Follow-up, 2026-10-10: fit-context reuse and current chain allocations

The source review found that the flat-direction safeguard was being
recomputed for every collapsed move. The context builder sanitized local
copies of the subject fits, but the move did not put those sanitized base fits
back into the chain state. The move now caches them when the stored fit list
has one entry per subject. It conditions a fresh copy on the current pooled
coefficients each iteration, and a proposal rebuild replaces the cached fits.
Fallback fits computed on demand are not persisted. A 20-context microcheck
on the frozen 93-subject checkpoint took 0.317 seconds before reuse and 0.026
seconds with the sanitized fits cached. This measures context construction,
not whole-fit ESS per second; it leaves the target and proposal scale
unchanged.

The allocator also no longer replaces each chain's collapsed-round choice
with the largest choice from another chain. A regression test confirms that
rates such as 1.125 and 0 remain distinct through rebuilding. Each chain's
own solver result already prices its movement and cost; a cross-chain maximum
could undo a zero allocation and spend time on a move that chain rejected.
The earlier local-versus-shared ESS comparison remains evidence for that
case, but it used the discarded 0.044 calibrated scale. On a fresh rebuild of
the frozen LBA checkpoint with the restored unit-scale proposal, all three
chains selected zero rounds, so this policy change did not alter that LBA run
and no ESS improvement is claimed for it. A current-source replay on the T30
`fullc` fixture, where collapsed rounds were allocated, is summarized below.
Current-source rebuilds on the corr15 and Forstmann `fullc` checkpoints chose
zero rounds for all chains; those targets therefore provide no local-versus-
shared continuation comparison.

The focused `move-allocation`, `collapsed-move`, and `production-options`
suites pass with the chain-local policy and cached-fit conditioning. All 64
subject-HMC assertions also pass with the original split-run test setup (one
development-only model-validation test is skipped). The cache test changes
pooled coefficients between contexts and checks that only the base fit is
reused.

### Current-source `fullc` check of chain-local round allocation

The archived T30 `fullc` checkpoint has 40 subjects and 13 group parameters.
A fresh current-source rebuild selected 7, 7.5, and 6.5 collapsed rounds per
transition for chains 1–3. The former shared-maximum rule would assign 7.5 to
all chains, adding 1.5 expected rounds per fit transition. I continued the
same rebuilt checkpoint under both policies for 200 transitions per chain,
charging each arm for the common proposal-rebuild cost. The second pair used a
new RNG stream and reversed arm order. Both arms used the same HMC allocation
within each pair.

| Pair | Policy | Total wall seconds | Median group-mean ESS/sec | Median group log-SD ESS/sec | Median all-subject ESS/sec | Minimum all-subject ESS/sec |
|---:|---|---:|---:|---:|---:|---:|
| 1 | Chain-local | 79.1 | 2.879 | 1.786 | 3.842 | 0.296 |
| 1 | Shared maximum | 86.3 | 1.976 | 2.486 | 3.733 | 0.383 |
| 2 | Chain-local | 82.3 | 4.271 | 3.697 | 4.793 | 1.120 |
| 2 | Shared maximum | 80.9 | 4.172 | 3.227 | 4.423 | 0.619 |

The local policy improved median all-subject ESS/sec in both pairs, by 3% and
8%, and improved its minimum in both pairs. Group-mean ESS/sec was higher in
both pairs; group log-SD ESS/sec fell in the first pair and rose in the
second. The wall-time result reversed: local allocation was 8% faster in pair
1 and 2% slower in pair 2. This supports retaining chain-local rates as a
reasonable way to avoid work a chain's own allocation rejected, but it does
not establish a general runtime or ESS guarantee. These are short runs from
one T30 checkpoint; some R-hats remained above 1.1 in the first pair.

### LBA surrogate-window screen, 2026-10-10 (inconclusive)

I compared a quadratic fit to the most recent 250 production draws with the
existing all-history fit, using a fresh current-source proposal metric and
holding the three `t0` group-mean chart rows fixed. At unit scale, 12
support-safe proposals per fit and seed gave a less-negative median log
correction for the recent window on three paired seeds (`-246`, `-262`,
`-235`) than for all history (`-308`, `-293`, `-326`). All six estimates
still implied negligible exact acceptance. In a separate eight-proposal
scale screen, the ordering reversed at unit scale and all-history had more
acceptance-weighted group movement at every tested scale. These samples are
too small and seed-sensitive to support changing the production fit window;
neither candidate made the unit-scale LBA move usable. No fit-window or scale
change was retained.

### LBA mode-anchored curvature screen, 2026-10-10 (rejected)

I tested the independent review's narrower suggestion: retain the empirical
quadratic curvature, but set its linear term so the surrogate gradient matches
the mode-centered Laplace likelihood fit at that fit's conditional mode under
the current group prior. I rebuilt the empirical fits and proposal metric from
the current-source window across all three saved chains. In the 36-coordinate
LBA chart, I froze the three `t0` group-mean coordinates and their three
diagonal-scale coordinates, then screened proposals from one chain at unit
scale. The matching empirical-fit arm produced two proposals within the
`t0` floor in 40 attempts; their median exact log correction was `-176.7` and
their mean expected acceptance was `5.1e-69`. The mode-anchored arm had no
support-safe proposal in 250 attempts. A separate current-state-gradient
hybrid screen also had none in 250 attempts. These are small proposal screens,
but they provide no basis for a production slope change.

Freezing the `t0` group coordinates did not freeze each subject's `t0` value:
the conditional Gaussian frame still shifts it through cross-parameter
likelihood curvature as other group coordinates move. This is one reason a
simple group-coordinate mask cannot guarantee support-safe subject proposals.
The exact correction on the few safe empirical proposals was also negligible,
so removing the `t0` group directions alone would not make this LBA collapsed
move useful. The current-source allocator assigns this checkpoint zero
collapsed rounds, which avoids paying exact likelihood cost for the unusable
move during transitions. Neither the mode-anchor nor the earlier group-chart
support mask was retained.

### Conditional transport for bounded subject coefficients, 2026-10-10

I replaced the proposed global `scale * C` response to support crossings with
a conditional transport. The move holds selected subject coefficients fixed
and transports the remaining coefficients under their Gaussian conditional.
Its inner group target is the collapsed surrogate multiplied by the Gaussian
marginal density of the fixed coefficients. The exact correction remains the
existing sum of likelihood-minus-surrogate differences. Finite-difference
checks verify the analytic group gradient, and a forward/reverse map check
verifies the subject transport. In production, the fixed set includes `t0`
coefficients (the observed-RT ceiling is subject-specific and absent from the
global bound matrix) plus coefficients whose hard model bounds are narrower
than their transforms. The original covariance metric is unchanged.

On the frozen 93-subject LBA snapshot, fixing only the three subject `t0`
coefficients removed the floor crossings in a 20-proposal unit-scale screen,
but expected exact acceptance was still only `4.9e-34`; there were no global
bound or RT-ceiling crossings. The unconditioned arm had a floor crossing in
19 of 20 proposals and expected acceptance `6.0e-97`. At scale `0.25`, the
`t0`-only arm still had expected acceptance `4.6e-5`. Thus removing t0
crossings is a valid support safeguard, but it does not resolve the LBA
surrogate correction error.

Fixing all eight coefficients selected by the stricter-hard-bound rule gave
expected acceptance around `8e-4` at unit scale and `0.207` at scale `0.25`
in separate 20-proposal screens. The latter had chart ESJD per second `0.54`;
scale `0.044` had higher acceptance (`0.73`) but lower chart ESJD per second
(`0.10`). A local Hessian metric did not improve on the quarter-scale screen.
These are proposal diagnostics, not ESS-per-second results. In particular,
they do not justify shipping a quarter-scale multiplier.

I also screened unit-scale `t0`-conditional proposals on the archived `fullc`
T30, corr15 and Forstmann checkpoints (12 proposals per arm, one chain each):

| Fixture | Unconditioned acceptance / chart ESJD/s | Conditional acceptance / chart ESJD/s |
|---|---:|---:|
| T30 | 0.128 / 7.75 | 0.169 / 9.63 |
| corr15 | 0.736 / 37.52 | 0.861 / 23.10 |
| Forstmann | 0.031 / 3.87 | 0.750 / 32.94 |

No candidate in these three screens crossed the model bounds, t0 floor, or
observed-RT ceiling. The conditional move helped T30 slightly and Forstmann
substantially, while its extra conditional work reduced corr15's movement
proxy. The broader hard-bound mask also fixed six DDM drift coefficients in
T30/corr15; those model bounds are legitimate constraints, but freezing them
can alter proposal efficiency. Allocation remains chain-local and can select
zero rounds when the move does not pay for its cost. These short screens do
not establish a package-wide ESS/sec gain.

The likelihood's existing `ok` mask reflects model-wide `minmax` bounds; it
does not construct a separate `t0 < min(rt_s)` bound for each subject. The
LBA kernel checks `rt - t0` per observation, and a pre-onset row need not make
the whole subject likelihood zero when other rows or the guess component
contribute. Therefore I did not add a per-subject hard cap to `design_model()`
or a whole-subject likelihood bypass. The conditional transport instead keeps
the current subject `t0` values unchanged during this group move.

Finally, when a fractional allocation draws zero rounds, the move now records
the iteration and RNG state and returns before rebuilding the conditional
context. The `collapsed-move` and `move-allocation` suites pass with this path.
The LBA bottleneck remains open: the support fix removes floor collisions, but
the exact correction is still too poor at the original scale, and the scale
that helped this checkpoint was not consistently best on the `fullc` controls.

### Follow-up: subject-HMC support blocks and the observed-RT ceiling

The existing collapsed-move selector identifies 11 support-sensitive LBA
coefficients in this snapshot: the `A`, `B`, and `sv` coefficients plus all
three `t0` coefficients. I compared one full-vector conditional HMC transition
with two successive blocks (those 11 coefficients, then the other 10) on four
subjects from the saved checkpoint. Each arm had 20 warmup transitions and 200
measured transitions; block step sizes adapted separately. The table reports
median coordinate bulk ESS per second over all coefficients, the selected
support-sensitive coefficients, the remainder, and `t0` alone. Time includes
the warmup. These are single-chain conditional-kernel screens, not full-fit
ESS estimates.

| Chain / subject | Kernel | Seconds | Median ESS/s | Sensitive ESS/s | Other ESS/s | t0 ESS/s | Aborts / 200 |
|---|---|---:|---:|---:|---:|---:|---:|
| 3 / 14 | Full | 16.90 | 1.09 | 0.46 | 2.68 | 0.16 | 1 |
| 3 / 14 | Split | 14.71 | 4.03 | 1.46 | 10.05 | 0.93 | 0 |
| 2 / 63 | Full | 7.04 | 5.53 | 3.86 | 5.68 | 9.18 | 17 |
| 2 / 63 | Split | 6.53 | 10.49 | 10.31 | 12.70 | 10.31 | 1 |
| 1 / 60 | Full | 6.42 | 13.71 | 13.29 | 13.82 | 7.37 | 0 |
| 1 / 60 | Split | 5.30 | 10.84 | 9.60 | 20.40 | 17.52 | 0 |
| 3 / 48 | Full | 8.07 | 11.91 | 8.40 | 14.64 | 6.79 | 0 |
| 3 / 48 | Split | 6.63 | 3.01 | 2.96 | 18.77 | 11.04 | 3 |

The split helped the two subjects with the most full-vector aborts, but its
overall ESS/s was mixed and it reduced movement in some sensitive coordinates.
This does not justify a default subject-HMC change. The diagnostic prototype
was run inline and left no HMC block code in `R/` or `src/`.

### Ridge regularization of the LBA quadratic surrogate, 2026-10-10

I screened whether the empirical quadratic's full set of interaction terms
was overfit. For each of the three chains in the saved 93-subject LBA
checkpoint, I fitted each subject's surrogate on the other two chains and
measured log-likelihood prediction on the held-out chain. A ridge penalty of
0.1 on standardized regression features lowered median per-subject held-out
RMSE from 1.436 to 1.388 and the 90th percentile from 2.166 to 2.034; it
improved RMSE for 88 of 93 subjects. This is evidence of modest predictive
overfit in the unregularized regression, but posterior-draw prediction is not
the same as accuracy at transported proposal states.

I therefore screened exact collapsed proposals before considering a code
change. Both fits used the same saved metric, the current eight-coordinate
hard-bound transport mask, and 20 proposals from the same checkpoint state.
The unregularized fit had mean exact acceptance `1.1e-5`, median log
correction `-18.27`, and acceptance-weighted group movement per second
`2.8e-6`. Ridge 0.1 had mean acceptance `1.7e-6`, median log correction
`-21.07`, and movement per second `4.2e-7`. The regularized fit predicted
held-out posterior likelihoods somewhat better but corrected actual
transported proposals substantially worse. I rejected ridge fitting and kept
the production surrogate unchanged. These are short proposal screens, not a
whole-fit ESS comparison. Reproduction scripts and outputs are
`benchmark_surrogate_ridge_cv.R`, `surrogate_ridge_cv_summary.csv`,
`surrogate_ridge_cv_subject.csv`, `benchmark_surrogate_ridge_exact_screen.R`,
and `surrogate_ridge_exact_screen_summary.csv`.

I also checked what an observed-RT upper bound would actually do in the current
likelihood. `do_bound()` checks the model's global `minmax` matrix; for LBA
that gives `t0 >= 0.05` and no finite upper bound. The LBA kernel separately
checks each observed row's `rt - t0`. If that is nonpositive, the winner
density is zero and the trial likelihood is floored; a loser with negative
accumulation time contributes survivor probability one. The remaining trials
are still evaluated. For chain 3, subject 14, the saved state had log
likelihood 4.73. Adding 0.5 to each of its three `t0` coefficients made 120
of 448 accumulator rows exceed their corresponding RT and changed the log
likelihood to -2593. Adding 1 made 368 rows exceed RT and gave -4583, still
above the HMC floor of -4789. Adding 2 made all 208 unique observations
impossible and reached that floor. Thus one `t0 >= min(RT)` crossing does not
by itself make the current subject likelihood equal to its floor; the other
observations continue contributing.

A short replay of the stored full-vector HMC kernel for this same subject
identified the abort source. At its saved step size (`eps = 0.363`), 16 of 20
diagnostic transitions aborted. Every floor-valued candidate inspected had
`t0` below the model-wide lower bound 0.05 (mapped values ranged from about
`2e-8` to 0.0257); none had a winning `t0 >= RT` violation. The model-wide
`do_bound()` mask rejected 416 accumulator rows in these candidates. This
short replay is not a production abort-rate estimate, but it shows that the
observed-RT ceiling was not the cause of these HMC aborts: the lower bound was.
Because each HMC transition uses one trajectory over all coefficients, that
lower-bound failure rejects the free-coordinate movement from the same
trajectory as well.

A scalar subject-specific cap is safe only for designs where the same `t0`
value applies to every observation and where no likelihood component supports
responses before `t0`. In the saved LBA design, `t0` varies by condition, so
the correct support check is observation-specific; the coefficients are
contrasts and cannot each be assigned `min(RT)` as an independent maximum.
The likelihood already avoids the expensive LBA density calculation on an
individual pre-onset row. Adding a whole-subject cap would change the current
floored likelihood for cases where other observations remain supported. No
such cap or bypass was added.

### HMC abort-rate context, 2026-10-10

The saved sampler's aggregate subject-HMC counters put the earlier 16-aborts-
in-20 diagnostic in context. Across its production checkpoint, the three
chains recorded 333/81,087, 415/73,093, and 581/69,543 HMC aborts (0.41%,
0.57%, and 0.84%; 0.59% pooled). These are all HMC abort causes, not just
support-bound failures. A separate 100-transition replay from chain 3,
subject 14, at its saved step size (`eps = 0.299`) had no aborts. That single
replay is not a rate estimate; it confirms that the earlier short replay was
state- and stream-sensitive. The aggregate counters do not support a default
HMC block split or a global step-size reduction as a package-wide fix. The
probe-only finite-difference fallback remains a possible robustness change,
but this audit did not establish a meaningful production cost to recover.

### Reuse the conditional Gaussian factorization, 2026-10-10

Equal-step microbenchmarks showed that the conditioned collapsed trajectory
spent roughly 3--5 times as long per trajectory as the unconditional path on
the archived `fullc` checkpoints. Those runs are not ESS comparisons: at a
shared step size, acceptance differed substantially for T30, and some arms
were measured only briefly. Source inspection found redundant per-subject
factorizations in the conditional target: it factored the full precision, the
free-coordinate precision, and the fixed-coordinate Schur complement at each
gradient evaluation.

The implementation now orders free coordinates before fixed coordinates and
uses one Cholesky factorization. Its leading block gives the free-coordinate
conditional covariance and mean adjustment; its trailing block factors the
fixed-coordinate marginal precision. This is an algebraic refactoring of the
same Gaussian terms, with no change to the conditioned target. To reduce
additional allocations, the sufficient-statistic accumulator now adds only
the free-coordinate covariance block instead of constructing a full zero-padded
matrix for every subject. The focused collapsed-move and move-allocation suites
passed after these changes: 779 assertions, no failures or warnings.

The historical all-hard-bound mask microbenchmark measured the cost of one C++ surrogate-gradient
evaluation on archived fixture states. It alternated 250 evaluations per arm
over six repetitions; the unconditioned arm uses the original full Gaussian
subject target and is a timing reference, not an equivalent support-safe
kernel.

| Fixture | Fixed subject coefficients | Conditioned, ms/evaluation | Unconditioned, ms/evaluation | Ratio |
|---|---:|---:|---:|---:|
| T30 | 7 | 0.132 | 0.048 | 2.75 |
| corr15 | 7 | 0.100 | 0.044 | 2.27 |
| Forstmann | 1 | 0.024 | 0.010 | 2.40 |
| LBA | 8 | 0.436 | 0.140 | 3.11 |

The block factorization removes redundant decompositions, but this historical
conditioned path made each gradient evaluation roughly 2.3--3.1 times as expensive as
the unconditioned reference. Absolute cost ranged from 0.024 to 0.436 ms per
evaluation. The benchmark does not establish whole-transition or ESS-per-
second gains.

I also screened 64 one-round proposals from each saved checkpoint using the
same metric, step size, starting state and per-proposal random seed across
three masks. The movement score is exact-acceptance-weighted, C-whitened group
ESJD per second; it includes measured context and state setup per proposal.
These are local kernel screens, not stationary ESS estimates.

| Fixture | Mask | Expected exact acceptance | Group ESJD/s | Any global-bound hit | Any `t0 >= RT` hit |
|---|---|---:|---:|---:|---:|
| T30 | Unconditioned | 0.378 | 1,989 | 0/64 | 0/64 |
| T30 | `t0` only | 0.585 | 2,899 | 0/64 | 0/64 |
| T30 | All hard bounds | 0.427 | 1,958 | 0/64 | 0/64 |
| corr15 | Unconditioned | 0.747 | 6,927 | 0/64 | 0/64 |
| corr15 | `t0` only | 0.883 | 7,408 | 0/64 | 0/64 |
| corr15 | All hard bounds | 0.984 | 7,234 | 0/64 | 0/64 |
| Forstmann | Unconditioned | 0.105 | 453 | 0/64 | 3/64 |
| Forstmann | `t0` only | 0.750 | 2,547 | 0/64 | 0/64 |
| Forstmann | All hard bounds | 0.750 | 1,211 | 0/64 | 0/64 |
| LBA | Unconditioned | 0 | 0 | 64/64 | 0/64 |
| LBA | `t0` only | 7.2e-8 | 3.1e-5 | 0/64 | 0/64 |
| LBA | All hard bounds | 0.0063 | 13.3 | 0/64 | 0/64 |

The package's own `do_bound()` semantics were used to count bound hits, and
the mapping needed for that audit ran outside the timed section. T30 and
corr15 had no global-bound crossings under any mask; freezing all their
hard-bounded coefficients did not prevent an observed crossing in this screen.
The `t0`-only mask improved ESJD/s in all three `fullc` fixtures, while adding
all hard-bound coefficients was mixed. For the LBA state, the all-hard mask
had only a 2.2% inner acceptance and zero realized exact moves in 64 attempts;
its positive expected score is too small to establish useful sampling. This
does not justify the all-hard-bound mask. Matched sampling continuations have
since compared all-hard, t0-only and unconditioned transport on ESS per wall
second.

The historical timing and mask screens are reproduced by
`benchmark_conditional_gradient_cost.R` and
`benchmark_conditional_mask_efficiency.R`; results are in
`conditional_gradient_cost_all_hard_bounds.csv` and
`conditional_mask_efficiency_historical.csv`.

### Matched mask continuations and current production choice, 2026-10-10

The earlier continuation note is now resolved. I ran three archived `fullc`
checkpoints for 200 new draws per chain in three chains, starting each arm
from the same checkpoint and resetting the same RNG seed. Each arm used one
collapsed round per transition and the archived metric and step size. Subject
HMC was disabled to isolate the collapsed-mask effect; other sampler moves
remained active. No scalar scale multiplier or exact-pilot calibration was
used. These short runs compare kernels, not final convergence.

| Fixture | Mask | Group-mean ESS/s (median) | `t0` group-mean ESS/s | Subject ESS/s (median) | Exact acceptance | Seconds | Max group-mean R-hat |
|---|---|---:|---:|---:|---:|---:|---:|
| T30 | Unconditioned | 0.940 | 1.182 | 1.098 | 0.217 | 117.6 | 1.052 |
| T30 | `t0` only | 1.260 | 0.203 | 1.410 | 0.395 | 127.2 | 1.099 |
| T30 | All hard bounds | 0.085 | 0.085 | 0.265 | 0.175 | 113.4 | 2.568 |
| corr15 | Unconditioned | 5.918 | 5.918 | 5.554 | 0.615 | 101.2 | 1.006 |
| corr15 | `t0` only | 5.496 | 5.506 | 5.207 | 0.764 | 107.8 | 1.004 |
| corr15 | All hard bounds | 5.321 | 5.871 | 5.187 | 0.895 | 104.4 | 1.007 |
| Forstmann | Unconditioned | 12.416 | 12.573 | 12.400 | 0.373 | 46.7 | 1.008 |
| Forstmann | `t0` only | 11.139 | 12.002 | 11.128 | 0.916 | 52.2 | 1.010 |

All-hard conditioning was clearly harmful on T30: it held seven coefficients
fixed instead of one `t0` coefficient and had poor R-hat after this short
run. The narrower `t0` mask is mixed against unconditioned transport. It
improved median group-mean and subject ESS/s on T30, but reduced `t0` group
mean ESS/s by about 83%. Unconditioned transport was about 7% faster on
corr15 and 10% faster on Forstmann by the median group and subject ESS/s,
despite lower exact acceptance. In Forstmann, conditioning raised exact
acceptance from 0.37 to 0.92 while reducing group and subject ESS/s. Acceptance
alone is not an efficiency measure.

The production mask now holds only `t0` subject coefficients fixed. This is a
support-specific compromise: it avoids tying unrelated bounded coordinates
into the collapsed conditional target while retaining the exact correction
for other bounds. The measurements do not show that this mask is faster for
every fit. The chain-local allocator still determines whether and how often
the collapsed move runs. Scale calibration remains removed because these
matched runs did not establish a package-wide ESS/s gain from it.

After narrowing the mask, a repeated C++ gradient microbenchmark (250
evaluations per arm, six repetitions) measured 0.162/0.044 ms for T30,
0.124/0.044 ms for corr15, 0.024/0.008 ms for Forstmann, and 0.540/0.140 ms
for LBA (conditioned/unconditioned). Thus the current conditional gradient
cost is 2.8--3.9 times the unconditioned reference, although its absolute
cost remains below 0.6 ms per evaluation.

The current gradient results are in `conditional_gradient_cost_current.csv`;
the historical all-hard results are in
`conditional_gradient_cost_all_hard_bounds.csv`. The local mask screen is in
`conditional_mask_efficiency_historical.csv`. Paired ESS outputs are in
`t0_vs_all_bounds_ess_t30.csv`, `t0_vs_all_bounds_ess_corr15.csv`,
`t0_vs_unconditioned_ess_t30.csv`, `t0_vs_unconditioned_ess_corr15.csv`, and
`t0_vs_unconditioned_ess_forstmann.csv`.

### Multi-start, all-parameter ESS and production allocation, 2026-10-10

The earlier matched-mask table summarized selected parameter families from
one archived start per fixture. I extended it to all 13 group means, all 13
group log-scales, and every subject coefficient, and repeated the comparison
from all three archived fullc starts for T30, corr15, and Forstmann. Each pair
started from the same three chains, used the same random seed, ran 200 draws
per chain, retained the archived covariance metric and step size, disabled
subject HMC, and forced one collapsed round per transition. This isolates the
mask effect with other sampler moves left active. Per-variable outputs are
in `t0_vs_unconditioned_all_parameter_ess_*_run*_variables.csv`.

The table reports median ratios across each fixture's three starts. Ratios
above one favor `t0` conditioning. `Sum ESS/s` sums bulk ESS per wall second
over the individual variables in each run before taking the median ratio.

| Fixture | All-variable median ESS/s | Sum ESS/s | 10th-percentile ESS/s | `t0` group-mean ESS/s |
|---|---:|---:|---:|---:|
| T30 | 1.335 | 1.280 | 1.039 | 0.225 |
| corr15 | 0.973 | 0.964 | 0.955 | 0.964 |
| Forstmann | 0.943 | 0.923 | 0.904 | 0.986 |

Across all nine paired starts, the geometric mean ratio was 1.068 for
all-variable median ESS/s, 1.052 for summed ESS/s, and 0.993 for 10th-percentile
ESS/s. The strong T30 gain outweighs the small losses on corr15 and Forstmann
under equal fixture weighting. The tradeoff is concentrated in T30 `t0` group
means: their median ESS/s ratio was 0.225, while corr15 and Forstmann were
near parity. This supports the current mask as a package compromise for
typical movement, but it does not make it a per-parameter or per-model
optimum. Acceptance was higher under conditioning in every fixture family;
it was not sufficient to predict ESS/s.

I also reran production allocation from archived sample histories with
subject HMC enabled and the rebuild counter reset. The allocator selected the
same rates for both masks in each paired start. For T30, the three pairs
received 1.71, 1.25, and 1.54 collapsed rounds per transition respectively.
The subject-HMC rates were 0.0258 for both arms on the first start and zero
for both arms on the other two. The all-variable median ESS/s ratios for the
`t0` arm were 1.047, 0.887, and 1.348; the 10th-percentile ratios were 0.661,
0.427, and 1.144. Thus the production result is variable across starts, and
the mask did not change HMC or collapsed-round allocation in these tests.
For corr15 and Forstmann, production allocation selected zero collapsed rounds
and zero subject-HMC rate for both masks; their resulting ESS/s differences
were within timing variation.

The production-mode runs are in
`t0_vs_unconditioned_all_parameter_ess_t30_run*_production.csv`,
`t0_vs_unconditioned_all_parameter_ess_corr15_run1_production.csv`, and
`t0_vs_unconditioned_all_parameter_ess_forstmann_run2_production.csv`.
The reproduction script now supports both controlled one-round comparisons
and allocator-rebuilt production comparisons through `EMC2_ESS_PRODUCTION`.
The evidence supports retaining the unscaled covariance metric and the
support-aware `t0` mask as the current general default, with the measured
per-model and per-parameter tradeoffs documented above.

### Conditional-gradient cost after factor reuse, 2026-10-10

The previous gradient timings above used a general conditioned implementation
that repeated matrix work at every gradient evaluation. I replaced it with a
rank-one update for one fixed coordinate and a rank-k conditional covariance
update for multiple fixed coordinates. Both use the Cholesky factor and
inverse-lower factor already computed for the collapsed likelihood. The
unused general implementation was removed. The test now checks the exact
target, finite-difference gradient, and transport for three fixed coordinates,
matching the three `t0` contrasts in the LBA case; the single-coordinate test
continues to cover the usual fixture path.

For a less timer-quantized measurement, this benchmark used 2,000 gradient
evaluations per arm and six alternating repetitions. It measures only the C++
surrogate gradient, not whole transitions or ESS per second.

| Fixture | Fixed coordinates | Conditioned, ms/evaluation | Unconditioned, ms/evaluation | Ratio |
|---|---:|---:|---:|---:|
| T30 | 1 | 0.0518 | 0.0428 | 1.21 |
| corr15 | 1 | 0.0508 | 0.0443 | 1.15 |
| Forstmann | 1 | 0.0105 | 0.0095 | 1.11 |
| LBA | 3 | 0.3673 | 0.1440 | 2.55 |

This reduces the one-coordinate overhead from the earlier 2.8--3.9x range to
about 11--21% on these fixtures. The three-coordinate LBA path is still more
expensive, though its measured cost fell from 0.552 to 0.367 ms/evaluation.
The saved LBA checkpoint receives no collapsed rounds from the current
chain-local allocator, so this microbenchmark does not imply a production
runtime change for that checkpoint.

The paired 200-draw-per-chain continuations after the rank-one optimization
still show why acceptance cannot select the default mask by itself:

| Fixture | Mask | Group-mean ESS/s | `t0` ESS/s | Subject ESS/s | Acceptance | Seconds |
|---|---|---:|---:|---:|---:|---:|
| T30 | Unconditioned | 0.925 | 1.163 | 1.080 | 0.217 | 119.5 |
| T30 | `t0` only | 1.352 | 0.260 | 1.471 | 0.395 | 124.6 |
| corr15 | Unconditioned | 5.929 | 5.929 | 5.565 | 0.615 | 101.0 |
| corr15 | `t0` only | 5.688 | 5.697 | 5.388 | 0.764 | 104.2 |
| Forstmann | Unconditioned | 12.396 | 12.552 | 12.380 | 0.373 | 46.8 |
| Forstmann | `t0` only | 11.516 | 12.409 | 11.505 | 0.916 | 50.5 |

Conditioning remains a mixed package-wide choice: it improves T30 group and
subject ESS/s but greatly reduces T30 `t0` ESS/s; unconditioned transport is
slightly more efficient on corr15 and Forstmann despite lower acceptance.
The timings do not justify a universal scale reduction or claim that the
`t0` mask maximizes ESS/s for every model. The rank-one and rank-k changes
reduce the cost of the current conditional implementation without changing
its target or proposal scale.

The updated microbenchmark is in `conditional_gradient_cost_current.csv`; the
pre-rank-k comparison is retained in `conditional_gradient_cost_rank1.csv`.
The matched continuation summaries are in
`t0_vs_unconditioned_ess_t30.csv`, `t0_vs_unconditioned_ess_corr15.csv`, and
`t0_vs_unconditioned_ess_forstmann.csv`.

### Support-mask selectivity check, 2026-10-10

I checked whether the blanket `t0` mask could be limited to subjects whose
current values are close to a support boundary. At the archived fullc states,
the minimum subject `t0` values were 0.208 s for T30, 0.230 s for corr15, and
0.147 s for Forstmann; none of those subjects had a minimum below 0.10 s. In
the LBA checkpoint, 56 of 93 subjects were below 0.10 s, with the minimum at
0.0509 s. This suggests that a lower-bound rule would leave those three fullc
fits unconditioned while selecting many LBA subjects.

I screened 32 unconditioned proposals at each fullc checkpoint, suppressing
exact likelihood evaluation and mapping the proposed subject values only.
There were no lower-bound or observed-RT crossings for T30 (31 proposals
reached the transport) or corr15 (32 proposals); Forstmann had two proposals
with an observed-RT crossing. An earlier 64-proposal LBA screen found support
risk spread across 38 subjects. Freezing subjects whose current minimum `t0`
was below 0.10 s removed 99.8% of held-out subject-level crossings, but two
proposals still had a crossing. These are local screens from one state per
fixture; they do not validate an adaptive subject mask or estimate stationary
ESS.

I did not implement that mask. The 0.10 s rule is tied to the LBA's 0.05 s
lower bound, misses the observed-RT risk seen in Forstmann, and has no
unit-invariant interpretation across models. A principled risk score would
need to use each model's design and subject-frame covariance. The current C++
target and R transport accept one shared set of fixed coordinates for every
subject; a varying mask would change the Gaussian conditional terms and
transport per subject. The selected subject set would also need to stay fixed
for the production kernel: choosing it anew from each current draw makes the
kernel state-dependent and would require a proposal-selection correction.
This needs correctness and ESS-per-second evidence before becoming package
behavior. In particular, current-source LBA allocation selects zero collapsed
rounds, so this screen does not show a production efficiency gain there.

This also answers the design-bound question for this path: `do_bound()` uses
model-wide bounds and does not infer an upper bound from RT. Its `ok` flag
checks the declared lower and upper model bounds; it does not include the
data-dependent `t0 < rt` condition. The local `bad` in `.rfun_LBA()` is for
simulation only: it marks trials where every accumulator finish time is
infinite. In the LBA likelihood path, `t0 >= rt` gives that
trial a zero raw density, but the likelihood then floors the per-trial log
likelihood at `min_ll` (normally `log(1e-10)`). It is therefore a severe
downweighting rather than a hard rejection. Marking the entire subject invalid
would also floor its other trials and change the target; a hard design-time
`t0` ceiling would remove parameter values that currently retain finite
likelihood, and contamination/censoring can further change the semantics.
`design()` could compute per-subject/per-condition RT summaries, but the
current global model bound cannot represent those as ordinary parameter
bounds. I did not add a dynamic RT ceiling or subject mask.

The allocator records collapsed movement as expected exact acceptance times
the squared chart jump, and records the measured CPU cost per collapsed round.
It solves for collapsed and subject-HMC rates together. Near-zero collapsed
acceptance therefore makes that move look unproductive, while its measured
cost makes the allocator less likely to schedule it at the next rebuild. A
rate already selected continues to pay its exact-likelihood cost until that
rebuild. This can shift the HMC rate indirectly through the joint allocation,
but it does not directly change HMC's own acceptance or step-size adaptation;
HMC movement and cost have separate measurements. The saved LBA checkpoint's
current production allocation selected zero collapsed rounds, so its rejected
collapsed proposals were not consuming production likelihood calls there.

### Kernel-mixture and focused t0-chart follow-up, 2026-10-10

I tested two ways to recover movement in the t0 group parameters without
shrinking every collapsed chart coordinate. On T30, I compared the existing
t0-conditioned kernel with state-independent random mixtures that selected an
unconditioned kernel on 20% or 10% of collapsed rounds. Each comparison used
the same starts, one collapsed round per transition, three chains, 200 draws
per chain, and no subject HMC. The 10% mixture raised t0 group-mean ESS/s from
0.148 to 0.431, while all-variable median ESS/s fell from 1.544 to 1.372 and
the 10th percentile fell from 0.610 to 0.577. The 20% mixture reached t0
group-mean ESS/s 0.394, all-variable median 1.268, and 10th percentile 0.488.
Thus the mixture recovered some t0 movement but reduced overall efficiency in
this short T30 continuation. Acceptance alone would have favored the mixtures
less clearly than the parameter-specific ESS results.

I also tested a focused chart update of the t0 group mean and log-scale. The
first prototype incorrectly took a submatrix of the full Cholesky factor;
that is a conditional metric, not the marginal metric for the selected chart
coordinates. I corrected the screen to form the selected covariance block
from the full covariance and factor that block. The corrected 200-draw
continuation still performed poorly: all-variable median ESS/s was 0.353,
10th-percentile ESS/s 0.111, t0 group-mean ESS/s 0.175, and maximum group-mean
R-hat 1.354. I discarded this block design and left no block-specific
collapsed implementation in package source.

These are short, single-fixture kernel comparisons, not a package-wide
sampler result. They do not justify adding a kernel-mixture option or a
focused t0 chart block. The package source contains no scale-search,
directional-search, mixture, or chart-block prototype; the temporary scripts
and CSVs for these screens were kept outside the repository and removed after
recording the results here.

### Fresh-stream repeat of the forced LBA split-scale move, 2026-10-10

To check whether a four-round split-t0 move could be worthwhile when the
allocator selected zero rounds, I repeated the forced comparison twice from
the same saved three-chain LBA start with fresh, paired RNG streams. Both runs
used current source, 200 draws per chain, subject HMC off, and the same
candidate metric: global scale 0.03125, t0-mean rows at half that scale,
t0-hierarchy rows at one quarter, with four collapsed rounds per iteration
on chain 2 only. The retained kernel holds the three subject t0 coefficients
fixed and uses the current 1.5 trajectory length. Candidate wall time includes
its 40-trajectory surrogate step-size tuning (0.035 and 0.036 seconds); the
starting metric construction is common to each pair.

| Fresh stream | Arm | Wall seconds | Expected exact acceptance | All-variable median ESS/s | 10th percentile ESS/s | Group-mean median ESS/s | t0-mean median ESS/s | Subject median ESS/s | Sum ESS/s | Max R-hat |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1, candidate first | No collapsed | 380.0 | — | 0.854 | 0.271 | 0.181 | 1.000 | 0.859 | 1,456 | 1.295 |
| 1, candidate first | Split t0 | 349.5 | 0.637 | 0.909 | 0.277 | 0.211 | 1.154 | 0.921 | 1,561 | 1.303 |
| 2, control first | No collapsed | 320.8 | — | 0.495 | 0.033 | 0.044 | 0.960 | 0.505 | 1,091 | 1.466 |
| 2, control first | Split t0 | 407.8 | 0.647 | 0.417 | 0.028 | 0.046 | 0.794 | 0.431 | 900 | 1.444 |

The wall-time result changed sign when arm order was reversed. The first
candidate was 8% faster; the second was 27% slower. The summed raw ESS was
1.4% lower for the candidate in the first pair and 4.9% higher in the second.
Across both pairs, the forced move used 8% more wall time for about 1% more
summed raw ESS, so combined summed ESS per wall second was about 6.5% lower.
The high expected exact acceptance (0.64--0.65) did not indicate a reliable
efficiency gain. The second pair's maximum R-hat exceeded 1.4 in both arms, and the first
pair's exceeded 1.29. These 200-draw continuations are too short to resolve
the remaining parameter-level differences, but they do show that 0.64
acceptance does not reliably buy better wall-time per ESS. The controlled
scale setting remains a rejected research candidate, and no scale or
allocation change is justified by these runs.

### Allocator objective review, 2026-10-10

The allocation solver minimizes the largest predicted coordinate-wise
`tau * cost`, which is approximately maximizing the minimum coordinate ESS
per second. It does not maximize median or summed ESS per second. If a
coordinate that an optional move cannot improve sets the maximum, the solver
can rationally choose zero rounds even when the move would improve many other
coordinates. The fixed-coordinate allocator test covers this floor behavior;
it is a consequence of the chosen objective, not an accidental scale tie
between parameters.

On the saved LBA state, the worst measured group coordinate had a one-round
collapsed factor near 0.998 and a shear factor near 0.95. The model therefore
prices collapsed rounds as ineffective for its bottleneck. Short forced
continuations sometimes improved group-mean ESS, but fresh paired repeats
showed lower combined summed ESS per second. This is not enough evidence to
replace the conservative objective with a median or sum objective: those
objectives make different tradeoffs and need a package-wide target and
replicated benchmarks. The current chain-local choices also do not jointly
optimize a shared parallel wall-time objective across chains. No allocator
policy change is retained from this review.

The source tree has no exact-pilot scale search, directional scale search,
mixture kernel, focused chart block, or data-derived hard RT cap. The
production `t0` safeguard remains conditional transport of subject `t0`
coefficients, while exact likelihood evaluation enforces the model-specific
support. A simple per-subject `t0 < min(rt)` hard bound would not be a safe
replacement for that behavior in designs with varying `t0`, censoring, or
contaminant likelihood components.

### Inner trajectory-count screen, 2026-10-10

I screened 1, 2, 3, 4, and 6 inner surrogate trajectories per exact correction
from the saved LBA state, with the current full covariance metric, trajectory
length 1.5, and subject `t0` coefficients fixed. Each setting generated 24
proposals in two repeats. Exact acceptance stayed between about `1e-22` and
`1e-31` in the first repeat and between `1e-21` and `1e-25` in the second;
acceptance-weighted group movement was negligible for every setting. Timing
was too noisy to distinguish trajectory counts, and changing the count did
not repair the surrogate correction. I retained the default of three; fewer
trajectories are not a fix for the LBA bottleneck, and more would add inner
work without evidence of better exact moves.

A separate screen that sums stored per-subject correction terms suggests
that correction error accumulates with the number of transported subjects.
That calculation is not a valid benchmark of a subject-blocked kernel: a
blocked kernel changes the conditional group target when the other subjects
are held fixed. Such a kernel would need a new per-subject target and exact
likelihood dispatch path, then matched ESS-per-second comparisons across the
`fullc` fixtures. I did not add that complexity based on the subset screen.
The current LBA checkpoint receives zero collapsed rounds from the allocator
at its final saved allocation, avoiding exact proposals once the move is
measured as unproductive.

### Repeated unmeasured allocation probes, 2026-10-10

I inspected the allocation history in the saved 93-subject LBA checkpoint.
Across 800 transitions, its three chains attempted 487, 491, and 484
collapsed rounds and accepted none. The allocation history shows repeated
sub-unit rates during burn and adaptation, followed by zero in the later
production blocks. Outside production, proposal rebuilds intentionally use
only the just-finished block because the surrogate and metric change. When a
block runs no collapsed moves, that next block has no move measurement and
the allocator returns to its optimistic prior factor; the move can then be
turned on again. So the final zero rate did not prevent substantial earlier
waste. At 93 subjects, each rejected round evaluates the exact subject
likelihood for all 93 subjects.

The move does not directly change subject-HMC acceptance or step-size
adaptation, but the shared allocator solves for HMC and collapsed rates
together. An optimistic unmeasured collapsed score can therefore temporarily
compete with HMC allocation. I initially proposed capping its first measured
window to ten exact corrections. That cap was provisional and was later
removed after the matched T30 `fullc` comparison below showed worse short-run
ESS and wall time.

### Follow-up: bounded probes versus retained movement evidence, 2026-10-10

I tested the ten-correction cap on the archived T30 `fullc` checkpoint, with
subject HMC disabled and a matched 200-transition continuation from the same
three-chain start. The checkpoint stores its rebuild schedule, so a session
option did not shorten the first run's proposal window. After setting the
checkpoint's stored schedule explicitly to two 100-transition windows, the
capped run allocated 0.1 collapsed rounds/transition in its first window,
made about ten proposals per chain, and then raised the next rates to 0.25--0.375
after measuring acceptance. The uncapped arm started at 0.625--0.75 and
continued at 0.625--1.375 after its first measurement.

| Arm | Wall seconds | Median group-mean ESS | `t0` group-mean ESS | First-window rate | Second-window rate |
|---|---:|---:|---:|---:|---:|
| Ten-proposal cap | 60.86 | 92.61 | 41.76 | 0.1 | 0.25--0.375 |
| Uncapped | 50.14 | 181.54 | 83.79 | 0.625--0.75 | 0.625--1.375 |

This is a short, single-start comparison, but it gives no support for a
package-wide cap: the capped arm took 21% longer and had about half the median
group and `t0` ESS. I removed the cap and its test rather than keep a policy
that can hold a useful unmeasured kernel below its productive rate.

The more targeted allocator failure is that an empty burn/adapt block erases a
previously measured collapsed score. The saved LBA histories show a repeated
pattern: a nonzero rate with the prior factor 0.5, a measured zero-movement
factor of 1, then another prior-factor allocation after a block with no moves.
I changed `.emc_collapsed_set_plan()` to retain the last measured block when a
same-stage rebuild has fewer than five new proposals. It does not pool
different nonempty burn/adapt blocks, and a stage change starts fresh; sample
stage measurements continue to use their existing pool. This preserves the
full initial allocation for an unmeasured useful move while preventing an
empty block from erasing a measured failure. A focused test covers both the
same-stage carry-forward and stage reset. The focused R test suite passes 783
assertions with this change.

The `bad` flag check also shows why adding an RT-derived upper parameter bound
at `design()` is not just a compute shortcut. The current LBA likelihood
returns a floored per-trial contribution when `t0` exceeds an observed RT. A
hard ceiling would change those finite contributions and thus the target.

### Unit-scale and subject-block follow-up, 2026-10-10

I directly checked whether conditioning the subject `t0` coefficients was
enough to use the original covariance factor on the LBA checkpoint. At unit
scale, the 20-proposal screen had no explicit-bound or RT-ceiling violations,
but its mean exact acceptance probability was only `4.9e-34`. At scale
`0.25`, it was `4.6e-5`. This separates support safety from surrogate accuracy:
fixing the `t0` values removes those boundary failures, but it does not make
the full 93-subject exact correction usable at the original transport size.
The scale reduction was responding to accumulated surrogate error as well as
support crossings.

I then screened a state-independent block of 16 subjects per collapsed round.
The other subjects' current coefficients were included in the conditional
group target, and only the selected subjects were transported and evaluated
by the exact likelihood. This is a valid conditional kernel, but its inner
target still depends on all subjects. On 16 local proposals per chain with
the original covariance factor, expected exact acceptance ranged from `0.48`
to `0.87`; retuned step sizes were only `0.09`–`0.10`, so the higher
acceptance also required more leapfrog steps.

For an end-to-end check, I randomized the 16-subject block each transition,
used the original covariance factor, retuned the step size once per chain,
disabled subject HMC, and continued the same saved three-chain checkpoint for
100 draws per chain. The no-collapsed arm took 149.9 seconds. The block arm
took 180.8 seconds plus 1.45 seconds of step-size tuning, accepted with mean
probability `0.468`, and moved on 48% of proposals. Yet median group-mean
ESS/sec fell from `0.431` to `0.123` (−71%), median ESS/sec across all group
and subject variables fell from `1.409` to `0.852` (−40%), and maximum group
mean R-hat rose from `1.077` to `1.179`. This is one short continuation and
not a precise efficiency estimate, but it gives no support for the added
conditional kernel. Its prototype, worker-index plumbing, and prototype-only
tests were removed from package source.

The package keeps conditional transport for bounded subject coefficients,
which prevents those coefficients from crossing their support during a
collapsed move. It does not apply a global scale multiplier or run the
subject-block prototype. The original-scale LBA move remains unusable after
support crossings are removed, and the current allocator assigns it zero
rounds on the saved checkpoint. No general package-wide scale or block policy
has passed an ESS-per-second comparison; the surviving changes address the
verified surrogate-fit and allocation bookkeeping errors without claiming to
solve the remaining LBA ridge.

### Surrogate-residual preflight screen, 2026-10-10 (rejected)

I tested whether stored exact subject likelihoods and the already-fitted
quadratic surrogates can provide a zero-new-likelihood warning before a
collapsed move is allocated. Across the saved `fullc` fixtures, the
in-sample total exact-minus-surrogate residual had SD 1.0--1.7; for the LBA
checkpoint it was 18--33. A leave-one-chain-out version showed the same broad
separation, but refitting all subject surrogates took 0.3--1.7 seconds and
would not be a cost-free rebuild check.

Neither simple score is reliable enough to drive allocation. Treating two
residual draws as independent predicted correction acceptance of about
0.64--0.71 for the fixtures, but also 0.51--0.52 for LBA, whose direct
support-conditioned proposal screens had expected acceptance near zero.
The importance-weight ESS ratio was below 0.01 for some fixture chains that
have productive proposals, so a cutoff on that score could suppress useful
moves. The observed residual spread is a useful audit clue for this set of
checkpoints, but it has not been validated as a package-wide predictor of
proposal acceptance or ESS per second. No preflight score, threshold, or
allocation policy was added to package code.
