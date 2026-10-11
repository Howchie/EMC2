# Subject HMC, shear rounds and the shared move allocation

Implements `plans/subject_hmc_plan_2026-10-07.md`, with the deviations listed
under Design. Code: `R/subject_hmc.R` (kernel, warmup, step-size updates),
`R/shear_move.R` (shear rounds), `R/move_allocation.R` (measurement and the
solver), wired into `run_stage()` (`R/sampling.R`), the worker pool
(`R/chain_pool.R`) and `add_proposals()` (`R/fitting.R`). Committed in
`e43fdc30` on `playground`. The collapsed transport is enabled by default and
can be disabled with `options(emc2.collapsed_move = FALSE)`.

## Design

- Iteration order: group Gibbs, particle step, ensemble group update, shear
  rounds, subject HMC (its own worker round), pooled parameters.
- Shear round: one Tjelmeland multiple-proposal step with 8 candidates along
  the orbit t -> T_t(x) of one group coordinate (a mean direction or a log
  group sd), subjects transported in a non-centred conditional frame, all
  candidates in one likelihood round.
- Subject HMC: whitened by the subject's likelihood fit plus the group prior;
  forward-difference gradient (p + 1 rows per gradient), corrected exactly by
  the Hamiltonian acceptance.
- Allocation: at each proposal rebuild every coordinate's lag-one
  autocorrelation over the rest of the iteration is measured as
  rho = 1 - ESJD / (2 * pooled variance across chains), a lower bound. One
  solver chooses expected shear rounds per group coordinate (cap 4) and HMC
  rates per subject to minimise the slowest coordinate's tau times CPU per
  iteration, buying nothing unless the predicted gain is at least 10%. x
  rounds multiply a coordinate's rho by g^x, with g the measured factor of
  one round.
- HMC warmup: in a subject's first burn block every eligible subject moves at
  rate 0.5 while its step size adapts by Robbins-Monro towards acceptance
  0.75, for 20 transitions. Its draws are discarded with burn. Later step
  sizes change only at rebuilds.

Deviations from the plan: HMC moved last so it is measured against the whole
of the rest of the iteration (measured against the particle step alone it
bought rows that the ensemble already provided); one allocation for both
moves in place of separate HMC rates; up-front calibration replaced by the
warmup; forward in place of central differences.

## Results

Default `fit()` stop rule (mu and alpha R-hat < 1.1, at most 20 sample
tries), cgroup CPU in seconds. `base` is the ensemble update alone; `full`
adds shear rounds and subject HMC. Neither arm enables the collapsed group
transport; only the separate `fullc` arm does. The reported `full` column
pools runs whose builds differ from the commit in one of: the saturating round
composition, up-front calibration in place of the warmup, or central
differences.

| case | base | full |
|---|---|---|
| forstmann (LNR) | 574, 567 | 561, 560 |
| corr15 (RDM, 30 subj) | 1465, 1463 | 1516-1579 (4 runs) |
| T30 (Wald, 40 subj) | fails both: mu R-hat 1.24, 1.15 | converges in all 12 runs: 10-19 tries, 2516-4038 s |
| lbaA (LBA, 20 subj) | 4811 (10 tries), 7193 (15 tries) | met in 2 of 4: 5404-6053 s; 2 runs hit 20+ tries |

These fixture results do not measure collapsed-transport acceptance or
efficiency. The reported HMC acceptance belongs to the subject-level HMC
kernel. Direct collapsed-move tests use small synthetic Gaussian likelihoods;
the large nonlinear LBA support case is covered by the checkpoint probes in
`plans/collapsed_transport_adaptation_2026-10-09/report.md`.

- T30 is the case the shear exists for: the ensemble update does not mix mu
  there, and shear mean rounds do.
- Tries on T30 and lbaA are held open by alpha R-hat. The slow alphas are
  coupled to the slow group log sds, and 85-97% of the rounds bought in
  production go to log-sd coordinates. A log-sd round moves only within one
  coordinate's conditional slice, which holds 5-30% of that coordinate's
  variance, so log-sd rounds help (min log-sd ESS roughly 3x over
  mean-only rounds) and then stop helping. The remaining slowness is this
  log-sd funnel, which needs a move on the covariance as a whole
  (`plans/lsd_funnel_handoff_2026-10-08.md`).
- Round composition: a saturating model (each coordinate's gain capped at
  its slice share, with a fixed per-round autocorrelation) was tried and
  removed. Over 5 T30 seeds per arm it matched the geometric composition on
  total CPU (3306 vs 3295 s) while hardcoding a saturation point that
  differs by model.
- Forward differences, benchmarked on the same states, step sizes and
  seeds: acceptance and accepted jump unchanged to 3 dp, half the rows,
  CPU per transition -7% (forstmann), -31% (T30), -40% (lbaA).

Driver: `run_arm.sh <lib> <case> <arm> <seed> <outdir> [cores]` runs `fit_arm.R` in a
cgroup scope; `summarise.R <dir>` tabulates a directory of runs. Build each
test library from a clean tree (`git archive`) so uncommitted work elsewhere
in the worktree is not picked up.

## Evaluation protocol

Single from-scratch fits do not rank sampler builds. Three fullc fits of the
same lbaT fixture under near-identical builds took 17.1k, 30.4k and 44.8k
CPU-seconds; two T30 seeds under one build took 8.9k and 5.1k. Convergence
time is a stopping time of a noisy R-hat and swings by 2-3x between seeds.
Continuations of 200-300 iterations are no better: on the N-back checkpoint,
arms at identical settings gave slowest-direction ESS between 13 and 73.
Comparisons are therefore made in two phases.

Phase A, convergence: `fit_arm.R`, the default stop rule, three seeds per
case. Reports whether the rule was met, the tries, and CPU to the stop.
This is the acceptance criterion (R-hat < 1.1 first), not the ranking
statistic. The seed-1 end state (`SAVE_EMC=1`) is the case's checkpoint.

Phase B, stationary efficiency: `continue_arm.R` continues the checkpoint
for a fixed number of production iterations, three seeds per arm, every arm
from the same checkpoint. The first 100 new iterations are discarded (one
block: a rebuild is forced at the start so each arm allocates under its own
kernel, and the first block's plan is the only one built from the
checkpoint's measurements). The score is ESS per cgroup CPU-second on
quantities fixed by the checkpoint: the group means and log sds,
standardised by their checkpoint sds and rotated to the directions of
largest lag-one autocorrelation (generalised eigenvectors of the lag-one
autocovariance against the covariance, chains pooled; cached in
`<ckpt>_dirs.rds`). The slowest of these directions is the headline, with
coordinate-wise mu, log-sd and alpha minima, split-half ESS as a
stationarity check, R-hat on the new draws, and the move totals. Runs are
long enough that the reference arm's slowest direction reaches ESS of
about 100 (relative error of the ESS estimate near 15%); on the N-back
checkpoint that is 1800 iterations.

Ranking rule: `summarise_cont.R <dir> <reference>` reports each arm's
median over seeds as a ratio to the reference and whether the seed ranges
are separated. A build is a gain only when the ratio exceeds 1 on the
slowest direction with separated ranges, and Phase A still meets the stop
rule.

Machine protocol: one fit at a time, 3 chains x 9 cores, pinned with
`taskset` (`queue.sh <jobs.txt>` runs a job list this way). The allocator
prices moves by wall-clock, so concurrent load changes plans, not only CPU
totals. A frozen copy of the installed library per launched run, since a
reinstall under a running fit changes the code its workers load.

Cases: T30 (Wald), lbaT (LBA ridges), and the N-back fixtures built by
`make_fixture_nback.R` from the user's fitted object (real data, written to
`FIXTURE_DIR` outside the repo): `nback` (93 subjects, p = 21, the minimum
target) and `nback31` (every third subject; same design and prior, a third
of the cost per iteration, used for screening before the full fixture).

## Fixture manifest

| Fixture | Source | md5 / spec |
|---|---|---|
| T30 | `fixtures/fixture_T30.rds`, copied from session scratch `gtest/fx_T30_s1.0.rds` (the earlier genuine-group-test fixture). Generator: `WorkingTests/generate_correlated_group_rdm_fixture.R` with A = 0 (Wald), `CORR_STRENGTH=1`, `CORR_SD_SCALE=1`, `CORR_N_SUBJECTS=40`, `CORR_TRIALS_PER_CELL=30` | `bbbadba905b9d5aeda9cc3972d994d2e` |
| corr15 | `benchmarks/mode-start/results_plausible/fixture_corr15.rds` (RDM, 30 subjects, 400 trials/cell) | `cb9c47fe214f52608e99b760522e30b5` |
| forstmann | package data, LNR `m ~ lM, s ~ 1, t0 ~ 1` | package data |
| lbaA | `fixtures/fixture_lbaA.rds` from `make_fixture.R` (below) | see below |
| lbaT | `fixtures/fixture_lbaT.rds` from `make_fixture_lbaT.R`: LBA `v ~ lM` with an exponential practice trend `v + w exp(-d_ed trial)`, `B ~ 1, A ~ 1, t0 ~ 1, sv ~ lM` (8 parameters), 80 subjects x 150 trials, group SD of the log rate 0.15, seed 20261011. Few trials per subject leave the subject posteriors on curved ridges and the group means identified by pooling. First fullc fit (seed 1): no convergence in 20 tries (1135 draws, 1547 s); `B`, `A`, `t0` means and log SDs at ESS 6-9, R-hat 1.4, trend parameters fine (ESS 44-128). | generated |
| lba71sim stuck state | `benchmarks/production-adaptation/results_shear/shear_v4/emc_lba71sim_doubling_shearK4_2.rds`, subject 11 / chain 1, 2000 production transitions, next rebuild 3200, `diagonal-gamma`, 20 subjects, p = 71 | `38c06da04e0c21fdda3539305f6c2fd9` |

### Identified free-A fixture (lbaA)

`make_fixture.R`: the lba71sim LBA design without the random factor F (17
parameters, 20 subjects). Screened with `screen_identification.R` (multi-start
BFGS from the truth and dispersed starts, observed information at the best
optimum with two finite-difference steps):

| variant | outcome |
|---|---|
| A mean 0.5, 200 trials/cell | rejected: 3/20 subjects with condition ~1e8, A SE 55-138 (log scale) |
| A mean 0.5, 500 trials/cell | rejected: max condition 9.9e8 |
| A mean 1.0, 500 trials/cell | rejected: subject 7 flat in A at its MLE (A SE 38) |
| A mean 1.0, 1000 trials/cell | 20/20 PD at 6 starts, but at 8 starts subject 7 again reached a flat ridge towards A -> 0 (condition 2.3e8) |

The generator therefore redraws any subject whose own multi-start MLE is not
identified (condition number <= 1e6 at both steps); redraw counts are saved in the
fixture. Locked: A mean 1.0, 1000 trials/cell, seed 20261007.

## Kernel correctness (stage 2)

`kernel_tests.R` (scratch copy of the package tests, larger samples):

- Gaussian force (central differences, as tested at the time): equals the
  analytic whitened gradient to 2e-11 (h = 1e-4); truncation error O(h^2),
  cancellation below 1e-6. The committed kernel uses forward differences,
  O(h); see Results for the matched benchmark.
- Energy error at fixed integration time: log-log slope 2.00 against eps.
- Reversibility: forward then momentum-reversed path returns to 2e-11.
- Invariance on a curved (banana) conditional from exact draws, 3 transitions
  each: KS p 0.90/0.54 (eps 0.5), 0.91/0.14 (eps 1.2); largest moment shift
  |z| 1.6/2.1. Broken control (acceptance without kinetic energy): KS p
  0.033/<0.001, shift |z| 4.8.
- Support boundary (truncated target): no draw leaves the support; KS p
  0.63/0.52; a start with a probe outside aborts without moving.
- Allocator vs dense rate-grid oracle, 60 random problems: never worse.

## Fixed-group replay: lba71sim subject 11 (stage 3)

Measured with the earlier kernel (central differences, HMC after the
particle step, up-front calibration). The step-size and escape findings
carry over; the calibration study below describes a scheme the warmup has
since replaced.

Group state fixed at chain 1's last draw; reference conditional from a
1400-transition HMC chain: A median -0.511, sd 0.167.

h sensitivity (`hsens.R`): relative force difference between adjacent h is
~1e-8 for h in [1e-5, 1e-4] at the stuck state, the Laplace mean, three
Laplace draws and a 2.5-sd outlier; 1e-2 -> 1e-3 differ by 7e-5 (truncation),
1e-6 -> 1e-7 by 6e-7 (cancellation). Default h = 1e-4.

Step-size grid (120 transitions each from the Laplace mean, L = ceil(pi/2/eps)):

| eps | acceptance | aborts | rows/transition | acc-weighted ESJD per 1000 rows |
|---|---|---|---|---|
| 0.08 | 0.99 | 0 | 3044 | 0.70 |
| 0.12 | 0.98 | 0 | 2078 | 1.00 |
| 0.17 | 0.96 | 0 | 1542 | 1.39 |
| 0.25 | 0.91 | 0 | 1106 | 1.80 |
| 0.35 | 0.33 | 17% | 851 | 0.90 |
| 0.50 | 0 | 68% | 644 | 0 |

The one-step bracket gives 0.172 from both the stuck state and the Laplace
mean (2 trials, 428 rows): within 1.3x of the grid optimum, on the safe side of
the acceptance cliff.

Escape from dispersed starts, 300 iterations, eps 0.2 (`escape.R`):

| start | particle only: A mean (sd) | particle + HMC each iteration: escape iteration, A mean (sd) |
|---|---|---|
| chain 1 stuck (A -0.215) | -0.266 (0.040), never crosses the median | 2, -0.538 (0.154) |
| chain 2 state (-0.700) | -0.535 (0.046) | 1, -0.517 (0.168) |
| chain 3 state (-0.623) | -0.749 (0.040) | 2, -0.537 (0.167) |
| Laplace + 2 sd draws | -0.243 (0.036), -0.417 (0.069), -0.243 (0.052) | 4, 7 (one started in the bulk); sd 0.147-0.179 |

The particle step alone keeps each start near where it began (sd 0.04-0.07
against 0.167). HMC every iteration reaches the bulk within 7 iterations from
every start and matches the reference width. HMC alone failed from one far-tail
start (acceptance 0 throughout), where the particle step's first move brought the
chain within HMC's reach: the composition matters. At rate 0.1, escape took
5-206 iterations.

Calibration (`calib_study.R`, from the stuck state):

- Step-size rule: the jitter-halves comparison with 40 attempts per window,
  started at the bracketed 0.172, moved to 0.258 (the grid optimum region)
  after 1-6 windows in all four runs. One run then took a second step to
  0.344, at the acceptance cliff (acceptance 0.33 there); the decreasing gain
  (x1.5, x1.33, x1.26, ...) moves it back down when the next window resolves
  the slope. The rule finds the right region but can overshoot a cliff by one
  step.
- Pilot batch size: the allocation's measure, the largest per-coordinate
  rho_H over 71 coordinates, is biased upward in small batches (the maximum
  of noisy estimates): 0.78 (sd 0.08) at 10 pilots, 0.68 (0.05) at 20, 0.57
  (0.09) at 40, while the median coordinate is 0.01-0.09. The bias makes the
  allocator pessimistic (about 1.3x more rate than needed at 20 vs 40), not
  starved; later scheduled and selected attempts add to the same evidence.
  20 pilots kept as the plan's provisional size.

## What validity is claimed (stage 4)

Fixed kernels. Within a block between two plan rebuilds every quantity the
HMC step uses (likelihood precision H_s, step size eps_s, rate r_s, finite
difference scale h) is fixed. For subject s, given the group state, the
kernel is (1 - r_s) I + r_s K_s, where K_s is a mixture over the
state-independent jitter U of Metropolis-Hastings kernels on (alpha_s, v),
v ~ N(0, I): the proposal is L leapfrog steps followed by a momentum flip,
an involution with unit Jacobian (every leapfrog sub-step is a shear, also
with a finite-difference force, because the force depends on position
only). The abort rule is a function of the set of positions visited, which a
path and its reversal share, so it multiplies the acceptance probability by
a symmetric indicator. K_s is therefore reversible for the subject's full
conditional, which the metric S S' = (H_s + V^-1)^-1 may depend on because
the group state is fixed during the subject update. A shear round is a
multiple-proposal Metropolis step on its orbit and leaves the joint
posterior invariant. Each sweep is a composition of pi-invariant updates.
The warmup adapts step sizes after every transition, so it is not a fixed
kernel; it runs only in a subject's first burn block, whose draws are
discarded. The package tests check reversibility, invariance on a curved
target from exact draws (with a broken control that fails), boundary
behaviour, shear invariance from exact draws (with a broken-Jacobian
control), and serial/static/queued identity under a forced allocation.

Adaptation. Plans change only at the proposal rebuilds: every block in burn
and adapt (discarded stages), and in production at 100, 200, 400, ...
transitions (or an explicit schedule; `emc2.sample_rebuilds = numeric(0)`
freezes the first production plan). This is the increasingly rare
adaptation of AirMCMC (Chimisov, Latuszynski and Roberts, 2018): the
adapted parameter at each rebuild is a function of the history only, and
the gaps between rebuilds grow. In production the likelihood fits are not
refreshed, so the adapted parameter of these moves is (eps_s, r_s) in the
compact set [0.02, 2] x [0, 1] per subject (eps is clipped to that range)
and expected rounds x_k in [0, 4] per group coordinate (with the shear step
scales), pilot work per rebuild is bounded (a subject first eligible outside
burn gets at most 12 bracketing trials and 20 pilot transitions), and the particle kernel with its group-prior
component is part of every sweep.

What is not established. The AirMCMC convergence theorems also need a
drift condition holding simultaneously (or locally simultaneously) over the
family of kernels the adaptation can select. That condition is not proved
here for the subject HMC family, nor has it been proved for the existing
production proposal refreshes, which rely on the same schedule. Compactness
of the tuning set, bounded pilot work and the retained particle kernel are
restrictions that such a proof would use; they do not by themselves
establish containment. The adaptive production method is therefore
experimental and off by default. The exactly valid control is a frozen
production plan (`emc2.sample_rebuilds = numeric(0)`), under which
production runs one fixed pi-invariant kernel.
