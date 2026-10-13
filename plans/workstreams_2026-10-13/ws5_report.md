# WS5 report: the ensemble update's cost, then its size

Branch `ws5/ensemble-cost` (worktree `/data/work/EMC2_ws5`), from `playground` at `0a5773dd`.
Private libraries under `/data/work/EMC2_ws5_libs/` (`base` = unchanged `0a5773dd`; `new*` = this branch).
Harness scripts in `/data/work/EMC2_ws5_libs/h/` (copied into `benchmarks/ensemble-cost/` on the branch).

Status (2026-10-13):

| Commit | Content |
|---|---|
| 7d875047 | Ensemble update in C++ (pool draws, stacked pools, Gibbs rounds), RNG-identical |
| ad3183b4 | Main group Gibbs step (standard, diagonal-gamma) in C++ |
| 873c734a | Iteration profile: `ensemble_elapsed`, `collapsed_elapsed`, `collapsed_moves` |
| c910a3c4, 141df019 | Between-block pool/rounds rule, two versions: **reverted** (07a5cf11, 37ba83dc) |
| 369ef06d | Proposal factors once per block; the inverse-Wishart draw's own precision factor; constant terms of the proposal density hoisted (RNG-identical) |
| 1bb98ee6 | Three ensemble rounds by default: **reverted** (1eeac033), ten stays |
| others | this report and `benchmarks/ensemble-cost/` |

Phase 1 is done and verified. Phase 2's rule is withdrawn (see "A between-block rule: tried and
withdrawn"): it valued the pool only by the group state's displacement inside one move and so
cut the pool to two members, losing a third of the worst coordinate's ESS per CPU-s on
Forstmann. The defaults stay M = 32, R = 10, now with measured support (see "Rounds" and "Pool
size: cost and sensitivity").

Tests on the final code (1eeac033, library `final2`, CPUs 23-31): full level test-ensemble-cpp,
-move, -pool-ll, -latent-output; default level test-subject-hmc, -shear-move,
-shear-standard-invariance, -collapsed-move, -move-allocation, -sampling-rejection,
-sample-proposal-schedule, -sampler-guards, -sample-store-inplace, -fit, -refit-alpha: no
failures, no errors (skips are the level-gated tests). nback31 was not run (its checkpoint build was stopped to release the shared
lock). Latent-factor hierarchies and standard ones with group designs keep the R rounds. The
`add_proposals()` hook the rule added is gone with it, so WS2's between-block code is untouched.

## Phase 1: what the move computes, and what was ported

### The arithmetic per iteration (step 2)

Inputs frozen for the block: per subject the pool proposal `q_s = (1 - w) N(mu_s, V_s) + w t_4(mu_s, c V_s)`,
`c = 1 + 3 / sqrt(p)`, `w = 0.25`. Per iteration, with `M` the pool size and `R` the rounds:

Before the subject step (`.emc_ensemble_draw_pools`, group RNG stream):
- per subject: `L_s = chol(V_s)` (lower), `M - 1` draws `x = mu_s + s_j L_s z_j`, `z_j ~ N(0, I_p)`,
  `s_j = 1` or `sqrt(c nu / chi2_nu)` with probability `w`; their densities `log q_s(x)` from the
  squared Mahalanobis distance `s_j^2 |z_j|^2` (no solve needed). RNG order per subject: the
  `p (M - 1)` normals, then `M - 1` uniforms, then one chi-square per tail draw.
- The `M - 1` likelihoods per subject ride on the particle round (workers).

The move (`.emc_ensemble_move`, group RNG stream), once per iteration:
- `k0_s ~ Uniform{1..M}` (position of the current state); stack `n M` rows; `log q_s` of the
  current state (one triangular solve per subject); `base_w = ll - log q` (group-independent).

Per round `r = 1..R`:
1. `alpha_k` = the selected row per subject (`p x n`).
2. Group Gibbs `h | alpha_k` (the package's conjugate update):
   - **standard**: `tmu ~ N((P0 + n V^-1)^-1 (P0 m0 + V^-1 sum_s alpha_s), (P0 + n V^-1)^-1)` via a
     Cholesky of the precision and `p` normals; residuals; per covariance block an inverse-Wishart
     `IW(v + d - 1 + n, 2 v diag(1/a_half) + S_block)` drawn as the inverse of a Wishart with scale
     `S^-1` through its Bartlett factor (`d` chi-squares with decreasing df, then `d (d - 1) / 2`
     normals filling the upper triangle column by column); per unblocked parameter one gamma;
     then `p` inverse-gamma `a_half` draws.
   - **diagonal-gamma**: `p` independent normals for `tmu` when the prior precision is diagonal
     (a joint draw otherwise), then `p` gammas for the precisions.
3. Subject-prior log density of every one of the `n M` rows under the new `h`:
   standard `-(1/2)(p log 2 pi + log|V| + |L^-1 (x - tmu)|^2)` with one Cholesky of `V` and one
   triangular solve of a `p x nM` matrix; diagonal-gamma the same with the diagonal.
4. Reselection: `k_s = argmax_j (base_w_j + log g_j + G_j)`, `G_j` iid Gumbel from `n M` uniforms in
   stack order (Gumbel-max is an exact categorical draw).
5. On the last round only, the weight ESS per subject from the same log weights.

After the rounds: the selected rows become the subject states and their likelihoods; the group
state is the last round's `h`.

**What is hoisted and what is not.** `base_w`, the stack and the current state's `log q` are
computed once (already so in R). Nothing else can be: every round changes `h`, so the subject
prior density of all `n M` rows and the covariance factorisation must be redone per round, and
the Gibbs draws depend on the selection. The redundant work in the R version was not arithmetic
but overhead: the scratch sample store (`.emc_ensemble_store`, `fill_samples`) rebuilt per round
only so that `gibbs_step()` could read the previous state; `riwish()`'s guard (`eigen` + `kappa`
+ two `solve`s per block per round); a `t()` and `matrix()` of the stack per round; and about
40 small R-level matrix calls per round. For the latent-factor hierarchies (factor,
infnt_factor, SEM) the group update has no cheap closed form of the same shape (it draws the
scores, loadings and precisions), so those keep the R rounds.

### Implementation (step 3)

`src/ensemble_move.cpp`, one translation unit:
- `emc_ensemble_draw_pools_cpp`: the pool draws, densities and the per-subject factors
  (`L`, `mu`, `logdet`, `c`) for every subject in one call.
- `emc_ensemble_stack_cpp`: the stacked pools, their log-likelihoods and `base_w`, with the
  current state's `log q` (type-independent; also feeds the R rounds).
- `emc_ensemble_rounds_cpp`: all `R` rounds for standard (no group designs) and diagonal-gamma.
- `emc_gibbs_standard_cpp`, `emc_gibbs_diag_gamma_cpp`: the same group update as the main
  iteration's Gibbs step (step 4); `gibbs_step_standard` routes to it when there are no group
  designs, `gibbs_step_diag_gamma` always. The R update with group designs remains as
  `.gibbs_step_standard_r`.

Random numbers: every draw goes through R's generator (`norm_rand`, `runif`, `rchisq`, `rgamma`)
in the order the R code consumed them, so `.Random.seed` after each call is bit-identical. The
one R-side draw kept in R is `sample.int(M, 1)` for `k0` (the `R_unif_index` path). What differs
is rounding only: LAPACK/Armadillo factorisations and summation order, and `-ffast-math`.

What the R rounds still do: the latent-factor hierarchies (`.emc_ensemble_rounds_r`, through
`gibbs_step()` and `.emc_ensemble_log_g()`), and standard hierarchies with group designs. These
were not profiled as significant and are not among the four fixtures.

Follow-up (369ef06d), after the external review: the pool proposals are frozen between
rebuilds, so `run_stage` computes their Cholesky factors and log determinants once
(`emc_ensemble_factors_cpp`) and passes them to every later pool draw of the block; the
inverse-Wishart draw returns its Bartlett factor `U = Z C` (upper triangular, the drawn
precision's Cholesky factor), so the precision block is `U'U` rather than a second inversion,
and the rounds' group density is evaluated block by block through those factors (no Cholesky
of the full covariance, no dense p x p triangular solve over all pool rows); the proposal
density's constant terms (two `lgamma` calls) are computed once per call instead of per pool
member. The Gumbel-max selection stays: replacing it would change the random numbers. Identity
(`sbs_ens2.R`, build `cpp2`, one chain, 9 cores, 200 iterations): Forstmann and T30 both 0
selection and 0 RNG mismatches; largest relative differences 1.2e-14 (Forstmann) and 1.4e-12
(T30) over pool draws, moves and 2200 Gibbs calls each.

Guard semantics: the R `riwish()` rejected a scale matrix by `eigen`/`kappa` thresholds
(`rcond < 1e-12`); the C++ path rejects when the symmetrised scale, its inverse or the Wishart
draw fails a Cholesky (`inv_sympd`). Both are error paths classified as numerical rejections by
`.emc_classify_failure()`; no accepted draw differs.

### Identity evidence (step 5)

Harness `sbs_ens.R`: a real sample-stage run on the new build where every pool draw, every
move and every main Gibbs step is also computed by the pure-R reference from the same inputs
and the same `.Random.seed`; the chain continues on the C++ result. Counted: selection
mismatches (the chosen pool rows differ), RNG-state mismatches after the call, and the largest
relative numeric difference. One chain, 4 cores, seed 1.

| Fixture | type | p | n | calls (draw / move) | selection mismatches | RNG mismatches | max rel diff tmu / tvar / a_half / wess |
|---|---|---|---|---|---|---|---|
| Forstmann | standard | 4 | 19 | 200 / 200 | 0 | 0 | 1.1e-14 / 1.9e-15 / 8.8e-15 / 2.2e-16 |
| T30 | standard | 13 | 40 | 600 / 600 | 0 | 0 | 1.1e-12 / 6.8e-14 / 1.0e-12 / 1.4e-15 |
| lba7r | diagonal-gamma | 20 | 106 | 1600 / 1600 | 0 | 0 | 1.3e-15 / 7.8e-16 / - / 2.2e-16 |

(`fit(stage = "sample", iter = n)` keeps trying the stop rule, hence 600 calls on T30 and 1600 on
lba7r for the requested 200 and 100.)

The main Gibbs step's port (second commit), same harness on build `head`, every `gibbs_step()`
call (the main one and the reference rounds') computed both ways:

| Fixture | Gibbs calls | RNG mismatches | max rel diff tmu / tvar / a_half |
|---|---|---|---|
| Forstmann | 2200 | 0 | 6.0e-15 / 2.4e-15 / 1.1e-14 |
| T30 | 2200 | 0 | 6.5e-13 / 4.7e-14 / 3.2e-12 |
| lba7r | 1100 | 0 | 1.1e-15 / 6.7e-16 / - |

Unit test `tests/testthat/test-ensemble-cpp.R`: pool draws, the stack's densities, the rounds
and the group updates against R references on random states for standard (unblocked, two
blocks, blocks mixed with gamma parameters) and diagonal-gamma (diagonal and correlated prior
precision): identical selections and `.Random.seed`, values within 1e-8; a control with one
extra round must change the RNG state.

### Timings (step 1 and step 5)

Paired per-part profiles (`prof_pair.sh`): base 0a5773dd against head (the two Phase 1
commits), one chain on 3 cores pinned to CPUs 23-31 with other work on the machine, (32,10)
throughout; head run twice. ms per iteration. Base's Gibbs column includes the ten Gibbs calls
inside the move's rounds; head's is the main step alone (the rounds' updates run inside the
move's C++). The stop rule ran a different number of iterations per run (T30 600 base, 500
head; lba7r 800 base, 1100 head), so whole-iteration wall times are indicative; the per-part
times are not affected.

| Fixture | move: base -> head | pool draws | Gibbs | iteration wall |
|---|---|---|---|---|
| Forstmann | 11.5 -> 1.7 / 0.8 | 2.3 -> 0.86 / 0.84 | 7.2 -> 0.21 / 0.09 | 47.5 -> 30.2 / 34.5 |
| T30 | 23.4 -> 1.9 / 2.6 | 5.6 -> 1.4 / 1.5 | 12.0 -> 0.16 / 0.19 | 164 -> 92 / 116 |
| lba7r | 26.0 -> 5.5 / 5.3 | 18.7 -> 7.8 / 7.1 | 2.1 -> 0.20 / 0.21 | 1487 -> 1306 / 1392 |

The ensemble's serial chain-process cost falls 5-10x. On Forstmann and T30 that is a large share
of the iteration; on lba7r the likelihood dominates and the wall difference is within load
noise.

## Phase 2: the pool size and the rounds

### What M and R buy (step 6)

Harness `sweep_mr.R`: on saved sample-stage states (8 states spread over chain 1's sample
draws, 4 replicates each), fresh pools from the chain's frozen proposals at pool size M, their
likelihoods, then the rounds chained one at a time up to 20 (identical to one call with R
rounds). After each round: the squared displacement of the group means and log sds from the
state the sweep started at, standardised by the posterior variance of the chain's draws and
averaged over coordinates (2 = an independent draw, on average over states); the share of
subjects whose selection differs from the start (`moved`); the fresh members' mean weight
relative to the carried state's (`eff`, the odds of displacing it); and the CPU of the
likelihood rows and of the rounds. Build `p2` (this branch). Single seed per state: indicative.

Note what the rounds are: round 1's group draw conditions on the subject states the sweep
received, so it is the same draw the main Gibbs step makes; the reselections' effect on the
group state appears from round 2 on.

**Forstmann** (standard, p = 4, 19 subjects, LNR): pools are efficient (`eff` 0.9 median
0.8) and the group decorrelates in one round.

| M | moved after 1 round | moved after 20 | jump mu, round 1 / 3 / 10 | jump log sd, round 1 / 3 / 10 | ll CPU per sweep |
|---|---|---|---|---|---|
| 2 | 0.41 | 0.41 | 1.6 / 1.7 / 2.5 | 1.7 / 1.5 / 2.3 | 3 ms |
| 4 | 0.67 | 0.70 | 2.3 / 1.8 / 1.8 | 1.4 / 1.7 / 1.6 | 3.5 ms |
| 8 | 0.82 | 0.84 | 2.4 / 1.9 / 2.0 | 1.9 / 1.7 / 1.7 | 4 ms |
| 32 | 0.97 | 0.97 | 1.7 / 2.0 / 1.8 | 1.2 / 2.0 / 1.8 | 9 ms |
| 128 | 0.98 | 0.99 | 2.2 / 2.1 / 2.3 | 1.4 / 2.1 / 2.2 | 27 ms |

The means are at the independence limit from the first (plain Gibbs) round at every M. The log
sds gain about 0.3-0.5 from the first reselection (round 1 to rounds 2-3) and nothing after;
the noise of 32 sweeps is about ±0.3, so no M dependence is resolvable here. The rounds cost
0.1 ms each (M = 32).

**T30** (standard, p = 13, 40 subjects, RDM): pools are inefficient (`eff` median 1e-4 to
1e-3; the mean 0.03 comes from a few subjects) and nothing the pool does moves the group.

| M | moved after 1 round | moved after 20 | jump mu, round 1 / 3 / 10 / 20 | jump log sd, round 1 / 3 / 10 / 20 | ll CPU per sweep |
|---|---|---|---|---|---|
| 2 | 0.006 | 0.007 | 0.67 / 0.59 / 0.41 / 0.59 | 0.35 / 0.33 / 0.39 / 0.40 | 8 ms |
| 8 | 0.027 | 0.048 | 0.49 / 0.75 / 0.49 / 0.69 | 0.34 / 0.39 / 0.46 / 0.42 | 13 ms |
| 32 | 0.087 | 0.20 | 0.69 / 0.60 / 0.64 / 0.69 | 0.31 / 0.34 / 0.42 / 0.34 | 29 ms |
| 128 | 0.20 | 0.43 | 0.60 / 0.59 / 0.62 / 0.64 | 0.40 / 0.38 / 0.43 / 0.46 | 95 ms |

The group jump is flat in M and in the rounds: it is the plain Gibbs conditional's width
(about 0.3-0.6 of the marginal) and the reselections add nothing to it, although at M = 128
and 20 rounds 43% of the subjects end with a fresh state. Those refreshed states lie within the
group conditional's reach. The rounds cost 0.2-0.5 ms each.

So, on these two fixtures, a pool larger than a handful of members and more than two or three
rounds buy no group movement inside one move. That was the wrong yardstick: the move's product
is the subject refresh (the `moved` columns), which grows with M on both fixtures and, where
the pools are inefficient (T30), with the rounds too (M = 32: 0.09 after one round, 0.20 after
twenty). Refreshed subject states change what the following iterations' group draws condition
on, which a single move's group displacement cannot show. The rule below was built on the
group displacement and failed for that reason.

### A between-block rule: tried and withdrawn

Two versions were implemented and reverted (c910a3c4, 141df019; reverts 07a5cf11, 37ba83dc).
Both set M and R between blocks to maximise the group means' and log sds' expected squared
jump per wall time, using the in-move group displacement by round, the fresh members' weight
odds e_s, and measured costs. Version 1 priced a pool member by the pool rows' share of the
particle round by row count and averaged displacement over coordinates; version 2 measured
the pool rows' worker CPU and maximised the slowest coordinate. Both share the flaw above: the
pool's benefit is the subject refresh share, which the objective never credits, so both drive
M to 2 whenever the in-move group gain is small (Forstmann: small because the means are already
independent; T30: small because the pools are inefficient).

Continuations from saved states with the first version (`cont_queue.sh`, build `p2`), 3 chains
x 3 cores pinned to 9 CPUs, one seed (indicative). CPU is the session cgroup's, which other
sessions share; wall times are clean.

| Fixture, 400 it | arm | worst coordinate | worst ESS | per CPU-s | lag-1 of worst | wall s |
|---|---|---|---|---|---|---|
| Forstmann | (32,10) | mu_m | 1086 | 11.7 | 0.03 | 35 |
| Forstmann | rule | lsd_t0 | 680 | 7.8 | 0.17 | 33 |
| Forstmann | (8,3) | lsd_t0 | 733 | 8.4 | 0.15 | 31 |
| Forstmann | (2,2) | lsd_t0 | 514 | 5.9 | 0.31 | 29 |
| T30 | (32,10) | lsd_B_Eneutral | 110 | 0.23 | 0.89 | 141 |
| T30 | rule | lsd_B_Eneutral | 165 | 0.34 | 0.90 | 132 |
| T30 | (8,3) | lsd_v_lMTRUE:Eaccuracy | 152 | 0.33 | 0.35 | 121 |
| T30 | (2,2) | lsd_B_Eaccuracy | 54 | 0.12 | 1.07 | 122 |

lba7r (200 iterations asked; the stop rule's tries gave 400 draws/chain): (32,10) worst
coordinate mu_v.g_lure, ESS 63, 0.011 per CPU-s, 0.08 per wall-s, lag-1 0.92, wall 791 s. The
first-version rule arm was stopped after 31 minutes (2.4x the fixed arm, no result; that rule is
superseded); (2,2) dropped as clearly worse.

On Forstmann the rule's cut to M = 2 slowed log sd(t0) (lag-1 0.03 -> 0.17) and cost a third
of the worst ESS per CPU-s for a 2 s saving in 35. On lba7r the rule arm was stopped after 31
minutes, 2.4x the fixed arm's time. An exact-length rerun (`diag_lba7r.R`, 300 iterations,
`max_tries = 1`) showed the per-iteration cost unchanged (2.2 s per iteration; a chain
iteration 0.9-1.2 s at M = 32 and 1.1-1.3 s at M = 2-6 after the cut): the long run was
`fit()` repeating 100-iteration tries (up to `max_tries = 20`) because the stop rule kept
failing with the smaller pool, i.e. a mixing loss, not a cost one.

### Rounds

Ten rounds stay the default. A three-round default (1bb98ee6) was committed on the strength
of the in-move group displacement and reverted (1eeac033) once continuations showed the loss.
Forstmann, M = 32, exactly 400 iterations, 3 chains x 3 cores, 4 seeds per arm (fresh chain
streams per seed; SE of the worst ESS about 20-30), build `cpp2`:

| R | worst ESS (log sd t0) | lag-1 of worst | alpha min ESS | CPU s |
|---|---|---|---|---|
| 3 | 865 | 0.14 | 875 | 92.7 |
| 5 | 1031 | 0.06 | 943 | 92.5 |
| 10 | 994 | 0.03 | 895 | 92.1 |
| 20 | 1077 | 0.03 | 922 | 92.0 |
| 40 | 1014 | 0.06 | 909 | 95.6 |

The benefit settles between three and five rounds; from five to forty the arms are within about
two SE, the lag-1 of the slowest coordinate bottoms out by ten, and forty starts to cost CPU.
The rounds are Gibbs on the group state and the selections given the pools, so as R grows the
move tends to an exact draw given the pools and the gain is bounded, while the serial cost
(about 0.1 ms per round on Forstmann at M = 32, 0.2-0.5 ms on T30) grows linearly. Ten keeps
a margin over the settling point at a cost of about 2% of a Forstmann iteration. Checkpoints
carry their recorded rounds in `production_options`.

### Pool size: cost and sensitivity

Benefit. The move's product is the share of subjects it refreshes, m(M) = mean_s (M-1) e_s /
(1 + (M-1) e_s), with e_s the fresh members' weight odds against the carried state, i.e. how
well the frozen proposal matches the subject's conditional posterior. It saturates once M - 1
exceeds about 1/e_s, so the useful M is set by the proposals' efficiency, not by the number of
subjects. Forstmann (e about 0.9), exactly 400 iterations, 3 seeds (2 for (32,10), (64,10)):

| M | refreshed per move | worst ESS, R = 3 | worst ESS, R = 10 | CPU s (R = 10) |
|---|---|---|---|---|
| 2 | 0.42 | 682 | 596 | 86.4 |
| 8 | 0.84 | 812 | 936 | 86.7 |
| 19 (= n) | 0.93 | 714 | 1003 | 88.8 |
| 32 | 0.96 | 859 | 957 | 93.0 |
| 64 | 0.98 | 866 | 1088 | 103.2 |

M = 2 is clearly worse; from M = 8 up the differences are within noise (replicate SD about
50). On T30 (e about 1e-3) even M = 128 refreshes 20-43% of the subjects (sweep above), and no
practical pool fixes proposals that poor.

Cost. Each iteration evaluates M - 1 fresh members per subject, one likelihood row each, in the
workers' particle round. Single-core timing of that call (`pool_cost.R`, rows from the
subject's own posterior draws), against the CPU of one chain iteration from the continuations:

| Fixture | n, p | per subject call: fixed + per row | pool CPU per iteration, M = 8 / 32 / n | share of iteration CPU, M = 32 / n |
|---|---|---|---|---|
| Forstmann | 19, 4 | 0.13 ms + 0.008 ms | 3.4 / 7.2 / 5.1 ms | 9% / 7% |
| T30 | 40, 13 | 0.30 ms + 0.017 ms | 16 / 32 / 37 ms | 8% / 10% |
| lba7r | 106, 20 | 0.85 ms + 0.19 ms | 0.21 / 0.69 / 2.16 s | 15% / 46% |

The cost is n x (fixed call + (M - 2) x row): linear in M per subject, and the fixed per-call
part dominates below about 15-20 rows on the cheap models (batching the pool rows into the
particle step's likelihood call would remove it). Wall time hid this cost at 3 cores per chain
(Forstmann wall flat at 20-22 s from M = 2 to 64) because the workers have slack; the CPU is
real. M = n is not a good rule: n does not enter what a subject gains, while the cost becomes
quadratic in n (lba7r: M = 106 would add about 1.5 s of CPU per chain iteration over M = 32).
M = 32 is past saturation on efficient proposals at a cost of 8-15% of the iteration on these
fixtures. A pool sized by the measured e_s (M - 1 a few times 1/median e_s, capped) is the
principled refinement, left for later; on poor proposals (T30) it would hit the cap and buy
little, so the cap is what would matter.
