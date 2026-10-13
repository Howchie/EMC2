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
| c910a3c4 | Between-block pool/rounds rule, first version (hook in `add_proposals()`: WS2) |
| 141df019 | Rule revised: measured pool-row CPU as the member price, slowest coordinate as the objective |
| others | this report and `benchmarks/ensemble-cost/` |

Phase 1 is done and verified. The first-version rule lost efficiency on Forstmann (see Check);
the revised rule is unit-tested only. Pending, queued in `h/queue_g.sh`
(= `benchmarks/ensemble-cost/ablation_forstmann.sh`): a Forstmann ablation of the two revisions
(3 seeds), then the test suites on build `final` (141df019). Whether the rule stays the
default, or (32,10) does, waits on that ablation; the rule is only to be run on T30 and lba7r
if Forstmann shows a gain. nback31 was not run (its checkpoint build was stopped to release the
shared lock). Latent-factor hierarchies and standard ones with group designs keep the R rounds.

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
rounds buy no group movement. The ensemble's value on Forstmann comes through the subject
refresh: the refreshed states let the next iteration's group draw be independent, which the
in-sweep displacement shows as the first-round jump already sitting at 2.

### The rule (step 7), first version

`R/ensemble_plan.R`, applied once per block where the pool proposal is rebuilt
(`add_proposals()`, one added call, which WS2 should know about). Objective: the group state's
expected squared displacement per unit of wall time; every quantity is measured in the block:
- inside the move (C++ returns the trajectory): the standardised displacement after each round
  from the state the move received, J_r; the move's gain G = J_R - J_1 (the reselections'
  displacement; J_1 is the Gibbs draw the main step also makes) and its increments dJ_r;
- the fresh members' weight odds e_s per subject, from the final weights;
- the chain's own lag-1 displacement of the same coordinates over the block, J_tot, and the
  iteration's wall time C, which set the exchange rate J_tot / C;
- the wall time of a round, c_round, and of a fresh member, c_member (the pool rows' share of
  the particle round's wall time, a pool row costing what a particle row costs).

Rounds: keep round r while dJ_r / c_round >= J_tot / C (its marginal displacement per cost
beats the chain's average); double when the last round still qualifies; floor 2 (one
reselection must be observable), cap 64. Pool: the reselections' displacement is taken to scale
with the share of subjects they refresh, m(x) = mean_s x e_s / (1 + x e_s) for x = M - 1 fresh
members (an independence sampler's displacement odds), G(x) = D m(x) with D = G / m(x_cur);
keep members while D m'(x) / c_member >= J_tot / C; double when the current pool still
qualifies; floor 1 fresh member (M = 2), cap 256. Both change only between blocks. The rule
has no dataset constants; the floor and caps are structural.

Unit tests (`test-ensemble-plan.R`): constructed block records with known answers (rounds
grow / stop after the last qualifying increment / floor; pool grows / stops at the last
qualifying member / floor when the reselections bring nothing / unchanged on a record too
short to price). Invariance (`test-ensemble-move.R`, full level): the joint posterior is
invariant when the pool and rounds change between blocks on a schedule
(2,2) -> (8,3) -> (24,6) -> (4,2).

Functional run (short Forstmann fit from scratch, 4 cores, build `dev`): at the first sample
rebuild the rule cut M from 32 to 2 and kept 9 rounds, then settled at M = 2, R = 4-5, with
the chain's lag-1 displacement unchanged (J_tot about 2 per coordinate before and after). The
measured prices: 31 fresh members were about 65% of the iteration's wall time on 4 cores and
bought G = 0.44 of displacement inside the move; a round about 5e-5 s.

### Check (step 8)

Continuations from saved states (`cont_queue.sh`, build `p2` = the first-version rule), 3
chains x 3 cores pinned to 9 CPUs, rebuilds every 100 transitions, one seed: indicative, and
only comparisons between arms are meaningful. Worst coordinate's ESS (bulk, over the group
means and log sds) and per CPU-s. The CPU figure is the session cgroup's, which other
sessions share, so it may include their work; wall times are clean. The harness now also
reports the run's own process-tree CPU (used for the ablation).

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

The first-version rule cut M to 2 on every chain of both fixtures at the first rebuild. On
Forstmann that cost a third of the worst coordinate's ESS per CPU-s: log sd(t0) slows (lag-1
0.03 -> 0.17) while the saving is 2 s of 35. Two faults:
- the price of a pool member was the pool rows' share of the particle round by row count. The
  likelihood is vectorised over rows with a per-call cost, and the round's wall time includes
  latency and idle cores, so that share overstated the members' cost several-fold (the
  functional run's "65% of the iteration" against a measured 2 s of 35);
- the objective averaged displacement over coordinates. The pool's gain sits on the log sds,
  and the slowest coordinate is a log sd; the average hid its loss behind the means, which sit
  at the independence limit from round 1.
On T30 the rule kept 6-10 rounds although the rounds' measured gain is noise around zero:
rounds are cheap enough that any positive noise passes. Harmless in cost.

### The rule, revised

Both faults fixed in `R/ensemble_plan.R`:
- the pool rows' cost is measured: each worker times the candidate likelihood of its
  subjects' pool rows (CPU), returned with the particle results (`extra_times` from
  `.emc_wpool_iter`, `extra_time` per subject on the fork path); the block's member price is
  that CPU spread over the chain's cores;
- the objective is the slowest coordinate's displacement per wall time, the next plan
  maximising, over pools up to twice the current and rounds up to the current,
  min_k (T_k + G_k(x', R') - G_k(x, R)) / (C + c_member (x' - x) + c_round (R' - R)),
  with T_k the chain's lag-1 displacement of coordinate k over the block, G_k(R') its
  reselection gain after R' rounds (negative counted as none) scaled by m(x')/m(x), the
  cheapest plan winning ties; the rounds double when the last round still raises the
  objective. Without a measured member price the pool is kept.
Tests (`test-ensemble-plan.R`): rounds grow / drop a round that adds nothing / floor with no
gain; the pool falls monotonically with the member's price and doubles when members are cheap;
a slow coordinate that the reselections move keeps the pool, while gains confined to already
fast coordinates cut it to one fresh member (the case the averaged objective got wrong); pool
CPU is spread over the chain's cores.

Check of the revision (queued, `ablation_forstmann.sh`): Forstmann continuations, 400
iterations, 3 seeds, one fit at a time on CPUs 23-31: fixed (32,10) and the first rule (build
`p2`), the price change alone (`abl_price`: 141df019 with the first-version objective), the
objective change alone (`abl_obj`: 141df019 with the row-share price), and both (`final` =
141df019). Results in `/data/work/EMC2_ws5_libs/abl/`. Not yet run.
