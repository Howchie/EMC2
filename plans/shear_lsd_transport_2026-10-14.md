# Shear moves along slow group log sds: a transport that survives non-Gaussian subjects

## Issue

Some group log-sd coordinates mix slowly in the canonical sampler, and the
shear move barely helps them. The clearest case is log sd(B_Espeed) on lbaE.
Shear rounds on it move less than 10% of its window spread per round, even
though the allocator funds hundreds of rounds per block.

**Where the orbit breaks.** A shear move along a log sd changes the group
covariance and transports every subject with it, using each subject's fixed
Gaussian likelihood frame (Laplace fit, `lik_laplace`). Along that orbit, the
log density is flat near the current state and then falls sharply, mostly on
the side where the sd grows. The fall comes from:

- **A few tail subjects at a time.** Growing the sd stretches subjects in the
  tail of the group distribution further out.
- **Their own likelihood, not a bound.** The likelihood is skewed in log
  threshold and falls off far faster than its Gaussian frame predicts.
  - Example: subject 39 is flat to +0.5 window sd, then loses 6, 30 and 113
    nats at +1, +1.5 and +2.
  - The hierarchical term changes by 1–2 nats, the wall Jacobian is 0, and
    t0 is far from its 0.05 floor.
- **A different subject on each draw.** The subject holding the orbit back
  changes from draw to draw. Total orbit spread in one chain was 54 at one draw
  and 3171 at the next.
- **Many small errors that add up.** Even after the cliff subjects are
  handled, every subject's small departure from its Gaussian frame contributes,
  and the sum over n subjects blocks long jumps. The same limit appeared on
  N-back's slow mean directions, and it is why the collapsed move fails at
  n = 93.

The aim is a transport (or a move) whose orbit along such a coordinate stays
flat enough to accept jumps of about 1 window sd, at a cost per round a
production fit can pay.

## What has been tried

### State-indexed curved frames (removed)

Curved frames were added in 07eb3cfd..142fc828. They were removed in 3efd260e,
restored, and removed again in 3b9c101d. Read `git show 3efd260e` for the code.

**Design.** Per block, an atlas of conditional-mode fits:
- **Nodes:** fits of selected subjects at group states along the slow
  coordinate's path.
- **Regression:** each subject's (lin, prec) regressed on the curve coordinates
  (ridge, linear or quadratic).
- **Gates:** leave-one-out, positive definiteness, and an orbit-spread
  validation.

**Trigger.** Up to 4 coordinates could be curved. A coordinate needed at least
5 rounds and ineffective rounds (g >= 0.9) in the last block. Mean directions
also needed a pilot gain of at least 2.

**Outcome.**
- lbaE (6 seeds): no significant difference; slowest-direction ESS/CPU paired
  t was about -0.5.
- T30 (3 seeds): neutral.
- lba7r (2 seeds per arm): a wash; seed-to-seed variance dominated.
- Direct cost: small (7.5 s across two rebuilds on lba7r).

**Why it failed:**
- **Selection.**
  - Only lsd(B_Espeed) ever triggered on lbaE. Most log sds got 2–4 rounds and
    could never be measured.
  - The mean-direction pilot gain was 1.1, below the required 2.
  - Subjects were chosen at one rebuild state, and rebuilds double
    (100, 200, 400, ... transitions). The binding subject usually had no nodes.
- **Staleness.** The group state moves a lot per iteration. An atlas fitted at
  the rebuild state is wrong one draw later.
- **Shape.** A fit linear in lsd only partly follows a skewed likelihood:
  subject 39 went from -113 to -11 nats, not to 0.

### Replacement variants, tried offline

All on lbaE, replaying one rebuild. The table gives the total orbit spread (sum
over subjects of the variance of the orbit log density over the ±2 window-sd
grid) across 8 draws of one chain:

| Frames | Total orbit spread | Verdict |
|---|---|---|
| Fixed Laplace frames (current) | 54 – 150,811 | baseline |
| Removed curved frames | 15 – 10,360 | weak |
| Every subject refitted (BFGS + Hessian) at every state, started at the fixed frame's predicted mode (deterministic in the group state) | **15 – 47** | the ceiling for Gaussian frames |
| One Newton step on the mode from the fixed frame | blows up (up to 1e9) | out |
| Damped Newton, block Hessian, 3 or 6 iterations | still ~1e5 on bad draws | out: updating the curvature is what matters |
| All-subject atlas at the rebuild state, interpolated in lsd | 15 at that state, 343–1550 a draw later | out: stale |
| Screened refit: refit only subjects whose Gaussian frame mispredicts the likelihood change (tolerance 0.5 / 0.1 / 0.02 nats) | 27–538 / 32–124 / 24–89 | at 0.02 it refits 26–59 of 80 subjects per state |

**What even the ceiling buys.** Mean acceptance of a jump along lsd(B_Espeed):

| Jump (window sds) | -1 | -0.5 | +0.5 | +1 |
|---|---|---|---|---|
| Fixed frames | 0.03 | 0.20 | 0.39 | 0.13 |
| Every subject refitted | 0.10 | 0.37 | 0.45 | 0.16 |

The refit removes the single-subject cliffs, but the remaining non-Gaussian
residuals summed over 80 subjects still block jumps of ±1. It also costs about
6 core-s per state for 80 subjects, roughly 5 lbaE iterations of CPU per
candidate state. **Conclusion: Gaussian per-subject frames, however well they
are fitted, top out at about 1.2–1.8x better acceptance on short jumps here.
More work on Gaussian frames is not promising.**

### Related history (see the memory notes and git log)

- Bound walls for shear and collapsed transport (4be15aa8, f6ce7d54). These
  are active, and they help on lbaE.
- The collapsed surrogate move: allocator-gated. It fails at large n because
  surrogate errors sum.
- On N-back's slow mean directions, every transport tried left orbit sd at
  0.07–0.35 posterior sd: Laplace, regression, empirical linear response,
  translation and exact mode refits.

## Candidate directions (not exhaustive)

- **Non-Gaussian per-subject transport along the move.** The move pushes each
  subject along a known path, so transport it through its own conditional
  distribution instead of a Gaussian frame, using a one-dimensional
  quantile/CDF map along that path. The Jacobian comes from the map. This
  targets the skew directly.
  - Emulate it on the replays before building anything: does the summed
    residual over subjects shrink enough for ±1 window-sd jumps?
  - Keep the exactness condition in mind: frames and maps must be functions
    of the group state, not of the subjects being transported.
- **Moves that do not transport every subject.** Examples: partial
  transport, or updating the group sd jointly with only the subjects that
  carry the orbit, the rest held. This avoids paying the n-subject residual sum.
- **Fix the trigger and allocation coupling.** Coordinates that get few rounds
  are never measured. Measure every slow coordinate.
- **Cheaper exact refits.** Only if a direction shows the refit ceiling is
  worth paying for: C++ optimiser, warm starts that stay deterministic in the
  group state, analytic gradients.
- **Model-side relief.** Linking bounds to transforms (see
  `plans/bound_transform_linking_2026-10-14.md`) removes wall cliffs, but those
  are not the lbaE problem. Check whether any target fixture's cliffs are
  bound-driven before relying on it.

Other approaches are welcome. Judge them on the same numbers: orbit spread,
jump acceptance at ±0.5 and ±1 window sd, and cost per round.

## Materials

- **Replay state:** `/data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds`.
  This is an lbaE emc after one forced production rebuild with the curve-capable
  build. It still carries `shear$curve` and `shear$atlas`, which the scripts
  use to find the coordinate.
- **Checkpoint:** `/data/work/EMC2_workstreams/fixtures/lbaE_fit.rds` (an lbaE
  sample-stage checkpoint). Other fixtures are in the same directory
  (lba7r_sample_stage, T30_seed1_emc, forstmann).
- **Diagnostic scripts** in `plans/shear_lsd_transport_2026-10-14/`:
  - `replay_save.R` (force a rebuild and save the emc).
  - `terms.R` / `cliff.R` (per-subject orbit terms and their decomposition).
  - `rotate.R` (which subjects bind, by draw).
  - `oracle.R` / `variants.R` (frame variants).
  - `accept.R` (jump acceptance).
  - `curvelog.R` (logs curve-selection decisions during a continuation).

  They were written against a build that still had curves (4824b8dc): they read
  `ch$shear$curve$l` for the coordinate and call `.emc_shear_fits_at`. On
  current HEAD, take the coordinate index directly (`p + l`) and use
  `ctx$packed` for fixed frames. Restore kernels after `readRDS`
  (`restore_custom_kernel_pointers(restore_duplicates(...))`).

## Constraints

- One canonical sampler: no options framework. Switches are acceptable only as
  scaffolding while testing, and are removed before commit.
- Moves must not idle a chain's workers or unbalance chains: price in-process
  work at n_cores, use C++ kernels for hot loops, and keep subject likelihood
  work on the claim queue.
- Every mechanism must demonstrably accept on the case it is designed for.
  Something that is never built is not "rejected"; it means the gate or the
  construction is wrong.
- Screen with short replays and continuations first: about 400-iteration
  continuations saving draws, judged by lag-k against checkpoint variance. Only
  then run stop-rule benchmarks, one at a time under
  `/data/work/EMC2_workstreams/bench_lock.sh`, 3 chains x 9 cores, at least 3
  seeds. Expect large seed-to-seed variance: lba7r sample-stage time varied 2x
  between seeds with the same build.
- Test on the fixtures the mechanism targets (lbaE; lba7r; T30 and N-back for
  the summed-residual limit), not only the easy ones.
- Code comments say what and why only: no history, numbers, fixture names or
  dates.
- Only the C++ likelihood path is authoritative.

## Deliverables

1. An offline emulation, on the replays, of the chosen direction, reported as
   orbit spread, ±0.5 and ±1 jump acceptance, and cost per state. Compare it
   with the fixed-frame and refit-ceiling rows above. Stop and report if it
   does not clearly beat the ceiling.
2. If it does: an implementation in the shear move with tests, then screening
   continuations, then stop-rule benchmarks against the current HEAD.
