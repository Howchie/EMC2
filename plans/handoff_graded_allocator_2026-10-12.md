# Handoff: graded move allocator, lbaT starvation, refresh (2026-10-12)

State of the `playground` branch at HEAD 65268ee3 (plus this document).
Scratch directory of session d88c8a86 (referred to as `$S` below):
`/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/d88c8a86-efae-415e-999b-147668d64e27/scratchpad`.

## What was found

1. **lbaT does not converge on any build**, and run-to-run variation on the
   same build is as large as the between-build difference (seven phase A
   fits, seed 1: R-hat 1.06 to 1.45; all hit the 20-try cap). The slow
   direction is the group log sd of A with the subjects' A: chains sit in
   the funnel neck (log sd(A) about -5.4) or the bulk (about -1.8) and
   cross between them over hundreds of iterations. One posterior, exact
   kernels, slow mixing.
2. **The cause was allocation starvation, not a missing move.** The shear
   rounds on log sd(A) accept and mix it wherever they run (chain 3 of the
   checkpoint: 1.9 rounds per iteration, 71% moved, lag-100 autocorrelation
   0.04). Chains 1 and 2 got 0.04 and 0.01 rounds per iteration. Five
   rounds taken in the wrong place read as a factor of 0.995, and a round
   was priced at 1.59 against a rest-of-iteration cost of 0.56 because the
   per-iteration frame build was divided by the few rounds. Nothing cleared
   the gain threshold, nothing was remeasured. The scaffold arm `force4` (4
   rounds per iteration on log sd(A), no rebuild) lifted its ESS from 6 to
   92 and proved the lever. The collapsed move is bought at zero on LBA
   because its acceptance is ~1e-12 (surrogate errors coherent across 80
   subjects); that is correct allocation and not tonight's target
   (plans/collapsed_move_lba_acceptance_2026-10-12.md, superseded section,
   and the two reviews under plans/collapsed_move_lba_review_2026-10-12/).
3. **Wall-clock "doubling" on lbaT was a load confound**: paired fits on an
   idle machine put the current build at 1282 s against 1469 s for the
   revert build, with equal CPU per iteration. N-back 31 has no revert-build
   reference fit; one was queued and killed, so that comparison is still
   open (`$S/bench/A_revert` empty).

## What was changed (all committed, tests green)

| commit | change |
|---|---|
| 58620bf4 | profiling: per-iteration rows survive multi-block fits (`.emc_iteration_record/_profile`) |
| 55b29f4d, f57c955d | allocator: a coordinate no move can lower does not bound the target; the collapsed move and a regime count as movers |
| c486378a | one total HMC rate under a regime (`r_min`), exact zero movement (no `1 - 1e-6` clamp), pilot gate on the total rate, zero-rate regime is identity, re-solve at shared particle rates, curve switch gates existing atlases |
| f6ce7d54 | wall transport: a subject outside its wall takes the mirror quantile map (was not a bijection; CODEX review) |
| 54dd5ccf, 65268ee3 | **graded allocator**: one cost table per move family (`w$fam`: cpu, own, iters, moves; `.emc_alloc_marginal`), marginal cost per move + fixed cost per iteration-with-moves paid once (with-versus-without, same movable set), evidence-weighted factors (`.EMC_ALLOC_ROUND_EVIDENCE = 20` moves equal the reference: shear vs its kind's pooled factor, collapsed vs the prior, subject HMC vs the subjects' pooled factor), log-sd rounds pooled across chains, evidence decay `.EMC_ALLOC_EVIDENCE_DECAY = 0.5` per production rebuild (shear pool, collapsed pool, HMC evidence), exploration floor `(20 - n)/B` and cap `max(20 - n, 2n)/B` for shear coordinates and the collapsed move |

Key code: R/move_allocation.R (header comment is current), R/shear_move.R
(`round_el`, pool decay), R/collapsed_move.R (`round_el`, pool decay),
R/subject_hmc.R (`.emc_hmc_decay_evidence`), R/sampling.R (cost call with
`moves = list(shear, collapsed, hmc)`). Tests: test-move-allocation.R (842),
test-collapsed-move.R (90), test-shear-move.R (323), test-subject-hmc.R (65),
test-profile-schema.R (487). The full suite was not run tonight.

## Results (lbaT seed 1, 1800-iteration continuations of one checkpoint)

Checkpoint `$S/bench/P_mech/lbaT_fullc_1_emc.rds` (current build phase A,
profiled); results `$S/bench/CL/<arm>_1.{rds,log}`; summarise with
`benchmarks/subject-hmc/summarise_cont.R $S/bench/CL def`. Libraries:
Rlib_bench_mech (old), mech2 (old + refresh), mech3 (floor/regime fixes),
mech4 (+ wall fix), mech6 (graded, unbounded), mech7 (graded, bounded = HEAD).

| arm | lib | wall s | CPU s | log sd(A) ESS | mean A ESS | slowest-dir ESS | per kCPU-s | per 1000 wall-s |
|---|---|---|---|---|---|---|---|---|
| def | mech | 998 | 8387 | 6 | 5 | 5 | 0.61 | 5.0 |
| static (atlas kept; = def) | mech | 1019 | 8462 | 6 | 5 | 5 | 0.60 | 4.9 |
| noreg (regimes off) | mech | 777 | 7367 | 6 | 5 | 5 | 0.66 | 6.4 |
| refresh | mech2 | 1032 | 7671 | 44 | 24 | 16 | 2.06 | 15.5 |
| fixed (floor fixes) | mech3 | 1017 | 8283 | 7 | 5 | 6 | 0.67 | 5.9 |
| force4 (scaffold) | mech3 | 1126 | 9413 | 92 | 8 | 7 | 0.78 | 6.2 |
| graded, unbounded | mech6 | 407 | 5140 | 8 | 6 | 5 | 0.94 | 12.3 |
| **graded, bounded (HEAD)** | mech7 | 500 | 5495 | 60 | 7 | 7 | 1.36 | 14.0 |
| graded + refresh | mech7 | 536 | 4918 | 15 | 4 | 4 | 0.86 | 7.5 |

Single seed each; lbaT continuation ESS is small and the halves disagree,
so these rank mechanisms only where they differ several-fold. N-back 31
(`$S/bench/C31`, 1800 iterations):

| arm | lib | wall s | CPU s | slowest-dir ESS | per kCPU-s | per 1000 wall-s |
|---|---|---|---|---|---|---|
| def (3 seeds) | mech | 923-945 | 16.8-17.2k | 143-161 | 8.5-9.6 | 15.2-17.2 |
| fixed (floor/regime/wall fixes) | mech4 | 961 | 17.4k | 208 | 12.0 | 21.6 |
| graded, bounded (HEAD) | mech7 | 1119 | 19.4k | 207 | 10.6 | 18.5 |

So on the easy fixture the graded build mixes like `fixed` but costs 15%
more CPU and 20% more wall than `fixed` (still better than def on both).
That extra is the exploration: every under-measured coordinate runs its
floor rounds, the collapsed move runs ~0.07 per iteration (108-134 moves
per chain) and 1000 HMC transitions per chain were tried; the profiled run
below is the place to see whether the frames built for those rounds are
the cost. Against the user's rule (a large gain on hard models is worth a
small cost on easy ones) this is 2x on lbaT for ~20% wall on N-back; the
exploration rate (`.EMC_ALLOC_ROUND_EVIDENCE / B`) is the knob if that is
too much, and a cheaper option is to exempt coordinates whose kind already
has ample pooled evidence from the floor. Phase A references: `$S/bench/P_mech`,
`$S/bench/P_revert`, `$S/lbaT*`, `$S/bench/A/nback31_fullc_1.log`.

Reading: the graded, bounded build buys log sd(A) rounds in every chain on
its own (0.16, 0.40, 0.68 per iteration), runs the collapsed move only at
its exploration trickle, and halves the wall of def at 35% less CPU; mean A
is now the slowest direction. Regimes bought full-rate HMC for all of chain
1's subjects for nothing on the old build; with the regime accounting fixed
the graded build chose regime 0.

## Refresh: wanted, not safe yet (build this)

`emc2.alloc_subject_floor = "refresh"` lets the cheaper of the particle step
and HMC refresh a subject (q_s + r_s >= 1) instead of keeping the particle
step as a floor. On the old allocator it was the best arm (3.4x def) because
the double-counted regime ran both kernels: HMC at rate 1 on top of particle
rates 0.3 to 0.5, i.e. a lot of subject movement. On the graded allocator
it set q = 0 for every subject in every chain and ran one HMC transition per
subject per iteration (144k per chain); the subjects met their own lag-one
target more cheaply, but the group moves, Gibbs on V and the shear rounds
live off subject movement, and mean A collapsed (ESS 4, R-hat 1.87).

What it needs: the particle step's factor on the subjects' deviations from
the group mean, `p_dev`, measured like `h_dev` (the window already holds the
particle jumps `w$P`; score them against `d_var` instead of `s_var`), and the
group-credit model (`.emc_alloc_regime`, `r_bar`, `g_floor`) extended so that
a plan's subject kernels (q_s, r_s) set the deviation autocorrelation
`d_rho_s' = d_rho_s (1 - q (1 - p_dev)) (1 - r (1 - h_dev))` relative to the
measured one, and the group coordinates' rho follow it as the regime
coupling does now. Then standing the particle step down is charged as
higher group autocorrelation and refresh can be judged honestly. Two
further points from the review: with independent Bernoulli selection
`q + r >= 1` guarantees a refresh only in expectation (couple the selection:
run HMC when the particle step did not), and the arm must be re-run on the
same checkpoint against graded-bounded.

## Idle time: attributed (two profiles, 2026-10-12 02:45)

Profiled graded-bounded continuation on the lbaT checkpoint (`$S/bench/CL/
gradedP_1.rds`, PROFILE=1, 1800 iterations, 505 s wall) and the review's
300-iteration probe (plans/graded_refresh_proposal_2026-10-12/) agree:

- Per-iteration wall by chain 0.174 / 0.199 / 0.264 s; the slowest chain's
  iterations sum to 476 of the 505 s, so chains 1 and 2 spend ~30% of the
  wall waiting at the block barrier for chain 3, which bought 2.5x the
  shear rounds (1.48 vs 0.58 per iteration) and more HMC. **Chain imbalance
  from per-chain plans is the largest idle term.** Work outside iterations
  (rebuilds, proposals, progress checks) is 29 s of 505.
- Within an iteration: Gibbs 0.008 s, particle step 0.083-0.096 s (its
  response wait equals the workers' elapsed, i.e. the workers are busy),
  shear CPU 0.08-0.24 s, HMC CPU 0.02-0.05 s. Worker busy share 0.44 / 0.44
  / 0.30 of (workers x wall). Shear setup (frames) is only 2-3% of the
  shear's reserved cost (review), so moving the frame build to the workers
  is NOT the first target. 78-86% of HMC's reserved cost is fixed: a round
  attempting 4-9 subjects per iteration activates every iteration with its
  latency and idle reservation.

So the order is: (1) balance chains, either by pricing the barrier (a
chain's plan cost is really the block's slowest chain) or by sharing shear
round rates across chains as HMC/particle rates are; (2) sparse-family
activation pricing (charge the fixed cost at the probability any round or
transition runs, `1 - prod(1 - x_k)` for independent rounds; group low-rate
rounds on a shared activation); (3) HMC round latency/right-sizing;
(4) only then worker-resident frame/transport work.

**Benchmark limitation (review, confirmed):** the checkpoint holds 2000
production transitions and no family cost fields, the harness forces one
rebuild, and the doubling schedule puts the next at 4000, beyond the 1800
iterations. Every continuation above therefore allocated once from the
checkpoint's legacy evidence with the family fixed costs at zero and never
repriced. `emc2.sample_rebuilds` accepts an explicit schedule: the
continuation should run a calibration segment (rebuilds at +100, +300)
before the scored window, with DISCARD covering it. Not yet done; it is the
first thing to change before trusting any pricing conclusion, including the
N-back cost figure.

Also from the review and fixed (bb7c3acb): the HMC deviation factors were
shrunk toward a reference built from ordinary subject factors; they now
pool against their own. Noted, not changed: `.emc_alloc_record_subjects`
measures the subjects just before the HMC round against the previous
iteration's final state, so the deviation statistic spans the previous HMC
round and this iteration's transports; any refresh model must anchor to
that measured schedule (the review's `d_new = d_old * F_new / F_old`), and
coupled selection (one uniform per subject: particle if U < q, HMC if
U >= 1 - r) is the scheduler that guarantees a refresh.

## Clean-up list (not done)

- Remove `emc2.shear_interleave` (no-op without HMC rates; its HMC tasks are
  not recorded, review P2) and `emc2.alloc_hmc_regimes` if the N-back and
  one more lbaT seed confirm regime 0 is chosen; remove `emc2.shear_walls`,
  `emc2.shear_curve` switches once their arms are settled (static arm now
  strips atlases via `arm_hook`; a re-run with atlases stripped was queued
  and killed, result absent).
- The collapsed move stays, allocator-gated.
- `interleave`, `nowall`, `p60` arms have no usable result; `static` with
  atlases stripped has none.
- Memory files for this work: project_emc2_graded_allocator.md,
  project_emc2_alloc_floor_veto.md, project_emc2_lbaT_lsdA_starvation.md,
  project_emc2_wall_bijection_fix.md, feedback_waiters_and_claims.md,
  feedback_no_walltime_regression.md.
