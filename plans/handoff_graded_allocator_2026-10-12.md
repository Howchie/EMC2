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
(`$S/bench/C31`): def seeds 8.5 to 9.6 per kCPU-s; `fixed` (mech4) 11.96,
no HMC bought, no regression; `graded` (mech7) was running at handoff time
(`$S/bench/C31/graded_1.log`). Phase A references: `$S/bench/P_mech`,
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

## Idle time: priced, not fixed, not yet attributed

CPU per wall-second (27 cores): def 8.4, refresh 7.4, graded bounded 11.0.
The shear frame build (80 subjects, in the chain process) and the collapsed
context are centralised work that idles the workers; they are now priced
as fixed costs so the allocator does not buy rounds they do not pay for,
but when it does the idling remains. Attribute it with a profiled
continuation: add `PROFILE=1` to continue_arm.R (set
`options(emc2.sampler_profile = TRUE)` and save
`EMC2:::.emc_iteration_profile()`; the rows now carry gibbs, particle, fill,
response_wait, shear and hmc per iteration) and run the graded build on the
lbaT checkpoint. If the frame build is a large share, move it to the
workers (same shape as the collapsed move's C++ inner chain).

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
