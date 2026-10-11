# EMC2 shear move review

## 1. Overhead

I profiled 40 iterations from the saved freeA15 warmup state, with 30 subjects, 14 group coordinates, and 15 workers. CPU time includes the parent and worker processes: parent-only timing misses much of the work. The package was loaded from an isolated temporary install; no project files were changed during the investigation.

| Run | Total CPU | CPU per iteration |
|---|---:|---:|
| Ensemble | 24.8 s | 0.620 s |
| Shear | 79.5 s | 1.987 s |

The paired difference was **1.367 CPU seconds per iteration**. About 0.19 s/iteration of that difference was extra ordinary particle-likelihood work as the two chains diverged. Within the shear step, the largest measured costs were about **0.94 worker CPU seconds** for exact likelihoods and **0.24 parent CPU seconds** in the shear move, including about 0.16 seconds dispatching worker calls.

The chain had only 27 usable subject fits, so the all-or-nothing delayed-acceptance gate disabled screening. Shear then made about **28 exact likelihood rounds per iteration**, evaluating all 30 subjects each time. Planning was negligible: a median **0.001 CPU seconds** per rebuild. State construction made about 29 candidate states per iteration, or roughly 870 subject frames; frame building used about **0.038 CPU seconds per iteration**, with Cholesky itself about 0.008 seconds. Reusing frames for mean-only proposals could trim some of this cost, but the exact likelihood rounds dominate.

A same-checkpoint, 50-iteration runtime comparison showed that allowing delayed acceptance with partial fits cut exact evaluations from 3,642 to 2,563 across three chains: **30% fewer**. The existing likelihood-delta code skips subjects without fits, and their exact likelihood terms remain in the final correction. A toy exact-posterior test of this partial-fit screen gave the expected target summaries. This is the clearest target-preserving overhead reduction; it leaves easy runs unchanged because shear is optional.

## 2. Seed 2 inconsistency

The saved production summary shows poor mixing: for example, group A had **Rhat 1.180, ESS_bulk 12.2**; subject A|21 had **Rhat 1.209, ESS_bulk 10.7**. The saved ensemble run had group A at **Rhat 1.024, ESS_bulk 56.4**.

I replayed 100 transitions from the saved seed-2 warmup checkpoint with `SAVE_EMC=1`. Group A’s three chain means moved from **−1.196, −1.808, −1.244** to **−1.146, −1.386, −1.130**. The second chain moved toward the others, while the replay still showed group-A **Rhat 1.16, ESS_bulk 13.4**. That supports slow between-chain movement along A; it does not establish a second mode or a permanently stuck chain.

Warmup plans selected all mean and scale directions: `lambda` stayed around 0.996–0.9999 and all scale weights were 1. Mean-direction acceptance was roughly 0.39–0.45. I found no evidence that directions were omitted, though these checks do not show whether the learned directions remain well aligned during production. Two replay chains had incomplete subject fits, which disabled delayed acceptance under the current gate. The measurements do not isolate a single root cause, and I did not find evidence sufficient to blame a wrong Laplace direction or the step-size adaptation.

## 3. Fix for the ridge

**Curved transport is the strongest next candidate for improving mixing, but the freeA15 evidence does not yet establish it as a fix.** Across three slow directions, measured path widths changed from 0.243 to 0.220, 0.206 to 0.285, and 0.301 to 0.370 for linear versus quadratic transport. The median chainwise improvement was about **1.23×**, with one chain getting worse. Residual widths were still 0.22–0.37 against a unit-width Gaussian step.

I also tested the proposal on an exact toy posterior with three subjects whose conditional means curve quadratically with the group parameter. I fitted the transport on direct exact posterior draws, froze it, and used a symmetric group random walk with the exact MH ratio. The map shifts subject parameters by the fitted conditional-mean change; the full transformation has unit Jacobian and an exact reverse move. Across four chains, curved transport raised acceptance from **0.424 to 0.744** and group `ess_bulk` from **607 to 4,643**; `rhat` was 1.00 and posterior means and standard deviations agreed with the exact reference. This verifies the construction and the toy case, not its performance on freeA15.

The README’s frozen-group comparison reports much better subject ESS when the group is held fixed, which supports cross-level coupling as the problem. A subject-only move would leave the group coordinate unchanged, so I would test a joint group-plus-subject transport first. Keep it opt-in until it shows consistent gains on held-out freeA15 chains.

## Ranked recommendations

1. **Enable partial-fit delayed acceptance for shear.** It has a target-preserving correction and cut exact calls by 30% in the short comparison. Confirm with a longer matched-checkpoint run; compare ESS and Rhat at matched CPU, and revert if savings or mixing gains disappear.
2. **Continue curved transport as an opt-in mixing experiment.** On held-out seed-2 checkpoints, require consistent path-width improvement across chains, then compare fixed-CPU ESS and Rhat. The current mixed path results are not enough to promote it.
3. **Reuse frames for mean-only proposals.** Check equality of proposed states and MH ratios over random moves, then repeat the short CPU profile. Drop the optimization if it adds complexity without a measurable speedup.

I would keep the ensemble move as the default. The current evidence supports optimizing optional shear and testing curved transport further, without imposing extra cost on easy problems.

Summaries used `posterior::ess_bulk` and `posterior::rhat`. Relevant sources: [investigation README](../benchmarks/group-move/results_shear_investigation/README.md), [shear move](../R/shear_move.R), [ensemble move](../R/ensemble_move.R), [new_particle()](../R/sampling.R), and [curved-path feasibility script](../benchmarks/group-move/results_shear_investigation/scripts/curve_path.R).
