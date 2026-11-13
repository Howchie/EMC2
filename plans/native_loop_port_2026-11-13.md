# Native sampling-loop port: status and remaining work

Goal: move the structure of the sampling loop to C++ while leaving every
sampling mechanism (particle step, ensemble, shear, collapsed, subject HMC,
Gibbs, allocation) and every draw unchanged. Each step is checked against the
R path on the same random stream.

Fixture: 20-subject REXG, three chains, one core per chain, default stop rules
with 1,000 sample iterations. DMC2 on the same data: 7.8 s.

## Done

| Piece | Where | Check |
|---|---|---|
| Particle step, all subjects in one call (plan, Laplace conditional, efficient-proposal shift, adaptation) | `src/native_particle.cpp` `emc2_native_particle_subjects`; `.emc_native_particle_batch` | R path draws, all stages |
| Subject HMC transition, all subjects in one call | `emc2_native_hmc_subjects`; `.emc_native_hmc` | task lists equal to `.emc_hmc_run`, accepted / rejected / aborted |
| Prepared likelihood arguments per (data, model) pair | `src/ll_fast_cache.cpp`; `calc_ll_manager` | values equal to the full path |
| Batched subject likelihoods for shear, collapsed and ensemble candidates | `emc2_ll_cached_subjects`; `.emc_ll_subjects_fast` | R path draws |

A subject the native code cannot serve (no stored arguments, non-finite
operand, likelihood error) runs on the R path from the stream its request
carries, so failures are reported exactly as before.

Timings (seconds):

| Stage | Before | After |
|---|---|---|
| preburn (50) | 1.35 | 1.06 |
| burn (300) | 12.8 | 9.4 |
| adapt (100) | 4.6 | 3.25 |
| sample (1,000) | 38.8 | 28.5 |
| total | 57.6 | 42.2 |

One iteration in-process: 31 ms to about 18 ms.

Second pass (after commit d77173d9):

| Piece | Where | Check |
|---|---|---|
| Parameter mapper kept with its context and reused across calls | `PtMapperContext::spare`, `pt_mapper_reload` in `src/particle_ll.cpp` | R path draws |
| Shared-state objects built only on the path that reads them | `calc_ll_oo` race branch, `unused_state<T>()` | R path draws |
| Shear and collapsed likelihood rows built in C++ | `emc_candidate_rows_cpp` | R path draws |

A one-row likelihood call: 31 us to 14.8 us (about 8 us inside `calc_ll_oo`,
2.3 us of it kernel). One iteration in-process: about 13.4 ms. Whole fit:
42.2 s to 37.2 s (preburn 1.1, burn 8.8, adapt 3.0, sample 24.3).

The shear and collapsed arithmetic was already in C++
(`emc_shear_candidates_cpp`, `emc_shear_frames_cpp`,
`collapsed_trajectories_cpp`); after the row builder their R side is a flat
tail with no line above 0.15 ms per iteration.

## Remaining, by size

Items 3 and 4 below are done to the extent the profile supports (second pass
above); the rebuilds are now the largest structural cost: in a 400-iteration
sample window, 3.9 of 11.2 s. Per rebuild: forked jobs about 0.58 s
(`.emc_rebuild_lapply`, four calls: shear orbit steps 0.27, efficient
proposals, Laplace), `.emc_alloc_solve` 36 calls of 6 ms in the master
(0.22 s), `.emc_empirical_covariance_ok` about 410 calls (0.14 s). Chain
blocks themselves have no measurable fork overhead (master wall equals the
slowest chain, 15-17 ms per iteration with three chains running).

1. **Proposal rebuilds between blocks**: about 1.3 s wall each, nine per fit
   (about 12 of 42 s). Per rebuild: allocation 0.40 s (the unit solve is 0.23 s
   of serial R in the master), shear orbit plans 0.28 s, efficient proposals
   0.22 s, chain proposals 0.11 s, Laplace check 0.10 s, collapsed plans
   0.09 s, Laplace refresh 0.24 s. Sample blocks already run straight to the
   next rebuild point.
2. **Kernel rows**: about 7 ms per iteration. The moves evaluate roughly 1,100
   likelihood rows per iteration, more than DMC2's whole 3.4 ms iteration, so
   the comparison that matters is effective samples per second.
3. **Per-call setup in `calc_ll_oo`**: about 27 us for each of about 120 calls
   per iteration (3.3 ms). Even with a parameter-mapping context,
   `make_pt_mapper` is 31% of a one-row call (pre-transform specs, table copy
   and clone, constant columns), and the race branch rebuilds its adapter and
   shared-state objects. Needs an evaluator kept per (data, model) pair.
4. **Shear and collapsed move logic in R**: about 3 ms per iteration.
5. **Gibbs, ensemble bookkeeping, allocation records**: about 1 ms.

## Notes

- The hierarchical burn stage is not reproducible run to run on the R path
  (one chain differs from the first iteration of a block), so burn is checked
  per request, not by whole-run comparison. Adapt and sample are reproducible.
- `native_particle.o` is built without fast-math and without contraction; the
  identity checks depend on it.
- Likelihood entries are keyed on the pointer pair because the same data is
  evaluated with two model objects (the slim worker model and the sampler's).

## Rebuild pass (uncommitted, after 625e9232)

Whole fit 35.9 -> 33.6 s (sample 23.3 -> 22.0 s per 1,000); each rebuild 1.35 -> 1.04 s.
Draws identical to 625e9232 on a 400-iteration sample continuation.

- Fork rounds were job work, not fork overhead (startup about 10 ms a round), so
  merging rounds was dropped.
- `.emc_empirical_covariance_ok`: one symmetric eigen solve instead of two plus an
  `all.equal` symmetry test; same decision on 36,558 test matrices.
- Shear orbit terms: hierarchical density computed per path state for all subjects,
  likelihood through the prepared-data path.
- Allocator: solves repeated within a rebuild are answered from a per-rebuild record
  (108 -> 68 solves over three rebuilds).
- Rebuild dispatcher: a job timed in the calling process decides serial against fork,
  because fork-measured job times run 5-8x above the in-process time for cheap jobs.
  Efficient-proposal conditionals now go through the dispatcher.

Per rebuild now: Laplace refresh 0.20, shear plans 0.20, allocator 0.34 (solves 0.16),
efficient proposals 0.13, fit checks 0.11, collapsed plans 0.09, chain proposals 0.04 s.

## Baseline against DMC2 (measured, not diagnosed)

Whole fits with shear, collapsed, subject HMC and the ensemble update all off, two seeds,
1,000 sample iterations, 25 particles in both packages:

| | sample wall | ESS mu (mu, sigma, tau) | ESS log sd | alpha median / min | move rate |
|---|---|---|---|---|---|
| DMC2 | 3.2-3.5 s | 2913 / 1014 / 1702 | 1750 / 292 / 551 | 1591 / 414 | 0.98-0.99 |
| EMC2, moves off | 9.5-10.3 s | 2501 / 262 / 1179 | 1767 / 53 / 394 | 1562 / 140 | 0.93 |

Resolved: DMC2's default is not Gibbs plus particles. It adds native scale moves and
location moves for standard hierarchies (`options(emc.scale_move = FALSE)` disables them;
`src/sampler/runner.cpp`). With them off DMC2 gives ESS mu 297, log sd 52-53, alpha min
101-131 in 2.4-2.5 s, the same mixing as our moves-off baseline at a quarter of the wall.
Its moves cost about 1 ms an iteration and buy 3.4x (mu) to 5.5x (log sd).

## Exact conjugate check (scratchpad `conj/`)

LNR with only six cell means free (200/200/20/20/2/2 trials a subject, 20 subjects, s and
t0 constant): both packages' likelihoods equal the Gaussian closed form to 1e-6 and the
priors are identical. Reference: 6 x 150,000 exact Gibbs draws (min ESS 3,425; its own
per-chain 5% point of the weakest log sd spans 0.6 sd, which bounds what the check resolves).
3 x 4,000 sample draws each:

| | sample wall | ESS weak log sd (5, 6) | ESS mu6 | worst abs z of means |
|---|---|---|---|---|
| DMC2 default | 27.6 s | 410 / 413 | 1655 | 2.2 |
| DMC2, moves off | 18.2 s | 50 / 49 | 164 | - |
| EMC2, moves off | 50.2 s | 117 / 96 | 262 | 2.9 |
| EMC2 full | 180.4 s | 4949 / 4873 | 8186 | 2.5 |

No arm departs from the reference beyond its noise. DMC2's kernel in adapt and sample is
exact by construction (local auxiliary-centred step or independent step, chosen at random).
