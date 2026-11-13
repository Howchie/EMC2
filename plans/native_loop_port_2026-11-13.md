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
