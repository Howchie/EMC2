# EMC2 architecture and efficiency audit

Status: live implementation plan, 2026-09-10. Architecture review, baseline measurements and evidence checks are complete. Recommendations remain open implementation work, with validation gates below. No package implementation changes are included.

Scope: ordinary BAwL, RDM and DDM fits to empirical data, condition-dependent parameters, hierarchical PMwG estimation, and their memory/parallel execution costs. FPE, GSL integration optimisation, correlated/quadrature-heavy variants, and logical rules are deliberately outside the optimisation programme. Their existing fallback behaviour must survive shared infrastructure changes. The obsolete R likelihoods are not numerical references for this work.

Revision inspected: `6280194719bcc423cfd922fdeedb7d7d7445ddbc`, package 3.4.0. Existing changes to plotting code, documentation and its test are unrelated and have been preserved.

## Decision and priorities

Keep PMwG as the default. There is no evidence from this audit that replacing the sampler wholesale would improve time to a trustworthy posterior. The largest *demonstrated* computational problem is the parameter mapper's 256-cell performance cliff in wide designs. The largest architectural opportunity is using the available cores when there are fewer subjects than cores, without repeatedly forking or compromising worker ownership. For ordinary narrow designs, the native likelihood already dominates: moving just the R sampler loop to C++ cannot produce a large speedup there.

Fix the worker lifecycle/protocol issues before extending parallel execution. These are prerequisites to dependable speed, not optional polish.

| Priority | Work | Evidence / expected benefit | Scope and cost |
|---|---|---|---|
| P0 | Make worker retirement, framing and deadlines complete | Source-confirmed gaps; prevents hangs, silent serial degradation and retained processes | Medium; independent of model numerics |
| P0 | Correct performance telemetry | `worker_max` is currently the maximum **subject** CPU time, not maximum worker workload | Small; needed to rank scheduler changes correctly |
| P1 | Remove the 256-cell mapper cliff; then stop storing scalar coefficients at trial resolution | Identical-result RDM A/B: 54.83 ms cell route versus 113.50 ms fallback; ~2.07x difference | Small/medium for dynamic scratch; larger for compact representation |
| P1 | Add bounded use of spare cores and reduce within-iteration tails | Eight subjects with a 32-core budget and default `r_cores=1` use only eight subject workers | Medium for scheduling; large for a native threaded backend |
| P1 | Reuse parameter-cell kernel preparation | Narrow-design likelihoods spend roughly 95% beyond the mapper; repeated BAwL normalisers and DDM parameter-only terms remain | Medium; model-specific, benchmark before rollout; speedup unmeasured |
| P2 | Make sample history append efficiently and retain suitable worker state across blocks | Whole-history `abind` at each block; worker/context construction repeats | Medium/large; best for long, wide hierarchical fits |
| **P1 (promoted)** | Compile an immutable likelihood context and slim the mixed-data path | Repeated materialisation, partition copies and validation survive outside the finite fast path; also a stated prerequisite of any native threading | Medium; especially small batches and datasets containing omissions |
| P2 | Fuse native proposal-density work where profiling warrants it | Repeated proposal slices, centring, matrix products and mixture arrays | Medium; narrow-model benefit is demonstrably limited |
| P2 experiment | Tune particles/proposals against posterior ESS per second | Potential multiplicative savings, but current particle-weight ESS is not posterior-chain ESS | Must validate inference and adaptation; no promised gain |
| P3 | Broader C++ sampler or alternative MCMC | No present evidence of a net win over targeted changes | Do not start a wholesale rewrite |

Priorities describe benefit on the affected workload. The 256-cell fix is high priority for large models, not a claim that every ordinary eight-parameter fit will become twice as fast.

Priority is not execution order. See "Commit sequence" below for the dependency-ordered commits, the shared mechanisms each one builds or consumes, and which items were merged because they revise the same structure.

## Architecture actually on the fitting path

```text
data + design(formulas, constants, transforms, bounds)
  -> make_emc / design_model / compress_dadm / subject splitting
  -> run_emc: proposal updates, chain blocks, diagnostics, concatenation, saves
     -> run_stages / run_stage: one chain's iteration loop
        -> group Gibbs update + covariance-factor cache
        -> persistent workers assigned subject partitions
           -> new_particle: proposal draws + likelihood + prior/mixture weights
              -> calc_ll_pooled / calc_ll_manager
                 -> calc_ll_oo
                    -> PtMapper / ParamTable / transforms / bounds
                    -> race or DDM native likelihood
        -> native writes into sample arrays
```

Source anchors, relative to the package root:

| Layer | Main locations | Existing optimisations to preserve |
|---|---|---|
| Preparation | `R/fitting.R:803`, `R/design.R:828`, `R/design.R:1371` | Trial/design compression, constants, dropped unobserved coefficients, cached likelihood metadata |
| Chain orchestration | `R/fitting.R:106`, `R/fitting.R:1106` | Chain completion markers and donated core budgets, protection for queued replacement chains |
| Group update | `R/variant_standard.R:228` | Cached design cross-products, Cholesky solves, covariance blocks |
| Particle step | `R/sampling.R:518`, `R/sampling.R:1013` | Block-constant subject factorisations, iteration-level group factors, epsilon rescaling, four-component sample-stage proposal mixture |
| Worker pool | `R/chain_pool.R` | Clean installed-package template, persistent processes, framed FIFOs, shared covariance file, LPT subject assignment, subject-local RNG streams, grow/recycle backoff |
| Likelihood boundary | `R/sampling.R:1641`, `R/oo_map.R:8` | A batch of particles per native call; compressed designs passed through without expansion |
| Parameter mapping | `src/particle_ll.cpp:256`, `src/ParamTable.h` | Mapping/transform plans, invariant columns, row-constant and design-cell paths, reusable bounds/materialisation buffers |
| Race evaluation | `src/particle_ll.cpp:904`, `src/model_BAwL.cpp`, `src/model_RDM_adapters.cpp` | Raw column pointers, winner/loser masks, reusable output buffers, log-space fallbacks, SIMD-aware reductions |
| DDM evaluation | `src/likelihood_ddm.cpp`, `src/model_DDM.h`, `src/ddm_functions_inline.h` | Raw parameter columns, finite/untruncated route, inline Wiener series, endpoint caching |
| History | `src/sample_store.cpp`, `R/sampling.R:1432`, `R/objects.R:516` | Owned arrays and native in-place slice writes; whole-history concatenation remains at block boundaries |

Do not propose “batch the likelihood”, “use column-major indexing”, “cache Cholesky factors”, “use persistent workers”, or “add a finite-trial fast path” as new work: all already exist. Full OpenMP threading is distinct from the existing `omp simd` loops.

## Measurements and their limits

Reproduction helper: [architecture-efficiency-bench.R](architecture-efficiency-bench.R). It prints the loaded/workspace DLL hashes for comparison, validates that direct and managed calls agree, and rejects a fixture whose likelihood is only the invalid-parameter floor.

Environment: Linux x86-64, AMD EPYC Genoa virtual CPU, 32 visible CPUs, R 4.6.1, OpenBLAS, `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`. Generated Makevars reports `legacy`, `march=native`, SIMD enabled. Installed and workspace DLL MD5 both equal `b1a8f056d897e155c097a60dd27aaa50`. Function bodies also match for `calc_ll_manager`, `new_particle`, `run_stage`, `.emc_wpool_iter`, `.emc_wpool_serve` and `.cache_ll_data_attrs`. This checks the actual available binary; it is not a fresh cross-platform rebuild certification.

Fixtures repeat one participant's empirical Forstmann data, add small RT jitter to defeat RT compression, and use plausible parameter clouds with 100 particles. BAwL uses positive leak and the ordinary normal launch without clocks. DDM uses its ordinary zero-`SZ`/zero-`st0` member. These are controlled execution measurements, not fitted/posterior performance or calibration results. Reported times are medians of three repeated batches. This is a shared host, so small differences should not be treated as stable improvements.

| Model / design | Sampled parameters | Augmented rows | Whole call ms | Prologue ms | Approx. prologue share |
|---|---:|---:|---:|---:|---:|
| BAwL, narrow | 8 | 4,000 | 40.58 | 1.75 | 4.3% |
| BAwL, 64 conditions | 70 | 4,000 | 44.67 | 5.42 | 12.1% |
| RDM, narrow | 7 | 4,000 | 37.00 | 1.75 | 4.7% |
| RDM, 64 conditions | 69 | 4,000 | 41.42 | 6.33 | 15.3% |
| DDM, narrow | 6 | 2,000 | 17.33 | 0.92 | 5.3% |
| DDM, 64 conditions | 68 | 2,000 | 18.92 | 2.50 | 13.2% |
| RDM, 256 conditions | 261 | 4,000 | 56.50 | 25.83 | 45.7% |
| RDM, 257 conditions | 262 | 4,000 | 116.50 | 80.83 | 69.4% |

`pt_prologue_oo` includes its own setup and a bounds-counting loop, while the whole call also includes adapter setup, masks and reductions. Therefore “whole minus prologue” is not a pure measurement of transcendental arithmetic. An infinitely fast mapper would improve the narrow BAwL fixture by only about 1.045x; its wide-design benefit can be much larger.

The 256/257-condition rows change the model width by one parameter. A stricter identical-result comparison appends an **unused** row to a 256-row compressed design, triggering the fallback without changing any mapped values, data or particles. This gave **54.83 ms for the cell route versus 113.50 ms for the fallback**, a 2.07x difference, with the complete likelihood vectors bit-identical (`max_ll_delta = 0`). This establishes avoidable overhead in the existing routes; it is not a claim that a new dynamic-scratch implementation has already been built or benchmarked. Repeated wide-model timings varied, especially the prologue, but the large cutoff effect persisted.

A second cutoff probe holds the sampled dimension at seven and varies a numeric covariate's distinct design rows: 256 rows gave 37.17 ms whole / 1.92 ms prologue; 257 gave 39.17 / 3.92 ms. Thus the fixed cutoff is already visible with a small parameter vector; it is much worse with wide design matrices.

For 2,000-trial RDM batches, the managed versus direct native call was approximately 0.75/0.67 ms for one particle, 2.33/2.17 ms for five, and 37.33/36.92 ms for 100. Caching the R wrapper alone is a small win at ordinary batch sizes. Avoid turning this into a large cache-lifecycle project ahead of the measured mapper problem.

A controlled sample-stage `new_particle` profile with all four proposal components, 200 trials, 100 nominal particles and fixed positive-definite proposal matrices spent about 90–95% of the time inside `new_particle` in `calc_ll_oo` (BAwL 2.26/2.37 s, RDM 1.97/2.11 s, DDM 1.00/1.11 s across 500 updates). Rprof samples at 10 ms here. Its surrounding `system.time` triggered extra GC outside the particle loop; that GC is not evidence of per-iteration overhead. The fixture holds proposal settings fixed and does not establish posterior mixing.

An existing hierarchical LBA pool smoke benchmark (24 subjects, 200 trials, 20 preburn iterations, four workers, 100 nominal particles, RT compression at 0.02 s) completed in 2.14 s with no reported pool error. Mean recorded iteration time was 37.3 ms; sends 0.6 ms; shared publication 0.1 ms. This supports keeping the shared-file broadcast rather than assuming it is the main bottleneck. It is a short preburn smoke test, not a BAwL/RDM/DDM ESS benchmark. The profiler's `worker_max` cannot be used for utilisation estimates as currently implemented.

## P0: make worker ownership and protocol recovery complete

These findings are from executable control flow. Existing worker-pool tests pass, including non-CRAN integration tests, but do not prove the missing cases below. No destructive fault injection was performed on other running fits.

1. **Worker ID overflow.** `.emc_wpool_serve` writes `as.raw(msg$w)`; `.emc_wpool_await_done` reads one byte. `.emc_wpool_ll` guards dynamic dispatch with `.EMC_WPOOL_MAX_DYN_WORKERS = 255`, but `.emc_wpool_iter` enables notifications whenever `pool$done` exists, without the equivalent guard. Direct R conversion of 256 or 257 produces raw zero with a warning. A sufficiently wide subject pool can therefore enter the “lost alignment” recovery path. Source: `R/chain_pool.R:614`, `:789-795`, `:990`, `:1300-1308`. Reachability: subject-pool width is bounded by `min(n_subjects, cores)`, so triggering this needs both 256+ subjects and 256+ cores. It is a latent defect, not a live hazard on the target hardware; it is carried free by the single wire-format revision (C4) and should not be scheduled as standalone urgent work.
2. **Deadlines cover only part of the exchange.** The subject pool times the wait for a completion token, then enters a blocking framed receive. A worker that announces completion and stalls part-way through a response can still hang it. Request writes, static likelihood receives and the dynamic likelihood queue's private `await()` also lack equivalent deadlines. A length header and short-read loops solve framing; they do not make I/O bounded. Source: `:77–146`, `:249`, `:980–1134`, `:1348`.
3. **Failed recycling can discard the owner.** On the clean backend, failed respawn returns a dead pool without the template, directory or job handles. `.emc_wpool_stop_workers` does not itself terminate/reap spawned jobs; that responsibility normally belongs to the template. Losing its handle removes the subsequent normal shutdown route. Source: `:668–739`. This is an ownership gap, not an assertion that every recycle failure leaks.
4. **Partial template spawn needs transactional cleanup.** The template accumulates newly forked workers in `made` and commits them to `jobs` only after the loop. A failure after some successful forks can leave those jobs outside the tracked list. Cleanup must own each child immediately. Source: `:810–850`.
5. **Fallback must retire pending work before recomputation.** A dead worker or invalid token marks the pool dead and recomputes outstanding partitions, but does not universally terminate/drain all remaining workers first. Late completion tokens can remain in the channel retained during recycling. Add generation and request identifiers, and retire the whole affected generation before accepting new work. A timeout followed by rerunning a truly nonterminating computation in the chain process is not a guaranteed recovery; preserve the last complete iteration and report the failure instead of hanging the parent.
6. **Count and classify rejected updates.** `safe_new_particle` catches every error and repeats the previous state (`R/sampling.R:878`); the group Gibbs step also catches all errors (`:683`). Intended numerical rejection can therefore look the same as a programming, context or infrastructure error. Record rejection reasons and rates without changing the numerical rejection policy; that half changes no behaviour and belongs with the telemetry schema (C1). Propagating clearly classified infrastructure/programming failures **does** change behaviour: repeating the previous state is what preserves the PMwG invariant under numerical rejection, so a wide fit that previously limped would now stop. Write the policy down before implementing it (C6) — which classified failures abort, whether a chain dies or is marked and continues, what the user sees, and what happens mid-block. Otherwise a fast run that repeatedly copies the old state may be mistaken for efficient sampling.

Implementation contract:

- One pool owner retains all handles even in degraded state; cleanup is idempotent and registered as resources are acquired. Release donated cores only after their workers have stopped.
- Use a framed completion record with a wide worker ID, pool generation, iteration/request ID and status. Keep writes within the platform's atomic pipe bound when sharing a FIFO, and accumulate partial reads correctly. Merely widening `readBin` to four bytes does not establish the protocol.
- Use readiness-driven I/O with monotonic deadlines for sends, headers and payloads. Track progress/heartbeats separately from whether the process is alive. Allow genuinely expensive valid requests; do not assume “20x a historical median can never time out legitimate work”.
- Reap known children with bounded graceful shutdown and escalation for owned jobs; do not signal arbitrary PIDs from a broad process-name search. Preserve failed-pool context until cleanup completes.
- A missing/incomplete result is never committed. Recompute recoverable work from the saved subject RNG stream only after retiring the old assignment; exactly one result updates the iteration.
- Keep safety limits on allocation/frames, but distinguish them from task-size limits. Segment legitimate large messages; validate lengths before allocating. Test both sides of 4 KiB, 8 KiB, 64 KiB, worker 255/256, and maximum supported sizes. The current 2 GiB frame guard is protection, not evidence that all larger work can be handled by a single frame.
- Check available R connections/file descriptors before growing: this pool consumes separate request/reply connections for each worker plus template/completion channels. Ordinary connection limits may prevent reaching 256 workers first; that does not make the one-byte protocol sound on installations with increased limits. Allocate a feasible worker count and queue the remaining work instead of discovering the limit part-way through startup.
- Add fault tests for interrupt/parent exit, death before/after notification, partial payload, partial startup/grow/recycle failure, repeated teardown, stale generation replies, file publication failure and connection exhaustion. Test with a small simulated worker-ID space rather than requiring 256 real processes.

Windows currently parallelises chains through `auto_mclapply`/PSOCK and explicitly rejects `cores_per_chain > 1`. The FIFO pool returns NULL there. Maintain that working route while extending portability; do not present Unix FIFOs or fork as cross-platform solutions.

## P0: repair the measurement contract

`.emc_wpool_compute_particle` records per-subject CPU times. `.emc_wpool_iter` reports `max(times)` as `worker_max` (`R/chain_pool.R:1434`), yet `WorkingTests/bench_worker_pool.R` labels it “slowest worker's own work”. With multiple subjects per worker this comparison understates worker cost, making `response_wait - worker_max` look like transport overhead.

Return both per-subject times and per-worker elapsed/CPU times, plus dispatch/finish timestamps. For the existing static partition, at least report `max(vapply(part, function(s) sum(times[s]), numeric(1)))` as the maximum assigned CPU workload. It still does not equal wall time under descheduling or nested work. Record effective particle counts, active workers, queue delay, recycle/startup time and fallbacks. Measure RSS/PSS or platform-appropriate private memory across the full process tree; summing forked RSS double-counts shared pages. Keep expensive allocation/serialization measurements opt-in.

## P1: remove mapper scaling cliffs, then compact its representation

`src/ParamTable.h:40` fixes `EMC2_PT_MAX_CELLS` at 256. `init_design_plan` admits the cell path only at or below this cap (`:383`); mapping, transforms and bounds all use corresponding stack buffers. Above it, mapping returns to coefficient-by-trial loops and transforms/bounds return to row resolution. This is a correct fallback with a serious performance cliff, not an out-of-bounds finding.

First change: allocate reusable, dynamically sized cell scratch outside the particle loop in **all** consumers (`ParamTable.h`, `transform_utils.cpp`). Select cell computation using actual reuse/cost and a memory budget, not a hidden fixed cell count. Preserve a general row route for continuous covariates with no reuse. Preserve the present coefficient accumulation order, self-referencing intercept semantics, pre-sum/post-sum split transforms, bounds exceptions and absent-cell representatives. Cover 255, 256, 257 and much larger cells in equivalence tests. Do not only raise the constant to 1,024 and move the problem.

Second change: for the no-trend path, represent sampled coefficients as scalars rather than repeating every coefficient down a full trial-length base column. `fill_from_particle_row_planned` currently fills those columns on every particle (`ParamTable.h:959`), even when their only purpose is forming a small number of model-parameter columns. Separate coefficient storage from natural parameter outputs. Compute outputs at cell resolution and retain existing contiguous model-parameter columns initially, preserving the raw-kernel interface. This removes trial-count times coefficient-count storage/traffic without forcing gather accesses into every kernel.

In the same change, add a model-agnostic **common-refinement** operation over design cells. Each `DesignEntry` already induces a partition of the trial index set: `expand_idx` maps trial to design row, so trials reading the same row share a cell (`ParamTable.h:355-393`). The mapper holds one partition per output column and has no operation to combine two of them. Add one: given a **set of output column indices**, return the partition in which trials `i` and `j` share a cell iff they agree on `expand_idx` for every column in the set, as a dense cell id per trial plus representatives.

This is set arithmetic on integers. `ParamTable` receives an index set and returns a partition; it learns nothing about what any parameter means, which is required — `src/col_registry.h` carries roughly forty model namespaces whose parameters have no vocabulary in common, and `DesignEntry` describes an output only as `out_idx`. Two properties matter. The refinement is bounded by `min(prod(n_cells), T)`, so a continuous covariate drives it to `T`, meaning no reuse: the operation must report that and defer to the general row route, which is the same budget decision this section already makes. And it depends only on design and data, never on parameter values, so it is computed once per dadm and cached rather than per particle.

Scalar coefficient storage and cell-resolution output are the same mechanism; implement them together (C8) rather than as two prototypes.

Computation granularity and storage granularity are separable, so this is not a choice between cell gathers and contiguous columns. The shape to adopt throughout is **compute at cells, materialise at rows**: evaluate the expensive parameter-only work once per cell, expand cells to rows into flat scratch (a pure copy with no dependent arithmetic, which vectorises), then run the consuming loop contiguously with SIMD intact and no gather in the hot path. The expansion pass is the price of keeping the hot loop contiguous, and it is usually worth it.

Tile the pre-pass rather than materialising whole trial-length arrays. Each derived quantity held at row resolution is a length-`T` array written and then read; several of them at once become extra memory streams and can turn an ALU saving into a bandwidth loss. Process trials in chunks sized to keep the scratch L1-resident, running pre-pass then consumer per chunk.

A generic sparse matrix conversion or BLAS multiply is not automatically faster than the existing small cell loops. Do not compute all particle × trial × parameter outputs at once.

Structural dependencies must govern skipping: a blocked proposal repeats unchanged coordinates in every proposal row, but `make_pt_mapper` identifies sampled coefficients from column **names**, not a varying-coordinate mask (`src/particle_ll.cpp:490`). Passing an explicit, verified update mask would allow unchanged mappings within a proposal batch to remain fixed. Invalidate between blocks and group/subject states as appropriate; never infer invariance from the first two equal floating-point values. Trend runtimes currently disable this optimisation for a reason and must retain the general route until their dependency contracts are explicit.

## P1: use spare capacity safely; native likelihood threads are a staged option

The current scheduler already reallocates cores across finishing chains and balances subjects with LPT. Remaining gaps:

- `.particle_core_budget(8, 32, 1, 32)` returns eight subject workers and one likelihood core. With multiple subjects and default inner width, spare cores are unused even if likelihoods are expensive. Donated cores also stop helping when all subjects have workers.
- `ctx$r_cores` is fixed for the pool lifetime. Explicit inner parallelism goes through per-call `auto_mclapply` in `calc_ll_manager`; there is a persistent particle likelihood pool specifically for the one-subject case, not a reusable nested pool for every subject.
- LPT assignments are fixed within an iteration. A worker that finishes its partition cannot take a remaining subject from another partition. Previous subject CPU cost is useful but does not predict every adaptive particle count or unusually expensive proposal cloud.

Implement a single bounded task budget. Start with dynamic whole-subject assignment only where measured tail imbalance exceeds dispatch costs; publish group state once per iteration and retain/reuse worker-local factors across its tasks. Preserve static LPT for cheap uniform work. For fewer subjects than useful cores, allow bounded particle tasks using the same scheduler or an explicitly bounded native executor, without each subject independently starting another full pool. Commit the group update only after the existing iteration barrier. Do not run later iterations with stale group draws.

Preserve subject-local L'Ecuyer streams for proposal generation, selection and adaptation. Deterministic likelihood tasks carry particle indices, draw no random numbers and write ordered result slots. Core allocation includes BLAS/native threads and useful parent computation; it is not enough to count only R workers. Inspect available affinity/quota/memory in defaults and respect the user's explicit budget. Set thread policy at supported process/library boundaries; environment variables changed after a BLAS library initialises are not a universal control mechanism.

A native executor is feasible without changing the public R API, but **not** by placing `omp parallel for` around `calc_ll_oo`'s existing particle loop. `PtMapper` allocates/mutates Rcpp objects; bounds return R vectors; trend plans can call R; contexts and some caches/debug counters are mutable. Build an immutable, owning data/design plan and private native workspaces first — that is C11, and it is a prerequisite, so filing it under P2 below inverts the dependency; it is promoted in the commit sequence. Then run a serial implementation against the current baseline, and only then parallelise independent particles with a portable library and a serial fallback. Audit every reachable Rmath fallback rather than assuming all double-valued kernels are thread-safe. Keep R allocation, RNG, callbacks and error/interrupt reporting on the main thread. RcppParallel's [thread-safety documentation](https://rcppcore.github.io/RcppParallel/) explicitly restricts worker access to the R/Rcpp API; its [portability guidance](https://rcppcore.github.io/RcppParallel/tbb.html) describes supported backends and the need for fallback.

Use bounded scratch per worker, not per particle, and memory-aware concurrency. Test Windows, macOS ARM/x86 where available, and Linux. Defer this larger migration if ordinary subject parallelism already saturates the user's machines and the measured gain cannot justify it.

## P1: eliminate repeated kernel preparation without changing formulas

For narrow designs, focus on actual numerical kernels rather than the wrapper:

- **BAwL:** `natural_normalizer` recomputes `pnorm(v/sv)` for repeated parameter cells (`src/model_LBA.h:120`, `:233`, `:292`). Prepare its exact denominator, validity/fallback decision, and reusable parameter-only geometry per particle and cell. Reuse for trials with identical complete relevant inputs. Keep the existing natural/log-space acceptance tests and `BAWL_DENOM_FLOOR`; retain the exact zero-leak and small-start-range limits. RT-dependent exponentials remain RT-dependent.
- **RDM:** `drdm_raw` and `prdm_raw` repeatedly form scale/geometry terms (`src/model_RDM_adapters.cpp:38-48`, `:59-69`). Their preparation is **identical** — same `inv_s`, same three scaled arguments — but each recomputes it on a complementary mask, density for winners and survivor for losers. One shared pre-pass serving both consumers therefore computes it once for all rows, roughly halving the parameter-only work independently of any cell structure, and removes the masked-row question entirely: a pass feeding both kernels has no wasted rows. Do this before, and separately from, any cell-resolution work, since it needs neither the refinement nor a reuse key. Note the arithmetic saved here is a division and a few multiplies against a transcendental density, so the ceiling is low; the structural tidiness is worth more than the cycles.
- **DDM:** the default raw path repeatedly scales `a`, `v`, `sv`; `dwiener_inline` recomputes parameter-only terms such as `log(a)`, and the long-time series evaluates `sin(k*pi*w)` for shared bias/response cells (`src/model_DDM.h:92`, `src/ddm_functions_inline.h:49`, `:175`). Cache exact reusable terms only where repeated use pays. Keep RT-dependent series-length selection, error tolerance, reflection by response and summation order. Do not truncate the series to a fixed length or replace sine with an unvalidated recurrence.

Reuse keys must come from the mapper's common-refinement operation above, taking the refinement over exactly the columns the reusable subexpression reads, rather than from a per-model key scheme. Different parameter formulas induce different cell partitions, so keying on any single parameter's cells is wrong whenever another column feeding the same subexpression varies independently of it; the refinement is the general statement of that condition, and it is why this work is sequenced after C8 rather than beside it. Declaring a narrower column set than the kernel receives is an optimisation, not a requirement — the full set is always correct. Exact complete parameter tuples remain the fallback key where structure gives none. Trends/continuous covariates require the general route when reuse disappears. Bound context growth and fall back safely.

This is a justified prototype target, not a measured speedup. Supporting prior: the BAwLcorr preburn fix already won ~2.6x by removing exactly this class of per-trial recomputation (per-trial `q_AB` and truncation-Z), which raises the prior on the BAwL item without making it a measured result here. Add counters or native profiles for time spent in these terms and their reuse rates, using the shared profile schema rather than ad-hoc counters; abandon each specialisation if setup/indirection outweighs saved arithmetic. Require a repeatable end-to-end benefit, not merely a faster isolated helper.

## P2: history, block setup and diagnostics

Native slice writes already avoid full-array copying every iteration. The remaining history issue is `concat_emc` (`R/objects.R:516`), which rebuilds sample arrays with `abind` every block. With a fixed block length and retained history this repeatedly copies increasing prefixes: total copying grows quadratically in the number of blocks. `run_emc` also subsets state, refreshes proposals, restarts chain processes/worker templates and optionally serialises the whole fit at block boundaries.

For a standard hierarchy, just `alpha` and `theta_var` require about `8 * I * (P*N + P^2)` bytes per chain, before copies and other arrays. At P=200, N=140, I=10,000 that is 5.44 GB per chain. This calculation is a capacity estimate, not a measured resident-memory result.

Use internal appendable chunks or explicitly owned growable capacity; materialise the familiar arrays at public API/checkpoint boundaries. A fixed-capacity internal array must distinguish allocated and committed iterations so diagnostics cannot consume uninitialised samples. Preserve names, nuisance stores, serialization/resume, and copy isolation when users keep earlier fit objects. Simply mutating previously returned arrays would break R expectations. Audit `extend_sampler`'s last-dimension identification if capacity differs from logical length.

Where block-start costs are material, keep clean worker templates/context alive within a fit and explicitly version updates to proposals, stages and transforms. Preserve coordinator barriers for cross-chain diagnostics/proposals. Do this only after P0 ownership fixes; longer-lived workers increase the consequences of lifecycle bugs. Choose recycle intervals from measured private-memory growth and restart cost with a hard memory fallback; never disable recycling blindly.

Profile `create_chain_proposals`, `create_eff_proposals`, `check_progress` and save/concatenate time as separate block costs. They are absent from the per-iteration timing. Do not change diagnostic frequency or ESS/Rhat criteria silently to improve a timing number. Incremental checkpoints may be useful for large fits, with atomic publication and a supported public resume format.

## P2: compiled context and mixed data

The finite/untruncated route requires *all* relevant trials to qualify. A dataset containing omissions or truncation takes the mixed path, though its finite trials are still batched; they have not all reverted to scalar integration. Existing shared state and endpoint grouping must be credited.

Remaining work includes `pt.materialize_reusable()` before ordinary mixed race evaluation (`src/particle_ll.cpp:1194`), and assignments that **copy** cached `std::vector` partitions into local vectors each particle (`:1807–1809`), despite the nearby comment saying “no copy”. Take const views/references with explicit lifetime, pass native parameter-column views into the mixed path, and retain shared reusable scratch. Keep trial classes separate so a small omission subset does not impose avoidable housekeeping on the finite majority. Preserve global trial validity and time/nogo/RACE rules; do not simply concatenate independently evaluated subsets when normalisation depends on their full accumulator block.

The R manager calls the model closure and cache validation each time, constructs a design list, and even `r_cores=1` normally enters the split/subset/lapply branch when there is more than one proposal (`R/sampling.R:1668`). A direct serial branch removes an unnecessary proposal matrix copy and list round-trip. This is a small, straightforward first patch, not the largest speedup.

A longer-lived immutable likelihood context can avoid repeatedly constructing `PtMapper`, adapter dispatch, transform/bound specifications and data masks. **This is promoted to a P1 prerequisite (C11):** the native-threading section requires it first, and the mixed-path fix above needs the same lifetime discipline, so the two are one commit rather than two independently invented ownership schemes. Data-independent compiled plans must be separate from per-call values, invariant masks and bound seeds. The present mapper caches values tied to the first particle; persisting that mutable state unchanged across calls would give wrong likelihoods. Use per-worker ownership, explicit invalidation for data/levels/design/constants/transforms/bounds/trends/guess-window changes, and rebuild after deserialization. Do not serialize an external pointer and expect it to work in a spawned process. Public low-level calls with mutable dadms still need validation; an immutable internal context can validate once.

Exact compression remains useful, but changing `rt_resolution` changes the represented data. Never recommend coarser RT rounding as a numerically neutral optimisation. For existing exact duplicate trials, precomputed multiplicities can replace a long expansion-index sum in the total-likelihood API; this changes reduction order and requires an explicit error budget. Keep expansion/order for pointwise likelihoods and information criteria. Its benefit is small unless duplicates are numerous and reductions are a measured cost.

## P2: proposal computation and statistical efficiency

`new_particle` repeatedly slices the same particle block for `fast_dmvnorm_rooti`, which centres and multiplies into fresh native matrices (`src/fast_dmvnorm.cpp`). A single native proposal-density/mixture entry point can reuse centred scratch, covariance factors, output arrays and index maps, returning the log mixture and prior density. Keep proposal generation/selection in R initially to preserve RNG order and make equivalence reviewable. Sample-stage chain-centred proposals share a covariance source; some transformed centred terms can be reused, but epsilon/mean changes must be applied correctly.

For multiple proposal blocks, the full prior path calls `fast_dmvnorm` again and factorises the same group covariance (`R/sampling.R:1193`). Cache its full factor once per group draw/worker as well as component factors, without replacing the full correlated prior by independent component priors. Retain current singular-matrix handling and distinguish proposal regularisation from the actual prior distribution. Do not cache prior densities across changing group draws.

A full native sampler migration would additionally require matching component bookkeeping, numerical rejection, nuisance/joint/group-design support, adaptation, RNG consumption and sample-store ownership. The narrow-model profile limits what moving only orchestration can save. Benchmark wide blocked fits before expanding the native scope.

`update_pm_settings` already adapts particle counts using particle-weight ESS after `gd_good`, with a minimum of 25 and a cap derived from initial particles (`R/sampling.R:1356`). That quantity measures concentration within one proposal cloud, not autocorrelation of the saved posterior chain. Evaluate particle counts, mixture quality, blocking and parameterisation by worst-relevant-parameter bulk/tail ESS per elapsed second and total time to the same convergence criteria. Include between-subject covariance and condition contrasts, not just population means or mean acceptance.

Any changed count/adaptation rule needs the PMwG invariant-target argument and a frozen-versus-adaptive validation strategy. Do not prune particles based on cheap preliminary scores, omit mixture density corrections, remove the retained particle, or casually change state-dependent adaptation. Reparameterisations can change the prior unless transformed with the correct Jacobian; compare the same intended model.

HMC/NUTS is not a drop-in acceleration: these likelihoods are double-valued C++ with support boundaries, floors and branch/series decisions, not an automatic-differentiation target. Finite-difference gradients multiply expensive likelihood calls with dimension. A separate differentiable backend could be researched for a restricted model, but there is no evidence here to recommend it as the package default. GPU execution, approximate surrogates and wholesale sampler replacement are also outside the actionable near-term plan.

## Numerical and portability acceptance gates

1. Same-build, pure scheduling/storage/cache changes should reproduce particle likelihood vectors and subject streams exactly. Differential tests must include repeated calls, reordered particles, chunk boundaries, worker growth/recycling/fallback and restarted fits.
2. Mapper-only changes must match mapped natural parameters and bounds against the existing independent `get_pars_c_batch_wrapper_oo` route. Test constants with duplicate names/last-writer precedence, split transforms, missing design cells, one-column dimensions, fixed versus sampled parameters, and trend fallback.
3. A pre-pass or hoist must adopt **one** existing operand order explicitly and state which. `src/model_RDM_adapters.cpp` already carries two spellings of the same quantity — `par[1]*inv_s + 0.5*par[2]*inv_s` on the scalar path (`:14`) versus `(B + 0.5*A)*inv_s` on the raw path (`:44`) — which are algebraically equal and not bit-identical, and which `legacy` fast-math may reassociate differently again. Whichever a hoist adopts, the other path shifts by an ulp and the bit-identical gate rejects it. Migrate the other path deliberately as its own reviewed change, or leave it untouched. This is the same class as the earlier cell-accumulation FMA mismatch.
4. Arithmetic changes must preserve support, floor and mixture semantics. Measure per-trial **and** total log-likelihood error over ordinary and adversarial grids; specify tolerances before accepting results and scale total error checks with trial count. Near-zero probabilities require log-space/absolute checks rather than unstable relative error. Independent high-precision/analytic tests and the existing C++ baseline are references; stale R likelihood code is not.
5. Cover RT near t0, zero/small A and leak, zero drift where supported, extreme tails, invalid proposals, NA versus NaN, ±Inf censoring, known/unknown responses, LT/UT endpoints, varying accumulator counts/inactive rows, contaminants/guessing and compressed/uncompressed expansion. Shared code must retain excluded models' fallbacks.
6. No new relaxed-math defaults, reduced precision, changed pnorm approximations, increased integration tolerances or RT binning. Existing `legacy` uses fast-math with finite/NaN protections; `compatibility` remains a necessary comparison. `-march=native` is useful for local builds but is not a portable distributed binary. Optional ISA specialisations need runtime dispatch and a baseline build, not unconditional AVX instructions.
7. Validate GCC/Clang/Rtools and macOS/Linux/Windows builds, serial execution without a native threading backend, multiple core budgets, BLAS oversubscription, allocation failure and connection limits. Existing portability is a constraint, not a speed tradeoff.
8. Compare full hierarchical runs using the same data, priors, convergence criteria and wall-clock boundaries. Include start-up, all stages, proposal updates, diagnostics, checkpointing and cleanup. Report peak total private memory as well as ESS/s. Short timing loops establish implementation costs, not posterior equivalence or convergence.

## Commit sequence

Each numbered item below is one reviewable commit. Dependencies are stated; commits sharing a
dependency set can proceed in parallel. Three shared-infrastructure commits land first because
five later sections would otherwise each grow their own version of the same mechanism.

Cross-cutting mechanisms deliberately built **once**:

| Mechanism | Built in | Consumed by | Would otherwise be reimplemented in |
|---|---|---|---|
| Owning, dynamically sized scratch with a per-call/per-worker lifetime | C7 | C8, C9, C11, C15 | mapper cells, kernel preparation, native workspaces, mixed-path buffers, proposal centring |
| Common refinement of design-cell partitions over a caller-supplied column set | C8 | C9, C10 | three separate model-specific reuse keys |
| Profile record schema and its single reporting path | C1 | C9, C13, C14, C3 | worker telemetry, rejection counting, kernel reuse rates, block costs |
| Immutable owning native context and explicit invalidation | C11 | C13 native option, C15 | mixed-path views, threading prerequisites, compiled-context work |
| Framed record fields (worker ID, generation, request ID, status) | C4 | C5, C13, C16 | two independent wire-format revisions |

### Foundation (blocks everything below)

**C1 — One profile record schema.** Define a single extensible profiling record with documented
field names and one reporting path used by both `WorkingTests/bench_worker_pool.R` and
`bench_large_fit_sampling.R`. Fix the measurement contract in the same commit: report maximum
assigned worker workload, `max(vapply(part, function(s) sum(times[s]), numeric(1)))`, under a name
distinct from the per-subject maximum, and keep both. Add dispatch/finish timestamps, per-worker
elapsed and CPU time, effective particle counts, active workers, queue delay, recycle/startup time,
fallback counts and tree-aware private memory (RSS/PSS, not summed forked RSS). Add rejection
counters and reason classification here, where they change no behaviour. Keep expensive
allocation/serialization measurements opt-in. No behaviour change; no dependencies.

**C2 — Equivalence and fault test harness.** The matrix in "Numerical and portability acceptance
gates" is demanded by nearly every commit below and must exist before the first of them, not inside
whichever lands first. Build the parameterised differential fixture (repeated calls, reordered
particles, chunk boundaries, worker growth/recycling/fallback, restarted fits), the mapper
comparison against `get_pars_c_batch_wrapper_oo` (255/256/257 and much larger cells, duplicate
constant names, split transforms, missing cells, one-column dimensions, fixed versus sampled, trend
fallback), the numerical grid (RT near t0, zero/small A and leak, zero drift, extreme tails, invalid
proposals, NA versus NaN, +/-Inf censoring, LT/UT endpoints, varying accumulator counts,
contaminants, compressed/uncompressed) and the fault cases (interrupt/parent exit, death
before/after notification, partial payload, partial startup/grow/recycle failure, repeated teardown,
stale generation replies, publication failure, connection exhaustion). Use a small simulated
worker-ID space rather than 256 real processes. No dependencies.

**C3 — Hierarchical regression benchmark harness.** Sequence items below and every "accept only
reproducible gains" gate need matched hierarchical runs: N=8/32/140+, P~8/64/250+, uneven trial
counts, multiple conditions, varying particle counts, multiple chains, and 1/2/4 plus larger core
budgets rather than all detected cores. Report worst-relevant-parameter bulk/tail ESS per elapsed
second, peak total private memory, and full wall clock including start-up, all stages, proposal
updates, diagnostics, checkpointing and cleanup. Consumes C1's schema. Depends on C1.

### P0 — worker ownership, protocol and measurement

**C4 — Revise the framed record once.** Findings 1, 2 and 5 all change the same wire format;
each revision re-runs the whole boundary test matrix, so revise it exactly once. Carry a wide worker
ID, pool generation, iteration/request ID and status together. Add readiness-driven I/O with
monotonic deadlines to sends, headers and payloads, including the dynamic likelihood queue's private
`await()` and the static likelihood receives, and track progress/heartbeats separately from
liveness. Keep writes within the platform's atomic pipe bound, segment legitimate large messages and
validate lengths before allocating. Test both sides of 4 KiB, 8 KiB, 64 KiB and simulated worker
255/256. Depends on C2.

*Reachability note for finding 1:* subject-pool width is bounded by `min(n_subjects, cores)`, so
reaching worker 256 needs both 256+ subjects and 256+ cores. The defect is real and correctly
diagnosed, but it is latent on the target hardware. It rides along free in C4; do not schedule it as
standalone urgent work.

**C5 — Ownership and generation retirement.** Findings 3 and 4 are one concern: cleanup must own
each child at the moment it is acquired. Retain template, directory, library and context handles
even in degraded state; make cleanup idempotent; release donated cores only after their workers
stop; reap owned jobs with bounded graceful shutdown and escalation, never by process-name search.
Add the retirement half of finding 5: retire the whole affected generation before accepting new
work, preserve the last complete iteration and report rather than hang, and let exactly one result
update the iteration. Generation retirement is meaningless without C4's generation field. Depends on
C4.

**C6 — Failure classification policy.** Counting lands in C1 and changes nothing. Propagation
changes behaviour and needs its policy written down before it is implemented: `safe_new_particle`
repeating the previous state (`R/sampling.R:878`) is what preserves the PMwG invariant under
numerical rejection, and the group Gibbs catch-all (`:683`) does the same. Specify which classified
failures abort, whether a chain dies or is marked and continues, what the user sees, and what
happens mid-block. This is the same mechanism that previously masked a diverging solve as
"efficient" sampling, so leaving the policy implicit repeats that. Depends on C1.

### P1 — mapper, kernels, scheduling

**C7 — Dynamic cell scratch; remove the 256 cliff.** Cover both consumers in one commit:
`ParamTable.h:730-732` (`post_val`, `pre_val`, `self_val`) and `transform_utils.cpp:124,195,263,283`
(`ok_cell`, `v`). Introduce here the owning, dynamically sized scratch with an explicit per-call
lifetime that C8, C9, C11 and C15 all reuse; this is the single scratch primitive, not a
mapper-local fix. Select cell computation from actual reuse/cost and a memory budget, not a hidden
fixed count, and preserve a general row route for continuous covariates with no reuse. Preserve
coefficient accumulation order, self-referencing intercept semantics, pre-sum/post-sum split
transforms, bounds exceptions and absent-cell representatives. Do not merely raise the constant.
`ParamTable.h` reaches roughly a dozen translation units transitively through `transform_utils.h`,
`likelihood_ddm.h` and `TrendEngine.h`; header dependency tracking already exists in
`src/Makevars.in`, but this is a further reason to revise the layout once rather than across two
commits. Gate: C2 equivalence bit-identical, plus a re-run of the identical-input cutoff benchmark.
Depends on C2, C7 has no other dependency.

**C8 — Joint cell partition, cell-resolution outputs and scalar coefficients.** Each `DesignEntry` induces a
partition of the trial index set via `expand_idx`, but the mapper holds one per output column and
cannot combine them (`src/ParamTable.h:355-393`). Add a model-agnostic common-refinement operation:
it takes a **set of output column indices** and returns the partition in which two trials share a
cell iff they agree on `expand_idx` for every column in the set. `ParamTable` stays agnostic — it
receives an index set, never a meaning. The subset choice is the only model-dependent part and
already lives in a model-specific translation unit; the safe default is every column the kernel is
handed, which each model already declares through its `emc2col::<model>` enum / p_types order, and a
narrower declared subset is an opt-in coarsening that is never required for correctness. Bound the
result by `min(prod(n_cells), T)` and defer to the row route when a covariate drives it to `T`.
Compute it once per dadm and cache it: it depends on design and data only, not on parameter values.
Without this one operation, every model that wants reuse builds its own key. In the same commit stop storing scalar coefficients at trial resolution:
`fill_from_particle_row_planned` currently runs `std::fill(col, col + T, val)` for coefficients whose
only purpose is forming a small number of model-parameter columns (`ParamTable.h:959-981`). Separate
coefficient storage from natural-parameter outputs, compute outputs at cell resolution and retain
the existing contiguous model-parameter columns so the raw-kernel interface is unchanged. Scalar
storage and cell-resolution output are the same mechanism; splitting them writes that path twice.
Do not compute all particle x trial x parameter outputs at once. Depends on C7.

**C9 — Kernel preparation: measure, then specialise.** First deliverable is instrumentation using
C1's schema: reuse rate and time share for BAwL `natural_normalizer` (`src/model_LBA.h:120,233,292`),
RDM scale/geometry in `drdm_raw`/`prdm_raw` (`src/model_RDM_adapters.cpp:30`), and the DDM
`a/s`, `v/s`, `sv/s` scaling (`src/model_DDM.h:102,121,131,169-184`), `std::log(a)`
(`src/ddm_functions_inline.h:156,191,220`) and `std::sin` in the series (`:63,159`). Then one commit
per model, each declaring the column set its reusable subexpression reads and consuming C8's
refinement and C7's scratch rather than its own key scheme;
abandon each individually if setup and indirection outweigh saved arithmetic. Preserve expression
order before any algebraic rearrangement, keep RT-dependent series-length selection, error
tolerance, reflection by response, summation order, the natural/log-space acceptance tests and
`BAWL_DENOM_FLOOR`. Supporting prior: the BAwLcorr preburn fix already won ~2.6x by removing exactly
this class of per-trial recomputation (per-trial `q_AB` and truncation-Z), which raises the prior on
the BAwL item specifically without making it a measured result here. Depends on C1, C8; C3 for the
end-to-end gate.

**C10 — Explicit update mask for blocked proposals.** `make_pt_mapper` identifies sampled
coefficients from column **names** (`src/particle_ll.cpp:490-500`), not a varying-coordinate mask,
so a blocked proposal that repeats unchanged coordinates still re-maps them. Pass an explicit,
verified mask and invalidate between blocks and group/subject states. Never infer invariance from
equal floating-point values. Trend runtimes retain the general route. Needs C8's cell-resolution
invalidation. Depends on C8.

**C11 — Immutable owning native context. Promoted out of P2.** The native-threading section says to
build this "first" and the mixed-data section needs the same lifetime discipline, so it is a
prerequisite rather than a later nicety; scheduling it as P2 inverts the dependency. Give
`PtMapper`, adapter dispatch, transform/bound specifications and data masks an explicit owner and
lifetime, separating data-independent compiled plans from per-call values, invariant masks and bound
seeds. The present mapper caches values tied to the first particle, so persisting that mutable state
unchanged across calls would give wrong likelihoods. Use per-worker ownership, explicit invalidation
for data/levels/design/constants/transforms/bounds/trends/guess-window changes, and rebuild after
deserialization; never serialize an external pointer and expect it to survive a spawned process.
Absorb the mixed-path fix here rather than separately: replace the `std::vector` partition
assignments at `src/particle_ll.cpp:1804-1809`, whose own comment already claims "no copy", with
const views of explicit lifetime, and address `pt.materialize_reusable()` before ordinary mixed race
evaluation (`:1194`). Both are "give native state an owner"; done apart, the mixed-path fix invents a
lifetime scheme this commit then replaces. Keep trial classes separate and preserve global trial
validity and time/nogo/RACE rules. Depends on C7, C2.

**C12 — Direct serial branch in `calc_ll_manager`.** `nrow(proposals) <= r_cores`
(`R/sampling.R:1668`) means even `r_cores=1` enters the split/subset/lapply branch whenever there is
more than one proposal. A direct serial branch removes a proposal matrix copy and a list round-trip.
Small, independent, safe to land at any point; explicitly not the largest speedup. Depends on C2.

**C13 — Bounded task budget and spare-core use.** `.particle_core_budget(8, 32, 1, 32)` returns
eight subject workers and one likelihood core, leaving 24 cores idle. Introduce one bounded task
budget: dynamic whole-subject assignment only where measured tail imbalance exceeds dispatch cost,
static LPT retained for cheap uniform work, group state published once per iteration, worker-local
factors reused across tasks, and the group update committed only after the existing barrier.
Preserve subject-local L'Ecuyer streams; deterministic likelihood tasks carry particle indices, draw
no random numbers and write ordered slots. Count BLAS/native threads and useful parent computation
in the allocation, and set thread policy at supported process/library boundaries rather than through
environment variables set after a BLAS library initialises. This commit cannot be ranked without
C1's corrected metric, and dynamic dispatch amplifies lifecycle bugs. Depends on C1, C3, C4, C5.

The native threaded executor remains a **staged option after** C11 and C13, not part of them: build
the serial implementation against the current baseline first, then parallelise independent particles
with a portable library and a serial fallback, auditing every reachable Rmath fallback and keeping R
allocation, RNG, callbacks and error/interrupt reporting on the main thread. Defer it entirely if
ordinary subject parallelism already saturates the target machines.

### P2 — history, proposals, statistical efficiency

**C14 — Appendable sample history.** `concat_emc` rebuilds sample arrays with `abind` every block
(`R/objects.R:523`, and again for the nuisance store at `:531`), so retained history copies growing
prefixes and total copying grows quadratically in blocks. Use internal appendable chunks or owned
growable capacity, distinguishing allocated from committed iterations so diagnostics cannot consume
uninitialised samples, and materialise the familiar arrays at public API and checkpoint boundaries.
Preserve names, nuisance stores, serialization/resume and copy isolation; never mutate previously
returned arrays. Audit `extend_sampler`'s last-dimension identification if capacity differs from
logical length. Depends on C2, C3.

**C15 — Native proposal-density entry point.** Reuse centred scratch (from C7), covariance factors,
output arrays and index maps behind one entry point returning the log mixture and prior density, and
cache the full group covariance factor once per group draw/worker alongside component factors
(`R/sampling.R:1193` refactorises it again for multiple blocks). Keep proposal generation and
selection in R initially to preserve RNG order. Retain singular-matrix handling, distinguish
proposal regularisation from the prior itself, and never cache prior densities across changing group
draws. Narrow-model benefit is demonstrably limited; gate on wide blocked fits. Depends on C7, C11,
C3.

**C16 — Retain worker templates and context across blocks.** Only after C4 and C5: longer-lived
workers increase the consequences of lifecycle bugs. Version updates to proposals, stages and
transforms explicitly, preserve coordinator barriers for cross-chain diagnostics, and choose recycle
intervals from measured private-memory growth and restart cost with a hard memory fallback. Profile
`create_chain_proposals`, `create_eff_proposals`, `check_progress` and save/concatenate as separate
block costs using C1's schema; they are absent from per-iteration timing today. Depends on C1, C5,
C14.

**C17 — Particles and adaptation experiment.** `update_pm_settings` adapts particle counts from
particle-weight ESS, `sum(weights)^2 / sum(weights^2)` (`R/sampling.R:1356`), which measures
concentration within one proposal cloud rather than autocorrelation of the saved chain. Evaluate
counts, mixture quality, blocking and parameterisation by worst-relevant-parameter bulk/tail ESS per
second and time to the same convergence criteria, including between-subject covariance and condition
contrasts. Any changed rule needs the PMwG invariant-target argument and frozen-versus-adaptive
validation. Do not prune particles on cheap preliminary scores, omit mixture density corrections,
remove the retained particle or casually change state-dependent adaptation. No promised gain.
Depends on C3.

### Stop rules

Abandon a candidate if a faster iteration loses its gain to worse mixing, setup cost, memory
pressure or repeated fallback. Accept C7-C10 only on bit-identical likelihood vectors, and C9's
specialisations only on repeatable end-to-end benefit rather than a faster isolated helper. Do not
make the native rewrite a prerequisite for the smaller demonstrated improvements, and do not start
a wholesale sampler rewrite (P3) at all.

Baseline checks completed: `test-proposal-cache.R`, `test-param-table-prologue.R`,
`test-particle-core-budget.R`, `test-sample-store-inplace.R`, `test-ll-data-cache.R`, and
`test-worker-pool.R` with `NOT_CRAN=true` passed. The intentional recycle-failure test emits its
degradation warning. The cache tests initially lacked the namespace parent when invoked individually
and passed when invoked in the correct namespace environment; that was a test invocation issue, not
a package defect. Cross-platform builds, new failure-injection cases, and full posterior ESS
comparisons have not been run in this audit.

### Reproduction commands

Run from the repository root against a matching installed package. Run timing workloads separately to avoid competing benchmarks; compare the printed DLL hashes and active build configuration first.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript plans/architecture-efficiency-bench.R
AUDIT_MODE=wide AUDIT_REPS=6 Rscript plans/architecture-efficiency-bench.R
AUDIT_MODE=particle Rscript plans/architecture-efficiency-bench.R
EMC_N=24 EMC_TRIALS=200 EMC_ITER=20 EMC_WORKERS=4 EMC_PARTICLES=100 EMC_PROFILE=true Rscript WorkingTests/bench_large_fit_sampling.R
```

For the focused baseline tests:

```r
library(EMC2)
Sys.setenv(NOT_CRAN = "true")
files <- c("test-proposal-cache.R", "test-param-table-prologue.R",
           "test-particle-core-budget.R", "test-sample-store-inplace.R",
           "test-ll-data-cache.R", "test-worker-pool.R")
for (f in files) {
  testthat::test_file(file.path("tests/testthat", f),
    env = new.env(parent = asNamespace("EMC2")),
    reporter = "summary", stop_on_failure = TRUE)
}
```
