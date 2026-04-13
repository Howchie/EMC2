# Integration Plan: BAwL + Refactored RDM into `dev-oo` (while preserving speed work)

## Goal
Integrate from `/data/work/EMC2_zachdev` (local `dev_extended`) into `/data/work/EMC2_dev_oo`:
1. **BAwL** (LBA with optional leak `k`).
2. **Refactored RDM** with canonical Wald/SWTN forms and optional between-trial drift variability (`sv`).

Constraint: keep `dev-oo` performance architecture intact (OO + SIMD + raw fast paths + fast pnorm + cens/trunc speedups).

## What I inspected

### `EMC2_dev_oo` optimization work (recent)
- `7e6855df` Optimize likelihood calculations/utilities.
- `1d3915f0` Significant censored-path speedups.
- `074d8f9e` Finite+censored path optimization.

Key optimized surfaces:
- `src/particle_ll.cpp`
  - model adapter resolver and fast-path dispatch (`resolve_race_model_adapter`, raw pointers, cached partitions).
  - all-finite untruncated fast path and SIMD reductions.
  - batched truncation normalizer path.
- `src/utils.h`
  - scalar adapters + raw-buffer adapters (`*_raw`) + batch log-survivor at scalar `t`.
- `src/wald_functions.h`
  - `pnorm_std()` with `USE_FAST_PNORM` path (`fast_norm_phi`).
- `src/model_LBA.h`
  - routes normal CDF calls through `pnorm_std()`.
- `src/Makevars`
  - `-O3 -march=native -ffast-math -fno-math-errno -DUSE_FAST_PNORM`.

### `EMC2_zachdev` (`dev_extended` + local modifications)
Target files inspected:
- `src/model_LBA.h`
  - adds leak kernels: `dleakyba_norm`, `pleakyba_norm`, and vectorized `_c` wrappers.
- `src/model_RDM.h`
  - introduces `dwald/pwald`, `dswtn/pswtn`, and `drdmswtn/prdmswtn` with `sv`.
- `src/gaussian.h`
  - helper math for Gaussian/bivariate CDF/integrals used by new RDM math.
- `src/utils.h` / `src/particle_ll.cpp`
  - adapters for BAwL and new RDM, but based on older likelihood plumbing than `dev-oo`.
- `src/composite_functions.h` / `src/utility_functions.h`
  - small helper additions (`clamp_pos`, `safe_log`, prob clamps), plus Gauss-quad convenience in utility file.

Note: there is **no `gaussian.cpp`** in this tree; implementation is in `gaussian.h`.

## Findings: compatibility map

## 1) `particle_ll.cpp` and `utils.h` from `dev_extended` must **not** replace `dev-oo`
`dev-oo` has a newer/longer fast-path architecture (raw-buffer, SIMD reduction loops, cached trial partitions, batched truncation handling). Replacing with `dev_extended` versions would drop important speed work.

## 2) `model_LBA.h` is mergeable with targeted edits
- BAwL math can be ported in.
- Current `dev_extended` LBA/BAwL calls `R::pnorm`; `dev-oo` LBA intentionally routes via `pnorm_std` (fast pnorm hook).
- Port should preserve the `pnorm_std` strategy for both standard LBA and BAwL where possible.

## 3) `model_RDM.h` needs split-path integration
`dev_extended` RDM is a broad rewrite. In `dev-oo`, RDM currently maps directly to `digt/pigt` and has fast raw adapters expecting 5 columns (`v,B,A,t0,s`).

If we blindly swap:
- we lose tuned fast behavior,
- and risk parameter-layout mismatch once `sv` is present.

Recommended: keep existing `digt/pigt` path as the default fast branch, then add an explicit SWTN branch only when needed.

## 4) `gaussian.h` can be added as a helper module
Add as a standalone header and include only from new SWTN/Wald code paths.

## 5) R model API currently diverges
- `dev-oo` RDM model has no `sv` parameter in `p_types`.
- `dev_extended` RDM includes `sv` always and routes through SWTN wrappers.

A compatibility strategy is required to avoid breaking existing fits/scripts.

## Recommended integration strategy

Guiding decision:
- `RDMSWTN` and `BAwL` are the target canonical supersets.
- Legacy `RDM` (`digt`) and `LBA` are temporary baselines for equivalence/performance gating only.

## Phase 0: Baseline and guardrails
1. Snapshot current `dev-oo` perf and likelihood snapshots.
2. Add focused benchmark scripts before merging (see Benchmark section).
3. Define acceptance gates for eventual replacement:
   - numerical equivalence (within tolerance) when `sv=0` (`RDMSWTN` vs legacy `RDM`) and `k=0` (`BAwL` vs `LBA`),
   - runtime parity or acceptable delta on hot paths.

## Phase 1: Import math kernels only
1. Port BAwL kernels from `model_LBA.h`:
   - `dleakyba_norm`, `pleakyba_norm`, `dleakyba_c`, `pleakyba_c`.
2. Add `gaussian.h` and only the pieces required by RDM SWTN.
3. Port RDM SWTN kernels into `model_RDM.h`:
   - `dwald/pwald`, `dswtn/pswtn`, `drdmswtn/prdmswtn`, and wrappers.
4. Keep compile/link setup from `dev-oo` (`gsl_bundled`, `Makevars` flags).

## Phase 2: Wire into `dev-oo` adapter architecture (do not regress fast path)
1. Add explicit model tags and adapters:
   - `BAwL` / `BAwLIO` -> new BAwL adapters.
   - `RDMSWTN` -> new SWTN/Wald adapters.
2. Add scalar adapters in `utils.h` for BAwL and RDMSWTN.
3. Add raw-buffer adapters as first-class hot paths:
   - `dleakyba_raw` / `pleakyba_raw`.
   - `drdmswtn_raw` / `prdmswtn_raw`.
4. Keep legacy `RDM` (`digt/pigt`) and `LBA` only as temporary reference paths during migration.

## Phase 3: R API and migration path
1. Introduce/keep `RDMSWTN()` as the canonical superset implementation and `BAwL()` as the canonical LBA superset.
2. Ensure parameterization supports exact reduction:
   - `RDMSWTN` with `sv=0` reduces to legacy `RDM` behavior.
   - `BAwL` with `k=0` reduces to `LBA` behavior.
3. Maintain `RDM` and `LBA` wrappers temporarily for compatibility and benchmarking gates.
4. Once gates are met, re-point wrappers (or deprecate legacy internals) so canonical path is `RDMSWTN`/`BAwL`.

## Phase 4: Optimize SWTN path specifically
1. Avoid repeated dynamic quadrature setup in hot loops.
   - pre-cache nodes/weights (20-point default) in static arrays/structures.
2. Fast-reduction special cases:
   - `sv == 0` and `A == 0` -> point Wald (or `digt0`-equivalent).
   - `sv == 0` and `A > 0` -> existing `digt/pigt` path where possible.
3. Route normal CDF calls through a unified wrapper that can use `USE_FAST_PNORM` where numerically acceptable.
4. Add raw-buffer implementations so SWTN can use `dev-oo` all-finite fast path machinery.

## Phase 5: Validation
1. Numerical equivalence tests vs `dev_extended` for BAwL and SWTN cases.
2. Existing likelihood snapshots in `dev-oo` must remain green.
3. Censoring/truncation and GNG behavior regression tests.
4. Stress tests at edge values (`k→0`, `sv→0`, tiny/large `t`, near-boundary CDFs).

## Wald vs `digt`: project stance

- Treat `dwald`/`RDMSWTN` as the canonical math and implementation target.
- `digt` is a temporary reference baseline for migration only.
- For `sv=0`, enforce numerical equivalence and optimize `dwald` hot path to parity with current `digt` performance.
- For `sv>0`, accept additional cost from drift-variability integration, then optimize that path as far as practical.

## Benchmark plan (required before final merge)
1. Compare old RDM (`digt`) vs new SWTN path at `sv=0` and `sv>0`.
2. Compare BAwL (`k=0` vs `k>0`) against existing LBA timings.
3. Scenarios:
   - all-finite no truncation (fast path),
   - truncation/censoring heavy,
   - multi-accumulator race (2, 4+ accumulators).
4. Metrics:
   - wall time per likelihood evaluation,
   - relative slowdown vs current `dev-oo` baseline,
   - numerical delta in log-likelihood.

## Concrete implementation order
1. Add helper math (`gaussian.h`) and BAwL/SWTN kernels in model headers.
2. Add new adapters in `utils.h` compatible with existing `dev-oo` interfaces.
3. Extend `resolve_race_model_adapter` only (do not replace core race likelihood flow).
4. Add R model constructors/wrappers (`BAwL`, and preferably `RDMSWTN`).
5. Add/adjust tests and benchmarks.
6. Only then consider replacing any legacy RDM internals for the `sv==0` case (likely unnecessary).

## Main risk points
- Breaking `dev-oo` raw fast path by importing older `utils.h`/`particle_ll.cpp` patterns.
- Parameter-column layout mismatches (`RDM` 5-col vs SWTN 6-col).
- Numeric stability changes if `pnorm` call paths are altered inconsistently.
- Startup/runtime overhead from quadrature node creation inside hot loops.

## Recommendation
Proceed with a **surgical merge** into current `dev-oo` architecture, but with a canonical-superset migration objective:
- make `RDMSWTN` and `BAwL` the optimized target paths,
- use legacy `RDM`/`LBA` only for equivalence and performance gates,
- then replace legacy internals once parity criteria are met.

## Execution checklist with acceptance gates

## A) Scaffolding and wiring
1. Create branch from `dev-oo` for migration work.
2. Add `src/gaussian.h` from `EMC2_zachdev` (header-only import).
3. Port BAwL kernels into [model_LBA.h](/data/work/EMC2_dev_oo/src/model_LBA.h):
   `dleakyba_norm`, `pleakyba_norm`, `dleakyba_c`, `pleakyba_c`.
4. Port RDMSWTN kernels into [model_RDM.h](/data/work/EMC2_dev_oo/src/model_RDM.h):
   `dwald/pwald`, `dswtn/pswtn`, `drdmswtn/prdmswtn`, vectorized wrappers.
5. Add adapters in [utils.h](/data/work/EMC2_dev_oo/src/utils.h):
   scalar + raw + log-survivor-at-t for `BAwL` and `RDMSWTN`.
6. Extend model resolver in [particle_ll.cpp](/data/work/EMC2_dev_oo/src/particle_ll.cpp):
   `BAwL`/`BAwLIO` and `RDMSWTN` routing into existing fast architecture.
7. Add R-side constructors/wrappers:
   `BAwL()` canonical, `RDMSWTN()` canonical, compatibility wrappers for `LBA()`/`RDM()` retained.

## B) Hot-path optimization tasks
1. Replace direct `R::pnorm` calls on hot paths with unified wrapper compatible with `pnorm_std` / `USE_FAST_PNORM`.
2. Pre-cache quadrature nodes/weights (20-point default) once per process, never per likelihood call.
3. Add degenerate-case fast reductions:
   - `RDMSWTN`: `sv==0,A==0` and `sv==0,A>0` specialized branches.
   - `BAwL`: `k==0` direct reduction to LBA branch.
4. Ensure raw-buffer adapters are used in all-finite untruncated race path.
5. Re-check SIMD reduction loops remain vectorizable (`#pragma omp simd` loops unaffected).

## C) Numerical equivalence tests
Test sets: small deterministic fixtures plus randomized grids.

1. `RDMSWTN(sv=0)` vs legacy `RDM` (`digt/pigt`) at matched parameterization.
2. `BAwL(k=0)` vs legacy `LBA`.
3. Include censoring/truncation and GNG where applicable.

Acceptance thresholds:
1. Per-accumulator PDF/CDF pointwise error:
   - central region (`1e-6 < p < 1-1e-6`): absolute error <= `1e-6`.
   - tail region: absolute error <= `1e-5`.
2. Trial log-likelihood difference: absolute error <= `1e-5` median, <= `1e-4` at 99th percentile.
3. Summed log-likelihood over 10k trials: absolute difference <= `1e-3`.
4. No systematic bias: mean signed LL difference magnitude <= `1e-5`.

## D) Performance benchmark suite
Run each benchmark with warm-up + 20 timed repetitions; report median and p90 runtime.

Scenarios:
1. All-finite, no truncation (primary hot path).
2. Mixed finite + censor/truncation.
3. 2-accumulator and 4-accumulator race setups.
4. `RDMSWTN` with `sv=0` and `sv>0`.
5. `BAwL` with `k=0` and `k>0`.

Acceptance thresholds for replacement readiness:
1. `RDMSWTN(sv=0)` vs legacy `RDM`:
   - all-finite path: <= 5% slower (target parity or faster).
   - censor/trunc mix: <= 10% slower.
2. `BAwL(k=0)` vs legacy `LBA`:
   - all-finite path: <= 5% slower (target parity or faster).
   - censor/trunc mix: <= 10% slower.
3. No memory-regression red flags (no repeated heap growth in long benchmark loops).

## E) Cutover criteria
1. All numerical gates in section C pass.
2. All performance gates in section D pass.
3. Existing `dev-oo` likelihood regression tests pass unchanged.
4. Then re-point wrappers so canonical internals are `RDMSWTN`/`BAwL`.
5. Keep legacy kernels only behind a temporary fallback flag for one release cycle; remove after stability window.

## F) Suggested implementation sequence (commit-level)
1. Commit 1: add math headers/kernels (`gaussian.h`, `model_LBA.h`, `model_RDM.h` additions).
2. Commit 2: adapter wiring in `utils.h` + `particle_ll.cpp` resolver integration.
3. Commit 3: raw fast-path optimizations and quadrature caching.
4. Commit 4: R API additions (`RDMSWTN`, `BAwL`) + compatibility wrappers.
5. Commit 5: equivalence tests.
6. Commit 6: benchmark harness + benchmark report artifact.
7. Commit 7: cutover (wrapper re-point / deprecation notes) once gates pass.
