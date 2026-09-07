# EMC2 optimisation sweep — 2026-09-06

Status: **implemented 2026-09-07**, everything except Commit 3 (see below).
Written against branch `playground` (`6dbc1f5b`, plus the uncommitted
FIFO-framing work in `R/chain_pool.R` and the `TrendEngine` changes). Every
number below was measured on this worktree, on this box (`legacy` build
profile: `-O3 -march=native -ffast-math -fno-finite-math-only -fno-math-errno
-DUSE_FAST_PNORM`). Single-process timings are **CPU time**, not elapsed.

**Correction to an earlier premise in this document.** The first pass recorded
this box as having "1 usable core", on the strength of `nproc`. That is wrong:
`nproc` returns 1 only because the harness exports `OMP_NUM_THREADS=1`. The
affinity mask is `0-31`, `nproc --all` is 32, and `parallel::detectCores()` --
which is what actually governs the worker pool -- is **32**. Every conclusion
that turned on the core count has been re-measured; the one that changed is the
second bullet of Commit 10, which is now done and is the largest single win in
the whole sweep on a wide model.

## Outcome

| benchmark | before | after | delta |
|---|---|---|---|
| LBA likelihood, 100 particles, 4 000 rows | 54.0 ms | 38.0 ms | **−29.7 %** |
| LBA likelihood, censored/mixed path, 3 412 rows | 46.9 ms | 32.5 ms | **−30.6 %** |
| end-to-end preburn fit, 6 subjects x 400 trials, 60 iters | 8.43 s | 5.92 s | **−29.8 %** |
| per-particle prologue's share of a likelihood call | 23.4 % | 3.7 % | 6.3x cheaper |
| `.cache_ll_data_attrs` rebuild, 20 000 trials | 290 ms | 1.5 ms | **194x**, and now linear |

And, added 2026-09-07 once the box's real core count was known (32 subjects x
200 trials, P = 56, 40 preburn iterations, **elapsed** per sampler iteration):

| workers | before | after | delta |
|---|---|---|---|
| 1 | 1066 ms | 1049 ms | neutral (see Commit 10) |
| 2 | 555 ms | 557 ms | neutral |
| 4 | 332 ms | 275 ms | **-17 %** |
| 8 | 176 ms | 152 ms | **-14 %** |
| 16 | 147 ms | 84 ms | **-43 %** |
| 32 | 217 ms | 87 ms | **-60 %** |

Note the shape of the *before* column: past 16 workers the old build gets
**slower** (147 -> 217 ms), because the per-iteration broadcast grows with the
worker count. The new build is monotone and flattens where the work does.

Correctness, at every step:

* Log-likelihoods **bit-identical** on both the all-finite and the mixed path.
* A 2-chain, 4-subject, 25-iteration hierarchical preburn fit is
  `identical()` between the baseline and the optimised build — every `alpha`,
  `theta_mu`, `theta_var` and `subj_ll`.
* Full suite at `EMC2_TEST_LEVEL=full NOT_CRAN=true`: the untouched baseline
  gives **failed=97, error=16**; the final build gives **failed=94, error=16**,
  a strict subset — **the sweep adds no failure anywhere.** Those 94/16 are
  pre-existing on this branch. Passing assertions 5 884 -> 5 927.
  (The three-assertion difference is two `test-worker-pool.R` core-reallocation
  tests that are flaky *inside the full suite* on both builds and pass on both
  when the file is run alone; they happened to pass in the final run.)
* The same suite under the new assertion build
  (`EMC2_EXTRA_CPPFLAGS=-DEMC2_PT_CONST_CHECK`), which re-scans every column
  whose row-constant or cell-constant claim is about to be acted on: same
  failure set, **zero** flag violations.
* The vectorised `.cache_ll_data_attrs` was checked attribute-for-attribute
  with `identical()` against a transcription of the loop it replaces, on seven
  data shapes (all-finite / 15 % / 50 % omissions, RACE, censored RACE, mixed
  2AFC-GNG with and without omissions).
* The incremental `accept_rate()` was checked against the original over every
  iteration of a 320-iteration run at 1, 3 and 40 subjects, warm and cold cache.
* The Commit 10 broadcast change (added 2026-09-07): a 2-chain, 8-subject,
  30-iteration preburn fit is `identical()` to the pre-change build at 1, 3 and
  8 workers, on both the `fork` and the `spawn` pool backend, and the two
  backends agree with each other.

## What is NOT done, and why

* **Commit 3** (adapt the particle count during burn) — dropped at the user's
  direction: adaptation deliberately happens *after* burn, before sampling, so
  the gate is intentional and the premise of that item was wrong.
* **Commit 13, last bullet** (geometric growth in `extend_obj`) — conditioned in
  the plan on a long fit actually being seen to thrash. Not observed. Code
  unchanged.
* **Commit 13, `.emc_cores_now()` caching** — implemented, then reverted: it
  changes an observable contract (a released core must be visible to the next
  caller) and made three worker-pool tests fail intermittently, for a saving of
  tens of microseconds per ~100 ms iteration. Details under Commit 13.
* **Commit 2's "free bonus"** (skip the kernel when a particle is entirely out
  of bounds) — now measured, and it does not pay: at proposal sd 0.1/0.5/1.0/2.0
  the fraction of particles with *no* in-bounds row is 0 % / 0 % / 0.2 % / 6.0 %.
  Real preburn clouds sit at the low end. Recorded so it is not re-attempted.
* **Packing the group covariance as an upper triangle** on the wire — measured
  (another 2.2x on the send) and **rejected**: `theta_var` is not exactly
  symmetric, so it would silently move every draw. See Commit 10.
* **Broadcasting the shared blob through a file rather than `n_workers` pipes**
  — implemented after the isolated measurement (4.5 ms against 41 ms for the
  covariance pipes on 32 workers). The master publishes one complete payload
  file per iteration and sends only its path; workers report missing or
  truncated files through the existing serial-recompute path. The path is
  removed after all replies, including failed-worker fallbacks, are settled.

---

## 0. Where the time actually goes

### 0.1 R vs C++

`Rprof` on a hierarchical LBA preburn (8 subjects x 200 trials, 7 sampled
parameters, 150 iterations, 1 core, worker pool active):

| bucket | share |
|---|---|
| `calc_ll_oo` (C++) | **92.4 %** |
| `new_particle` (R, self) | 0.9 % |
| everything else in R (draws, chol, dmvnorm, log-sum-exp, gibbs, fill) | ~6.7 % |

`gibbs_step_standard` was 0.8 % — consistent with the earlier finding recorded in
`project_emc2_worker_pool.md` that it is not worth optimising.

**Conclusion: R-side micro-optimisation of the sampler loop cannot buy more than
a few percent. The wins are inside `calc_ll_oo`.** The R-side items in this plan
are therefore filed under fragility/correctness or under paths that are not the
per-iteration hot loop (setup, `predict`, IC).

### 0.2 Inside `calc_ll_oo` — the mapping tax

Harness: single-subject LBA, `v ~ lM`, `B ~ E`, `A ~ 1`, `t0 ~ 1`,
`constants = c(sv = log(1))`, 2 000 unique trials (4 000 dadm rows),
100 particles, realistic parameter values (so the natural-space kernel path is
taken, not the log fallback). One `calc_ll_manager()` call = **50.7 ms**.

Phase timers compiled into `PtMapper::prepare()` (`src/particle_ll.cpp:356`):

| phase | per call | share |
|---|---|---|
| `fill_from_particle_row_planned` | 0.55 ms | 1.1 % |
| `map_from_designs` + `c_do_transform_pt` | 9.0 ms | **18.1 %** |
| `c_do_bound_pt_from` | 2.6 ms | **5.3 %** |
| model kernel + summation + `make_pt_mapper` | 37.8 ms | 75.6 % |

Cross-check: forcing `t0` so large that every row short-circuits inside
`dbawl_raw`/`pbawl_raw` leaves 16.6 ms of the 50.7 ms — the same ~33 % once the
kernel's own per-row early-exit loop is included.

I also confirmed the kernel is *not* falling back to the expensive log-space
path: instrumented branch counters gave 199 968 / 200 000 pdf rows and
200 000 / 200 000 cdf rows on the cheap natural path. The 75.6 % is genuine
kernel arithmetic.

### 0.3 Why the mapping tax is pure waste

The parameter columns only carry **design-cell** information, but every stage of
the per-particle prologue runs at **trial** resolution. For the benchmark design:

```
v              nrow=2   ncol=2   expand: 4000 rows -> 2 distinct
B              nrow=3   ncol=3   expand: 4000 rows -> 3 distinct
A, t0, eta, pContaminant, pGuess, sv   nrow=1      -> 1 distinct
```

Six of the eight parameters are intercept-only, yet `c_do_transform_pt`
(`src/transform_utils.cpp:156`) still calls `std::exp` 4 000 times per column per
particle, and `c_do_bound_pt_from` (`src/transform_utils.cpp:122`) still compares
4 000 rows per bound spec per particle. `objdump` on `transform_utils.o` shows a
plain undefined `exp` — no `_ZGV*` vector variant — so those are 4 000 scalar
libm calls, not a vectorised loop.

This is the single largest recoverable cost in the package, and it is
**model-independent**: every race model, the DDM, LogicalRules, SS and the
censored/truncated path all go through the same prologue.

### 0.4 Fixed per-call cost (checked, and it is fine)

`calc_ll_manager` cost vs particle count, same harness:

| particles | 1 | 2 | 5 | 10 | 25 | 100 |
|---|---|---|---|---|---|---|
| ms | 1.1 | 1.7 | 3.2 | 5.7 | 13.2 | 50.7 |

Linear fit: 0.60 ms fixed + 0.501 ms/particle. The fixed part (`make_pt_mapper`:
pretransform, `ParamTable` construction, design-plan build) is 1.2 % at 100
particles and ~5 % at a 4-way particle split. **`make_pt_mapper` caching is not
worth doing** — I had expected it to be a target and it is not.

### 0.5 The censored/truncated path

Same shape of harness with 15 % omissions and `UC = 2.0` (6 812 dadm rows,
`emc2_all_finite_trials = FALSE`, so the raw fast path is off and
`c_log_likelihood_race` runs): 95.3 ms per 100-particle call, i.e. 14 ns per
row-particle vs 12.7 ns on the all-finite path. `.is_valid_ll_cache` was
**0.6 %** of that call — the fix recorded in `project_emc2_ll_cache_quadratic.md`
holds; there is nothing left there.

---

## Commit 1 — Row-constant parameter columns (SIGNIFICANCE: HIGH)

**DONE.** `src/ParamTable.h` (`col_const`, `const_init`, `col_is_const`,
`DesignEntry::single_cell`), `src/transform_utils.cpp`, `src/particle_ll.cpp`.
Shipped with the assertion build described under "Risk" below, as
`EMC2_EXTRA_CPPFLAGS=-DEMC2_PT_CONST_CHECK` (a `./configure` pass-through added
for it), and with `tests/testthat/test-param-table-prologue.R`. Superseded in
part by Commit 2, which generalises it from one cell to n.

### What

Track, per base column of `ParamTable`, whether the column currently holds one
repeated value, and take an O(1) path in the three places that otherwise scan
`n_trials`:

* `map_from_designs` (`src/ParamTable.h:476`) — when the design entry reads a
  single design row for every trial *and* every coefficient column is itself
  row-constant, evaluate the linear combination once and `std::fill`.
* `c_do_transform_pt` (`src/transform_utils.cpp:156`) — one `exp`/`pnorm`, then
  `std::fill`.
* `c_do_bound_pt` / `c_do_bound_pt_from` (`src/transform_utils.cpp:89,122`) —
  one comparison; on failure write `false` across, on success skip the spec
  entirely.

The flags are maintained, not sniffed: `reset_base_to_zero` marks everything
constant, `fill_from_particle_row_planned` marks each column it scalar-fills,
and `map_from_designs` recomputes the flag for each output it writes. Invariant
columns are touched by none of those and correctly keep the flag they got on the
template particle.

`DesignEntry` gains `bool single_cell`, computed once in `init_design_plan`:
true when `expand_idx` is constant, or (uncompressed) when the design has one
row.

### Guard

Enable only when there is no `TrendRuntime` — the trend engine writes into `base`
outside this pipeline (`TrendRuntime::apply_base_for_op`). That is exactly the
condition the existing "planned" lane already requires
(`src/particle_ll.cpp:356`), so no new capability matrix. `set_column_by_name`
clears the flag defensively.

### Measured

| benchmark | before | after | delta |
|---|---|---|---|
| LBA likelihood, 100 particles, 4 000 rows | 50.7 ms | 43.7 ms | **−13.8 %** |
| LBA likelihood, censored/mixed path, 6 812 rows | 95.3 ms | 82.7 ms | **−13.2 %** |
| end-to-end preburn fit, 6 subjects x 400 trials, `v ~ lM*E`, `B ~ E*lR`, 60 iters | 7.98 s | 6.62 s | **−17.1 %** |

Log-likelihoods are **bit-identical** in every case
(`-572.7993 / -364.4686 / -752.0166` and `-5586.640 / -5893.988`).
`test-likelihoods.R`, `test-oo-particle-matrix.R` and `test-trend.R` pass; the
one failure in `test-map.R` (`.credint_label_mar` not found) reproduces on the
unpatched build and is unrelated.

### Risk

Low, and it fails loudly rather than silently: if a flag were ever wrongly set
the column would be filled from row 0, which changes the likelihood by a visible
amount rather than a rounding-level amount. Ship it with an assertion build
(`EMC2_PT_CONST_CHECK`) that re-scans the column and stops on disagreement, plus
a test that fits the same data with and without a dummy trend and compares
likelihoods.

### Ceiling

This recovers the intercept-only parameters. The `map + transform + bound`
budget is 23.4 % of a likelihood call, so ~10 points remain — that is Commit 2.

---

## Commit 2 — Design-cell resolution for the rest (SIGNIFICANCE: HIGH, larger effort)

**DONE, and it was the single largest item: it took the prologue from 10.0 % of
a likelihood call to 3.7 %.** Implemented not as the joint-cell-index scratch
table sketched below but as the natural generalisation of Commit 1 — a
per-column *cell* tag (`ParamTable::col_cell`, `DesignEntry::cell_usable /
n_cells / cell_rep`) saying "this column's distinct values are laid out by
design i's expand map". `map_from_designs` then evaluates n_cells values and
scatters, `c_do_transform_pt` takes n_cells transcendentals and scatters, and
the bound check settles the verdict in n_cells comparisons — skipping the spec
outright when they all pass. No scratch table, no joint index, and
row-constancy stays the n_cells == 1 corner of the same mechanism.

Two things worth knowing:

* The cell accumulation deliberately mirrors the general per-row loop's operand
  order (coefficient-major into a zeroed accumulator). The first version folded
  per cell instead and differed by one ulp on a continuous-covariate design,
  because the compiler contracted the two loops into FMAs differently. With the
  orders matched the two routes are bit-identical, which the new tests assert at
  `tolerance = 0`.
* Designs with more than `EMC2_PT_MAX_CELLS` (256) cells take the per-row path:
  that is the cap on the stack scratch, and by then the gather has little left
  to save.

The kernel-side prize below (hoisting each model's parameter-only
precomputation to cell resolution) is now *unlocked* but not taken — every
model would need its own change. The cell index it needs is on the table.

**Original sketch, kept for the record:**

Commit 1 is the `n_cells == 1` special case. The general statement is: a
parameter column takes at most `nrow(design)` distinct values, and `nrow(design)`
is typically 2–12 against `n_trials` in the thousands.

### Shape

1. Once per likelihood call, build the **joint cell index**: the interaction of
   every design's `expand` vector, giving `cell_of_row[0..T)` and `n_cells`.
   Fall back to the current path when `n_cells > n_trials / 4` (continuous
   covariates in a design blow the cell count up; that guard is what keeps this
   honest).
2. Per particle, run fill → map → transform → bounds on a scratch `ParamTable`
   with `n_cells` rows.
3. Scatter each `keep_names` column back to full length once
   (`base[r] = small[cell_of_row[r]]`), and expand the bound verdict.

Per parameter per particle this replaces ~6 full-length passes (two fills, K map
passes, one transcendental, one comparison) with one gather. Columns that are
row-constant stay on Commit 1's `std::fill`, which is faster than a gather.

### The bigger prize behind it

Once a cell index exists, **every kernel can hoist its parameter-only
precomputation to cell resolution.** Concretely, `ba_natural_pdf` and
`ba_natural_cdf` (`src/model_LBA.h:212,281`) both call `natural_normalizer`,
i.e. `pnorm_std(v/sv)`, once per row, and it depends on no time. That is one of
three `pnorm` calls in the LBA hot path. For the heavier families the payoff is
much larger — this is the same lever that gave 4x on BAwDD by memoising the
saturation Newton solve on `(c, rho, gamma)`
(`project_emc2_bawdd_alpha_unbounded.md`), generalised so every model gets it for
free instead of one at a time.

### Free bonus available at the same time

With bounds evaluated at cell resolution, an entirely out-of-bounds particle is
detected in O(n_params) instead of O(n_trials), and `lls[i]` can be set to
`n_expand * min_ll` **without running the kernel at all**. In preburn the
proposal cloud is wide and a large fraction of particles are out of bounds, so
this is worth measuring separately — instrument the out-of-bounds fraction per
stage before committing to a number.

---

## Commit 3 — Adapt the particle count during burn, not only sample (SIGNIFICANCE: HIGH, needs validation)

**NOT DONE — the premise is wrong.** Particle-count adaptation is deliberately
placed *after* burn and before sampling; the `length(pm_settings$mix) > 3` gate
is the design, not an oversight. Left exactly as it is.

**Original item, kept for the record:**

`update_pm_settings` (`R/sampling.R:1346`) only adapts `n_particles` when
`length(pm_settings$mix) > 3` (i.e. the `sample` stage) **and** `gd_good`.
Everywhere else the count is pinned at `round(particle_factor * sqrt(n_pars))`
(`R/fitting.R:236`), doubled in preburn — 224 particles for a 20-parameter model,
regardless of how well the proposal is doing.

Burn is routinely the longest and most expensive stage of a fit, and cost is
almost exactly linear in the particle count (§0.4). The ESS logic already in the
function is stage-agnostic; only the gate is not.

The count does not change the invariant distribution of this importance-MH
kernel — it changes efficiency only — so adapting it in burn/adapt is legitimate.
Proposed: same ESS rule, a higher floor in burn than in sample (say 50 rather
than 25) and no upward clamp relaxation, gated behind an option for one release
so it can be A/B'd.

**Validate before shipping**: fit the same data with a fixed and an adapted burn
particle count, and compare (a) wall-clock to a given `max_gd`, (b) ESS per CPU
second in the subsequent sample stage, (c) the posteriors themselves. This is the
one item in the plan whose size I have not measured; it is here because the
mechanism is already written and switched off, not because I have a number.

---

## Commit 4 — Pointer hoists in the per-particle prologue (SIGNIFICANCE: MEDIUM)

**DONE**, including both "not separately measured" items: `c_do_bound_pt_from`
now folds into a `PtMapper`-owned buffer (`c_do_bound_pt_from_into`) instead of
allocating and cloning per particle, and `build_plan` no longer puts a column in
both `plan_zero_base_idx` and `plan_fill_pairs`.

**Measured: −2.0 %, standalone.**

Three inner loops index Rcpp matrices with `operator()`, which performs a bounds
check and an `i + nrow*j` multiply per element:

* `map_from_designs` (`src/ParamTable.h:534,556`): `design(drow, j)` in both the
  split-transform and ordinary accumulation loops. Hoist
  `const double* dcol = &design(0, j);` and split the loop on
  `expand_idx.empty()` so the emptiness test leaves the inner loop too.
* `c_do_bound_pt` (`src/transform_utils.cpp:108`) and `c_do_bound_pt_from`
  (`src/transform_utils.cpp:146`): `base(i, col_idx)`. Hoist
  `const double* col = &base(0, col_idx);`.

Measured 50.7 ms → 49.8 ms on the LBA benchmark, identical likelihoods. It is
subsumed by Commit 1/2 for the constant columns but still applies to the
multi-cell ones, so land it either before or alongside them.

Two more allocations in the same region, not separately measured:

* `c_do_bound_pt` allocates an R `LogicalVector(nrows)` per particle and
  `c_do_bound_pt_from` additionally `Rcpp::clone`s the seed. Both could write
  into a reusable buffer owned by `PtMapper`.
* `build_plan` (`src/particle_ll.cpp:296`) puts every non-invariant base column
  into `plan_zero_base_idx` *and* every particle-backed one into
  `plan_fill_pairs`, so `fill_from_particle_row_planned`
  (`src/ParamTable.h:682`) memsets each sampled parameter's column to zero and
  then immediately overwrites all of it. Removing the overlap deletes ~7
  full-length memsets per particle. The whole fill phase measured 1.1 %, so this
  is worth a few tenths of a percent — take it because it is three lines, not
  because it matters.

---

## Commit 5 — Standard-normal `dnormP` fast path (SIGNIFICANCE: LOW–MEDIUM)

**DONE** (`src/wald_functions.h`). Re-measured in place at **−6.6 % on the LBA
likelihood** (44.1 -> 41.2 ms), well above the 1–2 % estimated from the
micro-benchmark below.

**Original estimate: −1 to −2 % on LBA; larger on the log-space families.**

`dnormP` (`src/wald_functions.h:61`) forwards unconditionally to `R::dnorm`, an
out-of-line call into Rmath that re-validates `mean`/`sd` on every invocation.
Micro-benchmark on this box: `R::dnorm` 10.5 ns, inline
`M_1_SQRT_2PI * exp(-0.5*x*x)` 4.9 ns, `R::pnorm` 31.3 ns.

There are 70 call sites and all but a handful pass `(0, 1)`. Rmath's own `x < 5`
branch *is* the inline expression, so:

```cpp
if (mean == 0.0 && sd == 1.0) {
  const double z2 = x * x;
  if (log) return -0.5 * z2 - LOG_SQRT_2PI;
  if (z2 < 25.0) return 0.39894228040143267794 * std::exp(-0.5 * z2);
}
return R::dnorm(x, mean, sd, log);   // |x| >= 5, NaN, non-standard
```

is exact where it fires and defers to Rmath's split-argument branch where that
branch actually matters. On the LBA benchmark this was worth 1–2 % (0.0507 →
0.0505–0.0499 s across runs); it is listed low because the LBA natural path calls
it only twice per row. It should be worth more in `bawl_corr_exact.h` (13 call
sites), `bawl_geometry.h` (10) and the `log = true` sites in `model_BAwD/F/R.cpp`
and `model_BTAwL.cpp`, which is where the family's log-space work lives. Measure
on a BAwL/BAwD fit before claiming a number there.

---

## Commit 6 — Mixed (censored/truncated) path per-particle overheads (SIGNIFICANCE: MEDIUM)

**DONE except the full column-pointer contract.** The `pC`/`pG` scans now hoist
their column pointer; `ll_unique` moved into `RaceSharedState`; and
`materialize_reusable()` re-copies only the columns that can change per particle
(`materialize_into_subset`), which is where most of the memcpy went — the
invariant columns hold the same natural-scale values for the whole call. Giving
`c_log_likelihood_race` a `const double* const*` contract instead of a
materialised matrix is left undone: it reaches into `correlated_likelihood.cpp`'s
quadrature paths as well, and the remaining copy is now the smaller half.

This path is the one the current research programme actually runs
(`project_emc2_cens_trunc2_divergence.md`, the omission work), and it is the one
that does *not* get the raw column-pointer treatment.

* **`pC` / `pG` "all zero" scans** (`src/particle_ll.cpp:1572,1605`). Both loop
  `for (i = 0; i < pars.nrow(); ++i)` with the bounds-checked `pars(i, col)`,
  and both break early only when a non-zero is found — so the *common* case
  (`pContaminant` pinned at 0) is the one that scans all `n_trials` rows, twice,
  per particle. `pContaminant` and `pGuess` are intercept-only in every design I
  have seen: with Commit 1's flag in place this becomes `pars(0, col) == 0`;
  without it, at minimum hoist the column pointer.
* **`materialize_reusable()` per particle** (`src/particle_ll.cpp:1155`). The
  mixed path takes a materialised `NumericMatrix pars`, so every particle pays a
  `keep_names.size() * n_trials * 8` byte memcpy (436 kB in the censored
  benchmark) that the all-finite path avoids entirely by handing the kernels
  `const double* const*` into `ParamTable::base`. Giving
  `c_log_likelihood_race` the same column-pointer contract as
  `c_log_likelihood_race`'s raw twin is the structural fix; it also removes the
  `RACE_mask` NA-poking loop (`src/particle_ll.cpp:1537`), which only exists
  because the copy is private.
* `std::vector<double> ll_unique(n_unique_trials, min_ll)` is heap-allocated per
  particle; it belongs in `RaceSharedState` next to `res_buf`.

---

## Commit 7 — `.cache_ll_data_attrs` build path is quadratic (SIGNIFICANCE: MEDIUM, setup-only)

**DONE.** Both loops vectorised (the RACE row loop and the per-unique-trial
loop). 20 000 trials: 290 ms -> 1.5 ms, and the growth is now linear. Verified
`identical()` attribute-for-attribute against a transcription of the old loop on
seven data shapes including RACE and mixed 2AFC/GNG.

`R/utils.R:382-404`. The per-unique-trial loop grows two index vectors with
`c(finite_rt_unique_trial_indices, j0)` / `c(other_unique_trial_indices, j0)`,
which reallocates and copies on every iteration.

Measured, `force_rebuild = TRUE`, 15 % omissions:

| trials | dadm rows | rebuild |
|---|---|---|
| 2 000 | 3 412 | 10.6 ms |
| 8 000 | 13 612 | 64.8 ms |
| 20 000 | 34 012 | **256.8 ms** |

4x the rows costs 6.1x the time, then 2.5x costs 4.0x. It is paid once per
subject at fit setup (`R/fitting.R:913`, `R/design.R:1535`) and again on every
`predict()` / `make_data()` — so a 100-subject censored study spends ~26 s of
pure setup here, and posterior-predictive loops pay it repeatedly.

The whole loop is vectorisable exactly the way `.is_valid_ll_cache` already was
(`R/utils.R:266-294` shows the pattern: `trial_of_row`, `acc_of_row`,
`matrix(..., nrow = n_lR)` + `colSums`). Also vectorise the scalar RACE loop at
`R/utils.R:339`, which is `for (i in seq_len(n_trials))` over a factor.

Not a sampling win — it never runs inside the iteration loop — but it is the
largest single R-side cost I found anywhere.

---

## Commit 8 — `sub_blocking()` averages one matrix N times (SIGNIFICANCE: MEDIUM — correctness)

**DONE.** Both indices fixed, the normaliser is now the explicit count rather
than `out[1,1]` (identical value, but immune to a NaN on the diagonal), and the
`cutree(k = n_blocks)` case the existing comment warned about is guarded:
a singleton group becomes its own block and `k` is clamped to the group size.

`R/fitting.R:547`:

```r
for(i in 1:length(covs)){
  cov_tmp <- covs[[1]]                  # should be covs[[i]]
  for(j in 1:length(cov_tmp)){
    out <- out + cov2cor(cov_tmp[[1]])  # should be cov_tmp[[j]]
  }
}
```

Both indices are stuck at 1, so `out` is `n_chains * n_subjects` copies of
**chain 1, subject 1**'s correlation matrix. The `hclust`/`cutree` that decides
which parameters share a proposal block therefore sees a single subject from a
single chain rather than the pooled structure.

This runs whenever `n_blocks > 1` — i.e. exactly the high-dimensional fits where
blocking is supposed to rescue the acceptance rate. It is a sampling-efficiency
bug dressed as a typo. Fix, then re-check that the resulting blocks change on a
model where subjects differ (they should).

While there: the comment on `R/fitting.R:563` already flags that `cutree(k =
n_blocks)` "could go wrong if one group has just one member" — worth a guard now
that the input will actually vary.

---

## Commit 9 — `extend_obj()` decides what to extend by numeric sniffing (SIGNIFICANCE: MEDIUM — fragility)

**DONE.** Replaced with the structural rule (`is_iteration_array()`: the final
dimension equals the iteration count before the extend). I checked every
`sample_store_*` — standard, factor, infnt_factor, SEM, diag_gamma and the base
alpha/subj_ll pair — and every element has `iters` as its last dimension, so the
guard was protecting nothing. `reject_sample_iteration()` carried the same
`d[1] != d[2]` guess and is fixed with it.

`R/sampling.R:1435`, inside `extend_obj()`):

```r
if(nrow(obj) == ncol(obj)){
  if(nrow(obj) > 1){
    if(mean(abs(abs(rowSums(obj/max(obj))) - abs(colSums(obj/max(obj))))) < .01) return(obj)
  }
}
```

`extend_sampler` `rapply`s this over the whole `samples` list every time a stage
block starts. The guard is trying to say "this is a fixed covariance-shaped
matrix, don't extend it", but it decides that by looking at the *numbers*:

* It only even looks when `nrow == ncol`, i.e. when the parameter count happens
  to equal the current iteration count — a condition every run passes through
  exactly once.
* `max(obj)` is not guarded. An all-zero or all-NA slice gives `NaN`, and
  `if (NaN < .01)` is an error, not a `FALSE` — the fit dies with "missing value
  where TRUE/FALSE needed" at a random iteration count.
* If it ever returns `TRUE` for a genuine sample array, that array silently stops
  growing and the next `fill_samples` writes out of range.

This is the same class of failure as the SIGPIPE race: a heuristic that is right
almost always, and whose failure mode is not "slower" but "dead run at scale".
Replace it with a structural rule — extend an array iff its last dimension equals
the current stored iteration count (`length(samples$stage)` before the update),
which is exact, cheap and cannot be fooled. I checked every `sample_store_*`
(`variant_standard`, `factor`, `infnt_factor`, `SEM`, `diag_gamma`) and found no
2-D square element that needs protecting, so the guard may simply be dead weight
— confirm before deleting rather than replacing.

---

## Commit 10 — Worker pool: a wedged worker hangs the fit forever (SIGNIFICANCE: MEDIUM — fragility)

**DONE (first and third bullets).** Rather than rewriting the transport to
non-blocking reads — the code that had just been fixed for a race — the
iteration path now reuses the completion FIFO the dynamic likelihood queue
already owns: `.emc_wpool_iter` sets `notify`/`w` on each message, and
`.emc_wpool_await_done()` waits on that non-blocking channel with the same
spin-then-back-off schedule, polling worker liveness and a deadline. On timeout
it names the worker and the subjects it holds, kills it, retires the pool and
recomputes those shares in the master — the path that raises the underlying
error properly. The deadline is 20x the median observed receive time, floored at
120 s (600 s cold), overridable with `options(emc2.worker_timeout=)`; `Inf`
disables it. Replies are now also consumed in completion order rather than
worker-index order. The send-ordering skew is documented in place.

`.emc_wpool_get_bytes` (`R/chain_pool.R:92`) does a blocking `readBin` on the
answer FIFO with no deadline, and `.emc_wpool_iter` (`R/chain_pool.R:1151`) waits
on each worker in turn. The recovery paths are good — a *dead* worker closes its
pipe, the read returns 0 bytes, the master recomputes that share serially — but
they all key on the worker **exiting**.

A worker that is alive and stuck produces no such signal. That is not
hypothetical here: `project_emc2_rlf_degenerate_stall.md` records a silent
sampler wedge from a single PDE solve diverging at `alpha -> 2`, and the FPE/RLF
solvers now sit on this path. The result is a fit that appears to be running with
a frozen progress bar and no way to tell which subject is responsible.

Proposal: give the receive loop a liveness poll — non-blocking reads with a short
sleep, checking `.emc_wpool_pid_alive()` and a per-iteration deadline derived
from a rolling median of `it$times`. On timeout, warn naming the worker and the
subjects it holds, terminate it, and fall back to the master's own
`.emc_wpool_compute()` for that share (which will raise the underlying error
properly, exactly as the `failed` branch already does).

Two smaller items in the same file, both now settled:

### The shared broadcast — **DONE 2026-09-07, and it was the big one**

The original note read: "`.emc_wpool_iter` embeds the full serialised `shared`
blob in every worker's message, so `group_chol` crosses `n_workers` pipes per
iteration... only worth doing if profiling shows `shared_serialize`/`send` above
a percent." Once the box's real core count was known, that profile could be
taken. On 32 workers with a 56-parameter model, `request_send` was **72.6 % of
the entire sampler iteration.**

Two facts the original note had wrong, both in the same direction:

* `.chol_factor()` keeps **the root and its inverse** next to the covariance, so
  the cache is *three* P x P matrices, not two. Sending the covariance alone is
  a 2.9x cut (77 162 -> 26 750 bytes at P = 56), not the 2x predicted.
* A blocking R FIFO moves about **16-27 MB/s**, not the GB/s a pipe suggests.
  Cost is linear in payload over the whole range that matters (1 kB -> 2.0 ms,
  25 kB -> 34 ms, 77 kB -> 153 ms for 32 workers), so the cut converts almost
  1:1 into wall time.

Implemented as: `.emc_wpool_iter` puts `group_var` + `idx_list` on the wire and
keeps its own built cache for the master-side recompute path;
`.emc_wpool_compute_particle` rebuilds the cache with `build_group_chol_cache()`
when only the covariance arrived. One `chol()` per component per worker per
iteration, against three P x P matrices per worker per iteration on the wire.

The trade only gets better with P: wire is O(P²) at ~2e7 B/s, `chol` is
O(P³)/3 at ~2e9 flop/s, and they do not cross until P is in the **thousands**.
There is no configuration in which this loses.

Measured, 32 subjects x 200 trials, P = 56, per iteration:

| | before | after |
|---|---|---|
| `shared_bytes` per worker | 77 162 | 26 750 |
| wire bytes per iteration | 2.48 MB | 0.94 MB |
| `request_send` | 136.6 ms | 39.9 ms |
| iteration total | 201 ms | 97 ms |

`response_wait` is now 55 % of the iteration and equals `worker_max` (58 ms) to
within a few ms — i.e. the pool is finally bound by the work rather than by the
broadcast.

A second, smaller item fell out of it: when the pool is dead (`cores_per_chain
= 1`, Windows, no `mkfifo`) the master used to serialise the blob every
iteration and throw it away, since it computes every share itself from the
decoded object. It now builds the wire copy only when there is a live worker.

**Correctness.** A 2-chain, 8-subject, 30-iteration preburn fit is `identical()`
to the pre-change build — `alpha`, `theta_mu`, `theta_var`, `subj_ll` — at 1, 3
and 8 workers, on **both** pool backends (`fork` and `spawn`), and fork and
spawn agree with each other. This is the property that matters: every worker
starts from the identical matrix and makes the identical LAPACK call, so its
factors are bit-for-bit the ones the master would have sent. Two
`test-worker-pool.R` tests pinned the old wire shape and were rewritten to pin
the new one; a third was added asserting that the wire payload is less than half
the cache's size *and* that the worker's reconstruction is `identical()` to the
master's object. `test-worker-pool.R` is 175/175 green, and the full suite is
unchanged at failed=94 / error=16 over the same 49 tests.

### Measured and rejected here too

* **Packing the covariance as an upper triangle** (25 kB -> 13 kB, and the pipe
  round trip 41 ms -> 19 ms on 32 workers). **Do not do this.** `theta_var` is
  *not* exactly symmetric — 0 of 16 stored iterations satisfy
  `identical(V, t(V))` — so rebuilding from a triangle silently changes the
  matrix, hence the factor, hence every draw. A guard on exact symmetry would
  simply never fire.

### File broadcast — **DONE 2026-09-07**

The shared covariance/index payload is written once to an atomic temporary-file
rename in the pool directory. Each worker request carries only that path, and
the master unlinks it after all replies or serial recomputations have finished.
Workers validate the file size and complete read before unserialising; a missing
or truncated file is a worker-local failure, so the existing master-side
recompute and pool-degradation logic remains the safety net. No file is created
for a dead or single-worker pool. The focused worker-pool suite covers the
round-trip, missing-file recovery, cache reconstruction, reallocation, and both
fork and spawn fit paths.

### Still on the table

* Sends are issued to all workers before any receive, and a send blocks until
  that worker drains the pipe. Worker `n` therefore starts after workers
  `1..n-1` have each taken their message, which slightly violates the
  simultaneous-start assumption behind the LPT partition. This mattered much
  more before the change above (the last of 32 workers started ~170 ms in, and
  `worker_max` was only 66 ms — the pool was starved, not busy). At 26 kB it is
  down to ~40 ms and the file broadcast above would remove it outright.
* `.EMC_WPOOL_WRITE_CHUNK` is 4096 (PIPE_BUF) so that a short write is
  effectively impossible, and `put_bytes()` errors rather than retries on one.
  Writing in 4 kB pieces is 13x slower than one `writeBin` **to /dev/null**
  (144 vs 1907 MB/s), which looks like an easy win and is not: with a real
  reader on the other end the chunk size makes no measurable difference
  (151.8 ms at 4 kB vs 155.0 ms at 1 MB on 32 workers), because the reader, not
  the writer, is the constraint. Leave the SIGPIPE hardening alone.

---

## Commit 11 — `src/` has no header dependency tracking (SIGNIFICANCE: MEDIUM — fragility)

**DONE**, and it is the reason the rest of this work could be trusted.
`./configure` probes `-MMD -MP` (recording `header-deps=` in the Makevars
metadata and the configure summary) and `src/Makevars.in` / `src/Makevars.win`
include the generated `.d` files.

One trap worth recording: R's shlib run supplies **no** explicit make goal, so
including the `.d` files makes the first target of the first `.d` the default
goal and the build stops after one object. `.DEFAULT_GOAL := all` before the
include pins it. Verified: touching `ParamTable.h` now rebuilds exactly the
seven translation units that include it, and nothing else.

R's default rules compile `foo.cpp -> foo.o` with **no** dependency on the
headers it includes. Editing `src/ParamTable.h` and re-running `R CMD INSTALL`
recompiles only the `.cpp` files whose timestamps changed; every other
translation unit keeps an object file built against the *old* struct layout. The
result is an ODR/layout mismatch and a segfault at a random address.

This is not theoretical — it happened during this sweep. The Commit 1 prototype
segfaulted at `address (nil)` and looked like a null-pointer bug in my new code;
it was a stale `.o`. `rm src/*.o` and a rebuild produced correct, identical
likelihoods with no source change. `project_emc2_rlf_graded_mesh.md` records the
same trap ("touch all 8 TUs that include `model_RLF.h` TRANSITIVELY or you get an
ODR segfault") and `project_emc2_roup_build_op.md` records it again ("rm src/*.o
or the header change is ignored").

Fix it once, in `src/Makevars.in`:

```make
PKG_CXXFLAGS += -MMD -MP
-include $(OBJECTS:.o=.d)
```

(and add `*.d` to `.gitignore` / the install `.Rbuildignore` cleanup). Failing
that, an explicit `$(OBJECTS): $(wildcard *.h)` line is blunt but correct.
Anything is better than a documented instruction to remember to delete object
files, because the failure mode is a benchmark that silently measures the wrong
binary or a crash that reads as a code bug.

---

## Commit 12 — Regression harness for the prologue (SIGNIFICANCE: MEDIUM)

**DONE.** `tests/testthat/test-param-table-prologue.R` (40 assertions) and
`WorkingTests/bench_likelihood_prologue.R`. Both hang off one new export,
`pt_prologue_oo()`, which runs exactly the lane `calc_ll_oo` uses and either
times it or hands back the parameters and bound verdicts it produced — so the
benchmark needs no instrumentation compiled into the sampler, and the tests can
pin the fast lane against `get_pars_c_batch_wrapper_oo`, an independent
implementation of the same mapping that takes none of the short cuts.

Covered: intercept-only, 2- and 3-level factors, interactions, no-intercept,
continuous covariates, compressed *and* expanded design matrices, split
(pre-sum) transforms, constants, per-row bound verdicts (one factor level out of
bounds must reject exactly its own trials), invariant bounds, an identity trend
against the un-trended lane, and particle independence (mapping one particle
alone must equal mapping it inside a batch).

Commits 1, 2, 4 and 6 all change the code that turns a particle row into
per-trial natural-scale parameters. That code has no direct test today: it is
covered only transitively, through likelihood values.

Add `tests/testthat/test-param-table-prologue.R`:

* For a matrix of designs (intercept-only; factor with 2–3 levels; interaction;
  continuous covariate; `uses_self`; split transform; `constants`; a trend), call
  `get_pars_c_batch_wrapper_oo` on a fixed particle matrix and snapshot the
  natural-scale columns. This pins mapping/transform/bounds independently of any
  model kernel.
* Assert the row-constant flag agrees with a brute-force re-scan
  (behind `EMC2_PT_CONST_CHECK`).
* Assert that a design with a continuous covariate is *not* treated as
  single-cell — that is the one input shape that would silently break Commit 1/2.
* A likelihood-equality test between a model fitted with and without a trend that
  multiplies by 1, so the trend-gated and un-gated lanes are checked against each
  other.

Also add a `WorkingTests/bench_likelihood_prologue.R` that reports the §0.2
phase split, so the next person can see in one command whether the mapping tax
has crept back.

---

## Commit 13 — Small R-side items (SIGNIFICANCE: LOW)

**DONE except two bullets.** `accept_rate()` is now incremental (it slides the
per-subject change counts instead of recomputing an n_subjects x 200 block:
~0.5 ms/iteration at 140 subjects, and the cache validates itself against the
array it was built from, so it can only ever be a speed-up);
`marginal_ll_from_grid()` uses a `pmax` fold; `merge_group_level()` uses
`matrix()`.

`.emc_cores_now()` was cached and then **reverted**. The reallocation contract
is that a core released by a finished sibling is visible to the *next* caller,
and the pool's own tests assert exactly that (`test-worker-pool.R:337,340,361`);
a 0.25 s cache made them fail — intermittently, which is worse. The listing is
tens of microseconds against a ~100 ms iteration, so there was nothing to win
either. If a networked tmpdir ever makes it expensive the fix is to put the
arena on a local disk, not to answer with a stale count. The comment in the code
now says so.

`extend_obj()`'s geometric growth is not done — see "What is NOT done" at the top.

Filed together because none of them individually clears 1 % of a fit, and §0.1
bounds what they can ever be worth.

* `accept_rate()` (`R/messaging.R:160`) slices a 3-D array
  (`samples$alpha[1, , start:end]`) and does an `n_subjects x 200` comparison
  **every iteration** just to feed the progress bar. Keep a running acceptance
  counter in `pm_settings` instead — `update_pm_settings` already computes the
  information. Only matters at high subject counts and only when
  `verboseProgress = TRUE`.
* `.emc_cores_now()` (`R/fitting.R:1105`) calls `list.files()` on the core-control
  directory once per sampler iteration per chain. Cache the count and re-read at
  most a few times a second; a `readdir` per iteration on a network `tmpdir`
  would be a real cost.
* `marginal_ll_from_grid()` (`R/sampling.R:146`) uses
  `apply(log_terms, 1L, max, na.rm = TRUE)` — an `apply` over rows of an
  `np x K` matrix. `matrixStats::rowMaxs` is not available (no dependency), but
  the same `pmax`-over-columns fold already used a few lines below in
  `new_particle` (`R/sampling.R:1204-1207`) is both faster and dependency-free.
* `merge_group_level()` (`R/sampling.R:1754`) builds
  `do.call(cbind, rep(list(c(tmu_nuis)), ncol(subj_mu)))` where
  `matrix(tmu_nuis, nrow = sum(is_nuisance), ncol = ncol(subj_mu))` does it
  without the intermediate list. Nuisance models only.
* `extend_obj()` allocates a fresh full-size array per stage block per sample
  element; with `saved + writes` in the thousands and 140 subjects this is the
  memory high-water mark of a long fit. Growing geometrically (double the
  capacity, track a fill pointer) would turn `n_blocks` reallocations into
  `log2(n)` of them. Worth doing only if a long fit is actually seen to thrash.

---

## Measured and rejected

Recording these so they are not re-attempted.

* **Skipping the kernel for an entirely out-of-bounds particle** (Commit 2's
  "free bonus"). Measured after the fact on a 400-trial LBA with `v ~ lM*E`,
  `B ~ E*lR`: the fraction of particles with no in-bounds row is 0 % at proposal
  sd 0.1 and 0.5, 0.2 % at 1.0 and 6.0 % at 2.0. Real preburn clouds sit at the
  low end, and the branch would have to reproduce the contaminant/guess mixing
  exactly to be safe. Not worth it.
* **Rewriting the worker transport to non-blocking reads** to get a receive
  deadline. The completion FIFO already in the pool gives the same thing without
  touching the framing code (see Commit 10).

* **Caching `make_pt_mapper` across likelihood calls.** Fixed per-call cost is
  0.60 ms against 0.501 ms per particle (§0.4) — 1.2 % at 100 particles. Not
  worth the staleness risk.
* **`dnormP` on the LBA path alone.** 1–2 %. Listed as Commit 5 only because it
  is three lines and should be worth more elsewhere; do not expect a headline
  number from it on ballistic models.
* **Chasing the log-space fallback in the BA family.** Instrumented counters show
  199 968/200 000 pdf rows and 200 000/200 000 cdf rows already take the cheap
  natural path at realistic parameters. The acceptance thresholds
  (`BA_ACCEPT_RAW`, `natural_cdf_safe`) are not costing anything in the regime
  fits actually visit.
* **`.is_valid_ll_cache` on the hot path.** 0.6 % of a censored likelihood call.
  Already fixed; nothing left.
* **`gibbs_step`.** 0.8 % of a preburn fit, confirming the earlier measurement.

---

## Suggested order

*(This is the order the work was actually done in, and it held up: Commit 11
first made every later benchmark trustworthy, and Commit 12 caught the one-ulp
FMA difference in Commit 2 the same hour it was introduced.)*

1. **Commit 11** (header deps) first — everything after it is a benchmark, and
   without it the benchmarks are not trustworthy.
2. **Commit 12** (prologue tests) second, so Commits 1/2/4/6 land against a net.
3. **Commit 1** (row-constant columns) — the measured 17 % end-to-end win.
4. **Commit 4** + the `build_plan` overlap fix, then **Commit 5**.
5. **Commits 8, 9, 10** — correctness and fragility, independent of the above.
6. **Commit 6**, then **Commit 7**.
7. **Commit 2** (cell resolution + per-cell kernel memoisation) as its own
   project once 1 and 12 are in; it is the largest remaining structural win and
   the one that generalises the BAwDD memoisation to every model.
8. **Commit 3** (burn-stage particle adaptation) whenever there is time to run
   the A/B properly. It is the only item here that could plausibly beat Commit 1,
   and the only one I have no number for.

---

## Reproducing the measurements

The scripts used are small and are worth keeping in `WorkingTests/`:

* Fit profile (§0.1): `make_emc` on replicated `forstmann`, `Rprof(interval =
  0.01)` around `run_emc(stage = "preburn", stop_criteria = list(iter = 150))`.
* Likelihood benchmark (§0.2–0.5): single subject, `type = "single"`,
  `rt_resolution = NULL`, `rt` jittered so compression does not collapse the
  rows, `EMC2:::calc_ll_manager()` in a loop. **Set realistic parameter values** —
  at `sampled_pars()`'s zeros, `t0 = exp(0) = 1 s` exceeds every RT, every row
  floors at `min_ll`, and the benchmark measures the wrong path (17.9 ms instead
  of 50.7 ms, and it is 3x too fast for the wrong reason).
* Phase split: temporary `clock_gettime(CLOCK_PROCESS_CPUTIME_ID)` timers around
  the three blocks of `PtMapper::prepare()`, exported via a
  `// [[Rcpp::export]]` accessor.
* A/B: install into two temp libraries and alternate runs
  (`R CMD INSTALL --library=...`), with `.libPaths()` set inside the script and
  the resolved namespace path printed — `library(EMC2)` otherwise picks up the
  stale system install (`project_emc2_verification_trap.md`).
* **Two temp libraries are not a fair A/B for an R-only change.** Chasing an
  apparent 8 % single-core gain from the Commit 10 broadcast change cost an hour:
  `.Rlib_optbase` 5.83/5.84/5.83/5.82 s against `.Rlib_optnew` 5.36/5.35/5.36/
  5.40 s, seven consecutive interleaved pairs, with **byte-identical `EMC2.so`**
  (same md5) and **bit-identical draws**. It is not real. Restoring the old code
  *inside one library* by `assignInNamespace()`-ing a `deparse`/`sub`/`eval`
  copy of the two functions gives 5.46/5.48 and 5.39/5.63 — no difference, same
  `alpha` checksum. Two separate `R CMD INSTALL` runs differ in lazy-load
  database and bytecode layout by several percent on whole-fit CPU. For an
  R-side change, patch the namespace in-process; keep the two-library A/B for
  C++ changes, where it is the only option.
* Diffing two installed builds: `deparse()` every function in both namespaces
  and `diff` the dumps. That is what established the two libraries above
  differed *only* in the intended functions, which is what made the timing
  discrepancy worth explaining rather than shrugging at.
* Pool profiling: `options(emc2.sampler_profile = TRUE)` puts a per-iteration
  data frame on `attr(samples, "sampler_profile")` with `shared_serialize`,
  `request_send`, `response_wait`, `worker_max`, `worker_sum`, `shared_bytes`,
  `private_bytes`, `wire_bytes` and `workers`. `worker_max` next to
  `response_wait` is the tell for whether a pool is bound by work or by wire.
* `perf` is unavailable on this box (`perf_event_paranoid = 4`) and there is no
  `gdb`, so all attribution above is from compiled-in counters and timers, not
  sampling.
