# FPE solver performance: audit, changes, and plan

Status of the Fokker–Planck solver family (ROU / ROUp / BOU) and the RLF
propagator, as of 2026-09-01, branch `playground`, baseline `b5c0ed25`.

Every number below is a measurement on this branch, not an estimate. Where
something is an estimate — the expected gains in the plan — it says so.

**Benchmark.** `forstmann`, subject 1, 810 trials, single-subject `make_emc`,
`rt_resolution = NULL`. AMD EPYC Genoa, AVX-512, one core.
**Timing.** `CLOCK_PROCESS_CPUTIME_ID` / `user.self` throughout, minimum of
repeated runs — this box is one core under shared load and elapsed time swings
~30%. Phase attribution inside the march is `rdtsc` against a measured TSC rate.

---

## 1. What landed

Files: `src/fpe_race.h`, `src/fpe_solver.h`, `src/fpe_bou.h`. All uncommitted.

### Root cause

The time march was **latency-bound on store-to-load forwarding** — not on
arithmetic and not on memory bandwidth. Both Thomas sweeps read `Q[i-1]` and
`Q[i+1]` back out of memory one iteration after storing them, which puts a
~6-cycle store-forward on top of an 8-cycle multiply-add dependency chain.

Two measurements rule out the alternatives:

* cost per cell-step was flat (2.0–2.8 ns) from `nx` = 64 to `nx` = 1024, so it
  is not bandwidth;
* the eight-lane path was only 2.06× faster *per lane* than the scalar one
  (863 µs vs 419 µs per key), so it is not arithmetic throughput.

Instrumented cycles per vector row per time step, `nx` = 512, 8 lanes, 2.45 GHz:

| Phase | Shipped | After |
|---|---:|---:|
| RHS assembly | 7.0 | — (fused) |
| Forward sweep | 14.1 | 7.5 (incl. RHS) |
| Back substitution | 12.1 | 4.0 |
| **Total** | **33.2** | **11.5** |

The back sweep's 4.0 is exactly one fnmadd latency, i.e. the floor.

### The three changes

All are pure evaluation-order rewrites of the same scheme.

1. **Carry both recurrences in registers.** The value the next row needs is
   already in a register; reading it back from memory is what costs.
   (`fpe_race.h` lane march, `FPE_Tri::solve` in `fpe_solver.h`, `fpe_bou.h`
   lane march.)
2. **Fuse the Crank–Nicolson right-hand side into the forward sweep.** A rolling
   three-register window of the pre-step `q` removes an entire pass over `Q` and
   the operator, and the `RHS` array with it. `q[i+1]` is read before `q[i]` is
   overwritten, so the in-place update is safe. Valid only when the two passes
   read the *same* operator, so a moving boundary keeps the split path in
   `fpe_race.h` (`fuse_rhs = !any_moving`); in `fpe_bou.h` the old and new
   operators live in different arrays, so it fuses unconditionally. New
   `FPE_Tri::solve_cn` covers the scalar path.
3. **Vectorise the moving-operator factorisation.** The collapsing-bound branch
   called `set_lane_factor` — a scalar per-lane loop over the mesh — at every
   time step. The ROUp branch already had a vectorised all-lane version; both
   now share `set_all_lane_factors`. This is the single largest item for
   collapsing bounds, where the operator never amortises.

### Results

Kernel level, 12 parameter keys, `t_max` = 1.5 s, minimum of three runs:

| Path | Shipped | After | Gain | Output |
|---|---:|---:|---:|---|
| ROU, fixed bound | 5.86 ms | 2.22 ms | 2.6× | bit-identical |
| ROU, collapsing bound | 22.42 ms | 7.70 ms | 2.9× | 1e−12 relative |
| ROUp, pulse drift | 10.77 ms | 7.59 ms | 1.4× | bit-identical |
| BOU, 8-lane batch | 21.14 ms | 15.48 ms | 1.4× | bit-identical |

End to end, through `design` → `make_emc` → the C++ likelihood:

| Model | Shipped | After | Gain | vs. analytic peer |
|---|---:|---:|---:|---|
| ROU | 3.70 ms/particle | 1.65 ms/particle | 2.25× | 19.5× → 8.7× RDM |
| BOU | 753 ms/particle | 514 ms/particle | 1.47× | 391× → 267× DDM |

Log-likelihoods identical to every printed digit in both cases. ROU carries
~0.25 ms/particle of non-solve overhead, so an infinitely fast solver would
still leave it at ~1.3× RDM.

The collapsing-bound 1e−12 is FMA contraction in the factorisation: it now runs
the same vectorised code the ROUp branch has always used. Everything else
reproduces the previous log-density and log-survivor grids exactly.

### Verification

Every kernel change was validated by a checksum over the full `log_pdf` and
`log_S` grids of a fixed key set, compiled against the old and new headers in
the same binary shape. "Bit-identical" above means all 15 printed digits agreed.

`test-rou.R`, `test-fpe-fht.R`, `test-fpe-bounded-ou.R`, `test-rlf-fht.R`,
`test-rlf-model.R` all pass with `NOT_CRAN=true` and `EMC2_TEST_LEVEL=full`.

The four `test-roup.R` errors are **pre-existing** — they fail identically on the
baseline library. They are a `design_model` scoping artefact of calling
`test_file` outside `test_check`, not a numerical failure.

### Rebuild trap

R's makefiles ignore header dependencies, and `fpebou::SolveCache` is reachable
from more translation units than the obvious ones. Rebuilding only the "obvious"
TUs after changing that struct produced a clean build and a segfault (classic
ODR mismatch). **Delete `src/*.o` after any change to these headers.**

---

## 2. Why BOU reads as 267× DDM

BOU pays **(design cells) × `n_sv` × `n_sz`** marches. Smith & Ratcliff anchor
the decay at the start point, so `drift(x) = (v + βz) − βx`: a start-point
quadrature node changes the *operator*, not just the initial condition. With
both `sv` and `SZ` free that is 49 solves per cell.

Same design, same data, varying only which variability parameters are free:

| Free parameters | Solves/cell | ms/particle | vs. DDM |
|---|---:|---:|---:|
| `sv`, `SZ`, `st0` | 49 | 516 | 268× |
| one of `sv`/`SZ`, + `st0` | 7 | 72 | 37× |
| `st0` only | 1 | 64 | 33× |
| none | 1 | 62 | 32× |
| DDM, same design | — | 1.93 | 1× |

Any earlier timing with `SZ` fixed would have been ~7× cheaper than the headline
number, which is the most likely explanation for a remembered figure much closer
to DDM.

**Note the bottom three rows.** They differ by 10 ms despite a 49× spread in
solve count. There is a **~55 ms/particle floor that does not shrink with `nx`** —
it survives at `nx` = 16. That is `nt_min` = 512 plus per-step, per-lane output
bookkeeping, and it is 28× DDM's entire likelihood before a single mesh cell is
touched. See step 3.

**What it is not:** a cache-lookup problem. `fpebou::SolveCache::find` was a
linear scan called ~500k times per particle, which looked like the obvious
culprit. Replacing it with a hash index made **no measurable difference** on a
six-cell design. The index is kept — it is the one part of the loop whose cost
grows with the design, and it mirrors what `fperace::SolveCache` already does —
but it is not where the time was. The comment in `fpe_bou.h` records this so
nobody re-derives it.

---

## 3. BOU's defaults over-resolve by about ten-fold

Against a converged reference (`nx` = 768, `dt` = 2.5e−4, 11×11×11 nodes), the
shipped defaults are accurate to 3.7e−6 nats/trial — roughly two orders tighter
than anything a fit can resolve. The dominant term is the time step:
`dt_target` = 5e−4 is about eight times finer than needed and costs 3× on its own.

| Setting | ms/particle | max abs error, nats/trial |
|---|---:|---:|
| shipped — 384 / 5e−4 / 7×7 | 731 | 3.7e−06 |
| `dt` → 2e−3 only | 251 | 5.2e−06 |
| `nx` → 192, `dt` → 2e−3 | 167 | 2.0e−05 |
| 192 / 2e−3 / 5×5 | 77 | 2.0e−05 |

With the march rewrite already in, the last row measures **62.5 ms/particle** —
**12×** the shipped default.

Held across five parameter regimes: strong leak, large `sv`/`SZ`, small and
large `a`, and β ≈ 0 against the DDM oracle. Worst error anywhere in that set was
4.1e−05 nats/trial.

---

## 4. Measured and rejected

Recorded with the measurement that killed each one, so the next pass does not
spend a day rediscovering them. Several were headline claims in the deleted
proposal document.

**A modal / spectral solver replacing the time march (claimed 15–35×).** Not
credible. The eight-lane path is already partly throughput-bound (4 → 8 lanes
costs 1.57×, not 1.0×), so wider vectors do not pay either. More fundamentally,
a truncated modal or shift-invert Krylov expansion resolves the *slow* modes and
needs the fast ones at small `t` — the density has to be right from ~10 ms.
Segmenting the horizon restores the fast modes but restores the step count with
them. The RLF propagator already *is* shift-invert Krylov, and that switch was
worth 1.3–3.2×, not 15–35×.

**Removing RLF horizon bucketing (claimed 2–3×).** Zero effect at α = 1.7: the
bucket is unoccupied and the split is a no-op, exactly as `model_RLF.h` predicts.
The code comment already documents a re-measurement after Richardson landed
showing the split *is* needed at α 1.1–1.5, where one slow trial otherwise
unresolves the boundary layer for every trial sharing its parameters.

**Toeplitz acceleration of the RLF operator build.** Already implemented.
`build_rlf_operator_graded` fills `band[k] = phi(k * uniform_h)` and indexes it
by face separation inside the uniform core, so the transcendentals are already
O(n) rather than O(n²). Only the assembly bookkeeping remains — see step 4.

**Richardson extrapolation in `nx` for ROU.** The scheme is cleanly second order,
so the extrapolation is valid — it just does not pay. A (192 + 288) pair costs
3.60 ms for the same mean error as a single `nx` = 384 at 3.22 ms. It *does*
improve the max (tail) error about 3×, so it is an accuracy tool, not a speed one.

**Cache-set staggering of the lane buffers.** At `nx` = 512 with 8 lanes every
buffer is exactly 32 kB — the L1 size — so all eight map to identical cache sets.
Padding them apart changed nothing measurable. Null result.

**Retuning `RLF_KRYLOV_MIN`.** 12% faster at α = 1.7, slightly slower at α = 1.3.
Parameter-dependent, so not a general win. The real target is the rung count —
see step 5.

**Retuning ROU's grid.** Unlike BOU there is almost no slack: the (`nx`, `nt`)
Pareto frontier puts the shipped (512, 256) close to optimal, and the best
available trade is ~1.6× for a 30% error increase. `dt` is already irrelevant
below 4e−3 because `nt_min` = 256 floors it.

---

## 5. Plan forward

Ordered by value per unit of risk.

### Step 1 — Re-balance the BOU defaults

*Expected ~9× on top of what landed. Touches `model_BOU.h` defaults only.*

Move `emc2.bou_nx` 384 → 192, `emc2.bou_dt` 5e−4 → 2e−3, and `n_sv` / `n_sz`
7 → 5. The options already exist, so this is a defaults change and a
documentation change, not a code change. The risk is entirely in whether five
scenarios generalise.

Tests needed before flipping the defaults:

1. **Parameter-net convergence.** A Latin hypercube of ~150 points over
   (`v`, `a`, `Z`, `sv`, `SZ`, `st0`, `beta`, `t0`) spanning the prior box, each
   scored against the converged reference. Gate on **max** |Δll|/trial, not mean
   — the failure mode is a corner, not a drift. Five scenarios put that at
   4.1e−05; a threshold of 1e−03 leaves margin and is still far below what a fit
   resolves.
2. **Independent oracle at β = 0.** At zero leak BOU *is* the DDM, so compare
   against Navarro–Fuss rather than against a finer BOU. This is the only check
   in the set that is not self-referential. Sweep `sv`, `SZ`, `st0` through their
   ranges; gate on max |Δlog f| over the central mass of each response.
3. **Grid-is-not-binding invariant.** The existing rationale is that the grid is
   set just fine enough that quadrature dominates. Re-establish it at the new
   node counts: at 5×5, refining (`nx`, `dt`) must move the likelihood *less*
   than going 5×5 → 7×7 does. If it does not, the grid is the binding term and
   `nx` = 192 is too coarse.
4. **Recovery equivalence.** Fit one simulated data set at old and new defaults;
   require posterior means to agree within MCSE. Expensive, but the only test
   that answers the question a user actually asks.
5. **Edge behaviour.** `SZ` at its cap (start range touching a boundary), `st0`
   large relative to `t0`, very large `beta`, `a` at both extremes. These are
   where a coarser mesh degenerates rather than merely losing accuracy.
6. **Timing guard.** A cheap regression assertion so the defaults cannot
   silently drift back.

### Step 2 — Move the BOU decay anchor to the midpoint

*Expected 5–7×, and more accurate rather than less.*

With `anchor = a/2` instead of `anchor = z`, the drift becomes `(v + βa/2) − βx`
— independent of the start point. The operator stops depending on `z`, so the
whole `n_sz` quadrature collapses into a **single march seeded with the uniform
start density**. That is an exact integration replacing a 5- or 7-node
approximation, so it is both cheaper and better: `fpe_seed` already handles
`z_lo ≠ z_hi` with exact cell-overlap averaging at t = 0, with no
frozen-coefficient warm-up at all.

Physically this moves the zero of the restoring force from the start point to
the midpoint of the decision interval, so accumulation is no longer pinned where
it began and decays through it.

Shape of the change — `anchor_at_z` already exists and is plumbed; it is
hardcoded `true` at `model_BOU.h:129`:

* flip that flag (and expose it as a model option, see test 6);
* drop `z` from `fpebou::Key`'s solve identity and add the start *interval*;
* replace `fpe_seed(m, p.z, p.z, …)` in `bou_solve` with the interval form;
* `bou_mix` loses its inner `iz` loop entirely.

Tests needed:

1. **Operator independence.** Assert directly that two keys differing only in
   `z` now build identical operators. The whole saving rests on this property,
   so it should be a unit test rather than an inference.
2. **Uniform seed equals the quadrature limit.** Under the *same* midpoint
   anchor, the single uniform-seeded march must agree with an `n_sz` = 31
   Gauss–Legendre sweep over point starts to solver accuracy. This isolates the
   seeding from the anchor change.
3. **β = 0 still reproduces the DDM with `SZ`.** At zero leak the anchor is
   irrelevant, so this must be unchanged from today — a good regression
   trip-wire for the key/seed rewiring.
4. **Quantify the model difference.** This is a different model, not just a
   faster one. Report |Δll|/trial between the two anchors across the prior box so
   the change is documented with a known size, and note it in NEWS.
5. **Recovery under the new anchor.** Confirm `beta` is still identified, and
   specifically check the `v`–`beta` correlation at fixed `a`: the effective
   drift now carries a common `βa/2` offset instead of a per-node `βz` one,
   which could tighten that ridge.
6. **Keep the old anchor reachable** as e.g. `BOU(anchor = "start")` so published
   fits stay reproducible.

### Step 3 — Kill BOU's resolution-independent floor

*~55 ms/particle today; becomes the binding cost once steps 1 and 2 land.*

Once the marches are cheap, what is left is per-step bookkeeping: at every one of
`nt` ≥ 512 steps, for every lane, BOU takes two `std::log` calls and five
`push_back`s, then stores a full time grid per entry. None of that scales with
`nx`, which is why the floor survives at `nx` = 16.

Two candidates, in order of appeal:

* Store the fluxes in the natural scale and take logs at lookup instead of at
  every grid point. Cheap, but only wins when queries are rarer than grid
  points — needs checking.
* Give BOU the sparse output path ROU already has, recording only at the times
  actually asked for. Harder: the `st0` convolution shifts the query times per
  row, so the sparse plan has to be built per design cell rather than per row.

Also worth checking whether `nt_min` = 512 is still the right floor once
`dt_target` moves.

### Step 4 — Trim the RLF operator assembly

*Operator build is 18% of the RLF solve; expect a few per cent back.*

The transcendental work is already Toeplitz-exploited. What remains is
bookkeeping: `build_rlf_operator_graded` allocates and zero-fills two dense
`nf × nf` scratch tables on every call — about 270 kB per build, twice per key
for the Richardson pair, ~6.5 MB per particle at n ≈ 130.

* `Xt` is entirely derivable and can go. It satisfies
  `Xv(p,q) = Gv(p,q) * (face[p] - face[q])` exactly, for both orderings, so one
  multiply replaces a ~135 kB table.
* `G` should move into the reusable `RLF_KrylovWork` scratch rather than being
  reallocated per call.

Neither changes an arithmetic result, so the test is a checksum on the assembled
operator plus the existing RLF suite.

The larger idea — keeping the operator in Toeplitz-plus-correction form so the
apply is O(n log n) — founders on the dense LU of (I − γL), which the
shift-invert Arnoldi needs and which is not Toeplitz. Worth a sketch only if
steps 1–3 land and RLF becomes the slowest model.

### Step 5 — Break RLF's two-eigensolve floor

*`rlf_reduce_modes` is 44% of the solve at ~2 rungs each; worth ~20%.*

RLF profile at α = 1.7, `nx` = 70, 12 keys (24 solves — Richardson pairs):

| Phase | Share |
|---|---:|
| `rlf_reduce_modes` (dhseqr + dtrevc per rung) | 44% |
| `rlf_arnoldi_extend` | 22% |
| `build_rlf_operator_graded` | 18% |
| shifted LU factor | 8% |

The Krylov search exits on *agreement between successive rungs*
(`have_previous && moved < RLF_KRYLOV_TOL`), so **two eigensolves are guaranteed
by construction** — measured at 1.96 rungs per solve, i.e. already at that floor.
Tightening or loosening `RLF_KRYLOV_TOL` cannot help; nor can retuning
`RLF_KRYLOV_MIN`, which only shifts where the two rungs sit.

The machinery for a single-rung exit already exists: `rlf_reduce_modes` computes
a genuine shift-and-invert defect as a modal curve and `modes.residual_at()`
evaluates it at any time. It is deliberately *not* a stopping rule — the comment
on `RLF_KRYLOV_RESID` says both diagnostics "are diagnostics rather than stopping
rules" — but it *is* trusted to gate the adaptive fast path
(`rlf_fast_path_safe`), so the question is open rather than settled.

The study that would settle it:

1. **Instrument every rung.** Over a parameter net in (α, `v`, `b`, `A`,
   horizon, `nx`), record at each rung both the residual over the probe sweep
   *and* that rung's true error against a converged reference (`m` = `m_cap`, or
   a dense `expm` at small n).
2. **Ask whether a threshold exists.** Plot true error against residual. The
   question is not "are they correlated" but "is there a residual threshold with
   **no false accepts**" — a single accepted-but-wrong rung is a silently wrong
   likelihood, which is worse than the 20%.
3. **Judge against the right yardstick.** The existing note records that the
   two-rung test is conservative by three orders: 1e−05 of movement leaves
   9e−06 of error against a spatial error near 1e−02. Any single-rung rule only
   has to beat the *spatial* error, not the current tolerance.
4. **Record the negative result either way.** If no clean threshold exists, that
   closes the item permanently and explains the two-rung design to whoever asks
   next.

---

## Not pursued

**Cross-particle key batching.** Gathering keys across the particle loop would
fill lanes exactly (12 keys currently chunk 8 + 4, and a 4-chunk costs 64% of an
8-chunk). Measured ceiling ~1.2–1.35× for a substantial restructuring of the
particle loop. Not worth it.

**Wider SIMD (16 or 32 lanes).** The 4 → 8 lane step already costs 1.57×, so the
path is partly throughput-bound on Zen 4's double-pumped AVX-512 and wider
vectors would not recover per-lane cost.
