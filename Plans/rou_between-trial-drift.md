# Between-trial drift variability (`sv`) for ROU

Status: **plan only, nothing implemented.** Written 2026-09-01 against branch `playground`
(`b5c0ed25`), measured on the build installed from this worktree (`src/EMC2.so`, 2026-08-31).

Target: normally-distributed between-trial variability in the accumulation rate of
`R/model_ROU.R`, with a `posdrift = TRUE/FALSE` option matching `RDMSWTN`
(`R/model_RDM.R:836`).

---

## 1. The structural fact that decides everything

`A` (start-point range) is free because it is an **initial condition**. The Fokker-Planck
equation is linear in `p`, so a uniform start is one solve with a smeared seed
(`fpe_seed`, `src/fpe_models.h:490`).

`sv` is a **parameter mixture**. `v` sits inside the spatial operator (`fpe::build_op`), so

    f(t) = INT f(t | v) phi(v; mu, sv) dv

is a mixture of *different PDEs*. Three escape routes were considered and all fail:

* **Augmented 2-D state `(x, v)`.** With no `v`-dynamics the generator is block-diagonal,
  i.e. literally `n_v` decoupled 1-D solves. The "2-D solve" *is* the quadrature.
* **Algebraic reduction.** `Y = X - v/k` turns the OU into a driftless mean-reverting
  process, but it moves `v` into *both* the start point and the boundary location
  (`b - v/k`), so each node remains a genuinely different first-passage problem. The key
  still needs four numbers.
* **Marginal Gaussian trick.** For `k = 0`, `X(t) = z + v t + s W(t)` with `v ~ N` is
  marginally Gaussian with inflated variance, but the *process* is not Markov, so no FPE
  represents it. This is exactly why LBA/RDM integrate over `v` analytically rather than
  solving one modified diffusion.

**Cost is therefore linear in the number of quadrature nodes, and the whole engineering
problem is minimising that node count.** Everything below follows from that.

---

## 2. Where it plugs in - one function, not the kernels

The race kernels (`drou_raw`, `prou_raw`, `rou_logS_at_t`, `drou_scalar`, `prou_scalar`),
the censoring / endpoint-query path, and `rou_prepare_rows` must **not** change. They all
consume a cache `Entry` holding `log_pdf` / `log_S` on a time grid. Make the mixture *be*
that Entry.

* `fperace::Key` gains `sv` and a `posdrift` flag. **The parent key count is unchanged**,
  so `rou_prepare_rows`' pass-1 linear scan over keys (`src/model_ROU.h:160`) does not
  degrade and `row_group` is untouched. This is the main argument for expanding inside
  `cache_get_batch` rather than in `rou_prepare_rows`.
* Inside `cache_get_batch` / `cache_get`, a missing key with `sv > 0` expands into `n`
  **transient child keys** (`v_j = v + sv * z_j`, `sv = 0`), which join the ordinary lane
  batch alongside every other key. The parent Entry is assembled as
  `sum_j w_j * pdf_j`, `sum_j w_j * surv_j`, then `safe_log`.
* Mix on `surv`, never on `1 - cdf`. `FPE_Result::surv` exists precisely to keep relative
  accuracy far into the tail (`src/fpe_solver.h`, `FPE_Result` comment); a positive-weight
  sum of positive quantities preserves that, `1 - cdf` throws it away.
* Children are **not** inserted into `C.index`; only the mixed parent is stored. Memory per
  entry is unchanged.

Two consequences worth stating plainly:

* `sv = 0` yields one node of weight 1, so results are **bit-identical to today** and the
  cost is unchanged. This is the no-regression property that makes the feature safe to ship
  on by default.
* Every downstream consumer works unmodified: contaminants, `pGuess`, truncation
  normalisers, the GSL integrands on the unknown-winner branch, and `EndpointQueryPlan`.

### 2.1 The one alignment problem

Lanes in `fpe_solve_batch_ou_lanes` run **independent time schedules**
(`src/fpe_race.h:339`), and `fpe_seed` picks a per-model `t_seed`. Siblings share a horizon,
hence `nt` and the block schedule - but with `A = 0` (point start) they get different
`t_seed`, so their grids do not line up for elementwise mixing.

`fpe_seed` already accepts `t_seed_force` (`src/fpe_models.h:493`) for exactly this reason
and **no caller currently uses it**. Force `t_seed = min` over the sibling set; the analytic
seed is valid at any sufficiently small time, so no lane is seeded later - and therefore
less accurately - than it would have been alone. Grids then match exactly and the mixture is
an exact elementwise weighted sum with no interpolation error. With `A > 0` the seed returns
`t0 = 0` for every lane anyway and this is a no-op.

---

## 3. The quadrature rule - the highest-leverage decision

Measured against a 40-node reference (ROU, `k = 1`, `b = 1.3`, `A = 0.3`, `t0 = 0.15`,
`t` from 0.25 s to 5 s so the far tail is covered). Entries are `max |d log f|`.

| rule | mu=2, sv=0.5 | mu=2, sv=1.0 | mu=1, sv=1.0 |
|---|---|---|---|
| Gauss-Legendre in `v` over +/-8 sd | n=21 -> 9e-8; **n=9 -> 2e-2** | n=9 -> 2e-2 | n=9 -> 6e-3 |
| Gauss-Legendre in probability space (`v = Phi^-1(u)`) | n=31 -> 2.6e-4 | - | - |
| plain Gauss-Hermite | n=7 -> 3e-5, then **stalls** | **stalls at 2e-2** | **stalls at 1.7e-1** |
| **half-range Gauss (rule for the truncated weight)** | **n=7 -> 2.8e-4** | **n=9 -> 1.7e-3** | **n=7 -> 1.2e-3** |

Survivor accuracy for the recommended rule, same setup (`max |d log S|`):

| nodes | mu=2, sv=0.25 | mu=2, sv=0.5 | mu=2, sv=1.0 | mu=1, sv=1.0 |
|---|---|---|---|---|
| 3 | 2.0e-2 | 4.7e-2 | 8.9e-2 | 6.1e-2 |
| 5 | 7.9e-5 | 3.8e-3 | 1.1e-2 | 6.2e-3 |
| 7 | 1.3e-6 | 5.4e-4 | 1.5e-3 | 6.5e-4 |
| 9 | 1.5e-8 | 1.7e-5 | 2.0e-4 | 6.9e-5 |
| 11 | 1.7e-10 | 5.8e-7 | 2.4e-5 | 7.3e-6 |

Three findings worth keeping:

* **Probability-space Gauss-Legendre is a trap.** `Phi^-1(u)` has unbounded derivatives at
  the endpoints, so it converges algebraically (~1/n^2). It is the obvious first thing to
  write and it is the worst option on the table.
* **Plain Gauss-Hermite stalls, and the plateau is truncation bias**, not quadrature error -
  GH integrates over `v < 0` whatever you do with the weights. It is fine only when
  `mu/sv >~ 3.5`.
* **A Gauss rule built for the truncated weight itself converges geometrically and is
  truncation-exact.** It needs 2-3x fewer nodes than anything else, which is 2-3x off the
  entire feature's runtime cost. This is the single most valuable decision in the plan.

### 3.1 Constructing the half-range rule

The standardised weight `phi(u) * 1{u > -a}` depends on the single shape parameter
`a = mu / sv`. Build the three-term recurrence by **discretised Stieltjes** on a 300-400
point Gauss-Legendre discretisation of the weight, then Golub-Welsch on the resulting Jacobi
matrix. `src/gh_quad.h::make_gh_rule` already does Golub-Welsch from a Jacobi matrix, so
this is a sibling of existing code, not a new technique in the repo.

Cost is an `n x n` symmetric eigensolve plus `O(nN)` inner products, ~20 us, against ~0.4 ms
**per node** of PDE solve - under 1%. Cache on `(n, a)` in the `SolveCache` and clear per
particle; the number of distinct `a` equals the number of design cells.

`posdrift = TRUE` guarantees `mu > 0` (the `v` transform stays `exp`), so `a > 0` always and
the recurrence is well conditioned. For `posdrift = FALSE` there is no truncation, so the
**existing cached `gh_rule(n)` is exactly right and free**.

---

## 4. Node-count policy

Anchor the target on the solver's own bias. At the shipped `nx = 512` the PDE carries
**3.9e-4** in `log f` (measured against `nx = 2048, dt = 1e-3`). Quadrature error below that
is wasted money.

Measured node requirement scales with `cv = sv / mu`:

| `cv` | nodes for <~1e-4 |
|---|---|
| 0 | 1 (exact passthrough) |
| <= 0.15 | 5 |
| <= 0.3 | 7 |
| <= 0.6 | 9 |
| <= 1.0 | 11 |

**Default `n = 9`**, exposed as `getOption("emc2.rou_sv_nodes", 9L)` and read in
`rou_configure_grid` alongside the existing `emc2.fpe_*` knobs. Fitted `sv/v` in RDM/LBA
usually lands at 0.3-0.5, which `n = 9` covers with margin.

**Adaptive `n` from `cv` is worth doing here specifically because convergence is
geometric.** Set the switch thresholds so both sides are under ~1e-5 and the induced
discontinuity in the log-likelihood is < 0.01 nats over 800 trials - far below sampler
noise. This matters in practice because `sv`'s default is `log(0)`, so burn-in draws near
zero cost nothing at all. Note this reasoning does *not* transfer to an algebraically
converging rule, where the jump would be visible; it is a property of the rule chosen in
section 3, not a general licence.

---

## 5. Measured cost

All CPU time (`proc.time()[["user.self"]]`, min/mean of repeated runs), one core,
AVX-512 present, load ~0.6, `nx = 512`, `dt = 4e-3`, `grade = 8`, `tgrade = 32`.

* Marginal cost of one extra distinct key in the batched regime: **0.42-0.57 ms**
  (horizon 1.6 s, 135 RTs per key).
* Lane batching buys ~2.2x per key going from 1 lane to 8 (1.10 ms -> 0.50 ms), **but
  realistic designs already saturate it**, so it does not absorb the node multiplier.
* Forstmann, single subject, `v~lM, k~1, B~E+lR, A~1, t0~1`, `s` fixed - 12 distinct keys
  over 1620 rows: **4.3 ms per particle likelihood** today.

| nodes | keys | per-particle ll | multiplier |
|---|---|---|---|
| 1 (`sv = 0`) | 12 | 4.3 ms | 1.0x |
| 5 | 60 | ~22 ms | 5.1x |
| 7 | 84 | ~31 ms | 7.2x |
| 9 | 108 | ~40 ms | 9.3x |

Small designs do better than the naive multiplier: a 2-key design goes 2 -> 18 keys for only
~3.3x wall time, because the extra nodes fill lanes that were previously padding. Multi-core
fits divide this by the worker pool.

**So yes, it is genuinely ~9x - but 9x is the floor, not a symptom of a bad design.** The
two things holding it at 9x rather than ~25x are the half-range rule (section 3) and the
`sv = 0` / small-`sv` adaptivity (section 4).

---

## 6. Things that do *not* buy the cost back

Checked and rejected, with the measurement that killed each:

* **Coarsening `nx` to pay for the nodes.** The hope was that mixing smooths the density
  enough to tolerate a coarser mesh. It does not. `max |d log f|` vs `nx = 2048, dt = 1e-3`:

  | `nx` | `sv = 0` | `sv = 0.5` (9 nodes) |
  |---|---|---|
  | 192 | 3.81e-3 | 3.65e-3 |
  | 256 | 2.22e-3 | 2.12e-3 |
  | 384 | 9.26e-4 | 8.74e-4 |
  | 512 | 3.92e-4 | 3.63e-4 |
  | 768 | 1.14e-4 | 1.15e-4 |

  Identical, because discretisation error is a systematic bias that survives positive-weight
  averaging. The trade is bad in that direction anyway: `nx` 512 -> 192 is 2.0x faster for
  10x worse.
* **`dt`.** `nt_min = 256` floors it; raising `dt_target` does nothing. (Same finding as the
  ROUp work.)
* **Sharing solves across design cells** via a common `v` grid. Only valid when cells share
  `(k, b, A, s, boundary)`, which `B ~ E + lR` breaks immediately, and it trades a per-cell
  Gauss rule for an interpolatory one - worse accuracy for a modest constant factor.
* **Sharing the tridiagonal factorisation across nodes.** `v` enters only `a0`, additively,
  so the face Peclet numbers shift by a constant - but the Scharfetter-Gummel weights are
  nonlinear in `P`, and even with central differencing `(I - c L(v))` would still need a
  factorisation per node. Nothing to reuse.
* **Reduced-basis / POD in `v`.** Building the basis costs solves and the solution manifold
  moves every MCMC step.

---

## 7. Correctness traps this will hit

1. **`fpe_x_lo_ou` is unsafe for negative drift** (`src/fpe_models.h:393`). It sizes the
   domain floor at `min(z_min, v/lambda)`. Currently unreachable because `rou_key` rejects
   `v <= 0`, but `posdrift = FALSE` makes negative nodes legal and `v = -0.5, k = 1e-3`
   gives `v/k = -500`: the domain explodes and all barrier resolution is lost. The bound
   should be the mean **at `t_max`**,
   `min(z_min, z_min * exp(-k t_max) + (v/k) * (1 - exp(-k t_max)))`,
   which is exact, finite, and reduces to the existing BM limit as `k -> 0`.
   **Fix this before enabling `posdrift = FALSE`.**
2. **`rou_key` returns false on `v <= 0`** (final line of `rou_key`, `src/fpe_race.h`). Must
   become conditional on `posdrift`.
3. **`defective_upper_tail`.** ROU declares non-defective (`src/race_dispatch.cpp:128`), and
   that stays correct: with `k > 0` the OU is positive recurrent (stationary
   `N(v/k, s^2/2k)`), so it hits `b` almost surely **even for negative `v`**. Unlike RDM,
   `posdrift = FALSE` here creates no atom at infinity; only `k = 0` exactly is defective.
   This makes `posdrift = FALSE` better behaved for ROU than for RDM. But a strongly
   negative node leaves large `S(t_max)` inside a finite horizon, so exercise the truncation
   normaliser path on it.
4. **No-flux far field with negative drift.** `F_0 = 0` reflects mass that physically should
   escape to `-inf`. Domain sizing keeps this near 1e-9, but validate at the extreme node of
   a wide `sv`.
5. **Sparse-vs-full path.** At forstmann scale the crossover
   (`num_chunk_queries * 8 <= num_chunk_steps`) selects the **full-grid** path, so the
   `t_seed` alignment of section 2.1 is load-bearing rather than an edge case. Children must
   also carry the parent's query-time set so the sparse path stays exact when it *is*
   selected.

---

## 8. Implementation checklist

### C++
* `src/fpe_race.h`
  * `Key` + `KeyHash` + `Key::finite()` gain `sv` and `posdrift`.
  * `rou_key` scales `sv` by `1/s` alongside `v`; relaxes the `v > 0` guard under
    `posdrift = FALSE`.
  * New `rou_drift_rule(v, sv, posdrift, n) -> {nodes, weights}` with an `(n, a)` cache.
  * New `rou_solve_mixture()`, used by **both** `cache_get` (scalar / GSL path) and
    `cache_get_batch` (row path), so the two cannot drift apart.
  * Child expansion + mixture assembly inside `cache_get_batch`; keep siblings adjacent so
    they land in one 8-lane chunk.
* `src/fpe_models.h` - fix `fpe_x_lo_ou` per trap 1; thread `t_seed_force` from the sibling
  group.
* `src/model_ROU.h` - `rou_configure_cache` reads `emc2.rou_sv_nodes`; `RouColIdx` gains
  `sv`; the collapse-column offset moves from the hard-coded `emc2col::rou::Binf` to a value
  stored on `SolveCache`.
* `src/col_registry.h` - `sv`-bearing ColSpecs (see decision 9.1).
* `src/race_dispatch.cpp` - `_SV` / `_IO` suffix handling, tested before the existing
  boundary suffixes so the variants compose.

### R
* `R/model_ROU.R`
  * `ROU(drift_variability = FALSE, posdrift = TRUE, ...)`.
  * `sv` p_type: `exp` transform, `[0, Inf)`, default `log(0)`, `exception sv = 0`.
  * `posdrift = FALSE` flips `v` to the natural scale per `.rdmswtn_v_scale`
    (`R/model_RDM.R:992`).
  * `.rou_cols` picks up `sv`; roxygen parameter table row and `@param`.
* `rROU` / `.rfun_ROU_R` - **no C++ change needed**: draw `v_i ~ N(v, sv)` (truncated at 0
  when `posdrift`) in R and pass the drawn vector to the existing `rou_hit_times_vec`.
* `make_emc` cost warning: a trend on `v` already makes every trial a distinct key; `sv`
  multiplies that by the node count, so the message should say so.

### Tests (`tests/testthat/test-rou.R`)
* `sv = 0` reproduces current values exactly (bit-identical, not `tolerance`).
* Mixture vs Monte Carlo simulation, both `posdrift` settings.
* Node-count convergence table from section 3 reproduced as a regression guard.
* `k = 0, sv > 0` against the analytic `RDMSWTN` kernel - a strong independent check, since
  that limit *is* RDMSWTN.
* Negative-drift domain sizing after the trap-1 fix.

---

## 9. Open decisions

### 9.1 `sv` column placement
`RDMSWTN` makes `sv` unconditional, which here would add a sampled parameter to *every*
existing ROU design and break stored fits and tests. Recommendation: make it **opt-in** via
`ROU(drift_variability = TRUE)` with a `_SV` `c_name` suffix and a second ColSpec per chart
(6 total), storing the collapse-column offset on `SolveCache` rather than hard-coding
`emc2col::rou::Binf`. Slightly more registry code, zero churn for existing users. The
alternative - match RDMSWTN and always carry `sv` - is simpler and more uniform if the churn
is acceptable.

### 9.2 `sv` under the non-rate charts
Cleanest is that `sv` is always variability in the **canonical `v`**, applied *after*
`rou_map_to_rate`. Under `equilibrium` that reads as variability in `theta` at fixed `t_k`;
under `curvature`, variability in the `v` implied by `tstar`. The alternative is to restrict
`drift_variability` to `parameterization = "rate"` initially.

### 9.3 Is `sv` identifiable here at all - check before spending the 9x
In a leaky accumulator, `k` already produces slow errors and a flattened right tail, which
is the job `sv` does in the RDM; `A` already covers fast errors. Expect `sv` and `k` to trade
off. A recovery study on `(v, sv, k, A)` is probably a better use of the first hour than the
implementation.

### 9.4 Out of scope for now
`ROUp` (`v_S`, `v_T`) takes the same mechanism but doubles the decision space (variability on
one channel or both, correlated or not). Do ROU first.
