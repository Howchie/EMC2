# WS1 allocator specification (from the code)

Source: `R/move_allocation.R` and `src/alloc_need.cpp` at `0a5773dd` (branch `ws1/allocator`). The consumers are in `R/sampling.R`, `R/shear_move.R`, `R/collapsed_move.R` and `R/subject_hmc.R`. Everything below was checked against the code. "Chosen" means a constant with no derivation in the code or its history.

## 1. Where it runs

| When | What | Code |
|---|---|---|
| Every iteration, every chain, stage ≠ preburn, `alloc_on` | records into `sampler$alloc$window` | `sampling.R` 1061, 1078, 1178, 1224 |
| Every rebuild (`add_proposals`, after the fit checks and the shear and collapsed plans) | measures, solves and writes the plan | `.emc_alloc_rebuild` (`fitting.R:506`) |
| Preburn | nothing: no window and no plan | `.emc_alloc_rebuild` returns at once |

`alloc_on = stage != "preburn" && (.emc_shear_active || .emc_collapsed_active || (HMC enabled && supported))`. In a default fit of a `standard` or `diagonal-gamma` model without nuisance or group design, all three are on.

The rebuild schedule:
- burn and adapt rebuild every block (B = 100 transitions);
- sample rebuilds at 100, 200, 400, ... transitions (`emc2.sample_rebuilds = "doubling"`), so B doubles: 100, 100, 200, 400, ...

The plan's outputs, each fixed for the block:
- `shear$plan$rounds`: expected rounds per group coordinate (2p);
- `collapsed$plan$rounds`: expected collapsed moves per iteration;
- `subject_hmc$subjects[[s]]$rate` (r_s) and `alloc$particle$rate` (q_s, shared across chains as the maximum).

How the rates are used:
- Rounds: `floor(x)` moves, plus one more with probability `x - floor(x)`. Shear draws this independently per coordinate; the collapsed move draws one count per iteration.
- q and r share one uniform per subject: the particle step runs if U < q, HMC if U ≥ 1 − r (`sampling.R:927-951`).

## 2. Measurement (`.emc_alloc_measure`, per chain j)

### Window
These are sums since the last rebuild:
- **Group, "rest of the iteration"**: from the previous stored state to just before the shear. Coordinates are `d = (mu, ½ log diag V)`, accumulated as `G += d dᵀ` with count `nG`.
- **Subjects**: from the previous stored alpha to just before HMC. `A_s += d dᵀ` and `D_s += dd ddᵀ`, where dd is the deviation jump (alpha jump minus the mu jump), with count `nA`.
- **Particle step jumps**: `P_s` and `nP_s`, from the steps that ran.
- **Costs**: `cpu_rest`; per family, `fam[[f]] = (cpu, own, iters, moves)`; HMC `cpu_hmc_s` and `hmc_att_s`.

Costs are wall seconds × the chain's cores (`.emc_reserved_cost`). Particle and HMC costs are worker CPU (marginal).

`stored = (no window or w$n < MIN_PAIRS = 10)`. This holds at the first rebuild after preburn, which then uses whole-iteration jumps between the last 200 stored draws as the "rest".

**Production** (stage = sample and the previous rebuild was also sample):
- group G and nG pool over every sample block (`alloc$pool`, undecayed);
- variances use every production draw of every chain.

Outside production, variances use every chain's last `w$n` draws.

### Group coordinates (2p)
- `z = (Dinv mu, lsd)`, where `Dinv = solve(shear$plan$dirs)`.
- `var_k`: pooled over chains and draws.
- `esjd_k = (T G Tᵀ)_kk / nG`.
- Rest-of-iteration rho, as a lower bound:

  `rho_lo(esjd, var, n) = clip01(1 − m (1 + 2√(2/n)))`, with `m = esjd / (2 var)`.

  The factor 2√(2/n) is derived: two standard errors of a mean of n squared Gaussian jumps, ignoring autocorrelation of the jumps.

### Shear factor `g_k` (if shear is active)
The input is `last = shear$plan$last`:
- in production, the pooled record `pool_k = 0.5 pool_{k-1} + block_k`;
- otherwise, the last block only.

Records are matched across rebuilds by direction (|cos| ≥ 0.9) and rescaled. Log-sd records also pool across chains in production.

The factor:
- `g_obs = clip01(1 − (E/n) / (2 var))`.
- `g_ref`: the round-weighted mean of `g_obs` over the measured coordinates of the same kind (mean directions or log sds). Failing that, either kind's; failing that, `SHEAR_PRIOR_G = 0.5`.
- `g = (n g_obs + n0 g_ref) / (n + n0)`, with `n0 = ROUND_EVIDENCE = 20`. An unmeasured coordinate gets `g_ref`.

Saturation `a` (`.emc_alloc_saturation`):
- Solve `D = (1 − g) Σ_k w_k (1 − a^k) / (1 − a)` for a ∈ [0, g] by 60-step bisection. Here `D = N / Σ K / (2 var)` (net squared jump per iteration) and w is the empirical distribution of rounds per iteration.
- With no iteration of ≥ 2 rounds, a = g (geometric).

### Collapsed factor `c_fac` (if collapsed is active)
- Records: in production, `pool = 0.5 pool + block`; otherwise the last block, or failing that the previous block of the same stage if it had ≥ 5 moves.
- `c_obs = clip01(1 − diag(T E Tᵀ)/n / (2 var))`.
- `c_fac = (n c_obs + n0 mean(c_obs)) / (n + n0)`, shrinking towards its own mean over coordinates. Unmeasured: 0.5 for every coordinate.
- Saturation is computed as for shear, with one K row for all coordinates.

### Subjects
- `W_s` is the HMC reference whitening of the subject's likelihood precision against Vbar (the mean group V over the window). Failing that, the identity.
- `s_var`: pooled variance of `W_s alpha_s` over the variance draws.
- `s_rho = rho_lo(diag(W A Wᵀ)/nA, s_var, nA)`, the rest including the particle step at the old rate.
- Particle factor (with ≥ 10 steps): `p = clip01(1 − diag(W P Wᵀ)/nP / (2 s_var))`, and `p_dev` the same against the deviation variance. Otherwise it is carried over from the last rebuild.
- `c_part = cpu_part / nP`.
- Rest without the particle step: `t = min(1, s_rho / (1 − q_old (1 − p)))`.
- Deviation statistics: `d_rho` comes from D_s against `d_var` (lower bound). `r_bar` is the mean of all `d_rho`.
- Group floor (perfect subject re-equilibration): `g_floor = clip(g_rho, (g_rho − r_bar)/(1 − r_bar))`.

### HMC factor `h_s` (`.emc_alloc_hmc_factor`)
- From the subject's evidence: the acceptance-weighted squared jumps, decayed 0.5 per production rebuild, and reset when eps changes.
- Requires `n ≥ min_evidence = 10` (option `emc2.subject_hmc_min_evidence`).
- `h = clip01(1 − diag(W M Wᵀ)/n / (2 s_var))`, shrunk as `(n h + 20 h_ref) / (n + 20)`, where `h_ref` is the n-weighted mean over subjects. `h_dev` is the same against `d_var`.

### Costs (`.emc_alloc_marginal`)
- Per family: `c = own / moves` if moves ≥ 5, else NA.
- `fixed = (cpu − own) / iters` if iters ≥ 10, else carried from the last rebuild (0 the first time).
- `C_rest = (cpu_rest − particle cpu of priced subjects) / n`.
- Guesses for unmeasured costs (chosen): a shear round costs `0.1 C_iter`, a collapsed move `0.1 C_iter`, an HMC transition `3 C_iter / n_subjects` (or its evidence's CPU per move). `C_rest` defaults to 1 if not finite.

## 3. Units (`.emc_alloc_units`)

- **Group coordinate k**: `shear(rho = g_rho, g, a, cost = c_round, floor, cap)` if shear is on, else `fixed(rho = g_rho)`.
  - Evidence bounds, given B transitions to the next rebuild: `floor = max(0, 20 − n_ev)/B` and `cap = min(4, max(20 − n_ev, 2 n_ev)/B)`.
  - `MAX_ROUNDS = 4` (chosen).
- **Subject s**: `subject(rho, p, c_p, h, c_h, d_rho, p_dev, h_dev, q_old, r_old)`.
  - `rho = t` and `p` are set only if the particle step is measured ("priced"). Otherwise `rho = s_rho`, `p = NULL`, and the particle cost stays in C_rest.
  - `h` is set only if the subject is eligible and has ≥ 10 HMC transitions. In the calibration pass ("optimistic"), an eligible subject without evidence gets h = 0.
- **Joint (collapsed)**: `idx` = every group coordinate, `g = c_fac`, `a = c_sat`, `cost = c_coll` (or 0.1 C_iter), `fixed = coll_fixed`, `floor` from the evidence bounds, and `cap = min(evidence cap, max(4, 2 × measured rate))`.

## 4. The solve

**Composition.** x expected moves with per-move factor g and saturation a give
`F(g, x, a) = (1 − f) F_k(m) + f F_k(m + 1)`, where m = ⌊x⌋, f = x − m, and `F_k(k) = 1 − c (1 − a^k)` with `c = (1 − g)/(1 − a)`. When a = g this is g^k.

**Subject kernel factor under coupled selection** (`.emc_alloc_dev_factor`):
`Φ(q, r; P, H) = (1 − q − r + b) + (q − b) P + (r − b) H + b P H`, with `b = max(0, q + r − 1)`.
- A subject's predicted rho is `t · Φ(q, r_total; p, h)`.
- Its floor ("level") is `t · p`, the particle step alone.

**`.emc_alloc_solve(units, C_rest, joint, movable, shear_fixed, hmc_fixed)`**:
1. Candidate collapsed rates: `xc ∈ {j_floor} ∪ grid`, where the grid is eighths up to 4, then halves up to the cap. Each candidate multiplies the group units' rho by `F(c_fac, xc, c_sat)`.
2. `movable` units are those some move can lower. Only movable units set the target range.
3. Targets: 400 log-spaced tau values between the best reachable and the unallocated maximum, plus every unit's own level.
4. For each target t, each unit's cheapest allocation reaching `rho ≤ t` (`.emc_alloc_unit_need`):
   - shear: `x = F⁻¹(g, a, t/rho)`, floored and capped;
   - subject: over q ∈ {1, 7/8, …, 0} (or the held `q_fixed`), the least r with `rho_i Φ(q, r) ≤ min(t, level_i)` on every coordinate (C++ `emc_alloc_subject_need_cpp`). The cheapest `q c_p + (r − r_min) c_h` wins.
5. Objective at each target: `tau(t) × (C_rest + Σ unit costs + xc c_coll + fixed_c min(1, xc) + activation)`.
   - Activation = `shear_fixed (1 − Π_k (1 − min(1, x_k))) + hmc_fixed (1 − Π_s (1 − min(1, r_s)))`.
6. The minimum over targets and xc is taken. It is replaced by the zero plan (floors only) unless `objective ≤ MIN_GAIN (0.9) × objective0`.

**`.emc_alloc_solve_units(u)`** wraps the solve in three layers:
1. If the shear fixed cost is > 0: solve with shear (frames paid by activation) and without (exploration floors kept, and their frame cost). With shear wins only if it is better *and* ≤ 0.9 × without's objective0.
2. If HMC is available:
   - "without": h removed, and the group credited at q = 1, r = 0 through the deviation ratio.
   - "with": the plain solve at the measured credit, then each regime r ∈ {0.5, 1}. Every eligible subject runs HMC at ≥ r, and the group is credited by `.emc_alloc_dev_ratio`.
   - The best is taken, then `settle`: if the chosen (q, r) move the deviations less than the credit assumed, re-solve "own" (credited at its rates, rates held) against "kept" (q = 1, credited, r held) and take the cheaper.
   - With beats without only if it is better *and* ≤ 0.9 × without's objective0.
3. Group credit: `rho' = f + (rho − f) R` for R ≤ 1. For R > 1 it is the saturating form `1 − (1 − rho) exp(−(rho − f)(R − 1)/(1 − rho))`. Here `R = Σ base Φ(q, r; p_dev, h_dev) / Σ base Φ(q_old, r_old; …)` and `base = d_rho / (1 − q_old (1 − p_dev))`.

**`.emc_alloc_rebuild`**:
1. HMC plans per chain.
2. Measure.
3. Optimistic pass: subjects the optimistic plan would use without evidence get an eps bracket and 20 pilot transitions, forked.
4. Each chain's plan.
5. q = max over chains. Chains whose own q was lower are re-solved with `q_fixed`.
6. Records written to `alloc$last` and `alloc$history`.

## 5. Constants

| Constant | Value | Status |
|---|---|---|
| `MAX_ROUNDS` | 4 per coordinate per iteration | chosen |
| `MIN_PAIRS` | 10 (window iterations; fixed-cost iterations; particle steps) | chosen |
| `SHEAR_PRIOR_G` | 0.5 (shear and collapsed prior factor) | chosen |
| `ROUND_EVIDENCE` n0 | 20: factor shrinkage weight, and the exploration target | derived rationale (relative SE √(2/n) ≈ 0.3 at n = 20), value chosen |
| `EVIDENCE_DECAY` | 0.5 per production rebuild (shear, collapsed, HMC) | chosen |
| `MIN_GAIN` | 0.9 (one in each of the solve, shear with/without and HMC with/without) | chosen |
| `HMC_REGIMES` | {0.5, 1} | chosen |
| `PARTICLE_RATES` | eighths | chosen grid |
| UCB factor | 1 + 2√(2/n) | derived (Gaussian jumps, independent) |
| tau grid | 400 points | chosen (resolution) |
| collapsed grid | eighths to 4, then halves | chosen |
| unmeasured cost guesses | 0.1 C_iter (round and collapsed), 3 C_iter / n (HMC) | chosen |
| marginal-cost minimum | 5 moves | chosen |
| collapsed rate cap | max(4, 2 × measured rate) | chosen |
| HMC `min_evidence`, `pilots` | 10, 20 | chosen (options) |
| stored-window draws | 200 | chosen |
| HMC warm-up rate, moves | 0.5, 20 (`subject_hmc.R`) | chosen |

## 6. Branches

**Default-on**
- shear, collapsed and HMC units;
- production pooling;
- evidence floors and caps;
- regimes (`emc2.alloc_hmc_regimes = TRUE`);
- settle;
- the shared-q second pass;
- the shear and HMC fixed-cost with/without layers;
- MIN_GAIN (×3);
- the optimistic calibration pass.

**Reachable in a default fit only in particular states**
- `stored` window: the first burn rebuild.
- `Dinv` failure fallback: a singular shear direction matrix (defensive).
- "fixed" group units: never with shear on, since `g_fac` is always finite. They occur only when `emc2.group_move ≠ "shear"`.

**Not reachable in a default fit (scaffolding or test-only)**
- `emc2.alloc_subject_floor = "refresh"` (the `refresh` paths in `.emc_alloc_solve` and the C++ kernel).
- `emc2.subject_hmc_rate` / `fixed_rate` and `emc2.subject_hmc_eps` / `fixed_eps` branches.
- `kind = "hmc"` units: only tests construct them (`.emc_alloc_as_subject`).
- `emc2.shear_interleave > 1` (`sampling.R`, not the allocator's).

## 7. Findings from writing this down (each is checked in the report)

- **S1: the evidence cap does not mean what it says in production.** The cap is "no faster than twice the rate the evidence was gathered at", computed as 2 n_ev / B. In production, n_ev is a 0.5-decayed pool while B doubles at every rebuild. For a steady rate x, n_ev ≈ x (B_k + B_{k−1}/2 + …) ≈ x B_k (1 + 1/4 + 1/16 …), so `2 n_ev / B_{k+1} ≈ 1.33 x`. A move's rate can grow by at most about a third per production rebuild, not double.
  - The floor has the same mismatch.
  - Only log-sd records pool across chains: three chains' evidence makes their cap ≈ 4x.
- **S2: the lower confidence bound is applied only to the rest of the iteration.** The move factors g, h and p are point estimates, so every prediction is biased low by the UCB. At nG = 200 the factor is 1.4; once production pooling holds thousands of pairs it is about 1.05.
- **S3: `pred` existed nowhere.** No record held the rho the plan predicts, so calibration was impossible. `.emc_alloc_predict` now attaches it to every solve, and the rebuild record stores `pred_group`, `pred_subject`, `pred0_group` and `pred0_subject`.
