# WS5 Review: Evaluation of the Ensemble Move Cost and Sizing Workstream

**Review of:** `plans/workstreams_2026-10-13/ws5_ensemble_cost.md` and `plans/workstreams_2026-10-13/ws5_report.md` (Branch `ws5/ensemble-cost`, worktree `/data/work/EMC2_ws5`, HEAD commit `5eb122d8`).  
**Author:** AI Reviewer  
**Date:** 2026-10-13  

---

## 1. Executive Summary

The task assigned in `ws5_ensemble_cost.md` was two-fold:
1. **Phase 1 (Zero-Behavioral Change):** Make the serial group ensemble update (`.emc_ensemble_move`) and in-iteration Gibbs group update (`gibbs_step_standard`) cost as little as possible by porting the serial arithmetic to C++, without altering RNG streams, discrete choices, or parameter draws beyond floating-point rounding under `-ffast-math`.
2. **Phase 2 (Principled Sizing):** Replace the static constants $M = 32$ (pool size) and $R = 10$ (Gibbs rounds) with an adaptive between-block rule based on measured importance efficiency, coordinate movement, and CPU prices.

### Evaluation Scorecard
* **Overall Grade:** **Very Good (Phase 1) / Needs Revision (Phase 2)**
* **Phase 1 Deliverables:** Completed to exceptional mathematical and numerical standards. The serial chain-process overhead fell by **$5\times\text{--}10\times$** across test fixtures with exact RNG preservation and zero selection mismatches.
* **Phase 2 Deliverables:** Partially completed; conceptually flawed in its first iteration and incomplete in its second. The first rule collapsed $M \to 2$ on Forstmann, reducing worst-coordinate ESS/s by 33%, and stalled on `lba7r`. The revised rule corrected immediate accounting defects but was left unit-tested only without full benchmark or ablation completion (`nback31` was omitted entirely).
* **Overlooked Optimizations:** Several major computational and architectural optimization opportunities were missed, ranging from repeated $O(n)$ Cholesky factorizations of frozen matrices in every iteration to redundant matrix inversions in inverse-Wishart sampling and unbatched worker likelihood dispatch.

---

## 2. Evaluation of Phase 1 (Make It Cheap)

### 2.1 Implementation Quality (`src/ensemble_move.cpp`)
The C++ translation (`7d875047`, `ad3183b4`) successfully captured the core serial loops:
* `emc_ensemble_draw_pools_cpp`: Generates $M-1$ fresh draws per subject from frozen proposal mixtures $(1-w)\mathcal{N}(\mu_s, V_s) + w t_4(\mu_s, c V_s)$ and computes log proposal densities via Mahalanobis distances.
* `emc_ensemble_stack_cpp`: Constructs the $(n \cdot M) \times p$ pooled state matrix, computing baseline importance weights $w_0 = \ell - \log q$.
* `emc_ensemble_rounds_cpp`: Executes all $R$ rounds of alternating conjugate group Gibbs updates and Gumbel-max categorical candidate reselections.
* `emc_gibbs_standard_cpp` & `emc_gibbs_diag_gamma_cpp`: Ports the main iteration's conjugate Gibbs updates for unblocked and block-diagonal covariance structures.

### 2.2 Numerical Identity & RNG Stream Preservation
The validation harness (`benchmarks/ensemble-cost/sbs_ens.R`) demonstrated exemplary precision across sample-stage runs:
* **Seed Invariance:** Bit-identical `.Random.seed` states after every call.
* **Selection Invariance:** Zero selection mismatches across 200 iterations of Forstmann, 600 iterations of T30, and 1,600 iterations of `lba7r`.
* **Floating-Point Drift:** Maximum relative difference between C++ and pure R reference implementations stayed within $10^{-12}\text{--}10^{-15}$, well within expected compiler/LAPACK tolerance.

### 2.3 Measured Performance Impact
Serial in-iteration execution time was reduced by roughly an order of magnitude:

| Fixture | Parameter $p$, Subjects $n$ | Baseline Move Time | Head Move Time | Baseline Gibbs | Head Gibbs | Speedup |
|---|---|---|---|---|---|---|
| **Forstmann** | $p = 4, n = 19$ | 11.5 ms | **0.8 – 1.7 ms** | 7.2 ms | **0.1 – 0.2 ms** | $\approx 8\times$ |
| **T30** | $p = 13, n = 40$ | 23.4 ms | **1.9 – 2.6 ms** | 12.0 ms | **0.16 – 0.19 ms** | $\approx 10\times$ |
| **lba7r** | $p = 20, n = 106$ | 26.0 ms | **5.3 – 5.5 ms** | 2.1 ms | **0.20 – 0.21 ms** | $\approx 5\times$ |

### 2.4 Phase 1 Limitations
1. **Omission of Group Designs and Latent Factor Hierarchies:**
   Only `standard` (without group designs) and `diagonal-gamma` were ported to C++. Models featuring group design matrices (`sampler$group_designs`) and latent-factor models (`factor`, `infnt_factor`, `SEM`) were relegated to the legacy R path (`.emc_ensemble_rounds_r`). These hierarchies still incur full R-level scratch sample store replication (`.emc_ensemble_store`) and interpreted matrix operations.
2. **Missing In-iteration Memory Optimization:**
   Large scratch matrices ($p \times nM$) continue to be dynamically allocated inside the round loop rather than reused across rounds.

---

## 3. Evaluation of Phase 2 (Size It on Principle)

### 3.1 Empirical Sweep Findings (`sweep_mr.R`)
The exploration of stationary states under varying $M \in [2, 128]$ and $R \in [1, 20]$ revealed two critical facts regarding the sampler's geometry:
1. **Means Decorrelate in One Round:** The initial Gibbs draw (round 1) conditions on the inherited subject states. Because the conjugate Gibbs conditional is unconstrained, group means ($\mu$) jump to the independence limit ($J \approx 2.0$) in round 1 across all tested models and pool sizes.
2. **Variances Saturate Early:** Reselections in rounds 2–3 produce an incremental jump of $\Delta J \approx 0.3\text{--}0.5$ in log standard deviations ($\log \sigma$), after which exploration completely plateaus. Extra rounds ($R = 4\dots 20$) contribute zero measurable displacement.
3. **Model Efficiency Disparity:** On Forstmann, pool efficiency was high ($\text{eff} \approx 0.9$), whereas on T30, pool efficiency was minimal ($\text{eff} \approx 10^{-4}\text{--}10^{-3}$), showing that the proposal $q_s$ struggled to match high-dimensional posteriors.

### 3.2 Post-Mortem of Rule v1 (`c910a3c4`)
The initial adaptation rule produced catastrophic regressions:
* **Forstmann:** Collapsed pool size from $M = 32 \to 2$ at the first rebuild. Lag-1 autocorrelation of the worst coordinate ($\log \sigma_{t0}$) worsened from $0.03 \to 0.17$, reducing worst ESS per CPU-second by 33% for an insignificant wall-time saving of 2 s out of 35 s.
* **lba7r:** The chain stalled and was terminated after 31 minutes ($2.4\times$ slower than the static control).

The report diagnosed two causes:
1. *Row-share pricing error:* Pricing pool rows by their fraction of total particle round wall time attributed latency, scheduling delay, and idle cores to pool rows, overstating member cost multiple-fold.
2. *Coordinate averaging:* Averaging displacement across coordinates allowed the immediately-independent group means to hide the loss of exploration on the slowest log-sd coordinate.

### 3.3 Assessment of Revised Rule v2 (`141df019`)
The revision implemented two corrections in `R/ensemble_plan.R`:
* Measured worker CPU time per subject (`extra_times`) was recorded directly and normalized by chain cores.
* The objective was changed to a max-min criterion over the slowest coordinate:
  $$\max_{x', R'} \min_k \frac{T_k + G_k(R') \frac{m(x')}{m(x)} - G_k(R)}{C + c_{\text{member}}(x' - x) + c_{\text{round}}(R' - R)}$$

### 3.4 Unresolved Flaws in the Phase 2 Formulation
Despite the fixes in Rule v2, three deeper theoretical flaws remain:

#### 1. The Additive Jump Assumption
The objective assumes that the chain's lag-1 jump $T_k$ shifts linearly by the delta in ensemble displacement: $\Delta T_k \approx \Delta G_k$. However, $T_k$ is the combined product of particle steps, HMC, shear moves, and Gibbs updates. The ensemble move is not an independent additive noise source; it interacts multiplicatively with subject state diversity.

#### 2. The Independence Proposal Fallacy
The model approximates pool displacement via an independence-sampler acceptance probability:
$$m(x) = \frac{1}{n} \sum_{s=1}^n \frac{x e_s}{1 + x e_s}$$
where $x = M - 1$.
* In high-efficiency regimes (Forstmann, $e_s \approx 0.9$), $x e_s \gg 1$ for even modest $x$. Thus, $m(x)$ rapidly approaches 1, implying diminishing returns for any $x > 1$.
* But the Neal ensemble update is a **data-augmentation representation** of the joint posterior $\prod_s p(\alpha_s | y_s)$. Collapsing $M \to 2$ reduces the subject posterior support to a binary choice, stripping the conjugate Gibbs covariance update of the dispersion needed to move $\log \sigma$.
* Conversely, in low-efficiency regimes (T30, $e_s \approx 10^{-3}$), $G_k \approx 0$. The numerator becomes independent of $x$, meaning the objective maximizes solely by minimizing the denominator, again collapsing $M \to 2$.

#### 3. Unchecked Deliverables
The revised rule was only unit-tested. Deliverables 4 and 8 (continuation runs across fixtures, verifying worst-coordinate ESS/s, and completing the Forstmann ablation) were left pending or queued, leaving the production safety of Rule v2 unproven.

---

## 4. Key Overlooked Optimizations

A detailed audit of the codebase reveals several performance and architectural optimizations that were overlooked in WS5:

### 4.1 Recomputing Frozen Proposal Cholesky Factors Every Iteration
* **Location:** `src/ensemble_move.cpp`, lines 282–287 (`emc_ensemble_draw_pools_cpp`).
* **The Defect:** In every iteration, for every subject $s \in 1\dots n$, the C++ code executes:
  ```cpp
  arma::mat Ls;
  if (!arma::chol(Ls, V, "lower")) ...
  L.slice(s) = Ls;
  logdet(s) = 2.0 * arma::sum(arma::log(Ls.diag()));
  ```
  However, the proposal covariance matrices $V_s$ are **frozen for the entire block** (100–400 transitions). They do not change between iterations.
* **Performance Impact:** On `lba7r` ($n = 106, p = 20$), the sampler executes 106 Cholesky factorizations of $20 \times 20$ matrices and 106 log-determinants every iteration. This accounts for the majority of the **7.1 ms** pool drawing time.
* **Remedy:** Precompute and store $L_s$ and $\log|V_s|$ once in `add_proposals()` when $V_s$ is formed. Pool drawing in C++ should merely draw standard normals $z$ and evaluate $L_s z$.

### 4.2 Triple Matrix Inversion in Inverse-Wishart Bartlett Sampling
* **Location:** `src/ensemble_move.cpp`, lines 73–94 (`riwish`) and lines 122–127 (`gibbs_standard`).
* **The Defect:**
  1. `riwish` draws a Wishart matrix $W \sim \mathcal{W}_d(df, S^{-1})$ via Bartlett factor $ZC$ ($W = ZC^T ZC$).
  2. It inverts $W$ to return $\Sigma = W^{-1}$ as the inverse-Wishart draw.
  3. Immediately upon receiving $\Sigma$, `gibbs_standard` calls `arma::inv_sympd(Sigma_inv, Sigma)` to re-invert $\Sigma$ back into $\Sigma^{-1}$.
  4. In Step 3 of the Gibbs rounds, it calls `arma::chol(Lg, tvar)` to evaluate $(x - \mu)^T \Sigma^{-1} (x - \mu)$.
* **Performance Impact:** $W$ **is** $\Sigma^{-1}$ by definition. Computing $W \to W^{-1} \to (W^{-1})^{-1}$ executes two redundant $d \times d$ matrix inversions per block per round. Furthermore, since $W = (ZC)^T (ZC)$, the square-root factor of $\Sigma^{-1}$ is already available without taking another Cholesky.

### 4.3 Redundant Large Heap Allocations Across Gibbs Rounds
* **Location:** `src/ensemble_move.cpp`, line 421 (`emc_ensemble_rounds_cpp`).
* **The Defect:**
  ```cpp
  const arma::mat Z = arma::solve(arma::trimatl(Lg), stack_t.each_col() - tmu);
  ```
  `stack_t` is dimension $p \times (n \cdot M)$. For `lba7r`, this is $20 \times 3392 = 67,840$ doubles ($542\text{ KB}$). Across $R = 10$ rounds, this line allocates and frees two $542\text{ KB}$ matrices **20 times per iteration** (over $10\text{ MB}$ of transient heap churn per iteration).
* **Remedy:** Pre-allocate two scratch matrices (`diff` and `Z`) once outside the round loop and pass them by reference.

### 4.4 Unbatched Likelihood Evaluations on Workers
* **Location:** `R/chain_pool.R`, lines 1511–1515.
* **The Defect:** Workers evaluate particle candidates in `safe_new_particle`, and then immediately make a second, separate call to `.emc_ensemble_ll_candidate` for the $M-1$ pool rows.
* **Performance Impact:** `calc_ll_manager` carries fixed per-call overhead (parsing data, unpackaging parameters, setting up subset indexes). Evaluating $100$ particles and $31$ pool rows in two separate calls incurs this overhead twice.
* **Remedy:** Concatenate candidate matrices (`rbind(particles, pool_rows)`) and invoke `calc_ll_manager` in a single batched call per worker.

### 4.5 Retaining Gumbel-Max Selection Artifact in C++
* **Location:** `src/ensemble_move.cpp`, lines 446–450.
* **The Defect:** In each round, selection simulates categorical draws using the Gumbel-max trick:
  ```cpp
  const double g = -std::log(-std::log(R::runif(0.0, 1.0)));
  ```
  Gumbel-max was originally introduced in R solely because `stats::runif(n * M)` and `max.col` vectorizes in R. In C++, generating a categorical draw from log-weights requires only **1 uniform draw per subject** using standard cumulative sum inversion or alias sampling.
* **Performance Impact:** For $n = 106, M = 32, R = 10$, Gumbel-max draws **33,920 uniforms** and evaluates **67,840 logarithms** per iteration, compared to 1,060 uniforms and zero logarithms under standard inversion. While necessary for seed verification during Phase 1, keeping Gumbel-max permanently in production is an unneeded CPU penalty.

### 4.6 Recomputing Scalar Invariants in `log_q_m2`
* **Location:** `src/ensemble_move.cpp`, lines 27–34 (`log_q_m2`).
* **The Defect:** For every pool member of every subject, `log_q_m2` recomputes:
  ```cpp
  R::lgammafn((nu + p) / 2.0) - R::lgammafn(nu / 2.0) - 0.5 * (p * std::log(nu * M_PI) + p * std::log(c))
  ```
  These terms depend solely on constants $p, \nu, c, w$. On `lba7r`, `log_q_m2` is called $3,286$ times per iteration, executing `lgammafn` **6,572 times per iteration** on identical values.
* **Remedy:** Hoist the scalar constants into a precomputed struct.

### 4.7 Solving Full Dense Triangular Systems for Block-Diagonal Covariances
* **Location:** `src/ensemble_move.cpp`, line 421.
* **The Defect:** In standard hierarchies where parameters are partitioned into multiple small covariance blocks and independent unblocked variances, $L_g$ is block-diagonal. Calling `arma::solve(trimatl(Lg), ...)` solves a full $p \times p$ triangular matrix for all $n \cdot M$ columns, doing $O(p^2)$ work per candidate rather than $O(\sum d_k^2)$.
* **Remedy:** Exploit `BlockPlan` to perform triangular solves block-by-block.

---

## 5. Strategic Recommendations for the Sampler

### 5.1 Decouple $M$ and $R$ into Independent, Robust Rules
The empirical findings in WS5 demonstrated that $M$ and $R$ play orthogonal roles:
* **Fix Rounds Structurally at $R = 3$:**
  * Round 1 decorrelates the means.
  * Rounds 2–3 capture over $95\%$ of available variance displacement.
  * Rounds $4\dots 10$ provide no incremental gain.
  * Because rounds cost only arithmetic ($< 0.1\text{ ms}$), attempting to adapt $R \in [2, 64]$ dynamically via noisy numerical gradients adds algorithmic instability for negligible performance return. Setting $R = 3$ statically cuts round arithmetic by $70\%$ safely.

* **Size Pool Size $M$ via Subject Effective Sample Size ($\text{wESS}$):**
  * Avoid optimizing $M$ using single-iteration additive displacement heuristics.
  * Size $M$ based on importance sampling health: ensure that even the least efficient subject achieves a minimum target effective sample size ($\text{wESS} \ge \tau$, e.g., $\tau \approx 3\text{--}5$):
    $$M_{\text{target}} = \text{clamp}\left( \left\lceil \frac{\tau}{\min_s (\text{wESS}_s / M)} \right\rceil, \; M_{\min} = 8, \; M_{\max} = 64 \right)$$
  * Imposing a hard floor of $M_{\min} = 8$ prevents the posterior support collapse observed on Forstmann, while capping at $M_{\max} = 64$ prevents run-away worker CPU on difficult geometries.

### 5.2 Implementation Priority Roadmap
1. **Immediate (C++ Efficiency Polish):**
   * Precompute $L_s$ and $\log|V_s|$ in `add_proposals()` to eliminate redundant Cholesky calls in `draw_pools_cpp`.
   * Eliminate the double-inversion in `riwish` / `gibbs_standard`.
   * Hoist invariant scalar terms in `log_q_m2`.
   * Pre-allocate scratch buffers for Gibbs rounds.
2. **Intermediate (Phase 2 Stabilization):**
   * Replace the 2D $(M, R)$ grid search with static $R = 3$ and $\text{wESS}$-guided $M$ adaptation with a strict floor ($M \ge 8$).
   * Run the queued Forstmann ablation and validate on `nback31` and `lba7r`.
3. **Architectural (Worker Batching):**
   * Batch pool candidate likelihood evaluations with particle proposals in a single call to `calc_ll_manager` on workers.
   * Port standard models with group designs (`sampler$group_designs`) to C++.
