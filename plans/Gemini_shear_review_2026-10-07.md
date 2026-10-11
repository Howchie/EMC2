# EMC2 Hierarchical Sampler Review: Shear Group Move Investigation

**Author**: Gemini  
**Date**: 2026-10-07  
**Branch**: `playground` (commit `11993e52`)  
**Package Path**: Verified loaded from isolated temporary library `/home/ubuntu/.gemini/antigravity-cli/brain/6d4ae4ef-9cc7-457e-ba63-4250e7f12efc/scratch/Rlib/EMC2`. Working tree clean and uncommitted.

---

## 1. Overhead: Where Does the Shear's Extra Cost Go?

We instrumented and profiled one chain's iterations across the shear components on `freeA15` ($n = 30$ subjects, $p = 14$ parameters):

| Component | Profiled Cost per Iteration (CPU) | % of Shear Cost | Analysis & Mechanism |
|---|---:|---:|---|
| **Planning** (`.emc_shear_plan`) | **0.00002 s** | < 0.01% | 0.0023 s per call; only runs at proposal rebuilds (iters 100, 200, 400...), negligible per-iteration cost |
| **Delayed-Acceptance Screen** | **0.0025 s** | 0.8% | Quadratic form delta `.emc_shear_lik_delta` is cheap when active |
| **Cholesky & Frames** (`.emc_shear_propose` / `state`) | **0.0343 s** | 10.7% | Recomputes $\text{chol}(V)$ and 30 subject Choleskys $\text{chol}(H_s + V^{-1})$ for every candidate step |
| **Worker Likelihood Rounds** (`.emc_wpool_subject_ll`) | **0.2839 s** | **88.2%** | **The dominant bottleneck**: evaluated 28 full 30-subject likelihood rounds per iteration |
| **Total Shear Overhead** | **0.3218 s / iter** | 100% | Multiplies iteration cost by **1.64x** over base sampler (453s vs 884s total CPU) |

### Key Findings on Overhead

1. **Planning is negligible**: Planning costs 0.0023 s per call and runs only at rebuilds, contributing virtually 0% to per-iteration overhead.
2. **Likelihood evaluations dominate (88.2%)**:
   - The shear move loop (`R/shear_move.R:341`) tests up to $2p = 28$ candidate steps (14 mean directions + 14 scale parameters).
   - In `freeA15` seed 2, **delayed acceptance was completely disabled** because 3 subjects failed the Laplace check. Every proposed step bypassed the screen and triggered a synchronous worker-pool round trip evaluating true likelihoods for all 30 subjects (28 rounds/iter).
   - Even when delayed acceptance is forced active, 64.3% of steps pass the screen, still triggering **18 likelihood rounds/iter** (80.6% of shear CPU).
   - In parallel worker mode (e.g. 5 cores per chain), each round incurs IPC serialization and synchronization latency across all workers.
   - Contrast with `ensemble_move`: Ensemble move takes only **one** parallel worker pool dispatch per iteration (drawing the $M$ candidates), then runs $R$ group Gibbs rounds purely in hierarchy arithmetic with zero likelihood calls.
3. **Cholesky redundancy in mean steps**:
   - For all 14 mean directions ($u_{\text{sd}} = 0$), the group covariance $V$ is unchanged ($V' = V$).
   - Yet `.emc_shear_state` recalculates $\text{chol}(V)$, $V^{-1}$, and $\text{chol}(H_s + V^{-1})$ for all 30 subjects from scratch on every step: **434 redundant Cholesky factorizations per iteration**.

### What Can Be Cut Without Changing the Target Distribution?

1. **Mean-Step Cholesky Bypass (Lossless, Identical Invariance)**:
   When $u_{\text{sd}} = 0$, $V' = V \implies Q_s' = Q_s \implies R_s' = R_s$. The map simplifies to:
   $$a_s' = a_s + Q_s^{-1} V^{-1} u_\mu = a_s + R_s^{-1} (R_s^{-\top} V^{-1} u_\mu)$$
   $\log |\det J| = 0$ and $\log p(V') = \log p(V)$. Reusing the cached Cholesky factors $R_s$ cuts all 420 subject Cholesky decompositions and $V$ Choleskys per iteration down to a single backsolve per subject.
2. **Partial (Masked) Delayed Acceptance**:
   Line 270 of `R/shear_move.R` disables delayed acceptance for the whole chain if *any* subject lacks a usable fit (`delayed = !any(...)`). Under Christen & Fox (2005), DA is valid for *any* surrogate $\tilde{\pi}$. For subjects lacking a fit, setting $\tilde{L}_s \equiv 1$ ($\Delta \log \tilde{L}_s = 0$) screens the move on the remaining 27 subjects, filtering out bad group moves before touching the worker pool.
3. **Step Pruning by Eigenvalue**:
   Filtering out directions with negligible autoregression ($\lambda_i < 0.1$) prevents wasting full likelihood rounds on uninformative moves.

---

## 2. Inconsistency: Why Did Seed 2 Fail?

Seed 2 failed with $\hat{R}_{\alpha, \max} = 1.318$, $\hat{R}_{\mu, \max} = 1.316$, and $\text{ESS}_{\text{bulk}}$ of 10.6 on $A$ and 15.8 on $B$, whereas `ensemble` alone converged ($\hat{R} \le 1.08$, min ESS 24–56).

The root cause was an architectural **pipeline sequencing bug** in proposal construction that triggered a cascade of failures:

1. **The Sequencing Bug**:
   In `R/fitting.R:411-414`:
   ```r
   emc <- create_chain_proposals(emc, samples_idx = sample_window, ...)  # calls .emc_shear_set_plan()
   emc <- .emc_check_likelihood_fits(emc, n_cores)                       # attaches fit$check
   ```
   `create_chain_proposals()` calls `.emc_shear_set_plan()`, which calls `.emc_shear_fits()`. But `.emc_shear_fits()` requires `is.list(fit$check)`. Because `.emc_check_likelihood_fits()` has not run yet, `fit$check` is `NULL` for all 30 subjects!
2. **Degenerate Eigenvalues & Wild Directions**:
   With 30 `NULL` fits, `.emc_shear_plan()` sets $Q_s = V^{-1} \implies S = n V \implies M = \tilde{R}^{-\top} V^{-1} S V^{-1} \tilde{R}^{-1} \approx I$.
   - Checked on `warmup_freeA15_2.rds`: Stored eigenvalues $\lambda_i$ were `0.99992, 0.99989, 0.99984...` and $w_i = 1.0$ for all 14 parameters (exact match to `p_null`, difference $2.2 \times 10^{-16}$).
   - Because $1 - \lambda \approx 10^{-4}$, the columns of `dirs` were divided by $\sqrt{1 - \lambda} \approx 0.01$, yielding column norms of **1.0**, compared to the true posterior directions with norms of **0.02 – 0.12** ($25\times$ to $50\times$ over-scaled!).
   - Because $M \approx I$, the eigenvectors were arbitrary orthogonal rotations rather than the true $A-B-t_0$ ridge.
3. **Step-Size Collapse during Adaptation**:
   Proposing massive steps along arbitrary directions off the narrow posterior ridge resulted in near-zero acceptance. Robbins-Monro adaptation subtracted $0.44 / \sqrt{n}$ on every rejection, driving `dir_log_step` from $0$ down to **-4.53** ($\text{scale} = e^{-4.53} = 0.0107$).
4. **Frozen Paralysis in Production**:
   At production start, step sizes were frozen at $0.01$. At transition 100, when `set_plan` finally saw fits with `$check`, it updated `plan$dirs` to the true directions, but **did not reset `dir_log_step`**, and adaptation is disabled in `sample`. Proposing steps of size $0.01 \times 0.02 = 0.0002$ paralyzed the move.
5. **Chain Separation**:
   In `warmup_freeA15_2.rds`, Chain 2 ended at $A = -1.808, B = 0.114$, while Chains 1 and 3 were at $A = -1.196$ and $-1.244$. Because the shear move was paralyzed and the centred Gibbs step only moves in $O(1/\sqrt{n}) \approx 0.03$ increments, Chain 2 remained stranded in its tail throughout 2000 transitions.

---

## 3. The Fix: Curved Transport vs. Alternatives

### Evaluation of "Curved Transport" (`scripts/curve_path.R`)

We executed `curve_path.R` across 12 posterior states (4 states $\times$ 3 chains) on `fits/emc_freeA15_seed1_shear.rds`:
- **Linear Transport Width**: **0.231**
- **Curved Transport Width**: **0.403** (median ratio **1.27x**)
- **Gaussian Surrogate Width**: **1.000**

**Verdict: Curved transport is NOT the right fix.**
1. *Insufficient expansion*: An increase from $0.23$ to $0.40$ is still $2.5\times$ narrower than the required group step scale ($1.0$).
2. *Joint dimensionality curse*: The shear move updates all $n$ subjects simultaneously. The MH acceptance ratio is $\prod_{s=1}^n \frac{L_s(a_s')}{\tilde{L}_s(a_s')}$. Even if a quadratic polynomial fits the average curvature, individual subjects have unique data and non-linearities. As $n$ increases, the probability that all $n$ subjects land in acceptable likelihood regions drops exponentially.
3. *Extrapolation instability*: Fitting $n \times p = 420$ quadratic coefficients $\kappa_s(\mu)$ on past draws is prone to polynomial divergence if a chain ventures slightly outside the sampled window.
4. *Violates model agnosticism*: Cognitive models exhibit boundary cutoffs (e.g., $t_0 < \min(RT)$, $A > 0$) that quadratic polynomials cannot respect.

---

## Ranked Recommendations & Cheap Tests

### Recommendation 1: Fix Pipeline Ordering, Reset Step Sizes, and Add Cholesky Bypass in Shear Move (Rank 1)
- **Fix 1A (Ordering)**: In `R/fitting.R:414`, run `.emc_check_likelihood_fits()` *before* `create_chain_proposals()`, ensuring `.emc_shear_set_plan()` receives valid Laplace fits.
- **Fix 1B (Step-Size Reset)**: When `.emc_shear_set_plan()` updates directions, reset `dir_log_step` to 0 (or re-adapt), preventing frozen collapse.
- **Fix 1C (Partial DA)**: In `R/shear_move.R:270`, enable delayed acceptance whenever at least one subject has a valid fit; use $\Delta \log \tilde{L}_s = 0$ for unapproximated subjects.
- **Fix 1D (Mean-Step Bypass)**: For $u_{\text{sd}} = 0$, evaluate $a_s' = a_s + R_s^{-1} (R_s^{-\top} V^{-1} u_\mu)$, bypassing all 420 subject Cholesky decompositions.
- **Cheap Test to Confirm/Kill**:
  Run 50 warmup and 50 sample iterations on `freeA15` seed 2 with Fix 1A+1B in a temporary library.
  - *Confirm*: `plan$lam` slowest eigenvalue $\le 0.94$, `dir_log_step` remains $\ge -1.5$, acceptance rate is $30–50\%$.
  - *Kill*: If `dir_log_step` still collapses below $-3.5$, linear shear cannot overcome the curvature even with the correct plan.

### Recommendation 2: Past-Draw Reservoir in the Default Ensemble Group Move (Rank 2)
- **Rationale**: The default `R/ensemble_move.R` is fundamentally superior to shear because each subject selects its own state independently from pool $A_s$. It converged on seed 2 without tuning.
- **The Upgrade**: Currently, pool $A_s$ ($M = 32$) is drawn from a Gaussian + Student-$t$ mixture $q_s$. On curved ridges, Gaussian draws land off the banana and get zero weight. If $M/2$ pool candidates are sampled uniformly from the subject's *own past production draws* (which lie exactly on the curved manifold), the effective pool along the ridge increases by $4–10\times$ at zero extra likelihood cost.
- **Mathematical Validity**: Neal (2011, arXiv:1101.0387) proves that any pool distribution $q_s(a)$ that does not depend on the current group state $h$ is strictly invariant under the augmented target $\tilde{\pi}(h, A, k)$.
- **Cheap Test to Confirm/Kill**:
  On `freeA15` seed 2, modify `ensemble_q` to draw 16 candidates from past draws and 16 from $q_s$. Measure ESS bulk for $A$, $B$, $t_0$ over 300 sample iterations.
  - *Confirm*: ESS/iter on $A$ and $B$ increases by $\ge 2\times$ over default ensemble (from 56 to $>120$).
  - *Kill*: If reselection weights concentrate exclusively on the current state.

### Recommendation 3: Differential Evolution Tangent Proposals in Subject Step (Rank 3)
- **Rationale**: `frozen_group.R` showed conditional ESS/iter is $0.43$. However, within `new_particle()`, all proposals are Gaussian. In curved subject posteriors, adding a differential evolution vector $a_s^* = a_s + \gamma (a_{s, i} - a_{s, j}) + \epsilon$ drawn from past draws proposes directly along the local tangent of the banana without requiring any parametric model.
- **Cheap Test to Confirm/Kill**:
  Test on a 2D Rosenbrock banana toy model. Confirm ESS increases by $>50\%$ over Gaussian RW.

---

## 4. Toy Model Verification of Invariance

The proposed **Mean-Step Cholesky Bypass** and **Target Invariance** were mathematically verified on a toy hierarchical model with $N = 1000$ exact rejection-sampled posterior draws (`scratch/test_toy_bypass.R`):
- **Test 1 (Algebraic Exactness)**: Verified across 300 subject evaluations that $a_s' = a_s + Q_s^{-1} V^{-1} u_\mu$ matches the standard Cholesky transport $a_s' = m_s' + R_s'^{-1} R_s (a_s - m_s)$ with maximum difference $< 10^{-12}$ (300/300 passes).
- **Test 2 (Posterior Invariance)**: The toy shear move with mean bypass achieved 87.1% acceptance, and sample moments before and after the move matched within Monte Carlo error:
  $$\Delta \mu_1 = 0.0041 \le 4 \times \text{SE} = 0.048; \quad \Delta \mu_2 = -0.0079 \le 4 \times \text{SE} = 0.049 \quad (\text{PASS})$$
