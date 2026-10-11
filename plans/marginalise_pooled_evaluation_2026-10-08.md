# Evaluation of the Marginalise and Pooled Parameter Architectures in EMC2

**Date:** 2026-10-08  
**Scope:** Architecture evaluation of `marginalise` and `pooled` parameters, computational overhead breakdown, cross-sampler interactions, and missed optimization opportunities.

---

## 1. Executive Summary

Both the `marginalise` and `pooled` architectures in **EMC2** are designed to solve fundamental identification, mixing, and parameter-parsimony challenges in cognitive models by reducing the effective parameter dimension:
- **`marginalise`** integrates out a subject-level parameter (currently primarily non-decision time $t_0$ in race models) with respect to a supplied prior $p(t_0 \mid \eta)$. This removes $t_0$ from particle selection and may reduce coupling with drift or threshold parameters; it does not guarantee that all posterior dependence disappears.
- **`pooled`** enforces zero between-subject variance ($\sigma^2_{\text{between}} = 0$) on designated parameters $\phi$, collapsing $S$ individual subject parameters into a single global vector shared across all subjects.

The two architectures add different kinds of work, but source inspection alone does not establish their wall-clock or ESS-per-second impact:
1. **`marginalise`** evaluates a subject likelihood at multiple quadrature and mode-finding nodes. Under current defaults this is roughly 15--25 likelihood rows per particle on common paths, with additional rows for repairs. Rows are batched into C++ calls, so these figures are not function-call counts.
2. **`pooled`** adds an eight-candidate pooled-parameter move by default and computes a finite-difference curvature grid periodically during non-sample stages. The group update still operates on the other, non-pooled parameters; only pooled coefficients leave the ordinary group hierarchy.
3. **Feature gates:** Marginalisation and pooled coefficients skip several optional sampler features. Pooled sampling also has its own shear transformation inside `pooled_move`; that is distinct from the general `shear_move`. These gates are real, while their end-to-end warmup and sampling cost has not been measured here.

`marginalise` and `pooled` cannot currently be combined. Treat performance claims below as hypotheses until measured on representative models, trial counts, subject counts, and worker configurations.

---

## 2. Side-by-Side Architectural Comparison

| Dimension | `marginalise` Architecture | `pooled` Architecture |
| :--- | :--- | :--- |
| **Target Parameter** | One subject-level coordinate, currently primarily $t_0$, integrated for likelihood evaluation and reconstructed as a posterior draw | Coefficients marked by `design(pooled = ...)`, shared across subjects |
| **Statistical Assumption** | Uses the supplied fixed prior $p(t_0 \mid \eta)$; its hyperparameters are not estimated by this marginalised update | The pooled coefficient has zero between-subject variance and its own prior |
| **Parameter Reduction** | Removes the coordinate from likelihood-based particle selection; a conditional draw is still stored in output | Removes subject-specific copies and their group hyperparameters for that coefficient |
| **Sampling Mechanism** | Adaptive composite Gauss-Legendre rule per particle, then a conditional node draw for the selected particle | `pooled_move`: default eight-candidate conditional importance sampling with a shear transformation |
| **Likelihood rows / subject / iteration** | With default $n_{nodes}=12$ and $N_p=100$, common warm/cold paths span roughly 1,500--2,500 rows; predictive paths vary and repairs can add work | About 100 particle rows + 8 pooled candidate rows at $N_p=100$, excluding periodic curvature rows |
| **Curvature / tuning cost** | Mode scans and refinements depend on warm/prediction path and repairs | $1+2D^2$ rows per subject on a curvature refresh: 339 for $D=13$, 883 for $D=21$; refreshed during non-sample adaptation, not sampling |
| **Group-level update** | Pins the marginalised coordinate's group mean and variance to the supplied prior values | Pooled coefficients are outside the ordinary group covariance; other group parameters keep their existing updates |
| **Likely optimization target** | Reuse mappings for parameters unchanged across quadrature nodes; measure node evaluation and mapping costs | Reduce curvature work or allocation/transport costs if profiling shows they matter |
| **Mixing impact** | Integrating $t_0$ removes it from particle-selection coordinates and may reduce dependence with other parameters; it does not guarantee all remaining dependence disappears | The shear transformation targets coupling between pooled and subject/group states; acceptance and ESS depend on the workload |

---

## 3. Deep Dive: The `marginalise` Architecture

### 3.1 Mechanics & Cost Profile
In race models (LBA, RDM, RDMSWTN), response-time support constrains non-decision time $t_0$, which can be correlated with race speeds and thresholds. The conditional width and strength of this dependence vary with model and data and should be reported from representative workloads.
In [`src/marginal_likelihood.cpp`](file:///data/work/EMC2_dev_oo/src/marginal_likelihood.cpp), `calc_ll_oo_marginal_core` integrates $t_0$ over the feasible interval $[\log(\text{bound}), \log(\min(RT))]$:
1. **Pilot scan:** 7 points spanning the interval across all particles ($7 \times N_p$).
2. **Parabolic refinement:** 2 rounds of 3 points ($6 \times N_p$) to resolve the peak curvature.
3. **Quadrature:** 12 composite Gauss-Legendre nodes ($12 \times N_p$); up to two nodes are assigned to tail panels, depending on whether the fitted core panel leaves tail intervals.
4. **Verification & repair:** Retries under-resolved particles with full pilot + refine + eval ($25 \times N_{\text{repaired}}$).

Approximate likelihood rows per particle under the default settings:
- Cold start: $7 + 6 + 12 = 25$ evaluations per particle.
- Subset regression (`predict_mode = TRUE`): the row count depends on the number of varying predictors, subset size, and repair path. The current estimate of $\approx 19.6$ needs a derivation or measurement before use as a planning figure.
- Warm start (`warm_mode`, `warm_sd`): typically one three-point refinement round plus 12 quadrature rows, with more work when a hint needs another refinement or repair.

For $N_p = 100$ particles per subject, the cold path is 2,500 likelihood rows versus 100 rows for a standard likelihood step. This is a 25-fold row-count comparison for that path, not a wall-clock multiplier: quadrature rows are batched, and mapping, kernel, data size, repairs, and parallelism affect elapsed time. Report C++ invocation counts and measured wall time separately.

### 3.2 Potential Optimization Opportunities in `marginalise`

#### M1. Redundant Design Mapping & Nonlinear Transforms in C++
- **Location:** [`src/marginal_likelihood.cpp:56-97`](file:///data/work/EMC2_dev_oo/src/marginal_likelihood.cpp#L56-L97) (`marginal_eval_nodes`).
- **The Problem:** Quadrature nodes flatten into an expanded matrix `pm` of size $(N_p \cdot K) \times d$, which is passed to `calc_ll_oo` with `varying = R_NilValue`. In [`calc_ll_oo`](file:///data/work/EMC2_dev_oo/src/particle_ll.cpp#L637-L760), `ParamTable` re-applies design matrices (`m.table.apply_designs`) and evaluates nonlinear transforms (`m.table.apply_transforms`) across all trials for every single row.
- **The Inefficiency:** For a given particle, all parameters except $t_0$ are repeated across its quadrature-node rows. The current call has no varying mask, and the existing invariant-design optimization works across columns unchanged over the whole batch, not repeated particle-by-particle node blocks. This makes repeated design mapping and transformation a credible target. In supported race models, $t_0$ changes the response-time argument while other trial parameters may be reused.
- **Optimization:** Prototype a path that maps and transforms $t_0$-invariant parameters once per particle, then evaluates the race kernel over node-specific response-time arguments. Preserve each model's designs, transforms, bounds, trend handling, and response-time-bin behavior. The $3\times$--$5\times$ speedup is unverified; benchmark by model and trial count before retaining it.

#### M2. Frozen Sample-Stage Warm Start: Measure Before Changing
- **Location:** [`R/sampling.R:1789-1794`](file:///data/work/EMC2_dev_oo/R/sampling.R#L1789-L1794).
- **Observed behavior:** The code freezes both warm and prediction state in the sample stage:
  ```R
  if (!is.null(marginalise) && stage != "sample") {
    pm_settings[[i]]$marg_warm <- ...
    pm_settings[[i]]$marg_pred <- pred_next_state
  }
  ```
- **Potential cost:** If a frozen hint is repeatedly rejected during sampling, its backoff state is not retained and the same fallback work can recur. The frequency and cost of this case have not been measured.
- **Correctness constraint:** The finite quadrature rule is centered and shaped using the hint. The source treats the hint as part of the likelihood approximation and freezes it to keep the retained-stage approximate target fixed. Updating it during sampling can change the approximate likelihood; it is not established that the update is target-neutral.
- **Next step:** Instrument warm/prediction use, fallback, and repair rates during retained sampling. Consider a sample-stage change only after showing that it preserves the intended target or bounding numerical error to the project's criteria. Safer alternatives include choosing a robust fixed hint before sampling or designing a node rule whose result is independent of the hint to the required accuracy.

#### M3. Reuse Refinement Evaluations in Quadrature
- **Location:** [`src/marginal_likelihood.cpp:338-382, 600-639`](file:///data/work/EMC2_dev_oo/src/marginal_likelihood.cpp#L338-L382).
- **Observed work:** During parabolic refinement, the mode search evaluates the integrand at $(m - h, m, m + h)$ (up to 6 rows), then the quadrature evaluates its node locations separately.
- **Optimization hypothesis:** Prototype a nested rule that can reuse prior evaluations. The current refinement points are not nodes of the recentered composite Gauss-Legendre rule, so reuse requires changing the rule and its error checks. A low-order rule centered at the mode can miss tail mass or a narrow/skewed peak, especially near truncated support. Validate against a high-order reference over narrow, broad, skewed, boundary-adjacent, and repaired cases before estimating savings. The $5$--$8$ row saving is not yet established.

#### M4. Ensemble Candidate Warm Start: Not Applicable Yet
- The current ensemble support check rejects samplers with marginalised coordinates, so `extra_ll` does not score ensemble candidates for a marginalised sampler. The proposed cold-start penalty is unreachable on the current path. Revisit it only alongside support for marginalised coordinates in `ensemble_move`, with candidate-specific correctness and warm-start handling.

#### M5. R/C++ Matrix Allocation & Transfer Overhead
- **Location:** [`R/sampling.R:1768-1805`](file:///data/work/EMC2_dev_oo/R/sampling.R#L1768-L1805), [`src/marginal_likelihood.cpp:763-788`](file:///data/work/EMC2_dev_oo/src/marginal_likelihood.cpp#L763-L788).
- **The Problem:** `calc_ll_oo_marginal_nodes` returns full `nodes` and `log_terms` matrices although R uses the selected row to draw the stored $t_0$ and uses marginal likelihoods for particle selection. These matrices are R-managed objects created in C++; describing the boundary as a C++-to-R copy overstates the current transfer.
- **Optimization:** To avoid returning the full matrices, either return one conditional $t_0$ draw per particle row (then keep the draw for the row R selects) or make a second call after R selects the row. The first option changes where RNG is consumed and must preserve per-subject stream behavior and reproducibility. This can reduce returned/retained matrix memory; it does not by itself remove internal node and term storage needed for integration and repair.

---

## 4. Deep Dive: The `pooled` Architecture

### 4.1 Mechanics & Cost Profile
Pooled parameters $\phi$ are shared by all subjects. In [`R/pooled.R`](file:///data/work/EMC2_dev_oo/R/pooled.R):
- At the likelihood boundary, [`.emc_pooled_columns`](file:///data/work/EMC2_dev_oo/R/pooled.R#L40) appends $\phi$ to every subject's proposal matrix.
- After subject particle and ensemble group updates, [`.emc_pooled_move`](file:///data/work/EMC2_dev_oo/R/pooled.R#L316) samples $\phi$ using conditional importance sampling (CIS, $K = 8$ candidates) centered at $z \sim N(\phi, S)$.
- To prevent ridge-stalling between $\phi$ and subject parameters $\alpha_s$, a shift $\delta$ carries all $\alpha_s$ and $\mu$ along a shear slope $c$: $\alpha_s' = \alpha_s + c \delta$.
- Evaluating $K = 8$ candidates requires **$8 \times S$ likelihood rows** on the ordinary pooled move. The candidate count is configurable.
- Estimating $c$ via [`.emc_pooled_curvature`](file:///data/work/EMC2_dev_oo/R/pooled.R#L240) scores the finite-difference grid during non-sample adaptation, with a refresh every 20 adaptation steps. Sampling-stage curvature is frozen.

### 4.2 Potential Optimization Opportunities in `pooled`

#### P1. Finite-Difference Cost in `.emc_pooled_curvature`
- **Location:** [`R/pooled.R:213-235, 240-267`](file:///data/work/EMC2_dev_oo/R/pooled.R#L213-L235).
- **Measured from source:** `.emc_pooled_fd_rows` constructs a full central-difference grid of dimension $D = d + p$:
  $$\text{Grid size} = 1 + 2D + 4\binom{D}{2} = 1 + 2D^2$$
  - For $D = 13$ ($d=12,p=1$): **339 likelihood rows per subject**.
  - For $D = 21$: **883 likelihood rows per subject**.
  - With 40 subjects, one refresh scores 13,560--35,320 rows across subjects. The worker pool batches subject rows; this is not the same number of C++ calls or serialized messages.
  - Amortized over 20 non-sample adaptation steps, these correspond to about 17 or 44 curvature rows per subject per step. With 100 particle rows and 8 pooled-candidate rows, that is roughly 125 or 152 total rows per subject per adaptation step under an equal-row-cost assumption.
- **The Inefficiency:** The shear formula is
- **Optimization:**
  1. Prototype a diagonal approximation to $Q_{aa}$ using axial evaluations, then compare it with the full Hessian. `eff_var[[s]]` is a fitted conditional proposal covariance; it is not automatically a Laplace precision and needs a justified conversion before substitution.
  2. Compute cross-derivatives $Q_{ap}$ with forward differences as a lower-cost, lower-order alternative:
     $$\frac{\partial^2 f}{\partial \alpha_i \partial \phi} \approx \frac{f(x_0 + he_i + he_\phi) - f(x_0 + he_i) - f(x_0 + he_\phi) + f(x_0)}{h^2}$$
     This requires one additional row per $(\alpha_i, \phi)$ pair when the base and axial values are reused. Check truncation and cancellation error at parameter-appropriate step sizes. The resulting row formula also drops off-diagonal $\alpha$--$\alpha$ and $\phi$--$\phi$ curvature; the 39-row example is for $p=1$. For $p>1$, explicitly estimate or justify omitting pooled-coordinate cross terms.
  3. Total rows drop from $1 + 2(d+p)^2$ to $1 + 2(d+p) + d \cdot p$.
     For $d = 12, p = 1$: **drops from 339 to 39 rows per subject** (an 8.7-fold reduction in rows for that curvature refresh). Measure impact on total adaptation time and pooled-move mixing; do not present it as an end-to-end speedup before benchmarking.

#### P2. Pooled Candidate Row Construction and Transport
- **Location:** [`R/pooled.R:290-314, 344-348`](file:///data/work/EMC2_dev_oo/R/pooled.R#L290-L314).
- **Observed work:** The master builds a list of $S$ distinct matrices:
  ```R
  rows <- lapply(seq_len(n), function(s) {
    cbind(sweep(shift, 2L, alpha[, s], "+"), cands)
  })
  ```
  and passes these subject-specific rows through the worker-pool input path. In queue mode, shared inputs are serialized once to a shared file; without queue mode, workers receive their assigned subject subsets. The current design does not send the full list separately to every worker.
- **The Inefficiency:** `shift` ($K \times d$) and `cands` ($K \times p$) are **identical for every subject**. The worker already receives subject $s$'s current state `alpha[, s]`.
- **Optimization hypothesis:** Transmit the common `shift` and `cands` payload and compose each subject's rows where it is scored. This could reduce master-side matrix construction and shared-input size, but still creates candidate rows for evaluation. Use serialization and worker timing data to decide whether the change is worthwhile.

#### P3. Proposal Matrix Allocations (`.emc_pooled_columns`)
- **Location:** [`R/pooled.R:40-56`](file:///data/work/EMC2_dev_oo/R/pooled.R#L40-L56), [`R/sampling.R:2171`](file:///data/work/EMC2_dev_oo/R/sampling.R#L2171).
- **Observed work:** At likelihood calls where pooled columns are absent, `.emc_pooled_columns` creates a repeated matrix `add` and calls `cbind(proposals, add)`.
- **Cost to measure:** This allocation happens at the likelihood boundary when pooled columns are absent. Its impact depends on how often the boundary is entered and on proposal size; the stated allocation rate has not been measured.
- **Optimization:** First check whether pooled columns can be added once at the outer particle-step boundary and reused. A direct fixed-column path in C++ is another option, but pooled coefficients still participate in design mapping and transforms. Treating them as ordinary `constants` requires a new interface contract; existing constant handling does not establish that substitution by itself.

#### P4. Candidate Count $K = 8$ vs One-Candidate Metropolis-Hastings
- **Location:** [`R/pooled.R:26-30, 324`](file:///data/work/EMC2_dev_oo/R/pooled.R#L26-L30).
- **Observed cost:** Eight CIS candidates cost $8 \times S$ likelihood rows on each ordinary pooled move.
- **Hypothesis:** A one-candidate random-walk Metropolis option could reduce pooled-move scoring from $8S$ to $S$ rows. That is an eight-fold reduction within the pooled move, not the whole iteration. At $N_p=100$, the stated ordinary total falls from 108 to 101 rows per subject (about a 6.5% reduction, before curvature rows and unequal row costs).
- **Evaluation:** Compare the current conditional-importance move and a correctly constructed symmetric random-walk proposal using ESS per second for pooled coefficients and correlated subject/group parameters. The current pooled adaptation target is 0.5; 0.44 is a candidate setting for a particular random-walk regime, not a universal target for this sheared joint move.

#### P5. Pooled Move Worker Round Trip
- **Location:** [`R/sampling.R:843-852, 982-995`](file:///data/work/EMC2_dev_oo/R/sampling.R#L843-L852).
- **The Problem:** The particle step and ensemble step are evaluated in one worker round trip (`.emc_wpool_iter(..., extra = ens_pools$rows)`). The pooled move then runs as a separate worker round trip immediately afterward.
- **Constraint and experiment:** Pooled candidates depend on the final subject state after particle and ensemble updates. Curvature refreshes also determine the pooled proposal covariance before candidates are generated. Removing a round trip therefore requires restructuring dependent updates or changing move cadence. Profile barrier time first; evaluate any cadence change as a sampler change with mixing and invariance checks.

---

## 5. Feature Compatibility and Potential Cost

A relevant compatibility cost is that the current gates skip several optional features for these sampler configurations:

```
Active Setting                   Disabled Feature                     Reference
---------------------------------------------------------------------------------------------
marginalise = 't0'      -->     burn_mode_start (mode init)          R/mode_start.R:550
marginalise = 't0'      -->     laplace_refresh                      R/mode_start.R:759
marginalise = 't0'      -->     subject_hmc (Hamiltonian MC)         R/subject_hmc.R:91
marginalise = 't0'      -->     shear_move (group shear)             R/shear_move.R:91
marginalise = 't0'      -->     ensemble_move (group ensemble)       R/ensemble_move.R:69

pooled = 't0'           -->     burn_mode_start (mode init)          R/mode_start.R:550
pooled = 't0'           -->     laplace_refresh                      R/mode_start.R:759
pooled = 't0'           -->     subject_hmc (Hamiltonian MC)         R/subject_hmc.R:93
pooled = 't0'           -->     shear_move (group shear)             R/shear_move.R:92
```

The two configurations are mutually exclusive in the current implementation. `pooled` still supports the ordinary ensemble group move; its own pooled-parameter shear is part of `pooled_move` and is distinct from `shear_move`.

1. **`burn_mode_start` & `laplace_refresh`:** These routines skip marginalised and pooled samplers. This can affect initialization and proposal refresh, but the number of extra warmup iterations is workload-dependent and has not been established.
2. **`subject_hmc`:** When enabled by option, HMC support rejects marginalised and pooled states. This removes that optional update path for these configurations.
3. **`shear_move`:** The general group shear is unavailable. Pooled coefficients retain their separate built-in shear; remaining group parameters use their supported group updates. Measure any warmup or ESS impact instead of inferring a fixed penalty.

---

## 6. Actionable Implementation Roadmap

### Phase 0: Establish Baselines and Correctness Criteria
1. Benchmark representative marginalised and pooled models over trial counts, subject counts, particle counts, and worker configurations. Record wall time and ESS per second for target parameters and correlated subject/group parameters.
2. Count marginalisation rows by warm, prediction, cold, and repair paths. Record C++ batch-call counts separately from row counts.
3. For pooled sampling, record particle-step time, pooled candidate time, curvature time, and worker serialization/dispatch time. Compare adaptation and retained-sampling stages separately.

### Phase 1: Prototype Measured Engineering Targets
1. **Hoist invariant mapping for marginalisation:** Prototype M1 and compare likelihood values and posterior draws against the existing path; retain it only if the measured gain is material.
2. **Reduce pooled curvature rows:** Prototype the P1 approximation behind a controlled option. Compare approximate shear and covariance with the full result across realistic states, then compare acceptance and ESS per second. Keep the full path as a fallback.
3. **Reduce pooled allocation/transport only if profiling supports it:** Prototype P2/P3 and report matrix allocation, shared serialization, and worker elapsed-time changes. Preserve design and transform behavior for pooled coefficients.

### Phase 2: Algorithm and Numerical Experiments
1. **Nested marginal quadrature:** Prototype M3 against a high-order reference and explicit integration error criteria, including truncated and repaired cases. Keep full-support tail handling.
2. **Sample-stage warm hint:** Pursue M2 only after demonstrating a fixed retained-stage target or establishing a suitable numerical error guarantee. A speed result alone is insufficient.
3. **One-candidate pooled move:** Evaluate P4 as an opt-in kernel against the current CIS move, using ESS per second and checking pooled plus subject/group marginals.
4. **Worker barrier reduction:** Pursue P5 only if profiling identifies synchronization as material. Preserve the current update ordering's target, or validate the revised transition sequence.

### Phase 3: Expand Feature Support
1. **Enable mode start and Laplace refresh for pooled models:** Define optimization over the joint subject/pooled state, then validate starts and refreshed metrics against the current initialization path.
2. **Support marginalised coordinates in `ensemble_move`:** Add support for the non-marginalised coordinates, with conditional likelihood scoring and tests that preserve the marginalised target. This also makes M4 applicable for future evaluation.
