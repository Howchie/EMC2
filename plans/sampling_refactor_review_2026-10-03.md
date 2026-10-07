# Sampling refactor review

**Review date:** 2026-10-03  
**Checkout:** `playground`, HEAD `9aaedb42557c40e3cb71a1e884cb277a61bbff54`, including the existing uncommitted worktree changes.  
**Scope:** particle conditional importance sampling; the historical group proposal/transport moves; the replacement ensemble update; optimizer seeded starts and Hessians; proposal means, covariance estimates and adaptation; worker execution, checkpoint continuation and the current pooled-parameter integration. This is a review, with no sampler implementation changes.

## 1. Assessment

The refactor has a defensible core. The repaired particle kernel and the frozen-proposal ensemble construction have explicit augmented-target arguments. Batching likelihoods, reusing Mahalanobis distances, removing likelihood weighting from posterior moments, and preserving chain-local checkpoint state are substantial improvements.

The full implementation does **not yet justify an unconditional claim of correct posterior sampling across all supported hierarchies and production configurations**. The principal obstacles are:

1. Incorrect full conditionals in inherited factor/SEM/infinite-factor hierarchy code that the new ensemble repeatedly invokes.
2. An unsupported diminishing-adaptation claim for production refreshes, particularly the permanently bounded, randomly refitted efficient proposals.
3. SEM reporting and bridge densities that disagree with the subject-conditional model used by the repaired sampler.
4. A pooled likelihood evaluator that silently converts arbitrary failures into zero likelihood.

Other confirmed issues affect backend reproducibility, proposal quality and diagnostics. They should be addressed, but a poor proposal covariance is not itself posterior bias when its actual density is used consistently.

Passing tests are useful evidence about the configurations exercised. They do not resolve the algebraic defects below or establish adaptive ergodicity.

### Priority and evidence key

- **P1:** resolve before claiming general correctness or releasing the affected functionality as validated.
- **P2:** required robustness, reproducibility, diagnostic or efficiency improvement.
- **P3:** useful hardening or measured optimization.
- **Confirmed:** reproduced behavior or direct full-conditional algebra.
- **Proof gap:** a correctness assertion lacks the required argument; this does not establish bias in every application.
- **Heuristic:** potentially useful, but the stated justification or calibration is incomplete.

| ID | Priority | Finding | Evidence / affected configurations |
|---|---|---|---|
| F1 | P1 | Free factor/SEM regressions omit fixed contributions | Confirmed; mixed free/fixed nonzero coefficients |
| F2 | P1 | Infinite-factor shrinkage conditionals use stale cumulative products | Confirmed algebra; at least three factors |
| F3 | P1 | SEM cycles require a determinant term absent from the conjugate update | Confirmed algebra; accepted cyclic structural specifications |
| F4 | P1 | Doubling refresh times do not establish diminishing adaptation | Proof gap; retained production with capped efficient refits |
| F5 | P1 | SEM bridge target and reported moments are inconsistent | Confirmed; nonzero covariate effects |
| F6 | P1 | Pooled likelihood failures are silently swallowed | Confirmed code path; programming/infrastructure failures |
| F7 | P2 | SEM efficient proposals differ between serial and direct worker paths | Reproduced; nonzero covariate effects and fitted regression |
| F8 | P2 | Laplace fallback misses infinite/failed diagnostics and can retain stale results | Reproduced gate; failed/poor Gaussian approximation |
| F9 | P2 | General covariance regularization depends strongly on parameter units | Reproduced; nearly singular, differently scaled coordinates |
| F10 | P2 | Proposal stabilization can leak into the actual prior density | Confirmed code path; invalid/near-singular group covariance |
| F11 | P2 | Independent-chain option still shares block construction; convergence defaults are weak | Confirmed code; multiple chains and subblocks |
| F12 | P2/P3 | Checkpoint and sample-store contracts need broader continuation coverage | Current fixes are sound; remaining helper fragility |
| F13 | P2/P3 | Student-t scaling explanation and local-move calibration overstate their theory | Confirmed distributional calculation; efficiency claims |

## 2. Scope, history and review method

The changes referred to as `group.proposals` were traced through the internal historical setting `emc2.group_move_proposals` and the group-transport implementation. The current default is the ensemble method; the transport implementation has been removed.

| Revision | Relevant change | Review conclusion |
|---|---|---|
| `1fd2e30f` | Particle invariance: auxiliary centers and multinomial allocation | Correct direction; retained-state and proposal-density arguments matter |
| `3447e715` | Separate independent and local conditional-IS updates | Valid composition under fixed settings |
| `a994a8a7` | Noncontiguous parameter groups and multiple blocks | Important indexing repair; correlated-prior tests remain necessary |
| `aac54dcb` | Generalized group transport with multiple proposals | Auxiliary-variable interpretation is defensible; no basis to condemn it just for using multiple candidates |
| `ac3e4d36` | Ensemble default; latent-score ordering and SEM subject-prior repairs | Ensemble structure is valid subject to correct hierarchy updates |
| `4dd79e2d` | Ensemble likelihoods share the subject worker trip | Appropriate amortization; preserve proposal/RNG contracts |
| `05da1855` | Remove legacy, AM and RAM transport moves | Reduces maintenance and theoretical surface area |
| `abe2c850` | Gumbel selection and ensemble cadence | Exact categorical selection; cadence itself is not an invariance problem |
| `4a446e24` | Optimizer/Hessian seeded burn starts | Useful initialization; local approximations require diagnostic fallback |
| `9aaedb42` | Batched Hessians, joint blocks, Laplace diagnostics and integration fixes | Substantial engineering improvement with the limitations below |
| Current worktree | Posterior moments, covariance shrinkage, dynamic efficient means, step tuning, checkpoint preservation and pooled parameters | Reviewed directly; several important residual issues |

The existing Graphify graph was used for initial relationships. It predates the newest changes, so conclusions use direct source inspection and historical `git show`, rather than assuming the graph is current. The main inspected paths were `R/sampling.R`, `R/ensemble_move.R`, `R/mode_start.R`, `R/fitting.R`, `R/statistics.R`, `R/chain_pool.R`, `R/objects.R`, the hierarchy variants and the current `R/pooled.R`, together with their tests and convergence/benchmark notes. Adjacent SEM prior, bridge and summary paths were inspected where they define or report the same target.

This covers the current worktree fixes and their interactions with the sampling history. It is not a certification of every unrelated edit in every linked worktree, every likelihood implementation, or every model in the package.

## 3. Mathematical assessment of the transitions

### 3.1 Particle conditional importance sampling

For fixed hierarchy state, a subject has target density proportional to `L(a) g(a | h)`. Retaining its current value alongside independent fresh draws from a fixed mixture `q`, then selecting with weights `L g / q`, is a conditional importance sampling transition. Plain multinomial allocation of component counts gives the appropriate mixture bag; imposing a minimum count per component generally changes the bag distribution and invalidates that simple denominator.

For a local component, introducing an auxiliary center `z ~ N(a, S)` gives extended target `pi(a) N(z | a, S)`. Conditional on `z`, drawing candidates from `N(z, S)` makes the Gaussian auxiliary factor cancel against the proposal. Running each local component as its own such transition, composed with the independent-mixture transition, is defensible. The original state-centered mixture should not be reinstated without its own correct augmented construction.

Sequential subblock updates must evaluate the full correlated subject prior, keep other coordinates at their latest values and update the corresponding shared likelihood contribution. The current code explicitly takes this approach. It deserves joint-moment and block-order tests; checking a single coordinate's mean and variance cannot validate all cross-block behavior.

### 3.2 Historical group proposals

The transport held standardized subject residuals fixed while changing group mean/covariance. For Gaussian subject priors, the transport Jacobian cancels the relevant change in subject-prior density; the transformed group covariance still requires its own Jacobian. A single symmetric proposal can therefore have an ordinary Metropolis acceptance ratio. The multiple-candidate auxiliary-center construction can instead use categorical selection on the transformed target. It is not invalid merely because it lacks the reverse candidate machinery of an unrelated multiple-try Metropolis algorithm.

The historical proposal geometry is less compelling than the invariance argument. Subtracting a covariance term from empirical group covariance is not a general exact law-of-total-variance identity in nonlinear, correlated hierarchies. Adaptive Metropolis or RAM can learn marginal group geometry while the residual-preserving move needs geometry conditional on those residuals. Small conjugate examples do not establish high-dimensional effectiveness. With these methods removed, the practical action is to retain an archived derivation and benchmark record, not reintroduce them without a demonstrated advantage.

The archived AM/RAM implementation deserves credit for its explicit discarded adaptation, frozen-settling/production states, configuration signatures and refusal to initialize adaptation inside an existing production segment. Those contracts are useful precedents for F4 and F12. RAM's acceptance-based geometry is more directly associated with the transported conditional target than AM's marginal empirical covariance. Its batched shadow directions, however, are evaluated from the same current state and then applied as successive factor updates; this is an adaptation heuristic rather than a literal sequence of ordinary RAM transitions at successively updated factors/states. Because it freezes before production, that distinction concerns adaptation quality and cost, rather than invalidating the frozen Metropolis transition. Include shadow-direction likelihood cost and the discarded settling period in any historical performance comparison.

### 3.3 Ensemble update

Write the hierarchical target as

```text
pi(h, a) ∝ p(h) ∏s Ls(as) gs(as | h).
```

For each subject, place the current state at a uniformly random position in a pool of size `M`, with all other members drawn independently from a frozen `qs`. The augmented target is

```text
pi_tilde(h, A, k) ∝ p(h) ∏s [Ls(As,ks) gs(As,ks | h) / M
                            × ∏j≠ks qs(As,j)].
```

Reselecting indices with weights `L g / q` and applying a hierarchy transition that preserves the conditional target given the selected subjects preserves this augmented target. A correct Gibbs **sweep** suffices; it need not be one exact joint draw of every group parameter. Fixed positive round counts and fixed cadence are compatible with invariance. The construction agrees with the ensemble framework in [Neal (2011)](https://arxiv.org/abs/1101.0387) and its hierarchical use in [Shestopaloff and Neal (2013)](https://arxiv.org/abs/1305.0320).

The code correctly carries the final selected subjects' likelihoods, uses the proposal density for the retained member, and uses exact Gumbel-max categorical selection. Generating fresh pools before the subject step is reasonable when proposals are frozen and the random draws are distinct. Combining likelihood evaluation into the worker trip is then an implementation optimization.

Two qualifications are essential: the hierarchy sweep must itself be correct (F1–F3), and proposal adaptation from retained history needs a separate adaptive-MCMC argument (F4). The augmented-target proof alone establishes neither.

### 3.4 Latent scores and partial collapse

The subject particle update integrates out latent scores. Redrawing scores before a subsequent group update conditions on them repairs the previously stale-score ordering. Conditional ensemble weights use scores as part of the state; marginal weights integrate them out and require refreshing them before they are next used. This ordering requirement follows the partially collapsed Gibbs framework of [van Dyk and Park (2008)](https://www.ma.imperial.ac.uk/~dvandyk/Research/08-jasa-pcg.pdf).

With marginal reselection, stored scores can be stale relative to the final selected subject state. That does not automatically invalidate the marginal draws of hierarchy and subjects when scores are refreshed before use. However, joint summaries involving the stored `(alpha, eta)` pair require a final conditional score refresh or an explicit limitation on that stored pair. Add a joint-latent check, particularly when the ensemble is disabled or marginal weights are selected.

### 3.5 Optimizer and Hessian seeding

Initialization need not itself preserve the posterior when subsequent warmup is discarded. BFGS, batched finite differences, support-aware retries and prior-whitened curvature regularization are reasonable ways to construct starts and proposals. Treating joint likelihood blocks as separable is correct only for likelihood curvature; a correlated prior can still contribute off-block precision. The current approach retains the prior contribution.

If conditional precision is `P + H_L`, widening the prior fourfold gives local proposal precision `P/4 + H_L`. The associated Newton shift from conditional mode `m` is approximately

```text
(P/4 + H_L)^(-1) (3/4) P (m - prior_mean).
```

That is exact for the local quadratic approximation, not necessarily the true widened-prior mode of a nonlinear likelihood. Naming and diagnostics should reflect this distinction. A finite regularized Hessian or successful optimizer status does not establish a global mode, tail coverage or covariance accuracy.

The mode-sharing threshold based on averaged approximate log mass is a warmup heuristic. Reusing donor curvature under different group priors is not an exact evidence comparison. Keep independently dispersed starts in validation: synchronized starts can give reassuring chain comparisons while all chains miss another posterior mode.

## 4. Required changes and their acceptance criteria

### F1. Subtract fixed contributions in free-coefficient regressions — P1, confirmed

**Locations:** `R/variant_factor.R:240–248`; `R/variant_SEM.R:381–392` and `413–455`. These defects are present in the historical hierarchy code as well; the ensemble does not create them, but repeats the affected transitions.

For `y = Z_free beta_free + Z_fixed beta_fixed + error`, the free-coefficient conditional uses response `y - Z_fixed beta_fixed`. The factor loading update instead regresses `alphatilde[, j]` directly on free score columns. The SEM loading/K update similarly uses `ytilde[, j]` without subtracting fixed loadings/K contributions.

This is relevant to ordinary identification constraints: a factor row can have an anchor loading fixed to one and another loading free. With correlated score columns, the omitted fixed term is absorbed into the free coefficient.

**Reproduction:** holding both score columns at `(-2, -1, 0, 1, 2)`, setting row 2's fixed second loading to one, free first loading to zero, and setting `y2 = eta2`, the correct free-loading conditional mean is zero. A controlled call to the actual factor Gibbs function requests a mean of **0.9477505**. The probe controls score/mu draws and observes the conditional mean, rather than claiming to establish full-chain bias numerically.

For SEM B/G rows with correlated disturbances, the conditional response must also subtract the row's fixed structural contribution. After updating a row, its residual must use **all** current coefficients. The present refresh `eta[, p] - Z_p %*% coef_sample` drops fixed coefficients and corrupts the residual used by later rows.

**Required change:** construct full row predictors, subtract fixed contributions before each free regression, and recompute residuals from the complete current B/G row. Explicitly align heterogeneous prior variances by coefficient: `matrix(prior$lambda_var, n_pars, n_factors)` does not generally repeat a factor-specific variance across each column as intended.

**Acceptance:** analytic conditional mean/covariance tests with two or more factors, correlated scores, nonzero fixed coefficients and heterogeneous variances; then independent posterior validation. A one-factor example or a comparison against the same faulty Gibbs implementation is insufficient.

### F2. Refresh cumulative shrinkage after each delta draw — P1, confirmed algebra

**Location:** `R/variant_infnt_factor.R:244–255`; inherited.

The full conditional for each `delta[h]` depends on the current cumulative products `tau = cumprod(delta)`. The code refreshes `tau` after drawing `delta[1]`, but not after drawing each later delta. With three or more factors, the next rate therefore uses an earlier delta value.

**Example:** with deltas `(1, 2, 3)`, if the second draw becomes four, the third cumulative product must become twelve. The loop still uses six. Its `1/delta[3]` factor does not remove dependence on the newly sampled second delta.

**Required change:** recompute `cumprod(delta)` after each delta update, or update the affected suffix exactly.

**Acceptance:** deterministic conditional-rate tests with at least three factors and forced distinct delta draws, followed by independent posterior checks. The current two-factor invariance fixture cannot exercise this error.

### F3. Restrict SEM to an acyclic structure or include its determinant — P1, confirmed algebra

**Locations:** `R/variant_SEM.R:403–459`; specification handling at `1158–1169` and the lower-level SEM settings path.

The structural model is `(I - B) eta = G x + disturbance`. Its score density includes `|det(I - B)|^n_subjects`. The Gaussian B/G regression update omits this term. For an acyclic structure, after topological ordering the determinant is one, so that omission is harmless. Rejecting self-loops alone does not establish acyclicity.

A two-factor cycle with both `B12` and `B21` free has determinant `1 - B12 B21`, which depends on the parameters being sampled. The accepted specification therefore need not have the conjugate B conditional used by the code.

**Required change:** validate a DAG across both high-level and direct settings, including nonzero fixed edges, and reject cycles; alternatively implement a correct nonconjugate B transition including determinant and invertibility constraints. State the supported model class clearly.

**Acceptance:** cycle rejection tests, topologically permuted DAG equivalence, and an independent SEM posterior reference for a model with nonzero structural and covariate effects. If cycles are supported, test the determinant contribution directly.

### F4. Establish production adaptation correctness — P1, proof gap

**Locations:** `R/fitting.R:331–409`, especially the claim at `334–340`; `R/objects.R:120–135`; efficient proposal fitting.

The comment states that changing proposals at `100, 200, 400, ...` iterations makes adaptation diminish. Sparse **deterministic** adaptation times do not themselves imply diminishing adaptation: changes can remain large along that infinite subsequence. The standard condition concerns the difference between successive transition kernels, not the density of refresh times. A separate containment/control condition is also needed in the usual framework. See [Roberts and Rosenthal (2007)](https://www2.stat.duke.edu/~scs/Courses/Stat376/Papers/AdaptiveMC/RobRosCouplingErgod2007.pdf).

The chain/ensemble expanding-window moments may become stable under suitable assumptions. The efficient proposal uses a fresh random subset capped at 750 draws, however, so its fitting noise need not vanish even with unlimited history. Smaller capped conditioning subsets can add noise. Reblocking introduces additional discrete changes.

**Reproduction:** repeatedly extracting/refitting efficient proposals from the same 10,000-draw Gaussian history gives standard deviations **0.05081963** for an estimated variance and **0.03395488** for a regression slope across 100 refits. This shows a permanent random fitting component; it does not measure a kernel total-variation distance or demonstrate posterior bias on that fixture.

**Required change:** the simplest defensible default is to freeze all proposal parameters, mixture/count settings, conditioning regressions and block partitions after discarded adaptation, then retain production from that fixed transition. If continued production adaptation is important, increase the fitting sample budget, damp changes with an explicit vanishing gain and provide the relevant joint-chain ergodicity argument and covariance/scale controls. Another justified adaptive framework could be used, but must be specified and proved rather than inferred from doubling times.

Track proposal epochs using a monotone transition counter, independent of trimming or thinning the stored history. Remove the current theorem claim until the implementation meets its assumptions. Finite early adaptation followed by freezing is compatible with this recommendation.

**Acceptance:** a default frozen-production contract with assertions covering every mutable proposal field, or a documented theorem-to-code mapping. Retain the Gaussian refresh experiment as a smoke check; it cannot replace that argument.

### F5. Unify SEM subject densities, bridge densities and reported moments — P1, confirmed

**Locations:** `R/variant_SEM.R:486–503`, `563–608`, `860–881`, and prior-derived covariance at `251–253`.

With `A = Lambda (I - B)^(-1)` and `C = K + A G`, the subject-conditional model is

```text
mean_s = mu + C x_s
V_cond = A Delta^(-1) A' + E.
```

The repaired particle/ensemble paths use this conditional model. The bridge group density instead evaluates all subjects with one population mean/covariance averaged over covariates. For observed covariates, that is a different density. Group information criteria likewise use common population moments rather than each subject's conditional density.

Population covariance, when that summary is wanted, is

```text
V_pop = V_cond + C Var(x) C'.
```

The current formula separately adds `A G Var(x) G' A'` and `K Var(x) K'`, omitting the two cross terms. With scalar effects `K = 0.3`, `A G = 0.8`, `Var(x) = 1`, the code's covariate contribution is **0.73**, versus **1.21**. This affects summaries and prior-derived quantities even if subject sampling uses the corrected conditional variance.

There is also a direct reporting error in `get_mu_implied()`: the `+ K %*% x_mu` expression starts a separate statement after the assignment, so it is evaluated and discarded. A two-parameter, two-draw probe with `mean(x) = 1`, `A G = 0.8`, `K = 0.3` returns **0.8** where **1.1** is required. The helper also drops dimensions in one-draw/one-parameter cases.

**Required change:** share one tested SEM conditional-moment helper across sampling, bridge likelihoods, prediction and subject-level information criteria. Compute population moments separately with all cross terms. Repair implied-mean assignment and use explicit `drop = FALSE`/matrix dimensions. Check variance ordering in SEM priors and bridge reconstruction for heterogeneous coefficient priors.

**Acceptance:** compare exact group log densities across all paths at the same parameter state and covariates; test nonzero K and G together, centered and uncentered covariates, one parameter/draw and correlated disturbance blocks. Validate marginal-likelihood results against an independent simple reference, not only finite-return checks.

### F6. Apply the normal likelihood failure policy to pooled moves — P1, confirmed

**Location:** current untracked `R/pooled.R:114–121`.

The evaluator catches **every** likelihood exception and returns `-Inf` for every candidate. It also converts every nonfinite result, including positive infinity, to `-Inf`, without checking the expected row count or classifying the failure. Consequently a programming error or broken likelihood backend can appear as ordinary out-of-support proposals and repeatedly retain the pooled state. Changes elsewhere in the hierarchy can prevent a whole-chain stall check from identifying this.

The auxiliary Gaussian pooled move itself is defensible with frozen covariance/shear: shifting all subjects and the group mean together preserves their hierarchical residuals, the shear has unit Jacobian, and the pooled/group-mean priors are included. The identified defect is the evaluator contract, not that transformation.

**Required change:** use the established checked likelihood policy, require exactly one result per input row, distinguish legitimate `-Inf` support values from NA/NaN/+Inf and exceptions, and retry candidates individually only for classified proposal-domain failures. Honor strict failure settings and expose pooled failure/stall counters. Abort programming and infrastructure errors.

**Acceptance:** pooled equivalents of the existing failure-policy tests: valid support failures, one bad row, malformed length, arbitrary exception, +Inf, NA, worker loss and strict/non-strict settings. Verify the same surviving candidates and density evaluations across backends.

### F7. Make SEM efficient conditioning identical across backends — P2, reproduced

**Locations:** `R/sampling.R:1277–1285`; `R/variant_SEM.R:550–560`; worker direct-state construction in `R/chain_pool.R`.

Efficient proposals fit a regression on the SEM group intercept `theta_mu`. The serial particle path advances that mean using `parameters$tmu`; the direct worker path advances it using the subject-specific prior mean, which includes covariate effects. Those are different regressors. Standard grouped designs already expose `conditional_mu`; SEM does not.

**Reproduction:** with fitted efficient mean `(1, 1)`, regression `B = 0.5 I`, conditioning origin zero, intercept `(1, 2)` and subject mean `(4, 6)`, identical RNG seeds yield a worker-minus-serial fresh-proposal shift of **(1.5, 2)**. The actual `new_particle()` paths were exercised while capturing likelihood input rows.

Each path can still be a valid fixed-proposal CIS transition because it uses its corresponding proposal density. This is a backend parity and fitting-specification defect, not direct evidence that both transitions target different posteriors.

**Required change:** either fit on the subject-conditional mean and store it as `conditional_mu`, or explicitly transmit the global intercept to direct workers. Use the same regressor in fitting and execution.

**Acceptance:** serial, split, queued and degraded-worker equality with nonzero K/G/B, a nonzero fitted regression and preserved subject RNG streams. The current standard-LNR parity fixture does not cover this.

### F8. Make Laplace diagnostics total, fresh and relevant to the actual proposal — P2, reproduced

**Locations:** `R/fitting.R:919–921`, `1013–1036`; `R/mode_start.R:425–448` and `675–716`.

The poor-approximation gate begins with `is.finite(pareto_k)`. This marks infinite k as not poor and also suppresses a high outside-support fraction when k is NA. For k `(0.2, Inf, NA, 0.2)` and outside fraction `(0, 0, 0.9, 0.9)`, the current result is **FALSE, FALSE, FALSE, TRUE**.

When fewer than two importance draws have finite weight, the diagnostic returns NULL. Refresh code updates the stored k/outside result only on success, so old successful diagnostics can remain attached to newly computed curvature. The reporting claim that k is refreshed for every subject is therefore too strong.

Moreover, a poor Laplace covariance can be replaced by an empirical covariance while the ensemble mean still comes from `burn_mode_mu`. That hybrid need not be appropriate, especially near a ridge or displaced mode. The diagnostic assesses the conditional Gaussian, whereas the executed burn proposal uses widened curvature and the ensemble may multiply its covariance by three.

**Required change:** return an explicit status for every attempt, including optimization/PSIS/support failure; treat +Inf k as poor and support failure independently of finite k. Replace stale results with current failure/unknown status. Choose mean and covariance together. Name the assessed Gaussian and, if using k to select the executed proposal, evaluate that proposal too or explain the proxy.

Finite-sample k estimates at 256 draws are noisy even if their population tail index is a property of the density pair. A fixed draw count should be evaluated empirically across dimensions and support geometries.

**Acceptance:** Inf/NA k, all/mostly-outside support, repeated success-then-failure refresh, empirical covariance fallback, and skewed/multimodal/nonlinear fixtures. Record optimizer stationarity, curvature regularization and actual proposal choice separately.

### F9. Regularize covariance in standardized coordinates — P2, reproduced

**Locations:** `R/statistics.R:580–607`; `R/fitting.R:824–846`, `939–981`, `1044–1046`. The standardized deterministic repair in `condMVN()` at `892–910` is a good existing model.

The general regularizer floors raw eigenvalues at `max_eigenvalue * 1e-6` or `1e-8`. For a singular matrix of ones, the regularized correlation is **0.999998**. Rescaling the second coordinate by one million, regularizing and transforming back gives **0.0009999985**. A change of units therefore almost removes the ridge. Raw fallback `0.5 I` has the same units problem.

**Required change:** use deterministic reference scales, regularize the resulting correlation/prior-whitened matrix and restore units. Specify prior-based scales for zero or unusable empirical variances. Avoid inferring those scales from the largest raw variance.

Correlation shrinkage preserves marginal scales and can improve finite-window proposals. Its temporal correction is a heuristic: the median lag-one AR(1) estimate of coordinate ESS is not generally the ESS of covariance products `x_i x_j`, and concatenating chains introduces nonphysical boundary lags. A slowly mixing direction can be concealed by the median. The independent-draw shrinkage derivation in [Schäfer and Strimmer (2005)](https://strimmerlab.github.io/publications/journals/shrinkcov2005.pdf) does not by itself establish this time-series extension.

**Acceptance:** unit-equivariance, singular ridge recovery, zero variance, one-dimensional cases, unequal chain means and slow covariance/product processes. Compare shrinkage choices by minimum relevant ESS per second, with chain boundaries handled explicitly. State the AR(1) correction as an approximation.

### F10. Separate proposal factorization from target-prior factorization — P2, confirmed code path

**Locations:** `R/sampling.R:1176–1191` and `1465–1489`.

`.chol_factor()` repairs failed covariance factorization with jitter and then a diagonal approximation. For a single parameter block, the resulting group proposal density `lp` is also used as `prior_density`. If the group covariance needed repair, the target-prior term silently becomes the repaired Gaussian rather than the specified prior. Multi-block handling has a separate full-prior path.

Normal well-conditioned group draws do not trigger this issue, and there is no evidence here of routine bias from it. The fallback nevertheless breaks the target contract exactly where correctness is most fragile.

**Required change:** distinguish strict target factors from regularized proposal factors. Validate actual group covariances once, reject/abort invalid target states with context, and never substitute a diagonal or jittered prior implicitly. Proposal regularization can remain, provided generation and denominator agree.

**Acceptance:** forced failed factorization and near-singular target tests that compare evaluated prior density against the stated model for one and multiple blocks. Test cache invalidation for changed covariance, partition and marginalization masks.

### F11. Complete independent-chain semantics and strengthen production diagnostics — P2

**Locations:** `R/fitting.R:653–665`, `769–799`, `360–365`, `1–25`, and `run_emc()`'s `rhat_version = "old"` default.

The independent-chain option builds own-chain empirical and efficient proposals and disables mode sharing. However, `sub_blocking(emc, n_blocks)` still averages covariance correlations across chains and assigns the shared result to every chain. This also occurs during production proposal refreshes. Thus the statement that these chains share nothing is false when subblocking is enabled.

Sharing a frozen pilot fit need not compromise the later target. Sharing adaptive production inputs couples the chain system and deserves a joint-system argument and an appropriate interpretation of diagnostics. Highly synchronized starts and pooled adaptation can weaken the practical evidence supplied by between-chain agreement.

Default production stopping checks alpha and mu with old R-hat and threshold 1.1. It omits important covariance, scale, loading, structural and pooled quantities. Mean burn R-hat can additionally hide a few seriously unconverged parameters. The rank-normalized split/folded R-hat and bulk/tail ESS approach is better supported for such diagnostics; see [Vehtari et al.](https://arxiv.org/abs/1903.08008).

**Required change:** freeze production partitions under F4; make partitions chain-local when independence is requested, or explicitly narrow that option's promise. Promote modern diagnostics with an approximately 1.01 production threshold, key hierarchy scales/correlations and latent/pooled quantities, plus bulk/tail ESS and task-relevant MCSE. Keep looser initialization criteria clearly separate.

**Acceptance:** independent mode with at least two chains and multiple subblocks; confirm own-chain proposals/partitions and distinct RNG streams. Include a fixture where alpha/mu pass but covariance or pooled effects remain poorly explored. Evaluate independently dispersed chains for multimodal starts.

### F12. Preserve and test the complete continuation state — P2/P3

**Locations:** `R/fitting.R:1577–1626`; static-field schema in `R/objects.R:1`; `R/ensemble_move.R:76–87`.

The checkpoint fix is well chosen: only identical immutable fields are shared, dynamic fields and proposal caches remain local, future unknown fields are preserved, and legacy checkpoints warn that historical adaptive state cannot be recovered. Retain this behavior.

Broaden regression coverage to multiple chains with different adaptive state and real production refresh boundaries. A single preburn continuation does not cover efficient regression attributes, ensemble proposals, proposal epochs, pooled shear, latent scores and updated RNG streams. Record effective sampler options/version with checkpoints: ensemble settings are currently read from process-global options during execution, so changed options can change continuation without a durable configuration record.

**Minor confirmed helper defect:** `.emc_ensemble_store()` guesses temporal arrays from their last dimension. When history length equals parameter dimension, it slices static `last_theta_var_inv`; a 3×3 matrix becomes 3×1 in the probe. **The existing standard fill immediately restores it to 3×3 before Gibbs, so this is not a demonstrated current ensemble crash.** Nevertheless the scratch helper violates the existing static-field contract and should use the central schema.

**Acceptance:** uninterrupted versus save/strip/restore continuation with two or three distinct chains, spanning the 100/200 refresh boundaries if those remain; compare samples, RNG states, proposal means/covariances, regression attributes, schedule counters and dynamic metadata. Exercise serial and queued execution, worker retirement/recovery, trimming/thinning, and pooled/latent configurations. Static matrices must survive extraction, scratch creation, rejection, concatenation and extension independently of history length.

### F13. Correct proposal-tail claims and calibrate local controllers — P2/P3

**Locations:** `R/ensemble_move.R:25–39`, `93–124`; local tuning in `R/sampling.R`; pooled adaptation in `R/pooled.R`.

The normal/t4 mixture generator and its density match. However, `c V` is the Student-t **scale matrix**, not its covariance. For four degrees of freedom its covariance is `2 c V`. If `M2` is squared Mahalanobis distance relative to V, then `M2/p` approaches the nondegenerate variable `4c / chi_square_4` as dimension grows. Its mean is `2c`, and its variance is infinite. It does not concentrate near `c p` as the current Gaussian-style typical-set explanation suggests.

Heavy tails can help stability and weight coverage. Tails lighter than the target do not automatically destroy invariance, and polynomial tails do not establish containment for a hierarchy whose means/scales move. Correct these comments and assess t scale/df/mixture weight on their actual radial distribution. Do not change the density merely to fit a comment.

The local-controller repairs—separate steps, current move statistics, and freezing production mixture/count/scale adaptation—are sensible. Targets inherited from 2.38 random-walk asymptotics are heuristics for a finite categorical CIS step, however. Move probability depends on particle count and the retained state. The same caution applies to a 0.234 pooled target with many candidates. Trace-based rescaling preserves an average step variance, not every direction after a covariance rotation.

**Acceptance:** radial-distribution tests and high-dimensional pool overlap measurements; controller comparisons using squared jump per likelihood cost, actual ESS per second and worst-subject behavior across particle counts. Keep covariance-layout resets deterministic and verify scales remain bounded during discarded adaptation.

## 5. Verification performed

### Environment and limits

Focused tests used R 4.6.1 with the current R source loaded through:

```r
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tests/testthat/helper-rng.R")
source("tests/testthat/helper-model-validation.R")
```

Development-only checks were enabled with `EMC2_TEST_LEVEL=full` and `NOT_CRAN=true` for the longer invariance tests. The pre-existing `src/EMC2.so` was used; no fresh native rebuild or package-wide check was performed. Results therefore validate the exercised worktree R paths with that binary, not a clean installed-package release.

### Existing focused tests

All executed expectations passed in these 16 test files:

| Area | Files |
|---|---|
| Proposal implementation | `test-proposal-density.R`, `test-proposal-cache.R`, `test-particle-adaptation.R` |
| Starts and scheduling | `test-burn-mode-start.R`, `test-sample-proposal-schedule.R` |
| Target smoke checks | `test-particle-invariance.R`, `test-ensemble-move.R`, `test-pooled.R` |
| Ensemble likelihood transport | `test-ensemble-pool-ll.R` |
| Failures and sampler guards | `test-sampler-guards.R`, `test-sampling-rejection.R`, `test-failure-policy.R` |
| Group/block integration | `test-par-groups.R`, `test-group-design-moments.R` |
| Backends and continuation | `test-audit-differential.R`, `test-worker-pool.R` |

The initially gated sample-schedule invariant check was subsequently run with full validation enabled and passed. The longer particle and ensemble invariant checks also completed successfully. A worker-pool test skipped the installed clean-process backend because that backend was not active. The rejection/worker fixtures emitted their expected failure-path warnings; these were not silently treated as release coverage for the skipped backend.

Important coverage limitations:

- The factor posterior invariance fixture uses one factor, so it does not expose F1's mixed anchor/free row.
- The infinite-factor fixture uses two factors, so it does not expose F2's stale product at the third update.
- The latent full-chain invariant loop covers factor and infinite-factor, not a complete SEM posterior run. The SEM check verifies conditional prior moments.
- Comparing ensemble results with the same hierarchy Gibbs implementation can preserve a shared defect in both reference and candidate.
- Standard-model backend parity does not establish SEM dynamic-regressor parity.
- Gaussian mean/variance agreement does not certify tails, multimodality, covariance products, irregular support or adaptive ergodicity.

### Independent probes

The companion [reproduction script](sampling_refactor_review_checks.R) records the focused demonstrations used here. Run it from the repository root with `Rscript plans/sampling_refactor_review_checks.R`. It uses scoped testthat mocks and prints observations; it does not change sampler source or form a new package test suite.

| Probe | Observed result | Interpretation |
|---|---|---|
| Factor row with fixed anchor and free loading | Requested free conditional mean 0.9477505; correct mean 0 | Incorrect conditional |
| SEM dynamic efficient mean | Worker fresh rows shifted by (1.5, 2) with identical RNG | Backend conditioning mismatch |
| Implied SEM mean | 0.8 returned; correct 1.1 | K contribution discarded |
| Covariate covariance contribution | 0.73 coded; correct 1.21 | Missing cross terms |
| Poor-Laplace gate | Inf k and NA/high-outside cases not classified poor | Failed fallback logic |
| Covariance under unit change | Correlation 0.999998 versus 0.0009999985 | Units-dependent regularization |
| Repeated efficient fits on fixed history | Nonzero SD for variance and slope | Bounded refit noise does not vanish |
| Ensemble static matrix scratch | 3×3 becomes 3×1, then fill repairs it | Helper contract defect; current repair acknowledged |

F2 and F3 are additionally established by their displayed conditional algebra, not by claiming a completed numerical posterior comparison.

## 6. Optimization opportunities, after correctness fixes

1. **Cache ensemble proposal factors at proposal rebuild.** `.emc_ensemble_draw_pools()` currently reconstructs `.emc_ensemble_q()` and Cholesky factors for every subject every active iteration, although proposal covariance is frozen between refreshes. Cache factor, log determinant, dimension and scale with clear invalidation for covariance/configuration changes. This removes recurring `O(n_subjects × p^3)` factorization work. Preserve RNG draw order.
2. **Combine covariance validation and factorization.** `.emc_empirical_covariance_ok()` performs positive-definiteness and conditioning checks repeatedly for current, old, mean and ensemble matrices. Return a validated standardized factor/condition status once per build and reuse it. Keep target inverse-Wishart validation stricter and separate from proposal repair.
3. **Measure optimizer refresh cost explicitly.** Warmup refreshes can optimize every subject/chain and score 256 diagnostic draws. Report optimization, finite-difference, diagnostic and transition costs separately. Skip redundant successful refreshes or use a measured warmup-only refresh policy where it improves total time to reliable samples. A smaller burn iteration count alone is not a speed result.
4. **Optimize in prior-whitened coordinates and reuse evaluations.** Curvature regularization already uses prior geometry; optimization still operates in raw coordinates. Whiten its optimization variables and cache repeated center/finite-difference evaluations where practical. Maintain full prior precision for joint models. Support-aware deterministic restarts are preferable to declaring a regularized saddle a usable mode.
5. **Tune ensemble size/rounds/cadence by cost.** Pool size pays for likelihoods; rounds primarily pay for hierarchy arithmetic. Compare these costs across cheap and expensive likelihoods, dimension and subject count. Use minimum subject/hierarchy ESS per second, not only average pool weight ESS or move rate. Cache the likelihoods without introducing stale-data/model keys.
6. **Report stage-local diagnostics.** Show ensemble movement, weight ESS, failures and actual jump sizes by stage/window, with lower quantiles and worst subjects. Separate independent/local CIS costs; their selection probabilities are not directly comparable. Include pooled coefficients in stopping and stall diagnostics.
7. **Use a corrected baseline and complete timing.** Benchmark against the repaired frozen particle sampler, with initialization and discarded stages included. Improvements against the biased pre-invariance sampler or only steady-state iteration time are insufficient evidence. Preserve fixture provenance, revision, dirty-worktree hash, options and saved RNG streams.

These are expected opportunities, not measured speedup claims from this review. Profile before accepting them, and retain exact density/RNG/backend checks when optimizing.

## 7. Validation additions and implementation order

### First: repair target and failure contracts

Implement F1–F3, F5 and F6 with analytic conditional/log-density tests. Freeze production by default under F4, or supply a complete alternative adaptive argument. Repair SEM backend conditioning and Laplace failure handling next. Regularize proposals in standardized coordinates and make target factorization strict.

### Second: verify integration and continuation

Add a test matrix crossing standard/diagonal/factor/infinite-factor/SEM with serial/queued/degraded workers, multiple subblocks, noncontiguous groups and correlated priors. Include documented unsupported nuisance/marginalized cases, so skipped ensemble support is intentional. Verify checkpoint continuation with distinct chains and real adaptive state. Normalize zero-dimensional covariate cases consistently; test both helper-produced empty matrices and supported direct settings rather than relying on accidental NULL behavior.

### Third: validate posteriors independently

Use known conjugate targets with correlated coordinates, an independent SEM implementation for nontrivial constraints/covariates, and simulation-based calibration for nonlinear models. Examine hierarchy variances, correlations and joint latent quantities as well as subject means. Include boundary, weak-identification, skewed and multimodal fixtures. Quantify Monte Carlo error before attributing discrepancies to bias.

### Fourth: establish practical benefit

Perform controlled ablations of ensemble, seeded mode starts, covariance shrinkage, efficient-mean tracking and local tuning. Use held-out fixtures and seeds. Report complete wall time to acceptable rank R-hat, bulk/tail ESS and MCSE, with the minimum relevant efficiency and uncertainty across runs. The existing short multi-seed burn comparisons are promising exploratory evidence; they are not a correctness proof or a comprehensive performance result.

## 8. Release recommendation

Retain the repaired CIS/ensemble architecture and the checkpoint state-preservation fix. Resolve the P1 target/failure defects and production-adaptation contract before presenting the full suite as mathematically validated. Treat mode seeding, shrinkage, tail scaling and local controllers as empirically assessed proposal strategies. Promote them on demonstrated total efficiency after independent posterior validation, with their assumptions and unsupported configurations recorded.

The strongest conclusion supported today is: **the central fixed-proposal construction is sound, the exercised tests pass, and several concrete hierarchy/reporting/failure defects plus an adaptive-correctness proof gap remain.**

## 9. Implementation status (2026-10-04)

| ID | Status | Notes |
|---|---|---|
| F1 | Fixed | Factor and SEM loading/K rows subtract fixed contributions; SEM B/G rows subtract fixed entries and refresh residuals from the complete row; per-column prior variances aligned (`byrow`), also in prior sampling and bridge. `test-hierarchy-conditionals.R` checks the requested conditional means against closed forms (two factors, correlated scores, heterogeneous variances, correlated disturbances); it fails on the previous code. |
| F2 | Fixed | `tauh` refreshed after each `delta[h]`; three-factor rate test. |
| F3 | Fixed (restricted) | `.sem_check_acyclic()` in `add_info_SEM()` and `make_sem_structure()` rejects any cycle among nonzero (free or fixed) B entries; recursive models only. |
| F4 | Re-argued, not frozen | The doubling schedule is adaptation at increasingly rare times (AirMCMC; Chimisov, Łatuszyński & Roberts 2018, arXiv:1801.09309), whose ergodicity/LLN/CLT results need growing gaps between adaptations plus a simultaneous drift condition on the kernel family — not vanishing kernel changes. The review is right that the old *comment* (diminishing adaptation) was the wrong theorem; the design itself has a standard justification, so non-vanishing refit noise is allowed. The drift condition is stated as an assumption. The schedule now counts transitions (`sample_transitions`), independent of thinning/trimming. |
| F5 | Fixed | `.sem_moments()` / `.sem_total_effects()` shared by Gibbs, prior sampling, `get_mu_implied()`, bridge and `group__IC_SEM()`; bridge and IC score subject-conditional densities. Bridge density checked against direct conditional densities; end-to-end SEM fit + IC + bridge smoke run. |
| F6 | Fixed | `.emc_pooled_rows_ll()` follows the failure policy: classified, row-wise retry, length check, NA/NaN/+Inf recorded, strict abort; counts travel back as an attribute and are recorded in the master under source "pooled" (stage summary + profile row). |
| F7 | Fixed | `filtered_samples_SEM()` provides `conditional_mu`; both paths condition on the subject prior mean. Probe shift is now (0, 0). |
| F8 | Fixed | Total diagnostic status, gate treats Inf/NA k as poor, failed refresh clears curvature and diagnostic, ensemble pool mean and covariance chosen as a pair. |
| F9 | Fixed; shrinkage removed | Regularisation in correlation coordinates with an explicit reference scale. The correlation shrinkage was then removed entirely: on lba71sim (p = 71) it applied intensities of 0.44-0.88, and the first production rebuild cut subject move rate 0.87 -> 0.39 (raw window covariance: 0.93; Laplace: 0.87). Its Frobenius-loss rationale does not fit an unscaled independence proposal. |
| F10 | Fixed | `.chol_factor()` flags repairs; a repaired group factor is never used as the single-block prior (numerical rejection instead). |
| F11 | Partial | Independent chains sub-block per chain. Diagnostic defaults (rank R-hat 1.01, more quantities) NOT changed: user decision, it lengthens every fit. |
| F12 | Partial | Scratch store respects `.EMC_STATIC_SAMPLE_FIELDS`. Broader continuation matrix and recorded sampler options still open. |
| F13 | Comment fixed | t component described as scale matrix (covariance 2cV), radial distribution stated; no density change. Controller calibration study open. |

Also fixed in passing: SEM prior sampling drew `K`/`G` columns with the wrong dimension and variance vector; `n_pars == 1` bridge `theta_mu` prior used `diag(scalar)`. Still open: SEM prior sampling uses inverse-Wishart df `max(a_d, d + 2)` while Gibbs/bridge use `a_d` (reporting only); stored `(alpha, eta)` pairs under `emc2.ensemble_latent_weight = "marginal"` (non-default).
