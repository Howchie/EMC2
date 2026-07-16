# Plan: addressing review feedback on correlated BAwL (commit 3852f10b)

Review target: branch `bawl-correlated-race`, commit `3852f10b` ("add correlated BAwL").
Reviewer verdict: the architecture (conditional ordinary race likelihood + 1-D
Gauss–Hermite integration over a shared factor) is right; five issues need work
before fitting the n-back data. The reviewer's addendum retracts the
two-free-loadings identifiability criticism *given* a one-rho-per-trial-type
design, but requires that the single coefficient be turned into two
role-specific signed row loadings.

## Current state (important context)

The working tree already contains **uncommitted** changes to
`src/particle_ll.cpp` (+157/−14) that:

1. **Fix review issue #1** — `log_den` is now initialised to `R_NegInf` for
   truncated trials and left at `0.0` (= log 1) for untruncated ones
   (`src/particle_ll.cpp:6338-6351`), and accumulation only happens for
   trials in `trunc_mask`.
2. Batch the truncation normaliser: conditional-on-node race survivors at
   uniform LT/UT are evaluated via the column-major `logS_at_t` callback with a
   per-trial scalar/GSL fallback, instead of one
   `get_trunc_normaliser_rowmajor_cpp` call per trial per node.

So workstream 1 below is "finish and lock in", not "implement from scratch".

Key code locations:

- Likelihood: `c_log_likelihood_bawl_correlated`, `src/particle_ll.cpp:6252`
  - per-row loading conversion `sign(rho)*sv*sqrt(|rho|)*z`,
    `sv*sqrt(1-|rho|)`: `src/particle_ll.cpp:6370-6387`
  - exact rho==0 routing to `c_log_likelihood_race`: `src/particle_ll.cpp:6289-6304`
- GH rule (fixed 10 nodes): `src/gh_quad.h`
- R simulator `rBAwL_corr`: `R/model_LBA.R:409-433`
- C++ simulator `rbawl_corr_cpp`: `src/model_rng.cpp:329-364`
- Model definition (rho p_type, pnorm→(−1,1) transform, Ttransform):
  `R/model_LBA.R:503-554`, `BAwLcorr` wrapper at `R/model_LBA.R:592`
- Tests: `tests/testthat/test-bawl-correlated.R`

---

## Workstream 1 — truncation denominator bug (DONE)

**Status:** fixed, uncommitted.

Completed:

1. `has_trunc_trial` is the single structural mask used for denominator
   initialisation and accumulation.
2. The independent truncation-identity regression test is in
   `tests/testthat/test-bawl-correlated.R`; it fails on bare `3852f10b` and
   passes with the fix (see the verification record below).
3. The changes remain in the working tree for review; no commit was created
   so the existing uncommitted work is not folded into an unsolicited Git
   history change.

## Workstream 2 — one correlation per trial type, role-specific signed loadings

**Problem (per reviewer's addendum):** the kernel converts every row
independently via `a_j = sign(rho_j)·sqrt(|rho_j|)`. Two consequences:

- If the *same signed* rho is mapped onto both active rows (e.g. `rho ~ TrialType`),
  the induced correlation is `a·a = |rho|`: the sign is silently
  unestimable even though the parameter is allowed to be negative.
- If each row gets its own free rho (as the current test's `rho ~ 0 + lR`
  does), only the product `sign(ρ_T ρ_N)·sqrt(|ρ_T ρ_N|)` is identified — a
  ridge, plus a global sign flip symmetry.

**Target parameterisation:** one direct correlation `r_c ∈ (−1,1)` per trial
type, converted centrally into role-specific loadings:

- correct/matching racer: `a = +sqrt(|r_c|)`
- incorrect/mismatching racer: `a = sign(r_c)·sqrt(|r_c|)`

so that `Corr(D_correct, D_incorrect) = r_c` exactly, both racers keep
conditional variance `sv²(1−|r_c|)`, and the PM racer keeps loading 0.

**Implementation choice:** do the role mapping in BAwL's `Ttransform`
(`R/model_LBA.R:529`), which receives `(pars, dadm)` and whose output feeds
*both* the C++ likelihood and `rfun` (the same mechanism already derives
`b = B + A` and `lambda_g`/`lambda_k`). This is inside the model definition,
so a user formula cannot bypass it, and it keeps the C++ kernel's existing
per-row `sign(rho)·sqrt(|rho|)` conversion unchanged:

```r
# inside Ttransform, when correlated:
#   correct row  -> |rho|   => kernel loading +sqrt(|rho|)
#   incorrect row -> rho    => kernel loading sign(rho)*sqrt(|rho|)
if (!is.null(dadm$lM)) pars[, "rho"] <- ifelse(dadm$lM, abs(pars[, "rho"]), pars[, "rho"])
```

Details:

1. Role indicator = `dadm$lM` (latent match). For the n-back application this
   gives exactly the reviewer's table: N trials — target accumulator is
   correct; L/F trials — non-target accumulator is correct. PM racer is held
   at `rho = 0` by a constant, unaffected by the mapping.
2. **Require `lM`** when `correlated = TRUE` and error with a clear message if
   the design lacks a `matchfun` — do not silently fall back to raw per-row
   loadings. (If a raw-loading escape hatch is ever wanted, make it an explicit
   argument, not a fallback.)
3. Verify both simulators inherit the mapping for free: `rfun` receives
   Ttransformed `pars`, and `rBAwL_corr` / `rbawl_corr_cpp` already apply
   `sign(rho)·sqrt(|rho|)` per row (`R/model_LBA.R:420-427`,
   `src/model_rng.cpp:352-357`). No C++ changes needed if the mapping lives in
   Ttransform. Confirm with the sign-sensitivity test (Workstream 4).
4. **`rho` must be shared within a trial.** Under the new semantics `rho` is a
   trial-level quantity (the correlation of the racing pair), so its design
   formula must not vary over `lR`: the same value sits on every participating
   row of a trial, and Ttransform's `lM` mapping is what differentiates the
   rows. Formulas like `rho ~ 0 + lR` (the current test) belong to the *old*
   per-accumulator-loading semantics and reintroduce exactly the ridge the
   reviewer flagged. Correct idioms: `rho ~ 1`, or `rho ~ TrialType` /
   `rho ~ S:Cond` for cell-specific correlations.
5. **Excluding the PM racer** can no longer go through `rho ~ 0 + lR` +
   constant. Use a derived row-level *participation mask* via
   `design(functions = ...)`, e.g.
   `coupled = function(d) factor(d$lR != "pm", c(FALSE, TRUE), c("no", "yes"))`,
   then `rho ~ 0 + coupled` (crossed with trial-type factors as needed) with
   `constants = c(rho_coupledno = 0)`. Both pair rows then share one free
   signed value; the PM row is pinned to 0. Equivalently, with the PM level as
   the factor reference: `rho ~ coupled` with `constants = c(rho = 0)` pins
   the intercept, so PM rows get 0 and the coupled rows carry the single free
   slope. (Constants are on the sampled scale; sampled 0 maps to correlation
   0 through the pnorm-based transform, so both idioms pin exact independence.)

   Note `coupled` and `lM` are different columns doing orthogonal jobs and
   must not be conflated. `coupled` says *who loads on the shared factor* and
   is trial-invariant (target and non-target always, PM never). `lM` says
   *who is the correct/reference racer on this trial* and is trial-varying;
   it drives only the sign convention in Ttransform. For n-back:

   | trial | row        | coupled | lM  | mapped rho |
   |-------|------------|---------|-----|------------|
   | N     | target     | yes     | T   | +\|ρ\|     |
   | N     | non-target | yes     | F   | ρ          |
   | N     | pm         | no      | F   | 0          |
   | L/F   | target     | yes     | F   | ρ          |
   | L/F   | non-target | yes     | T   | +\|ρ\|     |
   | L/F   | pm         | no      | F   | 0          |

   Using `lM` itself as the participation factor (`rho ~ 0 + lM` with
   `rho_lMFALSE = 0`) would be wrong: it would zero the *mismatching pair
   member* on every trial and kill the correlation entirely.
6. **Enforce, don't just document, the within-trial constraint.** The formula
   machinery will happily map any factor onto `rho`, including `lR`, so
   shared-within-trial is purely conventional unless checked. Because `rho`
   depends on parameters only through its design matrix, the check can be
   structural and one-time at dadm construction (make_emc): within each
   trial, the non-structurally-zero rows of the `rho` design matrix must be
   identical (≤ 1 distinct nonzero design row per trial). That guarantees
   within-trial equality for *every* parameter vector, catches `rho ~ 0 + lR`
   at model-build time with a clear error, and costs nothing in the sampling
   hot path. (A per-call Ttransform value check is the fallback if wiring a
   model-specific dadm check proves awkward.)
7. Update docs (`man/BAwLcorr.Rd`, roxygen in `R/model_LBA.R`): `rho` is the
   *pairwise correlation between the correct and incorrect racers' drifts per
   trial type* (of the underlying Gaussian model — see Workstream 3), not a
   per-accumulator loading; document the shared-within-trial requirement and
   the `coupled` exclusion idiom from points 4–5.
8. Update `tests/testthat/test-bawl-correlated.R`: replace `rho ~ 0 + lR`
   (which under the new semantics puts different correlations on rows of the
   same draw) with the `rho ~ 0 + coupled` + constant idiom from point 5, and
   add a test that `rho ~ 0 + lR` is rejected by the point-6 check.
9. Keep the existing `pnorm`-based (−1,1) transform, `exception = 0`, and the
   endpoint guard `fmax(residual, 1e-12)` — all still correct for a direct
   correlation parameter.

## Workstream 3 — posdrift semantics: adopt the jointly truncated correlated Gaussian

**Problem:** the current code implements the conditional-truncation hierarchy
(`Z ~ N(0,1)`, then `D_j | Z ~ TN₍₀,∞₎(μ_j(z), σ_j(z))`). Under
`posdrift = TRUE` this means `rho` changes each accumulator's *marginal* drift
distribution, so a fit improvement from `rho` can partly buy univariate shape
flexibility rather than dependence. The reviewer recommends the jointly
truncated correlated Gaussian instead, which nests standard positive-drift
BAwL cleanly.

**Decision:** adopt the jointly truncated model as the meaning of
`BAwLcorr(posdrift = TRUE)`. (No back-compat concern — nothing has been fit
with this model yet.) `posdrift = FALSE` is already exact and unchanged.

**Likelihood implementation** (in `c_log_likelihood_bawl_correlated`): with
`q_k(z) = Φ(μ_k(z)/σ_k(z))` per active accumulator, reweight nodes by
`Π_k q_k(z)` in numerator *and* denominator:

```
L(y) = ∫ φ(z) [Π q_k(z)] L⁺(y|z) dz  /  ∫ φ(z) [Π q_k(z)] Z⁺(z or 1) dz
```

- Per node `q` and unique trial `j`, compute
  `log_pos_j = Σ_{k active} pnorm(μ_k(z)/σ_k(z), log.p) `
  using the already node-shifted `pars_q` values — cheap (`n_active` Φ calls),
  and reuse it for both sums.
- Numerator accumulation becomes `log_w + log_pos_j + node_ll[i]` (note
  `log_pos` is per unique trial; apply through the `expand` map, same as
  `log_den` is applied today).
- Denominator: under this model **every** trial gets a denominator, not just
  truncated ones — `log_den[j] = logsumexp_q(log_w + log_pos_j + log_z_j)`
  with `log_z_j = 0` for untruncated trials and the conditional truncation
  survivor mass otherwise. This *simplifies* the Workstream-1 logic: init all
  `log_den` to `−Inf` when posdrift, drop the truncated/untruncated split, and
  the positivity normaliser cancels correctly against truncation.
- Skip the reweighting entirely when `posdrift = FALSE` (kernel selection
  already distinguishes BAwL vs BAwLIO via `c_name`), and note the rho==0
  fast path stays exact: at zero loading `Π q_k` is node-constant and cancels.
- Respect RACE masks / `isok` when forming the active set for `Π q_k`
  (excluded rows must not contribute a positivity factor).

**Simulator implementation** (both `rBAwL_corr` and `rbawl_corr_cpp`): replace
per-accumulator `rtnorm` draws with trialwise rejection, which samples the
jointly truncated Gaussian exactly:

1. draw `z ~ N(0,1)`;
2. draw all active drifts *untruncated* `N(μ_k(z), σ_k(z))`;
3. if any active drift ≤ 0, redraw the whole trial (including `z`).

Guard with a max-iteration cap and an informative error (acceptance can be
poor when a mean drift is far below zero — that is a parameter-region warning,
same spirit as standard posdrift LBA). `posdrift = FALSE` keeps plain normal
draws.

**Docs:** state explicitly in the roxygen that with `posdrift = TRUE` the
model is a correlated Gaussian drift vector conditioned on all active drifts
being positive; `rho` is the correlation of the *underlying untruncated*
Gaussian; marginals coincide with standard positive-drift BAwL.

## Workstream 4 — substantive tests (current tests establish plumbing only)

Add to `tests/testthat/test-bawl-correlated.R` (or a new
`test-bawl-correlated-numeric.R`):

a. **Independent nesting.** For all-zero rho:
   `ll_BAwLcorr == ll_BAwL` to machine precision (the code routes this through
   the ordinary path, so `expect_identical`-level agreement), covering: plain,
   truncated (LT/UT), omissions/censoring, clock processes (`lambda_g`,
   `lambda_k`), posdrift TRUE and FALSE.

b. **Quadrature reference.** Slow R reference: for a small hand-built `pars`
   matrix, integrate the conditional race likelihood over z with 40–80 GH
   nodes (or `integrate()`), including the positivity reweighting from
   Workstream 3, and compare trialwise log-likelihoods against the C++ path
   for rho ∈ {−.7, −.3, .3, .7} × posdrift {TRUE, FALSE}. Tolerance set by the
   Workstream-5 convergence results (the 10- vs 80-node gap), not machine eps.

c. **Truncation identity.** Verify numerically that
   `ll_truncated == ll_unconditional − log P(LT < T < UT)` with the
   denominator computed *independently* (R-side high-precision quadrature over
   z of the race survivor difference, positivity-reweighted under posdrift).
   This is the test that catches the `log_den = 0` bug; confirm it fails on
   bare `3852f10b`.

d. **Sign sensitivity / role mapping.** With the Workstream-2 mapping, check
   that `r_c = +.6` and `r_c = −.6` produce *different* trialwise likelihoods
   on the same data, and that flipping which accumulator is `lM`-matching
   flips the effective sign. (Guards against the "same signed rho on both
   rows ⇒ |rho|" failure mode returning via refactor.)

e. **Simulator dependence check.** For `posdrift = FALSE` (where drift
   correlation is exact), simulate large-n with the C++ simulator and verify
   the empirical drift correlation ≈ `r_c` — either by exposing drifts from
   the simulator for testing, or indirectly: check the trialwise ll of the
   simulated data is maximised near the generating `r_c` on a coarse rho grid.

f. **Dependence recovery** (slow; `skip_on_cran` + long-runner or an
   `inst/`-level script): simulate one direct correlation per trial type,
   fit, and check recovery for r ∈ {−.7, −.3, 0, .3, .7}. Do **not** attempt
   recovery with two free row-specific loadings — that is a ridge by design.

## Workstream 5 — GH node-count convergence

1. Replace the hardcoded 10-node table in `src/gh_quad.h` with runtime
   Gauss–Hermite rules for arbitrary n.  The bundled GSL copy lacks the fixed
   Hermite API, so the implementation uses an equivalent cached Golub–Welsch
   rule.  The production default is selected from the convergence study below
   rather than retaining the biased 10-node rule.
2. Make the node count controllable for testing (e.g. an R option /
   env var read once per call, or an argument threaded through the model
   list), so the convergence study runs without recompiles.
3. Convergence study (script, results recorded in this file or a notebook):
   `max_i |ℓ_i^(Q) − ℓ_i^(80)|` for Q ∈ {10, 20, 40} over the reviewer's
   grid — |r| ∈ {.2, .5, .8, .95}; low/near-zero mean error drift; large sv;
   posdrift on; truncated and omitted observations. Target ≤ 1e−5 per trial
   comfortable, ≤ 1e−4 acceptable *if smooth in parameters* (check a fine rho
   sweep for kinks, since MCMC cares about surface smoothness more than bias).
4. Set the production default from the results (plausibly 10 for |r| ≤ .8 but
   20+ near |r| → 1 where the residual SD collapses and the positive-truncated
   kernel gets sharp). If needed, pick the rule adaptively from
   `max |rho|` in the particle — cheap since rho is a parameter column.

## Sequencing

1. Commit Workstream 1 (bugfix + truncation-identity test 4c). Small,
   self-contained, unblocks everything else.
2. Workstream 2 (role-mapped single correlation) + tests 4a/4d — changes the
   parameter meaning, so do it before any reference tests bake in semantics.
3. Workstream 3 (jointly truncated posdrift model) — likelihood + both
   simulators + docs in one commit so semantics never straddle commits.
4. Workstream 4b/4e reference tests against the final semantics.
5. Workstream 5 convergence study; set default node count; then 4f recovery
   as the final end-to-end check before fitting the n-back data.

## Verification notes

- Verify via a fresh `R CMD INSTALL` into a temp `R_LIBS` library, not
  `library(EMC2)` against the stale system install (known trap in this repo).
- No GitHub Actions ran for `3852f10b`; run the full local testthat suite
  (including the LogicalRules censoring/truncation suites touched by the
  commit) before and after each workstream.

## Verification record (2026-07-16)

- Workstream 1 regression: the independent truncation oracle fails on bare
  `3852f10b`.  For the same one-trial harness it reports
  `ll_truncated = -0.678424172132`, `ll_unconditional = -0.941129025204`,
  while the independently integrated old-semantics window probability is
  `log P(window) = -0.463821058598`; the identity residual is `-0.201116205526`.
  The fixed implementation passes the corresponding independent test in
  `test-bawl-correlated.R`.
- The convergence grid covers regular, truncated, omitted, and near-zero
  mean/large-`sv` cases at `|rho| = {.2, .5, .8, .95}`.  Against Q80, the
  maximum absolute per-trial gaps were Q10 `8.773881`, Q20 `0.769456`, and
  Q40 `0.009301608`; through `|rho| <= .8`, Q40 was below `1e-4`.
  A fine Q200 sweep over rho = `.10, .15, ..., .95` was finite and had maximum
  absolute second difference `0.04169035`.
- The production rule is therefore adaptive rather than fixed at the old
  10-node rule: 40 nodes through `|rho| <= .8`, 80 through `.9`, and 200
  above `.9`; `EMC2_BAWLCORR_GH_N` still overrides this for testing.  At the
  hard `.95`, near-zero/large-`sv` point, Q200 differed from Q256 by less than
  `4e-7`.
- The bundled GSL copy exposes adaptive quadrature but not the fixed Hermite
  API named in the original proposal.  `src/gh_quad.h` consequently uses the
  equivalent Golub–Welsch Hermite construction with a per-node-count cache.
- Fresh installation and the focused correlated suite pass; the opt-in long
  recovery suite also passes.  The full local suite was run serially; its
  remaining failures are unrelated pre-existing tests (`map`, mixed GNG
  dispatch, stop-success/stopS, and wald-logspace), while the custom-likelihood
  regression passes after the guarded model-spec check in `design()`.
