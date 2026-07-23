# Increment 2 handoff — wire t0 marginalization into the sampler

Status at this checkpoint (2026-07-23, branch `bawl-correlated-race`):

- **Increment 1 DONE + validated** — `calc_ll_oo` gained a `marginalise` arg
  (`Nullable<List>` of `param, mu, sigma, n_nodes[, span, eps]`). `NULL` ⇒ byte-identical
  existing path (zero cost, `FAIL=0` on existing tests). Set ⇒ `calc_ll_oo_marginal_core`
  integrates one shared t0 out of the COMPLETE subject likelihood by Gauss-Legendre
  quadrature on the log-t0 axis, reusing the registered kernel per node (model-agnostic).
  Matches the R reference `WorkingTests/marginal_t0_lib.R::marginal_ll_t0` to 3.4e-13.
- **Increment 1b DONE** — PtMapper fusion measured not worth it (setup = 3.2% of a call).
  The adaptive route now uses a 7-node pilot scan followed by a recentered composite
  Gauss-Legendre rule with explicit tail panels. The production default is 40 final
  nodes (47 kernel calls including the pilot); 40 is below the old fixed GL(80) cost
  while meeting the convergence check for LBA, RDM, and RDMSWTN with positive `sv`.
  Exported
  `calc_ll_oo_marginal_nodes(...) -> list(nodes[K] sampled-scale, log_terms[np×K])` for
  reconstruct-at-storage; `softmax_k(log_terms)` reproduces the ll exactly and
  `E[t0|y]=Σ w_k exp(nodes_k)` recovers near truth.

The **C++ likelihood + reconstruction layer and the Increment 2 sampler wiring are
implemented**. The Stage 0 scratch check now covers both GNG/no-go and standard-race
configurations, including RDMSWTN with `sv > 0`.

## Decisions locked (do not re-litigate)

- **API:** `design(marginalise = "t0")`. Generally available on ANY race model — NO
  nogo/time gating (a plain LBA consumes it fine; unobserved-finish cases are handled by the
  same whole-model integral with zero special-casing).
- **η (fixed, Stage 1):** the EXISTING t0 prior mean/sd. No new η API.
- **"Use the posterior as sampled":** reconstruct t0 at sample-STORAGE time from the node
  weights and write it into the stored `alpha[t0]`. So `predict()`/`make_data()` stay
  byte-identical; t0 is a THIRD parameter category (not sampled, not constant).
- **Scope:** race models only (DDM excluded — not a race model).

## Touch points (files/lines are at this checkpoint; re-grep before editing)

1. **`design()` — carry the flag.** `R/design.R:81` (`design <- function(...)`). Add a
   `marginalise = NULL` arg; validate each named param is a real p_type with `~1` (single
   shared coord) and a finite lower bound in `model()$bound$minmax`. Store as
   `attr(design, "marginalise")`. Thread onto `dadm` like other design attrs (see
   `R/design.R:544` `attr(dadm,"sampled_p_names")`).

2. **`sampled_pars()` — exclude from the SAMPLED vector, only in the sampling path.**
   `sampled_pars.emc.design` at `R/design.R:1433` (and the `.emc`/`.emc.prior` methods:
   `R/s3_funcs.R:1276`, `R/priors.R:650`). The marginalised names must be dropped from the
   vector the sampler treats as `alpha`, BUT `make_data`/`predict` still need t0 as an
   ordinary parameter. Cleanest: keep `sampled_pars` returning the full set, and drop the
   marginalised coords where the SAMPLER builds `alpha` (step 3), OR add an
   `include_marginalised=TRUE/FALSE` switch and pass FALSE from the sampler only. Prefer the
   former (fewer call-site changes; verify all `sampled_pars` consumers).

3. **`new_particle` — propose without t0, evaluate via the marginal.** `R/sampling.R:334-433`.
   - Build the proposal MVN over `alpha` WITHOUT the marginalised coords (drop rows/cols of
     mu/Sigma for those names).
   - The shared-parameter likelihood call at `R/sampling.R:398` (`calc_ll_manager(...)`)
     must route through the marginal: pass `marginalise = list(param="t0", mu=<prior mean>,
     sigma=<prior sd>, n_nodes=...)` down to `calc_ll_oo`. Thread a `marginalise` arg through
     `calc_ll_manager` (`R/sampling.R:732`) → `calc_ll_oo(..., marginalise=)`. The proposal
     matrix still needs a t0 COLUMN (placeholder; its value is overwritten per node) — inject
     it (e.g. as a constant-like column) so the design mapping works.
   - η source: read the t0 prior's mean/sd from the prior object attached to the sampler
     (`pmwgs$prior` / `get_prior`-style). Convert to the sampled (log) scale used by `alpha`.

4. **Reconstruct-at-storage.** Where an accepted iteration is stored (`sample_store` /
   `fill_samples`, `R/sampling.R:619-657`): for the accepted θ, call
   `calc_ll_oo_marginal_nodes(...)`, form `w = softmax(log_terms[accepted,])`, draw one node
   per subject with R's RNG (`sample.int(K, 1, prob=w)`), and write `nodes[k]` into the
   stored `alpha[t0]` slot (sampled/log scale). Keep the t0 slot present in the samples
   array so `predict`/`fill_samples` shapes are unchanged. New subjects (no data): draw
   `t0 ~ p(t0|η)` from the fixed prior instead.

5. **Gibbs / group level.** `gibbs_step_standard` (`R/variant_standard.R:215-364`) derives
   hyperparameters from `alpha` residuals via conjugate IW. In Stage 1 t0's group level is
   NOT updated (η fixed) — ensure the marginalised coords are excluded from the IW residual
   path so they don't contaminate Sigma, OR (simpler) let the reconstructed t0 sit in `alpha`
   and be updated normally but with η overridden by the fixed prior. Decide by testing that
   the other parameters' Sigma is unaffected. (Stage 2 replaces this with a proper
   hierarchical-η variant + Metropolis on η.)

## Verification for Increment 2

- **Gating unchanged:** a model WITHOUT `marginalise` must be byte-identical + same timing.
  Run `tests/testthat/test-logicalrules-*.R` and a standard-race test; diff lls/timing.
- **Validation set (prove "any race model consumes it"):** plain LBA (no nogo/time) +
  standard-race GNG + LogicalRules GNG + an `lR=time` race. End-to-end `fit()` on simulated
  data with known t0; check (a) the t0↔race-speed pair no longer has a stuck direction
  (ESS/R-hat on the confounded pair vs the naive sampler), (b) θ and reconstructed t0 recover
  the generating values, (c) `predict()` runs unchanged and PPC looks right.
- **Reconstruction correctness:** posterior E[t0|y] per subject ≈ generating t0; the stored
  categorical draws give correct spread (not the RB mean — that understates within-draw
  uncertainty).
- Follow the stale-library discipline: install into a temp `R_LIBS`, check the harness path
  banner (see memory `project_emc2_verification_trap`).

## Key references

- C++: `src/particle_ll.cpp` — `calc_ll_oo_marginal_core` / `calc_ll_oo_marginal` /
  `calc_ll_oo_marginal_nodes` (search these names); `calc_ll_oo` marginalise branch near its
  top. `gl_quad.h::gl_get_rule`.
- R reference / Stage 0: `WorkingTests/marginal_t0_lib.R`, `WorkingTests/stage0_marginal_t0.R`.
- Plan: `~/.claude/plans/please-consider-the-below-magical-tulip.md` (full context + all
  decisions).
- Memory: `project_emc2_marginal_t0.md`.
