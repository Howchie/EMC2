# Handoff: hierarchical variance run instability

> **Resolved 2026-09-25.** The cause was a bias in the subject-level particle update, not mixing, RNG or the likelihood. See `plans/particle-sampler-fix-2026-09-25.md` (commits `1fd2e30f`, `3447e715`, `a994a8a7`). Refit affected models before comparing group variances or pWAIC.

Updated 2026-09-25. Investigation is open; no sampler fix has been made.

## Objective and current conclusion

Explain why two applied fits with the same data, prior, and tested likelihood give similar group means but substantially different between-subject variance posteriors. Determine whether the diagonal-gamma variance update, RNG handling, parallel execution, or convergence/stopping behavior is responsible. A compact saved-fit witness exists. A self-contained simulated fit that exhibits the failure has **not** yet been found.

The variance update appears algebraically correct and agrees with the saved conditional draws. The current sample-stage R-hat rule omits `sigma2`, but simply adding it at `max_gd = 1.1` would **not** have kept either applied fit running: both pass the actual stopping check at 1,000 draws. A continuation experiment then made the initially high-variance fit's SDs move toward the other fit's SDs across all 16 parameters, while its within-fit R-hat rose sharply. This supports a transient mixing/convergence failure, though its exact sampler mechanism is still unknown.

## Files and commands

- `WorkingTests/repro_hier_variance_saved_fits.R`: deterministic witness using the two saved project fits. Run from this repo with `Rscript WorkingTests/repro_hier_variance_saved_fits.R`. Set `REPRO_FITS_DIR` if the fit files are copied elsewhere. Requires installed EMC2 and the saved project files.
- `WorkingTests/extend_hier_variance_saved_fits.R`: rerunnable continuation experiment. Run with `Rscript WorkingTests/extend_hier_variance_saved_fits.R`; it takes several minutes and writes two summaries to `REPRO_OUT_DIR` (defaults to `/tmp`). It sets the parent RNG before each fit, so a rerun need not be bitwise identical to the exploratory continuation below.
- `WorkingTests/repro_hier_variance_run_instability.R`: original self-contained simulation search. **It uses the default `standard` inverse-Wishart prior**, whereas the applied fits use `type = "diagonal-gamma"`; therefore it is not a faithful reproduction of the applied prior case.
- `plans/hier-variance-run-instability-2026-09-25.md`: original investigation plan and findings added during this investigation.
- Applied files: `/data/work/PM/NirvanaHons_Nback/Fits/Match_Control_Exp1_rdm7s_mass_noA_t0R_vs2_mt_d42.RData` and `/data/work/PM/NirvanaHons_Nback/Fits/Match_Control_Exp1_rdm7s_mass_noA_t0R_UT0_vs2_mt_d42.RData`.

`Rscript WorkingTests/repro_hier_variance_saved_fits.R` passed and printed:

```text
Same prior values: TRUE
Same trial columns: TRUE
Largest tested log-likelihood gap: 0
Within-fit max sigma2 Rhat: 1.02, 1.011
Would stop at 1000 with sigma2 included: TRUE, TRUE
Six-chain max SD Rhat: 1.24854
v.d_match_LoadHigh SD medians: 0.737, 0.888
```

The likelihood check evaluates the **current** two C++ likelihood implementations at 24 saved subject parameter vectors, spanning four subjects and both fits. It does not prove the likelihoods were identical at fit time for every vector. The witness checks all trial data columns for all 42 subjects and all stored prior values; data object attributes differ. The models are RDMSWTN and RDMSWTN_UT with urgency fixed at zero. Both have 3 chains and 1,000 retained sample draws.

## Evidence collected

1. **Variance update:** `R/variant_diag_gamma.R:98-119` uses the expected inverse-gamma conditional. For each of 16 parameters in the saved fits, empirical mean `theta_var` differs by less than 0.8% from `(rate + RSS/2) / (shape + n_subjects/2 - 1)`. For `v.d_match_LoadHigh`, mean subject residual sums of squares are 25.049 versus 36.604, predicting variance means 0.582 versus 0.844; saved variance means are 0.580 versus 0.840. The sampled subject `alpha` dispersion is already different when the variance update runs. This argues against an arithmetic defect in that Gibbs draw; it does not rule out issues upstream of `alpha` or an unidentified sampler problem.
2. **Stopping:** `R/fitting.R:19-26` defaults the R-hat selection to `alpha,mu`. Calling actual `EMC2:::check_progress()` on each saved fit with `selection = c("alpha","mu","sigma2")`, `iter = 1000`, `max_gd = 1.1`, and old R-hat returns `done = TRUE` for both. Maximum all-selection R-hats are 1.048 and 1.039. Within-fit sigma2 R-hats are 1.020 and 1.011 (rank-normalized versions 1.015 and 1.009). Thus sample-stage inclusion of sigma2 alone is insufficient. Inclusion during burn could alter the path, which remains untested.
3. **Across fits:** pooling the six SD chains yields maximum R-hat 1.249; `v.d_match_LoadHigh` has roughly 1.185 and `v.M` roughly 1.219. The first retained SD draws already occupy separated regions; early and late retained medians are roughly stable. Thus separation predates the sample-stage endpoint. Many group means are close, but they are not exactly equal (for example `v.M` about -1.139 versus -1.345).
4. **RNG and parallelism:** within-fit chains have distinct streams and no evident aligned sigma-trace correlation. Repeated small same-seed RDMSWTN parallel fits with fixed settings were bitwise identical; LNR was too. `R/chain_pool.R:1765-1824` derives per-chain L'Ecuyer streams even when the caller uses Mersenne-Twister. Parallel and serial settings do not currently yield bitwise identical fits with the same seed, because `R/fitting.R:658-679` draws parent-side `sample(...)` proposals while serial and forked chain work advance the parent RNG differently. This is a reproducibility limitation across core configurations, not evidence by itself of biased posteriors.
5. **Synthetic search so far:** a 4-parameter hierarchical RDMSWTN with `diagonal-gamma`, 12 subjects, 30 trials per design cell, 500 sample draws, and seeds 1/2/3 gave close SD medians and pooled R-hat 1.002. A separate small nonlinear mixture stress test was also stable. These are negative controls, not reproductions. The original seven-run script uses the wrong prior for the applied case and has not demonstrated this exact failure.

Relevant additional code: proposal adaptation at `R/fitting.R:314-335`, parent-side proposal construction at `R/fitting.R:658-679`, and the sample/continuation stop handling at `R/fitting.R:117-120`. All applied subjects have one parameter proposal block (`attr(data, "components")` is all 1), so a partial-parameter varying-mask optimization does not appear relevant.

## Continuation experiment completed

From 01:40 to 01:48 UTC, the two saved fits were extended sequentially with `selection = c("alpha","mu","sigma2")`, `max_gd = 1.1`, `cores_for_chains = 3`, `cores_per_chain = 9`, `step_size = 100`, and `max_tries = 6`. Fit 1 reached 1,500 draws; fit 2 reached 1,600 because its R-hat increased and it exhausted `max_tries`. The original 1,000 draws and subsequent draws were compared separately. Exploratory summaries are `/tmp/hier_variance_extend_1.rds` and `_2.rds`; `/tmp` is ephemeral. The rerunnable script above recomputes them.

| Measure | Fit 1 | Fit 2 |
| --- | ---: | ---: |
| Sample draws after continuation | 1,500 | 1,600 |
| `v.d_match_LoadHigh` SD median, original draws | 0.737 | 0.889 |
| `v.d_match_LoadHigh` SD median, new draws only | 0.763 | 0.769 |
| Final within-fit max `sigma2` R-hat | 1.012 | 1.200 |

The maximum absolute gap between fits across all 16 SD medians fell from 0.153 in the original draws to 0.0081 in the new draws; the median absolute gap fell from 0.0174 to 0.0040. In fit 2, the all-selection R-hat rose from 1.039 at 1,000 draws to 1.237 at 1,600. This is direct evidence that its apparently converged 1,000-draw interval was misleading and that later draws moved toward the other fit. It does **not** prove the new segments fully converged or establish the precise upstream cause. Keep in mind the first continuation has 500 new draws and the second has 600.

The continuation call for each fit is:

```r
emc <- run_emc(emc, stage = "sample",
  stop_criteria = list(iter = 500, max_gd = 1.1,
    selection = c("alpha", "mu", "sigma2"), omit_mpsrf = TRUE),
  cores_for_chains = 3, cores_per_chain = 9,
  step_size = 100, max_tries = 6, verbose = TRUE)
```

In `run_emc`, a continuation request adds `stop_criteria$iter` to prior retained iterations. Thus `iter = 500` means 500 new draws, reaching 1,500 total. An earlier attempt used `iter = 1500` and was interrupted immediately after recognizing it would request 1,500 **additional** draws; it produced no results.

## Next decisions and experiments

1. Examine fit 2's 1,000-to-1,600-draw `alpha`, `mu`, and `theta_var` trajectories by chain to identify which subjects/parameters moved and whether the transition was abrupt. The exploratory continuation saved **summaries**, not the full extended `emc`, so rerun the script or adapt it to save full objects if these trajectories are needed.
2. Compare each fit's subject-level `alpha`, `mu`, and `theta_var` from the earliest available stage, especially after preburn, burn, and adaptation. Identify the first stage where the fits separate and whether each fit's three chains share an initialization basin. Examine joint likelihood and prior density along an interpolation between representative states, and compare chain-specific log posterior values.
3. Build a faithful small **diagonal-gamma** reproduction with the applied 16-parameter design and fewer subjects/trials, or subset the applied data while preserving the weakly identified dimensions. Test repeated seeds and both mathematically equivalent model forms. State explicitly if no small reproduction emerges.
4. Test sampler mechanics one at a time: proposal adaptation, full-block vs smaller-block updates, preburn/burn length, and initializations. Inspect stage-specific proposal covariance and acceptance. Only change sampler code once a diagnostic isolates the mechanism. A stricter stopping rule or longer fixed minimum sample phase may help, but the original sigma2 R-hats passed 1.1, so neither has yet been validated as a sufficient fix.
5. For RNG, test same-seed repeats under identical core settings separately from cross-core reproducibility. The former passed in small cases. Cross-core bitwise invariance remains incomplete and can be fixed separately if desired, but should not be presented as the cause of the variance discrepancy without a causal test.

## Workspace caution

The repository already had many modified/untracked files when this investigation began, including RDM and C++ model files. This investigation created `WorkingTests/repro_hier_variance_saved_fits.R` and added the sigma2 test section to the existing untracked plan. It did not edit the sampler implementation. Preserve the existing dirty worktree and review `git status --short` before editing.
