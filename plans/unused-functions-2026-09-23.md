# Unused R functions (2026-09-23)

**Status (2026-09-23):** "Superseded" group deleted; "Duplicated wrappers" `waic_*`/`loo_*` exported (man/waic_loo.Rd); `ess_*` and the "Legacy" group left in place.

Source: graphify call graph of `R/` + `src/`, then verified by hand. Each function
below is not exported, not registered as an S3 method, and its name appears
nowhere in `R/`, `src/`, `tests/`, `inst/`, `vignettes/` or `man/` except its
own definition. Anything reached via `switch()`, S3 dispatch or a pasted name
was checked and excluded. Out-of-package callers (`EMC2:::` from analysis
scripts) were NOT checked.

## Superseded (safe to delete)

| Function | Where | Replaced by |
|---|---|---|
| `rexGaussian` | R/model_SS.R:368 | `rREXG` closure inside `REXG()` (R/model_LNR.R:458): zero-truncated `rtexG` + `t0` + `ok` + timed guess. `rexGaussian` is untruncated, has no `t0` and ignores `ok`, so it doesn't match REXG's likelihood. |
| `rexG` | R/model_SS.R:343 | Only caller is `rexGaussian`. `rtexG` is the live sampler. |
| `rBTAwLSeparate`, `rBTAwLSustained` | R/model_rng.R:127,139 | BTAwL constructors call `.rfun_BTAwL(..., mode=)` directly (R/model_BTAwL.R:627-787). The siblings `rBTAwL`/`rBTAwLTransient` are still referenced. |
| `group_dist_infnt_factor`, `prior_dist_infnt_factor`, `get_all_pars_infnt_factor` | R/variant_infnt_factor.R:309-318 | Old IS2 hooks. No other variant still has its versions of these, and `define_variants.R` never dispatches to them. |
| `apply_forward_fill` | R/trend.R:611 | Forward fill now happens in the C++ kernels (`BaseKernel::forward_fill_missing_outputs`). |
| `get_unique_rows` | R/group_design.R:219 | Unused helper. |
| `add_recalculated_pars` | R/map.R:130 | Unused helper. |
| `add_Ffunctions` | R/make_data.R:1110 | Ffunctions are now applied in design.R:1151-1156. |
| `check_CR` | R/sampling.R:1851 | Unused. |

## Duplicated wrappers (delete, or keep for interactive use)

| Function | Where | Note |
|---|---|---|
| `waic_subject`, `waic_pooled`, `loo_subject`, `loo_pooled` | R/waic.R:203-226 | Thin wrappers over `waic_from_ll`/`loo_from_ll`. `compare()` (R/statistics.R) does the same thing inline and doesn't call them. Either export them as a public API or delete them. |
| `ess_mean`, `ess_median`, `ess_sd`, `ess_tail` | R/diagnostics.R:234-249 | `@noRd`, never called; `ess_summary` doesn't use them. |

## Legacy / unfinished (your call)

| Function | Where | Note |
|---|---|---|
| `plot_mcmc` | R/plotting.R:641 | roxygen deliberately commented out (`# #'`). |
| `plot_fit_choice` | R/plotting.R:95 | roxygen commented out, ~160 lines. |
| `order_pp`, `robust_density` | R/plotting.R:863, :17 | Helpers left over from old plotting code. |
| `compare_obs_vs_postpred` | R/plotting_ss.R:1129 | Stop-signal PPC helper, never wired in. |
| `compare_MLL` | R/statistics.R:820 | roxygen commented out. |
| `std_error_IS2` | R/statistics.R:377 | IS2 bootstrap SE. |
| `SBC_hierarchical`, `make_smooth` | R/SBC.R:302, :1041 | SBC code paths not reachable from `run_sbc`. |
| `check_prior` | R/priors.R:242 | Unused validation helper. |
| `dhalft` | R/statistics.R:501 | Half-t density. |
