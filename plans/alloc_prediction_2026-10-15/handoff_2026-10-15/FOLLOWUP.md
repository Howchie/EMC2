# Follow-up, 2026-10-15

Read [the findings](../eval_operator/FINDINGS.md) and [computed results](../eval_operator/RESULTS.md) before repeating the handoff's untested directions.

The allocation performance issue is still open. A correction to `.emc_alloc_eff_inv()` was implemented and tested: targets below kappa's minimum tau return an unreachable raw target. The scalar and batched pricing regression tests cover this, and the noisy stage-tau test pins its generator and uses a longer simulated reference. All four requested move test files passed at the full test level on that private build.

Around 07:00–07:02 UTC another session removed kappa and its tests and made further allocation changes. Those edits were retained. The tested version is preserved in `eval_operator/evaluated_allocation.R`; the experiments now source that snapshot for reproducibility. `kappa_floor_regression.R` tests its correction independently. The historical full test pass does not validate the subsequently edited allocator, and the kappa fix was not reintroduced into the current runtime files.

New offline evidence:

- Removing subject credit alone changes predicted cost by roughly -5% to +3% across the five saved PM states, with two states essentially unchanged. It is not a demonstrated solution to the spending regression.
- Kappa often lowers predicted cost partly by reducing shear spending.
- Optimizing each subject's own tau times a share of cost tends to buy HMC almost every iteration. A cheapest non-regressing subject policy changes spending much less.
- Subject conditional-frame coordinates also change the replay decisions only modestly; the prototype approximates instantaneous jump factors and does not implement the wall quantile maps.
- Joint covariance predictors improve PM lag-10 forecasts, but underpredict slow tau on several continuation arms and overpredict during T30 burn. Quadratic features and ordered family-kernel composition do not fix the cross-fixture failure.
- A joint family operator gives useful slow-direction attribution, including approximately 70% shear share on the instrumented PM arm, but this does not validate its tau or spending recommendations.
- The concurrently edited policy was captured separately as `latest_allocation.R` and replayed read-only. Predicted costs change by -5.3% to +4.5% against point across the five states, and by -1.7% at the original checkpoint. This is not evidence of a realized speedup or a validation of subsequent edits.

There were no new continuation fits, no claimed speedup, and no commits. Fixed HMC zero/one controls and three-seed end-to-end validation remain outstanding. The scripts and CSVs in `eval_operator/` preserve the rejected candidates and their results; none is wired into the sampler.
