# Ensemble group update benchmarks

Scripts behind the 2026-09-27/28 evaluation of the ensemble group update
(R/ensemble_move.R) against the sampler without a group update. Outputs go to
`results/` (git-ignored). One fit at a time, 3 chains x 9 cores.

* `fixture_fits.R` -- simulated RDM fits from
  `WorkingTests/generate_correlated_group_rdm_fixture.R` and a wide
  24-parameter design, both calibrated so every group SD is at least ~1.4x its
  per-subject posterior SE (otherwise no arm converges and nothing is
  learned). Cases: diag, corr, strong, n80, freeA, corr_blocked, diag_dg,
  tight, wide, forstmann, corr_factor, corr_infnt, strong_long.
  `Rscript fixture_fits.R fit <lib> <case> <arm> <seed>`, then `summarize`
  (`FIXTURE_ROOT` selects the results folder).
The `baseline`, `ensemble` and `ensemble_marg` arms of that evaluation were
set by sampler options that no longer exist; `fixture_fits.R` now runs the
package's one sampler and `<arm>` only labels the output. The earlier interweaving
group move (legacy / AM / RAM) and its benchmarks were removed after commit
ac3e4d36; see that commit for them.
