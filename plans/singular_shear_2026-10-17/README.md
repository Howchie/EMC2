# Singular shear adaptation investigation

Investigated the 726 `solve(): system is singular` warnings in
`/data/work/PM/NirvanaHons_Nback/Fits/lba7r_blk_gL_sepA_sepL_sepV_svL_mt.log`.
The warnings occur during adaptation between the first and second blocks.
The original adaptation states are absent from the retained checkpoint, so
the original coordinate and probe cannot be identified retrospectively.

## Mechanism and repair

`emc_shear_frames_cpp()` used nested Armadillo triangular solves with the
default approximate fallback. At extreme group SDs, the fallback can drop
directions and return finite, incorrect frame centres with `ok = TRUE`.
That is insufficient for the shear transport and its Jacobian.

The width search also allowed an arbitrarily large expansion from a nearly
flat quadratic fit and returned its last width even when the grid had not
resolved the density. A failed numerical calculation was either an immediate
failure or was treated as zero posterior density.

Changes:

- Use checked triangular solves with `no_approx`; a failed or nonfinite
  frame returns `ok = FALSE`. Ordinary solves are unchanged.
- Reject nonfinite group states before inversion/frame construction.
- Limit each grid expansion to four times its previous width.
- Retry a numerical grid failure at one quarter the width. Preserve actual
  negative-infinite log densities as zero probability, rather than treating
  them as numerical errors.
- Allow eight grid evaluations to accommodate gradual expansion and retry.
  An unresolved search returns `NA`, leaving the existing tuning unchanged.

No likelihood, prior, design, or model constraint was changed. The running
fit and default installed package were not modified. A separate build is at
`/tmp/emc2-shear-check` (temporary storage).

## Checkpoint replay

`replay.R` runs against either installation. The snapshot used here was
retained separately at `/tmp/emc2-shear-checkpoint.rds`; it contains sampling
states, not the original adaptation states. It is not included in the repo.

```sh
OPENBLAS_NUM_THREADS=1 Rscript plans/singular_shear_2026-10-17/replay.R /tmp/emc2-shear-checkpoint.rds
OPENBLAS_NUM_THREADS=1 R_LIBS=/tmp/emc2-shear-check Rscript plans/singular_shear_2026-10-17/replay.R /tmp/emc2-shear-checkpoint.rds
```

For a controlled extreme probe, multiply the `sv_LoadLow` group SD by
`exp(-48)` in a copy of chain 1's last covariance. The original implementation
emits 309 singular-solve warnings and returns `ok = TRUE`. Its first subject's
frame centre differs from direct triangular substitution by 9349.206 in
maximum absolute difference. The patched version emits zero warnings and
returns `ok = FALSE`.

Normal width searches on that same checkpoint agree to the printed precision:

| Group log-SD coordinate | Original width, window SDs | Patched width |
|---|---:|---:|
| `v.s_low` | 0.7725887 | 0.7725887 |
| `v.g_lure` | 1.131110 | 1.131110 |
| `sv_LoadLow` | 0.6052008 | 0.6052008 |

Starting the `v.s_low` search at 100 window SDs instead, the original returns
`NA`. The patched search narrows its grid and resolves a width of 0.7781363,
within 1% of the ordinary estimate. This stress test emits no solve warnings
with the patch.

These results establish a numerical failure mechanism consistent with the
log. They do not establish which parameter caused the historical burst or
prove that the statistical model is well identified or converged.

## Validation command

All three suites below passed, including the development posterior-invariance
checks. One existing queued-fit integration test in `test-shear-move.R` was
skipped by its `skip_on_cran()` guard. `git diff --check` also passed.

The targeted suites use the isolated installed package, with development
model validation enabled. Files run separately to avoid the repository's
parallel test runner starting workers against another package installation.

```sh
OPENBLAS_NUM_THREADS=1 R_LIBS=/tmp/emc2-shear-check EMC2_TEST_LEVEL=full Rscript - <<'RS'
library(EMC2)
library(testthat)
for (name in c("shear-move", "shear-standard-invariance", "collapsed-move")) {
  test_dir("tests/testthat", filter = paste0("^", name, "$"),
           reporter = "summary", package = "EMC2", load_package = "installed",
           stop_on_failure = TRUE)
}
RS
```
