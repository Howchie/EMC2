# Handoff: DMC2 overhead work

## Status

The work produced real, measured speedups, but it did not meet the broader goal of closing EMC2's EXG and LNR fit-time gaps to DMC2. The earlier claim that the optimization objective was complete was premature. Keep the remaining EXG kernel gap and the ordinary-fit LNR system overhead open.

## EXG: truncation does not explain the reported gap

The often-quoted 150 ns versus 33 ns per row-particle result is the original, pre-Step-4 optimization probe, not the current EMC2 kernel. Step 4a/4b subsequently reduced the EMC2 704-row, 256-particle raw call to about 11 ms, or 61 ns per row-particle. The archived DMC2 probe was about 6 ms, or 33 ns per row-particle. That suggests a residual gap around 1.8x, not 4.5–5x. It is still an exploratory comparison, not a matched kernel benchmark.

The probe uses `calc_ll_oo` in EMC2 and `calc_ll` in DMC2. Each package constructs and compresses its own fixture (704 versus 705 rows), and the harness times the complete native likelihood call including its setup. The REXG models also differ: EMC2 normalizes a lower-truncated distribution; DMC2 uses the untruncated ex-Gaussian. Compiler flags differ as well. Dividing by row-particle count does not remove those differences.

For the probe parameters (`mu = 0.4`, `sigma = 0.1`, `tau = 0.2`, lower bound zero), the ex-Gaussian probability below zero is approximately `3.2e-6`, so the truncation normalizer is about `0.9999968`. The density normalization itself changes the log density by only about three parts per million. The normalizer is cached per exact parameter tuple, not recomputed for each row. Truncation therefore cannot plausibly account for the measured 4.5–5x runtime gap on this fixture. That conclusion is an inference from the fixture parameters and the distribution formula; the benchmark did not separately time the normalization term.

The remaining difference has not been isolated. The archived DMC2 race kernel uses separate winner/loser gather, contiguous compute and scatter loops, a compile-time `LOG_P` specialization, and SIMD directives. EMC2's common one-tuple path prepares parameters and computes directly in its row loop. The true gather/compute/scatter SIMD attempt did generate vector instructions but regressed both kernel and fit time, so it was reverted. The current EMC2 prepared loop did not show compiler vectorization of the density/survivor work. These facts point to kernel structure, dispatch, and vectorization as better places to investigate than the truncation normalizer, but do not establish which one causes the residual 1.8x.

## What improved

- The repeated REXG truncation calculation cache cut a 704-row raw call from 28 ms to 14 ms (2x). A fixed-transition EXG fit fell from 16.53 s to 11.26 s (1.47x), with identical samples and RNG.
- Stable-tail and one-tuple preparation changes reduced the raw call further (12.5 ms to 11 ms); a six-seed fit gate measured 1.097x speedup (95% paired interval 1.052–1.144).
- A finite, nonuniform scalar REXG race route ran 1.35–1.44x faster in its fit fixture.
- PSIS warming reduced default EXG fit time about 24% and LNR about 3.7%, with identical draws/RNG.
- Conditional subject-fit reuse reduced EXG time 13.5% and passed its ESS gate, so it is enabled for REXG. LNR reuse was 8.1% faster but failed to establish ESS non-inferiority; it remains opt-in.
- Block-heavy persistent-worker fixtures had large gains, including 2x and 4.71x cases. Those gains did not carry to ordinary one-core-per-chain whole fits.

See `RESULTS.md` for complete conditions and raw result locations.

## Why the overall target was missed

1. The cross-package 4.5–5x microbenchmark was treated too much like an explanation of avoidable EMC2 work. It compared different model semantics, package entry points, compressed row counts and build flags. It motivated useful same-package optimization, but did not identify a single 5x removable cost.
2. Same-package EXG optimization removed demonstrably repeated work and produced fit-level wins, but stopped with a measurable raw-kernel gap. The SIMD attempt that might have closed more of it regressed and was reverted; no successful specialized/vectorized replacement followed.
3. LNR's kernel was already comparable in the raw probe (about 56 ns versus 54 ns), while the exploratory shared-data whole-fit reference remained roughly 105 s versus 55 s. The planned process/worker changes were not shown to remove that ordinary-fit overhead: the persistent-chain default-fit gate was 1.006x for LNR (95% interval 0.969–1.043), with higher system CPU. Fewer minor faults alone are not a fit-time win.
4. Several measured gains were on narrow cached or block-heavy fixtures. They validate those paths, not the ordinary LNR workload. Treating these gates as completion of the broader speed objective was a judgment error.

## Next work

1. Build a matched EXG kernel benchmark from one shared `rt`, parameter, mask, censoring and design buffer. Compare EMC2's untruncated density path (`dexg`) with DMC2's untruncated path first, then add EMC2 truncation as a separate treatment. Compare per-row outputs as well as timing. Time kernel-only and full-call paths separately.
2. Run both kernels with the same compiler/toolchain settings where possible. Capture cycles, instructions, branches, vectorization diagnostics and an annotated profile. The existing row-normalized wall-clock probe is too coarse to assign the residual cost.
3. Use that profile to target the actual residual (likely the per-row compute/control path or dispatch). Preserve EMC2's normalized density, stable log behavior, masks and likelihood floor; do not adopt DMC2's untruncated distribution to obtain a timing win.
4. For LNR, profile a normal full fit by phase and process, then demonstrate that a system change reduces that measured cost and improves whole-fit wall time. The block-heavy persistent-worker results are not a substitute. Re-run posterior/ESS checks for any change that alters mode reuse or transition scheduling.
5. Only call the performance objective met after the matched raw-kernel cause and ordinary LNR whole-fit gap have direct evidence of improvement.

## Pointers

- `src/model_EXG.cpp`: tuple discovery/cache and REXG raw density path.
- `src/exgaussian_functions.h`: stable ex-Gaussian density, CDF and truncation formulas.
- `benchmarks/dmc2_overhead/kernel_archive.R`: original package-specific kernel probe.
- `benchmarks/dmc2_overhead/RESULTS.md`: measured before/after gates, including reverted SIMD and the ordinary-fit persistent-chain gate.
- Archived DMC2 reference: commit `fe8341a5`, `src/model_exgaussian.cpp` and `src/exgaussian_functions.h`.
