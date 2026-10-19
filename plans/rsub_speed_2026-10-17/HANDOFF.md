# RSUB speed handoff (2026-10-17)

`RSUB(evidence = c("dense", "packets"))` was committed in ed8b48e9 on
`playground`. The math is in `Math/RDMSWTN_Derivation.tex`, section
"Subordinator Accumulators". This note records the model's cost compared with
LBA and RDM, and where that cost comes from.

The build tested here is installed only in `.lib/rigs`; the system library
does not have RSUB. Before running either script against another library, set
`EMC_LIB=<lib>`.

## Scripts

| Script | What it times |
|---|---|
| `bench_ll.R` | Full compiled likelihood (`calc_ll_manager`) for LBA, RDM, RSUB dense and RSUB packets. Uses the first forstmann subject, 100 particles, 15 interleaved repetitions, and reports medians. |
| `bench_kernel.R` | RSUB kernels alone (`drsub` / `prsub` on log scale, 1e5 rows) at `s` = 1, 0.3 and 0.1. |

Run both from the repo root, pinned to a single core when the machine is
loaded:

    taskset -c 5 Rscript plans/rsub_speed_2026-10-17/bench_ll.R
    Rscript plans/rsub_speed_2026-10-17/bench_kernel.R

## Results

These were measured with load average around 15, so the ratios are reliable
and the absolute times are noisy.

Full likelihood, in ns per accumulator-row:

| | LBA | RDM | RSUB dense | RSUB packets |
|---|---|---|---|---|
| A > 0 | 86 | 86 | 327 (3.8x) | 778 (9.0x) |
| A = 0 | 56 | 49 | 74 (1.3x) | 358 (6.4x) |

Each fit used `v` around 2 (3 for LBA), `B` = 1, `A` = 0.4 or 0, `t0` = 0.15,
and the scale parameter (`s` or `sv`) fixed at 1.

Kernels alone, in ns per row, with `v` = 1.3, `B` = 0.9, `A` = 0.5 or 0, and
`u` ~ U(0.1, 3):

| law | s | pdf A>0 | surv A>0 | pdf A=0 | surv A=0 |
|---|---|---|---|---|---|
| dense | 1 | ~500 | ~410 | ~95 | ~155 |
| dense | 0.3 | ~760 | ~590 | ~120 | ~240 |
| dense | 0.1 | ~510 | ~370 | ~105 | ~140 |
| packets | 1 | ~810 | ~880 | ~75 | ~500 |
| packets | 0.3 | ~1420 | ~1490 | ~85 | ~825 |
| packets | 0.1 | ~2660 | ~2700 | ~60 | ~1470 |

## Where the cost goes

**Dense, A = 0.** The density is one Mills-ratio evaluation, and the survivor
is `phi(alpha)[R(alpha) + R(beta)]`. Both are already near LBA/RDM speed.

**Dense, A > 0.** At each end of the start range the kernel evaluates both
complementary forms (J or 1-J, H or L) and keeps the one that does not cancel.
Each form needs Mills-ratio integrals: `log_mills_int` takes a direct
difference, or uses 8-point Gauss-Legendre when the two ends are close. A
further GL fallback over sqrt(d) runs only when the two ends agree to within a
few bits. Possible savings:

- Decide early versus late from `alpha` alone and skip computing the unused
  form.
- Share one `mills_seq` between the J and H calls at the same end, since the
  pdf and survivor are often requested together.

**Packets, any survivor or any A > 0.** These take a Skellam side sum
(`cp::side_sums`) of length about `sqrt(kstop^2 + 40 z) + 10`, with
`z = 2 sqrt(ab)`, `a = 2 v u / s^2` and `b = 2 d / s^2`. The cost therefore
grows like `1/s`. The Bessel ratios come from a backward recurrence over a
`thread_local` buffer. The second side is summed only when the first side
holds more than `1 - 1e-3` of the mass (`kSideFull`). The point density alone
is a closed-form `I0` (series below z = 25, Hankel above) and is cheap.

## Suggested next step if packets get used heavily

When `a + b` is large (more than a few hundred), replace the Skellam side sums
with a uniform asymptotic expansion. A saddlepoint (Lugannani-Rice) on the
Skellam cumulant generating function `K(t) = a(e^t - 1) + b(e^-t - 1)` gives
`P(Delta <= 0)` and, through the tilted means, `E[Delta^+]` and
`E[Delta^-]`. That would make the cost independent of `s`.

Validate any such change against the brute-force references in
`tests/testthat/test-rsub.R` (`rsub_pk_lF`, `rsub_pk_pdf`). Those references
must sum past the quantiles of both `max(a, b)` and `sqrt(ab)`. Also check
accuracy at about 1e-10 relative on log S out to log S of about -2000; the
current kernels reach about 1e-13 there.

## Checks to rerun after any kernel change

    NOT_CRAN=true R_LIBS=.lib/rigs Rscript -e 'library(testthat); suppressMessages(library(EMC2)); \
      source("tests/testthat/helper-rng.R"); test_file("tests/testthat/test-rsub.R")'

That file currently gives 353 passes. Also run `test-opr.R` (118 passes),
which shares the Hankel `I_nu` helper.

Reinstall the library before testing:

    R CMD INSTALL --no-docs -l .lib/rigs .

If only headers changed, delete `src/*.o` before reinstalling.

## Optimisation pass (2026-10-17, uncommitted, src/model_RSUB.h only)

Built into `.lib/rsubfast`; `.lib/rigs` still holds the ed8b48e9 build for A/B.
Full likelihood, ns per accumulator-row (old -> new): dense A>0 333 -> 86,
A=0 74 -> 43; packets A>0 790 -> 130, A=0 358 -> 105.  LBA is 86 / 56.

Dense:
- `mills_seq(beta)` formed once per end and passed into `mills_int` (was 3x).
- `H = L + (m - d)`: one Mills integral per end gives both stop-loss forms.
- `log_pdf` forms only the J form it will difference (`prefer_Jc`), trying
  the other only on cancellation; `log Phi(alpha)` replaced by `1 - Jc`.
- Linear-arithmetic fast path (`fast_pdf`, `fast_cdf_surv`) for |alpha| < 30;
  the log-scale path remains the fallback.
- `R(x)` for x < 4 from the Hart rational with the exponential cancelled;
  CF above 4 seeded at its fixed point with a finer depth schedule.
- `log_surv`/`log_cdf` skip quadratures and point CDFs only the other needs.

Packets:
- For min(a,b) <= 60, max <= 650: direct forward Poisson double sums
  (`near_direct`), both start-range ends in one interleaved pass, linear
  differencing.  No Bessel functions at all.
- Otherwise: division-free Miller recurrence with Horner accumulation and
  I_0 normalised by I_0 + 2 sum I_k = e^z (`side_sums`); old ratio form kept
  for z < 1e-30.
- `skel` forms only the quantities its caller needs.

Accuracy: 2e5-point random grid vs old build, max rel dlog 9.5e-12 (99.9%
<= 2.4e-13); 256-bit brute-force Skellam reference on 150 packets points:
pdf 8.3e-13, cdf 3.7e-13 (old 6.3e-14), surv 1.8e-15.  test-rsub 353/353,
test-opr 118/118.

Still slow: packets at small s (min(a,b) > 60) is Miller, O(sqrt(a+b))
steps (~1.6 us at s = 0.1).  The uniform asymptotic / saddlepoint idea above
is the next step there.  Dense at s <= 0.3 is bound by CF depth for
beta in [4, 8).

## Small-s pass (2026-10-17, commits 66fe81e0, 1e18a548)

- Dense: Mills Taylor table extended to [8, 32) (spacing 1/2), CF depth 9/7
  above.  Dense cost is now flat in s (kernel bench pdf_A 137/143/167 ns at
  s = 1/0.3/0.1); there is no small-s-specific headroom left.
- Packets, Miller regime (min mean > 60 or max > 650, z >= 40): `near_asym`,
  Hankel expansion of I_0/I_1 times the binomial series in t/sqrt(big),
  integrated against Gaussian-tail moments U_j from the Mills R, D, E
  (forward recurrence; its error is damped by (delta/sqrt(big))^j).
  Coefficients C_k binom(p_k, j) are tabulated, so each Hankel order is a dot
  product.  ~140 ns vs ~1000 ns for Miller, max error 2e-15 vs 256-bit
  (Miller was up to 2e-14 at z ~ 1e4).  Accurate for z >= 40 and big >= 40;
  below z = 40 Miller remains (cheap there).
- Direct sums templated on the sums needed (t0/t1 for densities, tE for
  stop-losses), per-sum ratio bounds, no r/(1-r) division, no out-of-line
  `std::fmin/fmax` (they are libm calls under -fno-finite-math-only).
- Kernel bench packets pdf_A/surv_A: s=1 197/203 -> 140/120, s=0.3
  460/473 -> 313/260, s=0.1 1593/~1600 -> 417/480.  Prior-like cost map:
  s fixed 235/230 -> 175/150, s sampled 325/355 -> 220/225 ns/row.
- `bench_kernel.R` reads `EMC_LIB` and defaults to `.lib/rigs` (the old
  build): set `EMC_LIB=.lib/rsubfast` or it silently times the old code.
