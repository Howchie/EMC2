# Race-model numerical-stability implementation plan

## Goal

Add numerical protection to the remaining race-model likelihoods without making
the entire likelihood log-space-only.

The implementation must keep two independent concepts separate:

1. **Output scale** — controlled by `log_out`. This determines whether the
   public/helper function returns a log value or a natural-scale value.
2. **Calculation scale** — selected internally. The normal case should use the
   existing natural-scale calculation for speed. A log-space calculation should
   be selected only when the natural calculation is unsafe because of overflow,
   underflow, tail saturation, or cancellation.

The final conversion is therefore allowed in either direction:

```text
natural fast path -> return natural value or log(natural value)
log fallback      -> return log value or exp(log value)
```

The fallback should not be selected merely because `log_out = true`; a log
output request and a log-space calculation are separate decisions.

This plan excludes the legacy normal RDM implementation and its `pigt`/`digt`
(`*ig`) formulas, as requested. `RDMSWTN` is the supported Wald-form
replacement.

## Current gaps

The recent LBA/BAwL implementation is the intended pattern: guarded natural
evaluation followed by a log-space fallback. The relevant shared wrappers are
in [src/model_LBA.h](/data/work/EMC2_dev_oo/src/model_LBA.h:509), with targeted
tests in [tests/testthat/test-lba-logspace.R](/data/work/EMC2_dev_oo/tests/testthat/test-lba-logspace.R:8).

The remaining gaps are:

| Model/path | Current problem | Priority |
|---|---|---:|
| RDMSWTN, `sv = 0`, no kill | `dwald_k0()` is evaluated naturally and then logged in [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:1041) | High |
| RDMSWTN, `sv = 0`, no kill survivor | `pwald_k0()` is evaluated naturally and passed through `log1p(-cdf)` in [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:1057) | High |
| RDMSWTN positive-drift quadrature | Natural PDF/CDF quadrature is logged only after integration in [src/model_RDM.h](/data/work/EMC2_dev_oo/src/model_RDM.h:1428) | High |
| RDMGBM SPV density | `exp(exp_factor_log) * integral_result` is formed before logging in [src/model_RDM.h](/data/work/EMC2_dev_oo/src/model_RDM.h:744) | High |
| LNR density | Current raw path has regressed to natural density -> `log()` in [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:534) | High |
| LNR survivor | Fast normal log tails become `-Inf` above `z = 37` in [src/wald_functions.h](/data/work/EMC2_dev_oo/src/wald_functions.h:51) | High |
| RGAMMA density/survivor | Natural `dgamma`/`pgamma` are used before taking logs in [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:600) | High |
| REXG density/survivor | Natural density round-trip, plus upper-tail saturation in [src/exgaussian_functions.h](/data/work/EMC2_dev_oo/src/exgaussian_functions.h:123) | High |
| Generic GSL/truncation paths | Integrands multiply natural PDFs and survivors in [src/particle_ll.cpp](/data/work/EMC2_dev_oo/src/particle_ll.cpp:3345) | Medium |
| SSEXG | Log APIs exist, but truncated Ex-Gaussian survivor paths still use natural probabilities | Medium |
| SSRDEX | Uses excluded legacy `*ig` formulas | Out of scope |

## Design rules

### 1. Preserve the existing public contracts

Every public/helper function must continue to honour `log_out` exactly as it
does now. Internal helpers should not infer calculation scale from `log_out`.

For each model, provide logically separate evaluators:

```cpp
// Fast natural evaluation. Returns false when the result is unsafe.
bool model_pdf_natural_fast(..., double& out);

// Authoritative log evaluation for exceptional cases.
double model_log_pdf_stable(...);

// Output adapter: calculation fallback and final output conversion.
double model_pdf(..., bool log_out);
```

The output adapter should use the natural result when the fast evaluator accepts
it, otherwise call the stable log evaluator and convert only at the boundary.
For log-returning race kernels, the natural fast result is converted with
`log(out)` only after it has passed finite/positive checks.

### 2. Use explicit acceptance guards

Natural evaluation should be accepted only when all of the following hold:

- the result is finite and has the expected sign;
- probabilities are in their valid range;
- a density is not close to overflowing before the final requested conversion;
- a probability difference has not suffered cancellation;
- a survivor has not saturated to zero when a finite log survivor is needed;
- all intermediate products and sums remain finite.

The guards should be cheap and model-specific. They should not call the log
implementation in ordinary central cases.

Centralize shared constants and predicates, including `log(DBL_MAX)`,
`log(DBL_MIN)`, probability saturation tolerances, and relative cancellation
thresholds. Distinguish a legitimate log-tail `-Inf` from an invalid/failed
calculation so callers do not incorrectly turn every tail into a generic
minimum likelihood.

### 3. Keep natural integrands natural when safe

GSL and Gauss-Legendre APIs consume natural-scale integrands. Their normal path
should remain natural-scale. If a natural integrand is unsafe, use a protected
fallback rather than converting every evaluation to logs.

For products such as

```text
winner PDF × loser survivors
```

the fallback should form the product in log space. If the integration API still
requires a natural value, integrate a shifted integrand:

```text
g_scaled(t) = exp(log_g(t) - log_scale)
log_integral = log_scale + log(integral(g_scaled))
```

The scale can be selected from a cheap pilot evaluation or a fixed adaptive
quadrature pre-pass. This fallback is only entered after the natural integrand
fails its guard.

### 4. Do not use `log(exp(x))` or `log(1 - p)` as stable implementations

Where a log result is available, call it directly. In particular:

- use direct log densities instead of natural density followed by `log()`;
- use direct upper-tail log probabilities instead of `log1p(-cdf)` after a
  lower-tail CDF has rounded to one;
- use `log_diff_exp` or signed-log arithmetic for differences of close terms;
- use `log1m_exp` only when its input is already a valid log probability.

## Shared numerical primitives

### A. Robust normal log tails

Move the LBA-specific extended tail logic into a shared helper near
`pnorm_std()` in [src/wald_functions.h](/data/work/EMC2_dev_oo/src/wald_functions.h:64).

The shared helper should:

- use the fast approximation in the central region;
- use the existing direct tail approximation for ordinary tails;
- continue with a Mills-ratio/continued-fraction approximation beyond the
  current `z > 37` cutoff instead of returning `-Inf` solely because the natural
  probability underflows;
- preserve R semantics for invalid inputs and non-log calls.

Update LBA/BAwL to use the shared helper so there is one tail implementation.
Add tests showing that ordinary results are unchanged and that log tails remain
finite beyond the natural underflow range.

### B. Stable natural/log adapters

Add small reusable adapters for the common pattern:

```cpp
natural fast calculation
  -> accept if safe
log calculation
  -> convert to requested output scale
```

The adapters should not allocate and should be usable by raw kernels, scalar
helpers, and integration fallbacks.

### C. Stable upper-tail adapters

Provide a model callback or helper for direct log survivor evaluation. The race
combiner should be able to request a survivor log without reconstructing it from
a saturated natural CDF.

## Model-by-model implementation

### 1. RDMSWTN/Wald

Files: [src/wald_functions.h](/data/work/EMC2_dev_oo/src/wald_functions.h:211),
[src/model_RDM.h](/data/work/EMC2_dev_oo/src/model_RDM.h:301),
[src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:1041).

#### 1.1 Stable k=0 Wald primitives

Add log-space counterparts for the supported Wald-form k=0 calculations:

- point-start density;
- uniform start-point density;
- point-start CDF;
- uniform start-point CDF/survivor.

These must be implemented from the Wald-form expressions, not by routing back
through the excluded legacy `digt`/`pigt` functions.

For `A = 0`, use direct log terms. For `A > 0`, use log normal intervals and
signed-log subtraction for the start-point integral. Add a midpoint/degenerate
range fallback when the interval is too small to resolve safely.

Then implement natural wrappers with acceptance checks. Update
`rdmswtn_k0_logpdf()` and `rdmswtn_k0_logsurv()` to use the stable log primitives
when the natural calculation is rejected.

#### 1.2 Other Wald-form branches

Audit and guard:

- `dwald()`/`pwald()` guess paths that currently reconstruct survivors from
  natural CDFs;
- `positive_trunc_swtn_density_k0()` and
  `positive_trunc_swtn_cdf_k0()` natural differences;
- `drdmswtn_joint_A_sv_density_*()` probability/density differences;
- positive-drift quadrature paths that currently integrate natural values and
  only then take `log()`.

The quadrature should have two modes selected by the node guards:

1. natural-node accumulation for normal cases;
2. log-scaled node accumulation if any node overflows, underflows too early, or
   produces an unsafe cancellation.

Keep the existing `log_out` behavior independent from this mode selection.

#### 1.3 Acceptance tests

Add tests for:

- `sv = 0`, `A = 0` and `A > 0`;
- very high-density early-time cases;
- far-upper-tail survivors;
- positive-drift quadrature fallback;
- guess and kill-clock combinations;
- agreement between natural fast results and `exp(log result)` in the safe
  region.

Keep the existing reduction tests in
[tests/testthat/test-rdmswtn-k.R](/data/work/EMC2_dev_oo/tests/testthat/test-rdmswtn-k.R:1).

### 2. RDMGBM

Files: [src/model_RDM.h](/data/work/EMC2_dev_oo/src/model_RDM.h:699) and
[src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:356).

#### 2.1 SPV density

Split the current SPV density into:

- a natural fast evaluator for the usual parameter range;
- a log evaluator using `Gstar(..., true)`, log Gaussian interval integrals,
  and signed-log arithmetic for the two-term numerator.

The natural path must reject non-finite `exp_factor`, non-finite intermediate
terms, sign failures, and cancellation. The log path must avoid materializing
`exp(exp_factor_log)` until the final output conversion.

#### 2.2 Guess/local-combination paths

Replace natural CDF-to-survivor reconstruction in the guess path with a guarded
survivor evaluator. Update `dgbm_local_combo()` so it can combine log terms when
the natural decision density or survivor is unsafe, while retaining its current
natural combination in ordinary cases.

Review the `t = Inf` Erlang-2 expressions for the same issue: if individual
terms can overflow even though the final probability is finite, calculate the
terms in log/signed-log form.

### 3. LNR

Files: [src/wald_functions.h](/data/work/EMC2_dev_oo/src/wald_functions.h:88) and
[src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:522).

#### 3.1 Density

Use a guarded natural density fast path. When it is unsafe, call
`dlnorm_std(..., true)` directly. The raw race kernel should never perform the
natural density -> `log()` round trip currently present in `dlnr_raw()`.

#### 3.2 Survivor

Keep the natural CDF fast path for scalar consumers, but add a direct log
survivor fallback using the shared robust normal-tail helper. Use that fallback
when the CDF is close enough to one that `1 - CDF` is unreliable, and in the
batch raw/log-survivor path when the normal score is beyond the fast tail range.

Update `lnr_logS_at_t()` and `plnr_raw()` to use the same protected survivor
implementation.

### 4. RGAMMA

Files: [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:581).

#### 4.1 Density

Use natural `R::dgamma(..., false)` only when it is finite and positive. On
failure, use `R::dgamma(..., true)` and convert at the boundary according to
`log_out`.

#### 4.2 Survivor

Use natural `pgamma(..., lower_tail = true, log_p = false)` for the normal case.
When it is near one or has saturated, use
`R::pgamma(..., lower_tail = false, log_p = true)` directly. Apply this in:

- `prgamma_raw()`;
- `rgamma_logS_at_t()`;
- scalar survivor helpers when their result is used in an extreme truncation or
  integration calculation.

### 5. REXG

Files: [src/exgaussian_functions.h](/data/work/EMC2_dev_oo/src/exgaussian_functions.h:16)
and [src/utils.h](/data/work/EMC2_dev_oo/src/utils.h:670).

#### 5.1 Density

`dexg()` already derives the density in log terms. Expose that calculation as
the fallback for `drexg_raw()` instead of evaluating `dexg(..., false)` and
logging afterward.

#### 5.2 Direct upper-tail log probability

Add a true upper-tail log implementation for `pexg()`. It must not:

1. exponentiate the lower-tail log CDF;
2. round it to one;
3. apply `log1m()`.

Use a natural fast upper-tail path when it is safe, then switch to the direct
log upper-tail expression near saturation. Thread this through `prexg_raw()`
and `rexg_logS_at_t()`.

#### 5.3 Truncated Ex-Gaussian paths

Audit `dtexg()`/`ptexg()` with finite truncation bounds. Their normalizer and
tail differences should use log differences when the natural subtraction fails,
while retaining the existing natural path for ordinary bounds.

### 6. Stop-signal models

Treat SSEXG as a follow-on application of the Ex-Gaussian work:

- use the protected truncated Ex-Gaussian density and survivor callbacks;
- retain natural stop-success integrands when safe;
- add the scaled-log integrand fallback if a stop-success product is unsafe.

Leave SSRDEX unchanged in this work because it depends on the explicitly
excluded legacy `*ig` formulas.

## Race-combiner and integration work

### 1. Finite known-RT path

The batch raw path already expects log density/log survivor outputs. After the
model fixes, verify that every supported model supplies those values without
natural round-trips.

Add a consistency check for each raw callback:

- finite central-case output;
- valid `-Inf` for a true zero probability;
- no `NaN` from a failed natural calculation;
- no accidental replacement of a finite extreme log value with `min_ll`.

### 2. Scalar/truncation/GSL path

Keep the current natural path in [src/particle_ll.cpp](/data/work/EMC2_dev_oo/src/particle_ll.cpp:3345)
for speed. Add a fallback only when:

- the winner PDF is non-finite or non-positive;
- a loser CDF is outside its valid range;
- a survivor has saturated to zero;
- the product becomes zero, infinite, or non-finite.

The fallback should evaluate the complete integrand as a log sum and use a
shifted natural integrand for GSL, rather than changing every integrand call to
log space.

Apply the same pattern to:

- race truncation normalizers;
- BAwL/RDMSWTN adaptive integrals;
- stop-success integrals in `ss_raw.h`;
- any remaining natural products in local guess/kill combinations.

## Testing plan

### Unit-level tests

For every model, test both output scales independently from calculation-space
selection:

- safe central cases: natural result equals `exp(log result)`;
- density overflow cases: log result remains finite and natural output follows
  the documented representability behavior;
- density underflow cases: natural result may be zero, while log result remains
  finite;
- upper-tail survivor cases: log survivor remains finite after natural CDF
  saturation would have produced zero;
- cancellation cases: no `NaN`, sign inversion, or spurious minimum likelihood.

Add model-specific cases for LNR, RGAMMA, REXG, RDMGBM, and RDMSWTN. Extend the
existing LBA/BAwL tail tests with shared normal-tail cases so all models test
the same `z > 37` behavior.

### Likelihood-level tests

Construct small races with:

- one extreme winner and ordinary losers;
- one ordinary winner and an extreme loser survivor;
- truncation bounds in both extreme tails;
- mixed guess/kill-clock paths;
- natural-scale integration requests.

Verify that the likelihood remains finite when the mathematically relevant log
terms are representable, and that ordinary likelihood values remain unchanged.

### Performance/regression tests

Benchmark representative central cases before and after the changes. The
expected result is:

- no material slowdown in ordinary natural-scale cases;
- extra cost only when an acceptance guard rejects the natural path;
- no unconditional calls to expensive log-space formulas in hot loops.

Add optional counters or a debug-only build hook for fallback frequency during
development, then remove or disable them in the release path.

## Recommended implementation order

1. Add shared normal-tail and natural/log adapter primitives.
2. Fix LNR, RGAMMA, and REXG raw density/survivor kernels.
3. Add Wald-form k=0 log primitives and repair RDMSWTN raw shortcuts.
4. Add guarded RDMSWTN positive-drift and joint `A`/`sv` fallbacks.
5. Repair the RDMGBM SPV density and local-combination paths.
6. Add scalar/GSL shifted-log fallbacks for race truncation and integration.
7. Apply the Ex-Gaussian protections to SSEXG.
8. Add extreme-case tests, run the existing race test suite, and benchmark the
   natural fast paths.

## Acceptance criteria

- `log_out` remains purely an output-scale choice.
- Natural calculations remain the default in ordinary parameter regions.
- Every supported model has a guarded log fallback for density overflow,
  survivor saturation, and unstable probability differences.
- Raw race kernels return stable log density/survivor values without natural
  round-trips.
- Natural-scale integrands still work and use log-space only when their natural
  calculation is rejected.
- Existing central-case values and benchmark behavior remain within agreed
  tolerances.
- Normal RDM and legacy `*ig` formulas remain out of scope.
