# Uniform guess contaminant (`pGuess`) — implementation plan

Status: **plan only, nothing implemented.** Written 2026-08-04.

## Goal

Add a second, RT-contributing contaminant process to the race and DDM families.
The existing `pContaminant` is a Bernoulli *omission* mixture: it scales the
likelihood by a failure probability and adds mass only at `rt == +Inf`. It does
nothing for observed RTs. The new `pGuess` adds the standard uniform-outlier
mixture (Ratcliff & Tuerlinckx 2002; HDDM's `w_outlier`), which contributes
directly to observed RT densities and gives fast/slow outliers a likelihood
floor.

`pContaminant` must not change behaviour. At `pGuess = 0` every likelihood is
bit-for-bit identical to the current implementation.

## Current state (traced)

### `pContaminant` is a post-hoc per-trial mixture

Applied *after* the trial likelihood is complete, including after truncation
renormalisation:

- R reference: `R/likelihood.R:366-376` —
  `p[isMiss] <- pc + (1-pc)*p; p[!isMiss] <- (1-pc)*p`, where `isMiss` is
  `rt == +Inf` only (a left-censored `-Inf` or missing `NA` rt does **not**
  pick up omission mass).
- C++ hot path: `src/particle_ll.cpp:8648`, the `apply_pC` lambda — deliberately
  factored so the expand and compressed branches cannot drift apart.
- Six further C++ sites replicate the same shift:
  - `:4485` race all-finite raw fast path
  - `:7449-7470` `c_log_likelihood_logicalrules`
  - `:9811, :9877` LogicalRules fast-const path (`fast_log1m_pc`)
  - `:9283` BAwLcorr
  - `:10747` RDMSWTNcorr

### Parameter declaration pattern

`pContaminant` is a trailing p_type with default `qnorm(0)`, transform `pnorm`,
bound `c(0.001, 0.999)`, exception `0`. Two mechanisms make a trailing nuisance
parameter inert:

- `design()` (`R/design.R:357-366`) turns any p_type absent from the formula
  into a constant at its default, and suppresses the "assumed constant" message
  for names outside `p_types_canonical`.
- `emc2col::validate_col_prefix` (`src/col_registry.h`) checks only the
  canonical *prefix* of the parameter column order, so trailing columns are free.

This is the extension pattern to reuse.

### Truncation / censoring windows

`LT`, `UT`, `LC`, `UC` are dadm columns (`add_bound_column_if_needed`,
`R/design.R:900`), part of the compression key (`compress_dadm`,
`R/design.R:699-712`), and read C++-side via `get_col_with_default`.

### DDM

`DDM()` / `DDMGNG()` (`R/model_DDM.R:145`) have **no** `pContaminant` at all.
The C++ kernel `c_log_likelihood_DDM_pt` (`:2019`) has a full
truncation + censoring path (`:2083` onwards) and an all-finite fast path
(`:2047`), but no `pc_col`.

### Simulation

`make_missing()` (`R/make_data.R:195-207`) applies the omission *after* the
truncation cut and *before* censoring. `make_data()` auto-populates
`TC$pContaminant` from the pars matrix at `R/make_data.R:571`.

## Design decisions

### 1. No rename

Keep `pContaminant` as the omission parameter. Renaming it breaks saved emc
objects, `TC` lists, direct `make_missing()` calls and ~15 man pages for no
functional gain. Add `pGuess` alongside it. `pC_Omission` / `pC_Guess` may be
accepted as aliases rewritten to the canonical names at the top of `design()`
if the symmetric naming is wanted; that is one line and independent of the rest.

### 2. Nested (stick-breaking) weights, not a simplex

```
P(omission) = pC
P(guess)    = (1 - pC) * pG
P(process)  = (1 - pC) * (1 - pG)
```

Both remain independent `pnorm` parameters on [0, 1] and can never sum above 1.
At `pG = 0` the arithmetic reduces exactly to the current expressions — this is
the property that protects `pContaminant`.

### 3. One scalar guess window per dadm, stored as an attribute

No new dadm columns and no per-row windows. `resolve_guess_window()` runs once
in `design_model()`; the result is `attr(dadm, "guess_window")`, a length-2
numeric.

Resolution order:

1. `TC$guess_window` if supplied by the user — wins outright. This is the
   exposed knob, and normally would not be touched.
2. Otherwise `LG = max(LT, LC)`, `UG = min(UC, UT)` from the existing columns.
3. If `UG` is still infinite: `UG = max(5, floor(max_finite_rt) + 1)` — 5 s
   unless the data run past it, in which case the next whole second strictly
   above the largest observed RT.

If the bound columns are not constant across rows (possible with subject-wise
bounds from `make_missing`), take the union — `min` of the lowers, `max` of the
uppers — and emit a message.

Rationale for the fallback: HDDM's fixed 5 s window is improper by their own
admission ("in practice, the outlier model is applied to all RTs, even those
larger than 5"). Keeping the density exactly 0 above `UG` while defaulting `UG`
to 5 would leave a 7 s outlier with no mixture protection at all — precisely the
trial the model exists to catch. Scaling up to the next integer above `max(rt)`
stays proper and covers every observed RT. For data that fit inside 5 s the
behaviour is identical to HDDM.

### 4. Expose the window, not `w_outlier`

`w_outlier = 0.1` only means "5 s" because HDDM has exactly two responses; on a
four-accumulator race the same number silently implies a 2.5 s window. Store
`guess_window`; report the implied `w_outlier = 1 / (n_resp * W)` in the design
print so the HDDM correspondence stays visible. `w_outlier` may be accepted as
an alternative spelling that converts on entry.

### 5. Guess kernel

Uniform in RT on `[LG, UG]`, uniform over responses:

```
log_g = -log(n_resp * (UG - LG))        // one double, computed once per dadm
```

`n_resp` excludes the `nogo` and `time` pseudo-levels — a guess is an overt
response, so it cannot be a withheld/nogo trial and cannot be a timeout.

Because the window is `[max(LT, LC), min(UC, UT)]` by construction, a guess can
never be censored or truncated away. That removes all interval-mass machinery:
no `log_guess_mass`, no censored-branch guess terms.

### 6. Per-trial mixture — two cases only

| trial | mixture |
|---|---|
| finite rt, R known | `log[(1-pC) * ((1-pG) * L_proc + pG * exp(log_g))]` |
| everything else (`+Inf`, `-Inf`, `NA`, withheld, nogo) | `log[pC * 1{rt == +Inf} + (1-pC) * (1-pG) * L_proc]` |

Proper: `pC + (1-pC) * [pG * 1 + (1-pG) * 1] = 1`.

The one branch that touches an unknown R with a finite rt
(`particle_ll.cpp:8608`, `log_min_density_rowmajor`) takes `log_g` **without**
the `n_resp` division — the same expression with one operation skipped, so the
density still integrates to 1 there. No R-unknown machinery is built beyond that.

### 7. Semantics: proportion among *retained* trials

Both components are renormalised on the truncation window before mixing, and
`make_missing` already contaminates after the truncation cut. So `pG` (like
`pContaminant`) is the guess proportion among **retained** trials, not among
generated trials. This must be stated in the docs — it is the kind of thing that
otherwise only surfaces as an SBC miscalibration.

### 8. No R likelihood changes

`log_likelihood_race_missing`, `log_likelihood_ddm` and `log_likelihood_ddmgng`
are no longer supported paths and are left untouched; `pContaminant` stays
exactly as it is in them.

Consequence: a model with no `c_name` dispatches to the R path
(`R/sampling.R:1094`) and would silently ignore `pGuess`, sampling a parameter
that does nothing. Add a guard in `design()`: if `pGuess` is free or non-zero
and `model()$c_name` is NULL, `stop()` with "pGuess requires a compiled
likelihood". One line; turns a silent wrong answer into a loud one.

### 9. No identifiability warning

`pGuess` and `pContaminant` serve different use cases and are not expected to be
used together. Using `pGuess` on data with no omissions is the intended use, not
a mistake, and must not warn.

## Change list

### C++

**New `src/contaminant_mixture.h`** — a single shared helper so the seven
application sites cannot drift apart (the existing code already carries a
comment about exactly this failure mode):

```cpp
inline double mix_contaminants(double ll_proc, double pC, double pG,
                               double log_g,        // -Inf on non-finite-rt trials
                               bool is_omission);
```

**`src/particle_ll.cpp`**

- Resolve `pg_col` beside the existing `pc_col` at `:4198`, `:4332`, `:4658`,
  `:7610`, `:8815`, `:10458`; add new resolution in the DDM branch at `:4045`.
- Read the scalar `log_g` once from `attr(dadm, "guess_window")` per likelihood
  call.
- Swap all application sites onto `mix_contaminants`:
  - `:4485` race all-finite fast path
  - `:7449-7470` LogicalRules kernel
  - `:8648` generic race `apply_pC`
  - `:9283` BAwLcorr
  - `:9811 / :9877` LogicalRules fast-const path
  - `:10747` RDMSWTNcorr
  - `:2047` (fast) and `:2355` (comprehensive) in the DDM kernel — new, since
    the DDM has no contaminant handling today.

### R

- **`R/utils.R`** — `add_nuisance_pars(p_types, transform, minmax, exception)`
  appending both `pContaminant` and `pGuess`. Refactor the model constructors
  onto it: `model_LBA.R`, `model_LNR.R`, `model_RDM.R`, `model_RLF.R`,
  `model_ROU.R`, `model_GOM.R`, `model_BAwD.R`, `model_BOU.R`.
- **`R/model_DDM.R:145`** — add both parameters to `DDM()` and `DDMGNG()`, and
  set `p_types_canonical` there (it currently has none) so the "assumed
  constant" message stays quiet.
- **`R/design.R`** — `resolve_guess_window()` called from `design_model()`,
  result attached as an attribute and carried through `compress_dadm()`'s
  attribute block (`:729-750`) so it survives contraction; plus the no-`c_name`
  guard from decision 8.
- **`R/make_data.R`** — `make_missing(..., pGuess = NULL)`: draw
  `Bernoulli((1-pC) * pG)` on non-omitted trials, then `rt ~ U(LG, UG)` and
  `R ~ sample(response levels minus nogo/time)`, placed beside the omission
  block at `:195` and before the censoring step. Add `pGuess` and
  `guess_window` to `check_missing()` defaults (`:245-268`) and to the pars
  auto-populate at `:571`.

Deliberately **not** changed: dadm columns, compression keys, R likelihood
functions, `model_rng.R` (keeping guess injection in `make_missing` makes it
model-agnostic, as the omission already is).

### Docs

`man/` regeneration for every affected model; `NEWS.md`; explicit note that
`pGuess` is the observed-RT contaminant and `pContaminant` remains the omission
one, together with the retained-trials semantics of decision 7.

## Tests

The R path is no longer an oracle, so:

1. **Regression, the one that protects `pContaminant`:** `pG = 0` reproduces the
   pre-change C++ likelihood bit-for-bit on every family. Compare against the
   currently built `EMC2.so`.
2. **Test-local reference:** a small R implementation of the two-case mixture,
   living in `tests/testthat/` and not shipped, checked against C++ on a dadm
   covering finite / `+Inf` / `-Inf` / `NA` / withheld trials.
3. **Propriety:** numeric quadrature confirming the mixture integrates to 1 over
   the window.
4. **HDDM parity:** two-response DDM with `guess_window = c(0, 5)` gives a
   per-response density of 0.1.
5. **Recovery:** MC recovery run with `pGuess` free.
6. Expand vs compressed agreement, mirroring
   `tests/testthat/test-pcontaminant-compressed.R`, which is the right template.

## Phasing

1. Helper header, `p_types` additions, `resolve_guess_window()`, and the DDM
   kernel. Smallest self-contained slice, and the one with the clean HDDM parity
   check.
2. Generic race kernel (`:8648`) and the race fast path (`:4485`).
3. The four specialised kernels: LogicalRules (both paths), BAwLcorr,
   RDMSWTNcorr.
4. Simulation (`make_missing`), recovery study, docs.

Estimate: ~half a day R, ~half a day C++, ~a day tests and docs. The risk is not
the mathematics but keeping the seven C++ sites in sync — hence the single
shared header.
