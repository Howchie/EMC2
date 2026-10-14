# Link declared parameter bounds to the sampling transform

## Goal

Make a parameter's declared bound (`bound$minmax`) part of the transform that
maps the sampled real line to the natural scale. The sampler then moves
smoothly over the whole support window and never reaches the bound. Values
reported to users must keep meaning what they say on the natural scale.

This is a modelling and interface change, independent of the sampler work on
`playground`. Scope it before building it, and expect the touch points below to
be incomplete.

## Why

Two mechanisms currently describe a parameter's support, and they are never
reconciled:

- **The transform** (`transform$func/lower/upper`, filled by `fill_transform`
  in `R/map.R`). For example, `exp` with lower `L` maps x to `L + exp(x)`, and
  `pnorm` maps x to `lower + (upper - lower) * pnorm(x)`. When a model gives no
  lower, `exp` defaults to 0.
- **The bound** (`bound$minmax` plus `bound$exception`, filled by `fill_bound`).
  It is checked strictly (`val > min && val < max`) on the natural scale, after
  mapping. Outside the bound, every trial in that design cell is floored at
  `min_ll` (`log(1e-10)`): a large but finite drop, not `-Inf`.

So for t0 in LBA, RDM, LNR and DDM (`exp`, no lower; bound `[0.05, Inf)`), the
sampled log t0 covers (0, Inf), but 0.05 is a hard cliff inside the
likelihood. Posterior mass piles up against it (seen on an N-back fit: about 17
of 93 subjects had LoadLow-target cells sitting on 0.05). Shear, collapsed,
HMC and particle moves all lose proposals to the cliff. The intended design is
that the transform keeps the sampler on the real line; for bounded parameters
we currently do not.

The terminology is confusing:
- The transform's `lower`/`upper` are really the mathematical support of the
  link.
- The bound is the support we actually want, and it must lie within the
  transform's window.

## Target behaviour

1. **Transform follows the bound.** When the transform family can represent a
   bounded side, that side of the declared bound becomes the transform's
   `lower`/`upper`: `exp` takes a finite lower bound, `pnorm` takes both. This
   must be general, not a t0 special case.
2. **User overrides.**
   - An explicit user `transform$lower`/`upper` wins. Validate that the declared
     bound lies within it.
   - A user change to the bound (`design(bound = ...)`) moves the transform with
     it, unless the user also set the transform explicitly.
3. **The strict bound test stays as a backstop** for everything the transform
   cannot guarantee (see the challenges below).
4. **Natural-scale meaning for users.** Anything a user reads or types should be
   on the natural scale by default, or explicitly labelled as the internal
   sampled scale.

## Challenges already identified

### It changes the model, not just the sampler

Priors are specified on the sampled scale, and the hierarchical group
distribution is normal on whatever scale is sampled.

- Today the group distribution is normal on log t0, truncated at log 0.05.
- With the shifted transform it is normal on log(t0 - 0.05). That is a
  different model, and its implied prior density goes to 0 at the bound.

Subjects whose posterior currently piles on the bound will get a visibly
different posterior. This is arguably more honest, but it must be documented.

### Sampled-scale numbers change meaning

`get_pars` defaults to `map = FALSE`, so summaries show sampled-scale values.
Priors, start values and `make_data` parameter vectors are also given on the
sampled scale. For example, -2.12 means t0 = 0.12 today but 0.17 with lower =
0.05. The user does not want to see shifted-log numbers. Plan for:

- Natural-scale reporting by default (or clear labelling of the raw scale).
- Natural-scale input for priors, start values and `make_data` values,
  converted under the hood.
- Default priors that keep a sensible natural-scale meaning on the new scale.

### A purely internal change of variables is not enough

A Jacobian change of variables would keep today's model, but only when the
bound applies to one sampled coordinate. Bounds apply per design cell, and a
cell is a linear combination of coefficients (for example `t0 ~ Load +
itemType`). The transform acts on each cell's linear predictor, so shifting the
transform does enforce per-cell bounds. A per-coefficient change of variables
cannot.

### Some sides cannot be absorbed by the transform family

Do not silently change a transform family (for example identity to pnorm): that
changes the geometry far more than the bound does. Keep the strict test on any
side the family cannot represent. Examples from `DDM()` (`R/model_DDM.R`):

- `a`: `exp`, bound `[0, 10]`. The upper bound of 10 cannot be absorbed.
- `st0`: `exp`, bound `[0, .5]`. The upper bound cannot be absorbed.
- `v`: `identity`, bound `[-20, 20]`. Neither side can be absorbed.
- `Z` and `SZ`: `pnorm`, bound `[.01, .99]`. Both sides can be absorbed.

### Exceptions sit outside the window

`bound$exception` permits exact values outside `minmax`, usually for parameters
fixed as constants:

- DDM: `t0 = 0`, `sv = 0`, `SZ = 0`, `st0 = 0`, `pGuess = 0`.
- LBA: `A = 0`.
- BTAwL: `pi = 0` or `1` (names may repeat; see `exception_hits`).

A constant given on the sampled scale passes through the transform, so a
shifted lower bound would turn `-Inf` into `L` instead of 0. Constants that hit
exception values must bypass the shifted transform, or exceptions must be
handled another way.

### Trends and other natural-scale paths can still cross the bound

- In the C++ likelihood (`src/particle_ll.cpp`, `prepare()`), the bound is
  checked after `update_pt_only`. Natural-scale trends are applied there
  (`TrendEngine` postmap specs, `src/TrendEngine.h`) and can push a parameter
  past its bound regardless of the transform.
- Bounds that depend on the data, such as t0 below the fastest RT, cannot be
  absorbed either. Those trials are floored by the kernel, which is correct
  model behaviour.

### Numerical floors

Many bounds are tiny numerical floors, for example 1e-4 for `exp` parameters
(LBA `sv`/`A`/`B`, BAwD, many others) and MRI's `sd` at 0.001. Shifting those
changes almost nothing, but it does change the sampled scale. Decide whether to
shift every finite bound or only those that matter: a uniform rule is simpler
and more predictable.

### Old fits and objects

Saved fits, priors and designs carry their filled transform. They must keep
mapping exactly as before (store the filled transform in the object; do not
re-derive it on load). Add a NEWS entry.

## Touch points (non-exhaustive)

Start here, but search the code base: other consumers of `bound`, `transform`
or the sampled scale will exist.

- **`R/map.R`**
  - `fill_transform`: where lower/upper are filled and defaulted, and where
    the link goes.
  - `fill_bound`.
  - `do_bound` (strict test).
  - `fix_bound` (`shrink2bound`, used by `make_data`, `R/trend.R`,
    `R/omission_decomposition.R`).
  - `do_pre_transform`.
- **`R/design.R`**
  - Where `fill_bound` and `fill_transform` are called (around lines 490 and
    1729).
  - The marginalise check (around line 73), which needs a finite lower bound
    from `bound$minmax`.
- **C++ path** (authoritative; the R likelihood path is outdated and not the
  reference):
  - `src/transform_utils.cpp`: `c_do_bound_pt`, `c_do_bound_pt_from_into`,
    transform specs.
  - `src/utility_functions.h`: `c_do_bound`.
  - `src/particle_ll.cpp`: bound plan, invariant bound specs.
  - `src/marginal_likelihood.cpp`: natural-scale bounds in marginalisation.
  - The transform specs passed from R.
- **Sampler code that reads bounds**
  - `R/shear_move.R`: `.emc_shear_support` / `.emc_shear_walls` turn
    `bound$minmax` plus the transform into per-subject half-spaces on the
    sampled scale.
  - The collapsed move: wall-aware transport.
  - The subject HMC bound-floor handling in `R/subject_hmc.R`.

  Once the transform absorbs a bound, most of those walls disappear. Check
  that this code still behaves correctly, and simplifies, rather than double
  counting.
- **Models**
  - Every model constructor's `transform` and `bound` (`R/model_*.R`).
  - Constructors with arguments (`LBA(posdrift)`, `GOM()`, `ROU()`, `SS*()`).
  - Helpers that add columns (`add_nuisance_pars`, `add_time_warp_par`).
  - An audit of all models was done before (each model's `p_types` default
    must lie strictly inside its own bound); reuse that approach.
- **User-facing scale**
  - `get_pars` (`map = FALSE` default), `summary.emc`, `credint`, and
    plotting that shows parameter values.
  - `prior()` / `R/priors.R` (default priors and how priors are specified).
  - `make_data` parameter vectors.
  - `mapped_pars`.
  - Start values.
- **Documentation**: roxygen for `design`, the model constructors and
  `prior`; NEWS.

## Constraints

- Only the C++ likelihood path is authoritative. Validate against it (and
  quadrature or Monte Carlo where needed), not against `calc_ll_R`.
- Code comments say what and why only: no history, benchmark numbers, fixture
  names or dates.
- Install into a private library (`R CMD INSTALL` with `R_LIBS=`/`EMC_LIB=`),
  check `find.package("EMC2")` reports that library, and
  `git checkout src/Makevars` afterwards.
- Run the full suite inside the package namespace:
  `TESTTHAT_PARALLEL=false NOT_CRAN=true` with
  `test_dir(".", package = "EMC2", load_package = "installed")`.
- Heavy jobs go through `/data/work/EMC2_workstreams/bench_lock.sh`. Stop them
  if the user starts their own fits.
- Do not touch uncommitted user work in the tree (`plans/`, `tools/`, or
  modified files you did not create).

## Deliverables

1. **A short design note**, written first, before any code: the rule for
   linking, how exceptions, trends, one-sided families and user overrides are
   handled, and how users read and write natural-scale values. Agree it with
   the user before implementing.
2. **The implementation, with tests:**
   - Natural-scale mapping is unchanged for old objects.
   - New objects map within the bound and never reach it.
   - Exceptions still work.
   - Overrides are respected and validated.
   - A trend that pushes past a bound is still caught by the strict test.
3. **Verification on a fit** whose posterior sat on a bound (t0 at 0.05): show
   the sampler no longer loses proposals to the cliff, and report how the
   natural-scale posterior changes.
4. **A NEWS entry** describing the model change and the natural-scale
   interface.
