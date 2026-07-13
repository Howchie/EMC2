# Raw-kernel column-order refactor plan

**Goal:** remove every hard-coded parameter-column index from the C++ likelihood hot
paths, align all batch kernels on a single pointer-array ABI, and eliminate the last
per-particle parameter-staging copies. After this refactor there is exactly one place
in the codebase that says "for LBA, `v` is column 0" — a registry header — and every
kernel family (RACE, DDM, SS) consumes parameters the same way.

## 1. Motivation / current state

Three different parameter-passing conventions coexist today:

| Family | ABI | Column indirection | Per-particle copy |
|---|---|---|---|
| RACE (LBA, RDM, LNR, REXG, RGAMMA, BAwL, RDMSWTN, GBM) | `const double* pars_cm` (contiguous col-major) | **none — offsets hard-coded in each kernel** | yes: `race_staging_buf` memcpy (fast path) or `materialize_reusable()` (LR/mixed/pw paths) |
| DDM | `const double* pars_cm` + `std::vector<int> p_idx` permutation | per-call permutation vector | no (reads ParamTable base directly) |
| SS (TEXG/RDEX) | `const double* const* cols` pointer array | fill-fns + `idx_tf`/`idx_gf` (hard-coded numbers, but localized) | no (pointers into ParamTable base) |

The SS path (`ss_raw.h`, `c_log_likelihood_ss_pt`) is the model to converge on:
`cols[j]` points at ParamTable base column `j`-of-p_types, resolved **by name** once per
`calc_ll_oo` call via `name_to_base_idx`; kernels index columns symbolically. Base
column *positions* are stable across particles (only values are refilled by
`prepare_particle`), so the pointer array never needs rebuilding per particle — the
current per-particle refresh at particle_ll.cpp:2757-2760 is unnecessary and can be
hoisted out of the particle loop.

Costs of the status quo:

1. **Silent coupling:** each kernel's hard-coded offsets implicitly assume the R-side
   `p_types` order for that model (e.g. `LBA: v,sv,B,A,t0`). Reordering `p_types` in a
   model file compiles clean and produces silently wrong likelihoods. Nothing validates
   the contract.
2. **Per-particle copies:** the race fast path memcpys `n_trials × n_par` doubles per
   particle into `race_staging_buf` purely to establish the column order the kernels
   expect (particle_ll.cpp:2926-2936). The LR / mixed / pw paths do the equivalent via
   `materialize_reusable()`. These are the last per-particle bulk copies on hot paths.
3. **Three ABIs** for the same job: harder to review, harder to extend (adding a model
   means picking a convention), and DDM's `p_idx` plumbs a `std::vector<int>&` through
   every call.

## 2. Target design

### 2.1 Single-source column registry — new header `src/col_registry.h`

Per-model namespaces with enums whose order **is** the canonical kernel-facing column
order (identical to today's R `p_types` order, so no numerical behavior changes):

```cpp
namespace emc2col {
  namespace lba     { enum { v, sv, B, A, t0, N_REQ }; }                      // threshold = B + A
  namespace rdm     { enum { v, B, A, t0, s, N_REQ }; }
  namespace lnr     { enum { m, s, t0, N_REQ }; }
  namespace rexg    { enum { mu, sigma, tau, N_REQ }; }
  namespace rgamma  { enum { shape, scale, N_REQ }; }                          // verify names in kernel
  namespace bawl    { enum { v, sv, B, A, t0, /*leak etc.*/ ..., N_REQ }; }    // read off dbawl_raw
  namespace rdmswtn { enum { v, B, A, t0, s_out, sv, mG, mK, omega, N_REQ_BASE }; }
  namespace rdmgbm  { enum { v, B, A, t0, s, mG, mK, omega, N_REQ_BASE }; }
  namespace ddm     { enum { v, a, sv, t0, st0, s, Z, SZ, N_REQ }; }
  namespace ss_texg { enum { mu, sigma, tau, muS, sigmaS, tauS, tf, gf, lb, lbS, N_REQ }; }
  namespace ss_rdex { enum { v, B, A, t0, s, muS, sigmaS, tauS, tf, gf, lbS, N_REQ }; }
}
```

(Exact member lists to be read off the current kernel bodies during implementation —
the enum must reproduce today's offsets verbatim. Optional trailing columns, e.g.
`omega` only for `_EMIX`, `pContaminant` appended last, stay *optional*: kernels only
dereference them when the ctx flags say they exist, same as today.)

Plus a loud-failure validator, run once per `calc_ll_oo`/`calc_ll_oo_pw` dispatch (not
per particle):

```cpp
// Checks that keep_names[0..n-1] == the registry's expected names for this model.
// Converts a silently-wrong p_types reordering into an immediate Rcpp::stop().
inline void validate_col_prefix(const Rcpp::CharacterVector& keep_names,
                                const char* const* expected, int n_expected,
                                const char* model_label);
```

Each registry namespace carries its `expected_names()` array next to the enum, so the
enum and the name list can't drift apart.

### 2.2 Unified batch ABI

```cpp
// utils.h — replaces the pars_cm versions
typedef void (*RaceRawFun)(const double* rt, const double* const* cols, int n_rows,
                           const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx_);

typedef void (*RaceLogSAtTFun)(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int n_par,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* logS_out);
```

Kernel bodies change mechanically:

```cpp
// before                                   // after
const double* v_ = pars_cm + 0 * n_rows;    const double* v_ = cols[emc2col::rdm::v];
```

`n_rows` stays (loop bound); the stride coupling between "matrix leading dimension"
and "number of rows to process" disappears, which is what lets callers pass **views of
any storage** (ParamTable base, a materialized matrix, or a compact GL gather buffer)
without copying into a contiguous col-major block first.

Scalar helpers (`RacePdf1Fun`/`RaceCdf1Fun`, row-major `const double* par`) are **not**
changed: GSL integrands want one accumulator's params contiguous. The row-major
*order* is the same registry order; gathers (`copy_par_row_colmajor`,
`fill_trial_buffers`, LR `parA/parB/parN` fills) either keep reading a materialized
matrix or are switched to read through `cols` (`cols[c][row]` instead of
`pars_cm[c*n + row]`).

### 2.3 Context indices become registry constants

`ContextForRaceModels::t0_index / mean_g_index / mean_k_index / erlang_omega_index`
keep their meaning (index into the cols array) but are assigned from the registry in
`resolve_race_model_adapter` instead of bare integers, e.g.
`out.ctx.t0_index = emc2col::lba::t0;`.

## 3. Full touch-point inventory

### 3.1 Typedefs and kernels — `src/utils.h`

* Typedefs: `RaceRawFun` (utils.h:32), `RaceLogSAtTFun` (utils.h:37). `RacePdf1Fun`/
  `RaceCdf1Fun` (28-29) unchanged.
* Kernels to convert (all in utils.h; each reads `pars_cm + k*n_rows` today):
  * RDM: `drdm_raw` :293, `prdm_raw` :314, `rdm_logS_at_t` :337
  * GBM: `drdmgbm_raw` :366, `prdmgbm_raw` :413, `rdmgbm_logS_at_t` :465
  * LNR: `dlnr_raw` :532, `plnr_raw` :549, `lnr_logS_at_t` :566
  * RGAMMA: `drgamma_raw` :610, `prgamma_raw` :631, `rgamma_logS_at_t` :653
  * REXG: `drexg_raw` :680, `prexg_raw` :699, `rexg_logS_at_t` :719
  * LBA: uses the BAwL raw adapters (`dbawl_raw`, `pbawl_raw`, `bawl_logS_at_t`)
    with `k=0` and clock means fixed off; no optional columns are read
  * BAwL: `dbawl_raw` :860, `pbawl_raw` :898, `bawl_logS_at_t` :938
  * RDMSWTN: `drdmswtn_raw` :1094, `prdmswtn_raw` :1154, `rdmswtn_logS_at_t` :1219
    (these also consult `ctx->mode_hint` at :1225/:1264/:1299/:1315 and read
    mG/mK/omega via ctx indices — those reads become `cols[ctx->mean_g_index]` etc.)
* Delete the layout comments ("// RDM: column layout v=0, B=1 ...") — the registry
  header supersedes them.

### 3.2 Adapter — `src/particle_ll.cpp` `resolve_race_model_adapter` (:474-590)

Replace numeric `t0_index/mean_g_index/mean_k_index/erlang_omega_index` assignments
with registry constants for all 8 model branches (RDMSWTN :505-508, GBM :520-523,
BAwL :531-534, LBA :549, RDM :560, REXG :567, LNR :574, RGAMMA :581).

### 3.3 calc_ll_oo race fast path (:2847-3068)

* **Delete** `race_staging_buf` + the per-particle memcpy (:2851, :2906, :2926-2936).
  Keep `race_base_col_order` (name→base resolution, :2895-2905) and build
  `std::vector<const double*> race_cols` **once before the particle loop**:
  `race_cols[j] = &param_table_template.base(0, race_base_col_order[j])`. Missing
  columns (`bidx < 0`) must disable the raw fast path (today they were silently left
  stale in the staging buffer — this tightens correctness).
* `pContaminant` read (:3042-3043) → `race_cols[fast_pc_staging_col][base]`.
* RDMSWTN mode-hint scan (:2941-2973): `sv_col/lambda_g_col/lambda_k_col` →
  `race_cols[emc2col::rdmswtn::sv]` etc.
* Global-kill `mean_k_ptr` (:3002) → `race_cols[adapter.ctx.mean_k_index]`.
* Kernel invocations :2979/:2985/:2992/:2996 pass `race_cols.data()`.
* Run `validate_col_prefix` once at dispatch.

### 3.4 Mixed path — `c_log_likelihood_race` (:4833-…)

This path **keeps** its materialized `pars` NumericMatrix because of the RACE NA-fill
(:4896-4898) which *writes* `NA_REAL` into rows of `pars` where `!RACE_mask`. That
write must never land in ParamTable base views (base columns for constants are not
refilled per particle, so the corruption would persist). Change is therefore minimal:

* Build a local `const double* cols_local[64]` from `pars` columns
  (`cols_local[j] = &pars(0, j)`), pass to the raw-kernel calls (:5405/:5410,
  :5458/:5461) and `logS_at_t` calls (:5538/:5549).
* Internal row gathers (`fill_trial_buffers` :4978, `par_buf`/`par_buf_inf` loops
  :5185/:5247/:5261/:5281/:5305, `copy_par_row_colmajor`/`row_equal_colmajor` uses)
  keep reading `pars.begin()` — unchanged semantics, no stride ambiguity since the
  matrix is the storage.
* Optional later optimization (out of scope here): audit whether the NA-fill is still
  load-bearing given the leading-rows RACE invariant, and if not, drop materialization
  on this path too.

### 3.5 LogicalRules path — `c_log_likelihood_logicalrules` (:4050-…) and callers

The LR function uses `pars` only via `pars.nrow()`, `pars.ncol()`, `pars.begin()`
(:4079-4087) and never writes to it → **drop materialization entirely**:

* Signature: `(const Rcpp::NumericMatrix& pars, …)` →
  `(const double* const* cols, int n_par, …)`; `n_trials` already comes from `shared`.
* Callers at :2800-2804 (calc_ll_oo) and the pw twin (:3215 area): stop calling
  `pt.materialize_reusable()`; resolve base pointers by name once (same recipe as SS
  :2737-2750) and pass the pointer array. This deletes the last per-particle bulk copy
  in the LR path.
* Inside the function, replace `pars_cm_ptr[c*n_trials + row]` reads:
  * main sweep kernel calls :4120/:4123 → pass `cols`
  * GL pre-pass column gathers :4208-4218 (`src = pars_cm_ptr + p*n_trials` →
    `src = cols[p]`); the compact role buffers `pars_nA/pars_A/pars_nB/pars_B` stay
    contiguous but are handed to kernels through small stack pointer arrays
    (`const double* colptr_nA[64]; colptr_nA[p] = pars_nA.data() + p*n_unique;`).
    Build the four pointer arrays once per particle right after the gather.
  * t0 hoist reads :4164-4166 → `cols[t0_col][inA]`
  * row copies/equality :4449-4466, :4704-4705 → pass `cols` to
    `copy_par_row_colmajor`/`row_equal_colmajor` (change those two helpers, :3553 and
    :3563, to take `const double* const*`; they have no other callers).
* `lr_detection_no_response_prob` and `logicalrules_detection_trial_ll` receive
  row-major per-accumulator buffers (`parA/parB/parN`) — unchanged, only the gathers
  feeding them change source.

### 3.6 `calc_ll_oo_pw` (pointwise twin)

Mirrors calc_ll_oo: LR call site (~:3215), race path via `c_log_likelihood_race`
(covered by 3.4), fast path if present, DDM at :3105-3116. Same edits as their
calc_ll_oo counterparts.

### 3.7 DDM unification — `src/model_DDM.h` + call sites

* `d_DDM_Wien_raw` (:52) / `p_DDM_Wien_raw` (:133): drop the
  `const std::vector<int>& p_idx` parameter and the `pars_cm + p_idx[k]*n_rows`
  arithmetic; take `const double* const* cols` and read
  `cols[emc2col::ddm::v]` … `cols[emc2col::ddm::SZ]`. `n_par` argument becomes unused
  → drop it.
* `c_log_likelihood_DDM_pt` (particle_ll.cpp:1599): parameter
  `(const double* pars_cm, …, const std::vector<int>& p_idx)` → `const double* const* cols`.
  All 15 internal kernel calls (:1618-1868) shed the `p_idx` argument.
* `init_ddm_shared_state` (:2514): today resolves `p_idx[k]` = base column of
  `ddm_names[k]` = {"v","a","sv","t0","st0","s","Z","SZ"}. Replace `ddm_names` with the
  registry `expected_names`, and store `std::vector<const double*> cols` (pointers into
  base) in the DDM shared state — built once, valid for all particles.
* Call sites :2689 and :3116 pass `shared.cols.data()` instead of
  `param_table_template.base.begin()` + `p_idx`.
* Legacy materialized-matrix call site :1988 (`c_log_likelihood_DDM_pt(pars.begin(),…)`)
  builds a local pointer array from the matrix columns permuted by the old name lookup.

### 3.8 SS unification — `src/ss_raw.h`, `src/model_SS_EXG.h`, particle_ll.cpp

Already on the target ABI; remaining work is symbolic indices + one hoist:

* `ss_raw.h`: `ss_texg_fill_acc_raw`/`ss_texg_fill_stop_raw` (:203) — replace bare
  `cols[0..9]` with `emc2col::ss_texg::*`; same for the RDEX twins (:242) with
  `emc2col::ss_rdex::*`; model definitions :221-224 and :266-269 (`idx_tf=6/gf=7`,
  `idx_tf=8/gf=9`) → registry constants.
* particle_ll.cpp `resolve_ss_adapter` (:2028-2029): `is_exg ? 6 : 8` → registry.
* Materialized-fallback stop lambdas (:205-209): `P(0,3..5,9)` / `P(0,5..7,10)` →
  registry constants (these index the same p_types order).
* `model_SS_EXG.h` `pars(i, 0/1/2/8)` sites (:47, :81, :111, :115, :166-169, :226,
  :258, :286-290, :360-362): registry constants.
* Hoist the per-particle `ss_cols[j] = &base(0, idx)` refresh (:2757-2760) above the
  particle loop — base column addresses are particle-invariant.
* No ABI change needed: SS is the reference implementation.

### 3.9 Entry points / R side

Only `calc_ll_oo` and `calc_ll_oo_pw` are exported (RcppExports.cpp:1638-1639); there
is no legacy `calc_ll` entry that feeds these kernels from R-built matrices. No R code
changes, no RcppExports regeneration (no exported signature changes).

## 4. Invariants to preserve

1. **Numerical identity.** Every phase must produce bit-identical `calc_ll_oo` /
   `calc_ll_oo_pw` outputs (checksum comparison) — this is a pure plumbing refactor.
2. **Mask/isok semantics**: `!mask[i]` → `out[i]` untouched; `!isok[i]` →
   `raw_log_zero`/`0.0`. Unchanged.
3. **RACE NA-fill** stays on the materialized-matrix mixed path (3.4); base views are
   never written through.
4. **Optional trailing columns**: kernels may only dereference `cols[k]` for k beyond
   the required prefix when the ctx flags guarantee existence (`kill_active`,
   `kill_shape==3`, `pc_col>=0`) — same discipline as today, now enforced at pointer
   granularity (a missing column is a `nullptr` slot, so a violation segfaults in
   testing instead of reading a neighboring column's data silently).
5. **Leading-rows RACE invariant** (active accumulators first within a trial block)
   still assumed; document it at the cols-array build site.

## 5. Migration phases (each: build → targeted tests → full suite → commit)

* **Phase A — registry + symbolic indices, no ABI change.** Add `col_registry.h`;
  convert `resolve_race_model_adapter`, `resolve_ss_adapter`, `ss_raw.h`,
  `model_SS_EXG.h`, stop lambdas, RDMSWTN mode-hint cols, DDM `ddm_names`, and add
  `validate_col_prefix` at the race/DDM/SS dispatch points. Hoist the SS pointer
  refresh. Zero behavior change; failures can only be validator false-positives, which
  is exactly the information we want early.
* **Phase B — race ABI flip (single atomic commit).** Change the two typedefs, all 24
  kernels, and the four call-site families (fast path, mixed, LR, pw). The compiler
  enforces completeness — anything missed fails to build. Verify with checksum
  benches (LBA/RDM fast path, RDMSWTN, GNG, LR truncated workload) against the
  pre-refactor commit, then the full race/LR/SS test files.
* **Phase C — DDM ABI flip.** As per 3.7. DDM tests + full NOT_CRAN suite +
  LR/plain benchmarks; push.

## 6. Expected wins

* Deletes the per-particle `n_trials × n_par` staging memcpy (fast path) and the
  per-particle `materialize_reusable` copy (LR path) — the last per-particle bulk
  copies on hot paths. (Measured impact will be modest — kernels dominate — but it is
  strictly-less work per particle and unblocks trend/pretransform paths that scale
  `n_par`.)
* One header defines every model's column order; `validate_col_prefix` turns the
  silent p_types-reorder failure mode into an immediate error naming the model and the
  mismatched column.
* One ABI across RACE/DDM/SS: new models copy one pattern; reviewers check one
  contract.
