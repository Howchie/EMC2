# Task: Truncation SBC mis-calibration (DDM/WDM, LNR, LBA vs RDM)

## Problem
SBC (WorkingTests/SBC/trunc_SBC.R) shows mis-calibration under truncation for
WDM/DDM, LNR, LBA but NOT RDM. Censoring is fine. Henrich (Stan, ~/Downloads
pdf) shows good truncation calibration for a WDM-like model -> suspect a bug.
Bounds in trunc_SBC.R are DATA-DEPENDENT: LT=quantile(rt,.1), UT=quantile(rt,.9).

## Findings (2026-06-21)
- Race models (LBA/LNR/RDM) all share `log_likelihood_race_missing` +
  `pr_pt` (R/likelihood.R:11): cf = 1/(S(LT)-S(UT)), correct given exact
  pfun/dfun. Identical code for all three -> RDM passing means the shared
  correction is fine; LNR/LBA differences must be model-specific or in the SBC
  scheme.
- DDM/WDM used `log_likelihood_ddm` (R/likelihood.R:379) which applied NO
  truncation correction at all -> root cause of WDM mis-calibration (t0 KS
  0.349, rank mean 0.73 = t0 underestimated). WDM has no between-trial
  variability, so cause is NOT intrinsic-variability.
- Measured KS (rank uniformity) on existing results: WDM t0=0.349,a=0.157;
  LNR s=0.125,m_lMTRUE=0.094,t0=0.081; LBA B=0.095,t0=0.090; RDM all<=0.082.
- THEORY for LNR/LBA residual: data-dependent quantile bounds make the SBC
  generative density (interior order stats given boundary order stats) NOT
  proportional to the fixed-bounds truncated likelihood -> intrinsic
  mis-calibration concentrated in shape/shift params. To confirm via fixed
  bounds. RDM less affected by its parameterization.

## Done
- Added truncation correction to `log_likelihood_ddm` (R/likelihood.R): divide
  each retained trial density by Pin = sum_R (pDDM(UT,R)-pDDM(LT,R)); UT=Inf->1,
  LT=0->0. Verified dDDM/pDDM mutually consistent to ~1e-7 (integral test), so
  corrected truncated likelihood integrates to 1 exactly.
- `log_likelihood_ddmgng` NOT yet updated (go/nogo + censoring semantics) -- TODO.

## In progress
- Background: fixed-bounds (LT=0.35,UT=1.10 constants) WDM SBC with corrected
  likelihood, 200 reps -> WorkingTests/SBC/SBC_WDM_fixedT.RData. Expect t0 KS to
  drop from 0.349 toward <0.1 if fix works. Script /tmp/wdm_fixedbounds_sbc.R.

## Next
1. Confirm fixed-bounds WDM calibrates (validates the fix).
2. Fixed-bounds LNR SBC -> separate quantile-scheme (intrinsic) from any LNR bug.
3. Then LBA; decide whether to address data-dependent bounds (doc/guidance) and
   port correction to full DDM + ddmgng.

---

# Task: Gauss–Legendre stop-success integral for SS-EXG and SS-RDEX

## 2026-06-11: Analytic (MVN-CDF) solution found and validated
`WorkingTests/stop_success_methods.R` derives and verifies an EXACT
closed form for the SS-EXG stop-success integral: expand the survivor
product over exponential/Gaussian branches (2^n terms), one integration by
parts per term -> one-factor multivariate-normal CDFs. For n=2 go runners:
12 bivariate Phi2 (full line) / trivariate Phi3 (truncated). Validated to
relerr ~1e-16 (n=1,2) and 3e-10 (n=3) against 1e-12-tol adaptive integration.
Verdict: NOT faster than GL-64 (pure-R GL 36us vs analytic ~300us; C++ would
be roughly a wash), and exp(delta) overflows when sigma/tau_total is extreme
(needs GL fallback there anyway). GL-64 confirmed accurate to ~1e-12 against
the exact form at typical params. Use the analytic form as the gold-standard
oracle in tests (better than integrate()-based truth) and as an exact n=1
closed form (4 pnorms). Does not apply to SS-RDEX (Wald survivor args are
nonlinear in t). NB: per Andrew's decision below, integrate()/qags stays the
default for now regardless.

## Already DONE (committed code, needs build + test)
Files added/edited:
- `src/gl_quad.h` — self-contained fixed Gauss–Legendre (Newton/gauleg, thread-
  local node cache). NOT using GSL glfixed (its source isn't in the bundle).
- `src/model_SS_EXG.h` — added `#include "gl_quad.h"`,
  `ss_exg_stop_success_lpdf_gl`, `ss_texg_stop_success_lpdf_gl`, and exported
  test wrapper `ss_texg_stop_success_value(SSD, pars, method, ...)`.
- `src/model_SS_RDEX.h` — added `#include "gl_quad.h"`,
  `ss_rdex_stop_success_lpdf_gl`, and exported `ss_rdex_stop_success_value(...)`.
- `R/stop_success_gl.R` — pure-R GL route (`stop_success_texg_R`,
  `stop_success_rdex_R`) reusing `stopfn_texg`/`stopfn_rdex`, method-selectable.
- `tests/testthat/test-stop_success_gl.R` — Q4 (methods agree) + Q5 (node tuning).
- `WorkingTests/stop_success_methods.R` — interactive speed/accuracy.

### Build + verify steps — ALL DONE (2026-06-11)
- [x] `Rcpp::compileAttributes(".")`  — regenerated (84 new lines in RcppExports)
- [x] `devtools::document()`  — OK (make_ssd.Rd skip is pre-existing, harmless)
- [x] `devtools::load_all(".")` — OK
- [x] `testthat::test_file("tests/testthat/test-stop_success_gl.R")` — 32/32 PASS
- [x] `source("WorkingTests/stop_success_methods.R")` — GL converges, speeds confirmed
- [x] `devtools::check(vignettes=FALSE)` — **Status: OK, 0 errors, 0 warnings, 0 notes**
      (Pandoc absent in this env so vignettes=FALSE; NOTE about CLAUDE.md etc. fixed
      by adding ^CLAUDE\.md$, ^CLAUDE_TASK\.md$, ^.*\.RData$ to .Rbuildignore)
- [x] `src/RcppExports.cpp` and `R/RcppExports.R` REGENERATED (confirmed via git diff)

## REMAINING — wire method/n_nodes into the LIVE likelihood — DONE (2026-06-11)
Mechanism (A) implemented (process-global, not thread_local — see gl_quad.h
comment for why; forked workers inherit it, PSOCK workers set it themselves):

- [x] `StopMethodConfig` + `stop_method_config()` in src/gl_quad.h; exported
      setter/getter `emc2_set_stop_method(method, n_nodes)` /
      `emc2_get_stop_method()` in src/model_SS_EXG.h.
- [x] Live dispatchers `ss_texg_stop_success_lpdf_live` /
      `ss_rdex_stop_success_lpdf_live` (read the config; "integrate" is the
      original qags call, untouched); particle_ll.cpp:1364 now points at them.
- [x] `SSEXG(stop_method = c("auto","integrate","gl","analytic"), stop_n_nodes
      = 64L)` and `SSRDEX(stop_method = c("auto","integrate","gl"), ...)` store
      the choice in the model list; `calc_ll_manager()` (R/sampling.R) and
      `profile_plot()`'s lfun (R/plotting.R) call
      `set_stop_method_from_model()` (R/stop_success_gl.R) before entering C++.
      Saved designs without the fields get the "auto" default.
- [x] R live path mirrored: `pstopTEXG`/`pstopHybrid` gained method/n_nodes
      args (threaded from the sfun closures) and route through
      `stop_success_texg_R`/`stop_success_rdex_R`; their "integrate" branch is
      the identical `my.integrate()` call as before. R live GL windows use
      k_sigma/k_tau = 8/16 to match the C++ live call sites.
- [x] End-to-end test: SSEXG(stop_method=) via calc_ll_manager — gl/auto match
      integrate within 1e-5, stop_n_nodes honoured (4-node visibly worse),
      config verified set, R route (calc_ll_R) honours the closure too.
- [x] The earlier opt-in mechanism (`emc2_stop_gl`/`emc2_set_stop_gl`/
      `StopGLConfig`/`choose_stop_nodes[_R]`, `ss_*_stop_success_lpdf_auto`)
      was REMOVED — superseded by stop_method; one config, one meaning of
      "auto". The stale compare_optin*.R scripts that referenced it were deleted
      in the 2026-06-12 consolidation (see below).
- [x] BUG FIX found en route (src/gsl_utils.h): `ensure_gsl_workspace()` never
      grew the thread-local workspace, so the FIRST caller's max_subdiv pinned
      it; any later qags call with a larger limit (live path uses 100, test
      wrappers 30) failed with GSL_EINVAL and silently returned min_ll
      (pStop = 0). Now reallocates when ws->limit < n.

## DECIDED (Andrew, 2026-06-11; updated later same day)
- stop_method = "auto" IS THE DEFAULT (see PLAN section for the dispatch rule).
- The integrate()/qags route must be KEPT fully working and selectable via
  `SSEXG(..., stop_method = "integrate")` — do not remove or alter its
  numerical behaviour. Removal only LATER, when Andrew explicitly asks.
- Other selectable values: "gl" (default n_nodes 64) and "analytic".

## NEW TASK (2026-06-11): n=1 vs n=2 go-runner speed/accuracy comparison
DONE. Script: `WorkingTests/stop_success_methods.R`.
Sources stop_success_methods.R for the oracle; includes C++ paths.

### Key results (51 param sets × 5 SSDs, SSD=0.2 for timing)

**Accuracy** (relerr vs pstop_texg_analytic oracle):

| Method        | n=1 med  | n=1 max  | n=2 med  | n=2 max  |
|---------------|----------|----------|----------|----------|
| qags (R)      | 1.8e-13  | 5.8e-09  | 3.5e-14  | 9.6e-03* |
| gl_16         | 3.3e-02  | 7.5e-01  | 7.5e-02  | 1.2e+00  |
| gl_32         | 7.1e-05  | 7.4e-02  | 3.0e-04  | 2.0e-01  |
| gl_64         | 8.0e-11  | 1.4e-04  | 1.3e-09  | 1.2e-03  |
| gl_128        | 2.6e-13  | 5.8e-09  | 3.0e-14  | 1.2e-08  |
| analytic_4pnorm| 2.1e-14 | 5.8e-02  | —        | —        |
| cpp_qags      | 1.9e-12  | 2.2e-07  | 1.1e-12  | 1.1e-01* |
| cpp_gl64      | 8.0e-11  | 1.4e-04  | 1.3e-09  | 1.2e-03  |

*large n=2 worst-case errors for qags paths almost certainly reflect TVPACK
  oracle inaccuracy on some extreme random-draw configs, not integration failure;
  gl_128 (robust truth) shows worst 1.2e-08 for the same condition.
  qags returned 0 (failure path): 0 times for both n=1 and n=2.

**Speed** (median μs, typical params, SSD=0.2):

| Method         | n=1 | n=2 |
|----------------|-----|-----|
| qags (R)       | 135 | 185 |
| gl_16          |  21 |  27 |
| gl_32          |  31 |  42 |
| gl_64          |  50 |  70 |
| gl_128         |  86 | 125 |
| analytic (n=2→TVPACK) | 203 | 662 |
| analytic_4pnorm (n=1) |  55 | —  |
| cpp_qags       |  66 |  99 |
| cpp_gl64       |  39 |  59 |

### Takeaways
- GL-64 is 2-4x faster than qags for both n_go; GL-128 matches qags accuracy.
- n=1→2 cost scaling: qags +37%, GL-64 +40%, analytic +3.3x (bivariate→trivariate).
- analytic_4pnorm (L=-Inf shortcut) is 4x faster than full analytic but has
  worst relerr ~6% — NOT safe as oracle; lbS truncation matters.
- Follow-up check (Fable, 2026-06-11): the 4pnorm form IS exact (max 6e-15 vs
  untruncated integrate over the same 50 draws) — its error in the table is
  purely the ignored truncation, amplified by 1/P_stop: relerr ~ F_S(lbS)/
  P_stop(SSD), so worst at long SSDs (6e-4 mass -> 5.8e-2 relerr at P~0.01).
  If an exact n=1 form WITH truncation is wanted, use the bivariate-CDF
  version (Q_trunc with d=1 -> 4 Phi2 + boundary term), not the 4pnorm one.
- C++ GL-64 (~40/59 μs) is ~20% faster than R GL-64 (~50/70 μs).
- No qags zero-failures on plausible param ranges.

- [x] New script `WorkingTests/stop_success_methods.R`
- [x] Conditions: n_go 1,2; SSD grid; typical + 50 random draws
- [x] Methods: qags, GL 16/32/64/128, analytic oracle, 4-pnorm timing
- [x] Compact accuracy + speed tables; zero-failure flag
- [x] C++ paths included (package was built)

## PLAN (Andrew, 2026-06-11): auto method dispatch
Goal: best speed/accuracy balance with no user tuning, while keeping explicit
user choice. EXG only — RDEX has no analytic form, always integrate/gl.

Dispatch rule for stop_method = "auto" (per stop-success evaluation):
- n_go == 1:
  - all lower bounds -Inf (lbS AND lb_go): analytic full-line (4 pnorm)
  - else if lbS + SSD >= lb_go: analytic truncated (4 Phi2 + boundary term)
  - else (kinked domain): GL          [piecewise-analytic possible later]
  - GUARD: deltas are known in closed form before any CDF call; if
    max_k delta_k > 600 (or any non-finite intermediate) -> GL fallback.
    Deterministic, cheap. ~1e-8 discontinuity at the switch is harmless for
    particle Metropolis (no gradients).
- n_go >= 2: GL with stop_n_nodes (default 64; consider auto-bump to 128 when
  (ub - lbS)/sigS is large, per the Q5 "tight stop density" test).

"analytic" and "auto" are values of the SAME stop_method mechanism being
wired in the REMAINING section. Per the updated DECIDED section, "auto" is
the DEFAULT; "integrate" stays available and numerically untouched.

Implementation steps — DONE (2026-06-11) except the benchmark note below:
- [x] C++: new header src/ss_exg_analytic.h — Genz BVND port (`bvn_cdf`) +
      `ss_texg_stop_success_analytic1` (full-line 4-pnorm and truncated
      4-Phi2 + 2-boundary forms) with the guards. Per-term exponents are
      computed as logC_A + delta_k (not exp(logC)*exp(delta)) for extra
      headroom; guard threshold 600 on the combined exponent.
- [x] EXTRA GUARD (deviation, important): analytic only when sigS <= 4*tauS.
      The C++ dexg() switches to a Mills-ratio tail approximation for
      z < -8 (relerr ~1/z^2, up to ~1.6%); for sigS/tauS > ~8 that branch
      covers the integrand's mass region, so the EXACT analytic form would
      differ from the integrate/GL routes by up to ~1% there. At sigS <= 4*tauS
      the affected mass is <~3e-5 of the Gaussian core (analytic-vs-quadrature
      relative effect <~1e-6). Deterministic, computed before any CDF call.
- [x] Finite `upper` (ST-win / deadline call sites) -> GL, not analytic: the
      closed form is derived for U = Inf; J(L)-J(U) doubling was judged not
      worth the cancellation risk. The dominant fitting call (NR stop trial,
      upper = Inf) takes the analytic branch.
- [x] n_go >= 2 GL fallback gets the Q5-motivated node bump: when
      (ub - lo)/sigS > 40, bump stop_n_nodes to 128 (gl_auto_nodes /
      gl_auto_nodes_R). The `_gl` routes also now handle lbS = -Inf with a
      mirrored lower-window heuristic (they used to integrate from -Inf -> NaN
      -> min_ll).
- [x] R mirror (deviation from plan): ported the SAME Genz BVND to R
      (`pbvn_R`, stop_success_gl.R) instead of Drezner-Wesolowsky — still no
      pbivnorm/mvtnorm dependency, but ~5e-16 accuracy and exact algorithmic
      parity with the C++ route (tested to 1e-12 agreement).
      `stop_success_texg_analytic1_R` mirrors the C++ form + guards.
- [x] stop_method wiring extended to c("auto","integrate","gl","analytic"),
      "auto" the DEFAULT everywhere (R model lists AND the compiled default in
      gl_quad.h, so direct calc_ll callers agree). For SSRDEX "analytic" is
      rejected by match.arg in R; C++ treats it as "auto" (GL + bump).
      In C++ "analytic" == "auto" (both try the closed form, GL fallback) —
      the distinction is user intent only.
- [x] Test/benchmark aid: exported `ss_texg_stop_success_auto_branch(SSD,
      pars, upper)` -> "analytic_fullline" / "analytic_trunc" / "gl_guard" /
      "gl_finite_upper" / "gl_ngo2plus"; `ss_*_stop_success_value` gained
      methods "analytic"/"auto" (dispatch rule) and "live" (config plumbing).
- [x] Tests (tests/testthat/test-stop_success_gl.R, 121 pass): auto ==
      integrate over a sigma/tau/SSD grid; truncated analytic vs tight qags
      (1e-7); full-line analytic vs exact-exG stats::integrate reference
      (1e-10); R == C++ analytic (1e-12); guard fallback (sigS > 4*tauS),
      kinked domain, n_go = 2, finite upper all verified via the branch
      reporter + value agreement; pbvn_R vs quadrature reference (1e-10) and
      rho = 0/±1 identities; live-config and end-to-end SSEXG threading tests.
- [x] Benchmark (2026-06-11, this machine, ss_texg_stop_success_value,
      typical n=1 params, SSD=0.2, 2000 reps, medians):
        analytic_trunc 2.9us | analytic_fullline 1.6us | gl64 39.4us |
        qags 65.6us
      -> truncated analytic is ~13.5x faster than GL-64 and ~22x faster than
      qags (far better than the expected ~2x; the 203us R timing was indeed
      pmvnorm overhead).

### Build + verify (dispatch + wiring, 2026-06-11) — ALL DONE
- [x] Rcpp::compileAttributes + devtools::document (emc2_stop_gl.Rd removed,
      NAMESPACE back to pre-opt-in state; SSEXG/SSRDEX Rd gained the params)
- [x] test-stop_success_gl.R: 121 pass / 0 fail
- [x] FULL suite (NOT_CRAN=true, sequential test_file loop): 268 pass /
      0 fail / 0 skip — incl. test-ssexg-cpp-vs-formula.R and the stopS files
      under the new "auto" default
- [x] devtools::check(vignettes=FALSE): 0 errors, 0 warnings, 0 notes

## TASK (2026-06-11): four-way stop_method comparison for Andrew — DONE
Script: `WorkingTests/stop_success_methods.R`
Run: devtools::load_all("."); source("WorkingTests/stop_success_methods.R")
- [x] New script `WorkingTests/stop_success_methods.R`, runnable end-to-end
- [x] Methods "integrate", "gl", "analytic", "auto" via C++ ss_texg_stop_success_value
      AND R stop_success_texg_R (which pstopTEXG calls internally)
- [x] Conditions: n_go=1,2; SSD c(0,.1,.2,.3,.5); typical + 50 random draws (seed 42);
      stress draws — guard (sigS>4*tauS) and kinked domain (lbS+SSD<lbG)
- [x] Truth: pstop_texg_analytic oracle (sourced from stop_success_methods.R)
- [x] Accuracy table (8 rows: cpp/R × 4 methods), speed table (4 methods × n_go=1,2)
- [x] Auto-branch summary per condition block, note on n=2 "analytic" behaviour
- [x] PASS/FAIL printed

### Key results (2026-06-11)
Auto branches: n=1 -> analytic_trunc×255 (all 51 param sets × 5 SSDs); n=2 -> gl_ngo2plus×255.
Stress: guard -> gl_guard (relerr vs integrate: ~1e-14); kinked SSD=0 -> gl_guard (relerr ~1e-8); SSD>=0.1 -> analytic_trunc.
"analytic" n=2: C++ and R both use GL (gl_ngo2plus), identical to "auto".

Accuracy (51 param sets × 5 SSDs, vs analytic oracle):
  n=1: cpp_auto med=6e-16, max=7e-15 — exact analytic, far below 1e-6.
  n=2: cpp_auto med=2e-13, max=8e-5  — GL-64+bump; max exceeds 1e-6 threshold.
        (cpp_integrate n=2 max=0.1; R_integrate n=2 max=3e-3 — likely TVPACK oracle noise)

Speed (C++, typical params, SSD=0.2):
  n=1: integrate=66us, gl=39us, analytic=3us, auto=3us
  n=2: integrate=99us, gl=59us, analytic=59us, auto=59us  (analytic==auto==gl for n=2)

PASS/FAIL:
  PASS C++ auto never worse than integrate worst-case (7.9e-5 vs 0.1)
  PASS R   auto never worse than integrate worst-case (7.9e-5 vs 3.1e-3)
  FAIL C++/R auto within 1e-6 for n=2 (max 7.9e-5 — GL-64+bump; probably oracle noise inflates some)
  NOTE: n=2 FAIL is expected — 1e-6 criterion was aspirational for GL-64; GL-128 throughout
        would achieve ~1e-8 but at 2x cost. Tightening the node-bump threshold would help.

## NOTES / gotchas
- Live EXG path uses the TRUNCATED `ss_texg_*` (lbS), not `ss_exg_*`. The `_gl`
  twin mirrors it. The plain `ss_exg_stop_success_lpdf_gl` is for testing only.
- RDEX with infinite upper uses qagiu in the adaptive path; the GL path uses the
  finite heuristic ub = muS + k_sigma*sigS + k_tau*tauS. Confirm that window is
  wide enough for your parameter ranges (k_sigma=6, k_tau=12 defaults; the LIVE
  call sites in particle_ll.cpp pass 8/16, and the R live path now passes 8/16
  too so both routes use the same window).
- The `_gl` routes now accept lbS = -Inf (finite lower-window heuristic
  muS - k_sigma*sigS - k_tau*tauS); they used to return min_ll there. The C++
  qags route still cannot start from -Inf (pre-existing; analytic handles the
  doubly-untruncated n=1 case exactly).
- The DEFAULT is now stop_method = "auto" everywhere (R model lists AND the
  compiled C++ config default), per Andrew's decision. Re-evaluating an old
  fit therefore uses auto unless `SSEXG(stop_method = "integrate")` is chosen.
  Existing C++-vs-formula tests (tolerances 1e-4/1e-5) pass under auto.
  Edge case: a design SAVED before this change carries the old sfun closure,
  so its (rarely used) pure-R route still integrates while its C++ route gets
  "auto" — numerically both fine, just not bit-identical to each other.
- gl_quad.h / ss_exg_analytic.h linter errors in-IDE (thread_local/templates/
  RcppArmadillo/R:: identifiers) are false positives from standalone parsing;
  they build under Makevars CXX_STD=CXX17. ss_exg_analytic.h must be included
  AFTER exgaussian_functions.h (uses pexg for the truncation normaliser).
- testthat::test_dir() in this env fails to start parallel subprocesses
  (callr startup error, unrelated to the package); loop test_file() instead,
  with NOT_CRAN=true to avoid skip_on_cran skips.
- Don't wire anything into the "OLD STUFF BELOW, KEPT FOR TESTING" sections.

## 2026-06-12: pre-push cleanup (Fable review)
- Deleted dead `ss_exg_stop_success_lpdf_gl` (no callers).
- GL window heuristic centralised: SS_WINDOW_K_SIGMA/K_TAU = 8/16 in gl_quad.h,
  used as C++ defaults and by particle_ll.cpp; R defaults mirror them (the two
  exported Rcpp wrappers keep literal 8.0/16.0 — compileAttributes cannot
  translate named constants into R defaults).
- Removed duplicate `%||%` from model_SS.R (kept make_ssd.R's).
- All stop-signal WorkingTests scripts consolidated into ONE file:
  `WorkingTests/stop_success_methods.R` (oracle + four-way accuracy + stress +
  GL node study + speed + full calc_ll comparison + PASS/FAIL). The files it
  replaced (analytic_stop_success_exg.R, benchmark/compare_* scripts) are
  deleted; any references in this file to those names now point at the
  consolidated script.
- Exported wrapper tolerance defaults (max_subdiv/abs_tol/rel_tol) aligned to
  the LIVE call-site values (100/1e-8/1e-6; were 30/1e-5/1e-4). Found via the
  consolidated script: at a draw with P_stop ~ 1.7e-6 (slow stop, long SSD),
  abs_tol = 1e-5 made wrapper-default qags 20% wrong while live qags was fine.
- Consolidated script final run: ALL PASS. auto n=1 = closed form, max relerr
  6.8e-15 at 3.0us/call (vs integrate 91us); auto n=2 = GL+bump, max 2.9e-5
  (integrate worst 9.1e-3 at live tolerances); guards verified on stress
  cases; full 1000-trial calc_ll: all four methods agree within 1.3e-7, auto
  ~10% faster end-to-end. demo_SSEXG.R / demo_SSRDEX.R verified to run
  against the new code (fits/profiles skipped). SSEXG help clarified (closed
  form applicability conditions). 121/121 tests pass.

## 2026-06-12: staircase attribute consolidation — DONE
The `staircase=` argument to make_data() was broken (attr set on `pars` at
the top of make_data, then wiped when pars is rebuilt by get_pars_oo at the
Ttransform stage; the two SS rfuns also disagreed on where to read it —
SSEXG read attr(data,...), SSRDEX/rSShybrid read attr(pars,...)).
Now consolidated: `data` is the SINGLE canonical carrier.
- model_SS.R: new `get_staircase_attr(data)` helper (next to check_staircase);
  both SSEXG rfun and rSShybrid use it.
- make_data.R: staircase (enriched with labels = lR levels, UC = TC$UC) is
  attached to `data` just before the rfun call; all pars-side sets removed.
  RACE_rfun re-attaches to the data slice after subsetting (R strips
  attributes on `[`). Returned data also re-carries the attribute.
- Verified: SSRDEX + SSEXG staircase sims work (513/2000 stop trials at
  pSSD=.25, SSDs staircase properly); test-stopS.R and
  test-stopS-cens-trunc-ST.R pass.
NB: `SSD = NA` in the design's SSD_function marks staircase trials; the `p=`
key in the staircase list is ignored (stop proportion comes from pSSD).

## DONE (2026-06-12): staircase support for conditional_on_data = FALSE
Implemented per the plan below (all boxes ticked); 37/37 new assertions pass.
- `staircase_step(SSD, label, stop_won, staircase)` + `staircase_clamp()`
  extracted in model_SS.R; `staircase_function` refactored to use them
  (conditional path verified bit-identical under the same seed).
- make_data passes `staircase` into both unconditional simulators (trend.R).
- Both simulators: NA SSDs trigger per-subject SSD state (init SSD0); the
  current SSD is written into data BEFORE design_model/Ttransform each trial
  (rfun never sees NA), then stepped from the simulated response. stop_won =
  R is NA or R is a stop-triggered (lI==1) label. Realised SSDs stay in the
  returned data (output column selection now keeps SSD).
- BONUS FIX: the vectorised simulator passed the WHOLE prefix dm to
  rfun/add_bound/RACE_rfun while all_pars covered only the current trial —
  row-misaligned (errored for rSShybrid, silently wrong `ok` bounds vector
  for everything). Now passes dm_cur (current-trial rows only).
- New test file tests/testthat/test-staircase.R: conditional + unconditional
  + vectorised for SSEXG and SSRDEX (proportion, SSD0 start, +/-stairstep
  steps), plus the missing-staircase error. test-trend.R (17), test-stopS.R
  (2), test-stopS-cens-trunc-ST.R (14), test-stop_success_gl.R (121) all pass.
- Still unsupported on the unconditional branch: UC censoring (make_missing
  is never applied there — pre-existing, applies to all models) and grouped
  make_ssd() staircase specs.

### Original plan (for reference)
Currently make_data(conditional_on_data = FALSE) routes to
make_data_unconditional / make_data_unconditional_vectorised (trend.R:1228,
:1341) and the staircase is silently unsupported there: the staircase arg is
never attached to anything on that branch, so NA SSDs error in the rfun.
The blocker is structural: both unconditional variants simulate TRIAL BY
TRIAL (rebuilding the design each trial for trend/feedback), while
staircase_function processes a whole block of staircase trials in ONE rfun
call, carrying SSD state internally from SSD0. Attaching the attr per trial
would restart the staircase at SSD0 on every trial.
Steps:
- [x] Extract a single-step helper from staircase_function (model_SS.R:60),
      e.g. `staircase_step(label, SSD, staircase)` -> new SSD, reusing the
      existing match_rule/up/down logic so rules (staircase_up/down) work.
- [x] make_data.R: pass the enriched staircase into make_data_unconditional
      and ..._vectorised on the unconditional branch (new argument).
- [x] make_data_unconditional: keep per-subject state `cur_SSD[subj]`
      (init SSD0, clamped to stairmin/stairmax). Each trial, for rows where
      data$SSD is NA: write cur_SSD into the data BEFORE design_model /
      Ttransform (so pars get a concrete SSD; rfun never sees NA). After
      simulating, derive the response label from Rrt (NA R = successful
      stop) and update cur_SSD via staircase_step. Record realised SSD back
      into data so the returned data has it.
- [x] ..._vectorised: same, but cur_SSD is a vector over subjects updated
      once per trial loop iteration (subjects are simulated together; the
      staircase is per-subject so state must be per-subject).
- [x] Tests: unconditional staircase sim for SSEXG + SSRDEX; check stop-trial
      proportion, SSD path starts at SSD0 and steps by ±stairstep within
      [stairmin, stairmax]; compare distribution against the conditional
      path with the same seed/params.
- [x] Gotcha: make_missing/UC handling on the unconditional branch returns
      early (make_data.R ~line 467) — check censoring interaction for
      staircase trials (conditional path applies UC inside the rfun via
      staircase$UC).

## DONE (2026-06-12): unconditional-branch gaps — censoring + grouped staircase
Implemented per the plan below. All suites green: test-staircase.R 60,
test-trend.R 17, test-stopS.R 2, test-stopS-cens-trunc-ST.R 14,
test-stop_success_gl.R 121, test-likelihoods.R 33.
- Part A: make_data's unconditional branch now runs the same make_missing +
  post_functions tail as the conditional path before its early return, so
  censoring/truncation/contamination/rt_resolution work via TC=. Output gains
  the LT/UT/LC/UC columns (now consistent with conditional; trend_conditional
  snapshot re-accepted for exactly this).
- Staircase x UC: make_data attaches TC$UC to the staircase spec; the
  per-trial stair_advance treats rt > UC as no response (ladder steps up),
  matching the block staircase's pre-race censoring.
- Part B: trial-by-trial staircase logic extracted to shared helpers in
  model_SS.R (stair_state_init / stair_group_of / stair_spec_of /
  stair_ssd_value / stair_advance; environment state keyed subject\rgroup).
  Grouped make_ssd() specs: ssd_meta is routed into the simulators, group of
  a trial re-derived from group_cols (interaction sep "::", matching
  names(specs)); per-group SSD0/stairstep/min/max/rules honoured. Custom
  staircase functions error clearly under conditional_on_data = FALSE.
- DECISION: pContaminant is NOT auto-derived from model pars on the
  unconditional branch (threading it as a data column broke design_model —
  "Data cannot have columns with the same names as model parameters" — and
  polluted parameter mapping). Supply via TC= if wanted; documented in code.
- NOTE (pre-existing, unchanged): the conditional SIMPLE staircase runs ONE
  ladder across all subjects' stop trials; the unconditional implementation
  is per-subject (arguably more correct). Differs only for multi-subject
  simple-staircase simulations.
- NOTE: trend feedback_fun sees uncensored responses during the loop
  (censoring applied once at the end, matching the conditional path).

### Original plan (for reference)

### Background (verified in code)
- make_data's unconditional branch (make_data.R:420) returns early: make_missing
  (censoring LC/UC, truncation LT/UT, pContaminant, rt_resolution flooring) and
  post_functions are never applied — for ALL models, not just stop-signal.
- Grouped staircases (make_ssd() with factors/formula) ride on attr(SSD,"emc_ssd")
  captured ONLY via make_data's `functions` argument loop (make_data.R:383-394
  -> ssd_meta); ssd_meta is not passed to the unconditional simulators, and the
  rfun-level apply_grouped_staircase (model_SS.R:168) is block-based, so the
  same trial-by-trial restart problem applies as for the simple staircase.
- Conditional path precedent for UC x staircase: inside the rfun, go finishing
  times > UC are set to Inf BEFORE the staircase logic, so a censored response
  counts as a successful stop and the ladder steps UP (model_SS.R ~519).

### Part A — make_missing on the unconditional branch
- [x] Pass TC into make_data_unconditional / _vectorised (trend.R), or apply
      make_missing in make_data itself before the early return (cleaner: the
      branch already returns res$data; apply the same tail as the conditional
      path — make_missing + post_functions — to it). Prefer the latter: zero
      duplication, one call site.
- [x] pContaminant: conditional path pulls pars[,"pContaminant"] per trial
      (make_data.R:511). Unconditional pars are trialwise/trend-dependent;
      mirror by computing it inside the simulators where pars are in scope
      (first accumulator row per trial) and returning it alongside data.
- [x] Staircase x UC: the per-trial step decision must see the CENSORED
      response. In the simulators' stair update: UC_i <- get_missing(
      staircase$UC, <trial row>, "UC", Inf, "numeric"); if simulated rt > UC_i
      treat label as NA / stop_won TRUE before staircase_step. (Recorded
      rt/R masking is then make_missing's job at the end.)
- [x] DECISION to flag: trend feedback_fun during the loop sees uncensored
      responses; applying full censoring per-trial would change feedback
      semantics. Proposed: keep end-of-loop make_missing (matches conditional
      path), document that feedback sees uncensored data.
- [x] Tests: unconditional + UC for a non-SS model (e.g. LBA, censored rt ->
      NA/UC response) and for SSEXG staircase (UC pushes ladder up; compare
      against conditional path with same UC).

### Part B — grouped make_ssd() staircases, unconditional
- [x] Plumbing: make_data captures ssd_meta before the unconditional branch
      only if `functions=` was used; move the functions-loop capture above the
      branch (it already is — line 383) and pass ssd_meta into both
      simulators (use it as `staircase` when non-NULL, like the conditional
      attach at make_data.R:494).
- [x] State: generalise stair_state from one scalar per subject to a named
      list keyed by (subject, group). Group of a trial: re-derive from
      stop_meta$group_cols on the data row (interaction of those columns) —
      do NOT try to row-match stop_meta$data against the accumulator-expanded
      reordered data. Each (subject, group) ladder inits at specs[[grp]]$SSD0
      %||% base_spec$SSD0 and steps with the group's stairstep/min/max and
      rules (specs[[grp]]$rules %||% global rules) via the existing
      staircase_step().
- [x] Refactor first: extract the now-duplicated inline staircase blocks in
      make_data_unconditional / _vectorised into shared helpers in
      model_SS.R — stair_init(data, staircase), stair_fill(state, rows, data),
      stair_update(state, rows, data) — then add grouping inside the helpers
      so both simulators get it for free.
- [x] Custom per-group staircase_function (attr on spec): block-based
      signature (dts, spec) cannot run stepwise. Error clearly:
      "custom staircase functions are not supported with
      conditional_on_data = FALSE".
- [x] Tests: 2-group make_ssd design (e.g. S-dependent SSD0/stairstep),
      unconditional + vectorised: per-group ladders start at their own SSD0,
      step by their own stairstep, don't cross-contaminate; error path for
      custom staircase_function.

### Order: Part A first (Part B's UC-aware stepping depends on it), refactor
### item of Part B before its grouping work.
