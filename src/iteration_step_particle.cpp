// The iteration's particle step (R/iteration_steps.R .emc_step_particle) when
// every subject steps in this process on the native kernel: the allocation's
// skips drawn from each subject's own stream, the subjects' steps and their
// proposal settings, the ensemble candidates' likelihoods, and the steps'
// jumps into the allocation window.
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <RcppArmadillo.h>
#include <cmath>
#include <string>
#include <vector>
#include "iteration_state.h"
#include "pool_exchange.h"

using namespace Rcpp;
using namespace emc_iter;

List emc2_native_particle_subjects(List statics, NumericMatrix alpha, NumericVector prev_ll,
                                   NumericMatrix population_mu, NumericMatrix group_var,
                                   List group, List pm_list, List seeds, LogicalVector skip,
                                   std::string stage, List tune_in);
List emc2_ll_cached_subjects(List rows, List data, SEXP model);
SEXP emc_alloc_record_particle_cpp(List P, IntegerVector nP, NumericVector cpu_part,
                                   NumericMatrix before, NumericMatrix after,
                                   LogicalVector ran, double cost, SEXP times);

namespace {

// The settings of a subject's block after a step, in the order
// update_pm_settings first creates them (.emc_native_particle_batch).
SEXP settings_after(SEXP pm_k, SEXP sel, SEXP out, const char* stage, bool adapt) {
  RObject pm = VECTOR_ELT(pm_k, 0);
  if (flag(lget(sel, "moved"))) {
    Shield<SEXP> entry(List::create(_["n"] = lget(sel, "n"), _["disp"] = lget(sel, "disp"),
                                    _["steps"] = lget(sel, "steps")));
    SEXP cur = lget(pm, "sel");
    Shield<SEXP> base(cur == R_NilValue ? Rf_allocVector(VECSXP, 0) : cur);
    Shield<SEXP> stages(lset(base, stage, entry));
    pm = lset(pm, "sel", stages);
  }
  if (adapt) {
    const bool layout = flag(lget(out, "layout_set"));
    pm = lset(pm, "iter", lget(out, "iter"));
    if (layout) {
      pm = lset(pm, "proposal_counts", lget(out, "proposal_counts"));
      pm = lset(pm, "acc_counts", lget(out, "acc_counts"));
      pm = lset(pm, "jump_counts", lget(out, "jump_counts"));
    }
    if (flag(lget(out, "weight_set"))) pm = lset(pm, "weight_ess", lget(out, "weight_ess"));
    if (layout) {
      pm = lset(pm, "step_layout", lget(out, "step_layout"));
      pm = lset(pm, "epsilon", lget(out, "epsilon"));
      pm = lset(pm, "mix", lget(out, "mix"));
    }
    if (flag(lget(out, "particles_set")))
      pm = lset(pm, "n_particles", lget(out, "n_particles"));
  }
  Shield<SEXP> res(Rf_shallow_duplicate(pm_k));
  SET_VECTOR_ELT(res, 0, pm);
  return res;
}

// Each subject's prior mean, parameters x subjects: a standard hierarchy's
// own, the prior mean of a fit without a group level, or the hierarchy
// type's rule; NULL for states the R step handles.
SEXP prior_means(SEXP pars, SEXP cx, SEXP pmwgs) {
  SEXP type = lget(pmwgs, "type"), alpha = lget(pars, "alpha");
  if (TYPEOF(type) != STRSXP || XLENGTH(type) != 1 || TYPEOF(alpha) != REALSXP ||
      !Rf_isMatrix(alpha))
    return R_NilValue;
  const std::string kind = CHAR(STRING_ELT(type, 0));
  const int p = Rf_nrows(alpha), n = Rf_ncols(alpha);
  SEXP subj_mu = lget(pars, "subj_mu");
  if (kind == "standard" && TYPEOF(subj_mu) == REALSXP && Rf_isMatrix(subj_mu) &&
      Rf_nrows(subj_mu) == p && Rf_ncols(subj_mu) == n)
    return subj_mu;
  SEXP tmu = lget(pars, "tmu");
  if (kind == "single" && lget(cx, "single") != R_NilValue && TYPEOF(tmu) == REALSXP &&
      XLENGTH(tmu) == p) {
    Shield<SEXP> out(Rf_allocMatrix(REALSXP, p, n));
    for (int s = 0; s < n; ++s)
      std::copy(REAL(tmu), REAL(tmu) + p, REAL(out) + static_cast<R_xlen_t>(s) * p);
    return out;
  }
  // The structural model's means given each subject's covariates, or one
  // mean for all; the factor and gamma hierarchies' group mean.
  const bool structural = kind == "SEM";
  if (structural && TYPEOF(subj_mu) == REALSXP && Rf_isMatrix(subj_mu) &&
      Rf_nrows(subj_mu) == p && Rf_ncols(subj_mu) >= n) {
    Shield<SEXP> out(Rf_allocMatrix(REALSXP, p, n));
    std::copy(REAL(subj_mu), REAL(subj_mu) + static_cast<R_xlen_t>(p) * n, REAL(out));
    return out;
  }
  SEXP common = structural ? subj_mu : tmu;
  if ((kind == "factor" || kind == "infnt_factor" || kind == "diagonal-gamma" ||
       (structural && !Rf_isMatrix(subj_mu))) &&
      TYPEOF(common) == REALSXP && XLENGTH(common) == p) {
    Shield<SEXP> out(Rf_allocMatrix(REALSXP, p, n));
    for (int s = 0; s < n; ++s)
      std::copy(REAL(common), REAL(common) + p, REAL(out) + static_cast<R_xlen_t>(s) * p);
    return out;
  }
  // Any other state by the hierarchy type's own rule (get_group_level).
  static Function means = package_function(".emc_steps_prior_means");
  return means(pars, type, n);
}

// The allocation's skips of the subjects' particle steps and the joint
// selection with subject HMC, each from the subject's own stream (advanced in
// streams2). Whether any subject may skip.
bool draw_skips(SEXP pmwgs, SEXP plan, bool alloc_on, int n, SEXP streams2,
                LogicalVector& skip, RObject& hmc_pick) {
  SEXP alloc = lget(pmwgs, "alloc");
  bool skipping = false;
  SEXP q_part = alloc_on ? lget(lget(alloc, "particle"), "rate") : R_NilValue;
  if (TYPEOF(q_part) == REALSXP && XLENGTH(q_part) == n) {
    for (int s = 0; s < n; ++s) if (REAL(q_part)[s] < 1) skipping = true;
  }
  if (!skipping) return false;
  std::vector<double> r_now(n, 0.0);
  if (TYPEOF(plan) == VECSXP)
    for (int s = 0; s < n && s < XLENGTH(plan); ++s)
      if (VECTOR_ELT(plan, s) != R_NilValue)
        r_now[s] = num(lget(VECTOR_ELT(plan, s), "rate"), 0.0);
  SEXP r_alloc = lget(lget(alloc, "last"), "rate");
  const bool have_alloc = TYPEOF(r_alloc) == REALSXP && XLENGTH(r_alloc) == n;
  for (int s = 0; s < n; ++s) {
    const double q = REAL(q_part)[s];
    if (!(q < 1)) continue;
    const double ra = have_alloc ? REAL(r_alloc)[s] : r_now[s];
    if (ra > 0 && !(r_now[s] > 0)) continue;
    RObject next;
    const double u = stream_unif(VECTOR_ELT(streams2, s), next);
    SET_VECTOR_ELT(streams2, s, next);
    skip[s] = u >= q;
    if (r_now[s] > 0 && r_now[s] < 1) {
      if (static_cast<SEXP>(hmc_pick) == R_NilValue) {
        LogicalVector fresh(n, NA_LOGICAL);
        hmc_pick = fresh;
      }
      LOGICAL(hmc_pick)[s] = u >= 1 - r_now[s];
    }
  }
  return true;
}

// The step with the subjects' particle steps in this process. FALSE, with the
// state untouched, when they need the worker pool or the full path.
bool particle_local(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), pars = var(st, "pars"), wpool = var(st, "wpool");
  SEXP consts = lget(cx, "particle"), ctx = lget(cx, "wpool_ctx");
  SEXP group_chol = var(st, "group_chol_it"), pm_settings = var(st, "pm_settings");
  SEXP streams = var(st, "streams"), ens_pools = var(st, "ens_pools");
  if (var(st, "pars_comb") != R_NilValue || var(st, "wpool_ctx") != R_NilValue ||
      var(st, "data") != R_NilValue || lget(pmwgs, "pooled") != R_NilValue)
    return false;
  if (TYPEOF(consts) != VECSXP || TYPEOF(pmwgs) != VECSXP || TYPEOF(wpool) != VECSXP ||
      flag(lget(wpool, "alive")) || num(lget(wpool, "n"), 0.0) != 1.0 ||
      lget(wpool, "retained") != R_NilValue || TYPEOF(pm_settings) != VECSXP ||
      TYPEOF(streams) != VECSXP)
    return false;
  const int n = static_cast<int>(XLENGTH(pm_settings));
  SEXP samples = lget(pmwgs, "samples");
  SEXP s_alpha = lget(samples, "alpha"), s_ll = lget(samples, "subj_ll");
  SEXP alpha0 = lget(pars, "alpha");
  Shield<SEXP> subj_mu(prior_means(pars, cx, pmwgs));
  SEXP gf = VECTOR_ELT(lget(group_chol, "f"), 0), gv = lget(group_chol, "ref");
  if (TYPEOF(alpha0) != REALSXP || !Rf_isMatrix(alpha0) || Rf_ncols(alpha0) != n ||
      XLENGTH(streams) != n || TYPEOF(s_alpha) != REALSXP || TYPEOF(s_ll) != REALSXP)
    return false;
  const int p = Rf_nrows(alpha0);
  SEXP dn = Rf_getAttrib(alpha0, R_DimNamesSymbol);
  if (dn == R_NilValue || VECTOR_ELT(dn, 0) == R_NilValue || TYPEOF(gf) != VECSXP ||
      flag(lget(gf, "repaired")) || TYPEOF(gv) != REALSXP || !Rf_isMatrix(gv) ||
      Rf_nrows(gv) != p || Rf_ncols(gv) != p || TYPEOF(subj_mu) != REALSXP ||
      !Rf_isMatrix(subj_mu) || Rf_nrows(subj_mu) != p || Rf_ncols(subj_mu) != n)
    return false;
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (prev < 0 || XLENGTH(s_alpha) < (prev + 1) * slice || XLENGTH(s_ll) < (prev + 1) * n)
    return false;
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  SEXP alloc = lget(pmwgs, "alloc"), w = lget(alloc, "window");
  if (alloc_on && (TYPEOF(lget(w, "P")) != VECSXP || TYPEOF(lget(w, "nP")) != INTSXP ||
                   TYPEOF(lget(w, "cpu_part")) != REALSXP))
    return false;
  SEXP extra = lget(ens_pools, "rows");
  if (ens_pools != R_NilValue) {
    if (TYPEOF(extra) != VECSXP || XLENGTH(extra) != n) return false;
    for (int s = 0; s < n; ++s) if (VECTOR_ELT(extra, s) == R_NilValue) return false;
  }
  Shield<SEXP> seed0(global_seed());

  Shield<SEXP> streams2(Rf_shallow_duplicate(streams));
  LogicalVector skip(n);
  RObject hmc_pick;
  draw_skips(pmwgs, var(st, "hmc_plan"), alloc_on, n, streams2, skip, hmc_pick);
  const Clock started = clock_now();

  // The subjects' states and prior means as .emc_wpool_iter sends them.
  Shield<SEXP> alpha(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(alpha0), REAL(alpha0) + slice, REAL(alpha));
  Rf_setAttrib(alpha, R_DimNamesSymbol, dn);
  Shield<SEXP> pop_mu(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(subj_mu), REAL(subj_mu) + slice, REAL(pop_mu));
  {
    Shield<SEXP> dn2(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn2, 0, VECTOR_ELT(dn, 0));
    Rf_setAttrib(pop_mu, R_DimNamesSymbol, dn2);
  }
  NumericVector prev_ll(n);
  std::copy(REAL(s_ll) + prev * n, REAL(s_ll) + (prev + 1) * n, prev_ll.begin());
  Shield<SEXP> subs(Rf_allocVector(INTSXP, n));
  for (int s = 0; s < n; ++s) INTEGER(subs)[s] = s + 1;

  const std::string stage = as<std::string>(lget(cx, "stage"));
  const bool adapt = stage != "sample";
  List res, extra_ll;
  bool ok = true;
  try {
    static Function statics_for = package_function(".emc_native_statics");
    Shield<SEXP> par_names(VECTOR_ELT(dn, 0));
    Shield<SEXP> statics(statics_for(subs, ctx, par_names));
    List group = List::create(_["R"] = lget(gf, "R"), _["rooti"] = lget(gf, "rooti"),
                              _["sld"] = lget(gf, "sum_log_diag"));
    res = emc2_native_particle_subjects(List(static_cast<SEXP>(statics)),
                                        NumericMatrix(static_cast<SEXP>(alpha)), prev_ll,
                                        NumericMatrix(static_cast<SEXP>(pop_mu)),
                                        NumericMatrix(gv), group, List(pm_settings),
                                        List(static_cast<SEXP>(streams2)), skip, stage,
                                        List(lget(consts, "tune_in")));
    if (ens_pools != R_NilValue)
      extra_ll = emc2_ll_cached_subjects(List(extra), List(lget(ctx, "data")),
                                         lget(ctx, "model"));
  } catch (std::exception&) {
    ok = false;
  }
  if (ok) {
    const LogicalVector done = res["done"];
    for (int s = 0; s < n; ++s) if (!done[s] && !skip[s]) ok = false;
    if (ens_pools != R_NilValue)
      for (int s = 0; ok && s < n; ++s) if (Rf_isNull(extra_ll[s])) ok = false;
  }
  if (!ok) {
    set_global_seed(seed0);
    return false;
  }

  const LogicalVector done = res["done"];
  const NumericMatrix r_alpha = res["alpha"];
  const NumericVector r_ll = res["ll"], r_cpu = res["cpu"];
  const List r_pm = res["pm"], r_sel = res["sel"], r_seeds = res["seeds"];
  Shield<SEXP> props(Rf_allocMatrix(REALSXP, p + 1, n));
  Shield<SEXP> pm2(Rf_shallow_duplicate(pm_settings));
  NumericVector times(n);
  double* pr = REAL(props);
  for (int s = 0; s < n; ++s) {
    double* col = pr + static_cast<R_xlen_t>(s) * (p + 1);
    if (done[s]) {
      for (int k = 0; k < p; ++k) col[k] = r_alpha(k, s);
      col[p] = r_ll[s];
      SET_VECTOR_ELT(pm2, s, settings_after(VECTOR_ELT(pm_settings, s), r_sel[s], r_pm[s],
                                            stage.c_str(), adapt));
      SET_VECTOR_ELT(streams2, s, r_seeds[s]);
      times[s] = r_cpu[s];
    } else {
      for (int k = 0; k < p; ++k) col[k] = REAL(alpha)[k + static_cast<R_xlen_t>(s) * p];
      col[p] = prev_ll[s];
    }
  }

  // The steps' jumps and costs into the allocation window.
  RObject sampler = pmwgs;
  if (alloc_on) {
    LogicalVector ran(n);
    bool any_ran = false;
    for (int s = 0; s < n; ++s) {
      ran[s] = !skip[s];
      any_ran = any_ran || ran[s];
    }
    if (any_ran) {
      NumericMatrix before(p, n);
      std::copy(REAL(s_alpha) + prev * slice, REAL(s_alpha) + (prev + 1) * slice,
                before.begin());
      Shield<SEXP> after(proposal_rows(props, p));
      const double cost = reserved_cost(clock_now().elapsed - started.elapsed,
                                        num(var(st, "cores_it"), 1.0));
      Shield<SEXP> r(emc_alloc_record_particle_cpp(
          List(lget(w, "P")), IntegerVector(lget(w, "nP")), NumericVector(lget(w, "cpu_part")),
          before, NumericMatrix(static_cast<SEXP>(after)), ran, cost, times));
      if (r == R_NilValue) {
        set_global_seed(seed0);
        return false;
      }
      RObject w2 = lset(w, "P", lget(r, "P"));
      w2 = lset(w2, "nP", lget(r, "nP"));
      w2 = lset(w2, "cpu_part", lget(r, "cpu_part"));
      Shield<SEXP> alloc2(lset(alloc, "window", w2));
      sampler = lset(pmwgs, "alloc", alloc2);
    }
  }

  // The subjects' relative costs for the next partition.
  SEXP cost0 = var(st, "wpool_cost");
  if (TYPEOF(cost0) == REALSXP && XLENGTH(cost0) == n) {
    double low = R_PosInf;
    for (int s = 0; s < n; ++s) if (!skip[s] && times[s] > 0 && times[s] < low) low = times[s];
    if (R_FINITE(low)) {
      Shield<SEXP> cost(Rf_duplicate(cost0));
      for (int s = 0; s < n; ++s)
        if (!skip[s]) REAL(cost)[s] = std::max(static_cast<double>(times[s]), low);
      set_var(st, "wpool_cost", cost);
    }
  }
  set_var(st, "wpool_part", List::create(subs));
  if (ens_pools != R_NilValue) set_var(st, "ens_pools", lset(ens_pools, "ll", extra_ll));
  set_global_seed(seed0);
  set_var(st, "pmwgs", sampler);
  set_var(st, "proposals", props);
  set_var(st, "pm_settings", pm2);
  set_var(st, "streams", streams2);
  set_var(st, "hmc_pick", hmc_pick);
  return true;
}

// x + y of two reject tallies.
SEXP tally_sum(SEXP x, SEXP y) {
  if (TYPEOF(x) == INTSXP && TYPEOF(y) == INTSXP && XLENGTH(x) == XLENGTH(y)) {
    Shield<SEXP> out(Rf_duplicate(x));
    for (R_xlen_t k = 0; k < XLENGTH(x); ++k) INTEGER(out)[k] += INTEGER(y)[k];
    return out;
  }
  static Function plus("+");
  return plus(x, y);
}

// The step with the subjects' particle steps on the worker pool's shared
// subject queue: the pool's upkeep, the round and its records. FALSE, with
// the state untouched, for states the R step handles; a pool that does not
// serve the queue leaves the round to .emc_step_particle_run().
bool particle_pool(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), wpool = var(st, "wpool");
  // With nuisance parameters the subjects see the merged group state; under
  // pooled parameters the pool's context carries the current data.
  SEXP pars = var(st, "pars_comb") != R_NilValue ? var(st, "pars_comb") : var(st, "pars");
  SEXP ctx = var(st, "wpool_ctx") != R_NilValue ? var(st, "wpool_ctx") : lget(cx, "wpool_ctx");
  SEXP group_chol = var(st, "group_chol_it");
  SEXP pm_settings = var(st, "pm_settings"), streams = var(st, "streams");
  SEXP ens_pools = var(st, "ens_pools"), cost0 = var(st, "wpool_cost");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(wpool) != VECSXP || TYPEOF(pm_settings) != VECSXP ||
      TYPEOF(streams) != VECSXP || var(st, "llpool") != R_NilValue ||
      flag(lget(cx, "sampler_profile")))
    return false;
  const int n = static_cast<int>(XLENGTH(pm_settings));
  SEXP samples = lget(pmwgs, "samples");
  SEXP s_alpha = lget(samples, "alpha"), s_ll = lget(samples, "subj_ll");
  SEXP alpha0 = lget(pars, "alpha");
  Shield<SEXP> subj_mu(prior_means(pars, cx, pmwgs));
  if (TYPEOF(alpha0) != REALSXP || !Rf_isMatrix(alpha0) || Rf_ncols(alpha0) != n ||
      XLENGTH(streams) != n || TYPEOF(s_alpha) != REALSXP || TYPEOF(s_ll) != REALSXP ||
      TYPEOF(cost0) != REALSXP || XLENGTH(cost0) != n ||
      Rf_asInteger(lget(pmwgs, "n_subjects")) != n)
    return false;
  const int p = Rf_nrows(alpha0);
  if (TYPEOF(subj_mu) != REALSXP || !Rf_isMatrix(subj_mu) || Rf_nrows(subj_mu) != p ||
      Rf_ncols(subj_mu) != n || Rf_asInteger(lget(ctx, "n_pars")) != p)
    return false;
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (prev < 0 || XLENGTH(s_alpha) < (prev + 1) * slice || XLENGTH(s_ll) < (prev + 1) * n)
    return false;
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  SEXP extra = lget(ens_pools, "rows");
  Shield<SEXP> seed0(global_seed());

  Shield<SEXP> streams2(Rf_shallow_duplicate(streams));
  LogicalVector skip(n);
  RObject hmc_pick;
  const bool skipping = draw_skips(pmwgs, var(st, "hmc_plan"), alloc_on, n, streams2, skip,
                                   hmc_pick);
  set_global_seed(seed0);
  const Clock started = clock_now();

  // The pool's upkeep: recycling on its schedule, growth onto freed cores.
  RObject pool = wpool;
  const bool persistent = flag(lget(cx, "persistent_workers"));
  const int at = persistent ? Rf_asInteger(var(st, "pool_iteration")) : i;
  const int every = Rf_asInteger(lget(cx, "recycle_every"));
  if (every != NA_INTEGER && every > 0 && at > 1 && Rf_asInteger(lget(pool, "n")) > 1 &&
      ((at - 1) % every) == 0) {
    static Function recycle = package_function(".emc_wpool_recycle");
    pool = recycle(pool, at, every, ctx);
  }
  const double chain_cores = emc_pool::cores_now(lget(cx, "core_ctl"),
                                                 Rf_asReal(lget(cx, "n_cores")));
  const int per_worker = std::max(1, Rf_asInteger(lget(lget(cx, "pool_budget"), "likelihood")));
  const int target = std::min(n, std::max(1, static_cast<int>(chain_cores) / per_worker));
  const int n_before = Rf_asInteger(lget(pool, "n"));
  if (target > n_before) {
    static Function grow = package_function(".emc_wpool_grow");
    pool = grow(pool, target, ctx);
  }
  if (persistent) pool = lset(pool, "iterations", var(st, "pool_iteration"));
  const int n_workers = Rf_asInteger(lget(pool, "n"));
  RObject part = var(st, "wpool_part");
  if (n_workers != n_before || skipping) {
    std::vector<double> cost(REAL(cost0), REAL(cost0) + n);
    if (skipping) for (int s = 0; s < n; ++s) if (skip[s]) cost[s] = 0.0;
    part = emc_pool::lpt_partition(cost.data(), n, n_workers);
  }

  // The subjects' states and prior means as .emc_wpool_iter sends them.
  SEXP dn = Rf_getAttrib(alpha0, R_DimNamesSymbol);
  Shield<SEXP> alpha(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(alpha0), REAL(alpha0) + slice, REAL(alpha));
  Shield<SEXP> pop_mu(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(subj_mu), REAL(subj_mu) + slice, REAL(pop_mu));
  if (dn != R_NilValue) {
    Rf_setAttrib(alpha, R_DimNamesSymbol, dn);
    if (VECTOR_ELT(dn, 0) != R_NilValue) {
      Shield<SEXP> dn2(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(dn2, 0, VECTOR_ELT(dn, 0));
      Rf_setAttrib(pop_mu, R_DimNamesSymbol, dn2);
    }
  }
  NumericVector prev_ll(n);
  std::copy(REAL(s_ll) + prev * n, REAL(s_ll) + (prev + 1) * n, prev_ll.begin());
  SEXP ll_dn = Rf_getAttrib(s_ll, R_DimNamesSymbol);
  if (ll_dn != R_NilValue && VECTOR_ELT(ll_dn, 0) != R_NilValue)
    prev_ll.attr("names") = VECTOR_ELT(ll_dn, 0);

  emc_pool::ParticleRequest request;
  request.ctx = ctx;
  request.part = part;
  request.alpha = alpha;
  request.population_mu = pop_mu;
  request.group_chol = group_chol;
  request.pm_settings = pm_settings;
  request.prev_ll = prev_ll;
  request.seeds = streams2;
  request.extra = extra;
  request.pooled = lget(lget(pmwgs, "pooled"), "value");
  request.skip = skipping ? static_cast<SEXP>(skip) : R_NilValue;
  request.n_pars = p;
  emc_pool::ParticleReply reply;
  if (!emc_pool::particle_queue(pool, request, reply)) {
    set_var(st, "wpool", pool);
    set_var(st, "wpool_part", part);
    static Function run = package_function(".emc_step_particle_run");
    run(st, cx, Rf_ScalarInteger(i), request.skip, hmc_pick, streams2, true);
    return true;
  }

  const NumericVector times(static_cast<SEXP>(reply.times));
  if (ens_pools != R_NilValue) set_var(st, "ens_pools", lset(ens_pools, "ll", reply.extra_ll));
  {
    SEXP tallies = var(st, "stage_rejects");
    Shield<SEXP> sum(tally_sum(lget(tallies, "particle"), reply.rejects));
    set_var(st, "stage_rejects", lset(tallies, "particle", sum));
    SEXP counts = lget(reply.ensemble_rejects, "counts");
    bool any = false;
    for (R_xlen_t k = 0; TYPEOF(counts) == INTSXP && k < XLENGTH(counts); ++k)
      if (INTEGER(counts)[k] >= 1) any = true;
    if (any || TYPEOF(counts) != INTSXP) {
      static Function add = package_function(".emc_reject_add");
      add(reply.ensemble_rejects, "ensemble");
    }
  }

  // The subjects' relative costs and the next round's partition.
  Shield<SEXP> cost(Rf_duplicate(cost0));
  {
    double low = R_PosInf;
    for (int s = 0; s < n; ++s) if (!skip[s] && times[s] > 0 && times[s] < low) low = times[s];
    if (R_FINITE(low))
      for (int s = 0; s < n; ++s)
        if (!skip[s]) REAL(cost)[s] = std::max(static_cast<double>(times[s]), low);
  }
  Shield<SEXP> next_part(emc_pool::lpt_partition(REAL(cost), n, Rf_asInteger(lget(pool, "n"))));

  // The steps' jumps and costs into the allocation window. The pass also
  // computes the ensemble candidates: the particle steps' share of it is
  // their share of the workers' CPU.
  RObject sampler = pmwgs;
  if (alloc_on) {
    double share = 1.0;
    if (reply.worker_cpu > 0) {
      long double total = 0.0L;
      for (int s = 0; s < n; ++s) total += times[s];
      share = std::min(1.0, static_cast<double>(total) / reply.worker_cpu);
    }
    LogicalVector ran(n);
    bool any_ran = false;
    for (int s = 0; s < n; ++s) {
      ran[s] = !skip[s];
      any_ran = any_ran || ran[s];
    }
    if (any_ran) {
      NumericMatrix before(p, n);
      std::copy(REAL(s_alpha) + prev * slice, REAL(s_alpha) + (prev + 1) * slice,
                before.begin());
      Shield<SEXP> after(proposal_rows(reply.props, p));
      const double cost_now = share * reserved_cost(clock_now().elapsed - started.elapsed,
                                                    num(var(st, "cores_it"), 1.0));
      SEXP alloc = lget(pmwgs, "alloc"), w = lget(alloc, "window");
      RObject r;
      if (TYPEOF(lget(w, "P")) == VECSXP && TYPEOF(lget(w, "nP")) == INTSXP &&
          TYPEOF(lget(w, "cpu_part")) == REALSXP)
        r = emc_alloc_record_particle_cpp(
            List(lget(w, "P")), IntegerVector(lget(w, "nP")),
            NumericVector(lget(w, "cpu_part")), before, NumericMatrix(static_cast<SEXP>(after)),
            ran, cost_now, times);
      if (!r.isNULL()) {
        RObject w2 = lset(w, "P", lget(r, "P"));
        w2 = lset(w2, "nP", lget(r, "nP"));
        w2 = lset(w2, "cpu_part", lget(r, "cpu_part"));
        Shield<SEXP> alloc2(lset(alloc, "window", w2));
        sampler = lset(pmwgs, "alloc", alloc2);
      } else {
        static Function record = package_function(".emc_alloc_record_particle");
        sampler = record(pmwgs, before, after, ran, cost_now, times);
      }
    }
  }

  set_var(st, "pmwgs", sampler);
  set_var(st, "proposals", reply.props);
  set_var(st, "pm_settings", reply.pm_settings);
  set_var(st, "streams", reply.seeds);
  set_var(st, "wpool", pool);
  set_var(st, "wpool_part", next_part);
  set_var(st, "wpool_cost", cost);
  set_var(st, "hmc_pick", hmc_pick);
  return true;
}

}  // namespace

namespace emc_iter {

// In this process where the subjects' steps run here, or on the worker
// pool's shared subject queue; else the draws here and the steps through
// .emc_step_particle_run().
bool step_particle(SEXP st, SEXP cx, int i) {
  if (particle_local(st, cx, i)) return true;
  if (particle_pool(st, cx, i)) return true;
  SEXP pmwgs = var(st, "pmwgs"), streams = var(st, "streams");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(streams) != VECSXP ||
      Rf_length(lget(pmwgs, "n_subjects")) != 1)
    return false;
  const int n = Rf_asInteger(lget(pmwgs, "n_subjects"));
  if (XLENGTH(streams) != n) return false;
  Shield<SEXP> seed0(global_seed());
  Shield<SEXP> streams2(Rf_shallow_duplicate(streams));
  LogicalVector skip(n);
  RObject hmc_pick;
  const bool skipping = draw_skips(pmwgs, var(st, "hmc_plan"), flag(lget(cx, "alloc_on")), n,
                                   streams2, skip, hmc_pick);
  set_global_seed(seed0);
  static Function run = package_function(".emc_step_particle_run");
  run(st, cx, Rf_ScalarInteger(i), skipping ? static_cast<SEXP>(skip) : R_NilValue, hmc_pick,
      streams2);
  return true;
}

}  // namespace emc_iter
