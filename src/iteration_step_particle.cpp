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
  if (TYPEOF(consts) != VECSXP || TYPEOF(pmwgs) != VECSXP || TYPEOF(wpool) != VECSXP ||
      flag(lget(wpool, "alive")) || num(lget(wpool, "n"), 0.0) != 1.0 ||
      lget(wpool, "retained") != R_NilValue || TYPEOF(pm_settings) != VECSXP ||
      TYPEOF(streams) != VECSXP)
    return false;
  const int n = static_cast<int>(XLENGTH(pm_settings));
  SEXP samples = lget(pmwgs, "samples");
  SEXP s_alpha = lget(samples, "alpha"), s_ll = lget(samples, "subj_ll");
  SEXP alpha0 = lget(pars, "alpha"), subj_mu = lget(pars, "subj_mu");
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

}  // namespace

namespace emc_iter {

// In this process where the subjects' steps run here; else the draws here
// and the steps through .emc_step_particle_run(), which holds the worker
// pool's transport.
bool step_particle(SEXP st, SEXP cx, int i) {
  if (particle_local(st, cx, i)) return true;
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
