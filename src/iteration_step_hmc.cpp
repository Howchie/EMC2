// The iteration's subject HMC step (R/iteration_steps.R .emc_step_hmc) for
// rounds that run in this process: the subjects' net moves into the
// allocation window, each subject's selection from its own stream, the
// selected transitions, and their results into the states, the plan and the
// records.
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include "iteration_state.h"

using namespace Rcpp;
using namespace emc_iter;

SEXP emc_alloc_record_subjects_cpp(List A, List D, NumericMatrix alpha, NumericVector old,
                                   SEXP shift);
SEXP emc_hmc_record_tasks_cpp(List st, List tasks);
List emc2_native_hmc_subjects(NumericMatrix alpha, NumericVector ll, NumericMatrix mu,
                              NumericMatrix var, List seeds, NumericVector eps_base,
                              List plans, List data, SEXP model, CharacterVector par_names);

namespace emc_iter {

SEXP global_seed() {
  return R_getVarEx(Rf_install(".Random.seed"), R_GlobalEnv, FALSE, R_NilValue);
}

void set_global_seed(SEXP seed) {
  if (seed != R_NilValue) Rf_defineVar(Rf_install(".Random.seed"), seed, R_GlobalEnv);
}

double stream_unif(SEXP seed, RObject& next) {
  set_global_seed(seed);
  GetRNGstate();
  const double u = unif_rand();
  PutRNGstate();
  next = global_seed();
  return u;
}

SEXP proposal_rows(SEXP proposals, int p) {
  const int n = Rf_ncols(proposals);
  Shield<SEXP> out(Rf_allocMatrix(REALSXP, p, n));
  const double* src = REAL(proposals);
  double* dst = REAL(out);
  for (int s = 0; s < n; ++s)
    for (int k = 0; k < p; ++k) dst[k + static_cast<R_xlen_t>(s) * p] =
      src[k + static_cast<R_xlen_t>(s) * (p + 1)];
  SEXP dn = Rf_getAttrib(proposals, R_DimNamesSymbol);
  if (dn != R_NilValue) {
    Shield<SEXP> dn2(Rf_allocVector(VECSXP, 2));
    SEXP rn = VECTOR_ELT(dn, 0);
    if (rn != R_NilValue) {
      Shield<SEXP> rn2(Rf_allocVector(STRSXP, p));
      for (int k = 0; k < p; ++k) SET_STRING_ELT(rn2, k, STRING_ELT(rn, k));
      SET_VECTOR_ELT(dn2, 0, rn2);
    }
    SET_VECTOR_ELT(dn2, 1, VECTOR_ELT(dn, 1));
    Rf_setAttrib(out, R_DimNamesSymbol, dn2);
  }
  return out;
}

bool option_on(const char* name) {
  SEXP v = Rf_GetOption1(Rf_install(name));
  return v == R_NilValue || flag(v);
}


bool step_hmc(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), pars = var(st, "pars"), proposals = var(st, "proposals");
  SEXP plan = var(st, "hmc_plan"), streams = var(st, "streams"), wpool = var(st, "wpool");
  if (TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) || TYPEOF(pmwgs) != VECSXP)
    return false;
  const int p = Rf_nrows(proposals) - 1, n = Rf_ncols(proposals);
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  SEXP samples = lget(pmwgs, "samples");
  RObject sampler = pmwgs;

  if (alloc_on) {
    SEXP alloc = lget(pmwgs, "alloc"), w = lget(alloc, "window");
    SEXP A = lget(w, "A"), D = lget(w, "D"), nA = lget(w, "nA");
    SEXP s_alpha = lget(samples, "alpha"), theta_mu = lget(samples, "theta_mu");
    SEXP tmu = lget(pars, "tmu"), idx = lget(samples, "idx");
    if (TYPEOF(A) != VECSXP || TYPEOF(D) != VECSXP || TYPEOF(nA) != INTSXP ||
        XLENGTH(nA) != 1 || TYPEOF(s_alpha) != REALSXP || TYPEOF(theta_mu) != REALSXP ||
        TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != p || Rf_length(idx) != 1)
      return false;
    const R_xlen_t ix = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
    const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
    if (ix < 0 || XLENGTH(s_alpha) < (ix + 1) * slice || XLENGTH(theta_mu) < (ix + 1) * p)
      return false;
    Shield<SEXP> alpha_new(proposal_rows(proposals, p));
    NumericVector old(slice), shift(p);
    std::copy(REAL(s_alpha) + ix * slice, REAL(s_alpha) + (ix + 1) * slice, old.begin());
    for (int k = 0; k < p; ++k) shift[k] = REAL(tmu)[k] - REAL(theta_mu)[k + ix * p];
    Shield<SEXP> r(emc_alloc_record_subjects_cpp(A, D, NumericMatrix(static_cast<SEXP>(alpha_new)), old, shift));
    if (r == R_NilValue) return false;
    RObject w2 = w;
    if (flag(lget(r, "finite"))) {
      w2 = lset(w2, "A", lget(r, "A"));
      w2 = lset(w2, "D", lget(r, "D"));
      w2 = lset(w2, "nA", Rf_ScalarInteger(INTEGER(nA)[0] + 1));
    }
    Shield<SEXP> alloc2(lset(alloc, "window", w2));
    sampler = lset(pmwgs, "alloc", alloc2);
  }

  if (plan == R_NilValue) {
    set_var(st, "pmwgs", sampler);
    set_var(st, "hmc_tasks", R_NilValue);
    set_var(st, "hmc_el", Rf_ScalarReal(0.0));
    return true;
  }

  SEXP budget = lget(cx, "pool_budget"), data = lget(cx, "data"),
       par_names = lget(cx, "par_names"), model = lget(pmwgs, "model");
  SEXP subj_mu = lget(pars, "subj_mu"), tvar = lget(pars, "tvar");
  SEXP hmc_state = lget(pmwgs, "subject_hmc");
  if (TYPEOF(plan) != VECSXP || XLENGTH(plan) != n || TYPEOF(streams) != VECSXP ||
      XLENGTH(streams) != n || TYPEOF(data) != VECSXP || XLENGTH(data) != n ||
      !flag(lget(hmc_state, "supported")))
    return false;
  // Whether the selected transitions run here, on the stored-argument
  // likelihood path; otherwise, and for a transition that path does not
  // serve, through .emc_hmc_dispatch() and so the worker pool.
  const bool local =
      !(flag(lget(wpool, "alive")) && num(lget(wpool, "n"), 1.0) > 1.0) &&
      !(num(lget(budget, "subject"), 1.0) > 1.0) &&
      !(num(lget(budget, "likelihood"), 1.0) > 1.0) &&
      option_on("emc2.native_particle_block") &&
      Rf_inherits(VECTOR_ELT(data, 0), "data.frame") &&
      Rf_getAttrib(VECTOR_ELT(data, 0), Rf_install("pooled_values")) == R_NilValue &&
      TYPEOF(par_names) == STRSXP && XLENGTH(par_names) == p &&
      TYPEOF(subj_mu) == REALSXP && Rf_isMatrix(subj_mu) && Rf_nrows(subj_mu) == p &&
      Rf_ncols(subj_mu) == n && TYPEOF(tvar) == REALSXP && Rf_isMatrix(tvar) &&
      Rf_nrows(tvar) == p && Rf_ncols(tvar) == p;

  const Clock started = clock_now();
  Shield<SEXP> seed0(global_seed());
  SEXP pick = var(st, "hmc_pick");
  Shield<SEXP> streams2(Rf_shallow_duplicate(streams));
  std::vector<int> subs;
  for (int s = 0; s < n; ++s) {
    SEXP ps = VECTOR_ELT(plan, s);
    if (ps == R_NilValue) continue;
    const double rate = num(lget(ps, "rate"), 0.0);
    if (!(rate > 0)) continue;
    bool go;
    if (rate >= 1) {
      go = true;
    } else if (TYPEOF(pick) == LGLSXP && XLENGTH(pick) > s && LOGICAL(pick)[s] != NA_LOGICAL) {
      go = LOGICAL(pick)[s] != 0;
    } else {
      RObject next;
      go = stream_unif(VECTOR_ELT(streams2, s), next) < rate;
      SET_VECTOR_ELT(streams2, s, next);
    }
    if (go) subs.push_back(s);
  }
  const int m = static_cast<int>(subs.size());
  Shield<SEXP> tasks(Rf_allocVector(VECSXP, n));
  RObject new_proposals = proposals, new_plan = plan, new_pool = wpool;
  double hmc_el = 0.0;
  if (m > 0) {
    // Each selected subject's task and stream after its transition.
    List res_tasks(m), res_seeds(m);
    bool ran = false;
    if (local) {
      NumericMatrix a(p, m), mu(p, m);
      NumericVector ll(m), eps(m);
      List seeds(m), plans(m), data_m(m);
      const double* pr = REAL(proposals);
      for (int k = 0; k < m; ++k) {
        const int s = subs[k];
        for (int q = 0; q < p; ++q) {
          a(q, k) = pr[q + static_cast<R_xlen_t>(s) * (p + 1)];
          mu(q, k) = REAL(subj_mu)[q + static_cast<R_xlen_t>(s) * p];
        }
        ll[k] = pr[p + static_cast<R_xlen_t>(s) * (p + 1)];
        SEXP ps = VECTOR_ELT(plan, s);
        eps[k] = num(lget(ps, "eps"), NA_REAL);
        seeds[k] = VECTOR_ELT(streams2, s);
        plans[k] = ps;
        data_m[k] = VECTOR_ELT(data, s);
      }
      List res;
      ran = true;
      try {
        res = emc2_native_hmc_subjects(a, ll, mu, NumericMatrix(tvar), seeds, eps, plans,
                                       data_m, model, CharacterVector(par_names));
      } catch (std::exception&) {
        ran = false;
      }
      if (ran) {
        const LogicalVector done = res["done"];
        for (int k = 0; k < m; ++k) if (!done[k]) ran = false;
      }
      if (ran) {
        res_tasks = as<List>(res["tasks"]);
        res_seeds = as<List>(res["seeds"]);
      }
    }
    if (!ran) {
      static Function dispatch = package_function(".emc_hmc_dispatch");
      IntegerVector subs1(m);
      for (int k = 0; k < m; ++k) subs1[k] = subs[k] + 1;
      const List sent = dispatch(subs1, plan, proposals, pars, lget(pmwgs, "type"), streams2,
                                 data, model, wpool, lget(budget, "subject"),
                                 lget(budget, "likelihood"), par_names, true);
      const List results = sent["results"];
      std::vector<int> at(n, -1);
      for (int k = 0; k < m; ++k) at[subs[k]] = k;
      int filled = 0;
      for (R_xlen_t r = 0; r < results.size(); ++r) {
        const List one = results[r];
        const IntegerVector r_subs = as<IntegerVector>(one["subs"]);
        const List reply = one["res"];
        const List r_tasks = reply["tasks"], r_seeds = reply["seeds"];
        if (r_tasks.size() != r_subs.size() || r_seeds.size() != r_subs.size())
          Rcpp::stop("subject HMC: a reply differs from its request");
        for (R_xlen_t k = 0; k < r_subs.size(); ++k) {
          const int s = r_subs[k] - 1;
          if (s < 0 || s >= n || at[s] < 0) Rcpp::stop("subject HMC: an unrequested subject");
          res_tasks[at[s]] = r_tasks[k];
          res_seeds[at[s]] = r_seeds[k];
          ++filled;
        }
      }
      if (filled != m) Rcpp::stop("subject HMC: a selected subject has no reply");
      SEXP pool = sent["pool"];
      if (pool != R_NilValue && lget(pool, "alive") != R_NilValue)
        new_pool = lset(wpool, "alive", lget(pool, "alive"));
    }
    bool any_accepted = false, any_warm = false;
    for (int k = 0; k < m; ++k) {
      if (flag(lget(res_tasks[k], "accepted"))) any_accepted = true;
      if (flag(lget(VECTOR_ELT(plan, subs[k]), "warmup"))) any_warm = true;
    }
    if (any_accepted) new_proposals = Rf_duplicate(proposals);
    if (any_warm) new_plan = Rf_shallow_duplicate(plan);
    for (int k = 0; k < m; ++k) {
      const int s = subs[k];
      RObject tk = res_tasks[k];
      SET_VECTOR_ELT(streams2, s, res_seeds[k]);
      if (flag(lget(tk, "accepted"))) {
        SEXP ta = lget(tk, "alpha");
        double* np = REAL(new_proposals);
        for (int q = 0; q < p; ++q)
          np[q + static_cast<R_xlen_t>(s) * (p + 1)] = REAL(ta)[q];
        np[p + static_cast<R_xlen_t>(s) * (p + 1)] = Rf_asReal(lget(tk, "ll"));
      }
      if (flag(lget(VECTOR_ELT(plan, s), "warmup"))) {
        static Function warmup_step = package_function(".emc_hmc_warmup_step");
        Shield<SEXP> ps(warmup_step(VECTOR_ELT(new_plan, s), tk));
        SET_VECTOR_ELT(new_plan, s, ps);
        tk = lset(tk, "eps_next", lget(ps, "eps"));
      }
      SET_VECTOR_ELT(tasks, s, tk);
    }
    const Clock ended = clock_now();
    const double master = ended.cpu - started.cpu;
    hmc_el = ended.elapsed - started.elapsed;
    std::vector<double> cpu(m);
    long double total = 0.0L;
    for (int k = 0; k < m; ++k) {
      SEXP tk = VECTOR_ELT(tasks, subs[k]);
      const double work = num(lget(tk, "cpu"), 0.0);
      cpu[k] = work + master / m;
      total += cpu[k];
      RObject t2 = lset(tk, "work", Rf_ScalarReal(work));
      t2 = lset(t2, "cpu", Rf_ScalarReal(cpu[k]));
      SET_VECTOR_ELT(tasks, subs[k], t2);
    }
    if (alloc_on) {
      const double round_cost = reserved_cost(hmc_el, num(var(st, "cores_it"), 1.0));
      const double tot = static_cast<double>(total);
      for (int k = 0; k < m; ++k) {
        const double share = tot > 0 ? cpu[k] / tot : 1.0 / m;
        SET_VECTOR_ELT(tasks, subs[k],
                       lset(VECTOR_ELT(tasks, subs[k]), "cost",
                            Rf_ScalarReal(round_cost * share)));
      }
    }
  } else {
    hmc_el = clock_now().elapsed - started.elapsed;
  }
  Shield<SEXP> recorded(emc_hmc_record_tasks_cpp(hmc_state, List(static_cast<SEXP>(tasks))));
  if (recorded == R_NilValue) {
    set_global_seed(seed0);
    return false;
  }
  sampler = lset(sampler, "subject_hmc", recorded);
  set_global_seed(seed0);
  set_var(st, "pmwgs", sampler);
  set_var(st, "proposals", new_proposals);
  set_var(st, "streams", streams2);
  set_var(st, "hmc_plan", new_plan);
  if (static_cast<SEXP>(new_pool) != wpool) set_var(st, "wpool", new_pool);
  set_var(st, "hmc_tasks", tasks);
  set_var(st, "hmc_el", Rf_ScalarReal(hmc_el));
  return true;
}

}  // namespace emc_iter
