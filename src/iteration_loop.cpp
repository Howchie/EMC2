// The sampler's iteration loop (R/iteration_steps.R): each iteration is the
// step sequence over the state environment, a step running here where its
// native form covers the state it finds and through its R function where it
// does not. A native step leaves the state it would have left in R, on the
// same random streams; it declines before changing anything.
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <RcppArmadillo.h>
#include <cmath>
#include <cstring>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#include <time.h>
#endif
#include "iteration_state.h"

using namespace Rcpp;

void emc_set_last_slice(NumericVector target, int j, NumericVector values);
void emc_set_particle_slice(NumericVector alpha, NumericVector subj_ll,
                            NumericMatrix proposals, int j, int n_pars);
SEXP emc_alloc_record_tasks_cpp(NumericVector cpu_hmc, IntegerVector hmc_att,
                                NumericVector rows_hmc, List tasks);

namespace emc_iter {

SEXP lget(SEXP x, const char* name) {
  if (TYPEOF(x) != VECSXP) return R_NilValue;
  SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
  if (nm == R_NilValue) return R_NilValue;
  const R_xlen_t n = XLENGTH(x);
  for (R_xlen_t i = 0; i < n; ++i)
    if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) return VECTOR_ELT(x, i);
  return R_NilValue;
}

SEXP lset(SEXP x, const char* name, SEXP v) {
  Shield<SEXP> keep(v);
  if (TYPEOF(x) != VECSXP) Rcpp::stop("iteration state: not a list");
  SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
  const R_xlen_t n = XLENGTH(x);
  if (nm != R_NilValue)
    for (R_xlen_t i = 0; i < n; ++i)
      if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) {
        Shield<SEXP> out(Rf_shallow_duplicate(x));
        SET_VECTOR_ELT(out, i, v);
        return out;
      }
  Shield<SEXP> out(Rf_allocVector(VECSXP, n + 1));
  Shield<SEXP> names(Rf_allocVector(STRSXP, n + 1));
  for (R_xlen_t i = 0; i < n; ++i) {
    SET_VECTOR_ELT(out, i, VECTOR_ELT(x, i));
    SET_STRING_ELT(names, i, nm == R_NilValue ? R_BlankString : STRING_ELT(nm, i));
  }
  SET_VECTOR_ELT(out, n, v);
  SET_STRING_ELT(names, n, Rf_mkChar(name));
  Rf_copyMostAttrib(x, out);
  Rf_setAttrib(out, R_NamesSymbol, names);
  return out;
}

SEXP var(SEXP env, const char* name) {
  return R_getVarEx(Rf_install(name), env, FALSE, R_NilValue);
}

void set_var(SEXP env, const char* name, SEXP v) {
  Rf_defineVar(Rf_install(name), v, env);
}

double num(SEXP x, double if_null) {
  if (x == R_NilValue || Rf_length(x) < 1) return if_null;
  return Rf_asReal(x);
}

bool flag(SEXP x) {
  return TYPEOF(x) == LGLSXP && Rf_length(x) == 1 && LOGICAL(x)[0] == TRUE;
}

namespace {

Clock r_proc_time() {
  static SEXP call = NULL;
  if (call == NULL) {
    call = Rf_lang1(Rf_install("proc.time"));
    R_PreserveObject(call);
  }
  Shield<SEXP> t(Rf_eval(call, R_BaseEnv));
  Clock c;
  c.cpu = REAL(t)[0] + REAL(t)[1];
  c.elapsed = REAL(t)[2];
  return c;
}

}  // namespace

#if defined(__unix__) || defined(__APPLE__)
// The system's clocks read as proc.time() reads them: elapsed time from R's
// own origin (found once against proc.time()), processor time in whole
// milliseconds per component.
Clock clock_now() {
  static bool have_origin = false;
  static double origin = 0.0;
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  const double wall = static_cast<double>(ts.tv_sec) + 1e-9 * static_cast<double>(ts.tv_nsec);
  if (!have_origin) {
    origin = r_proc_time().elapsed - wall;
    have_origin = true;
  }
  struct rusage self;
  getrusage(RUSAGE_SELF, &self);
  Clock c;
  c.elapsed = wall + origin;
  c.cpu = (static_cast<double>(self.ru_utime.tv_sec) + 1e-3 * (self.ru_utime.tv_usec / 1000)) +
          (static_cast<double>(self.ru_stime.tv_sec) + 1e-3 * (self.ru_stime.tv_usec / 1000));
  return c;
}
#else
Clock clock_now() { return r_proc_time(); }
#endif

double reserved_cost(double elapsed, double cores) {
  return std::max(1.0, cores) * std::max(0.0, elapsed);
}

SEXP iteration_index(SEXP cx, int i) {
  SEXP start = lget(cx, "start_iter");
  if (TYPEOF(start) == INTSXP) return Rf_ScalarInteger(INTEGER(start)[0] + i);
  return Rf_ScalarReal(Rf_asReal(start) + i);
}

}  // namespace emc_iter

namespace {

using namespace emc_iter;

bool real_scalar(SEXP x) { return TYPEOF(x) == REALSXP && XLENGTH(x) == 1; }
bool int_scalar(SEXP x) { return TYPEOF(x) == INTSXP && XLENGTH(x) == 1; }

// One move family's cost record with a move's cost added.
SEXP family_add(SEXP f, double cpu, double own, double n_moves) {
  Shield<SEXP> a(lset(f, "cpu", Rf_ScalarReal(REAL(lget(f, "cpu"))[0] + cpu)));
  Shield<SEXP> b(lset(a, "own", Rf_ScalarReal(REAL(lget(f, "own"))[0] + own)));
  Shield<SEXP> c(lset(b, "iters", Rf_ScalarInteger(INTEGER(lget(f, "iters"))[0] + 1)));
  return lset(c, "moves", Rf_ScalarInteger(INTEGER(lget(f, "moves"))[0] +
                                           static_cast<int>(std::nearbyint(n_moves))));
}

bool family_ok(SEXP f) {
  return real_scalar(lget(f, "cpu")) && real_scalar(lget(f, "own")) &&
         int_scalar(lget(f, "iters")) && int_scalar(lget(f, "moves"));
}

// .emc_step_store: the draw into the sample store (fill_samples for a
// standard hierarchy) and the iteration's costs into the allocation window
// (.emc_alloc_record_cost).
bool step_store(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), pars = var(st, "pars"), proposals = var(st, "proposals");
  SEXP samples = lget(pmwgs, "samples");
  SEXP a_half = lget(samples, "a_half"), theta_mu = lget(samples, "theta_mu"),
       theta_var = lget(samples, "theta_var"), alpha = lget(samples, "alpha"),
       subj_ll = lget(samples, "subj_ll");
  SEXP p_a_half = lget(pars, "a_half"), tmu = lget(pars, "tmu"), tvar = lget(pars, "tvar"),
       tvinv = lget(pars, "tvinv");
  if (TYPEOF(a_half) != REALSXP || TYPEOF(theta_mu) != REALSXP ||
      TYPEOF(theta_var) != REALSXP || TYPEOF(alpha) != REALSXP ||
      TYPEOF(subj_ll) != REALSXP || TYPEOF(p_a_half) != REALSXP ||
      TYPEOF(tmu) != REALSXP || TYPEOF(tvar) != REALSXP || tvinv == R_NilValue ||
      TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) ||
      lget(samples, "last_theta_var_inv") == R_NilValue ||
      lget(samples, "idx") == R_NilValue)
    return false;
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  SEXP alloc = lget(pmwgs, "alloc"), w = lget(alloc, "window"), fam = lget(w, "fam");
  SEXP tasks = var(st, "hmc_tasks");
  SEXP r_tasks = R_NilValue;
  if (alloc_on) {
    if (!int_scalar(lget(w, "n")) || !real_scalar(lget(w, "cpu_rest")) ||
        !family_ok(lget(fam, "shear")) || !family_ok(lget(fam, "collapsed")) ||
        !family_ok(lget(fam, "hmc")))
      return false;
    if (tasks != R_NilValue) {
      SEXP cpu_hmc = lget(w, "cpu_hmc"), hmc_att = lget(w, "hmc_att"),
           rows_hmc = lget(w, "rows_hmc");
      if (TYPEOF(tasks) != VECSXP || TYPEOF(cpu_hmc) != REALSXP ||
          TYPEOF(hmc_att) != INTSXP || TYPEOF(rows_hmc) != REALSXP)
        return false;
      r_tasks = emc_alloc_record_tasks_cpp(cpu_hmc, hmc_att, rows_hmc, tasks);
      if (r_tasks == R_NilValue) return false;
    }
  }
  Shield<SEXP> r_keep(r_tasks);

  Shield<SEXP> j(iteration_index(cx, i));
  const int jj = static_cast<int>(Rf_asReal(j));
  const int n_pars = Rf_asInteger(lget(pmwgs, "n_pars"));
  emc_set_last_slice(a_half, jj, p_a_half);
  emc_set_last_slice(theta_mu, jj, tmu);
  emc_set_last_slice(theta_var, jj, tvar);
  emc_set_particle_slice(alpha, subj_ll, proposals, jj, n_pars);
  Shield<SEXP> s1(lset(samples, "last_theta_var_inv", tvinv));
  Shield<SEXP> s2(lset(s1, "idx", j));
  RObject result = lset(pmwgs, "samples", s2);
  if (alloc_on) {
    const double cores = num(var(st, "cores_it"), 1.0);
    const double shear_el = num(var(st, "shear_el"), 0.0), coll_el = num(var(st, "coll_el"), 0.0),
                 hmc_el = num(var(st, "hmc_el"), 0.0);
    const double rest = reserved_cost(clock_now().elapsed - num(var(st, "el_iter"), 0.0) -
                                      shear_el - coll_el - hmc_el, cores);
    RObject wc = lset(w, "n", Rf_ScalarInteger(INTEGER(lget(w, "n"))[0] + 1));
    wc = lset(wc, "cpu_rest", Rf_ScalarReal(REAL(lget(w, "cpu_rest"))[0] + rest));
    Shield<SEXP> famc(Rf_shallow_duplicate(fam));
    double work = 0.0;
    int attempts = 0;
    if (r_tasks != R_NilValue) {
      wc = lset(wc, "cpu_hmc", lget(r_tasks, "cpu_hmc"));
      wc = lset(wc, "hmc_att", lget(r_tasks, "hmc_att"));
      wc = lset(wc, "rows_hmc", lget(r_tasks, "rows_hmc"));
      work = Rf_asReal(lget(r_tasks, "work"));
      attempts = Rf_asInteger(lget(r_tasks, "attempts"));
    }
    const char* fams[3] = {"shear", "collapsed", "hmc"};
    const double m_cpu[3] = {num(var(st, "shear_cpu"), 0.0), num(var(st, "coll_cpu"), 0.0),
                             reserved_cost(hmc_el, cores)};
    const double m_own[3] = {num(var(st, "shear_round_cpu"), 0.0),
                             num(var(st, "coll_round_cpu"), 0.0), work};
    const double m_n[3] = {num(var(st, "shear_rounds"), 0.0), num(var(st, "coll_moves"), 0.0),
                           static_cast<double>(attempts)};
    const R_xlen_t n_fam = XLENGTH(famc);
    SEXP fam_names = Rf_getAttrib(famc, R_NamesSymbol);
    for (int k = 0; k < 3; ++k) {
      if (!R_FINITE(m_n[k]) || m_n[k] <= 0) continue;
      for (R_xlen_t q = 0; q < n_fam; ++q)
        if (std::strcmp(CHAR(STRING_ELT(fam_names, q)), fams[k]) == 0)
          SET_VECTOR_ELT(famc, q, family_add(VECTOR_ELT(famc, q), m_cpu[k], m_own[k], m_n[k]));
    }
    wc = lset(wc, "fam", famc);
    Shield<SEXP> alloc2(lset(alloc, "window", wc));
    result = lset(result, "alloc", alloc2);
  }
  set_var(st, "pmwgs", result);
  return true;
}

typedef bool (*NativeStep)(SEXP, SEXP, int);

NativeStep native_step(const std::string& name) {
  if (name == "store") return step_store;
  if (name == "hmc") return emc_iter::step_hmc;
  if (name == "group") return emc_iter::step_group;
  if (name == "particle") return emc_iter::step_particle;
  if (name == "ensemble") return emc_iter::step_ensemble;
  if (name == "shear") return emc_iter::step_shear;
  if (name == "collapsed") return emc_iter::step_collapsed;
  return NULL;
}

}  // namespace

// Iterations `from` to cx$iter. steps: the R step functions in order, named;
// native: which of them may run here; rejects: the function that takes the
// iteration's rejected updates into the stage's counts.
// [[Rcpp::export]]
void emc_iteration_loop_cpp(Environment st, List cx, int from, List steps,
                            LogicalVector native, Function rejects) {
  const int iter = Rf_asInteger(emc_iter::lget(cx, "iter"));
  const int n_steps = steps.size();
  const CharacterVector names = steps.names();
  std::vector<NativeStep> fast(n_steps, static_cast<NativeStep>(NULL));
  std::vector<bool> is_store(n_steps, false);
  for (int k = 0; k < n_steps; ++k) {
    const std::string nm = as<std::string>(names[k]);
    if (native[k]) fast[k] = native_step(nm);
    is_store[k] = nm == "store";
  }
  // How often each step ran here, left in the state as `native_runs`.
  IntegerVector runs(n_steps);
  runs.names() = names;
  for (int i = from; i <= iter; ++i) {
    for (int k = 0; k < n_steps; ++k) {
      if (fast[k] != NULL && fast[k](st, cx, i)) {
        if (is_store[k]) rejects(st);
        ++runs[k];
        continue;
      }
      const Function step = steps[k];
      step(st, cx, i);
    }
    Rcpp::checkUserInterrupt();
  }
  emc_iter::set_var(st, "native_runs", runs);
}
