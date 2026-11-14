// The iteration's group step (R/iteration_steps.R .emc_step_group) for a
// standard hierarchy without group designs: the group draw given the
// subjects' last states, the ensemble pool members from the same stream, and
// the factor of the group covariance the particle steps use.
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <RcppArmadillo.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
#include <cmath>
#include "iteration_state.h"

using namespace Rcpp;
using namespace emc_iter;

Rcpp::List emc_gibbs_standard_cpp(const arma::mat& alpha, const arma::mat& tvinv0,
                                  const arma::vec& a_half0, const Rcpp::List& prior,
                                  const Rcpp::LogicalVector& is_blocked,
                                  const Rcpp::IntegerVector& par_group);
Rcpp::List emc_ensemble_factors_cpp(const Rcpp::List& var);
Rcpp::List emc_ensemble_draw_pools_cpp(const Rcpp::List& mu, const arma::cube& L,
                                       const arma::vec& logdet, int M1, double tail_weight,
                                       double tail_df, double tail_excess,
                                       const Rcpp::CharacterVector& par_names);

namespace emc_iter {

SEXP chol_factor(SEXP sigma) {
  int p = Rf_nrows(sigma);
  Shield<SEXP> R(Rf_duplicate(sigma));
  double* r = REAL(R);
  for (int j = 0; j < p; ++j)
    for (int i = j + 1; i < p; ++i) r[i + static_cast<R_xlen_t>(j) * p] = 0.0;
  int info = 0;
  F77_CALL(dpotrf)("U", &p, r, &p, &info FCONE);
  if (info != 0) return R_NilValue;
  for (int i = 0; i < p; ++i) if (r[i + static_cast<R_xlen_t>(i) * p] == 0.0) return R_NilValue;
  Shield<SEXP> rooti(Rf_allocMatrix(REALSXP, p, p));
  double* ri = REAL(rooti);
  std::fill(ri, ri + static_cast<R_xlen_t>(p) * p, 0.0);
  for (int i = 0; i < p; ++i) ri[i + static_cast<R_xlen_t>(i) * p] = 1.0;
  const double one = 1.0;
  F77_CALL(dtrsm)("L", "U", "N", "N", &p, &p, &one, r, &p, ri, &p FCONE FCONE FCONE FCONE);
  long double sld = 0.0L;
  for (int i = 0; i < p; ++i) sld += std::log(ri[i + static_cast<R_xlen_t>(i) * p]);
  return List::create(_["R"] = R, _["rooti"] = rooti,
                      _["sum_log_diag"] = static_cast<double>(sld), _["p"] = p,
                      _["repaired"] = false);
}

namespace {

// Wall-clock start of the iteration and the chain's cores, for the
// allocation's costs.
void iteration_costs(SEXP st, SEXP cx, const Clock& started) {
  if (flag(lget(cx, "persistent_workers")))
    set_var(st, "pool_iteration",
            Rf_ScalarInteger(Rf_asInteger(var(st, "pool_iteration")) + 1));
  if (!flag(lget(cx, "alloc_on"))) return;
  SEXP ctl = lget(cx, "core_ctl");
  double cores = num(lget(cx, "n_cores"), 1.0);
  if (ctl != R_NilValue) {
    const double total = num(lget(ctl, "total"), 1.0), slots = num(lget(ctl, "n_slots"), 1.0);
    cores = std::max(1.0, std::floor(total / slots));
  }
  set_var(st, "el_iter", Rf_ScalarReal(started.elapsed));
  set_var(st, "cores_it", Rf_ScalarReal(cores));
}

// The step of a fit without a group level (gibbs_step_single): the prior as
// the group state and the subjects' last states.
bool group_single(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), single = lget(cx, "single");
  SEXP tmu = lget(single, "tmu"), tvar = lget(single, "tvar"), rows = lget(single, "rows");
  SEXP samples = lget(pmwgs, "samples"), s_alpha = lget(samples, "alpha");
  SEXP dims = Rf_getAttrib(s_alpha, R_DimSymbol);
  SEXP seed = lget(lget(pmwgs, "rng"), "gibbs");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(s_alpha) != REALSXP || TYPEOF(dims) != INTSXP ||
      XLENGTH(dims) != 3 || TYPEOF(seed) != INTSXP || TYPEOF(tvar) != REALSXP ||
      !Rf_isMatrix(tvar) || TYPEOF(rows) != STRSXP || !flag(lget(cx, "chol_simple")) ||
      flag(lget(cx, "alloc_on")))
    return false;
  const int p = INTEGER(dims)[0], n = INTEGER(dims)[1];
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (prev < 0 || XLENGTH(s_alpha) < (prev + 1) * slice || XLENGTH(rows) != p ||
      Rf_nrows(tvar) != p || Rf_ncols(tvar) != p)
    return false;
  const Clock started = clock_now();
  RObject group_chol = var(st, "group_chol_it");
  if (static_cast<SEXP>(group_chol) == R_NilValue ||
      !R_compute_identical(lget(group_chol, "ref"), tvar, 0)) {
    Shield<SEXP> f(chol_factor(tvar));
    if (f == R_NilValue) return false;
    group_chol = List::create(_["idx_list"] = lget(cx, "idx_list_blk"), _["ref"] = tvar,
                              _["marginal_idx"] = lget(cx, "marginal_idx_blk"),
                              _["include_full"] = false, _["full"] = R_NilValue,
                              _["f"] = List::create(f));
  }
  Shield<SEXP> alpha(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(s_alpha) + prev * slice, REAL(s_alpha) + (prev + 1) * slice, REAL(alpha));
  Shield<SEXP> dn(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, rows);
  Rf_setAttrib(alpha, R_DimNamesSymbol, dn);
  set_global_seed(seed);
  set_var(st, "pars", List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["alpha"] = alpha));
  set_var(st, "ens_pools", R_NilValue);
  set_var(st, "group_chol_it", group_chol);
  iteration_costs(st, cx, started);
  return true;
}

}  // namespace

bool step_group(SEXP st, SEXP cx, int i) {
  if (lget(cx, "single") != R_NilValue) return group_single(st, cx, i);
  SEXP pmwgs = var(st, "pmwgs");
  SEXP gibbs = lget(cx, "gibbs"), pools = lget(cx, "pools");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(gibbs) != VECSXP || pools == R_NilValue ||
      !flag(lget(cx, "chol_simple")))
    return false;
  SEXP samples = lget(pmwgs, "samples"), rng = lget(pmwgs, "rng"), prior = lget(pmwgs, "prior");
  SEXP s_alpha = lget(samples, "alpha"), s_a_half = lget(samples, "a_half"),
       tvinv = lget(samples, "last_theta_var_inv"), idx = lget(samples, "idx");
  SEXP dims = Rf_getAttrib(s_alpha, R_DimSymbol);
  SEXP seed = lget(rng, "gibbs");
  if (TYPEOF(s_alpha) != REALSXP || TYPEOF(dims) != INTSXP || XLENGTH(dims) != 3 ||
      TYPEOF(s_a_half) != REALSXP || TYPEOF(tvinv) != REALSXP || !Rf_isMatrix(tvinv) ||
      Rf_length(idx) != 1 || TYPEOF(seed) != INTSXP || TYPEOF(prior) != VECSXP)
    return false;
  const int p = INTEGER(dims)[0], n = INTEGER(dims)[1];
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t at = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (prev < 0 || at < 0 || XLENGTH(s_alpha) < (prev + 1) * slice ||
      XLENGTH(s_a_half) < (at + 1) * p || Rf_nrows(tvinv) != p || Rf_ncols(tvinv) != p)
    return false;
  const Clock started = clock_now();

  // The subjects' last states, named as the store names them.
  Shield<SEXP> alpha_prev(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(s_alpha) + prev * slice, REAL(s_alpha) + (prev + 1) * slice, REAL(alpha_prev));
  SEXP dn = Rf_getAttrib(s_alpha, R_DimNamesSymbol);
  if (dn != R_NilValue) {
    Shield<SEXP> dn2(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn2, 0, VECTOR_ELT(dn, 0));
    SET_VECTOR_ELT(dn2, 1, VECTOR_ELT(dn, 1));
    SEXP dnn = Rf_getAttrib(dn, R_NamesSymbol);
    if (dnn != R_NilValue) {
      Shield<SEXP> dnn2(Rf_allocVector(STRSXP, 2));
      SET_STRING_ELT(dnn2, 0, STRING_ELT(dnn, 0));
      SET_STRING_ELT(dnn2, 1, STRING_ELT(dnn, 1));
      Rf_setAttrib(dn2, R_NamesSymbol, dnn2);
    }
    Rf_setAttrib(alpha_prev, R_DimNamesSymbol, dn2);
  }

  set_global_seed(seed);
  GetRNGstate();
  List res, drawn;
  RObject factors = var(st, "ens_factors");
  const bool pooled = TYPEOF(pools) == VECSXP;
  bool ok = true;
  try {
    const arma::mat a(REAL(alpha_prev), p, n, false, true);
    const arma::mat tv(REAL(tvinv), p, p, false, true);
    const arma::vec ah(REAL(s_a_half) + at * p, p);
    res = emc_gibbs_standard_cpp(a, tv, ah, List(prior), LogicalVector(lget(gibbs, "blocked")),
                                 IntegerVector(lget(gibbs, "group")));
    if (pooled) {
      if (static_cast<SEXP>(factors) == R_NilValue) {
        const List f = emc_ensemble_factors_cpp(List(lget(pools, "var")));
        factors = List::create(_["L"] = f["L"], _["logdet"] = f["logdet"]);
      }
      drawn = emc_ensemble_draw_pools_cpp(
          List(lget(pools, "mu")), as<arma::cube>(lget(factors, "L")),
          as<arma::vec>(lget(factors, "logdet")), Rf_asInteger(lget(pools, "size")),
          Rf_asReal(lget(pools, "tail_weight")), Rf_asReal(lget(pools, "tail_df")),
          Rf_asReal(lget(pools, "tail_excess")), CharacterVector(lget(pools, "par_names")));
    }
  } catch (std::exception&) {
    ok = false;
  }
  PutRNGstate();
  if (!ok) return false;
  Shield<SEXP> seed_after(global_seed());

  SEXP tmu = res["tmu"], tvar = res["tvar"];
  if (TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != p || TYPEOF(tvar) != REALSXP ||
      !Rf_isMatrix(tvar))
    return false;
  Shield<SEXP> subj_mu(Rf_allocMatrix(REALSXP, p, n));
  for (int s = 0; s < n; ++s)
    std::copy(REAL(tmu), REAL(tmu) + p, REAL(subj_mu) + static_cast<R_xlen_t>(s) * p);
  List pars = List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["tvinv"] = res["tvinv"],
                           _["a_half"] = res["a_half"], _["subj_mu"] = subj_mu,
                           _["alpha"] = alpha_prev);

  RObject group_chol = var(st, "group_chol_it");
  if (static_cast<SEXP>(group_chol) == R_NilValue ||
      !R_compute_identical(lget(group_chol, "ref"), tvar, 0)) {
    Shield<SEXP> f(chol_factor(tvar));
    if (f == R_NilValue) return false;
    group_chol = List::create(_["idx_list"] = lget(cx, "idx_list_blk"), _["ref"] = tvar,
                              _["marginal_idx"] = lget(cx, "marginal_idx_blk"),
                              _["include_full"] = false, _["full"] = R_NilValue,
                              _["f"] = List::create(f));
  }

  RObject ens_pools;
  if (pooled) {
    ens_pools = List::create(_["rows"] = drawn["rows"], _["lq"] = drawn["lq"],
                             _["ll"] = R_NilValue, _["L"] = lget(factors, "L"),
                             _["mu"] = drawn["mu"], _["logdet"] = lget(factors, "logdet"),
                             _["c"] = drawn["c"]);
    set_var(st, "ens_factors", factors);
  }
  Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
  set_var(st, "pmwgs", lset(pmwgs, "rng", rng2));
  set_var(st, "pars", pars);
  set_var(st, "ens_pools", ens_pools);
  set_var(st, "group_chol_it", group_chol);
  iteration_costs(st, cx, started);
  return true;
}

}  // namespace emc_iter
