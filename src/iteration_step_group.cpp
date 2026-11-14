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
#include <cstring>
#include <vector>
#include "iteration_state.h"

using namespace Rcpp;
using namespace emc_iter;

Rcpp::List emc_gibbs_standard_cpp(const arma::mat& alpha, const arma::mat& tvinv0,
                                  const arma::vec& a_half0, const Rcpp::List& prior,
                                  const Rcpp::LogicalVector& is_blocked,
                                  const Rcpp::IntegerVector& par_group);
Rcpp::List emc_gibbs_diag_gamma_cpp(const arma::mat& alpha, const arma::vec& tvinv0,
                                    const Rcpp::List& prior);
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

// The step of a standard hierarchy beside nuisance parameters without a
// group level: the hierarchy's draw on its own rows, the nuisance prior mean
// on the others, and the two merged for the subjects (merge_group_level).
bool group_nuisance(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), nz = lget(cx, "nuis");
  SEXP is = lget(nz, "is"), tmu_nuis = lget(nz, "tmu");
  SEXP samples = lget(pmwgs, "samples"), rng = lget(pmwgs, "rng"), prior = lget(pmwgs, "prior");
  SEXP s_alpha = lget(samples, "alpha"), s_a_half = lget(samples, "a_half"),
       tvinv = lget(samples, "last_theta_var_inv"), idx = lget(samples, "idx");
  SEXP dims = Rf_getAttrib(s_alpha, R_DimSymbol), dn = Rf_getAttrib(s_alpha, R_DimNamesSymbol);
  SEXP seed = lget(rng, "gibbs"), nuis = lget(pmwgs, "sampler_nuis");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(is) != LGLSXP || TYPEOF(tmu_nuis) != REALSXP ||
      TYPEOF(s_alpha) != REALSXP || TYPEOF(dims) != INTSXP || XLENGTH(dims) != 3 ||
      TYPEOF(s_a_half) != REALSXP || TYPEOF(tvinv) != REALSXP || !Rf_isMatrix(tvinv) ||
      Rf_length(idx) != 1 || TYPEOF(seed) != INTSXP || TYPEOF(prior) != VECSXP ||
      TYPEOF(nuis) != VECSXP || dn == R_NilValue || XLENGTH(dn) != 3 ||
      VECTOR_ELT(dn, 0) == R_NilValue)
    return false;
  const int p = INTEGER(dims)[0], n = INTEGER(dims)[1];
  if (XLENGTH(is) != p || n < 2 || p < 2) return false;
  std::vector<int> main_rows, nuis_rows;
  for (int k = 0; k < p; ++k) (LOGICAL(is)[k] == TRUE ? nuis_rows : main_rows).push_back(k);
  const int pm = static_cast<int>(main_rows.size()), pn = static_cast<int>(nuis_rows.size());
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t at = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (pm < 1 || pn < 1 || XLENGTH(tmu_nuis) != pn || prev < 0 || at < 0 ||
      XLENGTH(s_alpha) < (prev + 1) * slice || XLENGTH(s_a_half) < (at + 1) * pm ||
      Rf_nrows(tvinv) != pm || Rf_ncols(tvinv) != pm)
    return false;
  const Clock started = clock_now();

  // Every parameter's last states, and the hierarchy's rows of them.
  Shield<SEXP> alpha_all(Rf_allocMatrix(REALSXP, p, n));
  std::copy(REAL(s_alpha) + prev * slice, REAL(s_alpha) + (prev + 1) * slice, REAL(alpha_all));
  Shield<SEXP> dn_all(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn_all, 0, VECTOR_ELT(dn, 0));
  SET_VECTOR_ELT(dn_all, 1, VECTOR_ELT(dn, 1));
  Shield<SEXP> dnn2(Rf_allocVector(STRSXP, 2));
  SEXP dnn = Rf_getAttrib(dn, R_NamesSymbol);
  if (dnn != R_NilValue) {
    SET_STRING_ELT(dnn2, 0, STRING_ELT(dnn, 0));
    SET_STRING_ELT(dnn2, 1, STRING_ELT(dnn, 1));
    Rf_setAttrib(dn_all, R_NamesSymbol, dnn2);
  }
  Rf_setAttrib(alpha_all, R_DimNamesSymbol, dn_all);
  Shield<SEXP> alpha_main(Rf_allocMatrix(REALSXP, pm, n));
  for (int s = 0; s < n; ++s)
    for (int k = 0; k < pm; ++k)
      REAL(alpha_main)[k + static_cast<R_xlen_t>(s) * pm] =
          REAL(alpha_all)[main_rows[k] + static_cast<R_xlen_t>(s) * p];
  {
    Shield<SEXP> rows(Rf_allocVector(STRSXP, pm));
    for (int k = 0; k < pm; ++k)
      SET_STRING_ELT(rows, k, STRING_ELT(VECTOR_ELT(dn, 0), main_rows[k]));
    Shield<SEXP> dn_main(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn_main, 0, rows);
    SET_VECTOR_ELT(dn_main, 1, VECTOR_ELT(dn, 1));
    if (dnn != R_NilValue) Rf_setAttrib(dn_main, R_NamesSymbol, dnn2);
    Rf_setAttrib(alpha_main, R_DimNamesSymbol, dn_main);
  }

  set_global_seed(seed);
  GetRNGstate();
  List res;
  bool ok = true;
  try {
    const arma::mat a(REAL(alpha_main), pm, n, false, true);
    const arma::mat tv(REAL(tvinv), pm, pm, false, true);
    const arma::vec ah(REAL(s_a_half) + at * pm, pm);
    res = emc_gibbs_standard_cpp(a, tv, ah, List(prior), LogicalVector(lget(nz, "blocked")),
                                 IntegerVector(lget(nz, "group")));
  } catch (std::exception&) {
    ok = false;
  }
  PutRNGstate();
  if (!ok) return false;
  Shield<SEXP> seed_after(global_seed());
  SEXP tmu = res["tmu"], tvar = res["tvar"];
  if (TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != pm || TYPEOF(tvar) != REALSXP ||
      !Rf_isMatrix(tvar) || Rf_nrows(tvar) != pm)
    return false;
  Shield<SEXP> subj_mu(Rf_allocMatrix(REALSXP, pm, n));
  for (int s = 0; s < n; ++s)
    std::copy(REAL(tmu), REAL(tmu) + pm, REAL(subj_mu) + static_cast<R_xlen_t>(s) * pm);
  List pars = List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["tvinv"] = res["tvinv"],
                           _["a_half"] = res["a_half"], _["subj_mu"] = subj_mu,
                           _["alpha"] = alpha_main);

  // The merged state: the nuisance rows carry their prior mean and no
  // covariance.
  Shield<SEXP> tmu_c(Rf_allocVector(REALSXP, p)), tvar_c(Rf_allocMatrix(REALSXP, p, p)),
      mu_c(Rf_allocMatrix(REALSXP, p, n));
  std::fill(REAL(tvar_c), REAL(tvar_c) + static_cast<R_xlen_t>(p) * p, 0.0);
  for (int k = 0; k < pm; ++k) REAL(tmu_c)[main_rows[k]] = REAL(tmu)[k];
  for (int k = 0; k < pn; ++k) REAL(tmu_c)[nuis_rows[k]] = REAL(tmu_nuis)[k];
  for (int c = 0; c < pm; ++c)
    for (int r = 0; r < pm; ++r)
      REAL(tvar_c)[main_rows[r] + static_cast<R_xlen_t>(main_rows[c]) * p] =
          REAL(tvar)[r + static_cast<R_xlen_t>(c) * pm];
  for (int s = 0; s < n; ++s)
    for (int k = 0; k < p; ++k)
      REAL(mu_c)[k + static_cast<R_xlen_t>(s) * p] = REAL(tmu_c)[k];
  List pars_comb = List::create(_["tmu"] = tmu_c, _["tvar"] = tvar_c, _["subj_mu"] = mu_c,
                                _["alpha"] = alpha_all);

  RObject group_chol = var(st, "group_chol_it");
  if (static_cast<SEXP>(group_chol) == R_NilValue ||
      !R_compute_identical(lget(group_chol, "ref"), tvar_c, 0)) {
    static Function build = package_function("build_group_chol_cache");
    group_chol = build(tvar_c, lget(cx, "idx_list_blk"), lget(cx, "marginal_idx_blk"),
                       _["include_full"] = lget(cx, "include_full"));
  }

  // The nuisance sampler's store holds only its place in the chain.
  const R_xlen_t n_f = XLENGTH(nuis);
  SEXP f_names = Rf_getAttrib(nuis, R_NamesSymbol);
  std::vector<R_xlen_t> keep;
  for (R_xlen_t k = 0; k < n_f; ++k)
    if (std::strcmp(CHAR(STRING_ELT(f_names, k)), "samples") != 0) keep.push_back(k);
  Shield<SEXP> nuis2(Rf_allocVector(VECSXP, keep.size() + 1)),
      names2(Rf_allocVector(STRSXP, keep.size() + 1));
  for (size_t k = 0; k < keep.size(); ++k) {
    SET_VECTOR_ELT(nuis2, k, VECTOR_ELT(nuis, keep[k]));
    SET_STRING_ELT(names2, k, STRING_ELT(f_names, keep[k]));
  }
  SET_VECTOR_ELT(nuis2, keep.size(), List::create(_["idx"] = static_cast<SEXP>(j)));
  SET_STRING_ELT(names2, keep.size(), Rf_mkChar("samples"));
  Rf_setAttrib(nuis2, R_NamesSymbol, names2);

  Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
  RObject sampler = lset(pmwgs, "rng", rng2);
  sampler = lset(sampler, "sampler_nuis", nuis2);
  set_var(st, "pmwgs", sampler);
  set_var(st, "pars", pars);
  set_var(st, "pars_comb", pars_comb);
  set_var(st, "ens_pools", R_NilValue);
  set_var(st, "group_chol_it", group_chol);
  iteration_costs(st, cx, started);
  return true;
}

}  // namespace

bool step_group(SEXP st, SEXP cx, int i) {
  if (lget(cx, "single") != R_NilValue) return group_single(st, cx, i);
  if (lget(cx, "nuis") != R_NilValue) return group_nuisance(st, cx, i);
  SEXP pmwgs = var(st, "pmwgs");
  SEXP gibbs = lget(cx, "gibbs"), pools = lget(cx, "pools");
  // A diagonal-gamma hierarchy: independent variances with gamma precisions.
  SEXP dgamma = lget(cx, "dgamma");
  const bool dg = TYPEOF(dgamma) == VECSXP;
  if (TYPEOF(pmwgs) != VECSXP || (!dg && TYPEOF(gibbs) != VECSXP) || pools == R_NilValue ||
      !flag(lget(cx, "chol_simple")))
    return false;
  SEXP samples = lget(pmwgs, "samples"), rng = lget(pmwgs, "rng");
  SEXP prior = dg ? lget(dgamma, "prior") : lget(pmwgs, "prior");
  SEXP s_alpha = lget(samples, "alpha"), s_a_half = lget(samples, "a_half"),
       tvinv = lget(samples, "last_theta_var_inv"), idx = lget(samples, "idx");
  SEXP dims = Rf_getAttrib(s_alpha, R_DimSymbol);
  SEXP seed = lget(rng, "gibbs");
  if (TYPEOF(s_alpha) != REALSXP || TYPEOF(dims) != INTSXP || XLENGTH(dims) != 3 ||
      (!dg && TYPEOF(s_a_half) != REALSXP) || TYPEOF(tvinv) != REALSXP ||
      !Rf_isMatrix(tvinv) || Rf_length(idx) != 1 || TYPEOF(seed) != INTSXP ||
      TYPEOF(prior) != VECSXP)
    return false;
  const int p = INTEGER(dims)[0], n = INTEGER(dims)[1];
  Shield<SEXP> j(iteration_index(cx, i));
  const R_xlen_t prev = static_cast<R_xlen_t>(Rf_asReal(j)) - 2;
  const R_xlen_t at = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * n;
  if (prev < 0 || at < 0 || XLENGTH(s_alpha) < (prev + 1) * slice ||
      (!dg && XLENGTH(s_a_half) < (at + 1) * p) || Rf_nrows(tvinv) != p ||
      Rf_ncols(tvinv) != p)
    return false;
  if (dg && (TYPEOF(lget(dgamma, "rows")) != STRSXP || XLENGTH(lget(dgamma, "rows")) != p))
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
    if (dg) {
      const arma::vec prec = tv.diag();
      res = emc_gibbs_diag_gamma_cpp(a, prec, List(prior));
    } else {
      const arma::vec ah(REAL(s_a_half) + at * p, p);
      res = emc_gibbs_standard_cpp(a, tv, ah, List(prior),
                                   LogicalVector(lget(gibbs, "blocked")),
                                   IntegerVector(lget(gibbs, "group")));
    }
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

  RObject tmu_keep, tvar_keep;
  List pars;
  if (dg) {
    // gibbs_step_diag_gamma: the named means and the diagonal matrices of
    // the variances and precisions.
    SEXP mu = res["tmu"], prec = res["tvinv"];
    if (TYPEOF(mu) != REALSXP || XLENGTH(mu) != p || TYPEOF(prec) != REALSXP ||
        XLENGTH(prec) != p || Rf_getAttrib(prec, R_DimSymbol) != R_NilValue)
      return false;
    tmu_keep = Rf_duplicate(mu);
    Rf_setAttrib(tmu_keep, R_NamesSymbol, lget(dgamma, "rows"));
    Shield<SEXP> var_m(Rf_allocMatrix(REALSXP, p, p)), prec_m(Rf_allocMatrix(REALSXP, p, p));
    std::fill(REAL(var_m), REAL(var_m) + static_cast<R_xlen_t>(p) * p, 0.0);
    std::fill(REAL(prec_m), REAL(prec_m) + static_cast<R_xlen_t>(p) * p, 0.0);
    for (int k = 0; k < p; ++k) {
      volatile double v = 1.0 / REAL(prec)[k];
      REAL(var_m)[k + static_cast<R_xlen_t>(k) * p] = v;
      REAL(prec_m)[k + static_cast<R_xlen_t>(k) * p] = REAL(prec)[k];
    }
    tvar_keep = var_m;
    pars = List::create(_["tmu"] = tmu_keep, _["tvar"] = tvar_keep, _["tvinv"] = prec_m,
                        _["alpha"] = alpha_prev);
  } else {
    tmu_keep = res["tmu"];
    tvar_keep = res["tvar"];
  }
  SEXP tmu = tmu_keep, tvar = tvar_keep;
  if (TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != p || TYPEOF(tvar) != REALSXP ||
      !Rf_isMatrix(tvar))
    return false;
  if (!dg) {
    Shield<SEXP> subj_mu(Rf_allocMatrix(REALSXP, p, n));
    for (int s = 0; s < n; ++s)
      std::copy(REAL(tmu), REAL(tmu) + p, REAL(subj_mu) + static_cast<R_xlen_t>(s) * p);
    pars = List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["tvinv"] = res["tvinv"],
                        _["a_half"] = res["a_half"], _["subj_mu"] = subj_mu,
                        _["alpha"] = alpha_prev);
  }

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
