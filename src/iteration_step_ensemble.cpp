// The iteration's ensemble step (R/iteration_steps.R .emc_step_ensemble) for
// hierarchies whose Gibbs rounds run natively (standard, diagonal-gamma and
// the latent-factor ones): each subject's pool with its
// current state at a uniform position, the rounds on the group state and
// the pool positions, the selected states, and the group state's jump into
// the allocation window.
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <RcppArmadillo.h>
#include <R_ext/Random.h>
#include <cfloat>
#include <cmath>
#include <vector>
#include "iteration_latent.h"
#include "iteration_state.h"

using namespace Rcpp;
using namespace emc_iter;

Rcpp::List emc_ensemble_stack_cpp(const Rcpp::List& rows, const Rcpp::List& lq,
                                  const arma::mat& ll_new, const arma::mat& proposals,
                                  const Rcpp::IntegerVector& k0, const arma::cube& L,
                                  const arma::mat& mu, const arma::vec& logdet, double c,
                                  double tail_weight, double tail_df);
Rcpp::List emc_ensemble_rounds_cpp(const arma::mat& pool_stack, const arma::vec& base_w,
                                   const Rcpp::IntegerVector& k0, int M, int rounds,
                                   int type, const Rcpp::List& state, const Rcpp::List& prior,
                                   const Rcpp::LogicalVector& is_blocked,
                                   const Rcpp::IntegerVector& par_group);

SEXP emc_ensemble_rounds_latent_cpp(const Rcpp::NumericMatrix& pool_stack,
                                    const Rcpp::NumericVector& base_w,
                                    const Rcpp::IntegerVector& k0, int M, int rounds,
                                    const Rcpp::List& state, SEXP free_entry,
                                    const Rcpp::List& prior, SEXP sem);

namespace {

// as.numeric(x) shaped nr x nc (a plain vector with nc < 0), for a state the
// latent hierarchies' rounds read.
SEXP numeric_part(SEXP x, int nr, int nc) {
  const R_xlen_t len = nc < 0 ? nr : static_cast<R_xlen_t>(nr) * nc;
  if (TYPEOF(x) != REALSXP || XLENGTH(x) != len) throw emc_iter::LatentDecline();
  Shield<SEXP> out(nc < 0 ? Rf_allocVector(REALSXP, len) : Rf_allocMatrix(REALSXP, nr, nc));
  std::copy(REAL(x), REAL(x) + len, REAL(out));
  return out;
}

// The means as a one-row matrix with the parameters' names on its columns.
SEXP named_row(SEXP mu, SEXP names) {
  const int p = Rf_length(mu);
  Shield<SEXP> out(Rf_allocMatrix(REALSXP, 1, p));
  std::copy(REAL(mu), REAL(mu) + p, REAL(out));
  Shield<SEXP> dn(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 1, names);
  Rf_setAttrib(out, R_DimNamesSymbol, dn);
  return out;
}

}  // namespace

namespace emc_iter {

List latent_state(SEXP pars, int kind, int n, int p, int k, int n_cov) {
  if (kind == 0)
    return List::create(_["eta"] = numeric_part(lget(pars, "eta"), n, k),
                        _["lambda"] = numeric_part(lget(pars, "lambda_untransf"), p, k),
                        _["mu"] = numeric_part(lget(pars, "tmu"), p, -1),
                        _["sig_err_inv"] = numeric_part(lget(pars, "sig_err_inv"), p, -1),
                        _["psi_inv"] = numeric_part(lget(pars, "psi_inv"), k, -1));
  if (kind == 1)
    return List::create(_["psi"] = numeric_part(lget(pars, "psi"), p, k),
                        _["lambda"] = numeric_part(lget(pars, "lambda"), p, k),
                        _["mu"] = numeric_part(lget(pars, "tmu"), p, -1),
                        _["epsilon_inv"] = numeric_part(lget(pars, "epsilon_inv"), p, -1),
                        _["delta"] = numeric_part(lget(pars, "delta"), k, -1));
  return List::create(_["lambda"] = numeric_part(lget(pars, "lambda"), p, k),
                      _["B"] = numeric_part(lget(pars, "B"), k, k),
                      _["K"] = numeric_part(lget(pars, "K"), p, n_cov),
                      _["G"] = numeric_part(lget(pars, "G"), k, n_cov),
                      _["delta_inv"] = numeric_part(lget(pars, "delta_inv"), k, k),
                      _["mu"] = numeric_part(lget(pars, "tmu"), p, -1),
                      _["epsilon_inv"] = numeric_part(lget(pars, "epsilon_inv"), p, -1));
}

List latent_pars(const List& state, int kind, SEXP names, SEXP delta_form, SEXP eta_rows,
                 SEXP alpha) {
  if (kind == 0)
    return List::create(_["tmu"] = named_row(state["mu"], names), _["tvar"] = state["var"],
                        _["lambda_untransf"] = state["lambda"],
                        _["lambda"] = state["lambda_scaled"], _["eta"] = state["eta"],
                        _["sig_err_inv"] = state["sig_err_inv"],
                        _["psi_inv"] = state["psi_inv"], _["alpha"] = alpha);
  if (kind == 1) {
    // The factors' shrinkage keeps the store's form.
    SEXP values = state["delta"];
    if (TYPEOF(delta_form) != REALSXP || XLENGTH(delta_form) != XLENGTH(values))
      throw emc_iter::LatentDecline();
    Shield<SEXP> delta(Rf_duplicate(delta_form));
    std::copy(REAL(values), REAL(values) + XLENGTH(values), REAL(delta));
    RObject eta = state["eta"];
    if (eta_rows != R_NilValue) {
      eta = Rf_duplicate(eta);
      Shield<SEXP> dn(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(dn, 0, eta_rows);
      Rf_setAttrib(eta, R_DimNamesSymbol, dn);
    }
    return List::create(_["tmu"] = named_row(state["mu"], names), _["tvar"] = state["var"],
                        _["lambda"] = state["lambda"], _["eta"] = eta,
                        _["epsilon_inv"] = state["epsilon_inv"], _["psi"] = state["psi"],
                        _["alpha"] = alpha, _["delta"] = delta);
  }
  return List::create(_["tmu"] = state["tmu"], _["lambda"] = state["lambda"],
                      _["eta"] = state["eta"], _["B"] = state["B"], _["K"] = state["K"],
                      _["G"] = state["G"], _["epsilon_inv"] = state["epsilon_inv"],
                      _["delta_inv"] = state["delta_inv"], _["alpha"] = alpha,
                      _["subj_mu"] = state["subj_mu"], _["subj_var"] = state["subj_var"],
                      _["pop_mean"] = state["pop_mean"], _["tvar"] = state["tvar"]);
}

}  // namespace emc_iter

namespace {

// mean() of a numeric vector.
double r_mean(const double* x, int n) {
  long double s = 0.0L;
  for (int i = 0; i < n; ++i) s += x[i];
  s /= n;
  long double t = 0.0L;
  for (int i = 0; i < n; ++i) t += (x[i] - s);
  s += t / n;
  return static_cast<double>(s);
}

// .emc_alloc_lsd() of a p x p matrix stored at `v`.
void log_sd(const double* v, int p, double* out) {
  for (int k = 0; k < p; ++k) {
    double d = v[k + static_cast<R_xlen_t>(k) * p];
    if (d < DBL_MIN) d = DBL_MIN;
    volatile double l = std::log(d);
    out[k] = 0.5 * l;
  }
}

}  // namespace

namespace emc_iter {

bool step_ensemble(SEXP st, SEXP cx, int i) {
  SEXP pmwgs = var(st, "pmwgs"), pars = var(st, "pars"), proposals = var(st, "proposals");
  SEXP pools = var(st, "ens_pools"), gibbs = lget(cx, "gibbs"), cfg = lget(cx, "pools");
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  if (pools == R_NilValue) return !alloc_on;
  // A diagonal-gamma hierarchy's rounds draw independent precisions.
  SEXP dgamma = lget(cx, "dgamma");
  const bool dg = TYPEOF(dgamma) == VECSXP;
  // A latent-factor hierarchy's rounds run on its own group draw.
  SEXP latent = lget(cx, "latent");
  const bool lat = TYPEOF(latent) == VECSXP;
  if (TYPEOF(pmwgs) != VECSXP || (!dg && !lat && TYPEOF(gibbs) != VECSXP) ||
      TYPEOF(cfg) != VECSXP ||
      TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals))
    return false;
  const int p = Rf_nrows(proposals) - 1, n = Rf_ncols(proposals);
  SEXP rows = lget(pools, "rows"), ll = lget(pools, "ll");
  SEXP rng = lget(pmwgs, "rng"), seed = lget(rng, "gibbs");
  SEXP prior = dg ? lget(dgamma, "prior") : lget(pmwgs, "prior");
  SEXP samples = lget(pmwgs, "samples");
  if (TYPEOF(rows) != VECSXP || XLENGTH(rows) != n || TYPEOF(ll) != VECSXP ||
      XLENGTH(ll) != n || TYPEOF(seed) != INTSXP || !Rf_isMatrix(VECTOR_ELT(rows, 0)) ||
      Rf_asInteger(lget(pmwgs, "n_pars")) != p || Rf_asInteger(lget(pmwgs, "n_subjects")) != n)
    return false;
  const int M1 = Rf_nrows(VECTOR_ELT(rows, 0)), M = Rf_asInteger(lget(cfg, "size")) + 1;
  if (M1 + 1 != M) return false;
  NumericMatrix ll_new(M1, n);
  for (int s = 0; s < n; ++s) {
    SEXP l = VECTOR_ELT(ll, s);
    if (TYPEOF(l) != REALSXP || XLENGTH(l) != M1) return false;
    for (int k = 0; k < M1; ++k) {
      const double x = REAL(l)[k];
      if (ISNAN(x) || x == R_PosInf) return false;
      ll_new(k, s) = x;
    }
  }
  SEXP tvinv = lget(pars, "tvinv"), a_half = lget(pars, "a_half");
  if (!lat && (TYPEOF(tvinv) != REALSXP || (!dg && TYPEOF(a_half) != REALSXP))) return false;
  const int kind = lat ? Rf_asInteger(lget(latent, "kind")) : -1;
  SEXP par_names = lat ? lget(latent, "par_names") : R_NilValue;
  List start;
  if (lat) {
    if (kind < 0 || kind > 2 || TYPEOF(par_names) != STRSXP || XLENGTH(par_names) != p ||
        TYPEOF(pars) != VECSXP)
      return false;
    try {
      start = latent_state(pars, kind, n, p, Rf_asInteger(lget(latent, "k")),
                           Rf_asInteger(lget(latent, "n_cov")));
    } catch (LatentDecline&) {
      return false;
    }
  }
  if (dg && (!Rf_isMatrix(tvinv) || Rf_nrows(tvinv) != p || Rf_ncols(tvinv) != p ||
             TYPEOF(prior) != VECSXP || TYPEOF(lget(dgamma, "rows")) != STRSXP ||
             XLENGTH(lget(dgamma, "rows")) != p))
    return false;
  SEXP alloc = lget(pmwgs, "alloc"), w = lget(alloc, "window");
  SEXP G = lget(w, "G"), nG = lget(w, "nG"), idx = lget(samples, "idx");
  SEXP theta_mu = lget(samples, "theta_mu"), theta_var = lget(samples, "theta_var");
  if (alloc_on) {
    if (TYPEOF(G) != REALSXP || !Rf_isMatrix(G) || Rf_nrows(G) != 2 * p ||
        Rf_ncols(G) != 2 * p || TYPEOF(nG) != INTSXP || XLENGTH(nG) != 1 ||
        TYPEOF(theta_mu) != REALSXP || TYPEOF(theta_var) != REALSXP || Rf_length(idx) != 1)
      return false;
    const R_xlen_t at = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
    if (at < 0 || XLENGTH(theta_mu) < (at + 1) * p ||
        XLENGTH(theta_var) < (at + 1) * static_cast<R_xlen_t>(p) * p)
      return false;
  }

  set_global_seed(seed);
  GetRNGstate();
  IntegerVector k0(n);
  for (int s = 0; s < n; ++s) k0[s] = static_cast<int>(R_unif_index(M) + 1);
  List stacked, res;
  bool ok = true;
  try {
    const arma::mat lln(ll_new.begin(), M1, n, false, true);
    const arma::mat props(REAL(proposals), p + 1, n, false, true);
    stacked = emc_ensemble_stack_cpp(List(rows), List(lget(pools, "lq")), lln, props, k0,
                                     as<arma::cube>(lget(pools, "L")),
                                     as<arma::mat>(lget(pools, "mu")),
                                     as<arma::vec>(lget(pools, "logdet")),
                                     Rf_asReal(lget(pools, "c")),
                                     Rf_asReal(lget(cfg, "tail_weight")),
                                     Rf_asReal(lget(cfg, "tail_df")));
    if (lat) {
      SEXP out = emc_ensemble_rounds_latent_cpp(
          as<NumericMatrix>(stacked["pool_stack"]), as<NumericVector>(stacked["base_w"]), k0,
          M, Rf_asInteger(lget(cfg, "rounds")), start, lget(latent, "free"),
          List(lget(latent, "prior")), lget(latent, "sem"));
      if (out == R_NilValue) ok = false;
      else res = out;
    } else if (dg) {
      NumericVector prec(p);
      for (int k = 0; k < p; ++k) prec[k] = REAL(tvinv)[k + static_cast<R_xlen_t>(k) * p];
      res = emc_ensemble_rounds_cpp(as<arma::mat>(stacked["pool_stack"]),
                                    as<arma::vec>(stacked["base_w"]), k0, M,
                                    Rf_asInteger(lget(cfg, "rounds")), 1,
                                    List::create(_["tvinv"] = prec), List(prior),
                                    LogicalVector(0), IntegerVector(0));
    } else {
      NumericVector ah(XLENGTH(a_half));
      std::copy(REAL(a_half), REAL(a_half) + XLENGTH(a_half), ah.begin());
      const List state = List::create(_["tvinv"] = tvinv, _["a_half"] = ah);
      res = emc_ensemble_rounds_cpp(as<arma::mat>(stacked["pool_stack"]),
                                    as<arma::vec>(stacked["base_w"]), k0, M,
                                    Rf_asInteger(lget(cfg, "rounds")), 0, state, List(prior),
                                    LogicalVector(lget(gibbs, "blocked")),
                                    IntegerVector(lget(gibbs, "group")));
    }
  } catch (std::exception&) {
    ok = false;
  }
  PutRNGstate();
  if (!ok) return false;
  Shield<SEXP> seed_after(global_seed());

  const IntegerVector selected = res["selected"];
  const NumericVector wess = res["wess"];
  const NumericMatrix pool_stack = stacked["pool_stack"];
  const NumericVector ll_stack = stacked["ll_stack"];
  RObject tmu_keep, tvar_keep, tvinv_keep;
  if (lat) {
    const List state = res["state"];
    tmu_keep = kind == 2 ? state["tmu"] : state["mu"];
    tvar_keep = kind == 2 ? state["tvar"] : state["var"];
  } else {
    tmu_keep = res["tmu"];
  }
  if (lat) {
  } else if (dg) {
    // The named means and the diagonal matrices of the variances and
    // precisions, as the hierarchy's own draw gives them.
    SEXP prec = res["tvinv"];
    if (TYPEOF(tmu_keep) != REALSXP || XLENGTH(tmu_keep) != p || TYPEOF(prec) != REALSXP ||
        XLENGTH(prec) != p)
      return false;
    tmu_keep = Rf_duplicate(tmu_keep);
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
    tvinv_keep = prec_m;
  } else {
    tvar_keep = res["tvar"];
  }
  SEXP tmu = tmu_keep, tvar = tvar_keep;
  if (TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != p || TYPEOF(tvar) != REALSXP) return false;
  Shield<SEXP> props2(Rf_duplicate(proposals));
  double* pr = REAL(props2);
  int n_moved = 0;
  std::vector<double> moved(n);
  for (int s = 0; s < n; ++s) {
    const R_xlen_t row = static_cast<R_xlen_t>(s) * M + selected[s] - 1;
    for (int k = 0; k < p; ++k) pr[k + static_cast<R_xlen_t>(s) * (p + 1)] = pool_stack(row, k);
    pr[p + static_cast<R_xlen_t>(s) * (p + 1)] = ll_stack[row];
    moved[s] = selected[s] != k0[s] ? 1.0 : 0.0;
    n_moved += selected[s] != k0[s];
  }
  Shield<SEXP> subj_mu(Rf_allocMatrix(REALSXP, p, n));
  for (int s = 0; s < n; ++s)
    std::copy(REAL(tmu), REAL(tmu) + p, REAL(subj_mu) + static_cast<R_xlen_t>(s) * p);
  Shield<SEXP> alpha(proposal_rows(props2, p));
  List pars2;
  if (lat) {
    try {
      pars2 = latent_pars(List(res["state"]), kind, par_names, lget(latent, "delta_form"),
                          R_NilValue, alpha);
    } catch (LatentDecline&) {
      return false;
    }
  } else if (dg) {
    pars2 = List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["tvinv"] = tvinv_keep,
                         _["alpha"] = alpha);
  } else {
    pars2 = List::create(_["tmu"] = tmu, _["tvar"] = tvar, _["tvinv"] = res["tvinv"],
                         _["a_half"] = res["a_half"], _["subj_mu"] = subj_mu,
                         _["alpha"] = alpha);
  }

  // The move's running statistics.
  RObject stats = lget(pmwgs, "ensemble_stats");
  if (static_cast<SEXP>(stats) == R_NilValue || lget(stats, "subj_moved") == R_NilValue)
    stats = List::create(_["sweeps"] = 0, _["moved"] = 0.0, _["wess"] = 0.0,
                         _["subj_moved"] = NumericVector(n), _["subj_wess"] = NumericVector(n));
  SEXP sm = lget(stats, "subj_moved"), sw = lget(stats, "subj_wess");
  if (TYPEOF(lget(stats, "sweeps")) != INTSXP || TYPEOF(lget(stats, "moved")) != REALSXP ||
      TYPEOF(lget(stats, "wess")) != REALSXP || TYPEOF(sm) != REALSXP || XLENGTH(sm) != n ||
      TYPEOF(sw) != REALSXP || XLENGTH(sw) != n)
    return false;
  NumericVector sm2(n), sw2(n);
  for (int s = 0; s < n; ++s) {
    sm2[s] = REAL(sm)[s] + moved[s];
    volatile double share = wess[s] / M;
    sw2[s] = REAL(sw)[s] + share;
  }
  const double mean_moved = static_cast<double>(static_cast<long double>(n_moved) / n);
  volatile double mean_wess = r_mean(wess.begin(), n) / M;
  RObject stats2 = lset(stats, "sweeps", Rf_ScalarInteger(INTEGER(lget(stats, "sweeps"))[0] + 1));
  stats2 = lset(stats2, "moved", Rf_ScalarReal(REAL(lget(stats, "moved"))[0] + mean_moved));
  stats2 = lset(stats2, "wess", Rf_ScalarReal(REAL(lget(stats, "wess"))[0] + mean_wess));
  stats2 = lset(stats2, "subj_moved", sm2);
  stats2 = lset(stats2, "subj_wess", sw2);
  RObject sampler = lset(pmwgs, "ensemble_stats", stats2);
  {
    Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
    sampler = lset(sampler, "rng", rng2);
  }

  // The group state's jump from the stored one (.emc_alloc_record_group).
  if (alloc_on) {
    {
      const R_xlen_t at = static_cast<R_xlen_t>(Rf_asReal(idx)) - 1;
      const int q = 2 * p;
      std::vector<double> d(q), now(p), was(p);
      log_sd(REAL(tvar), p, now.data());
      log_sd(REAL(theta_var) + at * static_cast<R_xlen_t>(p) * p, p, was.data());
      bool finite = true;
      for (int k = 0; k < p; ++k) {
        d[k] = REAL(tmu)[k] - REAL(theta_mu)[k + at * p];
        d[p + k] = now[k] - was[k];
      }
      // Less the jump the shear and collapsed moves made within that span
      // (.emc_alloc_record_group): theirs, not the rest's.
      SEXP from = var(st, "shear_from");
      SEXP pmu = lget(pars, "tmu"), pvar = lget(pars, "tvar");
      if (TYPEOF(from) == REALSXP && XLENGTH(from) == q && TYPEOF(pmu) == REALSXP &&
          XLENGTH(pmu) == p && TYPEOF(pvar) == REALSXP &&
          XLENGTH(pvar) == static_cast<R_xlen_t>(p) * p) {
        std::vector<double> pre(p);
        log_sd(REAL(pvar), p, pre.data());
        for (int k = 0; k < p; ++k) {
          d[k] = d[k] - (REAL(pmu)[k] - REAL(from)[k]);
          d[p + k] = d[p + k] - (pre[k] - REAL(from)[p + k]);
        }
      }
      for (int k = 0; k < q; ++k) finite = finite && R_FINITE(d[k]);
      RObject w2 = w;
      if (finite) {
        Shield<SEXP> G2(Rf_duplicate(G));
        double* g = REAL(G2);
        for (int c = 0; c < q; ++c)
          for (int r = 0; r < q; ++r) {
            volatile double term = d[r] * d[c];
            g[r + static_cast<R_xlen_t>(c) * q] = g[r + static_cast<R_xlen_t>(c) * q] + term;
          }
        w2 = lset(w2, "G", G2);
        w2 = lset(w2, "nG", Rf_ScalarInteger(INTEGER(nG)[0] + 1));
      }
      Shield<SEXP> alloc2(lset(alloc, "window", w2));
      sampler = lset(sampler, "alloc", alloc2);
    }
  }
  set_var(st, "pmwgs", sampler);
  set_var(st, "pars", pars2);
  set_var(st, "proposals", props2);
  return true;
}

}  // namespace emc_iter
