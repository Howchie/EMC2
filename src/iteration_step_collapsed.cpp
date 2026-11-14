// The iteration's collapsed step (R/iteration_steps.R .emc_step_collapsed,
// R/collapsed_move.R .emc_collapsed_move and .emc_collapsed_round): the draw
// of the move count, each move's inner trajectories on the surrogate, the
// subjects' transport between the group states' frames, their likelihoods on
// the stored-argument path, the exact correction and the block's records.
// Products, factorisations and sums follow R's, so the step is the R step on
// the same random stream.
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
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>
#include "iteration_state.h"

using namespace Rcpp;
using namespace emc_iter;

Rcpp::List emc_shear_state_cpp(Rcpp::NumericVector mu, Rcpp::NumericMatrix V,
                               const arma::cube& prec, const arma::mat& lin, Rcpp::List spec);
Rcpp::List emc_shear_transport_cpp(const arma::mat& alpha,
                                   const arma::mat& m0, const arma::cube& R0,
                                   const arma::vec& logdet0,
                                   const arma::mat& m1, const arma::cube& R1,
                                   const arma::vec& logdet1,
                                   const Rcpp::LogicalVector& has,
                                   const arma::mat& U, const arma::vec& b,
                                   const arma::cube& Uperp);
Rcpp::List collapsed_trajectories_cpp(const arma::vec& phi0, const arma::mat& C,
                                      double eps, int n_traj, double len,
                                      const Rcpp::List& model);
Rcpp::NumericVector emc_collapsed_log_lhat_cpp(Rcpp::NumericMatrix alpha,
                                               Rcpp::NumericVector H,
                                               Rcpp::NumericMatrix lin,
                                               Rcpp::IntegerVector has);
List emc_candidate_rows_cpp(List cand, CharacterVector par_names);
SEXP emc2_ll_cached(SEXP proposals, SEXP dadm, SEXP model, SEXP varying);

namespace {

const int kInner = 3;
const double kLength = 1.5;

void zero_costs(SEXP st) {
  set_var(st, "coll_cpu", Rf_ScalarReal(0.0));
  set_var(st, "coll_el", Rf_ScalarReal(0.0));
  set_var(st, "coll_round_cpu", Rf_ScalarReal(0.0));
  set_var(st, "coll_moves", Rf_ScalarInteger(0));
}

bool all_finite(const double* x, R_xlen_t n) {
  for (R_xlen_t i = 0; i < n; ++i) if (!R_FINITE(x[i])) return false;
  return true;
}

// One block of the chart: the parameters it spans, and where its factor's
// lower triangle (column-major positions lt) sits in phi.
struct ChartBlock {
  std::vector<int> idx, lt, at;
  int d;
};
struct Chart {
  int p, dim;
  std::vector<ChartBlock> blocks;
  std::vector<int> gamma, gamma_at;
};

std::vector<int> ints(SEXP x) {
  Shield<SEXP> v(Rf_coerceVector(x, INTSXP));
  return std::vector<int>(INTEGER(v), INTEGER(v) + XLENGTH(v));
}

bool read_chart(SEXP chart, int p, Chart& out) {
  SEXP blocks = lget(chart, "blocks");
  if (TYPEOF(chart) != VECSXP || TYPEOF(blocks) != VECSXP || Rf_length(lget(chart, "dim")) != 1)
    return false;
  out.p = p;
  out.dim = Rf_asInteger(lget(chart, "dim"));
  for (R_xlen_t k = 0; k < XLENGTH(blocks); ++k) {
    SEXP b = VECTOR_ELT(blocks, k);
    ChartBlock cb;
    cb.idx = ints(lget(b, "idx"));
    cb.lt = ints(lget(b, "lt"));
    cb.at = ints(lget(b, "at"));
    cb.d = static_cast<int>(cb.idx.size());
    if (cb.lt.size() != cb.at.size()) return false;
    for (int v : cb.idx) if (v < 1 || v > p) return false;
    for (int v : cb.lt) if (v < 1 || v > cb.d * cb.d) return false;
    for (int v : cb.at) if (v < 1 || v > out.dim) return false;
    out.blocks.push_back(cb);
  }
  out.gamma = ints(lget(chart, "gamma"));
  out.gamma_at = ints(lget(chart, "gamma_at"));
  if (out.gamma.size() != out.gamma_at.size()) return false;
  for (int v : out.gamma) if (v < 1 || v > p) return false;
  for (int v : out.gamma_at) if (v < 1 || v > out.dim) return false;
  return true;
}

// .emc_collapsed_to_phi(). False when a block has no Cholesky factor.
bool to_phi(const double* mu, const double* V, const Chart& c, arma::vec& phi) {
  const int p = c.p;
  phi.zeros(c.dim);
  for (int k = 0; k < p; ++k) phi[k] = mu[k];
  for (const ChartBlock& b : c.blocks) {
    int d = b.d;
    std::vector<double> U(static_cast<size_t>(d) * d, 0.0), L(static_cast<size_t>(d) * d, 0.0);
    for (int cc = 0; cc < d; ++cc)
      for (int r = 0; r <= cc; ++r)
        U[r + static_cast<size_t>(cc) * d] =
          V[(b.idx[r] - 1) + static_cast<R_xlen_t>(b.idx[cc] - 1) * p];
    int info = 0;
    F77_CALL(dpotrf)("U", &d, U.data(), &d, &info FCONE);
    if (info != 0) return false;
    for (int cc = 0; cc < d; ++cc)
      for (int r = 0; r <= cc; ++r)
        L[cc + static_cast<size_t>(r) * d] = U[r + static_cast<size_t>(cc) * d];
    for (int k = 0; k < d; ++k) L[k + static_cast<size_t>(k) * d] = std::log(L[k + static_cast<size_t>(k) * d]);
    for (size_t j = 0; j < b.lt.size(); ++j) phi[b.at[j] - 1] = L[b.lt[j] - 1];
  }
  for (size_t j = 0; j < c.gamma.size(); ++j) {
    const int g = c.gamma[j] - 1;
    phi[c.gamma_at[j] - 1] = 0.5 * std::log(V[g + static_cast<R_xlen_t>(g) * p]);
  }
  return true;
}

// mu and V of .emc_collapsed_from_phi(). False where R's product would leave
// the BLAS path.
bool from_phi(const arma::vec& phi, const Chart& c, NumericVector& mu, NumericMatrix& V) {
  const int p = c.p;
  mu = NumericVector(p);
  V = NumericMatrix(p, p);
  for (int k = 0; k < p; ++k) mu[k] = phi[k];
  for (const ChartBlock& b : c.blocks) {
    int d = b.d;
    std::vector<double> L(static_cast<size_t>(d) * d, 0.0), S(static_cast<size_t>(d) * d, 0.0);
    for (size_t j = 0; j < b.lt.size(); ++j) L[b.lt[j] - 1] = phi[b.at[j] - 1];
    for (int k = 0; k < d; ++k) L[k + static_cast<size_t>(k) * d] = std::exp(L[k + static_cast<size_t>(k) * d]);
    if (!all_finite(L.data(), static_cast<R_xlen_t>(L.size()))) return false;
    const double one = 1.0, zero = 0.0;
    F77_CALL(dsyrk)("U", "N", &d, &d, &one, L.data(), &d, &zero, S.data(), &d FCONE FCONE);
    for (int cc = 0; cc < d; ++cc)
      for (int r = 0; r <= cc; ++r) {
        const double v = S[r + static_cast<size_t>(cc) * d];
        V(b.idx[r] - 1, b.idx[cc] - 1) = v;
        V(b.idx[cc] - 1, b.idx[r] - 1) = v;
      }
  }
  for (size_t j = 0; j < c.gamma.size(); ++j) {
    const int g = c.gamma[j] - 1;
    volatile double twice = 2 * phi[c.gamma_at[j] - 1];
    V(g, g) = std::exp(twice);
  }
  return true;
}

// .emc_alloc_lsd().
void lsd(SEXP V, std::vector<double>& out) {
  const int p = Rf_nrows(V);
  out.resize(p);
  for (int k = 0; k < p; ++k) {
    double d = REAL(V)[k + static_cast<R_xlen_t>(k) * p];
    if (d < DBL_MIN) d = DBL_MIN;
    out[k] = 0.5 * std::log(d);
  }
}

// c(mu1 - mu0, lsd(V1) - lsd(V0)).
void state_jump(SEXP mu1, SEXP V1, SEXP mu0, SEXP V0, std::vector<double>& out) {
  const int p = static_cast<int>(XLENGTH(mu1));
  std::vector<double> l1, l0;
  lsd(V1, l1);
  lsd(V0, l0);
  out.resize(2 * static_cast<size_t>(p));
  for (int k = 0; k < p; ++k) {
    out[k] = REAL(mu1)[k] - REAL(mu0)[k];
    out[p + k] = l1[k] - l0[k];
  }
}

// M + w * tcrossprod(x) in place.
void add_outer(double* M, double w, const std::vector<double>& x) {
  const size_t m = x.size();
  for (size_t c = 0; c < m; ++c)
    for (size_t r = 0; r < m; ++r) {
      volatile double prod = x[r] * x[c];
      volatile double term = w * prod;
      M[r + c * m] = M[r + c * m] + term;
    }
}

// A group state with the subjects' frames (.emc_shear_state()).
struct GroupState {
  RObject mu, V, Vi, m, R, logdet;
};

}  // namespace

namespace emc_iter {

bool log_prior_mu(SEXP mu, SEXP prior, double& out) {
  SEXP mean = lget(prior, "theta_mu_mean"), invar = lget(prior, "theta_mu_invar");
  int p = static_cast<int>(XLENGTH(mu));
  if (TYPEOF(mu) != REALSXP || TYPEOF(mean) != REALSXP || XLENGTH(mean) != p ||
      TYPEOF(invar) != REALSXP || !Rf_isMatrix(invar) || Rf_nrows(invar) != p ||
      Rf_ncols(invar) != p)
    return false;
  std::vector<double> d(p), w(p);
  for (int k = 0; k < p; ++k) d[k] = REAL(mu)[k] - REAL(mean)[k];
  if (!all_finite(d.data(), p) || !all_finite(REAL(invar), XLENGTH(invar))) return false;
  const double one = 1.0, zero = 0.0;
  int ione = 1;
  F77_CALL(dgemv)("N", &p, &p, &one, REAL(invar), &p, d.data(), &ione, &zero, w.data(),
                  &ione FCONE);
  if (!all_finite(w.data(), p)) return false;
  double z = 0.0;
  F77_CALL(dgemv)("T", &p, &ione, &one, d.data(), &p, w.data(), &ione, &zero, &z,
                  &ione FCONE);
  out = -0.5 * z;
  return true;
}

bool step_collapsed(SEXP st, SEXP cx, int i) {
  SEXP active = lget(cx, "shear_active");
  if (TYPEOF(active) != LGLSXP) return false;
  const std::string stage = as<std::string>(lget(cx, "stage"));
  if (!flag(active) || stage == "preburn") {
    zero_costs(st);
    return true;
  }
  SEXP pmwgs = var(st, "pmwgs");
  SEXP cst = lget(pmwgs, "collapsed"), plan = lget(cst, "plan"), rounds = lget(plan, "rounds");
  if (TYPEOF(plan) != VECSXP || !Rf_isMatrix(lget(plan, "C")) || Rf_length(rounds) != 1 ||
      !(Rf_asReal(rounds) > 0)) {
    if (TYPEOF(plan) == VECSXP && Rf_length(rounds) > 1) return false;
    zero_costs(st);
    return true;
  }
  SEXP rng = lget(pmwgs, "rng"), seed = lget(rng, "gibbs"), iterations = lget(cst, "iterations");
  SEXP pars = var(st, "pars"), proposals = var(st, "proposals"), wpool = var(st, "wpool");
  SEXP data = lget(cx, "data"), model = lget(pmwgs, "model"), prior = lget(pmwgs, "prior");
  SEXP par_names = lget(pmwgs, "par_names"), tmu = lget(pars, "tmu"), tvar = lget(pars, "tvar");
  SEXP blk = lget(plan, "block"), C = lget(plan, "C");
  const int p = Rf_asInteger(lget(pmwgs, "n_pars")), n = Rf_asInteger(lget(pmwgs, "n_subjects"));
  const int m = 2 * p;
  SEXP bE = lget(blk, "E"), bn = lget(blk, "n"), bacc = lget(blk, "acc"),
       bmoved = lget(blk, "moved"), bN = lget(blk, "N"), bK = lget(blk, "K");
  if (TYPEOF(seed) != INTSXP || TYPEOF(iterations) != INTSXP || XLENGTH(iterations) != 1 ||
      TYPEOF(lget(cst, "rounds")) != INTSXP || TYPEOF(lget(cst, "accepted")) != INTSXP ||
      flag(lget(wpool, "alive")) || !flag(lget(pmwgs, "ll_cache_validated")) ||
      lget(pmwgs, "marginalise") != R_NilValue || TYPEOF(data) != VECSXP || XLENGTH(data) != n ||
      Rf_getAttrib(VECTOR_ELT(data, 0), Rf_install("pooled_values")) != R_NilValue ||
      TYPEOF(par_names) != STRSXP || TYPEOF(C) != REALSXP ||
      Rf_length(lget(plan, "eps")) != 1 || Rf_length(lget(plan, "dim")) != 1 ||
      TYPEOF(blk) != VECSXP || TYPEOF(bE) != REALSXP || !Rf_isMatrix(bE) ||
      Rf_nrows(bE) != m || Rf_ncols(bE) != m || TYPEOF(bN) != REALSXP || !Rf_isMatrix(bN) ||
      Rf_nrows(bN) != m || Rf_ncols(bN) != m || TYPEOF(bn) != INTSXP || XLENGTH(bn) != 1 ||
      TYPEOF(bacc) != REALSXP || XLENGTH(bacc) != 1 || TYPEOF(bmoved) != INTSXP ||
      XLENGTH(bmoved) != 1 || TYPEOF(bK) != INTSXP || !Rf_isMatrix(bK) || Rf_nrows(bK) != 1 ||
      Rf_getAttrib(bK, R_DimNamesSymbol) != R_NilValue || lget(blk, "iters") == R_NilValue ||
      TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) ||
      Rf_nrows(proposals) != p + 1 || Rf_ncols(proposals) != n ||
      TYPEOF(tvar) != REALSXP || !Rf_isMatrix(tvar) || Rf_nrows(tvar) != p ||
      Rf_ncols(tvar) != p || TYPEOF(tmu) != REALSXP || XLENGTH(tmu) != p)
    return false;
  const Clock started = clock_now();
  static Function collapsed_context = package_function(".emc_collapsed_context");

  set_global_seed(seed);
  GetRNGstate();
  bool in_rng = true;
  auto leave_rng = [&]() { if (in_rng) { PutRNGstate(); in_rng = false; } };
  auto enter_rng = [&]() { if (!in_rng) { GetRNGstate(); in_rng = true; } };
  auto decline = [&]() { leave_rng(); return false; };

  const double r = Rf_asReal(rounds), whole = std::floor(r);
  const int n_rounds = static_cast<int>(whole) + (unif_rand() < r - whole ? 1 : 0);
  RObject cst2 = lset(cst, "iterations", Rf_ScalarInteger(INTEGER(iterations)[0] + 1));
  RObject pars2 = pars, props2 = proposals;
  double round_el = 0.0;
  int moves = 0;

  if (n_rounds > 0) {
    leave_rng();
    Shield<SEXP> ctx(collapsed_context(pmwgs));
    SEXP args = lget(ctx, "args"), spec = lget(ctx, "spec"), cpp = lget(ctx, "cpp");
    if (TYPEOF(args) != VECSXP || TYPEOF(spec) != VECSXP || TYPEOF(cpp) != VECSXP ||
        TYPEOF(lget(cpp, "has")) != INTSXP)
      return false;
    Chart chart;
    if (!read_chart(lget(ctx, "chart"), p, chart)) return false;
    if (flag(lget(ctx, "cacheable"))) cst2 = lset(cst2, "fits", lget(ctx, "base_fits"));
    const arma::cube cube = as<arma::cube>(lget(args, "cube"));
    const arma::mat lin = as<arma::mat>(lget(args, "lin"));
    const arma::mat U = as<arma::mat>(lget(args, "U"));
    const arma::vec wall_b = as<arma::vec>(lget(args, "b"));
    const arma::cube Uc = as<arma::cube>(lget(args, "Uc"));
    const LogicalVector has(lget(args, "has"));
    const bool walls = lget(ctx, "walls") != R_NilValue;
    if (walls && has.size() != n) return false;
    const arma::mat Cm = as<arma::mat>(C);
    const double eps = Rf_asReal(lget(plan, "eps"));
    const List cpp_model(cpp);
    const NumericVector lhat_H(lget(cpp, "H"));
    const NumericMatrix lhat_lin(lget(cpp, "lin"));
    const IntegerVector lhat_has(lget(cpp, "has"));
    const List spec_list(spec);

    // The state the moves start from.
    GroupState cur;
    bool have_cur = false;
    if (Rf_asReal(lget(plan, "dim")) == chart.dim) {
      NumericVector mu(p);
      std::copy(REAL(tmu), REAL(tmu) + p, mu.begin());
      const List out = emc_shear_state_cpp(mu, NumericMatrix(tvar), cube, lin, spec_list);
      double lp_mu = 0.0;
      if (flag(out["ok"])) {
        if (!log_prior_mu(mu, prior, lp_mu)) return false;
        if (R_FINITE(lp_mu + Rf_asReal(out["lp_v"]))) {
          cur.mu = mu; cur.V = tvar; cur.Vi = out["Vi"];
          cur.m = out["m"]; cur.R = out["R"]; cur.logdet = out["logdet"];
          have_cur = true;
        }
      }
    }
    if (have_cur) {
      const GroupState start = cur;
      NumericMatrix alpha(p, n);
      NumericVector ll(n);
      {
        const double* pr = REAL(proposals);
        for (int s = 0; s < n; ++s) {
          for (int k = 0; k < p; ++k) alpha(k, s) = pr[k + static_cast<R_xlen_t>(s) * (p + 1)];
          ll[s] = pr[p + static_cast<R_xlen_t>(s) * (p + 1)];
        }
      }
      NumericMatrix E = clone(NumericMatrix(bE)), N = clone(NumericMatrix(bN));
      int cnt = INTEGER(bn)[0], mv = INTEGER(bmoved)[0], evaluated = 0;
      double acc_sum = REAL(bacc)[0];
      bool any_moved = false;
      std::vector<double> jump(m);
      const Clock rounds_started = clock_now();
      enter_rng();
      for (int k = 0; k < n_rounds; ++k) {
        double acc = 0.0;
        bool moved = false, bail = false;
        std::fill(jump.begin(), jump.end(), 0.0);
        // One move; leaves early, unmoved, where the R round returns its state.
        [&]() {
          arma::vec phi0;
          if (!to_phi(REAL(cur.mu), REAL(cur.V), chart, phi0)) { bail = true; return; }
          const List inner = collapsed_trajectories_cpp(phi0, Cm, eps, kInner, kLength, cpp_model);
          if (!R_FINITE(as<double>(inner["lp"]))) return;
          const arma::vec phi = as<arma::vec>(inner["phi"]);
          bool same = true;
          for (arma::uword j = 0; j < phi.n_elem; ++j)
            if (!(phi[j] == phi0[j]) && !(ISNAN(phi[j]) && ISNAN(phi0[j]))) { same = false; break; }
          if (same) return;
          NumericVector mu1;
          NumericMatrix V1;
          if (!from_phi(phi, chart, mu1, V1)) { bail = true; return; }
          const List out = emc_shear_state_cpp(mu1, V1, cube, lin, spec_list);
          if (!flag(out["ok"])) return;
          double lp_mu = 0.0;
          if (!log_prior_mu(mu1, prior, lp_mu)) { bail = true; return; }
          if (!R_FINITE(lp_mu + Rf_asReal(out["lp_v"]))) return;
          GroupState nw;
          nw.mu = mu1; nw.V = V1; nw.Vi = out["Vi"];
          nw.m = out["m"]; nw.R = out["R"]; nw.logdet = out["logdet"];
          const arma::mat a0(alpha.begin(), p, n, false, true);
          const arma::mat m0 = as<arma::mat>(cur.m);
          const arma::cube R0 = as<arma::cube>(cur.R);
          const List tr = emc_shear_transport_cpp(
              a0, m0, R0, as<arma::vec>(cur.logdet), as<arma::mat>(nw.m),
              as<arma::cube>(nw.R), as<arma::vec>(nw.logdet), has, U, wall_b, Uc);
          const LogicalVector ok(static_cast<SEXP>(tr["ok"]));
          for (R_xlen_t s = 0; s < ok.size(); ++s) if (ok[s] != TRUE) return;
          const NumericMatrix alpha_new(static_cast<SEXP>(tr["alpha"]));
          if (!all_finite(alpha_new.begin(), XLENGTH(alpha_new))) return;
          std::vector<double> jac(n, 0.0);
          if (walls) {
            const NumericVector tj(static_cast<SEXP>(tr["jac"])), r1sq(static_cast<SEXP>(tr["r1sq"]));
            std::vector<double> dlt(p), rr(p);
            const double one = 1.0, zero = 0.0;
            int ione = 1, pp = p;
            for (int s = 0; s < n; ++s) {
              if (has[s] != TRUE) continue;
              for (int q = 0; q < p; ++q) dlt[q] = alpha(q, s) - m0(q, s);
              const double* Rs = R0.slice_memptr(s);
              if (!all_finite(dlt.data(), p) || !all_finite(Rs, static_cast<R_xlen_t>(p) * p))
                { bail = true; return; }
              F77_CALL(dgemv)("N", &pp, &pp, &one, Rs, &pp, dlt.data(), &ione, &zero, rr.data(),
                              &ione FCONE);
              long double cs = 0.0L;
              for (int q = 0; q < p; ++q) {
                volatile double sq = rr[q] * rr[q];
                cs += sq;
              }
              volatile double diff = r1sq[s] - static_cast<double>(cs);
              volatile double half = 0.5 * diff;
              jac[s] = tj[s] - half;
            }
          }
          const List rows = emc_candidate_rows_cpp(List::create(alpha_new), CharacterVector(par_names));
          NumericVector ll_new(n);
          for (int s = 0; s < n; ++s) {
            RObject v;
            try {
              v = emc2_ll_cached(rows[s], VECTOR_ELT(data, s), model, R_NilValue);
            } catch (...) {
              bail = true;
              return;
            }
            if (TYPEOF(v) != REALSXP || XLENGTH(v) != 1 || ANY_ATTRIB(v) || ISNAN(REAL(v)[0]) ||
                REAL(v)[0] == R_PosInf)
              { bail = true; return; }
            ll_new[s] = REAL(v)[0];
          }
          ++evaluated;
          const NumericVector lhat_new = emc_collapsed_log_lhat_cpp(alpha_new, lhat_H, lhat_lin, lhat_has);
          const NumericVector lhat_old = emc_collapsed_log_lhat_cpp(alpha, lhat_H, lhat_lin, lhat_has);
          long double s_new = 0.0L, s_old = 0.0L, s_jac = 0.0L;
          for (int s = 0; s < n; ++s) {
            volatile double dn = ll_new[s] - lhat_new[s];
            volatile double dold = ll[s] - lhat_old[s];
            s_new += dn;
            s_old += dold;
            s_jac += jac[s];
          }
          volatile double la = static_cast<double>(s_new) - static_cast<double>(s_old);
          la = la + static_cast<double>(s_jac);
          acc = (R_FINITE(la) && all_finite(ll_new.begin(), n)) ? std::min(1.0, std::exp(la)) : 0.0;
          state_jump(nw.mu, nw.V, cur.mu, cur.V, jump);
          if (unif_rand() < acc) {
            cur = nw;
            alpha = alpha_new;
            ll = ll_new;
            moved = true;
          }
        }();
        if (bail) return decline();
        add_outer(E.begin(), acc, jump);
        ++cnt;
        acc_sum = acc_sum + acc;
        mv += moved ? 1 : 0;
        any_moved = any_moved || moved;
      }
      round_el = clock_now().elapsed - rounds_started.elapsed;
      moves = n_rounds;
      std::vector<double> d;
      state_jump(cur.mu, cur.V, start.mu, start.V, d);
      RObject K = bK;
      if (all_finite(d.data(), m)) {
        add_outer(N.begin(), 1.0, d);
        const int cols = std::max(Rf_ncols(bK), n_rounds);
        IntegerMatrix K2(1, cols);
        std::copy(INTEGER(bK), INTEGER(bK) + XLENGTH(bK), K2.begin());
        K2(0, n_rounds - 1) += 1;
        K = K2;
      }
      RObject blk2 = lset(blk, "E", E);
      blk2 = lset(blk2, "n", Rf_ScalarInteger(cnt));
      blk2 = lset(blk2, "acc", Rf_ScalarReal(acc_sum));
      blk2 = lset(blk2, "moved", Rf_ScalarInteger(mv));
      blk2 = lset(blk2, "N", N);
      blk2 = lset(blk2, "K", K);
      RObject plan2 = lset(plan, "block", blk2);
      cst2 = lset(cst2, "plan", plan2);
      cst2 = lset(cst2, "rounds",
                  Rf_ScalarInteger(INTEGER(lget(cst, "rounds"))[0] + evaluated));
      if (any_moved) {
        cst2 = lset(cst2, "accepted",
                    Rf_ScalarInteger(INTEGER(lget(cst, "accepted"))[0] + 1));
        if (flag(lget(spec, "uses_a_half"))) {
          std::vector<double> ah;
          if (!draw_a_half(cur.V, spec, ah)) return decline();
          NumericVector a_half(ah.begin(), ah.end());
          SEXP nm = Rf_getAttrib(lget(pars, "a_half"), R_NamesSymbol);
          if (nm != R_NilValue) Rf_setAttrib(a_half, R_NamesSymbol, nm);
          pars2 = lset(pars2, "a_half", a_half);
        }
        Shield<SEXP> mu2(Rf_duplicate(cur.mu));
        Rf_setAttrib(mu2, R_NamesSymbol, Rf_getAttrib(tmu, R_NamesSymbol));
        Shield<SEXP> V2(Rf_duplicate(cur.V)), Vi2(Rf_duplicate(cur.Vi));
        Rf_setAttrib(V2, R_DimNamesSymbol, Rf_getAttrib(tvar, R_DimNamesSymbol));
        Rf_setAttrib(Vi2, R_DimNamesSymbol, Rf_getAttrib(tvar, R_DimNamesSymbol));
        pars2 = lset(pars2, "tmu", mu2);
        pars2 = lset(pars2, "tvar", V2);
        pars2 = lset(pars2, "tvinv", Vi2);
        SEXP subj_mu = lget(pars, "subj_mu");
        if (Rf_isMatrix(subj_mu)) {
          if (TYPEOF(subj_mu) != REALSXP) return decline();
          Shield<SEXP> sm(Rf_duplicate(subj_mu));
          for (R_xlen_t q = 0; q < XLENGTH(sm); ++q) REAL(sm)[q] = REAL(cur.mu)[q % p];
          pars2 = lset(pars2, "subj_mu", sm);
        }
        props2 = Rf_duplicate(proposals);
        double* pr = REAL(props2);
        for (int s = 0; s < n; ++s) {
          for (int q = 0; q < p; ++q) pr[q + static_cast<R_xlen_t>(s) * (p + 1)] = alpha(q, s);
          pr[p + static_cast<R_xlen_t>(s) * (p + 1)] = ll[s];
        }
        pars2 = lset(pars2, "alpha", alpha);
      }
    }
  }
  leave_rng();
  Shield<SEXP> seed_after(global_seed());
  Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
  RObject sampler = lset(pmwgs, "rng", rng2);
  sampler = lset(sampler, "collapsed", cst2);
  set_var(st, "pmwgs", sampler);
  set_var(st, "pars", pars2);
  set_var(st, "proposals", props2);
  const double el = clock_now().elapsed - started.elapsed;
  const double cores = num(var(st, "cores_it"), 1.0);
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  set_var(st, "coll_moves", Rf_ScalarInteger(moves));
  set_var(st, "coll_el", Rf_ScalarReal(el));
  set_var(st, "coll_cpu", Rf_ScalarReal(alloc_on ? reserved_cost(el, cores) : 0.0));
  set_var(st, "coll_round_cpu", Rf_ScalarReal(alloc_on ? reserved_cost(round_el, cores) : 0.0));
  return true;
}

}  // namespace emc_iter
