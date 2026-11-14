// One shear round along a log-sd direction in a single call
// (R/shear_move.R .emc_shear_round): the candidates, their subject
// likelihoods on the stored-argument path, the selection weights and the
// selection. Sums and the selection follow R's sum(), colSums() and
// sample.int(), so the round is the R round on the same random stream.
// [[Rcpp::depends(RcppArmadillo)]]
// The weights and sums are compared with R's bit for bit, so this file's
// arithmetic is strict IEEE whatever the package's flags.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <RcppArmadillo.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif
#include <chrono>
#include <cmath>
#include <cstring>
#include <vector>
#include "pool_exchange.h"
#include "subject_ll.h"

using namespace Rcpp;

Rcpp::List emc_shear_candidates_cpp(const arma::vec& t, const arma::vec& v,
                                    const arma::mat& alpha, const arma::vec& mu,
                                    const arma::mat& V,
                                    const arma::mat& m0, const arma::cube& R0,
                                    const arma::vec& logdet0,
                                    const arma::cube& prec, const arma::mat& lin,
                                    const Rcpp::LogicalVector& has,
                                    const arma::mat& U, const arma::vec& b,
                                    const arma::cube& Uperp,
                                    Rcpp::Nullable<Rcpp::List> spec);
List emc_candidate_rows_cpp(List cand, CharacterVector par_names);
SEXP emc2_ll_cached(SEXP proposals, SEXP dadm, SEXP model, SEXP varying);

namespace {

double seconds_since(const std::chrono::steady_clock::time_point& t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// sample.int(n, 1, prob = w) - 1 for weights that are finite, non-negative
// and not all zero.
int sample_one(const std::vector<double>& w) {
  const int n = static_cast<int>(w.size());
  std::vector<double> p(w);
  double sum = 0.0;
  for (int i = 0; i < n; ++i) sum += p[i];
  for (int i = 0; i < n; ++i) p[i] /= sum;
  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = i + 1;
  revsort(p.data(), perm.data(), n);
  const double rT = 1.0 * unif_rand();
  double mass = 0.0;
  int j = 0;
  for (j = 0; j < n - 1; ++j) {
    mass += p[j];
    if (rT <= mass) break;
  }
  return perm[j] - 1;
}

// Every subject's likelihood of its candidate rows on the stored-argument
// path, a column per candidate. False when a subject's likelihood is not
// served there.
bool stored_lls(const List& rows, const List& data, SEXP model, int n_ok,
                NumericMatrix& ll_new) {
  const int n = ll_new.nrow();
  for (int s = 0; s < n; ++s) {
    RObject ll;
    try {
      ll = emc2_ll_cached(rows[s], data[s], model, R_NilValue);
    } catch (...) {
      return false;
    }
    if (TYPEOF(ll) != REALSXP || XLENGTH(ll) != n_ok || ANY_ATTRIB(ll)) return false;
    const double* x = REAL(ll);
    for (int c = 0; c < n_ok; ++c) {
      if (ISNAN(x[c]) || x[c] == R_PosInf) return false;
      ll_new(s, c) = x[c];
    }
  }
  return true;
}

// The round's selection among staying (0) and its K candidates with steps t:
// the valid candidates `ok` (0-based) weigh by their likelihood columns
// against the current likelihoods plus their prior change. E is the
// selection-weighted squared step.
int select_candidate(const double* t, int K, const std::vector<int>& ok,
                     const std::vector<double>& prior_delta, const NumericMatrix& ll_new,
                     const NumericVector& ll_cur, double& E) {
  const int n = ll_new.nrow(), n_ok = static_cast<int>(ok.size());
  std::vector<double> lw(K, R_NegInf);
  if (n_ok > 0) {
    long double cur_acc = 0.0L;
    for (int s = 0; s < n; ++s) cur_acc += ll_cur[s];
    const double cur_sum = static_cast<double>(cur_acc);
    for (int c = 0; c < n_ok; ++c) {
      long double col = 0.0L;
      for (int s = 0; s < n; ++s) col += ll_new(s, c);
      const double v_ok = static_cast<double>(col) - cur_sum + prior_delta[ok[c]];
      lw[ok[c]] = R_FINITE(v_ok) ? v_ok : R_NegInf;
    }
  }
  std::vector<double> w(K + 1);
  w[0] = 0.0;
  for (int k = 0; k < K; ++k) w[k + 1] = lw[k];
  double top = w[0];
  for (int k = 1; k <= K; ++k) if (w[k] > top) top = w[k];
  // One scalar exp per weight, as R's exp(): a vectorised exp may differ in
  // the last bit.
  for (int k = 0; k <= K; ++k) {
    volatile double arg = w[k] - top;
    w[k] = std::exp(arg);
  }
  long double total = 0.0L;
  for (int k = 0; k <= K; ++k) total += w[k];
  const double total_d = static_cast<double>(total);
  for (int k = 0; k <= K; ++k) w[k] /= total_d;
  const int j = sample_one(w);
  long double e_acc = 0.0L;
  // Each term is rounded to double before it is added, as R's sum() of a
  // product vector.
  for (int k = 0; k < K; ++k) {
    volatile double term = t[k] * t[k];
    term = w[k + 1] * term;
    e_acc += term;
  }
  E = static_cast<double>(e_acc);
  return j;
}

}  // namespace

bool emc_subject_ll_in_process(SEXP pool, SEXP sampler) {
  auto get = [](SEXP x, const char* name) -> SEXP {
    if (TYPEOF(x) != VECSXP) return R_NilValue;
    SEXP names = Rf_getAttrib(x, R_NamesSymbol);
    for (R_xlen_t i = 0; i < Rf_xlength(names); ++i)
      if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(x, i);
    return R_NilValue;
  };
  auto is_true = [](SEXP x) {
    return TYPEOF(x) == LGLSXP && Rf_length(x) == 1 && LOGICAL(x)[0] == TRUE;
  };
  SEXP data = get(sampler, "data");
  return !is_true(get(pool, "alive")) && is_true(get(sampler, "ll_cache_validated")) &&
    get(sampler, "marginalise") == R_NilValue && TYPEOF(data) == VECSXP &&
    XLENGTH(data) > 0 && Rf_inherits(VECTOR_ELT(data, 0), "data.frame") &&
    Rf_getAttrib(VECTOR_ELT(data, 0), Rf_install("pooled_values")) == R_NilValue;
}

bool emc_subject_lls(const List& rows, const List& data, SEXP model, int n_cand,
                     NumericMatrix& ll, SubjectLl* src) {
  if (src == nullptr || src->in_process) {
    if (stored_lls(rows, data, model, n_cand, ll)) return true;
    if (src == nullptr) return false;
  }
  const R_xlen_t n = ll.nrow();
  // The shared subject queue's failure paths run R code, as the dispatch does.
  PutRNGstate();
  if (rows.size() == n &&
      emc_pool::subject_ll_queue(src->pool, rows, src->sampler, src->part, src->source, n_cand,
                                 ll.begin())) {
    GetRNGstate();
    return true;
  }
  Function dispatch(src->dispatch);
  const List got = dispatch(src->pool, rows, src->sampler, src->part,
                            Named("source") = std::string(src->source));
  GetRNGstate();
  src->pool = got["pool"];
  SEXP v = got["ll"];
  if (TYPEOF(v) != REALSXP || XLENGTH(v) != n * n_cand)
    Rcpp::stop("subject likelihoods: the round's values differ in shape");
  std::copy(REAL(v), REAL(v) + XLENGTH(v), ll.begin());
  return true;
}

// `v` has no mean part, so every candidate keeps mu and its mean-prior
// density lp_mu. Returns fallback = TRUE, before any random number is drawn,
// when a subject's likelihood is not served by the stored-argument path; the
// caller then runs the round in R. Otherwise: E, the selected candidate
// (0 = stay), whether any candidate was valid, the seconds spent in the
// likelihoods and, for a move, the selected state, subjects and likelihoods.
// With a likelihood source the round never falls back: the source serves
// whatever the stored-argument path does not. A `v` with a mean part takes
// each candidate's mean-prior density from `lp_mu_each`.
List shear_round_from(const arma::vec& t, const arma::vec& v,
                      const arma::mat& alpha, const arma::vec& mu,
                      const arma::mat& V, const arma::mat& m0,
                      const arma::cube& R0, const arma::vec& logdet0,
                      const arma::cube& prec, const arma::mat& lin,
                      const LogicalVector& has, const arma::mat& U,
                      const arma::vec& b, const arma::cube& Uperp,
                      List spec, double lp_mu, double cur_hier, double cur_lp,
                      NumericVector ll_cur, List data, SEXP model,
                      CharacterVector par_names, SubjectLl* src,
                      const double* lp_mu_each) {
  const int K = static_cast<int>(t.n_elem);
  const int p = static_cast<int>(alpha.n_rows), n = static_cast<int>(alpha.n_cols);
  if (static_cast<int>(v.n_elem) != 2 * p || data.size() != n || ll_cur.size() != n)
    Rcpp::stop("shear round: inputs differ in shape");
  List cands = emc_shear_candidates_cpp(t, v, alpha, mu, V, m0, R0, logdet0, prec, lin,
                                        has, U, b, Uperp, spec);
  const IntegerVector size = spec["size"];
  std::vector<int> ok;
  std::vector<double> prior_delta(K, R_NegInf), lp(K, NA_REAL);
  for (int k = 0; k < K; ++k) {
    const List o = cands[k];
    if (!as<bool>(o["ok"])) continue;
    const double lp_v = as<double>(o["lp_v"]);
    lp[k] = (lp_mu_each != nullptr ? lp_mu_each[k] : lp_mu) + lp_v;
    if (!R_FINITE(lp[k])) continue;
    long double acc = 0.0L;
    for (int i = 0; i < p; ++i) {
      volatile double term = (static_cast<double>(size[i]) + 1) * t[k];
      term = term * v[p + i];
      acc += term;
    }
    const double log_jac = static_cast<double>(acc) + as<double>(o["jac"]);
    prior_delta[k] = as<double>(o["hier"]) + lp[k] + log_jac - cur_hier - cur_lp;
    ok.push_back(k);
  }
  const int n_ok = static_cast<int>(ok.size());
  NumericMatrix ll_new(n, n_ok > 0 ? n_ok : 1);
  double ll_seconds = 0.0;
  if (n_ok > 0) {
    List moved(n_ok);
    for (int c = 0; c < n_ok; ++c) moved[c] = as<List>(cands[ok[c]])["alpha"];
    const List rows = emc_candidate_rows_cpp(moved, par_names);
    const auto ll_t0 = std::chrono::steady_clock::now();
    if (!emc_subject_lls(rows, data, model, n_ok, ll_new, src))
      return List::create(Named("fallback") = true);
    ll_seconds = seconds_since(ll_t0);
  }
  double E = 0.0;
  const int j = select_candidate(t.memptr(), K, ok, prior_delta, ll_new, ll_cur, E);
  List out = List::create(Named("fallback") = false,
                          Named("E") = E,
                          Named("j") = j, Named("evaluated") = n_ok > 0,
                          Named("ll_seconds") = ll_seconds);
  if (j > 0) {
    const List o = cands[j - 1];
    int col = 0;
    while (col < n_ok && ok[col] != j - 1) ++col;
    if (col == n_ok) Rcpp::stop("shear round: an invalid candidate was selected");
    NumericVector ll_sel(n);
    for (int s = 0; s < n; ++s) ll_sel[s] = ll_new(s, col);
    NumericVector mu1(p);
    for (int i = 0; i < p; ++i) {
      volatile double step = t[j - 1] * v[i];
      mu1[i] = mu[i] + step;
    }
    out["mu"] = mu1;
    out["V"] = o["V"];
    out["Vi"] = o["Vi"];
    out["chol_V"] = o["chol_V"];
    out["m"] = o["m"];
    out["R"] = o["R"];
    out["logdet"] = o["logdet"];
    out["lp_v"] = o["lp_v"];
    out["lp"] = lp[j - 1];
    out["alpha"] = o["alpha"];
    out["ll"] = ll_sel;
  }
  return out;
}

// [[Rcpp::export]]
List emc_shear_round_cpp(const arma::vec& t, const arma::vec& v,
                         const arma::mat& alpha, const arma::vec& mu,
                         const arma::mat& V, const arma::mat& m0,
                         const arma::cube& R0, const arma::vec& logdet0,
                         const arma::cube& prec, const arma::mat& lin,
                         const LogicalVector& has, const arma::mat& U,
                         const arma::vec& b, const arma::cube& Uperp,
                         List spec, double lp_mu, double cur_hier, double cur_lp,
                         NumericVector ll_cur, List data, SEXP model,
                         CharacterVector par_names) {
  return shear_round_from(t, v, alpha, mu, V, m0, R0, logdet0, prec, lin, has, U, b, Uperp,
                          spec, lp_mu, cur_hier, cur_lp, ll_cur, data, model, par_names,
                          nullptr, nullptr);
}

// .emc_shear_log_hier: the subjects' joint log density under the group
// state, by the same triangular solve and the same sums.
double shear_log_hier(const NumericMatrix& alpha, const NumericVector& mu,
                      const NumericMatrix& chol_V) {
  const int p = alpha.nrow(), n = alpha.ncol();
  std::vector<double> z(static_cast<size_t>(p) * n);
  for (int s = 0; s < n; ++s)
    for (int i = 0; i < p; ++i) z[i + static_cast<size_t>(s) * p] = alpha(i, s) - mu[i];
  const double one = 1.0;
  F77_CALL(dtrsm)("L", "U", "T", "N", &p, &n, &one, chol_V.begin(), &p, z.data(), &p
                  FCONE FCONE FCONE FCONE);
  long double ss = 0.0L;
  for (size_t i = 0; i < z.size(); ++i) {
    volatile double sq = z[i] * z[i];
    ss += sq;
  }
  long double sl = 0.0L;
  for (int i = 0; i < p; ++i) sl += std::log(chol_V(i, i));
  volatile double a = -0.5 * static_cast<double>(ss);
  volatile double b = n * static_cast<double>(sl);
  return a - b;
}

// One shear round along a mean direction d with no transport walls
// (R/shear_move.R .emc_shear_round with .emc_shear_line): V and every frame
// factor stay, each frame centre moves along Delta_s = Q_s^-1 V^-1 d, and the
// candidates' prior change is a quadratic in the step. Returns fallback =
// TRUE when a subject's likelihood is not served by the stored-argument path
// or a product would leave R's BLAS path. Otherwise E, the selected candidate
// (0 = stay), whether any candidate was valid and, for a move, the moved
// mean, frame centres, subjects and likelihoods.
List shear_mean_round(const NumericVector& t, const NumericVector& d,
                      const NumericMatrix& alpha, const NumericVector& mu,
                      const NumericMatrix& chol_V, const NumericMatrix& Vi,
                      const NumericMatrix& m0, const NumericVector& R0,
                      const NumericVector& prior_mean, const NumericMatrix& prior_invar,
                      const NumericVector& ll_cur, const List& data, SEXP model,
                      const CharacterVector& par_names, SubjectLl* src) {
  const int K = t.size();
  int p = alpha.nrow(), n = alpha.ncol();
  const double one = 1.0, zero = 0.0;
  int ione = 1;
  auto finite = [](const double* x, R_xlen_t len) {
    for (R_xlen_t i = 0; i < len; ++i) if (!R_FINITE(x[i])) return false;
    return true;
  };
  const List fallback = List::create(Named("fallback") = true);
  if (d.size() != p || mu.size() != p || prior_mean.size() != p || prior_invar.nrow() != p ||
      prior_invar.ncol() != p || XLENGTH(R0) != static_cast<R_xlen_t>(p) * p * n ||
      !finite(Vi.begin(), XLENGTH(Vi)) || !finite(prior_invar.begin(), XLENGTH(prior_invar)) ||
      !finite(d.begin(), p))
    return fallback;
  // Each centre's move for a mean step u: Q_s^-1 V^-1 u through the frame's
  // factor, lower then upper.
  std::vector<double> Rt(static_cast<size_t>(p) * p);
  auto centre_shift = [&](const double* u, std::vector<double>& out) {
    std::vector<double> Vi_u(p);
    F77_CALL(dgemv)("N", &p, &p, &one, Vi.begin(), &p, u, &ione, &zero, Vi_u.data(),
                    &ione FCONE);
    out.resize(static_cast<size_t>(p) * n);
    for (int s = 0; s < n; ++s) {
      const double* R = R0.begin() + static_cast<R_xlen_t>(s) * p * p;
      for (int c = 0; c < p; ++c)
        for (int r = 0; r < p; ++r) Rt[c + static_cast<size_t>(r) * p] = R[r + static_cast<size_t>(c) * p];
      double* x = out.data() + static_cast<size_t>(s) * p;
      std::copy(Vi_u.begin(), Vi_u.end(), x);
      F77_CALL(dtrsm)("L", "L", "N", "N", &p, &ione, &one, Rt.data(), &p, x, &p
                      FCONE FCONE FCONE FCONE);
      F77_CALL(dtrsm)("L", "U", "N", "N", &p, &ione, &one, R, &p, x, &p
                      FCONE FCONE FCONE FCONE);
    }
  };
  std::vector<double> Delta;
  centre_shift(d.begin(), Delta);
  std::vector<double> Z0(static_cast<size_t>(p) * n), Zd(static_cast<size_t>(p) * n);
  for (int s = 0; s < n; ++s)
    for (int i = 0; i < p; ++i) {
      const size_t at = i + static_cast<size_t>(s) * p;
      Z0[at] = alpha(i, s) - mu[i];
      Zd[at] = Delta[at] - d[i];
    }
  F77_CALL(dtrsm)("L", "U", "T", "N", &p, &n, &one, chol_V.begin(), &p, Z0.data(), &p
                  FCONE FCONE FCONE FCONE);
  F77_CALL(dtrsm)("L", "U", "T", "N", &p, &n, &one, chol_V.begin(), &p, Zd.data(), &p
                  FCONE FCONE FCONE FCONE);
  std::vector<double> P0d(p);
  F77_CALL(dgemv)("N", &p, &p, &one, prior_invar.begin(), &p, d.begin(), &ione, &zero,
                  P0d.data(), &ione FCONE);
  long double zz = 0.0L, dd = 0.0L, mp = 0.0L, dp = 0.0L;
  for (size_t i = 0; i < Z0.size(); ++i) {
    volatile double cross = Z0[i] * Zd[i];
    volatile double sq = Zd[i] * Zd[i];
    zz += cross;
    dd += sq;
  }
  for (int i = 0; i < p; ++i) {
    volatile double centred = mu[i] - prior_mean[i];
    volatile double a = centred * P0d[i];
    volatile double b = d[i] * P0d[i];
    mp += a;
    dp += b;
  }
  volatile double a1 = static_cast<double>(zz) + static_cast<double>(mp);
  volatile double a2 = static_cast<double>(dd) + static_cast<double>(dp);
  volatile double half_a2 = 0.5 * a2;
  std::vector<int> ok;
  std::vector<double> prior_delta(K);
  for (int k = 0; k < K; ++k) {
    volatile double lin_term = a1 * t[k];
    volatile double tsq = t[k] * t[k];
    volatile double quad = half_a2 * tsq;
    volatile double total = lin_term + quad;
    prior_delta[k] = -total;
    if (R_FINITE(prior_delta[k])) ok.push_back(k);
  }
  auto moved_alpha = [&](double step) {
    NumericMatrix a(p, n);
    for (R_xlen_t i = 0; i < XLENGTH(a); ++i) {
      volatile double shift = step * Delta[i];
      a[i] = alpha[i] + shift;
    }
    return a;
  };
  const int n_ok = static_cast<int>(ok.size());
  NumericMatrix ll_new(n, n_ok > 0 ? n_ok : 1);
  if (n_ok > 0) {
    List moved(n_ok);
    for (int c = 0; c < n_ok; ++c) moved[c] = moved_alpha(t[ok[c]]);
    const List rows = emc_candidate_rows_cpp(moved, par_names);
    if (!emc_subject_lls(rows, data, model, n_ok, ll_new, src)) return fallback;
  }
  double E = 0.0;
  const int j = select_candidate(t.begin(), K, ok, prior_delta, ll_new, ll_cur, E);
  List out = List::create(Named("fallback") = false, Named("E") = E, Named("j") = j,
                          Named("evaluated") = n_ok > 0);
  if (j > 0) {
    int col = 0;
    while (col < n_ok && ok[col] != j - 1) ++col;
    if (col == n_ok) Rcpp::stop("shear round: an invalid candidate was selected");
    const double step = t[j - 1];
    NumericVector u(p), mu1(p), ll_sel(n);
    for (int i = 0; i < p; ++i) {
      volatile double ui = step * d[i];
      u[i] = ui;
      mu1[i] = mu[i] + ui;
    }
    std::vector<double> shift;
    centre_shift(u.begin(), shift);
    NumericMatrix m1(p, n);
    for (R_xlen_t i = 0; i < XLENGTH(m1); ++i) m1[i] = m0[i] + shift[i];
    for (int s = 0; s < n; ++s) ll_sel[s] = ll_new(s, col);
    out["mu"] = mu1;
    out["m"] = m1;
    out["alpha"] = moved_alpha(step);
    out["ll"] = ll_sel;
  }
  return out;
}

// Consecutive shear rounds along log-sd coordinates `ks` (1-based, each
// above the number of parameters), each as emc_shear_round_cpp from the
// state the previous one left: the candidate steps are drawn here as
// .emc_shear_round draws them, and the step sizes adapt as .emc_shear_move
// adapts them. `done` rounds completed; when it is short of length(ks) the
// next round's likelihoods need the R path and `t_next` holds its steps,
// already drawn; the run also ends, without `t_next`, after a round that
// changed the group mean, whose prior density `lp_mu` no longer holds. Returns each completed round's net step, expected squared
// jump and whether it moved, the step sizes and counts, and the final state
// when any round moved.
List shear_rounds_from(IntegerVector ks, NumericVector log_step, IntegerVector blk_n,
                       bool adapt, NumericVector step_range, double move_target,
                       int n_candidates,
                       NumericMatrix alpha, NumericVector mu, NumericMatrix V,
                       NumericMatrix chol_V, NumericMatrix m0, NumericVector R0,
                       NumericVector logdet0, double cur_hier, double cur_lp,
                       double lp_mu, NumericVector ll,
                       const arma::cube& prec, const arma::mat& lin,
                       const LogicalVector& has, const arma::mat& U,
                       const arma::vec& b, const arma::cube& Uperp, List spec,
                       List data, SEXP model, CharacterVector par_names, SubjectLl* src) {
  const int n_rounds = ks.size();
  const int p = alpha.nrow();
  NumericVector steps = clone(log_step);
  IntegerVector counts = clone(blk_n);
  NumericVector t_used(n_rounds), E(n_rounds);
  LogicalVector moved(n_rounds);
  NumericVector t_next;
  int done = 0, evaluated = 0;
  double ll_seconds = 0.0;
  bool any_moved = false, mu_changed = false;
  RObject Vi, lp_v;
  NumericVector logdet_cur = logdet0;
  for (; done < n_rounds; ++done) {
    const int k = ks[done] - 1;
    if (k < p || k >= 2 * p) Rcpp::stop("shear rounds: not a log-sd coordinate");
    const double step = std::exp(steps[k]);
    NumericVector t(n_candidates);
    volatile double t_aux = step * norm_rand();
    for (int c = 0; c < n_candidates; ++c) {
      volatile double inc = step * norm_rand();
      t[c] = t_aux + inc;
    }
    arma::vec v(2 * p, arma::fill::zeros);
    v[k] = 1.0;
    const arma::ivec dims = {static_cast<arma::sword>(p), static_cast<arma::sword>(p),
                             static_cast<arma::sword>(alpha.ncol())};
    const List o = shear_round_from(
        as<arma::vec>(t), v, as<arma::mat>(alpha), as<arma::vec>(mu), as<arma::mat>(V),
        as<arma::mat>(m0), arma::cube(R0.begin(), p, p, alpha.ncol(), false, true),
        as<arma::vec>(logdet_cur), prec, lin, has, U, b, Uperp, spec, lp_mu, cur_hier,
        cur_lp, ll, data, model, par_names, src, nullptr);
    if (as<bool>(o["fallback"])) {
      t_next = t;
      break;
    }
    if (as<bool>(o["evaluated"])) ++evaluated;
    ll_seconds += as<double>(o["ll_seconds"]);
    const int j = as<int>(o["j"]);
    E[done] = as<double>(o["E"]);
    moved[done] = j > 0;
    t_used[done] = j > 0 ? t[j - 1] : 0.0;
    if (j > 0) {
      any_moved = true;
      const NumericVector mu_new = as<NumericVector>(o["mu"]);
      for (int i = 0; i < p; ++i)
        if (!(mu_new[i] == mu[i])) mu_changed = true;
      mu = mu_new;
      V = as<NumericMatrix>(o["V"]);
      Vi = o["Vi"];
      chol_V = as<NumericMatrix>(o["chol_V"]);
      m0 = as<NumericMatrix>(o["m"]);
      R0 = as<NumericVector>(o["R"]);
      logdet_cur = as<NumericVector>(o["logdet"]);
      lp_v = o["lp_v"];
      cur_lp = as<double>(o["lp"]);
      alpha = as<NumericMatrix>(o["alpha"]);
      ll = as<NumericVector>(o["ll"]);
      cur_hier = shear_log_hier(alpha, mu, chol_V);
    }
    counts[k] = counts[k] + 1;
    if (adapt) {
      volatile double gain = ((j > 0 ? 1 : 0) - move_target) / std::sqrt(static_cast<double>(counts[k]));
      double next = steps[k] + gain;
      if (next < step_range[0]) next = step_range[0];
      if (next > step_range[1]) next = step_range[1];
      steps[k] = next;
    }
    if (mu_changed) {
      ++done;
      break;
    }
  }
  List out = List::create(Named("done") = done, Named("t") = t_used, Named("E") = E,
                          Named("moved") = moved, Named("log_step") = steps,
                          Named("evaluated") = evaluated,
                          Named("ll_seconds") = ll_seconds, Named("t_next") = t_next,
                          Named("any_moved") = any_moved);
  if (any_moved) {
    out["mu"] = mu; out["V"] = V; out["Vi"] = Vi; out["chol_V"] = chol_V;
    out["m"] = m0; out["R"] = R0; out["logdet"] = logdet_cur; out["lp_v"] = lp_v;
    out["lp"] = cur_lp; out["alpha"] = alpha; out["ll"] = ll; out["cur_hier"] = cur_hier;
  }
  return out;
}

// [[Rcpp::export]]
List emc_shear_rounds_cpp(IntegerVector ks, NumericVector log_step, IntegerVector blk_n,
                          bool adapt, NumericVector step_range, double move_target,
                          int n_candidates,
                          NumericMatrix alpha, NumericVector mu, NumericMatrix V,
                          NumericMatrix chol_V, NumericMatrix m0, NumericVector R0,
                          NumericVector logdet0, double cur_hier, double cur_lp,
                          double lp_mu, NumericVector ll,
                          const arma::cube& prec, const arma::mat& lin,
                          const LogicalVector& has, const arma::mat& U,
                          const arma::vec& b, const arma::cube& Uperp, List spec,
                          List data, SEXP model, CharacterVector par_names) {
  return shear_rounds_from(ks, log_step, blk_n, adapt, step_range, move_target, n_candidates,
                           alpha, mu, V, chol_V, m0, R0, logdet0, cur_hier, cur_lp, lp_mu, ll,
                           prec, lin, has, U, b, Uperp, spec, data, model, par_names, nullptr);
}
