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
#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>
#include <vector>

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

}  // namespace

// `v` has no mean part, so every candidate keeps mu and its mean-prior
// density lp_mu. Returns fallback = TRUE, before any random number is drawn,
// when a subject's likelihood is not served by the stored-argument path; the
// caller then runs the round in R. Otherwise: E, the selected candidate
// (0 = stay), whether any candidate was valid, the seconds spent in the
// likelihoods and, for a move, the selected state, subjects and likelihoods.
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
    lp[k] = lp_mu + lp_v;
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
  std::vector<double> lw(K, R_NegInf);
  NumericMatrix ll_new(n, n_ok > 0 ? n_ok : 1);
  double ll_seconds = 0.0;
  if (n_ok > 0) {
    List moved(n_ok);
    for (int c = 0; c < n_ok; ++c) moved[c] = as<List>(cands[ok[c]])["alpha"];
    const List rows = emc_candidate_rows_cpp(moved, par_names);
    const auto ll_t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < n; ++s) {
      RObject ll;
      try {
        ll = emc2_ll_cached(rows[s], data[s], model, R_NilValue);
      } catch (...) {
        return List::create(Named("fallback") = true);
      }
      if (TYPEOF(ll) != REALSXP || XLENGTH(ll) != n_ok || ANY_ATTRIB(ll))
        return List::create(Named("fallback") = true);
      const double* x = REAL(ll);
      for (int c = 0; c < n_ok; ++c) {
        if (ISNAN(x[c]) || x[c] == R_PosInf) return List::create(Named("fallback") = true);
        ll_new(s, c) = x[c];
      }
    }
    ll_seconds = seconds_since(ll_t0);
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
  List out = List::create(Named("fallback") = false,
                          Named("E") = static_cast<double>(e_acc),
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
    for (int i = 0; i < p; ++i) mu1[i] = mu[i] + t[j - 1] * v[i];
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
