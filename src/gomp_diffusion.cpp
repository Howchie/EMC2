// R-facing vector helpers for the Gompertz growth-process race.  The sampled
// likelihood uses the raw adapters in model_GOM.h; this wrapper gives dGOM/pGOM
// the same cache-backed path for prediction and the R likelihood fallback.

#include "utility_functions.h"
#include "model_GOM_core.h"

using namespace Rcpp;

namespace {

fperace::BndSpec gomp_bnd_at_vec(int bkind, const NumericVector& Binf,
                                 const NumericVector& tau,
                                 const NumericVector& pw, int i) {
  fperace::BndSpec bs;
  bs.kind = bkind;
  if (bkind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = Binf[i];
  bs.tau = tau[i];
  bs.pw = (bkind == fpe::FPE_BND_WEIBULL) ? pw[i] : 0.0;
  return bs;
}

void check_gomp_vec_lengths(int n, const NumericVector& x, const char* what) {
  if (x.size() != n) stop("gomp_pdf_cdf_vec: %s must match length(rt).", what);
}

} // namespace

// [[Rcpp::export]]
Rcpp::List gomp_pdf_cdf_vec(NumericVector rt, NumericVector alpha,
                            NumericVector beta, NumericVector K,
                            NumericVector B, NumericVector A,
                            NumericVector t0, int bkind = 0,
                            NumericVector Binf = NumericVector::create(),
                            NumericVector tau = NumericVector::create(),
                            NumericVector pw = NumericVector::create(),
                            int nx = 512, double dt_target = 4e-3,
                            double grade = 8.0, double tgrade = 32.0) {
  const int n = rt.size();
  check_gomp_vec_lengths(n, alpha, "alpha");
  check_gomp_vec_lengths(n, beta, "beta");
  check_gomp_vec_lengths(n, K, "K");
  check_gomp_vec_lengths(n, B, "B");
  check_gomp_vec_lengths(n, A, "A");
  check_gomp_vec_lengths(n, t0, "t0");
  if (bkind != fpe::FPE_BND_FIXED &&
      (Binf.size() != n || tau.size() != n ||
       (bkind == fpe::FPE_BND_WEIBULL && pw.size() != n))) {
    stop("gomp_pdf_cdf_vec: collapsing boundary needs Binf/tau (and pw for Weibull) at length(rt).");
  }

  NumericVector pdf(n, 0.0), cdf(n, 0.0);
  fperace::SolveCache C;
  C.grid.nx = std::max(8, nx);
  C.grid.dt_target = std::max(1e-6, dt_target);
  C.grid.grade = grade;
  C.grid.tgrade = tgrade;
  C.bnd_kind = bkind;

  std::vector<fperace::Key> keys;
  std::vector<double> horizons;
  std::vector<int> group(n, -1);
  for (int i = 0; i < n; ++i) {
    const double tt = rt[i] - t0[i];
    if (!R_finite(tt) || !(tt > 0.0)) {
      if (R_finite(rt[i]) && R_finite(alpha[i]) && R_finite(beta[i]) &&
          R_finite(K[i]) && R_finite(B[i]) && R_finite(A[i]) &&
          R_finite(t0[i]) && alpha[i] > 0.0 && beta[i] > 0.0 &&
          K[i] > 0.0 && B[i] >= 0.0 && A[i] >= 0.0 && t0[i] >= 0.0 &&
          rt[i] == R_PosInf)
        cdf[i] = 1.0;
      continue;
    }
    fperace::Key p;
    if (!gomp_key(alpha[i], beta[i], K[i], B[i], A[i],
                  gomp_bnd_at_vec(bkind, Binf, tau, pw, i), p)) continue;
    int g = -1;
    for (size_t j = 0; j < keys.size(); ++j) {
      if (keys[j] == p) { g = static_cast<int>(j); break; }
    }
    if (g < 0) {
      keys.push_back(p);
      horizons.push_back(tt);
      g = static_cast<int>(keys.size()) - 1;
    } else if (tt > horizons[g]) {
      horizons[g] = tt;
    }
    group[i] = g;
  }

  std::vector<int> cache_idx(keys.size(), -1);
  for (size_t j = 0; j < keys.size(); ++j)
    cache_idx[j] = fperace::cache_get(C, keys[j], horizons[j]);

  for (int i = 0; i < n; ++i) {
    if (group[i] < 0) continue;
    const double tt = rt[i] - t0[i];
    const fperace::Entry& en = C.e[cache_idx[group[i]]];
    const double lp = fperace::entry_log_pdf(en, tt);
    const double ls = fperace::entry_log_S(en, tt);
    pdf[i] = (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
    cdf[i] = (ls >= 0.0) ? 0.0 :
             ((ls <= fperace::LOG_FLOOR) ? 1.0 : -std::expm1(ls));
  }
  return Rcpp::List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                            _["n_solves"] = static_cast<int>(C.e.size()));
}
