#include <Rcpp.h>
#include <cmath>
#include <limits>

// The least HMC rate a subject needs, for each particle rate qs[a] (rows) and
// target rs[t] (columns): the largest over coordinates i of the rate at which
// rho_i F(q, r) falls to the coordinate's target, where F, the kernels' joint
// factor under coupled selection, falls at slope 1 - h_i up to the knee
// r = 1 - q and at slope p_i (1 - h_i) beyond it. Coordinates already at
// their target contribute nothing; one the HMC cannot move (no h, or h_i >= 1)
// makes the target unreachable. NaN propagates as R's pmax does.
// [[Rcpp::export]]
Rcpp::NumericMatrix emc_alloc_subject_need_cpp(const Rcpp::NumericVector& rho,
                                               const Rcpp::NumericVector& p,
                                               Rcpp::Nullable<Rcpp::NumericVector> h_,
                                               const Rcpp::NumericVector& level,
                                               const Rcpp::NumericVector& qs,
                                               const Rcpp::NumericVector& rs,
                                               double r_min, bool refresh) {
  const int n_i = rho.size(), n_q = qs.size(), n_t = rs.size();
  const bool has_h = h_.isNotNull();
  Rcpp::NumericVector h = has_h ? Rcpp::NumericVector(h_) : Rcpp::NumericVector(0);
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::NumericMatrix out(n_q, n_t);
  for (int t = 0; t < n_t; ++t) {
    for (int a = 0; a < n_q; ++a) {
      const double q = qs[a];
      double rn = std::max(r_min, refresh ? 1.0 - q : 0.0);
      for (int i = 0; i < n_i; ++i) {
        const double K = refresh ? rs[t] : std::min(rs[t], level[i]);
        const double P = p[i];
        const double alone = 1.0 - q * (1.0 - P);
        if (!(rho[i] * alone > K * (1.0 + 1e-12))) continue;
        if (!has_h || h[i] >= 1.0) { rn = inf; continue; }
        const double H = h[i];
        const double tt = K / rho[i];
        const double r1 = (alone - tt) / (1.0 - H);
        const double knee = 1.0 - q;
        double val;
        if (std::isnan(r1)) val = r1;
        else if (r1 <= knee) val = r1;
        else val = knee + (q * P + (1.0 - q) * H - tt) / (P * (1.0 - H));
        if (std::isnan(rn) || std::isnan(val)) rn = std::numeric_limits<double>::quiet_NaN();
        else if (val > rn) rn = val;
      }
      out(a, t) = rn;
    }
  }
  return out;
}
