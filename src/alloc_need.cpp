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

// Every subject unit's cheapest allocation for each target rs[t], in one call:
// row u of x, q and cost holds unit u's HMC rate above its r_min, particle
// rate and cost (Inf where out of reach), as .emc_alloc_unit_need gives them
// one unit at a time. Units are packed: coordinates off[u]..off[u+1]-1 of rho,
// p and level (p is 1 where the unit has no particle factor), h where has_h;
// q_fixed and r_fixed are NA unless held. A unit without p or h runs its
// particle step (q = 1).
// [[Rcpp::export]]
Rcpp::List emc_alloc_subjects_need_cpp(const Rcpp::NumericVector& rho,
                                       const Rcpp::NumericVector& p,
                                       const Rcpp::NumericVector& h,
                                       const Rcpp::NumericVector& level,
                                       const Rcpp::IntegerVector& off,
                                       const Rcpp::LogicalVector& has_p,
                                       const Rcpp::LogicalVector& has_h,
                                       const Rcpp::NumericVector& q_fixed,
                                       const Rcpp::NumericVector& r_fixed,
                                       const Rcpp::NumericVector& r_min,
                                       const Rcpp::NumericVector& c_p,
                                       const Rcpp::NumericVector& c_h,
                                       const Rcpp::NumericVector& q_grid,
                                       const Rcpp::NumericVector& rs, bool refresh) {
  const int n_u = has_p.size(), n_t = rs.size();
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::NumericMatrix X(n_u, n_t), Q(n_u, n_t), C(n_u, n_t);
  std::vector<double> qs, rn;
  for (int u = 0; u < n_u; ++u) {
    const bool hp = has_p[u], hh = has_h[u];
    qs.clear();
    if (!hp || !hh) qs.push_back(1.0);
    else if (!Rcpp::NumericVector::is_na(q_fixed[u])) qs.push_back(q_fixed[u]);
    else qs.assign(q_grid.begin(), q_grid.end());
    const int n_q = qs.size();
    const double rmin = hh ? r_min[u] : 0.0;
    const bool held = hh && !Rcpp::NumericVector::is_na(r_fixed[u]);
    rn.assign(n_q, 0.0);
    // At targets at or above every coordinate's level each coordinate's
    // target is its level, so the allocation is the same at all of them:
    // computed at the first and copied.
    double top = -inf;
    for (int i = off[u]; i < off[u + 1]; ++i) top = std::max(top, level[i]);
    int first_top = -1;
    for (int t = 0; t < n_t; ++t) {
      if (!refresh && rs[t] >= top) {
        if (first_top >= 0) {
          X(u, t) = X(u, first_top);
          Q(u, t) = Q(u, first_top);
          C(u, t) = C(u, first_top);
          continue;
        }
        first_top = t;
      }
      for (int a = 0; a < n_q; ++a) {
        const double q = qs[a];
        double r = std::max(rmin, refresh ? 1.0 - q : 0.0);
        for (int i = off[u]; i < off[u + 1]; ++i) {
          const double K = refresh ? rs[t] : std::min(rs[t], level[i]);
          const double P = p[i];
          const double alone = 1.0 - q * (1.0 - P);
          if (!(rho[i] * alone > K * (1.0 + 1e-12))) continue;
          if (!hh || h[i] >= 1.0) { r = inf; continue; }
          const double H = h[i];
          const double tt = K / rho[i];
          const double r1 = (alone - tt) / (1.0 - H);
          const double knee = 1.0 - q;
          double val;
          if (std::isnan(r1)) val = r1;
          else if (r1 <= knee) val = r1;
          else val = knee + (q * P + (1.0 - q) * H - tt) / (P * (1.0 - H));
          if (std::isnan(r) || std::isnan(val)) r = std::numeric_limits<double>::quiet_NaN();
          else if (val > r) r = val;
        }
        if (r > 1 + 1e-9) r = inf;
        else if (std::isfinite(r)) r = std::min(r, 1.0);
        if (held) r = (r <= r_fixed[u] + 1e-9) ? r_fixed[u] : inf;
        rn[a] = r;
      }
      // The cheapest particle rate, the first (highest) of equal costs.
      int best = 0;
      double best_c = inf;
      bool any = false;
      for (int a = 0; a < n_q; ++a) {
        if (!std::isfinite(rn[a])) continue;
        const double c = qs[a] * c_p[u] + (rn[a] - rmin) * c_h[u];
        if (!any || c < best_c) { best = a; best_c = c; any = true; }
      }
      X(u, t) = rn[best] - rmin;
      Q(u, t) = qs[best];
      C(u, t) = any ? best_c : inf;
    }
  }
  return Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("q") = Q,
                            Rcpp::Named("cost") = C);
}
