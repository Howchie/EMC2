#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

// The least HMC rate a subject needs, for each particle rate qs[a] (rows) and
// target rs[t] (columns): the largest over coordinates i of the rate at which
// rho_i F(q, r) falls to the coordinate's target, where F, the kernels' joint
// factor under coupled selection, falls at slope 1 - h_i up to the knee
// r = 1 - q and at slope p_i (1 - h_i) beyond it. Coordinates already at
// their target contribute nothing; one the HMC cannot move (no h, or h_i >= 1)
// makes the target unreachable. NaN propagates as R's pmax does.
// [[Rcpp::export]]
Rcpp::NumericMatrix emc_alloc_subject_need_cpp(const Rcpp::NumericVector& rho,
                                               const Rcpp::NumericVector& b_share,
                                               const Rcpp::NumericVector& p,
                                               Rcpp::Nullable<Rcpp::NumericVector> h_,
                                               const Rcpp::NumericVector& level,
                                               const Rcpp::NumericVector& qs,
                                               const Rcpp::NumericVector& rs,
                                               double r_min) {
  const int n_i = rho.size(), n_q = qs.size(), n_t = rs.size();
  const bool has_h = h_.isNotNull();
  Rcpp::NumericVector h = has_h ? Rcpp::NumericVector(h_) : Rcpp::NumericVector(0);
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::NumericMatrix out(n_q, n_t);
  for (int t = 0; t < n_t; ++t) {
    for (int a = 0; a < n_q; ++a) {
      const double q = qs[a];
      double rn = r_min;
      for (int i = 0; i < n_i; ++i) {
        const double K = std::min(rs[t], level[i]);
        const double P = p[i];
        const double b = b_share[i];
        const double alone = 1.0 - q * (1.0 - P);
        const double base = b + (1.0 - b) * rho[i] * alone;
        if (!(base > K * (1.0 + 1e-12))) continue;
        if (!has_h || h[i] >= 1.0) { rn = inf; continue; }
        const double H = h[i];
        const double den = (1.0 - b) * rho[i];
        if (!(den > 0.0)) { rn = inf; continue; }
        const double tt = (K - b) / den;
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
                                       const Rcpp::NumericVector& b_share,
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
                                       const Rcpp::NumericVector& rs) {
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
      if (rs[t] >= top) {
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
        double r = rmin;
        for (int i = off[u]; i < off[u + 1]; ++i) {
          const double K = std::min(rs[t], level[i]);
          const double P = p[i];
          const double b = b_share[i];
          const double alone = 1.0 - q * (1.0 - P);
          const double base = b + (1.0 - b) * rho[i] * alone;
          if (!(base > K * (1.0 + 1e-12))) continue;
          if (!hh || h[i] >= 1.0) { r = inf; continue; }
          const double H = h[i];
          const double den = (1.0 - b) * rho[i];
          if (!(den > 0.0)) { r = inf; continue; }
          const double tt = (K - b) / den;
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

// A family's activation probability at each target (column of X): one minus
// the product over its units (rows) of one minus min(x + add, 1), a
// non-finite rate counting as one. The product is accumulated in long double,
// as R's prod() does, so the result equals the R expression exactly.
// [[Rcpp::export]]
Rcpp::NumericVector emc_alloc_act_term_cpp(const Rcpp::NumericMatrix& X,
                                           const Rcpp::NumericVector& add) {
  const int n_r = X.nrow(), n_t = X.ncol();
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::NumericVector out(n_t);
  for (int t = 0; t < n_t; ++t) {
    long double s = 1.0L;
    for (int i = 0; i < n_r; ++i) {
      double v = X(i, t) + add[i];
      if (std::isnan(v) || v == -inf || v > 1.0) v = 1.0;
      s *= static_cast<long double>(1.0 - v);
    }
    out[t] = 1.0 - static_cast<double>(s);
  }
  return out;
}

// emc_alloc_act_term_cpp with the rows given as a list of vectors.
// [[Rcpp::export]]
Rcpp::NumericVector emc_alloc_act_term_list_cpp(const Rcpp::List& xs,
                                                const Rcpp::NumericVector& add) {
  const int n_r = xs.size();
  if (!n_r) return Rcpp::NumericVector(0);
  std::vector<Rcpp::NumericVector> rows(n_r);
  for (int i = 0; i < n_r; ++i) rows[i] = Rcpp::as<Rcpp::NumericVector>(xs[i]);
  const int n_t = rows[0].size();
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::NumericVector out(n_t);
  for (int t = 0; t < n_t; ++t) {
    long double s = 1.0L;
    for (int i = 0; i < n_r; ++i) {
      double v = rows[i][t] + add[i];
      if (std::isnan(v) || v == -inf || v > 1.0) v = 1.0;
      s *= static_cast<long double>(1.0 - v);
    }
    out[t] = 1.0 - static_cast<double>(s);
  }
  return out;
}

// The shear pricing below reproduces R's arithmetic exactly: no fused
// multiply-adds or fast-math rewrites, which the package flags allow.
#pragma GCC push_options
#pragma GCC optimize("no-fast-math", "fp-contract=off")

// R's x^y, which is R_pow (std::pow differs in the last bit).
static inline double alloc_r_pow(double x, double y) { return R_pow(x, y); }

// .emc_alloc_factor_k for one k.
static inline double alloc_factor_k(double g, double a, double k) {
  if (a < g - 1e-12) return 1.0 - (1.0 - g) / (1.0 - a) * (1.0 - alloc_r_pow(a, k));
  return alloc_r_pow(g, k);
}

// .emc_alloc_factor_inverse for one target t < 1, given the move's g < 1.
static inline double alloc_factor_inverse(double g, double a, double t) {
  const double inf = std::numeric_limits<double>::infinity();
  double m;
  if (a >= g - 1e-12) {
    if (!(t > 0)) return inf;
    m = std::ceil(std::log(t) / std::log(g) - 1e-9);
  } else {
    const double c = (1.0 - g) / (1.0 - a);
    const double s = 1.0 - (1.0 - t) / c;
    if (!(s > 0)) return inf;
    m = a <= 0 ? 1.0 : std::ceil(std::log(s) / std::log(a) - 1e-9);
  }
  if (!(m >= 1.0)) m = 1.0;
  const double f0 = alloc_factor_k(g, a, m - 1.0), f1 = alloc_factor_k(g, a, m);
  double v = (f0 - t) / (f0 - f1);
  if (!std::isnan(v)) v = std::min(1.0, std::max(0.0, v));
  return (m - 1.0) + v;
}

// Every shear unit's rounds and cost at each target rs, as
// .emc_alloc_unit_need prices one (a list per unit: x, q, cost): the floor
// where the coordinate already meets the target, else the least rounds that
// reach it, raised to the floor and Inf beyond the cap.
// [[Rcpp::export]]
Rcpp::List emc_alloc_shears_need_cpp(const Rcpp::NumericVector& rho,
                                     const Rcpp::NumericVector& b_share,
                                     const Rcpp::NumericVector& g,
                                     const Rcpp::NumericVector& a,
                                     const Rcpp::NumericVector& cost,
                                     const Rcpp::NumericVector& cap,
                                     const Rcpp::NumericVector& floor,
                                     const Rcpp::NumericVector& rs) {
  const int n_u = rho.size(), n_t = rs.size();
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::List out(n_u);
  for (int i = 0; i < n_u; ++i) {
    const bool reach = g[i] < 1 - 1e-9;
    Rcpp::NumericVector X(n_t), C(n_t), Q(n_t, 1.0);
    for (int t = 0; t < n_t; ++t) {
      double x = floor[i];
      const double base = b_share[i] + (1.0 - b_share[i]) * rho[i];
      if (base > rs[t]) {
        const double den = (1.0 - b_share[i]) * rho[i];
        const double target = den > 0.0 ? (rs[t] - b_share[i]) / den : -inf;
        x = reach ? alloc_factor_inverse(g[i], a[i], target) : inf;
        if (!std::isnan(x) && x < floor[i]) x = floor[i];
        if (!std::isfinite(x) || x > cap[i] + 1e-9) x = inf;
      }
      X[t] = x;
      C[t] = std::isfinite(x) ? x * cost[i] : inf;
    }
    out[i] = Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("q") = Q,
                                Rcpp::Named("cost") = C);
  }
  return out;
}

// Every fixed unit's allocation at each target: nothing to buy, and out of
// reach wherever its rho (largest coordinate rho_max) is above the target.
// [[Rcpp::export]]
Rcpp::List emc_alloc_fixed_need_cpp(const Rcpp::NumericVector& rho_max,
                                    const Rcpp::NumericVector& rs) {
  const int n_u = rho_max.size(), n_t = rs.size();
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::List out(n_u);
  for (int i = 0; i < n_u; ++i) {
    Rcpp::NumericVector X(n_t), C(n_t), Q(n_t, 1.0);
    for (int t = 0; t < n_t; ++t) C[t] = rho_max[i] > rs[t] ? inf : 0.0;
    out[i] = Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("q") = Q,
                                Rcpp::Named("cost") = C);
  }
  return out;
}

#pragma GCC pop_options

namespace {

// u$name of a unit; NULL when absent.
SEXP unit_field(SEXP u, const char* name) {
  SEXP names = Rf_getAttrib(u, R_NamesSymbol);
  if (TYPEOF(u) != VECSXP || TYPEOF(names) != STRSXP) return R_NilValue;
  for (R_xlen_t i = 0; i < XLENGTH(u); ++i)
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(u, i);
  return R_NilValue;
}

// A unit's scalar field as a number, `if_null` when absent; false for a
// field that is not one number.
bool unit_scalar(SEXP u, const char* name, double if_null, double& out) {
  SEXP v = unit_field(u, name);
  if (v == R_NilValue) {
    out = if_null;
    return true;
  }
  if ((TYPEOF(v) != REALSXP && TYPEOF(v) != INTSXP && TYPEOF(v) != LGLSXP) || XLENGTH(v) != 1)
    return false;
  out = Rf_asReal(v);
  return true;
}

}  // namespace

// The subject units' cheapest allocations at each target, from the units
// themselves: emc_alloc_subjects_need_cpp() on their packed fields, each
// unit's row returned as its own list (x, q, cost). NULL for units whose
// fields are not in the expected form, which R then packs itself.
// [[Rcpp::export]]
SEXP emc_alloc_units_subjects_need_cpp(const Rcpp::List& units, const Rcpp::NumericVector& q_grid,
                                       const Rcpp::NumericVector& rs) {
  const int n = units.size();
  std::vector<double> rho, b_share, p, h, level;
  Rcpp::IntegerVector off(n + 1);
  Rcpp::LogicalVector has_p(n), has_h(n);
  Rcpp::NumericVector q_fixed(n), r_fixed(n), r_min(n), c_p(n), c_h(n);
  auto numeric = [](SEXP v) { return TYPEOF(v) == REALSXP || TYPEOF(v) == INTSXP; };
  for (int j = 0; j < n; ++j) {
    SEXP u = units[j];
    SEXP r = unit_field(u, "rho"), b = unit_field(u, "b_share"), pp = unit_field(u, "p"),
         hh = unit_field(u, "h");
    if (!numeric(r)) return R_NilValue;
    const R_xlen_t k = XLENGTH(r);
    if ((b != R_NilValue && (!numeric(b) || XLENGTH(b) != k)) ||
        (pp != R_NilValue && (!numeric(pp) || XLENGTH(pp) != k)) ||
        (hh != R_NilValue && (!numeric(hh) || XLENGTH(hh) != k)))
      return R_NilValue;
    double qf, rf, rm, cp, ch;
    if (!unit_scalar(u, "q_fixed", NA_REAL, qf) || !unit_scalar(u, "r_fixed", NA_REAL, rf) ||
        !unit_scalar(u, "r_min", 0.0, rm) || !unit_scalar(u, "c_p", 0.0, cp) ||
        !unit_scalar(u, "c_h", 0.0, ch))
      return R_NilValue;
    Rcpp::NumericVector rv = Rcpp::as<Rcpp::NumericVector>(r);
    Rcpp::NumericVector bv, pv, hv;
    if (b != R_NilValue) bv = Rcpp::as<Rcpp::NumericVector>(b);
    if (pp != R_NilValue) pv = Rcpp::as<Rcpp::NumericVector>(pp);
    if (hh != R_NilValue) hv = Rcpp::as<Rcpp::NumericVector>(hh);
    for (R_xlen_t i = 0; i < k; ++i) {
      const double f = pp == R_NilValue ? 1.0 : pv[i];
      // The unit's rho under its particle step alone.
      volatile double stepped = rv[i] * f;
      double base = stepped;
      if (b != R_NilValue) {
        volatile double rest = 1 - bv[i];
        volatile double scaled = rest * stepped;
        volatile double sum = bv[i] + scaled;
        base = sum;
      }
      rho.push_back(rv[i]);
      b_share.push_back(b == R_NilValue ? 0.0 : bv[i]);
      p.push_back(f);
      h.push_back(hh == R_NilValue ? 1.0 : hv[i]);
      level.push_back(base);
    }
    off[j + 1] = off[j] + static_cast<int>(k);
    has_p[j] = pp != R_NilValue;
    has_h[j] = hh != R_NilValue;
    q_fixed[j] = qf;
    r_fixed[j] = rf;
    r_min[j] = rm;
    c_p[j] = cp;
    c_h[j] = ch;
  }
  const Rcpp::List got = emc_alloc_subjects_need_cpp(
      Rcpp::NumericVector(rho.begin(), rho.end()),
      Rcpp::NumericVector(b_share.begin(), b_share.end()),
      Rcpp::NumericVector(p.begin(), p.end()), Rcpp::NumericVector(h.begin(), h.end()),
      Rcpp::NumericVector(level.begin(), level.end()), off, has_p, has_h, q_fixed, r_fixed, r_min,
      c_p, c_h, q_grid, rs);
  const Rcpp::NumericMatrix x = got["x"], q = got["q"], cost = got["cost"];
  const int n_t = x.ncol();
  Rcpp::List out(n);
  for (int j = 0; j < n; ++j) {
    Rcpp::NumericVector xj(n_t), qj(n_t), cj(n_t);
    for (int t = 0; t < n_t; ++t) {
      xj[t] = x(j, t);
      qj[t] = q(j, t);
      cj[t] = cost(j, t);
    }
    out[j] = Rcpp::List::create(Rcpp::_["x"] = xj, Rcpp::_["q"] = qj, Rcpp::_["cost"] = cj);
  }
  return out;
}
