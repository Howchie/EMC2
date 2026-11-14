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

// .emc_alloc_factor(): the factor on a coordinate's rho from x expected moves
// (floor(x) moves, and one more with probability x - floor(x)), given one
// move's factor g and the moves' correlation a. Plain numeric vectors of one
// length (or length one) without missing values; NULL otherwise, which the R
// form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_factor_cpp(SEXP g, SEXP x, SEXP a) {
  auto plain = [](SEXP v) {
    return TYPEOF(v) == REALSXP && XLENGTH(v) >= 1 && !Rf_isObject(v) &&
           Rf_getAttrib(v, R_NamesSymbol) == R_NilValue &&
           Rf_getAttrib(v, R_DimSymbol) == R_NilValue;
  };
  if (!plain(g) || !plain(x) || !plain(a)) return R_NilValue;
  const R_xlen_t ng = XLENGTH(g), nx = XLENGTH(x), na = XLENGTH(a);
  const R_xlen_t n = std::max(ng, std::max(nx, na));
  if ((ng != 1 && ng != n) || (nx != 1 && nx != n) || (na != 1 && na != n)) return R_NilValue;
  const double *gp = REAL(g), *xp = REAL(x), *ap = REAL(a);
  for (R_xlen_t i = 0; i < ng; ++i) if (ISNAN(gp[i])) return R_NilValue;
  for (R_xlen_t i = 0; i < nx; ++i) if (ISNAN(xp[i])) return R_NilValue;
  for (R_xlen_t i = 0; i < na; ++i) if (ISNAN(ap[i])) return R_NilValue;
  // k moves' factor: g^k, or the saturating form where the moves are less
  // correlated than one move's factor.
  auto moves = [](double gi, double ai, double k) {
    volatile double edge = gi - 1e-12;
    if (!(ai < edge)) return k == 2.0 ? gi * gi : R_pow(gi, k);
    volatile double num = 1 - gi;
    volatile double den = 1 - ai;
    volatile double ratio = num / den;
    volatile double pw = k == 2.0 ? ai * ai : R_pow(ai, k);
    volatile double rest = 1 - pw;
    volatile double prod = ratio * rest;
    volatile double out = 1 - prod;
    return static_cast<double>(out);
  };
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    const double gi = gp[ng == 1 ? 0 : i], xi = xp[nx == 1 ? 0 : i], ai = ap[na == 1 ? 0 : i];
    volatile double m = std::floor(xi);
    volatile double f = xi - m;
    volatile double keep = 1 - f;
    volatile double up = m + 1;
    volatile double lo = moves(gi, ai, m);
    volatile double hi = moves(gi, ai, up);
    volatile double t1 = keep * lo;
    volatile double t2 = f * hi;
    volatile double sum = t1 + t2;
    REAL(out)[i] = sum;
  }
  UNPROTECT(1);
  return out;
}

namespace {

// A unit's numeric field without names, dimensions or class, of length one
// or `k` (any length for k < 0); false otherwise or for missing values.
bool plain_field(SEXP v, R_xlen_t k) {
  if (TYPEOF(v) != REALSXP || Rf_isObject(v) || Rf_getAttrib(v, R_NamesSymbol) != R_NilValue ||
      Rf_getAttrib(v, R_DimSymbol) != R_NilValue || XLENGTH(v) < 1)
    return false;
  if (k >= 0 && XLENGTH(v) != 1 && XLENGTH(v) != k) return false;
  for (R_xlen_t i = 0; i < XLENGTH(v); ++i)
    if (ISNAN(REAL(v)[i])) return false;
  return true;
}

bool kind_is(SEXP u, const char* kind) {
  SEXP k = unit_field(u, "kind");
  return TYPEOF(k) == STRSXP && XLENGTH(k) == 1 && std::strcmp(CHAR(STRING_ELT(k, 0)), kind) == 0;
}

}  // namespace

// .emc_alloc_predict(): each unit's lag-one rho under a solved plan, given
// the plan's allocations x and particle rates q, and the collapsed move's
// factor `fac` on the units numbered `idx` (NULL when it does not run). NULL
// for units whose fields are not plain numbers, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_predict_cpp(const Rcpp::List& units, SEXP sol_x, SEXP sol_q, SEXP idx, SEXP fac) {
  const int n = units.size();
  if ((sol_x != R_NilValue && TYPEOF(sol_x) != REALSXP) ||
      (sol_q != R_NilValue && TYPEOF(sol_q) != REALSXP) || TYPEOF(idx) != INTSXP ||
      (fac != R_NilValue && (TYPEOF(fac) != REALSXP || XLENGTH(fac) != XLENGTH(idx))))
    return R_NilValue;
  const R_xlen_t nx = sol_x == R_NilValue ? 0 : XLENGTH(sol_x);
  const R_xlen_t nq = sol_q == R_NilValue ? 0 : XLENGTH(sol_q);
  Rcpp::List out(n);
  for (int i = 0; i < n; ++i) {
    SEXP u = units[i];
    SEXP r = unit_field(u, "rho"), b = unit_field(u, "b_share");
    if (!plain_field(r, -1)) return R_NilValue;
    const R_xlen_t k = XLENGTH(r);
    if (b != R_NilValue && !plain_field(b, k)) return R_NilValue;
    // The collapsed move's factor on the unit, if it moves it.
    double moved = NA_REAL;
    if (fac != R_NilValue)
      for (R_xlen_t j = 0; j < XLENGTH(idx); ++j)
        if (INTEGER(idx)[j] == i + 1) {
          moved = REAL(fac)[j];
          if (ISNAN(moved)) return R_NilValue;
          break;
        }
    double x = 0.0;
    if (nx > i && R_FINITE(REAL(sol_x)[i])) x = REAL(sol_x)[i];
    const bool hmc = kind_is(u, "hmc"), subject = hmc || kind_is(u, "subject");
    const bool shear = kind_is(u, "shear"), fixed = kind_is(u, "fixed");
    // The unit's own moves' factor on each coordinate.
    std::vector<double> own(k, 1.0);
    bool scaled = false;
    double tail = NA_REAL;
    Rcpp::RObject f_keep;
    if (shear) {
      SEXP g = unit_field(u, "g"), a = unit_field(u, "a"), rt = unit_field(u, "rho_t");
      if (!plain_field(g, k) || (a != R_NilValue && !plain_field(a, k))) return R_NilValue;
      Rcpp::NumericVector xv = Rcpp::NumericVector::create(x);
      f_keep = emc_alloc_factor_cpp(g, xv, a == R_NilValue ? g : a);
      SEXP f = f_keep;
      if (f == R_NilValue || (XLENGTH(f) != 1 && XLENGTH(f) != k)) return R_NilValue;
      for (R_xlen_t c = 0; c < k; ++c) own[c] = REAL(f)[XLENGTH(f) == 1 ? 0 : c];
      scaled = true;
      if (rt != R_NilValue) {
        if (TYPEOF(rt) != REALSXP || XLENGTH(rt) != 1) return R_NilValue;
        if (R_FINITE(REAL(rt)[0])) {
          if (XLENGTH(f) != 1 && k != XLENGTH(f)) return R_NilValue;
          tail = REAL(rt)[0];
        }
      }
    } else if (subject) {
      SEXP p = hmc ? R_NilValue : unit_field(u, "p"), h = unit_field(u, "h");
      SEXP rm = unit_field(u, "r_min");
      if ((p != R_NilValue && !plain_field(p, k)) || (h != R_NilValue && !plain_field(h, k)))
        return R_NilValue;
      double r_min = 0.0;
      if (rm != R_NilValue && !hmc) {
        if (!plain_field(rm, 1)) return R_NilValue;
        r_min = REAL(rm)[0];
      }
      const double q = (p == R_NilValue || nq <= i) ? 1.0 : REAL(sol_q)[i];
      if (ISNAN(q)) return R_NilValue;
      double rate = 0.0;
      if (h != R_NilValue) {
        volatile double total = x + r_min;
        rate = std::min(1.0, static_cast<double>(total));
      }
      volatile double excess = q + rate;
      volatile double over = excess - 1;
      const double both = std::max(0.0, static_cast<double>(over));
      volatile double a0 = 1 - q;
      volatile double a1 = a0 - rate;
      volatile double none = a1 + both;
      volatile double only_p = q - both;
      volatile double only_h = rate - both;
      for (R_xlen_t c = 0; c < k; ++c) {
        const double P = p == R_NilValue ? 1.0 : REAL(p)[XLENGTH(p) == 1 ? 0 : c];
        const double H = h == R_NilValue ? 1.0 : REAL(h)[XLENGTH(h) == 1 ? 0 : c];
        volatile double t1 = only_p * P;
        volatile double t2 = only_h * H;
        volatile double t3 = both * P;
        volatile double t4 = t3 * H;
        volatile double s1 = none + t1;
        volatile double s2 = s1 + t2;
        volatile double s3 = s2 + t4;
        own[c] = s3;
      }
      scaled = true;
    } else if (!fixed) {
      SEXP kind = unit_field(u, "kind");
      if (TYPEOF(kind) != STRSXP || XLENGTH(kind) != 1) return R_NilValue;
    }
    Rcpp::NumericVector res(k);
    for (R_xlen_t c = 0; c < k; ++c) {
      const double bc = b == R_NilValue ? 0.0 : REAL(b)[XLENGTH(b) == 1 ? 0 : c];
      volatile double rho = REAL(r)[c];
      if (!ISNAN(moved)) {
        volatile double m = rho * moved;
        rho = m;
      }
      volatile double rest = 1 - bc;
      volatile double v = rest * rho;
      if (scaled) {
        volatile double w = v * own[c];
        v = w;
      }
      volatile double sum = bc + v;
      double val = sum;
      if (!ISNAN(tail)) {
        volatile double t = tail * own[c];
        if (ISNAN(val) || ISNAN(static_cast<double>(t))) return R_NilValue;
        val = std::max(val, static_cast<double>(t));
      }
      res[c] = val;
    }
    out[i] = res;
  }
  return out;
}

// The shear units' needs with their tails (.emc_alloc_need_with_tail): where
// a unit carries a finite tail rho, the rounds both it and its tail need, the
// tail priced as a coordinate of its own without a between-chain share, and
// their cost. `need` holds the units' own needs (emc_alloc_shears_need_cpp)
// and the other vectors their fixed fields. NULL for a tail that is not one
// number or a need with missing values, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_shears_tail_need_cpp(const Rcpp::List& units, const Rcpp::List& need,
                                    const Rcpp::NumericVector& g, const Rcpp::NumericVector& a,
                                    const Rcpp::NumericVector& cost,
                                    const Rcpp::NumericVector& cap,
                                    const Rcpp::NumericVector& floor,
                                    const Rcpp::NumericVector& rs) {
  const int n = units.size(), n_t = rs.size();
  if (need.size() != n || g.size() != n || a.size() != n || cost.size() != n ||
      cap.size() != n || floor.size() != n)
    return R_NilValue;
  const double inf = std::numeric_limits<double>::infinity();
  Rcpp::List out(n);
  for (int i = 0; i < n; ++i) {
    SEXP rt = unit_field(units[i], "rho_t");
    out[i] = need[i];
    if (rt == R_NilValue) continue;
    if (TYPEOF(rt) != REALSXP || XLENGTH(rt) != 1) return R_NilValue;
    if (!R_FINITE(REAL(rt)[0])) continue;
    const Rcpp::List tail = emc_alloc_shears_need_cpp(
        Rcpp::NumericVector::create(REAL(rt)[0]), Rcpp::NumericVector::create(0.0),
        Rcpp::NumericVector::create(g[i]), Rcpp::NumericVector::create(a[i]),
        Rcpp::NumericVector::create(cost[i]), Rcpp::NumericVector::create(cap[i]),
        Rcpp::NumericVector::create(floor[i]), rs);
    const Rcpp::List own = need[i];
    const Rcpp::NumericVector x0 = own["x"], xt = Rcpp::as<Rcpp::List>(tail[0])["x"];
    if (x0.size() != n_t || xt.size() != n_t) return R_NilValue;
    Rcpp::NumericVector X(n_t), C(n_t);
    for (int t = 0; t < n_t; ++t) {
      if (ISNAN(x0[t]) || ISNAN(xt[t])) return R_NilValue;
      const double x = std::max(x0[t], xt[t]);
      X[t] = x;
      C[t] = std::isfinite(x) ? x * cost[i] : inf;
    }
    out[i] = Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("q") = own["q"],
                                Rcpp::Named("cost") = C);
  }
  return out;
}

namespace {

#pragma GCC push_options
#pragma GCC optimize("no-fast-math", "fp-contract=off")

// The mean net squared jump over twice the variance that moves of factor g
// and correlation x give, for iterations with k moves in proportions w[k-1].
// A function of its own: the file's floating-point settings do not reach a
// lambda's body.
double saturation_curve(const std::vector<double>& w, double g, double x) {
  long double s = 0.0;
  for (std::size_t j = 0; j < w.size(); ++j) {
    const double k = static_cast<double>(j + 1);
    volatile double pw = k == 2.0 ? x * x : R_pow(x, k);
    volatile double rest = 1 - pw;
    volatile double term = w[j] * rest;
    s += term;
  }
  volatile double total = static_cast<double>(s);
  volatile double left = 1 - g;
  volatile double scaled = left * total;
  volatile double den = 1 - x;
  return scaled / den;
}

}  // namespace

// .emc_alloc_saturation(): the moves' correlation a in [0, g] per coordinate
// from D, the mean net squared jump over twice the variance, and K, the
// iterations' move counts (one row per coordinate, or one row for all),
// solved by bisection. NULL for inputs with missing values or without rows
// for every coordinate, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_saturation_cpp(SEXP g, SEXP D, SEXP K) {
  if (TYPEOF(g) != REALSXP || TYPEOF(D) != REALSXP || !Rf_isMatrix(K) ||
      (TYPEOF(K) != INTSXP && TYPEOF(K) != REALSXP))
    return R_NilValue;
  const R_xlen_t n = XLENGTH(g), nd = XLENGTH(D);
  const int nr = Rf_nrows(K), nc = Rf_ncols(K);
  if (nr < 1 || nc < 1 || (nr != 1 && nr < n)) return R_NilValue;
  for (R_xlen_t i = 0; i < n; ++i)
    if (ISNAN(REAL(g)[i])) return R_NilValue;
  std::vector<double> counts(static_cast<std::size_t>(nr) * nc);
  for (R_xlen_t e = 0; e < XLENGTH(K); ++e) {
    if (TYPEOF(K) == INTSXP) {
      if (INTEGER(K)[e] == NA_INTEGER) return R_NilValue;
      counts[e] = INTEGER(K)[e];
    } else {
      if (ISNAN(REAL(K)[e])) return R_NilValue;
      counts[e] = REAL(K)[e];
    }
  }
  SEXP out = PROTECT(Rf_duplicate(g));
  std::vector<double> w(nc);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int row = nr == 1 ? 0 : static_cast<int>(i);
    long double total = 0.0, later = 0.0;
    for (int j = 0; j < nc; ++j) {
      const double v = counts[row + static_cast<std::size_t>(j) * nr];
      total += v;
      if (j > 0) later += v;
    }
    const double gi = REAL(g)[i];
    if (static_cast<double>(total) <= 0 || static_cast<double>(later) <= 0 || i >= nd ||
        !R_FINITE(REAL(D)[i]) || gi >= 1 - 1e-9)
      continue;
    const double Di = REAL(D)[i], sum_w = static_cast<double>(total);
    for (int j = 0; j < nc; ++j) w[j] = counts[row + static_cast<std::size_t>(j) * nr] / sum_w;
    auto h = [&](double x) { return saturation_curve(w, gi, x); };
    if (Di <= h(0.0)) {
      REAL(out)[i] = 0.0;
      continue;
    }
    if (Di >= h(gi)) continue;
    double lo = 0.0, hi = gi;
    for (int it = 0; it < 60; ++it) {
      const double mid = (lo + hi) / 2;
      if (h(mid) < Di) lo = mid; else hi = mid;
    }
    REAL(out)[i] = (lo + hi) / 2;
  }
  UNPROTECT(1);
  return out;
}

#pragma GCC pop_options
