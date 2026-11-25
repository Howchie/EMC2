#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

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

// The pricing below reproduces R's arithmetic exactly: no fused
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
// .emc_alloc_unit_need prices one (a list per unit: x, cost): the floor
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
    Rcpp::NumericVector X(n_t), C(n_t);
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
    out[i] = Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("cost") = C);
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
    Rcpp::NumericVector X(n_t), C(n_t);
    for (int t = 0; t < n_t; ++t) C[t] = rho_max[i] > rs[t] ? inf : 0.0;
    out[i] = Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("cost") = C);
  }
  return out;
}

namespace {

// u$name of a unit; NULL when absent.
SEXP unit_field(SEXP u, const char* name) {
  SEXP names = Rf_getAttrib(u, R_NamesSymbol);
  if (TYPEOF(u) != VECSXP || TYPEOF(names) != STRSXP) return R_NilValue;
  for (R_xlen_t i = 0; i < XLENGTH(u); ++i)
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(u, i);
  return R_NilValue;
}

}  // namespace

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
// the plan's allocations x and the collapsed move's factor `fac` on the units
// numbered `idx` (NULL when it does not run). NULL for units whose fields are
// not plain numbers, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_predict_cpp(const Rcpp::List& units, SEXP sol_x, SEXP idx, SEXP fac) {
  const int n = units.size();
  if ((sol_x != R_NilValue && TYPEOF(sol_x) != REALSXP) || TYPEOF(idx) != INTSXP ||
      (fac != R_NilValue && (TYPEOF(fac) != REALSXP || XLENGTH(fac) != XLENGTH(idx))))
    return R_NilValue;
  const R_xlen_t nx = sol_x == R_NilValue ? 0 : XLENGTH(sol_x);
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
    // The unit's own rounds' factor on each coordinate.
    std::vector<double> own(k, 1.0);
    bool scaled = false;
    Rcpp::RObject f_keep;
    if (kind_is(u, "shear")) {
      SEXP g = unit_field(u, "g"), a = unit_field(u, "a");
      if (!plain_field(g, k) || (a != R_NilValue && !plain_field(a, k))) return R_NilValue;
      Rcpp::NumericVector xv = Rcpp::NumericVector::create(x);
      f_keep = emc_alloc_factor_cpp(g, xv, a == R_NilValue ? g : a);
      SEXP f = f_keep;
      if (f == R_NilValue || (XLENGTH(f) != 1 && XLENGTH(f) != k)) return R_NilValue;
      for (R_xlen_t c = 0; c < k; ++c) own[c] = REAL(f)[XLENGTH(f) == 1 ? 0 : c];
      scaled = true;
    } else if (!kind_is(u, "fixed")) {
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
      res[c] = sum;
    }
    out[i] = res;
  }
  return out;
}

namespace {

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

namespace {

// A unit's field at coordinate c, a field of length one serving every
// coordinate.
inline double field_at(SEXP v, R_xlen_t c) { return REAL(v)[XLENGTH(v) == 1 ? 0 : c]; }

// b + (1 - b) rho, or rho without a between-chain share.
inline double with_share(SEXP b, R_xlen_t c, double rho) {
  if (b == R_NilValue) return rho;
  const double bc = field_at(b, c);
  volatile double rest = 1 - bc;
  volatile double v = rest * rho;
  volatile double sum = bc + v;
  return sum;
}

// A shear unit's factor at x rounds; false when it is left to the R form.
bool shear_factor(SEXP u, R_xlen_t k, double x, std::vector<double>& out) {
  SEXP g = unit_field(u, "g"), a = unit_field(u, "a");
  if (!plain_field(g, k) || (a != R_NilValue && !plain_field(a, k)) || ISNAN(x)) return false;
  Rcpp::NumericVector xv = Rcpp::NumericVector::create(x);
  Rcpp::RObject keep = emc_alloc_factor_cpp(g, xv, a == R_NilValue ? g : a);
  SEXP f = keep;
  if (f == R_NilValue || (XLENGTH(f) != 1 && XLENGTH(f) != k)) return false;
  out.assign(k, 0.0);
  for (R_xlen_t c = 0; c < k; ++c) out[c] = field_at(f, c);
  return true;
}

// A unit's optional field that must be one double when present.
bool one_double(SEXP v) { return TYPEOF(v) == REALSXP && XLENGTH(v) == 1 && !ISNAN(REAL(v)[0]); }

}  // namespace

// .emc_alloc_range_parts(): each unit's rho with nothing allocated (levels)
// and the least rho its own moves reach (best), the largest over its
// coordinates. NULL for units whose fields are not plain numbers of the
// coordinates' length or of length one, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_range_parts_cpp(const Rcpp::List& units, SEXP max_rounds) {
  if (!one_double(max_rounds) || Rf_getAttrib(units, R_NamesSymbol) != R_NilValue)
    return R_NilValue;
  const int n = units.size();
  Rcpp::NumericVector levels(n), best(n);
  std::vector<double> f_floor, f_cap;
  for (int i = 0; i < n; ++i) {
    SEXP u = units[i];
    SEXP r = unit_field(u, "rho"), b = unit_field(u, "b_share");
    if (!plain_field(r, -1)) return R_NilValue;
    const R_xlen_t k = XLENGTH(r);
    if (b != R_NilValue && !plain_field(b, k)) return R_NilValue;
    const bool shear = kind_is(u, "shear");
    if (!shear && !kind_is(u, "fixed")) return R_NilValue;
    double level = R_NegInf, least = R_NegInf;
    auto raise = [](double& top, double v) { if (v > top) top = v; };
    bool bad = false;
    auto seen = [&bad](double v) { if (ISNAN(v)) bad = true; return v; };
    if (shear) {
      SEXP fl = unit_field(u, "floor"), cap = unit_field(u, "cap");
      if ((fl != R_NilValue && !one_double(fl)) || (cap != R_NilValue && !one_double(cap)))
        return R_NilValue;
      const bool floored = fl != R_NilValue && REAL(fl)[0] > 0;
      if (floored && !shear_factor(u, k, REAL(fl)[0], f_floor)) return R_NilValue;
      if (!shear_factor(u, k, REAL(cap == R_NilValue ? max_rounds : cap)[0], f_cap))
        return R_NilValue;
      for (R_xlen_t c = 0; c < k; ++c) {
        const double rho = REAL(r)[c];
        volatile double at_floor = floored ? rho * f_floor[c] : rho;
        raise(level, seen(with_share(b, c, at_floor)));
        volatile double at_cap = rho * f_cap[c];
        raise(least, seen(with_share(b, c, at_cap)));
      }
    } else {
      for (R_xlen_t c = 0; c < k; ++c) raise(level, seen(with_share(b, c, REAL(r)[c])));
      least = level;
    }
    if (bad) return R_NilValue;
    levels[i] = level;
    best[i] = least;
  }
  return Rcpp::List::create(Rcpp::Named("levels") = levels, Rcpp::Named("best") = best);
}

namespace {

// .emc_alloc_row_var() of the columns [from, from + n) of `x` (p rows, column
// major): each row's sample variance, sums accumulated as R's rowMeans() and
// rowSums() do.
void row_var(const double* x, int p, R_xlen_t n, std::vector<double>& out) {
  out.assign(p, 0.0);
  for (int i = 0; i < p; ++i) {
    long double s = 0.0;
    for (R_xlen_t j = 0; j < n; ++j) s += x[i + j * p];
    s /= n;
    const double mean = static_cast<double>(s);
    long double ss = 0.0;
    for (R_xlen_t j = 0; j < n; ++j) {
      volatile double d = x[i + j * p] - mean;
      volatile double sq = d * d;
      ss += sq;
    }
    volatile double total = static_cast<double>(ss);
    volatile double den = static_cast<double>(n - 1);
    out[i] = total / den;
  }
}

}  // namespace

// .emc_alloc_within_share(): from the chains' draws of p coordinates (rows),
// each coordinate's mean within-chain variance, its variance over all draws,
// and the share of that between chains. Unnamed; NULL for draws that are not
// finite matrices of doubles with one set of row names, for a chain with
// fewer than two draws, or where a variance is not finite, which the R form
// handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_within_share_cpp(const Rcpp::List& chains) {
  const int n_c = chains.size();
  if (n_c < 1) return R_NilValue;
  int p = -1;
  R_xlen_t n_all = 0;
  SEXP first_names = R_NilValue;
  for (int c = 0; c < n_c; ++c) {
    SEXP x = chains[c];
    if (TYPEOF(x) != REALSXP || !Rf_isMatrix(x) || Rf_isObject(x)) return R_NilValue;
    if (p < 0) p = Rf_nrows(x);
    if (Rf_nrows(x) != p || p < 1 || Rf_ncols(x) < 2) return R_NilValue;
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    SEXP rn = dn == R_NilValue ? R_NilValue : VECTOR_ELT(dn, 0);
    if (c == 0) first_names = rn;
    else if (!R_compute_identical(rn, first_names, 16)) return R_NilValue;
    for (R_xlen_t e = 0; e < XLENGTH(x); ++e)
      if (!R_FINITE(REAL(x)[e])) return R_NilValue;
    n_all += Rf_ncols(x);
  }
  if (n_all < 2) return R_NilValue;
  std::vector<double> all(static_cast<std::size_t>(p) * n_all);
  R_xlen_t at = 0;
  for (int c = 0; c < n_c; ++c) {
    SEXP x = chains[c];
    std::copy(REAL(x), REAL(x) + XLENGTH(x), all.begin() + at);
    at += XLENGTH(x);
  }
  std::vector<double> total, one;
  row_var(all.data(), p, n_all, total);
  std::vector<long double> sum(p, 0.0L);
  int counted = 0;
  for (int c = 0; c < n_c; ++c) {
    SEXP x = chains[c];
    if (Rf_ncols(x) < 2) continue;
    row_var(REAL(x), p, Rf_ncols(x), one);
    for (int i = 0; i < p; ++i) {
      if (ISNAN(one[i])) return R_NilValue;
      sum[i] += one[i];
    }
    ++counted;
  }
  if (!counted) return R_NilValue;
  Rcpp::NumericVector var(p), tot(p), share(p);
  for (int i = 0; i < p; ++i) {
    const double within = static_cast<double>(sum[i] / counted);
    if (!R_FINITE(within) || !R_FINITE(total[i])) return R_NilValue;
    volatile double gap = total[i] - within;
    const double between = gap > 0 ? static_cast<double>(gap) : 0.0;
    double b;
    if (total[i] == 0) {
      b = 0.0;
    } else {
      volatile double ratio = between / total[i];
      b = ratio;
    }
    if (!R_FINITE(b)) return R_NilValue;
    var[i] = within;
    tot[i] = total[i];
    share[i] = b > 1 ? 1.0 : (b < 0 ? 0.0 : b);
  }
  return Rcpp::List::create(Rcpp::Named("var") = var, Rcpp::Named("total") = tot,
                            Rcpp::Named("b_share") = share);
}

// What .emc_alloc_solve_now() reads off its units before pricing them: each
// unit's exploration floor (x0) and cost per move (cost), and the shear units
// (numbered from one). NULL for named units or fields that are not one
// double, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_unit_summary_cpp(const Rcpp::List& units) {
  if (Rf_getAttrib(units, R_NamesSymbol) != R_NilValue) return R_NilValue;
  const int n = units.size();
  Rcpp::NumericVector x0(n), cost(n);
  std::vector<int> shear;
  auto one = [](SEXP v) {
    return TYPEOF(v) == REALSXP && XLENGTH(v) == 1 && !Rf_isObject(v);
  };
  for (int i = 0; i < n; ++i) {
    SEXP u = units[i];
    if (TYPEOF(u) != VECSXP) return R_NilValue;
    const bool is_shear = kind_is(u, "shear");
    SEXP fl = unit_field(u, "floor"), co = unit_field(u, "cost");
    if (is_shear && fl != R_NilValue) {
      if (!one(fl)) return R_NilValue;
      x0[i] = REAL(fl)[0];
    }
    if (co != R_NilValue) {
      if (!one(co)) return R_NilValue;
      cost[i] = REAL(co)[0];
    }
    if (is_shear) shear.push_back(i + 1);
  }
  return Rcpp::List::create(Rcpp::Named("x0") = x0, Rcpp::Named("cost") = cost,
                            Rcpp::Named("shear") = Rcpp::wrap(shear));
}

// One candidate collapsed rate's plan in .emc_alloc_solve_now(): the cost at
// each target, base + joint + fixed + the moved units' costs (need[[i]]$cost,
// summed in order from zero) + the shear's activation term (NULL without
// one), and the target with the least tau * cost. list(k, cost, objective) at
// that target, k zero where no objective is a number. NULL for vectors that
// are not doubles of the targets' length, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_candidate_cpp(SEXP base, double joint, double fixed, const Rcpp::List& need,
                             SEXP act_shear, SEXP tau) {
  if (TYPEOF(base) != REALSXP || TYPEOF(tau) != REALSXP) return R_NilValue;
  const R_xlen_t n = XLENGTH(tau);
  auto fits = [n](SEXP v) { return TYPEOF(v) == REALSXP && XLENGTH(v) == n; };
  if (n < 1 || !fits(base) || (act_shear != R_NilValue && !fits(act_shear)))
    return R_NilValue;
  std::vector<const double*> costs(need.size());
  for (R_xlen_t i = 0; i < need.size(); ++i) {
    SEXP c = unit_field(need[i], "cost");
    if (!fits(c)) return R_NilValue;
    costs[i] = REAL(c);
  }
  R_xlen_t best = -1;
  double best_obj = 0.0, best_cost = 0.0;
  for (R_xlen_t t = 0; t < n; ++t) {
    volatile double moved = 0.0;
    for (const double* c : costs) {
      volatile double next = moved + c[t];
      moved = next;
    }
    volatile double act = 0.0;
    if (act_shear != R_NilValue) {
      volatile double next = act + REAL(act_shear)[t];
      act = next;
    }
    volatile double c1 = REAL(base)[t] + joint;
    volatile double c2 = c1 + fixed;
    volatile double c3 = c2 + moved;
    volatile double cost = c3 + act;
    volatile double obj = REAL(tau)[t] * cost;
    const double o = obj;
    if (!ISNAN(o) && (best < 0 || o < best_obj)) {
      best = t;
      best_obj = o;
      best_cost = cost;
    }
  }
  return Rcpp::List::create(Rcpp::Named("k") = static_cast<int>(best + 1),
                            Rcpp::Named("cost") = best_cost,
                            Rcpp::Named("objective") = best_obj);
}

// .emc_alloc_rho(): the lag-one correlation 1 - esjd / (2 var) held within
// [0, 1], the variance floored at the least positive number. NULL for
// anything but two plain vectors of doubles of one length without missing
// values, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_rho_cpp(SEXP esjd, SEXP var) {
  auto plain = [](SEXP v) {
    return TYPEOF(v) == REALSXP && !Rf_isObject(v) &&
           Rf_getAttrib(v, R_DimSymbol) == R_NilValue;
  };
  if (!plain(esjd) || !plain(var) || XLENGTH(esjd) != XLENGTH(var) || XLENGTH(esjd) < 1)
    return R_NilValue;
  const R_xlen_t n = XLENGTH(esjd);
  const double tiny = std::numeric_limits<double>::min();
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t i = 0; i < n; ++i) {
    const double e = REAL(esjd)[i], v = REAL(var)[i];
    if (ISNAN(e) || ISNAN(v)) {
      UNPROTECT(1);
      return R_NilValue;
    }
    volatile double den = 2 * (v > tiny ? v : tiny);
    volatile double ratio = e / den;
    volatile double rho = 1 - ratio;
    const double r = rho;
    if (ISNAN(r)) {
      UNPROTECT(1);
      return R_NilValue;
    }
    REAL(out)[i] = r < 0 ? 0.0 : (r > 1 ? 1.0 : r);
  }
  UNPROTECT(1);
  return out;
}

// .emc_alloc_lsd_draws(): half the log of the diagonals of draws k of a
// covariance array, floored at the least positive number (parameters x
// draws, without names). NULL for anything but a cube of doubles with square
// slices, positions within it and no missing values, which the R form
// handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_lsd_draws_cpp(SEXP V, SEXP k) {
  SEXP dim = Rf_getAttrib(V, R_DimSymbol);
  if (TYPEOF(V) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 3 ||
      (TYPEOF(k) != INTSXP && TYPEOF(k) != REALSXP))
    return R_NilValue;
  const int p = INTEGER(dim)[0], n_draws = INTEGER(dim)[2];
  const R_xlen_t n = XLENGTH(k);
  if (p < 1 || INTEGER(dim)[1] != p || n < 1) return R_NilValue;
  const double tiny = std::numeric_limits<double>::min();
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * p;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, p, static_cast<int>(n)));
  for (R_xlen_t j = 0; j < n; ++j) {
    double at = TYPEOF(k) == INTSXP
      ? (INTEGER(k)[j] == NA_INTEGER ? NA_REAL : static_cast<double>(INTEGER(k)[j]))
      : REAL(k)[j];
    if (!(at >= 1) || !(at <= n_draws) || at != std::floor(at)) {
      UNPROTECT(1);
      return R_NilValue;
    }
    const R_xlen_t d = static_cast<R_xlen_t>(at) - 1;
    for (int i = 0; i < p; ++i) {
      const double v = REAL(V)[d * slice + i + static_cast<R_xlen_t>(i) * p];
      if (ISNAN(v)) {
        UNPROTECT(1);
        return R_NilValue;
      }
      volatile double l = std::log(v > tiny ? v : tiny);
      REAL(out)[j * p + i] = 0.5 * l;
    }
  }
  UNPROTECT(1);
  return out;
}

#pragma GCC pop_options
