#ifndef EMC2_DDM_VARIABILITY_H
#define EMC2_DDM_VARIABILITY_H

// Wiener first-passage density averaged over a uniform start point (sz) and a
// uniform non-decision time (st0); drift variability (sv) stays in closed form.
//
// Every rule below has a fixed number of nodes whose positions move smoothly
// with the parameters, so the log density is smooth in the parameters. An
// adaptive integrator is not: its decision to subdivide depends on the
// parameters, and the result jumps wherever that decision changes.
//
// Conventions follow dwiener: the density is that of absorption at the LOWER
// boundary, d is the start point as a fraction of a measured from that
// boundary, and v is the drift in that frame. An upper-boundary response maps
// to d = 1 - w and v = -v.

#include <algorithm>
#include <cmath>
#include "utility_functions.h"
#include "ddm_functions_inline.h"
#include "composite_functions.h"
#include "wald_functions.h"
#include "gl_quad.h"

namespace ddmvar {

// u / a^2 below which the start point is integrated in closed form. Below it
// the start-point integrand can be a narrow Gaussian at the edge of the start
// range, which a fixed rule misses; above it the Gaussian is wide enough for a
// fixed rule over the large-time series, which needs only a few terms there.
// Both are accurate far below the likelihood's other errors at the switch, so
// changing regime does not show in the density.
constexpr double SMALL_TIME_MAX = 0.1;
constexpr int N_W = 12;       // start-point nodes in the large-time regime
constexpr int N_U = 16;       // nodes over a non-decision-time window
constexpr int N_RAMP = 12;    // nodes on each sloped side of the binned kernel
constexpr int N_NARROW = 8;   // nodes for a Gaussian integral over a short interval
constexpr double NARROW_SPAN = 3.0;  // span x (1 + max |z|) below which an interval is short
constexpr double SERIES_ERR = 1e-12;
constexpr double TERM_CUTOFF = 25.0;  // drop series terms below exp(-25) of the largest
constexpr double TERM_FLOOR = 1.3887943864964021e-11;  // exp(-TERM_CUTOFF)
constexpr int MAX_TERMS = 40;

inline const GLRule& rule(int n) {
  static thread_local const GLRule* cached[64] = {nullptr};
  if (cached[n] == nullptr) cached[n] = &gl_get_rule(n);
  return *cached[n];
}

// Moments of the standard normal over [lo, hi] in log scale:
//   p0 = int phi,  p1l = int (z - lo) phi,  p1r = int (hi - z) phi.
// All three are non-negative, so a linear function of fixed sign integrates
// as a sum of two non-negative terms anchored at its smaller endpoint.
struct NormalMoments {
  double p0, p1l, p1r;
};

inline NormalMoments normal_moments_upper(double lo, double hi);

inline NormalMoments normal_moments(double lo, double hi) {
  const double span = hi - lo;
  const double reach = std::max(std::fabs(lo), std::fabs(hi));
  if (span * (1.0 + reach) < NARROW_SPAN) {
    // Short interval: log phi changes by a few units at most, so a fixed rule
    // is accurate and avoids the tail subtractions below.
    const GLRule& g = rule(N_NARROW);
    const double half = 0.5 * span, mid = lo + half;
    double lw[N_NARROW], lmax = R_NegInf;
    double z[N_NARROW];
    for (int j = 0; j < N_NARROW; ++j) {
      z[j] = mid + half * g.x[j];
      lw[j] = std::log(half * g.w[j]) + log_phi_std(z[j]);
      lmax = std::max(lmax, lw[j]);
    }
    double s0 = 0.0, sl = 0.0, sr = 0.0;
    for (int j = 0; j < N_NARROW; ++j) {
      const double e = std::exp(lw[j] - lmax);
      s0 += e;
      sl += (z[j] - lo) * e;
      sr += (hi - z[j]) * e;
    }
    return {lmax + std::log(s0), lmax + std::log(sl), lmax + std::log(sr)};
  }
  if (lo >= 0.0) return normal_moments_upper(lo, hi);
  if (hi <= 0.0) {
    const NormalMoments m = normal_moments_upper(-hi, -lo);
    return {m.p0, m.p1r, m.p1l};
  }
  // Straddles zero: neither tail is extreme, so the direct forms are stable.
  const double p0 = 1.0 - fast_norm_phi(-hi) - fast_norm_phi(lo);
  const double phl = std::exp(log_phi_std(lo)), phh = std::exp(log_phi_std(hi));
  return {std::log(p0), std::log(phl - phh - lo * p0),
          std::log(hi * p0 + phh - phl)};
}

// 0 <= lo < hi and the interval is long relative to the local tail scale.
// Written through the Mills ratio R and its decrement D = 1 - zR, with every
// term scaled by phi(lo) so that nothing underflows before the log.
inline NormalMoments normal_moments_upper(double lo, double hi) {
  const double r_lo = mills_ratio_std(lo), d_lo = mills_ratio_decrement(lo);
  double r_hi = 0.0, d_hi = 0.0, ratio = 0.0;
  if (R_FINITE(hi)) {
    ratio = std::exp(log_phi_std(hi) - log_phi_std(lo));
    r_hi = mills_ratio_std(hi);
    d_hi = mills_ratio_decrement(hi);
  }
  const double base = log_phi_std(lo);
  const double p0 = r_lo - ratio * r_hi;
  const double p1l = d_lo - ratio * (d_hi + (hi - lo) * r_hi);
  const double p1r = (hi - lo) * r_lo - d_lo + ratio * d_hi;
  return {base + std::log(p0), base + std::log(p1l), base + std::log(p1r)};
}

// log of the integral over d in [d_lo, d_hi] of the lower-boundary density at
// decision time u, using the small-time series term by term. With
// T = 1 + sv^2 u and u' = u / a^2, term k is (d + 2k) times a Gaussian in d
// with mean M_k = -a v u' - 2kT and variance T u', so each term integrates in
// closed form through normal_moments.
inline double log_start_integral_small(double u, double a, double v, double sv,
                                       double d_lo, double d_hi) {
  const double up = u / (a * a);
  const double T = 1.0 + sv * sv * u;
  const double S = std::sqrt(T * up);
  const double log_S = std::log(S);

  const double log_width = std::log(d_hi - d_lo);
  // The exponent left after the Gaussian is completed, with the e^{-v^2 u / 2T}
  // drift factor of the density already folded in.
  auto exponent = [&](int k) {
    return 2.0 * k * a * v + 2.0 * k * k * sv * sv * a * a;
  };
  // Upper bound on log |term k|: largest |d + 2k| times the interval length
  // times the Gaussian's largest value on it. Lets the loops stop without
  // evaluating a term they would discard.
  auto bound = [&](int k) {
    const double M = -a * v * up - 2.0 * k * T;
    const double lo = (d_lo - M) / S, hi = (d_hi - M) / S;
    const double nearest = lo > 0.0 ? lo : (hi < 0.0 ? hi : 0.0);
    const double reach = k >= 0 ? d_hi + 2.0 * k : -(d_lo + 2.0 * k);
    return exponent(k) + std::log(reach) + log_width + log_phi_std(nearest);
  };
  auto term = [&](int k) -> double {
    const double M = -a * v * up - 2.0 * k * T;
    const double lo = (d_lo - M) / S, hi = (d_hi - M) / S;
    const NormalMoments m = normal_moments(lo, hi);
    // (d + 2k) is positive on the interval for k >= 0 and negative for k < 0;
    // anchor the linear factor at its endpoint nearest zero.
    double log_j;
    if (k >= 0) {
      const double at_lo = d_lo + 2.0 * k;
      log_j = at_lo > 0.0 ? log_sum_exp(std::log(at_lo) + m.p0, log_S + m.p1l)
                          : log_S + m.p1l;
    } else {
      const double at_hi = -(d_hi + 2.0 * k);
      log_j = at_hi > 0.0 ? log_sum_exp(std::log(at_hi) + m.p0, log_S + m.p1r)
                          : log_S + m.p1r;
    }
    return exponent(k) + log_S + log_j;
  };

  // Terms shrink monotonically away from k = 0 and k = -1, the two images
  // nearest the start interval.
  double plus = term(0), minus = R_NegInf;
  double top = plus;
  for (int k = 1; k <= MAX_TERMS && bound(k) > top - TERM_CUTOFF; ++k) {
    const double t = term(k);
    plus = log_sum_exp(plus, t);
    top = std::max(top, t);
  }
  for (int k = -1; k >= -MAX_TERMS && bound(k) > top - TERM_CUTOFF; --k) {
    const double t = term(k);
    minus = log_sum_exp(minus, t);
    top = std::max(top, t);
  }
  if (!(plus > minus)) return R_NegInf;
  const double series = minus == R_NegInf ? plus
    : plus + std::log1p(-std::exp(minus - plus));
  // The constant factors of the density: a^-2 T^-1/2 (2 pi u'^3)^-1/2, with
  // the Gaussian completion's sqrt(2 pi) cancelling the series' 1/sqrt(2 pi).
  return series - 2.0 * std::log(a) - 0.5 * std::log(T) - 1.5 * std::log(up);
}

// Start-point rule for one trial: the fixed nodes the large-time regime uses
// and everything about them that does not depend on decision time, so each
// node's sin(k pi d) table is built once per trial rather than once per u.
// Without start-point variability the rule is the single point d.
constexpr int K_MAX = 12;  // large-time terms needed for u / a^2 >= SMALL_TIME_MAX
struct StartRule {
  double a, v, sv, d, d_half;
  int n;
  double node[N_W], weight[N_W];  // weights average: they sum to one
  double sine[K_MAX + 1][N_W];    // sin(k pi d_j), node index fastest
  double tilt[N_W];               // (sv^2 a^2 d^2 / 2 - a v d); divided by T it is the drift exponent
  // With sv = 0 the drift tilt e^{-a v d} does not depend on u, so the node
  // sum collapses to one coefficient per series term.
  double coef[K_MAX + 1], coef_shift;
};

inline void start_rule(StartRule& r, double a, double v, double sv, double d,
                       double d_half) {
  r.a = a; r.v = v; r.sv = sv; r.d = d; r.d_half = d_half;
  r.n = d_half > 0.0 ? N_W : 1;
  const GLRule& g = rule(N_W);
  for (int j = 0; j < r.n; ++j) {
    r.node[j] = r.n == 1 ? d : d + d_half * g.x[j];
    r.weight[j] = r.n == 1 ? 1.0 : 0.5 * g.w[j];
    r.tilt[j] = 0.5 * sv * sv * a * a * r.node[j] * r.node[j] - a * v * r.node[j];
    const double c2 = 2.0 * std::cos(M_PI * r.node[j]);
    r.sine[0][j] = 0.0;
    r.sine[1][j] = std::sin(M_PI * r.node[j]);
    for (int k = 2; k <= K_MAX; ++k)
      r.sine[k][j] = c2 * r.sine[k - 1][j] - r.sine[k - 2][j];
  }
  if (sv == 0.0) {
    r.coef_shift = R_NegInf;
    for (int j = 0; j < r.n; ++j) r.coef_shift = std::max(r.coef_shift, r.tilt[j]);
    for (int k = 1; k <= K_MAX; ++k) r.coef[k] = 0.0;
    for (int j = 0; j < r.n; ++j) {
      const double w = r.weight[j] * std::exp(r.tilt[j] - r.coef_shift);
      for (int k = 1; k <= K_MAX; ++k) r.coef[k] += w * r.sine[k][j];
    }
  }
}

// Large-time series f = pi exp(lg1(d)) sum_k k e^{-k^2 pi^2 u'/2} sin(k pi d),
// averaged over the start-point nodes. Returns the log of the average.
inline double log_start_average_large(double u, const StartRule& r) {
  const double a = r.a, v = r.v, sv = r.sv;
  const double up = u / (a * a);
  const double T = 1.0 + sv * sv * u;
  const double decay = 0.5 * M_PI * M_PI * up;
  // Weight of term k relative to term 1, k e^{-(k^2 - 1) decay}, by the
  // recurrence e^{-(k^2 - 1) decay} = e^{-((k-1)^2 - 1) decay} q^{2k - 1}.
  // Terms stop once k times the weight (sin(k pi d) ~ k sin(pi d) near a
  // boundary) is below e^-TERM_CUTOFF.
  const double q = std::exp(-decay), q2 = q * q;
  double weight[K_MAX + 1];
  int K = 1;
  weight[1] = 1.0;
  double ratio = 1.0, step = q;  // e^{-(k^2 - 1) decay} and q^{2k - 1}
  while (K < K_MAX) {
    const int k = K + 1;
    step *= q2;
    ratio *= step;
    if (k * k * ratio < TERM_FLOOR) break;
    weight[k] = k * ratio;
    K = k;
  }
  const double base = M_LNPI - decay - 2.0 * std::log(a) - 0.5 * std::log(T) -
    v * v * u / (2.0 * T);
  if (sv == 0.0) {
    double sum = 0.0;
    for (int k = 1; k <= K; ++k) sum += weight[k] * r.coef[k];
    return sum > 0.0 ? base + r.coef_shift + std::log(sum) : R_NegInf;
  }
  const int n = r.n;
  const double inv_T = 1.0 / T;
  double sum[N_W], expo[N_W];
  for (int j = 0; j < n; ++j) sum[j] = r.sine[1][j];
  for (int k = 2; k <= K; ++k) {
    const double wk = weight[k];
    #pragma omp simd
    for (int j = 0; j < n; ++j) sum[j] += wk * r.sine[k][j];
  }
  double top = R_NegInf;
  for (int j = 0; j < n; ++j) {
    expo[j] = r.tilt[j] * inv_T;
    top = std::max(top, expo[j]);
  }
  double acc = 0.0;
  #pragma omp simd reduction(+:acc)
  for (int j = 0; j < n; ++j)
    acc += r.weight[j] * std::max(sum[j], 0.0) * std::exp(expo[j] - top);
  return acc > 0.0 ? base + top + std::log(acc) : R_NegInf;
}

// log of the start-point-averaged lower-boundary density at decision time u.
inline double log_density_u(double u, const StartRule& r) {
  if (!(u > 0.0)) return R_NegInf;
  if (u >= SMALL_TIME_MAX * r.a * r.a) return log_start_average_large(u, r);
  if (r.d_half <= 0.0)
    return dwiener_inline(-u, r.a, r.v, r.d, r.sv, SERIES_ERR, 0, 1);
  return log_start_integral_small(u, r.a, r.v, r.sv, r.d - r.d_half,
                                  r.d + r.d_half) - std::log(2.0 * r.d_half);
}

// log of the integral of kernel(u) * density(u) over [u_lo, u_hi], with
// kernel linear in u. The rule runs in y = log(1 + sqrt(u) / c), c = ONSET_SCALE a:
// on a window short against a^2 this is the sqrt(u) rule, which spreads nodes
// over the leading edge; on a long one the density's mass sits within ~a^2 of
// onset and the logarithm keeps nodes there instead of across the empty tail.
// The density vanishes to all orders at u = 0, so clipping the panel there
// keeps the result smooth.
constexpr double ONSET_SCALE = 0.1;
constexpr int N_PANEL_MAX = 32;
template <typename Kernel>
inline double log_panel(double u_lo, double u_hi, int n, const Kernel& kernel,
                        const StartRule& r) {
  u_lo = std::max(u_lo, 0.0);
  if (!(u_hi > u_lo)) return R_NegInf;
  const GLRule& g = rule(n);
  const double c = ONSET_SCALE * r.a;
  const double y_lo = std::log1p(std::sqrt(u_lo) / c);
  const double y_hi = std::log1p(std::sqrt(u_hi) / c);
  const double half = 0.5 * (y_hi - y_lo), mid = y_lo + half;
  double term[N_PANEL_MAX], scale[N_PANEL_MAX], top = R_NegInf;
  for (int j = 0; j < n; ++j) {
    const double y = mid + half * g.x[j];
    const double z = c * std::expm1(y);
    const double u = z * z;
    scale[j] = 2.0 * z * c * std::exp(y) * half * g.w[j] * kernel(u);
    term[j] = scale[j] > 0.0 ? log_density_u(u, r) : R_NegInf;
    top = std::max(top, term[j]);
  }
  if (top == R_NegInf) return R_NegInf;
  double acc = 0.0;
  for (int j = 0; j < n; ++j)
    if (term[j] != R_NegInf) acc += scale[j] * std::exp(term[j] - top);
  return acc > 0.0 ? top + std::log(acc) : R_NegInf;
}

// The rule in the lower-boundary frame of response R (1 = lower, 2 = upper).
// sz is the absolute start-point range.
inline void start_rule_for(StartRule& r, int R, double a, double v, double Z,
                           double sz, double sv) {
  if (R == 1) start_rule(r, a, v, sv, Z, 0.5 * sz);
  else start_rule(r, a, -v, sv, 1.0 - Z, 0.5 * sz);
}

// log density at response time t for onset t0 + U(0, st0).
inline double log_density(double t, int R, double a, double v, double t0,
                          double Z, double sz, double sv, double st0) {
  const double u_hi = t - t0;
  if (!(u_hi > 0.0)) return R_NegInf;
  StartRule r;
  start_rule_for(r, R, a, v, Z, sz, sv);
  if (!(st0 > 0.0)) return log_density_u(u_hi, r);
  return log_panel(u_hi - st0, u_hi, N_U, [](double) { return 1.0; }, r) -
    std::log(st0);
}

// Binned likelihood for st0 > 0. The bin and the onset are both uniform, so
// over decision time u the integrand carries their convolution: a trapezoid
// whose sloped sides have width r = min(hi - lo, st0) and whose height is r.
// It equals a box of height r between the midpoints of the two sides plus one
// correction per side. For st0 >= hi - lo the box is r times the st0 window at
// the bin midpoint, i.e. the midpoint rule.
//
// A side's correction is either exact (two half-panels, either side of its
// midpoint) or its leading term from the density at the side's two corners:
// against a box edge at the side's midpoint, a rising side adds
// -r^3 f' / 24 = -r^2 (f(top corner) - f(bottom corner)) / 24 and a falling
// side the negative of that. What the corner form leaves is O(r^4) relative,
// except where the side lies on the density's leading edge, near u = 0. There
// the density rises like exp(-d_min^2 / 2u), whose time scale at u is
// u^2 / d_min^2 (unit diffusion), so the corner form holds once u is well
// beyond sqrt(r) d_min, and a bin when d_min is ~0. Each side takes the exact
// form within BIN_EXACT_NEAR such lengths of u = 0 and blends to the corner
// form over BIN_EXACT_BLEND more, so the result is C1 in the parameters.
// Beyond the leading edge the density decays no faster than
// e^{-(pi^2 / 2a^2 + v^2 / 2) u}; when a side spans more than a fraction of
// that decay length the corner form fails everywhere, so both sides blend to
// exact between r times the rate BIN_DECAY_LO and BIN_DECAY_HI.
constexpr double BIN_EXACT_NEAR = 0.5;
constexpr double BIN_EXACT_BLEND = 0.5;
constexpr double BIN_DECAY_LO = 0.5;
constexpr double BIN_DECAY_HI = 1.0;
constexpr int N_HALF = 6;       // nodes on the lower half of an exact side
constexpr int N_HALF_TOP = 12;  // and on its upper half, where an onset bin's
                                // density climbs steeply to the bin edge

inline double smoothstep(double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x * x * (3.0 - 2.0 * x);
}

// Weight of the exact form for a side whose midpoint is at decision time m.
inline double bin_side_weight(double m, double r, double a, double d_min) {
  const double scale = std::max(r, std::sqrt(r) * a * std::max(d_min, 0.0));
  return smoothstep((BIN_EXACT_NEAR + BIN_EXACT_BLEND - m / scale) / BIN_EXACT_BLEND);
}

// One signed term of the trapezoid sum, as log magnitude and sign.
struct SignedLog {
  double log_value, sign;
};

// Appends the correction of the side on [c, c + r] (rising when rising is
// true) with weight w on the exact form.
inline void bin_side(double c, double r, bool rising, double w,
                     const StartRule& rule_, SignedLog* terms, int& n) {
  if (!(c + r > 0.0)) return;  // the density is zero below u = 0
  const double m = c + 0.5 * r, e = c + r;
  const double sign = rising ? 1.0 : -1.0;
  if (w > 0.0) {
    const double lw = std::log(w);
    // Rising: (u - c) on [c, m], minus (e - u) on [m, e]; falling mirrors it.
    terms[n++] = {lw + log_panel(c, m, N_HALF, [&](double u) { return u - c; }, rule_),
                  sign};
    terms[n++] = {lw + log_panel(m, e, N_HALF_TOP, [&](double u) { return e - u; }, rule_),
                  -sign};
  }
  if (w < 1.0) {
    const double lc = std::log1p(-w) + 2.0 * std::log(r) - std::log(24.0);
    terms[n++] = {lc + log_density_u(e, rule_), -sign};
    terms[n++] = {lc + log_density_u(c, rule_), sign};
  }
}

// log of the response-time probability in [lo, hi] divided by the nominal
// width h, with weights w_low and w_up on the exact form of each side.
inline double log_bin_density_sides(double lo, double hi, double h,
                                    const StartRule& rule_, double t0,
                                    double st0, double w_low, double w_up) {
  const double width = hi - lo;
  if (!(width > 0.0) || !(st0 > 0.0)) return R_NegInf;
  const double r = std::min(width, st0);
  const double c_low = lo - t0 - st0, c_up = hi - t0 - r;
  SignedLog terms[9];
  int n = 0;
  terms[n++] = {std::log(r) + log_panel(c_low + 0.5 * r, c_up + 0.5 * r, N_U,
                                        [](double) { return 1.0; }, rule_), 1.0};
  bin_side(c_low, r, true, w_low, rule_, terms, n);
  bin_side(c_up, r, false, w_up, rule_, terms, n);
  double top = R_NegInf;
  for (int j = 0; j < n; ++j) top = std::max(top, terms[j].log_value);
  if (top == R_NegInf) return R_NegInf;
  double acc = 0.0;
  for (int j = 0; j < n; ++j)
    if (terms[j].log_value != R_NegInf)
      acc += terms[j].sign * std::exp(terms[j].log_value - top);
  // The sum is the trapezoid integral up to the corner form's error, so it
  // is positive wherever that form is used; the guard keeps a far-from-valid
  // call finite.
  if (!(acc > 0.0)) return terms[0].log_value - std::log(h * st0);
  return top + std::log(acc) - std::log(h * st0);
}

inline double log_bin_density(double lo, double hi, double h, int R, double a,
                              double v, double t0, double Z, double sz,
                              double sv, double st0) {
  const double width = hi - lo;
  if (!(width > 0.0) || !(st0 > 0.0)) return R_NegInf;
  StartRule rule_;
  start_rule_for(rule_, R, a, v, Z, sz, sv);
  const double r = std::min(width, st0);
  const double d_min = rule_.d - rule_.d_half;
  const double c_low = lo - t0 - st0, c_up = hi - t0 - r;
  const double decay = r * (0.5 * M_PI * M_PI / (a * a) + 0.5 * v * v);
  const double coarse = 1.0 - smoothstep((decay - BIN_DECAY_LO) / (BIN_DECAY_HI - BIN_DECAY_LO));
  return log_bin_density_sides(
    lo, hi, h, rule_, t0, st0,
    1.0 - coarse * (1.0 - bin_side_weight(c_low + 0.5 * r, r, a, d_min)),
    1.0 - coarse * (1.0 - bin_side_weight(c_up + 0.5 * r, r, a, d_min)));
}

// log P(T <= t, R) for onset t0 + U(0, st0): the decision-time density
// against P(onset <= t - u), which is 1 up to t - t0 - st0 and then falls
// linearly to 0 at t - t0. The body runs from 0, so it gets more nodes than
// an st0 window. Past CDF_HORIZON a^2 every series term has decayed by at
// least e^{-pi^2 CDF_HORIZON / 2} from its value at a^2, so the mass there is
// dropped and t = Inf gives the response probability.
constexpr int N_CDF = 24;
constexpr int N_CDF_RAMP = 16;
constexpr double CDF_HORIZON = 8.0;
inline double log_cdf(double t, int R, double a, double v, double t0,
                      double Z, double sz, double sv, double st0) {
  const double x = t - t0;
  if (!(x > 0.0)) return R_NegInf;
  StartRule r;
  start_rule_for(r, R, a, v, Z, sz, sv);
  const double horizon = CDF_HORIZON * a * a;
  if (!(st0 > 0.0))
    return log_panel(0.0, std::min(x, horizon), N_CDF, [](double) { return 1.0; }, r);
  const double u_ramp = x - st0;
  double acc = R_NegInf;
  if (u_ramp > 0.0)
    acc = log_panel(0.0, std::min(u_ramp, horizon), N_CDF,
                    [&](double) { return st0; }, r);
  const double lo = std::max(u_ramp, 0.0), hi = std::min(x, horizon);
  if (hi > lo)
    acc = log_sum_exp(acc, log_panel(lo, hi, N_CDF_RAMP,
                                     [&](double u) { return x - u; }, r));
  return acc - std::log(st0);
}

}  // namespace ddmvar

#endif
