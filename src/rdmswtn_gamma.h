#ifndef EMC2_RDMSWTN_GAMMA_H
#define EMC2_RDMSWTN_GAMMA_H

// Point-start Wald first passage with Gamma(k, rate tau) drift: log f, log F
// and log S without quadrature.  Derivations: Math/RDMSWTN_Derivation.tex,
// section "Gamma Drift Variability".  Notation:
//   p = tau s / sqrt(t), delta = d / (s sqrt(t)), z = p - delta, q = p / delta,
//   Y ~ Gamma(k, p) is the drift in units of s / sqrt(t),
//   G_nu(x) = J_nu(x) / Gamma(nu),  J_nu(x) = int_0^inf y^(nu-1) e^(-y^2/2 - x y) dy,
//   D_nu = E_{Gamma(nu,p)} phi(Y - delta) = p^nu e^(-delta^2/2) G_nu(z) / sqrt(2 pi),
//   f = (delta / t) D_k.
// Branches:
//   NB      q <= 1: negative-binomial mixture of inverse-gamma laws.
//   MOMENT  q > 1:  S = (1/p) sum_{m>=2} (1 - r^(m-1)) D_(k+m),  r = 1 - 2/q,
//                   F = 2 p^k e^(-delta^2/2) / sqrt(2 pi) sum_{n even} W_n I_n(delta),
//                   W_n = J_(k+n)(z) / Gamma(k),  I_n(x) = G_(n+1)(x).
//   NARROW  sd(Y) small: Taylor series of the conditional laws about E Y
//                   against the exact central moments of Gamma(k, p).
//   WIDE    g_Y slowly varying near delta: Taylor series against the
//                   moments of an exponentially modified Gaussian.
//   FAR_TAIL the smaller tail below e^LOG_NEGLIGIBLE: upper bounds.
// Every sequence of cylinder moments is run downward through the
// division-free recurrence G_nu = (nu + 1) G_(nu+2) + x G_(nu+1) (x >= 0,
// where G is the minimal solution), seeded at two orders with
// kappa = nu + y0^2 >= KAPPA_MIN by a uniform Laplace expansion of log G.

#include <cmath>
#include <algorithm>
#include <limits>
#include <Rcpp.h>

namespace rdmswtn_gamma {

// a * b + c, fused where the target has a native fma (a libm call otherwise)
inline double fmadd(double a, double b, double c) {
#ifdef FP_FAST_FMA
  return std::fma(a, b, c);
#else
  return a * b + c;
#endif
}

constexpr double LOG_SQRT_2PI = 0.91893853320467274178;
constexpr double INV_SQRT_2PI = 0.39894228040143267794;
constexpr double LOG2 = 0.69314718055994530942;
constexpr double KAPPA_MIN = 64.0;
constexpr double EPS = 1e-16;
constexpr double NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double NAN_D = std::numeric_limits<double>::quiet_NaN();
constexpr double BIG = 1e200, LBIG = 460.51701859880913680;
constexpr double MID = 1e100, LMID = 230.25850929940456840;

enum Branch { NB = 1, MOMENT = 2, NARROW = 3, NB_NOSKIP = 4, FAR_TAIL = 5, WIDE = 6 };  // NB_NOSKIP: testing only
enum Want { WANT_PDF = 1, WANT_CDF = 2, WANT_SURV = 4, WANT_ALL = 7 };

struct Out {
  double logf = NAN_D, logF = NAN_D, logS = NAN_D;
  int branch = 0, terms = 0;
};

// log G_nu(x) = phi(y0) - log Gamma(nu) + log sqrt(2 pi / kappa) + log(1 + sum_j c_j(alpha) / kappa^j),
// c_j polynomials of degree 2j in b = 1 - alpha in [0, 1/2] (x >= 0), j = 1..9;
// written out so the nine polynomials evaluate in parallel.
// Stirling remainder log Gamma(y) - [(y - 1/2) log y - y + log(2 pi)/2], y >= 12.
inline double stirling(double y) {
  const double r = 1.0 / y, r2 = r * r;
  return r * (1.0 / 12 + r2 * (-1.0 / 360 + r2 * (1.0 / 1260
    + r2 * (-1.0 / 1680 + r2 * (1.0 / 1188 + r2 * (-691.0 / 360360
    + r2 * (1.0 / 156)))))));
}

// log Gamma(a), a > 0: shift to y >= 12, then Stirling.
inline double lgam(double a) {
  if (a >= 12.0) return (a - 0.5) * std::log(a) - a + LOG_SQRT_2PI + stirling(a);
  double y = a, prod = 1.0;
  while (y < 12.0) { prod *= y; y += 1.0; }
  return (y - 0.5) * std::log(y) - y + LOG_SQRT_2PI + stirling(y) - std::log(prod);
}

// Laplace saddle y0 of y^nu exp(-y^2/2 - x y): y0 (y0 + x) = nu.
inline double saddle(double nu, double x) {
  return 2.0 * nu / (x + std::sqrt(x * x + 4.0 * nu));
}
inline double kappa_of(double nu, double x) {
  const double y0 = saddle(nu, x);
  return nu + y0 * y0;
}

// log G_nu(x) = log J_nu(x) - log Gamma(nu), x >= 0, kappa >= KAPPA_MIN.
inline double logG_asym(double nu, double x) {
  const double sq = std::sqrt(x * x + 4.0 * nu);
  const double y0 = 2.0 * nu / (x + sq);
  const double kap = nu + y0 * y0, alpha = nu / kap, b = 1.0 - alpha;
  const double ik = 1.0 / kap, b2 = b * b;
  const double c1 = fmadd(b, 0.08333333333333333, fmadd(0.8333333333333334, b2, 0.08333333333333333));
  const double c2 = fmadd(b, fmadd(-2.4305555555555554, b2, 0.006944444444444444), fmadd(fmadd(5.347222222222222, b2, 0.07291666666666667), b2, 0.003472222222222222));
  const double c3 = fmadd(b, fmadd(fmadd(-63.72106481481482, b2, -0.23022762345679013), b2, -0.008043981481481482), fmadd(fmadd(fmadd(65.6520061728395, b2, 13.32349537037037), b2, -0.013483796296296296), b2, -0.0026813271604938273));
  const double c4 = fmadd(b, fmadd(fmadd(fmadd(-1767.1331661522634, b2, -98.75155928497942), b2, -0.019079539609053498), b2, -0.0009178883744855967), fmadd(fmadd(fmadd(fmadd(1195.4136123971193, b2, 774.420122813786), b2, 1.0755228427211934), b2, -0.004305716306584362), b2, -0.00022947209362139917));
  const double c5 = fmadd(b, fmadd(fmadd(fmadd(fmadd(-57280.235594028636, b2, -10429.983768164866), b2, -7.912190252861368), b2, 0.03347539705504115), b2, 0.0039201961086003335), fmadd(fmadd(fmadd(fmadd(fmadd(28889.16229959705, b2, 38850.94240290638), b2, 922.8723023967979), b2, 0.10477942949459877), b2, 0.01158848425558495), b2, 0.0007840392217200666));
  const double c6 = fmadd(b, fmadd(fmadd(fmadd(fmadd(fmadd(-2164279.7422781456, b2, -839390.6601494031), b2, -10443.414382411293), b2, -0.22468145341703102), b2, 0.007353266418378622), b2, 0.0004183688255019515), fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(870687.2526406334, b2, 1996361.486411738), b2, 157385.89701878844), b2, 74.41242200086582), b2, -0.016895430570673806), b2, 0.00205232030554688), b2, 6.972813758365857e-05));
  const double c7 = fmadd(b, fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(-93961666.01413502, b2, -63425938.491317466), b2, -2651690.6143824356), b2, -835.1891459741347), b2, -0.023560550299951993), b2, -0.05139411285101596), b2, -0.0041451650614758575), fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(31448394.33942478, b2, 109800918.2850892), b2, 18787083.906380385), b2, 138867.84589876706), b2, 1.1603771297511232), b2, -0.13196901278482948), b2, -0.016528364142715683), b2, -0.0005921664373536939));
  const double c8 = fmadd(b, fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(-4620293268.367157, b2, -4825625076.4083185), b2, -443921411.6875257), b2, -2121971.418601673), b2, -13.35014090549088), b2, -0.054908564602015056), b2, -0.009439996026865943), b2, -0.00041374327266084735), fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(1324108436.6661975, b2, 6541659127.529597), b2, 1978677921.166638), b2, 49532291.08839781), b2, 10992.5038628171), b2, -2.636334640950171), b2, -0.01949478310885037), b2, -0.0023059695549890835), b2, -5.171790908260592e-05));
  const double c9 = fmadd(b, fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(-254118477470.18774, b2, -380291299168.5134), b2, -62588997619.865204), b2, -1018071731.6067933), b2, -166490.8863974251), b2, 33.400512931609676), b2, 1.0225149482733464), b2, 0.13969170737553988), b2, 0.007555488486048785), fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(63679807593.00213, b2, 422809892301.4613), b2, 200375756754.99963), b2, 11149375950.049349), b2, 36639157.011975564), b2, 111.92390108060582), b2, 2.5800549736761313), b2, 0.4223449870814658), b2, 0.03773865399843197), b2, 0.0008394987206720873));
  const double acc = ik * fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(fmadd(c9, ik, c8), ik, c7), ik, c6), ik, c5), ik, c4), ik, c3), ik, c2), ik, c1);
  // nu log(y0 / nu) + y0^2 / 2 + log(alpha) / 2 + log(1 + acc) - stirling(nu)
  return nu * std::log(2.0 / (x + sq)) + 0.5 * y0 * y0 + 0.5 * std::log(alpha) +
    std::log1p(acc) - stirling(nu);
}

// Number of integer steps above nu needed to reach kappa >= KAPPA_MIN.
inline int seed_lift(double nu, double x) {
  return kappa_of(nu, x) >= KAPPA_MIN ? 0 : static_cast<int>(std::ceil(KAPPA_MIN - nu));
}

// Downward G recurrence g_nu = (nu + 1) g_(nu+2) + x g_(nu+1) for
// nu = top - 1, ..., bottom, from scaled seeds g1 = G_(top+1), g0 = G_top.
// The values are kept relative to exp(*lscale).
struct GWalk {
  double g0, g1, lscale;
  double nu;   // order held in g0
  double x;
  inline void step() {
    const double gn = (nu) * g1 + x * g0;   // order nu - 1: (nu - 1 + 1) g_(nu+1) + x g_nu
    g1 = g0; g0 = gn; nu -= 1.0;
  }
  inline void rescale() {
    if (g0 > BIG) { g0 /= BIG; g1 /= BIG; lscale += LBIG; }
    else if (g0 < 1.0 / BIG && g0 > 0.0) { g0 *= BIG; g1 *= BIG; lscale -= LBIG; }
  }
  inline double log_g0() const { return lscale + std::log(g0); }
};

// Start a walk holding G_top in g0, with seeds lifted to kappa >= KAPPA_MIN.
inline GWalk walk_from(double top, double x) {
  const int lift = seed_lift(top, x);
  const double ns = top + lift;
  const double l0 = logG_asym(ns, x), l1 = logG_asym(ns + 1.0, x);
  GWalk w{1.0, std::exp(l1 - l0), l0, ns, x};
  for (int i = 0; i < lift; ++i) { w.step(); if ((i & 31) == 31) w.rescale(); }
  w.rescale();
  return w;
}

// log G_nu(x) for any nu > 0, x >= 0.
inline double logG(double nu, double x) {
  if (kappa_of(nu, x) >= KAPPA_MIN) return logG_asym(nu, x);
  return walk_from(nu, x).log_g0();
}

inline double log1m_exp(double a) {  // log(1 - e^a); a >= 0 (rounding above 0) gives -Inf
  if (a >= 0.0) return NEG_INF;
  return a > -LOG2 ? std::log(-std::expm1(a)) : std::log1p(-std::exp(a));
}

// ---------------------------------------------------------------- NB branch
// q <= 1: Gamma(k, tau) = sum_n w_n Gamma(k + n, d / s^2), w_n = NB(n; k, q),
// and each component is InvGamma(alpha_n = (k + n + 1) / 2, A = delta^2 / 2)
// in time.  With inc(a) = e^-A A^a / Gamma(a + 1) and u_n = w_n inc(alpha_n):
//   f = sum_n alpha_n u_n / t,  S = sum_n u_n P_n / inc_n,  F = sum_n w_n Q_n.
// Within a parity alpha steps by one: P(alpha) = inc(alpha) + P(alpha + 1)
// turns the S sum into one over inc and cumulative weights (nb_branch), and
// Q(alpha + 1) = Q(alpha) + inc(alpha) runs upward on the scale of inc.
constexpr double TAIL_SMALL = 0.01;   // a tail below this is computed directly
constexpr int NB_MAX_TERMS = 1000000;

// sum_j x^j / ((a+1)_j) = P(a, x) / inc(a), for x < a + 1.
inline double lower_ratio(double a, double x) {
  double term = 1.0, sum = 1.0, ap = a;
  for (int j = 0; j < 100000; ++j) {
    ap += 1.0; term *= x / ap; sum += term;
    if (term <= 1e-17 * sum) break;
  }
  return sum;
}

// Q(a, x) / inc(a) = a * (Legendre continued fraction), for x >= a + 1.
inline double upper_ratio(double a, double x) {
  constexpr double TINY = 1e-300;
  double b = x + 1.0 - a, c = 1.0 / TINY, d = 1.0 / b, h = d;
  for (int i = 1; i < 100000; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b; if (std::fabs(d) < TINY) d = TINY;
    c = b + an / c; if (std::fabs(c) < TINY) c = TINY;
    d = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) <= 1e-16) break;
  }
  return a * h;
}

// log inc(a) = log(A^a e^-A / Gamma(a + 1)); the direct form loses digits
// once a log A + A is large, where R's saddle-point form takes over.
inline double log_inc(double a, double A, double lA) {
  if (a * std::fabs(lA) + A < 200.0) return a * lA - A - lgam(a + 1.0);
  return R::dgamma(A, a + 1.0, 1.0, 1);
}

// log w_n + log inc(alpha_n), the log of the summand u_n (inc(a) is the
// Gamma(a + 1, 1) density at A).
inline double nb_log_u(double n, double k, double q, double A) {
  return R::dnbinom(n, k, q, 1) + R::dgamma(A, 0.5 * (k + n + 3.0), 1.0, 1);
}

// Start of the summed window.  When A is large the u_n peak sits near
// alpha_n ~ (1 - q)^2 A, far from n = 0.  Below a start n0 with alpha_(n0) < A,
// inc(alpha) increases in alpha and Q / inc = upper_ratio increases in alpha,
// so with H = NB mass below n0 (an incomplete beta):
//   sum_(m<n0) u_m <= H inc(alpha_(n0)),   sum_(m<n0) w_m Q_m <= H inc(alpha_(n0)) upper_ratio(alpha_(n0), A),
// and the head's share of S is H less that.  u_n is log-concave in n for
// n >= 2 (the curvature of log inc, -psi'(alpha + 1) / 4, outweighs the
// (1 - k) / ((n + 2)(n + k)) of the NB weights when k < 1), so with
// rho = u_(n0) / u_(n0 - 1) > 1 also sum_(m<n0) u_m <= u_(n0) / (rho - 1) + u_0 + u_1,
// which is the bound used when the NB head mass is not small.  Returns n0
// (0 for no skip), with the log NB head mass and the log bound on everything
// else the head holds.
inline int nb_window_start(double k, double q, double A, double& lhead, double& lbound,
                           double& ldrop) {
  lhead = NEG_INF; lbound = NEG_INF; ldrop = NEG_INF;
  if (A < 64.0) return 0;
  const double q1 = 1.0 - q;
  // peak: (1 - q)^2 A ((k + n)/(n + 1))^2 = alpha_n, by fixed point
  double n = 2.0 * q1 * q1 * A;
  for (int it = 0; it < 6; ++it) {
    const double r = (k + n) / (n + 1.0);
    n = std::max(0.0, 2.0 * q1 * q1 * A * r * r - k - 1.0);
  }
  if (n < 200.0) return 0;
  const double al = 0.5 * (k + n + 1.0);
  const double curv = k / ((k + n) * (n + 1.0)) + 0.25 / (al + 1.0);
  const double lpeak = nb_log_u(std::floor(n), k, q, A);
  double width = 9.0 / std::sqrt(curv) + 10.0;
  for (int attempt = 0; attempt < 3; ++attempt, width *= 2.0) {
    const double n0 = std::floor(n - width);
    if (n0 < 50.0) return 0;
    const double a0 = 0.5 * (k + n0 + 1.0);
    if (a0 > A - 1.0) continue;
    // NB head mass P(N <= n0 - 1); when w_n is log-concave (k >= 1) and still
    // rising at n0 it is at most w_(n0) / (rw - 1), and below the incomplete
    // beta's range it is only bounded (returned as ldrop, the S share dropped)
    double lHb = R_PosInf;
    if (k >= 1.0) {
      const double rw = q1 * (k + n0 - 1.0) / n0;
      if (rw > 1.0) lHb = R::dnbinom(n0, k, q, 1) - std::log(rw - 1.0);
    }
    const bool tiny_head = lHb < -650.0;
    // the NB mass at and above n0 is at most w_(n0) / (1 - up) once the ratio
    // q1 (k + m) / (m + 1) is below 1; when that is negligible the head is all of it
    const double up = std::max(q1, q1 * (k + n0) / (n0 + 1.0));
    const bool full_head = !tiny_head && up < 1.0 &&
      R::dnbinom(n0, k, q, 1) - std::log1p(-up) < -40.0;
    const double lH = tiny_head ? NEG_INF : (full_head ? 0.0 : R::pbeta(q, k, n0, 1, 1));
    const double lu0 = nb_log_u(n0, k, q, A), lrho = lu0 - nb_log_u(n0 - 1.0, k, q, A);
    double lsum = (tiny_head ? lHb : lH) + R::dgamma(A, a0 + 1.0, 1.0, 1);
    if (lrho > 0.0) {
      const double lgeo = lu0 - std::log(std::expm1(lrho));
      const double l01 = std::max(nb_log_u(0.0, k, q, A), nb_log_u(1.0, k, q, A)) + LOG2;
      const double hi = std::max(lgeo, l01);
      lsum = std::min(lsum, hi + std::log1p(std::exp(std::min(lgeo, l01) - hi)));
    }
    const double lb = lsum + std::log(a0 + 1.0 + upper_ratio(a0, A));
    if (lb < lpeak - 40.0) {
      lhead = lH; lbound = lb; ldrop = tiny_head ? lHb : NEG_INF;
      return static_cast<int>(n0);
    }
  }
  return 0;
}

inline void nb_branch(double t, double k, double q, double A, int want, Out& o,
                      bool force_f = false, bool allow_skip = true) {
  const double lA = std::log(A), lt = std::log(t);
  if (q >= 1.0 - 1e-15) {   // a single law: InvGamma((k + 1) / 2, A)
    const double al = 0.5 * (k + 1.0);
    o.logF = R::pgamma(A, al, 1.0, 0, 1);
    o.logS = R::pgamma(A, al, 1.0, 1, 1);
    o.logf = al * lA - A - lgam(al) - lt;
    o.terms = 1;
    return;
  }
  const double q1 = 1.0 - q, mean = k * q1 / q;
  // F directly when the bulk of the alpha_n sits below A
  const double ahat = 0.5 * (k + mean + 1.0);
  const bool f_first = force_f ||
    ((want & (WANT_CDF | WANT_SURV)) && A > ahat + 2.0 * std::sqrt(ahat) + 2.0);
  const bool want_S = !f_first && (want & (WANT_CDF | WANT_SURV));
  double lhead = NEG_INF, lbound = NEG_INF, ldrop = NEG_INF;
  const int n0 = allow_skip ? nb_window_start(k, q, A, lhead, lbound, ldrop) : 0;
  // Parity streams e = n & 1: u_n = w_n inc_n = wm im[e] exp(lw + lb[e]).
  // S: P(alpha) = inc(alpha) + P(alpha + 1) gives, with the same-parity
  // cumulative weights W_m = sum_(n0 <= n <= m) w_n (Wc[e] on wm's scale),
  //   sum_(n0 <= n <= N) w_n P_n = sum_(m <= N) inc_m W_m + sum_top (rho_top - 1) inc_top W_top,
  // rho = P / inc at the last index of each parity.
  double wm = 1.0;
  // log w_n0 and log inc(alpha) through R's saddle-point forms (no lgamma
  // differences, which lose digits at large orders)
  double lw = (n0 == 0) ? k * std::log(q) : R::dnbinom(n0, k, q, 1);
  double im[2], lb[2], al[2], Qm[2] = {0.0, 0.0}, Wc[2] = {0.0, 0.0}, tS[2] = {0.0, 0.0}, atop[2];
  for (int j = 0; j < 2; ++j) {
    const int e = (n0 + j) & 1;
    al[e] = 0.5 * (k + n0 + j + 1.0);
    lb[e] = log_inc(al[e], A, lA);
    im[e] = 1.0;
    if (f_first) Qm[e] = upper_ratio(al[e], A);   // Q / inc
  }
  double Lu = lw + std::max(lb[0], lb[1]), fac[2];
  double sU = 0.0, sf = 0.0, sF = 0.0, sS = 0.0;
  // fac[e] = exp(lw + lb[e] - Lu) <= 1: Lu, the scale of the sums, moves up
  // with the factors.  wm, im and Qm stay within MID^(+-1), so no product of
  // them overflows.
  auto rescale_sums = [&](double sc) {
    sU *= sc; sf *= sc; sF *= sc; sS *= sc; tS[0] *= sc; tS[1] *= sc;
  };
  auto refac = [&](int e) {
    const double l = lw + lb[e] - Lu;
    if (l > 0.0) {
      rescale_sums(std::exp(-l));
      Lu += l;
      fac[0] = std::exp(lw + lb[0] - Lu); fac[1] = std::exp(lw + lb[1] - Lu);
    } else {
      fac[e] = std::exp(l);
    }
  };
  refac(0); refac(1);
  const double stop_alpha = A + (f_first ? 2.0 * std::sqrt(A) + 8.0 : 1.0);
  int n = n0;
  for (;; ++n) {
    const int e = n & 1;
    const double u = wm * im[e] * fac[e];
    sU += u; sf += u * al[e];
    if (f_first) { sF += wm * Qm[e] * fac[e]; Qm[e] += im[e]; }
    if (want_S) { Wc[e] += wm; tS[e] = Wc[e] * im[e] * fac[e]; sS += tS[e]; atop[e] = al[e]; }
    im[e] *= A / (al[e] + 1.0);
    al[e] += 1.0;
    const double ratio = q1 * (k + n) / (n + 1.0);
    wm *= ratio;
    if (f_first ? (Qm[e] > MID || im[e] > MID) : im[e] > MID) {
      Qm[e] /= MID; im[e] /= MID; lb[e] += LMID; refac(e);
    } else if (im[e] < 1.0 / MID && (!f_first || Qm[e] < 1.0 / MID)) {
      Qm[e] *= MID; im[e] *= MID; lb[e] -= LMID; refac(e);
    }
    const double wbig = std::max(wm, std::max(Wc[0], Wc[1]));
    if (wbig > MID) {
      wm /= MID; Wc[0] /= MID; Wc[1] /= MID; lw += LMID; refac(0); refac(1);
    } else if (wbig < 1.0 / MID) {
      wm *= MID; Wc[0] *= MID; Wc[1] *= MID; lw -= LMID; refac(0); refac(1);
    }
    if (sU > BIG || sf > BIG || sF > BIG || sS > BIG) {
      rescale_sums(1.0 / BIG); Lu += LBIG; refac(0); refac(1);
    }
    if (!std::isfinite(sU) || !std::isfinite(sf) || !std::isfinite(sF) || !std::isfinite(sS) ||
        n - n0 >= NB_MAX_TERMS) {
      o.terms = -1;
      return;
    }
    const double amin = std::min(al[0], al[1]);
    // the bound below costs several divisions; screen with the current term
    if (amin > stop_alpha &&
        u * (amin + 2.0) <= 1e-14 * std::min(f_first ? sF : sU, sf)) {
      // u_(m+1)/u_m <= g for m > n (Wendel's bound on the half-step of inc)
      const double g = q1 * std::max(1.0, (k + n) / (n + 1.0)) *
        std::sqrt(A * (amin + 1.5)) / (amin + 1.0);
      if (g < 1.0) {
        const double c = 1.0 / (1.0 - A / (amin + 1.0));   // rho bound past n
        const double tail = u * g / (1.0 - g);
        bool ok = tail * c <= EPS * (f_first ? sF : sU) && tail * (amin + 2.0) / (1.0 - g) <= EPS * sf;
        if (ok && f_first) {
          // past here w_m Q_m = w_m - w_m P_m with the P part negligible, so the
          // rest of F is the NB mass beyond n: bounded geometrically, else
          // closed with the regularized incomplete beta I_(1-q)(n + 1, k)
          const double up = std::max(ratio, q1);
          const double lwt = std::log(wm) + lw - Lu - std::log1p(-up);
          if (lwt > std::log(EPS * sF)) {
            if (lwt + Lu > -650.0) {
              sF += std::exp(R::pbeta(q1, n + 1.0, k, 1, 1) - Lu);
            } else {
              // below the incomplete beta's range: sum the NB tail directly
              double wt = wm * std::exp(lw - Lu), tail_w = wt;   // w_(n+1)
              for (int m = n + 2; m < n + 2 + NB_MAX_TERMS; ++m) {
                wt *= q1 * (k + m - 1.0) / m;
                tail_w += wt;
                if (wt <= EPS * 1e-2 * (sF + tail_w)) break;
              }
              sF += tail_w;
            }
          }
        }
        if (ok) break;
      }
    }
  }
  // the skipped head must be negligible next to every sum it would feed
  auto head_ok = [&](double sum) { return n0 == 0 || lbound - Lu <= std::log(EPS * sum); };
  if (!head_ok(sf) || (f_first && !head_ok(sF))) {
    nb_branch(t, k, q, A, want, o, force_f, false);
    return;
  }
  o.terms = n - n0 + 1;
  o.logf = Lu + std::log(sf) - lt;
  if (!(want & (WANT_CDF | WANT_SURV))) return;
  if (f_first) {
    o.logF = std::min(0.0, Lu + std::log(sF));
    o.logS = log1m_exp(o.logF);
    return;
  }
  for (int e = 0; e < 2; ++e)
    if (tS[e] > 0.0) sS += (lower_ratio(atop[e], A) - 1.0) * tS[e];
  if (n0 > 0) {
    sS += std::exp(lhead - Lu);   // the head's NB mass, P = 1 - Q there
    if (!head_ok(sS) || ldrop - Lu > std::log(EPS * sS)) { nb_branch(t, k, q, A, want, o, force_f, false); return; }
  }
  o.logS = std::min(0.0, Lu + std::log(sS));
  o.logF = log1m_exp(o.logS);
  if ((want & WANT_CDF) && o.logF < std::log(TAIL_SMALL)) {
    Out o2;
    nb_branch(t, k, q, A, WANT_CDF, o2, true, allow_skip);
    o.logF = o2.logF;
    o.terms += o2.terms;
  }
}

// ------------------------------------------------------------ MOMENT branch
// S = (1/p) sum_{m>=2} (1 - r^(m-1)) D_(k+m).  D_(k+m) is log-concave in m with
// its peak near k + m = p delta and sd sqrt(p^2 + p delta); the tails beyond
// the summed window [mlo, mhi] are bounded geometrically from its end ratios.
// D obeys G's recurrence: D_nu = (nu + 1) D_(nu+2) / p^2 + z D_(nu+1) / p.
// With tail sums T_j = sum_(m>=j) D_(k+m) the weights telescope:
//   p S = T_mlo (1 - r^(mlo-2)) + (1 - r) r^(mlo-2) H,  H = sum_(j>=mlo) r^(j-mlo) T_j,
// so T and H accumulate in the same downward pass (H = T_j + r H).
constexpr int DENS_SWEEP_MAX = 48;   // walk on down to k for f, else seed G_k afresh

inline double moment_logS(double k, double p, double delta, double z, double q,
                          bool want_dens, double& lDk, int& terms) {
  // 1 - r and log r are formed from 2 / q directly: r rounds to 1 for large q
  const double omr = 2.0 / q, r = 1.0 - omr, lp = std::log(p);
  const double lr = (r > 0.0) ? std::log1p(-omr) : NAN_D;
  const double ip = 1.0 / p, ip2 = ip * ip, zp = z * ip;
  const double peak = p * delta - k, sd = std::sqrt(p * p + p * delta);
  double half = 8.0 * sd + 6.0;
  const double dconst = -0.5 * delta * delta - LOG_SQRT_2PI;
  for (int attempt = 0; attempt < 6; ++attempt, half *= 2.0) {
    const int mlo = std::max(2, static_cast<int>(std::floor(peak - half)));
    const int mhi = std::max(mlo + (6 << attempt), static_cast<int>(std::ceil(peak + half)));
    const double top = k + mhi;
    const int lift = seed_lift(top, z);
    const double ns = top + lift;
    const double l0 = ns * lp + logG_asym(ns, z), l1 = (ns + 1.0) * lp + logG_asym(ns + 1.0, z);
    double u0 = 1.0, u1 = std::exp(l1 - l0), ls = l0 + dconst, nu = ns;
    for (int i = 0; i < lift; ++i) {
      const double un = nu * ip2 * u1 + zp * u0;   // order nu - 1
      u1 = u0; u0 = un; nu -= 1.0;
      if (u0 > BIG) { u0 /= BIG; u1 /= BIG; ls += LBIG; }
    }
    // window, top down: u0 holds D_(k+mhi)
    const double d_top = u0, d_top1 = nu * ip2 * u1 + zp * u0;   // D_(k+mhi-1)
    double T = u0, H = u0;
    for (int m = mhi - 1; m >= mlo; --m) {
      const double un = nu * ip2 * u1 + zp * u0;
      u1 = u0; u0 = un; nu -= 1.0;
      T += u0;
      H = T + r * H;
      if (T > BIG) { u0 /= BIG; u1 /= BIG; T /= BIG; H /= BIG; ls += LBIG; }
    }
    const double d_lo = u0, d_lo1 = u1;   // D_(k+mlo), D_(k+mlo+1), same scale
    const double rp = std::pow(r, mlo - 2);
    const double head = (r > 0.0) ? -std::expm1((mlo - 2) * lr) : 1.0 - rp;
    const double sum = T * head + omr * rp * H;
    const double sc_top = std::exp(-(ls - (l0 + dconst)));  // top values' scale vs final
    const double g_hi = d_top / d_top1;
    const double g_lo = d_lo1 / d_lo;
    bool ok = g_hi < 1.0 && 2.0 * d_top * sc_top * g_hi / (1.0 - g_hi) <= EPS * sum;
    if (ok && mlo > 2) ok = g_lo > 1.0 && 2.0 * d_lo / (g_lo - 1.0) <= EPS * sum;
    if (!ok) continue;
    terms = mhi - mlo + 1 + lift;
    const double lS = std::log(sum) + ls - lp;
    if (want_dens) {
      if (mlo <= DENS_SWEEP_MAX) {
        for (int i = 0; i < mlo; ++i) {
          const double un = nu * ip2 * u1 + zp * u0;
          u1 = u0; u0 = un; nu -= 1.0;
          if (u0 < 1.0 / BIG) { u0 *= BIG; u1 *= BIG; ls -= LBIG; }
          else if (u0 > BIG) { u0 /= BIG; u1 /= BIG; ls += LBIG; }
        }
        lDk = ls + std::log(u0);
        terms += mlo;
      } else {
        lDk = k * lp + logG(k, z) + dconst;
      }
    }
    return lS;
  }
  terms = -1;
  return NAN_D;
}

// F = 2 p^k e^(-delta^2/2) / sqrt(2 pi) sum_{n even} W_n I_n(delta), W_n =
// J_(k+n)(z) / Gamma(k) (recurrence W_n = (W_(n+2) + z W_(n+1)) / (k + n)),
// I_n = G_(n+1)(delta) (I_n = (n + 2) I_(n+2) + delta I_(n+1)).
// The term ratio is about y0(k+n, z) / (y0(n+1, delta) + delta), which rises
// towards 1 like 1 - p / (2 sqrt(n)) (terms ~ e^(-p sqrt(n))), so the tail past
// a last ratio g is bounded by 2 / (1 - g)^2 times the last term rather than
// geometrically; the sum is ended where the saddle model has dropped that far
// below e^-40 of its maximum.
constexpr int F_MAX_TERMS = 20000;

inline int moment_F_length(double k, double z, double delta) {
  double prod = 1.0, lt = 0.0, lmax = 0.0, g_prev = 0.0, lt_prev = 0.0;
  for (int nn = 0; nn < F_MAX_TERMS; ++nn) {
    const double g = saddle(k + nn, z) / (saddle(nn + 1.0, delta) + delta);
    prod *= g;
    if ((nn & 7) == 7) {
      lt += std::log(prod); prod = 1.0;
      lmax = std::max(lmax, lt);
      if (g < 1.0 && lt < lmax - 41.0 + 2.0 * std::log1p(-g)) return nn + 8;
      // While the ratio rises no later block falls faster than this one, and
      // the stopping drop only grows: give up once that rate cannot reach it
      // within F_MAX_TERMS (small p, where the terms fall like e^(-p sqrt(n))).
      if (nn > 7 && g >= g_prev) {
        const double rate = lt_prev - lt;
        const double blocks = (F_MAX_TERMS - 1 - nn) / 8.0;
        const double need = lt - (lmax - 41.0 + (g < 1.0 ? 2.0 * std::log1p(-g) : 0.0));
        if (rate <= 0.0 || rate * blocks < need) return -1;
      }
      g_prev = g; lt_prev = lt;
    }
  }
  return -1;
}

// W and I both run downward, so the even terms are summed as they appear and
// the sum is rescaled with whichever factor is rescaled.  NaN when the series
// is longer than F_MAX_TERMS (the callers then use the survivor series).
inline double moment_logF(double k, double p, double delta, double z, double& lDk, int& terms) {
  terms = -1;
  int N = moment_F_length(k, z, delta);
  if (N < 0) return NAN_D;
  const double lp = std::log(p), lgk = lgam(k);
  for (int attempt = 0; attempt < 3 && N <= 2 * F_MAX_TERMS; ++attempt, N *= 2) {
    // W seeds at orders k + N, k + N + 1 (log W = log G + lgamma(order) - lgamma(k))
    const int liftW = seed_lift(k + N, z);
    const double nW = k + N + liftW;
    const double lw0 = logG_asym(nW, z) + lgam(nW) - lgk;
    double w1 = std::exp(logG_asym(nW + 1.0, z) + lgam(nW + 1.0) - lgk - lw0), w0 = 1.0, lsW = lw0, nuW = nW;
    for (int i = 0; i < liftW; ++i) {
      const double wn = (w1 + z * w0) / (nuW - 1.0);
      w1 = w0; w0 = wn; nuW -= 1.0;
      if (w0 > BIG) { w0 /= BIG; w1 /= BIG; lsW += LBIG; }
    }
    // I seeds at orders N + 1, N + 2 of G at delta (I_N, I_(N+1))
    const int liftI = seed_lift(N + 1.0, delta);
    const double nI = N + 1.0 + liftI;
    const double li0 = logG_asym(nI, delta);
    double v1 = std::exp(logG_asym(nI + 1.0, delta) - li0), v0 = 1.0, lsI = li0, nuI = nI;
    for (int i = 0; i < liftI; ++i) {
      const double vn = nuI * v1 + delta * v0;     // G_(nu-1) = nu G_(nu+1) + x G_nu
      v1 = v0; v0 = vn; nuI -= 1.0;
      if (v0 > BIG) { v0 /= BIG; v1 /= BIG; lsI += LBIG; }
    }
    // the top two even terms decide convergence; the top term is kept as a
    // log against lrs, the log of the factors the sum has been divided by
    // (it can fall below the double range relative to the final sum)
    double sum = 0.0, ltl = 0.0, g = NAN_D, lrs = 0.0;
    int seen = 0;
    auto add = [&](int n) {
      if (n & 1) return;
      const double tt = w0 * v0;
      sum += tt;
      if (seen == 0) ltl = std::log(tt) + lrs;
      else if (seen == 1) g = std::exp(ltl - lrs - std::log(tt));
      ++seen;
    };
    add(N);
    for (int i = N - 1; i >= 0; --i) {
      const double wn = (w1 + z * w0) / (k + i);
      w1 = w0; w0 = wn;
      const double vn = (i + 2.0) * v1 + delta * v0;
      v1 = v0; v0 = vn;
      // W rises with order and I falls, so either can leave the range
      if (w0 > BIG) { w0 /= BIG; w1 /= BIG; lsW += LBIG; sum /= BIG; lrs += LBIG; }
      else if (w0 < 1.0 / BIG) { w0 *= BIG; w1 *= BIG; lsW -= LBIG; sum *= BIG; lrs -= LBIG; }
      if (v0 > BIG) { v0 /= BIG; v1 /= BIG; lsI += LBIG; sum /= BIG; lrs += LBIG; }
      else if (v0 < 1.0 / BIG) { v0 *= BIG; v1 *= BIG; lsI -= LBIG; sum *= BIG; lrs -= LBIG; }
      add(i);
    }
    if (!(g < 1.0 && ltl - lrs + std::log(2.0 / ((1.0 - g) * (1.0 - g))) <= std::log(EPS * sum))) continue;
    terms = N + liftW + liftI;
    const double base = k * lp - 0.5 * delta * delta - LOG_SQRT_2PI;
    lDk = base + lsW + std::log(w0);
    return LOG2 + base + lsW + lsI + std::log(sum);
  }
  terms = -1;
  return NAN_D;
}

// ------------------------------------------------------------ NARROW branch
// E g(Y) = sum_n g^(n)(y0) mu_n / n!, mu_n the central moments of Gamma(k, p).
// With x = y0 - delta, w = y0 + delta, G(y) = phi(y - delta) R(y + delta):
// S_W' = -2 delta G, F_W = 1 - S_W; the Taylor coefficients of G are the
// Cauchy product of a_n (phi(x + h)) and b_n = (-1)^n I_n(w) (R(w + h)).
// The series behaves like a Poisson(lambda) sum over pairs of orders, lambda =
// (x sd(Y))^2 / 2, so its length is set from lambda (at most NT_MAX).
constexpr int NT_MAX = 320;
struct InvTable { double v[NT_MAX + 2]; constexpr InvTable() : v() { for (int i = 0; i < NT_MAX + 2; ++i) v[i] = 1.0 / (i + 1); } };
constexpr InvTable INV1{};   // INV1.v[n] = 1 / (n + 1)

inline bool narrow_branch(double t, double k, double p, double delta, int want, Out& o) {
  const double y0 = k / p, x = y0 - delta, w = y0 + delta;
  const double lam = 0.5 * x * x * k / (p * p);
  const double nt_want = 2.0 * (lam + 10.0 * std::sqrt(lam) + 20.0);
  if (!(nt_want <= NT_MAX)) return false;
  const int NT = std::max(72, static_cast<int>(nt_want));
  // I_n(w) = G_(n+1)(w), n = 0..NT+1, downward from seeds at order NT + 2;
  // brev[NT + 1 - j] = (-1)^j I_j so each Cauchy coefficient is a forward dot product
  alignas(32) double I[NT_MAX + 2], brev[NT_MAX + 2], a[NT_MAX + 2];
  {
    GWalk g = walk_from(NT + 2.0, w);     // g0 = G_(NT+2) = I_(NT+1)
    I[NT + 1] = g.g0;
    for (int n = NT; n >= 0; --n) {
      g.step(); I[n] = g.g0;
      if (g.g0 > BIG) {
        for (int j = n; j <= NT + 1; ++j) I[j] /= BIG;
        g.g0 /= BIG; g.g1 /= BIG; g.lscale += LBIG;
      }
    }
    const double s0 = std::exp(g.lscale);
    for (int n = 0; n <= NT + 1; ++n) { I[n] *= s0; brev[NT + 1 - n] = (n & 1) ? -I[n] : I[n]; }
    if (!(I[0] > 0.0) || !std::isfinite(I[0])) return false;
  }
  // every quantity below is proportional to phi(x); it is carried as 1 and
  // restored on the log scale, so |x| far out in the tail does not underflow
  const double lphx = -0.5 * x * x - LOG_SQRT_2PI;
  const double ip = 1.0 / p, kip2 = k * ip * ip, d2 = 2.0 * delta;
  a[0] = 1.0; a[1] = -x;
  double mu_m = 1.0, mu = 0.0;           // mu_(n-1), mu_n
  double dens = 1.0, C = 0.0;
  int nconv = -1, quiet = 0;
  for (int n = 1; n <= NT; ++n) {
    // mu_n from mu_(n-1), mu_(n-2): mu_(m+1) = m (mu_m / p + k mu_(m-1) / p^2)
    if (n >= 2) { const double mn = (n - 1) * (ip * mu + kip2 * mu_m); mu_m = mu; mu = mn; }
    if (n >= 2) a[n] = -(x * a[n - 1] + a[n - 2]) * INV1.v[n - 1];
    const double td = a[n] * mu;
    dens += td;
    double tc = 0.0;
    if (n >= 2) {
      const double* bb = brev + (NT + 2 - n);   // bb[i] = (-1)^(n-1-i) I_(n-1-i)
      double g = 0.0;
      for (int i = 0; i < n; ++i) g += a[i] * bb[i];
      tc = d2 * g * mu * INV1.v[n - 1];
      C += tc;
    }
    if (n >= 6 && std::fabs(td) <= EPS * std::fabs(dens) && std::fabs(tc) <= EPS * 1e-3 * std::max(1.0, std::fabs(C))) {
      if (++quiet >= 3) { nconv = n; break; }
    } else quiet = 0;
  }
  if (nconv < 0 || !(dens > 0.0)) return false;
  o.terms = nconv;
  o.logf = std::log(delta / t) + lphx + std::log(dens);
  // with R the Mills ratio: S_W(y0) / phi(x) = R(x) - I_0(w), F_W(y0) / phi(x) = R(-x) + I_0(w)
  auto mills = [](double y) { return std::exp(R::pnorm(-y, 0.0, 1.0, 1, 1) + 0.5 * y * y + LOG_SQRT_2PI); };
  if (x > 0.0) {
    double SW;
    if (4.0 * delta < x) {
      // R(x) - R(w) = sum_{n>=1} (2 delta)^n I_n(w)
      double acc = 0.0, h = 1.0;
      for (int n = 1; n <= NT + 1; ++n) { h *= d2; const double tt = h * I[n]; acc += tt; if (tt < EPS * acc) break; }
      SW = acc;
    } else {
      SW = mills(x) - I[0];
    }
    const double Ss = SW - C;
    if (!(Ss > 0.0)) return false;
    o.logS = lphx + std::log(Ss);
    if (!(o.logS < 0.0)) return false;
    o.logF = log1m_exp(o.logS);
  } else {
    const double Fs = mills(-x) + I[0] + C;
    if (!(Fs > 0.0)) return false;
    o.logF = lphx + std::log(Fs);
    if (!(o.logF < 0.0)) return false;
    o.logS = log1m_exp(o.logF);
  }
  (void)want;
  return true;
}

// ------------------------------------------------------------ WIDE branch
// With U ~ Exp(2 delta) and Z ~ N(0, 1) independent, S_W(y) = P(delta - U + Z > y)
// for y >= 0 (an exponentially modified Gaussian kernel), so with X = U - Z
//   S = P(Y + X < delta) = E G_Y(delta - X),   F = E Q_Y(delta - X),
//   f = (delta / t) E g_Y(delta + Z),
// G_Y, Q_Y, g_Y the CDF, survivor and density of Y.  When g_Y varies slowly on
// the unit scale around delta, Taylor series about delta against the moments of
// X and Z converge quickly:
//   S = G_Y(delta) + g_Y(delta) sum_(n>=1) (-1)^n d_(n-1) m_n,  F = Q_Y(delta) - (same sum),
//   f = (delta / t) g_Y(delta) sum_(n even) d_n E Z^n / n!,
// d_n = g_Y^(n)(delta) / g_Y(delta) from y g' = (k - 1 - p y) g:
//   d_(n+1) = [(k - 1 - n - p delta) d_n - p n d_(n-1)] / delta,
// and m_n = E X^n / n! from the cumulants of X (kappa_1 = 1/l, kappa_2 = 1/l^2 + 1,
// kappa_j = (j-1)!/l^j, l = 2 delta): m_n = (1/n) sum_(j<=n) e_j m_(n-j), e_j = l^-j + [j = 2].
// The expansions are asymptotic only through the branch point of g_Y at y = 0,
// at distance delta, whose share is below e^(-delta^2/2).
constexpr double WIDE_DELTA = 15.0, WIDE_SLOPE = 1.5, WIDE_CURV = 0.05;
constexpr int WIDE_MAX = 96;

// Tilted form for a steep g_Y at delta.  With psi1 = (k - 1)/delta - p the
// log-slope of g_Y there, g_Y(delta + h) = g_Y(delta) e^(psi1 h) phi(h),
// phi(h) = (1 + h/delta)^(k-1) e^(-(k-1) h/delta) = sum_j c_j h^j (c_1 = 0,
// c_j from the exponential of the log1p series).  Tilting A = -X = Z - U by
// e^(psi1 A) (Z -> N(psi1, 1), U -> Exp(2 delta + psi1), mass
// M = e^(psi1^2/2) 2 delta / (2 delta + psi1)) gives, for the tail on the side
// psi1 points away from,
//   S = g_Y(delta) M sum_j c_j j! b_j  (psi1 > 0),  F = -g_Y(delta) M sum_j c_j j! b_j  (psi1 < 0),
//   b_j = (a_j - b_(j-1)) / psi1,  a_j = E_tilt A^j / j!,
//   f = (delta / t) g_Y(delta) e^(psi1^2/2) sum_j c_j j! E (Z + psi1)^j / j!.
// The c_j are of order ((k-1)/(2 delta^2))^(j/2) and the moments of order
// (|psi1| + 1/|psi1| + 1)^j, so the series is short when their product is.
// The tail series is asymptotic through the branch point at y = 0: its terms
// bottom out near j ~ |psi1| delta at about e^(-|psi1| delta), so that product
// must be large.
constexpr double TILT_GATE = 2.0, TILT_REACH = 80.0;

inline bool tilt_branch(double t, double k, double p, double delta, double slope, int want, Out& o) {
  const double eps = std::fabs(k - 1.0) / (2.0 * delta * delta);
  const double reach = std::fabs(slope) + 1.0 / std::fabs(slope) + 1.0;
  if (!(eps * reach * reach <= TILT_GATE && 4.0 * (reach + 2.0) <= delta && 2.0 * delta + slope > 0.0 &&
        std::fabs(slope) * delta >= TILT_REACH))
    return false;
  const double lg = R::dgamma(delta, k, 1.0 / p, 1);
  if (!std::isfinite(lg)) return false;
  const double lam = 2.0 * delta + slope;
  const double lM = 0.5 * slope * slope + std::log(2.0 * delta / lam);
  // c_j of phi: log phi = sum_(m>=2) l_m h^m, l_m = (k-1) (-1)^(m+1) / (m delta^m)
  double c[WIDE_MAX + 1], l[WIDE_MAX + 1], a[WIDE_MAX + 1], e[WIDE_MAX + 1], z[WIDE_MAX + 1];
  c[0] = 1.0; l[0] = l[1] = 0.0;
  { double dm = 1.0 / delta;
    for (int m = 2; m <= WIDE_MAX; ++m) { dm /= delta; l[m] = (k - 1.0) * ((m & 1) ? 1.0 : -1.0) * dm / m; } }
  // cumulants of A under the tilt: kappa_1 = psi1 - 1/lam, kappa_2 = 1 + 1/lam^2,
  // kappa_j = (-1)^j (j-1)! / lam^j; e_j = kappa_j / (j-1)!
  { double il = 1.0;
    for (int j = 1; j <= WIDE_MAX; ++j) { il /= lam; e[j] = ((j & 1) ? -il : il); }
    e[1] += slope; e[2] += 1.0; }
  a[0] = 1.0; z[0] = 1.0;
  double b_prev = 0.0, jf = 1.0;              // b_(j-1), j!
  double ssum = 0.0, dsum = 0.0;
  int quiet_s = 0, quiet_d = 0, j = 0;
  bool s_done = !(want & (WANT_CDF | WANT_SURV)), d_done = !(want & WANT_PDF);
  for (; j <= WIDE_MAX && !(s_done && d_done); ++j) {
    if (j >= 1) {
      double cj = 0.0;
      for (int m = 2; m <= j; ++m) cj += m * l[m] * c[j - m];
      c[j] = cj / j;
      jf *= j;
      double aj = 0.0;
      for (int i = 1; i <= j; ++i) aj += e[i] * a[j - i];
      a[j] = aj / j;
      z[j] = (slope * z[j - 1] + (j >= 2 ? z[j - 2] : 0.0)) / j;
    }
    const double cjf = c[j] * jf;
    if (!s_done) {
      const double bj = (a[j] - b_prev) / slope;
      b_prev = bj;
      const double ts = cjf * bj;
      ssum += ts;
      if (j >= 4 && std::fabs(ts) <= EPS * std::fabs(ssum)) { if (++quiet_s >= 3) s_done = true; }
      else quiet_s = 0;
    }
    if (!d_done) {
      const double td = cjf * z[j];
      dsum += td;
      if (j >= 4 && std::fabs(td) <= EPS * std::fabs(dsum)) { if (++quiet_d >= 3) d_done = true; }
      else quiet_d = 0;
    }
    if (!std::isfinite(ssum) || !std::isfinite(dsum) || !std::isfinite(jf)) return false;
  }
  if (!(s_done && d_done)) return false;
  o.terms = j;
  if (want & WANT_PDF) {
    if (!(dsum > 0.0)) return false;
    o.logf = std::log(delta / t) + lg + 0.5 * slope * slope + std::log(dsum);
  }
  if (want & (WANT_CDF | WANT_SURV)) {
    const double tail = (slope > 0.0) ? ssum : -ssum;
    if (!(tail > 0.0)) return false;
    const double lt = lg + lM + std::log(tail);
    if (!(lt < 0.0)) return false;
    if (slope > 0.0) { o.logS = lt; o.logF = log1m_exp(lt); }
    else             { o.logF = lt; o.logS = log1m_exp(lt); }
  }
  return true;
}

inline bool wide_branch(double t, double k, double p, double delta, int want, Out& o) {
  const double slope = (k - 1.0) / delta - p;
  if (delta >= WIDE_DELTA && std::fabs(slope) > 1.0 && tilt_branch(t, k, p, delta, slope, want, o))
    return true;
  if (!(delta >= WIDE_DELTA && std::fabs(slope) <= WIDE_SLOPE &&
        std::fabs(k - 1.0) <= WIDE_CURV * delta * delta)) return false;
  const double lg = R::dgamma(delta, k, 1.0 / p, 1);
  if (!std::isfinite(lg)) return false;
  const double il = 0.5 / delta, kd = k - 1.0 - p * delta;
  double d_prev = 0.0, d = 1.0;        // d_(n-1), d_n
  double m[WIDE_MAX + 1], e[WIDE_MAX + 1];
  m[0] = 1.0;
  double el = 1.0;
  for (int j = 1; j <= WIDE_MAX; ++j) { el *= il; e[j] = el + (j == 2 ? 1.0 : 0.0); }
  double dens = 1.0, hz = 1.0;         // sum_(n even) d_n E Z^n / n!, E Z^n / n! for even n
  double cs = 0.0;                     // sum_(n>=1) (-1)^n d_(n-1) m_n
  int quiet_d = 0, quiet_c = 0, n = 1;
  bool dens_done = !(want & WANT_PDF), cdf_done = !(want & (WANT_CDF | WANT_SURV));
  for (; n <= WIDE_MAX && !(dens_done && cdf_done); ++n) {
    // d holds d_(n-1) here
    if (!cdf_done) {
      double mn = 0.0;
      for (int j = 1; j <= n; ++j) mn += e[j] * m[n - j];
      m[n] = mn / n;
      const double tc = ((n & 1) ? -d : d) * m[n];
      cs += tc;
      if (n >= 4 && std::fabs(tc) <= EPS * std::max(1.0, std::fabs(cs))) { if (++quiet_c >= 3) cdf_done = true; }
      else quiet_c = 0;
    }
    const double dn = (kd - (n - 1)) * d / delta - p * (n - 1) * d_prev / delta;   // d_n
    d_prev = d; d = dn;
    if (!dens_done && !(n & 1)) {
      hz *= 1.0 / (2.0 * (n / 2));     // E Z^n / n! = 1 / (2^(n/2) (n/2)!)
      const double td = d * hz;
      dens += td;
      if (n >= 4 && std::fabs(td) <= EPS * std::fabs(dens)) { if (++quiet_d >= 2) dens_done = true; }
      else quiet_d = 0;
    }
    if (!std::isfinite(d) || !std::isfinite(cs)) return false;
  }
  if (!(dens_done && cdf_done)) return false;
  o.terms = n;
  if (want & WANT_PDF) {
    if (!(dens > 0.0)) return false;
    o.logf = std::log(delta / t) + lg + std::log(dens);
  }
  if (want & (WANT_CDF | WANT_SURV)) {
    // S = e^lG + e^lg cs,  F = e^lQ - e^lg cs, each formed on the larger scale
    auto comb = [&](double la, double sgn) {
      const double v = sgn * cs;
      if (la >= lg) { const double r = 1.0 + v * std::exp(lg - la); return r > 0.0 ? la + std::log(r) : NAN_D; }
      const double r = std::exp(la - lg) + v;
      return r > 0.0 ? lg + std::log(r) : NAN_D;
    };
    // the smaller tail directly, the other as its complement
    const double lG = R::pgamma(delta, k, 1.0 / p, 1, 1), lQ = R::pgamma(delta, k, 1.0 / p, 0, 1);
    if (lG <= lQ) {
      const double lS = comb(lG, 1.0);
      if (!(lS < 0.0)) return false;
      o.logS = lS; o.logF = log1m_exp(lS);
    } else {
      const double lF = comb(lQ, -1.0);
      if (!(lF < 0.0)) return false;
      o.logF = lF; o.logS = log1m_exp(lF);
    }
  }
  return true;
}

// Far-tail screen.  F_W(y) <= 2 Phi(y - delta) and S_W(y) <= Phi(delta - y),
// and phi(y - delta) <= phi(0), so splitting the drift law at any c between
// delta and y0, with h = |c - delta|, bounds the small tail and the density by
//   F <= 2 [P(Y > c) + Phi(-h)]   (y0 < delta),   S <= P(Y < c) + Phi(-h)   (y0 > delta),
//   f <= (delta / t) [P(tail) phi(0) + phi(h)].
// These bounds are returned in place of the values when the small tail is
// below e^LOG_NEGLIGIBLE (the other tail is then 1 to double precision), when
// only the other tail is wanted and the small one is below e^LOG_UNIT (1 - it
// rounds to 1), or when only the density is wanted and its bound is below
// lf_floor (a caller's floor).
constexpr double LOG_NEGLIGIBLE = -1000.0, LOG_UNIT = -40.0;

inline bool far_tail(double t, double k, double p, double delta, int want, double lf_floor, Out& o) {
  const double y0 = k / p, sdY = std::sqrt(k) / p, gap = std::fabs(y0 - delta);
  const bool early = y0 < delta;
  const double lscale = std::log(delta / t) - LOG_SQRT_2PI;
  // the largest log value of the small tail (or of f / lscale) that may be bounded
  double thr = LOG_NEGLIGIBLE;
  if (want == (early ? WANT_SURV : WANT_CDF)) thr = LOG_UNIT;
  if (want == WANT_PDF) thr = std::max(thr, lf_floor - lscale);
  // the split c is free: the midpoint, or where a normal drift law's tail
  // meets the Wald one (best when the drift law is narrow)
  const double w_lo = sdY / (1.0 + sdY);   // c = delta + (1 - w)(y0 - delta)
  const double hmax = std::max(0.5, 1.0 - w_lo) * gap;
  if (0.5 * hmax * hmax < -thr) return false;
  double ltail = R_PosInf, lfb = R_PosInf;
  for (const double wt : {0.5, w_lo}) {
    const double h = (1.0 - wt) * gap, c = delta + (1.0 - wt) * (y0 - delta);
    const double lP = R::pgamma(c, k, 1.0 / p, early ? 0 : 1, 1);
    const double lPhi = R::pnorm(-h, 0.0, 1.0, 1, 1), lphi = -0.5 * h * h;
    const double hi = std::max(lP, lPhi), fhi = std::max(lP, lphi);
    ltail = std::min(ltail, (early ? LOG2 : 0.0) + hi + std::log1p(std::exp(std::min(lP, lPhi) - hi)));
    lfb = std::min(lfb, fhi + std::log1p(std::exp(std::min(lP, lphi) - fhi)));
    if (wt == w_lo) break;
  }
  const bool ok = ltail < LOG_NEGLIGIBLE ||
    (want == WANT_PDF ? lfb < thr : ltail < thr);
  if (!ok) return false;
  o.logf = lscale + lfb;
  if (early) { o.logF = ltail; o.logS = 0.0; }
  else       { o.logS = ltail; o.logF = 0.0; }
  o.terms = 0;
  return true;
}

// -------------------------------------------------------------------- driver
constexpr double NARROW_P = 16.0, NARROW_SD = 0.3;
constexpr double DENS_DIRECT_KMAX = 2000.0;

// force: 0 automatic, else the branch to use (NB needs q <= 1, MOMENT q > 1;
// NB_NOSKIP is NB summed from n = 0).
// lf_floor: with want = WANT_PDF, a log density known to be below it may be
// returned as an upper bound (callers that floor the log likelihood).
inline Out gamma_wald(double t, double d, double s, double k, double tau, int want = WANT_ALL,
                      int force = 0, double lf_floor = NEG_INF) {
  Out o;
  if (!(t > 0.0)) { o.logf = NEG_INF; o.logF = NEG_INF; o.logS = 0.0; return o; }
  const double rt = std::sqrt(t);
  const double p = tau * s / rt, delta = d / (s * rt), z = p - delta, q = p / delta;
  const double sdY = std::sqrt(k) / p, y0 = k / p;
  if (force == 0 && far_tail(t, k, p, delta, want, lf_floor, o)) { o.branch = FAR_TAIL; return o; }
  if (want == WANT_PDF && force == 0 && q > 1.0 && k <= DENS_DIRECT_KMAX) {
    // f = (delta / t) p^k e^(-delta^2/2) G_k(z) / sqrt(2 pi), z > 0; k log p
    // and log G_k cancel to O(k) absolute error, so large k (the sv -> 0
    // limit) goes through the branches below
    o.logf = std::log(delta / t) + k * std::log(p) - 0.5 * delta * delta - LOG_SQRT_2PI + logG(k, z);
    o.branch = MOMENT;
    return o;
  }
  if (force == WIDE || force == 0) {
    if (wide_branch(t, k, p, delta, want, o)) { o.branch = WIDE; return o; }
    if (force == WIDE) return o;
  }
  if (force == NARROW || (force == 0 && p >= NARROW_P && sdY <= NARROW_SD && std::fabs(y0 - delta) <= p / 3.0)) {
    if (narrow_branch(t, k, p, delta, want, o)) { o.branch = NARROW; return o; }
    if (force == NARROW) return o;
  }
  if (force == NB || force == NB_NOSKIP || (force == 0 && q <= 1.0)) {
    nb_branch(t, k, q, 0.5 * delta * delta, want, o, false, force != NB_NOSKIP);
    o.branch = NB;
    return o;
  }
  o.branch = MOMENT;
  const double ldt = std::log(delta / t);
  // F series when F is small: the gamma bulk sits well below delta (F >= 2 Phi(-delta)
  // at y = 0) and so does all but a small share of its (for small k, long)
  // upper tail; otherwise the S series
  const bool f_first = (y0 + 2.0 * sdY) < delta - 1.0 && delta > 2.6 &&
    R::pgamma(delta, k, 1.0 / p, 0, 1) < std::log(TAIL_SMALL);
  int t1 = 0, t2 = 0;
  double lDk = NAN_D, ld2;
  bool done = false;
  if (f_first) {
    const double lF = moment_logF(k, p, delta, z, lDk, t1);
    if (std::isfinite(lF)) {
      o.logF = std::min(0.0, lF);
      o.logf = ldt + lDk;
      o.logS = log1m_exp(lF);
      if ((want & WANT_SURV) && o.logS < std::log(TAIL_SMALL)) {
        const double lS = moment_logS(k, p, delta, z, q, false, ld2, t2);
        if (std::isfinite(lS)) o.logS = std::min(0.0, lS);
      }
      done = true;
    }
  }
  if (!done) {
    o.logS = moment_logS(k, p, delta, z, q, want & WANT_PDF, lDk, t1);
    if (o.logS > 0.0) o.logS = 0.0;   // NaN (no window closed) is kept
    o.logf = ldt + lDk;
    o.logF = log1m_exp(o.logS);
    if ((want & WANT_CDF) && o.logF < std::log(TAIL_SMALL)) {
      const double lF = moment_logF(k, p, delta, z, ld2, t2);
      if (std::isfinite(lF)) o.logF = std::min(0.0, lF);
    }
  }
  if ((want & WANT_PDF) && !std::isfinite(o.logf))
    o.logf = ldt + k * std::log(p) - 0.5 * delta * delta - LOG_SQRT_2PI + logG(k, z);
  o.terms = t1 + t2;
  return o;
}

}  // namespace rdmswtn_gamma
#endif
