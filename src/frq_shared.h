#ifndef EMC2_FRQ_SHARED_H
#define EMC2_FRQ_SHARED_H

// ============================================================================
// FRQ with a shared trial rate (FRQ(rate_sharing = "shared")).
//
// In counter coordinates accumulator i finishes by time x iff
// Y_i <= W_i phi_i(x), with Y_i ~ Gamma(K_i) its threshold epoch,
// phi_i = q_i / (1 - q_i) its operational clock and W_i ~ Gamma(beta, 1) its
// trial input-rate multiplier (src/model_FRQ.h).  The independent race draws
// one W_i per accumulator; the shared race draws one W per trial for all of
// them.  Every marginal is unchanged, which requires one beta = 1/cv^2 per
// trial; availability, clock, acceleration, t0, K and plo stay per
// accumulator because they only enter through phi_i.
//
// CONJUGACY.  Given W = w the accumulators are independent counters, and
// w^{beta-1} e^{-w} is conjugate to every Poisson-counter factor.  For winner
// w with threshold K_w and losers l,
//
//   f_w(x) = f_w^ind(x) E[ prod_l Q(K_l, G psi_l) ],
//   G ~ Gamma(beta + K_w, 1),   psi_l = phi_l(x) (1 - q_w(x)),
//
// where f_w^ind is the independent FRQ density and Q the regularised upper
// incomplete gamma.  The tilted rate is Gamma(beta + K_w, 1 + phi_w), and
// rescaling it to unit rate gives psi_l.  The joint survivor is the same
// expectation at G ~ Gamma(beta, 1) with psi_i = phi_i over all accumulators.
//
// One factor is a beta tail:  E[Q(K, G psi)] = I_{1/(1+psi)}(b, K).  So the
// two-accumulator density is closed form, the independent density times
// I_{1-z}(beta + K_w, K_l) in place of the independent survivor
// I_{1-q_l}(beta, K_l).
//
// Two factors (two losers, or the joint survivor of two accumulators) give a
// one-dimensional series: expanding P(K_2, y) = sum_n e^{-y} y^{K_2+n} /
// Gamma(K_2 + n + 1) tilts G once more per term, so
//
//   E[Q_1 Q_2] = E[Q_1] - sum_n w_n J_n,
//   w_n = NB(n; b, r) on K_2 + n,  r = psi_2 / (1 + psi_2),
//   J_n = I_x(b + K_2 + n, K_1),   x = 1 / (1 + psi_1 / (1 + psi_2)).
//
// w_n and J_n both have O(1) recurrences, sum_n w_n = I_r(K_2, b) is known,
// and J_n decreases in n, so (remaining weight) x J_N bounds the truncation
// error.  The J_n recurrence carries an absolute error of order eps J_0, so a
// result that does not clear it by a wide margin, and three or more factors,
// go to a quadrature over log G, whose
// integrand is concave for K >= 1 (a log-concave survivor composed with an
// exponential), so a mode-centred Gauss-Legendre rule converges fast.
//
// LOW STATE.  Each accumulator draws its threshold independently given W, so
// each loser factor is the mixture (1 - plo) Q(K, y) + plo Q(K_lo, y), and
// the expectation expands over the 2^m low/main assignments.  At K_lo = 1,
// Q(1, y) = e^{-y} is a further exponential tilt: it multiplies by
// (1 + s)^{-b} and divides the remaining psi by 1 + s, s the summed low psi,
// so a low loser removes a factor rather than adding one.
//
// OMEGA.  The excess criterion J = Y_{K-1} + (1 + omega) Y_1 is a
// Gamma(K + N) epoch with N ~ Geometric(1 / (1 + omega)), so under sharing it
// is a geometric series over the closed forms above (one-dimensional for the
// winner's own excess, two-dimensional when a loser has one too).  It is not
// implemented here; R/model_FRQ.R refuses omega with a shared rate.
//
// cv -> 0 (the Poisson branch) makes W degenerate, so the race factorises
// and the independent kernels are used.
// ============================================================================

#include "model_FRQ.h"
#include "gl_quad.h"
#include <algorithm>
#include <array>
#include <cmath>

// Series terms before the quadrature takes over.
constexpr int FRQSH_SERIES_MAX = 20000;
// Relative error allowed in a series result before it is recomputed by
// quadrature; the subtraction E[Q_1] - sum loses digits as the result falls.
constexpr double FRQSH_CANCEL_REL = 1e-8;
// Factors per expectation (2^m low/main assignments), i.e. racing
// accumulators per trial.
constexpr int FRQSH_MAX_FACTORS = 8;

// One loser factor: its main threshold, log psi (-Inf for an unstarted
// accumulator, whose factor is one) and its low state.
struct FrqShFactor {
  double K;
  double log_psi;
  const FrqLow* lo;
};

// log E[Q(K, G psi)], G ~ Gamma(b, 1), = log I_{1/(1+psi)}(b, K).
inline double frqsh_log_e1(double b, double K, double log_psi) {
  if (log_psi == R_NegInf) return 0.0;
  if (log_psi == R_PosInf) return R_NegInf;
  const double l1p = frq_log1pexp(log_psi);
  const double log_x = -l1p, log_omx = log_psi - l1p;
  return frq_log_ibeta(std::exp(log_x), std::exp(log_omx), log_x, log_omx,
                       b, K, lbeta_fast(b, K));
}

// log of the log-G integrand: Gamma(b) density of s = log G times the factors.
inline double frqsh_log_integrand(double s, double b, double lgb, int m,
                                  const double* K, const double* lgK,
                                  const double* log_psi) {
  const double G = std::exp(s);
  double out = b * s - G - lgb;
  for (int l = 0; l < m; ++l) {
    const double log_y = s + log_psi[l];
    double lq;
    gamma_log_pq_fast(K[l], std::exp(log_y), log_y, lgK[l], nullptr, &lq);
    out += lq;
  }
  return out;
}

// d/ds of the log integrand: b - G - sum_l y_l p_K(y_l) / Q(K, y_l).
inline double frqsh_log_integrand_ds(double s, double b, int m,
                                     const double* K, const double* lgK,
                                     const double* log_psi) {
  double out = b - std::exp(s);
  for (int l = 0; l < m; ++l) {
    const double log_y = s + log_psi[l];
    const double y = std::exp(log_y);
    double lq;
    gamma_log_pq_fast(K[l], y, log_y, lgK[l], nullptr, &lq);
    out -= std::exp(K[l] * log_y - y - lgK[l] - lq);
  }
  return out;
}

// log E[prod_l Q(K_l, G psi_l)] by quadrature over s = log G.  The rule is
// centred on the integrand's mode and extends until the log integrand has
// fallen 40 nats on each side.
inline double frqsh_log_quad(double b, int m, const double* K,
                             const double* log_psi) {
  std::array<double, FRQSH_MAX_FACTORS> lgK{};
  for (int l = 0; l < m; ++l) lgK[l] = lgamma_fast(K[l]);
  const double lgb = lgamma_fast(b);
  auto L = [&](double s) {
    return frqsh_log_integrand(s, b, lgb, m, K, lgK.data(), log_psi);
  };
  auto D = [&](double s) {
    return frqsh_log_integrand_ds(s, b, m, K, lgK.data(), log_psi);
  };
  // The derivative is decreasing and D(log b) <= 0; step down to a sign change.
  double hi = std::log(b), lo = hi - 1.0, step = 1.0;
  int guard = 0;
  while (D(lo) < 0.0 && guard++ < 200) { hi = lo; step *= 2.0; lo = hi - step; }
  // Illinois regula falsi on D.
  double flo = D(lo), fhi = D(hi);
  double mode = hi;
  if (flo > 0.0 && fhi < 0.0) {
    int side = 0;
    for (int it = 0; it < 60; ++it) {
      mode = (lo * fhi - hi * flo) / (fhi - flo);
      const double fm = D(mode);
      if (std::fabs(fm) < 1e-10 || hi - lo < 1e-10) break;
      if (fm > 0.0) {
        lo = mode; flo = fm;
        if (side == -1) fhi *= 0.5;
        side = -1;
      } else {
        hi = mode; fhi = fm;
        if (side == 1) flo *= 0.5;
        side = 1;
      }
    }
  }
  const double Lm = L(mode);
  if (!R_FINITE(Lm)) return R_NegInf;
  const double d = 1e-3;
  const double curv = (L(mode + d) - 2.0 * Lm + L(mode - d)) / (d * d);
  const double sd = (curv < 0.0 && R_FINITE(curv)) ? 1.0 / std::sqrt(-curv) : 1.0;
  auto reach = [&](double dir) {
    double step = sd;
    for (int it = 0; it < 60; ++it) {
      const double s = mode + dir * step;
      if (!(Lm - L(s) < 40.0)) return s;
      step *= 1.6;
    }
    return mode + dir * step;
  };
  const double a = reach(-1.0), c = reach(1.0);
  const GLRule& rule = gl_get_rule(32);
  double max_log = R_NegInf;
  std::array<double, 64> vals{};
  int k = 0;
  for (int panel = 0; panel < 2; ++panel) {
    const double p0 = panel == 0 ? a : mode, p1 = panel == 0 ? mode : c;
    const double half = 0.5 * (p1 - p0), mid = 0.5 * (p1 + p0);
    for (int i = 0; i < 32; ++i, ++k) {
      vals[k] = std::log(half * rule.w[i]) + L(mid + half * rule.x[i]);
      if (vals[k] > max_log) max_log = vals[k];
    }
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < 64; ++i) sum += std::exp(vals[i] - max_log);
  return max_log + std::log(sum);
}

// log E[Q(K1, G psi1) Q(K2, G psi2)] by the negative-binomial series in
// loser 2; false when it should be recomputed by quadrature.
inline bool frqsh_log_e2_series(double b, double K1, double log_psi1,
                                double K2, double log_psi2, double* out) {
  const double le1 = frqsh_log_e1(b, K1, log_psi1);
  if (!(le1 > -700.0)) return false;
  const double E1 = std::exp(le1);
  const double l1p2 = frq_log1pexp(log_psi2);
  const double log_r = log_psi2 - l1p2, log_omr = -l1p2;
  const double a = b + K2;
  const double lw0 = lgamma_fast(a) - lgamma_fast(b) - lgamma_fast(K2 + 1.0) +
    K2 * log_r + b * log_omr;
  if (!(lw0 > -700.0)) return false;
  const double log_tot = frq_log_ibeta(std::exp(log_r), std::exp(log_omr),
                                       log_r, log_omr, K2, b, lbeta_fast(K2, b));
  const double tot = std::exp(log_tot);
  const double log_psi1p = log_psi1 - l1p2;
  const double l1p1 = frq_log1pexp(log_psi1p);
  const double log_x = -l1p1, log_omx = log_psi1p - l1p1;
  const double x = std::exp(log_x);
  double J = std::exp(frq_log_ibeta(x, std::exp(log_omx), log_x, log_omx,
                                    a, K1, lbeta_fast(a, K1)));
  const double J0 = J;
  double T = std::exp(a * log_x + K1 * log_omx - std::log(a) -
                      lbeta_fast(a, K1));
  const double r = std::exp(log_r);
  double w = std::exp(lw0), sum = 0.0, cum = 0.0;
  const double tol = 1e-16 * E1;
  int n = 0;
  for (;; ++n) {
    sum += w * J;
    cum += w;
    if ((tot - cum) * J <= tol || J <= 0.0) break;
    if (n >= FRQSH_SERIES_MAX) return false;
    w *= r * (a + n) / (K2 + n + 1.0);
    J -= T;
    if (J < 0.0) J = 0.0;
    T *= x * (a + n + K1) / (a + n + 1.0);
  }
  const double val = E1 - sum;
  // The J recurrence carries an absolute error of order eps J_0 per step, and
  // J_0 can be far larger than E[Q_1]; the difference must clear that error.
  const double err = 2.2e-16 * (n + 10.0) * (J0 * tot + E1);
  if (!(val > err / FRQSH_CANCEL_REL)) return false;
  *out = std::log(val);
  return true;
}

// log E[prod_l Q(K_l, G psi_l)] for main-threshold factors only.
inline double frqsh_log_e_main(double b, int m, const double* K,
                               const double* log_psi) {
  if (m == 0) return 0.0;
  if (m == 1) return frqsh_log_e1(b, K[0], log_psi[0]);
  if (m == 2) {
    // Expand the factor with the smaller psi: fewer terms, smaller subtraction.
    const int i2 = log_psi[1] <= log_psi[0] ? 1 : 0, i1 = 1 - i2;
    double out;
    if (frqsh_log_e2_series(b, K[i1], log_psi[i1], K[i2], log_psi[i2], &out))
      return out;
  }
  return frqsh_log_quad(b, m, K, log_psi);
}

// log E[prod_l ((1 - plo_l) Q(K_l, G psi_l) + plo_l Q(K_lo, G psi_l))],
// G ~ Gamma(b, 1).
inline double frqsh_log_expect(double b, int m, const FrqShFactor* f) {
  std::array<int, FRQSH_MAX_FACTORS> idx{};
  int n = 0;
  for (int l = 0; l < m; ++l) {
    if (f[l].log_psi == R_NegInf) continue;   // unstarted: factor one
    if (f[l].log_psi == R_PosInf) return R_NegInf;
    if (n >= FRQSH_MAX_FACTORS)
      Rcpp::stop("FRQ shared rate: at most %d racing accumulators are supported.",
                 FRQSH_MAX_FACTORS);
    idx[n++] = l;
  }
  std::array<int, FRQSH_MAX_FACTORS> mix{};
  int n_mix = 0;
  for (int j = 0; j < n; ++j) if (f[idx[j]].lo->active) mix[n_mix++] = j;
  double total = R_NegInf;
  std::array<double, FRQSH_MAX_FACTORS> K{}, lp{};
  std::array<unsigned char, FRQSH_MAX_FACTORS> low{};
  for (int mask = 0; mask < (1 << n_mix); ++mask) {
    low.fill(0);
    double lw = 0.0;
    for (int j = 0; j < n_mix; ++j) {
      const FrqLow& L = *f[idx[mix[j]]].lo;
      if (mask & (1 << j)) { low[mix[j]] = 1; lw += L.log_plo; }
      else lw += L.log_1m_plo;
    }
    if (lw == R_NegInf) continue;
    // Low factors at K_lo = 1 are exponential tilts.
    double log_s = R_NegInf;
    for (int j = 0; j < n; ++j) {
      if (low[j] && f[idx[j]].lo->K_lo == 1.0)
        log_s = frq_logaddexp(log_s, f[idx[j]].log_psi);
    }
    const double l1ps = (log_s == R_NegInf) ? 0.0 : frq_log1pexp(log_s);
    int c = 0;
    for (int j = 0; j < n; ++j) {
      const FrqShFactor& fj = f[idx[j]];
      if (low[j] && fj.lo->K_lo == 1.0) continue;
      K[c] = low[j] ? fj.lo->K_lo : fj.K;
      lp[c] = fj.log_psi - l1ps;
      ++c;
    }
    const double term = lw - b * l1ps + frqsh_log_e_main(b, c, K.data(), lp.data());
    total = frq_logaddexp(total, term);
  }
  return total;
}

// One accumulator at one time point: its availability state on the
// operational clock, with the density prefactor log(p g) (plus the
// acceleration Jacobian a u) and log phi.
struct FrqShPoint {
  bool started = false;
  double log_q = R_NegInf;
  double log_omq = 0.0;
  double log_pg = R_NegInf;
  double log_phi = R_NegInf;
};

inline FrqShPoint frqsh_point(double u, const FrqPars& s, double accel) {
  FrqShPoint pt;
  if (!(u > 0.0)) return pt;
  pt.started = true;
  const double x = frq_accel_warp(u, accel);
  const FrqState z = frq_state(x, s);
  pt.log_q = z.log_q;
  pt.log_omq = z.log_omq;
  pt.log_pg = s.log_p + z.log_g + ((accel > 0.0 && R_FINITE(u)) ? accel * u : 0.0);
  pt.log_phi = (z.log_omq == R_NegInf) ? R_PosInf : z.log_q - z.log_omq;
  return pt;
}

// The independent FRQ density at threshold K, from its point: p g q^{K-1}
// (1 - q)^{beta-1} / B(K, beta).
inline double frqsh_log_dens_k(const FrqShPoint& pt, const FrqPars& s,
                               double K, double lbeta_k) {
  double lp = pt.log_pg - lbeta_k;
  if (K != 1.0) lp += (K - 1.0) * pt.log_q;
  if (s.beta_minus_1 != 0.0) lp += s.beta_minus_1 * pt.log_omq;
  return lp;
}

// Joint log density that accumulator w finishes first at its point, with the
// other started accumulators (sl, pl) still running.  All rows share beta.
inline double frqsh_log_cause(const FrqPars& sw, const FrqShPoint& pw, int m,
                              const FrqPars* const* sl,
                              const FrqShPoint* pl) {
  if (!pw.started || pw.log_pg == R_NegInf) return R_NegInf;
  std::array<FrqShFactor, FRQSH_MAX_FACTORS> f{};
  if (m >= FRQSH_MAX_FACTORS)
    Rcpp::stop("FRQ shared rate: at most %d racing accumulators are supported.",
               FRQSH_MAX_FACTORS);
  for (int l = 0; l < m; ++l) {
    f[l].K = sl[l]->alpha;
    f[l].lo = &sl[l]->lo;
    f[l].log_psi = pl[l].started ? pl[l].log_phi + pw.log_omq : R_NegInf;
  }
  const double beta = sw.beta;
  double out = R_NegInf;
  if (sw.lo.log_1m_plo > R_NegInf) {
    out = sw.lo.log_1m_plo + frqsh_log_dens_k(pw, sw, sw.alpha, sw.lbeta_val) +
      frqsh_log_expect(beta + sw.alpha, m, f.data());
  }
  if (sw.lo.active) {
    out = frq_logaddexp(out, sw.lo.log_plo +
      frqsh_log_dens_k(pw, sw, sw.lo.K_lo, sw.lo.lbeta_lo) +
      frqsh_log_expect(beta + sw.lo.K_lo, m, f.data()));
  }
  return out;
}

// Joint log survivor of the n accumulators at their points.
inline double frqsh_log_surv(int n, const FrqPars* const* s,
                             const FrqShPoint* p) {
  std::array<FrqShFactor, FRQSH_MAX_FACTORS> f{};
  if (n > FRQSH_MAX_FACTORS)
    Rcpp::stop("FRQ shared rate: at most %d racing accumulators are supported.",
               FRQSH_MAX_FACTORS);
  double beta = 0.0;
  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (!p[i].started) continue;
    beta = s[i]->beta;
    f[m].K = s[i]->alpha;
    f[m].lo = &s[i]->lo;
    f[m].log_psi = p[i].log_phi;
    ++m;
  }
  if (m == 0) return 0.0;
  return frqsh_log_expect(beta, m, f.data());
}

#endif  // EMC2_FRQ_SHARED_H
