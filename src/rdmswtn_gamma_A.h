#ifndef EMC2_RDMSWTN_GAMMA_A_H
#define EMC2_RDMSWTN_GAMMA_A_H

// Gamma-drift Wald with a uniform start point: distance D ~ U[B, B + A],
// drift V ~ Gamma(k, rate tau).  log f, log F and log S averaged over D in
// closed form, without quadrature.  Derivation: Math/RDMSWTN_Derivation.tex,
// subsection "Start-Point Variability" of "Gamma Drift Variability".
//
// Notation as in rdmswtn_gamma.h (p = tau s / sqrt(t), Y ~ Gamma(k, p)), with
// distances u in units of s sqrt(t) running over [uL, uH], Delta = uH - uL,
// W = Y + eps (eps ~ N(0, 1)) and, at each u,
//   D = E phi(Y - u)  (the density of W),  C = P(W <= u),  T = P(W > u),
//   Sf = E (u - W)_+,  Pf = E (W - u)_+,   F, S the point-start CDF and survivor.
// Then, with z = p - u and y0 = k / p,
//   t Delta f_A = E[W 1{uL < W <= uH}],  E[W 1{W <= u}] = u C - Sf,  E[W 1{W > u}] = u T + Pf,
//   Delta F_A = Pi(uL) - Pi(uH),   Pi = Pf + (D + (p/2 - u) F - z T) / (1 - k),
//   Delta S_A = Sg(uH) - Sg(uL),   Sg = Sf + (z C + (u - p/2) S + D) / (1 - k).
// Pi(u) = E (M - u)_+ for the running maximum M of the drifted process; the
// 1 / (1 - k) terms are the moments E[g(Y) / Y] after integrating by parts in
// Y.  The singularity at k = 1 is removable; there the result is interpolated
// between k = 1 -+ K_BAND.
//
// Gamma (+) normal pieces (C, T, Sf, Pf), on the side of the median of W where
// they are the small tails:
//   lower  C = (1/p) sum_{m>=1} D_(k+m),  Sf = (1/p^2) sum_{m>=2} (m - 1) D_(k+m);
//   upper  T = D sum_n w_n I_n,  Pf = D sum_n (n + 1) w_n I_(n+1),
//          w_n = J_(k+n)(z) / J_k(z),  I_n = G_(n+1)(u);
//   narrow Taylor series about y0 against the gamma central moments.
// The D_(k+m) and w_n recurrences run downward from asymptotic seeds when
// z >= 0 and upward from the point-kernel densities at shapes k and k + 1
// when z < 0; in both directions every term is a positive combination.

#include "rdmswtn_gamma.h"
#include <initializer_list>
#include <vector>

namespace rdmswtn_gamma {

constexpr double K_BAND = 2e-3;        // |1 - k| below this: interpolate in k
constexpr double SHORT_SPAN4 = 0.05;   // spans (Delta times the local scale) below
constexpr double SHORT_SPAN8 = 0.5;    // these use 4 or 8 Gauss-Legendre nodes
constexpr double PIMG_BAND = 0.1;      // |1 - k| below this: image term also from the series
constexpr int PIMG_MAX_TERMS = 2000;   // when it is at most this long
constexpr double LOG_T_SMALL = -6.907755278982137;   // log 1e-3
constexpr double T_SMALL_U = 3.090232306167813;      // Phi(-u) >= 1e-3 below this
constexpr double LOG_T_LADDER = -11.512925464970229;  // log 1e-5: T = 1 - C keeps ~11 digits above this
constexpr int FALLBACK_PANELS = 128;  // composite Gauss-Legendre when the closed form fails
constexpr double FALLBACK_TOL = 1e-11; // two panel counts agree to this (relative, in logs)
constexpr double COND_MAX = 1e5;       // closed form kept while it loses < 5 digits
constexpr int CONV_MAX_TERMS = 200000;
// log(1 - 1e-4): a tail formed as one minus this keeps all but four digits
constexpr double LOG_COMPLEMENT_OK = -1.0000500033334732e-4;

enum APath { A_POINT = 1, A_SHORT = 2, A_CLOSED = 3, A_BAND = 4, A_FALLBACK = 5, A_FAR = 6 };

struct AOut {
  double logf = NAN_D, logF = NAN_D, logS = NAN_D;
  int path = 0;
};

// Gamma (+) normal and point values at one distance, as logs.  On the lower
// side lC and lSf are exact (lT, lPf follow from complements); on the upper
// side lT and lPf.  lF, lS and lD come from the point kernel.
struct Endpoint {
  double u = 0.0, lD = NAN_D, lD1 = NAN_D, lF = NAN_D, lS = NAN_D;
  double lC = NAN_D, lSf = NAN_D, lT = NAN_D, lPf = NAN_D;
  double lPimg = NAN_D;   // log E[K / (2Y)] from the upper series, when it ran
  bool lower = true, both = false;
};

inline double log_mills(double x) {   // log R(x) = log Phi(-x) - log phi(x)
  return R::pnorm(-x, 0.0, 1.0, 1, 1) + 0.5 * x * x + LOG_SQRT_2PI;
}

// ------------------------------------------------------------ lower ladder
// A0 = sum_{m>=1} D_(k+m), A1 = sum_{m>=2} (m - 1) D_(k+m); returns false
// when the window does not close.
inline bool ladder_lower(double k, double p, double u, double lD0, double lD1,
                         double& lA0, double& lA1) {
  const double z = p - u, lp = std::log(p), ip = 1.0 / p, ip2 = ip * ip, zp = z * ip;
  if (z >= 0.0) {
    // D_(k+m) is log-concave in m with its peak near k + m = p u
    const double peak = p * u - k, sd = std::sqrt(p * p + p * u);
    double half = 8.0 * sd + 6.0;
    const double dconst = -0.5 * u * u - LOG_SQRT_2PI;
    for (int attempt = 0; attempt < 6; ++attempt, half *= 2.0) {
      const int mlo = std::max(1, static_cast<int>(std::floor(peak - half)));
      const int mhi = std::max(mlo + (8 << attempt), static_cast<int>(std::ceil(peak + half)));
      const double top = k + mhi;
      const int lift = seed_lift(top, z);
      const double ns = top + lift;
      const double l0 = ns * lp + logG_asym(ns, z), l1 = (ns + 1.0) * lp + logG_asym(ns + 1.0, z);
      double u0 = 1.0, u1 = std::exp(l1 - l0), ls = l0 + dconst, nu = ns;
      for (int i = 0; i < lift; ++i) {
        const double un = nu * ip2 * u1 + zp * u0;
        u1 = u0; u0 = un; nu -= 1.0;
        if (u0 > BIG) { u0 /= BIG; u1 /= BIG; ls += LBIG; }
      }
      const double ls_top = ls;
      const double d_top = u0, d_top1 = nu * ip2 * u1 + zp * u0;
      double a0 = u0, a1 = (mhi - 1.0) * u0;
      for (int m = mhi - 1; m >= mlo; --m) {
        const double un = nu * ip2 * u1 + zp * u0;
        u1 = u0; u0 = un; nu -= 1.0;
        a0 += u0; a1 += (m - 1.0) * u0;
        if (a0 > BIG) { u0 /= BIG; u1 /= BIG; a0 /= BIG; a1 /= BIG; ls += LBIG; }
      }
      const double sc_top = std::exp(ls_top - ls);
      const double g_hi = d_top / d_top1;
      if (!(g_hi < 1.0)) continue;
      const double tail_hi = d_top * sc_top * g_hi / (1.0 - g_hi);
      const double tail_hi1 = tail_hi * (mhi + 1.0 / (1.0 - g_hi));
      if (!(tail_hi <= EPS * a0 && tail_hi1 <= EPS * std::max(a1, a0))) continue;
      if (mlo > 1) {
        const double g_lo = u1 / u0;   // D_(k+mlo+1) / D_(k+mlo)
        if (!(g_lo > 1.0 && 2.0 * u0 / (g_lo - 1.0) * mlo <= EPS * std::max(a1, a0) &&
              2.0 * u0 / (g_lo - 1.0) <= EPS * a0)) continue;
      }
      if (!(a0 > 0.0)) return false;
      lA0 = std::log(a0) + ls;
      lA1 = a1 > 0.0 ? std::log(a1) + ls : NEG_INF;
      return true;
    }
    return false;
  }
  // z < 0: upward from D_k, D_(k+1); the terms fall once k + m passes p u
  const double p2 = p * p, pz = p * z;
  double d0 = 1.0, d1 = std::exp(lD1 - lD0), ls = lD0;
  double a0 = d1, a1 = 0.0;
  for (int m = 0; m < CONV_MAX_TERMS; ++m) {
    // d0 = D_(k+m), d1 = D_(k+m+1)  ->  D_(k+m+2)
    const double d2 = (p2 * d0 - pz * d1) / (k + m + 1.0);
    d0 = d1; d1 = d2;
    const double mm = m + 2.0;    // order of d1 above k
    a0 += d1; a1 += (mm - 1.0) * d1;
    if (d1 < 1.0 / BIG && a0 < 1.0 / BIG) {
      d0 *= BIG; d1 *= BIG; a0 *= BIG; a1 *= BIG; ls -= LBIG;
    }
    const double g = d1 / d0;
    if (k + mm > p * u && g < 1.0) {
      const double tail = d1 * g / (1.0 - g);
      if (tail <= EPS * a0 && tail * (mm + 1.0 / (1.0 - g)) <= EPS * std::max(a1, a0)) {
        if (!(a0 > 0.0)) return false;
        lA0 = std::log(a0) + ls;
        lA1 = a1 > 0.0 ? std::log(a1) + ls : NEG_INF;
        return true;
      }
    }
  }
  return false;
}

// ------------------------------------------------------------ upper series
// moment_F_length's scan with the term ratio's log, which varies slowly in n,
// taken at the middle of each block of 8 past the first LEN_EXACT terms; the
// margin keeps the length at or above the exact scan's.
constexpr int LEN_EXACT = 64, LEN_MARGIN = 16;
inline int upper_length(double k, double z, double u, int max_terms = F_MAX_TERMS) {
  double prod = 1.0, lt = 0.0, lmax = 0.0, g_prev = 0.0, lt_prev = 0.0;
  for (int nn = 7; nn < max_terms; nn += 8) {
    double g;
    if (nn < LEN_EXACT) {
      for (int m = nn - 7; m <= nn; ++m) { g = saddle(k + m, z) / (saddle(m + 1.0, u) + u); prod *= g; }
      lt += std::log(prod); prod = 1.0;
    } else {
      const double mid = nn - 3.5;
      lt += 8.0 * std::log(saddle(k + mid, z) / (saddle(mid + 1.0, u) + u));
      g = saddle(k + nn, z) / (saddle(nn + 1.0, u) + u);
    }
    lmax = std::max(lmax, lt);
    if (g < 1.0 && lt < lmax - 41.0 + 2.0 * std::log1p(-g)) return nn + 1 + LEN_MARGIN;
    if (nn > 7 && g >= g_prev) {
      const double rate = lt_prev - lt;
      const double blocks = (max_terms - 1 - nn) / 8.0;
      const double need = lt - (lmax - 41.0 + (g < 1.0 ? 2.0 * std::log1p(-g) : 0.0));
      if (rate <= 0.0 || rate * blocks < need) return -1;
    }
    g_prev = g; lt_prev = lt;
  }
  return -1;
}

// B0 = sum_n w_n I_n, B1 = sum_n (n + 1) w_n I_(n+1), B2 = sum_{n even} w_n I_(n+1),
// w_n = J_(k+n)(z) / J_k(z),
// I_n = G_(n+1)(u) (I_n = (n + 2) I_(n+2) + u I_(n+1), downward).
inline bool series_upper(double k, double p, double u, double lD0, double lD1,
                         double& lB0, double& lB1, double& lB2) {
  const double z = p - u;
  int N = upper_length(k, z, u);
  if (N < 0) return false;
  thread_local std::vector<double> Ibuf;
  for (int attempt = 0; attempt < 3 && N <= 2 * F_MAX_TERMS; ++attempt, N *= 2) {
    double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, last0 = NAN_D, last1 = NAN_D;
    double ls = 0.0;
    if (z >= 0.0) {
      // both downward: w from seeds at orders k + N, k + N + 1, I from N + 1, N + 2
      const int liftW = seed_lift(k + N, z);
      const double nW = k + N + liftW;
      // only ratios to w_0 are needed, so w carries no absolute scale; a
      // rescale of w rescales the sums with it, a rescale of I is recorded in lsI
      const double lw0 = logG_asym(nW, z) + lgam(nW);
      double w1 = std::exp(logG_asym(nW + 1.0, z) + lgam(nW + 1.0) - lw0), w0 = 1.0, nuW = nW;
      for (int i = 0; i < liftW; ++i) {
        const double wn = (w1 + z * w0) / (nuW - 1.0);
        w1 = w0; w0 = wn; nuW -= 1.0;
        if (w0 > BIG) { w0 /= BIG; w1 /= BIG; }
      }
      // w0 = J_(k+N), w1 = J_(k+N+1); I from G_(N+2) = I_(N+1)
      GWalk g = walk_from(N + 2.0, u);
      double i1 = g.g0;
      g.step();
      double i0 = g.g0;                   // I_N
      double lsI = g.lscale;
      sum0 = w0 * i0; sum1 = (N + 1.0) * w0 * i1; sum2 = (N & 1) ? 0.0 : w0 * i1;
      last0 = sum0; last1 = sum1;
      double prev0 = NAN_D, prev1 = NAN_D;
      for (int n = N - 1; n >= 0; --n) {
        const double wn = (w1 + z * w0) / (k + n);
        w1 = w0; w0 = wn;
        const double in = (n + 2.0) * i1 + u * i0;
        i1 = i0; i0 = in;
        const double t0 = w0 * i0, t1 = (n + 1.0) * w0 * i1;
        if (n == N - 1) { prev0 = t0; prev1 = t1; }
        sum0 += t0; sum1 += t1;
        if (!(n & 1)) sum2 += w0 * i1;
        // ls: log of the factor the sums have been divided by since the top term
        if (w0 > BIG) { w0 /= BIG; w1 /= BIG; sum0 /= BIG; sum1 /= BIG; sum2 /= BIG; ls += LBIG; }
        else if (w0 < 1.0 / BIG) { w0 *= BIG; w1 *= BIG; sum0 *= BIG; sum1 *= BIG; sum2 *= BIG; ls -= LBIG; }
        if (i0 > BIG) { i0 /= BIG; i1 /= BIG; sum0 /= BIG; sum1 /= BIG; sum2 /= BIG; ls += LBIG; lsI += LBIG; }
        else if (i0 < 1.0 / BIG) { i0 *= BIG; i1 *= BIG; sum0 *= BIG; sum1 *= BIG; sum2 *= BIG; ls -= LBIG; lsI -= LBIG; }
      }
      // tail beyond N: last ratio g, remainder <= 2 / (1 - g)^2 of the last term
      const double gr0 = last0 / prev0, gr1 = last1 / prev1;
      if (!(gr0 < 1.0 && gr1 < 1.0)) continue;
      const double lt0 = std::log(last0) - ls, lt1 = std::log(last1) - ls;
      if (!(lt0 + std::log(2.0 / ((1.0 - gr0) * (1.0 - gr0))) <= std::log(EPS * sum0))) continue;
      if (!(lt1 + std::log(2.0 / ((1.0 - gr1) * (1.0 - gr1))) <= std::log(EPS * sum1))) continue;
      const double lw = std::log(w0);   // J_k on the sums' w scale
      lB0 = std::log(sum0) - lw + lsI;
      lB1 = std::log(sum1) - lw + lsI;
      lB2 = std::log(sum2) - lw + lsI;
      return true;
    }
    // z < 0: the ratios I_(n+1) / I_n from a downward walk, then upward with
    // the w ratios rho_n = w_n / w_(n-1) = (k + n - 2) / rho_(n-1) - z,
    // rho_1 = (k / p) D_(k+1) / D_k; term_n = w_n I_n / I_0 = term_(n-1) rho_n rI_(n-1)
    Ibuf.resize(N + 1);
    double lI0;
    {
      GWalk g = walk_from(N + 2.0, u);    // g0 = G_(N+2) = I_(N+1)
      for (int n = N; n >= 0; --n) {
        g.step();                         // g0 = G_(n+1) = I_n, g1 = I_(n+1)
        Ibuf[n] = g.g1 / g.g0;
        g.rescale();
      }
      lI0 = g.log_g0();
    }
    double rho = (k / p) * std::exp(lD1 - lD0), term = 1.0, lsum = 0.0;
    sum0 = 1.0; sum1 = Ibuf[0]; sum2 = Ibuf[0];
    double pterm = 1.0;
    for (int n = 1; n <= N; ++n) {
      if (n >= 2) rho = (k + n - 2.0) / rho - z;
      pterm = term;
      term *= rho * Ibuf[n - 1];
      sum0 += term; sum1 += (n + 1.0) * term * Ibuf[n];
      if (!(n & 1)) sum2 += term * Ibuf[n];
      if (term > BIG) { term /= BIG; pterm /= BIG; sum0 /= BIG; sum1 /= BIG; sum2 /= BIG; lsum += LBIG; }
    }
    const double gr = term / pterm;
    if (!(gr < 1.0)) continue;
    const double bnd = term * 2.0 / ((1.0 - gr) * (1.0 - gr));
    if (!(bnd <= EPS * sum0 && bnd * (N + 1.0 / (1.0 - gr)) <= EPS * sum1)) continue;
    lB0 = std::log(sum0) + lsum + lI0;
    lB1 = std::log(sum1) + lsum + lI0;
    lB2 = std::log(sum2) + lsum + lI0;
    return true;
  }
  return false;
}

// ------------------------------------------------------------ narrow series
// E h(Y) = sum_n h^(n)(y0) mu_n / n! with x0 = y0 - u and phi(x0 + h) / phi(x0)
// = sum a_n h^n; in units of phi(x0):
//   D = sum a_n mu_n,  T = R(-x0) + sC,  C = R(x0) - sC,  sC = sum_{n>=2} a_(n-1) mu_n / n,
//   Pf = 1 + x0 R(-x0) + sS,  Sf = 1 - x0 R(x0) + sS,  sS = sum_{n>=2} a_(n-2) mu_n / (n (n-1)),
// with 1 - |x| R(|x|) = G_2(|x|).
inline bool narrow_conv(double k, double p, double u, Endpoint& e) {
  const double y0 = k / p, x = y0 - u;
  const double lam = 0.5 * x * x * k / (p * p);
  const double nt_want = 2.0 * (lam + 10.0 * std::sqrt(lam) + 20.0);
  if (!(nt_want <= NT_MAX)) return false;
  const int NT = std::max(72, static_cast<int>(nt_want));
  const double ip = 1.0 / p, kip2 = k * ip * ip;
  double am2 = 0.0, am1 = 1.0, an = -x;   // a_(n-2), a_(n-1), a_n at n = 1
  double mu_m = 1.0, mu = 0.0;            // mu_(n-1), mu_n
  double sD = 1.0, sC = 0.0, sS = 0.0;
  int quiet = 0;
  bool conv = false;
  for (int n = 1; n <= NT; ++n) {
    if (n >= 2) {
      const double mn = (n - 1) * (ip * mu + kip2 * mu_m); mu_m = mu; mu = mn;
      am2 = am1; am1 = an; an = -(x * am1 + am2) * INV1.v[n - 1];
    }
    const double td = an * mu;
    sD += td;
    double tc = 0.0, ts = 0.0;
    if (n >= 2) {
      tc = am1 * mu / n; sC += tc;
      ts = am2 * mu / (n * (n - 1.0)); sS += ts;
    }
    if (n >= 6 && std::fabs(td) <= EPS * std::fabs(sD) && std::fabs(tc) <= EPS * 1e-3 &&
        std::fabs(ts) <= EPS * 1e-3) {
      if (++quiet >= 3) { conv = true; break; }
    } else quiet = 0;
  }
  if (!conv || !(sD > 0.0)) return false;
  const double lphx = -0.5 * x * x - LOG_SQRT_2PI;
  const double ax = std::fabs(x);
  const double lRa = log_mills(ax), lG2 = logG(2.0, ax);
  // the small tails take the series correction relative to their own leading term
  double Csm, Sfsm, Tsm, Pfsm;
  if (x >= 0.0) {
    Csm = std::exp(lRa) - sC; Sfsm = std::exp(lG2) + sS;
    Tsm = std::exp(log_mills(-x)) + sC; Pfsm = 1.0 + x * std::exp(log_mills(-x)) + sS;
  } else {
    Tsm = std::exp(lRa) + sC; Pfsm = std::exp(lG2) + sS;
    Csm = std::exp(log_mills(x)) - sC; Sfsm = 1.0 - x * std::exp(log_mills(x)) + sS;
  }
  if (!(Csm > 0.0 && Tsm > 0.0 && Sfsm > 0.0 && Pfsm > 0.0)) return false;
  e.lD = lphx + std::log(sD);
  e.lC = lphx + std::log(Csm); e.lT = lphx + std::log(Tsm);
  e.lSf = lphx + std::log(Sfsm); e.lPf = lphx + std::log(Pfsm);
  e.both = true;
  return true;
}

// The expansions about u leave out the mass of Y near 0, which adds at most
// Phi(-u) <= phi(u) / u to T and phi(u) / u^2 to Pf; they are kept only while
// that is below e^LOG_BRANCH_OK of the results (it is not when k is tiny).
constexpr double LOG_BRANCH_OK = -29.933606208922594;   // log 1e-13
inline bool branch_negligible(double u, double lT, double lPf) {
  const double lphi = -0.5 * u * u - LOG_SQRT_2PI, lu = std::log(u);
  return lphi - lu < lT + LOG_BRANCH_OK && lphi - 2.0 * lu < lPf + LOG_BRANCH_OK;
}

// ------------------------------------------------------------ wide series
// Broad drift law at a large distance (as the point kernel's WIDE branch): with
// g_Y slowly varying on the unit scale around u, expand Q_Y(u - eps) and
// psi_Y(u - eps) = E (Y - u + eps)_+ about u against the normal moments,
//   T = Q_Y(u) - g_Y(u) sum_{n even>=2} d_(n-1) / (2^(n/2) (n/2)!),
//   Pf = psi_Y(u) + g_Y(u) sum_{n even>=2} d_(n-2) / (2^(n/2) (n/2)!),
// d_n = g_Y^(n)(u) / g_Y(u) by the recurrence in rdmswtn_gamma.h, and
// psi_Y(u) = (Q_Y(u) / p) [k Q(k + 1, p u) / Q(k, p u) - p u].
// Without the Wald's exponential part the branch point's share is e^(-u^2/2),
// so the series holds from a smaller distance than in the point kernel; it is
// ended relative to the sums, as the asymptotic floor allows.
constexpr double WIDE_CONV_U = 7.0, WIDE_CONV_SLOPE = 3.0, WIDE_CONV_TOL = 1e-14;

inline bool wide_conv(double k, double p, double u, Endpoint& e) {
  const double slope = (k - 1.0) / u - p;
  if (!(u >= WIDE_CONV_U && std::fabs(slope) <= WIDE_CONV_SLOPE && std::fabs(k - 1.0) <= WIDE_CURV * u * u))
    return false;
  const double lg = R::dgamma(u, k, 1.0 / p, 1), lQ = R::pgamma(u, k, 1.0 / p, 0, 1);
  const double lQ1 = R::pgamma(u, k + 1.0, 1.0 / p, 0, 1);
  if (!std::isfinite(lg) || !std::isfinite(lQ) || !std::isfinite(lQ1)) return false;
  const double x = p * u;
  const double psi_rel = (k * std::exp(lQ1 - lQ) - x) / p;   // psi_Y / Q_Y
  if (!(psi_rel > 0.0)) return false;
  // d_n = g^(n) / g: d_(n+1) = [(k - 1 - n - p u) d_n - p n d_(n-1)] / u
  double d[WIDE_MAX + 2];
  d[0] = 1.0; d[1] = slope;
  double cn = 1.0, sT = 0.0, sP = 0.0;
  const double qg = std::exp(lQ - lg), pg = qg * psi_rel;
  int quiet = 0;
  bool conv = false;
  for (int n = 2; n <= WIDE_MAX; n += 2) {
    d[n] = ((k - 1.0 - (n - 1) - x) * d[n - 1] - p * (n - 1) * d[n - 2]) / u;
    d[n + 1] = ((k - 1.0 - n - x) * d[n] - p * n * d[n - 1]) / u;
    cn /= n;                               // 1 / (2^(n/2) (n/2)!) = prod 1 / (2 j)
    const double tT = cn * d[n - 1], tP = cn * d[n - 2];
    sT += tT; sP += tP;
    if (n >= 6 && std::fabs(tT) <= WIDE_CONV_TOL * std::fabs(qg - sT) &&
        std::fabs(tP) <= WIDE_CONV_TOL * std::fabs(pg + sP)) {
      if (++quiet >= 2) { conv = true; break; }
    } else quiet = 0;
    if (!std::isfinite(sT) || !std::isfinite(sP)) return false;
  }
  if (!conv) return false;
  // in units of g: T / g = Q / g - sT, Pf / g = psi / g + sP
  const double Tg = qg - sT, Pg = pg + sP;
  if (!(Tg > 0.0 && Pg > 0.0)) return false;
  const double lT = lg + std::log(Tg), lPf = lg + std::log(Pg);
  if (!branch_negligible(u, lT, lPf)) return false;
  e.lT = lT;
  e.lPf = lPf;
  return true;
}

// ------------------------------------------------------------ tilted series
// Steep g_Y at u (as the point kernel's tilt_branch, without the exponential
// part): with psi = (k - 1)/u - p < 0, g_Y(u + h) = g_Y(u) e^(psi h) phi(h),
// phi(h) = sum_j c_j h^j, and
//   T = int g_Y(u + h) Phi(h) dh,  Pf = int g_Y(u + h) Psi(h) dh,  Psi(h) = h Phi(h) + phi(h),
// so with z_j = E (Z + psi)^j / j!, beta_j = (z_j + beta_(j-1)) / (-psi) and
// gamma_j = (beta_j + gamma_(j-1)) / (-psi),
//   T = g_Y(u) e^(psi^2/2) sum_j c_j j! beta_j,  Pf = g_Y(u) e^(psi^2/2) sum_j c_j j! gamma_j.
// The series is asymptotic: its terms bottom out near j ~ -psi u at about
// e^(psi u) of the sums.
constexpr double TILT_CONV_REACH = 30.0, TILT_CONV_TOL = 1e-13;

inline bool tilt_conv(double k, double p, double u, Endpoint& e) {
  const double psi = (k - 1.0) / u - p;
  if (!(psi < 0.0)) return false;
  const double eps = std::fabs(k - 1.0) / (2.0 * u * u);
  const double reach = -psi - 1.0 / psi + 1.0;
  if (!(eps * reach * reach <= TILT_GATE && 2.0 * (reach + 2.0) <= u && -psi * u >= TILT_CONV_REACH))
    return false;
  const double lg = R::dgamma(u, k, 1.0 / p, 1);
  if (!std::isfinite(lg)) return false;
  // c_j of phi: log phi = sum_(m>=2) l_m h^m, l_m = (k-1) (-1)^(m+1) / (m u^m)
  double c[WIDE_MAX + 1], l[WIDE_MAX + 1];
  c[0] = 1.0; l[0] = l[1] = 0.0;
  { double um = 1.0 / u;
    for (int m = 2; m <= WIDE_MAX; ++m) { um /= u; l[m] = (k - 1.0) * ((m & 1) ? 1.0 : -1.0) * um / m; } }
  const double ia = -1.0 / psi;
  double z2 = 0.0, z1 = 1.0, beta = 0.0, gam = 0.0, jf = 1.0, sT = 0.0, sP = 0.0;
  int quiet = 0;
  bool conv = false;
  for (int j = 0; j <= WIDE_MAX; ++j) {
    double zj = 1.0;
    if (j >= 1) {
      double cj = 0.0;
      for (int m = 2; m <= j; ++m) cj += m * l[m] * c[j - m];
      c[j] = cj / j;
      jf *= j;
      zj = (psi * z1 + z2) / j; z2 = z1; z1 = zj;
    }
    beta = (zj + beta) * ia; gam = (beta + gam) * ia;
    const double cjf = c[j] * jf, tT = cjf * beta, tP = cjf * gam;
    sT += tT; sP += tP;
    if (j >= 4 && std::fabs(tT) <= TILT_CONV_TOL * std::fabs(sT) && std::fabs(tP) <= TILT_CONV_TOL * std::fabs(sP)) {
      if (++quiet >= 2) { conv = true; break; }
    } else quiet = 0;
    if (!std::isfinite(sT) || !std::isfinite(sP) || !std::isfinite(jf)) return false;
  }
  if (!conv || !(sT > 0.0) || !(sP > 0.0)) return false;
  const double lm = lg + 0.5 * psi * psi, lT = lm + std::log(sT), lPf = lm + std::log(sP);
  if (!branch_negligible(u, lT, lPf)) return false;
  e.lT = lT;
  e.lPf = lPf;
  return true;
}

// Approximate median of W (Wilson-Hilferty for Y); only the far tails depend on it.
inline double median_W(double k, double p) {
  const double c = 1.0 - 1.0 / (9.0 * k);
  return c > 0.0 ? (k / p) * c * c * c : 0.0;
}

// Point values and gamma (+) normal pieces at distance d = u s sqrt(t).
inline bool endpoint(double t, double d, double s, double k, double tau, double p,
                     double u, bool lower, Endpoint& e, const Out* pt = nullptr) {
  // near k = 1 the image term is also wanted from the positive series, which
  // avoids the 1 / (1 - k) form when it converges
  const bool want_pimg = std::fabs(1.0 - k) < PIMG_BAND;
  e.u = u; e.lower = lower; e.both = false;
  if (u > 0.0) {
    const Out o = pt ? *pt : gamma_wald(t, d, s, k, tau, WANT_ALL);
    e.lF = o.logF; e.lS = o.logS;
    e.lD = o.logf + std::log(t / u);
  } else {
    e.lF = 0.0; e.lS = NEG_INF;
    e.lD = k * std::log(p) - LOG_SQRT_2PI + logG(k, p);
  }
  if (!std::isfinite(e.lD)) return false;
  const double z = p - u;
  // D at shape k + 1 seeds the upward recurrences (z < 0); the expansions do not use it
  auto need_D1 = [&]() {
    if (!(z < 0.0) || std::isfinite(e.lD1)) return true;
    const Out o1 = gamma_wald(t, d, s, k + 1.0, tau, WANT_PDF);
    e.lD1 = o1.logf + std::log(t / u);
    return std::isfinite(e.lD1);
  };
  const double sdY = std::sqrt(k) / p, y0 = k / p;
  if (p >= NARROW_P && sdY <= NARROW_SD && std::fabs(y0 - u) <= p / 3.0) {
    Endpoint en = e;
    if (narrow_conv(k, p, u, en)) {
      // keep the point kernel's density; the narrow series supplies C, T, Sf, Pf
      en.lD = e.lD;
      e = en;
      return true;
    }
  }
  // the preferred side, else the other one when the wanted tail is not small
  // (it is then formed as a complement).  The upper series is long for broad
  // laws, so it is used only when T may be small: T >= Phi(-u) and
  // T >= P(eps > -2) P(Y > u + 2).
  if (!lower && (u <= T_SMALL_U || std::log(0.977) + R::pgamma(u + 2.0, k, 1.0 / p, 0, 1) > LOG_T_SMALL))
    lower = true;
  // the expansions about u give T and Pf to full relative precision when they hold
  if (!lower && (wide_conv(k, p, u, e) || tilt_conv(k, p, u, e))) {
    e.lower = false;
    return std::isfinite(e.lT);
  }
  const double lp = std::log(p);
  // That bound is loose for broad laws (small k).  With z < 0 the lower ladder
  // runs about p u terms and the upper series about (41 / p)^2 (its terms fall
  // like e^(-p sqrt(n))); when the ladder is the shorter, its C gives T = 1 - C
  // directly, which keeps enough digits down to T = 1e-5.
  const double pu = p * u;
  if (!lower && z < 0.0 && pu + 8.0 * std::sqrt(pu) + 10.0 < (41.0 / p) * (41.0 / p)) {
    double lA0, lA1;
    if (!need_D1()) return false;
    if (ladder_lower(k, p, u, e.lD, e.lD1, lA0, lA1) && log1m_exp(lA0 - lp) > LOG_T_LADDER)
      lower = true;
  }
  for (int pass = 0; pass < 2; ++pass) {
    const bool lo = (pass == 0) ? lower : !lower;
    if (lo) {
      double lA0, lA1;
      if (!need_D1()) return false;
      if (!ladder_lower(k, p, u, e.lD, e.lD1, lA0, lA1)) continue;
      e.lC = lA0 - lp; e.lSf = lA1 - 2.0 * lp;
      if (pass == 1 && !(e.lC < LOG_COMPLEMENT_OK)) return false;
      double lB0, lB1, lB2;
      const int nser = want_pimg ? upper_length(k, z, u, PIMG_MAX_TERMS) : -1;
      if (nser >= 0 && nser <= PIMG_MAX_TERMS &&
          series_upper(k, p, u, e.lD, e.lD1, lB0, lB1, lB2))
        e.lPimg = e.lD + lB2;
    } else {
      double lB0, lB1, lB2;
      if (wide_conv(k, p, u, e) || tilt_conv(k, p, u, e)) {
        // the image term from the series only when that is short
      } else {
        if (!need_D1()) return false;
        if (!series_upper(k, p, u, e.lD, e.lD1, lB0, lB1, lB2)) continue;
        e.lT = e.lD + lB0; e.lPf = e.lD + lB1; e.lPimg = e.lD + lB2;
      }
      if (pass == 1 && !(e.lT < LOG_COMPLEMENT_OK)) return false;
    }
    e.lower = lo;
    return std::isfinite(lo ? e.lC : e.lT);
  }
  return false;
}

// Linear values at an endpoint relative to e^L: the exact side as given, the
// other side by complement (only used where those values are of order one).
struct Lin { double D, C, T, Sf, Pf, F, S, Pimg; };
inline Lin lin_at(const Endpoint& e, double L, double y0) {
  Lin v;
  v.D = std::exp(e.lD - L);
  v.F = std::exp(e.lF - L); v.S = std::exp(e.lS - L);
  v.Pimg = std::isfinite(e.lPimg) ? std::exp(e.lPimg - L) : NAN_D;
  if (e.both) {
    v.C = std::exp(e.lC - L); v.T = std::exp(e.lT - L);
    v.Sf = std::exp(e.lSf - L); v.Pf = std::exp(e.lPf - L);
  } else if (e.lower) {
    v.C = std::exp(e.lC - L); v.Sf = std::exp(e.lSf - L);
    v.T = std::exp(-L) - v.C; v.Pf = v.Sf + (y0 - e.u) * std::exp(-L);
  } else {
    v.T = std::exp(e.lT - L); v.Pf = std::exp(e.lPf - L);
    v.C = std::exp(-L) - v.T; v.Sf = v.Pf - (y0 - e.u) * std::exp(-L);
  }
  return v;
}

inline double lmax(std::initializer_list<double> xs) {
  double m = NEG_INF;
  for (double x : xs) if (x > m) m = x;
  return m;
}

// Closed form at shape k.  cond is the ratio of the summed term magnitudes to
// the smallest of the results formed (f, and F or S by difference), so about
// log10(cond) digits are lost to cancellation.
// ptB, when given, is the point kernel at (t, B) with WANT_ALL.
inline bool closed_A(double t, double B, double A, double s, double k, double tau,
                     AOut& r, double& cond, bool* used_ik = nullptr, const Out* ptB = nullptr) {
  const double rt = std::sqrt(t), ss = s * rt;
  const double p = tau * s / rt, uL = B / ss, uH = (B + A) / ss, Delta = A / ss;
  const double y0 = k / p, med = median_W(k, p), ik = 1.0 / (1.0 - k), aik = std::fabs(ik);
  const bool upper_all = uL > med;
  const bool lower_all = uH <= med;
  Endpoint eL, eH;
  if (!endpoint(t, B, s, k, tau, p, uL, !upper_all, eL, ptB)) return false;
  if (!endpoint(t, B + A, s, k, tau, p, uH, lower_all, eH)) return false;
  const double lDt = std::log(Delta * t), lDelta = std::log(Delta);
  // E[W 1{W <= u}], E[W 1{W > u}], Pi and Sg with their summed term magnitudes
  auto Lf = [&](const Lin& v, double u, double& m) { m += u * std::fabs(v.C) + std::fabs(v.Sf); return u * v.C - v.Sf; };
  auto Uf = [&](const Lin& v, double u, double& m) { m += u * std::fabs(v.T) + std::fabs(v.Pf); return u * v.T + v.Pf; };
  // Pi and Sg in the 1 / (1 - k) form, or with the series image term (Sg then
  // drops its constant p / (2 (1 - k)), so both endpoints must use it)
  auto Pi = [&](const Lin& v, double u, bool ser, double& m) {
    if (ser) { m += v.Pf + v.Pimg; return v.Pf + v.Pimg; }
    const double a = v.D, b = (0.5 * p - u) * v.F, c = (p - u) * v.T;
    m += std::fabs(v.Pf) + aik * (std::fabs(a) + std::fabs(b) + std::fabs(c));
    return v.Pf + ik * (a + b - c);
  };
  auto Sg = [&](const Lin& v, double u, bool ser, double& m) {
    if (ser) { m += std::fabs(v.Sf) + v.Pimg; return v.Sf + v.Pimg; }
    const double a = (p - u) * v.C, b = (u - 0.5 * p) * v.S, c = v.D;
    m += std::fabs(v.Sf) + aik * (std::fabs(a) + std::fabs(b) + std::fabs(c));
    return v.Sf + ik * (a + b + c);
  };
  // the better conditioned of the two forms of Pi(uL) - Pi(uH) or Sg(uH) - Sg(uL)
  bool ik_form = false;
  auto diff = [&](const Lin& a, const Lin& b, bool isPi, double& m) {
    const bool ser = !ISNAN(a.Pimg) && !ISNAN(b.Pimg);
    double best = NAN_D, mbest = R_PosInf;
    for (int form = 0; form < (ser ? 2 : 1); ++form) {
      double mm = 0.0;
      const double v = isPi ? Pi(a, uL, form == 1, mm) - Pi(b, uH, form == 1, mm)
                            : Sg(b, uH, form == 1, mm) - Sg(a, uL, form == 1, mm);
      if (v > 0.0 && mm / v < mbest) { best = v; mbest = mm / v; ik_form = (form == 0); }
    }
    m = ISNAN(best) ? R_PosInf : mbest * best;
    return best;
  };
  double mf = 0.0, mt = 0.0, nf, tail;
  bool tail_is_F;
  double L = 0.0;
  if (upper_all) {
    L = lmax({eL.lT, eL.lPf, eL.lD, eL.lF});
    const Lin a = lin_at(eL, L, y0), b = lin_at(eH, L, y0);
    nf = Uf(a, uL, mf) - Uf(b, uH, mf);
    tail = diff(a, b, true, mt);
    tail_is_F = true;
  } else if (lower_all) {
    L = lmax({eH.lC, eH.lSf, eH.lD, eH.lS});
    const Lin a = lin_at(eL, L, y0), b = lin_at(eH, L, y0);
    nf = Lf(b, uH, mf) - Lf(a, uL, mf);
    tail = diff(a, b, false, mt);
    tail_is_F = false;
  } else {
    const Lin a = lin_at(eL, 0.0, y0), b = lin_at(eH, 0.0, y0);
    mf = std::fabs(y0);
    nf = y0 - Lf(a, uL, mf) - Uf(b, uH, mf);
    double m1 = 0.0, m2 = 0.0;
    const double fa = diff(a, b, true, m1), sa = diff(a, b, false, m2);
    tail_is_F = !(fa > sa);
    tail = tail_is_F ? fa : sa;
    mt = tail_is_F ? m1 : m2;
  }
  if (!(nf > 0.0) || !(tail > 0.0)) return false;
  cond = std::max(mf / nf, mt / tail);
  if (used_ik) *used_ik = ik_form;
  r.logf = std::log(nf) + L - lDt;
  const double lt = std::min(0.0, std::log(tail) + L - lDelta);
  if (tail_is_F) { r.logF = lt; r.logS = log1m_exp(lt); }
  else           { r.logS = lt; r.logF = log1m_exp(lt); }
  return std::isfinite(r.logf) && !ISNAN(r.logF) && !ISNAN(r.logS);
}

inline double lse2(double a, double b) {
  const double m = std::max(a, b);
  if (m == NEG_INF) return NEG_INF;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// Gauss-Legendre average of the point kernel over [B, B + A] (n = 4 or 8 nodes).
inline AOut point_average(double t, double B, double A, double s, double k, double tau, int n) {
  static const double x4[] = {-0.86113631159405257522, -0.33998104358485626480,
    0.33998104358485626480, 0.86113631159405257522};
  static const double w4[] = {0.17392742256872692869, 0.32607257743127307131,
    0.32607257743127307131, 0.17392742256872692869};
  static const double x8[] = {-0.96028985649753623168, -0.79666647741362673959, -0.52553240991632898582,
    -0.18343464249564980494, 0.18343464249564980494, 0.52553240991632898582,
    0.79666647741362673959, 0.96028985649753623168};
  static const double w8[] = {0.05061426814518812957, 0.11119051722668723527, 0.15685332293894364367,
    0.18134189168918099149, 0.18134189168918099149, 0.15685332293894364367,
    0.11119051722668723527, 0.05061426814518812957};
  const double* xs = (n == 4) ? x4 : x8;
  const double* ws = (n == 4) ? w4 : w8;
  AOut r;
  r.logf = NEG_INF; r.logF = NEG_INF; r.logS = NEG_INF;
  for (int i = 0; i < n; ++i) {
    const double d = B + 0.5 * A * (1.0 + xs[i]), lw = std::log(ws[i]);
    const Out o = gamma_wald(t, d, s, k, tau, WANT_ALL);
    r.logf = lse2(r.logf, lw + o.logf);
    r.logF = lse2(r.logF, lw + o.logF);
    r.logS = lse2(r.logS, lw + o.logS);
  }
  r.logF = std::min(0.0, r.logF); r.logS = std::min(0.0, r.logS);
  return r;
}

// Closed form at k, or across the removable singularity at k = 1 by cubic
// interpolation in k from k = 1 -+ K_BAND, 1 -+ 2 K_BAND (the density, which has
// no 1 / (1 - k) term, at k itself).
inline bool closed_or_band(double t, double B, double A, double s, double k, double tau,
                           AOut& r, double& cond, const Out* ptB = nullptr) {
  bool used_ik = true;
  const bool ok = closed_A(t, B, A, s, k, tau, r, cond, &used_ik, ptB);
  if (std::fabs(1.0 - k) >= K_BAND || (ok && cond < COND_MAX && !used_ik)) return ok;
  const double kn[4] = {1.0 - 2.0 * K_BAND, 1.0 - K_BAND, 1.0 + K_BAND, 1.0 + 2.0 * K_BAND};
  AOut rn[4];
  double cmax = 0.0, c;
  if (!ok) return false;
  const AOut r0 = r;
  for (int i = 0; i < 4; ++i) {
    if (!closed_A(t, B, A, s, kn[i], tau, rn[i], c)) return false;
    cmax = std::max(cmax, c);
  }
  // interpolate the log of the smaller tail
  const bool useS = rn[1].logS < rn[1].logF;
  double v = 0.0;
  for (int i = 0; i < 4; ++i) {
    double w = 1.0;
    for (int j = 0; j < 4; ++j) if (j != i) w *= (k - kn[j]) / (kn[i] - kn[j]);
    v += w * (useS ? rn[i].logS : rn[i].logF);
  }
  v = std::min(0.0, v);
  r.logf = r0.logf;
  if (useS) { r.logS = v; r.logF = log1m_exp(v); }
  else      { r.logF = v; r.logS = log1m_exp(v); }
  cond = cmax;
  r.path = A_BAND;
  return true;
}

// Driver.  B >= 0 is the near distance, A >= 0 the width of the start range.
// want and lf_floor as for gamma_wald: a log density below lf_floor (with want
// = WANT_PDF) may come back as an upper bound still below it.
inline AOut gamma_wald_A(double t, double B, double A, double s, double k, double tau,
                         int want = WANT_ALL, double lf_floor = NEG_INF) {
  AOut r;
  if (!(t > 0.0)) { r.logf = NEG_INF; r.logF = NEG_INF; r.logS = 0.0; return r; }
  if (!(A > 0.0)) {
    const Out o = gamma_wald(t, B, s, k, tau, want, 0, lf_floor);
    r.logf = o.logf; r.logF = o.logF; r.logS = o.logS; r.path = A_POINT;
    return r;
  }
  const double rt = std::sqrt(t), ss = s * rt;
  const double p = tau * s / rt, uH = (B + A) / ss, Delta = A / ss;
  // Early times.  The first-passage time grows with the distance on every
  // path, so F(t | d) <= F(t | B) on [B, B + A]; and for each drift V the
  // start-averaged density (1/A) int (d / t) phi_ss(d - V t) dd is at most
  // (B + A) / (t A) P(V t + ss Z >= B) <= (B + A) / (t A) F(t | B), the Wald
  // CDF being at least Phi((V t - B) / ss).  Both bounds are returned when
  // they are negligible on the scale the caller asked for (as far_tail does).
  // the point values at B are kept for the closed form's near endpoint
  Out oB;
  bool have_B = false;
  if (B > 0.0 && k / p < B / ss) {
    oB = gamma_wald(t, B, s, k, tau, WANT_ALL);
    have_B = true;
    const Out& o = oB;
    const double lF = o.logF, lfb = lF + std::log((B + A) / (t * A));
    double thr = LOG_NEGLIGIBLE;
    if (want == WANT_SURV) thr = LOG_UNIT;
    if (want == WANT_PDF) thr = std::max(thr, lf_floor);
    const bool far = (lF < LOG_NEGLIGIBLE && lfb < LOG_NEGLIGIBLE) ||
      (want == WANT_PDF ? lfb < thr : (want == WANT_SURV && lF < thr));
    if (far) {
      r.logf = lfb; r.logF = lF; r.logS = 0.0; r.path = A_FAR;
      return r;
    }
  }
  // the span in units of the scale on which f, F and S vary with the distance
  const double span = Delta * (1.0 + uH + std::fabs(k / p - uH) + p);
  if (span < SHORT_SPAN4) {
    r = point_average(t, B, A, s, k, tau, 4);
    r.path = A_SHORT;
    return r;
  }
  double cond = R_PosInf;
  if (span >= SHORT_SPAN8 && closed_or_band(t, B, A, s, k, tau, r, cond, have_B ? &oB : nullptr) && cond < COND_MAX) {
    r.path = (std::fabs(1.0 - k) < K_BAND && r.path == A_BAND) ? A_BAND : A_CLOSED;
    return r;
  }
  if (span < SHORT_SPAN8) {
    r = point_average(t, B, A, s, k, tau, 8);
    r.path = A_SHORT;
    return r;
  }
  // The span overstates the panels needed when the drift law is broad, so the
  // count doubles from 2 until two counts agree, up to the span's count.
  const int pmax = std::min(FALLBACK_PANELS, std::max(2, static_cast<int>(std::ceil(span))));
  auto composite = [&](int panels) {
    AOut c;
    c.logf = NEG_INF; c.logF = NEG_INF; c.logS = NEG_INF;
    for (int j = 0; j < panels; ++j) {
      const AOut q = point_average(t, B + A * j / panels, A / panels, s, k, tau, 8);
      c.logf = lse2(c.logf, q.logf); c.logF = lse2(c.logF, q.logF); c.logS = lse2(c.logS, q.logS);
    }
    const double lpn = std::log(static_cast<double>(panels));
    c.logf -= lpn; c.logF = std::min(0.0, c.logF - lpn); c.logS = std::min(0.0, c.logS - lpn);
    return c;
  };
  auto same = [](double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) return a == b;
    return std::fabs(a - b) <= FALLBACK_TOL * std::max(1.0, std::fabs(a));
  };
  int panels = std::min(2, pmax);
  r = composite(panels);
  while (panels < pmax) {
    panels = std::min(2 * panels, pmax);
    const AOut q = composite(panels);
    const bool agree = same(q.logf, r.logf) && same(q.logF, r.logF) && same(q.logS, r.logS);
    r = q;
    if (agree) break;
  }
  r.path = A_FALLBACK;
  return r;
}

}  // namespace rdmswtn_gamma
#endif
