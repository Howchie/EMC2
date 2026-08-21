#ifndef model_rng_h
#define model_rng_h

// RNG simulation kernels for posterior prediction (rfun port).
//
// Pure scalar primitives only -- no Rcpp container types in the inner loops.
// All draws go through R's RNG (R::rnorm, R::runif, R::rgamma, R::qnorm, ...)
// so set.seed() reproducibility and single-threaded semantics are preserved.
// GetRNGstate/PutRNGstate is handled by Rcpp attributes (Rcpp::RNGScope) at
// the [[Rcpp::export]] boundary in model_rng.cpp -- callers of these
// primitives must be reached through an exported function.
//
// These kernels are distributionally equivalent to, but NOT stream-identical
// with, the R reference rfuns (.lba_rfun, rBAwL, rRDM/rWald, rRDMSWTN/rSWTN):
// the order and shape of RNG draws differs (msm::rtnorm uses internal
// rejection; rwaldt draws vectors per-subset rather than per-trial). See
// rfun_port_plan.md section 3.

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include "wald_functions.h"  // pnorm_std
#include "bawd_kernel.h"     // finite-rho scalar arithmetic, shared with the likelihood

// One-sided truncated normal N(mu, sd) on [lo, Inf). lo == R_NegInf means
// "no truncation" (plain rnorm). Uses inverse-CDF for modest truncation and
// Robert (1995) exponential rejection when (lo - mu)/sd is far in the tail,
// to avoid qnorm tail-cancellation. Distributionally equivalent to
// msm::rtnorm(1, mu, sd, lower = lo).
inline double rtnorm_lower_r(double mu, double sd, double lo) {
  if (!(sd > 0.0)) return mu;
  if (lo == R_NegInf) return R::norm_rand() * sd + mu;

  const double z_lo = (lo - mu) / sd;

  if (z_lo > 4.0) {
    // Robert (1995): shifted-exponential proposal with optimal rate alpha.
    const double alpha = (z_lo + std::sqrt(z_lo * z_lo + 4.0)) / 2.0;
    for (int iter = 0; iter < 1000; iter++) {
      const double e = R::exp_rand();
      const double z = z_lo + e / alpha;
      const double rho = std::exp(-0.5 * (z - alpha) * (z - alpha));
      if (R::unif_rand() <= rho) return mu + sd * z;
    }
    return mu + sd * z_lo;  // pathological fallback, essentially unreachable
  }

  const double p_lo = pnorm_std(z_lo, true, false);
  double p = p_lo + R::unif_rand() * (1.0 - p_lo);
  p = std::fmin(std::fmax(p, 1e-300), 1.0 - 1e-16);
  const double z = R::qnorm(p, 0.0, 1.0, 1, 0);
  return mu + sd * z;
}

// Wald/inverse-Gaussian first passage time: criterion k, rate l, diffusion s.
// Michael-Schucany-Haas sampler, matches R's rwaldt (model_RDM.R:33-63)
// distributionally (draw order differs).
inline double rwald_fpt_r(double k, double l, double s) {
  const double tiny = 1e-6;
  if (l <= tiny) {
    // Degenerate zero-drift case -> Levy distribution.
    const double c = (k / s) * (k / s);
    const double u = R::unif_rand();
    const double q = R::qnorm(1.0 - u / 2.0, 0.0, 1.0, 1, 0);
    return c / (q * q);
  }
  const double mu = k / l;
  const double lambda = (k / s) * (k / s);
  double y = R::norm_rand();
  y = y * y;
  double x0 = mu + mu * mu * y / (2.0 * lambda) -
              std::sqrt(4.0 * mu * lambda * y + mu * mu * y * y) * mu / (2.0 * lambda);
  if (!(x0 > 0.0)) x0 = mu;  // numerical edge case, R's x[x<0]<-max(x) has no scalar analogue
  const double z = R::unif_rand();
  const double test = mu / (mu + x0);
  return (z > test) ? (mu * mu / x0) : x0;
}

// Full rWald semantics for a single accumulator: start-point range,
// posdrift branching, and defective (never-hit) negative-drift handling.
// Matches R's rWald (model_RDM.R:30-99).
inline double rwald_acc_r(double B, double v, double A, double s, bool posdrift) {
  const double bs = B + R::unif_rand() * A;
  const bool pos = posdrift ? (v > 0.0) : (v >= 0.0);
  if (pos) return rwald_fpt_r(bs, v, s);
  if (!posdrift && v < 0.0) {
    const double p_hit = std::exp(2.0 * v * bs / (s * s));
    if (R::unif_rand() < p_hit) return rwald_fpt_r(bs, -v, s);
    return R_PosInf;
  }
  return R_PosInf;  // posdrift && v <= 0: defective, never hits
}

// Erlang clock draw shared by BAwL/RDMSWTN guess & kill timers
// (matches .rdmswtn_erlang_omega mixture semantics, model_RDM.R:791-800):
//   shape 1/2: rgamma(shape, rate = lambda)
//   shape 3 (mixed): with prob omega use shape 1 rate lambda,
//                    else shape 2 rate 2*lambda.
// lambda <= 0 or NA -> +Inf (clock disabled).
// NOTE: R::rgamma(shape, scale) takes SCALE not rate -- inverted explicitly.
inline double rerlang_clock_r(double lambda, int shape, double omega) {
  if (ISNAN(lambda) || !(lambda > 0.0)) return R_PosInf;
  int use_shape = shape;
  double rate = lambda;
  if (shape == 3) {
    if (R::unif_rand() <= omega) {
      use_shape = 1;
      rate = lambda;
    } else {
      use_shape = 2;
      rate = 2.0 * lambda;
    }
  }
  return R::rgamma((double)use_shape, 1.0 / rate);
}

// BAwD first passage: solve V q_rho(u) - ell c_{rho,gamma}(u) = d, where
// q_rho and c_{rho,gamma} integrate the selected base kernel.  We solve for
// the FIRST (rising-limb) crossing of d = b - z.
//
// The rho = Inf branch below retains the historical exponential arithmetic.
// Finite rho uses the shared scalar kernel and its closed-form gamma = 1
// inversion, with bracketed Newton on the rising limb for gamma < 1.
inline double bawd_qf_r(double u, double k) {
  if (k <= 1e-10) return u;
  return -std::expm1(-k * u) / k;
}

inline double bawd_cf_r(double u, double k, double gamma) {
  if (k <= 1e-10 || gamma <= 1e-12) return u;
  if (gamma >= 1.0 - 1e-12) return bawd_qf_r(u, k);
  return -std::expm1(-gamma * k * u) / (gamma * k);
}

// BAwF first-passage time: the smallest u > 0 with h_rho(u)[z + V u] = b, i.e.
// the smallest root of F(x) = b H_rho(x) - z - V x / k in the phase x = k u.
//
// Unlike every other ballistic model here, this one is NOT a function of the
// distance b - z alone: the fading multiplies the start point too, so b and z
// enter separately.
//
// F is convex (H'' > 0) with F(0) = b - z > 0, and the first crossing can only
// occur at or before the tangency phase x_sat(z) <= x_max.  So a crossing
// exists if and only if F(x_max) <= 0, and the bracket [0, x_max] both decides
// existence and contains the root -- which is also why no time beyond
// T_max = x_max / k is ever returned.  Keeping the bracket at x_max (<= 2)
// rather than at the unconstrained minimiser avoids overflowing H for large V.
inline double bawf_hit_time_r(double V, double b, double z, double k,
                              double rho) {
  if (!(b > z)) return 0.0;                   // start already at threshold
  if (ISNAN(V) || !(V > 0.0)) return R_PosInf;
  if (k <= 1e-10) return (b - z) / V;         // exact LBA limit
  const bool rho_inf = (rho == 0.0) || (rho > 0.0 && !R_FINITE(rho));
  if (!rho_inf && !(rho > 1.0)) return R_PosInf;
  const double x_max = rho_inf ? 1.0 : rho / (rho - 1.0);
  const double VK = V / k;
  auto Fx = [&](double x) {
    const double H = rho_inf ? std::exp(x)
                             : std::exp(rho * std::log1p(x / rho));
    return b * H - z - VK * x;
  };
  auto Fp = [&](double x) {
    const double Hp = rho_inf ? std::exp(x)
                              : std::exp((rho - 1.0) * std::log1p(x / rho));
    return b * Hp - VK;
  };
  if (Fx(x_max) > 0.0) return R_PosInf;       // peak never reaches threshold
  double lo = 0.0, hi = x_max, x = 0.5 * x_max;
  for (int it = 0; it < 100; ++it) {
    const double f = Fx(x);
    if (f > 0.0) lo = x; else hi = x;
    const double d1 = Fp(x);
    double xn = (d1 < 0.0 && R_FINITE(d1)) ? x - f / d1 : 0.5 * (lo + hi);
    if (!(xn > lo) || !(xn < hi) || !R_FINITE(xn)) xn = 0.5 * (lo + hi);
    const bool done = std::fabs(xn - x) <= 1e-14 * std::fmax(1.0, xn);
    x = xn;
    if (done) break;
  }
  return x / k;
}

// BAwR first-passage time: the smallest u > 0 with
// G(u) = V u - kappa u^(p+1)/(p+1) equal to the distance d = b - z.
//
// Unlike BAwF (where the fading multiplies the start point) the decay here is
// a pure function of elapsed time, so the crossing depends on b and z only
// through their difference and this takes `d` like the rest of the family.
//
// G rises from 0, peaks at the trajectory peak u_peak = (V/kappa)^(1/p), and
// falls thereafter, so the first crossing is in [0, u_peak] and exists if and
// only if G(u_peak) >= d.  That maximum is closed form:
//
//   G(u_peak) = u_peak V p / (p + 1),
//
// because kappa u_peak^p = V by construction.  G is strictly increasing and
// concave on the bracket, so a bracketed Newton is unconditionally safe.
inline double bawr_hit_time_r(double V, double d, double kappa, double pw) {
  if (!(d > 0.0)) return 0.0;                 // start already at threshold
  if (ISNAN(V) || !(V > 0.0)) return R_PosInf;
  if (!(pw > 0.0)) return R_PosInf;           // p = 0 is not in the family
  if (kappa <= 1e-10) return d / V;           // exact LBA limit
  const double u_peak = std::exp((std::log(V) - std::log(kappa)) / pw);
  // An unrepresentable peak means the decay is negligible over any reachable
  // time; the trajectory is a straight line there, which is the LBA answer.
  if (!R_FINITE(u_peak) || !(u_peak > 0.0)) return d / V;
  if (u_peak * V * pw / (pw + 1.0) < d) return R_PosInf;  // never reaches b
  auto Gx = [&](double u) {
    return V * u - kappa * std::exp((pw + 1.0) * std::log(u)) / (pw + 1.0);
  };
  auto Gp = [&](double u) {
    return V - kappa * std::exp(pw * std::log(u));
  };
  double lo = 0.0, hi = u_peak, u = 0.5 * u_peak;
  for (int it = 0; it < 100; ++it) {
    const double f = Gx(u) - d;
    if (f < 0.0) lo = u; else hi = u;
    const double d1 = Gp(u);
    double un = (d1 > 0.0 && R_FINITE(d1)) ? u - f / d1 : 0.5 * (lo + hi);
    if (!(un > lo) || !(un < hi) || !R_FINITE(un)) un = 0.5 * (lo + hi);
    const bool done = std::fabs(un - u) <= 1e-14 * std::fmax(1.0, un);
    u = un;
    if (done) break;
  }
  return u;
}

inline double bawd_hit_time_r(double V, double d, double k, double ell,
                              double gamma, double rho) {
  if (!(d > 0.0)) return 0.0;                 // start already at threshold
  if (ISNAN(V) || !(V > 0.0)) return R_PosInf;

  // FFI sentinel: rho == 0.0 (the compileAttributes-safe default) or a
  // positive Inf both select the exponential kernel -- see bawd_geometry()
  // in model_BAwD.h for the identical test.  Everything below this branch
  // is the pre-existing exponential body, untouched (R1: no shared
  // arithmetic with the finite-rho kernel).
  const bool rho_inf = (rho == 0.0) || (rho > 0.0 && !R_FINITE(rho));
  if (rho_inf) {
    if (k <= 1e-10) {                           // exact LBA limit, drift V - ell
      return (V > ell) ? d / (V - ell) : R_PosInf;
    }
    if (ell <= 1e-12) {                         // exact BAwL/leak-free inversion
      const double x = 1.0 - k * d / V;
      return (x > 0.0) ? -std::log(x) / k : R_PosInf;
    }
    if (gamma >= 1.0 - 1e-12) {
      if (!(V > ell)) return R_PosInf;
      const double x = 1.0 - k * d / (V - ell);
      return (x > 0.0) ? -std::log(x) / k : R_PosInf;
    }
    if (!(V > ell)) return R_PosInf;
    const double u_p = (gamma <= 1e-12)
      ? std::log(V / ell) / k
      : std::log(V / ell) / ((1.0 - gamma) * k);
    if (V * bawd_qf_r(u_p, k) - ell * bawd_cf_r(u_p, k, gamma) < d)
      return R_PosInf;

    double u = d / (V - ell);
    if (!(u > 0.0) || !R_FINITE(u)) u = 1e-12;
    for (int it = 0; it < 100; ++it) {
      const double f = V * bawd_qf_r(u, k) -
        ell * bawd_cf_r(u, k, gamma) - d;
      const double fp = V * std::exp(-k * u) -
        ell * std::exp(-gamma * k * u);
      if (!(fp > 0.0)) break;                   // at the peak: root is u_p
      double un = u - f / fp;
      if (!(un > 0.0)) un = 0.5 * u;
      if (un > u_p) un = u_p;
      const bool done = std::fabs(un - u) <= 1e-13 * std::fmax(1.0, un);
      u = un;
      if (done) break;
    }
    return u;
  }

  // Finite rho >= 1 (bawd_geometry already rejects rho in (-Inf, 1) and
  // NaN before this is reached from the likelihood side; guard here too
  // since the simulator can be called directly).
  if (!(rho >= 1.0)) return R_PosInf;
  const bool rho_one = std::fabs(rho - 1.0) <= 1e-12;
  if (k <= 1e-10)
    return (V > ell) ? d / (V - ell) : R_PosInf;

  if (ell <= 1e-12) {
    // With no clearance the finite-rho clock is V*kQ(x)/k.  Its
    // rho=1 member is unbounded; rho>1 has a finite total clock.
    if (rho_one) {
      return std::expm1(k * d / V) / k;
    }
    const double S = k * d * (rho - 1.0) / (V * rho);
    if (!(S < 1.0)) return R_PosInf;
    const double x = rho * (std::pow(1.0 - S,
                                      -1.0 / (rho - 1.0)) - 1.0);
    return x / k;
  }

  if (gamma >= 1.0 - 1e-12) {
    // Closed-form inversion of (V - ell) * kQ_rho(x) = k * d; no root-find.
    if (!(V > ell)) return R_PosInf;
    if (rho_one) {                              // kQ(x) = log(1 + x), unbounded
      const double x = std::expm1(k * d / (V - ell));
      return x / k;
    }
    const double S = k * d * (rho - 1.0) / ((V - ell) * rho);
    if (!(S < 1.0)) return R_PosInf;            // omission: never reaches d
    const double x = rho * (std::pow(1.0 - S, -1.0 / (rho - 1.0)) - 1.0);
    return x / k;
  }

  // gamma < 1: same omission-then-Newton shape as the exponential branch
  // above, but the peak and the trajectory both route through the shared
  // finite-rho kernel (bawd_kernel.h) instead of the exponential formulas.
  if (!(V > ell)) return R_PosInf;
  const double x_p = bawd_pk_peak_x(rho, gamma, std::log(V / ell));
  const double L_p = bawd_pk_log_tau(rho, x_p);
  const double height =
    (V * bawd_pk_kq(rho, L_p) - ell * bawd_pk_kr(rho, gamma, L_p)) / k;
  if (height < d) return R_PosInf;

  const double u_p = x_p / k;
  double u = d / (V - ell);
  if (!(u > 0.0) || !R_FINITE(u)) u = 1e-12;
  for (int it = 0; it < 100; ++it) {
    const double L = bawd_pk_log_tau(rho, k * u);
    const double f = V * bawd_pk_kq(rho, L) / k -
      ell * bawd_pk_kr(rho, gamma, L) / k - d;
    const double fp = V * std::exp(-rho * L) -
      ell * std::exp(-rho * gamma * L);
    if (!(fp > 0.0)) break;                   // at the peak: root is u_p
    double un = u - f / fp;
    if (!(un > 0.0)) un = 0.5 * u;
    if (un > u_p) un = u_p;
    const bool done = std::fabs(un - u) <= 1e-13 * std::fmax(1.0, un);
    u = un;
    if (done) break;
  }
  return u;
}

// BAwDp first passage.  The evidence trajectory is monotone only up to the
// universal freeze time u* = log(1/lambda)/k.  We invert the strictly
// increasing elementary clock m(u) on [0, u*] with a safeguarded Newton step;
// the likelihood itself never needs this inversion (its PDF/CDF are the
// closed-form LBA change of variables).
inline double bawdp_hit_time_r(double V, double d, double k, double lambda) {
  if (!(d > 0.0) || ISNAN(V) || !(V > 0.0) || !(k >= 0.0) ||
      !(lambda >= 0.0) || !(lambda < 1.0)) return R_PosInf;
  if (k <= 1e-10) {
    return d / (V * (1.0 - lambda));
  }
  const double target = d / V;
  if (!(target > 0.0)) return R_PosInf;
  if (lambda == 0.0) {
    const double x = 1.0 - k * target;
    return (x > 0.0) ? -std::log(x) / k : R_PosInf;
  }
  const double log_lambda = std::log(lambda);
  const double u_star = -log_lambda / k;
  const double m_max = (-std::expm1(log_lambda) + lambda * log_lambda) / k;
  if (!(target > 0.0) || !(target < m_max)) return R_PosInf;

  auto m_at = [k, lambda](double u) {
    return -std::expm1(-k * u) / k - lambda * u;
  };
  double lo = 0.0, hi = u_star;
  double u = target / (1.0 - lambda);
  if (!(u > lo) || !(u < hi) || !R_FINITE(u)) u = 0.5 * (lo + hi);
  for (int it = 0; it < 100; ++it) {
    const double f = m_at(u) - target;
    if (f > 0.0) hi = u; else lo = u;
    const double fp = std::exp(-k * u) - lambda;
    double un = (fp > 0.0) ? u - f / fp : 0.5 * (lo + hi);
    if (!(un > lo) || !(un < hi) || !R_FINITE(un)) un = 0.5 * (lo + hi);
    if (std::fabs(un - u) <= 1e-13 * std::fmax(1.0, un)) return un;
    u = un;
  }
  return u;
}

#endif
