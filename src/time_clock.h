#ifndef EMC2_TIME_CLOCK_H
#define EMC2_TIME_CLOCK_H

// Deterministic operational clocks q(x) for exact time-changed race models.
// A model evaluated on the clock has F(t) = F_base(q(x)) and
// f(t) = f_base(q(x)) q'(x), x = t - t0.  One header serves the likelihood
// wrappers (time_warp.cpp), the scalar R exports and the simulators, which
// invert the clock.  It includes nothing from Rcpp, so light translation units
// can share it.
//
//   CLOCK_POWER        ballistic warp (Math/ballistic-time.md), parameter eta:
//                      q = ((1+x)^omega - 1)/omega, omega = exp(eta)
//   CLOCK_LINEAR       urgency, RDMSWTN_UT(clock = "linear"), parameter u:
//                      q = x + u x^2/2,          q' = 1 + u x
//   CLOCK_EXPONENTIAL  urgency, RDMSWTN_UT(clock = "exponential"), parameter u:
//                      q = expm1(u x)/u,         q' = exp(u x)
//   CLOCK_EXHAUSTION   RDMSWTN_TT, parameter tau:
//                      q = x - x^2/(2 tau) on x < tau, frozen at tau/2 after;
//                      q' = 1 - x/tau, zero from x = tau on
//
// Every clock is the identity at one parameter value (eta = 0, u = 0,
// tau = +Inf), and that value short-circuits to the exact parent model.

#include <R.h>
#include <Rmath.h>
#include <cmath>

namespace emc2tw {

enum Clock : int {
  CLOCK_POWER = 0,
  CLOCK_LINEAR = 1,
  CLOCK_EXPONENTIAL = 2,
  CLOCK_EXHAUSTION = 3
};

// --- CLOCK_POWER ----------------------------------------------------------

// s = c_eta(u) = ((1+u)^omega - 1)/omega, omega = exp(eta).
// The log1p/expm1 forms avoid direct exponentiation's overflow and cancellation regions.
inline double fwd(double u, double eta) {
  if (eta == 0.0) return u;                  // exact parent, bitwise
  if (!(u > 0.0)) return u;                  // u <= 0 and NaN pass through
  if (!R_FINITE(u)) return u;                // c(+Inf) = +Inf for every omega > 0
  if (std::isnan(eta)) return R_NaN;         // malformed eta propagates
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return R_PosInf;     // omega -> Inf
  if (omega < 1e-12) return std::log1p(u);   // omega -> 0
  const double x = omega * std::log1p(u);
  if (x > 709.0) return R_PosInf;            // (1+u)^omega overflows
  return std::expm1(x) / omega;
}

// log c'_eta(u) = (omega - 1) log(1 + u).
inline double log_jac(double u, double eta) {
  if (eta == 0.0) return 0.0;
  if (!(u > 0.0) || !R_FINITE(u)) return 0.0;
  if (std::isnan(eta)) return R_NaN;
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return R_PosInf;
  return (omega - 1.0) * std::log1p(u);
}

// u = c_eta^{-1}(s) = (1 + omega s)^(1/omega) - 1.
inline double inv(double s, double eta) {
  if (eta == 0.0) return s;
  if (!(s > 0.0)) return s;
  if (!R_FINITE(s)) return s;
  if (std::isnan(eta)) return R_NaN;
  const double omega = std::exp(eta);
  if (!R_FINITE(omega)) return 0.0;          // omega -> Inf
  if (omega < 1e-12) return (s > 709.0) ? R_PosInf : std::expm1(s);
  const double os = omega * s;
  // When omega*s overflows, log1p(omega*s) == log(omega) + log(s)
  // to well past double precision.
  const double y = (R_FINITE(os) ? std::log1p(os)
                                 : (std::log(omega) + std::log(s))) / omega;
  if (y > 709.0) return R_PosInf;
  return std::expm1(y);
}

// --- The clock family -----------------------------------------------------
// x is decision time (x > 0, +Inf allowed) and p the clock parameter.  The
// forward map, log-Jacobian and inverse assume clock_is_valid(clock, p).

inline bool clock_is_identity(int clock, double p) {
  return clock == CLOCK_EXHAUSTION ? p == R_PosInf : p == 0.0;
}

// The power clock keeps its historical contract: every eta, NaN included,
// reaches the kernels (NaN propagates).  The RDMSWTN clocks reject what their
// R bounds already exclude, so a malformed row is a zero likelihood.
inline bool clock_is_valid(int clock, double p) {
  switch (clock) {
    case CLOCK_LINEAR:
    case CLOCK_EXPONENTIAL: return R_FINITE(p) && p >= 0.0;
    case CLOCK_EXHAUSTION:  return p > 0.0;  // +Inf allowed, NaN rejected
    default:                return true;
  }
}

inline double clock_fwd(int clock, double x, double p) {
  switch (clock) {
    case CLOCK_LINEAR:
      if (p == 0.0) return x;
      return x * (1.0 + 0.5 * p * x);
    case CLOCK_EXPONENTIAL:
      if (p == 0.0) return x;
      return std::expm1(p * x) / p;
    case CLOCK_EXHAUSTION:
      if (!R_FINITE(p)) return x;
      return (x >= p) ? 0.5 * p : x * (1.0 - 0.5 * x / p);
    default:
      return fwd(x, p);
  }
}

inline double clock_log_jac(int clock, double x, double p) {
  switch (clock) {
    case CLOCK_LINEAR:      return std::log1p(p * x);
    case CLOCK_EXPONENTIAL: return p * x;
    case CLOCK_EXHAUSTION:
      if (!R_FINITE(p)) return 0.0;
      return (x >= p) ? R_NegInf : std::log1p(-x / p);
    default:
      return log_jac(x, p);
  }
}

// Supremum of q: the operational-time budget.  Finite only for exhaustion,
// where it is the ordinary-process time at which the clock freezes.
inline double clock_budget(int clock, double p) {
  return (clock == CLOCK_EXHAUSTION && R_FINITE(p)) ? 0.5 * p : R_PosInf;
}

// x = q^{-1}(y); +Inf for y beyond the budget (the accumulator never
// finishes).  The quadratic roots use the cancellation-free forms.
inline double clock_inv(int clock, double y, double p) {
  if (!(y > 0.0) || clock_is_identity(clock, p)) return y;
  switch (clock) {
    case CLOCK_LINEAR: {
      if (!R_FINITE(y)) return y;
      const double py2 = 2.0 * p * y;
      if (!R_FINITE(py2)) return std::sqrt(2.0) * std::sqrt(y) / std::sqrt(p);
      return 2.0 * y / (1.0 + std::sqrt(1.0 + py2));
    }
    case CLOCK_EXPONENTIAL: {
      if (!R_FINITE(y)) return y;
      const double py = p * y;
      return (R_FINITE(py) ? std::log1p(py) : std::log(p) + std::log(y)) / p;
    }
    case CLOCK_EXHAUSTION: {
      if (y > 0.5 * p) return R_PosInf;
      return 2.0 * y / (1.0 + std::sqrt(std::fmax(0.0, 1.0 - 2.0 * y / p)));
    }
    default:
      return inv(y, p);
  }
}

}  // namespace emc2tw

#endif  // EMC2_TIME_CLOCK_H
