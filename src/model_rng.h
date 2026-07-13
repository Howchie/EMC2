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

#endif
