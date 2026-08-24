#ifndef rdm_h
#define rdm_h

// Shared RDM/SWTN/GBM model kernels and inline helpers.  Rcpp-exported
// wrappers are defined in model_RDM.cpp.
#define _USE_MATH_DEFINES
#include <cmath>
#include <Rcpp.h>
#include "gsl_utils.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "wald_functions.h"
#include "composite_functions.h"
#include "gaussian.h"
#include "gl_quad.h"
using namespace Rcpp;
// Exported definitions carry the Rcpp defaults in model_RDM.cpp.  Suppress
// defaults while that translation unit parses this public declaration block so
// its attributed declarations can own the defaults without a duplicate.
#ifdef RDM_NO_DEFAULT_ARGUMENTS
#define RDM_DEFAULT_ARGUMENT(value)
#else
#define RDM_DEFAULT_ARGUMENT(value) = value
#endif


static constexpr double RDM_Q_EPSILON = 1e-8;

template <typename DensityFn>
inline double integrate_density_gl20_finite(double t_upper, DensityFn&& density_fn);

// --------------------------------------------------------------------------
// Analytic helpers for SPV + Erlang-2 killed Wald CDF
//
// Computes:
//   J(q) = ∫_lo^hi exp(q*x) Phi(a*x + c) dx
//   K(q) = ∫_lo^hi x exp(q*x) Phi(a*x + c) dx
//
// Then:
//   ∫ exp[eta*(b-x)] Phi(a*x+c) dx
//     = exp(eta*b) J(-eta)
//
//   ∫ (b-x) exp[eta*(b-x)] Phi(a*x+c) dx
//     = exp(eta*b) { b J(-eta) - K(-eta) }
//
// Assumes a > 0 and usually b > x_hi.
// --------------------------------------------------------------------------

inline signed_log slog_int_exp_pnorm(
    double q, double a, double c,
    double x_lo, double x_hi
) {
  if (x_hi <= x_lo) return make_signed_log(R_NegInf, 0);

  const double y_lo = a * x_lo + c;
  const double y_hi = a * x_hi + c;

  // q = 0:
  // ∫ Phi(a*x+c) dx = [y Phi(y) + phi(y)] / a
  if (std::abs(q) <= RDM_Q_EPSILON) {
    auto G0 = [](double y) {
      return y * pnorm_std(y, true, false) +
             std::exp(-0.5 * y * y - LOG_SQRT_2PI);
    };

    const double val = (G0(y_hi) - G0(y_lo)) / a;
    return signed_from_real(val);
  }

  const double r = q / a;
  const double log_C = -q * c / a + 0.5 * r * r;

  // Numerator:
  // exp(q*x_hi) Phi(y_hi)
  // - exp(q*x_lo) Phi(y_lo)
  // - C [Phi(y_hi-r) - Phi(y_lo-r)]
  signed_log num = make_signed_log(R_NegInf, 0);

  const double log_A =
    q * x_hi + pnorm_std(y_hi, true, true);
  const double log_B =
    q * x_lo + pnorm_std(y_lo, true, true);
  const double log_D =
    log_C + log_pnorm_diff(y_lo - r, y_hi - r);

  num = signed_log_add(num, make_signed_log(log_A,  1));
  num = signed_log_sub(num, make_signed_log(log_B,  1));

  if (log_D != R_NegInf) {
    num = signed_log_sub(num, make_signed_log(log_D, 1));
  }

  if (num.sign == 0 || num.log_abs == R_NegInf) {
    return make_signed_log(R_NegInf, 0);
  }

  // divide by q
  return make_signed_log(
    num.log_abs - std::log(std::abs(q)),
    num.sign * (q > 0.0 ? 1 : -1)
  );
}

inline signed_log slog_int_x_exp_pnorm(
    double q, double a, double c,
    double x_lo, double x_hi
) {
  if (x_hi <= x_lo) return make_signed_log(R_NegInf, 0);

  const double y_lo = a * x_lo + c;
  const double y_hi = a * x_hi + c;

  // q = 0:
  // K = ∫ x Phi(a*x+c) dx
  // y = a*x+c
  // x = (y-c)/a
  //
  // ∫ Phi(y) dy = y Phi(y) + phi(y)
  // ∫ y Phi(y) dy = 0.5 * ((y^2 - 1) Phi(y) + y phi(y))
  if (std::abs(q) <= RDM_Q_EPSILON) {
    auto G0 = [](double y) {
      return y * pnorm_std(y, true, false) +
             std::exp(-0.5 * y * y - LOG_SQRT_2PI);
    };

    auto G1 = [](double y) {
      const double Phi = pnorm_std(y, true, false);
      const double phi = std::exp(-0.5 * y * y - LOG_SQRT_2PI);
      return 0.5 * ((y * y - 1.0) * Phi + y * phi);
    };

    const double val =
      ((G1(y_hi) - G1(y_lo)) -
       c * (G0(y_hi) - G0(y_lo))) / (a * a);

    return signed_from_real(val);
  }

  const double abs_q = std::abs(q);
  const double q2 = q * q;

  // Boundary term:
  // [ exp(q*x) * (x/q - 1/q^2) * Phi(a*x+c) ]_lo^hi
  auto boundary_piece = [&](double x, double y) -> signed_log {
    const double coef = x / q - 1.0 / q2;
    if (std::abs(coef) <= FPM_EPSILON) {
      return make_signed_log(R_NegInf, 0);
    }

    return make_signed_log(
      q * x + std::log(std::abs(coef)) + pnorm_std(y, true, true),
      coef > 0.0 ? 1 : -1
    );
  };

  signed_log upper = boundary_piece(x_hi, y_hi);
  signed_log lower = boundary_piece(x_lo, y_lo);
  signed_log out = signed_log_sub(upper, lower);

  // L0 = ∫ exp(q*x) phi(a*x+c) dx
  // L1 = ∫ x exp(q*x) phi(a*x+c) dx
  const double r = q / a;
  const double z_lo = y_lo - r;
  const double z_hi = y_hi - r;
  const double log_C = -q * c / a + 0.5 * r * r;

  const double log_dPhi = log_pnorm_diff(z_lo, z_hi);

  if (log_dPhi != R_NegInf) {
    // L0 = C/a * DeltaPhi
    signed_log L0 = make_signed_log(
      log_C - std::log(a) + log_dPhi,
      1
    );

    // bracket for L1:
    // (r-c) DeltaPhi + phi(z_lo) - phi(z_hi)
    signed_log bracket = make_signed_log(R_NegInf, 0);

    const double rc = r - c;
    if (std::abs(rc) > FPM_EPSILON) {
      bracket = signed_log_add(
        bracket,
        make_signed_log(
          std::log(std::abs(rc)) + log_dPhi,
          rc > 0.0 ? 1 : -1
        )
      );
    }

    bracket = signed_log_add(
      bracket,
      make_signed_log(Gstar(1.0, z_lo, true), 1)
    );

    bracket = signed_log_sub(
      bracket,
      make_signed_log(Gstar(1.0, z_hi, true), 1)
    );

    signed_log L1 = make_signed_log(R_NegInf, 0);
    if (bracket.sign != 0 && bracket.log_abs != R_NegInf) {
      L1 = make_signed_log(
        log_C - 2.0 * std::log(a) + bracket.log_abs,
        bracket.sign
      );
    }

    // K = boundary - (a/q) L1 + (a/q^2) L0
    if (L1.sign != 0 && L1.log_abs != R_NegInf) {
      const int s = -1 * (q > 0.0 ? 1 : -1) * L1.sign;
      out = signed_log_add(
        out,
        make_signed_log(
          std::log(a) - std::log(abs_q) + L1.log_abs,
          s
        )
      );
    }

    out = signed_log_add(
      out,
      make_signed_log(
        std::log(a) - 2.0 * std::log(abs_q) + L0.log_abs,
        1
      )
    );
  }

  return out;
}

inline signed_log slog_int_eta_pnorm(
    double eta, double b, double a, double c,
    double x_lo, double x_hi
) {
  signed_log J = slog_int_exp_pnorm(-eta, a, c, x_lo, x_hi);

  if (J.sign == 0 || J.log_abs == R_NegInf) {
    return J;
  }

  return make_signed_log(
    eta * b + J.log_abs,
    J.sign
  );
}

inline signed_log slog_int_d_eta_pnorm(
    double eta, double b, double a, double c,
    double x_lo, double x_hi
) {
  const double q = -eta;

  signed_log J = slog_int_exp_pnorm(q, a, c, x_lo, x_hi);
  signed_log K = slog_int_x_exp_pnorm(q, a, c, x_lo, x_hi);

  signed_log bJ = make_signed_log(R_NegInf, 0);

  if (J.sign != 0 && J.log_abs != R_NegInf && b > 0.0) {
    bJ = make_signed_log(std::log(b) + J.log_abs, J.sign);
  }

  // bJ - K
  signed_log inside = signed_log_sub(bJ, K);

  if (inside.sign == 0 || inside.log_abs == R_NegInf) {
    return make_signed_log(R_NegInf, 0);
  }

  return make_signed_log(
    eta * b + inside.log_abs,
    inside.sign
  );
}

// ==========================================================================
// RDMSWTN: Racing Diffusion Model with Shifted Wald / Truncated-Normal drift
// ==========================================================================
//
// Canonical Wald form: absorbing Brownian motion with drift mu, diffusion
// sigma, starting uniformly on [b-A, b] (so b is the upper threshold).
//
// Parameter hierarchy:
//   dwald / pwald   — Wald PDF/CDF with uniform start-point range A on [0,b].
//   dswtn / pswtn   — Wald + between-trial drift variability sv (SWTN), fixed threshold.
//   drdmswtn / prdmswtn — full model: SWTN + start-point variability A.
//
// All exported functions take raw time t and internally use t_eam = t - t0.
// Erlang processes run on raw time t; evidence kernels run on t_eam.
// The `s` (diffusion) parameter is passed explicitly through the Wald/SWTN
// kernels; callers should keep drift, threshold, and start variability on their
// physical scale.
// ==========================================================================

inline double log_wald_posdrift_hit_normalizer(bool posdrift, double mu, double sigma,
                                               double b, double x_lo, double x_hi) {
  if (!posdrift || mu >= 0.0 || sigma <= 0.0) return 0.0;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  const double eta  = 2.0 * mu / (sigma * sigma);
  if (span <= FPM_EPSILON) {
    return eta * (b - x_hi);
  }
  return eta * (b - x_hi) + log_expm1_ratio(eta * span);
}

double pwald(double t, double mu, double b, double A RDM_DEFAULT_ARGUMENT(0.0), double sigma RDM_DEFAULT_ARGUMENT(1.0),
             double t0 RDM_DEFAULT_ARGUMENT(0.0), double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0), bool log_out RDM_DEFAULT_ARGUMENT(false),
             int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true),
             double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));

inline double mix_erlang12(double log_e1, double log_e2, double omega, bool log_out) {
  const double w = std::fmax(0.0, std::fmin(1.0, omega));
  if (w >= 1.0) return log_out ? log_e1 : std::exp(log_e1);
  if (w <= 0.0) return log_out ? log_e2 : std::exp(log_e2);
  const double out = log_sum_exp(std::log(w) + log_e1, std::log1p(-w) + log_e2);
  return log_out ? out : std::exp(out);
}

// log S = log(1 - F) from a log-scale CDF.  A CDF computed in log space
// keeps far more resolution near zero than a natural CDF keeps near one, so
// guess-path survivors are reconstructed from log CDFs rather than from
// clamped natural probabilities.
inline double log_surv_from_log_cdf(double log_cdf) {
  if (log_cdf == R_NegInf) return 0.0;
  if (log_cdf >= 0.0) return R_NegInf;
  return log1m_exp(log_cdf);
}

double dwald(double t, double mu, double b, double A RDM_DEFAULT_ARGUMENT(0.0), double sigma RDM_DEFAULT_ARGUMENT(1.0),
             double t0 RDM_DEFAULT_ARGUMENT(0.0), double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0), bool log_out RDM_DEFAULT_ARGUMENT(false),
             int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true),
             double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));




// --------------------------------------------------------------------------
// GBM first-passage with physical-space start-point variability.
// Process: dX_t = mu X_t dt + sigma X_t dW_t, start X0 ~ Unif(1, 1 + A).
// Boundary: absorb at X_t = b (b > 0).
// In log-space this is a Wald kernel with drift (mu - 0.5*sigma^2), diffusion
// sigma, and log-start y in [log(1), log(1 + A)] weighted by exp(y).
// The continuous-density code assumes all starts are below threshold:
// 1 + A <= b. The RDMGBM wrapper uses b = 1 + B + A, so this is equivalent
// to B >= 0.
// --------------------------------------------------------------------------
double pgbm(double t, double mu, double b, double A, double sigma,
            double t0, double lambda_g, double lambda_k, bool log_out,
            int kill_shape, bool guess, double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));

double dgbm(double t, double mu, double b, double A RDM_DEFAULT_ARGUMENT(0.0), double sigma RDM_DEFAULT_ARGUMENT(1.0),
            double t0 RDM_DEFAULT_ARGUMENT(0.0), double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0), bool log_out RDM_DEFAULT_ARGUMENT(false),
            int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));


// --------------------------------------------------------------------------
// SWTN: Shifted Wald with Truncated-Normal drift variability.
// No start-point variability (fixed threshold = threshold param).
// sv  = SD of between-trial drift distribution N(mu_drift, sv^2) [truncated at 0].
// s   = within-trial diffusion coefficient.
// c   = lower truncation for drift (default 0 for positive-drift assumption).
// --------------------------------------------------------------------------
inline double log_swtn_hit_mass(double mu, double sv, double s, double r) {
  if (r <= 0.0) return 0.0;
  
  const double s2 = s * s;
  
  // Degenerate fixed-drift case
  if (sv <= 1e-10) {
    if (mu >= 0.0) return 0.0;
    return 2.0 * r * mu / s2;
  }
  
  const double sv2 = sv * sv;
  const double alpha = 2.0 * r / s2;
  
  // log P(D >= 0)
  const double log_pos = pnorm_std(mu / sv, true, true);
  
  // log E[exp(alpha D); D < 0]
  const double z = -(mu + alpha * sv2) / sv;
  const double log_neg =
    alpha * mu +
    0.5 * alpha * alpha * sv2 +
    pnorm_std(z, true, true);
  
  return log_sum_exp(log_pos, log_neg);
}

inline double log_swtn_spv_hit_mass_full(double mu, double sv, double s,
                                         double b, double A,
                                         int n_gauss_nodes = 20) {
  if (A <= 1e-10) {
    return log_swtn_hit_mass(mu, sv, s, b);
  }

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;
  
  const double center = b - 0.5 * A;
  const double half_width = 0.5 * A;
  
  double log_sum = R_NegInf;
  
  for (int j = 0; j < n_nodes; ++j) {
    const double r_j = center + half_width * nodes[j];
    const double log_H_j = log_swtn_hit_mass(mu, sv, s, r_j);
    
    if (weights[j] > 0.0 && log_H_j != R_NegInf) {
      log_sum = log_sum_exp(log_sum, std::log(weights[j]) + log_H_j);
    }
  }
  
  // Average over [b-A,b]: 0.5 * sum weights * H(r_j)
  return log_sum - M_LN2;
}

inline double dswtn_core(double t_adj, double mu_drift, double threshold,
                         double s, double t0, double sv, double lambda,
                         double log_norm, bool log_out, int kill_shape = 1, bool guess = false,
                         bool posdrift = true);

double pswtn(double t, double mu_drift, double threshold, double s RDM_DEFAULT_ARGUMENT(1.0),
             double t0 RDM_DEFAULT_ARGUMENT(0.0), double sv RDM_DEFAULT_ARGUMENT(0.0), double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0),
             bool log_out RDM_DEFAULT_ARGUMENT(false), int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true),
             double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));

inline double dswtn_positive_drift_quad(double t, double mu_drift, double threshold,
                                        double s, double t0, double sv,
                                        double lambda_g, double lambda_k,
                                        bool log_out, int kill_shape, bool guess);
inline double pswtn_positive_drift_quad(double t, double mu_drift, double threshold,
                                        double s, double t0, double sv,
                                        double lambda_g, double lambda_k,
                                        bool log_out, int kill_shape, bool guess);
inline double drdmswtn_positive_drift_quad(double t, double mu_drift, double b, double A,
                                           double s, double t0, double sv,
                                           double lambda_g, double lambda_k,
                                           bool log_out, int kill_shape, bool guess);
inline double prdmswtn_positive_drift_quad(double t, double mu_drift, double b, double A,
                                           double s, double t0, double sv,
                                           double lambda_g, double lambda_k,
                                           bool log_out, int kill_shape, bool guess);
inline double positive_trunc_swtn_density_k0(double t, double mu_drift,
                                             double threshold, double s,
                                             double t0, double sv,
                                             bool log_out);
inline double drdmswtn_joint_A_sv_density_postrunc(
    double t_adj, double mu, double b, double A,
    double s, double sv, bool log_out
);
inline double local_combo_response_pdf(double t, double f_decision, double F_decision,
                                       double lambda_g, double lambda_k,
                                       int kill_shape, bool log_out);
inline double local_combo_response_pdf_log(double t, double log_fD, double log_FD,
                                           double lambda_g, double lambda_k,
                                           int kill_shape, bool log_out);


// --------------------------------------------------------------------------
// SWTN CDF.
// --------------------------------------------------------------------------

// Retained as the SWTN-facing name; the dispatch itself lives in gaussian.h
// so BAwL's correlated kernels can share it.
inline double norm_cdf_2d_stable(double x, double y, double rho) {
  return norm_cdf_2d_hybrid(x, y, rho);
}

// Log-space counterpart of integrate_positive_drift_quad below: accumulates
// log(w_j) + log_kernel(drift_j) with log_sum_exp, so node under/overflow in
// the natural accumulation cannot zero out a representable mixture.  Only
// entered after the natural quadrature has been rejected; log_out remains a
// pure output-scale choice at the call sites.
template <typename LogKernelFn>
inline double integrate_positive_drift_quad_log(double mu_drift, double sv,
                                                LogKernelFn&& log_kernel_fn,
                                                int n_gauss_nodes = 20) {
  if (!(sv > 1e-10) || !emc2_isfinite(sv)) {
    return log_kernel_fn(mu_drift);
  }

  const double lower_p_raw = pnorm_std(-mu_drift / sv, true, false);
  const double upper_p = std::nextafter(1.0, 0.0);
  const double lower_p = std::fmax(0.0, std::fmin(lower_p_raw, upper_p));
  const double width = upper_p - lower_p;
  if (!(width > 0.0)) return R_NegInf;

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  double log_acc = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    if (!(gl.w[j] > 0.0)) continue;
    const double p = lower_p + 0.5 * width * (gl.x[j] + 1.0);
    const double drift = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double lk = log_kernel_fn(drift);
    if (lk == R_NegInf || ISNAN(lk)) continue;
    log_acc = log_sum_exp(log_acc, std::log(gl.w[j]) + lk);
  }
  return log_acc - M_LN2;
}

inline double positive_trunc_swtn_cdf_k0(double t, double mu_drift,
                                         double threshold, double s,
                                         double t0, double sv,
                                         bool log_out) {
  const double dt = t - t0;
  if (dt <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (!(s > 1e-10) || !(sv > 1e-10) || sv < 0.0) return NA_REAL;
  if (!emc2_isfinite(dt)) return log_out ? 0.0 : 1.0;

  const double log_denom = pnorm_std(mu_drift / sv, true, true);
  if (log_denom == R_NegInf) return log_out ? R_NegInf : 0.0;

  const double sv2 = sv * sv;
  const double s2 = s * s;
  const double Q = s2 + sv2 * dt;
  const double sqrt_Q = std::sqrt(Q);
  const double denom = std::sqrt(dt * Q);
  const double rho_mag = sv * std::sqrt(dt) / sqrt_Q;
  const double z0 = -mu_drift / sv;

  const double h1 = (mu_drift * dt - threshold) / denom;
  const double phi_h1 = pnorm_std(h1, true, false);
  const double biv1 = norm_cdf_2d_stable(z0, h1, -rho_mag);
  const double term1 = std::fmax(0.0, phi_h1 - biv1);

  const double alpha = 2.0 * threshold / s2;
  const double mu_p = mu_drift + alpha * sv2;
  const double h2 = (-mu_p * dt - threshold) / denom;
  const double z0_p = -mu_p / sv;
  const double phi_h2 = pnorm_std(h2, true, false);
  const double biv2 = norm_cdf_2d_stable(z0_p, h2, rho_mag);
  const double term2_core = std::fmax(0.0, phi_h2 - biv2);

  const double log_exp_term = alpha * mu_drift + 0.5 * alpha * alpha * sv2;
  double log_num = R_NegInf;
  if (term1 > 0.0) log_num = log_sum_exp(log_num, std::log(term1));
  if (term2_core > 0.0) log_num = log_sum_exp(log_num, log_exp_term + std::log(term2_core));
  if (log_num == R_NegInf) {
    // The bivariate-normal rectangle probabilities cancelled below natural
    // resolution.  Average the stable point-Wald log CDF over the truncated
    // drift law instead (the quadrature includes the posdrift normaliser).
    const double log_cdf_q = integrate_positive_drift_quad_log(
        mu_drift, sv, [&](double drift) {
          return wald_pt_log_cdf(dt, threshold / s, drift / s);
        });
    if (!(log_cdf_q > R_NegInf) || ISNAN(log_cdf_q))
      return log_out ? R_NegInf : 0.0;
    const double lc = std::fmin(log_cdf_q, 0.0);
    return log_out ? lc : std::exp(lc);
  }

  double log_cdf = log_num - log_denom;
  if (ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  if (log_cdf > 0.0 && log_cdf < 1e-10) log_cdf = 0.0;
  if (log_cdf > 0.0) log_cdf = 0.0;
  return log_out ? log_cdf : std::exp(log_cdf);
}

inline double positive_trunc_swtn_density_k0(double t, double mu_drift,
                                             double threshold, double s,
                                             double t0, double sv,
                                             bool log_out) {
  const double dt = t - t0;
  if (dt <= 1e-10 || threshold <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (!(s > 1e-10) || !(sv > 1e-10) || sv < 0.0) return NA_REAL;

  const double v = sv * sv;
  const double s2 = s * s;
  const double den_common = dt * v + s2;
  const double log_full =
    std::log(threshold) -
    0.5 * (std::log(M_PI) + M_LN2 + 3.0 * std::log(dt) + std::log(den_common)) -
    (std::pow(threshold - mu_drift * dt, 2.0)) / (2.0 * dt * den_common);

  const double mu_new = (threshold * v + mu_drift * s2) / den_common;
  const double sigma_new = std::sqrt(s2 * v / den_common);
  const double log_tail = pnorm_std(mu_new / sigma_new, true, true);
  const double log_denom = pnorm_std(mu_drift / sv, true, true);
  const double log_pdf = log_full + log_tail - log_denom;

  if (ISNAN(log_pdf) || log_pdf == R_NegInf) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

template <typename KernelFn>
inline double integrate_positive_drift_quad(double mu_drift, double sv,
                                            KernelFn&& kernel_fn,
                                            int n_gauss_nodes = 20) {
  if (!(sv > 1e-10) || !emc2_isfinite(sv)) {
    return kernel_fn(mu_drift);
  }

  const double lower_p_raw = pnorm_std(-mu_drift / sv, true, false);
  const double upper_p = std::nextafter(1.0, 0.0);
  const double lower_p = std::fmax(0.0, std::fmin(lower_p_raw, upper_p));
  const double width = upper_p - lower_p;
  if (!(width > 0.0)) return 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;

  double acc = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double p = lower_p + 0.5 * width * (nodes[j] + 1.0);
    const double drift = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    acc += weights[j] * kernel_fn(drift);
  }

  const double result = 0.5 * acc;
  return emc2_isfinite(result) ? result : 0.0;
}

inline double dswtn_positive_drift_quad(double t, double mu_drift, double threshold,
                                        double s, double t0, double sv,
                                        double lambda_g, double lambda_k,
                                        bool log_out, int kill_shape, bool guess) {
  const double lambda = guess ? lambda_g : lambda_k;
  if (lambda <= 1e-10) {
    return positive_trunc_swtn_density_k0(t, mu_drift, threshold, s, t0, sv, log_out);
  }
  if (guess && lambda_k <= 1e-10) {
    const double log_fD =
      positive_trunc_swtn_density_k0(t, mu_drift, threshold, s, t0, sv, true);
    const double log_FD =
      positive_trunc_swtn_cdf_k0(t, mu_drift, threshold, s, t0, sv, true);
    return local_combo_response_pdf_log(t, log_fD, log_FD,
                                        lambda_g, 0.0, kill_shape, log_out);
  }

  const auto kernel = [&](double drift) {
    return dwald(t, drift, threshold, 0.0, s, t0,
                 lambda_g, lambda_k,
                 false, kill_shape, guess, false);
  };
  const double out = integrate_positive_drift_quad(mu_drift, sv, kernel);
  if (out > 0.0 && emc2_isfinite(out)) return log_out ? std::log(out) : out;
  // Natural node accumulation under/overflowed: redo it in log space.
  const double log_pdf = integrate_positive_drift_quad_log(
      mu_drift, sv, [&](double drift) {
        return dwald(t, drift, threshold, 0.0, s, t0,
                     lambda_g, lambda_k, true, kill_shape, guess, false);
      });
  if (!(log_pdf > R_NegInf) || ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return return_from_log(log_pdf, log_out);
}

inline double pswtn_positive_drift_quad(double t, double mu_drift, double threshold,
                                        double s, double t0, double sv,
                                        double lambda_g, double lambda_k,
                                        bool log_out, int kill_shape, bool guess) {
  const double lambda = guess ? lambda_g : lambda_k;
  if (lambda <= 1e-10) {
    return positive_trunc_swtn_cdf_k0(t, mu_drift, threshold, s, t0, sv, log_out);
  }
  if (guess && lambda_k <= 1e-10) {
    auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    const double fw = std::max(0.0, std::min(1.0,
      positive_trunc_swtn_cdf_k0(t, mu_drift, threshold, s, t0, sv, false)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t, lambda_g, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  const auto kernel = [&](double drift) {
    return pwald(t, drift, threshold, 0.0, s, t0,
                 lambda_g, lambda_k,
                 false, kill_shape, guess, false);
  };
  double out = integrate_positive_drift_quad(mu_drift, sv, kernel);
  out = std::fmax(0.0, std::fmin(1.0, out));
  if (out > 0.0) return log_out ? std::log(out) : out;
  // Natural node accumulation underflowed: redo it in log space.
  const double log_cdf = integrate_positive_drift_quad_log(
      mu_drift, sv, [&](double drift) {
        return pwald(t, drift, threshold, 0.0, s, t0,
                     lambda_g, lambda_k, true, kill_shape, guess, false);
      });
  if (!(log_cdf > R_NegInf) || ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  return return_from_log(std::fmin(log_cdf, 0.0), log_out);
}

inline double drdmswtn_joint_A_sv_density_fullgauss(
    double t_adj, double mu, double b, double A,
    double s, double sv, bool posdrift,
    int n_gauss_nodes, bool log_out
) {
  const double r_lo = b - A;
  const double r_hi = b;
  const double s2 = s * s;
  const double sv2 = sv * sv;
  
  const double m = mu * t_adj;
  const double V = t_adj * (s2 + sv2 * t_adj);
  const double S = std::sqrt(V);
  
  const double z_lo = (r_lo - m) / S;
  const double z_hi = (r_hi - m) / S;
  
  const double log_hi = pnorm_std(z_hi, true, true);
  const double log_lo = pnorm_std(z_lo, true, true);
  const double Phi_diff = (z_hi > z_lo) ? std::exp(log_hi) * (-std::expm1(log_lo - log_hi)) : 0.0;
  
  const double phi_lo = std::exp(-0.5 * z_lo * z_lo - LOG_SQRT_2PI);
  const double phi_hi = std::exp(-0.5 * z_hi * z_hi - LOG_SQRT_2PI);
  
  double dens = (m * Phi_diff + S * (phi_lo - phi_hi)) / (A * t_adj);
  dens = std::fmax(0.0, dens);
  
  if (posdrift) {
    const double log_Hbar =
      log_swtn_spv_hit_mass_full(mu, sv, s, b, A, n_gauss_nodes);
    dens *= std::exp(-log_Hbar);
  }
  
  if (!emc2_isfinite(dens) || dens <= 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  
  return log_out ? std::log(dens) : dens;
}

inline double drdmswtn_joint_A_sv_density_postrunc(
    double t_adj, double mu, double b, double A,
    double s, double sv, bool log_out
) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (A <= 1e-10) {
    return positive_trunc_swtn_density_k0(t_adj, mu, b, s, 0.0, sv, log_out);
  }
  if (!(s > 1e-10) || !(sv > 1e-10) || !(A > 0.0)) return NA_REAL;

  const double r_lo = b - A;
  const double r_hi = b;
  if (!(r_hi > r_lo)) return log_out ? R_NegInf : 0.0;

  const double s2 = s * s;
  const double sv2 = sv * sv;
  const double mx = mu * t_adj;
  const double sD = std::sqrt(t_adj * (s2 + sv2 * t_adj));
  const double rho = sv * std::sqrt(t_adj) / std::sqrt(s2 + sv2 * t_adj);
  const double gamma = -mu / sv;
  const double z_plus = pnorm_std(mu / sv, true, false);
  if (!(z_plus > 0.0)) return log_out ? R_NegInf : 0.0;

  const double a_lo = (r_lo - mx) / sD;
  const double a_hi = (r_hi - mx) / sD;

  const double phi_lo = std::exp(-0.5 * a_lo * a_lo - LOG_SQRT_2PI);
  const double phi_hi = std::exp(-0.5 * a_hi * a_hi - LOG_SQRT_2PI);
  const double phi_gamma = std::exp(-0.5 * gamma * gamma - LOG_SQRT_2PI);

  const double sqrt_one_minus_rho2 = std::sqrt(std::fmax(1e-15, 1.0 - rho * rho));
  const auto G = [&](double y, double phi_y) {
    return -phi_y * pnorm_std((rho * y - gamma) / sqrt_one_minus_rho2, true, false) +
      rho * phi_gamma * pnorm_std((y - rho * gamma) / sqrt_one_minus_rho2, true, false);
  };

  const double p_rect =
    std::fmax(0.0, pnorm_std(a_hi, true, false) - pnorm_std(a_lo, true, false) -
      (norm_cdf_2d_stable(a_hi, gamma, rho) - norm_cdf_2d_stable(a_lo, gamma, rho)));
  const double m_rect = G(a_hi, phi_hi) - G(a_lo, phi_lo);

  double dens = (mx * p_rect + sD * m_rect) / (A * z_plus * t_adj);
  dens = std::fmax(0.0, dens);

  if (!emc2_isfinite(dens) || dens <= 0.0) {
    // The rectangle probabilities cancelled below natural resolution while
    // the true density is positive: rebuild the defective fixed-drift SPV
    // density mixture in log space (the quadrature includes the posdrift
    // normalisation over the truncated drift law).
    const double log_pdf = integrate_positive_drift_quad_log(
        mu, sv, [&](double drift) {
          return dwald_k0_log(t_adj, b, drift, A, s);
        });
    if (!(log_pdf > R_NegInf) || ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
    return return_from_log(log_pdf, log_out);
  }
  return log_out ? std::log(dens) : dens;
}

inline double prdmswtn_joint_A_sv_cdf_postrunc(
    double t, double mu, double b, double A,
    double s, double t0, double sv,
    int n_gauss_nodes, bool log_out
) {
  const double dt = t - t0;
  if (dt <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (A <= 1e-10) {
    return positive_trunc_swtn_cdf_k0(t, mu, b, s, t0, sv, log_out);
  }

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;

  const double center = b - 0.5 * A;
  const double half_width = 0.5 * A;
  double acc = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double thresh_j = center + half_width * nodes[j];
    acc += weights[j] *
      positive_trunc_swtn_cdf_k0(t, mu, thresh_j, s, t0, sv, false);
  }

  const double cdf = std::fmax(0.0, std::fmin(1.0, 0.5 * acc));
  if (cdf > 0.0) return log_out ? std::log(cdf) : cdf;

  // Every node underflowed on the natural scale: repeat the threshold
  // average in log space.
  double log_acc = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    if (!(weights[j] > 0.0)) continue;
    const double thresh_j = center + half_width * nodes[j];
    const double lc = positive_trunc_swtn_cdf_k0(t, mu, thresh_j, s, t0, sv, true);
    if (lc == R_NegInf || ISNAN(lc)) continue;
    log_acc = log_sum_exp(log_acc, std::log(weights[j]) + lc);
  }
  const double log_cdf = log_acc - M_LN2;
  if (!(log_cdf > R_NegInf) || ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  return return_from_log(std::fmin(log_cdf, 0.0), log_out);
}

inline double drdmswtn_positive_drift_quad(double t, double mu_drift, double b, double A,
                                           double s, double t0, double sv,
                                           double lambda_g, double lambda_k,
                                           bool log_out, int kill_shape, bool guess) {
  if (A < 1e-7) {
    return dswtn_positive_drift_quad(t, mu_drift, b, s, t0, sv,
                                     lambda_g, lambda_k, log_out, kill_shape, guess);
  }

  const auto kernel = [&](double drift) {
    // We already integrate over the positive-truncated drift law here, so the
    // fixed-drift kernel must stay defective and unnormalised.
    return dwald(t, drift, b, A, s, t0,
                 lambda_g, lambda_k,
                 false, kill_shape, guess, false);
  };
  const double out = integrate_positive_drift_quad(mu_drift, sv, kernel);
  if (out > 0.0 && emc2_isfinite(out)) return log_out ? std::log(out) : out;
  // Natural node accumulation under/overflowed: redo it in log space.
  const double log_pdf = integrate_positive_drift_quad_log(
      mu_drift, sv, [&](double drift) {
        return dwald(t, drift, b, A, s, t0,
                     lambda_g, lambda_k, true, kill_shape, guess, false);
      });
  if (!(log_pdf > R_NegInf) || ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return return_from_log(log_pdf, log_out);
}

inline double prdmswtn_positive_drift_quad(double t, double mu_drift, double b, double A,
                                           double s, double t0, double sv,
                                           double lambda_g, double lambda_k,
                                           bool log_out, int kill_shape, bool guess) {
  if (A < 1e-7) {
    return pswtn_positive_drift_quad(t, mu_drift, b, s, t0, sv,
                                     lambda_g, lambda_k, log_out, kill_shape, guess);
  }

  const auto kernel = [&](double drift) {
    // We already integrate over the positive-truncated drift law here, so the
    // fixed-drift kernel must stay defective and unnormalised.
    return pwald(t, drift, b, A, s, t0,
                 lambda_g, lambda_k,
                 false, kill_shape, guess, false);
  };
  double out = integrate_positive_drift_quad(mu_drift, sv, kernel);
  out = std::fmax(0.0, std::fmin(1.0, out));
  if (out > 0.0) return log_out ? std::log(out) : out;
  // Natural node accumulation underflowed: redo it in log space.
  const double log_cdf = integrate_positive_drift_quad_log(
      mu_drift, sv, [&](double drift) {
        return pwald(t, drift, b, A, s, t0,
                     lambda_g, lambda_k, true, kill_shape, guess, false);
      });
  if (!(log_cdf > R_NegInf) || ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  return return_from_log(std::fmin(log_cdf, 0.0), log_out);
}

double drdmswtn(double t, double mu_drift, double b, double A,
                double s RDM_DEFAULT_ARGUMENT(1.0), double t0 RDM_DEFAULT_ARGUMENT(0.0), double sv RDM_DEFAULT_ARGUMENT(0.0),
                double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0),
                int n_gauss_nodes RDM_DEFAULT_ARGUMENT(20), bool log_out RDM_DEFAULT_ARGUMENT(false),
                int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true),
                double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));
double prdmswtn(double t, double mu_drift, double b, double A,
                double s RDM_DEFAULT_ARGUMENT(1.0), double t0 RDM_DEFAULT_ARGUMENT(0.0), double sv RDM_DEFAULT_ARGUMENT(0.0),
                double lambda_g RDM_DEFAULT_ARGUMENT(0.0), double lambda_k RDM_DEFAULT_ARGUMENT(0.0),
                int n_gauss_nodes RDM_DEFAULT_ARGUMENT(20), bool log_out RDM_DEFAULT_ARGUMENT(false),
                int kill_shape RDM_DEFAULT_ARGUMENT(1), bool guess RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true),
                double erlang_omega RDM_DEFAULT_ARGUMENT(1.0));

template <typename DensityFn>
inline double integrate_density_gl20_finite(double t_upper, DensityFn&& density_fn) {
  if (!(t_upper > 0.0)) return 0.0;
  const GLRule& gl = gl_get_rule(20);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;
  double acc = 0.0;
  for (int j = 0; j < static_cast<int>(nodes.size()); ++j) {
    const double u = 0.5 * t_upper * (nodes[j] + 1.0);
    acc += weights[j] * density_fn(u);
  }
  return 0.5 * t_upper * acc;
}

template <typename DensityFn>
inline double integrate_density_gl20_infinite(double rate_scale, DensityFn&& density_fn) {
  const double scale = std::fmax(rate_scale, 1e-8);
  const GLRule& gl = gl_get_rule(20);
  const std::vector<double>& nodes = gl.x;
  const std::vector<double>& weights = gl.w;
  double acc = 0.0;
  for (int j = 0; j < static_cast<int>(nodes.size()); ++j) {
    const double q = 0.5 * (nodes[j] + 1.0);
    const double qq = std::fmin(1.0 - 1e-12, std::fmax(1e-15, q));
    const double t = -std::log1p(-qq) / scale;
    const double jac = 1.0 / (scale * (1.0 - qq));
    acc += weights[j] * density_fn(t) * jac;
  }
  return 0.5 * acc;
}

template <typename DensityFn>
inline double integrate_density_adaptive_finite(double t_upper, DensityFn&& density_fn) {
  if (!(t_upper > 0.0)) return 0.0;

  using Fn = std::decay_t<DensityFn>;
  Fn fn = std::forward<DensityFn>(density_fn);
  struct Adapter {
    const Fn* fn;
  } adapter{&fn};

  gsl_function F;
  F.function = +[](double x, void* p) -> double {
    const auto* adapter_ptr = static_cast<const Adapter*>(p);
    return (*(adapter_ptr->fn))(x);
  };
  F.params = &adapter;

  GslIntegrationControls ctl = default_gsl_controls();
  ctl.try_qng_first_finite = true;
  ctl.qag_key = GSL_INTEG_GAUSS21;
  ctl.rel_tol = 1e-6;

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws, ctl.retry_limit);

  double result = 0.0;
  double err = 0.0;
  int status = GSL_EFAILED;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  size_t neval = 0;
  status = gsl_integration_qng(&F, 0.0, t_upper, ctl.abs_tol, ctl.rel_tol,
                               &result, &err, &neval);
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qag(&F, 0.0, t_upper,
                                 ctl.abs_tol, ctl.rel_tol,
                                 ctl.limit, ctl.qag_key,
                                 workspace, &result, &err);
  }
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qags(&F, 0.0, t_upper,
                                  ctl.retry_abs_tol, ctl.retry_rel_tol,
                                  ctl.retry_limit, workspace, &result, &err);
  }

  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS || !emc2_isfinite(result) || result < 0.0) {
    return integrate_density_gl20_finite(t_upper, fn);
  }
  return result;
}

template <typename DensityFn>
inline double integrate_density_adaptive_infinite(double rate_scale, DensityFn&& density_fn) {
  using Fn = std::decay_t<DensityFn>;
  Fn fn = std::forward<DensityFn>(density_fn);
  struct Adapter {
    const Fn* fn;
  } adapter{&fn};

  gsl_function F;
  F.function = +[](double x, void* p) -> double {
    const auto* adapter_ptr = static_cast<const Adapter*>(p);
    return (*(adapter_ptr->fn))(x);
  };
  F.params = &adapter;

  GslIntegrationControls ctl = default_gsl_controls();
  ctl.rel_tol = 1e-6;
  if (rate_scale > 1e-8 && emc2_isfinite(rate_scale)) {
    ctl.abs_tol = std::min(ctl.abs_tol, 1e-10 / rate_scale);
  }

  static thread_local GslWorkspacePtr ws(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* workspace = ensure_gsl_workspace(ws, ctl.retry_limit);

  double result = 0.0;
  double err = 0.0;
  int status = GSL_EFAILED;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

  status = gsl_integration_qagiu(&F, 0.0, ctl.abs_tol, ctl.rel_tol,
                                 ctl.limit, workspace, &result, &err);
  if (status != GSL_SUCCESS || !emc2_isfinite(result)) {
    status = gsl_integration_qagiu(&F, 0.0, ctl.retry_abs_tol, ctl.retry_rel_tol,
                                   ctl.retry_limit, workspace, &result, &err);
  }

  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS || !emc2_isfinite(result) || result < 0.0) {
    return integrate_density_gl20_infinite(rate_scale, fn);
  }
  return result;
}

inline double local_combo_response_pdf(double t, double f_decision, double F_decision,
                                       double lambda_g, double lambda_k,
                                       int kill_shape, bool log_out) {
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;
  const bool use_guess = lambda_g > 0.0;
  const bool use_kill  = lambda_k > 0.0;
  if (!use_guess && !use_kill) {
    if (!(f_decision > 0.0) || !emc2_isfinite(f_decision)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(f_decision) : f_decision;
  }

  const double log_fD = (f_decision > 0.0 && emc2_isfinite(f_decision)) ? std::log(f_decision) : R_NegInf;
  const double log_sG = use_guess ? erlang_log_surv(t, lambda_g, kill_shape) : 0.0;
  const double log_sK = use_kill  ? erlang_log_surv(t, lambda_k, kill_shape) : 0.0;
  const double log_hit = log_fD + log_sG + log_sK;
  if (!use_guess) return log_out ? log_hit : std::exp(log_hit);

  const double F_clamped = std::fmax(0.0, std::fmin(1.0, F_decision));
  const double log_sD = std::log1p(-F_clamped);
  const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape);
  const double log_guess = log_fG + log_sK + log_sD;
  const double log_pdf = log_sum_exp(log_hit, log_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

// Log-scale variant: consumes a log decision density and a log decision CDF
// so far-tail components survive intact (no natural round trip before the
// mixture).  Natural callers keep local_combo_response_pdf above.
inline double local_combo_response_pdf_log(double t, double log_fD, double log_FD,
                                           double lambda_g, double lambda_k,
                                           int kill_shape, bool log_out) {
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;
  const bool use_guess = lambda_g > 0.0;
  const bool use_kill  = lambda_k > 0.0;
  if (ISNAN(log_fD)) log_fD = R_NegInf;
  if (!use_guess && !use_kill) {
    return log_out ? log_fD : std::exp(log_fD);
  }

  const double log_sG = use_guess ? erlang_log_surv(t, lambda_g, kill_shape) : 0.0;
  const double log_sK = use_kill  ? erlang_log_surv(t, lambda_k, kill_shape) : 0.0;
  const double log_hit = log_fD + log_sG + log_sK;
  if (!use_guess) return log_out ? log_hit : std::exp(log_hit);

  const double log_sD = log_surv_from_log_cdf(log_FD);
  const double log_fG = erlang_log_pdf(t, lambda_g, kill_shape);
  const double log_guess = log_fG + log_sK + log_sD;
  const double log_pdf = log_sum_exp(log_hit, log_guess);
  if (ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double local_combo_response_cdf_exp(double t, double F_decision,
                                           double F_decision_killed,
                                           double lambda_g, double lambda_k,
                                           bool log_out) {
  const double rate = lambda_g + lambda_k;
  if (!(t > 0.0) || !(rate > 0.0)) return log_out ? R_NegInf : 0.0;

  const double D = std::fmax(0.0, std::fmin(1.0, F_decision_killed));
  const double S = 1.0 - std::fmax(0.0, std::fmin(1.0, F_decision));
  const double tail = emc2_isfinite(t) ? std::exp(-rate * t) * S : 0.0;
  const double A0 = (1.0 - tail - D) / rate;
  double out = D + lambda_g * std::fmax(0.0, A0);
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline void exp_poly_moments_0_1_2(double upper, double rate,
                                   double& e0, double& e1, double& e2) {
  if (!(upper > 0.0)) {
    e0 = e1 = e2 = 0.0;
    return;
  }
  if (!emc2_isfinite(upper)) {
    e0 = 1.0 / rate;
    e1 = 1.0 / (rate * rate);
    e2 = 2.0 / (rate * rate * rate);
    return;
  }

  const double ru = rate * upper;
  if (std::abs(ru) < 1e-5) {
    const double u2 = upper * upper;
    const double u3 = u2 * upper;
    const double u4 = u3 * upper;
    const double u5 = u4 * upper;
    e0 = upper - 0.5 * rate * u2 + (rate * rate) * u3 / 6.0;
    e1 = 0.5 * u2 - rate * u3 / 3.0 + (rate * rate) * u4 / 8.0;
    e2 = u3 / 3.0 - rate * u4 / 4.0 + (rate * rate) * u5 / 10.0;
    return;
  }

  const double er = std::exp(-ru);
  const double r2 = rate * rate;
  const double r3 = r2 * rate;
  e0 = -std::expm1(-ru) / rate;
  e1 = (1.0 - er * (1.0 + ru)) / r2;
  e2 = (2.0 - er * (ru * ru + 2.0 * ru + 2.0)) / r3;
}

struct tilted_wald_moments_0_1_2 {
  bool ok;
  double F;
  double M0;
  double M1;
  double M2;
};

inline tilted_wald_moments_0_1_2 point_wald_tilted_moments_0_1_2(
    double upper, double distance, double drift, double sigma, double rate,
    bool posdrift = true) {
  tilted_wald_moments_0_1_2 out{false, 0.0, 0.0, 0.0, 0.0};
  if (!(distance > 0.0) || !(sigma > 0.0) || !(rate > 0.0)) return out;
  if (posdrift && drift <= 0.0) {
    out.ok = true;
    return out;
  }

  const double sig2 = sigma * sigma;
  if (!emc2_isfinite(upper)) {
    out.F = std::exp(log_wald_posdrift_hit_normalizer(
      posdrift, drift, sigma, distance, distance, distance));
  } else if (upper > 0.0) {
    out.F = std::fmax(0.0, std::fmin(1.0,
      pwald(upper, drift, distance, 0.0, sigma, 0.0,
            0.0, 0.0, false, 1, false, posdrift)));
  } else {
    out.ok = true;
    return out;
  }

  const double q = std::sqrt(drift * drift + 2.0 * sig2 * rate);
  const double beta = emc2_isfinite(upper) ? q * std::sqrt(upper) / sigma : 1.0;
  if (!(q > 1e-12) || (emc2_isfinite(upper) && std::abs(beta) < 1e-7)) {
    return out;
  }

  const double L = std::exp(distance * (drift - q) / sig2);
  if (!emc2_isfinite(upper)) {
    out.M0 = L;
    out.M1 = L * distance / q;
    out.M2 = L * (distance * distance / (q * q) +
                  distance * sig2 / (q * q * q));
    out.ok = true;
    return out;
  }

  const double sqrt_u = std::sqrt(upper);
  const double st = sigma * sqrt_u;
  const double d_minus = (q * upper - distance) / st;
  const double phi_minus = std::exp(-0.5 * d_minus * d_minus - LOG_SQRT_2PI);

  out.M0 = std::fmax(0.0, std::fmin(1.0,
    pwald(upper, drift, distance, 0.0, sigma, 0.0,
          0.0, rate, false, 1, false, posdrift)));
  const double erlang2_cdf = std::fmax(0.0, std::fmin(1.0,
    pwald(upper, drift, distance, 0.0, sigma, 0.0,
          0.0, rate, false, 2, false, posdrift)));
  out.M1 = std::fmax(0.0, (erlang2_cdf - out.M0) / rate);
  out.M2 = (distance * distance / (q * q)) * out.M0 +
           (sig2 / (q * q)) * out.M1 -
           L * (2.0 * distance * sigma * sqrt_u / (q * q)) * phi_minus;

  out.M0 = std::fmax(0.0, out.M0);
  out.M1 = std::fmax(0.0, out.M1);
  out.M2 = std::fmax(0.0, out.M2);
  out.ok = true;
  return out;
}

inline double primitive_K1(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  const double h = z + eta * S;
  const double E = std::exp(eta * c + 0.5 * eta * eta * S * S);
  const double m = c + eta * S * S;
  return -S * E * (m * pnorm_std(h, true, false) +
                   S * std::exp(-0.5 * h * h - LOG_SQRT_2PI));
}

inline double primitive_K2(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  const double h = z + eta * S;
  const double E = std::exp(eta * c + 0.5 * eta * eta * S * S);
  const double m = c + eta * S * S;
  const double phi_h = std::exp(-0.5 * h * h - LOG_SQRT_2PI);
  return -S * E * ((m * m + S * S) * pnorm_std(h, true, false) +
                   S * (2.0 * m - S * h) * phi_h);
}

inline double primitive_K3(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  const double h = z + eta * S;
  const double E = std::exp(eta * c + 0.5 * eta * eta * S * S);
  const double m = c + eta * S * S;
  const double phi_h = std::exp(-0.5 * h * h - LOG_SQRT_2PI);
  return -S * E * ((m * m * m + 3.0 * m * S * S) * pnorm_std(h, true, false) +
                   S * (3.0 * m * m - 3.0 * m * S * h +
                        S * S * (h * h + 2.0)) * phi_h);
}

inline double primitive_H0(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  const double h = z + eta * S;
  const double E = std::exp(eta * c + 0.5 * eta * eta * S * S);
  return (std::exp(eta * a) * pnorm_std(z, true, false) -
          E * pnorm_std(h, true, false)) / eta;
}

inline double primitive_H1(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  if (std::abs(eta) <= 1e-8) {
    return 0.5 * a * a * pnorm_std(z, true, false) +
           primitive_K2(0.0, c, a, S) / (2.0 * S);
  }
  return a * std::exp(eta * a) * pnorm_std(z, true, false) / eta -
         primitive_H0(eta, c, a, S) / eta +
         primitive_K1(eta, c, a, S) / (eta * S);
}

inline double primitive_H2(double eta, double c, double a, double S) {
  const double z = (c - a) / S;
  if (std::abs(eta) <= 1e-8) {
    return (a * a * a / 3.0) * pnorm_std(z, true, false) +
           primitive_K3(0.0, c, a, S) / (3.0 * S);
  }
  return a * a * std::exp(eta * a) * pnorm_std(z, true, false) / eta -
         2.0 * primitive_H1(eta, c, a, S) / eta +
         primitive_K2(eta, c, a, S) / (eta * S);
}

inline double primitive_exp_a1(double eta, double a) {
  if (std::abs(eta) <= 1e-8) return 0.5 * a * a;
  return std::exp(eta * a) * (a / eta - 1.0 / (eta * eta));
}

inline double primitive_exp_a2(double eta, double a) {
  if (std::abs(eta) <= 1e-8) return a * a * a / 3.0;
  return std::exp(eta * a) *
         (a * a / eta - 2.0 * a / (eta * eta) + 2.0 / (eta * eta * eta));
}

inline double primitive_spv_m2_finite(double endpoint, double upper,
                                      double drift, double sigma, double rate,
                                      bool gbm_weighted) {
  const double sig2 = sigma * sigma;
  const double q = std::sqrt(drift * drift + 2.0 * sig2 * rate);
  if (!(q > 0.0) || !(upper > 0.0)) return NA_REAL;
  const double S = sigma * std::sqrt(upper);
  const double eta_m = (drift - q) / sig2 - (gbm_weighted ? 1.0 : 0.0);
  const double eta_p = (drift + q) / sig2 - (gbm_weighted ? 1.0 : 0.0);
  const double c_m = q * upper;
  const double c_p = -q * upper;
  const double q2 = q * q;
  const double q3 = q2 * q;
  return (primitive_H2(eta_m, c_m, endpoint, S) +
          primitive_H2(eta_p, c_p, endpoint, S)) / q2 +
         (sig2 / q3) *
          (primitive_H1(eta_m, c_m, endpoint, S) -
           primitive_H1(eta_p, c_p, endpoint, S)) -
         (2.0 * S / q2) * primitive_K1(eta_m, c_m, endpoint, S);
}

inline double primitive_spv_m2_infinite(double endpoint,
                                        double drift, double sigma, double rate,
                                        bool gbm_weighted) {
  const double sig2 = sigma * sigma;
  const double q = std::sqrt(drift * drift + 2.0 * sig2 * rate);
  if (!(q > 0.0)) return NA_REAL;
  const double eta = (drift - q) / sig2 - (gbm_weighted ? 1.0 : 0.0);
  const double q2 = q * q;
  const double q3 = q2 * q;
  return primitive_exp_a2(eta, endpoint) / q2 +
         sig2 * primitive_exp_a1(eta, endpoint) / q3;
}

inline double spv_wald_m2(double upper, double b, double A,
                          double drift, double sigma, double rate) {
  if (A <= 1e-8) {
    return point_wald_tilted_moments_0_1_2(
      upper, b, drift, sigma, rate, false).M2;
  }
  const double lo = b - A;
  if (!(lo > 0.0)) return NA_REAL;
  const double hi = b;
  const double p_hi = emc2_isfinite(upper)
    ? primitive_spv_m2_finite(hi, upper, drift, sigma, rate, false)
    : primitive_spv_m2_infinite(hi, drift, sigma, rate, false);
  const double p_lo = emc2_isfinite(upper)
    ? primitive_spv_m2_finite(lo, upper, drift, sigma, rate, false)
    : primitive_spv_m2_infinite(lo, drift, sigma, rate, false);
  return (p_hi - p_lo) / A;
}

inline double gbm_spv_m2(double upper, double b, double A,
                         double drift, double sigma, double rate) {
  if (A <= 1e-8) {
    return point_wald_tilted_moments_0_1_2(
      upper, std::log(b), drift, sigma, rate, false).M2;
  }
  if (!(b > 1.0 + A)) return NA_REAL;
  const double lo = std::log(b / (1.0 + A));
  const double hi = std::log(b);
  const double p_hi = emc2_isfinite(upper)
    ? primitive_spv_m2_finite(hi, upper, drift, sigma, rate, true)
    : primitive_spv_m2_infinite(hi, drift, sigma, rate, true);
  const double p_lo = emc2_isfinite(upper)
    ? primitive_spv_m2_finite(lo, upper, drift, sigma, rate, true)
    : primitive_spv_m2_infinite(lo, drift, sigma, rate, true);
  return (b / A) * (p_hi - p_lo);
}

inline tilted_wald_moments_0_1_2 spv_wald_tilted_moments_0_1_2(
    double upper, double b, double A, double drift, double sigma, double rate,
    bool posdrift = true) {
  tilted_wald_moments_0_1_2 out{false, 0.0, 0.0, 0.0, 0.0};
  if (posdrift && drift <= 0.0) {
    out.ok = true;
    return out;
  }
  if (A <= 1e-8) {
    return point_wald_tilted_moments_0_1_2(
      upper, b, drift, sigma, rate, posdrift);
  }
  
  if (!emc2_isfinite(upper)) {
    out.F = std::exp(log_wald_posdrift_hit_normalizer(
      posdrift, drift, sigma, b, b - A, b));
    out.M0 = std::exp(pwald(
      R_PosInf, drift, b, A, sigma, 0.0, 0.0, rate, true, 1, false, posdrift));
    const double erlang2_cdf = std::exp(pwald(
      R_PosInf, drift, b, A, sigma, 0.0, 0.0, rate, true, 2, false, posdrift));
    out.M1 = std::fmax(0.0, (erlang2_cdf - out.M0) / rate);
    out.M2 = spv_wald_m2(upper, b, A, drift, sigma, rate);
    out.ok = emc2_isfinite(out.M2) && out.M2 >= 0.0;
    return out;
  }

  if (upper <= 0.0) {
    out.ok = true;
    return out;
  }
  if (!(b - A > 0.0)) return out;
  out.F = std::fmax(0.0, std::fmin(1.0,
    pwald(upper, drift, b, A, sigma, 0.0, 0.0, 0.0, false, 1, false, posdrift)));
  out.M0 = std::fmax(0.0, std::fmin(1.0,
    pwald(upper, drift, b, A, sigma, 0.0, 0.0, rate, false, 1, false, posdrift)));
  const double erlang2_cdf = std::fmax(0.0, std::fmin(1.0,
    pwald(upper, drift, b, A, sigma, 0.0, 0.0, rate, false, 2, false, posdrift)));
  out.M1 = std::fmax(0.0, (erlang2_cdf - out.M0) / rate);
  out.M2 = spv_wald_m2(upper, b, A, drift, sigma, rate);
  out.ok = emc2_isfinite(out.M2) && out.M2 >= 0.0;
  return out;
}

inline tilted_wald_moments_0_1_2 gbm_tilted_moments_0_1_2(
    double upper, double mu, double b, double A, double sigma, double rate) {
  tilted_wald_moments_0_1_2 out{false, 0.0, 0.0, 0.0, 0.0};
  const double drift = mu - 0.5 * sigma * sigma;
  if (A <= 1e-8) {
    return point_wald_tilted_moments_0_1_2(
      upper, std::log(b), drift, sigma, rate, false);
  }
  if (emc2_isfinite(upper) && upper <= 0.0) {
    out.ok = true;
    return out;
  }
  if (!(b > 1.0 + A)) return out;
  out.F = std::fmax(0.0, std::fmin(1.0,
    pgbm(upper, mu, b, A, sigma, 0.0, 0.0, 0.0, false, 1, false)));
  out.M0 = std::fmax(0.0, std::fmin(1.0,
    pgbm(upper, mu, b, A, sigma, 0.0, 0.0, rate, false, 1, false)));
  const double erlang2_cdf = std::fmax(0.0, std::fmin(1.0,
    pgbm(upper, mu, b, A, sigma, 0.0, 0.0, rate, false, 2, false)));
  out.M1 = std::fmax(0.0, (erlang2_cdf - out.M0) / rate);
  out.M2 = gbm_spv_m2(upper, b, A, drift, sigma, rate);
  out.ok = emc2_isfinite(out.M2) && out.M2 >= 0.0;
  return out;
}

inline double local_combo_cdf_erlang2_from_moments(
    double t, double t0, double lambda_g, double lambda_k,
    const tilted_wald_moments_0_1_2& m, bool log_out) {
  const double rate = lambda_g + lambda_k;
  double pre0, pre1, pre2;
  const double pre_upper = emc2_isfinite(t) ? std::fmin(t, std::fmax(0.0, t0))
                                            : std::fmax(0.0, t0);
  exp_poly_moments_0_1_2(pre_upper, rate, pre0, pre1, pre2);
  const double pre_guess = lambda_g * lambda_g * (pre1 + lambda_k * pre2);

  if (emc2_isfinite(t) && t <= t0) {
    const double out = std::fmax(0.0, std::fmin(1.0, pre_guess));
    return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
  }

  const double gk = lambda_g * lambda_k;
  const double t02 = t0 * t0;
  const double a0 = 1.0 + rate * t0 + gk * t02;
  const double a1 = rate + 2.0 * gk * t0;
  const double a2 = gk;
  const double exp_t0 = std::exp(-rate * t0);
  const double decision = exp_t0 * (a0 * m.M0 + a1 * m.M1 + a2 * m.M2);

  const double upper = emc2_isfinite(t) ? t - t0 : R_PosInf;
  double p0, p1, p2;
  exp_poly_moments_0_1_2(upper, rate, p0, p1, p2);

  const double r2 = rate * rate;
  const double r3 = r2 * rate;
  const double J0 = p0 * m.F - (m.F - m.M0) / rate;
  const double J1 = p1 * m.F - (m.F - m.M0 - rate * m.M1) / r2;
  const double J2 = p2 * m.F -
    (2.0 * m.F - 2.0 * m.M0 - 2.0 * rate * m.M1 - r2 * m.M2) / r3;

  const double c0 = t0 + lambda_k * t02;
  const double c1 = 1.0 + 2.0 * lambda_k * t0;
  const double c2 = lambda_k;
  const double post_surv_int =
    c0 * std::fmax(0.0, p0 - J0) +
    c1 * std::fmax(0.0, p1 - J1) +
    c2 * std::fmax(0.0, p2 - J2);
  const double post_guess = lambda_g * lambda_g * exp_t0 * post_surv_int;

  double out = pre_guess + decision + post_guess;
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline double drdmswtn_local_combo(double t, double mu_drift, double b, double A,
                                   double s, double t0, double sv,
                                   double lambda_g, double lambda_k,
                                   int n_gauss_nodes = 20,
                                   bool log_out = false, int kill_shape = 1, bool posdrift = true,
                                   double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = drdmswtn_local_combo(t, mu_drift, b, A, s, t0, sv, lambda_g, lambda_k, n_gauss_nodes, true, 1, posdrift, 1.0);
    const double e2 = drdmswtn_local_combo(t, mu_drift, b, A, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, n_gauss_nodes, true, 2, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // EAM density/CDF computed at EAM time; core functions now take raw t and
  // t0.  Both enter the mixture on the log scale so far-tail components are
  // not squeezed through a natural round trip first.
  const double log_fD = drdmswtn(t, mu_drift, b, A, s, t0, sv, 0.0, 0.0, n_gauss_nodes, true, kill_shape, false, posdrift, 1.0);
  const double log_FD = prdmswtn(t, mu_drift, b, A, s, t0, sv, 0.0, 0.0, n_gauss_nodes, true, kill_shape, false, posdrift, 1.0);
  return local_combo_response_pdf_log(t, log_fD, log_FD,
                                      lambda_g, lambda_k, kill_shape, log_out);
}

inline double prdmswtn_local_combo(double t, double mu_drift, double b, double A,
                                   double s, double t0, double sv,
                                   double lambda_g, double lambda_k,
                                   int n_gauss_nodes = 20,
                                   bool log_out = false, int kill_shape = 1, bool posdrift = true,
                                   double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = prdmswtn_local_combo(t, mu_drift, b, A, s, t0, sv, lambda_g, lambda_k, n_gauss_nodes, true, 1, posdrift, 1.0);
    const double e2 = prdmswtn_local_combo(t, mu_drift, b, A, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, n_gauss_nodes, true, 2, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;

  if (kill_shape <= 1) {
    const double rate = lambda_g + lambda_k;
    const double F_decision = prdmswtn(t, mu_drift, b, A, s, t0, sv,
                                       0.0, 0.0, n_gauss_nodes,
                                       false, 1, false, posdrift, 1.0);
    const double F_decision_killed = prdmswtn(t, mu_drift, b, A, s, t0, sv,
                                              0.0, rate, n_gauss_nodes,
                                              false, 1, false, posdrift, 1.0);
    return local_combo_response_cdf_exp(t, F_decision, F_decision_killed,
                                        lambda_g, lambda_k, log_out);
  }

  if (kill_shape == 2 && sv <= 1e-10) {
    const double upper = emc2_isfinite(t) ? t - t0 : R_PosInf;
    tilted_wald_moments_0_1_2 m =
      spv_wald_tilted_moments_0_1_2(
        upper, b, A, mu_drift, s, lambda_g + lambda_k, posdrift);
    if (m.ok) {
      return local_combo_cdf_erlang2_from_moments(
        t, t0, lambda_g, lambda_k, m, log_out);
    }
  }

  // Drift-variability cases still need the drift mixture of the second tilted
  // moment; keep raw-time integration there.
  const double rate_scale = lambda_g + lambda_k;
  const auto density_fn = [&](double u) {
    return drdmswtn_local_combo(u, mu_drift, b, A, s, t0, sv,
                                lambda_g, lambda_k, n_gauss_nodes, false, kill_shape, posdrift);
  };
  double out = emc2_isfinite(t)
    ? integrate_density_adaptive_finite(t, density_fn)
    : integrate_density_adaptive_infinite(rate_scale, density_fn);
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline double dgbm_local_combo(double t, double mu, double b, double A,
                               double sigma, double t0,
                               double lambda_g, double lambda_k,
                               bool log_out = false, int kill_shape = 1,
                               double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = dgbm_local_combo(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, 1.0);
    const double e2 = dgbm_local_combo(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // EAM computed with raw t and t0; combined on the log scale.
  const double log_fD = dgbm(t, mu, b, A, sigma, t0, 0.0, 0.0, true, 1, false, 1.0);
  const double log_FD = pgbm(t, mu, b, A, sigma, t0, 0.0, 0.0, true, 1, false, 1.0);
  return local_combo_response_pdf_log(t, log_fD, log_FD,
                                      lambda_g, lambda_k, kill_shape, log_out);
}

inline double pgbm_local_combo(double t, double mu, double b, double A,
                               double sigma, double t0,
                               double lambda_g, double lambda_k,
                               bool log_out = false, int kill_shape = 1,
                               double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = pgbm_local_combo(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, 1.0);
    const double e2 = pgbm_local_combo(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;

  if (kill_shape <= 1) {
    const double rate = lambda_g + lambda_k;
    const double F_decision = pgbm(t, mu, b, A, sigma, t0,
                                   0.0, 0.0, false, 1, false);
    const double F_decision_killed = pgbm(t, mu, b, A, sigma, t0,
                                          0.0, rate, false, 1, false);
    return local_combo_response_cdf_exp(t, F_decision, F_decision_killed,
                                        lambda_g, lambda_k, log_out);
  }

  if (kill_shape == 2) {
    const double upper = emc2_isfinite(t) ? t - t0 : R_PosInf;
    tilted_wald_moments_0_1_2 m =
      gbm_tilted_moments_0_1_2(upper, mu, b, A, sigma, lambda_g + lambda_k);
    if (m.ok) {
      return local_combo_cdf_erlang2_from_moments(
        t, t0, lambda_g, lambda_k, m, log_out);
    }
  }

  // Fallback for numerical edge cases.
  const double rate_scale = lambda_g + lambda_k;
  const auto density_fn = [&](double u) {
    return dgbm_local_combo(u, mu, b, A, sigma, t0, lambda_g, lambda_k, false, kill_shape);
  };
  double out = emc2_isfinite(t)
    ? integrate_density_adaptive_finite(t, density_fn)
    : integrate_density_adaptive_infinite(rate_scale, density_fn);
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline double dswtn_core(double t_adj, double mu_drift, double threshold,
                         double s, double t0, double sv, double lambda,
                         double log_norm, bool log_out, int kill_shape, bool guess, bool posdrift) {
  // t_adj is EAM time (rt - t0). t_raw = t_adj + t0 is raw rt for erlang.
  const double t_raw = t_adj + t0;

  // When EAM hasn't started, only erlang guess contributes (S_R = 1).
  if (t_adj <= 0.0) {
    if (!guess || lambda <= 0.0) return log_out ? R_NegInf : 0.0;
    const double log_fG = erlang_log_pdf(t_raw, lambda, kill_shape);
    return log_out ? log_fG : std::exp(log_fG);
  }

  const double v  = sv * sv;
  const double s2 = s * s;
  const double tv = t_adj * v;
  const double den_common     = tv + s2;
  const double log_den_common = std::log(den_common);

  const double term_log_threshold = std::log(threshold);
  const double term_log_denom     = -0.5 * (std::log(M_PI) + M_LN2 +
                                             3.0 * std::log(t_adj) + log_den_common);
  const double term_log_exp       = -(std::pow(threshold - mu_drift * t_adj, 2.0)) /
                                      (2.0 * t_adj * den_common);

  const double mu_new    = (threshold * v + mu_drift * s2) / den_common;
  const double sigma_new = std::sqrt(s2 * v / den_common);

  // EAM hit density multiplied by Erlang kill survival.
  const double log_f_hit = term_log_threshold + term_log_denom -
                         log_norm + term_log_exp +
                         erlang_log_surv(t_raw, lambda, kill_shape);

  if (!guess || lambda <= 0.0) {
    if (ISNAN(log_f_hit)) return log_out ? R_NegInf : 0.0;
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // Guess density: f_K(t_raw) * S_R(t_adj); EAM survivor uses EAM time
  const double log_sk = erlang_log_pdf(t_raw, lambda, kill_shape);
  const double log_sr = log_surv_from_log_cdf(
      pswtn(t_raw, mu_drift, threshold, s, t0, sv, 0.0, 0.0, true, 1, false, posdrift, 1.0));
  const double log_f_guess = log_sk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  if (ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double pswtn_killed_quad(double t_adj, double mu_drift, double threshold,
                                double s, double t0, double sv, double lambda,
                                int n_gauss_nodes, int kill_shape, bool guess, bool posdrift = true) {
  // t_adj is EAM time; t_adj + t0 is raw rt. Pass t0 to pwald so erlang uses raw time.
  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& gl_nodes   = gl.x;
  const std::vector<double>& gl_weights = gl.w;

  double integral = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);  // map [-1,1] -> [0,1]
    const double p = std::fmin(std::nextafter(1.0, 0.0), std::fmax(1e-15, u));
    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    // Integrate defective fixed-drift kernels across the full drift
    // distribution, then normalise outside if posdrift=true.
    integral += gl_weights[j] * pwald(t_adj + t0, drift_j, threshold, 0.0, s, t0, lambda, lambda, false, kill_shape, guess, false);
  }
  double out_val = 0.5 * integral;

  if (posdrift) {
    out_val *= std::exp(-log_swtn_hit_mass(mu_drift, sv, s, threshold));
  }
  return out_val;
}

inline double pswtn_killed_inf_quad(double threshold, double mu_drift,
                                    double sv, double s, double lambda,
                                    int n_gauss_nodes = 20,
                                    bool log_out = false, int kill_shape = 1,
                                    bool posdrift = true) {
  auto finish_log = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (s <= 1e-10 || sv < 0.0 || lambda < 0.0) return NA_REAL;

  const double s2 = s * s;

  if (sv <= 1e-10) {
    const double nu = std::sqrt(mu_drift * mu_drift + 2.0 * s2 * lambda);
    double log_hit = threshold * (mu_drift - nu) / s2;
    if (kill_shape >= 2 && nu > FPM_EPSILON)
      log_hit += std::log1p(lambda * threshold / nu);
    if (posdrift)
      log_hit -= log_swtn_hit_mass(mu_drift, sv, s, threshold);
    return finish_log(log_hit);
  }

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& gl_nodes   = gl.x;
  const std::vector<double>& gl_weights = gl.w;

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);

    const double p = std::fmin(std::nextafter(1.0, 0.0),
                               std::fmax(1e-15, u));

    const double drift_j =
      mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double nu_j = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
    double log_hit_j = threshold * (drift_j - nu_j) / s2;
    if (kill_shape >= 2 && nu_j > FPM_EPSILON)
      log_hit_j += std::log1p(lambda * threshold / nu_j);

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  const double log_mass = log_integral - M_LN2;
  if (!posdrift) return finish_log(log_mass);
  return finish_log(log_mass - log_swtn_hit_mass(mu_drift, sv, s, threshold));
}

inline double prdmswtn_killed_inf_quad(double b, double mu_drift, double A,
                                       double sv, double s, double t0, double lambda,
                                       int n_gauss_nodes = 20,
                                       bool log_out = false, int kill_shape = 1, bool posdrift = true) {
  auto finish_log = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (b <= 1e-10) return log_out ? 0.0 : 1.0;
  if (s <= 1e-10 || sv < 0.0 || A < 0.0 || lambda < 0.0) return NA_REAL;

  const double s2 = s * s;

  if (sv <= 1e-10) {
    double log_p;
    if (kill_shape >= 2) {
      // Use the defective fixed-drift kernel here; the SPV posdrift
      // normalisation is applied explicitly below when requested. pwald
      // applies the Erlang-2 t0 kill head-start internally (kill clock runs
      // from stimulus onset, not from evidence onset t0).
      log_p = pwald(R_PosInf, mu_drift, b, A, s, t0, lambda, lambda, true, kill_shape, false, false);
    } else {
      const double nu = std::sqrt(mu_drift * mu_drift + 2.0 * s2 * lambda);
      const double eta = (mu_drift - nu) / s2;
      // Closed-form defective eventual hit mass for the killed Wald (kill from
      // evidence onset), shifted by exp(-lambda*t0) for the exponential kill's
      // head-start over [0, t0] before evidence accumulation begins.
      log_p = log_wald_posdrift_hit_normalizer(true, 0.5 * eta * s2, s, b, 0.0, A) - lambda * t0;
    }
    if (posdrift) {
      log_p -= log_wald_posdrift_hit_normalizer(true, mu_drift, s, b, 0.0, A);
    }
    return finish_log(log_p);
  }

  const int n_nodes = std::max(1, n_gauss_nodes);
  const GLRule& gl = gl_get_rule(n_nodes);
  const std::vector<double>& gl_nodes   = gl.x;
  const std::vector<double>& gl_weights = gl.w;

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);
    
    const double p = std::fmin(std::nextafter(1.0, 0.0),
                               std::fmax(1e-15, u));
    
    const double drift_j =
      mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    double log_hit_j;
    if (kill_shape >= 2) {
      // Defective fixed-drift kernel inside the drift mixture; pwald applies
      // the Erlang-2 t0 kill head-start (kill runs from stimulus onset). The
      // SPV posdrift normaliser is applied once after integration if requested.
      log_hit_j = pwald(R_PosInf, drift_j, b, A, s, t0, lambda, lambda, true, kill_shape, false, false);
    } else {
      const double nu_j  = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
      const double eta_j = (drift_j - nu_j) / s2;
      // Closed-form defective eventual hit mass for the fixed-drift kernel,
      // shifted by exp(-lambda*t0) for the exponential kill's head-start.
      log_hit_j = log_wald_posdrift_hit_normalizer(true, 0.5 * eta_j * s2, s, b, 0.0, A) - lambda * t0;
    }

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  const double log_mass = log_integral - M_LN2;
  if (!posdrift) return finish_log(log_mass);
  return finish_log(log_mass - log_swtn_spv_hit_mass_full(mu_drift, sv, s, b, A, n_gauss_nodes));
}



// --------------------------------------------------------------------------
// RDMSWTN under the finite linear exhaustion clock
//
// x = t - t0, q(x) = x - x^2/(2*tau), q'(x) = 1 - x/tau,
// 0 < x < tau.  Once x reaches tau the operational-time budget is Q=tau/2,
// so the CDF freezes and the density is identically zero.
// --------------------------------------------------------------------------

inline double rdmswtn_tt_q(double x, double tau) {
  // tau = +Inf is the exact identity-clock limit.
  if (!R_FINITE(tau)) return x;
  return x * (1.0 - 0.5 * x / tau);
}


double drdmswtn_tt(double t, double mu_drift, double b, double A,
                   double s RDM_DEFAULT_ARGUMENT(1.0), double t0 RDM_DEFAULT_ARGUMENT(0.0), double sv RDM_DEFAULT_ARGUMENT(0.0), double tau RDM_DEFAULT_ARGUMENT(1.0),
                   bool log_out RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true));

double prdmswtn_tt(double t, double mu_drift, double b, double A,
                   double s RDM_DEFAULT_ARGUMENT(1.0), double t0 RDM_DEFAULT_ARGUMENT(0.0), double sv RDM_DEFAULT_ARGUMENT(0.0), double tau RDM_DEFAULT_ARGUMENT(1.0),
                   bool log_out RDM_DEFAULT_ARGUMENT(false), bool posdrift RDM_DEFAULT_ARGUMENT(true));



#undef RDM_DEFAULT_ARGUMENT







#endif
