#ifndef rdm_h
#define rdm_h

#define _USE_MATH_DEFINES
#include <cmath>
#include <RcppArmadillo.h>
#include "wald_functions.h"
#include "composite_functions.h"
#include "gaussian.h"

using namespace Rcpp;

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
  if (std::abs(q) <= FPM_EPSILON) {
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
  if (std::abs(q) <= FPM_EPSILON) {
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

// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false){
  int n = t.size();
  NumericVector pdf(n);
  for (int i = 0; i < n; i++){
    double ti = t[i] - t0[i];
    if (ti <= 0){
      pdf[i] = log_out? R_NegInf : 0.0;
    } else {
      pdf[i] = digt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
    }
  }
  return pdf;
}


// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false){
  int n = t.size();
  NumericVector cdf(n);
  for (int i = 0; i < n; i++){
    double ti = t[i] - t0[i];
    if (ti <= 0){
      cdf[i] = log_out? 0.0 : R_NegInf;
    } else {
      cdf[i] = pigt_impl(ti, B[i] + .5 * A[i], v[i], .5 * A[i], log_out);
    }
  }
  return cdf;
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
// All functions expect t_adj = t - t0 (non-decision time already removed).
// The `s` (diffusion) parameter is absorbed by the caller: pass s=1 once
// parameters have been pre-scaled by s.
// ==========================================================================

double pwald(double t, double b, double mu, double sigma, double A,
             double k, bool log_out, int kill_shape, bool guess);

// [[Rcpp::export]]
double dwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0, double k = 0.0,
             bool log_out = false, int kill_shape = 1, bool guess = false) {
  double x_lo = 0.0;
  double x_hi = x_lo + A;

  const double var = sigma * sigma * t;
  if (var <= FPM_EPSILON || t <= FPM_EPSILON || sigma <= 0.0 || k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }

  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  double log_f_hit = R_NegInf;

  // 1. Evidence-based hit density (f_R * S_K)
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;
    if (d <= 0.0) {
      if (t == 0.0) {
        return log_out ? R_PosInf : std::numeric_limits<double>::infinity();
      }
      return log_out ? R_NegInf : 0.0;
    }
    const double delta = d - mu * t;
    log_f_hit = std::log(d) - std::log(t) + Gstar(var, delta, true) + erlang_log_surv(t, k, kill_shape);
  } else {
    // SPV form
    const double mu_new = b - mu * t;
    const double log_pdf_hi = Gstar(var, x_hi - mu_new, true);
    const double log_pdf_lo = Gstar(var, x_lo - mu_new, true);
    const double log_cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi, true);

    signed_log term1 = make_signed_log(R_NegInf, 0);
    if (std::abs(mu) > FPM_EPSILON && log_cdf_integral != R_NegInf) {
      term1 = make_signed_log(std::log(std::abs(mu)) + std::log(t) + log_cdf_integral, mu > 0.0 ? 1 : -1);
    }
    signed_log pdf_hi = make_signed_log(log_pdf_hi, 1);
    signed_log pdf_lo = make_signed_log(log_pdf_lo, 1);
    signed_log pdf_diff = signed_log_sub(pdf_hi, pdf_lo);
    signed_log term2 = make_signed_log(R_NegInf, 0);
    if (pdf_diff.sign != 0) {
      term2 = make_signed_log(std::log(var) + pdf_diff.log_abs, pdf_diff.sign);
    }
    signed_log total = signed_log_add(term1, term2);

    if (total.sign > 0 && total.log_abs != R_NegInf) {
      log_f_hit = total.log_abs - std::log(t) - std::log(span) + erlang_log_surv(t, k, kill_shape);
    }
  }

  if (!guess || k <= 0.0) {
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // 2. Guess-based density (f_K * S_R)
  // S_R = 1 - pwald(t, ..., k=0)
  const double log_sk = erlang_log_pdf(t, k, kill_shape);
  const double log_sr = std::log1p(-std::max(0.0, std::min(1.0, pwald(t, b, mu, sigma, A, 0.0, false, 1, false))));
  const double log_f_guess = log_sk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

// --------------------------------------------------------------------------
// Wald / killed Wald CDF with uniform start-point variability.
// Returns sub-CDF:
//   P(T_wald <= t AND T_kill > T_wald)
// k = 0 gives standard Wald, nested.
// t = Inf returns total finite-response probability.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double pwald(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
             double k = 0.0, bool log_out = false, int kill_shape = 1, bool guess = false) {
  auto finish = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;

    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;

    return log_out ? log_p : std::exp(log_p);
  };

  if (sigma <= 0.0 || k < 0.0) {
    return NA_REAL;
  }

  if (t <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  if (guess && k > 0.0) {
    // Mixture CDF: F_i(t) = 1 - S_R,i(t) * S_K,i(t)
    const double fw = std::max(0.0, std::min(1.0, pwald(t, b, mu, sigma, A, 0.0, false, 1, false)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t, k, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  double x_lo = 0.0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);

  const double span = x_hi - x_lo;
  const double sig2 = sigma * sigma;

  // Fast path: no kill — use canonical Wald SPV for finite t.
  if (k <= 0.0) {
    if (emc2_isfinite(t)) {
      const double inv_s = 1.0 / sigma;
      const double cdf = pwald_k0(t, b * inv_s, mu * inv_s, A * inv_s);
      const double cl = std::max(0.0, std::min(1.0, cdf));
      if (log_out) return (cl <= 0.0) ? R_NegInf : std::log(cl);
      return cl;
    }
    // t = Inf with k=0: fall through — log_p_inf() handles it via nu = |mu|, eta1 correct.
  }

  const double nu = std::sqrt(mu * mu + 2.0 * sig2 * k);
  const double eta1 = (mu - nu) / sig2;

  auto log_p_inf = [&]() -> double {
    if (span <= FPM_EPSILON) {
      const double d = b - x_hi;
      if (d <= 0.0) return 0.0;
      const double log_p1 = eta1 * d;
      if (kill_shape <= 1 || nu <= FPM_EPSILON) return log_p1;
      // Erlang-2 eventual hit: exp(eta1*d) * (1 + lambda*d/nu)
      return log_p1 + std::log1p(k * d / nu);
    }

    if (kill_shape >= 2 && k > FPM_EPSILON) {
      // Analytic: (1/span) * int_{d_lo}^{d_hi} exp(eta1*d)*(1 + k*d/nu) dd
      // where d = b - x, d_lo = b - x_hi, d_hi = b - x_lo.
      if (b <= x_hi) return NA_REAL;
      const double d_lo = b - x_hi;
      const double d_hi = b - x_lo;
      double I0, I1;
      if (std::abs(eta1) <= FPM_EPSILON) {
        I0 = d_hi - d_lo;
        I1 = 0.5 * (d_hi * d_hi - d_lo * d_lo);
      } else {
        const double e_hi = std::exp(eta1 * d_hi);
        const double e_lo = std::exp(eta1 * d_lo);
        I0 = (e_hi - e_lo) / eta1;
        I1 = e_hi * (d_hi / eta1 - 1.0 / (eta1 * eta1)) -
             e_lo * (d_lo / eta1 - 1.0 / (eta1 * eta1));
      }
      const double p = (I0 + (k / nu) * I1) / span;
      if (p <= 0.0) return R_NegInf;
      if (p >= 1.0) return 0.0;
      return std::log(p);
    }

    if (std::abs(eta1) < FPM_EPSILON) {
      return 0.0;
    }

    // Erlang-1 SPV: (1/span) * int_0^A exp(eta1*(b-x)) dx
    // = exp(eta1*b)/span * int_0^A exp(q*x) dx,  q = -eta1 > 0
    // = exp(eta1*b) * (exp(q*A)-1) / (q*span)
    const double q = -eta1;
    if (q <= FPM_EPSILON) return 0.0;
    return eta1 * b +
           log_diff_exp(q * x_hi, q * x_lo) -  // log(exp(q*A)-1), x_hi=A, x_lo=0
           std::log(q) -
           std::log(span);
  };

  if (!emc2_isfinite(t)) {
    return finish(log_p_inf());
  }

  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  // Point-start case.
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;

    if (d <= 0.0) {
      return log_out ? 0.0 : 1.0;
    }

    const double log_prefactor = d * (mu - nu) / sig2;

    const double term1_arg = (nu * t - d) / st;
    const double term2_arg = (-nu * t - d) / st;

    const double log_cdf1 = pnorm_std(term1_arg, true, true);
    const double log_cdf2 = pnorm_std(term2_arg, true, true);

    const double log_exp_term = 2.0 * nu * d / sig2;

    const double log_cdf_nu =
      log_sum_exp(log_cdf1, log_exp_term + log_cdf2);

    if (kill_shape >= 2 && k > FPM_EPSILON && nu > FPM_EPSILON) {
      // Erlang-2: add lambda * M_W term.
      // M_W (unnorm) = Phi(z1) - exp_term*Phi(z2), which equals C^{-1} * dL/d(-lambda).
      const double log_b2 = log_exp_term + log_cdf2;
      const double log_mw = (log_cdf1 > log_b2)
                              ? log_diff_exp(log_cdf1, log_b2)
                              : R_NegInf;
      const double log_w = std::log(k) + std::log(d) - std::log(nu);  // log(lambda*d/nu)
      const double log_cdf_e2 = (log_mw > R_NegInf)
                                  ? log_sum_exp(log_cdf_nu, log_w + log_mw)
                                  : log_cdf_nu;
      return finish(log_prefactor + log_cdf_e2);
    }

    return finish(log_prefactor + log_cdf_nu);
  }

  if (kill_shape >= 2 && k > FPM_EPSILON) {
    // SPV + Erlang-2: analytic via slog_int_{eta,d_eta}_pnorm.
    // total = (T1 + T2 + (k/nu)*(D1-D2)) / span
    if (b <= x_hi) return NA_REAL;
    const double a    = 1.0 / st;
    const double c1   = (nu * t - b) / st;
    const double c2   = (-nu * t - b) / st;
    const double eta2 = (mu + nu) / sig2;
    signed_log T1 = slog_int_eta_pnorm(eta1, b, a, c1, x_lo, x_hi);
    signed_log T2 = slog_int_eta_pnorm(eta2, b, a, c2, x_lo, x_hi);
    signed_log D1 = slog_int_d_eta_pnorm(eta1, b, a, c1, x_lo, x_hi);
    signed_log D2 = slog_int_d_eta_pnorm(eta2, b, a, c2, x_lo, x_hi);
    signed_log total  = signed_log_add(T1, T2);
    signed_log Ddiff  = signed_log_sub(D1, D2);
    if (Ddiff.sign != 0 && Ddiff.log_abs != R_NegInf) {
      signed_log extra = make_signed_log(
        std::log(k) - std::log(nu) + Ddiff.log_abs, Ddiff.sign);
      total = signed_log_add(total, extra);
    }
    if (total.sign <= 0 || total.log_abs == R_NegInf)
      return log_out ? R_NegInf : 0.0;
    return finish(total.log_abs - std::log(span));
  }

  // SPV + Erlang-1: analytic.
  const double a = 1.0 / st;

  // First term:
  // exp[(mu - nu)(b-x)/sig2] *
  // Phi((x + nu*t - b) / st)
  const double c1 = (nu * t - b) / st;

  const double log_term1 =
    eta1 * b +
    log_integrate_exp_times_normal_cdf(-eta1, a, c1, x_lo, x_hi);

  // Second term:
  // exp[(mu + nu)(b-x)/sig2] *
  // Phi((x - nu*t - b) / st)
  const double eta2 = (mu + nu) / sig2;
  const double c2 = (-nu * t - b) / st;

  const double log_term2 =
    eta2 * b +
    log_integrate_exp_times_normal_cdf(-eta2, a, c2, x_lo, x_hi);

  const double log_cdf_val =
    log_sum_exp(log_term1, log_term2) - std::log(span);

  return finish(log_cdf_val);
}

// --------------------------------------------------------------------------
// GBM first-passage with physical-space start-point variability.
// Process: dX_t = mu X_t dt + sigma X_t dW_t, start X0 ~ Unif(1, 1 + A).
// Boundary: absorb at X_t = b (b > 0).
// In log-space this is a Wald kernel with drift (mu - 0.5*sigma^2), diffusion
// sigma, and start range [log(1), log(1 + A)] weighted by exp(x0).
// --------------------------------------------------------------------------
double pgbm(double t, double b, double mu, double sigma, double A,
            double k, bool log_out, int kill_shape, bool guess);

// [[Rcpp::export]]
double dgbm(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
            double k = 0.0, bool log_out = false, int kill_shape = 1, bool guess = false) {
  if (!(t > FPM_EPSILON) || !(sigma > FPM_EPSILON) || !(b > FPM_EPSILON) || !(A >= 0.0) || k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  if (A > b) return log_out ? R_NegInf : 0.0;

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double log_mu = mu - 0.5 * sigma * sigma;
  const double var = sigma * sigma * t;
  const double norm_const = A;

  double log_f_hit = R_NegInf;

  // 1. Evidence hit density (f_R * S_K)
  if (!(norm_const > FPM_EPSILON)) {
    const double d = log_b - x_hi;  // = log_b when A = 0
    if (d <= 0.0) {
      log_f_hit = (t == 0.0) ? std::numeric_limits<double>::infinity() : R_NegInf;
    } else {
      const double delta = d - log_mu * t;
      log_f_hit = std::log(d) - std::log(t) + Gstar(var, delta, true)
                  + erlang_log_surv(t, k, kill_shape);
    }
  } else {
    // SPV: normalization for start density on [1, 1 + A] mapped to log-space.
    const double mu_new_exp = log_b - log_mu * t + var;
    const double exp_factor_log = log_b - log_mu * t + 0.5 * var;
    const double pdf_hi = Gstar(var, x_hi - mu_new_exp);
    const double pdf_lo = Gstar(var, x_lo - mu_new_exp);
    const double cdf_integral = Gstar_Integral(var, mu_new_exp, x_lo, x_hi);

    const double integral_term1 = (log_mu * t - var) * cdf_integral;
    const double integral_term2 = var * (pdf_hi - pdf_lo);
    const double integral_result = integral_term1 + integral_term2;
    const double pdf_val = (std::exp(exp_factor_log) * integral_result) / (norm_const * t);
    if (pdf_val > 0.0 && std::isfinite(pdf_val)) {
       log_f_hit = std::log(pdf_val) + erlang_log_surv(t, k, kill_shape);
    }
  }

  if (!guess || k <= 0.0) {
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // 2. Guess density (f_K * S_R)
  const double log_sk = erlang_log_pdf(t, k, kill_shape);
  const double log_sr = std::log1p(-std::max(0.0, std::min(1.0, pgbm(t, b, mu, sigma, A, 0.0, false, 1, false))));
  const double log_f_guess = log_sk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}

// [[Rcpp::export]]
double pgbm(double t, double b, double mu, double sigma = 1.0, double A = 0.0,
            double k = 0.0, bool log_out = false, int kill_shape = 1, bool guess = false) {
  auto finish = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (!(b > FPM_EPSILON) || !(sigma > FPM_EPSILON) || !(A >= 0.0) || k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  if (A > b) return log_out ? R_NegInf : 0.0;

  if (guess && k > 0.0) {
    // Mixture CDF: F_i(t) = 1 - S_R,i(t) * S_K,i(t)
    const double fw = std::max(0.0, std::min(1.0, pgbm(t, b, mu, sigma, A, 0.0, false, 1, false)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t, k, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double alpha = mu - 0.5 * sigma * sigma;  // log-space drift
  const double sig2 = sigma * sigma;
  const double norm_const = A;

  // Tilted log-space drift for kill (= |alpha| when k = 0).
  const double ak     = std::sqrt(alpha * alpha + 2.0 * sig2 * k);
  const double c_eta  = (alpha - ak) / sig2;   // eta1 in log-space
  const double q      = 1.0 - c_eta;           // Jacobian-corrected coefficient (1 - eta1)
  const double k_log  = 2.0 * ak / sig2;       // 2*ak/sigma^2

  // ---- t = Inf: total finite-response probability ----
  auto log_p_inf = [&]() -> double {
    if (!(norm_const > FPM_EPSILON)) {
      // Point-start: P(hit before kill) = b^c_eta, Erlang-2 adds log1p correction.
      double log_hit = c_eta * log_b;
      if (kill_shape >= 2 && ak > FPM_EPSILON)
        log_hit += std::log1p(k * log_b / ak);
      return log_hit;
    }

    if (kill_shape >= 2) {
      // P(hit before Erlang-2 kill) =
      // exp(c_eta*log_b)/A * [I0 + (k/ak)*I1],
      // I0 = int exp(qy)dy, I1 = int (log_b-y)exp(qy)dy, y in [0, x_hi].
      const double log_prefactor = c_eta * log_b - std::log(norm_const);
      double I0, I1;
      if (std::abs(q) <= FPM_EPSILON) {
        I0 = x_hi;
        I1 = log_b * x_hi - 0.5 * x_hi * x_hi;
      } else {
        const double e_hi = std::exp(q * x_hi);
        I0 = (e_hi - 1.0) / q;
        I1 = (e_hi * (log_b / q - x_hi / q + 1.0 / (q * q)) -
              (log_b / q + 1.0 / (q * q)));
      }
      const double p = std::exp(log_prefactor) * (I0 + (k / ak) * I1);
      if (p <= 0.0) return R_NegInf;
      if (p >= 1.0) return 0.0;
      return std::log(p);
    }

    // Erlang-1 SPV analytic:
    // exp(c_eta * log_b) / A * integral_0^{log1p(A)} exp(q * y) dy
    const double log_prefactor = c_eta * log_b - std::log(norm_const);
    if (std::abs(q) < FPM_EPSILON)
      return log_prefactor + std::log(x_hi);   // q = 0: integral = x_hi
    if (q > 0.0)
      return log_prefactor + log_diff_exp(q * x_hi, 0.0) - std::log(q);
    return log_prefactor + log_diff_exp(0.0, q * x_hi) - std::log(-q);
  };

  if (!emc2_isfinite(t)) {
    return finish(log_p_inf());
  }

  if (t <= FPM_EPSILON) return log_out ? R_NegInf : 0.0;

  const double st = sigma * std::sqrt(t);
  if (st <= FPM_EPSILON) return log_out ? R_NegInf : 0.0;

  const double a  = 1.0 / st;
  const double c1 = (ak * t - log_b) / st;   // c for Phi in term1
  const double c2 = (-ak * t - log_b) / st;  // c for Phi in term2

  // ---- Point-start: A = 0 ----
  if (!(norm_const > FPM_EPSILON)) {
    const double d = log_b - x_hi;  // = log_b when A = 0

    if (d <= 0.0) return log_out ? 0.0 : 1.0;

    const double log_prefactor = c_eta * d;

    const double log_cdf1    = pnorm_std((ak * t - d) / st, true, true);
    const double log_cdf2    = pnorm_std((-ak * t - d) / st, true, true);
    const double log_exp_t   = k_log * d;
    const double log_cdf_ak  = log_sum_exp(log_cdf1, log_exp_t + log_cdf2);

    if (kill_shape >= 2 && k > FPM_EPSILON && ak > FPM_EPSILON) {
      const double log_b2  = log_exp_t + log_cdf2;
      const double log_mw  = (log_cdf1 > log_b2) ? log_diff_exp(log_cdf1, log_b2) : R_NegInf;
      const double log_w   = std::log(k) + std::log(d) - std::log(ak);
      const double log_e2  = (log_mw > R_NegInf) ? log_sum_exp(log_cdf_ak, log_w + log_mw)
                                                  : log_cdf_ak;
      return finish(log_prefactor + log_e2);
    }

    return finish(log_prefactor + log_cdf_ak);
  }

  // ---- SPV + Erlang-2 ----
  if (kill_shape >= 2) {
    // Analytic for the common case where all starts are below boundary (log_b > x_hi):
    // F_K(t) = exp(c_eta*log_b)/A * [T1 + T2 + (k/ak)(D1-D2)]
    // with
    // T1 = int exp(qy) Phi(a*y + c1) dy
    // T2 = exp(k_log*log_b) int exp((q-k_log)y) Phi(a*y + c2) dy
    // D1 = int (log_b-y) exp(qy) Phi(a*y + c1) dy
    // D2 = exp(k_log*log_b) int (log_b-y) exp((q-k_log)y) Phi(a*y + c2) dy
    signed_log J1 = slog_int_exp_pnorm(q, a, c1, x_lo, x_hi);
    signed_log J2 = slog_int_exp_pnorm(q - k_log, a, c2, x_lo, x_hi);
    signed_log K1 = slog_int_x_exp_pnorm(q, a, c1, x_lo, x_hi);
    signed_log K2 = slog_int_x_exp_pnorm(q - k_log, a, c2, x_lo, x_hi);

    signed_log T1 = J1;
    signed_log T2 = J2;
    if (T2.sign != 0 && T2.log_abs != R_NegInf)
      T2 = make_signed_log(T2.log_abs + k_log * log_b, T2.sign);

    signed_log bJ1 = make_signed_log(R_NegInf, 0);
    if (J1.sign != 0 && J1.log_abs != R_NegInf && log_b > 0.0)
      bJ1 = make_signed_log(std::log(log_b) + J1.log_abs, J1.sign);
    signed_log D1 = signed_log_sub(bJ1, K1);

    signed_log bJ2 = make_signed_log(R_NegInf, 0);
    if (J2.sign != 0 && J2.log_abs != R_NegInf && log_b > 0.0)
      bJ2 = make_signed_log(std::log(log_b) + J2.log_abs, J2.sign);
    signed_log D2 = signed_log_sub(bJ2, K2);
    if (D2.sign != 0 && D2.log_abs != R_NegInf)
      D2 = make_signed_log(D2.log_abs + k_log * log_b, D2.sign);

    signed_log total = signed_log_add(T1, T2);
    signed_log Ddiff = signed_log_sub(D1, D2);
    if (Ddiff.sign != 0 && Ddiff.log_abs != R_NegInf) {
      signed_log extra = make_signed_log(
        std::log(k) - std::log(ak) + Ddiff.log_abs, Ddiff.sign);
      total = signed_log_add(total, extra);
    }
    if (total.sign <= 0 || total.log_abs == R_NegInf)
      return log_out ? R_NegInf : 0.0;
    return finish(c_eta * log_b + total.log_abs - std::log(norm_const));
  }

  // ---- SPV + Erlang-1: analytic (log-stable) ----
  // F_K(t) = (exp(c_eta * log_b) / A) *
  //   [ ∫_0^{x_hi} exp(q * y) Phi(a*y + c1) dy
  //   + exp(k_log*log_b) * ∫_0^{x_hi} exp((q-k_log)*y) Phi(a*y + c2) dy ]
  const double log_term1 = c_eta * log_b +
    log_integrate_exp_times_normal_cdf(q, a, c1, x_lo, x_hi);
  const double log_term2 = (c_eta + k_log) * log_b +
    log_integrate_exp_times_normal_cdf(q - k_log, a, c2, x_lo, x_hi);

  const double log_cdf_val = log_sum_exp(log_term1, log_term2) - std::log(norm_const);

  return finish(log_cdf_val);
}

// --------------------------------------------------------------------------
// SWTN: Shifted Wald with Truncated-Normal drift variability.
// No start-point variability (fixed threshold = threshold param).
// sv  = SD of between-trial drift distribution N(mu_drift, sv^2) [truncated at 0].
// s   = within-trial diffusion coefficient.
// c   = lower truncation for drift (default 0 for positive-drift assumption).
// --------------------------------------------------------------------------
inline double dswtn_core(double t_adj, double threshold, double mu_drift,
                         double sv, double s, double lambda, double c,
                         double log_prob_gt_c, bool log_out, int kill_shape = 1, bool guess = false);

double pswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s, double lambda, double c,
             bool log_out, int kill_shape, bool guess);

inline Rcpp::List get_gl_nodes_weights(int n_gauss_nodes);

// [[Rcpp::export]]
double dswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double lambda=0.0, double c = 0.0,
             bool log_out = false, int kill_shape = 1, bool guess = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (sv < 0.0) return R_NaN;
  if (s <= 1e-10) return log_out ? R_NegInf : 0.0;

  // sv == 0: reduces to standard Wald
  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, threshold, mu_drift, s, 0.0, lambda, log_out, kill_shape, guess);
  }

  // Normalisation P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

  return dswtn_core(t_adj, threshold, mu_drift, sv, s, lambda, c, log_prob_gt_c, log_out, kill_shape, guess);
}

// --------------------------------------------------------------------------
// SWTN CDF.
// --------------------------------------------------------------------------
inline Rcpp::List get_gl_nodes_weights(int n_gauss_nodes) {
  const int n_nodes = std::max(1, n_gauss_nodes);
  return (n_nodes == 20) ? get_gl20()
                         : Rcpp::as<Rcpp::List>(gauss_quad(n_nodes, "legendre"));
}

double drdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s, double lambda, double c,
                int n_gauss_nodes, bool log_out, int kill_shape, bool guess);
double prdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s, double lambda, double c,
                int n_gauss_nodes, bool log_out, int kill_shape, bool guess);

template <typename DensityFn>
inline double integrate_density_gl20_finite(double t_upper, DensityFn&& density_fn) {
  if (!(t_upper > 0.0)) return 0.0;
  const Rcpp::List& gl = get_gl20();
  const Rcpp::NumericVector nodes = gl["nodes"];
  const Rcpp::NumericVector weights = gl["weights"];
  double acc = 0.0;
  for (int j = 0; j < nodes.size(); ++j) {
    const double u = 0.5 * t_upper * (nodes[j] + 1.0);
    acc += weights[j] * density_fn(u);
  }
  return 0.5 * t_upper * acc;
}

template <typename DensityFn>
inline double integrate_density_gl20_infinite(double rate_scale, DensityFn&& density_fn) {
  const double scale = std::fmax(rate_scale, 1e-8);
  const Rcpp::List& gl = get_gl20();
  const Rcpp::NumericVector nodes = gl["nodes"];
  const Rcpp::NumericVector weights = gl["weights"];
  double acc = 0.0;
  for (int j = 0; j < nodes.size(); ++j) {
    const double q = 0.5 * (nodes[j] + 1.0);
    const double qq = std::fmin(1.0 - 1e-12, std::fmax(1e-15, q));
    const double t = -std::log1p(-qq) / scale;
    const double jac = 1.0 / (scale * (1.0 - qq));
    acc += weights[j] * density_fn(t) * jac;
  }
  return 0.5 * acc;
}

inline double local_combo_response_pdf(double t, double f_decision, double F_decision,
                                       double lambda_g, double lambda_k,
                                       int kill_shape, bool log_out) {
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;
  const bool use_guess = lambda_g > 0.0;
  const bool use_kill  = lambda_k > 0.0;
  if (!use_guess && !use_kill) {
    if (!(f_decision > 0.0) || !std::isfinite(f_decision)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(f_decision) : f_decision;
  }

  const double log_fD = (f_decision > 0.0 && std::isfinite(f_decision)) ? std::log(f_decision) : R_NegInf;
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

inline double drdmswtn_local_combo(double t_adj, double b, double mu_drift, double A,
                                   double sv, double s, double c,
                                   double lambda_g, double lambda_k,
                                   int n_gauss_nodes = 20, bool log_out = false,
                                   int kill_shape = 1) {
  const double f_decision = drdmswtn(t_adj, b, mu_drift, A, sv, s, 0.0, c,
                                     n_gauss_nodes, false, 1, false);
  const double F_decision = prdmswtn(t_adj, b, mu_drift, A, sv, s, 0.0, c,
                                     n_gauss_nodes, false, 1, false);
  return local_combo_response_pdf(t_adj, f_decision, F_decision,
                                  lambda_g, lambda_k, kill_shape, log_out);
}

inline double prdmswtn_local_combo(double t_adj, double b, double mu_drift, double A,
                                   double sv, double s, double c,
                                   double lambda_g, double lambda_k,
                                   int n_gauss_nodes = 20, bool log_out = false,
                                   int kill_shape = 1) {
  if (!(t_adj > 0.0)) return log_out ? R_NegInf : 0.0;
  const double rate_scale = lambda_g + lambda_k;
  const auto density_fn = [&](double u) {
    return drdmswtn_local_combo(u, b, mu_drift, A, sv, s, c,
                                lambda_g, lambda_k, n_gauss_nodes, false, kill_shape);
  };
  double out = emc2_isfinite(t_adj)
    ? integrate_density_gl20_finite(t_adj, density_fn)
    : integrate_density_gl20_infinite(rate_scale, density_fn);
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline double dgbm_local_combo(double t, double b, double mu, double sigma, double A,
                               double lambda_g, double lambda_k,
                               bool log_out = false, int kill_shape = 1) {
  const double f_decision = dgbm(t, b, mu, sigma, A, 0.0, false, 1, false);
  const double F_decision = pgbm(t, b, mu, sigma, A, 0.0, false, 1, false);
  return local_combo_response_pdf(t, f_decision, F_decision,
                                  lambda_g, lambda_k, kill_shape, log_out);
}

inline double pgbm_local_combo(double t, double b, double mu, double sigma, double A,
                               double lambda_g, double lambda_k,
                               bool log_out = false, int kill_shape = 1) {
  if (!(t > 0.0)) return log_out ? R_NegInf : 0.0;
  const double rate_scale = lambda_g + lambda_k;
  const auto density_fn = [&](double u) {
    return dgbm_local_combo(u, b, mu, sigma, A, lambda_g, lambda_k, false, kill_shape);
  };
  double out = emc2_isfinite(t)
    ? integrate_density_gl20_finite(t, density_fn)
    : integrate_density_gl20_infinite(rate_scale, density_fn);
  out = std::fmax(0.0, std::fmin(1.0, out));
  return log_out ? ((out > 0.0) ? std::log(out) : R_NegInf) : out;
}

inline double dswtn_core(double t_adj, double threshold, double mu_drift,
                         double sv, double s, double lambda, double c,
                         double log_prob_gt_c, bool log_out, int kill_shape, bool guess) {
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
  const double term_log_int = pnorm_std((c - mu_new) / sigma_new, false, true);

  const double log_f_hit = term_log_threshold + term_log_denom +
                         (-log_prob_gt_c) + term_log_exp + term_log_int +
                         erlang_log_surv(t_adj, lambda, kill_shape);

  if (!guess || lambda <= 0.0) {
    if (ISNAN(log_f_hit)) return log_out ? R_NegInf : 0.0;
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // Guess density: f_K * S_R
  const double log_sk = erlang_log_pdf(t_adj, lambda, kill_shape);
  const double log_sr = std::log1p(-std::max(0.0, std::min(1.0, pswtn(t_adj, threshold, mu_drift, sv, s, 0.0, c, false, 1, false))));
  const double log_f_guess = log_sk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  if (ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
  return log_out ? log_pdf : std::exp(log_pdf);
}

inline double pswtn_killed_quad(double t_adj, double threshold, double mu_drift,
                                double sv, double s, double lambda, double c,
                                int n_gauss_nodes = 20, int kill_shape = 1) {
  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return 0.0;

  double integral = 0.0;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);  // map [-1,1] -> [0,1]
    const double p = alpha + (1.0 - alpha) * u;
    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    integral += gl_weights[j] * pwald(t_adj, threshold, drift_j, s, 0.0, lambda, false, kill_shape);
  }

  return 0.5 * integral;
}

inline double pswtn_killed_inf_quad(double threshold, double mu_drift,
                                    double sv, double s, double lambda, double c,
                                    int n_gauss_nodes = 20,
                                    bool log_out = false, int kill_shape = 1) {
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
    return finish_log(log_hit);
  }

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);
    double p = alpha + (1.0 - alpha) * u;
    p = std::fmax(0.0, std::fmin(std::nextafter(1.0, 0.0), p));

    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    const double nu_j = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
    double log_hit_j = threshold * (drift_j - nu_j) / s2;
    if (kill_shape >= 2 && nu_j > FPM_EPSILON)
      log_hit_j += std::log1p(lambda * threshold / nu_j);

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  return finish_log(log_integral - M_LN2);
}

inline double log_uniform_threshold_hit(double eta, double b, double A) {
  if (A <= 1e-10) {
    const double d = std::fmax(b, 0.0);
    return eta * d;
  }

  const double d_lo_raw = b - A;
  const double d_hi_raw = b;

  if (d_hi_raw <= 0.0) return 0.0;

  const double zero_width = std::fmax(0.0, std::fmin(d_hi_raw, 0.0) - d_lo_raw);
  const double d_lo = std::fmax(d_lo_raw, 0.0);
  const double d_hi = d_hi_raw;

  double log_pos_integral = R_NegInf;
  if (d_hi > d_lo) {
    if (std::abs(eta) < 1e-12) {
      log_pos_integral = std::log(d_hi - d_lo);
    } else if (eta < 0.0) {
      log_pos_integral = log_diff_exp(eta * d_lo, eta * d_hi) - std::log(-eta);
    } else {
      log_pos_integral = log_diff_exp(eta * d_hi, eta * d_lo) - std::log(eta);
    }
  }

  double log_total = log_pos_integral;
  if (zero_width > 0.0) {
    log_total = log_sum_exp(log_total, std::log(zero_width));
  }

  return log_total - std::log(A);
}

inline double prdmswtn_killed_inf_quad(double b, double mu_drift, double A,
                                       double sv, double s, double lambda, double c,
                                       int n_gauss_nodes = 20,
                                       bool log_out = false, int kill_shape = 1) {
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
    if (kill_shape >= 2) {
      // Erlang-2: pwald(Inf,...) handles SPV via GL20
      return finish_log(pwald(R_PosInf, b, mu_drift, s, A, lambda, true, kill_shape));
    }
    const double nu = std::sqrt(mu_drift * mu_drift + 2.0 * s2 * lambda);
    const double eta = (mu_drift - nu) / s2;
    return finish_log(log_uniform_threshold_hit(eta, b, A));
  }

  const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
  if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;

  const int n_nodes = std::max(1, n_gauss_nodes);
  Rcpp::List gl = get_gl_nodes_weights(n_nodes);
  const Rcpp::NumericVector gl_nodes   = gl["nodes"];
  const Rcpp::NumericVector gl_weights = gl["weights"];

  double log_integral = R_NegInf;
  for (int j = 0; j < n_nodes; ++j) {
    const double u = 0.5 * (gl_nodes[j] + 1.0);
    double p = alpha + (1.0 - alpha) * u;
    p = std::fmax(0.0, std::fmin(std::nextafter(1.0, 0.0), p));

    const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
    double log_hit_j;
    if (kill_shape >= 2) {
      // Erlang-2: pwald(Inf, ...) with SPV via GL20
      log_hit_j = pwald(R_PosInf, b, drift_j, s, A, lambda, true, kill_shape);
    } else {
      const double nu_j  = std::sqrt(drift_j * drift_j + 2.0 * s2 * lambda);
      const double eta_j = (drift_j - nu_j) / s2;
      log_hit_j = log_uniform_threshold_hit(eta_j, b, A);
    }

    if (gl_weights[j] > 0.0) {
      log_integral = log_sum_exp(log_integral, std::log(gl_weights[j]) + log_hit_j);
    }
  }

  return finish_log(log_integral - M_LN2);
}

// [[Rcpp::export]]
double pswtn(double t_adj, double threshold, double mu_drift,
             double sv, double s = 1.0, double lambda = 0.0, double c = 0.0,
             bool log_out = false, int kill_shape = 1, bool guess = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;
  if (sv < 0.0) return R_NaN;

  if (guess && lambda > 0.0) {
    auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    const double fw = std::max(0.0, std::min(1.0, pswtn(t_adj, threshold, mu_drift, sv, s, 0.0, c, false, 1, false)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t_adj, lambda, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  if (sv <= 1e-10) {
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, threshold, mu_drift, s, 0.0, lambda, log_out, kill_shape, guess);
  }

  if (lambda > 1e-10) {
    if (!emc2_isfinite(t_adj)) {
      return pswtn_killed_inf_quad(threshold, mu_drift, sv, s, lambda, c, 20, log_out, kill_shape);
    }
    const double cdf = pswtn_killed_quad(t_adj, threshold, mu_drift, sv, s, lambda, c, 20, kill_shape);
    if (!(cdf > 0.0)) return log_out ? R_NegInf : 0.0;
    if (cdf >= 1.0) return log_out ? 0.0 : 1.0;
    return log_out ? std::log(cdf) : cdf;
  }

  const double v   = sv * sv;
  const double s2  = s * s;
  const double tv  = t_adj * v;
  const double th2 = threshold * threshold;

  // P(xi > c) for xi ~ N(mu_drift, sv^2)
  const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
  if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;

  // Shared denominator term
  const double s_denom_new = std::sqrt(t_adj * (s2 + tv));
  const double rho_new     = (std::sqrt(t_adj) * sv) / std::sqrt(s2 + tv);

  // Term 1: bivariate normal
  const double h1 = (mu_drift * t_adj - threshold) / s_denom_new;
  const double k1 = (mu_drift - c) / sv;
  double term1;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h1) > 8.0 || std::fabs(k1) > 8.0)
    term1 = norm_cdf_2d(h1, k1, rho_new);
  else
    term1 = norm_cdf_2d_fast(h1, k1, rho_new);
  const double log_A = (term1 > 1e-300) ? std::log(term1) : R_NegInf;

  // Term 2: reflected bivariate normal
  const double log_exp_term = (2.0 * threshold * mu_drift / s2) +
                               (2.0 * th2 * v / (s2 * s2));
  const double mu_p  = mu_drift + 2.0 * threshold * v / s2;
  const double h2    = (-mu_p * t_adj - threshold) / s_denom_new;
  const double k2    = (mu_p - c) / sv;
  double term2;
  if (std::fabs(rho_new) > 0.9999 || std::fabs(h2) > 8.0 || std::fabs(k2) > 8.0)
    term2 = norm_cdf_2d(h2, k2, -rho_new);
  else
    term2 = norm_cdf_2d_fast(h2, k2, -rho_new);
  const double log_B = (term2 > 1e-300) ? (log_exp_term + std::log(term2)) : R_NegInf;

  const double log_numerator = log_sum_exp(log_A, log_B);
  const double log_cdf       = log_numerator - log_prob_gt_c;

  if (ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  double cdf = std::exp(log_cdf);
  if (cdf < 0.0) {
    return log_out ? R_NegInf : 0.0;
  } else if (cdf > 1.0){
    return log_out ? 0.0 : 1.0;
  }
  return log_out ? std::log(cdf) : cdf;
}

// --------------------------------------------------------------------------
// Full RDMSWTN: SWTN + start-point variability.
// b  = upper threshold; threshold ~ Unif(b-A, b).
// Uses pre-cached GL20 nodes for the quadrature over [b-A, b].
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double drdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double lambda = 0.0, double c = 0.0,
                int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1, bool guess = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? R_NegInf : 0.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (dwald handles A=0 internally)
    if (mu_drift <= 1e-10) return log_out ? R_NegInf : 0.0;
    return dwald(t_adj, b, mu_drift, s, A, lambda, log_out, kill_shape, guess);
  } else if (no_A && !no_sv) {
    // SWTN with fixed threshold b
    return dswtn(t_adj, b, mu_drift, sv, s, lambda, c, log_out, kill_shape, guess);
  } else {
    // Full model: integrate dswtn_core (mixture if guess=true) over threshold ~ Unif(b-A, b).
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = get_gl_nodes_weights(n_nodes);
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    const double log_prob_gt_c = pnorm_std((c - mu_drift) / sv, false, true);
    if (emc2_isinf(log_prob_gt_c) && log_prob_gt_c < 0.0) return log_out ? R_NegInf : 0.0;
    // Map nodes from [-1,1] to [b-A, b]
    const double center     = b - 0.5 * A;
    const double half_width = 0.5 * A;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double thresh_j = center + half_width * gl_nodes[j];
      integral += gl_weights[j] * dswtn_core(
        t_adj, thresh_j, mu_drift, sv, s, lambda, c, log_prob_gt_c, false, kill_shape, guess
      );
    }
    const double out_val = integral * 0.5;
    if (out_val < 0.0 || ISNAN(out_val)) return log_out ? R_NegInf : 0.0;
    return log_out ? std::log(out_val) : out_val;
  }
}

// [[Rcpp::export]]
double prdmswtn(double t_adj, double b, double mu_drift, double A,
                double sv, double s = 1.0, double lambda = 0.0, double c = 0.0,
                int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1, bool guess = false) {
  if (t_adj <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? 0.0 : 1.0;

  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv < 1e-7);

  if (!emc2_isfinite(t_adj) && lambda > 1e-10) {
    return prdmswtn_killed_inf_quad(b, mu_drift, A, sv, s, lambda, c, n_gauss_nodes, log_out, kill_shape);
  }

  if (guess && lambda > 0.0) {
     auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    const double fw = std::max(0.0, std::min(1.0, prdmswtn(t_adj, b, mu_drift, A, sv, s, 0.0, c, n_gauss_nodes, false, 1, false)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t_adj, lambda, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  if (no_sv) {
    // sv=0: standard Wald with start-point variability A (pwald handles A=0 internally)
    if (mu_drift <= 1e-7) return log_out ? R_NegInf : 0.0;
    return pwald(t_adj, b, mu_drift, s, A, lambda, log_out, kill_shape, guess);
  } else if (no_A && !no_sv) {
    return pswtn(t_adj, b, mu_drift, sv, s, lambda, c, log_out, kill_shape, guess);
  } else {
    const int n_nodes = std::max(1, n_gauss_nodes);
    Rcpp::List gl = get_gl_nodes_weights(n_nodes);
    const Rcpp::NumericVector gl_nodes   = gl["nodes"];
    const Rcpp::NumericVector gl_weights = gl["weights"];
    const double alpha = pnorm_std((c - mu_drift) / sv, true, false);
    if (!(alpha < 1.0)) return log_out ? R_NegInf : 0.0;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double u = 0.5 * (gl_nodes[j] + 1.0);
      const double p = alpha + (1.0 - alpha) * u;
      const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
      integral += gl_weights[j] * pwald(t_adj, b, drift_j, s, A, lambda, false, kill_shape, guess);
    }
    double out_val = 0.5 * integral;
    out_val = std::fmax(0.0, std::fmin(1.0, out_val));
    return log_out ? std::log(out_val) : out_val;
  }
}

// --------------------------------------------------------------------------
// Vectorised R-callable exports (match zachdev naming so R code is portable).
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector t0, NumericVector sv,
                       NumericVector s = 1.0, NumericVector c = 0.0,
                       NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    if (lg > 0.0 && lk > 0.0) {
      pdf[i] = drdmswtn_local_combo(t_adj, pick(b, i), pick(v, i), pick(A, i),
                                    pick(sv, i), pick(s, i), pick(c, i),
                                    lg, lk, n_gauss_nodes, log_out, kill_shape);
    } else if (lg > 0.0) {
      pdf[i] = drdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                        pick(sv, i), pick(s, i), lg, pick(c, i),
                        n_gauss_nodes, log_out, kill_shape, true);
    } else {
      pdf[i] = drdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                        pick(sv, i), pick(s, i), lk, pick(c, i),
                        n_gauss_nodes, log_out, kill_shape, false);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector t0, NumericVector sv,
                       NumericVector s = 1.0, NumericVector c = 0.0,
                       NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    if (lg > 0.0 && lk > 0.0) {
      cdf[i] = prdmswtn_local_combo(t_adj, pick(b, i), pick(v, i), pick(A, i),
                                    pick(sv, i), pick(s, i), pick(c, i),
                                    lg, lk, n_gauss_nodes, log_out, kill_shape);
    } else if (lg > 0.0) {
      cdf[i] = prdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                        pick(sv, i), pick(s, i), lg, pick(c, i),
                        n_gauss_nodes, log_out, kill_shape, true);
    } else {
      cdf[i] = prdmswtn(t_adj, pick(b, i), pick(v, i), pick(A, i),
                        pick(sv, i), pick(s, i), lk, pick(c, i),
                        n_gauss_nodes, log_out, kill_shape, false);
    }
  }
  return cdf;
}

// --------------------------------------------------------------------------
// GBM vectorized wrappers.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector dGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0, NumericVector s = 1.0,
                      NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                      bool log_out = false, int kill_shape = 1) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    if (lg > 0.0 && lk > 0.0) {
      pdf[i] = dgbm_local_combo(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                                lg, lk, log_out, kill_shape);
    } else if (lg > 0.0) {
      pdf[i] = dgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                    lg, log_out, kill_shape, true);
    } else {
      pdf[i] = dgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                    lk, log_out, kill_shape, false);
    }
  }
  return pdf;
}

// [[Rcpp::export]]
NumericVector pGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0, NumericVector s = 1.0,
                      NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                      bool log_out = false, int kill_shape = 1) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t_adj = t[i] - pick(t0, i);
    if (t_adj <= 0.0) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    if (lg > 0.0 && lk > 0.0) {
      cdf[i] = pgbm_local_combo(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                                lg, lk, log_out, kill_shape);
    } else if (lg > 0.0) {
      cdf[i] = pgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                    lg, log_out, kill_shape, true);
    } else {
      cdf[i] = pgbm(t_adj, pick(b, i), pick(v, i), pick(s, i), pick(A, i),
                    lk, log_out, kill_shape, false);
    }
  }
  return cdf;
}

// --------------------------------------------------------------------------
// Vectorised adapters for the EMC2 race-likelihood machinery.
// Column layout (from Ttransform + p_types reorder):
//   v=0, B=1, A=2, t0=3, s=4, sv=5  [pContaminant and b follow, ignored here]
// Scaling convention: absorb s inside kernel; pass s=1 to drdmswtn/prdmswtn.
// --------------------------------------------------------------------------

#endif
