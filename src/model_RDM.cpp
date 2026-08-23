// [[Rcpp::export]]
double pwald(double t, double mu, double b, double A = 0.0, double sigma = 1.0,
             double t0 = 0.0, double lambda_g = 0.0, double lambda_k = 0.0, bool log_out = false,
             int kill_shape = 1, bool guess = false, bool posdrift = true,
             double erlang_omega = 1.0);
// [[Rcpp::export]]
double dwald(double t, double mu, double b, double A = 0.0, double sigma = 1.0,
             double t0 = 0.0, double lambda_g = 0.0, double lambda_k = 0.0, bool log_out = false,
             int kill_shape = 1, bool guess = false, bool posdrift = true,
             double erlang_omega = 1.0);
#include "model_RDM.h"
double dwald(double t, double mu, double b, double A, double sigma,
             double t0, double lambda_g, double lambda_k, bool log_out,
             int kill_shape, bool guess, bool posdrift,
             double erlang_omega) {
  if (kill_shape == 3) {
    const double e1 = dwald(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, guess, posdrift, 1.0);
    const double e2 = dwald(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam is the EAM-adjusted time (rt - t0).
  // Erlang processes start at physical t=0 and use raw t; EAM uses t_eam.
  if (t <= 0.0 || sigma <= 0.0 || lambda_g < 0.0 || lambda_k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  if (posdrift && mu <= 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  const double t_eam = t - t0;
  // Combined local kill+guess uses the dedicated *_local_combo helpers.
  // This kernel handles a single active Erlang branch selected by `guess`.
  const double lambda = guess ? lambda_g : lambda_k;
  const double k_eff = (lambda <= 1e-8 ? 0.0 : lambda);

  // When rt < t0 the EAM hasn't started: f_EAM = 0, S_R = 1.
  // Only the erlang guess component contributes.
  if (t_eam <= FPM_EPSILON) {
    if (!guess || k_eff <= 0.0) return log_out ? R_NegInf : 0.0;
    const double log_fG = erlang_log_pdf(t, k_eff, kill_shape);
    return log_out ? log_fG : std::exp(log_fG);
  }

  double x_lo = 0.0;
  double x_hi = x_lo + A;

  const double var = sigma * sigma * t_eam;
  if (var <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  if (x_hi < x_lo) std::swap(x_lo, x_hi);
  const double span = x_hi - x_lo;
  double log_f_hit = R_NegInf;

  // 1. Evidence-based hit density (f_R * S_K)
  if (span <= FPM_EPSILON) {
    const double d = b - x_hi;
    if (d <= 0.0) {
      if (t_eam == 0.0) {
        return log_out ? R_PosInf : std::numeric_limits<double>::infinity();
      }
      return log_out ? R_NegInf : 0.0;
    }
    const double delta = d - mu * t_eam;
    log_f_hit = std::log(d) - std::log(t_eam) + Gstar(var, delta, true);
    if (guess || k_eff > 0.0) {
      log_f_hit += erlang_log_surv(t, k_eff, kill_shape);
    } 
  } else {
    // SPV form
    const double mu_new = b - mu * t_eam;
    const double log_pdf_hi = Gstar(var, x_hi - mu_new, true);
    const double log_pdf_lo = Gstar(var, x_lo - mu_new, true);
    const double log_cdf_integral = Gstar_Integral(var, mu_new, x_lo, x_hi, true);

    signed_log term1 = make_signed_log(R_NegInf, 0);
    if (std::abs(mu) > FPM_EPSILON && log_cdf_integral != R_NegInf) {
      term1 = make_signed_log(std::log(std::abs(mu)) + std::log(t_eam) + log_cdf_integral, mu > 0.0 ? 1 : -1);
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
      log_f_hit = total.log_abs - std::log(t_eam) - std::log(span);
      if (guess || k_eff > 0.0) {
        log_f_hit+= erlang_log_surv(t, k_eff, kill_shape);
      }
    }
  }

  if (!guess || k_eff <= 0.0) {
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // 2. Guess-based density (f_K * S_R)
  // Erlang uses raw t; EAM survivor uses t_eam (k=0 → no kill in survivor)
  const double log_sk = erlang_log_pdf(t, k_eff, kill_shape);
  const double log_sr = log_surv_from_log_cdf(
      pwald(t, mu, b, A, sigma, t0, 0.0, 0.0, true, 1, false, posdrift));
  const double log_f_guess = log_sk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}
// [[Rcpp::export]]
NumericVector dWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false) {
  int n = t.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res[i] = dwald(t[i], v[i], B[i] + A[i], A[i], 1.0, t0[i], 0.0, 0.0, log_out, 1, false);
  }
  return res;
}
// [[Rcpp::export]]
NumericVector pWald(NumericVector t, NumericVector v,
                    NumericVector B, NumericVector A, NumericVector t0,
                    bool log_out = false) {
  int n = t.size();
  NumericVector res(n);
  for (int i = 0; i < n; i++) {
    res[i] = pwald(t[i], v[i], B[i] + A[i], A[i], 1.0, t0[i], 0.0, 0.0, log_out, 1, false);
  }
  return res;
}
// --------------------------------------------------------------------------
// Wald / killed Wald CDF with uniform start-point variability.
// Returns sub-CDF:
//   P(T_wald <= t AND T_kill > T_wald)
// k = 0 gives standard Wald, nested.
// t = Inf returns total finite-response probability.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double pgbm(double t, double mu, double b, double A, double sigma,
            double t0, double lambda_g, double lambda_k, bool log_out,
            int kill_shape, bool guess, double erlang_omega) {
  if (kill_shape == 3) {
    const double e1 = pgbm(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, guess, 1.0);
    const double e2 = pgbm(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
  auto finish = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;
    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;
    return log_out ? log_p : std::exp(log_p);
  };

  if (t <= 0.0 || !(b > FPM_EPSILON) || !(sigma > FPM_EPSILON) || !(A >= 0.0) || lambda_g < 0.0 || lambda_k < 0.0) {
    return finish(R_NegInf);
  }
  if (1.0 + A > b + FPM_EPSILON) return finish(R_NegInf);

  const double t_eam = t - t0;
  const double lambda = guess ? lambda_g : lambda_k;
  const double k_eff = (lambda <= 1e-8 ? 0.0 : lambda);

  if (guess && k_eff > 0.0) {
    // Mixture CDF: F_i(t) = 1 - S_R,i(t_eam) * S_K,i(t)
    const double log_sk = erlang_log_surv(t, k_eff, kill_shape);
    if (t_eam <= FPM_EPSILON) {
      return finish(std::log1p(-std::exp(log_sk)));
    }
    const double fw = std::max(0.0, std::min(1.0, pgbm(t, mu, b, A, sigma, t0, 0.0, 0.0, false, 1, false)));
    const double log_sr = std::log1p(-fw);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  // No-guess path: return 0 if EAM hasn't started.
  if (t_eam <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  // For t0 > 0, the Erlang kill process starts at 0 while EAM starts at t0.
  // Analytic: P(T_EAM + t0 < T_kill AND T_EAM + t0 < t)
  if (!guess && k_eff > 1e-10 && t0 > 1e-10) {
    const double log_sk_t0 = -k_eff * t0;
    const double F1 = pgbm(t - t0, mu, b, A, sigma, 0.0, 0.0, k_eff, false, 1, false);
    if (kill_shape <= 1) {
      return finish(log_sk_t0 + std::log(std::fmax(0.0, F1)));
    }
    const double F2 = pgbm(t - t0, mu, b, A, sigma, 0.0, 0.0, k_eff, false, 2, false);
    // Same shifted Erlang-2 identity as pwald.
    const double val = std::exp(log_sk_t0) * (k_eff * t0 * F1 + F2);
    return finish(std::log(std::fmax(0.0, val)));
  }

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double alpha = mu - 0.5 * sigma * sigma;  // log-space drift
  const double sig2 = sigma * sigma;
  const double norm_const = A;

  // Tilted log-space drift for kill (= |alpha| when k = 0).
  const double ak     = std::sqrt(alpha * alpha + 2.0 * sig2 * k_eff);
  const double c_eta  = (alpha - ak) / sig2;   // eta1 in log-space
  const double q      = 1.0 - c_eta;           // Jacobian-corrected coefficient (1 - eta1)
  const double k_log  = 2.0 * ak / sig2;       // 2*ak/sigma^2

  // ---- t = Inf: total finite-response probability ----
  auto log_p_inf = [&]() -> double {
    if (!(norm_const > FPM_EPSILON)) {
      // Point-start: P(hit before kill) = b^c_eta, Erlang-2 adds log1p correction.
      double log_hit = c_eta * log_b;
      if (kill_shape >= 2 && ak > FPM_EPSILON)
        log_hit += std::log1p(k_eff * log_b / ak);
      return log_hit;
    }

    if (kill_shape >= 2) {
      // P(hit before Erlang-2 kill) =
      // exp(c_eta*log_b)/A * [I0 + (k/ak)*I1],
      // I0 = int exp(qy)dy, I1 = int (log_b-y)exp(qy)dy, y in [0, x_hi].
      const double log_prefactor = c_eta * log_b - std::log(norm_const);
      double I0, I1;
      if (std::abs(q) <= 1e-8) {
        I0 = x_hi;
        I1 = log_b * x_hi - 0.5 * x_hi * x_hi;
      } else {
        const double u = q * x_hi;
        I0 = std::expm1(u) / q;
        const double J = (std::exp(u) * (u - 1.0) + 1.0) / (q * q);  // int_0^x y exp(qy) dy
        I1 = log_b * I0 - J;
      }
      const double p = std::exp(log_prefactor) * (I0 + (k_eff / ak) * I1);
      if (p > 0.0 && emc2_isfinite(p)) {
        if (p >= 1.0) return 0.0;
        return std::log(p);
      }
      // exp(log_prefactor) underflowed or exp(u) overflowed while the finite
      // probability itself is representable: rebuild in signed-log space.
      if (std::abs(q) <= 1e-8) return R_NegInf;
      const double u = q * x_hi;
      const double log_I0 = ((u > 0.0) ? log_diff_exp(u, 0.0)
                                       : log_diff_exp(0.0, u)) -
                            std::log(std::fabs(q));
      // u = q * x_hi shares q's sign (x_hi > 0), so I0 = expm1(u)/q > 0.
      const signed_log sI0 = make_signed_log(log_I0, 1);
      // J = (exp(u)(u-1) + 1) / q^2
      const signed_log sJ_num = signed_log_add(signed_log_product(u - 1.0, u),
                                               make_signed_log(0.0, 1));
      signed_log sJ = sJ_num;
      if (sJ.sign != 0)
        sJ = make_signed_log(sJ.log_abs - 2.0 * std::log(std::fabs(q)), sJ.sign);
      signed_log sI1 = make_signed_log(R_NegInf, 0);
      if (sI0.sign != 0 && log_b > 0.0)
        sI1 = make_signed_log(std::log(log_b) + sI0.log_abs, sI0.sign);
      sI1 = signed_log_sub(sI1, sJ);
      signed_log total = sI0;
      if (sI1.sign != 0) {
        total = signed_log_add(total, make_signed_log(
            std::log(k_eff) - std::log(ak) + sI1.log_abs, sI1.sign));
      }
      if (total.sign <= 0 || total.log_abs == R_NegInf) return R_NegInf;
      const double log_p = log_prefactor + total.log_abs;
      return log_p > 0.0 ? 0.0 : log_p;
    }

    // Erlang-1 SPV analytic:
    // exp(c_eta * log_b) / A * integral_0^{log1p(A)} exp(q * y) dy
    const double log_prefactor = c_eta * log_b - std::log(norm_const);
    if (std::abs(q) < 1e-8)
      return log_prefactor + std::log(x_hi);   // q = 0: integral = x_hi
    if (q > 0.0)
      return log_prefactor + std::log(std::expm1(q * x_hi)) - std::log(q);
    return log_prefactor + std::log(-std::expm1(q * x_hi)) - std::log(-q);
  };

  if (!emc2_isfinite(t_eam)) {
    return finish(log_p_inf());
  }

  const double st = sigma * std::sqrt(t_eam);
  if (st <= FPM_EPSILON) return log_out ? R_NegInf : 0.0;

  const double a  = 1.0 / st;
  const double c1 = (ak * t_eam - log_b) / st;   // c for Phi in term1
  const double c2 = (-ak * t_eam - log_b) / st;  // c for Phi in term2

  // ---- Point-start: A = 0 ----
  if (!(norm_const > FPM_EPSILON)) {
    const double d = log_b - x_hi;  // = log_b when A = 0

    if (d <= 0.0) return log_out ? 0.0 : 1.0;

    const double log_prefactor = c_eta * d;

    const double log_cdf1    = pnorm_std((ak * t_eam - d) / st, true, true);
    const double log_cdf2    = pnorm_std((-ak * t_eam - d) / st, true, true);
    const double log_exp_t   = k_log * d;
    const double log_cdf_ak  = log_sum_exp(log_cdf1, log_exp_t + log_cdf2);

    if (kill_shape >= 2 && k_eff > FPM_EPSILON && ak > FPM_EPSILON) {
      const double log_b2  = log_exp_t + log_cdf2;
      const double log_mw  = (log_cdf1 > log_b2) ? log_diff_exp(log_cdf1, log_b2) : R_NegInf;
      const double log_w   = std::log(k_eff) + std::log(d) - std::log(ak);
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
        std::log(k_eff) - std::log(ak) + Ddiff.log_abs, Ddiff.sign);
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
  //
  // Natural-scale fast path. Both integrals close with the identity
  //   ∫ exp(k y) Phi(a y + c) dy
  //     = [exp(k y) Phi(a y + c) - E(k,c) Phi(a y + c - k/a)] / k,
  //   E(k,c) = exp(k^2/(2 a^2) - k c / a),
  // and the two shifted arguments coincide: with a = 1/st,
  // c1 - c2 = 2*ak*t_eam/st = k_log*st, so
  //   (a y + c2) - (q - k_log)/a = a y + c1 - q/a.
  // That leaves 6 natural Phi evaluations instead of the 8 log-scale ones
  // the generic helper would take, and drops the signed-log machinery.
  // Rejected (and handed to the signed-log path below) on overflow or on
  // cancellation past EMC2_NAT_REL_CANCEL.
  const double k2 = q - k_log;
  if (std::fabs(q) > 1e-8 && std::fabs(k2) > 1e-8) {
    const double sA_hi = a * x_hi + c1;
    const double sA_lo = a * x_lo + c1;
    const double sB_hi = a * x_hi + c2;
    const double sB_lo = a * x_lo + c2;
    const double sS_hi = sA_hi - q * st;   // shared shifted argument
    const double sS_lo = sA_lo - q * st;

    const double PhiA_hi = pnorm_std(sA_hi, true, false);
    const double PhiA_lo = pnorm_std(sA_lo, true, false);
    const double PhiB_hi = pnorm_std(sB_hi, true, false);
    const double PhiB_lo = pnorm_std(sB_lo, true, false);
    const double PhiS_hi = pnorm_std(sS_hi, true, false);
    const double PhiS_lo = pnorm_std(sS_lo, true, false);

    const double E1 = std::exp(0.5 * q * q * st * st - q * c1 * st);
    const double E2 = std::exp(0.5 * k2 * k2 * st * st - k2 * c2 * st);

    const double t1a =  std::exp(q * x_hi) * PhiA_hi / q;
    const double t1b = -std::exp(q * x_lo) * PhiA_lo / q;
    const double t1c = -E1 * PhiS_hi / q;
    const double t1d =  E1 * PhiS_lo / q;
    const double I1  = (t1a + t1b) + (t1c + t1d);

    const double t2a =  std::exp(k2 * x_hi) * PhiB_hi / k2;
    const double t2b = -std::exp(k2 * x_lo) * PhiB_lo / k2;
    const double t2c = -E2 * PhiS_hi / k2;
    const double t2d =  E2 * PhiS_lo / k2;
    const double I2  = (t2a + t2b) + (t2c + t2d);

    const double mag1 = std::fabs(t1a) + std::fabs(t1b) +
                        std::fabs(t1c) + std::fabs(t1d);
    const double mag2 = std::fabs(t2a) + std::fabs(t2b) +
                        std::fabs(t2c) + std::fabs(t2d);

    const double W1 = std::exp(c_eta * log_b);
    const double W2 = std::exp((c_eta + k_log) * log_b);
    const double cdf_val = (W1 * I1 + W2 * I2) / norm_const;

    const bool natural_ok =
        emc2_isfinite(mag1) && emc2_isfinite(mag2) &&
        emc2_isfinite(W1) && emc2_isfinite(W2) &&
        I1 > EMC2_NAT_REL_CANCEL * mag1 &&
        I2 > EMC2_NAT_REL_CANCEL * mag2 &&
        emc2_isfinite(cdf_val) && cdf_val > 0.0;

    if (natural_ok)
      return finish(std::log(cdf_val > 1.0 ? 1.0 : cdf_val));
  }

  const double log_term1 = c_eta * log_b +
    log_integrate_exp_times_normal_cdf(q, a, c1, x_lo, x_hi);
  const double log_term2 = (c_eta + k_log) * log_b +
    log_integrate_exp_times_normal_cdf(q - k_log, a, c2, x_lo, x_hi);

  const double log_cdf_val = log_sum_exp(log_term1, log_term2) - std::log(norm_const);

  return finish(log_cdf_val);
}
// [[Rcpp::export]]
double dgbm(double t, double mu, double b, double A = 0.0, double sigma = 1.0,
            double t0 = 0.0, double lambda_g = 0.0, double lambda_k = 0.0, bool log_out = false,
            int kill_shape = 1, bool guess = false, double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = dgbm(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, guess, 1.0);
    const double e2 = dgbm(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
  if (t <= 0.0 || !(sigma > FPM_EPSILON) || !(b > FPM_EPSILON) || !(A >= 0.0) || lambda_g < 0.0 || lambda_k < 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  if (1.0 + A > b + FPM_EPSILON) return log_out ? R_NegInf : 0.0;

  const double t_eam = t - t0;
  const double lambda = guess ? lambda_g : lambda_k;
  const double k_eff = (lambda <= 1e-8 ? 0.0 : lambda);

  // EAM hasn't started; only erlang guess contributes (S_R = 1).
  if (t_eam <= FPM_EPSILON) {
    if (!guess || k_eff <= 0.0) return log_out ? R_NegInf : 0.0;
    const double log_fG = erlang_log_pdf(t, k_eff, kill_shape);
    return log_out ? log_fG : std::exp(log_fG);
  }

  const double x_lo = 0.0;
  const double x_hi = std::log1p(A);
  const double log_b = std::log(b);
  const double log_mu = mu - 0.5 * sigma * sigma;
  const double var = sigma * sigma * t_eam;
  const double norm_const = A;

  double log_f_hit = R_NegInf;

  // 1. Evidence hit density (f_R * S_K); erlang survival uses raw t
  if (!(norm_const > FPM_EPSILON)) {
    const double d = log_b - x_hi;  // = log_b when A = 0
    if (d <= 0.0) {
      log_f_hit = (t_eam == 0.0) ? std::numeric_limits<double>::infinity() : R_NegInf;
    } else {
      const double delta = d - log_mu * t_eam;
      log_f_hit = std::log(d) - std::log(t_eam) + Gstar(var, delta, true)
                  + erlang_log_surv(t, k_eff, kill_shape);
    }
  } else {
    // SPV: normalization for start density on [1, 1 + A] mapped to log-space.
    // Natural fast path first; the log evaluator below is authoritative when
    // exp(exp_factor_log) under/overflows or the two-term numerator cancels.
    const double mu_new_exp = log_b - log_mu * t_eam + var;
    const double exp_factor_log = log_b - log_mu * t_eam + 0.5 * var;
    const double pdf_hi = Gstar(var, x_hi - mu_new_exp);
    const double pdf_lo = Gstar(var, x_lo - mu_new_exp);
    const double cdf_integral = Gstar_Integral(var, mu_new_exp, x_lo, x_hi);

    const double integral_term1 = (log_mu * t_eam - var) * cdf_integral;
    const double integral_term2 = var * (pdf_hi - pdf_lo);
    const double integral_result = integral_term1 + integral_term2;
    const bool natural_ok =
        emc2_isfinite(exp_factor_log) && emc2_isfinite(integral_result) &&
        integral_result > 0.0 &&
        integral_result > EMC2_NAT_REL_CANCEL *
            (std::fabs(integral_term1) + std::fabs(integral_term2));
    double pdf_val = 0.0;
    if (natural_ok) {
      pdf_val = (std::exp(exp_factor_log) * integral_result) / (norm_const * t_eam);
    }
    if (natural_ok && pdf_val > 0.0 && emc2_isfinite(pdf_val)) {
      log_f_hit = std::log(pdf_val) + erlang_log_surv(t, k_eff, kill_shape);
    } else {
      // Signed-log evaluation: exp(exp_factor_log) is never materialised and
      // the numerator difference is formed with signed-log arithmetic.
      const double log_pdf_hi = Gstar(var, x_hi - mu_new_exp, true);
      const double log_pdf_lo = Gstar(var, x_lo - mu_new_exp, true);
      const double log_cdf_integral = Gstar_Integral(var, mu_new_exp, x_lo, x_hi, true);
      const signed_log term1 =
          signed_log_product(log_mu * t_eam - var, log_cdf_integral);
      signed_log term2 = signed_log_sub(make_signed_log(log_pdf_hi, 1),
                                        make_signed_log(log_pdf_lo, 1));
      if (term2.sign != 0)
        term2 = make_signed_log(std::log(var) + term2.log_abs, term2.sign);
      const signed_log total = signed_log_add(term1, term2);
      if (total.sign > 0 && total.log_abs != R_NegInf && !ISNAN(total.log_abs)) {
        log_f_hit = exp_factor_log + total.log_abs -
                    std::log(norm_const) - std::log(t_eam) +
                    erlang_log_surv(t, k_eff, kill_shape);
      }
    }
  }

  if (!guess || k_eff <= 0.0) {
    return log_out ? log_f_hit : std::exp(log_f_hit);
  }

  // 2. Guess hit density (f_K * S_R); erlang pdf uses raw t, EAM survivor uses t_eam
  const double log_fk = erlang_log_pdf(t, k_eff, kill_shape);
  const double log_sr = log_surv_from_log_cdf(
      pgbm(t, mu, b, A, sigma, t0, 0.0, 0.0, true, 1, false, 1.0));
  const double log_f_guess = log_fk + log_sr;

  const double log_pdf = log_sum_exp(log_f_hit, log_f_guess);
  return log_out ? log_pdf : std::exp(log_pdf);
}
// [[Rcpp::export]]
double dswtn(double t, double mu_drift, double threshold, double s = 1.0,
             double t0 = 0.0, double sv = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
             bool log_out = false, int kill_shape = 1, bool guess = false, bool posdrift = true,
             double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = dswtn(t, mu_drift, threshold, s, t0, sv, lambda_g, lambda_k, true, 1, guess, posdrift, 1.0);
    const double e2 = dswtn(t, mu_drift, threshold, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  if (threshold <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (sv < 0.0) return R_NaN;
  if (s <= 1e-10) return log_out ? R_NegInf : 0.0;
  if (posdrift && sv <= 1e-10 && mu_drift <= 0.0)
    return log_out ? R_NegInf : 0.0;

  const double t_eam = t - t0;
  // Combined local kill+guess is handled one level up. This SWTN kernel sees
  // only a single active Erlang branch selected by `guess`.
  const double lambda = guess ? lambda_g : lambda_k;

  // sv == 0: reduces to standard Wald (which handles t_eam <= 0 and neg drift internally)
  if (sv <= 1e-10) {
    return dwald(t, mu_drift, threshold, 0.0, s, t0,
                 guess ? lambda : 0.0, guess ? 0.0 : lambda,
                 log_out, kill_shape, guess, posdrift);
  }

  // Handle t_eam <= 0: EAM hasn't started; only erlang guess contributes.
  if (t_eam <= 1e-10) {
    if (!guess || lambda <= 0.0) return log_out ? R_NegInf : 0.0;
    const double log_fG = erlang_log_pdf(t, lambda, kill_shape);
    return log_out ? log_fG : std::exp(log_fG);
  }

  if (posdrift && !guess && kill_shape <= 1 && lambda > 1e-10) {
    const double log_sk = erlang_log_surv(t, lambda, kill_shape);
    const double log_pdf0 = positive_trunc_swtn_density_k0(
      t, mu_drift, threshold, s, t0, sv, true
    );
    if (!R_FINITE(log_pdf0) || !R_FINITE(log_sk)) return log_out ? R_NegInf : 0.0;
    const double log_pdf = log_pdf0 + log_sk;
    return log_out ? log_pdf : std::exp(log_pdf);
  }

  if (posdrift) {
    return dswtn_positive_drift_quad(t, mu_drift, threshold, s, t0, sv,
                                     lambda_g, lambda_k, log_out, kill_shape, guess);
  }

  // log_norm = 0 means defective density; dswtn_core normalises by exp(log_norm).
  return dswtn_core(t_eam, mu_drift, threshold, s, t0, sv, lambda, 0.0,
                    log_out, kill_shape, guess, false);
}
// [[Rcpp::export]]
double pswtn(double t, double mu_drift, double threshold, double s = 1.0,
             double t0 = 0.0, double sv = 0.0, double lambda_g = 0.0, double lambda_k = 0.0,
             bool log_out = false, int kill_shape = 1, bool guess = false, bool posdrift = true,
             double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = pswtn(t, mu_drift, threshold, s, t0, sv, lambda_g, lambda_k, true, 1, guess, posdrift, 1.0);
    const double e2 = pswtn(t, mu_drift, threshold, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  const double lambda = guess ? lambda_g : lambda_k;
  if (sv < 0.0) return R_NaN;
  if (posdrift && sv <= 1e-10 && mu_drift <= 0.0)
    return log_out ? R_NegInf : 0.0;

  // t is raw rt; dt = t - t0 is EAM decision time.
  const double dt = t - t0;
  if (dt <= 1e-10) {
    // EAM not started; only erlang guess can contribute (S_R = 1).
    if (!guess || lambda <= 0.0) return log_out ? R_NegInf : 0.0;
    auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    return finish(std::log1p(-std::exp(erlang_log_surv(t, lambda, kill_shape))));
  }
  if (threshold <= 1e-10) return log_out ? 0.0 : 1.0;

  const double t_raw = t;  // raw rt for erlang

  if (guess && lambda > 0.0) {
    if (sv > 1e-10 && posdrift) {
      return pswtn_positive_drift_quad(t, mu_drift, threshold, s, t0, sv,
                                       lambda_g, lambda_k, log_out, kill_shape, true);
    }
    auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    const double fw = std::max(0.0, std::min(1.0, pswtn(t, mu_drift, threshold, s, t0, sv, 0.0, 0.0, false, 1, false, posdrift, 1.0)));
    const double log_sr = std::log1p(-fw);
    const double log_sk = erlang_log_surv(t_raw, lambda, kill_shape);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  if (!guess && lambda > 1e-10 && !emc2_isfinite(dt)) {
    return pswtn_killed_inf_quad(threshold, mu_drift, sv, s, lambda,
                                 20, log_out, kill_shape, posdrift);
  }

  if (sv <= 1e-10) {
    // Pass raw time and t0 to pwald so erlang uses physical time.
    return pwald(t_raw, mu_drift, threshold, 0.0, s, t0, lambda, lambda,
                 log_out, kill_shape, guess, posdrift);
  }

  if (posdrift) {
    return pswtn_positive_drift_quad(t, mu_drift, threshold, s, t0, sv,
                                     lambda_g, lambda_k, log_out, kill_shape, false);
  }

  if (lambda > 1e-10) {
    if (!emc2_isfinite(dt)) {
      return pswtn_killed_inf_quad(threshold, mu_drift, sv, s, lambda, 20, log_out, kill_shape, posdrift);
    }
    const double cdf = pswtn_killed_quad(dt, mu_drift, threshold, s, t0, sv, lambda, 20, kill_shape, guess, posdrift);
    if (!(cdf > 0.0)) return log_out ? R_NegInf : 0.0;
    if (cdf >= 1.0) return log_out ? 0.0 : 1.0;
    return log_out ? std::log(cdf) : cdf;
  }
  
  if (!emc2_isfinite(dt) && lambda <= 1e-10) {
    if (posdrift) return log_out ? 0.0 : 1.0;
    const double log_H = log_swtn_hit_mass(mu_drift, sv, s, threshold);
    return log_out ? log_H : std::exp(log_H);
  }
  
  const double sv2 = sv * sv;
  const double s2  = s * s;
  const double Q   = s2 + sv2 * dt;
  const double denom = std::sqrt(dt * Q);
  
  const double h1 = (mu_drift * dt - threshold) / denom;
  
  const double alpha = 2.0 * threshold / s2;
  const double mu_p = mu_drift + alpha * sv2;
  
  const double log_exp_term =
    alpha * mu_drift + 0.5 * alpha * alpha * sv2;
  
  const double h2 = (-mu_p * dt - threshold) / denom;
  
  const double log_A = pnorm_std(h1, true, true);
  const double log_B = log_exp_term + pnorm_std(h2, true, true);
  
  const double log_num = log_sum_exp(log_A, log_B);
  
  const double log_cdf = log_num;
  
  if (ISNAN(log_cdf)) return log_out ? R_NegInf : 0.0;
  if (log_cdf > 0.0){
    return log_out ? 0.0 : 1.0;
  }
  return log_out ? log_cdf : std::exp(log_cdf);
}
// --------------------------------------------------------------------------
// Full RDMSWTN: SWTN + start-point variability.
// b  = upper threshold; threshold ~ Unif(b-A, b).
// Uses pre-cached GL20 nodes for the quadrature over [b-A, b].
// --------------------------------------------------------------------------
// [[Rcpp::export]]
double rdmswtn_tt_qinv(double u, double tau) {
  if (!(tau > 0.0) || ISNAN(u) || u < 0.0) return NA_REAL;
  // Keep the inverse consistent with rdmswtn_tt_q() at tau = +Inf.
  if (!R_FINITE(tau)) return u;
  const double Q = 0.5 * tau;
  if (u > Q) return R_PosInf;
  if (u == 0.0) return 0.0;
  const double radicand = std::fmax(0.0, 1.0 - 2.0 * u / tau);
  return 2.0 * u / (1.0 + std::sqrt(radicand));
}
// [[Rcpp::export]]
double drdmswtn_tt(double t, double mu_drift, double b, double A,
                   double s = 1.0, double t0 = 0.0, double sv = 0.0,
                   double tau = 1.0, bool log_out = false,
                   bool posdrift = true) {
  if (!(tau > 0.0) || ISNAN(t) || ISNAN(t0))
    return log_out ? R_NegInf : 0.0;
  const double x = t - t0;
  if (!(x > 0.0) || (R_FINITE(tau) && x >= tau) || !R_FINITE(t))
    return log_out ? R_NegInf : 0.0;
  const double q = rdmswtn_tt_q(x, tau);
  const double log_f = drdmswtn(
      q, mu_drift, b, A, s, 0.0, sv, 0.0, 0.0, 20,
      true, 1, false, posdrift, 1.0);
  if (!R_FINITE(log_f)) return log_out ? R_NegInf : 0.0;
  const double log_jac = R_FINITE(tau) ? std::log1p(-x / tau) : 0.0;
  const double out = log_f + log_jac;
  return log_out ? out : std::exp(out);
}
// [[Rcpp::export]]
double prdmswtn_tt(double t, double mu_drift, double b, double A,
                   double s = 1.0, double t0 = 0.0, double sv = 0.0,
                   double tau = 1.0, bool log_out = false,
                   bool posdrift = true) {
  if (!(tau > 0.0) || ISNAN(t) || ISNAN(t0))
    return log_out ? R_NegInf : 0.0;
  const double x = t - t0;
  if (!(x > 0.0)) return log_out ? R_NegInf : 0.0;
  const double q = (R_FINITE(tau) && x >= tau)
    ? 0.5 * tau
    : rdmswtn_tt_q(x, tau);
  return prdmswtn(
      q, mu_drift, b, A, s, 0.0, sv, 0.0, 0.0, 20,
      log_out, 1, false, posdrift, 1.0);
}
// --------------------------------------------------------------------------
// Vectorised R-callable exports (match zachdev naming so R code is portable).
// --------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector dSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector s = 1.0, NumericVector t0 = 0.0,
                       NumericVector sv = 0.0,
                       NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1,
                       bool posdrift = true, NumericVector erlang_omega = 1.0) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t0_i = pick(t0, i);
    const double dt = t[i] - t0_i;
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    const bool erl = (lg > 0.0 || lk > 0.0);
    if (t[i] <= 0.0 || (dt <= 0.0 && !erl)) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    if (lg > 0.0 && lk > 0.0) {
      pdf[i] = drdmswtn_local_combo(t[i], pick(v, i), pick(b, i), pick(A, i),
                                    pick(s, i), t0_i, pick(sv, i),
                                    lg, lk, n_gauss_nodes, log_out, kill_shape, posdrift, omega);
    } else if (lg > 0.0) {
      pdf[i] = drdmswtn(t[i], pick(v, i), pick(b, i), pick(A, i),
                        pick(s, i), t0_i, pick(sv, i), lg, 0.0,
                        n_gauss_nodes, log_out, kill_shape, true, posdrift, omega);
    } else {
      pdf[i] = drdmswtn(t[i], pick(v, i), pick(b, i), pick(A, i),
                        pick(s, i), t0_i, pick(sv, i), 0.0, lk,
                        n_gauss_nodes, log_out, kill_shape, false, posdrift, omega);
    }
  }
  return pdf;
}
// [[Rcpp::export]]
NumericVector pSWTNspv(NumericVector t, NumericVector v, NumericVector b,
                       NumericVector A, NumericVector s = 1.0, NumericVector t0 = 0.0,
                       NumericVector sv = 0.0,
                       NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                       int n_gauss_nodes = 20, bool log_out = false, int kill_shape = 1,
                       bool posdrift = true, NumericVector erlang_omega = 1.0) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t0_i = pick(t0, i);
    const double dt = t[i] - t0_i;
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    const bool erl = (lg > 0.0 || lk > 0.0);
    if (t[i] <= 0.0 || (dt <= 0.0 && !erl)) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    if (lg > 0.0 && lk > 0.0) {
      cdf[i] = prdmswtn_local_combo(t[i], pick(v, i), pick(b, i), pick(A, i),
                                    pick(s, i), t0_i, pick(sv, i),
                                    lg, lk, n_gauss_nodes, log_out, kill_shape, posdrift, omega);
    } else if (lg > 0.0) {
      cdf[i] = prdmswtn(t[i], pick(v, i), pick(b, i), pick(A, i),
                        pick(s, i), t0_i, pick(sv, i), lg, 0.0,
                        n_gauss_nodes, log_out, kill_shape, true, posdrift, omega);
    } else {
      cdf[i] = prdmswtn(t[i], pick(v, i), pick(b, i), pick(A, i),
                        pick(s, i), t0_i, pick(sv, i), 0.0, lk,
                        n_gauss_nodes, log_out, kill_shape, false, posdrift, omega);
    }
  }
  return cdf;
}
// --------------------------------------------------------------------------
// GBM vectorized wrappers.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector dGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0 = 0.0, NumericVector s = 1.0,
                      NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                      bool log_out = false, int kill_shape = 1,
                      NumericVector erlang_omega = 1.0) {
  const int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t0_i = pick(t0, i);
    const double dt = t[i] - t0_i;
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    const bool erl = (lg > 0.0 || lk > 0.0);
    if (t[i] <= 0.0 || (dt <= 0.0 && !erl)) { pdf[i] = log_out ? R_NegInf : 0.0; continue; }
    if (lg > 0.0 && lk > 0.0) {
      pdf[i] = dgbm_local_combo(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                                t0_i, lg, lk, log_out, kill_shape, omega);
    } else if (lg > 0.0) {
      pdf[i] = dgbm(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                    t0_i, lg, 0.0, log_out, kill_shape, true, omega);
    } else {
      pdf[i] = dgbm(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                    t0_i, 0.0, lk, log_out, kill_shape, false, omega);
    }
  }
  return pdf;
}
// [[Rcpp::export]]
NumericVector pGBMspv(NumericVector t, NumericVector v, NumericVector b,
                      NumericVector A, NumericVector t0 = 0.0, NumericVector s = 1.0,
                      NumericVector lambda_g = 0.0, NumericVector lambda_k = 0.0,
                      bool log_out = false, int kill_shape = 1,
                      NumericVector erlang_omega = 1.0) {
  const int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    const double t0_i = pick(t0, i);
    const double dt = t[i] - t0_i;
    const double lg = pick(lambda_g, i);
    const double lk = pick(lambda_k, i);
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    const bool erl = (lg > 0.0 || lk > 0.0);
    if (t[i] <= 0.0 || (dt <= 0.0 && !erl)) { cdf[i] = log_out ? R_NegInf : 0.0; continue; }
    if (lg > 0.0 && lk > 0.0) {
      cdf[i] = pgbm_local_combo(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                                t0_i, lg, lk, log_out, kill_shape, omega);
    } else if (lg > 0.0) {
      cdf[i] = pgbm(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                    t0_i, lg, 0.0, log_out, kill_shape, true, omega);
    } else {
      cdf[i] = pgbm(t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
                    t0_i, 0.0, lk, log_out, kill_shape, false, omega);
    }
  }
  return cdf;
}
double pwald(double t, double mu, double b, double A, double sigma,
             double t0, double lambda_g, double lambda_k, bool log_out, int kill_shape, bool guess, bool posdrift,
             double erlang_omega) {
  if (kill_shape == 3) {
    const double e1 = pwald(t, mu, b, A, sigma, t0, lambda_g, lambda_k, true, 1, guess, posdrift, 1.0);
    const double e2 = pwald(t, mu, b, A, sigma, t0, 2.0 * lambda_g, 2.0 * lambda_k, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is the EAM-adjusted time.
  // Erlang processes use raw t; EAM uses t_eam.
  auto finish = [&](double log_p) {
    if (ISNAN(log_p)) return NA_REAL;

    if (log_p > 0.0 && log_p < 1e-10) log_p = 0.0;
    if (log_p > 0.0) log_p = 0.0;

    return log_out ? log_p : std::exp(log_p);
  };

  if (sigma <= 0.0 || lambda_g < 0.0 || lambda_k < 0.0) {
    return NA_REAL;
  }
  if (posdrift && mu <= 0.0) {
    return log_out ? R_NegInf : 0.0;
  }
  // Combined local kill+guess uses the dedicated *_local_combo helpers.
  // This kernel handles a single active Erlang branch selected by `guess`.
  const double lambda = guess ? lambda_g : lambda_k;
  const double k_eff = (lambda <= 1e-8 ? 0.0 : lambda);

  if (t <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  const double t_eam = t - t0;

  if (guess && k_eff > 0.0) {
    // Mixture CDF: F_i(t) = 1 - S_R,i(t_eam) * S_K,i(t)
    // When t_eam <= 0, EAM hasn't started: S_R = 1, CDF = 1 - S_K(t)
    const double log_sk = erlang_log_surv(t, k_eff, kill_shape);
    if (t_eam <= FPM_EPSILON) {
      return finish(std::log1p(-std::exp(log_sk)));
    }
    const double fw = std::max(0.0, std::min(1.0, pwald(t, mu, b, A, sigma, t0, 0.0, 0.0, false, 1, false, posdrift)));
    const double log_sr = std::log1p(-fw);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  // No-guess path: EAM CDF only; return 0 if EAM hasn't started.
  if (t_eam <= FPM_EPSILON) {
    return log_out ? R_NegInf : 0.0;
  }

  // For t0 > 0, the Erlang kill process starts at 0 while EAM starts at t0.
  // Analytic: P(T_EAM + t0 < T_kill AND T_EAM + t0 < t)
  // For Erlang-1 (Exp): S_K(u) = exp(-k*u) = exp(-k*t0) * exp(-k*v), where v = u-t0
  // For Erlang-2: S_K(u) = (1 + k*u)*exp(-k*u) = exp(-k*t0) * [ (1 + k*t0)*exp(-k*v) + k*v*exp(-k*v) ]
  if (!guess && k_eff > 1e-10 && t0 > 1e-10) {
    const double log_sk_t0 = -k_eff * t0; // Exp part of S_K(t0)
    const double F1 = pwald(t - t0, mu, b, A, sigma, 0.0, 0.0, k_eff, false, 1, false, posdrift);
    if (kill_shape <= 1) {
      return finish(log_sk_t0 + std::log(std::fmax(0.0, F1)));
    }
    const double F2 = pwald(t - t0, mu, b, A, sigma, 0.0, 0.0, k_eff, false, 2, false, posdrift);
    // Erlang-2 with shifted t0:
    // exp(-k*t0) * [(1 + k*t0) * B0 + k * B1], with F2 = B0 + k*B1.
    // So use exp(-k*t0) * [k*t0 * B0 + F2].
    const double val = std::exp(log_sk_t0) * (k_eff * t0 * F1 + F2);
    return finish(std::log(std::fmax(0.0, val)));
  }

  double x_lo = 0.0;
  double x_hi = x_lo + A;
  if (x_hi < x_lo) std::swap(x_lo, x_hi);

  const double span = x_hi - x_lo;
  const double sig2 = sigma * sigma;
  // Fast path: no kill — use canonical Wald SPV for finite t_eam.
  if (k_eff <= 0.0) {
    if (emc2_isfinite(t_eam)) {
      const double cdf = pwald_k0(t_eam, b, mu, A, sigma);
      const double cl = std::max(0.0, std::min(1.0, cdf));
      if (!log_out) return cl;
      if (cl > 0.0) return std::log(cl);
      // Natural CDF underflowed (early times / strong negative drift): the
      // log-output request still has a representable value in log space.
      return wald_k0_log_cdf(t_eam, b, mu, A, sigma);
    }
    // t_eam = Inf with k=0: fall through — log_p_inf() handles it.
  }

  const double nu = std::sqrt(mu * mu + 2.0 * sig2 * k_eff);
  const double eta1 = (mu - nu) / sig2;

  auto log_p_inf = [&]() -> double {
    if (span <= FPM_EPSILON) {
      const double d = b - x_hi;
      if (d <= 0.0) return 0.0;
      const double log_p1 = eta1 * d;
      if (kill_shape <= 1 || nu <= FPM_EPSILON) return log_p1;
      // Erlang-2 eventual hit: exp(eta1*d) * (1 + lambda*d/nu)
      return log_p1 + std::log1p(k_eff * d / nu);
    }

    if (kill_shape >= 2 && k_eff > FPM_EPSILON) {
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
      const double p = (I0 + (k_eff / nu) * I1) / span;
      if (p > 0.0 && p < 1.0) return std::log(p);
      if (p >= 1.0) return 0.0;
      // exp(eta1*d) underflowed (eta1 < 0 for k > 0): the eventual-hit mass
      // is still representable in log space.  Rebuild the same integral with
      // signed-log antiderivatives:
      //   I0 = (exp(eta1 d_lo) - exp(eta1 d_hi)) / (-eta1)
      //   I1 = [exp(eta1 d)(d/eta1 - 1/eta1^2)]_{d_lo}^{d_hi}.
      if (std::abs(eta1) <= FPM_EPSILON) return R_NegInf;
      const double log_I0 = log_diff_exp(eta1 * d_lo, eta1 * d_hi) -
                            std::log(-eta1);
      const signed_log a_hi = signed_log_product(
          d_hi / eta1 - 1.0 / (eta1 * eta1), eta1 * d_hi);
      const signed_log a_lo = signed_log_product(
          d_lo / eta1 - 1.0 / (eta1 * eta1), eta1 * d_lo);
      const signed_log sI1 = signed_log_sub(a_hi, a_lo);
      signed_log total = make_signed_log(log_I0, 1);
      if (sI1.sign != 0) {
        total = signed_log_add(total, make_signed_log(
            std::log(k_eff) - std::log(nu) + sI1.log_abs, sI1.sign));
      }
      if (total.sign <= 0 || total.log_abs == R_NegInf) return R_NegInf;
      const double log_p = total.log_abs - std::log(span);
      return log_p > 0.0 ? 0.0 : log_p;
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

  if (!emc2_isfinite(t_eam)) {
    return finish(log_p_inf());
  }

  const double st = sigma * std::sqrt(t_eam);
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

    const double term1_arg = (nu * t_eam - d) / st;
    const double term2_arg = (-nu * t_eam - d) / st;

    const double log_cdf1 = pnorm_std(term1_arg, true, true);
    const double log_cdf2 = pnorm_std(term2_arg, true, true);

    const double log_exp_term = 2.0 * nu * d / sig2;

    const double log_cdf_nu =
      log_sum_exp(log_cdf1, log_exp_term + log_cdf2);

    if (kill_shape >= 2 && k_eff > FPM_EPSILON && nu > FPM_EPSILON) {
      // Erlang-2: add lambda * M_W term.
      // M_W (unnorm) = Phi(z1) - exp_term*Phi(z2), which equals C^{-1} * dL/d(-lambda).
      const double log_b2 = log_exp_term + log_cdf2;
      const double log_mw = (log_cdf1 > log_b2)
                              ? log_diff_exp(log_cdf1, log_b2)
                              : R_NegInf;
      const double log_w = std::log(k_eff) + std::log(d) - std::log(nu);  // log(lambda*d/nu)
      const double log_cdf_e2 = (log_mw > R_NegInf)
                                  ? log_sum_exp(log_cdf_nu, log_w + log_mw)
                                  : log_cdf_nu;
      return finish(log_prefactor + log_cdf_e2);
    }

    return finish(log_prefactor + log_cdf_nu);
  }

  if (kill_shape >= 2 && k_eff > FPM_EPSILON) {
    // SPV + Erlang-2: analytic via slog_int_{eta,d_eta}_pnorm.
    // total = (T1 + T2 + (k/nu)*(D1-D2)) / span
    if (b <= x_hi) return NA_REAL;
    const double a    = 1.0 / st;
    const double c1   = (nu * t_eam - b) / st;
    const double c2   = (-nu * t_eam - b) / st;
    const double eta2 = (mu + nu) / sig2;
    signed_log T1 = slog_int_eta_pnorm(eta1, b, a, c1, x_lo, x_hi);
    signed_log T2 = slog_int_eta_pnorm(eta2, b, a, c2, x_lo, x_hi);
    signed_log D1 = slog_int_d_eta_pnorm(eta1, b, a, c1, x_lo, x_hi);
    signed_log D2 = slog_int_d_eta_pnorm(eta2, b, a, c2, x_lo, x_hi);
    signed_log total  = signed_log_add(T1, T2);
    signed_log Ddiff  = signed_log_sub(D1, D2);
    if (Ddiff.sign != 0 && Ddiff.log_abs != R_NegInf) {
      signed_log extra = make_signed_log(
        std::log(k_eff) - std::log(nu) + Ddiff.log_abs, Ddiff.sign);
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
  // Phi((x + nu*t_eam - b) / st)
  const double c1 = (nu * t_eam - b) / st;

  const double log_term1 =
    eta1 * b +
    log_integrate_exp_times_normal_cdf(-eta1, a, c1, x_lo, x_hi);

  // Second term:
  // exp[(mu + nu)(b-x)/sig2] *
  // Phi((x - nu*t_eam - b) / st)
  const double eta2 = (mu + nu) / sig2;
  const double c2 = (-nu * t_eam - b) / st;

  const double log_term2 =
    eta2 * b +
    log_integrate_exp_times_normal_cdf(-eta2, a, c2, x_lo, x_hi);

  const double log_cdf_val =
    log_sum_exp(log_term1, log_term2) - std::log(span);

  return finish(log_cdf_val);
}
// [[Rcpp::export]]
double drdmswtn(double t, double mu_drift, double b, double A,
                double s = 1.0, double t0 = 0.0, double sv = 0.0,
                double lambda_g = 0.0, double lambda_k = 0.0,
                int n_gauss_nodes = 20, bool log_out = false,
                int kill_shape = 1, bool guess = false, bool posdrift = true,
                double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = drdmswtn(t, mu_drift, b, A, s, t0, sv, lambda_g, lambda_k, n_gauss_nodes, true, 1, guess, posdrift, 1.0);
    const double e2 = drdmswtn(t, mu_drift, b, A, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, n_gauss_nodes, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? R_NegInf : 0.0;

  const double t_eam = t - t0;
  const double lambda = guess ? lambda_g : lambda_k;
  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv <= 1e-10);
  if (posdrift && no_sv && mu_drift <= 0.0)
    return log_out ? R_NegInf : 0.0;

  if (no_sv) {
    // sv=0: standard Wald under the caller's posdrift semantics.
    return dwald(t, mu_drift, b, A, s, t0, lambda_g, lambda_k,
                 log_out, kill_shape, guess, posdrift, erlang_omega);
  } else if (posdrift && !guess && kill_shape <= 1 && lambda > 1e-10) {
    const double log_sk = erlang_log_surv(t, lambda, kill_shape);
    double log_pdf0;
    if (no_A) {
      log_pdf0 = positive_trunc_swtn_density_k0(t, mu_drift, b, s, t0, sv, true);
    } else {
      log_pdf0 = drdmswtn_joint_A_sv_density_postrunc(t_eam, mu_drift, b, A, s, sv, true);
    }
    if (!R_FINITE(log_pdf0) || !R_FINITE(log_sk)) return log_out ? R_NegInf : 0.0;
    const double log_pdf = log_pdf0 + log_sk;
    return log_out ? log_pdf : std::exp(log_pdf);
  } else if (posdrift && !guess && lambda <= 1e-10) {
    if (no_A) {
      return positive_trunc_swtn_density_k0(t, mu_drift, b, s, t0, sv, log_out);
    }
    return drdmswtn_joint_A_sv_density_postrunc(t_eam, mu_drift, b, A, s, sv, log_out);
  } else if (posdrift) {
    return drdmswtn_positive_drift_quad(t, mu_drift, b, A, s, t0, sv,
                                        lambda_g, lambda_k,
                                        log_out, kill_shape, guess);
  } else if (no_A && !no_sv) {
    // SWTN with fixed threshold b (dswtn handles t_eam <= 0 via its t0 param).
    return dswtn(t, mu_drift, b, s, t0, sv,
                 guess ? lambda : 0.0, guess ? 0.0 : lambda,
                 log_out, kill_shape, guess, posdrift);
  } else {
    if (!guess && lambda <= 1e-10) {
      return drdmswtn_joint_A_sv_density_fullgauss(t_eam, mu_drift, b, A, s, sv, false, n_gauss_nodes, log_out);
    }

    // Full model: integrate defective dswtn_core over threshold ~ Unif(b-A, b),
    // then normalise by the SPV hit mass if posdrift=true.
    const int n_nodes = std::max(1, n_gauss_nodes);
    const GLRule& gl = gl_get_rule(n_nodes);
    const std::vector<double>& gl_nodes   = gl.x;
    const std::vector<double>& gl_weights = gl.w;
    const double center     = b - 0.5 * A;
    const double half_width = 0.5 * A;
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double thresh_j = center + half_width * gl_nodes[j];
      integral += gl_weights[j] * dswtn_core(
        t_eam, mu_drift, thresh_j, s, t0, sv, lambda, 0.0, false, kill_shape, guess, false
      );
    }
    double out_val = integral * 0.5;
    if (out_val > 0.0 && !ISNAN(out_val)) return log_out ? std::log(out_val) : out_val;
    if (out_val < 0.0 || ISNAN(out_val)) return log_out ? R_NegInf : 0.0;
    // All nodes underflowed on the natural scale: log-space accumulation.
    double log_acc = R_NegInf;
    for (int j = 0; j < n_nodes; ++j) {
      if (!(gl_weights[j] > 0.0)) continue;
      const double thresh_j = center + half_width * gl_nodes[j];
      const double lf = dswtn_core(t_eam, mu_drift, thresh_j, s, t0, sv,
                                   lambda, 0.0, true, kill_shape, guess, false);
      if (lf == R_NegInf || ISNAN(lf)) continue;
      log_acc = log_sum_exp(log_acc, std::log(gl_weights[j]) + lf);
    }
    const double log_pdf = log_acc - M_LN2;
    if (!(log_pdf > R_NegInf) || ISNAN(log_pdf)) return log_out ? R_NegInf : 0.0;
    return return_from_log(log_pdf, log_out);
  }
}
// [[Rcpp::export]]
double prdmswtn(double t, double mu_drift, double b, double A,
                double s = 1.0, double t0 = 0.0, double sv = 0.0,
                double lambda_g = 0.0, double lambda_k = 0.0,
                int n_gauss_nodes = 20, bool log_out = false,
                int kill_shape = 1, bool guess = false, bool posdrift = true,
                double erlang_omega = 1.0) {
  if (kill_shape == 3) {
    const double e1 = prdmswtn(t, mu_drift, b, A, s, t0, sv, lambda_g, lambda_k, n_gauss_nodes, true, 1, guess, posdrift, 1.0);
    const double e2 = prdmswtn(t, mu_drift, b, A, s, t0, sv, 2.0 * lambda_g, 2.0 * lambda_k, n_gauss_nodes, true, 2, guess, posdrift, 0.0);
    return mix_erlang12(e1, e2, erlang_omega, log_out);
  }
  // t is raw rt; t_eam = t - t0 is EAM-adjusted time. Erlang uses raw t.
  if (t <= 0.0) return log_out ? R_NegInf : 0.0;
  if (b <= 1e-7) return log_out ? 0.0 : 1.0;

  const double t_eam = t - t0;
  const double lambda = guess ? lambda_g : lambda_k;
  const bool no_A  = (A  < 1e-7);
  const bool no_sv = (sv <= 1e-10);
  if (posdrift && no_sv && mu_drift <= 0.0)
    return log_out ? R_NegInf : 0.0;

  if (posdrift && !guess && lambda <= 1e-10 && !no_sv) {
    return prdmswtn_joint_A_sv_cdf_postrunc(t, mu_drift, b, A, s, t0, sv,
                                            n_gauss_nodes, log_out);
  }

  if (!guess && lambda > 1e-10 && !emc2_isfinite(t_eam) && !no_sv) {
    return prdmswtn_killed_inf_quad(b, mu_drift, A, sv, s, t0, lambda,
                                    n_gauss_nodes, log_out, kill_shape, posdrift);
  }

  if (posdrift && !no_sv) {
    return prdmswtn_positive_drift_quad(t, mu_drift, b, A, s, t0, sv,
                                        lambda_g, lambda_k,
                                        log_out, kill_shape, guess);
  }

  if (guess && lambda > 0.0) {
    auto finish = [&](double log_p) {
      if (ISNAN(log_p)) return NA_REAL;
      if (log_p > 0.0) log_p = 0.0;
      return log_out ? log_p : std::exp(log_p);
    };
    // Erlang uses raw t; EAM survivor uses t_eam (via prdmswtn with t0=0 and t_eam as raw t)
    const double log_sk = erlang_log_surv(t, lambda, kill_shape);
    if (t_eam <= 1e-10) {
      // EAM hasn't started: S_R = 1, CDF = 1 - S_K(t)
      return finish(std::log1p(-std::exp(log_sk)));
    }
    const double fw = std::max(0.0, std::min(1.0, prdmswtn(t_eam, mu_drift, b, A, s, 0.0, sv, 0.0, 0.0, n_gauss_nodes, false, 1, false, false)));
    const double log_sr = std::log1p(-fw);
    return finish(std::log1p(-std::exp(log_sr + log_sk)));
  }

  // No-guess path: return 0 if EAM hasn't started.
  if (t_eam <= 1e-10) return log_out ? R_NegInf : 0.0;

  if (no_sv) {
    // sv=0: standard Wald under the caller's posdrift semantics.
    // Pass raw t and t0 to pwald so erlang inside uses physical time.
    return pwald(t, mu_drift, b, A, s, t0, lambda, lambda,
                 log_out, kill_shape, guess, posdrift);
  }

  if (!emc2_isfinite(t_eam) && lambda <= 1e-10) {
    // At t=Inf with no killing the defective mass is exactly the hit probability.
    const double log_mass = no_sv
        ? pwald(R_PosInf, mu_drift, b, A, s, 0.0, 0.0, 0.0,
                true, 1, false, false)
        : no_A
            ? log_swtn_hit_mass(mu_drift, sv, s, b)
            : log_swtn_spv_hit_mass_full(mu_drift, sv, s, b, A, n_gauss_nodes);
    return log_out ? log_mass : std::exp(log_mass);
  }

  if (no_A && !no_sv) {
    // pswtn receives raw t and uses t0 internally to form dt.
    return pswtn(t, mu_drift, b, s, t0, sv, lambda, lambda, log_out, kill_shape, guess, posdrift);
  } else {
    const int n_nodes = std::max(1, n_gauss_nodes);
    const GLRule& gl = gl_get_rule(n_nodes);
    const std::vector<double>& gl_nodes   = gl.x;
    const std::vector<double>& gl_weights = gl.w;
    
    double integral = 0.0;
    for (int j = 0; j < n_nodes; ++j) {
      const double u = 0.5 * (gl_nodes[j] + 1.0);
      
      const double p = std::fmin(std::nextafter(1.0, 0.0),
                                 std::fmax(1e-15, u));
      
      const double drift_j =
        mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
      
      integral += gl_weights[j] * pwald(
        t, drift_j, b, A, s, t0,
        lambda, lambda, false, kill_shape, guess,
        false // defective fixed-drift kernel inside mixture
      );
    }
    double out_val = 0.5 * integral;

    out_val = std::fmax(0.0, std::fmin(1.0, out_val));
    if (out_val > 0.0 || !log_out) return log_out ? std::log(out_val) : out_val;
    // Log output requested and every node underflowed: log-space accumulation.
    double log_acc = R_NegInf;
    for (int j = 0; j < n_nodes; ++j) {
      if (!(gl_weights[j] > 0.0)) continue;
      const double u = 0.5 * (gl_nodes[j] + 1.0);
      const double p = std::fmin(std::nextafter(1.0, 0.0), std::fmax(1e-15, u));
      const double drift_j = mu_drift + sv * R::qnorm(p, 0.0, 1.0, true, false);
      const double lc = pwald(t, drift_j, b, A, s, t0,
                              lambda, lambda, true, kill_shape, guess, false);
      if (lc == R_NegInf || ISNAN(lc)) continue;
      log_acc = log_sum_exp(log_acc, std::log(gl_weights[j]) + lc);
    }
    const double log_cdf = log_acc - M_LN2;
    if (!(log_cdf > R_NegInf) || ISNAN(log_cdf)) return R_NegInf;
    return std::fmin(log_cdf, 0.0);
  }
}
// [[Rcpp::export]]
NumericVector dRDMSWTN_TT_cpp(
    NumericVector t, NumericVector v, NumericVector b, NumericVector A,
    NumericVector s, NumericVector t0, NumericVector sv, NumericVector tau,
    bool log_out = false, bool posdrift = true) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    out[i] = drdmswtn_tt(
      t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
      pick(t0, i), pick(sv, i), pick(tau, i), log_out, posdrift);
  }
  return out;
}
// [[Rcpp::export]]
NumericVector pRDMSWTN_TT_cpp(
    NumericVector t, NumericVector v, NumericVector b, NumericVector A,
    NumericVector s, NumericVector t0, NumericVector sv, NumericVector tau,
    bool log_out = false, bool posdrift = true) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; ++i) {
    out[i] = prdmswtn_tt(
      t[i], pick(v, i), pick(b, i), pick(A, i), pick(s, i),
      pick(t0, i), pick(sv, i), pick(tau, i), log_out, posdrift);
  }
  return out;
}
