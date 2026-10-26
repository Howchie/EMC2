#include <Rcpp.h>
#include <complex>
#include "model_RDMSWTN_exwald.h"
#include "faddeeva.h"
#include "wald_functions.h"
#include "col_registry.h"

// ============================================================
// RDMSWTN with exponential hesitation (c_name RDMSWTN_EXWALD)
// Column layout: v=0 (drift), B=1, A=2, t0=3, s=4, tau=5 (mean hesitation).
// X = D + H, D the Wald first passage with distance uniform on [B, B + A],
// H ~ Exp(mean tau) independent.  No drift variability, no clocks.
//
// With lambda = 1/tau and kappa^2 = v^2 - 2 s^2 lambda, tilting the Wald
// density by e^(lambda u) and averaging over the distance a gives, in s = 1
// units with E(a) = -(a - v t)^2 / (2 t) and y(c) = (a - c t) / sqrt(2 t),
//   f_X = [G(B + A) - G(B)] / A,
//   G(a) = e^E { [(v + kappa) erfcx(y(kappa)) + (v - kappa) erfcx(y(-kappa))] / 4
//               - v erfcx(y(v)) / 2 },
// and for a point start f_X = lambda e^E [erfcx(y(kappa)) + erfcx(y(-kappa))] / 2.
// When kappa = i omega the erfcx pair is w(xi), conj(w(xi)) with
// xi = (omega t + i a) / sqrt(2 t) in the upper half plane, so G and f_X are
// real parts of one Faddeeva evaluation.  Both branches are the same analytic
// function of kappa^2.  Survivor and CDF follow from the density:
//   S_X = S_D + tau f_X,  F_X = F_D - tau f_X.
// ============================================================

namespace {

namespace ecol = emc2col::rdmswtn_exwald;

constexpr double kLogUnderflow = -746.0;

// log erfcx(y), all real y
inline double log_erfcx(double y) {
  static const double log_sqrt_2_over_pi = 0.5 * std::log(2.0 / M_PI);
  if (y >= 0.0) return log_sqrt_2_over_pi + std::log(mills_ratio_std(M_SQRT2 * y));
  const double r = std::sqrt(2.0 / M_PI) * mills_ratio_std(-M_SQRT2 * y);
  return y * y + std::log(2.0 - std::exp(-y * y) * r);
}

// log(e^E erfcx(y)) given delta = E + y^2 in a cancellation-free form.  For
// y < 0 both E and y^2 can exceed 1e40 when v / s is large, so their sum is
// taken from delta and erfcx(y) e^(-y^2) = erfc(y) = 1 + erf(-y).
inline double log_exp_erfcx(double E, double delta, double y) {
  if (y >= 0.0) return E + log_erfcx(y);
  return delta + std::log1p(std::erf(-y));
}

struct Tilt {
  double lam, log_lam;
  bool cplx;
  double kap, v_minus_kap;  // real branch
  double om;                // complex branch: kappa = i om
};

inline Tilt make_tilt(double v, double lam) {
  Tilt k{};
  k.lam = lam;
  k.log_lam = std::log(lam);
  const double k2 = v * v - 2.0 * lam;
  k.cplx = k2 < 0.0;
  if (k.cplx) {
    k.om = std::sqrt(-k2);
  } else {
    k.kap = std::sqrt(k2);
    k.v_minus_kap = 2.0 * lam / (v + k.kap);  // v > 0: no cancellation
  }
  return k;
}

// E + y(sgn kappa)^2 = (v - c)(2a - (v + c) t) / 2 with c = sgn kappa, real
// branch; v - kappa is taken from its cancellation-free form.
inline double kappa_delta(double t, double a, double v, const Tilt& k, int sgn) {
  const double vmc = (sgn > 0) ? k.v_minus_kap : v + k.kap;
  const double vpc = (sgn > 0) ? v + k.kap : k.v_minus_kap;
  return 0.5 * vmc * (2.0 * a - vpc * t);
}

// log f_X for a point start at distance a >= 0 (s = 1 units)
inline double point_log_pdf(double t, double a, double v, const Tilt& k) {
  const double st2 = std::sqrt(2.0 * t);
  const double E = -(a - v * t) * (a - v * t) / (2.0 * t);
  if (!k.cplx) {
    const double y1 = (a - k.kap * t) / st2, y2 = (a + k.kap * t) / st2;
    return k.log_lam - M_LN2 +
           log_sum_exp(log_exp_erfcx(E, kappa_delta(t, a, v, k, 1), y1),
                       log_exp_erfcx(E, kappa_delta(t, a, v, k, -1), y2));
  }
  const double re = faddeeva::w(std::complex<double>(k.om * t / st2, a / st2)).real();
  return (re > 0.0) ? k.log_lam + E + std::log(re) : R_NegInf;
}

// G(a) of the start-point average (s = 1 units); max_log collects the
// largest component for the cancellation test.
inline signed_log antideriv(double t, double a, double v, const Tilt& k, double& max_log) {
  const double st2 = std::sqrt(2.0 * t);
  const double E = -(a - v * t) * (a - v * t) / (2.0 * t);
  const signed_log neg =
      make_signed_log(std::log(0.5 * v) + log_exp_erfcx(E, 0.0, (a - v * t) / st2), -1);
  signed_log pos;
  if (!k.cplx) {
    const double l1 = std::log(0.25 * (v + k.kap)) +
                      log_exp_erfcx(E, kappa_delta(t, a, v, k, 1), (a - k.kap * t) / st2);
    const double l2 = std::log(0.25 * k.v_minus_kap) +
                      log_exp_erfcx(E, kappa_delta(t, a, v, k, -1), (a + k.kap * t) / st2);
    pos = make_signed_log(log_sum_exp(l1, l2), 1);
  } else {
    const std::complex<double> w = faddeeva::w(std::complex<double>(k.om * t / st2, a / st2));
    const double v_re = 0.5 * v * w.real(), om_im = 0.5 * k.om * w.imag();
    const double c = v_re - om_im;
    max_log = std::fmax(max_log, E + std::log(std::fmax(std::fabs(v_re), std::fabs(om_im))));
    pos = signed_from_real(0.0);
    if (c != 0.0) pos = make_signed_log(E + std::log(std::fabs(c)), c > 0.0 ? 1 : -1);
  }
  max_log = std::fmax(max_log, std::fmax(pos.log_abs, neg.log_abs));
  return signed_log_add(pos, neg);
}

// log f_X (s = 1 units), distance uniform on [B, B + A]
inline double exwald_log_pdf_s1(double t, double v, double B, double A, double lam) {
  const Tilt k = make_tilt(v, lam);
  if (A <= 1e-12) return point_log_pdf(t, B, v, k);
  double max_log = R_NegInf;
  const signed_log hi = antideriv(t, B + A, v, k, max_log);
  const signed_log lo = antideriv(t, B, v, k, max_log);
  const signed_log total = signed_log_sub(hi, lo);
  if (total.sign > 0 && !ISNAN(total.log_abs) &&
      total.log_abs > max_log + EMC2_LOG_CANCEL_MIN)
    return total.log_abs - std::log(A);
  // The bracket cancelled (narrow A, or lambda t small at early times): average
  // the positive point density over the distance by panelled quadrature.
  return wald_k0_log_avg_over_d(B, B + A, [&](double a) {
    return point_log_pdf(t, a, v, k);
  }) - std::log(A);
}

// Hesitation-free limit: the Wald kernels of RDMSWTN (s = 1 units).
inline double wald_log(double t, double v, double B, double A, char what) {
  const double b = B + A;
  if (what == 'f') {
    double pdf;
    if (dwald_k0_natural(t, b, v, A, pdf)) return (pdf > 0.0) ? std::log(pdf) : R_NegInf;
    return dwald_k0_log(t, b, v, A);
  }
  return (what == 'F') ? wald_k0_log_cdf(t, b, v, A) : wald_k0_log_surv(t, b, v, A);
}

// what: 'f' density, 'F' CDF, 'S' survivor of the decision time tt.
inline double exwald_log(double tt, double v, double B, double A, double s, double tau,
                         char what) {
  if (!(tt > 0.0)) return (what == 'S') ? 0.0 : R_NegInf;
  if (!(v > 0.0) || !(s > 0.0) || !(B >= 0.0) || !(A >= 0.0) || !(B + A > 0.0) ||
      !(tau >= 0.0)) return R_NaN;
  if (!R_FINITE(tt)) return (what == 'F') ? 0.0 : R_NegInf;
  const double inv_s = 1.0 / s;
  const double vs = v * inv_s, Bs = B * inv_s, As = A * inv_s;
  if (!(tau > 0.0) || !R_FINITE(1.0 / tau)) return wald_log(tt, vs, Bs, As, what);
  if (!R_FINITE(tau)) return (what == 'S') ? 0.0 : R_NegInf;
  const double lam = 1.0 / tau;
  const double lf = exwald_log_pdf_s1(tt, vs, Bs, As, lam);
  if (what == 'f') return lf;
  const double l_tf = std::log(tau) + lf;
  if (what == 'S') return std::fmin(log_sum_exp(wald_log(tt, vs, Bs, As, 'S'), l_tf), 0.0);
  const double lFD = wald_log(tt, vs, Bs, As, 'F');
  if (lFD == R_NegInf) return R_NegInf;
  const double ratio = l_tf - lFD;  // log(tau f_X / F_D) < 0
  if (ratio < -1e-6) return lFD + log1m_exp(ratio);
  // F_D - tau f_X lost its digits (lambda tt small):
  //   F_X = int_0^tt f_D(u) (1 - e^(-lambda (tt - u))) du <= F_D.
  return std::fmin(lFD, wald_k0_log_avg_over_d(0.0, tt, [&](double u) {
    if (!(u > 0.0) || !(u < tt)) return R_NegInf;
    return wald_log(u, vs, Bs, As, 'f') + log1m_exp(-lam * (tt - u));
  }));
}

}  // namespace

double exwald_log_density(double t, double v, double B, double A, double s, double tau) {
  return exwald_log(t, v, B, A, s, tau, 'f');
}
double exwald_log_cdf(double t, double v, double B, double A, double s, double tau) {
  return exwald_log(t, v, B, A, s, tau, 'F');
}
double exwald_log_surv(double t, double v, double B, double A, double s, double tau) {
  return exwald_log(t, v, B, A, s, tau, 'S');
}

double drdmswtn_exwald_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[ecol::v])) return 0.0;
  const double tt = t - par[ecol::t0];
  if (!(tt > 0.0)) return 0.0;
  const double l = exwald_log_density(tt, par[ecol::v], par[ecol::B], par[ecol::A],
                                      par[ecol::s], par[ecol::tau]);
  return (R_FINITE(l) && l > kLogUnderflow) ? std::exp(l) : 0.0;
}

double prdmswtn_exwald_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[ecol::v])) return 0.0;
  const double tt = t - par[ecol::t0];
  if (!(tt > 0.0)) return 0.0;
  const double l = exwald_log_cdf(tt, par[ecol::v], par[ecol::B], par[ecol::A],
                                  par[ecol::s], par[ecol::tau]);
  return R_FINITE(l) ? std::exp(l) : 0.0;
}

void drdmswtn_exwald_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_   = cols[ecol::v];
  const double* B_   = cols[ecol::B];
  const double* A_   = cols[ecol::A];
  const double* t0_  = cols[ecol::t0];
  const double* s_   = cols[ecol::s];
  const double* tau_ = cols[ecol::tau];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double l = exwald_log_density(tt, v_[i], B_[i], A_[i], s_[i], tau_[i]);
    out[i] = raw_log_value(l, min_ll, floor_raw);
  }
}

void prdmswtn_exwald_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_   = cols[ecol::v];
  const double* B_   = cols[ecol::B];
  const double* A_   = cols[ecol::A];
  const double* t0_  = cols[ecol::t0];
  const double* s_   = cols[ecol::s];
  const double* tau_ = cols[ecol::tau];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = 0.0; continue; }
    const double ls = exwald_log_surv(tt, v_[i], B_[i], A_[i], s_[i], tau_[i]);
    if (ISNAN(ls)) { out[i] = 0.0; continue; }
    out[i] = (ls == R_NegInf) ? raw_log_zero(min_ll, floor_raw) : std::fmin(ls, 0.0);
  }
}

void rdmswtn_exwald_logS_at_t(double t, const double* const* cols,
                              int /*n_rows_total*/, int n_lR, int /*n_par*/,
                              const int* trunc_mask, int n_unique_trials,
                              const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* v_   = cols[ecol::v];
  const double* B_   = cols[ecol::B];
  const double* A_   = cols[ecol::A];
  const double* t0_  = cols[ecol::t0];
  const double* s_   = cols[ecol::s];
  const double* tau_ = cols[ecol::tau];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int k = 0; k < n_lR; ++k) {
      const int r = start + k;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;
      const double ls = exwald_log_surv(tt, v_[r], B_[r], A_[r], s_[r], tau_[r]);
      if (ISNAN(ls) || ls == R_NegInf) { bad = true; break; }
      logS += std::fmin(ls, 0.0);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// R-level dfun/pfun.  t is RT (t0 is subtracted here).  lower_tail = FALSE
// gives the survivor.
// [[Rcpp::export]]
Rcpp::NumericVector dRDMSWTN_exwald_cpp(Rcpp::NumericVector t, Rcpp::NumericVector v,
                                        Rcpp::NumericVector B, Rcpp::NumericVector A,
                                        Rcpp::NumericVector t0, Rcpp::NumericVector s,
                                        Rcpp::NumericVector tau, bool log_out = false) {
  const R_xlen_t n = t.size();
  if (v.size() != n || B.size() != n || A.size() != n || t0.size() != n || s.size() != n ||
      tau.size() != n)
    Rcpp::stop("dRDMSWTN_exwald_cpp: arguments must have equal lengths.");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double l = exwald_log_density(t[i] - t0[i], v[i], B[i], A[i], s[i], tau[i]);
    out[i] = log_out ? l : std::exp(l);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector pRDMSWTN_exwald_cpp(Rcpp::NumericVector t, Rcpp::NumericVector v,
                                        Rcpp::NumericVector B, Rcpp::NumericVector A,
                                        Rcpp::NumericVector t0, Rcpp::NumericVector s,
                                        Rcpp::NumericVector tau, bool log_out = false,
                                        bool lower_tail = true) {
  const R_xlen_t n = t.size();
  if (v.size() != n || B.size() != n || A.size() != n || t0.size() != n || s.size() != n ||
      tau.size() != n)
    Rcpp::stop("pRDMSWTN_exwald_cpp: arguments must have equal lengths.");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double tt = t[i] - t0[i];
    const double l = lower_tail ? exwald_log_cdf(tt, v[i], B[i], A[i], s[i], tau[i])
                                : exwald_log_surv(tt, v[i], B[i], A[i], s[i], tau[i]);
    out[i] = log_out ? l : std::exp(l);
  }
  return out;
}
