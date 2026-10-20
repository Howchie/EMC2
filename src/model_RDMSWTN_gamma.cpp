#include <Rcpp.h>
#include "model_RDMSWTN_gamma.h"
#include "rdmswtn_gamma_A.h"
#include "wald_functions.h"
#include "col_registry.h"

// ============================================================
// RDMSWTN with Gamma drift (c_name RDMSWTN_GDRIFT)
// Column layout: v=0 (mean drift), B=1, A=2, t0=3, s=4, sv=5 (drift SD).
// V ~ Gamma(k, rate tau) with k = (v/sv)^2, tau = v/sv^2; distance to
// threshold uniform on [B, B + A] (A = 0 is a point start); no clocks.
// sv = 0 is the plain Wald.
// ============================================================

namespace {

namespace gcol = emc2col::rdmswtn_gamma;

// Drift SD below this fraction of the mean is the Wald limit: the gamma
// correction to the Wald law is O(cv^2).
constexpr double kWaldCv = 1e-6;
// exp() of a log below this is 0 in double, so a density that is returned on the
// natural scale may be bounded there.
constexpr double kLogUnderflow = -746.0;

inline bool wald_limit(double v, double sv) {
  return !emc2_isfinite(sv) || !(sv > kWaldCv * v);
}

// log f, log F or log S of the decision time tt (want: one of the Want bits).
// A point-start log density below lf_floor may come back as an upper bound
// still below it.
inline double gdrift_log(double tt, double v, double B, double A, double s, double sv,
                         int want, double lf_floor = R_NegInf) {
  using namespace rdmswtn_gamma;
  if (!(tt > 0.0)) return (want == WANT_SURV) ? 0.0 : R_NegInf;
  if (!(v > 0.0) || !(s > 0.0) || !(B >= 0.0) || !(A >= 0.0) || !(B + A > 0.0) ||
      !(sv >= 0.0)) return R_NaN;
  if (!R_FINITE(tt)) return (want == WANT_CDF) ? 0.0 : R_NegInf;
  const double inv_s = 1.0 / s;
  if (wald_limit(v, sv)) {
    const double mu = v * inv_s, bb = (B + A) * inv_s, aa = A * inv_s;
    if (want == WANT_PDF) {
      double pdf;
      if (dwald_k0_natural(tt, bb, mu, aa, pdf)) return (pdf > 0.0) ? std::log(pdf) : R_NegInf;
      return dwald_k0_log(tt, bb, mu, aa);
    }
    const double cdf = std::max(0.0, std::min(1.0, pwald_k0(tt, bb, mu, aa)));
    if (want == WANT_CDF) return std::log(cdf);
    if (cdf <= 0.0) return 0.0;
    if (cdf >= 1.0 - EMC2_CDF_SAT_MARGIN) return wald_k0_log_surv(tt, bb, mu, aa);
    return std::log1p(-cdf);
  }
  const double k = (v / sv) * (v / sv), tau = v / (sv * sv);
  if (A > 0.0) {
    const AOut o = gamma_wald_A(tt, B, A, s, k, tau, want, lf_floor);
    return (want == WANT_PDF) ? o.logf : (want == WANT_CDF ? o.logF : o.logS);
  }
  const Out o = gamma_wald(tt, B, s, k, tau, want, 0, lf_floor);
  return (want == WANT_PDF) ? o.logf : (want == WANT_CDF ? o.logF : o.logS);
}

}  // namespace

double gdrift_log_density(double t, double v, double B, double A, double s, double sv,
                          double lf_floor) {
  return gdrift_log(t, v, B, A, s, sv, rdmswtn_gamma::WANT_PDF, lf_floor);
}
double gdrift_log_cdf(double t, double v, double B, double A, double s, double sv) {
  return gdrift_log(t, v, B, A, s, sv, rdmswtn_gamma::WANT_CDF);
}
double gdrift_log_surv(double t, double v, double B, double A, double s, double sv) {
  return gdrift_log(t, v, B, A, s, sv, rdmswtn_gamma::WANT_SURV);
}

double drdmswtn_gdrift_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[gcol::v])) return 0.0;
  const double tt = t - par[gcol::t0];
  if (!(tt > 0.0)) return 0.0;
  const double l = gdrift_log_density(tt, par[gcol::v], par[gcol::B], par[gcol::A],
                                     par[gcol::s], par[gcol::sv], kLogUnderflow);
  return R_FINITE(l) ? std::exp(l) : 0.0;
}

double prdmswtn_gdrift_scalar(double t, const double* par, void* /*ctx_*/) {
  if (R_IsNA(par[gcol::v])) return 0.0;
  const double tt = t - par[gcol::t0];
  if (!(tt > 0.0)) return 0.0;
  const double l = gdrift_log_cdf(tt, par[gcol::v], par[gcol::B], par[gcol::A],
                                 par[gcol::s], par[gcol::sv]);
  return R_FINITE(l) ? std::exp(l) : 0.0;
}

void drdmswtn_gdrift_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[gcol::v];
  const double* B_  = cols[gcol::B];
  const double* A_  = cols[gcol::A];
  const double* t0_ = cols[gcol::t0];
  const double* s_  = cols[gcol::s];
  const double* sv_ = cols[gcol::sv];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    // a floored density needs no digits below the floor
    const double l = gdrift_log_density(tt, v_[i], B_[i], A_[i], s_[i], sv_[i],
                                       floor_raw ? min_ll : R_NegInf);
    out[i] = raw_log_value(l, min_ll, floor_raw);
  }
}

void prdmswtn_gdrift_raw(const double* rt, const double* const* cols, int n_rows,
                         const int* mask, const int* isok,
                         double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* v_  = cols[gcol::v];
  const double* B_  = cols[gcol::B];
  const double* A_  = cols[gcol::A];
  const double* t0_ = cols[gcol::t0];
  const double* s_  = cols[gcol::s];
  const double* sv_ = cols[gcol::sv];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (R_IsNA(v_[i]) || !isok[i]) { out[i] = 0.0; continue; }
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0)) { out[i] = 0.0; continue; }
    const double ls = gdrift_log_surv(tt, v_[i], B_[i], A_[i], s_[i], sv_[i]);
    if (ISNAN(ls)) { out[i] = 0.0; continue; }
    out[i] = (ls == R_NegInf) ? raw_log_zero(min_ll, floor_raw) : std::fmin(ls, 0.0);
  }
}

void rdmswtn_gdrift_logS_at_t(double t, const double* const* cols,
                              int /*n_rows_total*/, int n_lR, int /*n_par*/,
                              const int* trunc_mask, int n_unique_trials,
                              const int* isok_all, void* /*ctx_*/, double* logS_out) {
  const double* v_  = cols[gcol::v];
  const double* B_  = cols[gcol::B];
  const double* A_  = cols[gcol::A];
  const double* t0_ = cols[gcol::t0];
  const double* s_  = cols[gcol::s];
  const double* sv_ = cols[gcol::sv];
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
      const double ls = gdrift_log_surv(tt, v_[r], B_[r], A_[r], s_[r], sv_[r]);
      if (ISNAN(ls) || ls == R_NegInf) { bad = true; break; }
      logS += std::fmin(ls, 0.0);
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// R-level dfun/pfun.  Columns as in the kernel contract; t is RT (t0 is
// subtracted here).  lower_tail = FALSE gives the survivor.
// [[Rcpp::export]]
Rcpp::NumericVector dRDMSWTN_gdrift_cpp(Rcpp::NumericVector t, Rcpp::NumericVector v,
                                        Rcpp::NumericVector B, Rcpp::NumericVector A,
                                        Rcpp::NumericVector t0,
                                        Rcpp::NumericVector s, Rcpp::NumericVector sv,
                                        bool log_out = false) {
  const R_xlen_t n = t.size();
  if (v.size() != n || B.size() != n || A.size() != n || t0.size() != n || s.size() != n || sv.size() != n)
    Rcpp::stop("dRDMSWTN_gdrift_cpp: arguments must have equal lengths.");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double l = gdrift_log_density(t[i] - t0[i], v[i], B[i], A[i], s[i], sv[i]);
    out[i] = log_out ? l : std::exp(l);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector pRDMSWTN_gdrift_cpp(Rcpp::NumericVector t, Rcpp::NumericVector v,
                                        Rcpp::NumericVector B, Rcpp::NumericVector A,
                                        Rcpp::NumericVector t0,
                                        Rcpp::NumericVector s, Rcpp::NumericVector sv,
                                        bool log_out = false, bool lower_tail = true) {
  const R_xlen_t n = t.size();
  if (v.size() != n || B.size() != n || A.size() != n || t0.size() != n || s.size() != n || sv.size() != n)
    Rcpp::stop("pRDMSWTN_gdrift_cpp: arguments must have equal lengths.");
  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const double tt = t[i] - t0[i];
    const double l = lower_tail ? gdrift_log_cdf(tt, v[i], B[i], A[i], s[i], sv[i])
                                : gdrift_log_surv(tt, v[i], B[i], A[i], s[i], sv[i]);
    out[i] = log_out ? l : std::exp(l);
  }
  return out;
}
