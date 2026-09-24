#ifndef EMC2_MODEL_FRQ_FADE_H
#define EMC2_MODEL_FRQ_FADE_H

#include "model_FRQ.h"
// composite_functions.h consumes emc2_isfinite from utility_functions.h.
#include "utility_functions.h"
#include "composite_functions.h"

#include <array>
#include <cstdint>
#include <unordered_map>

struct FrqFadePars {
  double alpha = 0.0;
  double beta = 0.0;
  double lambda = 0.0;
  double log_lambda = 0.0;
  double kappa = 0.0;
  double cv2 = 0.0;
  double lbeta_val = 0.0;
  double alpha_minus_1 = 0.0;
  double beta_minus_1 = 0.0;
  FrqH hh;
  bool ok = false;
  // The infinite-time endpoint costs two incomplete-beta evaluations and is
  // needed only by u == Inf queries and the mapped summary, never by a finite
  // RT.  It is filled on first use by frqfade_ensure_endpoint(), so a row
  // that misses the memo (a trend or covariate on lambda or kappa) pays only
  // for what its own evaluator reads.
  mutable bool endpoint_ready = false;
  mutable double q_inf = 0.0;
  mutable double omq_inf = 1.0;
  mutable double log_h = R_NegInf;
  mutable double log_om_h = 0.0;
};

struct FrqFadeKey {
  std::array<std::uint64_t, 2> bits;
  bool operator==(const FrqFadeKey& other) const { return bits == other.bits; }
};

struct FrqFadeKeyHash {
  std::size_t operator()(const FrqFadeKey& key) const {
    std::uint64_t x = key.bits[0] + 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x ^= key.bits[1] + 0x94d049bb133111ebULL + (x << 6) + (x >> 2);
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<std::size_t>(x ^ (x >> 31));
  }
};

inline std::uint64_t frqfade_bits(double x) {
  std::uint64_t bits;
  std::memcpy(&bits, &x, sizeof(bits));
  return bits;
}

struct FrqFadeCache {
  std::unordered_map<FrqFadeKey, double, FrqFadeKeyHash> lbeta;

  FrqFadeCache() { lbeta.reserve(64); }

  double get_lbeta(double alpha, double beta) {
    const FrqFadeKey key{{frqfade_bits(alpha), frqfade_bits(beta)}};
    const auto it = lbeta.find(key);
    if (it != lbeta.end()) return it->second;
    const double value = R::lbeta(alpha, beta);
    lbeta.emplace(key, value);
    return value;
  }

  void clear() { lbeta.clear(); }
};

inline double frqfade_log1p_product(double a, double b, double c) {
  if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0)) return 0.0;
  const double lx = std::log(a) + std::log(b) + std::log(c);
  return log1p_exp(lx);
}

inline double frqfade_log_survival(double lambda, double cv2, double A) {
  if (!(A > 0.0)) return 0.0;
  if (cv2 == 0.0) return -lambda * A;
  const double lx = std::log(lambda) + std::log(cv2) + std::log(A);
  const double log_d = log1p_exp(lx);
  if (!R_FINITE(log_d)) return R_NegInf;
  // Keep the exponential-rate limit when cv2*A*lambda is tiny.
  if (lx < -18.0) {
    const double x = std::exp(std::log(lambda) + std::log(A));
    const double z = std::exp(lx);
    return (z > 0.0) ? -x * (std::log1p(z) / z) : -x;
  }
  return -log_d / cv2;
}

struct FrqFadeState {
  double q = 0.0;
  double omq = 1.0;
  double log_qprime = R_NegInf;
};

inline void frqfade_ensure_endpoint(const FrqFadePars& s);

inline FrqFadeState frqfade_state(double u, const FrqFadePars& s) {
  FrqFadeState z;
  if (!R_FINITE(u)) {
    frqfade_ensure_endpoint(s);
    z.q = s.q_inf;
    z.omq = s.omq_inf;
    return z;
  }
  const double k_u = s.kappa * u;
  const double A = (s.kappa == 0.0)
    ? u
    : -std::expm1(-k_u) / s.kappa;
  const double log_surv = frqfade_log_survival(s.lambda, s.cv2, A);
  z.q = -std::expm1(log_surv);
  z.omq = std::exp(log_surv);
  const double log_d = (s.cv2 == 0.0)
    ? 0.0
    : frqfade_log1p_product(s.lambda, s.cv2, A);
  z.log_qprime = s.log_lambda - k_u + log_surv - log_d;
  return z;
}

inline FrqLogBeta frqfade_log_beta_pair(const FrqFadeState& z,
                                        const FrqFadePars& s) {
  FrqLogBeta out;
  out.lz = frq_log_ibeta(z.q, z.omq, s.alpha, s.beta, s.lbeta_val);
  out.lomz = frq_log_ibeta(z.omq, z.q, s.beta, s.alpha, s.lbeta_val);
  return out;
}

// Fill the infinite-time endpoint (q_inf, h) on first use.  Only valid rows
// reach it; an invalid row keeps the NaN-free defaults and is rejected by the
// callers' `ok` checks first.
inline void frqfade_ensure_endpoint(const FrqFadePars& s) {
  if (s.endpoint_ready || !s.ok) return;
  s.endpoint_ready = true;
  if (s.kappa == 0.0) {
    s.q_inf = 1.0;
    s.omq_inf = 0.0;
    s.log_h = 0.0;
    s.log_om_h = R_NegInf;
    return;
  }
  const double A_inf = 1.0 / s.kappa;
  const double log_omq = frqfade_log_survival(s.lambda, s.cv2, A_inf);
  s.omq_inf = std::exp(log_omq);
  s.q_inf = -std::expm1(log_omq);
  const double lz = frq_log_ibeta(s.q_inf, s.omq_inf, s.alpha, s.beta,
                                  s.lbeta_val);
  const double lomz = frq_log_ibeta(s.omq_inf, s.q_inf, s.beta, s.alpha,
                                    s.lbeta_val);
  if (s.hh.active) {
    s.log_h = frq_h_log(lz - lomz, s.hh);
    s.log_om_h = frq_h_log(lomz - lz, s.hh);
  } else {
    s.log_h = lz;
    s.log_om_h = lomz;
  }
}

inline FrqFadePars frqfade_derive(double alpha, double beta, double lambda,
                                 double kappa, double delta, double cv_u,
                                 double lbeta_val) {
  FrqFadePars s;
  if (ISNAN(alpha) || ISNAN(beta) || ISNAN(lambda) || ISNAN(kappa) ||
      ISNAN(delta) || ISNAN(cv_u)) return s;
  if (!(alpha >= 1.0) || !R_FINITE(alpha) ||
      !(beta >= 1.0) || !R_FINITE(beta) ||
      !(lambda > 0.0) || !R_FINITE(lambda) ||
      !(kappa >= 0.0) || !R_FINITE(kappa) ||
      !(cv_u >= 0.0) || !R_FINITE(cv_u)) return s;
  s.hh = frq_h_make(delta);
  if (!s.hh.ok) return s;
  const double cv2 = cv_u * cv_u;
  if (!R_FINITE(cv2)) return s;
  s.alpha = alpha;
  s.beta = beta;
  s.lambda = lambda;
  s.log_lambda = std::log(lambda);
  s.kappa = kappa;
  s.cv2 = cv2;
  s.lbeta_val = lbeta_val;
  s.alpha_minus_1 = alpha - 1.0;
  s.beta_minus_1 = beta - 1.0;
  s.ok = true;
  return s;
}

struct FrqFadeMemo {
  FrqFadeCache* shared = nullptr;
  double alpha = R_NaN, beta = R_NaN, lambda = R_NaN;
  double kappa = R_NaN, delta = R_NaN, cv_u = R_NaN;
  double lbeta_val = R_NaN;
  FrqFadePars value;

  explicit FrqFadeMemo(FrqFadeCache* cache = nullptr) : shared(cache) {}

  const FrqFadePars& get(double a, double b, double l, double k,
                         double d, double cv) {
    if (a == alpha && b == beta && l == lambda && k == kappa &&
        d == delta && cv == cv_u) return value;
    if (!(a == alpha && b == beta)) {
      alpha = a;
      beta = b;
      lbeta_val = (shared != nullptr && a >= 1.0 && b >= 1.0 &&
                   R_FINITE(a) && R_FINITE(b))
        ? shared->get_lbeta(a, b) : R::lbeta(a, b);
    }
    lambda = l;
    kappa = k;
    delta = d;
    cv_u = cv;
    value = frqfade_derive(a, b, l, k, d, cv, lbeta_val);
    return value;
  }
};

inline double frqfade_log_pdf_dt(double u, const FrqFadePars& s) {
  if (!s.ok || !(u > 0.0) || !R_FINITE(u)) return R_NegInf;
  const FrqFadeState z = frqfade_state(u, s);
  double out = z.log_qprime - s.lbeta_val;
  out += frq_xlogy(s.alpha_minus_1, z.q);
  out += frq_xlogy(s.beta_minus_1, z.omq);
  if (s.hh.active) {
    const FrqLogBeta b = frqfade_log_beta_pair(z, s);
    out += frq_h_logderiv(b.lz, b.lomz, s.hh);
  }
  return ISNAN(out) ? R_NegInf : out;
}

inline double frqfade_log_cdf_dt(double u, const FrqFadePars& s) {
  if (!s.ok || ISNAN(u) || !(u > 0.0)) return R_NegInf;
  if (!R_FINITE(u)) {
    frqfade_ensure_endpoint(s);
    return s.log_h;
  }
  const FrqFadeState z = frqfade_state(u, s);
  double out;
  if (s.hh.active) {
    const FrqLogBeta b = frqfade_log_beta_pair(z, s);
    out = frq_h_log(b.lz - b.lomz, s.hh);
  } else {
    out = frq_log_ibeta(z.q, z.omq, s.alpha, s.beta, s.lbeta_val);
  }
  return ISNAN(out) ? R_NegInf : out;
}

inline double frqfade_log_surv_dt(double u, const FrqFadePars& s) {
  if (!s.ok) return R_NaN;
  if (ISNAN(u)) return R_NaN;
  if (!(u > 0.0)) return 0.0;
  if (!R_FINITE(u)) {
    frqfade_ensure_endpoint(s);
    return s.log_om_h;
  }
  const FrqFadeState z = frqfade_state(u, s);
  double out;
  if (s.hh.active) {
    const FrqLogBeta b = frqfade_log_beta_pair(z, s);
    out = frq_h_log(b.lomz - b.lz, s.hh);
  } else {
    out = frq_log_ibeta(z.omq, z.q, s.beta, s.alpha, s.lbeta_val);
  }
  return ISNAN(out) ? R_NaN : out;
}

inline double frqfade_pdf_natural_dt(double u, const FrqFadePars& s) {
  const double lp = frqfade_log_pdf_dt(u, s);
  return lp > R_NegInf ? std::exp(lp) : 0.0;
}

inline double frqfade_cdf_natural_dt(double u, const FrqFadePars& s) {
  const double lp = frqfade_log_cdf_dt(u, s);
  if (!(lp > R_NegInf)) return 0.0;
  const double out = std::exp(lp);
  return out > 1.0 ? 1.0 : out;
}

inline double frqfade_internal_clock(double p, const FrqFadePars& s) {
  if (!(p >= 0.0) || !(p < 1.0)) return R_PosInf;
  const double y = -std::log1p(-p);
  if (s.cv2 == 0.0) return y / s.lambda;
  const double z = s.cv2 * y;
  if (z > 709.782712893384) return R_PosInf;
  if (z == 0.0) return y / s.lambda;
  return (y / s.lambda) * (std::expm1(z) / z);
}

inline double frqfade_decision_time(double p, const FrqFadePars& s) {
  if (s.kappa == 0.0) {
    // Unlimited opportunity is proper: an exact p = 1 is a zero-probability
    // endpoint that finite RNGs can still return, so keep it finite rather
    // than creating an omission (as rfrq_cpp does at its p = 1 boundary).
    return frqfade_internal_clock(p >= 1.0 ? std::nextafter(1.0, 0.0) : p, s);
  }
  const double A = frqfade_internal_clock(p, s);
  const double kA = s.kappa * A;
  if (!(kA < 1.0)) return R_PosInf;
  return -std::log1p(-kA) / s.kappa;
}

double dfrqfade_scalar(double t, const double* par, void* ctx_);
double pfrqfade_scalar(double t, const double* par, void* ctx_);
void dfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void pfrqfade_raw(const double* rt, const double* const* cols, int n_rows,
                  const int* mask, const int* isok,
                  double* out, double min_ll, void* ctx_);
void frqfade_logS_at_t(double t, const double* const* cols,
                       int n_rows_total, int n_lR, int n_par,
                       const int* trunc_mask, int n_unique_trials,
                       const int* isok_all, void* ctx_, double* logS_out);

#endif  // EMC2_MODEL_FRQ_FADE_H
