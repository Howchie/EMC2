#ifndef model_GOM_core_h
#define model_GOM_core_h

// Shared, lightweight Gompertz-to-OU mapping.  This header is intentionally
// independent of utils.h so R-facing vector helpers do not pull in the full
// likelihood translation unit and duplicate its non-inline symbols.

#include "col_registry.h"
#include "fpe_race.h"

struct GompColIdx {
  int alpha, beta, K, B, A, t0;
};

inline GompColIdx gomp_col_idx() {
  return {emc2col::gompertz::alpha, emc2col::gompertz::beta,
          emc2col::gompertz::K, emc2col::gompertz::B,
          emc2col::gompertz::A, emc2col::gompertz::t0};
}

inline fperace::BndSpec gomp_bnd_row(int kind, const double* const* cols,
                                     int i) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = cols[emc2col::gompertz::Binf][i];
  bs.tau = cols[emc2col::gompertz::tau][i];
  bs.pw = (kind == fpe::FPE_BND_WEIBULL)
            ? cols[emc2col::gompertz::pw][i] : 0.0;
  return bs;
}

inline fperace::BndSpec gomp_bnd_par(int kind, const double* par) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = par[emc2col::gompertz::Binf];
  bs.tau = par[emc2col::gompertz::tau];
  bs.pw = (kind == fpe::FPE_BND_WEIBULL)
            ? par[emc2col::gompertz::pw] : 0.0;
  return bs;
}

// Build the log-space OU key from the physical Gompertz parameters.  The
// process starts Uniform(1, 1 + A) in X-space and has initial threshold
// 1 + B + A; A, B and all collapse parameters remain physical in the key.
inline bool gomp_key(double alpha, double beta, double K, double B, double A,
                     const fperace::BndSpec& bs, fperace::Key& out) {
  if (!(A >= 0.0) || !std::isfinite(A)) return false;
  const double AA = A;
  const double b0 = 1.0 + B + AA;
  if (!(alpha > 0.0) || !std::isfinite(alpha) ||
      !(beta > 0.0) || !std::isfinite(beta) ||
      !(K > 0.0) || !std::isfinite(K) ||
      !(B >= 0.0) || !std::isfinite(B) ||
      !(b0 > 1.0) || !std::isfinite(b0)) return false;

  out.model_kind = 1;
  out.log_state = true;
  out.sigma = beta;
  out.k = alpha;
  out.v = alpha * std::log(K) - 0.5 * beta * beta;
  out.b = b0;
  out.A = AA;
  out.zlo = 0.0;
  out.zhi = std::log1p(AA);
  out.bkind = bs.kind;
  if (bs.kind == fpe::FPE_BND_FIXED) {
    out.binf = 0.0;
    out.tau = 0.0;
    out.pw = 0.0;
  } else {
    if (!(bs.Binf > 0.0) || !std::isfinite(bs.Binf) ||
        !(bs.tau > 0.0) || !std::isfinite(bs.tau)) return false;
    out.binf = bs.Binf;
    out.tau = bs.tau;
    out.pw = (bs.kind == fpe::FPE_BND_WEIBULL) ? bs.pw : 0.0;
    if (bs.kind == fpe::FPE_BND_WEIBULL &&
        (!(bs.pw > 0.0) || !std::isfinite(bs.pw))) return false;
  }
  return out.finite() && out.b > 1.0 && out.zhi >= out.zlo;
}

#endif // model_GOM_core_h
