#ifndef EMC2_MODEL_PCOUNTER_H
#define EMC2_MODEL_PCOUNTER_H

// ============================================================================
// PCOUNTER -- gamma-distributed trialwise input with pure-birth self-excitation
// and a geometric threshold excess.
//
// CLOSED FORM.  Input arrivals follow a gamma(shape, rate) process whose rate
// is itself gamma-distributed (nu, sv^2); the first K = 2 + floor(k + 0.5)
// threshold excesses trigger a response.  All probabilities are computed on
// the log scale from exact finite sums (Stirling numbers of the second kind for
// the self-excitation term), so no quadrature is needed.
//
// SPLIT.  Only the tiny scalar predicates/helpers that every entry point
// shares stay inline here.  The heavier kernels (the rolling Stirling row,
// pcounter_log_* closed forms, and the race adapters) have their bodies in
// model_PCOUNTER.cpp, which is also the sole owner of the [[Rcpp::export]]
// entry points dpcounter/ppcounter and their recycle helpers.  This header is
// therefore safe to include from any number of translation units.
// ============================================================================

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// Variability/excess parameters (sv, gamma, omega) at or below this threshold
// activate the degenerate (zero-variability / non-self-exciting / geometric
// fixed-threshold) branches.
static constexpr double PC_EPS = 1e-8;

// ---------------------------------------------------------------------------
// Inline scalar predicates/helpers.
// ---------------------------------------------------------------------------

// The canonical event threshold K = 2 + floor(k + 0.5) shared by every C++
// entry point.  There is deliberately no model-level upper bound on K;
// this only rejects values that cannot be represented by the integer loop
// indices used by the exact finite-K formulas.
inline bool pcounter_k_supported(double k) {
  if (!R_FINITE(k) || k < 0.0) return false;
  const double z = std::floor(k + 0.5);
  // Tail formulas request rows through K+65; keep every int loop index
  // representable without imposing a fitted-model threshold cap.
  return z <= static_cast<double>(std::numeric_limits<int>::max() - 67);
}

inline bool pcounter_needs_stirling(double sv, double gamma) {
  return R_FINITE(sv) && R_FINITE(gamma) &&
         sv >= PC_EPS && gamma >= PC_EPS;
}

inline double pcounter_logsumexp(const std::vector<double>& x) {
  if (x.empty()) return R_NegInf;
  double m = R_NegInf;
  for (double z : x) if (z > m) m = z;
  if (!R_FINITE(m)) return R_NegInf;
  double s = 0.0;
  for (double z : x) if (R_FINITE(z)) s += std::exp(z - m);
  return m + std::log(s);
}

inline double pcounter_logdiffexp(double a, double b) {
  if (!R_FINITE(a)) return R_NegInf;
  if (!R_FINITE(b)) return a;
  if (b >= a) return R_NegInf;
  return a + std::log1p(-std::exp(b - a));
}

inline int pcounter_k_int(double k) {
  if (!pcounter_k_supported(k)) return 0;
  const double z = std::floor(k + 0.5);
  return static_cast<int>(z) + 2;
}

inline double pcounter_log_rising(double a, int n) {
  double out = 0.0;
  for (int j = 0; j < n; ++j) out += std::log(a + static_cast<double>(j));
  return out;
}

// ---------------------------------------------------------------------------
// Non-inline kernels.  Bodies live in model_PCOUNTER.cpp.
// ---------------------------------------------------------------------------

// A single rolling row of log Stirling numbers of the second kind, stored in
// O(n) rather than the former O(n^2) table.  Only forward-declared here; the
// heavy body belongs to model_PCOUNTER.cpp.
class PcounterStirlingRows;

void pcounter_log_lm(double a, double nu, double sv,
                     double& logL, double& logM);
double pcounter_log_h(int n, double t, double nu, double sv, double gamma,
                      PcounterStirlingRows& rows);
void pcounter_log_pi_phi(int n, double t, double nu, double sv, double gamma,
                         bool gamma_zero,
                         PcounterStirlingRows& rows,
                         double& logpi, double& logphi);
double pcounter_log_tail(int start, double t, double nu, double sv, double gamma,
                         double logr, bool gamma_zero,
                         PcounterStirlingRows& rows, bool phi);
double pcounter_log_fixed_cdf(int start, double t, double nu, double sv,
                              double gamma, bool gamma_zero,
                              PcounterStirlingRows& rows);
double pcounter_log_geom_cdf(int start, double t, double nu, double sv,
                             double gamma, double omega, bool gamma_zero,
                             PcounterStirlingRows& rows);
void pcounter_log_eval(double t, double nu, double sv, double gamma,
                       double k, double omega,
                       double& logf, double& logS, double& logF);

// ---------------------------------------------------------------------------
// Race-model adapters (scalar PDF/CDF, raw batch, and batch log-survivor).
// ---------------------------------------------------------------------------
double dpcounter_scalar(double t, const double* par, void* ctx_);
double ppcounter_scalar(double t, const double* par, void* ctx_);
void dpcounter_raw(const double* rt, const double* const* cols, int n_rows,
                   const int* mask, const int* isok, double* out,
                   double min_ll, void* ctx_);
void ppcounter_raw(const double* rt, const double* const* cols, int n_rows,
                   const int* mask, const int* isok, double* out,
                   double min_ll, void* ctx_);
void pcounter_logS_at_t(double t, const double* const* cols,
                        int n_rows_total, int n_lR, int n_par,
                        const int* trunc_mask, int n_unique_trials,
                        const int* isok_all, void* ctx_, double* logS_out);

// ---------------------------------------------------------------------------
// R-callable entry points and their recycle helpers.
// ---------------------------------------------------------------------------
void pcounter_check_len(const NumericVector& x, int n, const char* nm);
double pcounter_pick(const NumericVector& x, int i);
NumericVector dpcounter(NumericVector t, NumericVector nu, NumericVector sv,
                        NumericVector gamma, NumericVector k,
                        NumericVector omega, NumericVector t0,
                        bool log_out);
NumericVector ppcounter(NumericVector t, NumericVector nu, NumericVector sv,
                        NumericVector gamma, NumericVector k,
                        NumericVector omega, NumericVector t0,
                        bool lower_tail, bool log_out);

#endif // EMC2_MODEL_PCOUNTER_H
