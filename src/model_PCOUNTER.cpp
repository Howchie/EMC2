// Non-inline bodies and R-callable entry points for PCOUNTER.
//
// These live in their own translation unit rather than in model_PCOUNTER.h
// because the header is small (inline predicates + declarations) and because a
// header carrying [[Rcpp::export]] definitions can only be seen by one TU.  The
// race-path kernels and adapters keep their exact signatures, branch order,
// formulas and raw log-survivor flooring behaviour; only their home moves here.

#include "model_PCOUNTER.h"

#include "col_registry.h"
#include "race_contract.h"

namespace {

inline double pcorr_log_choose(int n, int k) {
  if (k < 0 || k > n) return R_NegInf;
  k = std::min(k, n - k);
  double out = 0.0;
  for (int i = 1; i <= k; ++i)
    out += std::log(static_cast<double>(n - k + i)) - std::log(static_cast<double>(i));
  return out;
}

inline double pcorr_log_rising(double x, int n) {
  if (n == 0) return 0.0;
  if (!(x > 0.0) || !R_FINITE(x)) return R_NegInf;
  double out = 0.0;
  for (int i = 0; i < n; ++i) out += std::log(x + static_cast<double>(i));
  return out;
}

inline double pcorr_logsumexp_n(const double* x, int n) {
  if (n <= 0) return R_NegInf;
  double m = R_NegInf;
  for (int i = 0; i < n; ++i) if (x[i] > m) m = x[i];
  if (!R_FINITE(m)) return R_NegInf;
  double s = 0.0;
  for (int i = 0; i < n; ++i) if (R_FINITE(x[i])) s += std::exp(x[i] - m);
  return m + std::log(s);
}

// log E[lambda_1^r lambda_2^q exp(-s1 lambda_1-s2 lambda_2)] for the
// common-CV shared-shape construction.  All summands are non-negative.
double pcorr_log_moment(int r, int q, double s1, double s2,
                        double nu1, double nu2, double cv, double rho) {
  if (cv <= PC_EPS) {
    return r * std::log(nu1) + q * std::log(nu2) - s1 * nu1 - s2 * nu2;
  }
  const double a = 1.0 / (cv * cv), a0 = rho * a;
  const double b1 = a / nu1, b2 = a / nu2;
  const double l0 = std::log1p(s1 / b1 + s2 / b2);
  const double l1 = std::log1p(s1 / b1), l2 = std::log1p(s2 / b2);
  const int n_terms = (r + 1) * (q + 1);
  double terms_fixed[64];
  std::vector<double> terms;
  if (n_terms <= 64)
    std::fill(terms_fixed, terms_fixed + n_terms, R_NegInf);
  if (n_terms > 64) terms.reserve(static_cast<size_t>(n_terms));
  for (int u = 0; u <= r; ++u) for (int v = 0; v <= q; ++v) {
    const int shared_power = u + v;
    // (0)_m is zero for m>0; spelling this out avoids lgamma(0).
    if ((a0 == 0.0 && shared_power > 0) ||
        (a0 == a && r - u > 0) || (a0 == a && q - v > 0)) continue;
    double z = pcorr_log_choose(r, u) + pcorr_log_choose(q, v) -
      r * std::log(b1) - q * std::log(b2);
    if (shared_power) z += pcorr_log_rising(a0, shared_power);
    z -= (a0 + shared_power) * l0;
    const int p1 = r - u, p2 = q - v;
    if (p1) z += pcorr_log_rising(a - a0, p1);
    if (p2) z += pcorr_log_rising(a - a0, p2);
    z -= (a - a0 + p1) * l1 + (a - a0 + p2) * l2;
    if (n_terms <= 64) terms_fixed[static_cast<size_t>(u * (q + 1) + v)] = z;
    else terms.push_back(z);
  }
  return n_terms <= 64 ? pcorr_logsumexp_n(terms_fixed, n_terms) :
    pcounter_logsumexp(terms);
}

double pcorr_log_density(double t, int w, const double* p1, const double* p2,
                         int rho_col) {
  const double nu1 = p1[emc2col::pcounter_corr::nu];
  const double nu2 = p2[emc2col::pcounter_corr::nu];
  const double cv = p1[emc2col::pcounter_corr::c];
  const double rho = p1[rho_col];
  const double x1 = t - p1[emc2col::pcounter_corr::t0];
  const double x2 = t - p2[emc2col::pcounter_corr::t0];
  const int k1 = pcounter_k_int(p1[emc2col::pcounter_corr::k]);
  const int k2 = pcounter_k_int(p2[emc2col::pcounter_corr::k]);
  if (w < 0 || w > 1 || !R_FINITE(t) || x1 <= 0.0 && w == 0 || x2 <= 0.0 && w == 1 ||
      !R_FINITE(nu1) || !R_FINITE(nu2) || nu1 <= 0 || nu2 <= 0 || cv < 0 ||
      !R_FINITE(cv) || !R_FINITE(rho) || rho < 0 || rho > 1 ||
      std::fabs(cv - p2[emc2col::pcounter_corr::c]) > 1e-12 ||
      std::fabs(rho - p2[rho_col]) > 1e-12 ||
      p1[emc2col::pcounter_corr::gamma] != 0 || p2[emc2col::pcounter_corr::gamma] != 0 ||
      p1[emc2col::pcounter_corr::omega] != 0 || p2[emc2col::pcounter_corr::omega] != 0 ||
      (rho > 0 && cv <= PC_EPS) || k1 <= 0 || k2 <= 0) return R_NegInf;
  const int kw = w == 0 ? k1 : k2, kl = w == 0 ? k2 : k1;
  const double xw = w == 0 ? x1 : x2, xl = w == 0 ? x2 : x1;
  double terms_fixed[64];
  std::vector<double> terms;
  if (xl <= 0.0) {
    return (kw - 1) * std::log(xw) - R::lgammafn(kw) +
      (w == 0 ? pcorr_log_moment(kw, 0, xw, 0, nu1, nu2, cv, rho)
              : pcorr_log_moment(0, kw, 0, xw, nu1, nu2, cv, rho));
  }
  if (kl > 64) terms.reserve(static_cast<size_t>(kl));
  for (int n = 0; n < kl; ++n) {
    double z = (kw - 1) * std::log(xw) - R::lgammafn(kw) +
      n * std::log(xl) - R::lgammafn(n + 1.0);
    z += w == 0 ? pcorr_log_moment(kw, n, x1, x2, nu1, nu2, cv, rho)
                : pcorr_log_moment(n, kw, x1, x2, nu1, nu2, cv, rho);
    if (kl <= 64) terms_fixed[n] = z;
    else terms.push_back(z);
  }
  return kl <= 64 ? pcorr_logsumexp_n(terms_fixed, kl) :
    pcounter_logsumexp(terms);
}
} // namespace

double c_log_likelihood_pcounter_corr(Rcpp::NumericMatrix pars,
                                      Rcpp::DataFrame dadm, int n_trials,
                                      Rcpp::LogicalVector winner,
                                      Rcpp::IntegerVector expand, double min_ll,
                                      const Rcpp::LogicalVector isok, int n_lR,
                                      int rho_col, Rcpp::NumericVector* trial_ll_out) {
  if (n_lR != 2 || n_trials % 2 || pars.nrow() != n_trials ||
      pars.ncol() <= rho_col || isok.size() != n_trials || winner.size() != n_trials)
    Rcpp::stop("PCOUNTERcorr requires exactly two parameter rows per trial.");
  if (dadm.containsElementNamed("RACE") || dadm.hasAttribute("pGuess") ||
      dadm.hasAttribute("pContaminant"))
    Rcpp::stop("PCOUNTERcorr does not support variable races or nuisance mixtures.");
  Rcpp::NumericVector rt = dadm["rt"];
  const int n_unique = n_trials / 2;
  std::vector<double> ll(static_cast<size_t>(n_unique), min_ll);
  std::vector<double> p1(static_cast<size_t>(pars.ncol()));
  std::vector<double> p2(static_cast<size_t>(pars.ncol()));
  for (int j = 0; j < n_unique; ++j) {
    const int at = 2 * j;
    if (!isok[at] || !isok[at + 1] || winner[at] == winner[at + 1]) continue;
    const int w = winner[at] ? 0 : 1;
    const double t = rt[at];
    for (int col = 0; col < pars.ncol(); ++col) {
      p1[static_cast<size_t>(col)] = pars(at, col);
      p2[static_cast<size_t>(col)] = pars(at + 1, col);
    }
    double z = pcorr_log_density(t, w, p1.data(), p2.data(), rho_col);
    ll[static_cast<size_t>(j)] = R_FINITE(z) && z > min_ll ? z : min_ll;
  }
  if (trial_ll_out != nullptr) {
    if (trial_ll_out->size() != expand.size()) Rcpp::stop("PCOUNTERcorr trial output size mismatch.");
    for (int j = 0; j < expand.size(); ++j) (*trial_ll_out)[j] = ll[expand[j] - 1];
  }
  double total = 0.0;
  for (int j = 0; j < expand.size(); ++j) total += ll[expand[j] - 1];
  return total;
}

// ============================================================================
// Rolling log Stirling numbers of the second kind.
// ============================================================================
class PcounterStirlingRows {
 public:
  PcounterStirlingRows() : n_(0), row_(1, 0.0) {}

  double at(int n, int j) {
    if (n < 0 || j < 0 || j > n) return R_NegInf;
    if (n < n_) {
      n_ = 0;
      row_.assign(1, 0.0);
    }
    while (n_ < n) {
      const int next = n_ + 1;
      row_.push_back(R_NegInf);
      for (int col = next; col >= 1; --col) {
        const double a = row_[static_cast<size_t>(col - 1)];
        const double log_b = (col < next && R_FINITE(row_[static_cast<size_t>(col)]))
          ? std::log(static_cast<double>(next - 1)) +
                row_[static_cast<size_t>(col)]
          : R_NegInf;
        row_[static_cast<size_t>(col)] =
          (a == R_NegInf) ? log_b :
          (log_b == R_NegInf) ? a :
          (a > log_b ? a + std::log1p(std::exp(log_b - a))
                      : log_b + std::log1p(std::exp(a - log_b)));
      }
      row_[0] = R_NegInf;
      n_ = next;
    }
    return row_[static_cast<size_t>(j)];
  }

 private:
  int n_;
  std::vector<double> row_;
};

// ============================================================================
// Closed-form kernels.
// ============================================================================

void pcounter_log_lm(double a, double nu, double sv,
                     double& logL, double& logM) {
  if (sv <= PC_EPS) {
    logL = -nu * a;
    logM = std::log(nu) + logL;
    return;
  }
  const double shape = nu * nu / (sv * sv);
  const double rate = nu / (sv * sv);
  logL = -shape * std::log1p(a / rate);
  logM = std::log(shape) - std::log(rate + a) + logL;
}
double pcounter_log_h(int n, double t, double nu, double sv, double gamma,
                      PcounterStirlingRows& rows) {
  double logL, logM;
  pcounter_log_lm(t, nu, sv, logL, logM);
  if (sv <= PC_EPS)
    return logL + pcounter_log_rising(nu / gamma, n);
  const double shape = nu * nu / (sv * sv);
  const double rate = nu / (sv * sv);
  std::vector<double> terms;
  std::vector<double> rising(static_cast<size_t>(n + 1), 0.0);
  for (int j = 1; j <= n; ++j)
    rising[static_cast<size_t>(j)] =
      rising[static_cast<size_t>(j - 1)] +
      std::log(shape + static_cast<double>(j - 1));
  terms.reserve(static_cast<size_t>(n + 1));
  for (int j = 0; j <= n; ++j) {
    const double c = rows.at(n, j);
    if (!R_FINITE(c)) { terms.push_back(R_NegInf); continue; }
    terms.push_back(c + rising[static_cast<size_t>(j)] - j * std::log(gamma) -
                       j * std::log(rate + t));
  }
  return logL + pcounter_logsumexp(terms);
}

void pcounter_log_pi_phi(int n, double t, double nu, double sv, double gamma,
                         bool gamma_zero,
                         PcounterStirlingRows& rows,
                         double& logpi, double& logphi) {
  if (gamma_zero) {
    if (sv <= PC_EPS) {
      logpi = -nu * t;
      if (n > 0) logpi += n * std::log(nu * t) - R::lgammafn(n + 1.0);
      logphi = std::log(nu) + logpi;
      return;
    }
    const double shape = nu * nu / (sv * sv);
    const double rate = nu / (sv * sv);
    double logL, logM;
    pcounter_log_lm(t, nu, sv, logL, logM);
    logpi = logL + pcounter_log_rising(shape, n) - R::lgammafn(n + 1.0);
    if (n > 0) logpi += n * std::log(t);
    logpi -= n * std::log(rate + t);
    logphi = logpi + std::log(shape + n) - std::log(rate + t);
    return;
  }
  const double q = -std::expm1(-gamma * t);
  const double logq = q > 0.0 ? std::log(q) : R_NegInf;
  logpi = n * logq - R::lgammafn(n + 1.0) +
          pcounter_log_h(n, t, nu, sv, gamma, rows);
  logphi = std::log(gamma) + n * logq - R::lgammafn(n + 1.0) +
           pcounter_log_h(n + 1, t, nu, sv, gamma, rows);
}

double pcounter_log_tail(int start, double t, double nu, double sv, double gamma,
                         double logr, bool gamma_zero,
                         PcounterStirlingRows& rows, bool phi) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back((phi ? lf : lp) + n * logr);
  }
  return pcounter_logsumexp(terms);
}

double pcounter_log_fixed_cdf(int start, double t, double nu, double sv,
                              double gamma, bool gamma_zero,
                              PcounterStirlingRows& rows) {
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back(lp);
  }
  return pcounter_logsumexp(terms);
}

double pcounter_log_geom_cdf(int start, double t, double nu, double sv,
                             double gamma, double omega, bool gamma_zero,
                             PcounterStirlingRows& rows) {
  const double logr = std::log(omega) - std::log1p(omega);
  std::vector<double> terms;
  terms.reserve(65);
  for (int n = start; n <= start + 64; ++n) {
    if (n == start) { terms.push_back(R_NegInf); continue; }
    double lp, lf;
    pcounter_log_pi_phi(n, t, nu, sv, gamma, gamma_zero, rows, lp, lf);
    terms.push_back(lp + std::log(-std::expm1((n - start) * logr)));
  }
  return pcounter_logsumexp(terms);
}

void pcounter_log_eval(double t, double nu, double sv, double gamma,
                       double k, double omega,
                       double& logf, double& logS, double& logF) {
  logf = R_NegInf; logS = 0.0; logF = R_NegInf;
  if (ISNAN(t) || !R_FINITE(nu) || !R_FINITE(sv) || !R_FINITE(gamma) ||
      !pcounter_k_supported(k) || !R_FINITE(omega) || nu <= 0.0 ||
      sv < 0.0 || gamma < 0.0 || omega < 0.0) return;
  if (R_PosInf == t) { logS = R_NegInf; logF = 0.0; return; }
  if (!R_FINITE(t) || t <= 0.0) return;

  const int kk = pcounter_k_int(k);
  const bool gamma_zero = gamma < PC_EPS;
  const bool sv_zero = sv < PC_EPS;
  const bool omega_zero = omega < PC_EPS;
  PcounterStirlingRows rows;
  if (omega_zero) {
    std::vector<double> probs;
    probs.reserve(static_cast<size_t>(kk));
    for (int n = 0; n < kk; ++n) {
      double lp, lf;
      pcounter_log_pi_phi(n, t, nu, sv_zero ? 0.0 : sv,
                          gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
      probs.push_back(lp);
    }
    logS = std::min(0.0, pcounter_logsumexp(probs));
    double lp, lf;
    pcounter_log_pi_phi(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                        gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
    logf = lf;
    logF = logS > -1e-7 ? pcounter_log_fixed_cdf(kk, t, nu,
      sv_zero ? 0.0 : sv, gamma_zero ? 0.0 : gamma, gamma_zero, rows)
      : (logS < 0.0 ? std::log(-std::expm1(logS)) : R_NegInf);
    return;
  }

  const double logp = -std::log1p(omega);
  const double logr = std::log(omega) + logp;
  const double p = std::exp(logp);
  double ar, logd;
  if (gamma_zero) {
    ar = p * t; logd = 0.0;
  } else {
    const double q = -std::expm1(-gamma * t);
    logd = std::log1p(-std::exp(logr) * q);
    ar = t + logd / gamma;
  }
  double logLar, logMar;
  pcounter_log_lm(ar, nu, sv_zero ? 0.0 : sv, logLar, logMar);
  const int nlow = kk - 2;
  std::vector<double> low_pi, low_pi_r, low_phi_r;
  if (nlow >= 0) {
    low_pi.reserve(static_cast<size_t>(nlow + 1));
    low_pi_r.reserve(static_cast<size_t>(nlow + 1));
    low_phi_r.reserve(static_cast<size_t>(nlow + 1));
    for (int n = 0; n <= nlow; ++n) {
      double lp, lf;
      pcounter_log_pi_phi(n, t, nu, sv_zero ? 0.0 : sv,
                          gamma_zero ? 0.0 : gamma, gamma_zero, rows, lp, lf);
      low_pi.push_back(lp);
      low_pi_r.push_back(lp + n * logr);
      low_phi_r.push_back(lf + n * logr);
    }
  }
  const double loglow = pcounter_logsumexp(low_pi);
  const double loglowr = pcounter_logsumexp(low_pi_r);
  double logtail = pcounter_logdiffexp(logLar, loglowr);
  if (nlow >= 0 && (loglowr > logLar - 1e-7 || !R_FINITE(logtail)))
    logtail = pcounter_log_tail(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                gamma_zero ? 0.0 : gamma, logr, gamma_zero,
                                rows, false);
  logS = std::min(0.0, pcounter_logsumexp(
      std::vector<double>{loglow, (1.0 - kk) * logr + logtail}));
  const double logA = logMar - (gamma_zero ? 0.0 : logd);
  const double loglowf = pcounter_logsumexp(low_phi_r);
  double logtail_flux = pcounter_logdiffexp(logA, loglowf);
  logf = logp + (1.0 - kk) * logr + logtail_flux;
  if (nlow >= 0 && (loglowf > logA - 1e-7 || !R_FINITE(logf))) {
    logtail_flux = pcounter_log_tail(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                     gamma_zero ? 0.0 : gamma, logr, gamma_zero,
                                     rows, true);
    logf = logp + (1.0 - kk) * logr + logtail_flux;
  }
  if (logS > -1e-7)
    logF = pcounter_log_geom_cdf(kk - 1, t, nu, sv_zero ? 0.0 : sv,
                                 gamma_zero ? 0.0 : gamma, omega, gamma_zero,
                                 rows);
  else logF = logS < 0.0 ? std::log(-std::expm1(logS)) : R_NegInf;
}

// ============================================================================
// Race-model adapters.
// ============================================================================

double dpcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j)
    if (!R_FINITE(par[j])) return 0.0;
  if (par[emc2col::pcounter::nu] <= 0.0 ||
      par[emc2col::pcounter::sv] < 0.0 ||
      par[emc2col::pcounter::gamma] < 0.0 ||
      par[emc2col::pcounter::k] < 0.0 ||
      par[emc2col::pcounter::omega] < 0.0 ||
      !pcounter_k_supported(par[emc2col::pcounter::k])) return 0.0;
  double lf, ls, lF;
  pcounter_log_eval(t - par[emc2col::pcounter::t0], par[0], par[1], par[2],
                    par[3], par[4], lf, ls, lF);
  return R_FINITE(lf) ? std::exp(lf) : 0.0;
}

double ppcounter_scalar(double t, const double* par, void* /*ctx_*/) {
  for (int j = 0; j < emc2col::pcounter::N_REQ; ++j)
    if (!R_FINITE(par[j])) return 0.0;
  if (par[0] <= 0.0 || par[1] < 0.0 || par[2] < 0.0 || par[3] < 0.0 ||
      par[4] < 0.0 || !pcounter_k_supported(par[3])) return 0.0;
  double lf, ls, lF;
  pcounter_log_eval(t - par[5], par[0], par[1], par[2], par[3], par[4],
                    lf, ls, lF);
  return lF == 0.0 ? 1.0 : (R_FINITE(lF) ? std::exp(lF) : 0.0);
}

void dpcounter_raw(const double* rt, const double* const* cols, int n_rows,
                   const int* mask, const int* isok, double* out,
                   double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* nu = cols[emc2col::pcounter::nu];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* kk = cols[emc2col::pcounter::k];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || !R_FINITE(nu[i]) || !R_FINITE(sv[i]) ||
        !R_FINITE(ga[i]) || !R_FINITE(kk[i]) || !R_FINITE(om[i]) ||
        !R_FINITE(t0[i]) || nu[i] <= 0.0 || sv[i] < 0.0 ||
        ga[i] < 0.0 || om[i] < 0.0 || !pcounter_k_supported(kk[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i],
                      lf, ls, lF);
    out[i] = raw_log_value(lf, min_ll, floor_raw);
  }
}

void ppcounter_raw(const double* rt, const double* const* cols, int n_rows,
                   const int* mask, const int* isok, double* out,
                   double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  const double* nu = cols[emc2col::pcounter::nu];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* kk = cols[emc2col::pcounter::k];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    if (!isok[i] || !R_FINITE(nu[i]) || !R_FINITE(sv[i]) ||
        !R_FINITE(ga[i]) || !R_FINITE(kk[i]) || !R_FINITE(om[i]) ||
        !R_FINITE(t0[i]) || nu[i] <= 0.0 || sv[i] < 0.0 ||
        ga[i] < 0.0 || om[i] < 0.0 || !pcounter_k_supported(kk[i])) {
      out[i] = raw_log_zero(min_ll, floor_raw); continue;
    }
    double lf, ls, lF;
    pcounter_log_eval(rt[i] - t0[i], nu[i], sv[i], ga[i], kk[i], om[i],
                      lf, ls, lF);
    // Log-survivors are NOT floored at min_ll: a loser accumulator legitimately
    // carries large negative log-survival, and clamping it flattens the
    // likelihood surface out in the tails (and would disagree with
    // pcounter_logS_at_t, which sums the unfloored value for truncation).
    // Only an exactly-zero survivor falls back to raw_log_zero, as in
    // prdm_raw / plnr_raw / prexg_raw.
    out[i] = R_FINITE(ls) ? std::fmin(ls, 0.0)
                          : raw_log_zero(min_ll, floor_raw);
  }
}

void pcounter_logS_at_t(double t, const double* const* cols,
                        int /*n_rows_total*/, int n_lR, int /*n_par*/,
                        const int* trunc_mask, int n_unique_trials,
                        const int* isok_all, void* ctx_, double* logS_out) {
  (void)ctx_;
  const double* kk = cols[emc2col::pcounter::k];
  const double* sv = cols[emc2col::pcounter::sv];
  const double* ga = cols[emc2col::pcounter::gamma];
  const double* om = cols[emc2col::pcounter::omega];
  const double* t0 = cols[emc2col::pcounter::t0];
  const double* nu = cols[emc2col::pcounter::nu];
  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    double sum = 0.0;
    for (int k = 0; k < n_lR; ++k) {
      const int r = j * n_lR + k;
      if (!isok_all[r] || !R_FINITE(nu[r]) || !R_FINITE(sv[r]) ||
          !R_FINITE(ga[r]) || !R_FINITE(kk[r]) || !R_FINITE(om[r]) ||
          !R_FINITE(t0[r]) || nu[r] <= 0.0 || sv[r] < 0.0 ||
          ga[r] < 0.0 || om[r] < 0.0 || !pcounter_k_supported(kk[r])) {
        sum = R_NegInf; break;
      }
      double lf, ls, lF;
      pcounter_log_eval(t - t0[r], nu[r], sv[r], ga[r], kk[r], om[r],
                        lf, ls, lF);
      sum += ls;
    }
    logS_out[j] = sum;
  }
}

// ============================================================================
// R-callable entry points.
// ============================================================================

// Only scalars and full-length vectors recycle.  A partial-length argument is
// an error rather than R-style wrap-around: the loops index directly, so
// silently accepting one would read past the end of the vector.
void pcounter_check_len(const NumericVector& x, int n, const char* nm) {
  if (x.size() != 1 && x.size() != n)
    Rcpp::stop("PCOUNTER: `%s` must be length 1 or length %d, not %d.",
               nm, n, static_cast<int>(x.size()));
}

double pcounter_pick(const NumericVector& x, int i) {
  return x.size() == 1 ? x[0] : x[i];
}

// R-facing scalar-vector wrappers use exactly the same kernel as the compiled
// race path.  Keeping these here removes the former hand-maintained R mirror.
// [[Rcpp::export]]
NumericVector dpcounter(NumericVector t, NumericVector nu, NumericVector sv,
                        NumericVector gamma, NumericVector k,
                        NumericVector omega, NumericVector t0,
                        bool log_out = false) {
  const int n = t.size();
  pcounter_check_len(nu, n, "nu"); pcounter_check_len(sv, n, "sv");
  pcounter_check_len(gamma, n, "gamma"); pcounter_check_len(k, n, "k");
  pcounter_check_len(omega, n, "omega"); pcounter_check_len(t0, n, "t0");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double t0i = pcounter_pick(t0, i);
    double lf = R_NegInf, ls = 0.0, lF = R_NegInf;
    if (R_FINITE(t0i)) {
      pcounter_log_eval(
        t[i] - t0i, pcounter_pick(nu, i), pcounter_pick(sv, i),
        pcounter_pick(gamma, i), pcounter_pick(k, i), pcounter_pick(omega, i),
        lf, ls, lF);
    }
    out[i] = log_out ? lf : (R_FINITE(lf) ? std::exp(lf) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector ppcounter(NumericVector t, NumericVector nu, NumericVector sv,
                        NumericVector gamma, NumericVector k,
                        NumericVector omega, NumericVector t0,
                        bool lower_tail = true, bool log_out = false) {
  const int n = t.size();
  pcounter_check_len(nu, n, "nu"); pcounter_check_len(sv, n, "sv");
  pcounter_check_len(gamma, n, "gamma"); pcounter_check_len(k, n, "k");
  pcounter_check_len(omega, n, "omega"); pcounter_check_len(t0, n, "t0");
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double t0i = pcounter_pick(t0, i);
    double lf = R_NegInf, ls = 0.0, lF = R_NegInf;
    if (R_FINITE(t0i)) {
      pcounter_log_eval(
        t[i] - t0i, pcounter_pick(nu, i), pcounter_pick(sv, i),
        pcounter_pick(gamma, i), pcounter_pick(k, i), pcounter_pick(omega, i),
        lf, ls, lF);
    }
    const double lp = lower_tail ? lF : ls;
    out[i] = log_out ? lp : (lp == 0.0 ? 1.0 :
                              (R_FINITE(lp) ? std::exp(lp) : 0.0));
  }
  return out;
}
