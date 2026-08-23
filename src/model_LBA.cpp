#include "model_LBA.h"
// CDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double pleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false,
                     int launch = 0) {
  // At infinite time k = 0 has the usual LBA limit; for k > 0 only drifts
  // above k*b can finish, and the m = 0 point limit inside the evaluators
  // retains that defective upper tail instead of returning one.
  return bawl_cdf_norm(t, A, b, v, sv, k, posdrift, log_out,
                       BAWL_DENOM_FLOOR, launch);
}
// PDF of the leaky ballistic accumulator.
// [[Rcpp::export]]
double dleakyba_norm(double t, double A, double b,
                     double v, double sv, double k,
                     bool posdrift = true, bool log_out = false,
                     int launch = 0) {
  return bawl_pdf_norm(t, A, b, v, sv, k, posdrift, log_out,
                       BAWL_DENOM_FLOOR, launch);
}
// [[Rcpp::export]]
NumericVector dkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0,
                             int launch = 0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    pdf[i] = dkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega, launch);
  }
  return pdf;
}
// [[Rcpp::export]]
NumericVector pkilledleakyba(NumericVector t,
                             NumericVector v, NumericVector b, NumericVector A,
                             NumericVector sv, NumericVector t0,
                             NumericVector k, NumericVector lambda_g, NumericVector lambda_k,
                             bool posdrift = true, bool log_out = false,
                             int kill_shape = 1, bool guess = false,
                             NumericVector erlang_omega = 1.0,
                             int launch = 0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++) {
    const double omega = (kill_shape <= 1) ? 1.0 :
                         (kill_shape == 2 ? 0.0 : pick(erlang_omega, i));
    cdf[i] = pkilledleakyba_norm(t[i], pick(v,i), pick(b,i), pick(A,i), pick(sv,i), pick(t0,i),
                                 pick(k,i), pick(lambda_g,i), pick(lambda_k,i),
                                 posdrift, log_out, kill_shape, guess, omega, launch);
  }
  return cdf;
}
// Vectorised R-callable wrappers (recycle scalar parameters).
// [[Rcpp::export]]
NumericVector dleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true, int launch = 0) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = dleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i),
                           posdrift, false, launch);
  return pdf;
}
// [[Rcpp::export]]
NumericVector pleakyba(NumericVector t,
                       NumericVector A, NumericVector b,
                       NumericVector v, NumericVector sv, NumericVector k,
                       bool posdrift = true, int launch = 0) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = pleakyba_norm(t[i], pick(A,i), pick(b,i), pick(v,i), pick(sv,i), pick(k,i),
                           posdrift, false, launch);
  return cdf;
}
// Standard LBA (exact k = 0 member with the legacy LBA normalizer floor),
// restored so the R-side dfun/pfun agree exactly with the C++ likelihood
// kernels, which also use LBA_DENOM_FLOOR for this model.
// [[Rcpp::export]]
NumericVector dlba(NumericVector t,
                   NumericVector A, NumericVector b,
                   NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector pdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    pdf[i] = lba_k0_pdf_norm(t[i], pick(A,i), pick(b,i), pick(v,i),
                             pick(sv,i), posdrift, log_out);
  return pdf;
}
// [[Rcpp::export]]
NumericVector plba(NumericVector t,
                   NumericVector A, NumericVector b,
                   NumericVector v, NumericVector sv,
                   bool posdrift = true, bool log_out = false) {
  int n = t.size();
  NumericVector cdf(n);
  auto pick = [](const NumericVector& vec, int i) -> double {
    return vec.size() == 1 ? vec[0] : vec[i];
  };
  for (int i = 0; i < n; i++)
    cdf[i] = lba_k0_cdf_norm(t[i], pick(A,i), pick(b,i), pick(v,i),
                             pick(sv,i), posdrift, log_out);
  return cdf;
}
// --------------------------------------------------------------------------
// R-callable entry points.  These bypass ContextForRaceModels entirely, so the
// launch distribution MUST be passed explicitly; R/model_BAwD.R derives both
// this argument and the c_name suffix from one `drift_distribution` value so
// dfun/pfun cannot silently disagree with the sampled likelihood.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false, double gamma = 0.0,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out, gamma, rho);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawd(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    NumericVector ell, int launch = 1, bool posdrift = true,
                    bool log_out = false, double gamma = 0.0,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), pick(ell, i), launch,
                           posdrift, log_out, gamma, rho);
  return out;
}

// [[Rcpp::export]]
double dbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false, double gamma = 0.0,
                  double rho = 0.0) {
  return bawd_pdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out,
                       gamma, rho);
}

// [[Rcpp::export]]
double pbawd_norm(double t, double A, double b, double p1, double p2, double k,
                  double ell, int launch = 1, bool posdrift = true,
                  bool log_out = false, double gamma = 0.0,
                  double rho = 0.0) {
  return bawd_cdf_norm(t, A, b, p1, p2, k, ell, launch, posdrift, log_out,
                       gamma, rho);
}

// Right endpoint of the supported decision-time window (Inf when ell = 0 or
// k = 0, where nothing saturates).  Exported because the flat-CDF contract and
// the support constraint are both stated in terms of it, and a test that
// recomputed it in R would not be testing the same quantity.
// [[Rcpp::export]]
double bawd_tmax(double A, double b, double k, double ell, double gamma = 0.0,
                 double rho = 0.0) {
  const BawdGeom g = bawd_geometry(A, b, k, ell, gamma, rho);
  return g.ok ? g.T_max : NA_REAL;
}

// Vectorised bawd_tmax, for the R Ttransform: it derives T_max per accumulator
// row, so a per-row .Call would dominate mapped_pars() and make_data().  Length
// follows A and the rest recycle, as in dbawd/pbawd above.
// [[Rcpp::export]]
NumericVector bawd_tmax_vec(NumericVector A, NumericVector b, NumericVector k,
                            NumericVector ell, NumericVector gamma = 0.0,
                            NumericVector rho = 0.0) {
  const int n = A.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawdGeom g = bawd_geometry(A[i], pick(b, i), pick(k, i),
                                     pick(ell, i), pick(gamma, i), pick(rho, i));
    out[i] = g.ok ? g.T_max : NA_REAL;
  }
  return out;
}

// Inverse of bawd_tmax_vec: the clearance rate implied by a sampled endpoint.
// The R Ttransform reports `ell` alongside `Tmax` so that the wrappers, the
// simulator and mapped_pars() all keep seeing the mechanistic parameter, and
// so that a fit in either chart can be read in the other.  Length follows
// Tmax; the rest recycle.
// [[Rcpp::export]]
NumericVector bawd_ell_vec(NumericVector Tmax, NumericVector b, NumericVector k,
                           NumericVector gamma = 0.0, NumericVector rho = 0.0) {
  const int n = Tmax.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawd_ell_from_tmax(pick(b, i), pick(k, i), Tmax[i],
                                pick(gamma, i), pick(rho, i));
  return out;
}

// log E[(V - v)_+] for log V ~ N(mu, sigma^2); exposed so the cancellation
// layer can be checked against a high-precision reference from R.
// [[Rcpp::export]]
double lognormal_stoploss_log(double v, double mu, double sigma) {
  return log_lognormal_stoploss(v, mu, sigma);
}
// log Lambda_m(v) for the lognormal negative-moment survivor integral; for
// BAwD's frozen branch m = r - 1.  The export is retained for independent
// checks of the frozen-start closed form.
// [[Rcpp::export]]
double lognormal_power_stoploss_log(double v, double mu, double sigma, double m) {
  return log_lognormal_power_stoploss(v, mu, sigma, m);
}

// log int_v^Inf log(w / ell) P(V >= w) dw.
// [[Rcpp::export]]
double lognormal_logratio_stoploss_log(double v, double mu, double sigma,
                                       double log_ell) {
  return log_lognormal_logratio_stoploss(v, mu, sigma, log_ell);
}
// --------------------------------------------------------------------------
// R-callable entry points.  As for BAwD these bypass ContextForRaceModels, so
// the launch distribution MUST be passed explicitly; R/model_BAwF.R derives
// both this argument and the c_name suffix from one `drift_distribution`
// value so dfun/pfun cannot silently disagree with the sampled likelihood.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawf(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    int launch = 1, bool posdrift = true, bool log_out = false,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawf_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), launch, posdrift, log_out,
                           rho);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawf(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector k,
                    int launch = 1, bool posdrift = true, bool log_out = false,
                    double rho = 0.0) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawf_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(k, i), launch, posdrift, log_out,
                           rho);
  return out;
}

// [[Rcpp::export]]
double dbawf_norm(double t, double A, double b, double p1, double p2, double k,
                  int launch = 1, bool posdrift = true, bool log_out = false,
                  double rho = 0.0) {
  return bawf_pdf_norm(t, A, b, p1, p2, k, launch, posdrift, log_out, rho);
}

// [[Rcpp::export]]
double pbawf_norm(double t, double A, double b, double p1, double p2, double k,
                  int launch = 1, bool posdrift = true, bool log_out = false,
                  double rho = 0.0) {
  return bawf_cdf_norm(t, A, b, p1, p2, k, launch, posdrift, log_out, rho);
}

// Right endpoint of the supported decision-time window: x_max / k, with
// x_max = 1 for the exponential kernel and rho/(rho - 1) for finite rho.  It
// is Inf only at k = 0.  Exported because the entire point of BAwF is that
// this quantity is free of b, and a test that recomputed it in R would not be
// testing the same code the likelihood uses.
// [[Rcpp::export]]
double bawf_tmax(double A, double b, double k, double rho = 0.0) {
  const BawfGeom g = bawf_geometry(A, b, k, rho);
  return g.ok ? g.T_max : NA_REAL;
}

// Vectorised bawf_tmax for the R Ttransform, which derives T_max per
// accumulator row; a per-row .Call would dominate mapped_pars()/make_data().
// [[Rcpp::export]]
NumericVector bawf_tmax_vec(NumericVector A, NumericVector b, NumericVector k,
                            NumericVector rho = 0.0) {
  const int n = A.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawfGeom g = bawf_geometry(A[i], pick(b, i), pick(k, i),
                                     pick(rho, i));
    out[i] = g.ok ? g.T_max : NA_REAL;
  }
  return out;
}

// The critical launch strength at the lowest start point: V_c(0) = k b H'(x_max),
// which is e k b for the exponential kernel.  Given the time scale set by
// T_max, this is the observable that identifies b, so it is exported for the
// same reason bawf_tmax is.
// [[Rcpp::export]]
NumericVector bawf_vcrit_vec(NumericVector A, NumericVector b, NumericVector k,
                             NumericVector rho = 0.0) {
  const int n = A.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawfGeom g = bawf_geometry(A[i], pick(b, i), pick(k, i),
                                     pick(rho, i));
    out[i] = g.ok ? g.V_c0 : NA_REAL;
  }
  return out;
}
// --------------------------------------------------------------------------
// R-callable entry points.  As for BAwD and BAwF these bypass
// ContextForRaceModels, so the launch distribution MUST be passed explicitly;
// R/model_BAwR.R derives both this argument and the c_name suffix from one
// `drift_distribution` value so dfun/pfun cannot silently disagree with the
// sampled likelihood.
// --------------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector dbawr(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector kappa,
                    NumericVector pw, int launch = 1, bool posdrift = true,
                    bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawr_pdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(kappa, i), pick(pw, i), launch,
                           posdrift, log_out);
  return out;
}

// [[Rcpp::export]]
NumericVector pbawr(NumericVector t, NumericVector A, NumericVector b,
                    NumericVector p1, NumericVector p2, NumericVector kappa,
                    NumericVector pw, int launch = 1, bool posdrift = true,
                    bool log_out = false) {
  const int n = t.size();
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawr_cdf_norm(t[i], pick(A, i), pick(b, i), pick(p1, i),
                           pick(p2, i), pick(kappa, i), pick(pw, i), launch,
                           posdrift, log_out);
  return out;
}

// [[Rcpp::export]]
double dbawr_norm(double t, double A, double b, double p1, double p2,
                  double kappa, double pw, int launch = 1,
                  bool posdrift = true, bool log_out = false) {
  return bawr_pdf_norm(t, A, b, p1, p2, kappa, pw, launch, posdrift, log_out);
}

// [[Rcpp::export]]
double pbawr_norm(double t, double A, double b, double p1, double p2,
                  double kappa, double pw, int launch = 1,
                  bool posdrift = true, bool log_out = false) {
  return bawr_cdf_norm(t, A, b, p1, p2, kappa, pw, launch, posdrift, log_out);
}

// Right endpoint of the supported decision-time window,
// T_max = [(p+1) b / (kappa p)]^(1/(p+1)); Inf only at kappa = 0.  Unlike
// BAwF's it DOES depend on b, which is the substantive difference between the
// two models, so it is exported and reported by Ttransform: a test or a
// summary that recomputed it in R would not be testing the same code the
// likelihood uses.
// [[Rcpp::export]]
double bawr_tmax(double A, double b, double kappa, double pw) {
  const BawrGeom g = bawr_geometry(A, b, kappa, pw);
  return g.ok ? g.T_max : NA_REAL;
}

// Vectorised bawr_tmax for the R Ttransform, which derives T_max per
// accumulator row; a per-row .Call would dominate mapped_pars()/make_data().
// [[Rcpp::export]]
NumericVector bawr_tmax_vec(NumericVector A, NumericVector b,
                            NumericVector kappa, NumericVector pw) {
  // Sized by the LONGEST argument, not by A: Ttransform passes full columns,
  // but a scalar A with vector b (a caution sweep) would otherwise silently
  // return a single value.
  const int n = std::max(std::max(A.size(), b.size()),
                         std::max(kappa.size(), pw.size()));
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawrGeom g = bawr_geometry(pick(A, i), pick(b, i), pick(kappa, i),
                                     pick(pw, i));
    out[i] = g.ok ? g.T_max : NA_REAL;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector bawr_kappa_vec(NumericVector Tmax, NumericVector b,
                             NumericVector pw) {
  const int n = std::max(Tmax.size(), std::max(b.size(), pw.size()));
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i)
    out[i] = bawr_kappa_from_tmax(pick(b, i), pick(pw, i), pick(Tmax, i));
  return out;
}

// The critical launch strength at the lowest start point,
// V_c(0) = kappa T_max^p.  Launches below it never reach the threshold, so
// this is the observable that sets the omission rate.
// [[Rcpp::export]]
NumericVector bawr_vcrit_vec(NumericVector A, NumericVector b,
                             NumericVector kappa, NumericVector pw) {
  // Sized by the longest argument; see bawr_tmax_vec().
  const int n = std::max(std::max(A.size(), b.size()),
                         std::max(kappa.size(), pw.size()));
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) -> double {
    return x.size() == 1 ? x[0] : x[i];
  };
  for (int i = 0; i < n; ++i) {
    const BawrGeom g = bawr_geometry(pick(A, i), pick(b, i), pick(kappa, i),
                                     pick(pw, i));
    out[i] = g.ok ? g.V_c0 : NA_REAL;
  }
  return out;
}
// [[Rcpp::export]]
NumericVector dbtawl_transient(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector tau, int launch = 1,
                     bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_pdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i),
                                    pick(p2,i), pick(k,i), pick(tau,i), launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawl_transient(NumericVector t, NumericVector A, NumericVector b,
                     NumericVector p1, NumericVector p2, NumericVector k,
                     NumericVector tau, int launch = 1,
                     bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_cdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i),
                                    pick(p2,i), pick(k,i), pick(tau,i), launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_transient_log_surv_vec(NumericVector t, NumericVector A,
                                 NumericVector b, NumericVector p1,
                                 NumericVector p2, NumericVector k,
                                 NumericVector tau, int launch = 1,
                                 bool posdrift = true) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau.size()});
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i)
    out[i] = btawl_log_surv(pick(t, i), pick(A, i), pick(b, i), pick(p1, i),
                            pick(p2, i), pick(k, i), pick(tau, i), launch,
                            posdrift);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_local_race_log_surv_vec(NumericVector t, NumericVector A,
                                     NumericVector b, NumericVector p1,
                                     NumericVector p2, NumericVector k,
                                     NumericVector tau_s, NumericVector tau_t,
                                     NumericVector pi, int launch = 1,
                                     bool posdrift = true) {
  const int n = std::max({t.size(), A.size(), b.size(), p1.size(), p2.size(),
                          k.size(), tau_s.size(), tau_t.size(), pi.size()});
  NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i)
    out[i] = btawl_local_race_log_surv(pick(t, i), pick(A, i), pick(b, i), pick(p1, i),
                                pick(p2, i), pick(k, i), pick(tau_s, i),
                                pick(tau_t, i), pick(pi, i), launch, posdrift);
  return out;
}

// [[Rcpp::export]]
NumericVector dbtawl_local_race(NumericVector t, NumericVector A, NumericVector b,
                        NumericVector p1, NumericVector p2, NumericVector k,
                        NumericVector tau_s, NumericVector tau_t,
                        NumericVector pi, int launch = 1,
                        bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_local_race_pdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                        pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                        launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pbtawl_local_race(NumericVector t, NumericVector A, NumericVector b,
                        NumericVector p1, NumericVector p2, NumericVector k,
                        NumericVector tau_s, NumericVector tau_t,
                        NumericVector pi, int launch = 1,
                        bool posdrift = true, bool log_out = false) {
  const int n = t.size(); NumericVector out(n);
  auto pick = [](const NumericVector& x, int i) { return x.size() == 1 ? x[0] : x[i]; };
  for (int i = 0; i < n; ++i) {
    const double lp = btawl_local_race_cdf_log(t[i], pick(A,i), pick(b,i), pick(p1,i), pick(p2,i),
                                        pick(k,i), pick(tau_s,i), pick(tau_t,i), pick(pi,i),
                                        launch, posdrift);
    out[i] = log_out ? lp : (lp > R_NegInf ? std::exp(lp) : 0.0);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_tmax_vec(NumericVector k, NumericVector tau) {
  const int n = std::max(k.size(), tau.size()); NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = btawl_tmax(k.size() == 1 ? k[0] : k[i],
                                                   tau.size() == 1 ? tau[0] : tau[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_tau_vec(NumericVector k, NumericVector Ttrans) {
  const int n = std::max(k.size(), Ttrans.size()); NumericVector out(n);
  for (int i = 0; i < n; ++i)
    out[i] = btawl_tau_from_ttrans(k.size() == 1 ? k[0] : k[i],
                                   Ttrans.size() == 1 ? Ttrans[0] : Ttrans[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector btawl_vcrit_vec(NumericVector k, NumericVector tau,
                              NumericVector b) {
  const int n = std::max({k.size(), tau.size(), b.size()}); NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    const double ki = k.size() == 1 ? k[0] : k[i];
    const double ti = tau.size() == 1 ? tau[0] : tau[i];
    const double bi = b.size() == 1 ? b[0] : b[i];
    const double tm = btawl_tmax(ki, ti);
    const double h = btawl_h(tm, ki, ti);
    out[i] = (h > 0.0 && emc2_isfinite(h)) ? bi / h : R_PosInf;
  }
  return out;
}
