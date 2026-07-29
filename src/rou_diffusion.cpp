// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// R entry points for the racing Ornstein-Uhlenbeck model (ROU).
//
// These back dROU / pROU / rROU in R/model_ROU.R, which serve make_data(),
// predict() and the plotting path.  They share fpe_race.h with the C++
// likelihood kernels -- the same key construction, the same solve cache, the
// same log-space interpolation -- so the R-side density and the sampled
// likelihood cannot drift apart.
//
// This TU deliberately does NOT include utils.h: the kernels there need the
// race context, and none of that is required to answer a plain vector of RTs.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <Rcpp.h>
#include "fpe_race.h"

using namespace Rcpp;

namespace {

fperace::FPE_Grid rou_grid(int nx, double dt_target, double grade, double tgrade) {
  fperace::FPE_Grid g;
  if (nx > 8) g.nx = nx;
  if (dt_target > 0.0) g.dt_target = dt_target;
  g.grade = grade;
  g.tgrade = tgrade;
  return g;
}

// Boundary spec for row i.  bkind == 0 (fixed) ignores the shape vectors
// entirely, so the fixed-bound calls need not supply them at all.
fperace::BndSpec rou_bnd_at(int bkind, const NumericVector& Binf,
                            const NumericVector& tau, const NumericVector& pw,
                            int i) {
  fperace::BndSpec bs;
  bs.kind = bkind;
  if (bkind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = Binf[i];
  bs.tau  = tau[i];
  bs.pw   = (bkind == fpe::FPE_BND_WEIBULL) ? pw[i] : 0.0;
  return bs;
}

void check_bnd_lengths(int bkind, const NumericVector& Binf,
                       const NumericVector& tau, const NumericVector& pw,
                       int n, const char* who) {
  if (bkind == fpe::FPE_BND_FIXED) return;
  if (Binf.size() != n || tau.size() != n ||
      (bkind == fpe::FPE_BND_WEIBULL && pw.size() != n)) {
    stop("%s: collapsing boundary needs Binf/tau (and pw for Weibull) at length(rt).", who);
  }
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Density and CDF for a single accumulator, one parameter set per RT.
//
// Grouped exactly as the batch kernels group: distinct parameter tuples become
// distinct solves, each taken to the largest time any of its rows needs.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rou_pdf_cdf_vec(NumericVector rt, NumericVector v, NumericVector k,
                           NumericVector B, NumericVector A, NumericVector t0,
                           NumericVector s, int nx = 384,
                           double dt_target = 2e-3, double grade = 8.0,
                           double tgrade = 32.0, int bkind = 0,
                           NumericVector Binf = NumericVector::create(),
                           NumericVector tau = NumericVector::create(),
                           NumericVector pw = NumericVector::create()) {
  const int n = rt.size();
  if (v.size() != n || k.size() != n || B.size() != n || A.size() != n ||
      t0.size() != n || s.size() != n) {
    stop("rou_pdf_cdf_vec: all parameter vectors must match length(rt).");
  }
  check_bnd_lengths(bkind, Binf, tau, pw, n, "rou_pdf_cdf_vec");

  NumericVector pdf(n, 0.0), cdf(n, 0.0);
  fperace::SolveCache C;
  C.grid = rou_grid(nx, dt_target, grade, tgrade);

  // Pass 1: group and find each group's horizon.
  std::vector<fperace::Key> keys;
  std::vector<double> horizon;
  std::vector<int> grp(n, -1);
  for (int i = 0; i < n; ++i) {
    const double tt = rt[i] - t0[i];
    if (!R_finite(tt) || tt <= 0.0) continue;
    fperace::Key p;
    if (!fperace::rou_key(v[i], k[i], B[i], A[i], s[i],
                          rou_bnd_at(bkind, Binf, tau, pw, i), p)) continue;
    int g = -1;
    for (size_t j = 0; j < keys.size(); ++j) {
      if (keys[j] == p) { g = static_cast<int>(j); break; }
    }
    if (g < 0) {
      keys.push_back(p);
      horizon.push_back(tt);
      g = static_cast<int>(keys.size()) - 1;
    } else if (tt > horizon[g]) {
      horizon[g] = tt;
    }
    grp[i] = g;
  }

  // Pass 2: solve.
  std::vector<int> idx(keys.size(), -1);
  for (size_t j = 0; j < keys.size(); ++j) {
    idx[j] = fperace::cache_get(C, keys[j], horizon[j]);
  }

  // Pass 3: interpolate.
  for (int i = 0; i < n; ++i) {
    if (grp[i] < 0) continue;
    const fperace::Entry& en = C.e[idx[grp[i]]];
    const double tt = rt[i] - t0[i];
    const double lp = fperace::entry_log_pdf(en, tt);
    pdf[i] = (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
    const double lS = fperace::entry_log_S(en, tt);
    cdf[i] = (lS >= 0.0) ? 0.0
                         : ((lS <= fperace::LOG_FLOOR) ? 1.0 : -std::expm1(lS));
  }

  return Rcpp::List::create(_["pdf"] = pdf, _["cdf"] = cdf,
                            _["n_solves"] = static_cast<int>(C.e.size()));
}

// ---------------------------------------------------------------------------
// Reference simulator, one draw per element.  Returns Inf for an accumulator
// that had not finished by t_max, which is the correct answer for a leaky
// accumulator whose asymptote v/k sits below threshold -- the race machinery
// treats Inf as "lost".
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector rou_hit_times_vec(NumericVector v, NumericVector k, NumericVector B,
                                NumericVector A, NumericVector s,
                                double dt = 1e-3, double t_max = 30.0,
                                int bkind = 0,
                                NumericVector Binf = NumericVector::create(),
                                NumericVector tau = NumericVector::create(),
                                NumericVector pw = NumericVector::create()) {
  const int n = v.size();
  if (k.size() != n || B.size() != n || A.size() != n || s.size() != n) {
    stop("rou_hit_times_vec: all parameter vectors must be the same length.");
  }
  if (!(dt > 0.0) || !(t_max > 0.0)) {
    stop("rou_hit_times_vec: dt and t_max must be positive.");
  }
  check_bnd_lengths(bkind, Binf, tau, pw, n, "rou_hit_times_vec");

  NumericVector out(n, R_PosInf);
  auto runif = []() { return ::unif_rand(); };
  auto rnorm = []() { return ::norm_rand(); };

  Rcpp::RNGScope scope;
  for (int i = 0; i < n; ++i) {
    if (!R_finite(v[i]) || !R_finite(B[i]) || !R_finite(s[i]) || s[i] <= 0.0) {
      continue;
    }
    const double kk = (R_finite(k[i]) && k[i] > 0.0) ? k[i] : 0.0;
    const double AA = (R_finite(A[i]) && A[i] > 0.0) ? A[i] : 0.0;
    // Same mapping as fperace::rou_key: b0 = B + A, asymptote Binf measured
    // from zero, so the bound may descend into the start-point range.
    fpe::FPE_Boundary bnd;
    if (bkind == fpe::FPE_BND_FIXED) {
      bnd.set_kind(fpe::FPE_BND_FIXED, B[i] + AA, B[i] + AA, 0.0, 0.0, false);
    } else {
      bnd.set_kind(bkind, B[i] + AA, Binf[i], tau[i],
                   (bkind == fpe::FPE_BND_WEIBULL) ? pw[i] : 0.0, false);
    }
    out[i] = fperace::rou_hit_time_bnd(v[i], kk, s[i], AA, bnd, dt, t_max,
                                       runif, rnorm);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Standard EMC2 exported aliases
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List drou_cpp(NumericVector rt, NumericVector v, NumericVector k,
                    NumericVector B, NumericVector A, NumericVector t0,
                    NumericVector s, int nx = 384,
                    double dt_target = 2e-3, double grade = 8.0,
                    double tgrade = 32.0, int bkind = 0,
                    NumericVector Binf = NumericVector::create(),
                    NumericVector tau = NumericVector::create(),
                    NumericVector pw = NumericVector::create()) {
  return rou_pdf_cdf_vec(rt, v, k, B, A, t0, s, nx, dt_target, grade, tgrade,
                         bkind, Binf, tau, pw);
}

// [[Rcpp::export]]
NumericVector rrou_hit_times_cpp(NumericVector v, NumericVector k, NumericVector B,
                                 NumericVector A, NumericVector s,
                                 double dt = 1e-3, double t_max = 30.0,
                                 int bkind = 0,
                                 NumericVector Binf = NumericVector::create(),
                                 NumericVector tau = NumericVector::create(),
                                 NumericVector pw = NumericVector::create()) {
  return rou_hit_times_vec(v, k, B, A, s, dt, t_max, bkind, Binf, tau, pw);
}

