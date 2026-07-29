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
  SEXP sparse = Rf_GetOption1(Rf_install("emc2.rou_sparse_output"));
  if (sparse != R_NilValue && Rf_length(sparse) > 0) {
    const int enabled = Rf_asLogical(sparse);
    if (enabled != NA_LOGICAL) C.sparse_raw_output = enabled;
  }

  // Pass 1: group and find each group's horizon.
  std::vector<fperace::Key> keys;
  std::vector<double> horizon;
  std::vector<std::vector<double>> query_times;
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
      query_times.push_back(std::vector<double>(1, tt));
      g = static_cast<int>(keys.size()) - 1;
    } else {
      if (tt > horizon[g]) horizon[g] = tt;
      query_times[g].push_back(tt);
    }
    grp[i] = g;
  }

  for (auto& times : query_times) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  // Pass 2: solve.  This public vector path can contain several distinct
  // parameter tuples just like the likelihood path, so use the same SIMD
  // cache fill rather than solving each tuple serially.
  std::vector<int> idx;
  fperace::cache_get_batch(
      C, keys, horizon, idx,
      C.sparse_raw_output ? &query_times : nullptr);

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

// ---------------------------------------------------------------------------
// Optimized Batched Race Simulator with Early Exit Across Accumulators
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rrou_cpp(NumericMatrix pars, CharacterVector lR_levels, LogicalVector ok,
                    SEXP kind_sexp = R_NilValue, double dt = 1e-3, double t_max = 30.0) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows % n_acc != 0) {
    stop("rrou_cpp: parameter rows must form complete accumulator blocks.");
  }
  if (ok.size() != n_rows) {
    stop("rrou_cpp: length(ok) must equal nrow(pars).");
  }
  if (!(dt > 0.0) || !R_finite(dt) ||
      !(t_max > 0.0) || !R_finite(t_max)) {
    stop("rrou_cpp: dt and t_max must be finite and positive.");
  }
  const int n_trials = n_rows / n_acc;

  // Map parameter column indices
  CharacterVector col_names = colnames(pars);
  int iv = -1, ik = -1, iB = -1, iA = -1, it0 = -1, is_ = -1;
  int iBinf = -1, itau = -1, ipw = -1;
  for (int j = 0; j < col_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(col_names[j]);
    if (nm == "v") iv = j;
    else if (nm == "k") ik = j;
    else if (nm == "B") iB = j;
    else if (nm == "A") iA = j;
    else if (nm == "t0") it0 = j;
    else if (nm == "s") is_ = j;
    else if (nm == "Binf") iBinf = j;
    else if (nm == "tau") itau = j;
    else if (nm == "pw") ipw = j;
  }
  if (iv < 0 || iB < 0 || it0 < 0) {
    stop("rrou_cpp: pars matrix must contain at least 'v', 'B', and 't0' columns.");
  }

  // Resolve boundary collapse kind
  int bkind = fpe::FPE_BND_FIXED;
  if (!Rf_isNull(kind_sexp)) {
    std::string kstr = Rcpp::as<std::string>(kind_sexp);
    if (kstr == "weibull") bkind = fpe::FPE_BND_WEIBULL;
    else if (kstr == "exponential") bkind = fpe::FPE_BND_EXPONENTIAL;
    else if (kstr == "linear") bkind = fpe::FPE_BND_LINEAR;
  }

  IntegerVector R_out(n_trials, NA_INTEGER);
  NumericVector rt_out(n_trials, R_PosInf);

  Rcpp::RNGScope scope;

  // Trial scratch has a fixed accumulator width.  Reuse it rather than paying
  // for ten small heap allocations on every simulated trial.
  std::vector<double> X(n_acc), b(n_acc), v_acc(n_acc), k_acc(n_acc);
  std::vector<double> t0_acc(n_acc), phi(n_acc), drift_gain(n_acc), sd(n_acc);
  std::vector<double> inv_2var_bb(n_acc);
  std::vector<fpe::FPE_Boundary> bnd(n_acc);
  std::vector<unsigned char> active(n_acc, 0);

  for (int j = 0; j < n_trials; ++j) {
    double win_rt = R_PosInf;
    int winner = -1;
    std::fill(active.begin(), active.end(), 0);

    double min_t0 = R_PosInf;

    for (int a = 0; a < n_acc; ++a) {
      int r = j * n_acc + a;
      if (!ok[r]) continue;

      double vv = pars(r, iv);
      double BB = pars(r, iB);
      if (!R_finite(vv) || !R_finite(BB)) continue;

      double kk = (ik >= 0 && R_finite(pars(r, ik)) && pars(r, ik) > 0.0) ? pars(r, ik) : 0.0;
      double AA = (iA >= 0 && R_finite(pars(r, iA)) && pars(r, iA) > 0.0) ? pars(r, iA) : 0.0;
      double ss = 1.0;
      if (is_ >= 0) {
        ss = pars(r, is_);
        if (!R_finite(ss) || !(ss > 0.0)) continue;
      }
      double tt0 = R_finite(pars(r, it0)) ? pars(r, it0) : 0.0;

      double Binf_val = (iBinf >= 0) ? pars(r, iBinf) : 0.5;
      double tau_val  = (itau >= 0)  ? pars(r, itau) : 1.0;
      double pw_val   = (ipw >= 0)   ? pars(r, ipw)  : 1.0;

      X[a] = (AA > 0.0) ? AA * ::unif_rand() : 0.0;
      if (bkind == fpe::FPE_BND_FIXED) {
        bnd[a].set_kind(fpe::FPE_BND_FIXED, BB + AA, BB + AA, 0.0, 0.0, false);
      } else {
        bnd[a].set_kind(bkind, BB + AA, Binf_val, tau_val,
                        (bkind == fpe::FPE_BND_WEIBULL) ? pw_val : 0.0, false);
      }
      b[a] = bnd[a].b(0.0);
      t0_acc[a] = tt0;
      v_acc[a] = vv;
      k_acc[a] = kk;

      if (tt0 < min_t0) min_t0 = tt0;

      if (X[a] >= b[a]) {
        if (tt0 < win_rt) {
          win_rt = tt0;
          winner = a + 1;
        }
        continue;
      }

      active[a] = 1;
      phi[a] = std::exp(-kk * dt);
      double m1 = -std::expm1(-kk * dt);
      double m2 = -std::expm1(-2.0 * kk * dt);
      drift_gain[a] = (kk > 1e-10) ? (m1 / kk) : dt;
      double var_val = (kk > 1e-10) ? (ss * ss * m2 / (2.0 * kk)) : (ss * ss * dt);
      sd[a] = std::sqrt(std::max(var_val, 0.0));
      inv_2var_bb[a] = 2.0 / (ss * ss * dt);
    }

    double t = 0.0;
    while (t < t_max) {
      if (t + min_t0 >= win_rt) break;

      bool any_active = false;
      for (int a = 0; a < n_acc; ++a) {
        if (!active[a]) continue;
        if (t + t0_acc[a] >= win_rt) {
          active[a] = 0;
          continue;
        }
        any_active = true;

        double X1 = X[a] * phi[a] + v_acc[a] * drift_gain[a] + sd[a] * ::norm_rand();
        double b1 = bnd[a].fixed ? b[a] : bnd[a].b(t + dt);
        if (X1 >= b1) {
          double d0 = b[a] - X[a], d1 = b1 - X1;
          double frac = d0 / std::max(d0 - d1, 1e-300);
          double hit_t = t + std::min(std::max(frac, 0.0), 1.0) * dt + t0_acc[a];
          if (hit_t < win_rt) {
            win_rt = hit_t;
            winner = a + 1;
          }
          active[a] = 0;
          continue;
        }
        double pc = std::exp(-(b[a] - X[a]) * (b1 - X1) * inv_2var_bb[a]);
        if (::unif_rand() < pc) {
          double hit_t = t + ::unif_rand() * dt + t0_acc[a];
          if (hit_t < win_rt) {
            win_rt = hit_t;
            winner = a + 1;
          }
          active[a] = 0;
          continue;
        }
        X[a] = X1;
        b[a] = b1;
      }
      if (!any_active) break;
      t += dt;
    }

    if (winner > 0) {
      R_out[j] = winner;
      rt_out[j] = win_rt;
    }
  }

  return Rcpp::List::create(Rcpp::Named("R") = R_out, Rcpp::Named("rt") = rt_out);
}
