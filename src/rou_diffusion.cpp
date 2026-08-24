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

inline double roup_logaddexp_local(double a, double b) {
  if (a <= fperace::LOG_FLOOR) return b;
  if (b <= fperace::LOG_FLOOR) return a;
  const double hi = std::max(a, b);
  return hi + std::log1p(std::exp(std::min(a, b) - hi));
}

inline bool roup_local_keys_cpp(double v_S, double v_T, double tau_S, double tau_T,
                                double k, double B, double A, double s,
                                const fperace::BndSpec& bs,
                                fperace::Key& ks, fperace::Key& kt,
                                bool& hs, bool& ht) {
  if (!R_finite(v_S) || !R_finite(v_T) || v_S < 0.0 || v_T < 0.0)
    return false;
  hs = v_S > fpe::ROUP_DRIFT_EPS;
  ht = v_T > fpe::ROUP_DRIFT_EPS;
  if (!hs && !ht) return false;
  if (hs && !fperace::roup_key(v_S, 0.0, tau_S, 0.0, k, B, A, s, bs, ks)) return false;
  if (ht && !fperace::roup_key(0.0, v_T, 0.0, tau_T, k, B, A, s, bs, kt)) return false;
  if (!hs) ks = kt;
  if (!ht) kt = ks;
  return true;
}

inline int roup_group_cpp(const fperace::Key& p, double tt,
                          std::vector<fperace::Key>& keys,
                          std::vector<double>& horizon,
                          std::vector<std::vector<double>>& query_times) {
  int g = -1;
  for (size_t j = 0; j < keys.size(); ++j) {
    if (keys[j] == p) { g = static_cast<int>(j); break; }
  }
  if (g < 0) {
    keys.push_back(p); horizon.push_back(tt);
    query_times.push_back(std::vector<double>(1, tt));
    return static_cast<int>(keys.size() - 1);
  }
  horizon[g] = std::max(horizon[g], tt);
  query_times[g].push_back(tt);
  return g;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Alternative parameterisations -> (v, k, s).
//
// The R-side dROU/pROU/rROU keep taking (v, k, s), so this is the one call that
// R makes to reach the alternatives, and the map itself stays in fpe_race.h
// where the likelihood kernels read it.  Re-implementing the algebra in R would
// be the obvious shortcut and exactly the drift this file exists to prevent.
//
// par_kind: 0 = rate (a pass-through), 1 = curvature (tstar, k, s),
// 2 = equilibrium (tk, theta, chi).
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List rou_to_rate_vec(int par_kind, NumericVector p1, NumericVector p2,
                           NumericVector p3, NumericVector B, NumericVector A) {
  const int n = p1.size();
  if (p2.size() != n || p3.size() != n || B.size() != n || A.size() != n) {
    stop("rou_to_rate_vec: all parameter vectors must be the same length.");
  }
  NumericVector v(n), k(n), s(n);
  for (int i = 0; i < n; ++i) {
    double vv = 0.0, kk = 0.0, ss = 0.0;
    fperace::rou_map_to_rate(par_kind, p1[i], p2[i], p3[i], B[i], A[i],
                             vv, kk, ss);
    v[i] = vv; k[i] = kk; s[i] = ss;
  }
  return Rcpp::List::create(_["v"] = v, _["k"] = k, _["s"] = s);
}

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
// that had not finished by the finite simulation horizon t_max; this is useful
// for treating deadline omissions as lost draws without asserting an intrinsic
// point mass at infinity.
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
                    SEXP kind_sexp = R_NilValue, double dt = 1e-3, double t_max = 30.0,
                    int par_kind = 0) {
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

  // Map parameter column indices.  The parameterisation renames only the three
  // columns that stand in for (v, k, s); B, A, t0 and the collapse columns are
  // the same under all of them.
  const char* n1 = "v";
  const char* n2 = "k";
  const char* n3 = "s";
  if (par_kind == fperace::ROU_PAR_CURVATURE) {
    n1 = "tstar"; n2 = "k"; n3 = "s";
  } else if (par_kind == fperace::ROU_PAR_EQUILIBRIUM) {
    n1 = "tk"; n2 = "theta"; n3 = "chi";
  }
  CharacterVector col_names = colnames(pars);
  int ip1 = -1, ip2 = -1, ip3 = -1, iB = -1, iA = -1, it0 = -1;
  int iBinf = -1, itau = -1, ipw = -1;
  for (int j = 0; j < col_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(col_names[j]);
    if (nm == n1) ip1 = j;
    else if (nm == n2) ip2 = j;
    else if (nm == n3) ip3 = j;
    else if (nm == "B") iB = j;
    else if (nm == "A") iA = j;
    else if (nm == "t0") it0 = j;
    else if (nm == "Binf") iBinf = j;
    else if (nm == "tau") itau = j;
    else if (nm == "pw") ipw = j;
  }
  if (ip1 < 0 || iB < 0 || it0 < 0) {
    stop("rrou_cpp: pars matrix must contain at least '%s', 'B', and 't0' columns.",
         n1);
  }
  // Under rate, k and s have defaults; the alternatives have no meaningful
  // default for their second and third columns, so require them.
  if (par_kind != fperace::ROU_PAR_RATE && (ip2 < 0 || ip3 < 0)) {
    stop("rrou_cpp: pars matrix must contain '%s' and '%s' columns.", n2, n3);
  }

  // Resolve boundary collapse kind
  int bkind = fpe::FPE_BND_FIXED;
  if (!Rf_isNull(kind_sexp)) {
    std::string kstr = Rcpp::as<std::string>(kind_sexp);
    if (kstr == "weibull") bkind = fpe::FPE_BND_WEIBULL;
    else if (kstr == "exponential") bkind = fpe::FPE_BND_EXPONENTIAL;
    else if (kstr == "linear_additive") bkind = fpe::FPE_BND_LINEAR_ADDITIVE;
    else if (kstr == "linear_multiplicative") bkind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
  }
  const bool is_fixed = (bkind == fpe::FPE_BND_FIXED);

  IntegerVector R_out(n_trials, NA_INTEGER);
  NumericVector rt_out(n_trials, R_PosInf);

  Rcpp::RNGScope scope;

  // Trial scratch has a fixed accumulator width.  Reuse it rather than paying
  // for ten small heap allocations on every simulated trial.
  std::vector<double> X(n_acc), b(n_acc), v_acc(n_acc), k_acc(n_acc);
  std::vector<double> t0_acc(n_acc), phi(n_acc), drift_step(n_acc), sd(n_acc);
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

      double BB = pars(r, iB);
      if (!R_finite(BB)) continue;
      double AA = (iA >= 0 && R_finite(pars(r, iA)) && pars(r, iA) > 0.0) ? pars(r, iA) : 0.0;

      double vv, kk, ss;
      if (par_kind == fperace::ROU_PAR_RATE) {
        vv = pars(r, ip1);
        if (!R_finite(vv)) continue;
        kk = (ip2 >= 0 && R_finite(pars(r, ip2)) && pars(r, ip2) > 0.0) ? pars(r, ip2) : 0.0;
        ss = 1.0;
        if (ip3 >= 0) {
          ss = pars(r, ip3);
          if (!R_finite(ss) || !(ss > 0.0)) continue;
        }
      } else {
        fperace::rou_map_to_rate(par_kind, pars(r, ip1), pars(r, ip2),
                                 pars(r, ip3), BB, AA, vv, kk, ss);
        if (!R_finite(vv) || !R_finite(kk) || !R_finite(ss) || !(ss > 0.0)) continue;
        if (kk < 0.0) kk = 0.0;
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
      double drift_gain = (kk > 1e-10) ? (m1 / kk) : dt;
      drift_step[a] = vv * drift_gain;
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

        double X1 = X[a] * phi[a] + drift_step[a] + sd[a] * ::norm_rand();
        double b1 = is_fixed ? b[a] : bnd[a].b(t + dt);
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

// ---------------------------------------------------------------------------
// Pulse-based Racing OU (ROUp)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List droup_cpp(NumericVector rt, NumericVector v_S, NumericVector v_T,
                     NumericVector tau_S, NumericVector tau_T, NumericVector k,
                     NumericVector B, NumericVector A, NumericVector t0,
                     NumericVector s, int nx = 384,
                     double dt_target = 2e-3, double grade = 8.0,
                     double tgrade = 32.0, int bkind = 0,
                     NumericVector Binf = NumericVector::create(),
                     NumericVector tau = NumericVector::create(),
                     NumericVector pw = NumericVector::create(), int pooling = 0) {
  const int n = rt.size();
  if (v_S.size() != n || v_T.size() != n || tau_S.size() != n || tau_T.size() != n ||
      k.size() != n || B.size() != n || A.size() != n || t0.size() != n || s.size() != n) {
    stop("droup_cpp: all parameter vectors must match length(rt).");
  }
  check_bnd_lengths(bkind, Binf, tau, pw, n, "droup_cpp");
  if (pooling != 0 && pooling != 1) stop("droup_cpp: pooling must be 0 (coactive) or 1 (local_race).");
  for (int i = 0; i < n; ++i) {
    if (R_finite(v_S[i]) && v_S[i] < 0.0) stop("droup_cpp: negative sustained drift is invalid.");
    if (R_finite(v_T[i]) && v_T[i] < 0.0) stop("droup_cpp: negative transient drift is invalid.");
    if (R_finite(v_S[i]) && R_finite(v_T[i]) &&
        v_S[i] <= fpe::ROUP_DRIFT_EPS && v_T[i] <= fpe::ROUP_DRIFT_EPS)
      stop("droup_cpp: both pulse channels are disabled.");
  }

  NumericVector pdf(n, 0.0), cdf(n, 0.0);
  fperace::SolveCache C;
  C.roup_local = pooling == 1;
  C.grid = rou_grid(nx, dt_target, grade, tgrade);
  SEXP sparse = Rf_GetOption1(Rf_install("emc2.rou_sparse_output"));
  if (sparse != R_NilValue && Rf_length(sparse) > 0) {
    const int enabled = Rf_asLogical(sparse);
    if (enabled != NA_LOGICAL) C.sparse_raw_output = enabled;
  }

  // Pass 1: group and find each group's max horizon.  Local mode maintains one
  // mapping per pulse channel; both mappings point into the same cache entries.
  std::vector<fperace::Key> keys, keys_s, keys_t;
  std::vector<double> horizon, horizon_s, horizon_t;
  std::vector<std::vector<double>> query_times, query_s, query_t;
  std::vector<int> grp(n, -1), grp_s(n, -1), grp_t(n, -1);
  for (int i = 0; i < n; ++i) {
    const double tt = rt[i] - t0[i];
    if (!R_finite(tt) || tt <= 0.0) continue;
    const fperace::BndSpec bs = rou_bnd_at(bkind, Binf, tau, pw, i);
    if (pooling == 1) {
      fperace::Key ks, kt; bool hs = false, ht = false;
      if (!roup_local_keys_cpp(v_S[i], v_T[i], tau_S[i], tau_T[i], k[i], B[i], A[i], s[i],
                               bs, ks, kt, hs, ht)) continue;
      if (hs) grp_s[i] = roup_group_cpp(ks, tt, keys_s, horizon_s, query_s);
      if (ht) grp_t[i] = roup_group_cpp(kt, tt, keys_t, horizon_t, query_t);
    } else {
      fperace::Key p;
      if (!fperace::roup_key(v_S[i], v_T[i], tau_S[i], tau_T[i], k[i], B[i], A[i], s[i], bs, p)) continue;
      grp[i] = roup_group_cpp(p, tt, keys, horizon, query_times);
    }
  }

  for (auto& times : query_times) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }
  for (auto* qv : {&query_s, &query_t}) for (auto& times : *qv) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  // Pass 2: solve keys to max required horizon.
  std::vector<int> idx, idx_s, idx_t;
  if (pooling == 1) {
    fperace::cache_get_batch(C, keys_s, horizon_s, idx_s,
                             C.sparse_raw_output ? &query_s : nullptr);
    fperace::cache_get_batch(C, keys_t, horizon_t, idx_t,
                             C.sparse_raw_output ? &query_t : nullptr);
  } else {
    fperace::cache_get_batch(C, keys, horizon, idx,
                             C.sparse_raw_output ? &query_times : nullptr);
  }

  // Pass 3: interpolate.
  for (int i = 0; i < n; ++i) {
    const double tt = rt[i] - t0[i];
    double lp = fperace::LOG_FLOOR, lS = 0.0;
    if (pooling == 1) {
      const bool hs = grp_s[i] >= 0, ht = grp_t[i] >= 0;
      if (!hs && !ht) continue;
      const double lss = hs ? fperace::entry_log_S(C.e[idx_s[grp_s[i]]], tt) : 0.0;
      const double lst = ht ? fperace::entry_log_S(C.e[idx_t[grp_t[i]]], tt) : 0.0;
      lS = lss + lst;
      if (hs && ht) {
        lp = roup_logaddexp_local(
          fperace::entry_log_pdf(C.e[idx_s[grp_s[i]]], tt) + lst,
          fperace::entry_log_pdf(C.e[idx_t[grp_t[i]]], tt) + lss);
      } else {
        lp = hs ? fperace::entry_log_pdf(C.e[idx_s[grp_s[i]]], tt)
                : fperace::entry_log_pdf(C.e[idx_t[grp_t[i]]], tt);
      }
    } else {
      if (grp[i] < 0) continue;
      const fperace::Entry& en = C.e[idx[grp[i]]];
      lp = fperace::entry_log_pdf(en, tt);
      lS = fperace::entry_log_S(en, tt);
    }
    pdf[i] = (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
    cdf[i] = (lS >= 0.0) ? 0.0 : ((lS <= fperace::LOG_FLOOR) ? 1.0 : -std::expm1(lS));
  }
  return Rcpp::List::create(_["pdf"] = pdf, _["cdf"] = cdf, _["n_solves"] = static_cast<int>(C.e.size()));
}

// Simulate one pulse trajectory.  Local ROUp calls this once per active
// channel, giving each channel its own start point and noise stream.
double simulate_roup_channel_cpp(double vS, double vT, double tauS, double tauT,
                                 double kk, double BB, double AA, double ss,
                                 double dt, double t_max, int bkind,
                                 double Binf, double tau, double pw) {
  if (!R_finite(vS) || !R_finite(vT) || vS < 0.0 || vT < 0.0 ||
      !R_finite(BB) || !R_finite(ss) || ss <= 0.0) return R_PosInf;
  const bool as = vS > fpe::ROUP_DRIFT_EPS, at = vT > fpe::ROUP_DRIFT_EPS;
  if (!as && !at) return R_PosInf;
  if ((as && (!(tauS > 0.0) || !R_finite(tauS))) ||
      (at && (!(tauT > 0.0) || !R_finite(tauT)))) return R_PosInf;
  AA = (R_finite(AA) && AA > 0.0) ? AA : 0.0;
  fpe::FPE_Boundary bnd;
  if (bkind == fpe::FPE_BND_FIXED)
    bnd.set_kind(fpe::FPE_BND_FIXED, BB + AA, BB + AA, 0.0, 0.0, false);
  else
    bnd.set_kind(bkind, BB + AA, Binf, tau,
                 (bkind == fpe::FPE_BND_WEIBULL) ? pw : 0.0, false);
  double X = (AA > 0.0) ? AA * ::unif_rand() : 0.0;
  double b = bnd.b(0.0);
  if (X >= b) return 0.0;
  kk = (R_finite(kk) && kk > 0.0) ? kk : 0.0;
  const double phi = std::exp(-kk * dt);
  const double m1 = -std::expm1(-kk * dt);
  const double m2 = -std::expm1(-2.0 * kk * dt);
  const double drift_gain = (kk > 1e-10) ? (m1 / kk) : dt;
  const double var_val = (kk > 1e-10) ? (ss * ss * m2 / (2.0 * kk)) : (ss * ss * dt);
  const double sd_val = std::sqrt(std::max(var_val, 0.0));
  const double inv_2var_bb = 2.0 / (ss * ss * dt);
  for (double t = 0.0; t < t_max; t += dt) {
    const double tm = t + 0.5 * dt;
    const double muS = as ? vS * (1.0 - std::exp(-tm / tauS)) : 0.0;
    const double muT = at ? vT * (tm / tauT) * std::exp(-tm / tauT) : 0.0;
    const double X1 = X * phi + (muS + muT) * drift_gain + sd_val * ::norm_rand();
    const double b1 = (bkind == fpe::FPE_BND_FIXED) ? b : bnd.b(t + dt);
    if (X1 >= b1) {
      const double d0 = b - X, d1 = b1 - X1;
      const double frac = d0 / std::max(d0 - d1, 1e-300);
      return t + std::min(std::max(frac, 0.0), 1.0) * dt;
    }
    const double pc = std::exp(-(b - X) * (b1 - X1) * inv_2var_bb);
    if (::unif_rand() < pc) return t + ::unif_rand() * dt;
    X = X1; b = b1;
  }
  return R_PosInf;
}

// [[Rcpp::export]]
Rcpp::List rroup_cpp(NumericMatrix pars, CharacterVector lR_levels, LogicalVector ok,
                     SEXP kind_sexp = R_NilValue, double dt = 1e-3, double t_max = 30.0,
                     int pooling = 0) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows % n_acc != 0) {
    stop("rroup_cpp: parameter rows must form complete accumulator blocks.");
  }
  if (ok.size() != n_rows) {
    stop("rroup_cpp: length(ok) must equal nrow(pars).");
  }
  if (!(dt > 0.0) || !R_finite(dt) || !(t_max > 0.0) || !R_finite(t_max)) {
    stop("rroup_cpp: dt and t_max must be finite and positive.");
  }
  if (pooling != 0 && pooling != 1) stop("rroup_cpp: pooling must be 0 (coactive) or 1 (local_race).");
  const int n_trials = n_rows / n_acc;

  CharacterVector col_names = colnames(pars);
  int iv_S = -1, iv_T = -1, itau_S = -1, itau_T = -1;
  int ik = -1, is = -1, iB = -1, iA = -1, it0 = -1;
  int iBinf = -1, itau = -1, ipw = -1;
  for (int j = 0; j < col_names.size(); ++j) {
    std::string nm = Rcpp::as<std::string>(col_names[j]);
    if (nm == "v_S") iv_S = j;
    else if (nm == "v_T") iv_T = j;
    else if (nm == "tau_S") itau_S = j;
    else if (nm == "tau_T") itau_T = j;
    else if (nm == "k") ik = j;
    else if (nm == "s") is = j;
    else if (nm == "B") iB = j;
    else if (nm == "A") iA = j;
    else if (nm == "t0") it0 = j;
    else if (nm == "Binf") iBinf = j;
    else if (nm == "tau") itau = j;
    else if (nm == "pw") ipw = j;
  }
  if (iv_S < 0 || iv_T < 0 || itau_S < 0 || itau_T < 0 || iB < 0 || it0 < 0) {
    stop("rroup_cpp: pars matrix must contain v_S, v_T, tau_S, tau_T, B, t0 columns.");
  }

  int bkind = fpe::FPE_BND_FIXED;
  if (!Rf_isNull(kind_sexp)) {
    std::string kstr = Rcpp::as<std::string>(kind_sexp);
    if (kstr == "weibull") bkind = fpe::FPE_BND_WEIBULL;
    else if (kstr == "exponential") bkind = fpe::FPE_BND_EXPONENTIAL;
    else if (kstr == "linear_additive") bkind = fpe::FPE_BND_LINEAR_ADDITIVE;
    else if (kstr == "linear_multiplicative") bkind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
  }
  const bool is_fixed = (bkind == fpe::FPE_BND_FIXED);

  IntegerVector R_out(n_trials, NA_INTEGER);
  NumericVector rt_out(n_trials, R_PosInf);

  Rcpp::RNGScope scope;

  std::vector<double> X(n_acc), b(n_acc), k_acc(n_acc), t0_acc(n_acc), phi(n_acc), sd(n_acc);
  std::vector<double> vS_acc(n_acc), vT_acc(n_acc), tauS_acc(n_acc), tauT_acc(n_acc);
  std::vector<double> inv_2var_bb(n_acc), drift_gain(n_acc);
  std::vector<fpe::FPE_Boundary> bnd(n_acc);
  std::vector<unsigned char> active(n_acc, 0);

  for (int j = 0; j < n_trials; ++j) {
    double win_rt = R_PosInf;
    int winner = -1;

    if (pooling == 1) {
      for (int a = 0; a < n_acc; ++a) {
        const int r = j * n_acc + a;
        if (!ok[r]) continue;
        const double BB = pars(r, iB);
        if (!R_finite(BB)) continue;
        const double AA = (iA >= 0 && R_finite(pars(r, iA)) && pars(r, iA) > 0.0) ? pars(r, iA) : 0.0;
        const double vvS = pars(r, iv_S), vvT = pars(r, iv_T);
        if (R_finite(vvS) && vvS < 0.0) stop("rroup_cpp: negative sustained drift is invalid.");
        if (R_finite(vvT) && vvT < 0.0) stop("rroup_cpp: negative transient drift is invalid.");
        if (R_finite(vvS) && R_finite(vvT) &&
            vvS <= fpe::ROUP_DRIFT_EPS && vvT <= fpe::ROUP_DRIFT_EPS)
          stop("rroup_cpp: both pulse channels are disabled.");
        const double ttauS = pars(r, itau_S), ttauT = pars(r, itau_T);
        if (R_finite(vvS) && vvS > fpe::ROUP_DRIFT_EPS &&
            (!(ttauS > 0.0) || !R_finite(ttauS)))
          stop("rroup_cpp: active sustained channel needs a finite positive tau_S.");
        if (R_finite(vvT) && vvT > fpe::ROUP_DRIFT_EPS &&
            (!(ttauT > 0.0) || !R_finite(ttauT)))
          stop("rroup_cpp: active transient channel needs a finite positive tau_T.");
        const double kk = (ik >= 0 && R_finite(pars(r, ik)) && pars(r, ik) > 0.0) ? pars(r, ik) : 0.0;
        const double ss = (is >= 0 && R_finite(pars(r, is)) && pars(r, is) > 0.0) ? pars(r, is) : 1.0;
        const double tt0 = R_finite(pars(r, it0)) ? pars(r, it0) : 0.0;
        const double bi = iBinf >= 0 ? pars(r, iBinf) : 0.5;
        const double tv = itau >= 0 ? pars(r, itau) : 1.0;
        const double pv = ipw >= 0 ? pars(r, ipw) : 1.0;
        double hit;
        if (vvS > fpe::ROUP_DRIFT_EPS && vvT > fpe::ROUP_DRIFT_EPS) {
          const double hs = simulate_roup_channel_cpp(vvS, 0.0, ttauS, 1.0, kk, BB, AA, ss,
                                                      dt, t_max, bkind, bi, tv, pv);
          const double ht = simulate_roup_channel_cpp(0.0, vvT, 1.0, ttauT, kk, BB, AA, ss,
                                                      dt, t_max, bkind, bi, tv, pv);
          hit = std::min(hs, ht);
        } else {
          hit = simulate_roup_channel_cpp(vvS, vvT, ttauS, ttauT, kk, BB, AA, ss,
                                          dt, t_max, bkind, bi, tv, pv);
        }
        if (R_finite(hit) && hit + tt0 < win_rt) {
          win_rt = hit + tt0;
          winner = a + 1;
        }
      }
      if (winner > 0) { R_out[j] = winner; rt_out[j] = win_rt; }
      continue;
    }
    std::fill(active.begin(), active.end(), 0);

    double min_t0 = R_PosInf;

    for (int a = 0; a < n_acc; ++a) {
      int r = j * n_acc + a;
      if (!ok[r]) continue;

      double BB = pars(r, iB);
      if (!R_finite(BB)) continue;
      double AA = (iA >= 0 && R_finite(pars(r, iA)) && pars(r, iA) > 0.0) ? pars(r, iA) : 0.0;
      double vvS = pars(r, iv_S);
      double vvT = pars(r, iv_T);
      if (R_finite(vvS) && vvS < 0.0) stop("rroup_cpp: negative sustained drift is invalid.");
      if (R_finite(vvT) && vvT < 0.0) stop("rroup_cpp: negative transient drift is invalid.");
      if (R_finite(vvS) && R_finite(vvT) && vvS <= fpe::ROUP_DRIFT_EPS && vvT <= fpe::ROUP_DRIFT_EPS)
        stop("rroup_cpp: both pulse channels are disabled.");
      double ttauS = pars(r, itau_S);
      double ttauT = pars(r, itau_T);
      if (R_finite(vvS) && vvS > fpe::ROUP_DRIFT_EPS &&
          (!(ttauS > 0.0) || !R_finite(ttauS)))
        stop("rroup_cpp: active sustained channel needs a finite positive tau_S.");
      if (R_finite(vvT) && vvT > fpe::ROUP_DRIFT_EPS &&
          (!(ttauT > 0.0) || !R_finite(ttauT)))
        stop("rroup_cpp: active transient channel needs a finite positive tau_T.");
      double kk = (ik >= 0 && R_finite(pars(r, ik)) && pars(r, ik) > 0.0) ? pars(r, ik) : 0.0;
      double ss = (is >= 0 && R_finite(pars(r, is)) && pars(r, is) > 0.0) ? pars(r, is) : 1.0;
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
      vS_acc[a] = vvS;
      vT_acc[a] = vvT;
      tauS_acc[a] = ttauS;
      tauT_acc[a] = ttauT;
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

        // Midpoint drift evaluation
        double tm = t + 0.5 * dt;
        double mu_S = vS_acc[a] * (tauS_acc[a] > 1e-12 ? (1.0 - std::exp(-tm / tauS_acc[a])) : 1.0);
        double mu_T = vT_acc[a] * (tauT_acc[a] > 1e-12 ? (tm / tauT_acc[a]) * std::exp(-tm / tauT_acc[a]) : 0.0);
        double d_step = (mu_S + mu_T) * drift_gain[a];

        double X1 = X[a] * phi[a] + d_step + sd[a] * ::norm_rand();
        double b1 = is_fixed ? b[a] : bnd[a].b(t + dt);
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

// [[Rcpp::export]]
NumericVector rroup_hit_times_cpp(NumericVector v_S, NumericVector v_T,
                                 NumericVector tau_S, NumericVector tau_T,
                                 NumericVector k, NumericVector B,
                                 NumericVector A, NumericVector s,
                                 double dt = 1e-3, double t_max = 30.0,
                                 int bkind = 0,
                                 NumericVector Binf = NumericVector::create(),
                                 NumericVector tau = NumericVector::create(),
                                 NumericVector pw = NumericVector::create(),
                                 int pooling = 0) {
  const int n = v_S.size();
  if (v_T.size() != n || tau_S.size() != n || tau_T.size() != n ||
      k.size() != n || B.size() != n || A.size() != n || s.size() != n) {
    stop("rroup_hit_times_cpp: all parameter vectors must be the same length.");
  }
  if (!(dt > 0.0) || !(t_max > 0.0)) stop("rroup_hit_times_cpp: dt and t_max must be positive.");
  if (pooling != 0 && pooling != 1) stop("rroup_hit_times_cpp: pooling must be 0 (coactive) or 1 (local_race).");
  check_bnd_lengths(bkind, Binf, tau, pw, n, "rroup_hit_times_cpp");
  for (int i = 0; i < n; ++i) {
    if (R_finite(v_S[i]) && v_S[i] < 0.0) stop("rroup_hit_times_cpp: negative sustained drift is invalid.");
    if (R_finite(v_T[i]) && v_T[i] < 0.0) stop("rroup_hit_times_cpp: negative transient drift is invalid.");
    if (R_finite(v_S[i]) && R_finite(v_T[i]) &&
        v_S[i] <= fpe::ROUP_DRIFT_EPS && v_T[i] <= fpe::ROUP_DRIFT_EPS)
      stop("rroup_hit_times_cpp: both pulse channels are disabled.");
    if (R_finite(v_S[i]) && v_S[i] > fpe::ROUP_DRIFT_EPS &&
        (!(tau_S[i] > 0.0) || !R_finite(tau_S[i])))
      stop("rroup_hit_times_cpp: active sustained channel needs a finite positive tau_S.");
    if (R_finite(v_T[i]) && v_T[i] > fpe::ROUP_DRIFT_EPS &&
        (!(tau_T[i] > 0.0) || !R_finite(tau_T[i])))
      stop("rroup_hit_times_cpp: active transient channel needs a finite positive tau_T.");
  }

  NumericVector out(n, R_PosInf);
  auto runif = []() { return ::unif_rand(); };
  auto rnorm = []() { return ::norm_rand(); };

  Rcpp::RNGScope scope;
  for (int i = 0; i < n; ++i) {
    if (pooling == 1) {
      if (v_S[i] > fpe::ROUP_DRIFT_EPS && v_T[i] > fpe::ROUP_DRIFT_EPS) {
        const double hs = simulate_roup_channel_cpp(v_S[i], 0.0, tau_S[i], 1.0,
                                                     k[i], B[i], A[i], s[i], dt, t_max,
                                                     bkind, bkind == fpe::FPE_BND_FIXED ? 0.0 : Binf[i],
                                                     bkind == fpe::FPE_BND_FIXED ? 0.0 : tau[i],
                                                     bkind == fpe::FPE_BND_WEIBULL ? pw[i] : 0.0);
        const double ht = simulate_roup_channel_cpp(0.0, v_T[i], 1.0, tau_T[i],
                                                     k[i], B[i], A[i], s[i], dt, t_max,
                                                     bkind, bkind == fpe::FPE_BND_FIXED ? 0.0 : Binf[i],
                                                     bkind == fpe::FPE_BND_FIXED ? 0.0 : tau[i],
                                                     bkind == fpe::FPE_BND_WEIBULL ? pw[i] : 0.0);
        out[i] = std::min(hs, ht);
      } else {
        out[i] = simulate_roup_channel_cpp(v_S[i], v_T[i], tau_S[i], tau_T[i], k[i], B[i], A[i], s[i],
                                            dt, t_max, bkind,
                                            bkind == fpe::FPE_BND_FIXED ? 0.0 : Binf[i],
                                            bkind == fpe::FPE_BND_FIXED ? 0.0 : tau[i],
                                            bkind == fpe::FPE_BND_WEIBULL ? pw[i] : 0.0);
      }
      continue;
    }
    if (!R_finite(v_S[i]) || !R_finite(v_T[i]) || !R_finite(B[i]) || !R_finite(s[i]) || s[i] <= 0.0) {
      continue;
    }
    const double kk = (R_finite(k[i]) && k[i] > 0.0) ? k[i] : 0.0;
    const double AA = (R_finite(A[i]) && A[i] > 0.0) ? A[i] : 0.0;
    fpe::FPE_Boundary bnd;
    if (bkind == fpe::FPE_BND_FIXED) {
      bnd.set_kind(fpe::FPE_BND_FIXED, B[i] + AA, B[i] + AA, 0.0, 0.0, false);
    } else {
      bnd.set_kind(bkind, B[i] + AA, Binf[i], tau[i],
                   (bkind == fpe::FPE_BND_WEIBULL) ? pw[i] : 0.0, false);
    }
    
    // Custom trajectory logic for ROUp hit time since FPE race logic fperace::rou_hit_time_bnd assumes constant drift
    double X = (AA > 0.0) ? AA * runif() : 0.0;
    double b = bnd.b(0.0);
    if (X >= b) { out[i] = 0.0; continue; }
    
    double phi = std::exp(-kk * dt);
    double m1 = -std::expm1(-kk * dt);
    double m2 = -std::expm1(-2.0 * kk * dt);
    double drift_gain = (kk > 1e-10) ? (m1 / kk) : dt;
    double var_val = (kk > 1e-10) ? (s[i] * s[i] * m2 / (2.0 * kk)) : (s[i] * s[i] * dt);
    double sd_val = std::sqrt(std::max(var_val, 0.0));
    double inv_2var_bb = 2.0 / (s[i] * s[i] * dt);
    
    double t = 0.0;
    bool hit = false;
    while (t < t_max) {
      double tm = t + 0.5 * dt;
      double mu_S = v_S[i] * (tau_S[i] > 1e-12 ? (1.0 - std::exp(-tm / tau_S[i])) : 1.0);
      double mu_T = v_T[i] * (tau_T[i] > 1e-12 ? (tm / tau_T[i]) * std::exp(-tm / tau_T[i]) : 0.0);
      double d_step = (mu_S + mu_T) * drift_gain;
      
      double X1 = X * phi + d_step + sd_val * rnorm();
      double b1 = (bkind == fpe::FPE_BND_FIXED) ? b : bnd.b(t + dt);
      if (X1 >= b1) {
        double d0 = b - X, d1 = b1 - X1;
        double frac = d0 / std::max(d0 - d1, 1e-300);
        out[i] = t + std::min(std::max(frac, 0.0), 1.0) * dt;
        hit = true;
        break;
      }
      double pc = std::exp(-(b - X) * (b1 - X1) * inv_2var_bb);
      if (runif() < pc) {
        out[i] = t + runif() * dt;
        hit = true;
        break;
      }
      X = X1;
      b = b1;
      t += dt;
    }
    if (!hit) out[i] = R_PosInf;
  }
  return out;
}
