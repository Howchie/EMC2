// C++ simulation kernels for posterior prediction, one entry point per race
// family (LBA, BAwL, RDM, RDMSWTN). See rfun_port_plan.md for the design and
// R/model_LBA.R (.lba_rfun, rBAwL) / R/model_RDM.R (rRDM, rWald, rSWTN, rRDMSWTN)
// for the reference semantics these kernels replicate distributionally.
//
// This is the sole translation unit including model_rng.h's RNG primitives
// combined with Rcpp container types; kept separate from particle_ll.cpp's
// TU (which owns model_LBA.h/model_RDM.h) to avoid duplicate-symbol issues
// with those headers' non-inline [[Rcpp::export]] functions.

#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <vector>
#include "model_rng.h"
#include "drift_factor.h"
// FRQ's (h, tau) -> (p, lambda) inversion.  Included rather than duplicated so
// the simulator and the likelihood cannot describe different models; the header
// is export-free precisely so it can be shared across translation units.
#include "model_FRQ.h"

using namespace Rcpp;

// Implemented in model_RDM.h's owning translation unit (particle_ll.cpp).
// The correlated simulator uses the exact same full RDMSWTN marginal CDF as
// the likelihood and inverts it numerically.
double prdmswtn(double t, double mu_drift, double b, double A,
                double s, double t0, double sv,
                double lambda_g, double lambda_k,
                int n_gauss_nodes, bool log_out,
                int kill_shape, bool guess, bool posdrift,
                double erlang_omega);
double rdmswtn_tt_qinv(double u, double tau);

namespace {

std::unordered_map<std::string, int> col_index_map(const Rcpp::NumericMatrix& pars) {
  Rcpp::CharacterVector cn = Rcpp::colnames(pars);
  std::unordered_map<std::string, int> m;
  for (int j = 0; j < cn.size(); j++) m[Rcpp::as<std::string>(cn[j])] = j;
  return m;
}

}  // namespace

// Shared-capacity LogicalRules simulator.  This is the C++ counterpart of
// .lr_capacity_finish_times() in R/make_data.R.  The shared additive drift
// shift loads only on the A/B target pair in an AB stimulus condition; every
// other accumulator keeps its ordinary independent LBA drift draw.
//
// The return value is a matrix of accumulator finishing times, with one row
// per trial and one column per lR level.  Logical-rule assembly is deliberately
// left to the R caller so the ordinary and detection rule paths retain the
// same response/tie/time-racer semantics.
// [[Rcpp::export]]
Rcpp::NumericMatrix logicalrules_capacity_finish_cpp(
    Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
    Rcpp::CharacterVector stimulus, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("logicalrules_capacity_finish_cpp: invalid accumulator/parameter dimensions.");
  }
  const int n_trials = n_rows / n_acc;
  if (stimulus.size() != n_trials) {
    Rcpp::stop("logicalrules_capacity_finish_cpp: stimulus must have one value per trial.");
  }

  const auto ci = col_index_map(pars);
  const int iv = ci.at("v");
  const int isv = ci.at("sv");
  const int ib = ci.at("b");
  const int iA = ci.at("A");
  const int it0 = ci.at("t0");
  const int ikappa = ci.at("kappa");
  const int itau = ci.at("tau");

  int a_role = -1, b_role = -1;
  for (int a = 0; a < n_acc; ++a) {
    const std::string role = Rcpp::as<std::string>(lR_levels[a]);
    if (role == "A") a_role = a;
    if (role == "B") b_role = a;
  }
  if (a_role < 0 || b_role < 0) {
    Rcpp::stop("logicalrules_capacity_finish_cpp: lR_levels must contain A and B.");
  }

  Rcpp::NumericMatrix out(n_trials, n_acc);
  Rcpp::colnames(out) = lR_levels;
  const int max_iter = 10000;

  auto normalized_condition = [](std::string condition) {
    if (condition == "none") condition = "NN";
    else if (condition == "A") condition = "AN";
    else if (condition == "B") condition = "NB";
    else if (condition == "BA" || condition == "A+B" || condition == "B+A") condition = "AB";
    return condition;
  };

  auto draw_finish_time = [&](int row, double drift) {
    const double dt = (pars(row, ib) - pars(row, iA) * R::unif_rand()) / drift;
    if (!R_FINITE(dt) || dt < 0.0) return R_PosInf;
    return pars(row, it0) + dt;
  };

  for (int tr = 0; tr < n_trials; ++tr) {
    const int base = tr * n_acc;
    const double kappa = pars(base + a_role, ikappa);
    const double tau = pars(base + a_role, itau);
    const int b_row = base + b_role;
    if (pars(b_row, ikappa) != kappa || pars(b_row, itau) != tau) {
      Rcpp::stop("LogicalRules capacity requires kappa and tau shared by the A and B rows within each trial.");
    }
    if (!R_FINITE(kappa) || !R_FINITE(tau) || tau < 0.0) {
      Rcpp::stop("LogicalRules capacity requires finite kappa and tau >= 0.");
    }

    const std::string condition = normalized_condition(
      Rcpp::as<std::string>(stimulus[tr]));
    const bool pair_active = condition == "AB" && (kappa != 0.0 || tau != 0.0);

    std::vector<double> drifts(static_cast<size_t>(n_acc), NA_REAL);
    if (pair_active) {
      for (int a = 0; a < n_acc; ++a) {
        if (a == a_role || a == b_role) continue;
        const int row = base + a;
        drifts[static_cast<size_t>(a)] = rtnorm_lower_r(
          pars(row, iv), pars(row, isv), posdrift ? 0.0 : R_NegInf);
      }

      bool accepted = false;
      for (int iter = 0; iter < max_iter; ++iter) {
        const double z = R::norm_rand();
        const double shift = kappa + tau * z;
        const int row_a = base + a_role;
        const int row_b = base + b_role;
        const double drift_a = R::rnorm(pars(row_a, iv) + shift, pars(row_a, isv));
        const double drift_b = R::rnorm(pars(row_b, iv) + shift, pars(row_b, isv));
        drifts[static_cast<size_t>(a_role)] = drift_a;
        drifts[static_cast<size_t>(b_role)] = drift_b;
        if (!posdrift || (drift_a > 0.0 && drift_b > 0.0)) {
          accepted = true;
          break;
        }
      }
      if (!accepted) {
        Rcpp::stop("LogicalRules capacity jointly positive drift rejection exceeded %d sweeps; check that the target drift means are not far below zero.", max_iter);
      }
    } else {
      for (int a = 0; a < n_acc; ++a) {
        const int row = base + a;
        drifts[static_cast<size_t>(a)] = rtnorm_lower_r(
          pars(row, iv), pars(row, isv), posdrift ? 0.0 : R_NegInf);
      }
    }

    for (int a = 0; a < n_acc; ++a) {
      out(tr, a) = draw_finish_time(base + a, drifts[static_cast<size_t>(a)]);
    }
  }
  return out;
}

namespace {

struct RaceOut {
  std::vector<int> R;        // 1-based winner accumulator, 0 = NA
  std::vector<double> rt;    // meaningful only where R != 0
  std::vector<int> best_i;   // 0-based winner accumulator index, -1 = NA (pre time-resample)
  std::vector<int> omitted;  // 1 = omission (NA/Inf), 0 = skipped/invalid (NA/NA) or valid winner
};

// Argmin-per-trial race resolution shared by all four families.
//   dt      : n_acc x n_trials, column-major (matches pars row order)
//   t0col   : per-row t0 to add to the winner's finish time, or nullptr if
//             dt already has t0 baked in (BAwL/RDMSWTN erlang-clock path)
//   ok_row  : if non-null, trial forced to NA when ok_row[first row] == 0
//             (matches .lba_rfun/rBAwL's `ok <- matrix(ok,nrow)[1,]` reduction;
//             rRDM/rRDMSWTN have no such reduction, pass nullptr)
RaceOut resolve_race(const std::vector<double>& dt, int n_acc, int n_trials,
                     const std::vector<double>* t0col, const std::vector<int>* ok_row) {
  RaceOut out;
  out.R.assign(n_trials, 0);
  out.rt.assign(n_trials, NA_REAL);
  out.best_i.assign(n_trials, -1);
  out.omitted.assign(n_trials, 0);
  for (int tr = 0; tr < n_trials; tr++) {
    const int base = tr * n_acc;
    if (ok_row && !(*ok_row)[base]) continue;
    double best = R_PosInf;
    int best_a = -1;
    for (int a = 0; a < n_acc; a++) {
      const double v = dt[base + a];
      if (v < best) { best = v; best_a = a; }
    }
    if (best_a < 0) {
      out.omitted[tr] = 1;  // all-Inf column -> omission coded as NA/Inf
      continue;
    }
    out.R[tr] = best_a + 1;
    out.best_i[tr] = best_a;
    out.rt[tr] = t0col ? ((*t0col)[base + best_a] + best) : best;
  }
  return out;
}

// Resample "time"-level winners uniformly from non-time/non-nogo levels
// (matches .apply_timed_guess_winner, R/utils.R:53-72). Mutates R in place.
// Returns whether an isTime column should be included (i.e. "time" is a
// level AND at least one trial's winner was "time"); isTime is always
// boolean (never NA) on this path.
bool resolve_time_level(std::vector<int>& R, std::vector<int>& isTime_out,
                        const Rcpp::CharacterVector& lR_levels) {
  const int n_acc = lR_levels.size();
  int time_code = -1, nogo_code = -1;
  for (int i = 0; i < n_acc; i++) {
    const std::string lv = Rcpp::as<std::string>(lR_levels[i]);
    if (lv == "time") time_code = i + 1;
    if (lv == "nogo") nogo_code = i + 1;
  }
  if (time_code < 0) return false;

  const int n = (int)R.size();
  bool any_time = false;
  for (int i = 0; i < n; i++) if (R[i] == time_code) { any_time = true; break; }
  if (!any_time) return false;

  std::vector<int> guess_levels;
  for (int i = 1; i <= n_acc; i++) if (i != time_code && i != nogo_code) guess_levels.push_back(i);

  isTime_out.assign(n, 0);
  for (int i = 0; i < n; i++) if (R[i] == time_code) isTime_out[i] = 1;
  if (guess_levels.empty()) return true;  // R left as "time" level, per R semantics

  for (int i = 0; i < n; i++) {
    if (R[i] == time_code) {
      int idx = (int)(R::unif_rand() * guess_levels.size());
      if (idx >= (int)guess_levels.size()) idx = (int)guess_levels.size() - 1;
      R[i] = guess_levels[idx];
    }
  }
  return true;
}

Rcpp::List pack_result(const std::vector<int>& R, const std::vector<double>& rt,
                       const std::vector<int>& omitted,
                       bool has_isTime, const std::vector<int>& isTime) {
  const int n = (int)R.size();
  Rcpp::IntegerVector Rv(n);
  Rcpp::NumericVector rtv(n);
  for (int i = 0; i < n; i++) {
    if (R[i] == 0) {
      Rv[i] = NA_INTEGER;
      rtv[i] = omitted[i] ? R_PosInf : NA_REAL;
    } else {
      Rv[i] = R[i];
      rtv[i] = rt[i];
    }
  }
  if (!has_isTime) {
    return Rcpp::List::create(Rcpp::Named("R") = Rv, Rcpp::Named("rt") = rtv,
                              Rcpp::Named("isTime") = R_NilValue);
  }
  Rcpp::LogicalVector itv(n);
  for (int i = 0; i < n; i++) itv[i] = (isTime[i] < 0) ? NA_LOGICAL : (Rboolean)(isTime[i] == 1);
  return Rcpp::List::create(Rcpp::Named("R") = Rv, Rcpp::Named("rt") = rtv,
                            Rcpp::Named("isTime") = itv);
}

}  // namespace

// pars columns: v, sv, b, A, t0. lR_levels/ok/pars rows are trial-major,
// accumulator-minor (accumulators contiguous within a trial) -- matches
// R's .lba_rfun (model_LBA.R:35-66).
// [[Rcpp::export]]
Rcpp::List rlba_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                    Rcpp::LogicalVector ok, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), isv = ci.at("sv"), ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; r++) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double lo = posdrift ? 0.0 : R_NegInf;
    const double drift = rtnorm_lower_r(pars(r, iv), pars(r, isv), lo);
    const double u = R::unif_rand();
    const double d = (pars(r, ib) - pars(r, iA) * u) / drift;
    dt[r] = (d < 0.0) ? R_PosInf : d;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// pars columns: v, B, A, t0 (+ optional s). RDM is always posdrift = TRUE
// (rRDM never exposes a posdrift toggle -- model_RDM.R:101-134).
// [[Rcpp::export]]
Rcpp::List rrdm_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                    Rcpp::LogicalVector ok) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), iB = ci.at("B"), iA = ci.at("A"), it0 = ci.at("t0");
  const int is_ = ci.count("s") ? ci.at("s") : -1;

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  for (int r = 0; r < n_rows; r++) {
    t0col[r] = pars(r, it0);
    if (!ok[r]) continue;
    double v = pars(r, iv), B = pars(r, iB), A = pars(r, iA);
    if (is_ >= 0) {
      const double s = pars(r, is_);
      v /= s; B /= s; A /= s;
    }
    if (B < 0.0) B = 0.0;
    if (A < 0.0) A = 0.0;
    dt[r] = rwald_acc_r(B, v, A, 1.0, true);
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, nullptr);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// Shared BAwL simulator.  The optional drift_override is supplied by the
// correlated-factor entry point; all clock, leak, and race handling is shared.
static Rcpp::List rbawl_cpp_impl(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                                 Rcpp::LogicalVector ok, bool posdrift, int erlang,
                                 bool guess, bool global,
                                 const std::vector<double>* drift_override,
                                 int launch = 0) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  // The lognormal launch strength keeps BAwD's naming convention: (mu, sigma)
  // occupy the (v, sv) slots and every downstream step is identical.
  const bool logn = (launch == 1);
  const int iv = ci.at(logn ? "mu" : "v"), isv = ci.at(logn ? "sigma" : "sv"),
            ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ik = ci.at("k"), ilg = ci.at("lambda_g"), ilk = ci.at("lambda_k");
  const int iomega = ci.count("omega") ? ci.at("omega") : -1;
  const double eps = 1e-10;

  auto omega_at = [&](int r) -> double {
    if (erlang <= 1) return 1.0;
    if (erlang == 2) return 0.0;
    if (iomega < 0) Rcpp::stop("mixed BAwL Erlang mode requires parameter column 'omega'.");
    return pars(r, iomega);
  };

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; r++) ok_row[r] = ok[r] ? 1 : 0;

  std::vector<double> tk_global(n_trials, R_PosInf);
  if (global) {
    for (int tr = 0; tr < n_trials; tr++) {
      const int base = tr * n_acc;
      const double lam0 = pars(base, ilk);
      for (int a = 1; a < n_acc; a++) {
        if (pars(base + a, ilk) != lam0)
          Rcpp::stop("global=TRUE requires lambda_k to be constant across accumulators");
      }
      if (!ISNAN(lam0) && lam0 > 0.0) tk_global[tr] = rerlang_clock_r(lam0, erlang, omega_at(base));
    }
  }

  // Core EAM step (leaky ballistic accumulator hit-time inversion), raw-time axis.
  for (int r = 0; r < n_rows; r++) {
    if (!ok[r]) continue;
    const double lo = posdrift ? 0.0 : R_NegInf;
    const double drift = drift_override
      ? (*drift_override)[static_cast<size_t>(r)]
      : (logn ? R::rlnorm(pars(r, iv), pars(r, isv))
              : rtnorm_lower_r(pars(r, iv), pars(r, isv), lo));
    // A lognormal launch strength is positive by construction.
    if (!logn && posdrift && !(drift > 0.0)) {
      dt[r] = R_PosInf;
      continue;
    }
    const double k = pars(r, ik);
    double d;
    if (k < eps) {
      const double u = R::unif_rand();
      d = (pars(r, ib) - pars(r, iA) * u) / drift;
    } else if (drift > k * pars(r, ib)) {
      const double u = R::unif_rand();
      const double num = drift - k * pars(r, ib);
      const double den = drift - k * pars(r, iA) * u;
      double ratio = num / den;
      ratio = std::fmin(std::fmax(ratio, std::numeric_limits<double>::min()), 1.0 - 1e-15);
      d = (-1.0 / k) * std::log(ratio);
    } else {
      d = R_PosInf;
    }
    if (d < 0.0) d = R_PosInf;
    dt[r] = d + pars(r, it0);
  }

  if (guess || !global) {
    std::vector<double> tg_local(n_rows, R_PosInf), tk_local(n_rows, R_PosInf);
    if (!global) {
      for (int r = 0; r < n_rows; r++) {
        if (!ok[r]) continue;
        const double lam = pars(r, ilk);
        if (lam > 0.0) tk_local[r] = rerlang_clock_r(lam, erlang, omega_at(r));
      }
    }
    if (guess) {
      for (int tr = 0; tr < n_trials; tr++) {
        const int base = tr * n_acc;
        for (int a = 0; a < n_acc; a++) {
          const int r = base + a;
          if (!ok[r]) continue;
          if (Rcpp::as<std::string>(lR_levels[a]) == "nogo") continue;
          const double lam = pars(r, ilg);
          if (lam > 0.0) tg_local[r] = rerlang_clock_r(lam, erlang, omega_at(r));
        }
      }
    }
    for (int r = 0; r < n_rows; r++) {
      if (guess && !global) {
        const double cand = std::fmin(dt[r], tg_local[r]);
        dt[r] = (cand < tk_local[r]) ? cand : R_PosInf;
      } else if (guess) {
        dt[r] = std::fmin(dt[r], tg_local[r]);
      } else {
        dt[r] = (dt[r] < tk_local[r]) ? dt[r] : R_PosInf;
      }
    }
  }

  if (global) {
    for (int tr = 0; tr < n_trials; tr++) {
      const int base = tr * n_acc;
      double m = R_PosInf;
      for (int a = 0; a < n_acc; a++) m = std::fmin(m, dt[base + a]);
      if (tk_global[tr] < m) for (int a = 0; a < n_acc; a++) dt[base + a] = R_PosInf;
    }
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// pars columns: v, sv, b, A, t0, k, lambda_g, lambda_k (+ optional omega), or
// mu, sigma in place of v, sv when launch = 1 (lognormal).
// Matches R's rBAwL (model_LBA.R:347-479).
// [[Rcpp::export]]
Rcpp::List rbawl_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                     Rcpp::LogicalVector ok, bool posdrift, int erlang, bool guess,
                     bool global, int launch = 0) {
  return rbawl_cpp_impl(pars, lR_levels, ok, posdrift, erlang, guess, global, nullptr,
                        launch);
}

// Correlated BAwL simulator.  The low-level entry point receives signed
// row-level factor-variance shares; BAwLcorr's R Ttransform derives those from
// one cell-level rho when accumulator roles are available.  rho == 0 gives an
// accumulator-specific draw only, so arbitrary races and independent PM/false-
// alarm racers remain supported.
// [[Rcpp::export]]
Rcpp::List rbawl_corr_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                          Rcpp::LogicalVector ok, bool posdrift, int erlang,
                          bool guess, bool global) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), isv = ci.at("sv"), irho = ci.at("rho");
  if (n_acc <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rbawl_corr_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) Rcpp::stop("rbawl_corr_cpp: ok has the wrong length.");

  std::vector<double> drifts(static_cast<size_t>(n_rows), R_PosInf);
  drift_factor_draw_correlated(pars, iv, isv, irho, ok, n_acc, posdrift,
                               drifts, "rbawl_corr_cpp");

  return rbawl_cpp_impl(pars, lR_levels, ok, posdrift, erlang,
                        guess, global, &drifts);
}

// BAwD simulator.  pars columns: b, A, t0, k, ell plus the launch pair, named
// (v, sv) for launch = 0 (truncated normal) and (mu, sigma) for launch = 1
// (lognormal).  `launch` must be the value BAwD() derived from its
// drift_distribution argument -- the same one that reaches the likelihood as
// ctx->bawd_launch -- or simulation and estimation describe different models.
// Matches R's rBAwD (R/model_BAwD.R).
// [[Rcpp::export]]
Rcpp::List rbawd_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                     Rcpp::LogicalVector ok, int launch, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbawd_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rbawd_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == 1);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawd_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawd_cpp: the normal launch requires columns 'v' and 'sv'.");
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ik = ci.at("k"), iell = ci.at("ell");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = logn
      ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
      : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf);
    const double z = pars(r, iA) * R::unif_rand();
    const double u = bawd_hit_time_r(V, pars(r, ib) - z, pars(r, ik),
                                     pars(r, iell));
    dt[r] = (R_FINITE(u) && u >= 0.0) ? u : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// BAwDp simulator.  The launch pair is (v, sv) for launch = 0 and
// (mu, sigma) for launch = 1.  The internal clock is monotone until its
// universal freeze time; the sampled LBA crossing is inverted on that clock.
// [[Rcpp::export]]
Rcpp::List rbawdp_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                      Rcpp::LogicalVector ok, int launch, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbawdp_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rbawdp_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == 1);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawdp_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawdp_cpp: the normal launch requires columns 'v' and 'sv'.");
  for (const char* nm : {"b", "A", "t0", "k", "lambda"})
    if (!ci.count(nm)) Rcpp::stop("rbawdp_cpp: missing parameter column '%s'.", nm);
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0");
  const int ik = ci.at("k"), ilambda = ci.at("lambda");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = logn
      ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
      : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf);
    const double z = pars(r, iA) * R::unif_rand();
    const double u = bawdp_hit_time_r(V, pars(r, ib) - z, pars(r, ik),
                                      pars(r, ilambda));
    dt[r] = (R_FINITE(u) && u >= 0.0) ? u : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// FRQ simulator.  pars columns: alpha, beta, h, tau, t0, delta.
//
// Exact, not approximate: the process is distributionally identical to drawing
// the latent quorum U ~ Beta(alpha, beta) and inverting the registration CDF at
// U/p.  That follows from F(x) = I_{q(x)}(alpha, beta) with q(x) = p G(x): the
// decision time is the quantile transform of a Beta(alpha, beta) draw, so
// U <= q(x) <=> T <= x, i.e. T = G^{-1}(U/p).  U > p means the reservoir
// saturates below the required quorum, so the accumulator NEVER terminates and
// the trial contributes +Inf -- the package's omission convention.  For integer
// alpha = K and beta = N - K + 1 this reproduces explicitly building N cues and
// waiting for the K-th.  Matches R's rFRQ (R/model_FRQ.R).
//
// With threshold variability (delta > 0) the latent quorum percentile is no
// longer Beta: F(x) = H(I_{q(x)}(alpha, beta)), so the quantile transform runs
// through H as well.  Draw V ~ U(0,1) and set U = qbeta(H^{-1}(V), alpha, beta);
// then U <= q(x) <=> V <= F(x) exactly as before, and the omission test U > p is
// unchanged because it is equivalent to V > h.  delta == 0 keeps the direct
// rbeta draw rather than routing through qbeta(runif()), which is the same
// distribution but a different RNG stream -- seeded FRQ simulations predating
// delta must keep reproducing.
// [[Rcpp::export]]
Rcpp::List rfrq_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                    Rcpp::LogicalVector ok) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rfrq_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rfrq_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  for (const char* nm : {"alpha", "beta", "h", "tau", "t0"})
    if (!ci.count(nm)) Rcpp::stop("rfrq_cpp: missing parameter column '%s'.", nm);
  const int ia = ci.at("alpha"), ib = ci.at("beta"), ih = ci.at("h"),
            it = ci.at("tau"), it0 = ci.at("t0");
  // Absent delta means the documented default of no threshold variability.
  // Safe to default here, unlike on the likelihood path: this lookup is BY NAME,
  // so a missing column cannot silently resolve to a different parameter.
  const int idl = ci.count("delta") ? ci.at("delta") : -1;

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const FrqPars s = frq_derive(pars(r, ia), pars(r, ib), pars(r, ih),
                                 pars(r, it), idl >= 0 ? pars(r, idl) : 0.0);
    if (!s.ok) continue;                       // invalid row: never finishes
    const double u = s.hh.active
      ? R::qbeta(frq_h_inv(R::unif_rand(), s.hh), s.alpha, s.beta, 1, 0)
      : R::rbeta(s.alpha, s.beta);
    if (!(u <= s.p)) continue;                 // quorum unreachable: omission
    const double d = -std::log1p(-u / s.p) / s.lambda;
    dt[r] = (R_FINITE(d) && d >= 0.0) ? d : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// pars columns: v, b, A, t0, sv, lambda_g, lambda_k (+ optional s, omega).
// erlang_type: "none" | "local_kill" | "global_kill" | "local_guess" | "local_kill_guess".
// Matches R's rRDMSWTN (model_RDM.R:872-1014) and rSWTN (model_RDM.R:489-532).
//
// `drift_override` supplies a pre-drawn rate per row, which is how the
// correlated-draw simulator injects its equicorrelated draws (drift_factor.h)
// without duplicating the Erlang clock, start-point and race machinery below.
// A null pointer keeps the ordinary independent per-row draw.
static Rcpp::List rrdmswtn_cpp_impl(Rcpp::NumericMatrix pars,
                                    Rcpp::CharacterVector lR_levels,
                                    Rcpp::LogicalVector ok, int erlang_shape,
                                    std::string erlang_type, bool posdrift,
                                    const std::vector<double>* drift_override) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"), isv = ci.at("sv"),
            ilg = ci.at("lambda_g"), ilk = ci.at("lambda_k");
  const int is_ = ci.count("s") ? ci.at("s") : -1;
  const int iomega = ci.count("omega") ? ci.at("omega") : -1;

  auto omega_at = [&](int r) -> double {
    if (erlang_shape <= 1) return 1.0;
    if (erlang_shape == 2) return 0.0;
    if (iomega < 0) Rcpp::stop("mixed RDMSWTN Erlang mode requires parameter column 'omega'.");
    return pars(r, iomega);
  };

  const bool guess = (erlang_type == "local_guess" || erlang_type == "local_kill_guess");
  const bool global = (erlang_type == "global_kill");
  const bool local_kill = (erlang_type == "local_kill" || erlang_type == "local_kill_guess");

  std::vector<double> dt(n_rows, R_PosInf);

  std::vector<double> tk_global(n_trials, R_PosInf);
  if (global) {
    for (int tr = 0; tr < n_trials; tr++) {
      const int base = tr * n_acc;
      const double lam0 = pars(base, ilk);
      for (int a = 1; a < n_acc; a++) {
        if (pars(base + a, ilk) != lam0)
          Rcpp::stop("global_kill requires lambda_k to be constant across accumulators");
      }
      // Global kill uses the plain (non-mixture) rgamma draw, matching R's
      // rRDMSWTN global branch (model_RDM.R:917) -- an existing asymmetry
      // vs. the local guess/kill draws, which do use the shape-3 mixture.
      if (!ISNAN(lam0) && lam0 > 0.0) tk_global[tr] = R::rgamma((double)erlang_shape, 1.0 / lam0);
    }
  }

  for (int r = 0; r < n_rows; r++) {
    if (!ok[r]) continue;
    double b = pars(r, ib); if (b < 0.0) b = 0.0;
    double A = pars(r, iA); if (A < 0.0) A = 0.0;
    const double s = (is_ >= 0) ? pars(r, is_) : 1.0;
    const double v = pars(r, iv), sv = pars(r, isv);
    double v_draw = v;
    if (drift_override != nullptr) {
      v_draw = (*drift_override)[static_cast<size_t>(r)];
    } else if (R_FINITE(sv) && sv > 1e-12) {
      if (posdrift) {
        const double lo = R::pnorm(0.0, v, sv, 1, 0);
        double u = lo + R::unif_rand() * (1.0 - lo);
        u = std::fmin(std::fmax(u, 1e-300), 1.0 - 1e-16);
        v_draw = R::qnorm(u, v, sv, 1, 0);
      } else {
        v_draw = R::rnorm(v, sv);
      }
    }
    const double d = rwald_acc_r(b - A, v_draw, A, s, posdrift);
    dt[r] = d + pars(r, it0);
  }

  std::vector<int> guess_win(n_rows, 0);
  if (guess || local_kill) {
    std::vector<double> tg_local(n_rows, R_PosInf), tk_local(n_rows, R_PosInf);
    if (guess) {
      for (int tr = 0; tr < n_trials; tr++) {
        const int base = tr * n_acc;
        for (int a = 0; a < n_acc; a++) {
          const int r = base + a;
          if (!ok[r]) continue;
          if (Rcpp::as<std::string>(lR_levels[a]) == "nogo") continue;
          const double lam = pars(r, ilg);
          if (lam > 0.0) tg_local[r] = rerlang_clock_r(lam, erlang_shape, omega_at(r));
        }
      }
    }
    if (local_kill) {
      for (int r = 0; r < n_rows; r++) {
        if (!ok[r]) continue;
        const double lam = pars(r, ilk);
        if (lam > 0.0) tk_local[r] = rerlang_clock_r(lam, erlang_shape, omega_at(r));
      }
    }
    for (int r = 0; r < n_rows; r++) {
      if (guess && local_kill) {
        const bool win = (tg_local[r] < dt[r]) && (tg_local[r] < tk_local[r]);
        guess_win[r] = win ? 1 : 0;
        const double cand = std::fmin(dt[r], tg_local[r]);
        dt[r] = (cand < tk_local[r]) ? cand : R_PosInf;
      } else if (guess) {
        guess_win[r] = (tg_local[r] < dt[r]) ? 1 : 0;
        dt[r] = std::fmin(dt[r], tg_local[r]);
      } else {
        dt[r] = (dt[r] < tk_local[r]) ? dt[r] : R_PosInf;
      }
    }
  }

  if (global) {
    for (int tr = 0; tr < n_trials; tr++) {
      const int base = tr * n_acc;
      double m = R_PosInf;
      for (int a = 0; a < n_acc; a++) m = std::fmin(m, dt[base + a]);
      if (tk_global[tr] < m) for (int a = 0; a < n_acc; a++) dt[base + a] = R_PosInf;
    }
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, nullptr);
  std::vector<int> isTime;
  bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);

  // guess overwrites isTime entirely with "did the guess clock decide the
  // winner" (indexed by the ORIGINAL winner, pre time-resample), discarding
  // any time-level semantics -- matches R/model_RDM.R:1009-1012 exactly.
  if (guess) {
    isTime.assign(n_trials, -1);  // NA sentinel for bad (all-Inf) trials
    for (int tr = 0; tr < n_trials; tr++) {
      if (res.best_i[tr] < 0) continue;
      isTime[tr] = guess_win[tr * n_acc + res.best_i[tr]];
    }
    has_isTime = true;
  }

  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// [[Rcpp::export]]
Rcpp::List rrdmswtn_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                        Rcpp::LogicalVector ok, int erlang_shape,
                        std::string erlang_type, bool posdrift) {
  return rrdmswtn_cpp_impl(pars, lR_levels, ok, erlang_shape, erlang_type,
                           posdrift, nullptr);
}

// Correlated *drift draws*: one shared latent factor per trial couples the
// rates, then every accumulator races independently on its own draw.  This is
// the same construction BAwLcorr simulates (rbawl_corr_cpp); the only
// difference is the accumulator the drawn rate is fed to.
// [[Rcpp::export]]
Rcpp::List rrdmswtn_drift_corr_cpp(Rcpp::NumericMatrix pars,
                                   Rcpp::CharacterVector lR_levels,
                                   Rcpp::LogicalVector ok, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rrdmswtn_drift_corr_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) {
    Rcpp::stop("rrdmswtn_drift_corr_cpp: ok has the wrong length.");
  }
  const auto ci = col_index_map(pars);
  std::vector<double> drifts(static_cast<size_t>(n_rows), R_PosInf);
  drift_factor_draw_correlated(pars, ci.at("v"), ci.at("sv"), ci.at("rho"), ok,
                               n_acc, posdrift, drifts,
                               "rrdmswtn_drift_corr_cpp");
  return rrdmswtn_cpp_impl(pars, lR_levels, ok, 1, "none", posdrift, &drifts);
}

namespace {

double rdmswtn_quantile_cpp(double u, double v, double b, double A, double s,
                            double t0, double sv, bool posdrift = true) {
  u = std::fmin(std::nextafter(1.0, 0.0),
                std::fmax(std::numeric_limits<double>::min(), u));
  auto cdf = [&](double t) {
    return prdmswtn(t, v, b, A, s, t0, sv, 0.0, 0.0, 20, false,
                    1, false, posdrift, 1.0);
  };
  if (!posdrift) {
    // Defective marginal: F saturates at the hit probability p = F(Inf) < 1.
    // Copula uniforms at or above that plateau are the atom at infinity, i.e.
    // this accumulator never finishes.  Without this the bracketing loop below
    // would run to its ceiling and abort on a perfectly valid draw.
    const double cap = cdf(R_PosInf);
    if (!R_FINITE(cap)) {
      Rcpp::stop("rrdmswtn_corr_cpp: non-finite marginal hit probability while inverting a finishing-time quantile.");
    }
    if (u >= cap) return R_PosInf;
  }
  double lo = 0.0;
  double hi = std::fmax(1.0, t0 + 1.0);
  double fhi = cdf(hi);
  for (int iter = 0; iter < 80 && (!R_FINITE(fhi) || fhi < u); ++iter) {
    hi *= 2.0;
    if (!R_FINITE(hi) || hi > 1e12) {
      Rcpp::stop("rrdmswtn_corr_cpp: could not bracket a finishing-time quantile; check the marginal parameters.");
    }
    fhi = cdf(hi);
  }
  if (!R_FINITE(fhi) || fhi < u) {
    Rcpp::stop("rrdmswtn_corr_cpp: could not bracket a finishing-time quantile; check the marginal parameters.");
  }
  for (int iter = 0; iter < 100; ++iter) {
    const double mid = lo + 0.5 * (hi - lo);
    const double fm = cdf(mid);
    if (!R_FINITE(fm)) {
      Rcpp::stop("rrdmswtn_corr_cpp: non-finite marginal CDF while inverting a finishing-time quantile.");
    }
    if (fm < u) lo = mid;
    else hi = mid;
    if (hi - lo <= 1e-10 * std::fmax(1.0, hi)) break;
  }
  return lo + 0.5 * (hi - lo);
}

}  // namespace

// Gaussian-copula RDMSWTN simulator. rho is a direct natural-scale pair
// correlation: exactly two active rows may share one signed nonzero value;
// every other active row is independent.
// [[Rcpp::export]]
Rcpp::List rrdmswtn_corr_cpp(Rcpp::NumericMatrix pars,
                             Rcpp::CharacterVector lR_levels,
                             Rcpp::LogicalVector ok,
                             bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rrdmswtn_corr_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) {
    Rcpp::stop("rrdmswtn_corr_cpp: ok has the wrong length.");
  }
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), ib = ci.at("b"), iA = ci.at("A");
  const int it0 = ci.at("t0"), isv = ci.at("sv"), irho = ci.at("rho");
  const int is_ = ci.count("s") ? ci.at("s") : -1;
  const int ilg = ci.at("lambda_g"), ilk = ci.at("lambda_k");

  bool any_nonzero = false;
  for (int r = 0; r < n_rows; ++r) {
    if (!ok[r]) continue;
    const double rho = pars(r, irho);
    if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
      Rcpp::stop("rrdmswtn_corr_cpp: rho must be finite and lie in [-1, 1].");
    }
    if (pars(r, ilg) != 0.0 || pars(r, ilk) != 0.0) {
      Rcpp::stop("rrdmswtn_corr_cpp: guess and kill clocks are not supported.");
    }
    if (std::fabs(rho) <= 1e-12) continue;
    any_nonzero = true;
    // See .check_rdmswtn_corr_io_sv() in R/model_RDM.R: the copula couples
    // defective finishing times, but sv > 0 under unrestricted drifts is
    // reserved for a correlated-drift model.
    if (!posdrift && !(pars(r, isv) <= 1e-12)) {
      Rcpp::stop("rrdmswtn_corr_cpp: posdrift = FALSE requires sv = 0 on the correlated rows.");
    }
  }
  if (!any_nonzero) {
    return rrdmswtn_cpp(pars, lR_levels, ok, 1, "none", posdrift);
  }

  std::vector<double> u(static_cast<size_t>(n_rows), 0.5);
  for (int r = 0; r < n_rows; ++r) {
    if (ok[r]) u[static_cast<size_t>(r)] = R::unif_rand();
  }
  for (int tr = 0; tr < n_trials; ++tr) {
    const int start = tr * n_acc;
    int pair[2] = {-1, -1};
    int n_pair = 0;
    double pair_rho = 0.0;
    for (int a = 0; a < n_acc; ++a) {
      const int r = start + a;
      if (!ok[r]) continue;
      const double rho = pars(r, irho);
      if (std::fabs(rho) <= 1e-12) continue;
      if (n_pair < 2) pair[n_pair] = r;
      ++n_pair;
      if (n_pair == 1) pair_rho = rho;
      else if (std::fabs(rho - pair_rho) > 1e-12) {
        Rcpp::stop("rrdmswtn_corr_cpp: the two participating rows must have the same signed nonzero rho.");
      }
    }
    if (n_pair > 2) {
      Rcpp::stop("rrdmswtn_corr_cpp: at most two active rows may have nonzero rho in a trial.");
    }
    if (n_pair == 2) {
      const double z1 = R::norm_rand();
      const double z2 = pair_rho * z1 +
        std::sqrt(std::fmax(0.0, 1.0 - pair_rho * pair_rho)) * R::norm_rand();
      u[static_cast<size_t>(pair[0])] = R::pnorm(z1, 0.0, 1.0, 1, 0);
      u[static_cast<size_t>(pair[1])] = R::pnorm(z2, 0.0, 1.0, 1, 0);
    }
  }

  std::vector<double> dt(static_cast<size_t>(n_rows), R_PosInf);
  for (int r = 0; r < n_rows; ++r) {
    if (!ok[r]) continue;
    const double s = (is_ >= 0) ? pars(r, is_) : 1.0;
    dt[static_cast<size_t>(r)] = rdmswtn_quantile_cpp(
      u[static_cast<size_t>(r)], pars(r, iv), pars(r, ib), pars(r, iA),
      s, pars(r, it0), pars(r, isv), posdrift
    );
  }
  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, nullptr);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// RDMSWTN under the linear exhaustion clock. Draw the ordinary operational
// finishing time, omit when it exceeds Q=tau/2, and otherwise invert q.
// See rrdmswtn_cpp_impl() for the drift_override contract.
static Rcpp::List rrdmswtn_tt_cpp_impl(Rcpp::NumericMatrix pars,
                                       Rcpp::CharacterVector lR_levels,
                                       Rcpp::LogicalVector ok, bool posdrift,
                                       const std::vector<double>* drift_override) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rrdmswtn_tt_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) {
    Rcpp::stop("rrdmswtn_tt_cpp: ok has the wrong length.");
  }
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), ib = ci.at("b"), iA = ci.at("A");
  const int it0 = ci.at("t0"), isv = ci.at("sv"), itau = ci.at("tau");
  const int is_ = ci.count("s") ? ci.at("s") : -1;

  std::vector<double> dt(static_cast<size_t>(n_rows), R_PosInf);
  for (int r = 0; r < n_rows; ++r) {
    if (!ok[r]) continue;
    const double tau = pars(r, itau);
    // tau = +Inf is the exact identity-clock limit: the ordinary
    // operational finish is returned on the physical-time axis.  Reject
    // NaN and negative/non-positive values, but keep the analytic limit
    // available to direct simulator callers.
    if (ISNAN(tau) || !(tau > 0.0)) {
      Rcpp::stop("rrdmswtn_tt_cpp: tau must be positive (or +Inf).");
    }
    const double s = (is_ >= 0) ? pars(r, is_) : 1.0;
    const double v = pars(r, iv), sv = pars(r, isv);
    double v_draw = v;
    if (drift_override != nullptr) {
      v_draw = (*drift_override)[static_cast<size_t>(r)];
    } else if (R_FINITE(sv) && sv > 1e-12) {
      v_draw = rtnorm_lower_r(v, sv, posdrift ? 0.0 : R_NegInf);
    }
    double b = std::fmax(0.0, pars(r, ib));
    double A = std::fmax(0.0, pars(r, iA));
    const double operational = rwald_acc_r(b - A, v_draw, A, s, posdrift);
    if (R_FINITE(operational) &&
        (!R_FINITE(tau) || operational <= 0.5 * tau)) {
      dt[static_cast<size_t>(r)] =
        pars(r, it0) + rdmswtn_tt_qinv(operational, tau);
    }
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, nullptr);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// [[Rcpp::export]]
Rcpp::List rrdmswtn_tt_cpp(Rcpp::NumericMatrix pars,
                           Rcpp::CharacterVector lR_levels,
                           Rcpp::LogicalVector ok, bool posdrift) {
  return rrdmswtn_tt_cpp_impl(pars, lR_levels, ok, posdrift, nullptr);
}

// Correlated drift draws under the exhaustion clock; see
// rrdmswtn_drift_corr_cpp().
// [[Rcpp::export]]
Rcpp::List rrdmswtn_tt_drift_corr_cpp(Rcpp::NumericMatrix pars,
                                      Rcpp::CharacterVector lR_levels,
                                      Rcpp::LogicalVector ok, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rrdmswtn_tt_drift_corr_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) {
    Rcpp::stop("rrdmswtn_tt_drift_corr_cpp: ok has the wrong length.");
  }
  const auto ci = col_index_map(pars);
  std::vector<double> drifts(static_cast<size_t>(n_rows), R_PosInf);
  drift_factor_draw_correlated(pars, ci.at("v"), ci.at("sv"), ci.at("rho"), ok,
                               n_acc, posdrift, drifts,
                               "rrdmswtn_tt_drift_corr_cpp");
  return rrdmswtn_tt_cpp_impl(pars, lR_levels, ok, posdrift, &drifts);
}

namespace {

double rdmswtn_operational_quantile_bounded(
    double u, double FQ, double Q, double v, double b, double A,
    double s, double sv, bool posdrift = true) {
  if (u >= FQ) return Q;
  double lo = 0.0, hi = Q;
  for (int iter = 0; iter < 100; ++iter) {
    const double mid = lo + 0.5 * (hi - lo);
    const double fm = prdmswtn(
        mid, v, b, A, s, 0.0, sv, 0.0, 0.0, 20, false,
        1, false, posdrift, 1.0);
    if (!R_FINITE(fm)) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: non-finite marginal CDF while inverting.");
    }
    if (fm < u) lo = mid;
    else hi = mid;
    if (hi - lo <= 1e-10 * std::fmax(1.0, hi)) break;
  }
  return lo + 0.5 * (hi - lo);
}

}  // namespace

// Gaussian-copula simulator for the exhausted marginals. Copula uniforms above
// the finite marginal plateau represent the atom at infinity.
// [[Rcpp::export]]
Rcpp::List rrdmswtn_tt_corr_cpp(Rcpp::NumericMatrix pars,
                                Rcpp::CharacterVector lR_levels,
                                Rcpp::LogicalVector ok,
                                bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0) {
    Rcpp::stop("rrdmswtn_tt_corr_cpp: invalid accumulator/parameter dimensions.");
  }
  if (ok.size() != n_rows) {
    Rcpp::stop("rrdmswtn_tt_corr_cpp: ok has the wrong length.");
  }
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), ib = ci.at("b"), iA = ci.at("A");
  const int it0 = ci.at("t0"), isv = ci.at("sv"), itau = ci.at("tau");
  const int irho = ci.at("rho");
  const int is_ = ci.count("s") ? ci.at("s") : -1;

  bool any_nonzero = false;
  for (int r = 0; r < n_rows; ++r) {
    if (!ok[r]) continue;
    const double rho = pars(r, irho);
    if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: rho must be finite and lie in [-1, 1].");
    }
    const double tau = pars(r, itau);
    if (ISNAN(tau) || !(tau > 0.0)) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: tau must be positive (or +Inf).");
    }
    if (std::fabs(rho) <= 1e-12) continue;
    any_nonzero = true;
    // sv > 0 with unrestricted drifts is reserved for a correlated-drift
    // model; see .check_rdmswtn_corr_io_sv() in R/model_RDM.R.
    if (!posdrift && !(pars(r, isv) <= 1e-12)) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: posdrift = FALSE requires sv = 0 on the correlated rows.");
    }
  }
  if (!any_nonzero) {
    return rrdmswtn_tt_cpp(pars, lR_levels, ok, posdrift);
  }

  std::vector<double> u(static_cast<size_t>(n_rows), 0.5);
  for (int r = 0; r < n_rows; ++r) {
    if (ok[r]) u[static_cast<size_t>(r)] = R::unif_rand();
  }
  for (int tr = 0; tr < n_trials; ++tr) {
    const int start = tr * n_acc;
    int pair[2] = {-1, -1};
    int n_pair = 0;
    double pair_rho = 0.0;
    for (int a = 0; a < n_acc; ++a) {
      const int r = start + a;
      if (!ok[r]) continue;
      const double rho = pars(r, irho);
      if (std::fabs(rho) <= 1e-12) continue;
      if (n_pair < 2) pair[n_pair] = r;
      ++n_pair;
      if (n_pair == 1) pair_rho = rho;
      else if (std::fabs(rho - pair_rho) > 1e-12) {
        Rcpp::stop("rrdmswtn_tt_corr_cpp: participating rows must share one signed nonzero rho.");
      }
    }
    if (n_pair > 2) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: at most two active rows may have nonzero rho.");
    }
    if (n_pair == 2) {
      const double z1 = R::norm_rand();
      const double z2 = pair_rho * z1 +
        std::sqrt(std::fmax(0.0, 1.0 - pair_rho * pair_rho)) * R::norm_rand();
      u[static_cast<size_t>(pair[0])] = pnorm_std(z1, true, false);
      u[static_cast<size_t>(pair[1])] = pnorm_std(z2, true, false);
    }
  }

  std::vector<double> dt(static_cast<size_t>(n_rows), R_PosInf);
  for (int r = 0; r < n_rows; ++r) {
    if (!ok[r]) continue;
    const double tau = pars(r, itau);
    const double ur = u[static_cast<size_t>(r)];
    if (!R_FINITE(tau)) {
      // The exhaustion clock becomes the identity at tau = +Inf.  In the
      // correlated path there is no finite plateau to invert against, so use
      // the ordinary RDMSWTN marginal quantile directly.
      const double s = (is_ >= 0) ? pars(r, is_) : 1.0;
      dt[static_cast<size_t>(r)] = rdmswtn_quantile_cpp(
        ur, pars(r, iv), pars(r, ib), pars(r, iA), s,
        pars(r, it0), pars(r, isv), posdrift);
      continue;
    }
    const double Q = 0.5 * tau;  // tau is finite on this branch
    const double s = (is_ >= 0) ? pars(r, is_) : 1.0;
    const double FQ = prdmswtn(
        Q, pars(r, iv), pars(r, ib), pars(r, iA), s, 0.0, pars(r, isv),
        0.0, 0.0, 20, false, 1, false, posdrift, 1.0);
    if (!R_FINITE(FQ)) {
      Rcpp::stop("rrdmswtn_tt_corr_cpp: non-finite marginal response probability.");
    }
    if (ur > FQ) continue;
    const double operational = rdmswtn_operational_quantile_bounded(
        ur, FQ, Q, pars(r, iv), pars(r, ib), pars(r, iA),
        s, pars(r, isv), posdrift);
    dt[static_cast<size_t>(r)] =
      pars(r, it0) + rdmswtn_tt_qinv(operational, tau);
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, nullptr);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}
