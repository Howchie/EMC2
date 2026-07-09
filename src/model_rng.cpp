// C++ simulation kernels for posterior prediction, one entry point per race
// family (LBA, BAwL, RDM, RDMSWTN). See rfun_port_plan.md for the design and
// R/model_LBA.R (rLBA, rBAwL) / R/model_RDM.R (rRDM, rWald, rSWTN, rRDMSWTN)
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

using namespace Rcpp;

namespace {

std::unordered_map<std::string, int> col_index_map(const Rcpp::NumericMatrix& pars) {
  Rcpp::CharacterVector cn = Rcpp::colnames(pars);
  std::unordered_map<std::string, int> m;
  for (int j = 0; j < cn.size(); j++) m[Rcpp::as<std::string>(cn[j])] = j;
  return m;
}

struct RaceOut {
  std::vector<int> R;        // 1-based winner accumulator, 0 = NA
  std::vector<double> rt;    // meaningful only where R != 0
  std::vector<int> best_i;   // 0-based winner accumulator index, -1 = NA (pre time-resample)
};

// Argmin-per-trial race resolution shared by all four families.
//   dt      : n_acc x n_trials, column-major (matches pars row order)
//   t0col   : per-row t0 to add to the winner's finish time, or nullptr if
//             dt already has t0 baked in (BAwL/RDMSWTN erlang-clock path)
//   ok_row  : if non-null, trial forced to NA when ok_row[first row] == 0
//             (matches rLBA/rBAwL's `ok <- matrix(ok,nrow)[1,]` reduction;
//             rRDM/rRDMSWTN have no such reduction, pass nullptr)
RaceOut resolve_race(const std::vector<double>& dt, int n_acc, int n_trials,
                     const std::vector<double>* t0col, const std::vector<int>* ok_row) {
  RaceOut out;
  out.R.assign(n_trials, 0);
  out.rt.assign(n_trials, NA_REAL);
  out.best_i.assign(n_trials, -1);
  for (int tr = 0; tr < n_trials; tr++) {
    const int base = tr * n_acc;
    if (ok_row && !(*ok_row)[base]) continue;
    double best = R_PosInf;
    int best_a = -1;
    for (int a = 0; a < n_acc; a++) {
      const double v = dt[base + a];
      if (v < best) { best = v; best_a = a; }
    }
    if (best_a < 0) continue;  // all-Inf column -> NA
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
                       bool has_isTime, const std::vector<int>& isTime) {
  const int n = (int)R.size();
  Rcpp::IntegerVector Rv(n);
  Rcpp::NumericVector rtv(n);
  for (int i = 0; i < n; i++) {
    if (R[i] == 0) {
      Rv[i] = NA_INTEGER;
      rtv[i] = NA_REAL;
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
// R's rLBA (model_LBA.R:135-166).
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
  return pack_result(res.R, res.rt, has_isTime, isTime);
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
  return pack_result(res.R, res.rt, has_isTime, isTime);
}

// pars columns: v, sv, b, A, t0, k, lambda_g, lambda_k (+ optional omega).
// Matches R's rBAwL (model_LBA.R:347-479).
// [[Rcpp::export]]
Rcpp::List rbawl_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                     Rcpp::LogicalVector ok, bool posdrift, int erlang, bool guess, bool global) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const int iv = ci.at("v"), isv = ci.at("sv"), ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
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
    const double drift = rtnorm_lower_r(pars(r, iv), pars(r, isv), lo);
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
  return pack_result(res.R, res.rt, has_isTime, isTime);
}

// pars columns: v, b, A, t0, sv, lambda_g, lambda_k (+ optional s, omega).
// erlang_type: "none" | "local_kill" | "global_kill" | "local_guess" | "local_kill_guess".
// Matches R's rRDMSWTN (model_RDM.R:872-1014) and rSWTN (model_RDM.R:489-532).
// [[Rcpp::export]]
Rcpp::List rrdmswtn_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                        Rcpp::LogicalVector ok, int erlang_shape, std::string erlang_type,
                        bool posdrift) {
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
    if (R_FINITE(sv) && sv > 1e-12) {
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

  return pack_result(res.R, res.rt, has_isTime, isTime);
}
