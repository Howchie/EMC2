// Compiled posterior-prediction simulators.

#include <Rcpp.h>
#include <string>
#include <unordered_map>
#include <vector>
#include "model_rng.h"
#include "drift_factor.h"
#include "model_FRQ.h"
#include "model_BTAwL.h"
#include "ddm_functions_inline.h"

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

// Draw V from the continuous split-lognormal launch (launch == 2).  At
// delta == 0 this deliberately keeps the exact rlnorm path so the split
// variant reproduces the ordinary lognormal simulator when its extra
// parameter is zero.  Distributionally matches R's .bawd_split_rlnorm
// (R/model_BAwD.R:46-62).
inline double rsplit_lognormal_r(double mu, double sigma, double delta) {
  if (delta == 0.0) return std::exp(mu + sigma * R::norm_rand());
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NaN;
  const bool left = R::unif_rand() < h.a;
  const double z = std::fabs(R::norm_rand());
  const double y = left ? h.c - z * h.sL : h.c + z * h.sR;
  return std::exp(y);
}

}  // namespace

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

template <typename F>
inline double btawl_bisect_root(F&& f, double hi) {
  if (!(hi > 0.0) || !R_FINITE(hi)) return R_PosInf;
  double lo = 0.0;
  for (int it = 0; it < 100; ++it) {
    const double mid = 0.5 * (lo + hi);
    const double fm = f(mid);
    if (ISNAN(fm)) return R_PosInf;
    if (fm >= 0.0) hi = mid; else lo = mid;
    if (hi - lo <= 1e-11) break;
  }
  return 0.5 * (lo + hi);
}

inline double btawl_transient_state(double t, double V, double z,
                                    double k, double tau) {
  return z * std::exp(-k * t) + V * btawl_h(t, k, tau);
}

inline double btawl_sustained_state(double t, double V, double z,
                                    double k, double tau) {
  return z * std::exp(-k * t) + V * btawl_hs(t, k, tau);
}

inline double btawl_hit_time_transient_cpp(double V, double z, double b,
                                           double k, double tau) {
  if (!R_FINITE(V) || V <= 0.0 || !R_FINITE(z) || !(tau > 0.0) ||
      !R_FINITE(b) || !R_FINITE(k) || k < 0.0)
    return (R_FINITE(z) && R_FINITE(b) && z >= b) ? 0.0 : R_PosInf;
  if (z >= b) return 0.0;

  auto state_minus_b = [&](double t) {
    return btawl_transient_state(t, V, z, k, tau) - b;
  };

  if (k <= BTAWL_K_EPS) {
    double hi = std::fmax(tau, 1.0);
    const double cap = 1e12 * std::fmax(tau, 1.0);
    for (int it = 0; it < 200 && state_minus_b(hi) < 0.0 && hi < cap; ++it)
      hi = std::fmin(2.0 * hi, cap);
    if (state_minus_b(hi) < 0.0) return R_PosInf;
    return btawl_bisect_root(state_minus_b, hi);
  }

  const double tm = btawl_tmax(k, tau);
  if (!(tm > 0.0) || !R_FINITE(tm) || state_minus_b(tm) < 0.0)
    return R_PosInf;
  return btawl_bisect_root(state_minus_b, tm);
}

inline double btawl_hit_time_sustained_cpp(double V, double z, double b,
                                           double k, double tau_s) {
  if (!R_FINITE(V) || V <= 0.0 || !R_FINITE(z) || !(tau_s > 0.0) ||
      !R_FINITE(b) || !R_FINITE(k) || k < 0.0)
    return (R_FINITE(z) && R_FINITE(b) && z >= b) ? 0.0 : R_PosInf;
  if (z >= b) return 0.0;
  if (k > BTAWL_K_EPS && V <= k * b) return R_PosInf;

  auto state_minus_b = [&](double t) {
    return btawl_sustained_state(t, V, z, k, tau_s) - b;
  };
  double hi = std::fmax(tau_s, (k > BTAWL_K_EPS) ? 1.0 / k : 1.0);
  if (!R_FINITE(hi) || !(hi > 0.0)) return R_PosInf;
  const double cap = 1e12 * std::fmax(tau_s, 1.0);
  for (int it = 0; it < 200 && state_minus_b(hi) < 0.0 && hi < cap; ++it)
    hi = std::fmin(2.0 * hi, cap);
  if (state_minus_b(hi) < 0.0) return R_PosInf;
  return btawl_bisect_root(state_minus_b, hi);
}

inline double btawl_tau_for_row(const Rcpp::NumericMatrix& pars,
                                const std::unordered_map<std::string, int>& ci,
                                int row, bool prefer_ttrans) {
  if (prefer_ttrans && ci.count("Ttrans")) {
    const double Ttrans = pars(row, ci.at("Ttrans"));
    if (ci.count("tau")) return pars(row, ci.at("tau"));
    if (ci.count("tau_t")) return pars(row, ci.at("tau_t"));
    return btawl_tau_from_ttrans(pars(row, ci.at("k")), Ttrans);
  }
  if (ci.count("tau")) return pars(row, ci.at("tau"));
  if (ci.count("tau_t")) return pars(row, ci.at("tau_t"));
  if (ci.count("Ttrans"))
    return btawl_tau_from_ttrans(pars(row, ci.at("k")), pars(row, ci.at("Ttrans")));
  return R_NaN;
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

namespace {

// Invert the defective lower-boundary CDF of a fixed-parameter bounded Wiener
// process. The upper boundary is obtained by reflecting drift and start point.
double ddm_wiener_quantile(double probability, double a, double v, double w) {
  if (!(probability > 0.0) || !(a > 0.0) || !(w > 0.0) || !(w < 1.0)) {
    return 0.0;
  }

  const double log_target = std::log(probability);
  double hi = std::fmax(1.0, a * a);
  if (std::fabs(v) > 1e-10) {
    hi = std::fmax(hi, 4.0 * a / std::fabs(v));
  }

  bool bracketed = false;
  for (int i = 0; i < 120; ++i) {
    const double log_cdf = pwiener_inline(hi, a, v, w, 5e-3, 0, 1);
    if (R_FINITE(log_cdf) && log_cdf >= log_target) {
      bracketed = true;
      break;
    }
    if (!(hi < 1e12)) break;
    hi *= 2.0;
  }
  if (!bracketed) return R_PosInf;

  double lo = 0.0;
  for (int i = 0; i < 70; ++i) {
    const double mid = 0.5 * (lo + hi);
    const double log_cdf = pwiener_inline(mid, a, v, w, 5e-3, 0, 1);
    if (R_FINITE(log_cdf) && log_cdf >= log_target) hi = mid;
    else lo = mid;
    if (hi - lo <= 1e-10 * std::fmax(1.0, hi)) break;
  }
  return 0.5 * (lo + hi);
}

double ddm_wiener_finish(double a, double v, double w, int& response) {
  if (!(a > 0.0) || !(w > 0.0) || !(w < 1.0) || !R_FINITE(v)) {
    return R_PosInf;
  }

  double p_lower = std::exp(ddm_logP(0, a, v, w));
  if (!R_FINITE(p_lower)) p_lower = (v < 0.0) ? 1.0 : 0.0;
  p_lower = std::fmin(1.0, std::fmax(0.0, p_lower));

  const bool lower = R::unif_rand() < p_lower;
  response = lower ? 1 : 2;
  const double p_response = lower ? p_lower : (1.0 - p_lower);
  if (!(p_response > 0.0)) return R_PosInf;

  double u = R::unif_rand();
  u = std::fmin(std::fmax(u, 1e-15), 1.0 - 1e-15);
  double target = u * p_response;
  const double endpoint_v = lower ? v : -v;
  const double endpoint_w = lower ? w : 1.0 - w;
  const double endpoint_p = std::exp(ddm_logP(0, a, endpoint_v, endpoint_w));
  if (R_FINITE(endpoint_p) && target >= endpoint_p) {
    target = std::nextafter(endpoint_p, 0.0);
  }
  return ddm_wiener_quantile(target, a, endpoint_v, endpoint_w);
}

}  // namespace

// Single-process bounded DDM simulator. `pars` is the trial-wise matrix after
// DDM::Ttransform(), so Z and SZ are on their natural relative scales.
// [[Rcpp::export]]
Rcpp::List rddm_cpp(Rcpp::NumericMatrix pars,
                    Rcpp::CharacterVector response_levels,
                    Rcpp::LogicalVector ok) {
  if (response_levels.size() != 2) {
    Rcpp::stop("rddm_cpp: the DDM requires exactly two response levels.");
  }
  const int n_rows = pars.nrow();
  if (n_rows <= 0 || ok.size() != n_rows) {
    Rcpp::stop("rddm_cpp: invalid parameter or ok dimensions.");
  }

  const auto ci = col_index_map(pars);
  const int iv = ci.at("v");
  const int ia = ci.at("a");
  const int isv = ci.at("sv");
  const int it0 = ci.at("t0");
  const int ist0 = ci.at("st0");
  const int is = ci.at("s");
  const int iZ = ci.at("Z");
  const int iSZ = ci.at("SZ");

  std::vector<int> response(static_cast<size_t>(n_rows), 0);
  std::vector<double> rt(static_cast<size_t>(n_rows), NA_REAL);
  std::vector<int> omitted(static_cast<size_t>(n_rows), 0);

  for (int i = 0; i < n_rows; ++i) {
    if (!ok[i]) continue;

    const double a = pars(i, ia);
    const double s = pars(i, is);
    const double sv = pars(i, isv);
    const double st0 = pars(i, ist0);
    const double Z = pars(i, iZ);
    const double SZ = pars(i, iSZ);
    const double t0_base = pars(i, it0);
    if (!R_FINITE(pars(i, iv)) || !R_FINITE(a) || !R_FINITE(s) ||
        !R_FINITE(sv) || !R_FINITE(st0) || !R_FINITE(t0_base) ||
        !R_FINITE(Z) || !R_FINITE(SZ) || !(a > 0.0) || !(s > 0.0) ||
        sv < 0.0 || st0 < 0.0 ||
        !(Z > 0.0) || !(Z < 1.0) || SZ < 0.0 ||
        !(Z - 0.5 * SZ > 0.0) || !(Z + 0.5 * SZ < 1.0)) {
      omitted[static_cast<size_t>(i)] = 1;
      continue;
    }

    const double v = (sv > 0.0) ?
      (pars(i, iv) + sv * R::norm_rand()) / s : pars(i, iv) / s;
    const double w = (SZ > 0.0) ?
      Z + SZ * (R::unif_rand() - 0.5) : Z;
    int boundary = 0;
    const double fpt = ddm_wiener_finish(a / s, v, w, boundary);
    if (!R_FINITE(fpt)) {
      omitted[static_cast<size_t>(i)] = 1;
      continue;
    }

    response[static_cast<size_t>(i)] = boundary;
    const double t0 = t0_base +
      ((st0 > 0.0) ? st0 * R::unif_rand() : 0.0);
    const double finish = t0 + fpt;
    if (!R_FINITE(finish)) {
      response[static_cast<size_t>(i)] = 0;
      omitted[static_cast<size_t>(i)] = 1;
      continue;
    }
    rt[static_cast<size_t>(i)] = finish;
  }

  return pack_result(response, rt, omitted, false, std::vector<int>());
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
  const bool logn = (launch == 1 || launch == 2);
  const bool split = (launch == 2);
  if (split && !ci.count("delta"))
    Rcpp::stop("rbawl_cpp: the split-lognormal launch requires a 'delta' column.");
  const int iv = ci.at(logn ? "mu" : "v"), isv = ci.at(logn ? "sigma" : "sv"),
            ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ik = ci.at("k"), ilg = ci.at("lambda_g"), ilk = ci.at("lambda_k");
  const int idelta = split ? ci.at("delta") : -1;
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
      : (split ? rsplit_lognormal_r(pars(r, iv), pars(r, isv), pars(r, idelta))
               : (logn ? R::rlnorm(pars(r, iv), pars(r, isv))
                       : rtnorm_lower_r(pars(r, iv), pars(r, isv), lo)));
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

// [[Rcpp::export]]
Rcpp::List rbta_wl_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                       Rcpp::LogicalVector ok, int mode, bool posdrift,
                       int launch, bool endpoint_chart = false) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbta_wl_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows)
    Rcpp::stop("rbta_wl_cpp: ok has the wrong length.");
  if (mode < 0 || mode > 2)
    Rcpp::stop("rbta_wl_cpp: mode must be 0 (transient), 1 (sustained), or 2 (full).");
  if (launch < BTAWL_LAUNCH_NORMAL || launch > BTAWL_LAUNCH_SPLITLOGNORMAL)
    Rcpp::stop("rbta_wl_cpp: launch must be 0 (normal), 1 (lognormal), or 2 (split-lognormal).");

  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == BTAWL_LAUNCH_LOGNORMAL ||
                     launch == BTAWL_LAUNCH_SPLITLOGNORMAL);
  const bool split = (launch == BTAWL_LAUNCH_SPLITLOGNORMAL);
  const char* p1_name = logn ? "mu" : "v";
  const char* p2_name = logn ? "sigma" : "sv";
  if (!ci.count(p1_name) || !ci.count(p2_name))
    Rcpp::stop("rbta_wl_cpp: missing launch columns '%s' and/or '%s'.", p1_name, p2_name);
  if (split && !ci.count("delta"))
    Rcpp::stop("rbta_wl_cpp: the split-lognormal launch requires a 'delta' column.");
  for (const char* nm : {"A", "t0", "k"})
    if (!ci.count(nm)) Rcpp::stop("rbta_wl_cpp: missing parameter column '%s'.", nm);
  if (!ci.count("b") && !ci.count("B"))
    Rcpp::stop("rbta_wl_cpp: missing threshold column 'b' (or upper-bound column 'B').");
  if (mode == 1 || mode == 2)
    if (!ci.count("tau_s")) Rcpp::stop("rbta_wl_cpp: sustained mode requires 'tau_s'.");
  if (mode == 2 && !ci.count("pi"))
    Rcpp::stop("rbta_wl_cpp: full local-race mode requires 'pi'.");

  const int ip1 = ci.at(p1_name), ip2 = ci.at(p2_name);
  const int idelta = split ? ci.at("delta") : -1;
  const int iA = ci.at("A"), it0 = ci.at("t0"), ik = ci.at("k");
  const bool has_b = ci.count("b") != 0;
  const int ib = has_b ? ci.at("b") : ci.at("B");
  const int itau_s = ci.count("tau_s") ? ci.at("tau_s") : -1;
  const int ipi = ci.count("pi") ? ci.at("pi") : -1;

  auto draw_launch = [&](double p1, double p2, double delta) {
    if (launch == BTAWL_LAUNCH_LOGNORMAL)
      return R::rlnorm(p1, p2);
    if (launch == BTAWL_LAUNCH_SPLITLOGNORMAL)
      return rsplit_lognormal_r(p1, p2, delta);
    return rtnorm_lower_r(p1, p2, posdrift ? 0.0 : R_NegInf);
  };

  std::vector<double> dt(static_cast<size_t>(n_rows), R_PosInf);
  std::vector<int> ok_row(static_cast<size_t>(n_rows), 0);
  for (int r = 0; r < n_rows; ++r) ok_row[static_cast<size_t>(r)] = (ok[r] == TRUE) ? 1 : 0;

  for (int r = 0; r < n_rows; ++r) {
    if (!ok_row[static_cast<size_t>(r)]) continue;
    const double A = pars(r, iA);
    const double b = pars(r, ib) + (has_b ? 0.0 : A);
    const double k = pars(r, ik);
    const double t0 = pars(r, it0);
    const double delta = split ? pars(r, idelta) : 0.0;
    if (!R_FINITE(A) || !R_FINITE(b) || !R_FINITE(k) || k < 0.0 || !R_FINITE(t0)) continue;

    if (mode == 0) {
      const double V = draw_launch(pars(r, ip1), pars(r, ip2), delta);
      const double z = A * R::unif_rand();
      const double tau = btawl_tau_for_row(pars, ci, r, endpoint_chart);
      const double hit = btawl_hit_time_transient_cpp(V, z, b, k, tau);
      dt[static_cast<size_t>(r)] = R_FINITE(hit) ? hit + t0 : R_PosInf;
      continue;
    }

    if (mode == 1) {
      const double V = draw_launch(pars(r, ip1), pars(r, ip2), delta);
      const double z = A * R::unif_rand();
      const double hit = btawl_hit_time_sustained_cpp(
        V, z, b, k, pars(r, itau_s));
      dt[static_cast<size_t>(r)] = R_FINITE(hit) ? hit + t0 : R_PosInf;
      continue;
    }

    const double pi = pars(r, ipi);
    if (!R_FINITE(pi)) continue;
    double V_T = 0.0, V_S = 0.0;
    double p1_T = 0.0, p2_T = 0.0, p1_S = 0.0, p2_S = 0.0;
    if (pi <= 1e-14) {
      p1_T = pars(r, ip1); p2_T = pars(r, ip2);
      V_T = draw_launch(p1_T, p2_T, delta);
    } else if (pi >= 1.0 - 1e-14) {
      p1_S = pars(r, ip1); p2_S = pars(r, ip2);
      V_S = draw_launch(p1_S, p2_S, delta);
    } else if (logn) {
      p1_T = pars(r, ip1) + std::log1p(-pi);
      p1_S = pars(r, ip1) + std::log(pi);
      p2_T = p2_S = pars(r, ip2);
    } else {
      p1_T = pars(r, ip1) * (1.0 - pi);
      p1_S = pars(r, ip1) * pi;
      p2_T = pars(r, ip2) * (1.0 - pi);
      p2_S = pars(r, ip2) * pi;
    }

    const double z_T = A * R::unif_rand();
    const double z_S = A * R::unif_rand();
    if (pi > 1e-14 && pi < 1.0 - 1e-14) {
      V_T = draw_launch(p1_T, p2_T, delta);
      V_S = draw_launch(p1_S, p2_S, delta);
    }

    const double tau_t = btawl_tau_for_row(pars, ci, r, endpoint_chart);
    const double t_T = btawl_hit_time_transient_cpp(V_T, z_T, b, k, tau_t);
    const double t_S = btawl_hit_time_sustained_cpp(V_S, z_S, b, k, pars(r, itau_s));
    const double hit = std::fmin(t_T, t_S);
    dt[static_cast<size_t>(r)] = R_FINITE(hit) ? hit + t0 : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, nullptr, &ok_row);
  for (int tr = 0; tr < n_trials; ++tr)
    if (!ok_row[static_cast<size_t>(tr * n_acc)]) res.omitted[static_cast<size_t>(tr)] = 1;
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
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
                     Rcpp::LogicalVector ok, int launch, bool posdrift,
                     double gamma = 0.0, double rho = 0.0) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbawd_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rbawd_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == 1 || launch == 2);
  const bool split = (launch == 2);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawd_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (split && !ci.count("delta"))
    Rcpp::stop("rbawd_cpp: the split-lognormal launch requires a 'delta' column.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawd_cpp: the normal launch requires columns 'v' and 'sv'.");
  // Lookup is by NAME, not by position: `b` is required, and passing `B`
  // instead used to fail as an unhandled std::map::at ("_Map_base::at").
  for (const char* nm : {"b", "A", "t0", "k", "ell"})
    if (!ci.count(nm)) Rcpp::stop("rbawd_cpp: missing parameter column '%s'.", nm);
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int idelta = split ? ci.at("delta") : -1;
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ik = ci.at("k"), iell = ci.at("ell");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = split
      ? rsplit_lognormal_r(pars(r, ip1), pars(r, ip2), pars(r, idelta))
      : (logn
          ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
          : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf));
    const double z = pars(r, iA) * R::unif_rand();
    const double u = bawd_hit_time_r(V, pars(r, ib) - z, pars(r, ik),
                                     pars(r, iell), gamma, rho);
    dt[r] = (R_FINITE(u) && u >= 0.0) ? u : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// BAwF simulator.  pars columns: b, A, t0, k plus the launch pair, named
// (v, sv) for launch = 0 (truncated normal) and (mu, sigma) for launch = 1
// (lognormal).  There is no clearance column: the fading rate k is the only
// dynamic parameter.  `launch` must be the value BAwF() derived from its
// drift_distribution argument -- the same one that reaches the likelihood as
// ctx->bawd_launch -- or simulation and estimation describe different models.
// Matches R's rBAwF (R/model_BAwF.R).
// [[Rcpp::export]]
Rcpp::List rbawf_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                     Rcpp::LogicalVector ok, int launch, bool posdrift,
                     double rho = 0.0) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbawf_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rbawf_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == 1 || launch == 2);
  const bool split = (launch == 2);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawf_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (split && !ci.count("delta"))
    Rcpp::stop("rbawf_cpp: the split-lognormal launch requires a 'delta' column.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawf_cpp: the normal launch requires columns 'v' and 'sv'.");
  // Lookup is by NAME, not by position: `b` is required, and passing `B`
  // instead used to fail as an unhandled std::map::at ("_Map_base::at").
  for (const char* nm : {"b", "A", "t0", "k"})
    if (!ci.count(nm)) Rcpp::stop("rbawf_cpp: missing parameter column '%s'.", nm);
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int idelta = split ? ci.at("delta") : -1;
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ik = ci.at("k");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = split
      ? rsplit_lognormal_r(pars(r, ip1), pars(r, ip2), pars(r, idelta))
      : (logn
          ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
          : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf));
    const double z = pars(r, iA) * R::unif_rand();
    // b and z are passed separately: the fading multiplies the start point.
    const double u = bawf_hit_time_r(V, pars(r, ib), z, pars(r, ik), rho);
    dt[r] = (R_FINITE(u) && u >= 0.0) ? u : R_PosInf;
  }

  RaceOut res = resolve_race(dt, n_acc, n_trials, &t0col, &ok_row);
  std::vector<int> isTime;
  const bool has_isTime = resolve_time_level(res.R, isTime, lR_levels);
  return pack_result(res.R, res.rt, res.omitted, has_isTime, isTime);
}

// BAwR simulator.  The launch pair is (v, sv) for launch = 0 and (mu, sigma)
// for launch = 1.  Decay is in physical time only, so the crossing needs the
// distance b - z rather than b and z separately (contrast rbawf_cpp).
// [[Rcpp::export]]
Rcpp::List rbawr_cpp(Rcpp::NumericMatrix pars, Rcpp::CharacterVector lR_levels,
                     Rcpp::LogicalVector ok, int launch, bool posdrift) {
  const int n_acc = lR_levels.size();
  const int n_rows = pars.nrow();
  if (n_acc <= 0 || n_rows <= 0 || n_rows % n_acc != 0)
    Rcpp::stop("rbawr_cpp: invalid accumulator/parameter dimensions.");
  if (ok.size() != n_rows) Rcpp::stop("rbawr_cpp: ok has the wrong length.");
  const int n_trials = n_rows / n_acc;
  const auto ci = col_index_map(pars);
  const bool logn = (launch == 1 || launch == 2);
  const bool split = (launch == 2);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawr_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (split && !ci.count("delta"))
    Rcpp::stop("rbawr_cpp: the split-lognormal launch requires a 'delta' column.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawr_cpp: the normal launch requires columns 'v' and 'sv'.");
  // Lookup is by NAME, not by position: `b` is required, and passing `B`
  // instead used to fail as an unhandled std::map::at ("_Map_base::at").
  for (const char* nm : {"b", "A", "t0", "kappa", "p"})
    if (!ci.count(nm)) Rcpp::stop("rbawr_cpp: missing parameter column '%s'.", nm);
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int idelta = split ? ci.at("delta") : -1;
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0"),
            ika = ci.at("kappa"), ipw = ci.at("p");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = split
      ? rsplit_lognormal_r(pars(r, ip1), pars(r, ip2), pars(r, idelta))
      : (logn
          ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
          : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf));
    const double z = pars(r, iA) * R::unif_rand();
    const double u = bawr_hit_time_r(V, pars(r, ib) - z, pars(r, ika),
                                     pars(r, ipw));
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
  const bool logn = (launch == 1 || launch == 2);
  const bool split = (launch == 2);
  if (logn && !(ci.count("mu") && ci.count("sigma")))
    Rcpp::stop("rbawdp_cpp: the lognormal launch requires columns 'mu' and 'sigma'.");
  if (split && !ci.count("delta"))
    Rcpp::stop("rbawdp_cpp: the split-lognormal launch requires a 'delta' column.");
  if (!logn && !(ci.count("v") && ci.count("sv")))
    Rcpp::stop("rbawdp_cpp: the normal launch requires columns 'v' and 'sv'.");
  for (const char* nm : {"b", "A", "t0", "k", "lambda"})
    if (!ci.count(nm)) Rcpp::stop("rbawdp_cpp: missing parameter column '%s'.", nm);
  const int ip1 = logn ? ci.at("mu") : ci.at("v");
  const int ip2 = logn ? ci.at("sigma") : ci.at("sv");
  const int idelta = split ? ci.at("delta") : -1;
  const int ib = ci.at("b"), iA = ci.at("A"), it0 = ci.at("t0");
  const int ik = ci.at("k"), ilambda = ci.at("lambda");

  std::vector<double> dt(n_rows, R_PosInf);
  std::vector<double> t0col(n_rows);
  std::vector<int> ok_row(n_rows);
  for (int r = 0; r < n_rows; ++r) {
    t0col[r] = pars(r, it0);
    ok_row[r] = ok[r] ? 1 : 0;
    if (!ok[r]) continue;
    const double V = split
      ? rsplit_lognormal_r(pars(r, ip1), pars(r, ip2), pars(r, idelta))
      : (logn
          ? std::exp(pars(r, ip1) + pars(r, ip2) * R::norm_rand())
          : rtnorm_lower_r(pars(r, ip1), pars(r, ip2), posdrift ? 0.0 : R_NegInf));
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
    if (!ok[r] || !R_FINITE(t0col[r])) {
      dt[r] = R_PosInf;
      continue;
    }
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
