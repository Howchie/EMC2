// FRQ(rate_sharing = "shared"): the trial likelihood of a race whose
// accumulators share one Gamma rate multiplier.  The kernels are in
// frq_shared.h; this file is the trial bookkeeping (RACE masks, rt bins,
// censoring, truncation, contaminants), on the same outline as the RDMSWTNcorr
// trial evaluator in correlated_likelihood.cpp.

#include <Rcpp.h>
#include "frq_shared.h"
#include "frq_shared_likelihood.h"
#include "ParamTable.h"
#include "composite_functions.h"
#include "contaminant_mixture.h"
#include "rt_bins.h"
#include "gl_quad.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace Rcpp;

namespace {

Rcpp::NumericVector frqsh_col(const Rcpp::DataFrame& df, const char* name,
                              double fill) {
  if (df.containsElementNamed(name)) return Rcpp::NumericVector(df[name]);
  return Rcpp::NumericVector(df.nrows(), fill);
}

inline std::uint64_t frqsh_bits(double x) {
  if (x == 0.0) return 0;
  std::uint64_t b = 0;
  std::memcpy(&b, &x, sizeof(b));
  return b;
}

// The accumulators racing on one unique trial, with their derived states.
struct FrqShTrial {
  int n = 0;
  int row[FRQSH_MAX_FACTORS];
  double t0[FRQSH_MAX_FACTORS];
  double accel[FRQSH_MAX_FACTORS];
  FrqPars pars[FRQSH_MAX_FACTORS];
  bool poisson = false;
  bool ok = true;
};

// 64-point Gauss-Legendre in t on [lo, hi] (hi may be +Inf), log scale.
template <typename LogF>
double frqsh_integrate_log(double lo, double hi, LogF f) {
  if (!(hi > lo) || lo == R_PosInf) return R_NegInf;
  const GLRule& rule = gl_get_rule(64);
  std::array<double, 64> vals{};
  double max_log = R_NegInf;
  for (int i = 0; i < 64; ++i) {
    const double x = rule.x[static_cast<size_t>(i)];
    double t, log_jac;
    if (hi == R_PosInf) {
      const double u = 0.5 * (x + 1.0), one = 1.0 - u;
      t = lo + u / one;
      log_jac = -M_LN2 - 2.0 * std::log(one);
    } else {
      t = 0.5 * (hi - lo) * x + 0.5 * (hi + lo);
      log_jac = std::log(0.5 * (hi - lo));
    }
    vals[static_cast<size_t>(i)] = log_jac +
      std::log(rule.w[static_cast<size_t>(i)]) + f(t);
    max_log = std::max(max_log, vals[static_cast<size_t>(i)]);
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (double v : vals) sum += std::exp(v - max_log);
  return max_log + std::log(sum);
}

class FrqShEval {
 public:
  explicit FrqShEval(const FrqShTrial& tr) : tr_(tr) {}

  // Joint log density that trial accumulator w finishes first at t.
  double cause(int w, double t) const {
    const double uw = t - tr_.t0[w];
    if (!(uw > 0.0) || !R_FINITE(uw)) return R_NegInf;
    if (tr_.poisson) {
      double out = frq_accel_log_pdf_dt(uw, tr_.pars[w], tr_.accel[w]);
      for (int i = 0; i < tr_.n && out > R_NegInf; ++i) {
        if (i == w) continue;
        const double ui = t - tr_.t0[i];
        if (ui > 0.0) out += frq_accel_log_surv_dt(ui, tr_.pars[i], tr_.accel[i]);
      }
      return ISNAN(out) ? R_NegInf : out;
    }
    const FrqShPoint pw = frqsh_point(uw, tr_.pars[w], tr_.accel[w]);
    std::array<const FrqPars*, FRQSH_MAX_FACTORS> sl{};
    std::array<FrqShPoint, FRQSH_MAX_FACTORS> pl{};
    int m = 0;
    for (int i = 0; i < tr_.n; ++i) {
      if (i == w) continue;
      sl[m] = &tr_.pars[i];
      pl[m] = frqsh_point(t - tr_.t0[i], tr_.pars[i], tr_.accel[i]);
      ++m;
    }
    const double out = frqsh_log_cause(tr_.pars[w], pw, m, sl.data(), pl.data());
    return ISNAN(out) ? R_NegInf : out;
  }

  // Joint log survivor at t; t = +Inf is the never-finish atom.
  double surv(double t) const {
    if (tr_.poisson) {
      double out = 0.0;
      for (int i = 0; i < tr_.n; ++i) {
        const double ui = t - tr_.t0[i];
        if (ui > 0.0) out += frq_accel_log_surv_dt(ui, tr_.pars[i], tr_.accel[i]);
      }
      return out;
    }
    std::array<const FrqPars*, FRQSH_MAX_FACTORS> s{};
    std::array<FrqShPoint, FRQSH_MAX_FACTORS> p{};
    for (int i = 0; i < tr_.n; ++i) {
      s[i] = &tr_.pars[i];
      p[i] = frqsh_point(t - tr_.t0[i], tr_.pars[i], tr_.accel[i]);
    }
    return frqsh_log_surv(tr_.n, s.data(), p.data());
  }

 private:
  const FrqShTrial& tr_;
};

}  // namespace

FrqSharedState build_frq_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names) {
  FrqSharedState s;
  if (n_lR <= 0 || n_trials <= 0 || n_trials % n_lR != 0)
    Rcpp::stop("FRQ shared rate: invalid race dimensions.");
  if (n_lR > FRQSH_MAX_FACTORS)
    Rcpp::stop("FRQ shared rate: at most %d accumulators are supported.",
               FRQSH_MAX_FACTORS);
  if (ctx.tw.par_index >= 0 || ctx.mix.p_index >= 0)
    Rcpp::stop("FRQ shared rate: time warps and state mixtures are not supported.");
  if (ctx.time_code > 0 || ctx.nogo_code > 0)
    Rcpp::stop("FRQ shared rate: 'time' and 'nogo' accumulators are not supported.");
  s.accel_chart = (ctx.frq_shared == 2);
  s.n_trials = n_trials;
  s.n_lR = n_lR;
  s.n_unique = n_trials / n_lR;
  s.n_out = expand.length() > 0 ? expand.length() : s.n_unique;
  s.n_par = keep_names.size();
  s.rt = dadm["rt"];
  s.rt_resolution = rt_bin_width(dadm);
  s.LT = frqsh_col(dadm, "LT", 0.0);
  s.UT = frqsh_col(dadm, "UT", R_PosInf);
  s.LC = frqsh_col(dadm, "LC", R_PosInf);
  s.UC = frqsh_col(dadm, "UC", 0.0);
  s.has_RACE = dadm.containsElementNamed("RACE");
  if (s.has_RACE && dadm.hasAttribute("RACE_nacc_by_row") &&
      dadm.hasAttribute("RACE_mask")) {
    s.race_nacc = dadm.attr("RACE_nacc_by_row");
    s.race_mask = dadm.attr("RACE_mask");
    if (s.race_nacc.size() != n_trials || s.race_mask.size() != n_trials)
      s.has_RACE = false;
  } else {
    s.has_RACE = false;
  }
  s.winner_row.assign(static_cast<size_t>(s.n_unique), -1);
  for (int r = 0; r < n_trials; ++r)
    if (winner[r]) s.winner_row[static_cast<size_t>(r / n_lR)] = r;
  s.cols.assign(static_cast<size_t>(s.n_par), nullptr);
  for (int c = 0; c < s.n_par; ++c) {
    const std::string nm = Rcpp::as<std::string>(keep_names[c]);
    s.cols[static_cast<size_t>(c)] = &table.base(0, table.base_index_for(nm));
    if (nm == "pContaminant") s.pc_col = c;
    else if (nm == "pGuess") s.pg_col = c;
    else if (nm == "omega") s.omega_col = c;
  }
  s.guess = resolve_guess_kernel(dadm);
  s.guess.pg_col = s.pg_col;
  s.valid = true;
  return s;
}

static double frq_shared_trial_loglik(
    FrqSharedState& s, int j, const FrqShTrial& tr, ContextForRaceModels* ctx,
    const Rcpp::IntegerVector& response, double min_ll) {
  if (!tr.ok) return min_ll;
  const int start = j * s.n_lR;
  const bool known = response[start] != NA_INTEGER;
  int wk = -1;
  if (known) {
    for (int i = 0; i < tr.n; ++i)
      if (tr.row[i] == s.winner_row[static_cast<size_t>(j)]) wk = i;
    if (wk < 0) return min_ll;
  }
  const double rt = s.rt[start];
  const double LT = s.LT[start], UT = s.UT[start];
  const double LC = s.LC[start], UC = s.UC[start];
  const FrqShEval ev(tr);

  auto finite_cause = [&](int w, double t) {
    if (w >= 0) return ev.cause(w, t);
    double total = R_NegInf;
    for (int i = 0; i < tr.n; ++i) total = log_sum_exp(total, ev.cause(i, t));
    return total;
  };
  auto interval_cause = [&](double lo, double hi, int w) {
    return frqsh_integrate_log(lo, hi, [&](double t) { return finite_cause(w, t); });
  };
  auto survivor = [&](double t) { return t == 0.0 ? 0.0 : ev.surv(t); };
  auto survivor_difference = [&](double lo, double hi) {
    const double a = survivor(lo);
    if (hi == R_PosInf) return a;
    return log_diff_exp(a, survivor(hi));
  };
  const bool keep_atom = ctx->retains_defective_atom(UT);
  auto upper_mass = [&](double lo) {
    return defective_upper_log_mass(survivor, lo, UT, keep_atom);
  };

  double value = R_NegInf;
  if (R_FINITE(rt) && rt > 0.0) {
    if (s.rt_resolution > 0.0) {
      const RtBin bin = rt_bin(rt, s.rt_resolution, LT, UT, LC, UC);
      double refine = 0.0;
      bool boundary = !(bin.hi > bin.lo);
      for (int i = 0; i < tr.n; ++i) {
        refine = std::max(refine, rt_bin_survivor_weight(bin, tr.t0[i]));
        boundary = boundary || rt_bin_near_support(bin, tr.t0[i]);
      }
      if (!boundary) {
        value = finite_cause(wk, rt) + std::log((bin.hi - bin.lo) / bin.h);
      } else {
        auto bin_cause = [&](int i) {
          return rt_bin_log_average(bin, tr.t0[i],
            [&](double t) { return ev.cause(i, t); }, refine);
        };
        if (known) value = bin_cause(wk);
        else for (int i = 0; i < tr.n; ++i) value = log_sum_exp(value, bin_cause(i));
      }
    } else {
      value = finite_cause(wk, rt);
    }
  } else if (rt == R_NegInf) {
    value = known ? interval_cause(LT, LC, wk) : survivor_difference(LT, LC);
  } else if (rt == R_PosInf) {
    value = known ? interval_cause(UC, UT, wk) : upper_mass(UC);
  } else if (Rcpp::NumericVector::is_na(rt)) {
    value = known
      ? log_sum_exp(interval_cause(LT, LC, wk), interval_cause(UC, UT, wk))
      : log_sum_exp(survivor_difference(LT, LC), upper_mass(UC));
  }

  if (LT != 0.0 || UT != R_PosInf) {
    // Keyed on every column the trial's states read, and the window.
    std::vector<std::uint64_t> key;
    key.reserve(static_cast<size_t>(tr.n) * emc2col::frq::N_REQ + 2);
    for (int i = 0; i < tr.n; ++i)
      for (int c = 0; c < emc2col::frq::N_REQ; ++c)
        key.push_back(frqsh_bits(s.cols[static_cast<size_t>(c)][tr.row[i]]));
    key.push_back(frqsh_bits(LT));
    key.push_back(frqsh_bits(UT));
    double log_z;
    const auto hit = s.z_cache.find(key);
    if (hit != s.z_cache.end()) {
      log_z = hit->second;
    } else {
      log_z = upper_mass(LT);
      s.z_cache.emplace(std::move(key), log_z);
    }
    if (!R_FINITE(log_z)) return min_ll;
    value -= log_z;
  }
  if (s.pc_col >= 0) {
    const double p_c = s.cols[static_cast<size_t>(s.pc_col)][start];
    if (!R_FINITE(p_c) || p_c < 0.0 || p_c >= 1.0) return min_ll;
  }
  if (s.pc_col >= 0 || (s.pg_col >= 0 && s.guess.active())) {
    const double pC = (s.pc_col >= 0) ? s.cols[static_cast<size_t>(s.pc_col)][start] : 0.0;
    const double pG = (s.pg_col >= 0 && s.guess.active())
                        ? s.cols[static_cast<size_t>(s.pg_col)][start] : 0.0;
    const GuessKernel guess = s.rt_resolution > 0.0 && R_FINITE(rt) ?
      rt_bin_guess(s.guess, rt_bin(rt, s.rt_resolution, LT, UT, LC, UC)) : s.guess;
    value = mix_contaminants_rt(value, pC, pG, guess, rt, known,
                                guess.response_is_eligible(response[start]));
  }
  return (!R_FINITE(value) || value < min_ll) ? min_ll : value;
}

double c_log_likelihood_frq_shared(
    FrqSharedState& s, const Rcpp::DataFrame& dadm, double min_ll,
    const Rcpp::LogicalVector& isok, const Rcpp::IntegerVector& expand,
    ContextForRaceModels* ctx, Rcpp::NumericVector* trial_ll_out) {
  if (ctx == nullptr || ctx->frq_shared == 0 || !s.valid)
    Rcpp::stop("FRQ shared rate: invalid likelihood context.");
  if (ctx->frq_cache) ctx->frq_cache->clear();
  s.z_cache.clear();
  namespace col = emc2col::frq;
  const double* K_ = s.cols[col::K];
  const double* v_ = s.cols[col::v];
  const double* cv_ = s.cols[col::sv];
  const double* h_ = s.cols[col::h];
  const double* t0_ = s.cols[col::t0];
  const double* pl_ = s.cols[col::plo];
  const double* ck_ = s.cols[col::clock];
  const double* om_ = s.omega_col >= 0 ? s.cols[static_cast<size_t>(s.omega_col)] : nullptr;
  // One memo per accumulator position: consecutive trials then revisit the
  // same availability state in each.
  std::vector<FrqMemo> memo(static_cast<size_t>(s.n_lR),
                            FrqMemo(ctx->frq_cache.get()));
  const Rcpp::IntegerVector response = dadm["R"];
  std::vector<double> unique_ll(static_cast<size_t>(s.n_unique), min_ll);
  FrqShTrial tr;
  for (int j = 0; j < s.n_unique; ++j) {
    const int start = j * s.n_lR;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
    tr.n = 0;
    tr.ok = true;
    double cv0 = R_NaN;
    for (int k = 0; k < n_lR_j; ++k) {
      const int r = start + k;
      if (s.has_RACE && !s.race_mask[r]) continue;
      if (!isok[r] || R_IsNA(K_[r]) || !R_FINITE(t0_[r])) { tr.ok = false; break; }
      if (om_ != nullptr && om_[r] != 0.0)
        Rcpp::stop("FRQ shared rate: omega must be 0.");
      if (tr.n == 0) cv0 = cv_[r];
      else if (cv_[r] != cv0)
        Rcpp::stop("FRQ shared rate: cv must be equal across the accumulators "
                   "of a trial (one shared rate has one CV).");
      const double accel = s.accel_chart ? ck_[r] : 0.0;
      if (s.accel_chart && !frq_accel_ok(accel)) { tr.ok = false; break; }
      const FrqPars& p = memo[static_cast<size_t>(k)].get_cv(
          K_[r], v_[r], cv_[r], h_[r], pl_[r], s.accel_chart ? 1.0 : ck_[r], 0.0);
      if (!p.ok) { tr.ok = false; break; }
      tr.row[tr.n] = r;
      tr.t0[tr.n] = t0_[r];
      tr.accel[tr.n] = accel;
      tr.pars[tr.n] = p;
      ++tr.n;
    }
    if (tr.ok && tr.n == 0) tr.ok = false;
    tr.poisson = tr.ok && tr.pars[0].poisson;
    unique_ll[static_cast<size_t>(j)] =
      frq_shared_trial_loglik(s, j, tr, ctx, response, min_ll);
  }
  double total = 0.0;
  if (expand.length() > 0) {
    for (int i = 0; i < s.n_out; ++i) {
      const double v = unique_ll[static_cast<size_t>(expand[i] - 1)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = v;
      total += v;
    }
  } else {
    for (int j = 0; j < s.n_unique; ++j) {
      const double v = unique_ll[static_cast<size_t>(j)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = v;
      total += v;
    }
  }
  return total;
}

// ---------------------------------------------------------------------------
// R entry points for one trial, used by the tests and by dFRQ_shared().
// pars holds one row per accumulator in the frq_cv / frq_accel column order.
// ---------------------------------------------------------------------------

static FrqShTrial frqsh_trial_from_matrix(const NumericMatrix& pars, bool accel_chart) {
  if (pars.ncol() < emc2col::frq::N_REQ)
    Rcpp::stop("frq_shared: pars needs the %d FRQ columns.", emc2col::frq::N_REQ);
  if (pars.nrow() > FRQSH_MAX_FACTORS || pars.nrow() < 1)
    Rcpp::stop("frq_shared: 1 to %d accumulators.", FRQSH_MAX_FACTORS);
  namespace col = emc2col::frq;
  FrqShTrial tr;
  tr.n = pars.nrow();
  for (int i = 0; i < tr.n; ++i) {
    if (pars(i, col::sv) != pars(0, col::sv))
      Rcpp::stop("frq_shared: cv must be equal across accumulators.");
    tr.row[i] = i;
    tr.t0[i] = pars(i, col::t0);
    tr.accel[i] = accel_chart ? pars(i, col::clock) : 0.0;
    const FrqAvail av = frq_derive_avail(pars(i, col::K), pars(i, col::sv),
                                         pars(i, col::h), pars(i, col::plo), 0.0);
    frq_copy_avail(tr.pars[i], av);
    frq_set_rate(tr.pars[i], av.ok, pars(i, col::v),
                 accel_chart ? 1.0 : pars(i, col::clock));
    if (!tr.pars[i].ok) tr.ok = false;
  }
  tr.poisson = tr.ok && tr.pars[0].poisson;
  return tr;
}

// [[Rcpp::export]]
NumericVector frq_shared_cause_cpp(NumericVector t, NumericMatrix pars,
                                   int winner, bool accel_chart = false) {
  const FrqShTrial tr = frqsh_trial_from_matrix(pars, accel_chart);
  if (winner < 1 || winner > tr.n) Rcpp::stop("frq_shared: bad winner index.");
  NumericVector out(t.size(), R_NegInf);
  if (!tr.ok) return out;
  const FrqShEval ev(tr);
  for (int i = 0; i < t.size(); ++i) out[i] = ev.cause(winner - 1, t[i]);
  return out;
}

// [[Rcpp::export]]
NumericVector frq_shared_surv_cpp(NumericVector t, NumericMatrix pars,
                                  bool accel_chart = false) {
  const FrqShTrial tr = frqsh_trial_from_matrix(pars, accel_chart);
  NumericVector out(t.size(), R_NaN);
  if (!tr.ok) return out;
  const FrqShEval ev(tr);
  for (int i = 0; i < t.size(); ++i) out[i] = ev.surv(t[i]);
  return out;
}

// log E[prod_l ((1 - plo_l) Q(K_l, G psi_l) + plo_l e^{-G psi_l})],
// G ~ Gamma(b, 1).  route: 0 = automatic, 1 = quadrature only.
// [[Rcpp::export]]
double frq_shared_expect_cpp(double b, NumericVector K, NumericVector log_psi,
                             NumericVector plo, int route = 0) {
  const int m = K.size();
  if (m > FRQSH_MAX_FACTORS || log_psi.size() != m || plo.size() != m)
    Rcpp::stop("frq_shared_expect_cpp: bad factor vectors.");
  if (route == 1) {
    std::array<double, FRQSH_MAX_FACTORS> k{}, lp{};
    for (int l = 0; l < m; ++l) { k[l] = K[l]; lp[l] = log_psi[l]; }
    return frqsh_log_quad(b, m, k.data(), lp.data());
  }
  std::array<FrqLow, FRQSH_MAX_FACTORS> lo{};
  std::array<FrqShFactor, FRQSH_MAX_FACTORS> f{};
  for (int l = 0; l < m; ++l) {
    lo[l] = frq_low_make(plo[l]);
    f[l].K = K[l];
    f[l].log_psi = log_psi[l];
    f[l].lo = &lo[l];
  }
  return frqsh_log_expect(b, m, f.data());
}
