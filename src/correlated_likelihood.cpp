#include "gh_quad.h"
#include "utility_functions.h"
#include "race_contract.h"
#include "correlated_likelihood.h"
#include "utils.h"
#include "ParamTable.h"
#include "gh_quad.h"
#include "gsl_utils.h"
#include "bawl_geometry.h"
#include "bawl_corr_exact.h"
#include "bawl_corr_counters.h"
#include "drift_factor.h"
#include "contaminant_mixture.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

// Defined in particle_ll.cpp; the correlated routes delegate ordinary trials
// and generic conditional nodes to the existing race implementation.
double c_log_likelihood_race(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    const int n_trials,
    Rcpp::LogicalVector winner,
    Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector isok,
    int n_lR,
    void* model_context_for_funcs,
    bool all_finite_trials,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,
    RaceSharedState* shared = nullptr,
    Rcpp::NumericVector* trial_ll_out = nullptr,
    bool apply_truncation_correction = true);

// Defined in race_integrands.cpp; used for the generic conditional truncation
// fallback when a node has a variable-accumulator or global-kill race.
double get_trunc_normaliser_rowmajor_cpp(
    const double* pars_rowmajor, const int* isok_int, RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1, double LT, double UT, int n_lR, int n_par,
    const GslIntegrationControls& gsl_ctl, void* model_specific_context,
    GslWorkspacePtr& workspace);

// Build the data-fixed correlated BAwL state once per likelihood call: lM
// role mapping, RACE masks, truncation windows, winner rows, unique-trial
// expansion, and direct ParamTable column pointers in keep_names order.
// The fast OO mapper intentionally stops at the design-matrix/natural-scale
// mapping and does not invoke the model's R Ttransform.  The role mapping
// resolved here mirrors BAwLcorr's role conversion for the likelihood path;
// R Ttransform remains the corresponding conversion for mapped parameter
// displays and simulation.
CorrDriftSharedState build_corr_drift_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names) {
  CorrDriftSharedState s;
  if (n_lR <= 0 || n_trials < 0 || n_trials % n_lR != 0) {
    Rcpp::stop("c_log_likelihood_corr_drift: invalid race dimensions.");
  }
  s.n_trials = n_trials;
  s.n_lR = n_lR;
  s.n_unique = n_trials / n_lR;
  s.n_out = (expand.length() > 0) ? expand.length() : s.n_unique;
  s.n_par = keep_names.size();
  s.rho_col = ctx.corr_drift_rho_index;
  if (s.rho_col < 0 || s.rho_col >= s.n_par) {
    Rcpp::stop("c_log_likelihood_corr_drift: rho column is outside pars.");
  }
  s.factor_cols.v = ctx.corr_drift_v_col;
  s.factor_cols.sv = ctx.corr_drift_sv_col;
  s.generic_only = ctx.corr_drift_generic_only;
  if (!s.factor_cols.valid() || s.factor_cols.v >= s.n_par ||
      s.factor_cols.sv >= s.n_par) {
    Rcpp::stop("c_log_likelihood_corr_drift: the drift mean/SD columns are "
               "outside pars; resolve_race_model_adapter() must set "
               "corr_drift_v_col and corr_drift_sv_col.");
  }

  // lM role mapping.
  bool lM_is_logical = false;
  Rcpp::LogicalVector lM_logical;
  Rcpp::IntegerVector lM_factor;
  int lM_true_code = -1;
  if (!dadm.containsElementNamed("lM")) {
    Rcpp::stop("c_log_likelihood_corr_drift: BAwLcorr requires lM from matchfun.");
  }
  SEXP lM_sexp = dadm["lM"];
  if (TYPEOF(lM_sexp) == LGLSXP && Rf_length(lM_sexp) == n_trials) {
    lM_logical = Rcpp::LogicalVector(lM_sexp);
    lM_is_logical = true;
  } else if (TYPEOF(lM_sexp) == INTSXP && Rf_length(lM_sexp) == n_trials) {
    lM_factor = Rcpp::IntegerVector(lM_sexp);
    Rcpp::CharacterVector levels = lM_factor.attr("levels");
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == "TRUE") {
        lM_true_code = i + 1;
        break;
      }
    }
    if (lM_true_code < 0 && levels.size() == 0) lM_true_code = 1;
  }
  if (!lM_is_logical && lM_true_code < 0) {
    Rcpp::stop("c_log_likelihood_corr_drift: lM must be a logical or TRUE/FALSE factor.");
  }
  s.role_correct.assign(static_cast<size_t>(n_trials), 0);
  for (int r = 0; r < n_trials; ++r) {
    const bool missing = lM_is_logical
      ? (lM_logical[r] == NA_LOGICAL)
      : (lM_factor[r] == NA_INTEGER);
    if (missing) {
      Rcpp::stop("c_log_likelihood_corr_drift: lM contains missing role values.");
    }
    const bool correct = lM_is_logical
      ? static_cast<bool>(lM_logical[r])
      : (lM_factor[r] == lM_true_code);
    s.role_correct[static_cast<size_t>(r)] = correct ? 1 : 0;
  }

  // RACE masks (attribute-cached variable-accumulator designs only).
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

  // Truncation windows.  emc2_all_finite_trials guarantees LT = 0/UT = Inf,
  // so these are all zero whenever the caller reports all-finite data.
  s.rt = dadm["rt"];
  s.LT = get_col_with_default(dadm, "LT", 0.0);
  s.UT = get_col_with_default(dadm, "UT", R_PosInf);
  s.LC = get_col_with_default(dadm, "LC", R_PosInf);
  s.UC = get_col_with_default(dadm, "UC", 0.0);
  s.has_trunc_trial.assign(static_cast<size_t>(s.n_unique), 0);
  for (int j = 0; j < s.n_unique; ++j) {
    const int start = j * n_lR;
    if (s.LT[start] != 0.0 || s.UT[start] != R_PosInf) {
      s.has_trunc_trial[static_cast<size_t>(j)] = 1;
      s.has_truncation = true;
    }
  }

  // Winner rows and unique-trial expansion map.
  s.winner_row.assign(static_cast<size_t>(s.n_unique), -1);
  for (int r = 0; r < n_trials; ++r) {
    if (winner[r]) s.winner_row[static_cast<size_t>(r / n_lR)] = r;
  }
  s.j_to_i.assign(static_cast<size_t>(s.n_unique), -1);
  if (expand.length() > 0) {
    for (int i = s.n_out - 1; i >= 0; --i)
      s.j_to_i[static_cast<size_t>(expand[i] - 1)] = i;
  } else {
    for (int j = 0; j < s.n_unique; ++j) s.j_to_i[static_cast<size_t>(j)] = j;
  }

  // Direct column pointers; base storage is refilled in place per particle.
  // base_index_for throws for unknown names exactly like materialisation.
  s.cols.assign(static_cast<size_t>(std::max(s.n_par, 16)), nullptr);
  for (int j = 0; j < s.n_par; ++j) {
    const std::string nm = Rcpp::as<std::string>(keep_names[j]);
    s.cols[static_cast<size_t>(j)] = &table.base(0, table.base_index_for(nm));
    if (nm == "pContaminant") s.pc_col = j;
    else if (nm == "pGuess") s.pg_col = j;
  }
  s.guess = resolve_guess_kernel(dadm);
  s.guess.pg_col = s.pg_col;

  s.effective_rho.assign(static_cast<size_t>(n_trials), R_NegInf);
  s.layout.assign(static_cast<size_t>(s.n_unique), BAwLCorrTrialLayout());
  s.valid = true;
  return s;
}

// Canonical per-particle trial classification (plan design rule C).  Counts
// only isok && RACE-active rows; an exact-zero loading is an independent
// singleton, not a loaded dimension.  This classifier is the sole authority
// for loaded-row discovery and the positivity dimension.
static void corr_drift_classify_particle(CorrDriftSharedState& s,
                                        const Rcpp::LogicalVector& isok,
                                        bool no_clock_eligible) {
  constexpr double kLoadingEps = 1e-14;
  const int n_lR = s.n_lR;
  for (int j = 0; j < s.n_unique; ++j) {
    BAwLCorrTrialLayout& L = s.layout[static_cast<size_t>(j)];
    L = BAwLCorrTrialLayout();
    const int start = j * n_lR;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : n_lR;
    L.winner_row = s.winner_row[static_cast<size_t>(j)];
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (s.has_RACE && !s.race_mask[row]) continue;
      ++L.n_active;
      const double rho = s.effective_rho[static_cast<size_t>(row)];
      if (!isok[row] || !R_FINITE(rho) || std::fabs(rho) > 1.0) {
        L.any_bad_row = true;
        continue;
      }
      if (std::fabs(rho) > kLoadingEps) {
        if (L.n_loaded < 2) L.loaded_row[L.n_loaded] = row;
        ++L.n_loaded;
        if (row == L.winner_row) L.winner_loaded = true;
      }
    }
    if (L.any_bad_row) {
      L.route = BAwLCorrRoute::invalid;
    } else if (L.n_loaded <= 1) {
      L.route = BAwLCorrRoute::ordinary;
    } else {
      // A two-loaded no-clock trial is the exact component route.  Trials
      // with a larger factor dimension retain the fused GH fallback; active
      // clocks retain the generic race semantics even when only two rows are
      // loaded.  Models without an affine-in-drift survivor have neither the
      // exact nor the fused kernel and always take the generic route.
      L.route = (no_clock_eligible && !s.generic_only)
        ? ((L.n_loaded == 2) ? BAwLCorrRoute::exact_pair
                             : BAwLCorrRoute::gh_no_clock)
        : BAwLCorrRoute::gh_generic_clock;
    }
  }
}

struct BAwLCorrPairData {
  int row1 = -1;
  int row2 = -1;
  double mu1 = 0.0;
  double mu2 = 0.0;
  double sd1 = 1.0;
  double sd2 = 1.0;
  double rho = 0.0;
  double normalizer = 1.0;
  bool positive = false;
};

static inline double bawl_corr_row_k(const double* const* cols, int row,
                                     const ContextForRaceModels* ctx) {
  return ctx->bawl_k_fixed_zero ? 0.0 : cols[emc2col::bawl::k][row];
}

static inline BAwLTimeGeometry bawl_corr_row_geometry(
    double t, const double* const* cols, int row,
    const ContextForRaceModels* ctx) {
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  return bawl_time_geometry(t, cols[emc2col::bawl::t0][row], A, B + A,
                            bawl_corr_row_k(cols, row, ctx));
}

static inline bool bawl_corr_make_pair_data(
    const BAwLCorrTrialLayout& layout, const CorrDriftSharedState& s,
    const double* const* cols, const ContextForRaceModels* ctx,
    BAwLCorrPairData& out) {
  if (layout.n_loaded != 2 || layout.loaded_row[0] < 0 ||
      layout.loaded_row[1] < 0) return false;
  out.row1 = layout.loaded_row[0];
  out.row2 = layout.loaded_row[1];
  out.mu1 = cols[emc2col::bawl::v][out.row1];
  out.mu2 = cols[emc2col::bawl::v][out.row2];
  out.sd1 = cols[emc2col::bawl::sv][out.row1];
  out.sd2 = cols[emc2col::bawl::sv][out.row2];
  const double r1 = s.effective_rho[static_cast<size_t>(out.row1)];
  const double r2 = s.effective_rho[static_cast<size_t>(out.row2)];
  if (!(out.sd1 > 0.0) || !(out.sd2 > 0.0) || !R_FINITE(out.mu1) ||
      !R_FINITE(out.mu2) || !R_FINITE(r1) || !R_FINITE(r2) ||
      std::fabs(r1) > 1.0 || std::fabs(r2) > 1.0) return false;
  out.rho = (r1 * r2 < 0.0 ? -1.0 : 1.0) *
    std::sqrt(std::fabs(r1 * r2));
  out.positive = ctx->use_posdrift;
  const double key[5] = {out.mu1, out.sd1, out.mu2, out.sd2, out.rho};
  for (const auto& e : s.qab_cache) {
    if (e.key[0] == key[0] && e.key[1] == key[1] && e.key[2] == key[2] &&
        e.key[3] == key[3] && e.key[4] == key[4]) {
      out.normalizer = e.value;
      return R_FINITE(out.normalizer) && out.normalizer > 0.0;
    }
  }
  out.normalizer = bawl_corr_pair_positive_normalizer(
      out.mu1, out.sd1, out.mu2, out.sd2, out.rho, out.positive);
  if (!R_FINITE(out.normalizer) || !(out.normalizer > 0.0)) return false;
  if (s.qab_cache.size() < 64) {
    s.qab_cache.push_back({{key[0], key[1], key[2], key[3], key[4]},
                           out.normalizer});
  }
  return R_FINITE(out.normalizer) && out.normalizer > 0.0;
}

static inline double bawl_corr_single_log_survival(
    double t, int row, const double* const* cols,
    const ContextForRaceModels* ctx) {
  const double t0 = cols[emc2col::bawl::t0][row];
  const double tau = t - t0;
  if (ISNAN(t) || ISNAN(t0) || !(t > 0.0) || !(tau > 0.0)) return 0.0;
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  const double v = cols[emc2col::bawl::v][row];
  const double sv = cols[emc2col::bawl::sv][row];
  const double k = bawl_corr_row_k(cols, row, ctx);
  const double log_cdf = log_ba_cdf(
      tau, A, B + A, v, sv, k, ctx->use_posdrift,
      ctx->bawl_k_fixed_zero ? LBA_DENOM_FLOOR : BAWL_DENOM_FLOOR);
  if (log_cdf == R_NegInf) return 0.0;
  if (log_cdf >= 0.0) return R_NegInf;
  return log1m_exp(log_cdf);
}

static inline double bawl_corr_single_log_density(
    double t, int row, const double* const* cols,
    const ContextForRaceModels* ctx) {
  if (!(t > 0.0) || !R_FINITE(t)) return R_NegInf;
  const double t0 = cols[emc2col::bawl::t0][row];
  const double tau = t - t0;
  if (!(tau > 0.0) || !R_FINITE(tau)) return R_NegInf;
  const double A = cols[emc2col::bawl::A][row];
  const double B = cols[emc2col::bawl::B][row];
  const double v = cols[emc2col::bawl::v][row];
  const double sv = cols[emc2col::bawl::sv][row];
  const double k = bawl_corr_row_k(cols, row, ctx);
  return log_ba_pdf(tau, A, B + A, v, sv, k, ctx->use_posdrift,
                    ctx->bawl_k_fixed_zero ? LBA_DENOM_FLOOR : BAWL_DENOM_FLOOR);
}

static inline double bawl_corr_log_pair_survival(
    const BAwLCorrPairData& pair, double t, const double* const* cols,
    const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const BAwLTimeGeometry g1 = bawl_corr_row_geometry(t, cols, pair.row1, ctx);
  const BAwLTimeGeometry g2 = bawl_corr_row_geometry(t, cols, pair.row2, ctx);
  const bool use_numeric = numeric || std::fabs(pair.rho) >= 1.0 - 1e-10;
  const BAwLCorrPairResult result = use_numeric
    ? bawl_corr_pair_survival_numeric(g1, g2, pair.mu1, pair.sd1, pair.mu2,
                                      pair.sd2, pair.rho, pair.positive,
                                      pair.normalizer)
    : bawl_corr_pair_survival_exact(g1, g2, pair.mu1, pair.sd1, pair.mu2,
                                    pair.sd2, pair.rho, pair.positive,
                                    pair.normalizer);
  if (status_out != nullptr) *status_out = result.status;
  if (!(result.value > 0.0) || !R_FINITE(result.value)) return R_NegInf;
  return std::log(result.value);
}

static inline double bawl_corr_log_pair_cause(
    const BAwLCorrPairData& pair, int winner_row, double t,
    const double* const* cols, const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const bool first = winner_row == pair.row1;
  const int rw = first ? pair.row1 : pair.row2;
  const int rl = first ? pair.row2 : pair.row1;
  const BAwLTimeGeometry gw = bawl_corr_row_geometry(t, cols, rw, ctx);
  const BAwLTimeGeometry gl = bawl_corr_row_geometry(t, cols, rl, ctx);
  const bool use_numeric = numeric || std::fabs(pair.rho) >= 1.0 - 1e-10;
  const BAwLCorrPairResult result = use_numeric
    ? bawl_corr_pair_cause_numeric(
        gw, gl, cols[emc2col::bawl::v][rw], cols[emc2col::bawl::sv][rw],
        cols[emc2col::bawl::v][rl], cols[emc2col::bawl::sv][rl], pair.rho,
        pair.positive, pair.normalizer)
    : bawl_corr_pair_cause_exact(
        gw, gl, cols[emc2col::bawl::v][rw], cols[emc2col::bawl::sv][rw],
        cols[emc2col::bawl::v][rl], cols[emc2col::bawl::sv][rl], pair.rho,
        pair.positive, pair.normalizer);
  if (status_out != nullptr) *status_out = result.status;
  if (!(result.value > 0.0) || !R_FINITE(result.value)) return R_NegInf;
  return std::log(result.value);
}

static inline bool bawl_corr_row_active(const CorrDriftSharedState& s, int row) {
  return !s.has_RACE || s.race_mask[row];
}

static inline double bawl_corr_log_component_survival(
    const CorrDriftSharedState& s, const BAwLCorrPairData& pair, int j,
    double t, const double* const* cols, const ContextForRaceModels* ctx,
    bool numeric, BAwLCorrMomentStatus* status_out = nullptr) {
  BAwLCorrMomentStatus pair_status = BAwLCorrMomentStatus::ok;
  double out = bawl_corr_log_pair_survival(pair, t, cols, ctx, numeric,
                                           &pair_status);
  if (status_out != nullptr) *status_out = pair_status;
  if (pair_status == BAwLCorrMomentStatus::unstable ||
      pair_status == BAwLCorrMomentStatus::invalid) return R_NaN;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!bawl_corr_row_active(s, row) || row == pair.row1 || row == pair.row2) continue;
    const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
    if (!R_FINITE(ls)) return R_NegInf;
    out += ls;
  }
  return out;
}

static inline double bawl_corr_log_component_cause(
    const CorrDriftSharedState& s, const BAwLCorrPairData& pair, int j,
    int winner_row, double t, const double* const* cols,
    const ContextForRaceModels* ctx, bool numeric,
    BAwLCorrMomentStatus* status_out = nullptr) {
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  if (winner_row == pair.row1 || winner_row == pair.row2) {
    BAwLCorrMomentStatus ps = BAwLCorrMomentStatus::ok;
    double out = bawl_corr_log_pair_cause(pair, winner_row, t, cols, ctx,
                                          numeric, &ps);
    if (status_out != nullptr) *status_out = ps;
    if (ps == BAwLCorrMomentStatus::unstable ||
        ps == BAwLCorrMomentStatus::invalid) return R_NaN;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!bawl_corr_row_active(s, row) || row == pair.row1 || row == pair.row2) continue;
      const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
      if (!R_FINITE(ls)) return R_NegInf;
      out += ls;
    }
    return out;
  }
  if (!bawl_corr_row_active(s, winner_row)) return R_NegInf;
  double out = bawl_corr_single_log_density(t, winner_row, cols, ctx);
  if (!R_FINITE(out)) return out;
  BAwLCorrMomentStatus ps = BAwLCorrMomentStatus::ok;
  const double lp = bawl_corr_log_pair_survival(pair, t, cols, ctx, numeric, &ps);
  if (status_out != nullptr) *status_out = ps;
  if (ps == BAwLCorrMomentStatus::unstable || ps == BAwLCorrMomentStatus::invalid)
    return R_NaN;
  if (!R_FINITE(lp)) return R_NegInf;
  out += lp;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!bawl_corr_row_active(s, row) || row == winner_row ||
        row == pair.row1 || row == pair.row2) continue;
    const double ls = bawl_corr_single_log_survival(t, row, cols, ctx);
    if (!R_FINITE(ls)) return R_NegInf;
    out += ls;
  }
  return out;
}

template <typename LogIntegrand>
static double bawl_corr_integrate_log_interval(double lo, double hi,
                                               LogIntegrand log_integrand) {
  if (!(hi > lo) && !(lo == 0.0 && hi == R_PosInf)) return R_NegInf;
  const GLRule& rule = gl_get_rule(64);
  double max_log = R_NegInf;
  std::array<double, 64> values{};
  for (int i = 0; i < 64; ++i) {
    const double x = rule.x[i];
    double t = 0.0;
    double log_jac = 0.0;
    if (hi == R_PosInf) {
      const double u = 0.5 * (x + 1.0);
      const double one = 1.0 - u;
      t = lo + u / one;
      log_jac = -std::log(2.0) - 2.0 * std::log(one);
    } else {
      t = 0.5 * (hi - lo) * x + 0.5 * (hi + lo);
      log_jac = std::log(0.5 * (hi - lo));
    }
    values[static_cast<size_t>(i)] = log_jac + std::log(rule.w[i]) +
      log_integrand(t);
    max_log = std::max(max_log, values[static_cast<size_t>(i)]);
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (int i = 0; i < 64; ++i)
    sum += std::exp(values[static_cast<size_t>(i)] - max_log);
  return max_log + std::log(sum);
}

struct BAwLCorrExactTrialResult {
  double log_likelihood = R_NegInf;
  BAwLCorrMomentStatus status = BAwLCorrMomentStatus::invalid;
  bool used_numeric = false;
  bool point_start = false;
};

static BAwLCorrExactTrialResult bawl_corr_exact_trial_loglik(
    const CorrDriftSharedState& s, const BAwLCorrTrialLayout& layout, int j,
    const double* const* cols, const ContextForRaceModels* ctx,
    const Rcpp::IntegerVector& response, const Rcpp::LogicalVector& winner,
    double min_ll, bool force_numeric = false) {
  BAwLCorrExactTrialResult out;
  BAwLCorrPairData pair;
  if (!bawl_corr_make_pair_data(layout, s, cols, ctx, pair)) return out;
  const int start = j * s.n_lR;
  const int response_code = response[start];
  const bool known = response_code != NA_INTEGER;
  const int known_row = s.winner_row[static_cast<size_t>(j)];
  const double rt = s.rt[start];
  const double LT = s.LT[start];
  const double UT = s.UT[start];
  const double LC = s.LC[start];
  const double UC = s.UC[start];
  bool numeric = force_numeric;
  BAwLCorrMomentStatus st = BAwLCorrMomentStatus::ok;
  auto absorb_status = [&](BAwLCorrMomentStatus value) {
    if (value == BAwLCorrMomentStatus::invalid) st = value;
    else if (value == BAwLCorrMomentStatus::unstable &&
             st != BAwLCorrMomentStatus::invalid) st = value;
  };
  auto component_survival_mode = [&](double t, bool numeric_survival,
                                     BAwLCorrMomentStatus* status) {
    BAwLCorrMomentStatus local = BAwLCorrMomentStatus::ok;
    const double value = bawl_corr_log_component_survival(
        s, pair, j, t, cols, ctx, numeric_survival, &local);
    if (status != nullptr) *status = local;
    absorb_status(local);
    return value;
  };
  auto component_survival = [&](double t, BAwLCorrMomentStatus* status) {
    return component_survival_mode(t, numeric, status);
  };
  auto component_cause = [&](int row, double t, BAwLCorrMomentStatus* status) {
    BAwLCorrMomentStatus local = BAwLCorrMomentStatus::ok;
    const double value = bawl_corr_log_component_cause(
        s, pair, j, row, t, cols, ctx, numeric, &local);
    if (status != nullptr) *status = local;
    absorb_status(local);
    return value;
  };

  auto finite_cause = [&](double t, int row, BAwLCorrMomentStatus* st) {
    if (row >= 0) return component_cause(row, t, st);
    double total = R_NegInf;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int r = start + k;
      if (!bawl_corr_row_active(s, r)) continue;
      total = log_sum_exp(total, component_cause(r, t, st));
    }
    return total;
  };
  auto interval_cause = [&](double lo, double hi, int row,
                            BAwLCorrMomentStatus* st) {
    return bawl_corr_integrate_log_interval(lo, hi, [&](double t) {
      return finite_cause(t, row, st);
    });
  };
  auto survival_difference = [&](double lo, double hi,
                                 BAwLCorrMomentStatus* st) {
    const double a = (lo == 0.0) ? 0.0 : component_survival(lo, st);
    if (hi == R_PosInf) return a;
    const double b = component_survival(hi, st);
    if (a == R_NegInf) return R_NegInf;
    const double z = log_diff_exp(a, b);
    return z;
  };
  bool denominator_requires_numeric = false;
  auto stable_survival_difference = [&](double lo, double hi,
                                        BAwLCorrMomentStatus* st) {
    // The exact rectangle formula obtains a survivor by subtracting four BVN
    // CDF corners.  In a rare positive orthant those corners can agree to all
    // available digits while leaving a small, positive residual that looks
    // valid.  Dividing by that residual as an LT/UT normaliser then rewards
    // the parameter point with an arbitrarily large finite likelihood.
    //
    // The deterministic pair route integrates one marginal drift against the
    // conditional survivor of the other.  It is independent of the corner
    // derivative algebra and serves as an independent check on the data-window
    // normaliser, which is cached by parameter cell below.  An exact event is
    // retained only with an exact denominator; selecting the numeric
    // denominator makes the dispatcher retry the whole trial numerically.
    // Both routes use the same positive-orthant normaliser, so q cancels
    // exactly in their log ratio.
    // The numeric integral is an independent check, but its marginal scan is
    // deliberately finite.  When the surviving mass sits beyond that scan
    // (high positive rho with very asymmetric drift means), it can
    // underestimate Z by hundreds of log units even though the exact
    // rectangle result is well resolved.  Conversely, the exact corner
    // subtraction can under-resolve central rare-positive-orthant cases that
    // the numeric route handles well.  Both known failures lose positive
    // mass, so retain the larger valid estimate rather than unconditionally
    // replacing the exact denominator with the numeric one.
    auto difference_for_mode = [&](bool use_numeric,
                                   BAwLCorrMomentStatus* mode_status) {
      BAwLCorrMomentStatus combined = BAwLCorrMomentStatus::ok;
      auto survivor = [&](double t) {
        BAwLCorrMomentStatus local = BAwLCorrMomentStatus::ok;
        const double value = bawl_corr_log_component_survival(
            s, pair, j, t, cols, ctx, use_numeric, &local);
        if (local == BAwLCorrMomentStatus::invalid) combined = local;
        else if (local == BAwLCorrMomentStatus::unstable &&
                 combined != BAwLCorrMomentStatus::invalid) combined = local;
        else if (local == BAwLCorrMomentStatus::zero_mass &&
                 combined == BAwLCorrMomentStatus::ok) combined = local;
        return value;
      };
      const double a = (lo == 0.0) ? 0.0 : survivor(lo);
      if (hi == R_PosInf) {
        *mode_status = combined;
        return a;
      }
      const double b = survivor(hi);
      *mode_status = combined;
      if (a == R_NegInf) return R_NegInf;
      return log_diff_exp(a, b);
    };

    auto valid_difference = [](double log_z, BAwLCorrMomentStatus status) {
      return R_FINITE(log_z) && status != BAwLCorrMomentStatus::unstable &&
        status != BAwLCorrMomentStatus::invalid;
    };
    if (numeric) {
      BAwLCorrMomentStatus numeric_status = BAwLCorrMomentStatus::ok;
      const double numeric_log_z = difference_for_mode(true, &numeric_status);
      const bool numeric_valid = valid_difference(numeric_log_z, numeric_status);
      if (numeric_valid) {
        *st = numeric_status;
        absorb_status(*st);
        return numeric_log_z;
      }
      *st = numeric_status == BAwLCorrMomentStatus::invalid
        ? BAwLCorrMomentStatus::invalid : BAwLCorrMomentStatus::unstable;
      absorb_status(*st);
      return R_NegInf;
    }

    BAwLCorrMomentStatus exact_status = BAwLCorrMomentStatus::ok;
    const double exact_log_z = difference_for_mode(false, &exact_status);
    const bool exact_valid = valid_difference(exact_log_z, exact_status);
    // Corner cancellation is an absolute-precision problem.  A resolved
    // normaliser above 1e-8 is far from that regime and does not need the
    // substantially more expensive 64-node independent cross-check.
    constexpr double kLogZCrosscheckThreshold = -18.420680743952367; // log(1e-8)
    if (exact_valid && exact_log_z > kLogZCrosscheckThreshold) {
      *st = exact_status;
      absorb_status(*st);
      return exact_log_z;
    }

    BAwLCorrMomentStatus numeric_status = BAwLCorrMomentStatus::ok;
    const double numeric_log_z = difference_for_mode(true, &numeric_status);
    const bool numeric_valid = valid_difference(numeric_log_z, numeric_status);
    if (exact_valid || numeric_valid) {
      constexpr double kLogZAgreementTol = 1e-8;
      const bool take_exact = exact_valid &&
        (!numeric_valid || numeric_log_z <= exact_log_z + kLogZAgreementTol);
      *st = take_exact ? exact_status : numeric_status;
      absorb_status(*st);
      // Do not combine an exact event density with a numeric denominator.
      // Ask the existing outer dispatcher to repeat the whole trial on the
      // numeric-pair route, where an unrepresentable event is conservatively
      // floored instead of becoming an artificial finite likelihood spike.
      denominator_requires_numeric = !take_exact;
      return take_exact ? exact_log_z : numeric_log_z;
    }
    *st = (exact_status == BAwLCorrMomentStatus::invalid ||
           numeric_status == BAwLCorrMomentStatus::invalid)
      ? BAwLCorrMomentStatus::invalid : BAwLCorrMomentStatus::unstable;
    absorb_status(*st);
    return R_NegInf;
  };

  double log_value = R_NegInf;
  if (R_FINITE(rt) && rt > 0.0) {
    log_value = finite_cause(rt, known ? known_row : -1, &st);
  } else if (rt == R_NegInf) {
    if (known) log_value = interval_cause(LT, LC, known_row, &st);
    else log_value = survival_difference(LT, LC, &st);
  } else if (rt == R_PosInf) {
    if (known) log_value = interval_cause(UC, UT, known_row, &st);
    else log_value = survival_difference(UC, UT, &st);
  } else if (Rcpp::NumericVector::is_na(rt)) {
    if (known) {
      log_value = log_sum_exp(interval_cause(LT, LC, known_row, &st),
                              interval_cause(UC, UT, known_row, &st));
    } else {
      log_value = log_sum_exp(survival_difference(LT, LC, &st),
                              survival_difference(UC, UT, &st));
    }
  }
  if (st == BAwLCorrMomentStatus::unstable || st == BAwLCorrMomentStatus::invalid) {
    out.status = st;
    return out;
  }

  if (s.has_trunc_trial[static_cast<size_t>(j)]) {
    // Z depends only on the parameter cell and the data-constant truncation
    // window, never on rt, so unique trials sharing a cell share Z.  Only
    // pure pair trials are cacheable: a singleton survivor's parameters
    // would otherwise need to join the key.
    const bool z_cacheable = !numeric && layout.n_active == 2;
    double zkey[15] = {0.0};
    const double* z_hit = nullptr;
    if (z_cacheable) {
      const double* t0c = cols[emc2col::bawl::t0];
      const double* Ac = cols[emc2col::bawl::A];
      const double* Bc = cols[emc2col::bawl::B];
      zkey[0] = t0c[pair.row1]; zkey[1] = Ac[pair.row1]; zkey[2] = Bc[pair.row1];
      zkey[3] = bawl_corr_row_k(cols, pair.row1, ctx);
      zkey[4] = pair.mu1; zkey[5] = pair.sd1;
      zkey[6] = t0c[pair.row2]; zkey[7] = Ac[pair.row2]; zkey[8] = Bc[pair.row2];
      zkey[9] = bawl_corr_row_k(cols, pair.row2, ctx);
      zkey[10] = pair.mu2; zkey[11] = pair.sd2;
      zkey[12] = pair.rho; zkey[13] = LT; zkey[14] = UT;
      for (const auto& e : s.z_cache) {
        if (std::memcmp(e.key, zkey, sizeof(zkey)) == 0) {
          z_hit = &e.log_z;
          break;
        }
      }
    }
    if (z_hit != nullptr) {
      log_value -= *z_hit;
    } else {
      BAwLCorrMomentStatus zst = BAwLCorrMomentStatus::ok;
      const double log_z = stable_survival_difference(LT, UT, &zst);
      if (zst == BAwLCorrMomentStatus::unstable || zst == BAwLCorrMomentStatus::invalid ||
          !R_FINITE(log_z)) {
        out.status = zst == BAwLCorrMomentStatus::ok
          ? BAwLCorrMomentStatus::unstable : zst;
        return out;
      }
      if (denominator_requires_numeric) {
        out.status = BAwLCorrMomentStatus::unstable;
        return out;
      }
      if (z_cacheable && s.z_cache.size() < 64) {
        CorrDriftSharedState::ZCacheEntry e;
        std::memcpy(e.key, zkey, sizeof(zkey));
        e.log_z = log_z;
        s.z_cache.push_back(e);
      }
      log_value -= log_z;
    }
  }
  // Contaminant mixture; see src/contaminant_mixture.h.
  if (s.pc_col >= 0 || (s.pg_col >= 0 && s.guess.active())) {
    const double pC = (s.pc_col >= 0) ? cols[static_cast<size_t>(s.pc_col)][start] : 0.0;
    const double pG = (s.pg_col >= 0 && s.guess.active())
                        ? cols[static_cast<size_t>(s.pg_col)][start] : 0.0;
    log_value = mix_contaminants_rt(log_value, pC, pG, s.guess, rt, known);
  }
  out.log_likelihood = (!R_FINITE(log_value) || log_value < min_ll)
    ? min_ll : log_value;
  out.status = BAwLCorrMomentStatus::ok;
  out.used_numeric = numeric;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int r = start + k;
    if (bawl_corr_row_active(s, r) &&
        cols[emc2col::bawl::A][r] <= BAWL_A_EPS) out.point_start = true;
  }
  (void)winner;
  return out;
}

// Correlated BAwL likelihood using batched Gauss-Hermite integration.
// The scan is recentered for narrow or tail-peaked trial integrands.
// Conditional on the shared factor, accumulators remain independent.
// Raw and general evaluators supply node integrands for the common and
// exceptional layouts.  Truncation is normalized as
//   log int p(data | z) phi(z) dz - log int Z(z) phi(z) dz.
double c_log_likelihood_corr_drift(
    CorrDriftSharedState& cshared,
    Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector isok,
    int n_lR,
    void* model_context_for_funcs,
    bool all_finite_trials,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t,
    RaceSharedState* shared,
    NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize) {
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  if (ctx == nullptr || !ctx->corr_drift_active || ctx->corr_drift_rho_index < 0 ||
      !cshared.valid || cshared.n_trials != n_trials || cshared.n_lR != n_lR) {
    Rcpp::stop("c_log_likelihood_corr_drift: invalid correlated BAwL context.");
  }

  const int n_unique = cshared.n_unique;
  const int n_out = cshared.n_out;
  if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
    Rcpp::stop("c_log_likelihood_corr_drift: trial_ll_out size mismatch.");
  }
  const int n_par = cshared.n_par;
  if (n_par > 64 || n_lR > 64) {
    Rcpp::stop("c_log_likelihood_corr_drift: at most 64 parameter columns and 64 accumulators are supported.");
  }

  // Role-mapped effective loadings; roles were resolved once in the shared
  // state, so per particle this is one pass over the raw rho column.
  const double* rho_col_ptr = cshared.cols[static_cast<size_t>(cshared.rho_col)];
  std::vector<double>& effective_rho = cshared.effective_rho;
  for (int r = 0; r < n_trials; ++r) {
    const double rho_cell = rho_col_ptr[r];
    effective_rho[static_cast<size_t>(r)] =
      cshared.role_correct[static_cast<size_t>(r)] ? std::fabs(rho_cell)
                                                   : rho_cell;
  }
  // Cell-keyed memoization is only valid within one particle: the ParamTable
  // base storage the column pointers view is refilled in place per particle.
  cshared.qab_cache.clear();
  cshared.z_cache.clear();

  const bool count_routes = bawl_corr_counters_enabled();
  bawl_corr_counters_active() = count_routes;
  BAwLCorrCounters& route_counters = bawl_corr_counters();

  // If all active loadings are zero, use the exact existing path.  Apart from
  // being cheaper, this keeps rho=0 bit-for-bit equivalent to independent
  // BAwL.  This is the only route below that consumes a materialised matrix
  // besides the generic-clock node evaluator.
  bool any_loading = false;
  for (int j = 0; j < cshared.n_unique && !any_loading; ++j) {
    const int start = j * n_lR;
    const int n_lR_j = cshared.has_RACE ? cshared.race_nacc[start] : n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int r = start + k;
      if (cshared.has_RACE && !cshared.race_mask[r]) continue;
      if (!isok[r]) continue;
      const double rho = effective_rho[static_cast<size_t>(r)];
      if (!R_FINITE(rho) || std::fabs(rho) > 1e-14) {
        any_loading = true;
        break;
      }
    }
  }
  if (!any_loading) {
    if (count_routes)
      route_counters.ordinary_zero_rho_trials += n_unique;
    return c_log_likelihood_race(materialize(), dadm, pdf1, cdf1, n_trials,
                                  winner, expand, min_ll, isok, n_lR,
                                  model_context_for_funcs, all_finite_trials,
                                  model_dfun_raw, model_pfun_raw, logS_at_t,
                                  shared, trial_ll_out, true);
  }

  const DriftFactorColumns& factor_cols = cshared.factor_cols;
  const bool has_RACE_col = cshared.has_RACE;
  const Rcpp::IntegerVector& RACE = cshared.race_nacc;
  const Rcpp::LogicalVector& RACE_mask = cshared.race_mask;
  const Rcpp::NumericVector& LT = cshared.LT;
  const Rcpp::NumericVector& UT = cshared.UT;
  const bool has_truncation = !all_finite_trials && cshared.has_truncation;
  const std::vector<unsigned char>& has_trunc_trial = cshared.has_trunc_trial;

  // Raw batched conditional evaluator eligibility for the common design:
  // all-finite RTs, no truncation, no global-kill/timer/nogo machinery
  // outside the kernels.  Everything else falls back to the generic race
  // evaluator below.  The fused evaluator is built on BAwLPreparedRow, so it
  // is available only to models with an affine-in-drift survivor.
  const bool use_fast_node_eval =
      !cshared.generic_only &&
      all_finite_trials && !has_truncation &&
      model_dfun_raw != nullptr && model_pfun_raw != nullptr &&
      !ctx->has_global_kill() && !ctx->kill_active && !ctx->gng &&
      ctx->time_code == -1 && ctx->nogo_code == -1 &&
      (!dadm.containsElementNamed("RACE") || has_RACE_col);

  // Classify each particle's trial layout.
  const bool no_clock_model = !ctx->gng && !ctx->kill_active &&
    !ctx->has_global_kill() && ctx->time_code == -1 && ctx->nogo_code == -1;
  // The exact pair component is data-fixed and pointer-based, so it remains
  // eligible even when the caller disables the all-finite GH hint (that hint
  // controls only the fused node evaluator).  This also lets nonfinite and
  // truncated no-clock trials use the same component path.
  corr_drift_classify_particle(cshared, isok, no_clock_model);
  if (count_routes) {
    for (int j = 0; j < n_unique; ++j) {
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (L.n_loaded == 0) ++route_counters.loaded_dimension_0;
      else if (L.n_loaded == 1) ++route_counters.loaded_dimension_1;
      else if (L.n_loaded == 2) ++route_counters.loaded_dimension_2;
      else ++route_counters.loaded_dimension_3plus;
      if (L.route == BAwLCorrRoute::ordinary) {
        if (L.n_loaded == 0) ++route_counters.ordinary_zero_rho_trials;
        else ++route_counters.ordinary_single_loaded_trials;
      }
    }
  }

  // Exact pair trials are evaluated before the shared-factor fallback.  The
  // arrays are also masks for the node evaluator, so mixed particles retain
  // per-trial routing rather than promoting an exact pair because another
  // trial needs GH.  `R` is data-fixed and is only used to distinguish a
  // known response from an unknown winner/censoring observation.
  std::vector<unsigned char> exact_done(static_cast<size_t>(n_unique), 0);
  std::vector<double> exact_ll(static_cast<size_t>(n_unique), min_ll);

  // A trial with at most one loaded racer has no cross-racer dependence left:
  // the shared factor can be integrated out analytically and the trial is
  // exactly the ordinary independent BAwL race.  Evaluate those trials once
  // through the existing race likelihood, then keep them out of every GH
  // node.  This is deliberately done at the unique-trial level so the
  // ordinary path retains its complete handling of omissions, censoring,
  // truncation, Erlang clocks, contaminants, RACE masks, and expansion.
  //
  // The ordinary call receives a mask containing only ordinary trials.  Rows
  // belonging to correlated trials are invalid for that call and their
  // returned values are ignored; this avoids constructing a second DataFrame
  // or changing the data-fixed expansion/index bookkeeping.
  bool any_ordinary_trial = false;
  Rcpp::LogicalVector ordinary_isok(n_trials, false);
  for (int j = 0; j < n_unique; ++j) {
    const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
    if (L.route != BAwLCorrRoute::ordinary) continue;
    any_ordinary_trial = true;
    const int start = j * n_lR;
    const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!has_RACE_col || RACE_mask[row]) ordinary_isok[row] = isok[row];
    }
  }
  if (any_ordinary_trial) {
    Rcpp::NumericMatrix pars_ordinary = materialize();
    Rcpp::NumericVector ordinary_ll(n_unique, min_ll);
    Rcpp::IntegerVector ordinary_expand(0);
    // c_log_likelihood_race updates a few context flags while inspecting the
    // active rows.  Isolate those updates from the correlated remainder,
    // whose GH evaluator must retain the original particle context.
    ContextForRaceModels ordinary_ctx = *ctx;
    c_log_likelihood_race(
        pars_ordinary, dadm, pdf1, cdf1, n_trials, winner, ordinary_expand,
        min_ll, ordinary_isok, n_lR, &ordinary_ctx, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, &ordinary_ll, true);
    for (int j = 0; j < n_unique; ++j) {
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (L.route != BAwLCorrRoute::ordinary) continue;
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = ordinary_ll[j];
    }
  }

  const Rcpp::IntegerVector response = dadm["R"];
  for (int j = 0; j < n_unique; ++j) {
    BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
    if (L.route != BAwLCorrRoute::exact_pair) continue;
    BAwLCorrExactTrialResult exact = bawl_corr_exact_trial_loglik(
        cshared, L, j, cshared.cols.data(), ctx, response, winner, min_ll, false);
    if (exact.status == BAwLCorrMomentStatus::unstable) {
      L.route = BAwLCorrRoute::numeric_pair;
      exact = bawl_corr_exact_trial_loglik(
          cshared, L, j, cshared.cols.data(), ctx, response, winner, min_ll, true);
    }
    if (exact.status == BAwLCorrMomentStatus::ok) {
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = exact.log_likelihood;
      if (count_routes) {
        if (L.route == BAwLCorrRoute::numeric_pair)
          ++route_counters.numeric_pair_trials;
        else {
          ++route_counters.exact_pair_trials;
          if (L.winner_loaded) ++route_counters.exact_pair_pair_winner_trials;
          else ++route_counters.exact_pair_independent_winner_trials;
          if (exact.point_start) ++route_counters.exact_pair_point_start_trials;
        }
      }
    } else {
      // Both the exact rectangle and its numeric fallback found this trial's
      // geometry too degenerate to trust (e.g. a near-zero residual SD driven
      // by a runaway sv/rho combination).  That verdict is a property of the
      // parameter point, not of the kernel: routing it into the generic GH
      // quadrature would silently hand a numerically pathological integrand
      // to code with no equivalent stability check, which can return a large
      // finite value instead of failing.  Floor it directly instead of
      // giving the sampler a route to exploit that hole.
      exact_done[static_cast<size_t>(j)] = 1;
      exact_ll[static_cast<size_t>(j)] = min_ll;
      if (count_routes) ++route_counters.unstable_pair_floored_trials;
    }
  }

  int n_fallback_trials = 0;
  for (int j = 0; j < n_unique; ++j)
    if (!exact_done[static_cast<size_t>(j)]) ++n_fallback_trials;
  if (n_fallback_trials == 0) {
    // The exact component has already applied truncation, contamination, and
    // the final floor.  Do not enter the GH machinery when every unique trial
    // is a pair: this is the common canonical path and should pay no fallback
    // node or prepared-geometry work.
    double total_ll = 0.0;
    if (expand.length() > 0) {
      for (int i = 0; i < n_out; ++i) {
        const double value = exact_ll[static_cast<size_t>(expand[i] - 1)];
        if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
        total_ll += value;
      }
    } else {
      for (int j = 0; j < n_unique; ++j) {
        const double value = exact_ll[static_cast<size_t>(j)];
        if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
        total_ll += value;
      }
    }
    return total_ll;
  }

  const bool joint_posdrift = ctx->use_posdrift;
  // Numerator/denominator integrals per unique trial; expanded rows share
  // their unique trial's value and are written out at the end.
  std::vector<double> log_num(static_cast<size_t>(n_unique), R_NegInf);
  // Joint positive-drift conditioning needs a denominator for every trial:
  // integral phi(z) prod_k q_k(z) Z+(z) dz.  Without joint conditioning the
  // untruncated denominator is exactly one and only truncated trials need it.
  std::vector<double> log_den(static_cast<size_t>(n_unique),
                              joint_posdrift ? R_NegInf : 0.0);
  if (has_truncation && !joint_posdrift) {
    for (int j = 0; j < n_unique; ++j) {
      if (has_trunc_trial[static_cast<size_t>(j)])
        log_den[static_cast<size_t>(j)] = R_NegInf;
    }
  }

  // Positive-drift denominators without truncation have a closed form: the
  // Phi argument of racer k is a_k + b_k z, so
  //   int phi(z) prod_k Phi(a_k + b_k z) dz
  //     = Phi_K(a_k / sqrt(1 + b_k^2); rho_kl = b_k b_l / sqrt((1+b_k^2)(1+b_l^2))),
  // the one-factor orthant probability.  K <= 2 covers the dominant designs
  // exactly (norm_cdf_2d); larger K and truncated trials keep the quadrature
  // denominator below.  At high |rho| this integrand is phi against a steep
  // sigmoid, which fixed and moment-matched Hermite rules resolve worst, so
  // the closed form is an accuracy fix, not just a shortcut.
  // The positivity dimension is the layout's loaded-row count: exact-zero
  // loadings contribute the same constant q to numerator and denominator and
  // are cancelled algebraically on both sides, so an independent PM racer
  // neither raises the closed-form dimension nor forces the quadrature
  // denominator sweep.
  std::vector<unsigned char> den_by_quadrature(static_cast<size_t>(n_unique), 0);
  {
    const double* v0c = cshared.cols[static_cast<size_t>(factor_cols.v)];
    const double* sv0c = cshared.cols[static_cast<size_t>(factor_cols.sv)];
    for (int j = 0; j < n_unique; ++j) {
      if (has_trunc_trial[static_cast<size_t>(j)]) {
        if (exact_done[static_cast<size_t>(j)]) continue;
        den_by_quadrature[static_cast<size_t>(j)] = 1;
        continue;
      }
      if (!joint_posdrift) continue;  // denominator is exactly one
      const BAwLCorrTrialLayout& Lj = cshared.layout[static_cast<size_t>(j)];
      if (exact_done[static_cast<size_t>(j)]) continue;
      if (Lj.any_bad_row) {
        log_den[static_cast<size_t>(j)] = R_NegInf;
        continue;
      }
      if (Lj.n_loaded > 2) {
        den_by_quadrature[static_cast<size_t>(j)] = 1;
        continue;
      }
      double h[2] = {0.0, 0.0};
      double b[2] = {0.0, 0.0};
      bool bad = false;
      for (int m = 0; m < Lj.n_loaded; ++m) {
        const int row = Lj.loaded_row[m];
        const double rho = effective_rho[static_cast<size_t>(row)];
        const double sv = sv0c[row];
        const double svq =
          sv * std::fmax(std::sqrt(std::fmax(0.0, 1.0 - std::fabs(rho))), 1e-12);
        const double mu = v0c[row];
        if (!(svq > 0.0) || !R_FINITE(mu)) { bad = true; break; }
        const double slope =
          ((rho < 0.0) ? -1.0 : 1.0) * sv * std::sqrt(std::fabs(rho));
        const double bk = slope / svq;
        const double hk = (mu / svq) / std::sqrt(1.0 + bk * bk);
        if (!R_FINITE(bk) || !R_FINITE(hk)) { bad = true; break; }
        h[m] = hk;
        b[m] = bk;
      }
      if (bad) {
        log_den[static_cast<size_t>(j)] = R_NegInf;
        continue;
      }
      if (Lj.n_loaded == 0) {
        log_den[static_cast<size_t>(j)] = 0.0;
      } else if (Lj.n_loaded == 1) {
        log_den[static_cast<size_t>(j)] = pnorm_log_direct(h[0], true);
      } else {
        const double rho12 =
          b[0] * b[1] / std::sqrt((1.0 + b[0] * b[0]) * (1.0 + b[1] * b[1]));
        const double p = norm_cdf_2d_hybrid(h[0], h[1], rho12);
        log_den[static_cast<size_t>(j)] =
          (R_FINITE(p) && p > 0.0) ? std::log(std::fmin(p, 1.0)) : R_NegInf;
      }
    }
  }

  // The truncation/positivity denominator is independent of the observed RT
  // and response.  Compression therefore leaves many unique trials with the
  // same denominator cell (the mapped parameter rows, active-racer layout,
  // and LT/UT window are identical).  Pick one representative per exact cell
  // and evaluate the GH denominator only for those representatives.  The
  // comparison is deliberately on mapped natural-scale values: it is safe for
  // arbitrary user design functions and also catches equal cells produced by
  // different design-matrix rows.
  std::vector<int> den_representative(static_cast<size_t>(n_unique), -1);
  std::vector<unsigned char> den_eval_mask(static_cast<size_t>(n_unique), 0);
  auto same_denominator_cell = [&](int a, int b) {
    const int sa = a * n_lR;
    const int sb = b * n_lR;
    if (LT[sa] != LT[sb] || UT[sa] != UT[sb]) return false;
    const int na = has_RACE_col ? RACE[sa] : n_lR;
    const int nb = has_RACE_col ? RACE[sb] : n_lR;
    if (na != nb) return false;
    for (int k = 0; k < na; ++k) {
      const int ra = sa + k;
      const int rb = sb + k;
      if (has_RACE_col && RACE_mask[ra] != RACE_mask[rb]) return false;
      if (isok[ra] != isok[rb] ||
          effective_rho[static_cast<size_t>(ra)] !=
            effective_rho[static_cast<size_t>(rb)]) return false;
      for (int c = 0; c < n_par; ++c) {
        if (cshared.cols[static_cast<size_t>(c)][ra] !=
            cshared.cols[static_cast<size_t>(c)][rb]) return false;
      }
    }
    return true;
  };
  int n_den_cells = 0;
  for (int j = 0; j < n_unique; ++j) {
    if (!den_by_quadrature[static_cast<size_t>(j)]) continue;
    int representative = j;
    for (int r = 0; r < j; ++r) {
      if (!den_eval_mask[static_cast<size_t>(r)]) continue;
      if (same_denominator_cell(j, r)) {
        representative = r;
        break;
      }
    }
    den_representative[static_cast<size_t>(j)] = representative;
    if (representative == j) {
      den_eval_mask[static_cast<size_t>(j)] = 1;
      ++n_den_cells;
    }
  }

  // Positivity reweighting is evaluated on the node-shifted Gaussian means
  // and residual SDs.  RACE-inactive rows do not participate in the product,
  // and zero-loading rows are cancelled: their constant q appears in both
  // the numerator and denominator integrands, so it is omitted from both.
  auto log_positive_trial = [&](const Rcpp::NumericMatrix& pars_q,
                                const Rcpp::LogicalVector& isok_q,
                                int j) -> double {
    const int start = j * n_lR;
    const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
    double out = 0.0;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (has_RACE_col && !RACE_mask[row]) continue;
      if (!isok_q[row]) return R_NegInf;
      const double sd = pars_q(row, factor_cols.sv);
      const double mu = pars_q(row, factor_cols.v);
      const double rho = effective_rho[static_cast<size_t>(row)];
      if (!(sd > 0.0) || !R_FINITE(mu)) return R_NegInf;
      if (std::fabs(rho) <= DRIFT_FACTOR_RHO_EPS) continue;  // cancelled constant
      const double log_q = drift_factor_log_positive_row(mu, sd);
      if (!R_FINITE(log_q)) return R_NegInf;
      out += log_q;
    }
    return out;
  };

  // Both node evaluators report per-unique-trial integrands; expanded output
  // rows share their unique trial's integral, so expansion happens once at
  // the end rather than per node.  j_to_i lets the generic evaluator read one
  // representative expanded row per unique trial.
  const std::vector<int>& j_to_i = cshared.j_to_i;

  // The generic-clock node evaluator is the only remaining consumer of a
  // materialised parameter matrix; produce it once per particle and only
  // when that route is actually selected.
  Rcpp::NumericMatrix pars_generic;
  if (!use_fast_node_eval) pars_generic = materialize();
  // Per-trial log(1 - pG) for the generic node evaluator, and the matching
  // suppression of pGuess inside the nested race call -- same split as the fast
  // path above: the constant factor is integrated, the guess mass is not.
  std::vector<double> generic_log1m_pg(static_cast<size_t>(n_unique), 0.0);
  if (!use_fast_node_eval && cshared.pg_col >= 0 && cshared.guess.active()) {
    for (int j = 0; j < n_unique; ++j) {
      const double pG = pars_generic(j * n_lR, cshared.pg_col);
      if (pG != 0.0) generic_log1m_pg[static_cast<size_t>(j)] = log1m(pG);
    }
    for (int r = 0; r < n_trials; ++r) pars_generic(r, cshared.pg_col) = 0.0;
  }

  // Step 4: the no-clock route consumes prepared BAwL geometry directly;
  // prepared conditional survivors retain scalar natural-branch semantics,
  // including the saturated-CDF boundary and separate q accounting.
  const double* fast_rt = nullptr;
  const double* fast_v0 = nullptr;
  const double* fast_sv0 = nullptr;
  const double* fast_B0 = nullptr;
  const double* fast_A0 = nullptr;
  const double* fast_t00 = nullptr;
  const double* fast_k0 = nullptr;
  std::vector<double> fast_res, fast_log1m_pc;
  std::vector<double> fast_const_ll;
  std::vector<BAwLPreparedRow> fast_prepared;
  std::vector<int> fast_winner, fast_loser, fast_const_winner, fast_const_loser, fast_ok;
  std::vector<int> fine_winner, fine_loser;
  bool use_fine_masks = false;
  bool has_fast_const = false;
  if (use_fast_node_eval) {
    fast_rt = cshared.rt.begin();
    fast_v0 = cshared.cols[static_cast<size_t>(emc2col::bawl::v)];
    fast_sv0 = cshared.cols[static_cast<size_t>(emc2col::bawl::sv)];
    fast_B0 = cshared.cols[static_cast<size_t>(emc2col::bawl::B)];
    fast_A0 = cshared.cols[static_cast<size_t>(emc2col::bawl::A)];
    fast_t00 = cshared.cols[static_cast<size_t>(emc2col::bawl::t0)];
    fast_k0 = ctx->bawl_k_fixed_zero
      ? nullptr : cshared.cols[static_cast<size_t>(emc2col::bawl::k)];
    fast_res.assign(static_cast<size_t>(n_trials), R_NegInf);
    fast_prepared.assign(static_cast<size_t>(n_trials), BAwLPreparedRow());
    fast_winner.resize(static_cast<size_t>(n_trials));
    fast_loser.resize(static_cast<size_t>(n_trials));
    fast_const_winner.assign(static_cast<size_t>(n_trials), 0);
    fast_const_loser.assign(static_cast<size_t>(n_trials), 0);
    fast_ok.resize(static_cast<size_t>(n_trials));
    for (int r = 0; r < n_trials; ++r) {
      if (exact_done[static_cast<size_t>(r / n_lR)]) {
        fast_ok[static_cast<size_t>(r)] = 0;
        fast_winner[static_cast<size_t>(r)] = 0;
        fast_loser[static_cast<size_t>(r)] = 0;
        fast_const_winner[static_cast<size_t>(r)] = 0;
        fast_const_loser[static_cast<size_t>(r)] = 0;
        continue;
      }
      const double rho = effective_rho[static_cast<size_t>(r)];
      const bool ok_r = isok[r] && R_FINITE(rho) && std::fabs(rho) <= 1.0;
      fast_ok[static_cast<size_t>(r)] = ok_r ? 1 : 0;
      const double magnitude = std::fabs(rho);
      const bool active = !has_RACE_col || RACE_mask[r];
      // A valid exact-zero loading is independent of the shared factor. Its
      // race contribution is evaluated once below and added to every node;
      // only nonzero-loading rows remain in the node masks.
      const bool independent = ok_r && magnitude <= DRIFT_FACTOR_RHO_EPS;
      if (active && independent) has_fast_const = true;
      fast_winner[static_cast<size_t>(r)] =
        (active && winner[r] && !independent) ? 1 : 0;
      fast_loser[static_cast<size_t>(r)] =
        (active && n_lR > 1 && !winner[r] && !independent) ? 1 : 0;
      fast_const_winner[static_cast<size_t>(r)] =
        (active && winner[r] && independent) ? 1 : 0;
      fast_const_loser[static_cast<size_t>(r)] =
        (active && n_lR > 1 && !winner[r] && independent) ? 1 : 0;

      if (active && ok_r && !independent) {
        // All fixed-time LBA/leak work is prepared once per particle.  The
        // node loop below only updates v_q and evaluates the normal endpoints
        // through BAwLPreparedRow; it never rebuilds tau, leak factors, or
        // Jacobians and never stages a conditional parameter table.
        const double kval = ctx->bawl_k_fixed_zero ? 0.0 : fast_k0[r];
        const BAwLTimeGeometry g = bawl_time_geometry(
            fast_rt[r], fast_t00[r], fast_A0[r], fast_B0[r] + fast_A0[r], kval);
        const DriftFactorLoading l = drift_factor_loading(fast_sv0[r], rho);
        fast_prepared[static_cast<size_t>(r)] =
            bawl_prepare_row(g, fast_v0[r], l.slope, l.sv_res);
        if (count_routes) ++route_counters.prepared_rows;
      }
    }

    fast_const_ll.assign(static_cast<size_t>(n_unique), 0.0);
    if (has_fast_const) {
      // Evaluate exact-zero-loading race terms once.  These are ordinary
      // singleton likelihood factors, so their own positive-drift
      // normalisers remain in the singleton kernel and are not part of the
      // shared-factor denominator.
      model_dfun_raw(fast_rt, cshared.cols.data(), n_trials,
                   fast_const_winner.data(), fast_ok.data(), fast_res.data(),
                   R_NegInf, ctx);
      for (int j = 0; j < n_unique; ++j) {
        const int start = j * n_lR;
        const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
        double value = 0.0;
        for (int k = 0; k < n_lR_j; ++k) {
          const int row = start + k;
          if (has_RACE_col && !RACE_mask[row]) continue;
          if (fast_const_winner[static_cast<size_t>(row)])
            value += fast_res[static_cast<size_t>(row)];
        }
        fast_const_ll[static_cast<size_t>(j)] = value;
      }
      if (n_lR > 1) {
        model_pfun_raw(fast_rt, cshared.cols.data(), n_trials,
                     fast_const_loser.data(), fast_ok.data(), fast_res.data(),
                     R_NegInf, ctx);
        for (int j = 0; j < n_unique; ++j) {
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          for (int k = 0; k < n_lR_j; ++k) {
            const int row = start + k;
            if (has_RACE_col && !RACE_mask[row]) continue;
            if (fast_const_loser[static_cast<size_t>(row)])
              fast_const_ll[static_cast<size_t>(j)] += fast_res[static_cast<size_t>(row)];
          }
        }
      }
    }

    // The generic evaluator inherits the contaminant shift from the race
    // path; with all-finite RTs that is a per-trial log(1 - pC) added to the
    // numerator integrand.
    //
    // pGuess is split in two here, because unlike pC's log(1 - pC) the guess is
    // NOT a constant factor and must not be integrated over z.  Only its
    // (1 - pG) down-weight -- which IS a constant factor, and therefore rides
    // through both the quadrature and the joint-positivity division unchanged --
    // goes into the integrand; the flat guess MASS is added once to the finished
    // trial likelihood in the accumulation loop below.  The split is exact.
    fast_log1m_pc.assign(static_cast<size_t>(n_unique), 0.0);
    if (cshared.pc_col >= 0) {
      const double* pc_ptr = cshared.cols[static_cast<size_t>(cshared.pc_col)];
      for (int j = 0; j < n_unique; ++j) {
        const double pC = pc_ptr[j * n_lR];
        if (pC != 0.0) fast_log1m_pc[static_cast<size_t>(j)] = std::log1p(-pC);
      }
    }
    if (cshared.pg_col >= 0 && cshared.guess.active()) {
      const double* pg_ptr = cshared.cols[static_cast<size_t>(cshared.pg_col)];
      for (int j = 0; j < n_unique; ++j) {
        const double pG = pg_ptr[j * n_lR];
        if (pG != 0.0) fast_log1m_pc[static_cast<size_t>(j)] += log1m(pG);
      }
    }
  }

  // Evaluators take one latent value per unique trial so the refined pass
  // can center each trial's rule on its own integrand; a shared node is just
  // a constant z_by_trial.
  auto evaluate_node_fast = [&](const double* z_by_trial, bool need_num,
                                std::vector<double>& node_num,
                                std::vector<double>& node_den) {
    node_num.assign(static_cast<size_t>(n_unique), R_NegInf);
    node_den.assign(static_cast<size_t>(n_unique),
                    joint_posdrift ? R_NegInf : 0.0);
    const std::vector<int>& winner_mask = use_fine_masks ? fine_winner : fast_winner;
    const std::vector<int>& loser_mask = use_fine_masks ? fine_loser : fast_loser;
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      const int start = j * n_lR;
      const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
      double s = has_fast_const ? fast_const_ll[static_cast<size_t>(j)] : 0.0;
      double lp = 0.0;
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (has_RACE_col && !RACE_mask[row]) continue;
        if (fast_const_winner[static_cast<size_t>(row)] ||
            fast_const_loser[static_cast<size_t>(row)]) continue;
        const double rho = effective_rho[static_cast<size_t>(row)];
        const bool loaded = fast_ok[static_cast<size_t>(row)] &&
          std::fabs(rho) > 1e-14;
        if (!loaded) {
          if (need_num && (winner_mask[static_cast<size_t>(row)] ||
                           loser_mask[static_cast<size_t>(row)]))
            s = R_NegInf;
          if (joint_posdrift) lp = R_NegInf;
          continue;
        }
        const BAwLPreparedRow& prepared = fast_prepared[static_cast<size_t>(row)];
        const double vq = prepared.v0 + prepared.slope * z_by_trial[j];
        if (joint_posdrift) {
          const double log_q = bawl_prepared_log_q(prepared, vq);
          if (!R_FINITE(log_q)) lp = R_NegInf;
          else if (lp != R_NegInf) lp += log_q;
        }
        if (need_num) {
          double log_component = R_NegInf;
          if (winner_mask[static_cast<size_t>(row)]) {
            log_component = bawl_prepared_log_pdf(prepared, vq);
          } else if (loser_mask[static_cast<size_t>(row)]) {
            log_component = bawl_prepared_log_survivor(
                prepared, vq, joint_posdrift);
          }
          s += log_component;
        }
      }
      if (need_num)
        node_num[static_cast<size_t>(j)] =
          // Loaded prepared components are already unnormalised: a winner
          // contributes f0 and a loser contributes q-F0.  The q product is
          // reserved for the positivity denominator, so adding lp here
          // would double-count every loaded loser.
          s + fast_log1m_pc[static_cast<size_t>(j)];
      if (joint_posdrift && den_eval_mask[static_cast<size_t>(j)])
        node_den[static_cast<size_t>(j)] = lp;
    }
  };

  // Generic conditional evaluator: one shared-race pass per node vector.
  // need_num/need_den let the refinement passes skip the parts they do not
  // consume (the denominator pass needs no race likelihood, only the
  // positivity product and truncation normalisers).
  auto evaluate_node = [&](const double* z_by_trial, bool need_num,
                           bool need_den, std::vector<double>& node_num,
                           std::vector<double>& node_den) {
    node_num.assign(static_cast<size_t>(n_unique), R_NegInf);
    node_den.assign(static_cast<size_t>(n_unique),
                    joint_posdrift ? R_NegInf : 0.0);
    Rcpp::NumericMatrix pars_q = Rcpp::clone(pars_generic);
    Rcpp::LogicalVector isok_q = Rcpp::clone(isok);

    for (int r = 0; r < n_trials; ++r) {
      const int j = r / n_lR;
      if (exact_done[static_cast<size_t>(j)] ||
          (!need_num && need_den && !den_eval_mask[static_cast<size_t>(j)]))
        isok_q[r] = false;
    }

    // Conditional factor model preserving each marginal sv; see
    // drift_factor.h for the construction and the SD floor.
    for (int r = 0; r < n_trials; ++r) {
      if (!isok_q[r]) continue;
      if (!drift_factor_apply_row(pars_q, factor_cols, r,
                                  effective_rho[static_cast<size_t>(r)],
                                  z_by_trial[r / n_lR])) {
        isok_q[r] = false;
      }
    }

    std::vector<double> log_pos(static_cast<size_t>(n_unique), 0.0);
    if (joint_posdrift) {
      for (int j = 0; j < n_unique; ++j) {
        if (!exact_done[static_cast<size_t>(j)])
          log_pos[static_cast<size_t>(j)] = log_positive_trial(pars_q, isok_q, j);
      }
    }

    ContextForRaceModels node_ctx = *ctx;
    if (need_num) {
      Rcpp::NumericVector node_ll(n_out);
      c_log_likelihood_race(pars_q, dadm, pdf1, cdf1, n_trials,
                            winner, expand, R_NegInf, isok_q, n_lR,
                            &node_ctx, all_finite_trials,
                            model_dfun_raw, model_pfun_raw, logS_at_t,
                            shared, &node_ll, false);

      // Expanded rows of one unique trial share the integrand; report once.
      for (int j = 0; j < n_unique; ++j) {
        if (exact_done[static_cast<size_t>(j)]) continue;
        const int i = j_to_i[static_cast<size_t>(j)];
        if (i < 0) continue;
        node_num[static_cast<size_t>(j)] =
          log_pos[static_cast<size_t>(j)] + node_ll[i] +
          generic_log1m_pg[static_cast<size_t>(j)];
      }
    }

    if (!need_den) return;
    if (has_truncation) {
      // The factor must remain inside the truncation normaliser:
      //   Z = int [S_race(LT | z) - S_race(UT | z)] phi(z) dz.
      // Do not replace this with an unconditional race normaliser.  However,
      // conditional on a GH node the accumulators are independent, so the
      // conditional race survivors can be evaluated in a batch.  The old
      // implementation repacked every trial into row-major storage and called
      // get_trunc_normaliser_rowmajor_cpp once per trial and node.
      GslIntegrationControls gsl_ctl = default_gsl_controls();
      gsl_ctl.try_qng_first_finite = true;
      gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
      gsl_ctl.rel_tol = 1e-5;
      static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
      GslWorkspacePtr& workspace = workspace_tls;

      std::vector<int> trunc_mask(static_cast<size_t>(n_unique), 0);
      std::vector<int> isok_q_int(static_cast<size_t>(n_trials), 0);
      double uniform_lt = 0.0;
      double uniform_ut = R_PosInf;
      bool have_trunc = false;
      bool uniform_lt_ok = true;
      bool uniform_ut_ok = true;
      for (int j = 0; j < n_unique; ++j) {
        const int start = j * n_lR;
        const double lt = LT[start];
        const double ut = UT[start];
        if (!has_trunc_trial[static_cast<size_t>(j)] ||
            !den_eval_mask[static_cast<size_t>(j)]) continue;
        trunc_mask[static_cast<size_t>(j)] = 1;
        if (!have_trunc) {
          uniform_lt = lt;
          uniform_ut = ut;
          have_trunc = true;
        } else {
          if (lt != uniform_lt) uniform_lt_ok = false;
          if (ut != uniform_ut) uniform_ut_ok = false;
        }
      }

      if (have_trunc) {
        for (int r = 0; r < n_trials; ++r) {
          isok_q_int[static_cast<size_t>(r)] = isok_q[r] ? 1 : 0;
        }

        std::vector<double> log_s_lt(static_cast<size_t>(n_unique), R_NegInf);
        std::vector<double> log_s_ut(static_cast<size_t>(n_unique), R_NegInf);
        std::array<const double*, 64> node_cols{};
        const double* pars_q_ptr = pars_q.begin();
        for (int c = 0; c < n_par; ++c) {
          node_cols[static_cast<size_t>(c)] =
            pars_q_ptr + static_cast<size_t>(c) * n_trials;
        }

        // Evaluate one conditional race survivor without constructing a
        // row-major trial buffer.  This is the fallback for variable-accumulator
        // RACE layouts, which the existing column-major batch callback cannot
        // represent because each trial can have a different active prefix.
        auto log_surv_trial = [&](double t, int j) -> double {
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          double log_s = 0.0;
          std::array<double, 64> par_row{};
          for (int k = 0; k < n_lR_j; ++k) {
            const int row = start + k;
            if (!isok_q[row] || (has_RACE_col && !RACE_mask[row])) return R_NegInf;
            for (int c = 0; c < n_par; ++c) par_row[static_cast<size_t>(c)] = pars_q(row, c);
            double cdf = cdf1(t, par_row.data(), &node_ctx);
            cdf = clamp_cdf01_race(cdf);
            const double log_s_k = safe_log1m_race(cdf);
            if (!emc2_isfinite(log_s_k)) return R_NegInf;
            log_s += log_s_k;
          }
          return log_s;
        };

        // Batch the common fixed-race case through the model-specific callback.
        // This is still conditional on the current GH node, so the shared
        // factor is fully represented in the final log-sum-exp below.
        auto log_surv_batch = [&](double t, std::vector<double>& out) {
          if (!has_RACE_col && logS_at_t != nullptr) {
            logS_at_t(t, node_cols.data(), n_trials, n_lR, n_par,
                      trunc_mask.data(), n_unique, isok_q_int.data(),
                      &node_ctx, out.data());
            return;
          }
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)]) {
              out[static_cast<size_t>(j)] = log_surv_trial(t, j);
            }
          }
        };

        if (uniform_lt_ok) {
          if (uniform_lt == 0.0) {
            for (int j = 0; j < n_unique; ++j) {
              if (trunc_mask[static_cast<size_t>(j)]) log_s_lt[static_cast<size_t>(j)] = 0.0;
            }
          } else {
            log_surv_batch(uniform_lt, log_s_lt);
          }
        } else {
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)]) {
              log_s_lt[static_cast<size_t>(j)] = log_surv_trial(LT[j * n_lR], j);
            }
          }
        }
        if (uniform_ut_ok && uniform_ut != R_PosInf) {
          log_surv_batch(uniform_ut, log_s_ut);
        } else if (!uniform_ut_ok) {
          for (int j = 0; j < n_unique; ++j) {
            if (trunc_mask[static_cast<size_t>(j)] && UT[j * n_lR] != R_PosInf) {
              log_s_ut[static_cast<size_t>(j)] = log_surv_trial(UT[j * n_lR], j);
            }
          }
        }

        std::vector<double> rowmajor(static_cast<size_t>(n_lR) * n_par);
        std::vector<int> ok_int(static_cast<size_t>(n_lR));
        const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
        for (int j = 0; j < n_unique; ++j) {
          if (!trunc_mask[static_cast<size_t>(j)]) continue;
          const int start = j * n_lR;
          const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
          const double lt = LT[start];
          const double ut = UT[start];
          double log_z = R_NegInf;
          // Under a global kill the batched survivors exclude the shared clock
          // (see get_trunc_normaliser_rowmajor_cpp), so always take the scalar
          // path, which corrects it or refuses.
          bool need_scalar_fallback = node_ctx.has_global_kill();

          if (!need_scalar_fallback && R_FINITE(log_s_lt[static_cast<size_t>(j)])) {
            if (ut == R_PosInf) {
              // For BAwL's defective tail, +Inf is part of the retained
              // [LT, Inf) window and must not be subtracted.
              log_z = log_s_lt[static_cast<size_t>(j)];
            } else {
              log_z = log_diff_exp(log_s_lt[static_cast<size_t>(j)],
                                   log_s_ut[static_cast<size_t>(j)]);
              need_scalar_fallback =
                !(R_FINITE(log_z) && log_z > log_prob_eps);
            }
          }

          if (need_scalar_fallback) {
            for (int k = 0; k < n_lR_j; ++k) {
              const int row = start + k;
              ok_int[static_cast<size_t>(k)] =
                (isok_q[row] && (!has_RACE_col || RACE_mask[row])) ? 1 : 0;
              for (int c = 0; c < n_par; ++c) {
                rowmajor[static_cast<size_t>(k) * n_par + c] = pars_q(row, c);
              }
            }
            log_z = get_trunc_normaliser_rowmajor_cpp(
                rowmajor.data(), ok_int.data(), pdf1, cdf1, lt, ut,
                n_lR_j, n_par, gsl_ctl, &node_ctx, workspace);
          }

          if (has_trunc_trial[static_cast<size_t>(j)]) {
            node_den[static_cast<size_t>(j)] =
              log_pos[static_cast<size_t>(j)] + log_z;
          }
        }
      }
    }

    if (joint_posdrift) {
      // Untruncated trials still require the positivity normaliser.  For a
      // truncated trial the block above supplied log_z; do not add a second
      // unit-mass term for it.
      for (int j = 0; j < n_unique; ++j) {
        if (has_trunc_trial[static_cast<size_t>(j)] ||
            !den_eval_mask[static_cast<size_t>(j)]) continue;
        node_den[static_cast<size_t>(j)] = log_pos[static_cast<size_t>(j)];
      }
    }
  };

  auto eval_selected = [&](const double* z_by_trial, bool need_num,
                           bool need_den, std::vector<double>& num,
                           std::vector<double>& den) {
    if (use_fast_node_eval) {
      evaluate_node_fast(z_by_trial, need_num, num, den);
      if (count_routes)
        route_counters.fused_node_evaluations += n_fallback_trials;
    } else {
      evaluate_node(z_by_trial, need_num, need_den, num, den);
    }
  };

  // Two batched passes of Gauss-Hermite nodes.
  //
  // Scan pass: a fixed shared rule whose node masses give each trial's
  // integrand mean and SD.  The scan is also used directly for broad,
  // central trials; only difficult trials pay for the recentered pass.
  //
  // Refined pass: a smaller rule re-centered per trial on those moments.
  // A weak-drift trial whose winner-implied peak sits in a tail or is narrow
  // gets its nodes placed inside that peak.  Numerator and denominator
  // integrands peak in different places (the denominator mass sits near the
  // prior mode), so each gets its own centering; the denominator refinement
  // needs no race pass.
  double max_abs_rho = 0.0;
  for (int r = 0; r < n_trials; ++r) {
    if (!isok[r]) continue;
    const double rho = effective_rho[static_cast<size_t>(r)];
    if (R_FINITE(rho) && std::fabs(rho) <= 1.0)
      max_abs_rho = std::fmax(max_abs_rho, std::fabs(rho));
  }
  // The difficult part of the integral is the narrow factor-conditioned
  // likelihood near |rho| -> 1. Moderate correlations are smooth enough for
  // a 10-node scan/refinement pair; retain the original 12-node tier for the
  // sharp cases. Environment overrides remain available for convergence
  // tests and explicit user tuning.
  //
  // A model on the generic route has no closed-form conditional: each node is
  // its own numerically integrated kernel, whose z-integrand is measurably
  // harder than BAwL's. Measured on a 24-trial RDMSWTN drift-factor design
  // against an independent 160-node reference, the sharp tier costs 1.3e-2 in
  // total log-likelihood at 12 nodes, 2.2e-3 at 16 and 6.2e-4 at 28, so the
  // sharp tier is raised there. The smooth tier already reaches 6e-4 at 10.
  const bool generic_kernel = cshared.generic_only;
  const int scan_default = (max_abs_rho <= 0.6) ? (generic_kernel ? 12 : 10)
                                                : (generic_kernel ? 28 : 12);
  const int fine_default = scan_default;
  const GHRule& scan_rule = gh_rule(
      emc2_quad_nodes("EMC2_BAWLCORR_SCAN_N", scan_default));
  const GHRule& fine_rule = gh_rule(
      emc2_quad_nodes("EMC2_BAWLCORR_FINE_N", fine_default));
  const int n_scan = static_cast<int>(scan_rule.x.size());
  const int n_fine = static_cast<int>(fine_rule.x.size());
  const double sqrt2 = std::sqrt(2.0);
  bool any_quad_den = false;
  for (int j = 0; j < n_unique; ++j) {
    if (den_eval_mask[static_cast<size_t>(j)]) { any_quad_den = true; break; }
  }

  if (count_routes) {
    long long n_quad_den = 0;
    for (int j = 0; j < n_unique; ++j)
      if (den_by_quadrature[static_cast<size_t>(j)]) ++n_quad_den;
    route_counters.den_quadrature_trials += n_quad_den;
    if (use_fast_node_eval) route_counters.gh_no_clock_trials += n_fallback_trials;
    else route_counters.gh_generic_clock_trials += n_fallback_trials;
    route_counters.scan_node_evaluations +=
      static_cast<long long>(n_scan) * n_fallback_trials;
  }

  std::vector<double> node_num;
  std::vector<double> node_den;
  std::vector<double> z_by_trial(static_cast<size_t>(n_unique));
  std::vector<double> z_scan(static_cast<size_t>(n_scan));
  std::vector<unsigned char> refine_num(static_cast<size_t>(n_unique), 0);
  for (int j = 0; j < n_unique; ++j)
    refine_num[static_cast<size_t>(j)] = exact_done[static_cast<size_t>(j)] ? 0 : 1;
  // Column-major (node fastest) log node masses from the scan pass.
  std::vector<double> scan_num(static_cast<size_t>(n_scan) * n_unique, R_NegInf);
  std::vector<double> scan_den(any_quad_den
      ? static_cast<size_t>(n_scan) * n_unique : 0, R_NegInf);
  for (int q = 0; q < n_scan; ++q) {
    const double z = sqrt2 * scan_rule.x[static_cast<size_t>(q)];
    z_scan[static_cast<size_t>(q)] = z;
    std::fill(z_by_trial.begin(), z_by_trial.end(), z);
    eval_selected(z_by_trial.data(), true, any_quad_den, node_num, node_den);
    const double log_w = std::log(gh_standard_normal_weight(scan_rule, q));
    for (int j = 0; j < n_unique; ++j) {
      scan_num[static_cast<size_t>(j) * n_scan + q] =
        log_w + node_num[static_cast<size_t>(j)];
      if (any_quad_den)
        scan_den[static_cast<size_t>(j) * n_scan + q] =
          log_w + node_den[static_cast<size_t>(j)];
    }
  }

  std::vector<AGHCenter> center_num(static_cast<size_t>(n_unique));
  std::vector<AGHCenter> center_den(
      any_quad_den ? static_cast<size_t>(n_unique) : 0);
  for (int j = 0; j < n_unique; ++j) {
    center_num[static_cast<size_t>(j)] = agh_center_from_scan(
        &scan_num[static_cast<size_t>(j) * n_scan], z_scan.data(), n_scan);
    if (any_quad_den)
      center_den[static_cast<size_t>(j)] = agh_center_from_scan(
          &scan_den[static_cast<size_t>(j) * n_scan], z_scan.data(), n_scan);
  }

  // A loaded known winner supplies an analytic proposal for the latent
  // factor.  Its marginal winner drift has density
  // (gamma0 + gamma1 V) N(V | mu, sv), restricted to the cause interval.  The
  // shared incomplete-normal moments turn that into the mean/variance of z;
  // survivor-only and unknown-winner integrands deliberately retain the scan
  // proposal because they have no winner tilt.
  std::vector<unsigned char> analytic_center_trial(static_cast<size_t>(n_unique), 0);
  if (use_fast_node_eval) {
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      const BAwLCorrTrialLayout& L = cshared.layout[static_cast<size_t>(j)];
      if (!L.winner_loaded || L.winner_row < 0) continue;
      const int row = L.winner_row;
      const BAwLTimeGeometry g = bawl_corr_row_geometry(
          cshared.rt[row], cshared.cols.data(), row, ctx);
      BAwLCorrRegion cause_region;
      if (!bawl_corr_cause_interval(g, joint_posdrift, cause_region)) continue;
      double p = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
      bawl_corr_univ_interval_moments3(
          cshared.cols[emc2col::bawl::v][row],
          cshared.cols[emc2col::bawl::sv][row], cause_region.lo,
          cause_region.hi, p, m1, m2, m3);
      const double den = g.gamma0 * p + g.gamma1 * m1;
      const double first = g.gamma0 * m1 + g.gamma1 * m2;
      const double second = g.gamma0 * m2 + g.gamma1 * m3;
      if (!(den > 0.0) || !R_FINITE(den) || !R_FINITE(first) ||
          !R_FINITE(second)) continue;
      const double mu = cshared.cols[emc2col::bawl::v][row];
      const double sd = cshared.cols[emc2col::bawl::sv][row];
      const double rho = effective_rho[static_cast<size_t>(row)];
      const double ev = first / den;
      const double vv = std::fmax(second / den - ev * ev, 0.0);
      const double loading = std::sqrt(std::fabs(rho));
      const double ez = (rho < 0.0 ? -1.0 : 1.0) * loading * (ev - mu) / sd;
      const double vz = 1.0 - std::fabs(rho) +
        std::fabs(rho) * vv / (sd * sd);
      if (R_FINITE(ez) && R_FINITE(vz) && vz > 0.0) {
        center_num[static_cast<size_t>(j)].mu = ez;
        center_num[static_cast<size_t>(j)].sigma =
          std::fmin(std::fmax(1.3 * std::sqrt(vz), 0.3), 1.25);
        analytic_center_trial[static_cast<size_t>(j)] = 1;
        if (count_routes) {
          ++route_counters.analytic_center_eligible_trials;
          ++route_counters.analytic_center_success_trials;
        }
      }
    }
  }

  auto fine_log_gauss_weights = [](const GHRule& rule) {
    std::vector<double> lgw(rule.x.size());
    for (size_t q = 0; q < rule.x.size(); ++q) {
      lgw[q] = std::log(rule.w[q]) + rule.x[q] * rule.x[q];
    }
    return lgw;
  };
  const std::vector<double> fine_lgw = fine_log_gauss_weights(fine_rule);

  std::vector<double> lw_by_trial(static_cast<size_t>(n_unique));
  // The raw scan is already a valid quadrature estimate when its mass is
  // broad and well inside the scan rule.  Keep the refinement masks separate
  // so easy trials do not re-run the BAwL kernels in the second pass; narrow
  // or tail-peaked trials retain the recentered pass for accuracy.
  if (use_fast_node_eval && !any_quad_den) {
    for (int j = 0; j < n_unique; ++j) {
      if (exact_done[static_cast<size_t>(j)]) continue;
      if (analytic_center_trial[static_cast<size_t>(j)]) {
        refine_num[static_cast<size_t>(j)] = 1;
        continue;
      }
      const AGHCenter& c = center_num[static_cast<size_t>(j)];
      const double* g = &scan_num[static_cast<size_t>(j) * n_scan];
      double peak = R_NegInf;
      int peak_i = -1;
      for (int q = 0; q < n_scan; ++q) {
        if (g[q] > peak) { peak = g[q]; peak_i = q; }
      }
      const bool at_scan_edge = peak_i == 0 || peak_i == n_scan - 1;
      // These thresholds are deliberately conservative: the normal scan is
      // reused only when it has a comfortably wide, central mass.  The fine
      // pass remains the fallback for the high-rho/narrow cases that motivated
      // the adaptive rule in the first place.
      const bool well_resolved = R_FINITE(peak) && !at_scan_edge &&
        std::fabs(c.mu) < 0.75 && c.sigma > 0.45;
      refine_num[static_cast<size_t>(j)] = well_resolved ? 0 : 1;
      if (!well_resolved) continue;
      double scan_ll = R_NegInf;
      for (int q = 0; q < n_scan; ++q)
        scan_ll = log_sum_exp(scan_ll, g[q]);
      log_num[static_cast<size_t>(j)] = scan_ll;
    }
    fine_winner = fast_winner;
    fine_loser = fast_loser;
    for (int j = 0; j < n_unique; ++j) {
      if (refine_num[static_cast<size_t>(j)]) continue;
      const int start = j * n_lR;
      const int n_lR_j = has_RACE_col ? RACE[start] : n_lR;
      for (int k = 0; k < n_lR_j; ++k) {
        fine_winner[static_cast<size_t>(start + k)] = 0;
        fine_loser[static_cast<size_t>(start + k)] = 0;
      }
    }
    use_fine_masks = true;
  }

  if (count_routes) {
    long long n_refined = 0;
    for (int j = 0; j < n_unique; ++j)
      if (refine_num[static_cast<size_t>(j)]) ++n_refined;
    route_counters.scan_refinement_trials += n_refined;
    route_counters.fine_node_evaluations +=
      static_cast<long long>(n_fine) * n_refined;
  }

  // The per-node raw evaluator consults these masks only during the optional
  // refinement pass; the scan still evaluates every active racer.
  for (int q = 0; q < n_fine; ++q) {
    const double x = sqrt2 * fine_rule.x[static_cast<size_t>(q)];
    for (int j = 0; j < n_unique; ++j) {
      const AGHCenter& c = center_num[static_cast<size_t>(j)];
      const double z = c.mu + c.sigma * x;
      z_by_trial[static_cast<size_t>(j)] = z;
      lw_by_trial[static_cast<size_t>(j)] =
        agh_log_weight(fine_lgw[static_cast<size_t>(q)], c.sigma, z);
    }
    eval_selected(z_by_trial.data(), true, false, node_num, node_den);
    for (int j = 0; j < n_unique; ++j) {
      if (!refine_num[static_cast<size_t>(j)]) continue;
      log_num[static_cast<size_t>(j)] = log_sum_exp(
          log_num[static_cast<size_t>(j)],
          lw_by_trial[static_cast<size_t>(j)] + node_num[static_cast<size_t>(j)]);
    }
  }

  if (any_quad_den) {
    // Quadrature denominators (truncated trials, or more than two racers)
    // carry the steep positivity sigmoid and the truncation window, so they
    // get a denser re-centered rule.  These evaluations skip the race pass
    // entirely and are cheap relative to the numerator sweeps.
    const GHRule& den_rule = gh_rule(64);
    const int n_den = static_cast<int>(den_rule.x.size());
    const std::vector<double> den_lgw = fine_log_gauss_weights(den_rule);
    if (count_routes) {
      long long n_quad_den = 0;
      for (int j = 0; j < n_unique; ++j)
        if (den_eval_mask[static_cast<size_t>(j)]) ++n_quad_den;
      route_counters.den_quadrature_node_evaluations +=
        static_cast<long long>(n_den) * n_quad_den;
    }
    for (int q = 0; q < n_den; ++q) {
      const double x = sqrt2 * den_rule.x[static_cast<size_t>(q)];
      for (int j = 0; j < n_unique; ++j) {
        const AGHCenter& c = center_den[static_cast<size_t>(j)];
        const double z = c.mu + c.sigma * x;
        z_by_trial[static_cast<size_t>(j)] = z;
        lw_by_trial[static_cast<size_t>(j)] =
          agh_log_weight(den_lgw[static_cast<size_t>(q)], c.sigma, z);
      }
      eval_selected(z_by_trial.data(), false, true, node_num, node_den);
      for (int j = 0; j < n_unique; ++j) {
        if (!den_eval_mask[static_cast<size_t>(j)]) continue;
        log_den[static_cast<size_t>(j)] = log_sum_exp(
            log_den[static_cast<size_t>(j)],
            lw_by_trial[static_cast<size_t>(j)] + node_den[static_cast<size_t>(j)]);
      }
    }
    if (n_den_cells < n_unique) {
      for (int j = 0; j < n_unique; ++j) {
        const int representative = den_representative[static_cast<size_t>(j)];
        if (representative >= 0 && representative != j) {
          log_den[static_cast<size_t>(j)] =
            log_den[static_cast<size_t>(representative)];
        }
      }
    }
  }

  // The flat guess mass for quadrature trials.  Exact trials already carry the
  // full mixture (bawl_corr_exact_trial_loglik); quadrature trials carry only
  // its (1 - pG) factor, so the mass is added here, where `value` is a finished
  // trial likelihood rather than an integrand.  The guess term picks up
  // log(1 - pC) because the integrand it is being summed with already has it:
  //   (1-pC)[(1-pG) L + pG g]  ==  [(1-pC)(1-pG) L] + pG (1-pC) g.
  // A non-finite RT gets log_g = -Inf, leaving just the (1 - pG) factor already
  // in the integrand -- which is the correct treatment for an omission.
  const bool guess_quad =
    cshared.pg_col >= 0 && cshared.guess.active();
  auto add_guess_mass = [&](double value, int j) -> double {
    if (!guess_quad || exact_done[static_cast<size_t>(j)]) return value;
    const int start = j * n_lR;
    const double pG = cshared.cols[static_cast<size_t>(cshared.pg_col)][start];
    if (pG <= 0.0) return value;
    const double rt_j = cshared.rt[start];
    if (!R_FINITE(rt_j) || rt_j <= 0.0) return value;
    const double pC = (cshared.pc_col >= 0)
      ? cshared.cols[static_cast<size_t>(cshared.pc_col)][start] : 0.0;
    return log_sum_exp(value,
                       std::log(pG) + log1m(pC) + cshared.guess.log_g);
  };

  double total_ll = 0.0;
  if (expand.length() > 0) {
    for (int i = 0; i < n_out; ++i) {
      const int j = expand[i] - 1;
      double value = exact_done[static_cast<size_t>(j)]
        ? exact_ll[static_cast<size_t>(j)]
        : add_guess_mass(log_num[static_cast<size_t>(j)] -
                         log_den[static_cast<size_t>(j)], j);
      if (!R_FINITE(value) || value < min_ll) value = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
      total_ll += value;
    }
  } else {
    for (int j = 0; j < n_unique; ++j) {
      double value = exact_done[static_cast<size_t>(j)]
        ? exact_ll[static_cast<size_t>(j)]
        : add_guess_mass(log_num[static_cast<size_t>(j)] -
                         log_den[static_cast<size_t>(j)], j);
      if (!R_FINITE(value) || value < min_ll) value = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
      total_ll += value;
    }
  }
  return total_ll;
}

RDMSWTNCorrSharedState build_rdmswtn_corr_shared_state(
    const Rcpp::DataFrame& dadm, int n_trials, int n_lR,
    const Rcpp::LogicalVector& winner, const Rcpp::IntegerVector& expand,
    const ContextForRaceModels& ctx, ParamTable& table,
    const Rcpp::CharacterVector& keep_names) {
  RDMSWTNCorrSharedState s;
  if (n_lR <= 0 || n_trials <= 0 || n_trials % n_lR != 0) {
    Rcpp::stop("RDMSWTNcorr: invalid race dimensions.");
  }
  s.n_trials = n_trials;
  s.n_lR = n_lR;
  s.n_unique = n_trials / n_lR;
  s.n_out = expand.length() > 0 ? expand.length() : s.n_unique;
  s.n_par = keep_names.size();
  s.rho_col = ctx.rdmswtn_rho_index;
  if (s.rho_col < 0 || s.rho_col >= s.n_par) {
    Rcpp::stop("RDMSWTNcorr: rho column is outside the parameter table.");
  }
  s.rt = dadm["rt"];
  s.LT = get_col_with_default(dadm, "LT", 0.0);
  s.UT = get_col_with_default(dadm, "UT", R_PosInf);
  s.LC = get_col_with_default(dadm, "LC", R_PosInf);
  s.UC = get_col_with_default(dadm, "UC", 0.0);
  s.lR_codes = Rcpp::IntegerVector(dadm["lR"]);

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
  for (int r = 0; r < n_trials; ++r) {
    if (winner[r]) s.winner_row[static_cast<size_t>(r / n_lR)] = r;
  }
  s.cols.assign(static_cast<size_t>(std::max(s.n_par, 16)), nullptr);
  for (int c = 0; c < s.n_par; ++c) {
    const std::string nm = Rcpp::as<std::string>(keep_names[c]);
    s.cols[static_cast<size_t>(c)] =
      &table.base(0, table.base_index_for(nm));
    if (nm == "pContaminant") s.pc_col = c;
    else if (nm == "pGuess") s.pg_col = c;
    else if (nm == "sv") s.sv_col = c;
  }
  s.guess = resolve_guess_kernel(dadm);
  s.guess.pg_col = s.pg_col;
  s.layout.assign(static_cast<size_t>(s.n_unique),
                  RDMSWTNCorrTrialLayout());
  s.valid = true;
  return s;
}

static inline bool rdmswtn_corr_row_active(
    const RDMSWTNCorrSharedState& s, int row) {
  return !s.has_RACE || s.race_mask[row];
}

static inline double rdmswtn_corr_eval_selected(
    double t, int row, const double* const* cols, int n_par,
    RacePdf1Fun fn, ContextForRaceModels* ctx) {
  if (n_par > 32) return R_NaN;
  std::array<double, 32> par{};
  for (int c = 0; c < n_par; ++c) par[static_cast<size_t>(c)] = cols[c][row];
  return fn(t, par.data(), ctx);
}

static inline double rdmswtn_corr_logcdf(
    double t, int row, const double* const* cols, int n_par,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  if (!(t > 0.0)) return R_NegInf;
  const double value = rdmswtn_corr_eval_selected(
      t, row, cols, n_par, cdf1, ctx);
  if (ISNAN(value)) return R_NaN;
  if (!(value > 0.0)) return R_NegInf;
  if (value >= 1.0) return 0.0;
  return std::log(value);
}

static inline double rdmswtn_corr_logdensity(
    double t, int row, const double* const* cols, int n_par,
    RacePdf1Fun pdf1, ContextForRaceModels* ctx) {
  if (!(t > 0.0) || !R_FINITE(t)) return R_NegInf;
  const double value = rdmswtn_corr_eval_selected(
      t, row, cols, n_par, pdf1, ctx);
  if (ISNAN(value)) return R_NaN;
  return value > 0.0 ? std::log(value) : R_NegInf;
}

static inline double rdmswtn_corr_z_from_logcdf(double log_f) {
  if (log_f == R_NegInf) return R_NegInf;
  if (log_f >= 0.0) return R_PosInf;
  return R::qnorm(log_f, 0.0, 1.0, 1, 1);
}

template <typename LogIntegrand>
static double rdmswtn_corr_integrate_log(double lo, double hi,
                                         LogIntegrand log_integrand) {
  if (!(hi > lo) || lo == R_PosInf) return R_NegInf;
  const GLRule& rule = gl_get_rule(64);
  std::array<double, 64> values{};
  double max_log = R_NegInf;
  for (int i = 0; i < 64; ++i) {
    const double x = rule.x[static_cast<size_t>(i)];
    double point, log_jac;
    if (hi == R_PosInf) {
      const double u = 0.5 * (x + 1.0);
      const double one = 1.0 - u;
      point = lo + u / one;
      log_jac = -M_LN2 - 2.0 * std::log(one);
    } else {
      point = 0.5 * (hi - lo) * x + 0.5 * (hi + lo);
      log_jac = std::log(0.5 * (hi - lo));
    }
    const double value = log_jac +
      std::log(rule.w[static_cast<size_t>(i)]) + log_integrand(point);
    values[static_cast<size_t>(i)] = value;
    max_log = std::max(max_log, value);
  }
  if (!R_FINITE(max_log)) return R_NegInf;
  double sum = 0.0;
  for (double value : values) sum += std::exp(value - max_log);
  return max_log + std::log(sum);
}

static double rdmswtn_corr_log_bvn_upper(double z1, double z2, double rho) {
  if (z1 == R_NegInf)
    return z2 == R_NegInf ? 0.0 : log_normal_upper_tail(z2);
  if (z2 == R_NegInf) return log_normal_upper_tail(z1);
  if (z1 == R_PosInf || z2 == R_PosInf) return R_NegInf;
  const double p = norm_cdf_2d_hybrid(-z1, -z2, rho);
  if (R_FINITE(p) && p > 1e-280) return std::log(std::fmin(1.0, p));

  // Direct conditional-normal integration is a tail fallback for probabilities
  // smaller than the natural-scale BVN routine can represent.
  if (z2 > z1) std::swap(z1, z2);
  const double sd = std::sqrt(std::fmax(1e-16, 1.0 - rho * rho));
  return rdmswtn_corr_integrate_log(z1, R_PosInf, [&](double x) {
    return -0.5 * x * x - 0.91893853320467274178 +
      log_normal_upper_tail((z2 - rho * x) / sd);
  });
}

static inline double rdmswtn_corr_log_pair_survival(
    double t, int row1, int row2, double rho, const double* const* cols,
    int n_par, RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const double lf1 = rdmswtn_corr_logcdf(t, row1, cols, n_par, cdf1, ctx);
  const double lf2 = rdmswtn_corr_logcdf(t, row2, cols, n_par, cdf1, ctx);
  if (ISNAN(lf1) || ISNAN(lf2)) return R_NaN;
  return rdmswtn_corr_log_bvn_upper(
      rdmswtn_corr_z_from_logcdf(lf1),
      rdmswtn_corr_z_from_logcdf(lf2), rho);
}

static inline double rdmswtn_corr_log_single_survival(
    double t, int row, const double* const* cols, int n_par,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const double lf = rdmswtn_corr_logcdf(t, row, cols, n_par, cdf1, ctx);
  if (ISNAN(lf)) return R_NaN;
  if (lf == R_NegInf) return 0.0;
  if (lf >= 0.0) return R_NegInf;
  return log1m_exp(lf);
}

static double rdmswtn_corr_log_component_survival(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout, double t, double rho,
    const double* const* cols, RaceCdf1Fun cdf1,
    ContextForRaceModels* ctx) {
  double out = rdmswtn_corr_log_pair_survival(
      t, layout.pair_row[0], layout.pair_row[1], rho, cols,
      s.n_par, cdf1, ctx);
  if (ISNAN(out)) return out;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row) ||
        row == layout.pair_row[0] || row == layout.pair_row[1]) continue;
    const double ls = rdmswtn_corr_log_single_survival(
        t, row, cols, s.n_par, cdf1, ctx);
    if (ISNAN(ls)) return ls;
    out += ls;
  }
  return out;
}

static double rdmswtn_corr_log_component_cause(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout, int winner_row, double t,
    double rho, const double* const* cols, RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1, ContextForRaceModels* ctx) {
  const int pair1 = layout.pair_row[0], pair2 = layout.pair_row[1];
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  if (!rdmswtn_corr_row_active(s, winner_row)) return R_NegInf;
  double out = rdmswtn_corr_logdensity(
      t, winner_row, cols, s.n_par, pdf1, ctx);
  if (!R_FINITE(out)) return out;

  if (winner_row == pair1 || winner_row == pair2) {
    const int loser = winner_row == pair1 ? pair2 : pair1;
    const double lfw = rdmswtn_corr_logcdf(
        t, winner_row, cols, s.n_par, cdf1, ctx);
    const double lfl = rdmswtn_corr_logcdf(
        t, loser, cols, s.n_par, cdf1, ctx);
    if (ISNAN(lfw) || ISNAN(lfl)) return R_NaN;
    const double zw = rdmswtn_corr_z_from_logcdf(lfw);
    const double zl = rdmswtn_corr_z_from_logcdf(lfl);
    // Far left tail: both CDFs can underflow to exactly zero while the winner
    // density is still positive (reachable with strongly negative drifts under
    // posdrift = FALSE).  (zl - rho * zw) is then Inf - Inf; the limit of the
    // conditional survivor is 1 for every |rho| < 1, so add nothing.
    if (!(zw == R_NegInf && zl == R_NegInf)) {
      const double sd = std::sqrt(std::fmax(1e-16, 1.0 - rho * rho));
      out += log_normal_upper_tail((zl - rho * zw) / sd);
    }
  } else {
    out += rdmswtn_corr_log_pair_survival(
        t, pair1, pair2, rho, cols, s.n_par, cdf1, ctx);
  }

  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row) || row == winner_row ||
        row == pair1 || row == pair2) continue;
    out += rdmswtn_corr_log_single_survival(
        t, row, cols, s.n_par, cdf1, ctx);
  }
  return out;
}

static double rdmswtn_corr_trial_loglik(
    const RDMSWTNCorrSharedState& s, int j,
    const RDMSWTNCorrTrialLayout& layout,
    ContextForRaceModels* ctx, const Rcpp::IntegerVector& response,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, double min_ll) {
  if (!layout.params_ok) return min_ll;
  const int start = j * s.n_lR;
  const int n_lR_j = s.has_RACE ? s.race_nacc[start] : s.n_lR;
  const double rho =
    s.cols[static_cast<size_t>(s.rho_col)][layout.pair_row[0]];
  const bool known = response[start] != NA_INTEGER;
  const int known_row = layout.winner_row;
  const double rt = s.rt[start];
  const double LT = s.LT[start], UT = s.UT[start];
  const double LC = s.LC[start], UC = s.UC[start];
  if (known && known_row < 0) return min_ll;

  int time_row = -1, nogo_row = -1, n_resp = 0;
  for (int k = 0; k < n_lR_j; ++k) {
    const int row = start + k;
    if (!rdmswtn_corr_row_active(s, row)) continue;
    const int code = s.lR_codes[row];
    if (code == ctx->time_code) time_row = row;
    else if (code == ctx->nogo_code) nogo_row = row;
    else ++n_resp;
  }

  auto single_cause = [&](int row, double t) {
    return rdmswtn_corr_log_component_cause(
        s, j, layout, row, t, rho, s.cols.data(), pdf1, cdf1, ctx);
  };
  auto finite_cause = [&](int row, double t) {
    if (row < 0) {
      double total = R_NegInf;
      for (int k = 0; k < n_lR_j; ++k) {
        const int active_row = start + k;
        if (!rdmswtn_corr_row_active(s, active_row)) continue;
        const int code = s.lR_codes[active_row];
        if (code == ctx->time_code || code == ctx->nogo_code) continue;
        double cause = single_cause(active_row, t);
        if (time_row >= 0 && n_resp > 0) {
          cause = log_sum_exp(
              cause, single_cause(time_row, t) - std::log(n_resp));
        }
        total = log_sum_exp(total, cause);
      }
      return total;
    }
    double value = single_cause(row, t);
    if (time_row >= 0 && time_row != row && n_resp > 0 &&
        s.lR_codes[row] != ctx->nogo_code) {
      value = log_sum_exp(value,
                          single_cause(time_row, t) - std::log(n_resp));
    }
    return value;
  };
  auto interval_cause = [&](double lo, double hi, int row) {
    return rdmswtn_corr_integrate_log(lo, hi, [&](double t) {
      return finite_cause(row, t);
    });
  };
  auto survivor = [&](double t) {
    if (t == 0.0) return 0.0;
    return rdmswtn_corr_log_component_survival(
        s, j, layout, t, rho, s.cols.data(), cdf1, ctx);
  };
  auto survivor_difference = [&](double lo, double hi) {
    const double a = lo == 0.0 ? 0.0 : survivor(lo);
    if (hi == R_PosInf) return a;
    return log_diff_exp(a, survivor(hi));
  };

  double value = R_NegInf;
  if (R_FINITE(rt) && rt > 0.0) {
    value = finite_cause(known ? known_row : -1, rt);
  } else if (rt == R_NegInf) {
    if (!known && nogo_row >= 0) {
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (!rdmswtn_corr_row_active(s, row) || row == nogo_row) continue;
        value = log_sum_exp(value, interval_cause(LT, LC, row));
      }
    } else {
      value = known ? interval_cause(LT, LC, known_row)
                    : survivor_difference(LT, LC);
    }
  } else if (rt == R_PosInf) {
    if (nogo_row >= 0) {
      value = log_sum_exp(interval_cause(LT, UC, nogo_row), survivor(UC));
    } else {
      value = known ? interval_cause(UC, UT, known_row)
                    : survivor_difference(UC, UT);
    }
  } else if (Rcpp::NumericVector::is_na(rt)) {
    if (known) {
      value = log_sum_exp(interval_cause(LT, LC, known_row),
                          interval_cause(UC, UT, known_row));
    } else {
      value = log_sum_exp(survivor_difference(LT, LC),
                          survivor_difference(UC, UT));
    }
  }

  if (LT != 0.0 || UT != R_PosInf) {
    const double log_z = survivor_difference(LT, UT);
    if (!R_FINITE(log_z)) return min_ll;
    value -= log_z;
  }
  // Contaminant mixture; see src/contaminant_mixture.h.
  if (s.pc_col >= 0) {
    const double p_c = s.cols[static_cast<size_t>(s.pc_col)][start];
    if (!R_FINITE(p_c) || p_c < 0.0 || p_c >= 1.0) return min_ll;
  }
  if (s.pc_col >= 0 || (s.pg_col >= 0 && s.guess.active())) {
    const double pC = (s.pc_col >= 0) ? s.cols[static_cast<size_t>(s.pc_col)][start] : 0.0;
    const double pG = (s.pg_col >= 0 && s.guess.active())
                        ? s.cols[static_cast<size_t>(s.pg_col)][start] : 0.0;
    value = mix_contaminants_rt(value, pC, pG, s.guess, rt, known);
  }
  return (!R_FINITE(value) || value < min_ll) ? min_ll : value;
}

double c_log_likelihood_rdmswtn_correlated(
    RDMSWTNCorrSharedState& s, Rcpp::DataFrame dadm,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, const int n_trials,
    LogicalVector winner, Rcpp::IntegerVector expand, double min_ll,
    const Rcpp::LogicalVector isok, int n_lR,
    void* model_context_for_funcs, bool all_finite_trials,
    RaceRawFun model_dfun_raw, RaceRawFun model_pfun_raw,
    RaceLogSAtTFun logS_at_t, RaceSharedState* shared,
    NumericVector* trial_ll_out,
    const std::function<Rcpp::NumericMatrix()>& materialize) {
  auto* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  if (ctx == nullptr || !ctx->rdmswtn_correlated ||
      ctx->rdmswtn_rho_index < 0 || !s.valid ||
      s.n_trials != n_trials || s.n_lR != n_lR) {
    Rcpp::stop("RDMSWTNcorr: invalid correlated likelihood context.");
  }
  const double* rho_col = s.cols[static_cast<size_t>(s.rho_col)];
  bool any_correlated = false;
  bool any_ordinary = false;
  Rcpp::LogicalVector ordinary_isok(n_trials, false);

  for (int j = 0; j < s.n_unique; ++j) {
    RDMSWTNCorrTrialLayout& layout = s.layout[static_cast<size_t>(j)];
    layout = RDMSWTNCorrTrialLayout();
    layout.winner_row = s.winner_row[static_cast<size_t>(j)];
    const int start = j * n_lR;
    const int n_lR_j = s.has_RACE ? s.race_nacc[start] : n_lR;
    double first_rho = 0.0;
    for (int k = 0; k < n_lR_j; ++k) {
      const int row = start + k;
      if (!rdmswtn_corr_row_active(s, row)) continue;
      const double rho = rho_col[row];
      if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
        Rcpp::stop("RDMSWTNcorr requires finite natural-scale rho values in [-1, 1].");
      }
      // posdrift = FALSE gives defective marginals; the copula is still exact
      // there (uniforms above the plateau are the atom at +Inf, so rho also
      // couples the intrinsic omissions).  sv > 0 on a correlated row is
      // reserved for a future correlated-drift model, so it is rejected rather
      // than silently given copula semantics.  See
      // .check_rdmswtn_corr_io_sv() in R/model_RDM.R.
      if (!ctx->use_posdrift && std::fabs(rho) > 1e-12) {
        const double sv_row = (s.sv_col >= 0)
          ? s.cols[static_cast<size_t>(s.sv_col)][row] : R_NaN;
        if (!(sv_row <= 1e-12)) {
          Rcpp::stop("RDMSWTNcorr with posdrift = FALSE requires sv = 0 on the correlated rows.");
        }
      }
      if (!isok[row]) layout.params_ok = false;
      if (std::fabs(rho) <= 1e-12) continue;
      if (layout.n_nonzero < 2)
        layout.pair_row[layout.n_nonzero] = row;
      ++layout.n_nonzero;
      if (layout.n_nonzero == 1) first_rho = rho;
      else if (std::fabs(rho - first_rho) > 1e-12) {
        Rcpp::stop("RDMSWTNcorr requires the two participating rows to have the same signed nonzero rho; all opted-out rows must have rho = 0.");
      }
    }
    if (layout.n_nonzero > 2) {
      Rcpp::stop("RDMSWTNcorr supports one pair per trial: at most two active rows may have nonzero rho, and all other rows must have rho = 0.");
    }
    layout.ordinary = layout.n_nonzero <= 1;
    if (layout.ordinary) {
      any_ordinary = true;
      for (int k = 0; k < n_lR_j; ++k) {
        const int row = start + k;
        if (rdmswtn_corr_row_active(s, row)) ordinary_isok[row] = isok[row];
      }
    } else {
      any_correlated = true;
    }
  }

  // Exact nesting and the cheapest route for an all-zero (or otherwise
  // unpaired) particle.
  if (!any_correlated) {
    return c_log_likelihood_race(
        materialize(), dadm, pdf1, cdf1, n_trials, winner, expand, min_ll,
        isok, n_lR, model_context_for_funcs, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, trial_ll_out, true);
  }

  std::vector<double> unique_ll(static_cast<size_t>(s.n_unique), min_ll);
  if (any_ordinary) {
    Rcpp::NumericVector ordinary_ll(s.n_unique, min_ll);
    Rcpp::IntegerVector no_expand(0);
    ContextForRaceModels ordinary_ctx = *ctx;
    ordinary_ctx.rdmswtn_correlated = false;
    c_log_likelihood_race(
        materialize(), dadm, pdf1, cdf1, n_trials, winner, no_expand, min_ll,
        ordinary_isok, n_lR, &ordinary_ctx, all_finite_trials,
        model_dfun_raw, model_pfun_raw, logS_at_t, shared, &ordinary_ll, true);
    for (int j = 0; j < s.n_unique; ++j) {
      if (s.layout[static_cast<size_t>(j)].ordinary)
        unique_ll[static_cast<size_t>(j)] = ordinary_ll[j];
    }
  }

  const Rcpp::IntegerVector response = dadm["R"];
  for (int j = 0; j < s.n_unique; ++j) {
    const RDMSWTNCorrTrialLayout& layout =
      s.layout[static_cast<size_t>(j)];
    if (!layout.ordinary) {
      unique_ll[static_cast<size_t>(j)] = rdmswtn_corr_trial_loglik(
          s, j, layout, ctx, response, pdf1, cdf1, min_ll);
    }
  }

  double total = 0.0;
  if (expand.length() > 0) {
    for (int i = 0; i < s.n_out; ++i) {
      const double value = unique_ll[static_cast<size_t>(expand[i] - 1)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = value;
      total += value;
    }
  } else {
    for (int j = 0; j < s.n_unique; ++j) {
      const double value = unique_ll[static_cast<size_t>(j)];
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = value;
      total += value;
    }
  }
  return total;
}
