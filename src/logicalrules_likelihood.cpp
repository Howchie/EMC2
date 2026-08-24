#include "utility_functions.h"
#include "logicalrules_likelihood.h"

#include "utils.h"
#include "model_LBA.h"
#include "gsl_utils.h"
#include "gh_quad.h"
#include "contaminant_mixture.h"
#include "lr_capacity_counters.h"
#include "bawl_geometry.h"
#include "bawl_corr_exact.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

LogicalRulesSharedState build_logicalrules_shared_state(const Rcpp::DataFrame& dadm,
                                                               int n_trials,
                                                               int n_acc,
                                                               bool capacity) {
  LogicalRulesSharedState out;
  out.n_trials = n_trials;
  out.n_acc = n_acc;

  if (n_trials <= 0 || n_acc <= 0 || (n_trials % n_acc) != 0) {
    Rcpp::stop("LogicalRules likelihood expects n_trials > 0 and n_trials %% n_acc == 0.");
  }
  if (!dadm.containsElementNamed("rt") ||
      !dadm.containsElementNamed("lR") ||
      !dadm.containsElementNamed("R") ||
      !dadm.containsElementNamed("LogicalRule")) {
    Rcpp::stop("LogicalRules likelihood requires dadm columns: rt, lR, R, LogicalRule.");
  }

  Rcpp::NumericVector rts = dadm["rt"];
  SEXP role_sexp = dadm["lR"];
  SEXP resp_sexp = dadm["R"];
  SEXP rule_sexp = dadm["LogicalRule"];
  Rcpp::NumericVector LT = get_col_with_default(dadm, "LT", 0.0);
  Rcpp::NumericVector UT = get_col_with_default(dadm, "UT", R_PosInf);
  Rcpp::NumericVector LC = get_col_with_default(dadm, "LC", 0.0);
  Rcpp::NumericVector UC = get_col_with_default(dadm, "UC", R_PosInf);
  const bool has_stim_col = dadm.containsElementNamed("S") ||
                            dadm.containsElementNamed("stimulus") ||
                            dadm.containsElementNamed("condition");
  SEXP stim_sexp = R_NilValue;
  bool stim_is_factor = false;
  Rcpp::IntegerVector stim_code;
  Rcpp::CharacterVector stim_chr;
  int stimNN = -1, stimA = -1, stimB = -1, stimAB = -1;
  if (has_stim_col) {
    if (dadm.containsElementNamed("S")) stim_sexp = dadm["S"];
    else if (dadm.containsElementNamed("stimulus")) stim_sexp = dadm["stimulus"];
    else stim_sexp = dadm["condition"];
    stim_is_factor = Rf_inherits(stim_sexp, "factor");
  }
  if (rts.size() != n_trials ||
      Rf_length(role_sexp) != n_trials ||
      Rf_length(resp_sexp) != n_trials ||
      Rf_length(rule_sexp) != n_trials) {
    Rcpp::stop("LogicalRules likelihood: dadm column lengths must match n_trials.");
  }

  const bool role_is_factor = Rf_inherits(role_sexp, "factor");
  const bool resp_is_factor = Rf_inherits(resp_sexp, "factor");
  const bool rule_is_factor = Rf_inherits(rule_sexp, "factor");

  Rcpp::IntegerVector role_code;
  Rcpp::IntegerVector resp_code;
  Rcpp::IntegerVector rule_code;
  Rcpp::CharacterVector role_chr;
  Rcpp::CharacterVector resp_chr;
  Rcpp::CharacterVector rule_chr;

  auto level_code = [](const Rcpp::CharacterVector& levels, const std::string& target) -> int {
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == target) return i + 1;
    }
    return -1;
  };

  int codeA = -1, codeB = -1, codenA = -1, codenB = -1, codeNogo = -1;
  int codeYes = -1, codeNo = -1;
  int codeNN = -1, codeAN = -1, codeNB = -1, codeAB = -1;
  int codeAND = -1, codeOR = -1, codeXOR = -1, codeID = -1, codeORDet = -1, codeORGng = -1;

  if (role_is_factor) {
    role_code = Rcpp::IntegerVector(role_sexp);
    Rcpp::CharacterVector lev = role_code.attr("levels");
    codeA = level_code(lev, "A");
    codeB = level_code(lev, "B");
    codenA = level_code(lev, "n_A");
    codenB = level_code(lev, "n_B");
    codeNogo = level_code(lev, "nogo");
    if (codeA < 0 || codeB < 0) {
      Rcpp::stop("LogicalRules likelihood: lR levels must include A and B.");
    }
  } else {
    role_chr = Rcpp::CharacterVector(role_sexp);
  }

  if (resp_is_factor) {
    resp_code = Rcpp::IntegerVector(resp_sexp);
    Rcpp::CharacterVector lev = resp_code.attr("levels");
    codeYes = level_code(lev, "yes");
    codeNo = level_code(lev, "no");
    codeNN = level_code(lev, "NN");
    codeAN = level_code(lev, "AN");
    codeNB = level_code(lev, "NB");
    codeAB = level_code(lev, "AB");
    if (codeYes < 0 && codeNo < 0 && codeNN < 0 && codeAN < 0 && codeNB < 0 && codeAB < 0) {
      Rcpp::stop("LogicalRules likelihood: response levels must include binary (yes/no) and/or ID (NN/AN/NB/AB) codes.");
    }
  } else {
    resp_chr = Rcpp::CharacterVector(resp_sexp);
  }

  if (rule_is_factor) {
    rule_code = Rcpp::IntegerVector(rule_sexp);
    Rcpp::CharacterVector lev = rule_code.attr("levels");
    codeAND = level_code(lev, "AND");
    codeOR = level_code(lev, "OR");
    codeXOR = level_code(lev, "XOR");
    codeID = level_code(lev, "ID");
    codeORDet = level_code(lev, "OR_DETECTION_ANALYTIC");
    codeORGng = level_code(lev, "OR_DETECTION_GNG");
    if (codeAND < 0 && codeOR < 0 && codeXOR < 0 && codeID < 0 &&
        codeORDet < 0 && codeORGng < 0) {
      Rcpp::stop("LogicalRules likelihood: LogicalRule levels must include AND, OR, XOR, ID, OR_DETECTION_ANALYTIC, or OR_DETECTION_GNG.");
    }
  } else {
    rule_chr = Rcpp::CharacterVector(rule_sexp);
  }

  if (has_stim_col && stim_is_factor) {
    stim_code = Rcpp::IntegerVector(stim_sexp);
    Rcpp::CharacterVector lev = stim_code.attr("levels");
    stimNN = level_code(lev, "NN");
    stimA = level_code(lev, "AN");
    if (stimA < 0) stimA = level_code(lev, "A");
    stimB = level_code(lev, "NB");
    if (stimB < 0) stimB = level_code(lev, "B");
    stimAB = level_code(lev, "AB");
  } else if (has_stim_col) {
    stim_chr = Rcpp::CharacterVector(stim_sexp);
  }

  out.n_unique_trials = n_trials / n_acc;
  out.idxA.assign(out.n_unique_trials, -1);
  out.idxB.assign(out.n_unique_trials, -1);
  out.idxnA.assign(out.n_unique_trials, -1);
  out.idxnB.assign(out.n_unique_trials, -1);
  out.idxNogo.assign(out.n_unique_trials, -1);
  out.rule_code.assign(out.n_unique_trials, 0);
  out.cond_code.assign(out.n_unique_trials, 0);
  out.resp_code.assign(out.n_unique_trials, 0);
  out.rt_unique.assign(out.n_unique_trials, NA_REAL);
  out.LT_unique.assign(out.n_unique_trials, 0.0);
  out.UT_unique.assign(out.n_unique_trials, R_PosInf);
  out.LC_unique.assign(out.n_unique_trials, 0.0);
  out.UC_unique.assign(out.n_unique_trials, R_PosInf);
  out.rt_by_row.assign(rts.begin(), rts.end());

  for (int j = 0; j < out.n_unique_trials; ++j) {
    const int start = j * n_acc;
    out.rt_unique[static_cast<size_t>(j)] = rts[start];
    out.LT_unique[static_cast<size_t>(j)] = LT[start];
    out.UT_unique[static_cast<size_t>(j)] = UT[start];
    out.LC_unique[static_cast<size_t>(j)] = LC[start];
    out.UC_unique[static_cast<size_t>(j)] = UC[start];

    int rcode = 0;
    if (rule_is_factor) {
      const int rc = rule_code[start];
      if (rc == codeOR) rcode = 1;
      else if (rc == codeAND) rcode = 2;
      else if (rc == codeXOR) rcode = 3;
      else if (rc == codeID) rcode = 4;
      else if (rc == codeORDet) rcode = 5;
      else if (rc == codeORGng) rcode = 6;
    } else {
      const std::string rv = Rcpp::as<std::string>(rule_chr[start]);
      if (rv == "OR") rcode = 1;
      else if (rv == "AND") rcode = 2;
      else if (rv == "XOR") rcode = 3;
      else if (rv == "ID") rcode = 4;
      else if (rv == "OR_DETECTION_ANALYTIC") rcode = 5;
      else if (rv == "OR_DETECTION_GNG") rcode = 6;
    }
    if (rcode == 0) {
      Rcpp::stop("LogicalRules likelihood: each unique trial must have a supported rule.");
    }
    out.rule_code[static_cast<size_t>(j)] = rcode;

    int resp = 0;
    if (resp_is_factor) {
      const int rc = resp_code[start];
      if (rc == NA_INTEGER) resp = 0;
      else if (rc == codeYes) resp = 1;
      else if (rc == codeNo) resp = 2;
      else if (rc == codeNN) resp = 3;
      else if (rc == codeAN) resp = 4;
      else if (rc == codeNB) resp = 5;
      else if (rc == codeAB) resp = 6;
    } else {
      if (resp_chr[start] == NA_STRING) resp = 0;
      else {
        const std::string rv = Rcpp::as<std::string>(resp_chr[start]);
        if (rv == "yes") resp = 1;
        else if (rv == "no" || rv == "nogo") resp = 2;
        else if (rv == "NN") resp = 3;
        else if (rv == "AN") resp = 4;
        else if (rv == "NB") resp = 5;
        else if (rv == "AB") resp = 6;
      }
    }
    if (resp == 0 && rcode <= 4) {
      Rcpp::stop("LogicalRules likelihood: response must be yes/no/NN/AN/NB/AB on each unique trial.");
    }
    out.resp_code[static_cast<size_t>(j)] = resp;

    // Stimulus condition codes are mandatory for the detection rules and for
    // the capacity model (which needs to identify the redundant-target AB
    // condition for every rule).  Ordinary choice designs may carry an
    // unrelated S factor; their unmapped levels are recorded as -1 and never
    // consulted.
    // OR_DETECTION_GNG (rcode 6) is the full four-horse OR task with a withheld
    // (rt = +Inf) outcome in place of the overt "no"; it does not consult the
    // stimulus condition (the design drives the drifts, exactly as for OR), so
    // only the analytic detector (rcode 5) and the capacity model require it.
    int cond = (rcode == 5 || capacity) ? -1 : 0;
    if (rcode == 5 || capacity) {
      if (!has_stim_col) {
        Rcpp::stop(capacity
          ? "LogicalRules capacity requires stimulus column `S` (or `stimulus`/`condition`)."
          : "LogicalRules detection rules require stimulus column `S` (or `stimulus`/`condition`).");
      }
    }
    if (has_stim_col) {
      if (stim_is_factor) {
        const int sc = stim_code[start];
        if (sc == NA_INTEGER || sc == stimNN) cond = 0;
        else if (sc == stimA) cond = 1;
        else if (sc == stimB) cond = 2;
        else if (sc == stimAB) cond = 3;
        else cond = -1;
      } else {
        if (stim_chr[start] == NA_STRING) cond = 0;
        else {
          const std::string sv = Rcpp::as<std::string>(stim_chr[start]);
          if (sv == "NN" || sv == "none") cond = 0;
          else if (sv == "AN" || sv == "A") cond = 1;
          else if (sv == "NB" || sv == "B") cond = 2;
          else if (sv == "AB" || sv == "BA" || sv == "A+B" || sv == "B+A") cond = 3;
          else cond = -1;
        }
      }
    }
    if ((rcode == 5 || capacity) && (cond < 0 || cond > 3)) {
      Rcpp::stop(capacity
        ? "LogicalRules capacity requires stimulus NN/AN/NB/AB (or A/B/AB) on every trial."
        : "LogicalRules detection rules require stimulus NN/AN/NB/AB (or A/B/AB).");
    }
    out.cond_code[static_cast<size_t>(j)] = cond;

    for (int k = 0; k < n_acc; ++k) {
      const int idx = start + k;
      if (role_is_factor) {
        const int rc = role_code[idx];
        if (rc == codeA) out.idxA[static_cast<size_t>(j)] = idx;
        else if (rc == codeB) out.idxB[static_cast<size_t>(j)] = idx;
        else if (rc == codenA) out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rc == codenB) out.idxnB[static_cast<size_t>(j)] = idx;
        else if (rc == codeNogo) out.idxNogo[static_cast<size_t>(j)] = idx;
      } else {
        const std::string rv = Rcpp::as<std::string>(role_chr[idx]);
        if (rv == "A") out.idxA[static_cast<size_t>(j)] = idx;
        else if (rv == "B") out.idxB[static_cast<size_t>(j)] = idx;
        else if (rv == "n_A") out.idxnA[static_cast<size_t>(j)] = idx;
        else if (rv == "n_B") out.idxnB[static_cast<size_t>(j)] = idx;
        else if (rv == "nogo") out.idxNogo[static_cast<size_t>(j)] = idx;
      }
    }
    if (out.idxNogo[static_cast<size_t>(j)] >= 0) out.has_nogo = true;

    if (rcode <= 4 || rcode == 6) {
      // OR/AND/XOR/ID and OR_DETECTION_GNG all use the four-horse subrace layout
      // (two target/nontarget channels).  GNG differs only in that the "no"
      // outcome is a withheld response (rt = +Inf) rather than an overt "no".
      if (out.idxA[static_cast<size_t>(j)] < 0 ||
          out.idxB[static_cast<size_t>(j)] < 0 ||
          out.idxnA[static_cast<size_t>(j)] < 0 ||
          out.idxnB[static_cast<size_t>(j)] < 0) {
        Rcpp::stop(rcode == 6
          ? "LogicalRules likelihood: OR_DETECTION_GNG trials require A, n_A, B, n_B accumulators."
          : "LogicalRules likelihood: OR/AND/XOR/ID trials require A, B, n_A, n_B accumulators.");
      }
    } else if (rcode == 5) {
      if (out.idxA[static_cast<size_t>(j)] < 0 || out.idxB[static_cast<size_t>(j)] < 0) {
        Rcpp::stop("LogicalRules likelihood: OR_DETECTION_ANALYTIC trials require A and B accumulators.");
      }
    }
  }

  const int n_unique_trials = out.n_unique_trials;
  out.cell_id.resize(n_unique_trials);
  std::vector<Rcpp::IntegerVector> dm_expands;
  if (dadm.hasAttribute("designs")) {
    Rcpp::List designs = dadm.attr("designs");
    for (int d = 0; d < designs.size(); ++d) {
      if (TYPEOF(designs[d]) == EXTPTRSXP) continue;
      Rcpp::RObject obj = designs[d];
      if (obj.hasAttribute("expand")) {
        dm_expands.push_back(obj.attr("expand"));
      }
    }
  }
  std::vector<int> acc_des_id(static_cast<size_t>(n_trials * n_acc), 0);
  std::vector<std::vector<int>> unique_profiles;
  for (int r = 0; r < n_trials * n_acc; ++r) {
    std::vector<int> profile;
    if (!dm_expands.empty()) {
      for (const auto& exp : dm_expands) {
        profile.push_back((r >= 0 && r < exp.size()) ? exp[r] : r);
      }
    } else {
      profile.push_back(r);
    }
    int found_p = -1;
    for (size_t u = 0; u < unique_profiles.size(); ++u) {
      if (unique_profiles[u] == profile) {
        found_p = static_cast<int>(u);
        break;
      }
    }
    if (found_p < 0) {
      found_p = static_cast<int>(unique_profiles.size());
      unique_profiles.push_back(profile);
    }
    acc_des_id[static_cast<size_t>(r)] = found_p;
  }
  auto get_des_idx = [&acc_des_id](int idx) -> int {
    if (idx >= 0 && idx < static_cast<int>(acc_des_id.size())) {
      return acc_des_id[static_cast<size_t>(idx)];
    }
    return idx;
  };
  for (int j = 0; j < n_unique_trials; ++j) {
    LogicalRulesCellKey k{
      get_des_idx(out.idxA[j]), get_des_idx(out.idxB[j]),
      get_des_idx(out.idxnA[j]), get_des_idx(out.idxnB[j]),
      get_des_idx(out.idxNogo[j]),
      out.rule_code[j], out.cond_code[j],
      out.LT_unique[j], out.UT_unique[j]
    };
    int found = -1;
    for (size_t c = 0; c < out.cells.size(); ++c) {
      if (out.cells[c] == k) {
        found = static_cast<int>(c);
        break;
      }
    }
    if (found < 0) {
      found = static_cast<int>(out.cells.size());
      out.cells.push_back(k);
    }
    out.cell_id[j] = found;
  }
  out.n_cells = static_cast<int>(out.cells.size());
  out.cell_log_z.assign(static_cast<size_t>(out.n_cells), NA_REAL);
  out.cell_log_den_cap.assign(static_cast<size_t>(out.n_cells), NA_REAL);

  out.valid = true;
  return out;
}

static inline void copy_par_row_colmajor(const double* const* cols,
                                         int n_par,
                                         int row_idx,
                                         double* out_row) {
  for (int c = 0; c < n_par; ++c) {
    out_row[c] = cols[c][row_idx];
  }
}

static inline bool row_equal_colmajor(const double* const* cols,
                                      int n_par,
                                      int row_a,
                                      int row_b) {
  for (int c = 0; c < n_par; ++c) {
    if (cols[c][row_a] != cols[c][row_b]) {
      return false;
    }
  }
  return true;
}

struct LogicalRulesScratch {
  std::vector<double> ll_unique;
  std::vector<double> logf_all;
  std::vector<double> logS_all;
  std::vector<int> all_mask;
  std::vector<int> isok_int_all;
  std::vector<double> GA_no_gl;
  std::vector<double> GB_no_gl;
  std::vector<unsigned char> ch_eq_vec;
  std::vector<unsigned char> pair_ok_A;
  std::vector<unsigned char> pair_ok_B;
  std::vector<double> gl_h_A;
  std::vector<double> gl_h_B;
  std::vector<double> gl_lo_A;
  std::vector<double> gl_lo_B;
  std::vector<double> pars_nA;
  std::vector<double> pars_A;
  std::vector<double> pars_nB;
  std::vector<double> pars_B;
  std::vector<int> isok_nA;
  std::vector<int> isok_A;
  std::vector<int> isok_nB;
  std::vector<int> isok_B;
  // Auxiliary GL sweeps (G_no at censor/truncation time points), slot-major
  // flat layout [slot * n_unique_trials + j]; see LrAuxSlot.
  std::vector<double> aux_t;
  std::vector<double> aux_hA;
  std::vector<double> aux_hB;
  std::vector<double> aux_GA;
  std::vector<double> aux_GB;
  std::vector<int> aux_mask;
  std::vector<int> isok_nA_base;
  std::vector<int> isok_A_base;
  std::vector<int> isok_nB_base;
  std::vector<int> isok_B_base;
  std::vector<int> all1;
  std::vector<double> tmpA1;
  std::vector<double> tmpA2;
  std::vector<double> tmpB1;
  std::vector<double> tmpB2;
  std::vector<double> parA;
  std::vector<double> parB;
  std::vector<double> parnA;
  std::vector<double> parnB;
  std::vector<double> pars_2buf;
};

// Auxiliary time-point slots for the logical-rules GL pre-pass: G_no is
// batch-evaluated at these per-trial times so the trial loop's truncation
// normaliser and censored branches need no per-trial GSL integration on the
// raw-kernel path.
enum LrAuxSlot {
  LR_AUX_LT = 0,   // truncation lower bound (Z), OR/AND only
  LR_AUX_UT,       // truncation upper bound (Z + upper-censor subtraction)
  LR_AUX_CLO,      // lower-censor window start max(0, LT)
  LR_AUX_CHI,      // lower-censor window end LC
  LR_AUX_UC,       // upper-censor bound, OR/AND only
  LR_AUX_N
};

static inline void resize_assign_double(std::vector<double>& x, size_t n, double value) {
  x.assign(n, value);
}

static inline void resize_assign_int(std::vector<int>& x, size_t n, int value) {
  x.assign(n, value);
}

static double local_race_helper(double t,
                                const double* par_target,
                                const double* par_nontarget,
                                int n_par,
                                ContextForRaceModels* ctx,
                                RacePdf1Fun pdf1,
                                RaceCdf1Fun cdf1,
                                const GslIntegrationControls& gsl_ctl,
                                gsl_integration_workspace* w,
                                double* pars_2buf,
                                int* isok_2buf) {
  if (!(t > 0.0) || !R_FINITE(t)) return 0.0;
  if (ctx == nullptr || pdf1 == nullptr || cdf1 == nullptr) return NA_REAL;

  const int t0_idx = ctx->t0_index;
  double t0_tgt = (t0_idx >= 0 && t0_idx < n_par) ? par_target[t0_idx] : 0.0;
  double t0_nt  = (t0_idx >= 0 && t0_idx < n_par) ? par_nontarget[t0_idx] : 0.0;
  if (t0_tgt < 0.0) t0_tgt = 0.0;
  if (t0_nt < 0.0) t0_nt = 0.0;

  if (t <= t0_tgt) return 0.0;

  double term1 = 0.0;
  const double split = std::min(t, t0_nt);
  if (split > t0_tgt) {
    double c_hi = cdf1(split, par_target, ctx);
    double c_lo = cdf1(t0_tgt, par_target, ctx);
    c_hi = clamp_cdf01_race(c_hi);
    c_lo = clamp_cdf01_race(c_lo);
    if (emc2_isnan(c_hi) || emc2_isnan(c_lo)) return NA_REAL;
    term1 = c_hi - c_lo;
    if (term1 < 0.0) term1 = 0.0;
    if (term1 > 1.0) term1 = 1.0;
  }

  const double low_limit = std::max(t0_tgt, t0_nt);
  if (t <= low_limit) return term1;

  const GLRule& gl = gl_get_rule(31);
  const double h = 0.5 * (t - low_limit);
  double sum_g = 0.0;
  for (size_t i = 0; i < gl.x.size(); ++i) {
    const double s = low_limit + h * (1.0 + gl.x[i]);
    const double f = pdf1(s, par_target, ctx);
    if (!R_FINITE(f) || f <= 0.0) continue;
    double F = clamp_cdf01_race(cdf1(s, par_nontarget, ctx));
    if (emc2_isnan(F)) return NA_REAL;
    const double S = std::max(0.0, 1.0 - F);
    sum_g += gl.w[i] * h * f * S;
  }

  double g = term1 + sum_g;
  if (g < 0.0) g = 0.0;
  if (g > 1.0) g = 1.0;
  return g;
}

// P(target strictly wins its two-accumulator channel by time t), i.e.
// int_0^t f_target(u) S_nontarget(u) du.  Unlike local_race_helper this
// accepts t = +Inf, integrating the tail with Gauss-Legendre quadrature.
// Returns a probability in [0, 1], or NA_REAL on integration failure.
static double lr_channel_yes_prob(double t,
                                  const double* par_target,
                                  const double* par_nontarget,
                                  int n_par,
                                  ContextForRaceModels* ctx,
                                  RacePdf1Fun pdf1,
                                  RaceCdf1Fun cdf1,
                                  const GslIntegrationControls& gsl_ctl,
                                  gsl_integration_workspace* w,
                                  double* pars_2buf,
                                  int* isok_2buf) {
  if (!(t > 0.0)) return 0.0;
  if (R_FINITE(t))
    return local_race_helper(t, par_target, par_nontarget, n_par, ctx,
                             pdf1, cdf1, gsl_ctl, w, pars_2buf, isok_2buf);
  if (ctx == nullptr || pdf1 == nullptr || cdf1 == nullptr)
    return NA_REAL;

  const int t0_idx = ctx->t0_index;
  double t0_tgt = (t0_idx >= 0 && t0_idx < n_par) ? par_target[t0_idx] : 0.0;
  double t0_nt  = (t0_idx >= 0 && t0_idx < n_par) ? par_nontarget[t0_idx] : 0.0;
  if (t0_tgt < 0.0) t0_tgt = 0.0;
  if (t0_nt < 0.0) t0_nt = 0.0;

  // Early window [t0_tgt, t0_nt] where only the target can finish (nontarget
  // survivor is exactly one): a plain CDF gap of the target.
  double term1 = 0.0;
  if (t0_nt > t0_tgt) {
    double c_hi = clamp_cdf01_race(cdf1(t0_nt, par_target, ctx));
    double c_lo = clamp_cdf01_race(cdf1(t0_tgt, par_target, ctx));
    if (emc2_isnan(c_hi) || emc2_isnan(c_lo)) return NA_REAL;
    term1 = c_hi - c_lo;
    if (term1 < 0.0) term1 = 0.0;
    if (term1 > 1.0) term1 = 1.0;
  }

  const double low_limit = std::max(t0_tgt, t0_nt);
  const GLRule& gl = gl_get_rule(31);
  double tail = 0.0;
  for (size_t i = 0; i < gl.x.size(); ++i) {
    const double u = 0.5 * (1.0 + gl.x[i]);
    const double om_u = 1.0 - u;
    if (om_u <= 1e-12) continue;
    const double s = low_limit + u / om_u;
    const double jacobian = 0.5 * gl.w[i] / (om_u * om_u);
    const double f = pdf1(s, par_target, ctx);
    if (!R_FINITE(f) || f <= 0.0) continue;
    double F = clamp_cdf01_race(cdf1(s, par_nontarget, ctx));
    if (emc2_isnan(F)) return NA_REAL;
    const double S = std::max(0.0, 1.0 - F);
    tail += jacobian * f * S;
  }

  double g = term1 + tail;
  if (g < 0.0) g = 0.0;
  if (g > 1.0) g = 1.0;
  return g;
}

static inline bool eval_pdf_cdf_race_scalar(
    bool use_raw_local,
    int idx,
    double t,
    const double* p_row,
    double min_ll,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    double& f_out,
    double& F_out) {
  if (use_raw_local) {
    const double lf = logf_all[static_cast<size_t>(idx)];
    const double ls = logS_all[static_cast<size_t>(idx)];
    if (!R_FINITE(lf) || !R_FINITE(ls)) return false;
    const double f = (lf <= min_ll) ? 0.0 : std::exp(lf);
    double S = (ls <= min_ll) ? 0.0 : std::exp(ls);
    if (!(S >= 0.0) || !R_FINITE(S)) return false;
    if (S > 1.0) S = 1.0;
    const double F = clamp_cdf01_race(1.0 - S);
    if (!R_FINITE(f) || !R_FINITE(F) || emc2_isnan(F)) return false;
    f_out = f;
    F_out = F;
    return true;
  }
  const double f = pdf1(t, p_row, model_ctx);
  if (!R_FINITE(f) || f < 0.0) return false;
  double F = cdf1(t, p_row, model_ctx);
  F = clamp_cdf01_race(F);
  if (!R_FINITE(F) || emc2_isnan(F)) return false;
  f_out = f;
  F_out = F;
  return true;
}

static inline double safe_surv_at_race_scalar(
    const double* p_row,
    double t_eval,
    double min_surv,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx) {
  double F = cdf1(t_eval, p_row, model_ctx);
  F = clamp_cdf01_race(F);
  if (!R_FINITE(F) || emc2_isnan(F)) return 0.0;
  return std::max(min_surv, 1.0 - F);
}

// ---------------------------------------------------------------------------
// Logical-rules channel-resolution probabilities at a fixed time t.
//
// A "channel" is the 2-accumulator sub-race between a target ("yes") detector
// and a nontarget ("no") detector. At time t the channel is in exactly one of
// three states, with G_yes + G_no + S_dec = 1:
//   G_yes(t) = P(target won the sub-race at some u <= t)     channel said yes
//   G_no(t)  = P(nontarget won at some u <= t)               channel said no
//   S_dec(t) = P(neither finished by t) = S_tgt(t)*S_ntgt(t) undecided
// ---------------------------------------------------------------------------
struct LrChannelState {
  double G_yes = 0.0;
  double G_no = 0.0;
  double S_dec = 1.0;   // t = 0 state: undecided with certainty
};

// P(no overt response by t) for the OR/AND/XOR/ID rules. An overt response is
// NOT the first accumulator event: OR responds at the first channel-yes or the
// second channel-no; AND at the second yes or the first no; XOR/ID only once
// BOTH channels have resolved. Derived from channel independence, and
// consistent with the finite-rt densities used in the main trial loop
// (e.g. d/dt[1 - Q_A*Q_B] = gA_yes*Q_B + gB_yes*Q_A with Q = G_no + S_dec).
static inline double lr_no_response_prob(int rule_code,
                                         const LrChannelState& A,
                                         const LrChannelState& B) {
  if (rule_code == 3 || rule_code == 4) {  // XOR/ID: waiting iff either channel undecided
    return 1.0 - (1.0 - A.S_dec) * (1.0 - B.S_dec);
  }
  if (rule_code == 1) {                    // OR: no yes yet, minus both-no (responded "no")
    const double Q_A = A.G_no + A.S_dec;
    const double Q_B = B.G_no + B.S_dec;
    return Q_A * Q_B - A.G_no * B.G_no;
  }
  // AND: no no yet, minus both-yes (responded "yes")
  const double Q_A = A.G_yes + A.S_dec;
  const double Q_B = B.G_yes + B.S_dec;
  return Q_A * Q_B - A.G_yes * B.G_yes;
}

// CDF of (response identity, RT <= t) for OR/AND/XOR/ID. resp_code follows
// LogicalRulesSharedState (0=unknown, 1=yes, 2=no, 3=NN, 4=AN, 5=NB, 6=AB);
// resp_code 0 marginalizes over identities. Returns NA_REAL for a response
// code that is impossible under the rule.
static inline double lr_response_cdf(int rule_code, int resp_code,
                                     const LrChannelState& A,
                                     const LrChannelState& B) {
  if (resp_code == 0) return 1.0 - lr_no_response_prob(rule_code, A, B);
  switch (rule_code) {
  case 1:  // OR: yes at first channel-yes; no once both channels said no
    if (resp_code == 1) return 1.0 - (A.G_no + A.S_dec) * (B.G_no + B.S_dec);
    if (resp_code == 2) return A.G_no * B.G_no;
    break;
  case 2:  // AND: yes once both channels said yes; no at first channel-no
    if (resp_code == 1) return A.G_yes * B.G_yes;
    if (resp_code == 2) return 1.0 - (A.G_yes + A.S_dec) * (B.G_yes + B.S_dec);
    break;
  case 3:  // XOR: yes = exactly one channel yes; both channels resolved
    if (resp_code == 1) return A.G_yes * B.G_no + A.G_no * B.G_yes;
    if (resp_code == 2) return A.G_yes * B.G_yes + A.G_no * B.G_no;
    break;
  case 4:  // ID: identity = the resolved (yes/no) pair
    if (resp_code == 6) return A.G_yes * B.G_yes;
    if (resp_code == 4) return A.G_yes * B.G_no;
    if (resp_code == 5) return A.G_no  * B.G_yes;
    if (resp_code == 3) return A.G_no  * B.G_no;
    break;
  }
  return NA_REAL;
}

// P(no overt response by t) for the legacy analytic-detection helper. The
// four-horse GNG route is handled by lr4_N_at below; it has no separate nogo
// accumulator. Used for upper-censor masses and the truncation normaliser.
// t must be finite; t <= 0 returns 1.
static double lr_detection_no_response_prob(
    int rule_code,
    int cond,
    double t,
    const double* parA,
    const double* parB,
    const double* parN,
    int n_par,
    double min_surv,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& pars_rowmajor_buf,
    bool& ok) {
  ok = true;
  if (!(t > 0.0)) return 1.0;
  if (!R_FINITE(t)) { ok = false; return 1.0; }
  if (cond == 0 && rule_code == 5) return 1.0;  // no detector: no response possible
  const double S_A = (cond == 1 || cond == 3)
    ? safe_surv_at_race_scalar(parA, t, min_surv, cdf1, model_ctx) : 1.0;
  const double S_B = (cond == 2 || cond == 3)
    ? safe_surv_at_race_scalar(parB, t, min_surv, cdf1, model_ctx) : 1.0;
  if (rule_code == 5) return S_A * S_B;

  // rule 6, cond 0: only the nogo accumulator races, so no overt response can
  // ever occur — S_N(t) + P(nogo won by t) = 1 identically.
  if (cond == 0) return 1.0;

  // rule 6: joint survivor of the active set, plus nogo-win mass on [0, t]
  const double S_N = safe_surv_at_race_scalar(parN, t, min_surv, cdf1, model_ctx);
  const double p_none = S_A * S_B * S_N;

  const int t0_idx = model_ctx->t0_index;
  const double t0_N = (t0_idx >= 0 && t0_idx < n_par) ? parN[t0_idx] : 0.0;
  const double lo = std::max(0.0, t0_N);
  const double h = (t > lo) ? 0.5 * (t - lo) : 0.0;
  double p_nogo_win = 0.0;
  if (h > 0.0) {
    const GLRule& gl = gl_get_rule(31);
    for (size_t i = 0; i < gl.x.size(); ++i) {
      const double s = lo + h * (1.0 + gl.x[i]);
      const double fN = pdf1(s, parN, model_ctx);
      if (!R_FINITE(fN) || fN <= 0.0) continue;
      double SA = 1.0, SB = 1.0;
      if (cond == 1 || cond == 3) SA = safe_surv_at_race_scalar(parA, s, min_surv, cdf1, model_ctx);
      if (cond == 2 || cond == 3) SB = safe_surv_at_race_scalar(parB, s, min_surv, cdf1, model_ctx);
      p_nogo_win += gl.w[i] * h * fN * SA * SB;
    }
  }

  return std::min(1.0, p_none + p_nogo_win);
}

static double logicalrules_detection_trial_ll(
    int rule_code,
    int cond,
    int resp,
    double t,
    double LTj,
    double UTj,
    double LCj,
    double UCj,
    int idxA,
    int idxB,
    int idxN,
    int n_par,
    double min_ll,
    double min_surv,
    double log_Z,
    bool use_raw_local,
    const std::vector<double>& logf_all,
    const std::vector<double>& logS_all,
    const double* parA,
    const double* parB,
    const double* parN,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    ContextForRaceModels* model_ctx,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& pars_rowmajor_buf) {
  const bool finite_rt = (R_FINITE(t) && t > 0.0);
  double p_j = 0.0;

  if (finite_rt) {
    if (resp != 1) return min_ll;
    // A detection stimulus activates only the detector(s) named by `cond`.
    // The fixed-role design still carries both A and B rows so that every
    // trial has the same layout, but an inactive row must not be evaluated:
    // besides being unnecessary, an invalid inactive parameter could
    // otherwise floor an otherwise valid A-only or B-only trial.
    const bool activeA = (cond == 1 || cond == 3);
    const bool activeB = (cond == 2 || cond == 3);
    if (!activeA && !activeB) return min_ll;
    double fA = 0.0, FA = 0.0, fB = 0.0, FB = 0.0;
    if (activeA && !eval_pdf_cdf_race_scalar(
          use_raw_local, idxA, t, parA, min_ll, pdf1, cdf1, model_ctx,
          logf_all, logS_all, fA, FA)) return min_ll;
    if (activeB && !eval_pdf_cdf_race_scalar(
          use_raw_local, idxB, t, parB, min_ll, pdf1, cdf1, model_ctx,
          logf_all, logS_all, fB, FB)) return min_ll;
    const double S_A = activeA ? std::max(min_surv, 1.0 - FA) : 1.0;
    const double S_B = activeB ? std::max(min_surv, 1.0 - FB) : 1.0;
    if (rule_code == 5) {
      if (cond == 1) p_j = fA;
      else if (cond == 2) p_j = fB;
      else p_j = fA * S_B + fB * S_A;
    } else {
      double fN = 0.0, FN = 0.0;
      if (!eval_pdf_cdf_race_scalar(
            use_raw_local, idxN, t, parN, min_ll, pdf1, cdf1, model_ctx,
            logf_all, logS_all, fN, FN)) return min_ll;
      const double S_N = std::max(min_surv, 1.0 - FN);
      if (cond == 1) p_j = fA * S_N;
      else if (cond == 2) p_j = fB * S_N;
      else p_j = (fA * S_B + fB * S_A) * S_N;
    }
  } else {
    // Upper censor (+Inf) is an omission: recorded resp may be missing (0) or
    // an explicit "no" (2). Lower censor (-Inf) is an overt detection response
    // that beat LC: resp may be missing (0) or "yes" (1) — the mass is the
    // same either way, as these rules have a single overt response type.
    const bool upper_censored = (t > 0.0);
    if (upper_censored ? (resp != 0 && resp != 2)
                       : (resp != 0 && resp != 1)) return min_ll;
    if (rule_code == 5) {
      // Upper censor (rt = +Inf): neither detector fired by UC. With a finite
      // upper truncation bound (and a proper upper tail) the response is known
      // to land in (UC, UT], so N(UT) is subtracted; the trial mass is then
      // normalised by Z = N(LT) - N(UT) in the caller.
      // Lower censor (rt = -Inf): at least one detector fired before LC.
      if (t > 0.0) {
        if (!R_FINITE(UCj))
          Rcpp::stop("LogicalRules OR_DETECTION_ANALYTIC requires finite UC for upper-censored trials.");
        bool okd = true;
        p_j = lr_detection_no_response_prob(rule_code, cond, UCj, parA, parB, parN,
                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                            gsl_ctl, w, pars_rowmajor_buf, okd);
        if (okd && R_FINITE(UTj) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT = lr_detection_no_response_prob(rule_code, cond, UTj, parA, parB, parN,
                                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                                            gsl_ctl, w, pars_rowmajor_buf, ok_ut);
          okd = okd && ok_ut;
          p_j = std::max(0.0, p_j - N_UT);
        }
        if (!okd) return min_ll;
      } else {
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, LCj);
        if (cond == 0) {
          p_j = 0.0;
        } else if (cond == 1) {
          const double S_A_lo = safe_surv_at_race_scalar(parA, lo, min_surv, cdf1, model_ctx);
          const double S_A_hi = safe_surv_at_race_scalar(parA, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_A_lo - S_A_hi);
        } else if (cond == 2) {
          const double S_B_lo = safe_surv_at_race_scalar(parB, lo, min_surv, cdf1, model_ctx);
          const double S_B_hi = safe_surv_at_race_scalar(parB, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_B_lo - S_B_hi);
        } else {
          const double S_A_lo = safe_surv_at_race_scalar(parA, lo, min_surv, cdf1, model_ctx);
          const double S_B_lo = safe_surv_at_race_scalar(parB, lo, min_surv, cdf1, model_ctx);
          const double S_A_hi = safe_surv_at_race_scalar(parA, hi, min_surv, cdf1, model_ctx);
          const double S_B_hi = safe_surv_at_race_scalar(parB, hi, min_surv, cdf1, model_ctx);
          p_j = std::max(0.0, S_A_lo * S_B_lo - S_A_hi * S_B_hi);
        }
      }
    } else {
      if (t > 0.0) {
        // Upper censor: no overt response by UC = nogo won by UC or nobody
        // fired by UC. The nogo-win mass now accumulates from 0 (not LT):
        // truncation conditioning is handled by the caller's normaliser Z,
        // and N(UT) is subtracted when the window's upper bound is finite.
        if (!R_FINITE(UCj))
          Rcpp::stop("LogicalRules OR_DETECTION_GNG requires finite UC for upper-censored trials.");
        bool okd = true;
        p_j = lr_detection_no_response_prob(rule_code, cond, UCj, parA, parB, parN,
                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                            gsl_ctl, w, pars_rowmajor_buf, okd);
        if (okd && R_FINITE(UTj) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT = lr_detection_no_response_prob(rule_code, cond, UTj, parA, parB, parN,
                                                            n_par, min_surv, pdf1, cdf1, model_ctx,
                                                            gsl_ctl, w, pars_rowmajor_buf, ok_ut);
          okd = okd && ok_ut;
          p_j = std::max(0.0, p_j - N_UT);
        }
        if (!okd) return min_ll;
      } else {
        // Lower censoring records an overt response, not a nogo finish.
        // Evaluate the response mass from the active go channels.
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, LCj);
        if (cond == 0 || !(hi > lo)) {
          p_j = 0.0;
        } else {
          const GLRule& gl = gl_get_rule(31);
          const double h = 0.5 * (hi - lo);
          double p_go_win = 0.0;
          for (size_t i = 0; i < gl.x.size(); ++i) {
            const double s = lo + h * (1.0 + gl.x[i]);
            const double SN = safe_surv_at_race_scalar(parN, s, min_surv, cdf1, model_ctx);
            if (cond == 1 || cond == 3) {
              const double fA = pdf1(s, parA, model_ctx);
              if (R_FINITE(fA) && fA > 0.0) {
                const double SB = (cond == 3) ? safe_surv_at_race_scalar(parB, s, min_surv, cdf1, model_ctx) : 1.0;
                p_go_win += gl.w[i] * h * fA * SB * SN;
              }
            }
            if (cond == 2 || cond == 3) {
              const double fB = pdf1(s, parB, model_ctx);
              if (R_FINITE(fB) && fB > 0.0) {
                const double SA = (cond == 3) ? safe_surv_at_race_scalar(parA, s, min_surv, cdf1, model_ctx) : 1.0;
                p_go_win += gl.w[i] * h * fB * SA * SN;
              }
            }
          }
          p_j = p_go_win;
        }
      }
    }
  }

  if (!(p_j > 0.0) || !R_FINITE(p_j)) return min_ll;
  const double ll = std::log(p_j) - log_Z;
  return (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
}

// ===========================================================================
// LogicalRules correlated-target capacity routes.
//
// On a redundant-target (AB) trial one latent standard-normal factor z adds
// a shared shift to both target drift means: V_i | z ~
// N(v_i + kappa + tau*z, sv_i^2) for the
// A and B target rows; nontarget/nogo racers never load on the factor.  The
// logical rule is assembled conditional on z from normalized conditional
// channel quantities and the complete trial mass is integrated over z with a
// scan/recentred Gauss-Hermite pair (gh_quad.h).  Positive-drift models
// weight each node by q_A(z)*q_B(z) and divide once by the closed-form
// bivariate orthant probability q_AB (bawl_corr_pair_positive_normalizer);
// independent racers keep their ordinary posdrift-normalized kernels, whose
// constant q factors cancel between numerator and denominator.
// LogicalRulesLBA is plain LBA, i.e. the exact k = 0 member of the BAwL
// geometry, so all conditional target endpoints come from BAwLPreparedRow.
// ===========================================================================

// GL order for the conditional channel integrals; matches the batched GL
// pre-pass order used by the ordinary route.
static constexpr int LRCAP_N_GL = 31;

// Node-invariant conditional drift view of one loaded target row:
// V | z ~ N(v0 + slope*z, sv^2) with v0 = v + kappa and slope = tau.
struct LrCapTargetView {
  double v0 = 0.0, slope = 0.0, sv = 1.0;
  double t0 = 0.0, A = 0.0, b = 0.0;
  bool valid = false;
};

static inline LrCapTargetView lrcap_target_view(const double* const* cols,
                                                int row, double kappa,
                                                double tau) {
  LrCapTargetView tv;
  const double v = cols[emc2col::lba::v][row];
  const double sv = cols[emc2col::lba::sv][row];
  const double B = cols[emc2col::lba::B][row];
  const double A = cols[emc2col::lba::A][row];
  tv.t0 = cols[emc2col::lba::t0][row];
  tv.A = A;
  tv.b = B + A;
  tv.v0 = v + kappa;
  tv.slope = tau;
  tv.sv = sv;
  tv.valid = R_FINITE(v) && sv > 0.0 && A >= 0.0 && tv.b > 0.0 && tv.b >= A &&
    R_FINITE(tv.v0) && R_FINITE(tv.slope) &&
    R_FINITE(tv.t0) && tv.t0 >= 0.0;
  return tv;
}

static inline BAwLPreparedRow lrcap_prepare_at(const LrCapTargetView& tv,
                                               double t) {
  return bawl_prepare_row(bawl_time_geometry(t, tv.t0, tv.A, tv.b, 0.0),
                          tv.v0, tv.slope, tv.sv);
}

// log q(z) = log Phi(v_q / sv): the target's positivity constant, shared by
// every evaluation time of the row (it depends only on the drift model).
static inline double lrcap_log_q(const LrCapTargetView& tv, double z,
                                 bool posdrift) {
  if (!posdrift) return 0.0;
  return pnorm_log_direct((tv.v0 + tv.slope * z) / tv.sv, true);
}

// Normalized conditional survivor of a prepared target row at latent z.  The
// caller supplies log_q once per (row, z); a numerically hopeless node (or an
// invalid prepared row) reports ok = false and the node contributes zero
// integrand mass rather than failing the trial.
static inline double lrcap_cond_survivor(const BAwLPreparedRow& prep, double z,
                                         double log_q, bool posdrift,
                                         bool& ok) {
  if (prep.status == BAwLTimeStatus::invalid) { ok = false; return 0.0; }
  const double vq = prep.v0 + prep.slope * z;
  const double ls = bawl_prepared_log_survivor(prep, vq, posdrift);
  if (ISNAN(ls)) { ok = false; return 0.0; }
  if (ls == R_NegInf) return 0.0;
  double S = std::exp(ls - log_q);
  if (!R_FINITE(S)) { ok = false; return 0.0; }
  if (S < 0.0) S = 0.0;
  if (S > 1.0) S = 1.0;
  return S;
}

// Normalized conditional density of a prepared target row at latent z.
static inline double lrcap_cond_pdf(const BAwLPreparedRow& prep, double z,
                                    double log_q, bool& ok) {
  if (prep.status == BAwLTimeStatus::invalid) { ok = false; return 0.0; }
  const double vq = prep.v0 + prep.slope * z;
  const double lf = bawl_prepared_log_pdf(prep, vq);
  if (ISNAN(lf)) { ok = false; return 0.0; }
  if (lf == R_NegInf) return 0.0;
  const double f = std::exp(lf - log_q);
  if (!R_FINITE(f) || f < 0.0) { ok = false; return 0.0; }
  return f;
}

// One channel (loaded target vs independent nontarget) prepared at one
// evaluation time.  Everything except the conditional target mean is
// node-invariant: the independent racer's density/survivor, the prepared
// target geometry at the endpoint, and — when the channel-no CDF is needed —
// the GL abscissa weights wt_i*h*f_ind(s_i) with per-abscissa prepared target
// geometry (plan phase 4's two mandatory caches).
struct LrCapChannelPoint {
  bool built = false;
  bool at_zero = false;              // tt <= 0: channel still undecided
  BAwLPreparedRow prep_t;
  double f_ind_t = 0.0;              // independent density at tt
  double S_ind_t = 1.0;              // independent survivor at tt
  int n_gl = 0;
  std::array<double, LRCAP_N_GL> gl_w{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prep{};
};

static void lrcap_build_channel_point(LrCapChannelPoint& cp,
                                      const LrCapTargetView& tv,
                                      const double* par_ind, double tt,
                                      bool need_G, bool need_f_ind,
                                      RacePdf1Fun pdf1, RaceCdf1Fun cdf1,
                                      ContextForRaceModels* ctx,
                                      double min_surv) {
  cp = LrCapChannelPoint();
  cp.built = true;
  if (!(tt > 0.0)) { cp.at_zero = true; return; }
  cp.prep_t = lrcap_prepare_at(tv, tt);
  cp.S_ind_t = safe_surv_at_race_scalar(par_ind, tt, min_surv, cdf1, ctx);
  if (need_f_ind) {
    const double f = pdf1(tt, par_ind, ctx);
    cp.f_ind_t = (R_FINITE(f) && f > 0.0) ? f : 0.0;
  }
  if (need_G) {
    const int t0_idx = ctx->t0_index;
    const double t0n = (t0_idx >= 0) ? par_ind[t0_idx] : 0.0;
    const double lo = std::max(0.0, t0n);
    const double h = (tt > lo) ? 0.5 * (tt - lo) : 0.0;
    if (h > 0.0) {
      const GLRule& gl = gl_get_rule(LRCAP_N_GL);
      cp.n_gl = LRCAP_N_GL;
      for (int i = 0; i < LRCAP_N_GL; ++i) {
        const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
        const double f = pdf1(s, par_ind, ctx);
        cp.gl_w[static_cast<size_t>(i)] =
          (R_FINITE(f) && f > 0.0) ? gl.w[static_cast<size_t>(i)] * h * f : 0.0;
        cp.gl_prep[static_cast<size_t>(i)] = lrcap_prepare_at(tv, s);
      }
    }
  }
}

// Complete conditional channel quantities at the point's time given latent z.
struct LrCapChanVals {
  double f_t = 0.0;   // conditional target density
  double S_t = 0.0;   // conditional target survivor
  double f_n = 0.0;   // independent density
  double S_n = 1.0;   // independent survivor
  double G_no = 0.0;  // P(nontarget won the channel by tt | z)
  bool ok = true;
};

static LrCapChanVals lrcap_chan_vals(const LrCapChannelPoint& cp, double z,
                                     double log_q, bool posdrift,
                                     bool need_pdf, bool count) {
  LrCapChanVals out;
  if (cp.at_zero) { out.S_t = 1.0; return out; }
  out.S_t = lrcap_cond_survivor(cp.prep_t, z, log_q, posdrift, out.ok);
  if (!out.ok) return out;
  out.f_n = cp.f_ind_t;
  out.S_n = cp.S_ind_t;
  if (need_pdf) {
    out.f_t = lrcap_cond_pdf(cp.prep_t, z, log_q, out.ok);
    if (!out.ok) return out;
  }
  double G_no = 0.0;
  for (int i = 0; i < cp.n_gl; ++i) {
    const double w = cp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    G_no += w * lrcap_cond_survivor(cp.gl_prep[static_cast<size_t>(i)], z,
                                    log_q, posdrift, out.ok);
    if (!out.ok) return out;
  }
  if (count && cp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += cp.n_gl;
  out.G_no = std::min(1.0, std::max(0.0, G_no));
  return out;
}

static inline LrChannelState lrcap_state_from_vals(const LrCapChanVals& v) {
  LrChannelState st;
  st.S_dec = v.S_t * v.S_n;
  st.G_no = v.G_no;
  st.G_yes = std::min(1.0, std::max(0.0, 1.0 - st.S_dec - st.G_no));
  return st;
}

// One detection evaluation point: prepared A/B target rows plus (for GNG) the
// nogo survivor and the z-dependent GL integral caches.  `gl_w` carries
// wt_i*h*f_N(s_i) for the nogo-win mass, or wt_i*h*S_N(s_i) for the go-win
// window mass, depending on which evaluator consumes the point.
struct LrCapDetPoint {
  bool built = false;
  bool at_zero = false;
  BAwLPreparedRow prepA, prepB;
  double S_N = 1.0;
  int n_gl = 0;
  std::array<double, LRCAP_N_GL> gl_w{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prepA{};
  std::array<BAwLPreparedRow, LRCAP_N_GL> gl_prepB{};
};

// N(tt) point for OR_DETECTION_ANALYTIC / OR_DETECTION_GNG: survivors at tt,
// and for GNG the nogo-win mass cache over [max(0, t0_N), tt].
static void lrcap_build_det_N_point(LrCapDetPoint& dp,
                                    const LrCapTargetView& tvA,
                                    const LrCapTargetView& tvB,
                                    const double* parN, bool has_nogo,
                                    double tt, bool need_nogo_win,
                                    RacePdf1Fun pdf1, RaceCdf1Fun cdf1,
                                    ContextForRaceModels* ctx,
                                    double min_surv) {
  dp = LrCapDetPoint();
  dp.built = true;
  if (!(tt > 0.0)) { dp.at_zero = true; return; }
  dp.prepA = lrcap_prepare_at(tvA, tt);
  dp.prepB = lrcap_prepare_at(tvB, tt);
  if (has_nogo) {
    dp.S_N = safe_surv_at_race_scalar(parN, tt, min_surv, cdf1, ctx);
    if (need_nogo_win) {
      const int t0_idx = ctx->t0_index;
      const double t0n = (t0_idx >= 0) ? parN[t0_idx] : 0.0;
      const double lo = std::max(0.0, t0n);
      const double h = (tt > lo) ? 0.5 * (tt - lo) : 0.0;
      if (h > 0.0) {
        const GLRule& gl = gl_get_rule(LRCAP_N_GL);
        dp.n_gl = LRCAP_N_GL;
        for (int i = 0; i < LRCAP_N_GL; ++i) {
          const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
          const double f = pdf1(s, parN, ctx);
          dp.gl_w[static_cast<size_t>(i)] =
            (R_FINITE(f) && f > 0.0) ? gl.w[static_cast<size_t>(i)] * h * f : 0.0;
          dp.gl_prepA[static_cast<size_t>(i)] = lrcap_prepare_at(tvA, s);
          dp.gl_prepB[static_cast<size_t>(i)] = lrcap_prepare_at(tvB, s);
        }
      }
    }
  }
}

// Go-win window point for the GNG lower-censor mass over [lo, hi]:
// integral of [f_A(s|z) S_B(s|z) + f_B(s|z) S_A(s|z)] * S_N(s) ds.
static void lrcap_build_det_gowin_point(LrCapDetPoint& dp,
                                        const LrCapTargetView& tvA,
                                        const LrCapTargetView& tvB,
                                        const double* parN, bool has_nogo,
                                        double lo, double hi,
                                        RaceCdf1Fun cdf1,
                                        ContextForRaceModels* ctx,
                                        double min_surv) {
  dp = LrCapDetPoint();
  dp.built = true;
  const double h = (hi > lo) ? 0.5 * (hi - lo) : 0.0;
  if (!(h > 0.0)) { dp.at_zero = true; return; }
  const GLRule& gl = gl_get_rule(LRCAP_N_GL);
  dp.n_gl = LRCAP_N_GL;
  for (int i = 0; i < LRCAP_N_GL; ++i) {
    const double s = lo + h * (1.0 + gl.x[static_cast<size_t>(i)]);
    const double SN = has_nogo
      ? safe_surv_at_race_scalar(parN, s, min_surv, cdf1, ctx) : 1.0;
    dp.gl_w[static_cast<size_t>(i)] = gl.w[static_cast<size_t>(i)] * h * SN;
    dp.gl_prepA[static_cast<size_t>(i)] = lrcap_prepare_at(tvA, s);
    dp.gl_prepB[static_cast<size_t>(i)] = lrcap_prepare_at(tvB, s);
  }
}

// Conditional detection no-response probability at a prepared N point.
static double lrcap_det_N_at(const LrCapDetPoint& dp, int rule_code, double z,
                             double lqA, double lqB, bool posdrift, bool& ok,
                             bool count) {
  if (dp.at_zero) return 1.0;
  const double SA = lrcap_cond_survivor(dp.prepA, z, lqA, posdrift, ok);
  if (!ok) return 1.0;
  const double SB = lrcap_cond_survivor(dp.prepB, z, lqB, posdrift, ok);
  if (!ok) return 1.0;
  if (rule_code == 5) return SA * SB;
  double p = SA * SB * dp.S_N;
  double nogo_win = 0.0;
  for (int i = 0; i < dp.n_gl; ++i) {
    const double w = dp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    const double SAi = lrcap_cond_survivor(dp.gl_prepA[static_cast<size_t>(i)],
                                           z, lqA, posdrift, ok);
    if (!ok) return 1.0;
    const double SBi = lrcap_cond_survivor(dp.gl_prepB[static_cast<size_t>(i)],
                                           z, lqB, posdrift, ok);
    if (!ok) return 1.0;
    nogo_win += w * SAi * SBi;
  }
  if (count && dp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += 2LL * dp.n_gl;
  return std::min(1.0, p + nogo_win);
}

// Conditional go-win window mass at a prepared go-win point.
static double lrcap_det_gowin_at(const LrCapDetPoint& dp, double z,
                                 double lqA, double lqB, bool posdrift,
                                 bool& ok, bool count) {
  if (dp.at_zero) return 0.0;
  double mass = 0.0;
  for (int i = 0; i < dp.n_gl; ++i) {
    const double w = dp.gl_w[static_cast<size_t>(i)];
    if (w == 0.0) continue;
    const double fAi = lrcap_cond_pdf(dp.gl_prepA[static_cast<size_t>(i)], z,
                                      lqA, ok);
    if (!ok) return 0.0;
    const double SAi = lrcap_cond_survivor(dp.gl_prepA[static_cast<size_t>(i)],
                                           z, lqA, posdrift, ok);
    if (!ok) return 0.0;
    const double fBi = lrcap_cond_pdf(dp.gl_prepB[static_cast<size_t>(i)], z,
                                      lqB, ok);
    if (!ok) return 0.0;
    const double SBi = lrcap_cond_survivor(dp.gl_prepB[static_cast<size_t>(i)],
                                           z, lqB, posdrift, ok);
    if (!ok) return 0.0;
    mass += w * (fAi * SBi + fBi * SAi);
  }
  if (count && dp.n_gl > 0)
    lr_capacity_counters().channel_gl_evaluations += 2LL * dp.n_gl;
  return mass;
}

// Shared-factor integration: a fixed scan pass locates each integrand's mass
// (and is the fallback value), then a recentred pass refines it
// (gh_quad.h; scan-based centring per the capacity plan, phase 5).
template <typename LogMass>
static double lrcap_integrate_factor(LogMass&& log_mass, int n_scan,
                                     int n_fine, bool count,
                                     bool bounded_mass = false) {
  // Censoring masses and truncation-window masses are probabilities.  A
  // quadrature approximation with a positive log mass is therefore a
  // numerical failure, not a plausible large likelihood.  This is important
  // for the recentered pass: a narrow positive-drift integrand can be missed
  // by the scan and then badly over-weighted by the refined rule.
  static constexpr double kLogProbabilityTol = 1e-8;
  const double sqrt2 = std::sqrt(2.0);
  const GHRule& scan_rule = gh_rule(n_scan);
  std::array<double, 256> g{};
  std::array<double, 256> zv{};
  double scan_lse = R_NegInf;
  for (int q = 0; q < n_scan; ++q) {
    const double z = sqrt2 * scan_rule.x[static_cast<size_t>(q)];
    zv[static_cast<size_t>(q)] = z;
    const double w = gh_standard_normal_weight(scan_rule, q);
    const double v = (R_FINITE(w) && w > 0.0)
      ? std::log(w) + log_mass(z) : R_NegInf;
    g[static_cast<size_t>(q)] = v;
    scan_lse = log_sum_exp(scan_lse, v);
  }
  if (count) lr_capacity_counters().factor_node_evaluations += n_scan;
  if (!R_FINITE(scan_lse)) return R_NegInf;
  const AGHCenter c = agh_center_from_scan(g.data(), zv.data(), n_scan);
  const GHRule& fine_rule = gh_rule(n_fine);
  double fine_lse = R_NegInf;
  for (int q = 0; q < n_fine; ++q) {
    const double w = fine_rule.w[static_cast<size_t>(q)];
    if (!(R_FINITE(w) && w > 0.0)) continue;
    const double lgw = std::log(w) +
      fine_rule.x[static_cast<size_t>(q)] * fine_rule.x[static_cast<size_t>(q)];
    const double z = c.mu + c.sigma * sqrt2 * fine_rule.x[static_cast<size_t>(q)];
    fine_lse = log_sum_exp(fine_lse,
                           agh_log_weight(lgw, c.sigma, z) + log_mass(z));
  }
  if (count) lr_capacity_counters().factor_node_evaluations += n_fine;
  if (bounded_mass) {
    const bool scan_ok = R_FINITE(scan_lse) &&
      scan_lse <= kLogProbabilityTol;
    const bool fine_ok = R_FINITE(fine_lse) &&
      fine_lse <= kLogProbabilityTol;
    if (fine_ok) return fine_lse;
    // The fixed scan is a deliberately conservative fallback when the
    // recentered rule violates the probability bound.  If both rules fail,
    // report an invalid mass so the caller floors this particle/trial rather
    // than exposing an artificial positive likelihood to the sampler.
    return scan_ok ? scan_lse : R_NegInf;
  }
  return R_FINITE(fine_lse) ? fine_lse : scan_lse;
}

// Complete per-trial likelihood for one capacity (AB, kappa/tau active)
// trial.  Mirrors the ordinary logical-rules branches conditional on the
// factor: the rule is assembled inside the integral, censoring masses and the
// truncation normaliser are z-integrated, and the final contract is
// log num - log den with den the closed-form q_AB (untruncated posdrift) or
// the z-integrated truncation window.
static double lrcap_trial_ll(
    int rule_code, int cond, int resp,
    double t, double LTj, double UTj, double LCj, double UCj,
    double kappa, double tau,
    const double* const* pars_cols, int n_par,
    int idxA, int idxB,
    const double* parnA, const double* parnB, const double* parN,
    bool has_channels, bool has_nogo,
    double min_ll, double min_surv,
    RacePdf1Fun pdf1, RaceCdf1Fun cdf1, ContextForRaceModels* ctx,
    int n_scan, int n_fine, bool count,
    int cell_idx, std::vector<double>* cell_log_den_cap,
    const GslIntegrationControls& gsl_ctl,
    gsl_integration_workspace* w,
    std::vector<double>& gng_scratch) {
  if (cond != 3) return min_ll;  // capacity trials are AB by construction
  const bool posdrift = ctx->use_posdrift;
  // GNG (rule 6) is the full four-horse OR task whose "no" outcome is a
  // withheld response (rt = +Inf); it shares OR's finite-RT go density but
  // routes its censor/withheld/truncation masses through the conditional
  // "no overt go by t" probability N_GNG(t | z) = Q_A(t|z) Q_B(t|z).
  const bool is_gng = (rule_code == 6);
  const LrCapTargetView tvA = lrcap_target_view(pars_cols, idxA, kappa, tau);
  const LrCapTargetView tvB = lrcap_target_view(pars_cols, idxB, kappa, tau);
  if (!tvA.valid || !tvB.valid) return min_ll;
  const bool finite_rt = R_FINITE(t) && t > 0.0;
  const bool upper_censored = !finite_rt && t == R_PosInf;
  const bool lower_censored = !finite_rt && t == R_NegInf;
  if (!finite_rt && !upper_censored && !lower_censored) return min_ll;

  // Closed-form joint positivity normalizer for the pair (untruncated
  // denominator); the quadrature identity int qA(z) qB(z) phi(z) dz is kept
  // only as a test cross-check.
  double log_qAB = 0.0;
  if (posdrift) {
    if (tau == 0.0) {
      log_qAB = pnorm_log_direct(tvA.v0 / tvA.sv, true) +
        pnorm_log_direct(tvB.v0 / tvB.sv, true);
    } else {
      const double sdA = std::sqrt(tvA.sv * tvA.sv + tvA.slope * tvA.slope);
      const double sdB = std::sqrt(tvB.sv * tvB.sv + tvB.slope * tvB.slope);
      const double rho = (tvA.slope * tvB.slope) / (sdA * sdB);
      const double qAB = bawl_corr_pair_positive_normalizer(
          tvA.v0, sdA, tvB.v0, sdB, rho, true);
      if (!R_FINITE(qAB) || !(qAB > 0.0)) return min_ll;
      log_qAB = std::log(std::fmin(qAB, 1.0));
    }
    if (!R_FINITE(log_qAB)) return min_ll;
  }

  const bool is_detection = (rule_code == 5);
  const bool need_G_for_N = (rule_code == 1 || rule_code == 2);
  const bool has_trunc = (LTj != 0.0 || R_FINITE(UTj));
  const bool subtract_UT = R_FINITE(UTj) && !ctx->defective_upper_tail;

  // The capacity denominator is independent of the observed RT. Cache it
  // by cell_idx, leaving the RT-specific event numerator on the per-trial path.
  double cached_log_den = R_NegInf;
  bool z_cache_hit = false;
  if (has_trunc && cell_log_den_cap != nullptr && cell_idx >= 0 &&
      cell_idx < static_cast<int>(cell_log_den_cap->size())) {
    if (!ISNAN((*cell_log_den_cap)[cell_idx])) {
      z_cache_hit = true;
      cached_log_den = (*cell_log_den_cap)[cell_idx];
      if (cached_log_den == R_NegInf) return min_ll;
    }
  }

  // Response gates mirror the ordinary branches; an incompatible recorded
  // response is zero mass for every factor value.
  if (is_detection) {
    if (finite_rt && resp != 1) return min_ll;
    if (upper_censored && resp != 0 && resp != 2) return min_ll;
    if (lower_censored && resp != 0 && resp != 1) return min_ll;
    if (upper_censored && !R_FINITE(UCj)) {
      Rcpp::stop("LogicalRules OR_DETECTION_ANALYTIC requires finite UC for upper-censored trials.");
    }
  } else if (is_gng) {
    // Overt (go) responses are "yes"; a withheld response is rt = +Inf with a
    // missing or "no" identity. UC may be +Inf (the withheld mass extends to
    // infinity when there is no upper deadline).
    if (finite_rt && resp != 1) return min_ll;
    if (upper_censored && resp != 0 && resp != 2) return min_ll;
    if (lower_censored && resp != 0 && resp != 1) return min_ll;
  } else {
    if (!has_channels) return min_ll;
    if (upper_censored && !R_FINITE(UCj)) return min_ll;
  }

  // ---- Node-invariant evaluation-point caches -----------------------------
  // Detection points.
  LrCapDetPoint det_t, det_uc, det_ut, det_lt, det_clo, det_chi, det_gowin;
  // Choice channel points.
  LrCapChannelPoint A_t, B_t, A_uc, B_uc, A_ut, B_ut, A_lt, B_lt,
    A_clo, B_clo, A_chi, B_chi;

  const double clo = std::max(0.0, LTj);
  const double chi = std::max(clo, LCj);

  if (is_detection) {
    if (finite_rt) {
      lrcap_build_det_N_point(det_t, tvA, tvB, parN, has_nogo, t, false,
                              pdf1, cdf1, ctx, min_surv);
    } else if (upper_censored) {
      lrcap_build_det_N_point(det_uc, tvA, tvB, parN, has_nogo, UCj,
                              rule_code == 6, pdf1, cdf1, ctx, min_surv);
    } else {
      if (rule_code == 5) {
        lrcap_build_det_N_point(det_clo, tvA, tvB, parN, false, clo, false,
                                pdf1, cdf1, ctx, min_surv);
        lrcap_build_det_N_point(det_chi, tvA, tvB, parN, false, chi, false,
                                pdf1, cdf1, ctx, min_surv);
      } else {
        lrcap_build_det_gowin_point(det_gowin, tvA, tvB, parN, has_nogo,
                                    clo, chi, cdf1, ctx, min_surv);
      }
    }
    const bool need_ut_for_event = upper_censored && subtract_UT;
    if (has_trunc && (!z_cache_hit || need_ut_for_event)) {
      if (!z_cache_hit && LTj > 0.0) {
        lrcap_build_det_N_point(det_lt, tvA, tvB, parN, has_nogo, LTj,
                                rule_code == 6, pdf1, cdf1, ctx, min_surv);
      }
      if (R_FINITE(UTj) && (!z_cache_hit || need_ut_for_event)) {
        lrcap_build_det_N_point(det_ut, tvA, tvB, parN, has_nogo, UTj,
                                rule_code == 6, pdf1, cdf1, ctx, min_surv);
      }
    }
  } else {
    // GNG builds only the finite-RT go-density channel points; its censor,
    // withheld and truncation masses are conditional "no overt go by t"
    // probabilities evaluated on the fly by gng_N_at (which handles UC = +Inf).
    if (finite_rt) {
      lrcap_build_channel_point(A_t, tvA, parnA, t, true, true,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_t, tvB, parnB, t, true, true,
                                pdf1, cdf1, ctx, min_surv);
    } else if (!is_gng && upper_censored) {
      lrcap_build_channel_point(A_uc, tvA, parnA, UCj, need_G_for_N, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_uc, tvB, parnB, UCj, need_G_for_N, false,
                                pdf1, cdf1, ctx, min_surv);
    } else if (!is_gng) {
      lrcap_build_channel_point(A_clo, tvA, parnA, clo, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_clo, tvB, parnB, clo, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(A_chi, tvA, parnA, chi, true, false,
                                pdf1, cdf1, ctx, min_surv);
      lrcap_build_channel_point(B_chi, tvB, parnB, chi, true, false,
                                pdf1, cdf1, ctx, min_surv);
    }
    const bool need_ut_for_event = upper_censored && subtract_UT;
    if (!is_gng && has_trunc && (!z_cache_hit || need_ut_for_event)) {
      if (!z_cache_hit && LTj > 0.0) {
        lrcap_build_channel_point(A_lt, tvA, parnA, LTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
        lrcap_build_channel_point(B_lt, tvB, parnB, LTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
      }
      if (R_FINITE(UTj) && (!z_cache_hit || need_ut_for_event)) {
        lrcap_build_channel_point(A_ut, tvA, parnA, UTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
        lrcap_build_channel_point(B_ut, tvB, parnB, UTj, need_G_for_N, false,
                                  pdf1, cdf1, ctx, min_surv);
      }
    }
  }

  // Conditional "no overt go response by tt" probability for GNG at latent z,
  // N_GNG(tt | z) = Q_A(tt|z) Q_B(tt|z), Q_X = 1 - GX_yes(tt|z).  The
  // conditional target drift mean is materialized per node (v -> v0 + slope*z)
  // and the channel win probability is integrated with the shared scalar GSL
  // helper, which also handles tt = +Inf (the withheld-to-infinity case).
  // See c_log_likelihood_logicalrules for why these are thread_local vectors
  // rather than fixed arrays, and why they are padded past n_par.
  const size_t lrcap_row_width = static_cast<size_t>(std::max(n_par, 64));
  static thread_local std::vector<double> bufA_row_tls;
  static thread_local std::vector<double> bufB_row_tls;
  bufA_row_tls.assign(lrcap_row_width, 0.0);
  bufB_row_tls.assign(lrcap_row_width, 0.0);
  double* bufA_row = bufA_row_tls.data();
  double* bufB_row = bufB_row_tls.data();
  int gng_isok2[2] = {1, 1};
  // gng_N_at is used both for the censor/withheld event masses (!finite_rt) and
  // for the truncation-window denominator (any finite-RT GNG trial with LT/UT),
  // so initialize its per-node buffers whenever the trial is GNG.
  if (is_gng) {
    for (int p = 0; p < n_par; ++p) {
      bufA_row[p] = pars_cols[p][idxA];
      bufB_row[p] = pars_cols[p][idxB];
    }
    gng_scratch.resize(static_cast<size_t>(2 * n_par));
  }
  auto gng_N_at = [&](double tt, double z, bool& ok) -> double {
    ok = true;
    if (!(tt > 0.0)) return 1.0;
    bufA_row[emc2col::lba::v] = tvA.v0 + tvA.slope * z;
    bufB_row[emc2col::lba::v] = tvB.v0 + tvB.slope * z;
    const double gA = lr_channel_yes_prob(tt, bufA_row, parnA, n_par, ctx,
                                          pdf1, cdf1, gsl_ctl, w,
                                          gng_scratch.data(), gng_isok2);
    const double gB = lr_channel_yes_prob(tt, bufB_row, parnB, n_par, ctx,
                                          pdf1, cdf1, gsl_ctl, w,
                                          gng_scratch.data(), gng_isok2);
    if (!R_FINITE(gA) || !R_FINITE(gB)) { ok = false; return 1.0; }
    const double QA = std::min(1.0, std::max(0.0, 1.0 - gA));
    const double QB = std::min(1.0, std::max(0.0, 1.0 - gB));
    return QA * QB;
  };

  // Conditional no-response probability for the choice rules at a prepared
  // channel-point pair (XOR/ID need only the channel survivors).
  auto choice_N_at = [&](const LrCapChannelPoint& cpA,
                         const LrCapChannelPoint& cpB, double z,
                         double lqA, double lqB, bool& ok) -> double {
    if (cpA.at_zero && cpB.at_zero) return 1.0;
    const LrCapChanVals va = lrcap_chan_vals(cpA, z, lqA, posdrift, false, count);
    if (!va.ok) { ok = false; return 1.0; }
    const LrCapChanVals vb = lrcap_chan_vals(cpB, z, lqB, posdrift, false, count);
    if (!vb.ok) { ok = false; return 1.0; }
    if (rule_code == 3 || rule_code == 4) {
      const double SdA = va.S_t * va.S_n;
      const double SdB = vb.S_t * vb.S_n;
      return 1.0 - (1.0 - SdA) * (1.0 - SdB);
    }
    return lr_no_response_prob(rule_code, lrcap_state_from_vals(va),
                               lrcap_state_from_vals(vb));
  };

  // ---- Conditional event mass (numerator integrand, includes the node
  // positivity weight qA(z)*qB(z)) ------------------------------------------
  auto log_event_mass = [&](double z) -> double {
    const double lqA = lrcap_log_q(tvA, z, posdrift);
    const double lqB = lrcap_log_q(tvB, z, posdrift);
    if (!R_FINITE(lqA) || !R_FINITE(lqB)) return R_NegInf;
    bool ok = true;
    double p = 0.0;
    if (is_detection) {
      if (finite_rt) {
        const double fA = lrcap_cond_pdf(det_t.prepA, z, lqA, ok);
        if (!ok) return R_NegInf;
        const double SA = lrcap_cond_survivor(det_t.prepA, z, lqA, posdrift, ok);
        if (!ok) return R_NegInf;
        const double fB = lrcap_cond_pdf(det_t.prepB, z, lqB, ok);
        if (!ok) return R_NegInf;
        const double SB = lrcap_cond_survivor(det_t.prepB, z, lqB, posdrift, ok);
        if (!ok) return R_NegInf;
        p = fA * std::max(min_surv, SB) + fB * std::max(min_surv, SA);
        if (rule_code == 6) p *= det_t.S_N;
      } else if (upper_censored) {
        p = lrcap_det_N_at(det_uc, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (ok && subtract_UT) {
          const double N_UT = lrcap_det_N_at(det_ut, rule_code, z, lqA, lqB,
                                             posdrift, ok, count);
          p = std::max(0.0, p - N_UT);
        }
        if (!ok) return R_NegInf;
      } else {
        if (rule_code == 5) {
          const double N_lo = lrcap_det_N_at(det_clo, 5, z, lqA, lqB,
                                             posdrift, ok, count);
          if (!ok) return R_NegInf;
          const double N_hi = lrcap_det_N_at(det_chi, 5, z, lqA, lqB,
                                             posdrift, ok, count);
          if (!ok) return R_NegInf;
          p = std::max(0.0, N_lo - N_hi);
        } else {
          p = lrcap_det_gowin_at(det_gowin, z, lqA, lqB, posdrift, ok, count);
          if (!ok) return R_NegInf;
        }
      }
    } else {
      if (finite_rt) {
        const LrCapChanVals va = lrcap_chan_vals(A_t, z, lqA, posdrift, true, count);
        if (!va.ok) return R_NegInf;
        const LrCapChanVals vb = lrcap_chan_vals(B_t, z, lqB, posdrift, true, count);
        if (!vb.ok) return R_NegInf;
        const double one_m_FA = std::max(min_surv, va.S_t);
        const double one_m_FB = std::max(min_surv, vb.S_t);
        const double one_m_FnA = std::max(min_surv, va.S_n);
        const double one_m_FnB = std::max(min_surv, vb.S_n);
        const double gA_yes = va.f_t * one_m_FnA;
        const double gB_yes = vb.f_t * one_m_FnB;
        const double gA_no = va.f_n * one_m_FA;
        const double gB_no = vb.f_n * one_m_FB;
        const double S_dec_A = one_m_FA * one_m_FnA;
        const double S_dec_B = one_m_FB * one_m_FnB;
        const double GA_no = va.G_no;
        const double GB_no = vb.G_no;
        if (rule_code == 1 || is_gng) {
          // GNG shares OR's go-response density (first channel-yes); its "no"
          // outcome is withheld (rt = +Inf), never a finite-RT event, so only
          // resp == 1 (gated above) contributes here.
          const double s_GA_yes = std::min(1.0, std::max(min_surv, GA_no + S_dec_A));
          const double s_GB_yes = std::min(1.0, std::max(min_surv, GB_no + S_dec_B));
          if (resp == 1) p = gA_yes * s_GB_yes + gB_yes * s_GA_yes;
          else if (resp == 2 && !is_gng) p = gA_no * GB_no + gB_no * GA_no;
        } else if (rule_code == 2) {
          const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
          const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
          const double s_GA_no = std::min(1.0, std::max(min_surv, GA_yes + S_dec_A));
          const double s_GB_no = std::min(1.0, std::max(min_surv, GB_yes + S_dec_B));
          if (resp == 1) p = gA_yes * GB_yes + gB_yes * GA_yes;
          else if (resp == 2) p = gA_no * s_GB_no + gB_no * s_GA_no;
        } else {
          const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
          const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
          const double dAB = gA_yes * GB_yes + gB_yes * GA_yes;
          const double dAN = gA_yes * GB_no + gB_no * GA_yes;
          const double dNB = gA_no * GB_yes + gB_yes * GA_no;
          const double dNN = gA_no * GB_no + gB_no * GA_no;
          if (rule_code == 4) {
            if (resp == 6) p = dAB;
            else if (resp == 4) p = dAN;
            else if (resp == 5) p = dNB;
            else if (resp == 3) p = dNN;
          } else {
            if (resp == 1) p = dAN + dNB;
            else if (resp == 2) p = dAB + dNN;
          }
        }
      } else if (is_gng) {
        // Withheld (rt = +Inf): no overt go by UC, N_GNG(UC|z) minus N_GNG(UT|z)
        // when the window has a finite upper bound; lower censor: go mass in
        // [clo, chi] = N_GNG(clo|z) - N_GNG(chi|z).
        if (upper_censored) {
          p = gng_N_at(UCj, z, ok);
          if (ok && subtract_UT) {
            const double N_UT = gng_N_at(UTj, z, ok);
            p = std::max(0.0, p - N_UT);
          }
          if (!ok) return R_NegInf;
        } else {
          const double N_lo = gng_N_at(clo, z, ok);
          if (!ok) return R_NegInf;
          const double N_hi = gng_N_at(chi, z, ok);
          if (!ok) return R_NegInf;
          p = std::max(0.0, N_lo - N_hi);
        }
      } else if (upper_censored) {
        p = choice_N_at(A_uc, B_uc, z, lqA, lqB, ok);
        if (ok && subtract_UT) {
          const double N_UT = choice_N_at(A_ut, B_ut, z, lqA, lqB, ok);
          p = std::max(0.0, p - N_UT);
        }
        if (!ok) return R_NegInf;
      } else {
        const LrCapChanVals va_lo = lrcap_chan_vals(A_clo, z, lqA, posdrift, false, count);
        const LrCapChanVals vb_lo = lrcap_chan_vals(B_clo, z, lqB, posdrift, false, count);
        const LrCapChanVals va_hi = lrcap_chan_vals(A_chi, z, lqA, posdrift, false, count);
        const LrCapChanVals vb_hi = lrcap_chan_vals(B_chi, z, lqB, posdrift, false, count);
        if (!va_lo.ok || !vb_lo.ok || !va_hi.ok || !vb_hi.ok) return R_NegInf;
        const double cdf_hi = lr_response_cdf(rule_code, resp,
                                              lrcap_state_from_vals(va_hi),
                                              lrcap_state_from_vals(vb_hi));
        const double cdf_lo = lr_response_cdf(rule_code, resp,
                                              lrcap_state_from_vals(va_lo),
                                              lrcap_state_from_vals(vb_lo));
        if (ISNAN(cdf_hi) || ISNAN(cdf_lo)) return R_NegInf;
        p = std::max(0.0, cdf_hi - cdf_lo);
      }
    }
    if (!(p > 0.0) || !R_FINITE(p)) return R_NegInf;
    return std::log(p) + lqA + lqB;
  };

  // ---- Conditional truncation-window mass (denominator integrand) ---------
  auto log_trunc_mass = [&](double z) -> double {
    const double lqA = lrcap_log_q(tvA, z, posdrift);
    const double lqB = lrcap_log_q(tvB, z, posdrift);
    if (!R_FINITE(lqA) || !R_FINITE(lqB)) return R_NegInf;
    bool ok = true;
    double N_LT = 1.0;
    double N_UT = 0.0;
    if (is_detection) {
      if (LTj > 0.0) {
        N_LT = lrcap_det_N_at(det_lt, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = lrcap_det_N_at(det_ut, rule_code, z, lqA, lqB, posdrift, ok, count);
        if (!ok) return R_NegInf;
      }
    } else if (is_gng) {
      if (LTj > 0.0) {
        N_LT = gng_N_at(LTj, z, ok);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = gng_N_at(UTj, z, ok);
        if (!ok) return R_NegInf;
      }
    } else {
      if (LTj > 0.0) {
        N_LT = choice_N_at(A_lt, B_lt, z, lqA, lqB, ok);
        if (!ok) return R_NegInf;
      }
      if (R_FINITE(UTj)) {
        N_UT = choice_N_at(A_ut, B_ut, z, lqA, lqB, ok);
        if (!ok) return R_NegInf;
      }
    }
    const double Z = std::min(1.0, N_LT - N_UT);
    if (!(Z > 0.0) || !R_FINITE(Z)) return R_NegInf;
    return std::log(Z) + lqA + lqB;
  };

  double log_num = R_NegInf;
  double log_den = 0.0;
  if (tau == 0.0) {
    // The factor is degenerate: one node at z = 0 is the exact integral, and
    // the positivity weight cancels exactly against the rho = 0 orthant.
    if (count) ++lr_capacity_counters().tau_zero_trials;
    log_num = log_event_mass(0.0);
    const double log_w0 = lrcap_log_q(tvA, 0.0, posdrift) +
      lrcap_log_q(tvB, 0.0, posdrift);
    log_den = has_trunc
      ? (z_cache_hit ? cached_log_den : log_trunc_mass(0.0))
      : log_w0;
  } else {
    // All censored event masses are probabilities; finite-RT event densities
    // are not bounded by one and must retain the ordinary adaptive result.
    log_num = lrcap_integrate_factor(log_event_mass, n_scan, n_fine, count,
                                     !finite_rt);
    log_den = has_trunc
      ? (z_cache_hit
          ? cached_log_den
          : lrcap_integrate_factor(log_trunc_mass, n_scan, n_fine, count,
                                   true))
      : log_qAB;
  }
  if (has_trunc && !z_cache_hit && cell_log_den_cap != nullptr && cell_idx >= 0 &&
      cell_idx < static_cast<int>(cell_log_den_cap->size())) {
    const bool valid_z = R_FINITE(log_den);
    (*cell_log_den_cap)[cell_idx] = valid_z ? log_den : R_NegInf;
  }
  if (!R_FINITE(log_den)) return min_ll;
  const double ll = log_num - log_den;
  if (!finite_rt && (!R_FINITE(log_num) || log_num > 1e-8 ||
                     log_den > 1e-8 || ll > 1e-8)) return min_ll;
  return (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
}

double c_log_likelihood_logicalrules(
    const double* const* pars_cols,
    int n_par,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& ok_params,
    int n_acc,
    ContextForRaceModels* model_ctx,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    RaceRawFun model_dfun_raw,
    RaceRawFun model_pfun_raw,
    const LogicalRulesSharedState& shared,
    Rcpp::NumericVector* trial_ll_out,
    int kappa_col,
    int tau_col,
    int pc_col,
    int pg_col,
    const GuessKernel* guess) {
  if (!shared.valid || n_acc <= 0 || pdf1 == nullptr || cdf1 == nullptr || model_ctx == nullptr) {
    Rcpp::stop("c_log_likelihood_logicalrules: invalid logical-rules configuration.");
  }
  const int n_trials = shared.n_trials;
  if (ok_params.size() != n_trials) {
    Rcpp::stop("c_log_likelihood_logicalrules: ok_params size must match shared data.");
  }
  if (n_trials == 0) return 0.0;
  // ParamTable columns are refilled in place for each particle, so no cache
  // entry may survive this likelihood call.
  shared.clear_particle_cache();
  if (model_ctx->fpe_cache) model_ctx->fpe_cache->new_particle();
  if (model_ctx->rlf_cache) model_ctx->rlf_cache->new_particle();
  // Keep scratch width at least 64 for optional trailing kernel columns.
  const size_t lr_row_width = static_cast<size_t>(std::max(n_par, 64));
  const bool capacity = kappa_col >= 0 || tau_col >= 0;
  if (capacity && (kappa_col < 0 || tau_col < 0 ||
                   kappa_col >= n_par || tau_col >= n_par)) {
    Rcpp::stop("c_log_likelihood_logicalrules: invalid LogicalRules capacity columns.");
  }
  const bool count_capacity = capacity && lr_capacity_counters_enabled();
  const int lrcap_scan_n = capacity
    ? emc2_quad_nodes("EMC2_LRCAP_SCAN_N", 12) : 0;
  const int lrcap_fine_n = capacity
    ? emc2_quad_nodes("EMC2_LRCAP_FINE_N", 12) : 0;

  static const double kMinSurv = 1e-300;
  const int n_unique_trials = shared.n_unique_trials;
  static thread_local LogicalRulesScratch scratch;
  resize_assign_double(scratch.ll_unique, static_cast<size_t>(n_unique_trials), min_ll);
  std::vector<double>& ll_unique = scratch.ll_unique;
  const bool use_raw_local =
    (model_dfun_raw != nullptr &&
     model_pfun_raw != nullptr &&
     static_cast<int>(shared.rt_by_row.size()) == n_trials);

  GslIntegrationControls gsl_ctl = default_gsl_controls();
  gsl_ctl.try_qng_first_finite = true;
  gsl_ctl.qag_key = GSL_INTEG_GAUSS21;
  gsl_ctl.rel_tol = 1e-5;  // 1e-4 let the adaptive rule's optimistic error estimate through

  bool any_ok = false;
  for (int i = 0; i < n_trials && !any_ok; ++i) any_ok = (bool)ok_params[i];
  if (!any_ok) {
    const int n_out = expand.length() > 0 ? (int)expand.length() : n_unique_trials;
    if (trial_ll_out != nullptr) for (int i = 0; i < n_out; ++i) (*trial_ll_out)[i] = min_ll;
    return static_cast<double>(n_out) * min_ll;
  }

  std::vector<double>& logf_all = scratch.logf_all;
  std::vector<double>& logS_all = scratch.logS_all;
  std::vector<int>& all_mask = scratch.all_mask;
  std::vector<int>& isok_int_all = scratch.isok_int_all;
  if (use_raw_local) {
    resize_assign_double(logf_all, static_cast<size_t>(n_trials), min_ll);
    resize_assign_double(logS_all, static_cast<size_t>(n_trials), min_ll);
    resize_assign_int(all_mask, static_cast<size_t>(n_trials), 1);
    isok_int_all.resize(static_cast<size_t>(n_trials));
    for (int i = 0; i < n_trials; ++i) isok_int_all[static_cast<size_t>(i)] = ok_params[i] ? 1 : 0;
    model_dfun_raw(shared.rt_by_row.data(), pars_cols, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logf_all.data(), min_ll, model_ctx);
    model_pfun_raw(shared.rt_by_row.data(), pars_cols, n_trials,
                   all_mask.data(), isok_int_all.data(),
                   logS_all.data(), min_ll, model_ctx);
  }

  // Gauss-Legendre batch integration.
  // GA_no_gl[j] = P(nontarget-A wins channel-A sub-race before RT_j)  [linear probability]
  // GB_no_gl[j] = same for channel B, or = GA_no_gl[j] when channels are equal.
  // OR/AND rules use GA_no_gl directly; XOR/ID rules derive GA_yes = 1 - GA_no_gl - S_dec_A.
  const bool use_gl_pass = use_raw_local;
  std::vector<double>& GA_no_gl = scratch.GA_no_gl;
  std::vector<double>& GB_no_gl = scratch.GB_no_gl;
  std::vector<unsigned char>& ch_eq_vec = scratch.ch_eq_vec;

  if (use_gl_pass) {
    const GLRule& gl_rule31 = gl_get_rule(31);   // shared cache, gl_quad.h

    // Per-trial GL parameters for GA_no = P(n_A wins before RT) and
    // GB_no = P(n_B wins before RT).  The integration lower bound is the
    // density accumulator's t0, not its opponent's t0; otherwise unequal t0
    // pairs lose early mass where the opponent survivor is exactly one.
    const int t0_col = model_ctx->t0_index;
    std::vector<double>& gl_h_A = scratch.gl_h_A;
    std::vector<double>& gl_h_B = scratch.gl_h_B;
    std::vector<double>& gl_lo_A = scratch.gl_lo_A;
    std::vector<double>& gl_lo_B = scratch.gl_lo_B;
    std::vector<unsigned char>& pair_ok_A = scratch.pair_ok_A;
    std::vector<unsigned char>& pair_ok_B = scratch.pair_ok_B;
    gl_h_A.resize(static_cast<size_t>(n_unique_trials));
    gl_h_B.resize(static_cast<size_t>(n_unique_trials));
    gl_lo_A.resize(static_cast<size_t>(n_unique_trials));
    gl_lo_B.resize(static_cast<size_t>(n_unique_trials));
    ch_eq_vec.resize(static_cast<size_t>(n_unique_trials));
    pair_ok_A.resize(static_cast<size_t>(n_unique_trials));
    pair_ok_B.resize(static_cast<size_t>(n_unique_trials));
    bool any_unequal = false;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int iA  = shared.idxA[j],  inA = shared.idxnA[j];
      const int iB  = shared.idxB[j],  inB = shared.idxnB[j];
      const double t0nA = (t0_col >= 0 && t0_col < n_par && inA >= 0)
                        ? pars_cols[t0_col][inA] : 0.0;
      const double t0nB = (t0_col >= 0 && t0_col < n_par && inB >= 0)
                        ? pars_cols[t0_col][inB] : 0.0;
      const double RT = shared.rt_unique[j];
      const double loA = std::max(0.0, t0nA);
      const double loB = std::max(0.0, t0nB);
      gl_lo_A[static_cast<size_t>(j)] = loA;
      gl_lo_B[static_cast<size_t>(j)] = loB;
      gl_h_A[static_cast<size_t>(j)] = (RT > loA) ? (RT - loA) * 0.5 : 0.0;
      gl_h_B[static_cast<size_t>(j)] = (RT > loB) ? (RT - loB) * 0.5 : 0.0;
      ch_eq_vec[j] = (inA >= 0 && inB >= 0) &&
                     row_equal_colmajor(pars_cols, n_par, iA, iB) &&
                     row_equal_colmajor(pars_cols, n_par, inA, inB);
      if (!ch_eq_vec[j]) any_unequal = true;
      pair_ok_A[static_cast<size_t>(j)] =
        (iA >= 0 && inA >= 0 && ok_params[iA] && ok_params[inA] &&
         RT > 0.0 && R_FINITE(RT) && gl_h_A[static_cast<size_t>(j)] > 0.0) ? 1 : 0;
      pair_ok_B[static_cast<size_t>(j)] =
        (iB >= 0 && inB >= 0 && ok_params[iB] && ok_params[inB] &&
         RT > 0.0 && R_FINITE(RT) && gl_h_B[static_cast<size_t>(j)] > 0.0) ? 1 : 0;
    }

    // Compact column-major parameter matrices (n_unique_trials rows) for each accumulator role.
    // Avoids iterating over the full n_trials parameter matrix with sparse masks.
    const size_t cpt = static_cast<size_t>(n_unique_trials * n_par);
    std::vector<double>& pars_nA = scratch.pars_nA;
    std::vector<double>& pars_A = scratch.pars_A;
    std::vector<double>& pars_nB = scratch.pars_nB;
    std::vector<double>& pars_B = scratch.pars_B;
    std::vector<int>& isok_nA = scratch.isok_nA;
    std::vector<int>& isok_A = scratch.isok_A;
    std::vector<int>& isok_nB = scratch.isok_nB;
    std::vector<int>& isok_B = scratch.isok_B;
    pars_nA.resize(cpt);
    pars_A.resize(cpt);
    isok_nA.resize(static_cast<size_t>(n_unique_trials));
    isok_A.resize(static_cast<size_t>(n_unique_trials));
    if (any_unequal) {
      pars_nB.resize(cpt);
      pars_B.resize(cpt);
      isok_nB.resize(static_cast<size_t>(n_unique_trials));
      isok_B.resize(static_cast<size_t>(n_unique_trials));
    }
    for (int p = 0; p < n_par; ++p) {
      const double* src = pars_cols[p];
      double* d_nA = pars_nA.data() + static_cast<size_t>(p) * n_unique_trials;
      double* d_A  = pars_A.data()  + static_cast<size_t>(p) * n_unique_trials;
      for (int j = 0; j < n_unique_trials; ++j) {
        const int inA_j = shared.idxnA[j];
        d_nA[j] = (inA_j >= 0) ? src[inA_j] : 0.0;
        d_A[j]  = src[shared.idxA[j]];
      }
      if (any_unequal) {
        double* d_nB = pars_nB.data() + static_cast<size_t>(p) * n_unique_trials;
        double* d_B  = pars_B.data()  + static_cast<size_t>(p) * n_unique_trials;
        for (int j = 0; j < n_unique_trials; ++j) {
          const int inB_j = shared.idxnB[j];
          d_nB[j] = (inB_j >= 0) ? src[inB_j] : 0.0;
          d_B[j]  = src[shared.idxB[j]];
        }
      }
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      isok_nA[j] = pair_ok_A[static_cast<size_t>(j)] ? isok_int_all[shared.idxnA[j]] : 0;
      isok_A[j]  = pair_ok_A[static_cast<size_t>(j)] ? isok_int_all[shared.idxA[j]] : 0;
      if (any_unequal) {
        isok_nB[j] = (!ch_eq_vec[j] && pair_ok_B[static_cast<size_t>(j)]) ? isok_int_all[shared.idxnB[j]] : 0;
        isok_B[j]  = (!ch_eq_vec[j] && pair_ok_B[static_cast<size_t>(j)]) ? isok_int_all[shared.idxB[j]] : 0;
      }
    }

    // Column-pointer views of the compact role matrices for the raw kernels.
    // Slots beyond n_par stay nullptr: kernels may fetch (not dereference)
    // optional trailing columns this model variant lacks.
    static thread_local std::vector<const double*> colsp_nA_tls, colsp_A_tls,
                                                   colsp_nB_tls, colsp_B_tls;
    colsp_nA_tls.assign(lr_row_width, nullptr);
    colsp_A_tls.assign(lr_row_width, nullptr);
    colsp_nB_tls.assign(lr_row_width, nullptr);
    colsp_B_tls.assign(lr_row_width, nullptr);
    const double** colsp_nA = colsp_nA_tls.data();
    const double** colsp_A  = colsp_A_tls.data();
    const double** colsp_nB = colsp_nB_tls.data();
    const double** colsp_B  = colsp_B_tls.data();
    for (int p = 0; p < n_par; ++p) {
      colsp_nA[p] = pars_nA.data() + static_cast<size_t>(p) * n_unique_trials;
      colsp_A[p]  = pars_A.data()  + static_cast<size_t>(p) * n_unique_trials;
      if (any_unequal) {
        colsp_nB[p] = pars_nB.data() + static_cast<size_t>(p) * n_unique_trials;
        colsp_B[p]  = pars_B.data()  + static_cast<size_t>(p) * n_unique_trials;
      }
    }

    resize_assign_double(GA_no_gl, static_cast<size_t>(n_unique_trials), 0.0);
    resize_assign_double(GB_no_gl, static_cast<size_t>(n_unique_trials), 0.0);

    std::vector<int>& all1 = scratch.all1;
    resize_assign_int(all1, static_cast<size_t>(n_unique_trials), 1);
    std::vector<double>& gl_rt = scratch.tmpA1;
    std::vector<double>& lf_nA = scratch.tmpA2;
    std::vector<double>& lS_A = scratch.tmpB1;
    std::vector<double>& lf_nB = scratch.tmpB2;
    std::vector<double>& lS_B = scratch.parA;
    gl_rt.resize(static_cast<size_t>(n_unique_trials));
    resize_assign_double(lf_nA, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lS_A, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lf_nB, static_cast<size_t>(n_unique_trials), min_ll);
    resize_assign_double(lS_B, static_cast<size_t>(n_unique_trials), min_ll);

    for (size_t k = 0; k < gl_rule31.x.size(); ++k) {
      const double xi = gl_rule31.x[k], wt = gl_rule31.w[k];
      #pragma omp simd
      for (int j = 0; j < n_unique_trials; ++j) {
        gl_rt[static_cast<size_t>(j)] = gl_lo_A[static_cast<size_t>(j)]
                                      + gl_h_A[static_cast<size_t>(j)] * (1.0 + xi);
      }
      model_dfun_raw(gl_rt.data(), colsp_nA, n_unique_trials,
                     all1.data(), isok_nA.data(), lf_nA.data(), min_ll, model_ctx);
      model_pfun_raw(gl_rt.data(), colsp_A,  n_unique_trials,
                     all1.data(), isok_A.data(),  lS_A.data(),  min_ll, model_ctx);
      #pragma omp simd
      for (int j = 0; j < n_unique_trials; ++j)
        GA_no_gl[j] += wt * gl_h_A[j] * std::exp(lf_nA[j] + lS_A[j]);
      if (any_unequal) {
        #pragma omp simd
        for (int j = 0; j < n_unique_trials; ++j) {
          gl_rt[static_cast<size_t>(j)] = gl_lo_B[static_cast<size_t>(j)]
                                        + gl_h_B[static_cast<size_t>(j)] * (1.0 + xi);
        }
        model_dfun_raw(gl_rt.data(), colsp_nB, n_unique_trials,
                       all1.data(), isok_nB.data(), lf_nB.data(), min_ll, model_ctx);
        model_pfun_raw(gl_rt.data(), colsp_B,  n_unique_trials,
                       all1.data(), isok_B.data(),  lS_B.data(),  min_ll, model_ctx);
        #pragma omp simd
        for (int j = 0; j < n_unique_trials; ++j)
          GB_no_gl[j] += wt * gl_h_B[j] * std::exp(lf_nB[j] + lS_B[j]);
      }
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      if (ch_eq_vec[static_cast<size_t>(j)]) GB_no_gl[static_cast<size_t>(j)] = GA_no_gl[static_cast<size_t>(j)];
    }

    // ---- Auxiliary GL sweeps: G_no at censoring/truncation time points ----
    // Same machinery as the RT sweep above, evaluated at the per-trial times
    // the trial loop needs: the truncation normaliser endpoints (LT, UT), the
    // lower-censor window ends, and UC for upper-censored OR/AND trials.
    // Values left NaN are recomputed by the scalar GSL helper in the trial
    // loop, so a slot missed here costs speed, never correctness.
    std::vector<double>& aux_t  = scratch.aux_t;
    std::vector<double>& aux_hA = scratch.aux_hA;
    std::vector<double>& aux_hB = scratch.aux_hB;
    std::vector<double>& aux_GA = scratch.aux_GA;
    std::vector<double>& aux_GB = scratch.aux_GB;
    std::vector<int>& aux_mask  = scratch.aux_mask;
    const size_t aux_len = static_cast<size_t>(LR_AUX_N) * n_unique_trials;
    aux_t.assign(aux_len, NA_REAL);
    aux_hA.assign(aux_len, 0.0);
    aux_hB.assign(aux_len, 0.0);
    aux_GA.assign(aux_len, NA_REAL);
    aux_GB.assign(aux_len, NA_REAL);
    aux_mask.assign(aux_len, 0);
    bool slot_any[LR_AUX_N] = {false, false, false, false, false};

    std::vector<int>& isok_nA_b = scratch.isok_nA_base;
    std::vector<int>& isok_A_b  = scratch.isok_A_base;
    std::vector<int>& isok_nB_b = scratch.isok_nB_base;
    std::vector<int>& isok_B_b  = scratch.isok_B_base;
    resize_assign_int(isok_nA_b, static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_A_b,  static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_nB_b, static_cast<size_t>(n_unique_trials), 0);
    resize_assign_int(isok_B_b,  static_cast<size_t>(n_unique_trials), 0);

    bool any_aux = false;
    for (int j = 0; j < n_unique_trials; ++j) {
      const int rc = shared.rule_code[static_cast<size_t>(j)];
      if (rc == 5) continue;   // analytic detection keeps its own route; GNG (6) uses four horses
      const int iA = shared.idxA[static_cast<size_t>(j)], inA = shared.idxnA[static_cast<size_t>(j)];
      const int iB = shared.idxB[static_cast<size_t>(j)], inB = shared.idxnB[static_cast<size_t>(j)];
      const bool base_ok =
        iA >= 0 && inA >= 0 && iB >= 0 && inB >= 0 &&
        ok_params[iA] && ok_params[inA] && ok_params[iB] && ok_params[inB];
      if (!base_ok) continue;
      isok_nA_b[static_cast<size_t>(j)] = 1;
      isok_A_b[static_cast<size_t>(j)]  = 1;
      if (!ch_eq_vec[static_cast<size_t>(j)]) {
        isok_nB_b[static_cast<size_t>(j)] = 1;
        isok_B_b[static_cast<size_t>(j)]  = 1;
      }

      const bool needs_G = (rc == 1 || rc == 2 || rc == 6);   // XOR/ID N(t) is survivor-only
      const double rtj = shared.rt_unique[static_cast<size_t>(j)];
      const double LTj = shared.LT_unique[static_cast<size_t>(j)];
      const double UTj = shared.UT_unique[static_cast<size_t>(j)];
      const bool truncated = (LTj != 0.0 || R_FINITE(UTj));

      double slot_t[LR_AUX_N] = {NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL};
      if (needs_G && truncated && LTj > 0.0) slot_t[LR_AUX_LT] = LTj;
      if (needs_G && R_FINITE(UTj))          slot_t[LR_AUX_UT] = UTj;
      if (rtj == R_NegInf) {   // lower censor: response CDFs need G for all rules <= 4
        const double lo = std::max(0.0, LTj);
        const double hi = std::max(lo, shared.LC_unique[static_cast<size_t>(j)]);
        if (lo > 0.0) slot_t[LR_AUX_CLO] = lo;
        if (hi > 0.0) slot_t[LR_AUX_CHI] = hi;
      } else if (rtj == R_PosInf && needs_G) {
        const double UCj = shared.UC_unique[static_cast<size_t>(j)];
        if (R_FINITE(UCj)) slot_t[LR_AUX_UC] = UCj;
      }
      for (int s = 0; s < LR_AUX_N; ++s) {
        if (ISNAN(slot_t[s])) continue;
        const size_t o = static_cast<size_t>(s) * n_unique_trials + j;
        aux_t[o] = slot_t[s];
        aux_hA[o] = (slot_t[s] > gl_lo_A[static_cast<size_t>(j)])
                  ? (slot_t[s] - gl_lo_A[static_cast<size_t>(j)]) * 0.5 : 0.0;
        aux_hB[o] = (slot_t[s] > gl_lo_B[static_cast<size_t>(j)])
                  ? (slot_t[s] - gl_lo_B[static_cast<size_t>(j)]) * 0.5 : 0.0;
        aux_mask[o] = 1;
        aux_GA[o] = 0.0;
        aux_GB[o] = 0.0;
        slot_any[s] = true;
        any_aux = true;
      }
    }

    if (any_aux) {
      for (int s = 0; s < LR_AUX_N; ++s) {
        if (!slot_any[s]) continue;
        const double* hA = aux_hA.data() + static_cast<size_t>(s) * n_unique_trials;
        const double* hB = aux_hB.data() + static_cast<size_t>(s) * n_unique_trials;
        const int* msk   = aux_mask.data() + static_cast<size_t>(s) * n_unique_trials;
        double* GA = aux_GA.data() + static_cast<size_t>(s) * n_unique_trials;
        double* GB = aux_GB.data() + static_cast<size_t>(s) * n_unique_trials;
        for (size_t k = 0; k < gl_rule31.x.size(); ++k) {
          const double xi = gl_rule31.x[k], wt = gl_rule31.w[k];
          #pragma omp simd
          for (int j = 0; j < n_unique_trials; ++j) {
            gl_rt[static_cast<size_t>(j)] = gl_lo_A[static_cast<size_t>(j)]
                                          + hA[j] * (1.0 + xi);
          }
          model_dfun_raw(gl_rt.data(), colsp_nA, n_unique_trials,
                         msk, isok_nA_b.data(), lf_nA.data(), min_ll, model_ctx);
          model_pfun_raw(gl_rt.data(), colsp_A,  n_unique_trials,
                         msk, isok_A_b.data(),  lS_A.data(),  min_ll, model_ctx);
          for (int j = 0; j < n_unique_trials; ++j) {
            if (msk[j]) GA[j] += wt * hA[j] * std::exp(lf_nA[j] + lS_A[j]);
          }
          if (any_unequal) {
            #pragma omp simd
            for (int j = 0; j < n_unique_trials; ++j) {
              gl_rt[static_cast<size_t>(j)] = gl_lo_B[static_cast<size_t>(j)]
                                            + hB[j] * (1.0 + xi);
            }
            model_dfun_raw(gl_rt.data(), colsp_nB, n_unique_trials,
                           msk, isok_nB_b.data(), lf_nB.data(), min_ll, model_ctx);
            model_pfun_raw(gl_rt.data(), colsp_B,  n_unique_trials,
                           msk, isok_B_b.data(),  lS_B.data(),  min_ll, model_ctx);
            for (int j = 0; j < n_unique_trials; ++j) {
              if (msk[j] && !ch_eq_vec[static_cast<size_t>(j)]) {
                GB[j] += wt * hB[j] * std::exp(lf_nB[j] + lS_B[j]);
              }
            }
          }
        }
        for (int j = 0; j < n_unique_trials; ++j) {
          if (msk[j] && ch_eq_vec[static_cast<size_t>(j)]) GB[j] = GA[j];
        }
      }
    }
  }

  static thread_local GslWorkspacePtr workspace_tls(nullptr, &gsl_integration_workspace_free);
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace_tls, gsl_ctl.retry_limit);
  std::vector<double>& parA = scratch.parA;
  std::vector<double>& parB = scratch.parB;
  std::vector<double>& parnA = scratch.parnA;
  std::vector<double>& parnB = scratch.parnB;
  std::vector<double>& pars_2buf = scratch.pars_2buf;
  std::vector<double>& parN = scratch.tmpA1;
  std::vector<double>& pars_3buf = scratch.tmpA2;
  parA.resize(static_cast<size_t>(n_par));
  parB.resize(static_cast<size_t>(n_par));
  parnA.resize(static_cast<size_t>(n_par));
  parnB.resize(static_cast<size_t>(n_par));
  parN.resize(static_cast<size_t>(n_par));
  pars_2buf.resize(static_cast<size_t>(2 * n_par));
  int isok_2buf[2] = {1, 1};

  for (int j = 0; j < n_unique_trials; ++j) {
    const int idxA = shared.idxA[static_cast<size_t>(j)];
    const int idxB = shared.idxB[static_cast<size_t>(j)];
    const int idxnA = shared.idxnA[static_cast<size_t>(j)];
    const int idxnB = shared.idxnB[static_cast<size_t>(j)];

    const int rule_code = shared.rule_code[static_cast<size_t>(j)];
    // GNG (rule 6) is the full four-horse OR task whose "no" outcome is a
    // withheld response (rt = +Inf); it shares OR's subrace machinery here.
    const bool is_gng = (rule_code == 6);
    // Parameter-validity gate. OR/AND/XOR/ID and GNG need all four accumulators;
    // the analytic detector (rule 5) needs only valid A and B. Without this,
    // censored branches would evaluate survivors/densities with out-of-bounds
    // parameters.
    const bool pars_valid = (rule_code <= 4 || is_gng)
      ? (idxA >= 0 && idxB >= 0 && idxnA >= 0 && idxnB >= 0 &&
         ok_params[idxA] && ok_params[idxB] && ok_params[idxnA] && ok_params[idxnB])
      : (idxA >= 0 && idxB >= 0 && ok_params[idxA] && ok_params[idxB]);
    if (!pars_valid) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      if (count_capacity) ++lr_capacity_counters().invalid_trials;
      continue;
    }

    const double t = shared.rt_unique[static_cast<size_t>(j)];
    copy_par_row_colmajor(pars_cols, n_par, idxA, parA.data());
    copy_par_row_colmajor(pars_cols, n_par, idxB, parB.data());
    if (idxnA >= 0) copy_par_row_colmajor(pars_cols, n_par, idxnA, parnA.data());
    if (idxnB >= 0) copy_par_row_colmajor(pars_cols, n_par, idxnB, parnB.data());

    const int cell_idx = shared.cell_id[static_cast<size_t>(j)];

    // Capacity applies only to non-degenerate redundant-target trials;
    // other trials use the general evaluator.
    if (capacity) {
      const double kappa = pars_cols[kappa_col][idxA];
      const double tau = pars_cols[tau_col][idxA];
      const bool shared_capacity =
        R_FINITE(kappa) && R_FINITE(tau) && tau >= 0.0 &&
        pars_cols[kappa_col][idxB] == kappa &&
        pars_cols[tau_col][idxB] == tau;
      if (!shared_capacity) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
        if (count_capacity) ++lr_capacity_counters().invalid_trials;
        continue;
      }
      const bool pair_active = shared.cond_code[static_cast<size_t>(j)] == 3 &&
        (kappa != 0.0 || tau != 0.0);
      if (pair_active) {
        if (count_capacity) {
          if (rule_code == 5) ++lr_capacity_counters().capacity_detection_trials;
          else ++lr_capacity_counters().capacity_choice_trials;
        }
        ll_unique[static_cast<size_t>(j)] = lrcap_trial_ll(
          rule_code,
          shared.cond_code[static_cast<size_t>(j)],
          shared.resp_code[static_cast<size_t>(j)],
          t,
          shared.LT_unique[static_cast<size_t>(j)],
          shared.UT_unique[static_cast<size_t>(j)],
          shared.LC_unique[static_cast<size_t>(j)],
          shared.UC_unique[static_cast<size_t>(j)],
          kappa, tau,
          pars_cols, n_par, idxA, idxB,
          parnA.data(), parnB.data(), parN.data(),
          rule_code <= 4 || is_gng, false,
          min_ll, kMinSurv,
          pdf1, cdf1, model_ctx,
          lrcap_scan_n, lrcap_fine_n,
          count_capacity, cell_idx, &shared.cell_log_den_cap,
          gsl_ctl, w, pars_3buf);
        continue;
      }
      if (count_capacity) ++lr_capacity_counters().ordinary_trials;
    }

    // Channel-parameter equality (A vs B), shared by the truncation
    // normaliser and the censoring branches below.
    const bool ch_eq_j = use_gl_pass ? ch_eq_vec[static_cast<size_t>(j)]
      : ((rule_code <= 4 || is_gng) &&
         row_equal_colmajor(pars_cols, n_par, idxA, idxB) &&
         row_equal_colmajor(pars_cols, n_par, idxnA, idxnB));

    // Sub-race win probability G_no at time tt for one channel. The GL
    // pre-pass batch-computes these at the aux slots (LT/UT/censor bounds);
    // slot values are used when present and matching tt, otherwise fall back
    // to the scalar GSL helper (non-raw path, or an uncovered combination).
    auto lr4_G_no_at = [&](double tt, int aux_slot, bool useB, bool& okf) -> double {
      okf = true;
      if (R_FINITE(tt) && use_gl_pass && aux_slot >= 0) {
        const size_t o = static_cast<size_t>(aux_slot) * n_unique_trials + j;
        const double g = (useB && !ch_eq_j) ? scratch.aux_GB[o] : scratch.aux_GA[o];
        if (!emc2_isnan(g) && !emc2_isnan(scratch.aux_t[o]) &&
            std::fabs(scratch.aux_t[o] - tt) < 1e-12) {
          return std::min(1.0, std::max(0.0, g));
        }
      }
      // tt = +Inf (GNG withheld / eventual N): P(nontarget eventually wins the
      // channel) via the QAGIU-capable helper; local_race_helper reports 0 for
      // an infinite horizon, so route infinite tt through lr_channel_yes_prob.
      const double g = R_FINITE(tt)
        ? local_race_helper(tt,
                            useB ? parnB.data() : parnA.data(),
                            useB ? parB.data()  : parA.data(),
                            n_par, model_ctx, pdf1, cdf1, gsl_ctl, w,
                            pars_2buf.data(), isok_2buf)
        : lr_channel_yes_prob(tt,
                            useB ? parnB.data() : parnA.data(),
                            useB ? parB.data()  : parA.data(),
                            n_par, model_ctx, pdf1, cdf1, gsl_ctl, w,
                            pars_2buf.data(), isok_2buf);
      okf = R_FINITE(g);
      return okf ? std::min(1.0, std::max(0.0, g)) : 0.0;
    };

    // Full channel state at tt (survivors are cheap scalar calls; G_no from
    // the batch/fallback above; G_yes by complement).
    auto lr4_state_at = [&](double tt, int aux_slot, bool useB, bool& okf) -> LrChannelState {
      LrChannelState s{0.0, 0.0, 1.0};
      okf = true;
      if (!(tt > 0.0)) return s;
      const double* pT = useB ? parB.data()  : parA.data();
      const double* pN = useB ? parnB.data() : parnA.data();
      const double S_dec = safe_surv_at_race_scalar(pT, tt, kMinSurv, cdf1, model_ctx)
                         * safe_surv_at_race_scalar(pN, tt, kMinSurv, cdf1, model_ctx);
      const double G_no = lr4_G_no_at(tt, aux_slot, useB, okf);
      if (!okf) return s;
      s.S_dec = S_dec;
      s.G_no  = G_no;
      s.G_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec - G_no));
      return s;
    };

    // N(t) = P(no overt rule response by t) for OR/AND/XOR/ID. XOR/ID need
    // only the channel survivors; OR/AND additionally need G_no.
    auto lr4_N_at = [&](double tt, int aux_slot, bool& okf) -> double {
      okf = true;
      if (!(tt > 0.0)) return 1.0;
      if (rule_code == 3 || rule_code == 4) {
        const double SdA =
          safe_surv_at_race_scalar(parA.data(),  tt, kMinSurv, cdf1, model_ctx) *
          safe_surv_at_race_scalar(parnA.data(), tt, kMinSurv, cdf1, model_ctx);
        const double SdB = ch_eq_j ? SdA
          : safe_surv_at_race_scalar(parB.data(),  tt, kMinSurv, cdf1, model_ctx) *
            safe_surv_at_race_scalar(parnB.data(), tt, kMinSurv, cdf1, model_ctx);
        return 1.0 - (1.0 - SdA) * (1.0 - SdB);
      }
      LrChannelState Ast = lr4_state_at(tt, aux_slot, false, okf);
      if (!okf) return 1.0;
      LrChannelState Bst = Ast;
      if (!ch_eq_j) {
        Bst = lr4_state_at(tt, aux_slot, true, okf);
        if (!okf) return 1.0;
      }
      if (is_gng) {
        // GNG: "no overt go response by t" = neither channel has said yes yet.
        // A completed channel-no is a withheld outcome here, not an overt
        // response, so (unlike OR) it stays inside N. tt may be +Inf (the
        // eventual withheld mass), which lr4_state_at handles via lr4_G_no_at.
        const double QA = std::min(1.0, std::max(0.0, Ast.G_no + Ast.S_dec));
        const double QB = std::min(1.0, std::max(0.0, Bst.G_no + Bst.S_dec));
        return QA * QB;
      }
      return lr_no_response_prob(rule_code, Ast, Bst);
    };

    // --- Truncation normaliser -------------------------------------------
    // Z = P(overt response RT in [LT, UT]) under the rule's outcome logic:
    // Z = N(LT) - N(UT), where the UT = Inf case is DEFINED as Z = N(LT).
    //
    // With a DEFECTIVE upper tail and a FINITE UT the eventual-no-response
    // mass N(Inf) is retained in the data -- make_missing() cuts only finite
    // RTs outside the window, and the rt = +Inf branch below correspondingly
    // does NOT subtract N(UT) from its numerator.  So N(Inf) has to be added
    // back here, or the retained withheld outcome would be scored against a
    // normaliser that excludes it.
    const double LTj_tr = shared.LT_unique[static_cast<size_t>(j)];
    const double UTj_tr = shared.UT_unique[static_cast<size_t>(j)];
    const bool has_trunc = (LTj_tr != 0.0 || R_FINITE(UTj_tr));
    double log_Z_j = 0.0;
    if (has_trunc) {
      if (!ISNAN(shared.cell_log_z[cell_idx])) {
        log_Z_j = shared.cell_log_z[cell_idx];
        if (log_Z_j == R_NegInf) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      } else {
        const int cond_j = shared.cond_code[static_cast<size_t>(j)];
        bool z_ok = true;
        double N_LT = 1.0;
        double N_UT = 0.0;
        if (rule_code == 5) {
          if (LTj_tr > 0.0) {
            N_LT = lr_detection_no_response_prob(rule_code, cond_j, LTj_tr,
                                                 parA.data(), parB.data(), parN.data(),
                                                 n_par, kMinSurv, pdf1, cdf1, model_ctx,
                                                 gsl_ctl, w, pars_3buf, z_ok);
          }
          if (z_ok && R_FINITE(UTj_tr)) {
            N_UT = lr_detection_no_response_prob(rule_code, cond_j, UTj_tr,
                                                  parA.data(), parB.data(), parN.data(),
                                                  n_par, kMinSurv, pdf1, cdf1, model_ctx,
                                                  gsl_ctl, w, pars_3buf, z_ok);
          }
        } else {
          if (LTj_tr > 0.0) N_LT = lr4_N_at(LTj_tr, LR_AUX_LT, z_ok);
          if (z_ok && R_FINITE(UTj_tr)) N_UT = lr4_N_at(UTj_tr, LR_AUX_UT, z_ok);
        }
        // Retained never-respond atom (zero unless the tail is defective).
        double N_Inf = 0.0;
        if (z_ok && R_FINITE(UTj_tr) && model_ctx->defective_upper_tail) {
          if (rule_code == 5) {
            // lr_detection_no_response_prob() has no infinite-horizon branch
            // (rule 6 needs a nogo-win quadrature over [t0_N, Inf)), so the
            // atom cannot be formed here.  Fail loudly rather than normalise
            // by a Z that is missing retained mass.
            Rcpp::stop("LogicalRules detection rules do not support a finite "
                       "UT with a defective upper tail: the retained "
                       "never-respond mass cannot be added to the truncation "
                       "normaliser. Remove UT or use a proper race model.");
          }
          N_Inf = lr4_N_at(R_PosInf, -1, z_ok);
        }
        const double Z = std::min(1.0, N_LT - N_UT + N_Inf);
        const bool valid_z = z_ok && R_FINITE(Z) && (Z > 1e-12);
        log_Z_j = valid_z ? std::log(Z) : R_NegInf;
        shared.cell_log_z[cell_idx] = log_Z_j;
        if (!valid_z) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
    }

    if (rule_code == 5) {
      ll_unique[static_cast<size_t>(j)] = logicalrules_detection_trial_ll(
        rule_code,
        shared.cond_code[static_cast<size_t>(j)],
        shared.resp_code[static_cast<size_t>(j)],
        t,
        shared.LT_unique[static_cast<size_t>(j)],
        shared.UT_unique[static_cast<size_t>(j)],
        shared.LC_unique[static_cast<size_t>(j)],
        shared.UC_unique[static_cast<size_t>(j)],
        idxA,
        idxB,
        shared.idxNogo[static_cast<size_t>(j)],
        n_par,
        min_ll,
        kMinSurv,
        log_Z_j,
        use_raw_local,
        logf_all,
        logS_all,
        parA.data(),
        parB.data(),
        parN.data(),
        pdf1,
        cdf1,
        model_ctx,
        gsl_ctl,
        w,
        pars_3buf);
      continue;
    }

    if (!(t > 0.0) || !R_FINITE(t)) {
      double p_j = 0.0;
      const int rule_code_nc = shared.rule_code[static_cast<size_t>(j)];
      if (!R_FINITE(t) && t > 0.0) {
        // Upper censor (+Inf): probability no overt response fires by UC.
        // Semantics depend on rule because each rule has a different "first
        // response" trigger (see lr_no_response_prob):
        //   XOR/ID : response fires when both channels resolve  → omit iff either channel undecided
        //   OR     : yes = first yes; no = second no           → omit = P(no yes by UC) − P(both no by UC)
        //   AND    : yes = second yes; no = first no           → omit = P(no no by UC) − P(both yes by UC)
        // With a finite upper truncation bound (and a proper upper tail) the
        // response is known to land in (UC, UT], so N(UT) is subtracted; the
        // trial is then normalised by Z below.
        // GNG: rt = +Inf is the withheld response — no overt go by UC, where
        // N_GNG(t) = Q_A(t) Q_B(t) (a completed channel-no is withheld, not an
        // overt response). UC may be +Inf (the withheld-to-infinity mass), so
        // GNG does not require a finite UC (unlike OR/AND/XOR/ID, whose rt=+Inf
        // omission is only defined against a finite censoring bound).
        const double UCj = shared.UC_unique[static_cast<size_t>(j)];
        bool ok_uc = true;
        if (is_gng) {
          p_j = lr4_N_at(UCj, R_FINITE(UCj) ? LR_AUX_UC : -1, ok_uc);
        } else {
          if (!R_FINITE(UCj)) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
          p_j = lr4_N_at(UCj, LR_AUX_UC, ok_uc);
        }
        if (ok_uc && R_FINITE(UTj_tr) && !model_ctx->defective_upper_tail) {
          bool ok_ut = true;
          const double N_UT_uc = lr4_N_at(UTj_tr, LR_AUX_UT, ok_ut);
          ok_uc = ok_uc && ok_ut;
          p_j = std::max(0.0, p_j - N_UT_uc);
        }
        if (!ok_uc) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
      } else if (!R_FINITE(t) && t < 0.0) {
        // Lower-censoring mass is the difference of response CDFs at the
        // window endpoints, rather than the first-event survivor difference.
        const double lo = std::max(0.0, shared.LT_unique[static_cast<size_t>(j)]);
        const double hi = std::max(lo,  shared.LC_unique[static_cast<size_t>(j)]);
        const int resp_lc = shared.resp_code[static_cast<size_t>(j)];
        if (is_gng) {
          // Lower censor: an overt go response fired in [lo, hi]. The go-CDF is
          // 1 - N_GNG, so the window mass is N_GNG(lo) - N_GNG(hi).
          bool ok_lc = true;
          const double N_lo = lr4_N_at(lo, LR_AUX_CLO, ok_lc);
          const double N_hi = ok_lc ? lr4_N_at(hi, LR_AUX_CHI, ok_lc) : 1.0;
          if (!ok_lc) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
          p_j = std::max(0.0, N_lo - N_hi);
          if (!(p_j > 0.0) || !R_FINITE(p_j)) {
            ll_unique[static_cast<size_t>(j)] = min_ll;
          } else {
            const double ll = std::log(p_j) - log_Z_j;
            ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
          }
          continue;
        }
        bool states_ok = true;
        LrChannelState A_lo = lr4_state_at(lo, LR_AUX_CLO, false, states_ok);
        LrChannelState A_hi, B_lo, B_hi;
        if (states_ok) A_hi = lr4_state_at(hi, LR_AUX_CHI, false, states_ok);
        if (states_ok) {
          if (ch_eq_j) {
            B_lo = A_lo;
            B_hi = A_hi;
          } else {
            B_lo = lr4_state_at(lo, LR_AUX_CLO, true, states_ok);
            if (states_ok) B_hi = lr4_state_at(hi, LR_AUX_CHI, true, states_ok);
          }
        }
        if (!states_ok) { ll_unique[static_cast<size_t>(j)] = min_ll; continue; }
        const double cdf_hi = lr_response_cdf(rule_code_nc, resp_lc, A_hi, B_hi);
        const double cdf_lo = lr_response_cdf(rule_code_nc, resp_lc, A_lo, B_lo);
        if (ISNAN(cdf_hi) || ISNAN(cdf_lo)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
        p_j = std::max(0.0, cdf_hi - cdf_lo);
      }
      if (!(p_j > 0.0) || !R_FINITE(p_j)) {
        ll_unique[static_cast<size_t>(j)] = min_ll;
      } else {
        const double ll = std::log(p_j) - log_Z_j;
        ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
      }
      continue;
    }

    double fA = 0.0, fB = 0.0, fnA = 0.0, fnB = 0.0;
    double FA = 0.0, FB = 0.0, FnA = 0.0, FnB = 0.0;
    if (!eval_pdf_cdf_race_scalar(use_raw_local, idxA, t, parA.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fA, FA) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxB, t, parB.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fB, FB) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxnA, t, parnA.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fnA, FnA) ||
        !eval_pdf_cdf_race_scalar(use_raw_local, idxnB, t, parnB.data(), min_ll, pdf1, cdf1, model_ctx, logf_all, logS_all, fnB, FnB)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    const double one_m_FA = std::max(kMinSurv, 1.0 - FA);
    const double one_m_FB = std::max(kMinSurv, 1.0 - FB);
    const double one_m_FnA = std::max(kMinSurv, 1.0 - FnA);
    const double one_m_FnB = std::max(kMinSurv, 1.0 - FnB);

    const double gA_yes = fA * one_m_FnA;
    const double gB_yes = fB * one_m_FnB;
    const double gA_no = fnA * one_m_FA;
    const double gB_no = fnB * one_m_FB;

    const double S_dec_A = one_m_FA * one_m_FnA;
    const double S_dec_B = one_m_FB * one_m_FnB;

    const int rule_code_lr = shared.rule_code[static_cast<size_t>(j)];
    const int resp_code = shared.resp_code[static_cast<size_t>(j)];
    // GNG shares OR's finite-RT go-response density; only its "no" outcome
    // differs (withheld, handled above), so a finite RT is always a "yes".
    const bool is_or_rule = (rule_code_lr == 1 || rule_code_lr == 6);
    const bool is_and_rule = (rule_code_lr == 2);
    const bool is_xor_rule = (rule_code_lr == 3);
    const bool is_id_rule = (rule_code_lr == 4);
    const bool channels_equal = use_gl_pass ? ch_eq_vec[static_cast<size_t>(j)]
        : (row_equal_colmajor(pars_cols, n_par, idxA, idxB) &&
           row_equal_colmajor(pars_cols, n_par, idxnA, idxnB));

    double p_j = 0.0;
    if (is_or_rule || is_and_rule) {
      double GA_no, GB_no;
      if (use_gl_pass) {
        GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        GB_no = channels_equal ? GA_no
                               : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
      } else {
        GA_no = local_race_helper(t, parnA.data(), parA.data(), n_par, model_ctx,
                                  pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_no = GA_no;
        if (!channels_equal)
          GB_no = local_race_helper(t, parnB.data(), parB.data(), n_par, model_ctx,
                                    pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_no) || !R_FINITE(GB_no)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double s_GA_yes = std::min(1.0, std::max(kMinSurv, GA_no + S_dec_A));
      const double s_GB_yes = std::min(1.0, std::max(kMinSurv, GB_no + S_dec_B));
      if (is_or_rule) {
        if (resp_code == 1) p_j = gA_yes * s_GB_yes + gB_yes * s_GA_yes; // yes
        else if (resp_code == 2 && rule_code_lr != 6) p_j = gA_no * GB_no + gB_no * GA_no; // no
      } else { // AND
        const double GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        const double GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
        const double s_GA_no = std::min(1.0, std::max(kMinSurv, GA_yes + S_dec_A));
        const double s_GB_no = std::min(1.0, std::max(kMinSurv, GB_yes + S_dec_B));
        if (resp_code == 1) p_j = gA_yes * GB_yes + gB_yes * GA_yes;      // yes
        else if (resp_code == 2) p_j = gA_no * s_GB_no + gB_no * s_GA_no; // no
      }
    } else if (is_xor_rule || is_id_rule) {
      double GA_yes, GB_yes;
      if (use_gl_pass) {
        // Derive GA_yes from GA_no_gl: GA_yes + GA_no + S_dec_A = 1.
        const double GA_no = std::max(0.0, std::min(GA_no_gl[static_cast<size_t>(j)], 1.0));
        const double GB_no = channels_equal ? GA_no
                                            : std::max(0.0, std::min(GB_no_gl[static_cast<size_t>(j)], 1.0));
        GA_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_no));
        GB_yes = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_no));
      } else {
        GA_yes = local_race_helper(t, parA.data(), parnA.data(), n_par, model_ctx,
                                   pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        GB_yes = GA_yes;
        if (!channels_equal)
          GB_yes = local_race_helper(t, parB.data(), parnB.data(), n_par, model_ctx,
                                     pdf1, cdf1, gsl_ctl, w, pars_2buf.data(), isok_2buf);
        if (!R_FINITE(GA_yes) || !R_FINITE(GB_yes)) {
          ll_unique[static_cast<size_t>(j)] = min_ll;
          continue;
        }
      }
      const double GA_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_A - GA_yes));
      const double GB_no = std::min(1.0, std::max(0.0, 1.0 - S_dec_B - GB_yes));
      const double dAB = gA_yes * GB_yes + gB_yes * GA_yes;
      const double dAN = gA_yes * GB_no + gB_no * GA_yes;
      const double dNB = gA_no * GB_yes + gB_yes * GA_no;
      const double dNN = gA_no * GB_no + gB_no * GA_no;

      if (is_id_rule) {
        if (resp_code == 6) p_j = dAB;
        else if (resp_code == 4) p_j = dAN;
        else if (resp_code == 5) p_j = dNB;
        else if (resp_code == 3) p_j = dNN;
      } else { // XOR
        if (resp_code == 1) p_j = dAN + dNB; // yes
        else if (resp_code == 2) p_j = dAB + dNN; // no
      }
    } else {
      ll_unique[static_cast<size_t>(j)] = min_ll;
      continue;
    }

    if (!(p_j > 0.0) || !R_FINITE(p_j)) {
      ll_unique[static_cast<size_t>(j)] = min_ll;
    } else {
      const double ll = std::log(p_j) - log_Z_j;
      ll_unique[static_cast<size_t>(j)] = (R_FINITE(ll) && ll > min_ll) ? ll : min_ll;
    }
  }

  // Contaminant mixture: pC for intrinsic omissions (+Inf RT), pGuess for the
  // uniform outlier on observed RTs.  Both read the first row of each trial.
  // See src/contaminant_mixture.h -- the arithmetic lives there so this and the
  // seven other application sites cannot drift apart.
  bool use_pC = (pc_col >= 0);
  if (use_pC) {
    bool all_zero = true;
    for (int j = 0; j < n_unique_trials; ++j) {
      if (pars_cols[pc_col][j * n_acc] != 0.0) { all_zero = false; break; }
    }
    if (all_zero) use_pC = false;
  }

  const GuessKernel gk = (guess != nullptr) ? *guess : GuessKernel();
  bool use_pG = (pg_col >= 0) && gk.active();
  if (use_pG) {
    bool all_zero = true;
    for (int j = 0; j < n_unique_trials; ++j) {
      if (pars_cols[pg_col][j * n_acc] != 0.0) { all_zero = false; break; }
    }
    if (all_zero) use_pG = false;
  }

  if (use_pC || use_pG) {
    for (int j = 0; j < n_unique_trials; ++j) {
      const size_t sj = static_cast<size_t>(j);
      const double pC = use_pC ? pars_cols[pc_col][sj * n_acc] : 0.0;
      const double pG = use_pG ? pars_cols[pg_col][sj * n_acc] : 0.0;
      // resp_code 0 is a missing response; every other code names a response.
      ll_unique[sj] = mix_contaminants_rt(ll_unique[sj], pC, pG, gk,
                                          shared.rt_unique[sj],
                                          shared.resp_code[sj] != 0);
    }
  }

  double sum_ll = 0.0;
  if (expand.length() > 0) {
    const int n_out = expand.length();
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (expand path).");
    }
    for (int i = 0; i < n_out; ++i) {
      double val = ll_unique[static_cast<size_t>(expand[i] - 1)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[i] = val;
      sum_ll += val;
    }
  } else {
    const int n_out = n_unique_trials;
    if (trial_ll_out != nullptr && trial_ll_out->size() != n_out) {
      Rcpp::stop("c_log_likelihood_logicalrules: trial_ll_out size mismatch (compressed path).");
    }
    for (int j = 0; j < n_unique_trials; ++j) {
      double val = ll_unique[static_cast<size_t>(j)];
      if (!R_FINITE(val) || val < min_ll) val = min_ll;
      if (trial_ll_out != nullptr) (*trial_ll_out)[j] = val;
      sum_ll += val;
    }
  }
  return sum_ll;
}

// Test/benchmark observability accessors for the LogicalRules correlated
// capacity route.  Counting is active only while EMC2_LRCAP_COUNTERS is set;
// the accessors themselves always return the current counters.
// [[Rcpp::export]]
Rcpp::List lr_capacity_counter_values() {
  const LrCapacityCounters& c = lr_capacity_counters();
  return Rcpp::List::create(
      Rcpp::Named("ordinary_trials") = (double)c.ordinary_trials,
      Rcpp::Named("capacity_detection_trials") = (double)c.capacity_detection_trials,
      Rcpp::Named("capacity_choice_trials") = (double)c.capacity_choice_trials,
      Rcpp::Named("invalid_trials") = (double)c.invalid_trials,
      Rcpp::Named("tau_zero_trials") = (double)c.tau_zero_trials,
      Rcpp::Named("factor_node_evaluations") = (double)c.factor_node_evaluations,
      Rcpp::Named("channel_gl_evaluations") = (double)c.channel_gl_evaluations);
}

// [[Rcpp::export]]
void lr_capacity_counters_reset() {
  lr_capacity_counters().reset();
}
