#ifndef CONTAMINANT_MIXTURE_H
#define CONTAMINANT_MIXTURE_H

// Shared contaminant mixture arithmetic.
//
// EMC2 carries two independent contaminant processes, combined as nested
// (stick-breaking) weights so that neither can push the other out of [0, 1]:
//
//   P(omission) = pC
//   P(guess)    = (1 - pC) * pG
//   P(process)  = (1 - pC) * (1 - pG)
//
// `pContaminant` (pC) is an *omission* mixture: it contributes mass only at
// rt == +Inf and does nothing for observed RTs.  `pGuess` (pG) is the standard
// uniform-outlier mixture (Ratcliff & Tuerlinckx 2002; HDDM's `w_outlier`):
// it contributes a flat density directly to observed RT densities, giving
// fast/slow outliers a likelihood floor.
//
// At pG == 0 every expression below reduces *exactly* to the pre-pGuess
// arithmetic, which is what protects `pContaminant`'s established behaviour.
//
// Eight call sites across particle_ll.cpp apply this mixture (the generic race
// kernel, the race all-finite fast path, both LogicalRules paths, BAwLcorr,
// RDMSWTNcorr, and both DDM paths).  They live here rather than being copied
// so they cannot drift apart again -- the file they came from already carried a
// comment about exactly that failure mode.

#include <Rcpp.h>
#include <cmath>
#include "composite_functions.h" // log1m, log_sum_exp

// The uniform guess kernel, resolved once per dadm rather than per trial.
//
// `resolve_guess_window()` (R/design.R) computes the window from the existing
// truncation/censoring columns and attaches it to the dadm; because the window
// is [max(LT, LC), min(UC, UT)] by construction, a guess can never be censored
// or truncated away.  That is what removes all interval-mass machinery from
// this file: there is no guess term in any censored branch.
struct GuessKernel {
  int    pg_col     = -1;         // pGuess column in the ParamTable, -1 when absent
  double log_g      = R_NegInf;   // -log(n_resp * (UG - LG)): guess density for a known R
  double log_g_any  = R_NegInf;   // -log(UG - LG): R marginalised out (n_resp skipped)
  double LG         = 0.0;
  double UG         = R_PosInf;

  bool active() const { return pg_col >= 0 && log_g > R_NegInf; }
};

// Read attr(dadm, "guess_window") / attr(dadm, "guess_n_resp") and turn them
// into the two log densities.  `pg_col` is filled in separately by whichever
// column walk the caller already performs for pContaminant.
inline GuessKernel resolve_guess_kernel(const Rcpp::DataFrame& dadm) {
  GuessKernel g;
  Rcpp::RObject win_obj = dadm.attr("guess_window");
  Rcpp::RObject nr_obj  = dadm.attr("guess_n_resp");
  if (win_obj.isNULL() || nr_obj.isNULL()) return g;
  Rcpp::NumericVector win(win_obj);
  if (win.size() != 2) return g;
  const double LG = win[0], UG = win[1];
  const double n_resp = Rcpp::as<double>(nr_obj);
  if (!R_FINITE(LG) || !R_FINITE(UG) || !(UG > LG) || !(n_resp >= 1.0)) return g;
  g.LG        = LG;
  g.UG        = UG;
  g.log_g_any = -std::log(UG - LG);
  g.log_g     = g.log_g_any - std::log(n_resp);
  return g;
}

// The per-trial mixture.  Only two cases exist:
//
//   finite rt, R known   log[(1-pC) * ((1-pG) * L_proc + pG * exp(log_g))]
//   everything else      log[pC * 1{rt == +Inf} + (1-pC) * (1-pG) * L_proc]
//
// The second row covers +Inf, -Inf, NA, withheld and nogo trials: a guess is an
// overt response, so it can be neither a timeout nor a withheld trial.  Callers
// signal that by passing log_g = R_NegInf.
//
// The mixture is proper: pC + (1-pC) * [pG * 1 + (1-pG) * 1] = 1.
//
// `ll_proc` must already be truncation-renormalised -- both components are
// renormalised on the truncation window before mixing, which makes pG the guess
// proportion among *retained* trials (matching how make_missing() contaminates
// after the truncation cut).
inline double mix_contaminants(double ll_proc, double pC, double pG,
                               double log_g, bool is_omission) {
  double v = ll_proc;
  if (pG > 0.0 && pG < 1.0) {
    // log_g == -Inf (not guess-eligible) still down-weights by (1-pG); the
    // guess mass for such a trial is carried by the eligible trials instead.
    v = (log_g == R_NegInf) ? (log1m(pG) + v)
                            : log_sum_exp(log1m(pG) + v, std::log(pG) + log_g);
  } else if (pG >= 1.0) {
    v = (log_g == R_NegInf) ? R_NegInf : log_g;
  }
  if (pC != 0.0) {
    const double log1m_pC = log1m(pC);
    v = is_omission ? log_sum_exp(std::log(pC), log1m_pC + v)
                    : log1m_pC + v;
  }
  return v;
}

// Convenience wrapper for the common case where the caller has rt and knows
// whether R is observed.  Picks log_g / log_g_any / -Inf per decision 6:
// an unknown R with a finite rt takes the window density *without* the n_resp
// division, so the density still integrates to 1 over that trial's support.
inline double mix_contaminants_rt(double ll_proc, double pC, double pG,
                                  const GuessKernel& g, double rt, bool R_known) {
  double log_g = R_NegInf;
  if (g.active() && R_FINITE(rt) && rt > 0.0) log_g = R_known ? g.log_g : g.log_g_any;
  return mix_contaminants(ll_proc, pC, pG, log_g, rt == R_PosInf);
}

#endif // CONTAMINANT_MIXTURE_H
