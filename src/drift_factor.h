#ifndef EMC2_DRIFT_FACTOR_H
#define EMC2_DRIFT_FACTOR_H

// Shared equicorrelated drift-factor construction for the correlated race
// models (BAwLcorr, RDMSWTNcorr with correlate = "drifts").
//
// Every one of those models draws its accumulator rates from
//
//     V_i ~ N(v_i, sv_i^2),   cor(V_i, V_j) = rho
//
// optionally conditioned on the positive orthant, and races independently
// given the drawn rates.  A single standard-normal latent factor z reproduces
// that equicorrelation exactly:
//
//     V_i | z ~ N(v_i + sign(rho) sv_i sqrt(|rho|) z,
//                 (sv_i sqrt(1 - |rho|))^2)
//
// so conditional on z the accumulators are independent and each one is just
// the model's ordinary scalar kernel evaluated at the shifted mean and shrunk
// SD.  Nothing here is specific to a particular accumulator: the caller
// supplies the column positions of v and sv, and the model's own likelihood
// does the rest.
//
// The sign trick makes a negative rho a genuine two-accumulator
// anti-correlation; with three or more loaded rows only rho >= 0 is a valid
// equicorrelation matrix, which is the callers' responsibility to enforce.
//
// This header is deliberately free of Rcpp objects apart from the matrix
// helper at the bottom, so the hot paths can keep using raw column pointers.

#include <cmath>
#include "wald_functions.h"   // pnorm_log_direct

// Residual SDs are floored strictly above zero so a hand-built parameter
// matrix sitting exactly on |rho| = 1 cannot hand a scalar kernel an
// undefined sv of zero.  The public rho bound stays just inside +/-1.
constexpr double DRIFT_FACTOR_SD_FLOOR = 1e-12;

// Below this magnitude a row is treated as carrying no loading at all: its
// conditional mean and SD are its marginal ones, and its positivity constant
// cancels between the numerator and denominator integrands.
constexpr double DRIFT_FACTOR_RHO_EPS = 1e-14;

struct DriftFactorColumns {
  int v = -1;
  int sv = -1;
  bool valid() const { return v >= 0 && sv >= 0; }
};

struct DriftFactorLoading {
  double slope = 0.0;   // sign(rho) * sv * sqrt(|rho|); multiplies z
  double sv_res = 0.0;  // sv * sqrt(1 - |rho|), floored
};

inline DriftFactorLoading drift_factor_loading(double sv, double rho) {
  const double magnitude = std::fabs(rho);
  DriftFactorLoading out;
  out.slope = ((rho < 0.0) ? -1.0 : 1.0) * sv * std::sqrt(magnitude);
  out.sv_res = sv * std::fmax(std::sqrt(std::fmax(0.0, 1.0 - magnitude)),
                              DRIFT_FACTOR_SD_FLOOR);
  return out;
}

inline double drift_factor_mean(double v, const DriftFactorLoading& l,
                                double z) {
  return v + l.slope * z;
}

// log P(V > 0 | z) for one row, on the node-shifted mean and residual SD.
// The product of these over a trial's loaded rows is the unnormalised
// positivity weight the posdrift denominator integrates.
inline double drift_factor_log_positive_row(double mu, double sd) {
  if (!(sd > 0.0) || !R_FINITE(mu)) return R_NegInf;
  return pnorm_log_direct(mu / sd, true);
}

// Rewrite one row of a materialised parameter matrix in place to its
// conditional-on-z values.  Returns false when rho is unusable, which the
// caller turns into an inactive row.
inline bool drift_factor_apply_row(Rcpp::NumericMatrix& pars,
                                   const DriftFactorColumns& cols,
                                   int row, double rho, double z) {
  if (!R_FINITE(rho) || std::fabs(rho) > 1.0) return false;
  const double sv = pars(row, cols.sv);
  const DriftFactorLoading l = drift_factor_loading(sv, rho);
  pars(row, cols.v) = drift_factor_mean(pars(row, cols.v), l, z);
  pars(row, cols.sv) = l.sv_res;
  return true;
}

// Draw one equicorrelated rate per accumulator row into `drifts`, rejecting
// jointly non-positive draws under `posdrift`.  Rows with ok = false keep
// whatever `drifts` already held, so a caller can pre-seed them.  This is the
// simulation counterpart of the conditional construction above: one shared z
// per trial, then an independent residual per row.
inline void drift_factor_draw_correlated(const Rcpp::NumericMatrix& pars,
                                         int iv, int isv, int irho,
                                         const Rcpp::LogicalVector& ok,
                                         int n_acc, bool posdrift,
                                         std::vector<double>& drifts,
                                         const char* caller) {
  const int n_rows = pars.nrow();
  const int n_trials = n_rows / n_acc;
  // A joint positive-orthant draw is rejection sampling, so drift means far
  // below zero can make acceptance arbitrarily rare; fail loudly rather than
  // spin.
  const int max_iter = 100000;
  for (int tr = 0; tr < n_trials; ++tr) {
    const int start = tr * n_acc;
    bool has_active = false;
    for (int a = 0; a < n_acc; ++a) has_active = has_active || static_cast<bool>(ok[start + a]);
    if (!has_active) continue;

    bool accepted = false;
    for (int iter = 0; iter < max_iter; ++iter) {
      const double z = R::norm_rand();
      bool positive = true;
      for (int a = 0; a < n_acc; ++a) {
        const int r = start + a;
        if (!ok[r]) continue;
        const double rho = pars(r, irho);
        if (!R_FINITE(rho) || std::fabs(rho) > 1.0) {
          Rcpp::stop("%s: rho must be finite and lie in [-1, 1].", caller);
        }
        const DriftFactorLoading l = drift_factor_loading(pars(r, isv), rho);
        const double draw =
          R::rnorm(drift_factor_mean(pars(r, iv), l, z), l.sv_res);
        drifts[static_cast<size_t>(r)] = draw;
        if (posdrift && !(draw > 0.0)) positive = false;
      }
      if (!posdrift || positive) {
        accepted = true;
        break;
      }
    }
    if (!accepted) {
      Rcpp::stop("%s: jointly positive drift rejection exceeded %d attempts; "
                 "check that the drift means are not far below zero.",
                 caller, max_iter);
    }
  }
}

#endif
