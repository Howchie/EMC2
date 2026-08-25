#ifndef timer_helpers_h
#define timer_helpers_h

#include <cmath>

#include "utility_functions.h"
#include "race_contract.h"

inline double erlang_omega_for_shape(int kill_shape, const double* par = nullptr,
                                     int omega_index = -1) {
  if (kill_shape <= 1) return 1.0;
  if (kill_shape == 2) return 0.0;
  if (par == nullptr || omega_index < 0) return 1.0;
  return std::fmax(0.0, std::fmin(1.0, par[omega_index]));
}

// Convert a timer mean to the Erlang rate parameter used by erlang_log_surv /
// erlang_log_pdf.  The C++ particle likelihood receives raw timer means (mG,
// mK) directly from the sampled parameter space; Ttransform is NOT applied on
// the C++ path.  This function performs the shape-dependent conversion:
//   Erlang-1 (exponential): rate = 1 / mean
//   Erlang-2:               rate = 2 / mean   (so that E[T] = 2/rate = mean)
//   EMIX (shape 3):         rate = 1 / mean   (each component is rescaled
//                                               inside erlang_log_surv for n=3)
inline double erlang_lambda_from_mean(double mean, int kill_shape) {
  if (!(mean > 0.0) || !emc2_isfinite(mean)) return 0.0;
  return ((kill_shape == 2) ? 2.0 : 1.0) / mean;
}

struct TimedLambdaDispatch {
  double lambda_g;
  double lambda_k;
  bool guess;
  bool use_combo;
};

inline TimedLambdaDispatch timed_lambda_dispatch(const ContextForRaceModels* ctx,
                                                 double lambda_g,
                                                 double lambda_k) {
  constexpr double kLamEps = 1e-12;
  const bool local_guess_only = ctx && ctx->is_local_guess;
  const bool local_kill_guess = ctx && ctx->is_local_kill_guess;
  const bool has_guess = lambda_g > kLamEps;
  const bool has_kill  = lambda_k > kLamEps;

  TimedLambdaDispatch out{0.0, 0.0, false, local_kill_guess && has_guess && has_kill};
  if (out.use_combo) return out;

  if (local_guess_only || (local_kill_guess && has_guess)) {
    out.lambda_g = lambda_g;
    out.guess = true;
  } else if (has_kill) {
    out.lambda_k = lambda_k;
  }
  return out;
}

#endif
