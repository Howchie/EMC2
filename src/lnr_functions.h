#ifndef lnr_functions_h
#define lnr_functions_h

#include <cmath>
#include <Rcpp.h>
#include "wald_functions.h"

using namespace Rcpp;

// LNR-specific lognormal helpers.
// Keep these separate from wald_functions.h so the LNR dependency is explicit.
inline double plnorm_std(double x, double meanlog, double sdlog,
                         bool lower_tail = true, bool log_p = false) {
  if (x <= 0.0) {
    if (log_p) return lower_tail ? R_NegInf : 0.0;
    return lower_tail ? 0.0 : 1.0;
  }
  if (!std::isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::plnorm(x, meanlog, sdlog, lower_tail, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return pnorm_std(z, lower_tail, log_p);
}

inline double dlnorm_std(double x, double meanlog, double sdlog,
                         bool log_p = false) {
  if (x <= 0.0) return log_p ? R_NegInf : 0.0;
  if (!std::isfinite(sdlog) || !(sdlog > 0.0)) {
    return R::dlnorm(x, meanlog, sdlog, log_p);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  const double log_pdf = -std::log(x) - std::log(sdlog) - 0.5 * z * z - LOG_SQRT_2PI;
  return log_p ? log_pdf : std::exp(log_pdf);
}

inline double lnorm_log_surv_std(double x, double meanlog, double sdlog) {
  if (x <= 0.0) return 0.0;
  if (!std::isfinite(sdlog) || !(sdlog > 0.0)) {
    const double cdf = R::plnorm(x, meanlog, sdlog, true, false);
    if (cdf >= 1.0) return R_NegInf;
    if (cdf <= 0.0) return 0.0;
    return std::log1p(-cdf);
  }
  const double z = (std::log(x) - meanlog) / sdlog;
  return pnorm_std(z, false, true);
}

#endif
