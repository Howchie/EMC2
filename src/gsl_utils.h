#ifndef GSL_UTILS_H
#define GSL_UTILS_H

#include <gsl/gsl_integration.h>
#include <memory>
#include <Rcpp.h>

// Unique pointer for GSL integration workspace with custom deleter
using GslWorkspacePtr = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;

// Ensure a GSL workspace is allocated in the provided pointer.
// Must GROW when a caller needs more subdivisions than the cached workspace
// holds: qags fails outright (GSL_EINVAL -> caller's min_ll path) when its
// limit argument exceeds ws->limit, so a small first allocation must not pin
// the size for the rest of the thread's lifetime.
inline gsl_integration_workspace* ensure_gsl_workspace(GslWorkspacePtr& ws, size_t n = 1000) {
  if (!ws || ws->limit < n) {
    ws.reset(gsl_integration_workspace_alloc(n));
    if (!ws) Rcpp::stop("Failed to allocate GSL integration workspace.");
  }
  return ws.get();
}

// Controls for GSL integration, including retry logic
struct GslIntegrationControls {
  double abs_tol;
  double rel_tol;
  size_t limit;
  double retry_abs_tol;
  double retry_rel_tol;
  size_t retry_limit;
};

// Default controls matching the existing logic
inline GslIntegrationControls default_gsl_controls() {
  return {
    1e-10, // abs_tol
    1e-5,  // rel_tol
    200,   // max subintervals
    0.0,   // retry_abs_tol
    1e-7,  // retry_rel_tol
    1000   // retry max subintervals
  };
}

#endif // GSL_UTILS_H
