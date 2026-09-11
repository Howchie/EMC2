#include "emc_scratch.h"

#include <Rcpp.h>

namespace emc {

namespace {
// Function-local static, so there is exactly one budget however many
// translation units include the header.
std::size_t& budget_storage() {
  // 8 MiB admits a few hundred thousand cells at the mapper's three doubles
  // each, which is far past any design anyone fits and still refuses a
  // pathological one rather than trying to allocate for it.
  static std::size_t bytes = 8u * 1024u * 1024u;
  return bytes;
}
}  // namespace

std::size_t cell_scratch_budget() { return budget_storage(); }

void set_cell_scratch_budget(std::size_t bytes) { budget_storage() = bytes; }

namespace {
bool& defer_storage() {
  static bool on = true;
  return on;
}
}  // namespace

bool defer_scalar_fill() { return defer_storage(); }

void set_defer_scalar_fill(bool on) { defer_storage() = on; }

}  // namespace emc

//' Per-design cell scratch budget, in bytes
//'
//' Reads the budget; with a value, sets it and returns the previous one. The
//' setter exists so a test can drive a design onto the general row route
//' without building one with hundreds of thousands of cells.
//'
//' @noRd
// [[Rcpp::export]]
double emc_pt_cell_budget(Rcpp::Nullable<Rcpp::NumericVector> bytes = R_NilValue) {
  const double previous = static_cast<double>(emc::cell_scratch_budget());
  if (bytes.isNotNull()) {
    Rcpp::NumericVector v(bytes);
    if (v.size() != 1 || !R_finite(v[0]) || v[0] < 0) {
      Rcpp::stop("the cell scratch budget must be one finite, non-negative number of bytes");
    }
    emc::set_cell_scratch_budget(static_cast<std::size_t>(v[0]));
  }
  return previous;
}

//' Defer the per-trial fill of scalar coefficients
//'
//' Reads the setting; with a value, sets it and returns the previous one.
//'
//' @noRd
// [[Rcpp::export]]
bool emc_pt_defer_scalars(Rcpp::Nullable<Rcpp::LogicalVector> on = R_NilValue) {
  const bool previous = emc::defer_scalar_fill();
  if (on.isNotNull()) {
    Rcpp::LogicalVector v(on);
    if (v.size() != 1 || Rcpp::LogicalVector::is_na(v[0])) {
      Rcpp::stop("emc_pt_defer_scalars() takes one TRUE or FALSE");
    }
    emc::set_defer_scalar_fill(v[0] != 0);
  }
  return previous;
}
