#include "emc_scratch.h"

#include <Rcpp.h>
#include <limits>

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
    const double max_size = static_cast<double>(
      std::numeric_limits<std::size_t>::max());
    const bool max_is_exact =
      std::numeric_limits<std::size_t>::digits <=
      std::numeric_limits<double>::digits;
    if (v.size() != 1 || !R_finite(v[0]) || v[0] < 0 ||
        (max_is_exact ? v[0] > max_size : v[0] >= max_size)) {
      Rcpp::stop("the cell scratch budget must be one finite, non-negative number of bytes representable by size_t");
    }
    emc::set_cell_scratch_budget(static_cast<std::size_t>(v[0]));
  }
  return previous;
}
