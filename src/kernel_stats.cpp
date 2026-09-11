#include "kernel_stats.h"

#include <Rcpp.h>

#include <algorithm>
#include <map>

namespace emc {

namespace {

bool& on_storage() {
  static bool on = false;
  return on;
}

std::map<std::string, KernelStat>& table() {
  static std::map<std::string, KernelStat> t;
  return t;
}

}  // namespace

bool kernel_stats_on() { return on_storage(); }

void set_kernel_stats(bool on) { on_storage() = on; }

void kernel_stats_record(const std::string& model, long long rows,
                         long long cells, double seconds) {
  KernelStat& s = table()[model];
  s.model = model;
  s.calls += 1;
  s.rows += rows;
  s.cells += (cells < 0) ? rows : cells;
  s.seconds += seconds;
}

std::vector<KernelStat> kernel_stats_read() {
  std::vector<KernelStat> out;
  out.reserve(table().size());
  for (std::map<std::string, KernelStat>::const_iterator it = table().begin();
       it != table().end(); ++it) {
    out.push_back(it->second);
  }
  return out;
}

void kernel_stats_reset() { table().clear(); }

std::vector<std::string> kernel_reuse_columns(const std::string& c_name) {
  std::vector<std::string> cols;
  // RDM and its variants: `inv_s = 1/s`, and the geometry `(B + A/2)/s`,
  // `v/s`, `A/(2s)` handed to the Wald primitives.  `t0` is deliberately not
  // here: it enters only through `rt - t0`, which is per trial whatever the
  // design says.
  if (c_name.compare(0, 3, "RDM") == 0) {
    cols.push_back("v");
    cols.push_back("B");
    cols.push_back("A");
    cols.push_back("s");
    return cols;
  }
  // LBA and the ballistic accumulators built on it: `natural_normalizer(v, sv)`
  // is a pnorm of v/sv and nothing else.
  if (c_name.compare(0, 3, "LBA") == 0 || c_name.compare(0, 4, "BAwL") == 0 ||
      c_name.compare(0, 4, "BAwD") == 0 || c_name.compare(0, 4, "BAwF") == 0 ||
      c_name.compare(0, 4, "BAwR") == 0) {
    cols.push_back("v");
    cols.push_back("sv");
    return cols;
  }
  // DDM: the scale divisions and `log(a)`.
  if (c_name.compare(0, 3, "DDM") == 0) {
    cols.push_back("a");
    cols.push_back("v");
    cols.push_back("sv");
    cols.push_back("s");
    return cols;
  }
  return cols;
}

}  // namespace emc

//' Kernel reuse instrumentation
//'
//' Reads the flag; with a value, sets it and returns the previous one.
//'
//' @noRd
// [[Rcpp::export]]
bool emc_kernel_stats(Rcpp::Nullable<Rcpp::LogicalVector> on = R_NilValue) {
  const bool previous = emc::kernel_stats_on();
  if (on.isNotNull()) {
    Rcpp::LogicalVector v(on);
    if (v.size() != 1 || Rcpp::LogicalVector::is_na(v[0])) {
      Rcpp::stop("emc_kernel_stats() takes one TRUE or FALSE");
    }
    emc::set_kernel_stats(v[0] != 0);
  }
  return previous;
}

//' Read the kernel reuse counters
//'
//' One row per model seen since the last reset.  `rows` is how many per-trial
//' evaluations the kernel performed and `cells` how many a cell-resolution
//' version would have performed, so `rows / cells` is the reuse available.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame emc_kernel_stats_read() {
  const std::vector<emc::KernelStat> s = emc::kernel_stats_read();
  const int n = static_cast<int>(s.size());
  Rcpp::CharacterVector model(n);
  Rcpp::NumericVector calls(n), rows(n), cells(n), seconds(n);
  for (int i = 0; i < n; ++i) {
    model[i] = s[i].model;
    calls[i] = static_cast<double>(s[i].calls);
    rows[i] = static_cast<double>(s[i].rows);
    cells[i] = static_cast<double>(s[i].cells);
    seconds[i] = s[i].seconds;
  }
  return Rcpp::DataFrame::create(
    Rcpp::_["model"] = model, Rcpp::_["calls"] = calls,
    Rcpp::_["rows"] = rows, Rcpp::_["cells"] = cells,
    Rcpp::_["seconds"] = seconds, Rcpp::_["stringsAsFactors"] = false);
}

//' @noRd
// [[Rcpp::export]]
void emc_kernel_stats_reset() { emc::kernel_stats_reset(); }

//' The columns a model's reusable subexpression reads
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::CharacterVector emc_kernel_reuse_columns(std::string c_name) {
  return Rcpp::wrap(emc::kernel_reuse_columns(c_name));
}
