#include "emc_scratch.h"

#include <Rcpp.h>
#include <atomic>
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

double& reuse_storage() {
  // Zero retains the pre-Stage-5 admission rule.
  static double ratio = 0.0;
  return ratio;
}

bool& defer_storage() {
  static bool on = true;
  return on;
}

struct CounterStorage {
  std::atomic<bool> enabled{false};
  std::atomic<std::uint64_t> values[13];
  CounterStorage() {
    for (std::size_t i = 0; i < 13; ++i) values[i].store(0);
  }
};

CounterStorage& counters() {
  static CounterStorage s;
  return s;
}

inline void count(std::size_t i) {
  CounterStorage& s = counters();
  if (s.enabled.load(std::memory_order_relaxed))
    s.values[i].fetch_add(1, std::memory_order_relaxed);
}
}  // namespace

std::size_t cell_scratch_budget() { return budget_storage(); }

void set_cell_scratch_budget(std::size_t bytes) { budget_storage() = bytes; }

double cell_min_reuse_ratio() { return reuse_storage(); }

void set_cell_min_reuse_ratio(double ratio) { reuse_storage() = ratio; }

bool defer_scalar_fill() { return defer_storage(); }

void set_defer_scalar_fill(bool on) { defer_storage() = on; }

bool mapper_stats_enabled() {
  return counters().enabled.load(std::memory_order_relaxed);
}

void set_mapper_stats_enabled(bool on) {
  counters().enabled.store(on, std::memory_order_relaxed);
}

void reset_mapper_stats() {
  CounterStorage& s = counters();
  for (std::size_t i = 0; i < 13; ++i)
    s.values[i].store(0, std::memory_order_relaxed);
}

MapperStatsSnapshot mapper_stats_snapshot() {
  CounterStorage& s = counters();
  MapperStatsSnapshot out;
  out.map_scalar = s.values[0].load(std::memory_order_relaxed);
  out.map_cell = s.values[1].load(std::memory_order_relaxed);
  out.map_row = s.values[2].load(std::memory_order_relaxed);
  out.transform_scalar = s.values[3].load(std::memory_order_relaxed);
  out.transform_cell = s.values[4].load(std::memory_order_relaxed);
  out.transform_row = s.values[5].load(std::memory_order_relaxed);
  out.bound_scalar = s.values[6].load(std::memory_order_relaxed);
  out.bound_cell = s.values[7].load(std::memory_order_relaxed);
  out.bound_row = s.values[8].load(std::memory_order_relaxed);
  out.cell_admit = s.values[9].load(std::memory_order_relaxed);
  out.cell_reject = s.values[10].load(std::memory_order_relaxed);
  out.cell_budget_reject = s.values[11].load(std::memory_order_relaxed);
  out.cell_reuse_reject = s.values[12].load(std::memory_order_relaxed);
  return out;
}

void mapper_count_map(MapperRoute route) {
  if (route <= MAPPER_ROW) count(static_cast<std::size_t>(route));
}

void mapper_count_transform(MapperRoute route) {
  if (route <= MAPPER_ROW) count(3u + static_cast<std::size_t>(route));
}

void mapper_count_bound(MapperRoute route) {
  if (route <= MAPPER_ROW) count(6u + static_cast<std::size_t>(route));
}

void mapper_count_cell_admit() { count(9); }

void mapper_count_cell_reject(bool budget, bool reuse) {
  CounterStorage& s = counters();
  if (!s.enabled.load(std::memory_order_relaxed)) return;
  s.values[10].fetch_add(1, std::memory_order_relaxed);
  if (budget) s.values[11].fetch_add(1, std::memory_order_relaxed);
  if (reuse) s.values[12].fetch_add(1, std::memory_order_relaxed);
}

}  // namespace emc
//' Minimum fraction of trials reused by an admitted cell design
//'
//' Zero retains the historical admission rule.  A value in [0, 1] requires
//' at least that fraction of trial rows to share design cells.
//'
//' @noRd
// [[Rcpp::export]]
double emc_pt_cell_min_reuse(Rcpp::Nullable<Rcpp::NumericVector> ratio = R_NilValue) {
  const double previous = emc::cell_min_reuse_ratio();
  if (ratio.isNotNull()) {
    Rcpp::NumericVector v(ratio.get());
    if (v.size() != 1 || !R_finite(v[0]) || v[0] < 0.0 || v[0] > 1.0) {
      Rcpp::stop("the minimum cell reuse ratio must be one number in [0, 1]");
    }
    emc::set_cell_min_reuse_ratio(v[0]);
  }
  return previous;
}

//' Enable, read, or reset mapper route counters
//'
//' The counters are disabled by default.  Passing enabled changes the switch;
//' reset clears counts, and the returned list is a snapshot.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List emc_pt_mapper_stats(
    Rcpp::Nullable<Rcpp::LogicalVector> enabled = R_NilValue,
    bool reset = false) {
  if (enabled.isNotNull()) {
    Rcpp::LogicalVector v(enabled.get());
    if (v.size() != 1 || Rcpp::LogicalVector::is_na(v[0]))
      Rcpp::stop("enabled must be one TRUE or FALSE");
    emc::set_mapper_stats_enabled(v[0] != 0);
  }
  if (reset) emc::reset_mapper_stats();
  const emc::MapperStatsSnapshot s = emc::mapper_stats_snapshot();
  return Rcpp::List::create(
    Rcpp::_["enabled"] = emc::mapper_stats_enabled(),
    Rcpp::_["map_scalar"] = static_cast<double>(s.map_scalar),
    Rcpp::_["map_cell"] = static_cast<double>(s.map_cell),
    Rcpp::_["map_row"] = static_cast<double>(s.map_row),
    Rcpp::_["transform_scalar"] = static_cast<double>(s.transform_scalar),
    Rcpp::_["transform_cell"] = static_cast<double>(s.transform_cell),
    Rcpp::_["transform_row"] = static_cast<double>(s.transform_row),
    Rcpp::_["bound_scalar"] = static_cast<double>(s.bound_scalar),
    Rcpp::_["bound_cell"] = static_cast<double>(s.bound_cell),
    Rcpp::_["bound_row"] = static_cast<double>(s.bound_row),
    Rcpp::_["cell_admit"] = static_cast<double>(s.cell_admit),
    Rcpp::_["cell_reject"] = static_cast<double>(s.cell_reject),
    Rcpp::_["cell_budget_reject"] = static_cast<double>(s.cell_budget_reject),
    Rcpp::_["cell_reuse_reject"] = static_cast<double>(s.cell_reuse_reject));
}

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
    Rcpp::NumericVector v(bytes.get());
    const double max_size = static_cast<double>(
      std::numeric_limits<std::size_t>::max());
    const bool max_is_exact =
      std::numeric_limits<std::size_t>::digits <=
      std::numeric_limits<double>::digits;
    if (v.size() != 1 || !R_finite(v[0]) || v[0] < 0 ||
        (max_is_exact ? v[0] > max_size : v[0] >= max_size)) {
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
    Rcpp::LogicalVector v(on.get());
    if (v.size() != 1 || Rcpp::LogicalVector::is_na(v[0])) {
      Rcpp::stop("emc_pt_defer_scalars() takes one TRUE or FALSE");
    }
    emc::set_defer_scalar_fill(v[0] != 0);
  }
  return previous;
}
