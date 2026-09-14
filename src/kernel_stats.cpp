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

KernelStat& pcounter_storage() {
  static KernelStat s;
  return s;
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

void pcounter_stats_record_raw(int route, long long rows, double seconds,
                               double preparation_seconds,
                               double rt_sum_seconds) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  if (route == PC_DENSITY) {
    s.dpcounter_calls += 1;
    s.dpcounter_rows += rows;
    s.dpcounter_seconds += seconds;
  } else if (route == PC_SURVIVOR) {
    s.ppcounter_calls += 1;
    s.ppcounter_rows += rows;
    s.ppcounter_seconds += seconds;
  } else {
    s.pcounter_logS_calls += 1;
    s.pcounter_logS_rows += rows;
    s.pcounter_logS_seconds += seconds;
  }
  s.preparation_seconds += preparation_seconds;
  s.rt_sum_seconds += rt_sum_seconds;
}

void pcounter_stats_record_branch(bool sv_zero, bool gamma_zero,
                                  bool omega_zero) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  s.sv_zero += sv_zero ? 1 : 0;
  s.gamma_zero += gamma_zero ? 1 : 0;
  s.omega_zero += omega_zero ? 1 : 0;
}

void pcounter_stats_record_exit(bool support) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  if (support) s.support_exits += 1;
  else s.invalid_exits += 1;
}

void pcounter_stats_record_fallback_tail() {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  s.fallback_tail += 1;
}

void pcounter_stats_record_k(int K) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  s.k_observations += 1;
  s.k_sum += K;
  if (s.k_observations == 1 || K < s.k_min) s.k_min = K;
  if (s.k_observations == 1 || K > s.k_max) s.k_max = K;
  if (K <= 3) ++s.k_0_3;
  else if (K <= 7) ++s.k_4_7;
  else if (K <= 15) ++s.k_8_15;
  else if (K <= 31) ++s.k_16_31;
  else if (K <= 63) ++s.k_32_63;
  else if (K <= 127) ++s.k_64_127;
  else if (K <= 255) ++s.k_128_255;
  else if (K <= 1023) ++s.k_256_1023;
  else ++s.k_1024_plus;
}

void pcounter_stats_add_stirling_terms(long long n) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  s.stirling_terms += n;
}

void pcounter_stats_add_rising_terms(long long n) {
  KernelStat& s = pcounter_storage();
  s.model = "PCOUNTER";
  s.rising_terms += n;
}

bool pcounter_stats_seen() {
  const KernelStat& s = pcounter_storage();
  return s.dpcounter_calls != 0 || s.ppcounter_calls != 0 ||
         s.pcounter_logS_calls != 0 || s.sv_zero != 0 ||
         s.gamma_zero != 0 || s.omega_zero != 0 ||
         s.fallback_tail != 0 || s.invalid_exits != 0 ||
         s.support_exits != 0 || s.k_observations != 0 ||
         s.stirling_terms != 0 || s.rising_terms != 0;
}

KernelStat pcounter_stats_snapshot() { return pcounter_storage(); }

void pcounter_stats_reset() { pcounter_storage() = KernelStat(); }

std::vector<KernelStat> kernel_stats_read() {
  std::vector<KernelStat> out;
  out.reserve(table().size() + (pcounter_stats_seen() ? 1U : 0U));
  for (std::map<std::string, KernelStat>::const_iterator it = table().begin();
       it != table().end(); ++it) {
    out.push_back(it->second);
  }
  if (pcounter_stats_seen()) {
    KernelStat pc = pcounter_stats_snapshot();
    for (std::vector<KernelStat>::iterator it = out.begin(); it != out.end(); ++it) {
      if (it->model == "PCOUNTER") {
        // Keep the ordinary particle/kernel timing row, while adding the
        // adapter-specific fields to the same model row.
        pc.calls = it->calls;
        pc.rows = it->rows;
        pc.cells = it->cells;
        pc.seconds = it->seconds;
        *it = pc;
        return out;
      }
    }
    out.push_back(pc);
  }
  return out;
}

void kernel_stats_reset() {
  table().clear();
  pcounter_stats_reset();
}


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

//' One row per model seen since the last reset.  `rows` is how many per-trial
//' evaluations the kernel performed and `cells` how many a cell-resolution
//' version would have performed, so `rows / cells` is the reuse available.
//' For PCOUNTER, the additional adapter, branch, K, term, and timing columns
//' are populated; they are zero on other model rows.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame emc_kernel_stats_read() {
  const std::vector<emc::KernelStat> s = emc::kernel_stats_read();
  const int n = static_cast<int>(s.size());
  Rcpp::CharacterVector model(n);
  Rcpp::NumericVector calls(n), rows(n), cells(n), seconds(n);
  Rcpp::NumericVector d_calls(n), d_rows(n), d_seconds(n);
  Rcpp::NumericVector p_calls(n), p_rows(n), p_seconds(n);
  Rcpp::NumericVector ls_calls(n), ls_rows(n), ls_seconds(n);
  Rcpp::NumericVector sv_zero(n), gamma_zero(n), omega_zero(n);
  Rcpp::NumericVector fallback_tail(n), invalid_exits(n), support_exits(n);
  Rcpp::NumericVector k_observations(n), k_sum(n), k_min(n), k_max(n);
  Rcpp::NumericVector k_0_3(n), k_4_7(n), k_8_15(n), k_16_31(n);
  Rcpp::NumericVector k_32_63(n), k_64_127(n), k_128_255(n);
  Rcpp::NumericVector k_256_1023(n), k_1024_plus(n);
  Rcpp::NumericVector stirling_terms(n), rising_terms(n);
  Rcpp::NumericVector preparation_seconds(n), rt_sum_seconds(n);
  for (int i = 0; i < n; ++i) {
    model[i] = s[i].model;
    calls[i] = static_cast<double>(s[i].calls);
    rows[i] = static_cast<double>(s[i].rows);
    cells[i] = static_cast<double>(s[i].cells);
    seconds[i] = s[i].seconds;
    d_calls[i] = static_cast<double>(s[i].dpcounter_calls);
    d_rows[i] = static_cast<double>(s[i].dpcounter_rows);
    d_seconds[i] = s[i].dpcounter_seconds;
    p_calls[i] = static_cast<double>(s[i].ppcounter_calls);
    p_rows[i] = static_cast<double>(s[i].ppcounter_rows);
    p_seconds[i] = s[i].ppcounter_seconds;
    ls_calls[i] = static_cast<double>(s[i].pcounter_logS_calls);
    ls_rows[i] = static_cast<double>(s[i].pcounter_logS_rows);
    ls_seconds[i] = s[i].pcounter_logS_seconds;
    sv_zero[i] = static_cast<double>(s[i].sv_zero);
    gamma_zero[i] = static_cast<double>(s[i].gamma_zero);
    omega_zero[i] = static_cast<double>(s[i].omega_zero);
    fallback_tail[i] = static_cast<double>(s[i].fallback_tail);
    invalid_exits[i] = static_cast<double>(s[i].invalid_exits);
    support_exits[i] = static_cast<double>(s[i].support_exits);
    k_observations[i] = static_cast<double>(s[i].k_observations);
    k_sum[i] = static_cast<double>(s[i].k_sum);
    k_min[i] = static_cast<double>(s[i].k_min);
    k_max[i] = static_cast<double>(s[i].k_max);
    k_0_3[i] = static_cast<double>(s[i].k_0_3);
    k_4_7[i] = static_cast<double>(s[i].k_4_7);
    k_8_15[i] = static_cast<double>(s[i].k_8_15);
    k_16_31[i] = static_cast<double>(s[i].k_16_31);
    k_32_63[i] = static_cast<double>(s[i].k_32_63);
    k_64_127[i] = static_cast<double>(s[i].k_64_127);
    k_128_255[i] = static_cast<double>(s[i].k_128_255);
    k_256_1023[i] = static_cast<double>(s[i].k_256_1023);
    k_1024_plus[i] = static_cast<double>(s[i].k_1024_plus);
    stirling_terms[i] = static_cast<double>(s[i].stirling_terms);
    rising_terms[i] = static_cast<double>(s[i].rising_terms);
    preparation_seconds[i] = s[i].preparation_seconds;
    rt_sum_seconds[i] = s[i].rt_sum_seconds;
  }
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::_["model"] = model, Rcpp::_["calls"] = calls,
    Rcpp::_["rows"] = rows, Rcpp::_["cells"] = cells,
    Rcpp::_["seconds"] = seconds, Rcpp::_["stringsAsFactors"] = false);
  out["dpcounter_calls"] = d_calls;
  out["dpcounter_rows"] = d_rows;
  out["dpcounter_seconds"] = d_seconds;
  out["ppcounter_calls"] = p_calls;
  out["ppcounter_rows"] = p_rows;
  out["ppcounter_seconds"] = p_seconds;
  out["pcounter_logS_calls"] = ls_calls;
  out["pcounter_logS_rows"] = ls_rows;
  out["pcounter_logS_seconds"] = ls_seconds;
  out["sv_zero"] = sv_zero;
  out["gamma_zero"] = gamma_zero;
  out["omega_zero"] = omega_zero;
  out["fallback_tail"] = fallback_tail;
  out["invalid_exits"] = invalid_exits;
  out["support_exits"] = support_exits;
  out["k_observations"] = k_observations;
  out["k_sum"] = k_sum;
  out["k_min"] = k_min;
  out["k_max"] = k_max;
  out["k_0_3"] = k_0_3;
  out["k_4_7"] = k_4_7;
  out["k_8_15"] = k_8_15;
  out["k_16_31"] = k_16_31;
  out["k_32_63"] = k_32_63;
  out["k_64_127"] = k_64_127;
  out["k_128_255"] = k_128_255;
  out["k_256_1023"] = k_256_1023;
  out["k_1024_plus"] = k_1024_plus;
  out["stirling_terms"] = stirling_terms;
  out["rising_terms"] = rising_terms;
  out["preparation_seconds"] = preparation_seconds;
  out["rt_sum_seconds"] = rt_sum_seconds;
  out.attr("class") = Rcpp::CharacterVector::create("data.frame");
  out.attr("row.names") = n == 0
    ? Rcpp::IntegerVector(0)
    : Rcpp::IntegerVector::create(NA_INTEGER, -n);
  return out;
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
