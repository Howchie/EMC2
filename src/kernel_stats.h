#ifndef EMC2_KERNEL_STATS_H
#define EMC2_KERNEL_STATS_H

// Is there kernel arithmetic worth reusing, and how much of the call is it?
//
// The audit names three candidates -- BAwL's `natural_normalizer`, RDM's
// scale/geometry in `drdm_raw`/`prdm_raw`, and the DDM's `a/s`, `v/s`, `sv/s`
// scaling with `std::log(a)` -- and asks for measurement before any of them is
// specialised, because a subexpression only pays to hoist if it is recomputed
// far more often than it changes.
//
// That ratio is not a guess: a subexpression that reads a known set of
// parameter columns is constant exactly on the joint partition of those
// columns, which is what C8's `joint_cells()` returns.  So the measurement is
//
//   rows    per-trial evaluations the kernel actually performs
//   cells   evaluations a cell-resolution version would perform
//
// and rows/cells is the reuse available.  A ratio near 1 means the arithmetic
// varies per trial and specialising it would add indirection for nothing; the
// plan's instruction to abandon individual models is meant to be acted on with
// this number in hand.
//
// `ns` is wall time inside the kernel, so the saving can be put against the
// call rather than reported as a ratio with no denominator.
//
// Off by default and gated at the dispatch boundary: when it is on, the caller
// accepts a clock read per particle and one cached partition lookup.

#include <string>
#include <vector>

namespace emc {

struct KernelStat {
  std::string model;
  long long calls = 0;   // kernel invocations, i.e. particles
  long long rows = 0;    // per-trial evaluations performed
  long long cells = 0;   // evaluations at cell resolution
  double seconds = 0.0;  // wall time inside the kernel

  // PCOUNTER's raw adapters have a deliberately separate set of counters.
  // They are carried on the model row returned by emc_kernel_stats_read(),
  // rather than allocating an R object from a hot kernel.
  long long dpcounter_calls = 0;
  long long dpcounter_rows = 0;
  double dpcounter_seconds = 0.0;
  long long ppcounter_calls = 0;
  long long ppcounter_rows = 0;
  double ppcounter_seconds = 0.0;
  long long pcounter_logS_calls = 0;
  long long pcounter_logS_rows = 0;
  double pcounter_logS_seconds = 0.0;
  long long sv_zero = 0;
  long long gamma_zero = 0;
  long long omega_zero = 0;
  long long fallback_tail = 0;
  long long invalid_exits = 0;
  long long support_exits = 0;
  long long k_observations = 0;
  long long k_sum = 0;
  long long k_min = 0;
  long long k_max = 0;
  long long k_0_3 = 0;
  long long k_4_7 = 0;
  long long k_8_15 = 0;
  long long k_16_31 = 0;
  long long k_32_63 = 0;
  long long k_64_127 = 0;
  long long k_128_255 = 0;
  long long k_256_1023 = 0;
  long long k_1024_plus = 0;
  long long stirling_terms = 0;
  long long rising_terms = 0;
  double preparation_seconds = 0.0;
  double rt_sum_seconds = 0.0;
};

// PCOUNTER instrumentation is updated with plain native counters.  These
// helpers are no-ops at call sites when kernel_stats_on() is false; in
// particular they never construct Rcpp objects.
void pcounter_stats_record_raw(int route, long long rows, double seconds,
                               double preparation_seconds,
                               double rt_sum_seconds);
void pcounter_stats_record_branch(bool sv_zero, bool gamma_zero,
                                  bool omega_zero);
void pcounter_stats_record_exit(bool support);
void pcounter_stats_record_fallback_tail();
void pcounter_stats_record_k(int K);
void pcounter_stats_add_stirling_terms(long long n);
void pcounter_stats_add_rising_terms(long long n);
bool pcounter_stats_seen();
KernelStat pcounter_stats_snapshot();
void pcounter_stats_reset();

// Adapter identifiers used by pcounter_stats_record_raw.
enum PcounterStatsRoute {
  PC_DENSITY = 0,
  PC_SURVIVOR = 1,
  PC_LOG_SURVIVOR_AT_T = 2
};

bool kernel_stats_on();
void set_kernel_stats(bool on);

// `cells < 0` means the partition was not available (a covariate drives it to
// one cell per trial, or the table could not answer); it is recorded as equal
// to `rows`, i.e. no reuse, rather than left out.
void kernel_stats_record(const std::string& model, long long rows,
                         long long cells, double seconds);

std::vector<KernelStat> kernel_stats_read();
void kernel_stats_reset();

// The parameter columns a model's reusable subexpression reads, by the model's
// `c_name`.  This is the only model-dependent part of the measurement, and it
// is a declaration rather than a heuristic: each entry says which columns the
// candidate arithmetic is a function of, and nothing else.  An unknown model
// returns an empty set and is not measured.
std::vector<std::string> kernel_reuse_columns(const std::string& c_name);

}  // namespace emc

#endif
