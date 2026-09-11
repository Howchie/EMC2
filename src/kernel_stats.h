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
