#ifndef EMC2_SCRATCH_H
#define EMC2_SCRATCH_H

// One owning, growable arena for per-call working memory.
//
// Why this exists.  The mapper's cell path kept its three accumulators as
// fixed C arrays sized by a constant, `EMC2_PT_MAX_CELLS = 256`, and the
// transform and bound lanes each kept their own.  The constant was not a
// modelling decision: it was how much would fit on the stack.  A design with
// 257 cells therefore fell off a cliff onto the per-row route -- 113.50 ms
// against 54.83 ms for identical results on the audit's RDM benchmark -- and
// wide designs are exactly where that hurts, because the number of design
// cells is what grows.
//
// Raising the constant would move the cliff, not remove it, and would put more
// on the stack in a function reached from a dozen translation units.  What is
// needed is scratch that is:
//
//   owning              it holds its own memory rather than borrowing the
//                       stack, so the size is a budget rather than a limit
//   dynamically sized   the size follows the design, not a constant
//   per-call            everything taken inside a frame is released when the
//                       frame ends, so nothing accumulates across calls
//   reused              the buffer persists between calls, so the common case
//                       costs a pointer bump rather than an allocation
//
// The lifetime is the part worth being strict about, so a frame is the only
// way in.  `ScratchFrame` marks the arena on construction and rewinds it on
// destruction; taking memory outside one is not expressible.
//
// The one rule: the arena may only grow when nothing is outstanding, because
// growing moves the buffer and would dangle every span already handed out.
// `reserve()` therefore succeeds only in an outermost frame, and a caller that
// is refused takes its slower general route rather than proceeding with a
// pointer it cannot trust.  A nested frame is refused for the same reason,
// which is safe rather than merely allowed: the fallback computes the same
// numbers.
//
// Not thread-safe, deliberately.  One arena belongs to one ParamTable and is
// used by one thread at a time; a threaded native backend (the audit's C11 and
// C15) needs one arena per worker, not a lock around this one.  Putting it on
// the table rather than in a global is what makes that a change of ownership
// rather than a rewrite.

#include <cstddef>
#include <cstdint>
#include <vector>
namespace emc {

class ScratchFrame;

class Scratch {
public:
  std::size_t capacity() const { return buf_.size(); }
  std::size_t used() const { return used_; }

private:
  friend class ScratchFrame;

  bool grow(std::size_t bytes) {
    if (bytes <= buf_.size()) return true;
    // Outstanding spans point into `buf_`; resizing may move it.
    if (used_ != 0) return false;
    buf_.resize(bytes);
    return true;
  }

  void* bump(std::size_t bytes, std::size_t align) {
    const std::size_t off = (used_ + align - 1) & ~(align - 1);
    if (bytes > buf_.size() || off > buf_.size() - bytes) return nullptr;
    used_ = off + bytes;
    return buf_.data() + off;
  }

  void rewind(std::size_t mark) { used_ = mark; }

  std::vector<unsigned char> buf_;
  std::size_t used_ = 0;
};

class ScratchFrame {
public:
  explicit ScratchFrame(Scratch& s) : s_(s), mark_(s.used()) {}
  ~ScratchFrame() { s_.rewind(mark_); }
  ScratchFrame(const ScratchFrame&) = delete;
  ScratchFrame& operator=(const ScratchFrame&) = delete;

  // Ask for room for everything this frame will take, before taking any of it.
  // FALSE means the caller must not use the arena in this frame.
  bool reserve(std::size_t bytes) { return s_.grow(mark_ + bytes); }

  // NULL when the reservation did not cover this; callers check once, after
  // taking everything, and fall back as a group.
  template <typename T>
  T* take(std::size_t n) {
    if (n == 0) return nullptr;
    return static_cast<T*>(s_.bump(n * sizeof(T), alignof(T)));
  }

private:
  Scratch& s_;
  std::size_t mark_;
};

// How much per-call cell scratch a single design entry may claim, in bytes.
// A budget, not a cell count: what it admits depends on how much each cell costs,
// so the mapper does not have to encode a number of cells anywhere.
// Settable so that the fallback route above the budget stays reachable in a
// test without building a design of a few hundred thousand cells.
std::size_t cell_scratch_budget();
void set_cell_scratch_budget(std::size_t bytes);

// Minimum fraction of trials that must be reused by a cell design.  Zero is
// the historical admission rule (any compressed design with fewer cells than
// trials is eligible), so merely exposing this knob does not change a route.
double cell_min_reuse_ratio();
void set_cell_min_reuse_ratio(double ratio);

// Whether a sampled coefficient's value is written once at row 0 and widened
// only if something reads it at row resolution, or written to every trial up
// front as it always was.  On by default; the switch exists so the two can be
// compared directly, which is how "they are the same numbers" is asserted and
// how the saving is measured.
bool defer_scalar_fill();
void set_defer_scalar_fill(bool on);

enum MapperRoute { MAPPER_SCALAR = 0, MAPPER_CELL = 1, MAPPER_ROW = 2 };

// These counters are deliberately plain native state.  Incrementing while
// disabled is a single branch and never allocates an R object; snapshots are
// only materialised by the diagnostic wrapper.
struct MapperStatsSnapshot {
  std::uint64_t map_scalar = 0;
  std::uint64_t map_cell = 0;
  std::uint64_t map_row = 0;
  std::uint64_t transform_scalar = 0;
  std::uint64_t transform_cell = 0;
  std::uint64_t transform_row = 0;
  std::uint64_t bound_scalar = 0;
  std::uint64_t bound_cell = 0;
  std::uint64_t bound_row = 0;
  std::uint64_t cell_admit = 0;
  std::uint64_t cell_reject = 0;
  std::uint64_t cell_budget_reject = 0;
  std::uint64_t cell_reuse_reject = 0;
};

bool mapper_stats_enabled();
void set_mapper_stats_enabled(bool on);
void reset_mapper_stats();
MapperStatsSnapshot mapper_stats_snapshot();
void mapper_count_map(MapperRoute route);
void mapper_count_transform(MapperRoute route);
void mapper_count_bound(MapperRoute route);
void mapper_count_cell_admit();
void mapper_count_cell_reject(bool budget, bool reuse);

}  // namespace emc

#endif
