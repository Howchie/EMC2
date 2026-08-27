#ifndef fpe_race_h
#define fpe_race_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Making a grid-valued PDE solver behave like a pointwise race kernel.
//
// EMC2's race contract is pointwise: every kernel is handed one row per
// (trial x accumulator) and asked for f(rt | par) or S(rt | par) at that row's
// own rt.  A Fokker-Planck march is the opposite shape -- it amortises O(M*nt)
// work over a WHOLE time grid and answers every rt at once, so evaluating it
// once per row would cost several hundred times more than necessary.
//
// The bridge is a solve cache keyed on the parameter tuple.  The number of
// solves then equals the number of DISTINCT parameter rows, not the number of
// trials: a forstmann-style design (v~lM, s~lM, B~E+lR, A~1, t0~1) has 12 of
// them covering 524 compressed rows, so a subject costs ~12 solves rather than
// ~524.  A trend model makes every trial distinct and degenerates to one solve
// per row -- correct, just slow; see the cost message in make_emc().
//
// The cache is also what makes CENSORING affordable.  The common unknown-winner
// path (log_surv_cm, particle_ll.cpp) needs scalar CDF values at the censoring
// and truncation bounds, on parameter rows the batch kernels have already
// solved.  Those calls become interpolations instead of solves.
//
// Everything here is model-agnostic except rou_solve(); adding the BM, GBM or
// Gompertz race means another solve function and another Key layout.
//
// NOTE ON t0: it is deliberately NOT part of the cache key.  t0 shifts the time
// axis (tt = rt - t0) and does not change the solve, so two design cells that
// differ only in t0 share one march.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <type_traits>
#include <R_ext/Arith.h>
#if defined(__x86_64__) || defined(_M_X64) || defined(__AVX2__)
#include <immintrin.h>
#endif
#include "fpe_models.h"


namespace fperace {

#if defined(__AVX512F__) && defined(__FMA__)
constexpr size_t OU_BATCH_LANES = 8;
template <size_t LANES> struct OUSimd;
template <> struct OUSimd<8> {
  using Vec = __m512d;
  static inline Vec load(const double* p) { return _mm512_loadu_pd(p); }
  static inline void store(double* p, Vec x) { _mm512_storeu_pd(p, x); }
  static inline Vec set1(double x) { return _mm512_set1_pd(x); }
  static inline Vec mul(Vec a, Vec b) { return _mm512_mul_pd(a, b); }
  static inline Vec div(Vec a, Vec b) { return _mm512_div_pd(a, b); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm512_fmadd_pd(a, b, c);
  }
  static inline Vec fnmadd(Vec a, Vec b, Vec c) {
    return _mm512_fnmadd_pd(a, b, c);
  }
};
template <> struct OUSimd<4> {
  using Vec = __m256d;
  static inline Vec load(const double* p) { return _mm256_loadu_pd(p); }
  static inline void store(double* p, Vec x) { _mm256_storeu_pd(p, x); }
  static inline Vec set1(double x) { return _mm256_set1_pd(x); }
  static inline Vec mul(Vec a, Vec b) { return _mm256_mul_pd(a, b); }
  static inline Vec div(Vec a, Vec b) { return _mm256_div_pd(a, b); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm256_fmadd_pd(a, b, c);
  }
  static inline Vec fnmadd(Vec a, Vec b, Vec c) {
    return _mm256_fnmadd_pd(a, b, c);
  }
};
#elif defined(__AVX2__) && defined(__FMA__)
constexpr size_t OU_BATCH_LANES = 4;
template <size_t LANES> struct OUSimd;
template <> struct OUSimd<4> {
  using Vec = __m256d;
  static inline Vec load(const double* p) { return _mm256_loadu_pd(p); }
  static inline void store(double* p, Vec x) { _mm256_storeu_pd(p, x); }
  static inline Vec set1(double x) { return _mm256_set1_pd(x); }
  static inline Vec mul(Vec a, Vec b) { return _mm256_mul_pd(a, b); }
  static inline Vec div(Vec a, Vec b) { return _mm256_div_pd(a, b); }
  static inline Vec fmadd(Vec a, Vec b, Vec c) {
    return _mm256_fmadd_pd(a, b, c);
  }
  static inline Vec fnmadd(Vec a, Vec b, Vec c) {
    return _mm256_fnmadd_pd(a, b, c);
  }
};
#else
// Keep cache chunking bounded on non-x86 and baseline x86 builds. The solver
// below uses scalar loops in this configuration.
constexpr size_t OU_BATCH_LANES = 4;
#endif

// Log of the smallest density/survivor worth representing.  Grids are stored in
// LOG space and interpolated there, because both quantities are close to
// exponential in t over a grid interval -- linear interpolation of log S is
// dramatically better than of S once S is small, and that is precisely the
// regime a race likelihood spends its time in.
constexpr double LOG_FLOOR = -700.0;

inline double safe_log(double x) {
  return (x > 0.0 && std::isfinite(x)) ? std::max(std::log(x), LOG_FLOOR)
                                       : LOG_FLOOR;
}

// ---------------------------------------------------------------------------
// Resolution.  dt_target rather than a fixed nt: t_max is set by the data (and
// is extended when a scalar/GSL call needs a longer horizon), so pinning nt
// would make dt -- and therefore the accuracy -- a function of the horizon, and
// a re-solve to a longer horizon would silently change the likelihood of rows
// that had already been evaluated.  Pinning dt keeps re-solves consistent.
// ---------------------------------------------------------------------------
// Defaults set by the step-10 sweep (race_ou_integration_plan.md section 8):
// nx is the binding constraint and dt is nearly free once the time grid is
// graded, so the budget goes to nx.  At (384, 2e-3) the likelihood carries
// 1.1e-3 nats/trial of discretisation bias against a converged reference, and a
// fit is within 2e-4 nats/trial of the reference fit.
struct FPE_Grid {
  int nx = 512;
  double dt_target = 4e-3;
  int nt_min = 256;
  int nt_max = 4096;
  double grade = fpe::FPE_GRADE;
  double tgrade = fpe::FPE_TGRADE;

  int nt_for(double t_max) const {
    const double want = std::ceil(t_max / std::max(dt_target, 1e-6));
    const int n = (want < static_cast<double>(nt_max)) ? static_cast<int>(want)
                                                       : nt_max;
    return std::max(n, nt_min);
  }
};

// ---------------------------------------------------------------------------
// Cache key.  Bit-exact comparison is the right test: these doubles come from
// the same ParamTable columns for every row in a design cell, so equal cells
// compare equal and unequal cells must not share a solve.
// ---------------------------------------------------------------------------
struct Key {
  double v = 0.0, k = 0.0, b = 0.0, A = 0.0;
  double v_T = 0.0, tau_S = 0.0, tau_T = 0.0;

  // ROU uses a state-rescaled OU (sigma == 1, linear state).  Gompertz uses
  // the same OU march after Y = log(X), but must keep the physical diffusion
  // and the physical-start interval because A is defined on the X scale.
  // Keeping these fields on the cache key lets both models share the cache
  // machinery without pretending that a log-space uniform start is uniform.
  // ROUp channel mask: bit 0 sustained, bit 1 transient.  Other model kinds
  // leave the default value untouched.
  unsigned char roup_active = 3;
  double sigma = 1.0;
  double zlo = 0.0, zhi = 0.0;
  int model_kind = 0;              // 0 = ROU, 1 = Gompertz, 2 = ROUp
  bool log_state = false;

  // Boundary.  For FPE_BND_FIXED the shape vectors are forced to zero by
  // rou_key(), so every fixed-bound row compares equal on them and a model that
  // never collapses behaves exactly as it did before collapse existed.
  int bkind = fpe::FPE_BND_FIXED;
  double binf = 0.0, tau = 0.0, pw = 0.0;

  bool operator==(const Key& o) const {
    return v == o.v && k == o.k && b == o.b && A == o.A &&
           v_T == o.v_T && tau_S == o.tau_S && tau_T == o.tau_T &&
           sigma == o.sigma && zlo == o.zlo && zhi == o.zhi &&
           model_kind == o.model_kind && log_state == o.log_state &&
           roup_active == o.roup_active &&
           bkind == o.bkind && binf == o.binf && tau == o.tau && pw == o.pw;
  }
  bool finite() const {
    return std::isfinite(v) && std::isfinite(k) && std::isfinite(b) &&
           std::isfinite(A) && std::isfinite(v_T) && std::isfinite(tau_S) &&
           std::isfinite(tau_T) &&
           std::isfinite(sigma) && std::isfinite(zlo) &&
           std::isfinite(zhi) && std::isfinite(binf) && std::isfinite(tau) &&
           std::isfinite(pw);
  }
};

struct KeyHash {
  size_t operator()(const Key& key) const noexcept {
    size_t h = std::hash<double>{}(key.v);
    auto mix = [&](double x) {
      const size_t hx = std::hash<double>{}(x);
      h ^= hx + static_cast<size_t>(0x9e3779b9U) + (h << 6) + (h >> 2);
    };
    mix(key.k);
    mix(key.b);
    mix(key.A);
    mix(key.v_T);
    mix(key.tau_S);
    mix(key.tau_T);
    mix(key.sigma);
    mix(static_cast<double>(key.roup_active));
    mix(key.zlo);
    mix(key.zhi);
    mix(static_cast<double>(key.model_kind));
    mix(static_cast<double>(key.log_state));
    mix(static_cast<double>(key.bkind));
    mix(key.binf);
    mix(key.tau);
    mix(key.pw);
    return h;
  }
};

// Collapse spec on the MODEL's scale, before s-scaling.  Binf is the asymptotic
// boundary measured from zero (unlike B, which is measured from the top of the
// start-point range).  The solver's FPE_Boundary knows no parameter names; this
// is the one place that maps them.
struct BndSpec {
  int kind = fpe::FPE_BND_FIXED;
  double Binf = 0.0, tau = 0.0, pw = 0.0;
};

// ---------------------------------------------------------------------------
// Parameterisations.
//
// The alternatives are restricted charts on the same SDE family.  The map below
// is the only place they exist: everything downstream -- the key, the cache, the
// march, the simulator -- sees (v, k, s) and cannot tell which chart produced
// them.  RATE remains a three-assignment pass-through.
//
// CURVATURE (tstar, k, s) uses the deterministic path from x = 0 to
// b = B + A.  tstar is its reference crossing time; k and s retain their
// physical units and can therefore be shared across conditions.  This chart
// covers the deterministic-crossing regime and the k = 0 Wiener limit.
//
// EQUILIBRIUM (tk, theta, chi) uses physical state units.  theta is the actual OU
// equilibrium v/k and chi = s*sqrt(tk) is the diffusion scale accumulated over one
// relaxation.  In Y = X/(scale) and u = t/tk the process has the usual OU form;
// unlike the old threshold-normalised chart, changing B while theta/chi are held
// fixed therefore changes the distance to the bound rather than rescaling the
// entire accumulator.  It covers finite positive leak (tk > 0).
//
// Both alternatives leave one state-scale redundancy when their noise parameter
// is free, so designs should normally fix B rather than s.  See R/model_ROU.R.
// ---------------------------------------------------------------------------
enum : int {
  ROU_PAR_RATE = 0,
  ROU_PAR_CURVATURE = 1,
  ROU_PAR_EQUILIBRIUM = 2
};

// ROUp keeps the sustained channel in the same (v_S, tau_S) coordinates in
// every chart. Its transient alternative replaces v_T by E_T = v_T * tau_T;
// the solver and cache continue to use the canonical rate v_T.
enum : int {
  ROUP_PAR_RATE = 0,
  ROUP_PAR_AREA = 1
};

// (p1, p2, p3) are the parameterisation's own three columns, in p_types order.
// Unrepresentable input becomes NaN rather than a silently substituted value,
// so rou_key()'s finiteness test rejects the row exactly as it does for a bad
// rate row.
inline void rou_map_to_rate(int par_kind, double p1, double p2, double p3,
                            double B, double A, double& v, double& k,
                            double& s) {
  if (par_kind == ROU_PAR_RATE) { v = p1; k = p2; s = p3; return; }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double AA = (A > 0.0 && std::isfinite(A)) ? A : 0.0;
  const double b = B + AA;
  if (!(b > 0.0) || !std::isfinite(b)) { v = k = s = nan; return; }

  if (par_kind == ROU_PAR_CURVATURE) {
    const double tstar = p1, kk = p2, ss = p3;
    if (!(tstar > 0.0) || !std::isfinite(tstar) ||
        !(kk >= 0.0) || !std::isfinite(kk) ||
        !(ss >= 0.0) || !std::isfinite(ss)) {
      v = k = s = nan;
      return;
    }
    k = kk;
    // The denominator is well behaved for small k*tstar; the explicit k = 0
    // branch keeps the Wiener race exact rather than merely approximate.
    v = (kk > 0.0) ? b * kk / -std::expm1(-kk * tstar) : b / tstar;
    s = ss;
    return;
  }

  // ROU_PAR_EQUILIBRIUM
  const double tk = p1, theta = p2, chi = p3;
  if (!(tk > 0.0) || !std::isfinite(tk) || !std::isfinite(theta) ||
      !(chi >= 0.0) || !std::isfinite(chi)) {
    v = k = s = nan;
    return;
  }
  k = 1.0 / tk;
  v = theta / tk;                           // theta = v/k (physical equilibrium)
  s = chi / std::sqrt(tk);                  // chi = s*sqrt(tk)
}

// Build the key for one accumulator's parameters.  s is scaled out exactly as
// the RDM does (v/s, B/s, A/s, with k unchanged -- k has units of 1/time and is
// invariant under the state rescaling Y = X/s), which is what makes the k = 0
// member reproduce the analytic Wald parameterisation term for term.  Returns
// false for a row this model cannot represent, which the callers turn into a
// zero density / full survivor.
//
// It lives here rather than with the kernels so that the C++ likelihood path
// and the R-facing dROU/pROU cannot drift apart in how they form a key.
inline bool rou_key(double v, double k, double B, double A, double s,
                    const BndSpec& bs, Key& out) {
  if (!(s > 0.0) || !std::isfinite(s)) return false;
  const double inv_s = 1.0 / s;
  out.v = v * inv_s;
  out.k = (k > 0.0 && std::isfinite(k)) ? k : 0.0;
  out.sigma = 1.0;
  out.zlo = 0.0;
  out.model_kind = 0;
  out.log_state = false;
  out.roup_active = 3;
  const double AA = (A > 0.0 && std::isfinite(A)) ? A * inv_s : 0.0;
  out.A = AA;
  out.b = B * inv_s + AA;                  // b = B + A, both already scaled

  // The collapse asymptote is Binf itself, measured from ZERO -- not Binf + A.
  // The bound may therefore descend into the start-point range [0, A], which is
  // the point: a bound that meets the start point is a forced response, and
  // that is a state the model should be able to express.  Binf = 0 is allowed
  // outright: the domain is [x_lo, b(t)] with x_lo < 0, so b(t) -> 0 does not
  // degenerate it.
  out.bkind = bs.kind;
  if (bs.kind == fpe::FPE_BND_FIXED) {
    out.binf = 0.0; out.tau = 0.0; out.pw = 0.0;   // canonical, so keys compare
  } else {
    if (!(bs.Binf >= 0.0) || !std::isfinite(bs.Binf) ||
        !(bs.tau > 0.0) || !std::isfinite(bs.tau)) return false;
    out.binf = bs.Binf * inv_s;
    out.tau = bs.tau;
    // Only the Weibull form has a shape exponent; zero it otherwise so that two
    // rows differing in an unused column still share one solve.
    out.pw = (bs.kind == fpe::FPE_BND_WEIBULL) ? bs.pw : 0.0;
    if (bs.kind == fpe::FPE_BND_WEIBULL &&
        (!(bs.pw > 0.0) || !std::isfinite(bs.pw))) return false;
  }
  return out.finite() && out.b > 0.0 && out.v > 0.0;
}

inline bool roup_key(double v_S, double v_T, double tau_S, double tau_T, double k, double B, double A, double s,
                     const BndSpec& bs, Key& out) {
  if (!(s > 0.0) || !std::isfinite(s)) return false;
  if (!std::isfinite(v_S) || !std::isfinite(v_T) || v_S < 0.0 || v_T < 0.0)
    return false;
  const bool active_S = v_S > fpe::ROUP_DRIFT_EPS;
  const bool active_T = v_T > fpe::ROUP_DRIFT_EPS;
  if (!active_S && !active_T) return false;
  if (active_S && (!(tau_S > 0.0) || !std::isfinite(tau_S))) return false;
  if (active_T && (!(tau_T > 0.0) || !std::isfinite(tau_T))) return false;
  const double inv_s = 1.0 / s;
  out.v = active_S ? v_S * inv_s : 0.0; // repurpose v as v_S
  out.v_T = active_T ? v_T * inv_s : 0.0;
  // Disabled channels are canonicalised so their time constants cannot split
  // cache entries or affect a solve.
  out.tau_S = active_S ? tau_S : 0.0;
  out.tau_T = active_T ? tau_T : 0.0;
  out.k = (k > 0.0 && std::isfinite(k)) ? k : 0.0;
  out.sigma = 1.0;
  out.zlo = 0.0;
  out.model_kind = 2; // ROUp
  out.log_state = false;
  out.roup_active = static_cast<unsigned char>((active_S ? 1 : 0) |
                                               (active_T ? 2 : 0));
  const double AA = (A > 0.0 && std::isfinite(A)) ? A * inv_s : 0.0;
  out.A = AA;
  out.b = B * inv_s + AA;                  // b = B + A, both already scaled

  // Boundary
  out.bkind = bs.kind;
  if (bs.kind == fpe::FPE_BND_FIXED) {
    out.binf = 0.0; out.tau = 0.0; out.pw = 0.0;   // canonical, so keys compare
  } else {
    if (!(bs.Binf >= 0.0) || !std::isfinite(bs.Binf) ||
        !(bs.tau > 0.0) || !std::isfinite(bs.tau)) return false;
    out.binf = bs.Binf * inv_s;
    out.tau = bs.tau;
    out.pw = (bs.kind == fpe::FPE_BND_WEIBULL) ? bs.pw : 0.0;
    if (bs.kind == fpe::FPE_BND_WEIBULL &&
        (!(bs.pw > 0.0) || !std::isfinite(bs.pw))) return false;
  }
  return out.finite() && out.b > 0.0;
}

inline bool roup_transient_to_rate(int par_kind, double transient, double tau_T,
                                   double& v_T) {
  if (!std::isfinite(transient) || transient < 0.0) return false;
  if (par_kind == ROUP_PAR_RATE) {
    v_T = transient;
    return true;
  }
  if (par_kind != ROUP_PAR_AREA) return false;
  if (transient <= fpe::ROUP_DRIFT_EPS) {
    v_T = 0.0;
    return true;
  }
  if (!(tau_T > 0.0) || !std::isfinite(tau_T)) return false;
  v_T = transient / tau_T;
  return std::isfinite(v_T);
}

// Map a raw ROUp chart to the canonical pulse-rate chart before state scaling
// and key construction. Keeping this map next to roup_key makes the sampled
// likelihood, scalar censoring path, and R-facing vector path share one algebra.
inline bool roup_key_par(int par_kind, double v_S, double transient,
                         double tau_S, double tau_T, double k, double B,
                         double A, double s, const BndSpec& bs, Key& out) {
  double v_T = 0.0;
  if (!roup_transient_to_rate(par_kind, transient, tau_T, v_T)) return false;
  return roup_key(v_S, v_T, tau_S, tau_T, k, B, A, s, bs, out);
}

// Fixed-boundary overload for ROUp.
inline bool roup_key(double v_S, double v_T, double tau_S, double tau_T, double k, double B, double A, double s, Key& out) {
  static const BndSpec fixed;
  return roup_key(v_S, v_T, tau_S, tau_T, k, B, A, s, fixed, out);
}

// Fixed-boundary overload -- the common case, and every pre-collapse call site.
inline bool rou_key(double v, double k, double B, double A, double s, Key& out) {
  static const BndSpec fixed;
  return rou_key(v, k, B, A, s, fixed, out);
}

// Map first, then key.  Every call site that has raw parameter columns should
// use this rather than mapping itself, so that a new parameterisation is one
// edit to rou_map_to_rate() and nothing else.
inline bool rou_key_par(int par_kind, double p1, double p2, double p3,
                        double B, double A, const BndSpec& bs, Key& out) {
  double v, k, s;
  rou_map_to_rate(par_kind, p1, p2, p3, B, A, v, k, s);
  return rou_key(v, k, B, A, s, bs, out);
}

struct GridBlock {
  double t_start = 0.0;
  double t_end = 0.0;
  double dt = 0.0;
  double inv_dt = 0.0;
  size_t start_idx = 0;
};

struct GridIndex {
  std::vector<GridBlock> blocks;

  void build(const std::vector<double>& tg) {
    blocks.clear();
    const size_t n = tg.size();
    if (n < 2) return;

    size_t start = 0;
    while (start < n - 1) {
      const double dt = tg[start + 1] - tg[start];
      size_t end = start + 1;
      while (end + 1 < n) {
        const double next_dt = tg[end + 1] - tg[end];
        if (std::abs(next_dt - dt) > 1e-7 * (dt + 1e-12)) break;
        ++end;
      }
      GridBlock blk;
      blk.t_start = tg[start];
      blk.t_end = tg[end];
      blk.dt = dt;
      blk.inv_dt = (dt > 0.0) ? (1.0 / dt) : 0.0;
      blk.start_idx = start;
      blocks.push_back(blk);
      start = end;
    }
  }

  inline size_t find_index(const std::vector<double>& tg, double t) const {
    const size_t n = tg.size();
    if (n == 0) return 0;
    if (t <= tg.front()) return 0;
    if (t >= tg.back()) return n - 1;

    for (const auto& blk : blocks) {
      if (t >= blk.t_start && t <= blk.t_end) {
        if (blk.dt <= 0.0) break;
        int idx = static_cast<int>(blk.start_idx) +
                  static_cast<int>((t - blk.t_start) * blk.inv_dt) + 1;
        idx = std::max(1, std::min(idx, static_cast<int>(n) - 1));
        while (idx < static_cast<int>(n) - 1 && tg[idx] < t) ++idx;
        while (idx > 1 && tg[idx - 1] > t) --idx;
        return static_cast<size_t>(idx);
      }
    }
    return static_cast<size_t>(
        std::lower_bound(tg.begin(), tg.end(), t) - tg.begin());
  }
};

struct Entry {
  Key key;
  double t_max = 0.0;
  std::vector<double> t, log_pdf, log_S;
  GridIndex grid_idx;
  // Full entries contain every solver time point and may interpolate arbitrary
  // scalar queries. Sparse entries contain exact answers only for the finite
  // RTs requested by a raw batch.
  bool complete_grid = true;
};

struct SolveCache {
  FPE_Grid grid;
  // Boundary form for this model instance, set once when the adapter is built.
  // It lives here rather than on ContextForRaceModels because it is a property
  // of the solve, and because every path that can reach a solve already holds
  // the cache.
  int bnd_kind = fpe::FPE_BND_FIXED;
  // Parameterisation, likewise set once.  It selects which columns the kernels
  // read and how they map onto (v, k, s); the solve itself is unaffected.
  int par_kind = ROU_PAR_RATE;
  // ROUp's independent-channel variant keeps two row-to-entry maps while
  // sharing the generic solve entries and exact-key index.
  bool roup_local = false;
  bool sparse_raw_output = true;
  bool prepared = false;
  size_t n_entries = 0;
  std::vector<Entry> e;
  std::vector<int> row_group;   // scratch: row -> index into e, or -1
  std::vector<int> row_group_s; // local ROUp sustained entry, or -1
  std::vector<int> row_group_t; // local ROUp transient entry, or -1
  std::unordered_map<Key, int, KeyHash> index;

  // Cleared once per particle.  Keys are exact, so a stale entry could never be
  // returned for the wrong parameters; the clear exists to bound memory, since
  // a run visits thousands of particles.  Memory allocations of internal
  // vectors in e are preserved across particles to avoid heap churn.
  void new_particle() {
    prepared = false;
    for (size_t i = 0; i < n_entries; ++i) {
      e[i].t.clear();
      e[i].log_pdf.clear();
      e[i].log_S.clear();
      e[i].grid_idx.blocks.clear();
      e[i].complete_grid = true;
      e[i].t_max = 0.0;
    }
    n_entries = 0;
    row_group.clear();
    row_group_s.clear();
    row_group_t.clear();
    index.clear();
  }
};

// ---------------------------------------------------------------------------
// One racing-OU march.
//
//   dX = (v - k X) dt + s dW,   X_0 ~ U(0, A),   absorbed at b = B + A.
//
// s is scaled out exactly as the RDM does (Y = X/s carries v/s, b/s, A/s and
// the same k), so the k = 0 member reproduces the analytic Wald parameterisation
// term for term.
// ---------------------------------------------------------------------------
inline void rou_solve(const Key& p, double t_max, const FPE_Grid& gr, Entry& out) {
  fpe::FPE_ModelOU m;
  m.v = p.v;
  m.lambda = p.k;
  m.sigma = p.sigma;
  m.bnd.set_kind(p.bkind, p.b,
                 (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                 p.tau, p.pw, p.log_state);
  const double zlo = p.log_state ? p.zlo : 0.0;
  const double zhi = p.log_state ? p.zhi : p.A;
  m.xlo = fpe::fpe_x_lo_ou(zlo, p.v, p.k, p.sigma, t_max);

  const fpe::FPE_Result r =
      fpe::fpe_run(m, zlo, zhi, t_max, gr.nx, gr.nt_for(t_max), gr.grade,
                   gr.tgrade);

  const size_t n = r.t.size();
  out.key = p;
  out.t_max = t_max;
  out.complete_grid = true;
  out.t = r.t;
  out.log_pdf.resize(n);
  out.log_S.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.log_pdf[i] = safe_log(r.pdf[i]);
    out.log_S[i] = safe_log(r.surv[i]);
  }
  out.grid_idx.build(out.t);
}

inline void roup_solve(const Key& p, double t_max, const FPE_Grid& gr, Entry& out) {
  fpe::FPE_ModelPulseOU m;
  m.v_S = p.v;
  m.v_T = p.v_T;
  m.active_S = (p.roup_active & 1u) != 0u;
  m.active_T = (p.roup_active & 2u) != 0u;
  m.tau_S = p.tau_S;
  m.tau_T = p.tau_T;
  m.lambda = p.k;
  m.sigma = p.sigma;
  m.bnd.set_kind(p.bkind, p.b,
                 (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                 p.tau, p.pw, p.log_state);
  const double zlo = p.log_state ? p.zlo : 0.0;
  const double zhi = p.log_state ? p.zhi : p.A;
  double v_min = std::min(0.0, std::min(p.v, p.v + p.v_T)); 
  m.xlo = fpe::fpe_x_lo_ou(zlo, v_min, p.k, p.sigma, t_max);

  const fpe::FPE_Result r =
      fpe::fpe_run(m, zlo, zhi, t_max, gr.nx, gr.nt_for(t_max), gr.grade,
                   gr.tgrade);

  const size_t n = r.t.size();
  out.key = p;
  out.t_max = t_max;
  out.complete_grid = true;
  out.t = r.t;
  out.log_pdf.resize(n);
  out.log_S.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.log_pdf[i] = safe_log(r.pdf[i]);
    out.log_S[i] = safe_log(r.surv[i]);
  }
  out.grid_idx.build(out.t);
}

// Find the entry for `p` covering at least t_need, solving or extending as
// required.  Returns an INDEX, not a pointer: solving can reallocate `e`.
inline bool entry_has_query(const Entry& entry, double query_time) {
  if (!std::isfinite(query_time) || !(query_time > 0.0)) return false;
  if (entry.complete_grid) return entry.t_max >= query_time;
  return std::binary_search(entry.t.begin(), entry.t.end(), query_time);
}

inline int cache_get(SolveCache& C, const Key& p, double t_need,
                     double query_time = R_NaN) {
  if (!(t_need > 0.0)) t_need = 1e-3;
  auto it = C.index.find(p);
  if (it != C.index.end()) {
    const size_t i = static_cast<size_t>(it->second);
    if (entry_has_query(C.e[i], query_time)) return it->second;
    if (C.e[i].complete_grid && C.e[i].t_max >= t_need)
      return it->second;
    const double solve_to = std::max(t_need, C.e[i].t_max);
    if (p.model_kind == 2) roup_solve(p, solve_to, C.grid, C.e[i]);
    else rou_solve(p, solve_to, C.grid, C.e[i]);
    return it->second;
  }
  if (C.n_entries < C.e.size()) {
    const size_t idx = C.n_entries++;
    if (p.model_kind == 2) roup_solve(p, t_need, C.grid, C.e[idx]);
    else rou_solve(p, t_need, C.grid, C.e[idx]);
    C.index[p] = static_cast<int>(idx);
    return static_cast<int>(idx);
  } else {
    C.e.push_back(Entry());
    const size_t idx = C.n_entries++;
    if (p.model_kind == 2) roup_solve(p, t_need, C.grid, C.e[idx]);
    else rou_solve(p, t_need, C.grid, C.e[idx]);
    C.index[p] = static_cast<int>(idx);
    return static_cast<int>(idx);
  }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Multi-solve 4-lane interleaved Thomas solver for fixed-boundary OU models.
// Interleaves up to 4 solves in a Struct-of-Arrays (SoA) layout so the serial
// recurrence dependency chain in Thomas forward sweep and back-substitution
// executes 4 independent SIMD lanes concurrently, filling the 12-cycle latency
// shadow of CPU execution ports.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[[deprecated("common-clock batching does not preserve scalar OU schedules")]]
inline std::vector<fpe::FPE_Result> fpe_solve_batch_ou_common_clock(
    const std::vector<fpe::FPE_ModelOU>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const fpe::FPE_Mesh& g, int nt, double tgrade = 1.0) {
  const size_t num_models = models.size();
  std::vector<fpe::FPE_Result> results(num_models);
  if (num_models == 0) return results;

  double max_t_max = 0.0;
  for (size_t l = 0; l < num_models; ++l) {
    if (t_max_vec[l] > max_t_max) max_t_max = t_max_vec[l];
  }
  if (!(max_t_max > 0.0)) max_t_max = 1e-3;

  const int M = g.M;
  const fpe::FPE_TimeSchedule sched = fpe::fpe_time_schedule(max_t_max, nt, tgrade);

  alignas(32) std::vector<double> Q(M * 4, 0.0);
  alignas(32) std::vector<double> RHS(M * 4, 0.0);
  alignas(32) std::vector<double> SUB(M * 4, 0.0);
  alignas(32) std::vector<double> DINV(M * 4, 0.0);
  alignas(32) std::vector<double> CPRIME(M * 4, 0.0);
  alignas(32) std::vector<double> OP_DIAG(M * 4, 0.0);
  alignas(32) std::vector<double> OP_SUB(M * 4, 0.0);
  alignas(32) std::vector<double> OP_SUP(M * 4, 0.0);
  alignas(32) double OP_D[4] = {0.0, 0.0, 0.0, 0.0};
  alignas(32) double surv_prev[4] = {1.0, 1.0, 1.0, 1.0};

  for (size_t l = 0; l < 4; ++l) {
    const size_t src = (l < num_models) ? l : (num_models - 1);
    fpe::FPE_Op op;
    fpe::build_op(models[src], 0.0, g, op);
    OP_D[l] = op.D;

    for (int i = 0; i < M; ++i) {
      Q[i * 4 + l] = q0_vec[src][i];
      OP_DIAG[i * 4 + l] = op.diag[i];
      OP_SUB[i * 4 + l] = op.sub[i];
      OP_SUP[i * 4 + l] = op.sup[i];
    }

    double mass = 0.0;
    for (int i = 0; i < M; ++i) mass += g.dx[i] * q0_vec[src][i];
    surv_prev[l] = std::min(1.0, std::max(0.0, mass));

    if (l < num_models) {
      results[l].t.push_back(t0_vec[l]);
      double g0 = fpe::flux_out(op, q0_vec[src], g);
      if (!(g0 > 0.0)) g0 = 0.0;
      results[l].pdf.push_back(g0);
      results[l].surv.push_back(surv_prev[l]);
    }
  }

  int n_steps = 0;
  for (size_t b = 0; b < sched.steps.size(); ++b) n_steps += sched.steps[b];
  const int n_rann = std::min(2, n_steps);

  double t = 0.0;
  int done = 0;

  for (size_t b = 0; b < sched.dt.size(); ++b) {
    const double dtb = sched.dt[b];

    for (size_t l = 0; l < 4; ++l) {
      const size_t src = (l < num_models) ? l : (num_models - 1);
      fpe::FPE_Op op;
      fpe::build_op(models[src], 0.0, g, op);
      fpe::FPE_Tri tri;
      tri.resize(M);
      for (int i = 0; i < M; ++i) {
        tri.sub[i] = -0.5 * dtb * op.sub[i];
        tri.dia[i] = 1.0 - 0.5 * dtb * op.diag[i];
        tri.sup[i] = -0.5 * dtb * op.sup[i];
      }
      tri.factor();
      for (int i = 0; i < M; ++i) {
        SUB[i * 4 + l] = tri.sub[i];
        DINV[i * 4 + l] = tri.dinv[i];
        CPRIME[i * 4 + l] = tri.cprime[i];
      }
    }

    auto do_batch_step = [&](double step, bool backward_euler) {
      const double c = backward_euler ? step : (0.5 * step);
      t += step;

      // 1. Shifted operator step
      for (int i = 0; i < M; ++i) {
        if (backward_euler) {
          for (int l = 0; l < 4; ++l) RHS[i * 4 + l] = Q[i * 4 + l];
        } else {
#if defined(__AVX2__) && defined(__FMA__)
          __m256d q_c = _mm256_loadu_pd(&Q[i * 4]);
          __m256d op_d = _mm256_loadu_pd(&OP_DIAG[i * 4]);
          __m256d vc = _mm256_fmadd_pd(_mm256_set1_pd(c), _mm256_mul_pd(op_d, q_c), q_c);
          if (i > 0) {
            __m256d q_m = _mm256_loadu_pd(&Q[(i - 1) * 4]);
            __m256d op_s = _mm256_loadu_pd(&OP_SUB[i * 4]);
            vc = _mm256_fmadd_pd(_mm256_set1_pd(c), _mm256_mul_pd(op_s, q_m), vc);
          }
          if (i < M - 1) {
            __m256d q_p = _mm256_loadu_pd(&Q[(i + 1) * 4]);
            __m256d op_sp = _mm256_loadu_pd(&OP_SUP[i * 4]);
            vc = _mm256_fmadd_pd(_mm256_set1_pd(c), _mm256_mul_pd(op_sp, q_p), vc);
          }
          _mm256_storeu_pd(&RHS[i * 4], vc);
#else
          for (int l = 0; l < 4; ++l) {
            double v = Q[i * 4 + l] + c * OP_DIAG[i * 4 + l] * Q[i * 4 + l];
            if (i > 0)     v += c * OP_SUB[i * 4 + l] * Q[(i - 1) * 4 + l];
            if (i < M - 1) v += c * OP_SUP[i * 4 + l] * Q[(i + 1) * 4 + l];
            RHS[i * 4 + l] = v;
          }
#endif
        }
      }

      // 2. Forward elimination
      for (int l = 0; l < 4; ++l) {
        Q[0 * 4 + l] = RHS[0 * 4 + l] * DINV[0 * 4 + l];
      }
      for (int i = 1; i < M; ++i) {
#if defined(__AVX2__) && defined(__FMA__)
        __m256d r = _mm256_loadu_pd(&RHS[i * 4]);
        __m256d s = _mm256_loadu_pd(&SUB[i * 4]);
        __m256d q_prev = _mm256_loadu_pd(&Q[(i - 1) * 4]);
        __m256d d = _mm256_loadu_pd(&DINV[i * 4]);
        __m256d q_curr = _mm256_mul_pd(_mm256_fnmadd_pd(s, q_prev, r), d);
        _mm256_storeu_pd(&Q[i * 4], q_curr);
#else
        for (int l = 0; l < 4; ++l) {
          Q[i * 4 + l] = (RHS[i * 4 + l] - SUB[i * 4 + l] * Q[(i - 1) * 4 + l]) * DINV[i * 4 + l];
        }
#endif
      }

      // 3. Back substitution
      for (int i = M - 2; i >= 0; --i) {
#if defined(__AVX2__) && defined(__FMA__)
        __m256d q_next = _mm256_loadu_pd(&Q[(i + 1) * 4]);
        __m256d cp = _mm256_loadu_pd(&CPRIME[i * 4]);
        __m256d q_curr = _mm256_loadu_pd(&Q[i * 4]);
        q_curr = _mm256_fnmadd_pd(cp, q_next, q_curr);
        _mm256_storeu_pd(&Q[i * 4], q_curr);
#else
        for (int l = 0; l < 4; ++l) {
          Q[i * 4 + l] -= CPRIME[i * 4 + l] * Q[(i + 1) * 4 + l];
        }
#endif
      }

      // 4. Output recording
      for (size_t l = 0; l < num_models; ++l) {
        if (t > t_max_vec[l] + 1e-9 && !results[l].t.empty()) continue;

        double mass = 0.0;
        for (int i = 0; i < M; ++i) mass += g.dx[i] * Q[i * 4 + l];
        double s_keep = std::min(1.0, std::max(0.0, mass));
        if (s_keep > surv_prev[l]) s_keep = surv_prev[l];
        surv_prev[l] = s_keep;

        double gt = OP_D[l] * (g.cA * Q[(M - 1) * 4 + l] - g.cB * Q[(M - 2) * 4 + l]);
        if (!(gt > 0.0)) gt = 0.0;

        results[l].t.push_back(t);
        results[l].pdf.push_back(gt);
        results[l].surv.push_back(s_keep);
      }
    };

    for (int k = 0; k < sched.steps[b]; ++k, ++done) {
      if (done < n_rann) {
        do_batch_step(0.5 * dtb, true);
        do_batch_step(0.5 * dtb, true);
      } else {
        do_batch_step(dtb, false);
      }
    }
  }

  return results;
}

struct OUQueryResult {
  std::vector<double> t, log_pdf, log_S;
};

inline double interp_log_pair(double w, double value_lo, double value_hi) {
  const double lo = safe_log(value_lo);
  const double hi = safe_log(value_hi);
  if (lo <= LOG_FLOOR) return hi;
  if (hi <= LOG_FLOOR) return lo;
  return lo + w * (hi - lo);
}

// Each lane follows the same schedule fpe_solve() would have built for that
// key, including a point-start seed time and its exact terminal horizon.
// Independent Thomas recurrences remain interleaved, but their time-step
// coefficients are allowed to differ. AVX-512 builds process eight lanes;
// AVX2 builds process four, and other builds use the same four-lane layout with
// scalar inner loops. When query_times is supplied, the march is unchanged but
// only old-path log-interpolated answers at those times are retained.
template <class Model, size_t LANES, bool SPARSE = false>
inline std::vector<fpe::FPE_Result> fpe_solve_batch_ou_lanes(
    const std::vector<Model>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const std::vector<int>& nt_vec,
    const fpe::FPE_Mesh& g, double tgrade = 1.0,
    const std::vector<std::vector<double>>* query_times = nullptr,
    std::vector<OUQueryResult>* query_results = nullptr) {
  const size_t num_models = models.size();
  std::vector<fpe::FPE_Result> results(num_models);
  if constexpr (SPARSE)
    query_results->assign(num_models, OUQueryResult());
  if (num_models == 0) return results;

  struct Action {
    double step;
    double lhs_c;
    double rhs_c;
  };
  struct QueryBracket {
    double pdf_lo = 0.0, pdf_hi = 0.0;
    double survivor_lo = 1.0, survivor_hi = 1.0;
    double w = 0.0;
    // 0 = below seed, 1 = interpolate, 2 = terminal endpoint.
    int mode = 0;
  };

  const int M = g.M;
  std::vector<std::vector<Action>> actions(LANES);
  size_t max_actions = 0;

  const size_t storage_size = static_cast<size_t>(M) * LANES;
  std::vector<double> Q(storage_size, 0.0);
  std::vector<double> RHS(storage_size, 0.0);
  std::vector<double> SUB(storage_size, 0.0);
  std::vector<double> DINV(storage_size, 1.0);
  std::vector<double> CPRIME(storage_size, 0.0);
  std::vector<double> OP_DIAG(storage_size, 0.0);
  std::vector<double> OP_SUB(storage_size, 0.0);
  std::vector<double> OP_SUP(storage_size, 0.0);
  alignas(64) double op_D[LANES] = {};
  alignas(64) double survivor_prev[LANES] = {};
  alignas(64) double lane_pdf[LANES] = {};
  alignas(64) double lane_time[LANES] = {};
  size_t query_pos[LANES] = {};
  std::vector<std::vector<QueryBracket>> query_brackets;
  if constexpr (SPARSE) query_brackets.resize(num_models);

  for (size_t l = 0; l < LANES; ++l) {
    const size_t src = (l < num_models) ? l : (num_models - 1);
    // The seed is already at t0, so the initial flux must use the operator at
    // that time as well. This is immaterial for fixed-boundary OU, but matters
    // for ROUp where the pulse makes the operator time dependent.
    const double initial_time = (l < num_models) ? t0_vec[src] : 0.0;
    fpe::FPE_OpView op(&OP_SUB[l], &OP_DIAG[l], &OP_SUP[l], LANES);
    fpe::build_op(models[src], initial_time, g, op);
    op_D[l] = op.D;
    for (int i = 0; i < M; ++i)
      Q[static_cast<size_t>(i) * LANES + l] = q0_vec[src][i];

    double mass = 0.0;
    for (int i = 0; i < M; ++i) mass += g.dx[i] * q0_vec[src][i];
    survivor_prev[l] = std::min(1.0, std::max(0.0, mass));
    lane_time[l] = (l < num_models) ? t0_vec[l] : 0.0;

    if (l >= num_models) continue;
    double g0 = op_D[l] *
        (g.cA * Q[static_cast<size_t>(M - 1) * LANES + l] -
         g.cB * Q[static_cast<size_t>(M - 2) * LANES + l]);
    if (!(g0 > 0.0)) g0 = 0.0;
    lane_pdf[l] = g0;
    if constexpr (SPARSE) {
      const std::vector<double>& qt = (*query_times)[l];
      query_brackets[l].resize(qt.size());
      while (query_pos[l] < qt.size() && qt[query_pos[l]] <= t0_vec[l]) {
        QueryBracket& bracket = query_brackets[l][query_pos[l]];
        bracket.survivor_lo = survivor_prev[l];
        bracket.survivor_hi = survivor_prev[l];
        bracket.mode = 0;
        ++query_pos[l];
      }
    } else {
      results[l].t.push_back(t0_vec[l]);
      results[l].pdf.push_back(g0);
      results[l].surv.push_back(survivor_prev[l]);
    }

    const double span = std::max(t_max_vec[l] - t0_vec[l], 1e-12);
    const int nt = (l < nt_vec.size()) ? nt_vec[l] : 1;
    const fpe::FPE_TimeSchedule schedule =
        fpe::fpe_time_schedule(span, nt, tgrade);
    int n_steps = 0;
    for (size_t b = 0; b < schedule.steps.size(); ++b)
      n_steps += schedule.steps[b];
    const int n_rann = std::min(2, n_steps);
    actions[l].reserve(static_cast<size_t>(n_steps + n_rann));
    int done = 0;
    for (size_t b = 0; b < schedule.dt.size(); ++b) {
      const double dtb = schedule.dt[b];
      for (int k = 0; k < schedule.steps[b]; ++k, ++done) {
        if (done < n_rann) {
          actions[l].push_back({0.5 * dtb, 0.5 * dtb, 0.0});
          actions[l].push_back({0.5 * dtb, 0.5 * dtb, 0.0});
        } else {
          actions[l].push_back({dtb, 0.5 * dtb, 0.5 * dtb});
        }
      }
    }
    if constexpr (!SPARSE) {
      results[l].t.reserve(actions[l].size() + 1);
      results[l].pdf.reserve(actions[l].size() + 1);
      results[l].surv.reserve(actions[l].size() + 1);
    }
    max_actions = std::max(max_actions, actions[l].size());
  }

  auto set_identity_factor = [&](int lane) {
    for (int i = 0; i < M; ++i) {
      const size_t o = static_cast<size_t>(i) * LANES + lane;
      SUB[o] = 0.0;
      DINV[o] = 1.0;
      CPRIME[o] = 0.0;
    }
  };
  auto set_lane_factor = [&](int lane, double lhs_c) {
    double cprime_prev = 0.0;
    for (int i = 0; i < M; ++i) {
      const size_t o = static_cast<size_t>(i) * LANES + lane;
      const double sub = -lhs_c * OP_SUB[o];
      const double dia = 1.0 - lhs_c * OP_DIAG[o];
      const double sup = -lhs_c * OP_SUP[o];
      double denom = dia - ((i > 0) ? sub * cprime_prev : 0.0);
      if (denom == 0.0) denom = std::numeric_limits<double>::min();
      const double dinv = 1.0 / denom;
      SUB[o] = sub;
      DINV[o] = dinv;
      CPRIME[o] = (i + 1 < M) ? sup * dinv : 0.0;
      cprime_prev = CPRIME[o];
    }
  };

  alignas(64) double last_lhs[LANES] = {};
  bool factor_live[LANES] = {};
  for (size_t l = 0; l < LANES; ++l)
    set_identity_factor(static_cast<int>(l));

  bool is_moving[LANES] = {};
  for (size_t l = 0; l < num_models; ++l) {
    is_moving[l] = !models[l].static_op();
  }

  for (size_t aidx = 0; aidx < max_actions; ++aidx) {
    alignas(64) double rhs_c[LANES] = {};
    bool active[LANES] = {};
    for (size_t l = 0; l < LANES; ++l) {
      active[l] = l < num_models &&
                  aidx < actions[l].size();
      if (active[l]) {
        const Action& a = actions[l][aidx];
        rhs_c[l] = a.rhs_c;
        if (!is_moving[l] &&
            (!factor_live[l] || last_lhs[l] != a.lhs_c)) {
          set_lane_factor(static_cast<int>(l), a.lhs_c);
          last_lhs[l] = a.lhs_c;
          factor_live[l] = true;
        }
      } else if (factor_live[l]) {
        set_identity_factor(static_cast<int>(l));
        factor_live[l] = false;
      }
    }

    for (int i = 0; i < M; ++i) {
#if (defined(__AVX512F__) || defined(__AVX2__)) && defined(__FMA__)
      using Simd = OUSimd<LANES>;
      using Vec = typename Simd::Vec;
      const size_t o = static_cast<size_t>(i) * LANES;
      const Vec c = Simd::load(rhs_c);
      const Vec qc = Simd::load(&Q[o]);
      const Vec od = Simd::load(&OP_DIAG[o]);
      Vec value = Simd::fmadd(c, Simd::mul(od, qc), qc);
      if (i > 0) {
        const Vec qm = Simd::load(&Q[static_cast<size_t>(i - 1) * LANES]);
        const Vec os = Simd::load(&OP_SUB[o]);
        value = Simd::fmadd(c, Simd::mul(os, qm), value);
      }
      if (i + 1 < M) {
        const Vec qp = Simd::load(&Q[static_cast<size_t>(i + 1) * LANES]);
        const Vec os = Simd::load(&OP_SUP[o]);
        value = Simd::fmadd(c, Simd::mul(os, qp), value);
      }
      Simd::store(&RHS[o], value);
#else
      for (size_t l = 0; l < LANES; ++l) {
        const size_t o = static_cast<size_t>(i) * LANES + l;
        double value = Q[o] + rhs_c[l] * OP_DIAG[o] * Q[o];
        if (i > 0)
          value += rhs_c[l] * OP_SUB[o] *
                   Q[static_cast<size_t>(i - 1) * LANES + l];
        if (i + 1 < M)
          value += rhs_c[l] * OP_SUP[o] *
                   Q[static_cast<size_t>(i + 1) * LANES + l];
        RHS[o] = value;
      }
#endif
    }

    if constexpr (std::is_same_v<Model, fpe::FPE_ModelPulseOU>) {
      // Every ROUp lane has a moving operator. Once the CN RHS has consumed
      // the old operator, build the new operators and factor all lane-local
      // tridiagonals together; the Thomas dependency is only along i, so the
      // independent lanes are safe to vectorise here as well as in the solve.
      alignas(64) double lhs_vec[LANES] = {};
      for (size_t l = 0; l < LANES; ++l) {
        if (active[l]) {
          const Action& a = actions[l][aidx];
          lhs_vec[l] = a.lhs_c;
          const double t_new = lane_time[l] + a.step;
          fpe::FPE_OpView op(&OP_SUB[l], &OP_DIAG[l], &OP_SUP[l], LANES);
          fpe::build_op(models[l], t_new, g, op);
          op_D[l] = op.D;
        }
      }
#if (defined(__AVX512F__) || defined(__AVX2__)) && defined(__FMA__)
      using Simd = OUSimd<LANES>;
      using Vec = typename Simd::Vec;
      Vec cprime_prev = Simd::set1(0.0);
      const Vec one = Simd::set1(1.0);
      const Vec zero = Simd::set1(0.0);
      const Vec c = Simd::load(lhs_vec);
      for (int i = 0; i < M; ++i) {
        const size_t o = static_cast<size_t>(i) * LANES;
        const Vec os = Simd::load(&OP_SUB[o]);
        const Vec od = Simd::load(&OP_DIAG[o]);
        const Vec ou = Simd::load(&OP_SUP[o]);
        const Vec sub = Simd::fnmadd(c, os, zero);
        const Vec dia = Simd::fnmadd(c, od, one);
        const Vec sup = Simd::fnmadd(c, ou, zero);
        const Vec denom = Simd::fnmadd(sub, cprime_prev, dia);
        const Vec dinv = Simd::div(one, denom);
        const Vec cp = Simd::mul(sup, dinv);
        Simd::store(&SUB[o], sub);
        Simd::store(&DINV[o], dinv);
        Simd::store(&CPRIME[o], cp);
        cprime_prev = cp;
      }
      for (size_t l = 0; l < LANES; ++l) {
        if (active[l]) {
          last_lhs[l] = lhs_vec[l];
          factor_live[l] = true;
        }
      }
#else
      for (size_t l = 0; l < LANES; ++l) {
        if (!active[l]) continue;
        set_lane_factor(static_cast<int>(l), lhs_vec[l]);
        last_lhs[l] = lhs_vec[l];
        factor_live[l] = true;
      }
#endif
    } else {
      // Crank--Nicolson uses the operator at the previous time on its RHS and
      // the operator at the new time in the implicit solve. Update moving lanes
      // only after the RHS has been assembled so the batched solver follows the
      // scalar solver's time-centred discretisation.
      for (size_t l = 0; l < LANES; ++l) {
        if (!active[l] || !is_moving[l]) continue;
        const Action& a = actions[l][aidx];
        const double t_new = lane_time[l] + a.step;
        fpe::FPE_OpView op(&OP_SUB[l], &OP_DIAG[l], &OP_SUP[l], LANES);
        fpe::build_op(models[l], t_new, g, op);
        op_D[l] = op.D;
        set_lane_factor(static_cast<int>(l), a.lhs_c);
        last_lhs[l] = a.lhs_c;
        factor_live[l] = true;
      }
    }

    for (size_t l = 0; l < LANES; ++l)
      Q[l] = RHS[l] * DINV[l];
    for (int i = 1; i < M; ++i) {
#if (defined(__AVX512F__) || defined(__AVX2__)) && defined(__FMA__)
      using Simd = OUSimd<LANES>;
      using Vec = typename Simd::Vec;
      const size_t o = static_cast<size_t>(i) * LANES;
      const Vec r = Simd::load(&RHS[o]);
      const Vec s = Simd::load(&SUB[o]);
      const Vec qp = Simd::load(&Q[static_cast<size_t>(i - 1) * LANES]);
      const Vec d = Simd::load(&DINV[o]);
      const Vec qc = Simd::mul(Simd::fnmadd(s, qp, r), d);
      Simd::store(&Q[o], qc);
#else
      for (size_t l = 0; l < LANES; ++l) {
        const size_t o = static_cast<size_t>(i) * LANES + l;
        Q[o] = (RHS[o] - SUB[o] *
                Q[static_cast<size_t>(i - 1) * LANES + l]) * DINV[o];
      }
#endif
    }

    alignas(64) double mass[LANES];
#if (defined(__AVX512F__) || defined(__AVX2__)) && defined(__FMA__)
    using Simd = OUSimd<LANES>;
    using Vec = typename Simd::Vec;
    Vec massv = Simd::mul(
        Simd::set1(g.dx[M - 1]),
        Simd::load(&Q[static_cast<size_t>(M - 1) * LANES]));
    for (int i = M - 2; i >= 0; --i) {
      const size_t o = static_cast<size_t>(i) * LANES;
      const Vec qn = Simd::load(&Q[static_cast<size_t>(i + 1) * LANES]);
      const Vec cp = Simd::load(&CPRIME[o]);
      Vec qc = Simd::load(&Q[o]);
      qc = Simd::fnmadd(cp, qn, qc);
      Simd::store(&Q[o], qc);
      massv = Simd::fmadd(Simd::set1(g.dx[i]), qc, massv);
    }
    Simd::store(mass, massv);
#else
    for (size_t l = 0; l < LANES; ++l)
      mass[l] = g.dx[M - 1] *
          Q[static_cast<size_t>(M - 1) * LANES + l];
    for (int i = M - 2; i >= 0; --i) {
      for (size_t l = 0; l < LANES; ++l) {
        const size_t o = static_cast<size_t>(i) * LANES + l;
        Q[o] -= CPRIME[o] *
            Q[static_cast<size_t>(i + 1) * LANES + l];
        mass[l] += g.dx[i] * Q[o];
      }
    }
#endif

    for (size_t l = 0; l < num_models; ++l) {
      if (!active[l]) continue;
      const double prev_time = lane_time[l];
      const double prev_pdf = lane_pdf[l];
      const double prev_survivor = survivor_prev[l];
      const double current_time = prev_time + actions[l][aidx].step;
      double s_keep = std::min(1.0, std::max(0.0, mass[l]));
      if (s_keep > prev_survivor) s_keep = prev_survivor;
      double gt = op_D[l] *
          (g.cA * Q[static_cast<size_t>(M - 1) * LANES + l] -
           g.cB * Q[static_cast<size_t>(M - 2) * LANES + l]);
      if (!(gt > 0.0)) gt = 0.0;

      if constexpr (SPARSE) {
        const std::vector<double>& qt = (*query_times)[l];
        const bool final_action = aidx + 1 == actions[l].size();
        while (query_pos[l] < qt.size() &&
               qt[query_pos[l]] <= current_time) {
          const double tq = qt[query_pos[l]];
          QueryBracket& bracket = query_brackets[l][query_pos[l]];
          if (final_action && tq >= current_time) {
            // interp_log() treats the final grid point as an endpoint rather
            // than as the upper half of an interpolation interval.
            bracket.pdf_hi = gt;
            bracket.survivor_hi = s_keep;
            bracket.mode = 2;
          } else {
            bracket.pdf_lo = prev_pdf;
            bracket.pdf_hi = gt;
            bracket.survivor_lo = prev_survivor;
            bracket.survivor_hi = s_keep;
            bracket.w = (tq - prev_time) / (current_time - prev_time);
            bracket.mode = 1;
          }
          ++query_pos[l];
        }
        if (final_action) {
          // Floating summation can leave the terminal grid point a few ulps
          // below its requested horizon. The old lookup returns the endpoint
          // for every such query, so do the same here.
          while (query_pos[l] < qt.size()) {
            QueryBracket& bracket = query_brackets[l][query_pos[l]];
            bracket.pdf_hi = gt;
            bracket.survivor_hi = s_keep;
            bracket.mode = 2;
            ++query_pos[l];
          }
        }
      } else {
        results[l].t.push_back(current_time);
        results[l].pdf.push_back(gt);
        results[l].surv.push_back(s_keep);
      }
      lane_time[l] = current_time;
      lane_pdf[l] = gt;
      survivor_prev[l] = s_keep;
    }
  }

  if constexpr (SPARSE) {
    for (size_t l = 0; l < num_models; ++l) {
      const std::vector<double>& qt = (*query_times)[l];
      OUQueryResult& qr = (*query_results)[l];
      qr.t = qt;
      qr.log_pdf.resize(qt.size());
      qr.log_S.resize(qt.size());
      for (size_t j = 0; j < qt.size(); ++j) {
        const QueryBracket& bracket = query_brackets[l][j];
        if (bracket.mode == 0) {
          qr.log_pdf[j] = LOG_FLOOR;
          qr.log_S[j] = safe_log(bracket.survivor_hi);
        } else if (bracket.mode == 2) {
          qr.log_pdf[j] = safe_log(bracket.pdf_hi);
          qr.log_S[j] = safe_log(bracket.survivor_hi);
        } else {
          qr.log_pdf[j] = interp_log_pair(
              bracket.w, bracket.pdf_lo, bracket.pdf_hi);
          qr.log_S[j] = interp_log_pair(
              bracket.w, bracket.survivor_lo, bracket.survivor_hi);
        }
      }
    }
  }

  return results;
}

inline std::vector<fpe::FPE_Result> fpe_solve_batch_ou(
    const std::vector<fpe::FPE_ModelOU>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const std::vector<int>& nt_vec,
    const fpe::FPE_Mesh& g, double tgrade = 1.0) {
#if defined(__AVX512F__) && defined(__FMA__)
  if (models.size() <= 4) {
    return fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 4, false>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
  }
  return fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 8, false>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
#else
  return fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 4, false>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
#endif
}

inline std::vector<OUQueryResult> fpe_solve_batch_ou_queries(
    const std::vector<fpe::FPE_ModelOU>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const std::vector<int>& nt_vec,
    const std::vector<std::vector<double>>& query_times,
    const fpe::FPE_Mesh& g, double tgrade = 1.0) {
  std::vector<OUQueryResult> result;
#if defined(__AVX512F__) && defined(__FMA__)
  if (models.size() <= 4) {
    fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 4, true>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
        &query_times, &result);
  } else {
    fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 8, true>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
        &query_times, &result);
  }
#else
  fpe_solve_batch_ou_lanes<fpe::FPE_ModelOU, 4, true>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
      &query_times, &result);
#endif
  return result;
}

inline std::vector<fpe::FPE_Result> fpe_solve_batch_roup(
    const std::vector<fpe::FPE_ModelPulseOU>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const std::vector<int>& nt_vec,
    const fpe::FPE_Mesh& g, double tgrade = 1.0) {
#if defined(__AVX512F__) && defined(__FMA__)
  if (models.size() <= 4) {
    return fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 4, false>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
  }
  return fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 8, false>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
#else
  return fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 4, false>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade);
#endif
}

inline std::vector<OUQueryResult> fpe_solve_batch_roup_queries(
    const std::vector<fpe::FPE_ModelPulseOU>& models,
    const std::vector<std::vector<double>>& q0_vec,
    const std::vector<double>& t0_vec,
    const std::vector<double>& t_max_vec,
    const std::vector<int>& nt_vec,
    const std::vector<std::vector<double>>& query_times,
    const fpe::FPE_Mesh& g, double tgrade = 1.0) {
  std::vector<OUQueryResult> result;
#if defined(__AVX512F__) && defined(__FMA__)
  if (models.size() <= 4) {
    fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 4, true>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
        &query_times, &result);
  } else {
    fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 8, true>(
        models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
        &query_times, &result);
  }
#else
  fpe_solve_batch_ou_lanes<fpe::FPE_ModelPulseOU, 4, true>(
      models, q0_vec, t0_vec, t_max_vec, nt_vec, g, tgrade,
      &query_times, &result);
#endif
  return result;
}

inline bool entry_has_queries(const Entry& entry,
                              const std::vector<double>& query_times) {
  if (entry.complete_grid) return true;
  size_t j = 0;
  for (double t : query_times) {
    while (j < entry.t.size() && entry.t[j] < t) ++j;
    if (j >= entry.t.size() || entry.t[j] != t) return false;
  }
  return true;
}

// Batched cache lookup: solves missing keys at the widest SIMD width available
// to this build (eight lanes with AVX-512, four otherwise).
inline void cache_get_batch(SolveCache& C, const std::vector<Key>& keys,
                            const std::vector<double>& horizons,
                            std::vector<int>& out_cache_indices,
                            const std::vector<std::vector<double>>*
                                query_times = nullptr) {
  const size_t num_keys = keys.size();
  out_cache_indices.assign(num_keys, -1);
  if (num_keys == 0) return;

  std::vector<size_t> missing_keys;
  for (size_t i = 0; i < num_keys; ++i) {
    const Key& p = keys[i];
    double t_need = horizons[i];
    if (!(t_need > 0.0)) t_need = 1e-3;

    int found = -1;
    for (size_t j = 0; j < C.n_entries; ++j) {
      if (C.e[j].key == p && C.e[j].t_max >= t_need &&
          (C.e[j].complete_grid ||
           (query_times != nullptr &&
            entry_has_queries(C.e[j], (*query_times)[i])))) {
        found = static_cast<int>(j);
        break;
      }
    }
    if (found >= 0) {
      out_cache_indices[i] = found;
    } else {
      missing_keys.push_back(i);
    }
  }

  if (missing_keys.empty()) return;

  std::vector<size_t> standard_keys;
  std::vector<size_t> pulse_keys;
  for (size_t k_idx : missing_keys) {
    if (keys[k_idx].model_kind == 0) {
      standard_keys.push_back(k_idx);
    } else if (keys[k_idx].model_kind == 2) {
      pulse_keys.push_back(k_idx);
    } else {
      out_cache_indices[k_idx] = cache_get(C, keys[k_idx], horizons[k_idx]);
    }
  }

  // The normalised spatial mesh depends only on the configured grid, not on a
  // key.  Build it once for the whole cache-fill rather than once per SIMD chunk.
  if (standard_keys.empty() && pulse_keys.empty()) return;
  fpe::FPE_Mesh mesh;
  mesh.build(C.grid.nx, C.grid.grade);

  auto process_batches = [&](std::vector<size_t> group, auto make_model,
                             auto solve_full, auto solve_sparse) {
    if (group.empty()) return;

    // Similar horizons finish after a similar number of lane-local actions,
    // which minimises padded work after shorter lanes have completed.
    std::stable_sort(group.begin(), group.end(),
                     [&](size_t a, size_t b) {
                       return horizons[a] < horizons[b];
                     });

    using Model = std::decay_t<decltype(
        make_model(keys[group[0]], horizons[group[0]]))>;
    for (size_t b = 0; b < group.size(); b += OU_BATCH_LANES) {
      const size_t chunk_size = std::min(OU_BATCH_LANES, group.size() - b);
      std::vector<Model> models(chunk_size);
      std::vector<std::vector<double>> q0_vec(chunk_size);
      std::vector<double> t0_vec(chunk_size);
      std::vector<double> t_max_vec(chunk_size);
      std::vector<int> nt_vec(chunk_size);
      std::vector<std::vector<double>> chunk_queries;
      if (query_times != nullptr) chunk_queries.resize(chunk_size);
      size_t num_chunk_queries = 0;
      size_t num_chunk_steps = 0;

      for (size_t l = 0; l < chunk_size; ++l) {
        const size_t k_idx = group[b + l];
        const Key& p = keys[k_idx];
        models[l] = make_model(p, horizons[k_idx]);
        const double zlo = p.log_state ? p.zlo : 0.0;
        const double zhi = p.log_state ? p.zhi : p.A;

        const double t0 = fpe::fpe_seed(models[l], zlo, zhi, mesh,
                                        horizons[k_idx], q0_vec[l]);
        t0_vec[l] = t0;
        t_max_vec[l] = horizons[k_idx];
        nt_vec[l] = C.grid.nt_for(horizons[k_idx]);
        num_chunk_steps += static_cast<size_t>(nt_vec[l]);
        if (query_times != nullptr) {
          chunk_queries[l] = (*query_times)[k_idx];
          num_chunk_queries += chunk_queries[l].size();
        }
      }

      if (chunk_size > 1) {
        std::vector<fpe::FPE_Result> batch_res;
        std::vector<OUQueryResult> query_res;
        // Sparse bookkeeping is cheaper only while requested times are a small
        // fraction of the march. Above this crossover, contiguous full-grid
        // writes win despite retaining more output.
        const bool use_sparse =
            query_times != nullptr &&
            num_chunk_queries * 8 <= num_chunk_steps;
        if (use_sparse) {
          query_res = solve_sparse(models, q0_vec, t0_vec, t_max_vec, nt_vec,
                                   chunk_queries, mesh, C.grid.tgrade);
        } else {
          batch_res = solve_full(models, q0_vec, t0_vec, t_max_vec, nt_vec,
                                 mesh, C.grid.tgrade);
        }

        for (size_t l = 0; l < chunk_size; ++l) {
          const size_t k_idx = group[b + l];
          const Key& p = keys[k_idx];

          Entry entry;
          entry.key = p;
          entry.t_max = horizons[k_idx];
          if (use_sparse) {
            entry.complete_grid = false;
            entry.t = std::move(query_res[l].t);
            entry.log_pdf = std::move(query_res[l].log_pdf);
            entry.log_S = std::move(query_res[l].log_S);
          } else {
            const fpe::FPE_Result& r = batch_res[l];
            const size_t n = r.t.size();
            entry.complete_grid = true;
            entry.t = r.t;
            entry.log_pdf.resize(n);
            entry.log_S.resize(n);
            for (size_t j = 0; j < n; ++j) {
              entry.log_pdf[j] = safe_log(r.pdf[j]);
              entry.log_S[j] = safe_log(r.surv[j]);
            }
            entry.grid_idx.build(entry.t);
          }

          int existing = -1;
          auto it = C.index.find(p);
          if (it != C.index.end()) existing = it->second;

          if (existing >= 0) {
            C.e[existing] = std::move(entry);
            out_cache_indices[k_idx] = existing;
          } else if (C.n_entries < C.e.size()) {
            const size_t idx = C.n_entries++;
            C.e[idx] = std::move(entry);
            C.index[p] = static_cast<int>(idx);
            out_cache_indices[k_idx] = static_cast<int>(idx);
          } else {
            C.e.push_back(std::move(entry));
            C.index[p] = static_cast<int>(C.e.size() - 1);
            out_cache_indices[k_idx] = static_cast<int>(C.e.size() - 1);
            C.n_entries = C.e.size();
          }
        }
      } else {
        const size_t k_idx = group[b];
        out_cache_indices[k_idx] = cache_get(C, keys[k_idx], horizons[k_idx]);
      }
    }
  };

  auto make_ou_model = [](const Key& p, double t_max) {
    fpe::FPE_ModelOU m;
    m.v = p.v;
    m.lambda = p.k;
    m.sigma = p.sigma;
    m.bnd.set_kind(p.bkind, p.b,
                   (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                   p.tau, p.pw, p.log_state);
    const double zlo = p.log_state ? p.zlo : 0.0;
    m.xlo = fpe::fpe_x_lo_ou(zlo, p.v, p.k, p.sigma, t_max);
    return m;
  };
  auto make_roup_model = [](const Key& p, double t_max) {
    fpe::FPE_ModelPulseOU m;
    m.v_S = p.v;
    m.v_T = p.v_T;
    m.active_S = (p.roup_active & 1u) != 0u;
    m.active_T = (p.roup_active & 2u) != 0u;
    m.tau_S = p.tau_S;
    m.tau_T = p.tau_T;
    m.lambda = p.k;
    m.sigma = p.sigma;
    m.bnd.set_kind(p.bkind, p.b,
                   (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                   p.tau, p.pw, p.log_state);
    const double zlo = p.log_state ? p.zlo : 0.0;
    const double v_min = std::min(0.0, std::min(p.v, p.v + p.v_T));
    m.xlo = fpe::fpe_x_lo_ou(zlo, v_min, p.k, p.sigma, t_max);
    return m;
  };
  auto solve_ou_full = [](const auto& models, const auto& q0,
                          const auto& t0, const auto& tmax,
                          const auto& nt, const fpe::FPE_Mesh& g,
                          double tgrade) {
    return fpe_solve_batch_ou(models, q0, t0, tmax, nt, g, tgrade);
  };
  auto solve_ou_sparse = [](const auto& models, const auto& q0,
                            const auto& t0, const auto& tmax,
                            const auto& nt, const auto& queries,
                            const fpe::FPE_Mesh& g, double tgrade) {
    return fpe_solve_batch_ou_queries(models, q0, t0, tmax, nt, queries, g,
                                      tgrade);
  };
  auto solve_roup_full = [](const auto& models, const auto& q0,
                            const auto& t0, const auto& tmax,
                            const auto& nt, const fpe::FPE_Mesh& g,
                            double tgrade) {
    return fpe_solve_batch_roup(models, q0, t0, tmax, nt, g, tgrade);
  };
  auto solve_roup_sparse = [](const auto& models, const auto& q0,
                              const auto& t0, const auto& tmax,
                              const auto& nt, const auto& queries,
                              const fpe::FPE_Mesh& g, double tgrade) {
    return fpe_solve_batch_roup_queries(models, q0, t0, tmax, nt, queries, g,
                                        tgrade);
  };

  process_batches(standard_keys, make_ou_model, solve_ou_full,
                  solve_ou_sparse);
  process_batches(pulse_keys, make_roup_model, solve_roup_full,
                  solve_roup_sparse);
}

// ---------------------------------------------------------------------------
// Interpolation, in log space.
//
// Below the grid the density is zero (the march is seeded at t_seed > 0, so the
// grid does not reach back to 0) but the survivor is not -- cache_get's grid
// carries the mass absorbed before t_seed.  This is the same asymmetry
// fpe::grid_lookup handles with zero_below_grid.
// ---------------------------------------------------------------------------
inline double interp_log(const Entry& en, double t, bool floor_below) {
  const auto& tg = en.t;
  const size_t n = tg.size();
  if (n == 0) return LOG_FLOOR;
  if (t <= tg.front()) return floor_below ? LOG_FLOOR : (floor_below ? en.log_pdf.front() : en.log_S.front());
  if (t >= tg.back()) return floor_below ? en.log_pdf.back() : en.log_S.back();

  const size_t j = en.grid_idx.find_index(tg, t);
  const double lo = floor_below ? en.log_pdf[j - 1] : en.log_S[j - 1];
  const double hi = floor_below ? en.log_pdf[j] : en.log_S[j];
  // A floored neighbour carries no information; do not average garbage into a
  // live value, take the live one.
  if (lo <= LOG_FLOOR) return hi;
  if (hi <= LOG_FLOOR) return lo;
  const double w = (t - tg[j - 1]) / (tg[j] - tg[j - 1]);
  return lo + w * (hi - lo);
}

inline double sparse_entry_log(const Entry& en, double tt, bool density) {
  const auto it = std::lower_bound(en.t.begin(), en.t.end(), tt);
  if (it == en.t.end() || *it != tt)
    throw std::runtime_error("ROU sparse cache queried at an unprepared time");
  const size_t i = static_cast<size_t>(it - en.t.begin());
  return density ? en.log_pdf[i] : en.log_S[i];
}

inline double entry_log_pdf(const Entry& en, double tt) {
  return en.complete_grid ? interp_log(en, tt, true)
                          : sparse_entry_log(en, tt, true);
}
inline double entry_log_S(const Entry& en, double tt) {
  return en.complete_grid ? interp_log(en, tt, false)
                          : sparse_entry_log(en, tt, false);
}

// ---------------------------------------------------------------------------
// Reference simulator: dX = (v - k X) dt + s dW, X_0 ~ U(0,A), absorbed at b.
//
// Written here rather than reusing simulate_ou_hit_times_bb() because that one
// is parameterised by (lambda, theta) and rejects lambda <= 0, so it cannot
// reach the k = 0 member at all -- which is the one case with an independent
// analytic answer to check against.
//
// Exact OU transition per step (no Euler discretisation of the drift), with the
// k -> 0 limits taken through expm1 so k = 0 is a regular value:
//   phi   = exp(-k dt)
//   mean  = X phi + v * (1 - phi)/k        -> X + v dt
//   var   = s^2 (1 - phi^2)/(2k)           -> s^2 dt
// plus the Brownian-bridge correction for a crossing that happens strictly
// inside a step, which is what stops the hit times being biased late:
//   P(cross | X_t, X_{t+dt}) = exp(-2 (b - X_t)(b - X_{t+dt}) / (s^2 dt))
// ---------------------------------------------------------------------------
// A COLLAPSING bound needs nothing new here beyond evaluating b at both ends of
// the step.  Over one step the barrier is treated as linear, and the crossing
// probability for a Brownian bridge against a linear barrier is the same
// expression with the two endpoint gaps -- so the correction stays exact to the
// order the rest of the step already is, rather than degrading to "compare
// against b at one end", which biases hit times late exactly when the bound is
// moving fastest.
template <class Unif01, class Norm01>
inline double rou_hit_time_bnd(double v, double k, double s, double A,
                               const fpe::FPE_Boundary& bnd, double dt,
                               double t_max, Unif01 runif, Norm01 rnorm) {
  double X = (A > 0.0) ? A * runif() : 0.0;
  double b = bnd.b(0.0);
  if (X >= b) return 0.0;

  const double phi = std::exp(-k * dt);
  const double m1 = -std::expm1(-k * dt);            // 1 - phi
  const double m2 = -std::expm1(-2.0 * k * dt);      // 1 - phi^2
  const double drift_gain = (k > 1e-10) ? (m1 / k) : dt;
  const double drift_step = v * drift_gain;
  const double var = (k > 1e-10) ? (s * s * m2 / (2.0 * k)) : (s * s * dt);
  const double sd = std::sqrt(std::max(var, 0.0));
  const double inv_2var_bb = 2.0 / (s * s * dt);     // local BB uses s^2 dt

  double t = 0.0;
  while (t < t_max) {
    const double X1 = X * phi + drift_step + sd * rnorm();
    const double b1 = bnd.fixed ? b : bnd.b(t + dt);
    t += dt;
    if (X1 >= b1) {
      // Where the path and the barrier, both linear over the step, meet.
      const double d0 = b - X, d1 = b1 - X1;
      const double frac = d0 / std::max(d0 - d1, 1e-300);
      return t - dt + std::min(std::max(frac, 0.0), 1.0) * dt;
    }
    // Did the path cross and come back inside the step?
    const double pc = std::exp(-(b - X) * (b1 - X1) * inv_2var_bb);
    if (runif() < pc) return t - dt + runif() * dt;
    X = X1;
    b = b1;
  }
  return std::numeric_limits<double>::infinity();  // never finished in horizon
}

// Fixed-boundary overload, for the k = 0 / k > 0 cross-validation.
template <class Unif01, class Norm01>
inline double rou_hit_time(double v, double k, double s, double A, double b,
                           double dt, double t_max, Unif01 runif, Norm01 rnorm) {
  fpe::FPE_Boundary bnd;
  bnd.set_kind(fpe::FPE_BND_FIXED, b, b, 0.0, 0.0, false);
  return rou_hit_time_bnd(v, k, s, A, bnd, dt, t_max, runif, rnorm);
}

} // namespace fperace

#endif // fpe_race_h
