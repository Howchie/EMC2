#ifndef fpe_bou_h
#define fpe_bou_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Bounded OU (Smith & Ratcliff 2004) -- solve cache and across-trial variability.
//
// This is to fpe::FPE_ModelBoundedOU what fpe_race.h is to fpe::FPE_ModelOU: it
// owns the mapping from a trial's parameters to a cached Fokker-Planck march and
// the interpolation back out, so the likelihood kernel never sees a solve.
//
// It is a good deal simpler than the race version in one respect and harder in
// another.  Simpler: there is no race to combine, because one march produces
// BOTH response densities -- the flux out of each absorbing barrier.  Harder:
// across-trial variability is genuinely expensive here, for a reason specific to
// this model:
//
//   * sv (drift, normal)      changes the drift    => new operator => new solve
//   * sz (start, uniform)     changes the START.  With the legacy start anchor
//                             this also changes the operator; with the midpoint
//                             anchor the operator is unchanged and the whole
//                             interval is seeded in one solve.
//   * st0 (non-decision)      a convolution on the finished density => free
//
// With the start anchor the cost is n_sv * n_sz solves per distinct parameter
// row.  The midpoint anchor instead uses one interval-seeded solve per drift
// node, so the start-point quadrature disappears exactly.
//
// Consequence for defaults: this model spends its budget very differently from
// the race models, and must NOT inherit fperace::FPE_Grid.  Two measurements
// (see bounded_ou_plan.md) drive that -- the error here is time-dominated rather
// than space-dominated, and grading the mesh actively hurts because a bounded
// domain has no far field to economise on.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <cfloat>
#include <unordered_map>
#if defined(__AVX512F__) || defined(__AVX2__)
#include <immintrin.h>
#endif
#include "fpe_models.h"
#include "gh_quad.h"
#include "gl_quad.h"

namespace fpebou {

// ---------------------------------------------------------------------------
// Grid settings.  Deliberately not fperace::FPE_Grid.
//
// nt is set from a target dt rather than fixed, so extending a solve to a longer
// horizon does not silently change the accuracy of rows already evaluated -- the
// same argument fperace::FPE_Grid makes, and the one place the two agree.
// ---------------------------------------------------------------------------
// Grid and quadrature are balanced against each other rather than each set as
// fine as it can afford.  Measured against the DDM oracle at beta = 0 with
// sv = 1, SZ = 0.3, st0 = 0.1 (worst |dpdf| over the RT range):
//
//   nodes   3x3      5x5      7x7      9x9        <- quadrature alone
//   err   4.1e-3   2.1e-4   3.4e-5   3.4e-5
//
//   at 7x7:  nx=512/dt=2.5e-4  3.4e-5  0.336 s
//            nx=384/dt=5e-4    6.1e-5  0.129 s
//            nx=256/dt=1e-3    1.4e-4  0.049 s
//
// At 5x5 the grid is irrelevant in the measured parameter net (every setting
// from nx=512 to nx=256 returns 2.05e-4, i.e. the error is all quadrature).
// The midpoint anchor removes the n_sz solve multiplier, so the defaults use
// five nodes in each quadrature and a 256-cell grid with a 1 ms target step.
struct Grid {
  int nx = 256;
  double dt_target = 1e-3;
  int nt_min = 512;
  int nt_max = 16384;
  double grade = fpe::FPE_GRADE_BOUNDED;   // 1.0: uniform, measured best
  // A shorter grading ratio avoids coarse late-time interpolation at the
  // balanced 1 ms target step without materially changing solve cost.
  double tgrade = 8.0;
  int n_sv = 5;                 // Gauss-Hermite nodes for drift variability
  int n_sz = 5;                 // Gauss-Legendre nodes for start-point variability
  int n_st0 = 7;                // Gauss-Legendre nodes for the t0 convolution

  int nt_for(double t_max) const {
    const int n = static_cast<int>(std::ceil(t_max / dt_target));
    return std::max(nt_min, std::min(nt_max, n));
  }
};

// ---------------------------------------------------------------------------
// Cache key.  One entry per (drift, leak, separation, start interval, anchor,
// boundary) tuple, compared BIT-EXACTLY: these values are read from the same
// ParamTable
// columns for every row in a design cell, so equal parameters really do produce
// identical doubles.  Quadrature nodes are folded in before the key is built, so
// a node is just another key -- the cache dedupes across nodes for free when two
// rows share one.
//
// s is scaled out exactly as the DDM does (v/s, a/s, z/s): the first-passage
// TIME of the rescaled process is the same random variable, so no Jacobian is
// needed on the density.
// ---------------------------------------------------------------------------
struct Key {
  double v = 0.0, beta = 0.0, a = 1.0;
  double z_lo = 0.5, z_hi = 0.5, anchor = 0.5;
  int bkind = fpe::FPE_BND_FIXED;
  double binf = 0.0, tau = 0.0, pw = 0.0;

  bool operator==(const Key& o) const {
    return v == o.v && beta == o.beta && a == o.a &&
           z_lo == o.z_lo && z_hi == o.z_hi &&
           anchor == o.anchor && bkind == o.bkind && binf == o.binf &&
           tau == o.tau && pw == o.pw;
  }
};

// Build a key from natural-scale parameters, scaling s out.  Returns false for
// a parameter set the solver cannot answer, so the caller floors the row rather
// than solving nonsense.
inline bool bou_key(double v, double beta, double a, double z_lo, double z_hi,
                    double anchor,
                    double s, int bkind, double binf, double tau, double pw,
                    Key& out) {
  if (!(s > 0.0) || !R_FINITE(s)) return false;
  if (!(a > 0.0) || !R_FINITE(a)) return false;
  if (!R_FINITE(v) || !R_FINITE(beta) || beta < 0.0) return false;
  if (!R_FINITE(z_lo) || !R_FINITE(z_hi) || z_lo > z_hi) return false;
  if (!(z_hi > 0.0) || !(z_lo < a)) return false;
  if (!R_FINITE(anchor)) return false;

  const double inv = 1.0 / s;
  out.v = v * inv;
  out.beta = beta;                 // a rate: invariant under the state rescale
  out.a = a * inv;
  out.z_lo = z_lo * inv;
  out.z_hi = z_hi * inv;
  out.anchor = anchor * inv;
  out.bkind = bkind;
  if (bkind == fpe::FPE_BND_FIXED) {
    // Canonicalise the shape fields so every fixed-bound row compares equal
    // regardless of what the (unused) collapse columns happen to hold.
    out.binf = 0.0; out.tau = 0.0; out.pw = 0.0;
  } else {
    out.binf = binf * inv; out.tau = tau; out.pw = pw;
  }
  return true;
}

// Hash for the solve index below.  Mirrors fperace::KeyHash: the fields are
// compared bit-exactly, so they are hashed bit-exactly too.
struct KeyHash {
  size_t operator()(const Key& key) const noexcept {
    size_t h = std::hash<double>{}(key.v);
    auto mix = [&](double x) {
      const size_t hx = std::hash<double>{}(x);
      h ^= hx + static_cast<size_t>(0x9e3779b9U) + (h << 6) + (h >> 2);
    };
    mix(key.beta);
    mix(key.a);
    mix(key.z_lo);
    mix(key.z_hi);
    mix(key.anchor);
    mix(static_cast<double>(key.bkind));
    mix(key.binf);
    mix(key.tau);
    mix(key.pw);
    return h;
  }
};

// ---------------------------------------------------------------------------
// A solved entry.  Holds BOTH responses; `lower` is R's first factor level and
// `upper` its second, matching dDDM/pDDM.
//
// Stored in log space for the same reason fperace::Entry does it: the consumer
// is a log-likelihood, and the defective cdfs are needed for truncation, where
// taking logs late loses the tail.
// ---------------------------------------------------------------------------
struct Entry {
  Key key;
  double t_max = 0.0;
  std::vector<double> t;
  std::vector<double> log_pdf_lo, log_pdf_up;
  std::vector<double> cdf_lo, cdf_up;   // natural scale: they are mixed linearly
  fperace::GridIndex grid_idx;
};

struct SolveCache {
  Grid grid;
  int bnd_kind = fpe::FPE_BND_FIXED;
  size_t n_entries = 0;
  std::vector<Entry> e;

  // Horizon every solve is taken to, set ONCE by the caller from the largest
  // decision time it will ask about.  This is not an optimisation detail: a
  // cache that grew its horizon on demand would re-solve every quadrature node
  // each time a later RT arrived, so walking a sorted RT vector would cost
  // O(n_rows) marches instead of O(1).  model_ROU.h:104-107 makes the same two-
  // pass argument -- the horizon of a group is only known once every row in it
  // has been seen.
  double t_horizon = 0.0;

  // Key -> entry slot.  find() was a LINEAR scan over the entries, and it is
  // called once per (drift x start-interval) node per (row x t0 node).  The
  // legacy start-anchor path has n_sv * n_sz entries per design cell; the
  // midpoint path has n_sv entries because its interval seed is exact.
  //
  // MEASURED: on forstmann with six parameter cells (294 entries) this index
  // made no measurable difference -- the per-row cost there is the quadrature
  // bookkeeping, not the lookup.  It is kept because the scan is the one part
  // of the loop whose cost grows quadratically with the design, which is the
  // same argument fperace::SolveCache made when it replaced its own linear scan
  // (see cache_get_batch); do not expect it to show up on a small design.
  std::unordered_map<Key, int, KeyHash> index;

  // Cleared once per particle.  Keys are exact, so a stale entry could never be
  // returned for different parameters -- but the parameters DO change every
  // particle, so keeping them would only grow the scan.
  void new_particle() { n_entries = 0; index.clear(); }

  Entry* find(const Key& k) {
    const auto it = index.find(k);
    if (it == index.end()) return nullptr;
    const size_t i = static_cast<size_t>(it->second);
    return (i < n_entries) ? &e[i] : nullptr;
  }
  // The key must be supplied here rather than left to whoever fills the entry:
  // the index is what makes find() O(1), and an entry that never reached it
  // would be re-solved on every lookup.
  Entry& push(const Key& k) {
    if (n_entries >= e.size()) e.resize(n_entries + 1);
    index[k] = static_cast<int>(n_entries);
    return e[n_entries++];
  }
  // Grow the store so the next `n` push() calls cannot reallocate, which would
  // invalidate pointers a batch is still holding.
  void reserve_for(size_t n) {
    if (n_entries + n > e.size()) e.resize(n_entries + n);
  }
};

// ---------------------------------------------------------------------------
// One march for one key.
// ---------------------------------------------------------------------------
// Build the solver model from a key.  `binf` is the ASYMPTOTIC SEPARATION, in
// the same units as `a`; set_separation turns the (a, a_inf) pair into the
// (midpoint, half-width) form the solver marches in.
inline fpe::FPE_ModelBoundedOU bou_model(const Key& p) {
  fpe::FPE_ModelBoundedOU m;
  m.v = p.v;
  m.beta = p.beta;
  m.anchor = p.anchor;
  m.sigma = 1.0;                 // s scaled out in bou_key
  m.set_separation(p.a, p.bkind, p.binf, p.tau, p.pw);
  return m;
}

inline void bou_solve(const Key& p, double t_max, const Grid& gr, Entry& out) {
  const fpe::FPE_ModelBoundedOU m = bou_model(p);

  fpe::FPE_Mesh g;
  g.build(gr.nx, gr.grade, fpe::FPE_ModelBoundedOU::symmetric_mesh);
  std::vector<double> q0;
  double absorbed_lower = 0.0;
  // The key carries either a point start (legacy start anchor) or the complete
  // uniform start interval (midpoint anchor).
  const double t0 = fpe::fpe_seed(m, p.z_lo, p.z_hi, g, t_max, q0,
                                  &absorbed_lower);
  const fpe::FPE_Result r =
    fpe::fpe_solve(m, q0, t0, t_max, g, gr.nt_for(t_max), gr.tgrade, absorbed_lower);

  out.key = p;
  out.t_max = t_max;
  out.t = r.t;
  out.cdf_lo = r.cdf_lower;
  out.cdf_up = r.cdf;
  const size_t n = r.t.size();
  out.log_pdf_lo.resize(n);
  out.log_pdf_up.resize(n);
  for (size_t i = 0; i < n; ++i) {
    out.log_pdf_lo[i] = (r.pdf_lower[i] > 0.0) ? std::log(r.pdf_lower[i]) : R_NegInf;
    out.log_pdf_up[i] = (r.pdf[i] > 0.0) ? std::log(r.pdf[i]) : R_NegInf;
  }
  out.grid_idx.build(out.t);
}

// ---------------------------------------------------------------------------
// Lane-batched march.
//
// The quadrature is the whole cost of this model, and its nodes are an unusually
// good fit for batching: they share nx, horizon and time schedule exactly, so
// they march on ONE clock.  That is what makes this simpler than the race
// version (fpe_race.h:597), whose accumulators each need their own schedule and
// retire at different times -- the reason its common-clock variant is
// deprecated.  Here a common clock is not an approximation, it is the truth.
//
// The layout is lane-interleaved (cell i, lane l at [i*LANES + l]) so the inner
// loops are contiguous over lanes and the compiler auto-vectorises them under
// -O3 -march=native.  That is deliberately not hand-written intrinsics: it gets
// most of the win for a fraction of the code, and every lane here runs the
// identical instruction sequence, which is the case auto-vectorisation handles
// well.  All bounds are fixed for now, so the operator is constant and the
// factorisation is done once per schedule block.
// ---------------------------------------------------------------------------
// Width follows the widest available register, as fpe_race.h:44-88 does.
// Measured on this batch (118 rows, sv+SZ+st0, 49 nodes, nx=512):
//   scalar 0.855 s   4 lanes 0.458 s   8 lanes 0.365 s
#if defined(__AVX512F__)
constexpr size_t BOU_LANES = 8;
template <size_t LANES> struct BOUSimd;
template <> struct BOUSimd<8> {
  using Vec = __m512d;
  static inline Vec load(const double* p) { return _mm512_loadu_pd(p); }
  static inline void store(double* p, Vec x) { _mm512_storeu_pd(p, x); }
  static inline Vec set1(double x) { return _mm512_set1_pd(x); }
  static inline Vec mul(Vec a, Vec b) { return _mm512_mul_pd(a, b); }
  static inline Vec sub(Vec a, Vec b) { return _mm512_sub_pd(a, b); }
  static inline Vec div(Vec a, Vec b) { return _mm512_div_pd(a, b); }
};
#elif defined(__AVX2__)
constexpr size_t BOU_LANES = 4;
template <size_t LANES> struct BOUSimd;
template <> struct BOUSimd<4> {
  using Vec = __m256d;
  static inline Vec load(const double* p) { return _mm256_loadu_pd(p); }
  static inline void store(double* p, Vec x) { _mm256_storeu_pd(p, x); }
  static inline Vec set1(double x) { return _mm256_set1_pd(x); }
  static inline Vec mul(Vec a, Vec b) { return _mm256_mul_pd(a, b); }
  static inline Vec sub(Vec a, Vec b) { return _mm256_sub_pd(a, b); }
  static inline Vec div(Vec a, Vec b) { return _mm256_div_pd(a, b); }
};
#else
constexpr size_t BOU_LANES = 4;
#endif

// ---------------------------------------------------------------------------
// Build the tridiagonal operator for every lane at time t, straight into the
// interleaved layout.
//
// This mirrors fpe::build_op exactly -- same stencil, same Scharfetter-Gummel
// face weights, same uniform-mesh Peclet recurrence -- with the lane loop
// innermost, so each row is a contiguous vector operation and no per-lane
// fpe::FPE_Op is materialised and then scattered.
//
// It exists for the COLLAPSING case.  With fixed bounds the operator is built
// once for the whole march and this is called exactly once; with bounds that
// move, it is called once per time step, alongside a fresh factorisation, and
// the two together are what a collapse costs over and above a fixed bound.
// ---------------------------------------------------------------------------
inline void bou_build_op_lanes(const fpe::FPE_ModelBoundedOU* models,
                               double t, const fpe::FPE_Mesh& g,
                               double* DIAG, double* SUB, double* SUP,
                               double* opD) {
  const int M = g.M;
  double a0[BOU_LANES], a1[BOU_LANES], invD[BOU_LANES], Dl[BOU_LANES];
  double Pu[BOU_LANES], eP[BOU_LANES], er[BOU_LANES], dP[BOU_LANES];

  // The uniform-mesh fast path is taken only if EVERY lane qualifies, so the
  // inner loops stay branch-free.  Lanes in a batch are quadrature nodes of one
  // trial, so they differ only in drift and start point and in practice either
  // all qualify or none do.
  bool uni = g.uniform;
  for (size_t l = 0; l < BOU_LANES; ++l) {
    const double Lb = models[l].length(t);
    const double Lp = models[l].length_prime(t);
    const double Bt = models[l].B() / Lb;
    const double D = 0.5 * Bt * Bt;
    Dl[l] = D; opD[l] = D; invD[l] = 1.0 / D;
    models[l].atil_affine(t, Lb, Lp, a0[l], a1[l]);
    Pu[l] = (a0[l] + a1[l] * g.h) * g.h * invD[l];
    dP[l] = a1[l] * g.h * g.h * invD[l];
    if (!(std::isfinite(Pu[l]) && std::isfinite(dP[l]) &&
          std::abs(dP[l]) * M < 20.0))
      uni = false;
    eP[l] = std::exp(Pu[l]);
    er[l] = std::exp(dP[l]);
  }

  double wm_p[BOU_LANES], wp_p[BOU_LANES], gf_p[BOU_LANES];
  for (size_t l = 0; l < BOU_LANES; ++l) wm_p[l] = wp_p[l] = gf_p[l] = 0.0;

  for (int i = 0; i < M; ++i) {
    const double rdx = g.rdx[i];
    const size_t o = static_cast<size_t>(i) * BOU_LANES;

    // Inflow face: interior for i > 0, the lower ABSORBING face at i == 0.
    if (i > 0) {
      for (size_t l = 0; l < BOU_LANES; ++l) {
        SUB[o + l]  = gf_p[l] * wm_p[l] * rdx;
        DIAG[o + l] = -gf_p[l] * wp_p[l] * rdx;
        SUP[o + l]  = 0.0;
      }
    } else {
      for (size_t l = 0; l < BOU_LANES; ++l) {
        SUB[l]  = 0.0;
        DIAG[l] = -Dl[l] * g.cL_A * rdx;
        SUP[l]  = Dl[l] * g.cL_B * rdx;
      }
    }

    // Outflow face: interior for i < M-1, the upper absorbing face at i == M-1.
    if (i < M - 1) {
      const int j = i + 1;
      const double rdc = g.rdc[j], pu = g.pu[j], pv = g.pv[j];
      if (uni) {
        for (size_t l = 0; l < BOU_LANES; ++l) {
          const double P = Pu[l];
          const double wp = fpe::bern_from_exp(P, eP[l]);
          Pu[l] += dP[l];
          eP[l] *= er[l];
          const double wm = fpe::bern_neg(wp, P);
          const double gf = Dl[l] * rdc;
          DIAG[o + l] -= gf * wm * rdx;
          SUP[o + l]  += gf * wp * rdx;
          wm_p[l] = wm; wp_p[l] = wp; gf_p[l] = gf;
        }
      } else {
        for (size_t l = 0; l < BOU_LANES; ++l) {
          const double P = (a0[l] * pu + a1[l] * pv) * invD[l];
          const double wp = fpe::bern(P);
          const double wm = fpe::bern_neg(wp, P);
          const double gf = Dl[l] * rdc;
          DIAG[o + l] -= gf * wm * rdx;
          SUP[o + l]  += gf * wp * rdx;
          wm_p[l] = wm; wp_p[l] = wp; gf_p[l] = gf;
        }
      }
    } else {
      for (size_t l = 0; l < BOU_LANES; ++l) {
        DIAG[o + l] -= Dl[l] * g.cA * rdx;
        SUB[o + l]  += Dl[l] * g.cB * rdx;
      }
    }
  }
}

inline void bou_solve_batch(const std::vector<Key>& keys, double t_max,
                            const Grid& gr, std::vector<Entry*>& out) {
  const size_t nk = keys.size();
  if (nk == 0) return;

  fpe::FPE_Mesh g;
  g.build(gr.nx, gr.grade, fpe::FPE_ModelBoundedOU::symmetric_mesh);
  const int M = g.M;

  for (size_t base = 0; base < nk; base += BOU_LANES) {
    const size_t n_lane = std::min(BOU_LANES, nk - base);

    // Pass 1: learn each lane's preferred seed time, then adopt the smallest.
    std::vector<fpe::FPE_ModelBoundedOU> models(BOU_LANES);
    std::vector<std::vector<double>> q0(BOU_LANES);
    std::vector<double> absorbed(BOU_LANES, 0.0);
    double t_seed = std::numeric_limits<double>::infinity();
    for (size_t l = 0; l < BOU_LANES; ++l) {
      const size_t src = base + std::min(l, n_lane - 1);
      models[l] = bou_model(keys[src]);
      const double ts = fpe::fpe_seed(models[l], keys[src].z_lo, keys[src].z_hi,
                                      g, t_max, q0[l], &absorbed[l]);
      if (ts < t_seed) t_seed = ts;
    }
    for (size_t l = 0; l < BOU_LANES; ++l) {
      const size_t src = base + std::min(l, n_lane - 1);
      fpe::fpe_seed(models[l], keys[src].z_lo, keys[src].z_hi, g, t_max, q0[l],
                    &absorbed[l], t_seed);
    }

    // A collapse in ANY lane forces the whole batch onto the time-varying path;
    // the batch marches on one clock and one factorisation buffer, so it cannot
    // be static for some lanes and not others.  In practice the boundary
    // parameters are shared across a trial's quadrature nodes, so this is
    // all-or-nothing anyway.
    bool stat = true;
    for (size_t l = 0; l < BOU_LANES; ++l)
      if (!models[l].static_op()) stat = false;

    const fpe::FPE_TimeSchedule sched =
      fpe::fpe_time_schedule(t_max - t_seed, gr.nt_for(t_max), gr.tgrade);

    const size_t SZ = static_cast<size_t>(M) * BOU_LANES;
    std::vector<double> Q(SZ), DIAG(SZ), SUB(SZ), SUP(SZ),
                        DINV(SZ), CPRIME(SZ);
    std::vector<double> opD(BOU_LANES, 0.0);
    // Second operator buffer, allocated only when the bounds actually move: the
    // Crank-Nicolson right-hand side needs the operator at the OLD time while
    // the left-hand side needs it at the new one.
    std::vector<double> DIAG2, SUB2, SUP2, opD2;
    double *Dg = DIAG.data(), *Sb = SUB.data(), *Sp = SUP.data(),
           *oD = opD.data();
    double *Dg2 = nullptr, *Sb2 = nullptr, *Sp2 = nullptr, *oD2 = nullptr;
    if (!stat) {
      DIAG2.resize(SZ); SUB2.resize(SZ); SUP2.resize(SZ);
      opD2.assign(BOU_LANES, 0.0);
      Dg2 = DIAG2.data(); Sb2 = SUB2.data(); Sp2 = SUP2.data();
      oD2 = opD2.data();
    }

    bou_build_op_lanes(models.data(), t_seed, g, Dg, Sb, Sp, oD);
    for (size_t l = 0; l < BOU_LANES; ++l)
      for (int i = 0; i < M; ++i)
        Q[static_cast<size_t>(i) * BOU_LANES + l] = q0[l][i];

    // Per-lane running state, mirroring fpe_solve exactly.
    std::vector<double> cdf_prev(BOU_LANES, 0.0), surv_prev(BOU_LANES, 1.0),
                        cdf_lo_prev(BOU_LANES, 0.0), cdf_lo_flux(BOU_LANES, 0.0),
                        cdf_up_prev(BOU_LANES, 0.0), gl_prev(BOU_LANES, 0.0);
    std::vector<Entry*> ent(BOU_LANES, nullptr);
    for (size_t l = 0; l < n_lane; ++l) ent[l] = out[base + l];

    for (size_t l = 0; l < n_lane; ++l) {
      double mass = 0.0;
      for (int i = 0; i < M; ++i)
        mass += g.dx[i] * Q[static_cast<size_t>(i) * BOU_LANES + l];
      const double cdf0 = std::min(1.0, std::max(0.0, 1.0 - mass));
      cdf_lo_prev[l] = std::min(cdf0, std::max(0.0, absorbed[l]));
      cdf_lo_flux[l] = cdf_lo_prev[l];
      cdf_up_prev[l] = cdf0 - cdf_lo_prev[l];
      cdf_prev[l] = cdf0;
      surv_prev[l] = std::min(1.0, std::max(0.0, mass));

      double gu = oD[l] * (g.cA * Q[static_cast<size_t>(M - 1) * BOU_LANES + l] -
                           g.cB * Q[static_cast<size_t>(M - 2) * BOU_LANES + l]);
      double glo = oD[l] * (g.cL_A * Q[l] - g.cL_B * Q[BOU_LANES + l]);
      if (!(gu > 0.0)) gu = 0.0;
      if (!(glo > 0.0)) glo = 0.0;
      gl_prev[l] = glo;

      Entry& e = *ent[l];
      e.key = keys[base + l]; e.t_max = t_max;
      e.t.clear(); e.log_pdf_lo.clear(); e.log_pdf_up.clear();
      e.cdf_lo.clear(); e.cdf_up.clear();
      e.t.push_back(t_seed);
      e.log_pdf_up.push_back(gu > 0.0 ? std::log(gu) : R_NegInf);
      e.log_pdf_lo.push_back(glo > 0.0 ? std::log(glo) : R_NegInf);
      e.cdf_lo.push_back(cdf_lo_prev[l]);
      e.cdf_up.push_back(cdf0 - cdf_lo_prev[l]);
    }

    // Factorise (I - c*L) for the current c, interleaved across lanes.  With
    // fixed bounds this runs once per schedule block; with collapsing bounds the
    // operator changes every step, so it runs once per step and there is nothing
    // left to amortise.
    auto factor = [&](double c, const double* dg, const double* sb,
                      const double* sp) {
#if defined(__AVX512F__) || defined(__AVX2__)
      // The Thomas dependency is along the spatial index, not across lanes.
      // Keep the lane-independent recurrence in registers and perform one
      // vector pass over all BOU_LANES rather than a scalar pass per lane.
      using Simd = BOUSimd<BOU_LANES>;
      using Vec = typename Simd::Vec;
      const Vec one = Simd::set1(1.0);
      const Vec cv = Simd::set1(c);
      const Vec neg_cv = Simd::set1(-c);
      const Vec d0 = Simd::sub(one,
                               Simd::mul(cv, Simd::load(dg)));
      Simd::store(DINV.data(), Simd::div(one, d0));
      for (int i = 1; i < M; ++i) {
        const size_t o = static_cast<size_t>(i) * BOU_LANES;
        const size_t p = static_cast<size_t>(i - 1) * BOU_LANES;
        const Vec cp = Simd::mul(
            Simd::mul(neg_cv, Simd::load(&sp[p])),
            Simd::load(&DINV[p]));
        Simd::store(&CPRIME[p], cp);
        const Vec d = Simd::sub(
            Simd::sub(one, Simd::mul(cv, Simd::load(&dg[o]))),
            Simd::mul(Simd::mul(neg_cv, Simd::load(&sb[o])), cp));
        Simd::store(&DINV[o], Simd::div(one, d));
      }
#else
      for (size_t l = 0; l < BOU_LANES; ++l) DINV[l] = 1.0 / (1.0 - c * dg[l]);
      for (int i = 1; i < M; ++i) {
        const size_t o = static_cast<size_t>(i) * BOU_LANES;
        const size_t p = static_cast<size_t>(i - 1) * BOU_LANES;
        for (size_t l = 0; l < BOU_LANES; ++l) {
          CPRIME[p + l] = -c * sp[p + l] * DINV[p + l];
          const double d = (1.0 - c * dg[o + l]) -
                           (-c * sb[o + l]) * CPRIME[p + l];
          DINV[o + l] = 1.0 / ((d != 0.0) ? d : DBL_MIN);
        }
      }
#endif
    };

    int n_steps = 0;
    for (size_t b = 0; b < sched.steps.size(); ++b) n_steps += sched.steps[b];
    const int n_rann = std::min(2, n_steps);

    double t = t_seed;
    int done = 0;
    std::vector<double> mass(BOU_LANES, 0.0);

    for (size_t b = 0; b < sched.dt.size(); ++b) {
      const double dtb = sched.dt[b];
      // Rannacher half-step (c = dt/2) and CN full step (c = dt/2) share c, so
      // one factorisation serves the whole block -- the same identity fpe_solve
      // relies on.  Only while the operator itself is constant, though.
      if (stat) factor(0.5 * dtb, Dg, Sb, Sp);
      for (int k = 0; k < sched.steps[b]; ++k, ++done) {
        const int n_sub = (done < n_rann) ? 2 : 1;
        for (int sub = 0; sub < n_sub; ++sub) {
          const bool be = (done < n_rann);
          const double step = be ? 0.5 * dtb : dtb;
          // The left-hand-side coefficient: c = step for a backward-Euler half
          // step, c = step/2 for Crank-Nicolson.  Both equal dtb/2, which is why
          // the fixed-bound path can share one factorisation.
          const double c = be ? step : 0.5 * step;

          // Advance the operator to the new time BEFORE the solve, then swap so
          // Dg/Sb/Sp are the new (left-hand-side) operator and Dg2/Sb2/Sp2 the
          // old one the Crank-Nicolson right-hand side still needs.
          if (!stat) {
            bou_build_op_lanes(models.data(), t + step, g, Dg2, Sb2, Sp2, oD2);
            std::swap(Dg, Dg2); std::swap(Sb, Sb2); std::swap(Sp, Sp2);
            std::swap(oD, oD2);
            factor(c, Dg, Sb, Sp);
          }
          const double* odg = stat ? Dg : Dg2;   // operator at the OLD time
          const double* osb = stat ? Sb : Sb2;
          const double* osp = stat ? Sp : Sp2;

          // Thomas, with the Crank-Nicolson right-hand side FUSED into the
          // forward sweep and both recurrences carried in lane-local arrays.
          //
          // Two separate costs are removed here.  The right-hand side used to
          // be a whole extra pass over Q and the operator, materialised into
          // RHS; the three q values a row needs (before the step) fit in a
          // rolling three-register window instead, and q_{i+1} is read before
          // q_i is overwritten so updating in place is safe.  And the sweeps
          // used to read Q[i-1] / Q[i+1] back one iteration after storing
          // them, which puts a store-to-load forward (~6 cycles) on top of the
          // 8-cycle multiply-add chain -- roughly doubling the cost of a sweep
          // that is already latency-bound.  Same arithmetic, same order.
          //
          // The right-hand side reads the operator at the OLD time
          // (odg/osb/osp) while the elimination reads the factorisation of the
          // new one (Sb/DINV/CPRIME), so fusing is valid for a moving boundary
          // as well as a fixed one -- they are different arrays.
          double y[BOU_LANES], xb[BOU_LANES], ms[BOU_LANES];
          if (be) {
            // Backward Euler: the right-hand side IS q, so there is nothing to
            // build -- read Q[i] and overwrite it in the same iteration.
            for (size_t l = 0; l < BOU_LANES; ++l) {
              y[l] = Q[l] * DINV[l];
              Q[l] = y[l];
            }
            for (int i = 1; i < M; ++i) {
              const size_t o = static_cast<size_t>(i) * BOU_LANES;
              for (size_t l = 0; l < BOU_LANES; ++l) {
                y[l] = (Q[o + l] + c * Sb[o + l] * y[l]) * DINV[o + l];
                Q[o + l] = y[l];
              }
            }
          } else {
            const double ch = 0.5 * step;
            double qm1[BOU_LANES], qcur[BOU_LANES], qnext[BOU_LANES];
            for (size_t l = 0; l < BOU_LANES; ++l) {
              qcur[l] = Q[l];
              qnext[l] = Q[BOU_LANES + l];
              const double r = qcur[l] + ch * odg[l] * qcur[l]
                                       + ch * osp[l] * qnext[l];
              y[l] = r * DINV[l];
              Q[l] = y[l];
            }
            for (int i = 1; i < M - 1; ++i) {
              const size_t o = static_cast<size_t>(i) * BOU_LANES;
              for (size_t l = 0; l < BOU_LANES; ++l) {
                qm1[l] = qcur[l];
                qcur[l] = qnext[l];
                qnext[l] = Q[o + BOU_LANES + l];
                const double r = qcur[l] + ch * odg[o + l] * qcur[l]
                                         + ch * osb[o + l] * qm1[l]
                                         + ch * osp[o + l] * qnext[l];
                y[l] = (r + c * Sb[o + l] * y[l]) * DINV[o + l];
                Q[o + l] = y[l];
              }
            }
            if (M > 1) {
              const size_t o = static_cast<size_t>(M - 1) * BOU_LANES;
              for (size_t l = 0; l < BOU_LANES; ++l) {
                qm1[l] = qcur[l];
                qcur[l] = qnext[l];
                const double r = qcur[l] + ch * odg[o + l] * qcur[l]
                                         + ch * osb[o + l] * qm1[l];
                y[l] = (r + c * Sb[o + l] * y[l]) * DINV[o + l];
                Q[o + l] = y[l];
              }
            }
          }
          for (size_t l = 0; l < BOU_LANES; ++l) {
            xb[l] = y[l];                       // == Q[(M-1)*BOU_LANES + l]
            ms[l] = g.dx[M - 1] * xb[l];
          }
          for (int i = M - 2; i >= 0; --i) {
            const size_t o = static_cast<size_t>(i) * BOU_LANES;
            for (size_t l = 0; l < BOU_LANES; ++l) {
              xb[l] = Q[o + l] - CPRIME[o + l] * xb[l];
              Q[o + l] = xb[l];
              ms[l] += g.dx[i] * xb[l];
            }
          }
          for (size_t l = 0; l < BOU_LANES; ++l) mass[l] = ms[l];

          t += step;
          for (size_t l = 0; l < n_lane; ++l) {
            double cdf_mass = std::min(1.0, std::max(0.0, 1.0 - mass[l]));
            if (cdf_mass < cdf_prev[l]) cdf_mass = cdf_prev[l];
            cdf_prev[l] = cdf_mass;
            double s_keep = std::min(1.0, std::max(0.0, mass[l]));
            if (s_keep > surv_prev[l]) s_keep = surv_prev[l];
            surv_prev[l] = s_keep;

            double gu = oD[l] *
              (g.cA * Q[static_cast<size_t>(M - 1) * BOU_LANES + l] -
               g.cB * Q[static_cast<size_t>(M - 2) * BOU_LANES + l]);
            double glo = oD[l] * (g.cL_A * Q[l] - g.cL_B * Q[BOU_LANES + l]);
            if (!(gu > 0.0)) gu = 0.0;
            if (!(glo > 0.0)) glo = 0.0;

            cdf_lo_flux[l] += 0.5 * step * (gl_prev[l] + glo);
            gl_prev[l] = glo;
            // Same two-sided window fpe_solve uses; see the comment there for
            // why capping at cdf_mass alone leaves the UPPER cdf non-monotone
            // on a kinked collapse.
            const double clo_hi = cdf_mass - cdf_up_prev[l];
            double clo = cdf_lo_flux[l];
            if (clo > clo_hi) clo = clo_hi;
            if (clo < cdf_lo_prev[l]) clo = cdf_lo_prev[l];
            cdf_lo_prev[l] = clo;
            cdf_up_prev[l] = cdf_mass - clo;

            Entry& e = *ent[l];
            e.t.push_back(t);
            e.log_pdf_up.push_back(gu > 0.0 ? std::log(gu) : R_NegInf);
            e.log_pdf_lo.push_back(glo > 0.0 ? std::log(glo) : R_NegInf);
            e.cdf_lo.push_back(clo);
            e.cdf_up.push_back(cdf_mass - clo);
          }
        }
      }
    }
    for (size_t l = 0; l < n_lane; ++l) ent[l]->grid_idx.build(ent[l]->t);
  }
}

inline Entry& bou_cache_get(SolveCache& C, const Key& p, double t_need) {
  // Always solve to the batch horizon, never to this row's need -- see
  // SolveCache::t_horizon.  The max() keeps a caller that forgot to set it
  // correct, merely slow.
  const double t_max = std::max(C.t_horizon, std::max(t_need, 1e-3));
  Entry* hit = C.find(p);
  if (hit != nullptr && hit->t_max >= t_need) return *hit;
  if (hit != nullptr) { bou_solve(p, t_max, C.grid, *hit); return *hit; }
  Entry& en = C.push(p);
  bou_solve(p, t_max, C.grid, en);
  return en;
}

// ---------------------------------------------------------------------------
// Interpolation.  Densities in log space, cdfs in the natural scale (they are
// about to be mixed linearly across quadrature nodes, and a defective cdf near
// zero is not the quantity whose relative accuracy matters).
// ---------------------------------------------------------------------------
struct Interp { double log_pdf_lo, log_pdf_up, cdf_lo, cdf_up; };

inline Interp bou_interp(const Entry& en, double t) {
  Interp o{R_NegInf, R_NegInf, 0.0, 0.0};
  const size_t n = en.t.size();
  if (n == 0) return o;
  if (t <= en.t.front()) {
    // Before the seed time the density is ~0, but the cdf is NOT: it carries the
    // mass absorbed during the analytic warm-up.  Same split fpe::grid_lookup
    // makes with zero_below_grid.
    o.cdf_lo = en.cdf_lo.front(); o.cdf_up = en.cdf_up.front();
    return o;
  }
  if (t >= en.t.back()) {
    o.log_pdf_lo = en.log_pdf_lo.back(); o.log_pdf_up = en.log_pdf_up.back();
    o.cdf_lo = en.cdf_lo.back(); o.cdf_up = en.cdf_up.back();
    return o;
  }
  const size_t j = en.grid_idx.find_index(en.t, t);
  const double t0 = en.t[j - 1], t1 = en.t[j];
  const double w = (t1 > t0) ? (t - t0) / (t1 - t0) : 0.0;

  // A floored neighbour carries no information -- take the live one rather than
  // interpolating toward -Inf (the rule fperace::interp_log uses).
  auto lin_log = [&](const std::vector<double>& lv) {
    const double a = lv[j - 1], b = lv[j];
    if (!R_FINITE(a)) return b;
    if (!R_FINITE(b)) return a;
    return a + w * (b - a);
  };
  o.log_pdf_lo = lin_log(en.log_pdf_lo);
  o.log_pdf_up = lin_log(en.log_pdf_up);
  o.cdf_lo = en.cdf_lo[j - 1] + w * (en.cdf_lo[j] - en.cdf_lo[j - 1]);
  o.cdf_up = en.cdf_up[j - 1] + w * (en.cdf_up[j] - en.cdf_up[j - 1]);
  return o;
}

// ---------------------------------------------------------------------------
// Across-trial variability.
//
// Returns the mixed DEFECTIVE density and cdf of BOTH responses at one time,
// integrating over drift (normal, sd sv), start point (uniform, width sz) and
// non-decision time (uniform, width st0).  Everything is mixed in the NATURAL
// scale, because a mixture is a sum -- taking logs is the caller's last step.
//
// The t0 integral is the outer loop: it only shifts the query time, so it never
// triggers a solve and costs one interpolation per node.  sv and sz are the
// expensive ones and are the inner loops, so their solves are shared across t0
// nodes through the cache.
// ---------------------------------------------------------------------------
struct BouMix { double d_lo, d_up, F_lo, F_up; };

inline BouMix bou_mix(SolveCache& C, double t_dec,
                      double v, double sv, double a, double z, double sz,
                      double s, double beta, bool anchor_at_z, double anchor_fix,
                      int bkind, double binf, double tau, double pw,
                      bool want_cdf) {
  BouMix acc{0.0, 0.0, 0.0, 0.0};
  const Grid& gr = C.grid;

  // Collect the quadrature's cache misses and solve them as one lane batch.
  // Solving them one at a time would leave the lanes empty, which is the whole
  // point of batching here -- the nodes share a mesh, horizon and clock.
  std::vector<Key> miss;
  std::vector<Entry*> miss_e;

  const bool do_sv = (sv > 0.0) && R_FINITE(sv);
  const bool do_sz = (sz > 0.0) && R_FINITE(sz);
  const int nv = do_sv ? std::max(1, gr.n_sv) : 1;
  // A midpoint-anchored operator is independent of the start point.  Seed
  // the entire uniform interval once instead of solving one point-start
  // operator per Gauss-Legendre node.
  const int nz = (anchor_at_z && do_sz) ? std::max(1, gr.n_sz) : 1;

  // gh_rule() rejects n < 2, so only build a rule when it is actually used.
  static const GHRule gh_dummy{};
  static const GLRule gl_dummy{};
  const GHRule& ghr = do_sv ? gh_rule(nv) : gh_dummy;
  const GLRule& glr = do_sz ? gl_get_rule(nz) : gl_dummy;

  for (int iv = 0; iv < nv; ++iv) {
    // The rule integrates exp(-x^2), so the standard-normal node is sqrt(2)*x
    // and the weight is divided by sqrt(pi) -- gh_standard_normal_weight does
    // the second half, the sqrt2 here does the first.
    const double wv = do_sv ? gh_standard_normal_weight(ghr, iv) : 1.0;
    const double vi = do_sv ? (v + sv * M_SQRT2 * ghr.x[iv]) : v;
    if (!(wv > 0.0)) continue;

    for (int iz = 0; iz < nz; ++iz) {
      double wz = 1.0, zi = z, z_lo = z, z_hi = z;
      if (anchor_at_z) {
        // Uniform on [z - sz/2, z + sz/2]; GL nodes on [-1,1] map by sz/2, and
        // the uniform density cancels against the interval width, leaving w/2.
        wz = do_sz ? 0.5 * glr.w[iz] : 1.0;
        zi = do_sz ? (z + 0.5 * sz * glr.x[iz]) : z;
        z_lo = z_hi = zi;
      } else if (do_sz) {
        // fpe_seed integrates this interval by exact cell overlap at t = 0.
        z_lo = z - 0.5 * sz;
        z_hi = z + 0.5 * sz;
      }
      if (!(wz > 0.0)) continue;
      if (anchor_at_z) {
        if (!(zi > 0.0) || !(zi < a)) continue;
      } else if (!(z_hi > 0.0) || !(z_lo < a)) {
        continue;
      }

      Key k;
      const double anc = anchor_at_z ? zi : anchor_fix;
      if (!bou_key(vi, beta, a, z_lo, z_hi, anc, s, bkind, binf, tau, pw, k))
        continue;

      // Two passes over the same node list: gather misses, solve them batched,
      // then read every node out of the cache.
      if (C.find(k) == nullptr) {
        bool queued = false;
        for (const Key& q : miss) if (q == k) { queued = true; break; }
        if (!queued) { miss.push_back(k); miss_e.push_back(nullptr); }
      }
    }
  }

  if (!miss.empty()) {
    const double t_max = std::max(C.t_horizon, std::max(t_dec, 1e-3));
    // Grow the entry store ONCE before taking any pointer into it: pushing one
    // at a time would reallocate mid-loop and leave every earlier pointer
    // dangling.
    C.reserve_for(miss.size());
    for (size_t i = 0; i < miss.size(); ++i) miss_e[i] = &C.push(miss[i]);
    bou_solve_batch(miss, t_max, gr, miss_e);
  }

  for (int iv = 0; iv < nv; ++iv) {
    const double wv = do_sv ? gh_standard_normal_weight(ghr, iv) : 1.0;
    const double vi = do_sv ? (v + sv * M_SQRT2 * ghr.x[iv]) : v;
    if (!(wv > 0.0)) continue;
    for (int iz = 0; iz < nz; ++iz) {
      double wz = 1.0, zi = z, z_lo = z, z_hi = z;
      if (anchor_at_z) {
        wz = do_sz ? 0.5 * glr.w[iz] : 1.0;
        zi = do_sz ? (z + 0.5 * sz * glr.x[iz]) : z;
        z_lo = z_hi = zi;
      } else if (do_sz) {
        z_lo = z - 0.5 * sz;
        z_hi = z + 0.5 * sz;
      }
      if (!(wz > 0.0)) continue;
      if (anchor_at_z) {
        if (!(zi > 0.0) || !(zi < a)) continue;
      } else if (!(z_hi > 0.0) || !(z_lo < a)) {
        continue;
      }

      Key k;
      const double anc = anchor_at_z ? zi : anchor_fix;
      if (!bou_key(vi, beta, a, z_lo, z_hi, anc, s, bkind, binf, tau, pw, k))
        continue;

      Entry* enp = C.find(k);
      if (enp == nullptr) continue;
      const Entry& en = *enp;
      const Interp it = bou_interp(en, t_dec);
      const double w = wv * wz;
      if (R_FINITE(it.log_pdf_lo)) acc.d_lo += w * std::exp(it.log_pdf_lo);
      if (R_FINITE(it.log_pdf_up)) acc.d_up += w * std::exp(it.log_pdf_up);
      if (want_cdf) { acc.F_lo += w * it.cdf_lo; acc.F_up += w * it.cdf_up; }
    }
  }
  return acc;
}

// Uniform t0 on [t0 - st0/2, t0 + st0/2], applied by shifting the decision time.
// For the DENSITY this is a convolution; for the CDF it is the same average, so
// one routine serves both.
inline BouMix bou_mix_t0(SolveCache& C, double rt,
                         double v, double sv, double a, double z, double sz,
                         double s, double beta, bool anchor_at_z, double anchor_fix,
                         double t0, double st0,
                         int bkind, double binf, double tau, double pw,
                         bool want_cdf) {
  const bool do_st0 = (st0 > 0.0) && R_FINITE(st0);
  if (!do_st0) {
    const double td = rt - t0;
    if (!(td > 0.0)) return BouMix{0.0, 0.0, 0.0, 0.0};
    return bou_mix(C, td, v, sv, a, z, sz, s, beta, anchor_at_z, anchor_fix,
                   bkind, binf, tau, pw, want_cdf);
  }
  const int n = std::max(1, C.grid.n_st0);
  const GLRule& glr = gl_get_rule(n);
  BouMix acc{0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < n; ++i) {
    // t0 ~ Uniform(t0, t0 + st0) -- the LOWER-EDGE convention, which is what the
    // package's DDM uses.  Verified against dDDM by reconstructing an st0 > 0
    // density from point-t0 densities: the lower-edge average reproduces it to
    // 2.2e-3 while the centred average [t0-st0/2, t0+st0/2] is off by 3.4e-1.
    // GL nodes live on [-1,1], so (1 + x)/2 maps them onto [0,1].
    const double t0i = t0 + st0 * 0.5 * (1.0 + glr.x[i]);
    const double w = 0.5 * glr.w[i];
    const double td = rt - t0i;
    if (!(td > 0.0)) continue;
    const BouMix m = bou_mix(C, td, v, sv, a, z, sz, s, beta, anchor_at_z,
                             anchor_fix, bkind, binf, tau, pw, want_cdf);
    acc.d_lo += w * m.d_lo; acc.d_up += w * m.d_up;
    if (want_cdf) { acc.F_lo += w * m.F_lo; acc.F_up += w * m.F_up; }
  }
  return acc;
}

} // namespace fpebou

#endif // fpe_bou_h
