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
#if defined(__x86_64__) || defined(_M_X64) || defined(__AVX2__)
#include <immintrin.h>
#endif
#include "fpe_models.h"

namespace fperace {

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

  // Boundary.  For FPE_BND_FIXED the shape fields are forced to zero by
  // rou_key(), so every fixed-bound row compares equal on them and a model that
  // never collapses behaves exactly as it did before collapse existed.
  int bkind = fpe::FPE_BND_FIXED;
  double binf = 0.0, tau = 0.0, pw = 0.0;

  bool operator==(const Key& o) const {
    return v == o.v && k == o.k && b == o.b && A == o.A &&
           bkind == o.bkind && binf == o.binf && tau == o.tau && pw == o.pw;
  }
  bool finite() const {
    return std::isfinite(v) && std::isfinite(k) && std::isfinite(b) &&
           std::isfinite(A) && std::isfinite(binf) && std::isfinite(tau) &&
           std::isfinite(pw);
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

// Fixed-boundary overload -- the common case, and every pre-collapse call site.
inline bool rou_key(double v, double k, double B, double A, double s, Key& out) {
  static const BndSpec fixed;
  return rou_key(v, k, B, A, s, fixed, out);
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
};

struct SolveCache {
  FPE_Grid grid;
  // Boundary form for this model instance, set once when the adapter is built.
  // It lives here rather than on ContextForRaceModels because it is a property
  // of the solve, and because every path that can reach a solve already holds
  // the cache.
  int bnd_kind = fpe::FPE_BND_FIXED;
  std::vector<Entry> e;
  std::vector<int> row_group;   // scratch: row -> index into e, or -1
  int prepared_n_rows = -1;

  // Cleared once per particle.  Keys are exact, so a stale entry could never be
  // returned for the wrong parameters; the clear exists to bound memory, since
  // a run visits thousands of particles.
  void new_particle() {
    e.clear();
    row_group.clear();
    prepared_n_rows = -1;
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
  m.sigma = 1.0;                                  // s already divided out
  m.bnd.set_kind(p.bkind, p.b,
                 (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                 p.tau, p.pw, false);
  m.xlo = fpe::fpe_x_lo_ou(0.0, p.v, p.k, 1.0, t_max);

  const fpe::FPE_Result r =
      fpe::fpe_run(m, 0.0, p.A, t_max, gr.nx, gr.nt_for(t_max), gr.grade,
                   gr.tgrade);

  const size_t n = r.t.size();
  out.key = p;
  out.t_max = t_max;
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
inline int cache_get(SolveCache& C, const Key& p, double t_need) {
  if (!(t_need > 0.0)) t_need = 1e-3;
  for (size_t i = 0; i < C.e.size(); ++i) {
    if (!(C.e[i].key == p)) continue;
    if (C.e[i].t_max >= t_need) return static_cast<int>(i);
    // Extend in place.  dt is pinned by FPE_Grid, so the values on the shared
    // part of the horizon move only at the discretisation-error level.
    rou_solve(p, t_need, C.grid, C.e[i]);
    return static_cast<int>(i);
  }
  C.e.push_back(Entry());
  rou_solve(p, t_need, C.grid, C.e.back());
  return static_cast<int>(C.e.size()) - 1;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Multi-solve 4-lane interleaved Thomas solver for fixed-boundary OU models.
// Interleaves up to 4 solves in a Struct-of-Arrays (SoA) layout so the serial
// recurrence dependency chain in Thomas forward sweep and back-substitution
// executes 4 independent SIMD lanes concurrently, filling the 12-cycle latency
// shadow of CPU execution ports.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline std::vector<fpe::FPE_Result> fpe_solve_batch_ou(
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
#if defined(__AVX2__)
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
#if defined(__AVX2__)
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
#if defined(__AVX2__)
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

// Batched cache lookup: solves missing keys in 4-lane SIMD batches
inline void cache_get_batch(SolveCache& C, const std::vector<Key>& keys,
                            const std::vector<double>& horizons,
                            std::vector<int>& out_cache_indices) {
  const size_t num_keys = keys.size();
  out_cache_indices.assign(num_keys, -1);
  if (num_keys == 0) return;

  std::vector<size_t> missing_keys;
  for (size_t i = 0; i < num_keys; ++i) {
    const Key& p = keys[i];
    double t_need = horizons[i];
    if (!(t_need > 0.0)) t_need = 1e-3;

    int found = -1;
    for (size_t j = 0; j < C.e.size(); ++j) {
      if (C.e[j].key == p && C.e[j].t_max >= t_need) {
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

  // Process missing keys in 4-lane batches
  for (size_t b = 0; b < missing_keys.size(); b += 4) {
    size_t chunk_size = std::min(static_cast<size_t>(4), missing_keys.size() - b);
    std::vector<fpe::FPE_ModelOU> models(chunk_size);
    std::vector<std::vector<double>> q0_vec(chunk_size);
    std::vector<double> t0_vec(chunk_size);
    std::vector<double> t_max_vec(chunk_size);

    fpe::FPE_Mesh mesh;
    mesh.build(C.grid.nx, C.grid.grade);

    bool all_fixed = true;
    for (size_t l = 0; l < chunk_size; ++l) {
      const size_t k_idx = missing_keys[b + l];
      const Key& p = keys[k_idx];
      if (p.bkind != fpe::FPE_BND_FIXED) all_fixed = false;

      models[l].v = p.v;
      models[l].lambda = p.k;
      models[l].sigma = 1.0;
      models[l].bnd.set_kind(p.bkind, p.b,
                             (p.bkind == fpe::FPE_BND_FIXED) ? p.b : p.binf,
                             p.tau, p.pw, false);
      models[l].xlo = fpe::fpe_x_lo_ou(0.0, p.v, p.k, 1.0, horizons[k_idx]);

      double t0 = fpe::fpe_seed(models[l], 0.0, p.A, mesh, horizons[k_idx], q0_vec[l]);
      t0_vec[l] = t0;
      t_max_vec[l] = horizons[k_idx];
    }

    if (all_fixed && chunk_size > 1) {
      double max_horiz = 0.0;
      for (size_t l = 0; l < chunk_size; ++l) {
        if (t_max_vec[l] > max_horiz) max_horiz = t_max_vec[l];
      }
      int nt = C.grid.nt_for(max_horiz);
      std::vector<fpe::FPE_Result> batch_res =
          fperace::fpe_solve_batch_ou(models, q0_vec, t0_vec, t_max_vec, mesh, nt, C.grid.tgrade);

      for (size_t l = 0; l < chunk_size; ++l) {
        const size_t k_idx = missing_keys[b + l];
        const Key& p = keys[k_idx];
        const fpe::FPE_Result& r = batch_res[l];
        const size_t n = r.t.size();

        Entry entry;
        entry.key = p;
        entry.t_max = horizons[k_idx];
        entry.t = r.t;
        entry.log_pdf.resize(n);
        entry.log_S.resize(n);
        for (size_t j = 0; j < n; ++j) {
          entry.log_pdf[j] = safe_log(r.pdf[j]);
          entry.log_S[j] = safe_log(r.surv[j]);
        }
        entry.grid_idx.build(entry.t);

        int existing = -1;
        for (size_t j = 0; j < C.e.size(); ++j) {
          if (C.e[j].key == p) { existing = static_cast<int>(j); break; }
        }
        if (existing >= 0) {
          C.e[existing] = std::move(entry);
          out_cache_indices[k_idx] = existing;
        } else {
          C.e.push_back(std::move(entry));
          out_cache_indices[k_idx] = static_cast<int>(C.e.size()) - 1;
        }
      }
    } else {
      for (size_t l = 0; l < chunk_size; ++l) {
        const size_t k_idx = missing_keys[b + l];
        out_cache_indices[k_idx] = cache_get(C, keys[k_idx], horizons[k_idx]);
      }
    }
  }
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

inline double entry_log_pdf(const Entry& en, double tt) {
  return interp_log(en, tt, true);
}
inline double entry_log_S(const Entry& en, double tt) {
  return interp_log(en, tt, false);
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
  const double var = (k > 1e-10) ? (s * s * m2 / (2.0 * k)) : (s * s * dt);
  const double sd = std::sqrt(std::max(var, 0.0));
  const double inv_2var_bb = 2.0 / (s * s * dt);     // local BB uses s^2 dt

  double t = 0.0;
  while (t < t_max) {
    const double X1 = X * phi + v * drift_gain + sd * rnorm();
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
