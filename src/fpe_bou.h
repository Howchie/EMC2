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
//   * sz (start, uniform)     changes the START, and because Smith & Ratcliff
//                             anchor the decay AT the start point the drift is
//                             xi - beta*(x - z) = (xi + beta*z) - beta*x, so z
//                             changes the OPERATOR too => new solve
//   * st0 (non-decision)      a convolution on the finished density => free
//
// So the cost is n_sv * n_sz solves per distinct parameter row, not one.  There
// is no way around it while the anchor is z: substituting y = x - z removes z
// from the drift but moves it into the domain [-z, a-z] instead, so the solve is
// still z-specific.  (An anchor that does not depend on z -- the midpoint -- WOULD
// collapse the whole uniform start into a single march, which is why
// FPE_ModelBoundedOU keeps `anchor` as a field.  Note that zero is not available:
// it is one of the response boundaries.)
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
struct Grid {
  int nx = 512;
  double dt_target = 5e-4;      // the binding constraint for this model
  int nt_min = 512;
  int nt_max = 16384;
  double grade = fpe::FPE_GRADE_BOUNDED;   // 1.0: uniform, measured best
  double tgrade = fpe::FPE_TGRADE;
  int n_sv = 7;                 // Gauss-Hermite nodes for drift variability
  int n_sz = 7;                 // Gauss-Legendre nodes for start-point variability
  int n_st0 = 7;                // Gauss-Legendre nodes for the t0 convolution

  int nt_for(double t_max) const {
    const int n = static_cast<int>(std::ceil(t_max / dt_target));
    return std::max(nt_min, std::min(nt_max, n));
  }
};

// ---------------------------------------------------------------------------
// Cache key.  One entry per (drift, leak, separation, start, anchor, boundary)
// tuple, compared BIT-EXACTLY: these values are read from the same ParamTable
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
  double v = 0.0, beta = 0.0, a = 1.0, z = 0.5, anchor = 0.5;
  int bkind = fpe::FPE_BND_FIXED;
  double binf = 0.0, tau = 0.0, pw = 0.0;

  bool operator==(const Key& o) const {
    return v == o.v && beta == o.beta && a == o.a && z == o.z &&
           anchor == o.anchor && bkind == o.bkind && binf == o.binf &&
           tau == o.tau && pw == o.pw;
  }
};

// Build a key from natural-scale parameters, scaling s out.  Returns false for
// a parameter set the solver cannot answer, so the caller floors the row rather
// than solving nonsense.
inline bool bou_key(double v, double beta, double a, double z, double anchor,
                    double s, int bkind, double binf, double tau, double pw,
                    Key& out) {
  if (!(s > 0.0) || !R_FINITE(s)) return false;
  if (!(a > 0.0) || !R_FINITE(a)) return false;
  if (!R_FINITE(v) || !R_FINITE(beta) || beta < 0.0) return false;
  if (!(z > 0.0) || !(z < a)) return false;

  const double inv = 1.0 / s;
  out.v = v * inv;
  out.beta = beta;                 // a rate: invariant under the state rescale
  out.a = a * inv;
  out.z = z * inv;
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

  // Cleared once per particle.  Keys are exact, so a stale entry could never be
  // returned for different parameters -- but the parameters DO change every
  // particle, so keeping them would only grow the scan.
  void new_particle() { n_entries = 0; }

  Entry* find(const Key& k) {
    for (size_t i = 0; i < n_entries; ++i)
      if (e[i].key == k) return &e[i];
    return nullptr;
  }
  Entry& push() {
    if (n_entries >= e.size()) e.resize(n_entries + 1);
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
inline void bou_solve(const Key& p, double t_max, const Grid& gr, Entry& out) {
  fpe::FPE_ModelBoundedOU m;
  m.v = p.v;
  m.beta = p.beta;
  m.anchor = p.anchor;
  m.sigma = 1.0;                 // s scaled out in bou_key
  m.xlo = 0.0;
  if (p.bkind == fpe::FPE_BND_FIXED) {
    m.bnd.set_kind(fpe::FPE_BND_FIXED, p.a, p.a, 0.0, 0.0, false);
  } else {
    m.bnd.set_kind(p.bkind, p.a, p.binf, p.tau, p.pw, false);
  }

  fpe::FPE_Mesh g;
  g.build(gr.nx, gr.grade, fpe::FPE_ModelBoundedOU::symmetric_mesh);
  std::vector<double> q0;
  double absorbed_lower = 0.0;
  // Point start: the sz spread is carried by the quadrature outside this
  // function, because with the decay anchored at z each start point needs its
  // own operator and so its own march.
  const double t0 = fpe::fpe_seed(m, p.z, p.z, g, t_max, q0, &absorbed_lower);
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
#else
constexpr size_t BOU_LANES = 4;
#endif

inline void bou_solve_batch(const std::vector<Key>& keys, double t_max,
                            const Grid& gr, std::vector<Entry*>& out) {
  const size_t nk = keys.size();
  if (nk == 0) return;

  fpe::FPE_Mesh g;
  g.build(gr.nx, gr.grade, fpe::FPE_ModelBoundedOU::symmetric_mesh);
  const int M = g.M;

  auto make_model = [&](const Key& p) {
    fpe::FPE_ModelBoundedOU m;
    m.v = p.v; m.beta = p.beta; m.anchor = p.anchor; m.sigma = 1.0; m.xlo = 0.0;
    if (p.bkind == fpe::FPE_BND_FIXED)
      m.bnd.set_kind(fpe::FPE_BND_FIXED, p.a, p.a, 0.0, 0.0, false);
    else
      m.bnd.set_kind(p.bkind, p.a, p.binf, p.tau, p.pw, false);
    return m;
  };

  for (size_t base = 0; base < nk; base += BOU_LANES) {
    const size_t n_lane = std::min(BOU_LANES, nk - base);

    // Pass 1: learn each lane's preferred seed time, then adopt the smallest.
    std::vector<fpe::FPE_ModelBoundedOU> models(BOU_LANES);
    std::vector<std::vector<double>> q0(BOU_LANES);
    std::vector<double> absorbed(BOU_LANES, 0.0);
    double t_seed = std::numeric_limits<double>::infinity();
    for (size_t l = 0; l < BOU_LANES; ++l) {
      const size_t src = base + std::min(l, n_lane - 1);
      models[l] = make_model(keys[src]);
      const double ts = fpe::fpe_seed(models[l], keys[src].z, keys[src].z, g,
                                      t_max, q0[l], &absorbed[l]);
      if (ts < t_seed) t_seed = ts;
    }
    for (size_t l = 0; l < BOU_LANES; ++l) {
      const size_t src = base + std::min(l, n_lane - 1);
      fpe::fpe_seed(models[l], keys[src].z, keys[src].z, g, t_max, q0[l],
                    &absorbed[l], t_seed);
    }

    const fpe::FPE_TimeSchedule sched =
      fpe::fpe_time_schedule(t_max - t_seed, gr.nt_for(t_max), gr.tgrade);

    const size_t SZ = static_cast<size_t>(M) * BOU_LANES;
    std::vector<double> Q(SZ), RHS(SZ), DIAG(SZ), SUB(SZ), SUP(SZ),
                        DINV(SZ), CPRIME(SZ);
    std::vector<double> opD(BOU_LANES, 0.0);

    for (size_t l = 0; l < BOU_LANES; ++l) {
      fpe::FPE_Op op;
      fpe::build_op(models[l], t_seed, g, op);
      opD[l] = op.D;
      for (int i = 0; i < M; ++i) {
        const size_t o = static_cast<size_t>(i) * BOU_LANES + l;
        Q[o] = q0[l][i]; DIAG[o] = op.diag[i];
        SUB[o] = op.sub[i]; SUP[o] = op.sup[i];
      }
    }

    // Per-lane running state, mirroring fpe_solve exactly.
    std::vector<double> cdf_prev(BOU_LANES, 0.0), surv_prev(BOU_LANES, 1.0),
                        cdf_lo_prev(BOU_LANES, 0.0), cdf_lo_flux(BOU_LANES, 0.0),
                        gl_prev(BOU_LANES, 0.0);
    std::vector<Entry*> ent(BOU_LANES, nullptr);
    for (size_t l = 0; l < n_lane; ++l) ent[l] = out[base + l];

    for (size_t l = 0; l < n_lane; ++l) {
      double mass = 0.0;
      for (int i = 0; i < M; ++i)
        mass += g.dx[i] * Q[static_cast<size_t>(i) * BOU_LANES + l];
      const double cdf0 = std::min(1.0, std::max(0.0, 1.0 - mass));
      cdf_lo_prev[l] = std::min(cdf0, std::max(0.0, absorbed[l]));
      cdf_lo_flux[l] = cdf_lo_prev[l];
      cdf_prev[l] = cdf0;
      surv_prev[l] = std::min(1.0, std::max(0.0, mass));

      double gu = opD[l] * (g.cA * Q[static_cast<size_t>(M - 1) * BOU_LANES + l] -
                            g.cB * Q[static_cast<size_t>(M - 2) * BOU_LANES + l]);
      double glo = opD[l] * (g.cL_A * Q[l] - g.cL_B * Q[BOU_LANES + l]);
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

    // Factorise (I - c*L) for the current c, interleaved across lanes.
    auto factor = [&](double c) {
      for (size_t l = 0; l < BOU_LANES; ++l) DINV[l] = 1.0 / (1.0 - c * DIAG[l]);
      for (int i = 1; i < M; ++i) {
        const size_t o = static_cast<size_t>(i) * BOU_LANES;
        const size_t p = static_cast<size_t>(i - 1) * BOU_LANES;
        for (size_t l = 0; l < BOU_LANES; ++l) {
          CPRIME[p + l] = -c * SUP[p + l] * DINV[p + l];
          const double d = (1.0 - c * DIAG[o + l]) -
                           (-c * SUB[o + l]) * CPRIME[p + l];
          DINV[o + l] = 1.0 / ((d != 0.0) ? d : DBL_MIN);
        }
      }
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
      // relies on.
      factor(0.5 * dtb);
      for (int k = 0; k < sched.steps[b]; ++k, ++done) {
        const int n_sub = (done < n_rann) ? 2 : 1;
        for (int sub = 0; sub < n_sub; ++sub) {
          const bool be = (done < n_rann);
          const double step = be ? 0.5 * dtb : dtb;

          if (be) {
            RHS = Q;
          } else {
            const double c = 0.5 * step;
            for (int i = 0; i < M; ++i) {
              const size_t o = static_cast<size_t>(i) * BOU_LANES;
              for (size_t l = 0; l < BOU_LANES; ++l) {
                double vv = Q[o + l] + c * DIAG[o + l] * Q[o + l];
                if (i > 0) vv += c * SUB[o + l] * Q[o - BOU_LANES + l];
                if (i < M - 1) vv += c * SUP[o + l] * Q[o + BOU_LANES + l];
                RHS[o + l] = vv;
              }
            }
          }

          // Thomas: forward substitution then back substitution, with the mass
          // accumulated on the way back so the sweep is single-pass.
          const double c = be ? step : 0.5 * step;
          for (size_t l = 0; l < BOU_LANES; ++l) Q[l] = RHS[l] * DINV[l];
          for (int i = 1; i < M; ++i) {
            const size_t o = static_cast<size_t>(i) * BOU_LANES;
            const size_t p = o - BOU_LANES;
            for (size_t l = 0; l < BOU_LANES; ++l)
              Q[o + l] = (RHS[o + l] + c * SUB[o + l] * Q[p + l]) * DINV[o + l];
          }
          for (size_t l = 0; l < BOU_LANES; ++l)
            mass[l] = g.dx[M - 1] * Q[static_cast<size_t>(M - 1) * BOU_LANES + l];
          for (int i = M - 2; i >= 0; --i) {
            const size_t o = static_cast<size_t>(i) * BOU_LANES;
            for (size_t l = 0; l < BOU_LANES; ++l) {
              Q[o + l] -= CPRIME[o + l] * Q[o + BOU_LANES + l];
              mass[l] += g.dx[i] * Q[o + l];
            }
          }

          t += step;
          for (size_t l = 0; l < n_lane; ++l) {
            double cdf_mass = std::min(1.0, std::max(0.0, 1.0 - mass[l]));
            if (cdf_mass < cdf_prev[l]) cdf_mass = cdf_prev[l];
            cdf_prev[l] = cdf_mass;
            double s_keep = std::min(1.0, std::max(0.0, mass[l]));
            if (s_keep > surv_prev[l]) s_keep = surv_prev[l];
            surv_prev[l] = s_keep;

            double gu = opD[l] *
              (g.cA * Q[static_cast<size_t>(M - 1) * BOU_LANES + l] -
               g.cB * Q[static_cast<size_t>(M - 2) * BOU_LANES + l]);
            double glo = opD[l] * (g.cL_A * Q[l] - g.cL_B * Q[BOU_LANES + l]);
            if (!(gu > 0.0)) gu = 0.0;
            if (!(glo > 0.0)) glo = 0.0;

            cdf_lo_flux[l] += 0.5 * step * (gl_prev[l] + glo);
            gl_prev[l] = glo;
            double clo = cdf_lo_flux[l];
            if (clo < cdf_lo_prev[l]) clo = cdf_lo_prev[l];
            if (clo > cdf_mass) clo = cdf_mass;
            cdf_lo_prev[l] = clo;

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
  Entry& en = C.push();
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
  const int nz = do_sz ? std::max(1, gr.n_sz) : 1;

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
      // Uniform on [z - sz/2, z + sz/2]; GL nodes on [-1,1] map by sz/2, and the
      // uniform density cancels against the interval width, leaving w/2.
      const double wz = do_sz ? 0.5 * glr.w[iz] : 1.0;
      const double zi = do_sz ? (z + 0.5 * sz * glr.x[iz]) : z;
      if (!(wz > 0.0)) continue;
      if (!(zi > 0.0) || !(zi < a)) continue;   // a start outside the barriers

      Key k;
      const double anc = anchor_at_z ? zi : anchor_fix;
      if (!bou_key(vi, beta, a, zi, anc, s, bkind, binf, tau, pw, k)) continue;

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
    for (size_t i = 0; i < miss.size(); ++i) miss_e[i] = &C.push();
    bou_solve_batch(miss, t_max, gr, miss_e);
  }

  for (int iv = 0; iv < nv; ++iv) {
    const double wv = do_sv ? gh_standard_normal_weight(ghr, iv) : 1.0;
    const double vi = do_sv ? (v + sv * M_SQRT2 * ghr.x[iv]) : v;
    if (!(wv > 0.0)) continue;
    for (int iz = 0; iz < nz; ++iz) {
      const double wz = do_sz ? 0.5 * glr.w[iz] : 1.0;
      const double zi = do_sz ? (z + 0.5 * sz * glr.x[iz]) : z;
      if (!(wz > 0.0)) continue;
      if (!(zi > 0.0) || !(zi < a)) continue;

      Key k;
      const double anc = anchor_at_z ? zi : anchor_fix;
      if (!bou_key(vi, beta, a, zi, anc, s, bkind, binf, tau, pw, k)) continue;

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
