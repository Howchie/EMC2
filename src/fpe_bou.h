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

      Entry& en = bou_cache_get(C, k, std::max(t_dec, 1e-3));
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
