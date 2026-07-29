#ifndef model_ROU_h
#define model_ROU_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ROU -- the racing Ornstein-Uhlenbeck (leaky accumulator) race kernels.
//
//   dX_i = (v_i - k_i X_i) dt + s_i dW_i,   X_i(0) ~ U(0, A_i),
//   accumulator i finishes when X_i reaches b_i = B_i + A_i.
//
// This is the leaky accumulator of Smith & Ratcliff (2004, appendix), i.e. an
// OU race.  It is NOT the leaky COMPETING accumulator: there is no lateral
// inhibition between accumulators and no rectification of the activation at
// zero.  The process is free to go negative; only upward boundary hits count.
//
// The five entry points below are the whole race contract (see utils.h).  All
// of them route through fperace::SolveCache, so the number of PDE solves is the
// number of distinct parameter tuples rather than the number of rows -- and a
// scalar cdf1() call from the censoring path is an interpolation into a solve
// the batch kernels already paid for.
//
// Column order per src/col_registry.h: v, k, B, A, t0, s.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "col_registry.h"
#include "fpe_race.h"

// ---------------------------------------------------------------------------
// Resolution knobs, read from R options once when the adapter is built.
//
// These MUST be honoured by the sampled likelihood, not only by the R-side
// dROU/pROU: without that there is no way to check that a fit has converged in
// (nx, dt), and no way to trade accuracy for speed on a real data set.  Read
// once per likelihood call, never per row.
// ---------------------------------------------------------------------------
inline void rou_configure_grid(fperace::FPE_Grid& g) {
  auto opt_num = [](const char* nm, double& dst, double lo) {
    SEXP v = Rf_GetOption1(Rf_install(nm));
    if (v == R_NilValue || Rf_length(v) < 1) return;
    const double x = Rf_asReal(v);
    if (R_finite(x) && x >= lo) dst = x;
  };
  double nx = g.nx, dt = g.dt_target, gr = g.grade, tg = g.tgrade;
  opt_num("emc2.fpe_nx", nx, 8.0);
  opt_num("emc2.fpe_dt", dt, 1e-6);
  opt_num("emc2.fpe_grade", gr, 1.0);
  opt_num("emc2.fpe_tgrade", tg, 1.0);
  g.nx = static_cast<int>(nx);
  g.dt_target = dt;
  g.grade = gr;
  g.tgrade = tg;
  // nt_max has to follow dt_target, or asking for a finer dt silently does
  // nothing on the long-horizon groups where it matters most.
  g.nt_max = std::max(g.nt_max, static_cast<int>(std::ceil(10.0 / dt)));
}

// ---------------------------------------------------------------------------
// Named parameters -> boundary form.  This is the ONLY place that knows the
// mapping; fpe::FPE_Boundary knows no parameter names and fperace::BndSpec
// carries only (kind, Binf, tau, pw).  Both accessors are gated on bnd_kind, so
// the optional columns are never dereferenced for a fixed-bound model -- where
// they do not exist and the pointers are null.
// ---------------------------------------------------------------------------
inline fperace::BndSpec rou_bnd_row(int kind, const double* const* cols, int i) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = cols[emc2col::rou::Binf][i];
  bs.tau  = cols[emc2col::rou::tau][i];
  bs.pw   = (kind == fpe::FPE_BND_WEIBULL) ? cols[emc2col::rou::pw][i] : 0.0;
  return bs;
}

inline fperace::BndSpec rou_bnd_par(int kind, const double* par) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = par[emc2col::rou::Binf];
  bs.tau  = par[emc2col::rou::tau];
  bs.pw   = (kind == fpe::FPE_BND_WEIBULL) ? par[emc2col::rou::pw] : 0.0;
  return bs;
}

inline fperace::SolveCache* rou_cache(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (ctx == nullptr) return nullptr;
  if (!ctx->fpe_cache) {
    ctx->fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_grid(ctx->fpe_cache->grid);
  }
  return ctx->fpe_cache.get();
}

// ---------------------------------------------------------------------------
// Group the masked rows by parameter tuple and make sure every group has a
// solve reaching its largest required time.  Writes C.row_group[i] = index into
// C.e, or -1 for rows that are excluded, invalid, or have tt <= 0.
//
// Two passes rather than one: the horizon of a group is only known once every
// row in it has been seen, and solving to a short horizon and then extending
// would do the expensive part twice.
// ---------------------------------------------------------------------------
inline void rou_prepare_rows(fperace::SolveCache& C, const double* rt,
                             const double* const* cols, int n_rows,
                             const int* mask, const int* isok) {
  const double* v_  = cols[emc2col::rou::v];
  const double* k_  = cols[emc2col::rou::k];
  const double* B_  = cols[emc2col::rou::B];
  const double* A_  = cols[emc2col::rou::A];
  const double* t0_ = cols[emc2col::rou::t0];
  const double* s_  = cols[emc2col::rou::s];

  C.row_group.assign(n_rows, -1);

  // Pass 1: collect distinct keys and the horizon each one needs.
  std::vector<fperace::Key> keys;
  std::vector<double> horizon;
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i] || !isok[i] || R_IsNA(v_[i])) continue;
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !emc2_isfinite(tt)) continue;
    fperace::Key p;
    if (!fperace::rou_key(v_[i], k_[i], B_[i], A_[i], s_[i],
                          rou_bnd_row(C.bnd_kind, cols, i), p)) continue;

    int g = -1;
    for (size_t j = 0; j < keys.size(); ++j) {
      if (keys[j] == p) { g = static_cast<int>(j); break; }
    }
    if (g < 0) {
      keys.push_back(p);
      horizon.push_back(tt);
      g = static_cast<int>(keys.size()) - 1;
    } else if (tt > horizon[g]) {
      horizon[g] = tt;
    }
    C.row_group[i] = g;
  }

  // Pass 2: solve, then rewrite the local group indices as cache indices.
  std::vector<int> to_cache(keys.size(), -1);
  for (size_t j = 0; j < keys.size(); ++j) {
    to_cache[j] = fperace::cache_get(C, keys[j], horizon[j]);
  }
  for (int i = 0; i < n_rows; ++i) {
    if (C.row_group[i] >= 0) C.row_group[i] = to_cache[C.row_group[i]];
  }
}

// ---------------------------------------------------------------------------
// Batch log-density at the winner rows.
// ---------------------------------------------------------------------------
inline void drou_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  rou_prepare_rows(*C, rt, cols, n_rows, mask, isok);

  // The solve is in the rescaled state Y = X/s, whose first-passage time is the
  // SAME random variable, so the density in t needs no Jacobian -- exactly as
  // the RDM divides its arguments by s and applies none.
  const double* t0_ = cols[emc2col::rou::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double lp = fperace::entry_log_pdf(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(lp, min_ll, floor_raw);
  }
}

// ---------------------------------------------------------------------------
// Batch log-SURVIVOR at the loser rows.  Note this is log S, not log CDF; the
// survivor comes from the solver's mass route (sum of positive cell masses), so
// it keeps relative accuracy far past the point where 1 - CDF has none.
// ---------------------------------------------------------------------------
inline void prou_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok,
                     double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  rou_prepare_rows(*C, rt, cols, n_rows, mask, isok);

  const double* t0_ = cols[emc2col::rou::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) {
      // Before t0, or an unusable row.  An accumulator that cannot have
      // finished has survived, exactly as the analytic kernels report.
      out[i] = 0.0;
      continue;
    }
    const double lS = fperace::entry_log_S(C->e[g], rt[i] - t0_[i]);
    out[i] = (lS >= 0.0) ? 0.0
                         : ((lS <= fperace::LOG_FLOOR)
                                ? raw_log_zero(min_ll, floor_raw)
                                : lS);
  }
}

// ---------------------------------------------------------------------------
// Log-survivor of a whole trial at one scalar t, for truncation normalisers.
// ---------------------------------------------------------------------------
inline void rou_logS_at_t(double t, const double* const* cols,
                          int n_rows_total, int n_lR, int /*n_par*/,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* ctx_, double* logS_out) {
  fperace::SolveCache* C = rou_cache(ctx_);
  const double* v_  = cols[emc2col::rou::v];
  const double* k_  = cols[emc2col::rou::k];
  const double* B_  = cols[emc2col::rou::B];
  const double* A_  = cols[emc2col::rou::A];
  const double* t0_ = cols[emc2col::rou::t0];
  const double* s_  = cols[emc2col::rou::s];

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int m = 0; m < n_lR && !bad; ++m) {
      const int r = start + m;
      if (r >= n_rows_total) break;
      if (!isok_all[r] || R_IsNA(v_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;                       // cannot have finished
      fperace::Key p;
      if (!fperace::rou_key(v_[r], k_[r], B_[r], A_[r], s_[r],
                            rou_bnd_row(C->bnd_kind, cols, r), p)) { bad = true; break; }
      const int g = fperace::cache_get(*C, p, tt);
      const double lS = fperace::entry_log_S(C->e[g], tt);
      if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
      if (lS < 0.0) logS += lS;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

// ---------------------------------------------------------------------------
// Scalar helpers, for the GSL integrands on the known-winner / active-nogo
// branches, and for the scalar CDF calls the censoring path makes at the
// LT / LC / UC / UT bounds.
//
// They are NOT free the way an analytic model's scalars are: a miss costs a
// whole march.  Two things keep that from mattering.  First, the parameter
// tuple is almost always one the batch kernels already solved, so the call is
// an interpolation.  Second, a genuine miss solves to a horizon comfortably
// beyond the requested time (rou_scalar_horizon), so a GSL sweep over an
// interval pays for at most one solve rather than one per integrand call.
// ---------------------------------------------------------------------------
inline double rou_scalar_horizon(double tt) {
  return (tt > 1.0) ? (2.0 * tt) : (tt + 1.0);
}

inline double drou_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::rou::v])) return 0.0;
  const double tt = t - par[emc2col::rou::t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) return 0.0;
  fperace::Key p;
  if (!fperace::rou_key(par[emc2col::rou::v], par[emc2col::rou::k],
                    par[emc2col::rou::B], par[emc2col::rou::A],
                    par[emc2col::rou::s],
                    rou_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, rou_scalar_horizon(tt));
  const double lp = fperace::entry_log_pdf(C->e[g], tt);
  return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
}

inline double prou_scalar(double t, const double* par, void* ctx_) {
  if (R_IsNA(par[emc2col::rou::v])) return 0.0;
  const double tt = t - par[emc2col::rou::t0];
  if (!(tt > 0.0)) return 0.0;
  fperace::SolveCache* C = rou_cache(ctx_);
  if (C == nullptr) return 0.0;
  fperace::Key p;
  if (!fperace::rou_key(par[emc2col::rou::v], par[emc2col::rou::k],
                    par[emc2col::rou::B], par[emc2col::rou::A],
                    par[emc2col::rou::s],
                    rou_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, rou_scalar_horizon(tt));
  const double lS = fperace::entry_log_S(C->e[g], tt);
  if (lS >= 0.0) return 0.0;
  if (lS <= fperace::LOG_FLOOR) return 1.0;
  return -std::expm1(lS);            // 1 - S, accurate for S near 1
}

#endif // model_ROU_h
