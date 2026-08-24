#ifndef model_ROUp_h
#define model_ROUp_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ROUp -- Racing Ornstein-Uhlenbeck with Smith (1995) time-varying pulse drift
//
// Each accumulator has time-varying drift:
//   mu(t) = v_S * (1 - exp(-t / tau_S)) + v_T * (t / tau_T) * exp(-t / tau_T)
//   dX_i = (mu_i(t) - k_i X_i) dt + s_i dW_i,  X_i(0) ~ U(0, A_i),
//   accumulator i finishes when X_i reaches b_i = B_i + A_i.
//
// Column order per src/col_registry.h: v_S, v_T, tau_S, tau_T, k, B, A, t0, s.
// Optional collapse columns: Binf, tau, pw.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include "col_registry.h"
#include "fpe_race.h"

inline void roup_configure_grid(fperace::FPE_Grid& g) {
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
  g.nt_max = std::max(g.nt_max, static_cast<int>(std::ceil(10.0 / dt)));
}

inline void roup_configure_cache(fperace::SolveCache& cache) {
  roup_configure_grid(cache.grid);
  SEXP sparse = Rf_GetOption1(Rf_install("emc2.rou_sparse_output"));
  if (sparse != R_NilValue && Rf_length(sparse) > 0) {
    const int enabled = Rf_asLogical(sparse);
    if (enabled != NA_LOGICAL) cache.sparse_raw_output = enabled;
  }
}

inline fperace::BndSpec roup_bnd_row(int kind, const double* const* cols, int i) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = cols[emc2col::roup::Binf][i];
  bs.tau  = cols[emc2col::roup::tau][i];
  bs.pw   = (kind == fpe::FPE_BND_WEIBULL) ? cols[emc2col::roup::pw][i] : 0.0;
  return bs;
}

inline fperace::BndSpec roup_bnd_par(int kind, const double* par) {
  fperace::BndSpec bs;
  bs.kind = kind;
  if (kind == fpe::FPE_BND_FIXED) return bs;
  bs.Binf = par[emc2col::roup::Binf];
  bs.tau  = par[emc2col::roup::tau];
  bs.pw   = (kind == fpe::FPE_BND_WEIBULL) ? par[emc2col::roup::pw] : 0.0;
  return bs;
}

inline double roup_logaddexp(double a, double b) {
  if (a <= fperace::LOG_FLOOR) return b;
  if (b <= fperace::LOG_FLOOR) return a;
  const double hi = std::max(a, b);
  return hi + std::log1p(std::exp(std::min(a, b) - hi));
}

// Build the one- or two-entry representation used by local ROUp.  Inactive
// channels are canonicalised by roup_key(), so their tau values are ignored.
inline bool roup_local_keys(double v_S, double v_T, double tau_S, double tau_T,
                            double k, double B, double A, double s,
                            const fperace::BndSpec& bs,
                            fperace::Key& key_s, fperace::Key& key_t,
                            bool& has_s, bool& has_t) {
  if (!std::isfinite(v_S) || !std::isfinite(v_T) || v_S < 0.0 || v_T < 0.0)
    return false;
  has_s = v_S > fpe::ROUP_DRIFT_EPS;
  has_t = v_T > fpe::ROUP_DRIFT_EPS;
  if (!has_s && !has_t) return false;
  if (has_s && !fperace::roup_key(v_S, 0.0, tau_S, 0.0, k, B, A, s, bs, key_s))
    return false;
  if (has_t && !fperace::roup_key(0.0, v_T, 0.0, tau_T, k, B, A, s, bs, key_t))
    return false;
  if (!has_s) key_s = key_t;
  if (!has_t) key_t = key_s;
  return true;
}

inline void roup_local_eval(const fperace::SolveCache& C, int gs, int gt,
                            bool has_s, bool has_t, double tt,
                            double& log_pdf, double& log_S) {
  const double lss = has_s ? fperace::entry_log_S(C.e[gs], tt) : 0.0;
  const double lst = has_t ? fperace::entry_log_S(C.e[gt], tt) : 0.0;
  log_S = lss + lst;
  if (has_s && has_t) {
    const double lps = fperace::entry_log_pdf(C.e[gs], tt);
    const double lpt = fperace::entry_log_pdf(C.e[gt], tt);
    log_pdf = roup_logaddexp(lps + lst, lpt + lss);
  } else {
    log_pdf = has_s ? fperace::entry_log_pdf(C.e[gs], tt)
                    : fperace::entry_log_pdf(C.e[gt], tt);
  }
}

inline fperace::SolveCache* roup_cache(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (ctx == nullptr) return nullptr;
  if (!ctx->fpe_cache) {
    ctx->fpe_cache = std::make_shared<fperace::SolveCache>();
    roup_configure_cache(*ctx->fpe_cache);
  }
  return ctx->fpe_cache.get();
}

inline void roup_prepare_rows(fperace::SolveCache& C, const double* rt,
                              const double* const* cols, int n_rows,
                              const int* isok) {
  if (C.prepared && (C.roup_local
      ? C.row_group_s.size() == static_cast<size_t>(n_rows)
      : C.row_group.size() == static_cast<size_t>(n_rows))) {
    return;
  }

  const double* v_S_   = cols[emc2col::roup::v_S];
  const double* v_T_   = cols[emc2col::roup::v_T];
  const double* tau_S_ = cols[emc2col::roup::tau_S];
  const double* tau_T_ = cols[emc2col::roup::tau_T];
  const double* k_     = cols[emc2col::roup::k];
  const double* B_     = cols[emc2col::roup::B];
  const double* A_     = cols[emc2col::roup::A];
  const double* t0_    = cols[emc2col::roup::t0];
  const double* s_     = cols[emc2col::roup::s];

  if (C.roup_local) {
    C.row_group_s.assign(n_rows, -1);
    C.row_group_t.assign(n_rows, -1);
    std::vector<fperace::Key> keys_s, keys_t;
    std::vector<double> horizon_s, horizon_t;
    std::vector<std::vector<double>> query_s, query_t;
    std::vector<int> row_s(n_rows, -1), row_t(n_rows, -1);

    auto add_key = [](const fperace::Key& p, double tt,
                      std::vector<fperace::Key>& keys,
                      std::vector<double>& horizons,
                      std::vector<std::vector<double>>& queries) {
      int g = -1;
      for (size_t j = 0; j < keys.size(); ++j) {
        if (keys[j] == p) { g = static_cast<int>(j); break; }
      }
      if (g < 0) {
        keys.push_back(p); horizons.push_back(tt);
        queries.push_back(std::vector<double>(1, tt));
        return static_cast<int>(keys.size() - 1);
      }
      horizons[g] = std::max(horizons[g], tt);
      queries[g].push_back(tt);
      return g;
    };

    for (int i = 0; i < n_rows; ++i) {
      if (!isok[i] || R_IsNA(v_S_[i])) continue;
      const double tt = rt[i] - t0_[i];
      if (!(tt > 0.0) || !emc2_isfinite(tt)) continue;
      fperace::Key ks, kt; bool hs = false, ht = false;
      if (!roup_local_keys(v_S_[i], v_T_[i], tau_S_[i], tau_T_[i], k_[i],
                           B_[i], A_[i], s_[i], roup_bnd_row(C.bnd_kind, cols, i),
                           ks, kt, hs, ht)) continue;
      if (hs) row_s[i] = add_key(ks, tt, keys_s, horizon_s, query_s);
      if (ht) row_t[i] = add_key(kt, tt, keys_t, horizon_t, query_t);
    }
    for (auto& q : query_s) {
      std::sort(q.begin(), q.end()); q.erase(std::unique(q.begin(), q.end()), q.end());
    }
    for (auto& q : query_t) {
      std::sort(q.begin(), q.end()); q.erase(std::unique(q.begin(), q.end()), q.end());
    }
    std::vector<int> idx_s, idx_t;
    fperace::cache_get_batch(C, keys_s, horizon_s, idx_s,
                             C.sparse_raw_output ? &query_s : nullptr);
    fperace::cache_get_batch(C, keys_t, horizon_t, idx_t,
                             C.sparse_raw_output ? &query_t : nullptr);
    for (int i = 0; i < n_rows; ++i) {
      if (row_s[i] >= 0) C.row_group_s[i] = idx_s[row_s[i]];
      if (row_t[i] >= 0) C.row_group_t[i] = idx_t[row_t[i]];
    }
    C.row_group.assign(n_rows, -1);
    C.prepared = true;
    return;
  }

  C.row_group.assign(n_rows, -1);

  // Pass 1: collect distinct keys and max horizon across all valid rows
  std::vector<fperace::Key> keys;
  std::vector<double> horizon;
  std::vector<std::vector<double>> query_times;
  std::vector<int> row_key_idx(n_rows, -1);

  for (int i = 0; i < n_rows; ++i) {
    if (!isok[i] || R_IsNA(v_S_[i])) continue;
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !emc2_isfinite(tt)) continue;
    fperace::Key p;
    if (!fperace::roup_key(v_S_[i], v_T_[i], tau_S_[i], tau_T_[i], k_[i], B_[i], A_[i], s_[i],
                           roup_bnd_row(C.bnd_kind, cols, i), p)) continue;

    int g = -1;
    for (size_t j = 0; j < keys.size(); ++j) {
      if (keys[j] == p) { g = static_cast<int>(j); break; }
    }
    if (g < 0) {
      keys.push_back(p);
      horizon.push_back(tt);
      query_times.push_back(std::vector<double>(1, tt));
      g = static_cast<int>(keys.size()) - 1;
    } else {
      if (tt > horizon[g]) horizon[g] = tt;
      query_times[g].push_back(tt);
    }
    row_key_idx[i] = g;
  }

  for (auto& times : query_times) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  // Pass 2: solve keys to max required horizon
  std::vector<int> to_cache(keys.size(), -1);
  fperace::cache_get_batch(
      C, keys, horizon, to_cache,
      C.sparse_raw_output ? &query_times : nullptr);
  for (int i = 0; i < n_rows; ++i) {
    if (row_key_idx[i] >= 0) C.row_group[i] = to_cache[row_key_idx[i]];
  }
  C.prepared = true;
}

inline void droup_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  roup_prepare_rows(*C, rt, cols, n_rows, isok);

  const double* t0_ = cols[emc2col::roup::t0];
  if (C->roup_local) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      const int gs = C->row_group_s[i], gt = C->row_group_t[i];
      if (gs < 0 && gt < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
      double lp = 0.0, ls = 0.0;
      roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0), gs >= 0, gt >= 0,
                      rt[i] - t0_[i], lp, ls);
      out[i] = raw_log_value(lp, min_ll, floor_raw);
    }
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) { out[i] = raw_log_zero(min_ll, floor_raw); continue; }
    const double lp = fperace::entry_log_pdf(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(lp, min_ll, floor_raw);
  }
}

inline void proup_raw(const double* rt, const double* const* cols, int n_rows,
                      const int* mask, const int* isok,
                      double* out, double min_ll, void* ctx_) {
  const bool floor_raw = raw_floor_log_lik(ctx_);
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  roup_prepare_rows(*C, rt, cols, n_rows, isok);

  const double* t0_ = cols[emc2col::roup::t0];
  if (C->roup_local) {
    for (int i = 0; i < n_rows; ++i) {
      if (!mask[i]) continue;
      const int gs = C->row_group_s[i], gt = C->row_group_t[i];
      if (gs < 0 && gt < 0) { out[i] = 0.0; continue; }
      double lp = 0.0, ls = 0.0;
      roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0), gs >= 0, gt >= 0,
                      rt[i] - t0_[i], lp, ls);
      out[i] = raw_log_value(ls, min_ll, floor_raw);
    }
    return;
  }
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int g = C->row_group[i];
    if (g < 0) {
      out[i] = 0.0;
      continue;
    }
    const double ls = fperace::entry_log_S(C->e[g], rt[i] - t0_[i]);
    out[i] = raw_log_value(ls, min_ll, floor_raw);
  }
}

inline void roup_logS_at_t(double t, const double* const* cols,
                           int n_rows_total, int n_lR, int /*n_par*/,
                           const int* trunc_mask, int n_unique_trials,
                           const int* isok_all, void* ctx_, double* logS_out) {
  fperace::SolveCache* C = roup_cache(ctx_);
  const double* v_S_   = cols[emc2col::roup::v_S];
  const double* v_T_   = cols[emc2col::roup::v_T];
  const double* tau_S_ = cols[emc2col::roup::tau_S];
  const double* tau_T_ = cols[emc2col::roup::tau_T];
  const double* k_     = cols[emc2col::roup::k];
  const double* B_     = cols[emc2col::roup::B];
  const double* A_     = cols[emc2col::roup::A];
  const double* t0_    = cols[emc2col::roup::t0];
  const double* s_     = cols[emc2col::roup::s];

  for (int j = 0; j < n_unique_trials; ++j) {
    if (!trunc_mask[j]) continue;
    const int start = j * n_lR;
    double logS = 0.0;
    bool bad = false;
    for (int m = 0; m < n_lR && !bad; ++m) {
      const int r = start + m;
      if (r >= n_rows_total) break;
      if (!isok_all[r] || R_IsNA(v_S_[r])) { bad = true; break; }
      const double tt = t - t0_[r];
      if (!(tt > 0.0)) continue;
      if (C->roup_local) {
        fperace::Key ks, kt; bool hs = false, ht = false;
        if (!roup_local_keys(v_S_[r], v_T_[r], tau_S_[r], tau_T_[r], k_[r],
                             B_[r], A_[r], s_[r], roup_bnd_row(C->bnd_kind, cols, r),
                             ks, kt, hs, ht)) {
          bad = true; break;
        }
        const int gs = hs ? fperace::cache_get(*C, ks, tt) : -1;
        const int gt = ht ? fperace::cache_get(*C, kt, tt) : -1;
        double lp = 0.0, lS = 0.0;
        roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0), hs, ht, tt, lp, lS);
        if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
        if (lS < 0.0) logS += lS;
        continue;
      }
      fperace::Key p;
      if (!fperace::roup_key(v_S_[r], v_T_[r], tau_S_[r], tau_T_[r], k_[r], B_[r],
                             A_[r], s_[r], roup_bnd_row(C->bnd_kind, cols, r), p)) {
        bad = true; break;
      }
      const int g = fperace::cache_get(*C, p, tt);
      const double lS = fperace::entry_log_S(C->e[g], tt);
      if (lS <= fperace::LOG_FLOOR) { bad = true; break; }
      if (lS < 0.0) logS += lS;
    }
    logS_out[j] = bad ? R_NegInf : logS;
  }
}

inline double roup_scalar_horizon(double tt) {
  return (tt > 1.0) ? (2.0 * tt) : (tt + 1.0);
}

inline double droup_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) return 0.0;
  if (R_IsNA(par[emc2col::roup::v_S])) return 0.0;
  const double tt = t - par[emc2col::roup::t0];
  if (!(tt > 0.0)) return 0.0;
  if (C->roup_local) {
    fperace::Key ks, kt; bool hs = false, ht = false;
    if (!roup_local_keys(par[emc2col::roup::v_S], par[emc2col::roup::v_T],
                         par[emc2col::roup::tau_S], par[emc2col::roup::tau_T],
                         par[emc2col::roup::k], par[emc2col::roup::B],
                         par[emc2col::roup::A], par[emc2col::roup::s],
                         roup_bnd_par(C->bnd_kind, par), ks, kt, hs, ht)) return 0.0;
    const int gs = hs ? fperace::cache_get(*C, ks, roup_scalar_horizon(tt)) : -1;
    const int gt = ht ? fperace::cache_get(*C, kt, roup_scalar_horizon(tt)) : -1;
    double lp = 0.0, lS = 0.0;
    roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0), hs, ht, tt, lp, lS);
    return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
  }
  fperace::Key p;
  if (!fperace::roup_key(par[emc2col::roup::v_S], par[emc2col::roup::v_T],
                         par[emc2col::roup::tau_S], par[emc2col::roup::tau_T],
                         par[emc2col::roup::k], par[emc2col::roup::B],
                         par[emc2col::roup::A], par[emc2col::roup::s],
                         roup_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, roup_scalar_horizon(tt));
  const double lp = fperace::entry_log_pdf(C->e[g], tt);
  return (lp <= fperace::LOG_FLOOR) ? 0.0 : std::exp(lp);
}

inline double proup_scalar(double t, const double* par, void* ctx_) {
  fperace::SolveCache* C = roup_cache(ctx_);
  if (C == nullptr) return 0.0;
  if (R_IsNA(par[emc2col::roup::v_S])) return 0.0;
  const double tt = t - par[emc2col::roup::t0];
  if (!(tt > 0.0)) return 0.0;
  if (C->roup_local) {
    fperace::Key ks, kt; bool hs = false, ht = false;
    if (!roup_local_keys(par[emc2col::roup::v_S], par[emc2col::roup::v_T],
                         par[emc2col::roup::tau_S], par[emc2col::roup::tau_T],
                         par[emc2col::roup::k], par[emc2col::roup::B],
                         par[emc2col::roup::A], par[emc2col::roup::s],
                         roup_bnd_par(C->bnd_kind, par), ks, kt, hs, ht)) return 0.0;
    const int gs = hs ? fperace::cache_get(*C, ks, roup_scalar_horizon(tt)) : -1;
    const int gt = ht ? fperace::cache_get(*C, kt, roup_scalar_horizon(tt)) : -1;
    double lp = 0.0, lS = 0.0;
    roup_local_eval(*C, std::max(gs, 0), std::max(gt, 0), hs, ht, tt, lp, lS);
    if (lS >= 0.0) return 0.0;
    if (lS <= fperace::LOG_FLOOR) return 1.0;
    return -std::expm1(lS);
  }
  fperace::Key p;
  if (!fperace::roup_key(par[emc2col::roup::v_S], par[emc2col::roup::v_T],
                         par[emc2col::roup::tau_S], par[emc2col::roup::tau_T],
                         par[emc2col::roup::k], par[emc2col::roup::B],
                         par[emc2col::roup::A], par[emc2col::roup::s],
                         roup_bnd_par(C->bnd_kind, par), p)) return 0.0;
  const int g = fperace::cache_get(*C, p, roup_scalar_horizon(tt));
  const double lS = fperace::entry_log_S(C->e[g], tt);
  if (lS >= 0.0) return 0.0;
  if (lS <= fperace::LOG_FLOOR) return 1.0;
  return -std::expm1(lS);
}

#endif // model_ROUp_h
