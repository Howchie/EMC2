#ifndef model_RLF_kernels_h
#define model_RLF_kernels_h

// EMC2 race-kernel adapter for the nonlocal RLF solver.  This header is included
// from utils.h after ContextForRaceModels and raw_log_* are defined.

inline SEXP rlf_option(const char* primary, const char* shared = nullptr) {
  SEXP value = Rf_GetOption1(Rf_install(primary));
  if (value == R_NilValue && shared != nullptr) {
    value = Rf_GetOption1(Rf_install(shared));
  }
  return value;
}

inline void rlf_configure_grid(rlf::Grid& grid) {
  auto option_number = [&](const char* primary, const char* shared,
                           double& destination, double lower) {
    SEXP value = rlf_option(primary, shared);
    if (value == R_NilValue || Rf_length(value) < 1) return;
    const double x = Rf_asReal(value);
    if (R_finite(x) && x >= lower) destination = x;
  };
  auto option_logical = [&](const char* name, bool& destination) {
    SEXP value = rlf_option(name);
    if (value == R_NilValue || Rf_length(value) < 1) return;
    const int x = Rf_asLogical(value);
    if (x != NA_LOGICAL) destination = x;
  };

  double nx = grid.nx;
  double dt = grid.dt_target;
  double tgrade = grid.tgrade;
  option_number("emc2.rlf_nx", "emc2.fpe_nx", nx, 30.0);
  option_number("emc2.rlf_dt", "emc2.fpe_dt", dt, 1e-8);
  option_number("emc2.rlf_tgrade", "emc2.fpe_tgrade", tgrade, 1.0);
  grid.nx = static_cast<int>(nx);
  grid.dt_target = dt;
  grid.tgrade = tgrade;
  grid.nt_max = std::max(
    grid.nt_max, static_cast<int>(std::ceil(10.0 / grid.dt_target)));
  option_logical("emc2.rlf_adaptive", grid.adaptive);
  option_logical("emc2.rlf_explicit_inverse", grid.explicit_inverse);
  option_logical("emc2.rlf_sparse_output", grid.sparse_output);
  option_logical("emc2.rlf_simd_batch", grid.simd_batch);
  option_logical("emc2.rlf_horizon_split", grid.horizon_split);
  // Global rather than a Grid field: build_rlf_operator sees only RLF_Model.
  // Written once here on the main thread, read-only from the worker marches.
  option_logical("emc2.rlf_force_centred", rlf::rlf_force_centred_drift);
}

inline rlf::SolveCache* rlf_cache(void* context) {
  auto* ctx = static_cast<ContextForRaceModels*>(context);
  if (ctx == nullptr) return nullptr;
  if (!ctx->rlf_cache) {
    ctx->rlf_cache = std::make_shared<rlf::SolveCache>();
    rlf_configure_grid(ctx->rlf_cache->grid);
  }
  return ctx->rlf_cache.get();
}

inline void rlf_prepare_rows(rlf::SolveCache& cache, const double* rt,
                             const double* const* cols, int n_rows,
                             const int* isok) {
  const double* v = cols[emc2col::rlf::v];
  const double* B = cols[emc2col::rlf::B];
  const double* A = cols[emc2col::rlf::A];
  const double* t0 = cols[emc2col::rlf::t0];
  const double* s = cols[emc2col::rlf::s];
  const double* alpha = cols[emc2col::rlf::alpha];

  cache.row_group.assign(n_rows, -1);
  std::vector<rlf::Key> keys;
  std::vector<double> horizons;
  std::vector<std::vector<double>> query_times;
  std::vector<int> row_key(n_rows, -1);
  std::unordered_map<rlf::Key, int, rlf::KeyHash> groups;

  for (int i = 0; i < n_rows; ++i) {
    if (!isok[i] || R_IsNA(v[i])) continue;
    const double tt = rt[i] - t0[i];
    if (!(tt > 0.0) || !emc2_isfinite(tt)) continue;
    rlf::Key key;
    if (!rlf::rlf_key(v[i], s[i], alpha[i], B[i], A[i], key)) continue;
    // Bucket on this row's own horizon, so the group a row lands in is the one
    // whose domain was sized for it.
    if (cache.grid.horizon_split) {
      key.bucket = rlf::rlf_horizon_bucket(key, tt);
    }
    const auto found = groups.find(key);
    int group = -1;
    if (found == groups.end()) {
      group = static_cast<int>(keys.size());
      groups.emplace(key, group);
      keys.push_back(key);
      horizons.push_back(tt);
      query_times.push_back(std::vector<double>(1, tt));
    } else {
      group = found->second;
      horizons[group] = std::max(horizons[group], tt);
      query_times[group].push_back(tt);
    }
    row_key[i] = group;
  }

  for (auto& times : query_times) {
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
  }

  std::vector<int> cache_index;
  rlf::cache_get_batch(
    cache, keys, horizons, cache_index,
    cache.grid.sparse_output ? &query_times : nullptr);
  for (int i = 0; i < n_rows; ++i) {
    if (row_key[i] >= 0) cache.row_group[i] = cache_index[row_key[i]];
  }
}

inline void drlf_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok, double* out,
                     double min_ll, void* context) {
  const bool floor_raw = raw_floor_log_lik(context);
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) {
    for (int i = 0; i < n_rows; ++i) {
      if (mask[i]) out[i] = raw_log_zero(min_ll, floor_raw);
    }
    return;
  }
  rlf_prepare_rows(*cache, rt, cols, n_rows, isok);
  const double* t0 = cols[emc2col::rlf::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int group = cache->row_group[i];
    if (group < 0) {
      out[i] = raw_log_zero(min_ll, floor_raw);
      continue;
    }
    out[i] = raw_log_value(
      rlf::entry_log_pdf(cache->entries[group], rt[i] - t0[i]),
      min_ll, floor_raw);
  }
}

inline void prlf_raw(const double* rt, const double* const* cols, int n_rows,
                     const int* mask, const int* isok, double* out,
                     double min_ll, void* context) {
  const bool floor_raw = raw_floor_log_lik(context);
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) {
    for (int i = 0; i < n_rows; ++i) if (mask[i]) out[i] = 0.0;
    return;
  }
  rlf_prepare_rows(*cache, rt, cols, n_rows, isok);
  const double* t0 = cols[emc2col::rlf::t0];
  for (int i = 0; i < n_rows; ++i) {
    if (!mask[i]) continue;
    const int group = cache->row_group[i];
    if (group < 0) {
      out[i] = 0.0;
      continue;
    }
    out[i] = raw_log_value(
      rlf::entry_log_S(cache->entries[group], rt[i] - t0[i]),
      min_ll, floor_raw);
  }
}

inline double rlf_scalar_horizon(double tt) {
  return tt > 1.0 ? 2.0 * tt : tt + 1.0;
}

// The scalar entry points share the raw path's cache, so they must bucket on
// the same rule; otherwise a single long-horizon censoring query would grow a
// bulk entry's march and undo the split for every row using it.
inline bool rlf_key_from_par(const double* par, double horizon,
                             const rlf::Grid& grid, rlf::Key& key) {
  if (!rlf::rlf_key(
        par[emc2col::rlf::v], par[emc2col::rlf::s],
        par[emc2col::rlf::alpha], par[emc2col::rlf::B],
        par[emc2col::rlf::A], key)) {
    return false;
  }
  if (grid.horizon_split) key.bucket = rlf::rlf_horizon_bucket(key, horizon);
  return true;
}

inline double drlf_scalar(double time, const double* par, void* context) {
  if (R_IsNA(par[emc2col::rlf::v])) return 0.0;
  const double tt = time - par[emc2col::rlf::t0];
  if (!(tt > 0.0)) return 0.0;
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) return 0.0;
  const double horizon = rlf_scalar_horizon(tt);
  rlf::Key key;
  if (!rlf_key_from_par(par, horizon, cache->grid, key)) return 0.0;
  const int group = rlf::cache_get(*cache, key, horizon);
  const double log_pdf =
    rlf::entry_log_pdf(cache->entries[group], tt);
  return log_pdf <= rlf::RLF_LOG_FLOOR ? 0.0 : std::exp(log_pdf);
}

inline double prlf_scalar(double time, const double* par, void* context) {
  if (R_IsNA(par[emc2col::rlf::v])) return 0.0;
  const double tt = time - par[emc2col::rlf::t0];
  if (!(tt > 0.0)) return 0.0;
  rlf::SolveCache* cache = rlf_cache(context);
  if (cache == nullptr) return 0.0;
  const double horizon = rlf_scalar_horizon(tt);
  rlf::Key key;
  if (!rlf_key_from_par(par, horizon, cache->grid, key)) return 0.0;
  const int group = rlf::cache_get(*cache, key, horizon);
  const double log_survivor =
    rlf::entry_log_S(cache->entries[group], tt);
  if (log_survivor >= 0.0) return 0.0;
  if (log_survivor <= rlf::RLF_LOG_FLOOR) return 1.0;
  return -std::expm1(log_survivor);
}

inline void rlf_logS_at_t(double time, const double* const* cols,
                          int n_rows_total, int n_lR, int,
                          const int* trunc_mask, int n_unique_trials,
                          const int* isok_all, void* context,
                          double* logS_out) {
  rlf::SolveCache* cache = rlf_cache(context);
  const double* v = cols[emc2col::rlf::v];
  const double* B = cols[emc2col::rlf::B];
  const double* A = cols[emc2col::rlf::A];
  const double* t0 = cols[emc2col::rlf::t0];
  const double* s = cols[emc2col::rlf::s];
  const double* alpha = cols[emc2col::rlf::alpha];

  for (int trial = 0; trial < n_unique_trials; ++trial) {
    if (!trunc_mask[trial]) continue;
    const int start = trial * n_lR;
    double total = 0.0;
    bool bad = false;
    for (int accumulator = 0; accumulator < n_lR && !bad; ++accumulator) {
      const int row = start + accumulator;
      if (row >= n_rows_total) break;
      if (!isok_all[row] || R_IsNA(v[row])) {
        bad = true;
        break;
      }
      const double tt = time - t0[row];
      if (!(tt > 0.0)) continue;
      rlf::Key key;
      if (!rlf::rlf_key(v[row], s[row], alpha[row], B[row], A[row], key)) {
        bad = true;
        break;
      }
      if (cache->grid.horizon_split) {
        key.bucket = rlf::rlf_horizon_bucket(key, tt);
      }
      const int group = rlf::cache_get(*cache, key, tt);
      const double log_survivor =
        rlf::entry_log_S(cache->entries[group], tt);
      if (log_survivor <= rlf::RLF_LOG_FLOOR) {
        bad = true;
      } else if (log_survivor < 0.0) {
        total += log_survivor;
      }
    }
    logS_out[trial] = bad ? R_NegInf : total;
  }
}

#endif
