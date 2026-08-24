#ifndef model_GOM_h
#define model_GOM_h

// Gompertz growth-process race kernels.
//
//   dX = alpha X log(K / X) dt + beta X dW,
//   X(0) ~ Uniform(1, 1 + A), absorbed at b = 1 + B + A.
//
// On Y = log(X), this is an OU process with
//   dY = (alpha log(K) - beta^2/2 - alpha Y) dt + beta dW.
// The cache stores that OU representation, while fpe_seed() uses the
// log_state flag to retain the physical-uniform start distribution.  This is
// important: A is a range on X, not a range on log(X).

#include "model_GOM_core.h"

inline fperace::SolveCache* gomp_cache(void* ctx_) {
  auto* ctx = static_cast<ContextForRaceModels*>(ctx_);
  if (ctx == nullptr) return nullptr;
  if (!ctx->fpe_cache) {
    ctx->fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_cache(*ctx->fpe_cache);
  }
  return ctx->fpe_cache.get();
}

inline void gomp_prepare_rows(fperace::SolveCache& C, const double* rt,
                              const double* const* cols, int n_rows,
                              const int* isok) {
  if (C.prepared && C.row_group.size() == static_cast<size_t>(n_rows)) return;

  const GompColIdx ci = gomp_col_idx();
  const double* alpha_ = cols[ci.alpha];
  const double* beta_ = cols[ci.beta];
  const double* K_ = cols[ci.K];
  const double* B_ = cols[ci.B];
  const double* A_ = cols[ci.A];
  const double* t0_ = cols[ci.t0];
  C.row_group.assign(n_rows, -1);

  std::vector<fperace::Key> keys;
  std::vector<double> horizons;
  std::vector<int> row_key_idx(n_rows, -1);
  for (int i = 0; i < n_rows; ++i) {
    if (!isok[i] || R_IsNA(alpha_[i])) continue;
    const double tt = rt[i] - t0_[i];
    if (!(tt > 0.0) || !emc2_isfinite(tt)) continue;
    fperace::Key p;
    if (!gomp_key(alpha_[i], beta_[i], K_[i], B_[i], A_[i],
                  gomp_bnd_row(C.bnd_kind, cols, i), p)) continue;
    int g = -1;
    for (size_t j = 0; j < keys.size(); ++j) {
      if (keys[j] == p) { g = static_cast<int>(j); break; }
    }
    if (g < 0) {
      keys.push_back(p);
      horizons.push_back(tt);
      g = static_cast<int>(keys.size()) - 1;
    } else if (tt > horizons[g]) {
      horizons[g] = tt;
    }
    row_key_idx[i] = g;
  }

  for (size_t j = 0; j < keys.size(); ++j) {
    const int g = fperace::cache_get(C, keys[j], horizons[j]);
    for (int i = 0; i < n_rows; ++i)
      if (row_key_idx[i] == static_cast<int>(j)) C.row_group[i] = g;
  }
  C.prepared = true;
}

void dgomp_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_);

void pgomp_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok, double* out,
               double min_ll, void* ctx_);

double dgomp_scalar(double t, const double* par, void* ctx_);

double pgomp_scalar(double t, const double* par, void* ctx_);

void gomp_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);


#endif // model_GOM_h
