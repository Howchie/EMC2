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
// Column order per src/col_registry.h: v_S, v_T/E_T, tau_S, tau_T, k, B, A,
// t0, s. The area chart supplies E_T = v_T * tau_T and is mapped to v_T before
// the FPE solve.
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

struct RoupColIdx {
  int v_S, transient, tau_S, tau_T, k, B, A, t0, s;
};

inline RoupColIdx roup_col_idx(int par_kind) {
  if (par_kind == fperace::ROUP_PAR_AREA) {
    return {emc2col::roup_area::v_S, emc2col::roup_area::E_T,
            emc2col::roup_area::tau_S, emc2col::roup_area::tau_T,
            emc2col::roup_area::k, emc2col::roup_area::B,
            emc2col::roup_area::A, emc2col::roup_area::t0,
            emc2col::roup_area::s};
  }
  return {emc2col::roup::v_S, emc2col::roup::v_T,
          emc2col::roup::tau_S, emc2col::roup::tau_T,
          emc2col::roup::k, emc2col::roup::B,
          emc2col::roup::A, emc2col::roup::t0,
          emc2col::roup::s};
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
  if (C.prepared && C.row_group.size() == static_cast<size_t>(n_rows)) {
    return;
  }

  const RoupColIdx ci = roup_col_idx(C.par_kind);
  const double* v_S_   = cols[ci.v_S];
  const double* transient_ = cols[ci.transient];
  const double* tau_S_ = cols[ci.tau_S];
  const double* tau_T_ = cols[ci.tau_T];
  const double* k_     = cols[ci.k];
  const double* B_     = cols[ci.B];
  const double* A_     = cols[ci.A];
  const double* t0_    = cols[ci.t0];
  const double* s_     = cols[ci.s];

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
    if (!fperace::roup_key_par(C.par_kind, v_S_[i], transient_[i], tau_S_[i],
                               tau_T_[i], k_[i], B_[i], A_[i], s_[i],
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

void droup_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);

void proup_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);

void roup_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);

inline double roup_scalar_horizon(double tt) {
  return (tt > 1.0) ? (2.0 * tt) : (tt + 1.0);
}

double droup_scalar(double t, const double* par, void* ctx_);

double proup_scalar(double t, const double* par, void* ctx_);

#endif // model_ROUp_h
