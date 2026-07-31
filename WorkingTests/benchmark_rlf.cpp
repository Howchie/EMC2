#include <Rcpp.h>
#include <chrono>
#include "../src/model_RLF.h"

using namespace Rcpp;

extern "C" {
void F77_NAME(dgetrs)(const char*, const int*, const int*, const double*,
                      const int*, const int*, double*, const int*, int*,
                      size_t);
}

template <class Function>
double elapsed_ms(Function&& function) {
  const auto start = std::chrono::steady_clock::now();
  function();
  const auto end = std::chrono::steady_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
Rcpp::List benchmark_rlf_optimizations(int nx = 100, int repeats = 6) {
  const int n_keys = 8;
  const double horizon = 2.0;
  std::vector<rlf::Key> keys(n_keys);
  std::vector<double> horizons(n_keys, horizon);
  std::vector<std::vector<double>> queries(
    n_keys, std::vector<double>{0.25, 0.5, 1.0, horizon});
  for (int i = 0; i < n_keys; ++i) {
    rlf::rlf_key(
      1.0 + 0.15 * i, 1.0, std::min(2.0, 1.3 + 0.1 * i),
      0.9 + 0.03 * (i % 2), 0.2 * (i % 2), keys[i]);
  }

  rlf::Grid grid;
  grid.nx = nx;
  grid.dt_target = 0.01;
  grid.tgrade = 1.0;
  grid.adaptive = false;
  grid.explicit_inverse = true;
  grid.sparse_output = true;
  // Keep the component benchmarks single-grid.  Richardson disables the SIMD
  // lane path and would make "inverse versus LU" include a second resolution.
  grid.richardson = false;

  double sink = 0.0;
  const double naive_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      for (int key = 0; key < n_keys; ++key) {
        rlf::Entry entry;
        rlf::rlf_cache_solve(
          keys[key], horizon, grid, &queries[key], entry);
        sink += entry.log_pdf.back();
      }
    }
  });

  rlf::SolveCache cached;
  cached.grid = grid;
  const double cached_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      for (int key = 0; key < n_keys; ++key) {
        const int index = rlf::cache_get(cached, keys[key], horizon);
        sink += cached.entries[index].log_pdf.back();
      }
    }
  });

  auto fill_batch = [&](bool simd) {
    rlf::SolveCache cache;
    cache.grid = grid;
    cache.grid.simd_batch = simd;
    std::vector<int> index;
    rlf::cache_get_batch(cache, keys, horizons, index, &queries);
    for (int i : index) sink += cache.entries[i].log_pdf.back();
  };
  const double serial_batch_ms = elapsed_ms([&]() { fill_batch(false); });
  const double simd_batch_ms = elapsed_ms([&]() { fill_batch(true); });

  rlf::RLF_Model model;
  model.v = 1.5;
  model.sigma = 1.0;
  model.alpha = 1.7;
  model.b0 = 1.2;
  model.z0 = 0.2;
  rlf::RLF_Result lu;
  const double lu_ms = elapsed_ms([&]() {
    lu = rlf::rlf_solve(model, horizon, nx, 200, false, 1.0, false);
  });
  rlf::RLF_Result inverse;
  const double inverse_ms = elapsed_ms([&]() {
    inverse =
      rlf::rlf_solve(model, horizon, nx, 200, false, 1.0, true);
  });
  const rlf::RLF_Comparison inverse_error =
    rlf::compare_rlf_results(lu, inverse, horizon);

  auto timed_pair = [&](int pair_nx, double dt, double ratio) {
    rlf::SolveCache pair_cache;
    pair_cache.grid = grid;
    pair_cache.grid.nx = pair_nx;
    pair_cache.grid.dt_target = dt;
    pair_cache.grid.richardson = true;
    pair_cache.grid.richardson_ratio = ratio;
    std::vector<int> index;
    const double ms = elapsed_ms([&]() {
      rlf::cache_get_batch(
        pair_cache, keys, horizons, index, &queries);
    });
    return ms;
  };
  const double former_default_ms = timed_pair(128, 0.008, 1.5);
  const double production_default_ms = timed_pair(160, 0.016, 1.25);

  return Rcpp::List::create(
    _["nx"] = nx,
    _["rows"] = repeats * n_keys,
    _["unique_keys"] = n_keys,
    _["naive_ms"] = naive_ms,
    _["cached_ms"] = cached_ms,
    _["cache_speedup"] = naive_ms / cached_ms,
    _["cache_solves"] = static_cast<int>(cached.solve_count),
    _["serial_batch_ms"] = serial_batch_ms,
    _["simd_batch_ms"] = simd_batch_ms,
    _["simd_speedup"] = serial_batch_ms / simd_batch_ms,
    _["lu_ms"] = lu_ms,
    _["inverse_ms"] = inverse_ms,
    _["inverse_speedup"] = lu_ms / inverse_ms,
    _["inverse_pdf_error"] = inverse_error.pdf,
    _["inverse_cdf_error"] = inverse_error.cdf,
    _["former_default_ms"] = former_default_ms,
    _["production_default_ms"] = production_default_ms,
    _["production_default_speedup"] =
      former_default_ms / production_default_ms,
    _["sink"] = sink);
}

// [[Rcpp::export]]
Rcpp::List benchmark_rlf_matvec(int nx = 128, int repeats = 20000) {
  rlf::RLF_Model model;
  model.v = 1.5;
  model.sigma = 1.0;
  model.alpha = 1.5;
  model.b0 = 2.0;
  model.z0 = 0.5;
  const double extent = rlf::rlf_lower_extent(model, 2.0, nx);
  const int n = nx - 1;
  const double h = (model.b0 + extent) / nx;
  const rlf::RLF_Operator op = rlf::build_rlf_operator(model, n, h);
  const std::vector<double> lhs = rlf::rlf_trbdf2_lhs(op, 0.008);
  rlf::RLF_DenseInverse inverse;
  if (!inverse.build(n, lhs)) Rcpp::stop("inverse failed");
  const std::vector<double> matrix =
    rlf::rlf_trbdf2_step_matrix(inverse.value, n);
  std::vector<double> x(n), manual(n), blas(n);
  for (int i = 0; i < n; ++i) x[i] = 1.0 / (1.0 + i);

  const double manual_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      for (int i = 0; i < n; ++i) {
        const double* __restrict row =
          &matrix[static_cast<size_t>(i) * n];
        double sum = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
        for (int j = 0; j < n; ++j) sum += row[j] * x[j];
        manual[i] = sum;
      }
      x.swap(manual);
    }
  });
  const std::vector<double> manual_final = x;
  for (int i = 0; i < n; ++i) x[i] = 1.0 / (1.0 + i);
  const double blas_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      rlf::rlf_dense_matvec(matrix, x, blas, n);
      x.swap(blas);
    }
  });
  double max_error = 0.0;
  for (int i = 0; i < n; ++i) {
    max_error = std::max(max_error, std::abs(x[i] - manual_final[i]));
  }
  return Rcpp::List::create(
    _["nx"] = nx,
    _["repeats"] = repeats,
    _["manual_ms"] = manual_ms,
    _["blas_ms"] = blas_ms,
    _["speedup"] = manual_ms / blas_ms,
    _["max_error"] = max_error);
}

// Compare the two ways to apply the collapsed TR--BDF2 update.  Forming
// 2a M^-2 - (a+b) M^-1 costs one dense matrix multiply per parameter key,
// after which every time step is one GEMV.  Applying M^-1 twice avoids that
// setup cost but needs two GEMVs per step.  The crossover depends on both nx
// and the number of time steps, so measure the complete per-key work rather
// than comparing isolated BLAS calls.
// [[Rcpp::export]]
Rcpp::List benchmark_rlf_propagator_strategy(
    int nx = 160, int steps = 125, int repeats = 20) {
  rlf::RLF_Model model;
  model.v = 1.5;
  model.sigma = 1.0;
  model.alpha = 1.3;
  model.b0 = 2.0;
  model.z0 = 0.5;
  const double extent = rlf::rlf_lower_extent(model, 2.0, nx);
  const int n = nx - 1;
  const double h = (model.b0 + extent) / nx;
  const rlf::RLF_Operator op = rlf::build_rlf_operator(model, n, h);
  const std::vector<double> lhs =
    rlf::rlf_trbdf2_lhs(op, 2.0 / std::max(1, steps));
  rlf::RLF_DenseInverse inverse;
  if (!inverse.build(n, lhs)) Rcpp::stop("inverse failed");

  std::vector<double> initial;
  rlf::rlf_initial_density(model, -extent, h, n, initial);
  std::vector<double> x, next, stage, second;
  double sink = 0.0;
  const double precomputed_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      const std::vector<double> matrix =
        rlf::rlf_trbdf2_step_matrix(inverse.value, n);
      x = initial;
      for (int step = 0; step < steps; ++step) {
        rlf::rlf_dense_matvec(matrix, x, next, n);
        x.swap(next);
      }
      sink += x[n / 2];
    }
  });
  const std::vector<double> precomputed_final = x;

  const double staged_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      x = initial;
      for (int step = 0; step < steps; ++step) {
        rlf::rlf_dense_matvec(inverse.value, x, stage, n);
        rlf::rlf_dense_matvec(inverse.value, stage, second, n);
        next.resize(n);
        for (int i = 0; i < n; ++i) {
          next[i] =
            2.0 * rlf::RLF_TRBDF2_A * second[i] -
            (rlf::RLF_TRBDF2_A + rlf::RLF_TRBDF2_B) * stage[i];
        }
        x.swap(next);
      }
      sink += x[n / 2];
    }
  });
  double max_error = 0.0;
  for (int i = 0; i < n; ++i) {
    max_error = std::max(
      max_error, std::abs(x[i] - precomputed_final[i]));
  }

  return Rcpp::List::create(
    _["nx"] = nx,
    _["steps"] = steps,
    _["repeats"] = repeats,
    _["precomputed_ms"] = precomputed_ms,
    _["staged_ms"] = staged_ms,
    _["staged_speedup"] = precomputed_ms / staged_ms,
    _["max_error"] = max_error,
    _["sink"] = sink);
}

// [[Rcpp::export]]
Rcpp::List benchmark_rlf_lapack_solve(
    int nx = 160, int steps = 125, int repeats = 20) {
  rlf::RLF_Model model;
  model.v = 1.5;
  model.sigma = 1.0;
  model.alpha = 1.3;
  model.b0 = 2.0;
  model.z0 = 0.5;
  const double extent = rlf::rlf_lower_extent(model, 2.0, nx);
  const int n = nx - 1;
  const double h = (model.b0 + extent) / nx;
  const rlf::RLF_Operator op = rlf::build_rlf_operator(model, n, h);
  const std::vector<double> lhs =
    rlf::rlf_trbdf2_lhs(op, 2.0 / std::max(1, steps));
  std::vector<double> initial;
  rlf::rlf_initial_density(model, -extent, h, n, initial);

  std::vector<double> x, next, stage;
  double sink = 0.0;
  const double precomputed_ms = elapsed_ms([&]() {
    for (int repeat = 0; repeat < repeats; ++repeat) {
      rlf::RLF_DenseInverse inverse;
      if (!inverse.build(n, lhs)) Rcpp::stop("inverse failed");
      const std::vector<double> matrix =
        rlf::rlf_trbdf2_step_matrix(inverse.value, n);
      x = initial;
      for (int startup = 0; startup < rlf::RLF_STARTUP_STEPS; ++startup) {
        rlf::rlf_dense_matvec(inverse.value, x, next, n);
        x.swap(next);
      }
      for (int step = 0; step < steps; ++step) {
        rlf::rlf_dense_matvec(matrix, x, next, n);
        x.swap(next);
      }
      sink += x[n / 2];
    }
  });
  const std::vector<double> precomputed_final = x;

  const double lapack_ms = elapsed_ms([&]() {
    const char no_transpose = 'N';
    const int nrhs = 1;
    const int lda = n;
    const int ldb = n;
    for (int repeat = 0; repeat < repeats; ++repeat) {
      std::vector<double> factor(static_cast<size_t>(n) * n);
      for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
          factor[static_cast<size_t>(col) * n + row] =
            lhs[static_cast<size_t>(row) * n + col];
        }
      }
      std::vector<int> pivot(n);
      int info = 0;
      F77_CALL(dgetrf)(
        &n, &n, factor.data(), &lda, pivot.data(), &info);
      if (info != 0) Rcpp::stop("factorization failed");
      auto solve = [&](std::vector<double>& value) {
        F77_CALL(dgetrs)(
          &no_transpose, &n, &nrhs, factor.data(), &lda, pivot.data(),
          value.data(), &ldb, &info, 1);
        if (info != 0) Rcpp::stop("solve failed");
      };

      x = initial;
      for (int startup = 0; startup < rlf::RLF_STARTUP_STEPS; ++startup) {
        next = x;
        solve(next);
        x.swap(next);
      }
      for (int step = 0; step < steps; ++step) {
        stage = x;
        solve(stage);
        next = stage;
        solve(next);
        for (int i = 0; i < n; ++i) {
          next[i] =
            2.0 * rlf::RLF_TRBDF2_A * next[i] -
            (rlf::RLF_TRBDF2_A + rlf::RLF_TRBDF2_B) * stage[i];
        }
        x.swap(next);
      }
      sink += x[n / 2];
    }
  });
  double max_error = 0.0;
  for (int i = 0; i < n; ++i) {
    max_error = std::max(
      max_error, std::abs(x[i] - precomputed_final[i]));
  }

  return Rcpp::List::create(
    _["nx"] = nx,
    _["steps"] = steps,
    _["repeats"] = repeats,
    _["precomputed_ms"] = precomputed_ms,
    _["lapack_ms"] = lapack_ms,
    _["lapack_speedup"] = precomputed_ms / lapack_ms,
    _["max_error"] = max_error,
    _["sink"] = sink);
}
