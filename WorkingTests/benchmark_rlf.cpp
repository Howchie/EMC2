#include <Rcpp.h>
#include <chrono>
#include "../src/model_RLF.h"

using namespace Rcpp;

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
    _["sink"] = sink);
}
