// Standalone equivalence check for the RLF matrix-free stage solver.
//
// The package selects between the dense inverse and the FFT/BiCGSTAB path by
// a compile-time grid-size threshold, so the R test suite cannot exercise both
// on one grid.  This driver builds the same operator and compares, at
// identical n:
//   * the matrix-free matvec against a dense M*x, and
//   * the preconditioned BiCGSTAB solve against LAPACK dgesv.
// Both agree to ~1e-13 relative, at five to seven iterations independent of n
// and alpha.
//
// Build and run (from the package root):
//   g++ -std=gnu++17 -O3 -march=native -ffast-math -I src \
//       -I"$(R RHOME)/include" -o /tmp/check_rlf_matrix_free \
//       WorkingTests/check_rlf_matrix_free.cpp \
//       $(R CMD config LAPACK_LIBS) $(R CMD config BLAS_LIBS) \
//       -L"$(R RHOME)/lib" -lR
//   /tmp/check_rlf_matrix_free <alpha> <n> <dt>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>

#define R_NO_REMAP
#include <R_ext/RS.h>
#include "model_RLF.h"
#include "rlf_toeplitz.h"

extern "C" {
  void F77_NAME(dgesv)(const int* n, const int* nrhs, double* a,
                       const int* lda, int* ipiv, double* b,
                       const int* ldb, int* info);
}

int main(int argc, char** argv) {
  const double alpha = argc > 1 ? atof(argv[1]) : 1.5;
  const int n = argc > 2 ? atoi(argv[2]) : 512;
  const double dt = argc > 3 ? atof(argv[3]) : 8e-3;

  rlf::RLF_Model m;
  m.v = 1.0; m.sigma = 1.0; m.alpha = alpha; m.b0 = 1.3; m.z0 = 0.3;

  const double x_lo = -8.0;
  const double h = (m.b0 - x_lo) / n;
  rlf::RLF_Operator op = rlf::build_rlf_operator(m, n, h);

  const double theta_dt = 0.2928932188134524 * dt;

  // Dense reference: M = I - theta*dt*L, column-major for LAPACK.
  std::vector<double> A(static_cast<size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      A[static_cast<size_t>(j) * n + i] =
        (i == j ? 1.0 : 0.0) - theta_dt * op.L[static_cast<size_t>(i) * n + j];

  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  std::vector<double> b(n), xd(n), xt(n, 0.0);
  for (int i = 0; i < n; ++i) { b[i] = u(rng); xd[i] = b[i]; }

  std::vector<int> ipiv(n);
  int nrhs = 1, info = 0;
  F77_CALL(dgesv)(&n, &nrhs, A.data(), &n, ipiv.data(), xd.data(), &n, &info);
  if (info != 0) { printf("dgesv info=%d\n", info); return 1; }

  rlf::RLF_ToeplitzSystem sys =
    rlf::rlf_build_toeplitz_system(op.L, n, theta_dt);

  // Structural check: dense M*x vs matrix-free M*x.
  std::vector<double> y_dense(n, 0.0), y_free(n);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < n; ++j) {
      const double mij = (i == j ? 1.0 : 0.0) -
        theta_dt * op.L[static_cast<size_t>(i) * n + j];
      acc += mij * b[j];
    }
    y_dense[i] = acc;
  }
  rlf::rlf_toeplitz_apply(sys, b.data(), y_free.data());
  double mv = 0.0, mvs = 0.0;
  for (int i = 0; i < n; ++i) {
    mv = std::max(mv, std::abs(y_dense[i] - y_free[i]));
    mvs = std::max(mvs, std::abs(y_dense[i]));
  }

  const int iters = rlf::rlf_toeplitz_solve(sys, b.data(), xt.data(), 1e-13, 300);
  double err = 0.0, scale = 0.0;
  for (int i = 0; i < n; ++i) {
    err = std::max(err, std::abs(xd[i] - xt[i]));
    scale = std::max(scale, std::abs(xd[i]));
  }
  printf("alpha=%.2f n=%4d upwind=%d  matvec rel=%.3e  solve iters=%3d rel=%.3e\n",
         alpha, n, (int)op.upwind_drift, mv / mvs, iters, err / scale);
  return 0;
}
