#ifndef rlf_toeplitz_h
#define rlf_toeplitz_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Matrix-free TR-BDF2 stage solves for the RLF generator.
//
// The nonlocal generator is dense, but only because the jump kernel couples
// every pair of cells: its entries depend on i-j alone.  Writing the operator
// as
//
//   L = Toeplitz(t) + diag(d),      t[0] := 0,  d[j] := L[j][j],
//
// costs nothing to extract (off-diagonal entries are already translation
// invariant, and every state-dependent term -- the censoring rate added back
// at the lower edge and the sealed advective face at x_lo -- lives on the
// diagonal) and buys an O(n log n) matvec via circulant embedding.
//
// The dense path pays O(n^3) once to invert M = I - (gamma/2) dt L and form
// the TR-BDF2 step matrix, then O(n^2) per step.  That setup dominates as soon
// as the grid is refined, which is exactly what heavy tails demand: the
// domain must widen like tol^(-1/alpha) while the absorbing boundary still
// needs a fine h.  Solving each stage iteratively instead removes the cubic
// term entirely, so the cost becomes O(iters * n log n) per stage.
//
// M is a small perturbation of a circulant, so a Strang circulant
// preconditioner clusters the spectrum and BiCGSTAB converges in a handful of
// iterations independently of n -- the standard result for fractional
// diffusion operators.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace rlf {

using rlf_cplx = std::complex<double>;

// ---------------------------------------------------------------------------
// Iterative radix-2 FFT.  Every transform here is a power of two by
// construction: the matvec pads to the circulant embedding length, and the
// preconditioner is deliberately built at the next power of two above n
// rather than at n itself, so no mixed-radix or Bluestein path is needed.
// ---------------------------------------------------------------------------
// Twiddle factors for one direction and transform length, split into real and
// imaginary planes.
//
// Computing them by the recurrence w *= step inside the butterfly is what a
// textbook FFT does, and it is doubly bad here: the recurrence serialises the
// inner loop so it cannot vectorise at all, and the repeated complex multiply
// accumulates rounding error that grows with the stage length.  Tabulating
// them once per size removes both, and the table is reused across every solve
// on the grid.
struct RLF_Twiddles {
  std::vector<double> re, im;
};

inline const RLF_Twiddles& rlf_twiddles(size_t n, bool inverse) {
  // Two cache slots, one per direction; a solve alternates forward/inverse at
  // a single size, so this never thrashes.
  static thread_local size_t cached_n[2] = {0, 0};
  static thread_local RLF_Twiddles cache[2];
  const int slot = inverse ? 1 : 0;
  RLF_Twiddles& tw = cache[slot];
  if (cached_n[slot] == n) return tw;

  // Stage `len` contributes len/2 factors, so the table holds
  // 1 + 2 + ... + n/2 = n - 1 entries laid out back to back.
  tw.re.resize(n);
  tw.im.resize(n);
  for (size_t len = 2, base = 0; len <= n; len <<= 1) {
    const size_t half = len / 2;
    const double theta =
      (inverse ? 2.0 : -2.0) * M_PI / static_cast<double>(len);
    for (size_t k = 0; k < half; ++k) {
      const double ang = theta * static_cast<double>(k);
      tw.re[base + k] = std::cos(ang);
      tw.im[base + k] = std::sin(ang);
    }
    base += half;
  }
  cached_n[slot] = n;
  return tw;
}

inline void rlf_fft(std::vector<rlf_cplx>& a, bool inverse) {
  const size_t n = a.size();
  if (n < 2) return;

  for (size_t i = 1, j = 0; i < n; ++i) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }

  // std::complex<double> is guaranteed to have the layout of double[2], so the
  // butterflies can run over raw pairs and vectorise.
  double* __restrict z = reinterpret_cast<double*>(a.data());
  const RLF_Twiddles& tw = rlf_twiddles(n, inverse);

  for (size_t len = 2, base = 0; len <= n; len <<= 1) {
    const size_t half = len / 2;
    const double* __restrict wr = tw.re.data() + base;
    const double* __restrict wi = tw.im.data() + base;
    for (size_t i = 0; i < n; i += len) {
      const size_t lo = 2 * i;
      const size_t hi = 2 * (i + half);
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
      for (size_t k = 0; k < half; ++k) {
        const double ur = z[lo + 2 * k];
        const double ui = z[lo + 2 * k + 1];
        const double xr = z[hi + 2 * k];
        const double xi = z[hi + 2 * k + 1];
        const double vr = xr * wr[k] - xi * wi[k];
        const double vi = xr * wi[k] + xi * wr[k];
        z[lo + 2 * k] = ur + vr;
        z[lo + 2 * k + 1] = ui + vi;
        z[hi + 2 * k] = ur - vr;
        z[hi + 2 * k + 1] = ui - vi;
      }
    }
    base += half;
  }

  if (inverse) {
    const double inv = 1.0 / static_cast<double>(n);
    const size_t total = 2 * n;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (size_t k = 0; k < total; ++k) z[k] *= inv;
  }
}

inline size_t rlf_next_pow2(size_t n) {
  size_t m = 1;
  while (m < n) m <<= 1;
  return m;
}

// ---------------------------------------------------------------------------
// M = I - theta*dt*L, held as a circulant embedding of the Toeplitz part plus
// an explicit diagonal, together with the eigenvalues of its Strang circulant
// preconditioner.  All three are rebuilt once per graded-time block.
// ---------------------------------------------------------------------------
struct RLF_ToeplitzSystem {
  int n = 0;
  size_t m = 0;                     // circulant embedding length (matvec)
  size_t mp = 0;                    // preconditioner circulant length
  std::vector<rlf_cplx> eig;        // FFT of the embedded first column
  std::vector<double> diag;         // diagonal of M
  std::vector<rlf_cplx> pre_eig;    // Strang preconditioner eigenvalues
  mutable std::vector<rlf_cplx> work;
  mutable std::vector<rlf_cplx> pwork;
};

// Build M = I - theta*dt*L directly from the dense generator.  Only the first
// row and first column are read, so this is O(n) rather than O(n^2): the
// interior of L is redundant once translation invariance is known.
inline RLF_ToeplitzSystem rlf_build_toeplitz_system(
    const std::vector<double>& L, int n, double theta_dt) {
  RLF_ToeplitzSystem sys;
  sys.n = n;
  sys.m = rlf_next_pow2(static_cast<size_t>(2 * n));

  // Row 0 (columns > 0) gives the super-diagonals; column 0 (rows > 0) gives
  // the sub-diagonals.  Entry (0,0) is excluded from both -- it carries the
  // sealed-face correction, which belongs to the diagonal.
  std::vector<rlf_cplx> col(sys.m, rlf_cplx(0.0, 0.0));
  for (int k = 1; k < n; ++k) {
    const double sub = -theta_dt * L[static_cast<size_t>(k) * n];        // (k,0)
    const double sup = -theta_dt * L[static_cast<size_t>(k)];            // (0,k)
    col[static_cast<size_t>(k)] = sub;
    col[sys.m - static_cast<size_t>(k)] = sup;
  }
  sys.eig = col;
  rlf_fft(sys.eig, false);

  sys.diag.resize(n);
  for (int j = 0; j < n; ++j) {
    sys.diag[j] = 1.0 - theta_dt * L[static_cast<size_t>(j) * n + j];
  }

  // Strang preconditioner.  Building it at the true size n would need an
  // arbitrary-length DFT, and Bluestein costs three padded transforms per
  // call against the matvec's one -- the preconditioner would dominate the
  // solve it is meant to accelerate.  Instead it is built at the next power
  // of two and applied to the zero-padded residual.  A preconditioner only
  // has to approximate M^-1, and the jump weights decay, so folding at N
  // rather than n clusters the spectrum just as well while keeping every
  // transform radix-2.
  double mean_diag = 0.0;
  for (int j = 0; j < n; ++j) mean_diag += sys.diag[j];
  mean_diag /= static_cast<double>(n);

  sys.mp = rlf_next_pow2(static_cast<size_t>(n));
  std::vector<rlf_cplx> pre(sys.mp, rlf_cplx(0.0, 0.0));
  pre[0] = mean_diag;
  // Stop short of mp/2 so the sub- and super-diagonal folds cannot land on
  // the same slot and overwrite each other.
  const int half = static_cast<int>((sys.mp - 1) / 2);
  const int reach = std::min(half, n - 1);
  for (int k = 1; k <= reach; ++k) {
    const double sub = -theta_dt * L[static_cast<size_t>(k) * n];
    const double sup = -theta_dt * L[static_cast<size_t>(k)];
    pre[static_cast<size_t>(k)] = sub;
    pre[sys.mp - static_cast<size_t>(k)] = sup;
  }
  sys.pre_eig = pre;
  rlf_fft(sys.pre_eig, false);

  sys.work.assign(sys.m, rlf_cplx(0.0, 0.0));
  sys.pwork.assign(sys.mp, rlf_cplx(0.0, 0.0));
  return sys;
}

// y = M x
inline void rlf_toeplitz_apply(const RLF_ToeplitzSystem& sys,
                               const double* x, double* y) {
  const int n = sys.n;
  std::fill(sys.work.begin(), sys.work.end(), rlf_cplx(0.0, 0.0));
  for (int i = 0; i < n; ++i) sys.work[static_cast<size_t>(i)] = x[i];

  rlf_fft(sys.work, false);
  // libstdc++'s complex operator* carries a NaN-rescue branch that defeats
  // vectorisation; over raw pairs this is a clean SIMD kernel.
  {
    double* __restrict w = reinterpret_cast<double*>(sys.work.data());
    const double* __restrict e =
      reinterpret_cast<const double*>(sys.eig.data());
    const size_t total = sys.m;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (size_t k = 0; k < total; ++k) {
      const double ar = w[2 * k], ai = w[2 * k + 1];
      const double br = e[2 * k], bi = e[2 * k + 1];
      w[2 * k] = ar * br - ai * bi;
      w[2 * k + 1] = ar * bi + ai * br;
    }
  }
  rlf_fft(sys.work, true);

  {
    const double* __restrict w =
      reinterpret_cast<const double*>(sys.work.data());
    const double* __restrict d = sys.diag.data();
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) y[i] = w[2 * i] + d[i] * x[i];
  }
}

// z = C^-1 r for the Strang circulant C.
inline void rlf_toeplitz_precondition(const RLF_ToeplitzSystem& sys,
                                      const double* r, double* z) {
  const int n = sys.n;
  std::fill(sys.pwork.begin(), sys.pwork.end(), rlf_cplx(0.0, 0.0));
  for (int i = 0; i < n; ++i) sys.pwork[static_cast<size_t>(i)] = r[i];

  rlf_fft(sys.pwork, false);
  for (size_t k = 0; k < sys.mp; ++k) {
    const rlf_cplx lam = sys.pre_eig[k];
    // A near-singular mode would amplify noise rather than precondition; the
    // preconditioner only has to be approximate, so such modes pass through.
    if (std::abs(lam) > 1e-13) sys.pwork[k] /= lam;
  }
  rlf_fft(sys.pwork, true);

  for (int i = 0; i < n; ++i) z[i] = sys.pwork[static_cast<size_t>(i)].real();
}

// Preconditioned BiCGSTAB.  Returns the iteration count, or -1 if the
// tolerance was not reached, which lets the caller fall back to the dense
// path rather than march on an unconverged stage.
inline int rlf_toeplitz_solve(const RLF_ToeplitzSystem& sys,
                              const double* b, double* x,
                              double tol = 1e-12, int max_iter = 200) {
  // The reductions and axpys below are the only non-FFT work in the solve;
  // -ffast-math plus ivdep lets them fold into SIMD lanes rather than running
  // as serial accumulations.
  const int n = sys.n;
  std::vector<double> rv(n), r0v_(n), pv(n), vv(n), sv(n), tv(n), phv(n), shv(n);
  double* __restrict r = rv.data();
  double* __restrict r0 = r0v_.data();
  double* __restrict p = pv.data();
  double* __restrict v = vv.data();
  double* __restrict s = sv.data();
  double* __restrict t = tv.data();
  double* __restrict ph = phv.data();
  double* __restrict sh = shv.data();

  rlf_toeplitz_apply(sys, x, v);
  double bnorm = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
  for (int i = 0; i < n; ++i) {
    r[i] = b[i] - v[i];
    r0[i] = r[i];
    p[i] = r[i];
    bnorm += b[i] * b[i];
  }
  bnorm = std::sqrt(bnorm);
  if (!(bnorm > 0.0)) {
    std::fill(x, x + n, 0.0);
    return 0;
  }
  const double target = tol * bnorm;

  double rho = 0.0;
  for (int i = 0; i < n; ++i) rho += r0[i] * r[i];

  for (int iter = 1; iter <= max_iter; ++iter) {
    rlf_toeplitz_precondition(sys, p, ph);
    rlf_toeplitz_apply(sys, ph, v);

    double r0v = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) r0v += r0[i] * v[i];
    if (!(std::abs(r0v) > 0.0)) return -1;
    const double alpha = rho / r0v;

    double snorm = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) {
      s[i] = r[i] - alpha * v[i];
      snorm += s[i] * s[i];
    }
    if (std::sqrt(snorm) < target) {
      for (int i = 0; i < n; ++i) x[i] += alpha * ph[i];
      return iter;
    }

    rlf_toeplitz_precondition(sys, s, sh);
    rlf_toeplitz_apply(sys, sh, t);

    double tt = 0.0, ts = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) {
      tt += t[i] * t[i];
      ts += t[i] * s[i];
    }
    if (!(tt > 0.0)) return -1;
    const double omega = ts / tt;

    double rnorm = 0.0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) {
      x[i] += alpha * ph[i] + omega * sh[i];
      r[i] = s[i] - omega * t[i];
      rnorm += r[i] * r[i];
    }
    if (std::sqrt(rnorm) < target) return iter;
    if (!(std::abs(omega) > 0.0)) return -1;

    double rho_new = 0.0;
    for (int i = 0; i < n; ++i) rho_new += r0[i] * r[i];
    if (!(std::abs(rho_new) > 0.0)) return -1;
    const double beta = (rho_new / rho) * (alpha / omega);
    rho = rho_new;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC ivdep
#endif
    for (int i = 0; i < n; ++i) {
      p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
  }
  return -1;
}

} // namespace rlf

#endif // rlf_toeplitz_h
