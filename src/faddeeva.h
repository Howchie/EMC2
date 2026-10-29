#ifndef EMC2_FADDEEVA_H
#define EMC2_FADDEEVA_H

// Faddeeva function w(z) = exp(-z^2) erfc(-i z) in the closed upper half
// plane, by Weideman's rational expansion (SIAM J. Numer. Anal. 31, 1994)
// with N = 40 terms: relative error below ~1e-12 in w and in Re w over
// Im z >= 0, including the far field and the region next to the real axis.
// The expansion coefficients are computed once, on first use.

#include <cmath>
#include <complex>

namespace faddeeva {

constexpr int WEIDEMAN_N = 40;

struct WeidemanCoefs {
  double L;
  double a[WEIDEMAN_N];  // a[n - 1] multiplies Z^(n - 1)
  WeidemanCoefs() {
    const int N = WEIDEMAN_N, M = 2 * N;
    L = std::sqrt(N / std::sqrt(2.0));
    for (int n = 1; n <= N; ++n) {
      double acc = 0.0;
      for (int k = -M + 1; k <= M - 1; ++k) {
        const double t = L * std::tan(0.5 * k * M_PI / M);
        acc += std::exp(-t * t) * (L * L + t * t) * std::cos(M_PI * k * n / M);
      }
      a[n - 1] = acc / (2.0 * M);
    }
  }
};

inline const WeidemanCoefs& weideman_coefs() {
  static const WeidemanCoefs c;
  return c;
}

// Requires Im z >= 0.
inline std::complex<double> w(std::complex<double> z) {
  const WeidemanCoefs& c = weideman_coefs();
  const std::complex<double> iz(-z.imag(), z.real());
  const std::complex<double> den = c.L - iz;
  const std::complex<double> Z = (c.L + iz) / den;
  std::complex<double> p = c.a[WEIDEMAN_N - 1];
  for (int n = WEIDEMAN_N - 2; n >= 0; --n) p = p * Z + c.a[n];
  return 2.0 * p / (den * den) + 1.0 / (std::sqrt(M_PI) * den);
}

// w at two points in one pass: the Horner chains are interleaved, so their
// latencies overlap.  Same arithmetic as w().
inline void w2(std::complex<double> z1, std::complex<double> z2,
               std::complex<double>& out1, std::complex<double>& out2) {
  const WeidemanCoefs& c = weideman_coefs();
  const std::complex<double> iz1(-z1.imag(), z1.real()), iz2(-z2.imag(), z2.real());
  const std::complex<double> den1 = c.L - iz1, den2 = c.L - iz2;
  const std::complex<double> Z1 = (c.L + iz1) / den1, Z2 = (c.L + iz2) / den2;
  std::complex<double> p1 = c.a[WEIDEMAN_N - 1], p2 = p1;
  for (int n = WEIDEMAN_N - 2; n >= 0; --n) {
    p1 = p1 * Z1 + c.a[n];
    p2 = p2 * Z2 + c.a[n];
  }
  out1 = 2.0 * p1 / (den1 * den1) + 1.0 / (std::sqrt(M_PI) * den1);
  out2 = 2.0 * p2 / (den2 * den2) + 1.0 / (std::sqrt(M_PI) * den2);
}

}  // namespace faddeeva

#endif
