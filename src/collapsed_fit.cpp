// The collapsed move's quadratic fit of a subject's log-likelihood to its
// window draws (.emc_collapsed_fit_subject, R/collapsed_move.R), with the
// routines and sums R's own functions use.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#pragma GCC optimize("fp-contract=off")
#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cfloat>
#include <cmath>
#include <vector>

#ifndef FCONE
#define FCONE
#endif

namespace {

typedef std::vector<double> Vec;

bool all_finite(const Vec& x) {
  for (double v : x)
    if (!std::isfinite(v)) return false;
  return true;
}

// x %*% y for finite operands as R's product forms it.
bool matprod(const Vec& x, int nrx, int ncx, const Vec& y, int nry, int ncy, Vec& z) {
  if (ncx != nry || nrx < 1 || ncx < 1 || ncy < 1 || !all_finite(x) || !all_finite(y))
    return false;
  z.assign(static_cast<std::size_t>(nrx) * ncy, 0.0);
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  if (ncy == 1)
    F77_CALL(dgemv)("N", &nrx, &ncx, &one, x.data(), &nrx, y.data(), &ione, &zero, z.data(),
                    &ione FCONE);
  else if (nrx == 1)
    F77_CALL(dgemv)("T", &nry, &ncy, &one, y.data(), &nry, x.data(), &ione, &zero, z.data(),
                    &ione FCONE);
  else
    F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one, x.data(), &nrx, y.data(), &nry, &zero,
                    z.data(), &nrx FCONE FCONE);
  return true;
}

// eigen(x, symmetric = TRUE): values in decreasing order and their vectors.
bool eigen_sym(const Vec& x, int n, Vec& values, Vec& vectors) {
  Vec a(x), w(n), z(static_cast<std::size_t>(n) * n);
  std::vector<int> isuppz(2 * static_cast<std::size_t>(n));
  int m = 0, il = 0, iu = 0, info = 0, lwork = -1, liwork = -1, iw0 = 0;
  double vl = 0.0, vu = 0.0, abstol = 0.0, w0 = 0.0;
  F77_CALL(dsyevr)("V", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.data(), &n, isuppz.data(), &w0, &lwork, &iw0, &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) return false;
  lwork = static_cast<int>(w0);
  liwork = iw0;
  Vec work(lwork);
  std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)("V", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.data(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) return false;
  values.resize(n);
  vectors.resize(z.size());
  for (int c = 0; c < n; ++c) {
    const int src = n - 1 - c;
    values[c] = w[src];
    std::copy(z.begin() + static_cast<std::size_t>(src) * n,
              z.begin() + static_cast<std::size_t>(src + 1) * n,
              vectors.begin() + static_cast<std::size_t>(c) * n);
  }
  return true;
}

// qr.coef(qr(X), y) for a full-rank X (n x q): false where the rank is short
// or the solve fails.
bool qr_coef(Vec X, int n, int q, Vec y, Vec& coef) {
  Vec qraux(q), work(2 * static_cast<std::size_t>(q));
  std::vector<int> pivot(q);
  for (int j = 0; j < q; ++j) pivot[j] = j + 1;
  double tol = 1e-7;
  int rank = 0, ny = 1, info = 0;
  F77_CALL(dqrdc2)(X.data(), &n, &n, &q, &tol, &rank, qraux.data(), pivot.data(), work.data());
  if (rank != q) return false;
  coef.assign(q, NA_REAL);
  F77_CALL(dqrcf)(X.data(), &n, &rank, qraux.data(), y.data(), &ny, coef.data(), &info);
  return info == 0;
}

}  // namespace

// The fit: the precision and linear term of the quadratic in the subject's
// state, from draws A (rows) and their log-likelihoods y. NULL where the R
// form gives no fit or takes a path of its own.
// [[Rcpp::export(rng = false)]]
SEXP emc_collapsed_fit_subject_cpp(const Rcpp::NumericMatrix& A_in, const Rcpp::NumericVector& y_in) {
  const int n_all = A_in.nrow(), p = A_in.ncol();
  if (p < 1 || y_in.size() != n_all) return R_NilValue;
  std::vector<int> rows;
  for (int i = 0; i < n_all; ++i) {
    bool ok = R_FINITE(y_in[i]);
    for (int j = 0; ok && j < p; ++j) ok = !ISNAN(A_in(i, j));
    if (ok) rows.push_back(i);
  }
  const int n = static_cast<int>(rows.size());
  const int m = p * (p + 1) / 2, q = 1 + p + m;
  if (n < 2 * q) return R_NilValue;
  Vec A(static_cast<std::size_t>(n) * p), y(n);
  for (int i = 0; i < n; ++i) {
    y[i] = y_in[rows[i]];
    for (int j = 0; j < p; ++j) A[i + static_cast<std::size_t>(j) * n] = A_in(rows[i], j);
  }
  if (!all_finite(A)) return R_NilValue;
  // colMeans() and cov().
  Vec abar(p), xm(p), C(static_cast<std::size_t>(p) * p);
  for (int j = 0; j < p; ++j) {
    const double* x = &A[static_cast<std::size_t>(j) * n];
    long double sum = 0.0;
    for (int k = 0; k < n; ++k) sum += x[k];
    abar[j] = static_cast<double>(sum / n);
    long double tmp = sum / n;
    if (R_FINITE(static_cast<double>(tmp))) {
      sum = 0.0;
      for (int k = 0; k < n; ++k) sum += (x[k] - tmp);
      tmp = tmp + sum / n;
    }
    xm[j] = static_cast<double>(tmp);
  }
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      const double *xx = &A[static_cast<std::size_t>(i) * n], *yy = &A[static_cast<std::size_t>(j) * n];
      long double sum = 0.0;
      // The deviations and their product at extended precision.
      const long double xxm = xm[i], yym = xm[j];
      for (int k = 0; k < n; ++k) sum += (xx[k] - xxm) * (yy[k] - yym);
      const double v = static_cast<double>(sum / (n - 1));
      C[j + static_cast<std::size_t>(i) * p] = v;
      C[i + static_cast<std::size_t>(j) * p] = v;
    }
  if (!all_finite(C)) return R_NilValue;
  Vec ev, V;
  if (!eigen_sym(C, p, ev, V)) return R_NilValue;
  double top = R_NegInf;
  for (double v : ev) {
    if (!std::isfinite(v)) return R_NilValue;
    top = std::max(top, v);
  }
  volatile double edge = 1e-12 * top;
  for (double v : ev)
    if (v <= edge) return R_NilValue;
  // The whitening map W (rows are scaled eigenvectors) and whitened draws Z.
  Vec W(static_cast<std::size_t>(p) * p), Wt(W.size());
  for (int i = 0; i < p; ++i) {
    volatile double s = std::sqrt(ev[i]);
    for (int j = 0; j < p; ++j) {
      volatile double w = V[j + static_cast<std::size_t>(i) * p] / s;
      W[i + static_cast<std::size_t>(j) * p] = w;
      Wt[j + static_cast<std::size_t>(i) * p] = w;
    }
  }
  Vec centred(A.size()), Z;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i) {
      volatile double d = A[i + static_cast<std::size_t>(j) * n] - abar[j];
      centred[i + static_cast<std::size_t>(j) * n] = d;
    }
  if (!matprod(centred, n, p, Wt, p, p, Z)) return R_NilValue;
  // The quadratic regression: intercept, linear terms, and the products of
  // each pair of the upper triangle in column order.
  std::vector<int> iu1, iu2;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i <= j; ++i) {
      iu1.push_back(i);
      iu2.push_back(j);
    }
  Vec X(static_cast<std::size_t>(n) * q);
  for (int i = 0; i < n; ++i) X[i] = 1.0;
  std::copy(Z.begin(), Z.end(), X.begin() + n);
  for (int c = 0; c < m; ++c)
    for (int i = 0; i < n; ++i) {
      volatile double prod = Z[i + static_cast<std::size_t>(iu1[c]) * n] *
                             Z[i + static_cast<std::size_t>(iu2[c]) * n];
      X[i + static_cast<std::size_t>(1 + p + c) * n] = prod;
    }
  Vec b;
  if (!qr_coef(X, n, q, y, b)) return R_NilValue;
  for (double v : b)
    if (ISNAN(v)) return R_NilValue;
  // The curvature: a product's coefficient is minus its entry off the
  // diagonal and minus half of it on the diagonal.
  Vec B(static_cast<std::size_t>(p) * p, 0.0);
  for (int c = 0; c < m; ++c) {
    const double cq = b[1 + p + c];
    if (iu1[c] == iu2[c]) {
      volatile double d = -2 * cq;
      B[iu1[c] + static_cast<std::size_t>(iu2[c]) * p] = d;
    } else {
      volatile double up = -cq + 0.0;
      volatile double low = 0.0 + -cq;
      B[iu1[c] + static_cast<std::size_t>(iu2[c]) * p] = up;
      B[iu2[c] + static_cast<std::size_t>(iu1[c]) * p] = low;
    }
  }
  if (!all_finite(B)) return R_NilValue;
  Vec bv, BV;
  if (!eigen_sym(B, p, bv, BV)) return R_NilValue;
  // Its positive semi-definite part.
  Vec right(B.size()), Bp;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) {
      volatile double v = std::max(bv[i], 0.0) * BV[j + static_cast<std::size_t>(i) * p];
      right[i + static_cast<std::size_t>(j) * p] = v;
    }
  if (!matprod(BV, p, p, right, p, p, Bp)) return R_NilValue;
  // The linear term refitted with the curvature held.
  Vec ZB;
  if (!matprod(Z, n, p, Bp, p, p, ZB)) return R_NilValue;
  Vec rhs(n), X2(static_cast<std::size_t>(n) * (p + 1));
  {
    std::vector<long double> rs(n, 0.0L);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) {
        volatile double prod = ZB[i + static_cast<std::size_t>(j) * n] * Z[i + static_cast<std::size_t>(j) * n];
        rs[i] += prod;
      }
    for (int i = 0; i < n; ++i) {
      volatile double half = 0.5 * static_cast<double>(rs[i]);
      volatile double v = y[i] + half;
      rhs[i] = v;
    }
  }
  for (int i = 0; i < n; ++i) X2[i] = 1.0;
  std::copy(Z.begin(), Z.end(), X2.begin() + n);
  Vec b2;
  if (!qr_coef(X2, n, p + 1, rhs, b2)) return R_NilValue;
  Vec g(b2.begin() + 1, b2.end());
  for (double v : g)
    if (ISNAN(v)) return R_NilValue;
  // Slopes only along the directions the curvature keeps.
  double biggest = DBL_MIN;
  for (double v : bv) biggest = std::max(biggest, std::fabs(v));
  volatile double tol = biggest * 1e-10;
  std::vector<int> keep;
  for (int i = 0; i < p; ++i)
    if (bv[i] > tol) keep.push_back(i);
  if (!keep.empty()) {
    const int kk = static_cast<int>(keep.size());
    Vec U(static_cast<std::size_t>(p) * kk), Ug(kk), back;
    for (int c = 0; c < kk; ++c)
      std::copy(BV.begin() + static_cast<std::size_t>(keep[c]) * p,
                BV.begin() + static_cast<std::size_t>(keep[c] + 1) * p,
                U.begin() + static_cast<std::size_t>(c) * p);
    if (!all_finite(U) || !all_finite(g)) return R_NilValue;
    const double one = 1.0, zero = 0.0;
    const int ione = 1, pp = p;
    F77_CALL(dgemv)("T", &pp, &kk, &one, U.data(), &pp, g.data(), &ione, &zero, Ug.data(),
                    &ione FCONE);
    if (!matprod(U, p, kk, Ug, kk, 1, back)) return R_NilValue;
    g = back;
  } else {
    std::fill(g.begin(), g.end(), 0.0);
  }
  Vec WtB, H, Ha, Wg;
  if (!matprod(Wt, p, p, Bp, p, p, WtB) || !matprod(WtB, p, p, W, p, p, H)) return R_NilValue;
  {
    Vec H0(H);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < p; ++i) {
        volatile double s = H0[i + static_cast<std::size_t>(j) * p] + H0[j + static_cast<std::size_t>(i) * p];
        H[i + static_cast<std::size_t>(j) * p] = s / 2;
      }
  }
  if (!matprod(H, p, p, abar, p, 1, Ha) || !matprod(Wt, p, p, g, p, 1, Wg)) return R_NilValue;
  Rcpp::NumericMatrix prec(p, p);
  Rcpp::NumericVector lin(p);
  std::copy(H.begin(), H.end(), prec.begin());
  for (int i = 0; i < p; ++i) {
    volatile double v = Ha[i] + Wg[i];
    lin[i] = v;
    if (!std::isfinite(lin[i])) return R_NilValue;
  }
  if (!all_finite(H)) return R_NilValue;
  return Rcpp::List::create(Rcpp::_["prec"] = prec, Rcpp::_["lin"] = lin);
}

namespace {

// LINPACK's estimate of an upper triangular matrix's reciprocal condition in
// the one norm (dtrco), which kappa() takes for a QR's triangular factor.
double upper_rcond(const Vec& t, int n) {
  auto T = [&](int i, int j) { return t[i + static_cast<std::size_t>(j) * n]; };
  const int ione = 1;
  double tnorm = 0.0;
  for (int j = 0; j < n; ++j) {
    const int l = j + 1;
    tnorm = std::max(tnorm, F77_CALL(dasum)(&l, &t[static_cast<std::size_t>(j) * n], &ione));
  }
  // The transposed system for a right-hand side of growing signs.
  double ek = 1.0;
  Vec z(n, 0.0);
  for (int k = 0; k < n; ++k) {
    if (z[k] != 0.0) ek = std::copysign(std::fabs(ek), -z[k]);
    if (std::fabs(ek - z[k]) > std::fabs(T(k, k))) {
      double s = std::fabs(T(k, k)) / std::fabs(ek - z[k]);
      F77_CALL(dscal)(&n, &s, z.data(), &ione);
      ek = s * ek;
    }
    double wk = ek - z[k], wkm = -ek - z[k];
    double s = std::fabs(wk), sm = std::fabs(wkm);
    if (T(k, k) != 0.0) {
      wk = wk / T(k, k);
      wkm = wkm / T(k, k);
    } else {
      wk = 1.0;
      wkm = 1.0;
    }
    if (k != n - 1) {
      for (int j = k + 1; j < n; ++j) {
        sm = sm + std::fabs(z[j] + wkm * T(k, j));
        z[j] = z[j] + wk * T(k, j);
        s = s + std::fabs(z[j]);
      }
      if (s < sm) {
        const double w = wkm - wk;
        wk = wkm;
        for (int j = k + 1; j < n; ++j) z[j] = z[j] + w * T(k, j);
      }
    }
    z[k] = wk;
  }
  double s = 1.0 / F77_CALL(dasum)(&n, z.data(), &ione);
  F77_CALL(dscal)(&n, &s, z.data(), &ione);
  double ynorm = 1.0;
  // The system itself.
  for (int kk = 1; kk <= n; ++kk) {
    const int k = n - kk;
    if (std::fabs(z[k]) > std::fabs(T(k, k))) {
      s = std::fabs(T(k, k)) / std::fabs(z[k]);
      F77_CALL(dscal)(&n, &s, z.data(), &ione);
      ynorm = s * ynorm;
    }
    if (T(k, k) != 0.0) z[k] = z[k] / T(k, k);
    if (T(k, k) == 0.0) z[k] = 1.0;
    if (kk < n) {
      const int len = n - kk;
      double w = -z[k];
      F77_CALL(daxpy)(&len, &w, &t[static_cast<std::size_t>(k) * n], &ione, z.data(), &ione);
    }
  }
  s = 1.0 / F77_CALL(dasum)(&n, z.data(), &ione);
  F77_CALL(dscal)(&n, &s, z.data(), &ione);
  ynorm = s * ynorm;
  return tnorm != 0.0 ? ynorm / tnorm : 0.0;
}

}  // namespace

// .emc_empirical_covariance_ok(): whether a square matrix is finite, exactly
// symmetric, has no eigenvalue under 1e-8 and a reciprocal condition (from
// the triangular factor of its QR, as kappa() estimates it) of at least
// 1e-12. NULL where a decomposition fails, which the R form reports.
// [[Rcpp::export(rng = false)]]
SEXP emc_covariance_ok_cpp(const Rcpp::NumericMatrix& S) {
  int n = S.nrow();
  if (n < 1 || S.ncol() != n) return Rf_ScalarLogical(FALSE);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      if (!R_FINITE(S(i, j)) || S(i, j) != S(j, i)) return Rf_ScalarLogical(FALSE);
  // The eigenvalues alone.
  {
    Vec a(S.begin(), S.end()), w(n);
    std::vector<int> isuppz(2 * static_cast<std::size_t>(n));
    int m = 0, il = 0, iu = 0, info = 0, lwork = -1, liwork = -1, iw0 = 0;
    double vl = 0.0, vu = 0.0, abstol = 0.0, w0 = 0.0, z = 0.0;
    F77_CALL(dsyevr)("N", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                     &z, &n, isuppz.data(), &w0, &lwork, &iw0, &liwork, &info FCONE FCONE FCONE);
    if (info != 0) return R_NilValue;
    lwork = static_cast<int>(w0);
    liwork = iw0;
    Vec work(lwork);
    std::vector<int> iwork(liwork);
    F77_CALL(dsyevr)("N", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                     &z, &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                     &info FCONE FCONE FCONE);
    if (info != 0) return R_NilValue;
    for (double v : w)
      if (!R_FINITE(v) || v < 1e-8) return Rf_ScalarLogical(FALSE);
  }
  // The condition estimate of the QR's triangular factor in the one norm.
  Vec X(S.begin(), S.end()), qraux(n), work(2 * static_cast<std::size_t>(n));
  std::vector<int> pivot(n);
  for (int j = 0; j < n; ++j) pivot[j] = j + 1;
  double tol = 1e-7;
  int rank = 0, info = 0;
  F77_CALL(dqrdc2)(X.data(), &n, &n, &n, &tol, &rank, qraux.data(), pivot.data(), work.data());
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) X[i + static_cast<std::size_t>(j) * n] = 0.0;
  volatile double rc = upper_rcond(X, n);
  volatile double kappa = 1 / rc;
  volatile double rcond = 1 / kappa;
  return Rf_ScalarLogical(R_FINITE(rcond) && rcond >= 1e-12);
}
