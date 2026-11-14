// The group draws of the factor hierarchies: parameter-expanded factor
// analysis with loadings fixed or free by entry (R/variant_factor.R
// gibbs_step_factor) and the shrinkage prior over an unbounded number of
// factors (R/variant_infnt_factor.R gibbs_step_infnt_factor), and the structural equation model's draw up to
// its structural coefficients (R/variant_SEM.R).
// Every product, inverse and eigendecomposition goes through the routine R's
// own operators call, in R's order, so the draw equals the R form's bit for
// bit on the same random stream. The draw declines where the R form would
// raise a condition or take a path this file does not mirror. A product
// feeding a sum is rounded first, as in R's arithmetic.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#pragma GCC optimize ("fp-contract=off")
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
#include <cfloat>
#include <cmath>
#include <vector>
#include <cstring>

#include "gibbs_factor.h"

namespace {

using emc_factor::Mat;

// The condition under which a step leaves the draw to the R form.
struct Decline {};

void require_finite(const Mat& x) {
  for (double v : x.v)
    if (!R_FINITE(v)) throw Decline();
}

Mat transpose(const Mat& x) {
  Mat z(x.c, x.r);
  for (int j = 0; j < x.c; ++j)
    for (int i = 0; i < x.r; ++i) z(j, i) = x(i, j);
  return z;
}

Mat diagonal(const std::vector<double>& d) {
  const int n = static_cast<int>(d.size());
  Mat z(n, n);
  for (int i = 0; i < n; ++i) z(i, i) = d[i];
  return z;
}

// x %*% y as R computes it for finite operands: a matrix-vector product when
// either result dimension is one, else the general product.
Mat mult(const Mat& x, const Mat& y) {
  if (x.c != y.r) throw Decline();
  // A product over an empty dimension is the zero matrix.
  if (x.r == 0 || x.c == 0 || y.c == 0) return Mat(x.r, y.c);
  require_finite(x);
  require_finite(y);
  Mat z(x.r, y.c);
  int nrx = x.r, ncx = x.c, nry = y.r, ncy = y.c, ione = 1;
  double one = 1.0, zero = 0.0;
  if (ncy == 1) {
    F77_CALL(dgemv)("N", &nrx, &ncx, &one, x.v.data(), &nrx, y.v.data(), &ione, &zero,
                    z.v.data(), &ione FCONE);
  } else if (nrx == 1) {
    F77_CALL(dgemv)("T", &nry, &ncy, &one, y.v.data(), &nry, x.v.data(), &ione, &zero,
                    z.v.data(), &ione FCONE);
  } else {
    F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one, x.v.data(), &nrx, y.v.data(), &nry,
                    &zero, z.v.data(), &nrx FCONE FCONE);
  }
  return z;
}

// solve(a): the inverse from the LU solve against the identity, refused
// where R reports a singular or ill-conditioned matrix.
Mat inverse(const Mat& a) {
  if (a.r != a.c || a.r == 0) throw Decline();
  int n = a.r, info = 0;
  Mat lu = a, b(n, n);
  for (int i = 0; i < n; ++i) b(i, i) = 1.0;
  std::vector<int> ipiv(n);
  F77_CALL(dgesv)(&n, &n, lu.v.data(), &n, ipiv.data(), b.v.data(), &n, &info);
  if (info != 0) throw Decline();
  double anorm = F77_CALL(dlange)("1", &n, &n, a.v.data(), &n, NULL FCONE);
  std::vector<double> work(4 * static_cast<size_t>(n));
  double rcond = 0.0;
  F77_CALL(dgecon)("1", &n, lu.v.data(), &n, &anorm, &rcond, work.data(), ipiv.data(),
                   &info FCONE);
  if (!(rcond >= DBL_EPSILON)) throw Decline();
  return b;
}

// all.equal()'s mean relative difference of two equal-length vectors, zero
// when they are equal.
double mean_difference(const std::vector<double>& target, const std::vector<double>& current,
                       double tolerance) {
  long double diff = 0.0, scale = 0.0;
  size_t n = 0;
  for (size_t i = 0; i < target.size(); ++i) {
    if (target[i] == current[i]) continue;
    diff += std::fabs(target[i] - current[i]);
    scale += std::fabs(target[i]);
    ++n;
  }
  if (n == 0) return 0.0;
  double xy = static_cast<double>(diff) / n;
  const double sc = static_cast<double>(scale) / n;
  if (R_FINITE(sc) && sc > tolerance) xy = xy / sc;
  return xy;
}

// Whether isSymmetric(x, tol) holds with a margin: a matrix near the
// tolerance is left to the R form.
void require_symmetric(const Mat& x) {
  const double tol = std::sqrt(DBL_EPSILON);
  const int n = x.r;
  std::vector<double> a, b;
  if (n > 1) {
    const int rows[4] = {0, 1, n - 2, n - 1};
    for (int i : rows) {
      a.clear();
      b.clear();
      for (int j = 0; j < n; ++j) {
        a.push_back(x(i, j));
        b.push_back(x(j, i));
      }
      if (!(mean_difference(a, b, 8 * tol) <= 4 * tol)) throw Decline();
    }
  }
  const Mat xt = transpose(x);
  if (!(mean_difference(x.v, xt.v, tol) <= 0.5 * tol)) throw Decline();
}

// The factor rmvnorm() multiplies standard normal rows by: from the
// eigendecomposition of the symmetric covariance, largest eigenvalue first.
Mat normal_factor(const Mat& sigma) {
  require_symmetric(sigma);
  require_finite(sigma);
  int n = sigma.r, info = 0, m = 0, il = 0, iu = 0, lwork = -1, liwork = -1, itmp = 0;
  double vl = 0.0, vu = 0.0, abstol = 0.0, tmp = 0.0;
  Mat a = sigma, z(n, n);
  std::vector<double> w(n);
  std::vector<int> isuppz(2 * static_cast<size_t>(n));
  F77_CALL(dsyevr)("V", "A", "L", &n, a.v.data(), &n, &vl, &vu, &il, &iu, &abstol, &m,
                   w.data(), z.v.data(), &n, isuppz.data(), &tmp, &lwork, &itmp, &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) throw Decline();
  lwork = static_cast<int>(tmp);
  liwork = itmp;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)("V", "A", "L", &n, a.v.data(), &n, &vl, &vu, &il, &iu, &abstol, &m,
                   w.data(), z.v.data(), &n, isuppz.data(), work.data(), &lwork, iwork.data(),
                   &liwork, &info FCONE FCONE FCONE);
  if (info != 0) throw Decline();
  std::vector<double> values(n);
  Mat vectors(n, n);
  for (int k = 0; k < n; ++k) {
    values[k] = w[n - 1 - k];
    for (int i = 0; i < n; ++i) vectors(i, k) = z(i, n - 1 - k);
  }
  const double floor = -std::sqrt(DBL_EPSILON) * std::fabs(values[0]);
  for (int k = 0; k < n; ++k)
    if (!(values[k] >= floor)) throw Decline();
  Mat scaled(n, n);
  for (int i = 0; i < n; ++i) {
    const double s = std::sqrt(values[i] > 0.0 ? values[i] : 0.0);
    for (int j = 0; j < n; ++j) scaled(i, j) = vectors(j, i) * s;
  }
  return transpose(mult(vectors, scaled));
}

// rmvnorm(n, mean, sigma) given sigma's factor: standard normals drawn row
// by row. A NULL mean is the zero vector.
Mat normal_rows(int n, const double* mean, const Mat& factor) {
  const int k = factor.r;
  Mat z(n, k);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < k; ++j) z(i, j) = norm_rand();
  Mat out = mult(z, factor);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < n; ++i) out(i, j) = out(i, j) + (mean == NULL ? 0.0 : mean[j]);
  return out;
}

// colSums().
std::vector<double> col_sums(const Mat& x) {
  std::vector<double> out(x.c);
  for (int j = 0; j < x.c; ++j) {
    long double s = 0.0;
    for (int i = 0; i < x.r; ++i) s += x(i, j);
    out[j] = static_cast<double>(s);
  }
  return out;
}

// rgamma(shape, rate = rate).
double gamma_draw(double shape, double rate) {
  const double scale = 1.0 / rate;
  const double x = Rf_rgamma(shape, scale);
  if (!R_FINITE(x)) throw Decline();
  return x;
}

void draw(emc_factor::State& s, const Mat& alpha, const emc_factor::Prior& prior,
          const std::vector<int>& free_entry) {
  const int n = alpha.r, p = alpha.c, k = s.eta.c;
  const Mat sig = diagonal(s.sig_err_inv), psi = diagonal(s.psi_inv);

  // The factor scores' conditional covariance and its normal factor.
  auto score_cov = [&](const Mat& lambda) {
    Mat m = mult(mult(transpose(lambda), sig), lambda);
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < k; ++i) m(i, j) = psi(i, j) + m(i, j);
    return inverse(m);
  };
  auto centred = [&](const std::vector<double>& mu) {
    Mat c(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) c(i, j) = alpha(i, j) - mu[j];
    return c;
  };
  auto scores = [&](const Mat& eta_sig, const Mat& resid) {
    const Mat mean = mult(mult(mult(eta_sig, transpose(s.lambda)), sig), transpose(resid));
    const Mat noise = normal_rows(n, NULL, normal_factor(eta_sig));
    for (int f = 0; f < k; ++f)
      for (int i = 0; i < n; ++i) s.eta(i, f) = mean(f, i) + noise(i, f);
  };

  // The scores given the subjects' states, before anything conditions on
  // them.
  const Mat eta_sig = score_cov(s.lambda);
  scores(eta_sig, centred(s.mu));

  // The group means.
  Mat prec(p, p);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) {
      volatile double scaled = n * sig(i, j);
      prec(i, j) = scaled + prior.mu_invar(i, j);
    }
  const Mat mu_sig = inverse(prec);
  {
    const Mat fitted = mult(s.eta, transpose(s.lambda));
    Mat resid(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) resid(i, j) = alpha(i, j) - fitted(i, j);
    Mat sums(p, 1), mean0(p, 1);
    sums.v = col_sums(resid);
    mean0.v = prior.mu_mean;
    const Mat a = mult(sig, sums), b = mult(prior.mu_invar, mean0);
    Mat rhs(p, 1);
    for (int i = 0; i < p; ++i) rhs(i, 0) = a(i, 0) + b(i, 0);
    const Mat mu_mu = mult(mu_sig, rhs);
    const Mat mu = normal_rows(1, mu_mu.v.data(), normal_factor(mu_sig));
    s.mu = mu.v;
  }
  const Mat tilde = centred(s.mu);

  // The scores given the new means.
  scores(score_cov(s.lambda), tilde);

  // The error precisions.
  {
    const Mat fitted = mult(s.eta, transpose(s.lambda));
    Mat sq(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) {
        const double d = tilde(i, j) - fitted(i, j);
        volatile double d2 = d * d;
        sq(i, j) = d2;
      }
    const std::vector<double> ss = col_sums(sq);
    for (int j = 0; j < p; ++j)
      s.sig_err_inv[j] = gamma_draw(prior.shape_err[j], prior.rate_err[j] + ss[j] / 2);
  }

  // The free loadings of each parameter, a regression on the scores net of
  // the fixed loadings.
  for (int j = 0; j < p; ++j) {
    std::vector<int> fr, fx;
    for (int f = 0; f < k; ++f) (free_entry[j + static_cast<size_t>(f) * p] ? fr : fx).push_back(f);
    if (fr.empty()) continue;
    const int kf = static_cast<int>(fr.size()), kx = static_cast<int>(fx.size());
    Mat eta_free(n, kf), y(n, 1);
    for (int f = 0; f < kf; ++f)
      for (int i = 0; i < n; ++i) eta_free(i, f) = s.eta(i, fr[f]);
    for (int i = 0; i < n; ++i) y(i, 0) = tilde(i, j);
    if (kx > 0) {
      Mat eta_fixed(n, kx), load(kx, 1);
      for (int f = 0; f < kx; ++f) {
        load(f, 0) = s.lambda(j, fx[f]);
        for (int i = 0; i < n; ++i) eta_fixed(i, f) = s.eta(i, fx[f]);
      }
      const Mat part = mult(eta_fixed, load);
      for (int i = 0; i < n; ++i) y(i, 0) = y(i, 0) - part(i, 0);
    }
    const Mat eta_free_t = transpose(eta_free);
    Mat m = mult(eta_free_t, eta_free);
    const double sj = s.sig_err_inv[j];
    for (int b = 0; b < kf; ++b)
      for (int a = 0; a < kf; ++a) {
        volatile double scaled = sj * m(a, b);
        m(a, b) = scaled + (a == b ? prior.lambda_invar[j] : 0.0);
      }
    const Mat lambda_sig = inverse(m);
    Mat scaled(kf, kf);
    for (size_t e = 0; e < scaled.v.size(); ++e) scaled.v[e] = lambda_sig.v[e] * sj;
    const Mat lambda_mu = mult(scaled, mult(eta_free_t, y));
    const Mat row = normal_rows(1, lambda_mu.v.data(), normal_factor(lambda_sig));
    for (int f = 0; f < kf; ++f) s.lambda(j, fr[f]) = row(0, f);
  }

  // The factor precisions.
  {
    Mat sq(n, k);
    for (size_t e = 0; e < sq.v.size(); ++e) {
      volatile double e2 = s.eta.v[e] * s.eta.v[e];
      sq.v[e] = e2;
    }
    const std::vector<double> ss = col_sums(sq);
    for (int f = 0; f < k; ++f)
      s.psi_inv[f] = gamma_draw(prior.shape_psi[f], prior.rate_psi[f] + ss[f] / 2);
  }

  // The implied group covariance and the loadings on unit-variance factors.
  s.var = mult(mult(s.lambda, inverse(diagonal(s.psi_inv))), transpose(s.lambda));
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i)
      s.var(i, j) = s.var(i, j) + (i == j ? 1.0 / s.sig_err_inv[i] : 0.0);
  std::vector<double> sd(k);
  for (int f = 0; f < k; ++f) sd[f] = std::sqrt(1.0 / s.psi_inv[f]);
  s.lambda_scaled = mult(s.lambda, diagonal(sd));
}

// chol(): the upper factor, refused where the matrix is not positive
// definite.
Mat chol_upper(const Mat& x) {
  if (x.r != x.c || x.r == 0) throw Decline();
  int n = x.r, info = 0;
  Mat u = x;
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) u(i, j) = 0.0;
  F77_CALL(dpotrf)("U", &n, u.v.data(), &n, &info FCONE);
  if (info != 0) throw Decline();
  return u;
}

// backsolve() and forwardsolve() on the triangle of `r`.
Mat tri_solve(const Mat& r, const Mat& b, bool upper) {
  if (r.r != r.c || b.r != r.r) throw Decline();
  int k = r.r, ncb = b.c;
  for (int i = 0; i < k; ++i)
    if (r(i, i) == 0.0) throw Decline();
  Mat z = b;
  double one = 1.0;
  F77_CALL(dtrsm)("L", upper ? "U" : "L", "N", "N", &k, &ncb, &one, r.v.data(), &k,
                  z.v.data(), &k FCONE FCONE FCONE FCONE);
  return z;
}

// qr.R(qr(x)) of a square matrix.
Mat qr_upper(const Mat& x) {
  if (x.r != x.c || x.r == 0) throw Decline();
  int n = x.r, p = x.c, rank = 0;
  double tol = 1e-07;
  Mat q = x;
  std::vector<double> qraux(p), work(2 * static_cast<size_t>(p));
  std::vector<int> pivot(p);
  for (int j = 0; j < p; ++j) pivot[j] = j + 1;
  F77_CALL(dqrdc2)(q.v.data(), &n, &n, &p, &tol, &rank, qraux.data(), pivot.data(),
                   work.data());
  for (int j = 0; j < p; ++j)
    for (int i = j + 1; i < n; ++i) q(i, j) = 0.0;
  return q;
}

// MASS::ginv() of a square matrix whose singular values all count.
Mat pseudo_inverse(const Mat& x) {
  if (x.r != x.c || x.r == 0) throw Decline();
  require_finite(x);
  int n = x.r, info = 0, lwork = -1;
  Mat a = x, u(n, n), vt(n, n);
  std::vector<double> d(n);
  std::vector<int> iwork(8 * static_cast<size_t>(n));
  double tmp = 0.0;
  F77_CALL(dgesdd)("S", &n, &n, a.v.data(), &n, d.data(), u.v.data(), &n, vt.v.data(), &n,
                   &tmp, &lwork, iwork.data(), &info FCONE);
  if (info != 0) throw Decline();
  lwork = static_cast<int>(tmp);
  std::vector<double> work(lwork);
  F77_CALL(dgesdd)("S", &n, &n, a.v.data(), &n, d.data(), u.v.data(), &n, vt.v.data(), &n,
                   work.data(), &lwork, iwork.data(), &info FCONE);
  if (info != 0) throw Decline();
  const double tiny = std::sqrt(DBL_EPSILON) * d[0];
  for (int i = 0; i < n; ++i)
    if (!(d[i] > (tiny > 0.0 ? tiny : 0.0))) throw Decline();
  Mat w(n, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) w(i, j) = (1.0 / d[i]) * u(j, i);
  return mult(transpose(vt), w);
}

// crossprod(x) (`of_rows` false) and tcrossprod(x).
Mat sym_product(const Mat& x, bool of_rows) {
  if (x.r == 0 || x.c == 0) throw Decline();
  require_finite(x);
  int nr = x.r, nc = x.c, n = of_rows ? nr : nc;
  double one = 1.0, zero = 0.0;
  Mat z(n, n);
  if (of_rows)
    F77_CALL(dsyrk)("U", "N", &nr, &nc, &one, x.v.data(), &nr, &zero, z.v.data(), &nr FCONE FCONE);
  else
    F77_CALL(dsyrk)("U", "T", &nc, &nr, &one, x.v.data(), &nr, &zero, z.v.data(), &nc FCONE FCONE);
  for (int i = 1; i < n; ++i)
    for (int j = 0; j < i; ++j) z(i, j) = z(j, i);
  return z;
}

// cumprod().
std::vector<double> cum_prod(const std::vector<double>& x) {
  std::vector<double> out(x.size());
  long double prod = 1.0;
  for (size_t i = 0; i < x.size(); ++i) {
    prod *= x[i];
    out[i] = static_cast<double>(prod);
  }
  return out;
}

void draw_infinite(emc_factor::InfiniteState& s, const Mat& alpha,
                   const emc_factor::InfinitePrior& prior) {
  const int n = alpha.r, p = alpha.c, K = s.lambda.c;
  std::vector<double> tauh = cum_prod(s.delta);
  const Mat plam = mult(s.psi, diagonal(tauh));
  const Mat weighted = mult(diagonal(s.epsilon_inv), s.lambda);

  // The scores' conditional covariance S S' through the factor S.
  Mat v1 = mult(transpose(weighted), s.lambda);
  for (int j = 0; j < K; ++j)
    for (int i = 0; i < K; ++i) v1(i, j) = (i == j ? 1.0 : 0.0) + v1(i, j);
  const Mat S = pseudo_inverse(qr_upper(chol_upper(v1)));
  const Mat cov = sym_product(S, true), St = transpose(S);
  auto centred = [&](const std::vector<double>& mu) {
    Mat c(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) c(i, j) = alpha(i, j) - mu[j];
    return c;
  };
  auto scores = [&](const Mat& resid) {
    const Mat mean = mult(mult(resid, weighted), cov);
    Mat z(n, K);
    for (double& e : z.v) e = norm_rand();
    const Mat noise = mult(z, St);
    for (size_t e = 0; e < s.eta.v.size(); ++e) s.eta.v[e] = mean.v[e] + noise.v[e];
  };
  auto residual = [&](const Mat& base) {
    const Mat fitted = mult(s.eta, transpose(s.lambda));
    Mat r(n, p);
    for (size_t e = 0; e < r.v.size(); ++e) r.v[e] = base.v[e] - fitted.v[e];
    return r;
  };

  // The scores given the subjects' states, before the means condition on
  // them.
  s.eta = Mat(n, K);
  scores(centred(s.mu));

  // The group means, independent given the scores.
  {
    const std::vector<double> sums = col_sums(residual(alpha));
    std::vector<double> mu_mu(p), mu_sig(p);
    for (int j = 0; j < p; ++j) {
      volatile double precision = n * s.epsilon_inv[j];
      mu_sig[j] = 1.0 / (precision + prior.mu_invar[j]);
      volatile double data = s.epsilon_inv[j] * sums[j];
      volatile double prior_part = prior.mu_invar[j] * prior.mu_mean[j];
      volatile double total = data + prior_part;
      mu_mu[j] = mu_sig[j] * total;
    }
    const Mat mu = normal_rows(1, mu_mu.data(), normal_factor(diagonal(mu_sig)));
    s.mu = mu.v;
  }
  const Mat tilde = centred(s.mu);

  // The scores given the new means.
  scores(tilde);

  // Each parameter's loadings.
  const Mat eta2 = sym_product(s.eta, false), eta_t = transpose(s.eta);
  for (int j = 0; j < p; ++j) {
    const double ej = s.epsilon_inv[j];
    Mat q(K, K), y(n, 1);
    for (int b = 0; b < K; ++b)
      for (int a = 0; a < K; ++a) {
        volatile double scaled = ej * eta2(a, b);
        q(a, b) = (a == b ? plam(j, a) : 0.0) + scaled;
      }
    for (int i = 0; i < n; ++i) y(i, 0) = tilde(i, j);
    Mat rhs = mult(eta_t, y);
    for (double& e : rhs.v) e = ej * e;
    const Mat upper = chol_upper(q), lower = transpose(upper);
    Mat z(K, 1);
    for (double& e : z.v) e = norm_rand();
    const Mat mean = tri_solve(upper, tri_solve(lower, rhs, false), true);
    const Mat noise = tri_solve(upper, z, true);
    for (int f = 0; f < K; ++f) s.lambda(j, f) = noise(f, 0) + mean(f, 0);
  }

  // The local shrinkage of each loading.
  for (int h = 0; h < K; ++h)
    for (int j = 0; j < p; ++j) {
      volatile double l2 = s.lambda(j, h) * s.lambda(j, h);
      volatile double scaled = tauh[h] * l2;
      s.psi(j, h) = gamma_draw((prior.df + 1) / 2, (prior.df + scaled) / 2);
    }

  // The factors' global shrinkage.
  {
    Mat m(p, K);
    for (size_t e = 0; e < m.v.size(); ++e) {
      volatile double l2 = s.lambda.v[e] * s.lambda.v[e];
      volatile double w = s.psi.v[e] * l2;
      m.v[e] = w;
    }
    const std::vector<double> sums = col_sums(m);
    auto tail = [&](int from) {
      long double total = 0.0;
      for (int h = from; h < K; ++h) {
        volatile double term = tauh[h] * sums[h];
        total += term;
      }
      return static_cast<double>(total);
    };
    {
      volatile double count = 0.5 * p;
      volatile double shape_part = count * K;
      volatile double half = 0.5 * (1 / s.delta[0]);
      volatile double rate_part = half * tail(0);
      s.delta[0] = gamma_draw(prior.ad1 + shape_part, prior.bd1 + rate_part);
    }
    tauh = cum_prod(s.delta);
    if (prior.has_ad2) {
      for (int h = 1; h < K; ++h) {
        volatile double count = 0.5 * p;
        volatile double shape_part = count * (K - h);
        volatile double half = 0.5 * (1 / s.delta[h]);
        volatile double rate_part = half * tail(h);
        s.delta[h] = gamma_draw(prior.ad2 + shape_part, prior.bd2 + rate_part);
        tauh = cum_prod(s.delta);
      }
    } else {
      for (int h = 1; h < K; ++h) s.delta[h] = 1.0;
    }
  }

  // The error precisions.
  {
    Mat sq = residual(tilde);
    for (double& e : sq.v) {
      volatile double e2 = e * e;
      e = e2;
    }
    const std::vector<double> ss = col_sums(sq);
    for (int j = 0; j < p; ++j) {
      volatile double half = .5 * ss[j];
      s.epsilon_inv[j] = gamma_draw(prior.shape_err[j], prior.rate_err[j] + half);
    }
  }

  // The implied group covariance.
  s.var = mult(s.lambda, transpose(s.lambda));
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i)
      s.var(i, j) = s.var(i, j) + (i == j ? 1.0 / s.epsilon_inv[i] : 0.0);
}

// crossprod(x, y) for a single-column y.
Mat cross_column(const Mat& x, const Mat& y) {
  if (y.c != 1 || y.r != x.r || x.r == 0 || x.c == 0) throw Decline();
  require_finite(x);
  require_finite(y);
  Mat z(x.c, 1);
  int nr = x.r, nc = x.c, ione = 1;
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemv)("T", &nr, &nc, &one, x.v.data(), &nr, y.v.data(), &ione, &zero, z.v.data(),
                  &ione FCONE);
  return z;
}

// The columns `idx` of x.
Mat columns(const Mat& x, const std::vector<int>& idx) {
  Mat z(x.r, static_cast<int>(idx.size()));
  for (size_t f = 0; f < idx.size(); ++f)
    for (int i = 0; i < x.r; ++i) z(i, static_cast<int>(f)) = x(i, idx[f]);
  return z;
}

// A regression row's draw: coefficients of the columns `fr` of `design` for
// the response net of the fixed columns' part, at residual precision `tau`
// and prior precisions `prec` (by design column). `coef` holds the row's
// current coefficients and takes the draw.
void regression_row(const Mat& design, Mat response, double tau, const std::vector<int>& fr,
                    const std::vector<int>& fx, const std::vector<double>& prec,
                    std::vector<double>& coef, bool scale_response) {
  const int n = design.r, kf = static_cast<int>(fr.size());
  if (!fx.empty()) {
    Mat fixed(static_cast<int>(fx.size()), 1);
    for (size_t f = 0; f < fx.size(); ++f) fixed(static_cast<int>(f), 0) = coef[fx[f]];
    const Mat part = mult(columns(design, fx), fixed);
    for (int i = 0; i < n; ++i) response(i, 0) = response(i, 0) - part(i, 0);
  }
  const Mat x = columns(design, fr);
  Mat m = sym_product(x, false);
  for (int b = 0; b < kf; ++b)
    for (int a = 0; a < kf; ++a) {
      volatile double scaled = tau * m(a, b);
      m(a, b) = scaled + (a == b ? prec[fr[a]] : 0.0);
    }
  const Mat sig = inverse(m);
  Mat rhs = cross_column(x, response);
  if (scale_response)
    for (double& e : rhs.v) e = tau * e;
  const Mat mean = mult(sig, rhs);
  const Mat row = normal_rows(1, mean.v.data(), normal_factor(sig));
  for (int f = 0; f < kf; ++f) coef[fr[f]] = row(0, f);
}

// The structural model's factor scores drawn given the subjects' states
// less their means (subjects x parameters), the covariates and the error
// precisions as a diagonal matrix (.sem_draw_eta).
void sem_scores(emc_factor::SemState& s, const Mat& tilde, const Mat& cov, const Mat& eps) {
  const int n = tilde.r, k = s.lambda.c;
  Mat ib(k, k);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) ib(i, j) = (i == j ? 1.0 : 0.0) - s.B(i, j);
  const Mat b0 = inverse(ib);
  const Mat psi0 = inverse(mult(mult(b0, inverse(s.delta_inv)), transpose(b0)));
  const Mat lt = transpose(s.lambda);
  Mat prec = mult(mult(lt, eps), s.lambda);
  for (size_t e = 0; e < prec.v.size(); ++e) prec.v[e] = psi0.v[e] + prec.v[e];
  const Mat eta_sig = inverse(prec);
  Mat centred = tilde;
  const Mat effects = mult(cov, transpose(s.K));
  for (size_t e = 0; e < centred.v.size(); ++e) centred.v[e] = centred.v[e] - effects.v[e];
  const Mat data = mult(mult(lt, eps), transpose(centred));
  const Mat structural = mult(mult(mult(psi0, b0), s.G), transpose(cov));
  Mat total(k, n);
  for (size_t e = 0; e < total.v.size(); ++e) total.v[e] = data.v[e] + structural.v[e];
  const Mat mean = mult(eta_sig, total);
  const Mat noise = normal_rows(n, NULL, normal_factor(eta_sig));
  if (s.eta.r != n || s.eta.c != k) s.eta = Mat(n, k);
  for (int f = 0; f < k; ++f)
    for (int i = 0; i < n; ++i) s.eta(i, f) = mean(f, i) + noise(i, f);
}

void draw_sem(emc_factor::SemState& s, const Mat& y, const Mat& cov,
              const emc_factor::SemPrior& prior) {
  const int n = y.r, p = y.c, k = s.lambda.c, c = cov.c;
  Mat eps = diagonal(s.epsilon_inv);
  auto centred = [&](const std::vector<double>& mu) {
    Mat z(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) z(i, j) = y(i, j) - mu[j];
    return z;
  };
  auto minus = [](Mat a, const Mat& b) {
    for (size_t e = 0; e < a.v.size(); ++e) a.v[e] = a.v[e] - b.v[e];
    return a;
  };
  auto scores = [&](const Mat& tilde) { sem_scores(s, tilde, cov, eps); };
  auto residual = [&](const Mat& base) {
    return minus(minus(base, mult(cov, transpose(s.K))), mult(s.eta, transpose(s.lambda)));
  };

  s.eta = Mat(n, k);
  scores(centred(s.mu));

  // The means.
  {
    Mat prec(p, p), prior_prec = diagonal(prior.mu_prec);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < p; ++i) {
        volatile double scaled = n * eps(i, j);
        prec(i, j) = scaled + prior_prec(i, j);
      }
    const Mat mu_sig = inverse(prec);
    Mat sums(p, 1), mean0(p, 1);
    sums.v = col_sums(residual(y));
    mean0.v = prior.mu_mean;
    const Mat a = mult(eps, sums), b = mult(prior_prec, mean0);
    Mat rhs(p, 1);
    for (int i = 0; i < p; ++i) rhs(i, 0) = a(i, 0) + b(i, 0);
    const Mat mu_mu = mult(mu_sig, rhs);
    s.mu = normal_rows(1, mu_mu.v.data(), normal_factor(mu_sig)).v;
  }
  const Mat tilde = centred(s.mu);
  scores(tilde);

  // The error precisions.
  {
    Mat sq = residual(tilde);
    for (double& e : sq.v) {
      volatile double e2 = e * e;
      e = e2;
    }
    const std::vector<double> ss = col_sums(sq);
    for (int j = 0; j < p; ++j) {
      volatile double half = 0.5 * ss[j];
      s.epsilon_inv[j] = gamma_draw(prior.shape_err[j], prior.rate_err[j] + half);
    }
  }

  // Each parameter's covariate effects and loadings.
  {
    Mat design(n, c + k);
    std::copy(cov.v.begin(), cov.v.end(), design.v.begin());
    std::copy(s.eta.v.begin(), s.eta.v.end(), design.v.begin() + cov.v.size());
    for (int j = 0; j < p; ++j) {
      std::vector<int> fr, fx;
      std::vector<double> coef(c + k);
      for (int f = 0; f < c; ++f) {
        (prior.free_K[j + static_cast<size_t>(f) * p] ? fr : fx).push_back(f);
        coef[f] = s.K(j, f);
      }
      for (int f = 0; f < k; ++f) {
        (prior.free_lambda[j + static_cast<size_t>(f) * p] ? fr : fx).push_back(c + f);
        coef[c + f] = s.lambda(j, f);
      }
      if (fr.empty()) continue;
      Mat response(n, 1);
      for (int i = 0; i < n; ++i) response(i, 0) = tilde(i, j);
      regression_row(design, response, s.epsilon_inv[j], fr, fx, prior.load_prec, coef, true);
      for (int f = 0; f < c; ++f) s.K(j, f) = coef[f];
      for (int f = 0; f < k; ++f) s.lambda(j, f) = coef[c + f];
    }
  }

  // Each factor's structural coefficients given the others' residuals.
  {
    Mat design(n, c + k);
    std::copy(cov.v.begin(), cov.v.end(), design.v.begin());
    std::copy(s.eta.v.begin(), s.eta.v.end(), design.v.begin() + cov.v.size());
    Mat resid = minus(minus(s.eta, mult(s.eta, transpose(s.B))), mult(cov, transpose(s.G)));
    for (int f = 0; f < k; ++f) {
      std::vector<int> fr, fx;
      std::vector<double> coef(c + k);
      for (int g = 0; g < c; ++g) {
        (prior.free_G[f + static_cast<size_t>(g) * k] ? fr : fx).push_back(g);
        coef[g] = s.G(f, g);
      }
      for (int g = 0; g < k; ++g) {
        (prior.free_B[f + static_cast<size_t>(g) * k] ? fr : fx).push_back(c + g);
        coef[c + g] = s.B(f, g);
      }
      if (fr.empty()) continue;
      const double tau = s.delta_inv(f, f);
      Mat others(n, 1);
      if (k > 1) {
        std::vector<int> rest;
        for (int g = 0; g < k; ++g)
          if (g != f) rest.push_back(g);
        Mat weights(k - 1, 1);
        for (int g = 0; g < k - 1; ++g) weights(g, 0) = s.delta_inv(rest[g], f);
        others = mult(columns(resid, rest), weights);
      }
      // The response is scaled by the factor's own precision after the
      // fixed coefficients' part is removed.
      Mat response(n, 1);
      for (int i = 0; i < n; ++i) response(i, 0) = s.eta(i, f);
      if (!fx.empty()) {
        Mat fixed(static_cast<int>(fx.size()), 1);
        for (size_t g = 0; g < fx.size(); ++g) fixed(static_cast<int>(g), 0) = coef[fx[g]];
        const Mat part = mult(columns(design, fx), fixed);
        for (int i = 0; i < n; ++i) response(i, 0) = response(i, 0) - part(i, 0);
      }
      for (int i = 0; i < n; ++i) {
        volatile double scaled = tau * response(i, 0);
        response(i, 0) = scaled + others(i, 0);
      }
      regression_row(design, response, tau, fr, std::vector<int>(), prior.struct_prec, coef,
                     false);
      for (int g = 0; g < c; ++g) s.G(f, g) = coef[g];
      for (int g = 0; g < k; ++g) s.B(f, g) = coef[c + g];
      Mat row(c + k, 1);
      row.v = coef;
      const Mat fitted = mult(design, row);
      for (int i = 0; i < n; ++i) resid(i, f) = s.eta(i, f) - fitted(i, 0);
    }
  }
}

// The reciprocal 1-norm condition number of a square matrix, zero where it
// is singular.
double reciprocal_condition(const Mat& a) {
  int n = a.r, info = 0;
  Mat lu = a;
  std::vector<int> ipiv(n);
  F77_CALL(dgetrf)(&n, &n, lu.v.data(), &n, ipiv.data(), &info);
  if (info != 0) return 0.0;
  double anorm = F77_CALL(dlange)("1", &n, &n, a.v.data(), &n, NULL FCONE);
  std::vector<double> work(4 * static_cast<size_t>(n));
  double rcond = 0.0;
  F77_CALL(dgecon)("1", &n, lu.v.data(), &n, &anorm, &rcond, work.data(), ipiv.data(),
                   &info FCONE);
  return info == 0 ? rcond : 0.0;
}

// riwish(v, scale): the inverse of a Wishart draw on the scale's inverse,
// from the Bartlett factor with its chi-square diagonal drawn before its
// normal upper triangle. A scale near the conditioning riwish() refuses is
// left to the R form.
Mat inverse_wishart(double v, const Mat& scale) {
  const int p = scale.r;
  if (scale.c != p || p < 2) throw Decline();
  require_finite(scale);
  Mat s(p, p);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) s(i, j) = (scale(i, j) + scale(j, i)) / 2;
  chol_upper(s);
  if (!(reciprocal_condition(s) > 1e-8)) throw Decline();
  const Mat s_inv = inverse(s);
  if (!R_FINITE(v) || v <= p - 1) throw Decline();
  const Mat cc = chol_upper(s_inv);
  Mat z(p, p);
  for (int i = 0; i < p; ++i) {
    const double chi = Rf_rchisq(v - static_cast<double>(i));
    if (!R_FINITE(chi)) throw Decline();
    z(i, i) = std::sqrt(chi);
  }
  for (int j = 1; j < p; ++j)
    for (int i = 0; i < j; ++i) z(i, j) = norm_rand();
  return inverse(sym_product(mult(z, cc), false));
}

Mat symmetrised(const Mat& x) {
  Mat z(x.r, x.c);
  for (int j = 0; j < x.c; ++j)
    for (int i = 0; i < x.r; ++i) z(i, j) = (x(i, j) + x(j, i)) / 2;
  return z;
}

void finish(emc_factor::SemState& s, const Mat& cov, const emc_factor::SemBlocks& b,
            emc_factor::SemMoments& m) {
  const int n = s.eta.r, k = s.eta.c, p = s.lambda.r, c = cov.c;
  const Mat own = mult(s.eta, transpose(s.B));
  const Mat driven = mult(cov, transpose(s.G));
  Mat resid(n, k);
  for (int f = 0; f < k; ++f)
    for (int i = 0; i < n; ++i) resid(i, f) = s.eta(i, f) - own(i, f) - driven(i, f);
  int n_blocks = 0;
  for (int f = 0; f < k; ++f) n_blocks = std::max(n_blocks, b.block[f] + 1);
  for (int g = 0; g < n_blocks; ++g) {
    std::vector<int> idx;
    for (int f = 0; f < k; ++f)
      if (b.block[f] == g) idx.push_back(f);
    const int d = static_cast<int>(idx.size());
    if (d == 0) throw Decline();
    const Mat cross = sym_product(columns(resid, idx), false);
    const double df = b.a_d + n;
    if (!(df > d - 1)) throw Decline();
    if (d > 1) {
      Mat scale(d, d);
      for (int j = 0; j < d; ++j)
        for (int i = 0; i < d; ++i) scale(i, j) = (i == j ? b.b_d : 0.0) + cross(i, j);
      const Mat precision = inverse(inverse_wishart(df, scale));
      for (int j = 0; j < d; ++j)
        for (int i = 0; i < d; ++i) s.delta_inv(idx[i], idx[j]) = precision(i, j);
    } else {
      volatile double half = 0.5 * cross(0, 0);
      s.delta_inv(idx[0], idx[0]) = gamma_draw(b.a_d + n / 2.0, b.b_d + half);
    }
  }
  // Loadings on the disturbances and the covariates' total effects.
  Mat structure(k, k);
  for (int j = 0; j < k; ++j)
    for (int i = 0; i < k; ++i) structure(i, j) = (i == j ? 1.0 : 0.0) - s.B(i, j);
  const Mat a = mult(s.lambda, inverse(structure));
  Mat total = mult(a, s.G);
  for (int j = 0; j < c; ++j)
    for (int i = 0; i < p; ++i) total(i, j) = s.K(i, j) + total(i, j);
  Mat v = mult(mult(a, inverse(s.delta_inv)), transpose(a));
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < p; ++i) v(i, j) = v(i, j) + (i == j ? 1 / s.epsilon_inv[i] : 0.0);
  m.subj_var = symmetrised(v);
  m.subj_mu = Mat(p, n);
  if (c > 0) {
    const Mat shift = mult(total, transpose(cov));
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < p; ++j) m.subj_mu(j, i) = s.mu[j] + shift(j, i);
    Mat x_mean(c, 1);
    x_mean.v = b.x_mean;
    const Mat mean_shift = mult(total, x_mean);
    m.pop_mean.resize(p);
    for (int j = 0; j < p; ++j) m.pop_mean[j] = s.mu[j] + mean_shift(j, 0);
    const Mat spread = mult(mult(total, b.x_var), transpose(total));
    Mat pv(p, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < p; ++i) pv(i, j) = m.subj_var(i, j) + spread(i, j);
    m.pop_var = symmetrised(pv);
  } else {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < p; ++j) m.subj_mu(j, i) = s.mu[j];
    m.pop_mean = s.mu;
    m.pop_var = m.subj_var;
  }
}

Mat as_mat(const Rcpp::NumericMatrix& x) {
  Mat z(x.nrow(), x.ncol());
  std::copy(x.begin(), x.end(), z.v.begin());
  return z;
}

Rcpp::NumericMatrix as_r(const Mat& x) {
  Rcpp::NumericMatrix z(x.r, x.c);
  std::copy(x.v.begin(), x.v.end(), z.begin());
  return z;
}

std::vector<double> as_std(const Rcpp::NumericVector& x) {
  return std::vector<double>(x.begin(), x.end());
}

emc_factor::State factor_state(const Rcpp::List& x) {
  emc_factor::State s;
  s.eta = as_mat(Rcpp::as<Rcpp::NumericMatrix>(x["eta"]));
  s.lambda = as_mat(Rcpp::as<Rcpp::NumericMatrix>(x["lambda"]));
  s.mu = as_std(x["mu"]);
  s.sig_err_inv = as_std(x["sig_err_inv"]);
  s.psi_inv = as_std(x["psi_inv"]);
  return s;
}

emc_factor::Prior factor_prior(const Rcpp::List& x) {
  emc_factor::Prior prior;
  prior.mu_invar = as_mat(Rcpp::as<Rcpp::NumericMatrix>(x["mu_invar"]));
  prior.mu_mean = as_std(x["mu_mean"]);
  prior.shape_err = as_std(x["shape_err"]);
  prior.rate_err = as_std(x["rate_err"]);
  prior.lambda_invar = as_std(x["lambda_invar"]);
  prior.shape_psi = as_std(x["shape_psi"]);
  prior.rate_psi = as_std(x["rate_psi"]);
  return prior;
}

// Which loadings are free; FALSE with a missing entry.
bool free_entries(const Rcpp::LogicalMatrix& x, std::vector<int>& out) {
  out.resize(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (x[i] == NA_LOGICAL) return false;
    out[i] = x[i] == TRUE;
  }
  return true;
}

Rcpp::List factor_result(const emc_factor::State& s) {
  return Rcpp::List::create(
      Rcpp::_["mu"] = Rcpp::NumericVector(s.mu.begin(), s.mu.end()),
      Rcpp::_["var"] = as_r(s.var), Rcpp::_["lambda"] = as_r(s.lambda),
      Rcpp::_["lambda_scaled"] = as_r(s.lambda_scaled), Rcpp::_["eta"] = as_r(s.eta),
      Rcpp::_["sig_err_inv"] = Rcpp::NumericVector(s.sig_err_inv.begin(), s.sig_err_inv.end()),
      Rcpp::_["psi_inv"] = Rcpp::NumericVector(s.psi_inv.begin(), s.psi_inv.end()));
}

emc_factor::InfiniteState infinite_state(const Rcpp::List& x) {
  emc_factor::InfiniteState s;
  s.psi = as_mat(Rcpp::as<Rcpp::NumericMatrix>(x["psi"]));
  s.lambda = as_mat(Rcpp::as<Rcpp::NumericMatrix>(x["lambda"]));
  s.mu = as_std(x["mu"]);
  s.epsilon_inv = as_std(x["epsilon_inv"]);
  s.delta = as_std(x["delta"]);
  return s;
}

emc_factor::InfinitePrior infinite_prior(const Rcpp::List& x) {
  emc_factor::InfinitePrior prior;
  prior.mu_invar = as_std(x["mu_invar"]);
  prior.mu_mean = as_std(x["mu_mean"]);
  prior.shape_err = as_std(x["shape_err"]);
  prior.rate_err = as_std(x["rate_err"]);
  prior.df = Rcpp::as<double>(x["df"]);
  prior.ad1 = Rcpp::as<double>(x["ad1"]);
  prior.bd1 = Rcpp::as<double>(x["bd1"]);
  SEXP ad2 = x["ad2"];
  prior.has_ad2 = ad2 != R_NilValue;
  prior.ad2 = prior.has_ad2 ? Rf_asReal(ad2) : 0.0;
  prior.bd2 = prior.has_ad2 ? Rcpp::as<double>(x["bd2"]) : 0.0;
  return prior;
}

Rcpp::List infinite_result(const emc_factor::InfiniteState& s) {
  return Rcpp::List::create(
      Rcpp::_["mu"] = Rcpp::NumericVector(s.mu.begin(), s.mu.end()),
      Rcpp::_["var"] = as_r(s.var), Rcpp::_["lambda"] = as_r(s.lambda),
      Rcpp::_["eta"] = as_r(s.eta),
      Rcpp::_["epsilon_inv"] = Rcpp::NumericVector(s.epsilon_inv.begin(), s.epsilon_inv.end()),
      Rcpp::_["psi"] = as_r(s.psi),
      Rcpp::_["delta"] = Rcpp::NumericVector(s.delta.begin(), s.delta.end()));
}

}  // namespace

namespace emc_factor {

bool gibbs(State& state, const Mat& alpha, const Prior& prior,
           const std::vector<int>& free_entry) {
  const int n = alpha.r, p = alpha.c, k = state.eta.c;
  if (n < 1 || p < 2 || k < 1 || state.eta.r != n || state.lambda.r != p ||
      state.lambda.c != k || static_cast<int>(state.mu.size()) != p ||
      static_cast<int>(state.sig_err_inv.size()) != p ||
      static_cast<int>(state.psi_inv.size()) != k || prior.mu_invar.r != p ||
      prior.mu_invar.c != p || static_cast<int>(prior.mu_mean.size()) != p ||
      static_cast<int>(prior.shape_err.size()) != p ||
      static_cast<int>(prior.rate_err.size()) != p ||
      static_cast<int>(prior.lambda_invar.size()) != p ||
      static_cast<int>(prior.shape_psi.size()) != k ||
      static_cast<int>(prior.rate_psi.size()) != k ||
      free_entry.size() != static_cast<size_t>(p) * k)
    return false;
  try {
    draw(state, alpha, prior, free_entry);
  } catch (Decline&) {
    return false;
  }
  return true;
}

bool gibbs_infinite(InfiniteState& state, const Mat& alpha, const InfinitePrior& prior) {
  const int n = alpha.r, p = alpha.c, K = state.lambda.c;
  auto len = [](const std::vector<double>& x) { return static_cast<int>(x.size()); };
  if (n < 1 || p < 2 || K < 1 || state.lambda.r != p || state.psi.r != p || state.psi.c != K ||
      len(state.mu) != p || len(state.epsilon_inv) != p || len(state.delta) != K ||
      len(prior.mu_invar) != p || len(prior.mu_mean) != p || len(prior.shape_err) != p ||
      len(prior.rate_err) != p)
    return false;
  try {
    draw_infinite(state, alpha, prior);
  } catch (Decline&) {
    return false;
  }
  return true;
}

bool gibbs_sem(SemState& state, const Mat& y, const Mat& cov, const SemPrior& prior) {
  const int n = y.r, p = y.c, k = state.lambda.c, c = cov.c;
  auto len = [](const std::vector<double>& x) { return static_cast<int>(x.size()); };
  auto cells = [](const std::vector<int>& x) { return x.size(); };
  if (n < 1 || p < 2 || k < 1 || cov.r != n || state.lambda.r != p || state.B.r != k ||
      state.B.c != k || state.K.r != p || state.K.c != c || state.G.r != k || state.G.c != c ||
      state.delta_inv.r != k || state.delta_inv.c != k || len(state.mu) != p ||
      len(state.epsilon_inv) != p || len(prior.mu_prec) != p || len(prior.mu_mean) != p ||
      len(prior.shape_err) != p || len(prior.rate_err) != p || len(prior.load_prec) != c + k ||
      len(prior.struct_prec) != c + k ||
      cells(prior.free_lambda) != static_cast<size_t>(p) * k ||
      cells(prior.free_K) != static_cast<size_t>(p) * c ||
      cells(prior.free_B) != static_cast<size_t>(k) * k ||
      cells(prior.free_G) != static_cast<size_t>(k) * c)
    return false;
  try {
    draw_sem(state, y, cov, prior);
  } catch (Decline&) {
    return false;
  }
  return true;
}

bool finish_sem(SemState& state, const Mat& cov, const SemBlocks& blocks,
                SemMoments& moments) {
  const int n = state.eta.r, k = state.eta.c, c = cov.c;
  if (n < 1 || k < 1 || cov.r != n || static_cast<int>(blocks.block.size()) != k ||
      static_cast<int>(blocks.x_mean.size()) != c || blocks.x_var.r != c ||
      blocks.x_var.c != c || !R_FINITE(blocks.a_d) || !R_FINITE(blocks.b_d))
    return false;
  for (int g : blocks.block)
    if (g < 0 || g >= k) return false;
  try {
    finish(state, cov, blocks, moments);
  } catch (Decline&) {
    return false;
  }
  return true;
}

}  // namespace emc_factor

// The structural equation model's state, prior and blocks from the lists
// emc_gibbs_sem_cpp() takes; FALSE where which coefficients are free is
// missing.
emc_factor::SemState sem_state(const Rcpp::List& state) {
  emc_factor::SemState s;
  s.lambda = as_mat(Rcpp::as<Rcpp::NumericMatrix>(state["lambda"]));
  s.B = as_mat(Rcpp::as<Rcpp::NumericMatrix>(state["B"]));
  s.K = as_mat(Rcpp::as<Rcpp::NumericMatrix>(state["K"]));
  s.G = as_mat(Rcpp::as<Rcpp::NumericMatrix>(state["G"]));
  s.delta_inv = as_mat(Rcpp::as<Rcpp::NumericMatrix>(state["delta_inv"]));
  s.mu = as_std(state["mu"]);
  s.epsilon_inv = as_std(state["epsilon_inv"]);
  return s;
}

bool sem_prior(const Rcpp::List& prior, const Rcpp::List& free, emc_factor::SemPrior& pr) {
  pr.mu_prec = as_std(prior["mu_prec"]);
  pr.mu_mean = as_std(prior["mu_mean"]);
  pr.shape_err = as_std(prior["shape_err"]);
  pr.rate_err = as_std(prior["rate_err"]);
  pr.load_prec = as_std(prior["load_prec"]);
  pr.struct_prec = as_std(prior["struct_prec"]);
  return free_entries(Rcpp::as<Rcpp::LogicalMatrix>(free["lambda"]), pr.free_lambda) &&
         free_entries(Rcpp::as<Rcpp::LogicalMatrix>(free["K"]), pr.free_K) &&
         free_entries(Rcpp::as<Rcpp::LogicalMatrix>(free["B"]), pr.free_B) &&
         free_entries(Rcpp::as<Rcpp::LogicalMatrix>(free["G"]), pr.free_G);
}

emc_factor::SemBlocks sem_blocks(const Rcpp::List& blocks) {
  emc_factor::SemBlocks bl;
  for (int g : Rcpp::as<Rcpp::IntegerVector>(blocks["block"])) bl.block.push_back(g - 1);
  bl.a_d = Rcpp::as<double>(blocks["a_d"]);
  bl.b_d = Rcpp::as<double>(blocks["b_d"]);
  bl.x_mean = as_std(blocks["x_mean"]);
  bl.x_var = as_mat(Rcpp::as<Rcpp::NumericMatrix>(blocks["x_var"]));
  return bl;
}

Rcpp::List sem_result(const emc_factor::SemState& s, const emc_factor::SemMoments& mom) {
  return Rcpp::List::create(
      Rcpp::_["tmu"] = Rcpp::NumericVector(s.mu.begin(), s.mu.end()),
      Rcpp::_["lambda"] = as_r(s.lambda), Rcpp::_["eta"] = as_r(s.eta),
      Rcpp::_["B"] = as_r(s.B), Rcpp::_["K"] = as_r(s.K), Rcpp::_["G"] = as_r(s.G),
      Rcpp::_["epsilon_inv"] = Rcpp::NumericVector(s.epsilon_inv.begin(), s.epsilon_inv.end()),
      Rcpp::_["delta_inv"] = as_r(s.delta_inv), Rcpp::_["subj_mu"] = as_r(mom.subj_mu),
      Rcpp::_["subj_var"] = as_r(mom.subj_var),
      Rcpp::_["pop_mean"] = Rcpp::NumericVector(mom.pop_mean.begin(), mom.pop_mean.end()),
      Rcpp::_["tvar"] = as_r(mom.pop_var));
}

// The structural equation model's group draw: `y` subjects x parameters, the
// covariates (no columns without any), the state as the store holds it,
// which coefficients are free, the prior's precisions at full length, and
// the factors' blocks (numbered from one) with the disturbance prior and the
// covariates' means and covariance. NULL when the draw is left to the R
// form.
// [[Rcpp::export]]
SEXP emc_gibbs_sem_cpp(const Rcpp::NumericMatrix& y, const Rcpp::NumericMatrix& covariates,
                       const Rcpp::List& state, const Rcpp::List& free,
                       const Rcpp::List& prior, const Rcpp::List& blocks) {
  emc_factor::SemState s = sem_state(state);
  emc_factor::SemPrior pr;
  if (!sem_prior(prior, free, pr)) return R_NilValue;
  const emc_factor::SemBlocks bl = sem_blocks(blocks);
  emc_factor::SemMoments mom;
  const Mat cov = as_mat(covariates);
  if (!emc_factor::gibbs_sem(s, as_mat(y), cov, pr) ||
      !emc_factor::finish_sem(s, cov, bl, mom))
    return R_NilValue;
  return sem_result(s, mom);
}

// .emc_ensemble_refresh_scores()'s draw of a latent hierarchy's factor
// scores given the subjects' states (parameters x subjects) and the group
// state as the iteration holds it: `kind` 0 factor, 1 shrinkage prior, 2
// structural (with its covariates). The scores (subjects x factors), or NULL
// when the draw is left to the R form.
// [[Rcpp::export]]
SEXP emc_refresh_scores_cpp(const Rcpp::NumericMatrix& alpha, int kind, const Rcpp::List& pars,
                            SEXP covariates) {
  const int p = alpha.nrow(), n = alpha.ncol();
  auto real = [&](const char* name) -> SEXP {
    SEXP names = pars.names();
    for (R_xlen_t i = 0; i < pars.size(); ++i)
      if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
        SEXP v = VECTOR_ELT(pars, i);
        return TYPEOF(v) == REALSXP ? v : R_NilValue;
      }
    return R_NilValue;
  };
  auto shaped = [&](SEXP v, int nr, int nc) {
    Mat m(nr, nc);
    if (v == R_NilValue || XLENGTH(v) != static_cast<R_xlen_t>(nr) * nc) throw Decline();
    std::copy(REAL(v), REAL(v) + XLENGTH(v), m.v.begin());
    return m;
  };
  auto vec = [&](SEXP v, int len) {
    if (v == R_NilValue || XLENGTH(v) != len) throw Decline();
    return std::vector<double>(REAL(v), REAL(v) + len);
  };
  if (p < 1 || n < 1 || kind < 0 || kind > 2 || TYPEOF(pars.names()) != STRSXP)
    return R_NilValue;
  try {
    SEXP load = real(kind == 0 ? "lambda_untransf" : "lambda");
    if (load == R_NilValue || XLENGTH(load) % p != 0) throw Decline();
    const int k = static_cast<int>(XLENGTH(load) / p);
    if (k < 1) throw Decline();
    const Mat loadings = shaped(load, p, k);
    const std::vector<double> mu = vec(real("tmu"), p);
    Mat residual(n, p);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i < n; ++i) residual(i, j) = alpha(j, i) - mu[j];
    if (kind == 2) {
      if (TYPEOF(covariates) != REALSXP || !Rf_isMatrix(covariates) ||
          Rf_nrows(covariates) != n)
        throw Decline();
      const int c = Rf_ncols(covariates);
      emc_factor::SemState s;
      s.lambda = loadings;
      s.K = shaped(real("K"), p, c);
      s.G = shaped(real("G"), k, c);
      s.B = shaped(real("B"), k, k);
      s.delta_inv = shaped(real("delta_inv"), k, k);
      sem_scores(s, residual, shaped(covariates, n, c), diagonal(vec(real("epsilon_inv"), p)));
      return as_r(s.eta);
    }
    const std::vector<double> error = vec(real(kind == 0 ? "sig_err_inv" : "epsilon_inv"), p);
    const std::vector<double> score =
        kind == 0 ? vec(real("psi_inv"), k) : std::vector<double>(k, 1.0);
    Mat weighted(p, k);
    for (int f = 0; f < k; ++f)
      for (int j = 0; j < p; ++j) weighted(j, f) = loadings(j, f) * error[j];
    // The two-matrix cross product: one general product, or the level-2
    // routine for a single factor.
    Mat prec(k, k);
    if (k == 1) {
      prec = cross_column(loadings, weighted);
    } else {
      require_finite(loadings);
      require_finite(weighted);
      int pp = p, kk = k;
      double one = 1.0, zero = 0.0;
      F77_CALL(dgemm)("T", "N", &kk, &kk, &pp, &one, loadings.v.data(), &pp, weighted.v.data(),
                      &pp, &zero, prec.v.data(), &kk FCONE FCONE);
    }
    for (int f = 0; f < k; ++f)
      for (int g = 0; g < k; ++g) prec(g, f) = (g == f ? score[f] : 0.0) + prec(g, f);
    const Mat covariance = inverse(prec);
    const Mat mean = mult(mult(residual, weighted), covariance);
    const Mat noise = normal_rows(n, NULL, normal_factor(covariance));
    Mat out(n, k);
    for (size_t e = 0; e < out.v.size(); ++e) out.v[e] = mean.v[e] + noise.v[e];
    return as_r(out);
  } catch (Decline&) {
    return R_NilValue;
  }
}

// mvtnorm::rmvnorm(n, mean, sigma) by its default method, on the same random
// stream: `mean` NULL for the zero vector. NULL where that function would
// raise a condition; the stream is then spent.
// [[Rcpp::export]]
SEXP emc_rmvnorm_cpp(int n, SEXP mean, const Rcpp::NumericMatrix& sigma) {
  const int k = sigma.nrow();
  if (n < 1 || k < 1 || sigma.ncol() != k ||
      (mean != R_NilValue && (TYPEOF(mean) != REALSXP || XLENGTH(mean) != k)))
    return R_NilValue;
  try {
    const Mat out = normal_rows(n, mean == R_NilValue ? NULL : REAL(mean),
                                normal_factor(as_mat(sigma)));
    return as_r(out);
  } catch (Decline&) {
    return R_NilValue;
  }
}

// solve(a) of a square matrix without dimnames; NULL where solve() would
// raise a condition.
// [[Rcpp::export(rng = false)]]
SEXP emc_solve_cpp(const Rcpp::NumericMatrix& a) {
  if (a.nrow() != a.ncol() || a.nrow() < 1) return R_NilValue;
  try {
    return as_r(inverse(as_mat(a)));
  } catch (Decline&) {
    return R_NilValue;
  }
}

// gibbs_step_factor()'s draw: `alpha` subjects x parameters, the state as
// the store holds it (eta, lambda, mu, sig_err_inv, psi_inv), the free
// loadings and the prior's vectors at full length. NULL when the draw is
// left to the R form, which then starts again from the stream's state before
// the call.
// [[Rcpp::export]]
SEXP emc_gibbs_factor_cpp(const Rcpp::NumericMatrix& alpha, const Rcpp::List& state,
                          const Rcpp::LogicalMatrix& free_entry, const Rcpp::List& prior) {
  emc_factor::State s = factor_state(state);
  std::vector<int> fr;
  if (!free_entries(free_entry, fr)) return R_NilValue;
  if (!emc_factor::gibbs(s, as_mat(alpha), factor_prior(prior), fr)) return R_NilValue;
  return factor_result(s);
}

// gibbs_step_infnt_factor()'s draw, likewise: the state is psi, lambda, mu,
// epsilon_inv and delta; the prior's `ad2` is NULL for a prior without it.
// [[Rcpp::export]]
SEXP emc_gibbs_infnt_factor_cpp(const Rcpp::NumericMatrix& alpha, const Rcpp::List& state,
                                const Rcpp::List& prior) {
  emc_factor::InfiniteState s = infinite_state(state);
  if (!emc_factor::gibbs_infinite(s, as_mat(alpha), infinite_prior(prior))) return R_NilValue;
  return infinite_result(s);
}

// The ensemble move's Gibbs rounds for a factor hierarchy (`free_entry` a
// matrix), the shrinkage-prior one (`free_entry` NULL) or a structural
// equation model (`sem` its covariates, free coefficients and blocks as
// emc_gibbs_sem_cpp() takes them), as .emc_ensemble_rounds_r() runs them:
// the group draw given the selected pool members, every member's density
// under the subject prior given the scores, and each subject's selection by
// the Gumbel-max trick, one uniform per pool row in stack order. Returns the
// final and the last conditioning selections, the group state and the pools'
// weight ESS under it; NULL when a round is left to the R form.
// [[Rcpp::export]]
SEXP emc_ensemble_rounds_latent_cpp(const Rcpp::NumericMatrix& pool_stack,
                                    const Rcpp::NumericVector& base_w,
                                    const Rcpp::IntegerVector& k0, int M, int rounds,
                                    const Rcpp::List& state, SEXP free_entry,
                                    const Rcpp::List& prior, SEXP sem) {
  const int n = k0.size(), p = pool_stack.ncol();
  const R_xlen_t nM = pool_stack.nrow();
  const bool structural = sem != R_NilValue;
  const bool infinite = !structural && free_entry == R_NilValue;
  if (n < 1 || M < 1 || rounds < 1 || nM != static_cast<R_xlen_t>(n) * M ||
      base_w.size() != nM)
    return R_NilValue;
  emc_factor::State fs;
  emc_factor::Prior fp;
  emc_factor::InfiniteState is;
  emc_factor::InfinitePrior ip;
  emc_factor::SemState ss;
  emc_factor::SemPrior sp;
  emc_factor::SemBlocks sb;
  emc_factor::SemMoments mom;
  Mat cov, cov_t;
  std::vector<int> fr;
  if (structural) {
    const Rcpp::List parts(sem);
    ss = sem_state(state);
    if (!sem_prior(prior, Rcpp::as<Rcpp::List>(parts["free"]), sp)) return R_NilValue;
    sb = sem_blocks(Rcpp::as<Rcpp::List>(parts["blocks"]));
    cov = as_mat(Rcpp::as<Rcpp::NumericMatrix>(parts["covariates"]));
    if (cov.r != n) return R_NilValue;
    cov_t = transpose(cov);
  } else if (infinite) {
    is = infinite_state(state);
    ip = infinite_prior(prior);
  } else {
    fs = factor_state(state);
    fp = factor_prior(prior);
    if (!free_entries(Rcpp::LogicalMatrix(free_entry), fr)) return R_NilValue;
  }
  std::vector<int> selected(n), conditioned(n);
  for (int s = 0; s < n; ++s) {
    if (k0[s] < 1 || k0[s] > M) return R_NilValue;
    selected[s] = k0[s] - 1;
  }
  Mat alpha(n, p);
  std::vector<double> lw(nM), wess(n, 0.0), dvar(p), sd(p);
  try {
    for (int r = 1; r <= rounds; ++r) {
      for (int s = 0; s < n; ++s)
        for (int i = 0; i < p; ++i)
          alpha(s, i) = pool_stack(static_cast<R_xlen_t>(s) * M + selected[s], i);
      conditioned = selected;
      const bool drawn =
          structural ? emc_factor::gibbs_sem(ss, alpha, cov, sp) &&
                           emc_factor::finish_sem(ss, cov, sb, mom)
                     : infinite ? emc_factor::gibbs_infinite(is, alpha, ip)
                                : emc_factor::gibbs(fs, alpha, fp, fr);
      if (!drawn) return R_NilValue;
      const Mat& lambda = structural ? ss.lambda : infinite ? is.lambda : fs.lambda;
      const Mat& eta = structural ? ss.eta : infinite ? is.eta : fs.eta;
      const std::vector<double>& mu = structural ? ss.mu : infinite ? is.mu : fs.mu;
      const std::vector<double>& precision =
          structural ? ss.epsilon_inv : infinite ? is.epsilon_inv : fs.sig_err_inv;
      // Each subject's prior mean given its scores (and its covariates);
      // independent errors.
      Mat mean = mult(lambda, transpose(eta));
      if (structural) {
        const Mat direct = mult(ss.K, cov_t);
        for (int s = 0; s < n; ++s)
          for (int i = 0; i < p; ++i) mean(i, s) = mu[i] + direct(i, s) + mean(i, s);
      } else {
        for (int s = 0; s < n; ++s)
          for (int i = 0; i < p; ++i) mean(i, s) = mu[i] + mean(i, s);
      }
      long double log_var = 0.0;
      for (int i = 0; i < p; ++i) {
        dvar[i] = 1 / precision[i];
        sd[i] = std::sqrt(dvar[i]);
        log_var += std::log(dvar[i]);
      }
      volatile double two_pi = 2 * M_PI;
      volatile double scaled = p * std::log(two_pi);
      const double norm = scaled + static_cast<double>(log_var);
      for (R_xlen_t row = 0; row < nM; ++row) {
        const int s = static_cast<int>(row / M);
        long double q = 0.0;
        for (int i = 0; i < p; ++i) {
          const double z = (pool_stack(row, i) - mean(i, s)) / sd[i];
          volatile double z2 = z * z;
          q += z2;
        }
        volatile double density = -0.5 * (norm + static_cast<double>(q));
        lw[row] = base_w[row] + density;
      }
      for (int s = 0; s < n; ++s) {
        const R_xlen_t off = static_cast<R_xlen_t>(s) * M;
        int best = 0;
        double best_v = 0.0;
        bool missing = false;
        for (int j = 0; j < M; ++j) {
          const double g = -std::log(-std::log(unif_rand()));
          const double val = lw[off + j] + g;
          if (ISNAN(val)) missing = true;
          if (j == 0 || best_v < val) {
            best_v = val;
            best = j;
          }
        }
        if (missing) {
          // The rest of the round's uniforms are not drawn: the R form
          // starts again from the stream's state before the rounds.
          return R_NilValue;
        }
        selected[s] = best;
      }
      if (r == rounds) {
        for (int s = 0; s < n; ++s) {
          const R_xlen_t off = static_cast<R_xlen_t>(s) * M;
          double mx = lw[off];
          for (int j = 1; j < M; ++j)
            if (lw[off + j] > mx) mx = lw[off + j];
          long double sw = 0.0, sw2 = 0.0;
          for (int j = 0; j < M; ++j) {
            const double w = std::exp(lw[off + j] - mx);
            volatile double w2 = w * w;
            sw += w;
            sw2 += w2;
          }
          const double a = static_cast<double>(sw), b = static_cast<double>(sw2);
          volatile double a2 = a * a;
          wess[s] = a2 / b;
        }
      }
    }
  } catch (Decline&) {
    return R_NilValue;
  }
  Rcpp::IntegerVector sel(n), cond(n);
  for (int s = 0; s < n; ++s) {
    sel[s] = selected[s] + 1;
    cond[s] = conditioned[s] + 1;
  }
  return Rcpp::List::create(Rcpp::_["selected"] = sel, Rcpp::_["conditioned"] = cond,
                            Rcpp::_["state"] = structural ? sem_result(ss, mom)
                                               : infinite ? infinite_result(is)
                                                          : factor_result(fs),
                            Rcpp::_["wess"] = Rcpp::NumericVector(wess.begin(), wess.end()));
}
