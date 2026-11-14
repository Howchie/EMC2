// The group draws of the factor hierarchies: parameter-expanded factor
// analysis with loadings fixed or free by entry (R/variant_factor.R
// gibbs_step_factor) and the shrinkage prior over an unbounded number of
// factors (R/variant_infnt_factor.R gibbs_step_infnt_factor).
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
  if (x.c != y.r || x.r == 0 || x.c == 0 || y.c == 0) throw Decline();
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

}  // namespace emc_factor

// gibbs_step_infnt_factor()'s draw: `alpha` subjects x parameters, the state
// as the store holds it and the prior's vectors at full length (`ad2` NULL
// for a prior without it). NULL when the draw is left to the R form.
// [[Rcpp::export]]
SEXP emc_gibbs_infnt_factor_cpp(const Rcpp::NumericMatrix& alpha,
                                const Rcpp::NumericMatrix& psi,
                                const Rcpp::NumericVector& epsilon_inv,
                                const Rcpp::NumericMatrix& lambda,
                                const Rcpp::NumericVector& mu,
                                const Rcpp::NumericVector& delta,
                                const Rcpp::NumericVector& mu_invar,
                                const Rcpp::NumericVector& mu_mean, double df, double ad1,
                                double bd1, SEXP ad2, double bd2,
                                const Rcpp::NumericVector& shape_err,
                                const Rcpp::NumericVector& rate_err) {
  emc_factor::InfiniteState s;
  s.psi = as_mat(psi);
  s.lambda = as_mat(lambda);
  s.mu = as_std(mu);
  s.epsilon_inv = as_std(epsilon_inv);
  s.delta = as_std(delta);
  emc_factor::InfinitePrior prior;
  prior.mu_invar = as_std(mu_invar);
  prior.mu_mean = as_std(mu_mean);
  prior.shape_err = as_std(shape_err);
  prior.rate_err = as_std(rate_err);
  prior.df = df;
  prior.ad1 = ad1;
  prior.bd1 = bd1;
  prior.has_ad2 = ad2 != R_NilValue;
  prior.ad2 = prior.has_ad2 ? Rf_asReal(ad2) : 0.0;
  prior.bd2 = bd2;
  if (!emc_factor::gibbs_infinite(s, as_mat(alpha), prior)) return R_NilValue;
  return Rcpp::List::create(
      Rcpp::_["mu"] = Rcpp::NumericVector(s.mu.begin(), s.mu.end()),
      Rcpp::_["var"] = as_r(s.var), Rcpp::_["lambda"] = as_r(s.lambda),
      Rcpp::_["eta"] = as_r(s.eta),
      Rcpp::_["epsilon_inv"] = Rcpp::NumericVector(s.epsilon_inv.begin(), s.epsilon_inv.end()),
      Rcpp::_["psi"] = as_r(s.psi),
      Rcpp::_["delta"] = Rcpp::NumericVector(s.delta.begin(), s.delta.end()));
}

// gibbs_step_factor()'s draw: `alpha` subjects x parameters, the state as the
// store holds it and the prior's vectors at full length. NULL when the draw
// is left to the R form, which then starts again from the stream's state
// before the call.
// [[Rcpp::export]]
SEXP emc_gibbs_factor_cpp(const Rcpp::NumericMatrix& alpha, const Rcpp::NumericMatrix& eta,
                          const Rcpp::NumericVector& psi_inv,
                          const Rcpp::NumericVector& sig_err_inv,
                          const Rcpp::NumericMatrix& lambda, const Rcpp::NumericVector& mu,
                          const Rcpp::LogicalMatrix& free_entry,
                          const Rcpp::NumericMatrix& mu_invar,
                          const Rcpp::NumericVector& mu_mean,
                          const Rcpp::NumericVector& shape_err,
                          const Rcpp::NumericVector& rate_err,
                          const Rcpp::NumericVector& lambda_invar,
                          const Rcpp::NumericVector& shape_psi,
                          const Rcpp::NumericVector& rate_psi) {
  emc_factor::State s;
  s.eta = as_mat(eta);
  s.lambda = as_mat(lambda);
  s.mu = as_std(mu);
  s.sig_err_inv = as_std(sig_err_inv);
  s.psi_inv = as_std(psi_inv);
  emc_factor::Prior prior;
  prior.mu_invar = as_mat(mu_invar);
  prior.mu_mean = as_std(mu_mean);
  prior.shape_err = as_std(shape_err);
  prior.rate_err = as_std(rate_err);
  prior.lambda_invar = as_std(lambda_invar);
  prior.shape_psi = as_std(shape_psi);
  prior.rate_psi = as_std(rate_psi);
  std::vector<int> fr(free_entry.size());
  for (R_xlen_t i = 0; i < free_entry.size(); ++i) {
    if (free_entry[i] == NA_LOGICAL) return R_NilValue;
    fr[i] = free_entry[i] == TRUE;
  }
  if (!emc_factor::gibbs(s, as_mat(alpha), prior, fr)) return R_NilValue;
  return Rcpp::List::create(
      Rcpp::_["mu"] = Rcpp::NumericVector(s.mu.begin(), s.mu.end()),
      Rcpp::_["var"] = as_r(s.var), Rcpp::_["lambda"] = as_r(s.lambda),
      Rcpp::_["lambda_scaled"] = as_r(s.lambda_scaled), Rcpp::_["eta"] = as_r(s.eta),
      Rcpp::_["sig_err_inv"] = Rcpp::NumericVector(s.sig_err_inv.begin(), s.sig_err_inv.end()),
      Rcpp::_["psi_inv"] = Rcpp::NumericVector(s.psi_inv.begin(), s.psi_inv.end()));
}
