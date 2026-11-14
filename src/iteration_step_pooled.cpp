// The iteration's pooled-coefficient move (.emc_step_pooled, R/pooled.R):
// candidates around a shifted point, each scored over all subjects, one drawn
// in proportion to its posterior weight.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#pragma GCC optimize("fp-contract=off")
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <vector>

#include "iteration_state.h"

#ifndef FCONE
#define FCONE
#endif

using namespace Rcpp;

namespace {

using emc_iter::flag;
using emc_iter::lget;
using emc_iter::lset;
using emc_iter::package_function;
using emc_iter::set_var;
using emc_iter::var;

bool all_finite(const double* x, R_xlen_t n) {
  for (R_xlen_t i = 0; i < n; ++i)
    if (!std::isfinite(x[i])) return false;
  return true;
}

// x %*% y as R's product forms it: the level-2 routine when either side is a
// vector, the level-3 one otherwise, and plain long-double sums when an
// operand is not finite.
void matprod(const double* x, int nrx, int ncx, const double* y, int nry, int ncy, double* z) {
  if (nrx == 0 || ncx == 0 || nry == 0 || ncy == 0) {
    std::fill(z, z + static_cast<R_xlen_t>(nrx) * ncy, 0.0);
    return;
  }
  if (!all_finite(x, static_cast<R_xlen_t>(nrx) * ncx) ||
      !all_finite(y, static_cast<R_xlen_t>(nry) * ncy)) {
    for (int i = 0; i < nrx; ++i)
      for (int k = 0; k < ncy; ++k) {
        long double sum = 0.0;
        for (int j = 0; j < ncx; ++j)
          sum += static_cast<long double>(x[i + static_cast<R_xlen_t>(j) * nrx]) *
                 y[j + static_cast<R_xlen_t>(k) * nry];
        z[i + static_cast<R_xlen_t>(k) * nrx] = static_cast<double>(sum);
      }
    return;
  }
  const double one = 1.0, zero = 0.0;
  const int ione = 1;
  if (ncy == 1)
    F77_CALL(dgemv)("N", &nrx, &ncx, &one, x, &nrx, y, &ione, &zero, z, &ione FCONE);
  else if (nrx == 1)
    F77_CALL(dgemv)("T", &nry, &ncy, &one, y, &nry, x, &ione, &zero, z, &ione FCONE);
  else
    F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z,
                    &nrx FCONE FCONE);
}

// Position of `name` among a vector's names; -1 when absent.
int name_index(SEXP x, const char* name) {
  SEXP nm = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(nm) != STRSXP) return -1;
  for (R_xlen_t i = 0; i < XLENGTH(nm); ++i)
    if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) return static_cast<int>(i);
  return -1;
}

// x with one element raised by `by`.
SEXP bumped(SEXP x, int at, double by) {
  Shield<SEXP> out(Rf_duplicate(x));
  REAL(out)[at] = REAL(out)[at] + by;
  return out;
}

// The curvature grid of one subject (.emc_pooled_fd_rows) at x0 with steps h:
// the centre, each coordinate's two axial points, a forward corner for each
// pair with a subject coordinate, and four central corners for each pair of
// pooled coordinates.
SEXP grid_rows(const std::vector<double>& x0, const std::vector<double>& h, int d, int p,
               SEXP dimnames) {
  const int D = d + p;
  const int M = 1 + 2 * D + (d * (d - 1)) / 2 + d * p + 2 * p * (p - 1);
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, M, D));
  double* o = REAL(out);
  auto base = [&](int r) { for (int q = 0; q < D; ++q) o[r + static_cast<R_xlen_t>(q) * M] = x0[q]; };
  auto move = [&](int r, int q, double sign) {
    volatile double v = sign > 0 ? x0[q] + h[q] : x0[q] - h[q];
    o[r + static_cast<R_xlen_t>(q) * M] = v;
  };
  int r = 0;
  base(r++);
  for (int i = 0; i < D; ++i) {
    base(r); move(r++, i, 1);
    base(r); move(r++, i, -1);
  }
  for (int i = 0; i < d; ++i)
    for (int j = i + 1; j < D; ++j) {
      base(r); move(r, i, 1); move(r++, j, 1);
    }
  for (int i = 0; i + 1 < p; ++i)
    for (int j = i + 1; j < p; ++j) {
      const int ii = d + i, jj = d + j;
      base(r); move(r, ii, 1); move(r++, jj, 1);
      base(r); move(r, ii, 1); move(r++, jj, -1);
      base(r); move(r, ii, -1); move(r++, jj, 1);
      base(r); move(r, ii, -1); move(r++, jj, -1);
    }
  Rf_setAttrib(out, R_DimNamesSymbol, dimnames);
  UNPROTECT(1);
  return out;
}

// One subject's negative Hessian from its grid values added into A
// (.emc_pooled_fd_hessian).
void add_hessian(const double* fs, const std::vector<double>& hs, int d, int p,
                 std::vector<double>& Q, std::vector<double>& A) {
  const int D = d + p;
  std::fill(Q.begin(), Q.end(), 0.0);
  volatile double two = 2, four = 4;
  for (int i = 0; i < D; ++i) {
    volatile double twice = two * fs[0];
    volatile double a = fs[2 * i + 1] - twice;
    volatile double b = a + fs[2 * i + 2];
    volatile double sq = hs[i] * hs[i];
    volatile double q = -b / sq;
    Q[i + static_cast<std::size_t>(i) * D] = q;
  }
  int k = 2 * D + 1;
  for (int i = 0; i < d; ++i)
    for (int j = i + 1; j < D; ++j) {
      volatile double a = fs[k] - fs[2 * i + 1];
      volatile double b = a - fs[2 * j + 1];
      volatile double c = b + fs[0];
      volatile double den = hs[i] * hs[j];
      volatile double q = -c / den;
      Q[i + static_cast<std::size_t>(j) * D] = q;
      Q[j + static_cast<std::size_t>(i) * D] = q;
      ++k;
    }
  for (int i = 0; i + 1 < p; ++i)
    for (int j = i + 1; j < p; ++j) {
      const int ii = d + i, jj = d + j;
      volatile double a = fs[k] - fs[k + 1];
      volatile double b = a - fs[k + 2];
      volatile double c = b + fs[k + 3];
      volatile double d1 = four * hs[ii];
      volatile double den = d1 * hs[jj];
      volatile double q = -c / den;
      Q[ii + static_cast<std::size_t>(jj) * D] = q;
      Q[jj + static_cast<std::size_t>(ii) * D] = q;
      k += 4;
    }
  for (std::size_t e = 0; e < A.size(); ++e) {
    volatile double v = A[e] + Q[e];
    A[e] = v;
  }
}

// (x + t(x)) / 2 in place.
void symmetrise(std::vector<double>& x, int n) {
  std::vector<double> y(x);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      volatile double s = y[i + static_cast<std::size_t>(j) * n] + y[j + static_cast<std::size_t>(i) * n];
      x[i + static_cast<std::size_t>(j) * n] = s / 2;
    }
}

// solve(a, b) for a square a (n x n) and b (n x m): false where solve()
// raises a condition.
bool lu_solve(std::vector<double> a, int n, std::vector<double>& b, int m) {
  if (n < 1 || m < 1) return false;
  std::vector<double> a0(a);
  std::vector<int> ipiv(n);
  int info = 0;
  F77_CALL(dgesv)(&n, &m, a.data(), &n, ipiv.data(), b.data(), &n, &info);
  if (info != 0) return false;
  const double anorm = F77_CALL(dlange)("1", &n, &n, a0.data(), &n, NULL FCONE);
  std::vector<double> work(4 * static_cast<std::size_t>(n));
  double rcond = 0.0;
  F77_CALL(dgecon)("1", &n, a.data(), &n, &anorm, &rcond, work.data(), ipiv.data(), &info FCONE);
  return rcond >= DBL_EPSILON;
}

// A symmetric matrix that is not positive definite, with its eigenvalues
// taken in absolute value and floored relative to the largest. False where
// the decomposition fails.
bool absolute_spectrum(std::vector<double>& A, int n) {
  std::vector<double> a(A), w(n), z(static_cast<std::size_t>(n) * n);
  std::vector<int> isuppz(2 * static_cast<std::size_t>(n));
  int m = 0, il = 0, iu = 0, info = 0, lwork = -1, liwork = -1, iw0 = 0;
  double vl = 0.0, vu = 0.0, abstol = 0.0, w0 = 0.0;
  F77_CALL(dsyevr)("V", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.data(), &n, isuppz.data(), &w0, &lwork, &iw0, &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) return false;
  lwork = static_cast<int>(w0);
  liwork = iw0;
  std::vector<double> work(lwork);
  std::vector<int> iwork(liwork);
  F77_CALL(dsyevr)("V", "A", "L", &n, a.data(), &n, &vl, &vu, &il, &iu, &abstol, &m, w.data(),
                   z.data(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0) return false;
  // Values and vectors in decreasing order, as eigen() returns them.
  std::vector<double> V(z.size()), lam(n);
  double top = 0.0;
  for (int c = 0; c < n; ++c) top = std::max(top, std::fabs(w[c]));
  volatile double floor_at = top * DBL_EPSILON;
  for (int c = 0; c < n; ++c) {
    const int src = n - 1 - c;
    lam[c] = std::max(std::fabs(w[src]), static_cast<double>(floor_at));
    std::copy(z.begin() + static_cast<std::size_t>(src) * n,
              z.begin() + static_cast<std::size_t>(src + 1) * n,
              V.begin() + static_cast<std::size_t>(c) * n);
  }
  std::vector<double> right(z.size());
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      volatile double v = lam[i] * V[j + static_cast<std::size_t>(i) * n];
      right[i + static_cast<std::size_t>(j) * n] = v;
    }
  matprod(V.data(), n, n, right.data(), n, n, A.data());
  return true;
}

// The shear slope (d x p) and step covariance (p x p) from the summed
// curvature A (.emc_pooled_curvature_from's algebra). False where a step of
// it fails, which the R form then reports.
bool shear_and_cov(std::vector<double> A, int d, int p, std::vector<double>& shear,
                   std::vector<double>& cov) {
  const int D = d + p;
  if (!all_finite(A.data(), static_cast<R_xlen_t>(A.size()))) return false;
  {
    std::vector<double> U(A);
    int info = 0, n = D;
    F77_CALL(dpotrf)("U", &n, U.data(), &n, &info FCONE);
    if (info != 0 && !absolute_spectrum(A, D)) return false;
  }
  std::vector<double> aa(static_cast<std::size_t>(d) * d), ap(static_cast<std::size_t>(d) * p),
      pp(static_cast<std::size_t>(p) * p);
  for (int j = 0; j < d; ++j)
    for (int i = 0; i < d; ++i) aa[i + static_cast<std::size_t>(j) * d] = A[i + static_cast<std::size_t>(j) * D];
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < d; ++i) ap[i + static_cast<std::size_t>(j) * d] = A[i + static_cast<std::size_t>(d + j) * D];
    for (int i = 0; i < p; ++i) pp[i + static_cast<std::size_t>(j) * p] = A[d + i + static_cast<std::size_t>(d + j) * D];
  }
  shear = ap;
  if (!lu_solve(aa, d, shear, p)) return false;
  for (double& v : shear) v = -v;
  // crossprod() of the cross block with the slope.
  std::vector<double> cross(static_cast<std::size_t>(p) * p);
  if (!all_finite(ap.data(), static_cast<R_xlen_t>(ap.size())) ||
      !all_finite(shear.data(), static_cast<R_xlen_t>(shear.size())))
    return false;
  {
    const double one = 1.0, zero = 0.0;
    const int ione = 1;
    if (p == 1)
      F77_CALL(dgemv)("T", &d, &p, &one, ap.data(), &d, shear.data(), &ione, &zero, cross.data(),
                      &ione FCONE);
    else
      F77_CALL(dgemm)("T", "N", &p, &p, &d, &one, ap.data(), &d, shear.data(), &d, &zero,
                      cross.data(), &p FCONE FCONE);
  }
  for (std::size_t e = 0; e < pp.size(); ++e) {
    volatile double v = pp[e] + cross[e];
    pp[e] = v;
  }
  symmetrise(pp, p);
  cov.assign(static_cast<std::size_t>(p) * p, 0.0);
  for (int i = 0; i < p; ++i) cov[i + static_cast<std::size_t>(i) * p] = 1.0;
  return lu_solve(pp, p, cov, p);
}

}  // namespace

namespace emc_iter {

bool step_pooled(SEXP st, SEXP cx, int) {
  SEXP pmwgs = var(st, "pmwgs"), proposals = var(st, "proposals"), pars = var(st, "pars");
  if (TYPEOF(pmwgs) != VECSXP || TYPEOF(pars) != VECSXP || flag(lget(cx, "sampler_profile")))
    return false;
  SEXP pst = lget(pmwgs, "pooled"), stage_s = lget(cx, "stage");
  if (TYPEOF(pst) != VECSXP || TYPEOF(stage_s) != STRSXP || XLENGTH(stage_s) != 1) return false;
  const char* stage = CHAR(STRING_ELT(stage_s, 0));
  const bool sampling = std::strcmp(stage, "sample") == 0;
  SEXP value = lget(pst, "value"), tries = lget(pst, "tries"), moves = lget(pst, "moves");
  SEXP n_adapt = lget(pst, "n_adapt"), prior_mean = lget(pst, "prior_mean");
  SEXP prior_var = lget(pst, "prior_var"), par_names = lget(pmwgs, "par_names");
  SEXP nm = Rf_getAttrib(value, R_NamesSymbol);
  const int n_pars = Rf_asInteger(lget(pmwgs, "n_pars"));
  const int n = Rf_asInteger(lget(pmwgs, "n_subjects"));
  if (TYPEOF(value) != REALSXP || TYPEOF(nm) != STRSXP || XLENGTH(value) < 1) return false;
  const int p = static_cast<int>(XLENGTH(value));
  // A block's first move outside sampling measures the curvature first.
  bool measure = false;
  SEXP blk = lget(pst, "block");
  if (!sampling) {
    SEXP cb = lget(pst, "curvature_block");
    const bool same = Rf_length(cb) == 1 && Rf_length(blk) == 1 && Rf_isNumeric(cb) &&
                      Rf_isNumeric(blk) && !ISNAN(Rf_asReal(cb)) && !ISNAN(Rf_asReal(blk)) &&
                      Rf_asReal(cb) == Rf_asReal(blk);
    measure = !same;
    if (measure && blk == R_NilValue) return false;
  }
  const int at_tries = TYPEOF(tries) == REALSXP ? name_index(tries, stage) : -1;
  const int at_moves = TYPEOF(moves) == REALSXP ? name_index(moves, stage) : -1;
  if (at_tries < 0 || at_moves < 0 ||
      (TYPEOF(n_adapt) != INTSXP && TYPEOF(n_adapt) != REALSXP) || XLENGTH(n_adapt) != 1 ||
      TYPEOF(prior_mean) != REALSXP || XLENGTH(prior_mean) != p ||
      TYPEOF(prior_var) != REALSXP || XLENGTH(prior_var) != p ||
      TYPEOF(par_names) != STRSXP || XLENGTH(par_names) != n_pars || n_pars < 1 || n < 1 ||
      TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) ||
      Rf_nrows(proposals) != n_pars + 1 || Rf_ncols(proposals) != n ||
      !all_finite(REAL(value), p) || !all_finite(REAL(proposals), XLENGTH(proposals)))
    return false;
  static Function score = package_function(".emc_pooled_score");

  const double* prop = REAL(proposals);
  const R_xlen_t rowlen = n_pars + 1;
  Shield<SEXP> alpha(Rf_allocMatrix(REALSXP, n_pars, n));
  Shield<SEXP> prev_ll(Rf_allocVector(REALSXP, n));
  for (int s = 0; s < n; ++s) {
    for (int q = 0; q < n_pars; ++q) REAL(alpha)[q + static_cast<R_xlen_t>(s) * n_pars] = prop[q + s * rowlen];
    REAL(prev_ll)[s] = prop[n_pars + s * rowlen];
  }
  {
    Shield<SEXP> dn(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn, 0, par_names);
    Rf_setAttrib(alpha, R_DimNamesSymbol, dn);
    SEXP pdn = Rf_getAttrib(proposals, R_DimNamesSymbol);
    if (TYPEOF(pdn) == VECSXP && XLENGTH(pdn) == 2 && VECTOR_ELT(pdn, 1) != R_NilValue)
      Rf_setAttrib(prev_ll, R_NamesSymbol, VECTOR_ELT(pdn, 1));
  }
  Shield<SEXP> col_names(Rf_allocVector(STRSXP, n_pars + p));
  for (int q = 0; q < n_pars; ++q) SET_STRING_ELT(col_names, q, STRING_ELT(par_names, q));
  for (int i = 0; i < p; ++i) SET_STRING_ELT(col_names, n_pars + i, STRING_ELT(nm, i));
  Shield<SEXP> row_dn(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(row_dn, 1, col_names);
  // The candidates' scores over the subjects on the pool, which may lose
  // workers on the way.
  SEXP wpool = var(st, "wpool");
  RObject pool = wpool;
  SEXP wpool_ctx = var(st, "wpool_ctx");
  if (wpool_ctx == R_NilValue) wpool_ctx = lget(cx, "wpool_ctx");
  auto scores_of = [&](SEXP rows_in) -> SEXP {
    Shield<SEXP> scored(score(pmwgs, rows_in, alpha, pool, wpool_ctx, var(st, "wpool_part"), pars,
                              var(st, "group_chol_it"), var(st, "pm_settings"), prev_ll,
                              var(st, "streams"), lget(cx, "n_cores"),
                              _["ll_cache_validated"] = true));
    SEXP alive = lget(scored, "alive");
    if (static_cast<SEXP>(pool) != R_NilValue && alive != R_NilValue)
      pool = lset(pool, "alive", alive);
    return lget(scored, "ll");
  };

  RObject measured;
  if (measure) {
    // .emc_pooled_curvature(): each subject's grid around its state and the
    // pooled value, scored, and the shear and step covariance its summed
    // curvature gives.
    const int D = n_pars + p;
    volatile double h0 = 1e-3;
    std::vector<std::vector<double>> steps(n, std::vector<double>(D));
    Shield<SEXP> grids(Rf_allocVector(VECSXP, n));
    std::vector<double> x0(D);
    for (int s = 0; s < n; ++s) {
      for (int q = 0; q < n_pars; ++q) x0[q] = REAL(alpha)[q + static_cast<R_xlen_t>(s) * n_pars];
      for (int q = 0; q < p; ++q) x0[n_pars + q] = REAL(value)[q];
      for (int q = 0; q < D; ++q) {
        volatile double v = h0 * std::max(1.0, std::fabs(x0[q]));
        steps[s][q] = v;
      }
      SET_VECTOR_ELT(grids, s, grid_rows(x0, steps[s], n_pars, p, row_dn));
    }
    const int M = Rf_nrows(VECTOR_ELT(grids, 0));
    Shield<SEXP> f(scores_of(grids));
    if (TYPEOF(f) != REALSXP || !Rf_isMatrix(f) || Rf_nrows(f) != M || Rf_ncols(f) != n)
      Rcpp::stop("pooled scores do not match the curvature grid");
    RObject pst_new = pst;
    if (all_finite(REAL(f), XLENGTH(f))) {
      static Function mu_prec = package_function(".emc_pooled_mu_prec");
      Shield<SEXP> P(mu_prec(lget(pmwgs, "prior"), n_pars));
      std::vector<double> A(static_cast<std::size_t>(D) * D, 0.0), Q(A.size());
      for (int s = 0; s < n; ++s)
        add_hessian(REAL(f) + static_cast<R_xlen_t>(s) * M, steps[s], n_pars, p, Q, A);
      std::vector<double> shear_v, cov_v;
      bool ok = TYPEOF(P) == REALSXP && Rf_isMatrix(P) && Rf_nrows(P) == n_pars &&
                Rf_ncols(P) == n_pars;
      if (ok) {
        volatile double unit = 1, none = 0;
        for (int j = 0; j < n_pars; ++j)
          for (int i = 0; i < n_pars; ++i) {
            volatile double v = A[i + static_cast<std::size_t>(j) * D] +
                                REAL(P)[i + static_cast<R_xlen_t>(j) * n_pars];
            A[i + static_cast<std::size_t>(j) * D] = v;
          }
        for (int j = 0; j < p; ++j)
          for (int i = 0; i < p; ++i) {
            volatile double add = i == j ? unit / REAL(prior_var)[i] : none;
            volatile double v = A[n_pars + i + static_cast<std::size_t>(n_pars + j) * D] + add;
            A[n_pars + i + static_cast<std::size_t>(n_pars + j) * D] = v;
          }
        symmetrise(A, D);
        ok = shear_and_cov(A, n_pars, p, shear_v, cov_v);
      }
      if (ok) {
        Shield<SEXP> sh(Rf_allocMatrix(REALSXP, n_pars, p)), cv(Rf_allocMatrix(REALSXP, p, p));
        std::copy(shear_v.begin(), shear_v.end(), REAL(sh));
        std::copy(cov_v.begin(), cov_v.end(), REAL(cv));
        pst_new = lset(pst_new, "shear", sh);
        pst_new = lset(pst_new, "cov", cv);
        if (Rf_asReal(n_adapt) == 0) pst_new = lset(pst_new, "log_eps", Shield<SEXP>(Rf_ScalarReal(0)));
      } else {
        // The R form of the algebra, which reports what failed.
        static Function from = package_function(".emc_pooled_curvature_from");
        Shield<SEXP> hs(Rf_allocVector(VECSXP, n));
        for (int s = 0; s < n; ++s)
          SET_VECTOR_ELT(hs, s, NumericVector(steps[s].begin(), steps[s].end()));
        pst_new = from(pmwgs, pst, f, hs, n_pars);
      }
    }
    measured = lset(pst_new, "curvature_block", blk);
    pst = measured;
  }

  SEXP cov = lget(pst, "cov"), log_eps = lget(pst, "log_eps"), shear = lget(pst, "shear");
  const bool sheared = shear != R_NilValue;
  if (TYPEOF(cov) != REALSXP || !Rf_isMatrix(cov) || Rf_nrows(cov) != p || Rf_ncols(cov) != p ||
      TYPEOF(log_eps) != REALSXP || XLENGTH(log_eps) != 1 ||
      (sheared && (TYPEOF(shear) != REALSXP || !Rf_isMatrix(shear) ||
                   Rf_nrows(shear) != n_pars || Rf_ncols(shear) != p)) ||
      !all_finite(REAL(cov), XLENGTH(cov)) || !all_finite(REAL(log_eps), 1))
    return false;
  SEXP tmu = lget(pars, "tmu"), subj_mu = lget(pars, "subj_mu");
  const bool subj_shift = sheared && TYPEOF(subj_mu) == REALSXP && Rf_isMatrix(subj_mu);
  if (sheared && (TYPEOF(tmu) != REALSXP || XLENGTH(tmu) < n_pars ||
                  !all_finite(REAL(shear), XLENGTH(shear)) ||
                  (Rf_isMatrix(subj_mu) && !subj_shift) ||
                  (subj_shift && Rf_nrows(subj_mu) != n_pars)))
    return false;
  static Function options = package_function(".emc_pooled_options");
  const int K = Rf_asInteger(lget(Shield<SEXP>(options()), "candidates"));
  const double target = Rf_asReal(lget(Shield<SEXP>(options()), "target"));
  if (K < 1) return false;

  // The proposal's lower factor.
  std::vector<double> L(static_cast<std::size_t>(p) * p, 0.0);
  {
    const double scale = std::exp(2 * REAL(log_eps)[0]);
    std::vector<double> U(static_cast<std::size_t>(p) * p, 0.0);
    for (int j = 0; j < p; ++j)
      for (int i = 0; i <= j; ++i) {
        volatile double v = scale * REAL(cov)[i + static_cast<R_xlen_t>(j) * p];
        U[i + static_cast<std::size_t>(j) * p] = v;
      }
    int info = 0;
    F77_CALL(dpotrf)("U", &p, U.data(), &p, &info FCONE);
    if (info != 0) return false;
    for (int j = 0; j < p; ++j)
      for (int i = 0; i <= j; ++i)
        L[j + static_cast<std::size_t>(i) * p] = U[i + static_cast<std::size_t>(j) * p];
  }
  // The mean's prior precision under a shear, which moves the group mean too.
  RObject prec;
  if (sheared) {
    static Function mu_prec = package_function(".emc_pooled_mu_prec");
    prec = mu_prec(lget(pmwgs, "prior"), n_pars);
    SEXP pm = lget(lget(pmwgs, "prior"), "theta_mu_mean");
    if (TYPEOF(prec) != REALSXP || !Rf_isMatrix(prec) || Rf_nrows(prec) != n_pars ||
        Rf_ncols(prec) != n_pars || TYPEOF(pm) != REALSXP || XLENGTH(pm) < n_pars)
      return false;
  }
  static Function set_data = package_function(".emc_pooled_set_data");
  static Function rejects = package_function(".emc_step_pooled_rejects");
  static Function sample_int("sample.int");

  SEXP seed = lget(lget(pmwgs, "rng"), "gibbs");
  if (seed != R_NilValue) set_global_seed(seed);
  std::vector<double> e(p), E(static_cast<std::size_t>(p) * K);
  GetRNGstate();
  for (int i = 0; i < p; ++i) e[i] = norm_rand();
  for (std::size_t i = 0; i < E.size(); ++i) E[i] = norm_rand();
  PutRNGstate();

  // The shifted centre and the candidates around it, one per row.
  std::vector<double> z(p), LE(static_cast<std::size_t>(p) * K);
  matprod(L.data(), p, p, e.data(), p, 1, z.data());
  for (int i = 0; i < p; ++i) {
    volatile double v = REAL(value)[i] + z[i];
    z[i] = v;
  }
  matprod(L.data(), p, p, E.data(), p, K, LE.data());
  std::vector<double> cands(static_cast<std::size_t>(K) * p), delta(cands.size());
  for (int k = 0; k < K; ++k)
    for (int i = 0; i < p; ++i) {
      volatile double v = z[i] + LE[i + static_cast<std::size_t>(k) * p];
      cands[k + static_cast<std::size_t>(i) * K] = v;
      volatile double d = v - REAL(value)[i];
      delta[k + static_cast<std::size_t>(i) * K] = d;
    }
  // What each candidate adds to every subject's state.
  std::vector<double> shift(static_cast<std::size_t>(K) * n_pars, 0.0);
  if (sheared) {
    std::vector<double> ts(static_cast<std::size_t>(p) * n_pars);
    for (int q = 0; q < n_pars; ++q)
      for (int i = 0; i < p; ++i)
        ts[i + static_cast<std::size_t>(q) * p] = REAL(shear)[q + static_cast<R_xlen_t>(i) * n_pars];
    matprod(delta.data(), K, p, ts.data(), p, n_pars, shift.data());
  }

  Shield<SEXP> rows(Rf_allocVector(VECSXP, n));
  for (int s = 0; s < n; ++s) {
    SEXP r = Rf_allocMatrix(REALSXP, K, n_pars + p);
    SET_VECTOR_ELT(rows, s, r);
    double* rp = REAL(r);
    for (int q = 0; q < n_pars; ++q)
      for (int k = 0; k < K; ++k) {
        volatile double v = shift[k + static_cast<std::size_t>(q) * K] +
                            REAL(alpha)[q + static_cast<R_xlen_t>(s) * n_pars];
        rp[k + static_cast<R_xlen_t>(q) * K] = v;
      }
    std::copy(cands.begin(), cands.end(), rp + static_cast<R_xlen_t>(n_pars) * K);
    Rf_setAttrib(r, R_DimNamesSymbol, row_dn);
  }

  Shield<SEXP> ll(scores_of(rows));
  if (TYPEOF(ll) != REALSXP || !Rf_isMatrix(ll) || Rf_nrows(ll) != K || Rf_ncols(ll) != n)
    Rcpp::stop("pooled scores do not match the candidates");

  // Log weights: the current value first, then the candidates.
  std::vector<double> lw(K + 1);
  for (int c = 0; c <= K; ++c) {
    long double sum = 0.0;
    for (int i = 0; i < p; ++i) {
      const double x = c == 0 ? REAL(value)[i] : cands[(c - 1) + static_cast<std::size_t>(i) * K];
      sum += Rf_dnorm4(x, REAL(prior_mean)[i], std::sqrt(REAL(prior_var)[i]), 1);
    }
    lw[c] = static_cast<double>(sum);
  }
  if (sheared) {
    const double* pm = REAL(lget(lget(pmwgs, "prior"), "theta_mu_mean"));
    std::vector<double> r(n_pars), Pr(n_pars);
    for (int c = 0; c <= K; ++c) {
      for (int q = 0; q < n_pars; ++q) {
        const double m = c == 0 ? 0.0 : shift[(c - 1) + static_cast<std::size_t>(q) * K];
        volatile double mu = REAL(tmu)[q] + m;
        volatile double d = mu - pm[q];
        r[q] = d;
      }
      matprod(REAL(prec), n_pars, n_pars, r.data(), n_pars, 1, Pr.data());
      long double sum = 0.0;
      for (int q = 0; q < n_pars; ++q) {
        volatile double t = r[q] * Pr[q];
        sum += t;
      }
      volatile double lpm = -0.5 * static_cast<double>(sum);
      volatile double both = lw[c] + lpm;
      lw[c] = both;
    }
  }
  {
    long double sum = 0.0;
    for (int s = 0; s < n; ++s) sum += REAL(prev_ll)[s];
    volatile double v = lw[0] + static_cast<double>(sum);
    lw[0] = v;
    std::vector<long double> rs(K, 0.0L);
    for (int s = 0; s < n; ++s)
      for (int k = 0; k < K; ++k) rs[k] += REAL(ll)[k + static_cast<R_xlen_t>(s) * K];
    for (int k = 0; k < K; ++k) {
      volatile double w = lw[k + 1] + static_cast<double>(rs[k]);
      lw[k + 1] = w;
    }
  }
  bool any = false;
  double top = R_NegInf;
  for (int c = 0; c <= K; ++c) {
    if (!std::isfinite(lw[c])) lw[c] = R_NegInf;
    else any = true;
    if (lw[c] > top) top = lw[c];
  }
  int pick = 1;
  if (any) {
    NumericVector prob(K + 1);
    for (int c = 0; c <= K; ++c) prob[c] = std::exp(lw[c] - top);
    pick = Rf_asInteger(Shield<SEXP>(sample_int(K + 1, 1, _["prob"] = prob)));
  }
  const bool moved = pick > 1;
  const int k = pick - 2;

  RObject pst2 = pst, props2 = proposals, pars2 = pars;
  if (moved) {
    Shield<SEXP> v2(Rf_duplicate(value));
    for (int i = 0; i < p; ++i) REAL(v2)[i] = cands[k + static_cast<std::size_t>(i) * K];
    pst2 = lset(pst2, "value", v2);
    Shield<SEXP> pr(Rf_duplicate(proposals));
    double* pp = REAL(pr);
    for (int s = 0; s < n; ++s) {
      for (int q = 0; q < n_pars; ++q) {
        volatile double v = REAL(alpha)[q + static_cast<R_xlen_t>(s) * n_pars] +
                            shift[k + static_cast<std::size_t>(q) * K];
        pp[q + s * rowlen] = v;
      }
      pp[n_pars + s * rowlen] = REAL(ll)[k + static_cast<R_xlen_t>(s) * K];
    }
    props2 = pr;
    if (sheared) {
      Shield<SEXP> tmu2(Rf_duplicate(tmu));
      for (int q = 0; q < n_pars; ++q) {
        volatile double v = REAL(tmu2)[q] + shift[k + static_cast<std::size_t>(q) * K];
        REAL(tmu2)[q] = v;
      }
      pars2 = lset(pars2, "tmu", tmu2);
      if (subj_shift) {
        Shield<SEXP> sm(Rf_duplicate(subj_mu));
        for (R_xlen_t i = 0; i < XLENGTH(sm); ++i) {
          volatile double v = REAL(sm)[i] + shift[k + static_cast<std::size_t>(i % n_pars) * K];
          REAL(sm)[i] = v;
        }
        pars2 = lset(pars2, "subj_mu", sm);
      }
      Shield<SEXP> a2(emc_iter::proposal_rows(pr, n_pars));
      pars2 = lset(pars2, "alpha", a2);
    }
  }
  pst2 = lset(pst2, "tries", Shield<SEXP>(bumped(tries, at_tries, 1.0)));
  pst2 = lset(pst2, "moves", Shield<SEXP>(bumped(moves, at_moves, moved ? 1.0 : 0.0)));
  if (!sampling) {
    double count;
    if (TYPEOF(n_adapt) == INTSXP) {
      const int c = INTEGER(n_adapt)[0] + 1;
      pst2 = lset(pst2, "n_adapt", Shield<SEXP>(Rf_ScalarInteger(c)));
      count = c;
    } else {
      count = REAL(n_adapt)[0] + 1;
      pst2 = lset(pst2, "n_adapt", Shield<SEXP>(Rf_ScalarReal(count)));
    }
    // Divisions as written: a constant divisor must not become a product.
    volatile double ten = 10, unit = 1;
    volatile double tenth = count / ten;
    volatile double base = 1 + tenth;
    volatile double power = R_pow(base, 0.6);
    volatile double gamma = unit / power;
    volatile double step = gamma * ((moved ? 1.0 : 0.0) - target);
    volatile double le = REAL(log_eps)[0] + step;
    Shield<SEXP> le2(Rf_duplicate(log_eps));
    REAL(le2)[0] = le;
    pst2 = lset(pst2, "log_eps", le2);
  }
  Shield<SEXP> data2(set_data(lget(pmwgs, "data"), lget(pst2, "value")));
  Shield<SEXP> seed_after(global_seed());
  SEXP rng = lget(pmwgs, "rng");
  Shield<SEXP> rng2(rng == R_NilValue
                        ? static_cast<SEXP>(List::create(_["gibbs"] = seed_after))
                        : lset(rng, "gibbs", seed_after));
  RObject pm2 = lset(pmwgs, "pooled", pst2);
  pm2 = lset(pm2, "data", data2);
  pm2 = lset(pm2, "rng", rng2);

  set_var(st, "proposals", props2);
  set_var(st, "pars", pars2);
  if (static_cast<SEXP>(pool) != wpool) set_var(st, "wpool", pool);
  set_var(st, "data", data2);
  set_var(st, "wpool_ctx",
          Shield<SEXP>(wpool_ctx == R_NilValue
                           ? static_cast<SEXP>(List::create(_["data"] = data2))
                           : lset(wpool_ctx, "data", data2)));
  set_var(st, "pmwgs", pm2);
  rejects(st);
  return true;
}

}  // namespace emc_iter
