// Bookkeeping of the stored Laplace checks (R/subject_fit_cache.R) and of the
// shear move's window (R/shear_move.R).
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>

namespace {

// A list with names and no other attribute; the position of `name` in it, or
// -1. False for anything else.
bool named_list(SEXP x, const char* name, R_xlen_t& at) {
  if (TYPEOF(x) != VECSXP || Rf_isObject(x)) return false;
  SEXP names = Rf_getAttrib(x, R_NamesSymbol);
  if (TYPEOF(names) != STRSXP) return false;
  if (Rf_getAttrib(x, R_DimSymbol) != R_NilValue) return false;
  at = -1;
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      if (at >= 0) return false;
      at = i;
    }
  return true;
}

}  // namespace

// .emc_laplace_check_store()'s entries after storing `check` under `key`:
// entries of the same key dropped (the last one's draws kept when none are
// given), the others' draws dropped, the new entry appended and the list cut
// to its last `max_entries`. NULL for entries that are not an unnamed list of
// plain named lists, which the R form handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_laplace_check_store_cpp(SEXP entries, SEXP key, SEXP check, SEXP draws,
                                 int max_entries) {
  if (entries != R_NilValue &&
      (TYPEOF(entries) != VECSXP || Rf_isObject(entries) ||
       Rf_getAttrib(entries, R_NamesSymbol) != R_NilValue ||
       Rf_getAttrib(entries, R_DimSymbol) != R_NilValue))
    return R_NilValue;
  if (max_entries < 1) return R_NilValue;
  const R_xlen_t n = entries == R_NilValue ? 0 : XLENGTH(entries);
  Rcpp::List kept(n);
  R_xlen_t n_kept = 0;
  Rcpp::RObject carried = draws;
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP e = VECTOR_ELT(entries, i);
    R_xlen_t at_key, at_draws;
    if (!named_list(e, "key", at_key) || !named_list(e, "draws", at_draws)) return R_NilValue;
    SEXP e_key = at_key < 0 ? R_NilValue : VECTOR_ELT(e, at_key);
    if (R_compute_identical(e_key, key, 16)) {
      if (draws == R_NilValue) carried = at_draws < 0 ? R_NilValue : VECTOR_ELT(e, at_draws);
      continue;
    }
    if (at_draws < 0) {
      kept[n_kept++] = e;
      continue;
    }
    const R_xlen_t m = XLENGTH(e);
    SEXP names = Rf_getAttrib(e, R_NamesSymbol);
    Rcpp::List bare(m - 1);
    Rcpp::CharacterVector bare_names(m - 1);
    for (R_xlen_t j = 0, k = 0; j < m; ++j) {
      if (j == at_draws) continue;
      bare[k] = VECTOR_ELT(e, j);
      bare_names[k++] = STRING_ELT(names, j);
    }
    bare.attr("names") = bare_names;
    kept[n_kept++] = bare;
  }
  Rcpp::List entry = carried.isNULL()
    ? Rcpp::List::create(Rcpp::Named("key") = key, Rcpp::Named("check") = check)
    : Rcpp::List::create(Rcpp::Named("key") = key, Rcpp::Named("check") = check,
                         Rcpp::Named("draws") = carried);
  const R_xlen_t total = n_kept + 1;
  const R_xlen_t first = total > max_entries ? total - max_entries : 0;
  Rcpp::List out(total - first);
  for (R_xlen_t i = first; i < total; ++i)
    out[i - first] = i < n_kept ? static_cast<SEXP>(kept[i]) : static_cast<SEXP>(entry);
  return out;
}

// .emc_shear_window_lsd(): half the log of the diagonals of the covariance
// draws numbered idx (parameters x draws). NULL for anything but a cube of
// doubles with square slices and positions within it, which the R form
// handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_window_lsd_cpp(SEXP theta_var, SEXP idx, int n_pars) {
  SEXP dim = Rf_getAttrib(theta_var, R_DimSymbol);
  if (TYPEOF(theta_var) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 3 ||
      (TYPEOF(idx) != INTSXP && TYPEOF(idx) != REALSXP))
    return R_NilValue;
  const int p = INTEGER(dim)[0], n_draws = INTEGER(dim)[2];
  if (p < 1 || INTEGER(dim)[1] != p || n_pars != p) return R_NilValue;
  const R_xlen_t n = XLENGTH(idx);
  if (n < 1) return R_NilValue;
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, p, static_cast<int>(n)));
  const double* v = REAL(theta_var);
  const R_xlen_t slice = static_cast<R_xlen_t>(p) * p;
  for (R_xlen_t j = 0; j < n; ++j) {
    double at;
    if (TYPEOF(idx) == INTSXP) {
      if (INTEGER(idx)[j] == NA_INTEGER) { UNPROTECT(1); return R_NilValue; }
      at = INTEGER(idx)[j];
    } else {
      at = REAL(idx)[j];
    }
    if (!(at >= 1) || !(at < n_draws + 1.0)) { UNPROTECT(1); return R_NilValue; }
    const R_xlen_t d = static_cast<R_xlen_t>(at) - 1;
    for (int i = 0; i < p; ++i) {
      volatile double l = std::log(v[d * slice + i + static_cast<R_xlen_t>(i) * p]);
      REAL(out)[j * p + i] = 0.5 * l;
    }
  }
  UNPROTECT(1);
  return out;
}

// One bounded parameter type's half-spaces for a subject
// (.emc_shear_support): per row of its design X, the normal u over the p
// sampled parameters (X's columns placed at cols, NA for constants) and the
// offsets of the finite limits in lim, shifted by the constant columns'
// contribution (their values in consts, in column order). NULL for a design
// that is not a matrix of doubles without missing values, which the R form
// handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_shear_half_spaces_cpp(SEXP X, SEXP cols, SEXP consts, SEXP lim, int p) {
  if (TYPEOF(X) != REALSXP || !Rf_isMatrix(X) || TYPEOF(cols) != INTSXP ||
      TYPEOF(consts) != REALSXP || TYPEOF(lim) != REALSXP || XLENGTH(lim) != 2 || p < 1)
    return R_NilValue;
  const int nr = Rf_nrows(X), nc = Rf_ncols(X);
  if (XLENGTH(cols) != nc) return R_NilValue;
  R_xlen_t n_const = 0;
  for (int j = 0; j < nc; ++j) {
    const int c = INTEGER(cols)[j];
    if (c == NA_INTEGER) ++n_const;
    else if (c < 1 || c > p) return R_NilValue;
  }
  if (XLENGTH(consts) != n_const) return R_NilValue;
  const double* x = REAL(X);
  for (R_xlen_t e = 0; e < XLENGTH(X); ++e)
    if (ISNAN(x[e])) return R_NilValue;
  for (R_xlen_t e = 0; e < n_const; ++e)
    if (ISNAN(REAL(consts)[e])) return R_NilValue;
  const double lo = REAL(lim)[0], hi = REAL(lim)[1];
  if (ISNAN(lo) || ISNAN(hi)) return R_NilValue;
  Rcpp::Shield<SEXP> keep(Rf_allocVector(VECSXP, 2 * static_cast<R_xlen_t>(nr)));
  R_xlen_t n_out = 0;
  auto half_space = [&](const std::vector<double>& u, double sign, double b) {
    Rcpp::NumericVector uv(p);
    for (int i = 0; i < p; ++i) uv[i] = sign < 0 ? -u[i] : u[i];
    SET_VECTOR_ELT(keep, n_out++,
                   Rcpp::List::create(Rcpp::Named("u") = uv, Rcpp::Named("b") = b));
  };
  std::vector<double> u(p);
  for (int r = 0; r < nr; ++r) {
    std::fill(u.begin(), u.end(), 0.0);
    long double total = 0.0;
    R_xlen_t k = 0;
    bool any = false;
    for (int j = 0; j < nc; ++j) {
      const double v = x[r + static_cast<R_xlen_t>(j) * nr];
      const int c = INTEGER(cols)[j];
      if (c == NA_INTEGER) {
        volatile double term = v * REAL(consts)[k++];
        total += term;
      } else {
        u[c - 1] = v;
      }
    }
    for (int i = 0; i < p; ++i) if (u[i] != 0) any = true;
    if (!any) continue;
    const double shift = n_const ? static_cast<double>(total) : 0.0;
    if (ISNAN(shift)) return R_NilValue;
    if (R_FINITE(lo)) {
      volatile double b = lo - shift;
      half_space(u, 1.0, b);
    }
    if (R_FINITE(hi)) {
      volatile double b = hi - shift;
      half_space(u, -1.0, -b);
    }
  }
  Rcpp::List out(n_out);
  for (R_xlen_t i = 0; i < n_out; ++i) out[i] = VECTOR_ELT(keep, i);
  return out;
}

// .emc_lpt_partition()'s assignment: subjects taken in `order` (numbered from
// one), each to the worker with the least load so far, the first of equals.
// A list of each worker's subjects. NULL for costs with missing values or an
// order that is not a permutation's worth of positions, which the R form
// handles.
// [[Rcpp::export(rng = false)]]
SEXP emc_lpt_partition_cpp(SEXP cost, SEXP order, int n_workers) {
  if (TYPEOF(cost) != REALSXP || TYPEOF(order) != INTSXP || n_workers < 2 ||
      XLENGTH(order) != XLENGTH(cost))
    return R_NilValue;
  const R_xlen_t n = XLENGTH(cost);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int s = INTEGER(order)[i];
    if (ISNAN(REAL(cost)[i]) || s == NA_INTEGER || s < 1 || s > n) return R_NilValue;
  }
  std::vector<std::vector<int> > mine(n_workers);
  std::vector<double> load(n_workers, 0.0);
  for (R_xlen_t i = 0; i < n; ++i) {
    const int s = INTEGER(order)[i];
    int w = 0;
    for (int k = 1; k < n_workers; ++k)
      if (load[k] < load[w]) w = k;
    mine[w].push_back(s);
    volatile double next = load[w] + REAL(cost)[s - 1];
    load[w] = next;
    if (ISNAN(load[w])) return R_NilValue;
  }
  Rcpp::List out(n_workers);
  for (int k = 0; k < n_workers; ++k) out[k] = Rcpp::wrap(mine[k]);
  return out;
}
