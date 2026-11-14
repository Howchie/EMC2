// The stored Laplace checks' bookkeeping (R/subject_fit_cache.R).
#include <Rcpp.h>
#include <cstring>

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
