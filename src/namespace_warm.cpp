// Which lazy-loaded bindings of an environment are still unforced.
//
// A package namespace binds each object to a promise that loads it from the
// lazy-load database on first use. A forked child that forces such a promise
// loads its own copy, which dies with the child, so a process that forks a
// child per block pays the same loads in every block. Reporting the bindings
// a child forced lets the parent force them once, and every later fork then
// inherits them.

#include <Rcpp.h>
#include <Rversion.h>

//' Lazy-load state of namespace bindings
//'
//' @param env Environment to inspect, normally the package namespace.
//' @param names Binding names to look up.
//' @return Logical vector, TRUE where the binding is a promise that has not
//'   been forced. Always FALSE before R 4.6, which has no public way to ask
//'   without forcing.
//' @noRd
// [[Rcpp::export]]
Rcpp::LogicalVector emc_ns_unforced(SEXP env, Rcpp::CharacterVector names) {
  const R_xlen_t n = names.size();
  Rcpp::LogicalVector out(n, false);
#if R_VERSION >= R_Version(4, 6, 0)
  if (TYPEOF(env) != ENVSXP) return out;
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP sym = Rf_installChar(STRING_ELT(names, i));
    out[i] = R_GetBindingType(sym, env) == R_BindingTypeDelayed;
  }
#else
  (void) env;
#endif
  return out;
}
