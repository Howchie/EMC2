// Per-iteration evidence records of the move allocation and the subject HMC
// (R/move_allocation.R, R/subject_hmc.R). Each function is the R loop it
// replaces, value for value: sums accumulate as R's, products are rounded
// before they are added, and results carry the attributes R's arithmetic
// would give them. A function returns NULL when an input is not of the
// shape it handles; the caller then runs the R loop.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <cstring>

using namespace Rcpp;

namespace {

bool finite_bits(double x) {
  uint64_t b;
  std::memcpy(&b, &x, sizeof(b));
  return (b & 0x7FF0000000000000ull) != 0x7FF0000000000000ull;
}

bool is_matrix(SEXP m, int p) {
  if (TYPEOF(m) != REALSXP) return false;
  SEXP dim = Rf_getAttrib(m, R_DimSymbol);
  return TYPEOF(dim) == INTSXP && XLENGTH(dim) == 2 && INTEGER(dim)[0] == p &&
    INTEGER(dim)[1] == p;
}

// m + scale * tcrossprod(d) as R computes it: every product rounded on its
// own, and the row names of `d` as dimnames when `m` has none.
SEXP add_outer(SEXP m, const double* d, int p, SEXP names, bool scaled, double scale) {
  SEXP out = PROTECT(Rf_duplicate(m));
  double* o = REAL(out);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < p; ++i) {
      volatile double term = d[i] * d[j];
      if (scaled) term = scale * term;
      o[i + j * p] = o[i + j * p] + term;
    }
  }
  if (names != R_NilValue && Rf_getAttrib(out, R_DimNamesSymbol) == R_NilValue) {
    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dn, 0, names);
    SET_VECTOR_ELT(dn, 1, names);
    Rf_setAttrib(out, R_DimNamesSymbol, dn);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

SEXP row_names(SEXP x) {
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  return TYPEOF(dn) == VECSXP && XLENGTH(dn) >= 1 ? VECTOR_ELT(dn, 0) : R_NilValue;
}

int field(SEXP lst, const char* name) {
  SEXP nm = Rf_getAttrib(lst, R_NamesSymbol);
  if (TYPEOF(nm) != STRSXP) return -1;
  for (R_xlen_t i = 0; i < XLENGTH(nm); ++i)
    if (std::strcmp(CHAR(STRING_ELT(nm, i)), name) == 0) return static_cast<int>(i);
  return -1;
}

SEXP get(SEXP lst, const char* name) {
  const int i = field(lst, name);
  return i < 0 ? R_NilValue : VECTOR_ELT(lst, i);
}

// A bare numeric scalar (no attributes), as a double.
bool plain_number(SEXP x, double& value, bool& integer) {
  if (ANY_ATTRIB(x) || XLENGTH(x) != 1) return false;
  if (TYPEOF(x) == REALSXP) { value = REAL(x)[0]; integer = false; return true; }
  if (TYPEOF(x) == INTSXP && INTEGER(x)[0] != NA_INTEGER) {
    value = INTEGER(x)[0]; integer = true; return true;
  }
  return false;
}

// cur + inc under R's arithmetic: integer only when both are.
SEXP plus(SEXP cur, double inc, bool inc_integer) {
  if (TYPEOF(cur) == INTSXP && inc_integer)
    return Rf_ScalarInteger(INTEGER(cur)[0] + static_cast<int>(inc));
  const double base = TYPEOF(cur) == INTSXP ? static_cast<double>(INTEGER(cur)[0]) : REAL(cur)[0];
  return Rf_ScalarReal(base + inc);
}

bool plain_scalar(SEXP x) {
  double v; bool i;
  return plain_number(x, v, i);
}

SEXP append_real(SEXP v, double x) {
  const R_xlen_t n = XLENGTH(v);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n + 1));
  if (n) std::memcpy(REAL(out), REAL(v), sizeof(double) * n);
  REAL(out)[n] = x;
  UNPROTECT(1);
  return out;
}

struct Task {
  bool present = false, has_disp = false, accepted = false, abort = false;
  bool rows_integer = false;
  double rows = 0.0, cpu = 0.0, acc_prob = 0.0, U = 0.0, score = 0.0;
  SEXP disp = R_NilValue, eps_next = R_NilValue;
};

// The fields .emc_hmc_record_tasks reads; false when one has another shape.
bool read_task(SEXP tk, int p, Task& t) {
  if (tk == R_NilValue) return true;
  if (TYPEOF(tk) != VECSXP) return false;
  t.present = true;
  bool integer = false;
  if (!plain_number(get(tk, "rows"), t.rows, t.rows_integer)) return false;
  if (!plain_number(get(tk, "cpu"), t.cpu, integer) || integer) return false;
  t.disp = get(tk, "disp");
  if (t.disp == R_NilValue) return true;
  t.has_disp = true;
  if (TYPEOF(t.disp) != REALSXP || XLENGTH(t.disp) != p ||
      Rf_getAttrib(t.disp, R_DimSymbol) != R_NilValue) return false;
  SEXP acc = get(tk, "accepted");
  t.accepted = TYPEOF(acc) == LGLSXP && XLENGTH(acc) == 1 && LOGICAL(acc)[0] == TRUE;
  SEXP status = get(tk, "status");
  t.abort = !(TYPEOF(status) == STRSXP && XLENGTH(status) == 1 &&
              !ANY_ATTRIB(status) &&
              std::strcmp(CHAR(STRING_ELT(status, 0)), "ok") == 0);
  if (!plain_number(get(tk, "acc_prob"), t.acc_prob, integer) || integer) return false;
  if (!plain_number(get(tk, "U"), t.U, integer) || integer) return false;
  if (!plain_number(get(tk, "score"), t.score, integer) || integer) return false;
  t.eps_next = get(tk, "eps_next");
  if (t.eps_next != R_NilValue && !(TYPEOF(t.eps_next) == REALSXP && plain_scalar(t.eps_next)))
    return false;
  return true;
}

const char* const kEvidenceNames[] = {"key", "n", "M", "accepted", "aborts", "rows",
                                      "cpu", "U", "score", "arows", "accp"};

SEXP new_evidence(int p, double eps) {
  SEXP ev = PROTECT(Rf_allocVector(VECSXP, 11));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 11));
  for (int i = 0; i < 11; ++i) SET_STRING_ELT(nm, i, Rf_mkChar(kEvidenceNames[i]));
  Rf_setAttrib(ev, R_NamesSymbol, nm);
  SET_VECTOR_ELT(ev, 0, Rf_ScalarReal(Rf_fprec(eps, 12.0)));
  SET_VECTOR_ELT(ev, 1, Rf_ScalarInteger(0));
  SEXP M = PROTECT(Rf_allocMatrix(REALSXP, p, p));
  std::memset(REAL(M), 0, sizeof(double) * p * p);
  SET_VECTOR_ELT(ev, 2, M);
  SET_VECTOR_ELT(ev, 3, Rf_ScalarInteger(0));
  SET_VECTOR_ELT(ev, 4, Rf_ScalarInteger(0));
  SET_VECTOR_ELT(ev, 5, Rf_ScalarReal(0.0));
  SET_VECTOR_ELT(ev, 6, Rf_ScalarReal(0.0));
  for (int i = 7; i < 11; ++i) SET_VECTOR_ELT(ev, i, Rf_allocVector(REALSXP, 0));
  UNPROTECT(3);
  return ev;
}

bool evidence_ok(SEXP ev, int p, int at[11]) {
  if (TYPEOF(ev) != VECSXP) return false;
  for (int i = 0; i < 11; ++i) {
    at[i] = field(ev, kEvidenceNames[i]);
    if (at[i] < 0) return false;
  }
  if (!is_matrix(VECTOR_ELT(ev, at[2]), p)) return false;
  for (int i : {1, 3, 4, 5, 6}) if (!plain_scalar(VECTOR_ELT(ev, at[i]))) return false;
  for (int i = 7; i < 11; ++i) {
    SEXP v = VECTOR_ELT(ev, at[i]);
    if (TYPEOF(v) != REALSXP || ANY_ATTRIB(v)) return false;
  }
  return true;
}

}  // namespace

// .emc_alloc_record_particle's loop: for each subject that ran and whose
// jump is finite, P + tcrossprod(jump), a count, and its cost (its worker
// time, or its share of `cost` when no worker time was measured).
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_record_particle_cpp(List P, IntegerVector nP, NumericVector cpu_part,
                                   NumericMatrix before, NumericMatrix after,
                                   LogicalVector ran, double cost, SEXP times) {
  const int p = before.nrow(), n = before.ncol();
  if (after.nrow() != p || after.ncol() != n || ran.size() != n || P.size() != n ||
      nP.size() != n || cpu_part.size() != n) return R_NilValue;
  if (times != R_NilValue && (TYPEOF(times) != REALSXP || XLENGTH(times) != n))
    return R_NilValue;
  std::vector<int> idx;
  for (int s = 0; s < n; ++s) {
    if (ran[s] == NA_LOGICAL) return R_NilValue;
    if (ran[s]) idx.push_back(s);
  }
  if (idx.empty()) return R_NilValue;
  const int m = static_cast<int>(idx.size());
  std::vector<double> t_ran(times == R_NilValue ? 0 : m);
  long double total = 0.0L;
  for (size_t k = 0; k < t_ran.size(); ++k) {
    const double x = REAL(times)[idx[k]];
    if (ISNAN(x)) return R_NilValue;
    t_ran[k] = x < 0.0 ? 0.0 : x;
    total += t_ran[k];
  }
  const double sum_t = static_cast<double>(total);
  const bool timed = sum_t > 0.0;
  for (int s : idx) if (!is_matrix(P[s], p)) return R_NilValue;
  List P_out = clone(P);
  IntegerVector nP_out = clone(nP);
  NumericVector cpu_out = clone(cpu_part);
  SEXP names = row_names(after);
  if (names == R_NilValue) names = row_names(before);
  std::vector<double> d(p);
  for (int k = 0; k < m; ++k) {
    const int s = idx[k];
    bool ok = true;
    for (int i = 0; i < p; ++i) {
      d[i] = after(i, s) - before(i, s);
      ok = ok && finite_bits(d[i]);
    }
    if (!ok) continue;
    SET_VECTOR_ELT(P_out, s, add_outer(P_out[s], d.data(), p, names, false, 1.0));
    nP_out[s] = nP_out[s] + 1;
    double add;
    if (timed) {
      add = t_ran[k];
    } else {
      volatile double share = 1.0 / m;
      volatile double part = cost * share;
      add = part;
    }
    cpu_out[s] = cpu_out[s] + add;
  }
  return List::create(Named("P") = P_out, Named("nP") = nP_out, Named("cpu_part") = cpu_out);
}

// .emc_alloc_record_subjects's loop: A + tcrossprod(jump) and
// D + tcrossprod(jump - shift) for every subject, or NULL when a jump is
// not finite. `shift` is NULL or one value per parameter.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_record_subjects_cpp(List A, List D, NumericMatrix alpha,
                                   NumericVector old, SEXP shift) {
  const int p = alpha.nrow(), n = alpha.ncol();
  if (old.size() != static_cast<R_xlen_t>(p) * n || A.size() != n || D.size() != n)
    return R_NilValue;
  const bool shifted = shift != R_NilValue;
  if (shifted && (TYPEOF(shift) != REALSXP || XLENGTH(shift) != p)) return R_NilValue;
  for (int s = 0; s < n; ++s)
    if (!is_matrix(A[s], p) || !is_matrix(D[s], p)) return R_NilValue;
  std::vector<double> d(static_cast<size_t>(p) * n);
  for (size_t i = 0; i < d.size(); ++i) {
    d[i] = alpha[i] - old[i];
    if (!finite_bits(d[i])) return List::create(Named("finite") = false);
  }
  SEXP names = row_names(alpha);
  List A_out(Rf_shallow_duplicate(A)), D_out(Rf_shallow_duplicate(D));
  std::vector<double> dd(p);
  for (int s = 0; s < n; ++s) {
    const double* ds = d.data() + static_cast<size_t>(s) * p;
    A_out[s] = add_outer(A[s], ds, p, names, false, 1.0);
    for (int i = 0; i < p; ++i) dd[i] = shifted ? ds[i] - REAL(shift)[i] : ds[i];
    D_out[s] = add_outer(D[s], dd.data(), p, names, false, 1.0);
  }
  return List::create(Named("finite") = true, Named("A") = A_out, Named("D") = D_out);
}

// The task loop of .emc_alloc_record_cost: each attempted subject's cost
// (its `work`, else `cost`, else `cpu`), attempt count and rows added to the
// window's vectors, with the cost total and number of attempts.
// [[Rcpp::export(rng = false)]]
SEXP emc_alloc_record_tasks_cpp(NumericVector cpu_hmc, IntegerVector hmc_att,
                                NumericVector rows_hmc, List tasks) {
  const int n = tasks.size();
  if (cpu_hmc.size() != n || hmc_att.size() != n || rows_hmc.size() != n)
    return R_NilValue;
  NumericVector cpu_out = clone(cpu_hmc), rows_out = clone(rows_hmc);
  IntegerVector att_out = clone(hmc_att);
  double work = 0.0;
  int attempts = 0;
  for (int s = 0; s < n; ++s) {
    SEXP tk = tasks[s];
    if (tk == R_NilValue) continue;
    if (TYPEOF(tk) != VECSXP) return R_NilValue;
    SEXP c = get(tk, "work");
    if (c == R_NilValue) c = get(tk, "cost");
    if (c == R_NilValue) c = get(tk, "cpu");
    double c_s, rows = 0.0;
    bool integer = false;
    if (!plain_number(c, c_s, integer) || integer) return R_NilValue;
    SEXP r = get(tk, "rows");
    if (r != R_NilValue && !plain_number(r, rows, integer)) return R_NilValue;
    cpu_out[s] = cpu_out[s] + c_s;
    att_out[s] = att_out[s] + 1;
    rows_out[s] = rows_out[s] + rows;
    work = work + c_s;
    ++attempts;
  }
  return List::create(Named("cpu_hmc") = cpu_out, Named("hmc_att") = att_out,
                      Named("rows_hmc") = rows_out, Named("work") = work,
                      Named("attempts") = attempts);
}

// .emc_hmc_record_tasks: one round's transitions added to the totals and to
// each subject's evidence, a warm-up move's next step size stored with its
// evidence key. `st` is the chain's HMC state.
// [[Rcpp::export(rng = false)]]
SEXP emc_hmc_record_tasks_cpp(List st, List tasks) {
  const int n = tasks.size();
  const int i_tot = field(st, "totals"), i_ev = field(st, "evidence");
  const int i_sub = field(st, "subjects"), i_p = field(st, "p");
  if (i_tot < 0 || i_ev < 0 || i_sub < 0 || i_p < 0) return R_NilValue;
  double p_value; bool p_integer;
  if (!plain_number(st[i_p], p_value, p_integer)) return R_NilValue;
  const int p = static_cast<int>(p_value);
  SEXP totals0 = st[i_tot], evidence0 = st[i_ev], subjects0 = st[i_sub];
  if (TYPEOF(totals0) != VECSXP || TYPEOF(evidence0) != VECSXP ||
      TYPEOF(subjects0) != VECSXP || XLENGTH(evidence0) < n || XLENGTH(subjects0) < n)
    return R_NilValue;
  const char* const total_names[] = {"attempts", "rows", "cpu", "accepted", "aborts"};
  int t_at[5];
  for (int i = 0; i < 5; ++i) {
    t_at[i] = field(totals0, total_names[i]);
    if (t_at[i] < 0 || !plain_scalar(VECTOR_ELT(totals0, t_at[i]))) return R_NilValue;
  }
  std::vector<Task> tk(n);
  for (int s = 0; s < n; ++s) {
    if (!read_task(tasks[s], p, tk[s])) return R_NilValue;
    if (!tk[s].has_disp) continue;
    SEXP ev = VECTOR_ELT(evidence0, s);
    int at[11];
    if (ev != R_NilValue && !evidence_ok(ev, p, at)) return R_NilValue;
    SEXP sj = VECTOR_ELT(subjects0, s);
    if (TYPEOF(sj) != VECSXP || field(sj, "eps") < 0) return R_NilValue;
    if (ev == R_NilValue && !(TYPEOF(get(sj, "eps")) == REALSXP && plain_scalar(get(sj, "eps"))))
      return R_NilValue;
  }
  List out(Rf_shallow_duplicate(st));
  List totals(Rf_shallow_duplicate(totals0));
  List evidence(Rf_shallow_duplicate(evidence0));
  List subjects(Rf_shallow_duplicate(subjects0));
  for (int s = 0; s < n; ++s) {
    const Task& t = tk[s];
    if (!t.present) continue;
    totals[t_at[0]] = plus(totals[t_at[0]], 1.0, true);
    totals[t_at[1]] = plus(totals[t_at[1]], t.rows, t.rows_integer);
    totals[t_at[2]] = plus(totals[t_at[2]], t.cpu, false);
    if (!t.has_disp) continue;
    totals[t_at[3]] = plus(totals[t_at[3]], t.accepted ? 1.0 : 0.0, true);
    totals[t_at[4]] = plus(totals[t_at[4]], t.abort ? 1.0 : 0.0, true);
    SEXP ev0 = evidence[s];
    List ev(ev0 == R_NilValue
              ? new_evidence(p, REAL(get(VECTOR_ELT(subjects, s), "eps"))[0])
              : Rf_shallow_duplicate(ev0));
    int at[11];
    for (int i = 0; i < 11; ++i) at[i] = field(ev, kEvidenceNames[i]);
    ev[at[1]] = plus(ev[at[1]], 1.0, true);
    ev[at[2]] = add_outer(ev[at[2]], REAL(t.disp), p, Rf_getAttrib(t.disp, R_NamesSymbol),
                          true, t.acc_prob);
    ev[at[3]] = plus(ev[at[3]], t.accepted ? 1.0 : 0.0, true);
    ev[at[4]] = plus(ev[at[4]], t.abort ? 1.0 : 0.0, true);
    ev[at[5]] = plus(ev[at[5]], t.rows, t.rows_integer);
    ev[at[6]] = plus(ev[at[6]], t.cpu, false);
    ev[at[7]] = append_real(ev[at[7]], t.U);
    ev[at[8]] = append_real(ev[at[8]], t.score);
    ev[at[9]] = append_real(ev[at[9]], t.rows);
    ev[at[10]] = append_real(ev[at[10]], t.acc_prob);
    if (t.eps_next != R_NilValue) {
      List sj(Rf_shallow_duplicate(VECTOR_ELT(subjects, s)));
      sj[field(sj, "eps")] = t.eps_next;
      subjects[s] = sj;
      ev[at[0]] = Rf_ScalarReal(Rf_fprec(REAL(t.eps_next)[0], 12.0));
    }
    evidence[s] = ev;
  }
  out[i_tot] = totals;
  out[i_ev] = evidence;
  out[i_sub] = subjects;
  return out;
}
