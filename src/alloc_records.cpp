// Per-iteration evidence record of the subject HMC (R/subject_hmc.R): the R
// loop it replaces, value for value: sums accumulate as R's, products are rounded
// before they are added, and results carry the attributes R's arithmetic
// would give them. It returns NULL when an input is not of the shape it
// handles; the caller then runs the R loop.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <cstring>

using namespace Rcpp;

namespace {

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
  bool present = false, moved = false, accepted = false, abort = false;
  bool rows_integer = false;
  double rows = 0.0, cpu = 0.0, acc_prob = 0.0, U = 0.0, score = 0.0;
  SEXP eps_next = R_NilValue;
};

// The fields .emc_hmc_record_tasks reads; false when one has another shape.
bool read_task(SEXP tk, int p, Task& t) {
  if (tk == R_NilValue) return true;
  if (TYPEOF(tk) != VECSXP) return false;
  t.present = true;
  bool integer = false;
  if (!plain_number(get(tk, "rows"), t.rows, t.rows_integer)) return false;
  if (!plain_number(get(tk, "cpu"), t.cpu, integer) || integer) return false;
  SEXP disp = get(tk, "disp");
  if (disp == R_NilValue) return true;
  t.moved = true;
  if (TYPEOF(disp) != REALSXP || XLENGTH(disp) != p ||
      Rf_getAttrib(disp, R_DimSymbol) != R_NilValue) return false;
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

// A subject's evidence at one step size: its key, the transitions, and per
// transition the jitter U, the acceptance-weighted squared jump (score), the
// likelihood rows and the acceptance probability.
constexpr int kEvidence = 6;
const char* const kEvidenceNames[kEvidence] = {"key", "n", "U", "score", "arows", "accp"};

SEXP new_evidence(double eps) {
  SEXP ev = PROTECT(Rf_allocVector(VECSXP, kEvidence));
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, kEvidence));
  for (int i = 0; i < kEvidence; ++i) SET_STRING_ELT(nm, i, Rf_mkChar(kEvidenceNames[i]));
  Rf_setAttrib(ev, R_NamesSymbol, nm);
  SET_VECTOR_ELT(ev, 0, Rf_ScalarReal(Rf_fprec(eps, 12.0)));
  SET_VECTOR_ELT(ev, 1, Rf_ScalarInteger(0));
  for (int i = 2; i < kEvidence; ++i) SET_VECTOR_ELT(ev, i, Rf_allocVector(REALSXP, 0));
  UNPROTECT(2);
  return ev;
}

bool evidence_ok(SEXP ev, int at[kEvidence]) {
  if (TYPEOF(ev) != VECSXP) return false;
  for (int i = 0; i < kEvidence; ++i) {
    at[i] = field(ev, kEvidenceNames[i]);
    if (at[i] < 0) return false;
  }
  if (!plain_scalar(VECTOR_ELT(ev, at[1]))) return false;
  for (int i = 2; i < kEvidence; ++i) {
    SEXP v = VECTOR_ELT(ev, at[i]);
    if (TYPEOF(v) != REALSXP || ANY_ATTRIB(v)) return false;
  }
  return true;
}

}  // namespace

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
    if (!tk[s].moved) continue;
    SEXP ev = VECTOR_ELT(evidence0, s);
    int at[kEvidence];
    if (ev != R_NilValue && !evidence_ok(ev, at)) return R_NilValue;
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
    if (!t.moved) continue;
    totals[t_at[3]] = plus(totals[t_at[3]], t.accepted ? 1.0 : 0.0, true);
    totals[t_at[4]] = plus(totals[t_at[4]], t.abort ? 1.0 : 0.0, true);
    SEXP ev0 = evidence[s];
    List ev(ev0 == R_NilValue
              ? new_evidence(REAL(get(VECTOR_ELT(subjects, s), "eps"))[0])
              : Rf_shallow_duplicate(ev0));
    int at[kEvidence];
    for (int i = 0; i < kEvidence; ++i) at[i] = field(ev, kEvidenceNames[i]);
    ev[at[1]] = plus(ev[at[1]], 1.0, true);
    ev[at[2]] = append_real(ev[at[2]], t.U);
    ev[at[3]] = append_real(ev[at[3]], t.score);
    ev[at[4]] = append_real(ev[at[4]], t.rows);
    ev[at[5]] = append_real(ev[at[5]], t.acc_prob);
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
