// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
                         Rcpp::Nullable<Rcpp::List> marginalise,
                         Rcpp::Nullable<Rcpp::LogicalVector> varying,
                         SEXP pt_mapper_context);

namespace {

// The prepared likelihood arguments of one (data, model) pair. Both objects
// are preserved and marked immutable, so their addresses identify them for
// as long as the entry lives.
struct LlEntry {
  SEXP dadm, model, args;
  std::vector<std::string> cols;
};

struct LlKey {
  SEXP dadm, model;
  bool operator==(const LlKey& o) const { return dadm == o.dadm && model == o.model; }
};

struct LlKeyHash {
  size_t operator()(const LlKey& k) const {
    const size_t a = std::hash<void*>()(static_cast<void*>(k.dadm));
    const size_t b = std::hash<void*>()(static_cast<void*>(k.model));
    return a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2));
  }
};

std::unordered_map<LlKey, LlEntry, LlKeyHash> ll_entries;
pid_t ll_entries_pid = 0;
const size_t kMaxEntries = 4096;

void release(LlEntry& e) {
  R_ReleaseObject(e.dadm);
  R_ReleaseObject(e.model);
  R_ReleaseObject(e.args);
}

void clear_entries() {
  for (auto& kv : ll_entries) release(kv.second);
  ll_entries.clear();
}

bool same_names(SEXP names, const std::vector<std::string>& cols) {
  if (TYPEOF(names) != STRSXP ||
      static_cast<size_t>(XLENGTH(names)) != cols.size()) return false;
  for (size_t j = 0; j < cols.size(); ++j)
    if (cols[j] != CHAR(STRING_ELT(names, j))) return false;
  return true;
}

SEXP column_names(SEXP x) {
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  if (TYPEOF(dn) != VECSXP || XLENGTH(dn) < 2) return R_NilValue;
  return VECTOR_ELT(dn, 1);
}

}  // namespace

// The log-likelihoods of `proposals` when this (data, model) pair has its
// arguments stored for these parameter columns; NULL otherwise (the caller
// then takes the full path and stores them). The likelihood draws no random
// numbers, so the generator state is not synchronised around the call.
// [[Rcpp::export(rng = false)]]
SEXP emc2_ll_cached(SEXP proposals, SEXP dadm, SEXP model, SEXP varying) {
  // Parameter-table contexts are not carried across a fork.
  const pid_t pid = getpid();
  if (pid != ll_entries_pid) {
    clear_entries();
    ll_entries_pid = pid;
    return R_NilValue;
  }
  auto it = ll_entries.find(LlKey{dadm, model});
  if (it == ll_entries.end()) return R_NilValue;
  const LlEntry& e = it->second;
  if (TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) ||
      !same_names(column_names(proposals), e.cols)) return R_NilValue;
  const List a(e.args);
  const SEXP trend = a["trend"], context = a["pt_context"];
  return calc_ll_oo(NumericMatrix(proposals), DataFrame(dadm),
                    as<NumericVector>(a["constants"]), as<List>(a["designs"]),
                    as<String>(a["type"]), as<List>(a["bounds"]),
                    as<List>(a["transforms"]), as<List>(a["pretransforms"]),
                    as<CharacterVector>(a["p_types"]), as<double>(a["min_ll"]),
                    Nullable<List>(trend), R_NilValue,
                    Nullable<LogicalVector>(varying), context);
}

// Whether emc2_ll_cached() would evaluate these columns for this pair.
bool emc2_ll_cache_has(SEXP dadm, SEXP model, SEXP cols) {
  if (getpid() != ll_entries_pid) return false;
  auto it = ll_entries.find(LlKey{dadm, model});
  return it != ll_entries.end() && same_names(cols, it->second.cols);
}

// emc2_ll_cached() for one matrix of rows per subject. An element stays
// NULL when the subject has no stored arguments, the call fails, or a value is
// NA or +Inf; the caller evaluates those subjects on the full path, which
// reports the failure.
// [[Rcpp::export(rng = false)]]
List emc2_ll_cached_subjects(List rows, List data, SEXP model) {
  const R_xlen_t n = rows.size();
  List out(n);
  if (data.size() != n) return out;
  for (R_xlen_t k = 0; k < n; ++k) {
    const SEXP r = rows[k];
    if (TYPEOF(r) != REALSXP || !Rf_isMatrix(r)) continue;
    try {
      const Shield<SEXP> ll(emc2_ll_cached(r, data[k], model, R_NilValue));
      if (TYPEOF(ll) != REALSXP || XLENGTH(ll) != Rf_nrows(r) ||
          ANY_ATTRIB(ll)) continue;
      const double* v = REAL(ll);
      bool valid = true;
      for (R_xlen_t i = 0; i < XLENGTH(ll); ++i)
        if (ISNAN(v[i]) || v[i] == R_PosInf) { valid = false; break; }
      if (valid) out[k] = ll;
    } catch (Rcpp::LongjumpException&) {
      throw;
    } catch (Rcpp::internal::InterruptedException&) {
      throw;
    } catch (...) {
    }
  }
  return out;
}

// [[Rcpp::export]]
void emc2_ll_cache_put(SEXP dadm, SEXP model, CharacterVector cols, List args) {
  const pid_t pid = getpid();
  if (pid != ll_entries_pid) {
    clear_entries();
    ll_entries_pid = pid;
  }
  if (ll_entries.size() >= kMaxEntries) clear_entries();
  const LlKey key{dadm, model};
  auto it = ll_entries.find(key);
  if (it != ll_entries.end()) {
    release(it->second);
    ll_entries.erase(it);
  }
  LlEntry e;
  e.dadm = dadm;
  e.model = model;
  e.args = args;
  e.cols = as<std::vector<std::string>>(cols);
  R_PreserveObject(e.dadm);
  R_PreserveObject(e.model);
  R_PreserveObject(e.args);
  MARK_NOT_MUTABLE(e.dadm);
  MARK_NOT_MUTABLE(e.model);
  ll_entries.emplace(key, e);
}

// [[Rcpp::export]]
void emc2_ll_cache_clear() {
  clear_entries();
}
