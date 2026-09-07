#pragma once

#include <Rcpp.h>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cstring>
#include "utility_types.h"
using Rcpp::_;

struct DesignEntry {
  bool valid;                    // does this design have a usable mapping?
  bool skip_self_intercept;      // from design_is_self_intercept
  int out_idx;                   // column index in base to write into
  std::vector<int> coef_idx;     // length K; -1 for unused
  bool uses_self;                // does any coef refer to out_idx?
  bool split_transform;          // transform(sum(pre terms)) + sum(post terms)
  TransformCode split_code;      // transform to apply to pre-sum
  double split_lower;
  double split_upper;
  std::vector<char> coef_in_pre_sum; // length K; whether this design column contributes to pre-sum
  // 0-based row map for a compressed design matrix; empty means identity.
  std::vector<int> expand_idx;
  // Every trial reads the SAME design row, so this design's output column is
  // row-constant whenever all of its coefficient columns are. Computed once in
  // init_design_plan; see the row-constant notes on ParamTable::col_const.
  bool single_cell = false;
  // Design-CELL structure: this design reads n_cells distinct rows, so its
  // output column takes at most n_cells distinct values laid out by
  // expand_idx.  cell_rep[c] is the first trial row that reads design row c
  // (-1 when no trial does).  single_cell is the n_cells == 1 case.
  bool cell_usable = false;
  int n_cells = 0;
  std::vector<int> cell_rep;
};

// Above this many design cells the gather stops paying for itself and the
// per-cell scratch would not fit on the stack; such designs take the plain
// per-row path.
static const int EMC2_PT_MAX_CELLS = 256;

// View-based ParamTable: one base matrix + active column indices
// Invariants:
// - base contains all parameters (never shrinks)
// - active_cols is a view into base (ordering may change)
// - name_to_base_idx maps ALL base columns, active or not
struct ParamTable {
  Rcpp::NumericMatrix base;     // underlying matrix (all columns)

  std::vector<DesignEntry> design_plan;   // cached once

  Rcpp::CharacterVector base_names; // colnames for base
  int n_trials = 0;

  std::vector<int> active_cols;  // indices into base/base_names
  std::vector<char> is_active;  // size = base.ncol()

  // fast look-up map
  std::unordered_map<std::string,int> name_to_base_idx;

  // ---------------------------------------------------------------------
  // Row-constant column tracking.
  //
  // A parameter column carries only DESIGN-CELL information, but the whole
  // per-particle prologue (map -> transform -> bounds) runs at TRIAL
  // resolution. For the common intercept-only parameter that means n_trials
  // scalar std::exp calls and n_trials bound comparisons per particle to
  // produce one distinct number. col_const[j] records "base column j currently
  // holds one repeated value", which lets those three stages collapse to O(1)
  // plus a std::fill.
  //
  // The flags are MAINTAINED, never sniffed: reset_base_to_zero marks
  // everything constant, the fill_* helpers mark each column they scalar-fill,
  // and map_from_designs recomputes the flag for every output it writes.
  // Columns touched by none of those (the invariant-parameter lane) correctly
  // keep the flag they were given on the template particle, because their
  // contents are likewise carried over.
  //
  // Tracking is opt-in (const_init) and is enabled only when nothing outside
  // this pipeline writes into `base` -- i.e. when there is no TrendRuntime,
  // whose apply_base_for_op does exactly that.
  std::vector<char> col_const;
  // Per base column: the design_plan index whose expand map lays this column's
  // distinct values out across trials, or -1 for no usable cell structure.
  // This is the general case of which col_const is the n_cells == 1 corner:
  // a column with n_cells distinct values needs n_cells transcendentals and
  // n_cells bound comparisons, not n_trials of each.
  std::vector<int> col_cell;
  bool track_const = false;

  void const_init() {
    col_const.assign(base.ncol(), 1);
    col_cell.assign(base.ncol(), -1);
    track_const = true;
  }
  inline bool col_is_const(int j) const {
    const bool flagged =
      track_const && j >= 0 && j < (int)col_const.size() && col_const[j];
#ifdef EMC2_PT_CONST_CHECK
    if (flagged) const_check_column(j);
#endif
    return flagged;
  }
  inline void const_mark(int j, bool is_const) {
    if (track_const && j >= 0 && j < (int)col_const.size()) {
      col_const[j] = is_const ? 1 : 0;
      col_cell[j] = -1;
    }
  }

  // The design entry whose cell map governs column j, or nullptr.  A
  // row-constant column answers nullptr: its own O(1) path is cheaper.
  inline const DesignEntry* col_cell_entry(int j) const {
    if (!track_const || j < 0 || j >= (int)col_cell.size()) return nullptr;
    if (col_const[j] || col_cell[j] < 0) return nullptr;
    const DesignEntry& e = design_plan[col_cell[j]];
#ifdef EMC2_PT_CONST_CHECK
    cell_check_column(j, e);
#endif
    return &e;
  }

#ifdef EMC2_PT_CONST_CHECK
  // Assertion build (-DEMC2_PT_CONST_CHECK): every time a row-constant claim is
  // about to be ACTED on, re-scan the column and stop on disagreement. Checking
  // at the point of use rather than at the end of the prologue matters: a wrong
  // flag makes the very next transform fill the column from row 0, after which
  // a whole-table scan would find it consistent and pass.
  void cell_check_column(int j, const DesignEntry& e) const {
    const double* col = &base(0, j);
    for (int r = 0; r < n_trials; ++r) {
      const int rep = e.cell_rep[e.expand_idx[r]];
      const bool same = (col[r] == col[rep]) ||
                        (Rcpp::NumericVector::is_na(col[r]) &&
                         Rcpp::NumericVector::is_na(col[rep]));
      if (!same) {
        Rcpp::stop("ParamTable cell check: column '%s' is flagged cell-constant "
                   "but row %d holds %g against its cell representative row %d's %g",
                   Rcpp::as<std::string>(base_names[j]).c_str(), r, col[r], rep,
                   col[rep]);
      }
    }
  }

  void const_check_column(int j) const {
    const double* col = &base(0, j);
    const double v0 = col[0];
    for (int r = 1; r < n_trials; ++r) {
      const bool same = (col[r] == v0) ||
                        (Rcpp::NumericVector::is_na(col[r]) &&
                         Rcpp::NumericVector::is_na(v0));
      if (!same) {
        Rcpp::stop("ParamTable const check: column '%s' is flagged row-constant "
                   "but row %d holds %g against row 0's %g",
                   Rcpp::as<std::string>(base_names[j]).c_str(), r, col[r], v0);
      }
    }
  }
#endif

  // keep track of parameters that have an intercept-only design
  std::vector<char> design_is_self_intercept;
  std::unordered_map<std::string, std::unordered_set<std::string>> pre_sum_term_map_;
  std::unordered_map<std::string, TransformCode> split_code_map_;
  std::unordered_map<std::string, std::pair<double, double>> split_bounds_map_;

  ParamTable() = default;

  ParamTable(Rcpp::NumericMatrix base_,
             Rcpp::CharacterVector names_)
    : base(base_), base_names(names_) {
    n_trials = base_.nrow();
    const int p = base_.ncol();
    active_cols.resize(p);
    std::iota(active_cols.begin(), active_cols.end(), 0);
    Rcpp::colnames(base) = base_names;
    rebuild_name_map();
  }

  void rebuild_name_map() {
    name_to_base_idx.clear();
    const int p = base_names.size();
    name_to_base_idx.reserve(p);
    for (int j = 0; j < p; ++j) {
      name_to_base_idx[ Rcpp::as<std::string>(base_names[j]) ] = j;
    }
  }

  int base_index_for(const std::string& nm) const {
    auto it = name_to_base_idx.find(nm);
    if (it == name_to_base_idx.end()) {
      Rcpp::stop("ParamTable: Unknown parameter '%s'", nm.c_str());
    }
    return it->second;
  }
  int n_params() const { return (int)active_cols.size(); }

  // Name of j-th active column
  Rcpp::String name_at(int j) const {
    return base_names[ active_cols[j] ];
  }

  // Column view by active index (no deep copy)
  Rcpp::NumericVector column_by_active_index(int j) const {
    return base(_, active_cols[j]);  // view
  }

  bool is_active_base_idx(int base_idx) const {
    return std::find(active_cols.begin(), active_cols.end(), base_idx) != active_cols.end();
  }

  // Column view by parameter name
  Rcpp::NumericVector column_by_name(const std::string& nm) const {
    int base_idx = base_index_for(nm); // O(1) from map
        return base(_, base_idx);  // view
  }

  // Assign into existing column (writes into base)
  void set_column_by_name(const std::string& nm,
                          const Rcpp::NumericVector& col) {
    int base_idx = base_index_for(nm);
        // Arbitrary external write: conservatively drop this column's
        // row-constant claim rather than inspect the incoming vector.
        const_mark(base_idx, false);
        base(_, base_idx) = col;
        return;
  }


  void set_transform_metadata(const Rcpp::List& transform) {
    pre_sum_term_map_.clear();
    split_code_map_.clear();
    split_bounds_map_.clear();

    if (!transform.containsElementNamed("pre_sum_terms")) return;
    Rcpp::List pre_sum_terms = transform["pre_sum_terms"];
    Rcpp::CharacterVector pre_sum_names = pre_sum_terms.names();
    for (int i = 0; i < pre_sum_terms.size(); ++i) {
      std::string param = Rcpp::as<std::string>(pre_sum_names[i]);
      Rcpp::CharacterVector term_names = pre_sum_terms[i];
      auto& cur_set = pre_sum_term_map_[param];
      for (int j = 0; j < term_names.size(); ++j) {
        cur_set.insert(Rcpp::as<std::string>(term_names[j]));
      }
    }

    if (transform.containsElementNamed("func")) {
      Rcpp::CharacterVector func = transform["func"];
      Rcpp::CharacterVector fnames = func.names();
      for (int i = 0; i < func.size(); ++i) {
        std::string name = Rcpp::as<std::string>(fnames[i]);
        std::string f = Rcpp::as<std::string>(func[i]);
        if (f == "exp") split_code_map_[name] = EXP;
        else if (f == "pnorm") split_code_map_[name] = PNORM;
        else split_code_map_[name] = IDENTITY;
      }
    }
    if (transform.containsElementNamed("lower")) {
      Rcpp::NumericVector lower = transform["lower"];
      Rcpp::CharacterVector lnames = lower.names();
      for (int i = 0; i < lower.size(); ++i) {
        split_bounds_map_[Rcpp::as<std::string>(lnames[i])].first = lower[i];
      }
    }
    if (transform.containsElementNamed("upper")) {
      Rcpp::NumericVector upper = transform["upper"];
      Rcpp::CharacterVector unames = upper.names();
      for (int i = 0; i < upper.size(); ++i) {
        split_bounds_map_[Rcpp::as<std::string>(unames[i])].second = upper[i];
      }
    }
  }

  void init_design_plan(const Rcpp::List& designs) {
    using namespace Rcpp;
    using std::string;

    const int n_params = designs.size();
    design_plan.clear();
    design_plan.resize(n_params);

    CharacterVector design_names = designs.names();
    if (design_names.size() != n_params) {
      stop("ParamTable::init_design_plan: designs must be a named list");
    }

    const int T = base.nrow();

    for (int i = 0; i < n_params; ++i) {
      DesignEntry& entry = design_plan[i];
      entry.valid = false;
      entry.split_transform = false;
      entry.split_code = IDENTITY;
      entry.split_lower = 0.0;
      entry.split_upper = 1.0;
      entry.expand_idx.clear();
      entry.single_cell = false;
      entry.cell_usable = false;
      entry.n_cells = 0;
      entry.cell_rep.clear();

      if (designs[i] == R_NilValue) continue;

      // self-intercept-only designs can be marked as skippable
      bool skip_self = !design_is_self_intercept.empty() &&
        i < (int)design_is_self_intercept.size() &&
        design_is_self_intercept[i];
      entry.skip_self_intercept = skip_self;

      const string out_name = as<string>(design_names[i]);
      auto out_it = name_to_base_idx.find(out_name);
      if (out_it == name_to_base_idx.end()) {
        // no matching output column; nothing to do
        continue;
      }
      entry.out_idx = out_it->second;

      NumericMatrix design = designs[i];
      IntegerVector expand;
      SEXP expand_attr = design.attr("expand");
      if (expand_attr != R_NilValue) expand = as<IntegerVector>(expand_attr);
      if (expand.size() == T) {
        entry.expand_idx.resize(T);
        for (int r = 0; r < T; ++r) {
          const int idx = expand[r];
          if (IntegerVector::is_na(idx) || idx < 1 || idx > design.nrow()) {
            stop("ParamTable::init_design_plan: invalid expand index for '%s'",
                 out_name.c_str());
          }
          entry.expand_idx[r] = idx - 1;
        }
      } else if (design.nrow() != T) {
        stop("ParamTable::init_design_plan: design for '%s' must have n_trials rows",
             out_name.c_str());
      }

      const int K = design.ncol();
      CharacterVector coef_names = colnames(design);

      entry.coef_idx.assign(K, -1);
      entry.coef_in_pre_sum.assign(K, 0);
      entry.uses_self = false;

      auto split_it = pre_sum_term_map_.find(out_name);
      if (split_it != pre_sum_term_map_.end()) {
        entry.split_transform = true;
        auto code_it = split_code_map_.find(out_name);
        if (code_it != split_code_map_.end()) entry.split_code = code_it->second;
        auto bound_it = split_bounds_map_.find(out_name);
        if (bound_it != split_bounds_map_.end()) {
          entry.split_lower = bound_it->second.first;
          entry.split_upper = bound_it->second.second;
        }
      }

      for (int j = 0; j < K; ++j) {
        string coef_name = as<string>(coef_names[j]);
        auto it = name_to_base_idx.find(coef_name);
        if (it == name_to_base_idx.end()) continue;

        int cidx = it->second;
        entry.coef_idx[j] = cidx;
        if (cidx == entry.out_idx) entry.uses_self = true;
        if (entry.split_transform && split_it->second.find(coef_name) != split_it->second.end()) {
          entry.coef_in_pre_sum[j] = 1;
        }
      }

      // Row-constant capability: does every trial read the same design row?
      // With a compressed design that is a constant expand vector; without one
      // it is a single-row design (or a single-trial table).
      if (!entry.expand_idx.empty()) {
        bool same = true;
        const int first = entry.expand_idx[0];
        for (int r = 1; r < T; ++r) {
          if (entry.expand_idx[r] != first) { same = false; break; }
        }
        entry.single_cell = same;
      } else {
        entry.single_cell = (design.nrow() == 1) || (T <= 1);
      }

      // Cell structure: worth exploiting only for a compressed design that
      // genuinely has fewer rows than trials, and only up to the scratch cap.
      if (!entry.single_cell && !entry.expand_idx.empty() &&
          design.nrow() < T && design.nrow() <= EMC2_PT_MAX_CELLS) {
        entry.n_cells = design.nrow();
        entry.cell_rep.assign(entry.n_cells, -1);
        for (int r = 0; r < T; ++r) {
          const int c = entry.expand_idx[r];
          if (entry.cell_rep[c] < 0) entry.cell_rep[c] = r;
        }
        entry.cell_usable = true;
      }

      entry.valid = true;
    }
  }

  static inline double apply_split_transform_scalar(double x, TransformCode code,
                                                    double lower, double upper) {
    switch (code) {
    case EXP:
      return lower + std::exp(x);
    case PNORM:
      return lower + (upper - lower) * R::pnorm(x, 0.0, 1.0, 1, 0);
    case IDENTITY:
    default:
      return x;
    }
  }

  Rcpp::NumericMatrix materialize_by_param_names(const Rcpp::CharacterVector& param_names) const {
    using namespace Rcpp;
    const int k = param_names.size();
    NumericMatrix out(n_trials, k);
    CharacterVector out_names(k);

    for (int j = 0; j < k; ++j) {
      std::string nm = as<std::string>(param_names[j]);
      int base_idx = base_index_for(nm);  // throws if unknown


      // Copy base(:, base_idx) → out(:, j)
      double* out_col        = &out(0, j);
      const double* base_col = &base(0, base_idx);
      for (int r = 0; r < n_trials; ++r) {
        out_col[r] = base_col[r];
      }
      out_names[j] = base_names[base_idx];
    }

    colnames(out) = out_names;
    return out;
  }

  // Refill a pre-allocated (n_trials x k) matrix with base columns in the
  // caller-resolved order (base_idx_order[j] = base column for out column j).
  // Per-particle twin of materialize_by_param_names: the caller allocates the
  // matrix and resolves names ONCE per likelihood call, so refilling per
  // particle is a plain memcpy with no R-heap allocation.
  // As materialize_into, but only for the listed output columns.  The mixed
  // (censored/truncated) race path re-materialises the whole parameter matrix
  // for every particle; the columns whose designs carry no sampled coefficient
  // hold the same natural-scale values for every particle, so after the first
  // fill they are pure memcpy waste.
  void materialize_into_subset(Rcpp::NumericMatrix& out,
                               const std::vector<int>& base_idx_order,
                               const std::vector<int>& which_out_cols) const {
    for (std::size_t k = 0; k < which_out_cols.size(); ++k) {
      const int j = which_out_cols[k];
      std::memcpy(&out(0, j), &base(0, base_idx_order[j]),
                  static_cast<size_t>(n_trials) * sizeof(double));
    }
  }

  void materialize_into(Rcpp::NumericMatrix& out,
                        const std::vector<int>& base_idx_order) const {
    const int k = (int)base_idx_order.size();
    for (int j = 0; j < k; ++j) {
      double* out_col = &out(0, j);
      const double* base_col = &base(0, base_idx_order[j]);
      std::memcpy(out_col, base_col, static_cast<size_t>(n_trials) * sizeof(double));
    }
  }

  // Materialize matrix (n_trials x n_active) with all columns
  Rcpp::NumericMatrix materialize() const {
    const int p = (int)active_cols.size();
    Rcpp::NumericMatrix out(n_trials, p);
    Rcpp::CharacterVector out_names(p);

    for (int j = 0; j < p; ++j) {
      int base_j = active_cols[j];
      // avoid Rcpp Sugar copies like out(_, j) = base(_, base_j);
      double* out_col = &out(0, j);
      const double* base_col = &base(0, base_j);
      for (int r = 0; r < n_trials; ++r) {
        out_col[r] = base_col[r];
      }
      out_names[j] = base_names[base_j];
    }
    Rcpp::colnames(out) = out_names;
    return out;
  }

  static ParamTable from_p_types(int n_trials,
                                 const Rcpp::CharacterVector& p_types) {
    const int p = p_types.size();
    Rcpp::NumericMatrix base(n_trials, p);
    base.fill(0.0);
    Rcpp::CharacterVector names = Rcpp::clone(p_types);
    return ParamTable(base, names);
  }

  static ParamTable from_p_vector_and_designs(const Rcpp::NumericVector& p_vector,
                                              const Rcpp::List& designs,
                                              int n_trials,
                                              const Rcpp::List& transforms = Rcpp::List::create()) {
    using std::string;
    using Rcpp::as;
    using Rcpp::CharacterVector;
    using Rcpp::NumericMatrix;
    using Rcpp::NumericVector;

    // 1) Collect all unique parameter names in a stable order:
    //    (a) names(p_vector)
    //    (b) names(designs)
    //    (c) colnames(designs[[i]])
    std::vector<string> names_vec;
    names_vec.reserve(p_vector.size() + designs.size() * 2);
    std::unordered_set<string> seen;

    // (a) names(p_vector)
    CharacterVector pv_names = p_vector.names();
    for (int i = 0; i < pv_names.size(); ++i) {
      string nm = as<string>(pv_names[i]);
      if (!seen.count(nm)) {
        seen.insert(nm);
        names_vec.push_back(nm);
      }
    }

    // (b) names(designs)
    CharacterVector design_names = designs.names();
    for (int i = 0; i < design_names.size(); ++i) {
      string nm = as<string>(design_names[i]);
      if (!seen.count(nm)) {
        seen.insert(nm);
        names_vec.push_back(nm);
      }
    }

    // (c) colnames of each design matrix
    for (int i = 0; i < designs.size(); ++i) {
      if (designs[i] == R_NilValue) continue;
      NumericMatrix dm = designs[i];
      CharacterVector cn = Rcpp::colnames(dm);
      for (int j = 0; j < cn.size(); ++j) {
        string nm = as<string>(cn[j]);
        if (!seen.count(nm)) {
          seen.insert(nm);
          names_vec.push_back(nm);
        }
      }
    }

    // 2) Build base matrix and names vector
    const int p = static_cast<int>(names_vec.size());
    NumericMatrix base(n_trials, p);
    base.fill(0.0);

    CharacterVector base_names(p);
    for (int j = 0; j < p; ++j) {
      base_names[j] = names_vec[j];
    }
    Rcpp::colnames(base) = base_names;

    // 3) Build map name -> value from p_vector
    std::unordered_map<string, double> pval;
    pval.reserve(p_vector.size());
    for (int i = 0; i < p_vector.size(); ++i) {
      string nm = as<string>(pv_names[i]);
      pval[nm] = p_vector[i];
    }

    // 4) Fill columns that appear in p_vector with that scalar value, others remain 0
    for (int j = 0; j < p; ++j) {
      string nm = names_vec[j];
      auto it = pval.find(nm);
      if (it != pval.end()) {
        double val = it->second;
        for (int r = 0; r < n_trials; ++r) {
          base(r, j) = val;
        }
      }
      // else: not in p_vector -> keep zeros
    }

    // 5) Construct ParamTable

    // 5) Construct ParamTable
    ParamTable pt(base, base_names);
    pt.set_transform_metadata(transforms);
    pt.n_trials = n_trials;  // (already set in ctor, but fine to be explicit)

    // cache "self-intercept-only" designs
    const int n_designs = designs.size();
    pt.design_is_self_intercept.assign(n_designs, 0);

    for (int i = 0; i < n_designs; ++i) {
      if (designs[i] == R_NilValue) continue;

      NumericMatrix dm = designs[i];

      Rcpp::IntegerVector expand;
      SEXP expand_attr = dm.attr("expand");
      if (expand_attr != R_NilValue) expand = Rcpp::as<Rcpp::IntegerVector>(expand_attr);
      if (dm.nrow() != n_trials && expand.size() != n_trials) continue;

      // Exactly 1 column
      if (dm.ncol() != 1) continue;

      CharacterVector cn = Rcpp::colnames(dm);
      if (cn.size() != 1) continue;

      // Column name must match the design/output name
      string out_name  = as<string>(design_names[i]);
      string coef_name = as<string>(cn[0]);
      if (coef_name != out_name) continue;

      // Column must be all ones
      bool all_ones = true;
      const bool compressed = expand.size() == n_trials;
      for (int r = 0; r < n_trials; ++r) {
        const int drow = compressed ? expand[r] - 1 : r;
        if (drow < 0 || drow >= dm.nrow() || dm(drow, 0) != 1.0) {
          all_ones = false;
          break;
        }
      }

      if (all_ones) {
        pt.design_is_self_intercept[i] = 1;
      }
    }

    return pt;
  }

  // Map parameters from designs into this table, ignoring trends.
  // Uses the *current* trial-by-trial values in ParamTable::base as
  // coefficients for design columns.
  //
  // designs: named list of design matrices (possibly compressed).
  //          Each matrix may have an "expand" attribute:
  //          - if present and length == n_trials: use it to expand rows
  //          - else, require nrow(design) == n_trials.
  // include_param: optional logical mask; if length == length(designs),
  //                decides which entries to map; if empty, all TRUE.
  void map_from_designs(const Rcpp::List& designs,
                        const Rcpp::LogicalVector& include_param = Rcpp::LogicalVector()) {
    using namespace Rcpp;
    using std::string;

    const int n_params = designs.size();
    if (n_params == 0) return;

    // If no include_param or wrong length, include all
    LogicalVector use =
      (include_param.size() == n_params)
      ? include_param
    : LogicalVector(n_params, true);

    // Lazy initialisation of the plan
    if (design_plan.size() != (size_t)n_params) {
      init_design_plan(designs);
    }

    const int T = n_trials;

    for (int i = 0; i < n_params; ++i) {
      if (!use[i]) continue;
      if (designs[i] == R_NilValue) continue;

      DesignEntry& entry = design_plan[i];
      if (!entry.valid) continue;

      // A self-intercept-only design is normally already represented by the
      // parameter value in `base`, so remapping it is redundant.  It is not
      // redundant for split transforms: the intercept is the pre-transform
      // sum and must still pass through the requested transform.
      if (entry.skip_self_intercept && !entry.split_transform) {
        continue;
      }

      const int out_idx = entry.out_idx;

      NumericMatrix design = designs[i];   // still a light handle
      const int K = design.ncol();

      // ---- Row-constant fast path -------------------------------------
      // If every trial reads the same design row and every coefficient column
      // is itself row-constant, the whole output column is one number: settle
      // it from a single representative row and fill. This is what removes the
      // n_trials-long accumulate passes for intercept-only parameters, which
      // are the majority of parameters in a typical design.
      if (track_const) {
        bool all_coef_const = true;
        for (int j = 0; j < K; ++j) {
          const int cidx = entry.coef_idx[j];
          if (cidx < 0) continue;
          if (!col_is_const(cidx)) { all_coef_const = false; break; }
        }
        const bool out_const = entry.single_cell && all_coef_const;
        col_const[out_idx] = out_const ? 1 : 0;
        col_cell[out_idx] = -1;   // re-established below if the cell path is taken
        if (out_const) {
          const int drow = entry.expand_idx.empty() ? 0 : entry.expand_idx[0];
          double* outc = &base(0, out_idx);
          // uses_self reads the pre-clear value of the output column, which is
          // constant here, so one saved scalar stands in for the whole copy.
          const double self_val = entry.uses_self ? outc[0] : 0.0;
          double pre = 0.0, post = 0.0;
          for (int j = 0; j < K; ++j) {
            const int cidx = entry.coef_idx[j];
            if (cidx < 0) continue;
            const double cv = (entry.uses_self && cidx == out_idx) ? self_val
                                                                  : base(0, cidx);
            const double contrib = cv * design(drow, j);
            if (entry.split_transform && entry.coef_in_pre_sum[j]) pre += contrib;
            else post += contrib;
          }
          const double v = entry.split_transform
            ? apply_split_transform_scalar(pre, entry.split_code,
                                           entry.split_lower, entry.split_upper) + post
            : post;
          std::fill(outc, outc + T, v);
          continue;
        }

        // ---- Cell fast path ---------------------------------------------
        // Not one value, but n_cells of them: evaluate the linear combination
        // once per design row and scatter.  Replaces K full-length
        // multiply-accumulate passes (plus the zeroing pass) with n_cells * K
        // scalar operations and one gather, and leaves the column tagged so
        // the transform and the bound check can work at cell resolution too.
        if (entry.cell_usable && all_coef_const) {
          const int n_cells = entry.n_cells;
          double post_val[EMC2_PT_MAX_CELLS];
          double pre_val[EMC2_PT_MAX_CELLS];
          double self_val[EMC2_PT_MAX_CELLS];
          double* outc = &base(0, out_idx);
          const int* ex = entry.expand_idx.data();
          for (int c = 0; c < n_cells; ++c) { post_val[c] = 0.0; pre_val[c] = 0.0; }
          // uses_self reads the pre-clear output column; it is governed by this
          // same map (only this design writes this column), so cell c's self
          // value is the one at cell c's representative row.  Snapshot them
          // before the scatter overwrites the column.
          if (entry.uses_self) {
            for (int c = 0; c < n_cells; ++c) {
              const int rep = entry.cell_rep[c];
              self_val[c] = (rep >= 0) ? outc[rep] : 0.0;
            }
          }
          // Loop order and accumulation order deliberately mirror the general
          // per-row path below, coefficient-major into a zeroed accumulator, so
          // the two routes fold the same products in the same sequence.
          for (int j = 0; j < K; ++j) {
            const int cidx = entry.coef_idx[j];
            if (cidx < 0) continue;
            const bool is_self = entry.uses_self && cidx == out_idx;
            const double cv = is_self ? 0.0 : base(0, cidx);
            double* acc = (entry.split_transform && entry.coef_in_pre_sum[j])
              ? pre_val : post_val;
            const double* dcol = &design(0, j);
            if (is_self) {
              for (int c = 0; c < n_cells; ++c) acc[c] += self_val[c] * dcol[c];
            } else {
              for (int c = 0; c < n_cells; ++c) acc[c] += cv * dcol[c];
            }
          }
          if (entry.split_transform) {
            for (int c = 0; c < n_cells; ++c) {
              post_val[c] = apply_split_transform_scalar(
                pre_val[c], entry.split_code, entry.split_lower,
                entry.split_upper) + post_val[c];
            }
          }
          for (int r = 0; r < T; ++r) outc[r] = post_val[ex[r]];
          col_cell[out_idx] = i;
          continue;
        }
      }

      // Preserve self column if needed
      std::vector<double> self_copy;
      double* out = &base(0, out_idx);

      if (entry.uses_self) {
        self_copy.resize(T);
        const double* src = out;
        std::copy(src, src + T, self_copy.begin());
      }

      // Clear output
      std::fill(out, out + T, 0.0);

      if (entry.split_transform) {
        std::vector<double> pre_acc(T, 0.0);
        std::vector<double> post_acc(T, 0.0);

        for (int j = 0; j < K; ++j) {
          int cidx = entry.coef_idx[j];
          if (cidx < 0) continue;

          const double* coef =
            (entry.uses_self && cidx == out_idx)
            ? self_copy.data()
              : &base(0, cidx);

          double* acc = entry.coef_in_pre_sum[j] ? pre_acc.data() : post_acc.data();
          // Hoist the design column pointer and lift the expand_idx.empty()
          // test out of the inner loop: design(drow, j) otherwise costs an
          // Rcpp bounds check and an i + nrow*j multiply per element.
          const double* dcol = &design(0, j);
          if (entry.expand_idx.empty()) {
            for (int r = 0; r < T; ++r) acc[r] += coef[r] * dcol[r];
          } else {
            const int* ex = entry.expand_idx.data();
            for (int r = 0; r < T; ++r) acc[r] += coef[r] * dcol[ex[r]];
          }
        }

        for (int r = 0; r < T; ++r) {
          out[r] = apply_split_transform_scalar(
            pre_acc[r], entry.split_code, entry.split_lower, entry.split_upper
          ) + post_acc[r];
        }
      } else {
        for (int j = 0; j < K; ++j) {
          int cidx = entry.coef_idx[j];
          if (cidx < 0) continue;

          const double* coef =
            (entry.uses_self && cidx == out_idx)
            ? self_copy.data()
              : &base(0, cidx);

          const double* dcol = &design(0, j);
          if (entry.expand_idx.empty()) {
            for (int r = 0; r < T; ++r) out[r] += coef[r] * dcol[r];
          } else {
            const int* ex = entry.expand_idx.data();
            for (int r = 0; r < T; ++r) out[r] += coef[r] * dcol[ex[r]];
          }
        }
      }
    }
  }

  std::unordered_set<std::string> split_transform_params() const {
    std::unordered_set<std::string> out;
    out.reserve(design_plan.size());
    for (const auto& entry : design_plan) {
      if (!entry.valid || !entry.split_transform) continue;
      out.insert(Rcpp::as<std::string>(base_names[entry.out_idx]));
    }
    return out;
  }


  // Zero the entire base matrix
  void reset_base_to_zero() {
    // Every column becomes one repeated value (zero).
    if (track_const) std::fill(col_const.begin(), col_const.end(), 1);
    const int n = n_trials;
    const int p = base.ncol();
    for (int j = 0; j < p; ++j) {
      double* col = &base(0, j);
      for (int r = 0; r < n; ++r) {
        col[r] = 0.0;
      }
    }
  }

  // Fill columns corresponding to names(p_vector) with that scalar value
  void fill_from_p_vector(const Rcpp::NumericVector& p_vector) {
    using std::string;
    using Rcpp::as;

    Rcpp::CharacterVector pv_names = p_vector.names();
    const int n_names = pv_names.size();

    for (int i = 0; i < n_names; ++i) {
      string nm = as<string>(pv_names[i]);
      auto it = name_to_base_idx.find(nm);
      if (it == name_to_base_idx.end()) {
        // p_vector may contain names that weren't used in designs, ignore them
        continue;
      }
      int j = it->second;
      double val = p_vector[i];
      const_mark(j, true);
      double* col = &base(0, j);
      for (int r = 0; r < n_trials; ++r) {
        col[r] = val;
      }
    }
  }

  // Or pass a matrix and the row to fill. Should be faster due to the pre-defined
  // pm_col_to_base_idx vector.
  // When invariant_base_indices is non-empty, those columns (already in natural
  // scale from the template particle) are preserved across the reset so the
  // invariant-parameter optimization can skip re-transforming them.
  void fill_from_particle_row(const Rcpp::NumericMatrix& particles,
                              int row,
                              const std::vector<int>& pm_col_to_base_idx,
                              const std::vector<int>& invariant_base_indices = {})
  {
    const int ncols = particles.ncol();
    const int n_trials = this->n_trials;

    if (invariant_base_indices.empty()) {
      reset_base_to_zero();
      for (int j = 0; j < ncols; ++j) {
        int base_idx = pm_col_to_base_idx[j];
        if (base_idx < 0) continue;
        double val = particles(row, j);
        const_mark(base_idx, true);
        double* col = &base(0, base_idx);
        for (int r = 0; r < n_trials; ++r) col[r] = val;
      }
    } else {
      // Save natural-scale invariant column values from the template (i=0),
      // reset, fill only non-invariant columns from the particle row (log scale),
      // then restore the invariant columns. update_pt_only will be called with
      // invariant skip, so only the non-invariant columns get exp-transformed.
      const int n_inv = static_cast<int>(invariant_base_indices.size());
      std::vector<double> saved(static_cast<size_t>(n_inv) * n_trials);
      std::vector<char> saved_const(static_cast<size_t>(n_inv), 0);
      for (int k = 0; k < n_inv; ++k) {
        const int bidx = invariant_base_indices[k];
        saved_const[k] = col_is_const(bidx) ? 1 : 0;
        for (int r = 0; r < n_trials; ++r)
          saved[static_cast<size_t>(k) * n_trials + r] = base(r, bidx);
      }
      reset_base_to_zero();
      // Build a fast-lookup set of invariant base indices to skip in the fill loop.
      std::unordered_set<int> inv_set(invariant_base_indices.begin(),
                                      invariant_base_indices.end());
      for (int j = 0; j < ncols; ++j) {
        int base_idx = pm_col_to_base_idx[j];
        if (base_idx < 0 || inv_set.count(base_idx)) continue;
        double val = particles(row, j);
        const_mark(base_idx, true);
        double* col = &base(0, base_idx);
        for (int r = 0; r < n_trials; ++r) col[r] = val;
      }
      // Restore natural-scale invariant values, and with them the row-constant
      // claims they carried: reset_base_to_zero above marked every column
      // constant, which is wrong for a restored multi-cell column.
      for (int k = 0; k < n_inv; ++k) {
        const int bidx = invariant_base_indices[k];
        const_mark(bidx, saved_const[k] != 0);
        for (int r = 0; r < n_trials; ++r)
          base(r, bidx) = saved[static_cast<size_t>(k) * n_trials + r];
      }
    }
  }

  // Planned twin of fill_from_particle_row.  The caller resolves, ONCE per
  // likelihood call, which base columns need zeroing and which (particle
  // column -> base column) pairs need filling; invariant columns appear in
  // neither list, so their natural-scale values simply survive untouched.
  // That removes the per-particle save / reset / restore round-trip and the
  // per-particle unordered_set<int> the general path builds.
  void fill_from_particle_row_planned(
      const Rcpp::NumericMatrix& particles,
      int row,
      const std::vector<int>& zero_base_idx,
      const std::vector<std::pair<int, int>>& fill_pm_to_base)
  {
    const int T = n_trials;
    // Both loops leave their column holding a single repeated value; invariant
    // columns appear in neither list and keep the flag (and the contents) they
    // were given on the template particle.
    for (std::size_t k = 0; k < zero_base_idx.size(); ++k) {
      const int b = zero_base_idx[k];
      const_mark(b, true);
      double* col = &base(0, b);
      std::fill(col, col + T, 0.0);
    }
    for (std::size_t k = 0; k < fill_pm_to_base.size(); ++k) {
      const double val = particles(row, fill_pm_to_base[k].first);
      const int b = fill_pm_to_base[k].second;
      const_mark(b, true);
      double* col = &base(0, b);
      std::fill(col, col + T, val);
    }
  }
};



std::unordered_set<std::string> param_names_excluding(const ParamTable& pt,
                                                      std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::CharacterVector names_excluding(const Rcpp::CharacterVector& names,
                                      std::initializer_list<const std::unordered_set<std::string>*> excludes);

Rcpp::NumericMatrix add_constants_columns(Rcpp::NumericMatrix p_matrix,
                                          Rcpp::NumericVector constants);
