#ifndef utility_h
#define utility_h
#pragma once
#include <Rcpp.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cassert>


using namespace Rcpp;
/*-------------------------------------------------------*
 *  numeric constants (inline - one definition per TU)
 *-------------------------------------------------------*/
inline constexpr double inv_sqrt2_pi      = 1.0 / std::sqrt(2.0 * M_PI);
inline constexpr double neg_inv_sqrt2     = -1.0 / std::sqrt(2.0);
inline constexpr double inv_sqrt2         = 1.0 / std::sqrt(2.0);
inline constexpr double sqrt2             = std::sqrt(2.0);
inline constexpr double twoPI             = 2.0 * M_PI;
inline constexpr double fourPI            = 4.0 * M_PI;
inline constexpr double half_log_twoPI    = 0.5 * std::log(twoPI);
inline constexpr double inv_fourPI        = 1.0 / fourPI;
inline constexpr double sqrt_twoPI        = std::sqrt(2.0 * M_PI);
inline constexpr double log_sqrt_twoPI    = std::log(sqrt_twoPI);
inline constexpr double minus_inv_twoPI   = -1.0 / twoPI;
inline constexpr double inv3sqrt2pi       = 1.0 / (3.0 * std::sqrt(2.0 * M_PI));
inline constexpr double asin_one          = std::asin(1.0);
inline constexpr double inv_four_asin_one = 1.0 / (4.0 * asin_one);
const double FPM_EPSILON = 1e-12;
LogicalVector contains(CharacterVector sv, std::string txt) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    res[i] = (sv[i] == txt);
  }
  return res;
}

NumericVector vector_pow(NumericVector x1, NumericVector x2){
  NumericVector out(x1.length());
  for(unsigned int i = 0; i < out.length(); i ++){
    out[i] = pow(x1[i], x2[i]);
  }
  return(out);
}

NumericVector pnorm_multiple(NumericVector x){
  NumericVector out(x.size());
  for(int i = 0; i < x.size(); i++){
    out[i] = R::pnorm(x[i], 0, 1, TRUE, FALSE);
  }
  return out;
}

LogicalVector contains_multiple(CharacterVector sv, CharacterVector inputs) {
  LogicalVector res(sv.size());
  for (int i = 0; i < sv.size(); i ++) {
    int k = 0;
    for (int j = 0; j < inputs.size(); j++){
      if (sv[i] == inputs[j]){
        k++;
      }
    }
    res[i] = k > 0;
  }
  return res;
}

NumericMatrix submat_rcpp_col(NumericMatrix X, LogicalVector condition) {
  int n = X.nrow();
  int k = X.ncol();

  if (condition.size() != k) {
    stop("Length of logical vector must match number of columns of the matrix.");
  }

  int to_keep = sum(condition);

  // If all columns match, just return X to avoid unnecessary copying
  if (to_keep == k) {
    return X;
  }

  NumericMatrix out(n, to_keep);

  double* x_ptr = REAL(X);
  double* out_ptr = REAL(out);

  // We'll keep track of the next column in 'out' to fill
  int out_col_index = 0;

  // Extract the original column names
  CharacterVector orig_colnames = colnames(X);
  CharacterVector new_colnames(to_keep);

  for (int col = 0; col < k; col++) {
    if (condition[col]) {
      double* x_col_start = x_ptr + col * n;
      double* out_col_start = out_ptr + out_col_index * n;

      std::copy(x_col_start, x_col_start + n, out_col_start);

      // Assign the column name
      new_colnames[out_col_index] = orig_colnames[col];

      out_col_index++;
    }
  }

  // Set new column names on output
  colnames(out) = new_colnames;

  return out;
}

NumericMatrix submat_rcpp_col_by_names(NumericMatrix X, CharacterVector cols) {
  int n = X.nrow();
  int k = X.ncol();
  int m = cols.size();

  // Get the original column names from X
  CharacterVector orig_colnames = colnames(X);

  // Create a vector to store the column indices to keep
  std::vector<int> indices;
  indices.reserve(m);

  // For each name in 'cols', find the corresponding column in X
  for (int i = 0; i < m; i++) {
    bool found = false;
    for (int j = 0; j < k; j++) {
      if (as<std::string>(cols[i]) == as<std::string>(orig_colnames[j])) {
        indices.push_back(j);
        found = true;
        break;
      }
    }
    if (!found) {
      stop("Column name " + as<std::string>(cols[i]) + " not found in matrix.");
    }
  }

  // Create the output matrix with the same number of rows but only m columns
  NumericMatrix out(n, m);

  double* x_ptr = REAL(X);
  double* out_ptr = REAL(out);

  // Copy each selected column into the output matrix in the order given by 'cols'
  for (int i = 0; i < m; i++) {
    int col_idx = indices[i];
    double* x_col_start = x_ptr + col_idx * n;
    double* out_col_start = out_ptr + i * n;
    std::copy(x_col_start, x_col_start + n, out_col_start);
  }

  // Set the column names of the output matrix to be the ones provided
  colnames(out) = cols;

  return out;
}

NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) {
  int n = X.nrow(), k = X.ncol();
  int to_keep = sum(condition);
  // If all rows match, just return X (this avoids copying)
  if (to_keep == n) {
    return X;
  }
  NumericMatrix out(to_keep, k);
  for (int i = 0, j = 0; i < n; i++) {
    if (condition[i]) {
      out(j, _) = X(i, _);
      j++;
    }
  }
  colnames(out) = colnames(X);
  return out;
}


NumericVector c_expand(NumericVector x1, IntegerVector expand){
  const int n_out = expand.length();
  NumericVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

LogicalVector c_bool_expand(LogicalVector x1, IntegerVector expand){
  const int n_out = expand.length();
  LogicalVector out(n_out);
  int curr_idx;
  for(int i = 0; i < n_out; i++){
    curr_idx = expand[i] - 1; //expand created in 1-based R
    out[i] = x1[curr_idx];
  }
  return(out);
}

NumericVector c_add_vectors(NumericVector x1, NumericVector x2){
  if(is_na(x2)[0] ){
    return(x1);
  }
  NumericVector output(x1.size() + x2.size());
  std::copy(x1.begin(), x1.end(), output.begin());
  std::copy(x2.begin(), x2.end(), output.begin() + x1.size());
  CharacterVector all_names(x1.size() + x2.size());
  CharacterVector x1_names = x1.names();
  CharacterVector x2_names = x2.names();
  std::copy(x1_names.begin(), x1_names.end(), all_names.begin());
  std::copy(x2_names.begin(), x2_names.end(), all_names.begin() + x1.size());
  output.names() = all_names;
  return output;
}

// [[Rcpp::export]]
CharacterVector c_add_charvectors(CharacterVector x, CharacterVector y) {
  // Create a new vector of length = length(x) + length(y)
  CharacterVector z(x.size() + y.size());
  // Copy x into z
  std::copy(x.begin(), x.end(), z.begin());
  // Copy y into z after x
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  return z;
}

// Custom hash for a vector<double>
struct RowHash {
  std::size_t operator()(const std::vector<double> &v) const {
    std::size_t seed = 0;
    std::hash<double> hash_double;
    for (double d : v) {
      // A standard hash combination approach (based on boost::hash_combine)
      seed ^= hash_double(d) + 0x9e3779b97f4a7c16ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

// Custom equality for a vector<double>
struct RowEqual {
  bool operator()(const std::vector<double> &a, const std::vector<double> &b) const {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }
};

Rcpp::LogicalVector duplicated_matrix(Rcpp::NumericMatrix x) {
  int n = x.nrow();
  int m = x.ncol();

  Rcpp::LogicalVector dup(n);
  dup.fill(false);

  // Pre-allocate storage to reduce reallocations
  std::vector<double> buffer(m);
  std::unordered_set<std::vector<double>, RowHash, RowEqual> seen;
  seen.reserve(n);

  for (int i = 0; i < n; i++) {
    // Copy row data into buffer
    for (int j = 0; j < m; j++) {
      buffer[j] = x(i, j);
    }

    // Check if we have seen this row before
    if (seen.find(buffer) != seen.end()) {
      dup[i] = true;
    } else {
      // Insert a copy of the current row
      seen.insert(std::vector<double>(buffer.begin(), buffer.end()));
    }
  }

  return dup;
}

IntegerVector cumsum_logical(LogicalVector x) {
  int n = x.size();
  IntegerVector out(n);
  int running_total = 0;

  for (int i = 0; i < n; i++) {
    // Add 1 if TRUE, else 0
    if (x[i]) {
      running_total += 1;
    }
    out[i] = running_total;
  }

  return out;
}

IntegerVector which_rcpp(LogicalVector x) {
  int n = x.size();
  int count = 0;
  // First pass: count how many TRUE
  for (int i = 0; i < n; i++) {
    if (x[i]) count++;
  }

  // Allocate the output vector
  IntegerVector out(count);

  // Second pass: fill the output with the indices of TRUE values
  for (int i = 0, j = 0; i < n; i++) {
    if (x[i]) {
      out[j] = i;
      j++;
    }
  }

  return out;
}

LogicalVector lr_all(LogicalVector ok, int n_side)
{
  const R_xlen_t n = ok.size();
  if (n % n_side)
    stop("Vector length is not a multiple of n_side");

  // Allocate without initializing for speed
  LogicalVector out(no_init(n));

  const int* x = LOGICAL(ok);   // input pointer
  int*       o = LOGICAL(out);  // output pointer

  for (R_xlen_t i = 0; i < n; i += n_side) {

    int state = TRUE;                 // optimistic

    // inspect one "column"
    for (int j = 0; j < n_side; ++j) {
      const int v = x[i + j];
      if (v == FALSE) {             // any FALSE trumps everything
        state = FALSE;
        break;
      }
      if (v == NA_LOGICAL)          // remember NA unless FALSE appears
        state = NA_LOGICAL;
    }

    // replicate result to each entry of the column
    std::fill(o + i, o + i + n_side, state);
  }
  return out;
}

// For do_bounds
struct BoundSpec {
  int col_idx;         // which column in 'pars'
  double min_val;
  double max_val;
  bool has_exception;
  double exception_val;
};


std::vector<BoundSpec> make_bound_specs(NumericMatrix minmax,
                                        CharacterVector minmax_colnames,
                                        NumericMatrix pars,
                                        List bound)
{
  // 1) Build a map from param-name -> column index in 'pars'
  CharacterVector pcolnames = colnames(pars);
  std::unordered_map<std::string, int> colMap;
  for (int j = 0; j < pcolnames.size(); j++) {
    colMap[ Rcpp::as<std::string>(pcolnames[j]) ] = j;
  }

  // 2) Build a map from param-name -> exception value
  bool has_exception = bound.containsElementNamed("exception") && !Rf_isNull(bound["exception"]);
  std::unordered_map<std::string, double> exceptionMap;
  if (has_exception) {
    NumericVector except_vec = bound["exception"];
    CharacterVector except_names = except_vec.names();
    for (int i = 0; i < (int)except_vec.size(); i++) {
      exceptionMap[ Rcpp::as<std::string>(except_names[i])] = except_vec[i];
    }
  }

  // 3) Create BoundSpec for each column in minmax
  int ncols = minmax_colnames.size();
  std::vector<BoundSpec> specs(ncols);
  for (int j = 0; j < ncols; j++) {
    std::string var_name = Rcpp::as<std::string>(minmax_colnames[j]);

    // Fill the struct
    BoundSpec s;
    s.col_idx     = colMap[var_name];
    s.min_val     = minmax(0, j);
    s.max_val     = minmax(1, j);

    auto it = exceptionMap.find(var_name);
    if (it != exceptionMap.end()) {
      s.has_exception = true;
      s.exception_val = it->second;
    } else {
      s.has_exception = false;
      s.exception_val = NA_REAL;  // or 0
    }
    specs[j] = s;
  }
  return specs;
}

// For transforms
enum TransformCode {
  IDENTITY = 0,
  EXP      = 1,
  PNORM    = 2
};

struct TransformSpec {
  int col_idx;        // which column in 'pars'
  TransformCode code; // e.g. EXP, PNORM, ...
  double lower;
  double upper;
};

std::vector<TransformSpec> make_transform_specs(NumericMatrix pars, List transform)
{
  // gather 'func', 'lower', 'upper'
  CharacterVector func_charvec = transform["func"];
  NumericVector lower_numvec   = transform["lower"];
  NumericVector upper_numvec   = transform["upper"];

  // Build a map param_name -> code
  std::unordered_map<std::string,TransformCode> codeMap;
  {
    // e.g. "param_name" -> "exp" or "pnorm" in func_charvec
    CharacterVector fnames = func_charvec.names();
    for (int i = 0; i < func_charvec.size(); i++) {
      std::string name = Rcpp::as<std::string>(fnames[i]);
      std::string f    = Rcpp::as<std::string>(func_charvec[i]);
      if (f == "exp") {
        codeMap[name] = EXP;
      } else if (f == "pnorm") {
        codeMap[name] = PNORM;
      } else {
        codeMap[name] = IDENTITY;
      }
    }
  }

  // Build param_name -> (lower, upper)
  std::unordered_map<std::string,std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lower_numvec.names();
    for (int i = 0; i < lower_numvec.size(); i++) {
      std::string nm = Rcpp::as<std::string>(ln[i]);
      boundMap[nm].first = lower_numvec[i];
    }
    CharacterVector un = upper_numvec.names();
    for (int i = 0; i < upper_numvec.size(); i++) {
      std::string nm = Rcpp::as<std::string>(un[i]);
      boundMap[nm].second = upper_numvec[i];
    }
  }

  // Now fill specs for each col in pars
  int ncol = pars.ncol();
  std::vector<TransformSpec> specs(ncol);

  CharacterVector cparnames = colnames(pars);
  for (int j = 0; j < ncol; j++) {
    std::string colname = Rcpp::as<std::string>(cparnames[j]);
    TransformSpec sp;
    sp.col_idx = j;
    sp.code    = codeMap[colname];
    auto it    = boundMap.find(colname);
    if (it != boundMap.end()) {
      sp.lower = it->second.first;
      sp.upper = it->second.second;
    } else {
      sp.lower = 0.0;  // or default
      sp.upper = 1.0;  // or default
    }
    specs[j] = sp;
  }

  return specs;
}


// For pretransform
enum PreTFCode { PTF_EXP = 1, PTF_PNORM = 2, PTF_NONE = 0 };

struct PreTransformSpec {
  int index;       // index in p_vector
  PreTFCode code;
  double lower;
  double upper;
  // Possibly store the original name if needed
};

std::vector<PreTransformSpec> make_pretransform_specs(NumericVector p_vector, List transform)
{
  // e.g. transform["func"], transform["lower"], transform["upper"]
  CharacterVector func   = transform["func"];
  NumericVector lowervec = transform["lower"];
  NumericVector uppervec = transform["upper"];

  // Build a map param_name -> code
  std::unordered_map<std::string,PreTFCode> codeMap;
  CharacterVector fnames = func.names();
  for (int i = 0; i < func.size(); i++) {
    std::string name = Rcpp::as<std::string>(fnames[i]);
    std::string f    = Rcpp::as<std::string>(func[i]);
    if (f == "exp") {
      codeMap[name] = PTF_EXP;
    } else if (f == "pnorm") {
      codeMap[name] = PTF_PNORM;
    } else {
      codeMap[name] = PTF_NONE;
    }
  }

  // Build a map param_name -> (lower, upper)
  std::unordered_map<std::string, std::pair<double,double>> boundMap;
  {
    CharacterVector ln = lowervec.names();
    for (int i = 0; i < lowervec.size(); i++) {
      boundMap[ Rcpp::as<std::string>(ln[i]) ].first = lowervec[i];
    }
    CharacterVector un = uppervec.names();
    for (int i = 0; i < uppervec.size(); i++) {
      boundMap[ Rcpp::as<std::string>(un[i]) ].second = uppervec[i];
    }
  }

  // Now create PreTransformSpec for each element in p_vector
  CharacterVector p_names = p_vector.names();
  int n = p_vector.size();
  std::vector<PreTransformSpec> specs(n);
  for (int i = 0; i < n; i++) {
    std::string pname = Rcpp::as<std::string>(p_names[i]);
    PreTransformSpec s;
    s.index = i;
    s.code = codeMap[pname];
    auto it = boundMap.find(pname);
    if (it != boundMap.end()) {
      s.lower = it->second.first;
      s.upper = it->second.second;
    } else {
      s.lower = 0.0;
      s.upper = 1.0;
    }
    specs[i] = s;
  }
  return specs;
}

//TODO check if we can delete this
static const Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm");
// Make pmvnorm available from R package, used for RDMSWTN
// [[Rcpp::export]]
NumericVector pmvnorm_cpp(NumericVector upper, NumericMatrix corr) {
  Function pmvnorm = mvtnorm["pmvnorm"];
  
  int n = upper.length();
  NumericVector lower(n, R_NegInf);      // Lower: -Inf for each dimension
  NumericVector mean(n, 0.0);            // Mean: 0 for each dimension

  NumericVector res = pmvnorm(_["lower"] = lower,
                              _["upper"] = upper,
                              _["mean"] = mean,
                              _["corr"] = corr);
  return res;
}
// 1-liner that loads the namespace if necessary and returns it
static const Rcpp::Environment statmodNS = Rcpp::Environment::namespace_env("statmod");

// Function object lives for the life-time of the DLL
static const Rcpp::Function gauss_quad = statmodNS["gauss.quad"];

static Rcpp::List gl = gauss_quad(20, "legendre");
static Rcpp::List gh = gauss_quad(5, "hermite");

// Numerically stable log(e^a + e^b)
double log_sum_exp(double a, double b) {
    if (std::isinf(a) && a < 0) return b;
    if (std::isinf(b) && b < 0) return a;
    if (a > b) {
        return a + std::log1p(std::exp(b - a)); // log1p(x) = log(1+x)
    } else {
        return b + std::log1p(std::exp(a - b));
    }
}

// Useful Distribution Functions
// Normal PDF
double gaussian_pdf(double x, double mean=0.0, double var = 1.0)
{
    if (std::isinf(x)) return 0.;
	if (var <= 0.0) return 0.0;
	const double inv_sqrt_var = 1.0/std::sqrt(var);
    const double z = (x - mean)*inv_sqrt_var;
    if (z < 5. && z > -5.) {
        return inv_sqrt_var*inv_sqrt2_pi * std::exp(-0.5 * z * z);
    }
    else {
        double z1 = std::floor(z * 0x1.0p16 + 0.5) * 0x1.0p-16;
        double z2 = z - z1;
        return inv_sqrt_var*inv_sqrt2_pi * (std::exp(-0.5 * z1 * z1) * std::exp((-0.5 * z2 - z1) * z2));
    }
}

// Normal CDF
/* Adapted from cephes' 'ndtr' */
inline double gaussian_cdf(double x, double mean=0.0, double var=1.0) {
    if (std::isinf(x)) return (x > 0.0) ? 1.0 : 0.0;
	if (var <= 0.0) return (x < mean) ? 0.0 : 1.0;
    const double z = (x - mean) / std::sqrt(var)  * inv_sqrt2;
    const double az = std::fabs(z);

    // Near zero: symmetric and well-conditioned
    if (az < inv_sqrt2) {                          
        return 0.5 + 0.5 * std::erf(z);
    }

    // Tails: use erfc with nonnegative argument, avoid subtractive cancellation
    if (z > 0.0) {
        return 1.0 - 0.5 * std::erfc(az);
    } else {
        return 0.5 * std::erfc(az);                // no subtraction from 1
    }
}

// Heat kernel G*(t, delta) = N(delta | 0, t)
inline double Gstar(double var, double delta) {
  if (var <= 0.0) return 0.0;
  return gaussian_pdf(delta, 0.0, var);
}

// CDF of heat kernel N(mean, t) at x
inline double Gstar_CDF(double var, double mean, double x) {
	if (var <= 0.0) return (x < mean) ? 0.0 : 1.0;
    return gaussian_cdf(x, mean, var);
}

// Definite integral ∫_{x_lo}^{x_hi} N(mean, t) dx
inline double Gstar_Integral(double t, double mean, double x_lo, double x_hi) {
    if (t <= 0.0) {
        // Dirac delta at mean; match your half-open convention (x_lo, x_hi]
        return (mean > x_lo && mean <= x_hi) ? 1.0 : 0.0;
    }
    if (x_hi <= x_lo) return 0.0;
    return Gstar_CDF(t, mean, x_hi) - Gstar_CDF(t, mean, x_lo);
}

/**
 * @brief Computes ∫[x_lo, x_hi] Φ(ax+c) dx
 * Where Φ is the standard Normal CDF.
 */
inline double integrate_gaussian_cdf(double a, double c, double x_lo, double x_hi) {
    if (std::abs(a) < FPM_EPSILON) {
        // Fallback: integral of a constant
        return Gstar_CDF(1.0, 0.0, c) * (x_hi - x_lo);
    }

    // Indefinite integral is (1/a) * [u*Φ(u) + φ(u)], where u = ax+c
    auto eval_integral = [&](double x){
        const double u = a * x + c;
        // Gstar_CDF is Φ(u), Gstar is φ(u) (with var=1)
        return (u * Gstar_CDF(1.0, 0.0, u) + Gstar(1.0, u)) / a;
    };
    
    return eval_integral(x_hi) - eval_integral(x_lo);
}

/**
 * @brief Computes ∫[x_lo, x_hi] exp(kx) * Φ(ax+c) dx
 */
inline double integrate_exp_times_normal_cdf(double k, double a, double c, double x_lo, double x_hi) {
    if (std::abs(k) < FPM_EPSILON) { // Fallback to non-exp integral
        return integrate_gaussian_cdf(a, c, x_lo, x_hi);
    }
    if (std::abs(a) < FPM_EPSILON) { // Fallback: integral of exp(kx) * C
        const double C = Gstar_CDF(1.0, 0.0, c);
        return C * (std::exp(k * x_hi) - std::exp(k * x_lo)) / k;
    }

    // Standard identity for this integral
    const double exp_factor = std::exp(k * k / (2.0 * a * a) - k * c / a);
    
    auto eval_integral = [&](double x){
        const double u = a * x + c;
        const double u_shifted = a * x + c - k / a;
        
        const double term1 = (std::exp(k * x) / k) * Gstar_CDF(1.0, 0.0, u);
        const double term2 = (1.0 / k) * exp_factor * Gstar_CDF(1.0, 0.0, u_shifted);
        
        return term1 - term2;
    };
    
    return eval_integral(x_hi) - eval_integral(x_lo);
}

double sqrt_pos(double x) {
  return std::sqrt(x > 0.0 ? x : 0.0);
}

double safe_exp(double x) {
  // avoid denormals/overflow
  if (x < -745.0) return 0.0;
  if (x >  709.0) return std::exp(709.0);
  return std::exp(x);
}
double safe_div(double num, double den) {
  if (!std::isfinite(num) || !std::isfinite(den) || den <= 0.0) return 0.0;
  return num / den;
}

// Exponential decay function
inline double exp_decay_scalar(double t, double x0, double xinf, double tau, double p) {
  if (tau <= 0.0 || p <= 0.0) return xinf;
  x0 = std::abs(x0);
  xinf = std::abs(xinf);
  const double amp = x0 - xinf;
  if (!std::isfinite(t) || t <= 0.0) {
    return xinf + amp;
  }
  double y = p * (std::log(t) - std::log(tau));
  double s = safe_exp(y);
  double exp_term = safe_exp(-s);
  return xinf + amp * exp_term;
}

// [[Rcpp::export]]
NumericVector exp_decay(NumericVector t, double x0, double xinf, double tau, double p) {
  if (tau <= 0.0 || p <= 0.0) stop("tau and p must be > 0.");
  int n = t.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double ti = t[i];
    if (NumericVector::is_na(ti)) {
      out[i] = NA_REAL;
    } else {
      out[i] = exp_decay_scalar(ti, x0, xinf, tau, p);
    }
  }
  return out;
}

// Bivariate Normal CDF Functions
// Functions distributed in this file taken from https://github.com/david-cortes/approxcdf
// Copyright 2022 David Cortes under BSD-3 license
constexpr static const double GL5_w[] = {
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889
};
constexpr static const double GL5_x[] = {
    0.9061798459386640, 0.5384693101056831, 0.0000000000000000
};

constexpr static const double GL6_w[] = {
    0.1713244923791704, 0.3607615730481386, 0.4679139345726910
};
constexpr static const double GL6_x[] = {
    0.9324695142031521, 0.6612093864662645, 0.2386191860831969
};

constexpr static const double GL8_w[] = {
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620
};
constexpr static const double GL8_x[] = {
    0.9602898564975363, 0.7966664774136267,
    0.5255324099163290, 0.1834346424956498
};

constexpr static const double GL12_w[] = {
    0.0471753363865118, 0.1069393259953184, 0.1600783285433462,
    0.2031674267230659, 0.2334925365383548, 0.2491470458134028
};
constexpr static const double GL12_x[] = {
    0.9815606342467192, 0.9041172563704749, 0.7699026741943047,
    0.5873179542866175, 0.3678314989981802, 0.1252334085114689
};

constexpr static const double GL16_w[] = {
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
    0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
    0.1826034150449236, 0.1894506104550685
};
constexpr static const double GL16_x[] = {
    0.9894009349916499, 0.9445750230732326, 0.8656312023878318,
    0.7554044083550030, 0.6178762444026438, 0.4580167776572274,
    0.2816035507792589, 0.0950125098376374
};

constexpr static const double GL20_w[] = {
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
    0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
    0.1316886384491766, 0.1420961093183820, 0.1491729864726037,
    0.1527533871307258
};
constexpr static const double GL20_x[] = {
    0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
    0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
    0.5108670019508271, 0.3737060887154195, 0.2277858511416451,
    0.0765265211334973
};

constexpr static const double GL24_w[] = {
    0.0123412297999872, 0.0285313886289337, 0.0442774388174198,
    0.0592985849154368, 0.0733464814110803, 0.0861901615319533,
    0.0976186521041139, 0.1074442701159656, 0.1155056680537256,
    0.1216704729278034, 0.1258374563468283, 0.1279381953467522
};
constexpr static const double GL24_x[] = {
    0.9951872199970213, 0.9747285559713095, 0.9382745520027328,
    0.8864155270044011, 0.8200019859739029, 0.7401241915785544,
    0.6480936519369755, 0.5454214713888396, 0.4337935076260451,
    0.3150426796961634, 0.1911188674736163, 0.0640568928626056
};

enum GLApprox {GL6=3, GL8=4, GL12=6, GL16=8, GL20=10, GL24=12}; 

constexpr static const double GL8_xp[] = {
    0.5 * GL8_x[0] + 0.5, 0.5 * GL8_x[1] + 0.5, 0.5 * GL8_x[2] + 0.5,
    0.5 * GL8_x[3] + 0.5
};
constexpr static const double GL8_xn[] = {
    -0.5 * GL8_x[0] + 0.5, -0.5 * GL8_x[1] + 0.5, -0.5 * GL8_x[2] + 0.5,
    -0.5 * GL8_x[3] + 0.5
};

constexpr static const double GL16_xp[] = {
    0.5 * GL16_x[0] + 0.5, 0.5 * GL16_x[1] + 0.5, 0.5 * GL16_x[2] + 0.5,
    0.5 * GL16_x[3] + 0.5, 0.5 * GL16_x[4] + 0.5, 0.5 * GL16_x[5] + 0.5,
    0.5 * GL16_x[6] + 0.5, 0.5 * GL16_x[7]
};
constexpr static const double GL16_xn[] = {
    -0.5 * GL16_x[0] + 0.5, -0.5 * GL16_x[1] + 0.5, -0.5 * GL16_x[2] + 0.5,
    -0.5 * GL16_x[3] + 0.5, -0.5 * GL16_x[4] + 0.5, -0.5 * GL16_x[5] + 0.5,
    -0.5 * GL16_x[6] + 0.5, -0.5 * GL16_x[7]
};

#ifndef _OPENMP
constexpr static const double GL5_w_div4pi[] = {
    GL5_w[0] / (4. * M_PI), GL5_w[1] / (4. * M_PI), GL5_w[2] / (4. * M_PI)
};
constexpr static const double GL5_xp[] = {
    0.5 * GL5_x[0] + 0.5, 0.5 * GL5_x[1] + 0.5, 0.5 * GL5_x[2] + 0.5
};
constexpr static const double GL5_xn[] = {
    -0.5 * GL5_x[0] + 0.5, -0.5 * GL5_x[1] + 0.5, -0.5 * GL5_x[2] + 0.5
};
#else
constexpr static const double GL8_w_div4pi[] = {
    GL8_w[0] / (4. * M_PI), GL8_w[1] / (4. * M_PI),
    GL8_w[2] / (4. * M_PI), GL8_w[3] / (4. * M_PI)
};
#endif

/* Bivariate normal CDF.
   Algorithm:  Drezner (1978) 5-point GL with the p-split
               refined by Drezner & Wesolowsky (1990);
               parameter cut-off as in West (2004).
   Expected accuracy ~ 1e-6.
*/
double norm_lcdf_2d_fast(double x1, double x2, double rho)
{
    double x12 = 0.5 * (x1*x1 + x2*x2);
    
    double out = 0;
    double r1, x3;
    if (std::fabs(rho) >= 0.7) {
        double r2 = 1. - rho*rho;
        double r3 = std::sqrt(r2);
        if (rho < 0) {
            x2 = -x2;
        }
        x3 = x1*x2;
        double x7 = std::exp(-0.5 * x3);
        if (r2) {
            double x6 = std::fabs(x1 - x2);
            double x5 = 0.5 * x6*x6;
            x6 /= r3;
            double aa = 0.5 - 0.125*x3;
            double ab = 3. - 2.*aa*x5;
            out = inv3sqrt2pi * (
                x6 * ab * gaussian_cdf(-x6) -
                std::exp(-x5/r2) * std::fma(aa, r2, ab) * inv_sqrt2_pi
            );
            double rr;
            double nr1, nrr, nr2;
            #ifndef _OPENMP
            #pragma GCC unroll 2
            for (int ix = 0; ix < 2; ix++) {
                r1 = r3 * GL5_xp[ix];
                rr = r1*r1;
                r2 = std::sqrt(1. - rr);

                nr1 = r3 * GL5_xn[ix];
                nrr = nr1*nr1;
                nr2 = std::sqrt(1. - nrr);
                
                out -= GL5_w_div4pi[ix] * (
                    std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr) +
                    std::exp(-x5/nrr) * (std::exp(-x3/(1. + nr2))/nr2/x7 - 1.- aa*nrr)
                );
            }
            r1 = r3 * GL5_xp[2];
            rr = r1*r1;
            r2 = std::sqrt(1. - rr);
            out -= GL5_w_div4pi[2] * std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr);
            #else
            #ifndef _MSC_VER
            #pragma omp simd
            #endif
            for (int ix = 0; ix < 4; ix++) {
                r1 = r3 * GL8_xp[ix];
                rr = r1*r1;
                r2 = std::sqrt(1. - rr);

                nr1 = r3 * GL8_xn[ix];
                nrr = nr1*nr1;
                nr2 = std::sqrt(1. - nrr);
                
                out -= GL8_w_div4pi[ix] * (
                    std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr) +
                    std::exp(-x5/nrr) * (std::exp(-x3/(1. + nr2))/nr2/x7 - 1.- aa*nrr)
                );
            }
            #endif
        }
        if (rho > 0) {
            out = std::fma(out, r3*x7, gaussian_cdf(-std::fmax(x1, x2)));
        }
        else {
            out = std::fmax(0., gaussian_cdf(-x1) - gaussian_cdf(-x2)) - out*r3*x7;
        }
        return out;
    }
    else {
        x3 = x1*x2;
        double rr2;
        double nr1, nrr2;
        #ifndef _OPENMP
        #pragma GCC unroll 2
        for (int ix = 0; ix < 2; ix++) {
            r1 = rho * GL5_xp[ix];
            rr2 = 1. - r1*r1;

            nr1 = rho * GL5_xn[ix];
            nrr2 = 1. - nr1*nr1;

            out += GL5_w_div4pi[ix] * (
                std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2) +
                std::exp((nr1*x3 - x12) / nrr2) / std::sqrt(nrr2)
            );
        }
        r1 = rho * GL5_xp[2];
        rr2 = 1. - r1*r1;
        out += GL5_w_div4pi[2] * std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2);
        #else
        #ifndef _MSC_VER
        #pragma omp simd
        #endif
        for (int ix = 0; ix < 4; ix++) {
            r1 = rho * GL8_xp[ix];
            rr2 = 1. - r1*r1;

            nr1 = rho * GL8_xn[ix];
            nrr2 = 1. - nr1*nr1;

            out += GL8_w_div4pi[ix] * (
                std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2) +
                std::exp((nr1*x3 - x12) / nrr2) / std::sqrt(nrr2)
            );
        }
        #endif
        return std::fma(out, rho, gaussian_cdf(-x1) * gaussian_cdf(-x2));
    }
}

double norm_cdf_2d_fast(double x1, double x2, double rho)
{
    return norm_lcdf_2d_fast(-x1, -x2, rho);
}

/* Tsay, Wen-Jen, and Peng-Hsuan Ke.
   "A simple approximation for the bivariate normal integral."
   Communications in Statistics-Simulation and Computation (2021): 1-14. */
constexpr const static double c1 = -1.0950081470333;
constexpr const static double c2 = -0.75651138383854;
double norm_cdf_2d_vfast(double x1, double x2, double rho)
{
    if (std::fabs(rho) <= std::numeric_limits<double>::epsilon()) {
        return gaussian_cdf(x1) * gaussian_cdf(x2);
    }

    double denom = std::sqrt(1 - rho * rho);
    double a = -rho / denom;
    double b = x1 / denom;
    double aq_plus_b = a*x2 + b;

    if (a > 0) {
        if (aq_plus_b >= 0) {
            double aa = a * a;
            double a_sq_c1 = aa*c1;
            double a_sq_c2 = aa*c2;
            double sqrt2b = sqrt2*b;
            double sqrt2x2 = sqrt2*x2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(1. - a_sq_c2);
            double twicea_sqrt_recpr_a_sq_c2 = 2.*a*sqrt_recpr_a_sq_c2;
            double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
            double t1 = a_sq_c1*c1 + 2.*b*b*c2;
            double t2 = 2.*sqrt2b*c1;
            double t3 = 4. - 4.*a_sq_c2;

            return
                0.5 * (std::erf(x2 / sqrt2) + std::erf(b / (sqrt2*a))) +
                temp
                    * std::exp((t1 - t2) / t3)
                    * (1. - std::erf((sqrt2b - a_sq_c1) / twicea_sqrt_recpr_a_sq_c2)) -
                temp
                    * std::exp((t1 + t2) / t3)
                    * (
                        std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 - a*c1) / (2.*sqrt_recpr_a_sq_c2)) +
                        std::erf((a_sq_c1 + sqrt2b) / twicea_sqrt_recpr_a_sq_c2)
                    );

        }
        else {
            double sqrt2b = sqrt2*b;
            double sqrt2x2 = sqrt2*x2;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;

            return
                (1. / (4. * sqrt_recpr_a_sq_c2)) *
                std::exp((a_c1*a_c1 - 2.*sqrt2b*c1 + 2*b*b*c2) / (4.*recpr_a_sq_c2)) *
                (1. + std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)));
        }
    }
    else {
        if (aq_plus_b >= 0) {
            double sqrt2b = sqrt2*b;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;
            double sqrt2_x2 = sqrt2*x2;

            return
                0.5 + 0.5 * std::erf(x2 / sqrt2) -
                (1. / (4. * sqrt_recpr_a_sq_c2)) *
                std::exp((a_c1*a_c1 + 2.*sqrt2b*c1 + 2.*b*b*c2) / (4.*recpr_a_sq_c2)) *
                (1. + std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 - a_c1) / (2.*sqrt_recpr_a_sq_c2)));
        }
        else {
            double sqrt2a = sqrt2*a;
            double sqrt2b = sqrt2*b;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;
            double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
            double t1 = a_c1*a_c1 + 2.*b*b*c2;
            double t2 = 2.*sqrt2b*c1;
            double t3 = 4.*recpr_a_sq_c2;
            double sqrt2_x2 = sqrt2*x2;

            return
                0.5 - 0.5 * std::erf(b / sqrt2a) -
                temp
                    * std::exp((t1 + t2) / t3)
                    * (1. - std::erf((sqrt2b + a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))) +
                temp
                    * std::exp((t1 - t2) / t3)
                    * (
                        std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)) +
                        std::erf((sqrt2b - a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))
                    );
        }
    }
} 

[[gnu::flatten]]
double norm_lcdf_2d(double x1, double x2, double rho)
{
    double abs_rho = std::fabs(rho);

    double out = 0;
    double hk;
    if (abs_rho < 0.925) {
        if (abs_rho > std::numeric_limits<double>::epsilon()) {
            hk = x1 * x2;
            double hs = 0.5 * (x1*x1 + x2*x2);
            double asr = std::asin(rho);
            double asr_half = 0.5 * asr;
            double sn1, sn2;

            int gl_dim;
            #ifndef _OPENMP
            const double* GL_wtable;
            const double* GL_xtable;
            #endif
            if (abs_rho < 0.3) {
                gl_dim = (int)GL6;
                #ifndef _OPENMP
                GL_wtable = GL6_w;
                GL_xtable = GL6_x;
                #endif
            }
            else if (abs_rho < 0.5) {
                gl_dim = (int)GL12;
                #ifndef _OPENMP
                GL_wtable = GL12_w;
                GL_xtable = GL12_x;
                #endif
            }
            else {
                gl_dim = (int)GL20;
                #ifndef _OPENMP
                GL_wtable = GL20_w;
                GL_xtable = GL20_x;
                #endif
            }
            
            #ifndef _OPENMP
            #pragma GCC unroll 3
            for (int ix = 0; ix < gl_dim; ix++) {
                sn1 = std::sin(asr_half * (1. + GL_xtable[ix]));
                sn2 = std::sin(asr_half * (1. - GL_xtable[ix]));
                out += GL_wtable[ix] * (
                    std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                    std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                );
            }
            #else
            switch (gl_dim) {
                case GL6: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 4; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL8_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL8_x[ix]));
                        out += GL8_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
                case GL12: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 8; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL16_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL16_x[ix]));
                        out += GL16_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
                case GL20: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 12; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL24_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL24_x[ix]));
                        out += GL24_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
            }
            #endif
            out *= asr / fourPI;
        }
        out = std::fma(gaussian_cdf(-x1), gaussian_cdf(-x2), out);
    }
    else {

        if (rho < 0) {
            x2 = -x2;
        }
        if (abs_rho < 1) {
            hk = x1 * x2;
            double as = std::fma(-rho, rho, 1.);
            double a = std::sqrt(as);
            double b;
            double bs = (x1 - x2) * (x1 - x2);
            double c = std::fma(-0.125, hk, 0.5);
            double d = std::fma(-0.0625, hk, 0.75);
            double asr = -0.5 * (hk + bs / as);
            double rfdbs = std::fma(-d, bs, 5.)*(1./15.);
            if (asr > -100.) {
                out = a * std::exp(asr) * (1. - c * (bs - as) * rfdbs + 0.2*c*d*as*as);
            }
            if (hk > -100.) {
                b = std::sqrt(bs);
                out -= std::exp(-0.5 * hk) * sqrt_twoPI * gaussian_cdf(-b / a) * b * (1. - c*bs*rfdbs);
            }
            a *= 0.5;
            double xs;
            double rs;
            double temp;

            #ifndef _OPENMP
            #pragma GCC unroll 10
            for (int ix = 0; ix < 10; ix++) {
                temp = a * (1. + GL20_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                if (asr > -100.) {
                    out += a * GL20_w[ix] * std::exp(asr) *
                           (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
                }

                temp = a * (1. - GL20_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                if (asr > -100.) {
                    out += a * GL20_w[ix] * std::exp(asr) *
                           (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
                }
            }
            #else
            #ifndef _MSC_VER
            #pragma omp simd
            #endif
            for (int ix = 0; ix < 12; ix++) {
                temp = a * (1. + GL24_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                out += a * GL24_w[ix] * std::exp(asr) *
                       (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));

                temp = a * (1. - GL24_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                out += a * GL24_w[ix] * std::exp(asr) *
                       (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
            }
            #endif
            out *= minus_inv_twoPI;
        }
        if (rho > 0) {
            out += gaussian_cdf(-std::fmax(x1, x2));
        }
        else {
            out = -out;
            if (x2 > x1) {
                if (x1 < 0) {
                    out += gaussian_cdf(x2) - gaussian_cdf(x1);
                }
                else {
                    out += gaussian_cdf(-x1) - gaussian_cdf(-x2);
                }
            }
        }
    }
    return out;
}

double norm_cdf_2d(double x1, double x2, double rho)
{
    return norm_lcdf_2d(-x1, -x2, rho);
}

// [[Rcpp::export]]
double pbvn_tsay(double h, double k, double rho) {
    return norm_cdf_2d_vfast(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_tvpack(double h, double k, double rho) {
    return norm_cdf_2d(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_drezner(double h, double k, double rho) {
    return norm_cdf_2d_fast(h, k, rho);
}

#endif


