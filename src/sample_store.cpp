#include <Rcpp.h>
#include <algorithm>

// Sample arrays are freshly allocated by extend_sampler() before a block is
// filled.  R's ordinary `[<-` duplicates those arrays once they are reached
// through a sampler list and function arguments.  These small writers update
// the owned storage directly while retaining the package's public matrix/array
// representation.

// [[Rcpp::export]]
Rcpp::List emc_clone_sample_store(Rcpp::List source) {
  // make_emc() initially replicates one sampler object across chains.  Clone
  // only the compact, one-iteration sample store before the first native write
  // so those chains do not retain shared array references.  Later block
  // extension allocates unique history arrays itself.
  return Rcpp::clone(source);
}

// [[Rcpp::export]]
void emc_set_last_slice(Rcpp::NumericVector target, int j,
                        Rcpp::NumericVector values) {
  if (j < 1) Rcpp::stop("sample index must be positive");
  Rcpp::IntegerVector dims = target.attr("dim");
  if (dims.size() < 2) Rcpp::stop("sample storage must be a matrix or array");

  R_xlen_t stride = 1;
  for (int k = 0; k < dims.size() - 1; ++k) stride *= dims[k];
  if (j > dims[dims.size() - 1]) Rcpp::stop("sample index is out of bounds");
  if (values.size() != stride) Rcpp::stop("sample slice has the wrong size");

  R_xlen_t offset = static_cast<R_xlen_t>(j - 1) * stride;
  std::copy(values.begin(), values.end(), target.begin() + offset);
}

// [[Rcpp::export]]
void emc_set_particle_slice(Rcpp::NumericVector alpha,
                            Rcpp::NumericVector subj_ll,
                            Rcpp::NumericMatrix proposals,
                            int j, int n_pars) {
  if (j < 1 || n_pars < 1) Rcpp::stop("invalid particle sample index");
  Rcpp::IntegerVector alpha_dims = alpha.attr("dim");
  Rcpp::IntegerVector ll_dims = subj_ll.attr("dim");
  if (alpha_dims.size() != 3 || ll_dims.size() != 2) {
    Rcpp::stop("particle sample storage has invalid dimensions");
  }
  int n_subjects = alpha_dims[1];
  if (alpha_dims[0] != n_pars || proposals.nrow() != n_pars + 1 ||
      proposals.ncol() != n_subjects || ll_dims[0] != n_subjects ||
      j > alpha_dims[2] || j > ll_dims[1]) {
    Rcpp::stop("particle sample slice has incompatible dimensions");
  }

  R_xlen_t alpha_offset = static_cast<R_xlen_t>(j - 1) * n_pars * n_subjects;
  for (int s = 0; s < n_subjects; ++s) {
    for (int p = 0; p < n_pars; ++p) {
      alpha[alpha_offset + p + static_cast<R_xlen_t>(s) * n_pars] = proposals(p, s);
    }
  }
  R_xlen_t ll_offset = static_cast<R_xlen_t>(j - 1) * n_subjects;
  for (int s = 0; s < n_subjects; ++s) {
    subj_ll[ll_offset + s] = proposals(n_pars, s);
  }
}

// [[Rcpp::export]]
void emc_copy_sample_prefix(Rcpp::NumericVector target,
                            Rcpp::NumericVector source) {
  if (source.size() > target.size()) {
    Rcpp::stop("sample prefix is larger than its destination");
  }
  std::copy(source.begin(), source.end(), target.begin());
}
