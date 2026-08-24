#ifndef utils_h
#define utils_h

#include <RcppArmadillo.h>
#include "exgaussian_functions.h"
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <limits>
#include "col_registry.h"
#include "fpe_race.h"
#include "fpe_bou.h"
#include "model_RLF.h"
#include "utility_functions.h"
#include "model_RDM.h"
#include "model_RDM_adapters.h"
#include "model_RDMSWTN.h"
#include "model_GBM.h"
#include "model_LBA.h"
#include "model_BAwL.h"
#include "model_BAwD.h"
#include "model_BAwF.h"
#include "model_BAwR.h"
#include "model_BTAwL.h"
#include "race_contract.h"

using namespace Rcpp;

#include "timer_helpers.h"

// Included here, not at the top: the ROU kernels need ContextForRaceModels and
// the raw_log_* helpers above.  (The cache TYPE they store in the context comes
// from fpe_race.h, which has no such dependency and is included at the top.)
#include "model_ROU.h"
#include "model_ROUp.h"
#include "model_GOM.h"
#include "model_RLF_kernels.h"
// Likewise: the BOU primitives need ContextForDDMModels defined above.
#include "model_BOU.h"
// FRQ needs nothing from this file -- it is pure Rmath -- but it carries
// [[Rcpp::export]] entry points, so it must be seen by exactly one translation
// unit and this header is included by exactly one (particle_ll.cpp).
#include "model_FRQ.h"
#include "model_PCOUNTER.h"


// Helper to safely get a column from a DataFrame with a default value if missing
// Also fills NA values with the default for backward compatibility.
inline Rcpp::NumericVector get_col_with_default(const Rcpp::DataFrame& df, const std::string& name, double default_val) {
  if (df.containsElementNamed(name.c_str())) {
    Rcpp::NumericVector col = df[name];
    // Check for NAs and replace if necessary
    bool has_na = false;
    for (int i = 0; i < col.size(); ++i) {
      if (Rcpp::NumericVector::is_na(col[i])) {
        has_na = true;
        break;
      }
    }
    if (has_na) {
      Rcpp::NumericVector res = Rcpp::clone(col);
      for (int i = 0; i < res.size(); ++i) {
        if (Rcpp::NumericVector::is_na(res[i])) {
          res[i] = default_val;
        }
      }
      return res;
    }
    return col;
  }
  return Rcpp::NumericVector(df.nrow(), default_val);
}

#endif
