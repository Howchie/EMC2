#include <RcppArmadillo.h>
#include "utility_functions.h"
#include "model_lnr.h"
#include "model_LBA.h"
#include "model_RDM.h"
#include "model_DDM.h"
#include "model_MRI.h"
#include "model_SS.h"
#include "prob_utils.h"
#include "trend.h"
#include "utils.h"
#include "numerical_integration.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h> // For GSL error handling
#include <cmath>
#include <string>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <iomanip>
#include <array>
#include <limits>
#include <cstdlib>
using namespace Rcpp;

// Numerically stable log(exp(a) - exp(b)) for a > b
inline double log_diff_exp(double a, double b) {
  if (a <= b) return R_NegInf; // Probability <= 0 or logic error
  if (b == R_NegInf) return a;
  return a + std::log1p(-std::exp(b - a));
}

// Forward declaration: stop-signal likelihood with deadline censoring (rt==Inf uses UC).
double c_log_likelihood_ss_cens_trunc(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::DataFrame& dadm,
    int n_trials,
    const Rcpp::LogicalVector& winner,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& isok,
    int n_lR,
    const std::string& ss_type);

// Forward declaration: logical rules (OR/AND) redundant target race likelihood.
double c_log_likelihood_logicalrules(
    Rcpp::NumericMatrix pars,
    Rcpp::DataFrame dadm,
    RacePdfFun model_dfun,
    RaceCdfFun model_pfun,
    const int n_trials,
    const Rcpp::IntegerVector expand,
    double min_ll,
    const Rcpp::LogicalVector ok_params,
    int n_acc,
    void* model_specific_context);

// [[Rcpp::export]]
Rcpp::NumericMatrix do_transform(Rcpp::NumericMatrix pars, Rcpp::List transform) {
  // Build the specs for these parameters
  std::vector<TransformSpec> specs = make_transform_specs(pars, transform);
  // Apply transformation in place and return
  return c_do_transform(pars, specs);
}



NumericMatrix c_map_p(NumericVector p_vector,
                      CharacterVector p_types,
                      List designs,
                      int n_trials,
                      DataFrame data,
                      List trend,
                      const std::vector<TransformSpec>& full_specs) {
  
  // Extract information about trends
  const bool has_trend = (trend.length() > 0);
  bool premap = false;
  bool pretransform = false;
  CharacterVector trend_names;
  if (has_trend) {
    trend_names = trend.names();
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "premap") premap = true;
      if (ph == "pretransform") pretransform = true;
    }
  }
  
  const int n_params = p_types.size();
  NumericMatrix pars(n_trials, n_params);
  colnames(pars) = p_types;
  
  // Prepare trend parameter columns when needed
  NumericMatrix trend_pars;
  LogicalVector trend_index(n_params, FALSE);
  CharacterVector trend_pnames;
  if (has_trend && (premap || pretransform)) {
    // Fill in trend columns first so that they can be used in premapped trend
    // This function also applies transformations to the trend parameters
    // to ensure real-lines support.
    // The pre-transform trends are also included here, they are used after map_p
    // but they need to be transformed already (before the other trend parameters)
    // are transformed.
    trend_pars = build_trend_columns_from_design(p_vector, p_types, designs, n_trials, trend, full_specs);
    trend_pnames = colnames(trend_pars);
    trend_index = contains_multiple(p_types, trend_pnames);
  }
  
  // Map non-trend parameters from designs, applying premap trends if requested
  for (int i = 0; i < n_params; i++) {
    if (trend_index[i] == TRUE) continue; // skip trend parameters here
    NumericMatrix cur_design = designs[i];
    CharacterVector cur_names = colnames(cur_design);
    for (int j = 0; j < cur_design.ncol(); j++) {
      String cur_name(cur_names[j]);
      NumericVector p_mult_design(n_trials, p_vector[cur_name]);
      if (has_trend && premap) {
        p_mult_design = apply_premap_trends(data, trend, trend_names, cur_name, p_mult_design, trend_pars, p_vector);
      }
      p_mult_design = p_mult_design * cur_design(_, j);
      LogicalVector bad = is_na(p_mult_design) | is_nan(p_mult_design);
      p_mult_design[bad] = 0;
      pars(_, i) = pars(_, i) + p_mult_design;
    }
  }
  
  // If using pretransform trends, copy the pre-transformed trend cols into pars by name
  if (has_trend && pretransform) {
    // Only fill columns for pretransform entries
    CharacterVector tf_names = collect_trend_param_names_phase(trend, "pretransform");
    NumericMatrix trend_pars_tf = (tf_names.size() > 0) ? submat_rcpp_col_by_names(trend_pars, tf_names) : NumericMatrix(n_trials, 0);
    fill_trend_columns_for_pretransform(pars, p_types, trend_pars_tf);
  }
  
  // If premap, trend parameter columns are not part of the final matrix
  if (has_trend && premap) {
    CharacterVector names_premap = collect_trend_param_names_phase(trend, "premap");
    if (names_premap.size() > 0) {
      pars = submat_rcpp_col(pars, !contains_multiple(p_types, names_premap));
    }
  }
  return pars;
}

NumericMatrix get_pars_matrix(NumericVector p_vector, NumericVector constants, const std::vector<PreTransformSpec>& p_specs,
                              CharacterVector p_types, List designs, int n_trials, DataFrame data, List trend,
                              const std::vector<TransformSpec>& full_specs){
  const bool has_trend = (trend.length() > 0);
  bool pretransform = false;
  bool posttransform = false;
  if (has_trend) {
    for (int i = 0; i < trend.size(); ++i) {
      List cur = trend[i];
      std::string ph = Rcpp::as<std::string>(cur["phase"]);
      if (ph == "pretransform") pretransform = true;
      if (ph == "posttransform") posttransform = true;
    }
  }
  NumericVector p_vector_updtd(clone(p_vector));
  p_vector_updtd = c_do_pre_transform(p_vector_updtd, p_specs);
  p_vector_updtd = c_add_vectors(p_vector_updtd, constants);
  NumericMatrix pars = c_map_p(p_vector_updtd, p_types, designs, n_trials, data, trend, full_specs);
  // // Check if pretransform trend applies
  if(pretransform){
    pars = prep_trend_phase(data, trend, pars, "pretransform");
  }
  std::vector<TransformSpec> t_specs = make_transform_specs_from_full(pars, p_types, full_specs);
  pars = c_do_transform(pars, t_specs);
  // Check if posttransform trend applies
  if(posttransform){
    // Build trend parameter columns once (transformed) and pass override
    NumericMatrix trend_pars_all = build_trend_columns_from_design(p_vector_updtd, p_types, designs, n_trials, trend, full_specs);
    CharacterVector names_post = collect_trend_param_names_phase(trend, "posttransform");
    NumericMatrix trend_pars_post = (names_post.size() > 0) ? submat_rcpp_col_by_names(trend_pars_all, names_post) : NumericMatrix(n_trials, 0);
    pars = prep_trend_phase_with_pars(data, trend, pars, "posttransform", trend_pars_post);
  }
  // ok is calculated afterwards and Ttransform applied in the function
  return(pars);
}

double c_log_likelihood_DDM(NumericMatrix pars, DataFrame data,
                            const int n_trials, IntegerVector expand,
                            double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector rts = data["rt"];
  IntegerVector R = data["R"];
  NumericVector lls(n_trials);
  NumericVector lls_exp(n_out);
  lls = d_DDM_Wien(rts, R, pars, is_ok);
  lls_exp = c_expand(lls, expand); // decompress
  // lls_exp = lls;
  lls_exp[is_na(lls_exp)] = min_ll;
  lls_exp[is_infinite(lls_exp)] = min_ll;
  lls_exp[lls_exp < min_ll] = min_ll;
  return(sum(lls_exp));
}

// ZH - this has been replaced by the more robust cens/log function
double c_log_likelihood_race(NumericMatrix pars, DataFrame data,
                             NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                             const int n_trials, LogicalVector winner, IntegerVector expand,
                             double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  NumericVector lR = data["lR"];
  NumericVector lds_exp(n_out);
  const int n_lR = unique(lR).length();
  if(sum(contains(data.names(), "RACE")) == 1){
    NumericVector NACC = data["RACE"];
    CharacterVector vals_NACC = NACC.attr("levels");
    for(int x = 0; x < pars.nrow(); x++){
      // subtract 1 because R is 1 coded
      if(lR[x] > atoi(vals_NACC[NACC[x]-1])){
        pars(x,0) = NA_REAL;
      }
    }
  }
  
  //    Rcpp::Rcout << "oldC" << std::endl;
  
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok)); //first for compressed
  lds[winner] = win;
  if(n_lR > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds)] = min_ll;
  
  if(n_lR > 1){
    NumericVector ll_out = lds[winner];
    NumericVector lds_los = lds[!winner];
    if(n_lR == 2){
      ll_out = ll_out + lds_los;
    } else{
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_lR -1), (z+1) * (n_lR -1) -1)]);
      }
    }
    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    ll_out = c_expand(ll_out, expand); // decompress
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    lds_exp = c_expand(lds, expand); // decompress
    return(sum(lds_exp));
  }
}

double first_value_or_default_df(DataFrame df, std::string col, double fallback) {
  if (!df.containsElementNamed(col.c_str())) return fallback;
  NumericVector v = df[col];
  if (v.size() == 0) return fallback;
  return v[0];
}

double c_log_likelihood_race_missing(NumericMatrix pars, DataFrame data,
                                     NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                                     NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
                                     const int n_trials, LogicalVector winner, IntegerVector expand,
                                     double min_ll, LogicalVector is_ok){
  const int n_out = expand.length();
  NumericVector lds(n_trials);
  NumericVector rts = data["rt"];
  // Keep as integer codes (avoids factor -> character coercion).
  IntegerVector R = data["R"];
  NumericVector lR = data["lR"];
  NumericVector pCont = pars(_ , pars.ncol() - 1);
  NumericVector lds_exp(n_out);
  int n_acc = unique(lR).length();
  if(sum(contains(data.names(), "NACC")) == 1){
    NumericVector NACC = data["NACC"];
    for(int x = 0; x < pars.nrow(); x++){
      if(lR[x] > NACC[x]){
        pars(x,0) = NA_REAL;
      }
    }
  }
  NumericVector win = log(dfun(rts, pars, winner, exp(min_ll), is_ok)); //first for compressed
  lds[winner] = win;
  if(n_acc > 1){
    NumericVector loss = log(1- pfun(rts, pars, !winner, exp(min_ll), is_ok)); //cdfs
    loss[is_na(loss)] = min_ll;
    loss[loss == log(1 - exp(min_ll))] = min_ll;
    lds[!winner] = loss;
  }
  lds[is_na(lds) | (winner & is_infinite(rts))] = min_ll;
  lds[(!winner) & (is_infinite(rts) | is_na(rts))] = 0;
  
  double LT = first_value_or_default_df(data,"LT",0);
  double UT = first_value_or_default_df(data,"UT",R_PosInf);
  bool dotrunc = !(LT == 0 && UT == R_PosInf);
  double LC = first_value_or_default_df(data,"LC",0);
  double UC = first_value_or_default_df(data,"UC",R_PosInf);
  
  if (!(LT == 0 && LC == 0 && !R_finite(UT) && !R_finite(UC))) {
    
    //  Environment global = Environment::global_env();
    //  global["LC"]  = LC;
    //  global["UC"]  = UC;
    //  global["LT"]  = LT;
    //  global["UT"]  = UT;
    
    //For single trial integration
    LogicalVector ok1(n_acc,1);
    
    // Response known
    // Fast
    LogicalVector neginf = rts == R_NegInf; // also sets NA to FALSE
    LogicalVector nortfast = neginf & !is_na(R) & is_ok;
    
    if (is_true(any(nortfast))) {
      
      //    Rcpp::Rcout << "nortfast" << std::endl;
      
      NumericMatrix mparsfast(sum(nortfast),pars.ncol());
      for (int i = 0, j = 0; i < nortfast.length(); i++) {
        if (nortfast[i]) {
          mparsfast(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnerfastvec = winner[nortfast];
      LogicalMatrix winnerfast(n_acc, winnerfastvec.length() / n_acc, winnerfastvec.begin());
      
      LogicalVector tofixfast = (winner & nortfast);
      NumericVector ldstofixfast(sum(tofixfast));
      
      for (int i = 0; i < sum(tofixfast); i++) {
        NumericMatrix pifast(n_acc, mparsfast.ncol());
        if (n_acc == 1) {
          pifast(0,_) = mparsfast(i,_);
        } else {
          for (int j = 0; j < n_acc; j++) {
            pifast(j,_) = mparsfast(i * n_acc + j,_);
          }
        }
        NumericVector tmp = f_integrate(pifast, winnerfast(_,i), dfun, pfun, exp(min_ll), ok1, LT, LC);
        ldstofixfast[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
      }
      lds[tofixfast] = ldstofixfast;
    }
    
    // Slow
    LogicalVector posinf = rts == R_PosInf; // also sets NA to FALSE
    LogicalVector nortslow = posinf & !is_na(R) & is_ok;
    
    if (is_true(any(nortslow))) {
      
      //    Rcpp::Rcout << "nortslow" << std::endl;
      
      
      NumericMatrix mparsslow(sum(nortslow),pars.ncol());
      for (int i = 0, j = 0; i < nortslow.length(); i++) {
        if (nortslow[i]) {
          mparsslow(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnerslowvec = winner[nortslow];
      LogicalMatrix winnerslow(n_acc, winnerslowvec.length() / n_acc, winnerslowvec.begin());
      
      LogicalVector tofixslow = (winner & nortslow);
      NumericVector ldstofixslow(sum(tofixslow));
      
      for (int i = 0; i < sum(tofixslow); i++) {
        NumericMatrix pislow(n_acc, mparsslow.ncol());
        if (n_acc == 1) {
          pislow(0,_) = mparsslow(i,_);
        } else {
          for (int j = 0; j < n_acc; j++) {
            pislow(j,_) = mparsslow(i * n_acc + j,_);
          }
        }
        NumericVector tmp = f_integrate_slow(pislow, winnerslow(_,i), dfun, pfun, exp(min_ll), ok1, UC, UT);
        ldstofixslow[i] = std::log(std::max(0.0, std::min(tmp[0], 1.0)));
      }
      lds[tofixslow] = ldstofixslow;
    }
    
    // No direction
    LogicalVector nortno = is_na(rts) & !is_na(R) & is_ok;
    
    if (is_true(any(nortno))) {
      
      //    Rcpp::Rcout << "nortno" << std::endl;
      
      NumericMatrix mparsno(sum(nortno), pars.ncol());
      
      for (int i = 0, j = 0; i < nortno.length(); i++) {
        if (nortno[i]) {
          mparsno(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnernovec = winner[nortno];
      LogicalMatrix winnerno(n_acc, winnernovec.length() / n_acc, winnernovec.begin());
      
      LogicalVector tofixno = (winner & nortno);
      NumericVector ldstofixno(sum(tofixno));
      
      for (int i = 0; i < sum(tofixno); i++) {
        NumericMatrix pino(n_acc, mparsno.ncol());
        if (n_acc  == 1) {
          pino(0,_) = mparsno(i,_);
        } else {
          for (int j = 0; j < n_acc; j++) {
            pino(j,_) = mparsno(i * n_acc + j,_);
          }
        }
        NumericVector tmpslow = f_integrate_slow(pino, winnerno(_,i), dfun, pfun, exp(min_ll), ok1, UC, UT);
        NumericVector tmpfast = f_integrate(pino, winnerno(_,i), dfun, pfun, exp(min_ll), ok1, LT, LC);
        double tmp = tmpslow[0] + tmpfast[0];
        
        ldstofixno[i] = std::log(std::max(0.0, std::min(tmp, 1.0)));
      }
      lds[tofixno] = ldstofixno;
    }
    
    // Response unknown
    // Fast
    LogicalVector nortfastu = (rts == R_NegInf) & is_na(R)  & is_ok;
    LogicalVector tofixfast(winner.length());
    
    if (is_true(any(nortfastu))) {
      
      //    Rcpp::Rcout << "nortfastu" << std::endl;
      
      
      NumericMatrix mpars(sum(nortfastu), pars.ncol());
      
      for (int i = 0, j = 0; i < nortfastu.length(); i++) {
        if (nortfastu[i]) {
          mpars(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnerfastuvec = winner[nortfastu];
      LogicalMatrix winnerfastu(n_acc, winnerfastuvec.length() / n_acc, winnerfastuvec.begin());
      
      tofixfast = (winner & nortfastu);
      NumericVector ldstofixfast(sum(tofixfast));
      
      for (int i = 0; i < sum(tofixfast); i++) {
        NumericMatrix pi(n_acc, mpars.ncol());
        
        for (int j = 0; j < n_acc; j++) {
          pi(j,_) = mpars(i * n_acc + j,_);
        }
        
        LogicalVector idx(n_acc, 0);
        idx[0] = 1;
        
        NumericVector pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), ok1, LT, LC);
        double p;
        if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
          p = NA_REAL;
        } else{
          p = std::max(0.0 ,std::min(pc[0],1.0));
        }
        
        double cf;
        if (p != 0 && !(LT==0 && UT==R_PosInf)) {
          cf = pr_pt(pi, pfun, LT, UT);
        } else {
          cf = 1;
        }
        
        if (!traits::is_na<REALSXP>(cf)) {
          p *= cf;
        }
        
        if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
          for (int j = 1; j < n_acc; j++) {
            idx.fill(0);
            idx[j] = 1;
            pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), ok1, LT, LC);
            if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
              p = NA_REAL;
              break;
            }
            if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
              cf = pr_pt(pi, pfun, LT, UT);
            } else{
              cf = 1;
            }
            if (!traits::is_na<REALSXP>(cf)) {
              p += pc[0] * cf;
            }
          }
        }
        double lp = std::log(p);
        if (!traits::is_na<REALSXP>(lp)) {
          ldstofixfast[i] = lp;
        } else{
          ldstofixfast[i] = R_NegInf;
        }
      }
      lds[tofixfast] = ldstofixfast;
    }
    
    // Slow
    LogicalVector nortslowu = (rts == R_PosInf) & is_na(R) & is_ok;
    LogicalVector tofixslow(winner.length());
    if (is_true(any(nortslowu))) {
      
      //    Rcpp::Rcout << "nortslowu" << std::endl;
      
      NumericMatrix mpars(sum(nortslowu), pars.ncol());
      for (int i = 0, j = 0; i < nortslowu.length(); i++) {
        if (nortslowu[i]) {
          mpars(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnerslowuvec = winner[nortslowu];
      LogicalMatrix winnerslowu(n_acc, winnerslowuvec.length() / n_acc, winnerslowuvec.begin());
      tofixslow = (winner & nortslowu);
      NumericVector ldstofixslow(sum(tofixslow));
      for (int i = 0; i < sum(tofixslow); i++) {
        NumericMatrix pi(n_acc, mpars.ncol());
        for (int j = 0; j < n_acc; j++) {
          pi(j,_) = mpars(i * n_acc + j,_);
        }
        
        LogicalVector idx(n_acc);
        idx[0] = 1;
        NumericVector pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), ok1, UC, UT);
        double p;
        if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
          p = NA_REAL;
        } else{
          p = std::max(0.0,std::min(pc[0],1.0));
        }
        
        double cf;
        if (p != 0 && !(LT==0 && UT==R_PosInf)) {
          cf = pr_pt(pi, pfun, LT, UT);
        } else {
          cf = 1;
        }
        
        if (!traits::is_na<REALSXP>(cf)) {
          p *= cf;
        }
        
        if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
          for (int j = 1; j < n_acc; j++) {
            idx.fill(0);
            idx[j] = 1;
            pc = f_integrate(pi, idx, dfun, pfun, exp(min_ll), ok1, UC, UT);
            if (pc[2] != 0 || traits::is_nan<REALSXP>(pc[0])) {
              p = NA_REAL;
              break;
            }
            if (pc[0] != 0.0 && !(LT == 0 && UT == R_PosInf)) {
              cf = pr_pt(pi, pfun, LT, UT);
            } else{
              cf = 1;
            }
            if (!traits::is_na<REALSXP>(cf)) {
              p += pc[0] * cf;
            }
          }
        }
        double lp = std::log(p);
        if (!traits::is_na<REALSXP>(lp)) {
          ldstofixslow[i] = lp;
        } else{
          ldstofixslow[i] = R_NegInf;
        }
      }
      lds[tofixslow] = ldstofixslow;
    }
    
    // No direction
    LogicalVector nortnou = is_na(rts) & is_na(R) & is_ok;
    
    if (is_true(any(nortnou))) {
      
      //    Rcpp::Rcout << "nortnou" << std::endl;
      
      NumericMatrix mpars(sum(nortnou), pars.ncol());
      
      for (int i = 0, j = 0; i < nortnou.length(); i++) {
        if (nortnou[i]) {
          mpars(j,_) = pars(i,_);
          j++;
        }
      }
      
      LogicalVector winnernouvec = winner[nortnou];
      LogicalMatrix winnernou(n_acc, winnernouvec.length() / n_acc, winnernouvec.begin());
      
      LogicalVector tofix = (winner & nortnou);
      NumericVector ldstofix(sum(tofix));
      
      for (int i = 0; i < sum(tofix); i++) {
        NumericMatrix pi(n_acc, mpars.ncol());
        
        for (int j = 0; j < n_acc; j++) {
          pi(j,_) = mpars(i * n_acc + j,_);
        }
        
        LogicalVector idx(n_acc);
        idx[0] = 1;
        
        double pc = pLU(pi, idx, dfun, pfun, exp(min_ll), ok1, LT, LC, UC, UT);
        double p;
        double cf;
        if (traits::is_na<REALSXP>(pc)) {
          p = NA_REAL;
        } else{
          if (pc != 0.0 && !(LT == 0 && UT == R_PosInf)) {
            cf = pr_pt(pi, pfun, LT, UT);
          } else{
            cf = 1;
          }
          
          if (!traits::is_na<REALSXP>(cf)) {
            p = pc*cf;
          } else {
            p = NA_REAL;
          }
          
          if (!traits::is_na<REALSXP>(p) && n_acc > 1) {
            for (int j = 1; j < n_acc; j++) {
              idx.fill(0);
              idx[j] = 1;
              
              pc = pLU(pi, idx, dfun, pfun, exp(min_ll), ok1, LT, LC, UC, UT);
              if (traits::is_na<REALSXP>(pc)) {
                p = NA_REAL;
                break;
              }
              if (pc != 0 && !(LT == 0 && UT == R_PosInf)) {
                cf = pr_pt(pi, pfun, LT, UT);
              } else {
                cf = 1;
              }
              
              if (traits::is_na<REALSXP>(cf)) {
                p = NA_REAL;
                break;
              } else{
                p += pc * cf;
              }
            }
          }
        }
        double lp = std::log(p);
        if (!traits::is_na<REALSXP>(lp)) {
          ldstofix[i] = lp;
        } else{
          ldstofix[i] = R_NegInf;
        }
      }
      lds[tofix] = ldstofix;
    }
    
    // Truncation where not censored or censored and response known
    LogicalVector unique_nort = data.attr("unique_nort");
    NumericVector uniquewinlike = lds[unique_nort & winner];
    LogicalVector ok = is_finite(uniquewinlike);
    IntegerVector uniquewinresp = R[unique_nort & winner];
    LogicalVector alreadyfixed = is_na(uniquewinresp);
    ok = ok & !alreadyfixed & is_ok;
    
    if (dotrunc & is_true(any(ok))) {
      
      //    Rcpp::Rcout << "dotrunc" << std::endl;
      
      IntegerVector expand_nort = data.attr("expand_nort");
      NumericMatrix tpars(sum(unique_nort),pars.ncol());
      for (int i=0, j=0; i < unique_nort.length(); i++) {
        if (unique_nort[i]) {
          tpars(j,_) = pars(i,_);
          j++;
        }
      }
      LogicalVector winnertruncvec = winner[unique_nort];
      LogicalMatrix winnertrunc(n_acc, winnertruncvec.length() / n_acc, winnertruncvec.begin());
      NumericVector cf = rep(NA_REAL, ok.length());
      
      for (int i = 0; i < ok.length(); i++) {
        if (ok[i]) {
          NumericMatrix pi(n_acc, tpars.ncol());
          
          for (int j = 0; j < n_acc; j++) {
            pi(j,_) = tpars(i * n_acc + j , _ );
          }
          
          cf[i] = pr_pt(pi, pfun, LT, UT);
        }
      }
      NumericVector cf_log = rep_each(log(cf), n_acc);
      NumericVector cf_exp = c_expand(cf_log, expand_nort);
      LogicalVector fix = winner & !is_na(cf_exp) & is_finite(cf_exp);
      if (is_true(any(fix))) {
        lds[fix] = lds[fix] + cf_exp[fix];
      }
      LogicalVector badfix = winner & (is_na(cf_exp) | is_infinite(cf_exp));
      if (all(!is_na(tofixfast))) {
        badfix = badfix & !tofixfast;
      }
      if (all(!is_na(tofixslow))) {
        badfix = badfix & !tofixslow;
      }
      if (is_true(any(badfix))) {
        lds[badfix] = R_NegInf;
      }
    }
    
  }
  
  // Non-process (contaminant) miss.
  LogicalVector isPCont = (pCont > 0);
  LogicalVector okwin = winner & is_ok;
  if (is_true(any(isPCont & okwin))) {
    
    //    Rcpp::Rcout << "isPCont" << std::endl;
    
    NumericVector p = lds[okwin];
    p = exp(p);
    NumericVector pc = pCont[okwin];
    IntegerVector Rwin = R[okwin];
    LogicalVector isMiss = is_na(Rwin);
    for (int i = 0; i < p.length(); i++) {
      if (isMiss[i]) {
        p[i] = pc[i] + (1  - pc[i]) * p[i];
      } else {
        p[i] = (1 - pc[i]) * p[i];
      }
    }
    NumericVector ldswinner = log(p);
    lds[okwin] = ldswinner;
  }
  
  
  
  if(n_acc > 1){
    NumericVector ll_out = lds[winner];
    NumericVector lds_los = lds[!winner];
    if(n_acc == 2){
      ll_out = ll_out + lds_los;
    } else{
      for(int z = 0; z < ll_out.length(); z++){
        ll_out[z] = ll_out[z] + sum(lds_los[seq( z * (n_acc -1), (z+1) * (n_acc -1) -1)]);
      }
    }
    ll_out[is_na(ll_out)] = min_ll;
    ll_out[is_infinite(ll_out)] = min_ll;
    ll_out[ll_out < min_ll] = min_ll;
    ll_out = c_expand(ll_out, expand); // decompress
    return(sum(ll_out));
  } else{
    lds_exp[is_na(lds_exp)] = min_ll;
    lds_exp[is_infinite(lds_exp)] = min_ll;
    lds_exp[lds_exp < min_ll] = min_ll;
    lds_exp = c_expand(lds, expand); // decompress
    return(sum(lds_exp));
  }
  
}

// [[Rcpp::export]]
NumericVector calc_ll(NumericMatrix p_matrix, DataFrame data, NumericVector constants,
                      List designs, String type, List bounds, List transforms, List pretransforms,
                      CharacterVector p_types, double min_ll, List trend){
  const int n_particles = p_matrix.nrow();
  const int n_trials = data.nrow();
  NumericVector lls(n_particles);
  NumericVector p_vector(p_matrix.ncol());
  CharacterVector p_names = colnames(p_matrix);
  NumericMatrix pars(n_trials, p_types.length());
  p_vector.names() = p_names;
  LogicalVector is_ok(n_trials);
  
  // Once (outside the main loop over particles):
  NumericMatrix minmax = bounds["minmax"];
  CharacterVector mm_names = colnames(minmax);
  std::vector<PreTransformSpec> p_specs;
  std::vector<BoundSpec> bound_specs;
  std::vector<TransformSpec> full_t_specs; // precomputed transform specs for p_types
  
  if(type == "DDM"){
    IntegerVector expand = data.attr("expand");
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      lls[i] = c_log_likelihood_DDM(pars, data, n_trials, expand, min_ll, is_ok);
    }
  } else if(type == "MRI" || type == "MRI_AR1"){
    int n_pars = p_types.length();
    NumericVector y = extract_y(data);
    for(int i = 0; i < n_particles; i++){
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      // Precompute specs
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      if(type == "MRI"){
        lls[i] = c_log_likelihood_MRI(pars, y, is_ok, n_trials, n_pars, min_ll);
      } else{
        lls[i] = c_log_likelihood_MRI_white(pars, y, is_ok, n_trials, n_pars, min_ll);
      }
    }
  } else if(type == "MLBA" || type == "MRDM" || type == "MLNR" || type == "OLBA" || type == "ORDM" || type == "OLNR") {
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Love me some good old ugly but fast c++ pointers
    NumericVector (*dfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector);
    if(type == "MLBA" || type == "OLBA"){
      dfun = dlba_c1;
      pfun = plba_c1;
    } else if(type == "MRDM"|| type == "ORDM"){
      dfun = drdm_c;
      pfun = prdm_c;
    } else {
      dfun = dlnr_c;
      pfun = plnr_c;
    }
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    for (int i = 0; i < n_particles; ++i) {
      p_vector = p_matrix(i, _);
      if(i == 0){
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);
      if(type == "MLBA" || type == "MRDM" || type == "MLNR") {
        lls[i] = c_log_likelihood_race_missing(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
      } else {
        lls[i] = c_log_likelihood_race(pars, data, dfun, pfun, n_trials, winner, expand, min_ll, is_ok);
      }
    }
  } else {
    // ZH I have adjusted a lot here. I made new race model pointers that help branch things like intrinsic omissions LBA, and the new GNG models
    // The "1" versions handle scalar likelihoods that do not call Rcpp functions,which speeds up the GSL integrations used in censoring, truncation, and GNG
    IntegerVector expand = data.attr("expand");
    LogicalVector winner = data["winner"];
    // Standardize incoming model type string from R
    // For race models: e.g. "LBA", "LBAIO", "LBAGNG", "RDM", "RDMGNG", "LNR", "LNRGNG"
    // For stop-signal models (C++): e.g. "SSexG", "SShybrid"
    std::string type_std = Rcpp::as<std::string>(Rcpp::wrap(type));
    const bool is_ss = (type_std.find("SS") != std::string::npos);
    const bool is_logicalrules = (type_std.find("LogicalRules") != std::string::npos);
    RacePdfFun model_dfun_ptr = nullptr;
    RaceCdfFun model_pfun_ptr = nullptr;
    RacePdf1Fun pdf1_ptr = nullptr;
    RaceCdf1Fun cdf1_ptr = nullptr;
    ContextForRaceModels current_model_ctx;
    current_model_ctx.min_lik_for_pdf = std::exp(min_ll);
    current_model_ctx.use_posdrift = true; // Default: LBA uses posdrift, RDM/LNR ignore this context field.
    current_model_ctx.gng = false; // Default: Do not use go/no-go likelihoods by default
    if (!is_ss) {
      // Determine adapter functions and specific LBA posdrift setting based on type_std
      if (type_std.find("LBA") != std::string::npos) {
        model_dfun_ptr = &lba_dfun_adapter;
        model_pfun_ptr = &lba_pfun_adapter;
        pdf1_ptr = &dlba_scalar;
        cdf1_ptr = &plba_scalar;
        // Check for the 'IO' (Implicit Omissions / no posdrift) flag in the original type_std
        if (type_std.find("IO") != std::string::npos) {
          current_model_ctx.use_posdrift = false;
        }
      } else if (type_std == "RDM") {
        model_dfun_ptr = &rdm_dfun_adapter;
        model_pfun_ptr = &rdm_pfun_adapter;
        pdf1_ptr = &drdm_scalar;
        cdf1_ptr = &prdm_scalar;
      } else if (type_std.find("LNR") != std::string::npos) {
        model_dfun_ptr = &lnr_dfun_adapter;
        model_pfun_ptr = &lnr_pfun_adapter;
        pdf1_ptr = &dlnr_scalar;
        cdf1_ptr = &plnr_scalar;
      } else {
        Rcpp::stop("Unsupported race model type string in calc_ll: " + type_std);
      }
      if (type_std.find("GNG") != std::string::npos) {
        current_model_ctx.gng = true;
      }
    }
    // Expose scalar PDF/CDF + parameter count in the context for any composite likelihoods.
    current_model_ctx.pdf1 = pdf1_ptr;
    current_model_ctx.cdf1 = cdf1_ptr;
    current_model_ctx.n_par = p_types.length();
    NumericVector lR = data["lR"];
    int n_lR = unique(lR).length();
    
    constexpr int kLlCacheVersion = 1;
    bool cache_ok = false;
    if (data.hasAttribute("emc2_ll_cache_version")) {
      int ver = Rcpp::as<int>(data.attr("emc2_ll_cache_version"));
      cache_ok = (ver == kLlCacheVersion);
    }

    const bool has_RACE_col = data.containsElementNamed("RACE");

    bool all_finite_trials = true;
    bool have_all_finite_attr = false;
    if (cache_ok && data.hasAttribute("emc2_all_finite_trials")) {
      Rcpp::LogicalVector v = data.attr("emc2_all_finite_trials");
      if (v.size() == 1 && v[0] != NA_LOGICAL) {
        all_finite_trials = v[0];
        have_all_finite_attr = true;
      }
    }

    const bool have_race_attrs = (!has_RACE_col) ||
      (data.hasAttribute("RACE_nacc_by_row") && data.hasAttribute("RACE_mask"));
    const bool have_finite_attrs = data.hasAttribute("finite_rt_mask") &&
      data.hasAttribute("finite_rt_unique_trial_indices") &&
      data.hasAttribute("other_unique_trial_indices");

    const bool need_prepare =
      !(cache_ok && have_all_finite_attr && have_race_attrs && (all_finite_trials || have_finite_attrs));

    if (need_prepare) {
      if (has_RACE_col && !have_race_attrs) {
        const Rcpp::IntegerVector lR_dadm = data["lR"];
        const Rcpp::IntegerVector race_idx = data["RACE"];
        const Rcpp::CharacterVector race_levels = race_idx.attr("levels");
        std::vector<int> nacc_by_level;
        nacc_by_level.reserve(static_cast<size_t>(race_levels.size()));
        for (int i = 0; i < race_levels.size(); ++i) {
          nacc_by_level.push_back(std::stoi(Rcpp::as<std::string>(race_levels[i])));
        }
        Rcpp::IntegerVector race_nacc_by_row(n_trials, n_lR);
        Rcpp::LogicalVector race_mask(n_trials, true);
        for (int row = 0; row < n_trials; ++row) {
          if (race_idx[row] == NA_INTEGER) continue;
          const int n_lR_this_trial = nacc_by_level[static_cast<size_t>(race_idx[row] - 1)];
          race_nacc_by_row[row] = n_lR_this_trial;
          if (lR_dadm[row] > n_lR_this_trial) race_mask[row] = false;
        }
        data.attr("RACE_nacc_by_row") = race_nacc_by_row;
        data.attr("RACE_mask") = race_mask;
      }

      if (!have_all_finite_attr) {
        // Data-only: determine whether every unique trial has a finite RT, known response, and is within [LT, UT].
        // Computing this once here avoids repeating an O(n_unique_trials) scan per particle.
        all_finite_trials = true;
        if (n_trials > 0) {
          if (n_lR <= 0 || (n_trials % n_lR) != 0) {
            all_finite_trials = false;
          } else {
            const int n_unique_trials = n_trials / n_lR;
            const Rcpp::NumericVector rts_dadm = data["rt"];
            const Rcpp::IntegerVector R_idxs_dadm = data["R"];
            for (int j = 0; j < n_unique_trials; ++j) {
              const int start_row_idx = j * n_lR;
              const double rt_j = rts_dadm[start_row_idx];
              const int R_j_idx = R_idxs_dadm[start_row_idx];
              if (!(R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER)) {
                all_finite_trials = false;
                break;
              }
            }
          }
        }
        data.attr("emc2_all_finite_trials") = all_finite_trials;
      }

      if (!all_finite_trials && n_trials > 0 && n_lR > 0 && (n_trials % n_lR) == 0) {
        if (!have_finite_attrs) {
          const int n_unique_trials = n_trials / n_lR;
          Rcpp::LogicalVector finite_rt_mask(n_trials, false);
          std::vector<int> finite_rt_unique_trial_indices;
          std::vector<int> other_unique_trial_indices;
          finite_rt_unique_trial_indices.reserve(n_unique_trials);
          other_unique_trial_indices.reserve(n_unique_trials);
          const Rcpp::NumericVector rts_dadm = data["rt"];
          const Rcpp::IntegerVector R_idxs_dadm = data["R"];
          Rcpp::IntegerVector race_nacc_by_row;
          if (has_RACE_col && data.hasAttribute("RACE_nacc_by_row")) {
            race_nacc_by_row = data.attr("RACE_nacc_by_row");
          }
          for (int j = 0; j < n_unique_trials; ++j) {
            const int start_row_idx = j * n_lR;
            const double rt_j = rts_dadm[start_row_idx];
            const int R_j_idx = R_idxs_dadm[start_row_idx];
            const int n_lR_j = (has_RACE_col && race_nacc_by_row.size() == n_trials)
              ? race_nacc_by_row[start_row_idx]
              : n_lR;
            if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx != NA_INTEGER) {
              finite_rt_unique_trial_indices.push_back(j);
              for (int k_acc = 0; k_acc < n_lR_j; ++k_acc) {
                finite_rt_mask[start_row_idx + k_acc] = true;
              }
            } else {
              other_unique_trial_indices.push_back(j);
            }
          }
          data.attr("finite_rt_mask") = finite_rt_mask;
          data.attr("finite_rt_unique_trial_indices") = Rcpp::IntegerVector(finite_rt_unique_trial_indices.begin(),
                                                                            finite_rt_unique_trial_indices.end());
          data.attr("other_unique_trial_indices") = Rcpp::IntegerVector(other_unique_trial_indices.begin(),
                                                                        other_unique_trial_indices.end());
        }
      }

      data.attr("emc2_ll_cache_version") = kLlCacheVersion;
    }
    
    for (int i = 0; i < n_particles; ++i) {
      p_vector = p_matrix(i, Rcpp::_);
      if (i == 0) {
        p_specs = make_pretransform_specs(p_vector, pretransforms);
        // Precompute transform specs for all p_types using a one-time dummy
        NumericMatrix dummy(1, p_types.size());
        colnames(dummy) = p_types;
        full_t_specs = make_transform_specs(dummy, transforms);
      }
      pars = get_pars_matrix(p_vector, constants, p_specs, p_types, designs, n_trials, data, trend, full_t_specs);
      if (i == 0) {                            // first particle only, just to get colnames
        bound_specs = make_bound_specs(minmax,mm_names,pars,bounds);
      }
      is_ok = c_do_bound(pars, bound_specs);
      is_ok = lr_all(is_ok, n_lR);

      // Optional debugging: print which integrator path is being used (qng vs qag fallback)
      // for LogicalRules likelihood. To enable:
      //   Sys.setenv(EMC2_LOGICALRULES_LOG_INTEGRATION="1")
      // This prints only for the first particle to avoid flooding output.
      bool log_integrator = false;
      if (const char* v = std::getenv("EMC2_LOGICALRULES_LOG_INTEGRATION")) {
        log_integrator = (v[0] != '\0' && v[0] != '0');
      }
      current_model_ctx.log_out = (log_integrator && (i == 0));
      if (is_ss) {
        lls[i] = c_log_likelihood_ss_cens_trunc(pars, data,
                                                n_trials, winner, expand,
                                                min_ll, is_ok, n_lR, type_std);
      } else if (is_logicalrules){
        lls[i] = c_log_likelihood_logicalrules(pars, data,
                                              model_dfun_ptr, model_pfun_ptr,
                                              n_trials, expand, min_ll, is_ok,
                                              n_lR, &current_model_ctx);
      }
      else {
        lls[i] = c_log_likelihood_race_cens_trunc(pars, data,
                                                  model_dfun_ptr, model_pfun_ptr,
                                                  pdf1_ptr, cdf1_ptr,
                                                  n_trials,
                                                  winner, expand, min_ll, is_ok, n_lR,
                                                  &current_model_ctx,
                                                  all_finite_trials);
      }
    }
  }
  return(lls);
}

// gsl adapter for integrals - uses scalar, Rcpp-independent functions for speed
double gsl_f_race_scalar(double t, void* p) {
  auto* P = static_cast<gsl_race_params_scalar*>(p);
  if (t <= 0.0) return 0.0;
  const int w = P->winner_idx0;
  if (w < 0 || w >= P->n_lR) return 0.0;
  if (!P->isok[w]) return 0.0;
  
  const double* par_w = P->pars + static_cast<size_t>(w) * P->n_par;
  double out = P->pdf1(t, par_w, P->ctx);
  if (!(out > 0.0) || !std::isfinite(out)) return 0.0;
  
  for (int j = 0; j < P->n_lR; ++j) {
    if (j == w) continue;
    if (!P->isok[j]) return 0.0;
    const double* par_j = P->pars + static_cast<size_t>(j) * P->n_par;
    const double cdf = P->cdf1(t, par_j, P->ctx);
    if (!(cdf > 0.0)) continue;
    if (cdf >= 1.0 || !std::isfinite(cdf)) return 0.0;
    out *= (1.0 - cdf);
    if (!(out > 0.0) || !std::isfinite(out)) return 0.0;
  }
  return out;
}


bool row_is_finite(const Rcpp::NumericMatrix& mat, int row) {
  // Convenience helper for sanity checks on parameter rows.
  for (int j = 0; j < mat.ncol(); ++j) {
    if (!R_finite(mat(row, j))) return false;
  }
  return true;
}

inline double log_survivor_rowmajor(double t,
                                    const double* pars_rowmajor,
                                    const int* isok_int,
                                    int n_lR,
                                    int n_par,
                                    RaceCdf1Fun cdf1,
                                    void* ctx) {
  // Log survivor of the race at time t:
  //   log S(t) = sum_k log(1 - F_k(t))
  //
  // "rowmajor" here means per-trial parameters for the n_lR accumulators are
  // packed as a contiguous buffer:
  //   pars_rowmajor[k * n_par + c]
  //
  // This avoids repeatedly indexing into an Rcpp::NumericMatrix inside tight
  // loops and matches the representation needed for the scalar/GSL integrands.
  if (t == R_PosInf) return R_NegInf;
  double logS = 0.0;
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!std::isfinite(ll)) return R_NegInf;
    logS += ll;
  }
  return logS;
}

inline double log_min_density_rowmajor(double t,
                                       const double* pars_rowmajor,
                                       const int* isok_int,
                                       int n_lR,
                                       int n_par,
                                       RacePdf1Fun pdf1,
                                       RaceCdf1Fun cdf1,
                                       void* ctx,
                                       std::vector<double>& logS_k) {
  // Log density of the minimum (no known winner) at time t:
  //   f_min(t) = sum_k f_k(t) * prod_{j != k} (1 - F_j(t))
  //
  // Work in log space:
  //   log f_min(t) = logsumexp_k [ log f_k(t) + sum_{j != k} log(1 - F_j(t)) ].
  //
  // logS_k is a per-accumulator scratch buffer (reused across calls) holding
  // log(1 - F_k(t)) for the current t, to avoid reallocations in the "other
  // trial" loop.
  if (!(t > 0.0) || !std::isfinite(t)) return R_NegInf;
  double logS_all = 0.0;
  if (logS_k.size() < static_cast<size_t>(n_lR)) {
    logS_k.assign(static_cast<size_t>(n_lR), R_NegInf);
  } else {
    std::fill(logS_k.begin(), logS_k.begin() + n_lR, R_NegInf);
  }
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_int[k]) return R_NegInf;
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    double cdf = cdf1(t, par_k, ctx);
    cdf = clamp_cdf01_race(cdf);
    const double ll = safe_log1m_race(cdf);
    if (!std::isfinite(ll)) return R_NegInf;
    logS_k[static_cast<size_t>(k)] = ll;
    logS_all += ll;
  }
  double out = R_NegInf;
  for (int k = 0; k < n_lR; ++k) {
    const double* par_k = pars_rowmajor + static_cast<size_t>(k) * n_par;
    const double pdf = pdf1(t, par_k, ctx);
    if (!(pdf > 0.0) || !std::isfinite(pdf)) continue;
    const double term = std::log(pdf) + (logS_all - logS_k[static_cast<size_t>(k)]);
    out = log_sum_exp(out, term);
  }
  return out;
}

using GslWorkspacePtr = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;

inline gsl_integration_workspace* ensure_gsl_workspace(GslWorkspacePtr& ws, size_t n = 1000) {
  // Lazily allocate the GSL workspace: most trials never need numerical
  // integration (finite RT + no truncation/censoring), so avoid paying for it.
  if (!ws) {
    ws.reset(gsl_integration_workspace_alloc(n));
    if (!ws) Rcpp::stop("Failed to allocate GSL integration workspace.");
  }
  return ws.get();
}

double integrate_for_kth_winner_rowmajor_cpp(
    int k_winner_idx, // 1-based
    const double* pars_rowmajor,
    const int* isok_int,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w);

// Numerical integration helper using GSL
double integrate_for_kth_winner_cpp(
    int k_winner_idx, // 1-based
    const Rcpp::NumericMatrix& p_all_acc,
    LogicalVector isok,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w) {
  
  if (low >= upp && !(low == 0 && upp == R_PosInf)) return R_NegInf;
  if (k_winner_idx < 1 || k_winner_idx > (n_lR_j)) {
    Rcpp::Rcerr << "Warning: Invalid k_winner_idx in integrate_for_kth_winner_cpp: " << k_winner_idx << std::endl;
    return R_NegInf;
  }
  
  if (w == nullptr) {
    Rcpp::stop("integrate_for_kth_winner_cpp: GSL workspace is null.");
  }
  
  std::vector<double> pars_rowmajor(static_cast<size_t>(n_lR_j) * n_par);
  for (int r = 0; r < n_lR_j; ++r) {
    for (int c = 0; c < n_par; ++c) {
      pars_rowmajor[static_cast<size_t>(r) * n_par + c] = p_all_acc(r, c);
    }
  }
  
  std::vector<int> isok_int(n_lR_j, 0);
  for (int r = 0; r < n_lR_j; ++r) isok_int[r] = isok[r] ? 1 : 0;
  
  return integrate_for_kth_winner_rowmajor_cpp(k_winner_idx,
                                               pars_rowmajor.data(),
                                               isok_int.data(),
                                               low,
                                               upp,
                                               pdf1,
                                               cdf1,
                                               n_lR_j,
                                               n_par,
                                               epsilon,
                                               model_specific_context,
                                               w);
}

double integrate_for_kth_winner_rowmajor_cpp(
    int k_winner_idx, // 1-based
    const double* pars_rowmajor,
    const int* isok_int,
    double low,
    double upp,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    int n_lR_j,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w) {
  
  // Integrate the k-th winner density over an interval [low, upp] for a single
  // unique trial, using rowmajor buffers (raw pointers). This is the fast path
  // used by truncation/censoring normalisers and by go/no-go branches; it avoids
  // Rcpp object traffic inside the GSL callback.
  if (low >= upp && !(low == 0 && upp == R_PosInf)) return R_NegInf;
  if (k_winner_idx < 1 || k_winner_idx > n_lR_j) return R_NegInf;
  if (w == nullptr) Rcpp::stop("integrate_for_kth_winner_rowmajor_cpp: GSL workspace is null.");
  
  gsl_function F;
  gsl_race_params_scalar params_struct;
  params_struct.pars = pars_rowmajor;
  params_struct.n_lR = n_lR_j;
  params_struct.n_par = n_par;
  params_struct.winner_idx0 = k_winner_idx - 1;
  params_struct.isok = isok_int;
  params_struct.pdf1 = pdf1;
  params_struct.cdf1 = cdf1;
  params_struct.ctx = model_specific_context;
  
  F.function = &gsl_f_race_scalar;
  F.params = &params_struct;
  
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status;
  double result = 0.0;
  double error = 0.0;
  if (upp == R_PosInf) {
    if (low < 0) low = 0; // QAGIU requires a >= 0
    if (low >= R_PosInf) {
      result = 0.0;
      status = GSL_SUCCESS;
    } else {
      status = gsl_integration_qagiu(&F, low, 0, epsilon, 1000, w, &result, &error);
    }
  } else {
    status = gsl_integration_qags(&F, low, upp, 0, epsilon, 1000, w, &result, &error);
  }
  
  gsl_set_error_handler(old_handler);
  if (status != GSL_SUCCESS) return R_NegInf;
  if (!(result > 0.0) || !R_finite(result)) return R_NegInf;
  return std::log(result);
}

double get_trunc_normaliser_rowmajor_cpp(const double* pars_rowmajor,
                                         const int* isok_int,
                                         RacePdf1Fun pdf1,
                                         RaceCdf1Fun cdf1,
                                         double LT,
                                         double UT,
                                         int n_lR,
                                         int n_par,
                                         double epsilon,
                                         void* model_specific_context,
                                         GslWorkspacePtr& workspace) {
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  const double logS_LT = log_survivor_rowmajor(LT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  if (!R_FINITE(logS_LT)) return R_NegInf;
  if (UT == R_PosInf) return logS_LT;
  
  const double logS_UT = log_survivor_rowmajor(UT, pars_rowmajor, isok_int, n_lR, n_par, cdf1, model_specific_context);
  double logP = log_diff_exp(logS_LT, logS_UT);
  if (R_FINITE(logP) && logP > log_prob_eps) return logP;
  
  gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
  double log_total = R_NegInf;
  for (int k_win = 1; k_win <= n_lR; ++k_win) {
    const double log_k = integrate_for_kth_winner_rowmajor_cpp(k_win,
                                                               pars_rowmajor,
                                                               isok_int,
                                                               LT,
                                                               UT,
                                                               pdf1,
                                                               cdf1,
                                                               n_lR,
                                                               n_par,
                                                               epsilon,
                                                               model_specific_context,
                                                               w);
    log_total = log_sum_exp(log_total, log_k);
  }
  if (R_FINITE(log_total) && log_total > log_prob_eps) return log_total;
  return R_NegInf;
}

// New common normaliser
double get_trunc_normaliser_cpp(
    const Rcpp::NumericMatrix& p_all_acc,
    RacePdf1Fun pdf1,
    RaceCdf1Fun cdf1,
    LogicalVector isok_trial,
    double LT, double UT,
    int n_lR,
    int n_par,
    double epsilon,
    void* model_specific_context,
    gsl_integration_workspace* w){
  // Calculate sum of log(1-F(LT)) and sum of log(1-F(UT))
  // S(t) = product(1 - F_i(t))
  // P(LT < T < UT) = S(LT) - S(UT)
  // log P = log_diff_exp(log S(LT), log S(UT))
  
  double logS_LT = 0.0;
  double logS_UT = 0.0;
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  std::vector<double> par_row(static_cast<size_t>(n_par));
  
  const bool ut_finite = (UT != R_PosInf);
  for (int k = 0; k < n_lR; ++k) {
    if (!isok_trial[k]) return R_NegInf;
    for (int c = 0; c < n_par; ++c) par_row[static_cast<size_t>(c)] = p_all_acc(k, c);
    
    double cdf_lt = cdf1(LT, par_row.data(), model_specific_context);
    if (!std::isfinite(cdf_lt)) return R_NegInf;
    if (cdf_lt >= 1.0) { logS_LT = R_NegInf; break; }
    if (cdf_lt > 0.0) {
      if (cdf_lt > 1.0 - 1e-15) cdf_lt = 1.0 - 1e-15;
      logS_LT += std::log1p(-cdf_lt);
    }
    
    if (!ut_finite) continue;
    double cdf_ut = cdf1(UT, par_row.data(), model_specific_context);
    if (!std::isfinite(cdf_ut)) return R_NegInf;
    if (cdf_ut >= 1.0) { logS_UT = R_NegInf; }
    else if (cdf_ut > 0.0) {
      if (cdf_ut > 1.0 - 1e-15) cdf_ut = 1.0 - 1e-15;
      logS_UT += std::log1p(-cdf_ut);
    }
  }
  
  if (logS_LT == R_NegInf) return R_NegInf;
  if (!ut_finite) return logS_LT;
  
  double logP = log_diff_exp(logS_LT, logS_UT);
  if (R_FINITE(logP) && logP > log_prob_eps) return logP;
  
  // Fallback: numerical survivor-product failed or underflowed; integrate over winners.
  if (w == nullptr) Rcpp::stop("get_trunc_normaliser_cpp: GSL workspace is null.");
  double log_total = R_NegInf;
  for (int k_win = 1; k_win <= n_lR; ++k_win) {
    double log_k = integrate_for_kth_winner_cpp(k_win,
                                                p_all_acc,
                                                isok_trial,
                                                LT,
                                                UT,
                                                pdf1,
                                                cdf1,
                                                n_lR,
                                                n_par,
                                                epsilon,
                                                model_specific_context,
                                                w);
    log_total = log_sum_exp(log_total, log_k);
  }
  if (R_FINITE(log_total) && log_total > log_prob_eps) return log_total;
  return R_NegInf;
}

// Main C++ function for censored/truncated race likelihood calculation
// This function is now the unified entry point for all race models (LBA, RDM, LNR),
// whether they are standard or explicitly handling censoring/truncation.
// It uses batching for finite RTs and iterative processing for others (censored/NA RTs).
double c_log_likelihood_race_cens_trunc(
    Rcpp::NumericMatrix pars,               // Parameters for one particle, covering all dadm rows for that particle
    Rcpp::DataFrame dadm,                   // Data for unique trial conditions, structured for all accumulators
    RacePdfFun model_dfun,                  // Pointer to the model's PDF adapter function
    RaceCdfFun model_pfun,                  // Pointer to the model's CDF adapter function
    RacePdf1Fun pdf1,                       // Scalar PDF (for GSL)
    RaceCdf1Fun cdf1,                       // Scalar CDF (for GSL)
    const int n_trials,
    LogicalVector winner,
    Rcpp::IntegerVector expand,  // Vector for expanding unique LLs to full trial count
    double min_ll,                          // Minimum log-likelihood value
    const Rcpp::LogicalVector isok,   // Parameter validity for each row in 'pars' matrix
    int n_lR,                              // Number of accumulators in the race (must be > 0 if data exists)
    void* model_context_for_funcs,          // Context for model_dfun/model_pfun (e.g., contains posdrift for LBA)
    bool all_finite_trials             // Data-only hint: all trials finite/in-bounds/known response
) {
  
  // AH  write results to global environment
  // Environment global = Environment::global_env();
  
  // Only allocate a GSL workspace if/when numerical integration is needed.
  GslWorkspacePtr workspace(nullptr, &gsl_integration_workspace_free);
  
  // Fetch censoring and truncation values from dadm columns. These are passed across all rows for ease of access, but should be identical at least at the subject level as they don't correspond to a data entry. Attributes probably a better fit here but clunky.
  Rcpp::NumericVector LT = dadm["LT"];
  Rcpp::NumericVector UT = dadm["UT"];
  Rcpp::NumericVector LC = dadm["LC"];
  Rcpp::NumericVector UC = dadm["UC"];
  double integration_epsilon = 1e-7; // Tolerance for GSL integration
  int n_lR_j = n_lR;
  Rcpp::NumericVector rts_dadm = dadm["rt"];
  Rcpp::IntegerVector R_idxs_dadm = dadm["R"];
  NumericVector lds(n_trials, min_ll); // initialise at min_ll
  // If a RACE column exists, set parameters of accumulators not present on a
  // given trial to NA so the density functions return zero for them. This
  // mirrors logic from the old c_log_likelihood_race implementation.
  const bool has_RACE_col = dadm.containsElementNamed("RACE");
  Rcpp::IntegerVector RACE;
  Rcpp::LogicalVector RACE_mask;
  if (has_RACE_col) {
    bool has_precomputed_race = false;
    if (dadm.hasAttribute("RACE_nacc_by_row") && dadm.hasAttribute("RACE_mask")) {
      RACE = dadm.attr("RACE_nacc_by_row");
      RACE_mask = dadm.attr("RACE_mask");
      has_precomputed_race = (RACE.size() == n_trials && RACE_mask.size() == n_trials);
    }
    if (!has_precomputed_race) {
      Rcpp::IntegerVector lR_dadm = dadm["lR"];
      RACE = Rcpp::IntegerVector(n_trials, n_lR);
      RACE_mask = Rcpp::LogicalVector(n_trials, true); // Mask for all dadm rows
      // factor codes (1-based) for each row
      Rcpp::IntegerVector race_idx = dadm["RACE"];
      // character levels ("2", "3", ...)
      Rcpp::CharacterVector race_levels = race_idx.attr("levels");
      std::vector<int> nacc_by_level;
      nacc_by_level.reserve(static_cast<size_t>(race_levels.size()));
      for (int i = 0; i < race_levels.size(); ++i) {
        nacc_by_level.push_back(std::stoi(Rcpp::as<std::string>(race_levels[i])));
      }
      for (int row = 0; row < pars.nrow(); ++row) {
        // how many accumulators for this trial
        if (race_idx[row] == NA_INTEGER) continue;
        int n_lR_this_trial = nacc_by_level[static_cast<size_t>(race_idx[row] - 1)];
        RACE[row] = n_lR_this_trial;
        // lR_dadm is the (1-based) index of *this* accumulator on the trial
        if (lR_dadm[row] > n_lR_this_trial) {
          // accumulator not present - blank its parameter row
          std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
          RACE_mask[row] = false;
        }
      }
    }
    for (int row = 0; row < pars.nrow(); ++row) {
      if (!RACE_mask[row]) std::fill(pars.row(row).begin(), pars.row(row).end(), NA_REAL);
    }
  }
  if (n_trials == 0) return 0.0; // No data, no likelihood
  
  if (n_lR <= 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: n_lR must be positive and correctly determined before this call.");
  if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_race_cens_trunc: dadm nrows not a multiple of n_lR.");
  
  // Here we check for a pC parameter corresponding to probability of contaminant OMISSION.
  // Index the column (if it exists)
  bool use_pC = false;
  int pc_col = -1;
  Rcpp::List dimnames = pars.attr("dimnames");
  Rcpp::CharacterVector colnames = as<Rcpp::CharacterVector>(dimnames[1]);
  for (int j = 0; j < colnames.size(); ++j) {
    if (as<std::string>(colnames[j]) == "pContaminant") {
      pc_col = j;
      use_pC = true;
      break;
    }
  }
  // Also check that pC is not all zeroes (i.e. a model that has pC optional but is not using it here)
  if (use_pC) {
    bool all_zero = true;
    for (int i = 0; i < pars.nrow(); ++i) {
      if (pars(i, pc_col) != 0.0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero) {
      use_pC = false;
    }
  }
  
  int n_unique_trials = n_trials / n_lR;
  Rcpp::NumericVector ll_unique(n_unique_trials, min_ll);
  Rcpp::NumericVector pC_values(n_unique_trials);
  if (use_pC) {
    for (int j = 0; j < n_unique_trials; ++j) {
      pC_values[j] = pars(j * n_lR, pc_col);
    }
  }
  
  // Parameter matrix and validity vector checks
  if (pars.nrow() != n_trials) {
    Rcpp::Rcout << "pars.nrow(): " << pars.nrow() << ", n_trials: " << n_trials << std::endl;
    Rcpp::stop("c_log_likelihood_race_cens_trunc: pars matrix dimensions do not match total dadm rows.");
  }
  if (isok.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match "
                 "pars matrix rows.");
  }
  if (winner.size() != pars.nrow()) {
    Rcpp::stop("c_log_likelihood_race_cens_trunc: isok size does not match pars matrix rows.");
  }
  const int n_par = pars.ncol();
  std::vector<double> pars_rowmajor_buffer(static_cast<size_t>(n_lR) * n_par);
  std::vector<int> isok_int_buffer(static_cast<size_t>(n_lR), 0);
  std::vector<double> logS_k_buffer(static_cast<size_t>(n_lR), R_NegInf);
  auto fill_trial_buffers = [&](int start_row_idx, int n_lR_j) {
    for (int k = 0; k < n_lR_j; ++k) {
      isok_int_buffer[k] = isok[start_row_idx + k] ? 1 : 0;
      for (int c = 0; c < n_par; ++c) {
        pars_rowmajor_buffer[static_cast<size_t>(k) * n_par + c] = pars(start_row_idx + k, c);
      }
    }
  };
  
  // Fast path hint: computed once per calc_ll call (data-only), to avoid re-scanning per particle.
  const bool use_full_finite_batch = all_finite_trials;
  Rcpp::LogicalVector finite_rt_mask;
  std::vector<int> finite_rt_unique_trial_indices;
  std::vector<int> other_unique_trial_indices;
  bool has_precomputed_finite = false;
  if (!use_full_finite_batch) {
    if (dadm.hasAttribute("finite_rt_mask") &&
        dadm.hasAttribute("finite_rt_unique_trial_indices") &&
        dadm.hasAttribute("other_unique_trial_indices")) {
      finite_rt_mask = dadm.attr("finite_rt_mask");
      Rcpp::IntegerVector finite_attr = dadm.attr("finite_rt_unique_trial_indices");
      Rcpp::IntegerVector other_attr = dadm.attr("other_unique_trial_indices");
      if (finite_rt_mask.size() == n_trials) {
        finite_rt_unique_trial_indices.assign(finite_attr.begin(), finite_attr.end());
        other_unique_trial_indices.assign(other_attr.begin(), other_attr.end());
        has_precomputed_finite = true;
      }
    }
    if (!has_precomputed_finite) {
      finite_rt_mask = Rcpp::LogicalVector(n_trials, false); // Mask for all dadm rows
    }
  }
  double log_Z_this = 0;  // Default inv_Z if no truncation. Should never be used but here as a precaution.
  // Scratch interval bounds for per-trial branches below. Keeping these as
  // mutable locals makes the censoring/truncation logic read like "set bounds,
  // then integrate / take survivor", rather than threading LTj/UTj/LCj/UCj
  // through every helper call.
  double lower_for_trial = 0;
  double upper_for_trial = R_PosInf;
  ContextForRaceModels* ctx = static_cast<ContextForRaceModels*>(model_context_for_funcs);
  ContextForRaceModels dense_ctx = *ctx;
  dense_ctx.min_lik_for_pdf = 0.0;
  void* dense_ctx_ptr = static_cast<void*>(&dense_ctx);
  bool gng=ctx->gng;
  // Batch process finite rt trials
  if (!use_full_finite_batch && !has_precomputed_finite) {
    finite_rt_unique_trial_indices.reserve(n_unique_trials);
    other_unique_trial_indices.reserve(n_unique_trials);
    // Categorize unique trials based on RT properties and parameter validity
    for (int j = 0; j < n_unique_trials; ++j) {
      int start_row_idx = j * n_lR; // Starting row in dadm/pars for this unique trial
      double rt_j = rts_dadm[start_row_idx]; // RT for this unique trial
      int R_j_idx = R_idxs_dadm[start_row_idx]; // Winner index (1-based from R factor)
      int n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
      // Criteria for batch processing: finite, positive RT, and known winner
      if (R_FINITE(rt_j) && rt_j > 0 && R_j_idx != NA_INTEGER) {
        finite_rt_unique_trial_indices.push_back(j); // Store unique trial index
        for(int k_acc = 0; k_acc < n_lR_j; ++k_acc) {
          int dadm_row_idx = start_row_idx + k_acc;
          finite_rt_mask[dadm_row_idx] = true;
        }
      }  else {
        other_unique_trial_indices.push_back(j); // All other cases (NA, Inf, -Inf, outside bounds, or unknown winner with finite RT)
      }
    }
  }
  
  const bool has_finite_batch = use_full_finite_batch || (finite_rt_unique_trial_indices.size() > 0);
  if (has_finite_batch) {
    Rcpp::LogicalVector idx_win;
    Rcpp::LogicalVector idx_loss;
    bool any_win = false;
    bool any_loss = false;
    if (use_full_finite_batch && !has_RACE_col) {
      idx_win = winner;
      any_win = true;
      if (n_lR > 1) {
        idx_loss = !winner;
        any_loss = true;
      }
    } else {
      idx_win = Rcpp::LogicalVector(n_trials, false);
      idx_loss = Rcpp::LogicalVector(n_trials, false);
      for (int i = 0; i < n_trials; ++i) {
        if (has_RACE_col && !RACE_mask[i]) continue;
        if (!use_full_finite_batch && !finite_rt_mask[i]) continue;
        if (winner[i]) { idx_win[i] = true; any_win = true; }
        else { idx_loss[i] = true; any_loss = true; }
      }
    }
    
    if (any_win) {
      Rcpp::NumericVector win_pdf = model_dfun(rts_dadm, pars, idx_win, isok, dense_ctx_ptr);
      int out_i = 0;
      for (int i = 0; i < n_trials; ++i) {
        if (!idx_win[i]) continue;
        const double pdf = win_pdf[out_i++];
        lds[i] = (pdf > 0.0 && std::isfinite(pdf)) ? std::log(pdf) : R_NegInf;
      }
    }
    
  if (n_lR > 1 && any_loss) {
      Rcpp::NumericVector loss_cdf = model_pfun(rts_dadm, pars, idx_loss, isok, dense_ctx_ptr);
      int out_i = 0;
      for (int i = 0; i < n_trials; ++i) {
        if (!idx_loss[i]) continue;
        double cdf = loss_cdf[out_i++];
        cdf = clamp_cdf01_race(cdf);
        lds[i] = safe_log1m_race(cdf);
      }
    }
    
    // Apply truncation correction and calculate log-likelihood for each trial in the batch
    const int n_finite_unique = use_full_finite_batch
    ? n_unique_trials
    : static_cast<int>(finite_rt_unique_trial_indices.size());
    for (int i = 0; i < n_finite_unique; ++i) {
      // When all unique trials are finite, iterate in order 0..n_unique_trials-1.
      // Otherwise, iterate over the precomputed subset of finite-RT unique trials.
      int unique_trial_idx = use_full_finite_batch
        ? i
        : finite_rt_unique_trial_indices[static_cast<size_t>(i)];
      int start_row_idx = unique_trial_idx * n_lR;
      n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
      log_Z_this = 0.0;
      const double LTj = LT[start_row_idx];
      const double UTj = UT[start_row_idx];
      if (LTj != 0.0 || UTj != R_PosInf) { // truncation active
        fill_trial_buffers(start_row_idx, n_lR_j);
        log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                       isok_int_buffer.data(),
                                                       pdf1,
                                                       cdf1,
                                                       LTj,
                                                       UTj,
                                                       n_lR_j,
                                                       n_par,
                                                       integration_epsilon,
                                                       model_context_for_funcs,
                                                       workspace);
      }
      double current_trial_ll_sum = 0.0;
      for(int k=0; k < n_lR_j; ++k) {
        int dadm_row_idx = start_row_idx + k; // Correct index into lds
        current_trial_ll_sum += lds[dadm_row_idx];
      }
      
      if ((LT[start_row_idx] != 0 || UT[start_row_idx] != R_PosInf) && (NumericVector::is_na(log_Z_this) || !R_FINITE(log_Z_this))) {
        // If truncation active but inv_Z is bad, probability is effectively zero. Could mean the probability of observing an untruncated RT is zero, which would be bad.
        ll_unique[unique_trial_idx] = min_ll;
        continue;
      }
      else if (LT[start_row_idx] != 0 || UT[start_row_idx] != R_PosInf) { // Truncation active and inv_Z is good
        current_trial_ll_sum -=  log_Z_this;
      }
      
      ll_unique[unique_trial_idx] = current_trial_ll_sum;
      ll_unique[unique_trial_idx] = std::max(min_ll, ll_unique[unique_trial_idx]);
    }
    
    //        global["llC"]  = ll_unique;
    
  }
  // --- Process other trials (Infinite RTs, NA RTs, or finite RTs outside truncation) ---
  // These trials require individual processing, often involving numerical integration for censored intervals.
  // TODO update Case 3 to handle go/no-go withheld response
  double current_ll_val;
  const double log_prob_eps = std::log(std::numeric_limits<double>::epsilon());
  for (size_t i = 0; i < other_unique_trial_indices.size(); ++i) {
    int unique_trial_idx = other_unique_trial_indices[i];
    int start_row_idx = unique_trial_idx*n_lR;
    log_Z_this = 0.0;
    n_lR_j = has_RACE_col ? RACE[start_row_idx] : n_lR;
    const double rt_j = rts_dadm[start_row_idx];
    const int R_j_idx = R_idxs_dadm[start_row_idx];
    const double LTj = LT[start_row_idx];
    const double UTj = UT[start_row_idx];
    const double LCj = LC[start_row_idx];
    const double UCj = UC[start_row_idx];
    const bool has_trunc = (LTj != 0.0 || UTj != R_PosInf);
    
    const bool needs_model = (rt_j == R_NegInf) || (rt_j == R_PosInf) || Rcpp::NumericVector::is_na(rt_j) ||
      (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER);
    if (!needs_model) {
      ll_unique[unique_trial_idx] = min_ll;
      continue;
    }
    
    fill_trial_buffers(start_row_idx, n_lR_j);
    
    auto integrate_interval = [&](int k_winner_1based, double low, double upp) -> double {
      gsl_integration_workspace* w = ensure_gsl_workspace(workspace);
      return integrate_for_kth_winner_rowmajor_cpp(k_winner_1based,
                                                   pars_rowmajor_buffer.data(),
                                                   isok_int_buffer.data(),
                                                   low,
                                                   upp,
                                                   pdf1,
                                                   cdf1,
                                                   n_lR_j,
                                                   n_par,
                                                   integration_epsilon,
                                                   model_context_for_funcs,
                                                   w);
    };
    
    auto log_survivor = [&](double t) -> double {
      return log_survivor_rowmajor(t,
                                   pars_rowmajor_buffer.data(),
                                   isok_int_buffer.data(),
                                   n_lR_j,
                                   n_par,
                                   cdf1,
                                   model_context_for_funcs);
    };
    
    current_ll_val = R_NegInf;
    if (rt_j == R_NegInf) {
      lower_for_trial = LTj;
      upper_for_trial = LCj;
      if (R_j_idx == NA_INTEGER) {
        const double logP = log_diff_exp(log_survivor(lower_for_trial), log_survivor(upper_for_trial));
        if (R_FINITE(logP) && logP > log_prob_eps) {
          current_ll_val = logP;
        } else {
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
          }
        }
      } else {
        current_ll_val = integrate_interval(R_j_idx, lower_for_trial, upper_for_trial);
      }
    } else if (rt_j == R_PosInf) {
      if (gng) {
        const double D = UCj;
        lower_for_trial = LTj;
        upper_for_trial = D;
        int k_nogo = -1; // 1-based
        int n_true = 0;
        for (int k = 0; k < n_lR_j; ++k) {
          if (winner[start_row_idx + k]) { n_true++; k_nogo = k + 1; }
        }
        if (n_true != 1) Rcpp::stop("No winner identified in go/no-go withheld response");
        const double logA = integrate_interval(k_nogo, lower_for_trial, upper_for_trial);
        const double logB = log_survivor(D);
        current_ll_val = log_sum_exp(logA, logB);
      } else {
        lower_for_trial = UCj;
        upper_for_trial = UTj;
        if (n_lR_j == 1) {
          const double cdf = cdf1(lower_for_trial, pars_rowmajor_buffer.data(), model_context_for_funcs);
          current_ll_val = safe_log1m_race(clamp_cdf01_race(cdf));
        } else if (R_j_idx == NA_INTEGER) {
          const double logP = log_diff_exp(log_survivor(lower_for_trial), log_survivor(upper_for_trial));
          if (R_FINITE(logP) && logP > log_prob_eps) {
            current_ll_val = logP;
          } else {
            for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
              current_ll_val = log_sum_exp(current_ll_val, integrate_interval(k_win, lower_for_trial, upper_for_trial));
            }
          }
        } else {
          current_ll_val = integrate_interval(R_j_idx, lower_for_trial, upper_for_trial);
        }
      }
    } else if (Rcpp::NumericVector::is_na(rt_j)) {
      const double lower1 = LTj;
      const double upper1 = LCj;
      const double lower2 = UCj;
      const double upper2 = UTj;
      if (R_j_idx != NA_INTEGER) {
        current_ll_val = log_sum_exp(integrate_interval(R_j_idx, lower1, upper1),
                                     integrate_interval(R_j_idx, lower2, upper2));
      } else {
        const double logP1 = log_diff_exp(log_survivor(lower1), log_survivor(upper1));
        const double logP2 = log_diff_exp(log_survivor(lower2), log_survivor(upper2));
        const double logPsum = log_sum_exp(logP1, logP2);
        if (R_FINITE(logPsum) && logPsum > log_prob_eps) {
          current_ll_val = logPsum;
        } else {
          for (int k_win = 1; k_win <= n_lR_j; ++k_win) {
            const double ll_L_k = integrate_interval(k_win, lower1, upper1);
            const double ll_U_k = integrate_interval(k_win, lower2, upper2);
            current_ll_val = log_sum_exp(current_ll_val, log_sum_exp(ll_L_k, ll_U_k));
          }
        }
      }
    } else if (R_FINITE(rt_j) && rt_j > 0.0 && R_j_idx == NA_INTEGER) {
      current_ll_val = log_min_density_rowmajor(rt_j,
                                                pars_rowmajor_buffer.data(),
                                                isok_int_buffer.data(),
                                                n_lR_j,
                                                n_par,
                                                pdf1,
                                                cdf1,
                                                model_context_for_funcs,
                                                logS_k_buffer);
    }
    
    if (current_ll_val > min_ll && has_trunc) {
      log_Z_this = get_trunc_normaliser_rowmajor_cpp(pars_rowmajor_buffer.data(),
                                                     isok_int_buffer.data(),
                                                     pdf1,
                                                     cdf1,
                                                     LTj,
                                                     UTj,
                                                     n_lR_j,
                                                     n_par,
                                                     integration_epsilon,
                                                     model_context_for_funcs,
                                                     workspace);
      if (!R_FINITE(log_Z_this)) current_ll_val = min_ll;
      else current_ll_val -= log_Z_this;
    }
    
    ll_unique[unique_trial_idx] = std::max(min_ll, current_ll_val);
  }
  
  //    global["llC1"]  = ll_unique;
  
  // --- Summation of log-likelihoods for all unique trials ---
  double total_ll = 0;
  if (expand.length() > 0) { // If an expansion vector is provided (e.g. from non-compressed dadm)
    if (use_pC) { // multiply all likelihoods by the appropriate pC adjustment if it's present
      for (int j = 0; j < n_unique_trials; ++j) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC); // log(1-pC)
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        
        // AH changed to match contamination signaled by +Inf
        if (rt_j == R_PosInf) {
          double term1 = std::log(pC);
          double term2 = log1m_pC + ll_unique[j];
          ll_unique[j] = log_sum_exp(term1, term2); // pC * pIO * pCens
        } else {
          ll_unique[j] = log1m_pC + ll_unique[j];    // pRT * 1-pC
        }
        
        
      }
    }
    for (int i = 0; i < expand.length(); ++i) {
      total_ll += ll_unique[expand[i] - 1];
    }
  } else { // Default: sum ll_unique directly (each unique trial counted once)
    for (int j = 0; j < n_unique_trials; ++j) {
      if (use_pC) {
        double pC = pC_values[j];
        double log1m_pC = std::log1p(-pC);
        int start_row_idx = j * n_lR;
        double rt_j = rts_dadm[start_row_idx];
        
        if (R_FINITE(rt_j)) {
          ll_unique[j] = log1m_pC + ll_unique[j];
        } else {
          double term1 = std::log(pC);
          double term2 = log1m_pC + ll_unique[j];
          ll_unique[j] = log_sum_exp(term1, term2);
        }
      }
      total_ll += ll_unique[j];
    }
  }
  
  //    global["llC1C"]  = ll_unique;
  
  return total_ll;
}

namespace {

struct LogicalRulesYesCdfParams {
  const double* par_target;      // length n_par
  const double* par_nontarget;   // length n_par
  ContextForRaceModels* ctx;     // provides pdf1/cdf1 + any model options (e.g., posdrift)
};

struct LogicalRulesIntegrationStats {
  int64_t qng_ok = 0;
  int64_t qng_fail = 0;
  int64_t qag_ok = 0;
  int64_t qag_fail = 0;
  int64_t qng_neval = 0;
};

inline double clamp01(double x) {
  if (!R_FINITE(x)) return NA_REAL;
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x;
}

double logicalrules_yes_cdf_integrand(double u, void* vp) {
  auto* p = static_cast<LogicalRulesYesCdfParams*>(vp);
  if (!p || !p->ctx || !p->ctx->pdf1 || !p->ctx->cdf1) return 0.0;
  if (!R_FINITE(u) || u <= 0.0) return 0.0;

  const double f = p->ctx->pdf1(u, p->par_target, p->ctx);
  if (!R_FINITE(f) || f <= 0.0) return 0.0;

  const double Fn = clamp01(p->ctx->cdf1(u, p->par_nontarget, p->ctx));
  if (!R_FINITE(Fn)) return 0.0;
  const double sn = 1.0 - Fn;
  return f * ((sn >= 0.0) ? sn : 0.0);
}

double integrate_logicalrules_yes_cdf(double t,
                                     const double* par_target,
                                     const double* par_nontarget,
                                     ContextForRaceModels* ctx,
                                     gsl_integration_workspace* w,
                                     double epsrel,
                                     LogicalRulesIntegrationStats* stats,
                                     bool try_qng) {
  if (!R_FINITE(t) || t <= 0.0) return 0.0;
  if (!ctx || !ctx->pdf1 || !ctx->cdf1) return NA_REAL;
  if (!w) return NA_REAL;

  // IMPORTANT: GSL's default error handler aborts the process on non-success
  // statuses (including qng "failed to reach tolerance"). We want to detect
  // failure and fall back to an adaptive integrator instead.
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  struct HandlerRestore {
    gsl_error_handler_t* h;
    ~HandlerRestore() { gsl_set_error_handler(h); }
  } restore{old_handler};

  LogicalRulesYesCdfParams p{par_target, par_nontarget, ctx};
  gsl_function F;
  F.function = &logicalrules_yes_cdf_integrand;
  F.params = &p;

  // Fast path: non-adaptive Gauss-Kronrod (usually much faster than qag).
  // Falls back to adaptive integration if qng fails.
  if (try_qng) {
    double result = NA_REAL;
    double abserr = NA_REAL;
    size_t neval = 0;
    const int status = gsl_integration_qng(&F,
                                          0.0,
                                          t,
                                          0.0,    // epsabs
                                          epsrel, // epsrel
                                          &result,
                                          &abserr,
                                          &neval);
    if (stats) stats->qng_neval += static_cast<int64_t>(neval);
    if (status == GSL_SUCCESS && R_FINITE(result)) {
      if (stats) stats->qng_ok += 1;
      if (result < 0.0) result = 0.0;
      return result;
    }
    if (stats) stats->qng_fail += 1;

    if (ctx && ctx->log_out) {
      static int printed = 0;
      int max_print = 5;
      if (const char* v = std::getenv("EMC2_LOGICALRULES_LOG_INTEGRATION_FAILS")) {
        char* endp = nullptr;
        const long vv = std::strtol(v, &endp, 10);
        if (endp != v && vv >= 0) max_print = static_cast<int>(vv);
      }
      if (printed < max_print) {
        ++printed;
        Rprintf("LogicalRules qng fail (status=%d: %s): t=%.6g epsrel=%.2g abserr=%.2g neval=%llu\n",
                status, gsl_strerror(status), t, epsrel, abserr,
                static_cast<unsigned long long>(neval));
      }
    }
  }

  double result = NA_REAL;
  double abserr = NA_REAL;
  const int status = gsl_integration_qag(&F,
                                        0.0,
                                        t,
                                        0.0,    // epsabs
                                        epsrel, // epsrel
                                        w->limit,
                                        GSL_INTEG_GAUSS21, // cheaper than 61-point; adaptive handles hard cases
                                        w,
                                        &result,
                                        &abserr);
  if (status == GSL_SUCCESS && R_FINITE(result)) {
    if (stats) stats->qag_ok += 1;
  } else {
    if (stats) stats->qag_fail += 1;
  }
  if (status != GSL_SUCCESS || !R_FINITE(result)) return NA_REAL;
  if (result < 0.0) result = 0.0;
  return result;
}

} // namespace

// This is the Bushmakin et al (2017) model
double c_log_likelihood_logicalrules(
  Rcpp::NumericMatrix pars,
  Rcpp::DataFrame dadm,
  RacePdfFun model_dfun,
  RaceCdfFun model_pfun,
  const int n_trials,
  const Rcpp::IntegerVector expand,
  double min_ll,
  const Rcpp::LogicalVector ok_params,
	int n_acc,
  void* model_specific_context) {

  if (n_trials == 0) return 0.0;
  if (n_acc <= 0) Rcpp::stop("c_log_likelihood_redundant_target_race: n_acc must be positive");
  if (n_trials % n_acc != 0) Rcpp::stop("c_log_likelihood_redundant_target_race: dadm rows must be multiple of n_acc");
  if (pars.nrow() != n_trials) Rcpp::stop("c_log_likelihood_redundant_target_race: pars nrow must match n_trials");
  if (ok_params.size() != n_trials) Rcpp::stop("c_log_likelihood_redundant_target_race: ok_params length must match n_trials");
  if (!dadm.containsElementNamed("rt") ||
      !dadm.containsElementNamed("lR") ||
      !dadm.containsElementNamed("R") ||
      !dadm.containsElementNamed("LogicalRule")) {
    Rcpp::stop("c_log_likelihood_redundant_target_race: dadm must contain columns rt, lR, R, LogicalRule");
  }

  Rcpp::NumericVector rts = dadm["rt"];
  SEXP role_sexp = dadm["lR"];
  SEXP resp_sexp = dadm["R"];
  SEXP rule_sexp = dadm["LogicalRule"];
  if (rts.size() != n_trials ||
      Rf_length(role_sexp) != n_trials ||
      Rf_length(resp_sexp) != n_trials ||
      Rf_length(rule_sexp) != n_trials) {
    Rcpp::stop("c_log_likelihood_redundant_target_race: dadm column lengths must match n_trials");
  }

  // Avoid repeated factor -> character coercion (very costly inside PMwG likelihood loops).
  const bool role_is_factor = Rf_inherits(role_sexp, "factor");
  const bool resp_is_factor = Rf_inherits(resp_sexp, "factor");
  const bool rule_is_factor = Rf_inherits(rule_sexp, "factor");

  Rcpp::IntegerVector role_code;
  Rcpp::IntegerVector resp_code;
  Rcpp::IntegerVector rule_code;
  Rcpp::CharacterVector role_chr;
  Rcpp::CharacterVector resp_chr;
  Rcpp::CharacterVector rule_chr;

  auto level_code = [](const Rcpp::CharacterVector& levels, const std::string& target) -> int {
    for (int i = 0; i < levels.size(); ++i) {
      if (Rcpp::as<std::string>(levels[i]) == target) return i + 1; // factor codes are 1-based
    }
    return -1;
  };

  int codeA = -1, codeB = -1, codenA = -1, codenB = -1;
  int codeYes = -1, codeNo = -1;
  int codeAND = -1, codeOR = -1;

  if (role_is_factor) {
    role_code = Rcpp::IntegerVector(role_sexp);
    const Rcpp::CharacterVector lev = role_code.attr("levels");
    codeA = level_code(lev, "A");
    codeB = level_code(lev, "B");
    codenA = level_code(lev, "n_A");
    codenB = level_code(lev, "n_B");
  } else {
    role_chr = Rcpp::CharacterVector(role_sexp);
  }

  if (resp_is_factor) {
    resp_code = Rcpp::IntegerVector(resp_sexp);
    const Rcpp::CharacterVector lev = resp_code.attr("levels");
    codeYes = level_code(lev, "yes");
    codeNo = level_code(lev, "no");
  } else {
    resp_chr = Rcpp::CharacterVector(resp_sexp);
  }

  if (rule_is_factor) {
    rule_code = Rcpp::IntegerVector(rule_sexp);
    const Rcpp::CharacterVector lev = rule_code.attr("levels");
    codeAND = level_code(lev, "AND");
    codeOR = level_code(lev, "OR");
  } else {
    rule_chr = Rcpp::CharacterVector(rule_sexp);
  }

  // No concept of "winner" here: request PDF/CDF for all accumulators, respecting ok_params bounds.
  Rcpp::LogicalVector winner_all(n_trials, true);
  Rcpp::NumericVector f_all = model_dfun(rts, pars, winner_all, ok_params, model_specific_context);
  Rcpp::NumericVector F_all = model_pfun(rts, pars, winner_all, ok_params, model_specific_context);
  const double* f_ptr = f_all.begin();
  const double* F_ptr = F_all.begin();

  const int n_unique_trials = n_trials / n_acc;
  Rcpp::NumericVector ll_unique(n_unique_trials, min_ll);
 	
	// Here we check for a rule-following parameter (p) corresponding to probability of processing only a single channel with probability q.
	// Index the column (if it exists)
	bool use_rulebreak = false;
  int rulebreak_col = -1;
	int q_col = -1;
  Rcpp::CharacterVector colnames;
  if (pars.hasAttribute("dimnames")) {
    Rcpp::List dimnames = pars.attr("dimnames");
    if (dimnames.size() >= 2 && dimnames[1] != R_NilValue) {
      colnames = as<Rcpp::CharacterVector>(dimnames[1]);
    }
  }
  if (colnames.size() > 0) {
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "p") {
        rulebreak_col = j;
        break;
      }
    }
  }
  if (rulebreak_col >= 0) {
    // If p is present but effectively hardcoded to 1.0 for all rows, disable rulebreaking.
    bool all_one = true;
    for (int i = 0; i < pars.nrow(); ++i) {
      const double p_val = pars(i, rulebreak_col);
      if (R_FINITE(p_val) && p_val != 1.0) {
        all_one = false;
        break;
      }
    }
    use_rulebreak = !all_one;
  }

	// If we're rulebreaking, find the index of column q, then also check that the rulebreak probability isn't hardcoded to zero
  if (use_rulebreak) {
    for (int j = 0; j < colnames.size(); ++j) {
      if (as<std::string>(colnames[j]) == "q") {
        q_col = j;
        break;
      }
    }
    // Require both p and q columns for rule-breaking mixture.
    if (q_col < 0) use_rulebreak = false;
  }

	
	// set up parameters
	double p_j = 0.0;
  constexpr double kMinSurv = 1e-12;
  auto clamp01_finite = [](double x) -> double {
    if (!R_FINITE(x)) return NA_REAL;
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    return x;
  };
  auto safe_surv_from_cdf = [&](double cdf) -> double {
    const double F = clamp01_finite(cdf);
    if (!R_FINITE(F)) return NA_REAL;
    const double s = std::fma(-F, 1.0, 1.0); // 1 - F
    return (s >= kMinSurv) ? s : kMinSurv;
  };
  auto safe_surv_from_prod_cdf = [&](double cdf1, double cdf2) -> double {
    const double F1 = clamp01_finite(cdf1);
    const double F2 = clamp01_finite(cdf2);
    if (!R_FINITE(F1) || !R_FINITE(F2)) return NA_REAL;
    const double s = std::fma(-F1, F2, 1.0); // 1 - F1*F2 (single-rounding)
    return (s >= kMinSurv) ? s : kMinSurv;
  };
	
  // Workspace for numerical integration (local-race CDF terms).
  GslWorkspacePtr workspace(nullptr, &gsl_integration_workspace_free);
  auto ensure_workspace = [&]() {
    if (!workspace) workspace.reset(gsl_integration_workspace_alloc(1000));
  };

  auto* ctx = static_cast<ContextForRaceModels*>(model_specific_context);
  const double epsrel = 1e-7;
  const int n_par = pars.ncol();
  LogicalRulesIntegrationStats stats;
  bool try_qng = true;

  // Reuse buffers to avoid per-trial allocations.
  std::vector<double> parA(static_cast<size_t>(n_par));
  std::vector<double> parB(static_cast<size_t>(n_par));
  std::vector<double> parnA(static_cast<size_t>(n_par));
  std::vector<double> parnB(static_cast<size_t>(n_par));

  for(int j=0; j<n_unique_trials; ++j){
        int start = j*n_acc;
		    p_j = 0.0;
        double fA = NA_REAL, fB = NA_REAL, fnA = NA_REAL, fnB = NA_REAL;
        double FA = NA_REAL, FB = NA_REAL, FnA = NA_REAL, FnB = NA_REAL;
        int idxA = -1, idxB = -1, idxnA = -1, idxnB = -1;
        for(int k=0;k<n_acc;++k){
            int idx = start+k;
            if (role_is_factor) {
              const int rc = role_code[idx];
              if (rc == codeA) { fA = f_ptr[idx]; FA = F_ptr[idx]; idxA = idx; }
              else if (rc == codeB) { fB = f_ptr[idx]; FB = F_ptr[idx]; idxB = idx; }
              else if (rc == codenA) { fnA = f_ptr[idx]; FnA = F_ptr[idx]; idxnA = idx; }
              else if (rc == codenB) { fnB = f_ptr[idx]; FnB = F_ptr[idx]; idxnB = idx; }
            } else {
              const Rcpp::String r = role_chr[idx];
              if (r == "A") { fA = f_ptr[idx]; FA = F_ptr[idx]; idxA = idx; }
              else if (r == "B") { fB = f_ptr[idx]; FB = F_ptr[idx]; idxB = idx; }
              else if (r == "n_A") { fnA = f_ptr[idx]; FnA = F_ptr[idx]; idxnA = idx; }
              else if (r == "n_B") { fnB = f_ptr[idx]; FnB = F_ptr[idx]; idxnB = idx; }
            }
        }

        // Guard: if any required pieces are missing/non-finite, return min_ll for this unique trial.
        if (!R_FINITE(fA) || !R_FINITE(fB) || !R_FINITE(fnA) || !R_FINITE(fnB) ||
            !R_FINITE(FA) || !R_FINITE(FB) || !R_FINITE(FnA) || !R_FINITE(FnB)) {
          ll_unique[j] = min_ll;
          continue;
        }
        if (idxA < 0 || idxB < 0 || idxnA < 0 || idxnB < 0) {
          ll_unique[j] = min_ll;
          continue;
        }
        if (fA < 0.0) fA = 0.0;
        if (fB < 0.0) fB = 0.0;
        if (fnA < 0.0) fnA = 0.0;
        if (fnB < 0.0) fnB = 0.0;
        FA = clamp_cdf01(FA);
        FB = clamp_cdf01(FB);
        FnA = clamp_cdf01(FnA);
        FnB = clamp_cdf01(FnB);
        if (!R_FINITE(FA) || !R_FINITE(FB) || !R_FINITE(FnA) || !R_FINITE(FnB)) {
          ll_unique[j] = min_ll;
          continue;
        }

        const double one_m_FB = safe_surv_from_cdf(FB);
        const double one_m_FA = safe_surv_from_cdf(FA);
        const double one_m_FnB = safe_surv_from_cdf(FnB);
        const double one_m_FnA = safe_surv_from_cdf(FnA);
        if (!R_FINITE(one_m_FB) || !R_FINITE(one_m_FA) || !R_FINITE(one_m_FnB) ||
            !R_FINITE(one_m_FnA)) {
          ll_unique[j] = min_ll;
          continue;
        }

        const bool rule_is_or = rule_is_factor ? (rule_code[start] == codeOR) : (Rcpp::String(rule_chr[start]) == "OR");
        const bool rule_is_and = rule_is_factor ? (rule_code[start] == codeAND) : (Rcpp::String(rule_chr[start]) == "AND");
        if (!(rule_is_or || rule_is_and)) {
          ll_unique[j] = min_ll;
          continue;
        }
        const bool resp_is_yes = resp_is_factor ? (resp_code[start] == codeYes) : (Rcpp::String(resp_chr[start]) == "yes");
        const bool resp_is_no = resp_is_factor ? (resp_code[start] == codeNo) : (Rcpp::String(resp_chr[start]) == "no");
        if (!(resp_is_yes || resp_is_no)) {
          ll_unique[j] = min_ll;
          continue;
        }

        // Parameters for local (two-accumulator) races, copied into contiguous row-major buffers.
        for (int c = 0; c < n_par; ++c) {
          parA[static_cast<size_t>(c)] = pars(idxA, c);
          parB[static_cast<size_t>(c)] = pars(idxB, c);
          parnA[static_cast<size_t>(c)] = pars(idxnA, c);
          parnB[static_cast<size_t>(c)] = pars(idxnB, c);
        }

        // Subdistribution CDFs for "target wins local race by time t" in each channel.
        // G_yes(t) = ∫_0^t f_target(u) * (1 - F_nontarget(u)) du
        ensure_workspace();
        if (!workspace || !ctx || !ctx->pdf1 || !ctx->cdf1) {
          ll_unique[j] = min_ll;
          continue;
        }
        const double t = rts[start];

        // Heuristic: if qng is repeatedly failing, stop trying it to avoid the "double cost"
        // (failed qng + then qag) on every subsequent integral within this likelihood eval.
        if (stats.qng_fail >= 50 && stats.qng_ok <= 1) {
          try_qng = false;
        }

        double GA_yes = integrate_logicalrules_yes_cdf(t, parA.data(), parnA.data(), ctx, workspace.get(), epsrel, &stats, try_qng);
        double GB_yes = integrate_logicalrules_yes_cdf(t, parB.data(), parnB.data(), ctx, workspace.get(), epsrel, &stats, try_qng);
        if (!R_FINITE(GA_yes) || !R_FINITE(GB_yes)) {
          ll_unique[j] = min_ll;
          continue;
        }

        // Channel decision-time CDF: min(target, nontarget) <= t
        const double GA_dec = clamp01_finite(1.0 - (one_m_FA * one_m_FnA));
        const double GB_dec = clamp01_finite(1.0 - (one_m_FB * one_m_FnB));
        if (!R_FINITE(GA_dec) || !R_FINITE(GB_dec)) {
          ll_unique[j] = min_ll;
          continue;
        }
        // Numerical guard: enforce 0 <= G_yes <= G_dec, then derive G_no = G_dec - G_yes.
        if (GA_yes < 0.0) GA_yes = 0.0;
        if (GB_yes < 0.0) GB_yes = 0.0;
        if (GA_yes > GA_dec) GA_yes = GA_dec;
        if (GB_yes > GB_dec) GB_yes = GB_dec;
        const double GA_no = GA_dec - GA_yes;
        const double GB_no = GB_dec - GB_yes;

        // Local winner densities at time t.
        const double gA_yes = fA * one_m_FnA;
        const double gB_yes = fB * one_m_FnB;
        const double gA_no = fnA * one_m_FA;
        const double gB_no = fnB * one_m_FB;

        // Rule-breaking mixture (single-channel processing) components, if enabled.
        bool apply_mix = false;
        double p_follow = NA_REAL;
        double q_chooseB = NA_REAL;
        if (use_rulebreak) {
          p_follow = clamp01_finite(pars(start, rulebreak_col));
          q_chooseB = clamp01_finite(pars(start, q_col));
          apply_mix = R_FINITE(p_follow) && R_FINITE(q_chooseB) && (p_follow < 1.0);
        }
        const double p_rulebreak_yes = (R_FINITE(q_chooseB))
          ? (q_chooseB * gB_yes + (1.0 - q_chooseB) * gA_yes)
          : NA_REAL;
        const double p_rulebreak_no = (R_FINITE(q_chooseB))
          ? (q_chooseB * gB_no + (1.0 - q_chooseB) * gA_no)
          : NA_REAL;

        if (rule_is_or) {
          if (resp_is_yes) {
            p_j = gA_yes * safe_surv_from_cdf(GB_yes) + gB_yes * safe_surv_from_cdf(GA_yes);
            if (apply_mix && R_FINITE(p_rulebreak_yes)) {
              p_j = p_follow * p_j + (1.0 - p_follow) * p_rulebreak_yes;
            }
          } else if (resp_is_no) {
            p_j = gA_no * GB_no + gB_no * GA_no;
            if (apply_mix && R_FINITE(p_rulebreak_no)) {
              p_j = p_follow * p_j + (1.0 - p_follow) * p_rulebreak_no;
            }
          }
        } else if (rule_is_and) {
          if (resp_is_no) {
            p_j = gA_no * safe_surv_from_cdf(GB_no) + gB_no * safe_surv_from_cdf(GA_no);
            if (apply_mix && R_FINITE(p_rulebreak_no)) {
              p_j = p_follow * p_j + (1.0 - p_follow) * p_rulebreak_no;
            }
          } else if (resp_is_yes) {
            p_j = gA_yes * GB_yes + gB_yes * GA_yes;
            if (apply_mix && R_FINITE(p_rulebreak_yes)) {
              p_j = p_follow * p_j + (1.0 - p_follow) * p_rulebreak_yes;
            }
          }
        }

        if (p_j <= 0.0 || !R_FINITE(p_j)) {
          ll_unique[j] = min_ll;
        } else {
          ll_unique[j] = std::log(p_j);
        }
    }

  double sum_ll = 0.0;
  if (expand.length() > 0) {
    Rcpp::NumericVector ll_exp = c_expand(ll_unique, expand);
    for (int i = 0; i < ll_exp.size(); ++i) {
      double val = ll_exp[i];
      if (!R_FINITE(val) || Rcpp::NumericVector::is_na(val) || val < min_ll) val = min_ll;
      sum_ll += val;
    }
  } else {
    for (int j = 0; j < ll_unique.size(); ++j) {
      double val = ll_unique[j];
      if (!R_FINITE(val) || Rcpp::NumericVector::is_na(val) || val < min_ll) val = min_ll;
      sum_ll += val;
    }
  }

  if (ctx && ctx->log_out) {
    Rprintf("LogicalRules integration summary: qng ok=%lld fail=%lld; qag ok=%lld fail=%lld; qng neval=%lld\n",
            static_cast<long long>(stats.qng_ok),
            static_cast<long long>(stats.qng_fail),
            static_cast<long long>(stats.qag_ok),
            static_cast<long long>(stats.qag_fail),
            static_cast<long long>(stats.qng_neval));
  }
  return sum_ll;
}

namespace {

enum class SsGoDist { EXG, RDM };
  
  inline int find_col_idx(const Rcpp::CharacterVector& cols, const std::string& name) {
    for (int i = 0; i < cols.size(); ++i) {
      if (Rcpp::as<std::string>(cols[i]) == name) return i;
    }
    return -1;
  }
  
  inline double safe_log1m(double p) {
    if (!std::isfinite(p)) return R_NegInf;
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return R_NegInf;
    if (p > 1.0 - 1e-15) p = 1.0 - 1e-15;
    return std::log1p(-p);
  }
  
  struct gsl_ss_stopwin_params {
    const double* go_pars; // row-major: n_go * n_go_par
    int n_go;
    int n_go_par;
    RaceCdf1Fun go_cdf1;
    void* go_ctx;
    double SSD;
    double muS;
    double sigmaS;
    double tauS;
  };
  
  double gsl_f_ss_stopwin(double s, void* p) {
    auto* P = static_cast<gsl_ss_stopwin_params*>(p);
    const double f_stop = ss_exgauss_pdf(s, P->muS, P->sigmaS, P->tauS);
    if (!(f_stop > 0.0) || !std::isfinite(f_stop)) return 0.0;
    
    const double t_abs = s + P->SSD;
    double logS = 0.0;
    for (int k = 0; k < P->n_go; ++k) {
      const double* par_k = P->go_pars + static_cast<size_t>(k) * P->n_go_par;
      double cdf = P->go_cdf1(t_abs, par_k, P->go_ctx);
      cdf = clamp_cdf01(cdf);
      const double ll = safe_log1m(cdf);
      if (!std::isfinite(ll)) return 0.0;
      logS += ll;
    }
    const double surv_prod = std::exp(logS);
    if (!(surv_prod > 0.0) || !std::isfinite(surv_prod)) return 0.0;
    return f_stop * surv_prod;
  }
  
  struct gsl_ss_stoptrunc_params {
    const double* go_pars; // row-major: n_go * n_go_par
    int n_go;
    int n_go_par;
    RacePdf1Fun go_pdf1;
    RaceCdf1Fun go_cdf1;
    void* go_ctx;
    double tf;
    double SSD;
    double muS;
    double sigmaS;
    double tauS;
  };
  
  double gsl_f_ss_stoptrunc(double t, void* p) {
    auto* P = static_cast<gsl_ss_stoptrunc_params*>(p);
    if (!(t > 0.0) || !std::isfinite(t)) return 0.0;
    
    double S_all = 1.0;
    double sum_pdf_over_S = 0.0;
    for (int k = 0; k < P->n_go; ++k) {
      const double* par_k = P->go_pars + static_cast<size_t>(k) * P->n_go_par;
      double cdf = P->go_cdf1(t, par_k, P->go_ctx);
      cdf = clamp_cdf01(cdf);
      if (cdf >= 1.0) return 0.0;
      const double S_k = 1.0 - cdf;
      if (!(S_k > 0.0) || !std::isfinite(S_k)) return 0.0;
      S_all *= S_k;
      if (!(S_all > 0.0) || !std::isfinite(S_all)) return 0.0;
      
      double pdf = P->go_pdf1(t, par_k, P->go_ctx);
      if (!(pdf > 0.0) || !std::isfinite(pdf)) continue;
      sum_pdf_over_S += pdf / S_k;
    }
    const double f_min_go = S_all * sum_pdf_over_S;
    if (!(f_min_go > 0.0) || !std::isfinite(f_min_go)) return 0.0;
    
    double S_stop = ss_stop_surv_abs(t, P->SSD, P->muS, P->sigmaS, P->tauS);
    if (!std::isfinite(S_stop)) S_stop = 0.0;
    S_stop = clamp_prob01(S_stop);
    
    const double mix = P->tf + (1.0 - P->tf) * S_stop;
    if (!(mix > 0.0) || !std::isfinite(mix)) return 0.0;
    return f_min_go * mix;
  }
  
  double integrate_ss_stopwin_prob(double upper_duration,
                                   const double* go_pars,
                                   int n_go,
                                   int n_go_par,
                                   RaceCdf1Fun go_cdf1,
                                   void* go_ctx,
                                   double SSD,
                                   double muS,
                                   double sigmaS,
                                   double tauS,
                                   double epsilon,
                                   gsl_integration_workspace* w) {
    if (w == nullptr) Rcpp::stop("integrate_ss_stopwin_prob: GSL workspace is null.");
    if (n_go <= 0) return 0.0;
    if (!R_finite(upper_duration)) upper_duration = R_PosInf;
    if (upper_duration <= 0.0) return 0.0;
    
    gsl_function F;
    gsl_ss_stopwin_params P;
    P.go_pars = go_pars;
    P.n_go = n_go;
    P.n_go_par = n_go_par;
    P.go_cdf1 = go_cdf1;
    P.go_ctx = go_ctx;
    P.SSD = SSD;
    P.muS = muS;
    P.sigmaS = sigmaS;
    P.tauS = tauS;
    F.function = &gsl_f_ss_stopwin;
    F.params = &P;
    
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status;
    double result = 0.0;
    double error = 0.0;
    
    if (upper_duration == R_PosInf) {
      status = gsl_integration_qagiu(&F, 0.0, 0, epsilon, 1000, w, &result, &error);
    } else {
      status = gsl_integration_qags(&F, 0.0, upper_duration, 0, epsilon, 1000, w, &result, &error);
    }
    
    gsl_set_error_handler(old_handler);
    if (status != GSL_SUCCESS) return 0.0;
    if (!(result > 0.0) || !std::isfinite(result)) return 0.0;
    if (result > 1.0) result = 1.0;
    return result;
  }
  
  double integrate_ss_stoptrunc_Z(double LT,
                                  double UT,
                                  const double* go_pars,
                                  int n_go,
                                  int n_go_par,
                                  RacePdf1Fun go_pdf1,
                                  RaceCdf1Fun go_cdf1,
                                  void* go_ctx,
                                  double tf,
                                  double SSD,
                                  double muS,
                                  double sigmaS,
                                  double tauS,
                                  double epsilon,
                                  gsl_integration_workspace* w) {
    if (w == nullptr) Rcpp::stop("integrate_ss_stoptrunc_Z: GSL workspace is null.");
    if (n_go <= 0) return 0.0;
    if (LT >= UT && !(LT == 0.0 && UT == R_PosInf)) return 0.0;
    
    gsl_function F;
    gsl_ss_stoptrunc_params P;
    P.go_pars = go_pars;
    P.n_go = n_go;
    P.n_go_par = n_go_par;
    P.go_pdf1 = go_pdf1;
    P.go_cdf1 = go_cdf1;
    P.go_ctx = go_ctx;
    P.tf = tf;
    P.SSD = SSD;
    P.muS = muS;
    P.sigmaS = sigmaS;
    P.tauS = tauS;
    F.function = &gsl_f_ss_stoptrunc;
    F.params = &P;
    
    double low = LT;
    double upp = UT;
    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
    int status;
    double result = 0.0;
    double error = 0.0;
    if (upp == R_PosInf) {
      if (low < 0.0) low = 0.0;
      if (low >= R_PosInf) {
        result = 0.0;
        status = GSL_SUCCESS;
      } else {
        status = gsl_integration_qagiu(&F, low, 0, epsilon, 1000, w, &result, &error);
      }
    } else {
      status = gsl_integration_qags(&F, low, upp, 0, epsilon, 1000, w, &result, &error);
    }
    gsl_set_error_handler(old_handler);
    if (status != GSL_SUCCESS) return 0.0;
    if (!(result > 0.0) || !std::isfinite(result)) return 0.0;
    return result;
  }
  
} // namespace

double c_log_likelihood_ss_cens_trunc(
    const Rcpp::NumericMatrix& pars,
    const Rcpp::DataFrame& dadm,
    int n_trials,
    const Rcpp::LogicalVector& winner,
    const Rcpp::IntegerVector& expand,
    double min_ll,
    const Rcpp::LogicalVector& isok,
    int n_lR,
    const std::string& ss_type) {
  
  if (n_trials == 0) return 0.0;
  if (n_lR <= 0) Rcpp::stop("c_log_likelihood_ss_cens_trunc: n_lR must be positive.");
  if (n_trials % n_lR != 0) Rcpp::stop("c_log_likelihood_ss_cens_trunc: dadm nrows not a multiple of n_lR.");
  
  SsGoDist go_dist;
  if (ss_type == "SSexG") {
    go_dist = SsGoDist::EXG;
  } else if (ss_type == "SShybrid") {
    go_dist = SsGoDist::RDM;
  } else {
    Rcpp::stop("Unsupported stop-signal model type string in calc_ll: " + ss_type);
  }
  
  const Rcpp::CharacterVector pcols = Rcpp::colnames(pars);
  const int col_gf = find_col_idx(pcols, "gf");
  const int col_tf = find_col_idx(pcols, "tf");
  const int col_muS = find_col_idx(pcols, "muS");
  const int col_sigmaS = find_col_idx(pcols, "sigmaS");
  const int col_tauS = find_col_idx(pcols, "tauS");
  if (col_gf < 0 || col_tf < 0 || col_muS < 0 || col_sigmaS < 0 || col_tauS < 0) {
    Rcpp::stop("c_log_likelihood_ss_cens_trunc: required stop-signal parameter columns not found in pars.");
  }
  
  std::array<int, 5> go_cols{};
  int n_go_par = 0;
  RacePdf1Fun go_pdf1 = nullptr;
  RaceCdf1Fun go_cdf1 = nullptr;
  void* go_ctx = nullptr;
  
  if (go_dist == SsGoDist::EXG) {
    const int col_mu = find_col_idx(pcols, "mu");
    const int col_sigma = find_col_idx(pcols, "sigma");
    const int col_tau = find_col_idx(pcols, "tau");
    if (col_mu < 0 || col_sigma < 0 || col_tau < 0) {
      Rcpp::stop("c_log_likelihood_ss_cens_trunc: SSexG requires go columns mu/sigma/tau.");
    }
    go_cols = {col_mu, col_sigma, col_tau, -1, -1};
    n_go_par = 3;
    go_pdf1 = &ss_exg_go_pdf1;
    go_cdf1 = &ss_exg_go_cdf1;
  } else {
    const int col_v = find_col_idx(pcols, "v");
    const int col_B = find_col_idx(pcols, "B");
    const int col_A = find_col_idx(pcols, "A");
    const int col_t0 = find_col_idx(pcols, "t0");
    const int col_s = find_col_idx(pcols, "s");
    if (col_v < 0 || col_B < 0 || col_A < 0 || col_t0 < 0 || col_s < 0) {
      Rcpp::stop("c_log_likelihood_ss_cens_trunc: SShybrid requires go columns v/B/A/t0/s.");
    }
    go_cols = {col_v, col_B, col_A, col_t0, col_s};
    n_go_par = 5;
    go_pdf1 = &drdm_scalar;
    go_cdf1 = &prdm_scalar;
  }
  
  if (!dadm.containsElementNamed("SSD")) {
    Rcpp::stop("c_log_likelihood_ss_cens_trunc: dadm must contain SSD column.");
  }
  const Rcpp::NumericVector RT = dadm["rt"];
  const Rcpp::NumericVector SSDv = dadm["SSD"];
  const bool has_LT = dadm.containsElementNamed("LT");
  const bool has_UT = dadm.containsElementNamed("UT");
  const Rcpp::NumericVector LT = has_LT ? Rcpp::NumericVector(dadm["LT"]) : Rcpp::NumericVector();
  const Rcpp::NumericVector UT = has_UT ? Rcpp::NumericVector(dadm["UT"]) : Rcpp::NumericVector();
  const bool has_UC = dadm.containsElementNamed("UC");
  const Rcpp::NumericVector UC = has_UC ? Rcpp::NumericVector(dadm["UC"]) : Rcpp::NumericVector();
  
  // Determine which response accumulators are GO (lI==2) vs stop-triggered (ST; lI==1).
  std::vector<int> go_offsets;
  std::vector<int> st_offsets;
  if (dadm.containsElementNamed("lI")) {
    const Rcpp::IntegerVector lI = dadm["lI"];
    for (int k = 0; k < n_lR; ++k) {
      const int code = lI[k];
      if (code == 2) go_offsets.push_back(k);
      else if (code == 1) st_offsets.push_back(k);
      else {
        Rcpp::stop("c_log_likelihood_ss_cens_trunc: unexpected lI code (expected 1=ST, 2=GO).");
      }
    }
  } else {
    // No lI provided: treat all accumulators as GO and none as ST.
    go_offsets.reserve(static_cast<size_t>(n_lR));
    for (int k = 0; k < n_lR; ++k) go_offsets.push_back(k);
  }
  const int n_go = static_cast<int>(go_offsets.size());
  const int n_st = static_cast<int>(st_offsets.size());
  
  auto workspace = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>(
    gsl_integration_workspace_alloc(1000),
    &gsl_integration_workspace_free);
  if (!workspace) Rcpp::stop("c_log_likelihood_ss_cens_trunc: Failed to allocate GSL integration workspace.");
  
  const double epsilon = 1e-7;
  std::unordered_map<std::string, double> cache_prob;
  std::unordered_map<std::string, double> cache_logZ;
  
  const int n_unique = n_trials / n_lR;
  Rcpp::NumericVector ll_unique(n_unique, min_ll);
  
  auto build_resp_pars = [&](int start_row, std::vector<double>& resp_pars_out) {
    resp_pars_out.resize(static_cast<size_t>(n_lR) * static_cast<size_t>(n_go_par));
    for (int k = 0; k < n_lR; ++k) {
      const int row = start_row + k;
      for (int j = 0; j < n_go_par; ++j) {
        resp_pars_out[static_cast<size_t>(k) * n_go_par + j] =
          pars(row, go_cols[static_cast<size_t>(j)]);
      }
    }
  };
  
  auto build_goonly_pars = [&](const std::vector<double>& resp_pars_in,
                               std::vector<double>& goonly_pars_out) {
    goonly_pars_out.resize(static_cast<size_t>(n_go) * static_cast<size_t>(n_go_par));
    for (int i = 0; i < n_go; ++i) {
      const int k0 = go_offsets[static_cast<size_t>(i)];
      for (int j = 0; j < n_go_par; ++j) {
        goonly_pars_out[static_cast<size_t>(i) * n_go_par + j] =
          resp_pars_in[static_cast<size_t>(k0) * n_go_par + j];
      }
    }
  };
  
  auto log_surv_prod = [&](double t, const std::vector<double>& pars_pack, int n_acc) -> double {
    double logS = 0.0;
    for (int k = 0; k < n_acc; ++k) {
      const double* par_k = pars_pack.data() + static_cast<size_t>(k) * n_go_par;
      double cdf = go_cdf1(t, par_k, go_ctx);
      cdf = clamp_cdf01(cdf);
      const double ll = safe_log1m(cdf);
      if (!std::isfinite(ll)) return R_NegInf;
      logS += ll;
    }
    return logS;
  };
  
  auto log_surv_prod_subset = [&](double t,
                                  const std::vector<double>& resp_pars,
                                  const std::vector<int>& offsets) -> double {
                                    double logS = 0.0;
                                    for (int off : offsets) {
                                      const double* par_k = resp_pars.data() + static_cast<size_t>(off) * n_go_par;
                                      double cdf = go_cdf1(t, par_k, go_ctx);
                                      cdf = clamp_cdf01(cdf);
                                      const double ll = safe_log1m(cdf);
                                      if (!std::isfinite(ll)) return R_NegInf;
                                      logS += ll;
                                    }
                                    return logS;
                                  };
  
  auto surv_prod = [&](double t, const std::vector<double>& pars_pack, int n_acc) -> double {
    const double logS = log_surv_prod(t, pars_pack, n_acc);
    if (!std::isfinite(logS)) return 0.0;
    const double out = std::exp(logS);
    return (std::isfinite(out) ? out : 0.0);
  };
  
  auto surv_prod_subset = [&](double t,
                              const std::vector<double>& resp_pars,
                              const std::vector<int>& offsets) -> double {
                                const double logS = log_surv_prod_subset(t, resp_pars, offsets);
                                if (!std::isfinite(logS)) return 0.0;
                                const double out = std::exp(logS);
                                return (std::isfinite(out) ? out : 0.0);
                              };
  
  auto log_joint_for_winner_subset = [&](double t,
                                         const std::vector<double>& resp_pars,
                                         const std::vector<int>& offsets,
                                         int winner_off) -> double {
                                           const double* par_w = resp_pars.data() + static_cast<size_t>(winner_off) * n_go_par;
                                           double pdf_w = go_pdf1(t, par_w, go_ctx);
                                           if (!(pdf_w > 0.0) || !std::isfinite(pdf_w)) return R_NegInf;
                                           double logv = std::log(pdf_w);
                                           for (int off : offsets) {
                                             if (off == winner_off) continue;
                                             const double* par_k = resp_pars.data() + static_cast<size_t>(off) * n_go_par;
                                             double cdf = go_cdf1(t, par_k, go_ctx);
                                             cdf = clamp_cdf01(cdf);
                                             const double ll = safe_log1m(cdf);
                                             if (!std::isfinite(ll)) return R_NegInf;
                                             logv += ll;
                                           }
                                           return logv;
                                         };
  
  auto log_race_subset = [&](double t,
                             const std::vector<double>& resp_pars,
                             const std::vector<int>& offsets,
                             int observed_winner_idx0) -> double {
                               if (offsets.empty()) return R_NegInf;
                               if (observed_winner_idx0 >= 0) {
                                 // Only valid if the observed winner is in this subset.
                                 bool in_subset = false;
                                 for (int off : offsets) if (off == observed_winner_idx0) { in_subset = true; break; }
                                 if (!in_subset) return R_NegInf;
                                 return log_joint_for_winner_subset(t, resp_pars, offsets, observed_winner_idx0);
                               }
                               // Unknown winner: sum over subset.
                               double out = R_NegInf;
                               for (int off : offsets) {
                                 out = log_sum_exp(out, log_joint_for_winner_subset(t, resp_pars, offsets, off));
                               }
                               return out;
                             };
  
  // Cache keys using only the go-only parameter pack plus stop parameters.
  auto make_stopwin_key = [&](double upper_dur,
                              double SSD,
                              double muS,
                              double sigmaS,
                              double tauS,
                              const std::vector<double>& goonly_pars) -> std::string {
                                std::ostringstream oss;
                                oss.setf(std::ios::scientific);
                                oss << std::setprecision(17);
                                oss << "SSP|" << upper_dur << '|' << SSD << '|' << muS << '|' << sigmaS << '|' << tauS << '|';
                                for (double v : goonly_pars) {
                                  if (Rcpp::NumericVector::is_na(v)) oss << "NA";
                                  else oss << v;
                                  oss << ',';
                                }
                                return oss.str();
                              };
  
  auto make_stoptrunc_key = [&](double LTj,
                                double UTj,
                                double tf,
                                double SSD,
                                double muS,
                                double sigmaS,
                                double tauS,
                                const std::vector<double>& goonly_pars) -> std::string {
                                  std::ostringstream oss;
                                  oss.setf(std::ios::scientific);
                                  oss << std::setprecision(17);
                                  oss << "SSZ|" << LTj << '|' << UTj << '|' << tf << '|' << SSD << '|' << muS << '|' << sigmaS << '|' << tauS << '|';
                                  for (double v : goonly_pars) {
                                    if (Rcpp::NumericVector::is_na(v)) oss << "NA";
                                    else oss << v;
                                    oss << ',';
                                  }
                                  return oss.str();
                                };
  
  auto make_stoptrunc_key_st = [&](double LTj,
                                   double UTj,
                                   double gf,
                                   double tf,
                                   double SSD,
                                   double muS,
                                   double sigmaS,
                                   double tauS,
                                   const std::vector<double>& resp_pars) -> std::string {
                                     std::ostringstream oss;
                                     oss.setf(std::ios::scientific);
                                     oss << std::setprecision(17);
                                     oss << "SSZST|" << LTj << '|' << UTj << '|' << gf << '|' << tf << '|' << SSD << '|' << muS << '|' << sigmaS << '|' << tauS << '|';
                                     for (double v : resp_pars) {
                                       if (Rcpp::NumericVector::is_na(v)) oss << "NA";
                                       else oss << v;
                                       oss << ',';
                                     }
                                     return oss.str();
                                   };
  
  for (int j = 0; j < n_unique; ++j) {
    const int start_row = j * n_lR;
    if (!isok[start_row]) {
      ll_unique[j] = min_ll;
      continue;
    }
    
    const double rt_j = RT[start_row];
    const double SSD = SSDv[start_row];
    const bool is_stop = R_finite(SSD);
    const double gf = pars(start_row, col_gf);
    const double tf = pars(start_row, col_tf);
    const double muS = pars(start_row, col_muS);
    const double sigmaS = pars(start_row, col_sigmaS);
    const double tauS = pars(start_row, col_tauS);
    
    std::vector<double> resp_pars;
    build_resp_pars(start_row, resp_pars);
    std::vector<double> goonly_pars;
    build_goonly_pars(resp_pars, goonly_pars);
    
    if (rt_j == R_PosInf) {
      if (!has_UC) Rcpp::stop("c_log_likelihood_ss_cens_trunc: rt contains Inf but dadm has no UC column.");
      const double UCj = UC[start_row];
      const double S_go_UC = (n_go > 0) ? surv_prod(UCj, goonly_pars, n_go) : 1.0;
      
      double p_nr = 0.0;
      if (!is_stop) {
        p_nr = gf + (1.0 - gf) * S_go_UC;
      } else {
        const double S_stop_UC = clamp_prob01(ss_stop_surv_abs(UCj, SSD, muS, sigmaS, tauS));
        const double upper_dur = UCj - SSD;
        
        double pStop_UC = 0.0;
        if (n_go > 0 && (1.0 - gf) > 0.0) {
          const std::string key = make_stopwin_key(upper_dur, SSD, muS, sigmaS, tauS, goonly_pars);
          auto hit = cache_prob.find(key);
          if (hit != cache_prob.end()) {
            pStop_UC = hit->second;
          } else {
            pStop_UC = integrate_ss_stopwin_prob(upper_dur,
                                                 goonly_pars.data(), n_go, n_go_par,
                                                 go_cdf1, go_ctx,
                                                 SSD, muS, sigmaS, tauS,
                                                 epsilon, workspace.get());
            cache_prob.emplace(key, pStop_UC);
          }
        } else if (n_go == 0) {
          pStop_UC = 1.0;
        }
        pStop_UC = clamp_prob01(pStop_UC);
        
        // ST survivor by UC if stop triggers (deadline precedes trigger => ST cannot respond yet).
        double S_st_UC = 1.0;
        if (n_st > 0) {
          double t_eff = UCj - SSD;
          if (!R_finite(t_eff) || t_eff <= 0.0) t_eff = 0.0;
          S_st_UC = surv_prod_subset(t_eff, resp_pars, st_offsets);
        }
        
        const double p_nr_if_tf = gf + (1.0 - gf) * S_go_UC;
        const double p_nr_if_trig = S_st_UC * (gf + (1.0 - gf) * (pStop_UC + S_go_UC * S_stop_UC));
        p_nr = tf * p_nr_if_tf + (1.0 - tf) * p_nr_if_trig;
      }
      
      ll_unique[j] = (p_nr > 0.0 && std::isfinite(p_nr)) ? std::max(min_ll, std::log(p_nr)) : min_ll;
      continue;
    }
    
    if (Rcpp::NumericVector::is_na(rt_j)) {
      double p_nr = 0.0;
      if (!is_stop) {
        p_nr = gf;
      } else {
        if (n_st > 0) {
          // With stop-triggered accumulators, a true intrinsic non-response requires both go failure and trigger failure.
          p_nr = gf * tf;
        } else {
          const std::string key = make_stopwin_key(R_PosInf, SSD, muS, sigmaS, tauS, goonly_pars);
          auto hit = cache_prob.find(key);
          double pStop = 0.0;
          if (hit != cache_prob.end()) {
            pStop = hit->second;
          } else {
            pStop = integrate_ss_stopwin_prob(R_PosInf,
                                              goonly_pars.data(), n_go, n_go_par,
                                              go_cdf1, go_ctx,
                                              SSD, muS, sigmaS, tauS,
                                              epsilon, workspace.get());
            cache_prob.emplace(key, pStop);
          }
          pStop = clamp_prob01(pStop);
          p_nr = gf + (1.0 - gf) * (1.0 - tf) * pStop;
        }
      }
      ll_unique[j] = (p_nr > 0.0 && std::isfinite(p_nr)) ? std::max(min_ll, std::log(p_nr)) : min_ll;
      continue;
    }
    
    if (!R_FINITE(rt_j) || !(rt_j > 0.0)) {
      ll_unique[j] = min_ll;
      continue;
    }
    
    const double LTj = has_LT ? LT[start_row] : 0.0;
    const double UTj = has_UT ? UT[start_row] : R_PosInf;
    if (rt_j < LTj || rt_j > UTj) {
      ll_unique[j] = min_ll;
      continue;
    }
    
    int observed_winner_idx0 = -1;
    int n_true = 0;
    for (int k = 0; k < n_lR; ++k) {
      if (winner[start_row + k] == TRUE) {
        n_true++;
        observed_winner_idx0 = k;
      }
    }
    if (n_true != 1) observed_winner_idx0 = -1; // unknown
    
    double log_like = R_NegInf;
    if (!is_stop) {
      // GO trial: only GO accumulators race.
      const double log_race_go = log_race_subset(rt_j, resp_pars, go_offsets, observed_winner_idx0);
      if (std::isfinite(log_race_go) && (1.0 - gf) > 0.0) {
        log_like = std::log1p(-gf) + log_race_go;
      }
    } else {
      // Stop trial: GO responses are a mixture over trigger failure; ST responses require triggering.
      bool winner_is_go = false;
      bool winner_is_st = false;
      if (observed_winner_idx0 >= 0) {
        for (int off : go_offsets) if (off == observed_winner_idx0) { winner_is_go = true; break; }
        for (int off : st_offsets) if (off == observed_winner_idx0) { winner_is_st = true; break; }
      }
      
      // If we can't classify the winner, fall back to summing GO+ST response masses.
      const bool unknown_winner = (observed_winner_idx0 < 0);
      
      // GO response mass (only GO racers, at absolute rt).
      double log_go_mass = R_NegInf;
      if (n_go > 0 && (unknown_winner || winner_is_go)) {
        const double log_race_go = log_race_subset(rt_j, resp_pars, go_offsets, unknown_winner ? -1 : observed_winner_idx0);
        if (std::isfinite(log_race_go) && (1.0 - gf) > 0.0) {
          const double log_base = std::log1p(-gf) + log_race_go;
          // tf branch: stop doesn't trigger => GO only
          if (tf > 0.0) log_go_mass = log_sum_exp(log_go_mass, std::log(tf) + log_base);
          // triggered branch: stop triggers => require stop survivor and ST survivors
          if ((1.0 - tf) > 0.0) {
            const double S_stop_rt = clamp_prob01(ss_stop_surv_abs(rt_j, SSD, muS, sigmaS, tauS));
            double logS_stop = (S_stop_rt > 0.0) ? std::log(S_stop_rt) : R_NegInf;
            double logS_st = 0.0;
            if (n_st > 0) {
              double t_eff = rt_j - SSD;
              if (!R_finite(t_eff) || t_eff <= 0.0) t_eff = 0.0;
              logS_st = log_surv_prod_subset(t_eff, resp_pars, st_offsets);
            }
            if (std::isfinite(logS_stop) && std::isfinite(logS_st)) {
              log_go_mass = log_sum_exp(log_go_mass, std::log1p(-tf) + log_base + logS_stop + logS_st);
            }
          }
        }
      }
      
      // ST response mass (only if stop triggers).
      double log_st_mass = R_NegInf;
      if (n_st > 0 && (unknown_winner || winner_is_st) && (1.0 - tf) > 0.0) {
        const double t_eff = rt_j - SSD;
        const double log_race_st = log_race_subset(t_eff, resp_pars, st_offsets, unknown_winner ? -1 : observed_winner_idx0);
        if (std::isfinite(log_race_st)) {
          // pStop up to rt (duration upper = rt-SSD), depends only on GO racers.
          double pStop_rt = 0.0;
          if (n_go == 0) {
            pStop_rt = 1.0;
          } else if ((1.0 - gf) > 0.0) {
            const double upper_dur = rt_j - SSD;
            const std::string key = make_stopwin_key(upper_dur, SSD, muS, sigmaS, tauS, goonly_pars);
            auto hit = cache_prob.find(key);
            if (hit != cache_prob.end()) {
              pStop_rt = hit->second;
            } else {
              pStop_rt = integrate_ss_stopwin_prob(upper_dur,
                                                   goonly_pars.data(), n_go, n_go_par,
                                                   go_cdf1, go_ctx,
                                                   SSD, muS, sigmaS, tauS,
                                                   epsilon, workspace.get());
              cache_prob.emplace(key, pStop_rt);
            }
          }
          pStop_rt = clamp_prob01(pStop_rt);
          const double log_pStop = (pStop_rt > 0.0) ? std::log(pStop_rt) : R_NegInf;
          const double log1m_pStop = safe_log1m(pStop_rt);
          
          // GO survivor by rt (only matters when stop doesn't beat go).
          const double logS_go_rt = (n_go > 0) ? log_surv_prod(rt_j, goonly_pars, n_go) : 0.0;
          
          // gf branch: only ST races
          if (gf > 0.0) {
            log_st_mass = log_sum_exp(log_st_mass, std::log(gf) + log_race_st);
          }
          // (1-gf) branch: mixture over pStop vs go surviving
          if ((1.0 - gf) > 0.0) {
            double log_mix = log_pStop;
            if (std::isfinite(log1m_pStop) && std::isfinite(logS_go_rt)) {
              log_mix = log_sum_exp(log_mix, log1m_pStop + logS_go_rt);
            }
            if (std::isfinite(log_mix)) {
              log_st_mass = log_sum_exp(log_st_mass, std::log1p(-gf) + log_race_st + log_mix);
            }
          }
          if (std::isfinite(log_st_mass)) {
            log_st_mass += std::log1p(-tf);
          }
        }
      }
      
      log_like = log_sum_exp(log_go_mass, log_st_mass);
    }
    
    if (LTj != 0.0 || UTj != R_PosInf) {
      double logZ = R_NegInf;
      if (!is_stop) {
        const double logS_lt = (n_go > 0) ? log_surv_prod(LTj, goonly_pars, n_go) : 0.0;
        const double logS_ut = (n_go > 0)
          ? ((UTj == R_PosInf) ? R_NegInf : log_surv_prod(UTj, goonly_pars, n_go))
          : 0.0;
        double log_diff = log_diff_exp(logS_lt, logS_ut);
        if (!std::isfinite(log_diff)) log_diff = std::log(std::numeric_limits<double>::min());
        logZ = ((1.0 - gf) > 0.0) ? (std::log1p(-gf) + log_diff) : R_NegInf;
      } else {
        if (n_st > 0) {
          // Z = P(response in (LT, UT)) = P(no response by LT) - P(no response by UT).
          // Requires pStop(t) evaluation (stop-win CDF) at LT and UT; cached via make_stopwin_key.
          const std::string key = make_stoptrunc_key_st(LTj, UTj, gf, tf, SSD, muS, sigmaS, tauS, resp_pars);
          auto hit = cache_logZ.find(key);
          if (hit != cache_logZ.end()) {
            logZ = hit->second;
          } else {
            auto p_no_by_t = [&](double t_abs) -> double {
              // GO survivor by absolute time.
              const double S_go = (n_go > 0) ? surv_prod(t_abs, goonly_pars, n_go) : 1.0;
              const double p_tf = gf + (1.0 - gf) * S_go;
              
              // Stop survivor by absolute time.
              double S_stop = ss_stop_surv_abs(t_abs, SSD, muS, sigmaS, tauS);
              if (!std::isfinite(S_stop)) S_stop = 0.0;
              S_stop = clamp_prob01(S_stop);
              
              // pStop(t): P(stop finishes before GO by t).
              double pStop = 0.0;
              if (n_go == 0) {
                pStop = 1.0;
              } else if ((1.0 - gf) > 0.0) {
                const double upper_dur = t_abs - SSD;
                const std::string keyP = make_stopwin_key(upper_dur, SSD, muS, sigmaS, tauS, goonly_pars);
                auto hitP = cache_prob.find(keyP);
                if (hitP != cache_prob.end()) {
                  pStop = hitP->second;
                } else {
                  pStop = integrate_ss_stopwin_prob(upper_dur,
                                                    goonly_pars.data(), n_go, n_go_par,
                                                    go_cdf1, go_ctx,
                                                    SSD, muS, sigmaS, tauS,
                                                    epsilon, workspace.get());
                  cache_prob.emplace(keyP, pStop);
                }
              }
              pStop = clamp_prob01(pStop);
              
              // ST survivor by absolute time (shifted by SSD; deadline before trigger => ST cannot respond yet).
              double S_st = 1.0;
              if (n_st > 0) {
                double t_eff = t_abs - SSD;
                if (!R_finite(t_eff) || t_eff <= 0.0) t_eff = 0.0;
                S_st = surv_prod_subset(t_eff, resp_pars, st_offsets);
              }
              
              const double p_trig_core = gf + (1.0 - gf) * (pStop + S_go * S_stop);
              const double p_trig = S_st * p_trig_core;
              
              const double p = tf * p_tf + (1.0 - tf) * p_trig;
              return clamp_prob01(p);
            };
            
            const double p_lt = p_no_by_t(LTj);
            const double p_ut = p_no_by_t(UTj);
            if (!(p_lt > p_ut) || !std::isfinite(p_lt) || !std::isfinite(p_ut)) {
              logZ = std::log(std::numeric_limits<double>::min());
            } else {
              double logZtmp = log_diff_exp(std::log(p_lt), std::log(p_ut));
              if (!std::isfinite(logZtmp)) logZtmp = std::log(std::numeric_limits<double>::min());
              logZ = logZtmp;
            }
            cache_logZ.emplace(key, logZ);
          }
        } else {
          const std::string key = make_stoptrunc_key(LTj, UTj, tf, SSD, muS, sigmaS, tauS, goonly_pars);
          auto hit = cache_logZ.find(key);
          if (hit != cache_logZ.end()) {
            logZ = hit->second;
          } else {
            const double Zcore = integrate_ss_stoptrunc_Z(LTj, UTj,
                                                          goonly_pars.data(), n_go, n_go_par,
                                                          go_pdf1, go_cdf1, go_ctx,
                                                          tf, SSD, muS, sigmaS, tauS,
                                                          epsilon, workspace.get());
            const double Z = (1.0 - gf) * std::max(Zcore, std::numeric_limits<double>::min());
            logZ = (Z > 0.0) ? std::log(Z) : R_NegInf;
            cache_logZ.emplace(key, logZ);
          }
        }
      }
      
      if (!std::isfinite(logZ)) {
        ll_unique[j] = min_ll;
        continue;
      }
      log_like -= logZ;
    }
    
    ll_unique[j] = std::isfinite(log_like) ? std::max(min_ll, log_like) : min_ll;
  }
  
  if (expand.length() > 0) {
    return sum(c_expand(ll_unique, expand));
  }
  double total = 0.0;
  for (int i = 0; i < ll_unique.size(); ++i) total += ll_unique[i];
  return total;
}
