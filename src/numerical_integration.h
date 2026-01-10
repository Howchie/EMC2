#ifndef numerical_integration_h
#define numerical_integration_h

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
using namespace Numer;

class race_f : public Func {
private:
  Rcpp::NumericMatrix pars;
  Rcpp::LogicalVector winner;
  Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector);
  Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector);
  double min_ll;
  Rcpp::LogicalVector is_ok;

public:
  race_f(Rcpp::NumericMatrix pars_,
        Rcpp::LogicalVector winner_,
        Rcpp::NumericVector (*dfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
        Rcpp::NumericVector (*pfun_)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
        double min_ll_,
        Rcpp::LogicalVector is_ok_) :
  pars(pars_),
  winner(winner_),
  dfun(dfun_),
  pfun(pfun_),
  min_ll(min_ll_),
  is_ok(is_ok_){}

  double operator()(const double& x) const
  {
    double accumulators = pars.nrow();
    Rcpp::NumericVector t(accumulators);
    t.fill(x);
    Rcpp::NumericVector d = dfun(t, pars, winner, min_ll, is_ok);
    double out = Rcpp::as<double>(d);

    if (accumulators > 1) {
      Rcpp::NumericVector p = 1 - pfun(t, pars, !winner, min_ll, is_ok);
      double prod_p = std::accumulate(p.begin(), p.end(), 1.0,
                                      std::multiplies<double>());
      out *= prod_p;
    }
    return out;
  }
};

Rcpp::NumericVector f_integrate(Rcpp::NumericMatrix pars,
                   Rcpp::LogicalVector winner,
                   Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
                   Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
                   double min_ll,
                   Rcpp::LogicalVector is_ok,
                   double lower,
                   double upper)
{

  if (R_finite(upper) && lower >= upper) {
    return Rcpp::NumericVector::create(0.0, 0.0, 0.0);
  }

  race_f f(pars, winner, dfun, pfun, min_ll, is_ok);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  Rcpp::NumericVector out{res, err_est, (double) err_code};
  return out;
}

Rcpp::NumericVector f_integrate_slow(Rcpp::NumericMatrix pars,
                               Rcpp::LogicalVector winner,
                               Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
                               Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
                               double min_ll,
                               Rcpp::LogicalVector is_ok,
                               double lower,
                               double upper)
{

  if (R_finite(upper) && lower >= upper) {
    return Rcpp::NumericVector::create(0.0, 0.0, 0.0);
  }

  race_f f(pars, winner, dfun, pfun, min_ll, is_ok);
  double err_est;
  int err_code;
  double res = integrate(f, lower, upper, err_est, err_code);
  if (err_code == 1 && upper == R_PosInf) {
    double err_est_hacky;
    int err_code_hacky;
    double res_hacky = integrate(f, lower, 10, err_est_hacky, err_code_hacky);
    Rcpp::NumericVector out{res_hacky, err_est_hacky, (double) err_code_hacky};
    return out;
  } else if (err_code > 0) {
    Rcpp::NumericVector out{min_ll, err_est, (double) err_code};
    return out;
  } else {
    Rcpp::NumericVector out{res, err_est, (double) err_code};
    return out;
  }
}


#include <Rcpp.h>
using namespace Rcpp;

inline double clamp01(double x){
  if (!R_finite(x)) return NA_REAL;
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0 - 1e-15;
  return x;
}

inline double surv_prod(const NumericVector& cdf){
  double prod = 1.0;
  for (int i = 0; i < cdf.size(); ++i){
    if (NumericVector::is_na(cdf[i])) return NA_REAL;
    const double c = clamp01(cdf[i]);
    if (!R_finite(c)) return NA_REAL;
    prod *= (1.0 - c);
    if (!(prod > 0.0)) return 0.0; // early zero
  }
  return prod;
}

// pr_pt: 1 / | prod(1-CDF(LT)) - prod(1-CDF(UT)) |, NA if undefined/zero
double pr_pt(NumericMatrix pars,
             NumericVector (*pfun)(NumericVector, NumericMatrix, LogicalVector, double, LogicalVector),
             double LT,
             double UT)
{
  const int n = pars.nrow();
  LogicalVector idx(n, true);     // use all rows/accumulators
  LogicalVector is_ok(n, true);   // assume all valid
  const double min_lik = 0.0;     // not used here; pass 0

  NumericVector cLT = pfun(NumericVector::create(LT), pars, idx, min_lik, is_ok);

  const double Slt = surv_prod(cLT);
  if (NumericVector::is_na(Slt) ) return NA_REAL;
  if (!R_finite(UT)) return 1/Slt;

  NumericVector cUT = pfun(NumericVector::create(UT), pars, idx, min_lik, is_ok);
  const double Sut = surv_prod(cUT);
  if (NumericVector::is_na(Slt) || NumericVector::is_na(Sut)) return NA_REAL;

  const double diff = std::fabs(Slt - Sut);
  if (!(diff > 0.0) || !R_finite(diff)) return NA_REAL;

  return 1.0 / diff;
}

double pLU(Rcpp::NumericMatrix pars,
           Rcpp::LogicalVector winner,
           Rcpp::NumericVector (*dfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
           Rcpp::NumericVector (*pfun)(Rcpp::NumericVector, Rcpp::NumericMatrix, Rcpp::LogicalVector, double, Rcpp::LogicalVector),
           double min_ll,
           Rcpp::LogicalVector is_ok,
           double LT,
           double LC,
           double UC,
           double UT) {
  Rcpp::NumericVector pL = f_integrate(pars, winner, dfun, pfun, min_ll, is_ok, LT, LC);
  if ((pL[2] != 0) || Rcpp::traits::is_nan<REALSXP>(pL[0])) {
    return NA_REAL;
  }
  Rcpp::NumericVector pU = f_integrate(pars, winner, dfun, pfun, min_ll, is_ok, UC, UT);
  if ((pU[2] != 0) || Rcpp::traits::is_nan<REALSXP>(pU[0])) {
    return NA_REAL;
  }
  double out = std::max(0.0, std::min(pL[0], 1.0)) + std::max(0.0, std::min(pU[0], 1.0));
  return out;
}

#endif
