
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include <vector>
// #include <mutex>
#include <R.h>
#include <gsl/gsl_errno.h>

//
// extern std::mutex mtx_R_CUI; // mutex for R_CheckUserInterrupt


/*------------------------------------------------*/

/* gsl_sf_result_struct and cheb_series_struct
  are copied from the GNU scientific library version 2.6 */
struct emc2_gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct emc2_gsl_sf_result_struct emc2_gsl_sf_result;

struct emc2_cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct emc2_cheb_series_struct emc2_cheb_series;


/*------------------------------------------------*/
struct my_params {
	double t;
	int low_or_up;
	double a;
	double v;
	double t0;
	double w;
	double sw;
	double sv;
	double st;
	double errorW;
	int K;
	int epsFLAG;
	double *val_ptr;
};

struct point {
	double x = 0;
	double h = 0;
	double dh = 0;
};

struct piece {
	double z = 0;
	double slope = 0;
	double absc = 0;
	double center = 0;
};

struct ars_archiv {
	std::vector<point> hstore;
	std::vector<piece> lowerstore;
	std::vector<piece> upperstore;
	double startstore;
	double scalestore;
	double normstore;
	std::vector<double> sstore;
};


/*------------------------------------------------*/

/*MATHEMATICAL CONSTANTS*/

#ifndef M_PISQ
#define M_PISQ        9.86960440108935861883449099987615113531
#endif

#ifndef LNNORM_MAX_X
#define LNNORM_MAX_X  38.0
#endif

#ifndef LNNORM_MIN_X
#define LNNORM_MIN_X  -1.00e9
#endif

#ifndef SQR2PI
#define SQR2PI        2.506628274631000502415765e0
#endif
  /*-------------------------*/

  /* The following mathematical constants
  are copied from GSL version 2.6 */

  // Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough. GPL 3

#ifndef GSL_DBL_EPSILON
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#endif

#ifndef M_LNPI
#define M_LNPI        1.14472988584940017414342735135      /* ln(pi) */
#endif
  /*-------------------------*/

  /* The following mathematical constants
  are copied from Rmath.h */

  // Copyright (C) 1998-2011  The R Core Team. GPL 2.1
  // Copyright (C) 2004       The R Foundation. GPL 2.1

#ifndef M_PI
#define M_PI          3.141592653589793238462643383280        /* pi */
#endif

#ifndef M_LN2
#define M_LN2         0.693147180559945309417232121458        /* ln(2) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI  0.572364942924700087071713675677        /* log(sqrt(pi)) == log(pi)/2 */
#endif
  /*-------------------------*/
/*------------------------------------------------*/



/* USEFUL FUNCTIONS */
double lnnorm(double);
double logsum(double, double);
double logdiff(double, double);
double rexp(double);
double lognormal(double);
double logMill(double);
void R_CheckUserInterruptGuarded();

/*------------------------------------------------*/



/* from the GNU scientific library version 2.6 */
double rat_eval(const double, const size_t, const double, const size_t, const double);
double small(double);
double intermediate(double);
double tail(double);
double emc2_gsl_cdf_ugaussian_Pinv(const double);

/* Forward declaration of our custom error reporter (avoids duplicate symbol with libgsl) */
void emc2_gsl_error(const char * reason, const char * file,
                    int line, int gsl_errno);

// Helpers to control libgsl error handling (print via Rprintf, no abort).
gsl_error_handler_t* emc2_gsl_set_error_handler_rprintf();
gsl_error_handler_t* emc2_gsl_set_error_handler_off();

/* EMC2_GSL_ERROR_VAL: call our error reporter, and return the given value */
#define EMC2_GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       emc2_gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)


/* modified gsl_ran_gaussian */
double onenorm();
/* modified gsl_ran_gaussian_tail */
double emc2_gsl_ran_gaussian_tail(const double, const double);
/* modified gsl_ran_ugaussian_tail */
double emc2_gsl_ran_ugaussian_tail(const double a);
int emc2_gsl_sf_erfc_e(double, emc2_gsl_sf_result *);
double emc2_gsl_sf_erfc(double);

/* gsl_error_handler_t comes from <gsl/gsl_errno.h> */

#define EVAL_RESULT(fn) \
   emc2_gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     EMC2_GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;

#define EVAL_DOUBLE(fn) \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     EMC2_GSL_ERROR_VAL(#fn, status, result); \
   } ; \
   return result;


/* -------------------------------------------------------- */




double oneuni();
double oneuniL();
double phi(double, double, double);
double lower_bound_var(double, double, double);
double coth(double);
double lower_bound_time(double, double, double);
double exp_mean(int, double, double, double);

int is_interruption();



#endif
