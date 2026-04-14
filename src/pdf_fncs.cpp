
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

#include "tools.h"
#include "pdf_fncs.h"
#include "ddm_functions_inline.h"
#include <cmath>


/* DENSITY */

/* calculate number of terms needed for short t */
double ks(double t, double w, double eps) {
    return ddm_ks(t, w, eps);
}

/* calculate number of terms needed for large t */
double kl(double q, double v, double w, double err) {
    return ddm_kl(q, v, w, err);
}

/* calculate terms of the sum for short t */
double logfs(double t, double w, int K) {
    // Note: the inline version uses the 'linear' variant for performance
    // but the old logfs used logsum/logdiff. For backward compatibility
    // with any callers of this non-inline function, we just route to 
    // the ddm_logfs_linear equivalent.
    return ddm_logfs_linear(t, w, K);
}

double logfs_linear(double t, double w, int K) {
    return ddm_logfs_linear(t, w, K);
}

/* calculate terms of the sum for large t */
double logfl(double q, double v, double w, int K) {
    return ddm_logfl_linear(q, v, w, K);
}

double logfl_linear(double q, double v, double w, int K) {
    return ddm_logfl_linear(q, v, w, K);
}

/* calculate density */
double dwiener(double q, double a, double vn, double wn, double sv, double err, int K, int epsFLAG) {
    return dwiener_inline(q, a, vn, wn, sv, err, K, epsFLAG);
}
