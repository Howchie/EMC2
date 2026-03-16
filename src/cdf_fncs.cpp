
// Chair of Social Psychology, University of Freiburg
// Authors: Christoph Klauer and Raphael Hartmann

//#include "cstdio"
#include "cdf_fncs.h"
#include "tools.h"
//#include <cmath>
//#include <limits>
//#include <Rinternals.h>


/* DISTRIBUTION */

/* P term in the distribution function */
double logP(int pm, double a, double v, double w) {
	const double em1 = 1.0-1.0e-6;
	if (pm == 1) { v = -v; w = 1.0 - w; }
	if (fabs(v) == 0.0) return log1p(-w);
	// if (fabs(v) < 1.0e-16) return log1p(-w);
	
	double prob;
	double e = (-2.0 * v * a * (1.0 - w));
	if (e < 0) {
		double tt = exp(e);
	  if (tt >= em1) return log1p(-w);
		tt = log1p(-tt) - logdiff(2 * v * a * w, e);
		prob = tt;
	}
	else {
		double tt = exp(-e);
	  if (tt >= em1) return log1p(-w);
	  tt = log1p(-tt) - log1p(-exp(2 * v * a));
		prob = tt;
	}
	return prob;
}

/* calculate number of terms needed for short t */
double Ks(double t, double v, double a, double w, double eps)
{
	double K1 = 0.5 * (fabs(v) / a * t - w);
	double arg = fmax(0, fmin(1, exp(v*a*w + v * v*t / 2 + (eps)) / 2));
	double K2 = (arg==0) ? INFINITY : (arg==1) ? -INFINITY : -sqrt(t) / 2 / a * emc2_gsl_cdf_ugaussian_Pinv(arg);
	return ceil(fmax(K1, K1 + K2));
}

/* calculate number of terms needed for large t */
double Kl(double t, double v, double a, double w, double err) {
	double api = a / M_PI, vsq = v * v;
	double sqrtL1 = sqrt(1 / t) * api;
	double sqrtL2 = sqrt(fmax(1.0, -2 / t * api * api * (err+std::log(M_PI*t / 2 * (vsq + (M_PI / a) * (M_PI / a))) + v * a * w + vsq * t / 2)));
	return ceil(fmax(sqrtL1, sqrtL2));
}

/* calculate terms of the sum for short t */
double logFs(double t, double v, double a, double w, int K)
{
	double fplus = -INFINITY, fminus = -INFINITY;
	double sqt = sqrt(t), temp = -v * a*w - v * v*t / 2;
	double vt = v * t;

	for (int k = K; k >= 0; k--)
	{
		double rj = a*(2 * k + w);
		double dj = lognormal(rj / sqt);
		double pos1 = dj + logMill((rj - vt) / sqt);
		double pos2 = dj + logMill((rj + vt) / sqt);
		fplus = logsum(logsum(pos1, pos2), fplus);
		rj = a*(2.0 * k + 2.0 - w);
		dj =  lognormal(rj / sqt);
		double neg1 = dj + logMill((rj - vt) / sqt);
		double neg2 = dj + logMill((rj + vt) / sqt);
		fminus = logsum(logsum(neg1, neg2), fminus);
	}

	return logdiff(fplus, fminus)+temp;
}

/* calculate terms of the sum for large t */
double logFl(double q, double v, double a, double w, int K)
{
	double fplus = -INFINITY, fminus = -INFINITY;
	double la = std::log(a), lv = std::log(fabs(v));
	double F = -INFINITY;
	for (int k = K; k >= 1; k--) {
		double temp0 = std::log(k * 1.0), temp1 = k * M_PI, temp2 = temp1 * w;
		double check = sin(temp2);
		if (check > 0) {
			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * (temp1 / a) * (temp1 / a) * q + std::log(check);
			fplus = logsum(temp, fplus);
		}
		else if (check < 0)
		{
			double temp = temp0 - logsum(2 * lv, 2 * (temp0 + M_LNPI - la)) - 0.5 * (temp1 / a) * (temp1 / a) * q + std::log(-check);
			fminus = logsum(temp, fminus);
		}
	}
	F = logdiff(fplus, fminus);
	return (F - v * a * w - 0.5 * v * v * q);
}

/* calculate distribution */
double pwiener(double q, double a, double v, double w, double err, int K, int epsFLAG) {
	//if (q == 0) return(GSL_NEGINF);
	double Kll, Kss, ans;
	if(!epsFLAG && K==0) {
		err = -27.63102; // exp(err) = 1.e-12
		epsFLAG = 1;
	}
	else if(!epsFLAG && K>0) err = -27.63102; // exp(err) = 1.e-12
	else if(epsFLAG) err = std::log(err);

	if (std::isinf(q)) return logP(0, a, v, w);

	Kss = Ks(q, v, a, w, err);
	Kll = Kl(q, v, a, w, err);
	double lg = M_LN2 + M_LNPI - 2.0 * std::log(a);

  if (3 * Kss < Kll /*|| fabs(v) < 1.0e-1*/) {
		if((epsFLAG && Kss<K) || !epsFLAG) Kss = K;
		ans = logFs(q, v, a, w, static_cast<int>(Kss));
		//Rprintf("short hier\n");
	}
	else {
		if((epsFLAG && Kll<K) || !epsFLAG) Kll = K;
		ans = logdiff(logP(0, a, v, w), lg + logFl(q, v, a, w, static_cast<int>(Kll)));
		//Rprintf("large\n");
	}
	return ans;
}