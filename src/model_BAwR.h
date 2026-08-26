#ifndef EMC2_MODEL_BAWR_H
#define EMC2_MODEL_BAWR_H

// ---------------------------------------------------------------------------
// BAwR: ballistic accumulation with a ramping clearance rate.
//
//   dX/du = V - kappa u^p,   p > 0,   z ~ U(0, A),   b = B + A
//   X(u)  = z + V u - kappa u^(p+1) / (p + 1)
//
// The momentary drive starts at the launch strength V and deteriorates
// deterministically with PHYSICAL time.  Unlike BAwL (leak, -k X) and BAwF
// (global fading, X = h(u)[z + V u]) there is no state dependence whatsoever
// and, unlike BAwD, no second opposing process with its own evidence scale:
// the whole decay law is one coefficient and one exponent.  Math/bawd.tex
// writes the coefficient `a`; it is `kappa` here so that formulas cannot
// confuse it with the start-point range `A`.
//
// Because each trial loses drive at the same ABSOLUTE rate, a strong launch
// stays productive longer -- the trajectory peaks at u = (V/kappa)^(1/p) --
// which is the mechanism the model is for: weak sensory representations
// expire quickly, strong ones keep supporting accumulation.
//
// The crossing law X(u) = b is
//
//   V*(u, z) = (b - z) / u + kappa u^p / (p + 1),
//
// U-shaped in u, so weak launches intrinsically omit.  Tangency (dV*/du = 0)
// is at u_sat(z) = [(p+1)(b-z) / (kappa p)]^(1/(p+1)), giving the hard right
// endpoint and the critical launch
//
//   T_max = [(p+1) b / (kappa p)]^(1/(p+1)),    V_c(z) = kappa u_sat(z)^p.
//
// NOTE the contrast with BAwF, whose endpoint is free of b.  Here T_max grows
// with b: the decay is in physical time and knows nothing about the
// threshold, so a more cautious accumulator simply gets longer before the
// drive expires.  That keeps B a pure caution parameter (raising it moves the
// deadline out) at the cost of the b-free endpoint.  kappa is the
// caution-free property of the stimulus representation; T_max is derived and
// reported by the R-side Ttransform.
//
// Both branches of the likelihood are closed form.  V* is affine in z, so the
// live branch is the usual stop-loss reduction, and the density is
//
//   A f(u) = int_{w_lo}^{w_hi} (w - kappa u^p) g(w) dw
//
// -- structurally identical to BAwF's with the affine offset c = kappa u^p in
// place of k b H'(k u), which is why this header reuses that machinery rather
// than restating it.  The frozen branch has |dz/dw| = kappa^(-1/p) w^(1/p), a
// positive-power survivor integral: closed for the lognormal launch at any p
// via log_lognormal_power_stoploss() with m = -(p+1)/p, and evaluated by the
// shared positive-integrand quadrature for the normal launch.
//
// Included by model_LBA.h after model_BAwF.h; reuses bawd_log_gl_split() and
// wald_functions.h's lognormal stop-loss primitives.
// ---------------------------------------------------------------------------

#include <cmath>
#include "model_BAwF.h"

constexpr double BAWR_K_EPS = BAWD_K_EPS;
constexpr double BAWR_A_EPS = BAWD_A_EPS;
constexpr double BAWR_MIN_SPAN = BAWD_MIN_SPAN;
constexpr double BAWR_MIN_LOG_GAP = BAWD_MIN_LOG_GAP;
constexpr double BAWR_LOG_BRACKET_MIN = BAWD_LOG_BRACKET_MIN;
constexpr double BAWR_DENOM_FLOOR = BAWD_DENOM_FLOOR;

// The launch selector is BAwD's, as for BAwF: one definition of
// "0 = truncated normal, 1 = lognormal, 2 = continuous split-lognormal"
// across the ballistic family.
constexpr int BAWR_LAUNCH_NORMAL = BAWD_LAUNCH_NORMAL;
constexpr int BAWR_LAUNCH_LOGNORMAL = BAWD_LAUNCH_LOGNORMAL;
constexpr int BAWR_LAUNCH_SPLITLOGNORMAL = BAWD_LAUNCH_SPLITLOGNORMAL;
// --------------------------------------------------------------------------
// Geometry (launch-distribution free)
// --------------------------------------------------------------------------

struct BawrGeom {
  bool ok = false;
  bool kappa_zero = false;    // exact LBA limit: no decay, no endpoint
  double b = 0.0, A = 0.0, kappa = 0.0, pw = 1.0;
  double log_kappa = R_NegInf;
  double T_max = R_PosInf;    // saturation time at z = 0: the hard endpoint
  double T_sat_A = R_PosInf;  // saturation time at z = A (<= T_max)
  double V_c0 = R_PosInf;     // kappa T_max^p: the omission boundary at z = 0
};

// u_sat(z) = [(p+1)(b-z) / (kappa p)]^(1/(p+1)), in logs so that the
// exponentiation is a single exp() rather than a pow() of a ratio of
// possibly very different magnitudes.
double bawr_sat_time(const BawrGeom& g, double b_minus_z);

// V_c(s) = kappa s^p, the launch strength exactly tangent at saturation time s.
double bawr_log_critical_launch(const BawrGeom& g, double s);
double bawr_critical_launch(const BawrGeom& g, double s);

BawrGeom bawr_geometry(double A, double b, double kappa, double pw);


// Per-time quantities.  `Z` splits the start-point range into the live part
// [0, Z] and the already-saturated (frozen) part (Z, A]; start points ABOVE Z
// are the frozen ones because u_sat(z) decreases in z.
struct BawrAtU {
  bool ok = false;
  bool saturated = false;   // u >= T_max: Z = 0, only frozen mass remains
  bool partial = false;     // Z < A: a frozen contribution exists
  double q = 0.0;           // dz -> dw Jacobian; V* is (b - z)/u + ... so q = u
  double Z = 0.0;
  double c = 0.0;           // kappa u^p: the affine offset in the density
  double w_hi = 0.0;        // V*(u, 0), required launch at the lowest start
  double w_lo = 0.0;        // V*(u, Z)
  double s_lo = 0.0;        // frozen-integral limits, in saturation-time
  double s_hi = 0.0;        // coordinates
};

BawrAtU bawr_at_u(const BawrGeom& g, double u);

// --------------------------------------------------------------------------
// Frozen contribution: the mass of start points that have already saturated.
//
//   Psi_fr = int_{z=Z}^{A} Gbar(V_c(z)) dz.
//
// At tangency z = b - kappa p s^(p+1)/(p+1) and w = kappa s^p, so
//
//   |dz/ds| = kappa p s^p          (saturation-time coordinates)
//   |dz/dw| = kappa^(-1/p) w^(1/p) (launch coordinates)
//
// The second form is closed for the lognormal launch at any p: it is the
// positive-power survivor integral int w^(1/p) Gbar dw, i.e.
// log_lognormal_power_stoploss() with -(m + 1) = 1/p.  Like BAwF's, it is a
// difference of monotone primitives and loses its digits when the
// critical-launch interval is short, so it falls back to the s-coordinate
// quadrature -- which needs no root solve per node and is what the normal
// launch always uses.
// --------------------------------------------------------------------------

double bawr_log_frozen_quad(const BawrGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift, double delta);

double bawr_log_frozen_normal(const BawrGeom& g, double s_lo,
                                     double s_hi, double v, double sv);

double bawr_log_frozen_logn(const BawrGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma, double delta = 0.0);

// --------------------------------------------------------------------------
// log CDF
//
//   A F(u) = q [C(w_lo) - C(w_hi)] + Psi_fr(u),
//
// with C the stop-loss price E[(V - w)_+] and q = u.  The live term is the
// same "affine in z" reduction BAwD, BAwF and the LBA use.
// --------------------------------------------------------------------------

double log_bawr_cdf_normal(double u, const BawrGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor);

double log_bawr_cdf_logn(double u, const BawrGeom& g, double mu,
                                double sigma, double delta = 0.0);

double bawr_log_frozen_surv_quad(const BawrGeom& g, double s_lo,
                                        double s_hi, double p1, double p2,
                                        bool logn, bool posdrift, double delta);

double bawr_log_frozen_surv_normal(const BawrGeom& g, double s_lo,
                                          double s_hi, double v, double sv,
                                          bool posdrift);

double bawr_log_frozen_surv_logn(const BawrGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma,
                                        double delta = 0.0);

double log_bawr_surv_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor);

double log_bawr_surv_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, double delta = 0.0);

// --------------------------------------------------------------------------
// log PDF
//
//   A f(u) = int_{w_lo}^{w_hi} (w - kappa u^p) g(w) dw,
//
// a partial expectation minus a probability, both closed form.  The integrand
// vanishes exactly at w_lo while Z = z*(u), so the live and frozen parts meet
// continuously and the endpoint decay is quadratic for A > 0 (linear for a
// point start, where there is no shrinking interval to supply the second
// factor) -- the same generic endpoint order as BAwD and BAwF.
// --------------------------------------------------------------------------

double log_bawr_pdf_normal(double u, const BawrGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor);

double log_bawr_pdf_logn(double u, const BawrGeom& g, double mu,
                                double sigma, double delta = 0.0);

// --------------------------------------------------------------------------
// Guarded natural-space CDF & PDF evaluators.  Same acceptance contract as
// BAwD and BAwF: a false return means "use the log path".
// --------------------------------------------------------------------------

bool bawr_natural_cdf_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &cdf);

bool bawr_natural_cdf_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, int accept_mode, double &cdf,
                                  double delta = 0.0);

bool bawr_natural_pdf_normal(double u, const BawrGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &pdf);

bool bawr_natural_pdf_logn(double u, const BawrGeom& g, double mu,
                                  double sigma, int accept_mode, double &pdf,
                                  double delta = 0.0);

bool ba_natural_cdf_bawr(double u, double A, double b, double p1,
                                double p2, double kappa, double pw, int launch,
                                bool posdrift, double denom_floor,
                                int accept_mode, double &cdf, double delta = 0.0);

bool ba_natural_pdf_bawr(double u, double A, double b, double p1,
                                double p2, double kappa, double pw, int launch,
                                bool posdrift, double denom_floor,
                                int accept_mode, double &pdf, double delta = 0.0);

// --------------------------------------------------------------------------
// Dispatch and output wrappers.  `p1`/`p2` are (v, sv) for the normal launch
// and (mu, sigma) for the lognormal one; they occupy the same kernel columns.
// The trailing `delta` selects the continuous split-lognormal launch
// (BAWR_LAUNCH_SPLITLOGNORMAL); delta = 0 reduces every path to the plain
// lognormal numbers exactly.
// --------------------------------------------------------------------------
double bawr_log_cdf(double u, double A, double b, double p1, double p2,
                           double kappa, double pw, int launch, bool posdrift,
                           double denom_floor = BAWR_DENOM_FLOOR,
                           double delta = 0.0);

double bawr_log_surv(double u, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            double denom_floor = BAWR_DENOM_FLOOR,
                            double delta = 0.0);

double bawr_log_pdf(double u, double A, double b, double p1, double p2,
                           double kappa, double pw, int launch, bool posdrift,
                           double denom_floor = BAWR_DENOM_FLOOR,
                           double delta = 0.0);

double bawr_cdf_norm(double t, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor = BAWR_DENOM_FLOOR,
                            double delta = 0.0);

double bawr_pdf_norm(double t, double A, double b, double p1, double p2,
                            double kappa, double pw, int launch, bool posdrift,
                            bool log_out,
                            double denom_floor = BAWR_DENOM_FLOOR,
                            double delta = 0.0);

// Natural-scale scalar evaluators for consumers that clamp to [0, 1] and
// tolerate tail saturation: truncation normalisers and GSL integrands.
double bawr_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double kappa, double pw,
                                      int launch, bool posdrift,
                                      double delta = 0.0);

double bawr_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double kappa, double pw,
                                      int launch, bool posdrift,
                                      double delta = 0.0);



// --------------------------------------------------------------------------
// Race-adapter entry points, shared with particle_ll.cpp's pointer dispatch.
// Bodies live in model_BAwR.cpp.
// --------------------------------------------------------------------------
double dbawr_scalar(double t, const double* par, void* ctx_);
double pbawr_scalar(double t, const double* par, void* ctx_);
void dbawr_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void pbawr_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void bawr_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);
#endif  // EMC2_MODEL_BAWR_H
