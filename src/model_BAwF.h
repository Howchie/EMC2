#ifndef EMC2_MODEL_BAWF_H
#define EMC2_MODEL_BAWF_H

// ---------------------------------------------------------------------------
// BAwF: ballistic accumulation with global fading of decision-relevant
// evidence.
//
//   X(u) = h_rho(u) [z + V u],   h_rho(u) = (1 + k u / rho)^(-rho),
//                                h_Inf(u) = exp(-k u),   z ~ U(0, A), b = B + A
//
// The latent trace z + V u accumulates ballistically; h_rho(u) is a global
// temporal attenuation of how strongly that trace remains available to the
// decision mechanism.  Equivalently dX/du = V h_rho(u) - k X (for rho = Inf),
// which reads as "drive decay plus leak" but has ONE primitive rate: the
// clearance term is k X, so it borrows its evidence scale from the state
// instead of carrying an independent one.  That is the whole point of the
// model -- BAwD's `ell` is an independent evidence scale and is what its
// (B, k, ell) ridge is made of.
//
// Writing x = k u and H_rho(x) = 1 / h_rho(x / k) = (1 + x / rho)^rho, the
// crossing law X(u) = b is
//
//   V*(u, z) = (b H_rho(k u) - z) / u,
//
// which falls and then rises, so weak launches intrinsically omit.  Tangency
// is at z = b G(x), G(x) = H(x) - x H'(x), and G decreases from G(0) = 1 to
// G(x_max) = 0 with
//
//   x_max = 1               (rho = Inf)          T_max = x_max / k
//   x_max = rho/(rho - 1)   (finite rho > 1)
//
// so the hard right endpoint T_max depends ONLY on k and the fixed shape
// index -- B cannot counterfeit it.  The critical launch at tangency is
// V_c(z) = k b H'(x_sat(z)), so V_c(0) = e k b for rho = Inf: given the time
// scale, the omission boundary identifies b.  That separation of observables
// is the reason this model exists.
//
// rho = 1 has no finite endpoint (V* decreases monotonically to k b) and is a
// qualitatively different member; BAwF() admits only rho > 1.
//
// Normal, lognormal, and Weibull launch strengths share this geometry. The selected
// launch distribution is supplied by ContextForRaceModels::bawd_launch, which
// BAwF shares with BAwD (the two never coexist in one adapter).
//
// This header is included by particle_ll.cpp through utils.h, after
// model_BAwD.h: it reuses that header's positive-integrand log quadrature
// (bawd_log_gl_split) and bawd_pk_log_tau, and wald_functions.h's lognormal
// stop-loss primitives.
// ---------------------------------------------------------------------------

#include <cmath>
#include "model_BAwD.h"

constexpr double BAWF_K_EPS = BAWD_K_EPS;
constexpr double BAWF_A_EPS = BAWD_A_EPS;
constexpr double BAWF_MIN_SPAN = BAWD_MIN_SPAN;
constexpr double BAWF_MIN_LOG_GAP = BAWD_MIN_LOG_GAP;
constexpr double BAWF_LOG_BRACKET_MIN = BAWD_LOG_BRACKET_MIN;
constexpr double BAWF_DENOM_FLOOR = BAWD_DENOM_FLOOR;

// The launch selector is BAwD's: identical meaning, and sharing it keeps one
// definition of "0 = truncated normal, 1 = lognormal, 3 = Weibull" for both models.
constexpr int BAWF_LAUNCH_NORMAL = BAWD_LAUNCH_NORMAL;
constexpr int BAWF_LAUNCH_LOGNORMAL = BAWD_LAUNCH_LOGNORMAL;
constexpr int BAWF_LAUNCH_SPLITLOGNORMAL = BAWD_LAUNCH_SPLITLOGNORMAL;
constexpr int BAWF_LAUNCH_WEIBULL = BAWD_LAUNCH_WEIBULL;
// --------------------------------------------------------------------------
// Geometry (launch-distribution free)
// --------------------------------------------------------------------------

struct BawfGeom {
  bool ok = false;
  bool k_zero = false;        // exact LBA limit: no fading, no endpoint
  bool rho_inf = true;
  double rho = R_PosInf;
  double b = 0.0, A = 0.0, k = 0.0;
  double y_0 = 1.0;           // x_max, the saturation phase at z = 0
  double y_A = 0.0;           // saturation phase at z = A (<= y_0)
  double T_max = R_PosInf;    // y_0 / k
  double T_sat_A = R_PosInf;  // y_A / k
  double V_c0 = R_PosInf;     // k b H'(y_0): the omission boundary at z = 0
};

// H(x) = (1 + x/rho)^rho, H'(x) = (1 + x/rho)^(rho-1),
// H''(x) = ((rho-1)/rho) (1 + x/rho)^(rho-2); all three are e^x at rho = Inf.
// Logs are used throughout because the finite-rho powers are otherwise a
// pow() per node of the frozen quadrature.
double bawf_log_H(const BawfGeom& g, double x);
double bawf_log_Hp(const BawfGeom& g, double x);
double bawf_log_Hpp(const BawfGeom& g, double x);

// G(x) = H - x H' = z*/b, written as a product so that it stays exact at the
// endpoint where H and x H' cancel: G = H'(x) (1 - x / x_max).
double bawf_zrel(const BawfGeom& g, double x);

// V_c(s) = k b H'(s), the launch strength that is exactly tangent at phase s.
double bawf_critical_launch(const BawfGeom& g, double s);

// Inverse of the above in log space: the phase whose critical launch is w,
// given log_ratio = log(w / (k b)).  Used to place the frozen quadrature's
// split point at the launch density's turnover.
double bawf_s_of_logratio(const BawfGeom& g, double log_ratio);

// Solve G(x) = c on [0, x_max] for c in [0, 1]; G is strictly decreasing there
// with G'(x) = -x H''(x), so a bracketed Newton is unconditionally safe.  Near
// c = 1 the root is O(sqrt(1 - c)) because G'(0) = 0, which is what the
// initial guess encodes.
double bawf_newton_x(const BawfGeom& g, double c);

BawfGeom bawf_geometry(double A, double b, double k, double rho);

// Per-time quantities.  `Z` splits the start-point range into the live part
// [0, Z] and the already-saturated (frozen) part (Z, A]; start points ABOVE
// Z are the frozen ones because x_sat(z) decreases in z.
struct BawfAtU {
  bool ok = false;
  bool saturated = false;   // u >= T_max: Z = 0, only frozen mass remains
  bool partial = false;     // Z < A: a frozen contribution exists
  double q = 0.0;           // dz -> dw Jacobian; V* is (b H - z)/u so q = u
  double Z = 0.0;
  double c = 0.0;           // k b H'(k u): the affine offset in the density
  double w_hi = 0.0;        // V*(u, 0), required launch at the lowest start
  double w_lo = 0.0;        // V*(u, Z)
  double s_lo = 0.0;        // frozen-integral limits, in phase coordinates
  double s_hi = 0.0;
};

BawfAtU bawf_at_u(const BawfGeom& g, double u);

// --------------------------------------------------------------------------
// Frozen contribution: the mass of start points that have already saturated.
//
// At tangency z = b G(s) and w = k b H'(s), so
//
//   |dz/dw| = s / k        (rho = Inf, where s = log(w / (k b)))
//   |dz/dw| = (rho/k) [(w/(k b))^(1/(rho-1)) - 1]     (finite rho)
//
// and Psi_fr = int_{w_a}^{w_b} |dz/dw| Gbar(w) dw.  For the lognormal launch
// both forms are closed: the first is the log-ratio survivor primitive
// log_lognormal_logratio_stoploss() with the reference k b, and the second is
// a positive-power survivor integral log_lognormal_power_stoploss() minus the
// plain stop-loss.  Both are differences of monotone primitives and lose
// their digits when the critical-launch interval concentrates at w = k b, so
// each falls back to the positive s-space quadrature below.
//
// In phase coordinates |dz/ds| = s b H''(s), which needs no root solve per
// node and is what the normal launch always uses.
// --------------------------------------------------------------------------

double bawf_log_frozen_quad(const BawfGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift);

double bawf_log_frozen_normal(const BawfGeom& g, double s_lo,
                                     double s_hi, double v, double sv);

double bawf_log_frozen_logn(const BawfGeom& g, double s_lo, double s_hi,
                                   double mu, double sigma,
                                   double delta = 0.0);
double bawf_log_frozen_weib(const BawfGeom& g, double s_lo, double s_hi,
                                   double shape, double mean);

// --------------------------------------------------------------------------
// log CDF
//
//   A F(u) = q [C(w_lo) - C(w_hi)] + Psi_fr(u),
//
// with C the stop-loss price E[(V - w)_+] and q = u.  The live term is the
// same "affine in z" reduction BAwD and the LBA use.
// --------------------------------------------------------------------------

double log_bawf_cdf_normal(double u, const BawfGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor);

double log_bawf_cdf_logn(double u, const BawfGeom& g, double mu,
                                double sigma, double delta = 0.0);

// Log survivor: the same start-point integral as the CDF, with the launch
// survivor replaced by the launch CDF.  This is deliberately independent of
// log1m_exp(log_cdf), which saturates as soon as the CDF rounds to one.
double bawf_log_frozen_surv_quad(const BawfGeom& g, double s_lo,
                                        double s_hi, double p1, double p2,
                                        bool logn, bool posdrift);

double bawf_log_frozen_surv_normal(const BawfGeom& g, double s_lo,
                                          double s_hi, double v, double sv,
                                          bool posdrift);

double bawf_log_frozen_surv_logn(const BawfGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma,
                                        double delta = 0.0);

double log_bawf_surv_normal(double u, const BawfGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor);

double log_bawf_surv_logn(double u, const BawfGeom& g, double mu,
                                  double sigma, double delta = 0.0);

// --------------------------------------------------------------------------
// log PDF
//
//   A f(u) = int_{w_lo}^{w_hi} (w - k b H'(k u)) g(w) dw,
//
// a partial expectation minus a probability, both closed form.  The integrand
// vanishes exactly at w_lo when Z = z*(u), so the live and frozen parts meet
// continuously and the endpoint decay is quadratic for A > 0 (linear for a
// point start, where there is no shrinking interval to supply the second
// factor).
// --------------------------------------------------------------------------

double log_bawf_pdf_normal(double u, const BawfGeom& g, double v,
                                  double sv, bool posdrift,
                                  double denom_floor);

double log_bawf_pdf_logn(double u, const BawfGeom& g, double mu,
                                double sigma, double delta = 0.0);

// --------------------------------------------------------------------------
// Guarded natural-space CDF & PDF evaluators for BAwF.  Same acceptance
// contract as BAwD: a false return means "use the log path".
// --------------------------------------------------------------------------

bool bawf_natural_cdf_normal(double u, const BawfGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &cdf);

bool bawf_natural_cdf_logn(double u, const BawfGeom& g, double mu,
                                  double sigma, int accept_mode, double &cdf,
                                  double delta = 0.0);

bool bawf_natural_pdf_normal(double u, const BawfGeom& g, double v,
                                    double sv, bool posdrift,
                                    double denom_floor, int accept_mode,
                                    double &pdf);

bool ba_natural_cdf_bawf(double u, double A, double b, double p1,
                                double p2, double k, int launch, bool posdrift,
                                double rho, double denom_floor,
                                int accept_mode, double &cdf,
                                double delta = 0.0);

bool ba_natural_pdf_bawf(double u, double A, double b, double p1,
                                double p2, double k, int launch, bool posdrift,
                                double rho, double denom_floor,
                                int accept_mode, double &pdf,
                                double delta = 0.0);

// --------------------------------------------------------------------------
// Dispatch and output wrappers.  `p1`/`p2` are (v, sv) for the normal launch
// and (mu, sigma) for the lognormal one; they occupy the same kernel columns.
// --------------------------------------------------------------------------
double bawf_log_cdf(double u, double A, double b, double p1, double p2,
                           double k, int launch, bool posdrift, double rho,
                           double denom_floor = BAWF_DENOM_FLOOR,
                           double delta = 0.0);

double bawf_log_surv(double u, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, double rho,
                            double denom_floor = BAWF_DENOM_FLOOR,
                            double delta = 0.0);

double bawf_log_pdf(double u, double A, double b, double p1, double p2,
                           double k, int launch, bool posdrift, double rho,
                           double denom_floor = BAWF_DENOM_FLOOR,
                           double delta = 0.0);

double bawf_cdf_norm(double t, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, bool log_out,
                            double rho,
                            double denom_floor = BAWF_DENOM_FLOOR,
                            double delta = 0.0);

double bawf_pdf_norm(double t, double A, double b, double p1, double p2,
                            double k, int launch, bool posdrift, bool log_out,
                            double rho,
                            double denom_floor = BAWF_DENOM_FLOOR,
                            double delta = 0.0);

// Natural-scale scalar evaluators for consumers that clamp to [0, 1] and
// tolerate tail saturation: truncation normalisers and GSL integrands.
double bawf_cdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, int launch,
                                      bool posdrift, double rho,
                                      double delta = 0.0);

double bawf_pdf_scalar_natural(double t, double A, double b, double p1,
                                      double p2, double k, int launch,
                                      bool posdrift, double rho,
                                      double delta = 0.0);



// --------------------------------------------------------------------------
// Race-adapter entry points, shared with particle_ll.cpp's pointer dispatch.
// Bodies live in model_BAwF.cpp.
// --------------------------------------------------------------------------
double dbawf_scalar(double t, const double* par, void* ctx_);
double pbawf_scalar(double t, const double* par, void* ctx_);
void dbawf_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void pbawf_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void bawf_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);
#endif  // EMC2_MODEL_BAWF_H
