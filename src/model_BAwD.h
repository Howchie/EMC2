#ifndef EMC2_MODEL_BAWD_H
#define EMC2_MODEL_BAWD_H

// ---------------------------------------------------------------------------
// BAwD: ballistic accumulator with a power-law base decay.
//
//   h_rho(u) = (1 + k u / rho)^(-rho), h_Inf(u) = exp(-k u)
//   dX/du = V h_rho(u) - ell h_rho(u)^gamma, z ~ U(0, A)
//   X(u) = z + V Q_rho(u) - ell R_rho(u),         b = B + A
//
// The fixed rho options are 1, 2, 4, and Inf; the fixed gamma options are
// 0, 1/2, 2/3, 3/4, and 1.  Neither is an estimated parameter.  For gamma < 1
// and k > 0 and ell > 0, finite-rho trajectories have a finite peak and a
// hard right endpoint T_max.  At gamma = 1, drive and clearance co-decay:
// X(u) = z + (V - ell) Q_rho(u), so there is no finite endpoint although weak
// launches still produce intrinsic omissions.
//
// Normal and lognormal launch strengths share this geometry.  The selected
// launch distribution is supplied by ContextForRaceModels::bawd_launch.
//
// This header is included by particle_ll.cpp directly and through utils.h.
// ---------------------------------------------------------------------------

#include <cmath>
#include "bawd_kernel.h"
#include "model_LBA.h"
#include "quad_templates.h"

constexpr double BAWD_K_EPS = 1e-10;
constexpr double BAWD_A_EPS = 1e-10;
constexpr double BAWD_ELL_EPS = 1e-12;
constexpr double BAWD_GAMMA_EPS = 1e-12;
// Below this standardized start-point span the midpoint-in-z limit is used
// instead.
constexpr double BAWD_MIN_SPAN = 1e-8;
// Relative width below which a positive endpoint difference is replaced by its
// midpoint limit (which is then accurate to the square of the width).
constexpr double BAWD_MIN_LOG_GAP = 1e-7;
// Trust a signed-log numerator while it retains this (log-scale) fraction of
// its largest term; shared with BAwL.
constexpr double BAWD_LOG_BRACKET_MIN = BAWL_LOG_BRACKET_MIN;
// BAwD reuses BAwL's normalizer floor so that the ell = 0, A = 0 member
// reproduces pleakyba() exactly rather than to within a floor difference.
constexpr double BAWD_DENOM_FLOOR = BAWL_DENOM_FLOOR;

// Launch-strength distribution selector (ContextForRaceModels::bawd_launch and
// the `launch` argument of the exported d/p functions -- keep them in sync).
constexpr int BAWD_LAUNCH_NORMAL = 0;
constexpr int BAWD_LAUNCH_LOGNORMAL = 1;

// --------------------------------------------------------------------------
// Geometry (launch-distribution free)
// --------------------------------------------------------------------------

// Row geometry: everything that depends on (A, b, k, ell) but not on u or on
// the launch distribution.  Two Newton solves per row, never per quadrature
// node, so the truncation and censoring paths that evaluate the CDF at several
struct BawdGeom {
  bool ok = false;
  bool k_zero = false;
  bool ell_zero = false;
  bool rho_inf = true;
  bool rho_one = false;
  double gamma = 0.0;
  bool gamma_zero = true;
  bool gamma_one = false;
  double omg = 1.0;
  double rho = R_PosInf;
  double m_shape = R_PosInf;
  double kq_inf = 1.0;
  double frozen_m = 0.0;      // alpha - 1 in the lognormal frozen term
  double frozen_alpha = 1.0;  // exponent alpha in the frozen Jacobian
  double b = 0.0, A = 0.0, k = 0.0, ell = 0.0;
  double y_A = 0.0;
  double y_0 = 0.0;
  double T_sat_A = R_PosInf;
  double T_max = R_PosInf;
  double V_c0 = R_PosInf;
};

// Per-time quantities.  `Z` splits the start-point range into the live part
// [0, Z] and the already-saturated (frozen) part (Z, A].
struct BawdAtU {
  bool ok = false;
  bool saturated = false;   // u >= T_max: Z = 0, only frozen mass remains
  bool partial = false;     // Z < A: a frozen contribution exists
  double q = 0.0;           // (1 - e^{-k u}) / k, = u at k = 0
  double E = 1.0;
  double log_E = 0.0;
  double E_g = 1.0;
  double log_E_g = 0.0;
  double E_rel = 1.0;
  double log_E_rel = 0.0;
  double Z = 0.0;
  double w_hi = 0.0;        // V*(u, 0), required launch at the lowest start
  double w_lo = 0.0;        // V*(u, Z)
  double s_lo = 0.0;        // frozen-integral exponent limits (y coordinates)
  double s_hi = 0.0;
};

// BAwD analytic core declarations.  Defaults are part of the public C++ API.
double bawd_em1my(double y);

double bawd_newton_y(double c);

double bawd_log_psi_prime(const BawdGeom& g, double x);
double bawd_log_h(const BawdGeom& g, double x);
double bawd_log_wrel(const BawdGeom& g, double x);
double bawd_kq(const BawdGeom& g, double x);
double bawd_kr(const BawdGeom& g, double x);
double bawd_psi(double s, const BawdGeom& g);
double bawd_newton_s(double c, const BawdGeom& g);
bool bawd_shape_flags(BawdGeom& g, double gamma, double rho);
BawdGeom bawd_geometry(double A, double b, double k, double ell,
                       double gamma, double rho = R_PosInf);
double bawd_critical_launch(const BawdGeom& g, double s);
double bawd_c_gamma(const BawdGeom& g, double u, double q);
BawdAtU bawd_at_u(const BawdGeom& g, double u);
double bawd_log_frozen_normal(const BawdGeom& g, double s_lo, double s_hi,
                              double v, double sv);
double bawd_log_frozen_logn_quad(const BawdGeom& g, double s_lo, double s_hi,
                                 double mu, double sigma);
double bawd_log_frozen_logn(const BawdGeom& g, double s_lo, double s_hi,
                            double mu, double sigma);
double log_bawd_cdf_normal(double u, const BawdGeom& g, double v, double sv,
                           bool posdrift, double denom_floor);
double log_bawd_cdf_logn(double u, const BawdGeom& g, double mu, double sigma);
double bawd_log_frozen_surv_normal(const BawdGeom& g, double s_lo, double s_hi,
                                   double v, double sv, bool posdrift,
                                   double denom_floor);
double bawd_log_frozen_surv_logn(const BawdGeom& g, double s_lo, double s_hi,
                                 double mu, double sigma);
double log_bawd_surv_normal(double u, const BawdGeom& g, double v, double sv,
                            bool posdrift, double denom_floor);
double log_bawd_surv_logn(double u, const BawdGeom& g, double mu, double sigma);
double log_bawd_pdf_normal(double u, const BawdGeom& g, double v, double sv,
                           bool posdrift, double denom_floor);
double log_bawd_pdf_logn(double u, const BawdGeom& g, double mu, double sigma);
bool bawd_natural_cdf_normal(double u, const BawdGeom& g, double v, double sv,
                             bool posdrift, double denom_floor, int accept_mode,
                             double &cdf);
bool bawd_natural_cdf_logn(double u, const BawdGeom& g, double mu, double sigma,
                           int accept_mode, double &cdf);
bool bawd_natural_pdf_normal(double u, const BawdGeom& g, double v, double sv,
                             bool posdrift, double denom_floor, int accept_mode,
                             double &pdf);
bool bawd_natural_pdf_logn(double u, const BawdGeom& g, double mu, double sigma,
                           int accept_mode, double &pdf);
bool ba_natural_cdf_bawd(double u, double A, double b, double p1, double p2,
                         double k, double ell, int launch, bool posdrift,
                         double gamma, double rho, double denom_floor,
                         int accept_mode, double &cdf);
bool ba_natural_pdf_bawd(double u, double A, double b, double p1, double p2,
                         double k, double ell, int launch, bool posdrift,
                         double gamma, double rho, double denom_floor,
                         int accept_mode, double &pdf);
double bawd_log_cdf(double u, double A, double b, double p1, double p2,
                    double k, double ell, int launch, bool posdrift,
                    double gamma, double rho,
                    double denom_floor = BAWD_DENOM_FLOOR);
double bawd_log_surv(double u, double A, double b, double p1, double p2,
                     double k, double ell, int launch, bool posdrift,
                     double gamma, double rho,
                     double denom_floor = BAWD_DENOM_FLOOR);
double bawd_log_pdf(double u, double A, double b, double p1, double p2,
                    double k, double ell, int launch, bool posdrift,
                    double gamma, double rho,
                    double denom_floor = BAWD_DENOM_FLOOR);
double bawd_cdf_norm(double t, double A, double b, double p1, double p2,
                     double k, double ell, int launch, bool posdrift,
                     bool log_out, double gamma, double rho,
                     double denom_floor = BAWD_DENOM_FLOOR);
double bawd_pdf_norm(double t, double A, double b, double p1, double p2,
                     double k, double ell, int launch, bool posdrift,
                     bool log_out, double gamma, double rho,
                     double denom_floor = BAWD_DENOM_FLOOR);
double bawd_cdf_scalar_natural(double t, double A, double b, double p1,
                               double p2, double k, double ell, int launch,
                               bool posdrift, double gamma, double rho);
double bawd_pdf_scalar_natural(double t, double A, double b, double p1,
                               double p2, double k, double ell, int launch,
                               bool posdrift, double gamma, double rho);

// BAwD race-model adapter entry points.
// Definitions live in model_BAwD.cpp so utils.h remains a
// declaration-only integration point for these adapters.
double dbawd_scalar(double t, const double* par, void* ctx_);
double pbawd_scalar(double t, const double* par, void* ctx_);
void dbawd_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void pbawd_raw(const double* rt, const double* const* cols, int n_rows,
               const int* mask, const int* isok,
               double* out, double min_ll, void* ctx_);
void bawd_logS_at_t(double t, const double* const* cols,
                    int n_rows_total, int n_lR, int n_par,
                    const int* trunc_mask, int n_unique_trials,
                    const int* isok_all, void* ctx_, double* logS_out);
// BAwDp race-model adapter entry points.
// Definitions live in model_BAwD.cpp so utils.h remains a
// declaration-only integration point for these adapters.
double dbawdp_scalar(double t, const double* par, void* ctx_);
double pbawdp_scalar(double t, const double* par, void* ctx_);
void dbawdp_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok,
                double* out, double min_ll, void* ctx_);
void pbawdp_raw(const double* rt, const double* const* cols, int n_rows,
                const int* mask, const int* isok,
                double* out, double min_ll, void* ctx_);
void bawdp_logS_at_t(double t, const double* const* cols,
                     int n_rows_total, int n_lR, int n_par,
                     const int* trunc_mask, int n_unique_trials,
                     const int* isok_all, void* ctx_, double* logS_out);


#endif  // EMC2_MODEL_BAWD_H
