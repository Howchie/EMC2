#ifndef EMC2_MODEL_BTAWL_H
#define EMC2_MODEL_BTAWL_H

// BTAwL: ballistic transient/sustained local race with leak.
//
//   transient: dX_T/du = V_T (u/tau_t) exp(-u/tau_t) - k X_T
//   sustained: dX_S/du = V_S (1 - exp(-u/tau_s)) - k X_S
//   X_T(0) = z_T, X_S(0) = z_S, z_T and z_S independently ~ U(0,A)
//
// The transient and sustained members are raced locally. The transient
// response is affine in the launch point and the required
// launch V*(u,z) is therefore affine in z on the live part of the CDF.  The
// live start-point average is evaluated in closed form (normal and lognormal
// launches).  The already-saturated start region is a one-dimensional
// quadrature in tangency time; there is no grid or PDE solve.
//
// This header is included by model_LBA.h, after the shared normal/lognormal
// primitives and the GSL helpers have been declared.

#include <cmath>
#include <functional>
#include <algorithm>
#include <limits>
#include "model_BAwD.h"

// Keep BTAwL's launch selector and numerical thresholds identical to BAwL.
constexpr int BTAWL_LAUNCH_NORMAL = BAWL_LAUNCH_NORMAL;
constexpr int BTAWL_LAUNCH_LOGNORMAL = BAWL_LAUNCH_LOGNORMAL;
constexpr double BTAWL_K_EPS = BAWL_K_EPS;
constexpr double BTAWL_A_EPS = BAWL_A_EPS;
constexpr double BTAWL_DENOM_FLOOR = BAWL_DENOM_FLOOR;

struct BtawlGeom {
  bool ok = false;
  bool k_zero = false;
  double A = 0.0, b = 0.0, k = 0.0, tau = 0.0;
  double t_max = R_PosInf;
  double h_max = 0.0;
  // Tangency time for the upper start point z = A.  Frozen quadratures use
  // this as their lower endpoint; keeping it with the geometry avoids solving
  // the same monotone equation for every evaluation.
  double s_lo = R_PosInf;
};

double btawl_h(double t, double k, double tau);

double btawl_g(double t, double tau);

double btawl_hp(double t, double k, double tau);

double btawl_tmax(double k, double tau);

double btawl_tangent_time(double z, const BtawlGeom& g);

// Invert the endpoint chart Ttrans = T_max(k, tau).  The endpoint is strictly
// increasing in tau; the bracket (0, Ttrans) is guaranteed by the transient
// geometry, so a safeguarded bisection is sufficient and has no branch
// ambiguity in the stiff k*tau regime.
double btawl_tau_from_ttrans(double k, double Ttrans);

BtawlGeom btawl_geometry(double A, double b, double k, double tau);

// Endpoint-chart geometry: the R-side transform has already solved
// Ttrans = T_max(k, tau).  Reuse that solved value instead of running the
// forward endpoint root finder on every density/CDF evaluation.
BtawlGeom btawl_geometry_from_ttrans(double A, double b, double k,
                                            double Ttrans,
                                            double tau_hint = R_NaN);

double btawl_vstar(double t, double z, const BtawlGeom& g);

// Derivative of V*(t,z).  It is negative on the live/rising limb.
double btawl_vstar_prime(double t, double z, const BtawlGeom& g);

// Start point at which the tangency time equals t.  Gamma starts at one,
// rises to its maximum at u = tau, and then decreases to zero at t_max.  The
// clamp below is therefore load-bearing: Gamma is not globally monotone, and
// a Newton solve seeded from that assumption can select the wrong branch.
double btawl_z_t(double t, const BtawlGeom& g);

double btawl_tangent_z_prime(double t, const BtawlGeom& g);

double btawl_normal_denom(double v, double sv, bool posdrift);

double btawl_surv(double w, double p1, double p2, int launch,
                         bool posdrift);

double btawl_pdf_v(double w, double p1, double p2, int launch,
                          bool posdrift);

double btawl_J(double x);

double btawl_live_cdf(double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift);

double btawl_live_pdf(double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift);

double btawl_log_frozen(bool survivor, double t, const BtawlGeom& g,
                               double p1, double p2,
                               int launch, bool posdrift);

double btawl_frozen_cdf(double t, double zlo, double zhi,
                               const BtawlGeom& g, double p1, double p2,
                               int launch, bool posdrift);

// Launch-density turnover on the descending branch of g(s) = (s/tau)e^{-s/tau}.
// The frozen quadrature is concentrated near this point in the stiff endpoint
// window, so an arithmetic midpoint is a particularly poor split.
double btawl_frozen_split(const BtawlGeom& g, double s_lo, double s_hi,
                                 double target);

double btawl_log_live(bool survivor, double t, double zlo, double zhi,
                             const BtawlGeom& g, double p1, double p2,
                             int launch, bool posdrift);

double btawl_log_eval(bool survivor, double t, const BtawlGeom& g,
                             double p1, double p2, int launch, bool posdrift);

double log_btawl_cdf_normal(double t, const BtawlGeom& g, double v,
                                   double sv, bool posdrift);
double log_btawl_cdf_logn(double t, const BtawlGeom& g, double mu,
                                 double sigma);
double log_btawl_surv_normal(double t, const BtawlGeom& g, double v,
                                    double sv, bool posdrift);
double log_btawl_surv_logn(double t, const BtawlGeom& g, double mu,
                                  double sigma);

double btawl_cdf(double t, double A, double b, double p1, double p2,
                        double k, double tau, int launch, bool posdrift);

double btawl_pdf(double t, double A, double b, double p1, double p2,
                        double k, double tau, int launch, bool posdrift);

double btawl_cdf_from_geom(double t, const BtawlGeom& g, double p1,
                                  double p2, int launch, bool posdrift);

bool btawl_natural_cdf_from_geom(double t, const BtawlGeom& g,
                                        double p1, double p2, int launch,
                                        bool posdrift, double& cdf);

double btawl_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                                  double p2, int launch, bool posdrift);

double btawl_log_launch_pdf(double w, double p1, double p2, int launch,
                                   bool posdrift);

bool btawl_natural_pdf_accepted(double t, const BtawlGeom& g,
                                        double p1, double p2, int launch,
                                        bool posdrift, double p_nat);

double btawl_log_pdf_from_geom(double t, const BtawlGeom& g, double p1,
                                      double p2, int launch, bool posdrift);

double btawl_log_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift);
double btawl_log_surv(double t, double A, double b, double p1, double p2,
                             double k, double tau, int launch, bool posdrift);
double btawl_log_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift);
double btawl_cdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift);
double btawl_pdf_log(double t, double A, double b, double p1, double p2,
                            double k, double tau, int launch, bool posdrift);

double btawl_cdf_chart(double t, double A, double b, double p1,
                              double p2, double k, double clear, int launch,
                              bool posdrift, bool endpoint_chart,
                              double tau_hint = R_NaN);
double btawl_pdf_chart(double t, double A, double b, double p1,
                              double p2, double k, double clear, int launch,
                              bool posdrift, bool endpoint_chart,
                              double tau_hint = R_NaN);
double btawl_log_surv_chart(double t, double A, double b, double p1,
                                   double p2, double k, double clear,
                                   int launch, bool posdrift,
                                   bool endpoint_chart,
                                   double tau_hint = R_NaN);
double btawl_log_pdf_chart(double t, double A, double b, double p1,
                                  double p2, double k, double clear,
                                  int launch, bool posdrift,
                                  bool endpoint_chart,
                                  double tau_hint = R_NaN);

// ---------------------------------------------------------------------------
// Sustained + transient local-race extension.
// ---------------------------------------------------------------------------

double btawl_hs(double t, double k, double tau);

double btawl_hs_p(double t, double k, double tau);

double btawl_sustained_cdf(double t, double A, double b, double p1, double p2,
                                  double k, double tau_s, int launch, bool posdrift);

double btawl_sustained_pdf(double t, double A, double b, double p1, double p2,
                                  double k, double tau_s, int launch, bool posdrift);

double btawl_local_race_cdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift);

double btawl_local_race_pdf(double t, double A, double b, double p1, double p2,
                            double k, double tau_s, double tau_t, double pi,
                            int launch, bool posdrift);

double btawl_local_race_cdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift);

double btawl_local_race_pdf_log(double t, double A, double b, double p1, double p2,
                                double k, double tau_s, double tau_t, double pi,
                                int launch, bool posdrift);

double btawl_sustained_log_surv(double t, double A, double b, double p1, double p2,
                                       double k, double tau_s, int launch, bool posdrift);

double btawl_local_race_log_surv(double t, double A, double b, double p1,
                                 double p2, double k, double tau_s,
                                 double tau_t, double pi, int launch,
                                 bool posdrift);

// BTAwL adapter entry points used by particle_ll.cpp.
double dbtawl_transient_scalar(double t, const double* par, void* ctx_);
double pbtawl_transient_scalar(double t, const double* par, void* ctx_);
void dbtawl_transient_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_);
void pbtawl_transient_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_);
void btawl_transient_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int n_par,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_,
                               double* logS_out);

double dbtawl_local_race_scalar(double t, const double* par, void* ctx_);
double pbtawl_local_race_scalar(double t, const double* par, void* ctx_);
void dbtawl_local_race_raw(const double* rt, const double* const* cols,
                           int n_rows, const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx_);
void pbtawl_local_race_raw(const double* rt, const double* const* cols,
                           int n_rows, const int* mask, const int* isok,
                           double* out, double min_ll, void* ctx_);
void btawl_local_race_logS_at_t(double t, const double* const* cols,
                                int n_rows_total, int n_lR, int n_par,
                                const int* trunc_mask, int n_unique_trials,
                                const int* isok_all, void* ctx_, double* out);

double dbtawl_sustained_scalar(double t, const double* par, void* ctx_);
double pbtawl_sustained_scalar(double t, const double* par, void* ctx_);
void dbtawl_sustained_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_);
void pbtawl_sustained_raw(const double* rt, const double* const* cols,
                          int n_rows, const int* mask, const int* isok,
                          double* out, double min_ll, void* ctx_);
void btawl_sustained_logS_at_t(double t, const double* const* cols,
                               int n_rows_total, int n_lR, int n_par,
                               const int* trunc_mask, int n_unique_trials,
                               const int* isok_all, void* ctx_, double* out);

#endif
