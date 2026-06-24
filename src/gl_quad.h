#ifndef gl_quad_h
#define gl_quad_h

// ---------------------------------------------------------------------------
// Self-contained fixed-order Gauss-Legendre quadrature.
//
// The bundled GSL subset (see gsl_bundled.c) compiles in ONLY the adaptive
// routines (QAGS/QAGIU); gsl_integration_glfixed is declared in the header but
// its source (glfixed.c) is NOT compiled, so it would be a link error. To avoid
// touching the GSL bundle / Makevars we provide our own fixed GL rule here.
//
// Nodes/weights on [-1,1] are computed once per node-count via Newton iteration
// on the Legendre polynomial (the classic "gauleg" algorithm) and cached
// thread-locally so repeated likelihood calls pay the setup cost only once.
//
// Integrand signature matches gsl_function: double f(double x, void* params),
// so the existing captureless integrand lambdas can be reused verbatim.
// ---------------------------------------------------------------------------

#include <vector>
#include <map>
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// Finite-window heuristic for the GL/auto stop-success routes:
// [muS - K_SIGMA*sigS - K_TAU*tauS, muS + K_SIGMA*sigS + K_TAU*tauS].
// Single source of truth: used as the C++ defaults, passed explicitly by the
// live call sites in particle_ll.cpp, and mirrored by the R routes
// (stop_success_gl.R) — keep all of them in sync via these constants.
constexpr double SS_WINDOW_K_SIGMA = 8.0;
constexpr double SS_WINDOW_K_TAU   = 16.0;

struct GLRule {
  std::vector<double> x;   // nodes on [-1, 1]
  std::vector<double> w;   // weights
};

// Compute Gauss-Legendre nodes/weights on [-1,1] for n points (Newton/gauleg).
inline void gl_compute_rule(int n, GLRule& rule) {
  rule.x.assign(n, 0.0);
  rule.w.assign(n, 0.0);
  const double eps = 1e-15;
  const int m = (n + 1) / 2;          // by symmetry only need half the roots
  for (int i = 0; i < m; ++i) {
    double z = std::cos(M_PI * (i + 0.75) / (n + 0.5));   // initial guess
    double z1, pp;
    do {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 0; j < n; ++j) {     // Legendre recurrence to degree n
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);   // derivative
      z1 = z;
      z  = z1 - p1 / pp;                          // Newton step
    } while (std::fabs(z - z1) > eps);
    rule.x[i]         = -z;
    rule.x[n - 1 - i] =  z;
    rule.w[i]         = 2.0 / ((1.0 - z * z) * pp * pp);
    rule.w[n - 1 - i] = rule.w[i];
  }
}

// Thread-local cache keyed by node count.
inline const GLRule& gl_get_rule(int n) {
  static thread_local std::map<int, GLRule> cache;
  auto it = cache.find(n);
  if (it == cache.end()) {
    GLRule r;
    gl_compute_rule(n, r);
    it = cache.emplace(n, std::move(r)).first;
  }
  return it->second;
}

// Integrate f(x, params) over [a, b] with an n-point Gauss-Legendre rule.
inline double gl_integrate(double (*f)(double, void*), void* params,
                           double a, double b, int n) {
  if (!(b > a)) return 0.0;
  const GLRule& r = gl_get_rule(n);
  const double c1 = 0.5 * (b - a);   // scale
  const double c2 = 0.5 * (b + a);   // shift
  double s = 0.0;
  for (int i = 0; i < n; ++i) {
    s += r.w[i] * f(c1 * r.x[i] + c2, params);
  }
  return c1 * s;
}

// ---------------------------------------------------------------------------
// Process-global stop_method configuration for the LIVE stop-success path.
//
// Set from R via emc2_set_stop_method() (see model_SS_EXG.h). SSEXG()/SSRDEX()
// store stop_method/stop_n_nodes in the model list and calc_ll_manager() calls
// the setter once per likelihood call before entering C++. The live
// dispatchers (ss_texg_stop_success_lpdf_live / ss_rdex_stop_success_lpdf_live
// in the model headers) read this config per stop-success evaluation.
//
// Process-global (NOT thread_local): set once per likelihood call, read many
// times. thread_local would default-construct on every worker thread and
// silently ignore the choice; a plain static propagates to all threads in the
// process. Forked (mclapply) workers inherit it from the parent; PSOCK workers
// each run calc_ll_manager() themselves, so it is set in every worker process.
// ---------------------------------------------------------------------------
// NB: emc2_get_stop_method() (model_SS_EXG.h) indexes its method_names array
// by these values — keep the two in sync if methods are added or reordered.
enum StopMethod {
  STOP_METHOD_AUTO      = 0,   // DEFAULT (Andrew, 2026-06-11)
  STOP_METHOD_INTEGRATE = 1,   // original adaptive qags/qagiu route, untouched
  STOP_METHOD_GL        = 2,   // fixed Gauss-Legendre with n_nodes
  STOP_METHOD_ANALYTIC  = 3    // EXG n_go==1 closed form (GL fallback otherwise)
};

struct StopMethodConfig {
  int method  = STOP_METHOD_AUTO;
  int n_nodes = 64;
};

inline StopMethodConfig& stop_method_config() {
  static StopMethodConfig cfg;
  return cfg;
}

// "auto" GL node-density bump. A fixed n-node rule under-resolves a sharp
// (near-Gaussian, small-tauS) stop density when the integration window spans
// many stop-sigma widths. The Gaussian core width is ~sigS and does NOT shrink
// as tauS -> 0, so we hold a minimum node DENSITY per sigS across the window
// rather than triggering on a single width ratio. The old rule
// ((ub-lo)/sigS > 40 -> 128) was blind to the small-tauS regime: its window is
// dominated by muS/sigS + 8 and the only tauS-dependent term, 16*tauS/sigS,
// SHRANK toward the trigger exactly as the peak got sharper, so it never fired
// where it was most needed (e.g. muS/sigS ~ 6, tauS << sigS -> ratio ~ 20).
//
//   want  = GL_NODES_PER_SIG * (ub - lo) / sigS               [target nodes]
//   n_eff = clamp(roundup(want, GL_NODE_STEP), n_nodes, GL_MAX_NODES)
//
// n_nodes is a floor (never reduced). Quantising up to a multiple of
// GL_NODE_STEP keeps only a handful of distinct rules in the thread-local
// gauleg cache (setup cost is per node-count). Deterministic and cheap;
// mirrored in R by gl_auto_nodes_R(). Tune GL_NODES_PER_SIG for the
// accuracy/speed trade-off (6 ~= 14 nodes across the Gaussian-core FWHM).
constexpr double GL_NODES_PER_SIG = 6.0;
constexpr int    GL_NODE_STEP     = 32;
constexpr int    GL_MAX_NODES     = 256;

inline int gl_auto_nodes(int n_nodes, double lo, double ub, double sigS) {
  if (!(sigS > 0.0) || !(ub > lo)) return n_nodes;
  const double want = GL_NODES_PER_SIG * (ub - lo) / sigS;
  if (!(want > n_nodes)) return n_nodes;
  int n_eff = static_cast<int>(std::ceil(want / GL_NODE_STEP)) * GL_NODE_STEP;
  if (n_eff < n_nodes)     n_eff = n_nodes;
  if (n_eff > GL_MAX_NODES) n_eff = GL_MAX_NODES;
  return n_eff;
}

#endif
