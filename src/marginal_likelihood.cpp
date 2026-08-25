#include "marginal_likelihood.h"

#include "gl_quad.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

// Per-node grid of a t0 marginalization: the quadrature node values (sampled
// log-t0 scale) and the unnormalized per-particle log-terms
//   ell_ik = log quadrature_weight_ik + log p(x_ik|eta) + log L_i(x_ik)
// that BOTH the marginal-ll reducer and the storage-time weight accessor share,
// so the node grid and the reconstruction weights are guaranteed identical.
// Nodes are PER PARTICLE (np x K): each proposal gets a rule centred on its own
// conditional t0 posterior, so one particle's rule can never be dragged off its
// peak by the rest of a scattered batch.
struct MarginalGrid {
  Rcpp::NumericMatrix nodes;      // np x K: per-particle log-t0 node values
  Rcpp::NumericMatrix log_terms;  // np x K: ell_ik (softmax over k => p(t0=x_ik|y,theta,eta))
  Rcpp::NumericVector mode;       // np: fitted conditional mode (warm-start carrier)
  Rcpp::NumericVector sd;         // np: fitted Laplace scale  (warm-start carrier)
  double warm_used = 0.0;         // fraction of rows retaining the warm hint
  int pred_used = 0;              // 1 when the subset-regression hint replaced it
  int repaired = 0;               // particles whose rule missed and was refitted
};

// A composite GL rule lets the t0 marginal retain support over its complete
// feasible interval while concentrating most nodes around the response-informed
// posterior bump.  Unlike simply truncating to a centred window, the two tail
// panels retain the integral's mass and therefore keep the approximation
// coherent even for a broad or weakly identified t0 posterior.
struct MarginalRule {
  std::vector<double> x;
  std::vector<double> w;
};

static void marginal_append_gl_panel(MarginalRule& out, double lo, double hi,
                                     int n) {
  if (n <= 0 || !(hi > lo)) return;
  const GLRule& rule = gl_get_rule(n);
  const double half = 0.5 * (hi - lo);
  const double mid = 0.5 * (hi + lo);
  for (int k = 0; k < n; ++k) {
    out.x.push_back(mid + half * rule.x[static_cast<size_t>(k)]);
    out.w.push_back(half * rule.w[static_cast<size_t>(k)]);
  }
}

// Evaluate all particles at their node columns in one kernel pass.
static NumericMatrix marginal_eval_nodes(
    const NumericMatrix& particle_matrix, const NumericMatrix& X, int t0col,
    DataFrame data, NumericVector constants, List designs, String type,
    List bounds, List transforms, List pretransforms, CharacterVector p_types,
    double min_ll, Rcpp::Nullable<Rcpp::List> trend) {
  const int np = particle_matrix.nrow();
  const int n_par = particle_matrix.ncol();
  const int m = X.ncol();
  NumericMatrix out(np, m);
  if (np == 0 || m == 0) return out;
  CharacterVector p_names = colnames(particle_matrix);
  // Cap the rows per call so a wide grid on an uncompressed data set cannot
  // blow the expanded parameter matrix up.
  const int max_rows = 4096;
  int cols_per_call = m;
  if (np * m > max_rows) cols_per_call = std::max(1, max_rows / np);
  NumericMatrix pm(np * std::min(cols_per_call, m), n_par);
  colnames(pm) = p_names;
  for (int k0 = 0; k0 < m; k0 += cols_per_call) {
    const int mk = std::min(cols_per_call, m - k0);
    if (pm.nrow() != np * mk) {
      pm = NumericMatrix(np * mk, n_par);
      colnames(pm) = p_names;
    }
    for (int k = 0; k < mk; ++k) {
      for (int i = 0; i < np; ++i) {
        const int r = k * np + i;
        for (int j = 0; j < n_par; ++j) pm(r, j) = particle_matrix(i, j);
        pm(r, t0col) = X(i, k0 + k);
      }
    }
    NumericVector ll = calc_ll_oo(pm, data, constants, designs, type, bounds,
                                  transforms, pretransforms, p_types, min_ll,
                                  trend, R_NilValue);
    for (int k = 0; k < mk; ++k)
      for (int i = 0; i < np; ++i) out(i, k0 + k) = ll[k * np + i];
  }
  return out;
}

// Mode and scale of the parabola through three points of the log integrand
// g(x) = log L(x) + log p(x|eta).  Returns false unless the fit is a proper
// interior maximum, in which case the caller keeps the full feasible interval
// rather than a Laplace window.
static bool marginal_parabola(double x1, double g1, double x2, double g2,
                              double x3, double g3, double* mode, double* sd) {
  if (!(R_FINITE(g1) && R_FINITE(g2) && R_FINITE(g3))) return false;
  if (!(x1 < x2 && x2 < x3)) return false;
  const double d12 = (g2 - g1) / (x2 - x1);
  const double d23 = (g3 - g2) / (x3 - x2);
  const double curv = 2.0 * (d23 - d12) / (x3 - x1);   // ~ g''(mode)
  if (!R_FINITE(curv) || curv >= 0.0) return false;
  const double num = (x2 - x1) * (x2 - x1) * (g2 - g3) -
                     (x2 - x3) * (x2 - x3) * (g2 - g1);
  const double den = (x2 - x1) * (g2 - g3) - (x2 - x3) * (g2 - g1);
  if (!R_FINITE(den) || den == 0.0) return false;
  const double xm = x2 - 0.5 * num / den;
  if (!R_FINITE(xm)) return false;
  const double s = 1.0 / std::sqrt(-curv);
  if (!R_FINITE(s) || !(s > 0.0)) return false;
  *mode = xm;
  *sd = s;
  return true;
}

// Ridge-stabilised least squares of y on [1 X] (X is n x p, column-major), by
// Cholesky on the normal equations.  p is a handful of sampled parameters and n
// a few dozen probe particles, so this is negligible next to one kernel pass.
// Returns false when the system is not positive definite (rank-deficient design,
// e.g. a component whose columns are all constant across the batch).
static bool marginal_lsfit(const std::vector<double>& X, const std::vector<double>& y,
                           int n, int p, std::vector<double>& beta) {
  const int q = p + 1;
  if (n < q + 2) return false;
  std::vector<double> A(static_cast<size_t>(q) * q, 0.0), b(static_cast<size_t>(q), 0.0);
  auto xcol = [&](int j, int i) { return j == 0 ? 1.0 : X[static_cast<size_t>(j - 1) * n + i]; };
  for (int j = 0; j < q; ++j) {
    for (int k = j; k < q; ++k) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i) sum += xcol(j, i) * xcol(k, i);
      A[static_cast<size_t>(j) * q + k] = A[static_cast<size_t>(k) * q + j] = sum;
    }
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += xcol(j, i) * y[static_cast<size_t>(i)];
    b[static_cast<size_t>(j)] = sum;
  }
  double tr = 0.0;
  for (int j = 0; j < q; ++j) tr += A[static_cast<size_t>(j) * q + j];
  const double ridge = 1e-8 * (tr / q) + 1e-12;
  for (int j = 0; j < q; ++j) A[static_cast<size_t>(j) * q + j] += ridge;
  // Cholesky, in place
  for (int j = 0; j < q; ++j) {
    double d = A[static_cast<size_t>(j) * q + j];
    for (int k = 0; k < j; ++k) d -= A[static_cast<size_t>(j) * q + k] * A[static_cast<size_t>(j) * q + k];
    if (!(d > 0.0) || !R_FINITE(d)) return false;
    d = std::sqrt(d);
    A[static_cast<size_t>(j) * q + j] = d;
    for (int i = j + 1; i < q; ++i) {
      double v = A[static_cast<size_t>(i) * q + j];
      for (int k = 0; k < j; ++k) v -= A[static_cast<size_t>(i) * q + k] * A[static_cast<size_t>(j) * q + k];
      A[static_cast<size_t>(i) * q + j] = v / d;
    }
  }
  beta.assign(static_cast<size_t>(q), 0.0);
  for (int i = 0; i < q; ++i) {            // forward
    double v = b[static_cast<size_t>(i)];
    for (int k = 0; k < i; ++k) v -= A[static_cast<size_t>(i) * q + k] * beta[static_cast<size_t>(k)];
    beta[static_cast<size_t>(i)] = v / A[static_cast<size_t>(i) * q + i];
  }
  for (int i = q - 1; i >= 0; --i) {       // back
    double v = beta[static_cast<size_t>(i)];
    for (int k = i + 1; k < q; ++k) v -= A[static_cast<size_t>(k) * q + i] * beta[static_cast<size_t>(k)];
    beta[static_cast<size_t>(i)] = v / A[static_cast<size_t>(i) * q + i];
  }
  for (int j = 0; j < q; ++j) if (!R_FINITE(beta[static_cast<size_t>(j)])) return false;
  return true;
}

// Coherent marginalization core (option 3): integrate ONE shared t0 out of the
// COMPLETE subject likelihood by a per-particle Laplace-centred composite
// Gauss-Legendre rule on the log-t0 (sampled) axis.  Given a warm start
// (`warm_mode`/`warm_sd`: the values fitted for this subject's accepted
// particle last iteration) a single probe round re-fits each particle at the
// right scale; otherwise a coarse pilot scan
// brackets each particle's own conditional mode, two parabolic refinements
// resolve its curvature (the conditional posterior narrows like 1/sqrt(n_trials)
// and can be orders of magnitude thinner than any prior-scaled window), and the
// core panel is then flanked by one tail node each side so the whole feasible
// interval stays integrated.  Reuses the registered kernel verbatim -- nodes
// overwrite the sampled t0 column and re-enter calc_ll_oo with marginalise off
// -- so it is model-agnostic (every race model, no per-model code).  The
// interval is clipped to where the likelihood is live: below t0's natural-scale
// lower bound c_do_bound floors the whole subject to n_trials*min_ll; above
// min(rt) the fastest response is infeasible.  Correctness reference is
// WorkingTests/marginal_t0_lib.R::marginal_ll_t0.
static MarginalGrid calc_ll_oo_marginal_core(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
    Rcpp::List marginalise) {
  const int np = particle_matrix.nrow();
  const int n_par = particle_matrix.ncol();

  // --- marginalization spec: p(t0|eta) = lognormal, log t0 ~ N(mu, sigma) ---
  const std::string mparam = Rcpp::as<std::string>(marginalise["param"]);
  const double mu    = Rcpp::as<double>(marginalise["mu"]);
  const double sigma = Rcpp::as<double>(marginalise["sigma"]);
  int    n_nodes = marginalise.containsElementNamed("n_nodes")
                     ? Rcpp::as<int>(marginalise["n_nodes"]) : 12;
  int    n_scan = marginalise.containsElementNamed("n_scan")
                     ? Rcpp::as<int>(marginalise["n_scan"]) : 7;
  int    n_refine = marginalise.containsElementNamed("n_refine")
                     ? Rcpp::as<int>(marginalise["n_refine"]) : 2;
  const bool adaptive = marginalise.containsElementNamed("adaptive")
                     ? Rcpp::as<bool>(marginalise["adaptive"]) : true;
  const double span = marginalise.containsElementNamed("span")
                     ? Rcpp::as<double>(marginalise["span"]) : 6.0;
  // Half-width of the core panel in Laplace SDs: 5 covers a Gaussian core to
  // ~1e-6 relative, and the flanking tail nodes carry the rest.
  const double window_sd = marginalise.containsElementNamed("window_sd")
                     ? Rcpp::as<double>(marginalise["window_sd"]) : 5.0;
  const double eps  = marginalise.containsElementNamed("eps")
                     ? Rcpp::as<double>(marginalise["eps"]) : 1e-3;
  // Warm start: the mode/scale fitted for this subject's accepted particle on
  // the previous iteration.  The conditional t0 posterior moves very little
  // between MCMC iterations, so a valid hint replaces the whole pilot scan (and
  // one refinement) with a single probe round at the right scale.
  const double warm_mode = marginalise.containsElementNamed("warm_mode")
                     ? Rcpp::as<double>(marginalise["warm_mode"]) : NA_REAL;
  const double warm_sd = marginalise.containsElementNamed("warm_sd")
                     ? Rcpp::as<double>(marginalise["warm_sd"]) : NA_REAL;
  // Probe at a few hint-SDs so the bracket still spans the peak when the batch
  // has drifted; the fit is re-scaled from those points anyway.
  const double warm_inflate = marginalise.containsElementNamed("warm_inflate")
                     ? Rcpp::as<double>(marginalise["warm_inflate"]) : 4.0;
  int n_refine_warm = marginalise.containsElementNamed("n_refine_warm")
                     ? Rcpp::as<int>(marginalise["n_refine_warm"]) : 2;
  // Within-iteration hint: pilot a 1/pred_frac subset, regress its modes on the
  // other sampled parameters, and bracket the rest of the batch from that.
  const bool predict_mode = marginalise.containsElementNamed("predict_mode")
                     ? Rcpp::as<bool>(marginalise["predict_mode"]) : true;
  int pred_frac = marginalise.containsElementNamed("pred_frac")
                     ? Rcpp::as<int>(marginalise["pred_frac"]) : 8;
  int n_refine_pred = marginalise.containsElementNamed("n_refine_pred")
                     ? Rcpp::as<int>(marginalise["n_refine_pred"]) : 2;
  if (n_nodes < 2) n_nodes = 2;
  if (n_scan < 3) n_scan = 3;
  if (n_refine < 0) n_refine = 0;
  if (n_refine_warm < 1) n_refine_warm = 1;
  if (n_refine_pred < 1) n_refine_pred = 1;
  if (pred_frac < 2) pred_frac = 2;

  // --- locate the sampled t0 coordinate in the proposal matrix ---
  CharacterVector p_names = colnames(particle_matrix);
  int t0col = -1;
  for (int j = 0; j < p_names.size(); ++j)
    if (Rcpp::as<std::string>(p_names[j]) == mparam) { t0col = j; break; }
  if (t0col < 0)
    Rcpp::stop("calc_ll_oo marginalise: '%s' is not a sampled column.", mparam);

  // --- lower clip: t0's natural-scale bound (minmax), on the log axis. ---
  double log_lo = R_NegInf;
  {
    NumericMatrix minmax = bounds["minmax"];
    CharacterVector mm_names = colnames(minmax);
    for (int j = 0; j < mm_names.size(); ++j)
      if (Rcpp::as<std::string>(mm_names[j]) == mparam) {
        const double lo = minmax(0, j);
        if (R_FINITE(lo) && lo > 0.0) log_lo = std::log(lo);
        break;
      }
  }
  // --- upper clip: feasibility ceiling log(min finite rt - eps). ---
  double log_ms = R_PosInf;
  {
    NumericVector rt = data["rt"];
    double mn = R_PosInf;
    for (int k = 0; k < rt.size(); ++k)
      if (R_FINITE(rt[k]) && rt[k] < mn) mn = rt[k];
    if (R_FINITE(mn)) {
      double top = mn - eps;
      if (top < std::numeric_limits<double>::min()) top = std::numeric_limits<double>::min();
      log_ms = std::log(top);
    }
  }

  const double lo = std::max(mu - span * sigma, log_lo);
  const double hi = std::min(mu + span * sigma, log_ms);

  MarginalGrid g;
  g.nodes = NumericMatrix(np, n_nodes);
  g.log_terms = NumericMatrix(np, n_nodes);
  g.mode = NumericVector(np, NA_REAL);
  g.sd = NumericVector(np, NA_REAL);
  if (!(hi > lo)) {                                  // collapsed interval => no mass
    std::fill(g.nodes.begin(), g.nodes.end(), NA_REAL);
    std::fill(g.log_terms.begin(), g.log_terms.end(), R_NegInf);
    return g;
  }
  if (np == 0) return g;

  // --- per-particle Laplace fit of g(x) = log L(x) + log p(x|eta) -----------
  const double def_mode = 0.5 * (lo + hi), def_sd = 0.25 * (hi - lo);
  std::vector<double> mode(static_cast<size_t>(np), def_mode);
  std::vector<double> sdev(static_cast<size_t>(np), def_sd);
  std::vector<char> fitted(static_cast<size_t>(np), 0);
  const double max_sd = 0.5 * (hi - lo);

  // One refinement round: probe g at (mode - h, mode, mode + h) with h the
  // current scale estimate and re-fit the parabola THERE.  A fit at the wrong
  // spacing badly mis-states the width of a razor-thin posterior (it shrinks
  // like 1/sqrt(n_trials)), so each round contracts the probe onto the peak.
  // Returns how many particles the round did not resolve CLEANLY: either the
  // parabola was not a proper maximum, or its vertex sits outside the middle of
  // the probe bracket, meaning the incoming centre/scale was an extrapolation
  // rather than a local fit.  The warm path uses that count to decide whether
  // its hint is still usable; the cold path ignores it (its own next round is
  // the correction).
  auto refine_on = [&](const NumericMatrix& PM, std::vector<double>& md,
                       std::vector<double>& sv, std::vector<char>& ft) -> int {
    const int nn = PM.nrow();
    NumericMatrix Xr(nn, 3);
    for (int i = 0; i < nn; ++i) {
      const double h = std::max(std::min(sv[static_cast<size_t>(i)], max_sd), 1e-8);
      double c = std::min(hi, std::max(lo, md[static_cast<size_t>(i)]));
      double a = c - h, b = c + h;
      if (a < lo) { a = lo; c = std::min(hi, a + h); b = std::min(hi, c + h); }
      if (b > hi) { b = hi; c = std::max(lo, b - h); a = std::max(lo, c - h); }
      Xr(i, 0) = a; Xr(i, 1) = c; Xr(i, 2) = b;
    }
    NumericMatrix gr = marginal_eval_nodes(PM, Xr, t0col, data, constants,
                                           designs, type, bounds, transforms,
                                           pretransforms, p_types, min_ll, trend);
    int bad = 0;
    for (int i = 0; i < nn; ++i) {
      double m_i, s_i;
      const double g1 = gr(i, 0) + R::dnorm(Xr(i, 0), mu, sigma, 1);
      const double g2 = gr(i, 1) + R::dnorm(Xr(i, 1), mu, sigma, 1);
      const double g3 = gr(i, 2) + R::dnorm(Xr(i, 2), mu, sigma, 1);
      ft[static_cast<size_t>(i)] = 0;
      if (marginal_parabola(Xr(i, 0), g1, Xr(i, 1), g2, Xr(i, 2), g3, &m_i, &s_i)) {
        // Let the window contract freely but grow only gradually, so one bad
        // fit cannot throw the rule back out to the pilot scale.
        const double cap = std::min(max_sd, 4.0 * sv[static_cast<size_t>(i)]);
        const double mc = std::min(hi, std::max(lo, m_i));
        const double dev = mc - Xr(i, 1);
        // A mode pinned on a boundary can never sit at the probe centre (the
        // bracket is shifted inwards to stay feasible), so a flush bracket
        // pointing further out is consistent, not a failed fit.
        const bool ok_local = std::fabs(dev) <= 0.5 * (Xr(i, 2) - Xr(i, 1)) ||
                              (Xr(i, 0) <= lo && dev < 0.0) ||
                              (Xr(i, 2) >= hi && dev > 0.0);
        if (!ok_local) ++bad;
        else ft[static_cast<size_t>(i)] = 1;
        md[static_cast<size_t>(i)] = mc;
        sv[static_cast<size_t>(i)] = std::min(cap, s_i);
      } else {
        ++bad;
      }
    }
    return bad;
  };
  auto refine_round = [&]() -> int { return refine_on(particle_matrix, mode, sdev, fitted); };

  // Coarse uniform scan of the whole feasible interval: brackets each particle's
  // own mode.  The log integrand is smooth and unimodal in log t0, so the
  // bracket is reliable even when the exponentiated peak is far thinner than the
  // scan spacing -- but it costs n_scan full passes over the batch, which is why
  // it is worth running on a subset and predicting the rest.
  auto pilot_on = [&](const NumericMatrix& PM, std::vector<double>& md,
                      std::vector<double>& sv, std::vector<char>& ft) {
    const int nn = PM.nrow();
    NumericMatrix Xs(nn, n_scan);
    std::vector<double> xs(static_cast<size_t>(n_scan));
    for (int k = 0; k < n_scan; ++k) {
      xs[static_cast<size_t>(k)] = lo + (hi - lo) * k / (n_scan - 1);
      for (int i = 0; i < nn; ++i) Xs(i, k) = xs[static_cast<size_t>(k)];
    }
    NumericMatrix gs = marginal_eval_nodes(PM, Xs, t0col, data, constants,
                                           designs, type, bounds, transforms,
                                           pretransforms, p_types, min_ll, trend);
    for (int k = 0; k < n_scan; ++k) {
      const double lp = R::dnorm(xs[static_cast<size_t>(k)], mu, sigma, 1);
      for (int i = 0; i < nn; ++i) gs(i, k) += lp;
    }
    for (int i = 0; i < nn; ++i) {
      int best = 0;
      double bv = R_NegInf;
      for (int k = 0; k < n_scan; ++k) if (gs(i, k) > bv) { bv = gs(i, k); best = k; }
      const int j = std::min(std::max(best, 1), n_scan - 2);
      double m_i, s_i;
      if (marginal_parabola(xs[static_cast<size_t>(j - 1)], gs(i, j - 1),
                            xs[static_cast<size_t>(j)],     gs(i, j),
                            xs[static_cast<size_t>(j + 1)], gs(i, j + 1), &m_i, &s_i)) {
        ft[static_cast<size_t>(i)] = 1;
        md[static_cast<size_t>(i)] = std::min(hi, std::max(lo, m_i));
        sv[static_cast<size_t>(i)] = std::min(max_sd, s_i);
      }
    }
  };

  // Particles whose window came from a HINT rather than from the pilot scan:
  // only those can be wrong in a way the pilot would fix, so only those are
  // worth verifying against the finished grid.
  std::vector<char> hinted(static_cast<size_t>(np), 0);
  if (adaptive && n_nodes >= 6) {
    bool warm_ok = R_FINITE(warm_mode) && R_FINITE(warm_sd) && warm_sd > 0.0 &&
                   warm_mode >= lo && warm_mode <= hi;
    if (warm_ok) {
      const double h0 = std::min(max_sd, std::max(warm_sd * warm_inflate, 1e-8));
      for (int i = 0; i < np; ++i) {
        mode[static_cast<size_t>(i)] = warm_mode;
        sdev[static_cast<size_t>(i)] = h0;
        fitted[static_cast<size_t>(i)] = 1;
      }
      // One probe round is enough when every particle's mode sits inside the
      // hinted bracket. Otherwise give the batch a second round at the rescaled
      // spacing -- unless the hint was so far off that most particles missed,
      // in which case chasing it costs more than the pilot it replaces.
      int bad = refine_round();
      if (bad > 0 && bad * 4 <= np) {
        for (int r = 1; r < n_refine_warm && bad > 0; ++r) bad = refine_round();
      }
      g.warm_used = static_cast<double>(np - bad) / static_cast<double>(np);
      if (bad == 0) {
        std::fill(hinted.begin(), hinted.end(), 1);
      } else if (bad * 4 <= np) {
        std::vector<int> bad_rows;
        for (int i = 0; i < np; ++i) {
          if (fitted[static_cast<size_t>(i)]) hinted[static_cast<size_t>(i)] = 1;
          else bad_rows.push_back(i);
        }
        const int nb = static_cast<int>(bad_rows.size());
        NumericMatrix Pb(nb, n_par);
        colnames(Pb) = p_names;
        for (int j = 0; j < nb; ++j)
          for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(bad_rows[static_cast<size_t>(j)], c);
        std::vector<double> mb(static_cast<size_t>(nb), def_mode);
        std::vector<double> sb(static_cast<size_t>(nb), def_sd);
        std::vector<char> fb(static_cast<size_t>(nb), 0);
        pilot_on(Pb, mb, sb, fb);
        refine_on(Pb, mb, sb, fb);
        for (int j = 0; j < nb; ++j) {
          const int i = bad_rows[static_cast<size_t>(j)];
          mode[static_cast<size_t>(i)] = mb[static_cast<size_t>(j)];
          sdev[static_cast<size_t>(i)] = sb[static_cast<size_t>(j)];
          fitted[static_cast<size_t>(i)] = fb[static_cast<size_t>(j)];
        }
      } else {
        warm_ok = false;
        g.warm_used = 0;
        std::fill(mode.begin(), mode.end(), def_mode);
        std::fill(sdev.begin(), sdev.end(), def_sd);
        std::fill(fitted.begin(), fitted.end(), 0);
      }
    }
    if (!warm_ok) {
      // Within-iteration hint: the conditional t0 mode is a smooth function of
      // the OTHER sampled parameters, so pilot-scan a small spread-out subset of
      // the batch, regress its modes on those parameters, and use the fit to
      // bracket every remaining particle.  Unlike the cross-iteration warm start
      // this hint is per particle and never stale -- it tracks the batch it was
      // built from -- and it replaces n_scan full passes with n_scan*k/np.
      bool pred_ok = false;
      if (predict_mode) {
        // Predictors: every sampled column that actually varies here (a constant
        // column, e.g. another component's parameters, carries no information
        // and would make the design rank-deficient).
        std::vector<int> cols;
        for (int j = 0; j < n_par; ++j) {
          if (j == t0col) continue;
          double mn = particle_matrix(0, j), mx = mn;
          for (int i = 1; i < np; ++i) {
            const double v = particle_matrix(i, j);
            if (v < mn) mn = v; else if (v > mx) mx = v;
          }
          if (mx - mn > 0.0) cols.push_back(j);
        }
        const int p_pred = static_cast<int>(cols.size());
        const int k_sub = std::max(p_pred + 6, np / pred_frac);
        // Only worth it when the subset is a small fraction of the batch.
        if (p_pred > 0 && np >= 2 * k_sub) {
          NumericMatrix Ps(k_sub, n_par);
          colnames(Ps) = p_names;
          std::vector<int> S(static_cast<size_t>(k_sub));
          for (int j = 0; j < k_sub; ++j) {
            // Evenly spaced: the batch is laid out in proposal-type blocks, so
            // a spread subset spans all of them.
            S[static_cast<size_t>(j)] = static_cast<int>(
              (static_cast<double>(j) * np) / k_sub);
            for (int c = 0; c < n_par; ++c)
              Ps(j, c) = particle_matrix(S[static_cast<size_t>(j)], c);
          }
          std::vector<double> ms(static_cast<size_t>(k_sub), def_mode);
          std::vector<double> ss(static_cast<size_t>(k_sub), def_sd);
          std::vector<char> fs(static_cast<size_t>(k_sub), 0);
          pilot_on(Ps, ms, ss, fs);
          refine_on(Ps, ms, ss, fs);
          // Regress the resolved subset modes on their parameters.
          std::vector<double> yv;
          std::vector<int> rows;
          for (int j = 0; j < k_sub; ++j)
            if (fs[static_cast<size_t>(j)]) { rows.push_back(j); yv.push_back(ms[static_cast<size_t>(j)]); }
          const int n_ok = static_cast<int>(rows.size());
          std::vector<double> Xv(static_cast<size_t>(n_ok) * p_pred);
          for (int c = 0; c < p_pred; ++c)
            for (int i = 0; i < n_ok; ++i)
              Xv[static_cast<size_t>(c) * n_ok + i] = Ps(rows[static_cast<size_t>(i)], cols[static_cast<size_t>(c)]);
          std::vector<double> beta;
          if (marginal_lsfit(Xv, yv, n_ok, p_pred, beta)) {
            // Residual scale sets the probe half-width: the bracket has to be
            // wide enough to contain the prediction error, and the refinement
            // rounds contract from there onto each particle's own peak.
            // Robust scale: the proposal batch is a MIXTURE, and its few very
            // wide draws are outliers to the local linear fit.  A plain RMS lets
            // them set the probe width for the entire batch (and then the whole
            // hint gets rejected as uninformative); a median absolute residual
            // sizes the probe for the bulk, and the verification pass below
            // catches the handful the regression could not place.
            std::vector<double> absres(static_cast<size_t>(n_ok));
            for (int i = 0; i < n_ok; ++i) {
              double f = beta[0];
              for (int c = 0; c < p_pred; ++c) f += beta[static_cast<size_t>(c + 1)] * Xv[static_cast<size_t>(c) * n_ok + i];
              absres[static_cast<size_t>(i)] = std::fabs(yv[static_cast<size_t>(i)] - f);
            }
            std::sort(absres.begin(), absres.end());
            const double rms = 1.4826 * absres[absres.size() / 2];
            std::vector<double> sorted;
            for (int j = 0; j < k_sub; ++j) if (fs[static_cast<size_t>(j)]) sorted.push_back(ss[static_cast<size_t>(j)]);
            std::sort(sorted.begin(), sorted.end());
            const double med_sd = sorted.empty() ? def_sd : sorted[sorted.size() / 2];
            const double h0 = std::min(max_sd, std::max(std::max(2.5 * rms, 4.0 * med_sd), 1e-8));
            // A prediction no sharper than the interval itself is worthless.
            if (R_FINITE(rms) && h0 < 0.5 * (hi - lo)) {
              for (int i = 0; i < np; ++i) {
                double f = beta[0];
                for (int c = 0; c < p_pred; ++c)
                  f += beta[static_cast<size_t>(c + 1)] * particle_matrix(i, cols[static_cast<size_t>(c)]);
                if (!R_FINITE(f)) { f = def_mode; }
                mode[static_cast<size_t>(i)] = std::min(hi, std::max(lo, f));
                sdev[static_cast<size_t>(i)] = h0;
                fitted[static_cast<size_t>(i)] = 1;
              }
              std::fill(hinted.begin(), hinted.end(), 1);
              // Keep the subset's own resolved fits: they are exact, not predicted.
              for (int j = 0; j < k_sub; ++j) if (fs[static_cast<size_t>(j)]) {
                hinted[static_cast<size_t>(S[static_cast<size_t>(j)])] = 0;
                const int i = S[static_cast<size_t>(j)];
                mode[static_cast<size_t>(i)] = ms[static_cast<size_t>(j)];
                sdev[static_cast<size_t>(i)] = ss[static_cast<size_t>(j)];
              }
              pred_ok = true;
              g.pred_used = 1;
            }
          }
        }
      }
      int bad = 0;
      if (pred_ok) {
        for (int r = 0; r < n_refine_pred; ++r)
          bad = refine_on(particle_matrix, mode, sdev, fitted);
        // The prediction is only a bracket; if it failed to bracket a
        // meaningful share of the batch, fall back to the full pilot.
        if (bad * 4 > np) {
          pred_ok = false;
          g.pred_used = 0;
          std::fill(hinted.begin(), hinted.end(), 0);
          std::fill(mode.begin(), mode.end(), def_mode);
          std::fill(sdev.begin(), sdev.end(), def_sd);
          std::fill(fitted.begin(), fitted.end(), 0);
        }
      }
      if (!pred_ok) {
        std::fill(hinted.begin(), hinted.end(), 0);
        pilot_on(particle_matrix, mode, sdev, fitted);
        for (int r = 0; r < n_refine; ++r) refine_round();
      }
    }
  }

  // --- per-particle composite rule: [lo, core) core [core, hi] --------------
  // Records whether the rule actually spends a node on each tail, which is what
  // makes the finished grid self-checking below.
  std::vector<char> has_left(static_cast<size_t>(np), 0), has_right(static_cast<size_t>(np), 0);
  std::vector<int> core_lo(static_cast<size_t>(np), 0), core_hi(static_cast<size_t>(np), n_nodes - 1);
  NumericMatrix logw(np, n_nodes);   // quadrature weights, so g = log_terms - logw
  auto build_rules = [&](const std::vector<int>& idxs) {
    for (size_t t = 0; t < idxs.size(); ++t) {
      const int i = idxs[t];
      double clo = lo, chi = hi;
      if (fitted[static_cast<size_t>(i)]) {
        const double half = window_sd * sdev[static_cast<size_t>(i)];
        clo = std::max(lo, mode[static_cast<size_t>(i)] - half);
        chi = std::min(hi, mode[static_cast<size_t>(i)] + half);
        if (!(chi > clo)) { clo = lo; chi = hi; }
        g.mode[i] = mode[static_cast<size_t>(i)];
        g.sd[i] = sdev[static_cast<size_t>(i)];
      } else {
        g.mode[i] = NA_REAL;
        g.sd[i] = NA_REAL;
      }
      const int n_left  = (clo > lo && n_nodes >= 4) ? 1 : 0;
      const int n_right = (chi < hi && n_nodes >= 4) ? 1 : 0;
      if (n_left == 0) clo = lo;
      if (n_right == 0) chi = hi;
      has_left[static_cast<size_t>(i)] = static_cast<char>(n_left);
      has_right[static_cast<size_t>(i)] = static_cast<char>(n_right);
      core_lo[static_cast<size_t>(i)] = n_left;
      core_hi[static_cast<size_t>(i)] = n_nodes - 1 - n_right;
      MarginalRule rule;
      marginal_append_gl_panel(rule, lo, clo, n_left);
      marginal_append_gl_panel(rule, clo, chi, n_nodes - n_left - n_right);
      marginal_append_gl_panel(rule, chi, hi, n_right);
      for (int k = 0; k < n_nodes; ++k) {
        const double xk = rule.x[static_cast<size_t>(k)];
        g.nodes(i, k) = xk;
        logw(i, k) = std::log(rule.w[static_cast<size_t>(k)]);
        g.log_terms(i, k) = logw(i, k) + R::dnorm(xk, mu, sigma, 1);
      }
    }
  };
  // Evaluate the kernel on a set of rows at their own nodes and add it in.
  auto eval_rows = [&](const std::vector<int>& idxs) {
    const int nb = static_cast<int>(idxs.size());
    if (nb == 0) return;
    NumericMatrix Pb(nb, n_par), Xb(nb, n_nodes);
    colnames(Pb) = p_names;
    for (int j = 0; j < nb; ++j) {
      for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(idxs[static_cast<size_t>(j)], c);
      for (int k = 0; k < n_nodes; ++k) Xb(j, k) = g.nodes(idxs[static_cast<size_t>(j)], k);
    }
    NumericMatrix lb = marginal_eval_nodes(Pb, Xb, t0col, data, constants, designs, type,
                                           bounds, transforms, pretransforms, p_types,
                                           min_ll, trend);
    for (int j = 0; j < nb; ++j)
      for (int k = 0; k < n_nodes; ++k)
        g.log_terms(idxs[static_cast<size_t>(j)], k) += lb(j, k);
  };

  std::vector<int> all_rows(static_cast<size_t>(np));
  for (int i = 0; i < np; ++i) all_rows[static_cast<size_t>(i)] = i;
  build_rules(all_rows);

  // --- one batched pass over the whole np x K grid --------------------------
  NumericMatrix llf = marginal_eval_nodes(particle_matrix, g.nodes, t0col, data, constants,
                                          designs, type, bounds, transforms,
                                          pretransforms, p_types, min_ll, trend);
  for (int k = 0; k < n_nodes; ++k)
    for (int i = 0; i < np; ++i) g.log_terms(i, k) += llf(i, k);

  // --- verify the finished rule, and repair whatever it missed --------------
  // The mode search can be fooled: a hinted bracket that never contained the
  // peak, or a parabola fitted on a shoulder, yields a window that is either
  // off the peak or far too wide for it.  The finished grid detects both for
  // free.  Reading the integrand back off its own nodes (g = log_terms - log w),
  // a rule that resolved its peak has an interior maximum whose curvature
  // implies a width comparable to the node spacing.  A maximum at an end node,
  // or an implied width far BELOW the spacing (the tell-tale of a wide window
  // straddling a razor-thin peak, which silently loses most of the mass), means
  // the rule is not integrating what it should.  Those particles -- and only
  // those -- are refitted from the full pilot scan and re-integrated, so a hint
  // can cost time but never accuracy.
  if (adaptive && n_nodes >= 6) {
    std::vector<int> hurt;
    for (int i = 0; i < np; ++i) {
      // A window that already came from the pilot scan has nothing to fall back
      // to: refitting it would repeat the identical deterministic search.
      if (!hinted[static_cast<size_t>(i)]) continue;
      int best = 0;
      double bv = R_NegInf;
      for (int k = 0; k < n_nodes; ++k) if (g.log_terms(i, k) > bv) { bv = g.log_terms(i, k); best = k; }
      if (!R_FINITE(bv)) continue;
      const int clo_i = core_lo[static_cast<size_t>(i)], chi_i = core_hi[static_cast<size_t>(i)];
      // Peak at an end node, or outside the core panel: the window is off-peak.
      if (best <= clo_i || best >= chi_i) { hurt.push_back(i); continue; }
      const double x1 = g.nodes(i, best - 1), x2 = g.nodes(i, best), x3 = g.nodes(i, best + 1);
      double m_i, s_i;
      if (!marginal_parabola(x1, g.log_terms(i, best - 1) - logw(i, best - 1),
                             x2, g.log_terms(i, best) - logw(i, best),
                             x3, g.log_terms(i, best + 1) - logw(i, best + 1), &m_i, &s_i)) {
        hurt.push_back(i);
        continue;
      }
      // Gauss-Legendre nodes cluster towards the panel edges, so on a CORRECT
      // rule (core = mode +/- window_sd sigma) the spacing either side of the
      // peak is already a couple of sigma; only a width far below that means the
      // nodes are straddling a peak they never resolved.
      if (s_i < 0.2 * 0.5 * (x3 - x1)) hurt.push_back(i);   // peak under-resolved
    }
    g.repaired = static_cast<int>(hurt.size());
    if (!hurt.empty()) {
      const int nb = g.repaired;
      NumericMatrix Pb(nb, n_par);
      colnames(Pb) = p_names;
      for (int j = 0; j < nb; ++j)
        for (int c = 0; c < n_par; ++c) Pb(j, c) = particle_matrix(hurt[static_cast<size_t>(j)], c);
      std::vector<double> mb(static_cast<size_t>(nb), def_mode);
      std::vector<double> sb(static_cast<size_t>(nb), def_sd);
      std::vector<char> fb(static_cast<size_t>(nb), 0);
      pilot_on(Pb, mb, sb, fb);
      for (int r = 0; r < n_refine; ++r) refine_on(Pb, mb, sb, fb);
      for (int j = 0; j < nb; ++j) {
        const int i = hurt[static_cast<size_t>(j)];
        mode[static_cast<size_t>(i)] = mb[static_cast<size_t>(j)];
        sdev[static_cast<size_t>(i)] = sb[static_cast<size_t>(j)];
        fitted[static_cast<size_t>(i)] = fb[static_cast<size_t>(j)];
      }
      build_rules(hurt);
      eval_rows(hurt);
    }
  }
  return g;
}

// Marginal log-likelihood: log-sum-exp of the shared core's per-node terms.
// (Not Rcpp-exported: internal helper, reached only via calc_ll_oo.)
NumericVector calc_ll_oo_marginal(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
    Rcpp::List marginalise) {
  MarginalGrid g = calc_ll_oo_marginal_core(particle_matrix, data, constants, designs, type,
      bounds, transforms, pretransforms, p_types, min_ll, trend, marginalise);
  const int np = g.log_terms.nrow(), nn = g.log_terms.ncol();
  NumericVector out(np);
  for (int i = 0; i < np; ++i) {
    double m = R_NegInf;
    for (int k = 0; k < nn; ++k) if (g.log_terms(i, k) > m) m = g.log_terms(i, k);
    if (!R_FINITE(m)) { out[i] = R_NegInf; continue; }
    double s = 0.0;
    for (int k = 0; k < nn; ++k) s += std::exp(g.log_terms(i, k) - m);
    out[i] = m + std::log(s);
  }
  return out;
}

// Storage-time accessor for the t0 reconstruction (reconstruct-at-storage design,
// plan "use the posterior as sampled"). Returns the quadrature node grid and the
// unnormalized per-particle log-terms; the sampler forms w = softmax_k(log_terms)
// and draws one node per subject to write into the stored alpha[t0], so the
// posterior carries a valid t0 draw and predict()/make_data() stay unchanged.
// Nodes are on the sampled (log-t0) scale, matching alpha.
// [[Rcpp::export]]
List calc_ll_oo_marginal_nodes(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::List marginalise,
    Rcpp::Nullable<Rcpp::List> trend = R_NilValue) {
  MarginalGrid g = calc_ll_oo_marginal_core(particle_matrix, data, constants, designs, type,
      bounds, transforms, pretransforms, p_types, min_ll, trend, marginalise);
  const int np = g.log_terms.nrow(), nn = g.log_terms.ncol();
  NumericVector ll(np);
  for (int i = 0; i < np; ++i) {
    double m = R_NegInf;
    for (int k = 0; k < nn; ++k) if (g.log_terms(i, k) > m) m = g.log_terms(i, k);
    if (!R_FINITE(m)) { ll[i] = R_NegInf; continue; }
    double s = 0.0;
    for (int k = 0; k < nn; ++k) s += std::exp(g.log_terms(i, k) - m);
    ll[i] = m + std::log(s);
  }
  return List::create(Rcpp::Named("nodes") = g.nodes,
                      Rcpp::Named("log_terms") = g.log_terms,
                      Rcpp::Named("ll") = ll,
                      Rcpp::Named("mode") = g.mode,
                      Rcpp::Named("sd") = g.sd,
                      Rcpp::Named("warm_used") = g.warm_used,
                      Rcpp::Named("pred_used") = g.pred_used,
                      Rcpp::Named("repaired") = g.repaired);
}
