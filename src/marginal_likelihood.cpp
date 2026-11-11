#include "gh_quad.h"
#include "marginal_likelihood.h"

#include "gl_quad.h"
#include "rt_bins.h"

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

// Per-node grid of a marginalization: the quadrature node values (sampled
// scale) and the unnormalized per-particle log-terms
//   ell_ik = log quadrature_weight_ik + log p(x_ik|eta) + log L_i(x_ik)
// that BOTH the marginal-ll reducer and the storage-time weight accessor share,
// so the node grid and the reconstruction weights are guaranteed identical.
// Nodes are PER PARTICLE (np x K): each proposal gets a rule centred on its own
// conditional posterior, so one particle's rule can never be dragged off its
// peak by the rest of a scattered batch.
struct MarginalGrid {
  Rcpp::NumericMatrix nodes;      // np x W: per-particle node values (NA where unused)
  Rcpp::NumericMatrix log_terms;  // np x W: log weight + log prior + log likelihood
  Rcpp::NumericVector mode;       // np: located peak of the log integrand
  Rcpp::NumericVector sd;         // np: its curvature scale
  Rcpp::IntegerVector rounds;     // np: kernel passes the particle needed
  int gl_used = 0;                // particles integrated by the edge rule
  int screened = 0;               // particles left on their first rule
};

// Nodes and weights of a Gauss-Legendre rule over one interval.
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
    const NumericMatrix& particle_matrix, const NumericMatrix& X, int pcol,
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
        pm(r, pcol) = X(i, k0 + k);
      }
    }
    // No varying mask: the quadrature rewrites the marginalised coordinate for
    // every node, so which of the caller's coordinates were blocked says
    // nothing about which of these vary.
    NumericVector ll = calc_ll_oo(pm, data, constants, designs, type, bounds,
                                  transforms, pretransforms, p_types, min_ll,
                                  trend, R_NilValue, R_NilValue, R_NilValue);
    for (int k = 0; k < mk; ++k)
      for (int i = 0; i < np; ++i) out(i, k0 + k) = ll[k * np + i];
  }
  return out;
}

// log sum_k exp(m(i, k)) over the first n columns of row i.
static double log_sum_exp_row(const NumericMatrix& m, int i, int n) {
  double mx = R_NegInf;
  for (int k = 0; k < n; ++k) if (m(i, k) > mx) mx = m(i, k);
  if (!R_FINITE(mx)) return mx;
  double acc = 0.0;
  for (int k = 0; k < n; ++k) acc += std::exp(m(i, k) - mx);
  return mx + std::log(acc);
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

// Column of the marginalised parameter in the proposal matrix.
static int marginal_param_col(const CharacterVector& p_names, const std::string& mparam) {
  for (int j = 0; j < p_names.size(); ++j)
    if (Rcpp::as<std::string>(p_names[j]) == mparam) return j;
  Rcpp::stop("calc_ll_oo marginalise: '%s' is not a sampled column.", mparam);
  return -1;
}

// A natural-scale value of parameter `p` on its sampled axis: the inverse of
// the model transform (natural = lower + exp(x), lower + (upper - lower) *
// Phi(x), or x).
static double marginal_to_sampled(const List& transforms, const std::string& p,
                                  double nat) {
  std::string f = "identity";
  double lw = 0.0, up = 1.0;
  auto named = [&](const char* field, double* out) {
    if (!transforms.containsElementNamed(field)) return;
    NumericVector v = transforms[field];
    if (!v.hasAttribute("names")) return;
    CharacterVector nm = v.names();
    for (int j = 0; j < v.size(); ++j)
      if (Rcpp::as<std::string>(nm[j]) == p) { *out = v[j]; return; }
  };
  if (transforms.containsElementNamed("func")) {
    CharacterVector fn = transforms["func"];
    if (fn.hasAttribute("names")) {
      CharacterVector nm = fn.names();
      for (int j = 0; j < fn.size(); ++j)
        if (Rcpp::as<std::string>(nm[j]) == p) f = Rcpp::as<std::string>(fn[j]);
    }
  }
  named("lower", &lw);
  named("upper", &up);
  if (f == "exp") return nat > lw ? std::log(nat - lw) : R_NegInf;
  if (f == "pnorm") {
    if (!(up > lw)) return NA_REAL;
    return R::qnorm((nat - lw) / (up - lw), 0.0, 1.0, 1, 0);
  }
  return nat;
}

// Feasible interval of the marginalised parameter on its sampled axis: its
// natural-scale bounds (outside them c_do_bound floors the whole subject to
// n_trials*min_ll), and for t0 also the fastest response, above which that
// response is infeasible. A local guess clock or an estimated uniform-guess
// mixture removes that ceiling. Otherwise the ceiling is the retained bin
// top, without discarding the final millisecond of its feasible support.
static void marginal_feasible(DataFrame data, NumericVector constants, String type,
                              List bounds, List transforms,
                              const CharacterVector& p_names,
                              const std::string& mparam, double eps,
                              double* lo_out, double* hi_out) {
  double nat_lo = R_NegInf, nat_hi = R_PosInf;
  {
    NumericMatrix minmax = bounds["minmax"];
    CharacterVector mm_names = colnames(minmax);
    for (int j = 0; j < mm_names.size(); ++j)
      if (Rcpp::as<std::string>(mm_names[j]) == mparam) {
        nat_lo = minmax(0, j);
        nat_hi = minmax(1, j);
        break;
      }
  }
  if (mparam == "t0") {
    NumericVector rt = data["rt"];
    double mn = R_PosInf;
    for (int k = 0; k < rt.size(); ++k)
      if (R_FINITE(rt[k]) && rt[k] < mn) mn = rt[k];
    const double h = rt_bin_width(data);
    if (h > 0.0) {
      RtBinCache bins; bins.init(data, 1);
      mn = R_PosInf;
      for (int k = 0; k < rt.size(); ++k)
        if (R_FINITE(rt[k])) mn = std::min(mn, bins.bins[k].hi);
    }
    bool guesses = false;
    for (int j = 0; j < p_names.size(); ++j)
      if (Rcpp::as<std::string>(p_names[j]) == "pGuess") guesses = true;
    if (constants.hasAttribute("names")) {
      CharacterVector cn = constants.attr("names");
      for (int j = 0; j < cn.size(); ++j)
        if (Rcpp::as<std::string>(cn[j]) == "pGuess" && constants[j] != R_NegInf)
          guesses = true;
    }
    const std::string type_name = Rcpp::as<std::string>(Rcpp::wrap(type));
    guesses = guesses || type_name.find("LOCAL_GUESS") != std::string::npos ||
      type_name.find("LOCAL_KILL_GUESS") != std::string::npos;
    if (R_FINITE(mn) && !guesses) {
      double top = h > 0.0 ? std::nextafter(mn, 0.0) : mn - eps;
      if (top < std::numeric_limits<double>::min()) top = std::numeric_limits<double>::min();
      nat_hi = std::min(nat_hi, top);
    }
  }
  double lo = R_FINITE(nat_lo) ? marginal_to_sampled(transforms, mparam, nat_lo) : R_NegInf;
  double hi = R_FINITE(nat_hi) ? marginal_to_sampled(transforms, mparam, nat_hi) : R_PosInf;
  if (ISNAN(lo)) lo = R_NegInf;
  if (ISNAN(hi)) hi = R_PosInf;
  *lo_out = lo;
  *hi_out = hi;
}

// Marginalization core: integrate one shared parameter (e.g. t0) out of the
// complete subject likelihood, per particle, on its sampled axis. The log
// integrand g(x) = log L(x) + log N(x | mu_i, sigma) is close to quadratic
// around its peak, so the integral is a K-node Gauss-Hermite rule centred on
// the peak at its curvature scale. Each particle starts from its predicted
// peak `centre` and scale `probe_sd` (default: its prior). A rule is kept
// when a quadratic fitted through its own nodes puts the peak within half a
// scale of the centre and the curvature scale within a quarter of the rule's;
// otherwise it is placed again at that fitted peak while the placements
// converge. Particles whose rules do not converge go to probe triples, which
// locate the peak for a final rule there; their steps are limited in length
// and kept inside a bracket on the peak. A peak the triples do not locate is
// integrated by Gauss-Legendre around the best node evaluated. A peak within
// reach of a hard edge of the feasible interval is truncated, so it is
// integrated by Gauss-Legendre up to the edge instead. Given `log_offset`, the rest of each
// particle's selection weight, particles that cannot be selected keep their
// first rule. Otherwise every particle is treated on its own, so its value is
// a fixed function of that particle and its start. The registered kernel is reused verbatim (nodes
// overwrite the sampled column), so every race model is covered.
static MarginalGrid calc_ll_oo_marginal_core(
    NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
    List designs, String type, List bounds, List transforms, List pretransforms,
    CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
    Rcpp::List marginalise) {
  const int np = particle_matrix.nrow();
  const int n_par = particle_matrix.ncol();

  // --- marginalization spec: sampled x ~ N(mu_i, sigma) for particle i ---
  // `mu` is one value per particle (the conditional group mean given the
  // particle's other parameters) or a single value shared by all of them.
  const std::string mparam = Rcpp::as<std::string>(marginalise["param"]);
  const NumericVector mu_in = marginalise["mu"];
  const double sigma = Rcpp::as<double>(marginalise["sigma"]);
  if (mu_in.size() != 1 && mu_in.size() != np)
    Rcpp::stop("calc_ll_oo marginalise: 'mu' needs one value or one per particle.");
  if (!R_FINITE(sigma) || !(sigma > 0.0))
    Rcpp::stop("calc_ll_oo marginalise: 'sigma' must be positive.");
  auto mu_of = [&](int i) { return mu_in[mu_in.size() == 1 ? 0 : i]; };
  int n_nodes = marginalise.containsElementNamed("n_nodes")
                  ? Rcpp::as<int>(marginalise["n_nodes"]) : 5;
  if (n_nodes < 3) n_nodes = 3;
  // Per-particle starting centres and scales, when the caller can predict
  // each peak; otherwise each rule starts at the particle's prior.
  const NumericVector centre_in = marginalise.containsElementNamed("centre")
                  ? NumericVector(marginalise["centre"]) : NumericVector(0);
  if (centre_in.size() != 0 && centre_in.size() != np)
    Rcpp::stop("calc_ll_oo marginalise: 'centre' needs one value per particle.");
  const NumericVector probe_in = marginalise.containsElementNamed("probe_sd")
                  ? NumericVector(marginalise["probe_sd"]) : NumericVector(0);
  if (probe_in.size() > 1 && probe_in.size() != np)
    Rcpp::stop("calc_ll_oo marginalise: 'probe_sd' needs one value or one per particle.");
  // The rest of each particle's log weight in the caller's selection, when
  // the caller selects among the particles by weight.
  const NumericVector offset_in = marginalise.containsElementNamed("log_offset")
                  ? NumericVector(marginalise["log_offset"]) : NumericVector(0);
  if (offset_in.size() != 0 && offset_in.size() != np)
    Rcpp::stop("calc_ll_oo marginalise: 'log_offset' needs one value per particle.");
  int max_rounds = marginalise.containsElementNamed("max_rounds")
                  ? Rcpp::as<int>(marginalise["max_rounds"]) : 12;
  if (max_rounds < 1) max_rounds = 1;
  // Half-width, in SDs, of the Gauss-Legendre window around a truncated or
  // unresolved peak.
  const double span = marginalise.containsElementNamed("span")
                  ? Rcpp::as<double>(marginalise["span"]) : 6.0;
  const double eps = marginalise.containsElementNamed("eps")
                  ? Rcpp::as<double>(marginalise["eps"]) : 1e-3;
  // Rule acceptance: the quadrature error of a Gauss-Hermite rule grows with
  // its centre's offset from the peak and its scale's mismatch with the
  // curvature; within these limits it stays a few thousandths of a nat.
  const double max_offset = 0.5, max_log_scale = std::log(1.25);

  CharacterVector p_names = colnames(particle_matrix);
  const int pcol = marginal_param_col(p_names, mparam);
  double x_lo, x_hi;
  marginal_feasible(data, constants, type, bounds, transforms, p_names, mparam, eps,
                    &x_lo, &x_hi);

  const GHRule& gh = gh_rule(n_nodes);
  std::vector<double> z(static_cast<size_t>(n_nodes)), lwz(static_cast<size_t>(n_nodes));
  double zmax = 0.0;
  for (int k = 0; k < n_nodes; ++k) {
    z[static_cast<size_t>(k)] = std::sqrt(2.0) * gh.x[static_cast<size_t>(k)];
    lwz[static_cast<size_t>(k)] = std::log(gh_standard_normal_weight(gh, k));
    zmax = std::max(zmax, std::fabs(z[static_cast<size_t>(k)]));
  }
  const int n_gl = 2 * n_nodes;   // a truncated peak is not Gaussian-shaped
  const int width = std::max(n_nodes, n_gl);

  MarginalGrid g;
  g.nodes = NumericMatrix(np, width);
  g.log_terms = NumericMatrix(np, width);
  std::fill(g.nodes.begin(), g.nodes.end(), NA_REAL);
  std::fill(g.log_terms.begin(), g.log_terms.end(), R_NegInf);
  g.mode = NumericVector(np, NA_REAL);
  g.sd = NumericVector(np, NA_REAL);
  g.rounds = IntegerVector(np, 0);
  if (np == 0) return g;

  if (!(x_hi > x_lo)) {
    // The data leave no feasible value at all: the likelihood is at its floor
    // everywhere, so the integral is the floored likelihood, taken at each
    // particle's prior centre as a single node.
    NumericMatrix X1(np, 1);
    for (int i = 0; i < np; ++i) X1(i, 0) = mu_of(i);
    NumericMatrix l1 = marginal_eval_nodes(particle_matrix, X1, pcol, data, constants,
                                           designs, type, bounds, transforms,
                                           pretransforms, p_types, min_ll, trend);
    for (int i = 0; i < np; ++i) {
      g.nodes(i, 0) = X1(i, 0);
      g.log_terms(i, 0) = l1(i, 0);
    }
    return g;
  }

  auto eval = [&](const std::vector<int>& rows, const NumericMatrix& X) {
    const int nb = static_cast<int>(rows.size());
    NumericMatrix Pb(nb, n_par);
    colnames(Pb) = p_names;
    for (int j = 0; j < nb; ++j)
      for (int q = 0; q < n_par; ++q) Pb(j, q) = particle_matrix(rows[static_cast<size_t>(j)], q);
    NumericMatrix out = marginal_eval_nodes(Pb, X, pcol, data, constants, designs, type,
                                            bounds, transforms, pretransforms, p_types,
                                            min_ll, trend);
    for (int j = 0; j < nb; ++j) {
      const double mi = mu_of(rows[static_cast<size_t>(j)]);
      for (int k = 0; k < X.ncol(); ++k) out(j, k) += R::dnorm(X(j, k), mi, sigma, 1);
    }
    return out;
  };
  auto clamp = [&](double x) { return std::min(x_hi, std::max(x_lo, x)); };
  auto rule_inside = [&](double c, double h) {
    return c - zmax * h > x_lo && c + zmax * h < x_hi;
  };
  auto scale_ok = [&](double h, double s) { return std::fabs(std::log(h / s)) <= max_log_scale; };

  std::vector<double> c(static_cast<size_t>(np)), h(static_cast<size_t>(np), sigma);
  std::vector<double> v(static_cast<size_t>(np), NA_REAL), sv(static_cast<size_t>(np), NA_REAL);
  std::vector<char> fitted(static_cast<size_t>(np), 0), done(static_cast<size_t>(np), 0);
  // The best node evaluated so far, and an interval known to hold the peak
  // of a unimodal integrand.
  std::vector<double> xb(static_cast<size_t>(np), NA_REAL), gb(static_cast<size_t>(np), R_NegInf);
  std::vector<double> br_lo(static_cast<size_t>(np), x_lo), br_hi(static_cast<size_t>(np), x_hi);
  auto see = [&](size_t ii, double x, double gx) {
    if (R_FINITE(gx) && gx > gb[ii]) { gb[ii] = gx; xb[ii] = x; }
  };
  for (int i = 0; i < np; ++i) {
    const double ci = centre_in.size() ? centre_in[i] : mu_of(i);
    c[static_cast<size_t>(i)] = clamp(R_FINITE(ci) ? ci : mu_of(i));
    // The posterior of a log-concave likelihood times the prior is no wider
    // than the prior.
    if (probe_in.size()) {
      const double hi = probe_in[probe_in.size() == 1 ? 0 : i];
      if (R_FINITE(hi) && hi > 0.0 && hi < sigma) h[static_cast<size_t>(i)] = hi;
    }
  }

  // --- 1. Gauss-Hermite rules checked by a quadratic through their own
  // nodes (least squares in d = (x - c) / h). A rule off its peak is placed
  // again at the peak that fit locates, for as long as each placement lands
  // closer to its peak than the last; the rest go on to probe triples. ---
  double S0 = 0, S2 = 0, S4 = 0;
  for (int k = 0; k < n_nodes; ++k) {
    const double d2 = z[static_cast<size_t>(k)] * z[static_cast<size_t>(k)];
    S0 += 1.0; S2 += d2; S4 += d2 * d2;
  }
  const double det = S0 * S4 - S2 * S2;
  std::vector<int> active(static_cast<size_t>(np));
  for (int i = 0; i < np; ++i) active[static_cast<size_t>(i)] = i;
  std::vector<double> last_offset(static_cast<size_t>(np), R_PosInf);
  std::vector<double> first_est(static_cast<size_t>(np), R_NegInf);
  std::vector<double> first_opt(static_cast<size_t>(np), R_PosInf);
  for (int r = 0; r < max_rounds && !active.empty(); ++r) {
    const int na = static_cast<int>(active.size());
    NumericMatrix X(na, n_nodes);
    for (int j = 0; j < na; ++j) {
      const size_t ii = static_cast<size_t>(active[static_cast<size_t>(j)]);
      for (int k = 0; k < n_nodes; ++k) X(j, k) = c[ii] + h[ii] * z[static_cast<size_t>(k)];
    }
    NumericMatrix gx = eval(active, X);
    std::vector<int> next;
    for (int j = 0; j < na; ++j) {
      const int i = active[static_cast<size_t>(j)];
      const size_t ii = static_cast<size_t>(i);
      g.rounds[i] += 1;
      double T0 = 0, T1 = 0, T2 = 0, best = R_NegInf, xbest = c[ii];
      bool finite = true;
      for (int k = 0; k < n_nodes; ++k) {
        const double gk = gx(j, k), d = z[static_cast<size_t>(k)];
        see(ii, X(j, k), gk);
        if (!R_FINITE(gk)) { finite = false; break; }
        T0 += gk; T1 += gk * d; T2 += gk * d * d;
        if (gk > best) { best = gk; xbest = X(j, k); }
      }
      // The nodes are symmetric about 0, so the odd moments vanish.
      const double b = finite ? T1 / S2 : NA_REAL;
      const double q = finite ? (S0 * T2 - S2 * T0) / det : NA_REAL;
      if (r == 0) {
        for (int k = 0; k < n_nodes; ++k) {
          g.nodes(i, k) = X(j, k);
          g.log_terms(i, k) = gx(j, k) + lwz[static_cast<size_t>(k)] -
                              R::dnorm(X(j, k), c[ii], h[ii], 1);
        }
        first_est[ii] = log_sum_exp_row(g.log_terms, i, n_nodes);
        // A misplaced rule can miss the peak, so the screen also takes the
        // integral of the fitted quadratic.
        first_opt[ii] = R_PosInf;
        if (finite && R_FINITE(q) && q < 0.0) {
          const double a0 = (S4 * T0 - S2 * T2) / det;
          const double lap = a0 - b * b / (4.0 * q) + 0.5 * std::log(2.0 * M_PI) +
                             std::log(h[ii]) - 0.5 * std::log(-2.0 * q);
          if (R_FINITE(lap)) first_opt[ii] = std::max(first_est[ii], lap);
        }
      }
      if (finite && R_FINITE(q) && q < 0.0) {
        const double vi = c[ii] - h[ii] * b / (2.0 * q);
        const double si = h[ii] / std::sqrt(-2.0 * q);
        if (R_FINITE(vi) && R_FINITE(si) && si > 0.0) {
          fitted[ii] = 1; v[ii] = vi; sv[ii] = std::min(si, sigma);
          const double offset = std::fabs(vi - c[ii]) / si;
          if (offset <= max_offset && scale_ok(h[ii], si) && rule_inside(c[ii], h[ii])) {
            done[ii] = 1;
            g.mode[i] = vi; g.sd[i] = si;
            for (int k = 0; k < n_nodes; ++k) {
              g.nodes(i, k) = X(j, k);
              g.log_terms(i, k) = gx(j, k) + lwz[static_cast<size_t>(k)] -
                                  R::dnorm(X(j, k), c[ii], h[ii], 1);
            }
            continue;
          }
          const bool closer = offset < last_offset[ii];
          last_offset[ii] = offset;
          c[ii] = clamp(vi); h[ii] = sv[ii];
          if (closer && rule_inside(c[ii], h[ii])) next.push_back(i);
          continue;
        }
      }
      c[ii] = clamp(xbest);
      h[ii] = std::min(2.0 * h[ii], sigma);
    }
    if (r == 0 && offset_in.size()) {
      // A particle whose weight, even at its optimistic value, is below the
      // best one's by more than the double-precision resolution of the batch's
      // summed weights cannot change the selection: it keeps its first rule.
      double ref = R_NegInf;
      for (int i = 0; i < np; ++i) {
        const double e = first_est[static_cast<size_t>(i)] + offset_in[i];
        if (R_FINITE(e) && e > ref) ref = e;
      }
      if (R_FINITE(ref)) {
        const double cut = ref + std::log(std::numeric_limits<double>::epsilon()) -
                           std::log(static_cast<double>(np));
        for (int i = 0; i < np; ++i) {
          const size_t ii = static_cast<size_t>(i);
          if (done[ii]) continue;
          const double o = first_opt[ii] + offset_in[i];
          if (!ISNAN(o) && o < cut) {
            done[ii] = 1;
            g.screened += 1;
            if (fitted[ii]) { g.mode[i] = v[ii]; g.sd[i] = sv[ii]; }
          }
        }
        std::vector<int> keep;
        for (int i : next) if (!done[static_cast<size_t>(i)]) keep.push_back(i);
        next.swap(keep);
      }
    }
    active.swap(next);
  }

  // --- 2. locate the remaining peaks with probe triples c +- sqrt(3) h ---
  active.clear();
  for (int i = 0; i < np; ++i) if (!done[static_cast<size_t>(i)]) active.push_back(i);
  std::vector<char> located(static_cast<size_t>(np), 0);
  const double r3 = std::sqrt(3.0);
  for (int r = 0; r < max_rounds && !active.empty(); ++r) {
    const int na = static_cast<int>(active.size());
    NumericMatrix X(na, 3);
    for (int j = 0; j < na; ++j) {
      const size_t i = static_cast<size_t>(active[static_cast<size_t>(j)]);
      const double hh = std::max(r3 * h[i], 1e-8);
      double m = clamp(c[i]);
      double a = m - hh, b = m + hh;
      if (a < x_lo) { a = x_lo; m = std::min(x_hi, a + hh); b = std::min(x_hi, m + hh); }
      if (b > x_hi) { b = x_hi; m = std::max(x_lo, b - hh); a = std::max(x_lo, m - hh); }
      X(j, 0) = a; X(j, 1) = m; X(j, 2) = b;
    }
    NumericMatrix gx = eval(active, X);
    std::vector<int> next;
    for (int j = 0; j < na; ++j) {
      const int i = active[static_cast<size_t>(j)];
      const size_t ii = static_cast<size_t>(i);
      g.rounds[i] += 1;
      for (int k = 0; k < 3; ++k) see(ii, X(j, k), gx(j, k));
      int kb = 0;
      for (int k = 1; k < 3; ++k) if (gx(j, k) > gx(j, kb)) kb = k;
      // Narrow the bracket: for a unimodal integrand the peak lies on the
      // uphill side of the triple's middle, or between its ends when the
      // middle is highest.
      if (R_FINITE(gx(j, kb))) {
        double lo = br_lo[ii], hi = br_hi[ii];
        if (kb == 0) hi = std::min(hi, X(j, 1));
        else if (kb == 2) lo = std::max(lo, X(j, 1));
        else { lo = std::max(lo, X(j, 0)); hi = std::min(hi, X(j, 2)); }
        if (!(hi > lo)) { lo = X(j, 0); hi = X(j, 2); }   // not unimodal: restart
        br_lo[ii] = lo; br_hi[ii] = hi;
      }
      // An uphill step goes at most four triple widths past the triple's
      // highest end, and never outside the bracket; once the middle is
      // highest the parabola's peak lies between the triple's ends.
      const double width = X(j, 2) - X(j, 0);
      const double far_lo = kb == 0 ? X(j, 0) - 4.0 * width : X(j, 0);
      const double far_hi = kb == 2 ? X(j, 2) + 4.0 * width : X(j, 2);
      const double h_probe = h[ii];
      double vi, si;
      if (marginal_parabola(X(j, 0), gx(j, 0), X(j, 1), gx(j, 1), X(j, 2), gx(j, 2), &vi, &si)) {
        si = std::min(si, sigma);
        const double vc = clamp(vi);
        const double dev = vc - X(j, 1);
        // A peak beyond an edge leaves the bracket flush against that edge.
        const bool pinned = (X(j, 0) <= x_lo && dev < 0.0) || (X(j, 2) >= x_hi && dev > 0.0);
        const bool centred = std::fabs(dev) <= max_offset * si;
        fitted[ii] = 1; v[ii] = vc; sv[ii] = si;
        if ((centred && scale_ok(h[ii], si)) || (pinned && scale_ok(h[ii], si))) {
          located[ii] = 1;
          continue;
        }
        // A bracket narrower than the centring tolerance locates the peak as
        // well as a centred rule would, though a kinked or truncated peak
        // never passes the curvature check.
        if (kb == 1 && br_hi[ii] - br_lo[ii] <= max_offset * si) {
          located[ii] = 1;
          v[ii] = std::min(br_hi[ii], std::max(br_lo[ii], vc));
          continue;
        }
        const bool reached = vc >= far_lo && vc <= far_hi;
        c[ii] = std::min(far_hi, std::max(far_lo, vc));
        h[ii] = reached ? si : std::min(2.0 * h[ii], sigma);
      } else if (kb == 1) {
        c[ii] = X(j, 1);
        h[ii] = 0.5 * h[ii];
      } else {
        c[ii] = kb == 0 ? far_lo : far_hi;
        h[ii] = std::min(2.0 * h[ii], sigma);
      }
      // A step out of the bracket returns to the best node evaluated, with a
      // narrower triple, so the search cannot cycle between the bracket's ends.
      const double cand = clamp(c[ii]);
      if ((cand < br_lo[ii] || cand > br_hi[ii]) && R_FINITE(xb[ii]) &&
          br_lo[ii] <= xb[ii] && xb[ii] <= br_hi[ii]) {
        c[ii] = xb[ii];
        h[ii] = 0.5 * h_probe;
      } else {
        c[ii] = std::min(br_hi[ii], std::max(br_lo[ii], cand));
      }
      next.push_back(i);
    }
    active.swap(next);
  }

  // --- 3. the final rule for every particle not settled by step 1 ---
  std::vector<int> gh_rows, gl_rows;
  std::vector<double> gl_a, gl_b;
  for (int i = 0; i < np; ++i) {
    const size_t ii = static_cast<size_t>(i);
    if (done[ii]) continue;
    if (located[ii] && rule_inside(v[ii], sv[ii])) {
      gh_rows.push_back(i);
    } else {
      // A peak near an edge, or one never resolved: Gauss-Legendre over the
      // feasible part of the window around it. An unresolved peak's last
      // fitted vertex was never checked, so the window is centred on the best
      // node evaluated and kept inside the bracket.
      const bool use_best = !located[ii] && R_FINITE(xb[ii]);
      const double cc = located[ii] ? v[ii] : (use_best ? xb[ii] : mu_of(i));
      const double ss = fitted[ii] ? sv[ii] : sigma;
      // A bracket that excludes the best node came from a multimodal
      // integrand and is not used.
      const bool bracketed = use_best && br_lo[ii] <= cc && cc <= br_hi[ii];
      const double lo = bracketed ? std::max(x_lo, br_lo[ii]) : x_lo;
      const double hi = bracketed ? std::min(x_hi, br_hi[ii]) : x_hi;
      double a = std::max(lo, cc - span * ss), b = std::min(hi, cc + span * ss);
      if (!(b > a)) {
        // The window lies wholly outside the feasible interval: integrate the
        // feasible window of the same width nearest to it, from the edge on
        // its side (finite, since the window lies beyond it).
        const double w = 2.0 * span * ss;
        if (cc > hi) { b = hi; a = std::max(lo, hi - w); }
        else { a = lo; b = std::min(hi, lo + w); }
      }
      gl_rows.push_back(i); gl_a.push_back(a); gl_b.push_back(b);
    }
  }
  if (!gh_rows.empty()) {
    const int nb = static_cast<int>(gh_rows.size());
    NumericMatrix X(nb, n_nodes);
    for (int j = 0; j < nb; ++j) {
      const size_t ii = static_cast<size_t>(gh_rows[static_cast<size_t>(j)]);
      for (int k = 0; k < n_nodes; ++k) X(j, k) = v[ii] + sv[ii] * z[static_cast<size_t>(k)];
    }
    NumericMatrix gx = eval(gh_rows, X);
    for (int j = 0; j < nb; ++j) {
      const int i = gh_rows[static_cast<size_t>(j)];
      const size_t ii = static_cast<size_t>(i);
      g.rounds[i] += 1;
      g.mode[i] = v[ii]; g.sd[i] = sv[ii];
      for (int k = 0; k < n_nodes; ++k) {
        g.nodes(i, k) = X(j, k);
        g.log_terms(i, k) = gx(j, k) + lwz[static_cast<size_t>(k)] -
                            R::dnorm(X(j, k), v[ii], sv[ii], 1);
      }
    }
  }
  if (!gl_rows.empty()) {
    const int nb = static_cast<int>(gl_rows.size());
    NumericMatrix X(nb, n_gl), LW(nb, n_gl);
    for (int j = 0; j < nb; ++j) {
      MarginalRule rule;
      marginal_append_gl_panel(rule, gl_a[static_cast<size_t>(j)], gl_b[static_cast<size_t>(j)], n_gl);
      for (int k = 0; k < n_gl; ++k) {
        X(j, k) = rule.x[static_cast<size_t>(k)];
        LW(j, k) = std::log(rule.w[static_cast<size_t>(k)]);
      }
    }
    NumericMatrix gx = eval(gl_rows, X);
    for (int j = 0; j < nb; ++j) {
      const int i = gl_rows[static_cast<size_t>(j)];
      const size_t ii = static_cast<size_t>(i);
      g.rounds[i] += 1;
      g.gl_used += 1;
      if (located[ii]) { g.mode[i] = v[ii]; g.sd[i] = sv[ii]; }
      else if (R_FINITE(xb[ii])) { g.mode[i] = xb[ii]; g.sd[i] = fitted[ii] ? sv[ii] : sigma; }
      for (int k = 0; k < n_gl; ++k) {
        g.nodes(i, k) = X(j, k);
        g.log_terms(i, k) = gx(j, k) + LW(j, k);
      }
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

// Storage-time accessor for the reconstruction (reconstruct-at-storage design,
// plan "use the posterior as sampled"). Returns the quadrature node grid and the
// unnormalized per-particle log-terms; the sampler forms w = softmax_k(log_terms)
// and draws one node per subject to write into the stored alpha, so the
// posterior carries a valid draw and predict()/make_data() stay unchanged.
// Nodes are on the sampled scale, matching alpha.
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
                      Rcpp::Named("rounds") = g.rounds,
                      Rcpp::Named("gl_used") = g.gl_used,
                      Rcpp::Named("screened") = g.screened);
}

