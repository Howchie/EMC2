// The iteration's shear step (R/iteration_steps.R .emc_step_shear,
// R/shear_move.R .emc_shear_move) when the subject likelihoods are computed
// in this process: the rounds allotted to each coordinate, runs of log-sd
// rounds on the native kernel, the step sizes and the block's records, and
// the moved state. A round along a mean direction runs through
// .emc_shear_round().
// [[Rcpp::depends(RcppArmadillo)]]
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math")
#endif
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <RcppArmadillo.h>
#include <R_ext/Lapack.h>
#include <R_ext/Random.h>
#ifndef FCONE
#define FCONE
#endif
#include <cmath>
#include <string>
#include <vector>
#include "iteration_state.h"
#include "subject_ll.h"

using namespace Rcpp;
using namespace emc_iter;

Rcpp::List emc_shear_state_cpp(Rcpp::NumericVector mu, Rcpp::NumericMatrix V,
                               const arma::cube& prec, const arma::mat& lin, Rcpp::List spec);
List shear_rounds_from(IntegerVector ks, NumericVector log_step, IntegerVector blk_n,
                       bool adapt, NumericVector step_range, double move_target,
                       int n_candidates, NumericMatrix alpha, NumericVector mu,
                       NumericMatrix V, NumericMatrix chol_V, NumericMatrix m0,
                       NumericVector R0, NumericVector logdet0, double cur_hier,
                       double cur_lp, double lp_mu, NumericVector ll,
                       const arma::cube& prec, const arma::mat& lin,
                       const LogicalVector& has, const arma::mat& U, const arma::vec& b,
                       const arma::cube& Uperp, List spec, List data, SEXP model,
                       CharacterVector par_names, SubjectLl* src);
double shear_log_hier(const NumericMatrix& alpha, const NumericVector& mu,
                      const NumericMatrix& chol_V);
List shear_round_from(const arma::vec& t, const arma::vec& v,
                      const arma::mat& alpha, const arma::vec& mu,
                      const arma::mat& V, const arma::mat& m0,
                      const arma::cube& R0, const arma::vec& logdet0,
                      const arma::cube& prec, const arma::mat& lin,
                      const LogicalVector& has, const arma::mat& U,
                      const arma::vec& b, const arma::cube& Uperp,
                      List spec, double lp_mu, double cur_hier, double cur_lp,
                      NumericVector ll_cur, List data, SEXP model,
                      CharacterVector par_names, SubjectLl* src,
                      const double* lp_mu_each, Rcpp::Nullable<Rcpp::List> curve);
Rcpp::List emc_shear_curve_fits_cpp(const arma::vec& mu, Rcpp::List curve, int n_subjects);
List shear_mean_round(const NumericVector& t, const NumericVector& d,
                      const NumericMatrix& alpha, const NumericVector& mu,
                      const NumericMatrix& chol_V, const NumericMatrix& Vi,
                      const NumericMatrix& m0, const NumericVector& R0,
                      const NumericVector& prior_mean, const NumericMatrix& prior_invar,
                      const NumericVector& ll_cur, const List& data, SEXP model,
                      const CharacterVector& par_names, SubjectLl* src);

namespace {

const int kCandidates = 8;
const double kMoveTarget = 0.68;

// The group state and the subjects as the rounds carry them (the R move's x).
struct ShearState {
  RObject mu, V, Vi, chol_V, m, R, logdet, lp_v, alpha, ll;
  double lp, cur_hier;
};

double r_sum(SEXP x) {
  long double s = 0.0L;
  for (R_xlen_t i = 0; i < XLENGTH(x); ++i) s += REAL(x)[i];
  return static_cast<double>(s);
}

// x as .emc_shear_round() takes it.
SEXP state_list(const ShearState& x, SEXP wpool, SEXP wpool_part) {
  List fm = List::create(_["m"] = x.m, _["R"] = x.R, _["logdet"] = x.logdet);
  List cur = List::create(_["mu"] = x.mu, _["V"] = x.V, _["Vi"] = x.Vi,
                          _["chol_V"] = x.chol_V, _["fm"] = fm,
                          _["logdet"] = r_sum(x.logdet), _["lp_v"] = x.lp_v, _["lp"] = x.lp);
  return List::create(_["cur"] = cur, _["cur_hier"] = x.cur_hier, _["alpha"] = x.alpha,
                      _["ll"] = x.ll, _["rounds"] = 0, _["build_el"] = 0.0, _["ll_el"] = 0.0,
                      _["wpool"] = wpool, _["wpool_part"] = wpool_part);
}

void state_from(SEXP xl, ShearState& x) {
  SEXP cur = lget(xl, "cur"), fm = lget(cur, "fm");
  x.mu = lget(cur, "mu"); x.V = lget(cur, "V"); x.Vi = lget(cur, "Vi");
  x.chol_V = lget(cur, "chol_V"); x.m = lget(fm, "m"); x.R = lget(fm, "R");
  x.logdet = lget(fm, "logdet"); x.lp_v = lget(cur, "lp_v");
  x.lp = Rf_asReal(lget(cur, "lp"));
  x.cur_hier = Rf_asReal(lget(xl, "cur_hier"));
  x.alpha = lget(xl, "alpha"); x.ll = lget(xl, "ll");
}

}  // namespace

namespace emc_iter {

// .emc_shear_draw_a_half(): the mixing scales given the covariance, drawn
// from the session's stream. False when the covariance's blocks do not
// factorise.
bool draw_a_half(SEXP V, SEXP spec, std::vector<double>& out) {
  const int p = Rf_nrows(V);
  const double* v = REAL(V);
  std::vector<double> prec(p, 0.0);
  SEXP g = lget(spec, "gamma_idx"), iw = lget(spec, "iw"), size = lget(spec, "size"),
       A = lget(spec, "A"), df = lget(spec, "v");
  if (TYPEOF(iw) != VECSXP || Rf_length(df) != 1 || Rf_length(size) != p ||
      TYPEOF(A) != REALSXP || XLENGTH(A) != p)
    return false;
  for (R_xlen_t k = 0; k < Rf_xlength(g); ++k) {
    const int at = (TYPEOF(g) == INTSXP ? INTEGER(g)[k] : static_cast<int>(REAL(g)[k])) - 1;
    prec[at] = 1 / v[at + static_cast<R_xlen_t>(at) * p];
  }
  for (R_xlen_t q = 0; q < XLENGTH(iw); ++q) {
    SEXP b = VECTOR_ELT(iw, q);
    if (TYPEOF(b) != INTSXP) return false;
    int nb = static_cast<int>(XLENGTH(b));
    std::vector<double> sub(static_cast<size_t>(nb) * nb, 0.0);
    for (int c = 0; c < nb; ++c)
      for (int r = 0; r <= c; ++r)
        sub[r + static_cast<size_t>(c) * nb] =
          v[(INTEGER(b)[r] - 1) + static_cast<R_xlen_t>(INTEGER(b)[c] - 1) * p];
    int info = 0;
    F77_CALL(dpotrf)("U", &nb, sub.data(), &nb, &info FCONE);
    if (info != 0) return false;
    F77_CALL(dpotri)("U", &nb, sub.data(), &nb, &info FCONE);
    if (info != 0) return false;
    for (int r = 0; r < nb; ++r) prec[INTEGER(b)[r] - 1] = sub[r + static_cast<size_t>(r) * nb];
  }
  const double nu = Rf_asReal(df);
  out.resize(p);
  for (int k = 0; k < p; ++k) {
    const double sz = TYPEOF(size) == INTSXP ? INTEGER(size)[k] : REAL(size)[k];
    volatile double shape = (nu + sz) / 2;
    volatile double a2 = REAL(A)[k] * REAL(A)[k];
    volatile double inv_a2 = 1 / a2;
    volatile double scaled = nu * prec[k];
    volatile double rate = scaled + inv_a2;
    volatile double scale = 1 / rate;
    out[k] = 1 / Rf_rgamma(shape, scale);
  }
  return true;
}

}  // namespace emc_iter

namespace {

void zero_costs(SEXP st) {
  set_var(st, "shear_cpu", Rf_ScalarReal(0.0));
  set_var(st, "shear_el", Rf_ScalarReal(0.0));
  set_var(st, "shear_round_cpu", Rf_ScalarReal(0.0));
  set_var(st, "shear_rounds", Rf_ScalarInteger(0));
}

}  // namespace

namespace emc_iter {

bool step_shear(SEXP st, SEXP cx, int i) {
  SEXP active = lget(cx, "shear_active");
  if (TYPEOF(active) != LGLSXP) return false;
  const std::string stage = as<std::string>(lget(cx, "stage"));
  SEXP pmwgs = var(st, "pmwgs"), pars = var(st, "pars"), proposals = var(st, "proposals");
  if (!flag(active) || stage == "preburn") {
    zero_costs(st);
    return true;
  }
  if (!flag(lget(cx, "moves_native"))) return false;
  SEXP sst = lget(pmwgs, "shear"), plan = lget(sst, "plan");
  SEXP dirs = lget(plan, "dirs"), rounds = lget(plan, "rounds"), log_step = lget(plan, "log_step");
  const int p = Rf_asInteger(lget(pmwgs, "n_pars")), n = Rf_asInteger(lget(pmwgs, "n_subjects"));
  const int q = 2 * p;
  if (TYPEOF(plan) != VECSXP || !Rf_isMatrix(dirs) || Rf_nrows(dirs) != p ||
      Rf_ncols(dirs) != p || TYPEOF(rounds) != REALSXP || XLENGTH(rounds) != q) {
    if (TYPEOF(plan) == VECSXP && TYPEOF(rounds) != REALSXP && rounds != R_NilValue) return false;
    zero_costs(st);
    return true;
  }
  bool any_rounds = false;
  for (int k = 0; k < q; ++k) if (REAL(rounds)[k] > 0) any_rounds = true;
  if (!any_rounds) {
    zero_costs(st);
    return true;
  }
  SEXP wpool = var(st, "wpool"), data = lget(pmwgs, "data"), model = lget(pmwgs, "model");
  SEXP par_names = lget(pmwgs, "par_names"), rng = lget(pmwgs, "rng"), seed = lget(rng, "gibbs");
  SEXP blk = lget(plan, "block"), prior = lget(pmwgs, "prior");
  SEXP bE = lget(blk, "E"), bn = lget(blk, "n"), bmoved = lget(blk, "moved"),
       bN = lget(blk, "N"), bK = lget(blk, "K");
  SEXP tvar = lget(pars, "tvar"), tmu = lget(pars, "tmu");
  if (!option_on("emc2.native_shear_round") ||
      flag(Rf_GetOption1(Rf_install("emc2.shear_delayed"))) ||
      TYPEOF(data) != VECSXP || XLENGTH(data) != n ||
      TYPEOF(par_names) != STRSXP || TYPEOF(seed) != INTSXP ||
      (lget(sst, "iterations") != R_NilValue &&
       (TYPEOF(lget(sst, "iterations")) != INTSXP || TYPEOF(lget(sst, "rounds")) != INTSXP)) ||
      TYPEOF(log_step) != REALSXP || XLENGTH(log_step) != q ||
      TYPEOF(bE) != REALSXP || XLENGTH(bE) != q || TYPEOF(bn) != INTSXP || XLENGTH(bn) != q ||
      TYPEOF(bmoved) != INTSXP || XLENGTH(bmoved) != q || TYPEOF(bN) != REALSXP ||
      XLENGTH(bN) != q || TYPEOF(bK) != INTSXP || !Rf_isMatrix(bK) || Rf_nrows(bK) != q ||
      Rf_getAttrib(bK, R_DimNamesSymbol) != R_NilValue || lget(blk, "iters") == R_NilValue ||
      TYPEOF(proposals) != REALSXP || !Rf_isMatrix(proposals) ||
      Rf_nrows(proposals) != p + 1 || Rf_ncols(proposals) != n ||
      TYPEOF(tvar) != REALSXP || !Rf_isMatrix(tvar) || TYPEOF(tmu) != REALSXP ||
      XLENGTH(tmu) != p)
    return false;
  // The move's counters start at zero.
  RObject sst_counted = sst;
  if (lget(sst, "iterations") == R_NilValue) {
    sst_counted = lset(sst_counted, "rounds", Rf_ScalarInteger(0));
    sst_counted = lset(sst_counted, "iterations", Rf_ScalarInteger(0));
    sst = sst_counted;
  }
  const Clock started = clock_now();
  static Function shear_context = package_function(".emc_shear_context");
  static Function r_log_prior_mu = package_function(".emc_shear_log_prior_mu");
  // The mean's log prior; R's own where its product leaves the BLAS path.
  auto prior_mu = [&](SEXP mu) {
    double v = 0.0;
    if (emc_iter::log_prior_mu(mu, prior, v)) return v;
    return Rf_asReal(r_log_prior_mu(mu, prior));
  };
  static Function shear_round = package_function(".emc_shear_round");
  static Function subject_ll = package_function(".emc_wpool_subject_ll");
  // The rounds' likelihoods: here on the stored-argument path, else through
  // the worker pool.
  SubjectLl lls;
  lls.pool = wpool;
  lls.sampler = pmwgs;
  lls.part = var(st, "wpool_part");
  lls.source = "shear";
  lls.in_process = emc_subject_ll_in_process(wpool, pmwgs);
  lls.dispatch = subject_ll;

  set_global_seed(seed);
  GetRNGstate();
  bool in_rng = true;
  // Leaving through an R call or early: R's stream is the one the draws left.
  auto leave_rng = [&]() { if (in_rng) { PutRNGstate(); in_rng = false; } };
  auto enter_rng = [&]() { if (!in_rng) { GetRNGstate(); in_rng = true; } };

  std::vector<int> n_rounds(q);
  int total = 0;
  for (int k = 0; k < q; ++k) {
    const double r = REAL(rounds)[k], whole = std::floor(r);
    n_rounds[k] = static_cast<int>(whole) + (unif_rand() < r - whole ? 1 : 0);
    total += n_rounds[k];
  }
  RObject sst2 = lset(sst, "iterations",
                      Rf_ScalarInteger(INTEGER(lget(sst, "iterations"))[0] + 1));
  // The move ends without rounds: the stream and the iteration count advance.
  auto finish_idle = [&]() {
    leave_rng();
    Shield<SEXP> seed_after(global_seed());
    Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
    RObject sampler = lset(pmwgs, "rng", rng2);
    sampler = lset(sampler, "shear", sst2);
    set_var(st, "pmwgs", sampler);
    zero_costs(st);
    const double el = clock_now().elapsed - started.elapsed;
    set_var(st, "shear_el", Rf_ScalarReal(el));
    if (flag(lget(cx, "alloc_on")))
      set_var(st, "shear_cpu", Rf_ScalarReal(reserved_cost(el, num(var(st, "cores_it"), 1.0))));
    return true;
  };
  if (total == 0) return finish_idle();

  leave_rng();
  SEXP held = emc_iter::move_context("shear", pmwgs);
  Shield<SEXP> ctx(held != R_NilValue ? held : static_cast<SEXP>(shear_context(pmwgs)));
  SEXP args = lget(ctx, "args"), spec = lget(ctx, "spec");
  if (TYPEOF(args) != VECSXP || TYPEOF(spec) != VECSXP) return false;
  arma::cube cube = as<arma::cube>(lget(args, "cube"));
  arma::mat lin = as<arma::mat>(lget(args, "lin"));
  // Curved frames: the fits follow the mean along the curve's coordinate, so
  // a run of log-sd rounds takes the curve's fits at the mean it starts
  // from (stale: a mean round moved since they were last taken).
  SEXP packed = lget(ctx, "packed");
  SEXP curve = TYPEOF(packed) == VECSXP ? lget(packed, "curve") : R_NilValue;
  const bool curved = TYPEOF(curve) == VECSXP;
  auto fits_at = [&](SEXP mu) {
    const List f = emc_shear_curve_fits_cpp(as<arma::vec>(mu), List(curve), n);
    cube = as<arma::cube>(f["cube"]);
    lin = as<arma::mat>(f["lin"]);
  };
  std::vector<char> on_curve(p, 0);
  if (curved) {
    SEXP row = lget(curve, "row");
    if (TYPEOF(row) != REALSXP || XLENGTH(row) != p) return false;
    double rr = 0.0;
    for (int r = 0; r < p; ++r) rr += REAL(row)[r] * REAL(row)[r];
    for (int k = 0; k < p; ++k) {
      double rd = 0.0, dd = 0.0;
      for (int r = 0; r < p; ++r) {
        const double d = REAL(dirs)[r + static_cast<R_xlen_t>(k) * p];
        rd += REAL(row)[r] * d;
        dd += d * d;
      }
      on_curve[k] = std::fabs(rd) > 1e-9 * std::sqrt(rr * dd);
    }
    fits_at(tmu);
  }
  bool stale = false;
  const arma::mat U = as<arma::mat>(lget(args, "U"));
  const arma::vec wall_b = as<arma::vec>(lget(args, "b"));
  const arma::cube Uc = as<arma::cube>(lget(args, "Uc"));
  const LogicalVector has(lget(args, "has"));

  ShearState x;
  {
    NumericVector mu(p);
    std::copy(REAL(tmu), REAL(tmu) + p, mu.begin());
    const List out = emc_shear_state_cpp(mu, NumericMatrix(tvar), cube, lin, List(spec));
    if (!flag(out["ok"])) return finish_idle();
    x.mu = mu; x.V = tvar; x.Vi = out["Vi"]; x.chol_V = out["chol_V"];
    x.m = out["m"]; x.R = out["R"]; x.logdet = out["logdet"]; x.lp_v = out["lp_v"];
    x.lp = prior_mu(mu) + Rf_asReal(out["lp_v"]);
    if (!R_FINITE(x.lp)) return finish_idle();
  }
  {
    NumericMatrix a(p, n);
    NumericVector ll(n);
    const double* pr = REAL(proposals);
    for (int s = 0; s < n; ++s) {
      for (int k = 0; k < p; ++k) a(k, s) = pr[k + static_cast<R_xlen_t>(s) * (p + 1)];
      ll[s] = pr[p + static_cast<R_xlen_t>(s) * (p + 1)];
    }
    x.alpha = a;
    x.ll = ll;
    x.cur_hier = shear_log_hier(a, NumericVector(static_cast<SEXP>(x.mu)),
                                NumericMatrix(static_cast<SEXP>(x.chol_V)));
  }
  const bool adapt = stage != "sample";
  NumericVector steps = clone(NumericVector(log_step));
  NumericVector E = clone(NumericVector(bE)), N = clone(NumericVector(bN));
  IntegerVector cnt = clone(IntegerVector(bn)), mv = clone(IntegerVector(bmoved));
  std::vector<double> net(q, 0.0);
  SEXP range = lget(cx, "shear_step_range");
  if (TYPEOF(range) != REALSXP || XLENGTH(range) != 2) return false;
  const double lo = REAL(range)[0], hi = REAL(range)[1];

  // The rounds in random order.
  enter_rng();
  std::vector<int> order(total);
  {
    std::vector<int> take, pool(total);
    take.reserve(total);
    for (int k = 0; k < q; ++k) for (int r = 0; r < n_rounds[k]; ++r) take.push_back(k + 1);
    for (int k = 0; k < total; ++k) pool[k] = k;
    int left = total;
    for (int k = 0; k < total; ++k) {
      const int pick = static_cast<int>(R_unif_index(left));
      order[k] = take[pool[pick]];
      pool[pick] = pool[--left];
    }
  }
  int evaluated = 0;
  bool any_moved = false;
  const Clock rounds_started = clock_now();
  // One round through .emc_shear_round(), with its steps when already drawn.
  auto r_round = [&](int k, SEXP t) {
    leave_rng();
    NumericVector v(q);
    if (k <= p) {
      for (int r = 0; r < p; ++r) v[r] = REAL(dirs)[r + static_cast<R_xlen_t>(k - 1) * p];
    } else {
      v[k - 1] = 1.0;
    }
    Shield<SEXP> xl(state_list(x, lls.pool, lls.part));
    Shield<SEXP> out(shear_round(xl, v, std::exp(steps[k - 1]), pmwgs, ctx, t));
    enter_rng();
    state_from(out, x);
    lls.pool = lget(out, "wpool");
    evaluated += Rf_asInteger(lget(out, "rounds"));
    const bool moved = flag(lget(out, "moved"));
    if (moved && k <= p) stale = true;
    net[k - 1] += Rf_asReal(lget(out, "t"));
    E[k - 1] = E[k - 1] + Rf_asReal(lget(out, "E"));
    cnt[k - 1] = cnt[k - 1] + 1;
    mv[k - 1] = mv[k - 1] + (moved ? 1 : 0);
    if (adapt) {
      volatile double gain = ((moved ? 1 : 0) - kMoveTarget) /
                             std::sqrt(static_cast<double>(cnt[k - 1]));
      steps[k - 1] = std::min(hi, std::max(lo, steps[k - 1] + gain));
    }
    any_moved = any_moved || moved;
  };
  // A mean-direction round on the line of its candidates; through
  // .emc_shear_round(), with the steps drawn here, when a likelihood needs
  // the full path.
  SEXP prior_mean = lget(prior, "theta_mu_mean"), prior_invar = lget(prior, "theta_mu_invar");
  const bool mean_line = lget(ctx, "walls") == R_NilValue && TYPEOF(prior_mean) == REALSXP &&
    XLENGTH(prior_mean) == p && TYPEOF(prior_invar) == REALSXP && Rf_isMatrix(prior_invar) &&
    Rf_nrows(prior_invar) == p && Rf_ncols(prior_invar) == p;
  auto mean_round = [&](int k) {
    const double step = std::exp(steps[k - 1]);
    NumericVector t(kCandidates), d(p);
    volatile double t_aux = step * norm_rand();
    for (int c = 0; c < kCandidates; ++c) {
      volatile double inc = step * norm_rand();
      t[c] = t_aux + inc;
    }
    for (int r = 0; r < p; ++r) d[r] = REAL(dirs)[r + static_cast<R_xlen_t>(k - 1) * p];
    const List o = shear_mean_round(
        t, d, NumericMatrix(static_cast<SEXP>(x.alpha)), NumericVector(static_cast<SEXP>(x.mu)),
        NumericMatrix(static_cast<SEXP>(x.chol_V)), NumericMatrix(static_cast<SEXP>(x.Vi)),
        NumericMatrix(static_cast<SEXP>(x.m)), NumericVector(static_cast<SEXP>(x.R)),
        NumericVector(prior_mean), NumericMatrix(prior_invar),
        NumericVector(static_cast<SEXP>(x.ll)), List(data), model, CharacterVector(par_names),
        &lls);
    if (as<bool>(o["fallback"])) {
      r_round(k, t);
      return;
    }
    if (as<bool>(o["evaluated"])) ++evaluated;
    const int j = as<int>(o["j"]);
    const bool moved = j > 0;
    if (moved) {
      stale = true;
      x.mu = o["mu"]; x.m = o["m"]; x.alpha = o["alpha"]; x.ll = o["ll"];
      leave_rng();
      x.lp = prior_mu(x.mu) + Rf_asReal(x.lp_v);
      enter_rng();
      x.cur_hier = shear_log_hier(NumericMatrix(static_cast<SEXP>(x.alpha)),
                                  NumericVector(static_cast<SEXP>(x.mu)),
                                  NumericMatrix(static_cast<SEXP>(x.chol_V)));
      net[k - 1] += t[j - 1];
    }
    E[k - 1] = E[k - 1] + as<double>(o["E"]);
    cnt[k - 1] = cnt[k - 1] + 1;
    mv[k - 1] = mv[k - 1] + (moved ? 1 : 0);
    if (adapt) {
      volatile double gain = ((moved ? 1 : 0) - kMoveTarget) /
                             std::sqrt(static_cast<double>(cnt[k - 1]));
      steps[k - 1] = std::min(hi, std::max(lo, steps[k - 1] + gain));
    }
    any_moved = any_moved || moved;
  };
  // A mean-direction round with transport walls, or along a curve's
  // coordinate: the general candidates, each with the prior density of its
  // own mean.
  auto walled_round = [&](int k) {
    const double step = std::exp(steps[k - 1]);
    NumericVector t(kCandidates), lp_each(kCandidates);
    volatile double t_aux = step * norm_rand();
    for (int c = 0; c < kCandidates; ++c) {
      volatile double inc = step * norm_rand();
      t[c] = t_aux + inc;
    }
    arma::vec v(q, arma::fill::zeros);
    for (int r = 0; r < p; ++r) v[r] = REAL(dirs)[r + static_cast<R_xlen_t>(k - 1) * p];
    bool native_prior = true;
    for (int c = 0; c < kCandidates && native_prior; ++c) {
      NumericVector mu_c(p);
      for (int r = 0; r < p; ++r) {
        volatile double shift = t[c] * v[r];
        mu_c[r] = REAL(x.mu)[r] + shift;
      }
      native_prior = emc_iter::log_prior_mu(mu_c, prior, lp_each[c]);
    }
    if (!native_prior) {
      r_round(k, t);
      return;
    }
    const List o = shear_round_from(
        as<arma::vec>(t), v, as<arma::mat>(x.alpha), as<arma::vec>(x.mu), as<arma::mat>(x.V),
        as<arma::mat>(x.m), as<arma::cube>(x.R), as<arma::vec>(x.logdet), cube, lin, has, U,
        wall_b, Uc, List(spec), 0.0, x.cur_hier, x.lp, NumericVector(static_cast<SEXP>(x.ll)),
        List(data), model, CharacterVector(par_names), &lls, lp_each.begin(),
        curved ? Rcpp::Nullable<List>(curve) : Rcpp::Nullable<List>(R_NilValue));
    if (as<bool>(o["evaluated"])) ++evaluated;
    const int j = as<int>(o["j"]);
    const bool moved = j > 0;
    if (moved) stale = true;
    if (moved) {
      SEXP V2 = o["V"];
      Shield<SEXP> Vn(Rf_duplicate(V2));
      Rf_setAttrib(Vn, R_DimNamesSymbol, Rf_getAttrib(x.V, R_DimNamesSymbol));
      x.mu = o["mu"]; x.V = Vn; x.Vi = o["Vi"]; x.chol_V = o["chol_V"];
      x.m = o["m"]; x.R = o["R"]; x.logdet = o["logdet"]; x.lp_v = o["lp_v"];
      x.lp = as<double>(o["lp"]);
      x.alpha = o["alpha"]; x.ll = o["ll"];
      x.cur_hier = shear_log_hier(NumericMatrix(static_cast<SEXP>(x.alpha)),
                                  NumericVector(static_cast<SEXP>(x.mu)),
                                  NumericMatrix(static_cast<SEXP>(x.chol_V)));
      net[k - 1] += t[j - 1];
    }
    E[k - 1] = E[k - 1] + as<double>(o["E"]);
    cnt[k - 1] = cnt[k - 1] + 1;
    mv[k - 1] = mv[k - 1] + (moved ? 1 : 0);
    if (adapt) {
      volatile double gain = ((moved ? 1 : 0) - kMoveTarget) /
                             std::sqrt(static_cast<double>(cnt[k - 1]));
      steps[k - 1] = std::min(hi, std::max(lo, steps[k - 1] + gain));
    }
    any_moved = any_moved || moved;
  };
  int at = 0;
  while (at < total) {
    const int k = order[at];
    if (k <= p) {
      if (curved && on_curve[k - 1]) walled_round(k);
      else if (mean_line) mean_round(k);
      else if (lget(ctx, "walls") != R_NilValue) walled_round(k);
      else r_round(k, R_NilValue);
      ++at;
      continue;
    }
    int last = at;
    while (last + 1 < total && order[last + 1] > p) ++last;
    IntegerVector ks(last - at + 1);
    for (int r = at; r <= last; ++r) ks[r - at] = order[r];
    leave_rng();
    const double lp_mu = prior_mu(x.mu);
    enter_rng();
    if (curved && stale) {
      fits_at(x.mu);
      stale = false;
    }
    const List o = shear_rounds_from(
        ks, steps, cnt, adapt, NumericVector(range), kMoveTarget, kCandidates,
        NumericMatrix(static_cast<SEXP>(x.alpha)), NumericVector(static_cast<SEXP>(x.mu)),
        NumericMatrix(static_cast<SEXP>(x.V)), NumericMatrix(static_cast<SEXP>(x.chol_V)),
        NumericMatrix(static_cast<SEXP>(x.m)), NumericVector(static_cast<SEXP>(x.R)),
        NumericVector(static_cast<SEXP>(x.logdet)), x.cur_hier, x.lp, lp_mu,
        NumericVector(static_cast<SEXP>(x.ll)), cube, lin, has, U, wall_b, Uc, List(spec),
        List(data), model, CharacterVector(par_names), &lls);
    evaluated += as<int>(o["evaluated"]);
    if (as<bool>(o["any_moved"])) {
      SEXP V2 = o["V"];
      Shield<SEXP> Vn(Rf_duplicate(V2));
      Rf_setAttrib(Vn, R_DimNamesSymbol, Rf_getAttrib(x.V, R_DimNamesSymbol));
      x.mu = o["mu"]; x.V = Vn; x.Vi = o["Vi"]; x.chol_V = o["chol_V"];
      x.m = o["m"]; x.R = o["R"]; x.logdet = o["logdet"]; x.lp_v = o["lp_v"];
      x.lp = as<double>(o["lp"]);
      x.alpha = o["alpha"]; x.ll = o["ll"];
      x.cur_hier = as<double>(o["cur_hier"]);
      any_moved = true;
    }
    const int done = as<int>(o["done"]);
    const NumericVector o_t = o["t"], o_E = o["E"], o_steps = o["log_step"];
    const LogicalVector o_moved = o["moved"];
    for (int r = 0; r < done; ++r) {
      const int kr = ks[r] - 1;
      net[kr] += o_t[r];
      E[kr] = E[kr] + o_E[r];
      cnt[kr] = cnt[kr] + 1;
      mv[kr] = mv[kr] + (o_moved[r] ? 1 : 0);
    }
    steps = clone(o_steps);
    at += done;
    SEXP t_next = o["t_next"];
    if (Rf_length(t_next) > 0) {
      r_round(order[at], t_next);
      ++at;
    } else if (done == 0) {
      leave_rng();
      return false;
    }
  }
  const double round_el = clock_now().elapsed - rounds_started.elapsed;

  // Each coordinate's net move over its rounds, by how many it had.
  int max_k = 0;
  for (int k = 0; k < q; ++k) {
    if (n_rounds[k] <= 0) continue;
    volatile double sq = net[k] * net[k];
    N[k] = N[k] + sq;
    max_k = std::max(max_k, n_rounds[k]);
  }
  const int k_cols = std::max(Rf_ncols(bK), max_k);
  IntegerMatrix K(q, k_cols);
  std::copy(INTEGER(bK), INTEGER(bK) + XLENGTH(bK), K.begin());
  for (int k = 0; k < q; ++k) if (n_rounds[k] > 0) K(k, n_rounds[k] - 1) += 1;
  RObject blk2 = lset(blk, "E", E);
  blk2 = lset(blk2, "n", cnt);
  blk2 = lset(blk2, "moved", mv);
  blk2 = lset(blk2, "N", N);
  blk2 = lset(blk2, "K", K);
  Rf_copyMostAttrib(log_step, steps);
  RObject plan2 = lset(plan, "log_step", steps);
  plan2 = lset(plan2, "block", blk2);
  sst2 = lset(sst2, "plan", plan2);
  sst2 = lset(sst2, "rounds", Rf_ScalarInteger(INTEGER(lget(sst, "rounds"))[0] + evaluated));

  RObject pars2 = pars;
  if (flag(lget(spec, "uses_a_half"))) {
    std::vector<double> ah;
    if (!draw_a_half(x.V, spec, ah)) {
      leave_rng();
      return false;
    }
    NumericVector a_half(ah.begin(), ah.end());
    SEXP nm = Rf_getAttrib(lget(pars, "a_half"), R_NamesSymbol);
    if (nm != R_NilValue) Rf_setAttrib(a_half, R_NamesSymbol, nm);
    pars2 = lset(pars2, "a_half", a_half);
  }
  leave_rng();
  Shield<SEXP> seed_after(global_seed());
  RObject sampler = pmwgs;
  {
    Shield<SEXP> rng2(lset(rng, "gibbs", seed_after));
    sampler = lset(sampler, "rng", rng2);
  }
  RObject props2 = proposals;
  if (any_moved) {
    Shield<SEXP> mu2(Rf_duplicate(x.mu));
    Rf_setAttrib(mu2, R_NamesSymbol, Rf_getAttrib(tmu, R_NamesSymbol));
    Shield<SEXP> V2(Rf_duplicate(x.V)), Vi2(Rf_duplicate(x.Vi));
    Rf_setAttrib(V2, R_DimNamesSymbol, Rf_getAttrib(tvar, R_DimNamesSymbol));
    Rf_setAttrib(Vi2, R_DimNamesSymbol, Rf_getAttrib(tvar, R_DimNamesSymbol));
    pars2 = lset(pars2, "tmu", mu2);
    pars2 = lset(pars2, "tvar", V2);
    pars2 = lset(pars2, "tvinv", Vi2);
    SEXP subj_mu = lget(pars, "subj_mu");
    if (Rf_isMatrix(subj_mu)) {
      if (TYPEOF(subj_mu) != REALSXP) return false;
      Shield<SEXP> sm(Rf_duplicate(subj_mu));
      for (R_xlen_t k = 0; k < XLENGTH(sm); ++k) REAL(sm)[k] = REAL(x.mu)[k % p];
      pars2 = lset(pars2, "subj_mu", sm);
    }
    props2 = Rf_duplicate(proposals);
    double* pr = REAL(props2);
    for (int s = 0; s < n; ++s) {
      for (int k = 0; k < p; ++k)
        pr[k + static_cast<R_xlen_t>(s) * (p + 1)] = REAL(x.alpha)[k + static_cast<R_xlen_t>(s) * p];
      pr[p + static_cast<R_xlen_t>(s) * (p + 1)] = REAL(x.ll)[s];
    }
    pars2 = lset(pars2, "alpha", x.alpha);
  }
  sampler = lset(sampler, "shear", sst2);
  set_var(st, "pmwgs", sampler);
  set_var(st, "pars", pars2);
  set_var(st, "proposals", props2);
  set_var(st, "wpool", lls.pool);
  const double el = clock_now().elapsed - started.elapsed;
  const double cores = num(var(st, "cores_it"), 1.0);
  const bool alloc_on = flag(lget(cx, "alloc_on"));
  set_var(st, "shear_rounds", Rf_ScalarInteger(evaluated));
  set_var(st, "shear_el", Rf_ScalarReal(el));
  set_var(st, "shear_cpu", Rf_ScalarReal(alloc_on ? reserved_cost(el, cores) : 0.0));
  set_var(st, "shear_round_cpu", Rf_ScalarReal(alloc_on ? reserved_cost(round_el, cores) : 0.0));
  return true;
}

}  // namespace emc_iter
