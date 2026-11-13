// [[Rcpp::depends(RcppArmadillo)]]
#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <RcppArmadillo.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <ctime>
#include <Rmath.h>
#include <cmath>
#include <algorithm>
#include <climits>
#include <limits>
#include <string>
#include <vector>

using namespace Rcpp;

arma::vec fast_dmvnorm_rooti(const arma::mat& x, const arma::rowvec& mean,
                             const arma::mat& rooti, double log_const);

NumericVector calc_ll_oo(NumericMatrix particle_matrix, DataFrame data, NumericVector constants,
                         List designs, String type, List bounds, List transforms, List pretransforms,
                         CharacterVector p_types, double min_ll, Rcpp::Nullable<Rcpp::List> trend,
                         Rcpp::Nullable<Rcpp::List> marginalise,
                         Rcpp::Nullable<Rcpp::LogicalVector> varying,
                         SEXP pt_mapper_context);

SEXP emc2_ll_cached(SEXP proposals, SEXP dadm, SEXP model, SEXP varying);
bool emc2_ll_cache_has(SEXP dadm, SEXP model, SEXP cols);

namespace {

// The block runs R's sampling arithmetic in R's order, so every draw, sum and
// selection matches new_particle() on the same random stream. A step that R
// would reject (non-finite weights, a likelihood error) ends the block at the
// start of its iteration; R then takes over from the returned stream.
struct BlockAbort {};

// `%*%` for finite operands (R's default matprod): dgemv when one side is a
// vector, dgemm otherwise.
void r_matprod(const double* x, int nrx, int ncx, const double* y, int nry,
               int ncy, double* z) {
  const char trans_n = 'N', trans_t = 'T';
  const double one = 1.0, zero = 0.0;
  const arma::blas_int ione = 1;
  arma::blas_int nrx_b = nrx, ncx_b = ncx, nry_b = nry, ncy_b = ncy;
  if (ncy == 1) {
    arma::blas::gemv<double>(&trans_n, &nrx_b, &ncx_b, &one, x, &nrx_b, y,
                             &ione, &zero, z, &ione);
  } else if (nrx == 1) {
    arma::blas::gemv<double>(&trans_t, &nry_b, &ncy_b, &one, y, &nry_b, x,
                             &ione, &zero, z, &ione);
  } else {
    arma::blas::gemm<double>(&trans_n, &trans_n, &nrx_b, &ncy_b, &ncx_b, &one,
                             x, &nrx_b, y, &nry_b, &zero, z, &nrx_b);
  }
}

// particle_draws(n, mu, R = R): rnorm(n * p) filled by column, times R, plus
// mu by column. Writes rows [row0, row0 + n) of the column-major `out`.
void draw_rows(int n, int p, const double* mu, const arma::mat& R,
               std::vector<double>& z, std::vector<double>& x, double* out,
               int out_rows, int row0) {
  z.resize(static_cast<size_t>(n) * p);
  x.resize(static_cast<size_t>(n) * p);
  for (size_t t = 0; t < z.size(); ++t) z[t] = norm_rand();
  r_matprod(z.data(), n, p, R.memptr(), p, p, x.data());
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i)
      out[(row0 + i) + static_cast<size_t>(j) * out_rows] = x[i + static_cast<size_t>(j) * n] + mu[j];
}

// sample(n, 1, prob = w): FixupProb then ProbSampleNoReplace. Returns the
// 0-based row.
int r_sample_one(std::vector<double> p) {
  const int n = static_cast<int>(p.size());
  double sum = 0.0;
  int npos = 0;
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(p[i]) || p[i] < 0.0) throw BlockAbort();
    if (p[i] > 0.0) { ++npos; sum += p[i]; }
  }
  if (npos == 0) throw BlockAbort();
  for (int i = 0; i < n; ++i) p[i] /= sum;
  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = i + 1;
  revsort(p.data(), perm.data(), n);
  const double rT = 1.0 * unif_rand();
  double mass = 0.0;
  int j = 0;
  for (j = 0; j < n - 1; ++j) {
    mass += p[j];
    if (rT <= mass) break;
  }
  return perm[j] - 1;
}

double r_max(const std::vector<double>& v) {
  // max(): NaN wins, then the largest value.
  double out = R_NegInf;
  bool nan = false;
  for (double value : v) {
    if (ISNAN(value)) nan = true;
    else if (value > out) out = value;
  }
  return nan ? R_NaN : out;
}

// R's sum() of doubles: a long double accumulator.
double r_sum(const std::vector<double>& v) {
  long double s = 0.0L;
  for (double value : v) s += value;
  return static_cast<double>(s);
}

// The adapted proposal settings (update_pm_settings), as R stores them.
struct PmState {
  double iter = 1.0;
  std::vector<double> proposal_counts, acc_counts, jump_counts, epsilon, mix;
  std::vector<std::string> step_layout;
  bool has_layout = false;
  double n_particles = 0.0;
  double weight_ess = 0.0;
  bool weight_set = false, layout_set = false, particles_set = false;
};

struct TuneArgs {
  double n0, step_target, mix_adapt, target_ESS, ESS_scale, max_particles;
};

// update_epsilon_continuous() for one component.
double update_epsilon(double epsilon, double acceptance, double target,
                      double iter, double d, double alpha_star, double damp,
                      double lo, double hi) {
  const double log_eps = std::log(epsilon);
  const double c_term =
    (1.0 - 1.0 / d) * std::sqrt(2.0 * M_PI) *
      std::exp(alpha_star * alpha_star / 2.0) / (2.0 * alpha_star) +
    1.0 / (d * target * (1.0 - target));
  const double step_size = c_term / std::max(damp, iter);
  const double diff_accept = acceptance - target;
  const double eps_new = std::exp(log_eps + step_size * diff_accept);
  return std::min(hi, std::max(eps_new, lo));
}

// update_pm_settings() with the layout new_particle() passes.
void update_pm(PmState& pm, const std::vector<double>& better,
               const std::vector<int>& numbers, double weight_ess,
               const std::vector<double>& jump, const TuneArgs& tune, int n_pars,
               const std::vector<std::string>& layout) {
  const int K = static_cast<int>(numbers.size());
  pm.iter = pm.iter + 1.0;
  if (pm.iter > tune.n0) {
    for (int k = 0; k < K; ++k) {
      pm.proposal_counts[k] = pm.proposal_counts[k] + numbers[k];
      pm.acc_counts[k] = pm.acc_counts[k] + better[k];
    }
    if (pm.jump_counts.size() != jump.size())
      pm.jump_counts.assign(jump.size(), 0.0);
    for (size_t k = 0; k < jump.size(); ++k)
      pm.jump_counts[k] = pm.jump_counts[k] + jump[k];

    // update_local_steps(better[-1], numbers[-1], ...)
    const int n_comp = K - 1;
    std::vector<double> eps(std::max<size_t>(n_comp, pm.epsilon.size()));
    for (size_t t = 0; t < eps.size(); ++t)
      eps[t] = pm.epsilon[t % pm.epsilon.size()];
    const double rw_scale = 2.38 / std::sqrt(static_cast<double>(n_pars));
    const bool old_ok = pm.has_layout &&
      static_cast<int>(pm.step_layout.size()) == n_comp;
    const double target = tune.step_target;
    const double alpha_star = -R::qnorm(target / 2.0, 0.0, 1.0, 1, 0);
    for (int c = 0; c < n_comp; ++c) {
      const bool local = !layout[c].empty();
      const bool fresh = local && (!old_ok || pm.step_layout[c] != layout[c]);
      if (fresh) {
        eps[c] = rw_scale;
      } else if (local && numbers[c + 1] > 0) {
        eps[c] = update_epsilon(eps[c],
                                better[c + 1] / static_cast<double>(numbers[c + 1]),
                                target, pm.iter, static_cast<double>(n_pars),
                                alpha_star, 100.0, rw_scale / 4.0, 5.0);
      }
    }
    pm.step_layout = layout;
    pm.has_layout = true;
    pm.layout_set = true;
    pm.epsilon = eps;

    if (pm.mix.size() > 2 && pm.jump_counts.size() == pm.mix.size()) {
      const size_t m = pm.mix.size();
      std::vector<double> per(m);
      for (size_t k = 0; k < m; ++k)
        per[k] = pm.proposal_counts[k] > 0
          ? pm.jump_counts[k] / pm.proposal_counts[k] : 0.0;
      const double per_sum = r_sum(per);
      if (per_sum > 0) {
        std::vector<double> nm(m);
        for (size_t k = 0; k < m; ++k) {
          const double perf = per[k] / per_sum;
          nm[k] = std::max((1.0 - tune.mix_adapt) * pm.mix[k] +
                             tune.mix_adapt * perf, 0.02);
        }
        const double nm_sum = r_sum(nm);
        for (size_t k = 0; k < m; ++k) pm.mix[k] = nm[k] / nm_sum;
      }
    }
  }
  pm.weight_ess = weight_ess;
  pm.weight_set = true;
}

// One subject's transition: the proposals' fixed parts.
struct Plan {
  int K = 0, p = 0;
  std::vector<arma::rowvec> mu;          // centres (unused where state-centred)
  std::vector<arma::mat> bR, brooti;     // unscaled factors per proposal
  std::vector<double> bsld;
  double log_const0 = 0.0, multiplier = 1.0;
  std::vector<std::vector<int>> step_ks;
  std::vector<bool> state_centred;
  arma::rowvec group_mu;
  std::vector<std::string> layout;
  bool have_jump = false;
  arma::mat jump_rooti;
};

struct LlArgs {
  DataFrame data;
  NumericVector constants;
  List designs;
  String type;
  List bounds, transforms, pretransforms;
  CharacterVector p_types;
  double min_ll;
  Nullable<List> trend;
  LogicalVector varying;
  SEXP pt_context;
  CharacterVector par_names;
  explicit LlArgs(const List& a)
    : data(a["data"]), constants(a["constants"]), designs(a["designs"]),
      type(as<String>(a["type"])), bounds(a["bounds"]),
      transforms(a["transforms"]), pretransforms(a["pretransforms"]),
      p_types(a["p_types"]), min_ll(as<double>(a["min_ll"])),
      trend(a["trend"]), varying(a["varying"]), pt_context(a["pt_context"]),
      par_names(a["par_names"]) {}
};

struct SelState {
  std::vector<double> n, disp;
  double steps = 0.0;
};

// Runs up to n_iter iterations from (state, state_ll) on the current stream
// (the caller holds GetRNGstate). Returns the number completed; on an abort
// `seed_at_abort` is the stream at the start of the aborted iteration and
// every setting is as that iteration found it.
int run_iterations(const Plan& plan, const LlArgs& ll, int n_iter,
                   std::vector<double>& state, double& state_ll, PmState& pm,
                   const TuneArgs& tune, bool gd_good, bool adapt,
                   SelState& sel, double* alpha_out, double* ll_out,
                   RObject& seed_at_abort) {
  const int K = plan.K, p = plan.p;
  const double p_dbl = static_cast<double>(p);
  std::vector<int> particle_numbers(K);
  std::vector<arma::mat> R(K), rooti(K);
  std::vector<double> log_consts(K);
  std::vector<double> z, x;
  int done = 0;
  for (int it = 0; it < n_iter; ++it) {
    // The stream at the start of this iteration, for R to resume from.
    PutRNGstate();
    RObject seed = Environment::global_env().get(".Random.seed");
    const SelState sel_before = sel;
    const PmState pm_before = pm;
    try {
      // The proposals at this iteration's step sizes: independence
      // components keep scale one.
      for (int k = 0; k < K; ++k) {
        // rep_len(c(1, epsilon), K), then 1 for independence components.
        const size_t e = k % (pm.epsilon.size() + 1);
        const double eps_k = (e == 0 || !plan.state_centred[k]) ? 1.0
          : pm.epsilon[e - 1];
        R[k] = plan.bR[k] * eps_k;
        rooti[k] = plan.brooti[k] / eps_k;
        log_consts[k] = (plan.bsld[k] - p_dbl * std::log(eps_k)) + plan.log_const0;
      }
      // rmultinom(1, n_particles * multiplier, mix) after FixupProb.
      const double size_d = pm.n_particles * plan.multiplier;
      if (!R_FINITE(size_d) || size_d < 0 || size_d >= INT_MAX) throw BlockAbort();
      const int size = static_cast<int>(size_d);
      std::vector<double> prob = pm.mix;
      {
        double sum = 0.0;
        int npos = 0;
        for (double v : prob) {
          if (!R_FINITE(v) || v < 0.0) throw BlockAbort();
          if (v > 0.0) { ++npos; sum += v; }
        }
        if (npos == 0) throw BlockAbort();
        for (double& v : prob) v /= sum;
      }
      rmultinom(size, prob.data(), K, particle_numbers.data());
      std::vector<double> better(K, 0.0), jump(K, 0.0);
      double weight_ess = 0.0;
      std::vector<double> ref = state;
      double ref_ll = state_ll;
      bool have_step_ref = false;
      double step_ref_ll = 0.0;
      for (const std::vector<int>& ks : plan.step_ks) {
        int n_part = 0;
        for (int k : ks) n_part += particle_numbers[k];
        if (n_part == 0) continue;
        std::vector<arma::rowvec> mus_i(K);
        for (int k : ks) {
          if (plan.state_centred[k]) {
            // Auxiliary centre z ~ N(ref, S_k).
            std::vector<double> centre(p);
            draw_rows(1, p, ref.data(), R[k], z, x, centre.data(), 1, 0);
            mus_i[k] = arma::rowvec(centre);
          } else {
            mus_i[k] = plan.mu[k];
          }
        }
        const int n_rows = 1 + n_part;
        arma::mat P(n_rows, p);
        for (int j = 0; j < p; ++j) P(0, j) = ref[j];
        int next_row = 1;
        for (int k : ks) {
          const int n_k = particle_numbers[k];
          if (n_k > 0)
            draw_rows(n_k, p, mus_i[k].memptr(), R[k], z, x, P.memptr(), n_rows, next_row);
          next_row += n_k;
        }

        // Likelihoods: after the first step the reference row is the
        // previous step's winner, already evaluated.
        const int first = have_step_ref ? 1 : 0;
        NumericMatrix pmat(n_rows - first, p);
        for (int j = 0; j < p; ++j)
          for (int i = first; i < n_rows; ++i) pmat(i - first, j) = P(i, j);
        colnames(pmat) = ll.par_names;
        NumericVector ll_eval;
        try {
          ll_eval = calc_ll_oo(pmat, ll.data, ll.constants, ll.designs, ll.type,
                               ll.bounds, ll.transforms, ll.pretransforms,
                               ll.p_types, ll.min_ll, ll.trend, R_NilValue,
                               ll.varying, ll.pt_context);
        } catch (Rcpp::LongjumpException&) {
          throw;
        } catch (Rcpp::internal::InterruptedException&) {
          throw;
        } catch (...) {
          throw BlockAbort();
        }
        if (ll_eval.size() != n_rows - first) throw BlockAbort();
        std::vector<double> lw(n_rows);
        if (have_step_ref) lw[0] = step_ref_ll;
        for (int i = first; i < n_rows; ++i) lw[i] = ll_eval[i - first];
        for (double& v : lw) if (!R_FINITE(v)) v = R_NegInf;

        const arma::vec lp = fast_dmvnorm_rooti(P, plan.group_mu, rooti[0], log_consts[0]);
        const int m = static_cast<int>(ks.size());
        std::vector<arma::vec> comps(m);
        for (int j = 0; j < m; ++j) {
          const int k = ks[j];
          const double lmix = std::log(pm.mix[k]);
          arma::vec dens = (k == 0) ? lp
            : fast_dmvnorm_rooti(P, mus_i[k], rooti[k], log_consts[k]);
          for (int i = 0; i < n_rows; ++i) dens[i] = lmix + dens[i];
          comps[j] = dens;
        }
        std::vector<double> lm(n_rows);
        for (int i = 0; i < n_rows; ++i) {
          double max_log = comps[0][i];
          for (int j = 1; j < m; ++j) {
            const double v = comps[j][i];
            // pmax(): NA if either is NA.
            if (ISNAN(max_log) || ISNAN(v)) max_log = ISNAN(max_log) ? max_log : v;
            else if (v > max_log) max_log = v;
          }
          long double acc = 0.0L;
          for (int j = 0; j < m; ++j) acc += std::exp(comps[j][i] - max_log);
          lm[i] = max_log + std::log(static_cast<double>(acc));
        }
        bool any_bad = false;
        double fin_min = R_PosInf;
        bool any_fin = false;
        for (double v : lm) {
          if (!R_FINITE(v)) any_bad = true;
          else { any_fin = true; if (v < fin_min) fin_min = v; }
        }
        if (any_bad) {
          const double repl = any_fin ? fin_min : -1e10;
          for (double& v : lm) if (!R_FINITE(v)) v = repl;
        }
        std::vector<double> aux(n_rows, 0.0);
        for (int k : ks) {
          if (!plan.state_centred[k]) continue;
          const arma::vec dens = fast_dmvnorm_rooti(P, mus_i[k], rooti[k], log_consts[k]);
          for (int i = 0; i < n_rows; ++i) aux[i] = aux[i] + dens[i];
        }
        std::vector<double> l(n_rows);
        for (int i = 0; i < n_rows; ++i) {
          const double lw_total = (lw[i] + ref_ll) - lw[0];
          l[i] = ((lw_total + lp[i]) + aux[i]) - lm[i];
        }
        const double l_max = r_max(l);
        std::vector<double> w(n_rows);
        for (int i = 0; i < n_rows; ++i) w[i] = std::exp(l[i] - l_max);
        const int pick = r_sample_one(w);
        sel.steps += 1;
        int comp = -1, offset = 1;
        for (int k : ks) {
          const int n_k = particle_numbers[k];
          if (pick >= offset && pick < offset + n_k) comp = k;
          offset += n_k;
        }
        if (comp >= 0) {
          sel.n[comp] = sel.n[comp] + 1;
          std::vector<double> d(p);
          long double ss = 0.0L;
          for (int j = 0; j < p; ++j) {
            d[j] = P(pick, j) - ref[j];
            ss += d[j] * d[j];
          }
          sel.disp[comp] = sel.disp[comp] + std::sqrt(static_cast<double>(ss));
          if (plan.have_jump) {
            std::vector<double> dj(p);
            r_matprod(d.data(), 1, p, plan.jump_rooti.memptr(), p, p, dj.data());
            long double sj = 0.0L;
            for (int j = 0; j < p; ++j) sj += dj[j] * dj[j];
            jump[comp] = jump[comp] + static_cast<double>(sj);
          }
        }
        offset = 1;
        for (int k : ks) {
          const int n_k = particle_numbers[k];
          if (n_k > 0) {
            int nb = 0;
            for (int i = offset; i < offset + n_k; ++i) if (w[i] > w[0]) ++nb;
            better[k] = nb;
            offset += n_k;
          }
        }
        // The independence step's weight ESS.
        if (std::find(ks.begin(), ks.end(), 0) != ks.end()) {
          std::vector<double> w2(n_rows);
          for (int i = 0; i < n_rows; ++i) w2[i] = w[i] * w[i];
          const double sw = r_sum(w);
          weight_ess = sw * sw / r_sum(w2);
        }
        for (int j = 0; j < p; ++j) ref[j] = P(pick, j);
        ref_ll = lw[pick];
        step_ref_ll = lw[pick];
        have_step_ref = true;
      }
      if (adapt) {
        update_pm(pm, better, particle_numbers, weight_ess, jump, tune, p, plan.layout);
        if (pm.iter > tune.n0 && pm.mix.size() > 3 && gd_good && weight_ess > 0) {
          const double scale = R_pow(tune.target_ESS / weight_ess, tune.ESS_scale);
          const double nn = nearbyint(pm.n_particles * scale);
          pm.n_particles = std::max(25.0, std::min(tune.max_particles, nn));
          pm.particles_set = true;
        }
      }
      state = ref;
      state_ll = ref_ll;
    } catch (BlockAbort&) {
      sel = sel_before;
      pm = pm_before;
      seed_at_abort = seed;
      break;
    }
    for (int j = 0; j < p; ++j) alpha_out[j + static_cast<size_t>(it) * p] = state[j];
    ll_out[it] = state_ll;
    done = it + 1;
  }
  return done;
}

PmState read_pm(const List& pm_in) {
  PmState pm;
  pm.iter = as<double>(pm_in["iter"]);
  pm.proposal_counts = as<std::vector<double>>(pm_in["proposal_counts"]);
  pm.acc_counts = as<std::vector<double>>(pm_in["acc_counts"]);
  pm.jump_counts = as<std::vector<double>>(pm_in["jump_counts"]);
  pm.epsilon = as<std::vector<double>>(pm_in["epsilon"]);
  pm.mix = as<std::vector<double>>(pm_in["mix"]);
  pm.n_particles = as<double>(pm_in["n_particles"]);
  if (!Rf_isNull(pm_in["step_layout"])) {
    pm.step_layout = as<std::vector<std::string>>(pm_in["step_layout"]);
    pm.has_layout = true;
  }
  return pm;
}

List write_pm(const PmState& pm) {
  return List::create(
    _["iter"] = pm.iter, _["proposal_counts"] = wrap(pm.proposal_counts),
    _["acc_counts"] = wrap(pm.acc_counts), _["jump_counts"] = wrap(pm.jump_counts),
    _["epsilon"] = wrap(pm.epsilon), _["mix"] = wrap(pm.mix),
    _["step_layout"] = wrap(pm.step_layout), _["n_particles"] = pm.n_particles,
    _["weight_ess"] = pm.weight_ess, _["weight_set"] = pm.weight_set,
    _["layout_set"] = pm.layout_set, _["particles_set"] = pm.particles_set);
}

TuneArgs read_tune(const List& t) {
  return TuneArgs{as<double>(t["n0"]), as<double>(t["step_target"]),
                  as<double>(t["mix_adapt"]), as<double>(t["target_ESS"]),
                  as<double>(t["ESS_scale"]), as<double>(t["max_particles"])};
}

// chol(): the upper factor by dpotrf on a copy with a zeroed lower triangle.
// False when the matrix is not positive definite.
bool r_chol(const arma::mat& A, arma::mat& out) {
  int n = static_cast<int>(A.n_rows);
  out = A;
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) out(i, j) = 0.0;
  int info = 0;
  F77_CALL(dpotrf)("U", &n, out.memptr(), &n, &info FCONE);
  return info == 0;
}

// chol2inv(): dpotri on the upper triangle, then mirrored.
bool r_chol2inv(const arma::mat& Rf, arma::mat& out) {
  int n = static_cast<int>(Rf.n_rows);
  out.set_size(n, n);
  out.zeros();
  for (int j = 0; j < n; ++j)
    for (int i = 0; i <= j; ++i) out(i, j) = Rf(i, j);
  int info = 0;
  F77_CALL(dpotri)("U", &n, out.memptr(), &n, &info FCONE);
  if (info != 0) return false;
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) out(i, j) = out(j, i);
  return true;
}

// backsolve(R, diag(p)).
bool r_backsolve_identity(const arma::mat& Rf, arma::mat& out) {
  int n = static_cast<int>(Rf.n_rows);
  for (int i = 0; i < n; ++i) if (Rf(i, i) == 0.0) return false;
  out.eye(n, n);
  const double one = 1.0;
  F77_CALL(dtrsm)("L", "U", "N", "N", &n, &n, &one,
                  const_cast<double*>(Rf.memptr()), &n, out.memptr(), &n
                  FCONE FCONE FCONE FCONE);
  return true;
}

// A named element that is present and not NULL.
bool has(const List& x, const char* name) {
  return x.containsElementNamed(name) && !Rf_isNull(x[name]);
}

bool all_finite(const arma::mat& m) {
  for (arma::uword i = 0; i < m.n_elem; ++i) if (!R_FINITE(m[i])) return false;
  return true;
}

// A subject's HMC transition leaves the native path (R then redoes it from
// the subject's stream) when an operand is not finite or a likelihood call
// fails.
struct HmcAbort {};

void need_finite(const double* x, int n) {
  for (int i = 0; i < n; ++i) if (!R_FINITE(x[i])) throw HmcAbort();
}

// The whitened frame of one subject (.emc_hmc_frame) and its likelihood.
struct HmcFrame {
  int p = 0;
  const double* mu = nullptr;
  const arma::mat* Vi = nullptr;
  arma::mat R, S;
  SEXP dadm = R_NilValue, model = R_NilValue;
  double h = 0.0, floor = 0.0;
  // Rows for a centre plus p probes, and for the probes alone.
  NumericMatrix x_centre, x_probes;
};

struct HmcPoint {
  bool ok = false;
  int rows = 0;
  std::vector<double> a, z, g;
  double ll = 0.0, lp = 0.0;
};

// .emc_hmc_eval(): log target, whitened gradient and likelihood at z.
void hmc_eval(const std::vector<double>& z, const HmcFrame& fr, bool have_centre,
              double ll_centre, HmcPoint& out) {
  const int p = fr.p;
  need_finite(z.data(), p);
  out.z = z;
  out.a.resize(p);
  std::vector<double> sz(p);
  r_matprod(fr.S.memptr(), p, p, z.data(), p, 1, sz.data());
  for (int j = 0; j < p; ++j) out.a[j] = fr.mu[j] + sz[j];
  const int off = have_centre ? 0 : 1;
  const NumericMatrix& X = have_centre ? fr.x_probes : fr.x_centre;
  const int rows = p + off;
  double* x = REAL(X);
  for (int j = 0; j < p; ++j) {
    if (off) x[static_cast<size_t>(j) * rows] = out.a[j];
    for (int i = 0; i < p; ++i)
      x[(off + i) + static_cast<size_t>(j) * rows] = out.a[j] + fr.S(j, i) * fr.h;
  }
  out.rows = rows;
  out.ok = false;
  SEXP ll_sexp;
  try {
    ll_sexp = emc2_ll_cached(X, fr.dadm, fr.model, R_NilValue);
  } catch (Rcpp::LongjumpException&) {
    throw;
  } catch (Rcpp::internal::InterruptedException&) {
    throw;
  } catch (...) {
    throw HmcAbort();
  }
  if (TYPEOF(ll_sexp) != REALSXP || XLENGTH(ll_sexp) != rows) throw HmcAbort();
  const Shield<SEXP> keep(ll_sexp);
  const double* ll = REAL(ll_sexp);
  for (int i = 0; i < rows; ++i)
    if (!R_FINITE(ll[i]) || ll[i] <= fr.floor) return;
  if (have_centre && !(R_FINITE(ll_centre) && ll_centre > fr.floor)) return;
  const double lc = have_centre ? ll_centre : ll[0];
  std::vector<double> d(p), vid(p), cp(p);
  for (int j = 0; j < p; ++j) d[j] = out.a[j] - fr.mu[j];
  need_finite(d.data(), p);
  r_matprod(fr.Vi->memptr(), p, p, d.data(), p, 1, vid.data());
  need_finite(vid.data(), p);
  {
    // crossprod(S, Vid).
    const char trans_t = 'T';
    const double one = 1.0, zero = 0.0;
    const arma::blas_int ione = 1;
    arma::blas_int pb = p;
    arma::blas::gemv<double>(&trans_t, &pb, &pb, &one, fr.S.memptr(), &pb,
                             vid.data(), &ione, &zero, cp.data(), &ione);
  }
  out.g.resize(p);
  for (int i = 0; i < p; ++i) out.g[i] = (ll[off + i] - lc) / fr.h - cp[i];
  long double q = 0.0L;
  for (int j = 0; j < p; ++j) q += d[j] * vid[j];
  out.ll = lc;
  out.lp = lc - 0.5 * static_cast<double>(q);
  out.ok = true;
}

double sum_squares(const std::vector<double>& v) {
  long double s = 0.0L;
  for (double x : v) s += x * x;
  return static_cast<double>(s);
}

// R's mean() of doubles: a long double mean with one correction pass.
double r_mean(const std::vector<double>& v) {
  const size_t n = v.size();
  long double s = 0.0L;
  for (double x : v) s += x;
  s /= n;
  long double t = 0.0L;
  for (double x : v) t += (x - s);
  s += t / n;
  return static_cast<double>(s);
}

}  // namespace

// Subject HMC transitions (.emc_hmc_run) for the subjects of one request, each
// from its own stream. `done` marks the subjects that ran here; the caller
// runs the others in R from the streams it sent.
// [[Rcpp::export]]
List emc2_native_hmc_subjects(NumericMatrix alpha, NumericVector ll,
                              NumericMatrix mu, NumericMatrix var, List seeds,
                              NumericVector eps_base, List plans, List data,
                              SEXP model, CharacterVector par_names) {
  const int n = alpha.ncol(), p = alpha.nrow();
  LogicalVector done(n);
  List tasks(n), seeds_out(n);
  if (mu.nrow() != p || var.nrow() != p || var.ncol() != p ||
      par_names.size() != p)
    return List::create(_["done"] = done, _["tasks"] = tasks, _["seeds"] = seeds_out);
  const arma::mat V(var.begin(), p, p, false);
  arma::mat Vc, Vi;
  if (!all_finite(V) || !r_chol(V, Vc) || !r_chol2inv(Vc, Vi))
    return List::create(_["done"] = done, _["tasks"] = tasks, _["seeds"] = seeds_out);
  Environment global = Environment::global_env();
  const SEXP alpha_dn = Rf_getAttrib(alpha, R_DimNamesSymbol);
  const SEXP alpha_names = (TYPEOF(alpha_dn) == VECSXP) ? VECTOR_ELT(alpha_dn, 0)
                                                        : R_NilValue;
  const SEXP ll_names = Rf_getAttrib(ll, R_NamesSymbol);
  const List x_dimnames = List::create(R_NilValue, par_names);

  for (int k = 0; k < n; ++k) {
    done[k] = false;
    if (Rf_isNull(plans[k]) || Rf_isNull(seeds[k])) continue;
    const std::clock_t c0 = std::clock();
    const List plan = plans[k];
    if (!has(plan, "prec") || !has(plan, "h") || !has(plan, "floor")) continue;
    const SEXP prec_sexp = plan["prec"];
    if (!Rf_isMatrix(prec_sexp) || TYPEOF(prec_sexp) != REALSXP ||
        Rf_nrows(prec_sexp) != p || Rf_ncols(prec_sexp) != p) continue;
    const arma::mat prec(REAL(prec_sexp), p, p, false);
    const double eps0 = eps_base[k];
    if (!all_finite(prec) || !R_FINITE(eps0) || !(eps0 > 0.0)) continue;

    HmcFrame fr;
    fr.p = p;
    fr.mu = &mu(0, k);
    fr.Vi = &Vi;
    fr.dadm = data[k];
    fr.model = model;
    fr.h = as<double>(plan["h"]);
    fr.floor = as<double>(plan["floor"]);
    if (!emc2_ll_cache_has(fr.dadm, fr.model, par_names)) continue;
    {
      arma::mat Q = prec + Vi, Qs(p, p);
      for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) Qs(i, j) = (Q(i, j) + Q(j, i)) / 2;
      if (!r_chol(Qs, fr.R) || !all_finite(fr.R) ||
          !r_backsolve_identity(fr.R, fr.S) || !all_finite(fr.S)) continue;
    }
    const bool have_W = has(plan, "W");
    arma::mat W;
    if (have_W) {
      const SEXP w_sexp = plan["W"];
      if (!Rf_isMatrix(w_sexp) || TYPEOF(w_sexp) != REALSXP ||
          Rf_ncols(w_sexp) != p) continue;
      W = arma::mat(REAL(w_sexp), Rf_nrows(w_sexp), p, false);
      if (!all_finite(W)) continue;
    }
    fr.x_centre = NumericMatrix(p + 1, p);
    fr.x_probes = NumericMatrix(p, p);
    fr.x_centre.attr("dimnames") = x_dimnames;
    fr.x_probes.attr("dimnames") = x_dimnames;

    std::vector<double> a0(p), mu_k(p);
    for (int j = 0; j < p; ++j) { a0[j] = alpha(j, k); mu_k[j] = mu(j, k); }
    const double ll0 = ll[k];

    global.assign(".Random.seed", seeds[k]);
    GetRNGstate();
    bool aborted = false;
    // Draws, in order: the jitter, p momenta, the acceptance uniform.
    const double U = R::runif(0.9, 1.1);
    std::vector<double> v0(p);
    for (int j = 0; j < p; ++j) v0[j] = R::rnorm(0.0, 1.0);
    const double u_acc = R::runif(0.0, 1.0);
    PutRNGstate();
    const double eps = eps0 * U;
    const double l_real = std::ceil((M_PI / 2) / eps);
    if (!(l_real < INT_MAX)) continue;
    const int L = std::max(1, static_cast<int>(l_real));

    int rows = 0, steps = 0;
    bool ok = false, accepted = false;
    double acc_prob = 0.0, energy = NA_REAL;
    HmcPoint cur, nw;
    std::vector<double> disp(p, 0.0);
    try {
      std::vector<double> d0(p), z(p), v(p);
      for (int j = 0; j < p; ++j) d0[j] = a0[j] - mu_k[j];
      need_finite(d0.data(), p);
      r_matprod(fr.R.memptr(), p, p, d0.data(), p, 1, z.data());
      hmc_eval(z, fr, true, ll0, cur);
      rows = cur.rows;
      if (cur.ok) {
        const double half = 0.5 * eps;
        for (int j = 0; j < p; ++j) v[j] = v0[j] + half * cur.g[j];
        z = cur.z;
        ok = true;
        for (int l = 1; l <= L; ++l) {
          for (int j = 0; j < p; ++j) z[j] = z[j] + eps * v[j];
          hmc_eval(z, fr, false, 0.0, nw);
          rows += nw.rows;
          steps = l;
          if (!nw.ok) { ok = false; break; }
          const double w = (l < L) ? eps : half;
          for (int j = 0; j < p; ++j) v[j] = v[j] + w * nw.g[j];
        }
        if (ok) {
          energy = (-nw.lp + 0.5 * sum_squares(v)) - (-cur.lp + 0.5 * sum_squares(v0));
          acc_prob = R_FINITE(energy) ? std::min(1.0, std::exp(-energy)) : 0.0;
          if (ISNAN(acc_prob)) throw HmcAbort();
          for (int j = 0; j < p; ++j) disp[j] = nw.a[j] - a0[j];
          accepted = u_acc < acc_prob;
        }
      }
    } catch (HmcAbort&) {
      aborted = true;
    }
    if (aborted) continue;

    double score = NA_REAL;
    if (have_W) {
      need_finite(disp.data(), p);
      std::vector<double> wd(W.n_rows);
      r_matprod(W.memptr(), static_cast<int>(W.n_rows), p, disp.data(), p, 1,
                wd.data());
      for (double& x : wd) x = x * x;
      score = acc_prob * r_mean(wd);
    }

    // The task as .emc_hmc_run() builds it, names included.
    NumericVector alpha_k(p), disp_k(p), ll_k(1), energy_k(1);
    if (accepted) {
      for (int j = 0; j < p; ++j) alpha_k[j] = nw.a[j];
      ll_k[0] = nw.ll;
    } else {
      for (int j = 0; j < p; ++j) alpha_k[j] = a0[j];
      if (!Rf_isNull(alpha_names)) alpha_k.attr("names") = alpha_names;
      ll_k[0] = ll0;
      if (!Rf_isNull(ll_names))
        ll_k.attr("names") = CharacterVector::create(String(STRING_ELT(ll_names, k)));
    }
    for (int j = 0; j < p; ++j) disp_k[j] = disp[j];
    if (ok && !Rf_isNull(alpha_names)) disp_k.attr("names") = alpha_names;
    energy_k[0] = energy;
    if (ok && !Rf_isNull(ll_names))
      energy_k.attr("names") = CharacterVector::create(String(STRING_ELT(ll_names, k)));
    List task = List::create(
      _["attempted"] = true, _["rows"] = rows,
      _["status"] = ok ? "ok" : "abort", _["alpha"] = alpha_k, _["ll"] = ll_k,
      _["accepted"] = accepted, _["acc_prob"] = acc_prob, _["disp"] = disp_k,
      _["U"] = U, _["eps"] = eps, _["L"] = L, _["steps"] = steps,
      _["energy"] = energy_k, _["score"] = score,
      _["cpu"] = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC);
    tasks[k] = task;
    seeds_out[k] = global.get(".Random.seed");
    done[k] = true;
  }
  return List::create(_["done"] = done, _["tasks"] = tasks, _["seeds"] = seeds_out);
}

// [[Rcpp::export]]
List emc2_native_particle_block(List plan_in, int n_iter, NumericVector start,
                                double start_ll, List ll_args,
                                NumericVector sel_n, NumericVector sel_disp,
                                double sel_steps, List pm_in, List tune_in,
                                bool adapt) {
  const List Mus = plan_in["Mus"];
  const List base_R = plan_in["base_R"], base_rooti = plan_in["base_rooti"];
  const NumericVector base_sld = plan_in["base_sld"];
  const List steps = plan_in["steps"];
  const LogicalVector state_centred = plan_in["state_centred"];
  const CharacterVector layout_in = plan_in["layout"];
  const RObject jump_rooti_in = plan_in["jump_rooti"];

  Plan plan;
  plan.p = start.size();
  plan.K = Mus.size();
  const int K = plan.K, p = plan.p;
  plan.log_const0 = as<double>(plan_in["const"]);
  plan.multiplier = as<double>(plan_in["multiplier"]);
  plan.group_mu = as<arma::rowvec>(plan_in["group_mu"]);
  const LlArgs ll(ll_args);

  PmState pm = read_pm(pm_in);
  if (static_cast<int>(pm.mix.size()) != K || pm.epsilon.empty() ||
      static_cast<int>(pm.proposal_counts.size()) != K ||
      static_cast<int>(pm.acc_counts.size()) != K ||
      layout_in.size() != K - 1)
    stop("native particle plan is inconsistent");
  const TuneArgs tune = read_tune(tune_in);
  const bool gd_good = as<bool>(pm_in["gd_good"]);
  plan.layout.resize(K - 1);
  for (int c = 0; c < K - 1; ++c) plan.layout[c] = as<std::string>(layout_in[c]);
  plan.mu.resize(K); plan.bR.resize(K); plan.brooti.resize(K);
  plan.bsld.assign(base_sld.begin(), base_sld.end());
  plan.state_centred.resize(K);
  for (int k = 0; k < K; ++k) {
    plan.bR[k] = as<arma::mat>(base_R[k]);
    plan.brooti[k] = as<arma::mat>(base_rooti[k]);
    plan.state_centred[k] = state_centred[k];
    if (!Rf_isNull(Mus[k])) plan.mu[k] = as<arma::rowvec>(Mus[k]);
  }
  plan.have_jump = !jump_rooti_in.isNULL();
  if (plan.have_jump) plan.jump_rooti = as<arma::mat>(jump_rooti_in);
  for (int s = 0; s < steps.size(); ++s) {
    IntegerVector ks = steps[s];
    std::vector<int> v;
    for (int k : ks) v.push_back(k - 1);
    plan.step_ks.push_back(v);
  }

  NumericMatrix alpha_out(p, n_iter);
  NumericVector ll_out(n_iter);
  SelState sel;
  sel.n.assign(sel_n.begin(), sel_n.end());
  sel.disp.assign(sel_disp.begin(), sel_disp.end());
  sel.steps = sel_steps;
  std::vector<double> state(start.begin(), start.end());
  double state_ll = start_ll;

  GetRNGstate();
  RObject seed_at_abort;
  const int done = run_iterations(plan, ll, n_iter, state, state_ll, pm, tune,
                                  gd_good, adapt, sel, alpha_out.begin(),
                                  ll_out.begin(), seed_at_abort);
  PutRNGstate();
  // On an abort the caller restores the stream to that iteration's start
  // (the exported wrapper writes the advanced stream back on return).
  return List::create(_["alpha"] = alpha_out, _["ll"] = ll_out,
                      _["done"] = done, _["sel_n"] = wrap(sel.n),
                      _["sel_disp"] = wrap(sel.disp), _["sel_steps"] = sel.steps,
                      _["pm"] = write_pm(pm),
                      _["seed"] = done < n_iter ? SEXP(seed_at_abort) : R_NilValue);
}

// One particle step for each subject of a multi-subject iteration. Each
// subject's plan is built here from the group state and its block-fixed
// parts, then stepped on the subject's own stream. A subject returns
// done = FALSE (its stream and settings untouched) when the R step must run:
// a factorisation fails or the step hits a case R handles.
//
// statics[[k]]: ll_args, chains_mu, chains (R, rooti, sld), eff_mu, eff_B,
// eff_given, eff (R, rooti, sld), lik_prec, lik_lin. group: R, rooti, sld.
// [[Rcpp::export]]
List emc2_native_particle_subjects(List statics, NumericMatrix alpha,
                                   NumericVector prev_ll,
                                   NumericMatrix population_mu,
                                   NumericMatrix group_var, List group,
                                   List pm_list, List seeds, LogicalVector skip,
                                   std::string stage, List tune_in) {
  const int n = alpha.ncol(), p = alpha.nrow();
  const bool adapt = stage != "sample";
  const TuneArgs tune = read_tune(tune_in);
  const arma::mat gR = as<arma::mat>(group["R"]), grooti = as<arma::mat>(group["rooti"]);
  const double gsld = as<double>(group["sld"]);
  const arma::mat V(group_var.begin(), p, p, false);
  const double log_const0 = -0.5 * p * std::log(2 * M_PI);
  // The group precision, shared by every subject's Laplace conditional.
  arma::mat Vi;
  bool have_Vi = false, tried_Vi = false;

  NumericMatrix alpha_out(p, n);
  NumericVector ll_out(n), cpu(n);
  LogicalVector done(n);
  List pm_out(n), seeds_out(n), sel_out(n);
  Environment global = Environment::global_env();

  for (int k = 0; k < n; ++k) {
    done[k] = false;
    if (skip[k] || Rf_isNull(statics[k])) continue;
    const std::clock_t c0 = std::clock();
    const List st = statics[k];
    const List pm_all = pm_list[k];
    const List pm_r = pm_all[0];
    const int mix_len = Rf_length(pm_r["mix"]);

    Plan plan;
    plan.p = p;
    plan.log_const0 = log_const0;
    plan.group_mu = arma::rowvec(p);
    for (int j = 0; j < p; ++j) plan.group_mu[j] = population_mu(j, k);

    // The efficient proposal at this group state.
    bool have_eff = false, ok = true;
    arma::rowvec eff_mu;
    arma::mat effR, effrooti;
    double effsld = 0.0;
    const bool want_eff = !(stage == "preburn" || (stage == "burn" && mix_len < 4));
    if (want_eff && !Rf_isNull(st["lik_prec"])) {
      if (!tried_Vi) {
        tried_Vi = true;
        arma::mat cv;
        have_Vi = all_finite(V) && r_chol(V, cv) && r_chol2inv(cv, Vi);
      }
      if (!have_Vi) continue;
      const arma::mat prec = as<arma::mat>(st["lik_prec"]);
      const arma::vec lin = as<arma::vec>(st["lik_lin"]);
      arma::mat Q = prec + Vi, Qs(p, p), Rq, Sig;
      for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) Qs(i, j) = (Q(i, j) + Q(j, i)) / 2;
      if (!all_finite(Qs) || !r_chol(Qs, Rq) || !r_chol2inv(Rq, Sig)) continue;
      std::vector<double> gm(p), t1(p), rhs(p), m(p);
      for (int j = 0; j < p; ++j) gm[j] = population_mu(j, k);
      r_matprod(Vi.memptr(), p, p, gm.data(), p, 1, t1.data());
      for (int j = 0; j < p; ++j) rhs[j] = lin[j] + t1[j];
      r_matprod(Sig.memptr(), p, p, rhs.data(), p, 1, m.data());
      arma::mat var(p, p);
      for (int j = 0; j < p; ++j)
        for (int i = 0; i < p; ++i) var(i, j) = (Sig(i, j) + Sig(j, i)) / 2;
      for (int j = 0; j < p; ++j) if (!R_FINITE(m[j])) ok = false;
      if (!ok || !all_finite(Sig)) continue;
      // .chol_factor() of the conditional covariance; a repaired factor is
      // left to R.
      if (!r_chol(var, effR) || !r_backsolve_identity(effR, effrooti)) continue;
      long double sl = 0.0L;
      for (int j = 0; j < p; ++j) sl += std::log(effrooti(j, j));
      effsld = static_cast<double>(sl);
      eff_mu = arma::rowvec(m);
      have_eff = true;
    } else if (want_eff && !Rf_isNull(st["eff_mu"])) {
      const NumericVector em = st["eff_mu"];
      eff_mu = arma::rowvec(p);
      for (int j = 0; j < p; ++j) eff_mu[j] = em[j];
      if (!Rf_isNull(st["eff_B"])) {
        // .emc_eff_mean(): the mean follows the group mean.
        const NumericMatrix B = st["eff_B"];
        const NumericVector given = st["eff_given"];
        std::vector<double> delta(p), shift(p);
        for (int j = 0; j < p; ++j) {
          delta[j] = population_mu(j, k) - given[j];
          if (!R_FINITE(delta[j])) delta[j] = 0.0;
        }
        r_matprod(B.begin(), p, p, delta.data(), p, 1, shift.data());
        bool fin = true;
        std::vector<double> out(p);
        for (int j = 0; j < p; ++j) {
          out[j] = em[j] + shift[j];
          if (!R_FINITE(out[j])) fin = false;
        }
        if (fin) for (int j = 0; j < p; ++j) eff_mu[j] = out[j];
      }
      if (Rf_isNull(st["eff"])) continue;
      const List ef = st["eff"];
      effR = as<arma::mat>(ef["R"]);
      effrooti = as<arma::mat>(ef["rooti"]);
      effsld = as<double>(ef["sld"]);
      have_eff = true;
    }

    // The stage's proposal layout (new_particle()).
    arma::mat cR, crooti;
    double csld = 0.0;
    arma::rowvec chains_mu;
    const bool have_chains = !Rf_isNull(st["chains"]);
    if (have_chains) {
      const List ch = st["chains"];
      cR = as<arma::mat>(ch["R"]);
      crooti = as<arma::mat>(ch["rooti"]);
      csld = as<double>(ch["sld"]);
    }
    if (!Rf_isNull(st["chains_mu"])) chains_mu = as<arma::rowvec>(st["chains_mu"]);
    // tags: 0 group, 1 chains, 2 eff.
    std::vector<int> tags;
    std::vector<int> centre;     // 0 group mean, 1 state, 2 chains mean, 3 eff mean
    if (stage == "preburn") {
      tags = {0, 0}; centre = {0, 1}; plan.multiplier = 2.0;
    } else if (stage == "burn" && mix_len < 4) {
      tags = {0, 0, 1}; centre = {0, 1, 1};
    } else if (stage == "adapt" && (!have_eff || mix_len < 4)) {
      tags = {0, 1, 1}; centre = {0, 1, 2};
    } else {
      tags = {0, 1, 1, 2}; centre = {0, 1, 2, 3};
    }
    const int K = static_cast<int>(tags.size());
    plan.K = K;
    plan.mu.resize(K); plan.bR.resize(K); plan.brooti.resize(K);
    plan.bsld.resize(K); plan.state_centred.resize(K);
    bool uses_chains = false;
    for (int j = 0; j < K; ++j) {
      plan.state_centred[j] = centre[j] == 1;
      if (tags[j] == 0) { plan.bR[j] = gR; plan.brooti[j] = grooti; plan.bsld[j] = gsld; }
      else if (tags[j] == 1) {
        if (!have_chains) { ok = false; break; }
        uses_chains = true;
        plan.bR[j] = cR; plan.brooti[j] = crooti; plan.bsld[j] = csld;
      } else {
        if (!have_eff) { ok = false; break; }
        plan.bR[j] = effR; plan.brooti[j] = effrooti; plan.bsld[j] = effsld;
      }
      if (centre[j] == 0) plan.mu[j] = plan.group_mu;
      else if (centre[j] == 2) {
        if (static_cast<int>(chains_mu.n_elem) != p) { ok = false; break; }
        plan.mu[j] = chains_mu;
      } else if (centre[j] == 3) plan.mu[j] = eff_mu;
      if (!plan.state_centred[j])
        for (int q = 0; q < p; ++q) if (!R_FINITE(plan.mu[j][q])) ok = false;
      if (!all_finite(plan.bR[j]) || !all_finite(plan.brooti[j]) ||
          !R_FINITE(plan.bsld[j])) ok = false;
    }
    if (!ok) continue;
    plan.have_jump = true;
    plan.jump_rooti = uses_chains ? crooti : grooti;
    static const char* tag_names[] = {"group", "chains", "eff"};
    plan.layout.resize(K - 1);
    for (int j = 1; j < K; ++j)
      plan.layout[j - 1] = plan.state_centred[j] ? tag_names[tags[j]] : "";
    std::vector<int> indep;
    for (int j = 0; j < K; ++j) if (!plan.state_centred[j]) indep.push_back(j);
    if (!indep.empty()) plan.step_ks.push_back(indep);
    for (int j = 0; j < K; ++j)
      if (plan.state_centred[j]) plan.step_ks.push_back(std::vector<int>{j});

    // The subject's settings, as R stores them.
    PmState pm;
    bool gd_good = false;
    SelState sel;
    try {
      if (mix_len != K || Rf_length(pm_r["epsilon"]) < 1 ||
          Rf_length(pm_r["n_particles"]) != 1) continue;
      pm.mix = as<std::vector<double>>(pm_r["mix"]);
      pm.epsilon = as<std::vector<double>>(pm_r["epsilon"]);
      pm.n_particles = as<double>(pm_r["n_particles"]);
      if (has(pm_r, "jump_counts"))
        pm.jump_counts = as<std::vector<double>>(pm_r["jump_counts"]);
      if (adapt) {
        if (!has(pm_r, "proposal_counts") || !has(pm_r, "acc_counts") ||
            !has(pm_r, "iter") || !has(pm_r, "gd_good") ||
            Rf_length(pm_r["proposal_counts"]) != K ||
            Rf_length(pm_r["acc_counts"]) != K || Rf_length(pm_r["iter"]) != 1 ||
            Rf_length(pm_r["gd_good"]) != 1) continue;
        pm.iter = as<double>(pm_r["iter"]);
        pm.proposal_counts = as<std::vector<double>>(pm_r["proposal_counts"]);
        pm.acc_counts = as<std::vector<double>>(pm_r["acc_counts"]);
        gd_good = as<bool>(pm_r["gd_good"]);
      } else {
        pm.proposal_counts.assign(K, 0.0);
        pm.acc_counts.assign(K, 0.0);
      }
      if (has(pm_r, "step_layout")) {
        pm.step_layout = as<std::vector<std::string>>(pm_r["step_layout"]);
        pm.has_layout = true;
      }
      sel.n.assign(K, 0.0);
      sel.disp.assign(K, 0.0);
      if (has(pm_r, "sel")) {
        const List sel_all = pm_r["sel"];
        if (has(sel_all, stage.c_str())) {
          const List sl = sel_all[stage];
          const NumericVector sn = sl["n"], sd = sl["disp"];
          for (int j = 0; j < K && j < sn.size(); ++j) sel.n[j] = sn[j];
          for (int j = 0; j < K && j < sd.size(); ++j) sel.disp[j] = sd[j];
          if (sn.size() > K || sd.size() > K) continue;
          sel.steps = as<double>(sl["steps"]);
        }
      }
    } catch (...) {
      continue;
    }
    if (!R_FINITE(prev_ll[k])) continue;

    const LlArgs ll(as<List>(st["ll_args"]));
    std::vector<double> state(p);
    for (int j = 0; j < p; ++j) state[j] = alpha(j, k);
    double state_ll = prev_ll[k];
    const double steps_before = sel.steps;
    global.assign(".Random.seed", seeds[k]);
    GetRNGstate();
    RObject seed_at_abort;
    const int n_done = run_iterations(plan, ll, 1, state, state_ll, pm, tune,
                                      gd_good, adapt, sel,
                                      &alpha_out(0, k), &ll_out[k], seed_at_abort);
    PutRNGstate();
    if (n_done != 1) continue;
    done[k] = true;
    seeds_out[k] = global.get(".Random.seed");
    pm_out[k] = write_pm(pm);
    sel_out[k] = List::create(_["n"] = wrap(sel.n), _["disp"] = wrap(sel.disp),
                              _["steps"] = sel.steps,
                              _["moved"] = sel.steps > steps_before);
    cpu[k] = static_cast<double>(std::clock() - c0) / CLOCKS_PER_SEC;
  }
  return List::create(_["alpha"] = alpha_out, _["ll"] = ll_out, _["done"] = done,
                      _["pm"] = pm_out, _["sel"] = sel_out,
                      _["seeds"] = seeds_out, _["cpu"] = cpu);
}
