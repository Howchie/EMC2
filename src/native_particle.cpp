// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R_ext/Utils.h>
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

}  // namespace

// [[Rcpp::export]]
List emc2_native_particle_block(List plan, int n_iter, NumericVector start,
                                double start_ll, List ll_args,
                                NumericVector sel_n, NumericVector sel_disp,
                                double sel_steps, List pm_in, List tune_in,
                                bool adapt) {
  const List Mus = plan["Mus"];
  const List base_R = plan["base_R"], base_rooti = plan["base_rooti"];
  const NumericVector base_sld = plan["base_sld"];
  const double log_const0 = as<double>(plan["const"]);
  const double multiplier = as<double>(plan["multiplier"]);
  const List steps = plan["steps"];
  const LogicalVector state_centred = plan["state_centred"];
  const arma::rowvec group_mu = as<arma::rowvec>(plan["group_mu"]);
  const CharacterVector par_names = ll_args["par_names"];
  const CharacterVector layout_in = plan["layout"];
  const RObject jump_rooti_in = plan["jump_rooti"];

  const DataFrame data = ll_args["data"];
  const NumericVector constants = ll_args["constants"];
  const List designs = ll_args["designs"];
  const String type = as<String>(ll_args["type"]);
  const List bounds = ll_args["bounds"], transforms = ll_args["transforms"],
    pretransforms = ll_args["pretransforms"];
  const CharacterVector p_types = ll_args["p_types"];
  const double min_ll = as<double>(ll_args["min_ll"]);
  const Nullable<List> trend(ll_args["trend"]);
  const LogicalVector varying = ll_args["varying"];
  SEXP pt_context = ll_args["pt_context"];

  const int p = start.size();
  const int K = Mus.size();
  const double p_dbl = static_cast<double>(p);

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
  if (static_cast<int>(pm.mix.size()) != K || pm.epsilon.empty() ||
      static_cast<int>(pm.proposal_counts.size()) != K ||
      static_cast<int>(pm.acc_counts.size()) != K ||
      layout_in.size() != K - 1)
    stop("native particle plan is inconsistent");
  TuneArgs tune{as<double>(tune_in["n0"]), as<double>(tune_in["step_target"]),
                as<double>(tune_in["mix_adapt"]), as<double>(tune_in["target_ESS"]),
                as<double>(tune_in["ESS_scale"]), as<double>(tune_in["max_particles"])};
  const bool gd_good = as<bool>(pm_in["gd_good"]);
  std::vector<std::string> layout(K - 1);
  for (int c = 0; c < K - 1; ++c) layout[c] = as<std::string>(layout_in[c]);

  std::vector<int> particle_numbers(K);
  std::vector<arma::mat> bR(K), brooti(K), R(K), rooti(K);
  std::vector<double> log_consts(K);
  std::vector<arma::rowvec> mu(K);
  for (int k = 0; k < K; ++k) {
    bR[k] = as<arma::mat>(base_R[k]);
    brooti[k] = as<arma::mat>(base_rooti[k]);
    if (!Rf_isNull(Mus[k])) mu[k] = as<arma::rowvec>(Mus[k]);
  }
  const bool have_jump = !jump_rooti_in.isNULL();
  arma::mat jump_rooti;
  if (have_jump) jump_rooti = as<arma::mat>(jump_rooti_in);
  std::vector<std::vector<int>> step_ks;
  for (int s = 0; s < steps.size(); ++s) {
    IntegerVector ks = steps[s];
    std::vector<int> v;
    for (int k : ks) v.push_back(k - 1);
    step_ks.push_back(v);
  }

  NumericMatrix alpha_out(p, n_iter);
  NumericVector ll_out(n_iter);
  NumericVector n_sel = clone(sel_n), disp = clone(sel_disp);
  double steps_done = sel_steps;
  std::vector<double> state(start.begin(), start.end());
  double state_ll = start_ll;

  GetRNGstate();
  RObject seed_at_abort;
  int done = 0;
  std::vector<double> z, x;
  for (int it = 0; it < n_iter; ++it) {
    // The stream at the start of this iteration, for R to resume from.
    PutRNGstate();
    RObject seed = Environment::global_env().get(".Random.seed");
    std::vector<double> sel_n_before(n_sel.begin(), n_sel.end());
    std::vector<double> disp_before(disp.begin(), disp.end());
    const double steps_before = steps_done;
    const PmState pm_before = pm;
    try {
      // The proposals at this iteration's step sizes: independence
      // components keep scale one.
      for (int k = 0; k < K; ++k) {
        // rep_len(c(1, epsilon), K), then 1 for independence components.
        const size_t e = k % (pm.epsilon.size() + 1);
        const double eps_k = (e == 0 || !state_centred[k]) ? 1.0
          : pm.epsilon[e - 1];
        R[k] = bR[k] * eps_k;
        rooti[k] = brooti[k] / eps_k;
        log_consts[k] = (base_sld[k] - p_dbl * std::log(eps_k)) + log_const0;
      }
      // rmultinom(1, n_particles * multiplier, mix) after FixupProb.
      const double size_d = pm.n_particles * multiplier;
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
      for (const std::vector<int>& ks : step_ks) {
        int n_part = 0;
        for (int k : ks) n_part += particle_numbers[k];
        if (n_part == 0) continue;
        std::vector<arma::rowvec> mus_i(K);
        for (int k : ks) {
          if (state_centred[k]) {
            // Auxiliary centre z ~ N(ref, S_k).
            std::vector<double> centre(p);
            draw_rows(1, p, ref.data(), R[k], z, x, centre.data(), 1, 0);
            mus_i[k] = arma::rowvec(centre);
          } else {
            mus_i[k] = mu[k];
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
        colnames(pmat) = par_names;
        NumericVector ll_eval;
        try {
          ll_eval = calc_ll_oo(pmat, data, constants, designs, type, bounds,
                               transforms, pretransforms, p_types, min_ll,
                               trend, R_NilValue, varying, pt_context);
        } catch (...) {
          throw BlockAbort();
        }
        if (ll_eval.size() != n_rows - first) throw BlockAbort();
        std::vector<double> lw(n_rows);
        if (have_step_ref) lw[0] = step_ref_ll;
        for (int i = first; i < n_rows; ++i) lw[i] = ll_eval[i - first];
        for (double& v : lw) if (!R_FINITE(v)) v = R_NegInf;

        const arma::vec lp = fast_dmvnorm_rooti(P, group_mu, rooti[0], log_consts[0]);
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
          if (!state_centred[k]) continue;
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
        steps_done += 1;
        int comp = -1, offset = 1;
        for (int k : ks) {
          const int n_k = particle_numbers[k];
          if (pick >= offset && pick < offset + n_k) comp = k;
          offset += n_k;
        }
        if (comp >= 0) {
          n_sel[comp] = n_sel[comp] + 1;
          std::vector<double> d(p);
          long double ss = 0.0L;
          for (int j = 0; j < p; ++j) {
            d[j] = P(pick, j) - ref[j];
            ss += d[j] * d[j];
          }
          disp[comp] = disp[comp] + std::sqrt(static_cast<double>(ss));
          if (have_jump) {
            std::vector<double> dj(p);
            r_matprod(d.data(), 1, p, jump_rooti.memptr(), p, p, dj.data());
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
        update_pm(pm, better, particle_numbers, weight_ess, jump, tune, p, layout);
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
      std::copy(sel_n_before.begin(), sel_n_before.end(), n_sel.begin());
      std::copy(disp_before.begin(), disp_before.end(), disp.begin());
      steps_done = steps_before;
      pm = pm_before;
      seed_at_abort = seed;
      break;
    }
    for (int j = 0; j < p; ++j) alpha_out(j, it) = state[j];
    ll_out[it] = state_ll;
    done = it + 1;
  }
  PutRNGstate();
  List pm_out = List::create(
    _["iter"] = pm.iter, _["proposal_counts"] = wrap(pm.proposal_counts),
    _["acc_counts"] = wrap(pm.acc_counts), _["jump_counts"] = wrap(pm.jump_counts),
    _["epsilon"] = wrap(pm.epsilon), _["mix"] = wrap(pm.mix),
    _["step_layout"] = wrap(pm.step_layout), _["n_particles"] = pm.n_particles,
    _["weight_ess"] = pm.weight_ess, _["weight_set"] = pm.weight_set,
    _["layout_set"] = pm.layout_set, _["particles_set"] = pm.particles_set);
  // On an abort the caller restores the stream to that iteration's start
  // (the exported wrapper writes the advanced stream back on return).
  return List::create(_["alpha"] = alpha_out, _["ll"] = ll_out,
                      _["done"] = done, _["sel_n"] = n_sel,
                      _["sel_disp"] = disp, _["sel_steps"] = steps_done,
                      _["pm"] = pm_out,
                      _["seed"] = done < n_iter ? SEXP(seed_at_abort) : R_NilValue);
}
