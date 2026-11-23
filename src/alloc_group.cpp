#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

// The group coordinates priced as one linear system (R/move_allocation.R,
// .emc_alloc_coupling): an iteration acts on linear functions of the
// coordinates as A = M_q ... M_1 B, B = diag(fac) A_rest, with
// M_k = I - theta_k e_k r_k', r_k = Q[k, ] / Q[k, k], and coordinate j's
// integrated autocorrelation time is 2 [(I - A)^-1 C]_jj / C_jj - 1.

namespace {

double group_factor_k(double g, double a, double k) {
  if (a < g - 1e-12) return 1.0 - (1.0 - g) / (1.0 - a) * (1.0 - std::pow(a, k));
  return std::pow(g, k);
}

// The factor of x expected rounds (.emc_alloc_factor).
double group_factor(double g, double x, double a) {
  const double m = std::floor(x), f = x - m;
  return (1.0 - f) * group_factor_k(g, a, m) + f * group_factor_k(g, a, m + 1.0);
}

struct GroupSystem {
  arma::mat C, R, B;
  std::vector<double> g, a;
  int q;

  bool moves(int k) const { return std::isfinite(g[k]) && g[k] < 1.0 - 1e-9; }

  double theta(int k, double x) const {
    if (!moves(k) || !(x > 0.0)) return 0.0;
    const double th = 1.0 - group_factor(g[k], x, a[k]);
    return std::min(1.0, std::max(0.0, th));
  }
};

// The system at a plan: tau of every coordinate, and what a change of one
// coordinate's theta does to them (A changes by -d u_k v_k').
struct GroupState {
  arma::vec tau;
  arma::mat U, V, NU, VN, VNC;
  arma::vec nc, den;
  std::vector<double> th;
  bool ok;
};

double tau_to_rho(double tau) {
  if (!std::isfinite(tau) || tau <= 0.0) return 1.0 - 1e-12;
  const double rho = (tau - 1.0) / (tau + 1.0);
  return std::min(1.0 - 1e-12, std::max(0.0, rho));
}

GroupState group_state(const GroupSystem& sys, const std::vector<double>& x, bool trials) {
  const int q = sys.q;
  GroupState st;
  st.ok = false;
  st.tau.set_size(q);
  st.tau.fill(std::numeric_limits<double>::infinity());
  std::vector<double>& th = st.th;
  th.resize(q);
  for (int k = 0; k < q; ++k) th[k] = sys.theta(k, x[k]);
  arma::mat S = sys.B;
  arma::mat& V = st.V;
  V.set_size(q, q);
  for (int k = 0; k < q; ++k) {
    V.row(k) = sys.R.row(k) * S;
    if (th[k] > 0.0) S.row(k) -= th[k] * V.row(k);
  }
  arma::mat N;
  if (!arma::inv(N, arma::eye(q, q) - S) || !N.is_finite()) return st;
  const arma::mat NC = N * sys.C;
  st.nc = NC.diag();
  for (int j = 0; j < q; ++j) st.tau[j] = 2.0 * st.nc[j] / sys.C(j, j) - 1.0;
  st.ok = true;
  if (!trials) return st;
  arma::mat P = arma::eye(q, q);
  arma::mat& U = st.U;
  U.set_size(q, q);
  for (int k = q - 1; k >= 0; --k) {
    U.col(k) = P.col(k);
    if (th[k] > 0.0) P -= th[k] * U.col(k) * sys.R.row(k);
  }
  st.NU = N * U;
  st.VN = V * N;
  st.VNC = st.VN * sys.C;
  st.den.set_size(q);
  for (int k = 0; k < q; ++k) st.den[k] = arma::dot(st.VN.row(k), U.col(k));
  return st;
}

// The state after coordinate k's theta rises by d: A loses d u_k v_k', the
// later coordinates' v and the earlier ones' u each change by a multiple of
// v_k and u_k, and (I - A)^-1 by a rank-one term, so every product the
// trials read is a rank-one update of itself.
void group_advance(const GroupSystem& sys, GroupState& st, int k, double d) {
  const int q = sys.q;
  const double s = d / (1.0 + d * st.den[k]);
  const arma::vec nu = st.NU.col(k);
  const arma::rowvec vn = st.VN.row(k), vnc = st.VNC.row(k);
  // c_j = [r_k' M_{k-1} ... M_{j+1}]_j for j < k; w_j = r_j' M_{j-1} ... M_{k+1} e_k for j > k.
  arma::rowvec c(q, arma::fill::zeros);
  arma::rowvec z = sys.R.row(k);
  for (int j = k - 1; j >= 0; --j) {
    c[j] = z[j];
    if (st.th[j] > 0.0) z -= st.th[j] * z[j] * sys.R.row(j);
  }
  arma::vec w(q, arma::fill::zeros), y(q, arma::fill::zeros);
  y[k] = 1.0;
  for (int j = k + 1; j < q; ++j) {
    w[j] = arma::dot(sys.R.row(j), y);
    if (st.th[j] > 0.0) y[j] -= st.th[j] * w[j];
  }
  const arma::rowvec vnU = vn * st.U;
  const arma::vec Vnu = st.V * nu;
  for (int j = 0; j < q; ++j) st.nc[j] -= s * nu[j] * vnc[j];
  st.NU -= nu * (d * c + s * (vnU - d * st.den[k] * c));
  const arma::vec left = s * Vnu + d * (1.0 - s * st.den[k]) * w;
  st.VN -= left * vn;
  st.VNC -= left * vnc;
  st.U -= d * st.U.col(k) * c;
  st.V -= d * w * st.V.row(k);
  st.th[k] += d;
  for (int j = 0; j < q; ++j) {
    st.den[j] = arma::dot(st.VN.row(j), st.U.col(j));
    st.tau[j] = 2.0 * st.nc[j] / sys.C(j, j) - 1.0;
  }
  st.ok = st.nc.is_finite() && st.den.is_finite();
}

GroupSystem group_system(const arma::mat& C, const arma::mat& Q, const arma::mat& Ar,
                         const Rcpp::NumericVector& fac, const Rcpp::NumericVector& g,
                         const Rcpp::NumericVector& a) {
  GroupSystem sys;
  sys.q = C.n_rows;
  sys.C = C;
  sys.R = Q;
  for (int k = 0; k < sys.q; ++k) sys.R.row(k) /= Q(k, k);
  sys.B = Ar;
  for (int k = 0; k < sys.q; ++k) sys.B.row(k) *= fac[k];
  sys.g.assign(g.begin(), g.end());
  sys.a.resize(sys.q);
  for (int k = 0; k < sys.q; ++k)
    sys.a[k] = std::isfinite(a[k]) ? a[k] : sys.g[k];
  return sys;
}

bool group_dims(const arma::mat& C, const arma::mat& Q, const arma::mat& Ar, int n_fac,
                int n_g, int n_a) {
  const int q = C.n_rows;
  return q > 0 && static_cast<int>(C.n_cols) == q && static_cast<int>(Q.n_rows) == q &&
    static_cast<int>(Q.n_cols) == q && static_cast<int>(Ar.n_rows) == q &&
    static_cast<int>(Ar.n_cols) == q && n_fac == q && n_g == q && n_a == q;
}

}  // namespace

// Each coordinate's equivalent lag-one rho, (tau - 1) / (tau + 1), under the
// plans in the columns of X (.emc_alloc_group_rho). A plan under which the
// system has no stationary answer reads as rho one.
// [[Rcpp::export]]
SEXP emc_alloc_group_rho_cpp(const arma::mat& C, const arma::mat& Q, const arma::mat& Ar,
                             const Rcpp::NumericVector& fac, const Rcpp::NumericVector& g,
                             const Rcpp::NumericVector& a, const Rcpp::NumericMatrix& X) {
  if (!group_dims(C, Q, Ar, fac.size(), g.size(), a.size()) ||
      X.nrow() != static_cast<int>(C.n_rows))
    return R_NilValue;
  const GroupSystem sys = group_system(C, Q, Ar, fac, g, a);
  const int q = sys.q, m = X.ncol();
  Rcpp::NumericMatrix out(q, m);
  std::vector<double> x(q);
  for (int c = 0; c < m; ++c) {
    for (int k = 0; k < q; ++k) x[k] = X(k, c);
    const GroupState st = group_state(sys, x, false);
    for (int k = 0; k < q; ++k) out(k, c) = tau_to_rho(st.tau[k]);
  }
  return out;
}

// The cheapest rounds that bring every movable coordinate's rho (floored by
// its between-chain share) to each target rs (ascending). One path from the
// floors: while a movable coordinate is above the loosest unmet target, the
// coordinate whose next step (an eighth of a round below one round, a
// quarter below two, a half beyond) takes most of the movable coordinates'
// excess autocorrelation time over that target per unit cost takes it; a
// step's effect on every coordinate is a rank-one update of the system.
// Returns x (coordinates by targets) and whether each target was reached; a
// target no step brings nearer is out of reach, as are those below it. The
// system is rebuilt in full every `rebuild` steps.
// [[Rcpp::export]]
SEXP emc_alloc_group_need_cpp(const arma::mat& C, const arma::mat& Q, const arma::mat& Ar,
                              const Rcpp::NumericVector& fac, const Rcpp::NumericVector& g,
                              const Rcpp::NumericVector& a, const Rcpp::NumericVector& floor,
                              const Rcpp::NumericVector& cap, const Rcpp::NumericVector& cost,
                              const Rcpp::NumericVector& b_share,
                              const Rcpp::LogicalVector& movable,
                              const Rcpp::NumericVector& rs, int rebuild = 16) {
  if (!group_dims(C, Q, Ar, fac.size(), g.size(), a.size())) return R_NilValue;
  const GroupSystem sys = group_system(C, Q, Ar, fac, g, a);
  const int q = sys.q, n_t = rs.size();
  if (floor.size() != q || cap.size() != q || cost.size() != q || b_share.size() != q ||
      movable.size() != q)
    return R_NilValue;
  std::vector<double> x(q);
  for (int k = 0; k < q; ++k) x[k] = std::min<double>(floor[k], cap[k]);
  auto floored_tau = [&](int j, double tau) {
    const double r = b_share[j] + (1.0 - b_share[j]) * tau_to_rho(tau);
    const double rr = std::min(r, 1.0 - 1e-12);
    return (1.0 + rr) / (1.0 - rr);
  };
  Rcpp::NumericMatrix X(q, n_t);
  Rcpp::LogicalVector reached(n_t);
  GroupState st = group_state(sys, x, true);
  std::vector<double> ft(q);
  auto refresh = [&]() {
    for (int j = 0; j < q; ++j) ft[j] = floored_tau(j, st.tau[j]);
  };
  refresh();
  bool stuck = !st.ok;
  const int max_steps = 64 * q + 64;
  int steps = 0;
  for (int t = n_t - 1; t >= 0; --t) {
    const double tt = (1.0 + std::min(rs[t], 1.0 - 1e-12)) / (1.0 - std::min(rs[t], 1.0 - 1e-12));
    while (!stuck) {
      double excess = 0.0;
      for (int j = 0; j < q; ++j)
        if (movable[j] && ft[j] > tt * (1.0 + 1e-12)) excess += ft[j] - tt;
      if (!(excess > 0.0)) break;
      if (++steps > max_steps) { stuck = true; break; }
      int best = -1;
      double best_gain = 0.0, best_step = 0.0, best_d = 0.0;
      for (int k = 0; k < q; ++k) {
        if (!sys.moves(k) || !(x[k] < cap[k] - 1e-12)) continue;
        const double unit = x[k] < 1.0 - 1e-12 ? 0.125 : (x[k] < 2.0 - 1e-12 ? 0.25 : 0.5);
        const double step = std::min(unit, cap[k] - x[k]);
        const double d = sys.theta(k, x[k] + step) - sys.theta(k, x[k]);
        if (!(d > 0.0)) continue;
        const double scale = d / (1.0 + d * st.den[k]);
        if (!std::isfinite(scale)) continue;
        double after = 0.0;
        for (int j = 0; j < q; ++j) {
          if (!movable[j]) continue;
          const double tau = 2.0 * (st.nc[j] - scale * st.NU(j, k) * st.VNC(k, j)) / sys.C(j, j) - 1.0;
          const double f = floored_tau(j, tau);
          if (f > tt) after += f - tt;
        }
        const double price = step * (cost[k] > 0.0 ? cost[k] : 1e-12);
        const double gain = (excess - after) / price;
        if (gain > best_gain) { best_gain = gain; best = k; best_step = step; best_d = d; }
      }
      if (best < 0 || !(best_gain * best_step > 1e-12 * excess)) { stuck = true; break; }
      x[best] += best_step;
      if (rebuild > 1 && steps % rebuild != 0) group_advance(sys, st, best, best_d);
      else st = group_state(sys, x, true);
      if (!st.ok) { stuck = true; break; }
      refresh();
    }
    bool met = !stuck;
    if (stuck) {
      met = st.ok;
      for (int j = 0; j < q && met; ++j)
        if (movable[j] && ft[j] > tt * (1.0 + 1e-12)) met = false;
    }
    reached[t] = met;
    for (int k = 0; k < q; ++k) X(k, t) = x[k];
  }
  return Rcpp::List::create(Rcpp::Named("x") = X, Rcpp::Named("reached") = reached);
}
