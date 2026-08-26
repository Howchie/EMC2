#ifndef EMC2_COL_REGISTRY_H
#define EMC2_COL_REGISTRY_H

#include <Rcpp.h>
#include <string>

// Single source of truth for the parameter-column order the C++ batch kernels
// expect. Each namespace's enum order MUST equal the R-side p_types order
// declared in the referenced model constructor; validate_col_prefix() enforces
// the contract once per likelihood call, turning a silent p_types reordering
// into an immediate error.
//
// Columns listed after N_REQ are optional: they exist only for some model
// variants, and kernels must gate every dereference of them on the
// corresponding context flag (kill_active, kill_shape == 3, pc_col >= 0, ...).

// `eta` (the operational-time warp of Math/ballistic-time.md) is deliberately
// NOT in any enum here: it is resolved BY NAME in configure_time_warp_context()
// (src/time_warp.cpp), like `rho`, so that adding it to a model's p_types needs
// no new ColSpec and cannot shift any positional index.
namespace emc2col {

// Column contracts use separate enum types to keep each model's layout
// explicit.  Convert only at the selection boundary when a kernel supports
// both the legacy and split layouts; this avoids -Wenum-compare diagnostics.
template <typename T, typename U>
constexpr int select_index(bool condition, T when_true, U when_false) noexcept {
  return condition ? static_cast<int>(when_true) : static_cast<int>(when_false);
}

struct ColSpec {
  const char* const* names;  // required column names, in kernel order
  int n_required;
  const char* label;         // model family, for error messages
};

// R/model_LBA.R — LBA / LBAIO (threshold b = B + A)
namespace lba {
  enum : int { v = 0, sv, B, A, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0"};
    return {n, N_REQ, "LBA"};
  }
}

// R/model_RDM.R — RDM
namespace rdm {
  enum : int { v = 0, B, A, t0, s, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s"};
    return {n, N_REQ, "RDM"};
  }
}

// R/model_RLF.R — RLF, the racing symmetric alpha-stable Lévy flight.
// The state parameters use the RDM convention b = B + A and X(0) ~ U(0, A);
// alpha controls the stable-noise tail and s is scaled out of the PDE key.
namespace rlf {
  enum : int { v = 0, B, A, t0, s, alpha, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s", "alpha"};
    return {n, N_REQ, "RLF"};
  }
}

// R/model_LNR.R — LNR
namespace lnr {
  enum : int { m = 0, s, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"m", "s", "t0"};
    return {n, N_REQ, "LNR"};
  }
}

// R/model_LNR.R — PCOUNTER (gamma-mixed, self-exciting Poisson counter race).
// The leading columns are the natural-scale parameters documented by
// PCOUNTER() in R/model_LNR.R.
namespace pcounter {
  enum : int { nu = 0, sv, gamma, k, omega, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"nu", "sv", "gamma", "k", "omega", "t0"};
    return {n, N_REQ, "PCOUNTER"};
  }
}

// R/model_LNR.R — REXG
namespace rexg {
  enum : int { mu = 0, sigma, tau, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "tau", "t0"};
    return {n, N_REQ, "REXG"};
  }
}

// R/model_LBA.R — BAwL (leaky ballistic accumulator with kill/guess clocks).
// mG/mK dereferences are gated on ctx->kill_active, omega on kill_shape == 3.
namespace bawl {
  enum : int { v = 0, sv, B, A, t0, k, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwL"};
  }
}
// Lognormal launch strength: same POSITIONS as bawl (mu occupies v's slot and
// sigma occupies sv's), so every kernel indexes through the bawl enum and only
// the names validate_col_prefix() insists on differ.  Mirrors bawd/bawd_logn.
namespace bawl_logn {
  enum : int { mu = 0, sigma, B, A, t0, k, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwL_LOGN"};
  }
}
namespace bawl_weib {
  enum : int { shape = 0, scale, B, A, t0, k, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwL_WEIB"};
  }
}
namespace bawlsplit {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwL_LOGN_SPLIT"};
  }
}

// R/model_BAwD.R — BAwD (ballistic accumulator with drive decay).  The two
// launch distributions differ only in the names of the first two columns, so
// the POSITIONS are deliberately identical and every kernel indexes them
// through the same enum; the spec selected in resolve_race_model_adapter() only
// changes what validate_col_prefix() insists the R p_types are called.
namespace bawd {
  enum : int { v = 0, sv, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k", "ell"};
    return {n, N_REQ, "BAwD"};
  }
}
namespace bawd_logn {
  // Same positions as bawd: mu occupies v's slot and sigma occupies sv's.
  enum : int { mu = 0, sigma, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k", "ell"};
    return {n, N_REQ, "BAwD_LOGN"};
  }
}
namespace bawd_weib {
  enum : int { shape = 0, scale, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k", "ell"};
    return {n, N_REQ, "BAwD_WEIB"};
  }
}
namespace bawdsplit {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, ell, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k", "ell"};
    return {n, N_REQ, "BAwD_LOGN_SPLIT"};
  }
}

// R/model_BAwF.R — BAwF (global fading of decision-relevant evidence).
// X(u) = h_rho(u) [z + V u]: there is no clearance parameter, so the required
// columns are BAwL's exactly.  It still needs its own namespaces, because
// sharing bawl::spec() would let a BAwF c_name validate against BAwL's kernel
// contract (which carries the optional mG/mK/omega kill columns BAwF has no
// meaning for).
namespace bawf {
  enum : int { v = 0, sv, B, A, t0, k, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwF"};
  }
}
namespace bawf_logn {
  // Same positions as bawf: mu occupies v's slot and sigma occupies sv's.
  enum : int { mu = 0, sigma, B, A, t0, k, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwF_LOGN"};
  }
}
namespace bawf_weib {
  enum : int { shape = 0, scale, B, A, t0, k, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwF_WEIB"};
  }
}
namespace bawfsplit {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k"};
    return {n, N_REQ, "BAwF_LOGN_SPLIT"};
  }
}

// R/model_BAwR.R — BAwR (ramping clearance, dX/du = V - kappa u^p).  The "R"
// is for the RAMP: the clearance rate kappa u^p grows with elapsed time.
// Deliberately not sharing bawf::spec(): the decay here needs two columns
// (coefficient and exponent) rather than one rate, and `kappa` has different
// units from BAwF's `k`.
namespace bawr {
  enum : int { v = 0, sv, B, A, t0, kappa, p, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "kappa", "p"};
    return {n, N_REQ, "BAwR"};
  }
}
namespace bawr_logn {
  // Same positions as bawr: mu occupies v's slot and sigma occupies sv's.
  enum : int { mu = 0, sigma, B, A, t0, kappa, p, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "kappa", "p"};
    return {n, N_REQ, "BAwR_LOGN"};
  }
}
namespace bawr_weib {
  enum : int { shape = 0, scale, B, A, t0, kappa, p, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "kappa", "p"};
    return {n, N_REQ, "BAwR_WEIB"};
  }
}
namespace bawrsplit {
  enum : int { mu = 0, sigma, delta, B, A, t0, kappa, p, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "kappa", "p"};
    return {n, N_REQ, "BAwR_LOGN_SPLIT"};
  }
}

// R/model_BTAwL.R — BTAwL pure transient member. As with the other ballistic
// launch models, the first two positions are either (v, sv) or (mu, sigma).
namespace btawl_transient {
  enum : int { v = 0, sv, B, A, t0, k, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k", "tau"};
    return {n, N_REQ, "BTAwL_RATE"};
  }
}
namespace btawl_transient_logn {
  enum : int { mu = 0, sigma, B, A, t0, k, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k", "tau"};
    return {n, N_REQ, "BTAwL_LOGN_RATE"};
  }
}
namespace btawl_transient_weib {
  enum : int { shape = 0, scale, B, A, t0, k, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k", "tau"};
    return {n, N_REQ, "BTAwL_WEIB_RATE"};
  }
}
namespace btawlsplit_transient {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k", "tau"};
    return {n, N_REQ, "BTAwL_LOGN_SPLIT_RATE"};
  }
}
namespace btawl_transient_split = btawlsplit_transient;

namespace btawl_local_race {
  // Full BTAwL local-race contract. pi remains the final position.
  enum : int { v = 0, sv, B, A, t0, k, tau_s, tau_t, pi, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k", "tau_s", "tau_t", "pi"};
    return {n, N_REQ, "BTAwL_RATE"};
  }
}
namespace btawl_local_race_logn {
  // Full BTAwL local-race contract. pi remains the final position.
  enum : int { mu = 0, sigma, B, A, t0, k, tau_s, tau_t, pi, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k", "tau_s", "tau_t", "pi"};
    return {n, N_REQ, "BTAwL_LOGN_RATE"};
  }
}
namespace btawl_local_race_weib {
  enum : int { shape = 0, scale, B, A, t0, k, tau_s, tau_t, pi, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k", "tau_s", "tau_t", "pi"};
    return {n, N_REQ, "BTAwL_WEIB_RATE"};
  }
}
namespace btawlsplit_local_race {
  // Full BTAwL local-race contract with split lognormal launch.
  enum : int { mu = 0, sigma, delta, B, A, t0, k, tau_s, tau_t, pi, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k", "tau_s", "tau_t", "pi"};
    return {n, N_REQ, "BTAwL_LOGN_SPLIT_RATE"};
  }
}
namespace btawl_local_race_split = btawlsplit_local_race;

// Pure sustained BTAwL wrapper.  The full BTAwL local race uses the nine-column
// contract above; this seven-column contract is only for the pi = 1 wrapper.
namespace btawl_sustained {
  enum : int { v = 0, sv, B, A, t0, k, tau_s, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k", "tau_s"};
    return {n, N_REQ, "BTAwL_SUSTAINED"};
  }
}
namespace btawl_sustained_logn {
  enum : int { mu = 0, sigma, B, A, t0, k, tau_s, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k", "tau_s"};
    return {n, N_REQ, "BTAwL_SUSTAINED_LOGN"};
  }
}
namespace btawl_sustained_weib {
  enum : int { shape = 0, scale, B, A, t0, k, tau_s, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k", "tau_s"};
    return {n, N_REQ, "BTAwL_SUSTAINED_WEIB"};
  }
}
namespace btawlsplit_sustained {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, tau_s, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k", "tau_s"};
    return {n, N_REQ, "BTAwL_SUSTAINED_LOGN_SPLIT"};
  }
}
namespace btawl_sustained_split = btawlsplit_sustained;

// R/model_BAwD.R — BAwDp (proportional-clearance drive clock).  As with BAwL,
// the launch pair occupies the first two positions; only the names differ for
// the lognormal launch.  `lambda` is the dimensionless clearance fraction.
namespace bawdp {
  enum : int { v = 0, sv, B, A, t0, k, lambda, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "sv", "B", "A", "t0", "k", "lambda"};
    return {n, N_REQ, "BAwDp"};
  }
}
namespace bawdp_logn {
  enum : int { mu = 0, sigma, B, A, t0, k, lambda, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "B", "A", "t0", "k", "lambda"};
    return {n, N_REQ, "BAwDp_LOGN"};
  }
}
namespace bawdp_weib {
  enum : int { shape = 0, scale, B, A, t0, k, lambda, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"shape", "scale", "B", "A", "t0", "k", "lambda"};
    return {n, N_REQ, "BAwDp_WEIB"};
  }
}
namespace bawdpsplit {
  enum : int { mu = 0, sigma, delta, B, A, t0, k, lambda, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "delta", "B", "A", "t0", "k", "lambda"};
    return {n, N_REQ, "BAwDp_LOGN_SPLIT"};
  }
}

// R/model_FRQ.R — FRQ (finite reservoir quorum).  alpha/beta are the Beta
// shapes of the latent quorum U ~ Beta(alpha, beta), i.e. the quorum size K
// and the residual redundancy N - K + 1.  The latter is called `beta` rather
// than `R` because `R` is a reserved data column name in EMC2
// (R/design.R:251).  h is the eventual completion probability and tau the
// conditional median decision time; the kernel inverts both to (p, lambda).
// delta is between-trial threshold variability (half-width of the uniform
// log-odds shift of the quorum percentile) and is REQUIRED rather than
// optional even though it defaults to zero: an optional trailing column is
// resolved positionally, so a design lacking it would silently read whatever
// parameter happened to land in slot 5 instead of erroring.
namespace frq {
  enum : int { alpha = 0, beta, h, tau, t0, delta, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"alpha", "beta", "h", "tau", "t0", "delta"};
    return {n, N_REQ, "FRQ"};
  }
}

// R/model_RDM.R — RDMGBM (geometric Brownian motion race with timers).
// mG/mK gated on ctx->kill_active, omega on kill_shape == 3.
namespace rdmgbm {
  enum : int { v = 0, B, A, t0, s, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s"};
    return {n, N_REQ, "RDMGBM"};
  }
}

// R/model_RDM.R — RDMSWTN (shifted Wald with timers, normal drift variability).
// mG/mK gated on ctx->kill_active, omega on kill_shape == 3.
namespace rdmswtn {
  enum : int { v = 0, B, A, t0, s, sv, N_REQ, mG = N_REQ, mK, omega };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s", "sv"};
    return {n, N_REQ, "RDMSWTN"};
  }
}

// R/model_RDM.R — RDMSWTN_TT (RDMSWTN under a finite linear exhaustion clock).
namespace rdmswtn_tt {
  enum : int { v = 0, B, A, t0, s, sv, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s", "sv", "tau"};
    return {n, N_REQ, "RDMSWTN_TT"};
  }
}

// R/model_ROU.R — ROU, the racing Ornstein-Uhlenbeck (leaky accumulator).
// Identical to RDM apart from the leak k, so an RDM design converts by adding
// k~1; k = 0 is the Wiener race, reached through the PDE rather than by
// dispatching to the Wald kernels.  Binf/tau/pw exist only for
// boundary_collapse != "fixed" and must be gated on SolveCache::bnd_kind before
// they are dereferenced; pw only for the Weibull form.
namespace rou {
  enum : int { v = 0, k, B, A, t0, s, N_REQ, Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"v", "k", "B", "A", "t0", "s"};
    return {n, N_REQ, "ROU"};
  }
}

// ROU(parameterization = "curvature"): (tstar, k, s) replace (v, k, s), with
// B, A, t0 and the optional collapse columns unchanged.  The optional columns
// land at the same indices as the rate layout because N_REQ is also 6.
namespace rou_curv {
  enum : int { tstar = 0, k, s, B, A, t0, N_REQ, Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"tstar", "k", "s", "B", "A", "t0"};
    return {n, N_REQ, "ROUCURV"};
  }
}

// ROU(parameterization = "equilibrium"): (tk, theta, chi) replace (v, k, s).  tk is
// the LEAK time constant 1/k; the boundary-collapse time constant keeps the
// name tau, which is why the leak one is not called tau here.
namespace rou_eq {
  enum : int { tk = 0, theta, chi, B, A, t0, N_REQ, Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"tk", "theta", "chi", "B", "A", "t0"};
    return {n, N_REQ, "ROUEQ"};
  }
}

// R/model_ROUp.R — ROUp, racing Ornstein-Uhlenbeck with Smith (1995) pulse drift.
namespace roup {
  enum : int { v_S = 0, v_T, tau_S, tau_T, k, B, A, t0, s, N_REQ, Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"v_S", "v_T", "tau_S", "tau_T", "k", "B", "A", "t0", "s"};
    return {n, N_REQ, "ROUp"};
  }
}

// R/model_ROUp.R parameterization = "area": E_T replaces v_T while the
// sustained channel and all boundary columns retain the rate-chart positions.
namespace roup_area {
  enum : int { v_S = 0, E_T, tau_S, tau_T, k, B, A, t0, s, N_REQ,
               Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"v_S", "E_T", "tau_S", "tau_T", "k", "B", "A", "t0", "s"};
    return {n, N_REQ, "ROUpAREA"};
  }
}

// R/model_GOM.R — Gompertz growth-process race.  The process is solved after
// Y = log(X), but alpha, beta and K remain on the physical Gompertz scale and
// A is the physical start-point range [1, 1 + A].
namespace gompertz {
  enum : int { alpha = 0, beta, K, B, A, t0, N_REQ,
               Binf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"alpha", "beta", "K", "B", "A", "t0"};
    return {n, N_REQ, "GOM"};
  }
}

// R/model_DDM.R — DDM
namespace ddm {
  enum : int { v = 0, a, sv, t0, st0, s, Z, SZ, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "a", "sv", "t0", "st0", "s", "Z", "SZ"};
    return {n, N_REQ, "DDM"};
  }
}

// R/model_BOU.R — BOU, the bounded Ornstein-Uhlenbeck (Smith & Ratcliff 2004).
// The DDM's columns plus the leak beta, in the same order, so a DDM design
// converts by adding beta~1 and beta = 0 recovers the DDM exactly -- reached
// through the Fokker-Planck solver rather than the Wiener series.
// aInf/tau/pw exist only for boundary_collapse != "fixed" and must be gated on
// ContextForDDMModels::bnd_kind before they are dereferenced, exactly as rou's
// Binf/tau/pw are; pw only for the Weibull form.  aInf is the ASYMPTOTIC
// SEPARATION, in the same units as a, and the two barriers close on the midpoint
// symmetrically.
namespace bou {
  enum : int { v = 0, a, sv, t0, st0, s, Z, SZ, beta, N_REQ,
               aInf = N_REQ, tau, pw };
  inline ColSpec spec() {
    static const char* n[] = {"v", "a", "sv", "t0", "st0", "s", "Z", "SZ",
                              "beta"};
    return {n, N_REQ, "BOU"};
  }
}

// R/model_SS.R — stop-signal truncated ex-Gaussian (SSEXG)
namespace ss_texg {
  enum : int { mu = 0, sigma, tau, muS, sigmaS, tauS, tf, gf, exg_lb, exgS_lb, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "tau", "muS", "sigmaS", "tauS",
                              "tf", "gf", "exg_lb", "exgS_lb"};
    return {n, N_REQ, "SSEXG"};
  }
}

// R/model_SS.R — stop-signal RDM go / ex-Gaussian stop (SSRDEX)
namespace ss_rdex {
  enum : int { v = 0, B, A, t0, s, muS, sigmaS, tauS, tf, gf, exgS_lb, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "B", "A", "t0", "s", "muS", "sigmaS", "tauS",
                              "tf", "gf", "exgS_lb"};
    return {n, N_REQ, "SSRDEX"};
  }
}

// Checks that col_names starts with the spec's required columns in order.
// Kernels index the leading columns positionally, so a mismatch here means the
// kernels would silently read the wrong parameters — stop loudly instead.
inline void validate_col_prefix(const Rcpp::CharacterVector& col_names,
                                const ColSpec& spec) {
  if (static_cast<int>(col_names.size()) < spec.n_required) {
    Rcpp::stop("%s kernels require %d leading parameter columns "
               "(see src/col_registry.h) but only %d were supplied.",
               spec.label, spec.n_required, static_cast<int>(col_names.size()));
  }
  for (int j = 0; j < spec.n_required; ++j) {
    const std::string got = Rcpp::as<std::string>(col_names[j]);
    if (got != spec.names[j]) {
      Rcpp::stop("%s kernels expect parameter column %d to be '%s' but got "
                 "'%s'; the R p_types order must match src/col_registry.h.",
                 spec.label, j + 1, spec.names[j], got.c_str());
    }
  }
}

} // namespace emc2col

#endif // EMC2_COL_REGISTRY_H
