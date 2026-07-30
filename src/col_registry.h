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

namespace emc2col {

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

// R/model_LNR.R — LNR
namespace lnr {
  enum : int { m = 0, s, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"m", "s", "t0"};
    return {n, N_REQ, "LNR"};
  }
}

// R/model_LNR.R — RGAMMA
namespace rgamma {
  enum : int { lambda = 0, shape, shift, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"lambda", "shape", "shift"};
    return {n, N_REQ, "RGAMMA"};
  }
}

// R/model_LNR.R — PCOUNTER (Poisson counter race).  Ratcliff & Smith (2004,
// Appendix, Eq. A10a/A10b): counter i accrues unit counts as a Poisson process
// with rate alpha until it reaches criterion K, so its first-passage time is
// Erlang(K, alpha) shifted by t0.  Same kernel positions as rgamma
// (rate, shape, shift) under counter-model names.
namespace pcounter {
  enum : int { alpha = 0, K, t0, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"alpha", "K", "t0"};
    return {n, N_REQ, "PCOUNTER"};
  }
}

// R/model_LNR.R — REXG
namespace rexg {
  enum : int { mu = 0, sigma, tau, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"mu", "sigma", "tau"};
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

// R/model_DDM.R — DDM
namespace ddm {
  enum : int { v = 0, a, sv, t0, st0, s, Z, SZ, N_REQ };
  inline ColSpec spec() {
    static const char* n[] = {"v", "a", "sv", "t0", "st0", "s", "Z", "SZ"};
    return {n, N_REQ, "DDM"};
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
