#include "gh_quad.h"
#include "utility_functions.h"
#include "race_dispatch.h"
#include "utils.h"
#include "model_lnr.h"
#include "model_EXG.h"
#include "col_registry.h"
#include "fpe_race.h"
#include "fpe_models.h"
#include "model_FRQ.h"
#include "model_PCOUNTER.h"
#include <cmath>
#include <memory>
#include <string>

RaceModelAdapter resolve_race_model_adapter(const std::string& type_std,
                                                          const std::string& caller) {
  RaceModelAdapter out;
  out.ctx.min_lik_for_pdf = std::exp(std::log(1e-10));
  out.ctx.use_posdrift = true;
  out.ctx.gng = false;
  out.ctx.kill_shape = (type_std.find("_EMIX") != std::string::npos) ? 3 :
                       ((type_std.find("_E2") != std::string::npos) ? 2 : 1);

  // Erlang process flag resolution
  out.ctx.is_local_guess = (type_std.find("_LOCAL_GUESS") != std::string::npos);
  out.ctx.is_global_kill = (type_std.find("_GLOBAL_KILL") != std::string::npos);
  out.ctx.is_local_kill  = (type_std.find("_LOCAL_KILL")  != std::string::npos);
  out.ctx.is_local_kill_guess = (type_std.find("_LOCAL_KILL_GUESS") != std::string::npos);

  // Backward compatibility: bare _GLOBAL maps to global_kill
  if (!out.ctx.is_local_guess && !out.ctx.is_global_kill && !out.ctx.is_local_kill && !out.ctx.is_local_kill_guess) {
    if (type_std.find("_GLOBAL") != std::string::npos) out.ctx.is_global_kill = true;
  }

  // A bare BAwL/RDM-with-timers type has timer columns in its parameter
  // contract, but no active Erlang process.  Keep that distinction in the
  // scalar and raw paths; treating the default mG/mK values as active clocks
  // changes the likelihood and needlessly evaluates the clock mixture for
  // every proposal.  The suffixed variants enable the process below.
  out.ctx.kill_active = out.ctx.is_local_guess || out.ctx.is_global_kill ||
                        out.ctx.is_local_kill || out.ctx.is_local_kill_guess;

  out.ctx.apply_lk_to_racers = !out.ctx.is_global_kill;

  if (out.ctx.is_global_kill) out.ctx.defective_upper_tail = true;

  if (type_std.find("RLF") != std::string::npos) {
    out.pdf1_ptr       = &drlf_scalar;
    out.cdf1_ptr       = &prlf_scalar;
    out.model_dfun_raw = &drlf_raw;
    out.model_pfun_raw = &prlf_raw;
    out.logS_at_t_ptr  = &rlf_logS_at_t;
    out.col_spec       = emc2col::rlf::spec();
    out.ctx.t0_index   = emc2col::rlf::t0;
    out.ctx.rlf_cache = std::make_shared<rlf::SolveCache>();
    rlf_configure_grid(out.ctx.rlf_cache->grid);
  } else if (type_std.find("GOM") != std::string::npos ||
             type_std.find("GOMP") != std::string::npos) {
    // Gompertz growth is an OU after Y = log(X).  The adapter keeps the
    // physical alpha/beta/K columns and gomp_key() performs the reduction; its
    // log_state flag makes the FPE seed integrate A uniformly on X rather than
    // incorrectly treating the transformed start range as uniform in Y.
    out.pdf1_ptr       = &dgomp_scalar;
    out.cdf1_ptr       = &pgomp_scalar;
    out.model_dfun_raw = &dgomp_raw;
    out.model_pfun_raw = &pgomp_raw;
    out.logS_at_t_ptr  = &gomp_logS_at_t;
    out.col_spec       = emc2col::gompertz::spec();
    out.ctx.t0_index   = emc2col::gompertz::t0;
    out.ctx.defective_upper_tail = false;
    out.ctx.fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_cache(*out.ctx.fpe_cache);
    if (type_std.find("_BWEIB") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_WEIBULL;
    else if (type_std.find("_BEXP") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_EXPONENTIAL;
    else if (type_std.find("_BLIN_MULT") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
    else if (type_std.find("_BLIN_ADD") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_ADDITIVE;
  } else if (type_std.find("ROUp") != std::string::npos) {
    out.pdf1_ptr       = &droup_scalar;
    out.cdf1_ptr       = &proup_scalar;
    out.model_dfun_raw = &droup_raw;
    out.model_pfun_raw = &proup_raw;
    out.logS_at_t_ptr  = &roup_logS_at_t;
    out.ctx.defective_upper_tail = out.ctx.is_global_kill;
    out.ctx.fpe_cache  = std::make_shared<fperace::SolveCache>();
    roup_configure_cache(*out.ctx.fpe_cache);
    // Local-race pooling suffix (ROUp(pooling = "local_race")): each marginal
    // kernel races its sustained-only and transient-only subraces
    // independently; the default coactive pooling shares one race.  The cache
    // flag carries the choice from the c_name to the FPE worker.
    out.ctx.fpe_cache->roup_local =
      (type_std.find("_LOCAL_RACE") != std::string::npos);
    if (type_std.find("ROUpAREA") != std::string::npos) {
      out.ctx.fpe_cache->par_kind = fperace::ROUP_PAR_AREA;
      out.col_spec       = emc2col::roup_area::spec();
      out.ctx.t0_index   = emc2col::roup_area::t0;
    } else {
      out.ctx.fpe_cache->par_kind = fperace::ROUP_PAR_RATE;
      out.col_spec       = emc2col::roup::spec();
      out.ctx.t0_index   = emc2col::roup::t0;
    }
    if (type_std.find("_BWEIB") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_WEIBULL;
    else if (type_std.find("_BEXP") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_EXPONENTIAL;
    else if (type_std.find("_BLIN_MULT") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
    else if (type_std.find("_BLIN_ADD") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_ADDITIVE;
  } else if (type_std.find("ROU") != std::string::npos) {
    // Ordered FIRST deliberately.  Dispatch here is by substring, so a key that
    // is a substring of a later one must be tested first; "ROU" collides with
    // nothing today, and testing it first is what keeps a future addition from
    // silently capturing it.
    out.pdf1_ptr       = &drou_scalar;
    out.cdf1_ptr       = &prou_scalar;
    out.model_dfun_raw = &drou_raw;
    out.model_pfun_raw = &prou_raw;
    out.logS_at_t_ptr  = &rou_logS_at_t;
    // With non-zero diffusion and a finite upper boundary, the OU hits the
    // boundary almost surely.  A subthreshold equilibrium creates long finite
    // survival and an approximately constant late hazard, not a point mass at
    // infinity.  Preserve a genuine global-kill defect if a suffixed variant
    // requests that clock.
    out.ctx.defective_upper_tail = out.ctx.is_global_kill;
    out.ctx.fpe_cache = std::make_shared<fperace::SolveCache>();
    rou_configure_cache(*out.ctx.fpe_cache);
    // Parameterisation, selected by ROU(parameterization=).  It changes only
    // which columns the kernels read and how they map onto (v, k, s); the
    // solve, the cache and the boundary handling are untouched.  Tested before
    // the boundary suffix because the two compose: ROUCURV_BEXP is both.
    if (type_std.find("ROUCURV") != std::string::npos) {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_CURVATURE;
      out.col_spec     = emc2col::rou_curv::spec();
      out.ctx.t0_index = emc2col::rou_curv::t0;
    } else if (type_std.find("ROUEQ") != std::string::npos) {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_EQUILIBRIUM;
      out.col_spec     = emc2col::rou_eq::spec();
      out.ctx.t0_index = emc2col::rou_eq::t0;
    } else {
      out.ctx.fpe_cache->par_kind = fperace::ROU_PAR_RATE;
      out.col_spec     = emc2col::rou::spec();
      out.ctx.t0_index = emc2col::rou::t0;
    }
    // Collapsing-bound variants, selected by ROU(boundary_collapse=).  The
    // suffix carries the FORM only; the shape parameters are ordinary optional
    // columns (Binf/tau/pw) that the design system estimates like any other.
    // Matched on the bare suffix rather than on "ROU_B*" so that the collapse
    // forms are available under every parameterisation; we are already inside
    // the ROU branch, so there is nothing else these can capture.
    if (type_std.find("_BWEIB") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_WEIBULL;
    else if (type_std.find("_BEXP") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_EXPONENTIAL;
    else if (type_std.find("_BLIN_MULT") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_MULTIPLICATIVE;
    else if (type_std.find("_BLIN_ADD") != std::string::npos)
      out.ctx.fpe_cache->bnd_kind = fpe::FPE_BND_LINEAR_ADDITIVE;
  } else if (type_std.find("RDMSWTN_TT") != std::string::npos) {
    // The time-changed model must precede generic RDMSWTN substring dispatch.
    out.pdf1_ptr       = &drdmswtn_tt_scalar;
    out.cdf1_ptr       = &prdmswtn_tt_scalar;
    out.model_dfun_raw = &drdmswtn_tt_raw;
    out.model_pfun_raw = &prdmswtn_tt_raw;
    out.logS_at_t_ptr  = &rdmswtn_tt_logS_at_t;
    out.col_spec       = emc2col::rdmswtn_tt::spec();
    out.ctx.t0_index   = emc2col::rdmswtn_tt::t0;
    out.ctx.defective_upper_tail = true;
    // "_CORRD" (correlated drift draws) contains "_CORR", so it must be
    // tested first; "_CORR" alone is the finishing-time copula.
    if (type_std.find("_CORRD") != std::string::npos) {
      out.ctx.corr_drift_active = true;
      out.ctx.corr_drift_v_col = emc2col::rdmswtn_tt::v;
      out.ctx.corr_drift_sv_col = emc2col::rdmswtn_tt::sv;
      out.ctx.corr_drift_generic_only = true;
    } else {
      out.ctx.rdmswtn_correlated =
        (type_std.find("_CORR") != std::string::npos);
    }
    if (type_std.find("_IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("RDMSWTN") != std::string::npos) {
    // Must be checked before "RDM" since "RDMSWTN" contains "RDM"
    out.pdf1_ptr       = &drdmswtn_scalar;
    out.cdf1_ptr       = &prdmswtn_scalar;
    out.model_dfun_raw = &drdmswtn_raw;
    out.model_pfun_raw = &prdmswtn_raw;
    out.logS_at_t_ptr  = &rdmswtn_logS_at_t;
    out.col_spec       = emc2col::rdmswtn::spec();
    out.ctx.t0_index   = emc2col::rdmswtn::t0;
    out.ctx.mean_g_index = emc2col::rdmswtn::mG;
    out.ctx.mean_k_index = emc2col::rdmswtn::mK;
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3) ? emc2col::rdmswtn::omega : -1;
    out.ctx.defective_upper_tail = true;
    // See the RDMSWTN_TT branch: "_CORRD" must be tested before "_CORR".
    if (type_std.find("_CORRD") != std::string::npos) {
      out.ctx.corr_drift_active = true;
      out.ctx.corr_drift_v_col = emc2col::rdmswtn::v;
      out.ctx.corr_drift_sv_col = emc2col::rdmswtn::sv;
      out.ctx.corr_drift_generic_only = true;
    } else {
      out.ctx.rdmswtn_correlated =
        (type_std.find("_CORR") != std::string::npos);
    }
    if (type_std.find("_IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
  } else if (type_std.find("GBM") != std::string::npos) {
    // Must be checked before "RDM" since "RDMGBM" contains "RDM"
    out.pdf1_ptr       = &drdmgbm_scalar;
    out.cdf1_ptr       = &prdmgbm_scalar;
    out.model_dfun_raw = &drdmgbm_raw;
    out.model_pfun_raw = &prdmgbm_raw;
    out.logS_at_t_ptr  = &rdmgbm_logS_at_t;
    out.col_spec       = emc2col::rdmgbm::spec();
    out.ctx.t0_index     = emc2col::rdmgbm::t0;
    out.ctx.mean_g_index = emc2col::rdmgbm::mG;
    out.ctx.mean_k_index = emc2col::rdmgbm::mK;
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3) ? emc2col::rdmgbm::omega : -1;
    out.ctx.defective_upper_tail = true;
  } else if (type_std.find("FRQ") != std::string::npos) {
    // Finite reservoir quorum.  "FRQ" is not a substring of any
    // other c_name and contains none, so its position among these branches is
    // free; it sits before BAwD only for readability.
    out.pdf1_ptr       = &dfrq_scalar;
    out.cdf1_ptr       = &pfrq_scalar;
    out.model_dfun_raw = &dfrq_raw;
    out.model_pfun_raw = &pfrq_raw;
    out.logS_at_t_ptr  = &frq_logS_at_t;
    out.col_spec       = emc2col::frq::spec();
    out.ctx.t0_index   = emc2col::frq::t0;
    // Always defective: the accumulator terminates only with probability
    // h = I_p(alpha, beta), and 1 - h of the mass sits at t = +Inf.  There is
    // no parameter setting that removes this, so the flag is unconditional.
    out.ctx.defective_upper_tail = true;
  } else if (type_std.find("BTAwL_SEPARATE") != std::string::npos) {
    out.pdf1_ptr       = &dbtawl_local_race_scalar;
    out.cdf1_ptr       = &pbtawl_local_race_scalar;
    out.model_dfun_raw = &dbtawl_local_race_raw;
    out.model_pfun_raw = &pbtawl_local_race_raw;
    out.logS_at_t_ptr  = &btawl_local_race_logS_at_t;
    const bool btawl_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool btawl_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool btawl_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = btawl_split
      ? emc2col::btawlsplit_local_race_separate::spec()
      : (btawl_weib ? emc2col::btawl_local_race_separate_weib::spec()
                    : (btawl_logn ? emc2col::btawl_local_race_separate_logn::spec()
                                  : emc2col::btawl_local_race_separate::spec()));
    out.ctx.t0_index = btawl_split ? int(emc2col::btawlsplit_local_race_separate::t0)
                                   : (btawl_weib ? int(emc2col::btawl_local_race_separate_weib::t0)
                                                 : (btawl_logn ? int(emc2col::btawl_local_race_separate_logn::t0)
                                                               : int(emc2col::btawl_local_race_separate::t0)));
    out.ctx.btawl_launch = btawl_split ? BTAWL_LAUNCH_SPLITLOGNORMAL
                         : (btawl_weib ? BTAWL_LAUNCH_WEIBULL
                                       : (btawl_logn ? BTAWL_LAUNCH_LOGNORMAL
                                                     : BTAWL_LAUNCH_NORMAL));
    out.ctx.btawl_separate = true;
    out.ctx.defective_upper_tail = true;
    if (!btawl_logn && type_std.find("_IO") != std::string::npos)
      out.ctx.use_posdrift = false;
    out.ctx.tw.supported = true;
  } else if (type_std.find("BTAwL_SUSTAINED") != std::string::npos) {
    out.pdf1_ptr       = &dbtawl_sustained_scalar;
    out.cdf1_ptr       = &pbtawl_sustained_scalar;
    out.model_dfun_raw = &dbtawl_sustained_raw;
    out.model_pfun_raw = &pbtawl_sustained_raw;
    out.logS_at_t_ptr  = &btawl_sustained_logS_at_t;
    const bool btawl_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool btawl_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool btawl_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = btawl_split
      ? emc2col::btawlsplit_sustained::spec()
      : (btawl_weib ? emc2col::btawl_sustained_weib::spec()
                    : (btawl_logn ? emc2col::btawl_sustained_logn::spec()
                                  : emc2col::btawl_sustained::spec()));
    out.ctx.t0_index = btawl_split ? int(emc2col::btawlsplit_sustained::t0)
                                   : int(emc2col::btawl_sustained::t0);
    out.ctx.btawl_launch = btawl_split ? BTAWL_LAUNCH_SPLITLOGNORMAL
                         : (btawl_weib ? BTAWL_LAUNCH_WEIBULL
                                       : (btawl_logn ? BTAWL_LAUNCH_LOGNORMAL
                                                     : BTAWL_LAUNCH_NORMAL));
    out.ctx.defective_upper_tail = true;
    if (!btawl_logn && type_std.find("_IO") != std::string::npos)
      out.ctx.use_posdrift = false;
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BTAwL_TRANSIENT") != std::string::npos) {
    out.pdf1_ptr       = &dbtawl_transient_scalar;
    out.cdf1_ptr       = &pbtawl_transient_scalar;
    out.model_dfun_raw = &dbtawl_transient_raw;
    out.model_pfun_raw = &pbtawl_transient_raw;
    out.logS_at_t_ptr  = &btawl_transient_logS_at_t;
    const bool btawl_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool btawl_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool btawl_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = btawl_split
      ? emc2col::btawlsplit_transient::spec()
      : (btawl_weib ? emc2col::btawl_transient_weib::spec()
                    : (btawl_logn ? emc2col::btawl_transient_logn::spec()
                                  : emc2col::btawl_transient::spec()));
    out.ctx.t0_index = btawl_split ? int(emc2col::btawlsplit_transient::t0)
                                   : int(emc2col::btawl_transient::t0);
    out.ctx.btawl_launch = btawl_split ? BTAWL_LAUNCH_SPLITLOGNORMAL
                         : (btawl_weib ? BTAWL_LAUNCH_WEIBULL
                                       : (btawl_logn ? BTAWL_LAUNCH_LOGNORMAL
                                                     : BTAWL_LAUNCH_NORMAL));
    out.ctx.defective_upper_tail = true;
    if (!btawl_logn && type_std.find("_IO") != std::string::npos)
      out.ctx.use_posdrift = false;
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BTAwL") != std::string::npos) {
    out.pdf1_ptr       = &dbtawl_local_race_scalar;
    out.cdf1_ptr       = &pbtawl_local_race_scalar;
    out.model_dfun_raw = &dbtawl_local_race_raw;
    out.model_pfun_raw = &pbtawl_local_race_raw;
    out.logS_at_t_ptr  = &btawl_local_race_logS_at_t;
    const bool btawl_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool btawl_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool btawl_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = btawl_split
      ? emc2col::btawlsplit_local_race::spec()
      : (btawl_weib ? emc2col::btawl_local_race_weib::spec()
                    : (btawl_logn ? emc2col::btawl_local_race_logn::spec()
                                  : emc2col::btawl_local_race::spec()));
    out.ctx.t0_index = btawl_split ? int(emc2col::btawlsplit_local_race::t0)
                                   : int(emc2col::btawl_local_race::t0);
    out.ctx.btawl_launch = btawl_split ? BTAWL_LAUNCH_SPLITLOGNORMAL
                         : (btawl_weib ? BTAWL_LAUNCH_WEIBULL
                                       : (btawl_logn ? BTAWL_LAUNCH_LOGNORMAL
                                                     : BTAWL_LAUNCH_NORMAL));
    out.ctx.defective_upper_tail = true;
    if (!btawl_logn && type_std.find("_IO") != std::string::npos)
      out.ctx.use_posdrift = false;
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BAwF") != std::string::npos) {
    // Global fading of the whole evidence trace: X(u) = h_rho(u)[z + V u].
    // "BAwF" contains and is contained by none of the other c_names, so its
    // position here is free; it sits with the other ballistic families.
    out.pdf1_ptr       = &dbawf_scalar;
    out.cdf1_ptr       = &pbawf_scalar;
    out.model_dfun_raw = &dbawf_raw;
    out.model_pfun_raw = &pbawf_raw;
    out.logS_at_t_ptr  = &bawf_logS_at_t;
    const bool bawf_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool bawf_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool bawf_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = bawf_split ? emc2col::bawfsplit::spec()
                  : (bawf_weib ? emc2col::bawf_weib::spec()
                               : (bawf_logn ? emc2col::bawf_logn::spec()
                                            : emc2col::bawf::spec()));
    out.ctx.t0_index = bawf_split ? int(emc2col::bawfsplit::t0)
                                  : int(emc2col::bawf::t0);
    // Shared with BAwD; see bawf_launch_of() in model_BAwF.cpp.
    out.ctx.bawd_launch = bawf_split ? BAWF_LAUNCH_SPLITLOGNORMAL
                        : (bawf_weib ? BAWF_LAUNCH_WEIBULL
                                     : (bawf_logn ? BAWF_LAUNCH_LOGNORMAL
                                                  : BAWF_LAUNCH_NORMAL));
    // Fixed fading-kernel shape from the c_name suffix; no suffix is the
    // exponential member.  rho = 1 has no finite endpoint and BAwF() refuses
    // it, so it emits no suffix here.
    out.ctx.bawd_rho =
      (type_std.find("_RHO2") != std::string::npos) ? 2.0 :
      ((type_std.find("_RHO4") != std::string::npos) ? 4.0 : R_PosInf);
    // Always defective: launches below V_c(z) never reach the threshold, and
    // no parameter setting removes that mass.
    out.ctx.defective_upper_tail = true;
    // posdrift is meaningless for the lognormal launch (V > 0 by
    // construction); BAwF() refuses posdrift = FALSE there rather than
    // silently ignoring it.
    if (!bawf_logn && !bawf_weib && type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BAwR") != std::string::npos) {
    // Ramping clearance: dX/du = V - kappa u^p.  "BAwR"
    // neither contains nor is contained by any other c_name, so its position
    // here is free; it sits with the other ballistic families.
    out.pdf1_ptr       = &dbawr_scalar;
    out.cdf1_ptr       = &pbawr_scalar;
    out.model_dfun_raw = &dbawr_raw;
    out.model_pfun_raw = &pbawr_raw;
    out.logS_at_t_ptr  = &bawr_logS_at_t;
    const bool bawr_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool bawr_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool bawr_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = bawr_split ? emc2col::bawrsplit::spec()
                  : (bawr_weib ? emc2col::bawr_weib::spec()
                               : (bawr_logn ? emc2col::bawr_logn::spec()
                                            : emc2col::bawr::spec()));
    out.ctx.t0_index = bawr_split ? int(emc2col::bawrsplit::t0)
                                  : int(emc2col::bawr::t0);
    // Shared with BAwD; see bawr_launch_of() in model_BAwR.cpp.  There is no rho:
    // the decay shape is the sampled exponent, not a fixed kernel index.
    out.ctx.bawd_launch = bawr_split ? BAWR_LAUNCH_SPLITLOGNORMAL
                        : (bawr_weib ? BAWR_LAUNCH_WEIBULL
                                     : (bawr_logn ? BAWR_LAUNCH_LOGNORMAL
                                                  : BAWR_LAUNCH_NORMAL));
    // Always defective: launches below V_c(z) never reach the threshold, and
    // no parameter setting removes that mass.
    out.ctx.defective_upper_tail = true;
    // posdrift is meaningless for the lognormal launch (V > 0 by
    // construction); BAwR() refuses posdrift = FALSE there rather than
    // silently ignoring it.
    if (!bawr_logn && !bawr_weib && type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BAwDp") != std::string::npos) {
    // BAwDp must precede the generic BAwD branch because its name contains the
    // string "BAwD".  It is an LBA evaluated at the closed-form internal clock
    // m(t), with m'(t) as the Jacobian and a finite frozen ceiling.
    out.pdf1_ptr       = &dbawdp_scalar;
    out.cdf1_ptr       = &pbawdp_scalar;
    out.model_dfun_raw = &dbawdp_raw;
    out.model_pfun_raw = &pbawdp_raw;
    out.logS_at_t_ptr  = &bawdp_logS_at_t;
    const bool bawdp_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool bawdp_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool bawdp_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec = bawdp_split ? emc2col::bawdpsplit::spec()
                  : (bawdp_weib ? emc2col::bawdp_weib::spec()
                                : (bawdp_logn ? emc2col::bawdp_logn::spec()
                                              : emc2col::bawdp::spec()));
    out.ctx.t0_index = bawdp_split ? int(emc2col::bawdpsplit::t0)
                                   : int(emc2col::bawdp::t0);
    out.ctx.bawl_launch = bawdp_split ? BAWL_LAUNCH_SPLITLOGNORMAL
                        : (bawdp_weib ? BAWL_LAUNCH_WEIBULL
                                      : (bawdp_logn ? BAWL_LAUNCH_LOGNORMAL
                                                    : BAWL_LAUNCH_NORMAL));
    out.ctx.defective_upper_tail = true;
    if (!bawdp_logn && !bawdp_weib && type_std.find("IO") != std::string::npos)
      out.ctx.use_posdrift = false;
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BAwD") != std::string::npos) {
    // Dispatch is by substring, and "BAwD" is a substring of nothing here and
    // contains neither "BAwL" nor "LBA", so placement relative to those is
    // safe either way; it sits next to BAwL for readability.
    out.pdf1_ptr       = &dbawd_scalar;
    out.cdf1_ptr       = &pbawd_scalar;
    out.model_dfun_raw = &dbawd_raw;
    out.model_pfun_raw = &pbawd_raw;
    out.logS_at_t_ptr  = &bawd_logS_at_t;
    // Normal and lognormal launches share column POSITIONS; split-lognormal
    // inserts delta immediately after sigma, shifting the rest by one.
    const bool bawd_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool bawd_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool bawd_weib = (type_std.find("_WEIB") != std::string::npos);
    out.ctx.t0_index = bawd_split ? int(emc2col::bawdsplit::t0) :
      (bawd_weib ? int(emc2col::bawd_weib::t0) : int(emc2col::bawd::t0));
    out.ctx.bawd_launch = bawd_split ? BAWD_LAUNCH_SPLITLOGNORMAL
                        : (bawd_weib ? BAWD_LAUNCH_WEIBULL
                                     : (bawd_logn ? BAWD_LAUNCH_LOGNORMAL
                                                  : BAWD_LAUNCH_NORMAL));
    // Fixed clearance exponent and base-kernel shape parsed from c_name.
    // Constant clearance and the exponential kernel emit no suffixes,
    // preserving existing BAwD routing.
    out.ctx.bawd_gamma =
      (type_std.find("_GAM100") != std::string::npos) ? 1.0 :
      ((type_std.find("_GAM34") != std::string::npos) ? 0.75 :
       ((type_std.find("_GAM23") != std::string::npos) ? (2.0 / 3.0) :
        ((type_std.find("_GAM12") != std::string::npos) ? 0.5 : 0.0)));
    // Slot 6 samples the clearance rate `ell` for every gamma/rho option; the
    // endpoint is derived on the R side and never stored in a design column.
    out.col_spec = bawd_split ? emc2col::bawdsplit::spec()
                  : (bawd_weib ? emc2col::bawd_weib::spec()
                               : (bawd_logn ? emc2col::bawd_logn::spec()
                                            : emc2col::bawd::spec()));
    // Fixed power-decay kernel parameter parsed from the c_name suffix.
    // Default (no suffix) is R_PosInf (exponential kernel).
    out.ctx.bawd_rho =
      (type_std.find("_RHO1") != std::string::npos) ? 1.0 :
      ((type_std.find("_RHO2") != std::string::npos) ? 2.0 :
       ((type_std.find("_RHO4") != std::string::npos) ? 4.0 : R_PosInf));
    // Always defective: weak launch strengths can miss the threshold in every
    // regime; co-decay removes the finite wall but not the omission mass.
    out.ctx.defective_upper_tail = true;
    // posdrift is meaningless for the lognormal launch (V > 0 by construction);
    // BAwD() refuses posdrift = FALSE there rather than silently ignoring it.
    if (!bawd_logn && !bawd_weib && type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
    out.ctx.tw.supported = true;   // operational-time warp (Math/ballistic-time.md)
  } else if (type_std.find("BAwL") != std::string::npos) {
    out.pdf1_ptr       = &dbawl_scalar;
    out.cdf1_ptr       = &pbawl_scalar;
    out.model_dfun_raw = &dbawl_raw;
    out.model_pfun_raw = &pbawl_raw;
    out.logS_at_t_ptr  = &bawl_logS_at_t;
    // The lognormal launch shares every column POSITION with the Gaussian one
    // (mu/sigma occupy v/sv); split-lognormal inserts delta after sigma, so
    // the spec, the launch flag, and the shifted indices all differ.
    const bool bawl_split = (type_std.find("_SPLIT") != std::string::npos);
    const bool bawl_logn = (type_std.find("_LOGN") != std::string::npos);
    const bool bawl_weib = (type_std.find("_WEIB") != std::string::npos);
    out.col_spec       = bawl_split ? emc2col::bawlsplit::spec()
                        : (bawl_weib ? emc2col::bawl_weib::spec()
                                     : (bawl_logn ? emc2col::bawl_logn::spec()
                                                  : emc2col::bawl::spec()));
    out.ctx.t0_index   = bawl_split ? int(emc2col::bawlsplit::t0) :
      (bawl_weib ? int(emc2col::bawl_weib::t0) : int(emc2col::bawl::t0));
    out.ctx.mean_g_index = bawl_split ? int(emc2col::bawlsplit::mG) :
      (bawl_weib ? int(emc2col::bawl_weib::mG) : int(emc2col::bawl::mG));
    out.ctx.mean_k_index = bawl_split ? int(emc2col::bawlsplit::mK) :
      (bawl_weib ? int(emc2col::bawl_weib::mK) : int(emc2col::bawl::mK));
    out.ctx.erlang_omega_index = (out.ctx.kill_shape == 3)
      ? (bawl_split ? int(emc2col::bawlsplit::omega) :
         (bawl_weib ? int(emc2col::bawl_weib::omega) : int(emc2col::bawl::omega))) : -1;
    out.ctx.bawl_launch = bawl_split ? BAWL_LAUNCH_SPLITLOGNORMAL
                        : (bawl_weib ? BAWL_LAUNCH_WEIBULL
                                     : (bawl_logn ? BAWL_LAUNCH_LOGNORMAL
                                                  : BAWL_LAUNCH_NORMAL));
    // The correlated-drift path is a one-factor decomposition of the *Gaussian*
    // drift vector, so it is not reachable with a lognormal launch; BAwL()
    // rejects that combination rather than silently ignoring one of them.
    out.ctx.corr_drift_active = (type_std.find("_CORR") != std::string::npos);
    out.ctx.corr_drift_v_col = emc2col::bawl::v;
    out.ctx.corr_drift_sv_col = emc2col::bawl::sv;
    // Leaky ballistic accumulators can have defective upper tails (never-finish
    // mass) even when posdrift=TRUE.
    out.ctx.defective_upper_tail = true;
    if (!bawl_logn && !bawl_weib && type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
    }
    if (type_std.find("_E2") != std::string::npos) out.ctx.kill_shape = 2;
    if (type_std.find("_EMIX") != std::string::npos) out.ctx.kill_shape = 3;
    // The Erlang kill/guess clocks run on RAW time and the correlated path has
    // its own inlined kernels, so neither composes with the operational-time
    // warp yet.  The R constructor already withholds `eta` from those variants;
    // this is the defensive half of the same contract.
    out.ctx.tw.supported = !out.ctx.kill_active && !out.ctx.corr_drift_active;
  } else if (type_std.find("LBA") != std::string::npos) {
    // Standard LBA is the exact k=0, no-clock member of the shared BAwL
    // family.  Keep the five-column LBA contract and force the optional BAwL
    // parameters off in the adapter rather than adding hidden columns.
    out.pdf1_ptr = &dbawl_scalar;
    out.cdf1_ptr = &pbawl_scalar;
    out.model_dfun_raw = &dbawl_raw;
    out.model_pfun_raw = &pbawl_raw;
    out.logS_at_t_ptr = &bawl_logS_at_t;
    out.col_spec     = emc2col::lba::spec();
    out.ctx.t0_index = emc2col::lba::t0;
    out.ctx.bawl_k_fixed_zero = true;
    out.ctx.bawl_clocks_fixed_off = true;
    out.ctx.defective_upper_tail = false;
    if (type_std.find("IO") != std::string::npos) {
      out.ctx.use_posdrift = false;
      out.ctx.defective_upper_tail = true;
    }
    // LogicalRulesLBA shares the "LBA" substring but has a separate compiled
    // time bookkeeping path; keep eta explicitly unsupported there.
    out.ctx.tw.supported = (type_std.find("LogicalRules") == std::string::npos);
  } else if (type_std.find("RDM") != std::string::npos) {
    out.pdf1_ptr = &drdm_scalar;
    out.cdf1_ptr = &prdm_scalar;
    out.model_dfun_raw = &drdm_raw;
    out.model_pfun_raw = &prdm_raw;
    out.logS_at_t_ptr = &rdm_logS_at_t;
    out.col_spec     = emc2col::rdm::spec();
    out.ctx.t0_index = emc2col::rdm::t0;
  } else if (type_std.find("REXG") != std::string::npos) {
    out.pdf1_ptr = &dexg_scalar;
    out.cdf1_ptr = &pexg_scalar;
    out.model_dfun_raw = &drexg_raw;
    out.model_pfun_raw = &prexg_raw;
    out.logS_at_t_ptr = &rexg_logS_at_t;
    out.col_spec = emc2col::rexg::spec();
    out.ctx.t0_index = emc2col::rexg::t0;
  } else if (type_std.find("LNR") != std::string::npos) {
    out.pdf1_ptr = &dlnr_scalar;
    out.cdf1_ptr = &plnr_scalar;
    out.model_dfun_raw = &dlnr_raw;
    out.model_pfun_raw = &plnr_raw;
    out.logS_at_t_ptr = &lnr_logS_at_t;
    out.col_spec     = emc2col::lnr::spec();
    out.ctx.t0_index = emc2col::lnr::t0;
  } else if (type_std.find("PCOUNTER") != std::string::npos) {
    out.pdf1_ptr = &dpcounter_scalar;
    out.cdf1_ptr = &ppcounter_scalar;
    out.model_dfun_raw = &dpcounter_raw;
    out.model_pfun_raw = &ppcounter_raw;
    out.logS_at_t_ptr = &pcounter_logS_at_t;
    out.col_spec     = emc2col::pcounter::spec();
    // PCOUNTER's t0 is an ordinary additive shift, so truncation/censoring
    // integration can skip its zero-density dead zone just like FRQ/LNR.
    out.ctx.t0_index = emc2col::pcounter::t0;
  } else {
    Rcpp::stop("Unsupported race model type string in %s: %s", caller.c_str(), type_std.c_str());
  }

  if (type_std.find("GNG") != std::string::npos) {
    out.ctx.gng = true;
  }
  return out;
}

void configure_corr_drift_context(RaceModelAdapter& adapter,
                                               const Rcpp::CharacterVector& keep_names,
                                               const std::string& caller) {
  if (!adapter.ctx.corr_drift_active) return;
  adapter.ctx.corr_drift_rho_index = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "rho") {
      adapter.ctx.corr_drift_rho_index = j;
      break;
    }
  }
  if (adapter.ctx.corr_drift_rho_index < 0) {
    Rcpp::stop("%s: a correlated-drift model requires a parameter column "
               "named 'rho'.", caller.c_str());
  }
}

void configure_rdmswtn_corr_context(
    RaceModelAdapter& adapter, const Rcpp::CharacterVector& keep_names,
    const std::string& caller) {
  if (!adapter.ctx.rdmswtn_correlated) return;
  adapter.ctx.rdmswtn_rho_index = -1;
  for (int j = 0; j < keep_names.size(); ++j) {
    if (Rcpp::as<std::string>(keep_names[j]) == "rho") {
      adapter.ctx.rdmswtn_rho_index = j;
      break;
    }
  }
  if (adapter.ctx.rdmswtn_rho_index < 0) {
    Rcpp::stop("%s: correlated RDMSWTN requires a parameter column named 'rho'.",
               caller.c_str());
  }
  if (adapter.ctx.kill_active) {
    Rcpp::stop("%s: correlated RDMSWTN does not support guess or kill clocks.",
               caller.c_str());
  }
}

bool is_stop_signal_type(const std::string& type_std) {
  return type_std == "SSEXG" || type_std == "SSRDEX";
}
