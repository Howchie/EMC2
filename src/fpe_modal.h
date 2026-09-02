#ifndef fpe_modal_h
#define fpe_modal_h

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modal propagator for a FIXED-boundary first-passage solve.
//
// When the boundary does not move, the finite-volume generator L built by
// fpe::build_op() is autonomous, and the whole solve is
//
//   q(t) = exp((t - t0) L) q0
//
// read through two linear functionals: the absorbing-face flux (the density)
// and the cell-mass sum (the survivor).  The Crank-Nicolson march in
// fpe_solver.h discretises the time axis to get there, at O(M * n_t); this
// module does not discretise time at all.  It projects the problem onto a small
// Krylov space, diagonalises the projection, and every requested time is then a
// sum of m real exponentials.
//
// Why this is possible here and not for a moving boundary: L's adjacent
// off-diagonals are both strictly positive -- L_{i,i-1} = D r_i bern(-P_i) dx_i^-1
// and L_{i-1,i} = D r_i bern(P_i) dx_{i-1}^-1, and bern > 0 everywhere -- so L is
// diagonally similar to a SYMMETRIC tridiagonal
//
//   S = W^-1 L W,   W_i / W_{i-1} = sqrt( L_{i,i-1} / L_{i-1,i} ),
//   S_{i,i-1} = S_{i-1,i} = sqrt( L_{i,i-1} L_{i-1,i} ).
//
// That buys three things at once.  Lanczos replaces Arnoldi, so the basis costs
// a three-term recurrence rather than a full Gram-Schmidt; the projected matrix
// is symmetric tridiagonal, so its eigendecomposition is a QL sweep rather than
// a Hessenberg QR plus back-transformation; and every eigenvalue is REAL, so
// the modal sums carry no complex arithmetic.  Compare model_RLF.h, which runs
// the same idea on a nonsymmetric dense generator and pays DHSEQR + DTREVC for
// it.
//
// The ratio has a closed form.  bern(-P)/bern(P) = e^P exactly, so
//
//   W_i / W_{i-1} = sqrt( dx_{i-1} / dx_i ) * exp(P_i / 2),
//
// i.e. log W is the discrete potential (1/sigma^2) * integral of the drift,
// which for OU is an inverted parabola centred on the equilibrium v/lambda.
// Its total SPAN over the domain is what decides whether the similarity is
// numerically usable, and it is reported so the caller can refuse the transform
// rather than return noise -- see fpe_modal_symmetrise().
//
// Shift-and-invert.  The Krylov space is built on (I - gamma S)^-1 rather than
// on S, which is one Thomas solve per vector and puts the resolution where the
// solution lives: the slow modes that survive to the response times, rather
// than the fast ones that have decayed by the time the first RT arrives.  Same
// device, same reasoning, as rlf_build_modes().
//
// Accuracy.  The modal answer carries NO time-discretisation error, so on a
// fixed boundary it is not an approximation to the march -- the march is an
// approximation to it.  Measured against a 64x-finer-in-time CN reference on
// the shipped grid, at m = 24 the modal route is 1e-5 to 1e-8 in |d log f| for
// t >= 0.15 s where the production march is 8e-5 to 9e-4.  See the ladder note
// on fpe_modal_build() for how m is chosen and what makes a solve refuse.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "fpe_solver.h"

namespace fpe {

// --- ladder and acceptance knobs ---------------------------------------------
// The ladder starts low and grows, because the cost of a rung is the QL sweep
// (O(m^2) rotations, the dominant per-key cost once the Lanczos itself is lane
// batched) and the Lanczos vectors are shared: the leading k x k block of the
// tridiagonal built at m IS the tridiagonal k steps would have produced, so a
// lower rung is a free re-read of work already done.
inline int FPE_MODAL_MIN = 16;
inline int FPE_MODAL_STEP = 6;
inline int FPE_MODAL_MAX = 48;
// The ladder is ADDITIVE, and deliberately so.  A rung costs O(m^2), so a
// geometric ladder looks like the obvious saving -- climbing 16, 22, 28, 34, 40
// to settle at 40 spends 96 us of QL where 16, 24, 36 spends 46.  It was tried
// and it loses, because the intermediate rungs are not waste: they are the
// chances to ACCEPT at a lower dimension, and the movement test is between
// CONSECUTIVE rungs, so widening the gaps also pushes the last comparison up
// against the cap.  Measured on an ordinary parameter net, growth factors of
// 1.5 and 1.75 turned 0.5% and 12% of keys from accepted into refused (they ran
// out of ladder before two rungs agreed) and made the mean solve slower, not
// faster, once the refused keys were marched.
inline int fpe_modal_next_rung(int m) { return m + FPE_MODAL_STEP; }
// Movement between successive rungs, in the currencies the likelihood reads it
// in: relative on the density against a denominator floored at FPE_MODAL_FLOOR
// of the peak, so the leading edge stays in the test without the deep tail
// dominating it, and absolute on log S, which is what the omission and
// truncation paths read.
//
// The tolerance is tight, and loosening it is not the saving it looks like.
// Measured over a random ROU net: at 1e-5 the worst accepted density error is
// 2.4e-2 and the mean solve is 127 us; at 1e-4 the worst is 0.80 and at 1e-3 it
// is 15.3 nats, for 121 and 106 us.  The mean dimension barely moves -- the
// tolerance is not buying vectors on ordinary keys, it is catching the handful
// of keys that need them.
inline double FPE_MODAL_TOL = 1e-5;
inline double FPE_MODAL_FLOOR = 1e-4;
// A floored test is an ABSOLUTE test wherever the density is below the floor,
// so on its own it will accept a rung whose leading-edge density is out by a
// large FACTOR as long as it is small in absolute terms.  Usually that is the
// right trade -- a density twelve orders below the peak is not what a fit turns
// on -- but a trial CAN land there, and then the error is nats of likelihood,
// not parts per million.  Measured on an ordinary parameter net with response
// times reaching 50 ms past the start, the floored test alone leaves a worst
// case of 2.0 nats in the earliest fifteen per cent of the horizon, against
// 0.24 for the march it is replacing.  So the movement in LOG density is
// required to be small as well, at a tolerance loose enough not to chase the
// noise floor.  Measured over that net, the worst early-time error falls to
// 1.7e-1 at 1e-1, 2.3e-2 at 1e-2 and 1.8e-3 at 1e-3, for a cost of about five
// per cent in mean solve time -- it does not move the accepted dimension, it
// sends the handful of keys that were being waved through to the march.
inline double FPE_MODAL_TOL_LOG = 1e-3;
// Shift, as a fraction of the horizon.  Small shifts resolve the slow end of
// the spectrum, which is the end that survives to a response time.
inline double FPE_MODAL_SHIFT = 0.02;
// Largest sum|term| / |sum| at which a modal value is still believed.  The
// expansion is a signed sum, so a value obtained by cancelling 15 digits is
// noise however well the Krylov space has converged; this is what separates
// "the density really is that small" from "the modal form cannot say".
inline double FPE_MODAL_CANCEL = 1e12;
// Modal coefficients can be large because the diagonal similarity is allowed
// to span hundreds of nats.  Evaluating coefficient * exp(lambda * t) as two
// operations therefore loses terms when the exponential underflows before the
// coefficient rescales it.  Form such tail products in log space and discard
// them only once each is ten nats below the likelihood's -700 log floor.
inline double FPE_MODAL_TERM_LOG_FLOOR = -710.0;
// Largest span of log W (nats) at which the diagonal similarity is attempted.
// The transform maps the seed by e^-logW and the functionals by e^+logW, so the
// span is exactly the dynamic range it costs, and past ~700 nats one of the two
// overflows outright.  It is set well inside that as a STRUCTURAL guard and
// nothing more: what actually decides whether a wide-span solve is usable is
// FPE_MODAL_CANCEL, measured on the answer rather than predicted from the
// transform.  Tightening the cap instead is worse on both counts -- over a
// random ROU net that reaches s = 0.3, a cap of 45 refused 33% of keys where
// the cancellation test refuses 20%, and the 13% it refused needlessly were
// solves the modal route was answering to 1e-5 where the march was out by 0.2
// nats.  Ordinary ROU parameters land at 7-31 nats; 67 (a convection-dominated
// cell at s = 0.2) is where the density starts to go and 560 (s = 0.05) is
// where both functionals do -- and the march is in trouble in both of those.
inline double FPE_MODAL_SPAN = 400.0;
// Most probes the convergence ladder is judged at.  The probe set is drawn from
// the times the caller actually wants, so it lands where the data is; the cap
// exists because a rung costs np * m exponentials and a key can carry thousands
// of distinct response times.  The ACCEPTED curve is then checked at every
// requested time before it is used, so the cap trades ladder cost for a final
// sweep, not for coverage.
inline int FPE_MODAL_PROBES = 24;

// Draw at most FPE_MODAL_PROBES probe offsets from `times` (absolute, sorted or
// not), keeping the first and last of those that lie above t0.
inline void fpe_modal_probes(const std::vector<double>& times, double t0,
                             std::vector<double>& probe) {
  probe.clear();
  std::vector<double> use;
  use.reserve(times.size());
  for (double t : times) if (t > t0) use.push_back(t - t0);
  if (use.empty()) return;
  std::sort(use.begin(), use.end());
  const int n = static_cast<int>(use.size());
  const int want = std::min(n, std::max(2, FPE_MODAL_PROBES));
  if (n <= want) { probe = use; return; }
  probe.reserve(want);
  for (int k = 0; k < want; ++k) {
    const int j = static_cast<int>(
      std::llround(static_cast<double>(k) * (n - 1) / (want - 1)));
    if (probe.empty() || use[j] > probe.back()) probe.push_back(use[j]);
  }
}

// Why a solve did not take the modal route, for the diagnostic counters.
enum : int {
  FPE_MODAL_OK = 0,
  FPE_MODAL_NOT_SYMMETRISABLE = 1,
  FPE_MODAL_SPAN_TOO_WIDE = 2,
  FPE_MODAL_SEED_UNUSABLE = 3,
  FPE_MODAL_QL_FAILED = 4,
  FPE_MODAL_NOT_CONVERGED = 5,
  FPE_MODAL_CANCELLING = 6
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Implicit-QL with Wilkinson shift on a symmetric tridiagonal, applying every
// Givens rotation to `nv` row vectors held in R (nv x n, row major) instead of
// accumulating the eigenvector matrix.
//
// Only three projections of the eigenbasis are ever wanted -- the two output
// functionals and e_1 -- so accumulating Z would be O(n^3) work to extract 3n
// numbers.  Rotating three vectors alongside the sweep is O(3 n^2) and is what
// keeps the reduction from dominating a lane-batched solve: measured at n = 32,
// 22 us here against 64 us for LAPACK's DSTEV, which forms all n eigenvectors.
//
//   d[0..n-1]  diagonal; overwritten with the eigenvalues, unsorted
//   e[0..n-2]  off-diagonal, e[i] between d[i] and d[i+1]; destroyed
//   R          on entry the nv vectors; on exit R[k*n + j] = v_k . z_j
//
// The caller's tridiagonal comes from Lanczos on (I - gamma S)^-1, whose
// spectrum lies in (0, 1], so no rotation can overflow the unscaled
// sqrt(f*f + g*g) and the branchy scaled pythag is not needed on what is
// already a latency-bound chain.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inline bool fpe_tql_project(int n, double* d, double* e, double* R, int nv) {
  if (n <= 0) return false;
  if (n == 1) return true;
  e[n - 1] = 0.0;
  const double tiny = std::numeric_limits<double>::epsilon();
  auto pythag = [](double a, double b) { return std::sqrt(a * a + b * b); };
  for (int l = 0; l < n; ++l) {
    int iter = 0, m = l;
    double r = 0.0;
    do {
      for (m = l; m + 1 < n; ++m) {
        const double dd = std::abs(d[m]) + std::abs(d[m + 1]);
        if (std::abs(e[m]) <= tiny * dd) break;
      }
      if (m == l) break;
      if (++iter > 50) return false;
      double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
      r = pythag(g, 1.0);
      g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));
      double s = 1.0, c = 1.0, p = 0.0;
      int i = m - 1;
      for (; i >= l; --i) {
        const double f = s * e[i];
        const double b = c * e[i];
        r = pythag(f, g);
        e[i + 1] = r;
        if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
        const double rinv = 1.0 / r;
        s = f * rinv;
        c = g * rinv;
        g = d[i + 1] - p;
        const double q = (d[i] - g) * s + 2.0 * c * b;
        p = s * q;
        d[i + 1] = g + p;
        g = c * q - b;
        for (int k = 0; k < nv; ++k) {
          double* row = R + static_cast<size_t>(k) * n;
          const double a1 = row[i + 1], a0 = row[i];
          row[i + 1] = s * a0 + c * a1;
          row[i]     = c * a0 - s * a1;
        }
      }
      if (r == 0.0 && i >= l) continue;
      d[l] -= p;
      e[l] = g;
      e[m] = 0.0;
    } while (true);
  }
  return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The finished propagator: two sums of m real exponentials, evaluable at any t.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct FPE_ModalCurve {
  int m = 0;
  int status = FPE_MODAL_NOT_CONVERGED;
  double t0 = 0.0;            // the seed time the exponentials are measured from
  double surv0 = 1.0;         // surviving mass at t0, the clamp for S(t)
  double span = 0.0;          // dynamic range of log W, nats (diagnostic)
  std::vector<double> lam;    // reduced-generator eigenvalues, all real and <= 0
  std::vector<double> cf, cs; // modal coefficients: absorbing flux, survivor

  bool ok() const { return status == FPE_MODAL_OK && m > 0; }

  static double scaled_term(double coefficient, double exponent) {
    if (coefficient == 0.0) return 0.0;
    // The ordinary product is both faster and more accurate away from the
    // underflow range.  In the tail, combine the exponents before calling exp
    // so a large modal coefficient can rescue a term that is still material.
    if (exponent >= -700.0) return coefficient * std::exp(exponent);
    const double log_term = std::log(std::abs(coefficient)) + exponent;
    if (log_term < FPE_MODAL_TERM_LOG_FLOOR) return 0.0;
    return std::copysign(std::exp(log_term), coefficient);
  }

  // Both functionals share one exponential sweep; a race likelihood always
  // wants the pair, and exp() is the whole cost of an evaluation.
  void at(double dt, double& pdf, double& surv) const {
    double f = 0.0, s = 0.0;
    for (int j = 0; j < m; ++j) {
      const double x = lam[j] * dt;
      f += scaled_term(cf[j], x);
      s += scaled_term(cs[j], x);
    }
    pdf = f;
    surv = s;
  }

  // Evaluate at MANY times at once.  `dt` must be sorted ascending.
  //
  // Mode-major rather than time-major, which is what makes this the cheap way
  // to answer a key that carries hundreds of distinct response times.  Two
  // things follow from the ordering:
  //
  //  * the inner loop runs over times at a fixed rate, so it is a flat
  //    exp() sweep -- and since every argument is non-positive it can go
  //    through exp_nonpos_v(), the same branch-free polynomial the lane-batched
  //    operator build uses, where the time-major loop would carry three
  //    accumulators and vectorise over nothing;
  //  * within a mode, the times past which BOTH functional contributions are
  //    safely below the output floor are found by one division.  This cutoff
  //    must include the coefficient: lambda alone is insufficient when the
  //    diagonal similarity has produced a large modal coefficient.
  void at_many(const std::vector<double>& dt, std::vector<double>& pdf,
               std::vector<double>& surv, std::vector<double>& amp) const {
    const size_t n = dt.size();
    pdf.assign(n, 0.0);
    surv.assign(n, 0.0);
    amp.assign(n, 0.0);
    if (n == 0) return;
    for (int j = 0; j < m; ++j) {
      const double lj = lam[j];
      const double a = cf[j], b = cs[j];
      const double coef = std::max(std::abs(a), std::abs(b));
      if (!(coef > 0.0)) continue;
      size_t stop = n;
      if (lj < 0.0) {
        const double t_lim =
          (FPE_MODAL_TERM_LOG_FLOOR - std::log(coef)) / lj;
        stop = static_cast<size_t>(
          std::upper_bound(dt.begin(), dt.end(), t_lim) - dt.begin());
      }
      // exp_nonpos_v deliberately flushes below -708 because its exponent-bit
      // scaling cannot create subnormals.  Large coefficients can make those
      // terms material, so only the safe prefix takes the vector polynomial;
      // scaled_term() handles the remaining tail without premature underflow.
      size_t vector_stop = stop;
      if (lj < 0.0) {
        const double t_lim = -708.0 / lj;
        vector_stop = std::min(
          stop, static_cast<size_t>(
            std::upper_bound(dt.begin(), dt.end(), t_lim) - dt.begin()));
      }
      size_t i = 0;
#if defined(__AVX512F__) && defined(__FMA__)
      {
        const __m512d lv = _mm512_set1_pd(lj), av = _mm512_set1_pd(a),
                      bv = _mm512_set1_pd(b);
        const __m512d sign = _mm512_castsi512_pd(
            _mm512_set1_epi64(0x7fffffffffffffffLL));
        for (; i + 8 <= vector_stop; i += 8) {
          const __m512d w =
            exp_nonpos_v(_mm512_mul_pd(lv, _mm512_loadu_pd(&dt[i])));
          const __m512d x = _mm512_mul_pd(av, w);
          _mm512_storeu_pd(&pdf[i], _mm512_add_pd(_mm512_loadu_pd(&pdf[i]), x));
          _mm512_storeu_pd(&amp[i], _mm512_add_pd(_mm512_loadu_pd(&amp[i]),
                                                 _mm512_and_pd(x, sign)));
          _mm512_storeu_pd(&surv[i],
                           _mm512_fmadd_pd(bv, w, _mm512_loadu_pd(&surv[i])));
        }
      }
#elif defined(__AVX2__) && defined(__FMA__)
      {
        const __m256d lv = _mm256_set1_pd(lj), av = _mm256_set1_pd(a),
                      bv = _mm256_set1_pd(b);
        const __m256d sign = _mm256_set1_pd(-0.0);
        for (; i + 4 <= vector_stop; i += 4) {
          const __m256d w =
            exp_nonpos_v(_mm256_mul_pd(lv, _mm256_loadu_pd(&dt[i])));
          const __m256d x = _mm256_mul_pd(av, w);
          _mm256_storeu_pd(&pdf[i], _mm256_add_pd(_mm256_loadu_pd(&pdf[i]), x));
          _mm256_storeu_pd(&amp[i], _mm256_add_pd(_mm256_loadu_pd(&amp[i]),
                                                 _mm256_andnot_pd(sign, x)));
          _mm256_storeu_pd(&surv[i],
                           _mm256_fmadd_pd(bv, w, _mm256_loadu_pd(&surv[i])));
        }
      }
#endif
      for (; i < stop; ++i) {
        const double exponent = lj * dt[i];
        const double x = scaled_term(a, exponent);
        pdf[i] += x;
        amp[i] += std::abs(x);
        surv[i] += scaled_term(b, exponent);
      }
    }
    for (size_t i = 0; i < n; ++i)
      amp[i] = (pdf[i] != 0.0) ? amp[i] / std::abs(pdf[i]) : HUGE_VAL;
  }

  // The same sweep, plus sum|term| / |sum| for the density: the cancellation
  // the log is being asked to survive.  Only the acceptance test wants the
  // third number, and it wants it at the same times as the first two, so it
  // rides along on the one set of exponentials rather than repeating them.
  void probe(double dt, double& pdf, double& surv, double& amp) const {
    double f = 0.0, s = 0.0, a = 0.0;
    for (int j = 0; j < m; ++j) {
      const double exponent = lam[j] * dt;
      const double x = scaled_term(cf[j], exponent);
      f += x;
      a += std::abs(x);
      s += scaled_term(cs[j], exponent);
    }
    pdf = f;
    surv = s;
    amp = (f != 0.0) ? a / std::abs(f) : HUGE_VAL;
  }
};

// Reusable scratch.  A cache holds one of these and every solve reuses it, so
// the M-sized buffers are allocated once per likelihood rather than per key.
struct FPE_ModalWork {
  std::vector<double> off, logw, phif, phis;   // M
  std::vector<double> v, vprev, w;             // M
  std::vector<double> sub, dinv, cprime;       // M
  std::vector<double> alpha, beta, gf, gs;     // m_max
  std::vector<double> qd, qe, qr;              // QL scratch
  std::vector<int> order;                      // mode permutation, slowest first
  // 0.5 * log(dx_{i-1} / dx_i) at face i.  Mesh-only, and every key in a batch
  // shares the mesh, so it is built once and reused -- see
  // fpe_modal_symmetrise() for why it is nearly the only logarithm the
  // transform needs.  Zero throughout on a uniform mesh.
  std::vector<double> mesh_half_log_ratio;
  const FPE_Mesh* mesh_cached = nullptr;
  int mesh_M = 0;

  void ensure_mesh(const FPE_Mesh& g) {
    if (mesh_cached == &g && mesh_M == g.M) return;
    mesh_half_log_ratio.assign(g.M, 0.0);
    if (!g.uniform) {
      for (int i = 1; i < g.M; ++i)
        mesh_half_log_ratio[i] = 0.5 * std::log(g.dx[i - 1] / g.dx[i]);
    }
    mesh_cached = &g;
    mesh_M = g.M;
  }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Symmetrise, and map the seed and the two functionals into the symmetric
// coordinates.
//
//   q(t) = W q~(t),   q~(t) = exp(t S) q~(0),   phi^T q = (W phi)^T q~
//
// so the seed is divided by W and the functionals multiplied by it.  log W is
// normalised at the PEAK OF THE SEED, which is the one choice that keeps
// q~(0) = O(q0): the seed is a narrow bump and everything else in the transform
// is measured relative to where it sits.
//
// Returns an FPE_MODAL_* status.  A span wider than FPE_MODAL_SPAN is refused
// here rather than allowed to produce a plausible-looking answer -- past that
// point the far field of q~ is amplified by e^span while the functional is
// damped by the same factor, and the products the modal sum has to form are
// cancellation noise.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// `work.v` is left holding the NORMALISED transformed seed, which is also the
// first Lanczos vector; fpe_modal_extend() picks it up from there.
//
// log W is accumulated from the face Peclet numbers rather than from
// log(L_{i,i-1} / L_{i-1,i}), which is the same number without the logarithm.
// At an interior face the two off-diagonals are
//   L_{i,i-1} = D r_i bern(-P_i) / dx_i,  L_{i-1,i} = D r_i bern(P_i) / dx_{i-1}
// and bern(-P)/bern(P) = e^P identically, so
//   log W_i - log W_{i-1} = P_i / 2 + log(dx_{i-1} / dx_i) / 2.
// The first term is affine in xi -- build_op() has to form it anyway -- and the
// second depends only on the mesh, which every key in a batch shares.  That
// takes M logarithms per key out of the setup, which with the paired
// exponentials below was measured as the largest single cost in a lane-batched
// solve after the QL reductions themselves.
//
// The two rows that TOUCH A BARRIER are not interior faces and the identity
// does not hold there: build_op() adds the one-sided quadratic closure weight
// (g.cB at the absorbing row, g.cL_B at the lower one when the model has it) to
// the same off-diagonal, so the ratio is no longer e^P.  Measured, using the
// identity at the absorbing row anyway misplaces log W by 0.18 nats -- a factor
// of 1.2 on the functional that reads the barrier.  Those two faces therefore
// take the exact logarithm, which costs two of them per solve rather than M.
template <class Model>
inline int fpe_modal_symmetrise(const Model& m, const FPE_Op& op,
                                const FPE_Mesh& g,
                                const std::vector<double>& q0,
                                FPE_ModalWork& work, double& beta0,
                                double& span, double& surv0) {
  const int M = g.M;
  work.ensure_mesh(g);
  work.off.assign(M, 0.0);
  work.logw.assign(M, 0.0);

  // t = 0 throughout: this module is only ever reached for a model whose
  // static_op() is true, and for those length(), length_prime() and
  // atil_affine() do not depend on t at all.
  const double Lb = m.length(0.0);
  const double invD = 1.0 / op.D;
  double a0 = 0.0, a1 = 0.0;
  m.atil_affine(0.0, Lb, m.length_prime(0.0), a0, a1);
  double lw = 0.0;
  for (int i = 1; i < M; ++i) {
    const double a = op.sub[i], b = op.sup[i - 1];
    if (!(a > 0.0) || !(b > 0.0) || !std::isfinite(a) || !std::isfinite(b)) {
      return FPE_MODAL_NOT_SYMMETRISABLE;
    }
    work.off[i] = std::sqrt(a * b);
    if (i == 1 || i == M - 1) {
      lw += 0.5 * std::log(a / b);
    } else {
      // pu/pv carry the face spacing and its position for BOTH mesh kinds -- on
      // a uniform mesh pu = h and pv = j h^2, so this is build_op()'s
      // arithmetic sequence written out rather than generated by recurrence.
      lw += 0.5 * (a0 * g.pu[i] + a1 * g.pv[i]) * invD +
            work.mesh_half_log_ratio[i];
    }
    work.logw[i] = lw;
  }

  int peak = 0;
  surv0 = 0.0;
  for (int i = 0; i < M; ++i) {
    surv0 += g.dx[i] * q0[i];
    if (q0[i] > q0[peak]) peak = i;
  }
  if (!(q0[peak] > 0.0)) return FPE_MODAL_SEED_UNUSABLE;

  const double centre = work.logw[peak];
  double lo = 0.0, hi = 0.0;
  for (int i = 0; i < M; ++i) {
    work.logw[i] -= centre;
    lo = std::min(lo, work.logw[i]);
    hi = std::max(hi, work.logw[i]);
  }
  span = hi - lo;
  if (!(span <= FPE_MODAL_SPAN)) return FPE_MODAL_SPAN_TOO_WIDE;

  work.v.resize(M);
  work.phis.resize(M);
  work.phif.assign(M, 0.0);
  beta0 = 0.0;
  // One exponential per cell, not two: the seed is divided by W and the
  // survivor functional multiplied by it, so the reciprocal serves both.
  for (int i = 0; i < M; ++i) {
    const double e = std::exp(work.logw[i]);
    work.phis[i] = g.dx[i] * e;
    work.v[i] = q0[i] / e;
    beta0 += work.v[i] * work.v[i];
  }
  beta0 = std::sqrt(beta0);
  if (!(beta0 > 0.0) || !std::isfinite(beta0)) return FPE_MODAL_SEED_UNUSABLE;
  const double rb = 1.0 / beta0;
  for (int i = 0; i < M; ++i) work.v[i] *= rb;
  work.phif[M - 1] =  op.D * g.cA * work.phis[M - 1] * g.rdx[M - 1];
  work.phif[M - 2] = -op.D * g.cB * work.phis[M - 2] * g.rdx[M - 2];
  for (int i = 0; i < M; ++i) {
    if (!std::isfinite(work.phis[i]) || !std::isfinite(work.v[i])) {
      return FPE_MODAL_SEED_UNUSABLE;
    }
  }
  if (!std::isfinite(work.phif[M - 1]) || !std::isfinite(work.phif[M - 2])) {
    return FPE_MODAL_SEED_UNUSABLE;
  }
  return FPE_MODAL_OK;
}

// Factorise I - gamma S, once per solve.  Same Thomas elimination as
// FPE_Tri::factor(), on the symmetrised operator.
inline void fpe_modal_factor(const FPE_Op& op, double gamma, int M,
                             FPE_ModalWork& work) {
  work.sub.resize(M);
  work.dinv.resize(M);
  work.cprime.resize(M);
  double cp = 0.0;
  for (int i = 0; i < M; ++i) {
    const double sub = (i > 0) ? -gamma * work.off[i] : 0.0;
    const double sup = (i + 1 < M) ? -gamma * work.off[i + 1] : 0.0;
    double den = (1.0 - gamma * op.diag[i]) - ((i > 0) ? sub * cp : 0.0);
    if (den == 0.0) den = std::numeric_limits<double>::min();
    const double di = 1.0 / den;
    work.sub[i] = sub;
    work.dinv[i] = di;
    cp = (i + 1 < M) ? sup * di : 0.0;
    work.cprime[i] = cp;
  }
}

// Extend the Lanczos basis from `built` vectors to `target`.  Returns the
// number actually built, which is smaller only on an invariant subspace -- and
// an invariant subspace is EXACT, so that is a reason to stop, not to refuse.
//
// No basis is stored and no reorthogonalisation is done.  The two go together:
// the only things wanted from the basis are its projections onto the two
// functionals, which are accumulated here as they are produced, so keeping
// M x m vectors around would exist solely to re-orthogonalise them.  For a
// matrix FUNCTION (as against an eigenvalue) the classical result is that loss
// of orthogonality is benign -- it happens in directions that have already
// converged -- and the ladder in fpe_modal_build() is an observable test on the
// answer, so a ghost mode that did move the answer would be caught by it and
// buy more vectors rather than being trusted.
inline int fpe_modal_extend(int M, int built, int target,
                            FPE_ModalWork& work) {
  work.w.resize(M);
  work.vprev.resize(M);
  if (built == 0) {
    std::fill(work.vprev.begin(), work.vprev.end(), 0.0);
    work.alpha.clear();
    work.beta.clear();
    work.gf.clear();
    work.gs.clear();
  }
  double bprev = work.beta.empty() ? 0.0 : work.beta.back();
  for (int j = built; j < target; ++j) {
    double pf = 0.0, ps = 0.0;
    for (int i = 0; i < M; ++i) {
      pf += work.v[i] * work.phif[i];
      ps += work.v[i] * work.phis[i];
    }
    work.gf.push_back(pf);
    work.gs.push_back(ps);

    // w = (I - gamma S)^-1 v
    double y = work.v[0] * work.dinv[0];
    work.w[0] = y;
    for (int i = 1; i < M; ++i) {
      y = (work.v[i] - work.sub[i] * y) * work.dinv[i];
      work.w[i] = y;
    }
    double z = work.w[M - 1];
    for (int i = M - 2; i >= 0; --i) {
      z = work.w[i] - work.cprime[i] * z;
      work.w[i] = z;
    }

    double al = 0.0;
    for (int i = 0; i < M; ++i) al += work.v[i] * work.w[i];
    work.alpha.push_back(al);
    double bn = 0.0;
    for (int i = 0; i < M; ++i) {
      const double u = work.w[i] - al * work.v[i] - bprev * work.vprev[i];
      work.w[i] = u;
      bn += u * u;
    }
    bn = std::sqrt(bn);
    if (!(bn > 1e-13) || !std::isfinite(bn)) return j + 1;   // invariant
    work.beta.push_back(bn);
    const double rb = 1.0 / bn;
    for (int i = 0; i < M; ++i) {
      work.vprev[i] = work.v[i];
      work.v[i] = work.w[i] * rb;
    }
    bprev = bn;
  }
  return target;
}

// Diagonalise the leading m x m block and turn it into the two modal curves.
// The leading block of a Lanczos tridiagonal IS the tridiagonal m steps would
// have produced, so a lower rung of the ladder costs one QL and no new vectors.
inline bool fpe_modal_reduce(FPE_ModalWork& work, int m, double beta0,
                             double gamma, FPE_ModalCurve& out) {
  if (m <= 0 || static_cast<int>(work.alpha.size()) < m) return false;
  work.qd.assign(work.alpha.begin(), work.alpha.begin() + m);
  work.qe.assign(m, 0.0);
  for (int j = 0; j + 1 < m; ++j) work.qe[j] = work.beta[j];
  // Rows: 0 = the flux functional, 1 = the survivor functional, 2 = e_1.
  work.qr.assign(static_cast<size_t>(3) * m, 0.0);
  for (int j = 0; j < m; ++j) {
    work.qr[j] = work.gf[j];
    work.qr[static_cast<size_t>(m) + j] = work.gs[j];
  }
  work.qr[static_cast<size_t>(2) * m] = 1.0;
  if (!fpe_tql_project(m, work.qd.data(), work.qe.data(),
                       work.qr.data(), 3)) {
    return false;
  }
  out.m = m;
  out.lam.resize(m);
  out.cf.resize(m);
  out.cs.resize(m);
  const double rg = 1.0 / gamma;
  for (int j = 0; j < m; ++j) {
    const double mu = work.qd[j];
    if (!(mu > 0.0)) return false;               // (I - gamma S) is an M-matrix
    out.lam[j] = (1.0 - 1.0 / mu) * rg;
    const double first = work.qr[static_cast<size_t>(2) * m + j];
    out.cf[j] = beta0 * work.qr[j] * first;
    out.cs[j] = beta0 * work.qr[static_cast<size_t>(m) + j] * first;
  }
  // Slowest mode first.  QL leaves the eigenvalues in whatever order deflation
  // produced; a stable ordering keeps the modal sums reproducible across
  // ladder rungs and accumulates the longest-lived contributions first.
  work.order.resize(m);
  for (int j = 0; j < m; ++j) work.order[j] = j;
  std::sort(work.order.begin(), work.order.end(),
            [&](int a, int b) { return out.lam[a] > out.lam[b]; });
  work.qd.assign(out.lam.begin(), out.lam.end());
  work.qe.assign(out.cf.begin(), out.cf.end());
  work.qr.assign(out.cs.begin(), out.cs.end());
  for (int j = 0; j < m; ++j) {
    const int o = work.order[j];
    out.lam[j] = work.qd[o];
    out.cf[j] = work.qe[o];
    out.cs[j] = work.qr[o];
  }
  return true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Build a propagator for one fixed-boundary model, choosing m by convergence.
//
// `probe` are the times, MEASURED FROM t0, that the ladder is judged at.  They
// should be the times the caller is actually going to ask for -- a modal form
// converged on the response times of a data set is converged where it matters,
// and the deep leading edge, where the expansion is a cancelling sum and no
// dimension settles it, is not something to spend vectors on.
//
// The two functionals are judged in different currencies, for the same reason
// they are in rlf_build_modes(): the density relatively but against a
// denominator floored at FPE_MODAL_FLOOR of the peak, so the leading edge stays
// in the test without the deep tail dominating it; the survivor absolutely in
// logs, because log S is what the omission and truncation paths read and it
// decays smoothly rather than collapsing into noise.
//
// A rung is refused outright, whatever its movement, if it produces a
// non-positive density or survivor at a probe where the density is above the
// floor, or if the sum that produced it amplified by more than
// FPE_MODAL_CANCEL.  Those are the two ways a modal expansion fails that more
// vectors do not fix.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// The rung test, factored out so the scalar and lane-batched drivers apply
// exactly the same acceptance rule.  Holds the previous rung's values.
struct FPE_ModalLadder {
  std::vector<double> prev_f, prev_s, fv, sv, av;
  bool have_prev = false;

  void reset(size_t np) {
    prev_f.assign(np, 0.0);
    prev_s.assign(np, 0.0);
    fv.assign(np, 0.0);
    sv.assign(np, 0.0);
    av.assign(np, 0.0);
    have_prev = false;
  }

  enum Verdict { ACCEPT, CLIMB, REFUSE };

  // `invariant` says the Lanczos space closed on itself (the rung is exact),
  // `at_cap` that there is nothing above this rung to climb to.
  Verdict test(const FPE_ModalCurve& rung, const std::vector<double>& probe,
               bool invariant, bool at_cap, int& status) {
    const size_t np = probe.size();
    double peak = 0.0;
    for (size_t k = 0; k < np; ++k) {
      rung.probe(probe[k], fv[k], sv[k], av[k]);
      peak = std::max(peak, fv[k]);
    }
    const double floor_f = std::max(peak * FPE_MODAL_FLOOR, 1e-300);

    // The density must be STRICTLY POSITIVE at every probe, not merely where it
    // carries the peak.  A probe is a time the caller is going to ask for, and
    // safe_log() of a non-positive density is LOG_FLOOR: a trial whose true log
    // density is -35 would enter the likelihood at -700.  The amplification --
    // sum|term| / |sum| -- is checked at every probe for the same reason: a
    // value can be positive and still have lost every digit it had.  Both are
    // reasons to CLIMB, not to give up: what fails at a low dimension is
    // usually the reconstruction rather than the transform, and the
    // cancellation pattern does change with m.  Measured, starting the ladder
    // at 10 rather than 16 turned 16% of an ordinary parameter net from
    // accepted into refused purely on first-rung failures.
    bool usable = true;
    for (size_t k = 0; k < np && usable; ++k) {
      if (!std::isfinite(fv[k]) || !std::isfinite(sv[k]) ||
          !(sv[k] > 0.0) || sv[k] > 1.0 + 1e-9 || !(fv[k] > 0.0)) usable = false;
      else if (av[k] > FPE_MODAL_CANCEL) usable = false;
    }

    if (usable) {
      double moved = HUGE_VAL;
      if (have_prev) {
        moved = 0.0;
        double moved_log = 0.0;
        for (size_t k = 0; k < np; ++k) {
          const double den = std::max(std::abs(fv[k]), floor_f);
          moved = std::max(moved, std::abs(fv[k] - prev_f[k]) / den);
          moved = std::max(moved, std::abs(std::log(sv[k] / prev_s[k])));
          if (prev_f[k] > 0.0)
            moved_log = std::max(moved_log, std::abs(std::log(fv[k] / prev_f[k])));
        }
        // Scale the log movement onto the same axis as the floored test so one
        // comparison decides both.
        moved = std::max(moved, moved_log * (FPE_MODAL_TOL / FPE_MODAL_TOL_LOG));
      }
      prev_f = fv;
      prev_s = sv;
      have_prev = true;
      if (invariant || moved < FPE_MODAL_TOL) { status = FPE_MODAL_OK; return ACCEPT; }
      if (at_cap) { status = FPE_MODAL_NOT_CONVERGED; return REFUSE; }
      return CLIMB;
    }

    have_prev = false;   // nothing usable to compare the next rung against
    if (invariant || at_cap) { status = FPE_MODAL_CANCELLING; return REFUSE; }
    return CLIMB;
  }
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Build a propagator for one fixed-boundary model, choosing m by convergence.
//
// `probe` are the times, MEASURED FROM t0, that the ladder is judged at.  They
// should be the times the caller is actually going to ask for -- a modal form
// converged on the response times of a data set is converged where it matters,
// and the deep leading edge, where the expansion is a cancelling sum and no
// dimension settles it, is not something to spend vectors on.
//
// The two functionals are judged in different currencies, for the same reason
// they are in rlf_build_modes(): the density relatively but against a
// denominator floored at FPE_MODAL_FLOOR of the peak, so the leading edge stays
// in the test without the deep tail dominating it; the survivor absolutely in
// logs, because log S is what the omission and truncation paths read and it
// decays smoothly rather than collapsing into noise.
//
// This is the scalar driver.  fperace::fpe_modal_batch_ou() runs the same
// ladder with the Lanczos itself interleaved across SIMD lanes; both go through
// FPE_ModalLadder, so they accept and refuse identically.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template <class Model>
inline FPE_ModalCurve fpe_modal_build(const Model& model, const FPE_Op& op,
                                      const FPE_Mesh& g,
                                      const std::vector<double>& q0,
                                      double t0, double horizon,
                                      const std::vector<double>& probe,
                                      FPE_ModalWork& work) {
  FPE_ModalCurve out;
  out.t0 = t0;
  const int M = g.M;
  double beta0 = 0.0;
  const int sym = fpe_modal_symmetrise(model, op, g, q0, work, beta0, out.span,
                                       out.surv0);
  if (sym != FPE_MODAL_OK) { out.status = sym; return out; }

  const double gamma = FPE_MODAL_SHIFT * std::max(horizon, 1e-9);
  fpe_modal_factor(op, gamma, M, work);

  const int m_cap = std::min(FPE_MODAL_MAX, M);
  FPE_ModalLadder ladder;
  ladder.reset(probe.size());
  int built = 0;
  int m = std::min(FPE_MODAL_MIN, m_cap);
  FPE_ModalCurve rung;
  rung.t0 = t0;
  rung.surv0 = out.surv0;
  rung.span = out.span;

  while (true) {
    const int grown = fpe_modal_extend(M, built, m, work);
    const bool invariant = grown < m;
    built = grown;
    m = grown;
    if (!fpe_modal_reduce(work, m, beta0, gamma, rung)) {
      out.status = FPE_MODAL_QL_FAILED;
      return out;
    }
    int status = FPE_MODAL_NOT_CONVERGED;
    const FPE_ModalLadder::Verdict v =
      ladder.test(rung, probe, invariant, m >= m_cap, status);
    if (v == FPE_ModalLadder::ACCEPT) { rung.status = FPE_MODAL_OK; return rung; }
    if (v == FPE_ModalLadder::REFUSE) {
      out.status = status;
      out.span = rung.span;
      return out;
    }
    m = std::min(fpe_modal_next_rung(m), m_cap);
  }
}

} // namespace fpe

#endif // fpe_modal_h
