#ifndef EMC2_RT_BINS_H
#define EMC2_RT_BINS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include "utility_functions.h"
#include "composite_functions.h"
#include "contaminant_mixture.h"
#include "gl_quad.h"

inline double rt_bin_width(const Rcpp::DataFrame& data) {
  if (!data.hasAttribute("rt_resolution")) return 0.0;
  const double h = Rcpp::as<double>(data.attr("rt_resolution"));
  return R_FINITE(h) && h > 0.0 ? h : 0.0;
}

struct RtBin {
  double lo, hi, h, edge;
};

// The stored RT is inside its grid bin. Keep the nominal width as the
// denominator even when observation bounds remove part of that bin.
inline RtBin rt_bin(double rt, double h, double LT = 0.0,
                    double UT = R_PosInf, double LC = 0.0,
                    double UC = R_PosInf) {
  const double edge = h * std::floor(rt / h);
  return {std::max(edge, std::max(LT, LC)),
          std::min(edge + h, std::min(UT, UC)), h, edge};
}

// Integrate unfloored log densities. The square-root coordinate removes the
// common onset singularity. A fixed low-order Gauss-Legendre rule in
// z = sqrt(t - onset) integrates polynomials up to degree 2n-1; floor-aware
// early exit skips refining negligible mass below the likelihood floor.
template <typename LogDensity>
double rt_bin_log_panel(double lo, double hi, const LogDensity& f,
                       double log_absolute_tolerance = -26.0, int depth = 0) {
  if (!(hi > lo)) return R_NegInf;
  auto rule = [&](int n) {
    const GLRule& g = gl_get_rule(n);
    double sum = R_NegInf;
    const double half = (hi - lo) / 2.0;
    for (int k = 0; k < n; ++k) {
      const double x = lo + half * (g.x[k] + 1.0);
      const double fx = f(x);
      if (R_FINITE(fx)) {
        sum = log_sum_exp(sum, std::log(half * g.w[k]) + fx);
      }
    }
    return sum;
  };
  const double a = rule(4), b = rule(8);
  // Floor-aware early exit: negligible probability mass well below the
  // likelihood floor (~1e-10) is clamped downstream.
  if (b < -26.0 && a < -26.0) return b;
  const double error = a > b ? log_diff_exp(a, b) : log_diff_exp(b, a);
  if (depth >= 3 || (a == R_NegInf && b == R_NegInf) ||
      error < log_absolute_tolerance ||
      (R_FINITE(a) && R_FINITE(b) && std::fabs(a - b) < 1e-6)) return b;
  const double mid = lo + (hi - lo) / 2.0;
  return log_sum_exp(rt_bin_log_panel(lo, mid, f, log_absolute_tolerance - M_LN2, depth + 1),
                     rt_bin_log_panel(mid, hi, f, log_absolute_tolerance - M_LN2, depth + 1));
}

inline double rt_bin_quadrature_weight(const RtBin& bin, double onset) {
  if (onset >= bin.hi) return 0.0;
  // If onset is inside the nominal bin, support is truncated: full quadrature weight.
  if (onset >= bin.edge) return 1.0;
  // If onset is below the bin: support is complete.
  // Smoothly blend to midpoint evaluation over at most 1 bin width (C1 cubic Hermite).
  const double u = (bin.edge - onset) / bin.h;
  if (u >= 1.0) return 0.0;
  const double x = 1.0 - u;
  return x * x * (3.0 - 2.0 * x);
}

inline bool rt_bin_near_support(const RtBin& bin, double onset) {
  return rt_bin_quadrature_weight(bin, onset) > 0.0;
}

// When a zero-start guess is present, t0 is an internal decision onset, not
// the onset of the complete density. Fade out beyond the bin too, since the
// remaining guess density must not jump from its integral to its midpoint.
inline double rt_bin_internal_onset_weight(const RtBin& bin, double onset) {
  if (onset < bin.hi) return rt_bin_quadrature_weight(bin, onset);
  const double u = (onset - bin.hi) / bin.h;
  if (u >= 1.0) return 0.0;
  const double x = 1.0 - u;
  return x * x * (3.0 - 2.0 * x);
}

// A loser's onset influences a generating cause through its survivor. Fade
// that influence at the upper edge too; its own density still gets the full
// onset rule in rt_bin_log_average().
inline double rt_bin_survivor_weight(const RtBin& bin, double onset) {
  const double qw = rt_bin_quadrature_weight(bin, onset);
  if (qw <= 0.0 || onset >= bin.hi) return 0.0;
  const double u = std::max(0.0, (bin.hi - onset) / bin.h);
  if (u >= 1.0) return qw;
  return qw * (u * u * (3.0 - 2.0 * u));
}

template <typename LogDensity>
double rt_bin_log_average(const RtBin& bin, double onset, const LogDensity& f,
                          double refinement = 0.0) {
  const double lo = std::max(bin.lo, onset);
  if (!(bin.hi > lo)) return R_NegInf;
  const double weight = std::max(refinement, rt_bin_quadrature_weight(bin, onset));
  auto midpoint = [&]() { return f(bin.lo + (bin.hi - bin.lo) / 2.0) +
      std::log((bin.hi - bin.lo) / bin.h); };
  if (weight == 0.0) return midpoint();
  const double zl = std::sqrt(lo - onset), zh = std::sqrt(bin.hi - onset);
  // Floor-aware early exit tolerance below min_ll (~1e-10)
  const double integral = rt_bin_log_panel(zl, zh, [&](double z) {
    return std::log(2.0 * z) + f(onset + z * z);
  }, std::log(bin.h) + std::log(1e-12)) - std::log(bin.h);
  return weight >= 1.0 ? integral : log_sum_exp(
    std::log(weight) + integral, std::log1p(-weight) + midpoint());
}

template <typename LogDensity>
double rt_bin_log_average_split(const RtBin& bin, double onset,
                                const LogDensity& f, double refinement,
                                const std::vector<double>& transitions) {
  const double lo = std::max(bin.lo, onset);
  if (!(bin.hi > lo)) return R_NegInf;
  const double weight = std::max(refinement, rt_bin_quadrature_weight(bin, onset));
  auto midpoint = [&]() { return f(bin.lo + (bin.hi - bin.lo) / 2.0) +
      std::log((bin.hi - bin.lo) / bin.h); };
  if (weight == 0.0) return midpoint();

  std::vector<double> zcuts{std::sqrt(lo - onset)};
  for (double t : transitions)
    if (R_FINITE(t) && t > lo && t < bin.hi)
      zcuts.push_back(std::sqrt(t - onset));
  zcuts.push_back(std::sqrt(bin.hi - onset));
  std::sort(zcuts.begin(), zcuts.end());
  zcuts.erase(std::unique(zcuts.begin(), zcuts.end(), [](double a, double b) {
    return std::fabs(a - b) < 1e-14;
  }), zcuts.end());

  double log_integral = R_NegInf;
  const double panel_tolerance = std::log(bin.h) + std::log(1e-12) -
    std::log(static_cast<double>(std::max<std::size_t>(1, zcuts.size() - 1)));
  for (size_t i = 1; i < zcuts.size(); ++i) {
    const double zl = zcuts[i - 1], zh = zcuts[i];
    const double panel = rt_bin_log_panel(zl, zh, [&](double z) {
      return std::log(2.0 * z) + f(onset + z * z);
    }, panel_tolerance);
    log_integral = log_sum_exp(log_integral, panel);
  }
  const double integral = log_integral - std::log(bin.h);
  return weight >= 1.0 ? integral : log_sum_exp(
    std::log(weight) + integral, std::log1p(-weight) + midpoint());
}

// Data-fixed intervals for the vectorised independent-race path. The normal
// kernel pass still evaluates just once at the stored midpoint; only trials
// near an accumulator onset need the scalar quadrature fallback.
struct RtBinCache {
  double h = 0.0;
  std::vector<RtBin> bins;
  std::vector<double> log_fraction;
  // Over the bins of finite positive rts: whether any is empty, whether any
  // is near an onset at zero, and the smallest onset that could reach one.
  bool any_empty = false;
  bool any_near_zero = false;
  double min_far_below = R_PosInf;
  void init(const Rcpp::DataFrame& data, int n_lR) {
    bins.clear();
    log_fraction.clear();
    any_empty = false;
    any_near_zero = false;
    min_far_below = R_PosInf;
    h = rt_bin_width(data);
    if (h == 0.0) return;
    Rcpp::NumericVector rt = data["rt"];
    auto column = [&](const char* name) -> const double* {
      return data.containsElementNamed(name) ? Rcpp::NumericVector(data[name]).begin() : nullptr;
    };
    const double* LT = column("LT"), *UT = column("UT");
    const double* LC = column("LC"), *UC = column("UC");
    bins.reserve(rt.size() / n_lR);
    for (int row = 0; row < rt.size(); row += n_lR) {
      bins.push_back(R_FINITE(rt[row]) ?
        rt_bin(rt[row], h, LT ? LT[row] : 0.0, UT ? UT[row] : R_PosInf,
               LC ? LC[row] : 0.0, UC ? UC[row] : R_PosInf) :
        RtBin{0.0, 0.0, h, 0.0});
      const RtBin& bin = bins.back();
      log_fraction.push_back(bin.hi > bin.lo ? std::log((bin.hi - bin.lo) / h) : R_NegInf);
      if (R_FINITE(rt[row]) && rt[row] > 0.0) {
        if (!(bin.hi > bin.lo)) any_empty = true;
        if (rt_bin_near_support(bin, 0.0)) any_near_zero = true;
        min_far_below = std::min(min_far_below, bin.edge - 2.0 * bin.h);
      }
    }
  }
  bool active() const { return h > 0.0; }
};

inline GuessKernel rt_bin_guess(const GuessKernel& guess, const RtBin& bin) {
  GuessKernel result = guess;
  const double lo = std::max(bin.lo, guess.LG);
  const double hi = std::min(bin.hi, guess.UG);
  const double fraction = hi > lo ? std::log((hi - lo) / bin.h) : R_NegInf;
  result.log_g += fraction;
  result.log_g_any += fraction;
  return result;
}

// rt_bin_guess into a caller-owned kernel that already holds a copy of
// `guess`: only the two densities change per bin, so per-trial loops avoid
// copying the kernel and its shared response mask.
inline const GuessKernel& rt_bin_guess_into(const GuessKernel& guess, const RtBin& bin,
                                            GuessKernel& out) {
  const double lo = std::max(bin.lo, guess.LG);
  const double hi = std::min(bin.hi, guess.UG);
  const double fraction = hi > lo ? std::log((hi - lo) / bin.h) : R_NegInf;
  out.log_g = guess.log_g + fraction;
  out.log_g_any = guess.log_g_any + fraction;
  return out;
}

#endif
