#ifndef EMC2_GH_QUAD_H
#define EMC2_GH_QUAD_H

// Fixed Gauss-Hermite rule for integrals against a standard normal latent
// factor.  The nodes/weights below are for exp(-x^2); callers use
// z = sqrt(2) * x and weight / sqrt(pi).

#include <array>
#include <cstddef>
#include <cmath>

struct GHRule10 {
  std::array<double, 10> x;
  std::array<double, 10> w;
};

inline const GHRule10& gh_rule10() {
  static const GHRule10 rule = {
    {{-3.4361591188377376033, -2.5327316742327897964,
      -1.7566836492998817735, -1.0366108297895136542,
      -0.3429013272237046088,  0.3429013272237046088,
       1.0366108297895136542,  1.7566836492998817735,
       2.5327316742327897964,  3.4361591188377376033}},
    {{0.00000764043285523262063, 0.00134364574678123269,
      0.03387439445548106314, 0.24013861108231468642,
      0.61086263373532579878, 0.61086263373532579878,
      0.24013861108231468642, 0.03387439445548106314,
      0.00134364574678123269, 0.00000764043285523262063}}
  };
  return rule;
}

inline double gh_standard_normal_weight10(int i) {
  return gh_rule10().w[static_cast<size_t>(i)] /
    std::sqrt(std::acos(-1.0));
}

#endif
