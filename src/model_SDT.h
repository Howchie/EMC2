#ifndef EMC2_MODEL_SDT_H
#define EMC2_MODEL_SDT_H

#include <Rcpp.h>
#include "wald_functions.h"
#include <algorithm>
#include <cmath>

/*
 * Hierarchical Unequal-Variance Signal Detection (hUVSD)
 * Based on Lages (2024), Psychonomic Bulletin & Review.
 * 
 * Parameters:
 *   d  - Sensitivity (distance between signal and noise means)
 *   c  - Bias (deviation from the midpoint between means)
 *   sd - Standard deviation of the signal distribution (noise SD is fixed at 1)
 * 
 * Identifiability:
 *   Lages (2024) shows that with hierarchical modeling, sd is identifiable 
 *   from binary data by exploiting variability across participants.
 */

inline double log_likelihood_huvsd_single(double d, double c, double sd, 
                                          bool is_signal, bool chosen_yes, 
                                          double min_ll) {
  // Equation 3 from the paper:
  // theta_h = Phi((+0.5 * d - c) / sd)
  // theta_f = Phi((-0.5 * d - c) / 1.0)
  
  double mu = is_signal ? (0.5 * d) : (-0.5 * d);
  double sigma = is_signal ? sd : 1.0;
  
  // P(Yes) = 1 - Phi((c - mu) / sigma) = Phi((mu - c) / sigma)
  double z = (mu - c) / sigma;
  
  double log_p;
  if (chosen_yes) {
    log_p = pnorm_std(z, true, true);
  } else {
    log_p = pnorm_std(z, false, true);
  }
  
  return std::max(log_p, min_ll);
}

#endif
