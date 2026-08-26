import re

with open('src/wald_functions.h', 'r') as f:
    code = f.read()

# I will systematically add overloads.

# 1. log_split_lognormal_density
density_orig = """inline double log_split_lognormal_density(double v, double mu, double sigma,
                                          double delta) {
  if (delta == 0.0) return dlnorm_std(v, mu, sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;"""

density_overload = """inline double log_split_lognormal_density(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return R_NegInf;
  const double x = std::log(v);
  const double s = (x <= h.c) ? h.sL : h.sR;
  const double mass = (x <= h.c) ? h.a : (1.0 - h.a);
  const double z = (x - h.c) / s;
  return std::log(2.0 * mass) - std::log(v * s * M_SQRT2PI) - 0.5 * z * z;
}

inline double log_split_lognormal_density(double v, double mu, double sigma,
                                          double delta) {
  if (delta == 0.0) return dlnorm_std(v, mu, sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_density(v, h);
}"""

# 2. log_split_lognormal_cdf
cdf_orig = """inline double log_split_lognormal_cdf(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return pnorm_log_direct((std::log(v) - mu) / sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  const double x = std::log(v);
  if (x <= h.c) {
    const double z = (x - h.c) / h.sL;
    return std::fmin(std::log(2.0 * h.a) + pnorm_log_direct(z, true), 0.0);
  }
  const double ls = log_split_lognormal_survivor(v, mu, sigma, delta);
  return (ls < 0.0) ? log1m_exp(ls) : R_NegInf;
}"""

cdf_overload = """inline double log_split_lognormal_survivor(double v, const split_lognormal_shape& h);

inline double log_split_lognormal_cdf(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return R_NegInf;
  const double x = std::log(v);
  if (x <= h.c) {
    const double z = (x - h.c) / h.sL;
    return std::fmin(std::log(2.0 * h.a) + pnorm_log_direct(z, true), 0.0);
  }
  const double z = (x - h.c) / h.sR;
  const double lp = std::log(2.0 * (1.0 - h.a)) + pnorm_log_direct(z, false);
  return (lp > R_NegInf) ? log1p_exp(std::log(2.0 * h.a - 1.0) - lp) + lp : R_NegInf;
}

inline double log_split_lognormal_cdf(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return pnorm_log_direct((std::log(v) - mu) / sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_cdf(v, h);
}"""

# 3. log_split_lognormal_survivor
surv_orig = """inline double log_split_lognormal_survivor(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return lnorm_log_surv_std(v, mu, sigma);
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);
}"""

surv_overload = """inline double split_lognormal_log_moment_tail(double v, const split_lognormal_shape& h, double r) {
  const double lpL = split_lognormal_log_side_tail(std::log(v), r, h, true);
  const double lpR = split_lognormal_log_side_tail(std::log(v), r, h, false);
  return log_sum_exp(lpL, lpR);
}

inline double split_lognormal_log_moment_tail(double v, double mu,
                                              double sigma, double delta,
                                              double r) {
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return split_lognormal_log_moment_tail(v, h, r);
}

inline double log_split_lognormal_survivor(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, h, 0.0);
}

inline double log_split_lognormal_survivor(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return lnorm_log_surv_std(v, mu, sigma);
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);
}"""

# Need to replace the moment tail function definition as well.
# We will do it properly with sed or a python replacer.

# 4. log_split_lognormal_first_partial_moment
first_orig = """inline double log_split_lognormal_first_partial_moment(double v, double mu,
                                                       double sigma,
                                                       double delta) {
  if (delta == 0.0) {
    if (!(v > 0.0)) return mu + 0.5 * sigma * sigma;
    const double x = (std::log(v) - mu) / sigma;
    return mu + 0.5 * sigma * sigma +
      pnorm_log_direct(x - sigma, false);
  }
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 1.0);
}"""

first_overload = """inline double log_split_lognormal_first_partial_moment(double v, const split_lognormal_shape& h) {
  return split_lognormal_log_moment_tail(v, h, 1.0);
}

inline double log_split_lognormal_first_partial_moment(double v, double mu,
                                                       double sigma,
                                                       double delta) {
  if (delta == 0.0) {
    if (!(v > 0.0)) return mu + 0.5 * sigma * sigma;
    const double x = (std::log(v) - mu) / sigma;
    return mu + 0.5 * sigma * sigma +
      pnorm_log_direct(x - sigma, false);
  }
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 1.0);
}"""

# 5. log_split_lognormal_stoploss
stoploss_orig = """inline double log_split_lognormal_stoploss(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return log_lognormal_stoploss(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  const double lm = split_lognormal_log_total_moment(mu, sigma, delta, 1.0);
  const double x = std::log(v);
  if (x > h.c) {
    const double lz = (x - h.c) / h.sR;
    if (lz > 1.5) {
      const double lr = log_mills_gap(lz, h.sR);
      if (emc2_isfinite(lr))
        return std::log(v) + std::log(2.0 * (1.0 - h.a)) +
          pnorm_log_direct(lz, false) + lr;
    }
  }
  const double lpL = split_lognormal_log_stoploss_side(x, h, true);
  const double lpR = split_lognormal_log_stoploss_side(x, h, false);
  return log_sum_exp(lpL, lpR);
}"""

stoploss_overload = """inline double split_lognormal_log_total_moment(const split_lognormal_shape& h, double r) {
  return split_lognormal_log_moment_tail(0.0, h, r);
}

inline double log_split_lognormal_stoploss(double v, const split_lognormal_shape& h) {
  const double x = std::log(v);
  if (x > h.c) {
    const double lz = (x - h.c) / h.sR;
    if (lz > 1.5) {
      const double lr = log_mills_gap(lz, h.sR);
      if (emc2_isfinite(lr))
        return std::log(v) + std::log(2.0 * (1.0 - h.a)) +
          pnorm_log_direct(lz, false) + lr;
    }
  }
  const double lpL = split_lognormal_log_stoploss_side(x, h, true);
  const double lpR = split_lognormal_log_stoploss_side(x, h, false);
  return log_sum_exp(lpL, lpR);
}

inline double log_split_lognormal_stoploss(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return log_lognormal_stoploss(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_stoploss(v, h);
}"""

# 6. log_split_lognormal_put
put_orig = """inline double log_split_lognormal_put(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return log_lognormal_put(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  const double x = std::log(v);
  const double lpL = split_lognormal_log_put_side(x, h, true);
  const double lpR = split_lognormal_log_put_side(x, h, false);
  return log_sum_exp(lpL, lpR);
}"""

put_overload = """inline double log_split_lognormal_put(double v, const split_lognormal_shape& h) {
  const double x = std::log(v);
  const double lpL = split_lognormal_log_put_side(x, h, true);
  const double lpR = split_lognormal_log_put_side(x, h, false);
  return log_sum_exp(lpL, lpR);
}

inline double log_split_lognormal_put(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return log_lognormal_put(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_put(v, h);
}"""

# 7. log_split_lognormal_power_stoploss
power_orig = """inline double log_split_lognormal_power_stoploss(double v, double mu,
                                                 double sigma, double delta,
                                                 double m) {
  if (delta == 0.0) return log_lognormal_power_stoploss(v, mu, sigma, m);
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  if (std::fabs(m) <= 1e-14) {
    // E[(Y-log v)+] is assembled from the two truncated first moments.
    split_lognormal_shape h;
    if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
    const double x = std::log(v);
    auto yside = [&](bool left) {
      const double s = left ? h.sL : h.sR;
      const double mass = left ? h.a : (1.0 - h.a);
      const double z = (x - h.c) / s;
      if (left) {
        const double p = pnorm_std(0.0, false, false) -
          pnorm_std(z, false, false);
        const double q = dnormP(z) - dnormP(0.0);
        const double val = (h.c - x) * p + s * q;
        return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
      }
      const double zz = (x > h.c) ? z : 0.0;
      const double p = pnorm_std(zz, false, false);
      const double q = dnormP(zz);
      const double val = (h.c - x) * p + s * q;
      return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
    };
    return log_sum_exp(yside(true), yside(false));
  }
  const double lm = split_lognormal_log_moment_tail(v, mu, sigma, delta, -m);
  const double ls = log_split_lognormal_survivor(v, mu, sigma, delta);
  const signed_log out = (m > 0.0)
    ? signed_log_sub(make_signed_log(std::log(v) * (-m) + ls, 1),
                     make_signed_log(lm, 1))
    : signed_log_sub(make_signed_log(lm, 1),
                     make_signed_log(std::log(v) * (-m) + ls, 1));
  return out.sign > 0 ? out.log_abs - std::log(std::fabs(m)) : R_NegInf;
}"""

power_overload = """inline double log_split_lognormal_power_stoploss(double v, const split_lognormal_shape& h, double m) {
  if (std::fabs(m) <= 1e-14) {
    const double x = std::log(v);
    auto yside = [&](bool left) {
      const double s = left ? h.sL : h.sR;
      const double mass = left ? h.a : (1.0 - h.a);
      const double z = (x - h.c) / s;
      if (left) {
        const double p = pnorm_std(0.0, false, false) -
          pnorm_std(z, false, false);
        const double q = dnormP(z) - dnormP(0.0);
        const double val = (h.c - x) * p + s * q;
        return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
      }
      const double zz = (x > h.c) ? z : 0.0;
      const double p = pnorm_std(zz, false, false);
      const double q = dnormP(zz);
      const double val = (h.c - x) * p + s * q;
      return val > 0.0 ? std::log(2.0 * mass * val) : R_NegInf;
    };
    return log_sum_exp(yside(true), yside(false));
  }
  const double lm = split_lognormal_log_moment_tail(v, h, -m);
  const double ls = log_split_lognormal_survivor(v, h);
  const signed_log out = (m > 0.0)
    ? signed_log_sub(make_signed_log(std::log(v) * (-m) + ls, 1),
                     make_signed_log(lm, 1))
    : signed_log_sub(make_signed_log(lm, 1),
                     make_signed_log(std::log(v) * (-m) + ls, 1));
  return out.sign > 0 ? out.log_abs - std::log(std::fabs(m)) : R_NegInf;
}

inline double log_split_lognormal_power_stoploss(double v, double mu,
                                                 double sigma, double delta,
                                                 double m) {
  if (delta == 0.0) return log_lognormal_power_stoploss(v, mu, sigma, m);
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_power_stoploss(v, h, m);
}"""

# Now write the replacements. First need to extract the original blocks from the file using regex,
# since the file content might be slightly different.

# Replace moment tail
code = re.sub(r'inline double split_lognormal_log_moment_tail\(double v, double mu,.*?return log_sum_exp\(lpL, lpR\);\n}', surv_overload, code, flags=re.DOTALL)
# The regex above will match `split_lognormal_log_moment_tail` down to `return log_sum_exp(lpL, lpR);\n}`

# The survivor is inside that chunk or separate?
# Let's just do it directly using replace_file_content if we can.
