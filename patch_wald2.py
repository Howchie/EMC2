import re

with open('src/wald_functions.h', 'r') as f:
    text = f.read()

def replace_regex(text, pattern, replacement):
    new_text, count = re.subn(pattern, replacement, text, flags=re.DOTALL)
    if count == 0:
        print("Failed to match:", pattern[:50])
    return new_text

# 1. split_lognormal_log_moment_tail
p1 = r'inline double split_lognormal_log_moment_tail\(double v, double mu,.*?return log_sum_exp\(split_lognormal_log_side_tail\(x, r, h, true\),\s*split_lognormal_log_side_tail\(x, r, h, false\)\);\n\}'
r1 = """inline double split_lognormal_log_moment_tail(double v, const split_lognormal_shape& h, double r) {
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
}"""
text = replace_regex(text, p1, r1)

# 2. split_lognormal_log_total_moment
p2 = r'inline double split_lognormal_log_total_moment\(double mu, double sigma,.*?return split_lognormal_log_moment_tail\(0\.0, mu, sigma, delta, r\);\n\}'
r2 = """inline double split_lognormal_log_total_moment(const split_lognormal_shape& h, double r) {
  return split_lognormal_log_moment_tail(0.0, h, r);
}

inline double split_lognormal_log_total_moment(double mu, double sigma,
                                               double delta, double r) {
  return split_lognormal_log_moment_tail(0.0, mu, sigma, delta, r);
}"""
text = replace_regex(text, p2, r2)

# 3. log_split_lognormal_survivor
p3 = r'inline double log_split_lognormal_survivor\(double v, double mu, double sigma,.*?return split_lognormal_log_moment_tail\(v, mu, sigma, delta, 0\.0\);\n\}'
r3 = """inline double log_split_lognormal_survivor(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, h, 0.0);
}

inline double log_split_lognormal_survivor(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return lnorm_log_surv_std(v, mu, sigma);
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);
}"""
text = replace_regex(text, p3, r3)

# 4. log_split_lognormal_first_partial_moment
p4 = r'inline double log_split_lognormal_first_partial_moment\(double v, double mu,.*?return split_lognormal_log_moment_tail\(v, mu, sigma, delta, 1\.0\);\n\}'
r4 = """inline double log_split_lognormal_first_partial_moment(double v, const split_lognormal_shape& h) {
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
text = replace_regex(text, p4, r4)

# 5. log_split_lognormal_stoploss
p5 = r'inline double log_split_lognormal_stoploss\(double v, double mu, double sigma,.*?return ll;\n\}'
r5 = """inline double log_split_lognormal_stoploss(double v, const split_lognormal_shape& h) {
  const double lm = split_lognormal_log_total_moment(h, 1.0);
  if (!(v > 0.0)) {
    if (v == 0.0) return lm;
    return log_sum_exp(lm, std::log(-v));
  }
  const double x = std::log(v);
  double ll = R_NegInf;
  if (x > h.c) {
    const double z = (x - h.c) / h.sR;
    const double r1 = mills_ratio_std(z - h.sR);
    const double r0 = mills_ratio_std(z);
    const double gap = r1 - r0;
    if (gap > 1e-12 * r1)
      ll = std::log(2.0 * (1.0 - h.a)) + std::log(v) +
        log_phi_std(z) + std::log(gap);
  }
  if (!(ll > R_NegInf))
    ll = log_sum_exp(split_lognormal_log_stoploss_side(x, v, h, true),
                     split_lognormal_log_stoploss_side(x, v, h, false));
  return ll;
}

inline double log_split_lognormal_stoploss(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return log_lognormal_stoploss(v, mu, sigma);
  if (!(sigma > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_stoploss(v, h);
}"""
text = replace_regex(text, p5, r5)

# 6. log_split_lognormal_put
p6 = r'inline double log_split_lognormal_put\(double v, double mu, double sigma,.*?return log_sum_exp\(side_put\(true\), side_put\(false\)\);\n\}'
r6 = """inline double log_split_lognormal_put(double v, const split_lognormal_shape& h) {
  const double x = std::log(v);
  auto side_put = [&](bool left) {
    const double lp = split_lognormal_log_side_lower(x, 0.0, h, left);
    const double lm = split_lognormal_log_side_lower(x, 1.0, h, left);
    if (!(lp > R_NegInf) || !(lm > R_NegInf)) return R_NegInf;
    const signed_log out = signed_log_sub(
      make_signed_log(std::log(v) + lp, 1), make_signed_log(lm, 1));
    return out.sign > 0 ? out.log_abs : R_NegInf;
  };
  return log_sum_exp(side_put(true), side_put(false));
}

inline double log_split_lognormal_put(double v, double mu, double sigma,
                                      double delta) {
  if (delta == 0.0) return log_lognormal_put(v, mu, sigma);
  if (!(sigma > 0.0) || !(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_put(v, h);
}"""
text = replace_regex(text, p6, r6)

# 7. log_split_lognormal_power_stoploss
p7 = r'inline double log_split_lognormal_power_stoploss\(double v, double mu,.*?return out\.sign > 0 \? out\.log_abs - std::log\(std::fabs\(m\)\) : R_NegInf;\n\}'
r7 = """inline double log_split_lognormal_power_stoploss(double v, const split_lognormal_shape& h, double m) {
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
text = replace_regex(text, p7, r7)

# 8. log_split_lognormal_density
p8 = r'inline double log_split_lognormal_density\(double v, double mu, double sigma,.*?0\.5 \* \(x - h\.c\) \* \(x - h\.c\) / \(s \* s\) - LOG_SQRT_2PI;\n\}'
r8 = """inline double log_split_lognormal_density(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return R_NegInf;
  const double x = std::log(v);
  const double s = (x <= h.c) ? h.sL : h.sR;
  const double mass = (x <= h.c) ? h.a : (1.0 - h.a);
  return std::log(2.0 * mass) - std::log(s) - std::log(v) -
    0.5 * (x - h.c) * (x - h.c) / (s * s) - LOG_SQRT_2PI;
}

inline double log_split_lognormal_density(double v, double mu, double sigma,
                                          double delta) {
  if (delta == 0.0) return dlnorm_std(v, mu, sigma, true);
  if (!(v > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (!split_lognormal_shape_params(mu, sigma, delta, h)) return R_NegInf;
  return log_split_lognormal_density(v, h);
}"""
text = replace_regex(text, p8, r8)

# 9. log_split_lognormal_cdf
p9 = r'inline double log_split_lognormal_cdf\(double v, double mu, double sigma,.*?return \(ls < 0\.0\) \? log1m_exp\(ls\) : R_NegInf;\n\}'
r9 = """inline double log_split_lognormal_cdf(double v, const split_lognormal_shape& h) {
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
text = replace_regex(text, p9, r9)

with open('src/wald_functions.h', 'w') as f:
    f.write(text)
print("Applied successfully.")
