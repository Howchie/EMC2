import re

with open('src/wald_functions.h', 'r') as f:
    text = f.read()

def replace_block(text, start_str, end_str, new_str):
    start = text.find(start_str)
    if start == -1:
        print("Not found:", start_str)
        return text
    end = text.find(end_str, start) + len(end_str)
    return text[:start] + new_str + text[end:]

# 1. split_lognormal_log_moment_tail
t_start = "inline double split_lognormal_log_moment_tail(double v, double mu,"
t_end = "return log_sum_exp(lpL, lpR);\n}"
t_new = """inline double split_lognormal_log_moment_tail(double v, const split_lognormal_shape& h, double r) {
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
text = replace_block(text, t_start, t_end, t_new)

# 2. split_lognormal_log_total_moment
tm_start = "inline double split_lognormal_log_total_moment(double mu, double sigma,"
tm_end = "return split_lognormal_log_moment_tail(0.0, mu, sigma, delta, r);\n}"
tm_new = """inline double split_lognormal_log_total_moment(const split_lognormal_shape& h, double r) {
  return split_lognormal_log_moment_tail(0.0, h, r);
}

inline double split_lognormal_log_total_moment(double mu, double sigma,
                                               double delta, double r) {
  return split_lognormal_log_moment_tail(0.0, mu, sigma, delta, r);
}"""
text = replace_block(text, tm_start, tm_end, tm_new)

# 3. log_split_lognormal_survivor
s_start = "inline double log_split_lognormal_survivor(double v, double mu, double sigma,"
s_end = "return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);\n}"
s_new = """inline double log_split_lognormal_survivor(double v, const split_lognormal_shape& h) {
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, h, 0.0);
}

inline double log_split_lognormal_survivor(double v, double mu, double sigma,
                                           double delta) {
  if (delta == 0.0) return lnorm_log_surv_std(v, mu, sigma);
  if (!(v > 0.0)) return 0.0;
  return split_lognormal_log_moment_tail(v, mu, sigma, delta, 0.0);
}"""
text = replace_block(text, s_start, s_end, s_new)

# 4. log_split_lognormal_stoploss
sl_start = "inline double log_split_lognormal_stoploss(double v, double mu, double sigma,"
sl_end = "return log_sum_exp(lpL, lpR);\n}"
sl_new = """inline double log_split_lognormal_stoploss(double v, const split_lognormal_shape& h) {
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
text = replace_block(text, sl_start, sl_end, sl_new)

# 5. log_split_lognormal_put
pt_start = "inline double log_split_lognormal_put(double v, double mu, double sigma,"
pt_end = "return log_sum_exp(lpL, lpR);\n}"
pt_new = """inline double log_split_lognormal_put(double v, const split_lognormal_shape& h) {
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
text = replace_block(text, pt_start, pt_end, pt_new)

# 6. log_split_lognormal_first_partial_moment
fm_start = "inline double log_split_lognormal_first_partial_moment(double v, double mu,"
fm_end = "return split_lognormal_log_moment_tail(v, mu, sigma, delta, 1.0);\n}"
fm_new = """inline double log_split_lognormal_first_partial_moment(double v, const split_lognormal_shape& h) {
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
text = replace_block(text, fm_start, fm_end, fm_new)

# 7. log_split_lognormal_power_stoploss
ps_start = "inline double log_split_lognormal_power_stoploss(double v, double mu,"
ps_end = "return out.sign > 0 ? out.log_abs - std::log(std::fabs(m)) : R_NegInf;\n}"
ps_new = """inline double log_split_lognormal_power_stoploss(double v, const split_lognormal_shape& h, double m) {
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
text = replace_block(text, ps_start, ps_end, ps_new)

# 8. log_split_lognormal_density
d_start = "inline double log_split_lognormal_density(double v, double mu, double sigma,"
d_end = "return std::log(2.0 * mass) - std::log(v * s * M_SQRT2PI) - 0.5 * z * z;\n}"
d_new = """inline double log_split_lognormal_density(double v, const split_lognormal_shape& h) {
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
text = replace_block(text, d_start, d_end, d_new)

# 9. log_split_lognormal_cdf
c_start = "inline double log_split_lognormal_cdf(double v, double mu, double sigma,"
c_end = "return (ls < 0.0) ? log1m_exp(ls) : R_NegInf;\n}"
c_new = """inline double log_split_lognormal_cdf(double v, const split_lognormal_shape& h) {
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
text = replace_block(text, c_start, c_end, c_new)

with open('src/wald_functions.h', 'w') as f:
    f.write(text)

