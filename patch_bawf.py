import re

with open('src/model_BAwF.cpp', 'r') as f:
    text = f.read()

# Replace bawf_log_frozen_quad
p1 = r'double bawf_log_frozen_quad\(const BawfGeom& g, double s_lo, double s_hi,.*?std::log\(g\.b\);\n\}'
r1 = """double bawf_log_frozen_quad(const BawfGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift, double delta) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (logn && delta != 0.0 && !split_lognormal_shape_params(p1, p2, delta, h))
    return R_NegInf;
  const double log_kb = std::log(g.k) + std::log(g.b);
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_jac = std::log(s) + bawf_log_Hpp(g, s);
    const double log_surv = logn
      ? (delta == 0.0
           ? pnorm_log_direct((log_kb + bawf_log_Hp(g, s) - p1) / p2, false)
           : log_split_lognormal_survivor(
               std::exp(log_kb + bawf_log_Hp(g, s)), h))
      : pnorm_log_direct((p1 - bawf_critical_launch(g, s)) / p2, true);
    return log_jac + log_surv;
  };
  (void)posdrift;
  const double mid = logn
    ? bawf_s_of_logratio(g, p1 - log_kb)
    : ((p1 > 0.0) ? bawf_s_of_logratio(g, std::log(p1) - log_kb) : s_lo);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log(g.b);
}"""
text = re.sub(p1, r1, text, flags=re.DOTALL)

# Replace bawf_log_frozen_surv_quad
p2 = r'double bawf_log_frozen_surv_quad\(const BawfGeom& g, double s_lo,.*?bawf_log_gl_split\(lf, s_lo, s_hi, mid, BAWD_GL_NODES\);\n\}'
r2 = """double bawf_log_frozen_surv_quad(const BawfGeom& g, double s_lo,
                                        double s_hi, double p1, double p2,
                                        bool logn, bool posdrift,
                                        double delta) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (logn && delta != 0.0 && !split_lognormal_shape_params(p1, p2, delta, h))
    return R_NegInf;
  const double log_kb = std::log(g.k) + std::log(g.b);
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_jac = std::log(s) + bawf_log_Hpp(g, s) + std::log(g.b);
    const double log_cdf = logn
      ? (delta == 0.0
           ? pnorm_log_direct((log_kb + bawf_log_Hp(g, s) - p1) / p2, true)
           : log_split_lognormal_cdf(
               std::exp(log_kb + bawf_log_Hp(g, s)), h))
      : (posdrift
           ? log_normal_cdf_positive_raw(bawf_critical_launch(g, s), p1, p2)
           : pnorm_log_direct((bawf_critical_launch(g, s) - p1) / p2, true));
    return log_jac + log_cdf;
  };
  const double mid = logn
    ? bawf_s_of_logratio(g, p1 - log_kb)
    : ((p1 > 0.0) ? bawf_s_of_logratio(g, std::log(p1) - log_kb) : s_lo);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES);
}"""
text = re.sub(p2, r2, text, flags=re.DOTALL)

with open('src/model_BAwF.cpp', 'w') as f:
    f.write(text)
print("BAwF patch applied")
