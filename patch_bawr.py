import re

with open('src/model_BAwR.cpp', 'r') as f:
    text = f.read()

p1 = r'double bawr_log_frozen_quad\(const BawrGeom& g, double s_lo, double s_hi,.*?BAWD_GL_NODES\);\n\}'
r1 = """double bawr_log_frozen_quad(const BawrGeom& g, double s_lo, double s_hi,
                                   double p1, double p2, bool logn,
                                   bool posdrift, double delta) {
  if (!(s_hi > s_lo) || !(p2 > 0.0)) return R_NegInf;
  split_lognormal_shape h;
  if (logn && delta != 0.0 && !split_lognormal_shape_params(p1, p2, delta, h))
    return R_NegInf;
  const auto lf = [&](double s) -> double {
    if (!(s > 0.0)) return R_NegInf;
    const double log_vc = bawr_log_critical_launch(g, s);
    const double log_jac = g.log_kappa + std::log(g.pw) + g.pw * std::log(s);
    const double log_surv = logn
      ? (delta == 0.0
           ? pnorm_log_direct((log_vc - p1) / p2, false)
           : log_split_lognormal_survivor(std::exp(log_vc), h))
      : pnorm_log_direct((p1 - std::exp(log_vc)) / p2, true);
    return log_jac + log_surv;
  };
  (void)posdrift;
  double mid = s_lo;
  const double log_ref = logn ? p1 : ((p1 > 0.0) ? std::log(p1) : R_NegInf);
  if (log_ref > R_NegInf) mid = std::exp((log_ref - g.log_kappa) / g.pw);
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES);
}"""
text = re.sub(p1, r1, text, flags=re.DOTALL)

with open('src/model_BAwR.cpp', 'w') as f:
    f.write(text)
print("BAwR patch applied")
