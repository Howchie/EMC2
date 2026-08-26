import re

with open('src/model_BAwD.cpp', 'r') as f:
    text = f.read()

# We can replace bawd_log_frozen_logn_quad
p1 = r'double bawd_log_frozen_logn_quad\(const BawdGeom& g, double s_lo,.*?return bawd_log_gl_split\(lf, s_lo, s_hi, mid, BAWD_GL_NODES\) \+\s*std::log1p\(-g\.gamma\) \+ std::log\(g\.ell\) - std::log\(g\.k\);\n\}'

r1 = """double bawd_log_frozen_logn_quad(const BawdGeom& g, double s_lo,
                                        double s_hi, double mu, double sigma,
                                        double delta) {
  split_lognormal_shape h;
  if (delta != 0.0 && !split_lognormal_shape_params(mu, sigma, delta, h))
    return R_NegInf;
  if (!g.rho_inf) {
    if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
    const double log_ell = std::log(g.ell);
    const auto lf = [&](double x) -> double {
      return bawd_log_psi_prime(g, x) +
        (delta == 0.0
          ? pnorm_log_direct((log_ell + bawd_log_wrel(g, x) - mu) / sigma, false)
          : log_split_lognormal_survivor(
              g.ell * std::exp(bawd_log_wrel(g, x)), h));
    };
    const double mid = bawd_pk_peak_x(g.rho, g.gamma, mu - log_ell);
    return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
      std::log(g.ell) - std::log(g.k);
  }
  if (!(s_hi > s_lo) || !(sigma > 0.0)) return R_NegInf;
  const double x_ell = (std::log(g.ell) - mu) / sigma;
  const auto lf = [&](double s) -> double {
    const double e1 = std::expm1(s);
    if (!(e1 > 0.0)) return R_NegInf;
    return std::log(e1) - g.gamma * s +
      (delta == 0.0
        ? pnorm_log_direct(x_ell + g.omg * s / sigma, false)
        : log_split_lognormal_survivor(
            g.ell * std::exp(g.omg * s), h));
  };
  const double mid = -sigma * x_ell / g.omg;
  return bawd_log_gl_split(lf, s_lo, s_hi, mid, BAWD_GL_NODES) +
    std::log1p(-g.gamma) + std::log(g.ell) - std::log(g.k);
}"""

# Do the replacement
new_text, n = re.subn(p1, r1, text, flags=re.DOTALL)
if n == 0:
    print("bawd_log_frozen_logn_quad replace failed")

with open('src/model_BAwD.cpp', 'w') as f:
    f.write(new_text)

