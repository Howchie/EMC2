import sympy as sp

D, v, s, k, b, a, lam = sp.symbols('D v s k b a lam', real=True, positive=True)

# Normal PDF
pdf = sp.exp(-(D - v)**2 / (2 * s**2)) / (sp.sqrt(2 * sp.pi) * s)

# Laplace transform term for BAwL
# e^{-lam T} where T = -1/k * ln((D - k*b)/(D - k*a))
exp_lam_T = ((D - k*b)/(D - k*a))**(lam/k)

# Integral over D
integral = sp.Integral(pdf * exp_lam_T, (D, k*b, sp.oo))

print("Is integrable?", integral.doit() != integral)
