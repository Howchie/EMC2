import math

def phi(beta, t):
    if t <= 0:
        return 0.0
    return math.exp(-beta*beta/(2*t)) / math.sqrt(2*math.pi*t)

def g_formula(t, beta):
    if t <= 0:
        return 0.0
    nu_t = -phi(beta, t)
    N = 2000
    u_max = math.sqrt(t)
    integral = 0.0
    for i in range(N):
        u_L = u_max * (1 - i / N)
        u_R = u_max * (1 - (i + 1) / N)
        u_M = 0.5 * (u_L + u_R)
        s_L = t - u_L * u_L
        s_R = t - u_R * u_R
        s_M = t - u_M * u_M
        nu_L = -phi(beta, s_L) if s_L > 0 else 0.0
        nu_R = -phi(beta, s_R) if s_R > 0 else 0.0
        nu_M = -phi(beta, s_M) if s_M > 0 else 0.0
        dt_L = t - s_L
        dt_R = t - s_R
        dt_M = t - s_M
        f_L = (nu_L - nu_t) / (math.sqrt(2 * math.pi) * dt_L * math.sqrt(dt_L)) if dt_L > 0 else 0.0
        f_R = (nu_R - nu_t) / (math.sqrt(2 * math.pi) * dt_R * math.sqrt(dt_R)) if dt_R > 0 else 0.0
        f_M = (nu_M - nu_t) / (math.sqrt(2 * math.pi) * dt_M * math.sqrt(dt_M)) if dt_M > 0 else 0.0
        panel_h_u = u_L - u_R
        integral += (panel_h_u / 6.0) * ((2.0 * u_L * f_L) + (8.0 * u_M * f_M) + (2.0 * u_R * f_R))
    term1 = 0.5 * integral
    term2 = -nu_t / math.sqrt(2 * math.pi * t)
    term3 = beta * math.exp(-beta * beta / (2 * t)) / (2 * math.sqrt(2 * math.pi) * (t ** 1.5))
    density = term1 + term2 + term3
    return density

def g_exact(t, beta):
    return beta * math.exp(-beta * beta / (2 * t)) / (math.sqrt(2 * math.pi) * (t ** 1.5))

beta = 1.0
for t in [0.1, 0.2, 0.5, 1.0]:
    approx = g_formula(t, beta)
    exact = g_exact(t, beta)
    print(t, approx, exact, approx / exact)
