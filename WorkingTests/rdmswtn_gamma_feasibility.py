"""Research checks only; no production kernels are changed.

Run: python WorkingTests/rdmswtn_gamma_feasibility.py
Requires numpy, scipy, mpmath. Results are written beside this script.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, minimize_scalar
from scipy.special import gammaln, log_ndtr, pbdv


def conditional(t, d, s, v, kind):
    c = s * math.sqrt(t)
    if kind == "f":
        return math.log(d / (s * math.sqrt(2 * math.pi))) - 1.5 * math.log(t) - (d-v*t)**2/(2*s*s*t)
    a = float(log_ndtr((v*t-d)/c))
    b = 2*d*v/(s*s) + float(log_ndtr(-(v*t+d)/c))
    if kind == "F":
        return float(np.logaddexp(a, b))
    a = float(log_ndtr((d-v*t)/c))
    # Stable difference; the tested domain has d > 0 and moderate v*t/c.
    if b >= a:
        raise ArithmeticError("conditional survivor cancellation")
    return a + math.log(-math.expm1(b-a))


def mixture(t, d, s, mean, k, kind):
    """Independent log-rate adaptive reference, centered on each integrand.

    Dropping 55 log units is an empirical truncation, not a certified bound.
    The log-rate Jacobian removes the Gamma singularity for shape < 1.
    """
    rate = k/mean
    norm = k*math.log(rate) - float(gammaln(k))

    def logfun(u):
        v = math.exp(u)
        return norm + k*u-rate*v + conditional(t, d, s, v, kind)

    mid = math.log(mean)
    opt = minimize_scalar(lambda u: -logfun(u), bounds=(mid-max(30,80/k), mid+15), method="bounded",
                          options={"xatol": 1e-12})
    u0, top = float(opt.x), -float(opt.fun)
    width = 1/math.sqrt(k)
    ends = []
    for sign in (-1, 1):
        w = width
        while logfun(u0+sign*w) > top-55:
            w *= 2
        ends.append(brentq(lambda u: logfun(u)-(top-55), *sorted((u0, u0+sign*w))))
    # Splitting at the mode prevents missing a narrow central peak.
    val = sum(quad(lambda u: math.exp(logfun(u)-top), a, b,
                   epsabs=1e-13, epsrel=2e-12, limit=200)[0]
              for a, b in ((ends[0],u0),(u0,ends[1])))
    return top+math.log(val)


def cylinder_logpdf(t, d, s, mean, k):
    with mp.workdps(60):
        t,d,s,mean,k = map(mp.mpf,(t,d,s,mean,k))
        rate = k/mean
        z = (rate-d/s**2)*s/mp.sqrt(t)
        return float(mp.log(d/(s*mp.sqrt(2*mp.pi))) - mp.mpf('1.5')*mp.log(t)
                     + k*mp.log(rate*s/mp.sqrt(t))-d*d/(2*s*s*t)
                     + z*z/4 + mp.log(mp.pcfd(-k,z)))


def laguerre(t, d, s, mean, k, kind, n):
    # Normalized Golub-Welsch weights avoid Gamma(k) overflow at large k.
    i = np.arange(n, dtype=float)
    x, vec = eigh_tridiagonal(2*i+k, np.sqrt(np.arange(1,n)*(np.arange(1,n)+k-1)))
    w = vec[0]**2
    return sum(wi*math.exp(conditional(t,d,s,xi*mean/k,kind)) for xi,wi in zip(x,w) if wi > 0)


def run():
    out = {"scope": "point distance; fixed diffusion; no clocks; deterministic research checks"}
    errors, raw_failures, cdf_checks = [], [], []
    cases = [(t,d,s,2.0,k) for k in (.25,1.,2.5,10.,100.)
             for t,d,s in ((.02,.2,1.),(.2,1.,1.),(1.,1.,1.),(5.,2.,.5))]
    for t,d,s,m,k in cases:
        ref = mixture(t,d,s,m,k,"f")
        exact = cylinder_logpdf(t,d,s,m,k)
        errors.append(abs(math.expm1(exact-ref)))
        z = (k/m-d/s**2)*s/math.sqrt(t)
        raw = float(pbdv(-k,z)[0])
        if not math.isfinite(raw) or raw <= 0:
            raw_failures.append({"case":[t,d,s,m,k],"z":z,"D":str(raw)})
        F,S = (math.exp(mixture(t,d,s,m,k,kind)) for kind in ("F","S"))
        h = t*1e-4
        side = "F" if F < .5 else "S"
        derivative = (math.exp(mixture(t+h,d,s,m,k,side))-math.exp(mixture(t-h,d,s,m,k,side)))/(2*h)
        if side == "S":
            derivative = -derivative
        cdf_checks.append({"case":[t,d,s,m,k],"F_plus_S_error":abs(F+S-1),
                           "derivative_relative_error":abs(derivative/math.exp(ref)-1)})
    out["density_cases"] = len(cases)
    out["max_cylinder_vs_quadrature_relative_error"] = max(errors)
    out["naive_scipy_cylinder_failures"] = raw_failures
    out["max_F_plus_S_error"] = max(x["F_plus_S_error"] for x in cdf_checks)
    out["max_central_derivative_relative_error"] = max(x["derivative_relative_error"] for x in cdf_checks)
    out["cdf_checks"] = cdf_checks
    out["prior_laguerre"] = []
    for k in (.25,1.,4.,100.):
        for t in (.2,2.,100.):
            t,d,s,m = t,1.,1.,2.
            refs = {kind:math.exp(mixture(t,d,s,m,k,kind)) for kind in ("f","F","S")}
            for n in (32,64,128):
                errs = {kind:abs(laguerre(t,d,s,m,k,kind,n)/refs[kind]-1) for kind in refs}
                out["prior_laguerre"].append({"k":k,"t":t,"n":n,"relative_errors":errs})
    out["tail_checks"] = []
    for k in (.5,1.,3.):
        d,s,m = 1.,1.,2.
        rate = k/m
        logC = (math.log(d)+(k-1)*math.log(s)+k*math.log(rate)+(k/2-1)*math.log(2)
                +float(gammaln(k/2))-float(gammaln(k))-.5*math.log(2*math.pi))
        for t in (1e4,1e6):
            lf,ls = (mixture(t,d,s,m,k,kind) for kind in ("f","S"))
            out["tail_checks"].append({"k":k,"t":t,
                "f_over_asymptotic":math.exp(lf-logC+(k+3)/2*math.log(t)),
                "S_over_asymptotic":math.exp(ls-logC-math.log(2/(k+1))+(k+1)/2*math.log(t))})
    # The exact fixed-rate limit approached by growing Gamma shape.
    out["constant_drift_limit"] = []
    for k in (100.,10000.,1e6):
        errs = {kind:math.expm1(mixture(.5,1.,1.,2.,k,kind)-conditional(.5,1.,1.,2.,kind))
                for kind in ("f","F","S")}
        out["constant_drift_limit"].append({"k":k,"relative_difference":errs})
    assert max(errors) < 1e-9
    assert out["max_F_plus_S_error"] < 1e-10
    assert out["max_central_derivative_relative_error"] < 1e-5
    path = Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:val for key,val in out.items() if key not in
                     ('cdf_checks','prior_laguerre','tail_checks')},indent=2))
    print('Worst prior-Laguerre relative error by node count:')
    for n in (32,64,128):
        print(n, max(max(x['relative_errors'].values()) for x in out['prior_laguerre'] if x['n']==n))
    print('Full results:',path)


if __name__ == '__main__':
    run()
