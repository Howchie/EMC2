"""Positive analytic gamma/Wald tails; research, not production dispatch.

No candidate below uses quadrature. ``mixture`` is an independent reference.
The density uses scaled fractional-order cylinder functions and backward
continued fractions rather than a raw large-order cylinder function.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.special import gammainc, gammaincc, gammaln, logsumexp, pbdv

from rdmswtn_gamma_feasibility import mixture


def inverse_gamma_mixture(t, d, s, rate, k, tol=2e-14):
    """Exact positive NB mixture when q <= 1; retained NB mass bounds error.

    This prototype bounds absolute tail error, not relative error of tiny F/S.
    """
    q = rate*s*s/d
    assert 0 < q <= 1
    if q == 1:
        return dict(F=float(gammaincc((k+1)/2, d*d/(2*s*s*t))),
                    S=float(gammainc((k+1)/2, d*d/(2*s*s*t))),
                    terms=1, omitted_mass_bound=0.)
    nmax = 100000
    lf, ls = [], []
    logw = k*math.log(q)
    x = d*d/(2*s*s*t)
    for n in range(nmax):
        alpha = (k+n+1)/2
        F, S = float(gammaincc(alpha, x)), float(gammainc(alpha, x))
        lf.append(logw+math.log(F) if F else -math.inf)
        ls.append(logw+math.log(S) if S else -math.inf)
        # For k<1 ratios increase towards 1-q; for k>1 they decrease.
        ratio = (1-q)*(k+n)/(n+1)
        upper_ratio = max(ratio,1-q)
        if n > k*(1-q)/q and upper_ratio < 1:
            bound = math.exp(logw)*ratio/(1-upper_ratio)
            if bound < tol:
                break
        logw += math.log(ratio)
    else:
        raise ArithmeticError('NB term budget')
    return dict(F=math.exp(float(logsumexp(lf))), S=math.exp(float(logsumexp(ls))),
                terms=n+1, omitted_mass_bound=bound)


def fractional_logmoment(rho, z):
    """log(J_rho(z)/Gamma(rho)), 0<rho<=1, z>=0; no raw large exponent.

    For z>20 the alternating expansion is stopped with its next-term bound.
    The finite sum is an approximation with a controlled analytic remainder.
    """
    if z <= 2:
        D = float(pbdv(-rho, z)[0])
        if D <= 0 or not math.isfinite(D):
            raise ArithmeticError('fractional cylinder seed')
        return z*z/4 + math.log(D)
    if z <= 20:
        # Some raw pbdv/hyperu implementations lose 1e-10 near z=6.
        # J_rho(z)=sum (zstar-z)^n J_(rho+n)(zstar)/n!, all positive.
        zstar=32.
        step=zstar-z
        upper_ratio=step/zstar
        count=math.ceil(math.log(1e-18*(1-upper_ratio))/math.log(upper_ratio))+10
        ratios=np.empty(count)
        end=rho+count+128
        R=2*end/(math.hypot(zstar,2*math.sqrt(end))+zstar)
        for j in range(count+127,-1,-1):
            R=(rho+j)/(zstar+R)
            if j<count:
                ratios[j]=R
        term=acc=1.
        for n in range(count):
            term*=step*ratios[n]/(n+1)
            acc+=term
        return fractional_logmoment(rho,zstar)+math.log(acc)
    term, acc = 1., 1.
    for n in range(100):
        nxt = -term*(rho+2*n)*(rho+2*n+1)/(2*(n+1)*z*z)
        if abs(nxt) < 2e-16*abs(acc):
            return -rho*math.log(z)+math.log(acc)
        if abs(nxt) >= abs(term):
            raise ArithmeticError('fractional asymptotic not accurate')
        acc += nxt
        term = nxt
    raise ArithmeticError('fractional seed budget')


def log_j_sequence(t, d, s, rate, k, count):
    """log E_{Gamma(k+j,rate)}[phi((d-Vt)/c)/c], j=0..count.

    J_{nu+1}/J_nu = nu/(z+J_{nu+2}/J_{nu+1}). The continued fraction
    is positive for z>=0. Depth doubling is an empirical convergence check,
    not a rigorous rounding/remainder certificate.
    """
    p = rate*s/math.sqrt(t)
    z = p-d/(s*math.sqrt(t))
    assert z >= -1e-12
    z = max(0., z)
    offset = max(0, math.ceil(k)-1)
    rho = k-offset
    top = offset+count
    nu = rho+np.arange(top+1)
    if z*math.sqrt(float(nu[-1])) <= .5:
        lm = (nu/2-1)*math.log(2)+gammaln(nu/2)-gammaln(nu)
        if z:
            term=np.ones_like(nu)
            acc=np.ones_like(nu)
            for n in range(40):
                term *= -z*math.sqrt(2)*np.exp(gammaln((nu+n+1)/2)-gammaln((nu+n)/2))/(n+1)
                acc += term
                if np.max(np.abs(term/acc)) < 2e-16:
                    break
            else:
                raise ArithmeticError('small-z series budget')
            lm += np.log(acc)
        depth = 0
    else:
        old = None
        depth = 128
        while depth <= 262144:
            end = rho+top+depth
            R = 2*end/(math.hypot(z, 2*math.sqrt(end))+z)
            ratios = np.empty(top)
            for j in range(top+depth-1, -1, -1):
                R = (rho+j)/(z+R)
                if j < top:
                    ratios[j] = R
            if old is not None and np.max(np.abs(ratios/old-1)) < 2e-15:
                break
            old = ratios.copy()
            depth *= 2
        else:
            raise ArithmeticError('continued fraction budget')
        # Recur on the already-normalised j values. Separately accumulating
        # nu*log(p) and log(J_nu/Gamma(nu)) loses digits at large orders.
        base = (-math.log(s*math.sqrt(2*math.pi*t))-d*d/(2*s*s*t)
                +rho*math.log(p)+fractional_logmoment(rho,z))
        logs=np.empty(top+1,dtype=np.longdouble)
        logs[0]=base
        increments=np.log(np.longdouble(p)*ratios.astype(np.longdouble)/nu[:-1].astype(np.longdouble))
        logs[1:]=np.longdouble(base)+np.cumsum(increments)
        return logs[offset:],depth
    logs = (-math.log(s*math.sqrt(2*math.pi*t))-d*d/(2*s*s*t)
            +nu*math.log(p)+lm)
    return logs[offset:], depth


def positive_survivor(t, d, s, rate, k, extra=0):
    """S=(t/rate) sum_{m>=2}(1-r**(m-1))*j_{k+m}, q>=1.

    The log-concavity of the j sequence gives a geometric truncation bound.
    CF/rounding errors still need separate certification in production.
    """
    q = rate*s*s/d
    assert q >= 1
    p = rate*s/math.sqrt(t)
    peak = rate*d/t
    count = max(80, int(max(0, peak-k)+14*math.sqrt(p*p+peak)+80))+extra
    logs, depth = log_j_sequence(t,d,s,rate,k,count+2)
    r = 1-2/q
    terms = []
    for m in range(2,count+1):
        if r >= 0:
            weight = -math.expm1((m-1)*math.log(r)) if r else 1.
        else:
            weight = 1-r**(m-1)
        if weight:
            terms.append(math.log(weight)+logs[m])
    logS=math.log(t/rate)+float(logsumexp(terms))
    ratio=math.exp(float(logs[count+2]-logs[count+1]))
    if ratio>=1:
        raise ArithmeticError('not past the j mode')
    bound=(math.log(t/rate)+math.log(1 if r>=0 else 2)
           +float(logs[count+1])-math.log1p(-ratio))
    return dict(logS=logS,terms=count-1,cf_depth=depth,
                logj0=float(logs[0]),relative_truncation_bound=math.exp(bound-logS))


def early_cdf(t,d,s,rate,k,tol=2e-14):
    """Finite analytic expansion with an exact-arithmetic next-term bound.

    Not a convergent infinite series. Arbitrary precision here evaluates
    negative-order incomplete gamma accurately; no integration is used.
    Return a failure status when the bound cannot certify this regime.
    """
    delta = rate-d/(s*s)
    if delta <= 0:
        return {'status':'outside_domain'}
    with mp.workdps(60):
        t,d,s,rate,k,delta = map(mp.mpf,(t,d,s,rate,k,delta))
        a=d*d/(2*s*s)
        scale=d*rate**k/(s*mp.sqrt(2*mp.pi)*delta**k)
        acc=mp.mpf(0)
        prev=mp.inf
        coefficient=mp.mpf(1)
        for n in range(200):
            term=scale*coefficient*a**(n-mp.mpf('.5'))*mp.gammainc(mp.mpf('.5')-n,a/t,mp.inf)
            # acc is the partial sum ending at n-1; term bounds its remainder.
            if acc>term and term<mp.mpf(tol)*(acc-term):
                return dict(status='bounded',logF=float(mp.log(acc)),terms=n,
                            relative_remainder_bound=float(term/(acc-term)))
            if n>5 and term>prev:
                return {'status':'bound_not_reached','terms':n}
            acc+=(-1)**n*term
            prev=term
            coefficient*= (k+2*n)*(k+2*n+1)/((n+1)*2*s*s*delta*delta)
        return {'status':'term_budget','terms':200}


def run():
    rows=[]
    for k in (.25,1.,2.5,10.,100.):
        for t in (.02,.2,1.,5.,100.):
            d=s=1.
            rate=k/2
            q=rate
            row=dict(k=k,t=t,q=q)
            refs={side:mixture(t,d,s,2.,k,side) for side in ('F','S','f')}
            if q<=1:
                ans=inverse_gamma_mixture(t,d,s,rate,k)
                row['NB']=ans
                row['F_absolute_error']=abs(ans['F']-math.exp(refs['F']))
                row['S_absolute_error']=abs(ans['S']-math.exp(refs['S']))
            if q>=1:
                ans=positive_survivor(t,d,s,rate,k)
                row['positive']=ans
                row['S_relative_error']=math.expm1(ans['logS']-refs['S'])
                row['pdf_relative_error']=math.expm1(math.log(d/t)+ans['logj0']-refs['f'])
                again=positive_survivor(t,d,s,rate,k,extra=200)
                row['extra_terms_logS_difference']=again['logS']-ans['logS']
                early=early_cdf(t,d,s,rate,k)
                if early['status']=='bounded':
                    early['relative_error']=math.expm1(early['logF']-refs['F'])
                row['early_CDF']=early
            rows.append(row)
    # q=1 is an exact inverse-gamma law for every real k; vary d and s too.
    junction=[]
    for k in (.3,2.5,100.):
        for d,s,t in ((.2,.7,.04),(2.,.5,3.)):
            rate=d/(s*s)
            exact=inverse_gamma_mixture(t,d,s,rate,k)
            ref=mixture(t,d,s,k/rate,k,'F')
            junction.append(dict(k=k,d=d,s=s,t=t,
                absolute_error=abs(exact['F']-math.exp(ref))))
    stress=[]
    for k in (.3,3.7,100.):
        for q in (1.0001,1.5,2.,20.):
            for A in (.1,1.,20.):
                d,s=1.3,.7
                rate=q*d/(s*s)
                t=d*d/(2*s*s*A)
                ans=positive_survivor(t,d,s,rate,k)
                ref=mixture(t,d,s,k/rate,k,'S')
                stress.append(dict(k=k,q=q,A=A,
                    relative_error=math.expm1(ans['logS']-ref),
                    terms=ans['terms'],relative_truncation_bound=ans['relative_truncation_bound']))
    result=dict(scope='Research; exact positive series, no production speed claim',
                cases=rows,inverse_gamma_identity=junction,stress=stress)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print('cases:',len(rows))
    print('NB maximum absolute error:',max(max(x.get('F_absolute_error',0),x.get('S_absolute_error',0)) for x in rows))
    print('positive maximum relative survivor error:',max(abs(x.get('S_relative_error',0)) for x in rows))
    print('scaled PDF maximum relative error:',max(abs(x.get('pdf_relative_error',0)) for x in rows))
    print('max terms:',max(x.get('positive',{}).get('terms',0) for x in rows))
    print('early CDF certified:',sum(x.get('early_CDF',{}).get('status')=='bounded' for x in rows))
    print('stress cases / maximum relative error:',len(stress),max(abs(x['relative_error']) for x in stress))
    for row in rows:
        if row['k']==100:
            print(row)
    assert max(abs(x.get('S_relative_error',0)) for x in rows)<2e-9
    assert max(abs(x.get('pdf_relative_error',0)) for x in rows)<2e-9
    assert max(max(x.get('F_absolute_error',0),x.get('S_absolute_error',0)) for x in rows)<1e-11
    assert max(x['absolute_error'] for x in junction)<1e-10
    assert max(abs(x['relative_error']) for x in stress)<2e-9
    assert max(x['relative_truncation_bound'] for x in stress)<1e-14
    for row in rows:
        early=row.get('early_CDF',{})
        if early.get('status')=='bounded':
            assert abs(early['relative_error'])<early['relative_remainder_bound']+1e-12


if __name__=='__main__':
    run()
