"""No-quadrature analytic survivor series: feasibility, not production code.

Run after/alongside rdmswtn_gamma_feasibility.py. Adaptive integration is used
only for an independent reference. The candidate itself is a recurrence.
"""
import json
import math
from pathlib import Path

import mpmath as mp
from scipy.special import gammaln
from rdmswtn_gamma_feasibility import mixture


def series(t, d, s, mean, k, digits=None):
    # Normalize the entire factor by J_k(0). Accumulate q_n = h_n t^-n/2.
    if digits:
        num, log, exp, gamma = mp.mpf, mp.log, mp.exp, mp.gamma
    else:
        num, log, exp = float, math.log, math.exp
    t,d,s,mean,k = map(num,(t,d,s,mean,k))
    rate=k/mean
    A=d*d/(2*s*s*t)
    B=(rate*s-d/s)/t**num('.5')
    ratio = (2**num('.5')*gamma((k+1)/2)/gamma(k/2) if digits else
             math.sqrt(2)*math.exp(float(gammaln((k+1)/2)-gammaln(k/2))))
    q=[num(1), -B*ratio]
    acc=q[0]/(k+1)+q[1]/(k+2)
    terms=[q[0]/(k+1),q[1]/(k+2)]
    quiet=0
    status='term_limit'
    for n in range(2000):
        prev=q[n-2] if n>=2 else num(0)
        nxt=-(((4*A-B*B)*n+2*A-k*B*B)*q[n]+(4*A*A-2*A*B*B)*prev)/((n+2)*(n+1))
        if not digits and not math.isfinite(nxt):
            return {'status':'overflow','terms':n+2}
        q.append(nxt)
        term=nxt/(k+n+3)
        acc+=term
        terms.append(term)
        tol=num('1e-60') if digits else num('1e-17')
        quiet=quiet+1 if abs(term)<tol*abs(acc) else 0
        if quiet>=20:
            status='converged_terms'
            break
    if not digits:
        acc=math.fsum(terms)
    if acc<=0:
        return {'status':'nonpositive_sum','terms':len(q)}
    if digits:
        lc=(log(2*d)+(k-1)*log(s)+k*log(rate)-log(2*mp.pi)/2
            -log(gamma(k))+(k/2-1)*log(2)+log(gamma(k/2)))
    else:
        lc=(log(2*d)+(k-1)*log(s)+k*log(rate)-log(2*math.pi)/2
            -float(gammaln(k))+(k/2-1)*log(2)+float(gammaln(k/2)))
    ls=lc-(k+1)/2*log(t)+log(acc)
    condition=sum(abs(x) for x in terms)/abs(acc)
    return {'status':status,'terms':len(q),'logS':float(ls),
            'log10_cancellation':float(log(condition)/log(10))}


def run():
    rows=[]
    for k in (.25,1.,2.5,10.,100.):
        for t in (.02,.2,1.,5.,100.):
            ref=mixture(t,1.,1.,2.,k,'S')
            row={'k':k,'t':t,'logS_reference':ref,'double':series(t,1.,1.,2.,k)}
            if 'logS' in row['double']:
                delta=row['double']['logS']-ref
                row['double']['relative_error']=math.expm1(delta) if delta<700 else 'overflow'
            # Avoid pretending that hundreds/thousands of lost digits are practical.
            if (k/2-1)**2/t<300:
                with mp.workdps(100):
                    row['100_digit']=series(t,1.,1.,2.,k,digits=100)
                if 'logS' in row['100_digit']:
                    row['100_digit']['relative_error']=math.expm1(row['100_digit']['logS']-ref)
            rows.append(row)
    # Explicit exponential-drift CDF (k=1). Its removable delta=0
    # singularity needs a separate branch in any production implementation.
    special=[]
    with mp.workdps(60):
        Phi=lambda x:mp.erfc(-x/mp.sqrt(2))/2
        for t0 in (.02,.2,1.,5.):
            t=mp.mpf(t0); d=s=mp.mpf(1); rate=mp.mpf('.5')
            c=s*mp.sqrt(t); h=d/c; delta=rate-2*d/s**2
            J=rate*mp.exp(-h*h/2)*mp.exp(((rate-d/s**2)*s/mp.sqrt(t))**2/4)*mp.pcfd(-1,(rate-d/s**2)*s/mp.sqrt(t))/(t*mp.sqrt(2*mp.pi))
            T=Phi(-h)+t*J/rate
            H=(rate*Phi(-h)-t*J)/delta
            lf=float(mp.log(T+H))
            special.append({'t':t0,'relative_error':math.expm1(lf-mixture(t0,1.,1.,2.,1.,'F'))})
    result={'series':rows,'exponential_drift_closed_cdf':special}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:
        print(row)
    print('Exponential-drift closed CDF:',special)


if __name__=='__main__':
    run()
