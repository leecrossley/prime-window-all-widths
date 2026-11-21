"""Separate enclosure of the finite and phase ranges at 151/100.

The trial centres below are from Andrew Odlyzko's zeta table
https://www-users.cse.umn.edu/~odlyzko/zeta_tables/ . Opposite signs of
the completed zeta function certify a zero in each interval, independently
of the nominal precision and indexing of the table.
"""
from bisect import bisect_right
from decimal import Decimal as D, ROUND_CEILING, ROUND_FLOOR, getcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import accumulate
from math import isqrt
from flint import acb, arb, ctx, fmpq

getcontext().prec=55
SCALE=10**23
CENTRES = (
    '14.134725142', '21.022039639', '25.010857580', '30.424876126',
    '32.935061588', '37.586178159', '40.918719012', '43.327073281',
    '48.005150881', '49.773832478', '52.970321478', '56.446247697',
    '59.347044003', '60.831778525', '65.112544048', '67.079810529',
    '69.546401711', '72.067157674', '75.704690699', '77.144840069',
    '79.337375020', '82.910380854', '84.735492981', '87.425274613',
)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def all_primes(n):
    """Trial division by previously established primes, independent of the release sieve."""
    out = [2]
    for p in range(3, n + 1, 2):
        for q in out:
            if q*q > p:
                out.append(p)
                break
            if p % q == 0:
                break
    return out


@lru_cache(maxsize=None)
def scaled_log(x):
    # Decimal.ln is correctly rounded. The two neighboring decimals enclose the exact log.
    center = D(x).ln()
    low = int((center.next_minus()*SCALE).to_integral_value(rounding=ROUND_FLOOR))
    high = int((center.next_plus()*SCALE).to_integral_value(rounding=ROUND_CEILING))
    return low, high


class PrimePrefixes:
    """Integer lower and upper prefix sums of log(p) and log(p)/n for n=p or n=p²."""
    def __init__(self, primes, exponent):
        self.n = [p**exponent for p in primes]
        lower, upper = zip(*(scaled_log(p) for p in primes))
        self.log_lo = [0, *accumulate(lower)]
        self.log_hi = [0, *accumulate(upper)]
        self.ratio_lo = [0, *accumulate(v//n for v, n in zip(lower, self.n))]
        self.ratio_hi = [0, *accumulate(-(-v//n) for v, n in zip(upper, self.n))]

    def sum_bound(self, u, a, c, upper):
        """Return a bound on a*c²*SCALE times the weighted sum at y=u/a²."""
        below, middle, end = (bisect_right(self.n, u//d) for d in (a*a, a*c, c*c))
        logs_pos, logs_neg, ratios_pos, ratios_neg = (
            (self.log_hi, self.log_lo, self.ratio_hi, self.ratio_lo) if upper
            else (self.log_lo, self.log_hi, self.ratio_lo, self.ratio_hi))
        # Direct evaluation of the two weight branches using signed interval arithmetic.
        left = a*a*c*(logs_pos[middle]-logs_pos[below]) - c*u*(ratios_neg[middle]-ratios_neg[below])
        right = a*u*(ratios_pos[end]-ratios_pos[middle]) - a*c*c*(logs_neg[end]-logs_neg[middle])
        return left + right


ctx.prec=144
a,c=151,100
A=a*a
b=arb(fmpq(a,c))
r=b.sqrt()


def arithmetic():
    ps=all_primes(228_010)
    direct=PrimePrefixes(ps,1)
    corners=sorted({A*100_000}|{u for p in ps for u in (A*p,a*c*p,c*c*p) if u<=A*100_000})
    scale=a*c*c*SCALE
    slope=(a-c)*scaled_log(D("1.51"))[0]
    g,u=min((slope*x-direct.sum_bound(x,a,c,upper=True),x) for x in corners)
    require(1000*g>59*scale,"direct finite lower bound")
    print(f"{len(corners):,} finite positive breakpoints; least certified grid deficit {D(g)/D(scale):.12f} at y={F(u,A)}.")

    sq=PrimePrefixes([p for p in ps if p<=151_000],2)
    lo,hi=A*100_000,A*10**10
    corners=sorted({lo,hi}|{x for p in sq.n for x in (p*A,p*a*c,p*c*c) if lo<x<hi})
    bounds=[sq.sum_bound(x,a,c,upper=False) for x in corners]
    for left,right,ql,qr in zip(corners,corners[1:],bounds,bounds[1:]):
        # The tangent taken at the square just below the LEFT end has rational coefficients.
        d=isqrt(left)
        for x,q in ((left,ql),(right,qr)):
            require(1000*2*a*d*q>115*scale*(x+d*d),"square tangent")
    print(f"{len(corners)-1:,} prime-square intervals with Q(y) > 0.115 sqrt(y).")


def xi_real(t):
    s=acb(arb(fmpq(1,2)),t)
    return (s*(s-1)/2 * (-s*arb.pi().log()/2).exp() * (s/2).gamma() * s.zeta()).real


def phases():
    N=24
    pairs=[]
    prev=arb(0)
    for row in CENTRES:
        centre=arb(row.strip())
        delta=arb("0.0000001")
        left,right=centre-delta,centre+delta
        vl,vr=xi_real(left),xi_real(right)
        require(left>prev and ((vl<0 and vr>0) or (vl>0 and vr<0)),"xi sign change")
        tau=left.union(right)
        u=tau*b.log()/2
        coeff=2*r*((r-1)**2+4*r*u.sin()**2)/(arb(fmpq(1,4))+tau*tau)
        pairs.append((tau,coeff))
        prev=right
    require(len(pairs)==N,"incomplete trial ordinates")
    K=sum(tau*w for tau,w in pairs)
    k=1000
    best=arb(100)
    for j in range(11*k,35*k):
        t=arb(fmpq(2*j+1,2*k))
        f=sum(2*w*(t*tau/2).cos()**2 for tau,w in pairs)
        lower=f-K/(2*k)
        require(lower>arb("0.020"),"phase gap")
        if lower<best:
            best=lower
    require(K < arb("3.2") and best > arb("0.0248"), "independent bounds")
    print("24 disjoint xi-sign intervals. Lipschitz constant below 3.2 and every mesh lower bound above 0.0248.")
    print("The phase gain is > 0.020 over the full closed interval [11,35].")


if __name__=="__main__":
    arithmetic()
    phases()
