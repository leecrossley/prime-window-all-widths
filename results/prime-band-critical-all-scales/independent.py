"""Independent prime, uniform phase and cube checks for the critical dilation interval."""
from pathlib import Path
from bisect import bisect_right
from decimal import Decimal as D, ROUND_CEILING, ROUND_FLOOR, getcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import accumulate
from math import isqrt
from flint import acb, arb, ctx, fmpq

getcontext().prec=55
SCALE=10**23
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


def other_finite():
    ps=all_primes(228010)
    pref=PrimePrefixes(ps,1)
    c=10000;count=0;worst=None
    for l in range(15000,15100):
        a=l+1
        knots=[p*a*c for p in ps[:bisect_right(ps,100000*a//c)]] + [100000*a*a]
        lnlo=scaled_log(D(l)/D(c))[0]
        for u in knots:
            val=l*(l-c)*u*lnlo-a*pref.sum_bound(u,a,c,True)
            if val<=0:raise ArithmeticError((l,u,val))
            if worst is None or val*worst[1]<worst[0]*a*a*c*c*SCALE:
                worst=(val,a*a*c*c*SCALE,F(l,c),F(u,a*a))
        count+=len(knots)
    require(D(worst[0])/D(worst[1])>D("0.058"),"prime floor")
    print(f"{count} endpoints in 100 dilation slices have positive gap above 0.058.")
    a,c=3,2;lo=a*a*100000;hi=a*a*10**10
    sq=PrimePrefixes([p for p in ps if p<=150000],2)
    knots=sorted({lo,hi}|{u for n in sq.n for u in (c*c*n,a*c*n,a*a*n) if lo<u<hi})
    values=[sq.sum_bound(u,a,c,False) for u in knots]
    for i,j,x,z in zip(knots,knots[1:],values,values[1:]):
        d=isqrt(i)
        for u,q in ((i,x),(j,z)):
            if 100*2*a*d*q<=11*a*c*c*SCALE*(u+d*d):raise ArithmeticError((i,j))
    print(f"{len(knots)-1} square intervals pass the independent 0.11 tangent bound.")


def other_phases():
    rows=[int(row) for row in Path(__file__).with_name("zero-centres.txt").read_text().splitlines()]
    require(len(rows)==64,"expected 64 trial centres")
    z=[]
    prev=arb(0)
    for raw in rows:
        mid=arb(fmpq(raw,10**9))
        left,right=mid-arb("0.0000001"),mid+arb("0.0000001")
        def f(t):
            s=acb(arb("0.5"),t)
            return (s*(s-1)/2*(-s*arb.pi().log()/2).exp()*(s/2).gamma()*s.zeta()).real
        a,b=f(left),f(right)
        require(left>prev and ((a<0 and b>0) or (a>0 and b<0)),"xi sign change")
        z.append(left.union(right))
        prev=right
    low,high=arb("1.50342"),arb("1.51")
    root=low.sqrt()
    midpoint=(low.log()+high.log())/2
    radius=(high.log()-low.log())/2
    weights=[]
    for tau in z:
        # A lower amplitude over the dilation interval from Lipschitz continuity of sin.
        w=(abs((tau*midpoint/2).sin())-tau*radius/2).lower().max(arb(0))
        coeff=2*root*((root-1)**2+4*root*w**2)/(arb("0.25")+tau*tau)
        weights.append(coeff.lower())
    lipschitz=sum(t*w for t,w in zip(z,weights))
    n=800
    for j in range(11*n,97*n):
        time=arb(fmpq(2*j+1,2*n))
        midgain=sum(2*w*(tau*time/2).cos()**2 for tau,w in zip(z,weights))
        require(midgain-lipschitz/(2*n)>arb("0.022"),"phase midpoint lower bound")
    b=low.union(high)
    require(arb("1.5")*(b.root(3)-1)*(b-b.root(3))>arb("0.077"),"cube drift")
    require(2*b.sqrt()*(b-1)<arb("1.26"),"square variation")
    require(arb("1.5")*(b-b.root(3)+b*b.root(3)-b.root(3)**2)<arb("1.18"),"cube variation")
    require(arb("0.077")-(arb("1.26")*48**2/arb(8).exp()+
                           arb("1.18")*32**2/arb(16).exp())/(8*arb.pi())>arb("0.035"),"infinite tail")
    print("64 certified sign crossings and 68800 Lipschitz cells give uniform gain > 0.022.")
    print("The first omitted power supplies more than 0.035 y^(1/3) after both prime-counting errors.")


if __name__=="__main__":
    other_finite()
    other_phases()
