"""Independent prime checks, xi sign intervals and a polynomial bound for the zero-density integrals."""
from bisect import bisect_right
from decimal import Decimal as D, ROUND_CEILING, ROUND_FLOOR, getcontext
from functools import lru_cache
from itertools import accumulate
from math import comb, factorial, isqrt
from pathlib import Path
from flint import acb, arb, ctx, fmpq

ctx.prec = 112
getcontext().prec = 55
SCALE = 10**23

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


def small():
    primes=all_primes(225000)
    direct=PrimePrefixes(primes,1)
    slope=scaled_log(D("1.5"))[0]
    centres=[6*p for p in primes if p<150000]+[900000]
    for u in centres:
        # 12*SCALE times b(b-1) y log b at y=u/9 is u*SCALE*log b.
        require(u*slope>direct.sum_bound(u,3,2,True),"direct main term")
    low,high=9*10**5,9*10**10
    squares=PrimePrefixes([p for p in primes if p<=150000],2)
    knots=sorted({low,high}|{d*z for z in squares.n for d in (4,6,9) if low<d*z<high})
    previous=squares.sum_bound(low,3,2,False)
    for left,right in zip(knots,knots[1:]):
        next_value=squares.sum_bound(right,3,2,False)
        m=isqrt(left)
        for u,q in ((left,previous),(right,next_value)):
            require(600*m*q>11*12*SCALE*(u+m*m),"integer square tangent")
        previous=next_value
    print(f"Trial-divided primes, {len(centres)} centres and {len(knots)-1} independent square tangents.")


def zero_data():
    centers=[int(v) for v in Path(__file__).with_name("zero-centres.txt").read_text().splitlines()]
    b=arb(fmpq(3,2))
    r=b.sqrt()
    prev=arb(0)
    coeffs=[]
    invsum=arb(0)
    def xi(t):
        z=acb(arb("0.5"),t)
        return (z*(z-1)/2 * (-z*arb.pi().log()/2).exp()*(z/2).gamma()*z.zeta()).real
    for trial in centers:
        t=arb(fmpq(trial,10**9))
        left,right=t-arb("0.0000001"),t+arb("0.0000001")
        a,c=xi(left),xi(right)
        require(prev<left and ((a<0 and c>0) or (a>0 and c<0)), "distinct xi zero")
        t=left.union(right)
        coeffs.append((t,2*r*(b+1-2*r*(t*b.log()).cos())/(arb("0.25")+t*t)))
        invsum+=1/t**2
        prev=right
    require(prev<5000 and arb(5000).zeta_nzeros()==len(centers)==4520,"all first 4520 zeros")
    # Brent, Platt and Trudgian, Corollary 1, rounded upwards.
    coefficient_mass=sum(a for t,a in coeffs)+2*r*(r+1)**2*(arb("0.023104993115419")-invsum)
    require(coefficient_mass<arb("0.127"),"mass from the unrestricted inverse-square sum")
    return coeffs[:160],coefficient_mass


def integrated_tail(exponent,log_power,degree,h,H):
    # Put t=H exp(v). Since log_power <= degree, replace
    # (1+v/h)^log_power by the integer power and integrate term by term.
    a=2-exponent
    factor=sum(arb(comb(degree,j)*factorial(j))/(h**j * a**(j+1)) for j in range(degree+1))
    return 2*H**(exponent-2)*h**log_power*factor


def off_line_terms():
    H=arb(3*10**12)
    h=H.log()
    # Round each KLN Table 1 constant upward. A common 4.1 bounds every B here.
    entries=[(80,"7"),(85,"9"),(86,"9.5"),(87,"10"),(88,"10.5"),(89,"11"),
             (90,"11.5"),(91,"12.1"),(92,"12.7"),(93,"13.3"),(94,"13.9"),
             (95,"14.6"),(96,"15.2"),(97,"15.9"),(98,"16.6"),(99,"17.3")]
    all_high=arb("0.00000000000149")
    bound=(1+(H/(2*arb.pi())).log())/(2*arb.pi()*H)+(arb("0.14")+arb("0.56")*h)/H**2
    require(bound<all_high,"inverse-square high tail")
    cuts=[]
    for percent,A in entries:
        sigma=arb(fmpq(percent,100))
        p=arb(fmpq(8,3))*(1-sigma)
        z=arb(A)*integrated_tail(p,5-2*sigma,4,h,H)
        z+=arb("4.1")*integrated_tail(arb(0),arb(2),2,h,H)
        cuts.append((sigma,z.min(all_high/2)))
    return all_high,cuts


def proof():
    modes,mass=zero_data()
    all_high,cuts=off_line_terms()
    r=arb("1.5").sqrt()
    def inverse_cost(t):
        rates=[((s-arb("0.5"))*t).exp() for s,q in cuts] + [(t/2).exp()]
        total=all_high+(rates[0]-1)*all_high/2
        return 15*(total+sum((v-u)*q for u,v,(_,q) in zip(rates,rates[1:],cuts)))
    check_t=28*arb(10).log()
    require(inverse_cost(check_t)<arb("0.0511"),"unknown zero bound")
    for x,floor in [(arb(10)**5,arb("0.10")), (arb(10)**10,arb("0.123"))]:
        require(2*r*(r-1)**2-r*x.log()**2/(8*arb.pi()*x.sqrt())>floor,"square bound")
    # A coarse endpoint upper bound proves the height condition over the full theta interval.
    require(5*(arb("1.5")*10**14).sqrt()<3*10**12,"partial-RH input")
    n=2500
    start=int((5*arb(10).log()*n).floor().unique_fmpz())
    end=int((check_t*n).ceil().unique_fmpz())
    logb=arb("1.5").log()
    limit20=20*arb(10).log()
    last_block=-1
    for j in range(start,end):
        block=j//50
        if block!=last_block:
            upper_time=arb(fmpq(50*(block+1),n))
            allowance=inverse_cost(upper_time)
            last_block=block
        lower_time=arb(fmpq(j,n))
        sq_floor=arb("0.123") if lower_time>limit20 else arb("0.10")
        target=mass+allowance-sq_floor
        time=arb(fmpq(2*j+1,2*n))+logb
        acc=arb(0)
        for tau,w in modes:
            # The change in this single term is at most w*tau/(2*n).
            acc+=(2*w*(tau*time/2).cos()**2-w*tau/(2*n)).max(arb(0))
            if acc>target:
                break
        else:
            raise ArithmeticError(("midpoint bound",j))
    print("4520 sign changes, the full inverse-square sum and polynomial density integrals passed.")
    print(f"{end-start} midpoint derivative enclosures confirm positivity through 10^28.")


if __name__=="__main__":
    small()
    proof()
