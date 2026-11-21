"""Interval certificates for the ordinary-prime tent at dilation 151/100.

Python 3.10+ with python-flint 0.9.0. All comparisons are strict comparisons
of enclosing Arb balls. Integer coordinates represent exact rational knots.
"""
from bisect import bisect_right
from math import isqrt
from flint import acb, arb, ctx, fmpq

ctx.prec = 112
A, C = 151, 100
DEN = A*A
B = arb(fmpq(A,C))
ROOT = B.sqrt()
LOG = B.log()
P0, P1, M = 100_000, 10**10, 21_000_000
ZERO_PAIRS = 40
CELLS_PER_UNIT = 1000


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def primes(limit):
    flag = bytearray(b"\x01")*(limit+1)
    flag[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit)+1):
        if flag[p]:
            flag[p*p::p] = bytes((limit-p*p)//p+1)
    return [p for p in range(2,limit+1) if flag[p]]


class Band:
    def __init__(self, ps, exponent):
        self.terms = [p**exponent for p in ps]
        self.logs = [arb(0)]
        self.ratios = [arb(0)]
        for p,n in zip(ps, self.terms):
            lp = arb(p).log()
            self.logs.append(self.logs[-1]+lp)
            self.ratios.append(self.ratios[-1]+lp/n)

    def at(self, u):
        # For y=u/A² the three boundaries y, by and b²y are all integer
        # comparisons. Endpoints shared by the branches have identical weight.
        i,j,k = (bisect_right(self.terms, u//v) for v in (A*A,A*C,C*C))
        y = arb(fmpq(u,DEN))
        left = B*((self.logs[j]-self.logs[i]) - y*(self.ratios[j]-self.ratios[i]))
        right = B*B*y*(self.ratios[k]-self.ratios[j]) - (self.logs[k]-self.logs[j])
        return left+right

    def knots(self, start, stop):
        return sorted({start,stop}|{u for n in self.terms for u in (n*C*C,n*A*C,n*A*A)
                                   if start<u<stop})


def finite_prime_ranges():
    ps = primes(A*A*P0//(C*C))
    direct = Band(ps,1)
    knots = direct.knots(0,DEN*P0)
    main = B*(B-1)*LOG
    for u in knots[1:]:
        check(main*arb(fmpq(u,DEN)) - direct.at(u) > arb("0.059"),
              f"ordinary-prime knot {u}/{DEN}")
    print(f"Ordinary primes: {len(knots)-1} positive knots through {P0}, each deficit > 0.059.")

    sq_limit = A*isqrt(P1)//C
    squares = Band([p for p in ps if p<=sq_limit],2)
    knots = squares.knots(DEN*P0,DEN*P1)
    left_value = squares.at(knots[0])
    for u,v in zip(knots,knots[1:]):
        right_value = squares.at(v)
        d = isqrt((u+v)//2)
        # An affine upper tangent for sqrt(w)/A is (w+d²)/(2 A d).
        for w,value in ((u,left_value),(v,right_value)):
            check(value > arb(fmpq(11*(w+d*d),100*2*A*d)),
                  f"prime-square interval {u}/{DEN} to {v}/{DEN}")
        left_value = right_value
    print(f"Prime squares: {len(knots)-1} intervals cover [{P0}, {P1}] with Q(y) > 0.11 sqrt(y).")


def spectral_coefficients():
    zeros = acb.zeta_zeros(1,ZERO_PAIRS)
    check(len(zeros)==ZERO_PAIRS,"zeta-zero count")
    pairs = []
    previous = arb(0)
    for rho in zeros:
        tau = rho.imag
        check(rho.real == arb(fmpq(1,2)) and tau>previous, "distinct critical-line zero")
        h = ROOT*(B+1-2*ROOT*(tau*LOG).cos())/(arb(fmpq(1,4))+tau*tau)
        check(h>0,"positive coefficient")
        pairs.append((tau,2*h))
        previous=tau
    return pairs


def phase_cover():
    pairs = spectral_coefficients()
    lo,hi = 11,35
    check((P0*B).log()>lo and (M*M*B).log()<hi,"logarithmic coverage")
    k = CELLS_PER_UNIT
    radius = arb(fmpq(1,2*k))
    for j in range(lo*k,hi*k):
        t = arb(fmpq(2*j+1,2*k),radius)
        gain = sum(h*(1+(tau*t).cos()) for tau,h in pairs)
        check(gain>arb("0.020"),f"phase cell {j}/{k}")
    print(f"Positive phase gain > 0.020 using {len(pairs)} certified zero pairs on {(hi-lo)*k} closed cells in [11,35].")


def analytic_bounds():
    # The single-scale ramp at b vanishes because b<2. This is L_b/sqrt(b).
    mass = ROOT*((B+1)*(B+1).log()-(B-1)*(B-1).log()+2*(B-1)*LOG
                 -2*B*arb(2).log()-(B-1)*(arb.const_euler()+arb.pi().log()))
    drift = 2*ROOT*(ROOT-1)**2
    variation = 2*ROOT*(B-1)
    def error(m):
        return variation*arb(m).log()**2/(8*arb.pi()*arb(m).sqrt())
    check(mass<arb("0.126"),"mass")
    check(drift-error(isqrt(P1))>arb("0.107"),"intermediate analytic square bound")
    check(drift-mass-error(M)>arb("0.00003"),"infinite tail")
    check(arb("0.110")-mass+arb("0.020")>0,"middle ranges join")
    check(arb("0.107")-mass+arb("0.020")>0,"last finite range join")
    print("L < 0.126. Q(y)/sqrt(y) > 0.107 for y >= 10^10 under RH.")
    print(f"Full normalized deficit > 0.00003 for y >= {M*M} under RH.")


if __name__ == "__main__":
    finite_prime_ranges()
    phase_cover()
    analytic_bounds()
