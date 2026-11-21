"""Independent rational proof of the transition bracket from the unshifted ramp.

Only Python's standard library is needed. Intervals below are pairs of exact
fractions. The trivial-zero series is summed directly, without its closed form.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt


def point(x):
    return Q(x),Q(x)


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def sub(a,b):
    return a[0]-b[1],a[1]-b[0]


def mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v),max(v)


def unit_log(x):
    z=(x-1)/(x+1)
    n=45
    t=2*sum(z**(2*k+1)/(2*k+1) for k in range(n))
    return t,t+2*z**(2*n+1)/((2*n+1)*(1-z*z))


@lru_cache(maxsize=None)
def logq(x):
    x=Q(x)
    if x <= 0:
        raise ValueError("positive logarithm argument required")
    j=0
    while x>=2:
        x/=2
        j+=1
    while x<1:
        x*=2
        j-=1
    return add(unit_log(x),mul(point(j),unit_log(Q(2))))


def atan_inv(d):
    z=Q(1,d)
    n=32
    partial=sum((-1)**k*z**(2*k+1)/(2*k+1) for k in range(n))
    following=partial+(-1)**n*z**(2*n+1)/(2*n+1)
    return min(partial,following),max(partial,following)


def constants():
    pi=sub(mul(point(16),atan_inv(5)),mul(point(4),atan_inv(239)))
    n=1024
    h=sum(Q(1,j) for j in range(1,n+1))
    approx=sub(point(h-Q(1,2*n)+Q(1,12*n**2)-Q(1,120*n**4)),logq(n))
    gamma=add(approx,(Q(0),Q(1,252*n**6)))
    lnpi=logq(pi[0])[0],logq(pi[1])[1]
    ln2pi=add(logq(2),lnpi)
    c0=sub(add(point(2),gamma),add(logq(2),ln2pi))
    return gamma,c0,ln2pi


def endpoint(b,gamma,c0,ln2pi):
    scale=10**35
    k=isqrt(b.numerator*scale**2//b.denominator)
    r=Q(k,scale),Q(k+1,scale)
    x=1/b**2
    n=78
    part=sum(x**j/(2*j*(2*j+1)) for j in range(1,n+1))
    T=part,part+x**(n+1)/((2*n+2)*(2*n+3)*(1-x))
    # For 1<b<=2 the von Mangoldt ramp is zero at b.
    S=add(sub(mul(point(b),sub(add(point(1),gamma),logq(b))),ln2pi),T)
    L_over_r=sub(mul(point(b+1),c0),mul(point(2),S))
    d=sub(r,point(1))
    return sub(mul(point(2),mul(d,d)),L_over_r)


if __name__=="__main__":
    gamma,c0,ln2pi=constants()
    left=endpoint(Q("1.50342897722316"),gamma,c0,ln2pi)
    right=endpoint(Q("1.50342897722317"),gamma,c0,ln2pi)
    if not (left[1]<-Q("1.28e-15") and right[0]>Q("2.59e-15")):
        raise ArithmeticError("direct ramp does not enclose the claimed signs")
    print("The 78-term trivial-zero series and exact rational interval arithmetic confirm both endpoint signs.")
