"""Verify the sign on every real scale through 100,000 for b=3/2."""


from bisect import bisect_right
from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt


if not __debug__:
    raise SystemExit("Run with Python assertions enabled.")

SCALE = 10**10


def log_unit_interval(t):
    """Enclose log(t), 1 <= t <= 2, by the atanh series and a geometric tail."""
    assert 1 <= t <= 2
    z = (t - 1) / (t + 1)
    n = 14
    low = 2 * sum(z ** (2*j + 1) / (2*j + 1) for j in range(n))
    return low, low + 2*z**(2*n+1) / ((2*n+1)*(1-z*z))


LOG2 = log_unit_interval(Q(2))


@lru_cache(maxsize=None)
def log_interval(t):
    t = Q(t)
    assert t > 0
    if t < 1:
        low, high = log_interval(1 / t)
        return -high, -low
    k = (t.numerator // t.denominator).bit_length() - 1
    lo, hi = log_unit_interval(t / 2**k)
    return k * LOG2[0] + lo, k * LOG2[1] + hi


def ceil_q(x):
    return -(-x.numerator // x.denominator)


def primes_to(bound):
    prime = bytearray([1]) * (bound + 1)
    prime[:2] = b'\0\0'
    for p in range(2, isqrt(bound) + 1):
        if prime[p]:
            prime[p*p::p] = bytes((bound-p*p)//p+1)
    return [p for p in range(2, bound+1) if prime[p]]


def small_intervals(primes, end, a=2, c=1):
    # At u = a²y the changes of slope for b = a/c occur at c²p, acp, a²p.
    n = a*a*end
    changes = [0]*(n+1)
    for p in primes[:bisect_right(primes,n//(c*c))]:
        coeff = ceil_q(SCALE*log_interval(p)[1]/p)
        for j, weight in ((c*c*p,a), (a*c*p,-a-c), (a*a*p,c)):
            if j <= n:
                changes[j] += weight*coeff
    slope = upper = 0
    main_slope = (a-c)*int(SCALE*log_interval(Q(a,c))[0])
    for j in range(1,n+1):
        upper += slope
        # The left side bounds a c² SCALE P_b(j/a²) from above.
        assert upper < j*main_slope, ("failed at b, u =", Q(a,c), j)
        slope += changes[j]
    print(f"At b={Q(a,c)}, certified strict sign at all {n} grid endpoints and by linearity for 0 < y <= {end}.")


def main():
    small_intervals(primes_to(225_000), 100_000, a=3, c=2)


if __name__ == "__main__":
    main()
