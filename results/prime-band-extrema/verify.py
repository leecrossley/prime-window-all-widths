"""Exact enclosures for the two dyadic spectral endpoints."""


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


def atan_inv(q, n):
    x = Q(1, q)
    a = sum((-1)**j * x**(2*j+1)/(2*j+1) for j in range(n))
    b = a + (-1)**n * x**(2*n+1)/(2*n+1)
    return min(a, b), max(a, b)


def fundamental_constants():
    """Enclose pi (Machin), Euler's constant (harmonic sums), and the zero mass."""
    a_lo, a_hi = atan_inv(5, 10)
    b_lo, b_hi = atan_inv(239, 4)
    pi_lo, pi_hi = 16*a_lo - 4*b_hi, 16*a_hi - 4*b_lo
    n = 100_000
    precision = 10**14
    harmonic_lo = Q(sum(precision//i for i in range(1, n+1)), precision)
    harmonic_hi = Q(sum((precision+i-1)//i for i in range(1, n+1)), precision)
    gamma_lo = harmonic_lo - log_interval(n+1)[1]
    gamma_hi = harmonic_hi - log_interval(n)[0]
    mass_lo = 2 + gamma_lo - log_interval(4*pi_hi)[1]
    mass_hi = 2 + gamma_hi - log_interval(4*pi_lo)[0]
    assert 0 < mass_lo <= mass_hi
    return pi_lo, pi_hi, gamma_lo, gamma_hi, mass_lo, mass_hi


def spectral(b):
    """Rational enclosures for B-L and B+L when the ramp R(b) vanishes."""
    assert 1 < b <= 2
    pi_lo, pi_hi, gamma_lo, gamma_hi, mass_lo, mass_hi = fundamental_constants()
    root = isqrt(b.numerator * SCALE**2 // b.denominator)
    r_lo, r_hi = Q(root, SCALE), Q(root+1, SCALE)
    assert 1 < r_lo and r_lo**2 <= b < r_hi**2
    log_lo, log_hi = log_interval(b)
    minus_lo, minus_hi = log_interval(b-1)
    plus_lo, plus_hi = log_interval(b+1)
    tail_lo = 1 + log_lo + ((b-1)*minus_lo - (b+1)*plus_hi)/2
    tail_hi = 1 + log_hi + ((b-1)*minus_hi - (b+1)*plus_lo)/2
    s_lo = b*(1 + gamma_lo - log_hi) - log_interval(2*pi_hi)[1] + tail_lo
    s_hi = b*(1 + gamma_hi - log_lo) - log_interval(2*pi_lo)[0] + tail_hi
    amplitude_lo = (b+1)*mass_lo - 2*s_hi
    amplitude_hi = (b+1)*mass_hi - 2*s_lo
    assert 0 < amplitude_lo <= amplitude_hi
    l_lo, l_hi = r_lo*amplitude_lo, r_hi*amplitude_hi
    drift_lo, drift_hi = 2*r_lo*(r_lo-1)**2, 2*r_hi*(r_hi-1)**2
    return drift_lo-l_hi, drift_hi-l_lo, drift_lo+l_lo, drift_hi+l_hi


def main():
    lo, hi, top_lo, top_hi = spectral(Q(2))
    assert Q(2198,10000) < lo <= hi < Q(2201,10000)
    assert Q(7504,10000) < top_lo <= top_hi < Q(7508,10000)
    print("0.2198 < B_2 - L_2 < 0.2201.")
    print("0.7504 < B_2 + L_2 < 0.7508.")


if __name__ == "__main__":
    main()
