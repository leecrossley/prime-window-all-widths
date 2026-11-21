"""Independent bounds from three signed ramps, concave chords and a trivial-zero series."""

import json
from bisect import bisect_right
from fractions import Fraction
from math import isqrt
from pathlib import Path

from flint import arb, ctx, fmpq

ctx.prec = 160
SCALE = 10**25
START, END = 10**5, 10**10
# (left width, right width, separation between power sum and spectral mass)
ROWS = [(151, 154, 127), (154, 165, 140), (165, 185, 195), (185, 227, 300)]


def check(condition, details):
    if not condition:
        raise ArithmeticError(details)


def sieve(limit):
    flags = bytearray([1]) * (limit+1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit)+1):
        if flags[p]:
            flags[p*p::p] = bytes((limit-p*p)//p+1)
    return [p for p in range(2, limit+1) if flags[p]]


def scaled(value):
    value = SCALE*value
    return int(value.lower().floor().unique_fmpz()), int(value.upper().ceil().unique_fmpz())


class RampSums:
    def __init__(self, pairs):
        self.nodes = []
        self.lower_log, self.upper_log = [0], [0]
        self.lower_harm, self.upper_harm = [0], [0]
        for n, p in pairs:
            self.nodes.append(n)
            log_p = arb(p).log()
            lo, hi = scaled(log_p)
            self.lower_log.append(self.lower_log[-1]+lo)
            self.upper_log.append(self.upper_log[-1]+hi)
            lo, hi = scaled(log_p/n)
            self.lower_harm.append(self.lower_harm[-1]+lo)
            self.upper_harm.append(self.upper_harm[-1]+hi)

    def directed(self, upper):
        # Bound b R(y) - (b+1) R(by) + R(b²y), where R(x)=x H(x)-theta(x).
        plus_h, plus_l = self.upper_harm, self.lower_log
        minus_h, minus_l = self.lower_harm, self.upper_log
        if not upper:
            plus_h, plus_l, minus_h, minus_l = minus_h, minus_l, plus_h, plus_l
        def value(z, d, a, c, i, j, k):
            first = a*c*(z*plus_h[i]-d*plus_l[i])
            middle = (a+c)*(a*z*minus_h[j]-c*d*minus_l[j])
            last = a*a*z*plus_h[k]-c*c*d*plus_l[k]
            return first-middle+last
        return value


def ordinary(primes):
    data = json.loads(Path(__file__).with_name("widths.json").read_text())
    c = data["denominator"]
    cells = data["intervals"]
    check(c == 10**5 and data["limit"] == START, "scale or denominator")
    check(cells[0][0] == 151000 and cells[-1][1] == 227000, "width endpoints")
    check(all(l < h for l, h in cells) and
          all(x[1] == y[0] for x, y in zip(cells, cells[1:])), "width cover")
    upper = RampSums((p, p) for p in primes).directed(upper=True)
    count = 0
    for low, a in cells:
        slope = low*(low-c)*scaled(arb(fmpq(low, c)).log())[0]
        middle = bisect_right(primes, START*a//c)
        for j in range(middle):
            p = primes[j]
            z = p*c
            i = bisect_right(primes, z//a)
            k = bisect_right(primes, p*a//c)
            if upper(z, a, a, c, i, j+1, k) >= slope*z:
                raise ArithmeticError(("three ramps", low, a, p))
        q = upper(START, 1, a, c, bisect_right(primes, START), middle,
                  bisect_right(primes, START*a*a//(c*c)))
        check(q < slope*START, "upper scale endpoint")
        count += middle+1
    check(count == data["endpoints"], "endpoint count")
    print(f"Three-ramp check: {count} positive endpoints across {len(cells)} width intervals.")


def middle(primes):
    pairs = []
    for p in primes:
        for n in (p*p, p*p*p):
            if 10000*n <= END*227**2:
                pairs.append((n, p))
    ramps = RampSums(sorted(pairs))
    lower = ramps.directed(upper=False)
    ns = ramps.nodes
    for a, _, h in ROWS:
        c = 100
        d = a*a
        start, end = d*START, d*END
        # Between outer support contacts the remaining slope jumps are downward.
        contacts = {start, end}
        for n in ns:
            for z in (c*c*n, d*n):
                if start < z < end:
                    contacts.add(z)
        factor2 = (h*c*c*a*SCALE)**2
        last = None
        for z in sorted(contacts):
            i = bisect_right(ns, z//d)
            j = bisect_right(ns, z//(a*c))
            k = bisect_right(ns, z//(c*c))
            value = 1000*lower(z, d, a, c, i, j, k)
            check(value > 0 and value**2 > factor2*z, ("chord endpoint", a, z))
            if last is not None:
                u, v = last
                check(2*v*value > factor2*(u+z), ("chord cross term", a, u, z))
            last = z, value
        print(f"Width {a}/100: {len(contacts)-1} lower chords clear the same square-root floor.")


def ball(rational):
    return arb(fmpq(rational.numerator, rational.denominator))


def atan_bounds(q):
    x = Fraction(1, q)
    even = sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(20))
    return ball(even).union(ball(even+x**41/41))


def constants():
    pi = 16*atan_bounds(5)-4*atan_bounds(239)
    n = 64
    harmonic = sum(Fraction(1, j) for j in range(1, n+1))
    base = harmonic-Fraction(1, 2*n)+Fraction(1, 12*n**2)-Fraction(1, 120*n**4)
    gamma = ball(base).union(ball(base+Fraction(1,252*n**6)))-arb(n).log()
    return gamma, pi


def spectral():
    gamma, pi = constants()
    c0 = 2+gamma-(4*pi).log()
    check(c0 < arb("0.0462"), "total inverse-square mass")
    for low, high, h in ROWS:
        # The width grid is independent of the adaptive closed-form calculation.
        # Start with intervals of length 1/1000 and halve any inconclusive one.
        for j in range(10*low, 10*high):
            pending = [(Fraction(j, 1000), Fraction(j+1, 1000))]
            while pending:
                left, right = pending.pop()
                b = ball(left).union(ball(right))
                ramp = (b/2-1)*arb(2).log() if left >= 2 else arb(0)
                x = 1/b**2
                power, partial = arb(1), arb(0)
                for k in range(1, 25):
                    power *= x
                    partial += power/(2*k*(2*k+1))
                tail = power*x/(50*51*(1-x))
                trivial = partial+arb(0).union(tail)
                s = ramp+b*(1+gamma-b.log())-(2*pi).log()+trivial
                r = b.sqrt()
                mass = r*((b+1)*c0-2*s)
                if mass < arb(fmpq(h,1000)) and 2*r*(r-1)**2-mass > arb("0.003"):
                    continue
                check(right-left > Fraction(1, 10**7), ("ramp series", left, right))
                mid = (left+right)/2
                pending.extend([(left, mid), (mid, right)])
    # Rational majorants for the fourth root, the error constant and the tail.
    q = Fraction
    check(q(123,100)**4 > q(227,100), "fourth root")
    check(2*q(223,100)*q(350,100)*q(462,10000) < q(73,100), "even powers")
    check(316**4 <= END and isqrt(END)**2 == END, "tail endpoints")
    check(q(3,1000)-q(73,31600)-q(1,(END-1)*isqrt(END)) > q(68,100000),
          "infinite scales")
    print("The unshifted von Mangoldt ramp and 24 trivial-zero terms confirm all spectral bounds.")


if __name__ == "__main__":
    ps = sieve(START*227**2//100**2)
    spectral()
    middle(ps)
    ordinary(ps)
