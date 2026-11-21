"""Check the ordinary-prime all-width inequalities, including real-parameter coverage."""

from bisect import bisect_right
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from functools import lru_cache
from itertools import accumulate
from math import isqrt
from pathlib import Path
import json

from flint import arb, ctx, fmpq

getcontext().prec = 65
ctx.prec = 144
UNIT = 10**24
SMALL = 100_000
LARGE = 10**10
# Width endpoints are hundredths. Each square-plus-cube floor is in thousandths.
BANDS = [(151, 154, 127), (154, 165, 140), (165, 185, 195), (185, 227, 300)]


def require(condition, label):
    if not condition:
        raise ArithmeticError(label)


def primes_through(limit):
    primes = [2]
    for candidate in range(3, limit + 1, 2):
        for divisor in primes:
            if divisor * divisor > candidate:
                primes.append(candidate)
                break
            if candidate % divisor == 0:
                break
    return primes


@lru_cache(maxsize=None)
def log_bounds(value):
    # Decimal.ln is correctly rounded. The adjacent values enclose the exact log.
    center = Decimal(value).ln()
    lo = (center.next_minus() * UNIT).to_integral_value(rounding=ROUND_FLOOR)
    hi = (center.next_plus() * UNIT).to_integral_value(rounding=ROUND_CEILING)
    # Multiplication by UNIT is an exact decimal exponent shift.
    return int(lo), int(hi)


class Window:
    def __init__(self, terms, upper):
        # terms is an increasing sequence of (n, p) with weight log(p) at n.
        self.nodes = [n for n, _ in terms]
        bounds = [log_bounds(p) for _, p in terms]
        log_lo, log_hi = zip(*bounds)
        inv_lo = [v // n for v, n in zip(log_lo, self.nodes)]
        inv_hi = [-(-v // n) for v, n in zip(log_hi, self.nodes)]
        lp, lm, hp, hm = (log_hi, log_lo, inv_hi, inv_lo) if upper else (
            log_lo, log_hi, inv_lo, inv_hi
        )
        self.lp, self.lm, self.hp, self.hm = (
            [0, *accumulate(values)] for values in (lp, lm, hp, hm)
        )

    def bound(self, z, d, a, c, i, j, k):
        """Directed bound multiplied by c² d UNIT at y=z/d, b=a/c."""
        left = a*c*(d*(self.lp[j]-self.lp[i]) - z*(self.hm[j]-self.hm[i]))
        right = a*a*z*(self.hp[k]-self.hp[j]) - c*c*d*(self.lm[k]-self.lm[j])
        return left + right


def small_scales(primes):
    data = json.loads(Path(__file__).with_name("widths.json").read_text())
    cells, c = data["intervals"], data["denominator"]
    require(c == 100000 and data["limit"] == SMALL, "small scale domain")
    require(cells[0][0] == 151000 and cells[-1][1] == 227000, "width endpoints")
    prime_sum = Window([(p, p) for p in primes], upper=True)
    count = 0
    previous = cells[0][0]
    for low, a in cells:
        require(previous == low < a, "width partition")
        previous = a
        slope = low*(low-c)*log_bounds(Decimal(low) / c)[0]
        middle = bisect_right(primes, SMALL*a//c)
        for j in range(middle):
            p = primes[j]
            z = p*c
            i = bisect_right(primes, z//a)
            k = bisect_right(primes, p*a//c)
            if slope*z <= prime_sum.bound(z, a, a, c, i, j+1, k):
                raise ArithmeticError(("prime centre", low, a, p))
        q = prime_sum.bound(SMALL, 1, a, c, bisect_right(primes, SMALL),
                            middle, bisect_right(primes, SMALL*a*a//(c*c)))
        require(slope*SMALL > q, "last scale")
        count += middle + 1
    require(count == data["endpoints"], "small scale count")
    print(f"Ordinary primes: {len(cells)} closed width intervals, {count} positive endpoints.")


def middle_scales(primes):
    limit = (227**2 * LARGE) // 100**2
    terms = sorted((p**k, p) for k in (2, 3) for p in primes if p**k <= limit)
    powers = Window(terms, upper=False)
    ns = powers.nodes
    for a, _, height in BANDS:
        c = 100
        d = a*a
        low, high = d*SMALL, d*LARGE
        knots = {low, high}
        for multiplier in (a*a, a*c, c*c):
            knots.update(n*multiplier for n in ns if low < n*multiplier < high)
        knots = sorted(knots)
        def at(z):
            i, j, k = (bisect_right(ns, z//v) for v in (a*a, a*c, c*c))
            return powers.bound(z, d, a, c, i, j, k)
        qleft = at(low)
        for left, right in zip(knots, knots[1:]):
            qright = at(right)
            m = isqrt(left)
            factor = height * a * c*c * UNIT
            for z, q in ((left, qleft), (right, qright)):
                require(2000*m*q > factor*(z+m*m), ("square-root tangent", a, z))
            qleft = qright
        print(f"Width {a}/100: squares and cubes exceed {height}/1000 times sqrt(y) on {len(knots)-1} affine intervals.")


def analytic_bounds():
    gamma = arb.const_euler()
    log_pi = arb.pi().log()
    log_two = arb(2).log()
    for low, high, height in BANDS:
        # Split at 2 so the sole prime-power ramp in this domain is exact.
        pending = [(1000*low, 1000*high)]
        if low < 200 < high:
            pending = [(1000*low, 200000), (200000, 1000*high)]
        while pending:
            left, right = pending.pop()
            b = arb(fmpq(left, 100000)).union(arb(fmpq(right, 100000)))
            r = b.sqrt()
            twice_ramp = (b-2)*log_two if left >= 200000 else arb(0)
            mass = r*((b+1)*(b+1).log() - (b-1)*(b-1).log() + 2*(b-1)*b.log()
                      - 2*b*log_two - (b-1)*(gamma+log_pi) - twice_ramp)
            drift = 2*r*(r-1)**2
            if mass < arb(fmpq(height, 1000)) and drift-mass > arb("0.003"):
                continue
            mid = (left+right)//2
            require(left < mid < right, ("spectral enclosure", left, right))
            pending.extend([(left, mid), (mid, right)])
    bmax = arb(fmpq(227,100))
    fourth_root = bmax.root(4)
    c0 = 2+gamma-(4*arb.pi()).log()
    error = 2*(fourth_root+1)*(fourth_root+bmax)*c0
    require(error < arb("0.73"), "even-power zero contribution")
    y = arb(LARGE)
    require(arb("0.003") - arb("0.73")/y.root(4) - 1/((y-1)*y.sqrt())
            > arb("0.00068"), "infinite-range floor")
    print("Four zero-mass ceilings hold. The drift margin exceeds 0.003 on [1.51,2.27].")
    print("At every y >= 10^10 the even-power bound leaves more than 0.00068 sqrt(y).")


if __name__ == "__main__":
    # This range includes all needed primes for both the direct and higher-power checks.
    ps = primes_through(227**2 * SMALL // 100**2)
    analytic_bounds()
    middle_scales(ps)
    small_scales(ps)
