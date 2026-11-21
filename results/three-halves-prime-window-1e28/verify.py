"""A finite-height, unconditional certificate for the three-halves prime window."""
from bisect import bisect_right
from math import isqrt
from flint import acb, arb, ctx, fmpq

ctx.prec = 104
b = arb(fmpq(3, 2))
r = b.sqrt()
H = arb(3 * 10**12)

def check(test, reason):
    if not test:
        raise ArithmeticError(reason)


def prime_table(limit):
    flags = bytearray([1]) * (limit + 1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = bytes((limit - p*p)//p + 1)
    return [n for n, flag in enumerate(flags) if flag]


class WeightedSums:
    def __init__(self, primes, k):
        self.values = [p**k for p in primes]
        self.logs = [arb(0)]
        self.ratios = [arb(0)]
        for p, value in zip(primes, self.values):
            logp = arb(p).log()
            self.logs.append(self.logs[-1] + logp)
            self.ratios.append(self.ratios[-1] + logp / value)

    def at(self, nine_y):
        """Sum at y = nine_y / 9, using integer comparisons at each corner."""
        i, j, h = (bisect_right(self.values, nine_y // d) for d in (9, 6, 4))
        y = arb(fmpq(nine_y, 9))
        return (b*(self.logs[j]-self.logs[i]) - b*y*(self.ratios[j]-self.ratios[i])
                + b*b*y*(self.ratios[h]-self.ratios[j]) - (self.logs[h]-self.logs[j]))


def small_ranges():
    primes = prime_table(225000)
    direct = WeightedSums(primes, 1)
    centers = [6*p for p in primes if 6*p < 900000] + [900000]
    for point in centers:
        check(b*(b-1)*b.log()*arb(fmpq(point, 9)) > direct.at(point), "ordinary prime centre")
    sq = WeightedSums([p for p in primes if p <= 150000], 2)
    low, high = 9*10**5, 9*10**10
    knots = sorted({low, high} | {u for z in sq.values for u in (4*z, 6*z, 9*z) if low < u < high})
    old = sq.at(low)
    for left, right in zip(knots, knots[1:]):
        new = sq.at(right)
        m = isqrt((left+right)//2)
        for u, value in ((left, old), (right, new)):
            check(value > arb(fmpq(11*(u+m*m), 600*m)), "square tangent")
        old = new
    print(f"{len(centers)} ordinary-prime centres and endpoints, {len(knots)-1} square intervals.")


def inverse_square_tail(T):
    """Brent, Platt and Trudgian, Example 1, for positive ordinates."""
    return ((1 + (T/(2*arb.pi())).log())/(2*arb.pi()*T)
            + (arb("0.14") + arb("0.56")*T.log())/T**2)


# Strict upward integer bounds for the two constants in KLN Table 1, allowing equality.
# The exponent sigma is given in hundredths. All tabulated k in these rows equal 1.
KLN = [(80,7,4), (85,9,4), (86,10,4), (87,10,4), (88,11,4), (89,11,4),
       (90,12,4), (91,13,4), (92,13,4), (93,14,4), (94,14,4), (95,15,4),
       (96,16,4), (97,16,4), (98,17,4), (99,18,4)]


def high_zero_bounds():
    logH = H.log()
    total = inverse_square_tail(H)
    tails = []
    for numerator, A, C in KLN:
        s = arb(fmpq(numerator, 100))
        p = arb(fmpq(8, 3))*(1-s)
        q = 5-2*s
        alpha = 2-p-q/logH
        beta = 2-2/logH
        check(alpha > 0 and beta > 0, "integral convergence")
        bound = (2*A*H**(p-2)*logH**q/alpha + 2*C*logH**2/(H**2*beta))
        tails.append((s, bound.min(total/2)))
    return total, tails


def high_cost(z, total, tails):
    # Every low-beta ordinate contributes the baseline; the positive excess
    # is bounded by a finite telescoping sum over cuts in beta.
    powers = [z**(s-arb("0.5")) for s, bound in tails] + [z.sqrt()]
    excess = (powers[0]-1)*total/2
    for k, (_, bound) in enumerate(tails):
        excess += (powers[k+1]-powers[k])*bound
    return 4*b*(b+1)*(total + excess)


def spectral_check():
    ordinates = acb.zeta_zeros(1,649)
    check(arb(1000).zeta_nzeros() == 649, "all ordinates below cutoff")
    previous=arb(0)
    modes=[]
    for zero in ordinates:
        t=zero.imag
        check(zero.real == arb("0.5") and previous < t and t < 1000, "ordered critical zero")
        a=2*r*((r-1)**2 + 4*r*(t*b.log()/2).sin()**2)/(arb("0.25")+t*t)
        modes.append((t,a))
        previous=t
    mass=sum(a for t,a in modes) + 2*r*(r+1)**2*inverse_square_tail(arb(1000))
    check(mass < arb("0.13121"), "verified-height mass")
    total, tails = high_zero_bounds()
    check(high_cost(arb(10)**28,total,tails) < arb("0.0505"), "outer zero error")
    drift=2*r*(r-1)**2
    variation=2*r*(b-1)
    for x, floor in ((arb(10)**5,arb("0.10")), (arb(10)**10,arb("0.123"))):
        check(drift-variation*x.log()**2/(8*arb.pi()*x.sqrt()) > floor,"prime-square floor")
    last_x=b*arb(10)**14
    check(arb("4.92")*(last_x/last_x.log()).sqrt() < H, "Büthe height condition")
    selected=modes[:128]
    factor=arb(10).log()/50
    cells=120
    logb=b.log()
    for j in range(250,1400):
        square_floor=arb("0.10") if j < 1000 else arb("0.123")
        error=high_cost(((j+1)*factor).exp(),total,tails)
        threshold=mass+error-square_floor
        for k in range(cells):
            left=(j*cells+k)*factor/cells+logb
            right=(j*cells+k+1)*factor/cells+logb
            time=left.union(right)
            gain=arb(0)
            for t,a in selected:
                gain+=a*(1+(t*time).cos()).max(arb(0))
                if gain>threshold:
                    break
            else:
                raise ArithmeticError(f"unresolved phase cell {j}/{k}")
    print("649 certified ordinates and 16 zero-density cuts bound all uncomputed zeros.")
    print(f"{1150*cells} closed logarithmic cells verify every scale from 10^5 to 10^28.")


if __name__ == "__main__":
    small_ranges()
    spectral_check()
