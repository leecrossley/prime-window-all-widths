"""Certify the numerical steps in the unique-dilation proof using Arb balls."""
from flint import arb, ctx

ctx.prec = 144
gamma = arb.const_euler()
pi = arb.pi()


def check(condition, message):
    if not condition:
        raise ArithmeticError(message)


def F(b):
    """(B_b-L_b)/sqrt(b) on 1<b<=2."""
    r = b.sqrt()
    return (2*(r-1)**2 - (b+1)*(b+1).log() + (b-1)*(b-1).log()
            - 2*(b-1)*b.log() + 2*b*arb(2).log() + (b-1)*(gamma+pi.log()))


def main():
    lo,hi=arb("1.50342897722316"),arb("1.50342897722317")
    low, high=F(lo),F(hi)
    check(low < -arb("1.28e-15") and high > arb("2.59e-15"), "root bracket")
    two=arb(2)
    r=two.sqrt()
    fp=1-r+gamma+(pi/3).log()
    p=3*r-10
    c=2+gamma-(4*pi).log()
    rough=2*((r-1)/(r+1))**2-c
    check(fp>arb("0.209") and p < -5 and rough>arb("0.0126"), "endpoint signs")
    print("F(1.50342897722316) < -1.28e-15; F(1.50342897722317) > 2.59e-15.")
    print("F'(2) > 0.209; P(sqrt(2)) < -5.")
    print("The coarse positive margin at 2 exceeds 0.0126 after division by sqrt(b)*(sqrt(b)+1)^2.")


if __name__ == "__main__":
    main()
