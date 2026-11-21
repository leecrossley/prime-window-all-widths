"""Two-sided logarithmic density bounds under RH and independence."""
from flint import arb, acb, ctx, fmpq

ctx.prec = 96
B = arb(fmpq(3,2))
R = B.sqrt()


def require(condition, detail):
    if not condition:
        raise ArithmeticError(detail)


def drift_and_mass():
    drift=2*R*(R-1)**2
    L=R*((B+1)*(B+1).log()-(B-1)*(B-1).log()+2*(B-1)*B.log()
         -2*B*arb(2).log()-(B-1)*(arb.const_euler()+arb.pi().log()))
    return drift,L


def weights(n):
    result=[]
    previous=arb(0)
    for j in range(1,n+1,500):
        batch=acb.zeta_zeros(j,min(500,n-j+1))
        for z in batch:
            t=z.imag
            require(z.real==arb(fmpq(1,2)) and t>previous,"distinct positive zero")
            h=2*R*((R-1)**2+4*R*(t*B.log()/2).sin()**2)/(arb(fmpq(1,4))+t*t)
            require(h>0,"positive coefficient")
            result.append(h)
            previous=t
    require(len(result)==n,"all requested ordinates")
    return result

def main():
    drift,L=drift_and_mass()
    a=weights(15000)
    d=L-drift
    tilt=1_170_000
    log_upper=tilt*d+sum((tilt*w).bessel_i(0,scaled=True).log() for w in a[:5000])
    require(log_upper < -1012*arb(10).log(),"upper density")
    budget=arb("0.0004")
    require(2*(L-sum(a))+budget<d,"open negative box")
    log_lower=arb(0)
    for w in a:
        ratio=(2*budget/(len(a)*w)).sqrt()/arb.pi()
        log_lower+=ratio.min(arb(1)).log()
    require(log_lower> -10000*arb(10).log(),"positive box volume")
    print("5000 certified pairs give density below 10^(-1012).")
    print("15000 certified pairs give density above 10^(-10000) with angle budget 0.0004.")


if __name__=="__main__":
    main()
