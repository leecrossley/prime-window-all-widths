"""Finite Fourier majorant for the 3/2 deficit under RH."""
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
    loss=L-drift
    require(0<loss and loss<arb("0.00163"),"distance to the lower endpoint")
    dimension,degree,tilt=500,500_000,150_000
    logarithm=tilt*loss
    for a in weights(dimension):
        v=tilt*a
        u=-arb(degree)*(-2*v/degree).expm1()/2
        require(u>0,"positive polynomial parameter")
        logarithm+=u.bessel_i(0,scaled=True).log()
    require(logarithm < -232*arb(10).log(),"logarithmic upper bound")
    print("500 isolated zero pairs. Degree 500000 and tilt 150000.")
    print("The certified logarithmic upper density is below 10^(-232).")


if __name__=="__main__":
    main()
