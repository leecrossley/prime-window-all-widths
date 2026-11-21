"""Cross-check the density calculation using individual xi sign intervals."""
from pathlib import Path
from flint import acb, arb, ctx, fmpq

ctx.prec=100
b=arb("1.5")
r=b.sqrt()


def require(condition,detail):
    if not condition:
        raise ArithmeticError(detail)


def original_ramp():
    g=arb.const_euler()
    c0=2+g-(4*arb.pi()).log()
    x=1/b**2
    n=64
    T=sum(x**j/(2*j*(2*j+1)) for j in range(1,n+1))
    tail=x**(n+1)/((2*n+2)*(2*n+3)*(1-x))
    S=(b*(1+g-b.log())-(2*arb.pi()).log()+T).union(b*(1+g-b.log())-(2*arb.pi()).log()+T+tail)
    return 2*r*(r-1)**2,r*((b+1)*c0-2*S)


def real_xi(t):
    s=acb(arb("0.5"),t)
    return (s*(s-1)/2 * (-s*arb.pi().log()/2).exp() * (s/2).gamma() * s.zeta()).real


def enclose_weights(n):
    file=Path(__file__).resolve().with_name("zero-centres.txt")
    centres=[int(s) for s in file.read_text().splitlines()]
    require(len(centres)==n,"number of centres")
    weights=[]
    previous=arb(0)
    step=arb("0.0000001")
    for k in centres:
        mid=arb(fmpq(k,10**9))
        left,right=mid-step,mid+step
        a,z=real_xi(left),real_xi(right)
        require(left>previous and ((a<0 and z>0) or (a>0 and z<0)),"xi sign change")
        t=left.union(right)
        w=2*r*(b+1-2*r*(t*b.log()).cos())/(arb("0.25")+t*t)
        require(w>0,"positive weight")
        weights.append(w)
        previous=right
    require((previous+arb("0.00001")).zeta_nzeros()==n,"zero count above all intervals")
    return weights


def scaled_bessel_upper(u):
    # Positive series exp(-u) I0(u). After K>=2u, consecutive ratios
    # are below 1/16, so 16/15 of the first omitted term bounds the tail.
    K=max(40,int((2*u).upper().floor().unique_fmpz())+1)
    w=(-u).exp()
    total=w
    x=(u/2)**2
    for j in range(1,K+1):
        w=w*x/(j*j)
        total+=w
    first_omitted=w*x/((K+1)**2)
    return total+arb(fmpq(16,15))*first_omitted

def main():
    drift,L=original_ramp()
    ws=enclose_weights(15000)
    tilt=1_170_000
    product=(tilt*(L-drift)).exp()
    for w in ws[:5000]:
        product*=scaled_bessel_upper(tilt*w)
    require(product<arb(10)**(-1012),"independent upper density")
    # Work directly with the angular box volume rather than its logarithm.
    prob=arb(1)
    e=arb("0.0004")
    require(e+2*(L-sum(ws))<L-drift,"negative throughout the box")
    for w in ws:
        prob*=((2*e/(len(ws)*w)).sqrt()/arb.pi()).min(arb(1))
    require(prob>arb(10)**(-10000),"independent lower density")
    print("15000 xi-sign intervals and a zero count establish the trial ordinates.")
    print("The direct box volume and positive Bessel series confirm both density bounds.")


if __name__=="__main__":
    main()
