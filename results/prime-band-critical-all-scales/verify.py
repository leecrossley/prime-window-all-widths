from bisect import bisect_right
from math import isqrt
from flint import acb,arb,ctx,fmpq

ctx.prec=112
Y0,Y1=10**5,10**10


def require(condition, message):
    if not condition: raise ArithmeticError(message)


def primes(n):
    flags=bytearray(b"\x01")*(n+1)
    flags[:2]=b"\x00\x00"
    for v in range(2,isqrt(n)+1):
        if flags[v]:flags[v*v::v]=bytes((n-v*v)//v+1)
    return [j for j in range(2,n+1) if flags[j]]


class Prefix:
    def __init__(self, ps, power):
        self.terms=[p**power for p in ps]
        self.L=[arb(0)]
        self.T=[arb(0)]
        for p,n in zip(ps,self.terms):
            v=arb(p).log()
            self.L.append(self.L[-1]+v)
            self.T.append(self.T[-1]+v/n)

    def band(self,num,den,a,c):
        y=arb(fmpq(num,den))
        b=arb(fmpq(a,c))
        i,j,k=(bisect_right(self.terms,q) for q in (num//den,num*a//(den*c),num*a*a//(den*c*c)))
        return (b*(self.L[j]-self.L[i])-b*y*(self.T[j]-self.T[i])
                +b*b*y*(self.T[k]-self.T[j])-(self.L[k]-self.L[j]))


def check_primes():
    ps=primes(228010)
    direct=Prefix(ps,1)
    count=0
    # Use 128 slices of [3/2,151/100]. The lower main term and upper
    # weight prove every interior parameter with a single comparison.
    c=12800
    for l in range(19200,19328):
        a=l+1
        b=arb(fmpq(l,c))
        slope=b*(b-1)*b.log()
        centers=[(p*c,a) for p in ps[:bisect_right(ps,Y0*a//c)]]
        for num,den in centers+[(Y0,1)]:
            require(slope*arb(fmpq(num,den))-direct.band(num,den,a,c)>arb("0.055"),"small prime centre")
        count+=len(centers)+1
    print(f"{count} prime centres and endpoints across 128 dilation slices. All finite endpoint deficits exceed 0.055.")
    # At any b>=3/2 the prime-square sum dominates its value at exactly 3/2.
    a,c=3,2
    sq=Prefix([p for p in ps if p<=150000],2)
    lo,hi=a*a*Y0,a*a*Y1
    knots=sorted({lo,hi}|{u for n in sq.terms for u in (c*c*n,a*c*n,a*a*n) if lo<u<hi})
    x=sq.band(knots[0],a*a,a,c)
    for left,right in zip(knots,knots[1:]):
        z=sq.band(right,a*a,a,c)
        m=isqrt((left+right)//2)
        for u,q in ((left,x),(right,z)):
            require(q>arb(fmpq(11*(u+m*m),100*2*a*m)),"square tangent")
        x=z
    print(f"Q(y)>0.11 sqrt(y) throughout {len(knots)-1} prime-square intervals.")


def analytic_and_phases():
    lo,hi=arb("1.50342"),arb("1.51")
    bs=lo.union(hi)
    r=bs.sqrt()
    drift2=2*r*(r-1)**2
    V2=2*r*(bs-1)
    k=bs.root(3)
    drift3=arb("1.5")*(k-1)*(bs-k)
    V3=arb("1.5")*(bs-k+bs*k-k*k)
    R=arb(96)
    pi=arb.pi()
    err2=V2*(R/2)**2*(-R/12).exp()/(8*pi)
    err3=V3*(R/3)**2*(-R/6).exp()/(8*pi)
    require(drift3-err2-err3>arb("0.035"),"uniform cube remainder")
    require(drift2<arb("0.13") and V2*arb(Y0).log()**2/(8*pi*arb(Y0).sqrt())<arb("0.021"),"middle constants")
    require((lo*Y0).log()>11 and hi.log()+R<97,"logarithmic coverage")
    require(arb("0.110")+arb("0.022")>arb("0.13"),"middle sums")
    zs=acb.zeta_zeros(1,64)
    coefficients=[]
    previous=arb(0)
    for z in zs:
        t=z.imag
        require(z.real==arb("0.5") and t>previous,"certified zeta zero")
        # h is monotone in sin² for fixed r. Take a uniform lower bound
        # allowing the whole range of b in the coefficient's phase.
        w=2*r*((r-1)**2+4*r*(t*bs.log()/2).sin()**2)/(arb("0.25")+t*t)
        weight=w.lower().max(arb(0))
        coefficients.append((t,weight))
        previous=t
    n=1000
    radius=arb(fmpq(1,2*n))
    for j in range(11*n,97*n):
        time=arb(fmpq(2*j+1,2*n),radius)
        g=sum(w*(1+(t*time).cos()) for t,w in coefficients)
        require(g>arb("0.022"),f"uniform phase cell {j}")
    print("64 zero pairs give phase gain > 0.022 on 86000 closed time cells for the whole dilation range.")
    print("The square loss is < 0.021 sqrt(y) beyond 10^10. At log(y)>=96 the net cube margin exceeds 0.035 y^(1/3).")


if __name__=="__main__":
    check_primes()
    analytic_and_phases()
