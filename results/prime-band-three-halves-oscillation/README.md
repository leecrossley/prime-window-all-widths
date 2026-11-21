# Conditional sign changes at dilation three halves

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410827.svg)](https://doi.org/10.5281/zenodo.24410827)

Assuming the Riemann hypothesis (RH) and that the distinct positive ordinates of the nontrivial zeta zeros are linearly independent over the rationals, the weighted sum of primes between $y$ and $2.25y$ at dilation $b=3/2$ does not stay on one side of its predicted size. It is below the prediction at arbitrarily large $y$ and above it at arbitrarily large $y$.

This contrasts with dilations above the transition $b_\ast\approx1.5034$ in the same family, where RH forces a shortfall at all sufficiently large $y$. It also contrasts with the finite check at this dilation, which finds a shortfall at every scale up to 100,000.

## Statement

Assume RH and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros. At $b=3/2$ the ordinary-prime deficit takes both signs arbitrarily far out.

## Setup

For $b>1$ and $y>0$:

- **Weight.** $W_b(t)$ equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere.
- **Prime sum.** $P_b(y)=\sum_{p\ {\rm prime}}(\log p)W_b(p/y)$, with natural logarithms. Its expected main term is $b(b-1)y\log b$.
- **Deficit.** $D_b(y)=b(b-1)y\log b-P_b(y)$.
- **Prime-square drift.** $B_b=2\sqrt b(\sqrt b-1)^2$, the asymptotic prime-square contribution.
- **Zero sum.** $L_b=\sum_\rho|H_b(\rho)|$ under RH, where $H_b(s)=(b^s-1)(b^s-b)/(s(s-1))$ is the Mellin transform of $W_b$.

## Proof

The general extremal formula gives $\limsup D_b(y)/\sqrt y=B_b+L_b$ under RH and $\liminf D_b(y)/\sqrt y=B_b-L_b$ with the additional independence hypothesis. [Derivation of the general extremal formula](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

The executable proves $B_b-L_b<-0.0015$ and $B_b+L_b>0$ at $b=3/2$. Both signs therefore occur beyond every finite bound.

Under RH alone this argument gives a lower bound of negative sign, which leaves the occurrence of negative values open.

## Verification

The written phase argument uses Kronecker's theorem for the negative extremum. The executable checks the signs of the two constants. Finding an explicit negative value remains open here.

```sh
bash verify-release.sh
```

This checks the [SHA-256 manifest](SHA256SUMS), reruns the arithmetic and compares its output with the [recorded result](evidence/check.txt). It needs Python 3.9 or later and only the standard library. The recorded run used Python 3.9.6 on arm64 macOS. [Latest release](https://github.com/leecrossley/prime-band-three-halves-oscillation/releases/latest)

## Related result

[Finite certificate at the same dilation](https://github.com/leecrossley/prime-band-three-halves-finite/tree/v1.0.0)
