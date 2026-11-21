# Eventual positivity at dilation 151/100

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410883.svg)](https://doi.org/10.5281/zenodo.24410883)

For the window from $y$ to $2.2801y$ (dilation $b=151/100$), the weighted sum of primes is eventually below its predicted size, provided the Riemann hypothesis (RH) holds. The shortfall, divided by $\sqrt y$, stays above $3/1000$ in the limit.

The margin is thin, and that is the point. This dilation sits only just above $3/2$, where the same family changes sign arbitrarily far out under RH and an independence hypothesis, so the two results show how sharply the behaviour switches.

## Statement

Assume RH and set $b=151/100$. Then $\liminf_{y\to\infty}D_b(y)/\sqrt y>3/1000$. In particular the prime sum on $(y,2.2801y)$ is eventually strictly below its main term.

## Setup

For $b>1$ and $y>0$:

- **Weight.** $W_b(t)$ equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere.
- **Prime sum.** $P_b(y)=\sum_{p\ {\rm prime}}(\log p)W_b(p/y)$, with natural logarithms. Its expected main term is $b(b-1)y\log b$.
- **Deficit.** $D_b(y)=b(b-1)y\log b-P_b(y)$.
- **Prime-square drift.** $B_b=2\sqrt b(\sqrt b-1)^2$, the asymptotic prime-square contribution.
- **Zero sum.** $L_b=\sum_\rho|H_b(\rho)|$ under RH, where $H_b(s)=(b^s-1)(b^s-b)/(s(s-1))$ is the Mellin transform of $W_b$.

## Proof

The general extremal bound gives $\liminf D_b(y)/\sqrt y\ge B_b-L_b$ under RH. [Derivation of the general extremal formula](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

Since $b\le2$, the ramp term $R(b)=\sum_n\Lambda(n)(b/n-1)_+$ of the general formula vanishes, so the closed formula uses only logarithms, $\pi$, Euler's constant and $\sqrt b$. The executable encloses these quantities with rational endpoints and proves $B_b-L_b>3/1000$. Numerically the endpoint is about $0.00312999$.

## Verification

The executable certifies the strict numerical inequality. The asymptotic implication uses the written explicit-formula argument. A later proof gives strict positivity for every real $y>0$ at this dilation under RH. [All-scale proof](https://github.com/leecrossley/prime-band-151-100-all-scales/blob/ce7d8f07b426e529c47c6748a5c919faf3d7d8cb/README.md)

```sh
bash verify-release.sh
```

This checks the [SHA-256 manifest](SHA256SUMS), reruns the arithmetic and compares its output with the [recorded result](evidence/check.txt). It needs Python 3.9 or later and only the standard library. The recorded run used Python 3.9.6 on arm64 macOS. [Latest release](https://github.com/leecrossley/prime-band-151-100/releases/latest)

## Related results

- [Positivity at every scale for this dilation](https://github.com/leecrossley/prime-band-151-100-all-scales/tree/v1.0.1)
- [Conditional oscillation at three halves](https://github.com/leecrossley/prime-band-three-halves-oscillation/tree/v1.0.1)
