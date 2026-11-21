# Extremal constants for a family of weighted prime sums

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410840.svg)](https://doi.org/10.5281/zenodo.24410840)

Under the Riemann hypothesis (RH), weighted sums of primes between $y$ and $b^2y$ tend to fall short of what the prime number theorem predicts. This repository works out how far the shortfall can swing. Under RH, once divided by $\sqrt y$, it is confined in the long run between two explicit constants, $B_b-L_b$ and $B_b+L_b$.

$B_b$ is a steady push from the prime squares, which a sum over primes leaves out. $L_b$ is the most the zeta zeros can push in either direction. Both have closed forms, so the sign of the computable number $B_b-L_b$ shows what is known at each dilation $b$. If $B_b-L_b>0$, RH gives a shortfall at all sufficiently large $y$. If it is negative, RH together with rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros gives both signs arbitrarily far out. At the unique critical width the normalised liminf is zero under RH and the independence hypothesis. A later proof gives $D_b(y)>0$ for every real $y>0$ at that width under RH alone. [Critical-width proof](https://github.com/leecrossley/prime-band-critical-all-scales/blob/99f0a80a6f332a43942344c470ab410eba0eb2f7/README.md)

## Statement

Fix $b>1$ and assume RH. Then

$$
 \limsup_{y\to\infty}\frac{D_b(y)}{\sqrt y}=B_b+L_b,\qquad
 \liminf_{y\to\infty}\frac{D_b(y)}{\sqrt y}\ge B_b-L_b.
$$

If the distinct positive ordinates of nontrivial zeta zeros are linearly independent over $\mathbb Q$, the liminf equals $B_b-L_b$.

Here, for $y>0$:

- $W_b(t)$ equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere.
- $P_b(y)=\sum_{p\ {\rm prime}}(\log p)W_b(p/y)$, with natural logarithms. Its expected main term is $b(b-1)y\log b$.
- $D_b(y)=b(b-1)y\log b-P_b(y)$ is the deficit.
- $B_b=2\sqrt b(\sqrt b-1)^2$ is the asymptotic prime-square contribution.
- $L_b=\sum_\rho|H_b(\rho)|$ over the nontrivial zeta zeros with multiplicity, where $H_b(s)=(b^s-1)(b^s-b)/(s(s-1))$ is the Mellin transform of $W_b$.

$L_b$ has a closed form in logarithms, $\pi$, Euler's constant and a finite sum over prime powers below $b$.

## Key numbers

| Quantity | Value |
|---|---|
| $B_2-L_2$ | about $0.2199792174$, certified between $0.2198$ and $0.2201$ |
| $B_2+L_2$ | about $0.7505835310$, certified between $0.7504$ and $0.7508$ |
| A first sufficient width for eventual bias | $b>1.84529052\ldots$ |
| Small-scale obstruction | for $1<b<\sqrt2$ the prime 2 alone exceeds the main term at $y=2/b$ |

At $b=2$, bramble_trace had already obtained the upper extremum by simultaneous recurrence. Their deficit at centre $2y$ is one third of ours. Converting the square-root normalisation multiplies their constant by $3\sqrt2$, taking $0.17691423488\ldots$ to $0.7505835310\ldots$. The formula above treats every $b>1$. [bramble_trace, revision 24](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24)

## How the proof goes

1. On the critical line $H_b(\rho)$ is $e^{i\tau\log b}$ times a positive number $h_\rho$. This has a separate Lean proof.
2. The explicit formula writes the deficit as a sum over zeta zeros plus prime-power terms.
3. Prime squares contribute $B_b\sqrt y$ in the limit, and higher powers are of smaller order.
4. This gives $D_b(y)/\sqrt y=B_b+\sum_\rho h_\rho e^{i\tau\log(by)}+o(1)$, an absolutely convergent series. Aligning the phases gives the limsup, and Kronecker's theorem gives the liminf under the linear independence hypothesis in the statement.
5. A contour integral evaluates $L_b=\sum_\rho h_\rho$ in closed form.

Every step is written out in [PROOF.md](PROOF.md).

## Verification

The executable certifies the dyadic numerical enclosures. The explicit formula, prime number theorem and simultaneous phase approximation are used in the written proof.

```sh
bash verify-release.sh
```

This checks the [SHA-256 manifest](SHA256SUMS), reruns the arithmetic and compares its output with the [recorded result](evidence/check.txt). It needs Python 3.9 or later and only the standard library. The recorded run used Python 3.9.6 on arm64 macOS. [Latest release](https://github.com/leecrossley/prime-band-extrema/releases/latest)

## Related results

- [Formal proof of the strip and phase identities](https://github.com/leecrossley/prime-band-mellin/tree/v1.0.1)
- [Eventual positivity at 151/100](https://github.com/leecrossley/prime-band-151-100/tree/v1.0.1)
- [Conditional oscillation at three halves](https://github.com/leecrossley/prime-band-three-halves-oscillation/tree/v1.0.1)
