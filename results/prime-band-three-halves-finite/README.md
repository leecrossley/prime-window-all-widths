# The three-halves prime deficit through 100,000

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410830.svg)](https://doi.org/10.5281/zenodo.24410830)

At dilation $b=3/2$, the weighted sum of primes between $y$ and $2.25y$ is strictly below its predicted size for every real $y$ up to 100,000. This is proved outright by exact arithmetic, with no hypothesis.

It is a statement about a finite range only. A companion result shows that, assuming the Riemann hypothesis and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros, the sign changes arbitrarily far out, so positivity on this range does not persist for all $y$ under those hypotheses.

## Statement

Set $b=3/2$. For every real $0<y\le100{,}000$, $D_b(y)>0$. This is an unconditional finite-range result.

## Setup

For $b>1$ and $y>0$:

- **Weight.** $W_b(t)$ equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere.
- **Prime sum.** $P_b(y)=\sum_{p\ {\rm prime}}(\log p)W_b(p/y)$, with natural logarithms. Its expected main term is $b(b-1)y\log b$.
- **Deficit.** $D_b(y)=b(b-1)y\log b-P_b(y)$.

## Proof

Write $b=a/c$ with $a=3$ and $c=2$, and put $u=a^2y=9y$.

1. For a prime $p$ the continuous weight $W_b(p/y)$ changes slope only at $c^2p$, $acp$ and $a^2p$ in $u$. After multiplying the weight by $ac^2$, its slope changes are respectively $a/p$, $-(a+c)/p$ and $c/p$.
2. An integer sieve through $225{,}000$ includes every contributing prime.
3. Coefficients $(\log p)/p$ are rounded upwards with rational log enclosures, before applying the nonnegative tent.
4. The main slope uses a lower enclosure of $(a-c)\log(a/c)$ with the same fixed-point scale.
5. Prefix summation evaluates these affine bounds at every integer $1\le u\le900{,}000$.

Each endpoint upper bound for the prime sum is strictly below the endpoint lower bound for the main term. Both bounds are affine between successive integers. At $u=0$ both vanish, so the strict inequality at $u=1$ covers the first open interval as well. The full assertion follows by interpolation.

## Verification

The executable checks all 900,000 positive endpoints, using integer and rational arithmetic. This covers every real scale in the stated interval.

```sh
bash verify-release.sh
```

This checks the [SHA-256 manifest](SHA256SUMS), reruns the arithmetic and compares its output with the [recorded result](evidence/check.txt). It needs Python 3.9 or later and only the standard library. The recorded run used Python 3.9.6 on arm64 macOS. [Latest release](https://github.com/leecrossley/prime-band-three-halves-finite/releases/latest)

## Related result

[Conditional oscillation at the same dilation](https://github.com/leecrossley/prime-band-three-halves-oscillation/tree/v1.0.1)
