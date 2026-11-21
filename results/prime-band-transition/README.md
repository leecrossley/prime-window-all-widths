# The transition in a family of weighted prime sums

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410825.svg)](https://doi.org/10.5281/zenodo.24410825)

Under the hypotheses stated below, weighted sums of primes over the window from $y$ to $b^2y$ behave in two different ways depending on the dilation $b$. This repository finds the single value where the behaviour switches, $b_\ast\approx1.50342897722316$, and encloses it rigorously to 14 decimal places.

Above $b_\ast$, the Riemann hypothesis (RH) implies that the primes fall short of the prediction at all sufficiently large $y$. Below it, RH together with rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros implies that both signs keep occurring. The threshold marks where the one-sided bias of the wider windows begins.

## Statement

There is a unique $b_\ast\in(1,2)$ at which the difference between the prime-square drift and the total absolute zero contribution vanishes. It satisfies

$$1.50342897722316 < b_\ast < 1.50342897722317.$$

| Dilation | Behaviour of $D_b(y)$ |
|---|---|
| $b>b_\ast$ | Under RH, $D_b(y)>0$ for all sufficiently large $y$. |
| $1<b<b_\ast$ | Under RH and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros, both signs occur arbitrarily far out. |
| $b=b_\ast$ | Under RH, $D_b(y)>0$ for every real $y>0$. With the independence hypothesis the normalised liminf is zero. [Critical-width proof](https://github.com/leecrossley/prime-band-critical-all-scales/blob/99f0a80a6f332a43942344c470ab410eba0eb2f7/README.md) |

For $b>1$ and $y>0$, let

$$
W_b(t)=
\begin{cases}
 b(1-1/t),&1 < t\le b,\\
 b^2/t-1,&b < t < b^2,\\
 0,&\text{otherwise},
\end{cases}
\qquad
D_b(y)=b(b-1)y\log b-\sum_{p\ {\rm prime}}\log p\thinspace W_b(p/y).
$$

The prime-square drift is $B_b=2\sqrt b(\sqrt b-1)^2$ and the total absolute zero contribution is $L_b=\sum_\rho|H_b(\rho)|$ over the nontrivial zeta zeros with multiplicity, where $H_b(s)=(b^s-1)(b^s-b)/(s(s-1))$ is the Mellin transform of $W_b$.

## How the proof goes

1. The sign of $B_b-L_b$ gives the asymptotic bounds in the table. Strict positivity at equality uses the later critical-width proof linked in that row.
2. For $1<b\le2$ the quotient $F(b)=(B_b-L_b)/\sqrt b$ has a closed form in logarithms, Euler's constant and $\pi$.
3. Calculus shows that $F$ decreases from 0 to a strict minimum and then strictly increases, so it has at most one zero in $(1,2)$.
4. Certified opposite signs of $F$ at the two ends of the displayed bracket show that this zero exists and lies in the bracket.
5. For $b\ge2$ a coefficient bound gives $B_b-L_b>0$ directly.

The derivation is in [PROOF.md](PROOF.md). The series, drift and phase argument behind step 1 are in [prime-band-extrema](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1).

## Prior work and scope

The dyadic instance, including an all-scale sign argument, appears in bramble_trace's same-day fixed-tent work. With their centre $x=2y$, our dyadic deficit equals three times theirs. Grenié and Molteni used the same centered exponential autocorrelation in class-group bounds. For $y\ge2$ this family is also a negative three-point difference of the integrated first Mertens error in Zhao's paper. These precedents supply the observable and much of the method. [Dyadic proof and certificate, revision 24](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24), [Grenié and Molteni, proof of Theorem 3.5](https://arxiv.org/pdf/1607.02430#page=6), [Zhao, Theorem 1 and Section 2](https://arxiv.org/html/2411.18903v2)

The contribution here is the unique transition over the full continuous dilation family and its certified enclosure. The rational-independence assumption applies only to the assertion of lower-endpoint equality and recurring negative values. Priority beyond the sources compared remains open.

## Reproduce the bounds

With Python 3.10 or newer, run the commands in a virtual environment.

```sh
python -m pip install -r requirements.txt
python verify.py
python verify_ramp.py
```

The first program uses Arb ball arithmetic on the displayed closed form. The second uses standard-library rational intervals, Machin's formula, the signed Euler-Maclaurin remainder for Euler's constant and 78 terms of the unshifted trivial-zero series with a geometric remainder. The calculus proof is in [PROOF.md](PROOF.md). [Signed remainder for the digamma expansion](https://dlmf.nist.gov/5.11#ii), [Verification output](evidence/check.txt), [File checksums](SHA256SUMS)
