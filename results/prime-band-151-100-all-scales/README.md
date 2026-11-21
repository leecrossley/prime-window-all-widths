# Positivity at every scale for dilation 151/100

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410846.svg)](https://doi.org/10.5281/zenodo.24410846)

At dilation $b=151/100$, the weighted sum of primes between $y$ and $2.2801y$ is strictly below its predicted size at every real scale $y>0$, from the first prime onwards, if the Riemann hypothesis (RH) holds. No independence assumption on zeta zeros is needed.

An earlier result gave this only for large enough $y$, without saying how large. Here four ranges, each with its own certificate, together cover every $y$.

## Statement

Assume RH. Then $D_{151/100}(y)>0$ for every real $y>0$. The proof also gives $D_{151/100}(y)>0.00003\sqrt y$ for $y\ge441{,}000{,}000{,}000{,}000$.

For $b>1$ and $y>0$,

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

[Earlier eventual bound](https://github.com/leecrossley/prime-band-151-100/tree/v1.0.1)

## Coverage of all real scales

Under RH the explicit formula bounds the deficit from below by the prime-square sum, minus the most the zeta zeros can remove, plus whatever a chosen finite set of zeros provably gives back. In units of $\sqrt y$ these are written $Q(y)/\sqrt y$, $L$ and the phase gain.

| Scale | How positivity is certified |
|---|---|
| $0<y\le10^5$ | Direct sum over primes. The deficit exceeds 0.059 at each of the 43,786 breakpoints, the upper endpoint included, and is affine between them. At zero it extends to zero. |
| $10^5\le y\le10^{10}$ | $Q(y)>0.110\sqrt y$ on 34,848 intervals between prime-square breakpoints. The loss is $L<0.126$ and the phase gain from 40 zeros exceeds 0.020. |
| $10^{10}\le y\le(21{,}000{,}000)^2$ | Schoenfeld's RH bound gives $Q(y)>0.107\sqrt y$. The same phase gain still covers the loss. |
| $y\ge(21{,}000{,}000)^2$ | Schoenfeld's RH bound for $Q(y)/\sqrt y$ exceeds the full loss $L$ by more than 0.00003, with every zero at its worst phase and no phase gain. |

Both lower bounds and the details of each certificate are in [PROOF.md](PROOF.md).

## Two certificates

Both programs check the direct prime sums for $0<y\le10^5$, the prime-square intervals on $[10^5,10^{10}]$ and the phase gain on $[11,35]$. Only `verify.py` checks $L<0.126$, $Q(y)>0.107\sqrt y$ for $y\ge10^{10}$, the tail margin above $0.00003$ and the joins $0.110-L+0.020>0$ and $0.107-L+0.020>0$.

- `verify.py` isolates 40 positive zeta zeros using FLINT and evaluates the gain on 24,000 closed intervals covering all of $[11,35]$ with rigorous ball arithmetic.
- `independent.py` uses 24 independently isolated zeros, verified by sign changes of the real completed zeta function around published trial centres, together with a derivative bound and rational integer bounds for the prime sums. It requires no external data file.

The explicit formula and partial summation are proved or referenced in the written proof.

With Python 3.10 or newer, run the commands in a virtual environment.

```sh
python -m pip install -r requirements.txt
python verify.py
python independent.py
```

[Verification output](evidence/check.txt), [file checksums](SHA256SUMS), [FLINT zeta zeros](https://python-flint.readthedocs.io/en/latest/acb.html#flint.acb.zeta_zeros)

## Prior work and scope

The dyadic instance, including an all-scale sign argument, appears in bramble_trace's same-day fixed-tent work. With their centre $x=2y$, our dyadic deficit equals three times theirs. Grenié and Molteni used the same centered exponential autocorrelation in class-group bounds. For $y\ge2$ this family is also a negative three-point difference of the integrated first Mertens error in Zhao's paper. These precedents supply the observable and much of the method. [Dyadic proof and certificate, revision 24](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24), [Grenié and Molteni, proof of Theorem 3.5](https://arxiv.org/pdf/1607.02430#page=6), [Zhao, Theorem 1 and Section 2](https://arxiv.org/html/2411.18903v2)

The contribution here is the all-scale conclusion at $151/100$, including the explicit coverage joining its finite and asymptotic ranges. I found no exact earlier result at this dilation in the sources compared.
