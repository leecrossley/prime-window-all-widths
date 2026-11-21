# The sharp width for a prime deficit at every scale

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410881.svg)](https://doi.org/10.5281/zenodo.24410881)

Assume the Riemann hypothesis (RH). Once the width parameter reaches $b_\ast=1.50342897722316\ldots$, the weighted sum of ordinary primes in $(y,b^2y)$ stays strictly below its main term for **every real scale $y>0$ and every real $b\ge b_\ast$**. This includes the critical width itself.

Below the cutoff, assume also that the distinct positive ordinates of the zeta zeros are linearly independent over the rationals. For each fixed $1<b<b_\ast$ the sum then passes above and below its main term arbitrarily far out. The threshold and this oscillation result were established in the [earlier transition theorem](https://github.com/leecrossley/prime-band-transition/tree/5fa6cfd68d6afa1451e45a56f6a88195c8b74dff).

## Which primes count

Every rational prime $p$ contributes $\log p$ times the weight below, with $t=p/y$. The weight vanishes at the two ends and reaches $b-1$ at $p=by$.

$$
W_b(t)=\begin{cases}
b(1-1/t),&1 < t\le b,\\
b^2/t-1,&b < t < b^2,\\
0,&\text{otherwise}.
\end{cases}
$$

The deficit is $D_b(y)=b(b-1)y\log b-\sum_p(\log p)W_b(p/y)$. For each fixed $b\ge b_\ast$, strict positivity at every scale is equivalent to RH. The converse uses the classical one-sign criterion for Mellin transforms. [Proof, including the converse](PROOF.md)

## How the missing widths are covered

Earlier work proved the all-scale inequality from $b_\ast$ through $1.51$. A positive sum of dilates expresses the deficit at $b^n$ in terms of deficits at $b$. This reduces the rest of the problem to $1.51\le b\le 2.27$. [Earlier critical-interval proof](https://github.com/leecrossley/prime-band-critical-all-scales/tree/ecb4538f96452cd25fff0d10fab52f71117d57ea)

| Scale on this width interval | Argument |
|---|---|
| $0<y\le 10^5$ | Direct prime counts and concavity. Two integer enclosures verify 39,689,083 positive endpoints across 2,375 width intervals. This part is unconditional. |
| $10^5\le y\le 10^{10}$ | Prime squares and cubes together exceed the maximum loss from all the zeros under RH. Two checks bound the sum between its breakpoints. |
| $y\ge 10^{10}$ | All even prime powers form a second von Mangoldt sum. Under RH its zero contribution costs less than $0.73y^{1/4}$. The deficit exceeds $0.00068\sqrt y$ throughout this range. |

## Prior work

Grenié and Molteni used the same centered kernel. The observable is a three-point difference of Zhao's integrated first Mertens error. Unser and Blu give the integer refinement filter for exponential splines. At width 2, bramble_trace had already proved positivity at every scale. The contribution here is the continuous all-width, all-scale theorem and the certified estimates that close the intervening range. The sources compared are listed with their formulas in the proof. [Grenié–Molteni](https://arxiv.org/pdf/1607.02430#page=6) · [Zhao](https://arxiv.org/html/2411.18903v2) · [Unser–Blu, §IV-D](https://bigwww.epfl.ch/publications/unser0503.pdf) · [bramble_trace, revision 24](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24)

## Reproduce the verification

With Python 3.10 or later in a virtual environment, run

```sh
python3 -m pip install -r requirements.txt
python3 verify-release.py
```

The command checks the file hashes and runs both certificates. The two programs implement the same written argument separately, so they share its analytic inputs and the width partition in `widths.json`. The explicit formula and the passage from the bounded width interval to every larger width are written proofs. [Recorded output](evidence/check.txt) · [SHA-256 manifest](SHA256SUMS) · [Release](https://github.com/leecrossley/prime-window-all-widths/releases/tag/v1.0.1)
