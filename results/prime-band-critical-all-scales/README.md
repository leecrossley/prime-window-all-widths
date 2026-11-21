# The critical prime window stays positive at every scale

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410842.svg)](https://doi.org/10.5281/zenodo.24410842)

Weighted sums of primes over a window from $y$ to $b^2y$ change behaviour at one exact dilation, $b_\ast\approx1.50342897722316$. Under the Riemann hypothesis (RH), for each fixed $b>b_\ast$ the sum stays below its predicted size for all sufficiently large $y$. Under RH and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros, for each fixed $1<b<b_\ast$ it is both above and below its predicted size at arbitrarily large $y$. Assuming RH, this repository shows that the weighted prime sum stays strictly below its predicted size at every real scale $y>0$, for every dilation from $b_\ast$ up to $151/100$, including $b_\ast$ itself.

The threshold is the hard case. There the steady push from prime squares exactly cancels the worst the zeta zeros can do, so the usual square-root argument leaves no margin. For $\log y\ge96$ the much smaller contribution of prime cubes, of size $y^{1/3}$, settles it, and finite certificates cover the smaller scales. Under RH and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros, the normalised lower limit at $b_\ast$ is nonetheless zero.

## Statement

Let $b_\ast$ be the unique dilation in $(1.50342897722316,1.50342897722317)$ where the square drift equals the full zero mass. Assuming RH, for every $b_\ast\le b\le151/100$ and **every real $y>0$** the ordinary-prime deficit $D_b(y)$ is strictly positive. For $\log y\ge96$ the stronger uniform estimate is

$$D_b(y)\ \ge\ (B_b-L_b)\sqrt y+0.035\thinspace y^{1/3}.$$

In particular it stays positive at $b=b_\ast$ itself. If the distinct positive zeta-zero ordinates are also linearly independent over the rationals, then $\liminf D_{b_\ast}(y)/\sqrt y=0$ nonetheless. [Definition, uniqueness and enclosure of the threshold](https://github.com/leecrossley/prime-band-transition/tree/v1.0.1)

For $b>1$, $W_b(t)$ is $b(1-1/t)$ on $1<t\le b$, $b^2/t-1$ on $b<t<b^2$, and zero elsewhere. Write $D_b(y)=b(b-1)y\log b-\sum_p(\log p)W_b(p/y)$ with $p$ running over ordinary primes. Under RH let $L_b$ denote the sum of the moduli of the nontrivial zero coefficients, including both signs and multiplicities, and $B_b=2\sqrt b(\sqrt b-1)^2$. The threshold theorem gives $B_b\ge L_b$ throughout the stated interval. [Zero expansion and mass](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

## How the proof goes

1. **Very large scales.** Under RH, for $\log y\ge96$, the explicit formula and Schoenfeld's RH bound applied to the prime squares and prime cubes give $D_b(y)\ge(B_b-L_b)\sqrt y+0.035\thinspace y^{1/3}$, uniformly for $1.50342\le b\le1.51$. On $b_\ast\le b\le151/100$ the threshold theorem gives $B_b\ge L_b$, which makes the right side positive.
2. **Everything below.** Three more ranges are each certified uniformly for $b_\ast\le b\le151/100$.

| Range of $y$ | Why the deficit is positive |
|---|---|
| $0<y\le10^5$ | Direct check of the prime sum over 128 dilation slices, at 1,778,200 positions in all, with margin above 0.055. The lower bound is concave between positions. |
| $10^5\le y\le10^{10}$ | Prime squares give $Q_2(y)/\sqrt y>0.11$ on 34,728 intervals, against $L_b<0.13$. Retained terms from 64 zeros add more than 0.022. |
| $10^{10}\le y\le e^{96}$ | Schoenfeld's bound loses less than 0.021 from the squares, and the same retained zeros add more than 0.022. Since $L_b\le B_b$, this leaves a positive margin. |

The full argument, including the constants $C_k$ and $V_k$ and the slicing in $b$, is in [PROOF.md](PROOF.md).

## Verification

With Python 3.10 or newer, run the following in a virtual environment. The zero centres used by the independent program were computed with FLINT and are checked again from sign changes. The certificates verify the finite inequalities in the written proof. The two programs implement the same written argument separately, so they share its analytic inputs, including the threshold theorem.

```sh
python -m pip install -r requirements.txt
python verify.py
python independent.py
```

[Recorded checks](evidence/check.txt), [checksums](SHA256SUMS)

## Related work

The centered exponential weight goes back to Grenié and Molteni. Zhao treats the underlying integrated Mertens error. The fixed dyadic case has a public all-scale certificate by bramble_trace, and the $151/100$ endpoint was previously checked separately. The assertion here covers the continuous interval down to its exact transition. The new positive margin at the boundary comes from prime cubes after the square drift and total zero mass cancel. [Grenié and Molteni](https://arxiv.org/pdf/1607.02430#page=6), [Zhao](https://arxiv.org/html/2411.18903v2), [Dyadic certificate](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24)
