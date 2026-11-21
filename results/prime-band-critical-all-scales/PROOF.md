# Proof: the critical prime window stays positive at every scale

This is the full argument behind the [summary](README.md).

## Statement

Let $b_\ast$ be the unique dilation in $(1.50342897722316,1.50342897722317)$ where the square drift equals the full zero mass. Assuming RH, for every $b_\ast\le b\le151/100$ and **every real $y>0$** the ordinary-prime deficit $D_b(y)$ is strictly positive. For $\log y\ge96$ the stronger uniform estimate is

$$D_b(y)\ \ge\ (B_b-L_b)\sqrt y+0.035\thinspace y^{1/3}.$$

In particular it stays positive at $b=b_\ast$ itself. If the distinct positive zeta-zero ordinates are also linearly independent over the rationals, then $\liminf D_{b_\ast}(y)/\sqrt y=0$ nonetheless. [Definition, uniqueness and enclosure of the threshold](https://github.com/leecrossley/prime-band-transition/tree/v1.0.1)

The liminf assertion is part of the cited threshold theorem. The argument below proves the positivity assertions.

## Setup

For $b>1$, $W_b(t)$ is $b(1-1/t)$ on $1<t\le b$, $b^2/t-1$ on $b<t<b^2$, and zero elsewhere. Write $D_b(y)=b(b-1)y\log b-\sum_p(\log p)W_b(p/y)$ with $p$ running over ordinary primes. Under RH let $L_b$ denote the sum of the moduli of the nontrivial zero coefficients, including both signs and multiplicities, and $B_b=2\sqrt b(\sqrt b-1)^2$. The threshold theorem gives $B_b\ge L_b$ throughout the stated interval. [Zero expansion and mass](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

## The cubes at equality

For $k=2,3$ put $Q_k(y)=\sum_p(\log p)W_b(p^k/y)$, $r_k=b^{1/k}$ and

$$C_k=\frac{k}{k-1}(r_k-1)(b-r_k),\qquad V_k=\frac{k}{k-1}(b-r_k+br_k-r_k^2).$$

Integration over $1\le v\le r_k^2$ gives $C_k=\int W_b(v^k)\thinspace dv$ and $V_k=\int v\thinspace|dW_b(v^k)|$. These may be checked by splitting at $r_k$. In particular $C_2=B_b$ and $C_3>0$. Schoenfeld's RH bound and integration by parts give, for $y^{1/k}\ge599$,

$$
Q_k(y)\ge C_k y^{1/k}
-\frac{V_k}{8\pi} y^{1/(2k)}(\log y/k)^2 .
$$

Here $\log^2t/\sqrt t$ decreases on the full integration interval. The weight vanishes at both endpoints. [Schoenfeld, Theorem 10, (6.3)](https://doi.org/10.1090/S0025-5718-1976-0457374-X)

In the explicit formula the trivial-zero terms in the deficit are positive for $y>1$. Retain squares and cubes, bound the nontrivial zero sum below by $-L_b\sqrt y$, and drop the other nonnegative prime-power contributions. At $\log y=t\ge96$ the net coefficient of $y^{1/3}$ is bounded below by

$$
C_3-\frac{V_2}{8\pi}(t/2)^2e^{-t/12}
    -\frac{V_3}{8\pi}(t/3)^2e^{-t/6}.
$$

Both losses decrease for $t\ge96$. Uniformly for $1.50342\le b\le1.51$, the certificate encloses $C_3>0.077$, $V_2<1.26$ and $V_3<1.18$. Substitution at 96 leaves more than 0.035. This proves the infinite range even when $B_b-L_b=0$.

## Every smaller scale

For each positive ordinate $\tau$ write $a_\tau=2\sqrt b(b+1-2\sqrt b\cos(\tau\log b))/(1/4+\tau^2)$. For a finite set $Z$ of distinct positive ordinates,

$$
D_b(y)/\sqrt y\ge Q_2(y)/\sqrt y-L_b+
      \sum_{\tau\in Z}a_\tau(1+\cos(\tau\log(by)))\qquad(y>1).
$$

This retains one conjugate pair for each selected ordinate. All unselected terms have already been allowed their worst sign in $-L_b$. [Finite phase argument at 151/100](https://github.com/leecrossley/prime-band-151-100-all-scales/tree/v1.0.1)

| Range of $y$ | Proof |
|---|---|
| $0<y\le10^5$ | Cover $[1.50,1.51]$ with 128 rational dilation intervals $[l,u]$. Both the main coefficient and $W_b(t)$ increase in $b$ for each fixed $t>1$. Use the main term at $l$ and the prime sum at $u$, and check the resulting lower bound at every centre $p/u$ and the upper endpoint. It exceeds 0.055 at all 1,778,200 tested positions, counting each slice separately. At zero it extends to zero. |
| $10^5\le y\le10^{10}$ | $Q_2(y)$ dominates its value at $b=3/2$. A tangent certificate on 34,728 intervals proves $Q_2(y)/\sqrt y>0.11$. In the target range $L_b\le B_b<0.13$. The retained phases contribute more than 0.022. |
| $10^{10}\le y\le e^{96}$ | The relative Schoenfeld loss in $Q_2/\sqrt y$ is less than 0.021, using $\sqrt y\ge10^5$. The same uniform phase contribution exceeds 0.022. |
| $y\ge e^{96}$ | The cube estimate above applies. |

For the first row the lower comparison is concave between successive prime centres. Each function $-W_u(p/y)$ has a downward slope jump at $p/u^2$ and at $p$, and its only upward jump is at $p/u$. An affine main term preserves concavity. The gaps between tested centres, including the first open interval, therefore follow from the endpoint bounds.

For the middle rows $11\le\log(by)\le97$. A rigorous lower bound on each of 64 coefficients, uniform over the whole dilation interval, reduces the calculation to one variable. The primary certificate encloses the resulting phase sum over 86,000 closed intervals. A second certificate uses the completed zeta function to prove 64 separate zero crossings, bounds the change of the coefficients with dilation and bounds the slope of the trigonometric sum between time midpoints.
