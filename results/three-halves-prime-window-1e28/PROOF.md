# Proof: the three-halves prime window through $10^{28}$

This is the full argument behind the [summary](README.md).

## Setup

For every real $0<y\le10^{28}$,

$$
\sum_{p\ {\rm prime}}(\log p)\thinspace W(p/y)<\tfrac34 y\log(3/2).
$$

Put $b=3/2$. The weight $W(t)$ equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere. Write $D(y)$ for the main term minus this ordinary-prime sum, and $Q(y)=\sum_p(\log p)W(p^2/y)$ for the square contribution. All logarithms are natural.

## Separating known zeros from the unknown ones

For $y>1$, the explicit formula for the continuous piecewise smooth weight gives

$$
D(y)=\sum_{\rho}y^\rho H_b(\rho)+\sum_{m\ge1}y^{-2m}H_b(-2m)
       +\sum_{k\ge2}\sum_p(\log p)W(p^k/y),\qquad
H_b(s)=\frac{(b^s-1)(b^s-b)}{s(s-1)} .
$$

The trivial-zero terms are positive and every prime-power weight is nonnegative. The nontrivial-zero series converges absolutely, including off the line, since $|H_b(\beta+i\tau)|\le 2b(b+1)/\tau^2$ for $0<\beta<1$. The formula follows by Mellin inversion of the von Mangoldt sum, or by smoothing the weight inside its support and taking the limit. [Explicit formula](https://people.maths.ox.ac.uk/greenbj/papers/primenumbers.pdf), [Mellin and phase identity](https://github.com/leecrossley/prime-window-all-widths/tree/main/results/prime-band-mellin)

Set $H=3\cdot10^{12}$. Every zero through this height is on the critical line by Platt and Trudgian. For each conjugate pair below $H$ the contribution to $D/\sqrt y$ is $a_\tau\cos(\tau\log(by))$, where

$$
a_\tau=\frac{2\sqrt b\thinspace(b+1-2\sqrt b\cos(\tau\log b))}{1/4+\tau^2}.
$$

Let $M$ be an upper bound for $\sum_{0<\tau\le H}a_\tau$. At a cutoff $T=1000$, sum the 649 coefficients below $T$ and bound the others by $2\sqrt b(\sqrt b+1)^2 V(T)$. Here

$$
V(T)=\frac{1+\log(T/(2\pi))}{2\pi T}+
       \frac{0.14+0.56\log T}{T^2}
$$

bounds $\sum_{\tau>T}\tau^{-2}$ counting **all** positive ordinates, with multiplicity. At the cutoff it follows from Brent, Platt and Trudgian, Lemma 1 and Example 1. Their half-weight convention at an endpoint only enlarges the upper bound for the strict tail. This gives $M<0.13121$. [Platt and Trudgian, Theorem 1](https://arxiv.org/abs/2004.09765), [Brent, Platt and Trudgian, Example 1](https://arxiv.org/pdf/2009.13791#page=9)

## The error above the verified height

For $\sigma>1/2$, set $S_\sigma=\sum_{\beta>\sigma,\ \tau>H}\tau^{-2}$ over zeros in the upper half-plane. There are no such zeros up to $H$, so partial summation gives $S_\sigma=2\int_H^\infty N(\sigma,t)t^{-3}dt$. Kadiri, Lumley and Ng give

$$
N(\sigma,t)\le A_\sigma t^p(\log t)^q+C_\sigma(\log t)^2,\qquad
p=\tfrac83(1-\sigma),\ q=5-2\sigma ,
$$

at the 16 cuts $0.80,0.85,0.86,\ldots,0.99$. The certificates round their Table 1 constants upward. With $h=\log H$ and $t=He^v$, use $(1+v/h)^q\le e^{qv/h}$. Thus

$$
U_\sigma=\min\left\lbrace\frac{V(H)}2,
\ \frac{2A_\sigma H^{p-2}h^q}{2-p-q/h}+
\frac{2C_\sigma H^{-2}h^2}{2-2/h}\right\rbrace\ \ge S_\sigma .
$$

The first bound uses reflection $\beta+i\tau\leftrightarrow1-\beta+i\tau$ at the same height. At most half of the mass lies strictly right of $1/2$. [Kadiri, Lumley and Ng, equation (1.7) and Table 1](https://www.cs.uleth.ca/~kadiri/articles/Explicit-zero-density-for-the-Riemann-zeta-function-JMAA-May2018.pdf)

Write the cuts as $s_1,\ldots,s_{16}$, append $s_{17}=1$, and set $f_j=y^{s_j-1/2}$. By charging the baseline 1 to every positive ordinate and the excess in each strip to those to its right, the absolute contribution above $H$ to $D(y)/\sqrt y$ is at most

$$
E(y)=4b(b+1)\left[V(H)+(f_1-1)\frac{V(H)}2+
             \sum_{j=1}^{16}(f_{j+1}-f_j)U_{s_j}\right].
$$

The factor includes both conjugate signs. Every increment increases for $y\ge1$. The main certificate gives $E(10^{28})<0.0505$ and uses smaller enclosures at each interval's right endpoint.

Now any finite set $Z$ of retained pairs below $H$ gives

$$
D(y)/\sqrt y\ \ge\ Q(y)/\sqrt y-M-E(y)+
   \sum_{\tau\in Z}a_\tau[1+\cos(\tau\log(by))].
$$

## The prime side and the continuous range

For $y\le10^5$, $D$ is concave between adjacent centres $p/b$. Indeed $-W(p/y)$ only jumps upward in slope at $p/b$ and its other two slope jumps are downward. Checking the 13,848 centres and the final endpoint gives strict positivity. At zero the extension equals zero, which covers the first open interval by concavity.

For $10^5\le y\le10^{10}$ the prime-square sum is affine between consecutive squared-prime breakpoints $y=p^2/b^2,p^2/b,p^2$. An affine upper tangent to $\sqrt y$ on each of 34,728 intervals gives $Q(y)>0.11\sqrt y$ there. For larger $y$, write $x=\sqrt y$, $B=2\sqrt b(\sqrt b-1)^2$ and $V_2=2\sqrt b(b-1)$. Integration by parts against $W(v^2)$ gives

$$
Q(y)/\sqrt y\ \ge B-\frac{V_2\log^2x}{8\pi\sqrt x}.
$$

Büthe's Theorem 2 applies to every $\vartheta(t)$ needed here. Indeed $10^5\le t\le b\thinspace10^{14}$ and $4.92\sqrt{t/\log t}<H$. The relative error decreases for $x\ge10^5$. These estimates give $Q(y)/\sqrt y>0.10$ on $[10^5,10^{28}]$ and $>0.123$ from $10^{20}$ onward. [Büthe, Theorem 2](https://arxiv.org/pdf/1410.7015#page=13)

The main program encloses the first 128 phase terms on 138,000 closed logarithmic intervals. On each interval their sum, with the square lower bound, exceeds $M+E$.
