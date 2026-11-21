# Proof: the transition in a family of weighted prime sums

This is the full argument behind the [summary](README.md).

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

## Derivation

Under RH, put $B_b=2\sqrt b(\sqrt b-1)^2$ and let $L_b$ be the sum over nontrivial zeta zeros, with multiplicities, of the positive moduli of the Mellin coefficients. The explicit formula yields an absolutely convergent trigonometric series for $D_b(y)/\sqrt y$ up to a term tending to zero. It gives limsup $B_b+L_b$, liminf at least $B_b-L_b$, and equality at the lower endpoint under the independence hypothesis. The sign classification therefore follows once the sign of $B_b-L_b$ is known. [Series, drift, zero mass and phase proof](https://github.com/leecrossley/prime-window-all-widths/tree/main/results/prime-band-extrema)

For $1<b\le2$ the von Mangoldt ramp at $b$ vanishes because its first possible prime-power term is at 2 and has zero weight there. Substituting this into the zero-mass formula gives

$$
\begin{aligned}
F(b)&=(B_b-L_b)/\sqrt b\\
&=2(\sqrt b-1)^2-(b+1)\log(b+1)+(b-1)\log(b-1)\\
&\quad -2(b-1)\log b+2b\log2+(b-1)(\gamma+\log\pi).
\end{aligned}
$$

Here $\gamma$ is Euler's constant. The continuous extension has $F(1)=0$. Differentiation gives

$$
\begin{aligned}
 F'(b)&=2/b-2/\sqrt b+\log\frac{4\pi e^\gamma(b-1)}{b^2(b+1)},\\
 b^2(b^2-1)F''(b)&=P(\sqrt b),\qquad
 P(r)=-2r^6+r^5+2r^2-r+2.
\end{aligned}
$$

For $r\ge1$, $P'(r)=-12r^5+5r^4+4r-1\le-3r^4-1<0$. Since $P(1)=2$ and $P(\sqrt2)=3\sqrt2-10<0$, $F'$ first increases and then decreases. Also $F'(b)\to-\infty$ as $b\downarrow1$, whereas $F'(2)=1-\sqrt2+\gamma+\log(\pi/3)>0$. It follows that $F'$ crosses zero exactly once. Thus $F$ decreases from 0 to a strict minimum and then strictly increases. The two certified opposite signs in the displayed bracket locate its unique interior zero.

For $b\ge2$ the coefficient bound $L_b\le C_0\sqrt b(\sqrt b+1)^2$ holds, where $C_0=2+\gamma-\log(4\pi)$. Hence it suffices that $2((\sqrt b-1)/(\sqrt b+1))^2>C_0$. The left side increases for $b>1$ and the strict inequality already holds at 2. This proves the assertion on the full half-line.
