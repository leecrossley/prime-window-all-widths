# Proof: how rare is a reversal in the three-halves prime window?

This is the full argument behind the [summary](README.md).

## Setup

Fix $b=3/2$. Give each ordinary prime $p$ weight $(\log p)W_b(p/y)$, where

$$
W_b(t)=\begin{cases}
 b(1-1/t)&1 < t\le b,\\
 b^2/t-1&b < t < b^2,\\
 0&\text{otherwise}.
\end{cases}
$$

Write $D_b(y)=b(b-1)y\log b-\sum_p(\log p)W_b(p/y)$. A reversal means $D_b(y)<0$. Its frequency up to $X>1$ is

$$
d(X)=\frac{1}{\log X}\int_1^X {\bf1}_{D_b(y)<0}\thinspace\frac{dy}{y}.
$$

All logarithms are natural.

## The finite test

Under RH the weighted explicit formula and the ordinary-prime square term give

$$
\frac{D_b(y)}{\sqrt y}=-\Delta+
  \sum_{\tau>0}a_\tau[1+\cos(\tau\log(by))]+o(1),
\qquad
a_\tau=\frac{2\sqrt b\thinspace[b+1-2\sqrt b\cos(\tau\log b)]}{1/4+\tau^2}.
$$

The sum counts conjugate pairs with multiplicity. The coefficients are positive and summable. Write $B=2\sqrt b(\sqrt b-1)^2$, $L=\sum_{\tau>0}a_\tau$ and $\Delta=L-B$. At $b=3/2$ we have $0<\Delta<0.00163$. A closed expression for this mass at $1<b\le2$ is

$$
L=\sqrt b\left[(b+1)\log(b+1)-(b-1)\log(b-1)+2(b-1)\log b
-2b\log 2-(b-1)(\gamma+\log\pi)\right],
$$

where $\gamma$ is Euler's constant. This absolute convergence lets us give every omitted coefficient its most adverse phase. [Derivation from the explicit formula](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

Let $\tau_1,\ldots,\tau_{500}$ be the first 500 positive ordinates and put $S(t)=\sum_{j=1}^{500}a_{\tau_j}(1+\cos(\tau_j t))$. For every $\varepsilon>0$, a reversal at sufficiently large $y$ gives $S(\log(by))\le\Delta+\varepsilon$. For $\lambda>0$ this bounds its indicator above by $e^{\lambda(\Delta+\varepsilon)}e^{-\lambda S(\log(by))}$. The average for the variable $t=\log(by)$ is the usual average over an interval of length $\log X$.

Take $d=500000$, $\lambda=150000$, $v_j=\lambda a_{\tau_j}$, $c_j=(1-e^{-2v_j/d})/2$ and $u_j=d c_j$. On $0\le q\le2$, convexity bounds $e^{-v_j q/d}$ by the chord $1-c_j q$. It follows that the nonnegative trigonometric polynomial

$$
P_j(\theta)=(1-c_j(1+\cos\theta))^d
$$

majorizes $e^{-v_j(1+\cos\theta)}$. It has Fourier degree at most $d$. Best and Trudgian established that no nonzero relation $\sum n_j\tau_j=0$ exists among these 500 ordinates for integer coefficients of absolute value at most $d$. In their Table 2 the bound for the first 500 is $5.54632\times10^5$. Consequently the time average of $\prod_j P_j(\tau_j t)$ is exactly the product of the constant Fourier coefficients. This step uses their published finite calculation. [Best and Trudgian, Appendix A, Table 2](https://arxiv.org/pdf/1209.3843#page=13)

Since $1-c_jq\le e^{-c_jq}$, each constant coefficient is at most $e^{-u_j}I_0(u_j)$, where $I_0(u)$ is the circle average of $e^{u\cos\theta}$. Taking the upper limit, then letting $\varepsilon\downarrow0$, gives

$$
\limsup_{X\to\infty}d(X)\ \le\
 e^{\lambda\Delta}\prod_{j=1}^{500} e^{-u_j}I_0(u_j)
 \ <\ 10^{-232}.
$$

This bound concerns logarithmic frequency. It gives no location for the first reversal.
