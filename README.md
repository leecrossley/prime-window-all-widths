# Two-sided density bounds for reversals of the three-halves prime window

Assume the Riemann hypothesis and rational linear independence of the distinct positive ordinates of the nontrivial zeta zeros. The $3/2$ weighted prime window overshoots its main term on a set with logarithmic density $\delta$ satisfying

$$10^{-10000}<\delta<10^{-1012}.$$

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

## Why the density exists

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

where $\gamma$ is Euler's constant. This absolute convergence lets us give every omitted coefficient its most adverse phase. [Derivation from the explicit formula](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.0)

Under the independence hypothesis, Kronecker's theorem gives uniform distribution for any finite vector of distinct phases modulo $2\pi$ as $\log(by)$ ranges over an interval of increasing length. Group repeated ordinates, assigning their full multiplicity to the coefficient. The positive coefficients are summable, hence finite sums approximate the full sum uniformly. Its limiting distribution is that of $-\Delta+\sum a_\tau(1+\cos\Theta_\tau)$, with independent uniform angles for the distinct ordinates and the indicated multiplicity in each coefficient. This distribution has no atom at zero. Conditioning on the other angles leaves one nonconstant cosine with a continuous distribution. Approximation of the sign indicator by continuous functions now gives $\delta=\lim_{X\to\infty}d(X)$, even with the $o(1)$ term in the explicit formula. [Rubinstein and Sarnak's prime-race framework](https://www.math.uwaterloo.ca/~mrubinst/publications/Chebyshev.pdf)

## Upper bound

The first 5,000 positive ordinates are distinct. On the reversal event their contribution $\sum_{j\le5000}a_j(1+\cos\Theta_j)$ is at most $\Delta$, since all later terms are nonnegative. Exponential Markov bounding with $\lambda=1170000$ gives

$$
\delta\le e^{\lambda\Delta}\prod_{j=1}^{5000}
                 e^{-\lambda a_j}I_0(\lambda a_j)<10^{-1012}.
$$

Here $I_0(u)$ is the circle average of $e^{u\cos\theta}$. The inequalities use certified enclosures for the ordinates, weights and Bessel functions.

## Lower bound

Take $N=15000$, $e=0.0004$, $A_N=\sum_{j\le N}a_j$, and impose the following condition on just the first $N$ angles,

$$
|\Theta_j-\pi|\le\phi_j=\min\left(\pi,\sqrt{2e/(Na_j)}\right)\quad(1\le j\le N),
$$

with angles in $[0,2\pi]$. Since $1-\cos v\le v^2/2$, these $N$ modes together contribute at most $e$. All remaining modes contribute at most $2(L-A_N)$, whatever their angles. The certificate proves $e+2(L-A_N)<\Delta$, so this whole box has negative normalized deficit. Its volume under independent uniform angles is

$$
\prod_{j=1}^N \frac{\phi_j}{\pi}>10^{-10000}.
$$

In particular the reversals occupy a positive proportion on the logarithmic scale. These bounds leave the location of the first reversal open.

## Prior work and checks

The existence argument belongs to the classical theory of prime races. Lamzouri studied large deviations for prime number races, and Zhao studied oscillation and bias for integrals of Mertens errors. The earlier result for this window establishes recurring conditional reversals. These two quantitative bounds apply to the specific $3/2$ ordinary-prime window. The literature checked yielded no earlier numerical interval for its reversal density. [Lamzouri](https://arxiv.org/abs/1103.0060) · [Zhao](https://arxiv.org/html/2411.18903v2) · [Earlier oscillation statement](https://github.com/leecrossley/prime-band-three-halves-oscillation/tree/v1.0.0)

The first program isolates 15,000 consecutive positive ordinates with FLINT and computes both bounds in logarithms. The second verifies opposite signs of the completed zeta function on both sides of every supplied centre, then confirms the total count. It uses the unshifted ramp with a bounded series remainder, evaluates the Bessel upper bounds by positive series, and multiplies the lower box volume directly. Both use ball arithmetic.

```sh
python -m pip install -r requirements.txt
python verify.py
python independent.py
```

[Recorded runs](evidence/check.txt) · [Checksums](SHA256SUMS)
