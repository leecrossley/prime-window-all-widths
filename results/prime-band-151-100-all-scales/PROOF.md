# Proof: positivity at every scale for dilation 151/100

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

## Two lower bounds

Set $b=151/100$, $r=\sqrt b$, $Q(y)=\sum_p(\log p)W_b(p^2/y)$ and $B=2r(r-1)^2$. For $\rho=1/2+i\tau$ define

$$
h_\tau=\frac{r(b+1-2r\cos(\tau\log b))}{1/4+\tau^2},\qquad
L=\sum_\rho h_\tau,
$$

counting multiplicity and both signs. The explicit formula for the prime-power sum gives, for $y>1$ and any finite collection $Z$ of distinct positive ordinates,

$$
\frac{D_b(y)}{\sqrt y}
\ \ge\ \frac{Q(y)}{\sqrt y}-L+
   \sum_{\tau\in Z}2h_\tau(1+\cos(\tau\log(by))).
$$

Indeed each conjugate pair contributes $2h_\tau\cos(\tau\log(by))$. Bound all other cosines below by $-1$. Omitted prime powers of exponent at least three and the trivial-zero terms have favourable sign. The Mellin transform is positive at each negative even integer by its positive integral. The series over the zeros converges absolutely. [Explicit formula and the phase of the moving tent](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)

Schoenfeld gives $|\vartheta(t)-t|\le \sqrt t\thinspace\log^2t/(8\pi)$ for $t\ge599$ under RH. For $m\ge599$ let $\epsilon_m=\log^2m/(8\pi\sqrt m)$. For $\sqrt y\ge m$ partial summation gives

$$
Q(y)/\sqrt y \ \ge\ B-2r(b-1)\epsilon_m.
$$

The integral of $W_b(v^2)$ from 1 to $b$ equals $B$. On the two pieces the error in replacing the prime integral by the ordinary integral is at most $\epsilon_m\int_1^b v|dW_b(v^2)|=2r(b-1)\epsilon_m$. The weight vanishes at both exterior endpoints, and $\log^2t/\sqrt t$ decreases for $t>e^4$. [Schoenfeld, Theorem 10, inequality (6.3)](https://doi.org/10.1090/S0025-5718-1976-0457374-X)

## Coverage of all real scales

| Scale | Certificate and resulting sign |
|---|---|
| $0<y\le10^5$ | Sum over ordinary primes directly. The deficit exceeds 0.059 on all 43,786 positive breakpoints and the upper endpoint taken together. Between consecutive points it is affine. At zero it extends to zero. |
| $10^5\le y\le10^{10}$ | On 34,848 intervals between prime-square breakpoints, $Q(y)>0.110\sqrt y$. The finite zero sum above exceeds 0.020 and $L<0.126$. |
| $10^{10}\le y\le(21{,}000{,}000)^2$ | Taking $m=10^5$ gives $Q(y)>0.107\sqrt y$. The same phase bound exceeds 0.020, which still covers the zero loss. |
| $y\ge(21{,}000{,}000)^2$ | Taking $m=21{,}000{,}000$ gives $B-L-2r(b-1)\epsilon_m>0.00003$. The empty collection of ordinates suffices in this range. |

For the finite prime checks use the integer coordinate $u=151^2y$. A summand at $n=p$ or $n=p^2$ changes formula only at $u=100^2n$, $151\cdot100n$ and $151^2n$. Prefix sums of $\log p$ and $(\log p)/n$ evaluate each exact piece. On a square interval, the affine function $(u+d^2)/(2\cdot151d)$ bounds $\sqrt y$ above for every integer $d>0$. Testing both endpoints against this tangent proves the square inequality throughout.

For the middle two ranges, $11\le\log(by)\le35$. The main certificate isolates 40 positive zeta zeros using FLINT. It evaluates the gain on 24,000 closed intervals covering all of $[11,35]$ with rigorous ball arithmetic. The second certificate uses 24 independently isolated zeros, direct signs of the real completed zeta function, a derivative bound and rational integer bounds for the prime sums. [FLINT zeta zeros](https://python-flint.readthedocs.io/en/latest/acb.html#flint.acb.zeta_zeros), [Independent check](independent.py)
