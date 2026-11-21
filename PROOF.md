# Proof of the all-width prime-window theorem

Use the natural logarithm, count each rational prime once, and set

$$
M_b=b(b-1)\log b,\qquad
D_b(y)=M_by-\sum_p(\log p)W_b(p/y).
$$

Here $b>1$, $y>0$ and $W_b(t)$ equals $b(1-1/t)$ on $1<t\le b$, $b^2/t-1$ on $b<t<b^2$, and zero elsewhere. Both exterior endpoints have weight zero. We prove under RH that $D_b(y)>0$ simultaneously for all $b\ge b_\ast$ and $y>0$. The root $b_\ast$ is defined below.

## The zero mass

With Mellin measure $t^{s-1}dt$ the transform of the weight is

$$
H_b(s)=\frac{(b^s-1)(b^s-b)}{s(s-1)} .
$$

Integration over the two branches gives this identity and the removable values at 0 and 1. In particular $H_b(1)=M_b$. For $y>1$ the explicit formula for the sum over all prime powers is

$$
\sum_n\Lambda(n)W_b(n/y)=M_by-\sum_\rho y^\rho H_b(\rho)
                                 -\sum_{j\ge1} y^{-2j}H_b(-2j).
$$

All zero sums in this proof count multiplicity and both signs of the ordinate. The kernel is continuous and piecewise smooth in log coordinates. Mellin inversion against $-\zeta'/\zeta$ gives the displayed poles. Equivalently one may smooth in log coordinates and pass to the limit. On RH the nontrivial-zero terms have a summable quadratic denominator, and the trivial-zero tail is geometric for $y>1$. At a negative real argument the transform is positive by its defining integral. [Explicit-formula details for this kernel](https://github.com/leecrossley/prime-window-all-widths/tree/main/results/prime-band-extrema) · [Ben Green, *Prime number theory*](https://people.maths.ox.ac.uk/greenbj/papers/primenumbers.pdf)

Write $\rho=1/2+i\tau$, $r=\sqrt b$, and define

$$
B_b=2r(r-1)^2,\qquad
L_b=\sum_\rho |H_b(\rho)|
    =\sum_\rho\frac{r(b+1-2r\cos(\tau\log b))}{1/4+\tau^2}.
$$

The phase of $H_b(\rho)$ is $e^{i\tau\log b}$. In the specified range $1<b\le 2.27$ let $R(b)=\max(b/2-1,0)\log 2$, where $x_+=\max(x,0)$. The single-scale von Mangoldt ramp gives

$$
\frac{L_b}{r}=(b+1)\log(b+1)-(b-1)\log(b-1)+2(b-1)\log b
       -2b\log2-(b-1)(\gamma+\log\pi)-2R(b).
$$

For the alternative calculation let $C_0=2+\gamma-\log(4\pi)=\sum_\rho 1/[\rho(1-\rho)]$. The completed zeta logarithmic derivative gives this constant. Define $S(b)=\sum_\rho b^\rho/[\rho(1-\rho)]$. Shifting the integral of $(-\zeta'/\zeta)(s)b^s/[s(s-1)]$ gives

$$
S(b)=R(b)+b(1+\gamma-\log b)-\log(2\pi)
             +\sum_{j\ge1}\frac{b^{-2j}}{2j(2j+1)},\qquad
L_b=r[(b+1)C_0-2S(b)] .
$$

Here 2 is the only prime power where the ramp can change slope below 2.27. Summing the elementary trivial-zero series gives the preceding closed form. The earlier transition theorem proves the unique zero $b_\ast\in(1,2)$ of $B_b-L_b$, encloses it in $(1.50342897722316,1.50342897722317)$ and proves $B_b-L_b<0$ below it. Under RH and rational independence of the **distinct** positive ordinates, finite simultaneous phase approximation then gives both signs arbitrarily late at each fixed $1<b<b_\ast$. At equality the normalized liminf is zero under those hypotheses. The earlier critical-interval theorem proves strict positivity at every scale for $b_\ast\le b\le1.51$, including equality. We use these results at their cited revisions. [Transition and oscillation](https://github.com/leecrossley/prime-window-all-widths/tree/main/results/prime-band-transition) · [Critical interval](https://github.com/leecrossley/prime-window-all-widths/tree/main/results/prime-band-critical-all-scales)

## Small scales

For fixed $t>0$, the branches of $W_b(t)$ meet continuously and increase with $b$ wherever they have positive weight. $M_b$ increases for $b>1$. For a closed width cell $[\ell,u]$ the function

$$
M_\ell y-\sum_p(\log p)W_u(p/y)
$$

therefore bounds $D_b(y)$ from below. As a function of $y$ it is concave between consecutive centres $p/u$. Indeed, the two outer contacts $p/u^2,p$ give downward slope jumps to the negative of a tent, and only its centre gives an upward jump. Check all positive centres through $10^5$ and the final endpoint. The lower comparison extends to zero at $y=0$, so its chord to the first positive checked endpoint is strictly positive for every intervening $y>0$.

The included partition covers $[1.51,2.27]$ with 2,375 adjoining closed rational intervals. The two verifiers check the partition itself and 39,689,083 centres and endpoints each. The first uses directed bounds on the two weight branches. The second computes three signed cumulative ramps with a separately sieved prime list. The small calculation needs no assumption on zeros. [Partition](widths.json) · [Certificate output](evidence/check.txt)

## Squares and cubes through ten billion

Set $Q_{k,b}(y)=\sum_p(\log p)W_b(p^k/y)$. In passing from the full von Mangoldt sum to ordinary primes, each higher power contributes positively to the deficit. The explicit formula therefore gives

$$
D_b(y)\ge Q_{2,b}(y)+Q_{3,b}(y)-L_b\sqrt y\qquad(y>1).
$$

For each row below, the two verifiers prove $L_b<h$ on the entire closed width cell and $Q_{2,\ell}(y)+Q_{3,\ell}(y)>h\sqrt y$ on the entire closed interval $10^5\le y\le10^{10}$. Monotonicity of the two prime-power weights in $b$ completes the comparison.

| Left width $\ell$ | Right width | $h$ | Affine intervals | Lower chords |
|---|---|---|---:|---:|
| 1.51 | 1.54 | 0.127 | 35,901 | 24,074 |
| 1.54 | 1.65 | 0.140 | 36,262 | 24,333 |
| 1.65 | 1.85 | 0.195 | 37,550 | 25,247 |
| 1.85 | 2.27 | 0.300 | 39,866 | 26,916 |

For $\ell=a/c$ and $n=p^2,p^3$, the first verifier splits at every $n/\ell^2,n/\ell,n$ in range. On each resulting interval the sum $Q$ is affine. Write $z=a^2y$ and choose the positive integer $m=\lfloor\sqrt{z_{\rm left}}\rfloor$. The affine function $(z+m^2)/(2am)$ majorizes $\sqrt y$ by the arithmetic-geometric mean inequality. Testing $Q$ against $h$ times this tangent at both endpoints covers the interval.

The second verifier splits only at $n/\ell^2,n$. Between them $Q$ is concave, since its remaining slope changes at tent centres are downward. If $U,V$ are positive lower bounds for $Q$ at endpoints $u,v$ in the original scale, its linear chord is a lower bound throughout. The checks are $U^2>h^2u$, $V^2>h^2v$, and $2UV>h^2(u+v)$. To see sufficiency, subtract $h^2((1-t)u+tv)$ from $((1-t)U+tV)^2$ and expand in $(1-t)^2, 2t(1-t), t^2$. All three coefficients are positive. Thus the two counts in the table correspond to different continuous enclosures.

The first spectral check encloses the closed expression for $L_b$ on subdivisions of each width cell. The second uses the unshifted ramp $S(b)$, sums 24 terms of its trivial-zero series and bounds the remainder by $x^{25}/[50\cdot51(1-x)]$, $x=b^{-2}$. It encloses $\pi$ using Machin's formula and $\gamma$ using the signed harmonic-number expansion at 64. Each verifier also obtains $B_b-L_b>0.003$ throughout $[1.51,2.27]$. [Remainder in the digamma expansion](https://dlmf.nist.gov/5.11#ii) · [Two implementations](evidence/check.txt)

## Large scales from all even powers

Write $x=\sqrt y$. The terms with even prime-power exponent form the single sum

$$
E_b(y)=\sum_n\Lambda(n) W_b((n/x)^2).
$$

Each $n=p^j$ contributes $\log p$ at $p^{2j}$ exactly once. The Mellin transform of $u\mapsto W_b(u^2)$ is $K_b(s)=H_b(s/2)/2$. Its pole main term in the explicit formula is $xK_b(1)=B_b\sqrt y$. On RH,

$$
|K_b(\rho)|=\frac{2|(b^{\rho/2}-1)(b^{\rho/2}-b)|}{|\rho|\thinspace|\rho-2|}
\le \frac{2(b^{1/4}+1)(b^{1/4}+b)}{|\rho|^2}.
$$

Since $|\rho-2|\ge|\rho|$ and $\sum_\rho|\rho|^{-2}=C_0$, the sum of these upper bounds is less than $0.73$ at $b=2.27$ and at every smaller width above 1. The nontrivial zeros thus cost at most $0.73 y^{1/4}$. The trivial zeros in $E_b$ have the opposite sign to the ones in the deficit. At $s=-j$ the explicit transform satisfies $0<H_b(-j)\le b/[j(j+1)]$. Their entire cost for $y>1$ is at most $b/[4(y-1)]<1/(y-1)$. Dropping the other nonnegative contributions from the deficit now gives

$$
D_b(y)\ge (B_b-L_b)\sqrt y-0.73y^{1/4}-\frac{1}{y-1}.
$$

The certified gap $0.003$ and substitution at $y=10^{10}$ leave more than $0.00068\sqrt y$. Both normalized losses decrease thereafter.

## Passing to every width

Let $T f(y)=f(by)$ and $R_p(y)=\sum_p(\log p)\max(y/p-1,0)$. A single prime contributes its tent to $(T-1)(T-b)R_p$, including all boundary cases. Acting on $y\log y$ gives $M_b y$. Thus $D_b=(T-1)(T-b)(y\log y-R_p(y))$. Factor $T^n-1$ and $T^n-b^n$ to get

$$
D_{b^n}(y)=\sum_{j,k=0}^{n-1} b^{n-1-k}D_b(b^{j+k}y),\qquad n\ge1 .
$$

Every coefficient is positive. For any $B\ge b_\ast$ take $n=\lfloor\log B/\log b_\ast\rfloor$. Then $1\le n\le\log B/\log b_\ast<n+1\le2n$, so $b=B^{1/n}$ belongs to $[b_\ast,b_\ast^2]$. The root enclosure gives $b_\ast^2<2.27$. The critical-interval theorem and the three scale ranges above cover this base interval. The identity therefore proves strict positivity at the original width $B$ for every $y>0$.

## Converse and comparisons

For any fixed $b>1$, even eventual nonnegativity implies RH. Put $U(s)=\sum_{p,k\ge2}(\log p)p^{-ks}$, holomorphic in $\Re s>1/2$ (the sum runs over primes $p$ and integers $k\ge2$). For $\Re s>1$, the transform of $D_b$ on $[1,\infty)$ is

$$
\frac{M_b}{s-1}-H_b(s)\left(-\frac{\zeta'(s)}{\zeta(s)}-U(s)\right)
     +\text{an entire function}.
$$

The entire correction comes from prime-window integrals over $0<y<1$, where only finitely many primes contribute and their supports stay away from 0. The pole at 1 cancels. There are no real zeta zeros in $(1/2,1)$ by the alternating eta series. Every zero with $1/2<\Re\rho<1$ gives a genuine pole, since $1<|b^\rho|<b$. For an eventually nonnegative function Landau's theorem forces a real singularity at its finite abscissa of convergence. Holomorphy at every real point above $1/2$ therefore excludes those nonreal poles. Functional-equation symmetry proves RH. This known direction supplies the equivalence stated in the README. [Kaczorowski–Pintz, introduction and Theorem 1](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1986.Kaczorowski_1.pdf) · [Earlier dyadic formalization](https://github.com/leecrossley/prime-window-rh-criterion/tree/f1e75967261918dec01883332a91ef072348d210)

Zhao's integrated first Mertens error $J$ differs from $R_p(y)-y\log y$ by an affine function for $y\ge2$. Thus $D_b(y)=-J(b^2y)+(b+1)J(by)-bJ(y)$ there. The sign of a cumulative error alone leaves its three-point difference undetermined. The centered expression $e^{u/2}W_b(be^u)=b e^{-|u|/2}-e^{|u|/2}$ for $|u|\le\log b$ is the exponential autocorrelation in Grenié and Molteni's Theorem 3.5. The two geometric sums in the dilation identity have the standard integer exponential-spline refinement form, with parameters $0,\log b$ up to reversal and normalization. [Zhao, Theorem 1 and §2](https://arxiv.org/html/2411.18903v2) · [Grenié–Molteni](https://arxiv.org/pdf/1607.02430#page=6) · [Unser–Blu, §IV-D, equations (39)–(40)](https://bigwww.epfl.ch/publications/unser0503.pdf)

bramble_trace established the dyadic sign on 20 November 2026. With their centre $x=2y$ the earlier normalized deficit is $D_2(y)/3$. Their later notes also contain a finite-height calculation for the fixed dyadic window. The earlier family papers establish eventual positivity above the transition and all-scale positivity through width $1.51$. The estimates above close the missing continuum in both variables. I found no earlier proof of that simultaneous assertion in the works compared. The analytic argument uses the ordinary zeta function. Alpöge's theta lifts motivated the investigation, and the finite identity for supplied adjacent rows is a separate result. [Dyadic proof, revision 24](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=24) · [Updated source, revision 35](https://rstation.ai/api/v1/posts/uw5cfb_J_Cj7XJO6/attachment?revision=35) · [Alpöge's announcement](https://alpo.ge/remark.pdf) · [Finite trace identity](https://github.com/leecrossley/theta-host-split-inert-race/blob/main/results/theta-prime-trace/ThetaPrimeTrace.lean#L66)
