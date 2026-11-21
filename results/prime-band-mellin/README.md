# The phase of the Mellin transform of a compact prime weight

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24410836.svg)](https://doi.org/10.5281/zenodo.24410836)

The tent weights in this family have a simple transform, $H_b(s)=(b^s-1)(b^s-b)/(s(s-1))$. This repository gives a Lean proof of two facts about it: it never vanishes for $0<\Re s<1$, and on the critical line it equals $b^s$ times a positive real number.

Of the companion repositories, prime-band-extrema cites this one, and it uses only the critical-line phase. Because the transform never vanishes in the strip, the weight cannot cancel the contribution of a zeta zero. In prime-band-extrema the exact phase lets the zero terms line up to give the limsup. The closed form there for the total zero contribution comes from the positive real factor together with a contour integral.

## Statement

For every real $b>1$:

- $H_b(s)$ is nonzero for $0<\Re s<1$. In Lean this is `ThetaTransfer.bandMellin_ne_zero_in_strip`.
- On the critical line $H_b(s)$ is $b^s$ times a strictly positive real factor. In Lean this is `ThetaTransfer.bandMellin_phase`.

Both are in [BandSpectrum.lean](BandSpectrum.lean). Written out for $\rho=1/2+i\tau$, the identity $b/b^\rho=\overline{b^\rho}$ gives

$$
 H_b(\rho)=e^{i\tau\log b}h_\rho,\qquad
 h_\rho=\frac{\sqrt b\thinspace(b+1-2\sqrt b\cos(\tau\log b))}{1/4+\tau^2}>0.
$$

Here $H_b$ is the Mellin transform of the weight $W_b(t)$, which equals $b(1-1/t)$ for $1<t\le b$, $b^2/t-1$ for $b<t<b^2$, and zero elsewhere.

## Proof

Integrating the two pieces against $t^{s-1}dt$ yields the displayed rational expression in $b^s$. On the open strip the modulus of $b^s$ lies strictly between $1$ and $b$, so neither numerator factor vanishes. For $\Re s=1/2$, put $z=b^s$, so $z\bar z=b$ and $s(s-1)=-|s|^2$. Consequently $(z-1)(z-b)=-z(b+1-2\Re z)$. The real multiplier of $z$ is positive since $b+1-2\Re z\ge(\sqrt b-1)^2>0$.

## Verification

The two assertions about the displayed expression are kernel-checked against Mathlib. The identification with the Mellin transform follows by integrating the two elementary pieces. The Lean statements are made away from the removable singularities at 0 and 1.

For a fresh build on x86-64 Linux with Git, curl, tar, zstd and Python 3.9 or later, run

```sh
bash verify-linux.sh
```

This checks the [SHA-256 manifest](SHA256SUMS), downloads the pinned dependencies and compares the resulting kernel and axiom summary with the [recorded run](evidence/check.txt).

This result imports Mathlib directly. With the pinned toolchain and Mathlib cache already installed, supply an empty directory as the source argument:

```sh
cd /path/to/pinned/mathlib4
lake env python3 /path/to/repository/verify.py /path/to/empty/directory /path/to/build
```

The build directory can begin empty. Reusing a dependency cache requires the same sources and compiler.

The recorded run used Lean 4.33.1 on arm64 macOS. It freshly elaborated this repository's module, used compiled Mathlib from the pinned revision, replayed the module with `leanchecker`, and inspected the named declarations' axioms. Compiler commit `819816b2e0a3bf405af45ae5c7af2491d8f5bee6` and Mathlib commit `0df444a360eaa60ab8c11dca51a86af692955474` fix the environment. [Lean release](https://github.com/leanprover/lean4/releases/tag/v4.33.1), [Mathlib revision](https://github.com/leanprover-community/mathlib4/tree/0df444a360eaa60ab8c11dca51a86af692955474), [Latest release](https://github.com/leecrossley/prime-band-mellin/releases/latest)

## Related result

[Derivation of the general extremal formula](https://github.com/leecrossley/prime-band-extrema/tree/v1.0.1)
