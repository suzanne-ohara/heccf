# hecqf-code

Code accompanying the paper : [O26]

> Suzanne O'Hara, "Heavenly elliptic curves over cubic fields,"
> to appear.

## Overview

This repository contains a script for a calculation that supports results appearing in [026].

* `balanced-filter.sage`: The balanced bound B(n,g) is shown to satisfy $B(3,1) < 1.80 \times 10^{16}$ in Corollary 5.3 of [MR26], but may be reduced to $B(3,1) \leq 73$ by a computational verification. The existence of a heavenly elliptic curve $E/K$ which is not balanced implies the Tate-Oort numbers $(j_{1}, j_{2})$ satisfy $j_{1} + j_{2} = e$ and the congruence (2.3):  
$$x := \tau_{e} - q^{j_{1}} - q^{j_{2}} \equiv 0 \pmod{\ell}$$.  
Here, $\mathfrak{p}$ is a prime of $K$ dividing $p$, a rational prime, $q = \mathbf{N}\mathfrak{p} = p^{f}$, and $\tau_{e}$ is the trace of Frobenius for $\mathfrak{p}^{e}$. The crucial observation is that (2.3) must hold for every $\mathfrak{p} \nmid \ell$. Taking $p = 3$ and factoring $x$ for every possible choice of $f$ and $\tau_{e}$, we obtain a finite set of possible *unbalanced* $\ell$. The script then checks the congruence (2.3) against primes $2 \leq p \leq 11$, removing any primes $\ell$ that fail to satisfy (2.3) for all pairs $(p,f)$ *with the same $(j_{1}, j_{2})$ for each $p$*.


## Usage

The balanced filter runs in a few seconds:

```bash
sage balanced-filter.sage
```


## Authors

* Suzanne O'Hara (Wesleyan University)

Comments welcome: seohara at wesleyan dot edu

## References



[BCP97] W. Bosma, J. Cannon, and C. Playoust, "The Magma algebra system. I.
The user language," *Journal of Symbolic Computation* **24** (1997),
pp. 235--265. Available at [http://magma.maths.usyd.edu.au/](http://magma.maths.usyd.edu.au/)

[Deu57] M. Deuring, "Die Zetafunktion einer algebraischen Kurve vom Geschlechte
Eins (iv)," *Nachrichten der Akademie der Wissenschaften in Göttingen*,
1957, pp. 55--80.

[DL15] H. Daniels and Á. Lozano-Robledo, "On the number of isomorphism classes of CM elliptic curves defined over a number field," *J. Number Theory* **157** (2015), pp. 367--396.

[K02] M. Kida, "Potential good reduction of elliptic curves," *J. Symbolic Computation* (2002) **34**, pp. 173--180. [doi:10.1006/jsco.2002.0555](https://doi.org/10.1006/jsco.2002.0555)

[Lan73] S. Lang, *Elliptic Functions*, Addison-Wesley, 1973.

[MR26] C. McLeman and C. Rasmussen, "Heavenly elliptic curves over quadratic
fields," *Res. number theory* (2026) **12**. [doi:10.1007/s40993-026-00757-8](https://doi.org/10.1007/s40993-026-00757-8)

[O26] S. O'Hara, "Heavenly elliptic curves over cubic fields," to appear, 2026.

[Sage] The Sage Developers, *SageMath, the Sage Mathematics Software System*.
Available at [https://www.sagemath.org/](https://www.sagemath.org/)
