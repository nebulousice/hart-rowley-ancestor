# Hart–Rowley Conjecture 1.3 in type A, and a type-B consequence

Author: Sukhraj Sachdeva. Date: 2026-11-27. This is the full proof behind the paper `paper.tex` / `paper.pdf` in this directory.

**Status.** Theorem A and Theorem B are **proved**: every step is argued by hand in this file, and the argument passed two rounds of review by AI referees (see the last section). Computations are **numerical** and are corroboration only: no step of the proof uses a computation. **Expert verification required.**

**Scope and credit.** This file makes no claim of priority or novelty. The conjecture and all definitions are due to S. B. Hart and P. J. Rowley ([arXiv:2502.00777](https://arxiv.org/abs/2502.00777); Arch. Math. 126 (2026) 231–237, [doi:10.1007/s00013-025-02220-z](https://doi.org/10.1007/s00013-025-02220-z)). The uniqueness of the ancestor (their Conj. 1.2), the fact that the ancestor is Bruhat-above every involution prefix, and the complementation step used in Lemma 2.6 are due to C. Gaetz and Y. Gao ([arXiv:2610.22857v1](https://arxiv.org/abs/2610.22857v1), unrefereed). Their proof of these facts (their Lemma 3.1, Theorem 1.1 and Corollary 1.2) is written out again, in our notation and with every step justified from textbook facts, in Lemma 2.6 and Appendix A; the argument is theirs. Conjecture 1.3 is **false in general**: S. Levinstein gave counterexamples in $H_3$, $F_4$, $E_6$, $E_7$ ([doi:10.5281/zenodo.24581942](https://doi.org/10.5281/zenodo.24581942)). What is shown here is that the conjecture does hold in type A.

## 0. Results

**Theorem A.** Let $N \ge 1$ and $w \in S_N$, with the Coxeter generators $s_1,\dots,s_{N-1}$ (adjacent transpositions). Then

```math
\ell_{\mathrm{inv}}(w) \le \mathrm{disp}(w) \le N-1 = \mathrm{rank}(S_N), \qquad \mathrm{disp}(w) = \max_{1\le i\le N} \lvert w(i)-i \rvert .
```

So Hart–Rowley Conjecture 1.3 holds for the Coxeter system of type $A_{N-1}$, for every $N$. Both inequalities are equalities for $w = s_1 s_2\cdots s_{N-1}$ (Proposition 4.6).

**Theorem B.** Let $n\ge 1$ and let $B_n$ be the group of signed permutations, with its Coxeter generators $s_0,\dots,s_{n-1}$, viewed inside $S_{2n}$ as in §5.1. For every $w\in B_n$ the involution length of $w$ in $B_n$ equals the involution length of $w$ in $S_{2n}$, and therefore

```math
\ell_{\mathrm{inv}}^{B_n}(w) = \ell_{\mathrm{inv}}^{S_{2n}}(w) \le \mathrm{disp}_{S_{2n}}(w) \le 2n-1 .
```

**Theorem B is not Conjecture 1.3 for type B.** The rank of $B_n$ is $n$, and the conjecture for $B_n$ asks for $\ell_{\mathrm{inv}}^{B_n}\le n$. The bound $2n-1$ is weaker for every $n\ge 2$. Conjecture 1.3 for $B_n$ stays **open** here. The same reduction is false for type D (Remark 5.6).

**Dependence.** Theorems A and B are stated without hypotheses. Their proofs use the textbook facts about Coxeter groups listed in §1.4 (E3) and one further statement, (E1) = Gaetz–Gao Cor. 1.2 (§1.4). (E1) is not assumed: it is proved in Appendix A from the same textbook facts, following Gaetz–Gao's proof. Nothing in Theorems A and B rests on the unrefereed preprint as an authority; it is credited as the source of the argument of Lemma 2.6 and Appendix A.

**Outline of the proof of Theorem A.** Let $a(w)$ be the ancestor and $\mathrm{rem}(w)=a(w)^{-1}w$. It is enough to show $\mathrm{disp}(\mathrm{rem}(w)) < \mathrm{disp}(w)$ for $w\ne e$. The remainder is a suffix of $w$, so it is Bruhat-below $w$, and displacement cannot grow in Bruhat order (Lemma 4.1). To see that it drops, take a position $p$ where $w$ has its largest displacement $D$ as an excedance, $w(p)=p+D$. Lemma K (Lemma 4.2) builds an explicit involution prefix $\iota$ of $w$ that moves the value $p+D$. By (E1) $\iota$ is Bruhat-below $a(w)$; by the complementation step (Lemma 2.6) $\mathrm{rem}(w)$ is Bruhat-below $\iota w$; and the condition "position $p$ holds the value $p+D$" is an up-set condition among permutations of displacement at most $D$. Since $\iota w$ does not satisfy it, neither does $\mathrm{rem}(w)$. Largest drops are reduced to largest excedances by the diagram automorphism $x\mapsto w_0xw_0$.

## 1. Definitions and conventions

### 1.1 The conjecture

Source: Hart–Rowley, [arXiv:2502.00777v1](https://arxiv.org/abs/2502.00777v1), §1; the journal version (Arch. Math. 126 (2026) 231–237) states Conjecture 1.3 in the same words. The definitions are restated in §1.2 in the notation of this file. The conjecture reads, verbatim: "Suppose Conjecture 1.2 holds, and let $W$ be a finite Coxeter group. For each $w \in W$, the involution length of $w$ is at most the rank of $W$." (Hart–Rowley, Conj. 1.3.)

Conjecture 1.2 (every non-identity element of a finite Coxeter group has a unique longest involution prefix) is now a theorem of Gaetz–Gao for all Coxeter groups, (E1) below, and is proved again in Appendix A. So the hypothesis of Conjecture 1.3 is satisfied and the involution length is defined for every element.

### 1.2 Coxeter groups

Let $(W,S)$ be a Coxeter system, $e$ its identity, $\ell$ its length function, and $T=\lbrace wsw^{-1} : w\in W, s\in S\rbrace$ its set of reflections. The rank is $\lvert S\rvert$. We use the basic properties of length ([BB, Ch. 1]): $\ell(x^{-1})=\ell(x)$, $\ell(xy)\le \ell(x)+\ell(y)$, and $\ell(xs)=\ell(x)\pm 1$, $\ell(sx)=\ell(x)\pm1$ for $s\in S$.

- **Prefix; right weak order.** $u\le_R w$ means $\ell(w)=\ell(u)+\ell(u^{-1}w)$. This is the source's notion of prefix: if $r_1\cdots r_k$ and $r_1\cdots r_n$ are reduced words for $u$ and $w$, then $r_{k+1}\cdots r_n$ is a word for $u^{-1}w$, so $\ell(u^{-1}w)\le n-k$, and $\ell(u^{-1}w)\ge \ell(w)-\ell(u)=n-k$ by subadditivity; conversely, if $\ell(w)=\ell(u)+\ell(u^{-1}w)$, then a reduced word for $u$ followed by a reduced word for $u^{-1}w$ is a word for $w$ of length $\ell(w)$, hence reduced. The source names this order the weak (left) Bruhat order; only the name differs.
- **Suffix.** $r$ is a suffix of $w$ if $\ell(w)=\ell(wr^{-1})+\ell(r)$.
- **Bruhat order.** $x\le_B y$ means there is a chain $x=x_0,x_1,\dots,x_m=y$ ($m\ge0$) with $x_{i}=x_{i-1}t_i$, $t_i\in T$ and $\ell(x_i)>\ell(x_{i-1})$ (the Bruhat order of [BB, Ch. 2]). Along a chain the length increases strictly, so $x\le_B y$ and $x\ne y$ force $\ell(x)<\ell(y)$.
- **Involution.** An element of order exactly 2. The identity is not an involution.
- **Ancestor.** For $w\ne e$, $\mathcal{IP}(w)$ is the set of involutions $u$ with $u\le_R w$. It is non-empty: the first letter of a reduced word of $w$ is in it. The ancestor $a(w)$ is the element of maximal length in $\mathcal{IP}(w)$; it is unique by (E1).
- **Remainder and involution length.** $\mathrm{rem}(w)=a(w)^{-1}w=a(w)w$ for $w\ne e$. Then $\ell(\mathrm{rem}(w))=\ell(w)-\ell(a(w))<\ell(w)$, and we set

```math
\ell_{\mathrm{inv}}(e)=0,\qquad \ell_{\mathrm{inv}}(w)=1+\ell_{\mathrm{inv}}(\mathrm{rem}(w))\quad (w\ne e).
```

This is the number $k$ of factors of the source's ancestor decomposition $w=a_1\cdots a_k$, $a_i=a(a_i\cdots a_k)$: if $w$ is an involution then $w\in\mathcal{IP}(w)$ is its own longest prefix, so $a(w)=w$, $\mathrm{rem}(w)=e$ and $\ell_{\mathrm{inv}}(w)=1$; otherwise $a_1=a(w)$ and $a_2\cdots a_k=\mathrm{rem}(w)$ is decomposed in the same way.

### 1.3 The symmetric group

$S_N$ is the group of permutations of $[N]=\lbrace 1,\dots,N\rbrace$, with $(uv)(i)=u(v(i))$. We write $w$ in one-line notation $w(1)w(2)\cdots w(N)$; "the value $b$ stands at position $j$" means $w(j)=b$. The generators are $s_i=(i\ \ i+1)$, $1\le i\le N-1$; $(S_N,\lbrace s_i\rbrace)$ is a Coxeter system of type $A_{N-1}$ and rank $N-1$, its reflections are the transpositions $(i\ j)$, and $\ell(w)$ is the number of inversions of $w$ ([BB, Ch. 1]).

Left multiplication acts on values and right multiplication on positions: the one-line notation of $(i\ j)w$ is that of $w$ with the **values** $i$ and $j$ exchanged, and that of $w(i\ j)$ is that of $w$ with the entries at the **positions** $i$ and $j$ exchanged.

For $w\in S_N$ put

```math
\mathrm{PInv}(w)=\lbrace (a,b) : a<b,\ w(a)>w(b)\rbrace ,\qquad
\mathrm{VInv}(w)=\lbrace (a,b) : a<b,\ w^{-1}(a)>w^{-1}(b)\rbrace .
```

$\mathrm{PInv}(w)$ is the set of inversions as pairs of positions. $\mathrm{VInv}(w)$ is the set of pairs of values $a<b$ such that $b$ stands to the left of $a$. By definition $\mathrm{VInv}(w)=\mathrm{PInv}(w^{-1})$. The map $(i,j)\mapsto (w(j),w(i))$ is a bijection from $\mathrm{PInv}(w)$ to $\mathrm{VInv}(w)$ with inverse $(a,b)\mapsto (w^{-1}(b),w^{-1}(a))$, so $\lvert \mathrm{PInv}(w)\rvert=\lvert \mathrm{VInv}(w)\rvert=\ell(w)$.

Further notation:

- $\mathrm{disp}(w)=\max_i \lvert w(i)-i\rvert$. A position $p$ with $w(p)-p=d>0$ is an *excedance of size $d$*; with $p-w(p)=d>0$ it is a *drop of size $d$*.
- $w_0$ is the permutation $i\mapsto N+1-i$, and $\varphi(x)=w_0xw_0$, so that $\varphi(x)(i)=N+1-x(N+1-i)$.
- Rank matrix: $x[i,k]=\mathrm{card}\lbrace j\le i : x(j)\ge k\rbrace$ for $i,k\in[N]$.

Everything in §2 uses only the order of the underlying set, so it applies to the permutations of any finite totally ordered set. This is used in §5.

### 1.4 Outside results used

**(E1) Gaetz–Gao, Cor. 1.2** ([arXiv:2610.22857v1](https://arxiv.org/abs/2610.22857v1)). Restated in our notation: let $W$ be any Coxeter group and $w\in W\setminus\lbrace e\rbrace$. Then

1. the interval $\lbrace u : u\le_R w\rbrace$ contains exactly one involution of maximal length, $a(w)$;
2. every involution $\iota$ with $\iota\le_R w$ satisfies $\iota\le_B a(w)$.

Their hypotheses are those stated: $W$ is an arbitrary Coxeter group (no finiteness or type assumption) and $w\ne e$; their $\le_R$ is the prefix order and their involutions exclude $e$. Here (E1) is applied with $W=S_N$ in §4 and with $W=S_{2n}$ in §5.

Standing of (E1). The preprint is unrefereed. **(E1) is proved in this file**: Appendix A writes out Gaetz–Gao's argument in our notation, filling in each step that the source leaves terse, from the textbook facts (E3) and Lemmas 2.4(2) and 2.6. None of these uses (E1), so there is no circularity. (E1) is therefore not an assumption of Theorems A and B, and the preprint is credited as the source of the argument, not relied on as an authority.

**(E2) Complementation step** (Gaetz–Gao, first paragraph of the proof of their Lemma 3.1): if $u,v\le_R w$ and $u\le_B v$, then $v^{-1}w\le_B u^{-1}w$. Valid in every Coxeter group. A full proof is given here as Lemma 2.6.

**(E3) Textbook facts**, from A. Björner and F. Brenti, *Combinatorics of Coxeter Groups*, GTM 231, Springer 2005, [doi:10.1007/3-540-27596-7](https://doi.org/10.1007/3-540-27596-7), cited as [BB] by chapter:

| Fact | Where in [BB] | Used in |
| --- | --- | --- |
| basic properties of length: $\ell(x^{-1})=\ell(x)$, $\ell(xy)\le\ell(x)+\ell(y)$, $\ell(xs)=\ell(x)\pm1$, $\ell(sx)=\ell(x)\pm1$ | Ch. 1 | §1.2, Lemmas 2.4, 2.6, Appendix A |
| $S_N$ is a Coxeter system of type $A_{N-1}$; $\ell$ = number of inversions; reflections = transpositions | Ch. 1 | §1.3 |
| definition of Bruhat order | Ch. 2 | §1.2 |
| subword property: if $s_1\cdots s_m$ is a reduced word for $y$, then $u\le_B y$ if and only if some subword $s_{i_1}\cdots s_{i_q}$ ($i_1<\dots<i_q$) is a reduced word for $u$ | Ch. 2 | Appendix A only |
| lifting property (one-sided form; the other side by inversion, see below) | Ch. 2 | Lemma 2.6, Appendix A |
| $B_n$ is a Coxeter system of type $B_n$ with $\ell_B=\mathrm{inv}+\mathrm{neg}+\mathrm{nsp}$ | Ch. 8 | §5 |
| $D_n$ is a Coxeter system of type $D_n$ with $\ell_D=\mathrm{inv}+\mathrm{nsp}$ | Ch. 8 | Remark 5.6 only |

Each fact is stated in full here, so the proofs do not depend on item numbers. The criteria usually quoted from [BB] for the weak order (the inversion-set criterion, Ch. 3) and for the Bruhat order of $S_N$ (the tableau criterion, Ch. 2) are **not** used as citations: the parts needed are proved in Lemmas 2.2 and 2.4.

**Lifting property**, in the form used. Let $u\le_B w$ and $s\in S$ with $\ell(ws)<\ell(w)$.

- (R1) If $\ell(us)<\ell(u)$, then $us\le_B ws$.
- (R2) If $\ell(us)>\ell(u)$, then $u\le_B ws$.
- (L1) Mirror image of (R1): if $u\le_B w$, $\ell(sw)<\ell(w)$ and $\ell(su)<\ell(u)$, then $su\le_B sw$.

The lifting property is stated in [BB, Ch. 2] for descents on one side. With right descents it reads: if $u<_B w$, $\ell(ws)<\ell(w)$ and $\ell(us)>\ell(u)$, then $u\le_B ws$ and $us\le_B w$. Either side suffices, because the other side follows by inversion. For instance, the left-handed form (if $u<_B w$, $\ell(sw)<\ell(w)$ and $\ell(su)>\ell(u)$, then $u\le_B sw$ and $su\le_B w$) follows from the right-handed one applied to $u^{-1},w^{-1}$, and conversely, because inversion preserves Bruhat order: a step $x\to xt$ with $\ell(xt)>\ell(x)$ becomes $x^{-1}\to tx^{-1}=x^{-1}(xtx^{-1})$ with $xtx^{-1}\in T$ and the same lengths. Also, inversion exchanges $\le_R$ and $\le_L$ (suffix order): $u\le_R w$ means $\ell(w)=\ell(u)+\ell(u^{-1}w)$, and $u^{-1}\le_L w^{-1}$ means $\ell(w^{-1})=\ell(w^{-1}u)+\ell(u^{-1})$, which is the same equation since $\ell(g^{-1})=\ell(g)$. (R2) is the right-handed statement, since $u\ne w$ when $s$ is a right descent of $w$ and not of $u$. For (R1): $us\le_B u$ because $u=(us)s$ is a step of the chain definition, so $us\le_B w$; $us\ne w$ since $\ell(us)<\ell(u)\le\ell(w)$; and $s$ is a right descent of $w$ and not of $us$; so the right-handed statement applied to $us<_B w$ gives $us\le_B ws$. (L1) is the same argument with the left-handed statement. (R1) and (L1) also hold when $u=w$, trivially.

### 1.5 Convention check against the source

The source's own worked example (type $A_6$) is reproduced by hand in §3.5 with the conventions above.

## 2. Preliminaries in the symmetric group

All statements in this section are proved here.

**Lemma 2.1 (inversions of a product).** Let $x,y\in S_N$ and let $i<j$ be positions. Write $\lbrace y(i),y(j)\rbrace=\lbrace a,b\rbrace$ with $a<b$. Then $(i,j)\in\mathrm{PInv}(xy)$ if and only if $(a,b)$ lies in exactly one of $\mathrm{PInv}(x)$, $\mathrm{VInv}(y)$. Consequently

```math
\ell(xy)=\ell(x)+\ell(y)-2\lvert \mathrm{PInv}(x)\cap\mathrm{VInv}(y)\rvert .
```

*Proof.* $(i,j)\in\mathrm{PInv}(xy)$ means $x(y(i))>x(y(j))$. The positions of $a,b$ in $y$ are $i,j$ in some order, and $(a,b)\in\mathrm{VInv}(y)$ means $y^{-1}(a)>y^{-1}(b)$, i.e. $y(i)=b$, $y(j)=a$. If $(a,b)\notin\mathrm{VInv}(y)$, then $y(i)=a$, $y(j)=b$ and the condition reads $x(a)>x(b)$, i.e. $(a,b)\in\mathrm{PInv}(x)$. If $(a,b)\in\mathrm{VInv}(y)$, the condition reads $x(b)>x(a)$, i.e. $(a,b)\notin\mathrm{PInv}(x)$. This proves the first statement. The map $\lbrace i,j\rbrace\mapsto\lbrace y(i),y(j)\rbrace$ is a bijection on 2-element subsets of $[N]$, so $\ell(xy)=\lvert\mathrm{PInv}(xy)\rvert$ is the size of the symmetric difference of $\mathrm{PInv}(x)$ and $\mathrm{VInv}(y)$, which is $\lvert\mathrm{PInv}(x)\rvert+\lvert\mathrm{VInv}(y)\rvert-2\lvert\mathrm{PInv}(x)\cap\mathrm{VInv}(y)\rvert$, and the two cardinalities are $\ell(x)$ and $\ell(y)$ (§1.3). ∎

Smallest cases: in $S_2$, $x=y=s_1$: both sets are $\lbrace(1,2)\rbrace$, the formula gives $1+1-2=0=\ell(e)$. In $S_3$, $x=s_1=213$, $y=s_2=132$: $\mathrm{PInv}(x)=\lbrace(1,2)\rbrace$, $\mathrm{VInv}(y)=\lbrace(2,3)\rbrace$, intersection empty, $\ell(xy)=\ell(231)=2$.

**Lemma 2.2 (prefix criterion).** For $u,w\in S_N$: $u\le_R w$ if and only if $\mathrm{VInv}(u)\subseteq\mathrm{VInv}(w)$. In particular an involution $\iota$ satisfies $\iota\le_R w$ if and only if $\mathrm{PInv}(\iota)\subseteq\mathrm{VInv}(w)$.

*Proof.* Apply Lemma 2.1 to $x=u^{-1}$, $y=w$, using $\mathrm{PInv}(u^{-1})=\mathrm{VInv}(u)$ and $\ell(u^{-1})=\ell(u)$:

```math
\ell(u^{-1}w)=\ell(u)+\ell(w)-2\lvert\mathrm{VInv}(u)\cap\mathrm{VInv}(w)\rvert .
```

So $\ell(u^{-1}w)=\ell(w)-\ell(u)$ holds if and only if $\lvert\mathrm{VInv}(u)\cap\mathrm{VInv}(w)\rvert=\ell(u)=\lvert\mathrm{VInv}(u)\rvert$, i.e. if and only if $\mathrm{VInv}(u)\subseteq\mathrm{VInv}(w)$. For an involution, $\mathrm{VInv}(\iota)=\mathrm{PInv}(\iota^{-1})=\mathrm{PInv}(\iota)$. ∎

(This is the type-A case of the general inversion-set criterion for the weak order, [BB, Ch. 3]. Not used.) Smallest case $N=2$: $s_1\le_R s_1$ and $\lbrace(1,2)\rbrace\subseteq\lbrace(1,2)\rbrace$; $s_1\not\le_R e$ and $\lbrace(1,2)\rbrace\not\subseteq\emptyset$.

**Lemma 2.3 (transpositions and length).** Let $x\in S_N$ and $i<j$. Then $\ell(x(i\ j))>\ell(x)$ if and only if $x(i)<x(j)$.

*Proof.* Let $y=x(i\ j)$, so $y(i)=x(j)$, $y(j)=x(i)$ and $y(c)=x(c)$ otherwise. Assume $x(i)<x(j)$ and compare inversions of $y$ and $x$ pair by pair. Pairs of positions disjoint from $\lbrace i,j\rbrace$ contribute equally. For $c<i$, the pairs $(c,i),(c,j)$ contribute $[x(c)>x(i)]+[x(c)>x(j)]$ to both (the two values compared with $x(c)$ are the same two values); likewise for $c>j$. For $i<c<j$ the pairs $(i,c),(c,j)$ contribute $[x(i)>x(c)]+[x(c)>x(j)]$ to $x$ and $[x(j)>x(c)]+[x(c)>x(i)]$ to $y$; since $x(j)>x(i)$ we have $[x(j)>x(c)]\ge[x(i)>x(c)]$ and $[x(c)>x(i)]\ge[x(c)>x(j)]$, so the contribution to $y$ is at least that to $x$. The pair $(i,j)$ contributes $0$ to $x$ and $1$ to $y$. Hence $\ell(y)\ge\ell(x)+1$. If instead $x(i)>x(j)$, then $y(i)<y(j)$ and $x=y(i\ j)$, so by the case just proved $\ell(x)\ge\ell(y)+1$. ∎

Smallest case $N=2$: $x=e$, $(i\ j)=(1\ 2)$: $x(1)<x(2)$ and $\ell(s_1)=1>0$; $x=s_1$: $x(1)>x(2)$ and $\ell(e)=0<1$.

**Lemma 2.4 (two consequences of the definition of Bruhat order).**

1. If $x\le_B y$ in $S_N$, then $x[i,k]\le y[i,k]$ for all $i,k\in[N]$.
2. In any Coxeter group, a suffix $r$ of $w$ satisfies $r\le_B w$.

*Proof.* (1) It suffices to treat one step $y=x(i\ j)$, $i<j$, $\ell(y)>\ell(x)$. By Lemma 2.3, $x(i)<x(j)$. For $m<i$ or $m\ge j$ the sets $\lbrace x(1),\dots,x(m)\rbrace$ and $\lbrace y(1),\dots,y(m)\rbrace$ are equal, so $y[m,k]=x[m,k]$. For $i\le m<j$ the second set is the first with $x(i)$ replaced by the larger value $x(j)$, so $y[m,k]-x[m,k]=[x(j)\ge k]-[x(i)\ge k]\ge0$.

(2) Let $w=cr$ with $\ell(w)=\ell(c)+\ell(r)$; induct on $\ell(c)$. If $c=e$ then $r=w$. Otherwise write $c=sc'$ with $s\in S$ and $\ell(c')=\ell(c)-1$ (first letter of a reduced word of $c$). Then $\ell(c'r)\le\ell(c')+\ell(r)=\ell(w)-1$ and $\ell(w)=\ell(s\cdot c'r)\le\ell(c'r)+1$, so $\ell(c'r)=\ell(w)-1=\ell(c')+\ell(r)$. By induction $r\le_B c'r$. Finally $w=s(c'r)=(c'r)t$ with $t=(c'r)^{-1}s(c'r)\in T$ and $\ell(w)>\ell(c'r)$, so $c'r\le_B w$. ∎

Only the implication (1) is needed; its converse is the non-trivial half of the tableau criterion ([BB, Ch. 2]) and is not used in the proofs. Smallest case of (1): $N=2$, $e\le_B s_1$, and the rank matrices are $e[1,1]=1,e[1,2]=0,e[2,1]=2,e[2,2]=1$ and $s_1[1,1]=1,s_1[1,2]=1,s_1[2,1]=2,s_1[2,2]=1$.

**Lemma 2.5 (the automorphism $\varphi$).** $\varphi(x)=w_0xw_0$ is a group automorphism of $S_N$ with $\varphi(s_i)=s_{N-i}$. It preserves $\ell$, $T$, involutions, $\le_R$ and $\le_B$, and $\mathrm{disp}(\varphi(x))=\mathrm{disp}(x)$. For $w\ne e$: $a(\varphi(w))=\varphi(a(w))$ and $\mathrm{rem}(\varphi(w))=\varphi(\mathrm{rem}(w))$. For every position $p$ and integer $d$: $x(p)=p-d$ if and only if $\varphi(x)(N+1-p)=(N+1-p)+d$.

*Proof.* Conjugation by $w_0$ sends the transposition $(i\ \ i+1)$ to $(w_0(i)\ \ w_0(i+1))=(N+1-i\ \ N-i)=s_{N-i}$. So $\varphi$ is an automorphism permuting $S$; it maps reduced words to reduced words, hence preserves $\ell$; it preserves $T$ (conjugates of $S$) and the order of elements. Since $\le_R$ and $\le_B$ are defined through the group law, $\ell$ and $T$, they are preserved. Therefore $\varphi$ maps $\mathcal{IP}(w)$ bijectively and length-preservingly onto $\mathcal{IP}(\varphi(w))$, which gives $a(\varphi(w))=\varphi(a(w))$ and then $\mathrm{rem}(\varphi(w))=\varphi(a(w))^{-1}\varphi(w)=\varphi(\mathrm{rem}(w))$. From $\varphi(x)(i)=N+1-x(N+1-i)$ we get $\varphi(x)(N+1-p)-(N+1-p)=p-x(p)$, which gives the last statement and, taking absolute values and the maximum over $p$, $\mathrm{disp}(\varphi(x))=\mathrm{disp}(x)$. ∎

Smallest cases: for $N=2$, $\varphi$ is the identity. For $N=3$, $\varphi(s_1)=s_2$ and $\varphi(231)=312$; in the table of §3.3, $a(231)=213$, $a(312)=132=\varphi(213)$, and $\mathrm{rem}(312)=213=\varphi(132)=\varphi(\mathrm{rem}(231))$.

**Lemma 2.6 (complementation; Gaetz–Gao, (E2)).** Let $(W,S)$ be a Coxeter system and $u,v,w\in W$ with $u\le_R w$, $v\le_R w$ and $u\le_B v$. Then $v^{-1}w\le_B u^{-1}w$.

*Proof* (the argument is Gaetz–Gao's; it is written out here so that the hypotheses can be checked). Induction on $\ell(v)$. If $\ell(v)=0$ then $v=e$, and $u\le_B e$ forces $\ell(u)\le0$, so $u=e$ and both sides equal $w$. Let $\ell(v)\ge1$ and choose $s\in S$ with $\ell(vs)=\ell(v)-1$. Put $a=v^{-1}w$ and $b=u^{-1}w$, so $\ell(w)=\ell(v)+\ell(a)=\ell(u)+\ell(b)$.

*Claim.* If $z\le_R w$ and $\ell(zs)=\ell(z)-1$, then $zs\le_R w$ and $\ell(s\cdot z^{-1}w)=\ell(z^{-1}w)+1$. Indeed $w=(zs)(s z^{-1}w)$, so $\ell(w)\le\ell(zs)+\ell(sz^{-1}w)\le(\ell(z)-1)+(\ell(z^{-1}w)+1)=\ell(w)$; both inequalities are equalities, which is the claim, with $(zs)^{-1}w=sz^{-1}w$.

By the claim for $z=v$: $vs\le_R w$, $(vs)^{-1}w=sa$ and $\ell(sa)=\ell(a)+1$.

*Case 1: $\ell(us)<\ell(u)$.* By the claim for $z=u$: $us\le_R w$, $(us)^{-1}w=sb$ and $\ell(sb)=\ell(b)+1$. By (R1), $us\le_B vs$. The induction hypothesis for $(us,vs,w)$, with $\ell(vs)<\ell(v)$, gives $sa\le_B sb$. Now $s$ is a left descent of both $sa$ and $sb$, so (L1) gives $a\le_B b$.

*Case 2: $\ell(us)>\ell(u)$.* By (R2), $u\le_B vs$. The induction hypothesis for $(u,vs,w)$ gives $sa\le_B b$. Also $a\le_B sa$, because $sa=a\cdot(a^{-1}sa)$ with $a^{-1}sa\in T$ and $\ell(sa)>\ell(a)$. Hence $a\le_B b$. ∎

Smallest case: $W=S_2$, $w=s_1$, $u=e$, $v=s_1$: $v^{-1}w=e\le_B s_1=u^{-1}w$.

**Lemma 2.7 (inversions of an involution; used only for the hand computations of §3.5 and Remark 5.6).** Let $\iota\in S_N$ be an involution.

1. If $(a\ b)$, $a<b$, is a 2-cycle of $\iota$, then $(a,b)\in\mathrm{PInv}(\iota)$.
2. If $(a\ b)$, $a<b$, is a 2-cycle of $\iota$ and $a<c<b$, then: $\iota(c)<b$ implies $(a,c)\in\mathrm{PInv}(\iota)$, and $\iota(c)>a$ implies $(c,b)\in\mathrm{PInv}(\iota)$.
3. If every 2-cycle of $\iota$ is of the form $(i\ \ i+1)$, then $\mathrm{PInv}(\iota)$ is the set of these pairs $(i,i+1)$.

*Proof.* (1) $\iota(a)=b>a=\iota(b)$. (2) $\iota(a)=b>\iota(c)$ in the first case, $\iota(c)>a=\iota(b)$ in the second. (3) Here $\iota(i)\le i+1$ and $\iota(j)\ge j-1$ for all $i,j$. If $i<j$ and $\iota(i)>\iota(j)$, then $i+1\ge\iota(i)>\iota(j)\ge j-1$, so $j=i+1$, and $\iota(i)>\iota(i+1)\ge i$ forces $\iota(i)=i+1$, i.e. $(i\ \ i+1)$ is a 2-cycle. Conversely each 2-cycle gives an inversion by (1). ∎

## 3. Small cases by hand

### 3.1 $S_1$ and $S_2$

$S_1=\lbrace e\rbrace$: $\ell_{\mathrm{inv}}=0=\mathrm{disp}=N-1$. $S_2$: $\ell_{\mathrm{inv}}(s_1)=1=\mathrm{disp}(s_1)=N-1$.

### 3.2 Procedure

By Lemma 2.2 the involution prefixes of $w$ are the involutions $\iota$ with $\mathrm{VInv}(\iota)\subseteq\mathrm{VInv}(w)$. So for each $w$ one lists $\mathrm{VInv}(w)$ (pairs $ab$ with $a<b$ and $b$ left of $a$), reads off the involution prefixes from the list of involutions, takes the longest one as $a(w)$, and obtains $\mathrm{rem}(w)=a(w)w$ by applying $a(w)$ to the values of $w$.

### 3.3 $S_3$

Involutions: $213$ with $\mathrm{VInv}=\lbrace12\rbrace$, $132$ with $\lbrace23\rbrace$, $321$ with $\lbrace12,13,23\rbrace$.

| $w$ | $\mathrm{VInv}(w)$ | involution prefixes | $a(w)$ | $\mathrm{rem}(w)$ | $\ell_{\mathrm{inv}}$ | $\mathrm{disp}(w)$ | $\mathrm{disp}(\mathrm{rem})$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 123 | – | – | – | – | 0 | 0 | – |
| 132 | 23 | 132 | 132 | 123 | 1 | 1 | 0 |
| 213 | 12 | 213 | 213 | 123 | 1 | 1 | 0 |
| 231 | 12 13 | 213 | 213 | 132 | 2 | 2 | 1 |
| 312 | 13 23 | 132 | 132 | 213 | 2 | 2 | 1 |
| 321 | 12 13 23 | 132 213 321 | 321 | 123 | 1 | 2 | 0 |

Distribution of $\ell_{\mathrm{inv}}$: $1,3,2$ for the values $0,1,2$; maximum $2=N-1$; $\ell_{\mathrm{inv}}\le\mathrm{disp}$ and $\mathrm{disp}(\mathrm{rem}(w))<\mathrm{disp}(w)$ in every row.

### 3.4 $S_4$

The nine involutions and their sets $\mathrm{VInv}$:

| $\iota$ | 2134 | 1324 | 1243 | 2143 | 3214 | 1432 | 3412 | 4231 | 4321 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $\mathrm{VInv}(\iota)$ | 12 | 23 | 34 | 12 34 | 12 13 23 | 23 24 34 | 13 14 23 24 | 12 13 14 24 34 | all six |
| $\ell(\iota)$ | 1 | 1 | 1 | 2 | 3 | 3 | 4 | 5 | 6 |

| $w$ | $\mathrm{VInv}(w)$ | involution prefixes | $a(w)$ | $\mathrm{rem}(w)$ | $\ell_{\mathrm{inv}}$ | $\mathrm{disp}(w)$ | $\mathrm{disp}(\mathrm{rem})$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1234 | – | – | – | – | 0 | 0 | – |
| 1243 | 34 | 1243 | 1243 | 1234 | 1 | 1 | 0 |
| 1324 | 23 | 1324 | 1324 | 1234 | 1 | 1 | 0 |
| 1342 | 23 24 | 1324 | 1324 | 1243 | 2 | 2 | 1 |
| 1423 | 24 34 | 1243 | 1243 | 1324 | 2 | 2 | 1 |
| 1432 | 23 24 34 | 1243 1324 1432 | 1432 | 1234 | 1 | 2 | 0 |
| 2134 | 12 | 2134 | 2134 | 1234 | 1 | 1 | 0 |
| 2143 | 12 34 | 1243 2134 2143 | 2143 | 1234 | 1 | 1 | 0 |
| 2314 | 12 13 | 2134 | 2134 | 1324 | 2 | 2 | 1 |
| 2341 | 12 13 14 | 2134 | 2134 | 1342 | 3 | 3 | 2 |
| 2413 | 12 14 34 | 1243 2134 2143 | 2143 | 1324 | 2 | 2 | 1 |
| 2431 | 12 13 14 34 | 1243 2134 2143 | 2143 | 1342 | 3 | 3 | 2 |
| 3124 | 13 23 | 1324 | 1324 | 2134 | 2 | 2 | 1 |
| 3142 | 13 23 24 | 1324 | 1324 | 2143 | 2 | 2 | 1 |
| 3214 | 12 13 23 | 1324 2134 3214 | 3214 | 1234 | 1 | 2 | 0 |
| 3241 | 12 13 14 23 | 1324 2134 3214 | 3214 | 1243 | 2 | 3 | 1 |
| 3412 | 13 14 23 24 | 1324 3412 | 3412 | 1234 | 1 | 2 | 0 |
| 3421 | 12 13 14 23 24 | 1324 2134 3214 3412 | 3412 | 1243 | 2 | 3 | 1 |
| 4123 | 14 24 34 | 1243 | 1243 | 3124 | 3 | 3 | 2 |
| 4132 | 14 23 24 34 | 1243 1324 1432 | 1432 | 2134 | 2 | 3 | 1 |
| 4213 | 12 14 24 34 | 1243 2134 2143 | 2143 | 3124 | 3 | 3 | 2 |
| 4231 | 12 13 14 24 34 | 1243 2134 2143 4231 | 4231 | 1234 | 1 | 3 | 0 |
| 4312 | 13 14 23 24 34 | 1243 1324 1432 3412 | 3412 | 2134 | 2 | 3 | 1 |
| 4321 | all six | all nine | 4321 | 1234 | 1 | 3 | 0 |

Distribution of $\ell_{\mathrm{inv}}$: $1,9,10,4$ for the values $0,\dots,3$ (the four elements with $\ell_{\mathrm{inv}}=3$ are $2341, 2431, 4123, 4213$). Maximum $3=N-1$. In every row the longest involution prefix is unique, $\ell_{\mathrm{inv}}(w)\le\mathrm{disp}(w)$ and $\mathrm{disp}(\mathrm{rem}(w))<\mathrm{disp}(w)$.

Status: proved (each row of both tables was derived by hand with the procedure above); numerical (the same tables are printed by `scripts/small_cases.py 6 tables`, archived in `scripts/out_small_cases_N6.txt`, and agree).

### 3.5 The source's worked example (convention check)

The source's example is $w=r_6r_3r_2r_1r_4r_5$ in type $A_6$, with ancestor decomposition $(r_3r_6)(r_2r_4)(r_1r_5)$. Take $r_i=s_i\in S_7$. Multiplying from the right end (start with $s_5=1234657$; then exchange the values $4,5$; then $1,2$; then $2,3$; then $3,4$; then $6,7$) gives $w=4125736$, with

```math
\mathrm{VInv}(w)=V=\lbrace (1,4),(2,4),(3,4),(3,5),(3,7),(6,7)\rbrace .
```

By Lemma 2.2 the involution prefixes of $w$ are the involutions $\iota$ with $\mathrm{PInv}(\iota)\subseteq V$. Let $\iota$ be one. By Lemma 2.7(1) each 2-cycle of $\iota$ is one of the six pairs of $V$. We exclude four of them with Lemma 2.7(2).

- $(1\ 4)$: for $c=2$ and $c=3$ the pair $(1,c)$ is not in $V$, so $\iota(c)>4$. Then $(2\ \iota(2))$ is a 2-cycle with $\iota(2)\ge5$, and $V$ contains no pair $(2,b)$ with $b\ge5$. Excluded.
- $(2\ 4)$: $(2,3)\notin V$, so $\iota(3)>4$, and the 2-cycle $(3\ \iota(3))$ is $(3\ 5)$ or $(3\ 7)$. Then $\iota(2)=4>3=\iota(\iota(3))$, so $(2,\iota(3))\in\mathrm{PInv}(\iota)$; it is $(2,5)$ or $(2,7)$, not in $V$. Excluded.
- $(3\ 5)$: take $c=4$. If $\iota(4)>3$, then $(4,5)\in\mathrm{PInv}(\iota)$, not in $V$. So $\iota(4)<3$ and $(\iota(4)\ 4)$ is the 2-cycle $(1\ 4)$ or $(2\ 4)$, both excluded. Excluded.
- $(3\ 7)$: take $c=4$. If $\iota(4)>3$, then $(4,7)\in\mathrm{PInv}(\iota)$, not in $V$. So $\iota(4)<3$, excluded as before.

So the 2-cycles of $\iota$ are among $(3\ 4)$ and $(6\ 7)$, and by Lemma 2.7(3) all three of $s_3$, $s_6$, $s_3s_6$ do have $\mathrm{PInv}\subseteq V$. Hence $\mathcal{IP}(w)=\lbrace s_3,s_6,s_3s_6\rbrace$, $a(w)=s_3s_6$, and $\mathrm{rem}(w)=3125647$ (exchange the values $3,4$ and $6,7$).

Next $\mathrm{VInv}(3125647)=\lbrace(1,3),(2,3),(4,5),(4,6)\rbrace$. A 2-cycle $(1\ 3)$ is excluded: $(1,2)$ is not in the set, so $\iota(2)>3$, and no pair $(2,b)$ with $b\ge4$ is in the set. A 2-cycle $(4\ 6)$ is excluded: if $\iota(5)>4$ then $(5,6)$ is an inversion, not in the set; so $\iota(5)<4$, and no pair $(b,5)$ with $b\le3$ is in the set. So the involution prefixes are $s_2$, $s_4$, $s_2s_4$ (Lemma 2.7(3)), the ancestor is $s_2s_4$, and the remainder is $2134657=s_1s_5$, an involution. This is the source's decomposition $(r_3r_6)(r_2r_4)(r_1r_5)$, with $\ell_{\mathrm{inv}}(w)=3$; also $\mathrm{disp}$ takes the values $3,2,1$ along $4125736$, $3125647$, $2134657$, in agreement with Corollary 4.5.

Status: proved (hand computation).

### 3.6 Lemma K and Proposition 4.3 replayed by hand

Two rows of the $S_4$ table (see §4 for the notation).

- $w=3412$, $D=2$, extremal excedances at $p=1$ ($v=3$) and $p=2$ ($v=4$). For $p=1$: $L=\lbrace3\rbrace$, $R=\lbrace1,2\rbrace$, $h(3)=1$, $h(2)=0$, so $c_0=2$, $J=\lbrace2,3\rbrace$ and $\iota=(2\ 3)=1324$, which is in the row's list of involution prefixes and moves $v=3$. For $p=2$: $L=\lbrace3,4\rbrace$, $R=\lbrace1,2\rbrace$, $h(4)=1,h(3)=2,h(2)=1,h(1)=0$, so $c_0=1$, $J=[1,4]$, $\iota=(1\ 3)(2\ 4)=3412=w$.
- $w=2431$, $D=3$; the only extremal entry is the drop $w(4)=1$. $\varphi(w)=4213$ has the extremal excedance $p=1$, $v=4$: $L=\lbrace4\rbrace$, $R=\lbrace1,2,3\rbrace$, $h(4)=1$, $h(3)=0$, $c_0=3$, $\iota=(3\ 4)=1243$, an involution prefix of $4213$ (table) with $\iota(4)=3$. The table gives $a(4213)=2143\ge_B 1243$ (one step $1243\to 1243\cdot(1\ 2)=2143$) and $\mathrm{rem}(4213)=3124$, whose first entry is $3\ne4$; and $\varphi(3124)=1342=\mathrm{rem}(2431)$ has no drop of size $3$, as Proposition 4.3 predicts.

## 4. Type A

### 4.1 Displacement and Bruhat order

**Lemma 4.1.** Let $d\ge1$. For $1\le i\le N-d$ put

```math
C_i^{d}=\lbrace x\in S_N : x[i,i+d]\ge1\rbrace=\lbrace x\in S_N : x(j)\ge i+d \text{ for some } j\le i\rbrace .
```

1. $C_i^d$ is an up-set of Bruhat order: $x\in C_i^d$ and $x\le_B y$ imply $y\in C_i^d$.
2. $x$ has an excedance of size at least $d$ if and only if $x\in C_i^d$ for some $i$.
3. If $\mathrm{disp}(x)\le d$ and $1\le i\le N-d$, then $x\in C_i^d$ if and only if $x(i)=i+d$.
4. $\mathrm{disp}$ is Bruhat-monotone: $x\le_B y$ implies $\mathrm{disp}(x)\le\mathrm{disp}(y)$.

*Proof.* (1) By Lemma 2.4(1), $y[i,i+d]\ge x[i,i+d]\ge1$.

(2) If $x(j)-j\ge d$, then $j\le N-d$ and $x\in C_j^d$ (take $i=j$). Conversely if $j\le i$ and $x(j)\ge i+d$, then $x(j)-j\ge i+d-j\ge d$.

(3) If $x(i)=i+d$ then $x\in C_i^d$ with $j=i$. Conversely let $j\le i$ with $x(j)\ge i+d$. Then $d\ge x(j)-j\ge i+d-j\ge d$, so all are equalities: $j=i$ and $x(i)=i+d$.

(4) Let $E_d$ be the set of permutations with an excedance of size at least $d$. By (2) and (1), $E_d$ is a union of up-sets, hence an up-set. By Lemma 2.5, $x$ has a drop of size at least $d$ if and only if $\varphi(x)\in E_d$, and $\varphi$ preserves Bruhat order, so $\varphi^{-1}(E_d)$ is an up-set. Hence $\lbrace x : \mathrm{disp}(x)\ge d\rbrace=E_d\cup\varphi^{-1}(E_d)$ is an up-set. If $x\le_B y$ and $d=\mathrm{disp}(x)\ge1$, then $y$ is in this set, so $\mathrm{disp}(y)\ge d$; if $\mathrm{disp}(x)=0$ the inequality holds because $\mathrm{disp}(y)\ge0$. ∎

Smallest case: $N=2$, $d=1$, $i=1$: $C_1^1=\lbrace s_1\rbrace$, an up-set; $e\le_B s_1$ and $\mathrm{disp}(e)=0\le1=\mathrm{disp}(s_1)$.

### 4.2 Lemma K

**Lemma 4.2 (Lemma K).** Let $w\in S_N$ with $D=\mathrm{disp}(w)\ge1$, and let $p$ be a position with $w(p)=p+D$. Put $v=w(p)$. Then there is an involution $\iota\in S_N$ with $\iota\le_R w$ and $\iota(v)<v$.

Explicitly: let $L=w(\lbrace1,\dots,p\rbrace)$ and $R=\lbrace1,\dots,v\rbrace\setminus L$. For $1\le c\le v$ put

```math
h(c)=\mathrm{card}\bigl(L\cap[c,v]\bigr)-\mathrm{card}\bigl(R\cap[c,v]\bigr).
```

Let $c_0$ be the largest $c$ with $h(c)=0$, $J=[c_0,v]$, and let $l_1<\dots<l_k$ and $r_1<\dots<r_k$ be the elements of $L\cap J$ and of $R\cap J$. Then $\iota$ is the product of the transpositions $(r_i\ l_i)$, $1\le i\le k$.

*Proof.*

**(a) $L\subseteq\lbrace1,\dots,v\rbrace$, $\lvert L\rvert=p$, $\lvert R\rvert=D$, $v\in L$.** Let $b>v$ be a value and $j$ its position. Then $b-j\le D$, so $j\ge b-D>v-D=p$. So the positions $1,\dots,p$ hold only values $\le v$: $L\subseteq\lbrace1,\dots,v\rbrace$. $\lvert L\rvert=p$ since $w$ is injective, hence $\lvert R\rvert=v-p=D$. $v=w(p)\in L$. Note that $\lbrace1,\dots,v\rbrace$ is the disjoint union of $L$ and $R$, and that an element of $\lbrace1,\dots,v\rbrace$ lies in $L$ if its position is $\le p$ and in $R$ if its position is $>p$.

**(b) Every $r\in R$ satisfies $r\ge p+1-D$.** The position $j$ of $r$ is $\ge p+1$, and $j-r\le D$.

**(c) $c_0$ exists, $c_0<v$, and $h(c)\ge1$ for $c_0<c\le v$.** We have $h(v)=1$ because $v\in L$, and $h(c-1)-h(c)=\pm1$ for $2\le c\le v$ because $c-1$ lies in exactly one of $L$, $R$. Let $c_1=\max(1,p+1-D)$. By (b) and $v\notin R$, $R\subseteq[c_1,v-1]$.

- If $p+1-D\ge1$, then $c_1=p+1-D$ and $[c_1,v]$ has $v-c_1+1=2D$ elements, of which $D$ are in $R$ and therefore $D$ in $L$. So $h(c_1)=0$.
- If $p+1-D<1$, then $c_1=1$ and $h(1)=\lvert L\rvert-\lvert R\rvert=p-D<0$. Since $h(v)=1>0$ and $h$ changes by $\pm1$ at each step, $h(c)=0$ for some $c$ with $1<c<v$.

In both cases $h$ has a zero in $[1,v]$; let $c_0$ be the largest. Then $c_0<v$ because $h(v)=1$. For $c_0<c\le v$, $h(c)\ne0$; since $h(v)=1$ and $h$ changes by $\pm1$, $h$ cannot become negative on $(c_0,v]$ without vanishing there, so $h(c)\ge1$.

**(d) $\lvert L\cap J\rvert=\lvert R\cap J\rvert=k\ge1$, $l_k=v$, and $r_i<l_i$ for all $i$.** The two cardinalities are equal because $h(c_0)=0$, and $k\ge1$, $l_k=v$ because $v=\max J$ lies in $L$. For $c_0\le c\le v-1$, taking complements inside $J$,

```math
\mathrm{card}\bigl(R\cap[c_0,c]\bigr)-\mathrm{card}\bigl(L\cap[c_0,c]\bigr)=\bigl(k-\mathrm{card}(R\cap[c+1,v])\bigr)-\bigl(k-\mathrm{card}(L\cap[c+1,v])\bigr)=h(c+1)\ge1 .
```

For $i<k$ apply this with $c=l_i\le v-1$: $\mathrm{card}(L\cap[c_0,l_i])=i$, so $\mathrm{card}(R\cap[c_0,l_i])\ge i+1$, hence $r_i\le l_i$, and $r_i\ne l_i$ since $L\cap R=\emptyset$. For $i=k$: $r_k\in J$ and $r_k\ne v$, so $r_k<v=l_k$.

**(e) $\iota$ is an involution with $\iota(v)<v$.** The $2k$ numbers $l_i,r_i$ are distinct, so $\iota$ is a product of $k\ge1$ disjoint transpositions, hence has order 2. $\iota(v)=\iota(l_k)=r_k<v$.

**(f) Every $(a,b)\in\mathrm{PInv}(\iota)$ has $a\in R\cap J$ and $b\in L\cap J$.** $\iota$ fixes every point outside $J$ and maps $J$ onto $J$. Let $a<b$ with $\iota(a)>\iota(b)$.

- If $a<c_0$: $\iota(a)=a$; if $b\in J$ then $\iota(b)\ge c_0>a$, and if $b\notin J$ then $\iota(b)=b>a$. No inversion.
- If $a>v$: then $b>v$ and both are fixed. No inversion.
- If $a\in J$ and $b\notin J$: then $b>v$, so $\iota(b)=b>v\ge\iota(a)$. No inversion.

So $a,b\in J$. If $a=l_i$, $b=l_j$ then $i<j$ and $\iota(a)=r_i<r_j=\iota(b)$: no inversion. If $a=r_i$, $b=r_j$ then $i<j$ and $\iota(a)=l_i<l_j=\iota(b)$: no inversion. If $a=l_i$, $b=r_j$: from $r_i<l_i<r_j$ we get $i<j$, so $\iota(a)=r_i<l_i<l_j=\iota(b)$: no inversion. The only remaining case is $a=r_j\in R$, $b=l_i\in L$.

**(g) $\iota\le_R w$.** Let $(a,b)\in\mathrm{PInv}(\iota)$. By (f), $a<b$, $a\in R$, $b\in L$. By (a) the position of $b$ is $\le p$ and the position of $a$ is $>p$, so $w^{-1}(a)>w^{-1}(b)$, i.e. $(a,b)\in\mathrm{VInv}(w)$. So $\mathrm{PInv}(\iota)\subseteq\mathrm{VInv}(w)$, and Lemma 2.2 gives $\iota\le_R w$. ∎

Smallest cases. $N=2$, $w=21$: $D=1$, $p=1$, $v=2$, $L=\lbrace2\rbrace$, $R=\lbrace1\rbrace$, $h(2)=1$, $h(1)=0$, $c_0=1$, $\iota=(1\ 2)=w$. $N=3$: $w=312$ ($p=1$, $v=3$, $L=\lbrace3\rbrace$, $R=\lbrace1,2\rbrace$, $h(3)=1$, $h(2)=0$, $c_0=2$, $\iota=(2\ 3)=132$, and indeed $a(312)=132$ in the table of §3); $w=321$ gives the same $L,R$ and $\iota=132\le_R 321$; $w=132$ ($p=2$, $v=3$, $L=\lbrace1,3\rbrace$, $R=\lbrace2\rbrace$, $\iota=(2\ 3)=w$); $w=213$ ($p=1$, $v=2$, $\iota=(1\ 2)=w$); $w=231$ has $D=2$ and no excedance of size 2 (its extremal entry is the drop $w(3)=1$, and $\varphi(231)=312$). The lemma does not involve $N$ beyond $v\le N$.

### 4.3 The remainder has smaller displacement

**Proposition 4.3.** Let $w\in S_N$, $w\ne e$, $D=\mathrm{disp}(w)$ and $r=\mathrm{rem}(w)$. Then $\mathrm{disp}(r)\le D-1$.

*Proof.* $D\ge1$ since $w\ne e$. Write $a=a(w)$.

*Step 1: $\mathrm{disp}(r)\le D$.* $w=a\cdot r$ with $\ell(w)=\ell(a)+\ell(r)$ because $a\le_R w$. So $r$ is a suffix of $w$, $r\le_B w$ by Lemma 2.4(2), and $\mathrm{disp}(r)\le\mathrm{disp}(w)$ by Lemma 4.1(4).

*Step 2: $r(p)\ne p+D$ for every position $p$.* This is statement $(\ast)$: *for every $w\in S_N\setminus\lbrace e\rbrace$ with $D=\mathrm{disp}(w)$, no position $p$ has $\mathrm{rem}(w)(p)=p+D$.* Suppose $r(p)=p+D$. Then $p\le N-D$ and $r\in C_p^D$ (take $j=i=p$ in the definition). Since $r\le_B w$ and $C_p^D$ is an up-set (Lemma 4.1(1)), $w\in C_p^D$; since $\mathrm{disp}(w)=D$, Lemma 4.1(3) gives $w(p)=p+D$. So Lemma 4.2 applies to $w$ and $p$ and gives an involution $\iota\le_R w$ with $\iota(p+D)\ne p+D$.

- By (E1)(2), applied to $W=S_N$, $w\ne e$ and the involution $\iota\le_R w$: $\iota\le_B a$.
- By Lemma 2.6, applied with its triple $(u,v,w)$ equal to $(\iota,a,w)$ (both $\iota\le_R w$ and $a\le_R w$, and $\iota\le_B a$): $a^{-1}w\le_B\iota^{-1}w$, i.e. $r\le_B\iota w$.
- $C_p^D$ is an up-set and $r\in C_p^D$, so $\iota w\in C_p^D$.
- $w=\iota\cdot(\iota w)$ with $\ell(w)=\ell(\iota)+\ell(\iota w)$ because $\iota\le_R w$ and $\iota^{-1}=\iota$. So $\iota w$ is a suffix of $w$; by Lemma 2.4(2) and Lemma 4.1(4), $\mathrm{disp}(\iota w)\le D$.
- By Lemma 4.1(3), $(\iota w)(p)=p+D$. But $(\iota w)(p)=\iota(w(p))=\iota(p+D)\ne p+D$. Contradiction.

The argument used nothing about $w$ except $w\ne e$, so it proves $(\ast)$ for every $w\in S_N\setminus\lbrace e\rbrace$.

*Step 3: $r(p)\ne p-D$ for every position $p$.* Suppose $r(p)=p-D$. Let $w'=\varphi(w)$. By Lemma 2.5, $w'\ne e$, $\mathrm{disp}(w')=D$, $\mathrm{rem}(w')=\varphi(r)$, and $\varphi(r)(N+1-p)=(N+1-p)+D$. This contradicts $(\ast)$ for $w'$ at the position $N+1-p$.

By Steps 1–3, $\lvert r(p)-p\rvert\le D-1$ for all $p$. ∎

Smallest cases: $N=2$: $w=s_1$, $r=e$, $\mathrm{disp}=0\le 0$. $N=3,4$: the last two columns of the tables in §3.

### 4.4 Theorem A

**Theorem 4.4 (Theorem A).** For every $N\ge1$ and $w\in S_N$: $\ell_{\mathrm{inv}}(w)\le\mathrm{disp}(w)\le N-1$.

*Proof.* The second inequality: $w(i),i\in[N]$, so $\lvert w(i)-i\rvert\le N-1$. The first, by induction on $\ell(w)$: for $w=e$ both sides are $0$. For $w\ne e$, $r=\mathrm{rem}(w)$ has $\ell(r)<\ell(w)$ (§1.2), so by induction and Proposition 4.3, $\ell_{\mathrm{inv}}(w)=1+\ell_{\mathrm{inv}}(r)\le1+\mathrm{disp}(r)\le1+(\mathrm{disp}(w)-1)=\mathrm{disp}(w)$. ∎

**Corollary 4.5.** Hart–Rowley Conjecture 1.3 holds for the Coxeter system $(S_N,\lbrace s_1,\dots,s_{N-1}\rbrace)$ of type $A_{N-1}$: $\ell_{\mathrm{inv}}(w)\le N-1=\mathrm{rank}$ for all $w$. Moreover, for $w\ne e$, along the ancestor decomposition $w=a_1\cdots a_k$ the remainders $w_i=a_i\cdots a_k$ satisfy $\mathrm{disp}(w_1)>\mathrm{disp}(w_2)>\dots>\mathrm{disp}(w_k)\ge1$.

*Proof.* The first statement is Theorem 4.4. The second is Proposition 4.3 applied to each $w_i$, $i<k$, with $w_{i+1}=\mathrm{rem}(w_i)$, and $w_k\ne e$. ∎

Status of 4.3–4.5: proved. They use (E1), which is proved in Appendix A; no outside hypothesis remains.

**Proposition 4.6 (the bound is attained; the Coxeter-element case is due to Hart–Rowley, §2 of the source).** Let $N\ge2$ and $c=s_1s_2\cdots s_{N-1}$. Then $\ell_{\mathrm{inv}}(c)=\mathrm{disp}(c)=N-1$.

*Proof.* For $1\le k\le N$ let $c_k=s_ks_{k+1}\cdots s_{N-1}$ (so $c_N=e$, $c_1=c$). We claim that $c_k$ is the permutation with one-line notation $1,\dots,k-1,k+1,k+2,\dots,N,k$. For $k=N$ this is $e$. If it holds for $k+1$, then $c_k=s_kc_{k+1}$ is obtained from $1,\dots,k,k+2,\dots,N,k+1$ by exchanging the values $k$ and $k+1$, which gives the claim for $k$. So $\mathrm{VInv}(c_k)=\lbrace(k,b) : k<b\le N\rbrace$ (the only value standing to the right of a larger value is $k$), and $\mathrm{disp}(c)=N-1$ from $c(N)=1$.

Let $k\le N-1$. First, $s_k\le_R c_k$ by Lemma 2.2, since $\mathrm{PInv}(s_k)=\lbrace(k,k+1)\rbrace\subseteq\mathrm{VInv}(c_k)$. Now let $\iota$ be any involution with $\iota\le_R c_k$, i.e. $\mathrm{PInv}(\iota)\subseteq\lbrace(k,b):b>k\rbrace$. So no pair $(i,j)$ with $i<j$ and $i\ne k$ is an inversion of $\iota$.

- For $i<k$: $\iota(i)<\iota(j)$ for all $j>i$. By induction on $i$, $\iota(i)=i$: if $\iota(i')=i'$ for all $i'<i$, then $\iota$ permutes $\lbrace i,\dots,N\rbrace$ and $\iota(i)$ is smaller than all other values taken there, so $\iota(i)=i$.
- Hence $\iota$ permutes $\lbrace k,\dots,N\rbrace$, and $\iota(k+1)<\iota(k+2)<\dots<\iota(N)$.
- If $\iota(k)=k$, then $\iota$ permutes $\lbrace k+1,\dots,N\rbrace$ increasingly, so it is the identity there and $\iota=e$, which is not an involution. So $\iota(k)=m>k$ and $\iota(m)=k$.
- The values taken on $\lbrace k+1,\dots,N\rbrace$ are $\lbrace k,\dots,N\rbrace\setminus\lbrace m\rbrace$, in increasing order, so $\iota(k+1)=k$. With $\iota(m)=k$ this gives $m=k+1$.
- Then $\iota$ permutes $\lbrace k+2,\dots,N\rbrace$ increasingly, so it is the identity there. Hence $\iota=s_k$.

So $\mathcal{IP}(c_k)=\lbrace s_k\rbrace$, $a(c_k)=s_k$ and $\mathrm{rem}(c_k)=s_kc_k=c_{k+1}$. Therefore $\ell_{\mathrm{inv}}(c_k)=1+\ell_{\mathrm{inv}}(c_{k+1})$ for $k\le N-1$, and $\ell_{\mathrm{inv}}(c)=\ell_{\mathrm{inv}}(c_1)=N-1$. ∎

Smallest case $N=2$: $c=s_1$, $\ell_{\mathrm{inv}}=1$. $N=4$: $c=2341$, row of §3, $\ell_{\mathrm{inv}}=3$.

## 5. Type B

### 5.1 Setup

Let $n\ge1$ and $I_n=\lbrace-n,\dots,-1,1,\dots,n\rbrace$ with its natural order. $S(I_n)$ is the group of permutations of $I_n$; the order-preserving bijection $\lambda:I_n\to[2n]$ ($\lambda(x)=n+x$ for $x>0$, $\lambda(x)=n+1+x$ for $x<0$) identifies it with $S_{2n}$, and all notions of §1.3 and §2 ($\ell_A$, $\mathrm{PInv}$, $\mathrm{VInv}$, $\le_R^A$, $\le_B$, rank matrices, $a_A$, $\mathrm{rem}_A$, $\ell_{\mathrm{inv}}^A$, $\mathrm{disp}_A$) are taken in $S(I_n)\cong S_{2n}$. The longest element of $S(I_n)$ is the order-reversing map $w_0(x)=-x$, and $\varphi(z)=w_0zw_0$.

```math
B_n=\lbrace w\in S(I_n) : w(-x)=-w(x)\text{ for all }x\rbrace=\lbrace w : \varphi(w)=w\rbrace .
```

With $s_0=(-1\ \ 1)$ and $s_i=(i\ \ i+1)(-i\ \ -i-1)$ for $1\le i\le n-1$, $(B_n,\lbrace s_0,\dots,s_{n-1}\rbrace)$ is a Coxeter system of type $B_n$ and rank $n$, and its length function is ([BB] Ch. 8)

```math
\ell_B(w)=\mathrm{inv}(w)+\mathrm{neg}(w)+\mathrm{nsp}(w),
```

where, for the window $(w(1),\dots,w(n))$: $\mathrm{inv}(w)=\mathrm{card}\lbrace 1\le i<j\le n : w(i)>w(j)\rbrace$, $\mathrm{neg}(w)=\mathrm{card}\lbrace i\in[n] : w(i)<0\rbrace$, $\mathrm{nsp}(w)=\mathrm{card}\lbrace 1\le i<j\le n : w(i)+w(j)<0\rbrace$. We write $\le_R^B$, $a_B$, $\mathrm{rem}_B$, $\ell_{\mathrm{inv}}^B$ for the notions of §1.2 in this Coxeter system. A pair $(a,b)$, $a<b$, of elements of $I_n$ is *symmetric* if $a=-b$; for a set $P$ of pairs, $F(P)$ is its set of symmetric pairs.

Two facts about $\mathrm{neg}$ for $z\in B_n$:

- (N1) $\mathrm{neg}(z)=\lvert F(\mathrm{PInv}(z))\rvert=\lvert F(\mathrm{VInv}(z))\rvert$. Indeed $(-c,c)\in\mathrm{PInv}(z)$ iff $z(-c)>z(c)$ iff $z(c)<0$, which gives the first equality; applied to $z^{-1}\in B_n$ it gives $\lvert F(\mathrm{VInv}(z))\rvert=\mathrm{neg}(z^{-1})$, and $\mathrm{neg}(z^{-1})=\mathrm{card}\lbrace c>0:z^{-1}(c)<0\rbrace=\mathrm{card}\lbrace d<0:z(d)>0\rbrace=\mathrm{card}\lbrace d>0:z(d)<0\rbrace=\mathrm{neg}(z)$ (substitute $d=z^{-1}(c)$, then use $z(-d)=-z(d)$).
- (N2) In the labelling by $[2n]$, $\mathrm{neg}(z)$ is the rank-matrix entry $z[n,n+1]$: this entry counts the negative positions $x$ with $z(x)>0$, and $\mathrm{card}\lbrace c>0 : z(-c)>0\rbrace=\mathrm{card}\lbrace c>0:z(c)<0\rbrace$.

### 5.2 Lengths and weak order

**Lemma 5.1.** For $w\in B_n$: $\ell_A(w)=2\ell_B(w)-\mathrm{neg}(w)$.

*Proof.* Count the inversions of $w$ as a permutation of $I_n$ by the signs of the two positions. Two positive positions $i<j$: inversion iff $w(i)>w(j)$; total $\mathrm{inv}(w)$. Two negative positions $-j<-i$ ($0<i<j$): inversion iff $w(-j)>w(-i)$ iff $w(i)>w(j)$; total $\mathrm{inv}(w)$. A negative position $-i$ and a positive position $j$: inversion iff $w(-i)>w(j)$ iff $w(i)+w(j)<0$. For $i=j$ this gives $\mathrm{neg}(w)$ pairs; for $i\ne j$ the condition is symmetric in $(i,j)$, so the ordered pairs $(i,j)$, $i\ne j$, give $2\mathrm{nsp}(w)$. Total: $2\mathrm{inv}+\mathrm{neg}+2\mathrm{nsp}=2\ell_B-\mathrm{neg}$. ∎

**Lemma 5.2.** For $x,y\in B_n$ let $P=\mathrm{PInv}(x)\cap\mathrm{VInv}(y)$. Then

```math
\ell_B(xy)=\ell_B(x)+\ell_B(y)-\lvert P\rvert-\lvert F(P)\rvert .
```

Consequently $\ell_B(xy)=\ell_B(x)+\ell_B(y)$ iff $P=\emptyset$ iff $\ell_A(xy)=\ell_A(x)+\ell_A(y)$; and for $u,w\in B_n$: $u\le_R^B w$ if and only if $u\le_R^A w$. (The right weak order of $B_n$ is the restriction of the right weak order of $S_{2n}$.)

*Proof.* By Lemma 2.1, $\ell_A(xy)=\ell_A(x)+\ell_A(y)-2\lvert P\rvert$. For $\mathrm{neg}(xy)$ use the first statement of Lemma 2.1: a pair of positions $i<j$ with $\lbrace y(i),y(j)\rbrace=\lbrace a<b\rbrace$ is in $\mathrm{PInv}(xy)$ iff $(a,b)$ is in exactly one of $\mathrm{PInv}(x)$, $\mathrm{VInv}(y)$. Since $y\in B_n$, the bijection $\lbrace i,j\rbrace\mapsto\lbrace y(i),y(j)\rbrace$ maps symmetric pairs to symmetric pairs ($\lbrace y(-c),y(c)\rbrace=\lbrace -y(c),y(c)\rbrace$), and so does its inverse, given by $y^{-1}\in B_n$. Therefore $\lvert F(\mathrm{PInv}(xy))\rvert$ is the size of the symmetric difference of $F(\mathrm{PInv}(x))$ and $F(\mathrm{VInv}(y))$, and by (N1)

```math
\mathrm{neg}(xy)=\mathrm{neg}(x)+\mathrm{neg}(y)-2\lvert F(P)\rvert .
```

By Lemma 5.1, $2\ell_B(xy)=\ell_A(xy)+\mathrm{neg}(xy)=2\ell_B(x)+2\ell_B(y)-2\lvert P\rvert-2\lvert F(P)\rvert$. This is the formula. Since $F(P)\subseteq P$, $\lvert P\rvert+\lvert F(P)\rvert=0$ iff $P=\emptyset$, and by Lemma 2.1 $P=\emptyset$ iff $\ell_A(xy)=\ell_A(x)+\ell_A(y)$. For the weak orders take $x=u$, $y=u^{-1}w$: $u\le_R^B w$ means $\ell_B(xy)=\ell_B(x)+\ell_B(y)$ and $u\le_R^A w$ means $\ell_A(xy)=\ell_A(x)+\ell_A(y)$. ∎

Smallest cases. $n=1$: $B_1=\lbrace e,s_0\rbrace=S(I_1)\cong S_2$; $\ell_A(s_0)=1=2\cdot1-1$. $n=2$, $x=y=s_0$ (in labels $1324$): $P=\lbrace(2,3)\rbrace$ in labels, i.e. the symmetric pair $(-1,1)$, so $\lvert P\rvert+\lvert F(P)\rvert=2$ and $\ell_B(e)=1+1-2$.

### 5.3 Ancestors

**Proposition 5.3.** Let $w\in B_n$, $w\ne e$. Then $a_A(w)\in B_n$ and $a_B(w)=a_A(w)$: the $\ell_B$-longest involution prefix of $w$ in $B_n$ is unique and equals the ancestor of $w$ in $S_{2n}$.

*Proof.* Let $a=a_A(w)$, the unique $\ell_A$-longest involution with $a\le_R^A w$ ((E1)(1) for $S_{2n}$). By Lemma 2.5, $\varphi(a)$ is an involution with $\varphi(a)\le_R^A\varphi(w)=w$ and $\ell_A(\varphi(a))=\ell_A(a)$. By uniqueness $\varphi(a)=a$, so $a\in B_n$. By Lemma 5.2, $a\le_R^B w$; so $a$ is an involution prefix of $w$ in $B_n$.

Let $\iota\in B_n$ be any involution with $\iota\le_R^B w$. By Lemma 5.2, $\iota\le_R^A w$. By (E1)(2) for $W=S_{2n}$: $\iota\le_B a$ in the Bruhat order of $S_{2n}$. Hence $\ell_A(\iota)\le\ell_A(a)$, with equality only if $\iota=a$ (§1.2), and $\mathrm{neg}(\iota)=\iota[n,n+1]\le a[n,n+1]=\mathrm{neg}(a)$ by (N2) and Lemma 2.4(1). By Lemma 5.1,

```math
2\ell_B(\iota)=\ell_A(\iota)+\mathrm{neg}(\iota)\le\ell_A(a)+\mathrm{neg}(a)=2\ell_B(a),
```

with equality only if $\ell_A(\iota)=\ell_A(a)$, i.e. only if $\iota=a$. So $a$ is the unique $\ell_B$-longest involution prefix of $w$ in $B_n$. ∎

(The uniqueness of $a_B(w)$ is also (E1)(1) for the Coxeter group $B_n$; the proof above does not need that.)

**Theorem 5.4 (Theorem B).** For every $w\in B_n$: $\ell_{\mathrm{inv}}^B(w)=\ell_{\mathrm{inv}}^A(w)\le\mathrm{disp}_A(w)\le 2n-1$. In terms of the window of $w$,

```math
\mathrm{disp}_A(w)=\max_{1\le i\le n}\ \delta_i,\qquad \delta_i=\lvert w(i)-i\rvert \text{ if } w(i)>0,\qquad \delta_i=i+\lvert w(i)\rvert-1 \text{ if } w(i)<0 .
```

*Proof.* Equality, by induction on $\ell_B(w)$: for $w=e$ both sides are $0$. For $w\ne e$, Proposition 5.3 gives $a_B(w)=a_A(w)=:a\in B_n$, so $\mathrm{rem}_B(w)=a^{-1}w=\mathrm{rem}_A(w)\in B_n$, with $\ell_B(\mathrm{rem}_B(w))<\ell_B(w)$. By induction $\ell_{\mathrm{inv}}^B(w)=1+\ell_{\mathrm{inv}}^B(a^{-1}w)=1+\ell_{\mathrm{inv}}^A(a^{-1}w)=\ell_{\mathrm{inv}}^A(w)$. The inequalities are Theorem 4.4 with $N=2n$. For the formula: the displacement at the position $\lambda(i)$ is $\lvert\lambda(w(i))-\lambda(i)\rvert$, which is $\lvert w(i)-i\rvert$ if $i,w(i)>0$ and $\lvert n+1+w(i)-n-i\rvert=i+\lvert w(i)\rvert-1$ if $i>0>w(i)$; the positions $\lambda(-i)$ give the same numbers because $w(-i)=-w(i)$ and $\lambda(-x)=2n+1-\lambda(x)$. ∎

Status of 5.1–5.4: proved. 5.3 and 5.4 use (E1) for $W=S_{2n}$, which is proved in Appendix A; no outside hypothesis remains.

**Smallest cases.** $n=1$: $B_1=S_2$. $n=2$: in labels, $B_2$ consists of the eight $w\in S_4$ with $w(5-i)=5-w(i)$: $1234$, $1324=s_0$, $2143=s_1$, $3142=s_0s_1$, $2413=s_1s_0$, $3412=s_0s_1s_0$, $4231=s_1s_0s_1$, $4321$. $B_2$ is dihedral of order 8: its two elements of length 2 have the single involution prefix $s_0$, resp. $s_1$, and involution length 2; the other five non-identity elements are involutions. The $S_4$ table of §3 gives the same: $a(3142)=1324$, $\mathrm{rem}=2143$, $\ell_{\mathrm{inv}}=2$; $a(2413)=2143$, $\mathrm{rem}=1324$, $\ell_{\mathrm{inv}}=2$; and $1324,2143,3412,4231,4321$ are involutions. So $\ell_{\mathrm{inv}}^B=\ell_{\mathrm{inv}}^A$ on $B_2$, with distribution $1,5,2$.

### 5.4 What Theorem B does not give

**Remark 5.5 (type B: the rank bound is open).** $\mathrm{rank}(B_n)=n$, and Conjecture 1.3 for $B_n$ asserts $\ell_{\mathrm{inv}}^B\le n$. Theorem 5.4 gives $2n-1$, and Theorem 4.4 cannot give more on its own: $\mathrm{disp}_A$ takes the value $2n-1$ on $B_n$ (for instance at $w=-\mathrm{id}$, where $\ell_{\mathrm{inv}}=1$). Numerical evidence for the rank bound (not used in any proof): the archived output `scripts/out_typeB_cases_n4.txt` gives $\max\ell_{\mathrm{inv}}^{B_n}=n$ for $n\le4$ (distributions $B_2$: 1, 5, 2; $B_3$: 1, 19, 20, 8; $B_4$: 1, 75, 156, 123, 29). No proof of the rank bound for $B_n$ is claimed here.

**Remark 5.6 (type D: the reduction fails).** $D_n=\lbrace w\in B_n:\mathrm{neg}(w)\text{ even}\rbrace$ is a Coxeter system with generators $s_0^D=(-2,-1,3,\dots,n)$ (window notation), $s_1,\dots,s_{n-1}$ and length $\ell_D=\mathrm{inv}+\mathrm{nsp}=\ell_B-\mathrm{neg}$ ([BB] Ch. 8). The computation of Lemma 5.2 then gives $\ell_D(xy)=\ell_D(x)+\ell_D(y)-(\lvert P\rvert-\lvert F(P)\rvert)$: a product in $D_n$ is $\ell_D$-additive exactly when $P$ consists of symmetric pairs, which is weaker than $P=\emptyset$. So for $u,w\in D_n$, $u\le_R^A w$ implies $u\le_R^D w$, but the converse fails: there are 32 pairs $(u,w)$ in $D_3$ and 788 in $D_4$ with $u\le_R^D w$ and not $u\le_R^A w$ (numerical; `scripts/out_typeB_cases_n4.txt`, key `prefix_in_W_not_in_A`). Moreover $D_n$ is not the fixed-point set of $\varphi$, and the ancestor in $S_{2n}$ of an element of $D_n$ need not lie in $D_n$, as the next example shows. So neither half of the argument of §5.3 is available.

Explicit example in $D_3$ (status: proved, by hand; numerical: `scripts/out_typeB_cases_n4.txt`, row `D3`, key `shortest_mismatch`). Let $w=s_2s_0^D$, with window $(-3,-1,2)$.

- In $D_3$: $\ell_D(w)=\mathrm{inv}(w)+\mathrm{nsp}(w)=0+2=2$, so every reduced word of $w$ has two letters, and its first letter $s$ satisfies $\ell_D(sw)=1$, i.e. $sw$ is a generator. Now $s_2w=s_0^D$, while $s_1w$ has window $(-3,-2,1)$ and $s_0^Dw$ has window $(-3,2,-1)$, neither of which is a generator. So $s_2s_0^D$ is the only reduced word, and the prefixes of $w$ are $e$, $s_2$, $w$. $w$ is not an involution ($w(w(1))=w(-3)=-2\ne1$). Hence $a_D(w)=s_2$, $\mathrm{rem}_D(w)=s_0^D$ and $\ell_{\mathrm{inv}}^{D}(w)=2$.
- In $S_6$: in labels $w=246135$ and $\mathrm{VInv}(w)=V=\lbrace(1,2),(1,4),(1,6),(3,4),(3,6),(5,6)\rbrace$. Let $\iota$ be an involution with $\mathrm{PInv}(\iota)\subseteq V$ and use Lemma 2.7. A 2-cycle $(1\ 4)$ or $(1\ 6)$ is excluded with $c=2$: $\iota(2)\ne1$, so $\iota(2)>1$ and $(2,4)$, resp. $(2,6)$, is an inversion, not in $V$. A 2-cycle $(3\ 6)$ is excluded with $c=4$: if $\iota(4)>3$ then $(4,6)$ is an inversion, not in $V$; so $\iota(4)<3$, and the 2-cycle $(\iota(4)\ 4)$ is $(1\ 4)$, excluded, or $(2\ 4)$, not in $V$. So the 2-cycles are among $(1\ 2),(3\ 4),(5\ 6)$, and by Lemma 2.7(3) every product of them is an involution prefix. Hence $a_A(w)=214365$. Its window is $(-1,3,2)$, with one negative entry, so $a_A(w)\in B_3\setminus D_3$.
- Continuing in $S_6$: $\mathrm{rem}_A(w)=135246$, with $\mathrm{VInv}=\lbrace(2,3),(2,5),(4,5)\rbrace$. A 2-cycle $(2\ 5)$ is excluded with $c=3$: if $\iota(3)>2$ then $(3,5)$ is an inversion, not in the set; so $\iota(3)=1$, and $(1,3)$ is not in the set. So the ancestor is $132546$ and the remainder is $124356=s_3^A$. Hence $\ell_{\mathrm{inv}}^{A}(w)=3\ne2=\ell_{\mathrm{inv}}^{D}(w)$.

Counts (numerical, exact): in $D_3$, 6 of the 23 non-identity elements have $a_D\ne a_A$ and 4 have $\ell_{\mathrm{inv}}^D\ne\ell_{\mathrm{inv}}^A$; in $D_4$, 55 of 191 and 33 (`scripts/out_typeB_cases_n4.txt`). Nothing is claimed about Conjecture 1.3 for type D.

## 6. Table of statements, status and evidence

**proved** = argued in full in this file. **numerical** = exact finite computation by one program, corroboration only; no proof step uses it. "N6", "N7" and "B4" are the archived outputs `scripts/out_small_cases_N6.txt`, `scripts/out_small_cases_N7.txt` and `scripts/out_typeB_cases_n4.txt`.

| Statement | Basis | Status | Numerical corroboration |
| --- | --- | --- | --- |
| (E1) Gaetz–Gao Cor. 1.2 | proved here (Theorem A.5; Gaetz–Gao's argument, credited) | proved | N7: `ancestor_ties` = 0, `gg1_fail` = 0 for $N\le7$ |
| (E3) textbook facts of §1.4 | cited ([BB], by chapter) | cited | – |
| A.1–A.5 | proved here (Gaetz–Gao's argument) | proved | – |
| Lemma 2.1, inversions of a product | proved here | proved | N6: `lemma21_length_formula_fail` = 0, all pairs, $N\le5$ |
| Lemma 2.2, prefix criterion | proved here | proved | N6: `lemma22_prefix_criterion_fail` = 0, all pairs, $N\le5$ |
| Lemma 2.3 | proved here | proved | – |
| Lemma 2.4, Bruhat order and rank matrices; suffixes | proved here | proved | N6: `lemma24_bruhat_implies_rank_fail` = 0 ($N\le6$, Bruhat order built from the definition), `suffix_not_bruhat_below_fail` = 0 |
| Lemma 2.5, the automorphism $\varphi$ | proved here (uses (E1)(1)) | proved | N7: `phi_equivariance_fail` = 0 |
| Lemma 2.6, complementation (Gaetz–Gao) | proved here (their argument) | proved | N6/N7: `gg2_fail` = 0 |
| Lemma 2.7 | proved here | proved | – |
| §3 small cases | by hand | proved | N6 (printed tables) |
| Lemma 4.1 | proved here | proved | N7: `disp_suffix_fail` = 0 |
| Lemma 4.2 (Lemma K) | proved here | proved | N7: `lemmaK_fail` = 0 on 2, 8, 38, 182, 1064, 7176 tests ($N=2,\dots,7$; one per position of largest displacement) |
| Proposition 4.3 | proved here | proved | N7: `disp_rem_not_smaller` = 0 ($N\le7$) |
| **Theorem 4.4 = Theorem A** | proved here | proved | N7: `linv_gt_disp` = 0, `linv_gt_rank` = 0 ($N\le7$) |
| Corollary 4.5 | from 4.3, 4.4 | proved | – |
| Proposition 4.6 (case due to Hart–Rowley) | proved here | proved | N7: `coxeter_element_linv`, `coxeter_chain_ok` |
| Lemmas 5.1, 5.2 | proved here | proved | B4: `b1_*` = 0; `b2_length_of_product_fail` = 0 ($n\le3$), `weak_order_mismatch` = 0 ($n\le4$) |
| Proposition 5.3 | proved here | proved | B4: `ancestor_mismatch_W_vs_A` = 0, $n\le4$ |
| **Theorem 5.4 = Theorem B** | proved here | proved | B4: `linv_mismatch_W_vs_A` = 0, `linvW_gt_dispA` = 0, $n\le4$ |
| Remark 5.5: Conjecture 1.3 for $B_n$ (bound $n$) | not addressed | open | B4, $n\le4$ only |
| Remark 5.6: the reduction fails for $D_n$ | example by hand; counts by script | proved (example); numerical (counts) | B4, rows `D3`, `D4` |

## 7. Computations

**No step of the proofs in §§2–5 and Appendix A depends on a computation.** The scripts in `scripts/` check the lemmas on small cases. They are a single implementation (Python 3.12.10, standard library only, exact integer arithmetic), so they are not certificates in the sense of two independent programs that agree. Exact commands and expected output are in `REPRODUCE.md`. Every failure counter in the three archived outputs is 0.

`small_cases.py` computes ancestors from the definition (length-additivity against every involution), tests Lemmas 2.1, 2.2 on all pairs for $N\le5$, builds Bruhat order from the chain definition for $N\le6$ and compares it with rank matrices, tests instances of (E1)(2) and Lemma 2.6, Lemma K with the replay of Step 2 of Proposition 4.3, and Proposition 4.6. For $N=7$ it uses the rank-matrix criterion ([BB, Ch. 2]) for Bruhat comparisons. `typeB_cases.py` computes $\ell_B$ and $\ell_D$ as word lengths by breadth-first search, and ancestors in $S_{2n}$ over all involutions of $S_{2n}$. The scripts are copied unedited; their comments refer to this file under its working name `PROOF_typeA.md`.

## 8. Remarks on the proof

1. **Direction of the conventions.** Left multiplication acts on values. "Prefix" is intrinsic to the Coxeter system, so the theorem does not depend on this choice; the tables of §3 do. The convention is checked against the source's own example in §3.5.
2. **Where each hypothesis of Lemma K is used.** $w(p)=p+D$ with $D=\mathrm{disp}(w)$ is used in (a) (no larger value to the left of $p$) and (b) (values of $R$ are not far to the left); both fail for a non-extremal excedance. The lemma handles excedances only; drops are transported by $\varphi$ in Step 3 of Proposition 4.3, which needs $a(\varphi(w))=\varphi(a(w))$ (Lemma 2.5).
3. **§3 and the example of Remark 5.6 are checks, not parts of the proof.** An error in the tables of §3, in the convention check of §3.5 or in the $D_3$ example would not affect the proofs of Theorems A or B.
4. **Type B.** Theorem B is a statement about $\ell_{\mathrm{inv}}$ in the Coxeter system $B_n$, but its bound is $2n-1$, not the rank $n$. It is not Conjecture 1.3 for type B.
5. **Formalization.** Nothing in this file is formalized in Lean. (E1), Lemma 2.6, Appendix A and the Bruhat order of a general Coxeter group are not in Mathlib as far as we know. The only Lean formalization connected with Conjecture 1.3 in this repository is the $H_3$ counterexample of the companion note.

## 9. Status

Disclosure: Human researcher assisted by AI systems during search, proof drafting, and verification. Internally audited and computationally verified, externally unrefereed.

**Expert verification required.** Verified so far with:

- Lean 4 and Mathlib: nothing in this file.
- Independent certificates: none; no proof step uses a computation. One exact implementation (two Python scripts) corroborates the results for $S_N$, $N\le7$, and for $B_n$, $D_n$, $n\le4$ (§7).
- AI referees, reading the full proof: round 1 (commit `4101231`), Claude Opus: minor issues (five), all fixed; Claude Fable: accept, no issues. In round 1, (E1) was still cited, not proved. Round 2 (commit `599cc84`), after Appendix A was added: Claude Opus, Claude Fable and Claude Sonnet each accept with no issues.

## Appendix A. Proof of (E1), after Gaetz–Gao

**Credit.** Everything in this appendix is the argument of C. Gaetz and Y. Gao, [arXiv:2610.22857v1](https://arxiv.org/abs/2610.22857v1), §3: their Lemma 3.1, Theorem 1.1 and Corollary 1.2. It is written out again here in our notation, in our own words, with each step that the source states briefly justified in full, so that Theorems A and B depend only on the textbook facts of §1.4. The first paragraph of their proof of Lemma 3.1 is Lemma 2.6 above. The only change of substance: the finiteness of $[e,x]_B$ used in their proof of Theorem 1.1 is replaced by a length argument (A.4), which needs nothing beyond §1.2. Status of every statement in this appendix: **proved**.

**Setting.** $(W,S)$ is an arbitrary Coxeter system. Besides $\le_R$ and $\le_B$ (§1.2) we use the **left weak order** (suffix order)

```math
u\le_L x \iff \ell(x)=\ell(xu^{-1})+\ell(u),
```

so "$u$ is a suffix of $x$" (§1.2) means $u\le_L x$. Inputs, all from §1.2, §1.4 and §2:

- (I1) basic length facts: $\ell(g^{-1})=\ell(g)$, $\ell(gh)\le\ell(g)+\ell(h)$, $\ell(gs)=\ell(g)\pm1$ for $s\in S$;
- (I2) $g\le_B h$ and $g\ne h$ imply $\ell(g)<\ell(h)$ (§1.2);
- (I3) the lifting property in the forms (R1), (R2) of §1.4;
- (I4) the subword property (§1.4, (E3));
- (I5) inversion preserves $\le_B$ and exchanges $\le_L$ and $\le_R$ (§1.4, paragraph on the lifting property);
- (I6) $u\le_L x$ implies $u\le_B x$ (Lemma 2.4(2)); hence, by (I5), $u\le_R x$ implies $u\le_B x$ (apply Lemma 2.4(2) to $u^{-1}\le_L x^{-1}$ and invert);
- (I7) complementation, Lemma 2.6: if $u,v\le_R x$ and $u\le_B v$, then $v^{-1}x\le_B u^{-1}x$.

For $x,y\in W$ put

```math
\mathcal{LB}(x,y)=\lbrace g : g\le_L x,\ g\le_B y\rbrace,\qquad
\mathcal{BR}(x,y)=\lbrace g : g\le_B x,\ g\le_R y\rbrace,\qquad
\mathcal{LR}(x,y)=\lbrace g : g\le_L x,\ g\le_R y\rbrace .
```

Each contains $e$. A *greatest element* of a set is an element of the set that is $\ge_B$ every element of it; it is unique if it exists, because $\le_B$ is a partial order.

### A.1 The downward Demazure operator

For $s\in S$ and $g\in W$ let $\delta_s(g)=gs$ if $\ell(gs)<\ell(g)$, and $\delta_s(g)=g$ otherwise. (Since $g$ and $gs$ differ by the reflection $s$, the shorter of the two is Bruhat-below the other, so $\delta_s(g)$ is the Bruhat-smaller of $g$ and $gs$; this is Gaetz–Gao's $\min_{\le_B}\lbrace g,gs\rbrace$.)

**Lemma A.1.** For $s\in S$ and $g,h\in W$:

1. $\delta_s(g)\le_B g$;
2. if $g\le_B h$, then $\delta_s(g)\le_B\delta_s(h)$.

*Proof.* (1) If $\delta_s(g)=gs$ with $\ell(gs)<\ell(g)$, then $g=(gs)s$ is one step of the chain definition, so $gs\le_B g$. Otherwise $\delta_s(g)=g$.

(2) If $\ell(hs)>\ell(h)$, then $\delta_s(h)=h$ and $\delta_s(g)\le_B g\le_B h$ by (1). Let $\ell(hs)<\ell(h)$, so $\delta_s(h)=hs$. If $\ell(gs)<\ell(g)$, then $\delta_s(g)=gs$ and (R1) gives $gs\le_B hs$. If $\ell(gs)>\ell(g)$, then $\delta_s(g)=g$ and (R2) gives $g\le_B hs$. ∎

Smallest case: $W=S_2$, $s=s_1$: $\delta_s(e)=e$, $\delta_s(s_1)=e$, and $e\le_B s_1$ gives $e\le_B e$.

### A.2 Sweeps

Fix a word $s_1s_2\cdots s_m$ in the generators. For a set of indices $I=\lbrace i_1<\dots<i_q\rbrace\subseteq\lbrace1,\dots,m\rbrace$ and $g\in W$, the **sweep** of $g$ along $I$ is

```math
\Delta_I(g)=\delta_{s_{i_1}}\bigl(\delta_{s_{i_2}}(\cdots\delta_{s_{i_q}}(g)\cdots)\bigr),
```

that is, the operators are applied in decreasing order of index, the one with the largest index first; $\Delta_\emptyset(g)=g$.

**Lemma A.2.**

1. Each $\Delta_I$ is Bruhat-monotone, and $\Delta_I(g)\le_B g$.
2. If $I\subseteq J$, then $\Delta_J(g)\le_B\Delta_I(g)$ for every $g$.

*Proof.* (1) A composition of monotone maps is monotone, and a composition of maps with $\delta_s(g)\le_B g$ has the same property (Lemma A.1 and transitivity).

(2) By induction on $\lvert J\setminus I\rvert$ it suffices to treat $J=I\cup\lbrace j\rbrace$ with $j\notin I$. Let $I_{>}$ and $I_{<}$ be the elements of $I$ larger and smaller than $j$. Then $\Delta_J=\Delta_{I_<}\circ\delta_{s_j}\circ\Delta_{I_>}$ and $\Delta_I=\Delta_{I_<}\circ\Delta_{I_>}$. Put $h=\Delta_{I_>}(g)$. By Lemma A.1(1), $\delta_{s_j}(h)\le_B h$, and by (1) applied to $\Delta_{I_<}$, $\Delta_{I_<}(\delta_{s_j}(h))\le_B\Delta_{I_<}(h)$, which is the claim. ∎

**Lemma A.2(3) (a sweep along a reduced suffix).** Let $w\le_L x$, and let $j_1<\dots<j_q$ be indices of the fixed word such that the subword $s_{j_1}\cdots s_{j_q}$ is a reduced word for $w$ (so $q=\ell(w)$). Then $\Delta_{\lbrace j_1,\dots,j_q\rbrace}(x)=xw^{-1}$, and each of the $q$ operators lowers the length by exactly one.

*Proof.* For $0\le p\le q$ let $w_p=s_{j_1}\cdots s_{j_p}$ ($w_0=e$, $w_q=w$) and let $g_p$ be the element obtained from $x$ after applying $\delta_{s_{j_q}},\dots,\delta_{s_{j_{p+1}}}$ (so $g_q=x$). We show by downward induction on $p$ that

```math
g_p=xw^{-1}w_p\qquad\text{and}\qquad \ell(g_p)=\ell(xw^{-1})+p .
```

For $p=q$: $g_q=x=xw^{-1}w$, and $\ell(x)=\ell(xw^{-1})+\ell(w)$ because $w\le_L x$. Let $1\le p\le q$ and assume the statement for $p$. Then $g_ps_{j_p}=xw^{-1}w_{p-1}$, so by (I1) $\ell(g_ps_{j_p})\le\ell(xw^{-1})+\ell(w_{p-1})\le\ell(xw^{-1})+p-1<\ell(g_p)$. Hence $s_{j_p}$ is a right descent of $g_p$, $g_{p-1}=\delta_{s_{j_p}}(g_p)=g_ps_{j_p}=xw^{-1}w_{p-1}$, and $\ell(g_{p-1})=\ell(g_p)-1=\ell(xw^{-1})+p-1$ by (I1). For $p=0$ this gives $g_0=xw^{-1}$. ∎

### A.3 The greatest elements of $\mathcal{LB}$ and $\mathcal{BR}$ (Gaetz–Gao, Lemma 3.1)

**Lemma A.3.** For all $x,y\in W$, $\mathcal{LB}(x,y)$ has a greatest element $\mathrm{lb}(x,y)$ and $\mathcal{BR}(x,y)$ has a greatest element $\mathrm{br}(x,y)$, and $\mathrm{br}(x,y)=\mathrm{lb}(y^{-1},x^{-1})^{-1}$.

*Proof.* Fix a reduced word $s_1\cdots s_m$ for $y$ ($m=\ell(y)$) and take sweeps along subsets of $\lbrace1,\dots,m\rbrace$ as in A.2.

*Construction.* Put $u_m=x$ and $u_{i-1}=\delta_{s_i}(u_i)$ for $i=m,m-1,\dots,1$, so that $u_0=\Delta_{\lbrace1,\dots,m\rbrace}(x)$. Let $i_1<\dots<i_r$ be the indices $i$ with $u_{i-1}\ne u_i$; at each of them $u_{i-1}=u_is_i$ and $\ell(u_{i-1})=\ell(u_i)-1$, and at every other index $u_{i-1}=u_i$. Since the operators act from the largest index down, $u_0=xs_{i_r}s_{i_{r-1}}\cdots s_{i_1}$. Put

```math
z=s_{i_1}s_{i_2}\cdots s_{i_r},\qquad\text{so that}\qquad u_0=xz^{-1},\quad \ell(u_0)=\ell(x)-r .
```

*$z\in\mathcal{LB}(x,y)$.* From $x=u_0z$ and (I1), $\ell(x)\le\ell(u_0)+\ell(z)\le(\ell(x)-r)+r=\ell(x)$, so both inequalities are equalities: $\ell(z)=r$, i.e. the word $s_{i_1}\cdots s_{i_r}$ is reduced, and $\ell(x)=\ell(xz^{-1})+\ell(z)$, i.e. $z\le_L x$. Since $s_{i_1}\cdots s_{i_r}$ is a reduced word for $z$ that is a subword of the reduced word $s_1\cdots s_m$ of $y$, the subword property (I4) gives $z\le_B y$.

*$z$ is greatest.* Let $w\in\mathcal{LB}(x,y)$. Since $w\le_B y$, the subword property (I4) gives indices $j_1<\dots<j_q$ such that $s_{j_1}\cdots s_{j_q}$ is a reduced word for $w$. Put $I=\lbrace j_1,\dots,j_q\rbrace$. Since $w\le_L x$, Lemma A.2(3) gives $\Delta_I(x)=xw^{-1}$, and Lemma A.2(2) with $I\subseteq\lbrace1,\dots,m\rbrace$ gives

```math
u_0=\Delta_{\lbrace1,\dots,m\rbrace}(x)\le_B\Delta_I(x)=xw^{-1}.
```

Both $u_0$ and $xw^{-1}$ are $\le_R x$: $\ell(x)=\ell(u_0)+\ell(u_0^{-1}x)$ because $u_0^{-1}x=z$ and $\ell(x)=\ell(u_0)+\ell(z)$; and $\ell(x)=\ell(xw^{-1})+\ell(w)$ with $(xw^{-1})^{-1}x=w$. So complementation (I7), applied with $(u,v,w)$ of Lemma 2.6 equal to $(u_0,xw^{-1},x)$, gives $(xw^{-1})^{-1}x\le_B u_0^{-1}x$, that is, $w\le_B z$. So $z$ is the greatest element $\mathrm{lb}(x,y)$ of $\mathcal{LB}(x,y)$.

*$\mathcal{BR}$.* By (I5), $g\le_B x$ iff $g^{-1}\le_B x^{-1}$, and $g\le_R y$ iff $g^{-1}\le_L y^{-1}$. So $g\mapsto g^{-1}$ maps $\mathcal{BR}(x,y)$ bijectively onto $\mathcal{LB}(y^{-1},x^{-1})$, and it preserves $\le_B$. Hence $\mathrm{lb}(y^{-1},x^{-1})^{-1}$ is the greatest element of $\mathcal{BR}(x,y)$. ∎

Smallest case: $W=S_3$, $x=y=s_1s_2=231$ (reduced word $s_1s_2$). Then $u_2=231$; $\delta_{s_2}(231)=213$ since $\ell(213)=1<2$; $\delta_{s_1}(213)=123$. So $r=2$, $z=s_1s_2=x$: indeed $x\le_L x$ and $x\le_B y$, so $\mathrm{lb}(x,x)=x$.

### A.4 The greatest element of $\mathcal{LR}$ (Gaetz–Gao, Theorem 1.1)

**Theorem A.4.** For all $x,y\in W$, $\mathcal{LR}(x,y)$ has a greatest element $\mathrm{lr}(x,y)$.

*Proof.* Define $z_0=x$ and, for $i\ge0$,

```math
z_{2i+1}=\mathrm{br}(z_{2i},y),\qquad z_{2i+2}=\mathrm{lb}(x,z_{2i+1}),
```

which exist by Lemma A.3.

*The sequence is non-increasing and becomes constant.* $z_{2i+1}\in\mathcal{BR}(z_{2i},y)$, so $z_{2i+1}\le_B z_{2i}$; $z_{2i+2}\in\mathcal{LB}(x,z_{2i+1})$, so $z_{2i+2}\le_B z_{2i+1}$. Hence $z_{j+1}\le_B z_j$ for every $j$, and by (I2) $\ell(z_{j+1})\le\ell(z_j)$, with equality only if $z_{j+1}=z_j$. A non-increasing sequence of non-negative integers is constant from some index $j_0$ on; so $z_{j+1}=z_j$ for all $j\ge j_0$, and the sequence is constant, equal to some $z$, from $j_0$ on. (Gaetz–Gao argue instead that $[e,x]_B$ is finite; the length argument avoids that fact.)

*$z\in\mathcal{LR}(x,y)$.* Choose $i$ with $2i\ge j_0$. Then $z=z_{2i+1}=\mathrm{br}(z_{2i},y)\in\mathcal{BR}(z,y)$, so $z\le_R y$; and $z=z_{2i+2}=\mathrm{lb}(x,z_{2i+1})\in\mathcal{LB}(x,z)$, so $z\le_L x$.

*$z$ is greatest.* Let $g\in\mathcal{LR}(x,y)$. We show $g\le_B z_j$ for all $j$ by induction on $j$. For $j=0$: $g\le_L x$, so $g\le_B x=z_0$ by (I6). If $g\le_B z_{2i}$, then, since $g\le_R y$, $g\in\mathcal{BR}(z_{2i},y)$, so $g\le_B\mathrm{br}(z_{2i},y)=z_{2i+1}$. If $g\le_B z_{2i+1}$, then, since $g\le_L x$, $g\in\mathcal{LB}(x,z_{2i+1})$, so $g\le_B\mathrm{lb}(x,z_{2i+1})=z_{2i+2}$. Hence $g\le_B z$. ∎

**Lemma A.4(2) (symmetry).** $g\in\mathcal{LR}(x,y)$ if and only if $g^{-1}\in\mathcal{LR}(y^{-1},x^{-1})$; hence $\mathrm{lr}(y^{-1},x^{-1})=\mathrm{lr}(x,y)^{-1}$.

*Proof.* By (I5), $g\le_L x$ iff $g^{-1}\le_R x^{-1}$, and $g\le_R y$ iff $g^{-1}\le_L y^{-1}$. So inversion maps $\mathcal{LR}(x,y)$ bijectively onto $\mathcal{LR}(y^{-1},x^{-1})$ and preserves $\le_B$, so it maps the greatest element to the greatest element. ∎

### A.5 (E1) (Gaetz–Gao, Corollary 1.2)

**Theorem A.5.** Let $w\in W\setminus\lbrace e\rbrace$ and $a=\mathrm{lr}(w^{-1},w)$. Then $a$ is an involution, $a\le_R w$, every involution $\iota\le_R w$ satisfies $\iota\le_B a$, and $a$ is the unique involution of maximal length in $\lbrace u : u\le_R w\rbrace$. In particular (E1)(1) and (E1)(2) hold, with $a(w)=a$.

*Proof.* Let $\mathcal{X}=\mathcal{LR}(w^{-1},w)=\lbrace g : g\le_L w^{-1},\ g\le_R w\rbrace$.

- *$a=a^{-1}$.* Lemma A.4(2) with $(x,y)=(w^{-1},w)$ gives $(y^{-1},x^{-1})=(w^{-1},w)$, so $\mathrm{lr}(w^{-1},w)=\mathrm{lr}(w^{-1},w)^{-1}$.
- *$a\ne e$.* Let $s$ be the first letter of a reduced word of $w$. Then $s\le_R w$ (§1.2, prefix), and by (I5) $s=s^{-1}\le_L w^{-1}$. So $s\in\mathcal{X}$ and $s\le_B a$ by Theorem A.4; by (I2), $\ell(a)\ge\ell(s)=1$.
- So $a^2=e$ and $a\ne e$: $a$ is an involution. Also $a\in\mathcal{X}$, so $a\le_R w$; thus $a\in\mathcal{IP}(w)$.
- *Domination.* Let $\iota$ be an involution with $\iota\le_R w$. By (I5), $\iota=\iota^{-1}\le_L w^{-1}$. So $\iota\in\mathcal{X}$, and $\iota\le_B a$ by Theorem A.4. This is (E1)(2).
- *Uniqueness.* For such $\iota$, (I2) gives $\ell(\iota)\le\ell(a)$, with equality only if $\iota=a$. So $a$ is the unique involution of maximal length that is $\le_R w$. This is (E1)(1). ∎

Smallest case: $W=S_3$, $w=321$. Every element is $\le_R w$ and $\le_L w^{-1}=w$, so $\mathcal{X}=S_3$ and $\mathrm{lr}(w,w)=321=w$, an involution, Bruhat-above the involutions $213$, $132$; this is the row $321$ of §3.3. For $w=231$: $\mathcal{X}=\lbrace g : g\le_L 312,\ g\le_R 231\rbrace=\lbrace e,213\rbrace$ ($312=s_2s_1$ and $231=s_1s_2$ each have only one reduced word, so the suffixes of $312$ are $e,s_1,312$ and the prefixes of $231$ are $e,s_1,231$), so $a=213$, the row $231$ of §3.3.

**Where each input is used.** (I3) in A.1; (I4) twice in A.3 (once each way); (I7), hence (I3) again through Lemma 2.6, in A.3; (I5) in A.3, A.4(2), A.5; (I6) in A.4; (I2) in A.4 and A.5. No other fact is used. In particular nothing in this appendix depends on $W$ being finite or of a particular type, and it is applied in this file with $W=S_N$ (§4) and $W=S_{2n}$ (§5).
