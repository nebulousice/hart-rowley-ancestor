# Hart–Rowley Conjecture 1.3 fails in $H_3$, $F_4$, $H_4$ and $E_6$: the full argument

*Write-up of 2026-11-27. **Expert verification required.** Verified so far with: Lean 4 (the H3 counterexample); two independent exact computer certificates; AI review of the type A proof (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 2: all three accept with zero issues). The type A theorem is not formalized in Lean. The counterexample argument in this document has not had a referee review.*

**Priority.** This is not the first refutation. S. Levinstein refuted Conjecture 1.3 first, in public, on 2026-11-27 ([Zenodo 10.5281/zenodo.24581942](https://doi.org/10.5281/zenodo.24581942); [s66lem/conjectures PR #15](https://github.com/s66lem/conjectures/pull/15)), with counterexamples in $H_3$, $F_4$ and $E_6$. This document is an independent verification. See [`LITERATURE_CHECK.md`](LITERATURE_CHECK.md).

## How to read the status labels

Every step below carries one status word.

| Word | Meaning here |
| --- | --- |
| **formalized** | Checked by the Lean 4 proof assistant's kernel: the build passes, with no `sorry` and only the standard axioms. |
| **certified** | An exact computation (integer arithmetic only), reproduced by an independent second program, with the output archived. |
| **numerical** | An exact computation by one program, or finite evidence, not independently reproduced. |
| **heuristic** | A plausible explanation, not a proof. |
| **proposed** | A written argument that has not yet passed independent review. |
| **open** | None of the above. |

Nothing in this document has the status "proved": that word is reserved for arguments that have passed independent review, and the arguments in this document have not been reviewed. (The type $A$ theorem recorded in Section 9 is proved; its proof was reviewed separately.)

## 0. What this document shows, and how

Hart and Rowley define a greedy way to write an element $w$ of a Coxeter group as a product of involutions. Their Conjecture 1.3 says that, in a finite Coxeter group, the number of factors never exceeds the rank.

We show that this fails.
- **Section 3** gives one explicit element of the group $H_3$ (rank 3) that needs 4 factors. The argument is short and done by hand. The same fact is also checked by two independent programs and by Lean.
- **Section 5** gives the exact maximum number of factors in $H_3$, $F_4$, $H_4$ and $E_6$, from exact computer enumeration.
- **Section 6** records numerical results for $E_7$ and $E_8$.
- **Section 9** records related results: the type $A$ theorem and its type $B$ consequence (proved, by a separately reviewed argument), and work in progress.

Sections 1, 2 and 4 set up the definitions and the standard facts used.

## 1. Definitions

### 1.1 Coxeter systems

A **Coxeter system** $(W,S)$ is a group $W$ with a finite generating set $S$ and a presentation

```math
W=\big\langle\ S \;\big|\; (st)^{m(s,t)}=e \text{ for all } s,t\in S\ \big\rangle ,
```

where $e$ is the identity, $m(s,s)=1$, and $m(s,t)=m(t,s)\in\{2,3,4,\dots\}\cup\{\infty\}$ for $s\ne t$ (no relation when $m=\infty$). In particular every $s\in S$ satisfies $s^2=e$. The **rank** of $W$ is $|S|$, the number of generators.

The **Coxeter graph** has vertex set $S$. Two vertices $s,t$ are joined by an edge when $m(s,t)\ge3$, and the edge carries the label $m(s,t)$ when $m(s,t)\ge4$. Unjoined generators commute.

### 1.2 Words, length, reduced words

A **word** is a finite sequence of generators. It **represents** the product of its letters. The **length** $\ell(w)$ of $w\in W$ is the least number of letters in a word representing $w$. A word of that least length is a **reduced word** for $w$. We write $\mathrm{Red}(w)$ for the set of reduced words of $w$.

Two facts follow at once from the definition.
- Any initial segment, and any final segment, of a reduced word is reduced. (If a segment could be shortened, so could the whole word.)
- Reversing a word for $w$ gives a word for $w^{-1}$, since every generator is its own inverse. So $\ell(w^{-1})=\ell(w)$, and $\mathrm{Red}(w^{-1})$ consists of the reversals of the words in $\mathrm{Red}(w)$.

### 1.3 Prefixes

Following Hart and Rowley [HR, §1], $u$ is a **prefix** of $w$ if some reduced word of $w$ begins with a reduced word of $u$.

**Lemma 1** (length criterion; status **proposed**, elementary). $u$ is a prefix of $w$ if and only if $\ell(u)+\ell(u^{-1}w)=\ell(w)$.

*Proof.* Suppose $\ell(u)+\ell(u^{-1}w)=\ell(w)$. Concatenate a reduced word of $u$ with a reduced word of $u^{-1}w$. The result represents $w$ and has $\ell(w)$ letters, so it is a reduced word of $w$ that begins with a reduced word of $u$.

Conversely, suppose a reduced word of $w$ splits as $ab$, where $a$ is a reduced word of $u$. Then $b$ represents $u^{-1}w$, so $\ell(u^{-1}w)\le\ell(w)-\ell(u)$. The reverse inequality $\ell(w)\le\ell(u)+\ell(u^{-1}w)$ holds because $w=u\cdot(u^{-1}w)$. $\square$

The relation "$u$ is a prefix of $w$" is the right weak order, written $u\le_R w$ in [BB, §3.1] and by Gaetz–Gao [GG]. Hart and Rowley call the same relation the weak (left) Bruhat order; only the name differs.

### 1.4 Involutions, ancestors, involution length

An **involution** is an element $u$ with $u\ne e$ and $u^2=e$.

For $w\ne e$, let $\mathcal{IP}(w)$ be the set of prefixes of $w$ that are involutions. It is not empty: the first letter of any reduced word of $w$ is a generator, hence an involution, and a prefix. Let $A(w)$ be the set of elements of $\mathcal{IP}(w)$ of largest length. Hart and Rowley say $w$ has the **ancestor property** if $A(w)$ has exactly one element. That element is the **ancestor** $a(w)$.

The **remainder** of $w$ is $\mathrm{rem}(w)=a(w)^{-1}w$. Because $a(w)$ is a prefix, Lemma 1 gives

```math
\ell(\mathrm{rem}(w))=\ell(w)-\ell(a(w))<\ell(w).
```

The **ancestor decomposition** of $w$ is $w=a_1a_2\cdots a_k$, where $a_1=a(w)$ and $a_2\cdots a_k$ is the ancestor decomposition of $\mathrm{rem}(w)$. The **involution length** is the number of factors:

```math
\ell_{\mathrm{inv}}(e)=0,\qquad \ell_{\mathrm{inv}}(w)=1+\ell_{\mathrm{inv}}(\mathrm{rem}(w))\quad(w\ne e).
```

The recursion stops because the length drops at every step. It is defined as long as every element met along the way has the ancestor property. This is the source's definition [HR, §1]: there, $w=a_1\cdots a_k$ with $a_i=a(a_i\cdots a_k)$ for $i<k$.

### 1.5 The conjectures

Paraphrasing [HR] (the wording of Conjecture 1.3 is the same in arXiv v1 and in the journal version):

- **Conjecture 1.2** (the Ancestor Conjecture). Every non-identity element of a finite Coxeter group has the ancestor property.
- **Conjecture 1.3.** Assume Conjecture 1.2. If $W$ is a finite Coxeter group, then $\ell_{\mathrm{inv}}(w)\le\mathrm{rank}(W)$ for every $w\in W$.

Gaetz and Gao proved Conjecture 1.2 for every Coxeter group, finite or infinite [GG, Cor. 1.2]. So the hypothesis of Conjecture 1.3 holds. Our counterexamples do not depend on their theorem: in each case we check the ancestor property directly (Sections 3 and 5).

### 1.6 Labelling of the groups

Generators are written $r_1,r_2,\dots$. We use these Coxeter graphs.
- $H_3$: $r_1 \overset{5}{-} r_2 - r_3$, that is, $m(r_1,r_2)=5$, $m(r_2,r_3)=3$, $m(r_1,r_3)=2$.
- $H_4$: $r_1 \overset{5}{-} r_2 - r_3 - r_4$.
- $F_4$: $r_1 - r_2 \overset{4}{-} r_3 - r_4$.
- $E_n$ ($n=6,7,8$), Bourbaki's labelling: the path $r_1 - r_3 - r_4 - \cdots - r_n$, with $r_2$ joined to $r_4$.

## 2. Two standard facts about reduced words

A **braid move** replaces, inside a word, a factor $sts\cdots$ of $m(s,t)$ alternating letters by the factor $tst\cdots$ of the same length. Here $s\ne t$ and $m(s,t)<\infty$. A braid move does not change the element represented, nor the number of letters.

In $H_3$ with our labelling there are exactly three kinds of braid move:

```math
r_1r_3\leftrightarrow r_3r_1,\qquad r_2r_3r_2\leftrightarrow r_3r_2r_3,\qquad r_1r_2r_1r_2r_1\leftrightarrow r_2r_1r_2r_1r_2 .
```

**Fact A** (Matsumoto's theorem; [BB, Thm. 3.3.1(ii)]). Any two reduced words of the same element are connected by a sequence of braid moves.

**Fact B** (Tits; [BB, Thm. 3.3.1(i)]). Any word can be turned into a reduced word for the same element by a sequence of braid moves and deletions of factors $ss$.

**Lemma 2** (reducedness test; status **proposed**, follows from Fact B). A word is reduced if and only if no word obtained from it by braid moves contains two equal adjacent letters.

*Proof.* Suppose the word is reduced. A braid-equivalent word with a factor $ss$ would give, after deleting $ss$, a shorter word for the same element. That is impossible.

Suppose the word is not reduced. By Fact B, some sequence of braid moves and deletions shortens it. Braid moves keep the number of letters, so at least one deletion occurs. The word just before the first deletion is obtained by braid moves alone, and it contains a factor $ss$. $\square$

**Lemma 3** (reduced words and involutions; status **proposed**, follows from Fact A). Let $a$ be a reduced word of $u\ne e$.
1. $\mathrm{Red}(u)$ is exactly the set of words obtained from $a$ by braid moves.
2. $u$ is an involution if and only if the reversal of $a$ lies in $\mathrm{Red}(u)$.

*Proof.* (1) Braid moves turn reduced words of $u$ into reduced words of $u$. Fact A gives the converse.

(2) The reversal of $a$ is a reduced word of $u^{-1}$ (Section 1.2). If $u=u^{-1}$, it is a reduced word of $u$. Conversely, if the reversal lies in $\mathrm{Red}(u)$, it represents both $u$ and $u^{-1}$, so $u=u^{-1}$. Then $u^2=e$, and $u\ne e$, so $u$ is an involution. $\square$

By definition, the prefixes of $w$ are the elements represented by the initial segments of the words in $\mathrm{Red}(w)$.

## 3. The counterexample in $H_3$

**Theorem 1.** In the Coxeter group $H_3$, the element

```math
w=r_2r_1r_2r_3r_1r_2r_1
```

has length 7 and ancestor decomposition

```math
w=(r_2r_1r_2)(r_1r_3)(r_2)(r_1).
```

Every element met in the recursion has the ancestor property. Hence $\ell_{\mathrm{inv}}(w)=4>3=\mathrm{rank}(H_3)$, and Conjecture 1.3 is false.

Status:
- **proposed** for the hand proof below;
- **certified** by two independent exact programs (Section 5);
- **formalized** in Lean 4, for a concrete permutation model of $H_3$ of order 120. The identification of that model with $H_3$ uses one cited fact, $|W(H_3)|=120$ (Section 7 and [`LEAN.md`](LEAN.md)).

*Plan of the proof.* We find all reduced words of $w$ with Lemma 2 and Fact A. This gives all prefixes of $w$. Lemma 3 tells us which prefixes are involutions. That gives the ancestor. We repeat the same steps on the remainder, three more times.

**Observation used throughout.** In every word that appears below, the only braid moves available are swaps of an adjacent pair $r_1,r_3$. Indeed, no such word contains $r_2r_3r_2$ or $r_3r_2r_3$, and the longest factor alternating between $r_1$ and $r_2$ has 4 letters ($r_2r_1r_2r_1$), fewer than $m(r_1,r_2)=5$. One checks this word by word in the tables below; every word in every listed class is shown.

**Step 1: $\ell(w)=7$, and $w$ has exactly two reduced words.** The word $r_2r_1r_2r_3r_1r_2r_1$ has one adjacent pair $r_3r_1$. Swapping it gives $r_2r_1r_2r_1r_3r_2r_1$. In that word the only adjacent $r_1,r_3$ pair swaps back. So the braid class of the word is

```math
\{\ r_2r_1r_2r_3r_1r_2r_1,\ \ r_2r_1r_2r_1r_3r_2r_1\ \}.
```

Neither word has two equal adjacent letters. By Lemma 2, both are reduced, so $\ell(w)=7$. By Lemma 3(1), they form $\mathrm{Red}(w)$.

**Step 2: the prefixes of $w$ and the ancestor $a(w)$.** The initial segments of the two reduced words give 9 distinct elements. For each one, the table lists its reduced words (the braid class of its segment, by Lemma 3(1)), the reversal of one of them, and the verdict of Lemma 3(2).

| Prefix $u$ | $\ell(u)$ | $\mathrm{Red}(u)$ | Reversal | Involution? |
| --- | --- | --- | --- | --- |
| $e$ | 0 | empty word | — | no ($u=e$) |
| $r_2$ | 1 | $r_2$ | $r_2$ | **yes** |
| $r_2r_1$ | 2 | $r_2r_1$ | $r_1r_2$ | no |
| $r_2r_1r_2$ | 3 | $r_2r_1r_2$ | $r_2r_1r_2$ | **yes** |
| $r_2r_1r_2r_1$ | 4 | $r_2r_1r_2r_1$ | $r_1r_2r_1r_2$ | no |
| $r_2r_1r_2r_3$ | 4 | $r_2r_1r_2r_3$ | $r_3r_2r_1r_2$ | no |
| $r_2r_1r_2r_1r_3$ | 5 | $r_2r_1r_2r_1r_3$, $r_2r_1r_2r_3r_1$ | $r_3r_1r_2r_1r_2$ | no |
| $r_2r_1r_2r_1r_3r_2$ | 6 | $r_2r_1r_2r_1r_3r_2$, $r_2r_1r_2r_3r_1r_2$ | $r_2r_3r_1r_2r_1r_2$ | no |
| $w$ | 7 | $r_2r_1r_2r_1r_3r_2r_1$, $r_2r_1r_2r_3r_1r_2r_1$ | $r_1r_2r_3r_1r_2r_1r_2$ | no |

In each "no" row with $u\ne e$, the reversal is not among the listed reduced words. So $\mathcal{IP}(w)=\{r_2,\ r_2r_1r_2\}$. The longest one, $r_2r_1r_2$, is unique. Hence $w$ has the ancestor property, with

```math
a(w)=r_2r_1r_2,\qquad \mathrm{rem}(w)=r_1r_3r_2r_1 .
```

The remainder is the final segment of the reduced word $r_2r_1r_2\cdot r_1r_3r_2r_1$, so it is reduced, of length 4.

**Step 3: the remainder $w_1=r_1r_3r_2r_1$.** Its braid class is $\{r_1r_3r_2r_1,\ r_3r_1r_2r_1\}$, so these are its reduced words. Its prefixes are the following 6 elements.

| Prefix $u$ | $\ell(u)$ | $\mathrm{Red}(u)$ | Reversal | Involution? |
| --- | --- | --- | --- | --- |
| $e$ | 0 | empty word | — | no |
| $r_1$ | 1 | $r_1$ | $r_1$ | **yes** |
| $r_3$ | 1 | $r_3$ | $r_3$ | **yes** |
| $r_1r_3$ | 2 | $r_1r_3$, $r_3r_1$ | $r_3r_1$ | **yes** |
| $r_1r_3r_2$ | 3 | $r_1r_3r_2$, $r_3r_1r_2$ | $r_2r_3r_1$ | no |
| $w_1$ | 4 | $r_1r_3r_2r_1$, $r_3r_1r_2r_1$ | $r_1r_2r_3r_1$ | no |

So $\mathcal{IP}(w_1)=\{r_1,\ r_3,\ r_1r_3\}$. The longest, $r_1r_3$, is unique. Hence $a(w_1)=r_1r_3$ and $\mathrm{rem}(w_1)=r_2r_1$.

**Step 4: $w_2=r_2r_1$.** Its only reduced word is $r_2r_1$. Its prefixes are $e$, $r_2$ and $r_2r_1$. Only $r_2$ is an involution (Step 2 shows $r_2r_1$ is not). So $a(w_2)=r_2$ and $\mathrm{rem}(w_2)=r_1$.

**Step 5: $w_3=r_1$.** It is an involution, so $a(w_3)=r_1$ and $\mathrm{rem}(w_3)=e$.

**Conclusion.** The recursion produced four ancestors, each unique:

```math
w=(r_2r_1r_2)(r_1r_3)(r_2)(r_1),\qquad \ell_{\mathrm{inv}}(w)=4>3. \qquad\square
```

The proof used only the Coxeter relations of $H_3$, Facts A and B, and Lemmas 1–3. It did not use Conjecture 1.2 or any computer output.

## 4. Parabolic subgroups

For $J\subseteq S$, the **standard parabolic subgroup** $W_J$ is the subgroup generated by $J$. The next lemma lets a counterexample in a smaller group pass to a larger one.

**Fact C** ([Hum, §5.5]). $(W_J,J)$ is a Coxeter system, and its length function is the restriction of $\ell$ to $W_J$.

**Lemma 4** (status **proposed**; standard). Let $w\in W_J$.
1. Every reduced word of $w$ in $(W,S)$ uses only letters from $J$.
2. The prefixes of $w$ in $(W,S)$ are the same elements as its prefixes in $(W_J,J)$, with the same lengths.
3. The ancestor decomposition of $w$, and $\ell_{\mathrm{inv}}(w)$, are the same in $W$ and in $W_J$.

Consequently, writing $M(W)=\max_{w\in W}\ell_{\mathrm{inv}}(w)$ when every non-identity element of $W$ has the ancestor property, we have $M(W)\ge M(W_J)$ whenever both sides are defined.

*Proof.* (1) By Fact C, a reduced word of $w$ in $(W_J,J)$ has $\ell(w)$ letters. So it is also a reduced word of $w$ in $W$, and it uses only letters from $J$. A braid move replaces an alternating factor in two letters $s,t$, which contains both $s$ and $t$, by another such factor. So braid moves do not change the set of letters used. By Fact A, all reduced words of $w$ in $W$ use the same set of letters, which lies in $J$.

(2) By (1) and Fact C, the reduced words of $w$ in $W$ and in $W_J$ are the same words. Prefixes are represented by their initial segments. They lie in $W_J$, where the two length functions agree.

(3) An element of $W_J$ is an involution in $W_J$ exactly when it is one in $W$. By (2), $\mathcal{IP}(w)$, $A(w)$ and $a(w)$ are the same in both groups. The remainder $a(w)^{-1}w$ lies in $W_J$. Induction on $\ell(w)$ finishes the proof. $\square$

Two consequences:
- $H_3$ is the parabolic subgroup of $H_4$ on $\{r_1,r_2,r_3\}$. So Theorem 1 gives only $M(H_4)\ge4=\mathrm{rank}(H_4)$. That is not a violation. The violation in $H_4$ (Section 5) needs its own computation.
- $E_7$ is the parabolic subgroup of $E_8$ on $\{r_1,\dots,r_7\}$. Levinstein used the analogous inclusion $E_6\subset E_7$ to get $M(E_7)\ge8$.

## 5. Exact computer verification: maxima in $H_3$, $F_4$, $H_4$ and $E_6$

**Proposition 5** (status **certified**). The following table holds. In each group, every element other than $e$ has the ancestor property, so $\ell_{\mathrm{inv}}$ is defined on the whole group.

| Group | Rank | Order | Involutions | $M(W)=\max\ell_{\mathrm{inv}}$ | Elements with $\ell_{\mathrm{inv}}>$ rank | Number of elements with $\ell_{\mathrm{inv}}=0,1,2,\dots$ |
| --- | --- | --- | --- | --- | --- | --- |
| $H_3$ | 3 | 120 | 31 | **4** | 3 | 1, 31, 52, 33, 3 |
| $F_4$ | 4 | 1152 | 139 | **5** | 50 | 1, 139, 374, 396, 192, 50 |
| $H_4$ | 4 | 14400 | 571 | **7** | 1746 | 1, 571, 2854, 5032, 4196, 1567, 173, 6 |
| $E_6$ | 6 | 51840 | 891 | **8** | 160 | 1, 891, 6664, 15236, 16850, 9622, 2416, 150, 10 |

Sample witnesses, as printed in the archived outputs:
- $H_3$: the three elements above the rank are $r_2r_1r_2r_1r_3r_2r_1$ (this is $w$ of Theorem 1), $r_2r_1r_3r_2r_1r_2r_3$ and $r_1r_2r_1r_2r_1r_3r_2r_1$, with decompositions $(r_2r_1r_2)(r_1r_3)(r_2)(r_1)$, $(r_2r_1r_3r_2)(r_1)(r_2)(r_3)$ and $(r_1r_2r_1r_2r_1)(r_3)(r_2)(r_1)$.
- $F_4$: $r_3r_4r_1r_2r_3r_2r_1$ has length 7 and decomposition $(r_1r_3)(r_2r_4)(r_3)(r_2)(r_1)$, so $\ell_{\mathrm{inv}}=5>4$.
- $H_4$: a shortest element above the rank (length 8) is $r_2r_1r_3r_2r_1r_2r_3r_4$.
- $E_6$: two shortest elements above the rank (length 11) are $r_1r_4r_3r_5r_4r_2r_6r_5r_4r_3r_1$ and $r_4r_3r_1r_6r_5r_4r_2r_3r_4r_5r_6$.

**How it was computed.** Two programs were written independently. They share no code with each other or with the exploratory scripts.

| | Primary program `cert/hr_tables.py` | Cross-check `cert/hr_tables_x.py` |
| --- | --- | --- |
| Model of the group | permutations of the root system, built exactly over $\mathbb{Z}[\varphi]$, $\varphi^2=\varphi+1$ | matrices of the Tits geometric representation, exact over $\mathbb{Z}[\varphi]$ ($H$ types) or $\mathbb{Z}[\sqrt2]$ ($F_4$); no roots |
| Length | number of positive roots sent to negative roots [Hum, §5.6] | breadth-first distance from $e$ in the Cayley graph (the definition) |
| Prefix test | inversion-set containment [BB, Prop. 3.1.3], self-tested against Lemma 1 on all 14,400 pairs of $H_3$ and all 1,327,104 pairs of $F_4$: 0 mismatches | walk down from $w$ by length-decreasing right multiplications by generators (the definition) |
| Involution test | the permutation squared is the identity | the matrix squared is the identity; cross-checked against word reversal: 0 mismatches |
| Ancestor | brute force over all involutions; uniqueness checked for every $w\ne e$ | first level below $w$ that contains an involution; uniqueness checked |

Both programs use integer arithmetic only. Neither uses floating point.

Coverage:
- $H_3$, $F_4$, $H_4$: both programs agree on every entry of the table.
- $E_6$: the primary program, together with the earlier exploratory implementation (`scripts/coxlib.py`, run through `scripts/rows.py`), agree on every entry. That exploratory implementation also models the group by root permutations, but it was written separately. The definition-level cross-check was too slow for $E_6$ and was stopped.

**Calibration against the source.** Both programs reproduce the source's worked example in $A_6$ (path graph $r_1-r_2-\cdots-r_6$), $w=r_6r_3r_2r_1r_4r_5$: the ancestor decomposition is $(r_3r_6)(r_2r_4)(r_1r_5)$, and the suffix version described in the source is $(r_3)(r_2r_4r_6)(r_1r_5)$. Swapping prefixes and suffixes would swap these two outputs. So the example fixes the side convention, and our reading matches the source's. Both programs also confirm the source's description of the ancestor of a **Coxeter element** (a product of all the generators, each used once, in some order): it is the product of the left descents of $w$ (the generators $s$ with $\ell(sw)<\ell(w)$). There were 0 failures on all 4 Coxeter elements of $H_3$ and all 8 of $F_4$ in both programs, and on all 8 of $H_4$ and all 32 of $E_6$ in the primary program.

**Other readings of the definition.** The certificates also compute three other greedy variants: the suffix version (peel the longest involution suffix on the right), and two alternating versions (prefix first, or suffix first). All three also exceed the rank. Their maxima are 4 in $H_3$ and 5 in $F_4$. So the failure does not come from a choice of side.

## 6. Numerical results: $E_7$ and $E_8$

These results are **numerical**: one whole-group program plus a single-element checker, not an independent second whole-group program.

- $E_7$ (order 2,903,040; rank 7). The exploratory program `scripts/big.py` finds $M(E_7)=9$, with 10,997 elements above the rank. The distribution of $\ell_{\mathrm{inv}}=0,1,\dots,9$ is 1, 10207, 139442, 531987, 919660, 802863, 385355, 102528, 9566, 1431. The ancestor is unique for every element.
- A witness. The element

```math
w_{7}=r_7r_5r_6r_4r_5r_3r_4r_1r_3r_2r_4r_5r_6r_7
```

has length 14 and ancestor decomposition

```math
(r_7r_5)(r_6r_4)(r_5r_3)(r_4r_1)(r_3r_2)(r_4)(r_5)(r_6)(r_7),
```

so $\ell_{\mathrm{inv}}(w_7)=9$. The single-element checker `scripts/check_word.py`, which uses a different method, gives the same decomposition.
- $E_8$ (rank 8). The checker gives $\ell_{\mathrm{inv}}(w_7)=9$ when $w_7$ is read in $E_8$, as Lemma 4 predicts. So $M(E_8)\ge9>8$. The exact value of $M(E_8)$ was not computed.

## 7. Formal verification in Lean (for $H_3$)

Status **formalized**, with one cited fact. Files: [`lean/Pipeline/HartRowley/Main.lean`](lean/Pipeline/HartRowley/Main.lean) and [`Data.lean`](lean/Pipeline/HartRowley/Data.lean). Details are in [`LEAN.md`](LEAN.md).

The Lean files define $G$, the group of permutations of the 30 roots of $H_3$ generated by three explicit permutations $s_1,s_2,s_3$. The Lean kernel checks, by evaluation (`decide +kernel`):
- $s_1,s_2,s_3$ satisfy the Coxeter relations of $H_3$;
- $G$ has exactly 120 elements, and its multiplication table is consistent with the permutations;
- word lengths, the list of involutions, and the prefix test of Lemma 1;
- the witness: length 7, exactly two involution prefixes ($r_2$ and $r_2r_1r_2$), the decomposition $(r_2r_1r_2)(r_3r_1)(r_2)(r_1)$, and $\ell_{\mathrm{inv}}=4$;
- the whole group: every element other than $e$ has a unique ancestor, the maximum of $\ell_{\mathrm{inv}}$ is 4, and exactly 3 elements exceed 3.

**Cited, not formalized:** $|W(H_3)|=120$ [Hum, §2.11]. With the checked relations, $r_i\mapsto s_i$ defines a surjective homomorphism $W(H_3)\to G$. Equal orders make it an isomorphism. So word length in $G$ equals the Coxeter length in $H_3$.

## 8. The source's computational evidence

Hart and Rowley report that Conjectures 1.2 and 1.3 were verified with Magma for every irreducible finite Coxeter group of order below 100,000, which includes $H_3$, $H_4$, $F_4$ and $E_6$ [HR, §3]. Proposition 5 contradicts this report for exactly these four groups. We do not know the cause.

**A possible explanation** (status **heuristic**: a guess about code we have not seen). Consider instead the minimal number of involution factors in a length-additive factorization (Section 9). That non-greedy count satisfies the rank bound in $H_3$, $F_4$, $H_4$ and $E_6$. A computation of that count, rather than of the greedy $\ell_{\mathrm{inv}}$, would have reported the conjecture as verified in these groups.

## 9. Related results and work in progress

Each item carries its current status. Items 1 and 2 are proved; their proof is separate from the counterexample argument above and is summarised in [`README.md`](README.md). Item 3 is work in progress.

Notation. $S_N$ is the symmetric group on $\{1,\dots,N\}$, the Coxeter group of type $A_{N-1}$ (rank $N-1$, generators the adjacent transpositions). The **displacement** of a permutation is $\mathrm{disp}(w)=\max_i|w(i)-i|$. A factorization $w=u_1\cdots u_k$ is **length-additive** if $\ell(w)=\ell(u_1)+\cdots+\ell(u_k)$. Let $m(w)$ be the least $k$ such that $w$ has a length-additive factorization into $k$ involutions, with $m(e)=0$. Since the ancestor decomposition is one such factorization, $m(w)\le\ell_{\mathrm{inv}}(w)$.

1. **Type $A$** (status **proved**: AI review, round 2, three of three referees accept with zero issues; also **numerical**). For every $w\in S_N$, $\ell_{\mathrm{inv}}(w)\le\mathrm{disp}(w)\le N-1$, which is Hart–Rowley Conjecture 1.3 in type $A$. The proof uses Gaetz–Gao [GG, Cor. 1.2], an unrefereed preprint, but does not assume it: an appendix re-proves it in full from textbook facts [BB] (the lifting and subword properties of Bruhat order), following Gaetz–Gao's argument, which is credited to them. Exact checks of all of $S_N$ for $N\le7$ agree; exploratory runs reach $S_{10}$, where the key combinatorial lemma had 0 failures on 4,850,718 extreme positions. These checks are corroboration only; no proof step uses them. Not formalized in Lean.
2. **Type $B$** (status **proved** for the reduction, same review; **numerical**, exact, for $B_n$ with $n\le4$ and for $D_3$, $D_4$). The weak order on $B_n$ is the restriction of the weak order on $S_{2n}$, so ancestors and involution lengths agree. With item 1 this gives $\ell_{\mathrm{inv}}\le 2n-1$ in $B_n$, not the conjectured $n$. Hart–Rowley Conjecture 1.3 for type $B$ is **open**. The reduction fails for type $D$ already in $D_3$: the signed permutation with window $(-3,-1,2)$ has involution length 2 in $D_3$ and 3 in $S_6$; $D_3$ has 6 ancestor and 4 involution-length mismatches, $D_4$ has 55 and 33. Type $D$ is **open**.
3. **The minimal count $m$** (status **numerical**, exact). The maximum of $m$ equals the rank in $H_3$ (3), $F_4$ (4), $H_4$ (4), $E_6$ (6) and $E_7$ (7, by two methods), in $A_n$ for $n\le9$ and in $B_n$ for $n\le8$. It is 3, 5, 5 in $D_4$, $D_5$, $D_6$, and 2 in the dihedral groups $I_2(k)$ for $3\le k\le20$. Whether $m(w)\le\mathrm{rank}(W)$ holds in every finite Coxeter group is **open**.

## References

- [BB] A. Björner and F. Brenti, *Combinatorics of Coxeter Groups*, Graduate Texts in Mathematics 231, Springer, 2005. [doi:10.1007/3-540-27596-7](https://doi.org/10.1007/3-540-27596-7)
- [GG] C. Gaetz and Y. Gao, *Separable elements and Bruhat orders*, preprint, 2026-10-16. [arXiv:2610.22857](https://arxiv.org/abs/2610.22857)
- [HR] S. B. Hart and P. J. Rowley, *A note on involution prefixes in Coxeter groups*, [arXiv:2502.00777](https://arxiv.org/abs/2502.00777) (2025); Arch. Math. (Basel) 126(3) (2026) 231–237, [doi:10.1007/s00013-025-02220-z](https://doi.org/10.1007/s00013-025-02220-z)
- [Hum] J. E. Humphreys, *Reflection Groups and Coxeter Groups*, Cambridge Studies in Advanced Mathematics 29, Cambridge University Press, 1990. [doi:10.1017/CBO9780511623646](https://doi.org/10.1017/CBO9780511623646)
- [Lev] S. Levinstein, *Counterexamples to the Hart–Rowley Involution-Length Conjecture in Finite Coxeter Groups*, Zenodo preprint, 2026-11-27. [doi:10.5281/zenodo.24581942](https://doi.org/10.5281/zenodo.24581942); also [s66lem/conjectures PR #15](https://github.com/s66lem/conjectures/pull/15)
