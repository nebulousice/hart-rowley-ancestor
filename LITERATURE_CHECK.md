# Literature and priority check

- **Date:** 2026-11-27 (three searches: 06:00–12:45 UTC, and 13:33–14:20 UTC).
- **Question:** has anyone else refuted, proved, or computed anything about Hart–Rowley Conjecture 1.3, and which of our values are already public?

## Conclusion first

- **We are not first.** Salem Levinstein refuted Conjecture 1.3 in public on 2026-11-27 at 08:30 UTC ([Zenodo 10.5281/zenodo.24581942](https://doi.org/10.5281/zenodo.24581942); [s66lem/conjectures PR #15](https://github.com/s66lem/conjectures/pull/15), merged 08:32 UTC). Priority for the refutation is Levinstein's.
- Levinstein's report gives counterexamples in $H_3$, $F_4$ and $E_6$, the maxima 4 ($H_3$), 5 ($F_4$) and 8 ($E_6$), the bound $\ge 8$ in $E_7$ (from $E_6\subset E_7$), and a Lean 4 proof for $H_3$.
- Our computations were made **independently**, without knowledge of Levinstein's report, in a private repository. They are offered as an independent confirmation, not as a priority claim.
- Not in Levinstein's report, and not found in any other source we searched: the maximum 7 in $H_4$ (certified), the maximum 9 in $E_7$ (numerical), and the bound $\ge 9$ in $E_8$ (numerical). This is a search result, not a claim of discovery.
- **Conjecture 1.2 is not ours either.** Gaetz and Gao proved it for every Coxeter group ([arXiv:2610.22857](https://arxiv.org/abs/2610.22857), 2026-10-16).

## Timeline

| Date (UTC) | Event | Source |
| --- | --- | --- |
| 2025-02-02 | Hart and Rowley post Conjectures 1.2 and 1.3, and report a Magma verification for all irreducible finite Coxeter groups of order below 100,000, which includes $H_3$, $H_4$, $F_4$, $E_6$ | [arXiv:2502.00777v1](https://arxiv.org/abs/2502.00777v1) |
| 2026-02-07 | Journal version online. Conjecture 1.3 is unchanged (pp. 231–232 read from the publisher's free preview). No erratum is registered | [Arch. Math. 126(3) 231–237](https://doi.org/10.1007/s00013-025-02220-z) |
| 2026-06-19 | Gimenez proves the ancestor property for fully commutative elements in any Coxeter system | [arXiv:2606.21104](https://arxiv.org/abs/2606.21104) |
| 2026-10-15 | The trureturing project refutes only the source's tentative remark that prefix and suffix involution lengths agree ($S_4$, $w=2431$: 3 versus 2). It states that this does not settle Conjectures 1.2 or 1.3 | [issue #17426](https://github.com/the-omega-institute/trureturing/issues/17426), [PR #17548](https://github.com/the-omega-institute/trureturing/pull/17548) |
| 2026-10-16 | Gaetz and Gao prove Conjecture 1.2 for every Coxeter group. The paper does not mention Conjecture 1.3 | [arXiv:2610.22857](https://arxiv.org/abs/2610.22857) |
| 2026-11-27, about 06:00 | Our exploratory scripts find elements of $H_3$, $F_4$, $H_4$, $E_6$, $E_7$ with involution length above the rank (private repository, unpublished; commit `fd69e12`, 06:06 UTC) | this repository's history |
| 2026-11-27, 08:30 | **Levinstein's Zenodo preprint is created; PR #15 merged at 08:32** | [Zenodo](https://doi.org/10.5281/zenodo.24581942), [PR #15](https://github.com/s66lem/conjectures/pull/15) |
| 2026-11-27, 12:40 | Our priority search finds Levinstein's report | this file |
| 2026-11-27, 13:06–13:26 | Our certificates and the Lean check are run and archived | [`cert/`](cert/) |
| 2026-11-27, 13:33–14:20 | Second search, for the remaining values: no source reports $M(H_4)=7$, $M(E_7)=9$ or $M(E_8)\ge9$, or anything on types $A$, $B$, $D$ or on the minimal-factorization version | this file |

**On the 06:00 entry.** It is listed only to keep the record accurate. A private, unpublished computation establishes no priority, and we claim none. The first public refutation is Levinstein's.

## Sources searched

| Source | What was searched | Finding |
| --- | --- | --- |
| Primary source, all versions | arXiv abs page (only v1 exists) and TeX source; Crossref metadata; the journal's free first-pages preview; zbMATH review | Statements as in [`PROOF.md`](PROOF.md) §1.5. Conjecture 1.3 is unchanged in the journal version. The journal's §3 (behind a paywall) could not be read; zbMATH and the unchanged reference list indicate it still reports the Magma check |
| C. Davis, `1337isnot1337/solved-conjectures`: all issues, PRs and branches | Full listings grepped for Hart, Rowley, ancestor, involution prefix/length, Coxeter, the arXiv ids, Levinstein | No hit (false positives only: "Ehrhart", git ancestors) |
| `PerAlexandersson/RealRooted` and SymCat | Full issue and PR listings; the open-problems page | No hit |
| `the-omega-institute/trureturing` | Issue and PR searches; problem dossier | Only the suffix-remark refutation (row above) |
| arXiv | Term, author and citer searches (OpenAlex, Semantic Scholar where available) | Gaetz–Gao (Conjecture 1.2); Gimenez; Hart–Kelsey–Rowley (motivation only). No paper on Conjecture 1.3 |
| Web and GitHub code search | Several queries on the conjecture, counterexample, $H_3$/$F_4$/$E_8$ | Only Levinstein's PR #15 and the trureturing files. Web summaries still call the conjecture open (index lag) |
| OEIS, MathOverflow, Open Problem Garden, `google-deepmind/formal-conjectures` | Keyword searches; full tree grep | Nothing relevant |
| MathDB | Entries #366919, #366920 (Conjecture 1.3), #408774 | #366920 still shows no progress after the refutation (catalog lag) |
| Zenodo | Phrase searches at 06:17 UTC and again at 12:40 UTC | 06:17: nothing. 12:40: exactly one record, Levinstein's |
| `s66lem/conjectures` (Levinstein) | All PRs, branches, the result folder (README, PROOF, AUDIT, paper, `verify.py`), the Lean file | See the comparison below. Nothing on $H_4$, $E_8$, classical types or a corrected statement; its audit says the $H_4$ and $E_8$ maxima were dropped as unverified |
| Hexagon, Prove2Me, formalpedia, other AI-group repositories | Site and tree searches | No hit |
| eulersolve, UnsolvedMath, Crab Research | Search endpoints unavailable | **Inconclusive** |

## Comparison with Levinstein's report

| | Levinstein ([Zenodo](https://doi.org/10.5281/zenodo.24581942), [PR #15](https://github.com/s66lem/conjectures/pull/15)) | This repository |
| --- | --- | --- |
| Public | Yes, 2026-11-27 08:30 UTC | No: our work was in a private repository; this repository was made public later on 2026-11-27 |
| Groups with a violation | $H_3$, $F_4$, $E_6$; $E_7$ by inclusion | $H_3$, $F_4$, $H_4$, $E_6$ (certified); $E_7$, $E_8$ (numerical) |
| Maxima | $H_3$ 4, $F_4$ 5, $E_6$ 8; $E_7\ge8$ | $H_3$ 4, $F_4$ 5, $H_4$ 7, $E_6$ 8 (certified); $E_7$ 9, $E_8\ge9$ (numerical) |
| $H_3$ labelling | Mathlib's: the edge labelled 5 joins $r_2$ and $r_3$ | ours: the edge labelled 5 joins $r_1$ and $r_2$ |
| $H_3$ witnesses | $r_2r_3r_1r_2r_3r_2r_1$ and $r_2r_3r_2r_3r_1r_2r_3$ (Mathlib labelling), length 7 | $r_2r_1r_2r_3r_1r_2r_1$ (our labelling), length 7; all three elements above the rank are listed |
| Exact computation | Python over $\mathbb{Q}(\sqrt5)$ (`verify.py`) | Two independent Python programs, integer arithmetic only |
| Lean 4 | $H_3$, for every abstract Coxeter system of type $H_3$ in Mathlib (987 lines; reported, not re-checked by us) | $H_3$, in a concrete permutation model; the identification with $H_3$ rests on the cited order 120. **Levinstein's statement is stronger** |
| Other readings of the definition | — | Suffix and alternating greedy variants also fail in $H_3$ and $F_4$ |
| Review | One-pass AI audit, reported as PASS | Counterexample write-up: no referee review yet. Type $A$ proof: AI review, round 2, three of three referees accept with zero issues |
| Types $A$, $B$, $D$; minimal factorizations | Not addressed | Type $A$: Conjecture 1.3 holds (proved). Type $B$: reduced to type $A$, giving the bound $2n-1$, not the rank. Type $D$ and minimal factorizations: work in progress (see [`README.md`](README.md)) |

## Gaetz–Gao

Gaetz and Gao ([arXiv:2610.22857](https://arxiv.org/abs/2610.22857), *Separable elements and Bruhat orders*, 2026-10-16, unrefereed preprint) prove that, in any Coxeter group, every $w\ne e$ has a unique longest involution prefix $a(w)$, and that every involution prefix of $w$ lies below $a(w)$ in Bruhat order. This is Conjecture 1.2 of Hart–Rowley, for all Coxeter groups. The paper states that ChatGPT produced the proof of its Theorem 1.1, on which this corollary rests. Our counterexamples do not depend on it: the ancestor property is checked directly in every group we use. Our type $A$ proof uses their Corollary 1.2, but does not assume it: it re-proves it in full from textbook facts (Björner–Brenti), following their argument and crediting it to them.
