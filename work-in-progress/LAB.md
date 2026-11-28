VERDICT C1 (2026-11-27): (1) type A ALIVE: Lemma K plus the per-condition up-set argument give l_inv(w) <= disp(w) <= rank (proposed; uses Gaetz–Gao Cor. 1.2 and Lemma 3.1); exhaustive S_3..S_10 with two implementations, 0 failures (numerical). (2) types B/D ALIVE, formula open: B_n ancestors equal type-A ancestors in S_2n (proposed; numerical B3..B6), so l_inv^B <= 2n-1; least admissible statistic G has max G = rank for B2..B8 and max G <= rank for D4..D7 (numerical); no closed formula; D8 not finished. (3) R4 ALIVE: no violation; exact max m = rank in H3, F4, H4, E6, E7 (numerical, two methods for E7); E8 not computed.

# Plan C1 lab notes (hart-rowley-ancestor, stage 3, re-scoped option C)

Explorer run of 2026-11-27. Everything below is working notes. Status words as in CLAUDE.md. Nothing is committed. All computations are exact integer computations (permutations, inversion counts, bitmasks); none uses floating point except the sign test of positive roots over Z[phi] in `npcox.py`, which is exact in effect because roots have bounded integer coefficients (the test is only used to split roots into positive and negative, and the group orders reproduce).

Conventions (as in BRIEF.md 1.2): u <=_R w iff l(w) = l(u) + l(u^-1 w) (u is a prefix of w); a(w) the longest involution prefix; rem(w) = a(w)^-1 w; l_inv(w) = 1 + l_inv(rem w). Permutations act on the left of values: for w in S_N, s_i w swaps the VALUES i, i+1, and w s_i swaps the POSITIONS i, i+1. disp(w) = max_i |w(i) - i|.

## 1. Type A: the rank bound

### 1.1 Statements

**Theorem A** (status: proposed; depends on Gaetz–Gao arXiv:2610.22857, unrefereed). For every w in S_N, l_inv(w) <= disp(w) <= N - 1 = rank(A_{N-1}). Equality l_inv = N-1 holds for the Coxeter element s_1 s_2 ... s_{N-1}.

The proof uses three outside facts:
- (GG1) Gaetz–Gao Cor. 1.2: every involution iota <=_R w satisfies iota <=_B a(w).
- (GG2) Gaetz–Gao, first paragraph of the proof of Lemma 3.1: if u, v <=_R w and u <=_B v, then v^-1 w <=_B u^-1 w. (A five-line induction with the lifting property; read and checked by hand in this run.)
- (BB) Björner–Brenti Thm 2.1.5 (rank-matrix criterion): x <=_B y in S_N iff x[i,k] <= y[i,k] for all i, k, where x[i,k] = #{j <= i : x(j) >= k}. Plus the subword property: a suffix of w is Bruhat-below w.

and one new combinatorial lemma, which needs no Gaetz–Gao:

**Lemma K** (status: proposed, with complete proof below; checked exhaustively for all w in S_N, N <= 10: 4,850,718 extremal pairs in S_10 alone, 0 failures). Let w in S_N, D = disp(w) >= 1, and let p be a position with w(p) = p + D; put v = w(p). Then there is an involution iota with iota <=_R w and iota(v) != v. Explicitly: L = w({1..p}), R = {1..v} \ L; walk downward from c = v with h(c) = #(L ∩ [c,v]) - #(R ∩ [c,v]); let c0 be the largest c with h(c) = 0 and J = [c0, v]; iota swaps the i-th smallest element of L ∩ J with the i-th smallest element of R ∩ J, for every i, and fixes all other values.

### 1.2 Proof of Lemma K (proposed)

(a) Every value b > v sits at a position > p: if b sits at position j, then b - j <= D, so j >= b - D > v - D = p. Hence L = w({1..p}) ⊆ {1..v}, |L| = p, and R = {1..v} \ L (the values <= v sitting right of p) has |R| = v - p = D.

(b) Every r in R satisfies r >= p + 1 - D: r sits at a position j >= p + 1, and j - r <= D.

(c) The walk returns to 0. h(v) = 1 since v is in L, and h changes by ±1 per step. Let c* = max(1, p + 1 - D). By (b), all D elements of R lie in [c*, v - 1]. If c* = p + 1 - D, the interval [c*, v] has v - c* + 1 = 2D slots, D of them in R, so h(c*) = D - D = 0. If c* = 1 (that is, p <= D), then h(1) = |L| - |R| = p - D <= 0. Either way h reaches 0; c0 exists and h(c) >= 1 for c0 < c <= v.

(d) Ballot property. Let l_1 < ... < l_k be L ∩ J and r_1 < ... < r_k be R ∩ J (equal numbers since h(c0) = 0, and k >= 1 since v is in J). For c0 <= c < v, #R∩[c0,c] - #L∩[c0,c] = h(c+1) >= 1 (complement within J). At c = l_i this gives #R∩[c0,l_i] >= i + 1 > i (or = i when l_i = v = l_k), so r_i < l_i for every i. In particular iota(v) = iota(l_k) = r_k < v.

(e) Inversions of iota. Let a < b with iota(a) > iota(b). Values outside J are fixed, and iota maps J into J, so a, b are both in J (if a < c0 then iota(a) = a < c0 <= iota(b); if b > v then iota(b) = b > v >= iota(a)). Inside J: two elements of L are not inverted (l_i -> r_i is increasing), two elements of R are not inverted (r_i -> l_i is increasing), and a pair l_i < r_j is not inverted because r_i < l_i < r_j forces i < j and then iota(l_i) = r_i < l_i < l_j = iota(r_j). So every inversion of iota is a pair (r, l) with r < l, r in R, l in L.

(f) iota <=_R w. For any x, x <=_R w iff T_L(x) ⊆ T_L(w) (Björner–Brenti Prop. 3.1.3), and for permutations the transposition (a b), a < b, lies in T_L(x) iff b stands left of a in x, i.e. iff (a, b) is an inversion of x^-1 (as a pair of positions). For the involution iota this is the inversion set of iota. By (e) every such pair is (r, l) with r < l, l in L (a position <= p) and r in R (a position > p): so l stands left of r in w and (r, l) is in T_L(w). Hence iota <=_R w. QED.

(The proof never uses N; the walk construction is exactly what `typeA_K_vec.py` implements, and its checks (k1)–(k4) test each of (c), (f) [as the length identity l(iota) + l(iota w) = l(w)], (d) and (e).)

### 1.3 Proof of Theorem A from Lemma K (proposed; uses GG1, GG2, BB)

Write r = rem(w), a = a(w), D = disp(w) >= 1. We show disp(r) < D; then induction on D gives l_inv(w) = 1 + l_inv(r) <= 1 + disp(r) <= D.

1. disp is Bruhat-monotone. For D' >= 1 put C_{i}^{D'} = {x : exists j <= i with x(j) >= i + D'} = {x : x[i, i+D'] >= 1}; by (BB) it is an up-set. x has an excedance >= D' iff x lies in some C_i^{D'} (take i = j). Drops: x has a drop >= D' iff phi(x) has an excedance >= D', where phi(x) = w0 x w0, and phi preserves Bruhat order. So {disp >= D'} is an up-set.
2. r is a suffix of w (w = a r, length-additive), so r <=_B w and disp(r) <= D.
3. No new extremal excedance. Suppose r(p) = p + D. Then r lies in C_p^D. On {disp <= D}, C_p^D = {x : x(p) = p + D} (j <= p and x(j) >= p + D >= j + D force j = p). Since r <=_B w and C_p^D is an up-set, w lies in C_p^D, so w(p) = p + D =: v. Lemma K gives an involution iota <=_R w with iota(v) != v. By (GG1) iota <=_B a; both are prefixes of w, so by (GG2) r = a^-1 w <=_B iota^-1 w = iota w. As r is in the up-set C_p^D, so is iota w; iota w is a suffix of w, so disp(iota w) <= D, hence (iota w)(p) = p + D, i.e. iota(v) = v. Contradiction.
4. No extremal drop. phi is a diagram automorphism of (S_N, S), so it preserves length, both weak orders, Bruhat order and involutions; hence a(phi(w)) = phi(a(w)) and rem(phi(w)) = phi(rem(w)). Also disp(phi(x)) = disp(x), and x(p) = p - D iff phi(x)(N+1-p) = (N+1-p) + D. Step 3 applied to phi(w) excludes r(p) = p - D.
5. So |r(p) - p| <= D - 1 for all p: disp(r) < D. QED.

Remarks.
- The argument does **not** need a single z that lowers all extremal entries at once (Lemma A′ of BRIEF P2). Each extremal entry is handled by its own involution prefix, and the up-set property of each single condition C_p^D transfers the conclusion to a(w). This is the "per-condition" form of lead L5; it applies to any statistic of the form f(x) = max_k [x in U_k] with each U_k a union of up-sets, provided each up-set can be escaped by one involution prefix.
- Lemma K is elementary and could be formalised in Lean (permutations of Fin N, inversion sets); GG1/GG2 and the rank-matrix criterion are not in Mathlib (see BRIEF P2 Lean note).
- Sharpness: disp is the least Bruhat-monotone upper bound of l_inv for A3–A6 (BRIEF L4), so Theorem A is the strongest statement of its kind.

### 1.4 Type A computations

| Check | Range | Implementation | Result | Status |
| --- | --- | --- | --- | --- |
| disp(rem w) < disp(w), all w != e | S_3 .. S_10 (all 3,628,800 of S_10) | I1 `typeA_dp.py` (prefix DP over Lehmer ranks) | 0 failures | numerical (exact) |
| same | S_3 .. S_10 | I2 `typeA_gg.py` (Gaetz–Gao lr sweeps, per element, no tables) | 0 failures [S_10: see 1.5] | numerical (exact) |
| a(w) identical in I1 and I2 | S_3 .. S_9 (S_10: see 1.5) | `compare_A.py` | 0 mismatches of 362,880 (S_9) | numerical (exact), two implementations |
| l_inv <= disp; l_inv <= N-1 | S_3 .. S_10 | I1, I2 | 0 failures; max l_inv = N-1 | numerical (exact) |
| a(w) is an involution, a prefix; ancestor ties | S_3 .. S_10 | I1 (DP), I2 | 0, 0, 0 | numerical (exact) |
| Lemma K: (k1) walk returns, (k2) iota prefix, (k3) iota(v) != v, (k4) inversion shape | all extremal pairs, S_3 .. S_10 (excedances and, via phi, drops) | `typeA_K_vec.py` (vectorised) and `check_K` in `typeA_gg.py` (pure python, S_3..S_8) | 0 failures; pairs: 8, 38, 182, 1064, 7176, 56004, 493050, 4850718 | numerical (exact); both codes agree on the pair counts |
| iota_p <=_B a(w) (GG1 instance) | S_3 .. S_7 | `typeA_gg.py` (tableau criterion) | 0 failures | numerical (exact) |
| random / banded permutations, N = 20..60 | see 1.5 | I2 | see 1.5 | numerical (exact) |

l_inv distributions reproduce BRIEF 1.3 exactly (A2..A9).

### 1.5 Large exhaustive and random runs

- S_10 exhaustive, I2 (Gaetz–Gao sweeps, 2335 s): a not involution 0, a not prefix 0, disp not decreasing 0, l_inv > disp 0, max l_inv 9, distribution 1, 9495, 202154, 842763, 1138044, 719549, 384272, 200818, 95572, 36132 (identical to I1 and to BRIEF/big.py). `compare_A.py`: **0 mismatches between I1 and I2 on all 3,628,800 ancestors** (and on S_3..S_9).
- Proof replay (`typeA_replay.py 7`): for every extremal pair of S_3..S_7 (7176 in S_7), iota_p <=_B a(w) (GG1), rem(w) <=_B iota_p w (GG2) and the conclusions of step 3 hold: 0 failures.
- Random uniform permutations (I2 + Lemma K, `typeA_gg.py random`): N = 20 (3000), 30 (2000), 40 (1000), 60 (300): 0 failures of (D), (I), (K); max l_inv 19, 29, 37, 47 (= N-1 for N = 20, 30).
- Banded permutations (disp <= b, random adjacent-transposition walks): (N, b) = (20, 4), (40, 6), (60, 10): 0 failures; max l_inv = 4, 6, 10 = b, so the bound l_inv <= disp is attained far below the rank.

## 2. Types B and D

### 2.1 Type B reduces to type A (status: proposed proof; numerical check B2..B6)

View B_n inside S_{2n} as the signed permutations of -n < ... < -1 < 1 < ... < n (the fixed points of conjugation by the longest element of S_{2n}).

- l_A(w) = 2 l_B(w) - neg(w), neg(w) = #{i > 0 : w(i) < 0} (long roots of B_n correspond to pairs of transpositions, short roots to single ones). Checked for every w in B2..B6 (`b_weak_restrict.py`, 0 failures).
- Weak orders agree: for u, w in B_n, u <=_R w in B_n iff u <=_R w in S_{2n}. Proof sketch (proposed): l(uv) = l(u) + l(v) - 2|C| with C = I_R(u) ∩ I_L(v), and the same identity holds for the count of short roots, neg(uv) = neg(u) + neg(v) - 2|C_short|. So A-additivity 2 l_B(uv) - neg(uv) = (2 l_B(u) - neg(u)) + (2 l_B(v) - neg(v)) reads 4|C| = 2|C_short|, forcing C = ∅, i.e. B-additivity; the converse is immediate. Checked on all pairs of B2..B4 and 250,000 pairs of B5, B6 (biased towards prefixes), 0 mismatches.
- Hence a_B(w) = a_A(w): a_A(w) is unique, so it is fixed by the automorphism, so it lies in B_n and is an involution prefix in B_n; by GG1 in S_{2n} it is Bruhat-above every involution prefix, Bruhat order on B_n is induced from S_{2n} (Björner–Brenti Cor. 8.1.9), and Bruhat order is graded by l_B, so it is the longest B-involution prefix. Checked directly: a_B = a_A and l_inv^B = l_inv^A for every element of B3..B6 (`bd_vs_A.py`, 0 mismatches). For D_n this is false (D4: 55 ancestor mismatches; D_n is not a fixed-point subgroup of S_{2n} in the same way).
- **Consequence (proposed, given Theorem A):** l_inv^{B_n}(w) = l_inv^{A_{2n-1}}(w) <= disp_{S_{2n}}(w) <= 2n - 1. This is the "partial result" of BRIEF P4, now with a proof route. The rank bound n needs a better statistic.

### 2.2 The least admissible statistic G

G = least function with G(e) = 0, Bruhat-monotone, G(w) >= 1 + G(rem w) (`npcox.py: envelope_G`, numpy, exact; agrees with `coxsp.py` on B2..B5). Any Bruhat-monotone statistic that decreases at every ancestor step is >= G, so **max G <= rank is necessary and sufficient for the "monotone statistic" route** (BRIEF P4 kill rule).

| Group | order | max l_inv | max G | rank | G - negmax (counts) | Status |
| --- | --- | --- | --- | --- | --- | --- |
| B5 | 3840 | 5 | 5 | 5 | 0: 3757, +1: 83 | numerical (exact) |
| B6 | 46080 | 6 | 6 | 6 | 0: 44470, +1: 1610 | numerical (exact) |
| B7 | 645120 | 7 | 7 | 7 | 0: 614374, +1: 30746 | numerical (exact) |
| B8 | 10321920 | 8 | 8 | 8 | 0: 9705213, +1: 614403, **+2: 2304** | numerical (exact) |
| D4 | 192 | 3 | 3 | 4 | | numerical (exact) |
| D5 | 1920 | 5 | 5 | 5 | | numerical (exact) |
| D6 | 23040 | 5 | 5 | 6 | | numerical (exact) |
| D7 | 322560 | 7 | 7 | 7 | | numerical (exact) |
| D8 | 5160960 | 7 | not computed (run stopped at the orchestrator's request, host overloaded) | 8 | | open |

So the route stays alive through B8 and D7, but G has no closed formula yet.

### 2.3 Candidates tested

| Candidate | Result | Evidence |
| --- | --- | --- |
| negmax(w) = max_i { abs(w(i) - i) if w(i) > 0; max(i, abs(w(i))) if w(i) < 0 } (signed one-line, Bruhat-monotone: it is a max over NE-quadrant conditions in S_{2n}) | **refuted** as a step-decreasing statistic: 2 step failures in B4, 36 in B5; also l_inv > negmax for 2 elements of B4 (e.g. w = (2,4,-3,1): l_inv = 4, negmax = 3) | `bd_cands.py` |
| G restricted to bigrassmannian elements, f*(x) = max G(beta) over bigrassmannian beta <= x | refuted in B4 (2 step failures); G is **not** of "max over single cells" form | `bd_explore.py` |
| G vs negmax | in B5..B7: G - negmax is 0 or 1; in B8 it reaches 2 (2304 elements), so no 'negmax + 0/1 correction' formula can be uniform in n; in D_n the difference takes values -1, 0, +1 | `bd_G.py` |
| BRIEF's earlier candidates (half displacement, unsigned + sign flag, S_2n displacement) | already refuted or too weak (BRIEF §4, L6) | — |

Obstructions (minimal elements of {G >= k} with negmax < k), B6 (`bd_exceptional.py`): (1,5,-3,2,4,6), (2,4,-3,1,5,6) [k=4]; (1,2,6,-4,3,5), (1,6,-4,2,3,5), (2,4,5,-3,1,6), (1,3,5,-4,2,6) [k=5]; (2,4,5,6,-3,1), (1,3,5,6,-4,2), (1,2,4,6,-5,3), (2,4,6,-5,-1,3) [k=6]. All but the last have length 2k - 1. D6 has a single one: (5,-4,-2,1,3,6) with G = 5. A formula for G must add +1 exactly above these patterns; the type-A proof template (one involution prefix per up-set) would then need a Lemma K for each pattern.

## 3. R4 survey (minimal number of involution factors in a length-additive factorization)

m(w) computed by two methods in `r4.py`: (A) exact Bellman–Ford relaxation over all involution prefixes (bitmask prefix test, self-tested against the length identity: 0 mismatches), (B) memoised exact recursion restricted to the elements with l_inv > rank. Both agree on every element they share.

| Group | rank | order | max l_inv | max m | elements with m > rank | method | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| I2(m), 3 <= m <= 20 | 2 | 2m | 2 | 2 | 0 | A (`out/r4_dihedral.log`; prefix self-test 0 mismatches) | numerical (exact) |
| H3 | 3 | 120 | 4 | 3 | 0 | A + B (CERT: 3) | numerical (exact); agrees with CERT.md |
| F4 | 4 | 1152 | 5 | 4 | 0 | A + B (CERT: 4) | numerical (exact); agrees with CERT.md |
| H4 | 4 | 14400 | 7 | 4 (distribution 1, 571, 6878, 6828, 122) | 0 | A + B (CERT: 4) | numerical (exact); agrees with CERT.md |
| E6 | 6 | 51840 | 8 | 6 | 0 | A + B (CERT: 6) | numerical (exact); agrees with CERT.md |
| E7 | 7 | 2903040 | 9 | **7** (distribution 1, 10207, 467170, 2126334, 294496, 4666, 162, 4) | 0 | A (5 sweeps, 123,130,795 prefix pairs) + B (exact m on all 10,997 elements with l_inv > 7: max 6; agreement A = B on 3000 of them) | numerical (exact), two methods agree |
| E8 | 8 | 696729600 | >= 9 | not computed | — | `e8_search.py` written (involutions by conjugacy closure of the w0(J), per-element R4 test) but stopped before output at the orchestrator's request; exhaustive search infeasible here (BRIEF P3) | open |
| A_n, n <= 9 | n | | n | n | 0 | m <= l_inv <= n (exhaustive l_inv) and m(s_1...s_n) = n; method A exact for A4..A6 (max m = n) | numerical (exact) |
| B_n, n <= 8 | n | | n | n | 0 | m <= l_inv <= n (exhaustive l_inv) and m(s_0...s_{n-1}) = n; method A exact for B4, B5 | numerical (exact) |
| D_n, n <= 8 | n | | 3,5,5,7,7 | D4 3, D5 5, D6 5 (method A); D7, D8 <= 7 | 0 | method A for D4..D6; m <= l_inv <= rank (exhaustive l_inv, BRIEF/big.py and npcox) for D7, D8 | numerical (exact) |

(Lower bounds m(c) = n for the path-order Coxeter elements of A_n, B_n: the word s_1 s_2 ... s_n has no commuting adjacent letters and no repeated letter, so it is the only reduced word; its prefixes s_1...s_k (k >= 2) are not involutions; so the only length-additive involution factorization is into single letters. Status: proposed (elementary).)

### 3.2 Scope-gate note: odd–even transposition sort (checked)
Odd–even transposition sort writes w in S_N as a length-additive product of its non-empty rounds; each round is a product of commuting adjacent transpositions, hence an involution. This gives m(w) <= N = rank + 1 in type A. `oddeven_check.py` (S_3..S_7, `out/oddeven_check.json`): the sort needs N rounds for some w (max rounds 3, 4, 5, 6, 7 = rank + 1; 1, 3, 24, 77, 868 elements need rank + 1 rounds), while max l_inv = rank in each case. So the observation (m <= rank + 1) is true but weaker by one than Theorem A (m <= l_inv <= rank, proposed; numerical through S_10). No conflict. Status: numerical (exact).

## 4. Commands run

All from `targets/hart-rowley-ancestor/plans/C1/` (python 3.12.10, numpy 2.5.5; python-flint and sympy not needed: every computation is integer/permutation arithmetic). Logs in `out/`.

| Command | Output | Result | Status |
| --- | --- | --- | --- |
| `python typeA_dp.py N`, N = 3..10 | stdout, `out/typeA_dp_9_10.log`, `out/A{N}_dp_anc.npy` | 0 failures (D), (I), (P); 0 ties | numerical (exact) |
| `python typeA_gg.py exhaustive N` (N = 3..8 with Lemma K, 9, 10 with `noK`) | `out/typeA_gg_8_10.log`, `out/A{N}_gg_anc.npy` | 0 failures | numerical (exact) |
| `python compare_A.py 3 4 5 6 7 8 9 10` | `out/compare_A.json` | 0 mismatches (S_10: 3,628,800 elements) | numerical (exact), two implementations agree |
| `python typeA_K_vec.py N`, N = 3..10 | stdout, `out/typeA_K_9_10.log` | Lemma K: 0 failures over 4,850,718 pairs (S_10) | numerical (exact) |
| `python typeA_replay.py 7` | stdout | proof replay, 0 failures | numerical (exact) |
| `python typeA_gg.py random N COUNT SEED [band]` (7 runs) | `out/typeA_random.log` | 0 failures | numerical (exact) |
| `python bd_vs_A.py B n` (n = 3..6), `D 4`, `D 5` | stdout | a_B = a_A: 0 mismatches; D: mismatches (expected) | numerical (exact) |
| `python b_weak_restrict.py` | stdout | 0 failures | numerical (exact) |
| `python bd_explore.py B n`, `python bd_cands.py B n` (n = 2..5) | stdout | negmax and f* refuted as step-decreasing | refuted (route), numerical |
| `python bd_G.py T n` (B5..B8, D4..D7; D8 stopped) | `out/bd_G_small.log`, `out/bd_G_B8.log`, `out/bd_G_D8.log` | max G = rank (B), <= rank (D) | numerical (exact) |
| `python bd_exceptional.py B 6`, `D 6` | stdout | obstruction list (section 2.3) | numerical (exact) |
| `python oddeven_check.py` | `out/oddeven_check.json` | sort needs up to rank + 1 rounds; l_inv <= rank | numerical (exact) |
| `python r4.py T n full|both` (H3, F4, H4, E6, E7, A4..A6, B4, B5, D4..D6, I2(3..20)) | `out/r4_*.log`, `out/R4_*_m.npy` | no R4 violation | numerical (exact) |
| `python e8_search.py 1 150` | none (stopped) | not run to completion | open |

Engines: `typeA_dp.py` (I1), `typeA_gg.py` (I2), `coxsp.py` (pure python, small groups), `npcox.py` (numpy, whole groups up to ~10^7). None of them imports code from `../../scripts/` or `../../cert/`. Gaetz–Gao source: `src_gg/` (arXiv e-print of 2610.22857v1, fetched 2026-11-27; `separable3.tex` read in full).

## 5. Next steps

1. **Type A to the Prover** (stage 4): Lemma K + section 1.3 is a complete argument modulo Gaetz–Gao Cor. 1.2 and the complementation step of their Lemma 3.1 (unrefereed preprint; the complementation step is short and could be re-proved inline). Certifier: archive the S_<=10 two-implementation agreement. Formaliser: Lemma K is pure permutation combinatorics and a reasonable Lean target; the Bruhat parts are not in Mathlib.
2. **Type B**: the problem is now "l_inv^A(w) <= n for centrally symmetric w in S_2n". Find a formula for G (data: G - negmax in {0,1}; the obstruction patterns of 2.3) and then a Lemma K per obstruction pattern; or find a different family of up-sets. Second Explorer run on P4 is justified (route alive, max G = rank through n = 8).
3. **Type D**: no reduction to type A (ancestors differ); needs its own statistic; max G <= rank through D8.
4. **R4**: holds everywhere computed; E8 exhaustively needs the C program of BRIEF P3 (not attempted).

### 3.1 Notes on R4
- In every exceptional group the maximum of m equals the rank exactly (H3 3, F4 4, H4 4, E6 6, E7 7), while the greedy l_inv exceeds it. The maximisers are not Coxeter elements in E6, E7 (Coxeter elements reach only 5 and 6); shortest E7 maximiser (0-based Bourbaki word) 6 5 4 3 2 1 3 4 0 2 3 1, length 12.
- Classical exact distributions (method A): D4 1, 43, 124, 24; D5 1, 155, 1080, 636, 46, 2; D6 1, 751, 10552, 10978, 736, 22; A5 1, 75, 414, 206, 22, 2; A6 1, 231, 2428, 2124, 228, 26, 2; B5 1, 311, 2290, 1194, 42, 2.
- R4 is the reading under which the source's Magma claim would have held (CERT.md, heuristic). It is now consistent with every group computed.
