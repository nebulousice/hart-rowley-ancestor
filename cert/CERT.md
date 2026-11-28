AGREEMENT, NO DEFINITIONAL MISMATCH. The source's definitions (arXiv:2502.00777v1 TeX) admit one reading, which reproduces the source's own A6 example and its Coxeter-element results. Under that reading, two new independent exact implementations agree with each other and with the Strategist on the witness and on the maxima. Conjecture 1.3 as stated fails in H3 (maximum involution length 4, 3 elements) and in F4 (maximum 5, 50 elements). Status proposed for the ledger: **certified** (HR-W, HR-T). The H3 statement is also kernel-checked in Lean, modulo the cited fact |W(H3)| = 120. (Certifier run 2026-11-27)

# Certification: Hart–Rowley Conjecture 1.3 (stage 5)

Target `hart-rowley-ancestor`. Certifier run on 2026-11-27, commit `d72b9ae3fc7c652864d9506640a3899a12394442`; the `cert/` files are uncommitted, and their sha256 hashes are in each `.out` header. Nothing was committed or pushed.

## Claim ids

| id | Claim |
| --- | --- |
| HR-D | Definitions audit: under the source's definitions, and under every alternative greedy reading tested, Conj. 1.3 fails in H3 and F4. |
| HR-W | In H3 (graph r1 –5– r2 – r3), w = r2r1r2r3r1r2r1 has length 7, involution prefixes {r2, r2r1r2}, ancestor decomposition (r2r1r2)(r3r1)(r2)(r1), and involution length 4 > 3. |
| HR-T | The maximum involution length and the number of elements above the rank: H3 4 (3 of 120), F4 5 (50 of 1152). Also H4 7 (1746 of 14400) and E6 8 (160 of 51840); see the coverage note below for E6. |

## (a) Definitions audit (HR-D)

**Source read.** The arXiv e-print of 2502.00777v1 was fetched on 2026-11-27 (gzip sha256 `95ec26a3106271e2a2d11575b327060c9d823dc897b24dc13794c458a2cf7c81`, `Involution-prefixes.tex` sha256 `85e2195ad1bd66fc15828e6725448d29df040776e78b02e421b298a72b35de05`, 169 lines) and read in full. The journal version (Arch. Math. 2026) is still unread, as recorded in GATE.md and BRIEF.md §2.

**The definitions are explicit and admit one reading:**

1. *Prefix* (§1, verbatim in GATE.md): u is a prefix of w if some reduced word r1…rn of w has r1…rk as a reduced word of u. The abstract restates this as "w = uv for some v in W such that ℓ(w) = ℓ(u) + ℓ(v)". Prefixes are therefore **left** factors.
2. *Ancestor*: the unique involution prefix of maximal length.
3. *Decomposition*: w = a1⋯ak with a_i = a(a_i⋯a_k). So a1 = a(w), and the recursion runs on a1⁻¹w = a2⋯ak: the ancestor is **peeled on the left**, and the recursion continues on the right-hand remainder. The involution length is k.
4. "Weak (left) Bruhat order" is only the name of the order whose lower intervals [1,w] are the prefix sets. This is Björner–Brenti's ≤_R, the same convention Gaetz–Gao use. The naming does not change the objects.
5. §3 of the source names the only alternative the authors consider: the suffix version, i.e. the same definitions with suffixes. It says that this version gives different decompositions, and that its bound would follow from the prefix bound because "the involution suffixes of any element w are the involution prefixes of w⁻¹". Hence suffix-ℓ_inv(w) = ℓ_inv(w⁻¹), and both versions have the same maximum over W.

**Calibration against the source's own examples.** Both implementations reproduce the following.
- The A6 example (§3, w = r6r3r2r1r4r5):
  - the prefix reading gives (r3r6)(r2r4)(r1r5), as in the source;
  - the suffix reading gives (r3)(r2r4r6)(r1r5), as in the source;
  - swapping the two readings would swap these outputs. This **pins the side convention**: the source's "ancestor decomposition" is our R1.
- Prop. 2.2 (a(w) = ∏ D_L(w) for Coxeter elements, with D the **left** descent set): 0 failures on all 4 Coxeter elements of H3 and all 8 of F4 (also all 8 of H4 and all 32 of E6 in the primary).
- Prop. 2.5 (ℓ_inv = path length for Coxeter elements): 0 failures (primary).

**Alternative readings, all computed by both implementations:**

| Reading | H3 max | H3 elements > 3 | F4 max | F4 elements > 4 |
| --- | --- | --- | --- | --- |
| R1: source (longest involution prefix, peeled on the left) | **4** | **3** | **5** | **50** |
| R2: suffix version (longest involution suffix, peeled on the right; = R1 on w⁻¹) | 4 | 3 | 5 | 50 |
| R3: alternating, prefix-ancestor first, then suffix-ancestor of the remainder, … | 4 | 3 | 5 | 54 |
| R3b: alternating, suffix first | 4 | 3 | 5 | 54 |
| R4: **non-greedy** minimum number of involutions in any length-additive factorization | 3 | 0 | 4 | 0 |

The primary also has, under R1/R2/R3/R4: H4 7/7/7/4 and E6 8/8/7/6.

**Verdict for HR-D:** no mismatch.
- The source's reading is unambiguous, and it reproduces the source's worked example. That example distinguishes prefixes from suffixes.
- Every greedy variant (R1–R3b) violates the rank bound in H3 and F4.
- The only reading tested that satisfies the bound is R4. R4 is **not** the source's definition: the source defines the decomposition greedily through ancestors and calls it unique.
- R4 is also ≤ rank in H4 and E6. A Magma computation of the *minimal* number of factors, rather than the greedy count, would therefore have "verified" Conj. 1.3 for F4, E6, H3 and H4. This is a possible explanation of the source's evidence claim. Status **heuristic**: a guess about unseen code, not a finding.
- On the witness, R4 gives 2. The unique minimal factorization is (r2)(r1r2r1r3r2r1), while the greedy count is 4.

## (b) Agreement (HR-W, HR-T)

The two implementations below share no code with each other or with `scripts/`.

| | Primary `hr_tables.py` | Cross-check `hr_tables_x.py` |
| --- | --- | --- |
| Group model | permutations of the root system, built over Z[φ] (Cartan-type, non-symmetric for m = 4) | matrices of the Tits geometric representation over Z[φ] (H) or Z[√2] (F); no roots |
| Length | inversion count \|Φ(w)\| | breadth-first word distance in the Cayley graph (the definition of ℓ) |
| Prefix test | Φ(u) ⊆ Φ(w) (Björner–Brenti 3.1.3), self-tested against ℓ(u)+ℓ(u⁻¹w)=ℓ(w) on all 14,400 pairs of H3 and 1,327,104 pairs of F4 (0 mismatches) | walk down from w by length-decreasing right multiplications (literal "prefix of a reduced word") |
| Involution test | p∘p = id on root permutations | M² = I on matrices, cross-checked against word reversal (0 mismatches) |
| Ancestor | brute force over all involutions, uniqueness checked | first level below w containing an involution, uniqueness checked |

**Results (identical in both runs, and identical to the Strategist's `scripts/out/`):**
- H3:
  - |W| = 120, 31 involutions, maximum length 15, ancestor uniqueness failures 0;
  - ℓ_inv distribution (0..4) = 1, 31, 52, 33, 3; max 4; 3 elements > 3;
  - the violators are r2r1r2r1r3r2r1 (= the witness), r2r1r3r2r1r2r3 and r1r2r1r2r1r3r2r1, with decompositions (r2r1r2)(r1r3)(r2)(r1), (r2r1r3r2)(r1)(r2)(r3) and (r1r2r1r2r1)(r3)(r2)(r1).
- Witness w = r2r1r2r3r1r2r1:
  - ℓ(w) = 7;
  - its 9 prefixes are exactly those listed in BRIEF.md §1.3 step 2;
  - its involution prefixes are exactly r2 and r2r1r2;
  - its decomposition is (r2r1r2)(r1r3)(r2)(r1), with factor lengths 3, 2, 1, 1; r1r3 = r3r1;
  - ℓ_inv = 4.
- F4:
  - |W| = 1152, 139 involutions, maximum length 24, uniqueness failures 0;
  - distribution (0..5) = 1, 139, 374, 396, 192, 50; max 5; 50 elements > 4;
  - witness r3r4r1r2r3r2r1 = (r1r3)(r2r4)(r3)(r2)(r1), ℓ_inv = 5.
- H4 (both implementations; cross-check in `hr_tables_x_h4.out`):
  - 14400 elements, 571 involutions;
  - distribution (0..7) = 1, 571, 2854, 5032, 4196, 1567, 173, 6; max 7; 1746 elements > 4.
- E6 (primary only):
  - 51840 elements, 891 involutions;
  - distribution (0..8) = 1, 891, 6664, 15236, 16850, 9622, 2416, 150, 10; max 8; 160 elements > 6.

**Agreement on the witness and on the H3 and F4 maxima: YES.** Disagreements: none.

**Coverage note.**
- H4 is cross-checked by both new implementations. The run also agrees with the Strategist's `scripts/big.py`.
- E6 has two implementations: the primary here and the Strategist's numpy engine `scripts/big.py` (already listed in the BRIEF), and their numbers agree exactly. My definition-level cross-check was too slow for E6: the walk down from elements near w0 visits very large levels. That run was stopped after more than 20 minutes and is not archived. E7 (10997 violators, M = 9) and the bound M(E8) ≥ 9 were **not** re-certified in this run. They rest on one whole-group implementation plus the single-element checker `scripts/check_word.py`, so they remain *numerical*.

## (c) Archived runs

All runs were made from a fresh `git clone` of the local repository at HEAD `d72b9ae`, with the `cert/` scripts copied in. The outputs were copied verbatim and never edited.

| File | Command | Result | Runtime |
| --- | --- | --- | --- |
| `hr_tables.out` | `bash run_cert.sh hr_tables.py hr_tables.out --big` (A6 example, H3, F4, H4, E6) | PASS | 206 s (host shared with other jobs) |
| `hr_tables_x.out` | `bash run_cert.sh hr_tables_x.py hr_tables_x.out` (A6 example, H3, F4) | PASS | 23 s |
| `hr_tables_x_h4.out` | `bash run_cert.sh hr_tables_x.py hr_tables_x_h4.out --h4` (adds H4) | PASS | 34 s |
| `lean_check.out` | `bash scripts/lean_check.sh HartRowley` (in the working repository; header added by script) | see the Lean section | see file |

Versions: Python 3.12.10 on Windows 11 (10.0.26200), standard library only (no numpy, sympy or flint); Lean 4 `v4.37.0-rc1`, with Mathlib per `lean/lake-manifest.json`.

Other files:
- `run_cert.sh` writes the headers.
- `hr_lean_gen.py` generates `lean/Pipeline/HartRowley/Data.lean`. That file holds untrusted hints only; every entry is re-verified in Lean.

## Certified constants (for the Consolidator)

- H3 (r1 –5– r2 – r3):
  - M(H3) = 4; #{w : ℓ_inv(w) > 3} = 3;
  - distribution 1, 31, 52, 33, 3;
  - witness r2r1r2r3r1r2r1, ℓ = 7, decomposition (r2r1r2)(r3r1)(r2)(r1).
- F4 (r1 – r2 –4– r3 – r4):
  - M(F4) = 5; #{ℓ_inv > 4} = 50;
  - distribution 1, 139, 374, 396, 192, 50;
  - witness r3r4r1r2r3r2r1 = (r1r3)(r2r4)(r3)(r2)(r1).
- H4: M = 7, 1746 elements > 4, distribution 1, 571, 2854, 5032, 4196, 1567, 173, 6. Certified by both new implementations.
- E6: M = 8, 160 elements > 6, distribution 1, 891, 6664, 15236, 16850, 9622, 2416, 150, 10. Certified by the primary here together with the Strategist's `big.py`.
- Not certified here: E7 (M = 9, 10997) and E8 (M ≥ 9).

## Lean (for the Formaliser)

**`scripts/lean_check.sh HartRowley`: PASS** (`lean_check.out`).
- Files: `lean/Pipeline/HartRowley/{Main,Data}.lean`, plus one import line in `lean/Pipeline.lean`.
- The proofs use `decide +kernel` only: no `sorry`, no `native_decide`.
- `#print axioms` reports [propext, Classical.choice] for every theorem.
- Timing:
  - the module alone (`lake build Pipeline.HartRowley.Main`) took 141 s;
  - the full `lean_check` took 47 min on the shared, heavily loaded host. Of that, "Built Pipeline" took 1650 s, and the axioms step re-elaborates Main.lean;
  - later builds are cached unless `Main.lean` or `Data.lean` change.
- The working repository's HEAD advanced to `5bb7344` during the run (commits by others); the archived header records this.

What the kernel checks:
- `relations_H3`: the three root permutations satisfy the H3 Coxeter relations, with exact orders.
- `table_ok`, `bfs_reaches_all`: the generated group has exactly 120 distinct elements, and the Cayley table is recomputed from the permutations.
- `lenData_ok`: lengths equal the BFS word distance.
- `involData_ok`: the involution list is exact.
- `witness_facts`: ℓ(w) = 7; exactly two involution prefixes, r2 and r2r1r2; the decomposition (r2r1r2)(r3r1)(r2)(r1); ℓ_inv = 4.
- `h3_table`: unique ancestor for every w ≠ e, max ℓ_inv = 4, exactly 3 elements > 3.
- `h3_counterexample`: the combined statement.

The definitions in Lean are literal: prefix ⇔ ℓ(u)+ℓ(u⁻¹w)=ℓ(w), with ℓ the word length.

**One cited fact, not formalized:** |W(H3)| = 120 (Humphreys §2.11). With the checked relations it makes r_i ↦ s_i an isomorphism W(H3) ≅ G, so word length in G is the Coxeter length. Upgrading this to Mathlib's `CoxeterMatrix.H₃` via `CoxeterSystem.lift` needs a proof that the lift is injective (for example, from the order of `H₃.Group`). Note that Mathlib's H₃ uses a different labelling of the 5-edge. This is left to the Formaliser.
