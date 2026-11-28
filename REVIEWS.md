# Review history

*Updated 2026-11-27.*

## Current status

- **The type $A$ proof and its type $B$ consequence passed review.** The process: two AI referees on different models review the full written proof blind. If they do not both accept with zero issues, the issues are fixed, and three referees (the two plus a third, on a third model) review the fixed version. At least two of the three must accept with zero issues, and any major issue blocks.
  - **Round 1** (commit `4101231` of the private working repository). Claude Opus checker: minor, 5 issues, about citations, the record, and how the dependence on Gaetz–Gao was handled; none affected the mathematics, and all 5 were fixed. Claude Fable adversary: accept, 0 issues. The main fix was an appendix that proves Gaetz–Gao Cor. 1.2 in full from textbook facts, so the proof no longer has an outside hypothesis.
  - **Round 2** (commit `599cc84`). Claude Opus checker, Claude Fable adversary and Claude Sonnet arbiter: all three accept with zero issues.
- **The counterexample write-up has not had a referee review.** This covers the hand proof of the $H_3$ case, the certificates and the Lean files; their checks are listed below.
- **Expert verification required.** Verified so far with: Lean 4 (the H3 counterexample); two independent exact computer certificates; AI review of the type A proof (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 2: all three accept with zero issues). The type A theorem is not formalized in Lean.

## Other checks (these are not referee reviews)

| Date | Check | By | Outcome |
| --- | --- | --- | --- |
| 2026-11-27 | Exploration: the source's first values, the counterexamples, a hand check of the $H_3$ case | Claude (AI) | Found the violations; flagged that the source's reported Magma verification does not reproduce |
| 2026-11-27 | Definitions audit and certification: is there a reading of the source's definitions under which the conjecture holds? Two independent exact programs | Claude (AI) | One reading only, fixed by the source's own $A_6$ example; every greedy variant fails in $H_3$ and $F_4$; both programs agree on $H_3$, $F_4$, $H_4$ ([`cert/CERT.md`](cert/CERT.md)) |
| 2026-11-27 | Lean 4 check of the $H_3$ case | Claude (AI) | `PASS`, standard axioms only ([`LEAN.md`](LEAN.md)) |
| 2026-11-27 | Priority searches (three runs) | Claude (AI) | Found Levinstein's earlier public refutation ([`LITERATURE_CHECK.md`](LITERATURE_CHECK.md)) |
| 2026-11-27 | Write-up | Claude (AI) | `paper.tex` compiled to [`paper.pdf`](paper.pdf) |

## Known weak points a referee should look at

1. **The hand proof of the $H_3$ case** ([`PROOF.md`](PROOF.md) §3) relies on listing complete braid classes by hand. Every class is short and written out in full, so this can be checked line by line.
2. **The $E_6$ row is certified by two separately written programs that both model the group by permutations of roots.** The definition-level cross-check (matrices, breadth-first search) was too slow for $E_6$ and was stopped. Our records name the second program as the exploratory `big.py`; the archived $E_6$ output of that code base is [`scripts/out/rows_2026-11-27.txt`](scripts/out/rows_2026-11-27.txt), produced by `rows.py` with `coxlib.py`.
3. **The $E_7$ and $E_8$ values are numerical.** They rest on one whole-group program (`big.py`, which identifies elements by hashed keys, with the final count checked against the known group order) and on a single-element checker.
4. **The Lean statement is about a concrete permutation group.** Its identification with $H_3$ uses the cited fact $\lvert W(H_3)\rvert=120$, and two short side arguments (step bound and closure) are not stated in Lean ([`LEAN.md`](LEAN.md)).
5. **The heuristic explanation** of the source's Magma claim (a non-greedy count would satisfy the bound) is a guess, not a finding.
6. **The type $A$ proof is not formalized in Lean**, and its full written text is not yet included in this folder ([`README.md`](README.md) gives the statement and the idea of the proof).
7. **Book references in the type $A$ proof.** The item numbers cited from Björner–Brenti, *Combinatorics of Coxeter Groups*, are not yet confirmed against the book. The referees noted this; the facts used are stated in full, and the proof does not depend on the numbering.
8. **The re-proof of Gaetz–Gao Cor. 1.2** follows Gaetz and Gao's argument. Their paper is an unrefereed preprint and states that ChatGPT produced the proof of its Theorem 1.1 ([`LITERATURE_CHECK.md`](LITERATURE_CHECK.md)). All three referees read the re-proof line by line in round 2 and found no gap.
