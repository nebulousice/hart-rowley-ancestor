# Reproducing the type A paper

Everything here is corroboration only. **No step of any proof in `paper.tex` or `PROOF.md` uses a computation.** The scripts are one exact implementation (pure Python, standard library, integer arithmetic), not an independent pair of certificates.

## 1. The paper

Requirements: [Tectonic](https://tectonic-typesetting.github.io/) (any recent release; it fetches the LaTeX packages it needs on first use). The source uses only pdfLaTeX-compatible packages (`amsart`, `geometry`, `fontenc`, `lmodern`, `amssymb`, `amsthm`, `booktabs`, `xurl`, `hyperref`), so `latexmk -pdf paper.tex` also works.

From the root of the result repository:

```
cd typeA
tectonic -X compile paper.tex
```

Expected: `paper.pdf`, 18 pages, with no errors, no undefined references or citations, and no overfull boxes. (Built on 2026-11-27 and rebuilt on 2026-11-28 with Tectonic; the only console line besides progress notes was a harmless `Fontconfig error: Cannot load default config file` on Windows.)

## 2. The corroboration scripts

Requirements: Python 3.12 (the archived runs used Python 3.12.10), standard library only.

From the root of the result repository:

```
cd typeA/scripts
python small_cases.py 6 tables > out_small_cases_N6.txt
python small_cases.py 7        > out_small_cases_N7.txt
python typeB_cases.py 4        > out_typeB_cases_n4.txt
```

Runtimes on a Windows 11 desktop (rerun 2026-11-28): `small_cases.py 6 tables` 8 s, `small_cases.py 7` 78 s, `typeB_cases.py 4` 12 s (the original runs of 2026-11-27 took 11 s, not timed, and 16 s).

Expected output: the three files must be **byte-identical** to the archived copies in `scripts/` (the scripts are deterministic). In the rerun of 2026-11-28 all three were identical. Each output file is one JSON line per case.

**`out_small_cases_N6.txt` and `out_small_cases_N7.txt`** (symmetric groups $S_N$, $N = 2, \dots, 6$ resp. $2, \dots, 7$; the N6 file also prints the $S_3$ and $S_4$ tables of the paper). Every key ending in `_fail`, and the keys `ancestor_ties`, `disp_rem_not_smaller`, `linv_gt_disp`, `linv_gt_rank`, must be `0`, and `coxeter_chain_ok` must be `true`. The last line of the N7 file is

```
{"N": 7, "order": 5040, "involutions": 231, "ancestor_ties": 0, "hist_linv": [1, 231, 1212, 1906, 1098, 438, 154], "max_linv": 6, "gg1_fail": 0, "suffix_not_bruhat_below_fail": 0, "disp_suffix_fail": 0, "gg2_instances": 467304, "gg2_fail": 0, "lemmaK_pairs": 7176, "lemmaK_fail": 0, "phi_equivariance_fail": 0, "disp_rem_not_smaller": 0, "linv_gt_disp": 0, "linv_gt_rank": 0, "coxeter_element_linv": 6, "coxeter_chain_ok": true}
```

The key `lemmaK_pairs` (2, 8, 38, 182, 1064, 7176 for $`N = 2, \dots, 7`$) counts tests of Lemma 4.2 of the paper, one for each position of largest displacement of each non-identity $`w`$: excedances $`w(p) = p + \mathrm{disp}(w)`$ directly, and drops $`w(p) = p - \mathrm{disp}(w)`$ through the symmetry $`\varphi`$. It is not the number of pairs with $`w(p) = p + \mathrm{disp}(w)`$ alone.

The other key names follow the numbering of `PROOF.md`, which differs from the paper: `lemma21_*` is paper Lemma 3.1 (inversions of a product), `lemma22_*` is paper Lemma 3.2 (prefix criterion), `lemma24_*` is paper Lemma 3.4 (Bruhat order and rank matrices), `ancestor_ties` checks the uniqueness in paper Theorem 2.11(1), `gg1_*` checks paper Theorem 2.11(2) (every involution prefix is Bruhat-below the ancestor), `gg2_*` checks the complementation step, paper Lemma 2.10, and `lemmaK_*` is paper Lemma 4.2.

Distributions of the involution length (`hist_linv`, counts for the values 0, 1, 2, ...): $S_2$ 1, 1; $S_3$ 1, 3, 2; $S_4$ 1, 9, 10, 4; $S_5$ 1, 25, 48, 36, 10; $S_6$ 1, 75, 240, 264, 104, 36; $S_7$ 1, 231, 1212, 1906, 1098, 438, 154. The maximum is $N-1$ in each case.

**`out_typeB_cases_n4.txt`** (types $B_1, \dots, B_4$ and $D_2, D_3, D_4$). For the $B_n$ rows every `_fail` and `_mismatch` key, and `linvA_gt_dispA`, `linvW_gt_dispA`, must be `0`; `hist_linv_W` is 1, 1 / 1, 5, 2 / 1, 19, 20, 8 / 1, 75, 156, 123, 29 for $n = 1, \dots, 4$. For type D the expected non-zero values are: `D3`: `weak_order_mismatch` 32, `ancestor_mismatch_W_vs_A` 6, `linv_mismatch_W_vs_A` 4, and `shortest_mismatch` with window `[-3, -1, 2]`, `linv_W` 2, `linv_A` 3; `D4`: 788, 55, 33.

## 3. Lean

Nothing in this directory is formalized in Lean. The only Lean development connected with Hart–Rowley Conjecture 1.3 in this repository is the $H_3$ counterexample of the companion note, under `../lean/` (the `lean/` directory at the root of the result repository, next to the companion note `../paper.tex`).
