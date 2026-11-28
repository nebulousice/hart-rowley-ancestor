# How to reproduce

All commands are run from the root of this repository unless a `cd` is shown. The archived outputs were produced on 2026-11-27 on Windows 11 (10.0.26200). Runtimes are wall-clock times on a shared, often heavily loaded machine; expect them to be shorter elsewhere.

## Versions

| Tool | Version used | Needed for |
| --- | --- | --- |
| Python | 3.12.10 | everything in `cert/`, `scripts/`, `work-in-progress/` |
| Python standard library only | — | the certificates in `cert/` (no numpy, sympy or flint) |
| NumPy | 2.5.5 | `scripts/big.py`, most of `work-in-progress/` |
| SymPy | not recorded | `scripts/second_impl.py` only (Todd–Coxeter enumeration) |
| Lean 4 (through `elan`) | `leanprover/lean4:v4.37.0-rc1` (from `lean/lean-toolchain`) | `lean/` |
| Mathlib | `8ec7bd8e5f837ab7e335374928e18341c4e9f385` (pinned in `lean/lake-manifest.json`) | only by the `lake` route below; no Lean file here imports it |
| bash, git, sha256sum | any | only to regenerate the archived headers with `cert/run_cert.sh` |

## 1. The certificates (the main result)

These two programs share no code. Both use integer arithmetic only.

```
cd cert
python hr_tables.py --big     # primary: A6 example, H3, F4, H4, E6   (206 s)
python hr_tables_x.py         # cross-check: A6 example, H3, F4        (23 s)
python hr_tables_x.py --h4    # cross-check, adding H4                 (34 s)
```

Expected last lines, exactly as in the archived outputs:

`python hr_tables.py --big` ([`cert/hr_tables.out`](cert/hr_tables.out)):
```
CERTIFIED: H3 max l_inv = 4 (3 elements > 3); F4 max l_inv = 5 (50 elements > 4)
CERTIFIED: H4 max l_inv = 7 (1746 elements > 4); E6 max l_inv = 8 (160 elements > 6)
PASS
```

`python hr_tables_x.py` ([`cert/hr_tables_x.out`](cert/hr_tables_x.out)):
```
CERTIFIED: H3 max l_inv = 4 (3 elements > 3); F4 max l_inv = 5 (50 elements > 4)
PASS
```

`python hr_tables_x.py --h4` ([`cert/hr_tables_x_h4.out`](cert/hr_tables_x_h4.out)):
```
CERTIFIED: H4 max l_inv = 7 (1746 elements > 4)
CERTIFIED: H3 max l_inv = 4 (3 elements > 3); F4 max l_inv = 5 (50 elements > 4)
PASS
```

Each program exits with a non-zero status, and prints no `PASS`, if any check fails. Other lines to compare: the witness line `H3 witness w = r2r1r2r3r1r2r1: l(w) = 7, decomposition ['r2r1r2', 'r1r3', 'r2', 'r1'], lengths [3, 2, 1, 1], l_inv = 4`, the `A6 example` line with `matches source ... True`, and `self-test prefix criterion ... mismatches = 0`.

**With the archived header** (command, date, commit, script hashes, versions, exit status, runtime). Run inside a git clone, and give a **new** output name so that the archived file is not overwritten:

```
cd cert
bash run_cert.sh hr_tables.py rerun_hr_tables.out --big
diff <(sed -n '/^----$/,/^----$/p' hr_tables.out) <(sed -n '/^----$/,/^----$/p' rerun_hr_tables.out)
```

The `diff` compares the bodies only; it should print nothing. The script hashes in the header should match those in the archived files:

```
07709f7c448cf0311e2321cbcff5e0e18bf357f7b2a18ac6f7b5adca36730fe4  hr_tables.py
939625ab44fa9788a05fd281c225ed84a140ad110bdcc724610e2e379c75b683  hr_tables_x.py
c1b55292378112c5f8d50e2f02772dfbeb37cbefe2d301191a168520ad5558ee  hr_lean_gen.py
62623851347d29339bd9e7fcf266c7ccd807949c53407aa7c759fb4dbc1b87d3  run_cert.sh
```

## 2. The exploratory scripts ($E_6$ second implementation; $E_7$ and $E_8$, numerical)

```
cd scripts
python rows.py           # all irreducible groups of order < 100,000, incl. the E6 row   (5.1 s)
python big.py E 7        # whole E7, numpy engine                                        (27.8 s)
python check_word.py     # single-element checks of the witnesses, incl. E7 and E8
python second_impl.py    # Todd-Coxeter implementation: A3, A4, B3, D4, H3, F4 (needs sympy)
```

Expected lines:
- `rows.py` ([`scripts/out/rows_2026-11-27.txt`](scripts/out/rows_2026-11-27.txt)):
  ```
  E6     51840  36  6    8         160            yes              [1, 891, 6664, 15236, 16850, 9622, 2416, 150, 10]
  ```
- `big.py E 7` (one JSON line; archived in [`scripts/out/big_runs_2026-11-27.jsonl`](scripts/out/big_runs_2026-11-27.jsonl)): `"order": 2903040`, `"ancestor_nonunique": 0`, `"max_linv": 9`, `"count_gt_rank": 10997`, `"hist_linv": [1, 10207, 139442, 531987, 919660, 802863, 385355, 102528, 9566, 1431]`. Timing fields will differ.
- `check_word.py` ([`scripts/out/check_word_2026-11-27.txt`](scripts/out/check_word_2026-11-27.txt)):
  ```
  E7  word r7 r5 r6 r4 r5 r3 r4 r1 r3 r2 r4 r5 r6 r7  length 14  l_inv = 9  (rank 7)  ancestor lengths [2, 2, 2, 2, 2, 1, 1, 1, 1]  #prefixes per step [27, 21, 16, 12, 8, 5, 4, 3, 2]
  E8  word r7 r5 r6 r4 r5 r3 r4 r1 r3 r2 r4 r5 r6 r7  length 14  l_inv = 9  (rank 8)  ancestor lengths [2, 2, 2, 2, 2, 1, 1, 1, 1]  #prefixes per step [27, 21, 16, 12, 8, 5, 4, 3, 2]
  ```
- `second_impl.py` ([`scripts/out/second_impl_2026-11-27.txt`](scripts/out/second_impl_2026-11-27.txt)): `H 3 120 15 {0: 1, 1: 31, 2: 52, 3: 33, 4: 3}` and `F 4 1152 24 {0: 1, 1: 139, 2: 374, 3: 396, 4: 192, 5: 50}`.

These outputs have no generated header. They were produced with Python 3.12.10 and NumPy 2.5.5.

## 3. Lean

The Lean files are [`lean/Pipeline/HartRowley/Main.lean`](lean/Pipeline/HartRowley/Main.lean) and [`Data.lean`](lean/Pipeline/HartRowley/Data.lean). They import no library, not even Mathlib.

**What was copied and what was added.** `lean/lean-toolchain`, `lean/lakefile.toml` and `lean/lake-manifest.json` are unmodified copies from the shared Lean project where the check was first run. That project's root module also imports unrelated targets, so it was not copied. Instead, [`lean/Pipeline.lean`](lean/Pipeline.lean) is a new, minimal root module that imports only `Pipeline.HartRowley.Main`. The `lakefile.toml` builds the library `Pipeline` from that root.

### 3a. With `lake` (as in the original check)

```
bash scripts/lean_check.sh HartRowley
```

This runs `lake build` in `lean/`, scans for `sorry`/`admit`, and checks the `#print axioms` lines. Expected last line:

```
PASS lean_check HartRowley (Lean (version 4.37.0-rc1, x86_64-w64-windows-gnu, commit 0b77e12952ad39077d80d8bd6b66d601f225c321, Release))
```

Notes:
- The copied `lakefile.toml` requires Mathlib, so the first `lake build` clones the packages pinned in `lake-manifest.json`. Since nothing imports Mathlib, none of it should need compiling. **This standalone `lake` route was not run during the write-up**; the original check ran inside the shared project ([`cert/lean_check.out`](cert/lean_check.out), 47 min in total on a loaded host; the module alone took 141 s).
- The script overwrites `lean/Pipeline/HartRowley/build.log` and `axioms.log`. The original record is preserved in [`cert/lean_check.out`](cert/lean_check.out); `git checkout` restores the logs.

### 3b. With plain `lean`, no Mathlib (run for this write-up)

```
cd lean
mkdir -p build/Pipeline/HartRowley
lean -o build/Pipeline/HartRowley/Data.olean Pipeline/HartRowley/Data.lean
LEAN_PATH=build lean -o build/Pipeline/HartRowley/Main.olean Pipeline/HartRowley/Main.lean
```

`elan` selects the toolchain from `lean/lean-toolchain`. Expected output: eight lines, one per theorem, each ending in `depends on axioms: [propext, Classical.choice]`:

```
'HartRowley.relations_H3' depends on axioms: [propext, Classical.choice]
'HartRowley.table_ok' depends on axioms: [propext, Classical.choice]
'HartRowley.bfs_reaches_all' depends on axioms: [propext, Classical.choice]
'HartRowley.lenData_ok' depends on axioms: [propext, Classical.choice]
'HartRowley.involData_ok' depends on axioms: [propext, Classical.choice]
'HartRowley.witness_facts' depends on axioms: [propext, Classical.choice]
'HartRowley.h3_table' depends on axioms: [propext, Classical.choice]
'HartRowley.h3_counterexample' depends on axioms: [propext, Classical.choice]
```

Archived run: [`lean/plain_lean_check_2026-11-27.log`](lean/plain_lean_check_2026-11-27.log) (928 s on a loaded host). Check that the file hashes match:

```
ff627893aeeae97d711528b34cfa2795c4d0d29ab5ebdfa49e48f30b92f68a94  lean/Pipeline/HartRowley/Main.lean
f468e40b6b0fbf03fefd9545275737212b6849bffc57bcba011c02cfb1eb0393  lean/Pipeline/HartRowley/Data.lean
```

**Regenerating the hint data** (optional). `Data.lean` is generated, and every entry is re-verified in `Main.lean`:

```
cd cert
python hr_lean_gen.py ../lean/Pipeline/HartRowley/Data.lean
```

## 4. Exploratory checks for types $A$, $B$, $D$ and the minimal count

[`work-in-progress/`](work-in-progress/) holds the exploratory programs for types $A$, $B$, $D$ and for the minimal number of involution factors. Its notes are [`work-in-progress/LAB.md`](work-in-progress/LAB.md), and its section 4 lists every command. The large `.npy` arrays written by these scripts are not included; they are regenerated by the commands. Main checks, with archived logs in [`work-in-progress/out/`](work-in-progress/out/):

```
cd work-in-progress
python typeA_K_vec.py 10      # key type-A lemma on S_10: "pairs": 4850718, "k1"..."k4": 0   (85 s)
python compare_A.py 3 4 5 6 7 8 9 10   # needs the .npy arrays from typeA_dp.py and typeA_gg.py first
python r4.py E 7 both         # minimal factorization count in E7: "full_max_m": 7         (about 52 min)
```

The type-$A$ theorem and the type-$B$ reduction are **proved** by a written argument that passed AI review; these runs are corroboration only, and no proof step uses them. The counts are **numerical**. The exact checks of $S_N$ for $N\le7$, of $B_n$ for $n\le4$ and of $D_3$, $D_4$ cited in [`README.md`](README.md) were made with scripts that are not yet included in this folder. See the research log in [`README.md`](README.md).
