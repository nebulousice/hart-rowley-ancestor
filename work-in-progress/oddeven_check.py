"""Scope-gate note: odd-even transposition sort writes w as a length-additive
product of at most N rounds, each round a product of commuting adjacent
transpositions (an involution), so m(w) <= N = rank + 1 in S_N.  Compare with
l_inv(w) <= disp(w) <= N - 1 (Theorem A) on all of S_N, N <= 7 (light check)."""
import itertools, json
import numpy as np
from typeA_gg import ancestor, compose, inverse

out = {}
for N in range(3, 8):
    P = np.array(list(itertools.permutations(range(N))), dtype=np.int16)
    rounds = []
    for w in P.tolist():
        x = list(w); r = 0; k = 0
        while x != sorted(x):
            sw = False
            for i in range(k % 2, N - 1, 2):
                if x[i] > x[i + 1]:
                    x[i], x[i + 1] = x[i + 1], x[i]; sw = True
            k += 1
            r += sw            # count non-empty rounds (each is an involution factor)
        rounds.append(r)
    linv = np.zeros(len(P), int); cur = P.copy(); alive = np.ones(len(P), bool); alive[0] = False
    while alive.any():
        idx = np.flatnonzero(alive); linv[idx] += 1
        a = ancestor(cur[idx]); cur[idx] = compose(inverse(a), cur[idx])
        alive = ~(cur == np.arange(N)[None, :]).all(axis=1)
    rounds = np.array(rounds)
    out['S%d' % N] = dict(max_oddeven_rounds=int(rounds.max()), rank=N - 1, max_linv=int(linv.max()),
                          oddeven_rounds_gt_rank=int((rounds > N - 1).sum()))
print(json.dumps(out))
