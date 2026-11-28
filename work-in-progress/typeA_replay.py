"""Replay of the proof of Theorem A on all of S_N (N <= 8): for every w and
every extremal position p (excedance or drop), with iota from Lemma K:
  (r1) iota <=_B a(w)                 [GG1 instance]
  (r2) rem(w) <=_B iota w             [GG2 instance]
  (r3) rem(w) <=_B w and disp(rem) <= D, (iota w)(p) != w(p)
Bruhat order by the tableau criterion; a(w) from Gaetz-Gao (typeA_gg)."""
import sys, json, itertools
import numpy as np
from typeA_gg import ancestor, check_K, iota_exc, bruhat_leq

out = {}
for N in range(3, int(sys.argv[1]) + 1):
    P = np.array(list(itertools.permutations(range(N))), dtype=np.int16)
    A = ancestor(P)
    cnt = r1 = r2 = r3 = 0
    w0 = lambda i: N - 1 - i
    for i in range(len(P)):
        w = tuple(int(x) for x in P[i]); a = tuple(int(x) for x in A[i])
        D = max(abs(w[j] - j) for j in range(N))
        if D == 0: continue
        rem = tuple(a[w[j]] for j in range(N))         # a^-1 w = a w
        for p in range(N):
            if w[p] == p + D:
                io = iota_exc(w, p)
            elif w[p] == p - D:
                phiw = tuple(w0(w[w0(j)]) for j in range(N))
                io2 = iota_exc(phiw, w0(p)); io = [w0(io2[w0(j)]) for j in range(N)]
            else:
                continue
            cnt += 1
            iw = tuple(io[w[j]] for j in range(N))
            r1 += not bruhat_leq(tuple(io), a)
            r2 += not bruhat_leq(rem, iw)
            r3 += not (bruhat_leq(rem, w) and max(abs(rem[j] - j) for j in range(N)) <= D and iw[p] != w[p])
    out['S%d' % N] = dict(extremal_pairs=cnt, r1_fail=r1, r2_fail=r2, r3_fail=r3)
print(json.dumps(out))
