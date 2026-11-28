"""Minimal elements (Bruhat) of {G >= k} with negmax < k in B_n / D_n:
the 'non-cellular' obstructions that make G differ from negmax.
Usage: python bd_exceptional.py B 6
"""
import sys, json
import numpy as np
from npcox import NpGroup
typ, n = sys.argv[1], int(sys.argv[2])
G = NpGroup(typ, n, 'signed'); linv = G.ancestors(); Gv = G.envelope_G()
pos = G.E[:, n:2 * n].astype(np.int64)
val = np.where(pos >= n, pos - n + 1, pos - n)
i = np.arange(1, n + 1)[None, :]
negmax = np.where(val > 0, np.abs(val - i), np.maximum(i, -val)).max(axis=1)
Tp = [G.E[t].astype(np.int64) for t in G.T]
cand = np.flatnonzero(Gv > negmax)
X = G.E[cand]
hasLowerSame = np.zeros(len(cand), bool)
for t in Tp:
    k = G.lookup(X[:, t])
    hasLowerSame |= (G.length[k] == G.length[cand] - 1) & (Gv[k] >= Gv[cand])
mins = cand[~hasLowerSame]
print(typ + str(n), 'excess elements', len(cand), 'minimal ones', len(mins))
for w in sorted(mins, key=lambda w: (Gv[w], G.length[w])):
    print('  G=%d negmax=%d len=%d linv=%d  w=%s' % (Gv[w], negmax[w], G.length[w], linv[w], tuple(val[w].tolist())))
