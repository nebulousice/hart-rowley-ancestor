"""Is the B_n (resp. D_n) ancestor equal to the type-A ancestor of the same
element viewed in S_{2n}?  (B_n = fixed points of conjugation by w0 in S_{2n}.)
a_B from the DP in coxsp.py; a_A from Gaetz-Gao's construction (typeA_gg.py).
Usage: python bd_vs_A.py B 5
"""
import sys, json
import numpy as np
from coxsp import Group
from typeA_gg import ancestor, compose, inverse

typ, n = sys.argv[1], int(sys.argv[2])
Gp = Group(typ, n)
Gp.ancestors()
E = Gp.E.astype(np.int16)
A = ancestor(E)
# a_B(w) = w rem(w)^-1
Rm = E[np.array(Gp.rem)]
AB = compose(E, inverse(Rm))
mism = int((AB != A).any(axis=1).sum())
# also: is l_inv^A == l_inv^B ?
linvA = np.zeros(Gp.order, int)
cur = E.copy(); alive = np.ones(Gp.order, bool); alive[0] = False
while alive.any():
    idx = np.flatnonzero(alive)
    linvA[idx] += 1
    a = ancestor(cur[idx]); cur[idx] = compose(inverse(a), cur[idx])
    alive = ~(cur == np.arange(2 * n)[None, :]).all(axis=1)
print(json.dumps(dict(type=typ + str(n), order=Gp.order, ancestor_mismatch_B_vs_A=mism,
                      linv_mismatch=int((linvA != np.array(Gp.linv)).sum()), max_linvA=int(linvA.max()),
                      max_linvW=max(Gp.linv))))
