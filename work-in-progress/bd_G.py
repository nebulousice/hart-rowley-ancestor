"""Least admissible statistic G (Bruhat-monotone, G(w) >= 1 + G(rem w)) for
B_n, D_n with the numpy engine; reports max G vs rank (P4 kill rule: max G >
rank means no Bruhat-monotone step-decreasing statistic bounded by the rank
exists), and how G compares with the candidate negmax.
Usage: python bd_G.py B 7
"""
import sys, json, time
import numpy as np
from npcox import NpGroup

typ, n = sys.argv[1], int(sys.argv[2])
t0 = time.time()
G = NpGroup(typ, n, 'signed')
linv = G.ancestors()
Gv = G.envelope_G()
# signed one-line: points coded +i -> n-1+i, -i -> n-i
pos = G.E[:, n:2 * n].astype(np.int64)           # images of +1..+n
val = np.where(pos >= n, pos - n + 1, pos - n)  # signed values
i = np.arange(1, n + 1)[None, :]
negmax = np.where(val > 0, np.abs(val - i), np.maximum(i, -val)).max(axis=1)
res = dict(type=typ + str(n), order=G.order, ties=G.ties, max_linv=int(linv.max()),
           max_G=int(Gv.max()), rank=n, hist_G=np.bincount(Gv).tolist(),
           G_minus_negmax=dict(zip(*[a.tolist() for a in np.unique(Gv - negmax, return_counts=True)])),
           linv_gt_G=int((linv > Gv).sum()), seconds=round(time.time() - t0, 1))
print(json.dumps(res), flush=True)
