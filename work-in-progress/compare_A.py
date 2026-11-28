"""Compare the ancestor tables of the two type-A implementations."""
import sys, json, numpy as np
out = {}
for N in map(int, sys.argv[1:]):
    a = np.load('out/A%d_dp_anc.npy' % N); b = np.load('out/A%d_gg_anc.npy' % N)
    out['S%d' % N] = dict(elements=len(a), mismatches=int((a != b).sum()))
print(json.dumps(out))
