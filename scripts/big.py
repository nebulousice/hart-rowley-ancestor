"""Vectorised (numpy) exhaustive computation of ancestors and involution
lengths for a whole finite Coxeter group.  Exploration-grade ("numerical" in
the pipeline's status words): one implementation, hashed element keys.  The
final element count is compared with the known group order, which rules out
hash collisions among the elements reached.

Usage:  python big.py TYPE RANK [roots|classical]
        e.g.  python big.py E 7        python big.py A 9 classical

Conventions as in coxlib.py: prefix order u <=_R w, right multiplication
table R[s][i] = index of w_i * s, elements sorted by length.
"""
import sys
import time
import json
import numpy as np
from coxlib import gens_perm

ORDERS = {('A', n): __import__('math').factorial(n + 1) for n in range(1, 12)}
ORDERS.update({('B', n): 2 ** n * __import__('math').factorial(n) for n in range(2, 11)})
ORDERS.update({('D', n): 2 ** (n - 1) * __import__('math').factorial(n) for n in range(3, 11)})
ORDERS.update({('E', 6): 51840, ('E', 7): 2903040, ('E', 8): 696729600,
               ('F', 4): 1152, ('G', 2): 12, ('H', 2): 10, ('H', 3): 120, ('H', 4): 14400})


def enumerate_group(typ, n, model='roots'):
    g = gens_perm(typ, n, model)
    deg = len(g[0])
    dt = np.uint8 if deg < 256 else np.uint16
    gens = np.array(g, dtype=dt)
    r = gens.shape[0]
    rng = np.random.default_rng(20261127)
    mult = rng.integers(1, 2 ** 62, size=deg, dtype=np.uint64) * np.uint64(2) + np.uint64(1)

    def keyof(P):
        return (P.astype(np.uint64) * mult[None, :]).sum(axis=1, dtype=np.uint64)

    F = np.arange(deg, dtype=dt)[None, :]
    K = keyof(F)
    Kprev = np.empty(0, np.uint64)
    offs = [0]
    Rparts = [[] for _ in range(r)]
    invparts = []
    offprev = 0
    while len(F):
        k = len(F)
        cands = [F[:, gens[s]] for s in range(r)]
        allk = np.concatenate([keyof(C) for C in cands])
        if len(Kprev):
            pos = np.searchsorted(Kprev, allk)
            pos[pos == len(Kprev)] = 0
            down = Kprev[pos] == allk
        else:
            pos = np.zeros(len(allk), np.int64)
            down = np.zeros(len(allk), bool)
        upk = allk[~down]
        Knext, first = np.unique(upk, return_index=True)
        allc = np.concatenate(cands)
        Fnext = allc[~down][first]
        del allc, cands
        off_cur = offs[-1]
        off_next = off_cur + k
        idx = np.empty(len(allk), np.int32)
        idx[down] = offprev + pos[down]
        idx[~down] = off_next + np.searchsorted(Knext, upk)
        for s in range(r):
            Rparts[s].append(idx[s * k:(s + 1) * k])
        Finv = np.argsort(F, axis=1).astype(dt)
        ik = keyof(Finv)
        p = np.searchsorted(K, ik)
        assert np.all(K[p] == ik), 'inverse not found on the same level'
        invparts.append((off_cur + p).astype(np.int32))
        offprev = off_cur
        offs.append(off_next)
        Kprev, K, F = K, Knext, Fnext
    R = [np.concatenate(Rparts[s]) for s in range(r)]
    inv = np.concatenate(invparts)
    offs = offs[:-1] if offs[-1] == offs[-2] else offs
    return r, R, inv, offs


def run(typ, n, model='roots'):
    t0 = time.time()
    r, R, inv, offs = enumerate_group(typ, n, model)
    order = len(inv)
    expected = ORDERS.get((typ, n))
    assert expected is None or order == expected, (order, expected)
    nlev = len(offs) - 1
    length = np.zeros(order, np.uint8)
    for L in range(nlev):
        length[offs[L]:offs[L + 1]] = L
    t1 = time.time()
    alen = np.zeros(order, np.int16)       # length of a(w)
    rem = np.zeros(order, np.int32)        # index of a(w)^-1 w
    linv = np.zeros(order, np.uint8)
    supp = np.zeros(order, np.int32)       # support bitmask
    nonunique = 0
    for L in range(1, nlev):
        lo, hi = offs[L], offs[L + 1]
        I = np.arange(lo, hi, dtype=np.int32)
        m = hi - lo
        cl = np.full((r, m), -1, np.int16)
        cr = np.zeros((r, m), np.int32)
        sp = np.zeros(m, np.int32)
        for s in range(r):
            k = R[s][lo:hi]
            d = length[k] < L
            kd = k[d]
            cl[s, d] = alen[kd]
            cr[s, d] = R[s][rem[kd]]          # a(ws)^-1 w = (a(ws)^-1 ws) s
            sp[d] = supp[kd] | (1 << s)
        supp[lo:hi] = sp
        best = cl.argmax(axis=0)
        bl = cl.max(axis=0)
        br = cr[best, np.arange(m)]
        # uniqueness of the longest involution prefix (Conj. 1.2): all
        # descents attaining the maximal length must give the same remainder
        tie = cl == bl[None, :]
        diff = tie & (cr != br[None, :])
        isinv = inv[lo:hi] == I
        nonunique += int((diff.any(axis=0) & ~isinv).sum())
        a = np.where(isinv, L, bl)
        rm = np.where(isinv, 0, br)
        alen[lo:hi] = a
        rem[lo:hi] = rm
        linv[lo:hi] = 1 + linv[rm]
    hist = np.bincount(linv)
    mx = int(linv.max())
    pop = np.zeros(order, np.uint8)
    for s in range(r):
        pop += ((supp >> s) & 1).astype(np.uint8)
    excess_supp = int((linv.astype(np.int16) - pop).max())
    noninv = (inv != np.arange(order)) & (np.arange(order) > 0)
    same_supp = int((supp[rem[noninv]] == supp[noninv]).sum())

    def word(i):
        out = []
        while i:
            for s in range(r):
                k = int(R[s][i])
                if length[k] < length[i]:
                    out.append(s)
                    i = k
                    break
        return out[::-1]

    def decomposition(i):
        out = []
        while i:
            j = int(rem[i])
            # a = w * rem^-1 : walk the inverse of rem's word on the right
            a = i
            for s in reversed(word(j)):
                a = int(R[s][a])
            out.append(word(a))
            i = j
        return out

    wit = int(np.flatnonzero(linv == mx)[0])     # shortest element of maximal l_inv
    res = {
        'type': typ + str(n), 'model': model, 'order': order, 'N': nlev - 1, 'rank': r,
        'involutions': int((inv == np.arange(order)).sum()) - 1,
        'ancestor_nonunique': nonunique,
        'max_linv': mx, 'count_gt_rank': int((linv > r).sum()),
        'hist_linv': hist.tolist(),
        'max_linv_minus_suppsize': excess_supp,
        'noninvolutions_with_supp(rem)==supp(w)': same_supp,
        'witness_word_0based': word(wit), 'witness_length': int(length[wit]),
        'witness_decomposition': decomposition(wit),
        'seconds_enum': round(t1 - t0, 1), 'seconds_total': round(time.time() - t0, 1),
        'numpy': np.__version__, 'python': sys.version.split()[0],
    }
    return res


if __name__ == '__main__':
    typ, n = sys.argv[1], int(sys.argv[2])
    model = sys.argv[3] if len(sys.argv) > 3 else 'roots'
    print(json.dumps(run(typ, n, model)), flush=True)
