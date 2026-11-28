"""Vectorised exhaustive test of Lemma K (type A), independent of any
ancestor computation.

Lemma K.  Let w in S_N, D = disp(w) >= 1, and p with w(p) = p + D (0-based
positions and values), v = w(p).  Put L = w({0..p}), R = {0..v} \\ L.
Walk down from c = v: h(c) = #(L & [c,v]) - #(R & [c,v]).  Let c0 be the
largest c with h(c) = 0, J = [c0, v].  Let iota swap the i-th smallest
element of L & J with the i-th smallest element of R & J (fixing all other
values).  Then iota is an involution, iota <=_R w, and iota(v) != v.

For drops (w(p) = p - D) the statement is applied to phi(w) = w0 w w0.
Here we check, for every w in S_N and every extremal excedance AND every
extremal drop (drops handled directly by phi), that:
  (k1) the walk returns to 0,
  (k2) l(iota) + l(iota w) = l(w)    [prefix, by definition of <=_R]
  (k3) iota(v) != v
  (k4) the inversions of iota are all of the form (r, l), r < l, r in R, l in L
       (the structural claim used in the proof)
Usage: python typeA_K_vec.py N
"""
import sys, time, json, itertools
import numpy as np


def length(P):
    N = P.shape[1]
    out = np.zeros(len(P), np.int64)
    for i in range(N):
        out += (P[:, i + 1:] < P[:, i:i + 1]).sum(axis=1)
    return out


def check(P):
    """P: (M x N) permutations, all with an extremal excedance at the
    position given in the companion array; returns failure counts."""
    pass


def run(N, CH=400000):
    t0 = time.time()
    allP = np.array(list(itertools.permutations(range(N))), dtype=np.int8)
    tot = dict(pairs=0, k1=0, k2=0, k3=0, k4=0)
    ar = np.arange(N)
    for lo in range(0, len(allP), CH):
        P0 = allP[lo:lo + CH].astype(np.int64)
        for kind in ('exc', 'drop'):
            if kind == 'exc':
                P = P0
            else:  # phi(w)(i) = N-1-w(N-1-i)
                P = (N - 1 - P0[:, ::-1])
            D = np.abs(P - ar[None, :]).max(axis=1)
            ext = (P - ar[None, :]) == D[:, None]
            ext &= (D[:, None] > 0)
            wi, pi = np.nonzero(ext)
            W = P[wi]
            m = len(wi)
            tot['pairs'] += m
            v = W[np.arange(m), pi]
            pos = np.empty_like(W)
            pos[np.arange(m)[:, None], W] = ar[None, :]      # pos[c] = position of value c
            inL = pos <= pi[:, None]                           # value c in L
            le_v = ar[None, :] <= v[:, None]
            step = np.where(inL, 1, -1) * le_v                 # contributions for c <= v
            # h(c) = sum_{c' = c..v} step(c')  -> reverse cumulative sum
            h = np.cumsum(step[:, ::-1], axis=1)[:, ::-1]
            zero = (h == 0) & le_v
            has = zero.any(axis=1)
            tot['k1'] += int((~has).sum())
            # largest c with h(c) = 0
            c0 = N - 1 - np.argmax(zero[:, ::-1], axis=1)
            inJ = (ar[None, :] >= c0[:, None]) & le_v
            Lj = inJ & inL
            Rj = inJ & ~inL
            rkL = np.cumsum(Lj, axis=1) * Lj                    # rank (1-based) within L&J
            rkR = np.cumsum(Rj, axis=1) * Rj
            iota = np.tile(ar, (m, 1))
            # pair i-th L with i-th R
            kmax = int(rkL.max()) if m else 0
            for k in range(1, kmax + 1):
                lk = (rkL == k)
                rk = (rkR == k)
                hasl = lk.any(axis=1)
                lval = np.argmax(lk, axis=1)
                rval = np.argmax(rk, axis=1)
                assert np.array_equal(hasl, rk.any(axis=1))
                rows = np.flatnonzero(hasl)
                iota[rows, lval[rows]] = rval[rows]
                iota[rows, rval[rows]] = lval[rows]
            IW = np.take_along_axis(iota, W, axis=1)          # (iota w)(i) = iota(w(i))
            k2 = length(iota) + length(IW) != length(W)
            tot['k2'] += int(k2.sum())
            tot['k3'] += int((iota[np.arange(m), v] == v).sum())
            # (k4) every inversion (a<b, iota(a) > iota(b)) has a in R, b in L
            bad4 = np.zeros(m, bool)
            for a in range(N):
                for b in range(a + 1, N):
                    invab = iota[:, a] > iota[:, b]
                    ok = (~inL[:, a]) & le_v[:, a] & inL[:, b]
                    bad4 |= invab & ~ok
            tot['k4'] += int(bad4.sum())
    tot.update(N=N, seconds=round(time.time() - t0, 1), numpy=np.__version__)
    return tot


if __name__ == '__main__':
    print(json.dumps(run(int(sys.argv[1]))), flush=True)
