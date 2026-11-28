"""Type A, implementation 1: exact dynamic programme over all of S_N.

Elements are indexed by their Lehmer rank (exact bijection, no hashing).
Conventions: permutations are 0-based one-line arrays x[0..N-1];
(xy)(i) = x[y[i]];  s_i swaps i, i+1;  right multiplication x s_i swaps the
entries in positions i, i+1;  l(x) = number of inversions.
u <=_R w  iff  l(w) = l(u) + l(u^-1 w)   (u is a prefix of w).

a(w) = longest involution prefix.  DP: if w is not an involution,
a(w) is a prefix of ws for some right descent s, so a(w) is the longest of
the a(ws), and rem(w) := a(w)^-1 w = rem(ws) s.  Ties with different
remainders are counted (Gaetz-Gao: there are none).

Checks, for every w != e:
  (D)  disp(rem w) < disp(w), disp(x) = max_i |x(i) - i|
  (I)  l_inv(w) <= disp(w) <= N-1
  (P)  a(w) := w rem(w)^-1 is an involution and l(a) + l(rem) = l(w)
Writes  out/A{N}_dp_anc.npy  (Lehmer rank of a(w), indexed by rank of w).

Usage: python typeA_dp.py N   (S_N, i.e. type A_{N-1})
"""
import sys, os, time, json, itertools, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def lehmer_rank(P):
    M, N = P.shape
    r = np.zeros(M, np.int64)
    for i in range(N):
        c = (P[:, i + 1:] < P[:, i:i + 1]).sum(axis=1)
        r += c * math.factorial(N - 1 - i)
    return r


def main(N):
    t0 = time.time()
    P = np.array(list(itertools.permutations(range(N))), dtype=np.int8)   # lex order = rank order
    M = len(P)
    assert M == math.factorial(N)
    assert np.array_equal(lehmer_rank(P[:: max(1, M // 1000)]), np.arange(0, M, max(1, M // 1000)))
    length = np.zeros(M, np.int16)
    for i in range(N):
        length += (P[:, i + 1:] < P[:, i:i + 1]).sum(axis=1).astype(np.int16)
    inv = np.empty_like(P)
    inv[np.arange(M)[:, None], P] = np.arange(N, dtype=np.int8)[None, :]
    invrank = lehmer_rank(inv)
    isinv = invrank == np.arange(M)
    R = []
    for s in range(N - 1):
        Q = P.copy()
        Q[:, [s, s + 1]] = Q[:, [s + 1, s]]
        R.append(lehmer_rank(Q).astype(np.int32))
    disp = np.abs(P - np.arange(N, dtype=np.int8)[None, :]).max(axis=1).astype(np.int16)
    t1 = time.time()
    order = np.argsort(length, kind='stable')
    lev = np.searchsorted(length[order], np.arange(length.max() + 2))
    alen = np.zeros(M, np.int16)
    rem = np.zeros(M, np.int64)
    linv = np.zeros(M, np.int16)
    ties = 0
    for L in range(1, length.max() + 1):
        I = order[lev[L]:lev[L + 1]]
        m = len(I)
        cl = np.full((N - 1, m), -1, np.int16)
        cr = np.zeros((N - 1, m), np.int64)
        for s in range(N - 1):
            k = R[s][I]
            d = P[I, s] > P[I, s + 1]                # right descent
            kd = k[d]
            cl[s, d] = alen[kd]
            cr[s, d] = R[s][rem[kd]]
        bl = cl.max(axis=0)
        best = cl.argmax(axis=0)
        br = cr[best, np.arange(m)]
        ii = isinv[I]
        ties += int(((cl == bl[None, :]) & (cr != br[None, :])).any(axis=0)[~ii].sum())
        alen[I] = np.where(ii, L, bl)
        rem[I] = np.where(ii, 0, br)
        linv[I] = 1 + linv[rem[I]]
    t2 = time.time()
    nz = np.arange(1, M)
    badD = int((disp[rem[nz]] >= disp[nz]).sum())
    badI = int((linv > disp).sum())
    # (P): a = w rem^-1, as permutations
    Rm = P[rem]
    Rinv = np.empty_like(Rm)
    Rinv[np.arange(M)[:, None], Rm] = np.arange(N, dtype=np.int8)[None, :]
    A = np.take_along_axis(P, Rinv.astype(np.int64), axis=1)      # A[i] = P[Rinv[i]]  -> (w rem^-1)(i)
    AA = np.take_along_axis(A, A.astype(np.int64), axis=1)
    notinv = int((AA != np.arange(N, dtype=np.int8)[None, :]).any(axis=1)[nz].sum())
    arank = lehmer_rank(A)
    badP = int((length[arank] + length[rem] != length).sum())
    badAlen = int((length[arank] != alen).sum())
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    np.save(os.path.join(HERE, 'out', 'A%d_dp_anc.npy' % N), arank.astype(np.int64))
    res = dict(N=N, type='A%d' % (N - 1), order=M, involutions_incl_e=int(isinv.sum()),
               ancestor_ties=ties, max_linv=int(linv.max()), max_disp=int(disp.max()),
               hist_linv=np.bincount(linv).tolist(),
               disp_not_decreasing=badD, linv_gt_disp=badI, a_not_involution=notinv,
               a_not_prefix=badP, alen_mismatch=badAlen,
               linv_eq_disp=int((linv == disp).sum()),
               sec_tables=round(t1 - t0, 1), sec_total=round(time.time() - t0, 1),
               numpy=np.__version__, python=sys.version.split()[0])
    print(json.dumps(res), flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))
