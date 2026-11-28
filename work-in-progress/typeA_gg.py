"""Type A, implementation 2: the ancestor computed one element at a time by
Gaetz-Gao's construction (arXiv:2610.22857, proofs of Lemma 3.1 and Thm 1.1):

  a(w) = lr(w^-1, w),
  lb(x, y): y = s_1...s_m reduced, u_m = x, u_{i-1} = delta_{s_i}(u_i),
            delta_s(g) = min_Bruhat(g, gs);  lb = u_0^-1 x
  br(x, y) = lb(y^-1, x^-1)^-1
  lr(x, y): z_0 = x, z_{2i+1} = br(z_{2i}, y), z_{2i+2} = lb(x, z_{2i+1}),
            until the sequence is stationary.

Vectorised over a batch of permutations of the same N; no tables, no shared
code with typeA_dp.py.  Conventions as in typeA_dp.py (0-based one-line,
(xy)(i) = x[y[i]], x s_i swaps positions i, i+1).

Also tests, for every w in the batch (pure combinatorics, no a(w) needed):
  (K)  Key Lemma: for every p with w(p) = p + D, D = disp(w), the explicit
       involution iota_p built by the walk construction (LAB.md, Lemma K)
       is an involution, is a prefix of w, and moves v = w(p);
       and the same for drops (w(p) = p - D) via phi(w) = w0 w w0.
and, using a(w):
  (D)  disp(a(w) w) < disp(w);   (I) l_inv(w) <= disp(w) (iterating)
  (B)  iota_p <=_Bruhat a(w)  (Gaetz-Gao Cor. 1.2; tableau criterion)

Usage:
  python typeA_gg.py exhaustive N        all of S_N; writes out/A{N}_gg_anc.npy
  python typeA_gg.py random N COUNT SEED [band]   random (or banded) permutations
"""
import sys, os, time, json, math, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def inverse(P):
    M, N = P.shape
    Q = np.empty_like(P)
    Q[np.arange(M)[:, None], P] = np.arange(N, dtype=P.dtype)[None, :]
    return Q


def compose(X, Y):          # (XY)(i) = X[Y[i]]
    return np.take_along_axis(X, Y.astype(np.int64), axis=1)


def length(P):
    N = P.shape[1]
    out = np.zeros(len(P), np.int64)
    for i in range(N):
        out += (P[:, i + 1:] < P[:, i:i + 1]).sum(axis=1)
    return out


def reduced_words(Y):
    """Return W (M x Lmax) with -1 padding, such that row i spells a reduced
    word s_{W[i,0]} s_{W[i,1]} ... for Y[i]."""
    M, N = Y.shape
    Lmax = int(length(Y).max()) if M else 0
    W = np.full((M, max(Lmax, 1)), -1, np.int16)
    Z = Y.copy()
    letters = []
    rows = np.arange(M)
    for k in range(Lmax):
        d = Z[:, :-1] > Z[:, 1:]
        has = d.any(axis=1)
        i = d.argmax(axis=1)
        r = rows[has]
        ii = i[has]
        a = Z[r, ii].copy()
        Z[r, ii] = Z[r, ii + 1]
        Z[r, ii + 1] = a
        col = np.full(M, -1, np.int16)
        col[has] = ii
        letters.append(col)
    # Y s_{d1} s_{d2} ... s_{dm} = e  =>  Y = s_{dm} ... s_{d1}
    lens = length(Y)
    L =np.array(letters).T if letters else np.zeros((M, 0), np.int16)   # M x Lmax, col k = d_{k+1}
    # reverse the non-padded part of each row
    idx = lens[:, None] - 1 - np.arange(max(Lmax, 1))[None, :]
    valid = idx >= 0
    idxc = np.where(valid, idx, 0)
    if Lmax:
        W = np.where(valid, np.take_along_axis(L, idxc, axis=1), -1).astype(np.int16)
    return W


def lb(X, Y):
    """lb(x, y) for each row."""
    W = reduced_words(Y)
    U = X.copy()
    M = len(U)
    rows = np.arange(M)
    for k in range(W.shape[1] - 1, -1, -1):
        s = W[:, k]
        act = s >= 0
        sc = np.where(act, s, 0).astype(np.int64)
        a = U[rows, sc]
        b = U[rows, sc + 1]
        sw = act & (a > b)              # s right descent of U: delta_s(U) = U s
        r = rows[sw]
        U[r, sc[sw]] = b[sw]
        U[r, sc[sw] + 1] = a[sw]
    return compose(inverse(U), X)


def br(X, Y):
    return inverse(lb(inverse(Y), inverse(X)))


def lr(X, Y):
    Z = X.copy()
    while True:
        Z1 = br(Z, Y)
        Z2 = lb(X, Z1)
        if np.array_equal(Z2, Z) and np.array_equal(Z1, Z):
            return Z
        Z = Z2


def ancestor(P):
    return lr(inverse(P), P)


def disp(P):
    N = P.shape[1]
    return np.abs(P.astype(np.int64) - np.arange(N)[None, :]).max(axis=1)


# ---------------------------------------------------------------- Lemma K
def iota_exc(w, p):
    """w: tuple (0-based), w[p] = p + D with D = disp(w).  Walk construction:
    L = w([0..p]) (values left of the cut after position p), R = [0..v] \\ L.
    Walk downwards from v; J = [c0, v] at the first return of #L - #R to 0;
    iota swaps the i-th smallest L with the i-th smallest R inside J."""
    N = len(w)
    v = w[p]
    Lset = set(w[:p + 1])
    h = 0
    c0 = None
    for c in range(v, -1, -1):
        h += 1 if c in Lset else -1
        if h == 0:
            c0 = c
            break
    assert c0 is not None, ('walk never returns', w, p)
    Ls = [c for c in range(c0, v + 1) if c in Lset]
    Rs = [c for c in range(c0, v + 1) if c not in Lset]
    io = list(range(N))
    for l, r in zip(Ls, Rs):
        io[l], io[r] = r, l
    return io


def inv_count(x):
    n = len(x)
    return sum(1 for i in range(n) for j in range(i + 1, n) if x[i] > x[j])


def check_K(w):
    """Returns (number of extremal positions checked, list of failures)."""
    N = len(w)
    D = max(abs(w[i] - i) for i in range(N))
    if D == 0:
        return 0, []
    fails = []
    cnt = 0
    lw = inv_count(w)
    w0 = lambda i: N - 1 - i
    phiw = tuple(w0(w[w0(i)]) for i in range(N))
    for p in range(N):
        for kind in ('exc', 'drop'):
            if kind == 'exc' and w[p] == p + D:
                io = iota_exc(w, p)
                v = w[p]
            elif kind == 'drop' and w[p] == p - D:
                q = w0(p)                         # phi(w)(q) = q + D
                io2 = iota_exc(phiw, q)
                io = [w0(io2[w0(i)]) for i in range(N)]    # phi(iota)
                v = w[p]
            else:
                continue
            cnt += 1
            ok_inv = all(io[io[i]] == i for i in range(N))
            iw = tuple(io[w[i]] for i in range(N))          # iota^-1 w = iota w
            ok_pre = inv_count(io) + inv_count(iw) == lw
            ok_mov = io[v] != v
            if not (ok_inv and ok_pre and ok_mov):
                fails.append((w, p, kind, tuple(io), ok_inv, ok_pre, ok_mov))
    return cnt, fails


def bruhat_leq(u, w):
    """Tableau / rank-matrix criterion (Bjorner-Brenti Thm 2.1.5)."""
    N = len(u)
    for i in range(N):
        a = sorted(u[:i + 1])
        b = sorted(w[:i + 1])
        if any(x > y for x, y in zip(a, b)):
            return False
    return True


def run(P, tag, kcheck=True, bcheck=False):
    t0 = time.time()
    M, N = P.shape
    out = {}
    A = np.empty_like(P)
    CH = 200000
    for lo in range(0, M, CH):
        A[lo:lo + CH] = ancestor(P[lo:lo + CH])
    t1 = time.time()
    nz = ~(P == np.arange(N, dtype=P.dtype)[None, :]).all(axis=1)
    AA = compose(A, A)
    out['a_not_involution'] = int(((AA != np.arange(N, dtype=P.dtype)[None, :]).any(axis=1) & nz).sum())
    Rm = compose(inverse(A), P)
    out['a_not_prefix'] = int((length(A) + length(Rm) != length(P)).sum())
    dP, dR = disp(P), disp(Rm)
    out['disp_not_decreasing'] = int(((dR >= dP) & nz).sum())
    # iterate to get l_inv
    linv = np.zeros(M, np.int64)
    cur = P.copy()
    alive = nz.copy()
    while alive.any():
        linv[alive] += 1
        idx = np.flatnonzero(alive)
        sub = cur[idx]
        a = np.empty_like(sub)
        for lo in range(0, len(sub), CH):
            a[lo:lo + CH] = ancestor(sub[lo:lo + CH])
        cur[idx] = compose(inverse(a), sub)
        alive = ~(cur == np.arange(N, dtype=P.dtype)[None, :]).all(axis=1)
    out['max_linv'] = int(linv.max())
    out['linv_gt_disp'] = int((linv > dP).sum())
    out['linv_gt_rank'] = int((linv > N - 1).sum())
    out['hist_linv'] = np.bincount(linv).tolist()
    t2 = time.time()
    if kcheck:
        cnt = 0
        fails = []
        bfail = 0
        for i in range(M):
            w = tuple(int(x) for x in P[i])
            c, f = check_K(w)
            cnt += c
            fails += f
            if bcheck and c:
                a = tuple(int(x) for x in A[i])
                N_ = len(w)
                D = max(abs(w[j] - j) for j in range(N_))
                for p in range(N_):
                    if w[p] == p + D:
                        if not bruhat_leq(iota_exc(w, p), a):
                            bfail += 1
        out['K_extremal_positions_checked'] = cnt
        out['K_failures'] = len(fails)
        out['K_first_failures'] = [str(f) for f in fails[:3]]
        if bcheck:
            out['iota_not_bruhat_below_a'] = bfail
    out.update(tag=tag, N=N, count=M, sec_anc=round(t1 - t0, 1), sec_total=round(time.time() - t0, 1),
               numpy=np.__version__, python=sys.version.split()[0])
    return A, out


if __name__ == '__main__':
    mode = sys.argv[1]
    N = int(sys.argv[2])
    if mode == 'exhaustive':
        P = np.array(list(itertools.permutations(range(N))), dtype=np.int8)
        kc = (sys.argv[3] != 'noK') if len(sys.argv) > 3 else True
        A, out = run(P, 'exhaustive', kcheck=kc, bcheck=(N <= 7))
        # Lehmer rank of a(w), to compare with typeA_dp.py
        r = np.zeros(len(A), np.int64)
        for i in range(N):
            r += (A[:, i + 1:] < A[:, i:i + 1]).sum(axis=1) * math.factorial(N - 1 - i)
        np.save(os.path.join(HERE, 'out', 'A%d_gg_anc.npy' % N), r)
        print(json.dumps(out), flush=True)
    else:
        cnt, seed = int(sys.argv[3]), int(sys.argv[4])
        band = int(sys.argv[5]) if len(sys.argv) > 5 else 0
        rng = np.random.default_rng(seed)
        if band:
            # banded permutations: products of random adjacent transpositions
            # with displacement capped at `band` (rejection of moves that exceed it)
            P = np.tile(np.arange(N, dtype=np.int16), (cnt, 1))
            rows = np.arange(cnt)
            for _ in range(6 * N * N):
                s = rng.integers(0, N - 1, size=cnt)
                Q = P.copy()
                a = Q[rows, s].copy()
                Q[rows, s] = Q[rows, s + 1]
                Q[rows, s + 1] = a
                ok = disp(Q) <= band
                P[ok] = Q[ok]
        else:
            P = np.array([rng.permutation(N) for _ in range(cnt)], dtype=np.int16)
        A, out = run(P, 'random' + ('_band%d' % band if band else ''), kcheck=True, bcheck=False)
        out['seed'] = seed
        out['max_disp'] = int(disp(P).max())
        print(json.dumps(out), flush=True)
