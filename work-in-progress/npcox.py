"""numpy engine for whole finite Coxeter groups (exploration; written for
plan C1, shares no code with ../../scripts).

Models (generators as permutations of a point set):
  roots  : any type (A..I), action on the root system; roots computed exactly
           over Z[phi] (phi = golden ratio) as integer pairs (a, b) = a + b*phi.
  signed : types B, D as signed permutations on 2n points (exact keys).
Right multiplication (x g)(c) = x[g[c]].  BFS by right multiplication, so
level = Coxeter length.  Keys: exact mixed radix for 'signed', 64-bit
random-weight hash for 'roots' (the element count is compared with the known
order, which rules out collisions among the elements reached).

Provides: enumeration, right multiplication tables R, inverse, lengths,
ancestors a(w) by the prefix DP (rem, linv, uniqueness ties), left inversion
masks (bitsets over positive roots, roots model), reflections and the
Bruhat envelope G (least Bruhat-monotone G with G(w) >= 1 + G(rem w)).
"""
import math
import numpy as np
from collections import deque

ORDERS = {}
for n in range(1, 13):
    ORDERS[('A', n)] = math.factorial(n + 1)
    ORDERS[('B', n)] = 2 ** n * math.factorial(n)
    ORDERS[('D', n)] = 2 ** (n - 1) * math.factorial(n)
ORDERS.update({('E', 6): 51840, ('E', 7): 2903040, ('E', 8): 696729600, ('F', 4): 1152,
               ('H', 3): 120, ('H', 4): 14400})

# ---------------------------------------------------------------- Z[phi]
def zadd(x, y): return (x[0] + y[0], x[1] + y[1])
def zsub(x, y): return (x[0] - y[0], x[1] - y[1])
def zmul(x, y):  # phi^2 = phi + 1
    a, b = x; c, d = y
    return (a * c + b * d, a * d + b * c + b * d)
Z0 = (0, 0)


def coxeter_matrix(typ, n):
    m = [[2] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = 1
    def e(i, j, v=3):
        m[i][j] = m[j][i] = v
    if typ == 'A':
        for i in range(n - 1): e(i, i + 1)
    elif typ == 'B':
        e(0, 1, 4)
        for i in range(1, n - 1): e(i, i + 1)
    elif typ == 'D':
        e(0, 2); e(1, 2)
        for i in range(2, n - 1): e(i, i + 1)
    elif typ == 'E':        # Bourbaki: 1-3-4-5-..., 2 attached to 4 (0-based 0-2-3-4..., 1-3)
        e(0, 2); e(1, 3)
        for i in range(2, n - 1): e(i, i + 1)
    elif typ == 'F':
        e(0, 1); e(1, 2, 4); e(2, 3)
    elif typ == 'H':
        e(0, 1, 5)
        for i in range(1, n - 1): e(i, i + 1)
    elif typ == 'I':        # I2(m): n is m
        return [[1, n], [n, 1]]
    return m


def cartan(typ, n):
    """Cartan-type matrix C[i][j] = <alpha_j, alpha_i^vee> over Z[phi]."""
    M = coxeter_matrix(typ, n)
    r = len(M)
    C = [[Z0] * r for _ in range(r)]
    for i in range(r):
        C[i][i] = (2, 0)
        for j in range(r):
            if i == j: continue
            m = M[i][j]
            if m == 2: continue
            if m == 3: C[i][j] = (-1, 0)
            elif m == 5: C[i][j] = (0, -1)
            elif m == 4:      # B/F: long-short, use (-2, -1) pattern
                pass
            elif m == 6:
                pass
            else:
                raise ValueError('use dihedral_perms for m=%d' % m)
    if typ == 'B':
        C[0][1], C[1][0] = (-2, 0), (-1, 0)
    if typ == 'F':
        C[1][2], C[2][1] = (-2, 0), (-1, 0)
    return C


def root_gens(typ, n):
    C = cartan(typ, n)
    r = len(C)
    simple = [tuple((1, 0) if j == i else Z0 for j in range(r)) for i in range(r)]
    def refl(i, b):
        c = Z0
        for j in range(r):
            if b[j] != Z0 and C[i][j] != Z0:
                c = zadd(c, zmul(b[j], C[i][j]))
        bb = list(b); bb[i] = zsub(bb[i], c); return tuple(bb)
    idx = {s: k for k, s in enumerate(simple)}
    roots = list(simple)
    q = deque(simple)
    while q:
        b = q.popleft()
        for i in range(r):
            c = refl(i, b)
            if c not in idx:
                idx[c] = len(roots); roots.append(c); q.append(c)
    def positive(b):
        return all((x[0] + x[1] * 1.618033988749895) >= -1e-9 for x in b) and any(x != Z0 for x in b)
    pos = [positive(b) for b in roots]
    gens = [tuple(idx[refl(i, b)] for b in roots) for i in range(r)]
    return gens, roots, pos


def dihedral_gens(m):
    """I2(m) acting on 2m points (the regular 2m-gon's edges/roots)."""
    # roots e^{i k pi/m}, k = 0..2m-1; simple roots k=0 and k=m+1? use reflections of the m-gon
    pts = 2 * m
    # reflection across axis at angle a: k -> (2a - k) mod 2m, with roots at angles k*pi/m
    s = tuple((m - k) % pts for k in range(pts))          # root alpha_0 at angle 0 -> angle pi (negated)
    t = tuple((m + 2 - k) % pts for k in range(pts))      # second simple root at angle (m+1)pi/m
    return [s, t]


def signed_gens(typ, n):
    pts = 2 * n
    code = lambda i: n - 1 + i if i > 0 else n + i
    def mk(f):
        P = [0] * pts
        for i in range(1, n + 1):
            P[code(i)] = code(f[i]); P[code(-i)] = code(-f[i])
        return tuple(P)
    g = []
    f = {i: i for i in range(1, n + 1)}
    if typ == 'B': f[1] = -1
    else: f[1], f[2] = -2, -1
    g.append(mk(f))
    for i in range(1, n):
        f = {k: k for k in range(1, n + 1)}; f[i], f[i + 1] = i + 1, i; g.append(mk(f))
    return g


class NpGroup:
    def __init__(self, typ, n, model='roots', verbose=False):
        self.typ, self.n, self.model = typ, n, model
        self.pos = None
        if model == 'signed':
            g = signed_gens(typ, n)
        elif typ == 'I':
            g = dihedral_gens(n)
            # positive roots for the simple roots at angles 0 and pi + pi/m: k = 0, m+1, ..., 2m-1
            self.pos = [k == 0 or k > n for k in range(2 * n)]
        else:
            g, self.roots, self.pos = root_gens(typ, n)
        self.gens = np.array(g, dtype=np.int16)
        self.r, self.deg = self.gens.shape
        dt = np.uint8 if self.deg < 256 else np.uint16
        self.dt = dt
        if model == 'signed':
            base = np.int64(self.deg)
            cols = np.arange(n, 2 * n)       # images of the positive points determine w
            w = base ** np.arange(n, dtype=np.int64)
            self.keyof = lambda P: (P[:, cols].astype(np.int64) * w[None, :]).sum(axis=1)
        else:
            rng = np.random.default_rng(271828)
            mult = rng.integers(1, 2 ** 62, size=self.deg, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
            self.keyof = lambda P: (P.astype(np.uint64) * mult[None, :]).sum(axis=1, dtype=np.uint64)
        self._enumerate(verbose)

    def _enumerate(self, verbose):
        gens = self.gens.astype(np.int64)
        F = np.arange(self.deg, dtype=self.dt)[None, :]
        levels = [F]
        keys = [self.keyof(F)]
        prevK = np.empty(0, keys[0].dtype)
        curK = keys[0]
        while True:
            cand = np.concatenate([F[:, gens[s]] for s in range(self.r)])
            ck = self.keyof(cand)
            u, first = np.unique(ck, return_index=True)
            keep = ~np.isin(u, prevK) & ~np.isin(u, curK)
            if not keep.any():
                break
            F = cand[first[keep]]
            prevK, curK = curK, u[keep]
            levels.append(F); keys.append(curK)
            if verbose: print('level', len(levels) - 1, len(F), flush=True)
        self.E = np.concatenate(levels)
        self.order = len(self.E)
        exp = ORDERS.get((self.typ, self.n))
        if self.typ == 'I': exp = 2 * self.n
        assert exp is None or exp == self.order, (self.order, exp)
        self.length = np.concatenate([np.full(len(l), k, np.int16) for k, l in enumerate(levels)])
        self.offs = np.cumsum([0] + [len(l) for l in levels])
        K = np.concatenate(keys)
        self.sortk = np.argsort(K)
        self.K = K[self.sortk]
        assert len(np.unique(self.K)) == self.order
        self.R = [self.lookup(self.E[:, gens[s]]) for s in range(self.r)]
        Einv = np.empty_like(self.E)
        Einv[np.arange(self.order)[:, None], self.E.astype(np.int64)] = np.arange(self.deg, dtype=self.dt)[None, :]
        self.inv = self.lookup(Einv)
        del Einv

    def lookup(self, P, chunk=2000000):
        out = np.empty(len(P), np.int32)
        for lo in range(0, len(P), chunk):
            k = self.keyof(P[lo:lo + chunk])
            p = np.searchsorted(self.K, k)
            p[p >= len(self.K)] = 0
            assert np.all(self.K[p] == k), 'element not found'
            out[lo:lo + chunk] = self.sortk[p]
        return out

    def ancestors(self):
        O, r = self.order, self.r
        alen = np.zeros(O, np.int16); rem = np.zeros(O, np.int32); linv = np.zeros(O, np.int16)
        isinv = self.inv == np.arange(O)
        ties = 0
        for L in range(1, len(self.offs) - 1):
            lo, hi = self.offs[L], self.offs[L + 1]
            m = hi - lo
            cl = np.full((r, m), -1, np.int16); cr = np.zeros((r, m), np.int32)
            for s in range(r):
                k = self.R[s][lo:hi]
                d = self.length[k] < L
                kd = k[d]
                cl[s, d] = alen[kd]
                cr[s, d] = self.R[s][rem[kd]]
            bl = cl.max(axis=0); best = cl.argmax(axis=0); br = cr[best, np.arange(m)]
            ii = isinv[lo:hi]
            ties += int(((cl == bl[None, :]) & (cr != br[None, :])).any(axis=0)[~ii].sum())
            alen[lo:hi] = np.where(ii, L, bl)
            rem[lo:hi] = np.where(ii, 0, br)
            linv[lo:hi] = 1 + linv[rem[lo:hi]]
        self.alen, self.rem, self.linv, self.ties, self.isinv = alen, rem, linv, ties, isinv
        return linv

    def word(self, i):
        out = []
        while i:
            for s in range(self.r):
                k = int(self.R[s][i])
                if self.length[k] < self.length[i]:
                    out.append(s); i = k; break
        return out[::-1]

    def reflections(self):
        # t = x s x^-1 for all x, s : as permutations
        T = set()
        Einv = self.E[self.inv]
        for s in range(self.r):
            # x s x^-1 (c) = x[ s[ xinv[c] ] ]
            g = self.gens[s].astype(np.int64)
            P = np.take_along_axis(self.E, g[Einv.astype(np.int64)], axis=1)
            T |= set(self.lookup(P).tolist())
        self.T = sorted(T)
        return self.T

    def envelope_G(self):
        """Least G: Bruhat-monotone, G(e)=0, G(w) >= 1 + G(rem w)."""
        if not hasattr(self, 'T'): self.reflections()
        G = np.zeros(self.order, np.int16)
        Tp = [self.E[t].astype(np.int64) for t in self.T]
        for L in range(1, len(self.offs) - 1):
            lo, hi = self.offs[L], self.offs[L + 1]
            X = self.E[lo:hi]
            g = 1 + G[self.rem[lo:hi]]
            for t in Tp:
                k = self.lookup(X[:, t])
                low = self.length[k] < L
                g = np.where(low, np.maximum(g, G[k]), g)
            G[lo:hi] = g
        self.G = G
        return G

    def left_inv_masks(self):
        """I_L(w) = {beta > 0 : w^-1 beta < 0} as Python-int-free uint64 words."""
        assert self.pos is not None
        posidx = np.flatnonzero(self.pos)
        self.posidx = posidx
        isneg = ~np.array(self.pos)
        nw = (len(posidx) + 63) // 64
        M = np.zeros((self.order, nw), np.uint64)
        Einv = self.E[self.inv]
        for b, beta in enumerate(posidx):
            neg = isneg[Einv[:, beta]]
            M[neg, b // 64] |= np.uint64(1) << np.uint64(b % 64)
        self.masks = M
        return M
