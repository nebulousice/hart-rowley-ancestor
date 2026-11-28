"""Small exact engine for finite Coxeter groups given by generator
permutations (exploration; up to ~10^6 elements).  Pure python + numpy.

Elements are permutations (tuples) of a point set; enumerated by BFS in
right multiplication, so the BFS level is the Coxeter length (exact).
Prefix order u <=_R w iff l(w) = l(u) + l(u^-1 w).

Models
  'A', n : S_{n+1} on {0..n}
  'B', n : signed permutations of {±1..±n}; points coded  +i -> n-1+i, -i -> n-i
           (so the points are ordered -n < ... < -1 < 1 < ... < n as 0..2n-1);
           s0 = sign change in position 1, s_i swaps positions i, i+1
  'D', n : s0' = swap positions 1,2 with both signs changed, s_i as in B
A signed permutation w is stored as its action on the 2n points:
   P[c] = code(w(point c)), with w(-i) = -w(i).  Right multiplication by a
generator acts on positions, as for permutations: (x g)(c) = x(g(c)).
"""
import numpy as np
from collections import deque


def code(n, i):          # signed integer i in ±[n] -> point 0..2n-1
    return n - 1 + i if i > 0 else n + i


def decode(n, c):
    return c - n + 1 if c >= n else c - n


def gens(typ, n):
    if typ == 'A':
        out = []
        for i in range(n):
            g = list(range(n + 1)); g[i], g[i + 1] = g[i + 1], g[i]; out.append(tuple(g))
        return out
    m = 2 * n
    def from_signed(f):     # f: dict i -> f(i) for i = 1..n
        P = [0] * m
        for i in range(1, n + 1):
            P[code(n, i)] = code(n, f[i]); P[code(n, -i)] = code(n, -f[i])
        return tuple(P)
    out = []
    if typ == 'B':
        f = {i: i for i in range(1, n + 1)}; f[1] = -1; out.append(from_signed(f))
    elif typ == 'D':
        f = {i: i for i in range(1, n + 1)}; f[1] = -2; f[2] = -1; out.append(from_signed(f))
    for i in range(1, n):
        f = {k: k for k in range(1, n + 1)}; f[i] = i + 1; f[i + 1] = i; out.append(from_signed(f))
    return out


class Group:
    def __init__(self, typ, n, tables=True):
        self.typ, self.n = typ, n
        self.gens = gens(typ, n)
        self.r = len(self.gens)
        e = tuple(range(len(self.gens[0])))
        self.elems = [e]
        self.index = {e: 0}
        self.length = [0]
        q = deque([0])
        while q:
            i = q.popleft()
            x = self.elems[i]
            for g in self.gens:
                y = tuple(x[c] for c in g)
                if y not in self.index:
                    self.index[y] = len(self.elems)
                    self.elems.append(y)
                    self.length.append(self.length[i] + 1)
                    q.append(self.index[y])
        self.order = len(self.elems)
        self.E = np.array(self.elems, dtype=np.int16)
        self.L = np.array(self.length)
        if tables:
            self.R = [[self.index[tuple(x[c] for c in g)] for x in self.elems] for g in self.gens]
            self.inv = [self.index[tuple(np.argsort(x).tolist())] for x in self.elems]

    def mul(self, i, j):
        x, y = self.elems[i], self.elems[j]
        return self.index[tuple(x[c] for c in y)]

    def ancestors(self):
        """a(w) via the DP; returns (alen, rem, linv), and counts ties."""
        O = self.order
        alen = [0] * O; rem = [0] * O; linv = [0] * O
        ties = 0
        for w in range(1, O):          # BFS order is by length
            if self.inv[w] == w:
                alen[w], rem[w] = self.length[w], 0
            else:
                best, br = -1, None
                cands = []
                for s in range(self.r):
                    k = self.R[s][w]
                    if self.length[k] < self.length[w]:
                        cands.append((alen[k], self.R[s][rem[k]]))
                best = max(c[0] for c in cands)
                rs = set(c[1] for c in cands if c[0] == best)
                ties += len(rs) > 1
                alen[w], rem[w] = best, rs.pop()
            linv[w] = 1 + linv[rem[w]]
        self.alen, self.rem, self.linv, self.ties = alen, rem, linv, ties
        return alen, rem, linv

    def reflections(self):
        T = set()
        for x in range(self.order):
            for s in range(self.r):
                T.add(self.mul(self.mul(x, self.R[s][0]), self.inv[x]))
        self.T = sorted(T)
        return self.T

    def signed(self, i):
        """Signed one-line notation [w(1), ..., w(n)] (types B, D)."""
        n = self.n
        x = self.elems[i]
        return tuple(decode(n, x[code(n, k)]) for k in range(1, n + 1))

    def descents(self, i):
        DR = [s for s in range(self.r) if self.length[self.R[s][i]] < self.length[i]]
        j = self.inv[i]
        DL = [s for s in range(self.r) if self.length[self.R[s][j]] < self.length[j]]
        return DL, DR

    def envelope_G(self):
        """Least G with G(e)=0, Bruhat-monotone, G(w) >= 1 + G(rem w)."""
        if not hasattr(self, 'T'):
            self.reflections()
        O = self.order
        G = [0] * O
        Tel = [self.elems[t] for t in self.T]
        for w in range(1, O):
            x = self.elems[w]
            g = 1 + G[self.rem[w]]
            lw = self.length[w]
            for t in Tel:
                k = self.index[tuple(x[c] for c in t)]
                if self.length[k] < lw and G[k] > g:
                    g = G[k]
            G[w] = g
        self.G = G
        return G
