"""Finite Coxeter groups as exact permutation groups, and the Hart-Rowley
ancestor decomposition (arXiv:2502.00777).

Conventions (fixed here, used everywhere in this target):
  * (W, S) finite Coxeter system, S = simple reflections, l = Coxeter length.
  * u is a PREFIX of w  iff  l(w) = l(u) + l(u^-1 w)   (u <=_R w, right weak
    order in the Bjorner-Brenti / Gaetz-Gao naming; Hart-Rowley call the same
    order "weak (left) Bruhat order").
  * a(w) = the involution prefix of w of maximal length (unique: Gaetz-Gao,
    arXiv:2610.22857, Cor. 1.2).  a(w) is defined for w != e.
  * ancestor decomposition w = a_1 a_2 ... a_k, a_i = a(a_i ... a_k);
    involution length l_inv(w) = k, l_inv(e) = 0.

Elements are permutations of a finite set on which W acts faithfully
(the roots, or a classical model), stored as bytes/tuples.  Everything is
exact integer arithmetic (Z[phi] for H3/H4).  No floating point.
"""
from collections import deque

# ----------------------------------------------------------------------------
# Z[phi] arithmetic, phi^2 = phi + 1; elements are pairs (a, b) = a + b*phi
# ----------------------------------------------------------------------------
def zmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)


def zsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


ZERO = (0, 0)


def cartan(typ, n):
    """Cartan-type matrix C with s_i(alpha_j) = alpha_j - C[i][j] alpha_i,
    entries in Z[phi] as pairs.  Node numbering: a path 0-1-...-(n-1) for
    A, B, F, H; see the code for D and E."""
    C = [[(2, 0) if i == j else ZERO for j in range(n)] for i in range(n)]

    def edge(i, j, cij=(-1, 0), cji=(-1, 0)):
        C[i][j] = cij
        C[j][i] = cji

    if typ == 'A':
        for i in range(n - 1):
            edge(i, i + 1)
    elif typ == 'B':
        # node 0 is the special one:  0 =4= 1 - 2 - ... - (n-1)
        for i in range(1, n - 1):
            edge(i, i + 1)
        if n >= 2:
            edge(0, 1, (-2, 0), (-1, 0))
    elif typ == 'D':
        # nodes 0 and 1 are the two fork leaves, both joined to node 2
        assert n >= 3
        edge(0, 2)
        edge(1, 2)
        for i in range(2, n - 1):
            edge(i, i + 1)
    elif typ == 'E':
        # path 0-2-3-4-...-(n-1), node 1 attached to node 3 (Bourbaki, 0-based)
        assert n in (6, 7, 8)
        edge(0, 2)
        edge(1, 3)
        for i in range(2, n - 1):
            edge(i, i + 1)
    elif typ == 'F':
        assert n == 4
        edge(0, 1)
        edge(1, 2, (-2, 0), (-1, 0))
        edge(2, 3)
    elif typ == 'G':
        assert n == 2
        edge(0, 1, (-3, 0), (-1, 0))
    elif typ == 'H':
        assert n in (2, 3, 4)
        edge(0, 1, (0, -1), (0, -1))      # m = 5 :  -phi
        for i in range(1, n - 1):
            edge(i, i + 1)
    else:
        raise ValueError(typ)
    return C


def coxeter_graph(typ, n):
    """Edges {i,j} of the Coxeter graph (non-commuting pairs)."""
    if typ == 'I':
        return [(0, 1)]
    C = cartan(typ, n)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if C[i][j] != ZERO]


def root_perms(typ, n):
    """Generators as permutations (tuples) of the root system."""
    C = cartan(typ, n)
    simple = [tuple((1, 0) if j == i else ZERO for j in range(n)) for i in range(n)]

    def refl(i, beta):
        c = ZERO
        for j in range(n):
            if beta[j] != ZERO and C[i][j] != ZERO:
                p = zmul(beta[j], C[i][j])
                c = (c[0] + p[0], c[1] + p[1])
        b = list(beta)
        b[i] = zsub(b[i], c)
        return tuple(b)

    roots = {r: k for k, r in enumerate(simple)}
    order = list(simple)
    q = deque(simple)
    while q:
        b = q.popleft()
        for i in range(n):
            c = refl(i, b)
            if c not in roots:
                roots[c] = len(order)
                order.append(c)
                q.append(c)
                if len(order) > 100000:
                    raise RuntimeError('root system not finite?')
    gens = [tuple(roots[refl(i, b)] for b in order) for i in range(n)]
    return gens


def classical_perms(typ, n):
    """Second, independent model for types A, B, D: (signed) permutations.
    Same node numbering as cartan().  Points: A_n acts on {0..n};
    B_n, D_n act on {+1..+n, -1..-n} coded as i-1 and n+i-1."""
    if typ == 'A':
        gens = []
        for i in range(n):
            g = list(range(n + 1))
            g[i], g[i + 1] = g[i + 1], g[i]
            gens.append(tuple(g))
        return gens
    pos = lambda i: i - 1          # +i, i = 1..n
    neg = lambda i: n + i - 1      # -i

    def swap(i, j):                # (i j)(-i -j)
        g = list(range(2 * n))
        g[pos(i)], g[pos(j)] = pos(j), pos(i)
        g[neg(i)], g[neg(j)] = neg(j), neg(i)
        return tuple(g)

    def nswap(i, j):               # (i -j)(j -i)
        g = list(range(2 * n))
        g[pos(i)], g[neg(j)] = neg(j), pos(i)
        g[pos(j)], g[neg(i)] = neg(i), pos(j)
        return tuple(g)

    if typ == 'B':
        g0 = list(range(2 * n))
        g0[pos(1)], g0[neg(1)] = neg(1), pos(1)
        return [tuple(g0)] + [swap(i, i + 1) for i in range(1, n)]
    if typ == 'D':
        return [nswap(1, 2)] + [swap(i, i + 1) for i in range(1, n)]
    raise ValueError(typ)


def gens_perm(typ, n, model='roots'):
    """Faithful permutation generators for the Coxeter group of the given type.
    I_2(m) is typ='I', n=m (dihedral of order 2m, acting on Z/m, m >= 3).
    model='roots' (action on the root system, any type) or 'classical'
    ((signed) permutations, types A, B, D only)."""
    if typ == 'I':
        m = n
        assert m >= 3
        s = tuple((-i) % m for i in range(m))
        t = tuple((1 - i) % m for i in range(m))
        return [s, t]
    if model == 'classical':
        return classical_perms(typ, n)
    return root_perms(typ, n)


class Cox:
    """Whole group, enumerated breadth-first by length.

    idx 0 is e.  R[s][i] = index of w_i * s,  Lm[s][i] = index of s * w_i.
    """

    def __init__(self, typ, n, need_left=True, model='roots'):
        self.typ, self.n = typ, n
        gens = gens_perm(typ, n, model)
        self.rank = len(gens)
        self.gens = gens
        deg = len(gens[0])
        self.deg = deg
        pack = bytes if deg < 256 else tuple
        e = pack(range(deg))
        index = {e: 0}
        elems = [e]
        length = [0]
        r = self.rank
        R = [[] for _ in range(r)]
        i = 0
        while i < len(elems):
            w = elems[i]
            for s in range(r):
                g = gens[s]
                ws = pack([w[g[j]] for j in range(deg)])   # (w s)(x) = w(s(x))
                k = index.get(ws)
                if k is None:
                    k = len(elems)
                    index[ws] = k
                    elems.append(ws)
                    length.append(length[i] + 1)
                R[s].append(k)
            i += 1
        self.elems, self.index, self.length, self.R = elems, index, length, R
        self.order = len(elems)
        self.N = max(length)
        inv = [0] * self.order
        for i, w in enumerate(elems):
            wi = [0] * deg
            for j in range(deg):
                wi[w[j]] = j
            inv[i] = index[pack(wi)]
        self.inv = inv
        # left multiplication: s*w = (w^-1 * s)^-1
        if need_left:
            self.Lm = [[inv[R[s][inv[i]]] for i in range(self.order)] for s in range(r)]
        self.pack = pack

    # -- basic helpers -------------------------------------------------------
    def mul(self, i, j):
        a, b = self.elems[i], self.elems[j]
        return self.index[self.pack([a[b[x]] for x in range(self.deg)])]

    def word(self, i):
        """One reduced word (list of generator indices) for element i."""
        out = []
        L, R = self.length, self.R
        while i:
            for s in range(self.rank):
                k = R[s][i]
                if L[k] < L[i]:
                    out.append(s)
                    i = k
                    break
        return out[::-1]

    def from_word(self, word):
        i = 0
        for s in word:
            i = self.R[s][i]
        return i

    def DR(self, i):
        return [s for s in range(self.rank) if self.length[self.R[s][i]] < self.length[i]]

    def DL(self, i):
        return [s for s in range(self.rank) if self.length[self.Lm[s][i]] < self.length[i]]

    def support_mask(self):
        """supp[i] = bitmask of generators in reduced words of w_i."""
        supp = [0] * self.order
        L, R = self.length, self.R
        for i in range(1, self.order):
            for s in range(self.rank):
                k = R[s][i]
                if L[k] < L[i]:
                    supp[i] = supp[k] | (1 << s)
                    break
        return supp

    # -- ancestors -----------------------------------------------------------
    def ancestors(self, check_unique=True):
        """DP over right weak order.
        anc[i]  = index of a(w_i) (0 for e),
        rem[i]  = index of a(w_i)^-1 w_i,
        linv[i] = involution length.
        Involution prefixes of w = {w if w^2=e} + union over right descents s
        of the involution prefixes of ws.  If check_unique, verify that the
        maximal-length involution prefix is unique (Conj. 1.2 / Gaetz-Gao) and
        raise if not."""
        n, L, R, inv = self.order, self.length, self.R, self.inv
        anc = [0] * n
        rem = [0] * n
        linv = [0] * n
        for i in range(1, n):                      # BFS order = nondecreasing length
            if inv[i] == i:
                anc[i] = i
                rem[i] = 0
                linv[i] = 1
                continue
            best, bl, bs = -1, -1, -1
            ties = set()
            for s in range(self.rank):
                k = R[s][i]
                if L[k] < L[i]:
                    a = anc[k]
                    la = L[a]
                    if la > bl:
                        best, bl, bs = a, la, s
                        ties = {a}
                    elif la == bl:
                        ties.add(a)
            if check_unique and len(ties) != 1:
                raise AssertionError(('ancestor not unique', self.typ, self.n, self.word(i)))
            anc[i] = best
            rem[i] = R[bs][rem[R[bs][i]]]          # a^-1 w = (a^-1 (w s)) s
            linv[i] = 1 + linv[rem[i]]
        self.anc, self.rem, self.linv = anc, rem, linv
        return anc, rem, linv

    def decomposition(self, i):
        out = []
        while i:
            out.append(self.anc[i])
            i = self.rem[i]
        return out


def ancestors_bruteforce(G):
    """Independent, definition-level implementation (O(|W| * #involutions)):
    for each w, list all involutions u with l(u) + l(u w) = l(w); take the
    longest; check uniqueness; iterate greedily.  Returns (anc, rem, linv)."""
    n, L, inv = G.order, G.length, G.inv
    invols = [u for u in range(1, n) if inv[u] == u]
    anc = [0] * n
    rem = [0] * n
    for w in range(1, n):
        best, bl, cnt, br = -1, -1, 0, -1
        for u in invols:
            if L[u] > L[w]:
                break                              # BFS order: lengths nondecreasing
            r = G.mul(u, w)                        # u^-1 w = u w
            if L[u] + L[r] == L[w]:
                if L[u] > bl:
                    best, bl, cnt, br = u, L[u], 1, r
                elif L[u] == bl:
                    cnt += 1
        assert cnt == 1, ('ancestor not unique (brute force)', G.word(w))
        anc[w], rem[w] = best, br
    linv = [0] * n
    for w in range(1, n):
        k, x = 0, w
        while x:
            x = rem[x]
            k += 1
        linv[w] = k
    return anc, rem, linv
