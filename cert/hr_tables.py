"""Primary certificate for target hart-rowley-ancestor (claims HR-W, HR-T, HR-D).

Method (written from scratch for the Certifier, no code shared with scripts/):
  * Root system built exactly over Z[phi] (phi^2 = phi + 1); crystallographic
    types only use the Z part.  Root coordinates in the basis of simple roots.
  * Group elements = permutations of the root set (faithful action).
  * Length l(w) = |Phi(w)|, Phi(w) = {beta > 0 : w^{-1} beta < 0} (left inversion set).
  * Prefix order: u <=_R w  (w = u v with l(w) = l(u) + l(v))  iff  Phi(u) subset Phi(w)
    (Bjorner-Brenti Prop. 3.1.3); additionally SELF-TESTED below against the
    length-additivity definition for every pair (u, w) in H3 and F4.
  * Ancestor a(w): brute force over ALL involutions u with Phi(u) subset Phi(w),
    maximal |Phi(u)|; uniqueness is checked for every w != e.
  * Involution length (source, arXiv:2502.00777v1): w = a(w) * rem, rem = a(w)^{-1} w,
    l_inv(w) = 1 + l_inv(rem), l_inv(e) = 0.

Alternative readings (definitions audit):
  R1  source reading: prefixes, ancestor peeled on the left.
  R2  suffix reading: longest involution suffix peeled on the right (= R1 on w^{-1}).
  R3  alternating: prefix-ancestor, then suffix-ancestor of the remainder, and so on.
  R3' alternating starting with the suffix side.
  R4  non-greedy minimum: least k with w = u_1...u_k length-additive, u_i involutions.

Decision path: integers only (Python int, exact sign test in Z[phi]).  No floats.
Usage:  python hr_tables.py            (A6 example, H3, F4)
        python hr_tables.py --big      (also H4, E6)
"""
import sys
import itertools
import platform
from collections import deque, Counter

# ---------------- Z[phi] arithmetic, x = (a, b) means a + b*phi ----------------
def zadd(x, y):
    return (x[0] + y[0], x[1] + y[1])

def zmul(x, y):
    a, b = x
    c, d = y
    return (a * c + b * d, a * d + b * c + b * d)

def zneg(x):
    return (-x[0], -x[1])

def zsign(x):
    # a + b*phi = (2a + b + b*sqrt5) / 2 ; sign of p + q*sqrt5 decided exactly
    a, b = x
    p, q = 2 * a + b, b
    if p == 0 and q == 0:
        return 0
    if p >= 0 and q >= 0:
        return 1
    if p <= 0 and q <= 0:
        return -1
    # opposite signs: compare p^2 with 5 q^2
    if p > 0:
        return 1 if p * p > 5 * q * q else -1
    return 1 if 5 * q * q > p * p else -1

ZERO, ONE, PHI = (0, 0), (1, 0), (0, 1)

# ---------------- Coxeter data ----------------
def coxeter_edges(typ, n):
    """Return dict {(i,j): m_ij} for m_ij >= 3, 0-based, in the brief's labelling."""
    e = {}
    def put(i, j, m):
        e[(i, j)] = m
        e[(j, i)] = m
    if typ == 'A':
        for i in range(n - 1):
            put(i, i + 1, 3)
    elif typ == 'H':          # r1 -5- r2 - r3 (- r4)
        put(0, 1, 5)
        for i in range(1, n - 1):
            put(i, i + 1, 3)
    elif typ == 'F':          # r1 - r2 -4- r3 - r4
        put(0, 1, 3); put(1, 2, 4); put(2, 3, 3)
    elif typ == 'B':          # r1 -4- r2 - ... - rn
        put(0, 1, 4)
        for i in range(1, n - 1):
            put(i, i + 1, 3)
    elif typ == 'E':          # Bourbaki: r1-r3-r4-...-rn, r2-r4
        put(0, 2, 3); put(1, 3, 3)
        for i in range(2, n - 1):
            put(i, i + 1, 3)
    else:
        raise ValueError(typ)
    return e

def cartan_entry(i, j, edges):
    """c(i,j) with s_i(alpha_j) = alpha_j - c(i,j) alpha_i ; c(i,j) c(j,i) = 4cos^2(pi/m)."""
    if i == j:
        return (2, 0)
    m = edges.get((i, j), 2)
    if m == 2:
        return ZERO
    if m == 3:
        return (-1, 0)
    if m == 4:
        return (-1, 0) if i < j else (-2, 0)
    if m == 5:
        return (0, -1)       # -phi, symmetric: phi^2 = 4cos^2(pi/5)
    raise ValueError(m)

class Group:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        self.edges = coxeter_edges(typ, n)
        C = [[cartan_entry(i, j, self.edges) for j in range(n)] for i in range(n)]
        self.C = C
        # roots
        def refl(i, beta):
            # s_i(beta) = beta - <beta, i> alpha_i, linear: <alpha_j, i> = c(i,j)
            coef = ZERO
            for j in range(n):
                coef = zadd(coef, zmul(beta[j], C[i][j]))
            nb = list(beta)
            nb[i] = zadd(nb[i], zneg(coef))
            return tuple(nb)
        simple = [tuple(ONE if k == j else ZERO for k in range(n)) for j in range(n)]
        roots = list(simple)
        index = {r: k for k, r in enumerate(roots)}
        q = deque(roots)
        while q:
            b = q.popleft()
            for i in range(n):
                c = refl(i, b)
                if c not in index:
                    index[c] = len(roots)
                    roots.append(c)
                    q.append(c)
        self.roots, self.rindex = roots, index
        def positive(b):
            for x in b:
                s = zsign(x)
                if s != 0:
                    return s > 0
            raise AssertionError('zero root')
        self.pos = [positive(b) for b in roots]
        # every root has all coefficients of one sign (sanity)
        for b in roots:
            sg = {zsign(x) for x in b} - {0}
            assert len(sg) == 1
        self.npos = sum(self.pos)
        self.posidx = [k for k in range(len(roots)) if self.pos[k]]
        self.bit = {k: 1 << t for t, k in enumerate(self.posidx)}
        self.gen = [tuple(index[refl(i, b)] for b in roots) for i in range(n)]
        self.simple_idx = list(range(n))
        # group by BFS on permutations
        ident = tuple(range(len(roots)))
        self.ident = ident
        elems = [ident]
        eidx = {ident: 0}
        q = deque([ident])
        while q:
            w = q.popleft()
            for g in self.gen:
                x = self.comp(w, g)
                if x not in eidx:
                    eidx[x] = len(elems)
                    elems.append(x)
                    q.append(x)
        self.elems, self.eidx = elems, eidx
        self.inv = [self.eidx[self.inverse(w)] for w in elems]
        self.mask = [self.phi_mask(w) for w in elems]
        self.length = [bin(m).count('1') for m in self.mask]
        self.invols = [k for k, w in enumerate(elems)
                       if w != ident and self.comp(w, w) == ident]
        self.order = sorted(range(len(elems)), key=lambda k: self.length[k])

    @staticmethod
    def comp(p, q):
        """p o q : apply q first, then p."""
        return tuple(p[k] for k in q)

    @staticmethod
    def inverse(p):
        r = [0] * len(p)
        for k, v in enumerate(p):
            r[v] = k
        return tuple(r)

    def phi_mask(self, w):
        wi = self.inverse(w)
        m = 0
        for k in self.posidx:
            if not self.pos[wi[k]]:
                m |= self.bit[k]
        return m

    def mul(self, a, b):
        return self.eidx[self.comp(self.elems[a], self.elems[b])]

    def from_word(self, word):  # word of 1-based generator labels, w = s_{i1} ... s_{ik}
        w = self.ident
        for i in word:
            w = self.comp(w, self.gen[i - 1])
        return self.eidx[w]

    def word(self, k):
        """Lexicographically first reduced word, by peeling left descents."""
        out = []
        while k != 0:
            for i in range(self.n):
                if self.mask[k] & self.bit[i]:      # alpha_i in Phi(w)  <=>  s_i left descent
                    out.append(i + 1)
                    k = self.mul(self.from_word([i + 1]), k)
                    break
        return out

    def wstr(self, k):
        return ''.join('r%d' % i for i in self.word(k)) or 'e'

    # ---------- ancestors ----------
    def ancestors(self):
        """a[w] for all w != e (brute force over all involutions), with uniqueness check."""
        a = [None] * len(self.elems)
        fails = 0
        for w in range(1, len(self.elems)):
            mw = self.mask[w]
            best, cand = -1, []
            for u in self.invols:
                mu = self.mask[u]
                if mu & ~mw == 0:
                    lu = self.length[u]
                    if lu > best:
                        best, cand = lu, [u]
                    elif lu == best:
                        cand.append(u)
            if len(cand) != 1:
                fails += 1
            a[w] = cand[0]
        self.anc = a
        return fails

    def rem(self, w):
        return self.mul(self.anc[w], w)       # a^{-1} w = a w

    def ilen_table(self):
        L = [0] * len(self.elems)
        for w in self.order:
            if w:
                L[w] = 1 + L[self.rem(w)]
        self.L = L
        return L

    def decomposition(self, w):
        out = []
        while w:
            out.append(self.anc[w])
            w = self.rem(w)
        return out

    def suffix_decomposition(self, w):
        """Longest involution suffix peeled on the right; returned left-to-right."""
        out = []
        while w:
            b = self.anc[self.inv[w]]          # longest involution suffix of w = a(w^{-1})
            out.append(b)
            w = self.mul(w, b)                 # w b^{-1} = w b
        return out[::-1]

    def alt_len(self, w, side):
        k = 0
        while w:
            if side == 0:
                w = self.rem(w)
            else:
                w = self.mul(w, self.anc[self.inv[w]])
            side ^= 1
            k += 1
        return k

    def min_table(self):
        M = [0] * len(self.elems)
        for w in self.order:
            if w == 0:
                continue
            mw = self.mask[w]
            best = None
            for u in self.invols:
                if self.mask[u] & ~mw == 0:
                    v = M[self.mul(u, w)]
                    if best is None or v + 1 < best:
                        best = v + 1
            M[w] = best
        return M

    def selftest_prefix_criterion(self):
        """Phi(u) subset Phi(w)  <=>  l(u) + l(u^{-1} w) = l(w), all pairs."""
        bad = 0
        N = len(self.elems)
        for u in range(N):
            ui = self.inv[u]
            for w in range(N):
                lhs = (self.mask[u] & ~self.mask[w]) == 0
                rhs = self.length[u] + self.length[self.mul(ui, w)] == self.length[w]
                if lhs != rhs:
                    bad += 1
        return bad

    def coxeter_element_check(self):
        """Source Prop. 2.2 (a(w) = prod D_L(w)) and Prop. 2.5 (l_inv = path length)."""
        n = self.n
        bad22 = bad25 = 0
        seen = set()
        for perm in itertools.permutations(range(1, n + 1)):
            w = self.from_word(list(perm))
            if w in seen:
                continue
            seen.add(w)
            DL = [i + 1 for i in range(n) if self.mask[w] & self.bit[i]]
            if self.anc[w] != self.from_word(DL):
                bad22 += 1
            # longest w-path: path in the Coxeter graph with increasing positions
            pos = {g: t for t, g in enumerate(perm)}
            best = {}
            for g in sorted(range(1, n + 1), key=lambda g: pos[g]):
                best[g] = 1 + max([best[h] for h in range(1, n + 1)
                                   if (g - 1, h - 1) in self.edges and pos[h] < pos[g]] or [0])
            # (a path in the DAG of edges oriented by position = a w-path)
            if self.L[w] != max(best.values()):
                bad25 += 1
        return len(seen), bad22, bad25


def run_type(typ, n, full_readings=True, selftest=False):
    G = Group(typ, n)
    fails = G.ancestors()
    L = G.ilen_table()
    N = len(G.elems)
    res = {}
    print('== %s%d ==' % (typ, n))
    print('  |W| = %d, #roots = %d, #positive roots = %d, #involutions = %d, max length = %d'
          % (N, len(G.roots), G.npos, len(G.invols), max(G.length)))
    if selftest:
        print('  self-test prefix criterion (all %d pairs): mismatches = %d' % (N * N, G.selftest_prefix_criterion()))
    print('  ancestor uniqueness failures (w != e): %d' % fails)
    dist = Counter(L)
    M = max(L)
    print('  R1 (source) distribution l_inv = 0..%d: %s' % (M, [dist[k] for k in range(M + 1)]))
    print('  R1 (source) max l_inv = %d ; elements with l_inv > rank (%d): %d' % (M, n, sum(1 for x in L if x > n)))
    res['R1'] = (M, sum(1 for x in L if x > n))
    viol = [w for w in G.order if L[w] > n]
    if viol:
        w0 = viol[0]
        print('  shortest violators (length %d): %s' % (G.length[w0],
              ', '.join(G.wstr(w) for w in viol if G.length[w] == G.length[w0])))
    if full_readings:
        L2 = [L[G.inv[w]] for w in range(N)]
        # R2 computed directly from suffix decomposition as well, to check L2
        L2d = [len(G.suffix_decomposition(w)) for w in range(N)]
        assert L2 == L2d
        L3 = [G.alt_len(w, 0) for w in range(N)]
        L3b = [G.alt_len(w, 1) for w in range(N)]
        L4 = G.min_table()
        for name, T in (('R2 suffix', L2), ('R3 alt(prefix first)', L3),
                        ('R3b alt(suffix first)', L3b), ('R4 min length-additive', L4)):
            print('  %s: max = %d ; elements > rank: %d' % (name, max(T), sum(1 for x in T if x > n)))
        res['R2'] = (max(L2), sum(1 for x in L2 if x > n))
        res['R3'] = (max(L3), sum(1 for x in L3 if x > n))
        res['R3b'] = (max(L3b), sum(1 for x in L3b if x > n))
        res['R4'] = (max(L4), sum(1 for x in L4 if x > n))
    ce, b22, b25 = G.coxeter_element_check()
    print('  Coxeter elements: %d distinct; Prop 2.2 failures %d; Prop 2.5 (path length) failures %d'
          % (ce, b22, b25))
    res['cox'] = (ce, b22, b25)
    res['fails'] = fails
    return G, res


def main():
    big = '--big' in sys.argv
    print('hr_tables.py (primary): python %s' % platform.python_version())
    ok = True

    # ---- the source's worked example, type A6 (r_i = (i, i+1)) ----
    G = Group('A', 6)
    G.ancestors(); G.ilen_table()
    w = G.from_word([6, 3, 2, 1, 4, 5])
    pd = [G.wstr(x) for x in G.decomposition(w)]
    sd = [G.wstr(x) for x in G.suffix_decomposition(w)]
    print('A6 example w = r6r3r2r1r4r5: ancestor decomposition', pd, ' suffix decomposition', sd)
    exp_p = [G.wstr(G.from_word(x)) for x in ([3, 6], [2, 4], [1, 5])]
    exp_s = [G.wstr(G.from_word(x)) for x in ([3], [2, 4, 6], [1, 5])]
    okA = (pd == exp_p and sd == exp_s)
    # does the swapped reading reproduce the source's example?  (it must not)
    print('  matches source (prefix -> (r3r6)(r2r4)(r1r5), suffix -> (r3)(r2r4r6)(r1r5)): %s' % okA)
    print('  swapped reading would give prefix=%s suffix=%s' % (sd, pd))
    ok &= okA

    # ---- H3 witness ----
    G, rH3 = run_type('H', 3, selftest=True)
    w = G.from_word([2, 1, 2, 3, 1, 2, 1])
    dec = G.decomposition(w)
    decw = [G.word(x) for x in dec]
    print('H3 witness w = r2r1r2r3r1r2r1: l(w) = %d, decomposition %s, lengths %s, l_inv = %d'
          % (G.length[w], ['r' + 'r'.join(map(str, x)) for x in decw], [G.length[x] for x in dec], G.L[w]))
    # prefixes and involution prefixes of w (by inversion set)
    pre = [u for u in range(len(G.elems)) if G.mask[u] & ~G.mask[w] == 0]
    ipre = [u for u in pre if u in set(G.invols)]
    print('  #prefixes = %d ; involution prefixes: %s' % (len(pre), sorted(G.wstr(u) for u in ipre)))
    expected = [G.from_word(x) for x in ([2, 1, 2], [3, 1], [2], [1])]
    okW = (dec == expected and G.L[w] == 4 and G.length[w] == 7)
    print('  witness decomposition equals (r2r1r2)(r3r1)(r2)(r1): %s' % okW)
    ok &= okW
    # reading R4 on the witness: all length-additive factorisations into min(#) involutions
    MT = G.min_table()
    def all_min(x):
        if x == 0:
            return [[]]
        out = []
        for u in G.invols:
            if G.mask[u] & ~G.mask[x] == 0:
                y = G.mul(u, x)
                if MT[y] == MT[x] - 1:
                    out += [[u] + rest for rest in all_min(y)]
        return out
    fac = all_min(w)
    print('  R4 on the witness: min #involutions = %d ; all minimal factorisations: %s'
          % (MT[w], ['(' + ')('.join(G.wstr(u) for u in f) + ')' for f in fac]))
    G, rF4 = run_type('F', 4, selftest=True)
    w = G.from_word([3, 4, 1, 2, 3, 2, 1])
    print('F4 witness w = r3r4r1r2r3r2r1: l(w) = %d, decomposition %s, l_inv = %d'
          % (G.length[w], [G.wstr(x) for x in G.decomposition(w)], G.L[w]))
    okF = G.L[w] == 5
    ok &= okF
    ok &= rH3['R1'] == (4, 3) and rF4['R1'] == (5, 50)
    ok &= rH3['fails'] == 0 and rF4['fails'] == 0
    ok &= rH3['cox'][1:] == (0, 0) and rF4['cox'][1:] == (0, 0)
    if big:
        G, rH4 = run_type('H', 4)
        G, rE6 = run_type('E', 6)
        ok &= rH4['R1'] == (7, 1746) and rE6['R1'] == (8, 160)
        ok &= rH4['fails'] == 0 and rE6['fails'] == 0
    print('CERTIFIED: H3 max l_inv = %d (%d elements > 3); F4 max l_inv = %d (%d elements > 4)'
          % (rH3['R1'] + rF4['R1']))
    if big:
        print('CERTIFIED: H4 max l_inv = %d (%d elements > 4); E6 max l_inv = %d (%d elements > 6)'
              % (rH4['R1'] + rE6['R1']))
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
