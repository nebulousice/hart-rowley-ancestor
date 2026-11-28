"""Independent cross-check for target hart-rowley-ancestor (claims HR-W, HR-T, HR-D).

Different method from hr_tables.py and from scripts/ (no roots, no inversion sets,
no shared code):
  * Group = image of the Tits (geometric) representation, matrices
    sigma_i = I - e_i (2B_i.)  with 2B_ij = 2cos-entries in Z[theta]:
      m=3 -> -1, m=4 -> -sqrt2 (theta^2 = 2), m=5 -> -phi (theta^2 = theta + 1).
    All entries are exact algebraic integers (pairs of Python ints).
  * Length = word metric: breadth-first distance in the Cayley graph from e
    (the definition of l).  Right and left multiplication tables by generators.
  * Prefixes of w (definition: w = u v with l(w) = l(u) + l(v)) are exactly the
    elements reached from w by successive length-decreasing RIGHT multiplications;
    suffixes likewise with LEFT multiplications.  The longest involution prefix is
    found by walking down from w level by level and stopping at the first level
    that contains an involution (M*M = I, M != I, checked on the matrix).
  * The remainder a^{-1} w is read off the letters stripped on the way down.
Decision path: integers only.  No floats.
Usage:  python hr_tables_x.py         (A6 example, H3, F4)
        python hr_tables_x.py --h4    (also H4)
        python hr_tables_x.py --big   (also H4, E6)
"""
import sys
import platform
import itertools
from collections import Counter


class Ring:
    """Z[theta] with theta^2 = p*theta + q; elements are pairs (a, b) = a + b*theta."""
    def __init__(self, p, q):
        self.p, self.q = p, q

    def mul(self, x, y):
        a, b = x
        c, d = y
        bd = b * d
        return (a * c + self.q * bd, a * d + b * c + self.p * bd)


def coxeter_matrix(typ, n):
    m = [[2] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = 1
    def link(i, j, v=3):
        m[i - 1][j - 1] = m[j - 1][i - 1] = v
    if typ == 'A':
        for i in range(1, n):
            link(i, i + 1)
    elif typ == 'H':
        link(1, 2, 5)
        for i in range(2, n):
            link(i, i + 1)
    elif typ == 'F':
        link(1, 2); link(2, 3, 4); link(3, 4)
    elif typ == 'E':
        link(1, 3); link(2, 4)
        for i in range(3, n):
            link(i, i + 1)
    return m


class TitsGroup:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        m = coxeter_matrix(typ, n)
        self.m = m
        ms = {m[i][j] for i in range(n) for j in range(n)}
        if 5 in ms:
            R = Ring(1, 1)
        elif 4 in ms:
            R = Ring(0, 2)
        else:
            R = Ring(0, 0)
        self.R = R
        def twoB(i, j):
            v = m[i][j]
            return {1: (2, 0), 2: (0, 0), 3: (-1, 0), 4: (0, -1), 5: (0, -1)}[v]
        # (sigma_j)_{j,k} = delta_jk - 2B_jk ; other rows are identity rows
        self.c = [[(1 - twoB(j, k)[0], -twoB(j, k)[1]) if j == k else
                   (-twoB(j, k)[0], -twoB(j, k)[1]) for k in range(n)] for j in range(n)]
        ident = tuple((1, 0) if r == s else (0, 0) for r in range(n) for s in range(n))
        self.I = ident
        mats = [ident]
        idx = {ident: 0}
        dist = [0]
        parent = [None]
        right = []
        frontier = [0]
        while frontier:
            nxt = []
            for x in frontier:
                row = []
                for j in range(n):
                    y = self.right_gen(mats[x], j)
                    k = idx.get(y)
                    if k is None:
                        k = len(mats)
                        idx[y] = k
                        mats.append(y)
                        dist.append(dist[x] + 1)
                        parent.append((x, j))
                        nxt.append(k)
                    row.append(k)
                right.append(row)
            frontier = nxt
        # rows were appended in BFS order = index order
        self.mats, self.idx, self.len, self.parent, self.right = mats, idx, dist, parent, right
        self.left = [[idx[self.left_gen(M, j)] for j in range(n)] for M in mats]
        self.is_inv = [x != 0 and self.square(mats[x]) == ident for x in range(len(mats))]
        self.N = len(mats)

    def right_gen(self, M, j):
        """M * sigma_j : column j -> -column j ; column k += c[j][k] * column j (k != j)."""
        n, R = self.n, self.R
        out = list(M)
        for r in range(n):
            mrj = M[r * n + j]
            if mrj == (0, 0):
                continue
            for k in range(n):
                cjk = self.c[j][k]
                if k == j:
                    out[r * n + k] = (-mrj[0], -mrj[1])
                elif cjk != (0, 0):
                    t = R.mul(mrj, cjk)
                    o = out[r * n + k]
                    out[r * n + k] = (o[0] + t[0], o[1] + t[1])
        return tuple(out)

    def left_gen(self, M, j):
        """sigma_j * M : only row j changes: row_j <- sum_k c[j][k] row_k."""
        n, R = self.n, self.R
        out = list(M)
        for s in range(n):
            acc0 = acc1 = 0
            for k in range(n):
                cjk = self.c[j][k]
                if cjk != (0, 0):
                    t = R.mul(cjk, M[k * n + s])
                    acc0 += t[0]
                    acc1 += t[1]
            out[j * n + s] = (acc0, acc1)
        return tuple(out)

    def square(self, M):
        n, R = self.n, self.R
        out = []
        for r in range(n):
            for s in range(n):
                a0 = a1 = 0
                for k in range(n):
                    t = R.mul(M[r * n + k], M[k * n + s])
                    a0 += t[0]
                    a1 += t[1]
                out.append((a0, a1))
        return tuple(out)

    # ---------- words and products via the Cayley graph ----------
    def word(self, x):
        w = []
        while x:
            x, j = self.parent[x]
            w.append(j)
        return w[::-1]

    def trace(self, x, word):
        for j in word:
            x = self.right[x][j]
        return x

    def elem(self, word1):
        return self.trace(0, [j - 1 for j in word1])

    def prod(self, x, y):
        return self.trace(x, self.word(y))

    def inverse(self, x):
        return self.trace(0, self.word(x)[::-1])

    def lexword(self, x):
        """Lexicographically first reduced word: repeatedly strip the smallest left descent."""
        out = []
        while x:
            for j in range(self.n):
                y = self.left[x][j]
                if self.len[y] < self.len[x]:
                    out.append(j + 1)
                    x = y
                    break
        return out

    def wstr(self, x):
        return ''.join('r%d' % j for j in self.lexword(x)) or 'e'

    # ---------- ancestors by walking down ----------
    def top_involutions(self, w, side):
        """side 'R': prefixes (strip on the right); side 'L': suffixes (strip on the left).
        Returns list of (involution, stripped letters in stripping order) at the highest level
        containing an involution."""
        tab = self.right if side == 'R' else self.left
        level = {w: []}
        while level:
            hits = [(x, p) for x, p in level.items() if self.is_inv[x]]
            if hits:
                return hits
            nxt = {}
            for x, p in level.items():
                lx = self.len[x]
                for j in range(self.n):
                    y = tab[x][j]
                    if self.len[y] < lx and y not in nxt:
                        nxt[y] = p + [j]
            level = nxt
        return []

    def anc_step(self, w, side):
        hits = self.top_involutions(w, side)
        a, stripped = hits[0]
        # side R: w = a * s_{jm} ... s_{j1}  (j1 stripped first) -> remainder word reversed
        # side L: w = s_{j1} ... s_{jm} * a  -> remainder (left part) = s_{j1} ... s_{jm}
        if side == 'R':
            rem = self.trace(0, stripped[::-1])
            assert self.trace(a, stripped[::-1]) == w
        else:
            rem = self.trace(0, stripped)
            assert self.trace(rem, self.word(a)) == w
        assert self.len[a] + self.len[rem] == self.len[w]
        return a, rem, len(hits)

    def tables(self, readings=True):
        order = sorted(range(self.N), key=lambda x: self.len[x])
        A = [None] * self.N
        REM = [None] * self.N
        nonuniq = 0
        for w in range(1, self.N):
            a, r, k = self.anc_step(w, 'R')
            A[w], REM[w] = a, r
            nonuniq += (k != 1)
        L = [0] * self.N
        for w in order:
            if w:
                L[w] = 1 + L[REM[w]]
        self.A, self.REM, self.L, self.nonuniq = A, REM, L, nonuniq
        out = {'R1': L}
        if readings:
            AS, REMS = [None] * self.N, [None] * self.N
            for w in range(1, self.N):
                b, r, k = self.anc_step(w, 'L')
                AS[w], REMS[w] = b, r
            LS = [0] * self.N
            for w in order:
                if w:
                    LS[w] = 1 + LS[REMS[w]]
            out['R2'] = LS
            self.AS, self.REMS, self.LS = AS, REMS, LS
            def alt(w, side):
                k = 0
                while w:
                    w = REM[w] if side == 0 else REMS[w]
                    side ^= 1
                    k += 1
                return k
            out['R3'] = [alt(w, 0) for w in range(self.N)]
            out['R3b'] = [alt(w, 1) for w in range(self.N)]
            # R4: minimum over all length-additive involution factorisations
            MN = [0] * self.N
            for w in order:
                if w == 0:
                    continue
                best = None
                seen = {w: []}
                stack = [w]
                while stack:
                    x = stack.pop()
                    if self.is_inv[x]:
                        r = self.trace(0, seen[x][::-1])
                        v = 1 + MN[r]
                        if best is None or v < best:
                            best = v
                    for j in range(self.n):
                        y = self.right[x][j]
                        if self.len[y] < self.len[x] and y not in seen:
                            seen[y] = seen[x] + [j]
                            stack.append(y)
                MN[w] = best
            out['R4'] = MN
        return out

    def decomposition(self, w):
        out = []
        while w:
            out.append(self.A[w])
            w = self.REM[w]
        return out

    def suffix_decomposition(self, w):
        out = []
        while w:
            out.append(self.AS[w])
            w = self.REMS[w]
        return out[::-1]


def report(G, out):
    n = G.n
    print('== %s%d (Tits matrices over %s) ==' % (G.typ, n,
          {(1, 1): 'Z[phi]', (0, 2): 'Z[sqrt2]', (0, 0): 'Z'}[(G.R.p, G.R.q)]))
    print('  |W| = %d, #involutions = %d, max length = %d'
          % (G.N, sum(G.is_inv), max(G.len)))
    inv_word = sum(1 for x in range(1, G.N) if (G.inverse(x) == x) != G.is_inv[x])
    print('  involution test (matrix square) vs word reversal: mismatches = %d' % inv_word)
    print('  ancestor uniqueness failures: %d' % G.nonuniq)
    L = out['R1']
    M = max(L)
    d = Counter(L)
    print('  R1 (source) distribution l_inv = 0..%d: %s' % (M, [d[k] for k in range(M + 1)]))
    print('  R1 (source) max l_inv = %d ; elements with l_inv > rank (%d): %d'
          % (M, n, sum(1 for x in L if x > n)))
    res = {}
    for key in ('R1', 'R2', 'R3', 'R3b', 'R4'):
        if key in out:
            T = out[key]
            res[key] = (max(T), sum(1 for x in T if x > n))
            if key != 'R1':
                print('  %s: max = %d ; elements > rank: %d' % (key, res[key][0], res[key][1]))
    # Coxeter elements: Prop 2.2 (a(w) = product of left descents)
    seen, bad = set(), 0
    for perm in itertools.permutations(range(1, n + 1)):
        w = G.elem(list(perm))
        if w in seen:
            continue
        seen.add(w)
        DL = [j + 1 for j in range(n) if G.len[G.left[w][j]] < G.len[w]]
        if G.A[w] != G.elem(DL):
            bad += 1
    print('  Coxeter elements: %d distinct; Prop 2.2 failures %d' % (len(seen), bad))
    res['cox'] = (len(seen), bad)
    return res


def main():
    big = '--big' in sys.argv
    h4only = '--h4' in sys.argv
    print('hr_tables_x.py (cross-check): python %s' % platform.python_version())
    ok = True
    G = TitsGroup('A', 6)
    G.tables(readings=True)
    w = G.elem([6, 3, 2, 1, 4, 5])
    pd = [G.wstr(x) for x in G.decomposition(w)]
    sd = [G.wstr(x) for x in G.suffix_decomposition(w)]
    exp_p = [G.wstr(G.elem(x)) for x in ([3, 6], [2, 4], [1, 5])]
    exp_s = [G.wstr(G.elem(x)) for x in ([3], [2, 4, 6], [1, 5])]
    print('A6 example r6r3r2r1r4r5: prefix-ancestor decomposition %s ; suffix decomposition %s' % (pd, sd))
    okA = pd == exp_p and sd == exp_s
    print('  matches source text: %s' % okA)
    ok &= okA

    G = TitsGroup('H', 3)
    rH3 = report(G, G.tables())
    w = G.elem([2, 1, 2, 3, 1, 2, 1])
    dec = G.decomposition(w)
    print('H3 witness r2r1r2r3r1r2r1: l = %d ; decomposition %s ; lengths %s ; l_inv = %d'
          % (G.len[w], [G.wstr(x) for x in dec], [G.len[x] for x in dec], G.L[w]))
    # explicit prefix list of w by walking down
    pre, stack = {w}, [w]
    while stack:
        x = stack.pop()
        for j in range(3):
            y = G.right[x][j]
            if G.len[y] < G.len[x] and y not in pre:
                pre.add(y)
                stack.append(y)
    print('  prefixes (%d): %s' % (len(pre), sorted((G.len[x], G.wstr(x)) for x in pre)))
    print('  involution prefixes: %s' % sorted(G.wstr(x) for x in pre if G.is_inv[x]))
    okW = (dec == [G.elem(x) for x in ([2, 1, 2], [3, 1], [2], [1])]) and G.L[w] == 4 and G.len[w] == 7
    print('  witness decomposition equals (r2r1r2)(r3r1)(r2)(r1): %s' % okW)
    ok &= okW
    viol = sorted((x for x in range(G.N) if G.L[x] > 3), key=lambda x: (G.len[x], G.lexword(x)))
    print('  all H3 violators: %s' % [(G.wstr(x), G.len[x], [G.wstr(a) for a in G.decomposition(x)]) for x in viol])

    G = TitsGroup('F', 4)
    rF4 = report(G, G.tables())
    w = G.elem([3, 4, 1, 2, 3, 2, 1])
    print('F4 witness r3r4r1r2r3r2r1: l = %d ; decomposition %s ; l_inv = %d'
          % (G.len[w], [G.wstr(x) for x in G.decomposition(w)], G.L[w]))
    ok &= G.L[w] == 5
    ok &= rH3['R1'] == (4, 3) and rF4['R1'] == (5, 50)
    ok &= rH3['cox'][1] == 0 and rF4['cox'][1] == 0
    if h4only:
        G = TitsGroup('H', 4)
        rH4 = report(G, G.tables(readings=False))
        ok &= rH4['R1'] == (7, 1746)
        print('CERTIFIED: H4 max l_inv = %d (%d elements > 4)' % rH4['R1'])
    if big:
        G = TitsGroup('H', 4)
        rH4 = report(G, G.tables(readings=False))
        G = TitsGroup('E', 6)
        rE6 = report(G, G.tables(readings=False))
        ok &= rH4['R1'] == (7, 1746) and rE6['R1'] == (8, 160)
    print('CERTIFIED: H3 max l_inv = %d (%d elements > 3); F4 max l_inv = %d (%d elements > 4)'
          % (rH3['R1'] + rF4['R1']))
    if big:
        print('CERTIFIED: H4 max l_inv = %d (%d elements > 4); E6 max l_inv = %d (%d elements > 6)'
              % (rH4['R1'] + rE6['R1']))
    print('PASS' if ok else 'FAIL')


if __name__ == '__main__':
    main()
