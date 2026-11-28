"""Second, independent implementation of the ancestor decomposition.

Independent of coxlib.py in every ingredient:
  * the group comes from the Coxeter PRESENTATION <s_i | s_i^2, (s_i s_j)^m_ij>
    by Todd-Coxeter coset enumeration over the trivial subgroup
    (sympy.combinatorics), not from a root system or a permutation model;
  * length = breadth-first distance in the resulting Cayley table;
  * the prefixes of w are found top-down (strip right descents), the
    involution test traces the word of u twice through the table, and the
    ancestor is taken straight from the definition (longest involution
    prefix, uniqueness asserted) -- no dynamic programme.

Usage:  python second_impl.py            (H3, F4, A3, A4, B3, D4 in full)
        python second_impl.py H 4 word   (check one element, word 0-based)
"""
import sys
from collections import Counter, deque
from sympy.combinatorics.free_groups import free_group
from sympy.combinatorics.fp_groups import FpGroup
from sympy.combinatorics.coset_table import coset_enumeration_r


def coxeter_matrix(typ, n):
    """m[i][j] written out by hand, independently of coxlib.cartan."""
    m = [[1 if i == j else 2 for j in range(n)] for i in range(n)]

    def e(i, j, v=3):
        m[i][j] = m[j][i] = v

    if typ == 'A':
        for i in range(n - 1):
            e(i, i + 1)
    elif typ == 'B':
        e(0, 1, 4)
        for i in range(1, n - 1):
            e(i, i + 1)
    elif typ == 'D':
        e(0, 2)
        e(1, 2)
        for i in range(2, n - 1):
            e(i, i + 1)
    elif typ == 'E':
        e(0, 2)
        e(1, 3)
        for i in range(2, n - 1):
            e(i, i + 1)
    elif typ == 'F':
        e(0, 1)
        e(1, 2, 4)
        e(2, 3)
    elif typ == 'H':
        e(0, 1, 5)
        for i in range(1, n - 1):
            e(i, i + 1)
    else:
        raise ValueError(typ)
    return m


class Table:
    def __init__(self, typ, n):
        m = coxeter_matrix(typ, n)
        F, *x = free_group(','.join('s%d' % i for i in range(n)))
        rels = [x[i] ** 2 for i in range(n)]
        rels += [(x[i] * x[j]) ** m[i][j] for i in range(n) for j in range(i + 1, n)]
        G = FpGroup(F, rels)
        C = coset_enumeration_r(G, [], max_cosets=4000000)
        C.compress()
        C.standardize()
        self.n = n
        self.R = [[row[2 * s] for row in C.table] for s in range(n)]   # w -> w*s
        self.order = len(C.table)
        # lengths and one reduced word per element (BFS from the identity)
        L = [-1] * self.order
        W = [None] * self.order
        L[0], W[0] = 0, ()
        q = deque([0])
        while q:
            i = q.popleft()
            for s in range(n):
                k = self.R[s][i]
                if L[k] < 0:
                    L[k] = L[i] + 1
                    W[k] = W[i] + (s,)
                    q.append(k)
        self.L, self.W = L, W

    def trace(self, start, word):
        for s in word:
            start = self.R[s][start]
        return start

    def is_involution(self, u):
        return u != 0 and self.trace(u, self.W[u]) == 0

    def prefixes(self, w):
        """dict u -> word v with u*v = w, l(u)+l(v)=l(w)."""
        out = {w: ()}
        q = deque([w])
        while q:
            x = q.popleft()
            for s in range(self.n):
                k = self.R[s][x]
                if self.L[k] < self.L[x] and k not in out:
                    out[k] = (s,) + out[x]
                    q.append(k)
        return out

    def ancestor(self, w):
        P = self.prefixes(w)
        inv = [u for u in P if self.is_involution(u)]
        mx = max(self.L[u] for u in inv)
        top = [u for u in inv if self.L[u] == mx]
        assert len(top) == 1, 'ancestor not unique'
        a = top[0]
        return a, self.trace(0, P[a])

    def decomposition(self, w):
        out = []
        while w:
            a, w = self.ancestor(w)
            out.append(list(self.W[a]))
        return out


def full(typ, n):
    T = Table(typ, n)
    memo = {0: 0}
    for w in sorted(range(1, T.order), key=lambda i: T.L[i]):
        a, r = T.ancestor(w)
        memo[w] = 1 + memo[r]
    c = Counter(memo.values())
    return T.order, max(T.L), dict(sorted(c.items()))


if __name__ == '__main__':
    if len(sys.argv) > 3:
        typ, n = sys.argv[1], int(sys.argv[2])
        word = [int(c) for c in sys.argv[3].split(',')]
        T = Table(typ, n)
        w = T.trace(0, word)
        d = T.decomposition(w)
        print(typ, n, 'order', T.order, 'word', word, 'length', T.L[w],
              'l_inv', len(d), 'decomposition', d)
    else:
        for typ, n in [('A', 3), ('A', 4), ('B', 3), ('D', 4), ('H', 3), ('F', 4)]:
            print(typ, n, *full(typ, n), flush=True)
