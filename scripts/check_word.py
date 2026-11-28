"""Check the ancestor decomposition of ONE element given by a word, without
enumerating the group.  Works in any finite type, including E7 and E8.

Method (different from coxlib.Cox and from big.py):
  * elements are permutations of the root system (exact, from coxlib.cartan);
  * length = number of positive roots sent to negative roots;
  * s is a right descent of w  iff  w(alpha_s) is negative;
  * the prefixes of w are found top-down by stripping right descents;
  * u is an involution iff u*u = identity (permutation composition);
  * the ancestor is the longest involution prefix (uniqueness asserted).

Usage:  python check_word.py TYPE RANK w0,w1,...     (0-based generator indices)
        python check_word.py                          (the witnesses quoted in BRIEF.md)
"""
import sys
from collections import deque
from coxlib import cartan, zmul, zsub, ZERO


def setup(typ, n):
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

    idx = {r: k for k, r in enumerate(simple)}
    roots = list(simple)
    q = deque(simple)
    while q:
        b = q.popleft()
        for i in range(n):
            c = refl(i, b)
            if c not in idx:
                idx[c] = len(roots)
                roots.append(c)
                q.append(c)
    gens = [tuple(idx[refl(i, b)] for b in roots) for i in range(n)]

    def positive(beta):
        # a root is positive iff its first nonzero simple-root coordinate is > 0;
        # in Z[phi], a + b*phi > 0 is decided exactly (phi = (1+sqrt5)/2)
        for x in beta:
            if x != ZERO:
                a, b = x
                # sign of a + b*phi = sign of (2a + b) + b*sqrt5
                u, v = 2 * a + b, b
                if u >= 0 and v >= 0:
                    return True
                if u <= 0 and v <= 0:
                    return False
                return (u * u > 5 * v * v) if u > 0 else (u * u < 5 * v * v)
        raise ValueError
    pos = [positive(b) for b in roots]
    return gens, pos


def decomposition(typ, n, word):
    gens, pos = setup(typ, n)
    deg = len(pos)
    ident = tuple(range(deg))

    def mul(a, b):
        return tuple(a[b[x]] for x in range(deg))

    def length(w):
        return sum(1 for x in range(deg) if pos[x] and not pos[w[x]])

    w = ident
    for s in word:
        w = mul(w, gens[s])
    out = []
    total = length(w)
    while w != ident:
        # prefixes: u -> word v with u v = w reduced
        pref = {w: ()}
        q = deque([w])
        while q:
            x = q.popleft()
            for s in range(n):
                if not pos[x[s]]:                 # roots are indexed with simple roots first
                    y = mul(x, gens[s])
                    if y not in pref:
                        pref[y] = (s,) + pref[x]
                        q.append(y)
        invs = [u for u in pref if u != ident and mul(u, u) == ident]
        mx = max(length(u) for u in invs)
        top = [u for u in invs if length(u) == mx]
        assert len(top) == 1, 'ancestor not unique'
        a = top[0]
        out.append((mx, len(pref)))
        r = ident
        for s in pref[a]:
            r = mul(r, gens[s])
        assert length(r) + mx == length(w)
        w = r
    return total, out


WITNESSES = [
    ('H', 3, [1, 0, 1, 2, 0, 1, 0]),
    ('F', 4, [2, 3, 0, 1, 2, 1, 0]),
    ('H', 4, [1, 2, 0, 1, 0, 1, 2, 3]),
    ('H', 4, [3, 1, 2, 0, 1, 0, 1, 2, 3, 2, 1, 0, 1, 2, 0, 1, 0, 1, 2, 3]),
    ('E', 6, [3, 4, 5, 0, 2, 3, 4, 1, 3, 2, 0]),
    ('E', 6, [2, 1, 3, 4, 5, 0, 2, 3, 4, 1, 3, 2, 0]),
    ('E', 7, [6, 4, 5, 3, 4, 2, 3, 0, 2, 1, 3, 4, 5, 6]),
    ('E', 8, [6, 4, 5, 3, 4, 2, 3, 0, 2, 1, 3, 4, 5, 6]),
]

if __name__ == '__main__':
    if len(sys.argv) > 3:
        jobs = [(sys.argv[1], int(sys.argv[2]), [int(c) for c in sys.argv[3].split(',')])]
    else:
        jobs = WITNESSES
    for typ, n, word in jobs:
        total, dec = decomposition(typ, n, word)
        print('%s%d  word %s  length %d  l_inv = %d  (rank %d)  ancestor lengths %s  #prefixes per step %s' % (
            typ, n, ' '.join('r%d' % (s + 1) for s in word), total, len(dec), n,
            [d[0] for d in dec], [d[1] for d in dec]), flush=True)
