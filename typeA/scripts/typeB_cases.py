"""Small-case checks for PROOF_typeA.md section 5 (hart-rowley-ancestor, HR-3).

B_n = signed permutations of {+-1..+-n}, viewed inside S_2n through the
order-preserving relabelling  -n < ... < -1 < 1 < ... < n  ->  1..2n.
Everything from the definitions, exact integers, pure Python, no shared code.

  B-length  : word length in s_0 = (-1,2,..,n), s_i = (i i+1), by breadth-first search
  D-length  : word length in s_0' = (-2,-1,3,..,n), s_i, by breadth-first search
  A-length  : number of inversions of the image in S_2n

Checks for B_n:
  (b1) l_B = inv + neg + nsp, and l_A = 2 l_B - neg                (Lemma 5.1)
  (b2) l_B(xy) = l_B(x) + l_B(y) - (|P| + |F(P)|)                  (Lemma 5.2, n <= 3)
  (b3) u <=_R w in B_n  iff  u <=_R w in S_2n, all pairs           (Lemma 5.2)
  (b4) a_B(w) (longest involution B-prefix, uniqueness tested) = a_A(w)
       (longest involution prefix in S_2n, over ALL involutions of S_2n)   (Prop. 5.3)
  (b5) l_inv^B = l_inv^A, l_inv <= disp in S_2n                    (Thm 5.4)
For D_n the same comparison is made and the mismatches are counted (Remark 5.6),
together with  (d1) l_D = inv + nsp, l_A = 2 l_D + neg  and
(d2) l_D(xy) = l_D(x) + l_D(y) - (|P| - |F(P)|)  (n <= 3).

Usage: python typeB_cases.py NMAX     (NMAX = 4 takes well under a minute)
"""
import sys, json, itertools
from collections import deque


def smul(u, v):
    """(u v)(i) = u(v(i)) for signed windows."""
    return tuple((u[abs(x) - 1] if x > 0 else -u[abs(x) - 1]) for x in v)


def sinv(u):
    r = [0] * len(u)
    for i, x in enumerate(u):
        r[abs(x) - 1] = (i + 1) if x > 0 else -(i + 1)
    return tuple(r)


def toA(w):
    """image in S_2n: positions and values -n..-1,1..n relabelled 1..2n."""
    n = len(w)
    lab = lambda x: n + x if x > 0 else n + 1 + x
    full = [-w[-i - 1] for i in range(-n, 0)] + list(w)     # w(-n),...,w(-1),w(1),...,w(n)
    return tuple(lab(x) for x in full)


def amul(u, v):
    return tuple(u[x - 1] for x in v)


def ainv(u):
    r = [0] * len(u)
    for i, x in enumerate(u):
        r[x - 1] = i + 1
    return tuple(r)


def alen(w):
    m = len(w)
    return sum(1 for i in range(m) for j in range(i + 1, m) if w[i] > w[j])


def pinv(w):
    m = len(w)
    return frozenset((a, b) for a in range(1, m + 1) for b in range(a + 1, m + 1) if w[a - 1] > w[b - 1])


def involutions(m):
    """all involutions of S_m (including e) as tuples."""
    out = []

    def rec(free, cur):
        if not free:
            out.append(tuple(cur)); return
        a = free[0]
        rec(free[1:], cur)                       # a fixed
        for b in free[1:]:
            cur[a - 1], cur[b - 1] = b, a
            rec([x for x in free[1:] if x != b], cur)
            cur[a - 1], cur[b - 1] = a, b
    rec(list(range(1, m + 1)), list(range(1, m + 1)))
    return out


def bfs(gens, e):
    L = {e: 0}; q = deque([e])
    while q:
        x = q.popleft()
        for s in gens:
            y = smul(x, s)
            if y not in L:
                L[y] = L[x] + 1; q.append(y)
    return L


def gens(typ, n):
    e = tuple(range(1, n + 1)); g = []
    if typ == 'B':
        g.append((-1,) + e[1:])
    elif n >= 2:
        g.append((-2, -1) + e[2:])
    for i in range(n - 1):
        s = list(e); s[i], s[i + 1] = s[i + 1], s[i]; g.append(tuple(s))
    return g, e


def run(typ, n):
    g, e = gens(typ, n)
    L = bfs(g, e)                         # Coxeter length in W = B_n or D_n
    W = list(L)
    A = {w: toA(w) for w in W}
    LA = {w: alen(A[w]) for w in W}
    eA = A[e]
    F = dict(type='%s%d' % (typ, n), order=len(W), rank=len(g))
    neg = {w: sum(1 for x in w if x < 0) for w in W}
    invw = lambda w: sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])
    nsp = lambda w: sum(1 for i in range(n) for j in range(i + 1, n) if w[i] + w[j] < 0)
    if typ == 'D':
        F['d1_length_formula_fail'] = sum(1 for w in W if L[w] != invw(w) + nsp(w))
        F['d1_lA_eq_2lD_plus_neg_fail'] = sum(1 for w in W if LA[w] != 2 * L[w] + neg[w])
    if typ == 'B':
        F['b1_length_formula_fail'] = sum(1 for w in W if L[w] != invw(w) + neg[w] + nsp(w))
        F['b1_lA_eq_2lB_minus_neg_fail'] = sum(1 for w in W if LA[w] != 2 * L[w] - neg[w])
        F['centralizer_fail'] = sum(1 for w in W if tuple(2 * n + 1 - A[w][2 * n - i] for i in range(1, 2 * n + 1)) != A[w])
    # weak orders, all pairs
    Winv = {w: sinv(w) for w in W}
    mismatch = only_W = only_A = npref = 0
    for u in W:
        ui = Winv[u]
        for w in W:
            v = smul(ui, w)
            pW = L[w] == L[u] + L[v]
            pA = LA[w] == LA[u] + LA[v]
            npref += pW
            mismatch += pW != pA
            only_W += pW and not pA
            only_A += pA and not pW
    F['prefix_pairs_in_W'] = npref
    F['weak_order_mismatch'] = mismatch
    F['prefix_in_W_not_in_A'] = only_W
    F['prefix_in_A_not_in_W'] = only_A
    if n <= 3:
        m = 2 * n
        fixed = lambda P: sum(1 for (a, b) in P if a + b == m + 1)
        sign = 1 if typ == 'B' else -1               # B: |P| + |F(P)|,  D: |P| - |F(P)|
        bad = 0
        for x in W:
            Px = pinv(A[x])
            for y in W:
                P = Px & pinv(ainv(A[y]))            # PInv(x) & VInv(y)
                bad += L[smul(x, y)] != L[x] + L[y] - (len(P) + sign * fixed(P))
        F['b2_length_of_product_fail' if typ == 'B' else 'd2_length_of_product_fail'] = bad
    # ancestors in W, from the definition
    invsW = [u for u in W if u != e and smul(u, u) == e]
    ancW = {}; ties = 0
    for w in W:
        if w == e: continue
        pref = [u for u in invsW if L[w] == L[u] + L[smul(u, w)]]
        mx = max(L[u] for u in pref)
        top = [u for u in pref if L[u] == mx]
        ties += len(top) != 1
        ancW[w] = top[0]
    F['ancestor_ties_in_W'] = ties
    linvW = {e: 0}
    for w in sorted(W, key=lambda x: L[x]):
        if w != e:
            linvW[w] = 1 + linvW[smul(ancW[w], w)]
    # ancestors in S_2n, from the definition, over all involutions of S_2n
    invsA = [u for u in involutions(2 * n) if u != eA]
    lenI = {u: alen(u) for u in invsA}
    cacheA = {}

    def ancA(x):
        if x not in cacheA:
            lx = alen(x)
            pref = [u for u in invsA if lx == lenI[u] + alen(amul(u, x))]
            mx = max(lenI[u] for u in pref)
            top = [u for u in pref if lenI[u] == mx]
            assert len(top) == 1
            cacheA[x] = top[0]
        return cacheA[x]

    def linvA(x):
        k = 0
        while x != eA:
            x = amul(ancA(x), x); k += 1
        return k
    amis = lmis = 0; first = None
    for w in W:
        if w == e: continue
        if ancA(A[w]) != A[ancW[w]]:
            amis += 1
            if first is None or L[w] < L[first]:
                first = w
    LinvA = {w: linvA(A[w]) for w in W}
    lmis = sum(1 for w in W if LinvA[w] != linvW[w])
    F['ancestor_mismatch_W_vs_A'] = amis
    F['linv_mismatch_W_vs_A'] = lmis
    histW = [0] * (max(linvW.values()) + 1)
    for w in W: histW[linvW[w]] += 1
    F['hist_linv_W'] = histW
    F['max_linv_W'] = max(linvW.values())
    F['max_linv_A'] = max(LinvA.values())
    dispA = lambda x: max(abs(v - (i + 1)) for i, v in enumerate(x))
    F['linvA_gt_dispA'] = sum(1 for w in W if LinvA[w] > dispA(A[w]))
    F['linvW_gt_dispA'] = sum(1 for w in W if linvW[w] > dispA(A[w]))
    if first is not None:
        F['shortest_mismatch'] = dict(w=first, image_in_S2n=A[first], a_W=ancW[first],
                                      a_W_image=A[ancW[first]], a_A=ancA(A[first]),
                                      linv_W=linvW[first], linv_A=LinvA[first])
    print(json.dumps(F))


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    print(json.dumps(dict(python=sys.version.split()[0], script='typeB_cases.py', nmax=nmax)))
    for n in range(1, nmax + 1):
        run('B', n)
    for n in range(2, nmax + 1):
        run('D', n)
