"""Small-case checks for PROOF_typeA.md (hart-rowley-ancestor, HR-2).

Everything is computed from the definitions, in exact integer arithmetic,
in pure Python (no numpy, no code shared with plans/C1 or scripts/).

Conventions (PROOF_typeA.md section 1): a permutation w of {1..N} is the
tuple (w(1), ..., w(N)); (u v)(i) = u(v(i)); length = number of inversions;
u <=_R w  iff  l(w) = l(u) + l(u^-1 w).

Usage:  python small_cases.py NMAX [tables]
  NMAX   : run every check for N = 2..NMAX (6 takes seconds, 7 about a minute)
  tables : also print the full S_3 and S_4 tables of section 3
Output: one JSON line per N with failure counts (all must be 0).
"""
import sys, json, itertools


def mul(u, v):
    return tuple(u[v[i] - 1] for i in range(len(v)))


def inv(u):
    r = [0] * len(u)
    for i, x in enumerate(u):
        r[x - 1] = i + 1
    return tuple(r)


def length(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def vinv(w):
    """value inversions: pairs (a, b), a < b, with b standing left of a."""
    pos = inv(w)
    n = len(w)
    return frozenset((a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)
                     if pos[a - 1] > pos[b - 1])


def pinv(w):
    """position inversions: pairs (a, b), a < b, with w(a) > w(b)."""
    n = len(w)
    return frozenset((a, b) for a in range(1, n + 1) for b in range(a + 1, n + 1)
                     if w[a - 1] > w[b - 1])


def disp(w):
    return max(abs(x - (i + 1)) for i, x in enumerate(w))


def rank_leq(x, y):
    """x[i,k] <= y[i,k] for all i, k  (x[i,k] = #{j <= i : x(j) >= k})."""
    n = len(x)
    for k in range(1, n + 1):
        cx = cy = 0
        for i in range(n):
            cx += x[i] >= k
            cy += y[i] >= k
            if cx > cy:
                return False
    return True


def bruhat_closure(perms, L):
    """Bruhat order from the definition: transitive closure of x -> x t,
    t a transposition, l(x t) = l(x) + 1.  Returns dict x -> set of y >= x."""
    n = len(perms[0])
    up = {x: set() for x in perms}
    for x in perms:
        for i in range(n):
            for j in range(i + 1, n):
                y = list(x); y[i], y[j] = y[j], y[i]; y = tuple(y)
                if L[y] == L[x] + 1:
                    up[x].add(y)
    above = {}
    for x in sorted(perms, key=lambda p: -L[p]):
        s = {x}
        for y in up[x]:
            s |= above[y]
        above[x] = s
    return above


def lemma_K(w, p):
    """The involution of Lemma K for an extremal excedance position p
    (1-based).  Returns (iota, data) or raises AssertionError."""
    n = len(w)
    D = disp(w)
    v = w[p - 1]
    assert D >= 1 and v == p + D
    Lset = set(w[:p])
    assert all(x <= v for x in Lset)                       # (a)
    Rset = set(range(1, v + 1)) - Lset
    assert len(Rset) == D and all(r >= p + 1 - D for r in Rset)   # (a), (b)
    h = 0; c0 = None
    for c in range(v, 0, -1):
        h += 1 if c in Lset else -1
        if h == 0:
            c0 = c
            break
    assert c0 is not None                                  # (c)
    J = range(c0, v + 1)
    ls = [c for c in J if c in Lset]
    rs = [c for c in J if c in Rset]
    assert len(ls) == len(rs) >= 1
    assert all(r < l for r, l in zip(rs, ls))              # (d)
    io = list(range(1, n + 1))
    for r, l in zip(rs, ls):
        io[r - 1] = l; io[l - 1] = r
    io = tuple(io)
    return io, (Lset, Rset, c0, ls, rs)


def phi(w):
    n = len(w)
    return tuple(n + 1 - w[n - i] for i in range(1, n + 1))


def run(N, tables=False):
    perms = [tuple(p) for p in itertools.permutations(range(1, N + 1))]
    e = tuple(range(1, N + 1))
    L = {w: length(w) for w in perms}
    invs = [u for u in perms if u != e and mul(u, u) == e]
    VI = {w: vinv(w) for w in perms}
    F = dict(N=N, order=len(perms), involutions=len(invs))

    # ancestors from the definition (length additivity, all involutions)
    anc = {}; ties = 0
    ip = {}
    for w in perms:
        if w == e:
            continue
        pref = [u for u in invs if L[u] + L[mul(inv(u), w)] == L[w]]
        ip[w] = pref
        m = max(L[u] for u in pref)
        top = [u for u in pref if L[u] == m]
        ties += len(top) != 1
        anc[w] = top[0]
    F['ancestor_ties'] = ties
    rem = {w: mul(anc[w], w) for w in anc}       # a^-1 w = a w
    linv = {e: 0}
    for w in sorted(perms, key=lambda p: L[p]):
        if w != e:
            linv[w] = 1 + linv[rem[w]]
    hist = [0] * N
    for w in perms:
        hist[linv[w]] += 1
    F['hist_linv'] = hist
    F['max_linv'] = max(linv.values())

    # Lemma 2.1 (length formula) and Lemma 2.2 (prefix criterion), all pairs
    if N <= 5:
        PI = {w: pinv(w) for w in perms}
        bad_len = bad_pref = 0
        for x in perms:
            for y in perms:
                bad_len += L[mul(x, y)] != len(PI[x] ^ VI[y])
                isp = L[y] == L[x] + L[mul(inv(x), y)]
                bad_pref += isp != (VI[x] <= VI[y])
        F['lemma21_length_formula_fail'] = bad_len
        F['lemma22_prefix_criterion_fail'] = bad_pref

    # Lemma 2.4: rank-matrix monotonicity against Bruhat order from the definition
    if N <= 6:
        above = bruhat_closure(perms, L)
        bad_easy = bad_conv = 0
        for x in perms:
            for y in perms:
                b = y in above[x]
                r = rank_leq(x, y)
                bad_easy += b and not r      # direction proved and used in the proof
                bad_conv += r and not b      # converse (BB Thm 2.1.5), not used
        F['lemma24_bruhat_implies_rank_fail'] = bad_easy
        F['rank_implies_bruhat_fail_not_used'] = bad_conv
        leqB = lambda x, y: y in above[x]
    else:
        leqB = rank_leq                      # BB Thm 2.1.5 (checked above for N <= 6)

    # suffix => Bruhat-below; disp Bruhat-monotone (Lemma 4.1) on the pairs used
    bad_suffix = bad_dispmono = 0
    # GG1 instances: every involution prefix is Bruhat-below a(w)
    bad_gg1 = 0
    for w in anc:
        a = anc[w]
        for u in ip[w]:
            bad_gg1 += not leqB(u, a)
            s = mul(u, w)
            bad_suffix += not leqB(s, w)
            bad_dispmono += disp(s) > disp(w)
    F['gg1_fail'] = bad_gg1
    F['suffix_not_bruhat_below_fail'] = bad_suffix
    F['disp_suffix_fail'] = bad_dispmono

    # GG2 instances: all u, v <=_R w with u <=_B v  (N <= 5: all prefixes; else involution prefixes)
    bad_gg2 = n_gg2 = 0
    if N <= 5:
        for w in perms:
            pref = [u for u in perms if VI[u] <= VI[w]]
            for u in pref:
                for v in pref:
                    if leqB(u, v):
                        n_gg2 += 1
                        bad_gg2 += not leqB(mul(inv(v), w), mul(inv(u), w))
    else:
        for w in anc:
            for u in ip[w]:
                for v in ip[w]:
                    if leqB(u, v):
                        n_gg2 += 1
                        bad_gg2 += not leqB(mul(v, w), mul(u, w))
    F['gg2_instances'] = n_gg2
    F['gg2_fail'] = bad_gg2

    # Lemma K on every extremal excedance; drops through phi
    pairs = bad_K = 0
    for w in anc:
        D = disp(w)
        for kind, ww in (('exc', w), ('drop', phi(w))):
            for p in range(1, N + 1):
                if ww[p - 1] == p + D:
                    pairs += 1
                    try:
                        io, (Ls, Rs, c0, ls, rs) = lemma_K(ww, p)
                        v = ww[p - 1]
                        ok = mul(io, io) == e and io != e
                        ok &= L[io] + L[mul(io, ww)] == L[ww]          # prefix, by definition
                        ok &= io[v - 1] < v
                        ok &= all(a in Rs and b in Ls for (a, b) in pinv(io))   # (f)
                        # Step 2 of Proposition 4.3, replayed
                        a = anc[ww]; r = rem[ww]; iw = mul(io, ww)
                        ok &= leqB(io, a) and leqB(r, iw)
                        ok &= r[p - 1] != p + D
                        bad_K += not ok
                    except AssertionError:
                        bad_K += 1
    F['lemmaK_pairs'] = pairs
    F['lemmaK_fail'] = bad_K

    # phi commutes with a and rem
    F['phi_equivariance_fail'] = sum(1 for w in anc if anc[phi(w)] != phi(anc[w]))

    # Theorem A
    F['disp_rem_not_smaller'] = sum(1 for w in anc if disp(rem[w]) >= disp(w))
    F['linv_gt_disp'] = sum(1 for w in perms if linv[w] > disp(w))
    F['linv_gt_rank'] = sum(1 for w in perms if linv[w] > N - 1)

    # Proposition 4.6: c = s_1 ... s_{N-1} = (2, 3, ..., N, 1)
    c = tuple(list(range(2, N + 1)) + [1])
    chain_ok = True; x = c; k = 1
    while x != e:
        sk = list(e); sk[k - 1], sk[k] = sk[k], sk[k - 1]
        chain_ok &= anc[x] == tuple(sk)
        x = rem[x]; k += 1
    F['coxeter_element_linv'] = linv[c]
    F['coxeter_chain_ok'] = chain_ok

    print(json.dumps(F))

    if tables and N in (3, 4):
        print('table S_%d: w | VInv(w) | involution prefixes | a(w) | rem(w) | l_inv | disp(w) | disp(rem)' % N)
        s = lambda t: ''.join(map(str, t))
        for w in perms:
            if w == e:
                print('  %s | {} | - | - | - | 0 | 0 | -' % s(w))
                continue
            vi = ' '.join('%d%d' % ab for ab in sorted(VI[w]))
            print('  %s | %s | %s | %s | %s | %d | %d | %d' % (
                s(w), vi, ' '.join(s(u) for u in sorted(ip[w], key=lambda u: L[u])),
                s(anc[w]), s(rem[w]), linv[w], disp(w), disp(rem[w])))


if __name__ == '__main__':
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    tables = len(sys.argv) > 2 and sys.argv[2] == 'tables'
    print(json.dumps(dict(python=sys.version.split()[0], script='small_cases.py', nmax=nmax)))
    for N in range(2, nmax + 1):
        run(N, tables)
