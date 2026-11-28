"""First values for Hart-Rowley Conjecture 1.3 (arXiv:2502.00777), reproduced
from the definitions.  Run:  python rows.py        (about 10 s)

Prints, for every irreducible finite Coxeter group of order < 100,000 (the
source's stated verification range) the distribution of the involution length
l_inv, its maximum, and the number of elements with l_inv > rank.  Also
re-checks the worked examples in the source and in trureturing PR #17548.

Conventions: see coxlib.py.  Exact arithmetic throughout.
"""
import itertools
import sys
import time
from collections import Counter
from coxlib import Cox, ancestors_bruteforce, coxeter_graph

GROUPS = ([('A', n) for n in range(1, 8)] + [('B', n) for n in range(2, 7)] +
          [('D', n) for n in range(4, 7)] +
          [('E', 6), ('F', 4), ('G', 2), ('H', 3), ('H', 4)] +
          [('I', m) for m in (5, 7, 8, 9, 10, 12)])
BRUTE = {('A', 2), ('A', 3), ('A', 4), ('B', 2), ('B', 3), ('D', 4), ('H', 3), ('F', 4), ('G', 2), ('I', 7)}


def name(typ, n):
    return 'I2(%d)' % n if typ == 'I' else '%s%d' % (typ, n)


def w1(word):
    return ' '.join('r%d' % (s + 1) for s in word)


def main():
    print('type   |W|     N  rank  max l_inv  #(l_inv>rank)  unique-ancestor  distribution of l_inv (0,1,2,...)')
    fails = []
    for typ, n in GROUPS:
        G = Cox(typ, n, need_left=False)
        anc, rem, linv = G.ancestors(check_unique=True)       # raises if Conj. 1.2 fails
        c = Counter(linv)
        hist = [c[k] for k in range(max(c) + 1)]
        bad = sum(v for k, v in c.items() if k > G.rank)
        tag = ''
        if (typ, n) in BRUTE:
            a2, r2, l2 = ancestors_bruteforce(G)
            assert (a2, r2, l2) == (anc, rem, linv)
            tag = '  [= brute force]'
        if typ in 'ABD' and n <= 6:
            H = Cox(typ, n, need_left=False, model='classical')
            assert Counter(H.ancestors()[2]) == c
            tag += '  [= classical model]'
        print('%-6s %-6d %-3d %-4d %-9d %-14d %-16s %s%s' % (
            name(typ, n), G.order, G.N, G.rank, max(linv), bad, 'yes', hist, tag))
        if bad:
            i = min((i for i in range(G.order) if linv[i] > G.rank), key=lambda i: G.length[i])
            fails.append((typ, n, G, i))
    print()
    print('Shortest elements with l_inv > rank (generators r1, r2, ... = nodes 0, 1, ... of coxlib.cartan):')
    for typ, n, G, i in fails:
        dec = ' | '.join(w1(G.word(a)) for a in G.decomposition(i))
        print('  %-3s w = %s   (length %d)   ancestor decomposition: %s   -> l_inv = %d > %d' % (
            name(typ, n), w1(G.word(i)), G.length[i], dec, G.linv[i], G.rank))

    print()
    print('Checks of worked examples:')
    # (1) trureturing PR #17548: S4, w = 2431 = r1 r3 r2 r3: prefix count 3, suffix count 2
    G = Cox('A', 3)
    G.ancestors()
    w = G.from_word([0, 2, 1, 2])
    print('  A3  w = r1 r3 r2 r3 (2431):  l_inv =', G.linv[w], ', suffix version (l_inv of w^-1) =', G.linv[G.inv[w]],
          ' expected 3, 2:', (G.linv[w], G.linv[G.inv[w]]) == (3, 2))
    # (2) Hart-Rowley Section 3: A6, r6 r3 r2 r1 r4 r5 -> (r3 r6)(r2 r4)(r1 r5); suffix (r3)(r2 r4 r6)(r1 r5)
    G = Cox('A', 6)
    G.ancestors()
    w = G.from_word([5, 2, 1, 0, 3, 4])
    dec = [sorted(G.word(a)) for a in G.decomposition(w)]
    sdec = [sorted(G.word(a)) for a in G.decomposition(G.inv[w])][::-1]
    print('  A6  w = r6 r3 r2 r1 r4 r5:  prefix decomposition', dec, ' suffix decomposition', sdec,
          ' as in the source:', dec == [[2, 5], [1, 3], [0, 4]] and sdec == [[2], [1, 3, 5], [0, 4]])
    # (3) Hart-Rowley Prop. 2.2 / 2.4: Coxeter elements: a(w) = product of D(w); l_inv = path length
    ok = True
    for typ, n in [('A', 4), ('A', 5), ('B', 4), ('D', 5), ('E', 6), ('F', 4), ('H', 3), ('H', 4)]:
        G = Cox(typ, n)
        G.ancestors()
        edges = coxeter_graph(typ, n)
        adj = {i: set() for i in range(G.rank)}
        for i, j in edges:
            adj[i].add(j)
            adj[j].add(i)
        seen = set()
        for perm in itertools.permutations(range(G.rank)):
            w = G.from_word(perm)
            if w in seen:
                continue
            seen.add(w)
            pos = {s: k for k, s in enumerate(perm)}
            # longest path r_1,...,r_m in the Coxeter graph with positions increasing
            best = {}
            for s in sorted(range(G.rank), key=lambda s: -pos[s]):
                best[s] = 1 + max([best[t] for t in adj[s] if pos[t] > pos[s]], default=0)
            ok &= (G.linv[w] == max(best.values()))
            ok &= (sorted(G.word(G.anc[w])) == sorted(G.DL(w)))
    print('  Coxeter elements (A4, A5, B4, D5, E6, F4, H3, H4): l_inv = path length and a(w) = prod D(w):', ok)
    return 0


if __name__ == '__main__':
    t = time.time()
    rc = main()
    print('\n(%.1f s; python %s)' % (time.time() - t, sys.version.split()[0]))
    sys.exit(rc)
