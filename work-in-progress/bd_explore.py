"""Types B, D exploration: the least admissible statistic G (Bruhat-monotone,
G(w) >= 1 + G(rem w)) restricted to bigrassmannian elements, and the
'join-irreducible envelope'  f*(x) = max{ G(beta) : beta bigrassmannian, beta <=_B x }.
Bruhat order: x >= beta computed by upward BFS through reflections (generic,
valid in every Coxeter group).
Usage: python bd_explore.py B 4
"""
import sys, json, time
from collections import Counter
from coxsp import Group


def upsets(Gp, betas):
    """For each beta, the boolean list [x >= beta]."""
    T = [Gp.elems[t] for t in Gp.T]
    out = {}
    for b in betas:
        seen = bytearray(Gp.order)
        seen[b] = 1
        stack = [b]
        while stack:
            x = stack.pop()
            ex = Gp.elems[x]
            lx = Gp.length[x]
            for t in T:
                k = Gp.index[tuple(ex[c] for c in t)]
                if not seen[k] and Gp.length[k] > lx:
                    seen[k] = 1
                    stack.append(k)
        out[b] = seen
    return out


def main(typ, n):
    t0 = time.time()
    Gp = Group(typ, n)
    Gp.ancestors()
    Gp.reflections()
    G = Gp.envelope_G()
    bigr = [w for w in range(1, Gp.order) if all(len(d) == 1 for d in Gp.descents(w))]
    ups = upsets(Gp, bigr)
    fstar = [0] * Gp.order
    for b in bigr:
        gb = G[b]
        s = ups[b]
        for x in range(Gp.order):
            if s[x] and gb > fstar[x]:
                fstar[x] = gb
    bad = sum(1 for w in range(1, Gp.order) if not fstar[Gp.rem[w]] < fstar[w])
    eqG = sum(1 for w in range(Gp.order) if fstar[w] == G[w])
    print('%s%d order %d ties %d max l_inv %d max G %d  #bigrassmannian %d  f*: max %d, step violations %d, f*==G on %d/%d  (%.1fs)' % (
        typ, n, Gp.order, Gp.ties, max(Gp.linv), max(G), len(bigr), max(fstar), bad, eqG, Gp.order, time.time() - t0), flush=True)
    rows = []
    for b in sorted(bigr, key=lambda b: (G[b], Gp.length[b])):
        DL, DR = Gp.descents(b)
        rows.append((G[b], Gp.length[b], Gp.signed(b) if typ != 'A' else Gp.elems[b], DL[0], DR[0], Gp.linv[b]))
    for r in rows:
        print('   G=%d len=%2d  w=%s  DL=%d DR=%d  linv=%d' % r)
    return Gp, G, fstar


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
