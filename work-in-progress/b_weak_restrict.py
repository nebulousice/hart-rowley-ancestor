"""B_n inside S_{2n} (signed permutations as permutations of -n..-1,1..n):
check (1) l_A(w) = 2 l_B(w) - neg(w) for all w, and (2) u <=_R w in B_n iff
u <=_R w in S_{2n}, for all pairs (n <= 4) / random pairs (n = 5, 6).
This is the input of the reduction a_B(w) = a_A(w) (LAB.md, section 2).
"""
import sys, json, random
from coxsp import Group

def invA(x):
    m = len(x); return sum(1 for i in range(m) for j in range(i + 1, m) if x[i] > x[j])

out = {}
for n in (2, 3, 4, 5, 6):
    G = Group('B', n)
    lA = [invA(x) for x in G.elems]
    neg = [sum(1 for v in G.signed(i) if v < 0) for i in range(G.order)]
    bad1 = sum(1 for i in range(G.order) if lA[i] != 2 * G.length[i] - neg[i])
    def uinv_w(u, w):
        return G.mul(G.inv[u], w)
    pairs = [(u, w) for u in range(G.order) for w in range(G.order)] if n <= 4 else \
            [(random.Random(k).randrange(G.order), random.Random(k + 10**6).randrange(G.order)) for k in range(200000)]
    # bias: also include genuine prefixes (u = prefix of w via random descents)
    rng = random.Random(7)
    for _ in range(50000 if n > 4 else 0):
        w = rng.randrange(G.order); u = w
        for _ in range(rng.randrange(0, 6)):
            d = [s for s in range(G.r) if G.length[G.R[s][u]] < G.length[u]]
            if not d: break
            u = G.R[rng.choice(d)][u]
        pairs.append((u, w))
    bad2 = 0; nprefB = 0
    for u, w in pairs:
        v = uinv_w(u, w)
        pB = G.length[w] == G.length[u] + G.length[v]
        pA = lA[w] == lA[u] + lA[v]
        nprefB += pB
        bad2 += (pB != pA)
    out['B%d' % n] = dict(length_formula_failures=bad1, pairs=len(pairs), B_prefix_pairs=nprefB, weak_order_mismatches=bad2)
print(json.dumps(out))
