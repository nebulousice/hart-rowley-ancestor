"""E8 without enumerating the group (696,729,600 elements; exhaustive search
is infeasible here).  Elements are permutations of the 240 roots.

1. All involutions of W(E8): closure of {w0(J) : J subset of S} under
   conjugation by the simple reflections (every involution is conjugate to
   some w0(J), Richardson 1982).  Stored with their left inversion bitmasks
   (120 positive roots -> two uint64 words) and lengths.
2. For an element w: involution prefixes = involutions iota with
   I_L(iota) subset of I_L(w); a(w) = the longest one (unique by Gaetz-Gao);
   l_inv by iteration; m(w) (least number of involution factors, length-
   additive) by memoised recursion.
3. Search: random elements, plus a greedy local search that tries to raise
   l_inv (multiply by a simple reflection on the right, keep if l_inv does
   not drop), plus the E7 witness family.  Every element visited gets the
   R4 test m(w) <= 8.

Usage: python e8_search.py SEED ROUNDS
"""
import sys, json, time
import numpy as np
from npcox import root_gens, coxeter_matrix

t0 = time.time()
gens, roots, pos = root_gens('E', 8)
gens = np.array(gens, dtype=np.int64)
pos = np.array(pos)
posidx = np.flatnonzero(pos)
assert len(roots) == 240 and len(posidx) == 120
rng0 = np.random.default_rng(31415)
MULT = rng0.integers(1, 2 ** 62, size=240, dtype=np.uint64) * np.uint64(2) + np.uint64(1)
E = np.arange(240, dtype=np.int64)


def key(p):
    return int((p.astype(np.uint64) * MULT).sum(dtype=np.uint64))


def mul(x, y):          # (xy)(c) = x[y[c]]
    return x[y]


def inv(x):
    q = np.empty_like(x); q[x] = np.arange(240); return q


def mask(x):
    """I_L(x) = {beta > 0 : x^-1 beta < 0} as an int (120 bits)."""
    xi = inv(x)
    neg = ~pos[xi[posidx]]
    return int(sum(1 << int(b) for b in np.flatnonzero(neg)))


def length(x):
    return int((~pos[x[posidx]]).sum())


# ---- involutions
def w0J(J):
    # longest element of W_J: BFS in W_J until no ascent
    x = E.copy()
    while True:
        up = False
        for s in J:
            y = mul(x, gens[s])
            if length(y) > length(x):
                x = y; up = True; break
        if not up:
            return x


invs = {}
frontier = []
import itertools
for k in range(1, 9):
    for J in itertools.combinations(range(8), k):
        x = w0J(J)
        kx = key(x)
        if kx not in invs:
            x = x.astype(np.int16); invs[kx] = x; frontier.append(x)
while frontier:
    nf = []
    for x in frontier:
        for s in range(8):
            y = gens[s][x[gens[s]]]          # s x s  (s = s^-1):  (s x s)(c) = s[x[s[c]]]
            ky = key(y)
            if ky not in invs:
                y = y.astype(np.int16); invs[ky] = y; nf.append(y)
    frontier = nf
INV = list(invs.values())
for x in INV[:50]:
    assert np.array_equal(x[x], E)
ninv = len(INV)
IM = np.zeros((ninv, 2), np.uint64)
IL = np.zeros(ninv, np.int64)
for j, x in enumerate(INV):
    m = mask(x)
    IM[j, 0] = np.uint64(m & (2 ** 64 - 1)); IM[j, 1] = np.uint64(m >> 64)
    IL[j] = length(x)
t1 = time.time()


def mask2(x):
    m = mask(x)
    return np.array([m & (2 ** 64 - 1), m >> 64], dtype=np.uint64)


def iprefixes(x):
    mw = mask2(x)
    ok = np.flatnonzero(((IM & mw[None, :]) == IM).all(axis=1))
    return ok


def ancestor(x):
    ok = iprefixes(x)
    if len(ok) == 0:
        return None
    L = IL[ok]
    best = ok[L == L.max()]
    assert len(best) == 1, 'ancestor not unique?!'
    return INV[best[0]]


def linv(x):
    k = 0
    while length(x) > 0:
        a = ancestor(x)
        x = mul(a, x)          # a^-1 x = a x
        k += 1
    return k


memo = {}


def m_le(x, k):
    lx = length(x)
    if lx == 0:
        return True
    if k <= 0:
        return False
    kk = (key(x), k)
    if kk in memo:
        return memo[kk]
    ok = iprefixes(x)
    res = False
    if len(ok) and IL[ok].max() == lx:        # x itself is an involution
        res = True
    else:
        for j in ok[np.argsort(-IL[ok])]:
            if m_le(mul(INV[j], x), k - 1):
                res = True; break
    memo[kk] = res
    return res


def from_word(word):
    x = E.copy()
    for s in word:
        x = mul(x, gens[s])
    return x


def main(seed, rounds):
    rng = np.random.default_rng(seed)
    stats = dict(involutions=ninv, sec_involutions=round(t1 - t0, 1))
    tested = 0; viol = []; best = (0, None); hist = {}
    # E7 witness (0-based Bourbaki) from BRIEF section 1.3, and the same word in E8
    wit7 = [6, 4, 5, 3, 4, 2, 3, 0, 2, 1, 3, 4, 5, 6]
    starts = [from_word(wit7)]
    for _ in range(rounds):
        starts.append(from_word(rng.integers(0, 8, size=int(rng.integers(10, 121))).tolist()))
    for x in starts:
        cur = x; lc = linv(cur)
        for step in range(40):                 # greedy local search on l_inv
            s = int(rng.integers(0, 8))
            y = mul(cur, gens[s])
            ly = linv(y)
            if ly >= lc:
                cur, lc = y, ly
            if lc >= 9:
                break
        for z in (x, cur):
            l = linv(z)
            hist[l] = hist.get(l, 0) + 1
            tested += 1
            if l > 8 and not m_le(z, 8):
                viol.append(z)
            if l > best[0]:
                best = (l, z)
    stats.update(tested=tested, hist_linv=dict(sorted(hist.items())), max_linv_found=best[0],
                 R4_violations=len(viol), seconds=round(time.time() - t0, 1), seed=seed)
    print(json.dumps(stats), flush=True)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
