"""R4 survey: m(w) = least k such that w = iota_1 ... iota_k with every iota_j
an involution and l(w) = sum l(iota_j)  (length-additive involution
factorization; m(e) = 0).  Corrected conjecture R4: max_W m <= rank(W).

Method A ('full'): exact m for every element by Bellman-Ford relaxation
  started from the greedy upper bound l_inv:  m(w) <- min(m(w), 1 + m(iota w))
  over all involution prefixes iota of w (prefix test: I_L(iota) subset of
  I_L(w), bitmasks over positive roots), until stationary.  At the fixed
  point m is exact (each value is realised by a factorization; and induction
  on length along an optimal factorization gives the lower bound).
Method B ('targeted'): exact decision m(w) <= k by memoised recursion, run
  only on the elements with l_inv(w) > rank (all others satisfy
  m <= l_inv <= rank), plus an exact computation of m on candidate maximisers
  (Coxeter elements) for the lower bound.

Usage: python r4.py TYPE N full|targeted|both
"""
import sys, json, time
import numpy as np
from npcox import NpGroup


def setup(typ, n):
    G = NpGroup(typ, n, 'roots')
    G.ancestors()
    G.left_inv_masks()
    invs = np.flatnonzero(G.isinv)
    invs = invs[invs != 0]
    G.invs = invs
    G.invmasks = G.masks[invs]
    G.invwords = [G.word(int(i)) for i in invs]
    return G


def selftest_prefix(G, sample=300, seed=1):
    """mask prefix test == length test on random pairs (iota, w)."""
    rng = np.random.default_rng(seed)
    bad = 0
    W = rng.integers(0, G.order, size=sample)
    for j, iota in enumerate(G.invs[:200]):
        mk = G.invmasks[j]
        test = ((G.masks[W] & mk[None, :]) == mk[None, :]).all(axis=1)
        # length test: l(w) = l(iota) + l(iota w); iota w via words
        for w, tv in zip(W, test):
            x = int(G.inv[w])
            for s in G.invwords[j]:
                x = int(G.R[s][x])
            r = int(G.inv[x])                    # r = iota^-1 w
            lt = G.length[w] == G.length[iota] + G.length[r]
            bad += int(lt != tv)
    return bad


def prefix_pairs(G, j):
    mk = G.invmasks[j]
    Wj = np.flatnonzero(((G.masks & mk[None, :]) == mk[None, :]).all(axis=1))
    x = G.inv[Wj]
    for s in G.invwords[j]:
        x = G.R[s][x]
    return Wj, G.inv[x]


def full(G):
    m = G.linv.astype(np.int16).copy()
    sweeps = 0
    npairs = 0
    while True:
        changed = 0
        for j in range(len(G.invs)):
            Wj, Rj = prefix_pairs(G, j)
            if sweeps == 0:
                npairs += len(Wj)
            c = 1 + m[Rj]
            better = c < m[Wj]
            if better.any():
                m[Wj[better]] = c[better]
                changed += int(better.sum())
        sweeps += 1
        if changed == 0:
            break
    return m, sweeps, npairs


class Targeted:
    def __init__(self, G):
        self.G = G
        self.memo = {}

    def iprefixes(self, w):
        G = self.G
        mw = G.masks[w]
        ok = np.flatnonzero(((G.invmasks & mw[None, :]) == G.invmasks).all(axis=1))
        out = []
        for j in ok:
            x = int(G.inv[w])
            for s in G.invwords[j]:
                x = int(G.R[s][x])
            out.append((int(G.invs[j]), int(G.inv[x])))
        # longest prefixes first (heuristic)
        out.sort(key=lambda p: -G.length[p[0]])
        return out

    def le(self, w, k):
        """m(w) <= k ?"""
        G = self.G
        if w == 0:
            return True
        if k <= 0:
            return False
        if G.linv[w] <= k or G.isinv[w]:
            return True
        key = (w, k)
        if key in self.memo:
            return self.memo[key]
        res = False
        for iota, r in self.iprefixes(w):
            if self.le(r, k - 1):
                res = True
                break
        self.memo[key] = res
        return res

    def exact(self, w):
        k = 0
        while not self.le(w, k):
            k += 1
        return k


def coxeter_elements(G):
    """All products of the simple reflections in some order (each once)."""
    import itertools
    seen = set()
    for perm in itertools.permutations(range(G.r)):
        x = 0
        for s in perm:
            x = int(G.R[s][x])
        seen.add(x)
    return sorted(seen)


def main(typ, n, mode):
    t0 = time.time()
    G = setup(typ, n)
    rank = G.r
    out = dict(type=typ + str(n), order=G.order, rank=rank, involutions=len(G.invs),
               max_linv=int(G.linv.max()), n_linv_gt_rank=int((G.linv > rank).sum()))
    if G.order <= 60000:
        out['prefix_selftest_mismatches'] = selftest_prefix(G)
    if mode in ('full', 'both'):
        m, sweeps, npairs = full(G)
        out.update(full_max_m=int(m.max()), full_hist_m=np.bincount(m).tolist(), full_sweeps=sweeps,
                   full_prefix_pairs=npairs, full_m_gt_rank=int((m > rank).sum()))
        wit = int(np.flatnonzero(m == m.max())[0])
        out['full_shortest_maximiser_word_0based'] = G.word(wit)
        np.save('out/R4_%s%d_m.npy' % (typ, n), m)
    if mode in ('targeted', 'both'):
        T = Targeted(G)
        cand = np.flatnonzero(G.linv > rank)
        viol = [int(w) for w in cand if not T.le(int(w), rank)]
        out['targeted_checked'] = len(cand)
        out['targeted_violations'] = len(viol)
        out['targeted_violation_words'] = [G.word(w) for w in viol[:5]]
        cox = coxeter_elements(G)
        cm = [T.exact(c) for c in cox]
        out['coxeter_elements'] = len(cox)
        out['coxeter_m_values'] = sorted(set(cm))
        # exact m on the targeted set (distribution), if small
        if len(cand) <= 20000:
            ex = [T.exact(int(w)) for w in cand]
            out['targeted_exact_m_hist_on_linv_gt_rank'] = np.bincount(ex).tolist() if ex else []
            out['targeted_max_m_on_checked_set_and_coxeter'] = max([max(cm)] + ex) if ex else max(cm)
        else:
            out['targeted_max_m_lower_bound'] = max(cm)
        if mode == 'both':
            out['agree_full_vs_targeted_on_cand'] = bool(all(T.exact(int(w)) == int(np.load('out/R4_%s%d_m.npy' % (typ, n))[w]) for w in cand[:3000]))
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
