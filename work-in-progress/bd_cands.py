"""Test candidate statistics f on B_n / D_n: (i) max f, (ii) f(rem w) < f(w)
for all w != e, (iii) l_inv <= f, (iv) Bruhat-monotone (checked on
reflection-covers x < x t), (v) comparison with G.
Usage: python bd_cands.py B 5
"""
import sys, json, time
from collections import Counter
from coxsp import Group


def f_neg_max(w):
    """positive entries: |w(i) - i|; negative entries: max(i, |w(i)|)."""
    out = 0
    for i, x in enumerate(w, 1):
        out = max(out, abs(x - i) if x > 0 else max(i, -x))
    return out


def f_neg_sum(w):
    out = 0
    for i, x in enumerate(w, 1):
        out = max(out, abs(x - i) if x > 0 else (i - x + 1) // 2)
    return out


CANDS = {'neg_max': f_neg_max}


def test(typ, n, cands=CANDS, G=None):
    t0 = time.time()
    Gp = Group(typ, n)
    Gp.ancestors()
    Gp.reflections()
    Tel = [Gp.elems[t] for t in Gp.T]
    S = [Gp.signed(i) for i in range(Gp.order)]
    res = {}
    for name, f in cands.items():
        F = [f(s) for s in S]
        step = [w for w in range(1, Gp.order) if not F[Gp.rem[w]] < F[w]]
        weak = sum(1 for w in range(Gp.order) if Gp.linv[w] > F[w])
        mono = 0
        for w in range(Gp.order):
            x = Gp.elems[w]
            for t in Tel:
                k = Gp.index[tuple(x[c] for c in t)]
                if Gp.length[k] < Gp.length[w] and F[k] > F[w]:
                    mono += 1
        res[name] = dict(max=max(F), step_violations=len(step), linv_gt_f=weak, monotone_violations=mono,
                         examples=[(S[w], S[Gp.rem[w]], F[w], F[Gp.rem[w]]) for w in step[:4]])
    print(typ + str(n), json.dumps(res), '(%.1fs)' % (time.time() - t0), flush=True)
    return Gp, res


if __name__ == '__main__':
    test(sys.argv[1], int(sys.argv[2]))
