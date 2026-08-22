"""Score affinity predictions on de-novo protein->small-molecule binders.

Joins reference/experimental_reference_ground_truth.csv with every predictions/*_predictions.csv
(by id) and evaluates:
  regression   on the 37 binders: does a method rank measured pKd? (Spearman/Pearson/RMSE + CI)
  specificity  on all 47: can a method tell the cognate ligand from the wrong one? (AUROC)
  combine      leave-one-out CV: do methods combine to beat the best single one?
"""

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOWER_BETTER = {"protenix-gpde"}                    # lower value = tighter (flip for rho/AUROC)
PK_METHODS = {"boltz-2", "nesso-1"}                 # predicted_affinity in pK units -> RMSE meaningful
COMBO = ("boltz-2", "nesso-1", "protenix-ligiptm", "clogp")


def _rd(path):
    return list(csv.DictReader(path.read_text().splitlines()))


def load():
    truth = {r["id"]: r for r in _rd(ROOT / "reference/experimental_reference_ground_truth.csv")}
    rows = {i: {"id": i, "is_binder": int(t["is_binder"]),
                "pKd": float(t["experimental_pKD"]) if t["experimental_pKD"] else None}
            for i, t in truth.items()}
    methods = []
    for p in sorted((ROOT / "predictions").glob("*_predictions.csv")):
        for r in _rd(p):
            m = r["method"]
            if m not in methods:
                methods.append(m)
            v = r["predicted_affinity"]
            rows[r["id"]][m] = float(v) if v.strip() else None
    return list(rows.values()), methods


def _rank(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for i, k in enumerate(o):
        r[k] = i
    return r


def spearman(xs, ys):
    rx, ry = _rank(xs), _rank(ys); n = len(xs)
    return 1 - 6 * sum((a - b) ** 2 for a, b in zip(rx, ry)) / (n * (n * n - 1))


def pearson(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5; sy = sum((y - my) ** 2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx and sy else float("nan")


def rmse(xs, ys):
    return (sum((x - y) ** 2 for x, y in zip(xs, ys)) / len(xs)) ** 0.5


def auroc(scores, labels):
    """P(score(pos) > score(neg)); labels 1=positive. Ties count 0.5."""
    pos = [s for s, l in zip(scores, labels) if l == 1]
    neg = [s for s, l in zip(scores, labels) if l == 0]
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def bootstrap_ci(xs, ys, stat, n_boot=2000, seed=0):
    import random
    rng = random.Random(seed); n = len(xs); out = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        try:
            out.append(stat([xs[i] for i in idx], [ys[i] for i in idx]))
        except (ZeroDivisionError, ValueError):
            pass
    out.sort()
    return out[int(0.025 * len(out))], out[int(0.975 * len(out))]


def regression(rows, methods):
    b = [r for r in rows if r["is_binder"] == 1 and r["pKd"] is not None]
    print(f"\n== Regression: rank measured pKd (binders, n={len(b)}) ==")
    print(f"{'method':20s}{'rho [95% CI]':>22s}{'pearson':>9s}{'RMSE':>7s}")
    for m in methods:
        v = [(r[m], r["pKd"]) for r in b if r.get(m) is not None]
        if len(v) < 3:
            continue
        xs, ys = map(list, zip(*v)); sgn = -1 if m in LOWER_BETTER else 1
        rho = sgn * spearman(xs, ys)
        lo, up = bootstrap_ci([sgn * x for x in xs], ys, spearman)
        rm = f"{rmse(xs, ys):.2f}" if m in PK_METHODS else "-"
        print(f"{m:20s}{f'{rho:+.2f} [{lo:+.2f},{up:+.2f}]':>22s}{sgn * pearson(xs, ys):>+9.2f}{rm:>7s}")


def specificity(rows, methods):
    print(f"\n== Specificity: cognate vs wrong ligand (AUROC, n={len(rows)}) ==")
    print(f"{'method':20s}{'AUROC [95% CI]':>22s}")
    for m in methods:
        v = [(r[m], r["is_binder"]) for r in rows if r.get(m) is not None]
        sc, lab = map(list, zip(*v)); sgn = -1 if m in LOWER_BETTER else 1
        a = auroc([sgn * x for x in sc], lab)
        lo, up = bootstrap_ci([sgn * x for x in sc], lab, auroc)
        print(f"{m:20s}{f'{a:.2f} [{lo:.2f},{up:.2f}]':>22s}")


def combine(rows):
    """Leave-one-out CV: do methods combine to beat the best single one? Small n, so LOO only."""
    import numpy as np
    from sklearn.linear_model import LinearRegression, LogisticRegression
    from sklearn.preprocessing import StandardScaler

    def loo(sub, y, model):
        X = np.array([[r[f] for f in COMBO] for r in sub]); y = np.array(y)
        pred = np.zeros(len(sub))
        for i in range(len(sub)):
            tr = [j for j in range(len(sub)) if j != i]
            sc = StandardScaler().fit(X[tr])
            m = model().fit(sc.transform(X[tr]), y[tr])
            pred[i] = (m.predict_proba(sc.transform(X[i:i + 1]))[0, 1]
                       if hasattr(m, "predict_proba") else m.predict(sc.transform(X[i:i + 1]))[0])
        return pred

    ok = [r for r in rows if all(r.get(f) is not None for f in COMBO)]
    b = [r for r in ok if r["is_binder"] == 1 and r["pKd"] is not None]
    print(f"\n== Combine ({'+'.join(COMBO)}), leave-one-out CV ==")
    p = loo(b, [r["pKd"] for r in b], LinearRegression)
    print(f"regression LOO Spearman rho = {spearman(list(p), [r['pKd'] for r in b]):+.2f}  (n={len(b)})")
    pc = loo(ok, [r["is_binder"] for r in ok], lambda: LogisticRegression(max_iter=1000))
    print(f"specificity LOO AUROC       = {auroc(list(pc), [r['is_binder'] for r in ok]):.2f}  (n={len(ok)})")


def main():
    argparse.ArgumentParser().parse_args()
    rows, methods = load()
    regression(rows, methods)
    specificity(rows, methods)
    combine(rows)


if __name__ == "__main__":
    main()
