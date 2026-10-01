"""Score affinity predictions on de-novo protein->small-molecule binders.

Joins reference/experimental_reference_ground_truth.csv with every predictions/*_predictions.csv
(by id) and evaluates:
  regression   on the binders: does a method rank measured pKd? (Spearman/Pearson/RMSE + CI)
  specificity  on all pairs: can a method tell the correct ligand from the wrong one? (AUROC)
  combine      leave-one-out CV: do methods combine to beat the best single one?
"""

import argparse
import contextlib
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOWER_BETTER = {"protenix-gpde", "boltz-2-pae-min", "boltz-2-pae-mean",
                "protenix-pae-min", "protenix-pae-mean",
                "rowan-gnina", "rowan-aevplig"}   # log10(Kd): lower = tighter (flip for rho/AUROC)
PK_METHODS = {"boltz-2", "nesso-1"}                 # predicted_affinity in pK units -> RMSE meaningful
COMBO = ("boltz-2", "nesso-1", "protenix-ligiptm", "clogp")


def _rd(path):
    return list(csv.DictReader(path.read_text().splitlines()))


def load():
    truth = {r["id"]: r for r in _rd(ROOT / "reference/experimental_reference_ground_truth.csv")}
    lig = {r["id"]: r["smiles"] for r in _rd(ROOT / "reference/system_reference.csv")}
    rows = {i: {"id": i, "is_binder": int(t["is_binder"]),
                "pKd": float(t["experimental_pKD"]) if t["experimental_pKD"] else None,
                "ligand": lig.get(i, "")}
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
    return 1 - 6 * sum((a - b) ** 2 for a, b in zip(rx, ry, strict=True)) / (n * (n * n - 1))


def pearson(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True))
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5; sy = sum((y - my) ** 2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx and sy else float("nan")


def rmse(xs, ys):
    return (sum((x - y) ** 2 for x, y in zip(xs, ys, strict=True)) / len(xs)) ** 0.5


def auroc(scores, labels):
    """P(score(pos) > score(neg)); labels 1=positive. Ties count 0.5."""
    pos = [s for s, l in zip(scores, labels, strict=True) if l == 1]
    neg = [s for s, l in zip(scores, labels, strict=True) if l == 0]
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def bootstrap_ci(xs, ys, stat, n_boot=2000, seed=0):
    import random
    rng = random.Random(seed); n = len(xs); out = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        with contextlib.suppress(ZeroDivisionError, ValueError):
            out.append(stat([xs[i] for i in idx], [ys[i] for i in idx]))
    out.sort()
    return out[int(0.025 * len(out))], out[int(0.975 * len(out))]


def regression(rows, methods):
    b = [r for r in rows if r["is_binder"] == 1 and r["pKd"] is not None]
    print(f"\n== Regression (POOLED; confounded by ligand composition, see within-ligand below): "
          f"rank pKd (binders, n={len(b)}) ==")
    print(f"{'method':20s}{'rho [95% CI]':>22s}{'pearson':>9s}{'RMSE':>7s}")
    for m in methods:
        v = [(r[m], r["pKd"]) for r in b if r.get(m) is not None]
        if len(v) < 3:
            continue
        xs, ys = map(list, zip(*v, strict=True)); sgn = -1 if m in LOWER_BETTER else 1
        rho = sgn * spearman(xs, ys)
        lo, up = bootstrap_ci([sgn * x for x in xs], ys, spearman)
        rm = f"{rmse(xs, ys):.2f}" if m in PK_METHODS else "-"
        print(f"{m:20s}{f'{rho:+.2f} [{lo:+.2f},{up:+.2f}]':>22s}{sgn * pearson(xs, ys):>+9.2f}{rm:>7s}")


def within_ligand_rho(rows, m, kmin=3):
    """Per-ligand Spearman(metric, pKd) among the de-novo proteins binding that ligand, then the
    median across ligand groups. Controls for ligand composition: a ligand-only property (cLogP, MW)
    is constant within a group, so it drops out. Returns (median, [per-group rho], n_groups)."""
    import collections
    import statistics
    sgn = -1 if m in LOWER_BETTER else 1
    groups = collections.defaultdict(list)
    for r in rows:
        if r["is_binder"] == 1 and r["pKd"] is not None and r.get(m) is not None:
            groups[r["ligand"]].append((sgn * r[m], r["pKd"]))
    per = [spearman([x for x, _ in g], [y for _, y in g])
           for g in groups.values() if len(g) >= kmin and len({x for x, _ in g}) > 1]
    if len(per) < 2:
        return None, per, len(per)
    return statistics.median(per), per, len(per)


def regression_within_ligand(rows, methods, kmin=3):
    import random
    import statistics
    print("\n== Regression (within-ligand): rank pKd among proteins binding the SAME ligand ==")
    print("   controls for ligand composition (a ligand-only score is constant per group, so it "
          f"drops out); median rho over ligands with >= {kmin} binders, direction-corrected")
    print(f"{'method':20s}{'median rho [95% CI]':>24s}{'n_lig':>7s}")
    for m in methods:
        med, per, n = within_ligand_rho(rows, m, kmin)
        if med is None:
            continue
        rng = random.Random(0); boots = []
        for _ in range(2000):
            boots.append(statistics.median([per[rng.randrange(n)] for _ in range(n)]))
        boots.sort()
        lo, hi = boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]
        print(f"{m:20s}{f'{med:+.2f} [{lo:+.2f},{hi:+.2f}]':>24s}{n:>7d}")


def specificity(rows, methods):
    print(f"\n== Specificity: correct vs wrong ligand (AUROC, n={len(rows)}) ==")
    print(f"{'method':20s}{'AUROC [95% CI]':>22s}")
    for m in methods:
        v = [(r[m], r["is_binder"]) for r in rows if r.get(m) is not None]
        sc, lab = map(list, zip(*v, strict=True)); sgn = -1 if m in LOWER_BETTER else 1
        a = auroc([sgn * x for x in sc], lab)
        lo, up = bootstrap_ci([sgn * x for x in sc], lab, auroc)
        print(f"{m:20s}{f'{a:.2f} [{lo:.2f},{up:.2f}]':>22s}")


def combine_loo(rows):
    """Leave-one-out CV predictions for the fixed COMBO feature set. Returns (reg_pred, reg_truth),
    (spec_pred, spec_truth). Small n, so LOO only."""
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
        return list(pred)

    ok = [r for r in rows if all(r.get(f) is not None for f in COMBO)]
    b = [r for r in ok if r["is_binder"] == 1 and r["pKd"] is not None]
    reg = (loo(b, [r["pKd"] for r in b], LinearRegression), [r["pKd"] for r in b])
    spec = (loo(ok, [r["is_binder"] for r in ok], lambda: LogisticRegression(max_iter=1000)),
            [r["is_binder"] for r in ok])
    return reg, spec


def combine(rows, best="boltz-2-iptm"):
    import random
    (rp, ry), (sp, sy) = combine_loo(rows)
    print(f"\n== Combine ({'+'.join(COMBO)}), leave-one-out CV ==")
    print(f"regression LOO Spearman rho = {spearman(rp, ry):+.2f}  (n={len(ry)})")
    print(f"specificity LOO AUROC       = {auroc(sp, sy):.2f}  (n={len(sy)})")
    # is the combination actually better than the best single metric? paired bootstrap of the difference.
    ok = [r for r in rows if all(r.get(f) is not None for f in COMBO)]
    b = [r for r in ok if r["is_binder"] == 1 and r["pKd"] is not None and r.get(best) is not None]
    ref = [r[best] for r in b]
    rng = random.Random(0); n = len(b); diffs = []
    for _ in range(2000):
        idx = [rng.randrange(n) for _ in range(n)]
        c = [rp[i] for i in idx]; s = [ref[i] for i in idx]; t = [ry[i] for i in idx]
        diffs.append(spearman(c, t) - spearman(s, t))
    diffs.sort()
    lo, hi = diffs[int(0.025 * len(diffs))], diffs[int(0.975 * len(diffs))]
    print(f"combination minus best single ({best}) rho, 95% CI = [{lo:+.2f}, {hi:+.2f}]  (spans 0 = not better)")


def main():
    argparse.ArgumentParser().parse_args()
    rows, methods = load()
    regression(rows, methods)
    regression_within_ligand(rows, methods)
    specificity(rows, methods)
    combine(rows)


if __name__ == "__main__":
    main()
