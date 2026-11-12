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
import math
from pathlib import Path

import numpy as np

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
    """Midranks: tied values share the average of the ranks they span.

    The obvious version -- sort, then number 0..n-1 -- breaks ties by original
    row order, which invents an ordering the data does not have. On this set
    that was not academic: `interface-ligburial` has 42 tied values in 60 and
    `clogp` is constant within a ligand by construction, so the old ranking
    reported a perfect correlation (rho +1.00 where scipy gives +0.77) for
    vectors that carry almost no ordering at all.
    """
    a = np.asarray(v, dtype=float)
    order = np.argsort(a, kind="mergesort")       # stable, so ties group together
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(len(a), dtype=float)
    # average the ranks within each run of equal values
    sa = a[order]
    i = 0
    while i < len(sa):
        j = i + 1
        while j < len(sa) and sa[j] == sa[i]:
            j += 1
        if j - i > 1:
            ranks[order[i:j]] = (i + j - 1) / 2.0
        i = j
    return ranks


def spearman(xs, ys):
    """Pearson correlation of the midranks.

    NOT the 1 - 6*sum(d^2)/(n(n^2-1)) shortcut, which is only valid when no
    value is tied. With ties it overstates the correlation, which is how the
    bug above went unnoticed: the formula still returns a plausible number.
    """
    rx, ry = _rank(xs), _rank(ys)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = math.sqrt(float(rx @ rx) * float(ry @ ry))
    return float(rx @ ry) / d if d else float("nan")


def pearson(xs, ys):
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True))
    sx = sum((x - mx) ** 2 for x in xs) ** 0.5; sy = sum((y - my) ** 2 for y in ys) ** 0.5
    return cov / (sx * sy) if sx and sy else float("nan")


def rmse(xs, ys):
    return (sum((x - y) ** 2 for x, y in zip(xs, ys, strict=True)) / len(xs)) ** 0.5


def auroc(scores, labels):
    """P(score(pos) > score(neg)); labels 1=positive. Ties count 0.5.

    Rank form of the Mann-Whitney statistic: O(n log n) instead of the O(P*N)
    double loop it replaces, and identical to it because the midranks already
    give ties their half-credit. On this set that is 49x11 pair comparisons per
    bootstrap resample, which is where the runtime was going.
    """
    s = np.asarray(scores, dtype=float)
    y = np.asarray(labels)
    npos = int((y == 1).sum())
    nneg = int((y == 0).sum())
    if not npos or not nneg:
        return None
    r = _rank(s) + 1.0                              # _rank is 0-based
    return float((r[y == 1].sum() - npos * (npos + 1) / 2.0) / (npos * nneg))


def _midranks_2d(a):
    """Midranks along axis 1 for a (B, n) array, vectorised over all B rows.

    Same convention as _rank, but without a Python loop per row: the whole
    bootstrap is ranked in one pass. Tie runs are found by comparing each sorted
    value to its neighbour and taking a running max of run starts forwards and a
    running min of run ends backwards, which gives every element the average of
    the positions its run spans.
    """
    b, n = a.shape
    order = np.argsort(a, axis=1, kind="stable")
    srt = np.take_along_axis(a, order, axis=1)
    pos = np.broadcast_to(np.arange(n, dtype=float), (b, n))

    new_run = np.ones((b, n), dtype=bool)
    new_run[:, 1:] = srt[:, 1:] != srt[:, :-1]
    start = np.maximum.accumulate(np.where(new_run, pos, -np.inf), axis=1)
    end_marker = np.ones((b, n), dtype=bool)
    end_marker[:, :-1] = srt[:, :-1] != srt[:, 1:]
    end = np.minimum.accumulate(np.where(end_marker, pos, np.inf)[:, ::-1], axis=1)[:, ::-1]

    srt_ranks = (start + end) / 2.0
    ranks = np.empty_like(srt_ranks)
    np.put_along_axis(ranks, order, srt_ranks, axis=1)
    return ranks


def batch_spearman(x_boot, y_boot):
    """Spearman for every row pair of two (B, n) arrays, in one pass."""
    rx = _midranks_2d(x_boot)
    ry = _midranks_2d(y_boot)
    rx = rx - rx.mean(axis=1, keepdims=True)
    ry = ry - ry.mean(axis=1, keepdims=True)
    num = (rx * ry).sum(axis=1)
    den = np.sqrt((rx * rx).sum(axis=1) * (ry * ry).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def batch_auroc(s_boot, y_boot):
    """AUROC for every row pair of a (B, n) score array and label array."""
    r = _midranks_2d(s_boot) + 1.0
    pos = y_boot == 1
    npos = pos.sum(axis=1).astype(float)
    nneg = (y_boot == 0).sum(axis=1).astype(float)
    rsum = np.where(pos, r, 0.0).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (rsum - npos * (npos + 1) / 2.0) / (npos * nneg)
    return np.where((npos > 0) & (nneg > 0), out, np.nan)


# stat -> its batched twin. bootstrap_ci uses these when it can, which turns
# 2000 Python calls per metric into one vectorised pass.
_BATCHED = {}


def bootstrap_ci(xs, ys, stat, n_boot=2000, seed=0):
    """Percentile CI over n_boot resamples.

    All resample indices are drawn in one numpy call rather than n_boot*n
    separate `random.randrange` calls, which was 5.9 million of them per run and
    most of the script's runtime. `stat` is still called per resample because
    the callers pass different statistics; the draw is what was expensive.
    """
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    n = len(x)
    idx = np.random.default_rng(seed).integers(0, n, size=(n_boot, n))
    batched = _BATCHED.get(stat)
    if batched is not None:
        vals = batched(x[idx], y[idx])
        out = np.sort(vals[np.isfinite(vals)])
    else:
        acc = []
        for row in idx:
            with contextlib.suppress(ZeroDivisionError, ValueError):
                acc.append(stat(x[row], y[row]))
        out = np.sort(np.asarray([v for v in acc if v is not None and np.isfinite(v)]))
    if not len(out):
        return float("nan"), float("nan")
    return float(out[int(0.025 * len(out))]), float(out[int(0.975 * len(out))])


_BATCHED.update({spearman: batch_spearman, auroc: batch_auroc})


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


def per_ligand_rho(rows, m, kmin=3):
    """{ligand: Spearman(metric, pKd)} over the de-novo proteins binding that ligand, plus the
    pair count. Keyed by ligand, not positional, because two metrics are not scored on the same
    ligand set -- pairing them by list position silently compares different ligands."""
    import collections
    sgn = -1 if m in LOWER_BETTER else 1
    groups = collections.defaultdict(list)
    for r in rows:
        if r["is_binder"] == 1 and r["pKd"] is not None and r.get(m) is not None:
            groups[r["ligand"]].append((sgn * r[m], r["pKd"]))
    used = {k: g for k, g in groups.items()
            if len(g) >= kmin and len({x for x, _ in g}) > 1}
    per = {k: spearman([x for x, _ in g], [y for _, y in g]) for k, g in used.items()}
    return per, sum(len(g) for g in used.values())


def within_ligand_rho(rows, m, kmin=3):
    """Per-ligand Spearman then the median across ligand groups. Controls for ligand composition:
    a ligand-only property (cLogP, MW) is constant within a group, so it drops out.

    Returns (median, {ligand: rho}, n_groups, n_pairs). n_pairs is the number of measured
    (protein, ligand) pairs that actually reach the median -- NOT the 60 pairs in the benchmark.
    Most ligands here have one or two binders and contribute nothing to a within-ligand rank."""
    import statistics
    per, n_pairs = per_ligand_rho(rows, m, kmin)
    if len(per) < 2:
        return None, per, len(per), n_pairs
    return statistics.median(per.values()), per, len(per), n_pairs


def regression_within_ligand(rows, methods, kmin=3):
    import random
    import statistics
    print("\n== Regression (within-ligand): rank pKd among proteins binding the SAME ligand ==")
    print("   controls for ligand composition (a ligand-only score is constant per group, so it "
          f"drops out); median rho over ligands with >= {kmin} binders, direction-corrected")
    print(f"{'method':20s}{'median rho [95% CI]':>24s}{'n_lig':>7s}{'n_pairs':>9s}")
    for m in methods:
        med, per, n, n_pairs = within_ligand_rho(rows, m, kmin)
        if med is None:
            continue
        vals = list(per.values())
        rng = random.Random(0); boots = []
        for _ in range(2000):
            boots.append(statistics.median([vals[rng.randrange(n)] for _ in range(n)]))
        boots.sort()
        lo, hi = boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]
        print(f"{m:20s}{f'{med:+.2f} [{lo:+.2f},{hi:+.2f}]':>24s}{n:>7d}{n_pairs:>9d}")
    print("   n_lig is the number of ligands the median is taken over, n_pairs the measured "
          "points behind it.\n   Neither is 60: a ligand with fewer than "
          f"{kmin} binders cannot support a rank and is dropped.")


def separability(rows, methods, kmin=3, n_boot=4000):
    """Which metrics are actually distinguishable from the best one, and which are not.

    A table of point estimates reads as a ranking even when the intervals swamp the gaps -- the
    same objection that sinks most published score comparisons, and it bites hard here: the
    within-ligand median is taken over single-digit numbers of ligands.

    Paired bootstrap resampling LIGANDS (the independent unit, not pairs), restricted to the
    ligands BOTH metrics scored, so the two medians always describe the same groups.
    """
    import random
    import statistics

    per_metric = {}
    for m in methods:
        per, _ = per_ligand_rho(rows, m, kmin)
        if len(per) >= 2:
            per_metric[m] = per
    if len(per_metric) < 2:
        return

    best = max(per_metric, key=lambda m: statistics.median(per_metric[m].values()))
    bper = per_metric[best]
    bmed = statistics.median(bper.values())
    print("\n== Separability: is anything actually worse than the best metric? ==")
    print(f"   paired bootstrap over the ligand groups BOTH metrics scored (best has {len(bper)})")
    print(f"   best = {best} (median rho {bmed:+.2f})")
    print(f"{'method':20s}{'median rho':>12s}{'shared':>8s}"
          f"{'delta vs best [95% CI]':>28s}{'separable?':>12s}")
    for m, per in sorted(per_metric.items(), key=lambda kv: -statistics.median(kv[1].values())):
        med = statistics.median(per.values())
        if m == best:
            print(f"{m:20s}{med:>+12.2f}{len(per):>8d}{'-':>28s}{'(best)':>12s}")
            continue
        shared = sorted(set(per) & set(bper))
        if len(shared) < 2:
            print(f"{m:20s}{med:>+12.2f}{len(shared):>8d}{'too few shared ligands':>28s}{'-':>12s}")
            continue
        bv = [bper[k] for k in shared]; mv = [per[k] for k in shared]
        rng = random.Random(0); k = len(shared); deltas = []
        for _ in range(n_boot):
            idx = [rng.randrange(k) for _ in range(k)]
            deltas.append(statistics.median([bv[i] for i in idx])
                          - statistics.median([mv[i] for i in idx]))
        deltas.sort()
        lo, hi = deltas[int(0.025 * n_boot)], deltas[int(0.975 * n_boot)]
        d = statistics.median(bv) - statistics.median(mv)
        print(f"{m:20s}{med:>+12.2f}{k:>8d}"
              f"{f'{d:+.2f} [{lo:+.2f},{hi:+.2f}]':>28s}{'yes' if lo > 0 else 'NO':>12s}")
    print("   'NO' = the interval spans zero: the point estimates differ, the metrics do not.")
    print(f"   Uncorrected for the {len(per_metric) - 1} comparisons made, so a lone 'yes' here is "
          "weak evidence; read the column as\n   'almost nothing in this table is separable', "
          "which is the honest summary of a 5-ligand median.")


def specificity(rows, methods):
    print(f"\n== Specificity: correct vs wrong ligand (AUROC, n={len(rows)}) ==")
    print(f"{'method':20s}{'AUROC [95% CI]':>22s}{'scored':>9s}")
    for m in methods:
        v = [(r[m], r["is_binder"]) for r in rows if r.get(m) is not None]
        if not v:
            print(f"{m:20s}{'no scored rows':>22s}{0:>9d}")
            continue
        sc, lab = map(list, zip(*v, strict=True)); sgn = -1 if m in LOWER_BETTER else 1
        # AUROC needs both classes. A partially-scored method (a Rowan run still
        # in flight, say) can cover only binders; that is reported, not skipped,
        # so a method never vanishes from the table without saying why.
        if len(set(lab)) < 2:
            only = "binders" if lab[0] else "non-binders"
            print(f"{m:20s}{f'n/a ({only} only)':>22s}{len(v):>9d}")
            continue
        a = auroc([sgn * x for x in sc], lab)
        lo, up = bootstrap_ci([sgn * x for x in sc], lab, auroc)
        ci = f"{a:.2f} [{lo:.2f},{up:.2f}]" if None not in (a, lo, up) else "n/a"
        print(f"{m:20s}{ci:>22s}{len(v):>9d}")


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
    separability(rows, methods)
    specificity(rows, methods)
    combine(rows)


if __name__ == "__main__":
    main()
