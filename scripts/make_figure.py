"""Benchmark figure: figures/benchmark.png, every metric, two panels.

Left panel ranks pKd within-ligand (median Spearman; dashed line = no correlation). Right panel tells
the correct ligand from the wrong one (AUROC over all pairs; dashed line = chance). Bars are 95%
bootstrap CIs. The leave-one-out combination is its own bar. Stats: score.py.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import (
    LOWER_BETTER,
    auroc,
    bootstrap_ci,
    combine_loo,
    load,
    per_ligand_rho,
    within_ligand_rho,
)

ROOT = Path(__file__).resolve().parents[1]
CAT = {"boltz-2": "aff", "boltz-2-pbind": "aff", "nesso-1": "aff", "nesso-1-pbind": "aff",
       "boltz-2-iptm": "prox", "boltz-2-pae-min": "prox", "boltz-2-pae-mean": "prox",
       "protenix-iptm": "prox", "protenix-ligiptm": "prox", "protenix-gpde": "prox",
       "protenix-ranking": "prox", "protenix-pae-min": "prox", "protenix-pae-mean": "prox",
       "protenix-ptm": "prox", "protenix-plddt-binder": "prox", "protenix-plddt-ligand": "prox",
       "clogp": "base", "molecular-weight": "base", "combination": "combo",
       "dtsfm-cosine": "seq",
       "rowan-sqm": "phys", "rowan-gnina": "dock", "rowan-aevplig": "dock",
       "interface-ncontacts": "iface", "interface-ligburial": "iface"}
COL = {"aff": "#35617a", "prox": "#c08a4a", "base": "#cbc3b8", "combo": "#8a5a72", "seq": "#7B4FA0",
       "phys": "#2e8b6b", "dock": "#a0522d", "iface": "#5a7d8c"}
REF, WHISK = "#8f8880", "#9b938a"
LABEL = {"aff": "dedicated affinity head", "prox": "structural proxy",
         "base": "trivial baseline", "combo": "combination (LOO)", "seq": "sequence-native FM",
         "phys": "physics (SQM)", "dock": "docking / ML rescore", "iface": "interface geometry"}
PRETTY = {"boltz-2": "Boltz-2", "nesso-1": "Nesso-1", "boltz-2-iptm": "Boltz-2 ipTM",
          "protenix-iptm": "Protenix ipTM", "clogp": "cLogP", "molecular-weight": "MW",
          "dtsfm-cosine": "dtSFM",
          "rowan-sqm": "SQM", "rowan-gnina": "GNINA", "rowan-aevplig": "AEV-PLIG",
          "interface-ncontacts": "iface contacts", "interface-ligburial": "iface burial",
          "protenix-ptm": "Protenix pTM", "protenix-plddt-binder": "Protenix pLDDT bind",
          "protenix-plddt-ligand": "Protenix pLDDT lig",
          "boltz-2-pbind": "Boltz-2 P(bind)", "nesso-1-pbind": "Nesso-1 P(bind)",
          "boltz-2-pae-min": "Boltz-2 iPAE min", "boltz-2-pae-mean": "Boltz-2 iPAE mean",
          "protenix-pae-min": "Protenix iPAE min", "protenix-pae-mean": "Protenix iPAE mean",
          "protenix-ligiptm": "Protenix lig-ipTM", "protenix-gpde": "Protenix gPDE",
          "protenix-ranking": "Protenix ranking"}


def not_separable(rows, methods, kmin=3, n_boot=4000):
    """How many metrics the top one is NOT distinguishable from (paired bootstrap over the shared
    ligands). Mirrors score.py's separability table so the figure cannot drift from the CLI."""
    import random
    import statistics
    per_metric = {m: per_ligand_rho(rows, m, kmin)[0] for m in methods}
    per_metric = {m: p for m, p in per_metric.items() if len(p) >= 2}
    best = max(per_metric, key=lambda m: statistics.median(per_metric[m].values()))
    bper = per_metric[best]
    n_ns = 0
    for m, per in per_metric.items():
        if m == best:
            continue
        shared = sorted(set(per) & set(bper))
        if len(shared) < 2:
            continue
        bv = [bper[k] for k in shared]; mv = [per[k] for k in shared]
        rng = random.Random(0); k = len(shared)
        d = sorted(statistics.median([bv[i] for i in idx]) - statistics.median([mv[i] for i in idx])
                   for idx in ([rng.randrange(k) for _ in range(k)] for _ in range(n_boot)))
        if d[int(0.025 * n_boot)] <= 0:
            n_ns += 1
    return n_ns, len(per_metric) - 1, len(bper)


def _wlig_ci(per, n):
    import random
    import statistics
    rng = random.Random(0)
    v = list(per.values()) if isinstance(per, dict) else list(per)
    boots = sorted(statistics.median([v[rng.randrange(n)] for _ in range(n)]) for _ in range(2000))
    return statistics.median(v), boots[int(0.025 * len(boots))], boots[int(0.975 * len(boots))]


def stats(rows, methods):
    # Pooled Spearman is omitted on purpose: across the set's many ligands it is confounded by
    # ligand composition (cLogP "wins" as an artifact), so only the composition-controlled axes are
    # shown. score.py still prints the pooled numbers for transparency.
    wlig, spec = [], []
    for m in methods:
        sgn = -1 if m in LOWER_BETTER else 1
        med, per, ng, npairs = within_ligand_rho(rows, m)
        if med is not None:
            wlig.append((PRETTY.get(m, m), *_wlig_ci(per, ng), CAT[m], f"{ng} lig / {npairs}"))
        sx, sl = zip(*[(sgn * r[m], r["is_binder"]) for r in rows if r.get(m) is not None], strict=True)
        spec.append((PRETTY.get(m, m), auroc(list(sx), list(sl)),
                     *bootstrap_ci(list(sx), list(sl), auroc), CAT[m], f"{len(sx)}"))
    (_, _), (sp, sy) = combine_loo(rows)
    spec.append(("combination", auroc(sp, sy), *bootstrap_ci(sp, sy, auroc), "combo", f"{len(sy)}"))
    return wlig, spec, not_separable(rows, methods)


def bars(ax, data, ref, ref_label, xlabel, xlim, n_header):
    data.sort(key=lambda t: t[1])
    n = len(data)
    ax.grid(axis="x", color="#ece8e2", lw=0.6, zorder=0)
    for y, (_lab, v, lo, hi, cat, nlab) in enumerate(data):
        ax.barh(y, v, color=COL[cat], height=0.64, zorder=3)
        ax.plot([lo, hi], [y, y], color=WHISK, lw=0.9, alpha=0.9, zorder=4)
        # the n behind each bar, in the plot's own x units so it never collides with the whisker
        ax.annotate(nlab, xy=(xlim[1], y), xytext=(4, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=7.5, color="#8a8178", annotation_clip=False)
    ax.annotate(n_header, xy=(xlim[1], n - 0.4), xytext=(4, 0), textcoords="offset points",
                ha="left", va="center", fontsize=7.5, color="#8a8178", fontstyle="italic",
                annotation_clip=False)
    ax.axvline(ref, ls="--", lw=1.1, color=REF, zorder=2)
    ax.annotate(ref_label, xy=(ref, n - 0.5), xytext=(0, 7), textcoords="offset points",
                ha="center", va="bottom", fontsize=8.5, color=REF, fontstyle="italic")
    ax.set_yticks(range(n)); ax.set_yticklabels([d[0] for d in data], fontsize=8.5)
    ax.set_ylim(-0.7, n + 0.3)
    ax.set_xlabel(xlabel, fontsize=10); ax.set_xlim(*xlim)
    ax.tick_params(axis="x", labelsize=9)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def render(wlig, spec, ns, out, height):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, height))
    # -1.05 not -0.85: several CIs reach rho = -1, and a clipped whisker understates the spread.
    bars(a, wlig, 0.0, "no correlation (rho = 0)", "within-ligand Spearman rho vs pKd",
         (-1.05, 1.0), "n  (ligands / pairs)")
    bars(b, spec, 0.5, "chance (AUROC = 0.5)", "AUROC: correct vs wrong ligand", (0.0, 1.0), "n")
    a.set_title("rank KD within a ligand  (composition controlled)", fontsize=10.5, pad=34)
    n_ns, n_cmp, n_lig = ns
    a.annotate(f"{n_ns} of {n_cmp} metrics are not separable from the top one\n"
               f"(paired bootstrap over the {n_lig} ligands with >=3 binders) -- read the order, "
               "not the ranking",
               xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 6), textcoords="offset points",
               ha="center", va="bottom", fontsize=8, color="#8a8178", fontstyle="italic")
    b.set_title("tell the right ligand from the wrong one", fontsize=10.5, pad=16)
    present = {d[4] for d in spec}
    cats = [c for c in ("aff", "prox", "phys", "dock", "iface", "seq", "base", "combo") if c in present]
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL[c]) for c in cats]
    fig.legend(handles, [LABEL[c] for c in cats], fontsize=8.5, loc="lower center",
               ncol=len(cats), frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Affinity predictors on de-novo protein-small-molecule binders", fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


def main():
    rows, methods = load()
    render(*stats(rows, methods), ROOT / "figures/benchmark.png", 8.8)


if __name__ == "__main__":
    main()
