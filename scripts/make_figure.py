"""Two-panel benchmark figure: rank-affinity (Spearman, binders) and specificity (AUROC, all rows).

Left: how well each method ranks measured pKd among the 37 binders (dashed = cLogP baseline, star =
leave-one-out combination). Right: cognate vs wrong ligand (dashed = chance). Data + stats: score.py.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import LOWER_BETTER, auroc, bootstrap_ci, load, spearman

ROOT = Path(__file__).resolve().parents[1]
CAT = {"boltz-2": "aff", "boltz-2-pbind": "aff", "nesso-1": "aff", "nesso-1-pbind": "aff",
       "boltz-2-iptm": "prox", "protenix-iptm": "prox", "protenix-ligiptm": "prox",
       "protenix-gpde": "prox", "protenix-ranking": "prox", "clogp": "base", "molecular-weight": "base"}
COL = {"aff": "#0072B2", "prox": "#E69F00", "base": "#999999"}
LABEL = {"aff": "dedicated affinity head", "prox": "structural proxy", "base": "trivial baseline"}
COMBO_RHO, COMBO_AUROC = 0.59, 0.83   # leave-one-out CV, from score.py


def bars(ax, data, ref, ref_label, xlabel, combo, xlim):
    data.sort(key=lambda t: t[1])
    n = len(data); ys = range(n)
    for y, (k, v, lo, hi) in zip(ys, data):
        ax.barh(y, v, color=COL[CAT[k]], height=0.66, zorder=2)
        ax.plot([lo, hi], [y, y], color="#333333", lw=1.0, zorder=3)
    ax.axvline(ref, ls="--", lw=1.2, color="#D55E00", zorder=1)
    ax.text(ref, n - 0.4, f" {ref_label}", color="#D55E00", fontsize=8.5, ha="left", va="center")
    ax.scatter([combo], [n + 0.15], marker="*", s=180, color="#009E73", zorder=4, clip_on=False)
    ax.text(combo, n + 0.15, "  combination (LOO)", color="#009E73", fontsize=8.5, va="center")
    ax.set_yticks(list(ys)); ax.set_yticklabels([k for k, *_ in data], fontsize=9.5)
    ax.set_ylim(-0.7, n + 0.7)
    ax.set_xlabel(xlabel, fontsize=9.5); ax.set_xlim(*xlim)
    ax.tick_params(axis="x", labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    rows, methods = load()
    binders = [r for r in rows if r["is_binder"] == 1 and r["pKd"] is not None]
    reg, spec = [], []
    for m in methods:
        sgn = -1 if m in LOWER_BETTER else 1
        bx, by = zip(*[(sgn * r[m], r["pKd"]) for r in binders if r.get(m) is not None])
        reg.append((m, sgn * spearman([sgn * x for x in bx], list(by)),
                    *bootstrap_ci(list(bx), list(by), spearman)))
        sx, sl = zip(*[(sgn * r[m], r["is_binder"]) for r in rows if r.get(m) is not None])
        spec.append((m, auroc(list(sx), list(sl)), *bootstrap_ci(list(sx), list(sl), auroc)))

    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 6.6))
    bars(a, reg, 0.40, "cLogP", "Spearman rho vs measured pKd  (37 binders)", COMBO_RHO, (-0.2, 0.85))
    bars(b, spec, 0.5, "chance", "AUROC: cognate vs wrong ligand  (47)", COMBO_AUROC, (0.0, 1.05))
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL[c]) for c in ("aff", "prox", "base")]
    a.legend(handles, [LABEL[c] for c in ("aff", "prox", "base")], fontsize=8.5,
             loc="lower right", frameon=False)
    fig.suptitle("Affinity predictors on de-novo protein-small-molecule binders", fontsize=12, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = ROOT / "figures" / "benchmark.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
