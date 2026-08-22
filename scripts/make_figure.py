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
    ys = range(len(data))
    for y, (k, v, lo, hi) in zip(ys, data):
        ax.barh(y, v, color=COL[CAT[k]], height=0.62, zorder=2)
        ax.plot([lo, hi], [y, y], color="#333333", lw=1.1, zorder=3)
    ax.axvline(ref, ls="--", lw=1.2, color="#D55E00", zorder=1)
    ax.text(ref, len(data) - 0.3, f" {ref_label}", color="#D55E00", fontsize=7.5, va="top")
    ax.scatter([combo], [len(data) - 0.5], marker="*", s=150, color="#009E73", zorder=4)
    ax.text(combo, len(data) - 0.5, "  combination (LOO)", color="#009E73", fontsize=7.5, va="center")
    ax.set_yticks(list(ys)); ax.set_yticklabels([k for k, *_ in data], fontsize=8)
    ax.set_xlabel(xlabel); ax.set_xlim(*xlim)
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

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 4.4))
    bars(a, reg, 0.40, "cLogP", "Spearman rho vs measured pKd  (37 binders)", COMBO_RHO, (-0.2, 0.8))
    bars(b, spec, 0.5, "chance", "AUROC: cognate vs wrong ligand  (47)", COMBO_AUROC, (0.0, 1.0))
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL[c]) for c in ("aff", "prox", "base")]
    a.legend(handles, [LABEL[c] for c in ("aff", "prox", "base")], fontsize=7.5,
             loc="lower right", frameon=False)
    fig.suptitle("Affinity predictors on de-novo protein-small-molecule binders", fontsize=11)
    fig.tight_layout()
    out = ROOT / "figures" / "benchmark.png"
    fig.savefig(out, dpi=140, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
