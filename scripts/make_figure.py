"""Two-panel benchmark figure. Left: how well each method ranks measured pKd among the 37 binders
(Spearman; dotted = cLogP baseline). Right: correct vs wrong ligand across all 47 (AUROC; dotted =
chance). The leave-one-out combination is shown as its own bar. Bars are 95% bootstrap CIs.

Shows one metric per model plus the baselines for readability; score.py reports all 11. Stats: score.py.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import auroc, bootstrap_ci, combine_loo, load, spearman

ROOT = Path(__file__).resolve().parents[1]
# (method, label, category) shown in the figure
SHOW = [("boltz-2", "Boltz-2", "aff"), ("nesso-1", "Nesso-1", "aff"),
        ("boltz-2-iptm", "Boltz-2 ipTM", "prox"), ("protenix-iptm", "Protenix ipTM", "prox"),
        ("clogp", "cLogP", "base"), ("molecular-weight", "MW", "base")]
COL = {"aff": "#35617a", "prox": "#c08a4a", "base": "#cbc3b8", "combo": "#8a5a72"}
REF, WHISK = "#8f8880", "#9b938a"
LABEL = {"aff": "dedicated affinity head", "prox": "structural proxy",
         "base": "trivial baseline", "combo": "combination (LOO)"}


def bars(ax, data, ref, ref_label, xlabel, xlim):
    data.sort(key=lambda t: t[1])
    n = len(data)
    for y, (lab, v, lo, hi, cat) in enumerate(data):
        ax.barh(y, v, color=COL[cat], height=0.62, zorder=3)
        ax.plot([lo, hi], [y, y], color=WHISK, lw=0.9, alpha=0.9, zorder=2)
    ax.axvline(ref, ls=":", lw=1.3, color=REF, zorder=1)
    ax.text(ref, n - 0.35, f" {ref_label}", color=REF, fontsize=9, ha="left", va="center")
    ax.set_yticks(range(n)); ax.set_yticklabels([d[0] for d in data], fontsize=10)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(xlabel, fontsize=10); ax.set_xlim(*xlim)
    ax.tick_params(axis="x", labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    rows, _ = load()
    binders = [r for r in rows if r["is_binder"] == 1 and r["pKd"] is not None]
    reg, spec = [], []
    for m, lab, cat in SHOW:
        bx, by = zip(*[(r[m], r["pKd"]) for r in binders if r.get(m) is not None])
        reg.append((lab, spearman(list(bx), list(by)), *bootstrap_ci(list(bx), list(by), spearman), cat))
        sx, sl = zip(*[(r[m], r["is_binder"]) for r in rows if r.get(m) is not None])
        spec.append((lab, auroc(list(sx), list(sl)), *bootstrap_ci(list(sx), list(sl), auroc), cat))
    (rp, ry), (sp, sy) = combine_loo(rows)
    reg.append(("combination", spearman(rp, ry), *bootstrap_ci(rp, ry, spearman), "combo"))
    spec.append(("combination", auroc(sp, sy), *bootstrap_ci(sp, sy, auroc), "combo"))

    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 5.2))
    bars(a, reg, 0.40, "cLogP", "Spearman rho vs measured pKd", (-0.1, 0.85))
    bars(b, spec, 0.5, "chance", "AUROC: correct vs wrong ligand", (0.0, 1.0))
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL[c]) for c in ("aff", "prox", "base", "combo")]
    a.legend(handles, [LABEL[c] for c in ("aff", "prox", "base", "combo")], fontsize=8.5,
             loc="lower right", frameon=False)
    fig.suptitle("Affinity predictors on de-novo protein-small-molecule binders", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = ROOT / "figures" / "benchmark.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
