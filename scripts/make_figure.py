"""Benchmark figures. Two per run:
  benchmark.png       select metrics (README face): affinity heads vs a proxy vs baselines + combination
  benchmark_full.png  every metric score.py computes

Left panel ranks measured pKd (Spearman among the binders; dotted = cLogP). Right panel tells the
correct ligand from the wrong one (AUROC over all pairs; dotted = chance). Bars are 95% bootstrap CIs.
The leave-one-out combination is its own bar. Stats: score.py.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from score import LOWER_BETTER, auroc, bootstrap_ci, combine_loo, load, spearman

ROOT = Path(__file__).resolve().parents[1]
CAT = {"boltz-2": "aff", "boltz-2-pbind": "aff", "nesso-1": "aff", "nesso-1-pbind": "aff",
       "boltz-2-iptm": "prox", "boltz-2-pae-min": "prox", "boltz-2-pae-mean": "prox",
       "protenix-iptm": "prox", "protenix-ligiptm": "prox", "protenix-gpde": "prox",
       "protenix-ranking": "prox", "protenix-pae-min": "prox", "protenix-pae-mean": "prox",
       "clogp": "base", "molecular-weight": "base", "combination": "combo",
       "dtsfm-cosine": "seq"}
COL = {"aff": "#35617a", "prox": "#c08a4a", "base": "#cbc3b8", "combo": "#8a5a72", "seq": "#7B4FA0"}
REF, WHISK = "#8f8880", "#9b938a"
LABEL = {"aff": "dedicated affinity head", "prox": "structural proxy",
         "base": "trivial baseline", "combo": "combination (LOO)", "seq": "sequence-native FM"}
SELECT = ["boltz-2", "nesso-1", "boltz-2-iptm", "protenix-iptm", "clogp", "molecular-weight"]
PRETTY = {"boltz-2": "Boltz-2", "nesso-1": "Nesso-1", "boltz-2-iptm": "Boltz-2 ipTM",
          "protenix-iptm": "Protenix ipTM", "clogp": "cLogP", "molecular-weight": "MW",
          "dtsfm-cosine": "dtSFM"}


def stats(rows, methods):
    binders = [r for r in rows if r["is_binder"] == 1 and r["pKd"] is not None]
    reg, spec = [], []
    for m in methods:
        sgn = -1 if m in LOWER_BETTER else 1
        bx, by = zip(*[(sgn * r[m], r["pKd"]) for r in binders if r.get(m) is not None], strict=True)
        reg.append((PRETTY.get(m, m), spearman(list(bx), list(by)),
                    *bootstrap_ci(list(bx), list(by), spearman), CAT[m]))
        sx, sl = zip(*[(sgn * r[m], r["is_binder"]) for r in rows if r.get(m) is not None], strict=True)
        spec.append((PRETTY.get(m, m), auroc(list(sx), list(sl)),
                     *bootstrap_ci(list(sx), list(sl), auroc), CAT[m]))
    (rp, ry), (sp, sy) = combine_loo(rows)
    reg.append(("combination", spearman(rp, ry), *bootstrap_ci(rp, ry, spearman), "combo"))
    spec.append(("combination", auroc(sp, sy), *bootstrap_ci(sp, sy, auroc), "combo"))
    return reg, spec


def bars(ax, data, ref, ref_label, xlabel, xlim):
    data.sort(key=lambda t: t[1])
    n = len(data)
    for y, (_lab, v, lo, hi, cat) in enumerate(data):
        ax.barh(y, v, color=COL[cat], height=0.62, zorder=3)
        ax.plot([lo, hi], [y, y], color=WHISK, lw=0.9, alpha=0.9, zorder=2)
    ax.axvline(ref, ls=":", lw=1.3, color=REF, zorder=1)
    ax.text(ref, n - 0.35, f" {ref_label}", color=REF, fontsize=8.5, ha="left", va="center")
    ax.set_yticks(range(n)); ax.set_yticklabels([d[0] for d in data], fontsize=9)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel(xlabel, fontsize=9.5); ax.set_xlim(*xlim)
    ax.tick_params(axis="x", labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def render(reg, spec, out, height):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, height))
    clogp_rho = next(v for lab, v, *_ in reg if lab == "cLogP")
    bars(a, reg, clogp_rho, "cLogP", "Spearman rho vs measured pKd", (-0.2, 0.85))
    bars(b, spec, 0.5, "chance", "AUROC: correct vs wrong ligand", (0.0, 1.0))
    present = {d[4] for d in reg}
    cats = [c for c in ("aff", "prox", "seq", "base", "combo") if c in present]
    handles = [plt.Rectangle((0, 0), 1, 1, color=COL[c]) for c in cats]
    fig.legend(handles, [LABEL[c] for c in cats], fontsize=8.5, loc="lower center",
               ncol=len(cats), frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Affinity predictors on de-novo protein-small-molecule binders", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


def main():
    rows, methods = load()
    render(*stats(rows, SELECT), ROOT / "figures/benchmark.png", 5.2)
    render(*stats(rows, methods), ROOT / "figures/benchmark_full.png", 8.4)


if __name__ == "__main__":
    main()
