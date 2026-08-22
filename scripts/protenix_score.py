"""Score designs on Protenix-v2, the independent held-out judge.

Its gpde and ranking_score (anchor-validated: a real binder scores gpde ~0.3, scrambled/random
~2.5) discriminate real binding from nonsense far better than a generic interface ipTM. Two modes:
  build  designs JSON -> a Protenix `pred` input JSON (proteinChain + ligand SMILES)
  parse  Protenix output dir -> iptm, ligand-specific ipTM (binder-ligand chain pair),
         gpde (lower better), ranking_score, written back into the designs JSON + a report.
"""

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

# maxAA above this = poly-X: a low-complexity sequence that games held-out ipTM (poly-A folds to a
# stable blob near the ligand and scores ~0.7) without being a designable protein. Disqualifier.
POLY_X = 0.35


def build(rows, out):
    entries = [{"name": f"d{i}", "sequences": [
        {"proteinChain": {"sequence": r["seq"], "count": 1, "id": ["A"]}},
        {"ligand": {"ligand": r["ligand"], "count": 1, "id": ["B"]}}]}
        for i, r in enumerate(rows)]
    Path(out).write_text(json.dumps(entries, indent=2))
    print(f"built {len(entries)} Protenix entries -> {out}", flush=True)


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def _maxaa(seq):
    return max(Counter(seq).values()) / len(seq) if seq else None


def parse(rows, outdir, seed):
    for i, r in enumerate(rows):
        r["maxaa"] = _maxaa(r.get("seq"))          # composition guard vs poly-X gaming
        # protenix writes N diffusion samples per seed; average across them
        hits = sorted(Path(outdir).glob(
            f"d{i}/seed_{seed}/predictions/d{i}_summary_confidence_sample_*.json"))
        if not hits:
            r["protenix_iptm"] = None
            print(f"  MISSING d{i}", flush=True)
            continue
        cs = [json.loads(h.read_text()) for h in hits]
        r["protenix_iptm"] = _mean([c.get("iptm") for c in cs])
        # build() writes exactly two chains, binder (A) first, ligand (B) second, so the binder
        # /ligand cell is [0][1]. Guard the shape so a reordered/extra-chain output fails loudly
        # instead of silently reading the wrong pair.
        r["protenix_lig_iptm"] = _mean(
            [c["chain_pair_iptm"][0][1] for c in cs
             if c.get("chain_pair_iptm") and len(c["chain_pair_iptm"]) == 2
             and len(c["chain_pair_iptm"][0]) == 2])
        r["protenix_gpde"] = _mean([c.get("gpde") for c in cs])
        r["protenix_ranking"] = _mean([c.get("ranking_score") for c in cs])

    def key(r):
        m = r.get("method", "all")                 # benchmark rows have no method; group them together
        return m if m.startswith("anchor") else f"{m}@{r.get('budget', 0)}"
    agg = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r.get("protenix_iptm") is None:
            continue
        for k in ("protenix_iptm", "protenix_lig_iptm", "protenix_gpde", "protenix_ranking", "maxaa"):
            if r.get(k) is not None:
                agg[key(r)][k].append(r[k])
    def fmt(a, k):
        return f"{sum(a[k]) / len(a[k]):.2f}" if a.get(k) else "  -"

    print("\n== Protenix-v2 held-out (mean per group) ==", flush=True)
    print(f"{'group':<20} {'iptm':>6} {'lig_iptm':>9} {'gpde':>6} {'rank':>6} {'maxAA':>6}  n", flush=True)
    for g in sorted(agg):
        a = agg[g]
        n = len(a.get("protenix_iptm", []))
        mean_maxaa = sum(a["maxaa"]) / len(a["maxaa"]) if a.get("maxaa") else None
        flag = "  poly-X (disqualified)" if mean_maxaa is not None and mean_maxaa > POLY_X else ""
        print(f"{g:<20} {fmt(a, 'protenix_iptm'):>6} {fmt(a, 'protenix_lig_iptm'):>9} "
              f"{fmt(a, 'protenix_gpde'):>6} {fmt(a, 'protenix_ranking'):>6} {fmt(a, 'maxaa'):>6}  {n}{flag}",
              flush=True)
    print(f"(maxAA = most-common-residue fraction; > {POLY_X} = poly-X, disqualified: a high held-out "
          "ipTM there is gaming, not binding)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["build", "parse"])
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--input-json", default="/tmp/protenix_in.json")
    ap.add_argument("--outdir", default="/tmp/protenix_out")
    ap.add_argument("--seed", default="101")
    a = ap.parse_args()
    rows = json.loads(Path(a.inp).read_text())
    if a.mode == "build":
        build(rows, a.input_json)
    else:
        parse(rows, a.outdir, a.seed)
        Path(a.inp).write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
