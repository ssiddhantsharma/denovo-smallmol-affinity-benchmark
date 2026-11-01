"""ipDAE for protein-small-molecule interfaces, from the Boltz-2 folds (no GPU, no network).

ipDAE (BindCraft2's `i_pDAE`) re-weights the interface PAE by a TM-score-like term whose distance
scale d0 grows with how many partners an anchor has, then takes the max over anchors and over both
directions. On protein-protein interfaces it is one of the stronger confidence scores. This is the
protein-ligand port, and it differs from the protein-protein one in two ways that matter:

1. CONTACT SET. The protein-protein version calls two residues in contact at CA-CA <= 8 A. A ligand
   has no CA, so contact here is residue-heavy-atom to ligand-heavy-atom <= 4.5 A -- the same
   definition `run_interface.py` already uses for ncontacts/burial, so the two agree on what the
   interface is. A CA-to-ligand-atom cutoff would be neither one thing nor the other.

2. d0 SATURATES, AND THAT IS THE INTERESTING PART. d0 = max(1.24*(max(n,19)-15)^(1/3) - 1.8, 1.0)
   over an anchor's partner count n. Solving d0 > 1 needs n >= 27 partners. A ligand atom contacting
   27 protein residues does not happen; a residue contacting 27 ligand heavy atoms is rare. So d0
   pins to 1.0 almost everywhere and the score collapses to max over anchors of mean 1/(1+PAE^2).
   That is still not a monotone transform of mean interface PAE -- the nonlinearity is applied
   before the average, and the max picks a single anchor -- so it can rank differently from
   `boltz-2-pae-mean`. But it is the same matrix read a slightly different way, so expect it to land
   near the existing iPAE columns rather than to add an independent axis. `--report` prints the
   fraction of anchors with d0 > 1 so this is measured rather than assumed.

Pairs Boltz-2's PAE with Boltz-2's OWN structure. The Protenix pose that run_interface.py reads is
a different model's coordinates and a different token order; crossing them would silently score one
model's confidence against another's geometry.

Writes predictions/boltz-2-ipdae_predictions.csv. Higher = more confident interface.
"""

import argparse
import csv
import glob
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
BOLTZ_OUT = ROOT / "_work" / "boltz_out" / "boltz_results_boltz_yaml" / "predictions"
CUT = 4.5


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def d0_from_partner_count(partner_count):
    """BindCraft2 `tm_score_distance_scale`, reproduced exactly."""
    return np.maximum(1.24 * (np.maximum(partner_count, 19.0) - 15.0) ** (1.0 / 3.0) - 1.8, 1.0)


def anchored_interface_tm(pae, contact):
    """Per-anchor score: mean TM-like term over that anchor's CONTACT partners only."""
    partner_counts = contact.sum(axis=-1).astype(float)
    d0 = d0_from_partner_count(partner_counts)[:, None]
    tm = 1.0 / (1.0 + (pae / d0) ** 2)
    return np.where(contact, tm, 0.0).sum(axis=-1) / np.maximum(partner_counts, 1.0)


def ipdae_from_block(pae_block, contact):
    """max over both directions of max over anchors, matching BindCraft2.

    `pae_block` is (n_protein_residues, n_ligand_atoms), already sliced out of the full matrix.
    Returns (score, frac_d0_above_1) or (None, None) when nothing is in contact.
    """
    pae_block = np.asarray(pae_block, dtype=float)
    contact = np.asarray(contact, dtype=bool)
    if pae_block.shape != contact.shape:
        raise ValueError(f"pae block {pae_block.shape} != contact {contact.shape}")
    if not contact.any():
        return None, None
    ab = anchored_interface_tm(pae_block, contact)
    ba = anchored_interface_tm(pae_block.T, contact.T)
    counts = np.concatenate([contact.sum(1), contact.sum(0)])
    counts = counts[counts > 0]
    frac = float((d0_from_partner_count(counts) > 1.0).mean())
    return float(max(ab.max(), ba.max())), frac


def split_chains(cif):
    """Protein residues (as gemmi residues) and ligand heavy-atom coords, Boltz-2 output.

    Same polymer-length rule as run_interface.py so both agree on which chain is the ligand.
    """
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    prot, lig = [], []
    for ch in st[0]:
        poly = ch.get_polymer()
        (prot if poly and len(poly) > 10 else lig).extend(ch)
    lig_xyz = np.array([[a.pos.x, a.pos.y, a.pos.z]
                        for r in lig for a in r if a.element.name != "H"])
    return prot, lig_xyz


def contact_matrix(prot_residues, lig_xyz, cut=CUT):
    """(n_residues, n_ligand_atoms) bool: any heavy atom of the residue within `cut` of that atom."""
    out = np.zeros((len(prot_residues), len(lig_xyz)), dtype=bool)
    for i, res in enumerate(prot_residues):
        xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"])
        if not len(xyz):
            continue
        d = np.sqrt(((xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1))
        out[i] = (d <= cut).any(axis=0)
    return out


def score_one(i, seq_len, report=False):
    hits = glob.glob(str(BOLTZ_OUT / f"c{i:02d}" / "pae*.npz"))
    cifs = sorted(glob.glob(str(BOLTZ_OUT / f"c{i:02d}" / "*.cif")))
    if not hits or not cifs:
        return None, None
    pae = np.load(hits[0])["pae"]
    prot, lig_xyz = split_chains(cifs[0])
    if not len(lig_xyz) or not len(prot):
        return None, None

    # The token layout (protein residues 0..L, then ligand atoms) is an assumption about Boltz-2's
    # output, and getting it wrong shifts every index silently. Check it instead of trusting it.
    n_tok = pae.shape[0]
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"c{i:02d}: PAE is {pae.shape}, expected square")
    if len(prot) != seq_len:
        raise ValueError(f"c{i:02d}: structure has {len(prot)} protein residues, "
                         f"sequence has {seq_len} -- token mapping would be wrong")
    if n_tok != seq_len + len(lig_xyz):
        raise ValueError(f"c{i:02d}: PAE is {n_tok} tokens but protein {seq_len} + ligand "
                         f"{len(lig_xyz)} = {seq_len + len(lig_xyz)} -- ligand may not be "
                         f"atom-tokenized as assumed; refusing to guess the offset")

    contact = contact_matrix(prot, lig_xyz)
    block = pae[:seq_len, seq_len:]
    score, frac = ipdae_from_block(block, contact)
    if report and score is not None:
        print(f"  c{i:02d}: ipdae={score:.4f}  contacts={int(contact.sum())}  "
              f"anchors with d0>1: {frac:.1%}")
    return score, frac


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="store_true", help="per-system diagnostics incl. d0 saturation")
    args = ap.parse_args()

    if not BOLTZ_OUT.exists():
        raise SystemExit(f"{BOLTZ_OUT} not found. ipDAE needs the Boltz-2 fold outputs "
                         f"(_work/ is gitignored); re-run scripts/run_predictors.py boltz first.")

    vals, fracs = {}, []
    for i, r in enumerate(refs()):
        score, frac = score_one(i, len(r["sequence"]), report=args.report)
        if score is None:
            print(f"  MISSING Boltz-2 pae/cif for {r['id']} (c{i:02d})", flush=True)
            continue
        vals[r["id"]] = score
        fracs.append(frac)

    with open(PRED / "boltz-2-ipdae_predictions.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "id", "predicted_affinity"])
        for r in refs():
            v = vals.get(r["id"])
            w.writerow(["boltz-2-ipdae", r["id"], "" if v is None else v])
    print(f"WROTE boltz-2-ipdae_predictions.csv ({len(vals)}/{len(refs())})", flush=True)
    if fracs:
        m = float(np.mean(fracs))
        print(f"anchors with d0 > 1: {m:.1%} on average.", flush=True)
        if m < 0.01:
            print("  d0 pinned at 1.0 throughout: on this interface size ipDAE is "
                  "max-anchor mean 1/(1+PAE^2), not a partner-count-scaled score.", flush=True)


if __name__ == "__main__":
    main()
