"""Score the benchmark with Rowan's binding-affinity workflow (SQM / GNINA / AEV-PLIG).

Rowan's pose-based scorers need a bound pose, so this reuses the co-folded complexes that
run_predictors.py produces under _work (re-run it first if _work was cleaned). Default pose
source is Protenix-v2 (the benchmark's strongest co-folder; best-ranked of its 5 samples);
--poses boltz uses the Boltz complexes instead. Each complex is split into the binder protein
(apo) and the ligand pose in the same coordinate frame, then submitted in Rowan's mode 2
(apo protein + external pose).

  # 1. poses (GPU): Protenix -> _work/px_out/d{i}/seed_101/predictions/d{i}_sample_*.cif
  PROTENIX_DIR=/path/to/Protenix PROTENIX_PY=/path/to/Protenix/.venv/bin/python \
    CUDA_VISIBLE_DEVICES=0 python scripts/run_predictors.py protenix

  # 2. dry run (no key, no API calls): split poses, print the submission plan
  python scripts/run_rowan.py

  # 3. submit (needs ROWAN_API_KEY): one workflow per (system, method)
  ROWAN_API_KEY=<your-rowan-key> python scripts/run_rowan.py --submit --methods sqm,gnina,aevplig \
      --max-credits 200

Writes predictions/rowan-{sqm,gnina,aevplig}_predictions.csv. Orientation: higher = tighter.
SQM is kcal/mol (negated on write); GNINA/AEV-PLIG are log10(M) (sign confirmed on the first
returned batch against known binders vs negatives before the leaderboard is trusted).
"""

import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
WORK = ROOT / "_work"
POSES = WORK / "rowan_poses"
BOLTZ_OUT = WORK / "boltz_out" / "boltz_results_boltz_yaml" / "predictions"
PX_OUT = WORK / "px_out"
SEED = 101

METHODS = {
    "sqm": ("rowan-sqm", "SinglePointEnergySettings", True),      # kcal/mol, lower = tighter
    "gnina": ("rowan-gnina", "GninaAffinitySettings", None),      # log10(M), sign TBD
    "aevplig": ("rowan-aevplig", "AEVPLIGAffinitySettings", None),
}


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def split_pose(cif, smiles, protein_pdb, ligand_sdf):
    """Co-folded complex -> apo protein PDB (chain A) + ligand SDF (bond orders from SMILES,
    coordinates from the pose). The co-folder emits ligand atoms in the input-SMILES order, so
    bonds are taken from the SMILES template and the pose coordinates transferred by index --
    no distance-based bond perception, which is unreliable on co-folded geometries. An element-
    order guard surfaces any system where that ordering does not hold. Returns None on success,
    an error string on failure (surfaced, never silently skipped)."""
    import gemmi
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]
    prot = gemmi.Structure()
    prot.add_model(gemmi.Model("1"))
    lig_atoms = []
    for ch in model:
        poly = ch.get_polymer()
        if poly and len(poly) > 10:
            nc = gemmi.Chain(ch.name)
            for r in ch:
                nc.add_residue(r)
            prot[0].add_chain(nc)
        else:
            for r in ch:
                lig_atoms.append(r)
    if not len(prot[0]) or not lig_atoms:
        return f"no protein/ligand chain in {cif.name}"
    prot.write_pdb(str(protein_pdb))

    heavy = [(a.element.name, a.pos.x, a.pos.y, a.pos.z)
             for r in lig_atoms for a in r if a.element.name != "H"]
    tmpl = Chem.MolFromSmiles(smiles)
    if tmpl is None:
        return f"bad SMILES template for {cif.name}"
    tmpl_el = [at.GetSymbol() for at in tmpl.GetAtoms()]
    if [e for e, *_ in heavy] != tmpl_el:
        return (f"ligand atom order/count mismatch for {cif.name} "
                f"(pose {len(heavy)} vs template {len(tmpl_el)}) - not index-transferable")
    conf = Chem.Conformer(tmpl.GetNumAtoms())
    for idx, (_, x, y, z) in enumerate(heavy):
        conf.SetAtomPosition(idx, Point3D(x, y, z))
    tmpl.AddConformer(conf, assignId=True)
    molH = Chem.AddHs(tmpl, addCoords=True)
    w = Chem.SDWriter(str(ligand_sdf))
    w.write(molH)
    w.close()
    return None


def _boltz_pose(i):
    cif = BOLTZ_OUT / f"c{i:02d}" / f"c{i:02d}_model_0.cif"
    return cif if cif.exists() else None


def _protenix_pose(i):
    """Top-ranked Protenix sample (5 diffusion samples/system; pick max ranking_score)."""
    d = PX_OUT / f"d{i}" / f"seed_{SEED}" / "predictions"
    best, best_rank = None, -1.0
    for cif in sorted(d.glob(f"d{i}_sample_*.cif")):
        n = cif.stem.split("_sample_")[1]
        try:
            rk = json.loads((d / f"d{i}_summary_confidence_sample_{n}.json").read_text()).get("ranking_score", -1.0)
        except (FileNotFoundError, ValueError):
            rk = -1.0
        if rk > best_rank:
            best, best_rank = cif, rk
    return best


def build_poses(mode):
    POSES.mkdir(parents=True, exist_ok=True)
    pick = _protenix_pose if mode == "protenix" else _boltz_pose
    ok, bad = [], []
    for i, r in enumerate(refs()):
        cif = pick(i)
        if cif is None or not cif.exists():
            bad.append((r["id"], f"missing {mode} pose (index d{i}/c{i:02d})"))
            continue
        pp = POSES / f"{r['id']}_protein.pdb"
        ls = POSES / f"{r['id']}_ligand.sdf"
        err = split_pose(cif, r["smiles"], pp, ls)
        (bad if err else ok).append((r["id"], err) if err else (r["id"], pp, ls))
    return ok, bad


def write_pred(method_col, values):
    """Write method,id,predicted_affinity for all systems. Systems scored this run (in `values`)
    get the new value; any others keep whatever is already on disk, so a targeted --only/--limit
    retry does not blank out the rest of the column."""
    path = PRED / f"{method_col}_predictions.csv"
    existing = {}
    if path.exists():
        existing = {r["id"]: r["predicted_affinity"]
                    for r in csv.DictReader(path.read_text().splitlines())}
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "id", "predicted_affinity"])
        for r in refs():
            if r["id"] in values:
                v = values[r["id"]]
                cell = "" if v is None else v
            else:
                cell = existing.get(r["id"], "")
            w.writerow([method_col, r["id"], cell])
    print(f"WROTE {method_col}_predictions.csv", flush=True)


def submit(methods, ok, max_credits, draft=False):
    import rowan

    settings = {"sqm": rowan.SinglePointEnergySettings, "gnina": rowan.GninaAffinitySettings,
                "aevplig": rowan.AEVPLIGAffinitySettings}
    try:  # folders are UI grouping only; some rowan-python builds mis-define the pydantic model
        folder = rowan.get_folder("denovo-smallmol-benchmark")
    except Exception as e:  # noqa: BLE001 - grouping is optional, never block the science
        print(f"  (folder grouping unavailable: {type(e).__name__}; submitting to account root)")
        folder = None

    # Workflows run concurrently server-side, so submit them all up front and collect afterwards
    # rather than blocking on each result in turn (183 affinity + 61 prep workflows would be hours
    # sequentially, ~an hour batched).
    for key in methods:
        col, _, sqm = METHODS[key]

        # Stage 1 (SQM only): submit every protein-prep so they run in parallel. SQM (MOZYME) needs
        # a protonated closed-shell protein; the co-folded PDB has no H, so geometry-based charge
        # inference yields odd-electron pockets. Protonate at pH 7.4 first.
        preps = {}
        if sqm and not draft:
            for sid, pp, _ in ok:
                try:
                    protein = rowan.upload_protein(sid, str(pp))
                    preps[sid] = rowan.submit_protein_preparation_workflow(
                        protein=protein.uuid, protonation_method="openmm", pH=7.4,
                        name=f"prep {sid}", folder=folder, max_credits=max_credits)
                except Exception as e:  # noqa: BLE001 - surface, keep the batch going
                    preps[sid] = e

        # Stage 2: submit every binding-affinity workflow (all run in parallel).
        jobs = []  # (sid, wf | None, err | None)
        for sid, pp, ls in ok:
            try:
                if sqm and not draft:
                    prep = preps.get(sid)
                    if isinstance(prep, Exception):
                        raise prep
                    target = prep.result().prepared_protein_uuid
                else:
                    target = rowan.upload_protein(sid, str(pp)).uuid
                ligand = next(iter(rowan.load_named_ligands(str(ls)).values()))
                wf = rowan.submit_binding_affinity_workflow(
                    protein=target, ligand_structures=[ligand],
                    binding_affinity_settings=settings[key](),
                    name=f"benchmark {col} {sid}", folder=folder,
                    max_credits=max_credits, is_draft=draft)
                jobs.append((sid, wf, None))
            except Exception as e:  # noqa: BLE001 - surface, keep the batch going
                jobs.append((sid, None, f"{type(e).__name__}: {str(e)[:140]}"))

        if draft:
            for sid, wf, err in jobs:
                print(f"  DRAFT {col} {sid}: {wf.uuid if wf else 'FAILED ' + err}", flush=True)
            continue

        # Stage 3: collect results (submitted together, so this drains as they finish).
        vals = {}
        for sid, wf, err in jobs:
            if wf is None:
                vals[sid] = None
                print(f"  {col} {sid}: FAILED ({err})", flush=True)
                continue
            try:
                score = wf.result().scores[0]
                raw = None if score is None else score.binding_affinity
                vals[sid] = None if raw is None else (-raw if sqm else raw)
                print(f"  {col} {sid}: {raw}", flush=True)
            except Exception as e:  # noqa: BLE001 - surface, keep the batch going
                vals[sid] = None
                print(f"  {col} {sid}: FAILED ({type(e).__name__}: {str(e)[:140]})", flush=True)
        write_pred(col, vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poses", choices=["protenix", "boltz"], default="protenix",
                    help="co-folder whose poses to score (protenix: best-ranked of 5 samples)")
    ap.add_argument("--methods", default="sqm,gnina,aevplig")
    ap.add_argument("--submit", action="store_true", help="actually call the Rowan API (needs ROWAN_API_KEY)")
    ap.add_argument("--max-credits", type=int, default=None, help="hard credit cap per workflow")
    ap.add_argument("--limit", type=int, default=None, help="use only the first N systems (smoke test)")
    ap.add_argument("--only", default=None, help="comma-separated system ids to run (targeted retry)")
    ap.add_argument("--draft", action="store_true",
                    help="with --submit: create is_draft workflows (validate the full path, no execution, no credits)")
    args = ap.parse_args()
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    ok, bad = build_poses(args.poses)
    print(f"\n{args.poses} poses: {len(ok)} ready, {len(bad)} failed")
    for b in bad:
        print(f"  FAILED {b[0]}: {b[1]}")
    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        ok = [r for r in ok if r[0] in want]
        print(f"--only: {len(ok)} system(s) {[r[0] for r in ok]}")
    if args.limit:
        ok = ok[:args.limit]
        print(f"--limit {args.limit}: using first {len(ok)} system(s) only")

    if not args.submit:
        n = len(ok) * len(methods)
        print(f"\nDRY RUN — no API calls. Would submit {n} workflows "
              f"({len(ok)} systems x {len(methods)} methods: {','.join(methods)}).")
        print("SQM dominates cost (QM opt + single-point); GNINA/AEV-PLIG are cheap ML rescoring.")
        print("Set --max-credits to cap each workflow, then re-run with --submit and ROWAN_API_KEY.")
        return
    if not os.environ.get("ROWAN_API_KEY"):
        raise SystemExit("--submit needs ROWAN_API_KEY in the environment")
    submit(methods, ok, args.max_credits, draft=args.draft)


if __name__ == "__main__":
    main()
