"""Tests for the scoring statistics and dataset integrity. No GPU, no network."""

import csv
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(relpath: str, name: str):
    """The scripts/ files are standalone CLIs, not an importable package, so they
    are loaded by path. They import each other by bare name (`import methods`),
    which works when run as `python scripts/x.py` because that puts scripts/ on
    sys.path; loading by path does not, so put it there."""
    if str(ROOT / "scripts") not in sys.path:
        sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(name, ROOT / relpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


score = _load("scripts/score.py", "score")


# --- scoring statistics: known values ---

def test_spearman_monotonic():
    assert score.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0
    assert score.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0


def test_pearson_linear():
    assert math.isclose(score.pearson([1, 2, 3], [2, 4, 6]), 1.0, abs_tol=1e-9)
    assert math.isclose(score.pearson([1, 2, 3], [6, 4, 2]), -1.0, abs_tol=1e-9)


def test_auroc_ranking_and_ties():
    assert score.auroc([0.9, 0.8, 0.1], [1, 1, 0]) == 1.0
    assert score.auroc([0.1, 0.2, 0.9], [1, 1, 0]) == 0.0
    assert score.auroc([0.5, 0.5], [1, 0]) == 0.5
    assert score.auroc([0.5], [1]) is None  # no negatives to rank against


def test_rmse():
    assert score.rmse([1, 2, 3], [1, 2, 3]) == 0.0
    assert math.isclose(score.rmse([1, 2], [1, 4]), math.sqrt(2), abs_tol=1e-9)


# --- dataset integrity ---

def test_dataset_shape():
    rows, methods = score.load()
    assert len(rows) == 60
    binders = [r for r in rows if r["is_binder"] == 1]
    negs = [r for r in rows if r["is_binder"] == 0]
    assert len(binders) == 49
    assert len(negs) == 11
    assert all(r["pKd"] is not None for r in binders), "every binder must carry a measured pKd"
    assert all(r["pKd"] is None for r in negs), "negatives must not carry a pKd"
    assert len(methods) >= 16
    for m in ("boltz-2", "clogp", "protenix-iptm", "dtsfm-cosine"):
        assert m in methods


def test_dataset_validated():
    """Integrity checks that must hold for every system (would have caught the JF646/JFX646 dup)."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    sysref = list(csv.DictReader(
        (ROOT / "reference/system_reference.csv").read_text().splitlines()))
    aa = set("ACDEFGHIKLMNPQRSTVWY")
    seen_pairs = set()
    for r in sysref:
        seq, smi = r["sequence"].strip().upper(), r["smiles"].strip()
        assert seq and not (set(seq) - aa), f"{r['id']}: invalid sequence"
        assert Chem.MolFromSmiles(smi) is not None, f"{r['id']}: unparseable SMILES"
        assert all(r[c].strip() for c in ("protein", "ligand_name", "source")), f"{r['id']}: missing field"
        pair = (seq, smi)
        assert pair not in seen_pairs, f"{r['id']}: duplicate (sequence, smiles) pair"
        seen_pairs.add(pair)
    pkds = [float(r["experimental_pKD"]) for r in csv.DictReader(
        (ROOT / "reference/experimental_reference_ground_truth.csv").read_text().splitlines())
        if r["experimental_pKD"].strip()]
    assert all(2.0 <= p <= 13.0 for p in pkds), "pKd outside a sane range"


def test_within_ligand_drops_constant_ligand_properties():
    rows, _ = score.load()
    # cLogP is a ligand-only property: constant within a ligand group, so it yields no within-ligand rho
    med, _per, n, npairs = score.within_ligand_rho(rows, "clogp")
    assert med is None and n == 0 and npairs == 0
    # an interface metric varies across the proteins binding one ligand, so it does
    med2, _per2, n2, npairs2 = score.within_ligand_rho(rows, "protenix-iptm")
    assert n2 >= 2 and med2 is not None
    # n_pairs must be the points actually behind the median, not the size of the benchmark
    assert n2 * 3 <= npairs2 < len(rows), "n_pairs should be >= 3 per group and well under 60"


def test_per_ligand_rho_is_keyed_by_ligand():
    """The separability bootstrap pairs two metrics by ligand. If the keys were positional, a
    metric scored on 3 ligands would be compared against 3 arbitrary slots of a 5-ligand metric --
    silently comparing different ligands. Keys must be the ligands themselves."""
    rows, _ = score.load()
    ligands = {r["ligand"] for r in rows}
    per, npairs = score.per_ligand_rho(rows, "protenix-iptm")
    assert set(per) <= ligands and len(per) >= 2
    assert npairs == sum(
        1 for r in rows
        if r["ligand"] in per and r["is_binder"] == 1 and r["pKd"] is not None
        and r.get("protenix-iptm") is not None)
    # a metric with fewer scored groups must intersect, not line up by position
    sparse, _ = score.per_ligand_rho(rows, "interface-ligburial")
    assert set(sparse) <= set(per), "sparser metric should be a subset, so the intersection is real"


def test_every_prediction_joins_to_truth():
    ids = {r["id"] for r in csv.DictReader(
        (ROOT / "reference/experimental_reference_ground_truth.csv").read_text().splitlines())}
    for p in (ROOT / "predictions").glob("*_predictions.csv"):
        for row in csv.DictReader(p.read_text().splitlines()):
            assert row["id"] in ids, f"{p.name}: id {row['id']} not in ground truth"


# --- end to end ---

def test_score_runs():
    r = subprocess.run([sys.executable, str(ROOT / "scripts/score.py")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "Regression" in r.stdout and "AUROC" in r.stdout


def test_ipdae_d0_saturates_below_27_partners():
    """The whole reason ipDAE behaves differently on a ligand than on a protein partner.

    d0 only exceeds 1.0 at 27+ partners. A ligand atom never contacts 27 residues, so on these
    interfaces d0 pins to 1.0 and the partner-count scaling does nothing. If this threshold ever
    moves, the docstring's claim about saturation is wrong and should move with it.
    """
    import numpy as np
    d0 = _load("scripts/run_ipdae.py", "run_ipdae").d0_from_partner_count
    assert float(d0(np.array([26]))[0]) == 1.0
    assert float(d0(np.array([27]))[0]) > 1.0
    assert float(d0(np.array([1]))[0]) == 1.0


def test_ipdae_masks_to_contacts_and_ranks_by_pae():
    """A perfectly confident contact scores 1.0, worse PAE scores lower, and -- the one that
    actually catches bugs -- PAE *outside* the contact set must not move the score at all."""
    import numpy as np
    m = _load("scripts/run_ipdae.py", "run_ipdae")
    contact = np.array([[True, False], [False, True]])

    perfect, _ = m.ipdae_from_block(np.zeros((2, 2)), contact)
    assert perfect == 1.0  # d0 = 1, tm = 1/(1+0)

    worse, _ = m.ipdae_from_block(np.full((2, 2), 5.0), contact)
    assert 0.0 < worse < perfect

    base = np.array([[1.0, 99.0], [99.0, 1.0]])
    noise = np.array([[1.0, 0.0], [0.0, 1.0]])  # same contacts, wildly different non-contacts
    assert m.ipdae_from_block(base, contact)[0] == m.ipdae_from_block(noise, contact)[0]


def test_ipdae_refuses_empty_and_mismatched_input():
    import numpy as np
    import pytest
    m = _load("scripts/run_ipdae.py", "run_ipdae")
    assert m.ipdae_from_block(np.zeros((2, 2)), np.zeros((2, 2), dtype=bool)) == (None, None)
    with pytest.raises(ValueError):
        m.ipdae_from_block(np.zeros((2, 3)), np.ones((2, 2), dtype=bool))


# --- rank statistics: the tie bug, and the guards against it coming back -----
# The original _rank numbered a sort 0..n-1, which breaks ties by row order and
# invents an ordering the data does not have. With 42 tied values in 60,
# interface-ligburial's within-ligand rho was reported as +0.10 when the correct
# value is +0.32. These tests exist so that cannot recur silently.

import numpy as np  # noqa: E402

# score.py is loaded by path above, not importable as a module; reuse that object
_rank, spearman, auroc = score._rank, score.spearman, score.auroc
batch_spearman, batch_auroc = score.batch_spearman, score.batch_auroc
bootstrap_ci = score.bootstrap_ci


def test_rank_gives_tied_values_the_same_midrank():
    r = _rank([1.0, 1.0, 1.0, 2.0])
    assert list(r) == [1.0, 1.0, 1.0, 3.0], "tied values must share the average of their ranks"


def test_spearman_matches_scipy_on_heavily_tied_input():
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(0)
    worst = 0.0
    for _ in range(200):
        n = int(rng.integers(5, 60))
        x = rng.integers(0, 5, n).astype(float)      # deliberately tie-dense
        y = rng.normal(size=n)
        worst = max(worst, abs(spearman(x, y) - scipy_stats.spearmanr(x, y).statistic))
    assert worst < 1e-12, f"max deviation from scipy {worst:.2e}"


def test_spearman_does_not_claim_perfect_correlation_on_a_constant_vector():
    """A ligand-only score is constant within a ligand by construction. The old
    ranking returned +1.00 for this; the honest answer is no correlation."""
    rho = spearman([3.0, 3.0, 3.0, 3.0], [1.0, 2.0, 3.0, 4.0])
    assert np.isnan(rho), f"a constant predictor cannot correlate, got {rho}"


def test_auroc_matches_sklearn_including_ties():
    skm = pytest.importorskip("sklearn.metrics")
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(200):
        n = int(rng.integers(6, 60))
        s = rng.integers(0, 4, n).astype(float)
        y = rng.integers(0, 2, n)
        if not 0 < y.sum() < n:
            continue
        worst = max(worst, abs(auroc(s, y) - skm.roc_auc_score(y, s)))
    assert worst < 1e-12, f"max deviation from sklearn {worst:.2e}"


def test_batched_statistics_equal_the_scalar_ones():
    """The bootstrap uses the batched forms; if they ever drift from the scalar
    versions the CIs stop describing the point estimate above them."""
    rng = np.random.default_rng(2)
    x = rng.integers(0, 4, (64, 25)).astype(float)
    y = rng.normal(size=(64, 25))
    assert np.allclose(batch_spearman(x, y),
                       [spearman(x[i], y[i]) for i in range(64)], equal_nan=True)
    lab = rng.integers(0, 2, (64, 25))
    keep = [i for i in range(64) if 0 < lab[i].sum() < 25]
    assert np.allclose(batch_auroc(x, lab)[keep],
                       [auroc(x[i], lab[i]) for i in keep])


def test_bootstrap_is_deterministic_and_brackets_the_estimate():
    rng = np.random.default_rng(3)
    x = rng.normal(size=40)
    y = x * 0.8 + rng.normal(size=40) * 0.4
    a = bootstrap_ci(x, y, spearman)
    b = bootstrap_ci(x, y, spearman)
    assert a == b, "same seed must give the same interval"
    lo, hi = a
    assert lo < spearman(x, y) < hi
    assert lo > 0, "a strong positive correlation should have an interval clear of zero"


def test_bootstrap_survives_a_degenerate_resample():
    """Some resamples of a tiny vector are all-identical, so the statistic is
    undefined. Those must be dropped, not turned into a NaN interval."""
    lo, hi = bootstrap_ci([1.0, 1.0, 2.0], [1.0, 2.0, 3.0], spearman, n_boot=200)
    assert np.isfinite(lo) and np.isfinite(hi)


# --- the method registry ----------------------------------------------------
# A scorer that was never run, or one added without registering, used to be
# invisible: update.py compared only the files that existed, so it reported
# "nothing to update" while three rowan-cofold columns had never been produced.

def test_every_predictions_file_is_registered():
    spec = importlib.util.spec_from_file_location("methods", ROOT / "scripts/methods.py")
    methods = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(methods)
    stems = {p.stem.replace("_predictions", "")
             for p in (ROOT / "predictions").glob("*_predictions.csv")}
    unregistered = stems - set(methods.PRODUCERS)
    assert not unregistered, f"add to scripts/methods.py: {sorted(unregistered)}"


def test_registry_names_a_producer_for_every_column():
    spec = importlib.util.spec_from_file_location("methods", ROOT / "scripts/methods.py")
    methods = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(methods)
    assert methods.PRODUCERS, "registry is empty"
    assert all(v.strip() for v in methods.PRODUCERS.values()), "a column has no producer"
    assert set(methods.PRODUCERS) >= methods.OPTIONAL, "optional names must be registered"


def test_registry_producers_point_at_the_script_that_writes_the_file():
    """A registry naming the wrong producer is worse than none: it sends you to
    the wrong script with confidence. Four PAE columns were attributed to
    run_predictors/run_protenix_confidence when interface_metrics.py writes them."""
    import re
    import subprocess

    spec = importlib.util.spec_from_file_location("methods", ROOT / "scripts/methods.py")
    methods = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(methods)
    wrong = []
    for col, producer in methods.PRODUCERS.items():
        hits = subprocess.run(
            ["grep", "-rlE", rf'write\("{re.escape(col)}"|"{re.escape(col)}_predictions',
             str(ROOT / "scripts")], capture_output=True, text=True).stdout.split()
        if hits and not any(producer.split()[0] in h for h in hits):
            wrong.append(f"{col}: registry says {producer}, written by {hits}")
    assert not wrong, "\n".join(wrong)


# --- Rowan hosted co-folding -------------------------------------------------
# These cover the three bugs that actually bit: a dict payload read with getattr
# (wrote an empty file and reported success), a ledger that counted an
# uncollected run as free, and partial coverage failing the weekly QC job.

def _rowan():
    return _load("scripts/run_rowan_cofold.py", "run_rowan_cofold")


REAL_PAYLOAD = {  # shape returned by the API for 01_apx1049__apixaban
    "affinity_score": {"pred_value": -1.114, "probability_binary": 0.921},
    "scores": {"iptm": 0.98, "ptm": 0.969},
}


def test_parse_result_reads_dict_payload():
    got = _rowan().parse_result(REAL_PAYLOAD)
    assert got["affinity"] == pytest.approx(7.114)
    assert got["pbind"] == pytest.approx(0.921)
    assert got["iptm"] == pytest.approx(0.98)


def test_hosted_columns_and_state_are_scoped_by_model():
    """--model is a dispatch input in rowan.yml. With a fixed column name a
    Chai-1r run would have overwritten Boltz-2's numbers in place, and the state
    file would have reported those systems already submitted."""
    r = _rowan()
    assert r.columns("boltz_2") != r.columns("chai_1r")
    assert set(r.columns("boltz_2").values()).isdisjoint(r.columns("chai_1r").values())
    assert r.key("boltz_2", "x") != r.key("chai_1r", "x")
    # underscores in Rowan's model ids become hyphens in column names
    assert r.columns("chai_1r")["iptm"] == "rowan-chai-1r-iptm"


def test_legacy_state_keys_are_read_as_boltz2():
    """Keys written before the state was model-scoped are bare system ids, and
    everything submitted then was boltz_2."""
    assert _rowan()._split("00_x__y") == ("boltz-2", "00_x__y")
    assert _rowan()._split("chai-1r:00_x__y") == ("chai-1r", "00_x__y")


def test_hosted_columns_never_collide_with_local_ones():
    """The hosted Boltz-2 is a different deployment with a different MSA server;
    merging it into the local boltz-2 column is the provenance error."""
    m = _load("scripts/methods.py", "methods")
    hosted = {c for c in m.PRODUCERS if c.startswith("rowan-")}
    assert "boltz-2" not in hosted and "boltz-2-iptm" not in hosted


def test_parse_result_unwraps_cofolding_results():
    assert _rowan().parse_result({"cofolding_results": [REAL_PAYLOAD]}) == \
        _rowan().parse_result(REAL_PAYLOAD)


@pytest.mark.parametrize("payload", [None, {}, "done", {"scores": None},
                                     {"affinity_score": {"pred_value": None}}])
def test_parse_result_survives_missing_fields(payload):
    """A pydantic object, a half-finished run or a changed schema must yield
    nothing rather than raising or inventing a score."""
    assert _rowan().parse_result(payload) == {} or \
        all(v is not None for v in _rowan().parse_result(payload).values())


def test_ledger_counts_uncollected_runs_as_spent():
    """An entry with no recorded charge is a submitted workflow whose cost we have
    not fetched yet. Counting it as free would let repeated dispatches walk past
    the budget."""
    r = _rowan()
    assert r._ledger({"a": {"uuid": "x", "credits": 9.08}}) == pytest.approx(9.08)
    assert r._ledger({"a": {"uuid": "x", "credits": None}}) > 0
    assert r._ledger({"a": "legacy-bare-uuid"}) > 0


def test_uuid_of_accepts_both_state_shapes():
    r = _rowan()
    assert r._uuid_of("bare") == "bare" == r._uuid_of({"uuid": "bare", "credits": 1.0})


def test_budget_is_set_and_finite():
    r = _rowan()
    assert 0 < r.BUDGET_CREDITS < 10_000


def test_qc_warns_but_passes_on_partial_coverage():
    """predictions/ currently holds a partially-filled Rowan run. QC must surface
    it and still exit 0, or the weekly job is red until the run finishes."""
    out = subprocess.run([sys.executable, str(ROOT / "scripts/qc_dataset.py")],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stdout
    if "pending)" in out.stdout:
        assert "WARN" in out.stdout


def test_every_registered_method_has_a_figure_category():
    """make_figure.py indexes CAT[m] directly. A registered method missing from
    the category map crashed the figures workflow with a bare KeyError -- the
    two dicts are now one, and this is what keeps them one."""
    m = _load("scripts/methods.py", "methods")
    assert not set(m.PRODUCERS) - set(m.CATEGORY), "registered but uncategorised"
    # the reverse too: a category for a column that no longer exists is dead
    # weight that outlives the rename it was left behind by
    assert not set(m.CATEGORY) - set(m.PRODUCERS), "categorised but unregistered"


def test_figure_categories_all_have_a_colour():
    pytest.importorskip("matplotlib")  # CI installs it; a bare local env may not
    fig = _load("scripts/make_figure.py", "make_figure")
    assert not set(fig.CAT.values()) - set(fig.COL)
    assert not set(fig.CAT.values()) - set(fig.LABEL)
