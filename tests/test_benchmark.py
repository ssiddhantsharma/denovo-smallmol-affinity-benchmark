"""Tests for the scoring statistics and dataset integrity. No GPU, no network."""

import csv
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_score():
    spec = importlib.util.spec_from_file_location("score", ROOT / "scripts/score.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


score = _load_score()


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
    assert len(rows) == 61
    binders = [r for r in rows if r["is_binder"] == 1]
    negs = [r for r in rows if r["is_binder"] == 0]
    assert len(binders) == 50
    assert len(negs) == 11
    assert all(r["pKd"] is not None for r in binders), "every binder must carry a measured pKd"
    assert len(methods) >= 16
    for m in ("boltz-2", "clogp", "protenix-iptm", "dtsfm-cosine"):
        assert m in methods


def test_within_ligand_drops_constant_ligand_properties():
    rows, _ = score.load()
    # cLogP is a ligand-only property: constant within a ligand group, so it yields no within-ligand rho
    med, _per, n = score.within_ligand_rho(rows, "clogp")
    assert med is None and n == 0
    # an interface metric varies across the proteins binding one ligand, so it does
    med2, _per2, n2 = score.within_ligand_rho(rows, "protenix-iptm")
    assert n2 >= 2 and med2 is not None


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
