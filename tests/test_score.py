"""Pure-logic tests for the scorer and a data-integrity check on the shipped predictions."""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import score


def test_spearman():
    assert score.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0
    assert score.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0


def test_pearson():
    assert abs(score.pearson([1, 2, 3], [2, 4, 6]) - 1.0) < 1e-9
    assert abs(score.pearson([1, 2, 3], [6, 4, 2]) + 1.0) < 1e-9


def test_rmse():
    assert score.rmse([1, 2, 3], [1, 2, 3]) == 0.0
    assert abs(score.rmse([0, 0], [3, 4]) - math.sqrt(12.5)) < 1e-9


def test_auroc():
    assert score.auroc([2, 2, 1, 1], [1, 1, 0, 0]) == 1.0       # perfect separation
    assert score.auroc([1, 1, 2, 2], [1, 1, 0, 0]) == 0.0       # perfectly wrong
    assert score.auroc([1, 1], [1, 1]) is None                  # no negatives


ROOT = Path(__file__).resolve().parents[1]


def _csv(path):
    return list(csv.DictReader(path.read_text().splitlines()))


def test_reference_integrity():
    gt = _csv(ROOT / "reference/experimental_reference_ground_truth.csv")
    sys_ref = _csv(ROOT / "reference/system_reference.csv")
    assert len(gt) == 47 and len(sys_ref) == 47
    assert sum(r["is_binder"] == "1" for r in gt) == 37
    assert sum(r["is_binder"] == "0" for r in gt) == 10
    ids = [r["id"] for r in gt]
    assert len(set(ids)) == len(ids)                            # unique ids
    assert {r["id"] for r in sys_ref} == set(ids)               # the two reference files align
    assert all(r["source_url"].strip() for r in sys_ref)        # every system keeps its reference


def test_predictions_strict_and_complete():
    ids = {r["id"] for r in _csv(ROOT / "reference/experimental_reference_ground_truth.csv")}
    files = list((ROOT / "predictions").glob("*_predictions.csv"))
    assert files
    for p in files:
        rows = _csv(p)
        assert list(rows[0].keys()) == ["method", "id", "predicted_affinity"], p.name  # strict 3-col
        assert {r["id"] for r in rows} == ids, p.name                                   # covers all ids
        assert all(r["predicted_affinity"].strip() for r in rows), p.name
